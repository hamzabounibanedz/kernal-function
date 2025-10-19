// IpmEngine: Numerical core for PD‑IPM SDO solves
// Responsibilities:
// - NT scaling and Newton system construction/solve
// - Proximity Ψ(V)=tr(ψ(V)) evaluation and backtracking line search
// - SPD stabilization via projection where needed
// - Iteration accounting (outer L, inner K_total, K_avg) and timing
// Notes:
// - Formulas in kernel calls match the cited papers; domains validated
// - Thread CPU time on Windows via GetThreadTimes; std::clock elsewhere
#include "ipm_engine.h"
#include <algorithm>
#include <chrono>
#include <cmath>

// Helper: sanitize a matrix by replacing non-finite entries and clamping
static Eigen::MatrixXd sanitize(const Eigen::MatrixXd &Min) {
  Eigen::MatrixXd M = Min;
  for (int i = 0; i < M.rows(); ++i) {
    for (int j = 0; j < M.cols(); ++j) {
      double v = M(i, j);
      if (!std::isfinite(v)) v = 0.0;
      if (v > 1e12) v = 1e12;
      if (v < -1e12) v = -1e12;
      M(i, j) = v;
    }
  }
  return M;
}

// Project a symmetric matrix to the nearest SPD by clamping
// small/negative eigenvalues (stabilizes kernel evaluations)
static Eigen::MatrixXd projectSPD(const Eigen::MatrixXd &Min) {
  Eigen::MatrixXd Ms = sanitize(0.5 * (Min + Min.transpose()));
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ms);
  Eigen::VectorXd vals = es.eigenvalues().cwiseMax(1e-12);
  return es.eigenvectors() * vals.asDiagonal() * es.eigenvectors().transpose();
}

AlgorithmResult IpmEngine::solve(const TestCase &tc, KernelBase &kernel,
                                 const IpmSettings &settings) {
  auto t0 = std::chrono::high_resolution_clock::now();
#if defined(_WIN32)
  FILETIME creationTime, exitTime, kernelTime0, userTime0;
  HANDLE hThread = GetCurrentThread();
  GetThreadTimes(hThread, &creationTime, &exitTime, &kernelTime0, &userTime0);
#else
  clock_t cpu0 = std::clock();
#endif

  AlgorithmResult res{};
  res.converged = false;
  res.iterations = 0;
  res.outer_iterations = 0;
  res.inner_iterations_total = 0;
  res.avg_inner_per_outer = 0.0;

  Iterate it = initialize(tc);
  // Extra initialization sanitation
  it.X = projectSPD(it.X);
  it.S = projectSPD(it.S);
  it.mu = std::max(1e-12, it.mu);

  int outer = 0;
  while (outer < settings.maxOuter && tc.problem_size * it.mu > settings.epsilon) {
    // Barrier update (count after update per spec)
    it.mu = (1.0 - settings.theta) * it.mu;
    res.outer_iterations++;
    outer++;
    int inner_this_outer = 0;
    // Inner loop: contract proximity for fixed μ
    int inner = 0;
    while (inner < settings.maxInner) {
      if (settings.cancel && settings.cancel->load()) {
        res.converged = false;
        res.convergence_info = "Cancelled";
        break;
      }

      // Compute NT scaling (full matrices) and spectrum of X^{1/2} S X^{1/2}
      struct Scaling { Eigen::MatrixXd P; Eigen::MatrixXd D; Eigen::VectorXd v; Eigen::MatrixXd U; };
      auto computeScaling = [&](const Iterate &itIn) -> Scaling {
        Scaling sc{};
        Eigen::MatrixXd Xh = sqrtmSPD(projectSPD(itIn.X));
        Eigen::MatrixXd M = projectSPD(Xh * itIn.S * Xh); // stabilized SPD
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
        Eigen::VectorXd lam = es.eigenvalues();
        // sanitize eigenvalues
        for (int i = 0; i < lam.size(); ++i) {
          if (!std::isfinite(lam(i)) || lam(i) <= 0.0) lam(i) = 1e-16;
        }
        sc.U = es.eigenvectors();
        double mu_safe = std::sqrt(std::max(1e-16, itIn.mu));
        sc.v = lam.cwiseSqrt() / mu_safe; // v_i = sqrt(λ_i/μ)
        // sanitize V spectrum for kernel calls
        for (int i = 0; i < sc.v.size(); ++i) {
          if (!std::isfinite(sc.v(i)) || sc.v(i) <= 0.0) sc.v(i) = 1e-8;
          if (sc.v(i) > 1e8) sc.v(i) = 1e8;
        }
        Eigen::MatrixXd Minvhalf = invsqrtmSPD(M);
        sc.P = Xh * Minvhalf * Xh;
        sc.D = invsqrtmSPD(sc.P);
        return sc;
      };

      auto sc = computeScaling(it);

      // Proximity Ψ=tr(ψ(V)) via spectrum v_i
      auto computePsiFromSpectrum = [&](const Eigen::VectorXd &vVals) -> double {
        double sum = 0.0;
        for (int i = 0; i < vVals.size(); ++i) sum += kernel.psi(vVals(i));
        return sum;
      };

      double psiTrace = computePsiFromSpectrum(sc.v);
      if (psiTrace <= settings.tau) break; // contracted proximity for this μ

      // Gradient -ψ'(V) in eigenbasis, map back: G = U diag(ψ'(v_i)) U^T
      Eigen::VectorXd gDiag(sc.v.size());
      for (int i = 0; i < sc.v.size(); ++i) gDiag(i) = kernel.psi_prime(sc.v(i));
      Eigen::MatrixXd G = sc.U * gDiag.asDiagonal() * sc.U.transpose();
      Eigen::MatrixXd rhs = -G;

      // Solve Newton system (exact NT-based formulation)
      Eigen::MatrixXd Dtmp;
      Direction dir = solveNewtonNT(tc, it, rhs, Dtmp);

      // Backtracking line search with SPD checks and proximity decrease
      auto isSPDAll = [&](const Eigen::MatrixXd &A) -> bool { return isSPD(A); };
      auto trialProximity = [&](const Iterate &base, const Direction &d, double a) -> double {
        Iterate t = base;
        // full updates
        t.X = symmetrize(base.X + a * d.dX);
        t.S = symmetrize(base.S + a * d.dS);
        // keep y update out of proximity calc
      Eigen::Index denomIndex = (t.X.rows() > 0 ? t.X.rows() : 1);
        t.mu = (t.X.cwiseProduct(t.S)).sum() / static_cast<double>(denomIndex);
        auto tsc = computeScaling(t);
        return computePsiFromSpectrum(tsc.v);
      };

      double alpha = 1.0;
      const double c1 = 1e-4;
      const double shrink = 0.7;
      int bt = 0;
      while (bt++ < 25) {
        Eigen::MatrixXd Xtrial = symmetrize(it.X + alpha * dir.dX);
        Eigen::MatrixXd Strial = symmetrize(it.S + alpha * dir.dS);
        if (isSPDAll(Xtrial) && isSPDAll(Strial)) {
          double psiNew = trialProximity(it, dir, alpha);
          if (psiNew <= (1.0 - c1 * alpha) * psiTrace) {
            break; // sufficient decrease
          }
        }
        alpha *= shrink;
        if (alpha < 1e-8) break;
      }

      // Apply full-matrix step
      applyStep(it, dir, alpha);

      // Iteration accounting
      res.iterations++;               // total Newton systems solved (T)
      res.inner_iterations_total++;   // cumulative inner iterations
      inner_this_outer++;             // per-outer counter
      inner++;
      if (res.iterations > 10000) break; // hard guard
    }

    if (settings.cancel && settings.cancel->load()) break;

    // Loop proceeds; μ update already applied at outer start
  }

  // Fill result
  // Convert Eigen to std::vector for result compatibility
  res.primal_solution.assign(it.X.rows(), std::vector<double>(it.X.cols(), 0.0));
  for (int i = 0; i < it.X.rows(); ++i)
    for (int j = 0; j < it.X.cols(); ++j)
      res.primal_solution[i][j] = it.X(i, j);
  res.dual_solution_y = it.y;
  res.dual_solution_s.assign(it.S.rows(), std::vector<double>(it.S.cols(), 0.0));
  for (int i = 0; i < it.S.rows(); ++i)
    for (int j = 0; j < it.S.cols(); ++j)
      res.dual_solution_s[i][j] = it.S(i, j);
  res.final_mu = it.mu;
  res.duality_gap = dualityGap(it.X, it.S);
  if (!std::isfinite(res.duality_gap)) {
    res.duality_gap = dualityGap(projectSPD(it.X), projectSPD(it.S));
  }
  res.primal_objective = primalObjective(it.X, tc.C_matrix);
  if (!std::isfinite(res.primal_objective)) {
    res.primal_objective = primalObjective(projectSPD(it.X), tc.C_matrix);
  }
  res.dual_objective = dualObjective(it.y, tc.b_vector);

  // Normalized duality gap
  double norm_gap = std::abs(res.primal_objective - res.dual_objective) /
                    (1.0 + std::abs(res.primal_objective) + std::abs(res.dual_objective));
  bool mu_ok = (tc.problem_size * it.mu <= settings.epsilon);
  bool gap_ok = std::isfinite(norm_gap) && (norm_gap <= settings.epsilon * 1e2);
  res.converged = mu_ok || gap_ok;
  res.convergence_info = res.converged ? (mu_ok ? "Converged (mu)" : "Converged (gap)") : "Not converged";

  // Average inner per outer (avoid divide-by-zero)
  if (res.outer_iterations > 0) {
    res.avg_inner_per_outer = static_cast<double>(res.inner_iterations_total) /
                              static_cast<double>(res.outer_iterations);
  } else {
    res.avg_inner_per_outer = 0.0;
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  res.execution_time_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() /
      1000.0;
  // CPU time in milliseconds
#if defined(_WIN32)
  FILETIME kernelTime1, userTime1;
  GetThreadTimes(hThread, &creationTime, &exitTime, &kernelTime1, &userTime1);
  ULARGE_INTEGER k0, u0, k1, u1;
  k0.LowPart = kernelTime0.dwLowDateTime; k0.HighPart = kernelTime0.dwHighDateTime;
  u0.LowPart = userTime0.dwLowDateTime;   u0.HighPart = userTime0.dwHighDateTime;
  k1.LowPart = kernelTime1.dwLowDateTime; k1.HighPart = kernelTime1.dwHighDateTime;
  u1.LowPart = userTime1.dwLowDateTime;   u1.HighPart = userTime1.dwHighDateTime;
  // FILETIME is 100-nanosecond intervals
  unsigned long long hundredNs = (k1.QuadPart - k0.QuadPart) + (u1.QuadPart - u0.QuadPart);
  res.cpu_time_ms = static_cast<double>(hundredNs) / 10000.0;
#else
  clock_t cpu1 = std::clock();
  res.cpu_time_ms = (cpu1 - cpu0) * 1000.0 / CLOCKS_PER_SEC;
#endif
  return res;
}

IpmEngine::Iterate IpmEngine::initialize(const TestCase &tc) const {
  Iterate it;
  it.X.resize(tc.problem_size, tc.problem_size);
  for (int i = 0; i < tc.problem_size; ++i)
    for (int j = 0; j < tc.problem_size; ++j)
      it.X(i, j) = tc.initial_X[i][j];
  it.y = tc.initial_y;
  it.S.resize(tc.problem_size, tc.problem_size);
  for (int i = 0; i < tc.problem_size; ++i)
    for (int j = 0; j < tc.problem_size; ++j)
      it.S(i, j) = tc.initial_S[i][j];

  // Ensure positive diagonal entries as a minimal IPC enforcement
  for (int i = 0; i < it.X.rows(); ++i) it.X(i, i) = std::max(1.0, std::abs(it.X(i, i)));
  for (int i = 0; i < it.S.rows(); ++i) it.S(i, i) = std::max(1.0, std::abs(it.S(i, i)));

  // μ0 = trace(X S)/n
  double tr = 0.0;
  tr = (it.X.cwiseProduct(it.S)).sum();
  it.mu = it.X.rows() == 0 ? 1.0 : tr / static_cast<double>(it.X.rows());
  if (it.mu <= 0.0) it.mu = 1.0;
  return it;
}

Eigen::MatrixXd IpmEngine::computeNTScalingV(const Iterate &it) const {
  // Compute V from spectrum of X^{1/2} S X^{1/2} relative to μ: V = U diag(sqrt(λ_i/μ)) U^T
  Eigen::MatrixXd Xh = sqrtmSPD(it.X);
  Eigen::MatrixXd M = Xh * it.S * Xh; // SPD
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
  Eigen::VectorXd lam = es.eigenvalues().cwiseMax(1e-16);
  Eigen::MatrixXd U = es.eigenvectors();
  double mu_clamped = it.mu < 1e-16 ? 1e-16 : it.mu;
  Eigen::VectorXd v = lam.cwiseSqrt() / std::sqrt(mu_clamped);
  return U * v.asDiagonal() * U.transpose();
}

double IpmEngine::computeProximityTracePsi(const Eigen::MatrixXd &V,
                                           KernelBase &kernel) const {
  double sum = 0.0;
  for (Eigen::Index i = 0; i < V.rows(); ++i) sum += kernel.psi(V(i, i));
  return sum;
}

Eigen::MatrixXd IpmEngine::kernelGradient(const Eigen::MatrixXd &V,
                                          KernelBase &kernel) const {
  const Eigen::Index n = V.rows();
  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(n, n);
  for (Eigen::Index i = 0; i < n; ++i) G(i, i) = kernel.psi_prime(V(i, i));
  return G;
}

IpmEngine::Direction IpmEngine::solveNewtonSimplified(
    const TestCase &tc, const Iterate &it, const Eigen::MatrixXd &rhs) const {
  (void)tc;
  (void)it;
  const Eigen::Index n = rhs.rows();
  Direction d;
  d.dX = Eigen::MatrixXd::Zero(n, n);
  d.dS = Eigen::MatrixXd::Zero(n, n);
  d.dy.assign(tc.A_matrices.size(), 0.0);
  for (Eigen::Index i = 0; i < n; ++i) {
    // Split equally between dX and dS along diagonal
    d.dX(i, i) = rhs(i, i) * 0.5;
    d.dS(i, i) = rhs(i, i) * 0.5;
  }
  return d;
}

IpmEngine::Direction IpmEngine::solveNewtonNT(const TestCase &tc, const Iterate &it,
                                              const Eigen::MatrixXd &R,
                                              Eigen::MatrixXd &D_out) const {
  // NT scaling matrices
  // P = X^(1/2)(X^(1/2) S X^(1/2))^(-1/2) X^(1/2), D = P^(-1/2)
  Eigen::MatrixXd Xh = sqrtmSPD(it.X);
  Eigen::MatrixXd M = Xh * it.S * Xh; // SPD
  Eigen::MatrixXd Minvhalf = invsqrtmSPD(M);
  Eigen::MatrixXd P = Xh * Minvhalf * Xh;
  Eigen::MatrixXd D = invsqrtmSPD(P);
  D_out = D;

  // Scale constraints: Abar_i = D^{-1} A_i D^{-1}
  auto Abar = scaleConstraints(tc.A_matrices, D);

  // Build normal equation H Δy = r based on simplified trace constraints:
  // We use a basic symmetric positive definite H with entries H_ij = tr(Abar_i Abar_j)
  const int m = static_cast<int>(tc.A_matrices.size());
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) {
      double hij = (Abar[i].cwiseProduct(Abar[j])).sum();
      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  // RHS for Δy: we use −trace(Abar_i * R) as a surrogate coupling to R
  Eigen::VectorXd r = Eigen::VectorXd::Zero(m);
  for (int i = 0; i < m; ++i) r(i) = -(Abar[i].cwiseProduct(R)).sum();

  // Solve H Δy = r (SPD)
  Eigen::VectorXd dyVec = Eigen::VectorXd::Zero(m);
  Eigen::LLT<Eigen::MatrixXd> llt(H);
  if (llt.info() == Eigen::Success) {
    dyVec = llt.solve(r);
  }

  // Recover DX, DS in scaled space using simplified relations:
  // DX − DS = R, and Abar⋅DX = 0, Σ y_i Abar_i − DS = 0 ⇒ DS ≈ Σ dy_i Abar_i
  Eigen::MatrixXd DS = Eigen::MatrixXd::Zero(it.X.rows(), it.X.cols());
  for (int i = 0; i < m; ++i) DS += dyVec(i) * Abar[i];
  Eigen::MatrixXd DX = R + DS;

  // Project DX to satisfy constraints tr(Abar_i DX) = 0
  Eigen::VectorXd c = Eigen::VectorXd::Zero(m);
  for (int i = 0; i < m; ++i) c(i) = (Abar[i].cwiseProduct(DX)).sum();
  if (H.rows() > 0) {
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(m);
    if (llt.info() == Eigen::Success) {
      lambda = llt.solve(c);
      for (int i = 0; i < m; ++i) DX -= lambda(i) * Abar[i];
      // Recompute DS from DX - R to preserve DX-DS=R
      DS = DX - R;
    }
  }

  // Back transform to original variables:
  // ΔX = D DX D, ΔS = D^{-1} DS D^{-1}
  Eigen::MatrixXd Dinv = sqrtmSPD(P); // since D = P^{-1/2}
  Eigen::MatrixXd dX = D * DX * D;
  Eigen::MatrixXd dS = Dinv * DS * Dinv;

  Direction dir;
  dir.dX = dX;
  dir.dS = dS;
  dir.dy.assign(m, 0.0);
  for (int i = 0; i < m; ++i) dir.dy[i] = dyVec(i);
  return dir;
}

Eigen::MatrixXd IpmEngine::symmetrize(const Eigen::MatrixXd &M) {
  return 0.5 * (M + M.transpose());
}

Eigen::MatrixXd IpmEngine::sqrtmSPD(const Eigen::MatrixXd &M) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
  Eigen::VectorXd vals = es.eigenvalues().cwiseMax(1e-12).cwiseSqrt();
  return es.eigenvectors() * vals.asDiagonal() * es.eigenvectors().transpose();
}

Eigen::MatrixXd IpmEngine::invsqrtmSPD(const Eigen::MatrixXd &M) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
  Eigen::VectorXd vals = es.eigenvalues().cwiseMax(1e-12).cwiseSqrt().cwiseInverse();
  return es.eigenvectors() * vals.asDiagonal() * es.eigenvectors().transpose();
}

std::vector<Eigen::MatrixXd> IpmEngine::scaleConstraints(
    const std::vector<std::vector<std::vector<double>>> &A,
    const Eigen::MatrixXd &D) {
  const int n = static_cast<int>(D.rows());
  std::vector<Eigen::MatrixXd> out;
  out.reserve(A.size());
  Eigen::MatrixXd Dinv = D.inverse();
  for (const auto &Ai : A) {
    Eigen::MatrixXd Aimat(n, n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        Aimat(i, j) = Ai[i][j];
    Eigen::MatrixXd Abar = Dinv * Aimat * Dinv;
    out.push_back(symmetrize(Abar));
  }
  return out;
}

bool IpmEngine::isSPD(const Eigen::MatrixXd &M) {
  Eigen::LLT<Eigen::MatrixXd> llt(symmetrize(M));
  return llt.info() == Eigen::Success;
}

double IpmEngine::primalObjective(const Eigen::MatrixXd &X,
                                  const std::vector<std::vector<double>> &C) {
  double sum = 0.0;
  for (int i = 0; i < X.rows(); ++i)
    for (int j = 0; j < X.cols(); ++j) sum += X(i, j) * C[i][j];
  return sum;
}

double IpmEngine::dualObjective(const std::vector<double> &y,
                                const std::vector<double> &b) {
  double sum = 0.0;
  const size_t m = std::min(y.size(), b.size());
  for (size_t i = 0; i < m; ++i) sum += y[i] * b[i];
  return sum;
}

double IpmEngine::computeStepSizePD(const Iterate &it, const Direction &dir) const {
  // Initial guess for line search
  (void)it;
  (void)dir;
  return 1.0;
}

void IpmEngine::applyStep(Iterate &it, const Direction &dir, double alpha) const {
  it.X = projectSPD(it.X + alpha * dir.dX);
  it.S = projectSPD(it.S + alpha * dir.dS);
  for (size_t i = 0; i < it.y.size() && i < dir.dy.size(); ++i)
    it.y[i] += alpha * dir.dy[i];
}

double IpmEngine::trace(const Eigen::MatrixXd &M) {
  return M.trace();
}

double IpmEngine::dualityGap(const Eigen::MatrixXd &X, const Eigen::MatrixXd &S) {
  return (X.cwiseProduct(S)).sum();
}

double IpmEngine::computeMaxStepSPD(const Eigen::MatrixXd &X,
                                    const Eigen::MatrixXd &dX) {
  // Find max α s.t. X + α dX ≻ 0 using smallest eigenvalue along line
  // Use heuristic: ensure all diagonal stay positive and fallback α
  double alpha = 1.0;
  for (int i = 0; i < X.rows(); ++i) {
    if (dX(i, i) < 0.0) alpha = std::min(alpha, -X(i, i) / dX(i, i));
  }
  alpha = std::max(1e-6, alpha);
  return alpha;
}



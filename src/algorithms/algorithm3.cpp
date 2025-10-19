#include "algorithm3.h"
#include "../kernels/parameterized_log_kernel.h"
#include "../utils/test_data.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

Algorithm3::Algorithm3()
    : AlgorithmBase(std::make_shared<ParameterizedLogKernel>()) {
  m_logKernel = std::make_shared<ParameterizedLogKernel>(1.0);
  // Set default parameters
  m_tolerance = 1e-8;
  m_theta = 0.5;
  m_tau = 3.0;
}

AlgorithmResult Algorithm3::solve(const TestCase &test_case) {
  auto start_time = std::chrono::high_resolution_clock::now();

  AlgorithmResult result;
  result.converged = false;
  result.iterations = 0;

  try {
    // Initialize the algorithm
    CurrentIterate current = initialize(test_case);

    // Outer loop: while n·μ ≥ ε do
    while (!checkConvergence(current.mu, test_case.tolerance,
                             test_case.problem_size)) {
      // μ ← (1−θ)·μ
      current.mu = (1.0 - test_case.theta) * current.mu;

      // V ← V / √(1−θ)
      updateScalingMatrix(current, test_case);

      // Inner loop: while Ψ(V) > τ do
      while (
          checkProximity(computeProximityMeasure(current.V), test_case.tau)) {
        // Compute kernel gradient -ψ′(V)
        std::vector<std::vector<double>> rhs_psi =
            computeKernelGradient(current.V);

        // Solve linear system for search direction
        SearchDirection direction =
            solveLinearSystem(test_case, current, rhs_psi);

        // Choose step size α (e.g. line-search)
        double alpha = computeStepSize(current, direction);

        // Update iterate
        updateIterate(current, direction, alpha);

        // Update V ← (1/√μ)·D⁻¹ X D⁻¹ = (1/√μ)·D S D
        updateVMatrix(current);

        result.iterations++;

        // Safety check to prevent infinite loops
        if (result.iterations > 1000) {
          result.convergence_info = "Maximum iterations exceeded";
          break;
        }
      }
    }

    // Compute final results
    result.converged = checkConvergence(current.mu, test_case.tolerance,
                                        test_case.problem_size);
    result.final_mu = current.mu;
    result.primal_solution = current.X;
    result.dual_solution_y = current.y;
    result.dual_solution_s = current.S;
    result.primal_objective =
        computeObjectiveValue(current.X, test_case.C_matrix);
    result.duality_gap = computeDualityGap(current.X, current.S);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
        end_time - start_time);
    result.execution_time_ms = duration.count() / 1000.0;

    if (result.converged) {
      result.convergence_info = "Algorithm converged successfully";
    }

  } catch (const std::exception &e) {
    result.convergence_info = std::string("Error: ") + e.what();
  }

  return result;
}

std::string Algorithm3::getName() const {
  return "Algorithm 3: PD IPM with Log-Kernel";
}

std::string Algorithm3::getDescription() const {
  return "Generic Primal-Dual Interior-Point Method for SDO using "
         "Parameterized Logarithmic "
         "kernel function ψ(t) = (t²-1-ln t)/2 + (e^{1/t^q-1}-1)/(2q), t>0, "
         "q≥1. "
         "Implements the algorithm from Section 2.2 of the third research "
         "paper.";
}

std::vector<TestCase> Algorithm3::getTestCases() const {
  return {TestDataProvider::getParameterizedLogTestCase()};
}

Algorithm3::CurrentIterate
Algorithm3::initialize(const TestCase &test_case) const {
  CurrentIterate current;
  int n = test_case.problem_size;

  // Initialize X ← Iₙ
  current.X = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    current.X[i][i] = 1.0;
  }

  // Initialize S ← Iₙ
  current.S = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    current.S[i][i] = 1.0;
  }

  // Initialize y
  current.y = test_case.initial_y;

  // Initialize μ ← 1
  current.mu = 1.0;

  // Initialize V ← Iₙ
  current.V = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    current.V[i][i] = 1.0;
  }

  return current;
}

double Algorithm3::computeProximityMeasure(
    const std::vector<std::vector<double>> &V) const {
  // Compute Ψ(V) = trace(ψ(V))
  double proximity = 0.0;
  for (size_t i = 0; i < V.size(); ++i) {
    for (size_t j = 0; j < V[i].size(); ++j) {
      if (i == j) { // Only diagonal elements for trace
        proximity += m_logKernel->psi(V[i][j]);
      }
    }
  }
  return proximity;
}

void Algorithm3::updateScalingMatrix(CurrentIterate &current,
                                     const TestCase &test_case) const {
  // V ← V / √(1−θ)
  double scale_factor = 1.0 / std::sqrt(1.0 - test_case.theta);
  for (auto &row : current.V) {
    for (double &val : row) {
      val *= scale_factor;
    }
  }
}

std::vector<std::vector<double>> Algorithm3::computeKernelGradient(
    const std::vector<std::vector<double>> &V) const {
  int n = V.size();
  std::vector<std::vector<double>> gradient(n, std::vector<double>(n, 0.0));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      gradient[i][j] = -m_logKernel->psi_prime(V[i][j]);
    }
  }

  return gradient;
}

Algorithm3::SearchDirection Algorithm3::solveLinearSystem(
    const TestCase &test_case, const CurrentIterate &current,
    const std::vector<std::vector<double>> &rhs_psi) const {
  (void)current; // Suppress unused parameter warning
  SearchDirection direction;
  int n = test_case.problem_size;
  int m = test_case.A_matrices.size();

  // Initialize search directions
  direction.delta_X =
      std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
  direction.delta_y.resize(m, 0.0);
  direction.delta_S =
      std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

  // Simplified implementation: solve ΔX + ΔS = −ψ′(V)
  // In practice, this would solve the complete system
  // tr(AᵢΔX) = 0, i=1…m
  // ∑ᵢ Δyᵢ Aᵢ + ΔS = 0
  // ΔX + ΔS = −ψ′(V)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      // For simplicity, distribute rhs_psi equally between ΔX and ΔS
      double delta_val = rhs_psi[i][j] / 2.0;
      direction.delta_X[i][j] = delta_val;
      direction.delta_S[i][j] = delta_val;
    }
  }

  return direction;
}

double Algorithm3::computeStepSize(const CurrentIterate &current,
                                   const SearchDirection &direction) const {
  // Simple line search: find maximum step size that maintains positive
  // definiteness
  double alpha = 1.0;

  // Check for X + αΔX > 0
  for (size_t i = 0; i < current.X.size(); ++i) {
    for (size_t j = 0; j < current.X[i].size(); ++j) {
      if (direction.delta_X[i][j] < 0) {
        alpha =
            std::min(alpha, -0.9 * current.X[i][j] / direction.delta_X[i][j]);
      }
    }
  }

  // Check for S + αΔS > 0
  for (size_t i = 0; i < current.S.size(); ++i) {
    for (size_t j = 0; j < current.S[i].size(); ++j) {
      if (direction.delta_S[i][j] < 0) {
        alpha =
            std::min(alpha, -0.9 * current.S[i][j] / direction.delta_S[i][j]);
      }
    }
  }

  return std::min(alpha, 0.95); // Conservative step size
}

void Algorithm3::updateIterate(CurrentIterate &current,
                               const SearchDirection &direction,
                               double alpha) const {
  // Update X ← X + α ΔX
  for (size_t i = 0; i < current.X.size(); ++i) {
    for (size_t j = 0; j < current.X[i].size(); ++j) {
      current.X[i][j] += alpha * direction.delta_X[i][j];
    }
  }

  // Update y ← y + α Δy
  for (size_t i = 0; i < current.y.size(); ++i) {
    current.y[i] += alpha * direction.delta_y[i];
  }

  // Update S ← S + α ΔS
  for (size_t i = 0; i < current.S.size(); ++i) {
    for (size_t j = 0; j < current.S[i].size(); ++j) {
      current.S[i][j] += alpha * direction.delta_S[i][j];
    }
  }
}

void Algorithm3::updateVMatrix(CurrentIterate &current) const {
  // Simplified implementation: V ← (1/√μ)·X·S
  // In practice, this would use more sophisticated scaling
  int n = current.X.size();
  std::vector<std::vector<double>> V(n, std::vector<double>(n, 0.0));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += current.X[i][k] * current.S[k][j];
      }
      V[i][j] = sum / std::sqrt(current.mu);
    }
  }

  current.V = V;
}

bool Algorithm3::checkConvergence(double mu, double epsilon, int n) const {
  // Converged when n·μ ≤ ε
  return n * mu <= epsilon;
}

bool Algorithm3::checkProximity(double proximity, double tau) const {
  return proximity > tau;
}

double Algorithm3::computeObjectiveValue(
    const std::vector<std::vector<double>> &X,
    const std::vector<std::vector<double>> &C) const {
  // Compute trace(C·X)
  double objective = 0.0;
  for (size_t i = 0; i < X.size() && i < C.size(); ++i) {
    for (size_t j = 0; j < X[i].size() && j < C[i].size(); ++j) {
      objective += C[i][j] * X[i][j];
    }
  }
  return objective;
}

double
Algorithm3::computeDualityGap(const std::vector<std::vector<double>> &X,
                              const std::vector<std::vector<double>> &S) const {
  // Compute trace(X·S)
  double gap = 0.0;
  for (size_t i = 0; i < X.size(); ++i) {
    for (size_t j = 0; j < X[i].size(); ++j) {
      gap += X[i][j] * S[i][j];
    }
  }
  return gap;
}

bool Algorithm3::isStrictlyFeasible(
    const std::vector<std::vector<double>> &X,
    const std::vector<std::vector<double>> &S) const {
  // Check if matrices are positive definite (simplified)
  for (const auto &row : X) {
    for (double val : row) {
      if (val <= 0)
        return false;
    }
  }
  for (const auto &row : S) {
    for (double val : row) {
      if (val <= 0)
        return false;
    }
  }
  return true;
}

double
Algorithm3::matrixTrace(const std::vector<std::vector<double>> &matrix) const {
  double trace = 0.0;
  for (size_t i = 0; i < matrix.size(); ++i) {
    if (i < matrix[i].size()) {
      trace += matrix[i][i];
    }
  }
  return trace;
}

std::vector<std::vector<double>>
Algorithm3::matrixAdd(const std::vector<std::vector<double>> &A,
                      const std::vector<std::vector<double>> &B) const {
  std::vector<std::vector<double>> result = A;
  for (size_t i = 0; i < result.size() && i < B.size(); ++i) {
    for (size_t j = 0; j < result[i].size() && j < B[i].size(); ++j) {
      result[i][j] += B[i][j];
    }
  }
  return result;
}

std::vector<std::vector<double>>
Algorithm3::matrixMultiply(const std::vector<std::vector<double>> &A,
                           const std::vector<std::vector<double>> &B) const {
  if (A.empty() || B.empty() || A[0].empty() || B[0].empty()) {
    return {};
  }

  int rows_A = A.size();
  int cols_A = A[0].size();
  int cols_B = B[0].size();

  std::vector<std::vector<double>> result(rows_A,
                                          std::vector<double>(cols_B, 0.0));

  for (int i = 0; i < rows_A; ++i) {
    for (int j = 0; j < cols_B; ++j) {
      for (int k = 0; k < cols_A; ++k) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  return result;
}
#include "algorithm4.h"
#include "../kernels/parametric_family_kernel.h"
#include "../utils/test_data.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

Algorithm4::Algorithm4()
    : AlgorithmBase(std::make_shared<ParametricFamilyKernel>()) {
  m_familyKernel = std::make_shared<ParametricFamilyKernel>(1.0);
  // Set default parameters
  m_tolerance = 1e-8;
  m_theta = 0.5;
  m_tau = 3.0;
}

AlgorithmResult Algorithm4::solve(const TestCase &test_case) {
  auto start_time = std::chrono::high_resolution_clock::now();

  AlgorithmResult result;
  result.converged = false;
  result.iterations = 0;

  try {
    // Initialize the algorithm
    CurrentIterate current = initialize(test_case);

    // Main algorithm loop: while nμ ≥ ε do
    while (!checkConvergence(current.mu, test_case.tolerance,
                             test_case.problem_size)) {
      // μ ← (1–θ)·μ
      current.mu = (1.0 - test_case.theta) * current.mu;

      // Inner loop: while Ψ(X,Z;μ) > τ do
      while (checkProximity(computeProximityMeasure(current), test_case.tau)) {
        // Compute kernel gradient -ψ′(V)
        std::vector<std::vector<double>> V = computeScalingMatrix(current);
        std::vector<std::vector<double>> rhs_psi = computeKernelGradient(V);

        // Solve the scaled NT-system for search direction
        SearchDirection direction =
            solveLinearSystem(test_case, current, rhs_psi);

        // Determine step-size α ∈ (0,1)
        double alpha = computeStepSize(current, direction);

        // Update (X,y,Z) ← (X,y,Z) + α (ΔX,Δy,ΔZ)
        updateIterate(current, direction, alpha);

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
    result.dual_solution_s = current.Z; // Note: using Z instead of S for CQSDO
    result.primal_objective =
        computeObjectiveValue(current.X, test_case.C_matrix);
    result.duality_gap = computeDualityGap(current.X, current.Z);

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

std::string Algorithm4::getName() const { return "Algorithm 4: PD for CQSDO"; }

std::string Algorithm4::getDescription() const {
  return "Generic Primal-Dual Algorithm for CQSDO using Parametric Family "
         "kernel function ψ(t) = (t²-1)/2 - (t^(-q)/q²-q+1)e^(1/(t-1)) + "
         "(1-q)/(q²-q+1), t>0, q≥1. "
         "Implements the algorithm from the fourth research paper.";
}

std::vector<TestCase> Algorithm4::getTestCases() const {
  return {TestDataProvider::getCQSDOProblem1TestCase(),
          TestDataProvider::getCQSDOProblem2TestCase()};
}

Algorithm4::CurrentIterate
Algorithm4::initialize(const TestCase &test_case) const {
  CurrentIterate current;

  // Initialize X, y, Z from test case
  current.X = test_case.initial_X;
  current.y = test_case.initial_y;
  current.Z = test_case.initial_S; // Use S as Z for CQSDO

  // Initialize μ₀ = 1 as specified in the algorithm
  current.mu = 1.0;

  return current;
}

double
Algorithm4::computeProximityMeasure(const CurrentIterate &current) const {
  // Compute Ψ(X,Z;μ) = trace(ψ(V)) where V is the scaling matrix
  std::vector<std::vector<double>> V = computeScalingMatrix(current);

  double proximity = 0.0;
  for (size_t i = 0; i < V.size(); ++i) {
    for (size_t j = 0; j < V[i].size(); ++j) {
      if (i == j) { // Only diagonal elements for trace
        proximity += m_familyKernel->psi(V[i][j]);
      }
    }
  }

  return proximity;
}

std::vector<std::vector<double>>
Algorithm4::computeScalingMatrix(const CurrentIterate &current) const {
  // Simplified implementation: V = (μ⁻¹·X·Z)¹ᐟ²
  // In practice, this would use more sophisticated scaling
  int n = current.X.size();
  std::vector<std::vector<double>> V(n, std::vector<double>(n, 0.0));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += current.X[i][k] * current.Z[k][j];
      }
      V[i][j] = std::sqrt(sum / current.mu);
    }
  }

  return V;
}

std::vector<std::vector<double>> Algorithm4::computeKernelGradient(
    const std::vector<std::vector<double>> &V) const {
  int n = V.size();
  std::vector<std::vector<double>> gradient(n, std::vector<double>(n, 0.0));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      gradient[i][j] = -m_familyKernel->psi_prime(V[i][j]);
    }
  }

  return gradient;
}

Algorithm4::SearchDirection Algorithm4::solveLinearSystem(
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
  direction.delta_Z =
      std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

  // Simplified implementation: solve ΔX + ΔZ = −ψ′(V)
  // In practice, this would solve the complete scaled NT-system
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      // For simplicity, distribute rhs_psi equally between ΔX and ΔZ
      double delta_val = rhs_psi[i][j] / 2.0;
      direction.delta_X[i][j] = delta_val;
      direction.delta_Z[i][j] = delta_val;
    }
  }

  return direction;
}

double Algorithm4::computeStepSize(const CurrentIterate &current,
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

  // Check for Z + αΔZ > 0
  for (size_t i = 0; i < current.Z.size(); ++i) {
    for (size_t j = 0; j < current.Z[i].size(); ++j) {
      if (direction.delta_Z[i][j] < 0) {
        alpha =
            std::min(alpha, -0.9 * current.Z[i][j] / direction.delta_Z[i][j]);
      }
    }
  }

  return std::min(alpha, 0.95); // Conservative step size
}

void Algorithm4::updateIterate(CurrentIterate &current,
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

  // Update Z ← Z + α ΔZ
  for (size_t i = 0; i < current.Z.size(); ++i) {
    for (size_t j = 0; j < current.Z[i].size(); ++j) {
      current.Z[i][j] += alpha * direction.delta_Z[i][j];
    }
  }
}

bool Algorithm4::checkConvergence(double mu, double epsilon, int n) const {
  return n * mu >= epsilon;
}

bool Algorithm4::checkProximity(double proximity, double tau) const {
  return proximity > tau;
}

double Algorithm4::computeObjectiveValue(
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
Algorithm4::computeDualityGap(const std::vector<std::vector<double>> &X,
                              const std::vector<std::vector<double>> &Z) const {
  // Compute trace(X·Z)
  double gap = 0.0;
  for (size_t i = 0; i < X.size(); ++i) {
    for (size_t j = 0; j < X[i].size(); ++j) {
      gap += X[i][j] * Z[i][j];
    }
  }
  return gap;
}

bool Algorithm4::isStrictlyFeasible(
    const std::vector<std::vector<double>> &X,
    const std::vector<std::vector<double>> &Z) const {
  // Check if matrices are positive definite (simplified)
  for (const auto &row : X) {
    for (double val : row) {
      if (val <= 0)
        return false;
    }
  }
  for (const auto &row : Z) {
    for (double val : row) {
      if (val <= 0)
        return false;
    }
  }
  return true;
}

double
Algorithm4::matrixTrace(const std::vector<std::vector<double>> &matrix) const {
  double trace = 0.0;
  for (size_t i = 0; i < matrix.size(); ++i) {
    if (i < matrix[i].size()) {
      trace += matrix[i][i];
    }
  }
  return trace;
}

std::vector<std::vector<double>>
Algorithm4::matrixAdd(const std::vector<std::vector<double>> &A,
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
Algorithm4::matrixMultiply(const std::vector<std::vector<double>> &A,
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
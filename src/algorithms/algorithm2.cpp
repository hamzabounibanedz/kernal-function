#include "algorithm2.h"
#include "../kernels/exponential_parametric_kernel.h"
#include "../utils/test_data.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

Algorithm2::Algorithm2()
    : AlgorithmBase(std::make_shared<ExponentialParametricKernel>()) {
  m_expKernel = std::make_shared<ExponentialParametricKernel>(1.0);
  // Set default parameters
  m_tolerance = 1e-8;
  m_theta = 0.4;
  m_tau = 1.0;
}

AlgorithmResult Algorithm2::solve(const TestCase &test_case) {
  auto start_time = std::chrono::high_resolution_clock::now();

  AlgorithmResult result;
  result.converged = false;
  result.iterations = 0;

  try {
    // Initialize the algorithm
    CurrentIterate current = initialize(test_case);

    // Main algorithm loop: while n·μ > ε do
    while (!checkConvergence(current.mu, test_case.tolerance,
                             test_case.problem_size)) {
      // μ ← (1–θ)·μ
      current.mu = (1.0 - test_case.theta) * current.mu;

      // Inner loop: while Ψ(X,S,μ) > τ do
      while (checkProximity(computeProximityMeasure(current), test_case.tau)) {
        // Compute scaling matrix V = (μ⁻¹·D⁻¹XSD)¹ᐟ²
        std::vector<std::vector<double>> V = computeScalingMatrix(current);

        // Compute kernel gradient -ψ′(V)
        std::vector<std::vector<double>> rhs_psi = computeKernelGradient(V);

        // Solve linear system for search direction
        SearchDirection direction =
            solveLinearSystem(test_case, current, V, rhs_psi);

        // Choose step-size α∈(0,1]
        double alpha = computeStepSize(current, direction);

        // Update iterate
        updateIterate(current, direction, alpha);

        // Update V ← (μ⁻¹·D⁻¹XSD)¹ᐟ²
        V = computeScalingMatrix(current);

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

std::string Algorithm2::getName() const {
  return "Algorithm 2: PD IPM for SDO";
}

std::string Algorithm2::getDescription() const {
  return "Generic Primal-Dual Interior-Point Method for Semidefinite "
         "Optimization "
         "using Exponential-Parametric kernel function ψ(t) = (t²-1)/2 - ∫₁ᵗ "
         "((e-1)/(e^x-1))^p dx, p≥1. "
         "Implements Algorithm 1 from the second research paper.";
}

std::vector<TestCase> Algorithm2::getTestCases() const {
  return {TestDataProvider::getExponentialParametricTestCase()};
}

Algorithm2::CurrentIterate
Algorithm2::initialize(const TestCase &test_case) const {
  CurrentIterate current;

  // Initialize X, y, S from test case with positive values
  current.X = test_case.initial_X;
  current.y = test_case.initial_y;
  current.S = test_case.initial_S;

  // Ensure positive diagonal elements for interior-point methods
  for (size_t i = 0; i < current.X.size(); ++i) {
    for (size_t j = 0; j < current.X[i].size(); ++j) {
      if (i == j) {
        // Ensure diagonal elements are positive and at least 1.0
        current.X[i][j] = std::max(1.0, std::abs(current.X[i][j]));
        current.S[i][j] = std::max(1.0, std::abs(current.S[i][j]));
      }
    }
  }

  // Initialize μ₀ = trace(X·S) / n
  double trace = 0.0;
  for (size_t i = 0; i < current.X.size(); ++i) {
    for (size_t j = 0; j < current.X[i].size(); ++j) {
      trace += current.X[i][j] * current.S[i][j];
    }
  }
  current.mu = trace / current.X.size();

  // Ensure μ is positive and reasonable
  if (current.mu <= 0.0) {
    current.mu = 1.0;
  }

  return current;
}

double
Algorithm2::computeProximityMeasure(const CurrentIterate &current) const {
  // Compute Ψ(X,S,μ) = trace(ψ(V)) where V = (μ⁻¹·D⁻¹XSD)¹ᐟ²
  std::vector<std::vector<double>> V = computeScalingMatrix(current);

  double proximity = 0.0;
  for (size_t i = 0; i < V.size(); ++i) {
    for (size_t j = 0; j < V[i].size(); ++j) {
      if (i == j) { // Only diagonal elements for trace
        proximity += m_expKernel->psi(V[i][j]);
      }
    }
  }

  return proximity;
}

std::vector<std::vector<double>>
Algorithm2::computeScalingMatrix(const CurrentIterate &current) const {
  // Compute V = (μ⁻¹·X·S)¹ᐟ²
  int n = static_cast<int>(current.X.size());
  std::vector<std::vector<double>> V(n, std::vector<double>(n, 0.0));

  // Ensure μ is positive
  if (current.mu <= 0.0) {
    throw std::runtime_error(
        "Invalid barrier parameter: μ=" + std::to_string(current.mu) +
        ". Value must be positive for interior-point methods.");
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += current.X[i][k] * current.S[k][j];
      }

      // Ensure the product is positive for sqrt
      if (sum <= 0.0) {
        // For diagonal elements, use direct product
        if (i == j) {
          double diagonal_product = current.X[i][i] * current.S[i][i];
          if (diagonal_product <= 0.0) {
            throw std::runtime_error(
                "Invalid diagonal product: X[" + std::to_string(i) + "][" +
                std::to_string(i) + "] * S[" + std::to_string(i) + "][" +
                std::to_string(i) + "] = " + std::to_string(diagonal_product) +
                ". Values must be positive for interior-point methods.");
          }
          V[i][j] = std::sqrt(diagonal_product / current.mu);
        } else {
          V[i][j] = 0.0; // Off-diagonal elements set to 0 for simplicity
        }
      } else {
        V[i][j] = std::sqrt(sum / current.mu);
      }

      // Ensure V[i][j] is positive (should be, but double-check)
      if (V[i][j] <= 0.0) {
        throw std::runtime_error(
            "Invalid scaling: V[" + std::to_string(i) + "][" +
            std::to_string(j) + "]=" + std::to_string(V[i][j]) +
            ". Scaling matrix must have positive elements.");
      }
    }
  }

  return V;
}

std::vector<std::vector<double>> Algorithm2::computeKernelGradient(
    const std::vector<std::vector<double>> &V) const {
  int n = static_cast<int>(V.size());
  std::vector<std::vector<double>> gradient(n, std::vector<double>(n, 0.0));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      gradient[i][j] = -m_expKernel->psi_prime(V[i][j]);
    }
  }

  return gradient;
}

Algorithm2::SearchDirection Algorithm2::solveLinearSystem(
    const TestCase &test_case, const CurrentIterate &current,
    const std::vector<std::vector<double>> &V,
    const std::vector<std::vector<double>> &rhs_psi) const {
  (void)current; // Suppress unused parameter warning
  (void)V;       // Suppress unused parameter warning
  SearchDirection direction;
  int n = test_case.problem_size;
  int m = static_cast<int>(test_case.A_matrices.size());

  // Initialize search directions
  direction.delta_X =
      std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
  direction.delta_y.resize(m, 0.0);
  direction.delta_S =
      std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

  // Simplified implementation: solve ΔX + ΔS = –ψ′(V)
  // In practice, this would solve the complete system
  // Āᵢ•ΔX=0, ∑ᵢΔyᵢĀᵢ+ΔS=0, ΔX+ΔS = –ψ′(V)
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

double Algorithm2::computeStepSize(const CurrentIterate &current,
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

void Algorithm2::updateIterate(CurrentIterate &current,
                               const SearchDirection &direction,
                               double alpha) const {
  // Update X ← X + α·ΔX
  for (size_t i = 0; i < current.X.size(); ++i) {
    for (size_t j = 0; j < current.X[i].size(); ++j) {
      current.X[i][j] += alpha * direction.delta_X[i][j];
    }
  }

  // Update y ← y + α·Δy
  for (size_t i = 0; i < current.y.size(); ++i) {
    current.y[i] += alpha * direction.delta_y[i];
  }

  // Update S ← S + α·ΔS
  for (size_t i = 0; i < current.S.size(); ++i) {
    for (size_t j = 0; j < current.S[i].size(); ++j) {
      current.S[i][j] += alpha * direction.delta_S[i][j];
    }
  }
}

bool Algorithm2::checkConvergence(double mu, double epsilon, int n) const {
  return n * mu <= epsilon;
}

bool Algorithm2::checkProximity(double proximity, double tau) const {
  return proximity > tau;
}

double Algorithm2::computeObjectiveValue(
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
Algorithm2::computeDualityGap(const std::vector<std::vector<double>> &X,
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

bool Algorithm2::isStrictlyFeasible(
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
Algorithm2::matrixTrace(const std::vector<std::vector<double>> &matrix) const {
  double trace = 0.0;
  for (size_t i = 0; i < matrix.size(); ++i) {
    if (i < matrix[i].size()) {
      trace += matrix[i][i];
    }
  }
  return trace;
}

std::vector<std::vector<double>>
Algorithm2::matrixAdd(const std::vector<std::vector<double>> &A,
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
Algorithm2::matrixMultiply(const std::vector<std::vector<double>> &A,
                           const std::vector<std::vector<double>> &B) const {
  if (A.empty() || B.empty() || A[0].empty() || B[0].empty()) {
    return {};
  }

  int rows_A = static_cast<int>(A.size());
  int cols_A = static_cast<int>(A[0].size());
  int cols_B = static_cast<int>(B[0].size());

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
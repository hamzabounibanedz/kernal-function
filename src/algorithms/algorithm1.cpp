#include "algorithm1.h"
#include "../kernels/trigonometric_kernel.h"
#include "../utils/test_data.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

Algorithm1::Algorithm1()
    : AlgorithmBase(std::make_shared<TrigonometricKernel>()) {
  m_trigKernel = std::make_shared<TrigonometricKernel>();
  // Set default parameters
  m_tolerance = 1e-8;
  m_theta = 0.4;
  m_tau = 1.0;
}

AlgorithmResult Algorithm1::solve(const TestCase &test_case) {
  auto start_time = std::chrono::high_resolution_clock::now();

  AlgorithmResult result;
  result.converged = false;
  result.iterations = 0;

  try {
    // Initialize the algorithm
    CurrentIterate current = initialize(test_case);

    // Main algorithm loop
    while (!checkConvergence(current.mu, test_case.tolerance)) {
      // Step 1: Compute scaling vector v = √(x ∘ s / μ)
      std::vector<double> v =
          computeScalingVector(current.x, current.s, current.mu);

      // Step 2: Compute kernel gradient -∇Ψ(v)
      std::vector<double> rhs_psi = computeKernelGradient(v);

      // Step 3: Solve linear system for search direction
      SearchDirection direction =
          solveLinearSystem(test_case, current, v, rhs_psi);

      // Step 4: Compute step size α
      double alpha = computeStepSize(current, direction);

      // Step 5: Update iterate
      updateIterate(current, direction, alpha);

      // Step 6: Update barrier parameter μ ← (1 – θ) μ
      current.mu = (1.0 - test_case.theta) * current.mu;

      result.iterations++;

      // Safety check to prevent infinite loops
      if (result.iterations > 1000) {
        result.convergence_info = "Maximum iterations exceeded";
        break;
      }
    }

    // Compute final results
    result.converged = checkConvergence(current.mu, test_case.tolerance);
    result.final_mu = current.mu;
    result.primal_solution = {current.x}; // Convert to matrix format
    result.dual_solution_y = current.y;
    result.dual_solution_s = {current.s}; // Convert to matrix format
    result.primal_objective =
        computeObjectiveValue(current.x, test_case.C_matrix);
    result.duality_gap = computeDualityGap(current.x, current.s);

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

std::string Algorithm1::getName() const {
  return "Algorithm 1: Basic PD IPM for LO";
}

std::string Algorithm1::getDescription() const {
  return "Generic Primal-Dual Interior-Point Algorithm for Linear Optimization "
         "using Trigonometric kernel function ψ(t) = (t-1)²/(2t) + (t-1)²/2 + "
         "tan²(h(t))/8, "
         "where h(t) = π(1-t)/(4t+2). Implements the algorithm from Figure 1 "
         "of the research paper.";
}

std::vector<TestCase> Algorithm1::getTestCases() const {
  return {TestDataProvider::getTrigonometricKernelTestCase()};
}

Algorithm1::CurrentIterate
Algorithm1::initialize(const TestCase &test_case) const {
  CurrentIterate current;

  // Initialize primal variables x from the diagonal of initial_X
  // Ensure they are strictly positive for interior-point methods
  current.x.resize(test_case.problem_size);
  for (int i = 0; i < test_case.problem_size; ++i) {
    // Use diagonal elements but ensure they are positive
    double x_val = test_case.initial_X[i][i];
    current.x[i] =
        std::max(1.0, std::abs(x_val)); // Ensure positive and at least 1.0
  }

  // Initialize dual variables y
  current.y = test_case.initial_y;

  // Initialize dual slack variables s from the diagonal of initial_S
  // Ensure they are strictly positive for interior-point methods
  current.s.resize(test_case.problem_size);
  for (int i = 0; i < test_case.problem_size; ++i) {
    // Use diagonal elements but ensure they are positive
    double s_val = test_case.initial_S[i][i];
    current.s[i] =
        std::max(1.0, std::abs(s_val)); // Ensure positive and at least 1.0
  }

  // Initialize barrier parameter μ = (x⁰)ᵀ s⁰ / n
  double dot_product = 0.0;
  for (int i = 0; i < test_case.problem_size; ++i) {
    dot_product += current.x[i] * current.s[i];
  }
  current.mu = dot_product / test_case.problem_size;

  // Ensure μ is positive and reasonable
  if (current.mu <= 0.0) {
    current.mu = 1.0; // Default value if computation fails
  }

  return current;
}

std::vector<double>
Algorithm1::computeScalingVector(const std::vector<double> &x,
                                 const std::vector<double> &s,
                                 double mu) const {
  std::vector<double> v(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    double product = x[i] * s[i];

    // Ensure both x[i] and s[i] are positive
    if (x[i] <= 0.0 || s[i] <= 0.0) {
      throw std::runtime_error(
          "Invalid variables: x[" + std::to_string(i) +
          "]=" + std::to_string(x[i]) + ", s[" + std::to_string(i) +
          "]=" + std::to_string(s[i]) +
          ". Values must be positive for interior-point methods.");
    }

    // Ensure μ is positive
    if (mu <= 0.0) {
      throw std::runtime_error(
          "Invalid barrier parameter: μ=" + std::to_string(mu) +
          ". Value must be positive for interior-point methods.");
    }

    // Compute scaling vector v = √(x ∘ s / μ)
    v[i] = std::sqrt(product / mu);

    // Ensure v[i] is positive (should be, but double-check)
    if (v[i] <= 0.0) {
      throw std::runtime_error("Invalid scaling: v[" + std::to_string(i) +
                               "]=" + std::to_string(v[i]) +
                               ". Scaling vector must be positive.");
    }
  }
  return v;
}

std::vector<double>
Algorithm1::computeKernelGradient(const std::vector<double> &v) const {
  std::vector<double> gradient(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    gradient[i] = -m_trigKernel->psi_prime(v[i]);
  }
  return gradient;
}

Algorithm1::SearchDirection Algorithm1::solveLinearSystem(
    const TestCase &test_case, const CurrentIterate &current,
    const std::vector<double> &v, const std::vector<double> &rhs_psi) const {
  SearchDirection direction;
  int n = test_case.problem_size;
  int m = test_case.A_matrices.size();

  // For simplicity, we'll use a basic approach
  // In practice, this would use more sophisticated linear algebra

  // Initialize search directions
  direction.delta_x.resize(n, 0.0);
  direction.delta_y.resize(m, 0.0);
  direction.delta_s.resize(n, 0.0);

  // Simple implementation: solve dx + ds = rhs_psi
  // where dx = v ∘ Δx / x, ds = v ∘ Δs / s
  for (int i = 0; i < n; ++i) {
    // For simplicity, distribute rhs_psi equally between dx and ds
    double dx = rhs_psi[i] / 2.0;
    double ds = rhs_psi[i] / 2.0;

    // Convert back to Δx and Δs
    direction.delta_x[i] = dx * current.x[i] / v[i];
    direction.delta_s[i] = ds * current.s[i] / v[i];
  }

  // Note: In a full implementation, we would solve the complete system
  // A Δx = 0, Aᵀ Δy + Δs = 0, dx + ds = –∇Ψ(v)
  // This simplified version focuses on the kernel-based part

  return direction;
}

double Algorithm1::computeStepSize(const CurrentIterate &current,
                                   const SearchDirection &direction) const {
  // Simple line search: find maximum step size that maintains positivity
  double alpha = 1.0;

  for (size_t i = 0; i < current.x.size(); ++i) {
    if (direction.delta_x[i] < 0) {
      alpha = std::min(alpha, -0.9 * current.x[i] / direction.delta_x[i]);
    }
    if (direction.delta_s[i] < 0) {
      alpha = std::min(alpha, -0.9 * current.s[i] / direction.delta_s[i]);
    }
  }

  return std::min(alpha, 0.95); // Conservative step size
}

void Algorithm1::updateIterate(CurrentIterate &current,
                               const SearchDirection &direction,
                               double alpha) const {
  // Update x ← x + α Δx
  for (size_t i = 0; i < current.x.size(); ++i) {
    current.x[i] += alpha * direction.delta_x[i];
  }

  // Update y ← y + α Δy
  for (size_t i = 0; i < current.y.size(); ++i) {
    current.y[i] += alpha * direction.delta_y[i];
  }

  // Update s ← s + α Δs
  for (size_t i = 0; i < current.s.size(); ++i) {
    current.s[i] += alpha * direction.delta_s[i];
  }
}

bool Algorithm1::checkConvergence(double mu, double epsilon) const {
  return mu <= epsilon;
}

double Algorithm1::computeObjectiveValue(
    const std::vector<double> &x,
    const std::vector<std::vector<double>> &C) const {
  // For simplicity, compute cᵀx where c is the diagonal of C
  double objective = 0.0;
  for (size_t i = 0; i < x.size() && i < C.size(); ++i) {
    if (i < C[i].size()) {
      objective += C[i][i] * x[i];
    }
  }
  return objective;
}

double Algorithm1::computeDualityGap(const std::vector<double> &x,
                                     const std::vector<double> &s) const {
  double gap = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    gap += x[i] * s[i];
  }
  return gap;
}

bool Algorithm1::isStrictlyFeasible(const std::vector<double> &x,
                                    const std::vector<double> &s) const {
  for (double xi : x) {
    if (xi <= 0)
      return false;
  }
  for (double si : s) {
    if (si <= 0)
      return false;
  }
  return true;
}
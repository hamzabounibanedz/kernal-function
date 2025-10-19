#include "algorithm_base.h"
#include <chrono>
#include <cmath>

AlgorithmBase::AlgorithmBase(std::shared_ptr<KernelBase> kernel)
    : m_kernel(kernel), m_tolerance(1e-8), m_theta(0.5), m_tau(3.0) {}

void AlgorithmBase::setParameters(double tolerance, double theta, double tau) {
  m_tolerance = tolerance;
  m_theta = theta;
  m_tau = tau;
}

bool AlgorithmBase::validateInput(const TestCase &test_case) const {
  // Basic validation
  if (test_case.problem_size <= 0)
    return false;
  if (test_case.A_matrices.empty())
    return false;
  if (test_case.b_vector.empty())
    return false;
  if (test_case.C_matrix.empty())
    return false;

  // Check dimensions consistency
  int n = test_case.problem_size;
  int m = static_cast<int>(test_case.A_matrices.size());

  // Check A matrices
  for (const auto &A : test_case.A_matrices) {
    if (A.size() != static_cast<size_t>(n) ||
        A[0].size() != static_cast<size_t>(n))
      return false;
  }

  // Check b vector
  if (test_case.b_vector.size() != static_cast<size_t>(m))
    return false;

  // Check C matrix
  if (test_case.C_matrix.size() != static_cast<size_t>(n) ||
      test_case.C_matrix[0].size() != static_cast<size_t>(n))
    return false;

  return true;
}

double AlgorithmBase::computeInitialMu(const TestCase &test_case) const {
  (void)test_case; // Suppress unused parameter warning
  // Simple initialization: μ₀ = 1.0 (as specified in most algorithms)
  return 1.0;
}

bool AlgorithmBase::checkConvergence(double mu, double duality_gap) const {
  return (mu <= m_tolerance) && (std::abs(duality_gap) <= m_tolerance);
}

double AlgorithmBase::matrixTrace(
    const std::vector<std::vector<double>> &matrix) const {
  double trace = 0.0;
  size_t size = std::min(matrix.size(), matrix[0].size());
  for (size_t i = 0; i < size; ++i) {
    trace += matrix[i][i];
  }
  return trace;
}

std::vector<std::vector<double>>
AlgorithmBase::matrixMultiply(const std::vector<std::vector<double>> &A,
                              const std::vector<std::vector<double>> &B) const {

  size_t rows_A = A.size();
  size_t cols_A = A[0].size();
  size_t cols_B = B[0].size();

  std::vector<std::vector<double>> result(rows_A,
                                          std::vector<double>(cols_B, 0.0));

  for (size_t i = 0; i < rows_A; ++i) {
    for (size_t j = 0; j < cols_B; ++j) {
      for (size_t k = 0; k < cols_A; ++k) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  return result;
}

std::vector<std::vector<double>>
AlgorithmBase::matrixAdd(const std::vector<std::vector<double>> &A,
                         const std::vector<std::vector<double>> &B) const {

  size_t rows = A.size();
  size_t cols = A[0].size();

  std::vector<std::vector<double>> result(rows, std::vector<double>(cols));

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      result[i][j] = A[i][j] + B[i][j];
    }
  }

  return result;
}

std::vector<std::vector<double>>
AlgorithmBase::matrixSubtract(const std::vector<std::vector<double>> &A,
                              const std::vector<std::vector<double>> &B) const {

  size_t rows = A.size();
  size_t cols = A[0].size();

  std::vector<std::vector<double>> result(rows, std::vector<double>(cols));

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      result[i][j] = A[i][j] - B[i][j];
    }
  }

  return result;
}

bool AlgorithmBase::isPositiveDefinite(
    const std::vector<std::vector<double>> &matrix) const {
  // Simple check: all diagonal elements positive (sufficient for many cases)
  for (size_t i = 0; i < matrix.size(); ++i) {
    if (matrix[i][i] <= 0.0)
      return false;
  }
  return true;
}
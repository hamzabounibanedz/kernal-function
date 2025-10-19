#include "utils/test_data.h"
#include <iomanip>
#include <iostream>


void testTestCaseData() {
  std::cout << "=== TEST CASE VERIFICATION ===\n\n";

  // Test each test case
  auto trigCase = TestDataProvider::getTrigonometricKernelTestCase();
  auto expCase = TestDataProvider::getExponentialParametricTestCase();
  auto logCase = TestDataProvider::getParameterizedLogTestCase();
  auto cqsdoCase = TestDataProvider::getCQSDOProblem1TestCase();
  auto bachirCase = TestDataProvider::getBachirTestCase();

  std::vector<std::pair<std::string, TestCase>> testCases = {
      {"Trigonometric Kernel", trigCase},
      {"Exponential Parametric", expCase},
      {"Parameterized Log", logCase},
      {"CQSDO Problem 1", cqsdoCase},
      {"Bachir Kernel", bachirCase}};

  for (const auto &[name, testCase] : testCases) {
    std::cout << "Test Case: " << name << "\n";
    std::cout << "  Name: " << testCase.name << "\n";
    std::cout << "  Paper: " << testCase.paper_reference << "\n";
    std::cout << "  Problem Size: n=" << testCase.problem_size
              << ", m=" << testCase.A_matrices.size() << "\n";
    std::cout << "  Parameters: ε=" << std::scientific << testCase.tolerance
              << ", θ=" << std::fixed << testCase.theta
              << ", τ=" << testCase.tau << "\n";
    std::cout << "  Expected Optimal Value: " << testCase.expected_optimal_value
              << "\n";

    // Check initial points
    std::cout << "  Initial X (diagonal): [";
    for (int i = 0; i < testCase.problem_size; ++i) {
      std::cout << testCase.initial_X[i][i];
      if (i < testCase.problem_size - 1)
        std::cout << ", ";
    }
    std::cout << "]\n";

    std::cout << "  Initial y: [";
    for (size_t i = 0; i < testCase.initial_y.size(); ++i) {
      std::cout << testCase.initial_y[i];
      if (i < testCase.initial_y.size() - 1)
        std::cout << ", ";
    }
    std::cout << "]\n";

    std::cout << "  Initial S (diagonal): [";
    for (int i = 0; i < testCase.problem_size; ++i) {
      std::cout << testCase.initial_S[i][i];
      if (i < testCase.problem_size - 1)
        std::cout << ", ";
    }
    std::cout << "]\n";

    // Check constraint matrices
    std::cout << "  Constraint Matrices A:\n";
    for (size_t i = 0; i < testCase.A_matrices.size(); ++i) {
      std::cout << "    A" << (i + 1) << " (first row): [";
      for (size_t j = 0; j < testCase.A_matrices[i][0].size(); ++j) {
        std::cout << testCase.A_matrices[i][0][j];
        if (j < testCase.A_matrices[i][0].size() - 1)
          std::cout << ", ";
      }
      std::cout << "]\n";
    }

    // Check b vector
    std::cout << "  b vector: [";
    for (size_t i = 0; i < testCase.b_vector.size(); ++i) {
      std::cout << testCase.b_vector[i];
      if (i < testCase.b_vector.size() - 1)
        std::cout << ", ";
    }
    std::cout << "]\n";

    // Check C matrix
    std::cout << "  C matrix (first row): [";
    for (size_t j = 0; j < testCase.C_matrix[0].size(); ++j) {
      std::cout << testCase.C_matrix[0][j];
      if (j < testCase.C_matrix[0].size() - 1)
        std::cout << ", ";
    }
    std::cout << "]\n";

    std::cout << "\n";
  }
}

void testKernelEvaluation() {
  std::cout << "=== KERNEL FUNCTION EVALUATION ===\n\n";

  // Test kernel functions at t = 1.5
  double t = 1.5;
  std::cout << "Testing kernel functions at t = " << t << ":\n\n";

  // Create kernel instances
  auto trig_kernel = std::make_shared<TrigonometricKernel>();
  auto exp_kernel = std::make_shared<ExponentialParametricKernel>(1.0);
  auto log_kernel = std::make_shared<ParameterizedLogKernel>(1.0);
  auto family_kernel = std::make_shared<ParametricFamilyKernel>(1.0);
  auto bachir_kernel = std::make_shared<BachirKernel>(0.5);

  std::vector<std::pair<std::string, std::shared_ptr<KernelBase>>> kernels = {
      {"Trigonometric", trig_kernel},
      {"Exponential-Parametric", exp_kernel},
      {"Parameterized Log", log_kernel},
      {"Parametric Family", family_kernel},
      {"Bachir φₘ(t)", bachir_kernel}};

  for (const auto &[name, kernel] : kernels) {
    try {
      auto [psi_val, psi_prime_val, psi_double_prime_val, _] =
          kernel->evaluateAll(t);

      std::cout << name << " Kernel:\n";
      std::cout << "  ψ(" << t << ") = " << std::setprecision(6) << psi_val
                << "\n";
      std::cout << "  ψ'(" << t << ") = " << std::setprecision(6)
                << psi_prime_val << "\n";
      std::cout << "  ψ''(" << t << ") = " << std::setprecision(6)
                << psi_double_prime_val << "\n";
      std::cout << "  Formula: " << kernel->getFormula() << "\n\n";

    } catch (const std::exception &e) {
      std::cout << name << " Kernel: ERROR - " << e.what() << "\n\n";
    }
  }
}

int main() {
  std::cout << "Kernel Function Test Verification\n";
  std::cout << "==================================\n\n";

  testTestCaseData();
  testKernelEvaluation();

  std::cout << "Verification complete!\n";
  return 0;
}
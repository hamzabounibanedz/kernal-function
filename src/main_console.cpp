#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "algorithms/algorithm1.h"
#include "algorithms/algorithm2.h"
#include "algorithms/algorithm3.h"
#include "algorithms/algorithm4.h"
#include "algorithms/algorithm5.h"
#include "kernels/bachir_kernel.h"
#include "kernels/exponential_parametric_kernel.h"
#include "kernels/parameterized_log_kernel.h"
#include "kernels/parametric_family_kernel.h"
#include "kernels/trigonometric_kernel.h"
#include "kernels/derbal20_param_log_kernel.h"
#include "kernels/ijnna23_fractional_kernel.h"
#include "engine/ipm_engine.h"
#include "utils/test_data.h"

void printHeader() {
  std::cout << "==========================================================\n";
  std::cout
      << "    Kernel Function Comparison Tool - Research Implementation\n";
  std::cout << "==========================================================\n\n";
}

void printAlgorithmHeader(const std::string &algorithmName) {
  std::cout << "\n" << std::string(50, '=') << "\n";
  std::cout << "  " << algorithmName << "\n";
  std::cout << std::string(50, '=') << "\n";
}

void printTestCaseInfo(const TestCase &testCase) {
  std::cout << "\nTest Case: " << testCase.name << "\n";
  std::cout << "Paper Reference: " << testCase.paper_reference << "\n";
  std::cout << "Algorithm Type: " << testCase.algorithm_type << "\n";
  std::cout << "Problem size: n=" << testCase.problem_size
            << ", m=" << testCase.A_matrices.size() << "\n";
  std::cout << "Parameters: ε=" << std::scientific << testCase.tolerance
            << ", θ=" << std::fixed << testCase.theta << ", τ=" << testCase.tau
            << "\n";
  std::cout << "Expected optimal value: " << testCase.expected_optimal_value
            << "\n";
  std::cout << "Notes: " << testCase.notes << "\n";

  // Display initial points
  std::cout << "\nInitial Points:\n";
  std::cout << "X₀ (diagonal): [";
  for (int i = 0; i < testCase.problem_size; ++i) {
    std::cout << testCase.initial_X[i][i];
    if (i < testCase.problem_size - 1)
      std::cout << ", ";
  }
  std::cout << "]\n";

  std::cout << "y₀: [";
  for (size_t i = 0; i < testCase.initial_y.size(); ++i) {
    std::cout << testCase.initial_y[i];
    if (i < testCase.initial_y.size() - 1)
      std::cout << ", ";
  }
  std::cout << "]\n";

  std::cout << "S₀ (diagonal): [";
  for (int i = 0; i < testCase.problem_size; ++i) {
    std::cout << testCase.initial_S[i][i];
    if (i < testCase.problem_size - 1)
      std::cout << ", ";
  }
  std::cout << "]\n";
}

void printAlgorithmResult(const AlgorithmResult &result,
                          const std::string &algorithmName) {
  std::cout << "\n" << algorithmName << " completed\n";
  std::cout << "   → " << result.iterations << " iterations, " << std::fixed
            << std::setprecision(2) << result.execution_time_ms
            << " ms, objective: " << std::setprecision(6)
            << result.primal_objective << "\n";

  if (!result.converged) {
    std::cout << "   → WARNING: Algorithm did not converge!\n";
    std::cout << "   → Convergence info: " << result.convergence_info << "\n";
  }

  std::cout << "   → Final μ: " << std::scientific << result.final_mu << "\n";
  std::cout << "   → Duality gap: " << std::scientific << result.duality_gap
            << "\n";
}

void testKernelFunctions() {
  std::cout << "PART 1: Kernel Function Evaluation\n";
  std::cout << "==================================\n\n";

  // Test all kernel functions at specific points
  std::vector<double> test_points = {0.5, 1.0, 1.5, 2.0, 3.0};

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
    std::cout << "Testing " << name << " Kernel:\n";
    std::cout << "Formula: " << kernel->getFormula() << "\n\n";

    std::cout << std::setw(8) << "t" << std::setw(15) << "ψ(t)" << std::setw(15)
              << "ψ'(t)" << std::setw(15) << "ψ''(t)" << std::setw(15)
              << "ψ'''(t)" << "\n";
    std::cout << std::string(70, '-') << "\n";

    for (double t : test_points) {
      try {
        auto [psi_val, psi_prime_val, psi_double_prime_val,
              psi_triple_prime_val] = kernel->evaluateAll(t);

        std::cout << std::setw(8) << std::fixed << std::setprecision(2) << t
                  << std::setw(15) << std::setprecision(6) << psi_val
                  << std::setw(15) << std::setprecision(6) << psi_prime_val
                  << std::setw(15) << std::setprecision(6)
                  << psi_double_prime_val << std::setw(15)
                  << std::setprecision(6) << psi_triple_prime_val << "\n";
      } catch (const std::exception &e) {
        std::cout << "Error at t=" << t << ": " << e.what() << "\n";
      }
    }
    std::cout << "\n";
  }
}

void executeAllAlgorithms() {
  std::cout << "\nPART 2: Algorithm Execution and Comparison\n";
  std::cout << "===========================================\n\n";

  std::cout
      << "Starting execution of all 5 algorithms with paper test cases...\n";
  std::cout << std::string(50, '=') << "\n";

  // Create algorithm instances
  auto alg1 = std::make_shared<Algorithm1>();
  auto alg2 = std::make_shared<Algorithm2>();
  auto alg3 = std::make_shared<Algorithm3>();
  auto alg4 = std::make_shared<Algorithm4>();
  auto alg5 = std::make_shared<Algorithm5>();

  // Get specific test cases for each algorithm
  auto trigTestCase = TestDataProvider::getTrigonometricKernelTestCase();
  auto expTestCase = TestDataProvider::getExponentialParametricTestCase();
  auto logTestCase = TestDataProvider::getParameterizedLogTestCase();
  auto cqsdoTestCase1 = TestDataProvider::getCQSDOProblem1TestCase();
  auto bachirTestCase = TestDataProvider::getBachirTestCase();

  std::vector<std::pair<std::string, std::shared_ptr<AlgorithmBase>>>
      algorithms = {{"Algorithm 1: Trigonometric Kernel", alg1},
                    {"Algorithm 2: Exponential Parametric", alg2},
                    {"Algorithm 3: Parameterized Log", alg3},
                    {"Algorithm 4: CQSDO Problem 1", alg4},
                    {"Algorithm 5: Bachir Kernel", alg5}};

  std::vector<TestCase> testCases = {trigTestCase, expTestCase, logTestCase,
                                     cqsdoTestCase1, bachirTestCase};

  // Execute each algorithm on its specific test case
  for (size_t i = 0; i < algorithms.size(); ++i) {
    const auto &[name, algorithm] = algorithms[i];
    const auto &testCase = testCases[i];

    printAlgorithmHeader(name);
    printTestCaseInfo(testCase);

    try {
      auto result = algorithm->solve(testCase);
      printAlgorithmResult(result, name);
    } catch (const std::exception &e) {
      std::cout << "\nERROR: " << e.what() << "\n";
    }
  }

  std::cout << "\n" << std::string(50, '=') << "\n";
  std::cout << "All algorithms completed successfully!\n";
  std::cout << "Results ready for analysis and comparison.\n";
  std::cout << std::string(50, '=') << "\n";
}

void executeAlgorithmWithMultipleThetas() {
  std::cout << "\nPART 3: Algorithm Performance with Different θ Values\n";
  std::cout << "=====================================================\n\n";

  // Test Algorithm 2 (Exponential Parametric) with different θ values
  auto alg2 = std::make_shared<Algorithm2>();
  auto testCase = TestDataProvider::getExponentialParametricTestCase();

  std::cout << "Testing Algorithm 2 with different θ values from paper:\n";
  std::cout << "θ ∈ {0.05, 0.4, 0.6, 0.95}\n\n";

  std::vector<double> thetaValues = {0.05, 0.4, 0.6, 0.95};

  for (double theta : thetaValues) {
    testCase.theta = theta;

    std::cout << "θ = " << theta << ":\n";
    printTestCaseInfo(testCase);

    try {
      auto result = alg2->solve(testCase);
      printAlgorithmResult(result, "Algorithm 2");
    } catch (const std::exception &e) {
      std::cout << "\nERROR: " << e.what() << "\n";
    }
    std::cout << "\n";
  }
}

void compareKernelFunctions() {
  std::cout << "\nPART 4: Kernel Function Comparison\n";
  std::cout << "===================================\n\n";

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

  // Compare kernel functions at t = 1.5 (typical value)
  double t = 1.5;
  std::cout << "Kernel Function Comparison at t = " << t << ":\n\n";

  std::cout << std::setw(20) << "Kernel Type" << std::setw(15) << "ψ(t)"
            << std::setw(15) << "ψ'(t)" << std::setw(15) << "ψ''(t)" << "\n";
  std::cout << std::string(70, '-') << "\n";

  for (const auto &[name, kernel] : kernels) {
    try {
      auto [psi_val, psi_prime_val, psi_double_prime_val, _] =
          kernel->evaluateAll(t);

      std::cout << std::setw(20) << name << std::setw(15)
                << std::setprecision(6) << psi_val << std::setw(15)
                << std::setprecision(6) << psi_prime_val << std::setw(15)
                << std::setprecision(6) << psi_double_prime_val << "\n";
    } catch (const std::exception &e) {
      std::cout << std::setw(20) << name << " Error: " << e.what() << "\n";
    }
  }
}

void runIndependentAlgorithmTests() {
  std::cout << "\nPART 5: Independent Algorithm Tests with Paper Examples\n";
  std::cout << "=======================================================\n\n";

  // Test each algorithm independently with its specific test case
  std::vector<std::pair<std::string, std::shared_ptr<AlgorithmBase>>>
      algorithms = {
          {"Algorithm 1: Trigonometric Kernel", std::make_shared<Algorithm1>()},
          {"Algorithm 2: Exponential Parametric",
           std::make_shared<Algorithm2>()},
          {"Algorithm 3: Parameterized Log", std::make_shared<Algorithm3>()},
          {"Algorithm 4: CQSDO Problem 1", std::make_shared<Algorithm4>()},
          {"Algorithm 5: Bachir Kernel", std::make_shared<Algorithm5>()}};

  std::vector<TestCase> testCases = {
      TestDataProvider::getTrigonometricKernelTestCase(),
      TestDataProvider::getExponentialParametricTestCase(),
      TestDataProvider::getParameterizedLogTestCase(),
      TestDataProvider::getCQSDOProblem1TestCase(),
      TestDataProvider::getBachirTestCase()};

  // Results storage for comparison
  struct TestResult {
    std::string algorithm_name;
    std::string test_case_name;
    double theta;
    int iterations;
    double execution_time_ms;
    double objective_value;
    bool converged;
    std::string status;
  };
  std::vector<TestResult> results;

  // Test each algorithm with different θ values
  std::vector<double> thetaValues = {0.1, 0.3, 0.5, 0.7, 0.9};

  for (size_t i = 0; i < algorithms.size(); ++i) {
    const auto &[name, algorithm] = algorithms[i];
    auto testCase = testCases[i];

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << name << "\n";
    std::cout << "  Test Case: " << testCase.name << "\n";
    std::cout << "  Paper: " << testCase.paper_reference << "\n";
    std::cout << std::string(60, '=') << "\n";

    for (double theta : thetaValues) {
      testCase.theta = theta;

      std::cout << "\nθ = " << theta << ":\n";

      try {
        auto result = algorithm->solve(testCase);

        TestResult testResult;
        testResult.algorithm_name = name;
        testResult.test_case_name = testCase.name;
        testResult.theta = theta;
        testResult.iterations = result.iterations;
        testResult.execution_time_ms = result.execution_time_ms;
        testResult.objective_value = result.primal_objective;
        testResult.converged = result.converged;
        testResult.status = result.convergence_info;

        results.push_back(testResult);

        std::cout << "   → " << result.iterations << " iterations, "
                  << std::fixed << std::setprecision(2)
                  << result.execution_time_ms
                  << " ms, objective: " << std::setprecision(6)
                  << result.primal_objective;

        if (!result.converged) {
          std::cout << " (DID NOT CONVERGE)";
        }
        std::cout << "\n";

      } catch (const std::exception &e) {
        std::cout << "   → ERROR: " << e.what() << "\n";
      }
    }
  }

  // Print comparison table
  std::cout << "\n" << std::string(80, '=') << "\n";
  std::cout << "COMPARISON TABLE: Iteration Counts and Execution Times\n";
  std::cout << std::string(80, '=') << "\n\n";

  std::cout << std::setw(25) << "Algorithm" << std::setw(15) << "θ"
            << std::setw(12) << "Iterations" << std::setw(15) << "Time (ms)"
            << std::setw(15) << "Objective" << std::setw(10) << "Status"
            << "\n";
  std::cout << std::string(95, '-') << "\n";

  for (const auto &result : results) {
    std::cout << std::setw(25) << result.algorithm_name.substr(0, 24)
              << std::setw(15) << std::fixed << std::setprecision(1)
              << result.theta << std::setw(12) << result.iterations
              << std::setw(15) << std::setprecision(2)
              << result.execution_time_ms << std::setw(15)
              << std::setprecision(6) << result.objective_value << std::setw(10)
              << (result.converged ? "OK" : "FAIL") << "\n";
  }

  // Find best performing algorithm for each θ
  std::cout << "\n" << std::string(60, '=') << "\n";
  std::cout << "BEST PERFORMANCE BY θ VALUE\n";
  std::cout << std::string(60, '=') << "\n\n";

  for (double theta : thetaValues) {
    std::cout << "θ = " << theta << ":\n";

    // Find minimum iterations for this theta
    int min_iterations = INT_MAX;
    std::string best_algorithm = "";

    for (const auto &result : results) {
      if (result.theta == theta && result.converged &&
          result.iterations < min_iterations) {
        min_iterations = result.iterations;
        best_algorithm = result.algorithm_name;
      }
    }

    if (min_iterations != INT_MAX) {
      std::cout << "   Best: " << best_algorithm << " (" << min_iterations
                << " iterations)\n";
    } else {
      std::cout << "   No algorithm converged\n";
    }
  }
}

int main() {
  printHeader();

  // Part 1: Test kernel functions
  testKernelFunctions();

  // Part 2: Execute all algorithms with REAL paper test cases
  std::cout << "\nPART 2: Algorithm Execution with REAL Paper Test Cases\n";
  std::cout << "=======================================================\n\n";

  // Create algorithms
  auto alg1 = std::make_shared<Algorithm1>();
  auto alg2 = std::make_shared<Algorithm2>();
  auto alg3 = std::make_shared<Algorithm3>();
  auto alg4 = std::make_shared<Algorithm4>();
  auto alg5 = std::make_shared<Algorithm5>();

  // Get REAL test cases from papers
  auto trigTestCase = TestDataProvider::getTrigonometricKernelTestCase();
  auto expTestCase = TestDataProvider::getExponentialParametricTestCase();
  auto logTestCase = TestDataProvider::getParameterizedLogTestCase();
  auto cqsdoTestCase1 = TestDataProvider::getCQSDOProblem1TestCase();
  auto bachirTestCase = TestDataProvider::getBachirTestCase();

  std::vector<std::pair<std::string, std::shared_ptr<AlgorithmBase>>>
      algorithms = {{"Algorithm 1: Trigonometric Kernel", alg1},
                    {"Algorithm 2: Exponential Parametric", alg2},
                    {"Algorithm 3: Parameterized Log", alg3},
                    {"Algorithm 4: CQSDO Problem 1", alg4},
                    {"Algorithm 5: Bachir Kernel", alg5}};

  std::vector<TestCase> testCases = {trigTestCase, expTestCase, logTestCase,
                                     cqsdoTestCase1, bachirTestCase};

  // Execute each algorithm on its REAL test case
  for (size_t i = 0; i < algorithms.size(); ++i) {
    const auto &[name, algorithm] = algorithms[i];
    const auto &testCase = testCases[i];

    printAlgorithmHeader(name);
    printTestCaseInfo(testCase);

    try {
      auto result = algorithm->solve(testCase);
      printAlgorithmResult(result, name);
    } catch (const std::exception &e) {
      std::cout << "\nERROR: " << e.what() << "\n";
    }
  }

  std::cout << "\n" << std::string(50, '=') << "\n";
  std::cout << "All algorithms completed with REAL test cases!\n";
  std::cout << "Results ready for analysis and comparison.\n";
  std::cout << std::string(50, '=') << "\n";

  // Part 3: Test with different θ values
  executeAlgorithmWithMultipleThetas();

  // Part 4: Compare kernel functions
  compareKernelFunctions();

  // Part 5: Independent algorithm tests with detailed comparison
  runIndependentAlgorithmTests();

  // Derbal20 vs IJNAA23 identical 5x5 test with IPM engine
  {
    std::cout << "\nDERBAL20 vs IJNAA23 on common 5x5 SDP (q=2.0)\n";
    TestCase tc = TestDataProvider::getParameterizedLogTestCase();
    IpmEngine engine;
    IpmSettings settings;
    settings.epsilon = tc.tolerance;
    settings.theta = tc.theta;
    settings.tau = tc.tau;
    Derbal20ParamLogKernel k1(2.0);
    Ijna23FractionalKernel k2(2.0);
    auto r1 = engine.solve(tc, k1, settings);
    auto r2 = engine.solve(tc, k2, settings);
    std::cout << "Derbal20: it=" << r1.iterations << ", obj=" << r1.primal_objective
              << ", conv=" << (r1.converged?"yes":"no") << "\n";
    std::cout << "IJNAA23: it=" << r2.iterations << ", obj=" << r2.primal_objective
              << ", conv=" << (r2.converged?"yes":"no") << "\n";
  }

  std::cout << "\n" << std::string(50, '=') << "\n";
  std::cout << "All tests completed!\n";
  std::cout << std::string(50, '=') << "\n";

  return 0;
}
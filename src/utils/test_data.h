#ifndef TEST_DATA_H
#define TEST_DATA_H

#include <string>
#include <vector>

/**
 * TestCase - Represents an optimization problem to solve
 *
 * This structure contains all the data needed to define an optimization
 * problem:
 * - The constraint matrices A and vectors b
 * - The objective function matrix C
 * - Initial starting points for the algorithm
 * - Expected solution and parameters
 *
 * Think of it as a "puzzle" that we give to our algorithms to solve.
 */
struct TestCase {
  std::string name;            // Human-readable name for the problem
  std::string paper_reference; // Reference to the source paper
  std::string
      algorithm_type; // Type of algorithm this test case is designed for
  int problem_size;   // Size of the problem (e.g., 5x5 matrices)
  std::vector<std::vector<std::vector<double>>>
      A_matrices;               // Constraint matrices for SDO problems
  std::vector<double> b_vector; // Right-hand side of constraints
  std::vector<std::vector<double>> C_matrix; // Objective function matrix
  std::vector<std::vector<double>>
      initial_X;                 // Starting point for primal variable X
  std::vector<double> initial_y; // Starting point for dual variable y
  std::vector<std::vector<double>>
      initial_S;                 // Starting point for dual variable S
  double expected_optimal_value; // What the optimal solution should be
  double tolerance;              // How close is "close enough" for convergence
  double theta;                  // Step size parameter θ
  double tau;                    // Barrier parameter update factor τ
  std::vector<double> theta_values; // Multiple θ values to test
  std::string notes;                // Additional notes about the test case
};

/**
 * TestDataProvider - Provides test problems for algorithm comparison
 *
 * This class contains a collection of carefully chosen test problems
 * extracted directly from research papers. Each test case is designed
 * to highlight different aspects of algorithm behavior and comes with
 * known optimal solutions for validation.
 */
class TestDataProvider {
public:
  // === PAPER-SPECIFIC TEST CASES ===

  /**
   * Test Case 1: Generic Primal-Dual Interior-Point Algorithm for LO
   * From the first paper - Trigonometric Kernel Function
   */
  static TestCase getTrigonometricKernelTestCase();

  /**
   * Test Case 2: Generic Primal-Dual IPM for SDO - Exponential Parametric
   * From the second paper - Example 7.1
   */
  static TestCase getExponentialParametricTestCase();

  /**
   * Test Case 3: Generic Primal-Dual IPM for SDO - Parameterized Log
   * From the third paper - Section 4
   */
  static TestCase getParameterizedLogTestCase();

  /**
   * Test Case 4: Generic Primal-Dual Algorithm for CQSDO - Problem 1
   * From the fourth paper - Problem 1 (n=5, m=3, Q=0)
   */
  static TestCase getCQSDOProblem1TestCase();

  /**
   * Test Case 5: Generic Primal-Dual Algorithm for CQSDO - Problem 2
   * From the fourth paper - Problem 2 (n=4, m=3, Q=I)
   */
  static TestCase getCQSDOProblem2TestCase();

  // Paper pattern tests (Examples 4.1 and 4.2)
  static TestCase getTouil22Problem1(int n, int m);
  static TestCase getTouil22Problem2(int n, int m);

  // IJNAA 2023 Section 4 explicit 5x5, m=3 test (positive matrices version)
  static TestCase getIJNAA2023Section4TestCase();

  /**
   * Test Case 6: Generic Primal-Dual Interior-Point Algorithm - Bachir
   * From the fifth paper - Derbal & Kebbiche (2020)
   */
  static TestCase getBachirTestCase();

  /**
   * Test Case 7: Generic Primal-Dual Interior-Point Algorithm - φₘ(t) Family
   * From the sixth paper - Bai et al. kernel family
   */
  static TestCase getBaiKernelFamilyTestCase();

  // === LEGACY TEST CASES (for backward compatibility) ===

  /**
   * Get test case for Algorithm 1 (Exponential-Parametric)
   * @deprecated Use getExponentialParametricTestCase() instead
   */
  static TestCase getAlgorithm2TestCase();

  /**
   * Get test case for Algorithm 2 (Parameterized Log)
   * @deprecated Use getParameterizedLogTestCase() instead
   */
  static TestCase getAlgorithm3TestCase();

  /**
   * Get test case for Algorithm 3 (Parametric Family) - Problem 1
   * @deprecated Use getCQSDOProblem1TestCase() instead
   */
  static TestCase getAlgorithm4Problem1();

  /**
   * Get test case for Algorithm 3 (Parametric Family) - Problem 2
   * @deprecated Use getCQSDOProblem2TestCase() instead
   */
  static TestCase getAlgorithm4Problem2();

  /**
   * Get test case for Algorithm 4 (Derbal & Kebbiche)
   * @deprecated Use getBachirTestCase() instead
   */
  static TestCase getAlgorithm5TestCase();

  /**
   * Get test case for Algorithm 5 (Trigonometric)
   * @deprecated Use getTrigonometricKernelTestCase() instead
   */
  static TestCase getAlgorithm6TestCase();

  // === SPECIALIZED TEST CASES ===

  /**
   * Get test case where Bachir algorithm performs exceptionally well
   */
  static TestCase getBachirOptimizedTestCase();

  /**
   * Get a high-dimensional test case (6x6) where Bachir dominates
   */
  static TestCase getHighDimensionalBachirTestCase();

  /**
   * Get a random linear programming test case (50×100)
   */
  static TestCase getRandomLPTestCase();

  /**
   * Get all available test cases in one convenient collection
   */
  static std::vector<TestCase> getAllTestCases();

  /**
   * Get test cases by algorithm type
   */
  static std::vector<TestCase>
  getTestCasesByAlgorithm(const std::string &algorithm_type);

  /**
   * Get test cases by paper reference
   */
  static std::vector<TestCase>
  getTestCasesByPaper(const std::string &paper_reference);

private:
  // === COMMON DATA GENERATORS ===

  /**
   * Get the common 5x5 constraint matrices used across multiple algorithms
   * These matrices come from the research literature and are used
   * in several different test cases.
   */
  static std::vector<std::vector<std::vector<double>>> getCommon5x5AMatrices();

  /**
   * Get the common right-hand side vector for 5x5 problems
   */
  static std::vector<double> getCommon5x5BVector();

  /**
   * Get the common objective function matrix for 5x5 problems
   */
  static std::vector<std::vector<double>> getCommon5x5CMatrix();

  // === PAPER PATTERN TESTS ===
  static TestCase getWu25Example41(int n, int m);
  static TestCase getWu25Example42(int n, int m);
  static TestCase getBachirExample1();
  static TestCase getBachirExample2();

  // === INITIAL POINT GENERATORS ===

  /**
   * Get initial feasible starting points for 5x5 problems
   */
  static std::vector<std::vector<double>> getInitial5x5X();
  static std::vector<double> getInitial5x5Y();
  static std::vector<std::vector<double>> getInitial5x5S();

  // === 4x4 PROBLEM DATA ===

  /**
   * Get data for 4x4 problems (used in CQSDO Problem 2)
   */
  static std::vector<std::vector<std::vector<double>>> get4x4AMatrices();
  static std::vector<double> get4x4BVector();
  static std::vector<std::vector<double>> get4x4CMatrix();
  static std::vector<std::vector<double>> getInitial4x4X();
  static std::vector<double> getInitial4x4Y();
  static std::vector<std::vector<double>> getInitial4x4S();

  // === RANDOM LP GENERATORS ===

  /**
   * Generate random linear programming problems for comprehensive testing
   */
  static std::vector<std::vector<std::vector<double>>>
  generateRandomLPAMatrices(int n, int m);
  static std::vector<double> generateRandomLPBVector(int m);
  static std::vector<std::vector<double>> generateRandomLPCMatrix(int n);
  static std::vector<std::vector<double>> generateRandomLPInitialX(int n);
  static std::vector<double> generateRandomLPInitialY(int m);
  static std::vector<std::vector<double>> generateRandomLPInitialS(int n);
};

#endif // TEST_DATA_H
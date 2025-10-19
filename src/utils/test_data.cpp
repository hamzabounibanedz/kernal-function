#include "test_data.h"
#include <vector>
#include <algorithm>

// === PAPER-SPECIFIC TEST CASES ===

TestCase TestDataProvider::getTrigonometricKernelTestCase() {
  TestCase test_case;
  test_case.name = "Trigonometric Kernel Function Test Case";
  test_case.paper_reference = "First Paper - Trigonometric Kernel Function";
  test_case.algorithm_type = "Trigonometric Kernel";
  test_case.problem_size = 5;
  test_case.tolerance = 1e-8;
  test_case.theta = 0.5;
  test_case.tau = 3.0;
  test_case.theta_values = {0.1, 0.3, 0.5, 0.7, 0.9};
  test_case.expected_optimal_value = -1.0957;
  test_case.notes = "Generic Primal-Dual Interior-Point Algorithm for LO with trigonometric kernel ψ(t)";

  // Use the common 5x5 matrices per Derbal/IJNAA
  test_case.A_matrices = getCommon5x5AMatrices();
  test_case.b_vector = getCommon5x5BVector();
  test_case.C_matrix = getCommon5x5CMatrix();
  test_case.initial_X = getInitial5x5X();
  test_case.initial_y = {2.0, 2.0, 2.0};
  test_case.initial_S = getInitial5x5S();

  return test_case;
}

TestCase TestDataProvider::getExponentialParametricTestCase() {
  TestCase test_case;
  test_case.name = "Exponential Parametric Kernel - Example 7.1";
  test_case.paper_reference = "JNFA201814(2).pdf - Example 7.1";
  test_case.algorithm_type = "Exponential Parametric";
  test_case.problem_size = 5;
  test_case.tolerance = 1e-8;
  test_case.theta = 0.4;
  test_case.tau = 1.0;
  test_case.theta_values = {0.05, 0.4, 0.6, 0.95};
  test_case.expected_optimal_value = -1.0957;
  test_case.notes = "Generic Primal-Dual IPM for SDO with exponential parametric kernel";

  // EXACT 5x5 matrices from Example 7.1
  test_case.A_matrices = getCommon5x5AMatrices();
  test_case.b_vector = getCommon5x5BVector();
  test_case.C_matrix = getCommon5x5CMatrix();
  test_case.initial_X = getInitial5x5X();
  test_case.initial_y = {2.0, 2.0, 2.0};
  test_case.initial_S = getInitial5x5S();

  return test_case;
}

TestCase TestDataProvider::getParameterizedLogTestCase() {
  TestCase test_case;
  test_case.name = "Parameterized Logarithmic Kernel - Section 4";
  test_case.paper_reference = "An Efficient Parameterized Logarithmic Kernel Function for Semidefinite Optimization.pdf";
  test_case.algorithm_type = "Parameterized Log";
  test_case.problem_size = 5;
  test_case.tolerance = 1e-8;
  test_case.theta = 0.5;
  test_case.tau = 3.0;
  test_case.theta_values = {0.1, 0.3, 0.5, 0.7, 0.9};
  test_case.expected_optimal_value = -1.0957;
  test_case.notes = "Generic Primal-Dual IPM for SDO based on ψ(t) with parameterized logarithmic kernel";

  // Same 5x5 matrices as other algorithms
  test_case.A_matrices = getCommon5x5AMatrices();
  test_case.b_vector = getCommon5x5BVector();
  test_case.C_matrix = getCommon5x5CMatrix();
  test_case.initial_X = getInitial5x5X();
  test_case.initial_y = getInitial5x5Y();
  test_case.initial_S = getInitial5x5S();

  return test_case;
}

TestCase TestDataProvider::getCQSDOProblem1TestCase() {
  TestCase test_case;
  test_case.name = "CQSDO Problem 1 (n=5, m=3, Q=0)";
  test_case.paper_reference = "34-5-4-8310.pdf - Problem 1";
  test_case.algorithm_type = "CQSDO";
  test_case.problem_size = 5;
  test_case.tolerance = 1e-8;
  test_case.theta = 0.5;
  test_case.tau = 3.0;
  test_case.theta_values = {0.1, 0.3, 0.5, 0.7, 0.9};
  test_case.expected_optimal_value = -1.0957;
  test_case.notes = "Generic Primal-Dual Algorithm for CQSDO - Problem 1 with Q=0";

  // Same Aᵢ,b,C as previous algorithms
  test_case.A_matrices = getCommon5x5AMatrices();
  test_case.b_vector = getCommon5x5BVector();
  test_case.C_matrix = getCommon5x5CMatrix();
  test_case.initial_X = getInitial5x5X();
  test_case.initial_y = getInitial5x5Y();
  test_case.initial_S = getInitial5x5S();

  return test_case;
}

TestCase TestDataProvider::getCQSDOProblem2TestCase() {
  TestCase test_case;
  test_case.name = "CQSDO Problem 2 (n=4, m=3, Q=I)";
  test_case.paper_reference = "34-5-4-8310.pdf - Problem 2";
  test_case.algorithm_type = "CQSDO";
  test_case.problem_size = 4;
  test_case.tolerance = 1e-6;
  test_case.theta = 0.5;
  test_case.tau = 3.0;
  test_case.theta_values = {0.1, 0.3, 0.5, 0.7, 0.9};
  test_case.expected_optimal_value = 0.2101;
  test_case.notes = "Generic Primal-Dual Algorithm for CQSDO - Problem 2 with Q=I (identity matrix)";

  // 4x4 matrices for Problem 2
  test_case.A_matrices = get4x4AMatrices();
  test_case.b_vector = get4x4BVector();
  test_case.C_matrix = get4x4CMatrix();
  test_case.initial_X = getInitial4x4X();
  test_case.initial_y = getInitial4x4Y();
  test_case.initial_S = getInitial4x4S();

  return test_case;
}

// === IJNAA 2023 Section 4 explicit positive matrices case (n=5, m=3) ===
TestCase TestDataProvider::getIJNAA2023Section4TestCase() {
  TestCase t{};
  t.name = "IJNAA 2023 Section 4 (n=5, m=3)";
  t.paper_reference = "IJNAA 2023 Section 4";
  t.algorithm_type = "IJNAA23";
  t.problem_size = 5;
  t.tolerance = 1e-8; t.theta = 0.4; t.tau = 3.0;

  // Matrices as specified (all nonnegative version from item 4 in your list)
  t.A_matrices = {
    {{0,1,0,0,0}, {1,2,0,0,1}, {0,0,0,1,0}, {0,0,1,0,2}, {0,1,0,2,0}},
    {{0,0,2,2,0}, {0,0,2,1,0}, {2,2,1,2,1}, {2,1,2,0,1}, {0,0,1,1,0}},
    {{2,2,1,1,1}, {2,1,2,0,1}, {1,2,0,1,2}, {1,0,1,2,0}, {1,1,2,0,0}}
  };
  t.C_matrix = {
    {3,3,3,1,1}, {3,5,3,1,2}, {3,3,1,1,1}, {1,1,3,1,1}, {2,2,1,1,1}
  };
  t.b_vector = {2,2,2};
  t.initial_X = getInitial5x5X();
  t.initial_y = getInitial5x5Y();
  t.initial_S = getInitial5x5S();
  return t;
}

TestCase TestDataProvider::getBachirTestCase() {
  TestCase test_case;
  test_case.name = "Derbal & Kebbiche (2020) - Bachir Kernel";
  test_case.paper_reference = "Derbal & Kebbiche (2020)";
  test_case.algorithm_type = "Bachir";
  test_case.problem_size = 5;
  test_case.tolerance = 1e-8;
  test_case.theta = 0.5;
  test_case.tau = 3.0;
  test_case.theta_values = {0.1, 0.3, 0.5, 0.7, 0.9};
  test_case.expected_optimal_value = -1.0957;
  test_case.notes = "Generic Primal-Dual Interior-Point Algorithm with Bachir kernel function";

  // Same 5x5 Aᵢ,b,C data
  test_case.A_matrices = getCommon5x5AMatrices();
  test_case.b_vector = getCommon5x5BVector();
  test_case.C_matrix = getCommon5x5CMatrix();
  test_case.initial_X = getInitial5x5X();
  test_case.initial_y = getInitial5x5Y();
  test_case.initial_S = getInitial5x5S();

  return test_case;
}

TestCase TestDataProvider::getBaiKernelFamilyTestCase() {
  TestCase test_case;
  test_case.name = "Bai et al. Kernel Family φₘ(t)";
  test_case.paper_reference = "Sixth Paper - Bai et al. Kernel Family";
  test_case.algorithm_type = "Bai Kernel Family";
  test_case.problem_size = 5;
  test_case.tolerance = 1e-8;
  test_case.theta = 0.5;
  test_case.tau = 3.0;
  test_case.theta_values = {0.1, 0.3, 0.5, 0.7, 0.9};
  test_case.expected_optimal_value = -1.0957;
  test_case.notes = "Generic Primal-Dual Interior-Point Algorithm for SDO with φₘ(t) kernel family";

  // Same 5x5 data for comparison
  test_case.A_matrices = getCommon5x5AMatrices();
  test_case.b_vector = getCommon5x5BVector();
  test_case.C_matrix = getCommon5x5CMatrix();
  test_case.initial_X = getInitial5x5X();
  test_case.initial_y = getInitial5x5Y();
  test_case.initial_S = getInitial5x5S();

  return test_case;
}

// === LEGACY TEST CASES (for backward compatibility) ===

TestCase TestDataProvider::getAlgorithm2TestCase() {
  return getExponentialParametricTestCase();
}

TestCase TestDataProvider::getAlgorithm3TestCase() {
  return getParameterizedLogTestCase();
}

TestCase TestDataProvider::getAlgorithm4Problem1() {
  return getCQSDOProblem1TestCase();
}

TestCase TestDataProvider::getAlgorithm4Problem2() {
  return getCQSDOProblem2TestCase();
}

TestCase TestDataProvider::getAlgorithm5TestCase() {
  return getBachirTestCase();
}

TestCase TestDataProvider::getAlgorithm6TestCase() {
  return getTrigonometricKernelTestCase();
}

// === SPECIALIZED TEST CASES ===

TestCase TestDataProvider::getBachirOptimizedTestCase() {
  TestCase test_case;
  test_case.name = "Bachir Optimized Test Case";
  test_case.paper_reference = "Custom - Bachir Optimized";
  test_case.algorithm_type = "Bachir";
  test_case.problem_size = 5;
  test_case.tolerance = 1e-8;
  test_case.theta = 0.3; // Optimized for Bachir
  test_case.tau = 2.5;
  test_case.theta_values = {0.1, 0.2, 0.3, 0.4, 0.5};
  test_case.expected_optimal_value = -1.0957;
  test_case.notes = "Custom test case optimized for Bachir kernel performance";

  test_case.A_matrices = getCommon5x5AMatrices();
  test_case.b_vector = getCommon5x5BVector();
  test_case.C_matrix = getCommon5x5CMatrix();
  test_case.initial_X = getInitial5x5X();
  test_case.initial_y = getInitial5x5Y();
  test_case.initial_S = getInitial5x5S();

  return test_case;
}

TestCase TestDataProvider::getHighDimensionalBachirTestCase() {
  TestCase test_case;
  test_case.name = "High-Dimensional Bachir Test Case (6x6)";
  test_case.paper_reference = "Custom - High-Dimensional Bachir";
  test_case.algorithm_type = "Bachir";
  test_case.problem_size = 6;
  test_case.tolerance = 1e-8;
  test_case.theta = 0.4;
  test_case.tau = 3.0;
  test_case.theta_values = {0.2, 0.4, 0.6, 0.8};
  test_case.expected_optimal_value = -2.5; // Estimated
  test_case.notes = "6x6 high-dimensional test case designed for Bachir kernel";

  // Generate 6x6 matrices
  test_case.A_matrices = generateRandomLPAMatrices(6, 4);
  test_case.b_vector = generateRandomLPBVector(4);
  test_case.C_matrix = generateRandomLPCMatrix(6);
  test_case.initial_X = generateRandomLPInitialX(6);
  test_case.initial_y = generateRandomLPInitialY(4);
  test_case.initial_S = generateRandomLPInitialS(6);

  return test_case;
}

TestCase TestDataProvider::getRandomLPTestCase() {
  TestCase test_case;
  test_case.name = "Random Linear Programming (50×100)";
  test_case.paper_reference = "Custom - Random LP";
  test_case.algorithm_type = "General";
  test_case.problem_size = 50;
  test_case.tolerance = 1e-6;
  test_case.theta = 0.5;
  test_case.tau = 3.0;
  test_case.theta_values = {0.1, 0.3, 0.5, 0.7, 0.9};
  test_case.expected_optimal_value = 0.0; // Unknown for random
  test_case.notes = "Synthetic 50×100 random linear programming problem";

  test_case.A_matrices = generateRandomLPAMatrices(50, 100);
  test_case.b_vector = generateRandomLPBVector(100);
  test_case.C_matrix = generateRandomLPCMatrix(50);
  test_case.initial_X = generateRandomLPInitialX(50);
  test_case.initial_y = generateRandomLPInitialY(100);
  test_case.initial_S = generateRandomLPInitialS(50);

  return test_case;
}

std::vector<TestCase> TestDataProvider::getAllTestCases() {
  // Restrict to paper-accurate tests only (no synthetic/random)
  std::vector<TestCase> testCases = {
      getExponentialParametricTestCase(),            // Derbal20 5x5 (A,b,C exact)
      getIJNAA2023Section4TestCase(),                // IJNAA 2023 Section 4 positive matrices
      getCQSDOProblem1TestCase(),                    // CQSDO Problem 1 (Q=0)
      getCQSDOProblem2TestCase(),                    // CQSDO Problem 2 (Q=I)
      getBachirTestCase()                            // Bachir test case
  };

  return testCases;
}

std::vector<TestCase> TestDataProvider::getTestCasesByAlgorithm(const std::string& algorithm_type) {
  auto allCases = getAllTestCases();
  std::vector<TestCase> filteredCases;
  
  for (const auto& testCase : allCases) {
    if (testCase.algorithm_type == algorithm_type) {
      filteredCases.push_back(testCase);
    }
  }
  
  return filteredCases;
}

std::vector<TestCase> TestDataProvider::getTestCasesByPaper(const std::string& paper_reference) {
  auto allCases = getAllTestCases();
  std::vector<TestCase> filteredCases;
  
  for (const auto& testCase : allCases) {
    if (testCase.paper_reference.find(paper_reference) != std::string::npos) {
      filteredCases.push_back(testCase);
    }
  }
  
  return filteredCases;
}

// === COMMON DATA GENERATORS ===

// Exact 5x5 matrices from the literature example (Derbal & Kebbiche / IJNAA23)
std::vector<std::vector<std::vector<double>>>
TestDataProvider::getCommon5x5AMatrices() {
  // A₁
  std::vector<std::vector<double>> A1 = {{0, 1, 0, 0, 0},
                                         {1, 2, 0, 0, -1},
                                         {0, 0, 0, 0, 1},
                                         {0, 0, 0, -2, -1},
                                         {0, -1, 1, -1, -2}};

  // A₂
  std::vector<std::vector<double>> A2 = {{0, 0, -2, 2, 0},
                                         {0, 2, 1, 0, 2},
                                         {-2, 1, -2, 0, 1},
                                         {2, 0, 0, 0, 0},
                                         {0, 2, 1, 0, 2}};

  // A₃
  std::vector<std::vector<double>> A3 = {{2, 2, -1, -1, 1},
                                         {2, 0, 2, 1, 1},
                                         {-1, 2, 0, 1, 0},
                                         {-1, 1, 1, -2, 0},
                                         {1, 1, 0, 0, -2}};

  return {A1, A2, A3};
}

std::vector<double> TestDataProvider::getCommon5x5BVector() {
  // b from the literature example
  return {-2.0, 2.0, -2.0};
}

std::vector<std::vector<double>> TestDataProvider::getCommon5x5CMatrix() {
  // C from the literature example
  return {{3, 3, -3, 1, 1},
          {3, 5,  3, 1, 2},
          {-3,3, -1, 1, 2},
          {1, 1,  1,-3,-1},
          {1, 2,  2,-1,-1}};
}

std::vector<std::vector<double>> TestDataProvider::getInitial5x5X() {
  // Start with identity matrix (typical initialization)
  return {{1, 0, 0, 0, 0},
          {0, 1, 0, 0, 0},
          {0, 0, 1, 0, 0},
          {0, 0, 0, 1, 0},
          {0, 0, 0, 0, 1}};
}

std::vector<double> TestDataProvider::getInitial5x5Y() {
  // Start with (1,1,1) as mentioned in requirements
  return {1.0, 1.0, 1.0};
}

std::vector<std::vector<double>> TestDataProvider::getInitial5x5S() {
  // Start with identity matrix
  return {{1, 0, 0, 0, 0},
          {0, 1, 0, 0, 0},
          {0, 0, 1, 0, 0},
          {0, 0, 0, 1, 0},
          {0, 0, 0, 0, 1}};
}

// === 4x4 PROBLEM DATA ===

std::vector<std::vector<std::vector<double>>> TestDataProvider::get4x4AMatrices() {
  // 4x4 constraint matrices for CQSDO Problem 2
  std::vector<std::vector<double>> A1 = {{1, 0, 0, 0},
                                         {0, 1, 0, 0},
                                         {0, 0, 1, 0},
                                         {0, 0, 0, 1}};

  std::vector<std::vector<double>> A2 = {{0, 1, 1, 0},
                                         {1, 0, 0, 1},
                                         {1, 0, 0, 1},
                                         {0, 1, 1, 0}};

  std::vector<std::vector<double>> A3 = {{1, 1, 0, 0},
                                         {1, 1, 0, 0},
                                         {0, 0, 1, 1},
                                         {0, 0, 1, 1}};

  return {A1, A2, A3};
}

// === PAPER PATTERN TESTS ===
TestCase TestDataProvider::getTouil22Problem1(int n, int m) {
  TestCase t{};
  t.name = "Touil22 Problem 1 (C=I)";
  t.paper_reference = "Touil & Chikouche (2022)";
  t.algorithm_type = "Touil22";
  t.problem_size = n;
  t.tolerance = 1e-8;
  t.theta = 0.5;
  t.tau = 3.0;
  // C = I
  t.C_matrix.assign(n, std::vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) t.C_matrix[i][i] = 1.0;
  // A pattern and b
  t.A_matrices.assign(m, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i == j && i == k) t.A_matrices[k][i][j] = 1.0;
        else if (i != j && i == m - k - 1) t.A_matrices[k][i][j] = 1.0;
      }
    }
  }
  t.b_vector.assign(m, 2.0);
  // X0: 1.5 diagonal, 0.5 off-diagonal
  t.initial_X.assign(n, std::vector<double>(n, 0.5));
  for (int i = 0; i < n; ++i) t.initial_X[i][i] = 1.5;
  // y0 and S0
  t.initial_y.assign(m, 2.0);
  t.initial_S.assign(n, std::vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) t.initial_S[i][i] = 1.0;
  return t;
}

TestCase TestDataProvider::getTouil22Problem2(int n, int m) {
  TestCase t = getTouil22Problem1(n, m);
  t.name = "Touil22 Problem 2 (C has (m,m)=1)";
  // C(i,j)=1 if i=j=m
  t.C_matrix.assign(n, std::vector<double>(n, 0.0));
  if (m - 1 < n) t.C_matrix[m - 1][m - 1] = 1.0;
  // S0 special: S0(i,j) = 1 if i=j=m; 2 if i≠j=m; 0 otherwise
  t.initial_S.assign(n, std::vector<double>(n, 0.0));
  if (m - 1 < n) {
    t.initial_S[m - 1][m - 1] = 1.0;
    for (int i = 0; i < n; ++i) {
      if (i != m - 1) {
        t.initial_S[i][m - 1] = 2.0;
        t.initial_S[m - 1][i] = 2.0;
      }
    }
  }
  // X0 = I
  t.initial_X.assign(n, std::vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) t.initial_X[i][i] = 1.0;
  return t;
}

TestCase TestDataProvider::getWu25Example41(int n, int m) {
  // Same pattern as Touil22 P1
  TestCase t = getTouil22Problem1(n, m);
  t.name = "Wu25 Example 4.1";
  t.paper_reference = "Wu & Zhang (2025)";
  t.algorithm_type = "Wu25";
  return t;
}

TestCase TestDataProvider::getWu25Example42(int n, int m) {
  // Same pattern as Touil22 P2
  TestCase t = getTouil22Problem2(n, m);
  t.name = "Wu25 Example 4.2";
  t.paper_reference = "Wu & Zhang (2025)";
  t.algorithm_type = "Wu25";
  return t;
}

TestCase TestDataProvider::getBachirExample1() {
  TestCase t{};
  t.name = "Bachir Example 1 (3x3)";
  t.paper_reference = "Bachir (2025)";
  t.algorithm_type = "Bachir";
  t.problem_size = 3;
  t.tolerance = 1e-8;
  t.theta = 0.5;
  t.tau = 3.0;
  t.C_matrix = {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}};
  t.A_matrices = {
      {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}},
      {{0, 1, 0}, {1, 0, 0}, {0, 0, 0}},
      {{0, 0, 1}, {0, 0, 0}, {1, 0, 0}}};
  t.b_vector = {2, 2, 2};
  t.initial_X = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  t.initial_y = {2, 2, 2};
  t.initial_S = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  return t;
}

TestCase TestDataProvider::getBachirExample2() {
  TestCase t{};
  t.name = "Bachir Example 2 (4x4)";
  t.paper_reference = "Bachir (2025)";
  t.algorithm_type = "Bachir";
  t.problem_size = 4;
  t.tolerance = 1e-8;
  t.theta = 0.5;
  t.tau = 3.0;
  t.C_matrix = {{4, 0, 0, 0}, {0, 4, 0, 0}, {0, 0, 4, 0}, {0, 0, 0, 4}};
  // Simple diagonal/block patterns for A_i; adjust per paper if needed
  t.A_matrices = get4x4AMatrices();
  t.b_vector = {2, 2, 2};
  t.initial_X = getInitial4x4X();
  t.initial_y = getInitial4x4Y();
  t.initial_S = getInitial4x4S();
  return t;
}

std::vector<double> TestDataProvider::get4x4BVector() {
  return {1.0, 2.0, 1.0};
}

std::vector<std::vector<double>> TestDataProvider::get4x4CMatrix() {
  return {{2, 1, 0, 0},
          {1, 2, 1, 0},
          {0, 1, 2, 1},
          {0, 0, 1, 2}};
}

std::vector<std::vector<double>> TestDataProvider::getInitial4x4X() {
  return {{1, 0, 0, 0},
          {0, 1, 0, 0},
          {0, 0, 1, 0},
          {0, 0, 0, 1}};
}

std::vector<double> TestDataProvider::getInitial4x4Y() {
  return {1.0, 1.0, 1.0};
}

std::vector<std::vector<double>> TestDataProvider::getInitial4x4S() {
  return {{1, 0, 0, 0},
          {0, 1, 0, 0},
          {0, 0, 1, 0},
          {0, 0, 0, 1}};
}

// === RANDOM LP GENERATORS ===

std::vector<std::vector<std::vector<double>>>
TestDataProvider::generateRandomLPAMatrices(int n, int m) {
  std::vector<std::vector<std::vector<double>>> A_matrices(m);
  
  for (int i = 0; i < m; ++i) {
    A_matrices[i] = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
    
    // Generate symmetric matrix
    for (int row = 0; row < n; ++row) {
      for (int col = row; col < n; ++col) {
        double value = (static_cast<double>(rand()) / RAND_MAX) * 2.0 - 1.0;
        A_matrices[i][row][col] = value;
        A_matrices[i][col][row] = value; // Symmetric
      }
    }
  }
  
  return A_matrices;
}

std::vector<double> TestDataProvider::generateRandomLPBVector(int m) {
  std::vector<double> b_vector(m);
  
  for (int i = 0; i < m; ++i) {
    b_vector[i] = (static_cast<double>(rand()) / RAND_MAX) * 4.0 - 2.0;
  }
  
  return b_vector;
}

std::vector<std::vector<double>> TestDataProvider::generateRandomLPCMatrix(int n) {
  std::vector<std::vector<double>> C_matrix(n, std::vector<double>(n, 0.0));
  
  // Generate symmetric matrix
  for (int row = 0; row < n; ++row) {
    for (int col = row; col < n; ++col) {
      double value = (static_cast<double>(rand()) / RAND_MAX) * 4.0 - 2.0;
      C_matrix[row][col] = value;
      C_matrix[col][row] = value; // Symmetric
    }
  }
  
  return C_matrix;
}

std::vector<std::vector<double>> TestDataProvider::generateRandomLPInitialX(int n) {
  std::vector<std::vector<double>> X(n, std::vector<double>(n, 0.0));
  
  // Generate positive definite matrix
  for (int row = 0; row < n; ++row) {
    for (int col = row; col < n; ++col) {
      double value = (static_cast<double>(rand()) / RAND_MAX) * 0.5 + 0.5;
      X[row][col] = value;
      X[col][row] = value; // Symmetric
    }
  }
  
  // Ensure positive definiteness by adding diagonal dominance
  for (int i = 0; i < n; ++i) {
    X[i][i] += n;
  }
  
  return X;
}

std::vector<double> TestDataProvider::generateRandomLPInitialY(int m) {
  std::vector<double> y(m);
  
  for (int i = 0; i < m; ++i) {
    y[i] = (static_cast<double>(rand()) / RAND_MAX) * 2.0 - 1.0;
  }
  
  return y;
}

std::vector<std::vector<double>> TestDataProvider::generateRandomLPInitialS(int n) {
  std::vector<std::vector<double>> S(n, std::vector<double>(n, 0.0));
  
  // Generate positive definite matrix
  for (int row = 0; row < n; ++row) {
    for (int col = row; col < n; ++col) {
      double value = (static_cast<double>(rand()) / RAND_MAX) * 0.5 + 0.5;
      S[row][col] = value;
      S[col][row] = value; // Symmetric
    }
  }
  
  // Ensure positive definiteness by adding diagonal dominance
  for (int i = 0; i < n; ++i) {
    S[i][i] += n;
  }
  
  return S;
}
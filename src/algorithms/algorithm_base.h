#ifndef ALGORITHM_BASE_H
#define ALGORITHM_BASE_H

#include "../kernels/kernel_base.h"
#include "../utils/test_data.h"
#include <memory>
#include <string>
#include <vector>

/**
 * AlgorithmResult - Stores the results of running an optimization algorithm
 *
 * When we run an algorithm to solve an optimization problem, we get back
 * various pieces of information about how it performed. This structure
 * organizes all that information in one place.
 */
struct AlgorithmResult {
  bool converged;           // Did the algorithm find a solution?
  int iterations;           // How many steps did it take?
  double final_mu;          // Final barrier parameter value
  double primal_objective;  // Value of the primal objective function
  double dual_objective;    // Value of the dual objective function
  double duality_gap;       // Gap between primal and dual solutions
  double execution_time_ms; // How long did it take (in milliseconds)?
  std::vector<std::vector<double>>
      primal_solution;                 // The solution X (or x) that was found
  std::vector<double> dual_solution_y; // The dual solution y
  std::vector<std::vector<double>>
      dual_solution_s;          // The dual solution S (or s)
  std::string convergence_info; // Any additional information about convergence
};

// TestCase struct is now defined in ../utils/test_data.h

/**
 * AlgorithmBase - The foundation for all optimization algorithms
 *
 * This is the base class that all our optimization algorithms inherit from.
 * It defines what every algorithm must be able to do:
 * - Solve optimization problems
 * - Provide information about itself
 * - Handle parameters and test cases
 *
 * Each specific algorithm (like Algorithm 1, Algorithm 2, etc.) inherits
 * from this class and implements the actual optimization logic in its own way.
 * This allows us to easily compare different algorithms using the same
 * interface.
 */
class AlgorithmBase {
public:
  // Constructor - creates an algorithm that uses a specific kernel function
  // The kernel function determines how the algorithm approaches the solution
  explicit AlgorithmBase(std::shared_ptr<KernelBase> kernel);

  // Destructor - cleans up when the algorithm is no longer needed
  virtual ~AlgorithmBase() = default;

  /**
   * Solve an optimization problem
   *
   * This is the main function - it takes a test problem and tries to find
   * the best solution using the algorithm's approach. The result tells us
   * whether it succeeded, how long it took, and what solution it found.
   *
   * @param test_case The optimization problem to solve
   * @return Detailed results of the algorithm's execution
   */
  virtual AlgorithmResult solve(const TestCase &test_case) = 0;

  /**
   * Get a human-readable name for this algorithm
   *
   * Used in the user interface to show which algorithm is being used.
   * For example: "Algorithm 1: Exponential-Parametric" or "Algorithm 2: Log
   * Kernel"
   *
   * @return A string describing the algorithm
   */
  virtual std::string getName() const = 0;

  /**
   * Get a detailed description of this algorithm
   *
   * This provides more information about how the algorithm works,
   * what makes it unique, and when it might be useful.
   *
   * @return A detailed description of the algorithm
   */
  virtual std::string getDescription() const = 0;

  /**
   * Get the kernel function used by this algorithm
   *
   * Each algorithm uses a specific kernel function that determines
   * how it approaches the optimization problem.
   *
   * @return A pointer to the kernel function
   */
  std::shared_ptr<KernelBase> getKernel() const { return m_kernel; }

  /**
   * Set the algorithm's parameters
   *
   * These parameters control how the algorithm behaves:
   * - tolerance: How close to the solution is "good enough"
   * - theta: Controls step size and convergence behavior
   * - tau: Influences barrier parameter updates (default is 1.0)
   *
   * @param tolerance Convergence tolerance
   * @param theta Step size parameter
   * @param tau Barrier parameter update factor
   */
  virtual void setParameters(double tolerance, double theta, double tau = 1.0);

  /**
   * Get test cases that work well with this algorithm
   *
   * Different algorithms might work better with different types of problems.
   * This function returns a set of test problems that are particularly
   * suitable for testing this algorithm's performance.
   *
   * @return A vector of test cases
   */
  virtual std::vector<TestCase> getTestCases() const = 0;

protected:
  // === ALGORITHM COMPONENTS ===

  // The kernel function that this algorithm uses
  std::shared_ptr<KernelBase> m_kernel;

  // === ALGORITHM PARAMETERS ===

  // How close to the solution is "good enough" (convergence criterion)
  double m_tolerance;

  // Controls step size and convergence behavior (0 < θ < 1)
  double m_theta;

  // Influences barrier parameter updates (typically τ = 1.0)
  double m_tau;

  // === HELPER FUNCTIONS ===

  /**
   * Check if the input problem is valid for this algorithm
   *
   * Before trying to solve a problem, we need to make sure it's well-formed
   * and suitable for this algorithm. This function performs those checks.
   *
   * @param test_case The problem to validate
   * @return true if the problem is valid, false otherwise
   */
  virtual bool validateInput(const TestCase &test_case) const;

  /**
   * Calculate the initial barrier parameter value
   *
   * The barrier parameter μ controls how close we are to the boundary
   * of the feasible region. We need to start with a reasonable value.
   *
   * @param test_case The problem we're solving
   * @return The initial value of μ
   */
  virtual double computeInitialMu(const TestCase &test_case) const;

  /**
   * Check if the algorithm has converged to a solution
   *
   * We need to know when to stop the algorithm. This function checks
   * if we're close enough to the solution to consider the problem solved.
   *
   * @param mu Current barrier parameter value
   * @param duality_gap Current gap between primal and dual solutions
   * @return true if converged, false if we should continue
   */
  virtual bool checkConvergence(double mu, double duality_gap) const;

  // === MATRIX UTILITY FUNCTIONS ===

  // These functions help with matrix operations that are common in optimization

  /**
   * Calculate the trace of a matrix (sum of diagonal elements)
   *
   * @param matrix The matrix to compute the trace of
   * @return The trace value
   */
  double matrixTrace(const std::vector<std::vector<double>> &matrix) const;

  /**
   * Multiply two matrices together
   *
   * @param A First matrix
   * @param B Second matrix
   * @return The product A * B
   */
  std::vector<std::vector<double>>
  matrixMultiply(const std::vector<std::vector<double>> &A,
                 const std::vector<std::vector<double>> &B) const;

  /**
   * Add two matrices together
   *
   * @param A First matrix
   * @param B Second matrix
   * @return The sum A + B
   */
  std::vector<std::vector<double>>
  matrixAdd(const std::vector<std::vector<double>> &A,
            const std::vector<std::vector<double>> &B) const;

  /**
   * Subtract one matrix from another
   *
   * @param A First matrix
   * @param B Second matrix
   * @return The difference A - B
   */
  std::vector<std::vector<double>>
  matrixSubtract(const std::vector<std::vector<double>> &A,
                 const std::vector<std::vector<double>> &B) const;

  /**
   * Check if a matrix is positive definite
   *
   * Positive definiteness is an important property in optimization.
   * This function checks if a matrix has this property.
   *
   * @param matrix The matrix to check
   * @return true if positive definite, false otherwise
   */
  bool isPositiveDefinite(const std::vector<std::vector<double>> &matrix) const;
};

#endif // ALGORITHM_BASE_H
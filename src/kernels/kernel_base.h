#ifndef KERNEL_BASE_H
#define KERNEL_BASE_H

#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

/**
 * KernelBase - The foundation for all kernel functions
 *
 * In optimization theory, a kernel function is a mathematical tool that helps
 * algorithms find solutions to complex problems. Think of it as a "guide" that
 * tells the algorithm how to move toward the best solution.
 *
 * This base class defines what every kernel function must be able to do:
 * - Calculate its value at any point
 * - Calculate its derivatives (how it changes)
 * - Provide information about itself
 *
 * Each specific kernel function (like exponential, logarithmic, etc.) inherits
 * from this class and implements these calculations in its own way.
 */
class KernelBase {
public:
  // Virtual destructor - needed because this is a base class
  // When we delete a kernel object, we need to make sure the right destructor
  // is called
  virtual ~KernelBase() = default;

  /**
   * Calculate the kernel function value ψ(t)
   *
   * This is the main function - it takes a number t and returns the value
   * of the kernel function at that point. The value t must be positive
   * because kernel functions are typically defined only for positive numbers.
   *
   * @param t The input value (must be positive)
   * @return The value of the kernel function ψ(t)
   */
  virtual double psi(double t) const = 0;

  /**
   * Calculate the first derivative ψ'(t)
   *
   * The first derivative tells us how fast the function is changing at point t.
   * This is crucial for optimization algorithms because they need to know
   * which direction to move to improve the solution.
   *
   * @param t The input value (must be positive)
   * @return The first derivative ψ'(t)
   */
  virtual double psi_prime(double t) const = 0;

  /**
   * Calculate the second derivative ψ''(t)
   *
   * The second derivative tells us about the curvature of the function.
   * It helps algorithms understand how the rate of change is itself changing,
   * which is important for determining step sizes and convergence behavior.
   *
   * @param t The input value (must be positive)
   * @return The second derivative ψ''(t)
   */
  virtual double psi_double_prime(double t) const = 0;

  /**
   * Calculate the third derivative ψ'''(t)
   *
   * The third derivative provides even more detailed information about
   * how the function behaves. While not always needed, it can help with
   * advanced optimization techniques and error analysis.
   *
   * @param t The input value (must be positive)
   * @return The third derivative ψ'''(t)
   */
  virtual double psi_triple_prime(double t) const = 0;

  /**
   * Get a human-readable name for this kernel function
   *
   * This is used in the user interface to show which kernel is being used.
   * For example: "Exponential-Parametric Kernel" or "Logarithmic Kernel"
   *
   * @return A string describing the kernel function
   */
  virtual std::string getName() const = 0;

  /**
   * Get the mathematical formula for this kernel function
   *
   * This returns a string representation of the mathematical formula,
   * useful for displaying in the user interface or documentation.
   * The format is similar to LaTeX notation.
   *
   * @return A string with the mathematical formula
   */
  virtual std::string getFormula() const = 0;

  /**
   * Get the current parameter values for this kernel
   *
   * Many kernel functions have parameters that can be adjusted
   * (like θ, τ, p, etc.). This function returns the current values
   * of those parameters.
   *
   * @return A vector containing the parameter values
   */
  virtual std::vector<double> getParameters() const = 0;

  /**
   * Set new parameter values for this kernel
   *
   * This allows users to adjust the behavior of the kernel function
   * by changing its parameters. For example, changing θ from 0.5 to 0.8
   * might make the algorithm converge faster or slower.
   *
   * @param params A vector containing the new parameter values
   */
  virtual void setParameters(const std::vector<double> &params) = 0;

  /**
   * Check if an input value is valid for this kernel function
   *
   * Most kernel functions only work with positive numbers, so this
   * function checks if the input t is acceptable. If not, the kernel
   * might return nonsense results or crash.
   *
   * @param t The input value to check
   * @return true if t is valid, false otherwise
   */
  virtual bool isValidInput(double t) const { return t > 0.0; }

  /**
   * Evaluate the kernel function over a range of values
   *
   * This is useful for plotting the kernel function or understanding
   * its behavior over a range of inputs. It returns pairs of (t, ψ(t))
   * values that can be used for visualization.
   *
   * @param t_min The minimum t value to evaluate
   * @param t_max The maximum t value to evaluate
   * @param num_points How many points to evaluate between min and max
   * @return A vector of (t, ψ(t)) pairs
   */
  std::vector<std::pair<double, double>>
  evaluateRange(double t_min, double t_max, int num_points) const;

  /**
   * Calculate all derivatives at once for efficiency
   *
   * Sometimes we need all the derivatives at the same point. Instead
   * of calling each derivative function separately, this function
   * calculates them all at once, which can be more efficient.
   *
   * @param t The input value
   * @return A tuple containing (ψ(t), ψ'(t), ψ''(t), ψ'''(t))
   */
  std::tuple<double, double, double, double> evaluateAll(double t) const;

protected:
  /**
   * Helper function to validate input and throw an error if invalid
   *
   * This is used internally by kernel functions to check their inputs
   * and provide helpful error messages if something goes wrong.
   *
   * @param t The input value to validate
   * @param function_name The name of the calling function (for error messages)
   */
  void validateInput(double t, const std::string &function_name) const;
};

#endif // KERNEL_BASE_H
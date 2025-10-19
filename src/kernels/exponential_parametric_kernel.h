#ifndef EXPONENTIAL_PARAMETRIC_KERNEL_H
#define EXPONENTIAL_PARAMETRIC_KERNEL_H

#include "kernel_base.h"
#include <cmath>

/**
 * Exponential-Parametric Kernel Function - Algorithm 1
 *
 * This kernel function uses an exponential decay approach with a parameter p
 * that you can adjust to control how the function behaves. The mathematical
 * formula is:
 *
 * ψ(t) = (t²-1)/2 - ∫₁ᵗ ((e-1)/(e^x-1))^p dx
 *
 * The parameter p (≥ 1) controls how fast the exponential decay happens.
 * When p = 1, you get a standard exponential decay. When p > 1, the decay
 * becomes more aggressive.
 *
 * This kernel is particularly useful for problems where you want smooth,
 * controlled convergence behavior. The exponential term helps the algorithm
 * approach the solution gradually without overshooting.
 */
class ExponentialParametricKernel : public KernelBase {
public:
  // Constructor - you can set the parameter p when creating the kernel
  // Default value is 1.0, which gives standard exponential behavior
  explicit ExponentialParametricKernel(double p = 1.0);

  // Destructor - cleans up when the kernel is no longer needed
  virtual ~ExponentialParametricKernel() = default;

  // === CORE KERNEL FUNCTION CALCULATIONS ===

  // Calculate the main kernel function value ψ(t)
  double psi(double t) const override;

  // Calculate the first derivative ψ'(t) = t - ((e-1)/(e^t-1))^p
  double psi_prime(double t) const override;

  // Calculate the second derivative ψ''(t)
  double psi_double_prime(double t) const override;

  // Calculate the third derivative ψ'''(t)
  double psi_triple_prime(double t) const override;

  // === METADATA AND INFORMATION ===

  // Get a human-readable name for this kernel
  std::string getName() const override;

  // Get the mathematical formula as a string
  std::string getFormula() const override;

  // Get the current parameter values (just p in this case)
  std::vector<double> getParameters() const override;

  // Set new parameter values (allows changing p during runtime)
  void setParameters(const std::vector<double> &params) override;

private:
  // The parameter p that controls the exponential decay behavior
  // p ≥ 1, where larger values make the decay more aggressive
  double m_p;

  /**
   * Calculate the integral term numerically
   *
   * The kernel function has an integral that can't be solved analytically
   * for most values of p, so we compute it numerically using numerical
   * integration techniques.
   *
   * @param t The upper limit of integration
   * @return The value of ∫₁ᵗ ((e-1)/(e^x-1))^p dx
   */
  double computeIntegral(double t) const;

  /**
   * Helper function g(t) = ((e-1)/(e^t-1))^p
   *
   * This is the function inside the integral. It appears in the derivatives
   * as well, so we compute it separately to avoid code duplication.
   *
   * @param t The input value
   * @return The value of g(t)
   */
  double g(double t) const;

  /**
   * First derivative of the helper function g'(t)
   *
   * This is needed for computing the second derivative of the kernel.
   *
   * @param t The input value
   * @return The first derivative g'(t)
   */
  double g_prime(double t) const;

  /**
   * Second derivative of the helper function g''(t)
   *
   * This is needed for computing the third derivative of the kernel.
   *
   * @param t The input value
   * @return The second derivative g''(t)
   */
  double g_double_prime(double t) const;

  // Mathematical constant e (Euler's number) for exponential calculations
  // Using a high-precision value to ensure accurate computations
  static constexpr double E = 2.71828182845904523536;
};

#endif // EXPONENTIAL_PARAMETRIC_KERNEL_H
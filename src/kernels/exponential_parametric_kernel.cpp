#include "exponential_parametric_kernel.h"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// Constructor - creates a new exponential-parametric kernel
// The parameter p controls how aggressive the exponential decay is
// We ensure p is at least 1.0 to maintain mathematical properties
ExponentialParametricKernel::ExponentialParametricKernel(double p)
    : m_p(std::max(1.0, p)) {}

// Calculate the main kernel function value ψ(t)
// This is the heart of the kernel - it combines a quadratic term with an
// integral
double ExponentialParametricKernel::psi(double t) const {
  // First, make sure the input is valid (t must be positive)
  validateInput(t, "psi");

  // The formula is: ψ(t) = (t²-1)/2 - ∫₁ᵗ ((e-1)/(e^x-1))^p dx
  // We break it into two parts for clarity:

  // Part 1: The quadratic term (t²-1)/2
  // This gives the basic shape of the kernel function
  double term1 = (t * t - 1.0) / 2.0;

  // Part 2: The integral term (computed numerically)
  // This adds the exponential behavior controlled by parameter p
  double term2 = computeIntegral(t);

  // Return the difference: quadratic term minus integral term
  return term1 - term2;
}

// Calculate the first derivative ψ'(t)
// This tells us how fast the function is changing at point t
double ExponentialParametricKernel::psi_prime(double t) const {
  validateInput(t, "psi_prime");

  // The derivative formula is: ψ'(t) = t - ((e-1)/(e^t-1))^p
  // This is much simpler than the main function because the derivative
  // of the integral term cancels out nicely
  return t - g(t);
}

// Calculate the second derivative ψ''(t)
// This tells us about the curvature of the function
double ExponentialParametricKernel::psi_double_prime(double t) const {
  validateInput(t, "psi_double_prime");

  // The second derivative is: ψ''(t) = 1 - g'(t)
  // The derivative of t is 1, and we subtract the derivative of g(t)
  return 1.0 - g_prime(t);
}

// Calculate the third derivative ψ'''(t)
// This provides even more detailed information about the function's behavior
double ExponentialParametricKernel::psi_triple_prime(double t) const {
  validateInput(t, "psi_triple_prime");

  // The third derivative is: ψ'''(t) = -g''(t)
  // The derivative of 1 is 0, so we just get the negative of g's second
  // derivative
  return -g_double_prime(t);
}

// Compute the integral term numerically using Simpson's rule
// This integral can't be solved analytically for most values of p,
// so we use numerical integration to get an approximation
double ExponentialParametricKernel::computeIntegral(double t) const {
  // If t is very close to 1, the integral is 0 (since we're integrating from 1
  // to 1)
  if (std::abs(t - 1.0) < 1e-10) {
    return 0.0;
  }

  // Use Simpson's rule for numerical integration
  // This method gives good accuracy with reasonable computational cost
  int n = 1000; // Number of intervals - more intervals = better accuracy
  double h = (t - 1.0) / n; // Step size
  double sum = 0.0;

  // Simpson's rule: ∫[a,b] f(x)dx ≈ (h/3) * [f(a) + 4f(x₁) + 2f(x₂) + 4f(x₃) +
  // ... + f(b)]
  for (int i = 0; i <= n; ++i) {
    double x = 1.0 + i * h; // Current x value

    // Weight depends on position: endpoints get weight 1, odd indices get
    // weight 4, even indices get weight 2
    double weight = (i == 0 || i == n) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);

    sum += weight * g(x); // Add weighted function value
  }

  // Apply the Simpson's rule formula
  return sum * h / 3.0;
}

// Helper function g(t) = ((e-1)/(e^t-1))^p
// This is the function inside the integral, and it appears in the derivatives
// too
double ExponentialParametricKernel::g(double t) const {
  // Calculate the components step by step for clarity
  double e_minus_1 = E - 1.0;         // (e-1) term
  double exp_t = std::exp(t);         // e^t
  double exp_t_minus_1 = exp_t - 1.0; // (e^t-1) term

  // Handle the case when t is very close to 0 (to avoid division by zero)
  if (exp_t_minus_1 < 1e-10) {
    // When t ≈ 0, e^t ≈ 1, so (e^t-1) ≈ 0
    // In this limit, g(t) approaches (e-1)^p
    return std::pow(e_minus_1, m_p);
  }

  // Calculate the ratio (e-1)/(e^t-1) and raise it to the power p
  double ratio = e_minus_1 / exp_t_minus_1;
  return std::pow(ratio, m_p);
}

// First derivative of g(t)
// This is needed for computing the second derivative of the kernel
double ExponentialParametricKernel::g_prime(double t) const {
  // The derivative formula is: g'(t) = -p * ((e-1)/(e^t-1))^p * e^t / (e^t-1)
  double exp_t = std::exp(t);
  double exp_t_minus_1 = exp_t - 1.0;

  // Handle the case when t is very close to 0
  if (exp_t_minus_1 < 1e-10) {
    return 0.0;
  }

  // Use the fact that g(t) = ((e-1)/(e^t-1))^p
  double g_val = g(t);

  // Apply the derivative formula
  return -m_p * g_val * exp_t / exp_t_minus_1;
}

// Second derivative of g(t)
// This is needed for computing the third derivative of the kernel
double ExponentialParametricKernel::g_double_prime(double t) const {
  // The second derivative has a more complex formula
  double exp_t = std::exp(t);
  double exp_t_minus_1 = exp_t - 1.0;

  // Handle the case when t is very close to 0
  if (exp_t_minus_1 < 1e-10) {
    return 0.0;
  }

  // Get the values we need
  double g_val = g(t);
  double g_prime_val = g_prime(t);

  // The second derivative has two terms
  double term1 = -m_p * g_prime_val * exp_t / exp_t_minus_1;
  double term2 =
      -m_p * g_val * exp_t / exp_t_minus_1 * (1.0 + exp_t / exp_t_minus_1);

  return term1 + term2;
}

// Return a human-readable name for this kernel
std::string ExponentialParametricKernel::getName() const {
  return "Exponential-Parametric Kernel";
}

// Return the mathematical formula as a string
std::string ExponentialParametricKernel::getFormula() const {
  return "ψ(t) = (t²-1)/2 - ∫₁ᵗ ((e-1)/(e^x-1))^p dx, p≥1";
}

// Get the current parameter values
// In this case, we only have one parameter: p
std::vector<double> ExponentialParametricKernel::getParameters() const {
  return {m_p};
}

// Set new parameter values
// This allows users to change the parameter p during runtime
void ExponentialParametricKernel::setParameters(
    const std::vector<double> &params) {
  if (!params.empty()) {
    // Ensure p is at least 1.0 to maintain mathematical properties
    m_p = std::max(1.0, params[0]);
  }
}
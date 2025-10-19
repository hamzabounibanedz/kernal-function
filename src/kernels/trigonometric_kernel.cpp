#include "trigonometric_kernel.h"
#include <cmath>
#include <string>
#include <vector>

double TrigonometricKernel::psi(double t) const {
  // EXACT formula: ψ(t) = (t-1)²/(2t) + (t-1)²/2 + tan²(h(t))/8
  validateInput(t, "psi");

  double t_minus_1 = t - 1.0;
  double term1 = (t_minus_1 * t_minus_1) / (2.0 * t);
  double term2 = (t_minus_1 * t_minus_1) / 2.0;

  double h_val = h(t);
  double tan_h = std::tan(h_val);
  double term3 = (tan_h * tan_h) / 8.0;

  return term1 + term2 + term3;
}

double TrigonometricKernel::psi_prime(double t) const {
  // EXACT formula: ψ'(t) = (2t³-t²-1)/(2t²) + [h'(t)tan(h(t))(1+tan²(h(t)))]/4
  validateInput(t, "psi_prime");

  double t2 = t * t;
  double t3 = t2 * t;

  double term1 = (2.0 * t3 - t2 - 1.0) / (2.0 * t2);

  double h_val = h(t);
  double h_prime_val = h_prime(t);
  double tan_h = std::tan(h_val);
  double tan_h_squared = tan_h * tan_h;
  double term2 = (h_prime_val * tan_h * (1.0 + tan_h_squared)) / 4.0;

  return term1 + term2;
}

double TrigonometricKernel::psi_double_prime(double t) const {
  // EXACT formula: ψ''(t) = (1+t³)/t³ + [(1+tan²(h(t)))/4][h''(t)tan(h(t)) +
  // h'(t)²(1+3tan²(h(t)))]
  validateInput(t, "psi_double_prime");

  double t3 = t * t * t;
  double term1 = (1.0 + t3) / t3;

  double h_val = h(t);
  double h_prime_val = h_prime(t);
  double h_double_prime_val = h_double_prime(t);
  double tan_h = std::tan(h_val);
  double tan_h_squared = tan_h * tan_h;
  double one_plus_tan2 = 1.0 + tan_h_squared;

  double bracket_term = h_double_prime_val * tan_h +
                        h_prime_val * h_prime_val * (1.0 + 3.0 * tan_h_squared);
  double term2 = (one_plus_tan2 / 4.0) * bracket_term;

  return term1 + term2;
}

double TrigonometricKernel::psi_triple_prime(double t) const {
  // EXACT formula: ψ'''(t) = -3/t⁴ + [(1+tan²(h(t)))/4]k(t)
  validateInput(t, "psi_triple_prime");

  double t4 = t * t * t * t;
  double term1 = -3.0 / t4;

  double h_val = h(t);
  double tan_h = std::tan(h_val);
  double tan_h_squared = tan_h * tan_h;
  double one_plus_tan2 = 1.0 + tan_h_squared;

  double k_val = k(t);
  double term2 = (one_plus_tan2 / 4.0) * k_val;

  return term1 + term2;
}

double TrigonometricKernel::h(double t) const {
  // EXACT formula: h(t) = π(1-t)/(4t+2)
  return PI * (1.0 - t) / (4.0 * t + 2.0);
}

double TrigonometricKernel::h_prime(double t) const {
  // EXACT formula: h'(t) = -6π/(2+4t)²
  double denominator = 2.0 + 4.0 * t;
  return -6.0 * PI / (denominator * denominator);
}

double TrigonometricKernel::h_double_prime(double t) const {
  // EXACT formula: h''(t) = 48π/(2+4t)³
  double denominator = 2.0 + 4.0 * t;
  return 48.0 * PI / (denominator * denominator * denominator);
}

double TrigonometricKernel::h_triple_prime(double t) const {
  // EXACT formula: h'''(t) = -576π/(2+4t)⁴
  double denominator = 2.0 + 4.0 * t;
  double denom_4th = denominator * denominator * denominator * denominator;
  return -576.0 * PI / denom_4th;
}

double TrigonometricKernel::k(double t) const {
  // EXACT formula: k(t) = 3h'(t)h''(t)(1+3tan²(h(t))) +
  // 4h'(t)³tan(h(t))(2+3tan²(h(t))) + h'''(t)tan(h(t))
  double h_val = h(t);
  double h_prime_val = h_prime(t);
  double h_double_prime_val = h_double_prime(t);
  double h_triple_prime_val = h_triple_prime(t);

  double tan_h = std::tan(h_val);
  double tan_h_squared = tan_h * tan_h;

  double term1 =
      3.0 * h_prime_val * h_double_prime_val * (1.0 + 3.0 * tan_h_squared);
  double term2 = 4.0 * h_prime_val * h_prime_val * h_prime_val * tan_h *
                 (2.0 + 3.0 * tan_h_squared);
  double term3 = h_triple_prime_val * tan_h;

  return term1 + term2 + term3;
}

std::string TrigonometricKernel::getName() const {
  return "Trigonometric Kernel";
}

std::string TrigonometricKernel::getFormula() const {
  return "ψ(t) = (t-1)²/(2t) + (t-1)²/2 + tan²(h(t))/8, h(t) = π(1-t)/(4t+2)";
}

std::vector<double> TrigonometricKernel::getParameters() const {
  return {}; // No parameters for this kernel
}

void TrigonometricKernel::setParameters(const std::vector<double> &params) {
  // No parameters to set for this kernel
  (void)params; // Suppress unused parameter warning
}
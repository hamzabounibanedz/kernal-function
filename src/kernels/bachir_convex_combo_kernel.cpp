#include "bachir_convex_combo_kernel.h"
#include <cmath>

double BachirConvexComboKernel::psi(double t) const {
  validateInput(t, "BachirConvexComboKernel::psi");
  return t * t + m_m * std::log(t) + (1.0 - m_m) * std::exp(1.0 / t - 1.0);
}

double BachirConvexComboKernel::psi_prime(double t) const {
  validateInput(t, "BachirConvexComboKernel::psi_prime");
  double a = 2.0 * t + m_m / t;
  double b = (1.0 - m_m) * std::exp(1.0 / t - 1.0) * (-1.0 / (t * t));
  return a + b;
}

double BachirConvexComboKernel::psi_double_prime(double t) const {
  validateInput(t, "BachirConvexComboKernel::psi_double_prime");
  double term1 = 2.0 - m_m / (t * t);
  double e = std::exp(1.0 / t - 1.0);
  double de = e * (-1.0 / (t * t));
  double dde = e * (2.0 / (t * t * t)) + de * (-1.0 / (t * t));
  double term2 = (1.0 - m_m) * dde;
  return term1 + term2;
}

double BachirConvexComboKernel::psi_triple_prime(double t) const {
  validateInput(t, "BachirConvexComboKernel::psi_triple_prime");
  double h = std::max(1e-5, t * 1e-3);
  return (psi_double_prime(t + h) - psi_double_prime(t - h)) / (2.0 * h);
}






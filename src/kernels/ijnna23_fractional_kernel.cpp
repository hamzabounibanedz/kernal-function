#include "ijnna23_fractional_kernel.h"
#include <cmath>

double Ijna23FractionalKernel::psi(double t) const {
  validateInput(t, "Ijna23FractionalKernel::psi");
  // FORMULA: ψ(t) = (t²−1)/2 + (1/2)·ln(t) + 2/(t^q + 1)
  // Term 1: (t²−1)/2
  double term1 = (t * t - 1.0) / 2.0;
  // Term 2: + (1/2)·ln(t)
  double term2 = 0.5 * std::log(t);
  // Term 3: + 2/(t^q + 1)
  double denom = std::pow(t, m_q) + 1.0;
  double term3 = 2.0 / denom;
  return term1 + term2 + term3;
}

double Ijna23FractionalKernel::psi_prime(double t) const {
  validateInput(t, "Ijna23FractionalKernel::psi_prime");
  // DERIVATIVE: ψ'(t) = t + 1/(2t) − 2q·t^{q−1}/(t^q + 1)²
  double denom = std::pow(t, m_q) + 1.0;
  double d_term = -2.0 * m_q * std::pow(t, m_q - 1.0) / (denom * denom);
  return t + 0.5 / t + d_term;
}

double Ijna23FractionalKernel::psi_double_prime(double t) const {
  validateInput(t, "Ijna23FractionalKernel::psi_double_prime");
  // DERIVATIVE: ψ''(t) = 1 − 1/(2t²) − 2q(q−1)t^{q−2}/(t^q+1)² + 4q² t^{2q−2}/(t^q+1)³
  double denom = std::pow(t, m_q) + 1.0;
  double termA = -2.0 * (m_q * (m_q - 1.0)) * std::pow(t, m_q - 2.0) / (denom * denom);
  double termB = 4.0 * m_q * m_q * std::pow(t, 2 * m_q - 2.0) / (denom * denom * denom);
  return 1.0 - 1.0 / (2.0 * t * t) + termA + termB;
}

double Ijna23FractionalKernel::psi_triple_prime(double t) const {
  validateInput(t, "Ijna23FractionalKernel::psi_triple_prime");
  // DERIVATIVE: complex higher-order terms; use numerical differentiation
  double h = std::max(1e-5, t * 1e-3);
  return (psi_double_prime(t + h) - psi_double_prime(t - h)) / (2.0 * h);
}






#include "wu_hyperbolic_kernel.h"
#include <cmath>

static inline double coth(double x) { return std::cosh(x) / std::sinh(x); }
static inline double csch(double x) { return 1.0 / std::sinh(x); }

double WuHyperbolicKernel::psi(double t) const {
  validateInput(t, "WuHyperbolicKernel::psi");
  // FORMULA: ψ(t) = t² + a·p·coth(p t)·sinh⁻²(p·coth(1))
  // where a = sinh⁻²(1)·coth(p+1)
  // Step 1: coth(1)
  double coth_1 = std::cosh(1.0) / std::sinh(1.0);
  // Step 2: p·coth(1)
  double p_coth_1 = m_p * coth_1;
  // Step 3: sinh(p·coth(1))
  double sinh_p_coth_1 = std::sinh(p_coth_1);
  // Step 4: sinh⁻²(p·coth(1))
  double sinh_inv_sq_p_coth_1 = 1.0 / (sinh_p_coth_1 * sinh_p_coth_1);
  // Step 5: a = sinh⁻²(1)·coth(p+1)
  double sinh_inv_sq_1 = 1.0 / (std::sinh(1.0) * std::sinh(1.0));
  double coth_p_plus_1 = std::cosh(m_p + 1.0) / std::sinh(m_p + 1.0);
  double a = sinh_inv_sq_1 * coth_p_plus_1;

  double term = a * m_p * coth(m_p * t) * sinh_inv_sq_p_coth_1;
  return t * t + term;
}

double WuHyperbolicKernel::psi_prime(double t) const {
  validateInput(t, "WuHyperbolicKernel::psi_prime");
  // DERIVATIVE: ψ'(t) = 2t − a·p²·csch²(p·t)·sinh⁻²(p·coth(1))
  // Recompute constants (could be cached if class state allowed)
  double coth_1 = std::cosh(1.0) / std::sinh(1.0);
  double sinh_inv_sq_1 = 1.0 / (std::sinh(1.0) * std::sinh(1.0));
  double coth_p_plus_1 = std::cosh(m_p + 1.0) / std::sinh(m_p + 1.0);
  double a = sinh_inv_sq_1 * coth_p_plus_1;
  double sinh_inv_sq_p_coth_1 = 1.0 / std::pow(std::sinh(m_p * coth_1), 2.0);
  // d/dt coth(p t) = −p·csch²(p t)
  double term_prime = a * m_p * (-m_p * std::pow(csch(m_p * t), 2.0)) * sinh_inv_sq_p_coth_1;
  return 2.0 * t + term_prime;
}

double WuHyperbolicKernel::psi_double_prime(double t) const {
  validateInput(t, "WuHyperbolicKernel::psi_double_prime");
  // DERIVATIVE: ψ''(t) = 2 + 2a·p³·coth(p t)·csch²(p t)·sinh⁻²(p·coth(1))
  double coth_1 = std::cosh(1.0) / std::sinh(1.0);
  double sinh_inv_sq_1 = 1.0 / (std::sinh(1.0) * std::sinh(1.0));
  double coth_p_plus_1 = std::cosh(m_p + 1.0) / std::sinh(m_p + 1.0);
  double a = sinh_inv_sq_1 * coth_p_plus_1;
  double sinh_inv_sq_p_coth_1 = 1.0 / std::pow(std::sinh(m_p * coth_1), 2.0);
  double term = 2.0 * a * std::pow(m_p, 3.0) * coth(m_p * t) * std::pow(csch(m_p * t), 2.0) * sinh_inv_sq_p_coth_1;
  return 2.0 + term;
}

double WuHyperbolicKernel::psi_triple_prime(double t) const {
  validateInput(t, "WuHyperbolicKernel::psi_triple_prime");
  // Use numerical differentiation for third derivative due to complexity
  double h = std::max(1e-6, t * 1e-4);
  return (psi_double_prime(t + h) - psi_double_prime(t - h)) / (2.0 * h);
}



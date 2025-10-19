#include "touil_hyperbolic_kernel.h"
#include <cmath>
#include <stdexcept>

static inline double coth(double x) { return std::cosh(x) / std::sinh(x); }
static inline double csch(double x) { return 1.0 / std::sinh(x); }

// (no file-scope constants needed)

double TouilHyperbolicKernel::psi(double t) const {
  validateInput(t, "TouilHyperbolicKernel::psi");
  // FORMULA: ψ(t) = t² + (1/2)·sinh⁻²(1)·coth(t) − coth(1)
  const double sinh_inv_sq_1 = 1.0 / (std::sinh(1.0) * std::sinh(1.0));
  const double coth_1 = std::cosh(1.0) / std::sinh(1.0);
  // Term 1: t²
  double term1 = t * t;
  // Term 2: (1/2)·sinh⁻²(1)·coth(t)
  double term2 = 0.5 * sinh_inv_sq_1 * coth(t);
  // Term 3: −coth(1) constant
  double term3 = -coth_1;
  return term1 + term2 + term3;
}

// numerical derivative helpers removed as we have closed-form derivatives

double TouilHyperbolicKernel::psi_prime(double t) const {
  validateInput(t, "TouilHyperbolicKernel::psi_prime");
  // DERIVATIVE: ψ'(t) = 2t − (1/2)·sinh⁻²(1)·csch²(t)
  const double sinh_inv_sq_1 = 1.0 / (std::sinh(1.0) * std::sinh(1.0));
  return 2.0 * t - 0.5 * sinh_inv_sq_1 * (csch(t) * csch(t));
}

double TouilHyperbolicKernel::psi_double_prime(double t) const {
  validateInput(t, "TouilHyperbolicKernel::psi_double_prime");
  // DERIVATIVE: ψ''(t) = 2 + sinh⁻²(1)·coth(t)·csch²(t)
  const double sinh_inv_sq_1 = 1.0 / (std::sinh(1.0) * std::sinh(1.0));
  return 2.0 + sinh_inv_sq_1 * coth(t) * (csch(t) * csch(t));
}

double TouilHyperbolicKernel::psi_triple_prime(double t) const {
  validateInput(t, "TouilHyperbolicKernel::psi_triple_prime");
  // DERIVATIVE: ψ'''(t) = sinh⁻²(1)·[−csch⁴(t) − 2·coth²(t)·csch²(t)]
  const double sinh_inv_sq_1 = 1.0 / (std::sinh(1.0) * std::sinh(1.0));
  double csch2 = csch(t) * csch(t);
  double coth2 = coth(t) * coth(t);
  return sinh_inv_sq_1 * (-(csch2 * csch2) - 2.0 * coth2 * csch2);
}



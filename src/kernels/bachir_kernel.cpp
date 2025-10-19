#include "bachir_kernel.h"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>


BachirKernel::BachirKernel(double m) : m_m(std::max(0.0, std::min(1.0, m))) {}

double BachirKernel::psi(double t) const {
  // Domain restriction from paper: t > 1 (singularity at t=1)
  if (!(t > 1.0)) {
    throw std::invalid_argument("Invalid input t=" + std::to_string(t) +
                                " in psi. Domain requires t>1 for BachirKernel.");
  }
  validateInput(t, "psi");

  // FORMULA: φₘ(t) = (t²−1)/2 − m·ln(t) + (1−m)·(e^{1/(t−1)} − 1)
  // Term 1: (t²−1)/2
  double term1 = (t * t - 1.0) / 2.0;
  // Term 2: − m·ln(t)
  double term2 = -m_m * std::log(t);
  // Term 3: (1−m)·(e^{1/(t−1)} − 1)
  double term3 = (1.0 - m_m) * (std::exp(1.0 / (t - 1.0)) - 1.0);

  return term1 + term2 + term3;
}

double BachirKernel::psi_prime(double t) const {
  if (!(t > 1.0)) {
    throw std::invalid_argument("Invalid input t=" + std::to_string(t) +
                                " in psi_prime. Domain requires t>1 for BachirKernel.");
  }
  validateInput(t, "psi_prime");

  // DERIVATIVE: φₘ'(t) = t − m/t − (1−m)·e^{1/(t−1)}/(t−1)²
  double term1 = t;
  double term2 = -m_m / t;
  double t_minus_1 = t - 1.0;
  double exp_term = std::exp(1.0 / t_minus_1);
  double term3 = -(1.0 - m_m) * exp_term / (t_minus_1 * t_minus_1);
  return term1 + term2 + term3;
}

double BachirKernel::psi_double_prime(double t) const {
  if (!(t > 1.0)) {
    throw std::invalid_argument("Invalid input t=" + std::to_string(t) +
                                " in psi_double_prime. Domain requires t>1 for BachirKernel.");
  }
  validateInput(t, "psi_double_prime");

  // DERIVATIVE: φₘ''(t) = 1 + m/t² + (1−m)·e^{1/(t−1)}·[2/(t−1)³ − 1/(t−1)⁴]
  double term1 = 1.0;
  double term2 = m_m / (t * t);
  double tm1 = t - 1.0;
  double e = std::exp(1.0 / tm1);
  double bracket = 2.0 / (tm1 * tm1 * tm1) - 1.0 / (tm1 * tm1 * tm1 * tm1);
  double term3 = (1.0 - m_m) * e * bracket;
  return term1 + term2 + term3;
}

double BachirKernel::psi_triple_prime(double t) const {
  if (!(t > 1.0)) {
    throw std::invalid_argument("Invalid input t=" + std::to_string(t) +
                                " in psi_triple_prime. Domain requires t>1 for BachirKernel.");
  }
  validateInput(t, "psi_triple_prime");

  // DERIVATIVE: φₘ'''(t) = −2m/t³ + (1−m)·e^{1/(t−1)}·[higher order terms]
  // Use numerical differentiation for exponential term; keep exact −2m/t³.
  double term1 = -2.0 * m_m / (t * t * t);
  double h = std::max(1e-6, t * 1e-4);
  double num = (psi_double_prime(t + h) - psi_double_prime(t - h)) / (2.0 * h);
  return term1 + (num - term1);
}

std::string BachirKernel::getName() const { return "Bachir φₘ(t) Family"; }

std::string BachirKernel::getFormula() const {
  return "φₘ(t) = (t²-1)/2 - m ln(t) + (1-m)(e^{1/(t-1)}-1), 0≤m≤1";
}

std::vector<double> BachirKernel::getParameters() const { return {m_m}; }

void BachirKernel::setParameters(const std::vector<double> &params) {
  if (!params.empty()) {
    m_m = std::max(0.0, std::min(1.0, params[0])); // Clamp to [0,1]
  }
}
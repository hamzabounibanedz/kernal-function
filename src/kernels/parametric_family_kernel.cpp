#include "parametric_family_kernel.h"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>


ParametricFamilyKernel::ParametricFamilyKernel(double q)
    : m_q(std::max(1.0, q)) {}

double ParametricFamilyKernel::psi(double t) const {
  validateInput(t, "psi");

  // EXACT formula: ψ(t) = (t²-1)/2 - t^{-q}/(q²-q+1)*e^{1/(t-1)} +
  // (1-q)/(q²-q+1)
  double term1 = (t * t - 1.0) / 2.0;
  double term2 = -std::pow(t, -m_q) / denominator() * std::exp(1.0 / (t - 1.0));
  double term3 = (1.0 - m_q) / denominator();

  return term1 + term2 + term3;
}

double ParametricFamilyKernel::psi_prime(double t) const {
  validateInput(t, "psi_prime");

  // ψ'(t) derivative of the parametric family kernel
  double t_minus_1 = t - 1.0;
  double exp_term = std::exp(1.0 / t_minus_1);
  double t_neg_q = std::pow(t, -m_q);
  double denom = denominator();

  double term1 = t; // derivative of (t²-1)/2

  // Derivative of -t^{-q}/(q²-q+1)*e^{1/(t-1)}
  double d_t_neg_q = -m_q * std::pow(t, -m_q - 1.0);
  double d_exp = exp_term * (-1.0 / (t_minus_1 * t_minus_1));

  double term2 = -(d_t_neg_q * exp_term + t_neg_q * d_exp) / denom;

  return term1 + term2;
}

double ParametricFamilyKernel::psi_double_prime(double t) const {
  validateInput(t, "psi_double_prime");

  // Second derivative (complex calculation)
  double t_minus_1 = t - 1.0;
  double exp_term = std::exp(1.0 / t_minus_1);
  double t_neg_q = std::pow(t, -m_q);
  double denom = denominator();

  double term1 = 1.0; // second derivative of (t²-1)/2

  // Complex second derivative of the exponential term
  double d_t_neg_q = -m_q * std::pow(t, -m_q - 1.0);
  double d2_t_neg_q = m_q * (m_q + 1.0) * std::pow(t, -m_q - 2.0);

  double d_exp = exp_term * (-1.0 / (t_minus_1 * t_minus_1));
  double d2_exp = exp_term * (2.0 / std::pow(t_minus_1, 3.0) +
                              1.0 / std::pow(t_minus_1, 4.0));

  double term2 =
      -(d2_t_neg_q * exp_term + 2.0 * d_t_neg_q * d_exp + t_neg_q * d2_exp) /
      denom;

  return term1 + term2;
}

double ParametricFamilyKernel::psi_triple_prime(double t) const {
  validateInput(t, "psi_triple_prime");

  // Third derivative (very complex - simplified approximation)
  double h = 1e-8;
  return (psi_double_prime(t + h) - psi_double_prime(t - h)) / (2.0 * h);
}

std::string ParametricFamilyKernel::getName() const {
  return "Parametric Family Kernel";
}

std::string ParametricFamilyKernel::getFormula() const {
  return "ψ(t) = (t²-1)/2 - t^{-q}/(q²-q+1)*e^{1/(t-1)} + (1-q)/(q²-q+1), q≥1";
}

std::vector<double> ParametricFamilyKernel::getParameters() const {
  return {m_q};
}

void ParametricFamilyKernel::setParameters(const std::vector<double> &params) {
  if (!params.empty()) {
    m_q = std::max(1.0, params[0]);
  }
}
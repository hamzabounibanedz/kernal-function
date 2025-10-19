#include "log_exponential_kernel.h"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>


LogExponentialKernel::LogExponentialKernel(double q) : m_q(std::max(1.0, q)) {}

double LogExponentialKernel::psi(double t) const {
  validateInput(t, "psi");

  // FORMULA (paper-exact): ψ(t) = (t²−1)/2 − (1/2)·ln(t) + (e^{1/t^q−1} − 1)/(2q)
  // Term 1: quadratic barrier (t²-1)/2
  double term1 = (t * t - 1.0) / 2.0;
  // Term 2: logarithmic barrier −(1/2)·ln(t)
  double term2 = -0.5 * std::log(t);
  // Term 3: exponential term (e^{1/t^q-1} − 1)/(2q)
  double exp_term_val = std::exp(1.0 / std::pow(t, m_q) - 1.0);
  double term3 = (exp_term_val - 1.0) / (2.0 * m_q);

  return term1 + term2 + term3;
}

double LogExponentialKernel::psi_prime(double t) const {
  validateInput(t, "psi_prime");

  // DERIVATIVE (paper-exact): ψ'(t) = t − 1/(2t) − e^{1/t^q−1}/(2·t^{q+1})
  double term1 = t;            // derivative of (t²-1)/2
  double term2 = -0.5 / t;     // derivative of −(1/2)·ln(t)
  // Exponential derivative part: (1/(2q))·e^{...}·(−q·t^{−(q+1)})
  double exp_term_val = std::exp(1.0 / std::pow(t, m_q) - 1.0);
  double derivative_exponent = -m_q / std::pow(t, m_q + 1.0);
  double term3 = exp_term_val * derivative_exponent / (2.0 * m_q);
  return term1 + term2 + term3;
}

double LogExponentialKernel::psi_double_prime(double t) const {
  validateInput(t, "psi_double_prime");

  // DERIVATIVE (paper-exact): ψ''(t) = 1 + 1/(2t²) + (1/(2q))·e^{1/t^q−1}·[q(q+1)·t^{−(q+2)} + q²·t^{−(2q+2)}]
  double term1 = 1.0 + 1.0 / (2.0 * t * t);
  double exp_term_val = std::exp(1.0 / std::pow(t, m_q) - 1.0);
  double first_deriv_exp = -m_q / std::pow(t, m_q + 1.0);
  double second_deriv_exp = m_q * (m_q + 1.0) / std::pow(t, m_q + 2.0);
  double term2 = (exp_term_val * (first_deriv_exp * first_deriv_exp + second_deriv_exp)) / (2.0 * m_q);
  return term1 + term2;
}

double LogExponentialKernel::psi_triple_prime(double t) const {
  validateInput(t, "psi_triple_prime");

  // DERIVATIVE (paper-exact): ψ'''(t) = −1/t³ + (1/(2q))·d³/dt³[e^{1/t^q−1}]
  double term1 = -1.0 / (t * t * t);
  // Use numerical differentiation for third derivative of exponential term
  double h = std::max(1e-6, t * 1e-4);
  double d2p_ph = psi_double_prime(t + h);
  double d2p_mh = psi_double_prime(t - h);
  double term2 = (d2p_ph - d2p_mh) / (2.0 * h) - (-1.0 / (t * t * t));
  return term1 + term2;
}

std::string LogExponentialKernel::getName() const {
  return "Log + Exponential Kernel";
}

std::string LogExponentialKernel::getFormula() const {
  return "ψ(t) = (t²-1-ln t)/2 + (e^{1/t^q-1}-1)/(2q), q≥1";
}

std::vector<double> LogExponentialKernel::getParameters() const {
  return {m_q};
}

void LogExponentialKernel::setParameters(const std::vector<double> &params) {
  if (!params.empty()) {
    m_q = std::max(1.0, params[0]);
  }
}
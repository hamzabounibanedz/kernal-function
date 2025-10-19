#include "parameterized_log_kernel.h"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>


ParameterizedLogKernel::ParameterizedLogKernel(double q)
    : m_q(std::max(1.0, q)) {}

double ParameterizedLogKernel::psi(double t) const {
  validateInput(t, "psi");

  // FORMULA (paper-exact): ψ(t) = (t²-1)/2 − (1/2)·ln(t) + (e^{1/t^q−1} − 1)/(2q)
  // Term 1: quadratic barrier (t²-1)/2
  double term1 = (t * t - 1.0) / 2.0;
  // Term 2: logarithmic barrier −(1/2)·ln(t)
  double term2 = -0.5 * std::log(t);
  // Term 3: exponential kernel term (e^{1/t^q-1} - 1)/(2q)
  double term3 = (exp_term(t) - 1.0) / (2.0 * m_q);

  return term1 + term2 + term3;
}

double ParameterizedLogKernel::psi_prime(double t) const {
  validateInput(t, "psi_prime");

  // DERIVATIVE (paper-exact): ψ'(t) = t − 1/(2t) − e^{1/t^q−1}/(2·t^{q+1})
  // Term 1: derivative of (t²-1)/2 → t
  // Term 2: derivative of −(1/2)·ln(t) → −1/(2t)
  double term1 = t;
  double term2 = -0.5 / t;
  // Term 3: (1/(2q))·d/dt[e^{1/t^q-1}] = (1/(2q))·e^{...}·(−q·t^{−(q+1)}) = − e^{...}/(2·t^{q+1})
  double term3 = exp_term_prime(t) / (2.0 * m_q);

  return term1 + term2 + term3;
}

double ParameterizedLogKernel::psi_double_prime(double t) const {
  validateInput(t, "psi_double_prime");

  // DERIVATIVE (paper-exact): ψ''(t) = 1 + 1/(2t²) + (1/(2q))·e^{1/t^q−1}·[q(q+1)·t^{−(q+2)} + q²·t^{−(2q+2)}]
  // Term 1: derivative of t → 1
  // Term 2: derivative of −1/(2t) → +1/(2t²)
  double term1 = 1.0 + 1.0 / (2.0 * t * t);
  // Term 3: (1/(2q))·e^{...}·((d(1/t^q))² + d²(1/t^q))
  double term2 = exp_term_double_prime(t) / (2.0 * m_q);

  return term1 + term2;
}

double ParameterizedLogKernel::psi_triple_prime(double t) const {
  validateInput(t, "psi_triple_prime");

  // DERIVATIVE (paper-exact): ψ'''(t) = −1/t³ + (1/(2q))·d³/dt³[e^{1/t^q−1}]
  // Term 1: derivative of 1/(2t²) → −1/t³
  double term1 = -1.0 / (t * t * t);
  // Term 2: exact third derivative of exponential component
  double term2 = exp_term_triple_prime(t) / (2.0 * m_q);

  return term1 + term2;
}

double ParameterizedLogKernel::exp_term(double t) const {
  // e^{1/t^q-1}
  double exponent = 1.0 / std::pow(t, m_q) - 1.0;
  return std::exp(exponent);
}

double ParameterizedLogKernel::exp_term_prime(double t) const {
  // d/dt[e^{1/t^q-1}] = e^{1/t^q-1} * d/dt[1/t^q] = e^{1/t^q-1} * (-q/t^{q+1})
  double exp_val = exp_term(t);
  double derivative_exponent = -m_q / std::pow(t, m_q + 1.0);
  return exp_val * derivative_exponent;
}

double ParameterizedLogKernel::exp_term_double_prime(double t) const {
  // Second derivative of e^{1/t^q-1}
  double exp_val = exp_term(t);
  double first_deriv_exp = -m_q / std::pow(t, m_q + 1.0);
  double second_deriv_exp = m_q * (m_q + 1.0) / std::pow(t, m_q + 2.0);

  return exp_val * (first_deriv_exp * first_deriv_exp + second_deriv_exp);
}

double ParameterizedLogKernel::exp_term_triple_prime(double t) const {
  // Third derivative of e^{1/t^q-1}
  double exp_val = exp_term(t);
  double first_deriv_exp = -m_q / std::pow(t, m_q + 1.0);
  double second_deriv_exp = m_q * (m_q + 1.0) / std::pow(t, m_q + 2.0);
  double third_deriv_exp =
      -m_q * (m_q + 1.0) * (m_q + 2.0) / std::pow(t, m_q + 3.0);

  return exp_val * (first_deriv_exp * first_deriv_exp * first_deriv_exp +
                    3.0 * first_deriv_exp * second_deriv_exp + third_deriv_exp);
}

std::string ParameterizedLogKernel::getName() const {
  return "Parameterized Log Kernel";
}

std::string ParameterizedLogKernel::getFormula() const {
  return "ψ(t) = (t²-1-ln t)/2 + (e^{1/t^q-1}-1)/(2q), q≥1";
}

std::vector<double> ParameterizedLogKernel::getParameters() const {
  return {m_q};
}

void ParameterizedLogKernel::setParameters(const std::vector<double> &params) {
  if (!params.empty()) {
    m_q = std::max(1.0, params[0]);
  }
}
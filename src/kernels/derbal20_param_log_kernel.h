#ifndef DERBAL20_PARAM_LOG_KERNEL_H
#define DERBAL20_PARAM_LOG_KERNEL_H

#include "kernel_base.h"
#include <algorithm>
#include <cmath>

// Derbal & Kebbiche (2020): ψ(t) = t^2 + log t + 2(e-1)/(t^q + 1)
class Derbal20ParamLogKernel : public KernelBase {
public:
  explicit Derbal20ParamLogKernel(double q = 2.0) : m_q(std::max(1.0, q)) {}
  ~Derbal20ParamLogKernel() override = default;

  double psi(double t) const override {
    validateInput(t, "Derbal20ParamLogKernel::psi");
    const double denom = std::pow(t, m_q) + 1.0;
    return t * t + std::log(t) + 2.0 * (E_MINUS_1) / denom;
  }

  double psi_prime(double t) const override {
    validateInput(t, "Derbal20ParamLogKernel::psi_prime");
    const double denom = std::pow(t, m_q) + 1.0;
    const double d_denom = m_q * std::pow(t, m_q - 1.0);
    const double d_term = -2.0 * (E_MINUS_1) * d_denom / (denom * denom);
    return 2.0 * t + 1.0 / t + d_term;
  }

  double psi_double_prime(double t) const override {
    validateInput(t, "Derbal20ParamLogKernel::psi_double_prime");
    const double denom = std::pow(t, m_q) + 1.0;
    const double d_denom = m_q * std::pow(t, m_q - 1.0);
    const double dd_denom = m_q * (m_q - 1.0) * std::pow(t, m_q - 2.0);
    const double dd_term = -2.0 * (E_MINUS_1) *
                           (dd_denom * denom - 2.0 * d_denom * d_denom) /
                           (denom * denom * denom);
    return 2.0 - 1.0 / (t * t) + dd_term;
  }

  double psi_triple_prime(double t) const override {
    validateInput(t, "Derbal20ParamLogKernel::psi_triple_prime");
    double h = std::max(1e-5, t * 1e-3);
    return (psi_double_prime(t + h) - psi_double_prime(t - h)) / (2.0 * h);
  }

  std::string getName() const override { return "Parameterized-Log [Derbal20]"; }
  std::string getFormula() const override {
    return "ψ(t)=t^2 + log t + 2(e-1)/(t^q+1)";
  }
  std::vector<double> getParameters() const override { return {m_q}; }
  void setParameters(const std::vector<double> &params) override {
    if (!params.empty()) {
      // clamp to UI range [1, 10]
      double v = params[0];
      if (v < 1.0) v = 1.0;
      if (v > 10.0) v = 10.0;
      m_q = v;
    }
  }

private:
  double m_q; // q>=1
  static constexpr double E_MINUS_1 = 2.71828182845904523536 - 1.0;
};

#endif // DERBAL20_PARAM_LOG_KERNEL_H



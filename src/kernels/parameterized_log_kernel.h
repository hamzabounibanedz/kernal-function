#ifndef PARAMETERIZED_LOG_KERNEL_H
#define PARAMETERIZED_LOG_KERNEL_H

#include "kernel_base.h"
#include <cmath>

/**
 * @brief Parameterized Log Kernel Function for Algorithm 3
 *
 * Kernel function from Equation (1):
 * ψ(t) = (t²-1-ln t)/2 + (e^{1/t^q-1}-1)/(2q), t>0, q≥1
 *
 * This kernel is used in PD IPM for SDO with Log-Kernel.
 */
class ParameterizedLogKernel : public KernelBase {
public:
  explicit ParameterizedLogKernel(double q = 1.0);
  virtual ~ParameterizedLogKernel() = default;

  // Core kernel function evaluation
  double psi(double t) const override;
  double psi_prime(double t) const override;
  double psi_double_prime(double t) const override;
  double psi_triple_prime(double t) const override;

  // Metadata
  std::string getName() const override;
  std::string getFormula() const override;
  std::vector<double> getParameters() const override;
  void setParameters(const std::vector<double> &params) override;

private:
  double m_q; // Parameter q ≥ 1

  /**
   * @brief Helper function for exponential term
   * @param t Input value
   * @return e^{1/t^q-1}
   */
  double exp_term(double t) const;

  /**
   * @brief First derivative of exponential term
   * @param t Input value
   * @return d/dt[e^{1/t^q-1}]
   */
  double exp_term_prime(double t) const;

  /**
   * @brief Second derivative of exponential term
   * @param t Input value
   * @return d²/dt²[e^{1/t^q-1}]
   */
  double exp_term_double_prime(double t) const;

  /**
   * @brief Third derivative of exponential term
   * @param t Input value
   * @return d³/dt³[e^{1/t^q-1}]
   */
  double exp_term_triple_prime(double t) const;
};

#endif // PARAMETERIZED_LOG_KERNEL_H
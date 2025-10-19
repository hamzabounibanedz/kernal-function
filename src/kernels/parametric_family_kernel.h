#ifndef PARAMETRIC_FAMILY_KERNEL_H
#define PARAMETRIC_FAMILY_KERNEL_H

#include "kernel_base.h"

/**
 * @brief Parametric Family Kernel Function for Algorithm 4 (CQSDO)
 *
 * Kernel function from Equation (1):
 * ψ(t) = (t²-1)/2 - t^{-q}/(q²-q+1)*e^{1/(t-1)} + (1-q)/(q²-q+1), q≥1
 */
class ParametricFamilyKernel : public KernelBase {
public:
  explicit ParametricFamilyKernel(double q = 1.0);
  virtual ~ParametricFamilyKernel() = default;

  double psi(double t) const override;
  double psi_prime(double t) const override;
  double psi_double_prime(double t) const override;
  double psi_triple_prime(double t) const override;

  std::string getName() const override;
  std::string getFormula() const override;
  std::vector<double> getParameters() const override;
  void setParameters(const std::vector<double> &params) override;

private:
  double m_q;
  double denominator() const { return m_q * m_q - m_q + 1.0; }
};

#endif // PARAMETRIC_FAMILY_KERNEL_H
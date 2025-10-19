#ifndef LOG_EXPONENTIAL_KERNEL_H
#define LOG_EXPONENTIAL_KERNEL_H

#include "kernel_base.h"

/**
 * @brief Log + Exponential Kernel Function for Algorithm 5 (Derbal & Kebbiche)
 *
 * ψ(t) = (t²-1-ln t)/2 + (e^{1/t^q-1}-1)/(2q), q≥1
 * Same as ParameterizedLogKernel but used in different algorithm context
 */
class LogExponentialKernel : public KernelBase {
public:
  explicit LogExponentialKernel(double q = 1.0);
  virtual ~LogExponentialKernel() = default;

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
};

#endif // LOG_EXPONENTIAL_KERNEL_H
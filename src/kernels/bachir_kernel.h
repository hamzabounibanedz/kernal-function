#ifndef BACHIR_KERNEL_H
#define BACHIR_KERNEL_H

#include "kernel_base.h"

/**
 * @brief Bachir φₘ(t) Family Kernel Function for Algorithm 6
 *
 * φₘ(t) = (t²-1)/2 - m ln t + (1-m)(e^{1/(t-1)}-1), 0≤m≤1, t>0
 * Recommended m ≈ 0.2–0.5
 */
class BachirKernel : public KernelBase {
public:
  explicit BachirKernel(double m = 0.3);
  virtual ~BachirKernel() = default;

  double psi(double t) const override;
  double psi_prime(double t) const override;
  double psi_double_prime(double t) const override;
  double psi_triple_prime(double t) const override;

  std::string getName() const override;
  std::string getFormula() const override;
  std::vector<double> getParameters() const override;
  void setParameters(const std::vector<double> &params) override;

private:
  double m_m; // Parameter m ∈ [0,1]
};

#endif // BACHIR_KERNEL_H
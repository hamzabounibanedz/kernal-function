#ifndef WU_HYPERBOLIC_KERNEL_H
#define WU_HYPERBOLIC_KERNEL_H

#include "kernel_base.h"

// Wu & Zhang (2025) hyperbolic kernel: φ(t)=t^2 + a * p*coth(p t) * sinh^{-2}(p*coth(1))
// with p>0 and a= sinh^{-2}(1)*coth(p+1)
class WuHyperbolicKernel : public KernelBase {
public:
  explicit WuHyperbolicKernel(double p = 1.0) : m_p(p) {}
  ~WuHyperbolicKernel() override = default;

  double psi(double t) const override;
  double psi_prime(double t) const override;
  double psi_double_prime(double t) const override;
  double psi_triple_prime(double t) const override;

  std::string getName() const override { return "Hyperbolic [Wu25]"; }
  std::string getFormula() const override {
    return "φ(t)=t^2 + a*p*coth(p t)*sinh^{-2}(p*coth(1)), a=sinh^{-2}(1)*coth(p+1)";
  }
  std::vector<double> getParameters() const override { return {m_p}; }
  void setParameters(const std::vector<double> &params) override {
    if (!params.empty()) {
      // clamp to UI range [0.5, 5]
      double v = params[0];
      if (v < 0.5) v = 0.5;
      if (v > 5.0) v = 5.0;
      m_p = v;
    }
  }

private:
  double m_p; // p>0
};

#endif // WU_HYPERBOLIC_KERNEL_H



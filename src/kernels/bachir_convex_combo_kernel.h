#ifndef BACHIR_CONVEX_COMBO_KERNEL_H
#define BACHIR_CONVEX_COMBO_KERNEL_H

#include "kernel_base.h"

// Bachir convex combination: φ(t)=t^2 + m*log t + (1-m)*exp(1/t - 1), 0<m<1
class BachirConvexComboKernel : public KernelBase {
public:
  explicit BachirConvexComboKernel(double m = 0.5) : m_m(m) {}
  ~BachirConvexComboKernel() override = default;

  double psi(double t) const override;
  double psi_prime(double t) const override;
  double psi_double_prime(double t) const override;
  double psi_triple_prime(double t) const override;

  std::string getName() const override { return "Convex Combo [Bachir25]"; }
  std::string getFormula() const override {
    return "φ(t)=t^2 + m*log t + (1-m)*exp(1/t - 1)";
  }
  std::vector<double> getParameters() const override { return {m_m}; }
  void setParameters(const std::vector<double> &params) override {
    if (!params.empty()) {
      // clamp to UI range [0.01, 0.99]
      double v = params[0];
      if (v < 0.01) v = 0.01;
      if (v > 0.99) v = 0.99;
      m_m = v;
    }
  }

private:
  double m_m; // in (0,1)
};

#endif // BACHIR_CONVEX_COMBO_KERNEL_H



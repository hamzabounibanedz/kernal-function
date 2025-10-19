#ifndef TOUIL_HYPERBOLIC_KERNEL_H
#define TOUIL_HYPERBOLIC_KERNEL_H

#include "kernel_base.h"
#include <string>

// Touil & Chikouche (2022) Hyperbolic barrier kernel (no parameters)
// φ(t) = t^2 + 0.5*sinh^{-2}(1)·coth(t) - coth(1), t>0
class TouilHyperbolicKernel : public KernelBase {
public:
  TouilHyperbolicKernel() = default;
  ~TouilHyperbolicKernel() override = default;

  double psi(double t) const override;
  double psi_prime(double t) const override;
  double psi_double_prime(double t) const override;
  double psi_triple_prime(double t) const override;

  std::string getName() const override { return "Novel-Hyperbolic [Touil22]"; }
  std::string getFormula() const override {
    return "φ(t)=t^2 + 0.5*sinh^{-2}(1)·coth(t) - coth(1)";
  }
  std::vector<double> getParameters() const override { return {}; }
  void setParameters(const std::vector<double> &) override {}

private:
};

#endif // TOUIL_HYPERBOLIC_KERNEL_H






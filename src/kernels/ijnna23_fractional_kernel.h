#ifndef IJNNA23_FRACTIONAL_KERNEL_H
#define IJNNA23_FRACTIONAL_KERNEL_H

#include "kernel_base.h"

// IJNAA 2023 fractional barrier variant: φ(t)=t^2 + log t + 2/(t^q+1)
class Ijna23FractionalKernel : public KernelBase {
public:
  explicit Ijna23FractionalKernel(double q = 2.0) : m_q(q) {}
  ~Ijna23FractionalKernel() override = default;

  double psi(double t) const override;
  double psi_prime(double t) const override;
  double psi_double_prime(double t) const override;
  double psi_triple_prime(double t) const override;

  std::string getName() const override { return "Fractional-Barrier [IJNAA23]"; }
  std::string getFormula() const override {
    return "φ(t)=t^2 + log t + 2/(t^q+1)";
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
  double m_q; // q>1
};

#endif // IJNNA23_FRACTIONAL_KERNEL_H



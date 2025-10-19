#ifndef TRIGONOMETRIC_KERNEL_H
#define TRIGONOMETRIC_KERNEL_H

#include "kernel_base.h"
#include <cmath>

/**
 * @brief Trigonometric Kernel Function for Algorithm 1
 *
 * EXACT implementation from requirements:
 * ψ(t) = (t-1)²/(2t) + (t-1)²/2 + tan²(h(t))/8
 * where h(t) = π(1-t)/(4t+2)
 *
 * First three derivatives:
 * ψ'(t) = (2t³-t²-1)/(2t²) + [h'(t)tan(h(t))(1+tan²(h(t)))]/4
 * ψ''(t) = (1+t³)/t³ + [(1+tan²(h(t)))/4][h''(t)tan(h(t)) +
 * h'(t)²(1+3tan²(h(t)))] ψ'''(t) = -3/t⁴ + [(1+tan²(h(t)))/4]k(t)
 *
 * where:
 * h'(t) = -6π/(2+4t)²
 * h''(t) = 48π/(2+4t)³
 * h'''(t) = -576π/(2+4t)⁴
 * k(t) = 3h'(t)h''(t)(1+3tan²(h(t))) + 4h'(t)³tan(h(t))(2+3tan²(h(t))) +
 * h'''(t)tan(h(t))
 */
class TrigonometricKernel : public KernelBase {
public:
  TrigonometricKernel() = default;
  virtual ~TrigonometricKernel() = default;

  // Core kernel function evaluation - EXACT formulas
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
  // Helper function h(t) = π(1-t)/(4t+2)
  double h(double t) const;

  // Derivatives of h(t) - EXACT formulas
  double h_prime(double t) const;        // h'(t) = -6π/(2+4t)²
  double h_double_prime(double t) const; // h''(t) = 48π/(2+4t)³
  double h_triple_prime(double t) const; // h'''(t) = -576π/(2+4t)⁴

  // Helper function k(t) for third derivative
  double k(double t) const;

  static constexpr double PI = 3.14159265358979323846;
};

#endif // TRIGONOMETRIC_KERNEL_H
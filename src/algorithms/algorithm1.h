#ifndef ALGORITHM1_H
#define ALGORITHM1_H

#include "../kernels/trigonometric_kernel.h"
#include "algorithm_base.h"
#include <chrono>

/**
 * Algorithm 1: Generic Primal-Dual Interior-Point Algorithm for LO
 *
 * From the first research paper - Figure 1
 *
 * Input: A ∈ ℝᵐ×ⁿ, b ∈ ℝᵐ, c ∈ ℝⁿ,
 *        strictly feasible (x⁰, y⁰, s⁰) > 0,
 *        parameters θ ∈ (0,1), tolerance ε > 0.
 *
 * Initialize: μ ← (x⁰)ᵀ s⁰ / n.
 *
 * While μ > ε do
 *   1. Let v = √(x ∘ s / μ).                 # element-wise
 *   2. Solve for (Δx, Δy, Δs):
 *        A Δx = 0,
 *        Aᵀ Δy + Δs = 0,
 *        dx + ds = –∇Ψ(v),
 *      where dx = v ∘ Δx / x, ds = v ∘ Δs / s.
 *   3. Choose step-size α > 0 (e.g. via line search or theoretical bound).
 *   4. Update x ← x + α Δx, y ← y + α Δy, s ← s + α Δs.
 *   5. Update μ ← (1 – θ) μ.                # large-update if θ constant
 * End while
 *
 * Output: Approximate primal-dual solution (x, y, s).
 */
class Algorithm1 : public AlgorithmBase {
public:
  Algorithm1();
  virtual ~Algorithm1() = default;

  AlgorithmResult solve(const TestCase &test_case) override;
  std::string getName() const override;
  std::string getDescription() const override;
  std::vector<TestCase> getTestCases() const override;

private:
  // Search direction structure
  struct SearchDirection {
    std::vector<double> delta_x;
    std::vector<double> delta_y;
    std::vector<double> delta_s;
  };

  // Current iterate structure
  struct CurrentIterate {
    std::vector<double> x; // Primal variables
    std::vector<double> y; // Dual variables
    std::vector<double> s; // Dual slack variables
    double mu;             // Barrier parameter
  };

  /**
   * Initialize the algorithm with starting point
   */
  CurrentIterate initialize(const TestCase &test_case) const;

  /**
   * Compute scaling vector v = √(x ∘ s / μ)
   */
  std::vector<double> computeScalingVector(const std::vector<double> &x,
                                           const std::vector<double> &s,
                                           double mu) const;

  /**
   * Compute kernel gradient -∇Ψ(v)
   */
  std::vector<double> computeKernelGradient(const std::vector<double> &v) const;

  /**
   * Solve the linear system for search direction
   * A Δx = 0
   * Aᵀ Δy + Δs = 0
   * dx + ds = –∇Ψ(v)
   * where dx = v ∘ Δx / x, ds = v ∘ Δs / s
   */
  SearchDirection solveLinearSystem(const TestCase &test_case,
                                    const CurrentIterate &current,
                                    const std::vector<double> &v,
                                    const std::vector<double> &rhs_psi) const;

  /**
   * Compute step size α using line search
   */
  double computeStepSize(const CurrentIterate &current,
                         const SearchDirection &direction) const;

  /**
   * Update the current iterate
   */
  void updateIterate(CurrentIterate &current, const SearchDirection &direction,
                     double alpha) const;

  /**
   * Check convergence criterion
   */
  bool checkConvergence(double mu, double epsilon) const;

  /**
   * Compute objective function value
   */
  double computeObjectiveValue(const std::vector<double> &x,
                               const std::vector<std::vector<double>> &C) const;

  /**
   * Compute duality gap
   */
  double computeDualityGap(const std::vector<double> &x,
                           const std::vector<double> &s) const;

  /**
   * Validate that current point is strictly feasible
   */
  bool isStrictlyFeasible(const std::vector<double> &x,
                          const std::vector<double> &s) const;

  std::shared_ptr<TrigonometricKernel> m_trigKernel;
};

#endif // ALGORITHM1_H
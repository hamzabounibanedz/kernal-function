#ifndef ALGORITHM4_H
#define ALGORITHM4_H

#include "../kernels/parametric_family_kernel.h"
#include "algorithm_base.h"
#include <chrono>

/**
 * Algorithm 4: Generic Primal-Dual Algorithm for CQSDO
 *
 * From the fourth research paper
 *
 * Inputs:
 *   • Threshold parameter τ (τ ≥ 1)
 *   • Accuracy parameter ε > 0
 *   • Barrier update parameter θ ∈ (0, 1)
 *   • Strictly feasible starting point (X₀,y₀,Z₀) with μ₀=1 and Ψ(X₀,Z₀;μ₀)≤τ
 *
 * Procedure:
 *   X ← X₀; y ← y₀; Z ← Z₀; μ ← μ₀
 *   while nμ ≥ ε do
 *     μ ← (1–θ)·μ
 *     while Ψ(X,Z;μ) > τ do
 *       1. Solve the scaled NT-system with kernel-based right-hand side
 * (equation (9)) to get (ΔX,Δy,ΔZ).
 *       2. Determine step-size α ∈ (0,1).
 *       3. Update (X,y,Z) ← (X,y,Z) + α (ΔX,Δy,ΔZ).
 *     end
 *   end
 */
class Algorithm4 : public AlgorithmBase {
public:
  Algorithm4();
  virtual ~Algorithm4() = default;

  AlgorithmResult solve(const TestCase &test_case) override;
  std::string getName() const override;
  std::string getDescription() const override;
  std::vector<TestCase> getTestCases() const override;

private:
  // Search direction structure
  struct SearchDirection {
    std::vector<std::vector<double>> delta_X;
    std::vector<double> delta_y;
    std::vector<std::vector<double>> delta_Z;
  };

  // Current iterate structure
  struct CurrentIterate {
    std::vector<std::vector<double>> X; // Primal variables
    std::vector<double> y;              // Dual variables
    std::vector<std::vector<double>>
        Z;     // Dual slack variables (Z instead of S for CQSDO)
    double mu; // Barrier parameter
  };

  /**
   * Initialize the algorithm with starting point
   */
  CurrentIterate initialize(const TestCase &test_case) const;

  /**
   * Compute proximity measure Ψ(X,Z;μ)
   */
  double computeProximityMeasure(const CurrentIterate &current) const;

  /**
   * Compute kernel gradient -ψ′(V)
   */
  std::vector<std::vector<double>>
  computeKernelGradient(const std::vector<std::vector<double>> &V) const;

  /**
   * Compute scaling matrix V
   */
  std::vector<std::vector<double>>
  computeScalingMatrix(const CurrentIterate &current) const;

  /**
   * Solve the scaled NT-system for search direction
   */
  SearchDirection
  solveLinearSystem(const TestCase &test_case, const CurrentIterate &current,
                    const std::vector<std::vector<double>> &rhs_psi) const;

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
  bool checkConvergence(double mu, double epsilon, int n) const;

  /**
   * Check proximity criterion
   */
  bool checkProximity(double proximity, double tau) const;

  /**
   * Compute objective function value
   */
  double computeObjectiveValue(const std::vector<std::vector<double>> &X,
                               const std::vector<std::vector<double>> &C) const;

  /**
   * Compute duality gap
   */
  double computeDualityGap(const std::vector<std::vector<double>> &X,
                           const std::vector<std::vector<double>> &Z) const;

  /**
   * Validate that current point is strictly feasible
   */
  bool isStrictlyFeasible(const std::vector<std::vector<double>> &X,
                          const std::vector<std::vector<double>> &Z) const;

  /**
   * Matrix trace operation
   */
  double matrixTrace(const std::vector<std::vector<double>> &matrix) const;

  /**
   * Matrix addition
   */
  std::vector<std::vector<double>>
  matrixAdd(const std::vector<std::vector<double>> &A,
            const std::vector<std::vector<double>> &B) const;

  /**
   * Matrix multiplication
   */
  std::vector<std::vector<double>>
  matrixMultiply(const std::vector<std::vector<double>> &A,
                 const std::vector<std::vector<double>> &B) const;

  std::shared_ptr<ParametricFamilyKernel> m_familyKernel;
};

#endif // ALGORITHM4_H
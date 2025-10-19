#ifndef ALGORITHM5_H
#define ALGORITHM5_H

#include "../kernels/bachir_kernel.h"
#include "algorithm_base.h"
#include <chrono>

/**
 * Algorithm 5: Generic Primal-Dual Interior-Point Algorithm for SDO
 *
 * From the fifth research paper - Section 2.2
 *
 * Input:
 *   • kernel function ψ(t)
 *   • threshold τ > 1
 *   • accuracy ε > 0
 *   • barrier-update θ ∈ (0,1)
 *
 * Initialize:
 *   X := I, S := I, μ := 1, V := I
 *
 * while n·μ ≥ ε do             ← outer iterations
 *   μ ← (1–θ)·μ
 *   V ← V / √(1–θ)
 *
 *   while Ψ(V) > τ do          ← inner iterations
 *     • Solve for (ΔX, Δy, ΔS) via the system
 *         Ai·ΔX = 0  (i=1…m)
 *         ∑ᵢ Δyᵢ Ai + ΔS = 0
 *         ΔX + ΔS = –ψ′(V)
 *     • Line-search to pick step α
 *     • X ← X + α ΔX,
 *       S ← S + α ΔS,
 *       y ← y + α Δy
 *     • V ← (1/√μ)·D⁻¹ X D⁻¹   (or equivalently using S)
 *   end inner
 * end outer
 */
class Algorithm5 : public AlgorithmBase {
public:
  Algorithm5();
  virtual ~Algorithm5() = default;

  AlgorithmResult solve(const TestCase &test_case) override;
  std::string getName() const override;
  std::string getDescription() const override;
  std::vector<TestCase> getTestCases() const override;

private:
  // Search direction structure
  struct SearchDirection {
    std::vector<std::vector<double>> delta_X;
    std::vector<double> delta_y;
    std::vector<std::vector<double>> delta_S;
  };

  // Current iterate structure
  struct CurrentIterate {
    std::vector<std::vector<double>> X; // Primal variables
    std::vector<double> y;              // Dual variables
    std::vector<std::vector<double>> S; // Dual slack variables
    double mu;                          // Barrier parameter
    std::vector<std::vector<double>> V; // Scaling matrix
  };

  /**
   * Initialize the algorithm with starting point
   */
  CurrentIterate initialize(const TestCase &test_case) const;

  /**
   * Compute proximity measure Ψ(V)
   */
  double
  computeProximityMeasure(const std::vector<std::vector<double>> &V) const;

  /**
   * Update scaling matrix V ← V / √(1–θ)
   */
  void updateScalingMatrix(CurrentIterate &current,
                           const TestCase &test_case) const;

  /**
   * Compute kernel gradient -ψ′(V)
   */
  std::vector<std::vector<double>>
  computeKernelGradient(const std::vector<std::vector<double>> &V) const;

  /**
   * Solve the linear system for search direction
   * Ai·ΔX = 0  (i=1…m)
   * ∑ᵢ Δyᵢ Ai + ΔS = 0
   * ΔX + ΔS = –ψ′(V)
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
   * Update V ← (1/√μ)·D⁻¹ X D⁻¹
   */
  void updateVMatrix(CurrentIterate &current) const;

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
                           const std::vector<std::vector<double>> &S) const;

  /**
   * Validate that current point is strictly feasible
   */
  bool isStrictlyFeasible(const std::vector<std::vector<double>> &X,
                          const std::vector<std::vector<double>> &S) const;

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

  std::shared_ptr<BachirKernel> m_bachirKernel;
};

#endif // ALGORITHM5_H
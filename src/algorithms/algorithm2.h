#ifndef ALGORITHM2_H
#define ALGORITHM2_H

#include "../kernels/exponential_parametric_kernel.h"
#include "algorithm_base.h"
#include <chrono>

/**
 * Algorithm 2: Generic Primal-Dual IPM for SDO
 *
 * From the second research paper - Algorithm 1
 *
 * Input:
 *   • threshold τ > 0
 *   • accuracy ε > 0
 *   • barrier-update parameter 0 < θ < 1
 *   • strictly feasible (X₀,y₀,S₀) with μ₀=1 and Ψ(X₀,S₀,μ₀)≤τ
 *
 * Begin
 *   X←X₀; S←S₀; μ←μ₀;
 *   while n·μ > ε do
 *     μ ← (1–θ)·μ;
 *     while Ψ(X,S,μ) > τ do
 *       solve
 *         Āᵢ•ΔX=0,  ∑ᵢΔyᵢĀᵢ+ΔS=0,
 *         ΔX+ΔS = –ψ′(V)
 *       for (ΔX,Δy,ΔS) and recover (ΔX,ΔS) via scaling;
 *       choose step-size α∈(0,1];
 *       X ← X + α·ΔX;
 *       y ← y + α·Δy;
 *       S ← S + α·ΔS;
 *       V ← (μ⁻¹·D⁻¹XSD)¹ᐟ²;
 *     end while
 *   end while
 * End
 */
class Algorithm2 : public AlgorithmBase {
public:
  Algorithm2();
  virtual ~Algorithm2() = default;

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
  };

  /**
   * Initialize the algorithm with starting point
   */
  CurrentIterate initialize(const TestCase &test_case) const;

  /**
   * Compute proximity measure Ψ(X,S,μ)
   */
  double computeProximityMeasure(const CurrentIterate &current) const;

  /**
   * Compute scaling matrix V = (μ⁻¹·D⁻¹XSD)¹ᐟ²
   */
  std::vector<std::vector<double>>
  computeScalingMatrix(const CurrentIterate &current) const;

  /**
   * Compute kernel gradient -ψ′(V)
   */
  std::vector<std::vector<double>>
  computeKernelGradient(const std::vector<std::vector<double>> &V) const;

  /**
   * Solve the linear system for search direction
   * Āᵢ•ΔX=0, ∑ᵢΔyᵢĀᵢ+ΔS=0, ΔX+ΔS = –ψ′(V)
   */
  SearchDirection
  solveLinearSystem(const TestCase &test_case, const CurrentIterate &current,
                    const std::vector<std::vector<double>> &V,
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

  std::shared_ptr<ExponentialParametricKernel> m_expKernel;
};

#endif // ALGORITHM2_H
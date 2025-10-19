// IpmEngine public API and types
#ifndef IPM_ENGINE_H
#define IPM_ENGINE_H

#include <atomic>
#include <vector>
// CPU time helpers
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <ctime>
#endif

#include <Eigen/Dense>

#include "../algorithms/algorithm_base.h"
#include "../kernels/kernel_base.h"
#include "../utils/test_data.h"

// Settings controlling the IPM solve
struct IpmSettings {
  double epsilon = 1e-8;          // convergence tolerance
  double theta = 0.1;             // barrier update factor (0,1)
  double tau = 1.0;               // proximity threshold
  int maxOuter = 200;             // max outer iterations
  int maxInner = 50;              // max inner iterations per μ
  std::atomic<bool> *cancel = nullptr; // optional cancellation flag
};

class IpmEngine {
public:
  struct Iterate {
    Eigen::MatrixXd X;            // primal variable
    std::vector<double> y;        // dual multipliers
    Eigen::MatrixXd S;            // dual slack
    double mu = 1.0;              // barrier parameter
  };

  struct Direction {
    Eigen::MatrixXd dX;
    std::vector<double> dy;
    Eigen::MatrixXd dS;
  };

  // Solve SDO using kernel-based NT directions
  AlgorithmResult solve(const TestCase &tc, KernelBase &kernel,
                        const IpmSettings &settings);

private:
  Iterate initialize(const TestCase &tc) const;

  Eigen::MatrixXd computeNTScalingV(const Iterate &it) const;
  double computeProximityTracePsi(const Eigen::MatrixXd &V,
                                  KernelBase &kernel) const;
  Eigen::MatrixXd kernelGradient(const Eigen::MatrixXd &V,
                                 KernelBase &kernel) const;

  Direction solveNewtonSimplified(const TestCase &tc, const Iterate &it,
                                  const Eigen::MatrixXd &rhs) const;
  Direction solveNewtonNT(const TestCase &tc, const Iterate &it,
                          const Eigen::MatrixXd &R,
                          Eigen::MatrixXd &D_out) const;

  static Eigen::MatrixXd symmetrize(const Eigen::MatrixXd &M);
  static Eigen::MatrixXd sqrtmSPD(const Eigen::MatrixXd &M);
  static Eigen::MatrixXd invsqrtmSPD(const Eigen::MatrixXd &M);

  static std::vector<Eigen::MatrixXd>
  scaleConstraints(const std::vector<std::vector<std::vector<double>>> &A,
                   const Eigen::MatrixXd &D);

  static bool isSPD(const Eigen::MatrixXd &M);

  static double primalObjective(const Eigen::MatrixXd &X,
                                const std::vector<std::vector<double>> &C);
  static double dualObjective(const std::vector<double> &y,
                              const std::vector<double> &b);

  double computeStepSizePD(const Iterate &it, const Direction &dir) const;
  void applyStep(Iterate &it, const Direction &dir, double alpha) const;

  static double trace(const Eigen::MatrixXd &M);
  static double dualityGap(const Eigen::MatrixXd &X,
                           const Eigen::MatrixXd &S);
  static double computeMaxStepSPD(const Eigen::MatrixXd &X,
                                  const Eigen::MatrixXd &dX);
};

#endif // IPM_ENGINE_H


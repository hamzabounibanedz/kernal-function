#include "kernel_base.h"
#include <cmath>
#include <stdexcept>
#include <tuple>

void KernelBase::validateInput(double t,
                               const std::string &function_name) const {
  if (!isValidInput(t)) {
    throw std::invalid_argument("Invalid input t=" + std::to_string(t) +
                                " in " + function_name +
                                ". Value must be positive.");
  }
}

std::vector<std::pair<double, double>>
KernelBase::evaluateRange(double t_min, double t_max, int num_points) const {
  if (t_min <= 0 || t_max <= 0 || t_min >= t_max || num_points <= 0) {
    throw std::invalid_argument(
        "Invalid range parameters for kernel evaluation");
  }

  std::vector<std::pair<double, double>> result;
  result.reserve(num_points);

  double step = (t_max - t_min) / (num_points - 1);

  for (int i = 0; i < num_points; ++i) {
    double t = t_min + i * step;
    result.emplace_back(t, psi(t));
  }

  return result;
}

std::tuple<double, double, double, double>
KernelBase::evaluateAll(double t) const {
  validateInput(t, "evaluateAll");
  return std::make_tuple(psi(t), psi_prime(t), psi_double_prime(t),
                         psi_triple_prime(t));
}
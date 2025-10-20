#include "kernel_manager.h"

#include "../kernels/touil_hyperbolic_kernel.h"
#include "../kernels/derbal20_param_log_kernel.h"
#include "../kernels/ijnna23_fractional_kernel.h"
#include "../kernels/wu_hyperbolic_kernel.h"
#include "../kernels/bachir_convex_combo_kernel.h"
#include "../kernels/parameterized_log_kernel.h"

KernelManager::KernelManager() {
  m_kernels.push_back(KernelInfo{
      "Touil & Chikouche (2022)",
      "Novel Hyperbolic Kernel ψ(t) = t² + ½sinh⁻²(coth t) - coth 1",
      {},
      []() { return std::make_shared<TouilHyperbolicKernel>(); }});

  m_kernels.push_back(KernelInfo{
      "Derbal & Kebbiche (2020)",
      "Parameterized Log Kernel ψ(t) = t² + log t + 2(e-1)/(t^q+1)",
      {KernelParameter{"q", 1.0, 10.0, 2.0}},
      []() { return std::make_shared<Derbal20ParamLogKernel>(2.0); }});

  m_kernels.push_back(KernelInfo{
      "IJNAA (2023)",
      "Fractional Barrier Kernel ψ(t) = t² + log t + 2/(t^q+1)",
      {KernelParameter{"q", 1.0, 10.0, 2.0}},
      []() { return std::make_shared<Ijna23FractionalKernel>(2.0); }});

  m_kernels.push_back(KernelInfo{
      "Wu & Zhang (2025)",
      "Hyperbolic Parameterized Kernel ψ(t) = t² + a·p·coth(pt)·sinh⁻²(p·coth(1))",
      {KernelParameter{"p", 0.5, 5.0, 1.0}},
      []() { return std::make_shared<WuHyperbolicKernel>(1.0); }});

  m_kernels.push_back(KernelInfo{
      "Bachir (2025)",
      "Convex Combination Kernel ψ(t) = t² + m·log t + (1-m)·exp(1/t-1)",
      {KernelParameter{"m", 0.01, 0.99, 0.5}},
      []() { return std::make_shared<BachirConvexComboKernel>(0.5); }});

  // Legacy Parameterized Log (default q=2.0)
  m_kernels.push_back(KernelInfo{
      "Parameterized Log (Legacy)", "ψ(t) with q parameter (legacy entry)",
      {KernelParameter{"q", 1.0, 10.0, 2.0}},
      []() { return std::make_shared<ParameterizedLogKernel>(2.0); }});
}

std::shared_ptr<KernelBase> KernelManager::createKernel(
    const KernelInfo &info, const std::vector<double> &paramValues) const {
  auto k = info.factory();
  if (!paramValues.empty()) k->setParameters(paramValues);
  return k;
}







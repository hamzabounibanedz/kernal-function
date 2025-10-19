#ifndef KERNEL_MANAGER_H
#define KERNEL_MANAGER_H

#include <functional>
#include <memory>
#include <vector>
#include <QString>

#include "../kernels/kernel_base.h"

struct KernelParameter {
  QString name;
  double minVal;
  double maxVal;
  double defaultVal;
};

struct KernelInfo {
  QString name;
  QString description;
  std::vector<KernelParameter> parameters;
  std::function<std::shared_ptr<KernelBase>()> factory;
};

class KernelManager {
public:
  KernelManager();

  const std::vector<KernelInfo> &getKernels() const { return m_kernels; }

  std::shared_ptr<KernelBase> createKernel(const KernelInfo &info,
                                           const std::vector<double> &paramValues) const;

private:
  std::vector<KernelInfo> m_kernels;
};

#endif // KERNEL_MANAGER_H



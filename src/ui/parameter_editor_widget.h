#ifndef PARAMETER_EDITOR_WIDGET_H
#define PARAMETER_EDITOR_WIDGET_H

#include <QDoubleSpinBox>
#include <QLabel>
#include <QSlider>
#include <QVBoxLayout>
#include <QWidget>
#include <vector>

#include "kernel_manager.h"

class ParameterEditorWidget : public QWidget {
  Q_OBJECT
public:
  explicit ParameterEditorWidget(QWidget *parent = nullptr);

  void displayParameters(const KernelInfo &info);
  void clearParameters();
  std::vector<double> currentValues() const;
  void setValues(const std::vector<double> &values);

signals:
  void valuesChanged(const std::vector<double> &values);

private:
  QVBoxLayout *m_layout;
  std::vector<QDoubleSpinBox *> m_spins;
  std::vector<QSlider *> m_sliders;
};

#endif // PARAMETER_EDITOR_WIDGET_H




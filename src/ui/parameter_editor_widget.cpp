#include "parameter_editor_widget.h"

#include <QHBoxLayout>

ParameterEditorWidget::ParameterEditorWidget(QWidget *parent) : QWidget(parent) {
  m_layout = new QVBoxLayout(this);
  m_layout->setSpacing(8);
  m_layout->setContentsMargins(0, 0, 0, 0);
}

void ParameterEditorWidget::clearParameters() {
  QLayoutItem *child;
  while ((child = m_layout->takeAt(0)) != nullptr) {
    if (child->widget()) child->widget()->deleteLater();
    delete child;
  }
  m_spins.clear();
  m_sliders.clear();
}

void ParameterEditorWidget::displayParameters(const KernelInfo &info) {
  clearParameters();
  for (const auto &p : info.parameters) {
    QWidget *row = new QWidget(this);
    QHBoxLayout *hl = new QHBoxLayout(row);
    hl->setSpacing(8);
    hl->setContentsMargins(0, 0, 0, 0);

    QLabel *label = new QLabel(p.name, row);
    QDoubleSpinBox *spin = new QDoubleSpinBox(row);
    spin->setRange(p.minVal, p.maxVal);
    spin->setValue(p.defaultVal);
    spin->setDecimals(3);
    spin->setSingleStep((p.maxVal - p.minVal) / 100.0);

    QSlider *slider = new QSlider(Qt::Horizontal, row);
    slider->setRange(static_cast<int>(p.minVal * 100), static_cast<int>(p.maxVal * 100));
    slider->setValue(static_cast<int>(p.defaultVal * 100));

    QObject::connect(slider, &QSlider::valueChanged, [this, spin](int v) {
      spin->setValue(v / 100.0);
      std::vector<double> vals; vals.reserve(m_spins.size());
      for (auto *s : m_spins) vals.push_back(s->value());
      emit valuesChanged(vals);
    });
    QObject::connect(spin, qOverload<double>(&QDoubleSpinBox::valueChanged), [this, slider](double v) {
      slider->setValue(static_cast<int>(v * 100));
      std::vector<double> vals; vals.reserve(m_spins.size());
      for (auto *s : m_spins) vals.push_back(s->value());
      emit valuesChanged(vals);
    });

    hl->addWidget(label);
    hl->addWidget(spin, 1);
    hl->addWidget(slider, 2);
    row->setLayout(hl);
    m_layout->addWidget(row);

    m_spins.push_back(spin);
    m_sliders.push_back(slider);
  }
  m_layout->addStretch();
}

std::vector<double> ParameterEditorWidget::currentValues() const {
  std::vector<double> vals;
  vals.reserve(m_spins.size());
  for (auto *s : m_spins) vals.push_back(s->value());
  return vals;
}

void ParameterEditorWidget::setValues(const std::vector<double> &values) {
  const int n = std::min<int>(static_cast<int>(values.size()), static_cast<int>(m_spins.size()));
  for (int i = 0; i < n; ++i) {
    m_spins[i]->setValue(values[i]);
    m_sliders[i]->setValue(static_cast<int>(values[i] * 100));
  }
}




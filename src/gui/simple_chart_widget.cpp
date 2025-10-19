#include "simple_chart_widget.h"
#include <QtCore/Qt>
#include <QtGui/QBrush>
#include <QtGui/QFont>
#include <QtGui/QPainter>
#include <QtGui/QPen>
#include <QtWidgets/QLabel>
#include <QtWidgets/QVBoxLayout>

SimpleChartWidget::SimpleChartWidget(QWidget *parent) : QWidget(parent) {
  setMinimumSize(400, 300);
  setStyleSheet(
      "background: qlineargradient(x1:0, y1:0, x2:0, y2:1,"
      "    stop:0 #1a1a2e, stop:1 #16213e); border: 2px solid #4a90e2; "
      "border-radius: 10px;");
}

void SimpleChartWidget::setData(const QVector<double> &values,
                                const QVector<QString> &labels) {
  m_values = values;
  m_labels = labels;
  update();
}

void SimpleChartWidget::setTitle(const QString &title) {
  m_title = title;
  update();
}

void SimpleChartWidget::setYAxisLabel(const QString &label) {
  m_yAxisLabel = label;
  update();
}

void SimpleChartWidget::paintEvent(QPaintEvent *event) {
  Q_UNUSED(event)

  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);

  // Set up the coordinate system
  int margin = 60;
  int chartWidth = width() - 2 * margin;
  int chartHeight = height() - 2 * margin;

  // Draw title
  if (!m_title.isEmpty()) {
    painter.setPen(QPen(QColor("#00d4ff")));
    painter.setFont(QFont("Segoe UI", 14, QFont::Bold));
    painter.drawText(rect(), Qt::AlignTop | Qt::AlignHCenter, m_title);
  }

  // Draw chart area
  QRect chartRect(margin, margin, chartWidth, chartHeight);

  // Draw axis
  drawAxis(painter);

  // Draw bars
  drawBarChart(painter);

  // Draw labels
  drawLabels(painter);
}

void SimpleChartWidget::drawAxis(QPainter &painter) {
  int margin = 60;
  int chartWidth = width() - 2 * margin;
  int chartHeight = height() - 2 * margin;

  painter.setPen(QPen(QColor("#4a90e2"), 2));

  // X-axis
  painter.drawLine(margin, height() - margin, width() - margin,
                   height() - margin);

  // Y-axis
  painter.drawLine(margin, margin, margin, height() - margin);

  // Y-axis label
  if (!m_yAxisLabel.isEmpty()) {
    painter.setPen(QPen(QColor("#00d4ff")));
    painter.setFont(QFont("Segoe UI", 10, QFont::Bold));
    painter.save();
    painter.translate(20, height() / 2);
    painter.rotate(-90);
    painter.drawText(QRect(-50, 0, 100, 20), Qt::AlignCenter, m_yAxisLabel);
    painter.restore();
  }
}

void SimpleChartWidget::drawBarChart(QPainter &painter) {
  if (m_values.isEmpty())
    return;

  int margin = 60;
  int chartWidth = width() - 2 * margin;
  int chartHeight = height() - 2 * margin;

  double maxValue = *std::max_element(m_values.begin(), m_values.end());
  if (maxValue <= 0)
    maxValue = 1.0;

  int barCount = m_values.size();
  int barWidth = chartWidth / (barCount + 1);
  int barSpacing = barWidth / 4;

  QVector<QColor> colors = {QColor("#00d4ff"), QColor("#4a90e2"),
                            QColor("#ff6b6b"), QColor("#feca57"),
                            QColor("#48dbfb"), QColor("#ff9ff3")};

  for (int i = 0; i < barCount; ++i) {
    double value = m_values[i];
    int barHeight = static_cast<int>((value / maxValue) * chartHeight);

    int x = margin + (i + 1) * barWidth - barWidth / 2;
    int y = height() - margin - barHeight;

    QRect barRect(x - barWidth / 2 + barSpacing, y, barWidth - 2 * barSpacing,
                  barHeight);

    // Draw bar with gradient
    QLinearGradient gradient(barRect.topLeft(), barRect.bottomLeft());
    gradient.setColorAt(0, colors[i % colors.size()]);
    gradient.setColorAt(1, colors[i % colors.size()].darker(120));

    painter.setBrush(QBrush(gradient));
    painter.setPen(QPen(colors[i % colors.size()].darker(150), 1));
    painter.drawRect(barRect);

    // Draw value on top of bar
    painter.setPen(QPen(QColor("#ffffff")));
    painter.setFont(QFont("Segoe UI", 8, QFont::Bold));
    painter.drawText(barRect, Qt::AlignCenter, QString::number(value, 'f', 2));
  }
}

void SimpleChartWidget::drawLabels(QPainter &painter) {
  if (m_labels.isEmpty())
    return;

  int margin = 60;
  int chartWidth = width() - 2 * margin;
  int chartHeight = height() - 2 * margin;

  int barCount = m_labels.size();
  int barWidth = chartWidth / (barCount + 1);

  painter.setPen(QPen(QColor("#00d4ff")));
  painter.setFont(QFont("Segoe UI", 9, QFont::Bold));

  for (int i = 0; i < barCount; ++i) {
    int x = margin + (i + 1) * barWidth - barWidth / 2;
    int y = height() - margin + 10;

    QString label = m_labels[i];
    if (label.length() > 15) {
      label = label.left(12) + "...";
    }

    painter.drawText(QRect(x - barWidth / 2, y, barWidth, 20), Qt::AlignCenter,
                     label);
  }
}
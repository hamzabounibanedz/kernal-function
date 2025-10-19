#ifndef SIMPLE_CHART_WIDGET_H
#define SIMPLE_CHART_WIDGET_H

#include <QPainter>
#include <QString>
#include <QVector>
#include <QWidget>

class SimpleChartWidget : public QWidget {
  Q_OBJECT

public:
  explicit SimpleChartWidget(QWidget *parent = nullptr);

  void setData(const QVector<double> &values, const QVector<QString> &labels);
  void setTitle(const QString &title);
  void setYAxisLabel(const QString &label);

protected:
  void paintEvent(QPaintEvent *event) override;

private:
  QVector<double> m_values;
  QVector<QString> m_labels;
  QString m_title;
  QString m_yAxisLabel;

  void drawBarChart(QPainter &painter);
  void drawAxis(QPainter &painter);
  void drawLabels(QPainter &painter);
};

#endif // SIMPLE_CHART_WIDGET_H
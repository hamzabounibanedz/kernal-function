#ifndef CPU_TIME_CHART_WIDGET_H
#define CPU_TIME_CHART_WIDGET_H

#include <QWidget>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QLogValueAxis>
#include <QtCharts/QCategoryAxis>
#include <QToolButton>
#include <QCheckBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMap>

class CpuTimeChartWidget : public QWidget {
  Q_OBJECT
public:
  explicit CpuTimeChartWidget(QWidget* parent=nullptr);

  // Append or update a kernel's CPU time (seconds)
  void upsertSample(const QString& kernelName, double seconds);
  void clearAll();

private slots:
  void onBarMode();
  void onLineMode();
  void onLogScaleToggled(bool);
  void onExportPng();

private:
  void rebuildAxes();
  void refresh();

  QChart* m_chart;
  QChartView* m_chartView;
  QBarSeries* m_barSeries;
  QBarSet* m_barSet;
  QLineSeries* m_lineSeries;
  QValueAxis* m_yLinear;
  QLogValueAxis* m_yLog;
  QCategoryAxis* m_xCats;
  QToolButton* m_btnBar;
  QToolButton* m_btnLine;
  QCheckBox* m_chkLog;
  QToolButton* m_btnExport;
  QMap<QString,double> m_nameToSeconds;
};

#endif // CPU_TIME_CHART_WIDGET_H




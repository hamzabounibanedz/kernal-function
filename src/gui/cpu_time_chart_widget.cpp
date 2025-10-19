#include "cpu_time_chart_widget.h"
#include <QFileDialog>
#include <QPixmap>

CpuTimeChartWidget::CpuTimeChartWidget(QWidget* parent)
  : QWidget(parent),
    m_chart(new QChart()),
    m_chartView(new QChartView(m_chart)),
    m_barSeries(new QBarSeries()),
    m_barSet(new QBarSet("CPU Time (s)")),
    m_lineSeries(new QLineSeries()),
    m_yLinear(new QValueAxis()),
    m_yLog(new QLogValueAxis()),
    m_xCats(new QCategoryAxis()),
    m_btnBar(new QToolButton()),
    m_btnLine(new QToolButton()),
    m_chkLog(new QCheckBox("Log scale")),
    m_btnExport(new QToolButton()) {

  // Header controls
  m_btnBar->setText("Bar");
  m_btnLine->setText("Line");
  m_btnBar->setCheckable(true);
  m_btnLine->setCheckable(true);
  m_btnBar->setChecked(false);
  m_btnExport->setText("Export");

  connect(m_btnBar, &QToolButton::clicked, this, &CpuTimeChartWidget::onBarMode);
  connect(m_btnLine, &QToolButton::clicked, this, &CpuTimeChartWidget::onLineMode);
  connect(m_chkLog, &QCheckBox::toggled, this, &CpuTimeChartWidget::onLogScaleToggled);
  connect(m_btnExport, &QToolButton::clicked, this, &CpuTimeChartWidget::onExportPng);

  // Chart setup
  m_barSeries->append(m_barSet);
  m_chart->addSeries(m_barSeries);
  m_chart->addSeries(m_lineSeries);
  m_chart->legend()->hide();
  m_chartView->setRenderHint(QPainter::Antialiasing);

  m_yLinear->setTitleText("Seconds");
  m_yLog->setTitleText("Seconds (log)");
  m_yLog->setBase(10.0);
  m_chart->addAxis(m_yLinear, Qt::AlignLeft);
  m_chart->addAxis(m_xCats, Qt::AlignBottom);

  m_barSeries->attachAxis(m_yLinear);
  m_barSeries->attachAxis(m_xCats);
  m_lineSeries->attachAxis(m_yLinear);
  m_lineSeries->attachAxis(m_xCats);

  // Layout
  auto header = new QWidget(this);
  auto headerLay = new QHBoxLayout(header);
  headerLay->setContentsMargins(0,0,0,0);
  headerLay->addWidget(m_btnBar);
  headerLay->addWidget(m_btnLine);
  headerLay->addSpacing(12);
  headerLay->addWidget(m_chkLog);
  headerLay->addStretch();
  headerLay->addWidget(m_btnExport);

  auto lay = new QVBoxLayout(this);
  lay->setContentsMargins(8,8,8,8);
  lay->setSpacing(8);
  lay->addWidget(header);
  lay->addWidget(m_chartView, 1);

  rebuildAxes();
  // Default view: Line + Log scale for better readability of large, non-harmonic values
  onLineMode();
  m_chkLog->setChecked(true);
  onLogScaleToggled(true);
}

void CpuTimeChartWidget::rebuildAxes() {
  // Rebuild categories in insertion order
  // QtCharts QCategoryAxis has no clear(); rebuild via replacement
  auto newCats = new QCategoryAxis();
  m_chart->removeAxis(m_xCats);
  m_xCats->deleteLater();
  m_xCats = newCats;
  m_chart->addAxis(m_xCats, Qt::AlignBottom);
  m_barSeries->attachAxis(m_xCats);
  m_lineSeries->attachAxis(m_xCats);

  m_barSet->remove(0, m_barSet->count());
  m_lineSeries->clear();

  int idx = 0;
  double maxY = 1.0;
  for (auto it = m_nameToSeconds.begin(); it != m_nameToSeconds.end(); ++it) {
    const QString name = it.key();
    const double s = it.value();
    m_xCats->append(name, idx);
    *m_barSet << s;
    m_lineSeries->append(idx, s);
    maxY = std::max(maxY, s);
    ++idx;
  }
  m_yLinear->setRange(0.0, maxY * 1.1);
  m_yLog->setRange(0.1, std::max(1.0, maxY * 1.1));
}

void CpuTimeChartWidget::refresh() {
  rebuildAxes();
}

void CpuTimeChartWidget::upsertSample(const QString& kernelName, double seconds) {
  m_nameToSeconds[kernelName] = seconds;
  refresh();
}

void CpuTimeChartWidget::clearAll() {
  m_nameToSeconds.clear();
  refresh();
}

void CpuTimeChartWidget::onBarMode() {
  m_btnBar->setChecked(true);
  m_btnLine->setChecked(false);
  m_barSeries->setVisible(true);
  m_lineSeries->setVisible(false);
}

void CpuTimeChartWidget::onLineMode() {
  m_btnBar->setChecked(false);
  m_btnLine->setChecked(true);
  m_barSeries->setVisible(false);
  m_lineSeries->setVisible(true);
}

void CpuTimeChartWidget::onLogScaleToggled(bool on) {
  // Switch y-axis
  if (on) {
    m_chart->removeAxis(m_yLinear);
    m_chart->addAxis(m_yLog, Qt::AlignLeft);
    m_barSeries->attachAxis(m_yLog);
    m_lineSeries->attachAxis(m_yLog);
  } else {
    m_chart->removeAxis(m_yLog);
    m_chart->addAxis(m_yLinear, Qt::AlignLeft);
    m_barSeries->attachAxis(m_yLinear);
    m_lineSeries->attachAxis(m_yLinear);
  }
}

void CpuTimeChartWidget::onExportPng() {
  QString file = QFileDialog::getSaveFileName(this, "Export CPU Time Chart", "cpu_time.png", "PNG Files (*.png)");
  if (file.isEmpty()) return;
  QPixmap pm = m_chartView->grab();
  pm.save(file, "PNG");
}




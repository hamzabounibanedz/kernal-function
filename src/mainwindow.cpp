#include "mainwindow.h"
#include "gui/algorithm_dialog.h"
#include "gui/matrix_input_widget.h"
#include "kernels/bachir_kernel.h"
#include "kernels/bachir_convex_combo_kernel.h"
#include "kernels/exponential_parametric_kernel.h"
#include "kernels/ijnna23_fractional_kernel.h"
#include "kernels/derbal20_param_log_kernel.h"
#include "kernels/log_exponential_kernel.h"
#include "kernels/parameterized_log_kernel.h"
#include "kernels/parametric_family_kernel.h"
#include "kernels/touil_hyperbolic_kernel.h"
#include "kernels/wu_hyperbolic_kernel.h"

#include "utils/test_data.h"

#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QStyle>
#include <QHeaderView>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QFileDialog>
#include <cmath>
#include <map>

// Persist per-kernel edited parameter values across selections in this session
static std::map<QString, std::vector<double>> g_kernelEditedValues;
#include <QProgressBar>
#include <QPushButton>
#include <QScreen>
#include <QScrollArea>
#include <QSplitter>
#include <QStatusBar>
#include <QTabWidget>
#include <QTableWidget>
#include <QTextEdit>
#include <QTimer>
#include <QVBoxLayout>
#include <memory>
#include <random>
#include <QFile>
#include <QTextStream>

std::shared_ptr<KernelBase> MainWindow::createSelectedKernel() const {
  // Prefer the new multi-selection list with per-kernel parameters
  if (m_kernelList) {
    auto items = m_kernelList->selectedItems();
    if (items.size() == 1) {
      QString name = items.first()->text();
      for (const auto &info : m_kernelManager.getKernels()) {
        if (info.name == name) {
          std::vector<double> values;
          auto kv = g_kernelEditedValues.find(name);
          if (kv != g_kernelEditedValues.end()) values = kv->second;
          else if (m_parameterEditor) values = m_parameterEditor->currentValues();
          if (values.empty()) { for (const auto &p : info.parameters) values.push_back(p.defaultVal); }
          return m_kernelManager.createKernel(info, values);
        }
      }
    }
  }

  // Fallback: first registered kernel with default parameters
  const auto &all = m_kernelManager.getKernels();
  if (!all.empty()) {
    const auto &info = all.front();
    std::vector<double> values; for (const auto &p : info.parameters) values.push_back(p.defaultVal);
    return m_kernelManager.createKernel(info, values);
  }

  // Ultimate fallback
  return std::make_shared<ParameterizedLogKernel>(2.0);
}

void MainWindow::onKernelChanged(int index) {
  // Update parameter label/range per kernel
  if (index == 0) { // Touil22 no params
    m_kernelParamLabel->setText("(no parameter)");
    m_kernelParamSpin->setEnabled(false);
  } else if (index == 1) { // Derbal20 q
    m_kernelParamLabel->setText("q (1..10)");
    m_kernelParamSpin->setRange(1.0, 10.0);
    m_kernelParamSpin->setValue(2.0);
    m_kernelParamSpin->setEnabled(true);
  } else if (index == 2) { // IJNAA23 q
    m_kernelParamLabel->setText("q (1..10)");
    m_kernelParamSpin->setRange(1.0, 10.0);
    m_kernelParamSpin->setValue(2.0);
    m_kernelParamSpin->setEnabled(true);
  } else if (index == 3) { // Wu25 p
    m_kernelParamLabel->setText("p (0.5..5)");
    m_kernelParamSpin->setRange(0.5, 5.0);
    m_kernelParamSpin->setValue(1.0);
    m_kernelParamSpin->setEnabled(true);
  } else if (index == 4) { // Bachir m
    m_kernelParamLabel->setText("m (0.01..0.99)");
    m_kernelParamSpin->setRange(0.01, 0.99);
    m_kernelParamSpin->setValue(0.5);
    m_kernelParamSpin->setEnabled(true);
  }
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), m_isRunning(false) {
  setupUI();
  initializeData();
  setupTableHeaders();

  // Set window properties
  setWindowTitle("SDP Kernel Benchmarking Suite v1.0");
  setWindowIcon(QIcon());
  resize(1400, 900);

  // Center window
  move(QApplication::primaryScreen()->geometry().center() -
       frameGeometry().center());
}

MainWindow::~MainWindow() = default;


void MainWindow::setupUI() {
  // Create central widget and main layout
  m_centralWidget = new QWidget(this);
  setCentralWidget(m_centralWidget);

  // Use app-level QSS; avoid local overrides that conflict with theme

  // Main layout with pixel-perfect spacing
  m_mainLayout = new QHBoxLayout(m_centralWidget);
  m_mainLayout->setSpacing(0);
  m_mainLayout->setContentsMargins(0, 0, 0, 0);

  // Create main splitter for responsive layout
  m_mainSplitter = new QSplitter(Qt::Horizontal, this);
  m_mainSplitter->setObjectName("MainSplitter");
  m_mainSplitter->setHandleWidth(4); // Splitter Handle: 4px wide
  m_mainSplitter->setChildrenCollapsible(false);

  m_mainLayout->addWidget(m_mainSplitter);

  setupMenuBar();
  setupInputPanel();
  setupResultsPanel();

  // Create light theme status bar
  statusBar()->showMessage("Ready - No kernels selected");
}

void MainWindow::setupMenuBar() {
  QMenu *fileMenu = menuBar()->addMenu("&File");
  fileMenu->addAction("&New");
  fileMenu->addAction("&Open");
  fileMenu->addAction("&Save");
  fileMenu->addAction("&Export", this, &MainWindow::exportResults);
  fileMenu->addSeparator();
  fileMenu->addAction("E&xit", this, &QWidget::close);

  QMenu *viewMenu = menuBar()->addMenu("&View");
  viewMenu->addAction("&Light Theme");
  viewMenu->addAction("&Dark Theme");
  viewMenu->addAction("&Reset Layout");

  QMenu *helpMenu = menuBar()->addMenu("&Help");
  helpMenu->addAction("&Documentation");
  helpMenu->addAction("&About", this, &MainWindow::showAlgorithmInfo);
}

void MainWindow::setupInputPanel() {
  m_leftTabs = new QTabWidget();

  // Benchmark Setup tab
  QWidget *setupTab = new QWidget();
  QVBoxLayout *setupLay = new QVBoxLayout(setupTab);
  QWidget *listsRow = new QWidget(); QHBoxLayout *h = new QHBoxLayout(listsRow);
  m_kernelList = new QListWidget(); m_kernelList->setSelectionMode(QAbstractItemView::ExtendedSelection);
  m_testList = new QListWidget(); m_testList->setSelectionMode(QAbstractItemView::ExtendedSelection);
  h->addWidget(m_kernelList, 1); h->addWidget(m_testList, 1);
  m_parameterEditor = new ParameterEditorWidget();
  setupLay->addWidget(listsRow); setupLay->addWidget(m_parameterEditor);
  m_leftTabs->addTab(setupTab, "Benchmark Setup");

  // Problem Editor tab
  m_matrixInputWidget = new MatrixInputWidget(this);
  m_leftTabs->addTab(m_matrixInputWidget, "Problem Editor");

  // Run Settings tab
  QWidget *runTab = new QWidget(); QVBoxLayout *runLay = new QVBoxLayout(runTab);
  QGroupBox *paramsGroup = new QGroupBox("Run Settings"); QFormLayout *form = new QFormLayout(paramsGroup);
  m_toleranceSpinBox = new QDoubleSpinBox(); m_toleranceSpinBox->setRange(1e-12, 1e-1); m_toleranceSpinBox->setDecimals(8); m_toleranceSpinBox->setValue(1e-8);
  m_thetaSpinBox = new QDoubleSpinBox(); m_thetaSpinBox->setRange(0.01, 0.99); m_thetaSpinBox->setDecimals(2); m_thetaSpinBox->setValue(0.1);
  m_tauSpinBox = new QDoubleSpinBox(); m_tauSpinBox->setRange(0.1, 10.0); m_tauSpinBox->setDecimals(2); m_tauSpinBox->setValue(1.0);
  form->addRow("Tolerance (ε)", m_toleranceSpinBox);
  form->addRow("Theta (θ)", m_thetaSpinBox);
  form->addRow("Tau (τ)", m_tauSpinBox);
  runLay->addWidget(paramsGroup);
  // Controls
  QWidget *ctrl = new QWidget(); QHBoxLayout *ctrlLay = new QHBoxLayout(ctrl);
  m_runAllButton = new QPushButton("Run All");
  m_runSingleButton = new QPushButton("Run One");
  m_cancelButton = new QPushButton("Cancel"); m_cancelButton->setEnabled(false);
  // Style properties for QSS
  m_runAllButton->setProperty("type", "primary");
  m_cancelButton->setProperty("type", "danger");
  ctrlLay->addWidget(m_runAllButton); ctrlLay->addWidget(m_runSingleButton); ctrlLay->addWidget(m_cancelButton); ctrlLay->addStretch();
  connect(m_runAllButton, &QPushButton::clicked, this, &MainWindow::runAllAlgorithms);
  connect(m_runSingleButton, &QPushButton::clicked, this, &MainWindow::runSingleAlgorithm);
  connect(m_cancelButton, &QPushButton::clicked, [this](){ m_cancel.store(true); m_statusLabel->setText("Cancelling..."); });
  runLay->addWidget(ctrl);
  m_progressBar = new QProgressBar(); m_progressBar->setVisible(false); runLay->addWidget(m_progressBar);
  m_statusLabel = new QLabel("Ready - No kernels selected");
  m_statusLabel->setObjectName("StatusBanner");
  m_statusLabel->setProperty("status", "idle");
  m_statusLabel->style()->unpolish(m_statusLabel); m_statusLabel->style()->polish(m_statusLabel);
  runLay->addWidget(m_statusLabel);
  m_progressText = new QLabel(""); runLay->addWidget(m_progressText);
  m_leftTabs->addTab(runTab, "Run Settings");

  // Populate kernels and tests
  for (const auto &k : m_kernelManager.getKernels()) m_kernelList->addItem(k.name);
  // Tests will be populated in initializeData()

  // Kernel parameter editor wiring
  connect(m_parameterEditor, &ParameterEditorWidget::valuesChanged, this, [this](const std::vector<double> &vals){
    auto items = m_kernelList->selectedItems();
    if (items.size() == 1) {
      g_kernelEditedValues[items.first()->text()] = vals;
    }
  });
  connect(m_kernelList, &QListWidget::itemSelectionChanged, this, [this]() {
    auto items = m_kernelList->selectedItems();
    if (items.size() == 1) {
      QString name = items.first()->text();
      for (const auto &info : m_kernelManager.getKernels()) {
        if (info.name == name) {
          m_parameterEditor->displayParameters(info);
          auto it = g_kernelEditedValues.find(name);
          if (it != g_kernelEditedValues.end()) m_parameterEditor->setValues(it->second);
          if (m_statusLabel) {
            int k = m_kernelList->selectedItems().size();
            int t = m_testList ? m_testList->selectedItems().size() : 0;
            m_statusLabel->setText(QString("Ready - %1 kernel%2 selected, %3 test%4 selected")
                                    .arg(k).arg(k==1?"":"s").arg(t).arg(t==1?"":"s"));
          }
          return;
        }
      }
    }
    m_parameterEditor->clearParameters();
    if (m_statusLabel) {
      int k = m_kernelList->selectedItems().size();
      int t = m_testList ? m_testList->selectedItems().size() : 0;
      m_statusLabel->setText(QString("Ready - %1 kernel%2 selected, %3 test%4 selected")
                              .arg(k).arg(k==1?"":"s").arg(t).arg(t==1?"":"s"));
    }
  });

  connect(m_testList, &QListWidget::itemSelectionChanged, this, [this]() {
    if (m_statusLabel) {
      int k = m_kernelList ? m_kernelList->selectedItems().size() : 0;
      int t = m_testList ? m_testList->selectedItems().size() : 0;
      m_statusLabel->setText(QString("Ready - %1 kernel%2 selected, %3 test%4 selected")
                              .arg(k).arg(k==1?"":"s").arg(t).arg(t==1?"":"s"));
    }
  });

  m_mainSplitter->addWidget(m_leftTabs);
}

void MainWindow::setupResultsPanel() {
  m_rightPanel = new QWidget();
  m_rightPanel->setObjectName("RightPanel");
  m_rightLayout = new QVBoxLayout(m_rightPanel);
  m_rightLayout->setSpacing(12);
  m_rightLayout->setContentsMargins(16, 16, 16, 16);

  m_resultsTabWidget = new QTabWidget();
  m_tableTab = new QWidget();
  m_tableLayout = new QVBoxLayout(m_tableTab);
  m_resultsTable = new QTableWidget();
  m_resultsTable->setAlternatingRowColors(true);
  m_resultsTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  m_resultsTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_resultsTable->verticalHeader()->setVisible(false);
  m_resultsTable->verticalHeader()->setDefaultSectionSize(32); // taller rows for readability
  m_resultsTable->setShowGrid(false);
  m_tableLayout->addWidget(m_resultsTable);

  QHBoxLayout *exportRow = new QHBoxLayout();
  exportRow->setSpacing(8);
  exportRow->setContentsMargins(0, 8, 0, 0);
  m_exportCsvButton = new QPushButton("📊 Export CSV");
  m_exportPngButton = new QPushButton("📈 Export PNG");
  m_exportPdfButton = new QPushButton("📄 Export PDF");
  connect(m_exportCsvButton, &QPushButton::clicked, this, &MainWindow::exportResults);
  connect(m_exportPngButton, &QPushButton::clicked, [this](){ QMessageBox::information(this, "Export PNG", "PNG export will render the table as an image in a future update."); });
  connect(m_exportPdfButton, &QPushButton::clicked, [this](){ QMessageBox::information(this, "Export PDF", "PDF export will generate a report in a future update."); });
  exportRow->addWidget(m_exportCsvButton);
  exportRow->addWidget(m_exportPngButton);
  exportRow->addWidget(m_exportPdfButton);
  exportRow->addStretch();
  m_tableLayout->addLayout(exportRow);

  m_resultsTabWidget->addTab(m_tableTab, "Results Table");

  // CPU Time chart tab (only if charts enabled at build time)
#ifdef CPU_CHARTS_ENABLED
  m_chartTab = new QWidget();
  m_chartLayout = new QVBoxLayout(m_chartTab);
  m_cpuChart = new CpuTimeChartWidget(m_chartTab);
  m_chartLayout->addWidget(static_cast<CpuTimeChartWidget*>(m_cpuChart));
  m_resultsTabWidget->addTab(m_chartTab, "CPU Time");
#endif

  m_logTab = new QWidget();
  m_logLayout = new QVBoxLayout(m_logTab);
  m_logTextEdit = new QTextEdit();
  m_logTextEdit->setObjectName("ExecutionLog");
  m_logLayout->addWidget(m_logTextEdit);
  m_resultsTabWidget->addTab(m_logTab, "Execution Log");

  m_rightLayout->addWidget(m_resultsTabWidget);
  m_mainSplitter->addWidget(m_rightPanel);
  m_mainSplitter->setSizes({420, 980});
}

void MainWindow::initializeData() {
  // Initialize legacy kernel list for dialogs (optional)
  m_kernels.clear();
  for (const auto &info : m_kernelManager.getKernels()) {
    auto k = info.factory();
    m_kernels.push_back(k);
  }

  // Load test cases and populate list
  m_testCases = TestDataProvider::getAllTestCases();
  if (m_testList) {
    m_testList->clear();
    for (const auto &tc : m_testCases) m_testList->addItem(QString::fromStdString(tc.name));
  }
}

void MainWindow::setupTableHeaders() {
  m_resultsTable->setColumnCount(7);
  QStringList headers = {"Algorithm", "Parameters", "Test Case", "Iterations", "CPU Time", "Duality Gap", "Status"};
  m_resultsTable->setHorizontalHeaderLabels(headers);
  m_resultsTable->horizontalHeader()->setStretchLastSection(true);
}

TestCase MainWindow::createTestCaseFromMatrixInput() {
  TestCase testCase;

  // Get problem size from MatrixInputWidget
  auto problemSize = m_matrixInputWidget->getProblemSize();
  testCase.problem_size = problemSize.first; // n
  int m = problemSize.second;                // m

  // Get matrix data from MatrixInputWidget
  testCase.A_matrices = m_matrixInputWidget->getAMatrices();
  testCase.b_vector = m_matrixInputWidget->getBVector();
  testCase.C_matrix = m_matrixInputWidget->getCMatrix();
  testCase.initial_X = m_matrixInputWidget->getInitialX();
  testCase.initial_y = m_matrixInputWidget->getInitialY();
  testCase.initial_S = m_matrixInputWidget->getInitialS();

  // Get parameters from spinboxes
  testCase.tolerance = m_toleranceSpinBox->value();
  testCase.theta = m_thetaSpinBox->value();
  testCase.tau = m_tauSpinBox->value();

  // Set a default name
  testCase.name = QString("Custom Problem (n=%1, m=%2)")
                      .arg(testCase.problem_size)
                      .arg(m)
                      .toStdString();
  testCase.expected_optimal_value = -1.0957; // Default value

  return testCase;
}

void MainWindow::runAllAlgorithms() {
  if (m_isRunning) return;

  auto kernels = getSelectedKernels();
  auto tests = getSelectedTests();
  if (kernels.empty() || tests.empty()) {
    QMessageBox::warning(this, "Run All", kernels.empty() ? "No kernel selected." : "No test selected.");
    m_statusLabel->setProperty("status", "error"); m_statusLabel->style()->unpolish(m_statusLabel); m_statusLabel->style()->polish(m_statusLabel);
    m_statusLabel->setText("Cannot run - please select at least one kernel and one test.");
    return;
  }
  int total = int(kernels.size() * tests.size());

  m_isRunning = true;
  m_runAllButton->setEnabled(false);
  m_runSingleButton->setEnabled(false);
  m_cancelButton->setEnabled(true);
  m_progressBar->setVisible(true);
  m_progressBar->setRange(0, std::max(1,total));
  m_progressBar->setValue(0);
  m_statusLabel->setText(QString("Ready - %1 kernels selected, %2 tests selected").arg(kernels.size()).arg(tests.size()));
  m_statusLabel->setProperty("status", "running"); m_statusLabel->style()->unpolish(m_statusLabel); m_statusLabel->style()->polish(m_statusLabel);
  // Immediate feedback in log
  m_logTextEdit->append(QString("Queued %1×%2 = %3 runs...").arg(kernels.size()).arg(tests.size()).arg(total));

  m_resultsTable->setRowCount(0);
  m_logTextEdit->clear();
  m_cancel.store(false);

  // Reset CPU chart at the start of a new batch
  #ifdef CPU_CHARTS_ENABLED
  if (m_cpuChart) m_cpuChart->clearAll();
  #endif

  m_future = QtConcurrent::run([this, kernels, tests, total]() {
    int completed = 0;
    for (const auto &k : kernels) {
      if (m_cancel.load()) break;
      for (const auto &tc : tests) {
        if (m_cancel.load()) break;
        QString runningMsg = QString("Running - %1 on %2...").arg(k.first).arg(QString::fromStdString(tc.name));
        QMetaObject::invokeMethod(this, [this, runningMsg, completed, total]() {
          m_statusLabel->setText(runningMsg);
          m_progressText->setText(formatProgressBar(completed, total));
          // Show first activity immediately
          if (m_logTextEdit && completed == 0) {
            m_logTextEdit->append("Started first run...");
            // Indeterminate progress until first completion
            m_progressBar->setRange(0, 0);
          }
        });

        // Use Run Settings from UI, not the test-case defaults
        IpmSettings settings; 
        settings.epsilon = m_toleranceSpinBox ? m_toleranceSpinBox->value() : tc.tolerance;
        settings.theta   = m_thetaSpinBox ? m_thetaSpinBox->value() : tc.theta;
        settings.tau     = m_tauSpinBox ? m_tauSpinBox->value() : tc.tau;
        settings.cancel  = &m_cancel;
        // Increase iteration budgets for convergence
        settings.maxOuter = 1000;
        settings.maxInner = 100;
        auto start = std::chrono::high_resolution_clock::now();
        AlgorithmResult result;
        try {
          // Guard against invalid starting data that can produce NaNs
          if (!std::isfinite(settings.epsilon) || settings.epsilon <= 0.0) throw std::runtime_error("Invalid epsilon");
          if (!std::isfinite(settings.theta) || settings.theta <= 0.0 || settings.theta >= 1.0) throw std::runtime_error("Invalid theta");
          if (!std::isfinite(settings.tau) || settings.tau <= 0.0) throw std::runtime_error("Invalid tau");
          // Validate initial matrices and vectors for finiteness
          auto finiteMatrix = [](const std::vector<std::vector<double>>& M){
            for (const auto &row : M) for (double v : row) if (!std::isfinite(v)) return false; return true; };
          auto finiteVector = [](const std::vector<double>& V){ for (double v : V) if (!std::isfinite(v)) return false; return true; };
          if (!finiteMatrix(tc.C_matrix) || !finiteVector(tc.b_vector)) throw std::runtime_error("Invalid C or b data");
          if (!tc.A_matrices.empty()) {
            for (const auto &A : tc.A_matrices) if (!finiteMatrix(A)) throw std::runtime_error("Invalid A matrix data");
          }
          if (!finiteMatrix(tc.initial_X) || !finiteMatrix(tc.initial_S) || !finiteVector(tc.initial_y)) throw std::runtime_error("Invalid initial X/S/y");

          // Log the settings used for this run
          // Log full run configuration including kernel parameters
          auto kernelParams = k.second->getParameters();
          QStringList kpStr; for (double v : kernelParams) kpStr << QString::number(v, 'g', 6);
          QMetaObject::invokeMethod(this, [this, settings, k, tc, kpStr]() {
            if (m_logTextEdit) m_logTextEdit->append(
              QString("Run config → Kernel: %1, Test: %2, ε=%3, θ=%4, τ=%5, params=[%6]")
                .arg(k.first)
                .arg(QString::fromStdString(tc.name))
                .arg(settings.epsilon, 0, 'g', 6)
                .arg(settings.theta, 0, 'g', 6)
                .arg(settings.tau, 0, 'g', 6)
                .arg(kpStr.join(", "))
            );
          });

          result = m_engine.solve(tc, *k.second, settings);
        } catch (const std::exception &e) {
          result.converged = false;
          result.convergence_info = std::string("Error: ") + e.what();
        } catch (...) {
          result.converged = false;
          result.convergence_info = std::string("Unknown error during solve");
        }
        auto end = std::chrono::high_resolution_clock::now();
        // Prefer precise CPU time if provided; else fall back to measured wall time
        double cpuSec = (result.cpu_time_ms > 0.0)
            ? (result.cpu_time_ms / 1000.0)
            : ((result.execution_time_ms > 0.0)
               ? (result.execution_time_ms / 1000.0)
               : (std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() / 1000.0));

        // Prepare kernel parameter string for logging
        auto kparams = k.second->getParameters();
        QStringList kparamsStr; for (double v : kparams) kparamsStr << QString::number(v, 'g', 6);
        QMetaObject::invokeMethod(this, [this, k, tc, result, cpuSec, completed, total, kparamsStr]() {
    int row = m_resultsTable->rowCount();
    m_resultsTable->insertRow(row);
          m_resultsTable->setItem(row, 0, new QTableWidgetItem(k.first));
          // Parameters column is embedded in name; also show explicit params if available
          QString params;
          auto ps = k.second->getParameters();
          if (!ps.empty()) {
            if (k.first.startsWith("Derbal")) params = QString("q=%1").arg(ps[0],0,'f',2);
            else if (k.first.startsWith("Bachir")) params = QString("m=%1").arg(ps[0],0,'f',2);
            else if (k.first.startsWith("Wu")) params = QString("p=%1").arg(ps[0],0,'f',2);
            else if (k.first.startsWith("IJNAA")) params = QString("q=%1").arg(ps[0],0,'f',2);
          } else params = "(none)";
          m_resultsTable->setItem(row, 1, new QTableWidgetItem(params));
          m_resultsTable->setItem(row, 2, new QTableWidgetItem(QString::fromStdString(tc.name)));
          m_resultsTable->setItem(row, 3, new QTableWidgetItem(QString::number(result.iterations)));
          m_resultsTable->setItem(row, 4, new QTableWidgetItem(QString::number(cpuSec, 'f', 3) + "s"));
          m_resultsTable->setItem(row, 5, new QTableWidgetItem(QString::number(result.duality_gap, 'e', 2)));
          QString status = result.converged ? "✅" : "Error - " + QString::fromStdString(result.convergence_info);
          m_resultsTable->setItem(row, 6, new QTableWidgetItem(status));

          // Update CPU time chart (Bachir results only, per request)
          #ifdef CPU_CHARTS_ENABLED
          if (m_cpuChart) {
            bool ok = false; QString tstr = m_resultsTable->item(row, 4)->text();
            double secs = tstr.endsWith("s") ? tstr.left(tstr.length()-1).toDouble(&ok) : tstr.toDouble(&ok);
            if (ok) {
              // Use kernel name as category so we always see data
              QString category = k.first;
              m_cpuChart->upsertSample(category, secs);
            }
          }
          #endif
          // Switch from indeterminate to determinate on first completion
          if (m_progressBar->minimum() == 0 && m_progressBar->maximum() == 0) {
            m_progressBar->setRange(0, std::max(1,total));
          }
          int next = completed + 1; m_progressBar->setValue(next);
          m_progressText->setText(formatProgressBar(next, total));
          if (m_logTextEdit) {
            m_logTextEdit->append(QString("Result → Kernel: %1, iters=%2, time=%3s, gap=%4, status=%5, params=[%6]")
              .arg(k.first)
              .arg(result.iterations)
              .arg(cpuSec, 0, 'f', 3)
              .arg(result.duality_gap, 0, 'e', 3)
              .arg(result.converged ? "Converged" : QString::fromStdString(result.convergence_info))
              .arg(kparamsStr.join(", "))
            );
            // Additional iteration statistics as requested
            m_logTextEdit->append(QString("   → Outer L=%1, Inner total=%2, K_avg=%3")
              .arg(result.outer_iterations)
              .arg(result.inner_iterations_total)
              .arg(result.avg_inner_per_outer, 0, 'f', 3));
          }
        });

        completed++;
      }
    }

    QMetaObject::invokeMethod(this, [this, total]() {
  m_isRunning = false;
  m_runAllButton->setEnabled(true);
  m_runSingleButton->setEnabled(true);
      m_cancelButton->setEnabled(false);
  m_progressBar->setVisible(false);
      // Completion summary
      int okCount = 0; for (int r=0; r<m_resultsTable->rowCount(); ++r) if (m_resultsTable->item(r,6) && m_resultsTable->item(r,6)->text().startsWith("✅")) okCount++;
      if (total>0) m_statusLabel->setText(QString("Complete - %1/%2 combinations successful").arg(okCount).arg(total));
      // Set status banner color based on outcome
      if (total > 0 && okCount == total) m_statusLabel->setProperty("status", "success");
      else m_statusLabel->setProperty("status", "error");
      m_statusLabel->style()->unpolish(m_statusLabel); m_statusLabel->style()->polish(m_statusLabel);
    });
  });
}

void MainWindow::runSingleAlgorithm() {
  if (m_isRunning)
    return;

  m_isRunning = true;
  m_runAllButton->setEnabled(false);
  m_runSingleButton->setEnabled(false);
  m_cancelButton->setEnabled(true);
  m_progressBar->setVisible(true);
  m_progressBar->setRange(0, 1);
  m_progressBar->setValue(0);
  m_statusLabel->setText("Running selected algorithm...");
  m_statusLabel->setProperty("status", "running"); m_statusLabel->style()->unpolish(m_statusLabel); m_statusLabel->style()->polish(m_statusLabel);

  m_logTextEdit->clear();
  m_logTextEdit->append(QString("Running single solve..."));
  m_logTextEdit->append("-------------------------------------------");

  // Create TestCase from current matrix input data
  m_currentTestCase = createTestCaseFromMatrixInput();

  // Log the problem being solved
  m_logTextEdit->append(
      QString("Problem: %1").arg(QString::fromStdString(m_currentTestCase.name)));
  m_logTextEdit->append(QString("Size: n=%1, m=%2")
                            .arg(m_currentTestCase.problem_size)
                            .arg(m_currentTestCase.A_matrices.size()));
  m_logTextEdit->append(QString("Parameters: ε=%1, θ=%2, τ=%3")
                            .arg(m_currentTestCase.tolerance)
                            .arg(m_currentTestCase.theta)
                            .arg(m_currentTestCase.tau));
  m_logTextEdit->append("");

  m_cancel.store(false);
  auto task = [this]() {
    auto kernel = createSelectedKernel();
    IpmSettings settings;
    settings.epsilon = m_currentTestCase.tolerance;
    settings.theta = m_currentTestCase.theta;
    settings.tau = m_currentTestCase.tau;
    settings.cancel = &m_cancel;
    AlgorithmResult result;
    try {
      result = m_engine.solve(m_currentTestCase, *kernel, settings);
  } catch (const std::exception &e) {
    result.converged = false;
    result.convergence_info = std::string("Error: ") + e.what();
  }

    QMetaObject::invokeMethod(this, [this, kernel, result]() {
      QString algorithmName = "Primal-Dual IPM (NT)";
      QString kernelName = QString::fromStdString(kernel->getName());
      // Prepare parameters string
      QStringList paramsStr; for (double v : kernel->getParameters()) paramsStr << QString::number(v, 'g', 6);
      QString params = paramsStr.isEmpty() ? "(none)" : QString("[%1]").arg(paramsStr.join(", "));

      // One-row results matching headers:
      // 0 Algorithm | 1 Parameters | 2 Test Case | 3 Iterations | 4 CPU Time | 5 Duality Gap | 6 Status
      m_resultsTable->setRowCount(1);
      m_resultsTable->setItem(0, 0, new QTableWidgetItem(algorithmName + " / " + kernelName));
      m_resultsTable->setItem(0, 1, new QTableWidgetItem(params));
      m_resultsTable->setItem(0, 2, new QTableWidgetItem(QString::fromStdString(m_currentTestCase.name)));
      m_resultsTable->setItem(0, 3, new QTableWidgetItem(QString::number(result.iterations)));
      // Prefer CPU seconds in table if available, else wall seconds, else 0
      double seconds = result.cpu_time_ms > 0.0 ? (result.cpu_time_ms / 1000.0)
                        : (result.execution_time_ms > 0.0 ? (result.execution_time_ms / 1000.0) : 0.0);
      m_resultsTable->setItem(0, 4, new QTableWidgetItem(QString::number(seconds, 'f', 3) + "s"));
      m_resultsTable->setItem(0, 5, new QTableWidgetItem(QString::number(result.duality_gap, 'e', 3)));
      QString status = result.converged ? "Converged ✅" : "Failed";
      if (!result.convergence_info.empty()) status += QString(" (%1)").arg(QString::fromStdString(result.convergence_info));
      m_resultsTable->setItem(0, 6, new QTableWidgetItem(status));

      m_progressBar->setValue(1);
      m_logTextEdit->append(QString("%1 completed").arg(algorithmName));
      m_logTextEdit->append(QString("   → %1 iterations, %2 s, gap=%3, objective: %4")
                              .arg(result.iterations)
                              .arg(seconds, 0, 'f', 3)
                              .arg(result.duality_gap, 0, 'e', 3)
                              .arg(result.primal_objective, 0, 'f', 6));
      if (!result.converged) {
        m_logTextEdit->append(QString("   → Status: %1").arg(QString::fromStdString(result.convergence_info)));
      }
      statusBar()->showMessage("Single algorithm execution completed", 3000);

      m_isRunning = false;
      m_runAllButton->setEnabled(true);
      m_runSingleButton->setEnabled(true);
      m_cancelButton->setEnabled(false);
      m_progressBar->setVisible(false);
      m_statusLabel->setText(m_cancel.load() ? "Cancelled" : (result.converged ? "Completed" : "Failed"));
      if (m_cancel.load()) m_statusLabel->setProperty("status", "error");
      else m_statusLabel->setProperty("status", result.converged ? "success" : "error");
      m_statusLabel->style()->unpolish(m_statusLabel); m_statusLabel->style()->polish(m_statusLabel);
    });
  };

  m_future = QtConcurrent::run(task);
}

void MainWindow::loadTestCase() {
  int selected = m_testCaseCombo->currentIndex();
  if (selected < 0 || selected >= static_cast<int>(m_testCases.size()))
    return;

  const TestCase &testCase = m_testCases[selected];

  QString info = QString("%1\n\n").arg(QString::fromStdString(testCase.name));
  info += QString("Problem Configuration:\n");
  info += QString("   • Size: n=%1, m=%2\n")
              .arg(testCase.problem_size)
              .arg(testCase.A_matrices.size());
  info += QString("   • Tolerance: ε=%1\n").arg(testCase.tolerance);
  info += QString("   • Theta: θ=%1\n").arg(testCase.theta);
  info += QString("   • Tau: τ=%1\n").arg(testCase.tau);
  info += QString("   • Expected optimal: %1\n\n")
              .arg(testCase.expected_optimal_value);

  // Display A₁ matrix
  if (!testCase.A_matrices.empty() && !testCase.A_matrices[0].empty()) {
    info += "A₁ Constraint Matrix:\n";
    for (const auto &row : testCase.A_matrices[0]) {
      info += "   [";
      for (size_t j = 0; j < row.size(); ++j) {
        info += QString("%1").arg(row[j], 4, 'f', 0);
        if (j < row.size() - 1)
          info += " ";
      }
      info += "]\n";
    }
    info += "\n";
  }

  // Display b vector
  if (!testCase.b_vector.empty()) {
    info += "b Vector: [";
    for (size_t i = 0; i < testCase.b_vector.size(); ++i) {
      info += QString::number(testCase.b_vector[i]);
      if (i < testCase.b_vector.size() - 1)
        info += ", ";
    }
    info += "]\n";
  }

  m_matrixDisplay->setPlainText(info);

  // Update parameter spinboxes
  m_toleranceSpinBox->setValue(testCase.tolerance);
  m_thetaSpinBox->setValue(testCase.theta);
  m_tauSpinBox->setValue(testCase.tau);

  // Update MatrixInputWidget with loaded data
  m_matrixInputWidget->setProblemSize(testCase.problem_size,
                                      static_cast<int>(testCase.A_matrices.size()));
  m_matrixInputWidget->setAMatrices(testCase.A_matrices);
  m_matrixInputWidget->setBVector(testCase.b_vector);
  m_matrixInputWidget->setCMatrix(testCase.C_matrix);
  m_matrixInputWidget->setInitialX(testCase.initial_X);
  m_matrixInputWidget->setInitialY(testCase.initial_y);
  m_matrixInputWidget->setInitialS(testCase.initial_S);

  statusBar()->showMessage(QString("Loaded: %1").arg(QString::fromStdString(testCase.name)), 3000);
}

void MainWindow::applyCustomMatrixData() {
  // Get data from MatrixInputWidget
  auto problemSize = m_matrixInputWidget->getProblemSize();
  int n = problemSize.first;
  int m = problemSize.second;

  auto aMatrices = m_matrixInputWidget->getAMatrices();
  auto bVector = m_matrixInputWidget->getBVector();
  auto cMatrix = m_matrixInputWidget->getCMatrix();

  // Get current parameter values from spinboxes
  double newTolerance = m_toleranceSpinBox->value();
  double newTheta = m_thetaSpinBox->value();
  double newTau = m_tauSpinBox->value();

  // Update matrix display with custom data
  QString info = QString("Custom Matrix Data Applied:\n\n");
  info += QString("Problem Configuration:\n");
  info += QString("   • Size: n=%1, m=%2\n").arg(n).arg(m);
  info += QString("   • Tolerance: ε=%1\n").arg(newTolerance);
  info += QString("   • Theta: θ=%1\n").arg(newTheta);
  info += QString("   • Tau: τ=%1\n").arg(newTau);
  info += QString("   • Expected optimal: -1.0957\n\n");

  // Display A₁ matrix
  if (!aMatrices.empty() && !aMatrices[0].empty()) {
    info += "A₁ Constraint Matrix:\n";
    for (const auto &row : aMatrices[0]) {
      info += "   [";
      for (size_t j = 0; j < row.size(); ++j) {
        info += QString("%1").arg(row[j], 4, 'f', 0);
        if (j < row.size() - 1)
          info += " ";
      }
      info += "]\n";
    }
    info += "\n";
  }

  // Display b vector
  if (!bVector.empty()) {
    info += "b Vector: [";
    for (size_t i = 0; i < bVector.size(); ++i) {
      info += QString::number(bVector[i]);
      if (i < bVector.size() - 1)
        info += ", ";
    }
    info += "]\n";
  }

  m_matrixDisplay->setPlainText(info);

  statusBar()->showMessage(
      QString("Custom Data Applied: n=%1, m=%2").arg(n).arg(m), 3000);
}

void MainWindow::onMatrixDataChanged() {
  // This signal is emitted when the user changes the matrix data in
  // MatrixInputWidget
  statusBar()->showMessage(
      "Matrix data modified - click 'Apply Custom Data' to use it", 2000);
}

void MainWindow::onProblemSizeChanged(int n, int m) {
  // This signal is emitted when the user changes the problem size in
  // MatrixInputWidget
  statusBar()->showMessage(
      QString("Problem size changed to n=%1, m=%2").arg(n).arg(m), 2000);
}

// simulateAlgorithmExecution removed
// executeNextAlgorithm removed

void MainWindow::onTableCellClicked(int row, int column) {
  Q_UNUSED(column)

  if (row < 0 || row >= m_resultsTable->rowCount())
    return;

  QTableWidgetItem *algorithmItem = m_resultsTable->item(row, 0);
  if (algorithmItem) {
    QString algorithmName = algorithmItem->text();
    statusBar()->showMessage(
        QString("Selected: %1 - Click Details for more information")
            .arg(algorithmName),
        3000);
  }
}

void MainWindow::showAlgorithmDialog(const QString &algorithmName) {
  // Create and show the algorithm dialog with dynamic content
  AlgorithmDialog *dialog = new AlgorithmDialog(this);

  // Extract algorithm index from the name (e.g., "Algorithm 2" -> index 1)
  int algorithmIndex = -1;
  if (algorithmName.contains("Algorithm")) {
    QString numStr = algorithmName.split(" ").last();
    algorithmIndex = numStr.toInt() - 1; // Convert to 0-based index
  }

  // If we have a valid algorithm index and corresponding kernel, use dynamic
  // content
  if (algorithmIndex >= 0 &&
      algorithmIndex < static_cast<int>(m_kernels.size())) {
    // Get the actual kernel function
    auto kernel = m_kernels[algorithmIndex];

    // Create dynamic algorithm name based on the kernel
    QString dynamicName = QString("Algorithm %1: %2")
                              .arg(algorithmIndex + 1)
                              .arg(QString::fromStdString(kernel->getName()));

    // Use the new dynamic method with actual kernel function
    dialog->showAlgorithmDetails(dynamicName, kernel);
  } else {
    // Fallback to the original name
    dialog->showAlgorithmDetails(algorithmName);
  }

  dialog->exec();
  dialog->deleteLater();
}

void MainWindow::updateCharts() {
  // This would update the charts - placeholder for now
  statusBar()->showMessage("Charts updated", 2000);
}

QString MainWindow::formatProgressBar(int completed, int total) const {
  if (total <= 0) total = 1;
  if (completed < 0) completed = 0;
  if (completed > total) completed = total;
  int percent = (completed * 100) / total;
  return QString("%1/%2 (%3%)").arg(completed).arg(total).arg(percent);
}

std::vector<std::pair<QString, std::shared_ptr<KernelBase>>> MainWindow::getSelectedKernels() const {
  std::vector<std::pair<QString, std::shared_ptr<KernelBase>>> out;
  if (!m_kernelList) return out;
  auto items = m_kernelList->selectedItems();
  // Parameter values: only apply editor values if single selection
  std::vector<double> editorVals;
  if (items.size() == 1) editorVals = m_parameterEditor ? m_parameterEditor->currentValues() : std::vector<double>{};
  for (auto *it : items) {
    QString name = it->text();
    for (const auto &info : m_kernelManager.getKernels()) {
      if (info.name == name) {
        // Prefer per-kernel persisted values if present; else editor (single select); else defaults
        std::vector<double> values;
        auto kv = g_kernelEditedValues.find(name);
        if (kv != g_kernelEditedValues.end()) values = kv->second;
        else values = editorVals;
        if (values.empty()) { for (const auto &p : info.parameters) values.push_back(p.defaultVal); }
        auto k = m_kernelManager.createKernel(info, values);
        out.emplace_back(name, k);
      }
    }
  }
  // If none selected, pick first kernel with defaults
  // Removed silent fallback; we now require explicit selection to avoid confusion
  return out;
}

std::vector<TestCase> MainWindow::getSelectedTests() const {
  std::vector<TestCase> out;
  if (!m_testList) return out;
  auto items = m_testList->selectedItems();
  if (items.empty()) return out;
  for (auto *it : items) {
    QString name = it->text();
    for (const auto &tc : m_testCases) {
      if (QString::fromStdString(tc.name) == name) { out.push_back(tc); break; }
    }
  }
  return out;
}

void MainWindow::exportResults() {
  // Export visible table to CSV
  QString fileName = QFileDialog::getSaveFileName(this, "Export Results to CSV", "results.csv", "CSV Files (*.csv)");
  if (fileName.isEmpty()) return;
  QFile file(fileName);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QMessageBox::warning(this, "Export Results", "Failed to open file for writing.");
    return;
  }
  QTextStream out(&file);
  // header
  QStringList headers;
  for (int c = 0; c < m_resultsTable->columnCount(); ++c) headers << m_resultsTable->horizontalHeaderItem(c)->text();
  out << headers.join(',') << "\n";
  // rows
  for (int r = 0; r < m_resultsTable->rowCount(); ++r) {
    QStringList row;
    for (int c = 0; c < m_resultsTable->columnCount(); ++c) {
      QTableWidgetItem *item = m_resultsTable->item(r, c);
      row << (item ? item->text() : "");
    }
    out << row.join(',') << "\n";
  }
  file.close();
  QMessageBox::information(this, "Export Results", "Exported results to " + fileName);
}

void MainWindow::showAlgorithmInfo() {
  QMessageBox::about(this, "About Kernel Function Comparison Tool",
                     "<b>Kernel Function Comparison Tool v1.0</b><br><br>"
                     "A comprehensive analysis tool for comparing six "
                     "different primal-dual interior-point methods:<br><br>"
                     "<b>Algorithms Implemented:</b><br>"
                     "• Algorithm 1: Generic PD IPM for Linear Optimization<br>"
                     "• Algorithm 2: PD IPM for Semidefinite Optimization<br>"
                     "• Algorithm 3: PD IPM with Log-Kernel<br>"
                     "• Algorithm 4: PD for Conic Quadratic SDO<br>"
                     "• Algorithm 5: Derbal & Kebbiche Method<br>"
                     "• Algorithm 6: Bachir φₘ(t) Family<br><br>"
                     "<b>Features:</b><br>"
                     "• Interactive matrix input and parameter control<br>"
                     "• Real-time algorithm execution and comparison<br>"
                     "• Advanced visualization with charts and graphs<br>"
                     "• Detailed algorithm explanations and formulas<br>"
                     "• Professional results export capabilities<br><br>"
                     "Built with Qt6 and modern C++17<br>"
                     "For research and educational purposes");
}
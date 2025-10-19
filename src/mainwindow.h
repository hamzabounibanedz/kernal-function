#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "engine/ipm_engine.h"
#include "utils/test_data.h"
#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMainWindow>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QSplitter>
#include <QTabWidget>
#include <QTableWidget>
#include <QTextEdit>
#include <QTimer>
#include <QVBoxLayout>
#include <QThread>
#include <QtConcurrent>
#include <QCheckBox>
#include <QSlider>
#include <QListWidget>
#ifdef CPU_CHARTS_ENABLED
#include "gui/cpu_time_chart_widget.h"
#endif
#include <atomic>
#include <memory>
#include <vector>

// Forward declarations - these are other classes we'll use
// We declare them here so the compiler knows they exist, but we don't include
// their full definitions yet
class ChartWidget;       // For displaying performance charts
class AlgorithmDialog;   // Dialog box for algorithm information
class MatrixInputWidget; // Widget for entering matrix data
struct TestCase;         // Represents a test problem
class KernelBase;        // Base class for kernel functions
class AlgorithmBase;     // Base class for optimization algorithms

// New UI helpers
#include "ui/kernel_manager.h"
#include "ui/parameter_editor_widget.h"

/**
 * MainWindow - This is the heart of our application
 *
 * This class creates and manages the entire user interface for comparing kernel
 * functions. It's organized into several panels:
 * - Left panel: Input controls and parameter settings
 * - Right panel: Results display and charts
 * - Bottom: Status and progress information
 *
 * The interface is designed to be intuitive - you can select test problems,
 * adjust parameters, run algorithms, and see the results all in one place.
 */
class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  // Constructor - sets up the entire user interface
  explicit MainWindow(QWidget *parent = nullptr);
  // Destructor - cleans up when the window is closed
  ~MainWindow();

private slots:
  // These are functions that get called when the user interacts with the
  // interface They're connected to buttons, menus, and other UI elements

  // Runs all available algorithms on the current test problem
  void runAllAlgorithms();
  // Runs just one specific algorithm
  void runSingleAlgorithm();
  // Loads a test case from the dropdown menu
  void loadTestCase();
  // Shows detailed information about an algorithm
  void showAlgorithmInfo();
  // Handles clicks on the results table
  void onTableCellClicked(int row, int column);
  // Updates the performance charts with new data
  void updateCharts();
  // Exports results to a file
  void exportResults();
  // Applies custom matrix data from MatrixInputWidget
  void applyCustomMatrixData();
  // Handles matrix data changes from MatrixInputWidget
  void onMatrixDataChanged();
  // Handles problem size changes from MatrixInputWidget
  void onProblemSizeChanged(int n, int m);
  // Shows a dialog with information about a specific algorithm
  void showAlgorithmDialog(const QString &algorithmName);
  // Handle kernel selection change
  void onKernelChanged(int index);

private:
  // These are helper functions that set up different parts of the interface

  // Creates the entire user interface layout
  void setupUI();
  // Sets up the menu bar at the top
  void setupMenuBar();
  // Creates the left panel with input controls
  void setupInputPanel();
  // Creates the control panel with buttons and progress bar
  void setupControlPanel();
  // Creates the right panel with results and charts
  void setupResultsPanel();
  // Loads test data and initializes algorithms
  void initializeData();
  // Sets up the headers for the results table
  void setupTableHeaders();
  // Creates a TestCase from MatrixInputWidget data
  TestCase createTestCaseFromMatrixInput();
  // Create kernel by selection
  std::shared_ptr<KernelBase> createSelectedKernel() const;

  // === USER INTERFACE COMPONENTS ===

  // Main layout components
  QWidget *m_centralWidget;  // The main widget that contains everything
  QHBoxLayout *m_mainLayout; // Horizontal layout for left/right panels
  QSplitter *m_mainSplitter; // Allows user to resize left/right panels

  // === LEFT PANEL - Input and Control ===
  QWidget *m_leftPanel;      // Container for the left side
  QVBoxLayout *m_leftLayout; // Vertical layout for left panel contents

  // Left tabs (Benchmark Setup / Problem Editor / Run Settings)
  QTabWidget *m_leftTabs {nullptr};
  // Benchmark Setup
  QListWidget *m_kernelList {nullptr};
  QListWidget *m_testList {nullptr};
  ParameterEditorWidget *m_parameterEditor {nullptr};

  // Input Section - where users select test problems and enter data
  QGroupBox *m_inputGroup;                // Group box with title "Input"
  QVBoxLayout *m_inputLayout;             // Layout for input controls
  QComboBox *m_testCaseCombo;             // Dropdown to select test problems
  QPushButton *m_loadTestButton;          // Button to load selected test case
  QTextEdit *m_matrixDisplay;             // Shows the current matrix data
  MatrixInputWidget *m_matrixInputWidget; // Custom widget for matrix input

  // Control Section - where users set parameters and run algorithms
  QGroupBox *m_controlGroup;          // Group box with title "Control"
  QVBoxLayout *m_controlLayout;       // Layout for control elements
  QDoubleSpinBox *m_toleranceSpinBox; // Input for convergence tolerance
  QDoubleSpinBox *m_thetaSpinBox;     // Input for θ parameter
  QDoubleSpinBox *m_tauSpinBox;       // Input for τ parameter
  QComboBox *m_kernelCombo;           // Legacy single kernel selection
  QLabel *m_kernelParamLabel;         // Legacy param label
  QDoubleSpinBox *m_kernelParamSpin;  // Legacy param value
  QPushButton *m_runAllButton;        // Button to run all algorithms
  QPushButton *m_runSingleButton;     // Button to run one algorithm
  QPushButton *m_cancelButton;        // Button to cancel running computations
  QProgressBar *m_progressBar;        // Shows algorithm execution progress
  QLabel *m_statusLabel;              // Shows current status messages
  QLabel *m_progressText;             // Detailed progress text as per spec

  // Kernel multi-selection per spec
  QGroupBox *m_kernelSelectionGroup;
  QCheckBox *m_chkTouil;
  QCheckBox *m_chkDerbal;
  QCheckBox *m_chkBachir;
  QCheckBox *m_chkWu;
  QCheckBox *m_chkIjnna;
  QDoubleSpinBox *m_qDerbal; QSlider *m_qDerbalSlider;
  QDoubleSpinBox *m_mBachir; QSlider *m_mBachirSlider;
  QDoubleSpinBox *m_pWu; QSlider *m_pWuSlider; QLabel *m_aWuLabel;
  QDoubleSpinBox *m_qIjnna; QSlider *m_qIjnnaSlider;

  // Test case selection per spec
  QGroupBox *m_testsSelectionGroup;
  QCheckBox *m_chkAllTests;
  QCheckBox *m_chkTouilP1; QCheckBox *m_chkTouilP2;
  QCheckBox *m_chkDerbal5x5;
  QCheckBox *m_chkBachirEx1; QCheckBox *m_chkBachirEx2;
  QCheckBox *m_chkWu41; QCheckBox *m_chkWu42;
  QCheckBox *m_chkIjnna5x5;

  // === RIGHT PANEL - Results and Charts ===
  QWidget *m_rightPanel;          // Container for the right side
  QVBoxLayout *m_rightLayout;     // Vertical layout for right panel contents
  QTabWidget *m_resultsTabWidget; // Tabbed interface for different result views

  // Results Table Tab - shows numerical results in a table format
  QWidget *m_tableTab;          // Tab widget for the table
  QVBoxLayout *m_tableLayout;   // Layout for table tab contents
  QTableWidget *m_resultsTable; // Table showing algorithm results
  QPushButton *m_exportButton;  // Legacy export button
  QPushButton *m_exportCsvButton;  // CSV
  QPushButton *m_exportPngButton;  // PNG
  QPushButton *m_exportPdfButton;  // PDF

  // CPU Time chart tab
  QWidget *m_chartTab {nullptr};
  QVBoxLayout *m_chartLayout {nullptr};
#ifdef CPU_CHARTS_ENABLED
  CpuTimeChartWidget *m_cpuChart {nullptr};
#else
  void *m_cpuChart {nullptr};
#endif

  // Charts Tab - shows performance charts and graphs
  QWidget *m_chartsTab;        // Tab widget for charts
  QVBoxLayout *m_chartsLayout; // Layout for charts tab contents
  ChartWidget *m_chartWidget;  // Custom widget for displaying charts

  // Log Tab - shows detailed execution log
  QWidget *m_logTab;        // Tab widget for log
  QVBoxLayout *m_logLayout; // Layout for log tab contents
  QTextEdit *m_logTextEdit; // Text area showing execution details

  // === DATA STORAGE ===

  // These store the actual data that the application works with
  std::vector<TestCase> m_testCases; // Available test problems
  std::vector<std::shared_ptr<KernelBase>>
      m_kernels; // Available kernel functions
  IpmEngine m_engine;                 // Core IPM engine
  KernelManager m_kernelManager;      // Data-driven kernel registry

  // === APPLICATION STATE ===

  // These track the current state of the application
  bool m_isRunning;           // Whether an algorithm is currently running
  TestCase m_currentTestCase; // Current test case being executed
  QFuture<void> m_future;     // Background run handle
  std::atomic<bool> m_cancel{false};

  // Helpers for new UI flow
  QString formatProgressBar(int completed, int total) const;
  std::vector<std::pair<QString, std::shared_ptr<KernelBase>>> getSelectedKernels() const;
  std::vector<TestCase> getSelectedTests() const;
};

#endif // MAINWINDOW_H
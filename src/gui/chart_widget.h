#ifndef CHART_WIDGET_H
#define CHART_WIDGET_H

#include "../kernels/kernel_base.h"
#include <QtCharts/QBarCategoryAxis>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <memory>
#include <vector>

QT_BEGIN_NAMESPACE
QT_END_NAMESPACE

using namespace QtCharts;

/**
 * ChartWidget - Visualizes algorithm performance and kernel functions
 *
 * This widget creates charts and graphs to help users understand how
 * different algorithms perform. It can show:
 * - How kernel functions behave over different input values
 * - How many iterations each algorithm took to find a solution
 * - How long each algorithm took to run
 * - How good the final solutions were
 *
 * The charts make it easy to compare different approaches and see
 * which algorithms work best for different types of problems.
 */
class ChartWidget : public QWidget {
  Q_OBJECT

public:
  // Constructor - creates the chart widget with controls
  explicit ChartWidget(QWidget *parent = nullptr);

  // Destructor - cleans up when the widget is destroyed
  ~ChartWidget();

  // === CHART DATA UPDATES ===

  /**
   * Update the kernel comparison chart
   *
   * This shows how different kernel functions behave over a range
   * of input values. It helps users understand the mathematical
   * properties of each kernel.
   *
   * @param kernels Vector of kernel functions to compare
   */
  void updateKernelComparison(
      const std::vector<std::shared_ptr<KernelBase>> &kernels);

  /**
   * Update the algorithm results chart
   *
   * This shows performance metrics for different algorithms:
   * - How many iterations each took
   * - How long each took to run
   * - How good the final solution was
   *
   * @param algorithms Names of the algorithms
   * @param iterations Number of iterations for each algorithm
   * @param times Execution time for each algorithm
   * @param objectives Final objective value for each algorithm
   */
  void updateAlgorithmResults(const std::vector<QString> &algorithms,
                              const std::vector<int> &iterations,
                              const std::vector<double> &times,
                              const std::vector<double> &objectives);

  /**
   * Update just the iteration comparison chart
   *
   * This focuses specifically on how many iterations each
   * algorithm needed to find a solution.
   *
   * @param algorithms Names of the algorithms
   * @param iterations Number of iterations for each algorithm
   */
  void updateIterationChart(const std::vector<QString> &algorithms,
                            const std::vector<int> &iterations);

public slots:
  // These functions are called when the user interacts with the interface

  /**
   * Handle chart type selection changes
   *
   * When the user selects a different chart type from the dropdown,
   * this function updates the display to show the selected chart.
   */
  void onChartTypeChanged();

  /**
   * Export the current chart to a file
   *
   * Allows users to save charts as images for reports or presentations.
   */
  void exportChart();

private:
  // === UI SETUP FUNCTIONS ===

  // Create the overall layout and controls
  void setupUI();

  // === CHART CREATION FUNCTIONS ===

  // Create a chart comparing different kernel functions
  void createKernelComparisonChart();

  // Create a bar chart comparing iteration counts
  void createIterationComparisonChart();

  // Create a bar chart comparing execution times
  void createTimeComparisonChart();

  // Create a bar chart comparing objective values
  void createObjectiveComparisonChart();

  // Apply consistent styling to charts
  void setupChartStyling(QChart *chart);

  // === USER INTERFACE COMPONENTS ===

  // Main layout that contains everything
  QVBoxLayout *m_mainLayout;

  // Horizontal layout for controls (chart type selector, export button)
  QHBoxLayout *m_controlLayout;

  // Title label for the chart
  QLabel *m_titleLabel;

  // Dropdown to select different chart types
  QComboBox *m_chartTypeCombo;

  // Button to export the chart as an image
  QPushButton *m_exportButton;

  // === CHART COMPONENTS ===

  // The main chart view that displays the chart
  QChartView *m_chartView;

  // The chart object that contains the data and axes
  QChart *m_chart;

  // Series for bar charts (used for comparing discrete values)
  QBarSeries *m_barSeries;

  // Series for line charts (used for continuous functions)
  QLineSeries *m_lineSeries;

  // Y-axis (vertical axis) for numerical values
  QValueAxis *m_axisY;

  // X-axis (horizontal axis) for categories or values
  QBarCategoryAxis *m_axisX;

  // === DATA STORAGE ===

  // Kernel functions for comparison charts
  std::vector<std::shared_ptr<KernelBase>> m_kernels;

  // Algorithm names for performance charts
  std::vector<QString> m_algorithmNames;

  // Performance data for different algorithms
  std::vector<int> m_iterationData;    // How many iterations each took
  std::vector<double> m_timeData;      // How long each took to run
  std::vector<double> m_objectiveData; // How good the final solution was

  // === CHART TYPE MANAGEMENT ===

  // Different types of charts we can display
  enum ChartType {
    KernelComparison,    // Compare kernel function behaviors
    IterationComparison, // Compare iteration counts
    TimeComparison,      // Compare execution times
    ObjectiveComparison  // Compare final objective values
  };

  // Which chart type is currently being displayed
  ChartType m_currentChartType;
};

#endif // CHART_WIDGET_H
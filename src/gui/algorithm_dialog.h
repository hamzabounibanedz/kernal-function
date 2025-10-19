#ifndef ALGORITHM_DIALOG_H
#define ALGORITHM_DIALOG_H

#include <QDialog>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QScrollArea>
#include <QString>
#include <QTabWidget>
#include <QTextEdit>
#include <QVBoxLayout>
#include <memory>

#include "kernels/kernel_base.h"

/**
 * @brief Dialog for displaying detailed algorithm descriptions and formulas
 *
 * This dialog shows when users click on algorithm names in the results table.
 * It provides comprehensive information about each algorithm including:
 * - Complete algorithm pseudocode
 * - Kernel function formulas with derivatives
 * - Theoretical background and complexity
 * - Numerical examples from papers
 * - Implementation details
 */
class AlgorithmDialog : public QDialog {
  Q_OBJECT

public:
  explicit AlgorithmDialog(QWidget *parent = nullptr);
  ~AlgorithmDialog();

  /**
   * Show algorithm details for specific algorithm
   * @param algorithmName Name of the algorithm (e.g., "Algorithm 1", "Bachir")
   */
  void showAlgorithmDetails(const QString &algorithmName);

  /**
   * Show algorithm details with actual kernel function for dynamic content
   * @param algorithmName Name of the algorithm
   * @param kernel The actual kernel function to get dynamic information from
   */
  void showAlgorithmDetails(const QString &algorithmName,
                            std::shared_ptr<KernelBase> kernel);

private slots:
  void onCopyFormula();
  void onExportDetails();
  void onRunExample();

private:
  void setupUI();
  void setupTabs();

  /**
   * @brief Get complete details for each algorithm matching user specifications
   */
  QString getAlgorithm1Details() const; // Trigonometric Kernel
  QString getAlgorithm2Details() const; // Exponential-Parametric Kernel
  QString getAlgorithm3Details() const; // Parameterized Log Kernel
  QString getAlgorithm4Details() const; // Parametric Family Kernel
  QString getAlgorithm5Details() const; // Derbal & Kebbiche
  QString getAlgorithm6Details() const; // Bachir φₘ(t) Family

  /**
   * @brief Format mathematical formulas with proper LaTeX-style notation
   */
  QString formatFormula(const QString &formula) const;
  QString formatAlgorithmSteps(const QString &steps) const;
  QString formatNumericalExample(const QString &example) const;

  // UI Components
  QVBoxLayout *m_mainLayout;
  QTabWidget *m_tabWidget;

  // Algorithm Tab
  QWidget *m_algorithmTab;
  QTextEdit *m_algorithmText;

  // Kernel Tab
  QWidget *m_kernelTab;
  QTextEdit *m_kernelText;

  // Example Tab
  QWidget *m_exampleTab;
  QTextEdit *m_exampleText;

  // Control buttons
  QHBoxLayout *m_buttonLayout;
  QPushButton *m_copyButton;
  QPushButton *m_exportButton;
  QPushButton *m_runExampleButton;
  QPushButton *m_closeButton;

  // Current algorithm
  QString m_currentAlgorithm;
};

#endif // ALGORITHM_DIALOG_H
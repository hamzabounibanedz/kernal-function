#ifndef MATRIX_INPUT_WIDGET_H
#define MATRIX_INPUT_WIDGET_H

#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

/**
 * @brief Widget for inputting matrices and vectors for optimization problems
 */
class MatrixInputWidget : public QWidget {
  Q_OBJECT

public:
  explicit MatrixInputWidget(QWidget *parent = nullptr);
  ~MatrixInputWidget();

  // Matrix data access
  std::vector<std::vector<std::vector<double>>> getAMatrices() const;
  std::vector<double> getBVector() const;
  std::vector<std::vector<double>> getCMatrix() const;
  std::vector<std::vector<double>> getInitialX() const;
  std::vector<double> getInitialY() const;
  std::vector<std::vector<double>> getInitialS() const;

  // Set data from test cases
  void
  setAMatrices(const std::vector<std::vector<std::vector<double>>> &matrices);
  void setBVector(const std::vector<double> &vector);
  void setCMatrix(const std::vector<std::vector<double>> &matrix);
  void setInitialX(const std::vector<std::vector<double>> &matrix);
  void setInitialY(const std::vector<double> &vector);
  void setInitialS(const std::vector<std::vector<double>> &matrix);

  // Problem size
  void setProblemSize(int n, int m);
  std::pair<int, int> getProblemSize() const;

signals:
  void dataChanged();
  void problemSizeChanged(int n, int m);

private slots:
  void onProblemSizeChanged();
  void onCellChanged(int row, int col);
  void clearAllData();
  void loadFromText();
  void exportToText();

private:
  void setupUI();
  void setupSizeControls();
  void setupMatrixTables();
  void setupTextInput();
  void resizeTables();
  void updateTableSizes();
  void updateConstraintTables();
  void setDefaultValues();

  // UI components
  QVBoxLayout *m_mainLayout;
  QTabWidget *m_tabWidget;

  // Problem size controls
  QHBoxLayout *m_sizeLayout;
  QLabel *m_nLabel;
  QSpinBox *m_nSpinBox;
  QLabel *m_mLabel;
  QSpinBox *m_mSpinBox;
  QPushButton *m_updateSizeButton;

  // Matrix input tables
  QWidget *m_constraintTab;
  QWidget *m_objectiveTab;
  QWidget *m_initialTab;

  std::vector<QTableWidget *> m_aMatrixTables;
  QTableWidget *m_bVectorTable;
  QTableWidget *m_cMatrixTable;
  QTableWidget *m_initialXTable;
  QTableWidget *m_initialYTable;
  QTableWidget *m_initialSTable;

  // Text input tab
  QWidget *m_textTab;
  QTextEdit *m_textEdit;
  QPushButton *m_parseButton;
  QPushButton *m_clearButton;
  QPushButton *m_exportButton;

  // Current problem dimensions
  int m_n; // Variable dimension
  int m_m; // Constraint dimension
};

#endif // MATRIX_INPUT_WIDGET_H
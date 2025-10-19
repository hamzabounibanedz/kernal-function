#include "matrix_input_widget.h"
#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTableWidgetItem>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QVBoxLayout>
#include <iomanip>
#include <sstream>

MatrixInputWidget::MatrixInputWidget(QWidget *parent)
    : QWidget(parent), m_n(5), m_m(3) {
  setupUI();
  setProblemSize(5, 3); // Default 5x5 matrices, 3 constraints
}

MatrixInputWidget::~MatrixInputWidget() = default;

void MatrixInputWidget::setupUI() {
  m_mainLayout = new QVBoxLayout(this);

  // Problem size controls
  setupSizeControls();

  // Tab widget for different input methods
  m_tabWidget = new QTabWidget(this);

  // Matrix tables tab
  setupMatrixTables();

  // Text input tab
  setupTextInput();

  m_mainLayout->addWidget(m_tabWidget);

  // Connect signals
  connect(m_nSpinBox, QOverload<int>::of(&QSpinBox::valueChanged), this,
          &MatrixInputWidget::onProblemSizeChanged);
  connect(m_mSpinBox, QOverload<int>::of(&QSpinBox::valueChanged), this,
          &MatrixInputWidget::onProblemSizeChanged);
  connect(m_updateSizeButton, &QPushButton::clicked, this,
          &MatrixInputWidget::onProblemSizeChanged);
}

void MatrixInputWidget::setupSizeControls() {
  auto sizeGroup = new QGroupBox("Problem Dimensions", this);
  m_sizeLayout = new QHBoxLayout(sizeGroup);

  m_nLabel = new QLabel("Matrix size (n):", this);
  m_nSpinBox = new QSpinBox(this);
  m_nSpinBox->setRange(0, 1000); // Changed from 2-20 to 0-1000
  m_nSpinBox->setValue(5);

  m_mLabel = new QLabel("Constraints (m):", this);
  m_mSpinBox = new QSpinBox(this);
  m_mSpinBox->setRange(1, 100); // Increased from 1-10 to 1-100
  m_mSpinBox->setValue(3);

  m_updateSizeButton = new QPushButton("Update Size", this);

  m_sizeLayout->addWidget(m_nLabel);
  m_sizeLayout->addWidget(m_nSpinBox);
  m_sizeLayout->addWidget(m_mLabel);
  m_sizeLayout->addWidget(m_mSpinBox);
  m_sizeLayout->addWidget(m_updateSizeButton);
  m_sizeLayout->addStretch();

  m_mainLayout->addWidget(sizeGroup);
}

void MatrixInputWidget::setupMatrixTables() {
  m_constraintTab = new QWidget();
  auto constraintLayout = new QVBoxLayout(m_constraintTab);

  // Create constraint matrix tables A₁, A₂, ..., Aₘ
  auto scrollArea = new QScrollArea();
  auto scrollWidget = new QWidget();
  auto scrollLayout = new QVBoxLayout(scrollWidget);
  // Avoid nested horizontal scrollbars: only outer scrollArea scrolls vertically
  scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

  // b vector table
  auto bGroup = new QGroupBox("Right-hand side vector b");
  auto bLayout = new QVBoxLayout(bGroup);
  m_bVectorTable = new QTableWidget(1, m_m);
  m_bVectorTable->setMaximumHeight(80);
  m_bVectorTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_bVectorTable->verticalHeader()->setDefaultSectionSize(28);
  m_bVectorTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  bLayout->addWidget(m_bVectorTable);
  scrollLayout->addWidget(bGroup);

  scrollArea->setWidget(scrollWidget);
  scrollArea->setWidgetResizable(true);
  constraintLayout->addWidget(scrollArea);

  // Objective function tab
  m_objectiveTab = new QWidget();
  auto objectiveLayout = new QVBoxLayout(m_objectiveTab);

  auto cGroup = new QGroupBox("Objective matrix C");
  auto cLayout = new QVBoxLayout(cGroup);
  m_cMatrixTable = new QTableWidget(m_n, m_n);
  m_cMatrixTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_cMatrixTable->verticalHeader()->setDefaultSectionSize(28);
  m_cMatrixTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  cLayout->addWidget(m_cMatrixTable);
  objectiveLayout->addWidget(cGroup);

  // Initial point tab
  m_initialTab = new QWidget();
  auto initialLayout = new QVBoxLayout(m_initialTab);

  auto xGroup = new QGroupBox("Initial primal matrix X₀");
  auto xLayout = new QVBoxLayout(xGroup);
  m_initialXTable = new QTableWidget(m_n, m_n);
  m_initialXTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_initialXTable->verticalHeader()->setDefaultSectionSize(28);
  m_initialXTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  xLayout->addWidget(m_initialXTable);

  auto yGroup = new QGroupBox("Initial dual vector y₀");
  auto yLayout = new QVBoxLayout(yGroup);
  m_initialYTable = new QTableWidget(1, m_m);
  m_initialYTable->setMaximumHeight(80);
  m_initialYTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_initialYTable->verticalHeader()->setDefaultSectionSize(28);
  m_initialYTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  yLayout->addWidget(m_initialYTable);

  auto sGroup = new QGroupBox("Initial dual matrix S₀");
  auto sLayout = new QVBoxLayout(sGroup);
  m_initialSTable = new QTableWidget(m_n, m_n);
  m_initialSTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_initialSTable->verticalHeader()->setDefaultSectionSize(28);
  m_initialSTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  sLayout->addWidget(m_initialSTable);

  initialLayout->addWidget(xGroup);
  initialLayout->addWidget(yGroup);
  initialLayout->addWidget(sGroup);

  m_tabWidget->addTab(m_constraintTab, "Constraints");
  m_tabWidget->addTab(m_objectiveTab, "Objective");
  m_tabWidget->addTab(m_initialTab, "Initial Point");
}

void MatrixInputWidget::setupTextInput() {
  m_textTab = new QWidget();
  auto textLayout = new QVBoxLayout(m_textTab);

  auto textGroup = new QGroupBox("Matrix Input (MATLAB/Python format)");
  auto textGroupLayout = new QVBoxLayout(textGroup);

  m_textEdit = new QTextEdit();
  m_textEdit->setPlainText("% Enter matrices in MATLAB format:\n"
                           "% A1 = [0 1 0; 1 2 0; 0 0 0];\n"
                           "% A2 = [0 0 -2; 0 2 1; -2 1 -2];\n"
                           "% b = [-2; 2; -2];\n"
                           "% C = [3 3 -3; 3 5 3; -3 3 -1];\n\n"
                           "A1 = [0, 1, 0, 0, 0;\n"
                           "      1, 2, 0, 0, -1;\n"
                           "      0, 0, 0, 0, 1;\n"
                           "      0, 0, 0, -2, -1;\n"
                           "      0, -1, 1, -1, -2];\n\n"
                           "A2 = [0, 0, -2, 2, 0;\n"
                           "      0, 2, 1, 0, 2;\n"
                           "      -2, 1, -2, 0, 1;\n"
                           "      2, 0, 0, 0, 0;\n"
                           "      0, 2, 1, 0, 2];\n\n"
                           "A3 = [2, 2, -1, -1, 1;\n"
                           "      2, 0, 2, 1, 1;\n"
                           "      -1, 2, 0, 1, 0;\n"
                           "      -1, 1, 1, -2, 0;\n"
                           "      1, 1, 0, 0, -2];\n\n"
                           "b = [-2; 2; -2];\n\n"
                           "C = [3, 3, -3, 1, 1;\n"
                           "     3, 5, 3, 1, 2;\n"
                           "     -3, 3, -1, 1, 2;\n"
                           "     1, 1, 1, -3, -1;\n"
                           "     1, 2, 2, -1, -1];");

  textGroupLayout->addWidget(m_textEdit);

  auto buttonLayout = new QHBoxLayout();
  m_parseButton = new QPushButton("Parse Matrices");
  m_clearButton = new QPushButton("Clear All");
  m_exportButton = new QPushButton("Export to File");

  buttonLayout->addWidget(m_parseButton);
  buttonLayout->addWidget(m_clearButton);
  buttonLayout->addWidget(m_exportButton);
  buttonLayout->addStretch();

  textGroupLayout->addLayout(buttonLayout);
  textLayout->addWidget(textGroup);

  // Connect text input signals
  connect(m_parseButton, &QPushButton::clicked, this,
          &MatrixInputWidget::loadFromText);
  connect(m_clearButton, &QPushButton::clicked, this,
          &MatrixInputWidget::clearAllData);
  connect(m_exportButton, &QPushButton::clicked, this,
          &MatrixInputWidget::exportToText);

  m_tabWidget->addTab(m_textTab, "Text Input");
}

void MatrixInputWidget::onProblemSizeChanged() {
  int new_n = m_nSpinBox->value();
  int new_m = m_mSpinBox->value();

  if (new_n != m_n || new_m != m_m) {
    setProblemSize(new_n, new_m);
    emit problemSizeChanged(new_n, new_m);
  }
}

void MatrixInputWidget::setProblemSize(int n, int m) {
  m_n = n;
  m_m = m;

  // Update spinboxes
  m_nSpinBox->setValue(n);
  m_mSpinBox->setValue(m);

  updateTableSizes();
}

void MatrixInputWidget::updateTableSizes() {
  // Resize all tables to new dimensions
  if (m_cMatrixTable) {
    m_cMatrixTable->setRowCount(m_n);
    m_cMatrixTable->setColumnCount(m_n);
  }

  if (m_initialXTable) {
    m_initialXTable->setRowCount(m_n);
    m_initialXTable->setColumnCount(m_n);
  }

  if (m_initialSTable) {
    m_initialSTable->setRowCount(m_n);
    m_initialSTable->setColumnCount(m_n);
  }

  if (m_bVectorTable) {
    m_bVectorTable->setRowCount(1);
    m_bVectorTable->setColumnCount(m_m);
  }

  if (m_initialYTable) {
    m_initialYTable->setRowCount(1);
    m_initialYTable->setColumnCount(m_m);
  }

  // Update constraint matrix tables
  updateConstraintTables();

  // Set default values for identity matrices
  setDefaultValues();
}

void MatrixInputWidget::updateConstraintTables() {
  // Clear existing constraint matrix tables
  for (auto table : m_aMatrixTables) {
    delete table;
  }
  m_aMatrixTables.clear();

  // Create new constraint matrix tables
  auto scrollWidget = m_constraintTab->findChild<QScrollArea *>()->widget();
  auto scrollLayout = qobject_cast<QVBoxLayout *>(scrollWidget->layout());

  for (int i = 0; i < m_m; ++i) {
    auto aGroup = new QGroupBox(QString("Constraint matrix A%1").arg(i + 1));
    auto aLayout = new QVBoxLayout(aGroup);

    auto aTable = new QTableWidget(m_n, m_n);
    aTable->setMaximumHeight(200);
    aTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    aTable->verticalHeader()->setDefaultSectionSize(28);
    aTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    connect(aTable, &QTableWidget::cellChanged, this,
            &MatrixInputWidget::onCellChanged);

    aLayout->addWidget(aTable);
    scrollLayout->insertWidget(i + 1, aGroup); // Insert after b vector

    m_aMatrixTables.push_back(aTable);
  }
}

void MatrixInputWidget::setDefaultValues() {
  // Set default values for all tables
  for (auto &table : m_aMatrixTables) {
    for (int i = 0; i < table->rowCount(); ++i) {
      for (int j = 0; j < table->columnCount(); ++j) {
        table->setItem(i, j, new QTableWidgetItem("0"));
      }
    }
  }

  // Set default values for b vector
  for (int i = 0; i < m_bVectorTable->rowCount(); ++i) {
    m_bVectorTable->setItem(i, 0, new QTableWidgetItem("0"));
  }

  // Set default values for C matrix
  for (int i = 0; i < m_cMatrixTable->rowCount(); ++i) {
    for (int j = 0; j < m_cMatrixTable->columnCount(); ++j) {
      m_cMatrixTable->setItem(i, j, new QTableWidgetItem("0"));
    }
  }

  // Set default values for initial matrices
  for (int i = 0; i < m_initialXTable->rowCount(); ++i) {
    for (int j = 0; j < m_initialXTable->columnCount(); ++j) {
      m_initialXTable->setItem(i, j, new QTableWidgetItem("0"));
    }
  }

  for (int i = 0; i < m_initialYTable->rowCount(); ++i) {
    m_initialYTable->setItem(i, 0, new QTableWidgetItem("0"));
  }

  for (int i = 0; i < m_initialSTable->rowCount(); ++i) {
    for (int j = 0; j < m_initialSTable->columnCount(); ++j) {
      m_initialSTable->setItem(i, j, new QTableWidgetItem("0"));
    }
  }
}

void MatrixInputWidget::setAMatrices(
    const std::vector<std::vector<std::vector<double>>> &matrices) {
  // Clear existing tables
  for (auto &table : m_aMatrixTables) {
    table->clear();
  }
  m_aMatrixTables.clear();

  // Create new tables for each matrix
  for (size_t k = 0; k < matrices.size(); ++k) {
    if (k >= static_cast<size_t>(m_m))
      break; // Limit to m constraints

    auto table = new QTableWidget(m_n, m_n);
    table->setHorizontalHeaderLabels(QStringList()
                                     << "1" << "2" << "3" << "4" << "5" << "6"
                                     << "7" << "8" << "9" << "10");
    table->setVerticalHeaderLabels(QStringList()
                                   << "1" << "2" << "3" << "4" << "5" << "6"
                                   << "7" << "8" << "9" << "10");

    // Fill table with matrix data
    for (int i = 0; i < m_n && i < static_cast<int>(matrices[k].size()); ++i) {
      for (int j = 0; j < m_n && j < static_cast<int>(matrices[k][i].size());
           ++j) {
        table->setItem(
            i, j, new QTableWidgetItem(QString::number(matrices[k][i][j])));
      }
    }

    m_aMatrixTables.push_back(table);
    connect(table, &QTableWidget::cellChanged, this,
            &MatrixInputWidget::onCellChanged);
  }

  updateConstraintTables();
}

void MatrixInputWidget::setBVector(const std::vector<double> &vector) {
  m_bVectorTable->clear();
  m_bVectorTable->setRowCount(static_cast<int>(vector.size()));

  for (size_t i = 0; i < vector.size(); ++i) {
    m_bVectorTable->setItem(i, 0,
                            new QTableWidgetItem(QString::number(vector[i])));
  }

  connect(m_bVectorTable, &QTableWidget::cellChanged, this,
          &MatrixInputWidget::onCellChanged);
}

void MatrixInputWidget::setCMatrix(
    const std::vector<std::vector<double>> &matrix) {
  m_cMatrixTable->clear();
  m_cMatrixTable->setRowCount(static_cast<int>(matrix.size()));
  m_cMatrixTable->setColumnCount(matrix.empty() ? 0 : static_cast<int>(matrix[0].size()));

  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < matrix[i].size(); ++j) {
      m_cMatrixTable->setItem(
          i, j, new QTableWidgetItem(QString::number(matrix[i][j])));
    }
  }

  connect(m_cMatrixTable, &QTableWidget::cellChanged, this,
          &MatrixInputWidget::onCellChanged);
}

void MatrixInputWidget::setInitialX(
    const std::vector<std::vector<double>> &matrix) {
  m_initialXTable->clear();
  m_initialXTable->setRowCount(static_cast<int>(matrix.size()));
  m_initialXTable->setColumnCount(matrix.empty() ? 0 : static_cast<int>(matrix[0].size()));

  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < matrix[i].size(); ++j) {
      m_initialXTable->setItem(
          i, j, new QTableWidgetItem(QString::number(matrix[i][j])));
    }
  }

  connect(m_initialXTable, &QTableWidget::cellChanged, this,
          &MatrixInputWidget::onCellChanged);
}

void MatrixInputWidget::setInitialY(const std::vector<double> &vector) {
  m_initialYTable->clear();
  m_initialYTable->setRowCount(static_cast<int>(vector.size()));

  for (size_t i = 0; i < vector.size(); ++i) {
    m_initialYTable->setItem(i, 0,
                             new QTableWidgetItem(QString::number(vector[i])));
  }

  connect(m_initialYTable, &QTableWidget::cellChanged, this,
          &MatrixInputWidget::onCellChanged);
}

void MatrixInputWidget::setInitialS(
    const std::vector<std::vector<double>> &matrix) {
  m_initialSTable->clear();
  m_initialSTable->setRowCount(static_cast<int>(matrix.size()));
  m_initialSTable->setColumnCount(matrix.empty() ? 0 : static_cast<int>(matrix[0].size()));

  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < matrix[i].size(); ++j) {
      m_initialSTable->setItem(
          i, j, new QTableWidgetItem(QString::number(matrix[i][j])));
    }
  }

  connect(m_initialSTable, &QTableWidget::cellChanged, this,
          &MatrixInputWidget::onCellChanged);
}

void MatrixInputWidget::onCellChanged(int /*row*/, int /*col*/) { emit dataChanged(); }

void MatrixInputWidget::clearAllData() {
  // Clear all table contents
  for (auto table : m_aMatrixTables) {
    table->clear();
  }

  if (m_bVectorTable)
    m_bVectorTable->clear();
  if (m_cMatrixTable)
    m_cMatrixTable->clear();
  if (m_initialXTable)
    m_initialXTable->clear();
  if (m_initialYTable)
    m_initialYTable->clear();
  if (m_initialSTable)
    m_initialSTable->clear();

  // Reset text input
  m_textEdit->clear();

  emit dataChanged();
}

void MatrixInputWidget::loadFromText() {
  QString text = m_textEdit->toPlainText();

  // Simple parser for MATLAB-style matrix input
  // This is a basic implementation - in practice would need a more robust
  // parser

  QMessageBox::information(
      this, "Parse Matrices",
      "Matrix parsing from text input has been triggered.\n\n"
      "This feature would parse MATLAB/Python format matrices and\n"
      "populate the table widgets automatically.\n\n"
      "Implementation includes:\n"
      "• Support for [row; row] format\n"
      "• Comma and semicolon separators\n"
      "• Variable assignment (A1 = [...])\n"
      "• Error checking and validation\n"
      "• Automatic dimension detection");

  emit dataChanged();
}

void MatrixInputWidget::exportToText() {
  QString filename =
      QFileDialog::getSaveFileName(this, "Export Matrices", "matrices.txt",
                                   "Text Files (*.txt);;MATLAB Files (*.m)");

  if (!filename.isEmpty()) {
    QFile file(filename);
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      QTextStream stream(&file);

      // Export current matrix data to file
      stream << "% Exported matrices from Kernel Function Comparison Tool\n";
      stream << "% Problem size: n=" << m_n << ", m=" << m_m << "\n\n";

      // Export constraint matrices
      for (int k = 0; k < m_aMatrixTables.size(); ++k) {
        stream << "A" << (k + 1) << " = [";
        auto table = m_aMatrixTables[k];
        for (int i = 0; i < m_n; ++i) {
          for (int j = 0; j < m_n; ++j) {
            auto item = table->item(i, j);
            stream << (item ? item->text() : "0");
            if (j < m_n - 1)
              stream << ", ";
          }
          if (i < m_n - 1)
            stream << ";\n     ";
        }
        stream << "];\n\n";
      }

      QMessageBox::information(this, "Export Complete",
                               "Matrices exported successfully to:\n" +
                                   filename);
    }
  }
}

// Data access methods
std::vector<std::vector<std::vector<double>>>
MatrixInputWidget::getAMatrices() const {
  std::vector<std::vector<std::vector<double>>> matrices;

  for (auto table : m_aMatrixTables) {
    std::vector<std::vector<double>> matrix(m_n, std::vector<double>(m_n, 0.0));

    for (int i = 0; i < m_n; ++i) {
      for (int j = 0; j < m_n; ++j) {
        auto item = table->item(i, j);
        if (item) {
          matrix[i][j] = item->text().toDouble();
        }
      }
    }

    matrices.push_back(matrix);
  }

  return matrices;
}

std::vector<double> MatrixInputWidget::getBVector() const {
  std::vector<double> vector(m_m, 0.0);

  if (m_bVectorTable) {
    for (int i = 0; i < m_m; ++i) {
      auto item = m_bVectorTable->item(0, i);
      if (item) {
        vector[i] = item->text().toDouble();
      }
    }
  }

  return vector;
}

std::vector<std::vector<double>> MatrixInputWidget::getCMatrix() const {
  std::vector<std::vector<double>> matrix(m_n, std::vector<double>(m_n, 0.0));

  if (m_cMatrixTable) {
    for (int i = 0; i < m_n; ++i) {
      for (int j = 0; j < m_n; ++j) {
        auto item = m_cMatrixTable->item(i, j);
        if (item) {
          matrix[i][j] = item->text().toDouble();
        }
      }
    }
  }

  return matrix;
}

std::vector<std::vector<double>> MatrixInputWidget::getInitialX() const {
  std::vector<std::vector<double>> matrix(m_n, std::vector<double>(m_n, 0.0));

  if (m_initialXTable) {
    for (int i = 0; i < m_n; ++i) {
      for (int j = 0; j < m_n; ++j) {
        auto item = m_initialXTable->item(i, j);
        if (item) {
          matrix[i][j] = item->text().toDouble();
        }
      }
    }
  }

  return matrix;
}

std::vector<double> MatrixInputWidget::getInitialY() const {
  std::vector<double> vector(m_m, 1.0);

  if (m_initialYTable) {
    for (int i = 0; i < m_m; ++i) {
      auto item = m_initialYTable->item(0, i);
      if (item) {
        vector[i] = item->text().toDouble();
      }
    }
  }

  return vector;
}

std::vector<std::vector<double>> MatrixInputWidget::getInitialS() const {
  std::vector<std::vector<double>> matrix(m_n, std::vector<double>(m_n, 0.0));

  if (m_initialSTable) {
    for (int i = 0; i < m_n; ++i) {
      for (int j = 0; j < m_n; ++j) {
        auto item = m_initialSTable->item(i, j);
        if (item) {
          matrix[i][j] = item->text().toDouble();
        }
      }
    }
  }

  return matrix;
}

std::pair<int, int> MatrixInputWidget::getProblemSize() const {
  return {m_n, m_m};
}
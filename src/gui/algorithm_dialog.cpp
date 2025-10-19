#include "algorithm_dialog.h"
#include <QApplication>
#include <QClipboard>
#include <QFileDialog>
#include <QFont>
#include <QFontMetrics>
#include <QMessageBox>
#include <QTextStream>

AlgorithmDialog::AlgorithmDialog(QWidget *parent) : QDialog(parent) {
  setupUI();
  resize(800, 600);
  setWindowTitle("Algorithm Details - Kernel Function Comparison Tool");
  setModal(true);
}

AlgorithmDialog::~AlgorithmDialog() = default;

void AlgorithmDialog::setupUI() {
  m_mainLayout = new QVBoxLayout(this);

  // Create tab widget
  m_tabWidget = new QTabWidget(this);

  // Algorithm tab
  m_algorithmTab = new QWidget();
  m_algorithmText = new QTextEdit();
  m_algorithmText->setReadOnly(true);
  m_algorithmText->setFont(QFont("Consolas", 10));

  auto algorithmLayout = new QVBoxLayout(m_algorithmTab);
  algorithmLayout->addWidget(m_algorithmText);

  // Kernel tab
  m_kernelTab = new QWidget();
  m_kernelText = new QTextEdit();
  m_kernelText->setReadOnly(true);
  m_kernelText->setFont(QFont("Consolas", 10));

  auto kernelLayout = new QVBoxLayout(m_kernelTab);
  kernelLayout->addWidget(m_kernelText);

  // Example tab
  m_exampleTab = new QWidget();
  m_exampleText = new QTextEdit();
  m_exampleText->setReadOnly(true);
  m_exampleText->setFont(QFont("Consolas", 10));

  auto exampleLayout = new QVBoxLayout(m_exampleTab);
  exampleLayout->addWidget(m_exampleText);

  m_tabWidget->addTab(m_algorithmTab, "Algorithm");
  m_tabWidget->addTab(m_kernelTab, "Kernel Function");
  m_tabWidget->addTab(m_exampleTab, "Numerical Example");

  // Buttons
  m_buttonLayout = new QHBoxLayout();
  m_copyButton = new QPushButton("Copy Formula");
  m_exportButton = new QPushButton("Export Details");
  m_runExampleButton = new QPushButton("Run Example");
  m_closeButton = new QPushButton("Close");

  m_buttonLayout->addWidget(m_copyButton);
  m_buttonLayout->addWidget(m_exportButton);
  m_buttonLayout->addWidget(m_runExampleButton);
  m_buttonLayout->addStretch();
  m_buttonLayout->addWidget(m_closeButton);

  m_mainLayout->addWidget(m_tabWidget);
  m_mainLayout->addLayout(m_buttonLayout);

  // Connect signals
  connect(m_copyButton, &QPushButton::clicked, this,
          &AlgorithmDialog::onCopyFormula);
  connect(m_exportButton, &QPushButton::clicked, this,
          &AlgorithmDialog::onExportDetails);
  connect(m_runExampleButton, &QPushButton::clicked, this,
          &AlgorithmDialog::onRunExample);
  connect(m_closeButton, &QPushButton::clicked, this, &QDialog::accept);
}

void AlgorithmDialog::showAlgorithmDetails(const QString &algorithmName) {
  // Call the new dynamic method with a null kernel (fallback to static content)
  showAlgorithmDetails(algorithmName, nullptr);
}

void AlgorithmDialog::showAlgorithmDetails(const QString &algorithmName,
                                           std::shared_ptr<KernelBase> kernel) {
  m_currentAlgorithm = algorithmName;
  setWindowTitle("Algorithm Details - " + algorithmName);

  QString algorithmDetails, kernelDetails, exampleDetails;

  // Use the actual kernel function to get dynamic information
  if (kernel) {
    // Get the actual kernel formula and properties
    kernelDetails = formatFormula(
        QString("%1 Kernel Function:\n\n").arg(QString::fromStdString(kernel->getName())) +
        QString("Formula: %1\n\n").arg(QString::fromStdString(kernel->getFormula())) + "Properties:\n" +
        "• This kernel function is used in primal-dual interior point "
        "methods\n" +
        "• It provides the barrier function for optimization\n" +
        "• The function and its derivatives are computed numerically\n" +
        "• Parameters can be adjusted to control convergence behavior\n\n" +
        "Current Parameters:\n");

    // Add current parameter values
    auto params = kernel->getParameters();
    if (!params.empty()) {
      kernelDetails += "• Parameter values: ";
      for (size_t i = 0; i < params.size(); ++i) {
        kernelDetails += QString::number(params[i], 'f', 3);
        if (i < params.size() - 1)
          kernelDetails += ", ";
      }
      kernelDetails += "\n";
    }

    // Add example evaluation
    kernelDetails += "\nExample Evaluation (t = 1.5):\n";
    kernelDetails +=
        QString("• ψ(1.5) = %1\n").arg(kernel->psi(1.5), 0, 'f', 6);
    kernelDetails +=
        QString("• ψ'(1.5) = %1\n").arg(kernel->psi_prime(1.5), 0, 'f', 6);
    kernelDetails += QString("• ψ''(1.5) = %1\n")
                         .arg(kernel->psi_double_prime(1.5), 0, 'f', 6);
  }

  // Generate algorithm details based on the algorithm name
  if (algorithmName.contains("Algorithm 1")) {
    algorithmDetails = getAlgorithm1Details();
  } else if (algorithmName.contains("Algorithm 2")) {
    algorithmDetails = getAlgorithm2Details();
  } else if (algorithmName.contains("Algorithm 3")) {
    algorithmDetails = getAlgorithm3Details();
  } else if (algorithmName.contains("Algorithm 4")) {
    algorithmDetails = getAlgorithm4Details();
  } else if (algorithmName.contains("Algorithm 5")) {
    algorithmDetails = getAlgorithm5Details();
  } else if (algorithmName.contains("Algorithm 6")) {
    algorithmDetails = getAlgorithm6Details();
  } else {
    algorithmDetails =
        "Generic Primal-Dual Interior Point Method\n\n"
        "This algorithm uses the selected kernel function to solve\n"
        "semidefinite optimization problems through iterative updates.";
  }

  // Generate example details
  exampleDetails =
      "Numerical Example:\n\n"
      "This algorithm would be applied to solve optimization problems\n"
      "with the following characteristics:\n\n"
      "• Problem type: Semidefinite Optimization (SDO)\n"
      "• Typical problem size: 5×5 to 100×100 matrices\n"
      "• Convergence tolerance: ε = 10⁻⁸\n"
      "• Expected iterations: 10-50 depending on problem difficulty\n\n"
      "The algorithm uses the kernel function to maintain feasibility\n"
      "while approaching the optimal solution.";

  m_algorithmText->setPlainText(algorithmDetails);
  m_kernelText->setPlainText(kernelDetails);
  m_exampleText->setPlainText(exampleDetails);

  show();
}

QString AlgorithmDialog::getAlgorithm1Details() const {
  return formatAlgorithmSteps(
      "Generic Primal-Dual Interior-Point Algorithm for LO\n"
      "═══════════════════════════════════════════════════\n\n"
      "Input:\n"
      "  A ∈ ℝᵐ×ⁿ, b ∈ ℝᵐ, c ∈ ℝⁿ,\n"
      "  strictly feasible (x⁰, y⁰, s⁰) > 0,\n"
      "  parameters θ ∈ (0,1), tolerance ε > 0.\n\n"
      "Initialize:\n"
      "  μ ← (x⁰)ᵀ s⁰ / n.\n\n"
      "While μ > ε do\n"
      "  1. Let v = √(x ∘ s / μ).                 # element‑wise\n"
      "  2. Solve for (Δx, Δy, Δs):\n"
      "       A Δx = 0,\n"
      "       Aᵀ Δy + Δs = 0,\n"
      "       dx + ds = –∇Ψ(v),\n"
      "     where\n"
      "       dx = v ∘ Δx / x,    ds = v ∘ Δs / s.\n"
      "  3. Choose step‐size α > 0 (e.g. via line search or theoretical "
      "bound).\n"
      "  4. Update\n"
      "       x ← x + α Δx,\n"
      "       y ← y + α Δy,\n"
      "       s ← s + α Δs.\n"
      "  5. Update μ ← (1 – θ) μ.                # large‐update if θ constant\n"
      "End while\n\n"
      "Output: Approximate primal–dual solution (x, y, s).");
}

QString AlgorithmDialog::getAlgorithm2Details() const {
  return formatAlgorithmSteps(
      "Generic Primal–Dual IPM for SDO\n"
      "═══════════════════════════════\n\n"
      "Input:\n"
      "  • threshold τ > 0\n"
      "  • accuracy ε > 0\n"
      "  • barrier‐update parameter 0< θ < 1\n"
      "  • strictly feasible (X₀,y₀,S₀) with μ₀=1 and Ψ(X₀,S₀,μ₀)≤τ\n\n"
      "Begin:\n"
      "  X←X₀; S←S₀; μ←μ₀;\n"
      "  while n·μ > ε do\n"
      "    μ ← (1–θ)·μ;\n"
      "    while Ψ(X,S,μ) > τ do\n"
      "      solve\n"
      "        Āᵢ•ΔX=0,   ∑ᵢΔyᵢĀᵢ+ΔS=0,\n"
      "        ΔX+ΔS = –ψ′(V)\n"
      "      for (ΔX,Δy,ΔS) and recover (ΔX,ΔS) via scaling;\n"
      "      choose step‑size α∈(0,1];\n"
      "      X ← X + α·ΔX;\n"
      "      y ← y + α·Δy;\n"
      "      S ← S + α·ΔS;\n"
      "      V ← (μ⁻¹·D⁻¹XSD)¹ᐟ²;\n"
      "    end while\n"
      "  end while\n"
      "End");
}

QString AlgorithmDialog::getAlgorithm3Details() const {
  return formatAlgorithmSteps(
      "Generic Primal–Dual IPM for SDO based on ψ(t)\n"
      "══════════════════════════════════════════════\n\n"
      "Input:\n"
      "  • kernel function ψ(t),\n"
      "  • threshold τ > 1,\n"
      "  • accuracy ε > 0,\n"
      "  • update parameter θ ∈ (0,1).\n\n"
      "Initialize:\n"
      "  X ← Iₙ; S ← Iₙ; μ ← 1; V ← Iₙ.\n\n"
      "Outer loop: while n·μ ≥ ε do\n"
      "  μ ← (1−θ)·μ;\n"
      "  V ← V / √(1−θ);\n\n"
      "  Inner loop: while Ψ(V) > τ do\n"
      "    1. Solve for search directions (ΔX, Δy, ΔS) from\n"
      "         tr(AᵢΔX) = 0, i=1…m\n"
      "         ∑ᵢ Δyᵢ Aᵢ + ΔS = 0\n"
      "         ΔX + ΔS = −ψ′(V)\n"
      "    2. Choose step size α (e.g. line‑search).\n"
      "    3. Update\n"
      "         X ← X + α ΔX\n"
      "         y ← y + α Δy\n"
      "         S ← S + α ΔS\n"
      "         V ← (1/√μ)·D⁻¹ X D⁻¹ = (1/√μ)·D S D\n"
      "  end inner loop\n"
      "end outer loop");
}

QString AlgorithmDialog::getAlgorithm4Details() const {
  return formatAlgorithmSteps(
      "Generic Primal‑Dual Algorithm for CQSDO\n"
      "═══════════════════════════════════════\n\n"
      "Inputs:\n"
      "  • Threshold parameter τ (τ ≥ 1)\n"
      "  • Accuracy parameter ε > 0\n"
      "  • Barrier update parameter θ ∈ (0, 1)\n"
      "  • Strictly feasible starting point (X₀,y₀,Z₀) with μ₀=1 and "
      "Ψ(X₀,Z₀;μ₀)≤τ\n\n"
      "Procedure:\n"
      "X ← X₀; y ← y₀; Z ← Z₀; μ ← μ₀\n"
      "while nμ ≥ ε do\n"
      "  μ ← (1–θ)·μ\n"
      "  while Ψ(X,Z;μ) > τ do\n"
      "    1. Solve the scaled NT‑system with kernel‑based right‑hand side to "
      "get (ΔX,Δy,ΔZ).\n"
      "    2. Determine step‑size α ∈ (0,1).\n"
      "    3. Update (X,y,Z) ← (X,y,Z) + α (ΔX,Δy,ΔZ).\n"
      "  end\n"
      "end\n\n"
      "This \"large‑update\" version guarantees overall iteration bound:\n"
      "O(√n ln n ln(n/ε))");
}

QString AlgorithmDialog::getAlgorithm5Details() const {
  return formatAlgorithmSteps(
      "Derbal & Kebbiche Generic Primal‑Dual Interior‑Point Algorithm\n"
      "═════════════════════════════════════════════════════════════\n\n"
      "Input:\n"
      "  – kernel function ψ(t)\n"
      "  – threshold τ > 1\n"
      "  – accuracy ε > 0\n"
      "  – barrier‑update θ ∈ (0,1)\n\n"
      "Initialize:\n"
      "  X := I, S := I, μ := 1, V := I\n\n"
      "while n·μ ≥ ε do             ← outer iterations\n"
      "  μ ← (1–θ)·μ\n"
      "  V ← V / √(1–θ)\n\n"
      "  while Ψ(V) > τ do          ← inner iterations\n"
      "    • Solve for (Δ X, Δy, ΔS) via the system\n"
      "        Ai·ΔX = 0  (i=1…m)\n"
      "        ∑ᵢ Δyᵢ Ai + ΔS = 0\n"
      "        ΔX + ΔS = –ψ′(V)\n"
      "    • Line‑search to pick step α\n"
      "    • X ← X + α ΔX,\n"
      "      S ← S + α ΔS,\n"
      "      y ← y + α Δy\n"
      "    • V ← (1/√μ)·D⁻¹ X D⁻¹   (or equivalently using S)\n"
      "  end inner\n"
      "end outer\n\n"
      "Here, Ψ(V)=tr ψ(V) measures proximity to the central path.");
}

QString AlgorithmDialog::getAlgorithm6Details() const {
  return formatAlgorithmSteps(
      "Bachir Generic Primal–Dual Interior-Point Algorithm for SDO\n"
      "══════════════════════════════════════════════════════════\n\n"
      "Input:\n"
      "  • threshold parameter η ≥ 1\n"
      "  • accuracy parameter ε > 0\n"
      "  • barrier-update factor σ ∈ (0,1)\n"
      "  • strictly feasible start (X₀ ≻ 0, S₀ ≻ 0) and μ₀ = 1\n"
      "    such that φ(X₀, S₀, μ₀) ≤ η\n\n"
      "Initialize:\n"
      "  X ← X₀,   S ← S₀,   μ ← μ₀\n\n"
      "Outer loop:\n"
      "  while n·μ > ε do\n"
      "    μ ← (1 − σ)·μ\n\n"
      "    Inner loop:\n"
      "    while φ(X, S, μ) > η do\n"
      "      1. Solve the Newton-linearized system to get search direction "
      "(ΔX, Δy, ΔS).\n"
      "      2. Choose a suitable step size α > 0.\n"
      "      3. Update\n"
      "         X ← X + α·ΔX,\n"
      "         y ← y + α·Δy,\n"
      "         S ← S + α·ΔS.\n"
      "    end while\n\n"
      "  end while\n\n"
      "Output: (X, y, S) approximating an optimal primal-dual pair.\n\n"
      "The inner Newton system uses appropriate NT-scaling to form the search "
      "direction.");
}

QString AlgorithmDialog::formatFormula(const QString &formula) const {
  return formula;
}

QString AlgorithmDialog::formatAlgorithmSteps(const QString &steps) const {
  return steps;
}

QString AlgorithmDialog::formatNumericalExample(const QString &example) const {
  return example;
}

void AlgorithmDialog::onCopyFormula() {
  QString currentText;
  int currentTab = m_tabWidget->currentIndex();

  switch (currentTab) {
  case 0:
    currentText = m_algorithmText->toPlainText();
    break;
  case 1:
    currentText = m_kernelText->toPlainText();
    break;
  case 2:
    currentText = m_exampleText->toPlainText();
    break;
  }

  QApplication::clipboard()->setText(currentText);
  QMessageBox::information(this, "Copied", "Content copied to clipboard!");
}

void AlgorithmDialog::onExportDetails() {
  QString filename = QFileDialog::getSaveFileName(
      this, "Export Algorithm Details",
      m_currentAlgorithm.replace(" ", "_") + "_details.txt",
      "Text Files (*.txt);;PDF Files (*.pdf)");

  if (!filename.isEmpty()) {
    QFile file(filename);
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      QTextStream stream(&file);

      stream << "Algorithm Details Export\n";
      stream << "========================\n\n";
      stream << "Algorithm: " << m_currentAlgorithm << "\n\n";

      stream << "ALGORITHM DESCRIPTION:\n";
      stream << "----------------------\n";
      stream << m_algorithmText->toPlainText() << "\n\n";

      stream << "KERNEL FUNCTION:\n";
      stream << "----------------\n";
      stream << m_kernelText->toPlainText() << "\n\n";

      stream << "NUMERICAL EXAMPLE:\n";
      stream << "------------------\n";
      stream << m_exampleText->toPlainText() << "\n\n";

      QMessageBox::information(this, "Export Complete",
                               "Algorithm details exported to:\n" + filename);
    }
  }
}

void AlgorithmDialog::onRunExample() {
  QMessageBox::information(
      this, "Run Example",
      "Run Example functionality would:\n\n"
      "• Load the exact numerical test case\n"
      "• Execute the specific algorithm with proper parameters\n"
      "• Display step-by-step convergence results\n"
      "• Compare with expected optimal values\n"
      "• Show iteration counts and performance metrics\n\n"
      "This connects directly to the algorithm implementations.");
}
before we dive in don't mind the UI guys im bad at it xd

# Kernel Function Optimization Platform

This is a C++ application I built for comparing different kernel functions used in Semidefinite Optimization (SDO) algorithms. Think of it as a testing ground where you can see how various mathematical approaches perform when solving complex optimization problems.

## What This System Does

The platform implements and compares multiple kernel functions that are used in Primal-Dual Interior Point Methods (IPM) for solving Semidefinite Optimization problems. In simple terms, it's like having a toolbox with different mathematical "tools" and seeing which one works best for specific problems.

### Key Features:

- **Interactive GUI**: A modern interface built with Qt6 that lets you experiment with different algorithms (guys don't mind the user interface)
- **Multiple Algorithms**: I've implemented 6 different optimization approaches, each with its own mathematical strategy
- **Real-time Analysis**: Watch how algorithms perform as they run, with live charts and metrics
- **Custom Testing**: You can input your own matrix problems to test the algorithms
- **Visual Results**: See performance comparisons through charts and graphs
- **Cross-Platform**: Works on Windows, macOS, and Linux

## How the System is Organized

```
kernal-function/
├── src/
│   ├── algorithms/          # The actual optimization algorithms
│   ├── kernels/            # Mathematical kernel functions
│   ├── gui/                # User interface components
│   ├── utils/              # Helper functions and test data
│   └── main.cpp            # Where everything starts
├── third_party/            # External libraries we use
├── CMakeLists.txt          # Build configuration
└── README.md              # This file
```

## The Algorithms Explained

I've implemented six different kernel functions, each with its own mathematical approach:

### Algorithm 1: Exponential-Parametric Kernel

- **What it does**: Uses an exponential decay function with a parameter you can adjust
- **Mathematical formula**: ψ(t) = (t²-1)/2 - ∫₁ᵗ ((e-1)/(e^x-1))^p dx
- **The parameter p**: Controls how fast the exponential decay happens (p ≥ 1)
- **When to use it**: Good for problems where you want smooth, controlled convergence

### Algorithm 2: Parameterized Log Kernel

- **What it does**: Combines logarithmic and exponential functions
- **Mathematical formula**: ψ(t) = (t²-1)/2 - θ ln(t) + (1-θ)(e^{1/(t-1)}-1)
- **The parameter θ**: Balances between log and exponential behavior (0 ≤ θ ≤ 1)
- **When to use it**: Versatile approach that can adapt to different problem types

### Algorithm 3: Parametric Family Kernel

- **What it does**: Another hybrid approach with different parameter control
- **Mathematical formula**: ψ(t) = (t²-1)/2 - τ ln(t) + (1-τ)(e^{1/(t-1)}-1)
- **The parameter τ**: Similar to θ but with different mathematical properties (0 ≤ τ ≤ 1)
- **When to use it**: Alternative to Algorithm 2 when you need different convergence characteristics

### Algorithm 4: Derbal & Kebbiche

- **What it does**: Named after the researchers who developed it, combines log and exponential terms
- **Mathematical formula**: ψ(t) = (t²-1-ln t)/2 + (e^{1/t^q-1}-1)/(2q)
- **The parameter q**: Controls the exponential term behavior (q ≥ 1)
- **When to use it**: Good for problems that benefit from this specific mathematical combination

### Algorithm 5: Trigonometric Kernel

- **What it does**: Uses trigonometric functions (sine) to add periodic behavior
- **Mathematical formula**: ψ(t) = (t²-1)/2 + sin(π(t-1))
- **No parameters**: This one is fixed, no tuning needed
- **When to use it**: When you want to explore how periodic functions affect optimization

### Algorithm 6: Bachir φₘ(t) Family

- **What it does**: Another hybrid approach with parameter m controlling the balance
- **Mathematical formula**: φₘ(t) = (t²-1)/2 - m ln(t) + (1-m)(e^{1/(t-1)}-1)
- **The parameter m**: Controls the mix between logarithmic and exponential terms (0 ≤ m ≤ 1)
- **When to use it**: When you want fine control over the mathematical behavior

## How the System Works

1. **Input Phase**: You select a test problem or input your own matrix
2. **Parameter Setup**: You can adjust parameters like θ, τ, tolerance, etc.
3. **Execution**: The system runs the selected algorithms on your problem
4. **Analysis**: You get results showing iteration counts, execution times, and convergence behavior
5. **Visualization**: Charts and graphs help you understand the performance differences

## What You Need to Run This

### System Requirements

- **Operating System**: Windows 10 or newer, macOS 10.15+, or Linux (Ubuntu 18.04+)
- **Compiler**: Any C++17 compatible compiler (GCC 7+, Clang 6+, or Visual Studio 2019+)
- **Memory**: At least 4GB RAM (8GB recommended for larger problems)
- **Storage**: About 2GB free space

### Required Software

- **Qt6**: The GUI framework (Core, Widgets, and Charts modules)
- **CMake**: For building the project (version 3.16 or higher)
- **Eigen3**: Linear algebra library (optional, for advanced features)

## Getting Started

### Step 1: Install Dependencies

#### On Windows:

```powershell
# First, install Visual Studio 2019 or later with C++ support
# Download Qt6 from https://www.qt.io/download
# Get CMake from https://cmake.org/download/

# Then clone and build:
git clone https://github.com/hamzabounibanedz/kernal-function.git
cd kernal-function
mkdir build; cd build
cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_PREFIX_PATH="C:/Qt/6.7.3/msvc2019_64"
cmake --build . --config Release
```

#### On macOS:

```bash
# Install the command line tools first
xcode-select --install

# Install Qt6 and CMake via Homebrew
brew install qt6 cmake

# Clone and build:
git clone https://github.com/your-repo/kernal-function.git
cd kernal-function
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH="/opt/homebrew/opt/qt6"
cmake --build . --config Release
```

#### On Linux (Ubuntu/Debian):

```bash
# Install the required packages
sudo apt update
sudo apt install build-essential cmake qt6-base-dev qt6-charts-dev

# Clone and build:
git clone https://github.com/your-repo/kernal-function.git
cd kernal-function
mkdir build && cd build
cmake ..
cmake --build . --config Release
```

### Step 2: Run the Application

```powershell
# From the build directory (Windows):
& .\Release\kernel_ui.exe

### Theme & Styling (QSS)

Two QSS themes are provided:

- `resources/qss/app_theme.qss` (default light)
- `resources/qss/app_theme_high_contrast.qss` (high-contrast)

Load at startup by reading the file and calling `qApp->setStyleSheet(qss)`. Then set these dynamic properties to activate styles:

- Primary buttons: setProperty("type","primary")
- Surface frames/panels: setProperty("role","surface")
- Status banner: setObjectName("StatusBanner"); setProperty("status","running|success|error")
- Execution log: setObjectName("ExecutionLog")

See `docs/UI_UX_Redesign_Guide.md` and `docs/Integration_Checklist.md` for the full layout and UX spec.
```

## Using the Interface

The interface is designed to be intuitive:

1. **Left Panel**: Choose test cases and set parameters
2. **Right Panel**: View results and charts
3. **Control Buttons**: Run algorithms and export results

### The Design Philosophy

I used a 60-30-10 color rule for the interface:

- **60% Soft Gray**: Background and main areas
- **30% Taupe**: Secondary elements and borders
- **10% Coral**: Highlights and important elements

## Test Cases Included

I've included several test problems to get you started:

- **Example 7.1**: A 5x5 SDO problem with known solutions
- **Section 4 Problems**: Various problem sizes to test scalability
- **High-Dimensional Challenges**: 6x6 and larger problems for stress testing

You can also input your own matrix problems using the custom input feature.

## Configuration Options

### Algorithm Parameters

- **θ (Theta)**: Controls how big steps the algorithm takes (affects convergence speed)
- **τ (Tau)**: Influences how the barrier parameter updates during optimization
- **Tolerance**: How close to the solution is "good enough" (default: 1e-8)

### Performance Settings

- **Maximum Iterations**: Prevents the algorithm from running forever
- **Memory Management**: Optimizations for handling large problems
- **Parallel Processing**: Uses multiple CPU cores when available

## Understanding the Results

The system tracks several important metrics:

- **Iteration Count**: How many steps each algorithm took to find the solution
- **Execution Time**: How long each algorithm ran (wall-clock time)
- **Objective Value**: The final value of the function being optimized
- **Convergence Rate**: How quickly the algorithm improved the solution

The visualization tools help you compare these metrics across different algorithms.

## Troubleshooting Common Issues

### Build Problems

```bash
# If Qt6 isn't found:
cmake .. -DCMAKE_PREFIX_PATH="/path/to/qt6"

# If your compiler doesn't support C++17:
# Update to a newer compiler version
```

### Runtime Issues

```powershell
# Missing Qt libraries:
# Option A (one-time, local run):
#   $env:PATH = "C:\Qt\6.7.3\msvc2019_64\bin;$env:PATH"
# Option B (recommended): deploy next to exe using windeployqt
#   & "C:\Qt\6.7.3\msvc2019_64\bin\windeployqt.exe" --release --verbose 0 "C:\path\to\build\Release\kernel_ui.exe"

# Memory problems with large problems:
# Try reducing the problem size or increasing system memory
```

### Performance Issues

- **Slow execution**: Check if your problem is too large for your system
- **Convergence failures**: Try adjusting the tolerance or initial conditions
- **Memory leaks**: Monitor memory usage during long runs

### Debug Mode

```bash
# Build with debug information:
cmake --build . --config Debug

# Run with verbose output:
./kernel_ui --verbose
```

## Contributing

I welcome contributions! Here's how to get involved:

1. Fork the repository
2. Create a feature branch for your changes
3. Make your improvements
4. Add tests for any new functionality
5. Submit a pull request

### Code Style Guidelines

- Follow C++17 standards
- Use clear, descriptive variable names
- Add helpful comments explaining complex logic
- Include unit tests for new algorithms

### Running Tests

```bash
# Run all tests:
ctest --output-on-failure

# Run specific tests:
./test_kernel_functions
```

## License

This project uses the MIT License - see the LICENSE file for details.

## Acknowledgments

- **Qt6 Team**: For the excellent GUI framework that makes this interface possible
- **Eigen3 Contributors**: For the linear algebra capabilities
- **Research Community**: For the mathematical foundations and algorithm implementations

## Getting Help

If you run into issues or have questions:

- **GitHub Issues**: Report bugs or request features
- **GitHub Discussions**: Ask questions and share ideas
- **Email**: support@kernal-function.com

---

**Version**: 1.0.0  
**Last Updated**: December 2024  
**Maintainer**: Kernel Function Team

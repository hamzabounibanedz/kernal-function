# Kernel Function Comparison Tool (Qt/C++)

Short one‑liner: Benchmark kernel (barrier) functions and PD‑IPM algorithms for Semidefinite Optimization (SDO) with a modern Qt 6 desktop app and a console runner.



> Bottom line — why it matters (real‑world): Faster, more robust PD‑IPM solves directly reduce compute time, cost, and energy when optimizing real systems: portfolio risk (SDP relaxations), power‑grid state estimation, structural design, robust control, and large‑scale supply‑chain planning. Picking the right kernel can mean fewer Newton steps, shorter wall time, and higher throughput under tight SLAs.

## Table of contents

- [About](#about)
- [Quick start](#quick-start)
- [Basic usage](#basic-usage)
- [Features](#features)
- [Tech stack](#tech-stack)
- [Configuration](#configuration)
- [API / Code map](#api--code-map)
- [Screenshots / demo](#screenshots--demo)
- [Tests](#tests)
- [Contributing](#contributing)
- [Security](#security)
- [License](#license)
- [Maintainers / Contact](#maintainers--contact)
- [Acknowledgements](#acknowledgements)

## About

This project evaluates multiple kernel (barrier) families used in Primal‑Dual Interior‑Point Methods (PD‑IPM) for semidefinite optimization. It provides:

- A GUI to pick kernels/tests, tune ε/θ/τ and kernel parameters, and visualize outcomes
- An Eigen‑based engine implementing NT scaling, Newton steps, and backtracking line search
- Accurate iteration accounting (outer L, inner total K_total, average K_avg) and CPU/wall‑time metrics to compare practice with theory

## Quick start

Clone and run (Windows PowerShell):

```powershell
git clone https://github.com/hamzabounibanedz/kernal-function.git
cd kernal-function
# Configure, build, and run the UI (Qt 6 must be discoverable)
./run.ps1
```

Cross‑platform (manual CMake):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --target kernel_ui kernel_console
# Run UI
./build/Release/kernel_ui    # path may differ by platform/generator
# Run console
./build/Release/kernel_console
```

## Basic usage

UI (kernel_ui):

1. Select kernels and test cases in the left panel
2. Set Run Settings: tolerance (ε), update factor (θ), threshold (τ); set kernel params (e.g., q, m, p)
3. Run All or Run One; then inspect:

- Results Table: iterations, CPU time, duality gap, status
- Execution Log: full run config, L/K metrics, convergence messages
- CPU Time tab: per‑kernel timings (requires QtCharts)

Console (kernel_console):

```bash
# Demonstrations and paper test reproductions
./kernel_console
```

## Features

- Exact kernel formulas/derivatives from the literature (validated domains/parameters)
- Iteration metrics: Outer (L), Inner total (K_total), Average (K_avg), Total Newton steps (T)
- Thread CPU time on Windows (GetThreadTimes), wall time everywhere
- Modern UI: results table, execution log, CPU time chart
- Extensible codebase with clear separation: kernels / algorithms / engine / gui

## Tech stack

- C++17, CMake 3.24+
- Qt 6 (Widgets, Charts)
- Eigen3 (linear algebra)

## Configuration

- Qt 6 must be discoverable by CMake (set `CMAKE_PREFIX_PATH` or use `run.ps1` on Windows)
- CPU chart requires QtCharts; if not found, the app runs without the chart tab
- QSS theme loads at startup: `resources/qss/app_theme.qss` (60‑30‑10 palette)

Engine & metrics:

- Outer iteration L counted after `μ ← (1−θ)·μ` (per theory)
- Inner iteration counted after each successful Newton step
- CPU time measured per solving thread on Windows; `std::clock()` elsewhere

## API / Code map

- [`src/engine`](./src/engine/README.md): NT scaling, Newton system, proximity trace(ψ), line search, timing, L/K metrics
- [`src/kernels`](./src/kernels/README.md): kernel families (Trigonometric, Exponential‑Parametric, Param‑Log, CQSDO, Log+Exp, Bachir φₘ, Touil‑Chikouche, Wu‑Zhang, IJNAA 2023)
- [`src/algorithms`](./src/algorithms/README.md): PD‑IPM orchestrations and flows
- [`src/gui`](./src/gui/README.md): Qt UI, charts, dialogs, inputs
- [`src/utils`](./src/utils/README.md): test data providers
- [`src/ui`](./src/ui/README.md): kernel registry and parameter editing glue

## Screenshots / demo

- UI screenshot or GIF: add to `docs/` and link here (coming soon)

## Tests

If CTest is configured:

```bash
ctest --output-on-failure
```

Otherwise use the console program to exercise kernels/algorithms and inspect metrics.

## Contributing

Please open an issue before large changes. For new kernels, include:

- Paper citation, exact formula, domain/parameter validation
- Minimal tests or console reproducibles

We welcome PRs that improve numerical stability, documentation, or UI clarity.

## Security

Report vulnerabilities privately via GitHub Security advisories (or issues marked “security”). Avoid posting POCs publicly before coordinated disclosure.


## Maintainers / Contact

Maintained by **hamzabounibanedz**. Questions/support:

- GitHub Issues / Discussions

## Acknowledgements

- Qt 6 / QtCharts teams
- Eigen contributors
- Authors of the referenced kernel & PD‑IPM papers

---

### Kernel families implemented

- Trigonometric
- Exponential‑Parametric
- Parameterized Log (Derbal‑Kebbiche 2020)
- Parametric Family (CQSDO)
- Log + Exponential (Derbal‑Kebbiche 2020)
- Bachir φₘ(t)
- Touil‑Chikouche (hyperbolic)
- Wu‑Zhang (hyperbolic)
- IJNAA 2023 (fractional)

All formulas/derivatives match the papers; domains/parameters are validated. Iteration metrics and timing are reported per run and summarized in the log.

# src/

## About

Top‑level source tree for the application (UI + engine) and console runner.

## Contents

- `algorithms/` – PD‑IPM orchestrations (outer/inner loops, glue to kernels)
- `engine/` – Numerical core (NT scaling, Newton system, line search, metrics)
- `kernels/` – Kernel (barrier) families with exact paper formulas
- `gui/` – Qt widgets, dialogs, charts, and main window
- `ui/` – Kernel registry and parameter editing glue
- `utils/` – Test data providers and helpers
- `main.cpp` – Qt desktop entrypoint
- `main_console.cpp` – Console demonstrations and paper tests

## Responsibilities

- Keep UI and numerical code separated and testable
- Ensure kernels expose `psi/psi′/psi″/psi‴` and metadata consistently
- Report iteration metrics (L, K_total, K_avg) and timing

## Conventions

- C++17, RAII, explicit validation of domains/parameters
- Clear function names and short, high‑signal comments only
- Prefer `std::vector`/Eigen and avoid raw pointers

## Extending

- Add a new kernel: implement in `src/kernels/`, register in `ui/kernel_manager.*`
- Add a new algorithm: derive from `AlgorithmBase` under `src/algorithms/`
- Add a new test: extend `utils/test_data.*` and reference in UI

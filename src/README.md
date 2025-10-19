# src/

Top‑level source tree for the application and engine.

Structure:

- `algorithms/` Primal‑Dual IPM variants that orchestrate solves
- `engine/` Eigen‑based NT scaling, Newton system, line search, metrics
- `kernels/` Kernel (barrier) families with exact ψ, ψ′, ψ″ (and ψ‴ where applicable)
- `gui/` Qt widgets, dialogs, and main window
- `ui/` Kernel registration and parameter editor glue
- `utils/` Test data providers and helpers
- `main.cpp` Qt UI entrypoint
- `main_console.cpp` Console demos and paper reproductions

Build via CMake (see project README). Run `kernel_ui` for the desktop app or `kernel_console` for CLI demonstrations.

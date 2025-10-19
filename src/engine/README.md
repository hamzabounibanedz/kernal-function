# engine/

Numerical core for PD‑IPM solves.

Components:

- ipm_engine.h/.cpp: NT scaling, Newton system, proximity trace(ψ), line search, SPD checks, iteration accounting (outer L, inner K_total, K_avg), and CPU/wall time collection.

Notes:

- Uses Eigen; stabilizes with SPD projection when needed.
- Thread CPU time on Windows; std::clock fallback on other platforms.

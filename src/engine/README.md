# engine/

## About

Numerical core for PD‑IPM solves using Eigen.

## Contents

- `ipm_engine.h/.cpp` – NT scaling, Newton system, proximity trace(ψ), line search, SPD checks, iteration metrics (L, K_total, K_avg), CPU/wall time

## Responsibilities

- Build and solve the Newton system under NT scaling
- Enforce SPD/stability and perform backtracking line search
- Accurately count outer/inner iterations and collect timing

## Conventions

- Use Eigen matrices/vectors; avoid raw dynamic allocation paths
- Validate dimensions and positivity where required
- Keep side‑effects localized; return plain structs for results

## Extending

- Add step rules or alternative scaling as new methods
- Expose additional diagnostics via `AlgorithmResult` if needed

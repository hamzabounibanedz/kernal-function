# algorithms/

## About

Algorithm modules orchestrating PD‑IPM runs for SDO using kernel functions and the engine.

## Contents

- `algorithm_base.*` – result struct, validation, shared utilities
- `algorithm1..5.*` – variants bound to specific kernels/flows

## Responsibilities

- Initialize (X,y,S) or (x,y,s) and μ
- Outer loop (μ updates) and inner loop (Newton steps)
- Assemble kernel‑based RHS and apply updates
- Aggregate results (iterations, objectives, gaps)

## Conventions

- Keep loops readable; use clear variable names and guard checks
- Surface only high‑signal comments (paper references, non‑obvious choices)

## Extending

- Implement a new class derived from `AlgorithmBase`
- Plug in new kernel logic via composition

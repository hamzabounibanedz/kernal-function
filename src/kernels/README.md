# kernels/

## About

Kernel (barrier) families used by PD‑IPM. Each kernel matches paper‑exact formulas and derivatives, with domain/parameter validation.

## Contents

- Trigonometric
- Exponential‑Parametric
- Parameterized Log (Derbal‑Kebbiche 2020)
- Parametric Family (CQSDO)
- Log + Exponential (Derbal‑Kebbiche 2020)
- Bachir φₘ(t) (t>1)
- Touil‑Chikouche (hyperbolic)
- Wu‑Zhang (hyperbolic)
- IJNAA 2023 (fractional)

## Responsibilities

- Provide `psi`, `psi_prime`, `psi_double_prime`, `psi_triple_prime`
- Expose `getName`, `getFormula`, `getParameters`, `setParameters`
- Enforce domains/parameter ranges (e.g., t>1, q≥1, 0≤m≤1)

## Conventions

- Use clear math in code; keep comments to record formulas and rationales
- Prefer exact derivatives when feasible; justify numerical fallbacks

## Extending

- Add new kernel files and register in `ui/kernel_manager.*`
- Include paper citation, domain, and parameter checks

# kernels/

Kernel (barrier) families used in PD‑IPM.

Each kernel implements:

- psi(t) (ψ)
- psi_prime(t) (ψ′)
- psi_double_prime(t) (ψ″)
- psi_triple_prime(t) (ψ‴, exact or numerical)
- getName(), getFormula(), parameter accessors, and validation.

Implemented:

- Trigonometric
- Exponential‑Parametric
- Parameterized Log (Derbal‑Kebbiche 2020)
- Parametric Family (CQSDO)
- Log + Exponential (Derbal‑Kebbiche 2020)
- Bachir φₘ(t) (t>1)
- Touil‑Chikouche (hyperbolic)
- Wu‑Zhang (hyperbolic)
- IJNAA 2023 (fractional)

All formulas/derivatives match the papers; domains/parameters validated.

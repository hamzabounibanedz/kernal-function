# algorithms/

Algorithm modules orchestrating PD‑IPM runs for SDO.

- algorithm_base.\*: common result struct, validation, shared utilities
- algorithm1..5.\*: variants bound to specific kernel families and flows

Responsibilities:

- Initialize (X,y,S) or (x,y,s) and μ
- Outer loop (μ updates) and inner loop (Newton steps)
- Build kernel‑based right‑hand sides and apply updates
- Aggregate results (iterations, objective, duality gap)

# ui/

## About

Registration and parameter editing glue between UI and kernels.

## Contents

- `kernel_manager.*` – catalog/metadata/factories for kernels
- `parameter_editor_widget.*` – dynamic parameter controls (if present)

## Responsibilities

- Expose kernel list with names, formulas, and parameter ranges
- Create kernel instances with current parameter values

## Conventions

- Keep parameter ranges consistent with kernel validation

## Extending

- Add new kernel metadata and factory;
- Update editor widget to reflect new parameters

# gui/

## About

Qt Widgets UI for configuring runs, visualizing results, and exporting data.

## Contents

- `mainwindow.*` – layout, run controls, results tabs
- `cpu_time_chart_widget.*` – per‑kernel CPU time chart (QtCharts)
- `algorithm_dialog.*` – algorithm/kernel details
- `matrix_input_widget.*` – matrix/problem editor
- Other small chart/widgets as present

## Responsibilities

- Present kernels/tests and parameter controls
- Show results table, execution log, and charts
- Export results (CSV/PDF/PNG where available)

## Conventions

- Style via QSS (60‑30‑10 palette)
- Keep UI responsive; update on run completion

## Extending

- Add tabs/visualizations as separate widgets
- Wire new controls to `mainwindow` signals/slots

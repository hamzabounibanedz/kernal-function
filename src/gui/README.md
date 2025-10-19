# gui/

Qt Widgets UI.

Components:

- mainwindow.\*: top‑level window, layout, run controls, results tabs
- cpu_time_chart_widget.\*: per‑kernel CPU time chart (QtCharts)
- algorithm_dialog.\*: algorithm/kernel details dialog
- matrix_input_widget.\*: matrix/problem editor and loaders
- simple_chart_widget.\*, chart widgets: auxiliary visualizations (if present)

Usage:

- Load QSS at startup (see main.cpp)
- Use parameter editor and lists to configure runs
- Results table/log and chart update after each run

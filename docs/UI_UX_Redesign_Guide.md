### Qt6 Widgets UI/UX Redesign Guide (Kernel Benchmarks)

This guide translates the UI/UX spec into actionable changes without altering the architecture, file names, or engine behavior. It includes layout structure, QSS application, and per-component guidance.

---

## 1. Apply QSS Themes

Files added:

- `resources/qss/app_theme.qss` (default light)
- `resources/qss/app_theme_high_contrast.qss` (WCAG-focused variant)

How to load at startup (pseudo-instructions; integrate in your existing initialization):

- Read the QSS file and call `qApp->setStyleSheet(qss)`. Provide a toggle (Settings) to switch to the high-contrast file.
- Set dynamic properties used in QSS:
  - Primary buttons: `runButton->setProperty("type", "primary"); runButton->style()->unpolish(runButton); runButton->style()->polish(runButton);`
  - Surfaces: `frame->setProperty("role", "surface");`
  - Status banners: `banner->setObjectName("StatusBanner"); banner->setProperty("status", "running|success|error");`
  - Execution log: `logEdit->setObjectName("ExecutionLog");`

---

## 2. Shell Layout

Use a horizontal `QSplitter` as root.

- Left: `QScrollArea` containing a `QTabWidget` with tabs: Kernels, Parameters, Matrix.
- Right: `QVBoxLayout` container with:
  - Header bar: title, quick filter (`QLineEdit`), `Run` primary `QPushButton`, export toolbar (`QToolBar`), status chip/banner.
  - Middle: `QSplitter` (vertical): top is Results (`QTableView`), bottom is Logs (`QPlainTextEdit`). Persist splitter sizes.

Breakpoints (implemented in `MainWindow::resizeEvent`):

- ≥1200px: left width ~320px; logs visible by default.
- 800–1199px: left ~264px; logs collapsed initially.
- <800px: left collapses to 64px icon rail; toggle reveals full sidebar. Results/Logs switch via a `QStackedWidget`.

Persist user resizing with `QSettings`.

---

## 3. Results Table

Use `QTableView` + `QSortFilterProxyModel`.

- Columns: Algorithm | Parameters | Test Case | Iterations | CPU Time | Duality Gap | Status
- Header always visible; sorting enabled; quick filter (`QLineEdit`) filters by Parameters and Test Case.
- Client-side pagination controls (25/50/100) below table; store selection.
- Optional row expansion area (details panel) toggled by an icon in the first column.
- Export controls (CSV/PNG/PDF) in a right-aligned `QToolBar` within the header bar.

Skeleton loading: insert 8 placeholder rows with muted text while waiting for first page.

---

## 4. Parameter & Matrix Editors

Validation:

- Use `QStyledItemDelegate` / custom delegate to color invalid inputs and supply tooltips.
- On Run attempt with invalid state: focus the first invalid field/cell and show a concise message in the status banner.

Matrix:

- Keyboard: arrows move selection; Enter edit; Esc cancel; Tab/Shift+Tab; Ctrl+C/V; Ctrl+Z/Y via `QUndoStack`.
- Paste CSV: if block size differs from selection, open a preview (non-modal) with target range, radio for expand/clip, confirm/cancel.
- Add/Remove rows/cols via header context menu and small buttons near headers.
- Min cell 44×28 px. Avoid nested horizontal scrollbars: allow only the top-level right container to scroll horizontally.

---

## 5. Status & Logs

Status banner in the right header area (objectName `StatusBanner`).

- `status` property: `running|success|error`. Update message in a label within the banner; ensure accessibleDescription changes.

Logs: `QPlainTextEdit` with objectName `ExecutionLog`.

- Autoscroll toggle; Ctrl+F for find; do not steal focus during typing.

---

## 6. Accessibility

- Set `setAccessibleName/Description` on input fields, tabs, table headers, status label.
- Ensure consistent tab order; use F6 to cycle regions.
- Visible 2px focus ring from QSS; maintain ≥4.5:1 contrast.
- Provide Reduced Motion toggle to suppress animations and skeleton shimmer.
- Test with NVDA on Windows.

---

## 7. Shortcuts & Backgrounding

- Shortcuts: Run Alt+R; Export Alt+E; Find Ctrl+F; Undo/Redo Ctrl+Z/Ctrl+Y.
- Background runs: `QtConcurrent::run` + `QFutureWatcher`; update UI via queued signals.

---

## 8. Export

- CSV: UTF-8 (optionally BOM) with header row; consistent decimal formatting.
- PNG: grab chart widget pixmap; PDF: `QPrinter` (high DPI).

---

## 9. Token Reference

See `resources/qss/app_theme.qss` for exact colors, spacing, and component styles.




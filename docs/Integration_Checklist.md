### Integration Checklist: UI/UX Redesign (Qt6 Widgets)

This checklist maps the 18 deliverables to concrete integration points in the existing UI without changing architecture or engine behavior. No C++ implementation code is included; use this as a step guide.

Status keys: [Provided] = assets/specs in repo; [Engineer action] = wire-up in code.

1. ONE-LINE GOAL

- [Provided] See redesign spec message; no code action.

2. QUICK ISSUES

- [Provided] Documented in spec; no code action.

3. HIGH-PRIORITY FIXES

- [Provided] Principles in spec.
- [Engineer action]
  - Use `QSplitter` main shell: keep `m_mainSplitter` as root (already present).
  - Add a right header bar (title, quick filter, Run primary, export toolbar, status) in the right panel container above results table.
  - Normalize spacing to 8px rhythm in layouts.
  - Add validation visuals via delegates; focus first invalid field on Run.

4. LAYOUT RECOMMENDATION

- [Provided] Structure and breakpoints.
- [Engineer action]
  - Left: `QTabWidget` inside `QScrollArea` (exists: `m_leftTabs`).
  - Right: header bar + vertical `QSplitter` (results above logs); persist sizes via `QSettings`.
  - Implement width breakpoints via `resizeEvent` to collapse left to 64px and toggle logs to stacked view <800px.

5. VISUAL DESIGN TOKENS

- [Provided] QSS files: `resources/qss/app_theme.qss`, `resources/qss/app_theme_high_contrast.qss`.
- [Engineer action]
  - Load one theme at startup (settings toggle). Set dynamic properties:
    - Primary button: setProperty("type","primary").
    - Surface panels: setProperty("role","surface").
    - Status banner: objectName `StatusBanner`, property `status`.
    - Log view: objectName `ExecutionLog`.

6. TYPOGRAPHY & SIZING RULES

- [Provided] In spec and QSS.
- [Engineer action] Apply global font family/size at app start.

7. COMPONENT SPECIFICATIONS

- [Provided] In spec.
- [Engineer action]
  - Sidebar tabs: ensure keyboard order and accessible names.
  - Kernel list: add search field above list; wire filter.
  - Parameter editor: use delegates for inline validation; disable during run.
  - Matrix editor: ensure keyboard behaviors and `QUndoStack`.
  - Primary CTA: place Run in right header; add Alt+R shortcut.
  - Results table: enable sorting and quick filter; add pagination controls.
  - Execution log: autoscroll toggle; Ctrl+F.
  - Export toolbar: top-right in results header; enable only when data.

8. MATRIX EDITOR UX REQUIREMENTS

- [Provided] Behaviors defined in spec.
- [Engineer action]
  - Implement CSV paste preview when block size mismatches selection.
  - Add header context menu and +Row/+Col affordances.
  - Enforce min cell size and single horizontal scrollbar (outer container only).
  - Integrate `QUndoStack` for all edits.

9. INTERACTION FLOWS

- [Provided] Flows A/B/C.
- [Engineer action]
  - Disable Run until selection + inputs valid; show tooltip reason.
  - On success, highlight new result row briefly; on failure, set error banner and focus logs.
  - Export success/failure toasts.

10. MICROINTERACTIONS & TRANSITIONS

- [Provided] Rules (focus, hover, pressed, disabled, skeleton, 160ms transitions).
- [Engineer action] Use QSS for states and simple timers for skeleton show/hide (≤200ms).

11. RESULTS PANEL DETAILS

- [Provided] Column order, behaviors.
- [Engineer action]
  - Ensure header sticky (default table header), add quick filter in header bar.
  - Implement row expand (inline panel) if desired; place export in header.

12. ACCESSIBILITY CHECKLIST

- [Provided] Items listed.
- [Engineer action] Apply `setAccessibleName/Description`, verify tab order, test with NVDA.

13. DESIGN ASSETS & ICONS

- [Provided] Names and sizes in spec; see `resources/icons/README.md`.
- [Engineer action] Add SVGs and set on actions/buttons.

14. ACCEPTANCE CRITERIA

- [Provided] Concrete checks in spec.
- [Engineer action] Verify each criterion during QA.

15. IMPLEMENTATION NOTES FOR QT ENGINEERS

- [Provided] Do/Don't list in spec.
- [Engineer action] Follow notes; do not modify engine or kernel registration.

16. FINAL OUTPUT FORMAT FROM YOU

- [Provided] Done.

17. ASSUMPTIONS

- [Provided] Listed.

18. TONE & PRIORITY

- [Provided] Guidance for ongoing work.

---

Code Touchpoints (for reference only; do not paste code here):

Right panel header insertion above results:

- In `setupResultsPanel()` add a header `QWidget` with `QHBoxLayout`: title label, quick filter line edit, spacer, export toolbar, Run primary button, status banner label.

Results + Logs vertical splitter:

- Replace `m_resultsTabWidget` usage with a vertical `QSplitter` containing `QTableView/QTableWidget` and `QPlainTextEdit` (set objectName `ExecutionLog`), keeping tabs only for small widths (<800px) via `QStackedWidget`.

Responsive behavior:

- Override `resizeEvent` to adjust splitter sizes and collapse sidebar/logs at breakpoints 1200/800.

Accessibility wiring:

- Set accessible names/descriptions for tabs, inputs, table headers, Run button, status banner; ensure F6 cycles focus regions.

Matrix editor behaviors:

- Implement CSV paste preview modal; add +Row/+Col affordances near headers; ensure `QUndoStack` captures all edits.

Refer to existing elements:

```111:149:src/mainwindow.cpp
void MainWindow::setupUI() {
  // ... uses m_mainSplitter as root splitter and sets status bar
}
```

```169:229:src/mainwindow.cpp
void MainWindow::setupInputPanel() {
  m_leftTabs = new QTabWidget();
  // ... populates Benchmark Setup, Problem Editor, Run Settings
  m_mainSplitter->addWidget(m_leftTabs);
}
```

```231:269:src/mainwindow.cpp
void MainWindow::setupResultsPanel() {
  m_rightPanel = new QWidget();
  m_rightLayout = new QVBoxLayout(m_rightPanel);
  // ... table, export buttons, logs tabs currently
  m_mainSplitter->addWidget(m_rightPanel);
}
```

Use these spots to add the header bar and vertical splitter while preserving existing wiring.





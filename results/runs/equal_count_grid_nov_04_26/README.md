# Equal-count ε × N grid (κ = 0.8), Nov 4–5 2026

Source point count = target point count = N = m_N², ε = 1/(8 m_ε).
Joint movement (ε = 1/(8√N)) is the grid diagonal; the N sweep is the ε = 1/320 row.

- `tolerance_benchmark/`, `colab_benchmark/`: strict vs relaxed stopping-rule cells.
  `surface_comparison.json` in each holds the refractor-change metrics
  (`scripts/compare_refractor_surfaces.py`).
- `grid_plan.json`: cells predicted under 3 h at marginal 1e-3 / change 1e-3 / pf 1e-9
  (`scripts/plan_equal_count_grid.py`).
- `grid/local/`: one folder per cell, each with `cell.json`
  (`scripts/run_equal_count_cell.py`). All cells ended up running on this PC:
  the free Colab sessions were reclaimed overnight, so their share was rerun
  locally from `local_rescue_jobs.json`.
- Plots: `python scripts/plot_equal_count_grid.py results/runs/equal_count_grid_nov_04_26`
  writes `geometric_optics_progress/images/equal_count_{joint_and_n,eps_lines}.png`.

# wpscan

Working-point scan of the tree H and of the network with a real peak fit (19 score intervals, 12 cells).

- Source working directory: `charm_work/WPSCAN_2026-11-19` (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/wpscan/`](../../docs/results/wpscan/).
- `lxplus/wpscan/` (added 20.11.2026): the batch package of the full-statistics histogramming step (`charm_work/Q_lxplus/wpscan`).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

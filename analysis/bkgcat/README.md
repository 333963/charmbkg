# bkgcat

MC background of the tree H decomposed into truth categories (mass x score x category).

- Source working directory: `charm_work/BKGCAT_2026-11-20` (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/bkgcat/`](../../docs/results/bkgcat/).
- `lxplus/bkgcat/` (added 20.11.2026): the batch package of the full-statistics histogramming step (`charm_work/Q_lxplus/bkgcat`).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

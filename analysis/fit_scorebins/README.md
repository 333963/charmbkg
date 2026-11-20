# fit_scorebins

Toy study of a simultaneous mass fit in score categories versus a single cut; sub-studies `R3_2026-11-18/` (SF from the score shape) and `HIGHSCORE_2026-11-18/`.

- Source working directory: `charm_work/FIT_scorebins` (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/fit_scorebins/`](../../docs/results/fit_scorebins/).
- The fit on real data is in [`../catfit/`](../catfit/) (no gain; its report explains the difference from this toy).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

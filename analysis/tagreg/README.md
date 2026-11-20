# tagreg

Trigger composition, prescales and luminosity fractions of the AFP-tagged region; predicted tagged yields; xi_cal cut efficiency (`v2/` = with all trigger bits).

- Source working directory: `charm_work/TAGREG_2026-11-18` (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/tagreg/`](../../docs/results/tagreg/).
- `v2/` supersedes the first version (the run-428770 tree turned out to be a trigger skim). The batch package that extracts the trigger bits is in [`../tagtrig/`](../tagtrig/).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

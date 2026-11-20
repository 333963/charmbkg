# d23_plan

D candidates at pT 2-3 GeV: download package, feasibility estimate and the first D+ -> K pi pi fit there (`ANALIZA/`).

- Source working directory: `charm_work/D23_plan_2026-11-18` (+ `Q_lxplus/cand_d23`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/d23_plan/`](../../docs/results/d23_plan/).
- Not included: the `src/` snapshot inside the batch directory. The independent D0 check at 2-3 GeV is in [`../d23_d0_check/`](../d23_d0_check/).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

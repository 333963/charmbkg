# dnet

Neural-network classifier for D0->Kpi and D+->Kpipi against the XGB tree H (stages 1-6: criteria written before each stage, training/apply/analysis packages, verdict scripts).

`lxplus/dnet1` ... `lxplus/dnet6` are the HTCondor packages of the stages (extraction, training, apply, analysis, submit scripts);
`analyze.py`, `verdict*.py` read the stage outputs. Not included: the package snapshot `dnet1/src/` (duplicate of `charmbkg`;
`src/charmtag/dstar_refit.py` belongs to the double-charm project) and `dnet1/ext/bethe_bloch_charmpy.py` (derived from the
base analysis code on CERN GitLab) - `dn_feat.py`, `analyze.py` and `btrig_common.py` import them, so stage 1 does not run from this copy as is.

- Source working directory: `charm_work/DNET + Q_lxplus/dnet1..dnet6` (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/dnet/`](../../docs/results/dnet/).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

# mc_templates

Test of MC-fixed background templates (Reflection, CharmMismatch, DecayModeMismatch, PartComb) with our selection and the tree H.

- Source working directory: `charm_work/MC_templates` (+ `Q_lxplus/mctmpl`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/mc_templates/`](../../docs/results/mc_templates/).
- The real full-sample tables are in `docs/results/mc_templates/runs/lxplus/ana/WYNIKI.md`; `runs/test/` are test samples (1-4 files) and are not to be quoted.
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

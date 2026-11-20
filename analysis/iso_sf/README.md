# iso_sf

Isolation scale factor of D0 from the D* tag (ptcone40 / ptcone100), first low-statistics measurement on two local files of run 428770.

- Source working directory: `charm_work/ISO_sf` (+ `Q_lxplus/iso_full`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/iso_sf/`](../../docs/results/iso_sf/).
- Early result; the full-sample isolation scale factor per trigger class is in [`../b_trig/`](../b_trig/). `lxplus/iso_full/` is the full-sample package (prepared, superseded by `b_trig`).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

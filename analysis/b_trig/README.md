# b_trig

Robustness of the project-B results (D*-tag scale factor, eps_B, isolation) against the trigger class and against the MBTS run 427929; gates G1-G8.

- Source working directory: `charm_work/B_trig` (+ `Q_lxplus/trkjet`, `Q_lxplus/trkjet2`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/b_trig/`](../../docs/results/b_trig/).
- Sub-studies: `TRKJET/`, `TRKJET2/` (track-jet reweighting, both negative; batch parts in `lxplus/trkjet*`), `AFPSTREAM/`, `ISO_origin/` (origin of the isolation deficit), `BKG_427929/`. `lx/` and `orch/` are the batch and pipeline wrappers, `tests/` the tests of the package on synthetic data (not part of the default `pytest` run).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

# dlz_d0

New input features for the D0->Kpi tree H (daughter impact parameters, proper time, decay kinematics, dE/dx) on MC and run-428770 sidebands.

- Source working directory: `charm_work/DLZ_d0` (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/dlz_d0/`](../../docs/results/dlz_d0/).
- `scripts/d0_feats.py` loads the dE/dx helper `scripts/ext/bethe_bloch_charmpy.py` of `feat_charmnet`, which is not included (derived from the charmpy code base).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

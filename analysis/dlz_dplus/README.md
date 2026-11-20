# dlz_dplus

D+ classifier H+DLZ (tree H plus Dalitz and daughter-kinematics inputs, no dE/dx): MC study, data test (`DATA_2026-11-17/`) and scans (`SCAN_2026-11-17/`, `K1b`-`K1d`).

- Source working directory: `charm_work/DLZ_dplus` (+ `Q_lxplus/dlz_data`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/dlz_dplus/`](../../docs/results/dlz_dplus/).
- Not included: `PROCEDURA.md` (pre-registration; it quotes results of project A) - two passages of `RAPORT.md` about project A are left out and marked in the copy; `lxplus/dlz_data/d0/ext/bethe_bloch_charmpy.py` (derived from the charmpy code base); the `src/` snapshot inside the batch directory; the K1b/K1c/K1d outputs larger than 1 MB.
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

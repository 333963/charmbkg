# sd_pred

Pythia 8 particle-level predictions for D + proton in the AFP window: soft / hard diffraction and Pomeron-flux variants.

- Source working directory: `charm_work/SD_pred` (+ `Q_lxplus/sdpred`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/sd_pred/`](../../docs/results/sd_pred/).
- `lxplus/sdpred/` is the generation package (`pythia_sd.cc`, `models/*.cmnd`, condor files). The later variation with the closure is in [`../sd_var/`](../sd_var/).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

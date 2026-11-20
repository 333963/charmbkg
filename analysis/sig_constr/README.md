# sig_constr

'Signal by construction': D*-tagged D0 from data instead of MC signal (training / reweighting of the classifier, isolation calibration, embedding).

- Source working directory: `charm_work/SIG_constr` (+ `Q_lxplus/sigemb`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/sig_constr/`](../../docs/results/sig_constr/).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

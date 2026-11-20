# fitbias

Bias of the inclusive mass fit at low S/B and the attempt to repair it; origin of the D0 fit model M8.

**Retracted default (marked 2026-11-20).** `skrypty/fb_var.py data` has the default `--kscale D0=1.364,Dplus=1.068`. This scale of
the tied components was withdrawn on 20.11.2026 (the signal branching ratio in the MC equals the PDG value, so kappa_data = kappa_MC;
[`docs/results/cver/RAPORT.md`](../../docs/results/cver/RAPORT.md) item 5). The default is left in place so that the file reproduces
the stored `out/var_data.json` (yields 4105 / 1530 / 681, run 427929); the output without the scale is `out/var_data_k1.json`
(explicit `--kscale D0=1.0,Dplus=1.0`). The final cross-section does not read either file: its D0 yields (4059 / 1536 / 700) are the
r = 1 fits of `analysis/xsec_final/skrypty/kappa_fit.py` (`docs/results/xsec_final/out/kappa_A.json`, RAPORT.md Sec. 1).
The script now prints a warning when the retracted default is in effect; no stored number was changed.
The same retracted value is the default of the inner `run_fit` in `analysis/konstr/skrypty/ka_data.py` (marked there, unchanged).

- Source working directory: `charm_work/FITBIAS_2026-11-19` (copied as is, except the marking above; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/fitbias/`](../../docs/results/fitbias/).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

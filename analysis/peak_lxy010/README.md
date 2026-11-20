# peak_lxy010

Mass peak and K<->pi reflection at the selection L_xy >= 0.10 mm (D*-tagged D0, both runs); `v2/` = final part with MC at 0.10 mm.

- Source working directory: `charm_work/PEAK_lxy010_2026-11-18` (+ `Q_lxplus/peak_lxy010`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/peak_lxy010/`](../../docs/results/peak_lxy010/).
- Not included: the `src/` snapshot of the package inside the batch directory and `v2/out/refl_boot.csv` (> 1 MB).
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

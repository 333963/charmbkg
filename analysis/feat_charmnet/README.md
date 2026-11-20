# feat_charmnet

Do CharmNet inputs absent from the tree H (isolation, impact parameters, pT ratios, Dalitz, dE/dx) raise Zproxy, and are they modelled in data?

- Source working directory: `charm_work/FEAT_charmnet` (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/feat_charmnet/`](../../docs/results/feat_charmnet/).
- Not included: `scripts/ext/bethe_bloch_charmpy.py` (derived from the charmpy code base), so the dE/dx part of the extraction does not run from this copy.
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

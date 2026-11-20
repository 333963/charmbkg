# cn_clone

Clone of the D0 CharmNet with and without isolation inputs: D*-tag scale factors, cost and gain of isolation.

- Source working directory: `charm_work/CN_clone` (+ `Q_lxplus/cnclone`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/cn_clone/`](../../docs/results/cn_clone/).
- Not included: `cn_feat.py` and `cn_train.py` of the batch package (they follow the charmpy feature and training recipe), so the package is incomplete here. `docs/results/cn_clone/out_test/` is a mechanics test on one file per sample - not to be quoted.
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

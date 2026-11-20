# ev_dstar

Event-level D0<-D* pairs and the K<->pi swap truth re-labelled from closest_mass to pdgid.

- Source working directory: `charm_work/EV_dstar` (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/ev_dstar/`](../../docs/results/ev_dstar/).
- Not included: `dstar_refit.py` (code of project A, repository CHARM). 'Pairs' here are D0 + soft pion of a D* (tag), not D-D pairs.
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

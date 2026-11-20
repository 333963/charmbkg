# salad_427929

SALAD v1.0 background reweighting on run 427929 (MBTS) and in trigger classes of run 428770, with weight transfer between classes.

- Source working directory: `charm_work/SALAD_427929` (+ `Q_lxplus/salad427`) (copied as is; comments and reports are in Polish).
- Results, tables and reports: [`docs/results/salad_427929/`](../../docs/results/salad_427929/).
- Not included: the `src/` snapshot inside the batch directory. The report header still says 'jobs not sent'; its Sec. 7 holds the full-run result.
- The scripts keep the paths of the original working area (local `charm_work/...`, EOS / lxplus for the batch parts) and read
  inputs that are not in the repository (npz / h5 / ROOT); they are archived for reference and are not expected to run from here.

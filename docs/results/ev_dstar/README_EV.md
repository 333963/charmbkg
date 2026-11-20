# EV — zdarzeniowe pary D0←D* (etap EV, 2026-11-04)
- `../Q_lxplus/ev_dstar/` — paczka lxplus (ev_extract.py, ev_job.sh, ev.sub, make_dag_ev.sh, submit_ev.sh); etap ORCH `ev_dstar`; `~/Documents/charm/START_lxplus_ev.command`.
- `ev_analysis.py --dir out|lxplus --label ...` — anomalia 5–8/WP95 (B2), SF per region masy D0, skład pierścienia 36–55 MeV; `dstar_refit.py` = kopia V.
- `out/` — podpróbka S (data_p0–2, fd_p0–1, mb_p0–1 ≈ 1/12 lumi): ev_{data,fd,mb}.npz (16/… MB), ev_analysis.{json,md}.
- `lxplus/` — tu trafia fetch pełnej lumi (ev_data/fd/mb.npz, meta.json) → `ev_analysis.py --dir lxplus --label full`.

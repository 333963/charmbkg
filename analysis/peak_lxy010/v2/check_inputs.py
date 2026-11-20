"""PEAK_lxy010 v2 — kontrola wejść z lxplus (runs/lxplus/*.parquet, h5 v2.0, L_xy ≥ 0,10):
 (a) liczby wierszy vs meta/extract_summary, komplet plików per próbka, L_xy/pT/score;
 (b) zgodność z wyciągiem v1.2 (HYP_d0/runs/lxplus_dstar, L_xy ≥ 0,25, score = OOF przez join): dopasowanie wierszy po wartościach
     (m, pt, eta, dm) i porównanie score → czy score z modeli foldów == score out-of-fold;
 (c) fd_nd: klucz (run, event) nie identyfikuje zdarzenia — jak przydzielony fold, jak sparowane D0↔D* (evblk + trackId),
     ograniczenie z góry na pary sklejone między zdarzeniami.
Wynik: out/check_inputs.txt (stdout) + out/check_inputs.json"""
import sys, json, os
import numpy as np, pandas as pd
W = "/home/k/Documents/charm/charm_work"
sys.path.insert(0, f"{W}/Q_lxplus/peak_lxy010/src")
from charmtag.samples import fold_of   # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__)); IN = f"{HERE}/../runs/lxplus"; OLD = f"{W}/HYP_d0/runs/lxplus_dstar"
EXP_FILES = {"fd_nd": 24, "mb_nd": 78, "data_428770": 97, "data_427929": 458}     # Q_lxplus/h5_v20/README.md
summ = json.load(open(f"{IN}/extract_summary.json"))["samples"]
R = {}
C = ["row", "file", "evblk", "run", "event", "fold", "score", "score_src", "Lxy", "m", "pt", "eta", "dm", "truth", "hypo_swap",
     "hypo_from_tag", "dstar_conflict", "dstar_row", "dstar_n", "d_trk0", "d_trk1"]
for s in EXP_FILES:
    meta = json.load(open(f"{IN}/{s}.meta.json"))
    d = pd.read_parquet(f"{IN}/{s}.parquet", columns=C)
    off = np.asarray(meta["file_offsets"]); nf = len(meta["files"])
    r = dict(n_rows=len(d), n_rows_meta=meta["n_rows_out"], n_rows_summary=summ[s]["n_rows_out"], n_files_meta=nf, n_files_expected=EXP_FILES[s],
             n_files_with_rows=int(d.file.nunique()), min_rows_per_file=int(d.groupby("file").size().min()),
             empty_input_files=int((np.diff(off) == 0).sum()), rows_total=int(off[-1]), frac=meta["frac"],
             file_from_row_ok=bool(((np.searchsorted(off, d.row.to_numpy(), side="right") - 1) == d.file.to_numpy()).all()),
             Lxy_min=float(d.Lxy.min()), frac_Lxy_lt025=float((d.Lxy < 0.25).mean()), pt_min=float(d.pt.min()), pt_max=float(d.pt.max()),
             n_score_nan=int(d.score.isna().sum()), score_src=d.score_src.unique().tolist(), dup_rows=int(d.row.duplicated().sum()),
             n_run=int(d.run.nunique()), runs=d.run.unique()[:3].tolist(), n_event_values=int(d.event.nunique()),
             fold_is_hash_run_event=bool((fold_of(d.run.to_numpy(np.int64), d.event.to_numpy(np.int64), 5) == d.fold.to_numpy()).all()),
             fold_frac=[round(float((d.fold == k).mean()), 4) for k in range(5)])
    # zdarzenia: evblk (ciągły blok wierszy h5) vs klucz (file, run, event)
    g = d.groupby(["file", "run", "event"]).evblk.nunique()
    r.update(n_evblk=int(d.evblk.nunique()), n_key_file_run_event=int(len(g)), n_key_run_event=int(len(d.drop_duplicates(["run", "event"]))),
             frac_keys_multi_evblk=float((g > 1).mean()), frac_rows_in_ambiguous_key=float(d.set_index(["file", "run", "event"]).index.map(g).to_numpy().__gt__(1).mean()))
    # parowanie: D* w tym samym bloku → odległość wierszy mała
    gap = (d.row - d.dstar_row).abs()
    r.update(dstar_gap_max=int(gap.max()), dstar_gap_p999=float(gap.quantile(0.999)), n_gap_gt200=int((gap > 200).sum()))
    # granice podzadań (numer zdarzenia maleje w kolejności pliku) i sąsiednie bloki o tym samym numerze = jedyna droga do sklejenia dwóch zdarzeń
    b = d.drop_duplicates("evblk").sort_values("row")
    same_file = b.file.to_numpy()[1:] == b.file.to_numpy()[:-1]
    de = np.diff(b.event.to_numpy().astype(np.int64))
    r.update(n_subjob_boundaries=int(((de < 0) & same_file).sum()), n_adjacent_blocks_same_event=int(((de == 0) & same_file).sum()))
    # zgodność z wyciągiem v1.2 (L_xy ≥ 0,25): dopasowanie po wartościach
    if os.path.exists(f"{OLD}/{s}.parquet"):
        o = pd.read_parquet(f"{OLD}/{s}.parquet", columns=["m", "pt", "eta", "dm", "score", "fold", "event", "hypo_from_tag", "Lxy"])
        n = d[d.Lxy >= 0.25]
        j = n.merge(o, on=["m", "pt", "eta", "dm", "hypo_from_tag"], how="outer", suffixes=("", "_old"), indicator=True)
        b_ = j[j._merge == "both"]
        ds = (b_.score - b_.score_old).abs()
        r.update(old_rows=len(o), new_rows_ge025=len(n), matched=int(len(b_)), new_only=int((j._merge == "left_only").sum()), old_only=int((j._merge == "right_only").sum()),
                 score_identical_frac=float((ds < 1e-6).mean()), score_maxdiff=float(ds.max()), fold_same_frac=float((b_.fold == b_.fold_old).mean()),
                 event_same_frac=float((b_.event == b_.event_old).mean()))
    if s in ("fd_nd", "mb_nd"):
        sig = d[np.isin(d.truth, [0, 1]) & (d.hypo_swap == 0) & (d.hypo_from_tag == 0) & (d.dstar_conflict == 0)]
        r.update(n_sig_tag=len(sig), n_sig_tag_ge025=int((sig.Lxy >= 0.25).sum()), n_sig_tag_slice=int((sig.Lxy < 0.25).sum()),
                 mean_score_sig_ge025=float(sig.score[sig.Lxy >= 0.25].mean()), mean_score_sig_slice=float(sig.score[sig.Lxy < 0.25].mean()))
    R[s] = r
    print("==", s); [print(f"  {k}: {v}") for k, v in r.items()]
# zbiór treningowy D0_H: czy fold = hash(run, event) i jaki zakres numerów
import h5py
with h5py.File("/home/k/Documents/charm/runs/lxplus/models/D0_H/oof.h5") as f:
    y, run, ev, fo = f["y"][:], f["run"][:], f["event"][:], f["fold"][:]
t = dict(n=len(y), n_sig=int((y == 1).sum()), sig_runs=np.unique(run[y == 1]).tolist()[:5], sig_event_min=int(ev[y == 1].min()), sig_event_max=int(ev[y == 1].max()),
         bkg_runs=np.unique(run[y == 0]).tolist()[:5], fold_is_hash=bool((fold_of(run.astype(np.int64), ev.astype(np.int64), 5) == fo).all()))
R["train_D0_H"] = t; print("== train", t)
json.dump(R, open(f"{HERE}/out/check_inputs.json", "w"), indent=1, default=str)

#!/usr/bin/env python3
"""DNET etap 4 — diagnostyka wolnego pionu D* (dopisek 19.11, tylko MC fd_nd, sama inferencja zapisanymi modelami foldów).
Sieć z torami otoczenia widzi wolny pion z D*→D0π w ΔR < 0,4 od D0 ZAWSZE w próbce kalibracyjnej (tag D*), a w próbce inkluzywnej tylko czasem.
Dla każdego modelu D0 (sieci i drzewo H jako kontrola), per bin pT, przy WP70 (próg = kwantyl score sygnału inkluzywnego fd_nd w oknie, jak w analizie):
  (a) ε D0 z próbki tagu D* (wiersze p_: para D0←D* jak w kalibracji, prawda D* i D0 sygnałowa; jeden wiersz na D0),
  (b) ε pozostałych D0 sygnałowych (wiersze a_ nietrafione przez żadną parę p_ z prawdą D*),
  (c) jak (a), ale z zamaskowanym tokenem wolnego pionu (tor otoczenia o pT zgodnym z p_pis_pt; pT tokenu = 0 → maska 0),
  (i) ε inkluzywna (wszystkie a_; = WP z konstrukcji).
W tabelach nie ma flagi prawdy „D0 z D*” → (b) zawiera też D0←D* bez zrekonstruowanego tagu: (a)/(b) jest OGRANICZENIEM DOLNYM efektu.
Wyjście: JSON. Użycie: dn4_pis.py --data DIR --env DIR --models DIR --out OUT.json [--tags REGEX] [--threads n]"""
import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dn4_lib as L4  # noqa: E402
import dn_feat as DF  # noqa: E402

SP, S = "D0", "fd_nd"
WPS = [0.7, 0.5, 0.8]
COLS = ["X", "Dd", "m", "pt", "eta", "phi", "uid", "truth", "swap", "w", "row", "ftag"]


def load(ddir, edir, pre, extra=()):
    z = np.load(os.path.join(ddir, f"dn_{S}.npz"), allow_pickle=False)
    k = f"{SP}__{pre}"
    ix, E = [], []
    for f in sorted(glob.glob(os.path.join(edir, f"dn2_{S}_*.npz"))):
        e = np.load(f, allow_pickle=False)
        if k + "idx" in e.files:
            ix.append(e[k + "idx"]); E.append(e[k + "E"])
    ix, E = np.concatenate(ix), np.concatenate(E)
    o = np.argsort(ix)
    ix, E = ix[o], E[o]
    T = {c: z[k + c][ix] for c in COLS + list(extra)}
    T["X"] = T["X"][:, [DF.ALLF[SP].index(f) for f in DF.H]]
    T["E"] = E.copy()
    T["n_full"] = int(len(z[k + "m"]))
    return T


def eff(pas, w):
    sw = float(w.sum())
    if sw <= 0:
        return None
    e = float((w * pas).sum() / sw)
    neff = sw ** 2 / float((w ** 2).sum())
    return [e, float(np.sqrt(max(e * (1 - e), 1e-9) / neff)), int(len(w))]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--env", required=True); ap.add_argument("--models", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--tags", default=""); ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    t0 = time.time()
    m0, win = DF.cfg["species"][SP]["m0"], DF.cfg["species"][SP]["signal_window"]
    A = load(a.data, a.env, "a_")
    P = load(a.data, a.env, "p_", extra=("pis_pt", "truth_ds"))
    sig = lambda T: np.isin(T["truth"], DF.SIGNAL_CATS) & (T["swap"] == 0) & (np.abs(T["m"] - m0) < win)  # noqa: E731
    # wolny pion wśród tokenów otoczenia: najbliższy w pT (float16 → tolerancja 1 %)
    Ept = P["E"][:, :, 0].astype(np.float64)
    pis = P["pis_pt"].astype(np.float64)
    rel = np.abs(Ept - pis[:, None]) / np.clip(pis[:, None], 1e-3, None)
    rel[Ept <= 0] = np.inf
    j = np.argmin(rel, 1)
    rows = np.arange(len(j))
    found = rel[rows, j] < 0.01
    dr = np.hypot(P["E"][rows, j, 1].astype(np.float64), P["E"][rows, j, 2].astype(np.float64))
    seen = found & (dr < 0.4) & (j < 12)                       # token faktycznie widziany przez sieć e04
    Em = P["E"].copy()
    Em[rows[found], j[found], :] = 0
    # próbka tagu: prawda D* i D0, jeden wiersz na D0 (pierwsza para)
    sp_ = sig(P) & np.isin(P["truth_ds"], DF.SIGNAL_CATS)
    key = lambda T: T["ftag"].astype(np.int64) * (1 << 40) + T["row"].astype(np.int64)  # noqa: E731
    kp = key(P)
    _, first = np.unique(np.where(sp_, kp, -1 - rows), return_index=True)
    one = np.zeros(len(kp), bool); one[first] = True
    tagsel = sp_ & one
    sa = sig(A)
    in_tag = np.isin(key(A), kp[sp_])
    res = {"n_a": int(sa.sum()), "n_a_tagged": int((sa & in_tag).sum()), "n_p_tag": int(tagsel.sum()), "rows_env": {"a": [len(A["m"]), A["n_full"]], "p": [len(P["m"]), P["n_full"]]},
           "pis_found": float(found[tagsel].mean()), "pis_seen_e04": float(seen[tagsel].mean()), "pis_dr_lt04": float((dr[tagsel & found] < 0.4).mean()),
           "note": "(b) = sygnał a_ bez zrekonstruowanej pary D* — zawiera D0←D* bez tagu; (a)/(b) to ograniczenie dolne efektu", "models": {}}
    print(f"sygnał a_ {res['n_a']} (w tym z tagiem {res['n_a_tagged']}), próbka tagu {res['n_p_tag']}; wolny pion znaleziony {res['pis_found']:.3f}, "
          f"widziany przez e04 {res['pis_seen_e04']:.3f}", flush=True)
    metas = {os.path.basename(os.path.dirname(f)): json.load(open(f)) for f in sorted(glob.glob(os.path.join(a.models, "*", "meta.json")))}
    for t, meta in metas.items():
        if meta["sp"] != SP or (a.tags and not re.search(a.tags, t)):
            continue
        meta, ms = L4.load_models(a.models, t)
        if meta["model"] == "xgb":
            sA = L4.score_rows("xgb", ms, A["uid"], A["X"]); sP = L4.score_rows("xgb", ms, P["uid"], P["X"]); sM = sP
        else:
            import torch
            torch.set_num_threads(max(1, a.threads))
            f = lambda T, E: L4.score_rows("tr", ms, T["uid"], T["X"], T["Dd"], E, T["eta"], T["phi"])  # noqa: E731
            sA, sP, sM = f(A, A["E"]), f(P, P["E"]), f(P, Em)
        R = {}
        for lo, hi in DF.BINS:
            ba = sa & (A["pt"] >= lo) & (A["pt"] < hi); bp = tagsel & (P["pt"] >= lo) & (P["pt"] < hi)
            if ba.sum() < 50 or bp.sum() < 20:
                continue
            o = {}
            for wp in WPS:
                thr = float(np.quantile(sA[ba], 1 - wp))
                e_i, e_a = eff(sA[ba] > thr, A["w"][ba]), eff(sP[bp] > thr, P["w"][bp])
                e_b, e_c = eff(sA[ba & ~in_tag] > thr, A["w"][ba & ~in_tag]), eff(sM[bp] > thr, P["w"][bp])
                e_at = eff(sA[ba & in_tag] > thr, A["w"][ba & in_tag])
                r_ab = [e_a[0] / e_b[0], e_a[0] / e_b[0] * float(np.hypot(e_a[1] / e_a[0], e_b[1] / e_b[0]))] if e_a and e_b and e_b[0] > 0 else None
                o[f"{wp:.1f}"] = {"thr": thr, "incl": e_i, "a_tag": e_a, "a_tag_in_a": e_at, "b_rest": e_b, "c_masked": e_c, "a_over_b": r_ab,
                                  "c_over_a": (e_c[0] / e_a[0] if e_c and e_a and e_a[0] > 0 else None),
                                  "d_score_mask_mean": float(np.mean(sM[bp] - sP[bp])), "tag_frac": float((ba & in_tag).sum() / ba.sum())}
            R[f"{lo}-{hi}"] = o
            x = o["0.7"]
            print(f"{t} {lo}-{hi}: ε(a tag D*) {x['a_tag'][0]:.3f}±{x['a_tag'][1]:.3f}, ε(b reszta) {x['b_rest'][0]:.3f}±{x['b_rest'][1]:.3f}, "
                  f"a/b {x['a_over_b'][0]:.3f}±{x['a_over_b'][1]:.3f}, ε(c maska) {x['c_masked'][0]:.3f}, udział tagu {x['tag_frac']:.3f}", flush=True)
        res["models"][t] = R
    res["sec"] = round(time.time() - t0, 1)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)
    print(f"zapisano {a.out}, {res['sec']}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

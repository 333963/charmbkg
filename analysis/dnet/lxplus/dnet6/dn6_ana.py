#!/usr/bin/env python3
"""DNET etap 6 — P1 dla sieci z maską wolnego pionu, z GOTOWYCH score OOF (bez inferencji; tylko MC fd_nd) → dnet6_p1.json.
Definicje jak dn4_pis.py (etap 4), próg WP = kwantyl score sygnału inkluzywnego fd_nd w oknie, per bin pT:
  (a) ε D0 z próbki tagu D* (wiersze p_: prawda D* i D0, poprawna hipoteza, okno masy; jeden wiersz na D0),
  (b) ε pozostałych D0 sygnałowych (wiersze a_ bez zrekonstruowanej pary D* z prawdą),
  (i) ε inkluzywna; podwójny stosunek DR = [ε_a/ε_b]_model / [ε_a/ε_b]_drzewo H (błędy dwumianowe MC, kwadratowo).
Dla sieci z maską (c) ≡ (a) z konstrukcji (maska jest w score). Opisowo: ε sygnału a_ osobno dla wierszy z usuniętym torem / bez (nrm).
Modele: wszystkie z --sc pasujące do --tags + drzewo odniesienia; pełna sieć etapu 2 (ta sama pula/ziarna) jako porównanie „bez maski”.
Użycie: dn6_ana.py --data DIR1/data --env DIR2/parts --sc DIR[,DIR…] --out dnet6_p1.json [--tags REGEX]"""
import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_Q = os.path.dirname(HERE)
for p in (os.environ.get("DN1") or os.path.join(_Q, "dnet1"), HERE):
    sys.path.insert(0, p)
import dn_feat as DF  # noqa: E402
import dn6_mask as MK  # noqa: E402

SP, S = "D0", "fd_nd"
WPS = [0.7, 0.5, 0.8]
REF = "D0.xgb.H.sb.n.s1"
COLS = ["Dd", "m", "pt", "eta", "phi", "uid", "truth", "swap", "w", "row", "ftag"]


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
    T["idx"], T["E"], T["n_full"] = ix, E, int(len(z[k + "m"]))
    return T


def eff(pas, w):
    sw = float(w.sum())
    if sw <= 0 or len(w) == 0:
        return None
    e = float((w * pas).sum() / sw)
    neff = sw ** 2 / float((w ** 2).sum())
    return [e, float(np.sqrt(max(e * (1 - e), 1e-9) / neff)), int(len(w))]


def ratio(x, y):
    if not x or not y or y[0] <= 0 or x[0] <= 0:
        return None
    r = x[0] / y[0]
    return [r, r * float(np.hypot(x[1] / x[0], y[1] / y[0]))]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--env", required=True); ap.add_argument("--sc", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--tags", default=r"^D0\.tr\.H\+e04(np)?\.mc\.d[0-9]+\.s[123]$")
    a = ap.parse_args(argv)
    t0 = time.time()
    m0, win = DF.cfg["species"][SP]["m0"], DF.cfg["species"][SP]["signal_window"]
    A = load(a.data, a.env, "a_")
    P = load(a.data, a.env, "p_", extra=("pis_pt", "truth_ds"))
    sig = lambda T: np.isin(T["truth"], DF.SIGNAL_CATS) & (T["swap"] == 0) & (np.abs(T["m"] - m0) < win)  # noqa: E731
    rows = np.arange(len(P["m"]))
    sp_ = sig(P) & np.isin(P["truth_ds"], DF.SIGNAL_CATS)
    key = lambda T: T["ftag"].astype(np.int64) * (1 << 40) + T["row"].astype(np.int64)  # noqa: E731
    kp = key(P)
    _, first = np.unique(np.where(sp_, kp, -1 - rows), return_index=True)
    one = np.zeros(len(kp), bool); one[first] = True
    tagsel = sp_ & one
    sa = sig(A)
    in_tag = np.isin(key(A), kp[sp_])
    # maska: liczniki (opisowo) — wolny pion w tokenach jak dn4_pis (pT zgodne w 1 %), czy usunięty
    Ept = P["E"][:, :, 0].astype(np.float64); pis = P["pis_pt"].astype(np.float64)
    rel = np.abs(Ept - pis[:, None]) / np.clip(pis[:, None], 1e-3, None); rel[Ept <= 0] = np.inf
    j = np.argmin(rel, 1); found = rel[rows, j] < 0.01
    dr = np.hypot(P["E"][rows, j, 1].astype(np.float64), P["E"][rows, j, 2].astype(np.float64))
    seen = found & (dr < 0.4) & (j < 12)
    vP = MK.veto(P["E"], P["Dd"], P["pt"], P["eta"], P["phi"], P["m"])
    vA = MK.veto(A["E"], A["Dd"], A["pt"], A["eta"], A["phi"], A["m"])
    in04 = lambda T: (T["E"][:, :12, 0] > 0) & (np.hypot(T["E"][:, :12, 1].astype(np.float64), T["E"][:, :12, 2].astype(np.float64)) < 0.4)  # noqa: E731
    nrmA = (vA[:, :12] & in04(A)).sum(1)
    res = {"generated": time.strftime("%FT%T"), "mask": {"dm_max": MK.DM_MAX, "sign": MK.SIGN}, "ref": REF, "wps": WPS,
           "n_a": int(sa.sum()), "n_a_tagged": int((sa & in_tag).sum()), "n_p_tag": int(tagsel.sum()),
           "pis_found": float(found[tagsel].mean()), "pis_seen_e04": float(seen[tagsel].mean()),
           "pis_removed_of_seen": float(vP[rows, j][tagsel & seen].mean()), "pis_left_visible_of_tag": float((seen & ~vP[rows, j])[tagsel].mean()),
           "sig_incl_frac_rows_masked": float((nrmA[sa] > 0).mean()), "sig_untagged_frac_rows_masked": float((nrmA[sa & ~in_tag] > 0).mean()),
           "sig_incl_frac_tracks_removed": float(nrmA[sa].sum() / max(in04(A)[sa].sum(), 1)),
           "note": "(b) zawiera D0←D* bez zrekonstruowanego tagu; DR sieć/drzewo: 1 = sieć nie rozróżnia D0←D* od pozostałych ponad kinematykę widzianą przez drzewo",
           "models": {}}
    print(f"sygnał a_ {res['n_a']} (z tagiem {res['n_a_tagged']}), próbka tagu {res['n_p_tag']}; π_s widziany w e04 {res['pis_seen_e04']:.3f}, usunięty maską "
          f"{res['pis_removed_of_seen']:.4f} widzianych, zostaje widoczny w {res['pis_left_visible_of_tag']:.4f} próbki tagu; wiersze sygnału z usuniętym torem "
          f"{res['sig_incl_frac_rows_masked']:.3f} (bez tagu {res['sig_untagged_frac_rows_masked']:.3f})", flush=True)
    scf = {}
    for d in a.sc.split(","):
        for f in sorted(glob.glob(os.path.join(d, "sc", "*.npz"))):
            t = os.path.basename(f)[:-4]
            if t == REF or re.search(a.tags, t):
                scf[t] = f                                                   # późniejszy katalog wygrywa (jak analyze.py)
    if REF not in scf:
        raise SystemExit(f"brak score odniesienia {REF}")
    SC = {}
    for t, f in scf.items():
        z = np.load(f)
        if f"{S}__a" not in z.files or f"{S}__p" not in z.files:
            print(f"{t}: brak score fd_nd a/p — pomijam", flush=True); continue
        SC[t] = (z[f"{S}__a"][A["idx"]], z[f"{S}__p"][P["idx"]], f)
    for t, (sA, sP, f) in SC.items():
        okA, okP = np.isfinite(sA), np.isfinite(sP)
        R = {}
        for lo, hi in DF.BINS:
            ba = sa & okA & (A["pt"] >= lo) & (A["pt"] < hi); bp = tagsel & okP & (P["pt"] >= lo) & (P["pt"] < hi)
            if ba.sum() < 50 or bp.sum() < 20:
                continue
            o = {}
            for wp in WPS:
                thr = float(np.quantile(sA[ba], 1 - wp))
                e_i, e_a = eff(sA[ba] > thr, A["w"][ba]), eff(sP[bp] > thr, P["w"][bp])
                e_b = eff(sA[ba & ~in_tag] > thr, A["w"][ba & ~in_tag])
                e_at = eff(sA[ba & in_tag] > thr, A["w"][ba & in_tag])
                e_m = eff(sA[ba & ~in_tag & (nrmA > 0)] > thr, A["w"][ba & ~in_tag & (nrmA > 0)])
                e_u = eff(sA[ba & ~in_tag & (nrmA == 0)] > thr, A["w"][ba & ~in_tag & (nrmA == 0)])
                o[f"{wp:.1f}"] = {"thr": thr, "incl": e_i, "a_tag": e_a, "a_tag_in_a": e_at, "b_rest": e_b, "a_over_b": ratio(e_a, e_b),
                                  "b_rest_masked_rows": e_m, "b_rest_unmasked_rows": e_u, "tag_frac": float((ba & in_tag).sum() / ba.sum())}
            R[f"{lo}-{hi}"] = o
        res["models"][t] = {"file": f, "bins": R}
    ref = res["models"][REF]["bins"]
    for t, M in res["models"].items():
        for b, o in M["bins"].items():
            for wp, x in o.items():
                y = ref.get(b, {}).get(wp)
                if x["a_over_b"] and y and y["a_over_b"]:
                    dr_ = x["a_over_b"][0] / y["a_over_b"][0]
                    x["double_ratio"] = [dr_, dr_ * float(np.hypot(x["a_over_b"][1] / x["a_over_b"][0], y["a_over_b"][1] / y["a_over_b"][0]))]
            x = o.get("0.7")
            if x and x.get("double_ratio"):
                print(f"{t} {b}: ε_a {x['a_tag'][0]:.3f}±{x['a_tag'][1]:.3f}, ε_b {x['b_rest'][0]:.3f}±{x['b_rest'][1]:.3f}, a/b {x['a_over_b'][0]:.3f}±{x['a_over_b'][1]:.3f}, "
                      f"DR {x['double_ratio'][0]:.4f}±{x['double_ratio'][1]:.4f}", flush=True)
    # bramka dla węzłów app: P1 „maska działa” dla ziarna 1 sieci z maską (KRYTERIA_etap6.md §2)
    g1 = [t for t in res["models"] if re.search(r"e04np\..*\.s1$", t)]
    gate = None
    if g1:
        v = [res["models"][g1[0]]["bins"].get(f"{lo}-{hi}", {}).get("0.7", {}).get("double_ratio") for lo, hi in DF.BINS]
        gate = bool(all(x is not None and abs(x[0] - 1) <= 0.03 for x in v))
    res["p1_gate_seed1"] = gate
    res["sec"] = round(time.time() - t0, 1)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)
    print(f"zapisano {a.out}, {res['sec']}s; bramka P1 (ziarno 1): {gate}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

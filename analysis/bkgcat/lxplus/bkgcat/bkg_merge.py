#!/usr/bin/env python3
"""BKGCAT — suma plików bkg_<NNNNNN>.npz (bkg_hist.py) → bkgcat.npz. Tylko numpy.
Wyjście: keys [nK] (kategorie), tab [N, 6] int64 (gatunek 0 = D0 / 1 = Dplus, bin pT 0..2, klucz, przedział score s, bin masy 1 MeV, n), sw [N] (suma wag),
  br{0,1,2}_key/_n (region: wszystko / |η| < 2,5 / pT ≥ 3 GeV ∧ |η| < 2,2; key = ((kod·8 + klasa braku)·2 + fromB)·4 + {421, 411, 413, inne}),
  mm_key/_n ((kod·4 + pdgidx)·400 + bin m_brak 5 MeV), dty10_mc, dty4_mc [5], thr|<gatunek>|<bin>, meta (json: pliki, statystyki dopasowań).
Użycie: bkg_merge.py --parts DIR --thr wps_hist.npz --out bkgcat.npz"""
import argparse
import glob
import json
import os
import sys

import numpy as np

TAG = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}


def agg(keys, vals):
    u, inv = np.unique(keys, axis=0, return_inverse=True)
    inv = np.ravel(inv)
    return u, [np.bincount(inv, weights=v, minlength=len(u)) for v in vals]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", required=True); ap.add_argument("--thr", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    fs = sorted(glob.glob(os.path.join(a.parts, "bkg_*.npz")))
    kid, tabs, sws, metas = {}, [], [], []
    br = {k: ([], []) for k in ("br0", "br1", "br2", "mm")}
    dty = {"dty10_mc": np.zeros(5, np.int64), "dty4_mc": np.zeros(5, np.int64)}
    codes = {"dty10_codes": set(), "dty4_codes": set()}
    other = 0
    for f in fs:
        z = np.load(f, allow_pickle=False)
        loc = np.array([kid.setdefault(str(k), len(kid)) for k in z["keys"]], np.int64)
        t = z["tab"].copy()
        if len(t):
            t[:, 2] = loc[t[:, 2]]
        tabs.append(t); sws.append(z["sw"]); metas.append(json.loads(str(z["meta"])))
        for k in br:
            br[k][0].append(z[k + "_key"]); br[k][1].append(z[k + "_n"])
        for k in dty:
            dty[k] += z[k]
        for k in codes:
            codes[k] |= set(z[k].tolist())
        other += int(z["pdg_other"][0])
    T = np.concatenate(tabs)
    u, (n, sw) = agg(T[:, :5], [T[:, 5].astype(np.float64), np.concatenate(sws)])
    out = {"keys": np.array(sorted(kid, key=kid.get)), "tab": np.concatenate([u, np.round(n).astype(np.int64)[:, None]], axis=1), "sw": sw}
    for k in br:
        uu, (nn,) = agg(np.concatenate(br[k][0]), [np.concatenate(br[k][1]).astype(np.float64)])
        out[k + "_key"], out[k + "_n"] = uu, np.round(nn).astype(np.int64)
    out.update(dty)
    for k in codes:
        out[k] = np.array(sorted(codes[k]), np.int64)
    Z = np.load(a.thr, allow_pickle=False)
    for sp, tag in TAG.items():
        for b in ("3-5", "5-8", "8-30"):
            out[f"thr|{sp}|{b}"] = Z[f"{sp}|{tag}|{b}|thr"]
    out["meta"] = np.array(json.dumps({"n_files": len(fs), "pdg_other": other, "eps_desc": json.loads(str(Z["meta"]))["eps_desc"], "tag": TAG, "parts": metas}))
    np.savez_compressed(a.out, **out)
    print(f"{len(fs)} plików, {len(kid)} kategorii, {len(out['tab'])} komórek, {int(out['tab'][:, 5].sum())} kandydatów → {a.out} ({os.path.getsize(a.out) / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

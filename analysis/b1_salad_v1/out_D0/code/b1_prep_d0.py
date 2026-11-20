#!/usr/bin/env python
"""B1-D0 prep: S_bkggen/sub/*.npz (surowe kolumny h5 v1.1, decayType 10 = D0) → kompaktowe npz z cechami H
(charmtag.features.compute), m, m_swap, pt, truth_category, parity. Jak b1_prep.py (D⁺), ale: zakres 1,75–2,00, pT 3–30,
weto refleksji |m_swap−m0| > 25 MeV WSZĘDZIE (S4), swapy 100/101 z closest_mass (mark_swaps w load_npz)."""
import os, sys, time
import numpy as np, pandas as pd
SRC = os.environ.get("CHARMTAG_SRC", "/work/b1d0/in/src")
sys.path.insert(0, SRC)
from charmtag import config as C
from charmtag.bkggen import data as D
from charmtag.kinematics import swapped_mass_sp

cfg = C.load(os.path.join(SRC, "config", "default.yaml"))
H = list(cfg["features"]["frozen7"]) + list(cfg["features"]["extra_H"])
UP = "/work/uploads/charm/charm_work/S_bkggen/sub"
OUT = "/work/b1d0/work/sub"; os.makedirs(OUT, exist_ok=True)
sp = cfg["species"]["D0"]; lo, hi = sp["mass_range"]; M0 = sp["m0"]; REFL = 0.025

SETS = {"mcbkg": ["mb_p0.npz", "mb_p1.npz", "mb_p2.npz", "mb_p3.npz"], "sig": ["fd_p0.npz", "fd_p1.npz"],
        "data": ["data_p0.npz", "data_p1.npz"]}
which = sys.argv[1:] or list(SETS)
for name in which:
    t0 = time.time(); parts = []; stats = {}
    for p in SETS[name]:
        L = D.load_npz(cfg, f"{UP}/{p}", "D0", tag=name)
        m = L["m"].to_numpy(float); pt = L["pt"].to_numpy(float)
        ms = swapped_mass_sp(L, sp); L["m_swap"] = ms
        n0 = len(L)
        keep = (m > lo) & (m < hi) & (pt >= 3) & (pt < 30)
        n_range = int(keep.sum())
        veto = np.abs(ms - M0) > REFL
        keep &= veto
        if name == "sig":
            keep &= L["truth"].isin([0, 1, 100, 101]).to_numpy()
        if name == "mcbkg":
            keep &= (L["truth"].to_numpy(float) >= 0)   # wszystko (0/1 sygnał w mb, 100/101 swapy, 2–7, ≥8 tło)
        stats[p] = dict(n_raw=n0, n_range=n_range, n_veto=int(keep.sum()), veto_eff_range=float(veto[(m > lo) & (m < hi) & (pt >= 3) & (pt < 30)].mean()))
        L = L[keep].reset_index(drop=True)
        ev = L["event"].to_numpy() if "event" in L else np.arange(len(L))
        out = {f: L[f].to_numpy(np.float32) for f in H}
        out.update(m=L["m"].to_numpy(np.float32), m_swap=L["m_swap"].to_numpy(np.float32), pt=L["pt"].to_numpy(np.float32),
                   truth=L["truth"].to_numpy(np.int16), parity=(ev % 2).astype(np.int8), w=L["weight"].to_numpy(np.float32),
                   src=np.full(len(L), SETS[name].index(p), np.int8))
        parts.append(pd.DataFrame(out)); del L
        print(name, p, stats[p], f"{time.time()-t0:.0f}s", flush=True)
    A = pd.concat(parts, ignore_index=True)
    np.savez(f"{OUT}/{name}.npz", **{c: A[c].to_numpy() for c in A.columns})
    tc = {int(k): int(v) for k, v in zip(*np.unique(A["truth"], return_counts=True))}
    print(name, len(A), "truth:", tc, "parity1:", float(A["parity"].mean()), f"{time.time()-t0:.0f}s", flush=True)

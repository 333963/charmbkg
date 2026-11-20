#!/usr/bin/env python
"""B1 prep: dp_*.npz (surowe kolumny h5 v1.1) → kompaktowe npz z cechami H (charmtag.features.compute), m, pt,
truth_category, parity. Weto D* i zakres mas jak charmtag (config/default.yaml)."""
import os, sys, time
import numpy as np, pandas as pd
SRC = os.environ.get("CHARMTAG_SRC", "/work/b1/in/AA_merge/src")
sys.path.insert(0, SRC)
from charmtag import config as C
from charmtag.bkggen import data as D
from charmtag.samples import dstar_veto_mask

cfg = C.load(os.path.join(SRC, "config", "default.yaml"))
H = list(cfg["features"]["frozen7"]) + list(cfg["features"]["extra_H"])
UP = "/work/uploads/charm/charm_work/S6_kand/sub"
OUT = "/work/b1/work/sub"; os.makedirs(OUT, exist_ok=True)
sp = cfg["species"]["Dplus"]; lo, hi = sp["mass_range"]

for name, pats in (("mcbkg", ["dp_mb_p0.npz", "dp_mb_p1.npz"]), ("sig", ["dp_fd_p0.npz", "dp_fd_p1.npz"]), ("data", ["dp_data_p0.npz"])):
    t0 = time.time()
    L = pd.concat([D.load_npz(cfg, f"{UP}/{p}", "Dplus", tag=name) for p in pats], ignore_index=True)
    m = L["m"].to_numpy(float)
    keep = (m > lo) & (m < hi) & dstar_veto_mask(L, sp)
    if name == "sig":
        keep &= L["truth"].isin([0, 1]).to_numpy()
    if name == "mcbkg":
        keep &= (L["truth"].to_numpy(float) >= 2)          # truth 2–7 zostają (wariant adv.), tło = ≥8
    L = L[keep].reset_index(drop=True)
    ev = L["event"].to_numpy() if "event" in L else np.arange(len(L))
    rng = np.random.default_rng(len(L))
    par = (ev % 2).astype(np.int8) if len(np.unique(ev)) > len(L) // 10 else (rng.random(len(L)) < 0.5).astype(np.int8)
    out = {f: L[f].to_numpy(np.float32) for f in H}
    out.update(m=L["m"].to_numpy(np.float32), pt=L["pt"].to_numpy(np.float32), truth=L["truth"].to_numpy(np.int16),
               parity=par, w=L["weight"].to_numpy(np.float32))
    np.savez(f"{OUT}/{name}.npz", **out)
    tc = {int(k): int(v) for k, v in zip(*np.unique(out["truth"], return_counts=True))}
    print(name, len(L), "parity_by_event" if par.dtype == np.int8 and len(np.unique(ev)) > len(L) // 10 else "parity_random",
          "truth:", tc, f"{time.time()-t0:.0f}s", flush=True)

"""DLZ_dplus/DATA: po treningu — progi WP (thr.json, z sygnału MC OOF), macierz cech sygnału MC (oof/S_feat.npz),
slim lokalnej ramki danych 428770 v1.2 (2 pliki, score OOF) w formacie apply.py → slim/data_428770loc.npz."""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dlzlib as D  # noqa: E402
from analyze import thresholds  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FR = "/home/k/Documents/charm/charm_work/FEAT_charmnet/frames"
FEATS = D.H + D.DLZ8
b = np.load(os.path.join(BASE, "oof", "base.npz"))
tags = [f[:-4] for f in sorted(os.listdir(os.path.join(BASE, "oof"))) if f.endswith(".npz") and "_s" in f and f[0] == "H"]
T = {}
sc = {}
for t in tags:
    z = np.load(os.path.join(BASE, "oof", t + ".npz"))
    sc[t] = z
    th = thresholds(z["S"], b["S_w"].astype(float), b["S_pt"].astype(float))
    T[t] = {f"{lo}-{hi}": {f"{wp:.2f}": th[(lo, hi, wp)] for wp in D.WPS} for lo, hi in D.BINS}
json.dump(T, open(os.path.join(BASE, "thr.json"), "w"), indent=1)
print(json.dumps(T, indent=0)[:600])
fd = pd.read_pickle(os.path.join(FR, "fd_nd_Dplus.pkl"))
fd = fd[fd["dstar_ok"].to_numpy(bool)].reset_index(drop=True)
S = fd[np.isin(fd["truth"].to_numpy(), [0, 1])]
assert len(S) == len(b["S_m"]) and np.allclose(S["m"].to_numpy(), b["S_m"])
np.savez_compressed(os.path.join(BASE, "oof", "S_feat.npz"), feat=S[FEATS].to_numpy(np.float32), names=np.array(FEATS))
# korelacja m z cechami i score w sygnale MC (sPlot: niezależność)
rho = {f: float(np.corrcoef(S["m"].to_numpy(float), np.nan_to_num(S[f].to_numpy(float)))[0, 1]) for f in FEATS}
for t in tags:
    rho["score_" + t] = float(np.corrcoef(b["S_m"].astype(float), sc[t]["S"])[0, 1])
json.dump(rho, open(os.path.join(BASE, "out_rho_signal_mc.json"), "w"), indent=1)
print({k: round(v, 3) for k, v in rho.items() if abs(v) > 0.03})
del fd
da = pd.read_pickle(os.path.join(FR, "data_428770_Dplus.pkl"))
da = da[da["dstar_ok"].to_numpy(bool)].reset_index(drop=True)
assert np.allclose(da["m"].to_numpy(), b["d_m"])
pt = da["pt"].to_numpy(float)
ib = D.ptbin(pt)
E1 = np.linspace(1.7, 2.2, 501)
hall = np.array([np.histogram(da["m"].to_numpy(float)[ib == i], E1)[0] for i in range(3)])
loose = np.zeros(len(da), bool)
for t in tags:
    if t.endswith("_s12345"):
        loose |= sc[t]["dloc"] >= np.array([T[t][f"{lo}-{hi}"]["0.95"] for lo, hi in D.BINS])[ib]
keepH = sc["H_s12345"]["dloc"] >= np.array([T["H_s12345"][f"{lo}-{hi}"]["0.95"] for lo, hi in D.BINS])[ib]
out = {"m": da["m"].to_numpy(np.float64)[loose], "pt": pt[loose].astype(np.float32), "run": da["run"].to_numpy(np.int64)[loose],
       "event": da["event"].to_numpy().astype(np.int64)[loose], "fold": b["d_fold"][loose]}
for t in tags:
    out["s_" + t] = sc[t]["dloc"][loose]
out["fi"] = np.flatnonzero(keepH[loose])
out["feat"] = da.loc[loose & keepH, FEATS].to_numpy(np.float32)
out["feat_names"] = np.array(FEATS)
out["n_files"] = 2
out["hall_any"] = hall
out["hall_mb"] = hall * 0
os.makedirs(os.path.join(BASE, "slim"), exist_ok=True)
np.savez_compressed(os.path.join(BASE, "slim", "data_428770loc.npz"), **out)
print("428770loc: po weto", len(da), "luźnych", int(loose.sum()), "z cechami", int((loose & keepH).sum()))

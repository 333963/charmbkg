"""FIT_scorebins — przygotowanie wejść (ŚLEPO: z danych tylko pasma boczne).

Dla gatunku i binu pT: komórki score wg kwantyli sygnału MC (fd, prawda 0/1), szablony masy sygnału (+ refleksja D0) per komórka,
tło: exp(pol2) dopasowane do pasm bocznych DANYCH per komórka i ekstrapolowane pod pik, N_S oczekiwane = S_MC(mb) × B_SB(dane)/B_SB(mb).
Wynik: out/inputs_<gatunek>.npz
"""
import os, sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fsb_common import EFF_EDGES, NCELL, PTB, xbasis, Model, init_beta

W = os.path.expanduser("~/Documents/charm/charm_work")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
CFG = {"D0": dict(m0=1.86484, lo=1.75, hi=2.00, nm=50, sb=0.075),
       "Dplus": dict(m0=1.86966, lo=1.72, hi=2.02, nm=60, sb=0.070)}


def load_D0():
    d = {}
    for s in ("data_428770", "fd_nd", "mb_nd"):
        z = np.load(f"{W}/B_trig/runs/lxplus/bt_{s}.npz", allow_pickle=True)
        d[s] = {k: z["all_" + k] for k in ("m", "pt", "score", "truth")}
    return d["data_428770"], d["fd_nd"], d["mb_nd"], "B_trig/runs/lxplus/bt_*.npz (all_*; dane i mb: 30 % zdarzeń)"


def load_Dplus():
    sc = np.load(f"{W}/DLZ_dplus/out/Dplus_mc/scores_H.npz")
    out = {}
    for s in ("data_428770", "fd_nd", "mb_nd"):
        L = pd.read_pickle(f"{W}/DLZ_dplus/frames/{s}_Dplus.pkl")[["m", "pt", "truth", "dstar_ok"]]
        L = L[L["dstar_ok"].to_numpy(bool)].reset_index(drop=True)
        out[s] = {k: L[k].to_numpy(np.float32) for k in ("m", "pt", "truth")}
    fd = out["fd_nd"]
    sig = np.isin(fd["truth"], [0, 1])
    assert sig.sum() == len(sc["S"]) and len(out["data_428770"]["m"]) == len(sc["dall"]), (sig.sum(), len(sc["S"]))
    fd = {k: v[sig] for k, v in fd.items()}
    fd["score"] = sc["S"]
    out["data_428770"]["score"] = sc["dall"]
    return out["data_428770"], fd, out["mb_nd"], "DLZ_dplus/frames/*_Dplus.pkl + out/Dplus_mc/scores_H.npz (H, OOF)"


def main(sp):
    c = CFG[sp]
    da, fd, mb, src = (load_D0 if sp == "D0" else load_Dplus)()
    edges = np.linspace(c["lo"], c["hi"], c["nm"] + 1)
    cen = 0.5 * (edges[1:] + edges[:-1])
    sbmask = (np.abs(cen - c["m0"]) > c["sb"]).astype(float)
    X = xbasis(c["nm"])
    res = {"X": X, "edges": edges, "sbmask": sbmask, "eff_edges": EFF_EDGES}
    info = {"src": src, "bins": {}}
    for ib, (lo, hi) in enumerate(PTB):
        ks = (fd["pt"] >= lo) & (fd["pt"] < hi) & np.isin(fd["truth"], [0, 1]) & (fd["m"] >= c["lo"]) & (fd["m"] < c["hi"])
        kr = (fd["pt"] >= lo) & (fd["pt"] < hi) & np.isin(fd["truth"], [100, 101]) & (fd["m"] >= c["lo"]) & (fd["m"] < c["hi"])
        thr = np.quantile(fd["score"][ks], 1 - EFF_EDGES)            # próg score dla każdej krawędzi wydajności
        thr[0], thr[-1] = -np.inf, np.inf
        cell = lambda s: np.clip(np.searchsorted(thr, s, side="right") - 1, 0, NCELL - 1)
        hs = np.zeros((NCELL, c["nm"])); hr = np.zeros_like(hs)
        np.add.at(hs, (cell(fd["score"][ks]), np.clip(np.digitize(fd["m"][ks], edges) - 1, 0, c["nm"] - 1)), 1)
        if kr.any():
            np.add.at(hr, (cell(fd["score"][kr]), np.clip(np.digitize(fd["m"][kr], edges) - 1, 0, c["nm"] - 1)), 1)
        f = hs.sum(1) / hs.sum()
        Us = hs / hs.sum(1, keepdims=True)
        Ur = hr / hs.sum(1, keepdims=True)
        # dane: TYLKO pasma boczne
        kd = (da["pt"] >= lo) & (da["pt"] < hi) & (np.abs(da["m"] - c["m0"]) > c["sb"]) & (da["m"] >= c["lo"]) & (da["m"] < c["hi"])
        hd = np.zeros((NCELL, c["nm"]))
        np.add.at(hd, (cell(da["score"][kd]), np.clip(np.digitize(da["m"][kd], edges) - 1, 0, c["nm"] - 1)), 1)
        hd *= sbmask[None, :]
        B = np.zeros_like(hd); chi = []
        for j in range(NCELL):
            m = Model(np.zeros((1, c["nm"])), X, "none", mask=sbmask[None, :])
            th0 = init_beta(hd[j] + 0.5, X, sbmask)
            th, _, _, ok, _ = m.fit(hd[j][None, :], th0)
            B[j] = np.exp(X @ th)
            chi.append(float((sbmask * (hd[j] - B[j]) ** 2 / np.maximum(B[j], 1e-9)).sum() / (sbmask.sum() - 3)))
            assert ok, (sp, ib, j)
        # normalizacja S: MC minbias (prawda 0/1) × stosunek pasm bocznych dane/MC (przed cięciem na score)
        kmb = (mb["pt"] >= lo) & (mb["pt"] < hi) & (mb["m"] >= c["lo"]) & (mb["m"] < c["hi"])
        s_mb = int((kmb & np.isin(mb["truth"], [0, 1])).sum())
        sb_mb = int((kmb & (np.abs(mb["m"] - c["m0"]) > c["sb"])).sum())
        sb_da = int(kd.sum())
        NS = s_mb * sb_da / sb_mb
        res[f"f_{ib}"], res[f"Us_{ib}"], res[f"Ur_{ib}"], res[f"B_{ib}"], res[f"NS_{ib}"], res[f"thr_{ib}"] = f, Us, Ur, B, NS, thr
        res[f"hd_{ib}"] = hd
        sr = np.abs(cen - c["m0"]) < 0.055
        info["bins"][f"{lo}-{hi}"] = dict(N_fd_sig=int(ks.sum()), N_fd_refl=int(kr.sum()), S_mb=s_mb, SB_mb=sb_mb, SB_data=sb_da,
                                         NS_exp=NS, B_SR_all=float(B[:, sr].sum()), chi2_sb_cells=[round(x, 2) for x in chi],
                                         B_SR_cells=[round(float(x), 1) for x in B[:, sr].sum(1)],
                                         thr=[float(x) for x in thr[1:-1]])
        print(sp, lo, hi, json.dumps({k: v for k, v in info["bins"][f"{lo}-{hi}"].items() if k not in ("thr",)}))
    np.savez_compressed(f"{OUT}/inputs_{sp}.npz", **res)
    json.dump(info, open(f"{OUT}/inputs_{sp}.json", "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])

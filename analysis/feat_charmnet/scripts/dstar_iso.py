"""FEAT_charmnet: izolacja PRAWDZIWEGO sygnału w danych — D*⁺→D⁰π_s (decayType 1), odjęcie tła z boków Δm,
vs sygnał D* w MC fd (truth 0/1) i vs tło. Okna jak charmbkg dstar_control: Δm sygnał 143,5–147,5 MeV, bok 150–165 MeV.
Δm w h5 v1.2 w GeV; Dmeson_mass = m(D*) → m(D⁰) = m − Δm, |m(D⁰) − m0| < 40 MeV. Użycie: python dstar_iso.py
"""
import glob
import json
import os

import h5py
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLS = ["Dmeson_decayType", "Dmeson_pt", "Dmeson_delta_mass", "Dmeson_mass", "Dmeson_ptcone40", "Dmeson_ptcone100",
        "Dmeson_truth_category"]


def read(sample):
    out = {c: [] for c in COLS}
    for p in sorted(glob.glob(os.path.join(BASE, "data", sample, "*_none.h5"))):
        a = h5py.File(p, "r")["events"].fields(COLS)[:]
        s = (a["Dmeson_decayType"] == 1) & (a["Dmeson_pt"] >= 3) & (a["Dmeson_pt"] < 30)
        for c in COLS:
            out[c].append(a[c][s])
    return {c: np.concatenate(v) for c, v in out.items()}


def main():
    res = {}
    D, M = read("data_428770"), read("fd_nd")
    bins = np.r_[np.linspace(0, 40, 161), 1e4]
    for lo, hi in ((3, 30), (3, 5), (5, 8), (8, 30)):
        r = {}
        for nm, X in (("data", D), ("fd", M)):
            pt, dm = X["Dmeson_pt"], X["Dmeson_delta_mass"]
            base = (pt >= lo) & (pt < hi) & (np.abs(X["Dmeson_mass"] - X["Dmeson_delta_mass"] - 1.86484) < 0.040)
            sr = base & (dm > 0.1435) & (dm < 0.1475)
            sb = base & (dm > 0.150) & (dm < 0.165)
            k = 4.0 / 15.0  # szerokości okien (płaskie tło w Δm — przybliżenie; tło rośnie ~√, zawyżenie ułamka tła)
            for c in ("Dmeson_ptcone40", "Dmeson_ptcone100"):
                hs, _ = np.histogram(X[c][sr], bins)
                hb, _ = np.histogram(X[c][sb], bins)
                sig = hs - k * hb
                if nm == "fd":
                    t = X["Dmeson_truth_category"]
                    hs_t, _ = np.histogram(X[c][sr & np.isin(t, [0, 1])], bins)
                    r[f"fd_truth_{c}_median"] = float(np.median(X[c][sr & np.isin(t, [0, 1])]))
                    r[f"fd_truth_{c}_frac_lt2"] = float((X[c][sr & np.isin(t, [0, 1])] < 2).mean())
                cdf = np.cumsum(sig) / max(sig.sum(), 1e-9)
                r[f"{nm}_{c}_median_sbsub"] = float(np.interp(0.5, cdf, bins[1:]))
                r[f"{nm}_{c}_frac_lt2_sbsub"] = float(sig[:8].sum() / max(sig.sum(), 1e-9))
                r[f"{nm}_{c}_bkg_median"] = float(np.median(X[c][sb])) if sb.any() else None
            r[f"{nm}_N_sr"], r[f"{nm}_N_sig_est"] = int(sr.sum()), float(sr.sum() - k * sb.sum())
        res[f"{lo}-{hi}"] = r
        print(lo, hi, json.dumps(r), flush=True)
    json.dump(res, open(os.path.join(BASE, "out", "dstar_iso.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

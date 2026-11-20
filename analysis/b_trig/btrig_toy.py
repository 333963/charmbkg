#!/usr/bin/env python3
"""B_trig — toy pass/fail Δm przy statystyce klasy triggera: obciążenie, pull, σ(ε) i moc testu SF95 = 0,83 vs 1.

Tło i sygnał ze stanu nominalnego EV (pełna lumi 428770, bin pT, okno ±55 MeV, WP95): liczby pass/fail w zakresie fitu z
`EV_dstar/lxplus/ev_data.npz`, S i ε z `ev_analysis_cm.json`; klasa = (f_S, f_B) — ułamki sygnału i tła nominału
(z btrig_sens: S_klasa/S_all, B_klasa/B_all) albo ręcznie. Kształt tła: empiryczny Δm danych (pass/fail osobno) z wyciętym
pikiem ±3 MeV i dopełnieniem interpolacją liniową. Kształt sygnału: 2 Gaussy z MC (shape_mc).
Użycie: btrig_toy.py --bin 5-8 --cls MB427929:0.25:0.08 --cls J12only:0.7:0.85 [--eps-true 0.79] [--ntoy 40] --out DIR
"""
import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import btrig_common as C  # noqa: E402

C.add_paths()
EV = os.path.join(C.CW, "EV_dstar", "lxplus")


def bkg_template(dm_mev):
    h, _ = np.histogram(dm_mev, C.E_DM)
    c = 0.5 * (C.E_DM[1:] + C.E_DM[:-1])
    pk = np.abs(c - 145.45) < 3.0
    h = h.astype(float)
    h[pk] = np.interp(c[pk], c[~pk], h[~pk])
    return h / h.sum()


def sig_template(shape):
    from scipy.stats import norm
    e = C.E_DM
    cdf = shape["f"] * norm.cdf(e, shape["mu"], shape["sc"]) + (1 - shape["f"]) * norm.cdf(e, shape["mu"], shape["sc"] * shape["k"])
    p = np.diff(cdf)
    return p / p.sum()


def _one(job):
    seed, S, eps, Bp, Bf, tp, tf, ts, shape = job
    rng = np.random.default_rng(seed)
    sp = rng.binomial(int(round(S)), eps)
    hp = rng.poisson(Bp * tp) + rng.multinomial(sp, ts)
    hf = rng.poisson(Bf * tf) + rng.multinomial(int(round(S)) - sp, ts)
    from charmtag import dstar_refit as DR
    try:
        r = DR.fit_passfail(hp, hf, C.E_DM, shape, interval=False, model_check=False)
        return float(r["eps"]), float(r["eps_err"]), bool(r["determined"])
    except Exception:  # noqa: BLE001
        return np.nan, np.nan, False


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default="5-8")
    ap.add_argument("--cls", action="append", default=[], help="nazwa:f_S:f_B (ułamki sygnału / tła nominału 428770 all)")
    ap.add_argument("--eps-true", type=float, default=None, help="domyślnie ε_dane nominału (SF95 0,83 → 0,79 w 5–8)")
    ap.add_argument("--eps-mc", type=float, default=0.947)
    ap.add_argument("--ntoy", type=int, default=40)
    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    js = json.load(open(os.path.join(EV, "ev_analysis_cm.json")))
    R = js["bins"][a.bin]
    shape = R["shape_mc"]
    nom = R["regions"]["nom"]
    S_nom = float(nom["N_sig_data"]["N"])
    eps_nom = float(nom["wp"]["0.95"]["eps_data"]) if a.eps_true is None else a.eps_true
    z = np.load(os.path.join(EV, "ev_data.npz"))
    lo, hi = C.BINS[a.bin]
    thr = js["thresholds"][a.bin]["0.95"]
    m = (z["pt"] >= lo) & (z["pt"] < hi) & (np.abs(z["m"] - js["m0"]) < C.NOM_WIN)
    dm = z["dm"][m] * 1e3
    pas = z["score"][m] >= thr
    tp, tf = bkg_template(dm[pas]), bkg_template(dm[~pas])
    Bp_nom = float(pas.sum()) - S_nom * eps_nom
    Bf_nom = float((~pas).sum()) - S_nom * (1 - eps_nom)
    ts = sig_template(shape)
    out = {"bin": a.bin, "S_nom": S_nom, "eps_true": eps_nom, "Bp_nom": Bp_nom, "Bf_nom": Bf_nom, "classes": {}}
    specs = a.cls or ["all:1:1"]
    with Pool(a.nproc) as P:
        for spec in specs:
            name, fs, fb = spec.split(":")
            fs, fb = float(fs), float(fb)
            jobs = [(1000 * i + 7, S_nom * fs, eps_nom, Bp_nom * fb, Bf_nom * fb, tp, tf, ts, shape) for i in range(a.ntoy)]
            res = np.array([x for x in P.map(_one, jobs)], dtype=object)
            e = res[:, 0].astype(float)
            er = res[:, 1].astype(float)
            det = res[:, 2].astype(bool)
            ok = np.isfinite(e) & det
            pull = (e[ok] - eps_nom) / er[ok]
            sf_err = float(np.std(e[ok]) / a.eps_mc) if ok.sum() > 2 else np.nan
            out["classes"][name] = {"f_S": fs, "f_B": fb, "S": S_nom * fs, "n_ok": int(ok.sum()), "n_toy": a.ntoy,
                                    "bias_eps": float(np.mean(e[ok]) - eps_nom) if ok.any() else np.nan,
                                    "rms_eps": float(np.std(e[ok])) if ok.any() else np.nan,
                                    "mean_err": float(np.mean(er[ok])) if ok.any() else np.nan,
                                    "pull_mean": float(np.mean(pull)) if ok.any() else np.nan,
                                    "pull_rms": float(np.std(pull)) if ok.any() else np.nan,
                                    "sigma_SF": sf_err, "Z_083_vs_1": 0.17 / sf_err if sf_err and np.isfinite(sf_err) else np.nan}
            print(name, out["classes"][name], flush=True)
    os.makedirs(a.out, exist_ok=True)
    json.dump(out, open(os.path.join(a.out, f"btrig_toy_{a.bin}.json"), "w"), indent=1, default=float)
    L = [f"# B_trig — toy pass/fail, bin {a.bin}, WP95 (S_nom = {S_nom:.0f}, ε = {eps_nom:.3f})", "",
         "| klasa | f_S | f_B | det/toy | obciążenie ε | rms ε | ⟨σ_HESSE⟩ | pull μ/σ | σ(SF) | Z(0,83 vs 1) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for k, v in out["classes"].items():
        L.append(f"| {k} | {v['f_S']} | {v['f_B']} | {v['n_ok']}/{v['n_toy']} | {v['bias_eps']:+.3f} | {v['rms_eps']:.3f} | {v['mean_err']:.3f} | "
                 f"{v['pull_mean']:+.2f}/{v['pull_rms']:.2f} | {v['sigma_SF']:.3f} | {v['Z_083_vs_1']:.1f} |")
    open(os.path.join(a.out, f"btrig_toy_{a.bin}.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    sys.exit(main())

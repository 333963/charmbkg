"""D23/ANALIZA krok 4: D⁰ inkluzywne i D* — tylko liczniki/istotność z histogramów pełnej statystyki (bez klasyfikatora).
Dodatkowo: domknięcie fitu D⁺ na fd_nd (N_fit vs prawda, udział klasy 'inne' pod pikiem) przy WP siatki.
Sygnał = szablon histogramowy z MC (klasa truth 0/1), tło gładkie. Zapis: out/d0_dstar.json.  Użycie: python s04_d0_dstar.py"""
import json
import os
import sys

import numpy as np
from iminuit import Minuit

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
CW = os.path.dirname(os.path.dirname(TOP))
IN = os.path.join(os.path.dirname(TOP), "runs", "lxplus")
sys.path.insert(0, os.path.join(CW, "DLZ_dplus", "SCAN_2026-11-17", "K1c", "scripts"))
import fastfit as FF  # noqa: E402

RUNS = {"427929": 1, "428770": 0}
MPI = 0.13957


def tfit(h, tmpl, x, order, kind="exp", t0=None, scale=0.1):
    """h = N*tmpl + nb*bkg(x); bkg: exp(pol_order(t)) albo próg u^a * exp(pol(u)). Zwraca N, err, chi2, ndf, valid, B pod pikiem (rdzeń 90 % szablonu)."""
    h = np.asarray(h, float)
    t = (x - (t0 if t0 is not None else x.mean())) / scale
    nm = ["N", "nb"] + [f"b{i+1}" for i in range(order)] + (["a"] if kind == "thr" else [])

    def bk(p):
        ex = sum(p[i] * t ** (i + 1) for i in range(order))
        f = np.exp(ex - ex.max())
        if kind == "thr":
            f = f * np.clip(x - MPI, 1e-9, None) ** p[order]
        return f / f.sum()

    def lam(p):
        return np.clip(p[0] * tmpl + p[1] * bk(p[2:]), 1e-9, None)

    def nll(*p):
        la = lam(p)
        return float(np.sum(la - h * np.log(la)))
    st = [0.0, h.sum()] + [0.0] * order + ([0.5] if kind == "thr" else [])
    m = Minuit(nll, *st, name=nm)
    m.errordef = 0.5
    m.limits["nb"] = (0, None)
    if kind == "thr":
        m.limits["a"] = (0, 5)
    m.migrad(ncall=40000)
    if not m.valid:
        m.migrad(ncall=40000)
    m.hesse()
    p = [m.values[n] for n in nm]
    la = lam(p)
    core = tmpl >= np.sort(tmpl)[::-1][np.searchsorted(np.cumsum(np.sort(tmpl)[::-1]), 0.9)]
    B = float((p[1] * bk(p[2:]))[core].sum())
    return dict(N=float(p[0]), err=float(m.errors["N"]), Z=float(p[0] / m.errors["N"]), chi2=float(((h - la) ** 2 / np.maximum(la, 1)).sum()), ndf=int(len(h) - m.nfit),
                valid=bool(m.valid), B_core=B, S_core=float(p[0] * tmpl[core].sum()), n_core_bins=int(core.sum()))


def main():
    out = {"d0": {}, "ds": {}, "dp_mc_closure": {}}
    hm = np.load(os.path.join(IN, "hist_fd_nd.npz"))
    # ---------------- D0: 5 MeV, 1,70–2,05
    e0 = hm["m0_edges"]
    x0 = 0.5 * (e0[1:] + e0[:-1]).reshape(70, 5).mean(1)
    t0 = hm["d0_h"][1].sum(0).reshape(70, 5).sum(1).astype(float)
    refl = hm["d0_h"][3].sum(0).reshape(70, 5).sum(1).astype(float)
    out["d0"]["mc"] = dict(n_sig=float(t0.sum()), n_inne=float(refl.sum()), rms_sig=float(np.sqrt(np.average((x0 - np.average(x0, weights=t0)) ** 2, weights=t0))),
                           peak=float(x0[np.argmax(t0)]))
    tm0 = t0 / t0.sum()
    print("D0 MC:", out["d0"]["mc"])
    # ---------------- D*: Δm 0,25 MeV, 0,140–0,165; |m(D0) − 1,865| < 40 MeV
    md, dm = hm["md_edges"], hm["dm_edges"]
    cmd = 0.5 * (md[1:] + md[:-1])
    selmd = np.abs(cmd - 1.865) < 0.040
    ndm = int(np.searchsorted(dm, 0.165 - 1e-9))
    xd = 0.5 * (dm[1:] + dm[:-1])[:ndm]
    ts = hm["ds_h"][0, 1].sum(0)[selmd].sum(0)[:ndm].astype(float)
    out["ds"]["mc"] = dict(n_sig=float(hm["ds_h"][0, 1].sum()), n_sig_sel=float(ts.sum()), peak=float(xd[np.argmax(ts)]), md_bins=[float(cmd[selmd].min()), float(cmd[selmd].max())])
    tms = ts / ts.sum()
    print("D* MC:", out["ds"]["mc"])
    for run, cls in RUNS.items():
        hd = np.load(os.path.join(IN, f"hist_data_{run}.npz"))
        out["d0"][run], out["ds"][run] = {}, {}
        for lab, sl in (("2-3", slice(0, 2)), ("2-2.5", slice(0, 1)), ("2.5-3", slice(1, 2))):
            h = hd["d0_h"][cls][sl].sum(0).reshape(70, 5).sum(1).astype(float)
            tm = hm["d0_h"][1][sl].sum(0).reshape(70, 5).sum(1).astype(float)
            tm = tm / tm.sum()
            r = {f"pol{o}": tfit(h, tm, x0, o, t0=1.865) for o in (3, 4, 5)}
            out["d0"][run][lab] = dict(n=float(h.sum()), **r)
            print(f"D0 {run} {lab}: n {h.sum():.0f}; " + "; ".join(f"pol{o}: N {r[f'pol{o}']['N']:.0f} ± {r[f'pol{o}']['err']:.0f} Z {r[f'pol{o}']['Z']:.1f} chi2 {r[f'pol{o}']['chi2']:.0f}/{r[f'pol{o}']['ndf']} S/B {r[f'pol{o}']['S_core']/r[f'pol{o}']['B_core']:.5f}"
                                                                 for o in (3, 4, 5)))
            hs = hd["ds_h"][0, cls][sl].sum(0)[selmd].sum(0)[:ndm].astype(float)
            tsl = hm["ds_h"][0, 1][sl].sum(0)[selmd].sum(0)[:ndm].astype(float)
            tsl = tsl / tsl.sum()
            r = {f"thr_pol{o}": tfit(hs, tsl, xd, o, kind="thr", t0=0.150, scale=0.01) for o in (1, 2, 3)}
            out["ds"][run][lab] = dict(n=float(hs.sum()), **r)
            print(f"D* {run} {lab}: n {hs.sum():.0f}; " + "; ".join(f"{k}: N {v['N']:.0f} ± {v['err']:.0f} Z {v['Z']:.1f} chi2 {v['chi2']:.0f}/{v['ndf']} S/B {v['S_core']/max(v['B_core'],1):.4f}" for k, v in r.items()))
    # ---------------- domknięcie D+ na fd_nd
    WP = json.load(open(os.path.join(TOP, "out", "wp.json")))
    H = hm["dp_h"]
    for w in WP["wp"]["grid"]:
        ie = w["ie"]
        hall = H[0][:, ie:].sum((0, 1)).reshape(100, 5).sum(1).astype(float)
        nt = float(H[1][:, ie:].sum())
        oth = H[3][:, ie:].sum((0, 1)).reshape(100, 5).sum(1).astype(float)
        bk = H[2][:, ie:].sum((0, 1)).reshape(100, 5).sum(1).astype(float)
        win = np.abs(FF.C - FF.M0) < FF.WIN
        r = FF.fit(hall, w["shape"], order=3, free_shape=False)
        rf = FF.fit(hall, w["shape"], order=3, free_shape=True)
        side = (np.abs(FF.C - FF.M0) > 0.07)
        out["dp_mc_closure"][str(w["wp"])] = dict(N_true=nt, N_fit=r["N"], N_err=r["N_err"], ratio=r["N"] / nt, ratio_err=r["N_err"] / nt, chi2=r["chi2"], ndf=r["ndf"],
                                                  N_free=rf["N"], dmu_free=rf["dmu"], sc_free=rf["sc"], inne_win=float(oth[win].sum()), tlo_win=float(bk[win].sum()),
                                                  inne_dens_win_over_side=float(oth[win].mean() / max(oth[side].mean(), 1e-9)), SB=r["SB"])
        print(f"D+ MC fd_nd WP{int(w['wp']*100)}: prawda {nt:.0f}, fit {r['N']:.0f} ± {r['N_err']:.0f} (iloraz {r['N']/nt:.3f} ± {r['N_err']/nt:.3f}; chi2 {r['chi2']:.0f}/{r['ndf']}); "
              f"'inne' w oknie {oth[win].sum():.0f} (gęstość okno/pasma {out['dp_mc_closure'][str(w['wp'])]['inne_dens_win_over_side']:.2f}), tło≥8 w oknie {bk[win].sum():.0f}; wolny: dmu {rf['dmu']*1e3:+.2f} sc {rf['sc']:.3f}")
    json.dump(out, open(os.path.join(TOP, "out", "d0_dstar.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

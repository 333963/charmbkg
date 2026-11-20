#!/usr/bin/env python3
"""XSEC_final (20.11.2026, projekt B) — rozstrzygniecie rozbieznosci kappa i koncowe uzyski D0 (M8, r = 1). DEFINICJA.md §2, §3. POST HOC.
Importuje bez zmian: FITBIAS fb_fit/fb_var, BKGCAT bc_m8 (patch), XSEC_atlas xa_fit. Dane: tylko inkluzywne piki 427929 (bez AFP).
  python -B kappa_fit.py A   -> out/kappa_A.json  (histogramy WPSCAN: warianty szablonu / r / modelu, bootstrap, domkniecie)
  python -B kappa_fit.py B   -> out/kappa_B.json  (tabele XSEC_atlas: skan kappa, punkty startu, ciecie eta, szablon BKGCAT)
  python -B kappa_fit.py C   -> out/kappa_altmin.json  (POST HOC: wielokrotne minima fitu nominalnego w skali szerokosci)"""
import json
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
CW = os.path.dirname(BASE)
sys.path.insert(0, os.path.join(CW, "BKGCAT_2026-11-20", "skrypty"))
sys.path.insert(0, os.path.join(CW, "XSEC_atlas_2026-11-20", "skrypty"))
import bc_m8 as M  # noqa: E402
V, F, L = M.V, M.F, M.L
BINS = F.BINS
SP = "D0"
RFB = 1.364


def fs(f):
    return {"N": f["N"][0], "err": (f.get("err") or [None])[0], "valid": f["valid"], "at_limit": f["at_limit"], "sc": f["sc"], "dmu": f["dmu"], "chi2": f["chi2"], "ndf": f["ndf"],
            "p": F.pval(f["chi2"], f["ndf"]), "fval": f["fval"]}


def part_A():
    M._init(L.IN)
    out = {}
    for bn in BINS:
        k0 = V.wcat(SP, bn, 0.7)
        D = V.wdata(k0, V.SRC["427929"])                       # [nrep+1, nfull]
        o = {"n_rep": int(D.shape[0] - 1), "kap_ref": k0["kap"]}
        kd = M.patch(k0, SP, bn, 0.7, None, "dir", 1.0)
        o["kap_dir"] = kd["kap"]
        for mode in ("ref", "dir"):
            for r in (RFB, 1.0):
                k = M.patch(k0, SP, bn, 0.7, None, mode, r)
                for m in ("M8", "M4", "M9", "M0"):
                    f = V.fit_set(D[0][None], [k], SP, m)
                    o[f"{m}|{mode}|{r}"] = fs(f)
        # nominal: M8 dir r=1 — bootstrap, normalizacja kappa
        k = M.patch(k0, SP, bn, 0.7, None, "dir", 1.0)
        f0 = V.fit_set(D[0][None], [k], SP, "M8")
        bN = []
        for i in range(1, D.shape[0]):
            fb = V.fit_set(D[i][None], [k], SP, "M8", start=f0["x"], hesse=False)
            if fb["valid"]:
                bN.append(fb["N"][0])
        o["boot"] = {"n_ok": len(bN), "rms": float(np.std(bN, ddof=1)), "mean": float(np.mean(bN))}
        for tag, s in (("k+5.8", 1.058), ("k-5.8", 0.942), ("k+20", 1.2), ("k-20", 0.8), ("k0", 0.0)):
            o[f"M8|dir|{tag}"] = fs(V.fit_set(D[0][None], [M.patch(k0, SP, bn, 0.7, None, "dir", s)], SP, "M8"))
        # skan r dla dir i ref (porownanie z B)
        o["scan_dir"] = {f"{r:.3f}": fs(V.fit_set(D[0][None], [M.patch(k0, SP, bn, 0.7, None, "dir", r)], SP, "M8")) for r in np.linspace(1.0, RFB, 5)}
        # domkniecie mb_nd
        H = V.wdata(k0, "mb|all")
        nt = float(V.wdata(k0, "mb|sig").sum())
        o["N_truth_mb"] = nt
        for m in ("M8", "M4", "M9"):
            f = V.fit_set(H[None], [M.patch(k0, SP, bn, 0.7, None, "dir", 1.0)], SP, m)
            o[f"closure|{m}|dir"] = {"b": f["N"][0] / nt - 1, "b_err": (f.get("err") or [float("nan")])[0] / nt, "valid": f["valid"], "p": F.pval(f["chi2"], f["ndf"])}
        out[bn] = o
        n = o["M8|dir|1.0"]
        print(f"A {bn}: M8 dir r=1 N = {n['N']:.1f} ± {n['err']:.1f} (boot {o['boot']['rms']:.1f}, n={o['boot']['n_ok']}) p = {n['p']:.3f}; ref r=1 {o['M8|ref|1.0']['N']:.1f}; ref rFB {o[f'M8|ref|{RFB}']['N']:.1f}; "
              f"dir rFB {o[f'M8|dir|{RFB}']['N']:.1f}; M4 {o['M4|dir|1.0']['N']:.1f}; M9 {o['M9|dir|1.0']['N']:.1f}; k±5.8 {o['M8|dir|k+5.8']['N']:.1f}/{o['M8|dir|k-5.8']['N']:.1f}; "
              f"domkn. {100 * o['closure|M8|dir']['b']:+.1f} ± {100 * o['closure|M8|dir']['b_err']:.1f} %", flush=True)
    json.dump(out, open(os.path.join(BASE, "out", "kappa_A.json"), "w"), indent=1, default=float)


def part_B():
    import xa_fit as X
    M._init(L.IN)
    D = X.load("data_427929", SP); D = {k: v[(D["cls"] & 3) == 3] for k, v in D.items()}
    MB, FD = X.load("mb_nd", SP), X.load("fd_nd", SP)
    out = {}
    for bn, lo, hi in X.X.BINS:
        o = {}
        for etac in (True, False):
            inb = lambda T: (T["pt"] >= lo) & (T["pt"] < hi) & ((np.abs(T["eta"]) < 2.2) if etac else True)  # noqa: E731
            d, m, f = inb(D) & X.sel_wp(D, SP, bn), inb(MB) & X.sel_wp(MB, SP, bn), inb(FD) & X.sel_wp(FD, SP, bn)
            comp = X.components(SP, FD, MB, f, m, bn)
            h = X.mh(D["m"][d], SP)
            e = {"n_cand": int(d.sum()), "kap": comp["kap"], "N_truth": comp["N_truth"], "frefl": comp["frefl"], "scan": {}}
            prev = None
            for r in np.linspace(1.0, RFB, 5):
                f0 = fs(X.dofit(h, comp, SP, kscale=float(r)))                       # start domyslny (jak XSEC_atlas)
                ft = X.dofit(h, comp, SP, kscale=float(r), start=prev) if prev is not None else None
                e["scan"][f"{r:.3f}"] = {"start0": f0, "start_prev": fs(ft) if ft else None}
                prev = (ft or X.dofit(h, comp, SP, kscale=float(r)))["x"]
            # punkt startu: rozwiazanie r=1 uzyte dla r=1,364 i odwrotnie
            f1 = X.dofit(h, comp, SP, kscale=1.0); f2 = X.dofit(h, comp, SP, kscale=RFB)
            e["r1364_start_from_r1"] = fs(X.dofit(h, comp, SP, kscale=RFB, start=f1["x"]))
            e["r1_start_from_r1364"] = fs(X.dofit(h, comp, SP, kscale=1.0, start=f2["x"]))
            # ten sam histogram, szablon i kappa BKGCAT (dir) / FITBIAS (ref)
            k0 = V.wcat(SP, bn, 0.7)
            for mode in ("dir", "ref"):
                for r in (1.0, RFB):
                    k = M.patch(k0, SP, bn, 0.7, None, mode, r)
                    c2 = dict(comp); c2["peakT"], c2["kap"] = k["peakT"], k["kap"]
                    e[f"tmpl_{mode}|{r}"] = fs(X.dofit(h, c2, SP))
            # histogram WPSCAN z komponentami XSEC_atlas (tylko bez ciecia eta ma sens)
            if not etac:
                Dw = V.wdata(k0, V.SRC["427929"])[0]
                hw = V.cutr(Dw, SP, F.SP[SP]["alt"])
                e["n_wps"] = float(hw.sum()); e["n_xa"] = float(h.sum()); e["max_abs_diff_hist"] = float(np.abs(hw - h).max())
                for r in (1.0, RFB):
                    e[f"wps_hist_xa_comp|{r}"] = fs(X.dofit(hw, comp, SP, kscale=r))
            o["eta22" if etac else "bez_eta"] = e
            s = e["scan"]
            print(f"B {bn} {'|eta|<2.2' if etac else 'bez eta'}: kand {e['n_cand']} kap {e['kap']:.3f}  " + "  ".join(f"r={r}: {v['start0']['N']:.1f} (sc {v['start0']['sc']:.3f}, dmu {1e3 * v['start0']['dmu']:+.2f}, fval {v['start0']['fval']:.2f}"
                  + (f"; prev {v['start_prev']['N']:.1f}, fval {v['start_prev']['fval']:.2f}" if v['start_prev'] else "") + ")" for r, v in s.items()), flush=True)
            print(f"     r=1,364 start z r=1: {e['r1364_start_from_r1']['N']:.1f} fval {e['r1364_start_from_r1']['fval']:.2f};  r=1 start z r=1,364: {e['r1_start_from_r1364']['N']:.1f} fval {e['r1_start_from_r1364']['fval']:.2f};  "
                  + "  ".join(f"{k}: {e[k]['N']:.1f}" for k in e if k.startswith("tmpl_") or k.startswith("wps_hist")), flush=True)
            if not etac:
                print(f"     hist WPSCAN {e['n_wps']:.0f} vs tabele {e['n_xa']:.0f}, max |roznica| w binie {e['max_abs_diff_hist']:.0f}")
        out[bn] = o
    json.dump(out, open(os.path.join(BASE, "out", "kappa_B.json"), "w"), indent=1, default=float)


def part_C():
    """POST HOC (po obejrzeniu part_B): minima fitu nominalnego (M8, dir, r = 1) z roznych punktow startu skali szerokosci. Alternatywa = wazne minimum z |dN/N| > 1 % i d(-2lnL) < 1."""
    M._init(L.IN)
    out = {}
    for bn in BINS:
        k0 = V.wcat(SP, bn, 0.7)
        D = V.wdata(k0, V.SRC["427929"])[0]
        k = M.patch(k0, SP, bn, 0.7, None, "dir", 1.0)
        f0 = V.fit_set(D[None], [k], SP, "M8")
        rows, alt = [], None
        for sc0 in (0.80, 0.85, 0.90, 0.945, 1.0, 1.05, 1.10, 1.20):
            x = list(f0["x"]); x[1] = sc0
            f = V.fit_set(D[None], [k], SP, "M8", start=x)
            r = dict(fs(f), sc0=sc0, dfval=f["fval"] - f0["fval"])
            rows.append(r)
            if f["valid"] and not f["at_limit"] and abs(f["N"][0] / f0["N"][0] - 1) > 0.01 and r["dfval"] < 1.0 and (alt is None or abs(r["N"] - f0["N"][0]) > abs(alt["N"] - f0["N"][0])):
                alt = r
        out[bn] = {"nom": fs(f0), "starts": rows, "altmin": alt}
        print(f"C {bn}: nom {f0['N'][0]:.1f} (sc {f0['sc']:.3f}); " + ("alt " + f"{alt['N']:.1f} (sc {alt['sc']:.3f}, d(-2lnL) {alt['dfval']:+.2f}, err {alt['err']})" if alt else "brak alternatywy"), flush=True)
    json.dump(out, open(os.path.join(BASE, "out", "kappa_altmin.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    {"A": part_A, "B": part_B, "C": part_C}[sys.argv[1]]()

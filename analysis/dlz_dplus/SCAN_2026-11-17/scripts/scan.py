"""SCAN_2026-11-17: R1 (skan cięcia N/ε vs WP) i R2 (pseudo-dane przy S/B danych) dla H i H+DLZ, D+ -> Kpipi.
Wejście: wyłącznie DATA_2026-11-17/runs/lxplus/ana/{res_RUN.json, rep_RUN.npz} (tylko odczyt). KRYTERIA.md obok.
Użycie: python scan.py [--ntoy 100] [--nproc 6]
"""
import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from iminuit import Minuit

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
DATA = os.path.join(os.path.dirname(OUT), "DATA_2026-11-17")
sys.path.insert(0, os.path.join(DATA, "scripts"))
import fitlib as FL  # noqa: E402

RUNS = ["427929", "428770"]
BINS = ["3-5", "5-8", "8-30"]
TAGS = {"H": "H_s12345", "DLZ": "HDLZ_s12345"}
WPS = [0.5, 0.7, 0.8, 0.9]          # skan (kryteria)
WPX = [0.5, 0.7, 0.8, 0.9, 0.95]    # + WP95 opisowo
C = FL.C
SB45 = np.abs(C - FL.M0) > 0.045
SB70 = np.abs(C - FL.M0) > 0.070


def key(t, b, wp):
    return f"{TAGS[t]}|{b}|{wp}"


def fit_bkg_only(h, mask, order, sig=None):
    """exp(pol-order) dopasowane do binów mask (Poisson ML); zwraca oczekiwane tło we wszystkich binach."""
    names = ["nb"] + [f"b{i+1}" for i in range(order)]

    def lam_of(p):
        return np.clip(p[0] * np.diff(FL.bkg_cdf(FL.E, p[1:])), 1e-9, None)

    def nll(*p):
        lam = lam_of(p) + (0.0 if sig is None else sig)   # sig: oczekiwane ogony sygnału (POST HOC, --tailsub), stałe
        return float(np.sum((lam - h * np.log(lam))[mask]))
    m = Minuit(nll, *([h.sum(), -0.5] + [0.0] * (order - 1)), name=names)
    m.errordef = 0.5
    m.limits["nb"] = (0, None)
    m.migrad(ncall=40000)
    if not m.valid:
        m.migrad(ncall=40000)
    lam = lam_of(list(m.values))
    tot = lam + (0.0 if sig is None else sig)
    pear = ((h - tot) ** 2 / np.maximum(tot, 1))[mask]
    return lam, bool(m.valid), float(pear.sum() / (mask.sum() - order - 1))


def toy_cell(a):
    """jedna komórka R2: Asimov + ntoy pseudo-eksperymentów, fit swobodny i o kształcie ustalonym."""
    (lab, sh, ntrue, bkg, ntoy, seed) = a
    sig = ntrue * np.diff(FL.sig_cdf(FL.E, sh, 0.0, 1.0))
    exp_ = sig + bkg
    rng = np.random.default_rng(seed)
    res = {"lab": lab, "N_true": float(ntrue), "B_win": float(bkg[~SB45].sum()), "SB_true": float(sig[~SB45].sum() / max(bkg[~SB45].sum(), 1e-9))}
    toys = rng.poisson(exp_, (ntoy, len(exp_))).astype(float)
    for fv, free in (("free", True), ("fix", False)):
        st0 = dict(N=ntrue, dmu=0.0, sc=1.0, nb=bkg.sum(), b1=-0.5, b2=0.0)
        ra = FL.fit(exp_, sh, free_shape=free, start=st0)
        start = [ra["N"], ra["dmu"] if free else 0.0, ra["sc"] if free else 1.0, ra["nb"]] + ra["b"]
        r = np.array([FL.fit_fast(toys[i], sh, start, free_shape=free) for i in range(ntoy)])
        ok = (r[:, 4] > 0) & np.isfinite(r[:, 0]) & (r[:, 1] > 0)
        x = r[ok, 0] / ntrue
        pull = (r[ok, 0] - ntrue) / r[ok, 1]
        res[fv] = dict(asimov=float(ra["N"] / ntrue), asimov_valid=bool(ra["valid"]), asimov_at_limit=bool(ra["at_limit"]), asimov_chi2=float(ra["chi2_ndf"] * (len(exp_) - (6 if free else 4))),
                       n_ok=int(ok.sum()), f_ok=float(ok.mean()),
                       b=float(x.mean()) if ok.any() else float("nan"), b_err=float(x.std(ddof=1) / np.sqrt(ok.sum())) if ok.sum() > 1 else float("nan"),
                       med=float(np.median(x)) if ok.any() else float("nan"), rms=float(x.std(ddof=1)) if ok.sum() > 1 else float("nan"),
                       pull_mean=float(pull.mean()) if ok.any() else float("nan"), pull_rms=float(pull.std(ddof=1)) if ok.sum() > 1 else float("nan"),
                       relerr_hesse=float(np.mean(r[ok, 1]) / ntrue) if ok.any() else float("nan"))
    return res


def gls_const(y, Cv):
    Ci = np.linalg.inv(Cv)
    one = np.ones(len(y))
    v = 1.0 / (one @ Ci @ one)
    c = v * (one @ Ci @ y)
    d = y - c
    return float(c), float(np.sqrt(v)), float(d @ Ci @ d)


def gls_slope(y, Cv, x):
    """y = a + g·(x − 0,7); zwraca s = g/a i błąd (propagacja)."""
    A = np.c_[np.ones(len(y)), np.asarray(x) - 0.7]
    Ci = np.linalg.inv(Cv)
    cov = np.linalg.inv(A.T @ Ci @ A)
    be = cov @ (A.T @ Ci @ y)
    s = be[1] / be[0]
    J = np.array([-be[1] / be[0] ** 2, 1 / be[0]])
    d = y - A @ be
    return float(s), float(np.sqrt(J @ cov @ J)), float(d @ Ci @ d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ntoy", type=int, default=100)
    ap.add_argument("--nproc", type=int, default=6)
    ap.add_argument("--tailsub", action="store_true", help="POST HOC: tło prawdy A/C dopasowane do SB z odjętymi oczekiwanymi ogonami sygnału S0·ε·kształt MC; zapis z sufiksem _posthoc")
    a = ap.parse_args()
    out = {"ntoy": a.ntoy, "wps": WPS, "cells": {}}
    jobs = []
    R = {}
    for run in RUNS:
        res = json.load(open(os.path.join(DATA, "runs", "lxplus", "ana", f"res_{run}.json")))
        rep = np.load(os.path.join(DATA, "runs", "lxplus", "ana", f"rep_{run}.npz"))
        R[run] = (res, rep)
        for b in BINS:
            lim = {t: {str(wp): bool(res["at_limit"][key(t, b, wp)]) for wp in WPX} for t in TAGS}
            fv = "fix" if any(lim[t][str(wp)] for t in TAGS for wp in WPS) else "free"
            cell = {"at_limit_free": lim, "fitvar": fv, "R1": {}, "R2": {}}
            out["cells"][f"{run}|{b}"] = cell
            for v in ("free", "fix"):
                F = res["fits"] if v == "free" else res["fits_fixshape"]
                pre = "" if v == "free" else "fix|"
                # macierz replik: kolumny (tag, wp) po WPX
                cols = [(t, wp) for t in TAGS for wp in WPX]
                Nrep = np.column_stack([rep[pre + key(t, b, wp)][:, 0] for t, wp in cols])
                okr = np.column_stack([(rep[pre + key(t, b, wp)][:, 4] > 0) & (rep[pre + key(t, b, wp)][:, 0] > 0) for t, wp in cols])
                eps = np.array([res["shapes"][key(t, b, wp)]["eps_mc"] for t, wp in cols])
                nmc = np.array([res["shapes"][key(t, b, wp)]["n_mc"] for t, wp in cols])
                N = np.array([F[key(t, b, wp)]["N"] for t, wp in cols])
                Ne = np.array([F[key(t, b, wp)]["N_err"] for t, wp in cols])
                d = dict(cols=[f"{t}|{wp}" for t, wp in cols], N=N.tolist(), N_err_hesse=Ne.tolist(), eps=eps.tolist(), y=(N / eps).tolist(),
                         chi2_ndf=[F[key(t, b, wp)]["chi2_ndf"] for t, wp in cols], valid=[bool(F[key(t, b, wp)]["valid"]) for t, wp in cols],
                         SB=[F[key(t, b, wp)]["SB"] for t, wp in cols], Z=[F[key(t, b, wp)]["Z"] for t, wp in cols],
                         boot_ok_frac=okr.mean(0).tolist(), N_err_boot=[float(np.std(Nrep[okr[:, i], i])) for i in range(len(cols))])
                cell["R1"][v] = d
                cell["R1"][v]["_Nrep"] = Nrep
                cell["R1"][v]["_ok"] = okr
            # R2: prawda
            for v in ("free", "fix"):
                pass
            Fc = res["fits"] if fv == "free" else res["fits_fixshape"]
            for vv in ("free", "fix"):
                Fv = res["fits"] if vv == "free" else res["fits_fixshape"]
                cell.setdefault("S0", {})[vv] = Fv[key("H", b, 0.7)]["N"] / res["shapes"][key("H", b, 0.7)]["eps_mc"]
            S0 = cell["S0"][fv]
            for t in TAGS:
                for wp in WPX:
                    k = key(t, b, wp)
                    h = np.asarray(res["hists"][k], float)
                    sh = res["shapes"][k]
                    sg = S0 * sh["eps_mc"] * np.diff(FL.sig_cdf(FL.E, sh, 0.0, 1.0)) if a.tailsub else None
                    bA, vA, cA = fit_bkg_only(h, SB45, 4, sg)
                    bC, vC, cC = fit_bkg_only(h, SB70, 4, sg)
                    f0 = Fc[k]
                    bP = f0["nb"] * np.diff(FL.bkg_cdf(FL.E, f0["b"]))
                    cell["R2"][f"{t}|{wp}|truthbkg"] = dict(A_valid=vA, A_chi2sb=cA, C_valid=vC, C_chi2sb=cC, B_win_A=float(bA[~SB45].sum()), B_win_C=float(bC[~SB45].sum()),
                                                            B_win_fit=float(bP[~SB45].sum()), pol2_chi2sb_data=float(f0["chi2_sb"] / (f0["n_sb"] - 3)))
                    for bv, bk in (("A", bA), ("C", bC), ("P", bP)):
                        jobs.append((f"{run}|{b}|{t}|{wp}|{bv}", sh, S0 * sh["eps_mc"], bk, a.ntoy, 20261117 + len(jobs)))
    print("komórek R2:", len(jobs), flush=True)
    with Pool(a.nproc) as p:
        for i, r in enumerate(p.imap_unordered(toy_cell, jobs, chunksize=1)):
            run, b, t, wp, bv = r["lab"].split("|")
            out["cells"][f"{run}|{b}"]["R2"][f"{t}|{wp}|{bv}"] = r
            if i % 20 == 0:
                print("toy", i, r["lab"], flush=True)
    # --- pochodne + kryteria
    for ck, cell in out["cells"].items():
        iw = [WPX.index(w) for w in WPS]
        for v in ("free", "fix"):
            d = cell["R1"][v]
            Nrep, okr = d.pop("_Nrep"), d.pop("_ok")
            eps, N = np.array(d["eps"]), np.array(d["N"])
            nW = len(WPX)
            sel = [i for i in iw] + [nW + i for i in iw]
            ok = okr[:, sel].all(1)
            d["n_rep_ok"] = int(ok.sum())
            Y = Nrep[ok][:, sel] / eps[sel]
            y = N[sel] / eps[sel]
            Cy = np.cov(Y.T)
            n4 = len(iw)
            der = {}
            for ti, t in enumerate(TAGS):
                sl = slice(ti * n4, (ti + 1) * n4)
                c, ce, chi = gls_const(y[sl], Cy[sl, sl])
                s, se, chil = gls_slope(y[sl], Cy[sl, sl], WPS)
                der[t] = dict(y=y[sl].tolist(), y_err=np.sqrt(np.diag(Cy[sl, sl])).tolist(), mean=c, mean_err=ce, chi2_flat=chi, slope=s, slope_err=se, chi2_lin=chil,
                              corr=np.corrcoef(Y[:, sl].T).tolist())
            Rrep = Y[:, n4:] / Y[:, :n4]
            Rv = y[n4:] / y[:n4]
            CR = np.cov(Rrep.T)
            c, ce, chi = gls_const(Rv, CR)
            der["R"] = dict(R=Rv.tolist(), R_err=np.sqrt(np.diag(CR)).tolist(), mean=c, mean_err=ce, chi2_const=chi)
            # kontrola błędu różnicy zagnieżdżonych: WP90 − WP70
            for ti, t in enumerate(TAGS):
                i70, i90 = ti * n4 + WPS.index(0.7), ti * n4 + WPS.index(0.9)
                sa, sb_ = np.std(Nrep[ok][:, sel][:, i90]), np.std(Nrep[ok][:, sel][:, i70])
                der[t]["nested_check_90_70"] = dict(sig_diff_boot=float(np.std(Nrep[ok][:, sel][:, i90] - Nrep[ok][:, sel][:, i70])), sig_diff_approx=float(np.sqrt(abs(sa**2 - sb_**2))),
                                                    sig_diff_uncorr=float(np.hypot(sa, sb_)))
            # WP95 opisowo
            der["wp95"] = {t: dict(y=float(N[ti * nW + 4] / eps[ti * nW + 4]), y_err=float(np.std(Nrep[okr[:, ti * nW + 4], ti * nW + 4]) / eps[ti * nW + 4])) for ti, t in enumerate(TAGS)}
            # korekta R2 (ten sam wariant fitu)
            bA = {t: np.array([cell["R2"][f"{t}|{wp}|A"][v]["b"] for wp in WPS]) for t in TAGS}
            bAe = {t: np.array([cell["R2"][f"{t}|{wp}|A"][v]["b_err"] for wp in WPS]) for t in TAGS}
            bC = {t: np.array([cell["R2"][f"{t}|{wp}|C"][v]["b"] for wp in WPS]) for t in TAGS}
            cor = {}
            for ti, t in enumerate(TAGS):
                sl = slice(ti * n4, (ti + 1) * n4)
                yc = y[sl] / bA[t]
                Cc = Cy[sl, sl] / np.outer(bA[t], bA[t]) + np.diag(yc**2 * ((bAe[t] / bA[t]) ** 2 + ((bA[t] - bC[t]) / bA[t]) ** 2))
                c, ce, chi = gls_const(yc, Cc)
                s, se, chil = gls_slope(yc, Cc, WPS)
                cor[t] = dict(y=yc.tolist(), y_err=np.sqrt(np.diag(Cc)).tolist(), mean=c, mean_err=ce, chi2_flat=chi, slope=s, slope_err=se)
            rA, rC = bA["DLZ"] / bA["H"], bC["DLZ"] / bC["H"]
            rAe = rA * np.hypot(bAe["DLZ"] / bA["DLZ"], bAe["H"] / bA["H"])
            Rc = Rv / rA
            CRc = CR / np.outer(rA, rA) + np.diag(Rc**2 * ((rAe / rA) ** 2 + ((rA - rC) / rA) ** 2))
            c, ce, chi = gls_const(Rc, CRc)
            cor["R"] = dict(R=Rc.tolist(), R_err=np.sqrt(np.diag(CRc)).tolist(), mean=c, mean_err=ce, chi2_const=chi, rA=rA.tolist(), rC=rC.tolist())
            cell["R1"][v]["der"] = der
            cell["R1"][v]["corr"] = cor
        v = cell["fitvar"]
        i7 = WPS.index(0.7)
        cor, der = cell["R1"][v]["corr"], cell["R1"][v]["der"]
        b7 = {t: cell["R2"][f"{t}|0.7|A"][v] for t in TAGS}
        sens = all(abs(b7[t]["b"] - 1) < 0.05 and b7[t]["f_ok"] >= 0.8 for t in TAGS)
        r7, r7e = cor["R"]["R"][i7], cor["R"]["R_err"][i7]
        cell["crit"] = dict(c_sensitive=bool(sens), b70={t: [b7[t]["b"], b7[t]["b_err"], b7[t]["f_ok"]] for t in TAGS}, c_decisive=bool(sens and r7e <= 0.10),
                            a1_R70=[r7, r7e], a1_nsig=float((r7 - 1) / r7e), a1_pass=bool(abs(r7 - 1) < 2 * r7e),
                            a2_cell=[cor["R"]["mean"], cor["R"]["mean_err"]],
                            b_dchi2=float(cor["DLZ"]["chi2_flat"] - cor["H"]["chi2_flat"]), b_pass=bool(cor["DLZ"]["chi2_flat"] - cor["H"]["chi2_flat"] <= 4),
                            raw_R70=[der["R"]["R"][i7], der["R"]["R_err"][i7]], raw_dchi2=float(der["DLZ"]["chi2_flat"] - der["H"]["chi2_flat"]))
    verdict = {}
    for run in RUNS:
        cs = {b: out["cells"][f"{run}|{b}"]["crit"] for b in BINS}
        m = np.array([cs[b]["a2_cell"][0] for b in BINS])
        e = np.array([cs[b]["a2_cell"][1] for b in BINS])
        w = 1 / e**2
        mean, err = float((w * m).sum() / w.sum()), float(1 / np.sqrt(w.sum()))
        a2 = bool(abs(mean - 1) < err or abs(mean - 1) < 0.05)
        sensitive = [b for b in BINS if cs[b]["c_sensitive"]]
        decisive = [b for b in BINS if cs[b]["c_decisive"]]
        fail = [b for b in sensitive if not (cs[b]["a1_pass"] and cs[b]["b_pass"])]
        if fail or (not a2 and sensitive):
            vd = "FAIL"
        elif len(decisive) == len(BINS):
            vd = "PASS"
        else:
            vd = "NIEROZSTRZYGNIETE"
        verdict[run] = dict(verdict=vd, a2_mean=[mean, err], a2_pass=a2, sensitive=sensitive, decisive=decisive, fail_cells=fail,
                            fitvar={b: out["cells"][f"{run}|{b}"]["fitvar"] for b in BINS}, cells=cs)
    out["verdict"] = verdict
    os.makedirs(os.path.join(OUT, "out"), exist_ok=True)
    suf = "_posthoc" if a.tailsub else ""
    out["tailsub"] = bool(a.tailsub)
    json.dump(out, open(os.path.join(OUT, "out", f"scan{suf}.json"), "w"), indent=1)
    json.dump(verdict, open(os.path.join(OUT, f"verdict{suf}.json"), "w"), indent=1)
    print("zapisano", flush=True)


if __name__ == "__main__":
    main()

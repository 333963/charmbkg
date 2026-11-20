"""K1d: obciążenie fitu masy D+ na MC minbias z pełnym tłem (KRYTERIA_K1d.md).
Warianty: a (MC wprost), b1_R / b0_R (tło MC przeważone do pasm bocznych danych runu R z/bez zależności od masy; sygnał skalowany do N_H danych).
Użycie: python k1d.py [--nrep 200] [--ntoy 200] [--nproc 7]   → out/k1d.json, out/weights.json"""
import argparse
import json
import os
import sys
import time
import multiprocessing as _mp

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

HERE = os.path.dirname(os.path.abspath(__file__))
K1D = os.path.dirname(HERE)
SCAN = os.path.dirname(K1D)
sys.path.insert(0, os.path.join(SCAN, "K1c", "scripts"))
import prep as PR  # noqa: E402  (tylko odczyt: load_data, stałe)
FF = PR.FF
Pool = _mp.get_context("fork").Pool
RUNS, BINS, BN = PR.RUNS, PR.BINS, PR.BN
CL = ["H", "DLZ"]
WPS = [0.5, 0.7, 0.8]
ORD = [2, 3, 4]
NB = FF.NB
M0 = FF.M0
GROUPS = ["sig"] + [str(i) for i in range(2, 16)]
G = {}


def logit(s):
    s = np.clip(s.astype(np.float64), 1e-12, 1 - 1e-12)
    return np.clip(np.log(s / (1 - s)), -12, 12)


def feats(m, zh, zd, lp, mass=True):
    t = (m - M0) / 0.1
    c = [zh, zd, zh ** 2, zd ** 2, zh * zd, lp, lp * zh, lp * zd]
    if mass:
        c += [t, t ** 2, t * zh, t * zd]
    return np.stack(c, 1)


def fit_logit(Xd, wd, Xm):
    """min Σ_d w_d log(1+e^{-θx}) + Σ_mc log(1+e^{θx}); zwraca funkcję log-wagi."""
    mu = Xm.mean(0)
    sd = Xm.std(0) + 1e-12
    A = np.c_[np.ones(len(Xd)), (Xd - mu) / sd]
    Bm = np.c_[np.ones(len(Xm)), (Xm - mu) / sd]
    nt = len(Xd) + len(Xm)

    def f(th):
        a, b = A @ th, Bm @ th
        L = np.sum(wd * np.logaddexp(0, -a)) + np.sum(np.logaddexp(0, b))
        g = -(A.T @ (wd * expit(-a))) + Bm.T @ expit(b)
        return L / nt, g / nt
    th0 = np.zeros(A.shape[1])
    th0[0] = np.log(max(wd.sum(), 1) / len(Xm))
    r = minimize(f, th0, jac=True, method="L-BFGS-B", options=dict(maxiter=500, maxfun=2000, ftol=1e-12, gtol=1e-8))
    th = r.x
    return (lambda X: np.c_[np.ones(len(X)), (X - mu) / sd] @ th), dict(theta=th.tolist(), mu=mu.tolist(), sd=sd.tolist(), ok=bool(r.success), nit=int(r.nit))


def chain(h, sh, free, st0=None):
    out, prev = {}, st0
    for o in ORD:
        r = FF.fit(h, sh, order=o, free_shape=free, start=prev)
        out[o] = r
        prev = r["par"]
    return out


def job(a):
    var, bi, ci, wi, ntoy, seed = a
    bn, c, wp = BN[bi], CL[ci], WPS[wi]
    P = G["P"]
    free = G["fitvar"][bn] == "free"
    sh = P["shapes"][f"{c}|{bn}|{wp}"]
    D = G["hist"][var]
    h = D["nom"][bi, ci, wi]
    hg = D["grp"][bi, ci, wi]          # (ngroups, NB)
    nt = float(D["nt_nom"][bi, ci, wi])
    hb, ntb = D["boot"][:, bi, ci, wi, :], D["nt_boot"][:, bi, ci, wi]
    win = np.abs(FF.C - M0) < FF.WIN
    out = dict(lab=f"{var}|{bn}|{c}|{wp}", N_true=nt, n_tot=float(h.sum()), S_win=float(hg[0][win].sum()), B_win=float((h - hg[0])[win].sum()),
               sig_frac_in_win=float(hg[0][win].sum() / max(hg[0].sum(), 1e-9)), fit={})
    nom = chain(h, sh, free, st0=[nt, 0.0, 1.0, h.sum(), -0.5])
    rng = np.random.default_rng(seed)
    toys = rng.poisson(np.clip(h, 0, None), (ntoy, NB)).astype(float)
    hbkg = h - hg[0]
    for o in ORD:
        f = nom[o]
        q = dict(N=f["N"], N_err=f["N_err"], b=f["N"] / nt, rel_true=f["N_err"] / nt, chi2=f["chi2"], ndf=f["ndf"], valid=f["valid"], at_limit=f["at_limit"],
                 dmu=f["dmu"], sc=f["sc"], SB_fit=f["SB"])
        x = np.array([FF.fit_fast(hb[i], sh, f["par"], order=o, free_shape=free) for i in range(hb.shape[0])])
        ok = (x[:, 4] > 0) & np.isfinite(x[:, 0]) & (x[:, 1] > 0)
        q["boot"] = dict(N=x[:, 0].tolist(), err=x[:, 1].tolist(), ok=ok.astype(int).tolist(), nt=ntb.tolist())
        y = np.array([FF.fit_fast(toys[i], sh, f["par"], order=o, free_shape=free) for i in range(ntoy)])
        oky = (y[:, 4] > 0) & np.isfinite(y[:, 0]) & (y[:, 1] > 0)
        if oky.sum() > 2:
            pl = (y[oky, 0] - nt) / y[oky, 1]
            q["toy"] = dict(f_ok=float(oky.mean()), mean=float(y[oky, 0].mean() / nt), rms=float(y[oky, 0].std(ddof=1) / nt), relerr=float(y[oky, 1].mean() / nt),
                            pull_mean=float(pl.mean()), pull_w=float(pl.std(ddof=1)), pull_w_err=float(pl.std(ddof=1) / np.sqrt(2 * oky.sum())))
        else:
            q["toy"] = dict(f_ok=float(oky.mean()))
        # leave-one-out wg prawdy
        loo = {}
        for gi, g in enumerate(GROUPS):
            if hg[gi].sum() <= 0:
                loo[g] = dict(n=0.0, n_win=0.0, dN=0.0)
                continue
            r = FF.fit(np.clip(h - hg[gi], 0, None), sh, order=o, free_shape=free, start=f["par"])
            loo[g] = dict(n=float(hg[gi].sum()), n_win=float(hg[gi][win].sum()), dN=float(f["N"] - r["N"]), valid=r["valid"])
        pc = hg[1:7].sum(0)
        r = FF.fit(np.clip(h - pc, 0, None), sh, order=o, free_shape=free, start=f["par"])
        loo["pc_2-7"] = dict(n=float(pc.sum()), n_win=float(pc[win].sum()), dN=float(f["N"] - r["N"]), valid=r["valid"])
        st = list(f["par"])
        st[0] = 10.0
        r = FF.fit(np.clip(hbkg, 0, None), sh, order=o, free_shape=False, start=st[:3] + st[3:])
        loo["bkg_only_fixshape"] = dict(n=float(hbkg.sum()), n_win=float(hbkg[win].sum()), dN=float(r["N"]), valid=r["valid"])
        r = FF.fit(np.clip(hg[7:].sum(0), 0, None), sh, order=o, free_shape=False, start=st)
        loo["comb_only_fixshape"] = dict(n=float(hg[7:].sum()), n_win=float(hg[7:].sum(0)[win].sum()), dN=float(r["N"]), valid=r["valid"])
        q["loo"] = loo
        out["fit"][str(o)] = q
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nrep", type=int, default=200)
    ap.add_argument("--ntoy", type=int, default=200)
    ap.add_argument("--nproc", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()
    P = json.load(open(os.path.join(SCAN, "K1c", "out", "prep.json")))
    NOM = json.load(open(os.path.join(SCAN, "K1c", "out", "k1c_nominal.json")))
    G["P"] = P
    fitvar = {bn: NOM["fitvar"][f"427929|{bn}"] for bn in BN}
    for bn in BN:
        assert NOM["fitvar"][f"428770|{bn}"] == fitvar[bn]
    G["fitvar"] = fitvar
    oof = os.path.join(SCAN, "cand", "train", "oof")
    b = np.load(os.path.join(oof, "base.npz"))
    m, pt, tr, key = b["mb_m"].astype(float), b["mb_pt"].astype(float), b["mb_truth"], b["mb_key"]
    assert np.all(b["mb_w"] == 1)
    sc = {"H": np.load(os.path.join(oof, "H_s12345.npz"))["mb"], "DLZ": np.load(os.path.join(oof, "HDLZ_s12345.npz"))["mb"]}
    thr = lambda c, bn, wp: P["thr"][f"{c}|{bn}|{wp}"]  # noqa: E731
    ibin = np.digitize(pt, [3, 5, 8, 30]) - 1
    inr = (m >= FF.LO) & (m < FF.HI)
    ue, einv = np.unique(key, return_inverse=True)
    # --- kandydaci MC w obszarze U per bin
    MC = []
    for bi, bn in enumerate(BN):
        u = np.flatnonzero((ibin == bi) & inr & ((sc["H"] >= thr("H", bn, 0.8)) | (sc["DLZ"] >= thr("DLZ", bn, 0.8))))
        d = dict(m=m[u], lp=np.log(pt[u]), zh=logit(sc["H"][u]), zd=logit(sc["DLZ"][u]), tr=tr[u], ev=einv[u],
                 mb=np.clip(((m[u] - FF.LO) / 0.005).astype(int), 0, NB - 1), sig=tr[u] <= 1)
        d["pass"] = {(ci, wi): (sc[c][u] >= thr(c, bn, wp)) for ci, c in enumerate(CL) for wi, wp in enumerate(WPS)}
        d["grp"] = np.where(d["sig"], 0, d["tr"] - 1)        # 0 = sygnał, 1..14 = kody 2..15
        d["nsig_nocut"] = int(((ibin == bi) & inr & (tr <= 1)).sum())
        MC.append(d)
    # --- wagi wariantów
    W = {"a": [np.ones(len(d["m"])) for d in MC]}
    winfo = {}
    for run in RUNS:
        dd = PR.load_data(run)
        dib = np.digitize(dd["pt"], [3, 5, 8, 30]) - 1
        for vn, mass in (("b1", True), ("b0", False)):
            W[f"{vn}_{run}"] = []
        for bi, bn in enumerate(BN):
            d = MC[bi]
            nH = NOM["nominal"][f"{run}|{bn}|H|0.7"][fitvar[bn]]["fits"]["3"]["N"]
            k = nH / float((d["sig"] & d["pass"][(0, 1)]).sum())
            du = (dib == bi) & ((dd["s"]["H"] >= thr("H", bn, 0.8)) | (dd["s"]["DLZ"] >= thr("DLZ", bn, 0.8))) & (dd["m"] >= FF.LO) & (dd["m"] < FF.HI)
            dsb = du & (np.abs(dd["m"] - M0) > 0.045)
            msb = np.abs(d["m"] - M0) > 0.045
            sB, sS = msb & ~d["sig"], msb & d["sig"]
            for vn, mass in (("b1", True), ("b0", False)):
                Xd = feats(dd["m"][dsb], logit(dd["s"]["H"][dsb]), logit(dd["s"]["DLZ"][dsb]), np.log(dd["pt"][dsb].astype(float)), mass)
                Xs = feats(d["m"][sS], d["zh"][sS], d["zd"][sS], d["lp"][sS], mass)
                Xm = feats(d["m"][sB], d["zh"][sB], d["zd"][sB], d["lp"][sB], mass)
                fn, info = fit_logit(np.r_[Xd, Xs], np.r_[np.ones(len(Xd)), -k * np.ones(len(Xs))], Xm)
                lw = fn(feats(d["m"], d["zh"], d["zd"], d["lp"], mass))
                w = np.where(d["sig"], k, np.exp(lw))
                W[f"{vn}_{run}"].append(w)
                wb = w[~d["sig"]]
                info.update(k_sig=k, n_data_sb=int(dsb.sum()), n_mc_sb=int(sB.sum()), w_mean=float(wb.mean()), w_p1=float(np.percentile(wb, 1)), w_p99=float(np.percentile(wb, 99)),
                            w_max=float(wb.max()), neff_frac=float(wb.sum() ** 2 / (wb ** 2).sum() / len(wb)))
                # domknięcie: dane vs MC ważone po cięciu (pasma i okno)
                cl = {}
                for ci, c in enumerate(CL):
                    for wi, wp in enumerate(WPS):
                        ps = d["pass"][(ci, wi)]
                        hm = np.bincount(d["mb"][ps], weights=w[ps], minlength=NB)
                        hv = np.bincount(d["mb"][ps], weights=w[ps] ** 2, minlength=NB)
                        hd = np.asarray(P["hists"][f"{run}|{c}|{bn}|{wp}"], float)
                        sbm = np.abs(FF.C - M0) > 0.045
                        chi = float(np.sum(((hd - hm) ** 2 / np.maximum(hd + hv, 1))[sbm]) / sbm.sum())
                        fd = NOM["nominal"][f"{run}|{bn}|{c}|{wp}"][fitvar[bn]]["fits"]["3"]
                        cl[f"{c}|{wp}"] = dict(sb_ratio=float(hm[sbm].sum() / hd[sbm].sum()), chi2_sb=chi, win_ratio=float(hm[~sbm].sum() / hd[~sbm].sum()),
                                               win_data=float(hd[~sbm].sum()), win_mc=float(hm[~sbm].sum()), win_mc_err=float(np.sqrt(hv[~sbm].sum())), SB_data_fit=fd["SB"], N_data_fit=fd["N"])
                info["closure"] = cl
                winfo[f"{vn}_{run}|{bn}"] = info
                print(f"wagi {vn}_{run} {bn}: k={k:.2f} <w>={info['w_mean']:.2f} p99={info['w_p99']:.1f} max={info['w_max']:.0f} neff/n={info['neff_frac']:.2f} ok={info['ok']} "
                      f"SB70 H/DLZ ratio {cl['H|0.7']['sb_ratio']:.3f}/{cl['DLZ|0.7']['sb_ratio']:.3f} chi2 {cl['H|0.7']['chi2_sb']:.2f}/{cl['DLZ|0.7']['chi2_sb']:.2f} "
                      f"win {cl['H|0.7']['win_ratio']:.3f}/{cl['DLZ|0.7']['win_ratio']:.3f}  ({time.time()-t0:.0f} s)", flush=True)
        del dd
    json.dump(winfo, open(os.path.join(K1D, "out", "weights.json"), "w"), indent=1)
    # --- histogramy nominalne, grupy prawdy, bootstrap
    rng = np.random.default_rng(20261118)
    BW = rng.poisson(1.0, (a.nrep, len(ue))).astype(np.float32)
    G["hist"] = {}
    for var in W:
        D = dict(nom=np.zeros((3, 2, 3, NB)), grp=np.zeros((3, 2, 3, len(GROUPS), NB)), nt_nom=np.zeros((3, 2, 3)),
                 boot=np.zeros((a.nrep, 3, 2, 3, NB)), nt_boot=np.zeros((a.nrep, 3, 2, 3)))
        for bi in range(3):
            d, w = MC[bi], W[var][bi]
            for (ci, wi), ps in d["pass"].items():
                mb_, w_, g_, ev_, sg = d["mb"][ps], w[ps], d["grp"][ps], d["ev"][ps], d["sig"][ps]
                D["nom"][bi, ci, wi] = np.bincount(mb_, weights=w_, minlength=NB)
                D["grp"][bi, ci, wi] = np.bincount(g_ * NB + mb_, weights=w_, minlength=len(GROUPS) * NB).reshape(len(GROUPS), NB)
                D["nt_nom"][bi, ci, wi] = w_[sg].sum()
                for r in range(a.nrep):
                    ww = BW[r, ev_] * w_
                    D["boot"][r, bi, ci, wi] = np.bincount(mb_, weights=ww, minlength=NB)
                    D["nt_boot"][r, bi, ci, wi] = ww[sg].sum()
        G["hist"][var] = D
    print(f"histogramy: {time.time()-t0:.0f} s", flush=True)
    jobs = [(var, bi, ci, wi, a.ntoy, 555000 + i) for i, (var, bi, ci, wi) in
            enumerate((v, bi, ci, wi) for v in W for bi in range(3) for ci in range(2) for wi in range(3))]
    jobs.sort(key=lambda j: WPS[j[3]] != 0.7)
    out = dict(nrep=a.nrep, ntoy=a.ntoy, fitvar=fitvar, variants=list(W), groups=GROUPS, cells={},
               mc=dict(n_events=int(len(ue)), nsig_nocut=[d["nsig_nocut"] for d in MC],
                       eps_mb={f"{c}|{BN[bi]}|{wp}": float((MC[bi]["sig"] & MC[bi]["pass"][(ci, wi)]).sum() / MC[bi]["nsig_nocut"])
                               for bi in range(3) for ci, c in enumerate(CL) for wi, wp in enumerate(WPS)}))
    with Pool(a.nproc) as p:
        for i, r in enumerate(p.imap_unordered(job, jobs, chunksize=1)):
            out["cells"][r["lab"]] = r
            if i % 10 == 9:
                print(f"  {i+1}/{len(jobs)} komórek, {time.time()-t0:.0f} s", flush=True)
                json.dump(out, open(os.path.join(K1D, "out", "k1d_partial.json"), "w"))
    json.dump(out, open(os.path.join(K1D, "out", "k1d.json"), "w"))
    print(f"koniec: {time.time()-t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()

"""LOWPT_calib: skan cięcia H dla D+ (dane 427929/428770) + wstrzyknięcie na MC mb z pełnym tłem. KRYTERIA.md.
Silnik fitu i przygotowanie: DLZ_dplus/SCAN_2026-11-17/K1c/scripts (tylko odczyt/import).
Użycie: python scan.py [--nrep 200] [--nproc 7]  -> out/scan.json"""
import argparse
import json
import os
import sys
import time
import multiprocessing as _mp

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
SCAN = os.path.join(os.path.dirname(TOP), "DLZ_dplus", "SCAN_2026-11-17")
sys.path.insert(0, os.path.join(SCAN, "K1c", "scripts"))
import prep as PR  # noqa: E402
FF = PR.FF
Pool = _mp.get_context("fork").Pool
RUNS, BINS, BN = PR.RUNS, PR.BINS, PR.BN
WPS = [0.5, 0.6, 0.7, 0.8, 0.9]
ORD = [2, 3, 4]
BORD = [3, 4]
NB = FF.NB
G = {}


def job(a):
    src, bi, wi = a
    bn, wp = BN[bi], WPS[wi]
    sh = G["P"]["shapes"][f"H|{bn}|{wp}"]
    free = G["fitvar"][bn] == "free"
    D = G["hist"][src]
    h = D["nom"][bi, wi]
    out = dict(lab=f"{src}|{bn}|{wp}", n=float(h.sum()), fit={}, boot={})
    prev = None
    for o in ORD:
        r = FF.fit(h, sh, order=o, free_shape=free, start=prev)
        prev = r["par"]
        out["fit"][str(o)] = {k: r[k] for k in ("N", "N_err", "chi2", "ndf", "valid", "at_limit", "SB", "dmu", "sc")}
        if o in BORD:
            hb = D["boot"][:, bi, wi, :]
            x = np.array([FF.fit_fast(hb[i], sh, r["par"], order=o, free_shape=free) for i in range(hb.shape[0])])
            out["boot"][str(o)] = dict(N=x[:, 0].tolist(), err=x[:, 1].tolist(), ok=((x[:, 4] > 0) & np.isfinite(x[:, 0]) & (x[:, 1] > 0)).astype(int).tolist())
    return out


def hists(mb, ev, nev, passes, w, BW):
    """passes: lista (bi, wi, maska). Zwraca nom (3,5,NB), boot (nrep,3,5,NB)."""
    nrep = BW.shape[0]
    nom = np.zeros((3, 5, NB))
    boot = np.zeros((nrep, 3, 5, NB))
    for bi, wi, ps in passes:
        m_, e_, w_ = mb[ps], ev[ps], w[ps]
        nom[bi, wi] = np.bincount(m_, weights=w_, minlength=NB)
        for r in range(nrep):
            boot[r, bi, wi] = np.bincount(m_, weights=BW[r, e_] * w_, minlength=NB)
    return nom, boot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nrep", type=int, default=200)
    ap.add_argument("--nproc", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()
    P = json.load(open(os.path.join(SCAN, "K1c", "out", "prep.json")))
    NOM = json.load(open(os.path.join(SCAN, "K1c", "out", "k1c_nominal.json")))
    G["P"] = P
    G["fitvar"] = {bn: NOM["fitvar"][f"427929|{bn}"] for bn in BN}
    thr = lambda bn, wp: P["thr"][f"H|{bn}|{wp}"]  # noqa: E731
    out = dict(nrep=a.nrep, fitvar=G["fitvar"], info={}, cells={}, eps={}, mc={})
    G["hist"] = {}
    rng = np.random.default_rng(20261118)

    # ---- eps MC (fd) i błąd stat. ilorazów zagnieżdżonych
    b = np.load(os.path.join(SCAN, "cand", "train", "oof", "base.npz"))
    sS = np.load(os.path.join(SCAN, "cand", "train", "oof", "H_s12345.npz"))["S"]
    Spt, Sw = b["S_pt"].astype(float), b["S_w"].astype(float)
    for (lo, hi), bn in zip(BINS, BN):
        inb = (Spt >= lo) & (Spt < hi)
        W = Sw[inb].sum()
        p90 = inb & (sS >= thr(bn, 0.9))
        W90 = Sw[p90].sum()
        for wp in WPS:
            ps = inb & (sS >= thr(bn, wp))
            r = Sw[ps].sum() / W90
            var = ((Sw[ps] ** 2).sum() * (1 - r) ** 2 + (Sw[p90 & ~ps] ** 2).sum() * r ** 2) / W90 ** 2
            out["eps"][f"{bn}|{wp}"] = dict(eps=float(Sw[ps].sum() / W), eps_prep=P["shapes"][f"H|{bn}|{wp}"]["eps_mc"], n=int(ps.sum()),
                                           ratio90=float(r), ratio90_relerr=float(np.sqrt(var) / r))
    print("eps MC gotowe", flush=True)

    # ---- dane
    for run in RUNS:
        d = PR.load_data(run)
        s, m, pt = d["s"]["H"], d["m"], d["pt"]
        out["info"][run] = dict(n=int(len(m)), pt_min=float(pt.min()), pt_max=float(pt.max()), sH_min=float(s.min()),
                                thr90=[thr(bn, 0.9) for bn in BN])
        ib = np.digitize(pt, [3, 5, 8, 30]) - 1
        ue, einv = np.unique(d["ev"], return_inverse=True)
        mbin = np.clip(((m - FF.LO) / 0.005).astype(int), 0, NB - 1)
        BW = rng.poisson(1.0, (a.nrep, len(ue))).astype(np.float32)
        passes = [(bi, wi, (ib == bi) & (s >= thr(bn, wp))) for bi, bn in enumerate(BN) for wi, wp in enumerate(WPS)]
        nom, boot = hists(mbin, einv, len(ue), passes, np.ones(len(m)), BW)
        same = all(np.array_equal(nom[bi, wi], np.asarray(P["hists"][f"{run}|H|{bn}|{wp}"], float)) for bi, bn in enumerate(BN) for wi, wp in enumerate(WPS))
        out["info"][run]["hist_equal_prep"] = bool(same)
        out["info"][run]["n_events"] = int(len(ue))
        G["hist"][run] = dict(nom=nom, boot=boot)
        print(f"dane {run}: n={len(m)} pt∈[{pt.min():.3f},{pt.max():.2f}] min s_H={s.min():.4f} hist==prep {same}  ({time.time()-t0:.0f} s)", flush=True)
        del d, BW

    # ---- MC mb: nominalne i wstrzyknięcie
    m, pt, tr, key = b["mb_m"].astype(float), b["mb_pt"].astype(float), b["mb_truth"], b["mb_key"]
    s = np.load(os.path.join(SCAN, "cand", "train", "oof", "H_s12345.npz"))["mb"]
    ib = np.digitize(pt, [3, 5, 8, 30]) - 1
    keep = (m >= FF.LO) & (m < FF.HI) & (ib >= 0) & (ib < 3)
    t90 = np.array([thr(bn, 0.9) for bn in BN])
    keep &= s >= t90[np.clip(ib, 0, 2)]
    m, pt, tr, key, s, ib = m[keep], pt[keep], tr[keep], key[keep], s[keep], ib[keep]
    sig = tr <= 1
    ue, einv = np.unique(key, return_inverse=True)
    mbin = np.clip(((m - FF.LO) / 0.005).astype(int), 0, NB - 1)
    BW = rng.poisson(1.0, (a.nrep, len(ue))).astype(np.float32)
    passes = [(bi, wi, (ib == bi) & (s >= thr(bn, wp))) for bi, bn in enumerate(BN) for wi, wp in enumerate(WPS)]
    w_inj = np.ones(len(m))
    for bi, bn in enumerate(BN):
        n70 = float((sig & (ib == bi) & (s >= thr(bn, 0.7))).sum())
        n90 = float((sig & (ib == bi)).sum())
        q = 0.9 * n70 / n90
        w = q * (n90 - n70) / (n70 * (1 - q))
        w_inj[sig & (ib == bi) & (s >= thr(bn, 0.7))] = w
        out["mc"][f"w_inj|{bn}"] = w
    for var, w in (("mc_nom", np.ones(len(m))), ("mc_inj", w_inj)):
        nom, boot = hists(mbin, einv, len(ue), passes, w, BW)
        G["hist"][var] = dict(nom=nom, boot=boot)
        # prawda
        nt_unw = np.zeros((3, 5))
        nt = np.zeros((3, 5))
        ntb = np.zeros((a.nrep, 3, 5))
        ntb_unw = np.zeros((a.nrep, 3, 5))
        for bi, wi, ps in passes:
            q = ps & sig
            nt_unw[bi, wi] = q.sum()
            nt[bi, wi] = w[q].sum()
            for r in range(a.nrep):
                ntb[r, bi, wi] = (BW[r, einv[q]] * w[q]).sum()
                ntb_unw[r, bi, wi] = BW[r, einv[q]].sum()
        out["mc"][var] = dict(nt=nt.tolist(), nt_unw=nt_unw.tolist(), ntb=ntb.tolist(), ntb_unw=ntb_unw.tolist())
    out["mc"]["n_events"] = int(len(ue))
    print(f"MC gotowe ({time.time()-t0:.0f} s)", flush=True)

    jobs = [(src, bi, wi) for src in G["hist"] for bi in range(3) for wi in range(5)]
    with Pool(a.nproc) as p:
        for i, r in enumerate(p.imap_unordered(job, jobs, chunksize=1)):
            out["cells"][r["lab"]] = r
            if i % 10 == 9:
                print(f"  {i+1}/{len(jobs)} komórek, {time.time()-t0:.0f} s", flush=True)
    json.dump(out, open(os.path.join(TOP, "out", "scan.json"), "w"))
    print(f"koniec: {time.time()-t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""REG428 r2 (20.11.2026, projekt B; lokalnie) — fity per region triggera runu 428770 (+ MB / MB∧J12 runu 427929). POST HOC.
  yields : uzyski WP70 drzewa H per klasa (D0: M0 i M8 z κ × 1,364, test p(χ²) M8 < 0,01 → M0; D+: M0); σ = max(HESSE, rms 30 replik) gdzie są repliki
  dstar  : SF drzewa z tagu D* per region (fit_passfail z tłem progowym, btrig_common.fit_eps; progi i ε_MC z dnet4_results.json)
  sbeff  : ułamki z odjęciem pasm (R1: bit J12 / j20; 427929 MB: bit J12) w binach analizy i drobnych binach pT
Maszyneria fitu: FITBIAS_2026-11-19/skrypty/{fb_fit,fb_var}.py (bez zmian). Wyjście: ../out/{yields,dstar,sbeff}.json"""
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(CW, "FITBIAS_2026-11-19", "skrypty"))
sys.path.insert(0, os.path.join(CW, "Q_lxplus", "dnet1"))
sys.path.insert(0, os.path.join(CW, "Q_lxplus", "dnet1", "src"))          # charmtag.dstar_refit w wersji użytej w dnet4 (przed L_charmtag z venv)
import fb_fit as F  # noqa: E402
import fb_var as V  # noqa: E402

RH = os.path.join(HERE, "..", "runs", "lxplus", "reg_hist.npz")
DS = os.path.join(HERE, "..", "runs", "lxplus", "dstar.npz")
OUT = os.path.join(HERE, "..", "out")
KSC = {"D0": 1.364, "Dplus": 1.068}
J70 = 5
CLS = {"428770": ["allGRL", "R1", "R2", "R3x", "rest", "J12", "j20AFP", "J12hi", "R3xhi", "R2hi", "notJ12"], "427929": ["MB", "MB_J12", "MB_notJ12"]}
P_REJ = 0.01
G = None


def _init():
    global G
    V._init()
    G = np.load(RH, allow_pickle=False)


def hist(sp, bn, run, cls):
    """→ D [1 + nrep, nfull5] (WP70 skumulowane; repliki gdy są)."""
    if cls == "notJ12":
        a, b = hist(sp, bn, run, "allGRL"), hist(sp, bn, run, "J12")
        return a - b
    if cls == "MB_notJ12":
        return hist(sp, bn, run, "MB") - hist(sp, bn, run, "MB_J12")
    k = f"{sp}|{bn}|{run}|{cls}"
    d = [G[k][J70 + 1:].sum(0)[None]]
    if k + "|rep70" in G.files:
        d.append(G[k + "|rep70"])
    return V.full5(np.concatenate(d).astype(np.float64))


def rms(x):
    x = [v for v in x if v is not None]
    return float(np.std(x, ddof=1)) if len(x) >= 5 else float("nan")


def t_yield(t):
    sp, bn, run, cls = t
    k = V.wcat(sp, bn, 0.7)
    k["kap"] *= KSC[sp]
    D = hist(sp, bn, run, cls)
    o = {"sp": sp, "bin": bn, "run": run, "cls": cls, "n_rows": float(D[0].sum()), "nrep": int(D.shape[0] - 1)}
    for m in (["M0", "M8"] if sp == "D0" else ["M0"]):
        try:
            f = V.fit_set(D[0][None], [k], sp, m)
        except Exception as ex:  # noqa: BLE001
            o[m] = {"valid": False, "reason": str(ex)}; continue
        r = V.sres(f)
        r["p"] = F.pval(r["chi2"], r["ndf"])
        bN = []
        for i in range(1, D.shape[0]):
            try:
                fb = V.fit_set(D[i][None], [k], sp, m, start=f["x"], hesse=False)
                bN.append(fb["N"][0] if fb["valid"] else None)
            except Exception:  # noqa: BLE001
                bN.append(None)
        r["boot_rms"] = rms(bN)
        r["stat"] = max(r["err"] or 0.0, r["boot_rms"] if math.isfinite(r["boot_rms"]) else 0.0)
        o[m] = r
    if sp == "D0":
        m8 = o.get("M8", {})
        o["chosen"] = "M8" if (m8.get("valid") and m8.get("p", 0) >= P_REJ) else "M0"
    else:
        o["chosen"] = "M0"
    return o


def run_yields(nproc=6):
    tasks = [(sp, bn, run, c) for sp in V.TAG for bn in F.BINS for run in CLS for c in CLS[run]]
    t0 = time.time()
    with Pool(nproc, initializer=_init) as p:
        out = p.map(t_yield, tasks, chunksize=1)
    json.dump(out, open(os.path.join(OUT, "yields.json"), "w"), indent=1)
    print(f"yields: {len(out)} komórek, {time.time() - t0:.0f}s")


# ------------------------------------------------------------------------------------------------ D*
def bit(b, i):
    return (b & np.uint16(1 << i)) > 0


def _fit_ds(job):
    import btrig_common as BC
    key, dp, df_, shape = job
    f = BC.fit_eps(dp, df_, shape, interval=True)
    f["eps_sym"] = float(BC.eps_sym_err(f))
    return key, f


def run_dstar(nproc=6):
    import btrig_common as BC
    BC.add_paths()
    from charmtag import dstar_refit as DR
    z = np.load(DS, allow_pickle=False)
    R = json.load(open(os.path.join(CW, "DNET", "runs", "lxplus4", "dnet4_results.json")))["species"]["D0"]
    t = R["ref"]
    ed = [3, 5, 8, 30]
    fdb = np.digitize(z["fd_nd|pt"], ed) - 1
    jobs, meta = [], {}
    for bi, bn in enumerate(F.BINS):
        shape = DR.shape_from_mc(z["fd_nd|dm"][fdb == bi].astype(np.float64) * 1e3)
        for wp in (0.7, 0.5):
            th = R["cfg"][t][f"{bn}|{wp}"]["thr"]
            em, eme, _ = R["sf_eps_mc"][f"{t}|{bn}|{wp}"]
            sfd = z["fd_nd|sc"][fdb == bi]
            meta[f"{bn}|{wp}"] = {"thr": th, "eps_mc": em, "eps_mc_err": eme, "eps_mc_here": float((sfd >= th).mean()), "n_mc": int(len(sfd))}
            for run in ("428770", "427929"):
                s = f"data_{run}"
                pt, dm, sc = z[s + "|pt"], z[s + "|dm"].astype(np.float64), z[s + "|sc"]
                inb = (np.digitize(pt, ed) - 1) == bi
                if run == "428770":
                    b, lb = z[s + "|bits"], z[s + "|lb"]
                    grl, r1, j20, j12 = bit(b, 0), bit(b, 1), bit(b, 2), bit(b, 3)
                    cl = {"all": np.ones(len(b), bool), "R1": grl & r1, "R2": grl & j20 & ~r1, "R3x": grl & j12 & ~r1 & ~j20, "rest": grl & ~r1 & ~j20 & ~j12,
                          "J12": grl & j12, "notJ12": grl & ~j12, "j20AFP": grl & j20, "notR3x": grl & ~(j12 & ~r1 & ~j20)}
                else:
                    c = z[s + "|cls"]
                    mb = ((c & 1) > 0) & ((c & 2) > 0)
                    cl = {"MB": mb, "MB_J12": mb & ((c & 4) > 0), "MB_notJ12": mb & ~((c & 4) > 0)}
                for cn, cm in cl.items():
                    md = cm & inb
                    jobs.append(((run, cn, bn, wp, int(md.sum())), dm[md & (sc >= th)], dm[md & (sc < th)], shape))
    t0 = time.time()
    with Pool(nproc) as p:
        fits = p.map(_fit_ds, jobs, chunksize=1)
    out = []
    for (run, cn, bn, wp, n), f in fits:
        m = meta[f"{bn}|{wp}"]
        sf = f["eps"] / m["eps_mc"] if np.isfinite(f["eps"]) else float("nan")
        sfe = float(np.hypot(f["eps_sym"] / m["eps_mc"], sf * m["eps_mc_err"] / m["eps_mc"])) if np.isfinite(sf) else float("nan")
        out.append({"run": run, "cls": cn, "bin": bn, "wp": wp, "n_rows": n, "SF": sf, "SF_err": sfe, **{k: (v if not isinstance(v, float) or np.isfinite(v) else None) for k, v in f.items()}})
    json.dump({"meta": meta, "fits": out}, open(os.path.join(OUT, "dstar.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(f"dstar: {len(out)} fitów, {time.time() - t0:.0f}s")


# ------------------------------------------------------------------------------------------------ odjęcie pasm
def sb_counts(h1, sp):
    """h1 [..., nm] (1 MeV) → (okno, pasma) z pasm symetrycznych o łącznej szerokości = okno (k = 1)."""
    c = F.SP[sp]
    lo = c["mr"][0]
    m0, w = c["m0"], c["win"]
    gap = 0.010
    i = lambda x: int(round((x - lo) / 0.001))  # noqa: E731
    win = h1[..., i(m0 - w):i(m0 + w)].sum(-1)
    sb = h1[..., i(m0 - w - gap - w / 2):i(m0 - w - gap)].sum(-1) + h1[..., i(m0 + w + gap):i(m0 + w + gap + w / 2)].sum(-1)
    return win.astype(float), 2.0 * sb.astype(float), 4.0 * sb.astype(float)     # okno, tło w oknie (k = 2), wariancja tła (k² × n)


def frac(hp, ha, sp):
    """ε = S_pass/S_all z odjęciem pasm; pass ⊂ all; błąd z niezależnych pass i fail."""
    wp, bp, vp = sb_counts(hp, sp)
    wa, ba, va = sb_counts(ha, sp)
    sp_, sa = wp - bp, wa - ba
    sf = sa - sp_
    vsp, vsf = wp + vp, (wa - wp) + (va - vp)
    with np.errstate(divide="ignore", invalid="ignore"):
        e = sp_ / sa
        err = np.sqrt((sf ** 2 * vsp + sp_ ** 2 * vsf)) / sa ** 2
    return {"S_pass": sp_.tolist(), "S_pass_err": np.sqrt(vsp).tolist(), "S_all": sa.tolist(), "S_all_err": np.sqrt(wa + va).tolist(), "eps": e.tolist(), "eps_err": err.tolist(),
            "n_pass_win": wp.tolist(), "n_all_win": wa.tolist()}


def run_sbeff():
    z = np.load(RH, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    out = {"fine_edges": meta["fine"], "bins": F.BINS, "res": {}}
    pairs = [("427929", "MB_J12", "MB"), ("428770", "R1_J12", "R1"), ("428770", "R1_j20", "R1"), ("428770", "R1hi_j20", "R1hi"), ("428770", "R1_J12tbp", "R1"),
             ("428770", "J12_AFPtbp", "J12"), ("428770", "J12_j20AFP", "J12_AFPtbp"), ("428770", "j20AFP_J12", "j20AFP"), ("428770", "J12", "allGRL")]
    for sp in V.TAG:
        for run, cp, ca in pairs:
            hp = np.stack([z[f"{sp}|{bn}|{run}|{cp}"][J70 + 1:].sum(0) for bn in F.BINS])
            ha = np.stack([z[f"{sp}|{bn}|{run}|{ca}"][J70 + 1:].sum(0) for bn in F.BINS])
            out["res"][f"{sp}|{run}|{cp}/{ca}|bins"] = frac(hp, ha, sp)
            out["res"][f"{sp}|{run}|{cp}/{ca}|fine"] = frac(z[f"{sp}|fine|{run}|{cp}"][:, 1], z[f"{sp}|fine|{run}|{ca}"][:, 1], sp)
            out["res"][f"{sp}|{run}|{cp}/{ca}|cand_all_mass"] = {"pass": hp.sum(-1).tolist(), "all": ha.sum(-1).tolist()}
    json.dump(out, open(os.path.join(OUT, "sbeff.json"), "w"), indent=1)
    print("sbeff OK")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for c in (sys.argv[1:] or ["yields", "dstar", "sbeff"]):
        {"yields": run_yields, "dstar": run_dstar, "sbeff": run_sbeff}[c]()

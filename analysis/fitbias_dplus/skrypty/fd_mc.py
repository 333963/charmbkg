#!/usr/bin/env python3
"""FITBIAS_Dplus — próby na MC (KRYTERIA §3.1–3.3, §7): domknięcie pełne, połówki (poza próbką), rozkład na kategorie, Asimov. Bez danych.
  fd_mc.py [--nproc 7] [--nrep 40]  → out/mc.json"""
import argparse
import json
import os
import sys
import time
import zlib
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fd_lib as D  # noqa: E402
F, V = D.F, D.V
NREP = 40
DEC = ["DMM", "CM_Ds", "CM_oth", "PC", "NoCharm", "CombB"]
MODELLED = {"pk": ("DMM", "CM_Ds", "CM_oth"), "pc": ("PC",)}


def reps(name, H, sigT, T, f0, seed, n):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        f = D.fitA(name, rng.poisson(H).astype(float), sigT, T, start=f0["x"], hesse=False)
        if f["valid"]:
            out.append(f["N"][0])
    return out


def closure(name, H, nt, sigT, T, seed, nrep):
    f = D.fitA(name, H, sigT, T)
    r = D.res(f)
    bN = reps(name, H, sigT, T, f, seed, nrep)
    rms = float(np.std(bN)) if len(bN) > 0.7 * nrep else None
    e = max(1.3 * (r["err"] or 0.0), rms or 0.0)
    r.update({"N_true": nt, "b": r["N"] / nt - 1, "b_err": e / nt, "rms_rep": rms, "n_rep": len(bN)})
    ff = D.fitA(name, H, sigT, T, shape="fix")
    r["b_fix"] = ff["N"][0] / nt - 1
    return r


def task(t):
    bn, name, nrep = t
    V._init()
    sigT, _ = D.sig_tmpl(bn)
    S = D.SET[name]
    cm = {w: D.comps(w, bn) for w in ("full", "A", "B")}
    T = {w: D.tmpls(cm[w]) for w in cm}
    sd = zlib.crc32(f"{bn}|{name}".encode())
    out = {"bin": bn, "model": name, "kap": {w: {k: T[w][k][1] for k in T[w]} for w in T}}
    # 1) domknięcie pełne (szablony w próbce)
    out["full"] = closure(name, cm["full"]["all"], float(cm["full"]["Sig"].sum()), sigT, T["full"], sd, nrep)
    # 2) połówki: szablon z jednej, pseudo-dane z drugiej
    for tw, dw in (("A", "B"), ("B", "A")):
        out[f"half_{tw}to{dw}"] = closure(name, cm[dw]["all"], float(cm[dw]["Sig"].sum()), sigT, T[tw], sd + 1, nrep)
    # 3a) Asimov: N_prawda × fd + rzeczywiste DMM + CM + PC + wygładzona kombinatoryka; w próbce i poza próbką
    for tag, tw, dw in (("in", "full", "full"), ("AtoB", "A", "B"), ("BtoA", "B", "A")):
        c = cm[dw]
        nt = float(c["Sig"].sum())
        H = nt * sigT.bins(D.FE, D.C["m0"]) + c["pk"] + c["PC"] + D.smooth(c["Comb"], sigT)
        f = D.fitA(name, H, sigT, T[tw])
        r = D.res(f)
        ch = c["pk"] + c["PC"]
        base = nt * sigT.bins(D.FE, D.C["m0"]) + D.smooth(c["Comb"], sigT)
        rng = np.random.default_rng(sd + 7)
        bN = []
        for _ in range(nrep):
            fr = D.fitA(name, base + rng.poisson(ch), sigT, T[tw], start=f["x"], hesse=False)
            if fr["valid"]:
                bN.append(fr["N"][0])
        r.update({"N_true": nt, "b": r["N"] / nt - 1, "b_err": (float(np.std(bN)) / nt if len(bN) > 0.7 * nrep else None)})
        out[f"asimov_{tag}"] = r
    # 3b) rozkład na kategorie (kształt stały; składowe modelowane szablonem: w próbce = 0 z konstrukcji → liczone poza próbką A→B, B→A)
    dec = {}
    for tag, tw, dw in (("in", "full", "full"), ("AtoB", "A", "B"), ("BtoA", "B", "A")):
        c = cm[dw]
        nt = float(c["Sig"].sum())
        sA = nt * sigT.bins(D.FE, D.C["m0"])
        tA = {k: nt * T[tw][k][1] * T[tw][k][0].bins(D.FE, D.C["m0"]) for k in S["tied"]}     # Asimov składowych związanych (z szablonu)
        o = {}
        names = DEC + ([k for k in c if k.startswith("PC|")] if tag == "in" else [])
        for cn in names:
            if cn.startswith("PC|") and c[cn].sum() < 150:
                continue
            H = sA + c[cn]
            for k in S["tied"]:
                H = H + tA[k]
            # składowa modelowana: zastępuje swoją część Asimova szablonu (rzeczywiste zdarzenia zamiast szablonu)
            mk = [k for k in S["tied"] if cn in MODELLED[k] or (cn.startswith("PC|") and k == "pc")]
            if mk:
                if tag == "in":
                    o[cn] = {"modelled": True}
                    continue
                k = mk[0]
                H = H - nt * T[tw][k][1] * _sub_shape(cm[tw], cn, k)        # zdejmij przewidywany przez szablon udział tej składowej; zostają rzeczywiste zdarzenia c[cn]
            f = D.fitA(name, H, sigT, T[tw], shape="fix", hesse=False)
            rng = np.random.default_rng(zlib.crc32(f"{bn}|{name}|{tag}|{cn}".encode()))
            bN = []
            for _ in range(nrep if not cn.startswith("PC|") else nrep // 2):
                fr = D.fitA(name, H - c[cn] + rng.poisson(c[cn]), sigT, T[tw], shape="fix", start=f["x"], hesse=False)
                if fr["valid"]:
                    bN.append(fr["N"][0])
            o[cn] = {"spur": (f["N"][0] - nt) / nt, "err": float(np.std(bN)) / nt if len(bN) > 5 else None, "n": float(c[cn].sum()), "modelled": bool(mk)}
        dec[tag] = o
    out["dec"] = dec
    return out


def _sub_shape(cm_t, cn, k):
    """Udział składowej cn w szablonie k (z próbki szablonu), jako kształt w binach znormalizowany tak, że Σ po składowych = szablon/N_szablonu."""
    tot = cm_t["pk"] if k == "pk" else cm_t["PC"]
    return cm_t[cn] / max(tot.sum(), 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nproc", type=int, default=7); ap.add_argument("--nrep", type=int, default=NREP)
    a = ap.parse_args()
    t0 = time.time()
    tasks = [(bn, m, a.nrep) for bn in D.BINS for m in D.SET]
    with Pool(a.nproc) as p:
        R = []
        for i, r in enumerate(p.imap_unordered(task, tasks)):
            R.append(r)
            print(f"[{time.time() - t0:5.0f}s] {i + 1}/{len(tasks)} {r['bin']} {r['model']}", flush=True)
    R.sort(key=lambda r: (D.BINS.index(r["bin"]), r["model"]))
    json.dump(R, open(os.path.join(D.BASE, "out", "mc.json"), "w"))
    print(f"zapisano out/mc.json, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()

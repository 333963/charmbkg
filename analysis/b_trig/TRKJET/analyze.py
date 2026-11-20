#!/usr/bin/env python3
"""TRKJET — analiza lokalna (projekt B, 17.11.2026). Hipoteza H: deficyt izolacji D0 (ptcone40 < 2 GeV) w klasach triggera dżetowego bierze się
stąd, że D0 leży w dżecie / tuż obok, czego MC (fd_nd) nie ma w tej proporcji.

Wejście: runs/lxplus/tj_<próbka>.npz (etap ORCH trkjet; pary z b_trig + dżety torowe anty-kT R=0,4, stożki, półkule — Q_lxplus/trkjet/tj_kernel.py).
Narzędzia fitu/SBS/przeważania: B_trig/ISO_origin (iso_common, s2_universal_rw) — te same co w ISO_origin/RAPORT.md.

(a) ułamek sygnału z D0 w dżecie / przy dżecie, per klasa triggera i run vs MC;
(b) SF izolacji osobno „w dżecie” i „poza dżetem” (fit Δm pass/fail przy N_sig ≥ 150; inaczej SBS z oznaczeniem);
(c) SF po przeważeniu MC w zmiennych dżetowych — czy domyka się JEDNOCZEŚNIE w 8:J12only i 9:MB;
(d) ε_iso(dane) w funkcji ΔR do dżetu — uniwersalność między runami/klasami (χ²).
Werdykt (§0 raportu) liczony z kryteriów K1–K4 automatycznie; wynik negatywny ma tę samą formę co pozytywny.

UWAGA o błędnym kole: ptcone40 to Σ pT torów w ΔR < 0,4, więc dżet zbudowany z tych samych torów (kolekcje A, B) jest z izolacją związany
z konstrukcji. Test H bez koła = kolekcja O (dżety tylko z torów POZA stożkiem 0,4) + półkula przeciwna + pierścień 0,4–1,0.
Schematy/definicje oznaczone „[koło]” służą jako sufit (co najwyżej tyle da się domknąć), nie jako dowód.

Użycie: .venv_h5/bin/python charm_work/B_trig/TRKJET/analyze.py [--in runs/lxplus] [--out out] [--report RAPORT.md] [--nproc 6]
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "ISO_origin"))
from iso_common import C, DR, chi2_const, classes, eff_w, fit_pf, sbs_w  # noqa: E402
import s2_universal_rw as S2  # noqa: E402

INF = 1e9
MWIN, STHR = 0.040, 0.6
CUT = ("ptcone40", 2.0)
NFIT = 150          # N_sig (SBS) od którego robimy fit pass/fail
NWEAK, NBAD = 150, 50
PTB = {"3-30": (3, 30), "3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}
CLS = ["9:MB", "9:MB&!j20&!J12", "9:MB&j20&!J12", "8:all", "8:J12only", "8:j20&!J12", "8:MBTSbased", "8:other"]
DISJ = ["9:MB", "8:J12only", "8:MBTSbased", "8:other"]
KEY2 = ["9:MB", "8:J12only"]     # klasy, w których H ma domknąć się jednocześnie

# (a) własności „D0 w dżecie / przy dżecie”: nazwa → (opis, funkcja maski, koło?)
PROPS = {
    "A5":  ("dżet zawierający D0 (kol. A) pT > 5 GeV i ΔR(D0, oś) < 0,4", lambda D: (D["a_jD_pt"] > 5) & (D["a_dr_jD"] < 0.4), True),
    "A10": ("jw., pT > 10 GeV", lambda D: (D["a_jD_pt"] > 10) & (D["a_dr_jD"] < 0.4), True),
    "A20": ("jw., pT > 20 GeV", lambda D: (D["a_jD_pt"] > 20) & (D["a_dr_jD"] < 0.4), True),
    "z<0.5": ("z = pT(D0)/pT(dżetu z D0) < 0,5", lambda D: D["z"] < 0.5, True),
    "B2":  ("ΔR(D0, dżet otoczenia kol. B, pT > 2 GeV) < 0,4", lambda D: D["b_drn2"] < 0.4, True),
    "B3":  ("jw., pT > 3 GeV", lambda D: D["b_drn3"] < 0.4, True),
    "B5":  ("jw., pT > 5 GeV", lambda D: D["b_drn5"] < 0.4, True),
    "B8":  ("jw., pT > 8 GeV", lambda D: D["b_drn8"] < 0.4, True),
    "O3":  ("dżet spoza stożka (kol. O) pT > 3 GeV w ΔR < 0,8", lambda D: D["o_ptmax08"] > 3, False),
    "O5":  ("jw., pT > 5 GeV", lambda D: D["o_ptmax08"] > 5, False),
    "O10": ("jw., pT > 10 GeV", lambda D: D["o_ptmax08"] > 10, False),
    "away10": ("Σ pT torów w półkuli przeciwnej > 10 GeV", lambda D: D["h_s_away"] > 10, False),
}
# (b) podziały w/poza dżetem
SPLITS = {"A10 [koło]": "A10", "B3 [koło]": "B3", "O5": "O5", "O10": "O10"}
# (c) schematy przeważania MC: (nazwa, [(zmienna, krawędzie)…], koło?)
PT4 = ("pt", [3, 5, 8, 12, 30])
SCHEMES = [
    ("brak", [], False),
    ("pT", [("pt", [3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 15, 20, 30])], False),
    ("pT×lead_pt×nSelBL (ISO_origin)", [PT4, ("e_lead_pt", [0, 5, 8, 12, INF]), ("e_nSelBL", [0, 35, 55, INF])], False),
    ("pT×O: pT dżetu w ΔR<0,8", [PT4, ("o_ptmax08", [0, 1e-6, 2, 3, 5, 8, 12, 20, INF])], False),
    ("pT×O: (ΔR, pT) najbliższego dżetu pT>3", [PT4, ("o_drn3", [0, 0.6, 0.8, 1.1, 1.6, INF]), ("o_ptn3", [0, 5, 8, 12, 20, INF])], False),
    ("pT×O: (ΔR, pT) najbliższego dżetu pT>5", [PT4, ("o_drn5", [0, 0.6, 0.8, 1.1, 1.6, INF]), ("o_ptn5", [0, 8, 12, 20, INF])], False),
    ("pT×Σ pT półkuli przeciwnej", [PT4, ("h_s_away", [0, 4, 7, 10, 15, 22, 35, INF])], False),
    ("pT×pierścień 0,4–1,0 (własny)", [PT4, ("ring", [0, 1, 2, 4, 6, 9, 13, 20, INF])], False),
    ("pT×O pT(ΔR<0,8)×nSelBL", [PT4, ("o_ptmax08", [0, 1e-6, 3, 6, 12, INF]), ("e_nSelBL", [0, 35, 55, INF])], False),
    ("pT×O pT(ΔR<0,8)×pierścień", [PT4, ("o_ptmax08", [0, 1e-6, 3, 6, 12, INF]), ("ring", [0, 2, 5, 10, INF])], False),
    ("pT×B: (ΔR, pT) najbliższego dżetu pT>3 [koło]", [PT4, ("b_drn3", [0, 0.2, 0.4, 0.7, 1.1, INF]), ("b_ptn3", [0, 5, 8, 12, 20, INF])], True),
    ("pT×pT dżetu z D0 (kol. A) [koło]", [PT4, ("a_jD_pt", [0, 4, 5, 6.5, 8, 10, 13, 17, 25, 40, INF])], True),
]
# (d) ε_iso(ΔR): zmienna → (opis, krawędzie, koło?)
DRV = {
    "o_drn5": ("ΔR(D0, najbliższy dżet kol. O pT > 5)", [0, 0.6, 0.8, 1.0, 1.4, 2.0, 3.2, 50, INF], False),
    "o_drn3": ("ΔR(D0, najbliższy dżet kol. O pT > 3)", [0, 0.6, 0.8, 1.0, 1.4, 2.0, 3.2, 50, INF], False),
    "b_drn5": ("ΔR(D0, najbliższy dżet kol. B pT > 5) [koło dla ΔR < 0,4]", [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 2.0, 3.2, 50, INF], True),
    "a_dr_j1": ("ΔR(D0, dżet wiodący kol. A)", [0, 0.1, 0.4, 0.8, 1.4, 2.4, 2.9, INF], True),
    "o_ptmax08": ("pT dżetu kol. O w ΔR < 0,8 [GeV]", [0, 1e-6, 3, 5, 8, 12, 20, INF], False),
    "h_s_away": ("Σ pT półkuli przeciwnej [GeV]", [0, 4, 7, 10, 15, 22, 35, INF], False),
}
DR_MAIN = ["9:MB", "8:J12only", "8:MBTSbased", "8:other"]


# ------------------------------------------------------------------------------------------------------------------ wczytanie
def load(indir, sample):
    p = os.path.join(indir, f"tj_{sample}.npz")
    z = np.load(p, allow_pickle=False)
    D = {k: z[k] for k in z.files if k != "meta"}
    meta = json.loads(str(z["meta"]))
    for k in list(D):
        if k.startswith("e_tp_"):
            D[k[2:]] = D[k]
    D["rel40"] = D["ptcone40"] / D["pt"]
    D["ann"] = D["ptcone100"] - D["ptcone40"]
    D["ring"] = D["c_s10"] - D["c_s04"]
    D["dok"] = (D["nmatch"] + D["nd_lo"] == 2) & (D["same_entry"] == 1)
    return D, meta


def base(D, lo=3, hi=30):
    return (np.abs(D["m"] - C.M0_D0) < MWIN) & (D["score"] > STHR) & (D["pt"] >= lo) & (D["pt"] < hi) & D["dok"]


def mc_sel(M, lo=3, hi=30):
    return base(M, lo, hi) & np.isin(M["truth"], [0, 1]) & (M["dm"] > 0.1435) & (M["dm"] < 0.1475)


def make_classes(D9, D8):
    c9 = classes(D9, 427929); c8 = classes(D8, 428770)
    t9 = lambda c: D9.get("tp_" + c, np.zeros(len(D9["m"]), np.int8)) == 1  # noqa: E731
    mb = c9["MB427929"]; j9 = t9("HLT_noalg_L1J12"); j20 = t9("HLT_j20_L1MBTS_2")
    return {"9:MB": (D9, mb), "9:MB&!j20&!J12": (D9, mb & ~j20 & ~j9), "9:MB&j20&!J12": (D9, mb & j20 & ~j9),
            "8:all": (D8, c8["all"]), "8:J12only": (D8, c8["J12only"]), "8:j20&!J12": (D8, c8["j20_notJ12"]),
            "8:MBTSbased": (D8, c8["MBTSbased"]), "8:other": (D8, c8["other"])}


def sbs_eff(D, m, flag):
    """ε (albo ułamek) sygnału z własnością flag wśród kandydatów m — odejmowanie pasm Δm (iso_common.sbs_w)."""
    if m.sum() < 20:
        return {"ncand": int(m.sum()), "S": 0.0, "S_err": 0.0, "eps": np.nan, "err": np.nan}
    w = sbs_w(D["dm"][m], D["pt"][m], flag[m])
    e, er, S, eS = eff_w(w, flag[m].astype(float))
    return {"ncand": int(m.sum()), "S": S, "S_err": eS, "eps": e, "err": er}


def mc_eff(flag, m):
    n = int(m.sum())
    if n == 0:
        return {"S": 0.0, "eps": np.nan, "err": np.nan}
    e = float(flag[m].mean())
    return {"S": float(n), "eps": e, "err": float(np.sqrt(max(e * (1 - e), 0.25 / n) / n))}


def quality(r):
    S, eS = r.get("S", 0.0), r.get("S_err", 0.0)
    if not np.isfinite(r.get("eps", np.nan)) or S < NBAD or (eS > 0 and S / eS < 3):
        return "niemiarodajny"
    return "słaby" if S < NWEAK else "ok"


def _fit(job):
    key, dp, df_, shape = job
    try:
        return key, fit_pf(dp, df_, shape, interval=False)
    except Exception as ex:  # noqa: BLE001
        return key, {"eps": np.nan, "err": np.nan, "N": 0.0, "N_err": 0.0, "det": False, "valid": False, "method": "fit", "reason": str(ex)}


def best(r):
    """Wynik do tabel: fit, gdy jest, zbieżny i określony; inaczej SBS."""
    f = r.get("fit")
    if f and f.get("det") and f.get("valid") and np.isfinite(f.get("eps", np.nan)) and np.isfinite(f.get("err", np.nan)) and f["err"] > 0:
        return f["eps"], f["err"], "fit"
    return r.get("eps", np.nan), r.get("err", np.nan), "SBS"


def ratio(a, ea, b, eb):
    if not (np.isfinite(a) and np.isfinite(b)) or b <= 0:
        return np.nan, np.nan
    return a / b, float(np.hypot(ea / b, a * eb / b ** 2))


# ------------------------------------------------------------------------------------------------------------------ formatowanie
def pm(v, e, nd=3):
    if v is None or not np.isfinite(v):
        return "–"
    return f"{v:.{nd}f} ± {e:.{nd}f}".replace(".", ",") if e is not None and np.isfinite(e) else f"{v:.{nd}f}".replace(".", ",")


def num(v, nd=2):
    return "–" if v is None or not np.isfinite(v) else f"{v:.{nd}f}".replace(".", ",")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", default=os.path.join(HERE, "runs", "lxplus"))
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    ap.add_argument("--report", default=os.path.join(HERE, "RAPORT.md"))
    ap.add_argument("--nproc", type=int, default=6)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    D9, m9 = load(a.indir, "data_427929"); D8, m8 = load(a.indir, "data_428770")
    M, mM = load(a.indir, "fd_nd"); MB, mMB = load(a.indir, "mb_nd")
    metas = {"data_427929": m9, "data_428770": m8, "fd_nd": mM, "mb_nd": mMB}
    S = make_classes(D9, D8)
    mm0 = mc_sel(M)
    shape = DR.shape_from_mc(M["dm"][mm0] * 1e3)
    pasM = M[CUT[0]] < CUT[1]
    for D in (D9, D8, M, MB):
        D["_pas"] = D[CUT[0]] < CUT[1]
    partial = any(m.get("n_files", 0) < n for m, n in ((m8, 97), (m9, 458), (mM, 24), (mMB, 78)))
    R = {"inputs": {k: {"n_pairs": v["n_pairs"], "n_files": v.get("n_files"), "frac_daughters_ok": v.get("frac_daughters_ok"),
                        "frac_same_entry": v.get("frac_same_entry"), "frac_pis_match": v.get("frac_pis_match"),
                        "ptcone40_own_vs_branch": v.get("ptcone40_own_vs_branch")} for k, v in metas.items()},
         "partial": partial, "jet": m8.get("jet"), "mc_n": int(mm0.sum()), "mc_eps": float(pasM[mm0].mean()) if mm0.any() else None,
         "classes": {}, "a": {}, "b": {}, "c": {}, "c_pt": {}, "d": {}, "d_chi2": {}, "valid": {}}
    # zgodność decyzji izolacji własnej z gałęzią
    for nm, D in (("data_427929", D9), ("data_428770", D8), ("fd_nd", M), ("mb_nd", MB)):
        ok = D["dok"]
        if ok.any():
            d = D["v_s40_p2"][ok] - D["ptcone40"][ok]
            R["valid"][nm] = {"exact": float(np.mean(np.abs(d) < 1e-3)), "lt_200MeV": float(np.mean(np.abs(d) < 0.2)),
                              "own_le_branch": float(np.mean(d < 1e-3)), "same_decision_2GeV": float(np.mean((D["v_s40_p2"][ok] < 2) == (D["ptcone40"][ok] < 2))),
                              "dok": float(ok.mean())}
    jobs = []
    # ---- klasy: ε_iso, SF bez wag (odtworzenie ISO_origin)
    for name in CLS:
        D, cm = S[name]
        for pb, (lo, hi) in PTB.items():
            m = base(D, lo, hi) & cm
            r = sbs_eff(D, m, D["_pas"])
            r["mc"] = mc_eff(pasM, mc_sel(M, lo, hi))
            R["classes"][f"{name}|{pb}"] = r
            if r["S"] >= NFIT:
                jobs.append((("classes", f"{name}|{pb}"), D["dm"][m & D["_pas"]], D["dm"][m & ~D["_pas"]], shape))
    # ---- (a) ułamki
    for pn, (desc, fn, circ) in PROPS.items():
        fM = fn(M); fMB = fn(MB)
        R["a"][pn] = {"desc": desc, "circ": circ, "MC fd_nd": mc_eff(fM, mm0), "MC mb_nd": mc_eff(fMB, mc_sel(MB)), "cls": {}}
        for name in CLS:
            D, cm = S[name]
            R["a"][pn]["cls"][name] = sbs_eff(D, base(D) & cm, fn(D))
    # ---- (b) SF w/poza dżetem
    for sn, pn in SPLITS.items():
        fn = PROPS[pn][1]
        fM = fn(M)
        for side, neg in (("w", False), ("poza", True)):
            mcm = mm0 & (fM != neg)
            for name in CLS:
                D, cm = S[name]
                f = fn(D)
                m = base(D) & cm & (f != neg)
                r = sbs_eff(D, m, D["_pas"])
                r["mc"] = mc_eff(pasM, mcm)
                R["b"][f"{sn}|{side}|{name}"] = r
                if r["S"] >= NFIT:
                    jobs.append((("b", f"{sn}|{side}|{name}"), D["dm"][m & D["_pas"]], D["dm"][m & ~D["_pas"]], shape))
    # ---- (c) przeważenie
    for sname, scheme, circ in SCHEMES:
        for name in CLS:
            D, cm = S[name]
            for pb, (lo, hi) in PTB.items():
                if pb != "3-30" and name not in ("9:MB", "8:J12only", "8:MBTSbased"):
                    continue
                md = base(D, lo, hi) & cm; mm = mc_sel(M, lo, hi)
                if md.sum() < 50 or mm.sum() < 50:
                    continue
                try:
                    r = S2.sf_rw(D, md, D["_pas"], M, mm, pasM, scheme)
                except Exception as ex:  # noqa: BLE001
                    r = None; print("sf_rw", sname, name, pb, ex)
                if r:
                    (R["c"] if pb == "3-30" else R["c_pt"])[f"{sname}|{name}" + ("" if pb == "3-30" else f"|{pb}")] = r
    for sname, scheme, circ in SCHEMES:
        v = [R["c"].get(f"{sname}|{n}") for n in DISJ]
        v = [x for x in v if x and x["S"] > NBAD and x["err"] > 0]
        chi = float(sum(((x["SF"] - 1) / x["err"]) ** 2 for x in v))
        R["c"][f"{sname}|chi2"] = {"chi2": chi, "ndf": len(v), "circ": circ}
    # ---- (d) ε_iso(zmienna)
    for var, (desc, ed, circ) in DRV.items():
        R["d"][var] = {"desc": desc, "edges": ed, "circ": circ, "cls": {}, "MC": []}
        for lo, hi in zip(ed[:-1], ed[1:]):
            R["d"][var]["MC"].append(mc_eff(pasM, mm0 & (M[var] >= lo) & (M[var] < hi)))
        for name in DR_MAIN:
            D, cm = S[name]
            b = base(D) & cm
            R["d"][var]["cls"][name] = [sbs_eff(D, b & (D[var] >= lo) & (D[var] < hi), D["_pas"]) for lo, hi in zip(ed[:-1], ed[1:])]
        cl = R["d"][var]["cls"]
        for x, y in (("9:MB", "8:J12only"), ("9:MB", "8:MBTSbased"), ("8:J12only", "8:MBTSbased"), ("8:J12only", "8:other")):
            R["d_chi2"][f"{var}|{x}|{y}"] = S2.chi2_pair(cl[x], cl[y])
        for x in DR_MAIN:
            R["d_chi2"][f"{var}|{x}|MC"] = S2.chi2_pair(cl[x], R["d"][var]["MC"])
    print(f"fity: {len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    if jobs:
        with Pool(a.nproc) as P:
            for (sec, key), f in P.imap_unordered(_fit, [(k, dp, df_, sh) for k, dp, df_, sh in jobs]):
                R[sec][key]["fit"] = f

    # ------------------------------------------------------------------------------------------------------------ kryteria
    def cls_eps(name, pb="3-30"):
        return best(R["classes"][f"{name}|{pb}"])

    K = {}
    # K1: D0 częściej w dżecie w J12only niż w 9:MB i MC (bez koła: O5; z kołem: A10)
    for pn in ("O5", "A10"):
        x, y, mc = R["a"][pn]["cls"]["8:J12only"], R["a"][pn]["cls"]["9:MB"], R["a"][pn]["MC fd_nd"]
        d = x["eps"] - y["eps"]; e = np.hypot(x["err"], y["err"])
        K[f"K1_{pn}"] = {"J12only": x["eps"], "J12only_err": x["err"], "MB9": y["eps"], "MB9_err": y["err"], "MC": mc["eps"], "MC_err": mc["err"],
                         "nsigma": float(d / e) if np.isfinite(d) and e > 0 else None,
                         "pass": bool(np.isfinite(d) and e > 0 and d / e > 3)}
    # K2: SF „poza dżetem” (O5) zgodny z 1 w J12only i 9:MB; SF „w dżecie” też (H mówi: przy tej samej topologii MC opisuje izolację)
    for sn in ("O5", "O10"):
        o = {}
        for side in ("w", "poza"):
            for name in KEY2:
                r = R["b"][f"{sn}|{side}|{name}"]
                e, er, meth = best(r)
                sf, sfe = ratio(e, er, r["mc"]["eps"], r["mc"]["err"])
                o[f"{side}|{name}"] = {"SF": sf, "err": sfe, "q": quality(r), "pull": float((sf - 1) / sfe) if np.isfinite(sf) and sfe > 0 else None}
        ok = [v for v in o.values() if v["pull"] is not None and v["q"] != "niemiarodajny"]
        K[f"K2_{sn}"] = {"cells": o, "n_usable": len(ok), "pass": bool(len(ok) >= 3 and all(abs(v["pull"]) < 2 for v in ok)),
                         "max_abs_pull": max((abs(v["pull"]) for v in ok), default=None)}
    # K3: schemat bez koła domyka J12only i 9:MB jednocześnie (|SF−1| < 2σ w obu) i poprawia χ² względem ISO_origin
    ref = R["c"]["pT×lead_pt×nSelBL (ISO_origin)|chi2"]
    cand = []
    for sname, scheme, circ in SCHEMES:
        x = [R["c"].get(f"{sname}|{n}") for n in KEY2]
        if not all(x):
            continue
        pulls = [(v["SF"] - 1) / v["err"] for v in x]
        ch = R["c"][f"{sname}|chi2"]
        cand.append({"scheme": sname, "circ": circ, "SF_MB9": x[0]["SF"], "err_MB9": x[0]["err"], "SF_J12": x[1]["SF"], "err_J12": x[1]["err"],
                     "pulls": pulls, "closes": bool(all(abs(p) < 2 for p in pulls)), "chi2": ch["chi2"], "ndf": ch["ndf"],
                     "fallback_J12": x[1]["frac_fallback"], "dist_J12": abs(x[1]["SF"] - 1)})
    nc = [c for c in cand if not c["circ"] and c["scheme"] not in ("brak", "pT", "pT×lead_pt×nSelBL (ISO_origin)")]
    bst = min(nc, key=lambda c: c["chi2"]) if nc else None
    refc = next((c for c in cand if c["scheme"].startswith("pT×lead_pt×nSelBL")), None)
    K["K3"] = {"best": bst, "ref": refc, "ref_chi2": ref, "ceiling": min((c for c in cand if c["circ"]), key=lambda c: c["chi2"], default=None),
               "pass": bool(bst and bst["closes"]),
               "improves": bool(bst and refc and bst["dist_J12"] < 0.5 * refc["dist_J12"])}
    # K4: ε_iso(ΔR do dżetu O) uniwersalne między 9:MB i J12only (p > 0,01), ndf ≥ 3
    k4 = {v: R["d_chi2"][f"{v}|9:MB|8:J12only"] for v in ("o_drn5", "o_drn3", "o_ptmax08")}
    us = [x for x in k4.values() if x["ndf"] >= 3]
    K["K4"] = {"tests": k4, "pass": bool(us and all(x["p"] is not None and x["p"] > 0.01 for x in us)), "n_usable": len(us)}
    usable = K["K3"]["best"] is not None and R["classes"]["8:J12only|3-30"]["S"] >= NWEAK and R["classes"]["9:MB|3-30"]["S"] >= NWEAK
    if not usable:
        verdict = "NIEROZSTRZYGNIĘTE — za mała statystyka (bieg testowy / niepełne wejście)"
    elif K["K3"]["pass"] and K["K4"]["pass"] and K["K1_O5"]["pass"]:
        verdict = "H POTWIERDZONA — topologia dżetowa (zmienne bez koła) domyka SF w J12only i 9:MB jednocześnie, ε_iso(ΔR) uniwersalne"
    elif K["K3"]["pass"] or (K["K3"]["improves"] and K["K1_O5"]["pass"]):
        verdict = "H CZĘŚCIOWO — zmienne dżetowe zmniejszają resztę co najmniej 2×, ale nie wszystkie kryteria spełnione"
    else:
        verdict = "H ODRZUCONA — przeważenie w zmiennych dżetowych (bez koła) nie domyka J12only lepiej niż pT×lead_pt×nSelBL"
    R["K"] = K; R["verdict"] = verdict
    json.dump(R, open(os.path.join(a.out, "trkjet.json"), "w"), indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))

    # ------------------------------------------------------------------------------------------------------------ raport
    L = []
    w = L.append
    w("# TRKJET — czy deficyt izolacji D⁰ w trygerze dżetowym to „D⁰ w dżecie”? (projekt B)\n")
    w(f"Wygenerowane: `analyze.py`, wejście `{os.path.relpath(a.indir, HERE)}`, {time.strftime('%Y-%m-%d %H:%M')}. Liczby: `{os.path.relpath(a.out, HERE)}/trkjet.json`.")
    if partial:
        w("\n**UWAGA: wejście NIEPEŁNE (bieg testowy na kilku plikach) — liczby służą tylko sprawdzeniu, że analiza działa; werdykt nie jest wynikiem.**")
    w("\n## 0. Werdykt\n")
    w(f"**{verdict}.**\n")
    w("| kryterium | co sprawdza | liczby | wynik |\n|---|---|---|---|")
    for pn in ("O5", "A10"):
        k = K[f"K1_{pn}"]
        dsc = {"O5": "dżet spoza stożka (kol. O) pT > 5 GeV w ΔR < 0,8", "A10": "dżet zawierający D⁰ (kol. A) pT > 10 GeV"}[pn]
        w(f"| K1 ({pn}{' [koło]' if PROPS[pn][2] else ''}) | ułamek sygnału: {dsc}; J12only > 9:MB o > 3σ | J12only {pm(k['J12only'], k['J12only_err'])}, "
          f"9:MB {pm(k['MB9'], k['MB9_err'])}, MC {pm(k['MC'], k['MC_err'])}; różnica {num(k['nsigma'], 1)}σ | {'TAK' if k['pass'] else 'NIE'} |")
    for sn in ("O5", "O10"):
        k = K[f"K2_{sn}"]
        cells = "; ".join(f"{c}: {pm(v['SF'], v['err'], 2)}{'' if v['q'] == 'ok' else ' (' + v['q'] + ')'}" for c, v in k["cells"].items())
        w(f"| K2 ({sn}) | SF_iso osobno w/poza dżetem zgodny z 1 (|pull| < 2) w 9:MB i J12only | {cells} | {'TAK' if k['pass'] else 'NIE'} (miarodajnych {k['n_usable']}/4) |")
    k = K["K3"]
    if k["best"]:
        b_, r_ = k["best"], k["ref"]
        w(f"| K3 | przeważenie MC bez koła domyka 9:MB i J12only jednocześnie (|SF−1| < 2σ) | najlepszy: „{b_['scheme']}”: 9:MB {pm(b_['SF_MB9'], b_['err_MB9'])}, "
          f"J12only {pm(b_['SF_J12'], b_['err_J12'])}, χ²/ndf {num(b_['chi2'], 1)}/{b_['ndf']}; odniesienie pT×lead_pt×nSelBL: J12only "
          f"{pm(r_['SF_J12'], r_['err_J12']) if r_ else '–'}, χ² {num(k['ref_chi2']['chi2'], 1)}/{k['ref_chi2']['ndf']} | "
          f"{'TAK' if k['pass'] else 'NIE'}{' (reszta ≥ 2× mniejsza)' if k['improves'] and not k['pass'] else ''} |")
        if k["ceiling"]:
            c_ = k["ceiling"]
            w(f"| K3 sufit [koło] | to samo ze zmienną z torów stożka | „{c_['scheme']}”: 9:MB {pm(c_['SF_MB9'], c_['err_MB9'])}, J12only {pm(c_['SF_J12'], c_['err_J12'])}, "
              f"χ² {num(c_['chi2'], 1)}/{c_['ndf']} | (informacyjnie) |")
    k = K["K4"]
    w("| K4 | ε_iso(dane) w funkcji ΔR/pT dżetu kol. O: 9:MB vs J12only, p > 0,01 | " +
      "; ".join(f"{v}: χ²/ndf {num(x['chi2'], 1)}/{x['ndf']}" + (f" (p = {x['p']:.3g})" if x['p'] is not None else "") for v, x in k["tests"].items()) +
      f" | {'TAK' if k['pass'] else 'NIE'} (testów z ndf ≥ 3: {k['n_usable']}) |")
    w("\nJak czytać: H wymaga K3 i K4 (przy tej samej topologii dżetowej MC opisuje izolację i dane różnych klas są zgodne). K1 sam nie wystarcza "
      "(tryger dżetowy z definicji wybiera zdarzenia z dżetem). Jeżeli K3 = NIE przy K1 = TAK: D⁰ rzeczywiście leży częściej przy dżecie, "
      "ale to nie tłumaczy deficytu — różnica siedzi w samym stożku (fragmentacja wokół D⁰ w danych ≠ MC przy tym samym otoczeniu).")
    w("\n## 1. Metoda i zastrzeżenia\n")
    j = R["jet"] or {}
    w(f"- Pary D⁰←D* i score H wzięte z `b_trig` (`bt_<próbka>.npz`), bez ponownego liczenia. Selekcja: |m(Kπ) − m₀| < {MWIN * 1e3:.0f} MeV, score_H > {STHR}, "
      f"pT 3–30 GeV, komplet córek wycięty. 427929: MBTS_1 ∧ GRL. Cięcie izolacji: {CUT[0]} < {CUT[1]} GeV (gałąź ntupla).")
    w(f"- Dżety: {j.get('algo', 'anty-kT')}, R = {j.get('R')}, tory pT > {j.get('pt_min')} GeV, |η| < {j.get('eta_max')}, masa pionu. **fastjet nie jest dostępny** w `q_env` "
      "(LCG_110_swan) → własna implementacja w numba (`tj_kernel.py`, test jednostkowy w pliku).")
    w("- Kolekcje: **A** — D⁰ jako jedna pseudocząstka zamiast córek (dżet z D⁰, z); **B** — tory bez córek D⁰ i wolnego pionu; "
      "**O** — jak B, ale tylko tory z ΔR(tor, D⁰) > 0,4. Zmienne A i B zawierają tory stożka izolacji → związek z ptcone40 z konstrukcji (**[koło]**); "
      "test hipotezy opiera się na O, półkuli przeciwnej i pierścieniu 0,4–1,0.")
    w("- ε i ułamki: odejmowanie pasm Δm (SBS, `ISO_origin/iso_common.sbs_w`); fit jednoczesny pass/fail (`charmtag.dstar_refit.fit_passfail`, HESSE) "
      f"przy N_sig(SBS) ≥ {NFIT}. Oznaczenia: „słaby” N_sig < {NWEAK}, „niemiarodajny” N_sig < {NBAD} albo N_sig/σ < 3.")
    w(f"- MC: fd_nd, truth ∈ {{0,1}}, Δm 143,5–147,5 MeV, ta sama selekcja; N = {R['mc_n']}, ε_MC = {num(R['mc_eps'], 3)}. mb_nd tylko w (a).")
    w("- Przeważenie: SF = Σw·pass / Σw·ε_MC(komórka) (`s2_universal_rw.sf_rw`), komórki z n_MC < 25 cofają się do mniejszej liczby wymiarów "
      "(kolumna „cof.” = udział takiego sygnału). **Zastrzeżenie:** fd_nd ma mało D⁰ przy twardych dżetach — duże „cof.” oznacza, że schemat nie jest testem.")
    w("\n**Zgodność własnego ptcone z gałęzią `ptcone40`** (Σ pT wszystkich torów drzewa `Tracks` w ΔR < 0,4 bez córek D⁰, wolny pion wliczony):\n")
    w("| próbka | par | komplet córek | identyczne (< 1 MeV) | |Δ| < 0,2 GeV | własne ≤ gałąź | ta sama decyzja (< 2 GeV) |\n|---|---|---|---|---|---|---|")
    for nm, v in R["valid"].items():
        w(f"| {nm} | {R['inputs'][nm]['n_pairs']} | {num(v['dok'], 4)} | {num(v['exact'], 3)} | {num(v['lt_200MeV'], 3)} | {num(v['own_le_branch'], 4)} | {num(v['same_decision_2GeV'], 4)} |")
    w("\nRóżnica jest jednostronna (własne ≤ gałąź, kwanty 0,1–0,2 GeV): gałąź liczy tory od pT > 0,1 GeV, drzewo `Tracks` ma pT > 0,2 GeV. "
      "Gałąź zawiera wolny pion D* (wersja bez pionu zgadza się w < 1 %).")
    w("\n## 2. Odtworzenie punktu wyjścia — ε_iso i SF bez wag per klasa\n")
    w("| klasa | N_sig (SBS) | ε SBS | ε fit | ε_MC | SF | jakość |\n|---|---|---|---|---|---|---|")
    for name in CLS:
        r = R["classes"][f"{name}|3-30"]
        e, er, meth = best(r)
        sf, sfe = ratio(e, er, r["mc"]["eps"], r["mc"]["err"])
        f = r.get("fit")
        w(f"| {name} | {pm(r['S'], r['S_err'], 0)} | {pm(r['eps'], r['err'])} | {pm(f['eps'], f['err']) if f else '–'} | {num(r['mc']['eps'], 3)} | {pm(sf, sfe, 2)} ({meth}) | {quality(r)} |")
    w("\nOczekiwane z ISO_origin (pełna statystyka): 9:MB 0,98 ± 0,04; 8:J12only 0,46 ± 0,03; 8:MBTSbased 0,62 ± 0,07; 8:other 0,80 ± 0,09.")
    w("\n## 3. (a) Ułamek sygnału D⁰ w dżecie / przy dżecie\n")
    w("| własność | " + " | ".join(CLS) + " | MC fd_nd | MC mb_nd |\n|---|" + "---|" * (len(CLS) + 2))
    for pn, v in R["a"].items():
        w(f"| **{pn}**{' [koło]' if v['circ'] else ''}: {v['desc']} | " +
          " | ".join(pm(v["cls"][n]["eps"], v["cls"][n]["err"]) + ("" if quality(v["cls"][n]) == "ok" else " *") for n in CLS) +
          f" | {pm(v['MC fd_nd']['eps'], v['MC fd_nd']['err'])} | {pm(v['MC mb_nd']['eps'], v['MC mb_nd']['err'])} |")
    w("\n\\* N_sig < 150 albo niemiarodajny. Ułamki są dla sygnału (wagi SBS), nie dla wszystkich kandydatów.")
    w("\n## 4. (b) SF izolacji osobno „w dżecie” i „poza dżetem”\n")
    w("SF = ε_dane / ε_MC w tej samej kategorii (MC bez wag). Dla definicji [koło] ε „w dżecie” jest małe z konstrukcji — patrz tylko na „poza”.\n")
    w("| podział | strona | klasa | N_sig | ε dane | metoda | ε MC (N) | SF | jakość |\n|---|---|---|---|---|---|---|---|---|")
    for sn in SPLITS:
        for side in ("w", "poza"):
            for name in CLS:
                r = R["b"][f"{sn}|{side}|{name}"]
                e, er, meth = best(r)
                sf, sfe = ratio(e, er, r["mc"]["eps"], r["mc"]["err"])
                w(f"| {sn} | {side} | {name} | {pm(r['S'], r['S_err'], 0)} | {pm(e, er)} | {meth} | {num(r['mc']['eps'], 3)} ({int(r['mc']['S'])}) | {pm(sf, sfe, 2)} | {quality(r)} |")
    w("\n## 5. (c) SF po przeważeniu MC w zmiennych dżetowych\n")
    w("| schemat | " + " | ".join(CLS) + " | χ²/ndf vs 1 (4 klasy rozłączne) | cof. J12only |\n|---|" + "---|" * (len(CLS) + 2))
    for sname, scheme, circ in SCHEMES:
        row = []
        for n in CLS:
            x = R["c"].get(f"{sname}|{n}")
            row.append(pm(x["SF"], x["err"]) + ("" if x["S"] >= NWEAK else " *") if x else "–")
        ch = R["c"][f"{sname}|chi2"]; xj = R["c"].get(f"{sname}|8:J12only")
        w(f"| {sname} | " + " | ".join(row) + f" | {num(ch['chi2'], 1)}/{ch['ndf']} | {num(xj['frac_fallback'], 2) if xj else '–'} |")
    w("\n\\* N_sig < 150. Klasy rozłączne do χ²: " + ", ".join(DISJ) + ".")
    w("\nReszta w binach pT(D⁰) (ISO_origin: J12only 0,97 / 0,72 / 0,50 po pT×lead_pt×nSelBL):\n")
    w("| schemat | klasa | 3–5 | 5–8 | 8–30 |\n|---|---|---|---|---|")
    for sname, scheme, circ in SCHEMES:
        for n in ("9:MB", "8:J12only", "8:MBTSbased"):
            xs = [R["c_pt"].get(f"{sname}|{n}|{pb}") for pb in ("3-5", "5-8", "8-30")]
            if any(xs):
                w(f"| {sname} | {n} | " + " | ".join(pm(x["SF"], x["err"], 2) + ("" if x["S"] >= NWEAK else " *") if x else "–" for x in xs) + " |")
    w("\n## 6. (d) ε_iso w danych w funkcji ΔR do dżetu — uniwersalność\n")
    for var, v in R["d"].items():
        ed = v["edges"]
        lab = [("brak dżetu" if lo >= 50 else f"{lo:g}–{hi:g}" if hi < INF else f"≥ {lo:g}") for lo, hi in zip(ed[:-1], ed[1:])]
        w(f"**{v['desc']}** (`{var}`)\n")
        w("| próbka | " + " | ".join(lab) + " |\n|---|" + "---|" * len(lab))
        for n in DR_MAIN:
            w(f"| {n} | " + " | ".join((pm(r["eps"], r["err"], 2) + f" ({r['S']:.0f})") if np.isfinite(r["eps"]) and r["S"] > 25 else "–" for r in v["cls"][n]) + " |")
        w("| MC fd_nd | " + " | ".join((pm(r["eps"], r["err"], 2) + f" ({r['S']:.0f})") if r["S"] > 25 else "–" for r in v["MC"]) + " |")
        ch = [f"{x} vs {y}: {num(c['chi2'], 1)}/{c['ndf']}" + (f" (p = {c['p']:.2g})" if c["p"] is not None else "")
              for (vv, x, y), c in ((k.split("|"), c) for k, c in R["d_chi2"].items()) if vv == var]
        w("\nχ²/ndf (biny z N_sig > 25 w obu): " + "; ".join(ch) + ".\n")
    w("W nawiasach N_sig (SBS; MC: liczba kandydatów). „brak dżetu” = w zdarzeniu nie ma dżetu kolekcji powyżej progu.")
    w("\n## 7. Zastrzeżenia\n")
    w("- Dżety torowe (tylko naładowane, pT > 0,5 GeV) ≠ dżety kalorymetryczne trygera L1_J12/HLT_j20; pT dżetu torowego ≈ 0,6 pT pełnego, z dużym rozrzutem.")
    w("- Kolekcja O obcina dżet, w którym leży D⁰ (tory w ΔR < 0,4 usunięte): „dżet O w ΔR < 0,8” to reszta dżetu wokół D⁰ albo dżet sąsiedni.")
    w("- MC nie ma trygera; w fd_nd twarde dżety są rzadkie → przy dużych pT dżetu komórki MC puste (kolumna „cof.”). Domknięcie przy dużym „cof.” nie jest dowodem.")
    w("- SBS zakłada kształt tła Δm niezależny od zmiennej wewnątrz grupy (bin pT × pass/fail); przy podziale w/poza dżetem kontrola = fit (kolumna „metoda”).")
    w("- Pary fałszywe w MC (D⁰ i D* z różnych wpisów o tym samym eventNumber; ułamek w tabeli wejść `same_entry`) są odrzucone.")
    open(a.report, "w").write("\n".join(L) + "\n")
    print(f"werdykt: {verdict}\nraport: {a.report}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    sys.exit(main())

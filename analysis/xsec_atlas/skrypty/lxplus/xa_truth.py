#!/usr/bin/env python3
"""XSEC_atlas — prawda generatora mb_nd z ROOT v2.0 (TruthDmesons) z PRAWDZIWYM trybem rozpadu (lista corek), niezaleznie od `decayType`.
Na wiersz prawdziwego D0 / D+ / D*+ (pT >= 2 GeV, |eta| < 2,5): pdg, pt [GeV], eta, phi, decayType, kod topologii (liczby K-, K+, pi+, pi-, pi0, D0 —
znaki wzgledem znaku pdg rodzica), klasa brakujacej masy (0: m_D - m_wid < 30 MeV; 1: m_brak < 100 MeV; 2: K0; 3: eta; 4: inne), fromBottom;
dla D*: to samo dla corki D0 odszukanej wsrod TruthDmesons tego samego zdarzenia (zgodnosc pT, eta corki).
Klasyfikacja topologii i brakujacej masy: funkcja truth_tables z ../bkgcat/bkg_hist.py przepisana 1:1 (tu z dolaczeniem corki D0).
Uzycie: xa_truth.py --rawdir DIR --k K --n N --out OUT.npz"""
import argparse, glob, os, sys, time
import numpy as np

MPI, MK, MPI0 = 139.57039, 493.677, 134.9768


def one(rootfile):
    import awkward as ak
    import uproot
    t = uproot.open(rootfile)["CharmAnalysis"]
    tb = ["pdgId", "decayType", "pt", "eta", "phi", "m", "fromBottom", "daughterPdgId", "daughterPt", "daughterEta", "daughterPhi"]
    A = t.arrays(["TruthDmesons_" + b for b in tb])
    g = lambda b: A["TruthDmesons_" + b]  # noqa: E731
    cnt = ak.to_numpy(ak.num(g("pdgId")))
    entry = np.repeat(np.arange(len(cnt), dtype=np.int64), cnt)
    fl = lambda b: ak.to_numpy(ak.flatten(g(b))).astype(np.float64)  # noqa: E731
    pdg = ak.to_numpy(ak.flatten(g("pdgId"))).astype(np.int64)
    pt, eta, phi, m = fl("pt"), fl("eta"), fl("phi"), fl("m")
    did, dp, de, dph = (ak.flatten(g(b), axis=1) for b in ("daughterPdgId", "daughterPt", "daughterEta", "daughterPhi"))
    sg = ak.Array(np.sign(pdg))
    d = did * sg
    c = lambda mask: np.minimum(ak.to_numpy(ak.sum(mask, axis=1)), 9).astype(np.int64)  # noqa: E731
    nself = (np.abs(pdg) == 421).astype(np.int64)
    isd0 = abs(did) == 421
    nd0 = np.clip(ak.to_numpy(ak.sum(isd0, axis=1)) - nself, 0, 9)
    code = c(d == -321) + 10 * c(d == 321) + 100 * c(d == 211) + 1000 * c(d == -211) + 10000 * c(abs(did) == 111) + 100000 * nd0
    vis = ~isd0
    dm = ak.where(abs(did) == 321, MK, ak.where(abs(did) == 211, MPI, MPI0))
    vp, ve, vph, vm = dp[vis], de[vis], dph[vis], dm[vis]
    px, py, pz = vp * np.cos(vph), vp * np.sin(vph), vp * np.sinh(ve)
    E = np.sqrt(px * px + py * py + pz * pz + vm * vm)
    s = lambda x: ak.to_numpy(ak.sum(x, axis=1)).astype(np.float64)  # noqa: E731
    Ev, pxv, pyv, pzv = s(E), s(px), s(py), s(pz)
    mvis = np.sqrt(np.maximum(Ev ** 2 - pxv ** 2 - pyv ** 2 - pzv ** 2, 0.0))
    Px, Py, Pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(np.clip(eta, -20, 20))
    ED = np.sqrt(Px ** 2 + Py ** 2 + Pz ** 2 + m ** 2)
    mm2 = (ED - Ev) ** 2 - (Px - pxv) ** 2 - (Py - pyv) ** 2 - (Pz - pzv) ** 2
    mmiss = np.sign(mm2) * np.sqrt(np.abs(mm2))
    mc = np.full(len(pdg), 4, np.int64)
    mc[mmiss < 100] = 1
    mc[np.abs(mmiss - 497.611) < 3] = 2
    mc[np.abs(mmiss - 547.862) < 3] = 3
    mc[(m - mvis) < 30] = 0
    # corka D0 mezonu D*: pierwsza corka |pdg| = 421 -> (pt, eta); dopasowanie do wiersza prawdy D0 w tym samym zdarzeniu
    f0 = lambda x: ak.to_numpy(ak.fill_none(ak.firsts(x[isd0], axis=1), -1.0)).astype(np.float64)  # noqa: E731
    d0pt, d0eta = f0(dp), f0(de)
    dty = ak.to_numpy(ak.flatten(g("decayType"))).astype(np.int64)
    d0code = np.full(len(pdg), -1, np.int64); d0mc = np.full(len(pdg), -1, np.int64); d0dty = np.full(len(pdg), -1, np.int64)
    i0 = np.flatnonzero(np.abs(pdg) == 421)
    o = np.argsort(entry[i0], kind="mergesort"); i0 = i0[o]; e0 = entry[i0]
    ids = np.flatnonzero((np.abs(pdg) == 413) & (d0pt > 0))
    lo, hi = np.searchsorted(e0, entry[ids], "left"), np.searchsorted(e0, entry[ids], "right")
    for j, a, b in zip(ids, lo, hi):
        if b > a:
            cand = i0[a:b]
            dd = np.abs(pt[cand] / d0pt[j] - 1) + np.abs(eta[cand] - d0eta[j])
            k = cand[np.argmin(dd)]
            if dd.min() < 1e-3:
                d0code[j], d0mc[j], d0dty[j] = code[k], mc[k], dty[k]
    keep = (pt >= 2000.0) & (np.abs(eta) < 2.5) & np.isin(np.abs(pdg), (421, 411, 413))
    return {"entry": entry[keep], "pdg": pdg[keep].astype(np.int32), "pt": (pt[keep] / 1e3).astype(np.float32), "eta": eta[keep].astype(np.float32),
            "phi": phi[keep].astype(np.float32), "dty": dty[keep].astype(np.int16), "code": code[keep].astype(np.int32), "mc": mc[keep].astype(np.int8),
            "fromB": (ak.to_numpy(ak.flatten(g("fromBottom"))) > 0)[keep].astype(np.int8), "d0code": d0code[keep].astype(np.int32), "d0mc": d0mc[keep].astype(np.int8),
            "d0dty": d0dty[keep].astype(np.int16), "mmiss": mmiss[keep].astype(np.float32)}, len(cnt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rawdir", required=True); ap.add_argument("--k", type=int, default=0); ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--out", required=True); ap.add_argument("--max-files", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    files = [f for i, f in enumerate(sorted(glob.glob(os.path.join(a.rawdir, "*.root")))) if i % a.n == a.k]
    if a.max_files:
        files = files[:a.max_files]
    parts, nev, names = {}, 0, []
    for i, f in enumerate(files):
        T, n = one(f)
        T["fidx"] = np.full(len(T["pdg"]), i * a.n + a.k, np.int16)
        nev += n; names.append(os.path.basename(f))
        for k, v in T.items():
            parts.setdefault(k, []).append(v)
        print(f"[{a.k}/{a.n}] {os.path.basename(f)} zdarzen {n} wierszy {len(T['pdg'])} {time.time() - t0:.0f}s", flush=True)
    out = {k: np.concatenate(v) for k, v in parts.items()}
    out["nev"] = np.array([nev]); out["files"] = np.array(names)
    np.savez_compressed(a.out, **out)
    print("zapisano", a.out, len(out["pdg"]), nev)


if __name__ == "__main__":
    sys.exit(main())

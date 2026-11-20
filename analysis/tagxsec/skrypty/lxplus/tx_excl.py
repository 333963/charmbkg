#!/usr/bin/env python3
"""TAGXSEC (20.11.2026, po uwadze z BKGCAT / CVER) — prawdziwy tryb rozpadu D z córek TruthDmesons (definicja BKGCAT, bkg_hist.truth_tables):
klasa braku mc = 0 (m_D − m_wid < 30 MeV: ekskluzywny), 1 (m_brak < 100 MeV: radiacyjny), 2 (K0), 3 (η), 4 (inne). Tryb sygnału K∓π± / K∓π±π± = decayType 10 / 4 ∧ mc ∈ {0, 1}.
Wiersze: D0/D+ z decayType sygnału, pT ≥ 2,5 GeV, |η| < 2,7 (jak tabela gen w tx_mc.py). Użycie: tx_excl.py <próbka> OUT.npz [nproc=2]"""
import glob, re, sys, time
import numpy as np, uproot, awkward as ak
from multiprocessing import Pool
S, OUT = sys.argv[1:3]; NP = int(sys.argv[3]) if len(sys.argv) > 3 else 2
MK, MPI, MPI0 = 493.677, 139.57039, 134.9768
T0 = time.time()

def one(f):
    fid = int(re.search(r"\._(\d{6})\.", f).group(1))
    t = uproot.open(f)["CharmAnalysis"]; keys = set(t.keys())
    def k(b):
        for c in ("TruthDmesons_" + b, "TruthDmesons_" + b + "_NOSYS"):
            if c in keys: return c
        raise KeyError(b)
    A = t.arrays([k(b) for b in ("pdgId", "decayType", "pt", "eta", "phi", "m", "daughterPdgId", "daughterPt", "daughterEta", "daughterPhi")])
    g = lambda b: A[k(b)]
    cnt = ak.to_numpy(ak.num(g("pdgId"))); entry = np.repeat(np.arange(len(cnt), dtype=np.int64), cnt)
    fl = lambda b: ak.to_numpy(ak.flatten(g(b))).astype(np.float64)
    pdg = ak.to_numpy(ak.flatten(g("pdgId"))).astype(np.int64); pt, eta, phi, m = fl("pt"), fl("eta"), fl("phi"), fl("m")
    dty = ak.to_numpy(ak.flatten(g("decayType"))).astype(np.int64)
    did, dp, de, dph = (ak.flatten(g(b), axis=1) for b in ("daughterPdgId", "daughterPt", "daughterEta", "daughterPhi"))
    vis = abs(did) != 421
    dm = ak.where(abs(did) == 321, MK, ak.where(abs(did) == 211, MPI, MPI0))
    dp, de, dph, dm = dp[vis], de[vis], dph[vis], dm[vis]
    px, py, pz = dp * np.cos(dph), dp * np.sin(dph), dp * np.sinh(de)
    E = np.sqrt(px * px + py * py + pz * pz + dm * dm)
    s = lambda x: ak.to_numpy(ak.sum(x, axis=1)).astype(np.float64)
    Ev, pxv, pyv, pzv = s(E), s(px), s(py), s(pz)
    mvis = np.sqrt(np.maximum(Ev ** 2 - pxv ** 2 - pyv ** 2 - pzv ** 2, 0.0))
    Px, Py, Pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(np.clip(eta, -20, 20))
    ED = np.sqrt(Px ** 2 + Py ** 2 + Pz ** 2 + m ** 2)
    mm2 = (ED - Ev) ** 2 - (Px - pxv) ** 2 - (Py - pyv) ** 2 - (Pz - pzv) ** 2
    mmiss = np.sign(mm2) * np.sqrt(np.abs(mm2))
    mc = np.full(len(pdg), 4, np.int8); mc[mmiss < 100] = 1; mc[np.abs(mmiss - 497.611) < 3] = 2; mc[np.abs(mmiss - 547.862) < 3] = 3; mc[(m - mvis) < 30] = 0
    nd = ak.to_numpy(ak.num(did)).astype(np.int8)
    want = np.where(np.abs(pdg) == 421, 10, 4)
    sel = np.isin(np.abs(pdg), (421, 411)) & (dty == want) & (pt / 1000.0 >= 2.5) & (np.abs(eta) < 2.7)
    return fid, entry[sel], pdg[sel].astype(np.int32), (pt[sel] / 1000.0).astype(np.float32), eta[sel].astype(np.float32), mc[sel], nd[sel], len(cnt)

if __name__ == "__main__":
    fs = sorted(glob.glob(f"/eos/user/k/klyczak/charm_raw_v20/{S}/*.root"))
    o = {k: [] for k in ("fid", "entry", "pdg", "pt", "eta", "mc", "nd")}; nev = {}
    with Pool(NP) as p:
        for r in p.imap(one, fs):
            fid = r[0]; n = len(r[1]); o["fid"].append(np.full(n, fid, np.int32))
            for kk, v in zip(("entry", "pdg", "pt", "eta", "mc", "nd"), r[1:7]): o[kk].append(v)
            nev[fid] = r[7]; print(f"[{time.time()-T0:.0f}s] {fid}: {n} wierszy, {r[7]} zd.", flush=True)
    o = {k: np.concatenate(v) for k, v in o.items()}
    np.savez_compressed(OUT, files=np.array(sorted(nev)), **o)
    for p_, nm in ((421, "D0"), (411, "D+")):
        q = (np.abs(o["pdg"]) == p_) & (o["pt"] >= 3) & (np.abs(o["eta"]) < 2.2)
        print(S, nm, "decayType sygnału, pT≥3, |η|<2,2:", int(q.sum()), "klasy 0..4:", np.bincount(o["mc"][q], minlength=5).tolist(), "ekskl.+rad. =", round(float((o["mc"][q] <= 1).mean()), 4))

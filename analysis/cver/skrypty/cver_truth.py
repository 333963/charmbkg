#!/usr/bin/env python3
"""CVER (20.11.2026, projekt B) — niezalezna weryfikacja definicji C w XSEC: tablica prawdy TruthDmesons mb_nd (ROOT v2.0).
Dla kazdego prawdziwego D0 / D+ z pT >= 3 GeV, |eta| < 2.5: decayType, nDaughters, lista corek (wlasne kodowanie), m - m_wid, m_brak,
oraz kandydaci reco przypisani przez TruthDmesons_recoMatchIndices (INNE lacze niz w BKGCAT, ktore szlo przez DMesons_truthPt).
Tylko MC, bez AFP, bez danych. Uzycie: cver_truth.py OUT.npz plik1.root [plik2.root ...]   (jeden proces)"""
import sys
import collections
import numpy as np
import awkward as ak
import uproot

WANT = {421: 10, 411: 4}
MASS = {321: 493.677, 211: 139.57039, 111: 134.9768, 22: 0.0, 11: 0.511, 13: 105.658, 12: 0.0, 14: 0.0, 16: 0.0,
        310: 497.611, 130: 497.611, 311: 497.611, 221: 547.862, 2212: 938.272, 2112: 939.565}


def one(path):
    t = uproot.open(path)["CharmAnalysis"]
    tb = ["pdgId", "decayType", "pt", "eta", "phi", "m", "fromBottom", "nDaughters", "daughterPdgId", "daughterPt", "daughterEta", "daughterPhi",
          "recoMatchIndices", "recoMatchCategories", "recoMatchMass", "recoMatchPt"]
    A = t.arrays(["TruthDmesons_" + b for b in tb] + ["DMesons_decayType_NOSYS", "DMesons_m_NOSYS", "DMesons_pt_NOSYS", "DMesons_fitOutput__Lxy_NOSYS",
                                                       "DMesons_TruthMatchCategory_NOSYS"])
    g = lambda b: A["TruthDmesons_" + b]  # noqa: E731
    nev = len(A)
    pdg = g("pdgId")
    sel = ((abs(pdg) == 421) | (abs(pdg) == 411)) & (g("pt") >= 3000.0) & (abs(g("eta")) < 2.5)
    cnt_all = collections.Counter(ak.to_numpy(abs(ak.flatten(pdg))).tolist())
    ev = ak.local_index(pdg, axis=0)
    ev = ak.to_numpy(ak.flatten(ak.broadcast_arrays(ev, pdg)[0][sel]))
    f = lambda b: ak.to_numpy(ak.flatten(g(b)[sel]))  # noqa: E731
    P, pt, eta, phi, m = f("pdgId").astype(np.int64), f("pt").astype(np.float64), f("eta").astype(np.float64), f("phi").astype(np.float64), f("m").astype(np.float64)
    did = ak.flatten(g("daughterPdgId")[sel], axis=1)
    dpt, deta, dphi = (ak.values_astype(ak.flatten(g(b)[sel], axis=1), np.float64) for b in ("daughterPt", "daughterEta", "daughterPhi"))
    # lista zawiera takze sam mezon (pdgId = pdgId rodzica) — usuwany; liczba takich wpisow zapisana w `nself`
    notself = did != ak.Array(P)
    nself = ak.to_numpy(ak.sum(~notself, axis=1)).astype(np.int32)
    did, dpt, deta, dphi = did[notself], dpt[notself], deta[notself], dphi[notself]
    sg = ak.Array(np.sign(P))
    d = did * sg
    c = lambda mask: np.minimum(ak.to_numpy(ak.sum(mask, axis=1)), 9).astype(np.int64)  # noqa: E731
    known = (abs(did) == 321) | (abs(did) == 211) | (did == 111)
    # kod: K-(wg D), K+, pi+, pi-, pi0, inne
    code = c(d == -321) + 10 * c(d == 321) + 100 * c(d == 211) + 1000 * c(d == -211) + 10000 * c(did == 111) + 100000 * c(~known)
    oth = collections.Counter(ak.to_numpy(abs(ak.flatten(did[~known]))).tolist())
    dm = ak.zeros_like(dpt)
    for k, v in MASS.items():
        dm = ak.where(abs(did) == k, v, dm)
    px, py, pz = dpt * np.cos(dphi), dpt * np.sin(dphi), dpt * np.sinh(deta)
    E = np.sqrt(px * px + py * py + pz * pz + dm * dm)
    s = lambda x: ak.to_numpy(ak.sum(x, axis=1)).astype(np.float64)  # noqa: E731
    Ev, pxv, pyv, pzv = s(E), s(px), s(py), s(pz)
    mvis = np.sqrt(np.maximum(Ev ** 2 - pxv ** 2 - pyv ** 2 - pzv ** 2, 0.0))
    Px, Py, Pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    ED = np.sqrt(Px ** 2 + Py ** 2 + Pz ** 2 + m ** 2)
    mm2 = (ED - Ev) ** 2 - (Px - pxv) ** 2 - (Py - pyv) ** 2 - (Pz - pzv) ** 2
    emiss = ED - Ev
    want = np.where(np.abs(P) == 421, 10, 4)
    T = dict(pdg=P.astype(np.int32), pt=(pt / 1000).astype(np.float32), eta=eta.astype(np.float32), mD=m.astype(np.float32), dty=f("decayType").astype(np.int32),
             fromB=f("fromBottom").astype(np.int8), nD=f("nDaughters").astype(np.int32), nlist=ak.to_numpy(ak.num(did)).astype(np.int32),
             code=code.astype(np.int32), nself=nself, dmvis=(m - mvis).astype(np.float32), mm2=mm2.astype(np.float64), emiss=emiss.astype(np.float32))
    # reco przypisane przez recoMatchIndices
    idx = ak.flatten(g("recoMatchIndices")[sel], axis=1)
    rcat = ak.flatten(g("recoMatchCategories")[sel], axis=1)
    rmass = ak.flatten(g("recoMatchMass")[sel], axis=1)
    rpt = ak.flatten(g("recoMatchPt")[sel], axis=1)
    nk = ak.to_numpy(ak.num(idx))
    trow = np.repeat(np.arange(len(P)), nk)
    fi = ak.to_numpy(ak.flatten(idx)).astype(np.int64)
    nreco = ak.to_numpy(ak.num(A["DMesons_decayType_NOSYS"]))
    roff = np.concatenate([[0], np.cumsum(nreco)])[:-1]
    evk = ev[trow]
    valid = (fi >= 0) & (fi < nreco[evk])
    gi = np.where(valid, fi + roff[evk], 0)
    fr = lambda b: ak.to_numpy(ak.flatten(A[b]))  # noqa: E731
    rdt, rm, rp, rl, rc = (fr("DMesons_" + b + "_NOSYS") for b in ("decayType", "m", "pt", "fitOutput__Lxy", "TruthMatchCategory"))
    if len(rdt) == 0:
        rdt = np.zeros(1, np.int32); rm = np.zeros(1, np.float32); rp = np.zeros(1, np.float32); rl = np.zeros(1, np.float32); rc = np.zeros(1, np.int32)
    keep = valid & (rdt[gi] == want[trow])
    R = dict(r_trow=trow[keep].astype(np.int64), r_m=(rm[gi][keep] / 1000).astype(np.float32), r_pt=(rp[gi][keep] / 1000).astype(np.float32),
             r_lxy=rl[gi][keep].astype(np.float32), r_cat=rc[gi][keep].astype(np.int32),
             r_cat2=ak.to_numpy(ak.flatten(rcat))[keep].astype(np.int32), r_m2=(ak.to_numpy(ak.flatten(rmass))[keep] / 1000).astype(np.float32))
    st = dict(nev=nev, n_invalid_idx=int((~valid).sum()), n_idx=int(len(fi)))
    return T, R, cnt_all, oth, st


def main():
    out, files = sys.argv[1], sys.argv[2:]
    Ts, Rs, ca, co, off, nev, bad, nidx = [], [], collections.Counter(), collections.Counter(), 0, 0, 0, 0
    for p in files:
        T, R, c1, c2, st = one(p)
        R["r_trow"] = R["r_trow"] + off
        off += len(T["pdg"]); nev += st["nev"]; bad += st["n_invalid_idx"]; nidx += st["n_idx"]
        Ts.append(T); Rs.append(R); ca.update(c1); co.update(c2)
        print("OK", p.split("._")[-1][:6], st, len(T["pdg"]), flush=True)
    o = {k: np.concatenate([t[k] for t in Ts]) for k in Ts[0]}
    o.update({k: np.concatenate([r[k] for r in Rs]) for k in Rs[0]})
    o["nev"] = np.array([nev]); o["nfiles"] = np.array([len(files)]); o["bad_idx"] = np.array([bad, nidx])
    o["all_pdg"] = np.array(sorted(ca.items()), np.int64).reshape(-1, 2); o["other_dau"] = np.array(sorted(co.items()), np.int64).reshape(-1, 2)
    np.savez_compressed(out, **o)
    print("DONE", out, nev, off, flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""p1: zwarta tabela per (zdarzenie, strona) z ROOT v2.0: tory AFP per stacja, protony rek., (MC) proton prawdziwy, (dane) bity L1/HLT.
Bez mas kandydatów D, bez histogramów masy. Tylko odczyt EOS. Limit czasu (deadline) — przetwarza równomiernie rozłożony podzbiór plików.
Użycie: p1_extract.py <próbka> <out.npz> [max_s=480] [nproc=4] [co_ile_plik=1]
Wiersz zapisywany, gdy po stronie jest ≥1 tor AFP, ≥1 proton rek. albo (MC) proton prawdziwy z ξ<0,25. Liczniki wszystkich zdarzeń osobno.
Stacje: 0 = A far, 1 = A near, 2 = C near, 3 = C far. Strona 0 = A, 1 = C."""
import sys, glob, time, numpy as np, uproot, awkward as ak
from multiprocessing import Pool
S = sys.argv[1]; OUT = sys.argv[2]; MAXS = float(sys.argv[3]) if len(sys.argv) > 3 else 480; NP = int(sys.argv[4]) if len(sys.argv) > 4 else 4
STEP = int(sys.argv[5]) if len(sys.argv) > 5 else 1
T0 = time.time(); ISMC = not S.startswith("data")
TRIG = ["trigTBP_L1_AFP_A", "trigTBP_L1_AFP_C", "trigTBP_L1_AFP_A_OR_C", "trigTBP_L1_AFP_FSA_BGRP12", "trigTBP_L1_AFP_NSA_BGRP12",
        "trigTBP_L1_AFP_FSC_BGRP12", "trigTBP_L1_AFP_NSC_BGRP12", "trigTBP_L1_MBTS_1", "trigTBP_L1_MBTS_2",
        "trigPassed_HLT_mb_mbts_L1MBTS_1", "trigPassed_HLT_noalg_mb_L1MBTS_1", "trigPassed_HLT_noalg_L1AFP_A_OR_C", "trigPassed_HLT_noalg_L1J12"]
NEAR = {0: 1, 1: 2}; FAR = {0: 0, 1: 3}

def first(a, n, default=-99.):
    return ak.to_numpy(ak.fill_none(ak.pad_none(a, n, clip=True), default))

def one(fn):
    if time.time() - T0 > MAXS: return None
    try:
        t = uproot.open(fn)["CharmAnalysis"]
    except Exception as e:
        return ("ERR", fn, str(e))
    keys = set(t.keys())
    def k(b):
        for c in (b, b + "_NOSYS", "EventInfo_" + b, "EventInfo_" + b + "_NOSYS"):
            if c in keys: return c
        return None
    br = {x: k(x) for x in ["AfpTracks_station", "AfpTracks_x", "AfpTracks_y", "AfpTracks_nClusters", "AfpProtons_xi", "AfpProtons_side", "AfpProtons_px",
                            "AfpProtons_py", "AfpProtons_chi2", "nPileupVertices", "xi_cal_A", "xi_cal_C", "lumiBlock"]}
    if ISMC:
        br.update({x: k(x) for x in ["TruthProtons_px", "TruthProtons_py", "TruthProtons_pz", "TruthProtons_xi", "TruthProtons_side", "TruthDmesons_pdgId",
                                     "TruthDmesons_pt", "TruthDmesons_eta"]})
    else:
        br.update({x: k(x) for x in TRIG + ["passGRL_MinBiasAFP"]})
    br = {a: b for a, b in br.items() if b}
    a = t.arrays(list(br.values()), library="ak")
    g = lambda x: a[br[x]]
    n = len(a); rows = []
    cnt = {"n_ev": n}
    npu = ak.to_numpy(g("nPileupVertices")); lb = ak.to_numpy(g("lumiBlock")) if "lumiBlock" in br else np.zeros(n)
    xc = [ak.to_numpy(g("xi_cal_A")), ak.to_numpy(g("xi_cal_C"))]
    if ISMC:
        pid = abs(g("TruthDmesons_pdgId")); pt = g("TruthDmesons_pt"); eta = abs(g("TruthDmesons_eta"))
        sc = 1e-3 if ak.max(pt) > 1000 else 1.
        fid = (eta < 2.2) & (pt * sc >= 3) & (pt * sc < 30)
        dfl = ak.to_numpy(ak.any(fid & (pid == 421), axis=1)) * 1 + ak.to_numpy(ak.any(fid & (pid == 411), axis=1)) * 2 + ak.to_numpy(ak.any(fid & (pid == 413), axis=1)) * 4
        fid2 = (eta < 2.2) & (pt * sc >= 2) & (pt * sc < 30)
        dfl = dfl + ak.to_numpy(ak.any(fid2 & ((pid == 421) | (pid == 411) | (pid == 413)), axis=1)) * 8
        extra = dfl; grl = np.ones(n, bool)
    else:
        extra = np.zeros(n, np.int32)
        for i, x in enumerate(TRIG):
            if x in br: extra |= (ak.to_numpy(g(x)).astype(np.int32) != 0) << i
        grl = ak.to_numpy(g("passGRL_MinBiasAFP")) != 0
        cnt["n_grl"] = int(grl.sum()); cnt["n_grl_npu0"] = int((grl & (npu == 0)).sum())
        for i, x in enumerate(TRIG): cnt["grl_" + x] = int((grl & ((extra >> i) & 1 == 1)).sum())
    st = g("AfpTracks_station")
    for s in (0, 1):
        mn = st == NEAR[s]; mf = st == FAR[s]
        nn = ak.to_numpy(ak.sum(mn, axis=1)); nf = ak.to_numpy(ak.sum(mf, axis=1))
        xn = first(g("AfpTracks_x")[mn], 1)[:, 0]; yn = first(g("AfpTracks_y")[mn], 1)[:, 0]; cn = first(g("AfpTracks_nClusters")[mn], 1, -1)[:, 0]
        xf = first(g("AfpTracks_x")[mf], 1)[:, 0]; yf = first(g("AfpTracks_y")[mf], 1)[:, 0]; cf = first(g("AfpTracks_nClusters")[mf], 1, -1)[:, 0]
        mp = g("AfpProtons_side") == s
        npr = ak.to_numpy(ak.sum(mp, axis=1))
        xi = g("AfpProtons_xi")[mp]; o = ak.argsort(xi, axis=1)
        xi2 = first(xi[o], 2, -1.)
        ppx = first(g("AfpProtons_px")[mp][o], 1)[:, 0]; ppy = first(g("AfpProtons_py")[mp][o], 1)[:, 0]; chi = first(g("AfpProtons_chi2")[mp][o], 1)[:, 0]
        nwin = ak.to_numpy(ak.sum((xi >= 0.03) & (xi <= 0.10), axis=1))
        cols = [np.full(n, s), lb, npu, grl * 1, extra, xc[s], xc[1 - s], nn, nf, xn, yn, cn, xf, yf, cf, npr, nwin, xi2[:, 0], xi2[:, 1], ppx, ppy, chi]
        keep = (nn > 0) | (nf > 0) | (npr > 0)
        if ISMC:
            ms = g("TruthProtons_side") == s
            txi = g("TruthProtons_xi")[ms]; o = ak.argsort(txi, axis=1)
            nt = ak.to_numpy(ak.sum(ms, axis=1)); ntw = ak.to_numpy(ak.sum((txi >= 0.03) & (txi <= 0.10), axis=1))
            t1 = first(txi[o], 2, -1.)
            tpx = first(g("TruthProtons_px")[ms][o], 1)[:, 0]; tpy = first(g("TruthProtons_py")[ms][o], 1)[:, 0]; tpz = first(g("TruthProtons_pz")[ms][o], 1)[:, 0]
            cols += [nt, ntw, t1[:, 0], t1[:, 1], tpx, tpy, tpz]
            keep |= (nt > 0) & (t1[:, 0] < 0.25)
            cnt[f"n_truth_any_{s}"] = int((nt > 0).sum())
        cnt[f"n_keep_{s}"] = int(keep.sum())
        rows.append(np.stack([np.asarray(c, np.float32)[keep] for c in cols], axis=1))
    return ("OK", fn, cnt, np.concatenate(rows))

if __name__ == "__main__":
    fs = sorted(glob.glob(f"/eos/user/k/klyczak/charm_raw_v20/{S}/*dpt3to150*.root"))[::STEP]
    # kolejność: równomiernie po liście (bit-reversal-like), żeby przerwanie dało podzbiór rozłożony po runie
    idx = sorted(range(len(fs)), key=lambda i: (i * 0.6180339887) % 1.0); fs = [fs[i] for i in idx]
    tot = {}; R = []; done = []; err = []
    with Pool(NP) as p:
        for r in p.imap_unordered(one, fs):
            if r is None: continue
            if r[0] == "ERR": err.append(r[1:]); continue
            done.append(r[1].split("/")[-1])
            for a, b in r[2].items(): tot[a] = tot.get(a, 0) + b
            R.append(r[3])
    names = "side lb npu grl extra xical_same xical_opp n_near n_far x_near y_near ncl_near x_far y_far ncl_far n_prot n_win xi1 xi2 px py chi2".split()
    if ISMC: names += "nt nt_win txi1 txi2 tpx tpy tpz".split()
    np.savez_compressed(OUT, rows=np.concatenate(R) if R else np.zeros((0, len(names)), np.float32), names=np.array(names), files=np.array(done),
                        cnt_keys=np.array(list(tot.keys())), cnt_vals=np.array(list(tot.values()), np.int64), trig=np.array(TRIG), nfiles_all=len(fs))
    print(S, "pliki", len(done), "/", len(fs), "bledy", len(err), err[:2], "czas %.0f s" % (time.time() - T0)); print(tot)

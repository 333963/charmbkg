#!/usr/bin/env python3
"""TAGXSEC (20.11.2026) — tabele do C_D: zrekonstruowany prawdziwy sygnał (wiersze app + evtable) i prawdziwe D (gen) z krotnością zdarzenia.
  tx_mc.py mb_sd  APPDIR OUT.npz            gen z ROOT v2.0 (TruthDmesons, TruthProtons), złączenie po (plik, wpis)
  tx_mc.py mb_nd  APPDIR OUT.npz GEN.npz    gen z przesłanej tabeli truthD_mb (mod sygnału), krotność z evtable po (plik, wpis)
Sygnał reco = ok & truth∈{0,1} & (D0: swap≠1) & masa w zakresie & 3≤pT<30 — jak licznik C w XSEC (wps_hist.py). s = liczba progów ≤ score (WP70 ⇔ s ≥ 6)."""
import glob, json, os, re, sys, time
import numpy as np, h5py
S, APP, OUT = sys.argv[1:4]
EOS = "/eos/user/k/klyczak"
BINS = [(3, 5), (5, 8), (8, 30)]
SP = {"D0": {"mr": (1.70, 2.05), "tag": "D0.xgb.H.sb.n.s1", "pdg": 421, "dt": 10}, "Dplus": {"mr": (1.70, 2.20), "tag": "Dplus.xgb.H.sb.n.s1", "pdg": 411, "dt": 4}}
W = np.load(f"{EOS}/charm_q/wpscan/wps_hist.npz", allow_pickle=False)
thr = {(sp, b): W[f"{sp}|{c['tag']}|{b[0]}-{b[1]}|thr"] for sp, c in SP.items() for b in BINS}
t0 = time.time()
EVF = ["event_uid", "event_entry", "event_number", "nSelBL", "nTracksPt500Eta25", "afp_tag_A", "afp_tag_C", "nPileupVertices"]
ev = {}                                                    # id pliku → evtable posortowana po wpisie
key2id = {}
for f in sorted(glob.glob(f"{EOS}/charm_h5_v20/parts/{S}/{S}_*_none.h5")):
    fid = int(re.search(r"_(\d{6})_none\.h5$", f).group(1))
    with h5py.File(f, "r") as h:
        e = h["evtable"].fields(EVF)[:]
    e = e[np.argsort(e["event_entry"])]
    assert np.array_equal(e["event_entry"], np.arange(len(e))), f"{f}: wpisy nieciągłe"
    ev[fid] = e
    k = np.unique(e["event_uid"] >> np.uint64(22)); assert len(k) == 1, f
    key2id[int(k[0])] = fid
print(f"[{time.time()-t0:.0f}s] evtable: {len(ev)} plików, {sum(len(e) for e in ev.values())} zdarzeń", flush=True)
# ---- prawda protonu per zdarzenie (tylko mb_sd) i gen
tp = {}
G = {k: [] for k in ("fid", "entry", "pdg", "pt", "eta", "fromB", "mode")}
if S != "mb_nd":
    import uproot, awkward as ak
    for f in sorted(glob.glob(f"{EOS}/charm_raw_v20/{S}/*.root")):
        fid = int(re.search(r"\._(\d{6})\.", f).group(1))
        t = uproot.open(f)["CharmAnalysis"]; keys = set(t.keys())
        def k(b):
            for c in (b, b + "_NOSYS", "EventInfo_" + b, "EventInfo_" + b + "_NOSYS"):
                if c in keys: return c
            raise KeyError(b)
        br = {x: k(x) for x in ["TruthProtons_xi", "TruthProtons_side", "TruthDmesons_pdgId", "TruthDmesons_pt", "TruthDmesons_eta", "TruthDmesons_decayType",
                                "TruthDmesons_fromBottom", "eventNumber"]}
        a = t.arrays(list(br.values()), library="ak"); g = lambda x: a[br[x]]
        n = len(a)
        if fid not in ev or len(ev[fid]) != n or not np.array_equal(ak.to_numpy(g("eventNumber")).astype(np.uint64), ev[fid]["event_number"].astype(np.uint64)):
            print(f"UWAGA {f}: brak zgodności z evtable (fid {fid}, n {n}) — pomijam", flush=True); continue
        xi, sd = g("TruthProtons_xi"), g("TruthProtons_side")
        w = (xi >= 0.03) & (xi <= 0.10)
        tp[fid] = np.stack([ak.to_numpy(ak.sum(w & (sd == 0), axis=1)), ak.to_numpy(ak.sum(w & (sd == 1), axis=1)),
                            ak.to_numpy(ak.sum(sd == 0, axis=1)), ak.to_numpy(ak.sum(sd == 1, axis=1))], 1).astype(np.int8)
        pid = g("TruthDmesons_pdgId"); pt = g("TruthDmesons_pt") / 1000.0; eta = g("TruthDmesons_eta"); dt = g("TruthDmesons_decayType")
        sel = ((abs(pid) == 421) | (abs(pid) == 411)) & (pt >= 2.5) & (abs(eta) < 2.7)
        want = ak.where(abs(pid) == 421, 10, 4)
        cnt = ak.to_numpy(ak.num(pid)); evi = np.repeat(np.arange(n), cnt); fs = ak.to_numpy(ak.flatten(sel))
        G["fid"].append(np.full(fs.sum(), fid)); G["entry"].append(evi[fs]); G["pdg"].append(ak.to_numpy(ak.flatten(pid))[fs]); G["pt"].append(ak.to_numpy(ak.flatten(pt))[fs])
        G["eta"].append(ak.to_numpy(ak.flatten(eta))[fs]); G["fromB"].append(ak.to_numpy(ak.flatten(g("TruthDmesons_fromBottom")))[fs])
        G["mode"].append(ak.to_numpy(ak.flatten(dt == want))[fs])
        print(f"[{time.time()-t0:.0f}s] {os.path.basename(f)}: {n} zd., D {int(fs.sum())}", flush=True)
    G = {k: np.concatenate(v) for k, v in G.items()}
else:
    z = np.load(sys.argv[4], allow_pickle=False)
    G = {"fid": z["tag"].astype(np.int64), "entry": z["entry"], "pdg": z["pdg"], "pt": z["pt"], "eta": z["eta"], "fromB": z["fromB"], "mode": z["mode"]}
    evn = z["evnum"]
def evcols(fid, entry, evnum=None):
    n = len(fid); o = {"nsel": np.zeros(n, np.int32), "ntrk": np.zeros(n, np.int32), "tagA": np.zeros(n, np.int8), "tagC": np.zeros(n, np.int8), "npu": np.zeros(n, np.int16),
                       "pA": np.full(n, -1, np.int8), "pC": np.full(n, -1, np.int8), "nA": np.full(n, -1, np.int8), "nC": np.full(n, -1, np.int8), "okf": np.zeros(n, bool)}
    for f in np.unique(fid):
        m = fid == f
        if int(f) not in ev: continue
        e = ev[int(f)]; ix = entry[m]
        good = ix < len(e)
        if evnum is not None:
            good &= e["event_number"][np.minimum(ix, len(e) - 1)].astype(np.uint64) == evnum[m].astype(np.uint64)
        r = e[np.minimum(ix, len(e) - 1)]
        o["nsel"][m] = r["nSelBL"]; o["ntrk"][m] = r["nTracksPt500Eta25"]; o["tagA"][m] = r["afp_tag_A"] > 0; o["tagC"][m] = r["afp_tag_C"] > 0; o["npu"][m] = r["nPileupVertices"]
        okf = good.copy()
        if tp:
            if int(f) in tp:
                q = tp[int(f)][np.minimum(ix, len(e) - 1)]; o["pA"][m], o["pC"][m], o["nA"][m], o["nC"][m] = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
            else:
                okf[:] = False
        o["okf"][m] = okf
    return o
out = {}
c = evcols(G["fid"], G["entry"], evn if S == "mb_nd" else None)
for k, v in {**G, **c}.items(): out["g_" + k] = v
print(f"[{time.time()-t0:.0f}s] gen: {len(G['pt'])} wierszy, złączone {int(c['okf'].sum())}", flush=True)
# ---- reco
R = {k: [] for k in ("sp", "m", "pt", "lg", "s", "truth", "fid", "entry")}
nf = 0
for f in sorted(glob.glob(os.path.join(APP, f"app_{S}_*.npz"))):
    z = np.load(f, allow_pickle=False); tags = json.loads(str(z["meta"]))["tags"]; nf += 1
    for isp, (sp, cc) in enumerate(SP.items()):
        if f"{sp}__m" not in z.files: continue
        jt = tags[sp].index(cc["tag"])
        uid, m, pt, ok, truth, swap = (z[f"{sp}__{x}"] for x in ("uid", "m", "pt", "ok", "truth", "swap"))
        lg = z[f"{sp}__lg"][:, jt].astype(np.float64)
        sel = ok & np.isin(truth, (0, 1)) & (m >= cc["mr"][0]) & (m < cc["mr"][1]) & (pt >= 3) & (pt < 30)
        if sp == "D0": sel &= swap != 1
        uid, m, pt, lg, truth = uid[sel], m[sel], pt[sel], lg[sel], truth[sel]
        s = np.zeros(len(m), np.int8)
        for b in BINS:
            inb = (pt >= b[0]) & (pt < b[1]); s[inb] = np.searchsorted(thr[(sp, b)], lg[inb], side="right")
        key = (uid >> np.uint64(22)).astype(np.int64); fid = np.array([key2id[int(x)] for x in key], np.int64)
        ent = np.zeros(len(uid), np.int64)
        for ff in np.unique(fid):
            mm = fid == ff; e = ev[int(ff)]; o = np.argsort(e["event_uid"]); ix = o[np.searchsorted(e["event_uid"][o], uid[mm])]
            assert np.array_equal(e["event_uid"][ix], uid[mm]); ent[mm] = e["event_entry"][ix]
        for k, v in (("sp", np.full(len(m), isp, np.int8)), ("m", m), ("pt", pt), ("lg", lg.astype(np.float32)), ("s", s), ("truth", truth), ("fid", fid), ("entry", ent)):
            R[k].append(v)
R = {k: np.concatenate(v) for k, v in R.items()}
c = evcols(R["fid"], R["entry"])
for k, v in {**R, **c}.items(): out["r_" + k] = v
out["meta"] = np.array(json.dumps({"sample": S, "app_files": nf, "h5_files": len(ev), "root_files": len(tp), "n_events": int(sum(len(e) for e in ev.values()))}))
np.savez_compressed(OUT, **out)
for isp, sp in enumerate(SP):
    for b in BINS:
        q = (R["sp"] == isp) & (R["pt"] >= b[0]) & (R["pt"] < b[1])
        g = (abs(G["pdg"]) == SP[sp]["pdg"]) & (G["mode"] == 1) & (G["pt"] >= b[0]) & (G["pt"] < b[1]) & (abs(G["eta"]) < 2.2)
        print(f"{S} {sp} {b}: reco sygnał {int(q.sum())}, WP70 {int((q & (R['s'] >= 6)).sum())}, gen (mod sygnału) {int(g.sum())}")
print(f"[{time.time()-t0:.0f}s] zapisano {OUT}", flush=True)

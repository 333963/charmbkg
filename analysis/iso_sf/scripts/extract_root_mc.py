"""ISO_sf: MC sygnał D*→D0(Kπ)π z lokalnych ntupli v9 (mc23 minbias ND Dpt3to150, ntuple produkcji v9) — mają n_trk (nTracksPt500Eta25),
tory (lead_pt) i ptcone. Wiersz D* (decayType 1) z TruthMatchCategory ∈ {0,1}; izolacja i pT z wiersza D0 (D0Index, decayType 10).
Wyjście: out/root_mc_dstar.parquet (jednostki GeV)."""
import glob, os, sys
import awkward as ak, numpy as np, pandas as pd, uproot

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = "/home/k/Documents/charm/mc23_minbiasND_Dpt3to150GeV_v9"
BR = ["EventInfo_nTracksPt500Eta25_NOSYS", "DMesons_decayType_NOSYS", "DMesons_D0Index_NOSYS", "DMesons_TruthMatchCategory_NOSYS",
      "DMesons_pt_NOSYS", "DMesons_m_NOSYS", "DMesons_DeltaMass_NOSYS", "DMesons_ptcone40_NOSYS", "DMesons_ptcone100_NOSYS",
      "DMesons_daughterInfo__trackId_NOSYS", "DMesons_fitOutput__Lxy_NOSYS", "EventInfo_weight_mc_NOSYS"]


def one(path, fi):
    f = uproot.open(path)
    a = f["CharmAnalysis"].arrays(BR)
    tr = f["Tracks"].arrays(["Tracks_pt_NOSYS", "Tracks_eta_NOSYS"])
    lead = ak.max(ak.where(abs(tr.Tracks_eta_NOSYS) < 2.5, tr.Tracks_pt_NOSYS, 0.0), axis=1)
    lead = np.asarray(ak.fill_none(lead, 0.0)) / 1e3
    ncand = np.asarray(ak.num(a.DMesons_decayType_NOSYS))
    tid = ak.flatten(a.DMesons_daughterInfo__trackId_NOSYS, axis=2)
    ntrkid = np.array([len(set(x)) for x in ak.to_list(tid)], dtype=np.int32)
    dt, cat = a.DMesons_decayType_NOSYS, a.DMesons_TruthMatchCategory_NOSYS
    sel = (dt == 1)
    ev = ak.broadcast_arrays(ak.local_index(dt, axis=0), dt)[0]
    d0 = a.DMesons_D0Index_NOSYS
    E = np.asarray(ak.flatten(ev[sel])); I0 = np.asarray(ak.flatten(d0[sel]))
    out = pd.DataFrame({"file": fi, "ev": E, "cat": np.asarray(ak.flatten(cat[sel])),
                        "dm": np.asarray(ak.flatten(a.DMesons_DeltaMass_NOSYS[sel])) / 1e3,
                        "pt_dstar": np.asarray(ak.flatten(a.DMesons_pt_NOSYS[sel])) / 1e3})
    # wiersz D0
    off = np.r_[0, np.cumsum(ncand)[:-1]]
    g = off[E] + I0
    flat = lambda b: np.asarray(ak.flatten(a[b]))
    out["dt0"] = flat("DMesons_decayType_NOSYS")[g]
    out["pt"] = flat("DMesons_pt_NOSYS")[g] / 1e3
    out["m"] = flat("DMesons_m_NOSYS")[g] / 1e3
    out["cat0"] = flat("DMesons_TruthMatchCategory_NOSYS")[g]
    out["ptcone40"] = flat("DMesons_ptcone40_NOSYS")[g] / 1e3
    out["ptcone100"] = flat("DMesons_ptcone100_NOSYS")[g] / 1e3
    out["Lxy"] = flat("DMesons_fitOutput__Lxy_NOSYS")[g]
    out["n_trk"] = np.asarray(a.EventInfo_nTracksPt500Eta25_NOSYS)[E]
    out["lead_pt"] = lead[E]
    out["n_cand"], out["n_trkid"] = ncand[E], ntrkid[E]
    out["w_mc"] = np.asarray(a.EventInfo_weight_mc_NOSYS)[E]
    out = out[(out.dm < 0.170) & (out.pt >= 3) & (out.pt < 30) & out.cat.isin([0, 1])]
    # inkluzywne prawdziwe D0 (decayType 10, kat. 0/1) — kontrola zależności od n_trk z większą statystyką
    s0 = (dt == 10) & ((cat == 0) | (cat == 1)) & (a.DMesons_pt_NOSYS >= 3000) & (a.DMesons_pt_NOSYS < 30000)
    E0 = np.asarray(ak.flatten(ev[s0]))
    inc = pd.DataFrame({"file": fi, "ev": E0, "pt": np.asarray(ak.flatten(a.DMesons_pt_NOSYS[s0])) / 1e3,
                        "ptcone40": np.asarray(ak.flatten(a.DMesons_ptcone40_NOSYS[s0])) / 1e3,
                        "ptcone100": np.asarray(ak.flatten(a.DMesons_ptcone100_NOSYS[s0])) / 1e3})
    inc["n_trk"] = np.asarray(a.EventInfo_nTracksPt500Eta25_NOSYS)[E0]; inc["lead_pt"] = lead[E0]
    inc["n_cand"], inc["n_trkid"] = ncand[E0], ntrkid[E0]
    return out, inc


if __name__ == "__main__":
    fs = sorted(glob.glob(f"{SRC}/*.root"))
    parts, incs = [], []
    for i, p in enumerate(fs):
        o, n = one(p, i); parts.append(o); incs.append(n); print(i, os.path.basename(p), len(o), len(n), flush=True)
    pd.concat(parts, ignore_index=True).to_parquet(os.path.join(BASE, "out", "root_mc_dstar.parquet"))
    pd.concat(incs, ignore_index=True).to_parquet(os.path.join(BASE, "out", "root_mc_d0incl.parquet"))

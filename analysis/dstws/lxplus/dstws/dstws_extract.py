#!/usr/bin/env python3
"""DSTWS ekstrakcja (projekt B, 17.11.2026): tło D* przez konstrukcję — kandydat D0→Kπ (h5 v2.0, pełna selekcja wierzchołka, score H)
× tor z drzewa `Tracks` ROOT v2.0 → kombinacje RS / WS / ROT / MIX (dstws_kernel) w oknie Δm < --dm-max.

D0: jak D0 par B_trig (`btrig_extract.py`): decayType D0, L_xy ≥ --lxy-min (0,25 mm), |m(Kπ) − m0| < 75 MeV, score H modelem nominalnym
(models/D0_H, charmtag z ev_dstar/src) — ale WSZYSTKIE takie D0 (B_trig trzyma tylko D0 z nominalnym partnerem D* oraz 30 % „all_*”).
Pary nominalne (odniesienie RS): D* z ntupla jak w B_trig (0 < Δm < 0,165, L_xy(D*) ≥ lxy-min, zgodność m/pT/η z D0), złączenie po
event_entry (wpis drzewa; B_trig łączył po (run, event) — w MC niejednoznaczne). Tablica n_* = wszystkie pary nominalne z flagą n_found.
Wyjście npz: tablica D0 (d_*; tylko D0 z ≥ 1 zapisaną kombinacją), kombinacje (c_*; Δm < --row-dm), pary nominalne (n_*),
histogramy H_dm[typ, bin pT, okno m, score/0,02, Δm/0,25 MeV] dla Δm < --dm-max, H_d0[bin pT, okno m, score] (mianowniki), meta (json).
Bity d_trig: 1 J12, 2 MBTS_ANY, 4 mb_mbts_L1MBTS_1, 8 sptrk_pt4_L1MBTS_2, 16 j20_L1MBTS_2, 32 GRL.
Użycie: dstws_extract.py --sample S --models DIR --h5dir DIR --rawdir DIR --part k --nparts N --out OUT.npz [--max-files n]
"""
import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dstws_kernel as K  # noqa: E402

M_WIN = 0.075
NOM_WIN = 0.055
DM_MAX_NOM = 0.165
PT_BINS = [3.0, 5.0, 8.0, 30.0]
SC_N = 50
DM_E = np.arange(0.1375, 0.1700001, 0.00025)
MBTS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_noalg_mb_L1MBTS_1",
        "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1"]
EXTRA = ["event_entry", "event_number", "lumiBlock", "nSelBL", "passGRL_MinBiasAFP", "Dmeson_hypo_swap", "Dmeson_eta", "Dmeson_phi", "Dmeson_pt",
         "Dmeson_mass", "Dmeson_delta_mass", "Dmeson_truth_category", "Dmeson_soft_pion_d0", "Dmeson_Lxy", "trigPassed_HLT_noalg_L1J12"] + \
        ["trigPassed_" + c for c in MBTS] + \
        [f"Dmeson_daughter{i}_{v}" for i in (0, 1, 2) for v in ("pt", "eta", "phi", "charge", "pdgId", "truthPdgId", "truthParentPdgId")]


def read_h5(path, cols, dt_col, keep_dt, chunk, max_rows=0):
    import h5py
    parts = []
    with h5py.File(path, "r") as f:
        ds = f["events"]
        cols = [c for c in cols if c in ds.dtype.names]
        n = ds.shape[0] if not max_rows else min(ds.shape[0], max_rows)
        for a in range(0, n, chunk):
            b = min(n, a + chunk)
            x = ds.fields(cols)[a:b]
            k = np.isin(x[dt_col], keep_dt)
            df = pd.DataFrame({c: x[c][k] for c in cols})
            df["_local_row"] = np.arange(a, b)[k]
            parts.append(df)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def csr(idx, n):
    """Posortowane idx (0..n−1) → offsety CSR."""
    return np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=n))]).astype(np.int64)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--models", required=True)
    ap.add_argument("--h5dir", required=True); ap.add_argument("--rawdir", required=True)
    ap.add_argument("--part", type=int, default=0); ap.add_argument("--nparts", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-files", type=int, default=0, help="test: tylko pierwsze n plików części")
    ap.add_argument("--max-rows", type=int, default=0, help="test: tylko pierwsze N wierszy h5")
    ap.add_argument("--lxy-min", type=float, default=0.25)
    ap.add_argument("--dm-max", type=float, default=0.170); ap.add_argument("--row-dm", type=float, default=0.160)
    ap.add_argument("--h5-chunk", type=int, default=400_000); ap.add_argument("--step", type=int, default=20_000)
    a = ap.parse_args(argv)
    import awkward as ak
    import uproot
    from charmtag import features as F
    from charmtag.io import resolve_columns
    from charmtag.model import Ensemble
    from charmtag.samples import mark_swaps

    t0 = time.time()
    cfg = yaml.safe_load(open(os.path.join(a.models, "config_used.yaml")))
    meta_m = json.load(open(os.path.join(a.models, "meta.json")))
    sp = cfg["species"]["D0"]; dc = cfg["dstar_control"]; m0 = float(sp["m0"]); units = cfg.get("units")
    DT0, DT1 = int(sp["decay_type"]), int(dc["dstar_decay_type"])
    ens = Ensemble(a.models)
    is_mc = not a.sample.startswith("data_")
    files_all = sorted(os.path.basename(p) for p in glob.glob(os.path.join(a.h5dir, f"{a.sample}_*_none.h5")))
    mine = files_all[a.part::a.nparts]
    if a.max_files:
        mine = mine[:a.max_files]
    print(f"[{a.sample}] część {a.part}/{a.nparts}: {len(mine)} z {len(files_all)} plików h5, numba={K.HAVE_NUMBA}", flush=True)
    H_dm = np.zeros((4, 3, 2, SC_N, len(DM_E) - 1), np.int32)
    H_d0 = np.zeros((3, 2, SC_N), np.int64)
    Dp, Cp, Np = [], [], []
    nd_tot = 0
    stat = {"files": [], "per_file": {}, "n_event_mismatch": 0}
    for fn in mine:
        tf = time.time()
        mm = re.match(r"(.+)_(\d{6})_none\.h5$", fn)
        tag = mm.group(2)
        rf = [f for f in glob.glob(os.path.join(a.rawdir, f"*._{tag}.*.root")) if "dpt3to150" in f]
        if len(rf) != 1:
            raise RuntimeError(f"{a.sample} {tag}: plików ROOT dpt3to150 {len(rf)} (oczekiwany 1)")
        path = os.path.join(a.h5dir, fn)
        import h5py
        with h5py.File(path, "r") as f:
            cols_all = list(f["events"].dtype.names)
        cm2 = resolve_columns(cols_all, cfg["columns"], 2); cm3 = resolve_columns(cols_all, cfg["columns"], 3)
        dt_col = cm3.get("decay_type")
        need = sorted({c for c in (F.raw_columns(cm3, 3) + [dt_col] + EXTRA) if c})
        raw = read_h5(path, need, dt_col, [DT0, DT1], a.h5_chunk, a.max_rows)
        stat["files"].append(fn)
        fs = {"n_d0_sel": 0, "n_d0_win": 0, "n_nom": 0, "n_comb": 0}
        stat["per_file"][fn] = fs
        r0 = raw[raw[dt_col] == DT0].reset_index(drop=True)
        r1 = raw[raw[dt_col] == DT1].reset_index(drop=True)
        del raw
        if len(r0) == 0:
            continue
        L = F.to_logical(r0, cm2, 2)
        L["weight"] = L["weight"].fillna(1.0)
        ok = L["Lxy"].to_numpy(float) >= a.lxy_min
        L, r0 = L[ok].reset_index(drop=True), r0[ok].reset_index(drop=True)
        if len(L) == 0:
            continue
        L = F.compute(L, sp, units)
        mark_swaps(L, sp, cfg)
        L["score"], L["fold"] = ens.score(L)
        fs["n_d0_sel"] = int(len(L))
        m = L["m"].to_numpy(float)
        win = np.abs(m - m0) < M_WIN
        L, r0 = L[win].reset_index(drop=True), r0[win].reset_index(drop=True)
        nd = len(L); fs["n_d0_win"] = int(nd)
        if nd == 0:
            continue
        g = lambda c: r0[c].to_numpy()  # noqa: E731
        u = 1e-3 if np.nanmedian(g("Dmeson_daughter0_pt")) > 100 else 1.0
        ent = g("event_entry").astype(np.int64)
        kis0 = np.abs(g("Dmeson_daughter0_pdgId")) == 321
        pick = lambda v, k0: np.where(k0, g(f"Dmeson_daughter0_{v}"), g(f"Dmeson_daughter1_{v}")).astype(np.float64)  # noqa: E731
        d_pt = L["pt"].to_numpy(float); d_m = L["m"].to_numpy(float); d_sc = L["score"].to_numpy(float)
        if np.nanmax(np.abs(d_pt / g("Dmeson_pt") - 1)) > 1e-3 or np.nanmax(np.abs(d_m - g("Dmeson_mass"))) > 1e-4:
            raise RuntimeError(f"{fn}: pT/m logiczne ≠ Dmeson_pt/Dmeson_mass (jednostki?)")
        DD = np.column_stack([d_pt, g("Dmeson_eta"), g("Dmeson_phi"), d_m, pick("charge", kis0),
                              pick("pt", kis0) * u, pick("eta", kis0), pick("phi", kis0),
                              pick("pt", ~kis0) * u, pick("eta", ~kis0), pick("phi", ~kis0)]).astype(np.float64)
        tb = lambda c: (g("trigPassed_" + c) == 1) if "trigPassed_" + c in r0 else np.zeros(nd, bool)  # noqa: E731
        mb_any = np.zeros(nd, bool)
        for c in MBTS:
            mb_any |= tb(c)
        trig = (tb("HLT_noalg_L1J12") * 1 + mb_any * 2 + tb("HLT_mb_mbts_L1MBTS_1") * 4 + tb("HLT_mb_sptrk_pt4_L1MBTS_2") * 8
                + tb("HLT_j20_L1MBTS_2") * 16 + (g("passGRL_MinBiasAFP") == 1) * 32).astype(np.uint8)
        ib = np.digitize(d_pt, PT_BINS) - 1
        iw = (np.abs(d_m - m0) < NOM_WIN).astype(np.int64)
        isc = np.clip((d_sc * SC_N).astype(np.int64), 0, SC_N - 1)
        okb = (ib >= 0) & (ib < 3)
        np.add.at(H_d0, (ib[okb], iw[okb], isc[okb]), 1)
        # --- pary nominalne (jak B_trig, klucz = wpis drzewa)
        NPa = np.zeros((0, 3)); n_d0 = np.zeros(0, np.int64); nom = None
        TPe = np.zeros(0, np.int64); TPa = np.zeros((0, 3))
        if len(r1):
            ds = F.to_logical(r1, cm3, 3)
            e1 = r1["event_entry"].to_numpy(np.int64)
            dsel = ((ds["dm"] > 0) & (ds["dm"] < DM_MAX_NOM) & (ds["Lxy"].to_numpy(float) >= a.lxy_min)).to_numpy()
            dd = pd.DataFrame({"entry": e1, "m0_from": (ds["m"] - ds["dm"]).to_numpy(float), "pt_ds": ds["pt"].to_numpy(float),
                               "eta_ds": ds["eta"].to_numpy(float), "dm": ds["dm"].to_numpy(float), "truth_ds": ds["truth"].to_numpy(float),
                               "row_ds": r1["_local_row"].to_numpy(), "p_pt": r1["Dmeson_daughter2_pt"].to_numpy(float) * u,
                               "p_eta": r1["Dmeson_daughter2_eta"].to_numpy(float), "p_phi": r1["Dmeson_daughter2_phi"].to_numpy(float),
                               "p_q": r1["Dmeson_daughter2_charge"].to_numpy(float), "sp_d0": r1["Dmeson_soft_pion_d0"].to_numpy(float)})[dsel]
            Lk = pd.DataFrame({"entry": ent, "m": d_m, "pt": d_pt, "eta": L["eta"].to_numpy(float), "_i": np.arange(nd)})
            j = Lk.merge(dd, on="entry")
            okj = (np.abs(j["m"] - j["m0_from"]) < dc["match_tol_m"]) & (j["pt"] <= j["pt_ds"] * 1.001) & (j["pt"] > 0.8 * j["pt_ds"]) & \
                  (np.abs(j["eta"] - j["eta_ds"]) < 0.1)
            nom = j[okj].sort_values("_i", kind="stable").reset_index(drop=True)
            n_d0 = nom["_i"].to_numpy(np.int64)
            NPa = nom[["p_pt", "p_eta", "p_phi"]].to_numpy(np.float64)
            if is_mc:
                tp = (np.abs(r1["Dmeson_daughter2_truthPdgId"].to_numpy()) == 211) & (np.abs(r1["Dmeson_daughter2_truthParentPdgId"].to_numpy()) == 413)
                T = pd.DataFrame({"entry": e1[tp], "pt": (r1["Dmeson_daughter2_pt"].to_numpy(float) * u)[tp],
                                  "eta": r1["Dmeson_daughter2_eta"].to_numpy(float)[tp], "phi": r1["Dmeson_daughter2_phi"].to_numpy(float)[tp]})
                T = T.drop_duplicates().sort_values("entry", kind="stable")
                TPe = T["entry"].to_numpy(np.int64); TPa = T[["pt", "eta", "phi"]].to_numpy(np.float64)
        fs["n_nom"] = int(len(n_d0))
        noff = csr(n_d0, nd)
        n_found = np.zeros(len(n_d0), np.int8); n_dm_own = np.full(len(n_d0), np.nan, np.float32)
        nmatch = np.zeros(nd, np.int8); ntrk_d = np.zeros(nd, np.int32)
        cc = []
        with uproot.open(rf[0]) as R:
            T_, C_ = R["Tracks"], R["CharmAnalysis"]
            if T_.num_entries != C_.num_entries:
                raise RuntimeError(f"{rf[0]}: Tracks {T_.num_entries} ≠ CharmAnalysis {C_.num_entries}")
            if ent.max() >= T_.num_entries:
                raise RuntimeError(f"{rf[0]}: wpis {ent.max()} ≥ {T_.num_entries}")
            order = np.argsort(ent, kind="stable"); ents = ent[order]; ue = np.unique(ents)
            for s in range(0, T_.num_entries, a.step):
                e = min(T_.num_entries, s + a.step)
                j0, j1 = np.searchsorted(ue, [s, e])
                if j1 == j0:
                    continue
                uec = ue[j0:j1]
                tr = T_.arrays(["Tracks_pt_NOSYS", "Tracks_eta_NOSYS", "Tracks_phi_NOSYS", "Tracks_charge_NOSYS"], entry_start=s, entry_stop=e)[uec - s]
                cnt = ak.to_numpy(ak.num(tr["Tracks_pt_NOSYS"]))
                off = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64)
                tpt = ak.to_numpy(ak.flatten(tr["Tracks_pt_NOSYS"])).astype(np.float64) * 1e-3
                teta = ak.to_numpy(ak.flatten(tr["Tracks_eta_NOSYS"])).astype(np.float64)
                tphi = ak.to_numpy(ak.flatten(tr["Tracks_phi_NOSYS"])).astype(np.float64)
                tq = ak.to_numpy(ak.flatten(tr["Tracks_charge_NOSYS"])).astype(np.float64)
                del tr
                evn = C_.arrays(["EventInfo_eventNumber"], entry_start=s, entry_stop=e, library="np")["EventInfo_eventNumber"][uec - s]
                p0, p1 = np.searchsorted(ents, [s, e])
                pi = order[p0:p1]                                        # D0 tej porcji (indeksy pliku), rosnąco po wpisie
                pev = np.searchsorted(uec, ent[pi])
                if not is_mc:   # MC: event_number w h5 ma offset job << 32 → porównanie tylko w danych (MC: młodsze 32 bity)
                    stat["n_event_mismatch"] += int(np.sum(evn[pev].astype(np.int64) != g("event_number")[pi].astype(np.int64)))
                else:
                    stat["n_event_mismatch"] += int(np.sum((evn[pev].astype(np.int64) & 0xFFFFFFFF) != (g("event_number")[pi].astype(np.int64) & 0xFFFFFFFF)))
                # CSR par nominalnych i prawdziwych pionów w kolejności porcji
                cn = noff[pi + 1] - noff[pi]
                noff_c = np.concatenate([[0], np.cumsum(cn)]).astype(np.int64)
                nidx = (np.repeat(noff[pi] - noff_c[:-1], cn) + np.arange(int(cn.sum()))).astype(np.int64)
                NP_c = NPa[nidx] if len(nidx) else np.zeros((0, 3))
                tsel = np.nonzero(np.isin(TPe, uec))[0]                 # prawdziwe wolne piony w zdarzeniach tej porcji (TPe posortowane)
                toff_c = csr(np.searchsorted(uec, TPe[tsel]), len(uec)) if len(tsel) else np.zeros(len(uec) + 1, np.int64)
                TP_c = TPa[tsel] if len(tsel) else np.zeros((0, 3))
                pm = K.mix_partners(cnt)
                c_d0, typ, kk, dm, ppt, peta, dr, fl, cnom, nm = K.run(tpt, teta, tphi, tq, off, pev, DD[pi], noff_c, NP_c, toff_c, TP_c, pm,
                                                                         dm_max=a.dm_max)
                nmatch[pi] = nm; ntrk_d[pi] = cnt[pev]
                gi = pi[c_d0]                                           # indeks D0 w pliku
                hb = okb[gi]
                idm = np.clip(np.digitize(dm, DM_E) - 1, 0, len(DM_E) - 2)
                inr = (dm >= DM_E[0]) & (dm < DM_E[-1]) & hb
                np.add.at(H_dm, (typ[inr].astype(np.int64), ib[gi][inr], iw[gi][inr], isc[gi][inr], idm[inr]), 1)
                hasn = cnom >= 0
                gn = nidx[cnom[hasn]]
                n_found[gn] = np.where(typ[hasn] == K.T_RS, 1, 2); n_dm_own[gn] = dm[hasn]
                dmn = np.full(len(dm), np.nan, np.float32); tds = np.full(len(dm), -1, np.int16)
                if nom is not None and hasn.any():
                    dmn[hasn] = nom["dm"].to_numpy(np.float32)[gn]; tds[hasn] = nom["truth_ds"].to_numpy(np.float32)[gn].astype(np.int16)
                kr = dm < a.row_dm
                cc.append(pd.DataFrame({"d0": gi[kr], "typ": typ[kr], "k": kk[kr], "dm": dm[kr], "pis_pt": ppt[kr], "pis_eta": peta[kr],
                                        "dr": dr[kr], "flag": fl[kr], "dm_nom": dmn[kr], "truth_ds": tds[kr]}))
        Cf = pd.concat(cc, ignore_index=True) if cc else pd.DataFrame({"d0": np.zeros(0, np.int64)})
        fs["n_comb"] = int(len(Cf))
        fs["frac_daughters_matched"] = float(np.mean(nmatch == 2))
        used = np.unique(Cf["d0"].to_numpy(np.int64)) if len(Cf) else np.zeros(0, np.int64)
        remap = np.full(nd, -1, np.int64); remap[used] = np.arange(len(used)) + nd_tot
        if len(Cf):
            Cf["d0"] = remap[Cf["d0"].to_numpy(np.int64)]
            Cp.append(Cf)
        hs = g("Dmeson_hypo_swap") if "Dmeson_hypo_swap" in r0 else np.full(nd, -1)
        Dp.append(pd.DataFrame({"file": np.full(len(used), int(tag), np.int32), "row": r0["_local_row"].to_numpy()[used], "entry": ent[used],
                                "run": L["run"].to_numpy()[used], "event": L["event"].to_numpy()[used], "m": d_m[used], "pt": d_pt[used],
                                "eta": DD[used, 1], "score": d_sc[used], "truth": L["truth"].to_numpy(float)[used], "weight": L["weight"].to_numpy(float)[used],
                                "hypo_swap": hs[used], "qK": DD[used, 4], "trig": trig[used], "lb": g("lumiBlock")[used], "nSelBL": g("nSelBL")[used],
                                "ntrk": ntrk_d[used], "nmatch": nmatch[used]}))
        nd_tot += len(used)
        if nom is not None and len(nom):
            i_ = n_d0
            Np.append(pd.DataFrame({"file": np.full(len(i_), int(tag), np.int32), "row": r0["_local_row"].to_numpy()[i_], "row_ds": nom["row_ds"].to_numpy(),
                                    "m": d_m[i_], "pt": d_pt[i_], "score": d_sc[i_], "truth": L["truth"].to_numpy(float)[i_], "trig": trig[i_],
                                    "dm": nom["dm"].to_numpy(), "pis_pt": nom["p_pt"].to_numpy(), "pis_q_rs": (nom["p_q"].to_numpy() * DD[i_, 4] < 0),
                                    "truth_ds": nom["truth_ds"].to_numpy(), "sp_d0": nom["sp_d0"].to_numpy(), "found": n_found, "dm_own": n_dm_own,
                                    "nmatch": nmatch[i_]}))
        fs["sec"] = round(time.time() - tf, 1)
        print(f"[{a.sample}] {fn}: D0 sel {fs['n_d0_sel']} okno {nd}, pary nom {fs['n_nom']} (znalezione RS {int((n_found == 1).sum())}), "
              f"kombinacje {len(Cf)}, córki dopasowane {fs['frac_daughters_matched']:.4f}, {fs['sec']}s", flush=True)
    if stat["n_event_mismatch"]:
        raise RuntimeError(f"eventNumber ROOT ≠ h5: {stat['n_event_mismatch']} (złe wyrównanie Tracks/CharmAnalysis/h5)")
    out = {}
    I8 = ("typ", "k", "flag", "hypo_swap", "qK", "found", "nmatch"); I16 = ("truth_ds", "nSelBL"); I32 = ("file", "row", "row_ds", "d0", "lb", "ntrk")
    I64 = ("entry", "run", "event"); U8 = ("trig",)

    def pack(prefix, parts, empty_cols):
        X = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame({c: [] for c in empty_cols})
        for c in X.columns:
            v = X[c].to_numpy()
            dt = np.int8 if c in I8 else np.int16 if c in I16 else np.int32 if c in I32 else np.int64 if c in I64 else np.uint8 if c in U8 \
                else bool if c == "pis_q_rs" else np.float32
            out[prefix + c] = np.nan_to_num(v, nan=-1).astype(dt) if dt in (np.int8, np.int16) else v.astype(dt)
        return len(X)
    nD = pack("d_", Dp, ["file", "row", "entry", "run", "event", "m", "pt", "eta", "score", "truth", "weight", "hypo_swap", "qK", "trig", "lb",
                         "nSelBL", "ntrk", "nmatch"])
    nC = pack("c_", Cp, ["d0", "typ", "k", "dm", "pis_pt", "pis_eta", "dr", "flag", "dm_nom", "truth_ds"])
    nN = pack("n_", Np, ["file", "row", "row_ds", "m", "pt", "score", "truth", "trig", "dm", "pis_pt", "pis_q_rs", "truth_ds", "sp_d0", "found",
                         "dm_own", "nmatch"])
    out["H_dm"] = H_dm; out["H_d0"] = H_d0
    stat.update({"sample": a.sample, "part": a.part, "nparts": a.nparts, "n_files_total": len(files_all), "n_d0": int(nD), "n_comb": int(nC),
                 "n_nom": int(nN), "numba": K.HAVE_NUMBA, "seconds": round(time.time() - t0, 1), "m0": m0, "m_win": M_WIN, "nom_win": NOM_WIN,
                 "lxy_min": a.lxy_min, "dm_max": a.dm_max, "row_dm": a.row_dm, "dm_max_nom": DM_MAX_NOM, "pt_bins": PT_BINS, "sc_n": SC_N,
                 "dm_edges": [float(DM_E[0]), float(DM_E[-1]), len(DM_E) - 1], "nrot": K.NROT, "nmix": K.NMIX, "match_dr": K.MATCH_DR,
                 "match_rel": K.MATCH_REL, "types": ["RS", "WS", "ROT", "MIX"], "thresholds": meta_m.get("thresholds"),
                 "trig_bits": {"J12": 1, "MBTS_ANY": 2, "mbts_L1MBTS_1": 4, "sptrk_pt4_L1MBTS_2": 8, "j20_L1MBTS_2": 16, "GRL": 32},
                 "flag_bits": {"nominal_pair": 1, "true_dstar_pion": 2}})
    out["meta"] = np.array(json.dumps(stat))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.sample}] zapisano {a.out}: D0 {nD}, kombinacji {nC}, par nominalnych {nN}, {os.path.getsize(a.out) / 1e6:.1f} MB, "
          f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    sys.exit(main())

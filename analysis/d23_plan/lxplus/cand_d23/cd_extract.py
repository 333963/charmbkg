#!/usr/bin/env python3
"""CAND_D23 (18.11.2026, projekt B) — kandydaci D⁺→Kππ per kandydat z h5 v2.0 produkcji dpt2to3 (D pT 2–3 GeV) + score H modeli
foldów D⁺ z dlz_data (H_s12345; TRENOWANE NA pT ≥ 3 GeV → tu POZA ZAKRESEM TRENINGU) + histogramy/liczniki D⁰ i D*.

  cd_extract.py --sample data_428770 --tag 0 --out-dir WD --dp-models .../dlz_data/train/models/H_s12345 plik1.h5 plik2.h5 ...

Wyjście w --out-dir (część <próbka>_<tag>):
  dp_<s>_<tag>.parquet   D⁺ (decayType 4, m ∈ species.Dplus.mass_range, pT ∈ [pt-lo, pt-hi), L_xy ≥ lxy-min), z preskalą (niżej)
  hist_<s>_<tag>.npz     histogramy PEŁNEJ statystyki (bez preskali): dp_h[klasa, bin pT, bin score, bin m] po wecie D*,
                         d0_h[klasa, bin pT, bin m] (D⁰, obie hipotezy), ds_h[RS/WS, klasa, bin pT, bin m(D0), bin Δm] (D*)
  part_<s>_<tag>.json    liczniki, zakres pT, pliki (zapisywany na końcu = znacznik kompletu)
Preskala wierszy D⁺ (deterministyczna, hasz zdarzenia `ps_u` ∈ [0, 10⁶); wiersz zostaje gdy ps_u < p·10⁶; waga w_ps = 1/p):
  ps_tier 3: MC, prawda sygnału (truth 0/1)                p = 1
  ps_tier 2: score ≥ --sc-full                             p = 1
  ps_tier 1: |m − m0| < --win i score < sc-full            p = --ps-win
  ps_tier 0: reszta                                        p = --ps-out
Hasz: dane (run, event); MC (10⁶ + nr pliku, wpis) — w fd_nd (run, event) nie identyfikuje zdarzenia.
Score: `score` = model foldu wg wzorca dlz_data (dane: fold_of(run, event); MC: fold_of(1, nr_pliku·10⁷ + blok) — formalnie jak w
treningu, ale numeracja plików d23 ≠ d150, więc to NIE jest fold treningowy; bez znaczenia, bo kandydaci pT < 3 nie byli w treningu),
`score_avg`/`score_sd` = średnia i rozrzut 5 foldów. Okno pT: charmbkg.samples.in_range z env CHARMBKG_PT_RANGE (kopia w src/).
"""
import argparse
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("CHARMBKG_SRC") or os.path.join(HERE, "src"))

TRIG = ["HLT_mb_mbts_L1MBTS_1", "HLT_noalg_L1J12", "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1",
        "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C", "HLT_j20_L1MBTS_2",
        "HLT_noalg_L1AFP_A_OR_C", "HLT_noalg_L1RD0_FILLED"]
TBP = ["L1_MBTS_1", "L1_MBTS_2", "L1_J12", "L1_AFP_A_OR_C", "L1_RD0_FILLED"]
EV = {"lb": ("lumiBlock", np.int32), "entry": ("event_entry", np.int64), "nSelBL": ("nSelBL", np.int16), "mu": ("average_mu", np.float32),
      "mu_act": ("actual_mu", np.float32), "npu": ("nPileupVertices", np.int16), "npv": ("nPrimaryVertices", np.int16),
      "ntrk": ("nTracksPt500Eta25", np.int16), "grl": ("passGRL_MinBiasAFP", np.int8)}
EV.update({"trig_" + t: ("trigPassed_" + t, np.int8) for t in TRIG})
EV.update({"tbp_" + t: ("trigTBP_" + t, np.int8) for t in TBP})
TRUTHX = {"D_tpdg": ("Dmeson_truthPdgId", np.int32), "n_truth_ok": ("Dmeson_n_daughter_truth_ok", np.int8)}
QCOL = ["Dmeson_daughter0_charge", "Dmeson_daughter2_charge"]
M_EDGES = np.linspace(1.70, 2.20, 501)                                   # D⁺: 1 MeV
S_EDGES = np.r_[np.arange(0, 90) * 0.01, 0.9 + np.arange(0, 51) * 0.002]  # score: 0,01 do 0,9; 0,002 do 1 (140 binów)
M0_EDGES = np.linspace(1.70, 2.05, 351)                                  # D⁰: 1 MeV (okno produkcji v2.0)
MD_EDGES = np.linspace(1.70, 2.05, 71)                                   # D*: m(D0) = m − Δm, 5 MeV
DM_EDGES = np.linspace(0.140, 0.200, 241)                                # D*: Δm, 0,25 MeV
NCLS = 4
CLS_DATA = ["all", "GRL & HLT_mb_mbts_L1MBTS_1", "GRL & HLT_noalg_L1J12", "GRL"]
CLS_MC = ["all", "truth 0/1 (sygnal; D0: bez hypo_swap=1)", "truth >= 8 (tlo)", "inne (truth 2-7; D0: takze hypo_swap=1)"]


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def ibin(x, edges):
    return np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)


def mix_u(a, b):
    from charmbkg.hyp import _mix64
    with np.errstate(over="ignore"):
        x = _mix64((np.asarray(a, np.int64).astype(np.uint64) * np.uint64(0xD6E8FEB86659FD93)) ^ np.asarray(b, np.int64).astype(np.uint64))
    return (x % np.uint64(1_000_000)).astype(np.int32)


def load_models(mdir):
    import xgboost as xgb
    meta = json.load(open(os.path.join(mdir, "meta.json"))) if os.path.exists(os.path.join(mdir, "meta.json")) else {}
    out = []
    for k in range(5):
        b = xgb.Booster()
        b.load_model(os.path.join(mdir, f"fold{k}.json"))
        bi = b.attr("best_iteration")
        bi = int(bi) if bi is not None else (int(meta["best_it"][k]) if meta.get("best_it") else b.num_boosted_rounds() - 1)
        out.append((b, bi))
    return out, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--dp-models", required=True, help="katalog z fold{0..4}.json (+ meta.json) modelu H D⁺ (dlz_data/train/models/H_s12345)")
    ap.add_argument("--pt-lo", type=float, default=2.0)
    ap.add_argument("--pt-hi", type=float, default=3.0)
    ap.add_argument("--lxy-min", type=float, default=0.25)
    ap.add_argument("--sc-full", type=float, default=0.93)
    ap.add_argument("--win", type=float, default=0.06)
    ap.add_argument("--ps-win", type=float, default=0.05)
    ap.add_argument("--ps-out", type=float, default=0.01)
    ap.add_argument("--chunk", type=int, default=500_000)
    ap.add_argument("--max-rows", type=int, default=0, help="test: tylko pierwsze N wierszy każdego pliku")
    ap.add_argument("--nthread", type=int, default=2)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    os.environ["CHARMBKG_PT_RANGE"] = f"{a.pt_lo},{a.pt_hi}"      # obejście twardego pT ≥ 3 (src/charmbkg/samples.py: pt_range)
    import h5py
    import pyarrow as pa
    import pyarrow.parquet as pq
    import xgboost as xgb
    from charmbkg import config as C, features as F
    from charmbkg.hyp import event_blocks
    from charmbkg.samples import colmap, dstar_veto_mask, fold_of, in_range, pt_range
    assert pt_range() == (a.pt_lo, a.pt_hi), pt_range()
    cfg = C.load()
    SP = {**cfg["species"]["Dplus"], "n_daughters": 3}
    S0 = cfg["species"]["D0"]
    feats = F.feature_list(cfg, "H")
    models, mmeta = load_models(a.dp_models)
    if mmeta.get("features") and list(mmeta["features"]) != feats:
        raise SystemExit(f"cechy modelu ≠ cechy H paczki: {mmeta['features']} vs {feats}")
    for b, _ in models:
        b.set_param({"nthread": a.nthread})
        if b.num_features() != len(feats):
            raise SystemExit(f"model: {b.num_features()} cech, oczekiwane {len(feats)}")
    is_mc = not a.sample.startswith("data")
    PT_EDGES = np.r_[np.arange(a.pt_lo, a.pt_hi - 1e-9, 0.5), a.pt_hi]
    npt = len(PT_EDGES) - 1
    H_dp = np.zeros((NCLS, npt, len(S_EDGES) - 1, len(M_EDGES) - 1), np.int64)
    H_d0 = np.zeros((NCLS, npt, len(M0_EDGES) - 1), np.int64)
    H_ds = np.zeros((2, NCLS, npt, len(MD_EDGES) - 1, len(DM_EDGES) - 1), np.int64)
    cnt = {k: 0 for k in ("rows", "dp_dt", "dp_pt_outside", "dp_sel", "dp_veto", "dp_kept", "dp_score_nan", "dp_mc_sig",
                          "d0_dt", "d0_pt_outside", "d0_sel", "d0_mc_sig", "d0_mc_swap", "ds_dt", "ds_pt_outside", "ds_sel", "ds_rs",
                          "ds_rs_dm_signal", "ds_ws_dm_signal", "ds_mc_sig")}
    cnt.update({f"dp_tier{i}": 0 for i in range(4)})
    ptmm = [np.inf, -np.inf]
    missing_ev = set()
    os.makedirs(a.out_dir, exist_ok=True)
    base = f"{a.sample}_{a.tag}"
    out_pq = os.path.join(a.out_dir, f"dp_{base}.parquet")
    W = None
    fidx_list = []

    def classes(ev, truth, swap=None):
        c = np.zeros((NCLS, len(truth)), bool)
        c[0] = True
        if is_mc:
            sig = np.isin(truth, cfg["signal"]["truth_values"])
            if swap is not None:
                sig &= swap != 1
            c[1] = sig
            c[2] = truth >= 8
            c[3] = ~c[1] & ~c[2]
        else:
            g = ev["grl"] == 1
            c[1] = g & (ev["trig_HLT_mb_mbts_L1MBTS_1"] == 1)
            c[2] = g & (ev["trig_HLT_noalg_L1J12"] == 1)
            c[3] = g
        return c

    for path in a.files:
        mfi = re.search(r"_(\d{6})_", os.path.basename(path))
        fidx = int(mfi.group(1)) if mfi else len(fidx_list)
        fidx_list.append(fidx)
        cm = colmap(cfg, path, 3)
        with h5py.File(path, "r") as f:
            ds = f[cfg["io"]["dataset"]]
            names = set(ds.dtype.names)
            evmap = {k: v[0] for k, v in {**EV, **TRUTHX}.items() if v[0] in names}
            missing_ev |= {v[0] for k, v in {**EV, **TRUTHX}.items() if v[0] not in names}
            hs = cm.get("hypo_swap")
            cols = list(dict.fromkeys(F.raw_columns(cm, 3) + list(evmap.values()) + [c for c in QCOL if c in names]))
            cols = [c for c in cols if c in names]
            n = ds.shape[0] if not a.max_rows else min(ds.shape[0], a.max_rows)
            prev = None
            for start in range(0, n, a.chunk):
                stop = min(n, start + a.chunk)
                arr = ds.fields(cols)[start:stop]
                cnt["rows"] += len(arr)
                g = lambda c: np.asarray(arr[c])  # noqa: E731
                run = g(cm["run"]).astype(np.int64) if cm.has("run") else np.zeros(len(arr), np.int64)
                evn = g(cm["event"]).astype(np.int64)
                if "entry" in evmap:
                    evb, prev = event_blocks(np.zeros(len(arr), np.int64), g(evmap["entry"]).astype(np.int64), 0, prev)
                    evid = g(evmap["entry"]).astype(np.int64)
                else:
                    evb, prev = event_blocks(run, evn, 0, prev)
                    evid = evb
                u = mix_u(1_000_000 + fidx, evid) if is_mc else mix_u(run, evn)
                dt, pt, m = g(cm["decay_type"]), g(cm["pt"]).astype(np.float64), g(cm["m"]).astype(np.float64)
                truth = g(cm["truth"]).astype(np.int32) if cm.has("truth") else np.zeros(len(arr), np.int32)
                inpt = (pt >= a.pt_lo) & (pt < a.pt_hi)
                evc = {k: (g(evmap[k]) if k in evmap else np.full(len(arr), -1)) for k in ("grl", "trig_HLT_mb_mbts_L1MBTS_1", "trig_HLT_noalg_L1J12")}
                # ---------------- D⁰ (liczniki + histogram masy, obie hipotezy) i D* (liczniki + histogram m(D0) × Δm)
                s0 = dt == int(S0["decay_type"])
                cnt["d0_dt"] += int(s0.sum()); cnt["d0_pt_outside"] += int((s0 & ~inpt).sum())
                s0 &= inpt & (m >= M0_EDGES[0]) & (m < M0_EDGES[-1])
                if s0.any():
                    i0 = np.flatnonzero(s0)
                    sw = g(hs)[i0] if hs and hs in names else None
                    cl = classes({k: v[i0] for k, v in evc.items()}, truth[i0], sw)
                    flat = ibin(pt[i0], PT_EDGES) * (len(M0_EDGES) - 1) + ibin(m[i0], M0_EDGES)
                    for c in range(NCLS):
                        H_d0[c] += np.bincount(flat[cl[c]], minlength=H_d0[c].size).reshape(H_d0[c].shape)
                    cnt["d0_sel"] += len(i0)
                    if is_mc:
                        cnt["d0_mc_sig"] += int(cl[1].sum())
                        cnt["d0_mc_swap"] += int((np.isin(truth[i0], cfg["signal"]["truth_values"]) & (sw == 1)).sum()) if sw is not None else 0
                s1 = dt == 1
                cnt["ds_dt"] += int(s1.sum()); cnt["ds_pt_outside"] += int((s1 & ~inpt).sum())
                if s1.any() and cm.has("dm"):
                    dm = g(cm["dm"]).astype(np.float64)
                    md = m - dm
                    s1 &= inpt & (dm >= DM_EDGES[0]) & (dm < DM_EDGES[-1]) & (md >= MD_EDGES[0]) & (md < MD_EDGES[-1])
                    i1 = np.flatnonzero(s1)
                    if len(i1):
                        if all(c in names for c in QCOL):
                            rs = g(QCOL[0])[i1] == -g(QCOL[1])[i1]
                        else:
                            rs = np.ones(len(i1), bool)
                        cl = classes({k: v[i1] for k, v in evc.items()}, truth[i1])
                        shp = H_ds[0, 0].shape
                        flat = (ibin(pt[i1], PT_EDGES) * shp[1] + ibin(md[i1], MD_EDGES)) * shp[2] + ibin(dm[i1], DM_EDGES)
                        for c in range(NCLS):
                            for r, sel in ((0, rs), (1, ~rs)):
                                H_ds[r, c] += np.bincount(flat[cl[c] & sel], minlength=H_ds[r, c].size).reshape(shp)
                        dsig = (dm[i1] >= 0.1435) & (dm[i1] < 0.1475) & (np.abs(md[i1] - S0["m0"]) < 0.04)
                        cnt["ds_sel"] += len(i1); cnt["ds_rs"] += int(rs.sum())
                        cnt["ds_rs_dm_signal"] += int((rs & dsig).sum()); cnt["ds_ws_dm_signal"] += int((~rs & dsig).sum())
                        if is_mc:
                            cnt["ds_mc_sig"] += int(cl[1].sum())
                # ---------------- D⁺ per kandydat
                sp = dt == int(SP["decay_type"])
                cnt["dp_dt"] += int(sp.sum()); cnt["dp_pt_outside"] += int((sp & ~inpt).sum())
                idx = np.flatnonzero(sp & inpt & (m >= SP["mass_range"][0]) & (m < SP["mass_range"][1]))
                if not len(idx):
                    continue
                raw = pd.DataFrame({c: np.asarray(arr[c])[idx] for c in cols})
                L = F.to_logical(raw, cm, 3)
                if not cm.has("run"):
                    L["run"] = 0
                L["weight"] = L["weight"].fillna(1.0)
                L = F.compute(L, SP, cfg.get("units"))
                L["dstar_ok"] = dstar_veto_mask(L, SP)
                keep = in_range(L, SP).to_numpy() & (L["Lxy"].to_numpy() >= a.lxy_min)
                L, raw, idx = L[keep].reset_index(drop=True), raw[keep].reset_index(drop=True), idx[keep]
                if not len(L):
                    continue
                nn = len(L)
                mm, pp = L["m"].to_numpy(np.float64), L["pt"].to_numpy(np.float64)
                if pp.min() < a.pt_lo or pp.max() >= a.pt_hi:
                    raise SystemExit(f"BŁĄD: pT poza [{a.pt_lo}, {a.pt_hi}): {pp.min()} … {pp.max()}")
                ptmm = [min(ptmm[0], float(pp.min())), max(ptmm[1], float(pp.max()))]
                X = L[feats].to_numpy(np.float32)
                dmat = xgb.DMatrix(X)
                P = np.stack([b.predict(dmat, iteration_range=(0, bi + 1)) for b, bi in models])
                if is_mc:
                    fold = fold_of(np.full(nn, 1), fidx * 10_000_000 + evb[idx], 5)
                else:
                    fold = fold_of(run[idx], evn[idx], 5)
                score = P[fold.astype(np.int64), np.arange(nn)].astype(np.float32)
                cnt["dp_score_nan"] += int(np.isnan(score).sum())
                ev = {k: (raw[v[0]].to_numpy() if k in evmap else np.full(nn, -1)) for k, v in {**EV, **TRUTHX}.items()}
                tr = truth[idx]
                sig = np.isin(tr, cfg["signal"]["truth_values"]) if is_mc else np.zeros(nn, bool)
                ok = L["dstar_ok"].to_numpy(bool)
                cl = classes(ev, tr)
                shp = H_dp[0].shape
                flat = (ibin(pp, PT_EDGES) * shp[1] + ibin(score, S_EDGES)) * shp[2] + ibin(mm, M_EDGES)
                for c in range(NCLS):
                    H_dp[c] += np.bincount(flat[cl[c] & ok], minlength=H_dp[c].size).reshape(shp)
                tier = np.zeros(nn, np.int8)
                tier[np.abs(mm - SP["m0"]) < a.win] = 1
                tier[score >= a.sc_full] = 2
                tier[sig] = 3
                p = np.array([a.ps_out, a.ps_win, 1.0, 1.0])[tier]
                uu = u[idx]
                kp = uu < np.round(p * 1_000_000)
                cnt["dp_sel"] += nn; cnt["dp_veto"] += int((~ok).sum()); cnt["dp_kept"] += int(kp.sum()); cnt["dp_mc_sig"] += int(sig.sum())
                for i in range(4):
                    cnt[f"dp_tier{i}"] += int((kp & (tier == i)).sum())
                if not kp.any():
                    continue
                O = {"run": run[idx].astype(np.int32), "lb": ev["lb"].astype(np.int32), "event": evn[idx], "entry": ev["entry"].astype(np.int64),
                     "file_idx": np.full(nn, fidx, np.int32), "row": (idx + start).astype(np.int32), "evblk": evb[idx].astype(np.int32),
                     "m": mm, "pt": pp.astype(np.float32), "eta": L["eta"].to_numpy(np.float32), "phi": L["phi"].to_numpy(np.float32),
                     "charge": np.nan_to_num(L["charge"].to_numpy(float)).astype(np.int8),
                     "m_kpi1": L["m_kpi1"].to_numpy(np.float32), "m_kpi2": L["m_kpi2"].to_numpy(np.float32), "dstar_veto": (~ok).astype(np.int8),
                     "score": score, "score_avg": P.mean(0).astype(np.float32), "score_sd": P.std(0).astype(np.float32), "fold": fold.astype(np.int8),
                     "truth": tr.astype(np.int16), "is_sig": sig.astype(np.int8) if is_mc else np.full(nn, -1, np.int8),
                     "weight": L["weight"].to_numpy(np.float32)}
                for k, v in TRUTHX.items():
                    O[k] = ev[k].astype(v[1])
                for k, v in EV.items():
                    if k not in O:
                        O[k] = np.nan_to_num(ev[k].astype(np.float64), nan=-1).astype(v[1])
                for c in feats:
                    O[c] = X[:, feats.index(c)]
                O["ps_tier"] = tier
                O["ps_u"] = uu
                O["w_ps"] = (1.0 / p).astype(np.float32)
                t = pa.table({k: v[kp] for k, v in O.items()})
                if W is None:
                    codec = next((c for c in ("zstd", "gzip", "snappy") if pa.Codec.is_available(c)), "none")
                    W = pq.ParquetWriter(out_pq + ".tmp", t.schema, compression=codec)
                elif t.schema != W.schema:
                    t = t.cast(W.schema)
                W.write_table(t)
        log(f"{os.path.basename(path)}: D+ w oknie {cnt['dp_sel']}, zapisane {cnt['dp_kept']}; D0 {cnt['d0_sel']}; D* {cnt['ds_sel']}")
    if W is not None:
        W.close()
        os.replace(out_pq + ".tmp", out_pq)
    np.savez_compressed(os.path.join(a.out_dir, f"hist_{base}.npz"), dp_h=H_dp, d0_h=H_d0, ds_h=H_ds, m_edges=M_EDGES, s_edges=S_EDGES,
                        pt_edges=PT_EDGES, m0_edges=M0_EDGES, md_edges=MD_EDGES, dm_edges=DM_EDGES,
                        cls=np.array(CLS_MC if is_mc else CLS_DATA), ds_sign=np.array(["RS q(K)=-q(pi_s)", "WS"]))
    meta = {"sample": a.sample, "tag": a.tag, "is_mc": is_mc, "files": [os.path.basename(q) for q in a.files], "file_idx": fidx_list,
            "pt_range": [a.pt_lo, a.pt_hi], "dp_pt_min": None if not np.isfinite(ptmm[0]) else ptmm[0],
            "dp_pt_max": None if not np.isfinite(ptmm[1]) else ptmm[1], "counts": cnt,
            "dp_bytes": os.path.getsize(out_pq) if os.path.exists(out_pq) else 0, "missing_columns": sorted(missing_ev),
            "settings": {"lxy_min": a.lxy_min, "sc_full": a.sc_full, "win": a.win, "ps_win": a.ps_win, "ps_out": a.ps_out, "max_rows": a.max_rows},
            "model": {"dir": a.dp_models, "best_it": [bi for _, bi in models], "features": feats,
                      "note": "H_s12345 z dlz_data trenowany na pT >= 3 GeV; score w binie 2-3 GeV = ekstrapolacja poza zakres treningu"},
            "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
    with open(os.path.join(a.out_dir, f"part_{base}.json"), "w") as f:
        json.dump(meta, f, indent=1)
    log(f"OK {base}: {json.dumps(cnt)} pT D+ [{ptmm[0]:.4f}, {ptmm[1]:.4f}] brak kolumn: {sorted(missing_ev)}")
    if cnt["dp_score_nan"]:
        raise SystemExit("BŁĄD: score NaN")


if __name__ == "__main__":
    main()

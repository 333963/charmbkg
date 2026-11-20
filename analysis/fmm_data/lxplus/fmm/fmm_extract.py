#!/usr/bin/env python3
"""FMM — wyciąg par bliźniaczych D0 (K↔π) z h5 v2.0 (compound `events`) do pomiaru/ograniczenia udziału MisMatch* pod pikiem
(projekt B, otwarty punkt METHODS_NOTE §7 / B1). Logika wczytywania i score H = kopia `B_trig/btrig_extract.py` (ten sam model
nominalny `models/D0_H`, charmtag 2.0.2 z `ev_dstar/src`).

Wiersze: decayType 10, L_xy ≥ --lxy-min (0,25 mm = nominał v1.2), 3 ≤ pT < 30 GeV, m ∈ [--m-lo, --m-hi] (domyślnie pełne okno v2.0
1,70–2,05 GeV), podpróbka całych zdarzeń --frac (hash (run, event) w danych, event_uid w MC → bliźniaki zawsze razem).
Kolumny (npz):
  klucz    run, event, uid (event_uid), file_idx, row
  kandydat m, m_swap (masa TEJ SAMEJ pary torów z zamienionym przypisaniem K↔π, analitycznie — istnieje zawsze), pt, eta, lxy,
           score (H, OOF po foldach jak w apply), fold, k_idx (która córka jest K w tej hipotezie), q0, q1
  bliźniak n_grp (ile wierszy D0 dzieli tę parę torów po cięciach L_xy/pT/okno), tw_m, tw_score (masa/score wiersza-bliźniaka; NaN gdy
           bliźniak nie jest zapisany, tj. m_swap poza oknem ntupla), tw_row
  cechy H  f_<cecha> (float16; wszystkie cechy modelu)
  PID/tory dedx0/1, npts0/1, dpt0/1, deta0/1, dphi0/1, dd00/1, dz00/1
  prawda   truth (Dmeson_truth_category; 0/1 sygnał, 2–7 częściowo dopasowany charm, 8–11 MisMatch*, 12 OtherNoB, 13–15 Other*B),
           hypo_swap (TYLKO MC: 1 = hipoteza K↔π odwrotna do truthPdgId córek, 0 = zgodna, −1 = para nie jest K+π / dane),
           D_tpdg, tpdg0/1, tpar0/1 (pdgId prawdy córek i ich rodziców), n_tok
  zdarzenie tp_<chain> (int8; MC = −1), grl, nSelBL, lead_pt [GeV], npu, mu, w
Uwaga: grupowanie bliźniaków w porcjach --chunk; para rozcięta granicą porcji gubi tw_* (≲ 1e-5) — obserwablą nominalną jest m_swap.
Użycie: fmm_extract.py --models DIR --out OUT.npz --tag data|fd|mb [--frac 0.3] [--truth-max 7] FILE.h5 [...]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import yaml

TRIG = ["HLT_j20_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C",
        "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_noalg_L1AFP_A_OR_C", "HLT_noalg_L1J12", "HLT_noalg_L1RD0_FILLED", "HLT_noalg_mb_L1MBTS_1",
        "HLT_noalg_mb_L1MBTS_2"]
EV = {"grl": "passGRL_MinBiasAFP", "nSelBL": "nSelBL", "lead_pt": "leadTrackPtBL", "npu": "nPileupVertices", "mu": "actual_mu",
      "uid": "event_uid"}
RAWX = {"k_pdg0": "Dmeson_daughter0_pdgId", "q0": "Dmeson_daughter0_charge", "q1": "Dmeson_daughter1_charge",
        "trk0": "Dmeson_daughter0_trackId", "trk1": "Dmeson_daughter1_trackId",
        "tpdg0": "Dmeson_daughter0_truthPdgId", "tpdg1": "Dmeson_daughter1_truthPdgId",
        "tpar0": "Dmeson_daughter0_truthParentPdgId", "tpar1": "Dmeson_daughter1_truthParentPdgId",
        "D_tpdg": "Dmeson_truthPdgId", "hypo_swap": "Dmeson_hypo_swap", "n_tok": "Dmeson_n_daughter_truth_ok",
        "truth_raw": "Dmeson_truth_category",
        "dedx0": "Dmeson_daughter0_dEdx", "dedx1": "Dmeson_daughter1_dEdx", "npts0": "Dmeson_daughter0_npts_dEdx",
        "npts1": "Dmeson_daughter1_npts_dEdx", "dpt0": "Dmeson_daughter0_pt", "dpt1": "Dmeson_daughter1_pt",
        "deta0": "Dmeson_daughter0_eta", "deta1": "Dmeson_daughter1_eta", "dphi0": "Dmeson_daughter0_phi", "dphi1": "Dmeson_daughter1_phi",
        "dd00": "Dmeson_daughter0_d0", "dd01": "Dmeson_daughter1_d0", "dz00": "Dmeson_daughter0_z0SinTheta",
        "dz01": "Dmeson_daughter1_z0SinTheta"}
DTYPE = {"run": np.int32, "event": np.int64, "uid": np.uint64, "file_idx": np.int16, "row": np.int32, "tw_row": np.int32, "fold": np.int8,
         "k_idx": np.int8, "q0": np.int8, "q1": np.int8, "n_grp": np.int8, "truth": np.int16, "hypo_swap": np.int8, "D_tpdg": np.int32,
         "tpdg0": np.int32, "tpdg1": np.int32, "tpar0": np.int32, "tpar1": np.int32, "n_tok": np.int8, "npts0": np.int8, "npts1": np.int8,
         "grl": np.int8, "nSelBL": np.int16, "npu": np.int8, "m": np.float32, "m_swap": np.float32, "tw_m": np.float32,
         "tw_score": np.float32, "score": np.float32, "pt": np.float32, "eta": np.float32, "lxy": np.float32, "w": np.float32,
         "dpt0": np.float32, "dpt1": np.float32, "deta0": np.float32, "deta1": np.float32, "dphi0": np.float32, "dphi1": np.float32,
         "lead_pt": np.float32, "mu": np.float32}
MK, MPI = 0.493677, 0.13957039


def keep_hash(a, b, frac):
    """Deterministyczna podpróbka zdarzeń (splitmix64) — ta sama decyzja w każdym pliku/porcji (jak B_trig)."""
    if frac >= 1.0:
        return np.ones(len(a), bool)
    with np.errstate(over="ignore"):
        x = (np.asarray(a).astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15)) ^ np.asarray(b).astype(np.uint64)
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        x = x ^ (x >> np.uint64(31))
    return (x >> np.uint64(11)).astype(np.float64) / float(1 << 53) < frac


def swap_mass(pt0, eta0, phi0, pt1, eta1, phi1, k_is0):
    """m(Kπ) pary torów z przypisaniem ODWROTNYM do hipotezy wiersza (k_is0: córka 0 jest K w hipotezie wiersza). GeV."""
    def p4(pt, eta, phi, m):
        px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
        return np.sqrt(px * px + py * py + pz * pz + m * m), px, py, pz
    a = p4(pt0, eta0, phi0, np.where(k_is0, MPI, MK))
    b = p4(pt1, eta1, phi1, np.where(k_is0, MK, MPI))
    return np.sqrt(np.maximum((a[0] + b[0]) ** 2 - sum((a[i] + b[i]) ** 2 for i in (1, 2, 3)), 0.0))


def _chunks(path, cols, chunk, dt_col, dt, max_rows):
    import h5py
    with h5py.File(path, "r") as f:
        ds = f["events"]
        names = set(ds.dtype.names)
        cols = [c for c in cols if c in names]
        n = ds.shape[0] if not max_rows else min(ds.shape[0], max_rows)
        for a in range(0, n, chunk):
            b = min(n, a + chunk)
            x = ds.fields(cols)[a:b]
            df = pd.DataFrame({c: x[c] for c in cols})
            df["_row"] = np.arange(a, b)
            yield df[df[dt_col].to_numpy() == dt].reset_index(drop=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", required=True, help="data|fd|mb")
    ap.add_argument("--chunk", type=int, default=500_000)
    ap.add_argument("--m-lo", type=float, default=1.70)
    ap.add_argument("--m-hi", type=float, default=2.05)
    ap.add_argument("--lxy-min", type=float, default=0.25)
    ap.add_argument("--frac", type=float, default=1.0, help="podpróbka całych zdarzeń")
    ap.add_argument("--truth-max", type=int, default=-1, help="≥0: zachowaj tylko truth ≤ N (fd: 7 → sygnał + częściowo dopasowany charm)")
    ap.add_argument("--max-rows", type=int, default=0, help="test: tylko pierwsze N wierszy pliku")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args(argv)

    from charmtag import features as F
    from charmtag.io import resolve_columns
    from charmtag.model import Ensemble

    cfg = yaml.safe_load(open(os.path.join(a.models, "config_used.yaml")))
    meta = json.load(open(os.path.join(a.models, "meta.json")))
    sp = cfg["species"]["D0"]
    units = cfg.get("units") or {}
    upt = float(units.get("d_pt", 1.0))
    ens = Ensemble(a.models)
    feats = list(ens.feats)
    parts, nfile, t0 = [], {}, time.time()
    for fi, path in enumerate(a.files):
        import h5py
        with h5py.File(path, "r") as f:
            names = list(f["events"].dtype.names)
        cm = resolve_columns(names, cfg["columns"], 2)
        dt_col = cm.get("decay_type")
        miss = [c for c in list(RAWX.values()) + list(EV.values()) if c not in names]
        if miss:
            sys.exit(f"BŁĄD: {path}: brak kolumn {miss} (to nie jest h5 v2.0?)")
        trig = {c: f"trigPassed_{c}" for c in TRIG if f"trigPassed_{c}" in names}
        need = sorted({c for c in F.raw_columns(cm, 2) + [dt_col] + list(RAWX.values()) + list(EV.values()) + list(trig.values()) if c})
        n_in = n_out = n_tw = 0
        for r0 in _chunks(path, need, a.chunk, dt_col, int(sp["decay_type"]), a.max_rows):
            if len(r0) == 0:
                continue
            L = F.to_logical(r0, cm, 2)
            L["weight"] = L["weight"].fillna(1.0)
            n_in += len(L)
            m, pt = L["m"].to_numpy(float), L["pt"].to_numpy(float)
            uid = r0[EV["uid"]].to_numpy().astype(np.uint64)
            ok = (L["Lxy"].to_numpy(float) >= a.lxy_min) & (pt >= 3) & (pt < 30) & (m >= a.m_lo) & (m <= a.m_hi)
            if a.tag == "data":
                ok &= keep_hash(L["run"].to_numpy(), L["event"].to_numpy(), a.frac)
            else:
                ok &= keep_hash(uid, uid >> np.uint64(22), a.frac)
            if not ok.any():
                continue
            L, r0 = L[ok].reset_index(drop=True), r0[ok].reset_index(drop=True)
            L = F.compute(L, sp, units)
            score, fold = ens.score(L)
            O = pd.DataFrame({"run": L["run"].to_numpy(), "event": L["event"].to_numpy(), "file_idx": fi, "row": r0["_row"].to_numpy(),
                              "m": L["m"].to_numpy(float), "pt": L["pt"].to_numpy(float), "eta": L["eta"].to_numpy(float),
                              "lxy": L["Lxy"].to_numpy(float), "score": np.asarray(score, float), "fold": np.asarray(fold),
                              "w": L["weight"].to_numpy(float), "truth": r0[RAWX["truth_raw"]].to_numpy()})
            for k, c in RAWX.items():
                if k not in ("truth_raw", "k_pdg0", "trk0", "trk1"):
                    O[k] = r0[c].to_numpy()
            O["dpt0"] = O["dpt0"].to_numpy(float) * upt
            O["dpt1"] = O["dpt1"].to_numpy(float) * upt
            k0 = np.abs(r0[RAWX["k_pdg0"]].to_numpy()) == 321
            O["k_idx"] = np.where(k0, 0, 1)
            O["m_swap"] = swap_mass(O["dpt0"].to_numpy(float), O["deta0"].to_numpy(float), O["dphi0"].to_numpy(float),
                                    O["dpt1"].to_numpy(float), O["deta1"].to_numpy(float), O["dphi1"].to_numpy(float), k0)
            for k, c in EV.items():
                O[k] = r0[c].to_numpy()
            for c, col in trig.items():
                O[f"tp_{c}"] = r0[col].to_numpy().astype(np.int8)
            for c in feats:
                O[f"f_{c}"] = L[c].to_numpy(float).astype(np.float16)
            # bliźniaki: ta sama (nieuporządkowana) para torów w tym samym zdarzeniu
            ta, tb = r0[RAWX["trk0"]].to_numpy().astype(np.int64), r0[RAWX["trk1"]].to_numpy().astype(np.int64)
            key = pd.DataFrame({"u": uid[ok], "lo": np.minimum(ta, tb), "hi": np.maximum(ta, tb)})
            g = key.groupby(["u", "lo", "hi"], sort=False)
            gid, n = g.ngroup().to_numpy(), g["u"].transform("size").to_numpy()
            tw = np.full(len(O), -1, np.int64)
            idx = np.flatnonzero(n == 2)
            order = idx[np.argsort(gid[idx], kind="stable")]
            tw[order[0::2]], tw[order[1::2]] = order[1::2], order[0::2]
            has = tw >= 0
            twc = np.clip(tw, 0, None)
            O["n_grp"] = np.minimum(n, 127)
            O["tw_m"] = np.where(has, O["m"].to_numpy()[twc], np.nan)
            O["tw_score"] = np.where(has, O["score"].to_numpy()[twc], np.nan)
            O["tw_row"] = np.where(has, O["row"].to_numpy()[twc], -1)
            if a.truth_max >= 0:
                keep = (O["truth"] <= a.truth_max).to_numpy()
                O, has = O[keep], has[keep]
            n_out += len(O)
            n_tw += int(has.sum())
            parts.append(O)
        nfile[os.path.basename(path)] = [n_in, n_out]
        print(f"[{a.tag}] {os.path.basename(path)}: D0 {n_in} → {n_out} (z bliźniakiem w oknie {n_tw})  {time.time() - t0:.0f}s", flush=True)

    P = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["run", "event", "m", "n_grp"])
    out = {}
    for c in P.columns:
        v = P[c].to_numpy()
        if c in DTYPE:
            out[c] = v.astype(DTYPE[c])
        elif c.startswith("tp_"):
            out[c] = v.astype(np.int8)
        else:
            out[c] = v.astype(np.float16)
    out["meta"] = np.array(json.dumps({"tag": a.tag, "files": [os.path.basename(p) for p in a.files], "rows_per_file": nfile,
                                       "features": feats, "thresholds": meta.get("thresholds"), "m0": float(sp["m0"]),
                                       "window": [a.m_lo, a.m_hi], "lxy_min": a.lxy_min, "frac": a.frac, "truth_max": a.truth_max,
                                       "trig": TRIG, "n": int(len(P)), "fmm": "1.0"}))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.tag}] zapisano {a.out}: {len(P)} wierszy, {len(out) - 1} kolumn, {os.path.getsize(a.out) / 1e6:.1f} MB, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""B_trig — ekstrakcja z h5 (v1.2 albo v2.0, compound `events`) do analizy odporności na trigger. Kopia logiki
`Q_lxplus/ev_dstar/ev_extract.py` (pary D0←D*, score H modelem nominalnym, progi WP) + dodatki B_trig:

  1. pary: + ptcone40/ptcone100 wiersza D0 (izolacja, gdy kolumny są w h5) i kolumny zdarzeniowe e_<nazwa> (e_tp_<chain>,
     e_nSelBL, e_lead_pt, e_actual_mu, e_grl, e_npu) wzięte z wiersza D0 — jeśli są w h5 (MC: jedyny bezpieczny klucz);
  2. `all_*` (tło do ε_B): wszystkie D0 w oknie m ∈ [--all-lo, --all-hi] (domyślnie 1,75–2,00 GeV = okno v1.2; v2.0 ma 1,70–2,05,
     tniemy do 1,75–2,00, żeby SB były jak w nominale), Z KLUCZEM zdarzenia (all_run, all_event, all_file) — EV go nie zapisywał;
     podpróbka zdarzeń --all-frac (hash (run, event), deterministyczna);
  3. tabela zdarzeń `ev_*` (jeden wiersz na zdarzenie z D0/D*): run, event, file_idx, tp_/ps_<chain>, nSelBL, lead_pt, …
     (pusta poza run/event, gdy h5 nie ma tych kolumn → wtedy klasy z `trg_root_extract.py` + `--event-table`).

Nazwy kolumn trigger: tp_/ps_, trigPassed_/trigPS_, EventInfo_trigPassed_* (btrig_common.find_trig_cols); lead_pt: jeżeli mediana
> 200 → MeV (ROOT) → /1000.
Wyjście: npz jak EV (ev_analysis.py działa bez zmian) + nowe klucze. Użycie:
  btrig_extract.py --models DIR --out OUT.npz --tag data|fd|mb [--all-frac 0.3] [--chunk 500000] FILE.h5 [...]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from btrig_common import EVENT_EXTRA, TRIG_CHAINS, find_col, find_trig_cols  # noqa: E402

DM_MAX = 0.165
M_WIN = 0.075
ISO_COLS = {"ptcone40": ["Dmeson_ptcone40", "ptcone40"], "ptcone100": ["Dmeson_ptcone100", "ptcone100"]}


def keep_hash(run, event, frac):
    """Deterministyczna podpróbka zdarzeń (splitmix64 na (run, event)) — ta sama decyzja w każdym pliku/porcji."""
    if frac >= 1.0:
        return np.ones(len(run), bool)
    x = (np.asarray(run, np.uint64) * np.uint64(0x9E3779B97F4A7C15)) ^ np.asarray(event, np.uint64)
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    x = x ^ (x >> np.uint64(31))
    return (x >> np.uint64(11)).astype(np.float64) / float(1 << 53) < frac


def event_cols(names):
    """{nazwa logiczna: kolumna h5} dla kolumn zdarzeniowych (trigger + EVENT_EXTRA) obecnych w pliku."""
    m = {}
    for k, c in EVENT_EXTRA.items():
        b = find_col(names, c)
        if b:
            m[k] = b
    for c, v in find_trig_cols(names, TRIG_CHAINS).items():
        if v["pass"]:
            m[f"tp_{c}"] = v["pass"]
        if v["ps"]:
            m[f"ps_{c}"] = v["ps"]
    return m


def _load_h5_chunks(path, cols, chunk, dt_col, keep_dt):
    import h5py
    with h5py.File(path, "r") as f:
        ds = f["events"]
        names = set(ds.dtype.names)
        cols = [c for c in cols if c in names]
        n = ds.shape[0]
        for a in range(0, n, chunk):
            b = min(n, a + chunk)
            x = ds.fields(cols)[a:b]
            df = pd.DataFrame({c: x[c] for c in cols})
            df["_local_row"] = np.arange(a, b)
            df = df[np.isin(df[dt_col].to_numpy(), keep_dt)]
            yield df


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", required=True, help="katalog modelu D0_H (config_used.yaml, meta.json, fold*.json)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", required=True, help="data|fd|mb")
    ap.add_argument("--chunk", type=int, default=500_000)
    ap.add_argument("--all-lo", type=float, default=1.75)
    ap.add_argument("--all-hi", type=float, default=2.00)
    ap.add_argument("--all-frac", type=float, default=1.0, help="podpróbka zdarzeń dla all_* (dane: 0,3 wystarcza do ε_B)")
    ap.add_argument("--lxy-min", type=float, default=0.25,
                    help="L_xy(D0) i L_xy(D*) ≥ [mm]; 0,25 = v1.1/v1.2 (h5_v20 ma 0,10 — bez cięcia v2.0 ≠ v1.2; README h5_v20)")
    ap.add_argument("--max-rows", type=int, default=0, help="test: tylko pierwsze N wierszy pliku")
    ap.add_argument("files", nargs="+")
    args = ap.parse_args(argv)

    from charmtag import features as F
    from charmtag.io import resolve_columns
    from charmtag.kinematics import swapped_mass_sp
    from charmtag.model import Ensemble
    from charmtag.samples import mark_swaps

    cfg = yaml.safe_load(open(os.path.join(args.models, "config_used.yaml")))
    meta = json.load(open(os.path.join(args.models, "meta.json")))
    sp = cfg["species"]["D0"]
    dc = cfg["dstar_control"]
    m0 = float(sp["m0"])
    units = cfg.get("units")
    ens = Ensemble(args.models)
    feats = list(ens.feats)
    out_parts, all_parts, ev_parts = [], [], []
    evmap_all = {}
    t0 = time.time()
    for fi, path in enumerate(args.files):
        import h5py
        with h5py.File(path, "r") as f:
            cols_all = list(f["events"].dtype.names)
        cm2 = resolve_columns(cols_all, cfg["columns"], 2)
        cm3 = resolve_columns(cols_all, cfg["columns"], 3)
        dt_col = cm3.get("decay_type")
        run_col, evt_col = cm3.get("run"), cm3.get("event")
        evmap = event_cols(cols_all)
        evmap_all.update(evmap)
        iso = {k: find_col(cols_all, c) for k, c in ISO_COLS.items()}
        iso = {k: v for k, v in iso.items() if v}
        need = sorted({c for c in (F.raw_columns(cm3, 3) + [dt_col] + list(evmap.values()) + list(iso.values())) if c})
        chunks = _load_h5_chunks(path, need, args.chunk if not args.max_rows else min(args.chunk, args.max_rows), dt_col,
                                 [int(sp["decay_type"]), int(dc["dstar_decay_type"])])
        n_d0 = n_ds = n_pairs = 0
        for ic, raw in enumerate(chunks):
            if args.max_rows and ic * args.chunk >= args.max_rows:
                break
            # tabela zdarzeń (pierwszy wiersz każdego (run, event) w porcji; duplikaty między porcjami usuwane na końcu)
            E = raw[[run_col, evt_col] + list(evmap.values())].drop_duplicates([run_col, evt_col])
            E = E.rename(columns={run_col: "run", evt_col: "event", **{v: k for k, v in evmap.items()}})
            E["file_idx"] = fi
            ev_parts.append(E)
            r0 = raw[raw[dt_col] == int(sp["decay_type"])].reset_index(drop=True)
            r1 = raw[raw[dt_col] == int(dc["dstar_decay_type"])].reset_index(drop=True)
            if len(r0) == 0:
                continue
            L = F.to_logical(r0, cm2, 2)
            L["weight"] = L["weight"].fillna(1.0)
            lxy_ok = L["Lxy"].to_numpy(float) >= args.lxy_min
            L, r0 = L[lxy_ok].reset_index(drop=True), r0[lxy_ok].reset_index(drop=True)
            if len(L) == 0:
                continue
            L = F.compute(L, sp, units)
            L["row"] = r0["_local_row"].to_numpy()
            L["file_idx"] = fi
            for k, c in iso.items():
                L[k] = r0[c].to_numpy(np.float32) * float((units or {}).get("pt", 1.0))
            for k, c in evmap.items():          # kolumny zdarzeniowe niesione na wierszu (MC: brak unikalnego (run, event))
                L[f"e_{k}"] = r0[c].to_numpy()
            ecols = [f"e_{k}" for k in evmap]
            mark_swaps(L, sp, cfg)
            L["m_swap"] = swapped_mass_sp(L, sp)
            L["score"], L["fold"] = ens.score(L)
            n_d0 += len(L)
            m = L["m"].to_numpy()
            aw = (m >= args.all_lo) & (m <= args.all_hi) & keep_hash(L["run"].to_numpy(), L["event"].to_numpy(), args.all_frac)
            A = L.loc[aw, ["run", "event", "file_idx", "row", "score", "pt", "m", "truth", "weight"] + list(iso) + ecols]
            all_parts.append(A)
            if len(r1) == 0:
                continue
            win = np.abs(m - m0) < M_WIN
            ds = F.to_logical(r1, cm3, 3)
            ds = ds[(ds["dm"] > 0) & (ds["dm"] < DM_MAX) & (ds["Lxy"].to_numpy(float) >= args.lxy_min)]
            n_ds += len(ds)
            ds = ds.assign(m0_from=ds["m"] - ds["dm"], row_ds=r1["_local_row"].to_numpy()[ds.index], truth_ds=ds["truth"],
                           pt_ds=ds["pt"], eta_ds=ds["eta"], pis_pt=ds["d_pt2"] * float((units or {}).get("d_pt", 1.0)))[
                ["run", "event", "m0_from", "pt_ds", "eta_ds", "dm", "row_ds", "truth_ds", "pis_pt"]]
            Lk = L.loc[win, ["run", "event", "m", "pt", "eta"]].assign(_i=np.nonzero(win)[0])
            j = Lk.merge(ds, on=["run", "event"])
            ok = (np.abs(j["m"] - j["m0_from"]) < dc["match_tol_m"]) & (j["pt"] <= j["pt_ds"] * 1.001) & \
                 (j["pt"] > 0.8 * j["pt_ds"]) & (np.abs(j["eta"] - j["eta_ds"]) < 0.1)
            j = j[ok].reset_index(drop=True)
            if len(j) == 0:
                continue
            base = L.iloc[j["_i"].to_numpy()].reset_index(drop=True)
            keep = ["run", "event", "file_idx", "row", "m", "pt", "eta", "phi", "score", "fold", "truth", "weight", "m_swap"] \
                + list(iso) + ecols + feats
            part = pd.concat([base[keep], j[["row_ds", "dm", "truth_ds", "pt_ds", "eta_ds", "pis_pt"]]], axis=1)
            out_parts.append(part)
            n_pairs += len(part)
        print(f"[{args.tag}] {os.path.basename(path)}: D0 {n_d0}  D* {n_ds}  pary {n_pairs}  {time.time() - t0:.0f}s", flush=True)

    ints = ("run", "event", "file_idx", "row", "row_ds", "fold")
    P = pd.concat(out_parts, ignore_index=True) if out_parts else pd.DataFrame(columns=["run", "event", "dm", "score", "pt", "m"])
    def _np(X, c):
        if c in ints:
            return X[c].to_numpy(dtype=np.int64)
        if c.startswith("e_tp_") or c == "e_grl":
            return X[c].to_numpy(dtype=np.int8)
        return X[c].to_numpy(dtype=np.float32)
    out = {c: _np(P, c) for c in P.columns}
    A = pd.concat(all_parts, ignore_index=True) if all_parts else pd.DataFrame(columns=["run", "event", "file_idx", "row", "score",
                                                                                         "pt", "m", "truth", "weight"])
    for c in A.columns:
        out[f"all_{c if c != 'weight' else 'w'}"] = _np(A, c)
    EV = pd.concat(ev_parts, ignore_index=True) if ev_parts else pd.DataFrame(columns=["run", "event", "file_idx"])
    EV = EV.drop_duplicates(["file_idx", "run", "event"] if args.tag != "data" else ["run", "event"]).reset_index(drop=True)
    if "lead_pt" in EV and len(EV) and np.nanmedian(EV["lead_pt"].to_numpy(float)) > 200:
        EV["lead_pt"] = EV["lead_pt"] / 1000.0      # MeV (ROOT/charmpp v2.0) → GeV
    for c in EV.columns:
        v = EV[c].to_numpy()
        out[f"ev_{c}"] = v.astype(np.int8) if c.startswith("tp_") or c == "grl" else \
            v.astype(np.int64) if c in ("run", "event", "file_idx", "nSelBL", "npu", "lb", "ntrk500") else v.astype(np.float32)
    out["meta"] = np.array(json.dumps({"tag": args.tag, "files": [os.path.basename(p) for p in args.files], "features": feats,
                                       "thresholds": meta.get("thresholds"), "m0": m0, "dm_max": DM_MAX, "m_win": M_WIN,
                                       "all_window": [args.all_lo, args.all_hi], "all_frac": args.all_frac,
                                       "event_columns": sorted(evmap_all), "n_pairs": int(len(P)), "n_all": int(len(A)),
                                       "n_events": int(len(EV)), "btrig": "1.0"}))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez_compressed(args.out, **out)
    print(f"[{args.tag}] zapisano {args.out}: pary {len(P)}, D0 all {len(A)}, zdarzenia {len(EV)}, "
          f"kolumny zdarzeniowe {sorted(evmap_all)}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())

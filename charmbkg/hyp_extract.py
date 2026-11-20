"""HYP — ekstrakcja kandydatów D0 (obie hipotezy K/π) z h5 v1.2 do kompaktowego parquet (lxplus; jeden job = jedna próbka).
Wywołanie: `python3 -m charmbkg.hyp_extract` albo pkg/hyp_extract/extract_d0hyp.py.

  python3 -m charmbkg.hyp_extract --sample fd_nd --h5 '/eos/.../charm_h5_v12/parts/fd_nd/*_none.h5' --out /eos/.../fd_nd.parquet \
      --scores /eos/.../runs/lxplus_v11/scores/D0_H/signal_mc.h5 --models /eos/.../runs/lxplus_v11/models/D0_H [--frac 1] [--max-gb 1.5]

Wiersze: decayType 10, m w `species.D0.mass_range` i pT 3–30 (= `samples.in_range`, jak `charmbkg apply`). `row` = wiersz
GLOBALNY w posortowanej liście plików (jak `charmbkg.io.open_events` / `config.sample_paths`: glob sortowany, bez *.part).
Cechy H: `charmbkg.features.compute` (te same co w treningu D0_H). Score nominalny D0_H:
  --score-mode join  — dołącz po `row` z scores/D0_H/<klucz>.h5 (kontrola: m identyczne bit w bit w każdym wierszu, offsety
                       plików i nazwy plików zgodne); niezgodność = błąd,
  --score-mode model — licz modelem foldów `--models` (charmbkg.model.Ensemble; fold po (run, event) jak w apply),
  --score-mode auto  — join, a gdy niemożliwy/niezgodny → model (domyślnie); none — bez score (NaN).
Podpróbkowanie: --frac (całe zdarzenia, hash splitmix(run, event) → pary hipotez razem); --max-gb obniża frac tak, by
szacowany parquet ≤ limit (przebieg liczący 3 kolumny). DRYRUN=1 / --dry-run: tylko plan (pliki, wiersze, frac), bez zapisu.
T4 (09.11): tag D*⁺→D0π⁺_soft (`charmbkg.dstar_tag`, --dstar-tag auto|off): kolumny dstar_tag, dm, soft_pi_charge, hypo_from_tag
(0 zgodny z zapachem / 1 refleksja / −1 brak tagu), dstar_n, dstar_conflict, dstar_row, pis_pt — dopasowanie po trackId córek
w tym samym zdarzeniu; poprawna hipoteza K/π także na DANYCH. Statystyki (i na MC macierz hypo_from_tag × hypo_swap) w meta["dstar"].
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd

BYTES_PER_ROW = 200          # parquet zstd, ~65 kolumn (pomiar na toy: ~175 B/wiersz + T4 ~8 kolumn, głównie stałe) — zapas
EXTRA_V12 = {                # v1.2 → nazwy w parquet (brak w pliku → kolumna pominięta)
    "d_tpdg": ["Dmeson_daughter{i}_truthPdgId"], "d_tparent": ["Dmeson_daughter{i}_truthParentPdgId"],
    "d_charge": ["Dmeson_daughter{i}_charge"], "d_pdg": ["Dmeson_daughter{i}_pdgId"],
    "D_tpdg": ["Dmeson_truthPdgId"], "n_truth_ok": ["Dmeson_n_daughter_truth_ok"],
}
KEEP = ["sample", "row", "file", "evblk", "run", "event", "weight", "m", "pt", "eta", "phi", "truth", "hypo_swap", "score", "score_src", "fold",
        "cos_theta_star", "dstar_tag", "dm", "soft_pi_charge", "hypo_from_tag", "dstar_n", "dstar_conflict", "dstar_row", "pis_pt",
        ] + [f"{k}{i}" for k in ("d_pt", "d_eta", "d_phi", "d_d0", "d_z0", "d_dedx", "d_npts", "d_trk") for i in (0, 1)]
SAMPLE_KEY = {"fd_nd": "signal_mc", "mb_nd": "incl_mc", "data_428770": "data_none"}


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def files_of(spec: str) -> list[str]:
    out = []
    for it in spec.split(","):
        it = os.path.expanduser(it.strip())
        out += [g for g in sorted(glob.glob(it)) if not g.endswith(".part")] if any(c in it for c in "*?[") else [it]
    return list(dict.fromkeys(out))


def keep_event(run, event, frac: float) -> np.ndarray:
    if frac >= 1:
        return np.ones(len(event), bool)
    from charmbkg.hyp import _mix64
    x = _mix64((np.asarray(run, np.int64).astype(np.uint64) * np.uint64(0xD6E8FEB86659FD93)) ^ np.asarray(event, np.int64).astype(np.uint64))
    return (x % np.uint64(1_000_000)).astype(np.int64) < int(round(frac * 1_000_000))


def count_selected(files, cm, ds, sp, chunk):
    from charmbkg.io import open_events
    from charmbkg.samples import pt_range
    n = 0
    lo, hi = sp["mass_range"]
    plo, phi = pt_range()
    with open_events(files, ds) as E:
        for a in range(0, E.shape[0], chunk):
            x = E.fields([cm["decay_type"], cm["m"], cm["pt"]])[a:min(E.shape[0], a + chunk)]
            dt, m, pt = (np.asarray(x[c]) for c in (cm["decay_type"], cm["m"], cm["pt"]))
            n += int(((dt == 10) & (m >= lo) & (m < hi) & (pt >= plo) & (pt < phi)).sum())
    return n


class ScoreJoin:
    """scores/D0_H/<klucz>.h5 (charmbkg.store, kolumna-na-dataset; `row` rosnąco jak w apply)."""

    def __init__(self, path, files, offsets):
        import h5py
        self.f = h5py.File(path, "r")
        self.row = self.f["row"][:].astype(np.int64)
        if np.any(np.diff(self.row) < 0):
            raise ValueError("scores: `row` nie rosnące")
        self.problems = []
        src = self.f.attrs.get("source")
        if src is not None:
            src = src.decode() if isinstance(src, bytes) else str(src)
            names = [os.path.basename(s) for s in src.split(";")]
            if names != [os.path.basename(q) for q in files]:
                self.problems.append(f"lista plików scores ({len(names)}) ≠ lista h5 v1.2 ({len(files)}) lub inna kolejność/nazwy")
        fo = self.f.attrs.get("file_offsets")
        if fo is not None and json.loads(fo) != [int(x) for x in offsets]:
            self.problems.append("file_offsets scores ≠ offsety h5 v1.2 (inne pliki/wiersze)")
        self.checked = self.bad = self.missing = 0

    def get(self, rows, m):
        i = np.searchsorted(self.row, rows)
        i = np.clip(i, 0, len(self.row) - 1)
        hit = self.row[i] == rows
        out = np.full(len(rows), np.nan, np.float32)
        fold = np.full(len(rows), -1, np.int8)
        if hit.any():
            a, b = int(i[hit].min()), int(i[hit].max()) + 1      # wiersze porcji są ciągłe w scores → jeden odczyt wycinka
            sc = self.f["score"][a:b]
            mm = self.f["m"][a:b]
            fd = self.f["fold"][a:b] if "fold" in self.f else np.full(b - a, -1)
            j = i[hit] - a
            out[hit], fold[hit] = sc[j], fd[j]
            d = np.abs(mm[j].astype(np.float64) - np.asarray(m, np.float64)[hit])
            self.checked += int(hit.sum())
            self.bad += int((d > 1e-9).sum())
        self.missing += int((~hit).sum())
        return out, fold

    def close(self):
        self.f.close()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True, help="fd_nd | mb_nd | data_428770 (prefiks nazw)")
    ap.add_argument("--h5", required=True, help="glob/lista plików h5 v1.2 (po przecinku)")
    ap.add_argument("--out", required=True, help="plik .parquet (obok: <out>.meta.json)")
    ap.add_argument("--config", default=None, help="config charmbkg (domyślnie pakietowy default.yaml)")
    ap.add_argument("--scores", default=None)
    ap.add_argument("--models", default=None)
    ap.add_argument("--score-mode", default="auto", choices=["auto", "join", "model", "none"])
    ap.add_argument("--frac", type=float, default=1.0)
    ap.add_argument("--max-gb", type=float, default=0.0, help="limit rozmiaru parquet [GB]; 0 = bez limitu")
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--max-rows", type=int, default=0, help="test: tylko pierwsze N wierszy")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--dstar-tag", default="auto", choices=["auto", "off"], help="T4: tag D* po trackId (auto = gdy kolumny v1.2 są)")
    ap.add_argument("--dstar-dm-max", type=float, default=0.165, help="T4: górna granica Δm D* [GeV] dla tagu")
    ap.add_argument("--dstar-only", action="store_true", help="T4: zapisz tylko wiersze z dstar_tag=1 (pełna lumi, mały parquet)")
    a = ap.parse_args(argv)
    dry = a.dry_run or os.environ.get("DRYRUN") == "1"

    from charmbkg import config as C
    from charmbkg import features as F
    from charmbkg.io import list_columns, n_rows, open_events, resolve_columns
    from charmbkg.samples import colmap, fold_of, in_range, pt_range as _pt_range
    cfg = C.load(a.config)
    sp = cfg["species"]["D0"]
    ds = cfg["io"]["dataset"]
    files = files_of(a.h5)
    if not files or not all(os.path.exists(q) for q in files):
        raise SystemExit(f"brak plików: {a.h5} → {files[:3]}")
    offsets = np.concatenate([[0], np.cumsum([n_rows(q, ds) for q in files])]).astype(np.int64)
    ntot = int(offsets[-1]) if not a.max_rows else min(int(offsets[-1]), a.max_rows)
    cm = colmap(cfg, files[0], 2)
    for need in ("decay_type", "m", "pt", "eta", "event"):
        if not cm.has(need):
            raise SystemExit(f"brak kolumny {need} ({cfg['columns'][need]})")
    ex = resolve_columns(list_columns(files[0], ds), EXTRA_V12, 2)
    v12 = cm.has("hypo_swap") and cm.has("d_trk0")
    log(f"{a.sample}: {len(files)} plików, {ntot} wierszy; v1.2 (hypo_swap, trackId): {v12}; dodatkowe: "
        f"{[k for k, v in ex.mapping.items() if v]}")
    # --- T4: kolumny D* (córki 0..2: trackId, ładunek, pT π_soft)
    from charmbkg import dstar_tag as DT
    dx = resolve_columns(list_columns(files[0], ds), {**DT.DSTAR_V12, "ds_pt": cfg["columns"]["d_pt"]}, 3)
    tag_on = a.dstar_tag != "off" and cm.has("d_trk0") and cm.has("d_trk1") and cm.has("dm") and dx.has("ds_trk0") \
        and dx.has("ds_trk1") and (dx.has("ds_q2") or cm.has("charge"))
    dstar_cols = [c for c in [dx.get("ds_trk0"), dx.get("ds_trk1"), dx.get("ds_trk2"), dx.get("ds_q0"), dx.get("ds_q1"),
                              dx.get("ds_q2"), dx.get("ds_pt2"), cm.get("dm"), cm.get("charge")] if c] if tag_on else []
    log(f"T4 tag D*: {'TAK' if tag_on else 'NIE (brak kolumn v1.2 lub --dstar-tag off)'}; kolumny {dstar_cols}")
    if a.dstar_only and not tag_on:
        raise SystemExit("--dstar-only bez tagu D* (brak kolumn v1.2 trackId/Δm albo --dstar-tag off)")
    frac = float(a.frac)
    n_sel = None
    if a.max_gb > 0 and not a.dstar_only:   # --dstar-only: wiersze otagowane to ~1 % okna → limit nie obniża frac
        n_sel = count_selected(files, cm, ds, sp, a.chunk)
        cap = a.max_gb * 1e9 / BYTES_PER_ROW
        if n_sel * frac > cap:
            log(f"limit {a.max_gb} GB: {n_sel} wierszy w oknie × frac {frac} > {cap:.0f} → frac = {cap / n_sel:.4f}")
            frac = cap / n_sel
    meta = {"sample": a.sample, "files": [os.path.basename(q) for q in files], "file_offsets": offsets.tolist(), "rows_total": ntot,
            "n_selected_before_frac": n_sel, "frac": frac, "frac_requested": a.frac, "v12": bool(v12), "mass_range": sp["mass_range"],
            "pt_range": list(_pt_range()), "score_mode_requested": a.score_mode,
            "dstar": {"enabled": bool(tag_on), "dm_max": a.dstar_dm_max, "match": "evblk + {trackId0, trackId1}",
                      "only_tagged": bool(a.dstar_only)}}
    if a.sample.startswith("data"):
        meta["data_frac"] = frac
    # --- score: plan
    sj, ens, mode = None, None, a.score_mode
    if mode in ("auto", "join"):
        if a.scores and os.path.exists(a.scores):
            sj = ScoreJoin(a.scores, files, offsets)
            if sj.problems:
                msg = "; ".join(sj.problems)
                if mode == "join":
                    raise SystemExit(f"score join niemożliwy: {msg}")
                log("score join → fallback model:", msg)
                sj.close()
                sj, mode = None, "model"
            else:
                mode = "join"
        elif mode == "join":
            raise SystemExit(f"brak pliku scores {a.scores}")
        else:
            mode = "model"
    if mode == "model":
        if a.models and os.path.exists(os.path.join(a.models, "meta.json")):
            from charmbkg.model import Ensemble
            ens = Ensemble(a.models)
        else:
            log("UWAGA: brak modeli D0_H → score = NaN")
            mode = "none"
    meta["score_mode"] = mode
    if a.models and os.path.exists(os.path.join(a.models, "meta.json")):
        meta["thresholds"] = json.load(open(os.path.join(a.models, "meta.json"))).get("thresholds")
    log(f"score: {mode}; frac={frac:.4f}")
    if dry:
        print(json.dumps(meta, indent=1)[:3000])
        return 0
    # --- przebieg
    import pyarrow as pa
    import pyarrow.parquet as pq
    feats = F.feature_list(cfg, "H")
    raw_cols = list(dict.fromkeys(F.raw_columns(cm, 2) + [v for v in ex.mapping.values() if v] + dstar_cols))
    tmp = a.out + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    W, nout, nchk, prev = None, 0, 0, None
    from charmbkg.hyp import event_blocks
    hs_count = {}
    ds_carry, ds_stat, ds_conf, ds_conf_sw = None, {}, {}, {}
    dm_sig = (cfg.get("dstar_control") or {}).get("dm_signal", [0.1435, 0.1475])
    upt = float((cfg.get("units") or {}).get("d_pt", 1.0))
    with open_events(files, ds) as E:
        for start in range(0, ntot, a.chunk):
            stop = min(ntot, start + a.chunk)
            x = E.fields(raw_cols)[start:stop]
            raw = pd.DataFrame({c: np.asarray(x[c]) for c in raw_cols})
            fi = np.searchsorted(offsets, np.arange(start, stop), side="right") - 1        # indeks pliku (wiersz globalny)
            rr = raw[cm["run"]].to_numpy() if cm.has("run") else np.zeros(len(raw), np.int64)
            evb, prev = event_blocks(rr.astype(np.int64) * 4096 + fi, raw[cm["event"]].to_numpy(), 0, prev)  # blok = (plik, run, event) ciągły
            raw = raw.assign(_evblk=evb, _file=fi)
            T = None
            if tag_on:   # T4: tabela D* tej porcji (+ D* ostatniego zdarzenia poprzedniej porcji)
                dsr = raw[raw[cm["decay_type"]].to_numpy() == DT.DSTAR_DT]
                g = lambda k: dsr[dx[k]].to_numpy() if dx.has(k) else None  # noqa: E731
                qs = g("ds_q2")
                qc = dsr[cm["charge"]].to_numpy() if cm.has("charge") else None
                if qs is None:
                    qs = qc
                elif qc is not None:
                    qs = np.where(np.abs(qs) == 1, qs, qc)
                Tn = DT.dstar_table(dsr["_evblk"].to_numpy(), dsr.index.to_numpy() + start, g("ds_trk0"), g("ds_trk1"), g("ds_trk2"),
                                    dsr[cm["dm"]].to_numpy(), qs, g("ds_q0"), g("ds_q1"),
                                    None if g("ds_pt2") is None else g("ds_pt2") * upt, a.dstar_dm_max)
                DT.merge_counts(ds_stat, {"n_dstar_rows": len(dsr), "n_dstar_ok": len(Tn),
                                          "n_dstar_k_consistent": int((Tn["ds_q0"].to_numpy() == -Tn["ds_qs"].to_numpy()).sum()),
                                          "n_dstar_q0_known": int(np.isfinite(Tn["ds_q0"].to_numpy()).sum()),
                                          "n_dstar_dm_signal": int(((Tn["ds_dm"] >= dm_sig[0]) & (Tn["ds_dm"] < dm_sig[1])).sum())})
                T = Tn if ds_carry is None or not len(ds_carry) else pd.concat([ds_carry, Tn], ignore_index=True)
                ds_carry = DT.carry_last_event(T, evb[-1] if len(evb) else None)
            raw = raw[raw[cm["decay_type"]].to_numpy() == 10]
            if not len(raw):
                continue
            L = F.to_logical(raw.reset_index(drop=True), cm, 2)
            if not cm.has("run"):
                L["run"] = 0
            L["weight"] = L["weight"].fillna(1.0)
            L = F.compute(L, {**sp, "n_daughters": 2}, cfg.get("units"))
            L["row"] = raw.index.to_numpy() + start
            L["evblk"] = raw["_evblk"].to_numpy()
            L["file"] = raw["_file"].to_numpy().astype(np.int16)
            for k, v in ex.mapping.items():
                if v:
                    L[k] = raw[v].to_numpy()
            sel = in_range(L, sp).to_numpy() & keep_event(L["run"].fillna(0).to_numpy(), L["event"].to_numpy(), frac)
            L = L[sel].reset_index(drop=True)
            if not len(L):
                continue
            if tag_on:   # T4
                tg = DT.tag_d0(L["evblk"].to_numpy(), L["d_trk0"].to_numpy(), L["d_trk1"].to_numpy(), T,
                               L["d_charge0"].to_numpy() if "d_charge0" in L else None)
            else:
                tg = DT.empty_tags(len(L))
            for k, v in tg.items():
                L[k] = v
            t1 = tg["dstar_tag"] == 1
            DT.merge_counts(ds_stat, {"n_rows_tagged": int(t1.sum()), "n_rows_conflict": int(tg["dstar_conflict"].sum()),
                                      "n_rows_tag_ok": int((tg["hypo_from_tag"] == 0).sum()),
                                      "n_rows_tag_refl": int((tg["hypo_from_tag"] == 1).sum()),
                                      "n_rows_tag_dm_signal": int((t1 & (tg["dm"] >= dm_sig[0]) & (tg["dm"] < dm_sig[1])).sum())})
            if "hypo_swap" in L:   # MC: tag vs prawda
                DT.merge_counts(ds_conf, DT.confusion(tg["hypo_from_tag"], L["hypo_swap"].to_numpy()))
                DT.merge_counts(ds_conf_sw, DT.confusion(tg["hypo_from_tag"], L["hypo_swap"].to_numpy(),
                                                         (tg["dm"] >= dm_sig[0]) & (tg["dm"] < dm_sig[1])))
            if a.dstar_only:
                L = L[t1].reset_index(drop=True)
                if not len(L):
                    continue
            if mode == "join":
                L["score"], L["fold"] = sj.get(L["row"].to_numpy(np.int64), L["m"].to_numpy())
            elif mode == "model":
                L["score"], L["fold"] = ens.score(L)
            else:
                L["score"], L["fold"] = np.float32(np.nan), np.int8(-1)
            if mode != "model" and "fold" in L and (L["fold"] < 0).all():
                L["fold"] = fold_of(L["run"].fillna(0).astype(np.int64), L["event"].astype(np.int64), 5)
            L["score_src"] = {"join": 1, "model": 2}.get(mode, 0)
            L["sample"] = a.sample
            if "hypo_swap" not in L:
                L["hypo_swap"] = -1
            for v, c in zip(*np.unique(L["hypo_swap"].to_numpy(), return_counts=True)):
                hs_count[int(v)] = hs_count.get(int(v), 0) + int(c)
            cols = [c for c in KEEP if c in L] + [c for c in L.columns if c not in KEEP and (c in feats or c in ex.mapping
                                                                                            or c[:-1] in ex.mapping)]
            O = L[list(dict.fromkeys(cols))].copy()
            for c in O.columns:
                if O[c].dtype == np.float64 and c != "m":
                    O[c] = O[c].astype(np.float32)
                elif O[c].dtype == np.int64 and c not in ("row", "event", "evblk", "d_trk0", "d_trk1", "dstar_row"):
                    O[c] = O[c].astype(np.int32)
            t = pa.Table.from_pandas(O, preserve_index=False)
            if W is None:
                # 09.11: pyarrow z LCG_110_swan nie ma zstd (ArrowNotImplementedError) → pierwszy dostępny kodek
                codec = next((c for c in ("zstd", "gzip", "snappy") if pa.Codec.is_available(c)), "none")
                W = pq.ParquetWriter(tmp, t.schema, compression=codec)
            elif t.schema != W.schema:
                t = t.cast(W.schema)
            W.write_table(t)
            nout += len(O)
            nchk += 1
            log(f"  [{stop}/{ntot}] +{len(O)} (razem {nout})")
    if sj is not None:
        meta["score_join"] = {"checked": sj.checked, "mass_mismatch": sj.bad, "row_missing": sj.missing}
        sj.close()
        if sj.bad or sj.missing:
            if W is not None:
                W.close()
            os.remove(tmp)
            raise SystemExit(f"score join: {sj.bad} niezgodnych mas, {sj.missing} brakujących wierszy na {sj.checked} → "
                             "to nie te same wiersze co scores (uruchom z --score-mode model)")
    if W is None:
        raise SystemExit("0 wierszy po selekcji")
    W.close()
    os.replace(tmp, a.out)
    meta["dstar"].update({"counts": ds_stat, "confusion_vs_hypo_swap": ds_conf, "confusion_vs_hypo_swap_dm_signal": ds_conf_sw,
                          "dm_signal": dm_sig})
    meta.update({"n_rows_out": nout, "bytes": os.path.getsize(a.out), "hypo_swap_counts": hs_count,
                 "columns": pq.read_schema(a.out).names, "time": time.strftime("%Y-%m-%dT%H:%M:%S")})
    with open(a.out[:-len(".parquet")] + ".meta.json" if a.out.endswith(".parquet") else a.out + ".meta.json", "w") as f:
        json.dump(meta, f, indent=1)
    log(f"OK {a.out}: {nout} wierszy, {meta['bytes'] / 1e6:.1f} MB, hypo_swap {hs_count}; D* tag {ds_stat}; tag×prawda {ds_conf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

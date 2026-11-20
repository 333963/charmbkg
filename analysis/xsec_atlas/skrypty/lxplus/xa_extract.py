#!/usr/bin/env python3
"""XSEC_atlas (20.11.2026, projekt B) — tabele kandydatow inkluzywnych D0 / D+ ze zmiennymi ciec (selekcja cieciowa w stylu ATLAS, KRYTERIA.md cz. B).
Jedna grupa plikow h5 v2.0 jednej probki -> jeden npz z tabelami <gat>__t_* (format INCAL `table`, bez kolumn ROOT):
  wiersz zapisany, gdy ok (D+: weto D*) i [ preselekcja luzna (nadzbior ciec i ich wariantow)  LUB  score drzewa H >= prog WP70 binu pT
  (D0: L_xy >= 0,25 mm) LUB (tylko MC) prawda sygnalu 0/1 ].
Tylko import kodu INCAL / dnet1 (bez zmian w nich). Bez AFP, tylko pojedyncze D.
Uzycie: xa_extract.py --sample S --k K --n N --h5dir DIR --models DIR --out OUT.npz [--max-files n]"""
import argparse, glob, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("INC") or os.path.join(os.path.dirname(HERE), "incal"))
import inc_extract as IE  # noqa: E402
import pandas as pd  # noqa: E402

DF, cfg = IE.DF, IE.cfg
LOOSE = {"D0": {"lxy": 0.33, "ptmin": 0.6}, "Dplus": {"lxy": 0.75, "ptmin": 0.6}}
LXY_READ = {"D0": 0.25, "Dplus": 0.10}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--k", type=int, default=0); ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--h5dir", required=True); ap.add_argument("--models", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--max-files", type=int, default=0)
    a = ap.parse_args(argv)
    t0 = time.time()
    is_mc = not a.sample.startswith("data")
    files = IE.group(sorted(glob.glob(os.path.join(a.h5dir, a.sample, "*_none.h5"))), a.k, a.n)
    if a.max_files:
        files = files[:a.max_files]
    if not files:
        raise SystemExit("brak plikow h5")
    models = {sp: IE.load_models(a.models, t) for sp, t in IE.TAGS.items()}
    parts, stat = {}, []
    for path in files:
        tf = time.time()
        st = {"file": os.path.basename(path)}
        for spn in ("D0", "Dplus"):
            nd = int(cfg["species"][spn]["n_daughters"])
            L, _ = IE.read_species(path, spn, 3, 30, LXY_READ[spn])
            n = len(L)
            X = np.stack([L[c].to_numpy(np.float64) for c in DF.H], axis=1).astype(np.float32) if n else np.zeros((0, len(DF.H)), np.float32)
            L["score"] = IE.score_rows(models[spn], L["uid"].to_numpy(np.uint64), X) if n else np.zeros(0, np.float32)
            pt, sc = L["pt"].to_numpy(float), L["score"].to_numpy(float)
            thr = np.full(n, np.inf)
            for name, lo, hi in IE.PTB["d150"]:
                thr[(pt >= lo) & (pt < hi)] = IE.THR[spn][name][0]
            ptmin = np.min(np.stack([L[f"d_pt{i}"].to_numpy(float) for i in range(nd)], 1), 1) if n else np.zeros(0)
            lxy = L["Lxy"].to_numpy(float)
            loose = (lxy > LOOSE[spn]["lxy"]) & (ptmin > LOOSE[spn]["ptmin"])
            cls70 = sc >= thr
            keep = L["ok"].to_numpy(bool) & (loose | cls70)
            if is_mc:
                keep |= np.isin(L["truth"].to_numpy(float), DF.SIGNAL_CATS)
            Lw = L[keep]
            M = pd.DataFrame({c: np.full(len(Lw), np.nan, np.float32) for c in IE.RCOLS + ["anc"]})
            T = IE.table(Lw, M, spn, is_mc, f"{spn}__t_")
            for c in IE.RCOLS + ["anc", "pvz_h5", "dm"]:
                T.pop(f"{spn}__t_{c}", None)
            for i in range(nd):
                T.pop(f"{spn}__t_d{i}_z0", None)
            st[spn] = {"rows": int(n), "ok": int(L["ok"].sum()), "loose": int((L["ok"].to_numpy(bool) & loose).sum()), "cls70": int((L["ok"].to_numpy(bool) & cls70).sum()), "kept": int(keep.sum())}
            for k, v in T.items():
                parts.setdefault(k, []).append(v)
        st["sec"] = round(time.time() - tf, 1)
        stat.append(st)
        print(f"[{a.sample} {a.k}/{a.n}] {st}", flush=True)
    out = {k: np.concatenate(v) for k, v in parts.items()}
    out["meta"] = np.array(json.dumps({"sample": a.sample, "k": a.k, "n": a.n, "files": stat, "loose": LOOSE, "lxy_read": LXY_READ, "thr": IE.THR,
                                       "tags": IE.TAGS, "mrange": IE.MRANGE, "xsatlas": "1.0"}))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.sample} {a.k}/{a.n}] zapisano {a.out} ({os.path.getsize(a.out) / 1e6:.1f} MB): { {k: len(v) for k, v in out.items() if k.endswith('_m')} }  {time.time() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

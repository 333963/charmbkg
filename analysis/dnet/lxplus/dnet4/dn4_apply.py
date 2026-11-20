#!/usr/bin/env python3
"""DNET etap 4 — zastosowanie zapisanych modeli foldów (dn4_train.py) do JEDNEJ próbki: tabela dn_<próbka>.npz (dn4_extract.py) + tory
otoczenia dn2_<próbka>_*.npz (../dnet2/dn2_extract.py --species D0,Dplus) → zwarty plik wierszy:
  <gatunek>__uid (u8), m, pt (f4), ok (bool; D+: weto D*), cls (u1: bit0 GRL, bit1 HLT_mb_mbts_L1MBTS_1, bit2 HLT_noalg_L1J12, bit3 dowolny
  łańcuch MBTS), truth (i2), swap (i1), w (f4), m_swap (f4; D0), lg [n, nmodeli] (logit score, float16), meta (json: tags per gatunek).
Score = model foldu fold_of(event_uid) (dn4_lib.score_rows). Bez selekcji AFP, bez tagu; tylko pojedyncze D.
Użycie: dn4_apply.py --sample S --data DIR --env DIR --models DIR --out OUT.npz [--species D0,Dplus] [--tags REGEX] [--threads n]"""
import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dn4_lib as L4  # noqa: E402
import dn_feat as DF  # noqa: E402
import dn2_train as D2T  # noqa: E402


def load_one(ddir, edir, s, sp):
    """Wiersze a_ próbki s z torami otoczenia (jak dn2_train.load dla jednej próbki) + kolumny zdarzenia/trygera."""
    z = np.load(os.path.join(ddir, f"dn_{s}.npz"), allow_pickle=False)
    k = f"{sp}__a_"
    if k + "m" not in z.files:
        return None
    ix, E = [], []
    for f in sorted(glob.glob(os.path.join(edir, f"dn2_{s}_*.npz"))):
        e = np.load(f, allow_pickle=False)
        if k + "idx" in e.files:
            ix.append(e[k + "idx"]); E.append(e[k + "E"])
    if not ix:
        return None
    ix, E = np.concatenate(ix), np.concatenate(E)
    o = np.argsort(ix)
    ix, E = ix[o], E[o]
    n_full = len(z[k + "m"])
    hc = [DF.ALLF[sp].index(f) for f in DF.H]
    T = {c: z[k + c][ix] for c in D2T.KEYS if k + c in z.files}
    T["X"] = T["X"][:, hc]
    T["E"] = E
    n = len(ix)
    tp = lambda c: (z[k + "e_tp_" + c][ix] == 1) if k + "e_tp_" + c in z.files else np.zeros(n, bool)  # noqa: E731
    grl = (z[k + "e_grl"][ix] == 1) if k + "e_grl" in z.files else np.ones(n, bool)
    mbts = np.zeros(n, bool)
    for c in DF.MBTS_CHAINS:
        mbts |= tp(c)
    T["cls"] = (grl * 1 + tp("HLT_mb_mbts_L1MBTS_1") * 2 + tp(DF.J12) * 4 + mbts * 8).astype(np.uint8)
    return T, n_full


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--data", required=True); ap.add_argument("--env", required=True)
    ap.add_argument("--models", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--species", default="D0,Dplus"); ap.add_argument("--tags", default="")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("DN_CPUS", "4")))
    a = ap.parse_args(argv)
    t0 = time.time()
    metas = {}
    for f in sorted(glob.glob(os.path.join(a.models, "*", "meta.json"))):
        t = os.path.basename(os.path.dirname(f))
        if a.tags and not re.search(a.tags, t):
            continue
        metas[t] = json.load(open(f))
    if not metas:
        raise SystemExit(f"brak modeli w {a.models} (filtr „{a.tags}”)")
    out, stat = {}, {"sample": a.sample, "tags": {}, "rows": {}, "sec": {}}
    for sp in a.species.split(","):
        r = load_one(a.data, a.env, a.sample, sp)
        tags = sorted(t for t, m in metas.items() if m["sp"] == sp)
        stat["tags"][sp] = tags
        if r is None or not tags:
            print(f"[{a.sample}] {sp}: brak wierszy albo modeli — pomijam", flush=True)
            stat["tags"][sp] = []
            continue
        T, n_full = r
        n = len(T["m"])
        if n != n_full:
            print(f"UWAGA [{a.sample}] {sp}: tory otoczenia dla {n} z {n_full} wierszy", flush=True)
        lg = np.zeros((n, len(tags)), np.float16)
        for j, t in enumerate(tags):
            tm = time.time()
            meta, ms = L4.load_models(a.models, t)
            if meta["model"] == "xgb":
                s_ = L4.score_rows("xgb", ms, T["uid"], T["X"])
            else:
                import torch
                torch.set_num_threads(max(1, a.threads))
                s_ = L4.score_rows("tr", ms, T["uid"], T["X"], T["Dd"], T["E"], T["eta"], T["phi"])
            if not np.isfinite(s_).all():
                raise RuntimeError(f"{t}: score NaN")
            lg[:, j] = L4.to_lg16(s_)
            stat["sec"][t] = round(time.time() - tm, 1)
            print(f"[{a.sample}] {t}: {n} wierszy, {time.time() - tm:.0f}s ({n / max(time.time() - tm, 1e-3):.0f}/s)", flush=True)
        k = sp + "__"
        out.update({k + "uid": T["uid"].astype(np.uint64), k + "m": T["m"].astype(np.float32), k + "pt": T["pt"].astype(np.float32),
                    k + "ok": T["ok"].astype(bool), k + "cls": T["cls"], k + "truth": T["truth"].astype(np.int16), k + "swap": T["swap"].astype(np.int8),
                    k + "w": T["w"].astype(np.float32), k + "lg": lg})
        if sp == "D0":
            out[k + "m_swap"] = T["m_swap"].astype(np.float32)
        stat["rows"][sp] = [int(n), int(n_full)]
    stat["total_sec"] = round(time.time() - t0, 1)
    out["meta"] = np.array(json.dumps(stat))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.sample}] zapisano {a.out}: {stat['rows']}, {time.time() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

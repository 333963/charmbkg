#!/usr/bin/env python3
"""sd_merge.py — suma części (npz z sd_reduce) per model → <out>/<model>.json. Użycie: sd_merge.py <SD_OUT> <out_dir> [model ...]
Per histogram: s = Σ po jobach, s2 = Σ kwadratów (→ błąd stat. z rozrzutu między jobami, uwzględnia korelacje D w zdarzeniu).
σ(hist) [mb] = s · sigma_mb / n_ok  (sigma_mb = σ_gen całego biegu, średnia ważona n_ok; dla doHard typ 3/4 = σ po wecie)."""
import glob, json, os, sys
import numpy as np


def merge(parts):
    S, S2 = {}, {}; nok = 0; sw = 0.0; se2 = 0.0; proc = {}; edges = None; meta0 = None; nj = 0; noks = []
    for p in parts:
        try: z = np.load(p, allow_pickle=False)
        except Exception as e: print("  POMIJAM", p, e); continue
        m = json.loads(str(z["__meta__"])); edges = edges or json.loads(str(z["__edges__"])); meta0 = meta0 or m
        nj += 1; n = m["n_ok"]; nok += n; noks.append(n); sw += m["sigma_gen_mb_total"] * n; se2 += (m["sigma_err_mb_total"] * n) ** 2
        for k, v in m["proc"].items():
            q = proc.setdefault(k, {"n": 0, "tagA": 0, "tagC": 0, "tagAC": 0, "hxi": [0] * 40})
            for a in ("n", "tagA", "tagC", "tagAC"): q[a] += v[a]
            q["hxi"] = [a + b for a, b in zip(q["hxi"], v["hxi"])]
        for k in z.files:
            if k.startswith("__"): continue
            v = z[k]
            if k not in S: S[k] = np.zeros_like(v); S2[k] = np.zeros_like(v)
            S[k] += v; S2[k] += v * v
    if not nj: return None
    return {"njobs": nj, "n_ok": nok, "n_ok_jobs": noks, "sigma_mb": sw / nok, "sigma_err_mb": se2 ** 0.5 / nok, "meta_first": {k: meta0[k] for k in meta0 if k != "proc"},
            "proc": proc, "edges": edges, "hists": {k: {"s": S[k].tolist(), "s2": S2[k].tolist()} for k in sorted(S)}}


def main():
    base, out = sys.argv[1], sys.argv[2]
    models = sys.argv[3:] or sorted(d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d)) and d != "merged")
    os.makedirs(out, exist_ok=True); summ = {}
    for mname in models:
        r = merge(sorted(glob.glob(os.path.join(base, mname, "part_*.npz"))))
        if r is None: print(f"{mname}: brak części"); continue
        r["model"] = mname
        json.dump(r, open(os.path.join(out, mname + ".json"), "w"))
        nt = sum(sum(h["s"]) for k, h in r["hists"].items() if k.endswith("|pt_tag"))
        summ[mname] = {"njobs": r["njobs"], "n_ok": r["n_ok"], "sigma_mb": r["sigma_mb"], "n_D_tag": nt}
        print(f"{mname}: {r['njobs']} części, {r['n_ok']} zdarzeń, σ = {r['sigma_mb']:.4g} mb, D fid. z tagiem {nt:.0f}")
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""EV + v1.2: przepisuje prawdę swapu K↔π w ev_fd.npz (v1.1, closest_mass) na Dmeson_hypo_swap z wyciągu HYP v1.2
(fd_nd.parquet; wiersze v1.2 = v1.1, H2). Join po (plik, wiersz lokalny); kontrola zgodności masy. Wyjście: --out/ev_fd.npz
(+ dowiązania ev_data/ev_mb). Progi WP (all_*) bez zmian — definicja WP jak w nominale. Użycie:
  relabel_pdgid.py --ev lxplus --hyp ../HYP_d0/runs/lxplus --out pdgid
"""
import argparse, ast, json, os
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--ev", required=True); ap.add_argument("--hyp", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
z = np.load(os.path.join(a.ev, "ev_fd.npz"), allow_pickle=False)
d = {k: z[k] for k in z.files}
meta = json.loads(str(z["meta"]))
groups = meta["groups"]; groups = ast.literal_eval(groups) if isinstance(groups, str) else groups
hm = json.load(open(os.path.join(a.hyp, "fd_nd.meta.json")))
fidx = {f: i for i, f in enumerate(hm["files"])}; off = np.asarray(hm["file_offsets"], np.int64)
gf = {g["group"]: g["files"] for g in groups}
fi = d["file_idx"]
pq_file = np.array([fidx[gf[x // 1000][x % 1000]] for x in fi], np.int64)
grow = off[pq_file] + d["row"]
P = pd.read_parquet(os.path.join(a.hyp, "fd_nd.parquet"), columns=["row", "m", "truth", "hypo_swap"]).set_index("row")
J = P.reindex(grow)
found = J["m"].notna().to_numpy()
dm = np.abs(J["m"].to_numpy() - d["m"].astype(float))
ok = found & (dm < 1e-5)
base = np.where(d["truth"] >= 100, d["truth"] - 100, d["truth"])
hs = np.where(ok, J["hypo_swap"].fillna(-1).to_numpy(), -1)
new = base.copy()
sw = np.isin(base, [0, 1]) & (hs == 1)
new[sw] += 100
stat = {"n": int(len(fi)), "found": int(found.sum()), "mass_ok": int(ok.sum()),
        "truth01_rows": int(np.isin(base, [0, 1]).sum()),
        "changed_sig_to_swap": int(((d["truth"] < 100) & (new >= 100)).sum()),
        "changed_swap_to_sig": int(((d["truth"] >= 100) & (new < 100)).sum()),
        "truth01_no_match": int((np.isin(base, [0, 1]) & ~ok).sum())}
print(json.dumps(stat, indent=1))
os.makedirs(a.out, exist_ok=True)
d["truth"] = new.astype(np.float32)
meta["truth_swap_source"] = "pdgid (hypo_swap v1.2 via HYP fd_nd.parquet)"; meta["relabel"] = stat
d["meta"] = np.array(json.dumps(meta))
np.savez(os.path.join(a.out, "ev_fd.npz"), **d)
for f in ("ev_data.npz", "ev_mb.npz"):
    t = os.path.join(a.out, f)
    if not os.path.exists(t):
        os.symlink(os.path.abspath(os.path.join(a.ev, f)), t)
json.dump(stat, open(os.path.join(a.out, "relabel.json"), "w"), indent=1)

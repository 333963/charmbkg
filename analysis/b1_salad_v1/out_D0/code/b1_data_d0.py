"""B1-D0: produkcyjny bieg SALAD v1.0 na danych D0 (data_p0+p1 po wecie refleksji, SB+pasmo) × MC tła (mb_p0–p3 truth≥8) → artefakt wag.
Jak b1_data.py (D⁺). Użycie: python3 b1_data_d0.py [hgb|sk_mlp] [data_src: 0 = data_p0, 01 = p0+p1]"""
import os, sys, json, time, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.environ.get("CHARMTAG_SRC", "/work/b1d0/in/src")); warnings.filterwarnings("ignore")
from charmtag.bkggen import salad_v1 as SV
cfg = SV.load_cfg(); kind = sys.argv[1] if len(sys.argv) > 1 else cfg["kind"]; dsrc = sys.argv[2] if len(sys.argv) > 2 else "0"
load = lambda n: (lambda z: pd.DataFrame({k: z[k] for k in z.files}))(np.load(f"sub/{n}.npz"))
D, B = load("data"), load("mcbkg"); B["row"] = np.arange(len(B))
D = D[D["src"].isin([int(c) for c in dsrc])].reset_index(drop=True)
B = B[B["src"].isin([0, 1])].reset_index(drop=True)        # jak S4: mb_p0+p1
sp = SV.species_params(cfg, "D0"); rD = SV.regions(D["m"], sp)
t0 = time.time()
p = SV.run(cfg, "D0", D[rD != "SR"], B, "out", sample=f"data_p{dsrc}_{kind}", kind=kind)
w, meta = SV.load_weights(p)
print(kind, p, os.path.getsize(p) // 1024, "kB", f"{time.time()-t0:.0f}s")
print(json.dumps({k: v for k, v in meta.items() if k not in ("feats", "cfg")}, indent=1))

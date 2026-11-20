"""B1: produkcyjny bieg SALAD v1.0 na danych (data_p0 D⁺ po wecie D*, SB+pasmo) × MC tła (mb_p0+p1 truth≥8) → artefakt wag."""
import os, sys, json, time, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.environ.get("CHARMTAG_SRC", "/work/b1/in/AA_merge/src")); warnings.filterwarnings("ignore")
from charmtag.bkggen import salad_v1 as SV
cfg = SV.load_cfg(); kind = sys.argv[1] if len(sys.argv) > 1 else cfg["kind"]
load = lambda n: (lambda z: pd.DataFrame({k: z[k] for k in z.files}))(np.load(f"sub/{n}.npz"))
D, B = load("data"), load("mcbkg"); B["row"] = np.arange(len(B))
sp = SV.species_params(cfg, "Dplus"); rD = SV.regions(D["m"], sp)
t0 = time.time()
p = SV.run(cfg, "Dplus", D[rD != "SR"], B, "out", sample=f"data_p0_{kind}", kind=kind)
w, meta = SV.load_weights(p)
print(kind, p, os.path.getsize(p) // 1024, "kB", f"{time.time()-t0:.0f}s")
print(json.dumps({k: v for k, v in meta.items() if k not in ("feats", "cfg")}, indent=1))

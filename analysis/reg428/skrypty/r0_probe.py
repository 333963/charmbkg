#!/usr/bin/env python3
"""REG428 r0 — sonda wejść (tylko odczyt): klucze tabel dn/sc/app, gałęzie trygera w ROOT v2.0, zgodność event_uid."""
import glob, json, sys, zipfile
import numpy as np
E = "/eos/user/k/klyczak"
def keys(f):
    with zipfile.ZipFile(f) as z:
        return [n[:-4] for n in z.namelist()]
k = keys(f"{E}/charm_q/dnet/data/dn_data_428770.npz")
print("dn keys D0 p_:", [x for x in k if x.startswith("D0__p_")])
print("dn keys other prefixes:", sorted({x.split("__")[0] + "__" + x.split("__")[1][:2] for x in k if "__" in x}))
for d in ("dnet", "dnet4"):
    f = f"{E}/charm_q/{d}/sc/D0.xgb.H.sb.n.s1.npz"
    try: print(d, "sc keys:", keys(f))
    except Exception as e: print(d, "sc:", e)
f = sorted(glob.glob(f"{E}/charm_q/dnet4/app/d150/app_data_428770_*.npz"))[0]
z = np.load(f); print(f, z.files); print(json.loads(str(z["meta"])))
u = z["D0__uid"]; print("uid task", np.unique(u >> np.uint64(36)), "fnum", np.unique((u >> np.uint64(22)) & np.uint64(0x3FFF))[:10], "n", len(u))
import uproot
rf = sorted(glob.glob(f"{E}/charm_raw_v20/data_428770/*.root"))
print("root files", len(rf))
t = uproot.open(rf[0])
print(t.keys()[:10])
tr = t[[k for k in t.keys() if "Charm" in k or "tree" in k.lower() or "reco" in k.lower()][0]] if len(t.keys()) > 1 else t[t.keys()[0]]
print(tr.name, tr.num_entries)
bk = tr.keys()
print([b for b in bk if "trigPassed" in b and ("AFP" in b or "J12" in b or "j20" in b)])
print([b for b in bk if b.startswith("EventInfo_") and "trig" not in b][:60])

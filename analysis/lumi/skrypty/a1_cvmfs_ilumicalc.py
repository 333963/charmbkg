"""a1: odczyt oficjalnych plików ilumicalc (cvmfs GoodRunsLists/data22_13p6TeV/2023{1005,1024}) -> świetlność per LB (csv w tabele/)."""
import sys, uproot, numpy as np, pandas as pd
out = sys.argv[1]
for f in sys.argv[2:]:
    r = uproot.open(f); print("\n#####", f.split("surowe/")[-1])
    t = r["LumiMetaData"]
    print(" branches:", [f"{k}:{b.typename}" for k, b in t.items()])
    num = [k for k, b in t.items() if "string" not in b.typename.lower() and "char" not in b.typename.lower()]
    a = pd.DataFrame({k: t[k].array(library="np") for k in num})
    for k in t.keys():
        if k not in num:
            try: print("  ", k, "unikalne:", sorted(set(map(str, t[k].array(library="np"))))[:5])
            except Exception as e: print("  ", k, "nieczytelne", type(e).__name__)
    tag = f.split("/")[-2] + "__" + f.split("/")[-1].replace(".root", "")
    a.to_csv(f"{out}/{tag}__LumiMetaData.csv", index=False)
    print(a.head(3).T.to_string())
    lc = [c for c in a.columns if "lumi" in c.lower()]
    for run, g in a.groupby("RunNbr"):
        print("  run", run, "n_LB", len(g), "LB", int(g.LBStart.min()), int(g.LBStart.max()), {c: float(g[c].sum()) for c in lc + ["LiveTime"] if c in g})

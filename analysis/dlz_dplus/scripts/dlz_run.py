"""DLZ_dplus: uruchomienie skryptów FEAT_charmnet (bez ich edycji) z nowymi wariantami D+ i wyjściem do DLZ_dplus/out.
Nowe warianty (GROUPS['Dplus']):
  ALL_noiso_nopid = dalitz + ptratio + maxdR + absimp   (ALL_noiso bez dE/dx)
  dalitz_ptratio  = dalitz + ptratio
Użycie: python dlz_run.py study  --test mc|data --variants ...
        python dlz_run.py sculpt mc|data
        python dlz_run.py transfer
        python dlz_run.py modelling
"""
import json
import os
import sys

FEAT = "/home/k/Documents/charm/charm_work/FEAT_charmnet/scripts"
DLZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, FEAT)
import study  # noqa: E402

study.BASE = DLZ  # frames -> symlink do FEAT_charmnet/frames; out -> DLZ_dplus/out
g = study.GROUPS["Dplus"]
g["ALL_noiso_nopid"] = g["dalitz"] + g["ptratio"] + g["maxdR"] + g["absimp"]
g["dalitz_ptratio"] = g["dalitz"] + g["ptratio"]

import numpy as np  # noqa: E402

cmd = sys.argv[1]
if cmd == "study":
    sys.argv = ["study.py", "--species", "Dplus", "--cap", "300000", "--nboot", "200"] + sys.argv[2:]
    study.main()
elif cmd in ("sculpt", "transfer"):
    import sculpt
    import transfer
    sculpt.BASE = transfer.BASE = DLZ
    if cmd == "sculpt":
        sculpt.main("Dplus", sys.argv[2])
    else:
        transfer.main("Dplus")
elif cmd == "modelling":
    import modelling as M
    from charmbkg.samples import reweight
    sp = "Dplus"
    S, Bmc, Bda, dall, mb = study.load(sp)
    sel = mb["ok"].to_numpy() & (mb["dm0"].to_numpy() > 0.110)
    MM = mb[sel].reset_index(drop=True)
    rng = np.random.default_rng(3)  # ten sam podzbiór co modelling.py
    n = 300000
    D = Bda.iloc[np.sort(rng.choice(len(Bda), min(n, len(Bda)), replace=False))].reset_index(drop=True)
    MM = MM.iloc[np.sort(rng.choice(len(MM), min(n, len(MM)), replace=False))].reset_index(drop=True)
    D["weight"] = 1.0
    wm = reweight(MM, D, study.cfg["training"]["reweight_bins"])
    res = {"n_data": len(D), "n_mc": len(MM), "c2st": {}}
    H = study.H
    sets = [("H", H), ("H+dalitz", H + g["dalitz"]), ("H+dalitz_ptratio", H + g["dalitz_ptratio"]),
            ("H+ALL_noiso_nopid", H + g["ALL_noiso_nopid"])]
    for seed in (0, 1, 2):
        for nm, feats in sets:
            v = M.c2st(np.nan_to_num(D[feats].to_numpy(float), nan=-9), np.nan_to_num(MM[feats].to_numpy(float), nan=-9), wm, seed=seed)
            res["c2st"].setdefault(nm, []).append(v)
            print(nm, seed, round(v, 4), flush=True)
    json.dump(res, open(os.path.join(DLZ, "out", "modelling_Dplus.json"), "w"), indent=1)

"""DLZ_d0: skrypty FEAT_charmnet (bez edycji) z nowymi grupami cech D0, wyjście do DLZ_d0/out, ramki DLZ_d0/frames (d0_feats.py).
Użycie: python d0_run.py study --test mc|data --variants ...   |  sculpt mc|data  |  transfer  |  modelling
XGB wymuszone na 3 wątki (limit rdzeni). Score H, H_seed2, H+ALL_noiso, H+ptratio skopiowane z FEAT jako cache (te same wiersze i repliki).
"""
import json
import os
import shutil
import sys

FEAT = os.environ.get("FEAT_DIR", "/home/k/Documents/charm/charm_work/FEAT_charmnet")
D0 = os.environ.get("D0_BASE") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.environ.get("FEAT_SCRIPTS", os.path.join(FEAT, "scripts")))
import xgboost as xgb  # noqa: E402

_XC = xgb.XGBClassifier


MODE = os.environ.get("D0_MODE", "full")      # screen: lżejszy XGB (lr 0,15; early stopping 30), osobny katalog DLZ_d0/screen, bez cache FEAT
NJ = int(os.environ.get("D0_NJ", "3"))


def _XGB3(**kw):
    kw["n_jobs"] = NJ
    return _XC(**kw)


xgb.XGBClassifier = _XGB3
import study  # noqa: E402

BASE = D0 if MODE == "full" else os.path.join(D0, "screen")
study.BASE = BASE
if MODE == "screen":
    os.makedirs(os.path.join(BASE, "out"), exist_ok=True)
    if not os.path.exists(os.path.join(BASE, "frames")):
        os.symlink(os.path.join(D0, "frames"), os.path.join(BASE, "frames"))
    study.cfg["training"]["xgb"].update(learning_rate=0.15, n_estimators=600, early_stopping_rounds=30)
study.xgb.XGBClassifier = _XGB3
g = study.GROUPS["D0"]
KIN2 = ["ctsK_nf", "ptAsym_nf"]
KINPT = ["dMinPt_nf", "dMaxPt_nf", "dMinPtRatio_cn", "dMaxPtRatio_cn"]
IP = ["d0prod_nf", "d0sigK_nf", "d0sigPi_nf", "d0sigMin_nf", "d0sigMax_nf", "d0chi2_nf", "dz0_nf", "absImpact_cn", "d0AbsZ0_cn", "d1AbsZ0_cn"]
TIME = ["ctau_nf"]
PID = ["bbK_K_nf", "bbK_pi_nf", "bbPi_pi_nf", "bbPi_K_nf", "pPi_nf", "nptsPi_nf", "dedxPi_nf"]
NEW = {"d0prod": ["d0prod_nf"], "ip": IP, "time": TIME, "kin2": KIN2, "kinall": KIN2 + KINPT + ["dR_nf"], "pidbb": PID,
       "NOKIN": IP + TIME, "NOKIN_pid": IP + TIME + PID, "ALLnew": IP + TIME + PID + KIN2 + KINPT + ["dR_nf"],
       "pidshift_ctrl": ["dedxPiSh_nf"]}   # kontrola artefaktu przesunięcia dE/dx MC (tylko modelling)
g.update(NEW)
import numpy as np  # noqa: E402

for t in (("mc", "data") if MODE == "full" else ()):
    src, dst = os.path.join(FEAT, "out", f"D0_{t}"), os.path.join(D0, "out", f"D0_{t}")
    os.makedirs(dst, exist_ok=True)
    for v in ("H", "H_seed2", "H+ALL_noiso", "H+ptratio"):
        for ext in ("npz", "json"):
            a, b = os.path.join(src, f"scores_{v}.{ext}"), os.path.join(dst, f"scores_{v}.{ext}")
            if os.path.exists(a) and not os.path.exists(b):
                shutil.copy2(a, b)

cmd = sys.argv[1] if len(sys.argv) > 1 else "none"
if cmd == "study":
    sys.argv = ["study.py", "--species", "D0", "--cap", "300000" if MODE == "full" else "150000", "--nboot", "200"] + sys.argv[2:]
    study.main()
elif cmd in ("sculpt", "transfer"):
    import sculpt
    import transfer
    sculpt.BASE = transfer.BASE = BASE
    if not os.path.exists(os.path.join(BASE, "out", "dstar_iso_fit.json")):
        shutil.copy2(os.environ.get("DSTAR_ISO_FIT") or os.path.join(FEAT, "out", "dstar_iso_fit.json"), os.path.join(BASE, "out", "dstar_iso_fit.json"))
    if cmd == "sculpt":
        sculpt.main("D0", sys.argv[2])
    else:
        transfer.main("D0")
elif cmd == "modelling":
    import modelling as M
    from charmbkg.samples import reweight
    sp = "D0"
    m0 = study.cfg["species"][sp]["m0"]
    S, Bmc, Bda, dall, mb = study.load(sp)
    sel = mb["ok"].to_numpy() & (mb["dm0"].to_numpy() > 0.100) & (np.abs(mb["m_swap"].to_numpy() - m0) > 0.025)
    MM = mb[sel].reset_index(drop=True)
    rng = np.random.default_rng(3)
    n = 300000
    D = Bda.iloc[np.sort(rng.choice(len(Bda), min(n, len(Bda)), replace=False))].reset_index(drop=True)
    MM = MM.iloc[np.sort(rng.choice(len(MM), min(n, len(MM)), replace=False))].reset_index(drop=True)
    D["weight"] = 1.0
    wm = reweight(MM, D, study.cfg["training"]["reweight_bins"])
    res = {"n_data": len(D), "n_mc": len(MM), "c2st": {}, "c2st1d": {}}
    H = study.H
    names = sys.argv[2].split(",") if len(sys.argv) > 2 else list(NEW)
    sets = [("H", H)] + [("H+" + k, H + NEW[k]) for k in names]
    for seed in (0, 1, 2):
        for nm, feats in sets:
            v = M.c2st(np.nan_to_num(D[feats].to_numpy(float), nan=-9), np.nan_to_num(MM[feats].to_numpy(float), nan=-9), wm, seed=seed)
            res["c2st"].setdefault(nm, []).append(v)
            print(nm, seed, round(v, 4), flush=True)
    for f in sorted(set(sum(NEW.values(), []))):
        res["c2st1d"][f] = M.c2st(np.nan_to_num(D[[f]].to_numpy(float), nan=-9), np.nan_to_num(MM[[f]].to_numpy(float), nan=-9), wm)
        res["c2st1d"][f + "|ks"] = M.wks(np.nan_to_num(D[f].to_numpy(float), nan=-9), np.ones(len(D)), np.nan_to_num(MM[f].to_numpy(float), nan=-9), wm)
        print(f, round(res["c2st1d"][f], 4), round(res["c2st1d"][f + "|ks"], 4), flush=True)
    json.dump(res, open(os.path.join(BASE, "out", "modelling_D0.json"), "w"), indent=1)

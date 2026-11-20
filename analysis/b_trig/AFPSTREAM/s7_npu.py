#!/usr/bin/env python3
"""Sprawdzenie: klasy AFP z wetem nPileupVertices == 0 (weto stosowane w głównym łańcuchu analizy) — tło f i SF sygnału (SBS).
Wynik: out/s7.json."""
import numpy as np
import afp_common as A
from iso_common import load_pairs, sbs_w, eff_w, C
import s2_universal_rw as U
from s1_iso import samples, bkg_sel
M = load_pairs("bt_fd_nd"); mm = U.mc_sel(M, 0.040, 0.6); emc = float((M["ptcone40"][mm] < 2).mean())
B = load_pairs("bt_mb_nd"); bm = bkg_sel(B, np.ones(len(B["m"]), bool), 3, 30)
D8, D9, S = samples()
res = {}
for name in ["8:AFPstream", "8:AFPsptrk4", "8:AFPnoalg_any", "8:AFP&J12", "8:J12only!A", "8:MBTS!A", "8:none", "9:MB427929"]:
    D, cm = S[name]
    for lab, sel in (("npu=0", D["e_npu"] == 0), ("npu>=1", D["e_npu"] >= 1)):
        c = cm & sel; pas = D["ptcone40"] < 2
        bb = bkg_sel(D, c, 3, 30); md = U.base(D, 0.040, 0.6) & c
        r = {"frac_pairs": float(c.sum() / max(cm.sum(), 1)), "n_bkg": int(bb.sum())}
        if bb.sum() > 30:
            f = float(pas[bb].mean()); r["f"] = f; r["f_err"] = float(np.sqrt(f * (1 - f) / bb.sum()))
            eps_c, _, _ = U.mc_eps_for(D, bb, B, bm, B["ptcone40"] < 2, U.SCHEMES["pT×lead_pt×nSelBL"]); r["R_rw"] = f / float(eps_c.mean())
        if md.sum() > 20:
            w = sbs_w(D["dm"][md], D["pt"][md], pas[md]); e, er, T, eT = eff_w(w, pas[md].astype(float))
            r.update({"S": T, "S_err": eT, "SF": e / emc, "SF_err": er / emc})
            rw = U.sf_rw(D, md, pas, M, mm, M["ptcone40"] < 2, U.SCHEMES["pT×lead_pt×nSelBL"])
            if rw:
                r["SF_rw"] = rw["SF"]; r["SF_rw_err"] = rw["err"]
        res[f"{name}|{lab}"] = r
        print(name, lab, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
A.jdump(res, f"{A.OUT}/s7.json")

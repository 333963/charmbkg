#!/usr/bin/env python3
"""TAGGED — kontrolka K2: D*+ → D0 π (kandydaci slim afp_ext, bez klasyfikatora), fit Δm. Tryby: prep (bez okna sygnału w zdarzeniach z tagiem) | unblind.
Wyjście: out/dst_prep.json, out/dst_unblind.json (+ histogramy do rysunków w out/dst_hist.npz)."""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tg_lib as L  # noqa: E402

AM = os.path.join(L.W, "AFP_mix", "runs", "lxplus")
MPI = 0.13957
LO, HI, BWD = 0.13945, 0.15995, 0.0005     # Δm w afp_ext skwantowane co 0,1 MeV → bin 0,5 MeV z krawędziami między wartościami (odstępstwo od 0,25 MeV, przed odczytem)
WLO, WHI = 0.14345, 0.14745
E = np.arange(LO, HI + BWD / 2, BWD)
CTR = 0.5 * (E[1:] + E[:-1])
X = np.maximum(CTR - MPI, 1e-5)
WM = (CTR > WLO) & (CTR < WHI)
ITEMS = {"Dst|427929|MBTS": ("427929", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_mbts_L1MBTS_1"), "Dst|428770|R1": ("428770", "HLT_noalg_L1AFP_A_OR_C", None)}


def load(run, trig_tag, trig_incl, blind):
    import pyarrow.parquet as pq
    cols = ["event_number", "passGRL_MinBiasAFP", "nPileupVertices", "xi_cal_A", "xi_cal_C", "afp_xi1_A", "afp_xi1_C", "afp_tag_A", "afp_tag_C"]
    tc = sorted({t for t in (trig_tag, trig_incl) if t})
    ev = pq.read_table(os.path.join(AM, f"data_{run}_events.parquet"), columns=cols + ["trigPassed_" + t for t in tc]).to_pandas()
    ev = ev[(ev.passGRL_MinBiasAFP > 0) & (ev.nPileupVertices == 0)]
    a, c = ev.afp_tag_A.to_numpy() > 0, ev.afp_tag_C.to_numpy() > 0
    xok = (~a | (ev.xi_cal_A.to_numpy() <= ev.afp_xi1_A.to_numpy() + 0.02)) & (~c | (ev.xi_cal_C.to_numpy() <= ev.afp_xi1_C.to_numpy() + 0.02))
    t2 = (a | c) & xok & (ev["trigPassed_" + trig_tag].to_numpy() > 0)
    inc = (ev["trigPassed_" + trig_incl].to_numpy() > 0) if trig_incl else np.ones(len(ev), bool)
    en = ev.event_number.to_numpy()
    o = np.argsort(en); en, t2, inc = en[o], t2[o], inc[o]
    ca = pq.read_table(os.path.join(AM, f"data_{run}_cands.parquet"), columns=["event_number", "Dmeson_species", "Dmeson_delta_mass", "Dmeson_pt"]).to_pandas()
    ca = ca[(ca.Dmeson_species == 1) & (ca.Dmeson_pt >= 3) & (ca.Dmeson_pt < 30)]
    ix = np.clip(np.searchsorted(en, ca.event_number.to_numpy()), 0, len(en) - 1)
    ok = en[ix] == ca.event_number.to_numpy()
    dm = ca.Dmeson_delta_mass.to_numpy()
    mi, mt = ok & inc[ix], ok & t2[ix]
    if blind:
        mt = mt & ~((dm > WLO) & (dm < WHI))
    return np.histogram(dm[mi], E)[0].astype(float), np.histogram(dm[mt], E)[0].astype(float)


def gcdf(x, mu, s):
    from scipy.special import erf
    return 0.5 * (1 + erf((x - mu) / (s * math.sqrt(2))))


def sshape(mu, s1, s2, fr):
    return fr * np.diff(gcdf(E, mu, s1)) + (1 - fr) * np.diff(gcdf(E, mu, s2))


def incl_fit(H):
    from iminuit import Minuit

    def model(p):
        b = np.power(X, p[6]) * np.exp(np.clip(p[7] * X * 50 + p[8] * (X * 50) ** 2, -50, 50))
        return np.maximum(p[0] * sshape(p[2], p[3], p[4], p[5]) + p[1] * b / b.sum(), 1e-9)

    def cost(p):
        mu = model(p)
        return 2.0 * float(np.sum(mu - H) + np.sum(np.where(H > 0, H * np.log(np.where(H > 0, H, 1.0) / mu), 0.0)))
    m = Minuit(cost, [0.002 * H.sum(), H.sum(), 0.14543, 0.0004, 0.0010, 0.6, 0.5, 0.0, 0.0], name=["NS", "NB", "mu", "s1", "s2", "fr", "a", "b", "c"])
    m.errordef = 1.0
    m.limits["NS"] = (0, H.sum()); m.limits["NB"] = (0, 2 * H.sum()); m.limits["mu"] = (0.1450, 0.1458); m.limits["s1"] = (0.0002, 0.0008); m.limits["s2"] = (0.0008, 0.0025)
    m.limits["fr"] = (0.2, 1.0); m.limits["a"] = (0.0, 3.0); m.limits["b"] = (-20, 20); m.limits["c"] = (-20, 20)
    m.migrad(ncall=20000); m.hesse()
    v = np.array(m.values)
    mu = model(v)
    return {"N": float(v[0]), "err": float(m.errors["NS"]), "par": [float(t) for t in v], "valid": bool(m.valid), "chi2": float(np.sum((H - mu) ** 2 / mu)), "ndf": len(H) - 9}, mu


def main(mode):
    from multiprocessing import Pool
    out, hs = {}, {}
    prep = json.load(open(os.path.join(L.OUT, "dst_prep.json"))) if mode == "unblind" else {}
    for key, (run, tt, ti) in ITEMS.items():
        Hi, Ht = load(run, tt, ti, blind=(mode == "prep"))
        fi, mui = incl_fit(Hi)
        S = sshape(*fi["par"][2:6])
        r = {"N_incl": fi["N"], "N_incl_err": fi["err"], "incl_valid": fi["valid"], "incl_chi2": fi["chi2"], "incl_ndf": fi["ndf"], "incl_par": fi["par"],
             "nsb_tag": float(Ht[~WM].sum()), "nsb_incl": float(Hi[~WM].sum())}
        r["P_tag_SB"] = r["nsb_tag"] / r["nsb_incl"]
        if mode == "prep":
            t = L.TFit(S, len(Ht), 1, mask=~WM, bkind="pow", x=X)
            vb = t.fit(Ht, fixNS=0.0)[0]
            r["bkg_sb"] = [float(q) for q in vb]
            npred = max(round(fi["N"] * r["P_tag_SB"]), 1)
            r["toys"] = []
            for p in (0, npred):
                with Pool(7) as pool:
                    T = np.concatenate(pool.map(L.toy_worker, [(S, len(Ht), 1, vb, ~WM, float(p), 43, 1000 * (c + 1) + p + (7 if run == "428770" else 0), "pow", X) for c in range(7)]))
                r["toys"].append(L.toy_summary(T, p))
        else:
            t = L.TFit(S, len(Ht), 1, bkind="pow", x=X)
            r["fit"] = t.full(Ht)
            r["n_window"] = float(Ht[WM].sum())
            hs[key.replace("|", "__") + "__tag"] = Ht; hs[key.replace("|", "__") + "__mu"] = t.mu(np.array(r["fit"]["x"])); hs[key.replace("|", "__") + "__bkg"] = t.mu(np.array([0.0] + r["fit"]["x"][1:]))
            hs[key.replace("|", "__") + "__incl"] = Hi; hs[key.replace("|", "__") + "__inclmu"] = mui
            r["toys"] = prep.get(key, {}).get("toys")
        out[key] = r
        print(key, {k: v for k, v in r.items() if k not in ("incl_par", "bkg_sb")}, flush=True)
    json.dump(out, open(os.path.join(L.OUT, f"dst_{mode}.json"), "w"), indent=1)
    if hs:
        np.savez_compressed(os.path.join(L.OUT, "dst_hist.npz"), edges=E, **hs)


if __name__ == "__main__":
    if sys.argv[1] == "unblind" and not os.path.exists(os.path.join(L.OUT, "UNBLINDED")):
        sys.exit("STOP: uruchamiane tylko z t2_unblind.py")
    main(sys.argv[1])

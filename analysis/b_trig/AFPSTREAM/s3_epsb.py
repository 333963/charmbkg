#!/usr/bin/env python3
"""Krok 3: ε_B pod pikiem (exp(pol2) z SB, jak BKG_427929/bkg_common) per klasa trygera, dane/MC (truth ≥ 2 i ≥ 8, zliczenie w SR),
oraz stosunek dane/dane w samych pasmach bocznych (zliczenie) klasa/427929 i klasa/J12only. Bootstrap Poissona 300 replik.
Wynik: out/s3.json."""
import numpy as np
import afp_common as A
import bkg_common as K

NREP = 300
rng = np.random.default_rng(428770)
WPI = [0, 2, 4]                                   # WP95, WP80, WP50
thr = K.thresholds()
z8 = np.load(f"{A.R}/bt_data_428770.npz", allow_pickle=True)
z9 = np.load(f"{A.R}/bt_data_427929.npz", allow_pickle=True)
zm = np.load(f"{A.R}/bt_mb_nd.npz", allow_pickle=True)
c8 = A.classes8(A.bits(z8, "all_e_")); c8["AFPnoalg_any"] = c8["_an"]
c9 = A.classes9(A.bits(z9, "all_e_"), z9["all_e_grl"] == 1)
SETS = {}
for k in A.DISJ8 + ["AFPonly", "AFPsptrk4_any", "AFPnoalg_any", "AFPstream", "J12only", "MBTSbased", "other", "all"]:
    SETS["8:" + k] = (z8, c8[k])
for k in ["MB427929", "MB9&!j20&!J12", "MB9&j20&!J12", "MB9&AFPnoalg"]:
    SETS["9:" + k] = (z9, c9[k])
tm = zm["all_truth"]
MC = {"mc_t2": tm >= 2, "mc_t8": tm >= 8, "mc_all": np.ones(len(tm), bool)}
res = {"thr": thr, "rows": {}, "mc": {}}
H = {}
for name, (z, cm) in list(SETS.items()) + [(k, (zm, v)) for k, v in MC.items()]:
    for b, (lo, hi) in K.BINS.items():
        s = cm & (z["all_pt"] >= lo) & (z["all_pt"] < hi)
        H[(name, b)] = K.cat_hist(z["all_m"][s], z["all_score"][s], thr[b])


def sb_frac(h):
    ha = h.sum(0)
    return np.array([h[j + 1:, K.SBm].sum() / max(ha[K.SBm].sum(), 1e-12) for j in range(5)])


def rep(h, method):
    out = np.zeros((NREP, 5))
    for i in range(NREP):
        out[i] = K.all_eps(K.boot(h, rng), method) if method != "sbfrac" else sb_frac(K.boot(h, rng))
    return out


for b in K.BINS:
    for mk in MC:
        h = H[(mk, b)]
        e = K.all_eps(h, "true"); rp = rep(h, "true")
        res["mc"][f"{mk}|{b}"] = {"eps_sr": e, "err": rp.std(0), "sbfrac": sb_frac(h), "sbfrac_err": rep(h, "sbfrac").std(0), "n_sr": int(h.sum(0)[K.SRm].sum())}
for name in SETS:
    for b in K.BINS:
        h = H[(name, b)]
        ha = h.sum(0)
        r = {"n": int(ha.sum()), "n_sb": int(ha[K.SBm].sum()), "n_sb_pass": [int(h[j + 1:, K.SBm].sum()) for j in range(5)]}
        if r["n_sb"] >= 200:
            e = K.all_eps(h, "pol2"); rp = rep(h, "pol2")
            r["pol2"] = e; r["pol2_err"] = np.nanstd(rp, 0)
            r["sym"] = K.all_eps(h, "sym"); r["lr"] = K.all_eps(h, "lr"); r["pol1"] = K.all_eps(h, "pol1")
            r["sbfrac"] = sb_frac(h); r["sbfrac_err"] = rep(h, "sbfrac").std(0)
            for mk in ("mc_t2", "mc_t8"):
                em = res["mc"][f"{mk}|{b}"]
                r["R_" + mk] = e / em["eps_sr"]
                r["R_" + mk + "_stat"] = r["R_" + mk] * np.hypot(r["pol2_err"] / e, em["err"] / em["eps_sr"])
            r["syst_rel"] = np.max([np.abs(r[v] / e - 1) for v in ("sym", "lr", "pol1")], axis=0)
            # like-for-like w SB: frakcja pass w SB dane / MC inkluzywne
            mm = res["mc"][f"mc_all|{b}"]
            r["Rsb_mc"] = r["sbfrac"] / mm["sbfrac"]
            r["Rsb_mc_err"] = r["Rsb_mc"] * np.hypot(r["sbfrac_err"] / np.maximum(r["sbfrac"], 1e-12), mm["sbfrac_err"] / mm["sbfrac"])
        res["rows"][f"{name}|{b}"] = r
A.jdump(res, f"{A.OUT}/s3.json")
for b in K.BINS:
    print("==", b, "MC t2", np.round(res["mc"][f"mc_t2|{b}"]["eps_sr"][WPI], 4))
    for name in SETS:
        r = res["rows"][f"{name}|{b}"]
        if "pol2" not in r:
            print(f"{name:20s} n_sb {r['n_sb']}"); continue
        print(f"{name:20s} n_sb {r['n_sb']:8d} R_t2 " + " ".join(f"{r['R_mc_t2'][j]:.3f}±{r['R_mc_t2_stat'][j]:.3f}(s{r['syst_rel'][j]:.2f})" for j in WPI)
              + " | Rsb " + " ".join(f"{r['Rsb_mc'][j]:.3f}±{r['Rsb_mc_err'][j]:.3f}" for j in WPI))

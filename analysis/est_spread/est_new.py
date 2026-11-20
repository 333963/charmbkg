"""EST_spread — NOWE liczby (17.11.2026): zestaw estymatorów ε_B pod pikiem D0 na tych samych próbkach bt_*.npz
(B_trig/runs/lxplus, tylko odczyt; inkluzywne D0 `all_*`, 30 % zdarzeń, v2.0), per run × bin pT × WP.

Estymatory (definicje i kod: B_trig/BKG_427929/bkg_common.py, importowane bez zmian):
  pasma boczne danych: pol2 = ML exp(pol2) [nominał], pol1 = ML exp(pol1), pol3 = ML exp(pol3), sym = zliczenie |dm| w [75,125] MeV,
                       lr = liniowa interpolacja ε(lewe), ε(prawe) do m0, far = ML exp(pol2) tylko na |dm| >= 90 MeV (dalekie pasma)
  MC: mcraw_t2 / mcraw_t8 = zliczenie w SR (truth >= 2 / >= 8), mcV1_t2 / mcV2_t2 = MC truth >= 2 w SR ważone wagami
      z B_trig/BKG_427929/out/s3_w_<run>.npz (V1: bez score; V2: ze score = przeważenie typu SALAD, uczone w pasmach bocznych)
Domknięcie na MC: te same estymatory pasm bocznych na tle MC (truth >= 2 oraz >= 8) wobec zliczenia tła w SR (prawda).
Błędy: bootstrap Poissona (dane: histogram kategoria×masa; MC: wagi Poissona per kandydat), NREP replik; pełna kowariancja między estymatorami.
Wyjście: out/est_new.json, out/est_new_rep.npz
"""
import json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "B_trig", "BKG_427929"))
import bkg_common as C  # noqa: E402

OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
NREP = 300
WPS = C.WPS
SB_EST = ["pol2", "pol1", "pol3", "sym", "lr", "far"]
FARm = np.abs(C.CEN - C.M0) >= 0.090


def eps_far(ha, hp):
    a, b = C.pred_sr(C.glm_pois(ha, 2, FARm)), C.pred_sr(C.glm_pois(hp, 2, FARm))
    return b / a if a > 0 else np.nan


def sb_all(H):
    """→ tablica [estymator, WP]"""
    ha = H.sum(0)
    out = []
    for e in SB_EST:
        f = eps_far if e == "far" else C.METHODS[e]
        out.append([f(ha, H[j + 1:].sum(0)) for j in range(H.shape[0] - 1)])
    return np.array(out)


def true_sr(H):
    ha = H.sum(0)
    return np.array([C.eps_true_sr(ha, H[j + 1:].sum(0)) for j in range(H.shape[0] - 1)])


t0 = time.time()
thr = C.thresholds()
rng = np.random.default_rng(20261117)
res = {"meta": {"nrep": NREP, "wps": WPS, "sb_est": SB_EST, "thr": thr,
                "regions": {"SR_half": C.SR_HALF, "SB_min": C.SB_MIN, "SB_max": C.SB_MAX, "window": [C.LO, C.HI], "far_min": 0.090}},
       "data": {}, "mc_closure": {}, "mc_est": {}, "N": {}}
reps = {}

# ---------- dane
for run in ("d427", "d428"):
    D = C.load_all(run)
    for b, (lo, hi) in C.BINS.items():
        s = (D["pt"] >= lo) & (D["pt"] < hi)
        H = C.cat_hist(D["m"][s], D["score"][s], thr[b])
        v = sb_all(H)
        R = np.array([sb_all(C.boot(H, rng)) for _ in range(NREP)])          # [rep, est, wp]
        res["data"][f"{run}|{b}"] = {e: {"val": v[i].tolist(), "err": np.nanstd(R[:, i], 0).tolist()} for i, e in enumerate(SB_EST)}
        reps[f"{run}|{b}"] = R
        res["N"][f"{run}|{b}"] = {"n": int(s.sum()), "n_sb": int(H.sum(0)[C.SBm].sum()), "n_sr": int(H.sum(0)[C.SRm].sum()),
                                   "n_sb_pass": [int(H[j + 1:, C.SBm].sum()) for j in range(5)]}
        print(run, b, "pol2", np.round(v[0], 4), f"{time.time()-t0:.0f}s", flush=True)
    del D

# ---------- MC
mb = C.load_all("mb")
sr = np.abs(mb["m"] - C.M0) < C.SR_HALF
i_sr = np.flatnonzero(sr)
W = {r: np.load(os.path.join(HERE, "..", "B_trig", "BKG_427929", "out", f"s3_w_{r}.npz")) for r in ("d427", "d428")}
for r in W:
    assert len(W[r]["V2 + score|sr"]) == len(i_sr), (len(W[r]["V2 + score|sr"]), len(i_sr))
for b, (lo, hi) in C.BINS.items():
    T = np.array(thr[b])
    inb = (mb["pt"] >= lo) & (mb["pt"] < hi)
    # domknięcie interpolacji na tle MC
    for tn, tmin in (("t2", 2), ("t8", 8)):
        s = inb & (mb["truth"] >= tmin)
        H = C.cat_hist(mb["m"][s], mb["score"][s], thr[b])
        v, tv = sb_all(H), true_sr(H)
        Rr, Rt = [], []
        for _ in range(NREP):
            Hb = C.boot(H, rng)
            Rr.append(sb_all(Hb)); Rt.append(true_sr(Hb))
        Rr, Rt = np.array(Rr), np.array(Rt)
        ratio = Rr / Rt[:, None, :]
        res["mc_closure"][f"{tn}|{b}"] = {"true": {"val": tv.tolist(), "err": Rt.std(0).tolist()},
                                           **{e: {"val": v[i].tolist(), "err": np.nanstd(Rr[:, i], 0).tolist(),
                                                  "ratio": (v[i] / tv).tolist(), "ratio_err": np.nanstd(ratio[:, i], 0).tolist()}
                                              for i, e in enumerate(SB_EST)}}
        reps[f"mcclos_{tn}|{b}"] = ratio
        res["N"][f"mc_{tn}|{b}"] = {"n": int(s.sum()), "n_sb": int(H.sum(0)[C.SBm].sum()), "n_sr": int(H.sum(0)[C.SRm].sum()),
                                     "n_sb_pass": [int(H[j + 1:, C.SBm].sum()) for j in range(5)],
                                     "n_sr_pass": [int(H[j + 1:, C.SRm].sum()) for j in range(5)]}
        print("MC closure", tn, b, "pol2/true", np.round(v[0] / tv, 3), flush=True)
    # estymatory z MC w SR (surowe i ważone)
    pt_sr, tr_sr, sc_sr = mb["pt"][i_sr], mb["truth"][i_sr], mb["score"][i_sr]
    sb_ = (pt_sr >= lo) & (pt_sr < hi)
    s2, s8 = sb_ & (tr_sr >= 2), sb_ & (tr_sr >= 8)
    P2 = (sc_sr[s2][:, None] >= T[None, :]).astype(np.float32)
    is8 = (tr_sr[s2] >= 8)
    names = ["mcraw_t2", "mcraw_t8", "mcV1_t2_d427", "mcV2_t2_d427", "mcV1_t2_d428", "mcV2_t2_d428"]
    ws = [np.ones(s2.sum(), np.float32), is8.astype(np.float32),
          W["d427"]["V1 bez score|sr"][s2], W["d427"]["V2 + score|sr"][s2], W["d428"]["V1 bez score|sr"][s2], W["d428"]["V2 + score|sr"][s2]]
    val = np.array([(w[:, None] * P2).sum(0) / w.sum() for w in ws])
    R = np.zeros((NREP, len(ws), 5))
    for k in range(NREP):
        kk = rng.poisson(1.0, s2.sum()).astype(np.float32)
        for i, w in enumerate(ws):
            kw = kk * w
            R[k, i] = (kw @ P2) / kw.sum()
    res["mc_est"][b] = {n: {"val": val[i].tolist(), "err": R[:, i].std(0).tolist(),
                            "neff_over_n": float(ws[i].sum() ** 2 / (ws[i] ** 2).sum() / max((ws[i] > 0).sum(), 1))} for i, n in enumerate(names)}
    reps[f"mcest|{b}"] = R
    print("MC est", b, {n: np.round(val[i], 4).tolist() for i, n in enumerate(names)}, flush=True)

res["meta"]["cpu_s"] = time.time() - t0
json.dump(res, open(os.path.join(OUT, "est_new.json"), "w"), indent=1)
np.savez_compressed(os.path.join(OUT, "est_new_rep.npz"), **reps)
print("OK", f"{time.time()-t0:.0f}s")

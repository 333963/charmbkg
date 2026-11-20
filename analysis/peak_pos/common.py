"""PEAK_pos: wspólne — próbki, progi WP, biny."""
import numpy as np, pandas as pd
OUT = "/home/k/Documents/charm/charm_work/PEAK_pos/out"
M0 = {"D0": 1.86484, "Dplus": 1.86966}
THR = {  # progi WP modelu analizy (models/D0_H, models/Dplus_H meta.json)
    "D0": {"3-5": {0.8: 0.636467695236206, 0.7: 0.7913140654563904, 0.5: 0.9381963610649109},
           "5-8": {0.8: 0.696060061454773, 0.7: 0.8363182544708252, 0.5: 0.9542829394340515},
           "8-30": {0.8: 0.6751556992530823, 0.7: 0.818087100982666, 0.5: 0.9405738711357117}},
    "Dplus": {"3-5": {0.8: 0.84221351146698, 0.7: 0.9323874711990356, 0.5: 0.9836874008178711},
              "5-8": {0.8: 0.8692589402198792, 0.7: 0.9428622722625732, 0.5: 0.984621524810791},
              "8-30": {0.8: 0.8230705857276917, 0.7: 0.920244038105011, 0.5: 0.9780784249305725}}}
PTB = [(3, 5), (5, 8), (8, 30)]
VETO = 0.025  # D0: weto refleksji |m_swap − m0| > 25 MeV


def thr_own(score, pt, wp, w=None):
    """progi WP z własnego sygnału MC (dla MC FEAT z innym modelem)"""
    w = np.ones(len(score)) if w is None else w
    t = {}
    for lo, hi in PTB:
        a = (pt >= lo) & (pt < hi)
        s, ww = score[a], w[a]
        o = np.argsort(-s); c = np.cumsum(ww[o]) / ww.sum()
        t[f"{lo}-{hi}"] = {wp: s[o][min(np.searchsorted(c, wp), len(o) - 1)]}
    return t


def pass_wp(score, pt, thr, wp):
    ok = np.zeros(len(score), bool)
    for lo, hi in PTB:
        ok |= (pt >= lo) & (pt < hi) & (score >= thr[f"{lo}-{hi}"][wp])
    return ok


def load(sp, tag, wp=0.7, veto=True):
    """zwraca DataFrame po cięciu WP (i weto refleksji dla D0, gdy veto=True). tag: data | mc (D0: HYP fd; D+: FEAT fd) | mbmc (D+)"""
    if sp == "D0":
        L = pd.read_parquet(f"{OUT}/D0_{tag}30.parquet")
        thr = THR["D0"]
    elif tag == "mc":
        L = pd.read_parquet(f"{OUT}/Dplus_mc.parquet")      # FEAT fd (własny model H → własne progi)
        L["q"] = L["charge"].astype(np.int8)
        thr = thr_own(L.score.to_numpy(), L.pt.to_numpy(), wp)
    else:
        L = pd.read_parquet(f"{OUT}/Dplus_{tag}30.parquet")
        thr = THR["Dplus"]
    ok = pass_wp(L.score.to_numpy(), L.pt.to_numpy(), thr, wp)
    if sp == "D0" and veto:
        ok &= np.abs(L.m_swap.to_numpy() - M0["D0"]) > VETO
    return L[ok].reset_index(drop=True)

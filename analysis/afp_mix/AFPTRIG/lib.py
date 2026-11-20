#!/usr/bin/env python3
"""AFPTRIG/lib.py — metoda kształtu (styl CMS-TOTEM 2102.06945 / ALFA 2212.00664) dla tła przypadkowego tagu AFP
w strumieniu trygera AFP. Projekt B, WSTĘPNE. Nic tu nie liczy uzysków sygnału D (żadnych fitów masy).

Estymator (per strona s, komórka c = blok LB × klasa nSelBL):
    N_acc(SR) = Σ_c N_CR,c(dane) · R_c,   R_c = T_c(Δ ≤ δ) / T_c(Δ > δ + m),   Δ = ξ_cal,s − ξ_AFP,s
T_c = szablon z mieszania: ξ_cal,s zdarzenia centralnego (pula centralna, bez protonu po stronie s) minus ξ protonu
przypadkowego (pula protonów) z tej samej komórki LB. Liczone dokładnie (wszystkie pary) przez searchsorted.
"""
from __future__ import annotations
import json, os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
INDIR = os.path.join(os.path.dirname(HERE), "runs", "lxplus")
CACHE = os.path.join(HERE, "cache")
RUNS = (427929, 428770)
SIDES = ("A", "C")
XI_LO, XI_HI, DELTA = 0.030, 0.100, 0.02
M_NOM = 0.025
AFP_CHAINS = ("HLT_noalg_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C")
SPECIES = {"Dplus": dict(code=0, var="Dmeson_mass", peak=(1.84, 1.90), sb=((1.75, 1.80), (1.94, 2.00))),
           "D0": dict(code=3, var="Dmeson_mass", peak=(1.835, 1.895), sb=((1.75, 1.80), (1.93, 2.00))),
           "Dstar": dict(code=1, var="Dmeson_delta_mass", peak=(0.1440, 0.1470), sb=((0.150, 0.160),))}
PT_BINS = ((3.0, 5.0), (5.0, 8.0), (8.0, 30.0))
CATS = [f"{sp}|{lo:g}-{hi:g}|{w}" for sp in SPECIES for lo, hi in PT_BINS for w in ("peak", "sb")]
NSEL_EDGES = np.array([0, 5, 10, 15, 20, 25, 30, 40, 50, 80, 10000])


def load(run, force=False):
    """Zdarzenia po GRL ∧ npu == 0 z ≥ 1 zapisanym łańcuchem (AFP lub nie-AFP) → dict tablic; macierz liczby kandydatów per kategoria."""
    os.makedirs(CACHE, exist_ok=True)
    fn = os.path.join(CACHE, f"ev_{run}.npz")
    if os.path.exists(fn) and not force:
        z = np.load(fn, allow_pickle=False)
        return {k: z[k] for k in z.files}
    ev = pd.read_parquet(os.path.join(INDIR, f"data_{run}_events.parquet"))
    tc = [c for c in ev.columns if c.startswith("trigPassed_")]
    n0 = len(ev)
    P = {c[len("trigPassed_"):]: ev[c].to_numpy() == 1 for c in tc}
    afp = np.zeros(n0, bool); non = np.zeros(n0, bool)
    for k, v in P.items():
        if "L1AFP" in k: afp |= v
        else: non |= v
    grl = ev["passGRL_MinBiasAFP"].to_numpy() == 1
    npu0 = ev["nPileupVertices"].to_numpy() == 0
    info = {"n_events": n0, "n_grl": int(grl.sum()), "n_grl_npu0": int((grl & npu0).sum()), "n_none": int((grl & ~afp & ~non).sum()),
            "n_afp_any_grl": int((grl & afp).sum()), "n_afp_only_grl": int((grl & afp & ~non).sum()), "n_nonafp_grl": int((grl & non).sum()),
            "n_afp_any_grl_npu0": int((grl & afp & npu0).sum())}
    ca = pd.read_parquet(os.path.join(INDIR, f"data_{run}_cands.parquet"), columns=["event_number", "Dmeson_species", "Dmeson_mass", "Dmeson_delta_mass", "Dmeson_pt"])
    crow = pd.Index(ev["event_number"].to_numpy()).get_indexer(ca["event_number"].to_numpy())
    assert (crow >= 0).all()
    spec = ca["Dmeson_species"].to_numpy(); pt = ca["Dmeson_pt"].to_numpy()
    csel = np.isin(spec, [s["code"] for s in SPECIES.values()])
    has = np.zeros(n0, bool); has[crow[csel]] = True
    keep = grl & npu0 & (afp | non)
    rowk = np.full(n0, -1, np.int64); rowk[keep] = np.arange(keep.sum())
    nk = int(keep.sum())
    cnt = np.zeros((nk, len(CATS)), np.uint16)
    ci = 0
    for sp, d in SPECIES.items():
        ms = spec == d["code"]; v = ca[d["var"]].to_numpy()
        pk = ms & (v >= d["peak"][0]) & (v <= d["peak"][1])
        sb = np.zeros(len(ca), bool)
        for lo, hi in d["sb"]:
            sb |= ms & (v >= lo) & (v <= hi)
        for lo, hi in PT_BINS:
            mp = (pt >= lo) & (pt < hi)
            for wm in (pk, sb):
                r = rowk[crow[wm & mp]]; r = r[r >= 0]
                cnt[:, ci] = np.bincount(r, minlength=nk).astype(np.uint16)
                ci += 1
    out = {"lb": ev["lumiBlock"].to_numpy()[keep].astype(np.int32), "mu": ev["actual_mu"].to_numpy()[keep].astype(np.float32),
           "nsel": ev["nSelBL"].to_numpy()[keep].astype(np.int32), "hasD": has[keep], "cnt": cnt,
           "afp": afp[keep], "non": non[keep], "evno": ev["event_number"].to_numpy()[keep]}
    for s in SIDES:
        out[f"xical_{s}"] = ev[f"xi_cal_{s}"].to_numpy()[keep].astype(np.float64)
        out[f"np_{s}"] = ev[f"afp_nProtons_{s}"].to_numpy()[keep].astype(np.int16)
        out[f"xi1_{s}"] = ev[f"afp_xi1_{s}"].to_numpy()[keep].astype(np.float64)
        out[f"tag_{s}"] = ev[f"afp_tag_{s}"].to_numpy()[keep] == 1
    for k, v in P.items():
        out["T_" + k] = v[keep]
    for c in ev.columns:
        if c.startswith("trigTBP_"):
            out[c] = ev[c].to_numpy()[keep] == 1
    out["info"] = np.array(json.dumps(info))
    np.savez(fn, **out)
    return out


def lb_groups(lb, ref_mask, n):
    """Bloki LB o równej liczbie zdarzeń odniesienia (ref_mask). → indeks grupy per zdarzenie (0..n−1)."""
    if n <= 1:
        return np.zeros(len(lb), np.int32)
    e = np.unique(np.quantile(lb[ref_mask], np.linspace(0, 1, n + 1)[1:-1]))
    return np.searchsorted(e, lb, side="right").astype(np.int32)


def pair_fracs(xc, wc, xp, wp, cuts_lo, cuts_hi):
    """Dla zdarzeń centralnych (ξ_cal xc, wagi wc) i protonów (ξ xp, wagi wp): Σ wag par z lo < Δ ≤ hi dla każdej pary (lo, hi),
    Δ = xc_i − xp_j. → (S[k] sumy, a[i,k] udział per zdarzenie centralne Σ_j wp_j 1[...], b[j,k] per proton Σ_i wc_i 1[...])."""
    o = np.argsort(xp); xps = xp[o]; cw = np.concatenate([[0.0], np.cumsum(wp[o])])
    oc = np.argsort(xc); xcs = xc[oc]; ccw = np.concatenate([[0.0], np.cumsum(wc[oc])])
    K = len(cuts_lo)
    a = np.zeros((len(xc), K)); b = np.zeros((len(xp), K))
    for k, (lo, hi) in enumerate(zip(cuts_lo, cuts_hi)):
        # lo < xc − xp ≤ hi  ⇔  xc − hi ≤ xp < xc − lo
        i1 = np.searchsorted(xps, xc - hi, side="left"); i2 = np.searchsorted(xps, xc - lo, side="left")
        a[:, k] = cw[i2] - cw[i1]
        # per proton: xp + lo < xc ≤ xp + hi
        j1 = np.searchsorted(xcs, xp + lo, side="right"); j2 = np.searchsorted(xcs, xp + hi, side="right")
        b[:, k] = ccw[j2] - ccw[j1]
    S = (a * wc[:, None]).sum(0)
    return S, a, b


def transfer(xc, wc, xp, wp, sr=(-np.inf, DELTA), cr=(DELTA + M_NOM, np.inf)):
    """R = T(SR)/T(CR) i błąd stat. (estymator ilorazowy: zdarzenia centralne ⊕ protony). → (R, σ_R, T_SR, T_CR)."""
    if len(xc) == 0 or len(xp) == 0:
        return np.nan, np.nan, 0.0, 0.0
    S, a, b = pair_fracs(xc, wc, xp, wp, [sr[0], cr[0]], [sr[1], cr[1]])
    if S[1] <= 0:
        return np.nan, np.nan, S[0], S[1]
    R = S[0] / S[1]
    vc = ((wc * (a[:, 0] - R * a[:, 1])) ** 2).sum() / S[1] ** 2
    vp = ((wp * (b[:, 0] - R * b[:, 1])) ** 2).sum() / S[1] ** 2
    return R, float(np.sqrt(vc + vp)), S[0], S[1]


def template_hist(xc, wc, xp, wp, edges):
    """Histogram Δ par (suma wag) w binach edges."""
    S, _, _ = pair_fracs(xc, wc, xp, wp, edges[:-1], edges[1:])
    return S


class Est:
    """Wynik: N_acc = Σ_c N_CR,c R_c z błędami (Poisson N_CR ⊕ stat szablonu)."""
    def __init__(self):
        self.n = 0.0; self.v_cr = 0.0; self.v_t = 0.0; self.ncr = 0.0; self.nmiss = 0.0; self.cells = []

    def add(self, ncr, n2cr, R, sR):
        """ncr = Σ wag danych w CR komórki, n2cr = Σ wag² (Poisson)."""
        if ncr == 0:
            return
        if not np.isfinite(R):
            self.nmiss += ncr; return
        self.n += ncr * R; self.v_cr += n2cr * R * R; self.v_t += (ncr * sR) ** 2; self.ncr += ncr
        self.cells.append((ncr, R, sR))

    @property
    def err(self):
        return float(np.sqrt(self.v_cr + self.v_t))


def estimate(d_delta, d_w, d_cell, c_x, c_w, c_cell, p_x, p_w, p_cell, m=M_NOM, sr=(-np.inf, DELTA), cr_hi=np.inf, cr_lo=None, pcell_of=None):
    """Dane w komórkach d_cell (Δ, wagi), pula centralna (c_x, c_w, c_cell), protony (p_x, p_w, p_cell — komórki protonów;
    pcell_of[c] mapuje komórkę danych/centralną na komórkę protonów, domyślnie tożsamość). CR = (cr_lo lub δ+m, cr_hi]."""
    lo = DELTA + m if cr_lo is None else cr_lo
    E = Est()
    incr = (d_delta > lo) & (d_delta <= cr_hi)
    cells = np.unique(d_cell[incr])
    for c in cells:
        md = incr & (d_cell == c)
        mc = c_cell == c
        pc = c if pcell_of is None else pcell_of[c]
        mp = p_cell == pc
        R, sR, _, _ = transfer(c_x[mc], c_w[mc], p_x[mp], p_w[mp], sr=sr, cr=(lo, cr_hi))
        E.add(d_w[md].sum(), (d_w[md] ** 2).sum(), R, sR)
    return E

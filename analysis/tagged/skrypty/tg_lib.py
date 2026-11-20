#!/usr/bin/env python3
"""TAGGED — wspólne: wczytanie tabel, selekcje (T2, regiony), kształt sygnału z fitu inkluzywnego, fit tagowany (Poisson, kształt stały),
profil wiarygodności, zabawki. Kryteria: ../KRYTERIA.md. Okno sygnału jest maskowane przy wczytaniu, dopóki blind=True."""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
W = os.path.dirname(TOP)
sys.path.insert(0, os.path.join(W, "FITBIAS_2026-11-19", "skrypty"))
import fb_fit as F  # noqa: E402
import fb_var as V  # noqa: E402

RUNS = os.environ.get("TG_RUNS", os.path.join(TOP, "runs", "interactive"))   # przebieg interaktywny lxplus 19.11 23:28–23:35 (job etapu ORCH `tagged` czekał w kolejce; ten sam kod)
OUT = os.environ.get("TG_OUT", os.path.join(TOP, "out"))
BITS = os.path.join(W, "TAGREG_2026-11-18", "runs", "lxplus")
BINS = ["3-5", "5-8", "8-30"]
SPN = ["D0", "Dplus"]
MODEL = {"D0": "M8", "Dplus": "M0"}
KSC = {"D0": 1.364, "Dplus": 1.068}
WIN = {"D0": (1.86484, 0.055), "Dplus": (1.86966, 0.045)}
RNG = {"D0": (1.70, 2.05), "Dplus": (1.70, 2.20)}
ALT = {"D0": (1.75, 2.00), "Dplus": (1.75, 2.10)}
DMU_S = {"D0": 0.0008, "Dplus": 0.0005}
TRIG = ["HLT_j20_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1MBTS_2",
        "HLT_noalg_L1AFP_A_OR_C", "HLT_noalg_L1J12", "HLT_noalg_L1RD0_FILLED", "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2"]
EPS = [round(0.95 - 0.05 * i, 2) for i in range(18)]
LUMI = {"427929": 0.4648, "R1": 0.4687, "R2": 27.65, "R3": 29.41}
FACC = {"427929": (0.0087, 0.0044), "428770": (0.048, 0.014), "AFP": (0.05, 0.03)}
PRED = {"D0": {"S": (0.00365, 0.0018, 0.0073), "M": (0.00194, 0.00077, 0.0077)}, "Dplus": {"S": (0.00237, 0.0012, 0.0047), "M": (0.00194, 0.00077, 0.0077)},
        "Dst": {"S": (0.0023, 0.0011, 0.0045), "M": (0.00194, 0.00077, 0.0077)}}


def tb(i):
    return 1 << TRIG.index(i)


def jwp(e):
    return EPS.index(e)


# ------------------------------------------------------------------------------------------------ wczytanie
def load_incl(run):
    z = np.load(os.path.join(RUNS, f"tg_incl_{run}.npz"), allow_pickle=False)
    return {k: z[k] for k in z.files if k != "meta"}, json.loads(str(z["meta"]))


def incl_h(I, sp, bn, cls, e):
    """hist 1 MeV pełnego zakresu, skumulowany WP ε."""
    return I[f"{sp}|{bn}|{cls}"][jwp(e) + 1:].sum(0).astype(np.float64)


def load_cand(run, blind=True):
    """Tabela kandydatów w zdarzeniach T0. blind=True: wiersze z oknem sygnału USUWANE przed zwróceniem."""
    z = np.load(os.path.join(RUNS, f"tg_cand_{run}.npz"), allow_pickle=False)
    C = {k: z[k] for k in z.files if k != "meta"}
    if blind:
        keep = np.ones(len(C["m"]), bool)
        for isp, sp in enumerate(SPN):
            m0, w = WIN[sp]
            keep &= ~((C["sp"] == isp) & (np.abs(C["m"] - m0) < w))
        C = {k: v[keep] for k, v in C.items()}
    a, c = C["tag_A"], C["tag_C"]
    C["T1"] = C["npu"] == 0
    C["xok"] = (~a | (C["xical_A"] <= C["xi_A"] + 0.02)) & (~c | (C["xical_C"] <= C["xi_C"] + 0.02))
    C["T2"] = C["T1"] & C["xok"]
    dA, dC = C["xical_A"] - C["xi_A"], C["xical_C"] - C["xi_C"]
    C["CR"] = C["T1"] & (~a | ((C["xical_A"] <= 0.12) & (dA > 0.045))) & (~c | ((C["xical_C"] <= 0.12) & (dC > 0.045)))
    C["base"] = C["ok"] & (C["ipt"] >= 0)
    if run == "428770":
        import pyarrow.parquet as pq
        ch = json.load(open(os.path.join(BITS, f"chains_{run}.json")))["chains"]
        t = pq.read_table(os.path.join(BITS, f"bits_{run}.parquet"), columns=["event_number", "trig_w0", "trig_w1"]).to_pandas()
        j = ch.index("HLT_j20_L1AFP_A_OR_C")
        w = t["trig_w0" if j < 64 else "trig_w1"].to_numpy().astype(np.uint64)
        j20 = ((w >> np.uint64(j % 64)) & np.uint64(1)) > 0
        en = t["event_number"].to_numpy().astype(np.uint64)
        o = np.argsort(en)
        ix = np.clip(np.searchsorted(en[o], C["event_number"]), 0, len(en) - 1)
        found = en[o][ix] == C["event_number"]
        C["found_bits"] = found
        C["j20afp"] = found & j20[o][ix]
    return C


def region_mask(C, run, reg):
    t = C["trig"]
    if run == "427929":
        mb = (t & tb("HLT_mb_mbts_L1MBTS_1")) > 0
        if reg == "MBTS":
            return mb
        if reg == "MBAFP":
            return mb & ((t & tb("HLT_noalg_L1AFP_A_OR_C")) > 0)
        if reg == "MBTS_A":
            return mb & C["tag_A"] & ~C["tag_C"]
        if reg == "MBTS_C":
            return mb & C["tag_C"] & ~C["tag_A"]
    else:
        r1 = (t & tb("HLT_noalg_L1AFP_A_OR_C")) > 0
        r2 = C["j20afp"] & ~r1
        j12 = (t & tb("HLT_noalg_L1J12")) > 0
        return {"R1": r1, "R2": r2, "R3": j12 & ~r1 & ~r2, "J12": j12}[reg]
    raise KeyError(reg)


INCL_CLS = {("427929", "MBTS"): "MB_npu0", ("427929", "MBAFP"): "MB_npu0", ("427929", "MBTS_A"): "MB_npu0", ("427929", "MBTS_C"): "MB_npu0",
            ("428770", "R1"): "all_npu0", ("428770", "R2"): "all_npu0", ("428770", "R3"): "all_npu0", ("428770", "J12"): "J12_npu0"}


def edges5(sp, rng=None):
    lo, hi = rng or RNG[sp]
    return np.arange(lo, hi + F.BW / 2, F.BW)


def hist5(m, sp, rng=None):
    return np.histogram(m, edges5(sp, rng))[0].astype(np.float64)


def winmask(sp, rng=None):
    e = edges5(sp, rng)
    ctr = 0.5 * (e[1:] + e[:-1])
    m0, w = WIN[sp]
    return np.abs(ctr - m0) < w


# ------------------------------------------------------------------------------------------------ fit inkluzywny → kształt
def incl_fit(I, sp, cls, e):
    """Fity inkluzywne per bin pT modelem nominalnym (D0 M8, D+ M0). → lista dict(N, err, dmu, sc, cat)."""
    if V.Z is None:
        V._init()
    out = []
    for bn in BINS:
        cat = V.wcat(sp, bn, e)
        cat["kap"] *= KSC[sp]
        H5 = V.full5(incl_h(I, sp, bn, cls, e))[None, :]
        f = V.fit_set(H5, [cat], sp, MODEL[sp])
        out.append({"bin": bn, "N": f["N"][0], "err": (f.get("err") or [float("nan")])[0], "dmu": f["dmu"], "sc": f["sc"], "valid": f["valid"], "chi2": f["chi2"], "ndf": f["ndf"],
                    "frefl": cat["frefl"], "kap": cat["kap"], "cat": cat})
    return out


def sig_shape(sp, fits, rng=None, w=None, dsc=1.0, ddmu=0.0, frs=1.0, ks=1.0, bins=None):
    """Kształt łączny Σ_b w_b S_b w binach 5 MeV zakresu rng; N_S = współczynnik przy unormowanym szablonie sygnału."""
    e = edges5(sp, rng)
    m0 = F.SP[sp]["m0"]
    use = [i for i, f in enumerate(fits) if bins is None or f["bin"] in bins]
    if w is None:
        w = np.array([max(fits[i]["N"], 0.0) for i in use])
    w = np.asarray(w, float) / np.sum(w)
    S = np.zeros(len(e) - 1)
    for wi, i in zip(w, use):
        f = fits[i]
        c = f["cat"]
        dmu = f["dmu"] + ddmu
        s = c["sigT"].bins(e, m0, dmu, f["sc"] * dsc)
        if c["refT"] is not None:
            s = s + c["frefl"] * frs * c["refT"].bins(e, m0, dmu, 1.0)
        if MODEL[sp] == "M8":
            s = s + f["kap"] * ks * c["peakT"].bins(e, m0, dmu, 1.0)
        S += wi * s
    return S


# ------------------------------------------------------------------------------------------------ fit tagowany
class TFit:
    """μ_i = N_S·S_i + N_B·b_i/Σb, b = exp(Czebyszew n). mask = biny użyte (None = wszystkie)."""

    def __init__(self, S, nb, nbkg, mask=None, bkind="exp", x=None):
        self.S = np.asarray(S, float)
        self.n = nbkg
        self.use = np.ones(nb, bool) if mask is None else mask
        self.x = np.linspace(-1 + 1.0 / nb, 1 - 1.0 / nb, nb) if x is None else x
        self.bkind = bkind

    def bsh(self, a):
        if self.bkind == "pow":                                   # D*: (Δm − m_π)^a, x = Δm − m_π
            sh = np.power(np.maximum(self.x, 1e-9), a[0])
        else:
            sh = F.bshape(self.x, a, "exp")
        return sh / sh[self.use].sum()

    def mu(self, p):
        return np.maximum(p[0] * self.S + p[1] * self.bsh(p[2:]), 1e-9)

    def nll2(self, p, H):
        mu = self.mu(p)[self.use]
        n = H[self.use]
        return 2.0 * float(np.sum(mu - n) + np.sum(np.where(n > 0, n * np.log(np.where(n > 0, n, 1.0) / mu), 0.0)))

    def fit(self, H, fixNS=None, start=None):
        from iminuit import Minuit
        ntot = float(H[self.use].sum())
        a0 = [0.5] if self.bkind == "pow" else [0.0] * self.n
        p0 = [0.0 if fixNS is None else fixNS, max(ntot, 1.0)] + a0
        if start is not None:
            p0 = list(start)
            if fixNS is not None:
                p0[0] = fixNS
        m = Minuit(lambda p: self.nll2(p, H), p0)
        m.errordef = 1.0
        m.limits[0] = (-3 * math.sqrt(ntot + 1) - 20, ntot + 20)
        m.limits[1] = (0.0, 2 * ntot + 20)
        for k in range(2, len(p0)):
            m.limits[k] = (0.0, 5.0) if self.bkind == "pow" else (-30.0, 30.0)
        if fixNS is not None:
            m.fixed[0] = True
        m.migrad(ncall=5000)
        if not m.valid:
            m.simplex(ncall=5000); m.migrad(ncall=10000)
        return np.array(m.values), float(m.fval), bool(m.valid), m

    def full(self, H, cl95=True):
        """N̂, przedziały z profilu (68 %, 95 %), q0, Z, HESSE."""
        from scipy.optimize import brentq
        v, f0, ok, m = self.fit(H)
        try:
            m.hesse(); he = float(m.errors[0])
        except Exception:  # noqa: BLE001
            he = float("nan")
        if not np.isfinite(he) or he <= 0:
            he = math.sqrt(max(H[self.use].sum(), 1.0))
        cache = {}

        def q(ns):
            if ns not in cache:
                cache[ns] = self.fit(H, fixNS=ns, start=v)[1] - f0
            return cache[ns]

        def cross(level, sgn):
            a, step = v[0], he * (1.0 if level < 2 else 2.0)
            b = a + sgn * step
            for _ in range(12):
                if q(b) > level:
                    break
                b += sgn * step
            else:
                return float("nan")
            try:
                return float(brentq(lambda t: q(t) - level, min(a, b), max(a, b), xtol=1e-3 * he + 1e-6))
            except Exception:  # noqa: BLE001
                return float("nan")
        r = {"N": float(v[0]), "NB": float(v[1]), "a": [float(t) for t in v[2:]], "hesse": he, "valid": ok, "fval": f0,
             "lo68": cross(1.0, -1), "hi68": cross(1.0, +1)}
        if cl95:
            r["lo95"], r["hi95"] = cross(3.84, -1), cross(3.84, +1)
        q0 = max(q(0.0), 0.0) if v[0] > 0 else 0.0
        r["q0"], r["Z"] = float(q0), float(math.sqrt(q0))
        mu = self.mu(v)
        r["chi2"] = float(np.sum(((H - mu) ** 2 / mu)[self.use])); r["ndf"] = int(self.use.sum() - len(v))
        r["x"] = [float(t) for t in v]
        return r


def bkg_order(S, Hsb, mask_sb):
    """Rząd tła z samych pasm: n = 2 gdy 2ΔlnL(2 vs 1) > 3,84. → (n, q, parametry tła n)."""
    nb = len(Hsb)
    r = {}
    for n in (1, 2):
        t = TFit(S, nb, n, mask=mask_sb)
        v, f, ok, _ = t.fit(Hsb, fixNS=0.0)
        r[n] = (v, f)
    dq = r[1][1] - r[2][1]
    n = 2 if dq > 3.84 else 1
    return n, float(dq), r[n][0]


def toy_worker(args):
    """Zabawki: tło z fitu pasm (ekstrapolowane pod okno) + wstrzyknięcie nS. → [N̂, lo68, hi68, hesse, q0]"""
    S, nb, n, vbkg, mask_sb, ninj, ntoy, seed, bkind, x = args
    rng = np.random.default_rng(seed)
    tsb = TFit(S, nb, n, mask=mask_sb, bkind=bkind, x=x)
    b = tsb.bsh(vbkg[2:]) * vbkg[1]                                # unormowane do pasm → pełny zakres
    t = TFit(S, nb, n, bkind=bkind, x=x)
    out = []
    for _ in range(ntoy):
        H = rng.poisson(np.maximum(b + ninj * S, 0.0)).astype(float)
        r = t.full(H, cl95=False)
        out.append([r["N"], r["lo68"], r["hi68"], r["hesse"], r["q0"]])
    return np.array(out)


def toy_summary(T, ninj):
    N, lo, hi, he, q0 = T.T
    sig = np.where(N > ninj, N - lo, hi - N)
    pull = (N - ninj) / np.where(np.isfinite(sig) & (sig > 0), sig, he)
    cov = np.mean((lo <= ninj) & (ninj <= hi))
    return {"ninj": ninj, "n": int(len(N)), "bias": float(np.mean(N) - ninj), "bias_err": float(np.std(N) / math.sqrt(len(N))), "rmsN": float(np.std(N)),
            "pull_mean": float(np.nanmean(pull)), "pull_rms": float(np.nanstd(pull)), "cov68": float(cov), "medZ": float(np.median(np.sqrt(q0))),
            "fracZ3": float(np.mean(np.sqrt(q0) >= 3)), "nan": int(np.sum(~np.isfinite(lo) | ~np.isfinite(hi)))}

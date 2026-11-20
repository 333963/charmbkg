#!/usr/bin/env python3
"""Wydajność weta wierzchołków `nPileupVertices == 0` na próbce INKLUZYWNEJ (bez tagu protonu AFP) i na MC.
Wejście (tylko odczyt): B_trig/runs/lxplus/{bt,trg}_*.npz (h5 v2.0 → etap ORCH btrig, 17.11). Wyjście: out/veto.json, out/tables.md.
Nie używa żadnej zmiennej AFP (protony, ξ). Klasa trygera „other” (ani J12, ani MBTS; głównie łańcuchy L1AFP) — tylko ułamek zdarzeń-tła,
bez uzysku D*. Zmienne zdarzeniowe brane z wiersza kandydata (e_*, all_e_*): w fd_nd klucz (file_idx, event) tabeli ev_ nie zgadza się
z tabelami kandydatów (v0d_explore.py), więc tabela ev_ służy tylko do zliczeń łącznych."""
import sys, os, json, time
sys.dont_write_bytecode = True
import numpy as np
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CW, "B_trig", "ISO_origin"))
import iso_common as I  # noqa  (sbs_w: wagi odejmowania tła z pasm Δm; C.classify: klasy trygera)
C = I.C
R = os.path.join(CW, "B_trig", "runs", "lxplus")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(20261118)
NB = 300
MWIN, STHR = 0.040, 0.6                      # |m(Kπ) − m_D0| < 40 MeV, H > 0,6 (jak B_trig/ISO_origin „base”)
PT_BINS = [(3, 5), (5, 8), (8, 30)]
NSEL_E = [0, 10, 20, 30, 40, 50, 65, 80, 1000]
LXY_E = [0.25, 0.4, 0.6, 1.0, 2.0, 1e3]
MU_E = {428770: [0, 0.025, 0.035, 0.045, 0.055, 0.065, 0.08, 0.12, 0.3], 427929: [0, 0.003, 0.004, 0.005, 0.007, 0.012, 0.02, 0.03]}
LB_E = {428770: [0, 100, 250, 290, 400], 427929: [119, 170, 220, 270, 325]}
res = {}
md = []


def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); md.append(s)


def binom(k, n):
    if n <= 0:
        return (np.nan, np.nan, 0)
    p = k / n
    v = p * (1 - p) if 0 < k < n else 1.0 / (n + 1)        # k = 0 lub n: ~1/n zamiast zera
    return (float(p), float(np.sqrt(v / n)), int(n))


def ev_frac(passed, mask=None):
    """Ułamek zdarzeń przechodzących weto; błąd dwumianowy (równoważny bootstrapowi po zdarzeniach dla wag jednostkowych — test w logu)."""
    if mask is not None:
        passed = passed[mask]
    return binom(int(passed.sum()), len(passed))


def boot_eps(uid, w, pas, raw=False):
    """ε = Σ w·pass / Σ w po kandydatach; bootstrap Poissona po ZDARZENIACH (wagi SBS stałe)."""
    if len(w) == 0 or w.sum() < 25:
        return (np.nan, np.nan, float(w.sum()) if len(w) else 0.0)
    u, inv = np.unique(uid, return_inverse=True)
    W = np.bincount(inv, w, len(u)); Wp = np.bincount(inv, w * pas, len(u))
    e0 = Wp.sum() / W.sum()
    g = RNG.poisson(1.0, (NB, len(u)))
    eb = (g @ Wp) / np.maximum(g @ W, 1e-9)
    if raw:
        return e0, eb
    return (float(e0), float(eb.std(ddof=1)), float(W.sum()))


def fmt(t):
    e, s, n = t
    return "–" if not np.isfinite(e) else f"{100*e:.2f} ± {100*s:.2f}"


def ff(f):
    return [f"{100*f['f']:.2f} ± {100*f['f_err']:.2f}", f"{f['eps_v']:.3f} ± {f['eps_v_err']:.3f}", f"{f['rho']:.2f}", f"{f['chi2']:.1f}/{f['ndf']}"]


def uid_of(fi, ev, ns, lp):
    """Identyfikator zdarzenia odporny na powtarzające się numery zdarzeń w MC (fd_nd: 1…10000 w każdym pliku)."""
    return (fi.astype(np.int64) << 44) ^ (ev.astype(np.int64) << 12) ^ (ns.astype(np.int64) << 2) ^ (np.round(lp * 1000).astype(np.int64) << 24)


def tab(z, pre, extra=()):
    d = {"uid": uid_of(z[pre + "file_idx"], z[pre + "event"], z[pre + "e_nSelBL"], z[pre + "e_lead_pt"])}
    for k in ("npu", "actual_mu", "lb", "grl", "nSelBL"):
        d[k] = z[pre + "e_" + k]
    for k in z.files:
        if k.startswith(pre + "e_tp_"):
            d[k[len(pre) + 2:]] = z[k]
    for k in extra:
        d[k] = z[pre + k]
    d["pas"] = d["npu"] == 0
    return d


def pois_fit(mu, k, n):
    """k_i ~ Bin(n_i, (1−f)·exp(−ε_v·μ_i)); zwraca f, ε_v, błędy (hesjan numeryczny), χ²/ndf."""
    mu, k, n = map(lambda a: np.asarray(a, float), (mu, k, n))

    def nll(p):
        q = np.clip((1 - p[0]) * np.exp(-p[1] * mu), 1e-9, 1 - 1e-9)
        return -(k * np.log(q) + (n - k) * np.log(1 - q)).sum()
    r = minimize(nll, [0.01, 0.8], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 4000})
    p = r.x; h = np.array([1e-4, 1e-3]); H = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            ei, ej = np.eye(2)[i] * h[i], np.eye(2)[j] * h[j]
            H[i, j] = (nll(p + ei + ej) - nll(p + ei - ej) - nll(p - ei + ej) + nll(p - ei - ej)) / (4 * h[i] * h[j])
    try:
        cov = np.linalg.inv(H); err = np.sqrt(np.maximum(np.diag(cov), 0)); rho = cov[0, 1] / max(err[0] * err[1], 1e-30)
    except np.linalg.LinAlgError:
        err = np.array([np.nan, np.nan]); rho = np.nan
    q = (1 - p[0]) * np.exp(-p[1] * mu); chi2 = (((k - n * q) ** 2) / np.maximum(n * q * (1 - q), 1e-9)).sum()
    return {"f": float(p[0]), "f_err": float(err[0]), "eps_v": float(p[1]), "eps_v_err": float(err[1]), "rho": float(rho),
            "chi2": float(chi2), "ndf": int(len(mu) - 2)}


def mu_bins_fit(mu, pas, nb=16):
    ed = np.unique(np.quantile(mu, np.linspace(0, 1, nb + 1))); ed[-1] += 1e-9
    ib = np.clip(np.digitize(mu, ed) - 1, 0, len(ed) - 2)
    n = np.bincount(ib, minlength=len(ed) - 1); k = np.bincount(ib, pas, len(ed) - 1); m = np.bincount(ib, mu, len(ed) - 1) / np.maximum(n, 1)
    ok = n > 50
    return pois_fit(m[ok], k[ok], n[ok])


def table(title, rows, header):
    P(f"\n**{title}**\n"); P("| " + " | ".join(header) + " |"); P("|" + "---|" * len(header))
    for r in rows:
        P("| " + " | ".join(str(x) for x in r) + " |")



#!/usr/bin/env python3
"""DNET etap 4 — metoda R: względny współczynnik skali sieć/drzewo z INKLUZYWNEGO piku D0 / D+ (kryteria: DNET/KRYTERIA_etap4.md).

  R = [N_S(sieć, WP)/ε_MC(sieć)] / [N_S(drzewo H, WP)/ε_MC(H)]   na TEJ SAMEJ próbce kandydatów → SF_sieć = R × SF_H (SF_H z tagu D*).

Kandydaci dzieleni na kategorie ROZŁĄCZNE: b = przechodzi oba, n = tylko sieć, t = tylko drzewo. Jeden fit łączny trzech widm masy
(rozszerzony Poisson, biny 5 MeV): w każdej kategorii N_S·[szablon sygnału MC + f_refl·szablon odbicia K↔π (D0; f_refl z MC)] + tło
exp(wielomian Czebyszewa); wspólne przesunięcie i skala szerokości piku (dmu, sc). R = [(N_b+N_n)/(e_b+e_n)] / [(N_b+N_t)/(e_b+e_t)],
e_c = udział sygnału MC (fd_nd) w kategorii; błąd z macierzy kowariancji (kategorie rozłączne ⇒ statystycznie niezależne) oraz z
bootstrapu Poissona po ZDARZENIACH (wielu kandydatów na zdarzenie); σ = max(HESSE, bootstrap).

  dn4_ana.py r --app DIR/app/d150 --out dnet4_R.json [--nboot 100] [--force-data]
      1) domknięcie na MC: pseudo-dane = mb_nd (wszystkie kategorie prawdy), ε i szablony z fd_nd → R_fit wobec R_prawda (liczenie prawdy);
      2) dane (428770: all, J12only, reszta; 427929: MB) — TYLKO dla (gatunek, model, bin, WP) z zaliczonym domknięciem (bramka).
  dn4_ana.py d23 --app DIR/app/d23 --r dnet4_R.json --out dnet4_d23.json      bin 2–3 GeV: WP z MC + pasm (okno zaślepione), potem fit piku
  dn4_ana.py report --res dnet4_R.json [--d23 dnet4_d23.json]                 tabele markdown (stdout)
  dn4_ana.py selftest [--dir KATALOG]                                         zabawki ze znanym R (numpy + iminuit; bez torch)
Pliki app: app_<próbka>_<k>.npz z dn4_apply.py (logit float16, progi liczone na tych samych wartościach). Bez AFP, tylko pojedyncze D."""
import argparse
import datetime
import glob
import json
import math
import os
import re
import sys
import time

import numpy as np

BINS = [(3, 5), (5, 8), (8, 30)]
WPS = [0.7, 0.5, 0.8]                       # 0,7 = nominalny (bootstrap + warianty); 0,5 i 0,8 opisowo (HESSE)
SP = {"D0": {"m0": 1.86484, "win": 0.055, "fit": (1.75, 2.00), "alt": (1.70, 2.05), "nbkg": 2, "refl": True},
      "Dplus": {"m0": 1.86966, "win": 0.045, "fit": (1.70, 2.20), "alt": (1.75, 2.10), "nbkg": 3, "refl": False}}
def ext_range(sp):
    """Zakres mas ekstrakcji w węzłach app (dn4_extract: DN4_MR_<gatunek> albo cfg mass_range) — zakres szablonów MC."""
    v = os.environ.get("DN4_MR_" + sp, "")
    return tuple(float(x) for x in v.split(",")) if v else {"D0": (1.70, 2.05), "Dplus": (1.70, 2.20)}[sp]


BW = 0.005                                  # szerokość binu masy [GeV]
GRID = 0.002                                # siatka dystrybuanty szablonu MC [GeV] (interpolacja monotoniczna PCHIP → gładka pochodna)
DMU_MAX, SC_LIM = 0.006, (0.7, 1.5)
REF = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
CLOS_NSIG, CLOS_MAXERR, CLOS_MAXBIAS = 2.5, 0.05, 0.04   # bramka domknięcia: |R_fit − R_prawda| ≤ 2,5σ i ≤ 0,04 oraz σ ≤ 0,05
D23_WPS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
D23_SCALE = 0.6                             # oczekiwany stosunek uzysku 2–3 / 3–5 GeV (D23_plan/RAPORT.md §3)


def bname(b):
    return f"{b[0]}-{b[1]}"


# ----------------------------------------------------------------------------------------------------------- wejście
def load_app(appdir, sample, sp, need=("uid", "m", "pt", "ok", "cls", "truth", "swap", "w", "lg")):
    fs = sorted(glob.glob(os.path.join(appdir, f"app_{sample}_*.npz")))
    cols, tags = {}, None
    for f in fs:
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        t = meta["tags"].get(sp) or []
        if not t or f"{sp}__m" not in z.files:
            continue
        if tags is None:
            tags = t
        elif t != tags:
            raise RuntimeError(f"{f}: inne modele ({t}) niż w poprzednich częściach ({tags})")
        for c in need:
            cols.setdefault(c, []).append(z[f"{sp}__{c}"])
    if tags is None:
        return None
    T = {c: np.concatenate(v) for c, v in cols.items()}
    T["tags"], T["nfiles"] = tags, len(fs)
    return T


def classes(sample, cls):
    grl, mb1, j12, mbts = (cls & 1) > 0, (cls & 2) > 0, (cls & 4) > 0, (cls & 8) > 0
    if sample == "data_427929":
        return {"MB427929": mb1 & grl}
    if sample == "data_428770":
        return {"all": np.ones(len(cls), bool), "J12only": j12 & ~mbts, "notJ12only": ~(j12 & ~mbts)}
    return {"all": np.ones(len(cls), bool)}


# ----------------------------------------------------------------------------------------------------------- szablony i fit
class Tmpl:
    """Szablon masy z MC jako dystrybuanta (siatka 2 MeV, interpolacja monotoniczna PCHIP); całka po całym zakresie ekstrakcji MC = 1."""

    def __init__(self, m, w, lo, hi):
        from scipy.interpolate import PchipInterpolator
        self.g = np.arange(lo - 0.02, hi + 0.02 + GRID / 2, GRID)
        h, _ = np.histogram(m, self.g, weights=w)
        tot = h.sum()
        self.n = float(tot)
        c = np.concatenate([[0.0], np.cumsum(h) / tot]) if tot > 0 else np.linspace(0, 1, len(self.g))
        self.f = PchipInterpolator(self.g, c)

    def bins(self, edges, m0, dmu=0.0, sc=1.0):
        x = np.clip((edges - m0 - dmu) / sc + m0, self.g[0], self.g[-1])
        return np.diff(self.f(x))


def cheb(x, a):
    """Σ_k a_k T_k(x), k = 1..len(a)."""
    t0, t1, s = np.ones_like(x), x, np.zeros_like(x)
    for k, ak in enumerate(a):
        s = s + ak * t1
        t0, t1 = t1, 2 * x * t1 - t0
    return s


def bshape(x, a, kind="exp"):
    pol = cheb(x, a)
    return np.exp(np.clip(pol, -50, 50)) if kind == "exp" else np.maximum(1.0 + pol, 1e-9)


def fit_cats(H, sigT, refT, frefl, edges, m0, nbkg, kind="exp", fix_shape=False, start=None, blind=None, win=0.045):
    """Fit łączny kategorii. H [ncat, nb]; sigT/refT: listy szablonów (refT[c] None → bez odbicia); frefl [ncat].
    blind = maska binów wyłączonych (fit samych pasm: N_S = 0, kształt stały). → dict (N, cov_N, NB, dmu, sc, valid, chi2, ndf, x)."""
    from iminuit import Minuit
    ncat, nb = H.shape
    ctr = 0.5 * (edges[1:] + edges[:-1])
    x = 2 * (ctr - edges[0]) / (edges[-1] - edges[0]) - 1
    use = np.ones(nb, bool) if blind is None else ~blind
    npc = 2 + nbkg

    def model(p):
        dmu, sc = p[0], p[1]
        mu = np.empty((ncat, nb))
        for c in range(ncat):
            q = p[2 + c * npc: 2 + (c + 1) * npc]
            s = sigT[c].bins(edges, m0, dmu, sc)
            if refT[c] is not None and frefl[c] > 0:
                s = s + frefl[c] * refT[c].bins(edges, m0, dmu, 1.0)
            sh = bshape(x, q[2:], kind)
            mu[c] = q[0] * s + q[1] * sh / sh[use].sum()
        return np.maximum(mu, 1e-9)

    def cost(p):
        mu = model(p)[:, use]
        n = H[:, use]
        return 2.0 * float(np.sum(mu - n) + np.sum(np.where(n > 0, n * np.log(np.where(n > 0, n, 1.0) / mu), 0.0)))
    names = ["dmu", "sc"]
    p0 = [0.0, 1.0]
    lim = [(-DMU_MAX, DMU_MAX), SC_LIM]
    inwin = np.abs(ctr - m0) < 0.03
    for c in range(ncat):
        ntot = float(H[c][use].sum())
        nwin = float(H[c][inwin & use].sum())
        nsb = float(H[c][~inwin & use].sum()) * inwin.sum() / max((~inwin & use).sum(), 1)
        s0 = 0.0 if blind is not None else max(nwin - nsb, 0.0) * 1.2
        names += [f"NS{c}", f"NB{c}"] + [f"a{c}_{k + 1}" for k in range(nbkg)]
        p0 += [s0, max(ntot - s0, 1.0)] + [0.0] * nbkg
        lim += [(-5 * math.sqrt(nwin + 1) - 10, 1.5 * ntot + 20), (0.0, 2.0 * ntot + 20)] + [(-30.0, 30.0)] * nbkg
    if start is not None and len(start) == len(p0):
        p0 = list(start)
    m = Minuit(cost, p0, name=names)
    m.errordef = 1.0
    m.strategy = 1
    for nme, l in zip(names, lim):
        m.limits[nme] = l
    m.fixed["dmu"] = True; m.fixed["sc"] = True
    if blind is not None:
        for c in range(ncat):
            m.values[f"NS{c}"] = 0.0; m.fixed[f"NS{c}"] = True
    def good():
        fm = m.fmin
        return fm is not None and not fm.has_reached_call_limit and (m.valid or fm.edm < 0.02) and not fm.is_above_max_edm
    m.migrad(ncall=20000)
    if not fix_shape and blind is None:
        m.fixed["dmu"] = False; m.fixed["sc"] = False
        m.migrad(ncall=20000)
    for it in range(3):                         # ratunek: simplex → migrad (potem strategia 2)
        if good():
            break
        if it == 2:
            m.strategy = 2
        m.simplex(ncall=20000)
        m.migrad(ncall=40000)
    cov = None
    try:
        m.hesse()
        fm = m.fmin
        if m.covariance is not None and fm.has_covariance and fm.has_posdef_covar and not fm.has_made_posdef_covar:
            cov = np.array(m.covariance)
    except Exception:  # noqa: BLE001
        cov = None
    v = np.array(m.values)
    iN = [2 + c * npc for c in range(ncat)]
    mu = model(v)
    ndf = int(use.sum() * ncat - m.nfit)
    pear = float(np.sum(((H - mu) ** 2 / mu)[:, use]))
    out = {"N": [float(v[i]) for i in iN], "NB": [float(v[i + 1]) for i in iN], "dmu": float(v[0]), "sc": float(v[1]),
           "valid": bool(good() and cov is not None), "accurate": bool(m.accurate), "edm": float(m.fmin.edm), "fval": float(m.fval), "chi2": pear, "ndf": ndf, "x": v.tolist(),
           "at_limit": bool((not m.fixed["dmu"]) and (abs(v[0]) > 0.98 * DMU_MAX or v[1] < SC_LIM[0] * 1.01 or v[1] > SC_LIM[1] * 0.99))}
    if cov is not None:
        out["cov_N"] = cov[np.ix_(iN, iN)].tolist()
        out["dmu_err"], out["sc_err"] = float(math.sqrt(max(cov[0, 0], 0))), float(math.sqrt(max(cov[1, 1], 0)))
    wmask = np.abs(ctr - m0) < win
    out["bkg_win"] = []
    for i in iN:
        sh = bshape(x, v[i + 2:i + 2 + nbkg], kind)
        out["bkg_win"].append(float(v[i + 1] * sh[wmask].sum() / sh[use].sum()))
    return out


def r_of(N, cov, e):
    """R, σ_R (z kowariancji N), uzyski, zysk w ficie. N = [N_b, N_n, N_t], e = [e_b, e_n, e_t]."""
    Nn, Nt = N[0] + N[1], N[0] + N[2]
    en, et = e[0] + e[1], e[0] + e[2]
    if Nn <= 0 or Nt <= 0 or en <= 0 or et <= 0:
        return None
    R = (Nn / en) / (Nt / et)
    o = {"R": float(R), "N_net": float(Nn), "N_tree": float(Nt), "y_net": float(Nn / en), "y_tree": float(Nt / et)}
    if cov is not None:
        C = np.asarray(cov)
        g = np.array([R * (1 / Nn - 1 / Nt), R / Nn, -R / Nt])
        o["R_err"] = float(math.sqrt(max(g @ C @ g, 0)))
        sn = math.sqrt(max(C[0, 0] + C[1, 1] + 2 * C[0, 1], 0)); st = math.sqrt(max(C[0, 0] + C[2, 2] + 2 * C[0, 2], 0))
        o.update({"N_net_err": sn, "N_tree_err": st, "relerr_net": sn / Nn, "relerr_tree": st / Nt,
                  "G_fit": float((Nn / sn) / (Nt / st)) if sn > 0 and st > 0 else None})
    return o


class MC:
    """Składniki z sygnału MC (fd_nd) dla (gatunek, model, odniesienie, bin, WP): progi, udziały kategorii, szablony."""

    def __init__(self, sp, F, tag, ref, pb, wp, lo_hi=None):
        c = SP[sp]
        self.sp, self.ok = sp, False
        lo, hi = lo_hi if lo_hi is not None else ext_range(sp)
        self.mc_range = (lo, hi)
        jn, jr = F["tags"].index(tag), F["tags"].index(ref)
        inb = (F["pt"] >= pb[0]) & (F["pt"] < pb[1]) & F["ok"] & np.isin(F["truth"], (0, 1))
        sig = inb & ((F["swap"] != 1) if sp == "D0" else True)
        sw = sig & (np.abs(F["m"] - c["m0"]) < c["win"])
        if sw.sum() < 200:
            return
        self.thr = (float(np.quantile(F["lg"][sw, jn].astype(np.float64), 1 - wp)), float(np.quantile(F["lg"][sw, jr].astype(np.float64), 1 - wp)))
        cat = self.cat(F["lg"][:, jn], F["lg"][:, jr])
        w = F["w"].astype(np.float64)
        W = float(w[sig].sum())
        self.k = [float(w[sig & (cat == j)].sum()) for j in range(3)]
        self.kn = [int((sig & (cat == j)).sum()) for j in range(3)]
        self.e = [k / W for k in self.k]
        self.n_sig = int(sig.sum())
        self.eps_net, self.eps_tree = self.e[0] + self.e[1], self.e[0] + self.e[2]
        k0, k1, k2 = self.kn
        self.eps_ratio_relerr = math.sqrt(k1 / max(k0 + k1, 1) ** 2 + k2 / max(k0 + k2, 1) ** 2)
        allT = Tmpl(F["m"][sig], w[sig], lo, hi)
        self.sigT = [Tmpl(F["m"][sig & (cat == j)], w[sig & (cat == j)], lo, hi) if self.kn[j] >= 300 else allT for j in range(3)]
        self.sigAny = {"net": Tmpl(F["m"][sig & (cat < 2)], w[sig & (cat < 2)], lo, hi), "tree": Tmpl(F["m"][sig & ((cat == 0) | (cat == 2))], w[sig & ((cat == 0) | (cat == 2))], lo, hi)}
        self.refT, self.frefl = [None] * 3, [0.0] * 3
        self.refAny, self.freflAny = {"net": None, "tree": None}, {"net": 0.0, "tree": 0.0}
        if c["refl"]:
            rf = inb & (F["swap"] == 1)
            if rf.sum() >= 300:
                allR = Tmpl(F["m"][rf], w[rf], lo, hi)
                for j in range(3):
                    s_ = rf & (cat == j)
                    self.frefl[j] = float(w[s_].sum() / max(self.k[j], 1e-9))
                    self.refT[j] = Tmpl(F["m"][s_], w[s_], lo, hi) if s_.sum() >= 300 else allR
                for nm, s_ in (("net", rf & (cat < 2)), ("tree", rf & ((cat == 0) | (cat == 2)))):
                    self.freflAny[nm] = float(w[s_].sum() / max(self.k[0] + self.k[1 if nm == "net" else 2], 1e-9))
                    self.refAny[nm] = Tmpl(F["m"][s_], w[s_], lo, hi) if s_.sum() >= 300 else allR
        self.ok = min(self.kn[0] + self.kn[1], self.kn[0] + self.kn[2]) >= 200

    def cat(self, lgn, lgr):
        pn, pr = lgn >= self.thr[0], lgr >= self.thr[1]
        return np.where(pn & pr, 0, np.where(pn, 1, np.where(pr, 2, 3))).astype(np.int8)


def hist_cats(m, cat, edges, w=None):
    nb = len(edges) - 1
    ib = np.floor((m - edges[0]) / BW + 1e-9).astype(np.int64)
    ok = (ib >= 0) & (ib < nb) & (cat < 3)
    idx = cat[ok].astype(np.int64) * nb + ib[ok]
    return np.bincount(idx, weights=None if w is None else w[ok], minlength=3 * nb).reshape(3, nb).astype(np.float64), ok, idx


def measure(sp, mc, m, cat, uid, nboot, rng, full=True):
    """Pomiar R na wierszach (m, cat, uid) jednej klasy i binu pT. full=True: bootstrap po zdarzeniach + warianty (tło, zakres, kształt)."""
    c = SP[sp]
    lo, hi = c["fit"]
    edges = np.arange(lo, hi + BW / 2, BW)
    H, ok, idx = hist_cats(m, cat, edges)
    f = fit_cats(H, mc.sigT, mc.refT, mc.frefl, edges, c["m0"], c["nbkg"])
    o = {"n_rows": [int(x) for x in H.sum(1)], "fit": {k: f[k] for k in ("N", "NB", "dmu", "sc", "valid", "accurate", "chi2", "ndf", "at_limit", "bkg_win")},
         "eps": mc.e, "eps_net": mc.eps_net, "eps_tree": mc.eps_tree, "thr": mc.thr, "eps_ratio_relerr": mc.eps_ratio_relerr, "determined": False}
    r = r_of(f["N"], f.get("cov_N"), mc.e) if f["valid"] else None
    if r is None or "R_err" not in r:
        return o
    o.update(r)
    o["R_err_hesse"] = r["R_err"]
    if not full:
        o["determined"] = bool(f["valid"] and not f["at_limit"])
        return o
    # bootstrap po zdarzeniach
    if nboot > 0:
        u, inv = np.unique(uid[ok], return_inverse=True)
        Rb, nok = [], 0
        for _ in range(nboot):
            wv = rng.poisson(1.0, len(u)).astype(np.float64)[inv]
            Hb = np.bincount(idx, weights=wv, minlength=H.size).reshape(H.shape)
            fb = fit_cats(Hb, mc.sigT, mc.refT, mc.frefl, edges, c["m0"], c["nbkg"], start=f["x"])
            rb = r_of(fb["N"], None, mc.e) if fb["valid"] else None
            if rb is not None:
                Rb.append(rb["R"]); nok += 1
        o["boot"] = {"n": nboot, "ok": nok, "rms": float(np.std(Rb)) if nok > 5 else None, "mean": float(np.mean(Rb)) if nok else None,
                     "cand_per_event": float(ok.sum() / max(len(u), 1))}
        if nok > 5:
            o["R_err"] = float(max(o["R_err_hesse"], np.std(Rb)))
    # warianty: tło (rząd ±1, wielomian bez exp), zakres fitu, kształt stały
    var = {}
    for nm, kw in ((f"exp_pol{c['nbkg'] - 1}", dict(nbkg=c["nbkg"] - 1)), (f"exp_pol{c['nbkg'] + 1}", dict(nbkg=c["nbkg"] + 1)),
                   (f"pol{c['nbkg']}", dict(nbkg=c["nbkg"], kind="pol")), ("shape_fix", dict(nbkg=c["nbkg"], fix_shape=True))):
        fv = fit_cats(H, mc.sigT, mc.refT, mc.frefl, edges, c["m0"], **kw)
        rv = r_of(fv["N"], fv.get("cov_N"), mc.e) if fv["valid"] else None
        var[nm] = {"R": rv["R"] if rv else None, "chi2": fv["chi2"], "ndf": fv["ndf"], "valid": fv["valid"]}
    alo, ahi = max(c["alt"][0], mc.mc_range[0]), min(c["alt"][1], mc.mc_range[1])
    if (alo, ahi) != (lo, hi) and len(m) and m.min() <= alo + 2 * BW and m.max() >= ahi - 2 * BW:
        e2 = np.arange(alo, ahi + BW / 2, BW)
        H2, _, _ = hist_cats(m, cat, e2)
        fv = fit_cats(H2, mc.sigT, mc.refT, mc.frefl, e2, c["m0"], c["nbkg"])
        rv = r_of(fv["N"], fv.get("cov_N"), mc.e) if fv["valid"] else None
        var[f"range_{alo:.2f}_{ahi:.2f}"] = {"R": rv["R"] if rv else None, "chi2": fv["chi2"], "ndf": fv["ndf"], "valid": fv["valid"]}
    o["var"] = var
    dv = [abs(v["R"] - o["R"]) for k, v in var.items() if v["R"] is not None and k != "shape_fix"]
    o["syst_bkg"] = float(max(dv)) if dv else None
    o["syst_shape"] = abs(var["shape_fix"]["R"] - o["R"]) if var["shape_fix"]["R"] is not None else None
    o["syst_mcstat"] = o["R"] * mc.eps_ratio_relerr
    o["R_err_tot"] = float(math.sqrt(o["R_err"] ** 2 + (o["syst_bkg"] or 0) ** 2 + (o["syst_shape"] or 0) ** 2 + o["syst_mcstat"] ** 2))
    o["n_var_ok"] = int(len(dv))
    o["determined"] = bool(f["valid"] and not f["at_limit"] and len(dv) >= 3)
    return o


# ----------------------------------------------------------------------------------------------------------- tryb r
def run_r(a):
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    rng = np.random.default_rng(4)
    res = {"generated": datetime.datetime.now().isoformat(timespec="seconds"), "nboot": a.nboot, "app": a.app, "force_data": bool(a.force_data),
           "cfg": {"SP": SP, "BW": BW, "WPS": WPS, "closure_gate": [CLOS_NSIG, CLOS_MAXERR, CLOS_MAXBIAS]}, "closure": {}, "data": {}, "blocked": [], "tags": {}}
    for sp in a.species.split(","):
        F = load_app(a.app, "fd_nd", sp)
        M = load_app(a.app, "mb_nd", sp)
        if F is None or M is None:
            log(f"{sp}: brak fd_nd/mb_nd w {a.app} — pomijam gatunek")
            continue
        ref = REF[sp]
        tags = [t for t in F["tags"] if t != ref and (not a.tags or re.search(a.tags, t))]
        if ref not in F["tags"] or not tags:
            log(f"{sp}: brak odniesienia {ref} albo sieci (są: {F['tags']})")
            continue
        res["tags"][sp] = {"ref": ref, "nets": tags, "n_fd": int(len(F["m"])), "n_mb": int(len(M["m"])), "files": {"fd_nd": F["nfiles"], "mb_nd": M["nfiles"]}}
        mcr = ext_range(sp)
        MCs, gate = {}, {}
        # 1) domknięcie na MC
        C = res["closure"].setdefault(sp, {})
        Mok = M["ok"] & (M["truth"] < 100)
        for t in tags:
            for b in BINS:
                for wp in WPS:
                    mc = MC(sp, F, t, ref, b, wp, mcr)
                    key = (t, bname(b), wp)
                    if not mc.ok:
                        C.setdefault(t, {}).setdefault(bname(b), {})[str(wp)] = {"pass": False, "why": "za mało sygnału MC"}
                        gate[key] = False
                        continue
                    MCs[key] = mc
                    inb = Mok & (M["pt"] >= b[0]) & (M["pt"] < b[1])
                    cat = mc.cat(M["lg"][inb, M["tags"].index(t)], M["lg"][inb, M["tags"].index(ref)])
                    o = measure(sp, mc, M["m"][inb].astype(np.float64), cat, M["uid"][inb], a.nboot if wp == WPS[0] else 0, rng, full=(wp == WPS[0]))
                    tr = M["truth"][inb]
                    ts = np.isin(tr, (0, 1)) & ((M["swap"][inb] != 1) if sp == "D0" else True)
                    Tn = [int((ts & (cat == j)).sum()) for j in range(3)]
                    rt = r_of(Tn, None, mc.e)
                    o["truth_N"] = Tn
                    if rt is not None:
                        o["R_truth"] = rt["R"]
                        o["R_truth_err"] = rt["R"] * math.sqrt(Tn[1] / max(Tn[0] + Tn[1], 1) ** 2 + Tn[2] / max(Tn[0] + Tn[2], 1) ** 2 + mc.eps_ratio_relerr ** 2)
                    if "R" in o and rt is not None:
                        o["bias"] = o["R"] - rt["R"]
                        o["pass"] = bool(o["determined"] or (wp != WPS[0] and o["fit"]["valid"])) and abs(o["bias"]) <= CLOS_NSIG * o["R_err"] and abs(o["bias"]) <= CLOS_MAXBIAS and o["R_err"] <= CLOS_MAXERR
                        o["N_fit_over_truth"] = [o["fit"]["N"][j] / Tn[j] if Tn[j] > 0 else None for j in range(3)]
                    else:
                        o["pass"] = False
                    gate[key] = bool(o["pass"])
                    C.setdefault(t, {}).setdefault(bname(b), {})[str(wp)] = o
                    log(f"domknięcie {t} {bname(b)} WP{int(wp * 100)}: R_fit {o.get('R', float('nan')):.4f} ± {o.get('R_err', float('nan')):.4f}, "
                        f"R_prawda {o.get('R_truth', float('nan')):.4f} → {'OK' if o['pass'] else 'NIE'}")
        del M
        # 2) dane — tylko po zaliczonym domknięciu
        D = res["data"].setdefault(sp, {})
        for sample in ("data_427929", "data_428770"):
            T = load_app(a.app, sample, sp, need=("uid", "m", "pt", "ok", "cls", "lg"))
            if T is None:
                log(f"{sp} {sample}: brak plików app — pomijam")
                continue
            if T["tags"] != F["tags"]:
                raise RuntimeError(f"{sp} {sample}: modele {T['tags']} ≠ MC {F['tags']}")
            log(f"{sp} {sample}: {len(T['m'])} wierszy, {T['nfiles']} plików")
            res["tags"][sp]["files"][sample] = T["nfiles"]
            for cn, cm in classes(sample, T["cls"]).items():
                for t in tags:
                    for b in BINS:
                        for wp in WPS:
                            key = (t, bname(b), wp)
                            kk = f"{sample[5:]}|{cn}"
                            if key not in MCs or (not gate.get(key) and not a.force_data):
                                res["blocked"].append(f"{sp}|{t}|{kk}|{bname(b)}|{wp}")
                                continue
                            mc = MCs[key]
                            inb = cm & T["ok"] & (T["pt"] >= b[0]) & (T["pt"] < b[1])
                            cat = mc.cat(T["lg"][inb, T["tags"].index(t)], T["lg"][inb, T["tags"].index(ref)])
                            o = measure(sp, mc, T["m"][inb].astype(np.float64), cat, T["uid"][inb], a.nboot if wp == WPS[0] else 0, rng, full=(wp == WPS[0]))
                            o["n_cand_bin"] = int(inb.sum())
                            D.setdefault(t, {}).setdefault(kk, {}).setdefault(bname(b), {})[str(wp)] = o
                            log(f"dane {t} {kk} {bname(b)} WP{int(wp * 100)}: R {o.get('R', float('nan')):.4f} ± {o.get('R_err', float('nan')):.4f} "
                                f"(tło {o.get('syst_bkg')}, kształt {o.get('syst_shape')}), G_fit {o.get('G_fit')}")
            del T
        json.dump(res, open(a.out, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    json.dump(res, open(a.out, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log(f"zapisano {a.out}; zablokowanych przez domknięcie: {len(res['blocked'])}")
    return 0


# ----------------------------------------------------------------------------------------------------------- tryb d23 (bin 2–3 GeV)
def run_d23(a):
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    rng = np.random.default_rng(23)
    R = json.load(open(a.r)) if a.r and os.path.isfile(a.r) else {"data": {}}
    res = {"generated": datetime.datetime.now().isoformat(timespec="seconds"), "nboot": a.nboot, "app": a.app, "species": {}}
    b = (2, 3)
    CLS = {"data_427929": "MB427929", "data_428770": "all"}
    for sp in a.species.split(","):
        c = SP[sp]
        F = load_app(a.app, "fd_nd", sp)
        if F is None:
            log(f"{sp}: brak fd_nd (2–3 GeV) — pomijam")
            continue
        lo, hi = c["fit"]
        edges = np.arange(lo, hi + BW / 2, BW)
        ctr = 0.5 * (edges[1:] + edges[:-1])
        mcr = ext_range(sp)
        inb = (F["pt"] >= b[0]) & (F["pt"] < b[1]) & F["ok"] & np.isin(F["truth"], (0, 1))
        sig = inb & ((F["swap"] != 1) if sp == "D0" else True)
        rf = inb & (F["swap"] == 1) if c["refl"] else np.zeros(len(sig), bool)
        w = F["w"].astype(np.float64)
        sw = sig & (np.abs(F["m"] - c["m0"]) < c["win"])
        O = res["species"].setdefault(sp, {"n_sig_mc": int(sig.sum()), "tags": F["tags"], "runs": {}})
        if sw.sum() < 500:
            O["why"] = "za mało sygnału MC 2–3 GeV"
            continue
        ref = REF[sp]
        nets = [t for t in F["tags"] if t != ref]
        net1 = next((t for t in nets if t.endswith(".mc.d5.s1") or ".mc." in t), nets[0] if nets else None)
        use_tags = [t for t in (ref, net1) if t in F["tags"]]
        mcw = {}
        for t in use_tags:
            j = F["tags"].index(t)
            for wp in D23_WPS:
                th = float(np.quantile(F["lg"][sw, j].astype(np.float64), 1 - wp))
                ps, pr = sig & (F["lg"][:, j] >= th), rf & (F["lg"][:, j] >= th)
                if ps.sum() < 300:
                    continue
                sT = Tmpl(F["m"][ps], w[ps], *mcr)
                mcw[(t, wp)] = {"thr": th, "eps": float(w[ps].sum() / w[sig].sum()), "sigT": sT, "n": int(ps.sum()),
                                "refT": Tmpl(F["m"][pr], w[pr], *mcr) if pr.sum() >= 300 else None, "frefl": float(w[pr].sum() / w[ps].sum()) if pr.sum() >= 300 else 0.0,
                                "fwin": float(sT.bins(np.array([c["m0"] - c["win"], c["m0"] + c["win"]]), c["m0"])[0])}
        for sample, cn in CLS.items():
            T = load_app(a.app, sample, sp, need=("uid", "m", "pt", "ok", "cls", "lg"))
            if T is None:
                continue
            if T["tags"] != F["tags"]:
                raise RuntimeError(f"{sp} {sample}: modele {T['tags']} ≠ MC {F['tags']}")
            cm = classes(sample, T["cls"])[cn] & T["ok"] & (T["pt"] >= b[0]) & (T["pt"] < b[1])
            run = sample[5:]
            # oczekiwany sygnał: 0,6 × uzysk drzewa w 3–5 GeV (ten sam run i klasa) z dnet4_R.json
            y35 = None
            for t_, v in (R.get("data", {}).get(sp, {}) or {}).items():
                x = v.get(f"{run}|{cn}", {}).get("3-5", {}).get("0.7", {})
                if x.get("determined") and x.get("y_tree"):
                    y35 = x["y_tree"]
                    break
            blind = np.abs(ctr - c["m0"]) < c["win"] + 0.015 - 1e-9
            OR = O["runs"].setdefault(f"{run}|{cn}", {"n_cand": int(cm.sum()), "y_tree_3_5": y35, "nfiles": T["nfiles"], "nfiles_fd": F["nfiles"], "tags": {}})
            log(f"{sp} {run}|{cn}: {int(cm.sum())} kandydatów 2–3 GeV, y(3–5) = {y35}")
            for t in use_tags:
                j = T["tags"].index(t)
                scan = {}
                for wp in D23_WPS:
                    if (t, wp) not in mcw:
                        continue
                    q = mcw[(t, wp)]
                    sel = cm & (T["lg"][:, j] >= q["thr"])
                    H = np.histogram(T["m"][sel], edges)[0].astype(np.float64)[None, :]
                    fb = fit_cats(H, [q["sigT"]], [q["refT"]], [q["frefl"]], edges, c["m0"], c["nbkg"], blind=blind, win=c["win"])
                    B = fb["bkg_win"][0]
                    S = D23_SCALE * y35 * q["eps"] * q["fwin"] if y35 else None
                    scan[wp] = {"thr": q["thr"], "eps": q["eps"], "B_win": B, "S_exp": S, "valid_sb": fb["valid"],
                                "Z_exp": (S / math.sqrt(S + B) if S and S + B > 0 else None), "fom_eps_sqrtB": q["eps"] / math.sqrt(max(B, 1.0)), "n_sel": int(sel.sum())}
                # wybór WP (bez okna sygnału): max Z_exp; bez y(3–5) → max ε/√B; remis < 2 % → luźniejszy
                fom = {wp: (s["Z_exp"] if y35 else s["fom_eps_sqrtB"]) for wp, s in scan.items() if s["valid_sb"]}
                if not fom:
                    OR["tags"][t] = {"scan": {str(k): v for k, v in scan.items()}, "why": "fit pasm nieudany"}
                    continue
                best = max(fom.values())
                wp0 = max(wp for wp, z in fom.items() if z >= 0.98 * best)
                out = {"scan": {str(k): v for k, v in scan.items()}, "wp": wp0, "rule": "Z_exp" if y35 else "eps/sqrtB", "fits": {}}
                for wp in [x for x in (round(wp0 + 0.1, 1), wp0, round(wp0 - 0.1, 1)) if (t, x) in mcw]:
                    q = mcw[(t, wp)]
                    sel = cm & (T["lg"][:, j] >= q["thr"])
                    msel, usel = T["m"][sel].astype(np.float64), T["uid"][sel]
                    H = np.histogram(msel, edges)[0].astype(np.float64)[None, :]
                    fo = {}
                    for nm, kw in (("nom", dict(nbkg=c["nbkg"], fix_shape=True)), (f"exp_pol{c['nbkg'] - 1}", dict(nbkg=c["nbkg"] - 1, fix_shape=True)),
                                   (f"exp_pol{c['nbkg'] + 1}", dict(nbkg=c["nbkg"] + 1, fix_shape=True)), ("free_shape", dict(nbkg=c["nbkg"]))):
                        f = fit_cats(H, [q["sigT"]], [q["refT"]], [q["frefl"]], edges, c["m0"], win=c["win"], **kw)
                        e_ = math.sqrt(max(f["cov_N"][0][0], 0)) if f.get("cov_N") else None
                        fo[nm] = {"N": f["N"][0], "N_err": e_, "Z": (f["N"][0] / e_ if e_ else None), "valid": f["valid"], "chi2": f["chi2"], "ndf": f["ndf"],
                                  "dmu": f["dmu"], "sc": f["sc"], "dmu_err": f.get("dmu_err"), "sc_err": f.get("sc_err"), "B_win": f["bkg_win"][0], "x": f["x"]}
                    if wp == wp0 and a.nboot > 0 and fo["nom"]["valid"]:
                        u, inv = np.unique(usel, return_inverse=True)
                        ib = np.floor((msel - edges[0]) / BW + 1e-9).astype(np.int64)
                        okb = (ib >= 0) & (ib < len(ctr))
                        Nb = []
                        for _ in range(a.nboot):
                            wv = rng.poisson(1.0, len(u)).astype(np.float64)[inv]
                            Hb = np.bincount(ib[okb], weights=wv[okb], minlength=len(ctr))[None, :]
                            fbt = fit_cats(Hb, [q["sigT"]], [q["refT"]], [q["frefl"]], edges, c["m0"], c["nbkg"], fix_shape=True, start=fo["nom"]["x"])
                            if fbt["valid"]:
                                Nb.append(fbt["N"][0])
                        if len(Nb) > 5:
                            rms = float(np.std(Nb))
                            fo["nom"]["boot_rms"], fo["nom"]["boot_ok"] = rms, len(Nb)
                            fo["nom"]["N_err"] = max(fo["nom"]["N_err"] or 0, rms)
                            fo["nom"]["Z"] = fo["nom"]["N"] / fo["nom"]["N_err"]
                    for v_ in fo.values():
                        v_.pop("x", None)
                    fo["eps"], fo["y"] = q["eps"], fo["nom"]["N"] / q["eps"]
                    fo["S_over_B"] = fo["nom"]["N"] * q["fwin"] / max(fo["nom"]["B_win"], 1e-9)
                    out["fits"][str(wp)] = fo
                OR["tags"][t] = out
                z = out["fits"].get(str(wp0), {}).get("nom", {})
                log(f"  {t}: WP{int(wp0 * 100)} ({out['rule']}), N = {z.get('N')} ± {z.get('N_err')}, Z = {z.get('Z')}")
            del T
    json.dump(res, open(a.out, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log(f"zapisano {a.out}")
    return 0


# ----------------------------------------------------------------------------------------------------------- raport
def pm(x, e=None, nd=3):
    if x is None or not np.isfinite(x):
        return "–"
    return f"{x:.{nd}f}" + (f" ± {e:.{nd}f}" if e is not None and np.isfinite(e) else "")


def report(a):
    R = json.load(open(a.res))
    L = [f"Wygenerowano {R['generated']} (`dn4_ana.py report`, {a.res}); bootstrap {R['nboot']} replik; bramka domknięcia: |R_fit − R_prawda| ≤ "
         f"{R['cfg']['closure_gate'][0]}σ i ≤ {R['cfg']['closure_gate'][2]}, σ ≤ {R['cfg']['closure_gate'][1]}.", ""]
    for sp, C in R["closure"].items():
        L += [f"### {sp} — domknięcie na MC (pseudo-dane mb_nd, ε i szablony z fd_nd), WP70", "",
              "| model | bin | R_fit | R_prawda (ε z fd_nd) | obciążenie | N_fit/N_prawda (b / n / t) | syst. tła | χ²/ndf | bramka |", "|---|---|---|---|---|---|---|---|---|"]
        for t, v in C.items():
            for b, w in v.items():
                o = w.get("0.7", {})
                nf = o.get("N_fit_over_truth") or [None] * 3
                L.append(f"| `{t}` | {b} | {pm(o.get('R'), o.get('R_err'))} | {pm(o.get('R_truth'), o.get('R_truth_err'))} | {pm(o.get('bias'))} | "
                         + " / ".join(pm(x, nd=2) for x in nf) + f" | {pm(o.get('syst_bkg'))} | {pm(o.get('fit', {}).get('chi2'), nd=0)}/{o.get('fit', {}).get('ndf', '–')} | "
                         f"{'OK' if o.get('pass') else 'NIE'} |")
        L.append("")
    for sp, D in R["data"].items():
        L += [f"### {sp} — dane: R = [N/ε]_sieć / [N/ε]_drzewo, WP70", "",
              "| model | run:klasa | bin | R ± stat | syst. tło / kształt / MC | σ całk. | WP50 | WP80 | G_fit | σ_N/N sieć / drzewo | boot/HESSE | określony |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for t, v in D.items():
            for k, w in v.items():
                for b, x in w.items():
                    o = x.get("0.7", {})
                    bt = o.get("boot", {}).get("rms")
                    L.append(f"| `{t}` | {k.replace('|', ':')} | {b} | {pm(o.get('R'), o.get('R_err'))} | {pm(o.get('syst_bkg'))} / {pm(o.get('syst_shape'))} / {pm(o.get('syst_mcstat'))} | "
                             f"{pm(o.get('R_err_tot'))} | {pm(x.get('0.5', {}).get('R'), x.get('0.5', {}).get('R_err'))} | {pm(x.get('0.8', {}).get('R'), x.get('0.8', {}).get('R_err'))} | "
                             f"{pm(o.get('G_fit'), nd=2)} | {pm(o.get('relerr_net'))} / {pm(o.get('relerr_tree'))} | "
                             f"{pm(bt / o['R_err_hesse'] if bt and o.get('R_err_hesse') else None, nd=2)} | {'tak' if o.get('determined') else 'NIE'} |")
        L.append("")
    if R.get("blocked"):
        L += [f"Zablokowane przez domknięcie (brak pomiaru na danych): {len(R['blocked'])} — " + "; ".join(R["blocked"][:12]) + (" …" if len(R["blocked"]) > 12 else ""), ""]
    if a.d23 and os.path.isfile(a.d23):
        D = json.load(open(a.d23))
        L += ["### Bin 2–3 GeV (WP wybrany z MC + pasm, okno zaślepione przy wyborze)", "",
              "| gatunek | run:klasa | klasyfikator | WP | ε_S MC | N_S | Z = N/σ | Z (tło rząd −1 / +1) | S/B | Z_oczek |", "|---|---|---|---|---|---|---|---|---|---|"]
        for sp, O in D["species"].items():
            for k, v in O.get("runs", {}).items():
                for t, o in v["tags"].items():
                    if "wp" not in o:
                        L.append(f"| {sp} | {k.replace('|', ':')} | `{t}` | – | – | – | – | – | – | {o.get('why')} |")
                        continue
                    f = o["fits"].get(str(o["wp"]), {})
                    alt = [f.get(n, {}).get("Z") for n in f if n.startswith("exp_pol")]
                    L.append(f"| {sp} | {k.replace('|', ':')} | `{t}` | {int(o['wp'] * 100)} | {pm(f.get('eps'))} | {pm(f.get('nom', {}).get('N'), f.get('nom', {}).get('N_err'), 0)} | "
                             f"{pm(f.get('nom', {}).get('Z'), nd=1)} | " + " / ".join(pm(x, nd=1) for x in alt) + f" | {pm(f.get('S_over_B'))} | {pm(o['scan'][str(o['wp'])].get('Z_exp'), nd=1)} |")
        L.append("")
    print("\n".join(L))
    return 0


# ----------------------------------------------------------------------------------------------------------- zabawki
def _toy_rows(rng, sp, n_sig, n_bkg, kind, shift_net=0.0, pt=(3.0, 30.0), mr=None):
    """Wiersze jak z dn4_apply: 2 modele (drzewo, sieć), sygnał Gauss (+ odbicie D0), tło exp. shift_net: przesunięcie logitu sieci sygnału
    (emulacja ε_dane ≠ ε_MC). Zwraca (tabela, prawda_sygnału)."""
    c = SP[sp]
    lo, hi = mr or c["alt"] if sp == "D0" else (mr or c["fit"])
    nr = n_sig if c["refl"] else 0
    n = n_sig + nr + n_bkg
    is_s = np.arange(n) < n_sig
    is_r = (np.arange(n) >= n_sig) & (np.arange(n) < n_sig + nr)
    m = np.empty(n)
    m[is_s] = c["m0"] + np.where(rng.random(n_sig) < 0.8, rng.normal(0, 0.016, n_sig), rng.normal(0, 0.035, n_sig))
    if nr:
        m[is_r] = c["m0"] + rng.normal(0.0, 0.09, nr)
    u = rng.random(n_bkg)
    k = 3.0
    m[~(is_s | is_r)] = lo - np.log(1 - u * (1 - np.exp(-k * (hi - lo)))) / k
    z1, z2 = rng.normal(size=n), rng.normal(size=n)
    lt = np.where(is_s, 2.0, np.where(is_r, 0.8, 0.0)) + z1
    ln = np.where(is_s, 3.0 + shift_net, np.where(is_r, 1.0, 0.0)) + 0.6 * z1 + 0.8 * z2 * 1.2
    ptv = np.clip(pt[0] + rng.exponential(2.5, n), pt[0], pt[1] - 1e-3)
    keep = (m >= lo) & (m < hi)
    ev = rng.integers(1, max(n // 2, 2), n).astype(np.uint64) * np.uint64(7919)          # ~2 kandydatów na zdarzenie
    truth = np.where(is_s | is_r, 0, 8).astype(np.int16) if kind != "data" else np.full(n, -1, np.int16)
    swap = np.where(is_r, 1, 0).astype(np.int8)
    if kind == "fd":
        keep &= is_s | is_r
    T = {"uid": ev[keep], "m": m[keep].astype(np.float32), "pt": ptv[keep].astype(np.float32), "ok": np.ones(keep.sum(), bool),
         "cls": np.where(rng.random(keep.sum()) < 0.5, 1 + 2 + 8, 1 + 4).astype(np.uint8), "truth": truth[keep], "swap": swap[keep],
         "w": np.ones(keep.sum(), np.float32), "lg": np.stack([lt[keep], ln[keep]], 1).astype(np.float16)}
    return T, is_s[keep]


def _toy_write(d, sample, parts):
    os.makedirs(d, exist_ok=True)
    out = {}
    tags = {}
    for sp, T in parts.items():
        for k, v in T.items():
            out[f"{sp}__{k}"] = v
        tags[sp] = [REF[sp], f"{sp}.tr.H+e04.mc.d5.s1"]
    out["meta"] = np.array(json.dumps({"sample": sample, "tags": tags}))
    np.savez(os.path.join(d, f"app_{sample}_0.npz"), **out)


def selftest(a):
    import tempfile
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    ctx = None
    d = a.dir
    if not d:
        ctx = tempfile.TemporaryDirectory(prefix="dn4_ana_"); d = ctx.name
    try:
        rng = np.random.default_rng(11)
        d150, d23 = os.path.join(d, "app", "d150"), os.path.join(d, "app", "d23")
        truth = {}
        SH = {"data_427929": -0.15, "data_428770": -0.5}                       # przesunięcie logitu sieci na sygnale „danych” → R < 1
        for sample, kind, ns, nbk in (("fd_nd", "fd", 150_000, 0), ("mb_nd", "mb", 60_000, 1_200_000), ("data_427929", "data", 40_000, 1_500_000),
                                       ("data_428770", "data", 60_000, 4_000_000)):
            parts = {}
            for sp in ("D0", "Dplus"):
                T, s = _toy_rows(rng, sp, ns, nbk, kind, SH.get(sample, 0.0))
                parts[sp] = T; truth[(sample, sp)] = (T, s)
            _toy_write(d150, sample, parts)
        ns = argparse.Namespace(app=d150, out=os.path.join(d, "R.json"), nboot=a.nboot, species="D0,Dplus", tags="", force_data=False)
        run_r(ns)
        R = json.load(open(ns.out))
        bad, npull, pulls = [], 0, []
        for sp in ("D0", "Dplus"):
            t = f"{sp}.tr.H+e04.mc.d5.s1"
            F = load_app(d150, "fd_nd", sp)
            for b in BINS:
                o = R["closure"][sp][t][bname(b)]["0.7"]
                if not o.get("pass"):
                    bad.append(f"domknięcie {sp} {bname(b)}: {o.get('R')} ± {o.get('R_err')} vs {o.get('R_truth')}")
                pulls.append((o["R"] - o["R_truth"]) / o["R_err"])
                mc = MC(sp, F, t, REF[sp], b, 0.7)
                for sample, kk in (("data_427929", "427929|MB427929"), ("data_428770", "428770|all")):
                    T, s = truth[(sample, sp)]
                    cm = classes(sample, T["cls"])[kk.split("|")[1]]
                    inb = cm & (T["pt"] >= b[0]) & (T["pt"] < b[1]) & s
                    cat = mc.cat(T["lg"][inb, 1], T["lg"][inb, 0])
                    Rt = r_of([int((cat == j).sum()) for j in range(3)], None, mc.e)["R"]
                    x = R["data"][sp][t][kk][bname(b)]["0.7"]
                    pull = (x["R"] - Rt) / x["R_err"]
                    pulls.append(pull); npull += 1
                    log(f"zabawka {sp} {kk} {bname(b)}: R_fit {x['R']:.4f} ± {x['R_err']:.4f} (HESSE {x['R_err_hesse']:.4f}, boot {x.get('boot', {}).get('rms')}), "
                        f"R_prawda {Rt:.4f}, pull {pull:+.2f}, syst tło {x['syst_bkg']:.4f}, G_fit {x['G_fit']:.2f}")
                    if abs(pull) > 3.5 or not x["determined"]:
                        bad.append(f"dane {sp} {kk} {bname(b)}: pull {pull:.2f}, determined {x['determined']}")
        pulls = np.array(pulls)
        log(f"pulle R (n = {len(pulls)}): średnia {pulls.mean():+.2f}, rms {pulls.std():.2f}, max |pull| {np.abs(pulls).max():.2f}; zablokowane: {len(R['blocked'])}")
        if abs(pulls.mean()) > 3 / math.sqrt(len(pulls)) or pulls.std() > 1.8:
            bad.append(f"rozkład pulli: średnia {pulls.mean():.2f}, rms {pulls.std():.2f}")
        # bin 2–3 GeV
        for sample, kind, ns_, nbk in (("fd_nd", "fd", 40_000, 0), ("data_427929", "data", 3_000, 1_500_000), ("data_428770", "data", 4_000, 3_000_000)):
            parts = {}
            for sp in ("D0", "Dplus"):
                parts[sp], _ = _toy_rows(rng, sp, ns_, nbk, kind, 0.0, pt=(2.0, 3.0))
            _toy_write(d23, sample, parts)
        n23 = argparse.Namespace(app=d23, r=ns.out, out=os.path.join(d, "d23.json"), nboot=min(a.nboot, 20), species="D0,Dplus")
        run_d23(n23)
        D = json.load(open(n23.out))
        for sp in ("D0", "Dplus"):
            for k, v in D["species"][sp]["runs"].items():
                for t, o in v["tags"].items():
                    f = o["fits"][str(o["wp"])]
                    ntrue = (3_000 if "427929" in k else 4_000) * (0.5 if "MB" in k else 1.0)
                    pull = (f["y"] - ntrue) / (f["nom"]["N_err"] / f["eps"])
                    log(f"zabawka 2–3 {sp} {k} {t}: WP{int(o['wp'] * 100)}, y = N/ε = {f['y']:.0f} (prawda ≈ {ntrue:.0f}), Z {f['nom']['Z']:.1f}, pull {pull:+.2f}")
                    if abs(pull) > 4:
                        bad.append(f"2–3 {sp} {k} {t}: pull {pull:.2f}")
        report(argparse.Namespace(res=ns.out, d23=n23.out))
        if bad:
            log("BŁĄD selftest: " + "; ".join(bad))
            return 1
        log("selftest OK: domknięcie R na zabawkach MC, R na „danych” z przesuniętą efektywnością sieci zgodny z prawdą, bin 2–3 działa")
    finally:
        if ctx is not None:
            ctx.cleanup()
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["r", "d23", "report", "selftest"])
    ap.add_argument("--app"); ap.add_argument("--out"); ap.add_argument("--r", default=""); ap.add_argument("--res"); ap.add_argument("--d23", default="")
    ap.add_argument("--nboot", type=int, default=int(os.environ.get("DN4_NBOOT", "100")))
    ap.add_argument("--species", default="D0,Dplus"); ap.add_argument("--tags", default="")
    ap.add_argument("--force-data", action="store_true"); ap.add_argument("--dir", default="")
    a = ap.parse_args(argv)
    return {"r": run_r, "d23": run_d23, "report": report, "selftest": selftest}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())

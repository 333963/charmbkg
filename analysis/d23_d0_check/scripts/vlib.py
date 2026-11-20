"""Wspolne: wczytanie zredukowanych plikow, progi WP z MC (jak dn4_ana.run_d23), fit (fit_cats/Tmpl z dn4_ana — import tylko do odczytu)."""
import os, sys, json, math
import numpy as np
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(BASE, "..", "Q_lxplus", "dnet4"))
import dn4_ana as A  # noqa
M0, WIN = 1.86484, 0.055
TAGS = ["D0.xgb.H.sb.n.s1", "D0.tr.H+e04.mc.d5.s1"]   # kolumny lg: 0 drzewo, 1 siec
NM = ["drzewo", "siec"]
MCR = (1.70, 2.05)

def load(name):
    z = np.load(os.path.join(BASE, "in", name + ".npz"))
    return {k: z[k] for k in z.files}

def mc_wp(F, j, wp):
    """prog = kwantyl sygnalu MC (poprawna hipoteza) w oknie; eps = ulamek calego sygnalu (jak etap 4)."""
    sig = F["cat"] == 0
    sw = sig & (np.abs(F["m"] - M0) < WIN)
    th = float(np.quantile(F["lg"][sw, j].astype(np.float16).astype(np.float64), 1 - wp))
    return th

def templ(F, j, lo, hi=np.inf):
    """szablony sygnalu/odbicia i eps dla wycinka score [lo, hi)."""
    w = F["w"].astype(np.float64)
    sig, rf = F["cat"] == 0, F["cat"] == 1
    s = (F["lg"][:, j] >= lo) & (F["lg"][:, j] < hi)
    ps, pr = sig & s, rf & s
    sT = A.Tmpl(F["m"][ps], w[ps], *MCR)
    rT = A.Tmpl(F["m"][pr], w[pr], *MCR) if pr.sum() >= 300 else None
    return dict(sigT=sT, refT=rT, frefl=float(w[pr].sum() / w[ps].sum()) if rT is not None else 0.0, eps=float(w[ps].sum() / w[sig].sum()),
                eps_err=math.sqrt(ps.sum()) / sig.sum(), n=int(ps.sum()),
                fwin=float(sT.bins(np.array([M0 - WIN, M0 + WIN]), M0)[0]))

def fit(m, q, rng=(1.75, 2.00), nbkg=2, kind="exp", free=False, frefl=None, w=None, start=None):
    edges = np.arange(rng[0], rng[1] + A.BW / 2, A.BW)
    H = np.histogram(m, edges, weights=w)[0].astype(np.float64)[None, :]
    fr = q["frefl"] if frefl is None else frefl
    f = A.fit_cats(H, [q["sigT"]], [q["refT"]], [fr], edges, M0, nbkg, kind=kind, fix_shape=not free, start=start)
    e = math.sqrt(max(f["cov_N"][0][0], 0)) if f.get("cov_N") else float("nan")
    return dict(N=f["N"][0], err=e, Z=f["N"][0] / e if e else float("nan"), chi2=f["chi2"], ndf=f["ndf"], dmu=f["dmu"], sc=f["sc"], dmu_err=f.get("dmu_err"), sc_err=f.get("sc_err"),
                valid=f["valid"], B=f["bkg_win"][0], x=f["x"], fval=f["fval"], H=H[0], edges=edges, at_limit=f["at_limit"])

def boot(m, uid, q, nb=100, seed=1, **kw):
    rng = np.random.default_rng(seed)
    u, inv = np.unique(uid, return_inverse=True)
    f0 = fit(m, q, **kw)
    out = []
    for _ in range(nb):
        wv = rng.poisson(1.0, len(u)).astype(np.float64)[inv]
        f = fit(m, q, w=wv, start=f0["x"], **kw)
        if f["valid"]:
            out.append(f["N"])
    return float(np.std(out)), len(out)


def fit2(H, edges, q, nbkg=2, kind="exp", free_shape=False, frefl="fix", fixN=None, gauss=None):
    """Niezalezna implementacja fitu (binned Poisson): N*[sig(dmu,sc) + frefl*refl] + NB*tlo. frefl: 'fix' | 'free' | liczba.
    gauss=(mu0,sig0): sygnal = Gauss z wolnym mu, sigma zamiast szablonu. fixN: ustalone N (np. 0 dla testu LR)."""
    from iminuit import Minuit
    from scipy.stats import norm
    ctr = 0.5 * (edges[1:] + edges[:-1]); x = 2 * (ctr - edges[0]) / (edges[-1] - edges[0]) - 1
    refl = q["refT"].bins(edges, M0) if q["refT"] is not None else np.zeros(len(ctr))
    fr0 = q["frefl"] if frefl in ("fix", "free") else float(frefl)

    def comps(p):
        N, NB, a, b, fr = p[0], p[1], p[2], p[3], p[4]
        s = np.diff(norm.cdf(edges, a, b)) if gauss else q["sigT"].bins(edges, M0, a, b)
        sh = A.bshape(x, p[5:], kind)
        return N * s, N * fr * refl, NB * sh / sh.sum()

    def cost(p):
        mu = np.maximum(sum(comps(p)), 1e-9)
        return 2.0 * float(np.sum(mu - H) + np.sum(np.where(H > 0, H * np.log(np.where(H > 0, H, 1.0) / mu), 0.0)))
    inw = np.abs(ctr - M0) < 0.03
    s0 = max(H[inw].sum() - H[~inw].sum() * inw.sum() / (~inw).sum(), 0) * 1.2
    names = ["N", "NB", "a", "b", "fr"] + [f"c{k}" for k in range(nbkg)]
    p0 = [s0 if fixN is None else fixN, H.sum(), gauss[0] if gauss else 0.0, gauss[1] if gauss else 1.0, fr0] + [0.0] * nbkg
    m = Minuit(cost, p0, name=names); m.errordef = 1.0
    m.limits["NB"] = (0, None); m.limits["fr"] = (0, 5)
    m.limits["a"] = (1.84, 1.89) if gauss else (-0.012, 0.012); m.limits["b"] = (0.008, 0.05) if gauss else (0.6, 1.6)
    m.fixed["fr"] = frefl != "free"
    m.fixed["a"] = True; m.fixed["b"] = True
    if fixN is not None:
        m.fixed["N"] = True
    m.migrad(ncall=20000)
    if (free_shape or gauss) and fixN is None:
        m.fixed["a"] = False; m.fixed["b"] = False
        m.migrad(ncall=20000)
    if not m.valid:
        m.simplex(); m.migrad(ncall=40000)
    m.hesse()
    v = np.array(m.values); e = np.array(m.errors)
    c = comps(v)
    wm = np.abs(ctr - M0) < WIN
    return dict(N=v[0], err=e[0], a=v[2], a_err=e[2], b=v[3], b_err=e[3], fr=v[4], fr_err=e[4], fval=m.fval, valid=bool(m.valid), chi2=float(np.sum((H - sum(c)) ** 2 / np.maximum(sum(c), 1e-9))),
                ndf=len(H) - m.nfit, comps=c, B=float(c[2][wm].sum()), Rw=float(c[1][wm].sum()), Sw=float(c[0][wm].sum()), ctr=ctr, x=v)


def hist(m, rng=(1.75, 2.00), w=None):
    edges = np.arange(rng[0], rng[1] + A.BW / 2, A.BW)
    return np.histogram(m, edges, weights=w)[0].astype(np.float64), edges


def boot2(m, uid, q, nb=100, seed=3, rng=(1.75, 2.00), keys=("N",), **kw):
    r = np.random.default_rng(seed)
    u, inv = np.unique(uid, return_inverse=True)
    out = []
    for _ in range(nb):
        wv = r.poisson(1.0, len(u)).astype(np.float64)[inv]
        H, e = hist(m, rng, wv)
        f = fit2(H, e, q, **kw)
        if f["valid"]:
            out.append([f[k] for k in keys])
    out = np.array(out)
    return out.std(0), len(out)

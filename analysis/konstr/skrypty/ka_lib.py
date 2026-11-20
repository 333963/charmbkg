"""KONSTR A — pomiar odbicia K↔π z bliźniaków (m_swap) wg KRYTERIA §A."""
import numpy as np
M0 = 1.86484
BINS = ["3-5", "5-8", "8-30"]
EPS = [round(0.95 - 0.05 * i, 2) for i in range(18)]
SE = np.arange(1.60, 2.15 + 1e-9, 0.005)
SC = 0.5 * (SE[1:] + SE[:-1])
WINM = np.abs(SC - 1.865) < 0.040          # 1,825–1,905 (16 binów)
SBM = (np.abs(SC - 1.865) > 0.040) & (np.abs(SC - 1.865) < 0.080)   # 2 × 8 binów
SL = 5                                      # 25 MeV
SLC = 1.70 + 0.025 * (np.arange(14) + 0.5)
USE = np.abs(SLC - 1.865) > 0.05 + 0.0124   # przedziały rozłączne ze skrzyżowaniem |m − m0| < 50 MeV


def cum(G, e):
    return G[EPS.index(e) + 1:].sum(0).astype(float)


def raw(G2):
    """G2 [70, 110] → (W − SB) per przedział m 25 MeV, wariancja Poissona."""
    g = G2.reshape(14, SL, 110).sum(1)
    w, s = g[:, WINM].sum(1), g[:, SBM].sum(1)
    return w - s, w + s


def eff(T1):
    """ε_okno − ε_pasma z rozkładu m_swap wierszy odbicia MC (t1), per przedział m i łącznie."""
    p = T1.sum(0)
    return float((p[WINM].sum() - p[SBM].sum()) / p.sum())


def tmpl_slices(T1m, s=1.0):
    """szablon odbicia fd (projekcja m, 5 MeV, 1,70–2,05) → ułamki per przedział 25 MeV; s = skala szerokości wokół m0 (interpolacja dystrybuanty)."""
    e = np.arange(1.70, 2.05 + 1e-9, 0.005)
    c = np.concatenate([[0], np.cumsum(T1m)]) / T1m.sum()
    x = np.clip((e - M0) / s + M0, e[0], e[-1])
    cs = np.interp(x, e, c)
    return np.diff(cs[::SL]), np.diff(cs)


def measure(G2, T1fd, free_s=False):
    """→ A_R (liczba wierszy odbicia w 1,70–2,05), błąd, s_R, χ², ndf, punkty."""
    d, v = raw(G2)
    e = eff(T1fd)
    y, ey = d / e, np.sqrt(np.maximum(v, 1.0)) / e
    T1m = T1fd.sum(1)
    best = None
    for s in (np.arange(0.6, 1.61, 0.02) if free_s else [1.0]):
        f, _ = tmpl_slices(T1m, s)
        u = USE
        A = np.sum(y[u] * f[u] / ey[u] ** 2) / np.sum(f[u] ** 2 / ey[u] ** 2)
        c2 = float(np.sum(((y[u] - A * f[u]) / ey[u]) ** 2))
        if best is None or c2 < best[3]:
            best = (float(A), float(1 / np.sqrt(np.sum(f[u] ** 2 / ey[u] ** 2))), float(s), c2)
    return {"A": best[0], "eA": best[1], "s": best[2], "chi2": best[3], "ndf": int(USE.sum()) - 1 - int(free_s), "y": y, "ey": ey, "eff": e}


# ---- estymator E2 (odstępstwo, po niezaliczeniu testu MC przez E1 = okno − pasma; ustalony przed danymi):
# per przedział m: ważony fit liniowy rozkładu m_swap w |m_swap − 1,865| < 120 MeV modelem [kształt sygnału w m_swap z fd (wiersze odbicia)] + wielomian 2.
def slices2(G2, T1fd, half=0.12, order=2):
    g = G2.reshape(14, SL, 110).sum(1)
    sel = np.abs(SC - 1.865) < half
    x = (SC[sel] - 1.865) / half
    p = T1fd.sum(0)
    X = np.stack([p[sel] / p.sum()] + [x ** k for k in range(order + 1)], 1)
    y, ey = np.zeros(14), np.zeros(14)
    for i in range(14):
        n = g[i, sel]
        w = 1 / np.maximum(n, 1.0)
        for _ in range(3):
            A = X.T @ (X * w[:, None])
            c = np.linalg.solve(A, X.T @ (n * w))
            w = 1 / np.maximum(X @ c, 1.0)
        y[i], ey[i] = c[0], np.sqrt(np.linalg.inv(A)[0, 0])
    return y, ey


def measure2(G2, T1fd, s=1.0):
    y, ey = slices2(G2, T1fd)
    f, _ = tmpl_slices(T1fd.sum(1), s)
    u = USE
    A = np.sum(y[u] * f[u] / ey[u] ** 2) / np.sum(f[u] ** 2 / ey[u] ** 2)
    eA = 1 / np.sqrt(np.sum(f[u] ** 2 / ey[u] ** 2))
    c2 = float(np.sum(((y[u] - A * f[u]) / ey[u]) ** 2))
    return {"A": float(A), "eA": float(eA), "chi2": c2, "ndf": int(u.sum()) - 1, "y": y, "ey": ey, "f": f}

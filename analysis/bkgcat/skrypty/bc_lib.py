#!/usr/bin/env python3
"""BKGCAT (20.11.2026, projekt B; POST HOC względem FITBIAS) — wspólne: wczytanie bkgcat.npz (Q_lxplus/bkgcat), grupy kategorii prawdy, wybór przedziałów score.
Importuje FITBIAS/skrypty/fb_fit.py bez zmian. Tylko MC mb_nd + (w bc_m8.py) inkluzywne piki danych z WPSCAN. Bez AFP."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "..")
CW = os.path.join(BASE, "..")
sys.path.insert(0, os.path.join(CW, "FITBIAS_2026-11-19", "skrypty"))
import fb_fit as F  # noqa: E402

IN = os.path.join(BASE, "in", "bkgcat.npz")
WPF = os.path.join(CW, "WPSCAN_2026-11-19", "runs", "lxplus", "wps_hist.npz")
TAG = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
SPS = ["D0", "Dplus"]
BINS = F.BINS
GORD = ["Sig", "Refl", "DMM", "CM", "PC", "MB", "ML", "OtherB", "NoCharm"]
GNAME = {"Sig": "sygnał", "Refl": "odbicie K↔π", "DMM": "DecayModeMismatch", "CM": "CharmMismatch", "PC": "PartComb", "MB": "Comb: tory z jednego hadronu b",
         "ML": "Comb: tory z jednego lekkiego rezonansu", "OtherB": "Comb: inne z b", "NoCharm": "Comb: bez powabu (OtherNoB)"}


def group(key):
    for p, g in (("Sig", "Sig"), ("Refl", "Refl"), ("DMM", "DMM"), ("CM", "CM"), ("PC|", "PC"), ("MB|", "MB"), ("ML|", "ML"), ("OtherNoB", "NoCharm"), ("Other", "OtherB")):
        if key.startswith(p):
            return g
    return "?"


class BK:
    def __init__(self, fn=IN):
        z = np.load(fn, allow_pickle=False)
        self.z = z
        self.keys = [str(k) for k in z["keys"]]
        self.tab, self.sw = z["tab"], z["sw"]
        self.meta = json.loads(str(z["meta"]))
        self.eps = self.meta["eps_desc"]
        self._d = {}

    def dense(self, sp, bn):
        """→ {klucz: [19, nm] float} (liczności; wagi ≡ 1 sprawdzone w bc_check)."""
        k = (sp, bn)
        if k not in self._d:
            lo, hi = F.SP[sp]["mr"]
            nm = int(round((hi - lo) / 0.001))
            t = self.tab[(self.tab[:, 0] == SPS.index(sp)) & (self.tab[:, 1] == BINS.index(bn))]
            o = {}
            for kid in np.unique(t[:, 2]):
                r = t[t[:, 2] == kid]
                h = np.zeros((19, nm))
                np.add.at(h, (r[:, 3], r[:, 4]), r[:, 5])
                o[self.keys[kid]] = h
            self._d[k] = o
        return self._d[k]

    def jj(self, elo, ehi=None):
        return self.eps.index(elo), (None if ehi is None else self.eps.index(ehi))

    def cat(self, sp, bn, elo, ehi=None):
        """→ {klucz: hist 1 MeV [nm]} dla kategorii score ε_lo ≥ … > ε_hi (ehi=None: skumulowana)."""
        jlo, jhi = self.jj(elo, ehi)
        return {k: h[jlo + 1:(None if jhi is None else jhi + 1)].sum(0) for k, h in self.dense(sp, bn).items()}

    def groups(self, sp, bn, elo, ehi=None):
        c = self.cat(sp, bn, elo, ehi)
        o = {g: 0.0 for g in GORD}
        for k, h in c.items():
            o[group(k)] = o[group(k)] + h
        z0 = np.zeros(next(iter(c.values())).shape)
        return {g: (v if isinstance(v, np.ndarray) else z0.copy()) for g, v in o.items()}


def h5(h1):
    return h1.reshape(h1.shape[:-1] + (h1.shape[-1] // 5, 5)).sum(-1).astype(np.float64)


def full_edges(sp):
    lo, hi = F.SP[sp]["mr"]
    return np.arange(lo, hi + F.BW / 2, F.BW)


def cutr(h, sp, rng):
    lo = F.SP[sp]["mr"][0]
    i0, i1 = int(round((rng[0] - lo) / F.BW)), int(round((rng[1] - lo) / F.BW))
    return h[..., i0:i1]


def win_mask(sp, rng):
    e = np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
    ctr = 0.5 * (e[1:] + e[:-1])
    return np.abs(ctr - F.SP[sp]["m0"]) < F.SP[sp]["win"], e


def pm(x, e, d=1, pct=False):
    s = 100.0 if pct else 1.0
    return f"{x * s:+.{d}f} ± {e * s:.{d}f}" if e is not None else f"{x * s:+.{d}f}"

"""charmbkg.bkggen.api — JEDEN interfejs generatorów tła (etap AA, scalenie S/S2/S4/S6).

Kontrakt:
    gen = make_generator(name, feats, sp, **kw)      # name: s_copula | s_flow | copula | flow | salad | candrot
    gen.fit(SB)                                       # SB: DataFrame z kolumnami cech + "m" (+ "pt"); SALAD: fit(SB, mc=MC_SB)
    X   = gen.sample(m)                               # m: tablica mas w SR → DataFrame z kolumnami FEATS_H + "m"
    gen.feats_out                                     # lista kolumn wyjściowych (zawsze FEATS_H z flow.py S2)

Baza cech: FEATS_GEN + derive_H (S2) — więzy kinematyczne spełnione ściśle. Adaptery:
  * S2Copula / S2Flow   — CopulaBkg / FlowBkg z flow.py (S2): generują FEATS_GEN, H przez derive_H.
  * SCopula / SFlow     — BkgGenerator z gen.py (S): generuje cechy H wprost (łamie więzy; zostawione dla
                          porównywalności z liczbami S1/S7). Wyjście przycinane do FEATS_H.
  * Salad               — SaladReweighter (S4): nie generuje, zwraca PRAWDZIWE zdarzenia MC w SR z wagą `w`
                          (sample(m) ignoruje m poza wyborem regionu; wymaga fit(SB_data, mc=MC)).
  * Candrot             — candrot.make (S6): obiekty kontrolne z kandydatów (rot_pi|rot_rand|flip_eta|swap|mix),
                          sample(m) zwraca kontrole w oknie |m−m0|<sr_half (kształt masy, nie generator cech).
Nic z istniejących modułów nie jest zmieniane — api.py tylko opakowuje.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import flow as F2

FEATS_H = list(F2.FEATS_H)
FEATS_GEN = list(F2.FEATS_GEN)

_REGISTRY = {}


def register(name):
    def deco(cls):
        _REGISTRY[name] = cls
        cls.name = name
        return cls
    return deco


def available():
    return sorted(_REGISTRY)


def make_generator(name, feats=None, sp=None, species="D0", **kw):
    if name not in _REGISTRY:
        raise KeyError(f"generator '{name}' nieznany; dostępne: {available()}")
    sp = sp or dict(F2.SPECIES[species])
    return _REGISTRY[name](feats=feats, sp=sp, species=species, **kw)


class Generator:
    """Klasa bazowa. Podklasy: _fit(SB, **kw) i _sample(m) → ndarray (len(m), len(feats_out))."""

    feats_out = FEATS_H

    def __init__(self, feats=None, sp=None, species="D0", seed=0, **kw):
        self.feats = list(feats) if feats is not None else list(FEATS_H)
        self.sp, self.species, self.seed, self.kw = sp, species, seed, kw
        self.fitted = False

    # --- API
    def fit(self, SB: pd.DataFrame, **kw):
        if "m" not in SB:
            raise ValueError("SB musi mieć kolumnę 'm'")
        self._fit(SB, **kw)
        self.fitted = True
        return self

    def sample(self, m, **kw) -> pd.DataFrame:
        if not self.fitted:
            raise RuntimeError("najpierw fit(SB)")
        m = np.asarray(m, float)
        out = self._sample(m, **kw)
        if isinstance(out, pd.DataFrame):
            df = out
        else:
            df = pd.DataFrame(np.asarray(out), columns=self.feats_out)
            df["m"] = m
        if "pt" not in df and "log_pt" in df:
            df["pt"] = np.expm1(df["log_pt"])
        return df

    def sample_sr(self, n, rng=None, **kw):
        """Masy z modelu tła w SR (MassModel S2) + cechy."""
        mm = getattr(self, "mass_model", None)
        if mm is None:
            raise RuntimeError("generator nie ma modelu masy — podaj m do sample(m)")
        return self.sample(mm.sample(n, region="SR", rng=rng), **kw)

    # --- do nadpisania
    def _fit(self, SB, **kw):
        raise NotImplementedError

    def _sample(self, m, **kw):
        raise NotImplementedError


# ----------------------------------------------------------------------------- S2 (derive_H)
class _S2Base(Generator):
    gen_feats = FEATS_GEN

    def _fit(self, SB, epochs=None, **kw):
        missing = [f for f in self.gen_feats if f not in SB]
        if missing:
            raise ValueError(f"brak cech generowanych (FEATS_GEN) w SB: {missing}; użyj bkggen.prep")
        X = SB[self.gen_feats].to_numpy(float)
        m = SB["m"].to_numpy(float)
        self.mass_model = F2.MassModel(self.sp).fit(m)
        self._fit_model(X, m, epochs=epochs, **kw)

    def _sample(self, m, **kw):
        G = self.model.sample(m)
        return F2.derive_H(np.asarray(G), m, gen=self.gen_feats)


@register("copula")
class S2Copula(_S2Base):
    def _fit_model(self, X, m, epochs=None, deg=1, **kw):
        self.model = F2.CopulaBkg(self.gen_feats, self.sp, deg=deg, seed=self.seed).fit(X, m)


@register("flow")
class S2Flow(_S2Base):
    def _fit_model(self, X, m, epochs=None, **kw):
        fkw = {k: v for k, v in self.kw.items() if k in ("transforms", "hidden", "bins", "device")}
        self.model = F2.FlowBkg(self.gen_feats, self.sp, seed=self.seed, **fkw)
        self.model.fit(X, m, **({"epochs": epochs} if epochs else {}), **kw)


# ----------------------------------------------------------------------------- S (gen.py, cechy H wprost)
class _SBase(Generator):
    model_name = "copula"

    def __init__(self, *a, **kw):
        import warnings
        warnings.warn(f"generator '{getattr(self, 'name', self.model_name)}' (S, cechy H wprost) jest przestarzały — "
                      f"użyj '{self.model_name}' (S2, derive_H) albo 'salad'", DeprecationWarning, stacklevel=3)
        super().__init__(*a, **kw)

    def _fit(self, SB, **kw):
        from .gen import BkgGenerator
        feats = [f for f in self.feats if f in SB]
        lo, hi = self.sp["mass_range"]
        sbm = self.sp["sb_min"]
        m0 = self.sp["m0"]
        self._bg = BkgGenerator(feats, m0, lambda m: np.abs(np.asarray(m) - m0) > sbm, (lo, hi),
                                model=self.model_name, seed=self.seed, **self.kw)
        self._bg.fit(SB, **kw)
        self._feats_fit = feats
        self.feats_out = feats
        self.mass_model = F2.MassModel(self.sp).fit(SB["m"].to_numpy(float))

    def _sample(self, m, **kw):
        df = self._bg.sample(len(m), (self.sp["m0"] - self.sp["sr_half"], self.sp["m0"] + self.sp["sr_half"]), m=m)
        return df[self._feats_fit].to_numpy(float)


@register("s_copula")
class SCopula(_SBase):
    """PRZESTARZAŁE: pierwsza wersja kopuły (cechy H generowane wprost → łamie więzy kinematyczne; closure C2ST 0,96).
    Zostawione do odtworzenia liczb z docu S. Używaj `copula` (derive_H)."""
    model_name = "copula"


@register("s_flow")
class SFlow(_SBase):
    """PRZESTARZAŁE: pierwsza wersja flow (cechy H wprost). Używaj `flow` (derive_H) albo `salad`."""
    model_name = "flow"


# ----------------------------------------------------------------------------- S4 SALAD
@register("salad")
class Salad(Generator):
    """fit(SB_data, mc=MC_all) → sample(m): prawdziwe zdarzenia MC z |m−m0|<sr_half (lub najbliższe podanym m)
    z wagą `w` = p/(1−p). kind: sk_mlp (domyślnie, sklearn, bez torcha) | hgb | hgb_nom | mlp (torch, S4).
    `kind` można podać w make_generator("salad", kind=...) albo w fit(..., kind=...)."""

    def _fit(self, SB, mc=None, kind=None, epochs=12, **kw):
        from .salad8 import SaladReweighter8
        kind = kind or self.kw.get("kind", "sk_mlp")
        if mc is None:
            raise ValueError("Salad.fit wymaga mc=DataFrame (MC tła, cały zakres mas)")
        feats = [f for f in self.feats if f in SB and f in mc]
        self._feats_fit = feats
        self.feats_out = feats
        sb_mc = np.abs(mc["m"].to_numpy(float) - self.sp["m0"]) > self.sp["sb_min"]
        self.rw = SaladReweighter8(feats, self.sp, kind=kind, seed=self.seed).fit(
            mc.loc[sb_mc, feats].to_numpy(float), mc.loc[sb_mc, "m"].to_numpy(float),
            SB[feats].to_numpy(float), SB["m"].to_numpy(float), epochs=epochs)
        self.mc = mc.reset_index(drop=True)
        self.mass_model = F2.MassModel(self.sp).fit(SB["m"].to_numpy(float))

    def _sample(self, m, **kw):
        mm = self.mc["m"].to_numpy(float)
        sr = np.abs(mm - self.sp["m0"]) < self.sp["sr_half"]
        rng = np.random.default_rng(self.seed)
        idx = np.flatnonzero(sr)
        if len(idx) == 0:
            raise RuntimeError("MC nie pokrywa SR")
        pick = rng.choice(idx, size=len(m), replace=len(idx) < len(m))
        df = self.mc.loc[pick, self._feats_fit + ["m"]].reset_index(drop=True)
        df["w"] = self.rw.weights(df[self._feats_fit].to_numpy(float), df["m"].to_numpy(float)) \
            if hasattr(self.rw, "weights") else self.rw.predict_w(df[self._feats_fit].to_numpy(float), df["m"].to_numpy(float))
        return df


# ----------------------------------------------------------------------------- S6 candrot
@register("candrot")
class Candrot(Generator):
    """fit(SB) zapamiętuje kandydatów źródłowych (ramka charmbkg `L`, nie cechy!); sample(m) → kontrole w oknie."""

    def _fit(self, SB, method="rot_rand", n_rep=3, vertex_mode="keep", **kw):
        self.L, self.method, self.n_rep, self.vmode, self.ckw = SB, method, n_rep, vertex_mode, kw

    def _sample(self, m, **kw):
        from . import candrot as R
        out = R.make(self.L, self.sp, self.species, self.method, n_rep=self.n_rep, seed=self.seed,
                     vertex_mode=self.vmode, **self.ckw)
        df = out[0] if isinstance(out, tuple) else out
        if "m" in df:
            keep = np.abs(df["m"].to_numpy(float) - self.sp["m0"]) < self.sp["sr_half"]
            df = df[keep].reset_index(drop=True)
        return df


try:  # B1: SALAD v1.0 (charmbkg/bkggen/salad_v1.py) rejestruje generator "salad_v1" — addytywnie, bez wpływu na resztę
    from . import salad_v1 as _salad_v1  # noqa: F401,E402
except Exception as _e:  # pragma: no cover
    import logging as _lg; _lg.getLogger("charmbkg.bkggen.api").debug("salad_v1 niedostępny: %s", _e)

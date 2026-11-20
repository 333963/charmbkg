"""AFP-MIX v2: tło przypadkowego tagu AFP (WP3) metodą mieszania protonów (event mixing).

Model. Zdarzenie z tagiem = prawdziwa dyfrakcja ⊕ niezależny proton (pile-up bez wierzchołka, halo, …). Obiekt tła:
kandydaci D zdarzenia docelowego („biorcy”) + zawartość AFP innego zdarzenia („dawcy”) z tej samej warstwy.

v2 (2026-11-17) — zmiany względem v1:

A. Estymator nominalny („A”): dawcy z wetem ξ_cal (ξ_cal > ξ_hi + δ + margines po stronie), **warstwy (run, LB) × klasa
   aktywności** (przedziały nSelBL × „z D/bez D”, opcjonalnie klasa trygera). Region A∨C = suma mieszań z pul
   jednostronnych (pula dwustronna ma inną aktywność niż biorcy). Częstość protonu przypadkowego zależy od aktywności,
   a pula ma inną aktywność niż biorcy z małym ξ_cal — warstwowanie usuwa to
   z konstrukcji (zamiast ograniczać pulę do zdarzeń z D). Hierarchia warstw: (run, LB, klasa) → (run, przedział μ, klasa)
   → (run, klasa) → (run); na każdym poziomie pełna pula poziomu. Overlay „union” liczony na jawnych parach
   (biorca, otagowany dawca): biorca z własnym protonem po stronie nie dostaje tagu przypadkowego (v1: poprawka
   iteracyjna „własnego okna” — usunięta). Błąd stat.: bootstrap Poissona po otagowanych dawcach (wspólne wagi dla
   wszystkich wariantów → błąd RÓŻNICY wariantów). Ekstrapolacja pula → obszar użycia: profil R(ξ_cal) = obs/pred
   w puli, liniowy w ln ξ_cal, wartość κ w ⟨ln ξ_cal⟩ obszaru użycia; bez korekty, √((κ−1)² + σ_κ²) jako systematyka.
B. Drugie estymatory (reguła „połowa różnicy”, CMS-TOTEM 2002.12146): (i) absolutne P(tag) z próbki quasi-minbias
   o pełnej aktywności (weto ξ_cal po OBU stronach, nSelBL ≥ 20) w przedziałach μ, bez warstw aktywności — styl ALFA
   1911.00453; (ii) normalizacja w obszarze kontrolnym Δ = ξ_cal − ξ_AFP > δ + margines u biorców spoza puli.
C. Bramki: G0 wejścia, G1 jednorodność puli w ξ_cal (Poisson LR), G2 domknięcie na zbiorach rozłącznych (klasy
   trygera), G5 kształt (Δ, ξ) w obszarze kontrolnym (Poisson LR, p z pseudoeksperymentów z bootstrapem dawców),
   G7 zgodność A–B(ii), G8 zgodność A–B(i), G9 obszar „dwuprotonowy” (prawdziwy tag po drugiej stronie).
   I3 (μ), I6 (η*) — tylko informacyjne. G4 z v1 usunięta. Opisy „co testuje / czego NIE”: ``GATE_DOC``.

Konwencje charmpp (gałąź minimum-bias): tag = dokładnie jeden proton po stronie i w oknie 0,030 ≤ ξ ≤ 0,100
(``tag_def='nside'``, dane h5 v2.0) albo dokładnie jeden w oknie (``'inwindow'``, toy); region `tag` = A lub C;
cięcie ξ_cal ≤ ξ1_AFP + 0,02 na każdej otagowanej stronie (ξ1 = najmniejsze ξ po stronie).
Zakres: zdarzenia z trygerów nie-AFP (W2) — zob. docs/afp_mix_method.md.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd
from scipy import stats

SIDES = ("A", "C")

# ----------------------------------------------------------------------------------------------------------------------
# mapowanie nazw: wewnętrzna → zewnętrzna (h5 v2.0/v2.1). Jedyne miejsce do zmiany przy innej konwersji.
# ----------------------------------------------------------------------------------------------------------------------
COLMAP: dict[str, str] = {
    "run": "run_number",
    "event": "event_number",
    "lb": "lumiBlock",
    "mu": "actual_mu",
    "npu": "nPileupVertices",
    "nsel": "nSelBL",
    "np_A": "afp_nProtons_A", "np_C": "afp_nProtons_C",
    "xi1_A": "afp_xi1_A", "xi1_C": "afp_xi1_C",
    "xi2_A": "afp_xi2_A", "xi2_C": "afp_xi2_C",
    "tag_A": "afp_tag_A", "tag_C": "afp_tag_C",
    "xical_A": "xi_cal_A", "xical_C": "xi_cal_C",
    # opcjonalne flagi zdarzenia (v2): zdarzenie z ≥ 1 kandydatem D, klasa trygera (int), próbka B(i)
    "hasD": "has_D", "trig": "trig_class", "b1": "b1_sample",
    # kandydaci (opcjonalne; tylko do szablonu masy i asymetrii η)
    "m": "Dmeson_mass",
    "pt": "Dmeson_pt",
    "eta": "Dmeson_eta",
}
EVENT_REQUIRED = ("run", "event", "lb", "mu", "np_A", "np_C", "xi1_A", "xi1_C", "xi2_A", "xi2_C", "xical_A", "xical_C")
EVENT_KEY = ("run", "event")


def standardize(df: pd.DataFrame, colmap: dict[str, str] | None = None, required=EVENT_REQUIRED) -> pd.DataFrame:
    """Zwraca kopię z nazwami wewnętrznymi (tylko kolumny z mapy, które są w ``df``). Brak wymaganej → KeyError."""
    cm = dict(COLMAP, **(colmap or {}))
    out = {k: df[v].to_numpy() for k, v in cm.items() if v in df.columns}
    miss = [f"{k} ({cm[k]})" for k in required if k not in out]
    if miss:
        raise KeyError(f"afp_mix: brak kolumn {miss}")
    return pd.DataFrame(out, index=df.index)


def events_from_candidates(cands: pd.DataFrame) -> pd.DataFrame:
    """Tabela zdarzeń z tabeli kandydatów (kolumny zdarzeniowe są powielone per kandydat): pierwszy wiersz zdarzenia."""
    cols = [c for c in cands.columns if c not in ("m", "pt", "eta")]
    return cands[cols].drop_duplicates(list(EVENT_KEY)).reset_index(drop=True)


def ekey(df) -> np.ndarray:
    """Klucz zdarzenia int64 = run·2⁴⁰ + event (event < 2⁴⁰)."""
    return df["run"].to_numpy().astype(np.int64) * (1 << 40) + df["event"].to_numpy().astype(np.int64)


# ----------------------------------------------------------------------------------------------------------------------
# konfiguracja
# ----------------------------------------------------------------------------------------------------------------------
@dataclass
class MixConfig:
    xi_lo: float = 0.030            # okno tagu (charmpp minbias)
    xi_hi: float = 0.100
    delta: float = 0.02             # ξ_cal ≤ ξ_AFP + delta
    tag_def: str = "inwindow"       # 'nside' (dane h5 v2.0: dokładnie 1 proton po stronie, w oknie) | 'inwindow' (toy)
    donor_mode: str = "xical_veto"  # 'xical_veto' (nominał) | 'all' (naiwne; tylko test) | 'mask'
    veto_margin: float = 0.02       # dawca: ξ_cal > ξ_hi + delta + veto_margin
    pool_sides: str = "region"      # 'region' (weto po stronach regionu) | 'both' (pełna aktywność: obie strony)
    donor_mask_col: str | None = None
    strata: str = "lb"              # 'lb': (run, LB, klasa) → (run, μ, klasa) → (run, klasa) → (run) | 'mu': od poziomu μ
    mu_edges: tuple = field(default_factory=lambda: tuple(np.geomspace(0.002, 0.2, 21)))
    act_col: str = "nsel"           # zmienna aktywności (v2); warstwy = przedziały act_edges
    act_edges: tuple = (10, 15, 20, 25, 30, 40, 50, 80)   # () → bez warstw aktywności (≈ v1)
    cls_extra: tuple = ("hasD",)    # dodatkowe kolumny klas (int < 64): 'hasD' (zdarzenie z kandydatem D), opcjonalnie 'trig'
    n_min_pool: int = 50            # min. liczba dawców w komórce; mniej → następny poziom hierarchii
    use_tag_column: bool = True     # dla okna nominalnego użyj afp_tag_* z ntupli
    cr_margin: float = 0.025        # obszar kontrolny: Δ = ξ_cal − ξ1 > delta + cr_margin (CMS-TOTEM: 0,048)
    cr_xmax: float | None = None    # obszar kontrolny: ξ_cal biorcy ≤ cr_xmax; domyślnie ξ_hi + δ (górny kres obszaru użycia)
    max_pairs: int = 400_000_000

    @property
    def veto_thr(self) -> float:
        return self.xi_hi + self.delta + self.veto_margin

    @property
    def cr_x(self) -> float:
        return self.xi_hi + self.delta if self.cr_xmax is None else self.cr_xmax


# ----------------------------------------------------------------------------------------------------------------------
# tag i cięcie
# ----------------------------------------------------------------------------------------------------------------------
def n_in_window(np_, xi1, xi2, lo, hi):
    np_ = np.asarray(np_)
    return ((np_ >= 1) & (xi1 >= lo) & (xi1 <= hi)).astype(np.int8) + ((np_ >= 2) & (xi2 >= lo) & (xi2 <= hi)).astype(np.int8)


def tag_window(np_, xi1, xi2, lo, hi, tag_def="inwindow"):
    """Tag z podsumowania (n, ξ1 ≤ ξ2). 'inwindow': dokładnie jeden proton w [lo, hi] (n ≥ 3 i ξ2 ≤ hi → trzeciego nie
    znamy → brak tagu). 'nside': dokładnie jeden proton po stronie i w oknie. Zwraca (tag, ≥ 1 znany w oknie)."""
    np_ = np.asarray(np_)
    n_in = n_in_window(np_, xi1, xi2, lo, hi)
    amb = (np_ >= 3) & (xi2 <= hi)
    anyin = (n_in >= 1) | (amb & (xi2 >= lo))
    if tag_def == "nside":
        return (np_ == 1) & (n_in == 1), anyin
    return (n_in == 1) & ~amb, anyin


def _nominal(cfg, window):
    return window is None or (abs(window[0] - cfg.xi_lo) < 1e-12 and abs(window[1] - cfg.xi_hi) < 1e-12)


def side_tags(ev: pd.DataFrame, cfg: MixConfig, window=None) -> dict:
    """{s: (tag_s, any_in_window_s)} dla okna (lo, hi); dla okna nominalnego tag z kolumny afp_tag_* (o ile jest)."""
    lo, hi = window or (cfg.xi_lo, cfg.xi_hi)
    out = {}
    for s in SIDES:
        t, anyin = tag_window(ev[f"np_{s}"].to_numpy(), ev[f"xi1_{s}"].to_numpy(), ev[f"xi2_{s}"].to_numpy(), lo, hi, cfg.tag_def)
        if _nominal(cfg, window) and cfg.use_tag_column and f"tag_{s}" in ev:
            t = ev[f"tag_{s}"].to_numpy().astype(bool)
        out[s] = (t, anyin)
    return out


def region_pass(tA, tC, cutA, cutC, region: str):
    """Region charmpp: A | C | AorC (podwójne zostają), cięcie ξ_cal na każdej otagowanej stronie."""
    if region == "A":
        return tA & cutA
    if region == "C":
        return tC & cutC
    return (tA | tC) & (~tA | cutA) & (~tC | cutC)


def observed_region(ev: pd.DataFrame, cfg: MixConfig, region="AorC", window=None, xical_cut=True):
    """Obserwowany tag w regionie (z własnych protonów zdarzenia)."""
    tg = side_tags(ev, cfg, window)
    cut = {s: (ev[f"xical_{s}"].to_numpy() <= ev[f"xi1_{s}"].to_numpy() + cfg.delta) if xical_cut else
           np.ones(len(ev), bool) for s in SIDES}
    return region_pass(tg["A"][0], tg["C"][0], cut["A"], cut["C"], region)


def region_sides(region: str):
    return ("A", "C") if region == "AorC" else (region,)


# ----------------------------------------------------------------------------------------------------------------------
# pula dawców, klasy, warstwy
# ----------------------------------------------------------------------------------------------------------------------
def donor_pool(ev: pd.DataFrame, cfg: MixConfig, region="AorC", margin=None) -> np.ndarray:
    m = np.ones(len(ev), bool)
    if cfg.donor_mode == "xical_veto":
        thr = cfg.xi_hi + cfg.delta + (cfg.veto_margin if margin is None else margin)
        for s in (SIDES if cfg.pool_sides == "both" else region_sides(region)):
            m &= ev[f"xical_{s}"].to_numpy() > thr
    elif cfg.donor_mode == "mask":
        m &= ev[cfg.donor_mask_col].to_numpy().astype(bool)
    elif cfg.donor_mode != "all":
        raise ValueError(cfg.donor_mode)
    return m


def event_class(ev: pd.DataFrame, cfg: MixConfig) -> np.ndarray:
    """Klasa aktywności (int64): przedział ``act_col`` wg ``act_edges`` (+ kolumny ``cls_extra``)."""
    c = np.zeros(len(ev), np.int64)
    if len(cfg.act_edges):
        c = np.digitize(ev[cfg.act_col].to_numpy(), np.asarray(cfg.act_edges)).astype(np.int64)
    for col in cfg.cls_extra:
        if col not in ev:
            raise KeyError(f"afp_mix: klasa warstw wymaga kolumny '{col}' u biorców i dawców (MixConfig.cls_extra)")
        c = c * 64 + ev[col].to_numpy().astype(np.int64)
    return c


def _mu_bin(mu, edges):
    return np.clip(np.digitize(mu, edges), 0, len(edges))


def _level_keys(ev, cls, cfg: MixConfig, runs):
    r = np.searchsorted(runs, ev["run"].to_numpy()).astype(np.int64)
    mb = _mu_bin(ev["mu"].to_numpy(), np.asarray(cfg.mu_edges)).astype(np.int64)
    K = 100_000
    out = []
    if cfg.strata == "lb":
        out.append(("lb", (r * 1_000_000 + ev["lb"].to_numpy().astype(np.int64)) * K + cls))
    out += [("mu", (r * 64 + mb) * K + cls), ("run_cls", r * K + cls), ("run", r)]
    return out


def strata_keys(targets: pd.DataFrame, donors: pd.DataFrame, cfg: MixConfig):
    """(v1; używane przez f̄) Klucze warstw bez klas: (run, LB); LB z < n_min_pool dawców → (run, przedział μ) → (μ)."""
    def k_lb(d):
        return d["run"].to_numpy().astype(np.int64) * 100_000 + d["lb"].to_numpy().astype(np.int64)

    def k_mu(d, with_run=True):
        b = _mu_bin(d["mu"].to_numpy(), np.asarray(cfg.mu_edges)).astype(np.int64)
        r = d["run"].to_numpy().astype(np.int64) if with_run else np.zeros(len(d), np.int64)
        return -(r * 1000 + b + 1) - (0 if with_run else 10**12)

    kt, kd = (k_lb(targets), k_lb(donors)) if cfg.strata == "lb" else (k_mu(targets), k_mu(donors))
    for level in ("run_mu", "mu") if cfg.strata == "lb" else ("mu",):
        cnt = pd.Series(kd).value_counts()
        poor_d = cnt.reindex(kd).fillna(0).to_numpy() < cfg.n_min_pool
        poor_t = cnt.reindex(kt).fillna(0).to_numpy() < cfg.n_min_pool
        if not (poor_d.any() or poor_t.any()):
            break
        wr = level == "run_mu"
        kd = np.where(poor_d, k_mu(donors, wr), kd)
        kt = np.where(poor_t, k_mu(targets, wr), kt)
    return kt, kd


# ----------------------------------------------------------------------------------------------------------------------
# bootstrap dawców (wspólny dla wariantów: wagi przypisane do zdarzenia)
# ----------------------------------------------------------------------------------------------------------------------
class Boot:
    """Wagi Poissona(1) per zdarzenie-dawca (klucz ``ekey``), B replik. Te same wagi we wszystkich wariantach puli →
    błąd różnicy wariantów. Fluktuacje mianownika (liczby dawców w komórce, n ≫ k) pominięte: względnie ~ częstość tagu."""

    def __init__(self, keys, B=200, seed=0):
        self.keys = np.unique(np.asarray(keys, np.int64))
        self.B = B
        self.W = np.random.default_rng(seed).poisson(1.0, (len(self.keys), B)).astype(np.float32)

    def get(self, keys):
        keys = np.asarray(keys, np.int64)
        i = np.searchsorted(self.keys, keys)
        i[i == len(self.keys)] = 0
        if len(keys) and not np.array_equal(self.keys[i], keys):
            raise KeyError("Boot: dawca spoza tabeli wag")
        return self.W[i]

    @classmethod
    def for_events(cls, donors, cfg, B=200, seed=0):
        tg = side_tags(donors, cfg)
        return cls(ekey(donors)[tg["A"][0] | tg["C"][0]], B, seed)


# ----------------------------------------------------------------------------------------------------------------------
# rdzeń: pary (biorca, otagowany dawca) w komórkach warstw
# ----------------------------------------------------------------------------------------------------------------------
_TCOLS = ("np_A", "np_C", "xi1_A", "xi1_C", "xi2_A", "xi2_C", "xical_A", "xical_C")


class MixSum:
    """Suma mieszań po stronach (pule jednostronne) dla tych samych biorców i tej samej logiki regionu — region A∨C
    = mieszanie strony A + mieszanie strony C (tag przypadkowy po obu stronach naraz: O(p²), pominięty)."""

    def __init__(self, parts):
        self.parts = list(parts)
        m = self.parts[0]
        self.region, self.cfg, self.nt, self.own, self.own_tag, self.T = m.region, m.cfg, m.nt, m.own, m.own_tag, m.T
        self.p_tag = np.sum([x.p_tag for x in self.parts], 0)
        lv = np.array([x.level for x in self.parts])
        self.level = np.where(lv.min(0) < 0, -1, lv.max(0))
        self.n_pool = np.min([x.n_pool for x in self.parts], 0)
        self.ntd = sum(x.ntd for x in self.parts)

    def p_acc(self):
        return np.sum([x.p_acc() for x in self.parts], 0)

    def total(self, mask=None, tw=None):
        return float(sum(x.total(tw=tw) for x in self.parts))

    def boot(self, boot, mask=None, tw=None):
        return np.sum([x.boot(boot, tw=tw) for x in self.parts], 0)


class MixResult:
    """Wynik mieszania. Per biorca: ``p_tag`` (P(dawca komórki ma tag w regionie), z wykluczeniem siebie), ``n_pool``,
    ``level`` (poziom hierarchii; −1 = brak puli), ``q[s]`` (ułamek dawców komórki z protonem blokującym po stronie).
    Pary (tylko biorcy z ξ_cal ≤ próg puli po stronie tagu dawcy): ``pt`` (biorca), ``pd`` (indeks w tabeli otagowanych
    dawców ``TD``), ``pw`` (1/n_komórki, z poprawką 1/(1 − q) — biorcy, którzy już mają proton przypadkowy, wypadają
    z sumy). Overlay „union”: ``tag[s]``, ``dtag[s]`` (tag od dawcy), ``dl[s]`` (Δ_s = ξ_cal,s(biorca) − ξ1,s(union)), ``acc`` (region przechodzi
    tylko dzięki protonowi dawcy)."""

    def __init__(self, region, cfg, window, T, own_by, own_tag, p_tag, n_pool, level, q, pt, pd_, pw, TD, logic=None):
        self.pool_region, self.cfg, self.window = region, cfg, window
        self.T, self.own_by, self.own_tag = T, own_by, own_tag
        self.nt = len(p_tag)
        self.p_tag, self.n_pool, self.level, self.q = p_tag, n_pool, level, q
        self.pt, self.pd, self.TD = pt, pd_, TD
        self.ntd = len(TD["key"])
        lo, hi = window or (cfg.xi_lo, cfg.xi_hi)
        self.tag, self.dtag, self.dl = {}, {}, {}
        w = pw.astype(np.float64)
        ov = region_sides(region)                    # strony, po których nakładana jest zawartość AFP dawcy
        for s in SIDES:
            ntg = T[f"np_{s}"][pt]
            nd = TD[f"np_{s}"][pd_] if s in ov else np.zeros(len(pt), np.int32)
            nind = TD[f"nin_{s}"][pd_] if s in ov else 0
            x1 = np.minimum(np.where(ntg >= 1, T[f"xi1_{s}"][pt], np.inf), np.where(nd >= 1, TD[f"xi1_{s}"][pd_], np.inf))
            if cfg.tag_def == "nside":
                tg = (ntg + nd == 1) & (x1 >= lo) & (x1 <= hi)
            else:
                tg = (T[f"nin_{s}"][pt] + nind) == 1
            self.tag[s] = tg                                        # tag union (czyjkolwiek proton)
            self.dtag[s] = (tg & TD[f"tag_{s}"][pd_] & ~own_tag[s][pt]) if s in ov else np.zeros(len(pt), bool)   # tag od dawcy
            self.dl[s] = (T[f"xical_{s}"][pt] - x1).astype(np.float32)
            w = np.where(self.dtag[s], w / np.maximum(1 - q[s][pt], 0.5), w)
        self.pw = w
        self._set_logic(logic or region)

    def _set_logic(self, logic):
        """Logika regionu dla „acc”: para liczy się, gdy region ``logic`` przechodzi w overlayu, a bez dawcy nie."""
        cut = {s: self.dl[s] <= self.cfg.delta for s in SIDES}
        self.region = logic
        self.own = self.own_by[logic]
        self.acc = region_pass(self.tag["A"], self.tag["C"], cut["A"], cut["C"], logic) & ~self.own[self.pt]

    def relogic(self, logic):
        import copy
        m = copy.copy(self)
        m._set_logic(logic)
        return m

    # ---- sumy i bootstrap
    def p_acc(self, mask=None):
        m = self.acc if mask is None else mask
        return np.bincount(self.pt[m], self.pw[m], self.nt)

    def total(self, mask=None, tw=None):
        m = self.acc if mask is None else mask
        w = self.pw[m] if tw is None else self.pw[m] * tw[self.pt[m]]
        return float(w.sum())

    def donor_contrib(self, mask=None, tw=None):
        m = self.acc if mask is None else mask
        w = self.pw[m] if tw is None else self.pw[m] * tw[self.pt[m]]
        return np.bincount(self.pd[m], w, self.ntd)

    def boot(self, boot: Boot, mask=None, tw=None, contrib=None):
        """Repliki sumy Σ pw (po parach z maski, z wagami biorców tw) — wektor (B,)."""
        c = self.donor_contrib(mask, tw) if contrib is None else contrib
        nz = np.flatnonzero(c)
        if not len(nz):
            return np.zeros(boot.B)
        return c[nz] @ boot.get(self.TD["key"][nz]).astype(np.float64)

    def boot_p_tag(self, boot: Boot):
        """Repliki Σ_e p_tag(e) (wszyscy biorcy, wagi 1)."""
        return self.boot(boot, contrib=self.TD["tsum"])

    def xi_of_pairs(self, side):
        return self.TD[f"xi1_{side}"][self.pd]


def mix(targets: pd.DataFrame, donors: pd.DataFrame, cfg: MixConfig, region="AorC", window=None, margin=None,
        donor_sel=None, pairs=True, logic=None) -> MixResult:
    """Mieszanie w komórkach warstw (hierarchia poziomów, pełna pula na każdym poziomie). ``region`` wyznacza pulę
    (weto ξ_cal po stronach regionu), tag dawcy i strony, po których nakładana jest zawartość AFP dawcy; ``logic``
    (domyślnie = region) — region, którego przejście definiuje parę „acc” (np. pula strony A, logika A∨C).
    ``donor_sel`` — dodatkowa maska dawców. ``pairs=False`` → tylko p_tag (np. biorcy = zdarzenia puli)."""
    pool = donor_pool(donors, cfg, region, margin)
    if donor_sel is not None:
        pool = pool & np.asarray(donor_sel, bool)
    Dp = donors.iloc[np.flatnonzero(pool)]
    nt, ndn = len(targets), len(Dp)
    lo, hi = window or (cfg.xi_lo, cfg.xi_hi)
    tg = side_tags(Dp, cfg, window)
    tA, tC = tg["A"][0], tg["C"][0]
    tagged = (tA | tC) if region == "AorC" else tg[region][0]
    Dnp = {c: Dp[c].to_numpy() for c in _TCOLS}
    block = {}
    for s in SIDES:
        block[s] = (Dnp[f"np_{s}"] >= 1) if cfg.tag_def == "nside" else (n_in_window(Dnp[f"np_{s}"], Dnp[f"xi1_{s}"], Dnp[f"xi2_{s}"], lo, hi) >= 1)
    T = {c: targets[c].to_numpy() for c in _TCOLS}
    for s in SIDES:
        T[f"nin_{s}"] = n_in_window(T[f"np_{s}"], T[f"xi1_{s}"], T[f"xi2_{s}"], lo, hi)
    own_by = {r: observed_region(targets, cfg, r, window) for r in ("A", "C", "AorC")}
    ot = side_tags(targets, cfg, window)
    own_tag = {s: ot[s][0] for s in SIDES}
    runs = np.unique(np.concatenate([targets["run"].to_numpy(), Dp["run"].to_numpy()]))
    KT = _level_keys(targets, event_class(targets, cfg), cfg, runs)
    KD = _level_keys(Dp, event_class(Dp, cfg), cfg, runs)
    # wykluczenie samego siebie: wiersz biorcy w puli dawców
    dk = ekey(Dp)
    if ndn:
        o = np.argsort(dk, kind="stable")
        pos = np.clip(np.searchsorted(dk[o], ekey(targets)), 0, ndn - 1)
        self_d = np.where(dk[o][pos] == ekey(targets), o[pos], -1)
    else:
        self_d = np.full(nt, -1)
    tag_idx = np.flatnonzero(tagged)
    pos_in_tag = np.full(ndn, -1, np.int64)
    pos_in_tag[tag_idx] = np.arange(len(tag_idx))
    tsum = np.zeros(len(tag_idx))
    p_tag, n_pool, level = np.full(nt, np.nan), np.zeros(nt), np.full(nt, -1, np.int8)
    q = {s: np.zeros(nt) for s in SIDES}
    sides = region_sides(region)
    xmax = cfg.veto_thr if cfg.donor_mode == "xical_veto" else np.inf
    low = np.zeros(nt, bool)
    for s in sides:
        low |= T[f"xical_{s}"] <= xmax
    PT, PD, PW = [], [], []
    npairs = 0
    rem = np.arange(nt)
    for L, ((_, kt_all), (_, kd)) in enumerate(zip(KT, KD)):
        if not len(rem) or not ndn:
            break
        last = L == len(KT) - 1
        uk, inv, cnt = np.unique(kd, return_inverse=True, return_counts=True)
        kt = kt_all[rem]
        ps = np.clip(np.searchsorted(uk, kt), 0, len(uk) - 1)
        n = np.where(uk[ps] == kt, cnt[ps], 0)
        take = n >= (1 if last else cfg.n_min_pool)
        it, cell, nn = rem[take], ps[take], n[take].astype(np.float64)
        rem = rem[~take]
        if not len(it):
            continue
        tcell = inv[tag_idx]
        kcell = np.bincount(tcell, minlength=len(uk))
        sd = self_d[it]
        inp = sd >= 0
        n_eff = np.maximum(nn - inp, 1.0)
        k_eff = kcell[cell] - np.where(inp, tagged[np.maximum(sd, 0)], 0)
        p_tag[it], n_pool[it], level[it] = k_eff / n_eff, nn - inp, L
        for s in SIDES:
            q[s][it] = (np.bincount(inv, block[s], len(uk)) / cnt)[cell]
        tsum += np.bincount(cell, 1.0 / n_eff, len(uk))[tcell]
        st = inp & tagged[np.maximum(sd, 0)]                 # biorca sam jest otagowanym dawcą: bez wkładu do siebie
        if st.any():
            np.subtract.at(tsum, pos_in_tag[sd[st]], 1.0 / n_eff[st])
        if not pairs:
            continue
        o2 = np.argsort(tcell, kind="stable")
        tsorted = tag_idx[o2]
        start = np.cumsum(kcell) - kcell
        sel = low[it] & (kcell[cell] > 0)
        itl, cl, ne = it[sel], cell[sel], n_eff[sel]
        kk = kcell[cl]
        csum = np.cumsum(kk)
        CH = 20_000_000
        i0 = 0
        while i0 < len(itl):
            i1 = int(np.searchsorted(csum, (csum[i0 - 1] if i0 else 0) + CH, side="right"))
            i1 = max(i1, i0 + 1)
            k_ = kk[i0:i1]
            tot = int(k_.sum())
            rt = np.repeat(itl[i0:i1], k_)
            offs = np.arange(tot) - np.repeat(np.cumsum(k_) - k_, k_)
            rd = tsorted[np.repeat(start[cl[i0:i1]], k_) + offs]
            keep = rd != self_d[rt]
            ok = np.zeros(tot, bool)       # przydatna para: strona tagu dawcy z ξ_cal biorcy ≤ próg
            for s, tS in (("A", tA), ("C", tC)):
                if s in sides:
                    ok |= tS[rd] & (T[f"xical_{s}"][rt] <= xmax)
            keep &= ok
            PT.append(rt[keep].astype(np.int32))
            PD.append(pos_in_tag[rd[keep]].astype(np.int32))
            PW.append((1.0 / np.repeat(ne[i0:i1], k_))[keep].astype(np.float32))
            npairs += int(keep.sum())
            if npairs > cfg.max_pairs:
                raise MemoryError(f"afp_mix.mix: > {cfg.max_pairs} par (poziom {L}); zwiększ rozdrobnienie warstw")
            i0 = i1
    cat = (lambda a, dt: np.concatenate(a) if a else np.zeros(0, dt))
    TD = {c: Dnp[c][tag_idx] for c in _TCOLS}
    for s in SIDES:
        TD[f"nin_{s}"] = n_in_window(TD[f"np_{s}"], TD[f"xi1_{s}"], TD[f"xi2_{s}"], lo, hi)
    TD.update(tag_A=tA[tag_idx], tag_C=tC[tag_idx], key=dk[tag_idx], tsum=tsum)
    return MixResult(region, cfg, window, T, own_by, own_tag, p_tag, n_pool, level, q,
                     cat(PT, np.int32), cat(PD, np.int32), cat(PW, np.float32), TD, logic)


def mix_region(targets, donors, cfg: MixConfig, region="AorC", **kw):
    """Mieszanie dla regionu z pul jednostronnych: {strona: MixResult(pula strony, logika regionu)} i ich suma."""
    parts = {s: mix(targets, donors, cfg, s, logic=region, **kw) for s in region_sides(region)}
    return parts, (MixSum(parts.values()) if len(parts) > 1 else next(iter(parts.values())))


# ----------------------------------------------------------------------------------------------------------------------
# estymata A i pochodne
# ----------------------------------------------------------------------------------------------------------------------
def estimate(mx: MixResult, boot: Boot | None = None, tw=None, fbar=0.0) -> dict:
    """N_acc = Σ par „acc” (overlay union; bez poprawek iteracyjnych), N_tag = obserwowane w regionie, f_acc, ε_ξcal,
    N_diff = (N_tag − N_acc)/(1 − f̄). ``tw`` — wagi biorców (np. liczba kandydatów w oknie masy)."""
    w = np.ones(mx.nt) if tw is None else np.asarray(tw, float)
    n_tag = float(w[mx.own].sum())
    n_acc = mx.total(tw=None if tw is None else w)
    e_acc = float(mx.boot(boot, tw=None if tw is None else w).std(ddof=1)) if boot is not None else float("nan")
    own_nc = region_pass(mx.own_tag["A"], mx.own_tag["C"], True, True, mx.region)
    n_ptag = float(np.nansum((w * mx.p_tag)[~own_nc]))
    f_acc = n_acc / n_tag if n_tag > 0 else float("nan")
    e_tag2 = float((w[mx.own] ** 2).sum())
    e_f = f_acc * math.sqrt((e_acc / n_acc) ** 2 + e_tag2 / n_tag ** 2) if (n_tag > 0 and n_acc > 0 and np.isfinite(e_acc)) else float("nan")
    fb = fbar if np.isfinite(fbar) else 0.0
    return {"region": mx.region, "n_targets": int((w > 0).sum()), "n_tag": n_tag, "n_acc": n_acc, "n_acc_err_stat": e_acc,
            "f_acc": f_acc, "f_acc_err_stat": e_f, "eps_xical_acc": n_acc / n_ptag if n_ptag > 0 else float("nan"),
            "fbar": fbar, "n_diff": (n_tag - n_acc) / (1 - fb), "n_diff_err": math.sqrt(e_tag2 + (0 if not np.isfinite(e_acc) else e_acc ** 2)) / (1 - fb),
            "frac_no_pool": float(np.mean(mx.level < 0)), "n_tagged_donors": int(mx.ntd)}


def spoil_probability(targets: pd.DataFrame, donors: pd.DataFrame, cfg: MixConfig, region="AorC", margin=None, chunk=2048) -> np.ndarray:
    """f̄(e) = P(protony dawcy NAŁOŻONE na protony biorcy e zepsują jego region) — dla biorców z tagiem (≈ dyfrakcja).
    Overlay 'union' po stronie wg ``cfg.tag_def``; cięcie ξ_cal na ξ1. Warstwy (run, LB) bez klas aktywności (f̄ ~ 10⁻³)."""
    pool = donor_pool(donors, cfg, region, margin)
    D = donors.loc[pool].reset_index(drop=True)
    kt, kd = strata_keys(targets, D, replace(cfg, n_min_pool=max(cfg.n_min_pool, 200)))
    lo, hi, dl = cfg.xi_lo, cfg.xi_hi, cfg.delta
    anyp = (D["np_A"].to_numpy() >= 1) | (D["np_C"].to_numpy() >= 1)
    a = {}
    for s in SIDES:
        a[s] = (targets[f"np_{s}"].to_numpy().astype(np.int32), np.where(targets[f"np_{s}"].to_numpy() >= 1, targets[f"xi1_{s}"].to_numpy(), np.inf),
                targets[f"xical_{s}"].to_numpy(), D[f"np_{s}"].to_numpy().astype(np.int32),
                np.where(D[f"np_{s}"].to_numpy() >= 1, D[f"xi1_{s}"].to_numpy(), np.inf),
                n_in_window(targets[f"np_{s}"].to_numpy(), targets[f"xi1_{s}"].to_numpy(), targets[f"xi2_{s}"].to_numpy(), lo, hi).astype(np.int32),
                n_in_window(D[f"np_{s}"].to_numpy(), D[f"xi1_{s}"].to_numpy(), D[f"xi2_{s}"].to_numpy(), lo, hi).astype(np.int32))
    base = observed_region(targets, cfg, region)
    q = np.full(len(targets), np.nan)
    ot = pd.Series(np.arange(len(targets))).groupby(kt).indices
    od = pd.Series(np.arange(len(D))).groupby(kd).indices
    for key, it in ot.items():
        idd = od.get(key)
        if idd is None:
            continue
        sel = idd[anyp[idd]]
        if not len(sel):
            q[it] = 0.0
            continue
        for c0 in range(0, len(it), chunk):
            ic = it[c0:c0 + chunk]
            t, cut = {}, {}
            for s in SIDES:
                npe, x1e, xc, npd, x1d, nie, nid = a[s]
                x1 = np.minimum(x1e[ic, None], x1d[None, sel])
                if cfg.tag_def == "nside":
                    t[s] = ((npe[ic, None] + npd[None, sel]) == 1) & (x1 >= lo) & (x1 <= hi)
                else:
                    t[s] = (nie[ic, None] + nid[None, sel]) == 1
                cut[s] = xc[ic, None] <= x1 + dl
            ok = region_pass(t["A"], t["C"], cut["A"], cut["C"], region)
            q[ic] = (~ok).sum(1) / len(idd)
    q[~base] = np.nan
    return q


def acc_template(cands: pd.DataFrame, targets: pd.DataFrame, p_acc, edges, windows: dict | None = None):
    """Szablon masy przypadkowego tagu: każdy kandydat z wagą p_acc swojego zdarzenia (wektor per biorca)."""
    pw = pd.Series(np.asarray(p_acc), index=pd.Index(ekey(targets)))
    pw = pw[~pw.index.duplicated()]
    w = pw.reindex(ekey(cands)).fillna(0).to_numpy()
    h, _ = np.histogram(cands["m"].to_numpy(), edges, weights=w)
    out = {"hist": h, "edges": np.asarray(edges), "w": w}
    if windows:
        m = cands["m"].to_numpy()
        out["mean_p_acc"] = {k: float(w[(m >= lo) & (m <= hi)].mean()) for k, (lo, hi) in windows.items()}
    return out


def xi_distribution(mx: MixResult, side: str, edges, after_cut=True):
    """Rozkład ξ_AFP przypadkowych tagów po stronie (z par mieszania; po cięciu: pary „acc”)."""
    m = mx.dtag[side] & (mx.acc if after_cut else True)
    h, _ = np.histogram(mx.xi_of_pairs(side)[m], edges, weights=mx.pw[m])
    return h


# ----------------------------------------------------------------------------------------------------------------------
# narzędzia statystyczne
# ----------------------------------------------------------------------------------------------------------------------
def binned_fraction(x, y, w=None, edges=None):
    """Frakcja ⟨y⟩ w przedziałach x (y ∈ [0,1], wagi w), błąd: binomialny z n_eff."""
    w = np.ones(len(x)) if w is None else np.asarray(w, float)
    b = np.digitize(x, edges) - 1
    nb = len(edges) - 1
    ok = (b >= 0) & (b < nb)
    sw = np.bincount(b[ok], w[ok], nb)
    sy = np.bincount(b[ok], w[ok] * y[ok], nb)
    sw2 = np.bincount(b[ok], w[ok] ** 2, nb)
    sx = np.bincount(b[ok], w[ok] * x[ok], nb)
    with np.errstate(invalid="ignore", divide="ignore"):
        f = sy / sw
        neff = sw ** 2 / sw2
        err = np.sqrt(np.clip(f * (1 - f), 1e-12, None) / neff)
        xc = sx / sw
    return xc, f, err


def linfit(x, y, e):
    """Ważone LS y = a + b x → (a, b, cov, χ², ndf)."""
    x, y, e = np.asarray(x, float), np.asarray(y, float), np.asarray(e, float)
    ok = np.isfinite(y) & np.isfinite(e) & (e > 0) & np.isfinite(x)
    X = np.stack([np.ones(ok.sum()), x[ok]], 1)
    Wt = 1 / e[ok] ** 2
    cov = np.linalg.inv(X.T @ (X * Wt[:, None]))
    a, b = cov @ (X.T @ (Wt * y[ok]))
    chi2 = float(((y[ok] - a - b * x[ok]) ** 2 * Wt).sum())
    return float(a), float(b), cov, chi2, int(ok.sum() - 2)


def poisson_lr(obs, exp):
    """Poissonowski iloraz wiarygodności (Baker–Cousins): 2 Σ [e − o + o ln(o/e)]; puste przedziały danych bez artefaktu."""
    o, e = np.asarray(obs, float), np.maximum(np.asarray(exp, float), 1e-300)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(o > 0, o * np.log(o / e), 0.0)
    return float(2 * (e - o + t).sum())


def merge_groups(exp, min_exp=5.0):
    """Grupy sąsiednich przedziałów (od ogona) o oczekiwanej liczbie ≥ min_exp → wektor etykiet grup."""
    exp = np.asarray(exp, float)
    lab = np.zeros(len(exp), int)
    g, acc = 0, 0.0
    for i in range(len(exp) - 1, -1, -1):
        lab[i] = g
        acc += exp[i]
        if acc >= min_exp:
            g, acc = g + 1, 0.0
    if acc > 0 and g > 0:        # niedomknięta ostatnia grupa → dołącz do poprzedniej
        lab[lab == g] = g - 1
    return lab.max() - lab


def count_z(n_obs, n_pred, e_pred=0.0):
    """Istotność N_obs względem przewidywania: mała liczba oczekiwana (< 25) → dokładny ogon Poissona (z niepewnością
    przewidywania przez rozmycie Gamma ≈ ujemny dwumianowy), inaczej (N_obs − N_pred)/√(N_pred + σ²_pred). Znak: obs − pred."""
    if not (n_pred > 0):
        return 0.0
    if n_pred < 25:
        r = (e_pred / n_pred) ** 2 if e_pred > 0 else 0.0
        if r > 1e-6:      # NB: mean μ, var μ + r μ²
            size = 1 / r
            dist = stats.nbinom(size, size / (size + n_pred))
        else:
            dist = stats.poisson(n_pred)
        p_hi, p_lo = dist.sf(n_obs - 1), dist.cdf(n_obs)
        return float(stats.norm.isf(p_hi)) if p_hi < p_lo else -float(stats.norm.isf(p_lo))
    return float((n_obs - n_pred) / math.sqrt(n_pred + e_pred ** 2))


def eta_star(eta, side):
    """η* > 0 = półkula układu X (z dala od protonu); DIFF_model: tag A → η* = −η, tag C → η* = +η."""
    return np.where(np.asarray(side) == "A", -np.asarray(eta), np.asarray(eta))


def eta_asymmetry(eta_s, w=None):
    w = np.ones(len(eta_s)) if w is None else np.asarray(w, float)
    sp, sm = w[eta_s > 0].sum(), w[eta_s < 0].sum()
    n = sp + sm
    A = (sp - sm) / n
    neff = n ** 2 / (w ** 2).sum()
    return float(A), float(math.sqrt(max(1 - A * A, 1e-12) / neff))


# ----------------------------------------------------------------------------------------------------------------------
# bramki v2
# ----------------------------------------------------------------------------------------------------------------------
@dataclass
class GateConfig:
    z_max: float = 2.5          # |z| dla bramek zgodności
    p_min: float = 0.01         # p dla G1 i G5 (nominalny; v1 zmieniło próg G5 post hoc na 0,002 — cofnięte)
    flat_edges: tuple = (0.2, 0.3, 0.5, 1.0)   # wewnętrzne brzegi przedziałów ξ_cal puli (G1, ekstrapolacja)
    cr_delta_edges: tuple = (0.055, 0.065, 0.075)   # wewnętrzne brzegi Δ w obszarze kontrolnym (G5)
    cr_xi_edges: tuple = (0.04, 0.05, 0.06)      # wewnętrzne brzegi ξ_AFP w obszarze kontrolnym (G5)
    cr_margins: tuple = (0.01, 0.025, 0.04, 0.06)  # skan dolnej granicy Δ − δ (diagnostyka przecieku własnych protonów)
    margins: tuple = (0.0, 0.05, 0.10)   # warianty marginesu weta dawcy
    act_edges_coarse: tuple = (20, 40)   # wariant: zgrubne klasy aktywności; także klasy B(i)
    b1_mu_bins: int = 6                  # B(i): liczba kwantylowych przedziałów μ
    b1_nsel_min: int = 20                # B(i): „pełna aktywność” — nSelBL ≥ b1_nsel_min (oraz weto ξ_cal po obu stronach)
    n_boot: int = 200
    n_pseudo: int = 2000                 # pseudoeksperymenty G5
    mu_edges: tuple = (0.0, 0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.2)   # I3
    a_diff: float = 0.20                 # I6: asymetria η* prawdziwej dyfrakcji (DIFF_model)
    a_diff_err: float = 0.02
    min_exp_gate: float = 5.0            # bramka zliczeniowa z N_pred < min_exp_gate → N/A


GATE_DOC = {
    "G0_inputs": ("pokrycie warstw (≥ 99 % biorców w komórce zachowującej klasę aktywności: LB×klasa → μ×klasa → run×klasa; 0 bez puli), brak zdarzeń wyłącznie z trygera AFP",
                  "poprawności fizycznej; zakresu poza trygerami nie-AFP (W2) — tam metoda nie obowiązuje"),
    "G1_pool_homogeneity": ("czy po warstwowaniu częstość tagu w puli nie zależy od ξ_cal dawcy POWYŻEJ progu (Poisson LR obs vs pred w przedziałach ξ_cal)",
                            "obszaru poniżej progu (obszar użycia); zależności, które znoszą się w przedziale"),
    "G2_disjoint_closure": ("domknięcie na zbiorach ROZŁĄCZNYCH w obszarze anty-ξ_cal: tagi zdarzeń puli z najliczniejszej klasy trygera vs przewidywanie z dawców pozostałych klas (te same warstwy LB × aktywność × D)",
                            "obszaru użycia (ξ_cal małe); zależności od ξ_cal; przy jednej klasie trygera N/A (v1: G2 była w ~90 % tożsamością)"),
    "I2_hasD_vs_noD": ("(informacyjna) skład puli W1: tagi zdarzeń z D w obszarze anty vs przewidywanie z dawców BEZ D, per klasa aktywności — uzasadnia klasę „z D/bez D” w warstwach",
                       "niczego o estymacie nominalnej (dawcy bez D nie są w niej używani) — nie wchodzi do werdyktu"),
    "G5_cr_shape": ("kształt (Δ, ξ_AFP) tagów w obszarze kontrolnym Δ > δ + 0,025 u biorców SPOZA puli (ξ_cal ≤ 0,12) vs mieszanie; Poisson LR, p z pseudoeksperymentów z bootstrapem dawców",
                    "normalizacji (to G7); obszaru Δ ≤ δ; przy < ~50 zdarzeniach nie ma mocy"),
    "G7_methods_A_vs_Bii": ("GŁÓWNA: częstość przypadkowych najbliżej obszaru użycia — N_obs/N_pred w obszarze kontrolnym (ξ_cal 0,075–0,12, Δ > 0,045; „przesuwane weto”) = stosunek metod B(ii)/A",
                            "obszaru ξ_cal < 0,075 (tam prawdziwe tagi są nieodróżnialne); przeciek własnych protonów przez rozdzielczość ξ_cal zawyża R (skan marginesu)"),
    "G8_methods_A_vs_Bi": ("GŁÓWNA: zgodność z drugą metodą — absolutne P(tag) z próbki quasi-minbias o pełnej aktywności (weto ξ_cal po obu stronach, nSelBL ≥ 20) w przedziałach μ, bez warstw aktywności (styl ALFA: pełna niezależność)",
                           "obszaru użycia (obie pule mają ξ_cal > próg); przy tej samej próbce trygerowej pule częściowo się pokrywają"),
    "G9_opposite_side": ("obszar „dwuprotonowy”: zdarzenia z prawdziwym tagiem (po cięciu) po stronie przeciwnej i ξ_cal > próg po stronie badanej — częstość tagu po stronie badanej vs mieszanie (czy zdarzenia dyfrakcyjne mają tę samą częstość przypadkowych)",
                         "zależności od ξ_cal po stronie badanej; CD z protonem poza oknem nie szkodzi, ale CD z rozdzielczością ξ_cal zawyża obs"),
    "I3_mu_slope": ("(informacyjna) nachylenie f_tag(μ) dane vs mieszanie", "niczego przy σ_b(dane) ≫ b_mix — nie wchodzi do werdyktu"),
    "I6_eta_asymmetry": ("(informacyjna) f_acc z asymetrii η* kandydatów przy założonym A_diff", "A_diff z Pythii nie pasuje do kandydatów inkluzywnych — nie wchodzi do werdyktu"),
}
VERDICT_GATES = ("G0_inputs", "G1_pool_homogeneity", "G2_disjoint_closure", "G5_cr_shape", "G7_methods_A_vs_Bii",
                 "G8_methods_A_vs_Bi", "G9_opposite_side")


def _gate(name, status, **kw):
    if isinstance(status, (bool, np.bool_)):
        status = "PASS" if status else "FAIL"
    return {"gate": name, "status": status, "tests": GATE_DOC[name][0], "not_tests": GATE_DOC[name][1], **kw}


def _zgate(name, n_obs, n_pred, e_pred, gc: GateConfig, **kw):
    z = count_z(n_obs, n_pred, e_pred)
    r = n_obs / n_pred if n_pred > 0 else float("nan")
    er = math.sqrt(max(n_pred, 0) + e_pred ** 2) / n_pred if n_pred > 0 else float("nan")
    st = "N/A" if not (n_pred >= gc.min_exp_gate) else ("PASS" if abs(z) < gc.z_max else "FAIL")
    return _gate(name, st, n_obs=float(n_obs), n_pred=float(n_pred), n_pred_err=float(e_pred), ratio=r, ratio_err=er, z=z, **kw)


def gate_g0_inputs(mx: MixResult, n_afp_only=0):
    lv = mx.level
    fr = {int(k): float(np.mean(lv == k)) for k in np.unique(lv)}
    n_lv = 4 if mx.cfg.strata == "lb" else 3
    fine = float(np.mean((lv >= 0) & (lv <= n_lv - 2)))     # poziomy zachowujące klasę aktywności: LB×klasa, μ×klasa, run×klasa
    return _gate("G0_inputs", fine >= 0.99 and fr.get(-1, 0.0) == 0.0 and n_afp_only == 0, coverage=fine, level_fractions=fr,
                 frac_no_pool=fr.get(-1, 0.0), n_afp_trigger_only=int(n_afp_only), median_pool=float(np.median(mx.n_pool)))


def xical_profile(mxp: MixResult, x, tagged, gc: GateConfig, cfg: MixConfig, u_use=None):
    """Profil częstości w puli: biorcy = zdarzenia puli (p_tag z warstw, bez siebie). Przedziały x = ξ_cal powyżej progu:
    k_j (tagi), pred_j = Σ p_tag. Test jednorodności (Poisson LR, pred unormowane do Σk) i fit liniowy R(u), u = ln ξ_cal;
    κ = R(u_use)/R(ū_puli) — przewidywana zmiana częstości przy przejściu do obszaru użycia."""
    edges = np.concatenate([[cfg.veto_thr], [e for e in gc.flat_edges if e > cfg.veto_thr], [np.inf]])
    b = np.digitize(x, edges) - 1
    ok = (b >= 0) & np.isfinite(mxp.p_tag)
    nb = len(edges) - 1
    k = np.bincount(b[ok], tagged[ok], nb).astype(float)
    pred = np.bincount(b[ok], mxp.p_tag[ok], nb)
    n = np.bincount(b[ok], minlength=nb)
    u = np.bincount(b[ok], np.log(x[ok]), nb) / np.maximum(n, 1)
    exp = pred * (k.sum() / pred.sum()) if pred.sum() > 0 else pred
    lab = merge_groups(exp, 5.0)
    ng = lab.max() + 1
    lr = poisson_lr(np.bincount(lab, k, ng), np.bincount(lab, exp, ng))
    ndf = int(ng - 1)
    p = float(stats.chi2.sf(lr, ndf)) if ndf > 0 else 1.0
    out = {"edges": edges.tolist(), "k": k.tolist(), "pred": pred.tolist(), "n": n.tolist(), "u": u.tolist(), "lr": lr, "ndf": ndf, "p": p,
           "rate": float(k.sum() / max(n.sum(), 1))}
    R = k / np.maximum(exp, 1e-300)
    eR = np.sqrt(np.maximum(exp, 1e-300)) / np.maximum(exp, 1e-300)
    good = exp > 0
    if good.sum() >= 3:
        ubar = float(np.average(u[good], weights=exp[good]))
        a, s, cov, chi2, nd = linfit(u[good] - ubar, R[good], eR[good])
        out.update(slope=s, slope_err=math.sqrt(cov[1, 1]), u_pool=ubar, chi2_lin=chi2, ndf_lin=nd)
        if u_use is not None and np.isfinite(u_use):
            du = u_use - ubar
            out.update(u_use=float(u_use), kappa=1 + s * du / a, kappa_err=abs(du) * math.sqrt(cov[1, 1]) / a)
    return out


def control_region(mx: MixResult, side: str, cfg: MixConfig, dmin=None, tw=None):
    """Obszar kontrolny po stronie (styl CMS-TOTEM / „przesuwane weto”): tag po stronie, ξ_cal ≤ cr_x = ξ_hi + δ (biorca
    spoza puli, w zakresie ξ_cal obszaru użycia) i Δ = ξ_cal − ξ1 > δ + cr_margin (własny proton niemożliwy z dokładnością do rozdzielczości).
    → (maska obserwowanych biorców, maska par przewidywanych)."""
    dmin = cfg.delta + cfg.cr_margin if dmin is None else dmin
    xc, x1 = mx.T[f"xical_{side}"], mx.T[f"xi1_{side}"]
    obs = mx.own_tag[side] & (xc <= cfg.cr_x) & ((xc - x1) > dmin)
    pm = mx.dtag[side] & (xc[mx.pt] <= cfg.cr_x) & (mx.dl[side] > dmin)
    return obs, pm


def cr_counts(mxs: dict, cfg: MixConfig, boot: Boot | None, dmin=None):
    """Zliczenia obszaru kontrolnego dla stron w ``mxs`` ({strona: MixResult regionu tej strony})."""
    out = {"n_obs": 0.0, "n_pred": 0.0, "sides": {}}
    reps = np.zeros(boot.B) if boot is not None else None
    for s, mx in mxs.items():
        obs, pm = control_region(mx, s, cfg, dmin)
        no, npd = float(obs.sum()), mx.total(pm)
        rb = mx.boot(boot, pm) if boot is not None else None
        out["sides"][s] = {"n_obs": no, "n_pred": npd, "n_pred_err": float(rb.std(ddof=1)) if rb is not None else float("nan")}
        out["n_obs"] += no
        out["n_pred"] += npd
        if rb is not None:
            reps += rb
    out["n_pred_err"] = float(reps.std(ddof=1)) if reps is not None else 0.0
    out["reps"] = reps
    return out


def gate_g5_cr_shape(mxs: dict, cfg: MixConfig, gc: GateConfig, boot: Boot, rng=7):
    """G5: kształt (Δ, ξ_AFP) w obszarze kontrolnym, strony razem (wspólna normalizacja do N_obs). Statystyka: Poisson LR
    rozkładów brzegowych Δ i ξ per strona (po połączeniu przedziałów do oczekiwanej ≥ 5). p-wartość: pseudoeksperymenty —
    „prawda” = replika bootstrapu dawców (niepewność i korelacje mieszania), dane ~ Poisson, test względem nominału."""
    de = np.concatenate([[cfg.delta + cfg.cr_margin], gc.cr_delta_edges, [np.inf]])
    xe = np.concatenate([[cfg.xi_lo], gc.cr_xi_edges, [np.inf]])
    nd, nx = len(de) - 1, len(xe) - 1
    obs2, nom2, rep2 = [], [], []
    for s, mx in mxs.items():
        o, pm = control_region(mx, s, cfg)
        b_o = (np.digitize((mx.T[f"xical_{s}"] - mx.T[f"xi1_{s}"])[o], de) - 1) * nx + (np.digitize(mx.T[f"xi1_{s}"][o], xe) - 1)
        obs2.append(np.bincount(np.clip(b_o, 0, nd * nx - 1), minlength=nd * nx).astype(float))
        idx = np.flatnonzero(pm)
        xi_p = mx.T[f"xical_{s}"][mx.pt[idx]] - mx.dl[s][idx]      # ξ1 union = ξ dawcy
        b_p = np.clip((np.digitize(mx.dl[s][idx], de) - 1) * nx + (np.digitize(xi_p, xe) - 1), 0, nd * nx - 1)
        nom2.append(np.bincount(b_p, mx.pw[idx], nd * nx))
        C = np.zeros((mx.ntd, nd * nx))
        np.add.at(C, (mx.pd[idx], b_p), mx.pw[idx])
        nz = np.flatnonzero(C.sum(1))
        rep2.append(boot.get(mx.TD["key"][nz]).astype(np.float64).T @ C[nz] if len(nz) else np.zeros((boot.B, nd * nx)))
    obs2, nom2, rep2 = np.concatenate(obs2), np.concatenate(nom2), np.concatenate(rep2, 1)
    n_obs, n_nom = obs2.sum(), nom2.sum()
    if n_obs < 10 or n_nom <= 0:
        return _gate("G5_cr_shape", "N/A", n_obs=float(n_obs), n_pred=float(n_nom), note="za mało zdarzeń w obszarze kontrolnym")
    ns = len(mxs)
    # rozkłady brzegowe per strona: macierz rzutowania (bin 2D → bin brzegowy), potem łączenie do exp ≥ 5
    proj = []
    for i in range(ns):
        for kind in ("d", "x"):
            n1 = nd if kind == "d" else nx
            P = np.zeros((ns * nd * nx, n1))
            for j in range(nd * nx):
                P[i * nd * nx + j, (j // nx) if kind == "d" else (j % nx)] = 1
            exp1 = (nom2 * n_obs / n_nom) @ P
            lab = merge_groups(exp1, 5.0)
            Mg = np.zeros((n1, lab.max() + 1))
            Mg[np.arange(n1), lab] = 1
            proj.append(P @ Mg)
    P = np.concatenate(proj, 1)
    nbins = P.shape[1]

    def T(o, e):
        e = e * (o.sum() / max(e.sum(), 1e-300))
        return poisson_lr(o @ P, e @ P)

    t_obs = T(obs2, nom2)
    g = np.random.default_rng(rng)
    tp = np.empty(gc.n_pseudo)
    for i in range(gc.n_pseudo):
        truth = np.maximum(rep2[g.integers(0, rep2.shape[0])], 0)
        truth = truth * (n_obs / max(truth.sum(), 1e-300))
        tp[i] = T(g.poisson(truth).astype(float), nom2)
    p = float((1 + np.sum(tp >= t_obs)) / (1 + gc.n_pseudo))
    return _gate("G5_cr_shape", p >= gc.p_min, lr=t_obs, n_bins=int(nbins), p=p, lr_median_null=float(np.median(tp)), n_obs=float(n_obs),
                 n_pred=float(n_nom), obs_delta=[(obs2[i * nd * nx:(i + 1) * nd * nx].reshape(nd, nx).sum(1)).tolist() for i in range(ns)],
                 pred_delta=[((nom2 * n_obs / n_nom)[i * nd * nx:(i + 1) * nd * nx].reshape(nd, nx).sum(1)).tolist() for i in range(ns)],
                 obs_xi=[(obs2[i * nd * nx:(i + 1) * nd * nx].reshape(nd, nx).sum(0)).tolist() for i in range(ns)],
                 pred_xi=[((nom2 * n_obs / n_nom)[i * nd * nx:(i + 1) * nd * nx].reshape(nd, nx).sum(0)).tolist() for i in range(ns)])


def info_i3_mu(targets, mx: MixResult, gc: GateConfig):
    """I3 (informacyjna): f_tag(μ) = a + b μ po pełnej selekcji — nachylenie danych vs mieszania."""
    try:
        mu = targets["mu"].to_numpy()
        edges = np.asarray(gc.mu_edges)
        x1, f1, e1 = binned_fraction(mu, mx.own.astype(float), edges=edges)
        x2, f2, _ = binned_fraction(mu, mx.p_acc(), edges=edges)
        a1, b1, c1, *_ = linfit(x1, f1, e1)
        a2, b2, c2, *_ = linfit(x2, f2, np.full_like(f2, 1e-6))
        return _gate("I3_mu_slope", "INFO", slope_data=b1, slope_data_err=math.sqrt(c1[1, 1]), slope_mix=b2,
                     z_slope=(b1 - b2) / math.sqrt(c1[1, 1]), intercept_data=a1, intercept_mix=a2)
    except Exception as ex:   # np. μ stałe → fit niemożliwy
        return _gate("I3_mu_slope", "INFO", error=repr(ex))


def info_i6_eta(cands, targets, cfg: MixConfig, gc: GateConfig, f_mix):
    """I6 (informacyjna): asymetria η* kandydatów w zdarzeniach z jednym tagiem (po cięciu): f_acc(η) = 1 − A_obs/A_diff."""
    obs_A, obs_C = observed_region(targets, cfg, "A"), observed_region(targets, cfg, "C")
    single = obs_A ^ obs_C
    sd = pd.Series(np.where(obs_A, "A", "C")[single], index=pd.Index(ekey(targets)[single]))
    cs = sd.reindex(ekey(cands)).to_numpy()
    sel = pd.notna(cs)
    if sel.sum() < 20:
        return _gate("I6_eta_asymmetry", "INFO", n_cand=int(sel.sum()))
    A, eA = eta_asymmetry(eta_star(cands["eta"].to_numpy()[sel], cs[sel]))
    f_eta = 1 - A / gc.a_diff
    ef = math.hypot(eA / gc.a_diff, A * gc.a_diff_err / gc.a_diff ** 2)
    return _gate("I6_eta_asymmetry", "INFO", A_obs=A, A_err=eA, f_acc_eta=f_eta, f_acc_eta_err=ef, f_acc_mix=f_mix,
                 z=(f_eta - f_mix) / ef, n_cand=int(sel.sum()))


# ----------------------------------------------------------------------------------------------------------------------
# procedura v2
# ----------------------------------------------------------------------------------------------------------------------
def b1_config(cfg: MixConfig, gc: GateConfig, donors, b1_mask) -> MixConfig:
    """Konfiguracja metody B(i) (styl ALFA): pełna aktywność (weto ξ_cal po OBU stronach; maska próbki zawiera też
    nSelBL ≥ b1_nsel_min), absolutne P(tag) w kwantylowych przedziałach μ próbki — BEZ klas aktywności i bez klasy D."""
    mu = donors["mu"].to_numpy()[np.asarray(b1_mask, bool)]
    ed = np.unique(np.quantile(mu, np.linspace(0, 1, gc.b1_mu_bins + 1))[1:-1]) if len(mu) else np.array([])
    return replace(cfg, strata="mu", mu_edges=tuple(ed), act_edges=(), cls_extra=(), pool_sides="both", n_min_pool=max(cfg.n_min_pool, 200))


def _class_z(obs_c, pred_c, k_c):
    """χ² = Σ z_j² po klasach (z_j z ogona Poissona/NB; błąd przewidywania ~ 1/√k otagowanych dawców klasy)."""
    chi2, nd, rows = 0.0, 0, []
    for c in sorted(pred_c):
        o, pr, k = obs_c.get(c, 0.0), pred_c[c], k_c.get(c, 0.0)
        if pr <= 0:
            continue
        z = count_z(o, pr, pr / math.sqrt(max(k, 1.0)))
        z = float(np.clip(z, -8, 8))
        chi2, nd = chi2 + z * z, nd + 1
        rows.append({"cls": int(c), "n_obs": float(o), "n_pred": float(pr), "k_donors": float(k), "z": z})
    return chi2, nd, rows


def run_procedure(cands, events, cfg: MixConfig | None = None, gc: GateConfig | None = None, region="AorC", targets=None,
                  n_afp_only=0, rng=0, log=None, keep=False) -> dict:
    """Pełna procedura v2 (docs/afp_mix_method.md). ``events`` = dawcy (nazwy wewnętrzne; kolumny opcjonalne: ``hasD`` —
    domyślnie wyznaczana z ``cands``/``targets``; ``trig`` — klasa trygera (G2, wariant puli); ``b1`` — maska próbki B(i)).
    ``targets`` = biorcy (domyślnie zdarzenia kandydatów). Region A∨C liczony jako suma mieszań z pul jednostronnych.
    Zwraca estymatę A z rozbiciem błędu (stat, pula, ekstrapolacja, połowa różnicy metod), B(i), B(ii), bramki, werdykt.
    ``keep=True`` → obiekty MixResult w ``_obj``."""
    cfg = cfg or MixConfig()
    gc = gc or GateConfig()
    log = log or (lambda *_: None)
    if targets is None:
        targets = events_from_candidates(cands)
    events = events.reset_index(drop=True)
    targets = targets.reset_index(drop=True)
    if "hasD" not in events:
        events = events.assign(hasD=np.isin(ekey(events), np.unique(ekey(cands)) if cands is not None else ekey(targets)))
    if "hasD" not in targets:
        targets = targets.assign(hasD=True)
    hasD = events["hasD"].to_numpy().astype(bool)
    boot = Boot.for_events(events, cfg, gc.n_boot, rng)
    sides = region_sides(region)
    cfg_noD = replace(cfg, cls_extra=tuple(c for c in cfg.cls_extra if c != "hasD"))     # klasy bez podziału na z D / bez D

    # ---- A: nominał (pule jednostronne, logika regionu)
    ms, m0 = mix_region(targets, events, cfg, region)
    est = {}
    for reg, m in [(region, m0)] + ([(s, ms[s].relogic(s)) for s in sides] if region == "AorC" else []):
        fb = float("nan")
        if m.own.any():
            qv = spoil_probability(targets.loc[m.own].reset_index(drop=True), events, cfg, reg)
            fb = float(np.nanmean(qv)) if np.isfinite(qv).any() else float("nan")
        est[reg] = estimate(m, boot, fbar=fb)
        log(f"A[{reg}]: N_tag {est[reg]['n_tag']:.0f}, N_acc {est[reg]['n_acc']:.1f} ± {est[reg]['n_acc_err_stat']:.1f}")
    e0 = est[region]
    n0 = e0["n_acc"]
    rep0 = m0.boot(boot)

    def variant(c_v, **kw):
        _, mv = mix_region(targets, events, c_v, region, **kw)
        nv, rv = mv.total(), mv.boot(boot)
        d = nv / n0 - 1 if n0 else float("nan")
        with np.errstate(invalid="ignore", divide="ignore"):
            ed = float(np.nanstd(rv / np.where(rep0 > 0, rep0, np.nan) - 1, ddof=1)) if n0 else float("nan")
        return {"n_acc": nv, "ratio": nv / n0 if n0 else float("nan"), "delta": d, "delta_err": ed,
                "significant": max(abs(d) - ed, 0.0) if np.isfinite(d) and np.isfinite(ed) else float("nan"),   # istotna część różnicy
                "n_tagged_donors": int(mv.ntd), "frac_no_pool": float(np.mean(mv.level < 0))}

    # ---- warianty puli (pełna statystyka, wspólny bootstrap → błąd różnicy)
    vdef = [(f"margin_{mg:g}", cfg, dict(margin=mg)) for mg in gc.margins]
    vdef.append(("strata_mu" if cfg.strata == "lb" else "strata_lb", replace(cfg, strata="mu" if cfg.strata == "lb" else "lb"), {}))
    vdef.append(("act_coarse", replace(cfg, act_edges=tuple(gc.act_edges_coarse)), {}))
    if "trig" in events and "trig" in targets and "trig" not in cfg.cls_extra:
        vdef.append(("plus_trig_class", replace(cfg, cls_extra=cfg.cls_extra + ("trig",)), {}))
    variants = {}
    for name, c_v, kw in vdef:
        variants[name] = variant(c_v, **kw)
        log(f"wariant {name}: ×{variants[name]['ratio']:.3f} ± {variants[name]['delta_err']:.3f}")
    s_pool = max([v["significant"] for v in variants.values() if np.isfinite(v["significant"])] or [0.0])
    # informacyjne (NIE w systematyce): pula z dawcami bez D; bez warstw aktywności i bez klasy D (≈ v1)
    info_var = {"donors_D_and_noD": variant(cfg_noD), "v1_like_LB_only": variant(replace(cfg, act_edges=(), cls_extra=()))}
    for k, v in info_var.items():
        log(f"informacyjnie {k}: ×{v['ratio']:.3f} ± {v['delta_err']:.3f}")

    # ---- profil puli w ξ_cal per strona (klasy używane przez biorców): G1, ekstrapolacja, G9, G2
    used_cls = np.unique(event_class(targets, cfg))
    in_used = np.isin(event_class(events, cfg), used_cls)
    prof, g9, g2 = {}, {"n_obs": 0.0, "n_pred": 0.0, "sides": {}}, {"n_obs": 0.0, "n_pred": 0.0, "sides": {}}
    g2rep = np.zeros(boot.B)
    lr_sum, ndf_sum = 0.0, 0
    kap_num = kap_den = kap_var = 0.0
    tcl = None
    if "trig" in events and "trig" not in cfg.cls_extra:
        tv, tn = np.unique(events["trig"].to_numpy()[in_used], return_counts=True)
        tcl = tv[np.argmax(tn)] if len(tv) > 1 else None
    for s in sides:
        o = "C" if s == "A" else "A"
        ps = donor_pool(events, cfg, s) & in_used
        P = events.loc[ps].reset_index(drop=True)
        mp = mix(P, events, cfg, s, pairs=False)
        tag_s = mp.own_tag[s]
        msk = ms[s].acc & ms[s].dtag[s]
        wsum = ms[s].pw[msk].sum()
        u_use = float((ms[s].pw[msk] * np.log(np.maximum(ms[s].T[f"xical_{s}"][ms[s].pt[msk]], 1e-6))).sum() / wsum) if wsum > 0 else float("nan")
        prof[s] = xical_profile(mp, P[f"xical_{s}"].to_numpy(), tag_s, gc, cfg, u_use)
        lr_sum, ndf_sum = lr_sum + prof[s]["lr"], ndf_sum + prof[s]["ndf"]
        if "kappa" in prof[s]:
            kap_num += wsum * prof[s]["kappa"]
            kap_den += wsum
            kap_var += (wsum * prof[s]["kappa_err"]) ** 2
        # G9: prawdziwy tag (po cięciu) po stronie przeciwnej
        sel = mp.own_tag[o] & (P[f"xical_{o}"].to_numpy() <= P[f"xi1_{o}"].to_numpy() + cfg.delta) & np.isfinite(mp.p_tag)
        g9["sides"][s] = {"n_events": int(sel.sum()), "n_obs": float(tag_s[sel].sum()), "n_pred": float(mp.p_tag[sel].sum())}
        g9["n_obs"] += g9["sides"][s]["n_obs"]
        g9["n_pred"] += g9["sides"][s]["n_pred"]
        del mp
        # G2: zbiory rozłączne — zdarzenia puli z najliczniejszej klasy trygera vs dawcy z pozostałych klas
        if tcl is not None:
            is_t = P["trig"].to_numpy() == tcl
            if is_t.any() and (~is_t).any():
                m2 = mix(P.loc[is_t].reset_index(drop=True), events, cfg, s, donor_sel=events["trig"].to_numpy() != tcl, pairs=False)
                okp = np.isfinite(m2.p_tag)
                no, npd = float(m2.own_tag[s][okp].sum()), float(m2.p_tag[okp].sum())
                g2["sides"][s] = {"n_targets": int(is_t.sum()), "n_obs": no, "n_pred": npd, "frac_no_pool": float(np.mean(~okp))}
                g2["n_obs"] += no
                g2["n_pred"] += npd
                g2rep += m2.boot_p_tag(boot)
                del m2
        del P
    kappa = kap_num / kap_den if kap_den > 0 else float("nan")
    kappa_err = math.sqrt(kap_var) / kap_den if kap_den > 0 else float("nan")
    s_ext = math.sqrt((kappa - 1) ** 2 + kappa_err ** 2) if np.isfinite(kappa) else float("nan")
    p_g1 = float(stats.chi2.sf(lr_sum, ndf_sum)) if ndf_sum > 0 else 1.0
    n_pool_tags = sum(sum(prof[s]["k"]) for s in sides)
    g2err = float(g2rep.std(ddof=1))

    # ---- I2 (informacyjna): biorcy w obszarze anty (z D) vs dawcy bez D, per klasa aktywności — uzasadnienie klasy „hasD”
    i2 = None
    if "hasD" in cfg.cls_extra and (~hasD).any():
        oc, pc, kc = {}, {}, {}
        for s in sides:
            anti = donor_pool(targets, cfg, s)
            Ta = targets.loc[anti].reset_index(drop=True)
            m2 = mix(Ta, events, cfg_noD, s, donor_sel=~hasD, pairs=False)
            okp = np.isfinite(m2.p_tag)
            ca_ = event_class(Ta, replace(cfg_noD, cls_extra=()))
            for c in np.unique(ca_):
                mk = (ca_ == c) & okp
                oc[c] = oc.get(c, 0.0) + float(m2.own_tag[s][mk].sum())
                pc[c] = pc.get(c, 0.0) + float(m2.p_tag[mk].sum())
            Dn_ = events.loc[donor_pool(events, cfg, s) & ~hasD]
            cd_ = event_class(Dn_, replace(cfg_noD, cls_extra=()))
            tgd = side_tags(Dn_, cfg)[s][0]
            for c in np.unique(cd_):
                kc[c] = kc.get(c, 0.0) + float(tgd[cd_ == c].sum())
            del m2
        chi2, nd, rows = _class_z(oc, pc, kc)
        i2 = _gate("I2_hasD_vs_noD", "INFO", chi2=chi2, ndf=nd, p=float(stats.chi2.sf(chi2, nd)) if nd else 1.0,
                   n_obs=float(sum(oc.values())), n_pred=float(sum(pc.values())), per_class=rows, act_edges=list(cfg.act_edges))

    # ---- B(ii): normalizacja w obszarze kontrolnym (strony regionu)
    cr = cr_counts(ms, cfg, boot)
    R_cr = cr["n_obs"] / cr["n_pred"] if cr["n_pred"] > 0 else float("nan")
    cr_scan = []
    for mg in gc.cr_margins:
        c = cr_counts(ms, cfg, boot, dmin=cfg.delta + mg)
        cr_scan.append({"margin": mg, "n_obs": c["n_obs"], "n_pred": c["n_pred"], "ratio": c["n_obs"] / c["n_pred"] if c["n_pred"] > 0 else float("nan"),
                        "ratio_err": math.sqrt(c["n_pred"] + c["n_pred_err"] ** 2) / c["n_pred"] if c["n_pred"] > 0 else float("nan")})
    b2 = {"available": bool(cr["n_pred"] >= gc.min_exp_gate), "ratio": R_cr, "n_obs_cr": cr["n_obs"], "n_pred_cr": cr["n_pred"]}
    if cr["n_pred"] > 0:
        b2["n_acc"] = n0 * R_cr
        # błąd: Poisson N_obs(CR) ⊕ bootstrap wspólny (N_A i N_pred(CR) skorelowane przez dawców)
        with np.errstate(invalid="ignore", divide="ignore"):
            rb = rep0 * cr["n_obs"] / np.where(cr["reps"] > 0, cr["reps"], np.nan)
        e_b = float(np.nanstd(rb, ddof=1))
        b2["n_acc_err"] = math.hypot(e_b, b2["n_acc"] / math.sqrt(max(cr["n_obs"], 1.0)))
        b2["diff_rel"] = R_cr - 1
        b2["diff_rel_err"] = math.sqrt(cr["n_pred"] + cr["n_pred_err"] ** 2) / cr["n_pred"]

    # ---- B(i): druga pula (quasi-minbias, pełna aktywność, absolutne P(tag) w przedziałach μ)
    b1 = {"available": False}
    if "b1" in events and events["b1"].to_numpy().astype(bool).any():
        bm = events["b1"].to_numpy().astype(bool)
        if cfg.act_col in events:
            bm = bm & (events[cfg.act_col].to_numpy() >= gc.b1_nsel_min)
        both = donor_pool(events, replace(cfg, pool_sides="both"), region)
        cb = b1_config(cfg, gc, events, bm & both)
        _, mb = mix_region(targets, events, cb, region, donor_sel=bm)
        nb, rbb = mb.total(), mb.boot(boot)
        ed = float(np.std(rbb - rep0, ddof=1))
        tgb = side_tags(events.loc[bm & both], cfg)
        ntd = int((tgb["A"][0] | tgb["C"][0]).sum())
        b1 = {"available": bool(ntd >= 10 and nb > 0), "n_acc": nb, "n_acc_err": float(rbb.std(ddof=1)), "n_tagged_donors": ntd,
              "n_pool": int((bm & both).sum()), "diff_rel": nb / n0 - 1 if n0 else float("nan"),
              "diff_rel_err": ed / n0 if n0 else float("nan"), "frac_no_pool": float(np.mean(mb.level < 0)),
              "mu_edges": list(cb.mu_edges), "nsel_min": gc.b1_nsel_min}
        log(f"B(i): N_acc {nb:.1f} ± {b1['n_acc_err']:.1f} (otagowanych dawców {ntd})")
        del mb

    # ---- bramki
    gates = [gate_g0_inputs(m0, n_afp_only),
             _gate("G1_pool_homogeneity", p_g1 >= gc.p_min, lr=lr_sum, ndf=ndf_sum, p=p_g1, n_pool_tags=float(n_pool_tags),
                   per_side={s: {k: prof[s].get(k) for k in ("lr", "ndf", "p", "rate", "slope", "slope_err", "kappa", "kappa_err", "k", "pred", "edges")} for s in sides}),
             _zgate("G2_disjoint_closure", g2["n_obs"], g2["n_pred"], g2err, gc, per_side=g2["sides"], trig_class=int(tcl)) if g2["n_pred"] > 0 else
             _gate("G2_disjoint_closure", "N/A", note="brak kolumny klasy trygera albo jedna klasa"),
             gate_g5_cr_shape(ms, cfg, gc, boot),
             _zgate("G7_methods_A_vs_Bii", cr["n_obs"], cr["n_pred"], cr["n_pred_err"], gc, per_side=cr["sides"], margin_scan=cr_scan)]
    if b1["available"]:
        z8 = b1["diff_rel"] / b1["diff_rel_err"] if b1["diff_rel_err"] > 0 else 0.0
        gates.append(_gate("G8_methods_A_vs_Bi", abs(z8) < gc.z_max, n_acc_A=n0, n_acc_Bi=b1["n_acc"], ratio=b1["n_acc"] / n0,
                           ratio_err=b1["diff_rel_err"], z=z8, n_tagged_donors=b1["n_tagged_donors"]))
    else:
        gates.append(_gate("G8_methods_A_vs_Bi", "N/A", note="brak próbki B(i) (kolumna b1) lub < 10 otagowanych dawców"))
    gates.append(_zgate("G9_opposite_side", g9["n_obs"], g9["n_pred"], g9["n_pred"] / math.sqrt(max(n_pool_tags, 1.0)), gc, per_side=g9["sides"]))
    if i2 is not None:
        gates.append(i2)
    gates.append(info_i3_mu(targets, m0, gc))
    if cands is not None and "eta" in cands:
        gates.append(info_i6_eta(cands, targets, cfg, gc, e0["f_acc"]))
    failed = [g["gate"] for g in gates if g["gate"] in VERDICT_GATES and g["status"] == "FAIL"]
    na = [g["gate"] for g in gates if g["gate"] in VERDICT_GATES and g["status"] == "N/A"]
    verdict = "PASS" if not failed else "FAIL"

    # ---- rozbicie błędu
    s_stat = e0["n_acc_err_stat"] / n0 if n0 else float("nan")
    halves = {k: 0.5 * abs(b["diff_rel"]) for k, b in (("Bi", b1), ("Bii", b2)) if b.get("available") and np.isfinite(b.get("diff_rel", np.nan))}
    s_half = max(halves.values()) if halves else float("nan")
    parts = [s_stat, s_pool, s_ext, s_half]
    syst = {"stat": s_stat, "pool": s_pool, "extrapolation": s_ext, "half_diff_methods": s_half, "half_diff_by_method": halves,
            "total_rel": math.sqrt(sum(p * p for p in parts if np.isfinite(p))), "kappa": kappa, "kappa_err": kappa_err}
    out = {"estimate": e0, "estimates_by_region": est, "variants": variants, "info_variants": info_var, "profile": prof,
           "method_Bi": b1, "method_Bii": {k: v for k, v in b2.items()}, "syst_rel": syst, "gates": gates, "verdict": verdict,
           "verdict_failed": failed, "verdict_na": na}
    if keep:
        mx = {region: m0}
        if region == "AorC":
            mx.update({s: ms[s].relogic(s) for s in sides})
        out["_obj"] = {"mx": mx, "boot": boot, "hasD": hasD, "sides": ms}
    return out

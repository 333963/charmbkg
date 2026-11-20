#!/usr/bin/env python3
"""TRKJET2 — analiza lokalna (projekt B, 17.11.2026). Trzy testy z literatury (RECON_lit_B2 T2, T3, T8) deficytu izolacji D0←D*
(ptcone40 < 2 GeV) w klasie trygera dżetowego 428770/J12only (SF 0,46 ± 0,03 przy 0,98 ± 0,04 w 427929/MBTS):

 T1  STOŻKI PROSTOPADŁE / LOSOWE — czy nadwyżka aktywności w stożku D0 (dane J12 − MC) to piedestał zdarzenia (ta sama nadwyżka
     w stożkach R = 0,4 obróconych o ±90° w φ), czy aktywność skorelowana z kierunkiem D0;
 T2  PODPRÓBKA NIEOBCIĄŻONA TRYGEREM — twardy dżet torowy po stronie przeciwnej (|Δφ| > 2π/3, pT > X) i brak dżetu pT > Y w półkuli D0;
     dodatkowo bity ortogonalne (HLT_mb_*, trigTBP_L1_J12) w obu runach;
 T3  ZASTĘPNIK ENERGII NEUTRALNEJ (spekulatywne) — r_n = ξ_trk/ξ_cal jako zmienna przeważania.

Kryteria rozstrzygnięcia (stałe KRYTERIA poniżej) zapisane PRZED obejrzeniem danych; raport drukuje je dosłownie, werdykt liczony automatycznie,
wynik negatywny ma tę samą formę co pozytywny.

Wejście: runs/lxplus/tj2_<próbka>.npz (etap ORCH trkjet2; Q_lxplus/trkjet2). Narzędzia fitu/SBS/przeważania: B_trig/ISO_origin (te same co TRKJET).
Użycie: .venv_h5/bin/python charm_work/B_trig/TRKJET2/analyze.py [--in runs/lxplus] [--out out] [--report RAPORT.md] [--nproc 6] [--criteria-only]
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "ISO_origin"))
from iso_common import C, DR, chi2_const, classes, eff_w, fit_pf, sbs_w  # noqa: E402  (na poziomie modułu: procesy Pool muszą je widzieć)
import s2_universal_rw as S2  # noqa: E402

INF = 1e9
SQRTS = 13600.0
MWIN, STHR = 0.040, 0.6
CUT = 2.0
NFIT = 150
NWEAK, NBAD = 150, 50
NMC_OK, NMC_MIN = 100, 25
CLS = ["9:MB", "8:all", "8:J12only", "8:MBTSbased", "8:other", "8:notJ12"]
DISJ = ["9:MB", "8:J12only", "8:MBTSbased", "8:other"]
KEY2 = ["9:MB", "8:J12only"]
PT12 = [3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 15, 20, 30]
PT4 = ("pt", [3, 5, 8, 12, 30])
PT3 = ("pt", [3, 5, 8, 30])
S_PT = [("pt", PT12)]
S_REF = [PT4, ("e_lead_pt", [0, 5, 8, 12, INF]), ("e_nSelBL", [0, 35, 55, INF])]
REFNAME = "pT×lead_pt×nSelBL (ISO_origin)"
PERP_E = [0, 0.3, 0.6, 1.0, 1.5, 2.2, 3.2, 5.0, INF]
# T1: schematy przeważania MC w zmiennych stożków prostopadłych/losowego (żaden nie używa torów stożka D0)
T1_SCHEMES = [
    ("brak", []),
    ("pT", S_PT),
    (REFNAME, S_REF),
    ("pT×perp_avg", [PT4, ("perp_avg", PERP_E)]),
    ("pT×perp_max×perp_min", [PT4, ("perp_max", [0, 0.6, 1.2, 2.2, 4.0, INF]), ("perp_min", [0, 0.3, 0.8, 1.6, INF])]),
    ("pT×perp_avg×nSelBL", [PT4, ("perp_avg", [0, 0.6, 1.2, 2.2, 4.0, INF]), ("e_nSelBL", [0, 35, 55, INF])]),
    ("pT×perp_avg×lead_pt", [PT4, ("perp_avg", [0, 0.6, 1.2, 2.2, 4.0, INF]), ("e_lead_pt", [0, 5, 8, 12, INF])]),
    ("pT×stożek losowy", [PT4, ("rnd", PERP_E)]),
]
T1_PERP = [s for s, _ in T1_SCHEMES if "perp" in s or "losowy" in s]
# T2: skan
T2_X = [0, 5, 8, 10, 12, 15]
T2_Y = [INF, 5, 3]
T2_VETO = {"O": "o_tow_pt", "B [koło]": "b_tow_pt"}
T2_CLS = ["8:J12only", "8:all", "9:MB"]
T2_PRIM_Y, T2_PRIM_V = 5, "O"
AWAY_E = [0, 3, 5, 8, 12, 20, INF]
# T3
RN_E = [0, 0.05, 0.08, 0.12, 0.17, 0.25, INF]
U_E = [0, 0.2, 0.4, 0.6, 0.8, 1.0001]
SOFT_MB = ["HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_1_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_mbts_L1MBTS_1_1",
           "HLT_mb_mbts_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_2", "HLT_noalg_L1MBTS_2", "HLT_mb_sptrk_L1RD0_FILLED",
           "HLT_noalg_L1RD0_FILLED", "HLT_noalg_mb_L1RD0_FILLED", "HLT_mb_mbts_L1RD0_FILLED", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_sptrk_pt2_L1RD0_FILLED"]
HARD_MB = [f"HLT_mb_sptrk_pt{p}_L1{l}" for p in (4, 6, 8) for l in ("MBTS_2", "MBTS_1", "RD0_FILLED")]

KRYTERIA = """Zapisane 17.11.2026 przed uruchomieniem ekstrakcji `trkjet2` na pełnych danych (stałe w `analyze.py`, drukowane dosłownie). Wspólne:
cięcie izolacji ptcone40 < 2 GeV (gałąź ntupla), selekcja i klasy jak w TRKJET; ε z fitu Δm pass/fail przy N_sig(SBS) ≥ 150, inaczej SBS
(„słaby” N_sig < 150, „niemiarodajny” N_sig < 50 lub N_sig/σ < 3 — nie wchodzi do werdyktu). „Domknięcie” = |SF − 1| < 2σ JEDNOCZEŚNIE
w 8:J12only i 9:MB. Odniesienie: pT×lead_pt×nSelBL (TRKJET: J12only 0,761 ± 0,040). Schemat z udziałem cofnięć („cof.”) > 0,30 w J12only
nie liczy się jako domknięcie. Bieg na niepełnym wejściu albo N_sig < 150 w klasie kluczowej → „NIEROZSTRZYGNIĘTE”.

**T1 — stożki prostopadłe / losowe (piedestał czy aktywność skorelowana z D⁰).** Zmienne: `perp1` = Σ pT torów (pT > 0,2 GeV, ta sama definicja co
własny stożek D⁰) w JEDNYM stożku R = 0,4 o osi (η_D⁰, φ_D⁰ ± π/2), znak wg parzystości numeru zdarzenia; `perp_avg` = średnia z obu;
`rnd` = stożek o losowym φ (≥ 0,8 od D⁰ i od osi wiodącego dżetu torowego).
- K1.1 (czy tryger podnosi piedestał): SF_perp = ε(perp1 < 2 GeV; dane) / ε(MC przeważone w pT). TAK, gdy SF_perp(J12only) < 1 o > 3σ.
  Kontrola: |SF_perp(9:MB) − 1| < 2σ (inaczej MC źle opisuje zdarzenie podkładowe niezależnie od trygera — zaznaczyć).
- K1.2 (udział piedestału w nadwyżce): F_ped = [⟨perp_avg⟩_dane − ⟨perp_avg⟩_MC] / [⟨stożek D⁰⟩_dane − ⟨stożek D⁰⟩_MC] w J12only
  (średnie sygnału, MC przeważone w pT). F_ped ≥ 0,7 → piedestał; F_ped ≤ 0,3 → aktywność skorelowana z D⁰; pomiędzy → mieszane.
  Liczone tylko, gdy mianownik > 3σ.
- K1.3 (domknięcie): (a) SF po odjęciu piedestału — warunek (własny stożek − perp_avg) < 2 GeV w danych i MC (MC przeważone w pT);
  (b) SF nominalny po przeważeniu MC w pT×perp (schematy `T1_PERP`; zmienne nie używają torów stożka D⁰ → brak koła).
  TAK, gdy (a) lub najlepszy (najmniejsze χ² na 4 klasach rozłącznych) schemat (b) domyka. „Częściowo”, gdy reszta |SF − 1| w J12only
  spada co najmniej 2× względem odniesienia przy 9:MB nadal w 2σ.
- Werdykt T1: K1.3 TAK → „PIEDESTAŁ tłumaczy deficyt”; K1.3 częściowo → „CZĘŚCIOWO piedestał”; K1.3 NIE → „NIE piedestał” z dopiskiem
  wg K1.1/K1.2 (piedestał podniesiony, lecz nie tłumaczy / piedestał niepodniesiony → nadwyżka skorelowana z D⁰).

**T2 — podpróbka nieobciążona trygerem.** Topologia T(X, Y): max pT dżetu torowego (kol. B) w |Δφ(dżet, D⁰)| > 2π/3 większy od X
oraz brak dżetu pT > Y w półkuli D⁰ (|Δφ| < π/2). Weto w wersji **O** (dżety z torów poza stożkiem D⁰ — bez koła; wersja główna)
i **B** (z torami stożka — [koło], informacyjnie). MC z tym samym wymaganiem; n_MC raportowane; n_MC < 100 = „MC słabe”, n_MC < 25 → brak SF.
- Punkt główny: Y = 5 GeV (O), największe X ∈ {5, 8, 10, 12, 15} z N_sig(J12only) ≥ 150 i n_MC ≥ 100. Brak takiego X → „NIEROZSTRZYGNIĘTE”.
- K2.1: SF_T = ε_fit(dane, T) / ε_MC(T, bez wag) w punkcie głównym; porównanie z SF inkluzywnym J12only (ta sama metoda).
  TAK („otoczenie D⁰ nieobciążone”), gdy |SF_T − 1| < 2σ oraz SF_T − SF_inkl > 2σ. CZĘŚCIOWO, gdy SF_T − SF_inkl > 2σ, ale SF_T < 1 o > 2σ.
  NIE, gdy SF_T − SF_inkl ≤ 2σ (deficyt zostaje także wtedy, gdy tryger tłumaczy dżet po przeciwnej stronie).
- K2.2 (informacyjnie): SF_T w funkcji X przy Y = 5 — czy rośnie z X; ε(dane) w binach pT dżetu strony przeciwnej.
- K2.3 (bity ortogonalne; każda komórka liczy się przy N_sig ≥ 150): (i) 428770, zdarzenia z HLT minbias bez wymagania twardości
  (lista `SOFT_MB`) i trigTBP_L1_J12 = 0 — SF zgodny z 1 (|pull| < 2) → warunki runu nie są przyczyną; (ii) zdarzenia zapisane trygerem
  minbias (428770 SOFT_MB albo 9:MB), lecz z trigTBP_L1_J12 = 1 — SF < 1 o > 2σ → deficyt idzie za samą obecnością dżetu L1, nie za tym,
  który tryger zapisał zdarzenie. Uwaga zapisana z góry: trigTBP_L1_MBTS_* jest ustawiony w prawie każdym zdarzeniu z D* — sam w sobie
  NIE definiuje próbki nieobciążonej (zdarzenie zapisał J12); raportowany tylko ułamek.

**T3 — zastępnik energii neutralnej (SPEKULATYWNE).** r_n = ξ_trk/ξ_cal, ξ_trk = Σ pT·(e^{+η} + e^{−η})/√s po wszystkich torach drzewa,
ξ_cal = ξ_cal_A + ξ_cal_C (suma stron — niezależna od przypisania A/C); f_n = 1 − r_n.
- K3.0 (czy MC nadaje się do przeważania wprost): ⟨r_n⟩_MC / ⟨r_n⟩_9:MB ∈ [0,67; 1,5]. NIE → schemat „pT×r_n” jest nieważny, decyduje
  schemat kwantylowy „pT×u” (u = dystrybuanta r_n w klasie odniesienia tego samego runu: 9:MB dla 427929, 8:notJ12 dla 428770; MC: własna).
- K3.1: ⟨r_n⟩ sygnału w J12only różni się od 8:notJ12 (ten sam run, te same warunki pile-upu) o > 3σ.
- K3.2: ΔSF = SF(schemat decydujący) − SF(pT) w J12only. „Zmienia”, gdy |ΔSF| ≥ 0,05 i |ΔSF| > 2σ; „domyka” wg definicji wspólnej.
  Dodatkowo dane→dane w tym samym runie: ε(J12only) / ε przewidziane z 8:notJ12 w komórkach (pT, r_n) względem samych komórek pT.
- Werdykt T3: domyka → „f_n TŁUMACZY (spekulatywne — wymaga potwierdzenia)”; zmienia, nie domyka → „f_n MA WPŁYW, nie domyka”;
  inaczej → „f_n BEZ WPŁYWU”.

**Łącznie:** hipoteza „deficyt = obciążenie otoczenia D⁰ przez tryger dżetowy, usuwalne zmienną zdarzenia” jest potwierdzona, gdy domyka
T1 albo T3, i wsparta, gdy T2 = TAK. Jeżeli T1 = NIE, T2 = NIE i T3 bez wpływu: deficyt siedzi w samym stożku D⁰ przy dowolnie dobranym
otoczeniu → nie da się go zdjąć przeważeniem MC minbias; pozostaje ε z danych per strumień (RECON_lit_B2 T1) albo embedding (T7)."""




# ------------------------------------------------------------------------------------------------------------------ wczytanie
def load(indir, sample):
    z = np.load(os.path.join(indir, f"tj2_{sample}.npz"), allow_pickle=False)
    D = {k: z[k] for k in z.files if k != "meta"}
    meta = json.loads(str(z["meta"]))
    for k in list(D):
        if k.startswith("e_tp_"):
            D[k[2:]] = D[k]
    n = len(D["m"])
    D["dok"] = (D["nmatch"] + D["nd_lo"] == 2) & (D["same_entry"] == 1)
    D["_pas"] = D["ptcone40"] < CUT
    D["cone"] = D["v_s40_p2"].astype(np.float64)
    D["perp_avg"] = 0.5 * (D["p_s2_p"] + D["p_s2_m"]).astype(np.float64)
    D["perp_max"] = np.maximum(D["p_s2_p"], D["p_s2_m"]); D["perp_min"] = np.minimum(D["p_s2_p"], D["p_s2_m"])
    D["perp1"] = np.where(D["event"].astype(np.int64) % 2 == 0, D["p_s2_p"], D["p_s2_m"])
    D["rnd"] = np.where(D["r_ok"] == 1, D["r_s2"], D["perp1"])          # brak dopuszczalnego φ (rzadkie) → stożek prostopadły
    D["cone_sub"] = D["cone"] - D["perp_avg"]
    D["_pas_sub"] = D["cone_sub"] < CUT
    D["_pas_own"] = D["cone"] < CUT
    D["xi_trk"] = (D["xt_p2_pos"] + D["xt_p2_neg"]).astype(np.float64) / SQRTS
    if "xi_cal_A" in D and "xi_cal_C" in D:
        xc = D["xi_cal_A"].astype(np.float64) + D["xi_cal_C"].astype(np.float64)
    else:
        xc = np.full(n, np.nan)
    D["xi_cal"] = xc
    with np.errstate(divide="ignore", invalid="ignore"):
        D["rn"] = np.where(np.isfinite(xc) & (xc > 0), D["xi_trk"] / xc, np.nan)
    D["rn_ok"] = np.isfinite(D["rn"])
    D["rn"] = np.where(D["rn_ok"], D["rn"], 0.0)
    return D, meta


def base(D, lo=3, hi=30):
    return (np.abs(D["m"] - C.M0_D0) < MWIN) & (D["score"] > STHR) & (D["pt"] >= lo) & (D["pt"] < hi) & D["dok"]


def mc_sel(M, lo=3, hi=30):
    return base(M, lo, hi) & np.isin(M["truth"], [0, 1]) & (M["dm"] > 0.1435) & (M["dm"] < 0.1475)


def xbit(D, name, pref="x_tp_"):
    n = len(D["m"])
    for k in (pref + name, ("tp_" + name) if pref == "x_tp_" else None):
        if k and k in D:
            return D[k] == 1
    return np.zeros(n, bool)


def make_classes(D9, D8):
    c9 = classes(D9, 427929); c8 = classes(D8, 428770)
    S = {"9:MB": (D9, c9["MB427929"]), "8:all": (D8, c8["all"]), "8:J12only": (D8, c8["J12only"]), "8:MBTSbased": (D8, c8["MBTSbased"]),
         "8:other": (D8, c8["other"]), "8:notJ12": (D8, c8["all"] & c8["notJ12"])}
    X = {}                                              # klasy bitów ortogonalnych (T2 K2.3)
    for tag, D, cm in (("8", D8, c8["all"]), ("9", D9, c9["MB427929"])):
        soft = np.zeros(len(D["m"]), bool); hard = np.zeros(len(D["m"]), bool)
        for c in SOFT_MB:
            soft |= xbit(D, c)
        for c in HARD_MB:
            hard |= xbit(D, c)
        tj = xbit(D, "L1_J12", "x_tbp_"); have = "x_tbp_L1_J12" in D
        tm = xbit(D, "L1_MBTS_1", "x_tbp_")
        if tag == "8":
            X["8:J12only ∧ TBP_MBTS_1"] = (D, c8["J12only"] & tm)
            X["8:softMB (niezależnie od J12)"] = (D, cm & soft)
            X["8:softMB ∧ TBP_J12=0"] = (D, cm & soft & ~tj & have)
            X["8:softMB ∧ TBP_J12=1"] = (D, cm & soft & tj)
            X["8:sptrk_pt≥4 ∧ TBP_J12=0"] = (D, cm & hard & ~tj & have)
            X["8:sptrk_pt≥4 ∧ TBP_J12=1"] = (D, cm & hard & tj)
            X["8:all ∧ TBP_J12=0"] = (D, cm & ~tj & have)
            X["8:all ∧ TBP_J12=1"] = (D, cm & tj)
        else:
            X["9:MB ∧ TBP_J12=0"] = (D, cm & ~tj & have)
            X["9:MB ∧ TBP_J12=1"] = (D, cm & tj)
            X["9:MB ∧ sptrk_pt≥4 ∧ TBP_J12=0"] = (D, cm & hard & ~tj & have)
    return S, X


# ------------------------------------------------------------------------------------------------------------------ estymatory
def sbs_eff(D, m, flag):
    if m.sum() < 20:
        return {"ncand": int(m.sum()), "S": 0.0, "S_err": 0.0, "eps": np.nan, "err": np.nan}
    w = sbs_w(D["dm"][m], D["pt"][m], flag[m])
    e, er, S, eS = eff_w(w, flag[m].astype(float))
    return {"ncand": int(m.sum()), "S": S, "S_err": eS, "eps": e, "err": er}


def sig_w(D, m):
    """Wagi sygnału (odejmowanie pasm Δm), grupy tła wg pass/fail izolacji nominalnej."""
    return sbs_w(D["dm"][m], D["pt"][m], D["_pas"][m])


def wmean(w, x):
    T = w.sum()
    if T <= 0 or len(w) < 20:
        return np.nan, np.nan
    mu = float((w * x).sum() / T)
    return mu, float(np.sqrt((w ** 2 * (x - mu) ** 2).sum()) / T)


def mc_pt_w(w, ptd, ptm):
    """Wagi MC odtwarzające widmo pT sygnału danych (biny PT12; puste/ujemne biny danych → 0)."""
    ed = np.asarray(PT12, float)
    hd = np.maximum(np.histogram(ptd, ed, weights=w)[0], 0.0); hm = np.histogram(ptm, ed)[0].astype(float)
    r = np.where(hm > 0, hd / np.maximum(hm, 1), 0.0)
    return r[np.clip(np.digitize(ptm, ed) - 1, 0, len(ed) - 2)]


def mc_eff(flag, m):
    n = int(m.sum())
    if n == 0:
        return {"S": 0.0, "eps": np.nan, "err": np.nan}
    e = float(flag[m].mean())
    return {"S": float(n), "eps": e, "err": float(np.sqrt(max(e * (1 - e), 0.25 / n) / n))}


def quality(r):
    S, eS = r.get("S", 0.0), r.get("S_err", 0.0)
    if not np.isfinite(r.get("eps", np.nan)) or S < NBAD or (eS > 0 and S / eS < 3):
        return "niemiarodajny"
    return "słaby" if S < NWEAK else "ok"


def _fit(job):
    key, dp, df_, shape = job
    try:
        return key, fit_pf(dp, df_, shape, interval=False)
    except Exception as ex:  # noqa: BLE001
        return key, {"eps": np.nan, "err": np.nan, "N": 0.0, "N_err": 0.0, "det": False, "valid": False, "method": "fit", "reason": str(ex)}


def best(r):
    f = r.get("fit")
    if f and f.get("det") and f.get("valid") and np.isfinite(f.get("eps", np.nan)) and np.isfinite(f.get("err", np.nan)) and f["err"] > 0:
        return f["eps"], f["err"], "fit"
    return r.get("eps", np.nan), r.get("err", np.nan), "SBS"


def ratio(a, ea, b, eb):
    if a is None or b is None or not (np.isfinite(a) and np.isfinite(b)) or b <= 0:
        return np.nan, np.nan
    return a / b, float(np.hypot(ea / b, a * eb / b ** 2))


def pull1(sf, e):
    return float((sf - 1) / e) if sf is not None and np.isfinite(sf) and e and np.isfinite(e) and e > 0 else None


def sf_dd(Dt, mt, Dr, mr, scheme, pas="_pas", min_s=30.0):
    """Dane→dane: SF = ε(cel) / ε przewidziane z klasy odniesienia w komórkach `scheme` (wagi SBS po obu stronach; komórki odniesienia
    z Σw < min_s cofają się do mniejszej liczby wymiarów). Błąd: delta (cel) ⊕ stat. komórek odniesienia."""
    if mt.sum() < 50 or mr.sum() < 50:
        return None
    wt = sbs_w(Dt["dm"][mt], Dt["pt"][mt], Dt[pas][mt]); p = Dt[pas][mt].astype(float)
    wr = sbs_w(Dr["dm"][mr], Dr["pt"][mr], Dr[pas][mr]); pr = Dr[pas][mr].astype(float)
    if wt.sum() <= 0 or wr.sum() <= 0:
        return None
    Xt = {v: Dt[v][mt] for v, _ in scheme} | {"pt": Dt["pt"][mt]}; Xr = {v: Dr[v][mr] for v, _ in scheme} | {"pt": Dr["pt"][mr]}
    n = len(wt)
    eps = np.full(n, np.nan); var = np.zeros(n); cid = np.zeros(n, np.int64); fb = np.zeros(n, bool); todo = np.ones(n, bool)
    for upto in range(len(scheme), -1, -1):
        if upto == 0:
            e0 = (wr * pr).sum() / wr.sum(); v0 = (wr ** 2 * (pr - e0) ** 2).sum() / wr.sum() ** 2
            eps[todo] = e0; var[todo] = v0; cid[todo] = -1; fb[todo] = len(scheme) > 0
            break
        nc = int(np.prod([len(ed) - 1 for _, ed in scheme[:upto]]))
        ir = S2.cell_index(Xr, scheme, upto); it = S2.cell_index(Xt, scheme, upto)
        Sr = np.bincount(ir, weights=wr, minlength=nc); Pr = np.bincount(ir, weights=wr * pr, minlength=nc)
        e = Pr / np.where(Sr > 0, Sr, 1.0)
        v = np.bincount(ir, weights=wr ** 2 * (pr - e[ir]) ** 2, minlength=nc) / np.where(Sr > 0, Sr, 1.0) ** 2
        good = Sr[it] >= min_s
        take = todo & good
        eps[take] = e[it[take]]; var[take] = v[it[take]]; cid[take] = it[take] + upto * 10_000_000
        if upto < len(scheme):
            fb[take] = True
        todo &= ~good
        if not todo.any():
            break
    A, B = (wt * p).sum(), (wt * eps).sum()
    if B <= 0:
        return None
    sf = A / B
    v_t = (wt ** 2 * (p - sf * eps) ** 2).sum() / B ** 2
    u, inv = np.unique(cid, return_inverse=True)
    wc = np.bincount(inv, weights=wt); vc = np.zeros(len(u)); vc[inv] = var
    v_r = sf ** 2 * (wc ** 2 * vc).sum() / B ** 2
    return {"SF": float(sf), "err": float(np.sqrt(v_t + v_r)), "eps_data": float(A / wt.sum()), "eps_pred": float(B / wt.sum()),
            "S": float(wt.sum()), "S_ref": float(wr.sum()), "frac_fallback": float((wt * fb).sum() / wt.sum())}


def wcdf(xref, wref):
    """Dystrybuanta ważona (wagi SBS mogą być ujemne → monotonizacja). Zwraca funkcję x → u ∈ [0, 1]."""
    o = np.argsort(xref)
    xs = xref[o]; c = np.maximum.accumulate(np.cumsum(wref[o]))
    c = c / c[-1] if c[-1] > 0 else np.linspace(0, 1, len(c))
    return lambda x: np.clip(np.interp(x, xs, c), 0.0, 1.0)


# ------------------------------------------------------------------------------------------------------------------ formatowanie
def pm(v, e, nd=3):
    if v is None or not np.isfinite(v):
        return "–"
    return f"{v:.{nd}f} ± {e:.{nd}f}".replace(".", ",") if e is not None and np.isfinite(e) else f"{v:.{nd}f}".replace(".", ",")


def num(v, nd=2):
    return "–" if v is None or not np.isfinite(v) else f"{v:.{nd}f}".replace(".", ",")


def tak(b):
    return "–" if b is None else ("TAK" if b else "NIE")


def header(L, partial, indir, out, criteria_only=False):
    w = L.append
    w("# TRKJET2 — trzy testy z literatury deficytu izolacji D⁰ w trygerze dżetowym (projekt B)\n")
    if criteria_only:
        w(f"Stan: **kryteria zapisane {time.strftime('%Y-%m-%d %H:%M')}, wyników jeszcze nie ma** (etap ORCH `trkjet2` niewysłany / w toku). "
          "Po pobraniu `runs/lxplus/tj2_*.npz`: `.venv_h5/bin/python charm_work/B_trig/TRKJET2/analyze.py` nadpisze ten plik (sekcja „Kryteria” bez zmian).")
    else:
        w(f"Wygenerowane: `analyze.py`, wejście `{os.path.relpath(indir, HERE)}`, {time.strftime('%Y-%m-%d %H:%M')}. Liczby: `{os.path.relpath(out, HERE)}/trkjet2.json`.")
        if partial:
            w("\n**UWAGA: wejście NIEPEŁNE (bieg testowy na kilku plikach) — liczby służą tylko sprawdzeniu, że analiza działa; werdykty nie są wynikiem.**")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", default=os.path.join(HERE, "runs", "lxplus"))
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    ap.add_argument("--report", default=os.path.join(HERE, "RAPORT.md"))
    ap.add_argument("--nproc", type=int, default=6)
    ap.add_argument("--nfit", type=int, default=0, help="TEST: próg N_sig dla fitu pass/fail (kryteria wymagają 150 — nie zmieniać w produkcji)")
    ap.add_argument("--criteria-only", action="store_true", help="zapisz RAPORT.md tylko z sekcją Kryteria (przed danymi)")
    a = ap.parse_args(argv)
    if a.criteria_only:
        L = []
        header(L, False, a.indir, a.out, True)
        L.append("\n## Kryteria (zapisane przed obejrzeniem danych)\n"); L.append(KRYTERIA)
        open(a.report, "w").write("\n".join(L) + "\n")
        print("kryteria →", a.report)
        return 0
    if a.nfit:
        globals()["NFIT"] = a.nfit
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    D9, m9 = load(a.indir, "data_427929"); D8, m8 = load(a.indir, "data_428770")
    M, mM = load(a.indir, "fd_nd"); MB, mMB = load(a.indir, "mb_nd")
    metas = {"data_427929": m9, "data_428770": m8, "fd_nd": mM, "mb_nd": mMB}
    S, XC = make_classes(D9, D8)
    mm0 = mc_sel(M)
    shape = DR.shape_from_mc(M["dm"][mm0] * 1e3)
    partial = any(m.get("n_files", 0) < n for m, n in ((m8, 97), (m9, 458), (mM, 24), (mMB, 78)))
    R = {"partial": partial, "mc_n": int(mm0.sum()),
         "inputs": {k: {"n_pairs": v["n_pairs"], "n_files": v.get("n_files"), "frac_rnd_ok": v.get("frac_rnd_ok"), "frac_xi_nan": v.get("frac_xi_nan"),
                        "frac_daughters_ok": v.get("frac_daughters_ok")} for k, v in metas.items()},
         "cls": {}, "t1": {"means": {}, "perp": {}, "sub": {}, "rw": {}, "dd": {}, "diff": {}}, "t2": {"grid": {}, "away": {}, "bits": {}},
         "t3": {"means": {}, "rw": {}, "dd": {}, "diff": {}}}
    jobs = []

    def add_fit(sec, key, D, m, flag, r):
        if r["S"] >= NFIT:
            jobs.append(((sec, key), D["dm"][m & flag], D["dm"][m & ~flag], shape))

    def at(sec):
        o = R
        for k in sec:
            o = o[k]
        return o

    # ---- klasy: odtworzenie punktu wyjścia + wagi sygnału / wagi MC w pT (cache)
    W, WM = {}, {}
    ptM = M["pt"][mm0]
    for name in CLS:
        D, cm = S[name]
        m = base(D) & cm
        r = sbs_eff(D, m, D["_pas"]); r["mc"] = mc_eff(M["_pas"], mm0)
        R["cls"][name] = r
        add_fit(("cls",), name, D, m, D["_pas"], r)
        if m.sum() >= 50:
            W[name] = (m, sig_w(D, m)); WM[name] = mc_pt_w(W[name][1], D["pt"][m], ptM)

    def mc_mean(name, x):
        return wmean(WM[name], x[mm0].astype(float)) if name in WM else (np.nan, np.nan)

    # ================================================================================================================ T1
    for name in CLS:
        if name not in W:
            continue
        D, cm = S[name]; m, w = W[name]
        o = {"S": float(w.sum())}
        for v in ("cone", "perp_avg", "p_s2_p", "p_s2_m", "rnd", "cone_sub", "p_n2_p"):
            d = wmean(w, D[v][m].astype(float)); mc = mc_mean(name, M[v]); mc0 = wmean(np.ones(int(mm0.sum())), M[v][mm0].astype(float))
            o[v] = {"data": d[0], "data_err": d[1], "mc": mc[0], "mc_err": mc[1], "mc_unw": mc0[0]}
        ec = o["cone"]["data"] - o["cone"]["mc"]; eec = np.hypot(o["cone"]["data_err"], o["cone"]["mc_err"])
        ep = o["perp_avg"]["data"] - o["perp_avg"]["mc"]; eep = np.hypot(o["perp_avg"]["data_err"], o["perp_avg"]["mc_err"])
        o["excess_cone"], o["excess_cone_err"], o["excess_perp"], o["excess_perp_err"] = ec, eec, ep, eep
        if np.isfinite(ec) and eec > 0 and ec / eec > 3:
            o["F_ped"], o["F_ped_err"] = ratio(ep, eep, ec, eec)
        R["t1"]["means"][name] = o
        for cut in (1.0, 2.0):                                         # ε(perp1 < cut)
            fl = D["perp1"] < cut
            r = sbs_eff(D, m, fl)
            mc = mc_mean(name, (M["perp1"] < cut).astype(float))
            r["mc"] = {"eps": mc[0], "err": mc[1], "S": float(mm0.sum())}
            R["t1"]["perp"][f"{name}|{cut:g}"] = r
            if cut == CUT:
                add_fit(("t1", "perp"), f"{name}|{cut:g}", D, m, fl, r)
        for tag, pas in (("nominalna (gałąź)", "_pas"), ("własny stożek", "_pas_own"), ("własny − perp_avg", "_pas_sub")):
            r = sbs_eff(D, m, D[pas])
            mc = mc_mean(name, M[pas].astype(float)); r["mc"] = {"eps": mc[0], "err": mc[1], "S": float(mm0.sum())}
            r["mc_unw"] = mc_eff(M[pas], mm0)
            R["t1"]["sub"][f"{tag}|{name}"] = r
            add_fit(("t1", "sub"), f"{tag}|{name}", D, m, D[pas], r)
        for sname, scheme in T1_SCHEMES:
            try:
                r = S2.sf_rw(D, m, D["_pas"], M, mm0, M["_pas"], scheme)
            except Exception as ex:  # noqa: BLE001
                r = None; print("sf_rw", sname, name, ex)
            if r:
                R["t1"]["rw"][f"{sname}|{name}"] = r
        for sname, scheme in (("pT", S_PT), ("pT×lead_pt×nSelBL", S_REF)):       # odjęcie piedestału + przeważenie kinematyczne
            try:
                r = S2.sf_rw(D, m, D["_pas_sub"], M, mm0, M["_pas_sub"], scheme)
            except Exception:  # noqa: BLE001
                r = None
            if r:
                R["t1"]["rw"][f"[odjęcie] {sname}|{name}"] = r
    allsch = [s for s, _ in T1_SCHEMES] + ["[odjęcie] pT", "[odjęcie] pT×lead_pt×nSelBL"]
    for sname in allsch:
        v = [R["t1"]["rw"].get(f"{sname}|{n}") for n in DISJ]
        v = [x for x in v if x and x["S"] > NBAD and x["err"] > 0]
        R["t1"]["rw"][f"{sname}|chi2"] = {"chi2": float(sum(((x["SF"] - 1) / x["err"]) ** 2 for x in v)), "ndf": len(v)}
    # dane→dane: 9:MB jako odniesienie
    if "9:MB" in W:
        Dr, mr = S["9:MB"][0], W["9:MB"][0]
        for name in CLS[1:]:
            if name not in W:
                continue
            D, m = S[name][0], W[name][0]
            for sname, scheme in (("pT", [PT3]), ("pT×perp_avg", [PT3, ("perp_avg", [0, 0.6, 1.2, 2.2, INF])]),
                                  ("pT×r_n", [PT3, ("rn", [0, 0.08, 0.14, INF])])):
                r = sf_dd(D, m, Dr, mr, scheme)
                if r:
                    R["t1"]["dd"][f"{sname}|{name}"] = r
    # ε_iso w binach perp_avg
    R["t1"]["diff"] = {"edges": PERP_E, "cls": {}, "MC": [mc_eff(M["_pas"], mm0 & (M["perp_avg"] >= lo) & (M["perp_avg"] < hi)) for lo, hi in zip(PERP_E[:-1], PERP_E[1:])]}
    for name in DISJ:
        D, cm = S[name]; b = base(D) & cm
        R["t1"]["diff"]["cls"][name] = [sbs_eff(D, b & (D["perp_avg"] >= lo) & (D["perp_avg"] < hi), D["_pas"]) for lo, hi in zip(PERP_E[:-1], PERP_E[1:])]
    R["t1"]["diff"]["chi2"] = {f"{x}|{y}": S2.chi2_pair(R["t1"]["diff"]["cls"][x], R["t1"]["diff"]["cls"][y] if y != "MC" else R["t1"]["diff"]["MC"])
                               for x, y in (("9:MB", "8:J12only"), ("9:MB", "MC"), ("8:J12only", "MC"), ("8:J12only", "8:MBTSbased"))}

    # ================================================================================================================ T2
    for vname, vcol in T2_VETO.items():
        for X in T2_X:
            for Y in T2_Y:
                tM = mm0 & (M["b_away_pt"] > X) & (M[vcol] < Y)
                for name in T2_CLS:
                    D, cm = S[name]
                    m = base(D) & cm & (D["b_away_pt"] > X) & (D[vcol] < Y)
                    key = f"{vname}|{X}|{'-' if Y >= INF else Y}|{name}"
                    r = sbs_eff(D, m, D["_pas"]); r["mc"] = mc_eff(M["_pas"], tM)
                    if m.sum() >= 50 and tM.sum() >= NMC_MIN:
                        try:
                            r["rw_pt"] = S2.sf_rw(D, m, D["_pas"], M, tM, M["_pas"], [PT4])
                        except Exception:  # noqa: BLE001
                            pass
                    R["t2"]["grid"][key] = r
                    if vname == "O" or name == "8:J12only":
                        add_fit(("t2", "grid"), key, D, m, D["_pas"], r)
    for name in T2_CLS:                                                # ε w binach pT dżetu strony przeciwnej, weto O Y = 5
        D, cm = S[name]; b = base(D) & cm & (D["o_tow_pt"] < T2_PRIM_Y)
        R["t2"]["away"][name] = [sbs_eff(D, b & (D["b_away_pt"] >= lo) & (D["b_away_pt"] < hi), D["_pas"]) for lo, hi in zip(AWAY_E[:-1], AWAY_E[1:])]
    R["t2"]["away"]["MC"] = [mc_eff(M["_pas"], mm0 & (M["o_tow_pt"] < T2_PRIM_Y) & (M["b_away_pt"] >= lo) & (M["b_away_pt"] < hi)) for lo, hi in zip(AWAY_E[:-1], AWAY_E[1:])]
    for xn, (D, cm) in XC.items():
        m = base(D) & cm
        r = sbs_eff(D, m, D["_pas"]); r["mc"] = mc_eff(M["_pas"], mm0)
        if m.sum() >= 50:
            try:
                r["rw_pt"] = S2.sf_rw(D, m, D["_pas"], M, mm0, M["_pas"], S_PT)
            except Exception:  # noqa: BLE001
                pass
        R["t2"]["bits"][xn] = r
        add_fit(("t2", "bits"), xn, D, m, D["_pas"], r)
    j12 = S["8:J12only"]
    R["t2"]["tbp_mbts_frac"] = sbs_eff(j12[0], base(j12[0]) & j12[1], xbit(j12[0], "L1_MBTS_1", "x_tbp_"))
    R["t2"]["have_tbp"] = {"8": "x_tbp_L1_J12" in D8, "9": "x_tbp_L1_J12" in D9, "8_mbts": "x_tbp_L1_MBTS_1" in D8}

    # ================================================================================================================ T3
    okM = mm0 & M["rn_ok"]
    for name in CLS:
        if name not in W:
            continue
        D, cm = S[name]; m, w = W[name]
        ok = D["rn_ok"][m]
        o = {"frac_ok": float(w[ok].sum() / w.sum())}
        for v in ("rn", "xi_cal", "xi_trk"):
            d = wmean(w[ok], D[v][m][ok].astype(float)); o[v] = {"data": d[0], "data_err": d[1]}
        R["t3"]["means"][name] = o
    o = {"frac_ok": float(okM.sum() / max(mm0.sum(), 1))}
    for v in ("rn", "xi_cal", "xi_trk"):
        d = wmean(np.ones(int(okM.sum())), M[v][okM].astype(float)); o[v] = {"data": d[0], "data_err": d[1]}
    R["t3"]["means"]["MC fd_nd"] = o
    # kwantyl u względem klasy odniesienia tego samego runu
    for D in (D9, D8, M):
        D["u"] = np.full(len(D["m"]), 0.5)
    for ref, D in (("9:MB", D9), ("8:notJ12", D8)):
        if ref in W:
            m, w = W[ref]; ok = D["rn_ok"][m]
            if ok.sum() > 50 and w[ok].sum() > 0:
                D["u"] = wcdf(D["rn"][m][ok], w[ok])(D["rn"])
    if okM.sum() > 50:
        M["u"] = wcdf(M["rn"][okM], np.ones(int(okM.sum())))(M["rn"])
    T3_SCHEMES = [("pT", S_PT), ("pT4", [PT4]), ("pT×r_n (wprost)", [PT4, ("rn", RN_E)]), ("pT×u (kwantyl r_n)", [PT4, ("u", U_E)]),
                  ("pT×u×nSelBL", [PT4, ("u", [0, 0.33, 0.67, 1.0001]), ("e_nSelBL", [0, 35, 55, INF])]),
                  ("pT×u×perp_avg", [PT4, ("u", [0, 0.33, 0.67, 1.0001]), ("perp_avg", [0, 0.6, 1.2, 2.2, INF])])]
    for name in CLS:
        D, cm = S[name]
        m = base(D) & cm & D["rn_ok"]
        if m.sum() < 50 or okM.sum() < 50:
            continue
        for sname, scheme in T3_SCHEMES:
            try:
                r = S2.sf_rw(D, m, D["_pas"], M, okM, M["_pas"], scheme)
            except Exception as ex:  # noqa: BLE001
                r = None; print("sf_rw T3", sname, name, ex)
            if r:
                R["t3"]["rw"][f"{sname}|{name}"] = r
    if "8:notJ12" in W:                                               # dane→dane w tym samym runie
        Dr, cr = S["8:notJ12"]; mr = base(Dr) & cr & Dr["rn_ok"]
        for name in ("8:J12only",):
            D, cm = S[name]; m = base(D) & cm & D["rn_ok"]
            for sname, scheme in (("pT", [PT3]), ("pT×r_n", [PT3, ("u", [0, 0.33, 0.67, 1.0001])]), ("pT×perp_avg", [PT3, ("perp_avg", [0, 0.6, 1.2, 2.2, INF])])):
                r = sf_dd(D, m, Dr, mr, scheme)
                if r:
                    R["t3"]["dd"][f"{sname}|{name}"] = r
    R["t3"]["diff"] = {"edges": U_E, "cls": {}}
    for name in ("9:MB", "8:J12only", "8:notJ12"):
        D, cm = S[name]; b = base(D) & cm & D["rn_ok"]
        R["t3"]["diff"]["cls"][name] = [sbs_eff(D, b & (D["u"] >= lo) & (D["u"] < hi), D["_pas"]) for lo, hi in zip(U_E[:-1], U_E[1:])]
    R["t3"]["diff"]["MC"] = [mc_eff(M["_pas"], okM & (M["u"] >= lo) & (M["u"] < hi)) for lo, hi in zip(U_E[:-1], U_E[1:])]
    R["t3"]["diff"]["flat"] = {}
    for name, rows in R["t3"]["diff"]["cls"].items():
        g = [x for x in rows if np.isfinite(x["eps"]) and x["S"] > 25 and x["err"] > 0]
        R["t3"]["diff"]["flat"][name] = chi2_const([x["eps"] for x in g], [x["err"] for x in g]) if len(g) >= 2 else None

    print(f"fity: {len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    if jobs:
        with Pool(a.nproc) as P:
            for (sec, key), f in P.imap_unordered(_fit, jobs):
                at(sec)[key]["fit"] = f

    # ================================================================================================================ kryteria
    K = {}
    usable = (not partial) and all(R["cls"][n]["S"] >= NWEAK for n in KEY2)

    def sf_cls(r, mc_key="mc"):
        e, er, meth = best(r)
        sf, sfe = ratio(e, er, r[mc_key]["eps"], r[mc_key]["err"])
        return sf, sfe, meth

    # --- T1
    k = {}
    for n in KEY2:
        r = R["t1"]["perp"].get(f"{n}|{CUT:g}")
        if r:
            sf, sfe, meth = sf_cls(r)
            k[n] = {"SF_perp": sf, "err": sfe, "method": meth, "pull": pull1(sf, sfe), "eps_data": best(r)[0], "eps_mc": r["mc"]["eps"]}
    p8 = k.get("8:J12only", {}).get("pull"); p9 = k.get("9:MB", {}).get("pull")
    K["K1.1"] = {"cells": k, "pass": None if p8 is None else bool(p8 < -3), "sanity_MB": None if p9 is None else bool(abs(p9) < 2)}
    mj = R["t1"]["means"].get("8:J12only", {})
    F, Fe = mj.get("F_ped"), mj.get("F_ped_err")
    K["K1.2"] = {"F_ped": F, "err": Fe, "excess_cone": mj.get("excess_cone"), "excess_cone_err": mj.get("excess_cone_err"),
                 "excess_perp": mj.get("excess_perp"), "excess_perp_err": mj.get("excess_perp_err"),
                 "class": None if F is None else ("piedestał" if F >= 0.7 else "skorelowana z D⁰" if F <= 0.3 else "mieszane")}
    cand = []
    for sname in allsch:
        x = [R["t1"]["rw"].get(f"{sname}|{n}") for n in KEY2]
        if not all(x):
            continue
        pulls = [(v["SF"] - 1) / v["err"] for v in x]
        ch = R["t1"]["rw"][f"{sname}|chi2"]
        cand.append({"scheme": sname, "SF_MB9": x[0]["SF"], "err_MB9": x[0]["err"], "SF_J12": x[1]["SF"], "err_J12": x[1]["err"], "pulls": pulls,
                     "fallback_J12": x[1]["frac_fallback"], "closes": bool(all(abs(p) < 2 for p in pulls) and x[1]["frac_fallback"] <= 0.30),
                     "chi2": ch["chi2"], "ndf": ch["ndf"], "dist_J12": abs(x[1]["SF"] - 1), "mb_ok": bool(abs(pulls[0]) < 2)})
    refc = next((c for c in cand if c["scheme"] == REFNAME), None)
    perp_c = [c for c in cand if c["scheme"] in T1_PERP]
    sub_c = next((c for c in cand if c["scheme"] == "[odjęcie] pT"), None)
    bst = min(perp_c, key=lambda c: c["chi2"]) if perp_c else None
    closes = bool((bst and bst["closes"]) or (sub_c and sub_c["closes"]))
    part = bool(refc and any(c and c["mb_ok"] and c["dist_J12"] < 0.5 * refc["dist_J12"] for c in (bst, sub_c)))
    K["K1.3"] = {"best_rw": bst, "subtraction": sub_c, "ref": refc, "pass": closes, "partial": part and not closes}
    if not usable or bst is None:
        v1 = "NIEROZSTRZYGNIĘTE — niepełne wejście albo za mała statystyka"
    elif closes:
        v1 = "PIEDESTAŁ tłumaczy deficyt — domknięcie J12only i 9:MB jednocześnie"
    elif part:
        v1 = "CZĘŚCIOWO piedestał — reszta w J12only co najmniej 2× mniejsza niż w odniesieniu, bez domknięcia"
    else:
        tail = ("piedestał w J12 jest podniesiony, lecz nie tłumaczy deficytu" if K["K1.1"]["pass"] else "piedestał niepodniesiony")
        v1 = f"NIE piedestał — {tail}; nadwyżka w stożku D⁰: {K['K1.2']['class'] or 'udział piedestału nieokreślony'}"
    K["T1"] = v1
    # --- T2
    inc = R["cls"]["8:J12only"]; sfi, sfie, _ = sf_cls(inc)
    scan = []
    for X in T2_X:
        r = R["t2"]["grid"].get(f"{T2_PRIM_V}|{X}|{T2_PRIM_Y}|8:J12only")
        if not r:
            continue
        sf, sfe, meth = sf_cls(r)
        scan.append({"X": X, "N_sig": r["S"], "N_err": r["S_err"], "n_mc": int(r["mc"]["S"]), "eps": best(r)[0], "eps_err": best(r)[1], "eps_mc": r["mc"]["eps"],
                     "SF": sf, "err": sfe, "method": meth, "SF_rw_pt": (r.get("rw_pt") or {}).get("SF"), "err_rw_pt": (r.get("rw_pt") or {}).get("err"),
                     "ok": bool(r["S"] >= NFIT and r["mc"]["S"] >= NMC_OK and meth == "fit")})
    prim = max((s for s in scan if s["ok"] and s["X"] >= 5), key=lambda s: s["X"], default=None)
    k = {"scan": scan, "primary": prim, "SF_incl": sfi, "SF_incl_err": sfie}
    if prim:
        d = prim["SF"] - sfi; de = float(np.hypot(prim["err"], sfie))
        k.update({"diff": d, "diff_err": de, "returns": bool(d > 2 * de), "consistent_1": bool(abs(prim["SF"] - 1) < 2 * prim["err"])})
        k["pass"] = bool(k["returns"] and k["consistent_1"]); k["partial"] = bool(k["returns"] and not k["consistent_1"] and prim["SF"] < 1)
    K["K2.1"] = k
    g = [s for s in scan if s["ok"]]
    K["K2.2"] = {"rising": bool(len(g) >= 2 and g[-1]["SF"] > g[0]["SF"]), "n_points": len(g),
                 "slope_sigma": float((g[-1]["SF"] - g[0]["SF"]) / np.hypot(g[-1]["err"], g[0]["err"])) if len(g) >= 2 else None}
    kb = {}
    for xn, r in R["t2"]["bits"].items():
        sf, sfe, meth = sf_cls(r)
        kb[xn] = {"N_sig": r["S"], "SF": sf, "err": sfe, "method": meth, "pull": pull1(sf, sfe), "q": quality(r),
                  "SF_rw_pt": (r.get("rw_pt") or {}).get("SF"), "err_rw_pt": (r.get("rw_pt") or {}).get("err")}
    c_i = kb.get("8:softMB ∧ TBP_J12=0"); c_ii = [kb.get("8:softMB ∧ TBP_J12=1"), kb.get("9:MB ∧ TBP_J12=1")]
    c_ii = [c for c in c_ii if c and c["N_sig"] >= NFIT and c["pull"] is not None]
    K["K2.3"] = {"cells": kb, "i": None if not c_i or c_i["N_sig"] < NFIT or c_i["pull"] is None else bool(abs(c_i["pull"]) < 2),
                 "ii": None if not c_ii else bool(all(c["pull"] < -2 for c in c_ii)), "ii_n": len(c_ii),
                 "tbp_mbts_frac": R["t2"]["tbp_mbts_frac"]["eps"], "tbp_mbts_frac_err": R["t2"]["tbp_mbts_frac"]["err"]}
    if not usable or not prim:
        v2 = "NIEROZSTRZYGNIĘTE — brak punktu głównego z N_sig ≥ 150 i n_MC ≥ 100 (albo niepełne wejście)"
    elif k["pass"]:
        v2 = "TAK — w podpróbce „dżet po przeciwnej stronie” SF wraca do 1: otoczenie D⁰ nieobciążone"
    elif k["partial"]:
        v2 = "CZĘŚCIOWO — SF w podpróbce istotnie wyższy niż inkluzywny, ale nadal < 1"
    else:
        v2 = "NIE — deficyt zostaje także wtedy, gdy tryger tłumaczy dżet po przeciwnej stronie"
    K["T2"] = v2
    # --- T3
    mm_, m9_ = R["t3"]["means"].get("MC fd_nd", {}).get("rn", {}), R["t3"]["means"].get("9:MB", {}).get("rn", {})
    rr = mm_.get("data") / m9_["data"] if m9_.get("data") and np.isfinite(m9_.get("data", np.nan)) and mm_.get("data") is not None else None
    compat = None if rr is None or not np.isfinite(rr) else bool(0.67 <= rr <= 1.5)
    K["K3.0"] = {"ratio_mc_over_mb": rr, "pass": compat}
    a_, b_ = R["t3"]["means"].get("8:J12only", {}).get("rn", {}), R["t3"]["means"].get("8:notJ12", {}).get("rn", {})
    ns = None
    if a_.get("data") is not None and b_.get("data") is not None and np.isfinite(a_["data"]) and np.isfinite(b_["data"]):
        ns = float((a_["data"] - b_["data"]) / np.hypot(a_["data_err"], b_["data_err"]))
    K["K3.1"] = {"J12only": a_.get("data"), "J12only_err": a_.get("data_err"), "notJ12": b_.get("data"), "notJ12_err": b_.get("data_err"),
                 "nsigma": ns, "pass": None if ns is None else bool(abs(ns) > 3)}
    dec = "pT×r_n (wprost)" if compat else "pT×u (kwantyl r_n)"
    x0 = R["t3"]["rw"].get("pT4|8:J12only"); x1 = R["t3"]["rw"].get(f"{dec}|8:J12only"); y1 = R["t3"]["rw"].get(f"{dec}|9:MB")
    k = {"scheme": dec, "base": x0, "J12": x1, "MB9": y1}
    if x0 and x1:
        d = x1["SF"] - x0["SF"]
        k.update({"dSF": d, "changes": bool(abs(d) >= 0.05 and abs(d) > 2 * x1["err"]),
                  "closes": bool(y1 and abs(x1["SF"] - 1) < 2 * x1["err"] and abs(y1["SF"] - 1) < 2 * y1["err"] and x1["frac_fallback"] <= 0.30)})
    K["K3.2"] = k
    if not usable or not (x0 and x1):
        v3 = "NIEROZSTRZYGNIĘTE — niepełne wejście albo brak ξ_cal"
    elif k["closes"]:
        v3 = "f_n TŁUMACZY (spekulatywne — wymaga potwierdzenia) — domknięcie J12only i 9:MB"
    elif k["changes"]:
        v3 = "f_n MA WPŁYW, nie domyka"
    else:
        v3 = "f_n BEZ WPŁYWU — przeważenie w r_n nie zmienia SF w J12only"
    K["T3"] = v3
    if not usable:
        vj = "NIEROZSTRZYGNIĘTE (bieg testowy / za mała statystyka)"
    elif K["K1.3"]["pass"] or k.get("closes"):
        vj = "POTWIERDZONA: deficyt usuwalny zmienną zdarzenia (" + ", ".join(t for t, c in (("T1", K["K1.3"]["pass"]), ("T3", k.get("closes"))) if c) + \
             ")" + ("; wsparta przez T2" if K["K2.1"].get("pass") else "")
    elif K["K2.1"].get("pass"):
        vj = "WSPARTA przez T2 (podpróbka nieobciążona działa), ale żadna zmienna zdarzenia nie domyka"
    else:
        vj = "ODRZUCONA: ani piedestał, ani podpróbka away-side, ani f_n nie zdejmują deficytu → ε z danych per strumień (T1 lit.) albo embedding (T7 lit.)"
    K["joint"] = vj
    R["K"] = K
    json.dump(R, open(os.path.join(a.out, "trkjet2.json"), "w"), indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))

    # ================================================================================================================ raport
    L = []
    w = L.append
    header(L, partial, a.indir, a.out)
    w("\n## 0. Werdykty\n")
    w("| test | werdykt |\n|---|---|")
    w(f"| T1 stożki prostopadłe / losowe | **{K['T1']}** |")
    w(f"| T2 podpróbka nieobciążona trygerem | **{K['T2']}** |")
    w(f"| T3 zastępnik energii neutralnej (spekulatywne) | **{K['T3']}** |")
    w(f"| łącznie: „deficyt = obciążenie otoczenia D⁰ usuwalne zmienną zdarzenia” | **{vj}** |")
    w("\n| kryterium | liczby | wynik |\n|---|---|---|")
    c = K["K1.1"]["cells"]
    w("| K1.1 SF_perp = ε(perp1 < 2 GeV) dane/MC | " + "; ".join(f"{n}: ε {num(v['eps_data'], 3)} / {num(v['eps_mc'], 3)} → {pm(v['SF_perp'], v['err'])} ({v['method']})" for n, v in c.items()) +
      f" | piedestał podniesiony w J12: {tak(K['K1.1']['pass'])}; kontrola 9:MB zgodna z 1: {tak(K['K1.1']['sanity_MB'])} |")
    k = K["K1.2"]
    w(f"| K1.2 udział piedestału F_ped (J12only) | nadwyżka dane − MC: stożek D⁰ {pm(k['excess_cone'], k['excess_cone_err'], 2)} GeV, stożek prostopadły "
      f"{pm(k['excess_perp'], k['excess_perp_err'], 2)} GeV; F_ped = {pm(k['F_ped'], k['err'], 2)} | {k['class'] or 'nieokreślone (mianownik < 3σ)'} |")
    k = K["K1.3"]

    def cstr(c_):
        return "–" if not c_ else (f"„{c_['scheme']}”: 9:MB {pm(c_['SF_MB9'], c_['err_MB9'])}, J12only {pm(c_['SF_J12'], c_['err_J12'])}, "
                                   f"χ²/ndf {num(c_['chi2'], 1)}/{c_['ndf']}, cof. {num(c_['fallback_J12'], 2)}")
    w(f"| K1.3 domknięcie | odjęcie piedestału: {cstr(k['subtraction'])}; najlepsze przeważenie: {cstr(k['best_rw'])}; odniesienie: {cstr(k['ref'])} | "
      f"{tak(k['pass'])}{' (częściowo)' if k['partial'] else ''} |")
    k = K["K2.1"]; p_ = k["primary"]
    w("| K2.1 SF w podpróbce T(X, Y = 5; weto O) | " + (f"punkt główny X = {p_['X']} GeV: N_sig {pm(p_['N_sig'], p_['N_err'], 0)}, n_MC {p_['n_mc']}, ε {pm(p_['eps'], p_['eps_err'])} / "
      f"{num(p_['eps_mc'], 3)} → SF_T {pm(p_['SF'], p_['err'])} (po pT: {pm(p_['SF_rw_pt'], p_['err_rw_pt'])}); inkluzywny {pm(k['SF_incl'], k['SF_incl_err'])}; "
      f"różnica {pm(k['diff'], k['diff_err'])}" if p_ else f"brak punktu z N_sig ≥ {NFIT} i n_MC ≥ {NMC_OK}; inkluzywny {pm(k['SF_incl'], k['SF_incl_err'])}") +
      f" | {tak(k.get('pass')) if p_ else 'nierozstrzygnięte'}{' (częściowo)' if k.get('partial') else ''} |")
    k = K["K2.2"]
    w(f"| K2.2 trend SF_T(X) | punktów miarodajnych {k['n_points']}; ostatni − pierwszy: {num(k['slope_sigma'], 1)}σ | rośnie: {tak(k['rising']) if k['n_points'] >= 2 else '–'} (informacyjnie) |")
    k = K["K2.3"]

    def bstr(n):
        v = k["cells"].get(n)
        return f"{n}: " + ("–" if not v else f"N_sig {num(v['N_sig'], 0)}, SF {pm(v['SF'], v['err'], 2)}" + ("" if v["q"] == "ok" else f" ({v['q']})"))
    w(f"| K2.3 bity ortogonalne | (i) {bstr('8:softMB ∧ TBP_J12=0')}; (ii) {bstr('8:softMB ∧ TBP_J12=1')}; {bstr('9:MB ∧ TBP_J12=1')}; "
      f"ułamek J12only z TBP_MBTS_1: {pm(k['tbp_mbts_frac'], k['tbp_mbts_frac_err'])} | (i) warunki runu w porządku: {tak(k['i'])}; (ii) deficyt idzie za dżetem L1: {tak(k['ii'])} |")
    k = K["K3.0"]
    w(f"| K3.0 MC nadaje się do przeważania w r_n wprost | ⟨r_n⟩_MC / ⟨r_n⟩_9:MB = {num(k['ratio_mc_over_mb'], 2)} | {tak(k['pass'])} |")
    k = K["K3.1"]
    w(f"| K3.1 ⟨r_n⟩ J12only vs 8:notJ12 | {pm(k['J12only'], k['J12only_err'], 4)} vs {pm(k['notJ12'], k['notJ12_err'], 4)}; różnica {num(k['nsigma'], 1)}σ | {tak(k['pass'])} |")
    k = K["K3.2"]
    w(f"| K3.2 SF po przeważeniu „{k['scheme']}” | J12only: {pm(k['J12']['SF'], k['J12']['err']) if k['J12'] else '–'} (sam pT: {pm(k['base']['SF'], k['base']['err']) if k['base'] else '–'}; "
      f"ΔSF = {num(k.get('dSF'), 3)}); 9:MB: {pm(k['MB9']['SF'], k['MB9']['err']) if k['MB9'] else '–'} | zmienia: {tak(k.get('changes'))}; domyka: {tak(k.get('closes'))} |")

    w("\n## 1. Kryteria (zapisane przed obejrzeniem danych)\n")
    w(KRYTERIA)

    w("\n## 2. Wejście i metoda\n")
    w("| próbka | par | plików | komplet córek | stożek losowy znaleziony | ξ_cal brak |\n|---|---|---|---|---|---|")
    for nm, v in R["inputs"].items():
        w(f"| {nm} | {v['n_pairs']} | {v['n_files']} | {num(v['frac_daughters_ok'], 4)} | {num(v['frac_rnd_ok'], 4)} | {num(v['frac_xi_nan'], 4)} |")
    w(f"\n- Pary, selekcja, klasy, fit i SBS jak w TRKJET (`../TRKJET/RAPORT.md` §1). MC: fd_nd, truth ∈ {{0,1}}, N = {R['mc_n']}; „MC przeważone w pT” = wagi w 12 binach pT "
      "do widma sygnału danej klasy.")
    w("- Stożki prostopadłe i losowy: te same tory i ta sama definicja co własny stożek D⁰ (`v_s40_p2`: wszystkie tory drzewa `Tracks`, pT > 0,2 GeV, bez córek D⁰; "
      "zgodność własnego stożka z gałęzią `ptcone40`: ta sama decyzja < 2 GeV w ≥ 99 %, TRKJET §1). Stożek losowy: pierwsza z 16 liczb losowych φ spełniająca warunki; "
      "gdy żadna — użyty stożek prostopadły.")
    w("- Punkt wyjścia (odtworzenie):\n")
    w("| klasa | N_sig (SBS) | ε SBS | ε fit | ε_MC | SF | jakość |\n|---|---|---|---|---|---|---|")
    for name in CLS:
        r = R["cls"][name]; sf, sfe, meth = sf_cls(r); f = r.get("fit")
        w(f"| {name} | {pm(r['S'], r['S_err'], 0)} | {pm(r['eps'], r['err'])} | {pm(f['eps'], f['err']) if f else '–'} | {num(r['mc']['eps'], 3)} | {pm(sf, sfe, 2)} ({meth}) | {quality(r)} |")

    w("\n## 3. T1 — stożki prostopadłe / losowe\n")
    w("**3.1 Średnia Σ pT sygnału [GeV]: stożek D⁰, stożki prostopadłe, losowy (MC przeważone w pT danej klasy)**\n")
    w("| klasa | stożek D⁰: dane | MC | perp_avg: dane | MC | perp +π/2 | perp −π/2 | losowy: dane | MC | nadwyżka stożek D⁰ | nadwyżka perp | F_ped |\n|" + "---|" * 12)
    for name, o in R["t1"]["means"].items():
        w(f"| {name} | {pm(o['cone']['data'], o['cone']['data_err'], 2)} | {num(o['cone']['mc'])} | {pm(o['perp_avg']['data'], o['perp_avg']['data_err'], 2)} | {num(o['perp_avg']['mc'])} | "
          f"{num(o['p_s2_p']['data'])} | {num(o['p_s2_m']['data'])} | {pm(o['rnd']['data'], o['rnd']['data_err'], 2)} | {num(o['rnd']['mc'])} | "
          f"{pm(o['excess_cone'], o['excess_cone_err'], 2)} | {pm(o['excess_perp'], o['excess_perp_err'], 2)} | {pm(o.get('F_ped'), o.get('F_ped_err'), 2)} |")
    w("\n**3.2 ε(perp1 < próg) — „izolacja” stożka prostopadłego, dane vs MC (przeważone w pT)**\n")
    w("| klasa | próg [GeV] | N_sig | ε dane | metoda | ε MC | SF_perp | jakość |\n|---|---|---|---|---|---|---|---|")
    for key, r in R["t1"]["perp"].items():
        name, cut = key.split("|"); e, er, meth = best(r); sf, sfe, _ = sf_cls(r)
        w(f"| {name} | {cut} | {pm(r['S'], r['S_err'], 0)} | {pm(e, er)} | {meth} | {num(r['mc']['eps'], 3)} | {pm(sf, sfe)} | {quality(r)} |")
    w("\n**3.3 SF izolacji: nominalna, własny stożek, po odjęciu piedestału (fit pass/fail; MC przeważone w pT)**\n")
    w("| definicja pass | klasa | N_sig | ε dane | metoda | ε MC (pT) | ε MC (bez wag) | SF (MC pT) | jakość |\n|---|---|---|---|---|---|---|---|---|")
    for key, r in R["t1"]["sub"].items():
        tag, name = key.split("|"); e, er, meth = best(r); sf, sfe, _ = sf_cls(r)
        w(f"| {tag} < 2 GeV | {name} | {pm(r['S'], r['S_err'], 0)} | {pm(e, er)} | {meth} | {num(r['mc']['eps'], 3)} | {num(r['mc_unw']['eps'], 3)} | {pm(sf, sfe)} | {quality(r)} |")
    w("\n**3.4 SF po przeważeniu MC (izolacja nominalna; „[odjęcie]” = warunek własny − perp_avg < 2 GeV)**\n")
    w("| schemat | " + " | ".join(CLS) + " | χ²/ndf vs 1 (4 klasy rozłączne) | cof. J12only |\n|---|" + "---|" * (len(CLS) + 2))
    for sname in allsch:
        row = []
        for n in CLS:
            x = R["t1"]["rw"].get(f"{sname}|{n}")
            row.append(pm(x["SF"], x["err"]) + ("" if x["S"] >= NWEAK else " *") if x else "–")
        ch = R["t1"]["rw"][f"{sname}|chi2"]; xj = R["t1"]["rw"].get(f"{sname}|8:J12only")
        w(f"| {sname} | " + " | ".join(row) + f" | {num(ch['chi2'], 1)}/{ch['ndf']} | {num(xj['frac_fallback'], 2) if xj else '–'} |")
    w("\n\\* N_sig < 150. Klasy rozłączne do χ²: " + ", ".join(DISJ) + ".")
    w("\n**3.5 Dane→dane: ε(klasa) / ε przewidziane z 9:MB w komórkach (bez MC)**\n")
    w("| komórki | klasa | ε dane | ε przewidziane z 9:MB | stosunek | cof. |\n|---|---|---|---|---|---|")
    for key, r in R["t1"]["dd"].items():
        sname, name = key.split("|")
        w(f"| {sname} | {name} | {num(r['eps_data'], 3)} | {num(r['eps_pred'], 3)} | {pm(r['SF'], r['err'])} | {num(r['frac_fallback'], 2)} |")
    d = R["t1"]["diff"]; ed = d["edges"]
    lab = [f"{lo:g}–{hi:g}" if hi < INF else f"≥ {lo:g}" for lo, hi in zip(ed[:-1], ed[1:])]
    w("\n**3.6 ε_iso w binach perp_avg [GeV]**\n")
    w("| próbka | " + " | ".join(lab) + " |\n|---|" + "---|" * len(lab))
    for n in DISJ:
        w(f"| {n} | " + " | ".join((pm(r["eps"], r["err"], 2) + f" ({r['S']:.0f})") if np.isfinite(r["eps"]) and r["S"] > 25 else "–" for r in d["cls"][n]) + " |")
    w("| MC fd_nd | " + " | ".join((pm(r["eps"], r["err"], 2) + f" ({r['S']:.0f})") if r["S"] > 25 else "–" for r in d["MC"]) + " |")
    w("\nχ²/ndf: " + "; ".join(f"{k_.replace('|', ' vs ')}: {num(c_['chi2'], 1)}/{c_['ndf']}" + (f" (p = {c_['p']:.2g})" if c_["p"] is not None else "") for k_, c_ in d["chi2"].items()) + ".")

    w("\n## 4. T2 — podpróbka nieobciążona trygerem\n")
    w("**4.1 Skan T(X, Y)** — X: min pT dżetu torowego w |Δφ| > 2π/3 (0 = bez wymagania); Y: weto na dżet w półkuli D⁰ („-” = bez weta). SF = ε(dane)/ε_MC(T) bez wag; „SF po pT” = MC przeważone w pT wewnątrz T.\n")
    w("| weto | X | Y | klasa | N_sig | ε dane | metoda | n_MC | ε MC | SF | SF po pT | uwagi |\n|" + "---|" * 12)
    for key, r in R["t2"]["grid"].items():
        vname, X, Y, name = key.split("|")
        if vname != "O" and name != "8:J12only":
            continue
        e, er, meth = best(r); sf, sfe, _ = sf_cls(r)
        nmc = int(r["mc"]["S"]); rw = r.get("rw_pt") or {}
        notes = [q for q in (quality(r) if quality(r) != "ok" else "", "MC słabe" if NMC_MIN <= nmc < NMC_OK else "", "brak MC" if nmc < NMC_MIN else "") if q]
        w(f"| {vname} | {X} | {Y} | {name} | {pm(r['S'], r['S_err'], 0)} | {pm(e, er)} | {meth} | {nmc} | {num(r['mc']['eps'], 3)} | "
          f"{pm(sf, sfe) if nmc >= NMC_MIN else '–'} | {pm(rw.get('SF'), rw.get('err'))} | {', '.join(notes)} |")
    lab = [f"{lo:g}–{hi:g}" if hi < INF else f"≥ {lo:g}" for lo, hi in zip(AWAY_E[:-1], AWAY_E[1:])]
    w("\n**4.2 ε_iso w binach pT dżetu strony przeciwnej [GeV] (weto O, Y = 5 GeV)**\n")
    w("| próbka | " + " | ".join(lab) + " |\n|---|" + "---|" * len(lab))
    for n in T2_CLS:
        w(f"| {n} | " + " | ".join((pm(r["eps"], r["err"], 2) + f" ({r['S']:.0f})") if np.isfinite(r["eps"]) and r["S"] > 25 else "–" for r in R["t2"]["away"][n]) + " |")
    w("| MC fd_nd | " + " | ".join((pm(r["eps"], r["err"], 2) + f" ({r['S']:.0f})") if r["S"] > 25 else "–" for r in R["t2"]["away"]["MC"]) + " |")
    hv = R["t2"]["have_tbp"]
    w(f"\n**4.3 Bity ortogonalne** (kolumny `x_tbp_L1_J12`: 428770 {'jest' if hv['8'] else 'BRAK'}, 427929 {'jest' if hv['9'] else 'BRAK'}; `x_tbp_L1_MBTS_1`: {'jest' if hv['8_mbts'] else 'BRAK'}). "
      "softMB = " + ", ".join(f"`{c}`" for c in SOFT_MB) + ". SF względem MC inkluzywnego (bez wag) i przeważonego w pT.\n")
    w("| klasa | N_sig | ε dane | metoda | SF (bez wag) | SF po pT | jakość |\n|---|---|---|---|---|---|---|")
    for xn, r in R["t2"]["bits"].items():
        e, er, meth = best(r); sf, sfe, _ = sf_cls(r); rw = r.get("rw_pt") or {}
        w(f"| {xn} | {pm(r['S'], r['S_err'], 0)} | {pm(e, er)} | {meth} | {pm(sf, sfe, 2)} | {pm(rw.get('SF'), rw.get('err'), 2)} | {quality(r)} |")

    w("\n## 5. T3 — zastępnik energii neutralnej (spekulatywne)\n")
    w("**5.1 Średnie sygnału** (r_n = ξ_trk/ξ_cal; ξ_cal = A + C; ξ_trk z torów pT > 0,2 GeV)\n")
    w("| próbka | ułamek z ξ_cal > 0 | ⟨ξ_cal⟩ | ⟨ξ_trk⟩ | ⟨r_n⟩ |\n|---|---|---|---|---|")
    for n, o in R["t3"]["means"].items():
        w(f"| {n} | {num(o['frac_ok'], 3)} | {pm(o['xi_cal']['data'], o['xi_cal']['data_err'], 4)} | {pm(o['xi_trk']['data'], o['xi_trk']['data_err'], 4)} | {pm(o['rn']['data'], o['rn']['data_err'], 4)} |")
    w("\n**5.2 SF po przeważeniu MC** (u = kwantyl r_n względem klasy odniesienia tego samego runu; MC: własny)\n")
    w("| schemat | " + " | ".join(CLS) + " | cof. J12only |\n|---|" + "---|" * (len(CLS) + 1))
    for sname, _ in T3_SCHEMES:
        row = []
        for n in CLS:
            x = R["t3"]["rw"].get(f"{sname}|{n}")
            row.append(pm(x["SF"], x["err"]) + ("" if x["S"] >= NWEAK else " *") if x else "–")
        xj = R["t3"]["rw"].get(f"{sname}|8:J12only")
        w(f"| {sname} | " + " | ".join(row) + f" | {num(xj['frac_fallback'], 2) if xj else '–'} |")
    w("\n**5.3 Dane→dane w 428770: ε(J12only) / ε przewidziane z 8:notJ12**\n")
    w("| komórki | ε J12only | ε przewidziane | stosunek | cof. |\n|---|---|---|---|---|")
    for key, r in R["t3"]["dd"].items():
        w(f"| {key.split('|')[0]} | {num(r['eps_data'], 3)} | {num(r['eps_pred'], 3)} | {pm(r['SF'], r['err'])} | {num(r['frac_fallback'], 2)} |")
    lab = [f"{lo:g}–{min(hi, 1):g}" for lo, hi in zip(U_E[:-1], U_E[1:])]
    w("\n**5.4 ε_iso w binach u (kwantyl r_n; małe u = mało torów względem kalorymetru = duża frakcja neutralna)**\n")
    w("| próbka | " + " | ".join(lab) + " | χ² stałej |\n|---|" + "---|" * (len(lab) + 1))
    for n, rows in R["t3"]["diff"]["cls"].items():
        fl = R["t3"]["diff"]["flat"].get(n)
        w(f"| {n} | " + " | ".join((pm(r["eps"], r["err"], 2) + f" ({r['S']:.0f})") if np.isfinite(r["eps"]) and r["S"] > 25 else "–" for r in rows) +
          f" | {(num(fl['chi2'], 1) + '/' + str(fl['ndf'])) if fl else '–'} |")
    w("| MC fd_nd | " + " | ".join((pm(r["eps"], r["err"], 2) + f" ({r['S']:.0f})") if r["S"] > 25 else "–" for r in R["t3"]["diff"]["MC"]) + " | |")

    w("\n## 6. Zastrzeżenia\n")
    w("- Stożki prostopadłe mierzą aktywność przy tej samej η, 90° od D⁰; w zdarzeniach dwudżetowych (D⁰ w dżecie, drugi dżet naprzeciw) to obszar poprzeczny — "
      "piedestał „czysty”. Stożek losowy może trafić w stronę przeciwną (dżet odrzutu) → zwykle ≥ perp.")
    w("- T2: dżety torowe ≠ dżety kalorymetryczne L1; „twardy dżet torowy naprzeciw” nie gwarantuje, że to on odpalił L1_J12 (frakcja neutralna, próg ~12 GeV EM). "
      "MC minbias ma mało takich topologii — patrz n_MC; SF przy n_MC < 100 nie wchodzi do werdyktu.")
    w("- T2 weto B [koło]: brak dżetu przy D⁰ zbudowanego m.in. z torów stożka bezpośrednio podnosi ε — w danych i w MC, ale niekoniecznie tak samo.")
    w("- T2 K2.3: trigTBP_L1_MBTS_* nie wybiera próbki nieobciążonej (zdarzenie i tak zapisał J12); nieobciążone są tylko zdarzenia zapisane łańcuchem minbias.")
    w("- T3: ξ_cal w MC ma inną skalę niż w danych (test ROOT 17.11: mediana ξ_cal_A 0,0045 w fd_nd vs 0,05 w 427929 i 0,14 w 428770), a w 428770 zawiera pile-up "
      "→ r_n nie jest czystą frakcją naładowaną; schemat kwantylowy zakłada, że ranga r_n w MC odpowiada randze w klasie odniesienia; 8:notJ12 nie jest próbką minbias. "
      "Wynik T3 w obie strony jest wskazówką, nie dowodem.")
    w("- Przeważanie: komórki z n_MC < 25 cofają się do mniejszej liczby wymiarów („cof.”); schemat z cof. > 0,30 w J12only nie jest testem.")
    w("- SBS zakłada kształt tła Δm niezależny od zmiennej wewnątrz grupy (bin pT × pass/fail).")
    open(a.report, "w").write("\n".join(L) + "\n")
    print(f"T1: {K['T1']}\nT2: {K['T2']}\nT3: {K['T3']}\nłącznie: {vj}\nraport: {a.report}  ({time.time() - t0:.0f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

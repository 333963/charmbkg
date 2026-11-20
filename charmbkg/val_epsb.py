"""B4 — estymatory ε_B pod pikiem (SR) dla ``charmtag validate`` (projekt A).

Addytywne, domyślnie wyłączone (``validate.eps_b_nominal: legacy``). Włączenie: ``validate.eps_b_nominal: pol2``
(nominał = ε_B(SR) z fitu exp(pol2) do SB danych ekstrapolowanego do okna — S8 ``data_eff``, VAL ``pol2_closure``)
albo ``salad`` (SALAD v1.0, ``bkggen/salad_v1``, wagi w na MC w SR; gdy brak wag → fallback ``pol2`` z adnotacją).

Proxy SB: **symetryczne** (obie strony, |m−m0| ∈ [sb_min, sb_max]) zamiast asymetrycznych (1,70–1,80 ∪ 1,94–2,20 dla D⁺),
bo ε_B tła D⁺ spada z m ×2,7 (X/S8/VAL): asym zaniża ε_B pod pikiem ~20 %, sym zawyża 6–14 %, pol2/SALAD trafia w 3–10 %.
Stare pola ``eps_regions`` (SB_asym, lewy, prawy…) nie są ruszane — ``nominal_block`` dopisuje obok nich blok
``eps_B_SR`` z polami ``nominal``/``method`` i stosunkami proxy/nominał oraz adnotacją ``legacy``.

Wszystko w GeV (jak ``charmbkg.constants``). Nie importuje nic z charmbkg poza numpy — można wołać z ``pipeline.cmd_validate``
albo samodzielnie na tablicach (m, maska pass per WP).
"""

from __future__ import annotations

import math

import numpy as np

LEGACY_NOTE = (
    "legacy: eps_regions.SB_asym / lewy / prawy = stare proxy (asymetryczne SB); dla D+ zanizone ~20 % wzgledem eps_B pod pikiem "
    "(S8/VAL). Nominal B4: pol2 (fit exp(pol2) SB->okno) lub SALAD v1.0; proxy = SB symetryczne."
)

LEGEND = {
    "nominal": "eps_B(SR) wg `method`: pol2 = N_pass_pred(SR)/N_all_pred(SR) z fitu exp(pol2) do SB danych "
    "(sygnal poza fitem); salad = sum(w*pass)/sum(w) na MC w SR z wagami SALAD v1.0",
    "sb_sym": "eps_B w symetrycznych SB |m-m0| in [sb_min, sb_max] (obie strony)",
    "sb_asym_legacy": "eps_B w asymetrycznych SB (stare proxy charmbkg) — tylko do porownania",
    "rel_sym_vs_nom": "sb_sym / nominal (S8: 1.06-1.14 dla D+)",
    "rel_asym_vs_nom": "sb_asym_legacy / nominal (S8: ~0.78-0.83 dla D+)",
    "stat_rel": "przyblizony blad wzgledny nominalu (Poisson na przewidywanych zliczeniach SB)",
}


# ----------------------------------------------------------------------------------------------------------------------
# fit exp(pol2) do SB i ekstrapolacja do okna
# ----------------------------------------------------------------------------------------------------------------------
def _hist(m: np.ndarray, lo: float, hi: float, nbins: int):
    counts, edges = np.histogram(m, bins=nbins, range=(lo, hi))
    centers = 0.5 * (edges[1:] + edges[:-1])
    return counts.astype(float), centers, edges


def fit_exp_pol2_sb(
    m: np.ndarray,
    m0: float,
    sr_half: float,
    lo: float,
    hi: float,
    nbins: int = 100,
    sb_exclude_half: float | None = None,
) -> dict:
    """Fit log N(m) = a + b t + c t² (t = m − m0) na binach SB (|m−m0| ≥ sb_exclude_half, domyślnie = sr_half),
    wagami Poissona; zwraca przewidywaną liczbę zdarzeń w oknie |m−m0| < sr_half oraz parametry.

    Zwraca dict: n_pred_sr, n_sb, coef (a, b, c), chi2_ndf (na binach SB), ok.
    """
    excl = sr_half if sb_exclude_half is None else sb_exclude_half
    counts, centers, _edges = _hist(np.asarray(m, float), lo, hi, nbins)
    t = centers - m0
    sb = np.abs(t) >= excl
    use = sb & (counts > 0)
    if use.sum() < 6:
        return {"n_pred_sr": float("nan"), "n_sb": float(counts[sb].sum()), "coef": None, "chi2_ndf": float("nan"), "ok": False}
    y = np.log(counts[use])
    w = np.sqrt(counts[use])  # σ(log n) = 1/√n → waga = √n
    coef = np.polyfit(t[use], y, 2, w=w)  # [c, b, a]
    pred = np.exp(np.polyval(coef, t))
    sr = np.abs(t) < sr_half
    resid = (counts[use] - pred[use]) ** 2 / np.maximum(pred[use], 1.0)
    ndf = max(int(use.sum()) - 3, 1)
    return {
        "n_pred_sr": float(pred[sr].sum()),
        "n_sb": float(counts[sb].sum()),
        "coef": [float(coef[2]), float(coef[1]), float(coef[0])],
        "chi2_ndf": float(resid.sum() / ndf),
        "ok": bool(np.isfinite(pred).all() and resid.sum() / ndf < 5.0),
    }


def eps_b_sr_pol2(
    m: np.ndarray,
    passed: np.ndarray,
    m0: float,
    sr_half: float,
    lo: float,
    hi: float,
    nbins: int = 100,
    sb_exclude_half: float | None = None,
) -> dict:
    """ε_B(SR) = N_pass_pred(SR) / N_all_pred(SR); oba z fitu exp(pol2) do SB (po i przed cięciem). Sygnał poza fitem."""
    m = np.asarray(m, float)
    passed = np.asarray(passed, bool)
    f_all = fit_exp_pol2_sb(m, m0, sr_half, lo, hi, nbins, sb_exclude_half)
    f_pass = fit_exp_pol2_sb(m[passed], m0, sr_half, lo, hi, nbins, sb_exclude_half)
    ok = f_all["ok"] and f_pass["ok"] and f_all["n_pred_sr"] > 0
    eps = f_pass["n_pred_sr"] / f_all["n_pred_sr"] if ok else float("nan")
    # przybliżony błąd: Poisson na zliczeniach SB użytych w fitach (dominuje strona „pass”)
    stat_rel = math.sqrt(1.0 / max(f_pass["n_sb"], 1.0) + 1.0 / max(f_all["n_sb"], 1.0)) if ok else float("nan")
    return {
        "eps": float(eps),
        "stat_rel": float(stat_rel),
        "n_pred_sr_all": f_all["n_pred_sr"],
        "n_pred_sr_pass": f_pass["n_pred_sr"],
        "chi2_ndf_all": f_all["chi2_ndf"],
        "chi2_ndf_pass": f_pass["chi2_ndf"],
        "ok": bool(ok),
    }


# ----------------------------------------------------------------------------------------------------------------------
# proxy SB: symetryczne / asymetryczne (legacy)
# ----------------------------------------------------------------------------------------------------------------------
def eps_b_sb_sym(m: np.ndarray, passed: np.ndarray, m0: float, sb_min: float, sb_max: float) -> dict:
    d = np.abs(np.asarray(m, float) - m0)
    sel = (d >= sb_min) & (d <= sb_max)
    n = int(sel.sum())
    k = int((sel & np.asarray(passed, bool)).sum())
    return {"eps": k / n if n else float("nan"), "n": n, "n_pass": k}


def eps_b_sb_asym(m: np.ndarray, passed: np.ndarray, lo: float, left_max: float, right_min: float, hi: float) -> dict:
    mm = np.asarray(m, float)
    sel = ((mm >= lo) & (mm < left_max)) | ((mm > right_min) & (mm <= hi))
    n = int(sel.sum())
    k = int((sel & np.asarray(passed, bool)).sum())
    return {"eps": k / n if n else float("nan"), "n": n, "n_pass": k}


def eps_b_salad(mc_sr_mask: np.ndarray, mc_passed: np.ndarray, w: np.ndarray) -> dict:
    """ε_B(SR) z MC ważonego SALAD v1.0: Σ w·pass / Σ w w SR."""
    sr = np.asarray(mc_sr_mask, bool)
    w = np.asarray(w, float)
    p = np.asarray(mc_passed, bool)
    sw = w[sr].sum()
    swp = w[sr & p].sum()
    n_eff = (w[sr].sum() ** 2 / max((w[sr] ** 2).sum(), 1e-12)) if sr.any() else 0.0
    return {"eps": swp / sw if sw > 0 else float("nan"), "sum_w": float(sw), "n_eff": float(n_eff)}


# ----------------------------------------------------------------------------------------------------------------------
# blok do validation/<sp>.json
# ----------------------------------------------------------------------------------------------------------------------
def nominal_block(
    m_data: np.ndarray,
    passed_by_wp: dict,
    m0: float,
    regions: dict,
    method: str = "pol2",
    salad: dict | None = None,
    nbins: int = 100,
) -> dict:
    """Blok ``eps_B_SR`` per WP. ``regions``: {lo, hi, sr_half, sb_min, sb_max, left_max, right_min} w GeV
    (D⁺: lo 1.70, hi 2.20, sr_half 0.045…0.055, sb_min 0.070, sb_max 0.170, left_max 1.80, right_min 1.94).
    ``salad``: {wp: {"mc_sr": mask, "mc_passed": mask, "w": wagi}} — gdy ``method == 'salad'``.
    """
    out = {"method": method, "legacy_note": LEGACY_NOTE, "legend": LEGEND, "wps": {}}
    r = regions
    for wp, passed in passed_by_wp.items():
        pol2 = eps_b_sr_pol2(m_data, passed, m0, r["sr_half"], r["lo"], r["hi"], nbins, r.get("sb_exclude_half"))
        sym = eps_b_sb_sym(m_data, passed, m0, r["sb_min"], r["sb_max"])
        asym = eps_b_sb_asym(m_data, passed, r["lo"], r["left_max"], r["right_min"], r["hi"])
        row = {
            "pol2": pol2,
            "sb_sym": sym,
            "sb_asym_legacy": asym,
        }
        used = "pol2"
        nominal = pol2["eps"]
        stat_rel = pol2["stat_rel"]
        if method == "salad":
            if salad and wp in salad:
                s = eps_b_salad(salad[wp]["mc_sr"], salad[wp]["mc_passed"], salad[wp]["w"])
                row["salad"] = s
                used = "salad"
                nominal = s["eps"]
                stat_rel = 1.0 / math.sqrt(max(s["n_eff"] * max(nominal, 1e-9), 1.0))
            else:
                row["salad"] = {"skipped": "brak wag SALAD dla tego WP — fallback pol2"}
        row["nominal"] = float(nominal)
        row["method_used"] = used
        row["stat_rel"] = float(stat_rel)
        row["rel_sym_vs_nom"] = float(sym["eps"] / nominal) if nominal and np.isfinite(nominal) else float("nan")
        row["rel_asym_vs_nom"] = float(asym["eps"] / nominal) if nominal and np.isfinite(nominal) else float("nan")
        out["wps"][str(wp)] = row
    return out


def default_regions(species: str) -> dict:
    """Okna jak w S8/B1 (GeV). D⁺: SR |Δm|<45–55 MeV → tu 0.050; SB sym 70–170 MeV; asym 1.70–1.80 ∪ 1.94–2.20.
    D0: SR 55 MeV, pasmo 55–75, SB >75 MeV (do 1.75–2.00); asym = 1.75–1.80 ∪ 1.93–2.00 (przybliżenie)."""
    if species.lower().startswith("dp"):
        return {"lo": 1.70, "hi": 2.20, "sr_half": 0.050, "sb_min": 0.070, "sb_max": 0.170, "left_max": 1.80, "right_min": 1.94}
    return {"lo": 1.75, "hi": 2.00, "sr_half": 0.055, "sb_min": 0.075, "sb_max": 0.125, "left_max": 1.80, "right_min": 1.93}

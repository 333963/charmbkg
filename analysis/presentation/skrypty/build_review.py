#!/usr/bin/env python3
"""Sklada PRZEGLAD.html i PRZEGLAD.md (ta sama tresc, po angielsku) z listy rysunkow. Zwykly lokalny HTML ze sciezkami wzglednymi; nic nie publikuje."""
import html
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TITLE = "Project B — final results in figures (open charm in ATLAS minimum-bias data, 2022, √s = 13.6 TeV)"
INTRO = ("Runs 427929 (clean minimum bias, L = 0.4648 ± 0.0103 nb⁻¹) and 428770 (jet-trigger dominated); "
         "D⁰→Kπ, D⁺→Kππ, D*±; classifier: tree H at WP70. All results are preliminary and post hoc. Every number on a figure is read by the script from the source file named in RAPORT.md. "
         "Layout of the 2D map (F4) follows the ROOT '_margins' plots: map, marginal distributions, nominal-mass lines and a residual panel.")
FIGS = [
    ("F1a", "fig/F1a_mass_fit_D0_M8.png", "D⁰ mass spectra with the nominal fit (model M8)",
     "The D⁰→Kπ peak is extracted in three pT bins of run 427929 with yields of 4059 ± 468, 1536 ± 161 and 700 ± 70 and acceptable fit quality (p = 0.24, 0.22, 0.07).",
     "Broad components tied to the signal (K↔π reflection, other charm decays) are taken from simulation; the fit has a second minimum at 3–5 and 5–8 GeV (4.6 % and 4.4 %, carried as a systematic). "
     "On simulation the fit underestimates the yield at 3–5 GeV by 22 ± 10 %; this is an uncertainty, not a correction."),
    ("F1b", "fig/F1b_mass_fit_Dplus_A2.png", "D⁺ mass spectra with the nominal fit (model A2)",
     "The D⁺→Kππ peak is clean above 5 GeV; yields are 3404 ± 292, 1630 ± 102 and 651 ± 53.",
     "Model A2 was chosen by ranking and post hoc; no setting meets the closure criterion on simulation at 5–8 GeV (+9.2 ± 6.2 %), carried as an uncertainty. "
     "At 3–5 GeV the fit quality is poor (p = 0.03) and the yield depends on the fitted width by up to 13 %."),
    ("F2", "fig/F2_cross_sections.png", "Inclusive cross-sections of D⁰, D⁺ and D*± per pT bin",
     "The cross-sections lie about 0.6 of Pythia 8 A3 above 5 GeV and inside or just above the FONLL band; D⁰ at 5–8 GeV is 3.0σ below CMS at 13 TeV after matching the definitions.",
     "Preliminary and post hoc; the tension with CMS is not explained. The ALICE comparison uses reference values that were not re-verified (8 % error assumed). "
     "Scale factors of D⁺ and of D⁰ at 3–5 GeV rest on a transfer from D⁰ above 5 GeV; tracking, vertex and luminosity-calibration uncertainties are taken from the literature or assumed."),
    ("F3", "fig/F3_ratios.png", "Ratios D⁺/D⁰ and D*⁺/D⁰",
     "D⁺/D⁰ = 0.50, 0.59, 0.49 and D*⁺/D⁰ = 0.48, 0.61 agree with all references within 1.4σ.",
     "The precision (24–33 %) does not separate 0.40 from 0.47. Both ratios sit 0.2–1.4σ above the measured references, which would correspond to a D⁰ cross-section low by 15–30 %."),
    ("F4", "fig/F4_map_mKpi_dm.png", "m(Kπ) × Δm map of D* candidates in data",
     "The D*± → D⁰π signal is visible as a correlated excess at the crossing of the nominal D⁰ mass and Δm = 145.4 MeV, without any classifier.",
     "Illustration only: no fit is drawn and no yield is taken from this map; σ(D*±) comes from a separate Δm fit. Selection as in that fit but without the D⁰ mass window; pT(D*) 5–30 GeV."),
    ("F5", "fig/F5_isolation_sf.png", "Isolation scale factor: minimum bias against the jet-triggered run",
     "The isolation efficiency of true D⁰ agrees with simulation in minimum-bias data (SF = 0.98 ± 0.04) but not in run 428770 (0.51 ± 0.02), so an MC-derived isolation efficiency is usable only for the minimum-bias sample.",
     "Measured on D⁰ tagged by a D*, whose soft pion lies inside the cone: on simulation the 427929 value carries over to inclusive D⁰ within 0.5 %, the 428770 one is model dependent. Statistical errors only; D⁺ not tested."),
    ("F6", "fig/F6_trigger.png", "Trigger: L1_J12 efficiency for events with a D and trigger regions of run 428770",
     "L1_J12 selects only 2–12 % of the events with a D and 70–87 % of the 428770 yield comes from the L1_J12-only region, so run 428770 gives no second absolute cross-section point.",
     "The efficiency is known only from run 427929 and with 16–71 % relative error; 12–24 % of the yield comes from events of none of the three chains, with unknown efficiency. "
     "D⁰ at 8–30 GeV is fitted with another model than in run 427929."),
    ("F7", "fig/F7_uncertainty_budget.png", "Uncertainty budget of σ(D⁰) and σ(D⁺)",
     "The measurement is systematics limited (totals 16–29 %): fit bias and fit model dominate for D⁰, the classifier scale factor at 5–8 GeV, and literature tracking uncertainties for D⁺.",
     "The fit bias is the closure on simulation (|b| ⊕ σ_b) and is not corrected. Tracking, impact-parameter, vertex and generator rows are from ATLAS publications, not measured on these data; the 2.2 % luminosity calibration is an assumption."),
    ("F8", "fig/F8_background_composition.png", "Background composition under the peak from MC truth",
     "At WP70 the D⁰ background is mostly charm-free combinatorics at low pT and becomes reflection- and charm-related at high pT, while the D⁺ background above 5 GeV is dominated by partially combinatorial charm.",
     "Simulation only (minimum-bias MC, unweighted), post hoc with respect to the fit-model choice; the data composition is tested only indirectly through the fits."),
    ("F9", "fig/F9_network.png", "Neural network: soft-pion effect and gain against the tree (sideband proxy and peak fit)",
     "Masking soft-pion candidates removes the D*-specific efficiency of the network (double ratio 1.37–1.45 → 1.014–1.017); its gain over the tree from the inclusive peak fit is then 0.99 / 1.10 / 1.09 in minimum-bias run 427929 (sideband proxy 1.04 / 1.05 / 1.17), below the pre-set 1.15, and 1.13 / 1.25 / 1.17 in the jet-triggered run 428770 (proxy 1.25 / 1.40 / 1.42).",
     "The sideband proxy overestimates the gain (ε_S from simulation). The peak-fit gain is from one training seed and has no error, so differences of about 0.05 are not significant. In run 428770 the gain sits in the L1_J12-only class and comes with a network yield, after the MC efficiency, 9–14 % below that of the tree; the closure of this yield ratio on simulation is not passed in any bin, so the ratio is descriptive. "
     "The D*-tag scale factor of the masked network is not stable between seeds at 5–8 GeV. The tree H stays nominal."),
]
TAG = [
    ("T1", "local_tagged/T1_tagged_D0_mass.png", "D⁰ mass spectrum in events with an AFP proton tag",
     "In run 427929 a D⁰ signal of 54.7 +19.4/−18.7 (stat) ± 12.6 (shape) candidates is indicated at Z = 3.05; a background shape fixed from the untagged sample gives the same result (54.1, Z = 3.05).",
     "An indication at the 3σ boundary, not an observation. The yield assumes that the broad components tied to the signal have the same proportions as inclusively; a counting check (window − sidebands) gives about 30. The variant is post hoc."),
    ("T2", "local_tagged/T2_tagged_R_sigma.png", "Tagged fraction R^fid and fiducial tagged cross-section against predictions",
     "R^fid = 0.99 +0.35/−0.34 ± 0.28 % and σ_tag^fid = 13.6 +4.8/−4.7 ± 4.6 µb are compatible with the Pythia envelopes and of the same order as the converted literature values.",
     "Post hoc with respect to the read-out of the tagged peaks. The window − sidebands variant (0.54 %, 7.4 µb) is a one-sided model systematic not included in ± syst; 64 % of the efficiency denominator comes from 3–5 GeV, where the classifier scale factor is a scenario. Literature values are order-of-magnitude conversions."),
]


def md():
    L = [f"# {TITLE}", "", INTRO, ""]
    for k, p, t, c, cav in FIGS:
        L += [f"## {k}. {t}", "", f"**{c}**", "", f"![{k}]({p})", "", f"*Caveats:* {cav}", ""]
    L += ["---", "", "# AFP-tagged measurement", "", "Measurement with the AFP proton tag (run 427929). Figures in `local_tagged/`.", ""]
    for k, p, t, c, cav in TAG:
        L += [f"## {k}. {t}", "", f"**{c}**", "", f"![{k}]({p})", "", f"*Caveats:* {cav}", ""]
    return "\n".join(L)


def htm():
    e = html.escape
    css = ("body{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;max-width:1180px;margin:2rem auto;padding:0 1rem;color:#0b0b0b;background:#fcfcfb;line-height:1.5}"
           "h1{font-size:1.5rem}h2{font-size:1.15rem;margin-top:2.6rem;border-top:1px solid #d9d8d2;padding-top:1rem}"
           "img{max-width:100%;height:auto;border:1px solid #d9d8d2;background:#fff}.c{font-weight:600}.cav{color:#52514e;font-size:.93rem}.intro{color:#52514e}"
           ".local{margin-top:4rem;padding:.6rem 1rem;border:2px solid #d03b3b;color:#d03b3b;font-size:1.4rem}a{color:#1c5cab}")
    H = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>{e(TITLE)}</title><style>{css}</style></head><body>",
         f"<h1>{e(TITLE)}</h1><p class='intro'>{e(INTRO)}</p>"]

    def block(k, p, t, c, cav):
        return (f"<h2 id='{k}'>{k}. {e(t)}</h2><p class='c'>{e(c)}</p><a href='{p}'><img src='{p}' alt='{e(k + ': ' + t)}'></a>"
                f"<p class='cav'><b>Caveats:</b> {e(cav)} (<a href='{p[:-4]}.pdf'>PDF</a>)</p>")
    H += [block(*f) for f in FIGS]
    H.append("<h1 class='local'>AFP-tagged measurement</h1><p class='intro'>Measurement with the AFP proton tag (run 427929). Figures in <code>local_tagged/</code>.</p>")
    H += [block(*f) for f in TAG]
    H.append("</body></html>")
    return "\n".join(H)


for k, p, *_ in FIGS + TAG:
    assert os.path.exists(os.path.join(BASE, p)) and os.path.exists(os.path.join(BASE, p[:-4] + ".pdf")), p
open(os.path.join(BASE, "PRZEGLAD.md"), "w").write(md())
open(os.path.join(BASE, "PRZEGLAD.html"), "w").write(htm())
print("zapisano PRZEGLAD.html, PRZEGLAD.md;", len(FIGS) + len(TAG), "rysunkow")

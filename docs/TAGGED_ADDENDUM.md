<!-- Not part of the methods note.
     Written 2026-11-20. Every number: reports under ~/Documents/charm/charm_work/ named in the src comments.
     Audit of 20.11 (charm_work/AUDYT_noty_2026-11-20): σ_tag^fid, R^fid and the budget recomputed by script with the final inclusive inputs (no new fit; N_tag unchanged). Backup: TAGGED_ADDENDUM.md.bak_przed_2011audyt.
     Closing update of 20.11 (charm_work/ZAMKNIECIE_2026-11-20): inclusive inputs of version v2 (larger error of the D⁰ scale factor above 5 GeV; quoted results unchanged), cross-check with a constructed background shape (Sec. 3), D⁺ remark. Backup: TAGGED_ADDENDUM.md.bak_przed_2011v2.
     Last additions of 20.11 (charm_work/ZAMKNIECIE_2026-11-20, last section): predictions recomputed after the closure of the Pythia variation and compared with the measurement (Sec. 5), with the statement on their independence; figures (Sec. 8). Backup: TAGGED_ADDENDUM.md.bak_przed_2011ost. -->

## Addendum: D mesons with a forward proton in AFP (run 427929)

This addendum collects the proton-tagged results of 19–20 November. They are deliberately kept out of `METHODS_NOTE_DRAFT.md`, whose Sec. 6.11 gives predictions only.

### 1 Status of each element

| element | status |
|---|---|
| tagged yields and tagged fraction (TAGGED, 19.11) | criteria pre-registered, written before the single unblinding of the tagged peaks (19.11, 23:51); the D⁺ result failed two of its control criteria and is descriptive |
| fiducial cross-section and R^fid (TAGXSEC, 20.11) | **post hoc**: criteria were written before the computation, but N_tag was already known, and the definition of the efficiency denominator was changed after a first pass (exclusive K∓π± mode); the values below are the arithmetic update of the audit of 20.11 with the final inclusive inputs |
| AFP acceptance and efficiency (AFP_akceptancja, 20.11) | separate study, does not use the tagged D yields |
| cross-check with a constructed background shape (KONSTR, 20.11) | **post hoc, after the reading of the tagged peaks**; a variant and a check, it does not replace the nominal |
| comparison with Pythia | post hoc; envelope taken from the generator study of 19.11 (range 2–30 GeV); the completion of that study on 20.11 leaves the envelope unchanged (× 0.40 to × 3.98) |
| predicted yields (SD_var, 19–20.11) | computed by a separate study that did not open the tagged result; two recorded incidents (Sec. 5) |
<!-- src: charm_work/TAGGED_2026-11-19/KRYTERIA.md, RAPORT.md; charm_work/TAGXSEC_2026-11-20/KRYTERIA.md, RAPORT.md; charm_work/AFP_akceptancja_2026-11-20/PLAN.md; charm_work/KONSTR_2026-11-20/KRYTERIA.md -->

### 2 Definitions

**Tag (detector level).** Exactly one reconstructed AFP proton with 0.030 ≤ ξ ≤ 0.100, on side A or C; pile-up vertex veto; ξ_cal ≤ ξ_AFP + 0.02. Run 427929, `HLT_mb_mbts_L1MBTS_1`, GRL, L = 0.4648 nb⁻¹; tree H at WP70.

**Fiducial definition (particle level).** D⁰ + anti-D⁰, prompt and non-prompt, 3 ≤ pT < 30 GeV, |η| < 2.2, in events with a proton of 0.030 ≤ ξ ≤ 0.100 on either side; no requirement on t or on the production process. The ξ_cal requirement and the vertex veto are corrected for and are not part of the definition.

σ_tag^fid = (N_tag − N_acc) / (C̄_tag · B · L · ⟨A·ε_AFP⟩ · ε_trig · ε_ξcal · ε_npu), and R^fid = σ_tag^fid / σ_incl in the same D⁰ phase space.
<!-- src: charm_work/TAGXSEC_2026-11-20/akapit_do_noty.md; RAPORT.md -->

### 3 Yields

| sample | species | N_tag (stat, syst) | significance | remark |
|---|---|---|---|---|
| 427929, MBTS | D⁰ | 54.7 −18.7 +19.4 ± 12.6 | Z = 3.05 (toys: p = 0.0015) | nominal: fit of the tagged peak |
| 427929, MBTS | D⁰ | 29.9 ± 12.5 | – | variant: signal window minus extrapolated sidebands, biased low by construction |
| 427929, MBTS | D⁺ | 29.1 −10.8 +11.5 ± 8.8 | Z = 2.9 | control criteria not met: descriptive only |
| 428770, R2 | D⁰ / D⁺ | 278 ± 49 ± 48 / 185 ± 30 ± 24 | Z = 6.0 / 6.5 | trigger-biased; no absolute normalisation (note, Sec. 6.11) |
| 428770, R1 / R3 | D⁰ | 0.9 ± 16 / 19.9 | – | no signal |

Tagged fraction of true D in run 427929: 0.88 −0.30 +0.31 ± 0.22 % (D⁰) and 0.50 −0.19 +0.20 ± 0.15 % (D⁺, descriptive). The D⁺ fraction was computed with the inclusive D⁺ yields of the old fit model; with the yields of the new model (note, Sec. 5.14: 3404 / 1630 / 651 instead of 3547 / 1651 / 687, a sum lower by 3.4 %) the denominator falls and the fraction would rise by about 3 % of its value, far inside its error. It was not recomputed, and D⁺ stays descriptive.
<!-- src: charm_work/TAGGED_2026-11-19/RAPORT.md; charm_work/TAGGED_2026-11-19/tabele; inclusive D⁺ yields: charm_work/FITBIAS_Dplus_2026-11-20/out/dplus_yields.json -->

**Cross-check with a constructed background shape (D⁰; post hoc, after the reading).** The combinatorial background shape of the tagged fit was fixed to that of the untagged sample of the same selection (2.1·10⁵ candidates), leaving only the signal and background normalisations free. Tagged events have a much lower track multiplicity (median n_sel^BL 18 against 38), but the sideband mass shapes agree without reweighting (χ²/ndf = 7.9/8), and reweighting in multiplicity and pT changes the yield by less than 0.5 events. On 3 × 200 pseudo-tagged samples drawn from the untagged candidates the procedure is unbiased (mean pull between −0.07 and −0.01, 68 % coverage 0.69–0.72), with an r.m.s. 3–6 % smaller than for the free background.

| fit of the tagged D⁰ peak | N_tag | Z |
|---|---|---|
| nominal (free background shape) | 54.7 −18.7 +19.4 | 3.05 |
| constructed background shape | 54.1 −18.3 +19.0 | 3.05 |
| constructed, fit range 1.75–2.00 GeV | 46.0 −19.7 +20.4 | 2.39 |
| constructed, broad signal-associated components moved into the background (three variants) | 44.3 / 46.6 / 52.8 | 2.45 / 3.10 / 3.02 |
| window minus extrapolated sidebands (variant of the table above) | 29.9 ± 12.5 | – |

Reading. (i) The free and the constructed background give the same yield and the same significance, so the freedom of the background shape was not the origin of the model spread. (ii) All models with the constructed background shape give 44–54. The value of about 30 from the window-minus-sidebands count is not reproduced by any of them, not even when all broad components are assigned to the background (52.8); it comes from the free extrapolation of the tagged sidebands alone, which puts 228 background candidates into the window against 187–222 in the fits. (iii) Consequence for the variant R^fid = 0.54 % of Sec. 4: it is kept as quoted, because it was defined before the reading and this check is post hoc, but it is no longer supported as an alternative model of the same standing as the nominal. It is the low edge given by an estimator known to be biased low, and the one-sided −46 % of Sec. 6 is a conservative envelope; within the constructed-background models the downward spread is 19 % (54 → 44). (iv) Shape systematic as a variant: with the constructed background it is 11.0 events (fit range 8.1, treatment of the broad components 7.5, reweighting 0.5), against the nominal 12.6 with the 30–55 spread not covered. The nominal yield, its systematic uncertainty of 12.6 and all numbers of Secs. 4 and 6 are unchanged. Not tested by this check: the assumption that the broad components have the same proportions in tagged events as inclusively (the check only measures its effect, 54 against 44–47), and the sideband shape test has little power.
<!-- src: charm_work/KONSTR_2026-11-20/NOTA_tekst_EN.md (part for Sec. 6); charm_work/KONSTR_2026-11-20/RAPORT.md §B; charm_work/KONSTR_2026-11-20/tabele/{B_prep,B_wyniki,B_extra}.md -->

### 4 Result (D⁰, run 427929)

| quantity | nominal | variant (window − sidebands) |
|---|---|---|
| R^fid = σ_tag^fid / σ_incl | **0.99 −0.34 +0.35 (stat) ± 0.29 (syst) %** | 0.54 ± 0.23 ± 0.09 % |
| σ_tag^fid | 13.9 −4.8 +5.0 (stat) ± 4.6 (syst) µb | 7.6 ± 3.2 ± 1.8 µb |

Inputs (audit of 20.11): the final inclusive values of the note, Sec. 6.13: σ_incl(D⁰, 3–30 GeV) = 1028 + 278 + 100 = 1406 µb, scale factors 0.967 ± 0.033 / 1.005 ± 0.135 / 0.965 ± 0.109 (version v2 of the inclusive computation: above 5 GeV the spread of the Δm estimators is added to the earlier ± 0.126 and ± 0.097; the central values, σ_incl and therefore σ_tag^fid and R^fid do not change, and their systematic uncertainties change only in the next digit, 4.61 → 4.62 µb), B = 3.936 % (without the doubly Cabibbo-suppressed mode, as in the denominator of C). With the inputs of the first computation (1375 µb, scale factor 1 ± 0.10 at 3–5 GeV, B = 3.951 %) the values were R^fid = 0.99 −0.34 +0.35 ± 0.28 % and σ_tag^fid = 13.6 −4.7 +4.8 ± 4.6 µb (variant 7.4 ± 3.1 ± 1.9 µb). R^fid does not change (0.986 → 0.988 %): B and L cancel exactly; the scale factor at 3–5 GeV moves numerator and denominator together (σ_tag^fid +2.1 %, residual in R −0.3 %); the new inclusive yields move σ_incl by −0.6 %. σ_tag^fid rises by 2.5 % (+2.1 % from the scale factor, +0.4 % from B).

R^fid is the main result: the branching fraction, the luminosity and the definition of the efficiency denominator cancel exactly, the classifier scale factor and the inclusive fit bias up to the difference of the tagged and inclusive pT spectra (≤ 1.3 %). The difference between the nominal and the variant is a separate one-sided model uncertainty (−46 %) and is not part of the systematic uncertainty quoted; after the cross-check of Sec. 3 it is to be read as a conservative envelope, since no fit with a background shape reproduces the variant.

Efficiency: C̄_tag = Σ_b g_b C_b SF_b ρ_b = 0.236 (0.241 with the scale factor 1 at 3–5 GeV), with C_b from the non-diffractive MC (exclusive-mode denominator), g_b the pT spectrum of D⁰ in single-diffractive MC with a proton in the window, and ρ = C(tagged)/C(ND) = 0.96 ± 0.12 for the lower multiplicity of tagged events (uncertainty = non-closure on single-diffractive MC, direct estimate 1.08 ± 0.09 from 98 reconstructed D⁰). About 63 % of C̄_tag comes from the 3–5 GeV bin.

Pythia 8 A3 (soft diffraction) gives R^fid = 0.35 % with an envelope of 0.14–1.40 %, a hard-diffraction scenario 0.66 % (0.58–1.41 %); both are compatible with the result (+1.4σ and +0.7σ nominal; +0.8σ and −0.5σ variant), which does not discriminate between them. For D⁺, R^fid = 0.55 −0.20 +0.22 ± 0.20 % is descriptive only; it was not recomputed with the new inclusive D⁺ yields (σ_incl(D⁺) lower by 4.0 / 1.3 / 5.2 %, which would raise it by a few per cent of its value).
<!-- src: charm_work/TAGXSEC_2026-11-20/akapit_do_noty.md; charm_work/TAGXSEC_2026-11-20/RAPORT.md; charm_work/TAGXSEC_2026-11-20/tabele/wyniki_excl.md, out/tagxsec_excl.json; Pythia envelope: charm_work/SD_var_2026-11-19/RAPORT.md; update with final inputs: charm_work/XSEC_final_2026-11-20/out/final.json, charm_work/AUDYT_noty_2026-11-20/skrypty/przeliczenia.py, out/przeliczenia.md §3; version v2: charm_work/XSEC_final_2026-11-20/v2/out/final.json, charm_work/ZAMKNIECIE_2026-11-20/out/przeliczenia.md §3 -->

### 5 AFP acceptance and efficiency

| factor | value | source of the number |
|---|---|---|
| geometric acceptance A in 0.030 ≤ ξ ≤ 0.100 | 0.97 ± 0.03 | acceptance study |
| reconstruction efficiency ε_AFP | 0.95 ± 0.02 | tag-and-probe: 0.968 (side A), 0.965 (side C) |
| ⟨A·ε_AFP⟩ | 0.92 ± 0.05 (0.922 used) | includes alignment and optics, ± 3.4 % |
| ε_trig (MBTS for events with a D) | 1.000 | note, Sec. 6.11 |
| ε_ξcal | 1.000 +0/−0.005 | |
| ε_npu (vertex veto) | 0.985 | note, Sec. 6.7 |

Efficiency of the L1 AFP trigger items (not used for the MBTS-triggered result): 0.904–0.915 on side A and 0.434–0.458 on side C. The factor 0.55 for "AFP acceptance and efficiency" in the earlier predictions of Sec. 6.11 of the note was an assumption and is superseded by 0.92.

**Measurement against prediction.** The predictions were recomputed on 20.11 after the completion of the Pythia variation (note, Sec. 6.11): N_tag = Σ over pT bins of the inclusive yield fitted in data × fraction, the fraction being R(pT) × 0.92 × 1.00 with the model envelope (generator scenario) or the tagged share of the sideband candidates with a factor 0.5–2 (sideband scenario). They replace the numbers of the audit of 20.11, which were the earlier predictions multiplied by 0.92/0.55 and by ratios of yields. The nominal row is WP70, the working point of the measurement; this choice was made on 20.11, after the reading of the tagged peaks (at WP80: D⁰ 26 (10–104) and 30 (15–60), D⁺ 21 (8–79) and 16 (8–32)); the fractions do not depend on the working point.

| sample | species | predicted, generator scenario (model envelope) | predicted, sideband scenario (× 0.5–2) | measured N_tag (stat, syst) | measured over predicted (generator / sidebands) |
|---|---|---|---|---|---|
| 427929, MBTS, WP70 | D⁰ | 20.4 (8.2–81) | 23.5 (12–47) | 54.7 −18.7 +19.4 ± 12.6 (variant 29.9 ± 12.5) | 2.7 / 2.3 (+1.5σ / +1.4σ) |
| 427929, MBTS, WP70 | D⁺ | 18.6 (7.1–71) | 14.2 (7.1–28) | 29.1 −10.8 +11.5 ± 8.8 (descriptive) | 1.6 / 2.0 (+0.8σ / +1.1σ) |
| 427929, MBTS | D* | 2.8 (1.1–11) | 2.1 (1.0–4) | not quoted in this addendum | – |
| 428770, R1, WP70 | D⁰ | 13 (5–52) | 15 (7–30) | 0.9 −16.0 +16.7 ± 11.1 | −0.6σ / −0.7σ |
| 428770, R1, WP70 | D⁺ | 12 (5–45) | 9 (5–18) | 24.8 −11.0 +11.6 ± 7.1 | 2.1 / 2.8 (+1.0σ / +1.2σ) |
| tagged fraction, 427929 | D⁰ | 0.32 % (0.13–1.28 %) | 0.37 % (0.19–0.74 %) | 0.88 −0.30 +0.31 ± 0.22 % | 2.7 / 2.4 |
| tagged fraction, 427929 | D⁺ | 0.32 % (0.12–1.20 %) | 0.24 % (0.12–0.48 %) | 0.50 −0.19 +0.20 ± 0.15 % (descriptive) | 1.6 / 2.1 |

Reading. The deviations in σ use the error of the measurement only (statistical on the side of the prediction, with the systematic one added in quadrature); the ranges of the predictions are spreads of models and assumptions, not confidence intervals, and predicted accidental tags (0.1–0.6) are negligible. The nominal D⁰ yield of run 427929 is above both central predictions by a factor of 2.3–2.7, which is 1.4–1.5σ; it lies inside the model envelope (8.2–81) and 0.3σ above the upper edge of the sideband range, and the window-minus-sidebands variant (29.9 ± 12.5) agrees with both central values. D⁺ is within 1.2σ of the predictions in both samples; D⁰ in R1 shows no signal and is compatible with its prediction. The comparison does not discriminate between the scenarios, and it does not support a rate above the model envelope. The predicted yields are not corrected for a possible difference of the classifier efficiency between diffractive and other events (ρ = 0.96 ± 0.12 of Sec. 4).

**Independence of the predictions (statement recorded by the prediction study).** The prediction study did not open the results of the tagged measurement (outputs, tables, figures, report); of the tagged study it read only the criteria file with the definitions. The method and the first numbers date from 19.11, 23:49 CET, two minutes before the unblinding of the tagged peaks; the final run is of 20.11, 02:01 CET. The changes in between come from three inputs outside the tagged region: the missing parts of the colour-reconnection variation (envelope changed by less than 0.3 %), the measured A·ε_AFP = 0.92 ± 0.05 in place of the assumed 0.55 (the acceptance study does not use the tagged D yields, but its plan records that its author had seen the title of the note with the tagged result), and the removal of a branching-fraction scale of the peaking component in the inclusive D⁰ fit (+0.1 % and −0.2 % on the prediction). Two incidents were recorded. (1) On 20.11 at about 00:05 CET the prediction study read a note of the inclusive cross-section study with D yields in events carrying the `HLT_noalg_L1AFP_A_OR_C` bit, without an offline tag (D⁰ 45 ± 30 / 9 ± 10 / 5 ± 4 and D⁺ 12 ± 14 / 0 ± 6 / 5 ± 6 per pT bin), and the probability of that bit for events with a D (0.41 ± 0.22 % in run 428770, about 0.5 ± 0.3 % in run 427929); these numbers are related to the R1 prediction, were not used, and nothing was adjusted to them (the R1 computation was already written; afterwards only the ratio of the AFP efficiencies, 0.94/0.95, was added to it, −1 %). (2) On 20.11 at about 02:20 CET, after the numbers had been frozen (02:01) and the report written, a file search displayed the name of a note of the tagged study with the result in its title (R^fid 0.99 %, variant 0.54 %, σ 13.6 µb, ρ 0.96 ± 0.12); the file was not opened and no prediction was changed afterwards (file time stamps 02:01).
<!-- src: charm_work/AFP_akceptancja_2026-11-20/PLAN.md; charm_work/AFP_akceptancja_2026-11-20/tabele (T1–T3 read; T4–T9 not read for this addendum); charm_work/TAGXSEC_2026-11-20/RAPORT.md; predictions and statement of independence: charm_work/SD_var_2026-11-19/RAPORT_domkniecie.md §0, §4, §5, out/pred_tables.md §1, §2; R1 yields: charm_work/TAGGED_2026-11-19/RAPORT.md (table of yields) -->

### 6 Uncertainty budget (D⁰, relative)

| source | % | enters |
|---|---|---|
| statistical | −34 / +36 | R and σ |
| tagged-fit shape | 23.2 | R and σ |
| ρ (multiplicity dependence of the efficiency) | 12.5 | R and σ |
| inclusive yield: statistics (larger of bootstrap and 1.3 × HESSE), background model | 8.7, 5.5 | R |
| inclusive yield: second minimum of the fit, peaking template | 3.5, 0.4 | R |
| AFP acceptance, efficiency, alignment | 5.1 | R and σ |
| inclusive fit bias | 18.6 | σ only |
| classifier scale factor | 4.4 | σ only |
| **systematic total** | **29.1 (R), 33.2 (σ)** | |
| model of the wide components (nominal against variant) | −46, one-sided | quoted separately |
Before the update with the final inclusive inputs the rows were: inclusive statistics 8.4, fit bias 19.4, scale factor 7.3, totals 28.8 (R) and 34.1 (σ); the rows for the second minimum and the peaking template are new. With the inclusive inputs of version v2 the scale-factor row is 4.4 instead of 4.2 and the totals stay 29.1 (R) and 33.2 (σ). Variant of the shape row (post hoc, Sec. 3): 11.0 instead of 12.6 events, i.e. 20.3 instead of 23.2 %, which would give totals of 26.8 (R) and 31.3 (σ); not adopted. The statistical row is −34.3 / +35.6 % of N_tag − N_acc (the earlier "+35" was a rounding slip).
<!-- src: charm_work/TAGXSEC_2026-11-20/RAPORT.md (budget), out/tagxsec_excl.json; charm_work/TAGXSEC_2026-11-20/akapit_do_noty.md; update: charm_work/XSEC_final_2026-11-20/out/final.json (rel per bin), charm_work/AUDYT_noty_2026-11-20/out/przeliczenia.md §3; v2: charm_work/ZAMKNIECIE_2026-11-20/out/przeliczenia.md §3; shape variant: charm_work/KONSTR_2026-11-20/RAPORT.md §B -->

### 7 Caveats

- One run, about 55 tagged D⁰; significance about 3σ. Nothing here is differential.
- **Updated to the final inclusive computation (audit of 20.11).** σ_tag^fid, R^fid and the budget of Secs. 4 and 6 use the final inclusive values of the note (Sec. 6.13: 1028 + 278 + 100 µb; scale factor 0.97 ± 0.03 at 3–5 GeV, Sec. 6.8); the update is arithmetic on the stored inputs (N_tag, g_b, C_b, ρ_b unchanged), not a new fit. Not updated: the weights of the fit-bias row (inclusive yields of 19.11) and ρ. The ± 0.08 for the D⁺ scale factor at 3–5 GeV, adopted in the note on 20.11, changes nothing here (D⁰ only); the larger error of the D⁰ scale factor above 5 GeV changes the scale-factor row only (4.2 → 4.4 %).
- The cross-check with a constructed background shape (Sec. 3) is post hoc. It supports the nominal yield and makes the window-minus-sidebands variant an estimator-specific low edge, but it does not change any quoted number.
- The Pythia envelope of Sec. 4 is carried over from 2–30 GeV. The completion of the generator study on 20.11 (missing parts of the colour-reconnection variation) leaves it unchanged (upper factor 3.97 → 3.98, less than 0.3 %), so the values of Sec. 4 were not recomputed. The measurement-against-prediction table of Sec. 5 is post hoc as a comparison; the choice of WP70 as its nominal row was made after the reading.
- An extrapolation to the total single-diffractive cross-section (factor about 3.9) is not recommended and not given.
- The D⁺ result failed its control criteria; no D⁺ cross-section is quoted.
- Run 428770 shows tagged peaks with higher significance but has a trigger bias, an accidental-tag fraction of 5–6 % and no absolute normalisation; it is not used for a cross-section.
- Comparison with theory predictions from the literature is deliberately absent: the source numbers have not been checked against the papers.
<!-- src: charm_work/TAGXSEC_2026-11-20/RAPORT.md, akapit_do_noty.md (Uwagi); charm_work/XSEC_final_2026-11-20/tabele/final.md F1; charm_work/SFABS_wynik_2026-11-20/RAPORT.md §0; charm_work/REG428_2026-11-20/RAPORT.md §5 -->

### 8 Figures (local, not in the note)

- T1: `~/Documents/charm/charm_work/PREZENTACJA_2026-11-20/local_tagged/T1_tagged_D0_mass.png` — D⁰ mass spectrum in events with an AFP proton tag, nominal fit and the fit with the constructed background shape (Sec. 3).
- T2: `~/Documents/charm/charm_work/PREZENTACJA_2026-11-20/local_tagged/T2_tagged_R_sigma.png` — R^fid and σ_tag^fid against the Pythia envelopes. The figure is drawn from the first computation (σ_tag^fid = 13.6 µb, systematic of R^fid 0.28 %), not from the updated values of Sec. 4 (13.9 µb, 0.29 %).
<!-- src: charm_work/PREZENTACJA_2026-11-20/RAPORT.md; skrypty/t1_t2_tagged.py (input TAGXSEC_2026-11-20/out/tagxsec_excl.json), skrypty/build_review.py -->

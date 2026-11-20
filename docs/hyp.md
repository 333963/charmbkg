# HYP — hipotezy D⁰ K↔π i arbiter PID („obiekty sygnałopodobne”)

Stan: 13.11.2026. Kod: `charmbkg.hyp` (arbiter, bloki zdarzeń, prawda z `hypo_swap`), `charmbkg.hyp_extract` (ekstrakcja do parquet,
CLI `charmbkg hyp-extract`), `charmbkg.hyp_toy` (toy dE/dx), `charmbkg.dstar_tag` (tag D*⁺→D⁰π⁺_soft na danych).
Notatki robocze (poza repo): `charm_work/HYP_d0/` (`report/hyp_report.md`, `t4/t4.md`, `t4/t4_results.md` §5–8).

## Idea
Ta sama para torów D⁰ daje dwa kandydaty: K⁻π⁺ i π⁻K⁺ (zamienione przypisanie mas). Fałszywa hipoteza jest „tłem przez konstrukcję”:
ma te same tory, wierzchołek i kinematykę co sygnał, różni się tylko przypisaniem PID. Arbiter wybiera hipotezę **wyłącznie z PID**
(dE/dx, npts, log p, pulle Bethe-Blocha, Δχ²) — bez α, cosθ*, pT, które odtwarzają masę i rzeźbią tło.
Klasyfikator antysymetryczny (P_A + P_B = 1), k-fold OOF po zdarzeniu.

## Wyniki
| co | wynik | źródło |
|---|---|---|
| toy: ε wyboru poprawnej hipotezy | arbiter 54,8 % (92 % dla p_min 0,7–1 GeV), hash 49,7 %, closest_mass 83 % **ale rzeźbi tło (Z = +168)**, pid_kin rzeźbi (Z = +72) | HYP toy |
| arbiter v1 (cechy par) | rzeźbi tło (T1) | t4 §6 |
| arbiter v2 (`mode: pertrack`) na danych z tagiem D* | **SF = 0,977 ± 0,041** (dane = MC w 0,6σ); v2 ≠ losowy tylko 1,9σ; nie rzeźbi (|Z| < 3 w 24 testach) | t4 §5–6 (po `fit_passfail`, 13.11) |
| po poprawce r(m) (10.11) | v2 ≈ hash | t4 §7–8 |
| pojedyncze D⁰ | keep_both najlepsze we wszystkich binach (σ_rel o 7–15 % niższe niż v2) | t4 §8 |
| pary | hash → v2 bez mierzalnego zysku σ_μ (najwyżej przy niskim pT nogi) | t4 §8 |
| R/S refleksji @WP70 | dane/MC 0,96 ± 0,04 (nowy `fit_passfail`, 13.11) → syst. szablonu refleksji można zawęzić z ±10 % do ~±5 % | t4 §4 |

## Decyzje
- Pojedyncze D⁰: `signal.swap_resolution: keep_both` (bez arbitra).
- Pary (projekt A, charmtag): `pairs.d0_arbitration: hash`. Przełączenie na v2 bezpieczne, ale opłaca się tylko przy okazji kolejnego pełnego biegu par.

## Otwarte
- Czy fit par (fit2d) z wyborem nogi hashem ma ten sam efekt r(m) w tle kombinatorycznym — niesprawdzone; closure par na MC.
- Kalibracja PID per tor na torach z tagu D* (research K2).
- Zawężenie systematyki szablonu refleksji do ~±5 % (przenieść do A jako zmianę systematyki).

# PREZENTACJA — rysunki końcowych wyników projektu B (20.11.2026)

Wszystko lokalnie; skrypty tylko czytają wejścia (uruchamiane z `python -B`).
Układ rysunku F4 wg wzorca `~/Downloads/plots_doubleD/*_margins.png` (mapa 2D + rozkłady brzegowe + linie mas nominalnych + panel reszt „corr. excess”); plików wzorca nie zmieniano.
Styl wspólny: `skrypty/pstyle.py` (jedna paleta kategoryczna sprawdzona walidatorem CVD, dane zawsze czarne, √s, run, L).
Każda liczba na rysunku jest wczytywana ze źródła; tam, gdzie fit jest odtwarzany, asercja porównuje go z zapisanym wynikiem.

Odtworzenie: `cd skrypty && ~/Documents/charm/.venv_h5/bin/python -B <skrypt>.py` (f1_fits, f2_f3_xsec, f4_map, f5_f9, t1_t2_tagged, build_review).

## Rysunki → źródła

| rysunek | plik (PNG + PDF) | skrypt | źródło liczb |
|---|---|---|---|
| F1a D⁰, fit M8, 3 biny | `fig/F1a_mass_fit_D0_M8` | `f1_fits.py` | histogramy `WPSCAN_2026-11-19/runs/lxplus/wps_hist.npz`; kod fitu `FITBIAS_2026-11-19/skrypty/{fb_fit,fb_var}.py` + `BKGCAT_2026-11-20/skrypty/bc_m8.py` (szablon „dir”, r = 1) importowany bez zmian; N, błąd stat., p(χ²) z `XSEC_final_2026-11-20/v2/out/final.json` (asercja: fit odtworzony = final.json) |
| F1b D⁺, fit A2, 3 biny | `fig/F1b_mass_fit_Dplus_A2` | `f1_fits.py` | jw. + `FITBIAS_Dplus_2026-11-20/skrypty/fd_lib.py` (A2), szablony `BKGCAT_2026-11-20/in/bkgcat.npz`; N, stat, p z `final.json` |
| F2 przekroje D⁰, D⁺, D*± | `fig/F2_cross_sections` | `f2_f3_xsec.py` | wyłącznie `XSEC_final_2026-11-20/v2/out/final.json` (`cells`, `cmp`, `dstar`); ALICE oznaczone „post hoc, unverified” |
| F3 stosunki D⁺/D⁰, D*⁺/D⁰ | `fig/F3_ratios` | `f2_f3_xsec.py` | `final.json` (`ratio`, `dstar`) |
| F4 mapa m(Kπ) × Δm | `fig/F4_map_mKpi_dm` | `f4_map.py` | lokalne tabele par `INCAL_2026-11-19/runs/lxplus/parts/inc_d150_data_427929_*.npz` (`D0__pair_*`), wczytanie `XSEC_atlas_2026-11-20/skrypty/xa_lib.py`; selekcja jak `xa_dstar.sel` bez okna masy D⁰, pT(D*) 5–30 GeV; 3138 par; bez fitu |
| F5 izolacja SF | `fig/F5_isolation_sf` | `f5_f9.py f5` | `B_trig/runs/lxplus/out/btrig_iso.json` (ptcone40 < 2 GeV; 0,98 ± 0,04 vs 0,51 ± 0,02 = liczby noty §6.2.2) |
| F6 tryger | `fig/F6_trigger` | `f5_f9.py f6` | `REG428_2026-11-20/out/summary.json` (`eps`, `yields`; = tabele T1, T3) |
| F7 budżet niepewności | `fig/F7_uncertainty_budget` | `f5_f9.py f7` | `final.json` (`cells[*].rel`, `L_rel`; asercja: suma w kwadracie = `tot_rel`) |
| F8 skład tła (MC) | `fig/F8_background_composition` | `f5_f9.py f8` | `BKGCAT_2026-11-20/out/comp.json` (WP70, okno; asercja: suma kategorii = tło) |
| F9 sieć | `fig/F9_network` | `f5_f9.py f9` | `DNET/runs/lxplus6/dnet6_p1.json` (DR), `dnet6_results.json` (G proxy), `dnet6_R.json` (G_fit z fitu piku, ziarno 1, i R runu 428770 — dodane 20.11; asercje: próg 1,15 w ≥ 2 z 3 binów niespełniony w 427929, spełniony w 428770, bramka domknięcia MC niezaliczona); zgodne z `WERDYKT_etap6.md` §4, §5 |
| T1 widmo D⁰ z tagiem | `local_tagged/T1_tagged_D0_mass` | `t1_t2_tagged.py` | nominał: `TAGGED_2026-11-19/out/{hists.npz,unblind.json}` (bez nowego fitu); wariant: `KONSTR_2026-11-20/skrypty/kb_tag.py` importowany bez zmian, fit odtworzony (asercja z `out/kb_window.json`) |
| T2 R^fid, σ_tag^fid | `local_tagged/T2_tagged_R_sigma` | `t1_t2_tagged.py` | `TAGXSEC_2026-11-20/out/tagxsec_excl.json` (`wynik`, `pred`); zakresy literatury parsowane z tabeli `TAGXSEC_2026-11-20/RAPORT.md` §5 (nie ma ich w JSON) |

Przegląd: `PRZEGLAD.html` (zwykły HTML, ścieżki względne) i `PRZEGLAD.md` — ta sama treść, generowane przez `skrypty/build_review.py`; sekcja tagowana na końcu pod nagłówkiem „AFP-tagged measurement”.

## Czego nie narysowano / ograniczenia

- **F1:** fity rysowane w dwóch plikach (D⁰, D⁺) zamiast jednego — 3 biny × 3 panele (widmo, dane − tło kombinatoryczne ze składowymi, pulle) nie mieszczą się czytelnie na jednym rysunku. Błąd N na rysunku to max(bootstrap, 1,3 × HESSE) z `final.json`, nie HESSE fitu.
- **F2:** pasmo FONLL jest tylko dla części prompt (tak w `final.json`), więc w górnym panelu leży poniżej danych z definicji; porównanie ilościowe jest w dolnym panelu (część prompt danych / FONLL). CMS tylko w dolnym panelu (inna definicja; dane przeliczone). D*± 3–5 GeV: brak pomiaru. Dla D*± w `final.json` nie ma FONLL ani ALICE.
- **F4:** zakres tabel par to m(Kπ) 1790–1940 MeV, Δm < 156 MeV — mapa nie sięga dalej; pT(D*) 5–30 GeV (wybór prezentacyjny: biny, w których mierzono σ(D*)). Ekstrakcji na lxplus nie uruchamiano.
- **F5:** „skutek” pokazany jako ε(dane) vs ε(MC) i czynnik 0,51 na rysunku; osobnego rachunku wpływu na przekrój nie ma w źródłach (dla 428770 przekroju bezwzględnego się nie podaje).
- **F6:** udziały regionów liczone względem fitu wszystkich zdarzeń; suma rozłącznych regionów to 96–102 % (podane nad słupkami). R1 < 1 % — niewidoczny.
- **F9 (aktualizacja 20.11, ostatnie dopisanie; kopie `*.bak_przed_2011ost`):** obok proxy narysowany zysk z fitu piku G_fit (`dnet6_R.json`; jedno ziarno, bez błędu): 427929 0,99 / 1,10 / 1,09, 428770 1,13 / 1,25 / 1,17, z linią progu 1,15; R jest opisowe (bramka domknięcia na MC niezaliczona 3/3). Opis w `PRZEGLAD.{html,md}` zmieniony (`build_review.py`). Poprzedni stan: rysunek pokazywał wyłącznie proxy.
- **Histogramy WPSCAN w przedziałach score:** użyte tylko jako wejście fitów przy WP70 (F1); osobnego rysunku skanu WP nie robiono (nie było na liście F1–F9; skan WP jest w nocie §5.12 liczbowo).
- **T2:** przewidywanie bezwzględne Pythii = R × σ_incl Pythii (`py_incl_ub`); dla scenariusza twardego i Siddikova–Schmidta wartości bezwzględnej w źródłach nie ma („not available”).
- Rysunków noty `docs/figures/` (F1, F4, F6, F8–F10, F13–F17) nie dublowano: dotyczą metod (zamknięcia, SF vs WP, klasy trygera szczegółowo); F5 jest jednopanelowym skrótem F14.

## Uwaga do wzorca (projekt A, bez analizy)

W `h_matched_d0_dplus_margins.png` w masie kandydata D⁺ widać drugi pik przy ~2010 MeV — odpowiada masie D*⁺ (D*⁺ → D⁰π⁺ → Kππ zrekonstruowane jako kandydat D⁺); do sprawdzenia w projekcie A.

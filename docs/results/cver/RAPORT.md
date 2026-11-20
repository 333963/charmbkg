# CVER — weryfikacja definicji C = A·ε w XSEC względem trybu rozpadu (projekt B, 2026-11-20)

Niezależne sprawdzenie sporu XSEC ↔ BKGCAT. Tylko MC mb_nd (78 plików ROOT v2.0, 19 220 986 zdarzeń) + odczyt plików XSEC/BKGCAT/WPSCAN.
Bez AFP, bez condora; XSEC / FITBIAS / BKGCAT / nota bez zmian.
Skrypty: `skrypty/cver_truth.py` (lxplus, interaktywnie, 6 procesów × ~2,5 min, `/tmp/klyczak_cver`), `skrypty/cver_ana.py` (lokalnie). Liczby: `out/cver.json`, `out/tabele.md`.

## 0. Werdykt

**Rację ma BKGCAT. Przekroje XSEC są zawyżone: σ_popr = σ_XSEC × 0,734 ± 0,002 (D⁰) i × 0,939 ± 0,001 (D⁺); D⁺/D⁰ × 1,279 ± 0,004.** Czynnik nie zależy od pT.
Mianownik C liczy prawdziwe D po `TruthDmesons_decayType` (10 / 4), a ten znacznik jest topologiczny (córki naładowane + π⁰ = dokładnie K∓π±[π±], dowolne inne neutralne);
licznik C i uzysk z fitu to wyłącznie ekskluzywne K⁻π⁺(γ) / K⁻π⁺π⁺(γ); mnoży się to przez BR_PDG trybu ekskluzywnego. Licznik i mianownik nie są w tym samym trybie.

## 1. Wzór faktycznie zaimplementowany w XSEC (`xsec.py` w. 181–237)

σ = N_fit / (C × SF × BR_PDG × L), C = `N_truth` / `N_gen_sig`

- **Mianownik** `N_gen_sig`: `GEN_mcstat/runs/lxplus_v20/truthD_mb.parquet`, kolumna `mode == 1`, |pdg| = 421 / 411, pT prawdy w binie, |η| < 2,2.
  `mode` powstaje w `Q_lxplus/gen_norm/effmap_gn.py` w. 62: `TruthDmesons_decayType == {421: 10, 411: 4}`. Żadnego warunku na pdgId córek ani na masę widzialną.
- **Licznik** `N_truth`: `WPSCAN_2026-11-19/wps_fit.json` (`wps_fit.py` w. 98) = suma histogramu `mb|sig` z `wps_hist.py` w. 127–128: kandydaci mb_nd z `ok`, score ≥ WP70,
  `Dmeson_truth_category` ∈ {0, 1}, D⁰: `swap ≠ 1`, masa w zakresie histogramu (D⁰ 1,70–2,05, D⁺ 1,70–2,20 GeV), bin pT reco, bez cięcia η prawdy.
- **BR**: D⁰ 3,936 + 0,0148 (DCS) = 3,951 %, D⁺ 9,31 % (PDG, tryb ekskluzywny z FSR).
- „BR w MC” 5,39 % / 9,95 % w XSEC = N(`mode == 1`)/N(wszystkie D) — czyli ułamek znacznika topologicznego, nie BR trybu.
- Pythia A3 w XSEC: `N_gen_all` (wszystkie D⁰ / D⁺ w obszarze, dowolny rozpad) × σ_ND / N_zdarzeń — **spójne**, niezależne od `decayType`. To samo dla stosunku D⁺/D⁰ Pythii.

## 2. Prawda generatora (własny skrypt; klasyfikacja po liście córek i bilansie czteropędu, bez użycia `decayType`)

Lista `TruthDmesons_daughterPdgId` zawiera wyłącznie K±, π±, π⁰ (plus w 66 % wierszy wpis samego mezonu — usuwany); fotonów, K⁰, η, leptonów nie ma.
`decayType` 10 ⇔ lista = {K⁻, π⁺} (33 950 / 33 950), `decayType` 4 ⇔ {K⁻, π⁺, π⁺} (32 824 / 32 824). „Ekskluzywny” = m_D − m(Kπ[π]) < 1 MeV (bez fotonu) albo brakujący czteropęd o |m_brak| < 100 MeV (FSR).

| | bin pT | N_gen wszystkie | N_gen `decayType` (= XSEC, zgodne co do sztuki) | ekskl. bez γ | ekskl. z γ | + K⁰ | + η | + inne (η′, φ, 2 neutralne) | f_ekskl [%] | BR_MC z FSR [%] | bez FSR [%] | PDG [%] | MC/PDG |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D⁰ | 3–5 | 413 204 | 22 306 | 15 223 | 1143 | 2730 | 1644 | 1566 | 73,37 ± 0,30 | 3,961 ± 0,031 | 3,684 | 3,951 | 1,003 ± 0,008 |
| D⁰ | 5–8 | 162 404 | 8706 | 5942 | 452 | 1013 | 648 | 651 | 73,44 ± 0,47 | 3,937 ± 0,049 | 3,659 | 3,951 | 0,997 ± 0,012 |
| D⁰ | 8–30 | 54 168 | 2909 | 1965 | 178 | 319 | 211 | 236 | 73,67 ± 0,82 | 3,956 ± 0,085 | 3,628 | 3,951 | 1,001 ± 0,022 |
| D⁰ | ≥ 3 | 630 213 | 33 950 | 23 150 | 1775 | 4063 | 2507 | 2455 | **73,42 ± 0,24** | 3,955 ± 0,025 | 3,673 | 3,951 | 1,001 ± 0,006 |
| D⁺ | 3–5 | 215 142 | 21 458 | 18 962 | 1179 | 1234 | 1 | 82 | 93,86 ± 0,16 | 9,362 ± 0,066 | 8,814 | 9,31 | 1,006 ± 0,007 |
| D⁺ | 5–8 | 85 738 | 8440 | 7486 | 437 | 490 | 1 | 26 | 93,87 ± 0,26 | 9,241 ± 0,104 | 8,731 | 9,31 | 0,993 ± 0,011 |
| D⁺ | 8–30 | 28 917 | 2913 | 2562 | 172 | 170 | 0 | 9 | 93,86 ± 0,44 | 9,455 ± 0,181 | 8,860 | 9,31 | 1,016 ± 0,019 |
| D⁺ | ≥ 3 | 330 015 | 32 824 | 29 021 | 1790 | 1894 | 2 | 117 | **93,87 ± 0,13** | 9,336 ± 0,053 | 8,794 | 9,31 | 1,003 ± 0,006 |

- BR trybu sygnału w mb_nd zgadza się z PDG w granicach 0,6 % (stat.) — teza XSEC „BR w MC ×1,364 / ×1,068” jest nieprawdziwa; 1,364 = 1/0,734 to udział rozpadów z niewidocznymi neutralnymi w znaczniku.
- FSR (foton > 1 MeV w układzie D): 7,1 % (D⁰), 5,8 % (D⁺) rozpadów ekskluzywnych; PDG zawiera FSR, więc właściwa jest kolumna „z FSR”.
- PDG: `BR_frag_2026-11-18/RAPORT.md` §1 i notatka bazy wiedzy „BR i frakcje fragmentacji powabu…” (PDG 2026: 3,936 ± 0,030 %, DCS 0,0148 %, 9,31 ± 0,14 %).
- Drobiazg: ekskluzywne DCS K⁺π⁻ (87 wierszy, 0,0138 %) mają `decayType` −1, nie 10 — komentarz w `xsec.py` w. 42 („bez względu na ładunek”) jest błędny; BR z DCS (3,951) przy mianowniku bez DCS to niespójność 0,4 % (uzysk w danych przy `keep_both` DCS zawiera). Pominięte w czynniku.

## 3. Co jest w liczniku i w uzysku z fitu

- Licznik WP70 (tabela BKGCAT `in/bkgcat.npz`, sumy = `N_truth` XSEC co do sztuki: 3653 / 1919 / 848 i 3245 / 1827 / 843): **100 % trybu ekskluzywnego** (z FSR), 0 kandydatów z brakującą masą.
- Niezależne łącze (`TruthDmesons_recoMatchIndices`, poziom preselekcji, bez score): kandydaci przypisani do nieekskluzywnych wierszy `decayType` 10: 83 na 19 845 (0,42 %), płasko w 1,70–2,04 GeV (nie pik; K i π z K⁻π⁺K⁰ mają m(Kπ) ≤ 1,37 GeV — to przypisania z obcym torem); D⁺: 2 na 8462. Przy WP70 w łączu BKGCAT: 0. Górna granica wpływu na czynnik: 0,4 % (D⁰).
- Rozpady nieekskluzywne ze znacznika (K⁻π⁺K⁰, K⁻π⁺η, K⁻π⁺η′…) leżą poniżej 1,37 GeV, poza zakresem fitu; nie są ani sygnałem, ani składową M8. Składowe pikujące M8 (DMM: K⁻π⁺π⁰, π⁺π⁻π⁰, K⁻ℓν, KK, ππ; CM) mają `truth_category` 6/7, 2/3 — poza licznikiem, w fitcie osobny szablon o stałym κ.
- Szablon sygnału w fitcie = fd_nd `truth` ∈ {0, 1} → N_fit to wyłącznie prawdziwe K⁻π⁺(γ) / K⁻π⁺π⁺(γ). **Licznik C i N_fit są zgodne; niezgodny jest tylko mianownik.**

## 4. C poprawne i równoważność dwóch sformułowań

Def. 1: σ = N / (C_popr × BR_PDG × L), C_popr = N_reco / N_gen(ekskl.). Def. 2: σ = N / (C′ × L), C′ = N_reco / N_gen(wszystkie D).

| | bin | C_XSEC [%] | C_popr [%] | C_popr/C_XSEC | C′ [10⁻³] | C_popr × BR_PDG [10⁻³] | C′/(C_popr·BR) | σ_popr/σ_XSEC (def. 1) | (def. 2) |
|---|---|---|---|---|---|---|---|---|---|
| D⁰ | 3–5 | 16,38 | 22,32 | 1,363 ± 0,006 | 8,841 | 8,818 | 1,003 | 0,7337 ± 0,0030 | 0,732 ± 0,005 |
| D⁰ | 5–8 | 22,04 | 30,01 | 1,362 ± 0,009 | 11,816 | 11,857 | 0,997 | 0,7344 ± 0,0047 | 0,737 ± 0,008 |
| D⁰ | 8–30 | 29,15 | 39,57 | 1,357 ± 0,015 | 15,655 | 15,634 | 1,001 | 0,7367 ± 0,0082 | 0,736 ± 0,013 |
| D⁺ | 3–5 | 15,12 | 16,11 | 1,065 ± 0,002 | 15,083 | 15,000 | 1,006 | 0,9386 ± 0,0016 | 0,933 ± 0,006 |
| D⁺ | 5–8 | 21,65 | 23,06 | 1,065 ± 0,003 | 21,309 | 21,468 | 0,993 | 0,9387 ± 0,0026 | 0,946 ± 0,010 |
| D⁺ | 8–30 | 28,94 | 30,83 | 1,066 ± 0,005 | 29,152 | 28,706 | 1,016 | 0,9386 ± 0,0044 | 0,924 ± 0,016 |

Def. 1 i 2 różnią się tylko o BR_MC/BR_PDG (≤ 1,6 %, w granicach statystyki MC) — dają to samo. XSEC zrobił hybrydę: mianownik „po znaczniku topologicznym” (ani tryb, ani wszystkie D) × BR_PDG trybu ekskluzywnego.
Zalecany czynnik (wspólny dla binów, def. 1, błąd dwumianowy MC): **D⁰ 0,7342 ± 0,0024; D⁺ 0,9387 ± 0,0013**. Błąd stat. MC samego C bez zmian (1,5–2,9 %).

## 5. Wielkości poprawione (427929; wartości XSEC × czynnik, błędy przeskalowane proporcjonalnie; nic nie przeliczano w XSEC)

| | bin | σ XSEC [µb] | σ popr. [µb] | dane/Pythia A3: było → jest |
|---|---|---|---|---|
| D⁰ | 3–5 (SF = 1) | 1365 | 1002 (SFABS: 1411 → 1036 ± 117 ± 277 ± 23) | 1,11 → 0,81 (SFABS 1,15 → 0,84) |
| D⁰ | 5–8 | 376 ± 31 ± 57 ± 8 | **276 ± 23 ± 42 ± 6** | 0,78 → 0,57 |
| D⁰ | 8–30 | 132 ± 11 ± 24 ± 3 | **97 ± 8 ± 18 ± 2** | 0,82 → 0,61 |
| D⁺ | 3–5 (SF = 1) | 542 | 509 (SFABS: 571 → 536 ± 42 ± 56 ± 12) | 0,85 → 0,79 |
| D⁺ | 5–8 | 176 | 165,5 (SFABS 178 → 167) | 0,69 → 0,65 |
| D⁺ | 8–30 | 54,8 | 51,5 | 0,64 → 0,60 |

D⁺/D⁰ (SF⁺ = SF⁰): 0,397 / 0,466 / 0,431 → **0,508 / 0,596 / 0,549** (× 1,279 ± 0,006 / 0,009 / 0,015); Pythia A3: 0,521 / 0,528 / 0,534. SFABS 3–5 GeV: 0,405 → 0,518.
Po poprawce D⁺/D⁰ zgadza się z Pythią, a dane/Pythia są jednakowe dla obu gatunków (0,57–0,65 przy ≥ 5 GeV) — wcześniejszy „niski” stosunek D⁺/D⁰ był tym artefaktem.

Dotknięte (wszystko, co zawiera C_XSEC × BR_PDG):
1. `XSEC out/tabele.md` T1, **T1b**, T2 (D⁺/D⁰), T5 (dane na zdarzenie / Pythia: D⁰ 1,24 / 0,87 / 0,91 → 0,91 / 0,64 / 0,67; D⁺ 0,94 / 0,77 / 0,71 → 0,88 / 0,72 / 0,67), T6 (σ_AFP regionu R1 — ta sama C; stosunki σ_AFP/σ_inkl tego samego gatunku bez zmian), `WYNIKI.md` §1, §2, §5, §7, §8.3–8.5, `DEFINICJA.md` §2–3 (akapit „BR w tym MC ≠ PDG”).
2. Porównania: CMS (D⁰ 5–8: 0,72 → **0,53 ± 0,11**, 8–30: 0,97 → 0,71 ± 0,16; D⁺ 0,93 / 0,90 → 0,87 / 0,84); FONLL (D⁰ 1,52 / 1,46 → 1,12 / 1,07; D⁺ 1,79 / 1,52 → 1,68 / 1,43); Pythia (wyżej). **Uwaga: D⁰ 5–8 GeV wypada teraz ~3σ poniżej CMS** — poprawka jest pewna, więc napięcie wskazuje na coś innego (SF, preselekcja dane/MC, obciążenie fitu), nie na błąd tej poprawki.
3. Scenariusz D⁺ (T1b) i `SFABS_wynik_2026-11-20` (RAPORT §„Przekrój”, `akapit_do_noty.md`: 571 → 536, 1411 → 1036, 0,41 → 0,52).
4. Przewidywania uzysków (także tagowanych): liczone jako σ_XSEC × C_XSEC × BR × L — bez zmian (błąd się skraca); liczone z Pythii/FONLL (wszystkie D) × C_XSEC × BR_PDG — zaniżone ×0,734 (D⁰), ×0,939 (D⁺); liczone z „sygnału MC na zdarzenie” przeskalowanego o 1/1,364 — błędnie zaniżone (skalowanie niepotrzebne).
5. **κ × 1,364 (D⁰) / × 1,068 (D⁺) w FITBIAS/M8** (`fb_var.py --kscale`): artefakt — BR sygnału w MC = PDG (pkt 2), więc κ_dane = κ_MC, **r = 1**. Potwierdzam rozumowanie BKGCAT; wpływu na uzyski sam nie przefitowywałem — wg `BKGCAT tabele/m8.md` (r = 1 + szablon wprost) 427929: −1,1 / +0,4 / +2,8 %, 428770: −1,3 / −0,3 / +2,3 % (sama skala: do +3,5 % i +8,4 % przy 8–30). Uzyski D⁰ w `fitbias_override.json` mają κ × 1,364 → dodatkowa zmiana σ(D⁰) o te ≤ 3 % (osobno od czynnika 0,734). Zdanie FITBIAS „efekt DMM w danych −7 % zamiast −5 %” odpada.
6. Wyjaśnienie „B7: D*/zdarzenie mb_nd/A3 = 1,36” i „łańcuch D* 3,70 % wobec 2,665 %” jako BR generatora — nieprawdziwe; to ten sam znacznik topologiczny. Do sprawdzenia przez właściciela `GEN_mcstat`/`gen_norm` (projekt A; nie ruszałem): mapy wydajności z mianownikiem `mode == 1` są wewnętrznie spójne (licznik/mianownik × liczba prawdziwych D z tym samym znacznikiem), ale każde mnożenie przez BR_PDG albo porównanie z samodzielną A3 wymaga czynnika 0,734 / 0,939.

Nie dotyczy: stosunki binów pT tego samego gatunku, stosunek runów, SF z tagu D* (DNET), świetlność, frakcja nieprompt, σ_Pythia.

## 6. Pewność i zastrzeżenia

- Werdykt: **bardzo wysoka pewność.** N_gen z `decayType` odtwarza mianowniki XSEC co do sztuki (22 306 / 8706 / 2909; 21 458 / 8440 / 2913); klasyfikacja ekskluzywna zrobiona dwiema niezależnymi drogami (BKGCAT: 73,62 / 94,01 % przy |η| < 2,5; tu: 73,42 / 93,87 % w obszarze fiducjalnym), piki brakującej masy przy m(K⁰) i m(η) jednoznaczne; BR ekskluzywne = PDG w 0,6 %.
- Czynnik: błąd stat. 0,3 % (D⁰), 0,15 % (D⁺); systematycznie ≤ 0,4 % (kandydaci z nieekskluzywnych wierszy na poziomie preselekcji; DCS).
- Nie sprawdzano: przefitowania M8 z r = 1 (liczby BKGCAT), fd_nd (kształt sygnału — nie wpływa na C), 428770 (przekroju bezwzględnego XSEC nie podaje; σ × L skaluje się tak samo).
- Pierwsza wersja `cver_truth.py` nie usuwała wpisu samego mezonu z listy córek (wynik bezużyteczny, zachowany w `in/v0_blad_self/`); poprawiona wersja = `in/part_grp_*.npz`.

## 7. Pozostałości

- lxplus9103 `/tmp/klyczak_cver` (skrypt, listy plików, `part_grp_*.npz` v0 i `v1/`, logi; ~55 MB).
- Lokalnie `CVER_2026-11-20/{in,out,skrypty}`; `in/v0_blad_self/` do skasowania za zgodą.

# SFABS_wynik — bezwzględny SF (dane/MC) wydajności drzewa H przy WP70: D⁺ (3–5 / 5–8 / 8–30 GeV) i D⁰ 3–5 GeV (projekt B)

20.11.2026, ~01:15 CET. Kryteria: `KRYTERIA.md` (zapisane 19.11 ~23:55, przed wszystkimi fitami danych tego zadania). **Całość post hoc** względem etapów 4/5a, LOWPT, WPSCAN, WPTIGHT, FITBIAS, INCAL, dnet4.
Praca wyłącznie lokalna: bez jobów, bez lxplus, bez zmian w ORCH/XSEC/nocie. Bez AFP; tylko piki inkluzywne i tag D*.
Wejścia: `WPSCAN_2026-11-19/runs/lxplus/wps_hist.npz` (A), tabele par D⁰←D* i sygnału MC z `INCAL_2026-11-19/runs/lxplus/parts` (B, C; pełna statystyka obu runów — nowa ekstrakcja nie była potrzebna), moduły `FITBIAS/skrypty/fb_fit.py, fb_var.py` i `INCAL/scripts/inc_lib.py` (import, bez zmian).
Tabele: `tabele/{a_tab,b_tab,c_tab,d_tab}.md`; liczby: `out/final.json`, `out/a_summary.json`; akapit do noty: `akapit_do_noty.md`.

## 0. Wynik

| gatunek | bin [GeV] | A drabina (bezpośrednia) | B tag D* 2D (bezpośrednia) | C przeniesienie z D⁰ ≥ 5 GeV | zgodność | **SF(WP70) łączny** |
|---|---|---|---|---|---|---|
| D⁺ | 3–5 | 0,74 ± 0,17 (tylko 427929; σ > 0,15 → „ograniczenie”) | — | 0,968 ± 0,029 ± 0,012 ± 0,039 | 1,3σ | **0,95 ± 0,05** |
| D⁺ | 5–8 | 1,17 ± 0,15 | — | 0,968 ± 0,029 ± 0,012 ± 0,045 | 1,3σ | **0,99 ± 0,05** |
| D⁺ | 8–30 | 1,09 ± 0,12 | — | 0,968 ± 0,029 ± 0,012 ± 0,048 | 0,9σ | **0,99 ± 0,05** |
| D⁰ | 3–5 | nie przechodzi bramek (opisowo 1,00 ± 0,27) | 0,96 ± 0,16 ± 0,08 (427929; 428770 nie przechodzi) | 0,968 ± 0,029 ± 0,012 ± 0,012 | 0,03σ | **0,97 ± 0,03** |

(C: ± stat D* ± model fitu ± przeniesienie.) **Rekomendacja: zamiast scenariusza „SF = 1 ± 10 %” wstawić wartości z ostatniej kolumny.**
Trzeba przy tym wiedzieć, na czym stoi błąd ±0,05: wartość łączna ma wagę C 0,81–0,97, a C jest **warunkowe** (różnica dane–MC typu rozdzielczość / dE/dx). Metody bezpośrednie w binie docelowym mają ±0,12–0,18 i same nie dałyby nic lepszego niż scenariusz.
Jedyny bezpośredni test założenia C dla D⁺: drabina uśredniona po trzech binach = **1,03 ± 0,08** wobec SF_ref 0,968 ± 0,029 (+0,7σ). Wariant ostrożny (gdyby nie ufać przeniesieniu między gatunkami): D⁺ ± 0,08.

Przekrój 427929 (własna kopia rachunku; D + anty-D, prompt + nieprompt, |η| < 2,2; wejście `XSEC_2026-11-19/out/xsec_results.json` w stanie z 19.11 23:49, D⁰ z modelem M8):

| | σ przy SF = 1 [µb] | **σ [µb] ± stat ± syst ± L** | razem | było (scen.) | dane/Pythia A3 |
|---|---|---|---|---|---|
| D⁰ 3–5 | 1365 | **1411 ± 159 ± 378 ± 32** | 29,2 % | 30,7 % | 1,15 ± 0,33 |
| D⁺ 3–5 | 542 | **571 ± 45 ± 60 ± 13** | 13,2 % | 15,8 % | 0,89 ± 0,12 |
| D⁺ 5–8 | 176 | 178 ± 8 ± 25 ± 4 | 14,9 % | 17,2 % | 0,70 ± 0,10 |
| D⁺ 8–30 | 55 | 55 ± 5 ± 6 ± 1 | 14,4 % | 16,8 % | 0,65 ± 0,09 |

**D⁺/D⁰ (3–5 GeV) = 0,405 ± 0,056 (stat) ± 0,115 (syst)** (było 0,397 przy SF⁺ = SF⁰); SF⁺/SF⁰ = 0,981 ± 0,042 (wspólny SF_ref się skraca); błąd razem 31,6 % (było 32,9 %) — stosunek ogranicza fit D⁰, nie SF.
Zysk dla D⁺: błąd całkowity 15,8 → 13,2 % (3–5 GeV); dla D⁰ 3–5 GeV SF przestaje być istotnym składnikiem (3,4 % wobec 25 % obciążenia fitu).

## 1. Metoda A — drabina (`tabele/a_tab.md`)
SF70 = SF(L) × SF_rel(70 | L), L = 0,90; fit łączny kategorii [WP70, ·) i [WP90, WP70); D⁺ model M0, D⁰ M8. Szczebel dolny: stałe κ (stosunek nieefektywności dane/MC) + ograniczenie κ ∈ [0,5; 2] (SF(90) ∈ [0,90; 1,05]).
- **Poprawka na obciążenie fitu b_L (mb_nd, pełne tło):** D⁺ 1,063 ± 0,132 / 0,928 ± 0,056 / 0,910 ± 0,058 (błąd = max(rozrzut Poissona, 1,3 × HESSE)); D⁰ 0,63 ± 0,11 / 0,89 ± 0,09 / 1,10 ± 0,10 — D⁰ 3–5: fit kategorii luźnej zawyża sygnał ×2,6.
- **G_A1 (wstrzyknięcia na mb_nd):** D⁺ odzyskuje SF70 w ≤ 0,010 we wszystkich binach (np. 3–5: 0,911 / 0,902 i 0,880 / 0,870) — spełnione; D⁰ 3–5 nie (0,801 wobec 0,863); D⁰ 5–8, 8–30 spełnione.
  Zastrzeżenie: oba wstrzyknięcia spełniają założenie stałego κ, więc G_A1 testuje fit, **nie** szczebel dolny; ten jest ograniczony tylko przez σ_szcz (0,034–0,050).
- **G_A2 (tło z pasm danych + sygnał MC, 40 zabawek):** D⁺ 428770 3–5 nie przechodzi (obciążenie +4,0 / +4,8 % > 3 %) → w 3–5 GeV zostaje sam 427929. Pozostałe komórki D⁺: pull −0,19…+0,17, rms 0,83–1,05.
- **Dane (L = 0,90):** D⁺ 427929 3–5: r = 0,846 ± 0,100, po poprawce 0,796, κ̂ = 1,59, SF70 = 0,741 ± 0,112 (stat) ± 0,117 (b) ± 0,034 (szczebel) = ± 0,165; p(χ²) = 0,014 (na granicy).
  5–8: 1,173 ± 0,147 (oba runy, zgodne p = 0,89); 8–30: 1,085 ± 0,117 (p = 0,28). Kontrole L = 0,95 / 0,80: 0,78 ± 0,22 / 0,97 ± 0,23 (3–5); 1,16 / 1,06 (5–8); 1,11 / 1,06 (8–30) — spójne.
- **Dlaczego błąd jest duży:** założenie stałego κ wzmacnia błąd SF_rel ×1,2–1,4, a poprawka b jest znana tylko na 6–13 % (jedna realizacja mb_nd o statystyce ≈ 427929). Wg kryterium G_A4 komórka 3–5 GeV to „ograniczenie”, nie pomiar.
- D⁰ (opisowo): 3–5 nie przechodzi G_A1 (1,00 ± 0,27); 5–8: 1,19 ± 0,29 (427929); 8–30: 0,93 ± 0,17 (427929); 428770 odpada na G_A3 (κ̂ na brzegu albo p(χ²) = 0).

## 2. Metoda B — tag D*, fit 2D (m_Kπ, Δm) (`tabele/b_tab.md`)
Model: S (D* prawdziwe; szablony m i Δm z fd_nd, przesunięcia i skale wolne), R (prawdziwe D⁰ × przypadkowy pion: kształt m sygnału × próg), C (exp × próg (Δm − m_π)^p·e^(−cΔm)), P (pik w Δm, gładkie w m); ε parametrem; siatka 10 MeV × 0,2 MeV; pass + fail łącznie.
**Decyzja podjęta na MC przed danymi:** ε_R = ε (składowa R to prawdziwe D⁰ — na MC ε_R/ε = 1,05 ± 0,05 mb_nd, 1,01 fd_nd). Ze swobodnym R w kategorii „fail” błąd ε przy 3–5 GeV był 0,25 (jak w dnet4), z wiązaniem 0,12; wariant ε_R = (1 ± 0,05) ε zmienia SF o ≤ 0,01.
- **G_B1 (mb_nd):** ε_fit − ε_prawda = +0,038 ± 0,140 / +0,017 ± 0,083 / −0,014 ± 0,074 — spełnione (pulle +0,27 / +0,20 / −0,19).
- **G_B2 (czystość danych):** przechodzą 427929 3–5 (obc. −0,006 / +0,024), 427929 8–30, 428770 5–8, 428770 8–30; **nie przechodzą 427929 5–8 (obc. +0,049 ± 0,016; jedna zabawka z pull rms 6,4) i 428770 3–5 (+0,104 ± 0,022)**.
- **Dane:** D⁰ 3–5 GeV, 427929: N_D* = 1081 ± 165, ε = 0,674, **SF = 0,962 ± 0,162 (bootstrap; HESSE 0,141) ± 0,077 (model)**, p(χ²) = 0,58 — „pomiar” wg G_B4 (σ ≤ 0,25). Plan przewidywał ±0,15–0,20 — potwierdzone.
  428770 3–5: 0,741 ± 0,108 ± 0,190, p(χ²) = 0,002, warianty 0,55–0,89 — odrzucone (G_B2, G_B3); średniej runów nie ma.
  ≥ 5 GeV: 427929 8–30: 1,074 ± 0,106; 428770 5–8: 0,906 ± 0,111; 428770 8–30: **0,960 ± 0,033 ± 0,008** (dnet4: 0,968 ± 0,043).
- Uwaga: w 427929 fit 2D daje SF wyższy o ~0,10 niż dnet4 na tych samych danych (5–8: 1,10 / 1,01; 8–30: 1,07 / 0,97); w 5–8 zabawki pokazują obciążenie fitu 2D +0,05–0,07 w SF, w 8–30 +0,03 (w granicach bramki). Rząd wielkości różnicy metod: ≤ 0,1 przy małej statystyce, 0,01 w 428770 8–30.

## 3. Metoda C — przeniesienie (`tabele/c_tab.md`)
- **SF_ref = 0,968 ± 0,029 (stat) ± 0,012 (model)**: 4 komórki D⁰ ≥ 5 GeV, χ²/ndf = 1,45/3 (p = 0,69). Komórka 427929 5–8 wzięta z dnet4 (1,005 ± 0,126), bo fit 2D nie przeszedł G_B2 — zgodnie z p. C4 kryteriów; bez niej 0,965 ± 0,030; same liczby dnet4: 0,962 ± 0,036. Waga 428770 8–30: 0,80.
- **Czułość (fd_nd, sygnał z prawdą; modele rozmycia T i F z INCAL):** przy s = 1,10 SF(WP70) = 0,983–0,989 (D⁰), 0,960–0,990 (D⁺); dE/dx × 0,95: 0,989–0,994; rozmycie dE/dx 5 %: 0,997–0,999.
  δ_T = 0,040 / 0,047 / 0,049 (D⁺ 3–5 / 5–8 / 8–30), 0,012 (D⁰ 3–5). Dla D⁺ dominuje niejednoznaczność traktowania migracji przez próg L_xy (dla D⁺ kandydatów wchodzących spod progu 0,5 mm nie da się odtworzyć) razem z różnicą modeli T/F; samo dE/dx ≤ 0,006.
- Informacyjnie: gdyby skala rozdzielczości w danych wynosiła 1,05, SF celu zmieniłby się o +0,002…−0,016 — SF_ref pochłania większość tego efektu.

## 4. Połączenie (`tabele/d_tab.md`)
Wszystkie pary metod zgodne w ≤ 1,3σ → średnia ważona (metody niezależne). Żadna komórka nie wymagała reguły „niezgodne > 2σ”. Zgodność trzech binów drabiny D⁺ między sobą: χ² = 4,2/2 (p = 0,12).

## 5. Odstępstwa od kryteriów i zastrzeżenia
1. Fit 2D: wiązanie ε_R = ε i warianty ε_R = (1 ± 0,05) ε dodane po pierwszym przebiegu na MC, przed otwarciem danych (dozwolone p. B2); wariant „swobodne R” raportowany osobno i **nie** wliczony do systematyki modelu (dane 427929 3–5: 1,006 zamiast 0,962).
2. Test χ² fitu 2D przy 8–30 GeV w 427929 liczony na grubszej siatce (dodane po fitach danych — poprawka techniczna, bez wpływu na ε).
3. SF_ref: błąd modelu fitu (0,012) dodany do błędu C, choć kryteria mówiły o samej statystyce — w stronę ostrożną.
4. Bootstrap w B po grupach kandydatów (`uid >> 8`), nie ściśle po zdarzeniach; w A — repliki po zdarzeniach z WPSCAN (30). Liczba replik/zabawek 30–40: błędy rms znane na ±12 %.
5. C nie widzi efektów innych niż rozdzielczość / dE/dx (INCAL C5); zakres skali 1,00–1,10 pochodzi z diagnostyki INCAL (nie wyniku) i literatury. Bezpośredni test dla D⁺: ±0,08.
6. SF_ref jest zdominowany przez run dżetowy 428770 (8–30 GeV); zgodność z 427929 w tym binie 1,0σ (1,074 ± 0,106 wobec 0,960 ± 0,033). Stosowany do przekroju 427929.
7. Drabina: poprawka b z jednej realizacji mb_nd (skład minimum-bias); dla 428770 nie jest reprezentatywna. D⁺ 427929 3–5: p(χ²) fitu 0,014.
8. Przekrój: pozostałe składniki budżetu przepisane z XSEC bez zmian (obciążenie fitu D⁰ 25 %, model tła itd.); XSEC jest równolegle zmieniany — liczby σ przy SF = 1 odpowiadają plikowi z 19.11 23:49.
9. Jedno ziarno drzewa; bin 2–3 GeV poza zakresem. Liczba porównań jak w kryteriach (36 + 6 + 4 + 4), bez poprawki na wielokrotność.

## 6. Pliki i pozostałości
`KRYTERIA.md`, `RAPORT.md`, `akapit_do_noty.md`, `skrypty/{sfa_ladder.py, sfa_report.py, sfb_cache.py, sfb_fit2d.py, sfc_sens.py, sfd_final.py}`, `tabele/{a_tab,b_tab,c_tab,d_tab}.md`,
`out/{a_mc,a_data,a_toys,a_summary,b_mc,b_mc_tieR,b_data,b_toys,b_gof_fix,c_sens,final}.json` + logi (`b_mc.json` = `b_mc_tieR.json`: model nominalny).
Pozostałość poza projektem: skrót tabel par (4 pliki `pairs_*.npz`, ~40 MB) w katalogu tymczasowym — odtwarzalny przez `sfb_cache.py`. Na lxplus/EOS nic nie powstało.

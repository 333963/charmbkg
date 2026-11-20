# ISO_pis — wolny pion D* w zmiennych izolacyjnych D⁰ (19.11.2026, projekt B, tylko MC)

Pytanie: czy π_s z D*→D⁰π wchodzi do ptcone kandydata D⁰ i czy przez to ε_iso / SF izolacji z tagu D* (nota Sec. 6.2) różni się od inkluzywnej.

## Krótko

1. **Tak, π_s jest w stożku prawie zawsze.** Gałąź `Dmeson_ptcone40/100` nie wyklucza π_s (wyklucza tylko córki D⁰; gałąź − pT(π_s) nigdy < 0).
   π_s leży w ΔR < 0,4 u 95 / 99,8 / 99,8 % D⁰ z tagiem (3–5 / 5–8 / 8–30 GeV; brakujące 5 % w 3–5 to π_s z pT < 0,2 GeV, których nie ma w drzewie `Tracks`),
   ΔR mediana 0,12 / 0,08 / 0,04. Wkład do ptcone40: mediana 0,31 / 0,47 / 0,80 GeV (q10–q90: 0,23–0,41 / 0,34–0,63 / 0,56–1,28), czyli mediana 37 / 36 / 32 %
   ptcone40; u 35 / 29 / 19 % D⁰ π_s jest jedynym torem w stożku. Próg pT > 0,5 GeV (definicja `pc40` etapu 3) π_s przechodzi u 0,3 / 40 / 97 %.
2. **Efekt samego π_s (ε_c/ε_a) jest duży i stabilny między próbkami**: ptcone40 < 2 GeV: +5 / +12 / +31 % (fd_nd), +4 / +15 / +28 % (mb_nd);
   cięcia względne i ptcone100: +2…+11 %.
3. **Netto (ε tagu / ε inkluzywna) efekt jest mniejszy, bo D⁰ z tagiem bez π_s są w MC bardziej izolowane niż reszta** (ε_c/ε_b = 1,06 / 1,09 / 1,13 w fd_nd):
   ptcone40 < 2 GeV: ε_a/ε_i = 1,006 ± 0,004 / 0,984 ± 0,006 / **0,903 ± 0,014** (fd_nd), 0,99 ± 0,02 / 0,96 ± 0,04 / 0,87 ± 0,09 (mb_nd).
   Próg 3 % przekroczony tylko w 8–30 GeV dla ptcone40 < 2 GeV (ptcone40 < 4 GeV: −2,1 %) oraz dla stożków 100 (fd_nd: +1…+8 %, mb_nd: −4…−13 % — znak zależy od próbki).
   ptcone40/pT < 0,5: zgodne w 0,6 % (fd_nd).
4. **SF (dane/MC)**: przy SF ≈ 1 (427929) efekt drugiego rzędu jest ≤ 0,5 % (fd_nd) — zastrzeżenie nieistotne wobec błędu ±0,04. Przy SF 0,5–0,8 (428770)
   SF_inkl/SF_tag zależy od modelu przesunięcia aktywności: −2…+8 % przy SF 0,76, −10…+27 % przy SF 0,51; znaku nie da się ustalić z MC.

## Wejście i definicje

- Tabele DNET etapu 1 `/eos/user/k/klyczak/charm_q/dnet/data/dn_{fd_nd,mb_nd}.npz` (h5 v2.0, L_xy ≥ 0,25 mm, pT 3–30 GeV; wiersze `a_` = D⁰ inkluzywne,
  `p_` = pary D⁰←D* z Δm < 0,165 GeV, `p_pis_pt` = pT trzeciej córki D*) + tory otoczenia etapu 2 `/eos/.../charm_q/dnet2/parts/dn2_*.npz`
  (20 torów o największym pT w ΔR < 1,0 bez córek D⁰; `Tracks`: pT > 0,2 GeV). Tylko odczyt.
- ptcone = gałąź ntupla (`ptcone40_iso`, `ptcone100_iso` w `X`): Σ pT torów w ΔR < 0,4 / 1,0 bez córek D⁰, tory od ≈ 0,1 GeV (B_trig/TRKJET/RAPORT.md),
  π_s wliczony. Własna suma z tokenów zgadza się z gałęzią w 1 MeV u 89 % sygnału, |Δ| < 0,2 GeV u 99 %.
- Sygnał: truth ∈ {0, 1}, swap = 0, |m − m₀| < 55 MeV. Bez cięcia score_H (nota 6.2: score_H > 0,6) i bez selekcji AFP.
- (a) próbka tagu: wiersze `p_` z prawdą D* i D⁰, jeden na D⁰. (b) reszta: sygnał `a_` bez pary z prawdą D* (zawiera D⁰←D* bez zrekonstruowanego tagu;
  10–13 % b ma w stożku tor zgodny z π_s wg maski etapu 6). (c) = (a) z ptcone − pT(π_s), gdy π_s (tor o pT zgodnym w 1 %) leży w stożku. (i) inkluzywne.
  am/bm: to samo po odjęciu torów maski bez prawdy (Δm < 0,155 GeV, ładunek −q_K); maska usuwa 99,6 % π_s leżących w stożku.
- Statystyka: fd_nd 60 500 tag / 195 637 inkl.; mb_nd 1426 tag (100 % zdarzeń) / 3375 inkl. (tory otoczenia dla 40 % zdarzeń). Udział tagu w inkluzywnych:
  0,29–0,34 (fd_nd), 0,16–0,18 (mb_nd). Błędy dwumianowe z n_eff; (a) i (c) są skorelowane (te same wiersze), błąd c/a zawyżony.

Pełne tabele (6 cięć × 4 biny × 2 próbki): `tables.md`; liczby: `iso_pis.json`.

## Wydajności, ptcone40 < 2 GeV

| próbka | pT | ε_a tag | ε_b reszta | ε_c tag bez π_s | ε_i inkl. | a/b | a/i | c/a | c/b |
|---|---|---|---|---|---|---|---|---|---|
| fd_nd | 3–5 | 0,765 ± 0,002 | 0,759 ± 0,002 | 0,803 ± 0,002 | 0,761 | 1,008 ± 0,004 | 1,006 ± 0,004 | 1,049 | 1,058 ± 0,004 |
| fd_nd | 5–8 | 0,630 ± 0,003 | 0,645 ± 0,002 | 0,703 ± 0,003 | 0,640 | 0,976 ± 0,006 | 0,984 ± 0,006 | 1,116 | 1,090 ± 0,007 |
| fd_nd | 8–30 | 0,400 ± 0,005 | 0,465 ± 0,004 | 0,524 ± 0,005 | 0,443 | 0,861 ± 0,014 | 0,903 ± 0,014 | 1,310 | 1,128 ± 0,015 |
| fd_nd | 3–30 | 0,671 ± 0,002 | 0,689 ± 0,001 | 0,732 ± 0,002 | 0,684 | 0,973 ± 0,003 | 0,981 ± 0,003 | 1,091 | 1,062 ± 0,003 |
| mb_nd | 3–5 | 0,791 ± 0,015 | 0,794 ± 0,010 | 0,820 ± 0,014 | 0,795 | 0,995 ± 0,022 | 0,994 ± 0,022 | 1,037 | 1,032 ± 0,022 |
| mb_nd | 5–8 | 0,647 ± 0,023 | 0,694 ± 0,016 | 0,742 ± 0,021 | 0,678 | 0,933 ± 0,039 | 0,955 ± 0,040 | 1,146 | 1,069 ± 0,039 |
| mb_nd | 8–30 | 0,409 ± 0,035 | 0,486 ± 0,026 | 0,525 ± 0,035 | 0,472 | 0,842 ± 0,085 | 0,866 ± 0,086 | 1,284 | 1,081 ± 0,094 |
| mb_nd | 3–30 | 0,693 ± 0,012 | 0,727 ± 0,008 | 0,755 ± 0,011 | 0,719 | 0,954 ± 0,020 | 0,963 ± 0,020 | 1,089 | 1,039 ± 0,020 |

Pozostałe cięcia, a/b (3–5 / 5–8 / 8–30), fd_nd | mb_nd:

| cięcie | fd_nd a/b | mb_nd a/b | fd_nd c/a |
|---|---|---|---|
| ptcone40 < 4 GeV | 1,005 / 1,004 / 0,969 (± 0,002–0,009) | 0,99 / 0,98 / 0,93 (± 0,01–0,06) | 1,013 / 1,034 / 1,115 |
| ptcone40/pT < 1,0 | 1,008 / 1,006 / 0,998 (± 0,002–0,004) | 0,99 / 0,98 / 0,99 (± 0,01–0,02) | 1,013 / 1,010 / 1,009 |
| ptcone40/pT < 0,5 | 1,009 / 0,999 / 0,996 (± 0,004–0,007) | 0,99 / 0,96 / 0,95 (± 0,02–0,04) | 1,054 / 1,055 / 1,060 |
| ptcone100 < 8 GeV | 1,055 / 1,042 / 1,010 (± 0,004–0,012) | 0,95 / 0,91 / 0,87 (± 0,02–0,07) | 1,020 / 1,042 / 1,105 |
| ptcone100/pT < 1,0 | 1,110 / 1,053 / 1,028 (± 0,009) | 0,99 / 0,92 / 0,90 (± 0,05) | 1,077 / 1,072 / 1,047 |

Odpowiednik podwójnego stosunku P1 (odniesienie = drzewo H bez otoczenia, a/b = 1,004 / 1,011 / 1,022): dla ptcone40 < 2 GeV
1,004 / 0,966 / 0,842 (fd_nd) — wobec 1,37 / 1,45 / 1,44 dla sieci z torami otoczenia. Kierunek przeciwny (π_s psuje izolację, a sieci pomaga) i efekt
kilkukrotnie mniejszy.

## Czego ten test nie wyjaśnia

- **ε_c/ε_b ≠ 1.** Po odjęciu π_s D⁰ z tagiem są w fd_nd bardziej izolowane niż reszta o 6–13 % (ptcone40 < 2 GeV) i 8–12 % (ptcone100 < 8 GeV).
  Nie jest to widmo pT w binie (po przeważeniu do widma b: a/b 1,011 / 0,978 / 0,862), nie nieoznakowane D* w b (c/b2 = 1,04 / 1,07 / 1,11) i nie udział
  hadronów b (tag 15–17 %, reszta 14–16 %; c/b osobno dla prompt 1,06 / 1,09 / 1,14). Przyczyna nieustalona.
- **fd_nd i mb_nd różnią się dla stożków 100**: a/b dla ptcone100 < 8 GeV (3–30) 1,039 ± 0,003 wobec 0,927 ± 0,019. Różnica próbki tagu od reszty
  poza π_s zależy więc od próbki MC (skład produkcji / zdarzenia), a nie od pionu. Dla stożka 0,4 obie próbki są zgodne w błędach mb_nd.
- Wniosek z pkt 3 „Krótko” (kompensacja) jest więc własnością MC, nie konstrukcji: część pionowa (c/a) jest pewna, część c/b nie.

## SF: efekt drugiego rzędu

Model: „dane” = MC, w którym aktywność spoza π_s (u = ptcone − π_s; w próbce inkluzywnej u = ptcone − tory maski) jest przesunięta tak, by SF próbki tagu
miał zadaną wartość; π_s bez zmian (tak samo modelowany). Dwa przesunięcia: mnożnik u → k·u (puste stożki zostają puste) i piedestał u → u + δ (wszystkie
kandydaty). Wynik: SF_inkl/SF_tag, ptcone40 < 2 GeV, fd_nd:

| pT | model | SF_tag 0,98 | 0,95 | 0,90 | 0,76 | 0,51 |
|---|---|---|---|---|---|---|
| 3–5 | k·u / u + δ | 0,999 / 1,001 | 0,997 / 1,002 | 0,994 / 1,006 | 0,981 / 1,020 | 0,914 / 1,052 |
| 5–8 | k·u / u + δ | 0,998 / 1,001 | 0,996 / 1,004 | 0,990 / 1,010 | 0,975 / 1,035 | 0,922 / 1,129 |
| 8–30 | k·u / u + δ | 0,996 / 1,003 | 0,994 / 1,012 | 0,990 / 1,026 | 0,969 / 1,078 | 0,905 / 1,266 |
| 3–30 | k·u / u + δ | 0,999 / 1,002 | 0,997 / 1,005 | 0,993 / 1,011 | 0,979 / 1,035 | 0,916 / 1,128 |

- SF_tag ≈ 0,98 (427929): |SF_inkl/SF_tag − 1| ≤ 0,5 % we wszystkich binach i cięciach w fd_nd (mb_nd: ≤ 2 %, w granicach szumu; 8–30 GeV ma 198 zdarzeń tagu).
- SF_tag 0,76–0,51 (428770 i jego klasy): rozrzut między modelami 4–11 % przy 0,76 i 14–36 % przy 0,51, o przeciwnych znakach. Dla ptcone40/pT < 0,5:
  0,98–1,03 przy 0,76, 0,92–1,12 przy 0,51; dla ptcone100: 0,97–1,02 i 0,94–1,05.
- Wkład pierwszego rzędu w próbce inkluzywnej (wydajność rekonstrukcji π_s inna w danych): 10 % utraty π_s zmienia ε_i o +0,2 / +0,5 / +1,1 % — pomijalne.
- **Ograniczenie:** to oszacowanie rozmiaru, nie pomiar. Zakłada, że przesunięcie aktywności jest takie samo dla D⁰ z D* i pozostałych oraz że da się je
  opisać jednym parametrem; nota 6.2.3–6.2.4a pokazuje, że w 428770 żaden jednoparametrowy opis nie domyka izolacji. Rzetelnie rozstrzyga tylko pomiar
  w danych: SF na ptcone po odjęciu toru maski (Δm < 0,155 GeV, ładunek −q_K; definicja bez prawdy, usuwa 99,6 % π_s) w próbce tagu.

## Wniosek

- Zastrzeżenie jest realne co do mechanizmu (π_s jest w ptcone u ≥ 95 % D⁰ z tagiem i obniża ε(ptcone40 < 2 GeV) o 5–31 %), ale dla SF w 427929
  (SF ≈ 1) nie ma znaczenia: < 0,5 % wobec błędu 4 %.
- ε_MC z próbki tagu nie nadaje się jako ε inkluzywna dla ptcone40 < 2 GeV w 8–30 GeV (−10 %) — nota używa ε_MC z tagu tylko w mianowniku SF, więc liczby
  SF to nie zmienia; ε inkluzywną trzeba brać z sygnału inkluzywnego MC.
- Dla 428770 (SF 0,4–0,8) przeniesienie SF z tagu na inkluzywne D⁰ ma dodatkową niepewność rzędu ±10 % przy SF ≈ 0,5 (do +27 % w 8–30 GeV), ±4 % przy
  SF ≈ 0,76, o nieustalonym znaku. Nie zmienia to wniosku Sec. 6.2 (MC nie opisuje izolacji w 428770).
- D⁺: nie badane (brak odpowiednika π_s w stożku; przeniesienie SF na D⁺ pozostaje założeniem).

## Proponowane zdanie do Sec. 6.2 (zamiast ogólnego zastrzeżenia)

> **D\* tag and the soft pion (MC test, 19.11).** The ptcone branches exclude only the D⁰ daughters, so the soft pion of the D\* lies inside the R = 0.4 cone
> for 95–100 % of the tagged D⁰ and contributes a median of 0.31 / 0.47 / 0.80 GeV (pT 3–5 / 5–8 / 8–30 GeV). On MC it lowers ε(ptcone40 < 2 GeV) of the
> tagged sample by 5 / 10 / 24 %; because tagged D⁰ are otherwise more isolated than the rest, the tagged efficiency differs from the inclusive one by
> +0.6 / −1.6 / −9.7 % (forced-decay MC; −1 / −4 / −13 % with ±2–9 % in minimum-bias MC), and by less than 1 % for ptcone40/pT < 0.5. For the scale
> factor the pion cancels to first order: shifting the non-pion activity in MC until the tagged SF equals 0.98 changes the inclusive SF by less than 0.5 %,
> so the 427929 scale factors apply to inclusive D⁰. For the 428770 values the transfer is model dependent: −2 to +8 % at SF = 0.76 and −10 to +27 % at
> SF = 0.51, depending on whether the excess activity is modelled as a scale or as a pedestal. The 428770 numbers of Secs. 6.2.1–6.2.4 therefore remain
> scale factors for D⁰ from D\*; the conclusion that the MC does not describe isolation in that run is unaffected. D⁺ was not tested.

## Pliki

- `iso_pis.py` — obliczenie (numpy; uruchomione interaktywnie na lxplus9103, 75 s, bez condora), `mk_tables.py` → `tables.md`, `iso_pis.json`, `run.log`.
- Na lxplus pozostawiony katalog `/tmp/klyczak_isopis` (lxplus9103: `iso_pis.py`, `iso_pis.json`, `run.log`; < 1 MB). Nic nie zapisano na EOS/AFS.
- Nie zmieniano: ORCH, etapy, nota, podsumowanie.

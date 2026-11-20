# K1d — obciążenie fitu masy D⁺ (H vs H+DLZ) na MC minbias z pełnym tłem (18.11.2026, projekt B)

**Werdykt wg `../KRYTERIA_K1d.md` (zapisane przed pierwszym fitem): obciążenie H+DLZ „w normie” NIE wykazane — wszystkie trzy warianty formalnie POZA NORMĄ.**
C1 (|b − 1| zgodne z 5 % w 1σ) nie przechodzi dla H+DLZ w 3–5 GeV (b = 1,09 ± 0,04) i 5–8 GeV (1,12 ± 0,06) — o włos (przekroczenie 0,001 i 0,009) i tylko na ~2σ od 1;
H przechodzi wszędzie, ale w 5–8 GeV ma tę samą wartość centralną (1,11 ± 0,07) i przechodzi wyłącznie większym błędem. C2 (Δ = DLZ − H zgodne z 0 w 2σ) przechodzi wszędzie przy WP70.
Żadna komórka nie ma mocy (σ_b = 4–8 % > 5 %): MC minbias ma statystykę ≈ runu 427929, więc test nie może być ostrzejszy niż sam pomiar. H zostaje nominałem.

Pliki: `tables.md` (wszystko), `posthoc.md`, `verdict.json`, `out/k1d.json`, `out/weights.json`, `scripts/{k1d,report,posthoc_fake}.py`, `logs/`. Liczenie 25 min (limit 40).

## 1. Wejście (zadanie 1)
- `base.npz` mb: 10 585 186 kandydatów, 2 768 510 zdarzeń, wagi 1. Prawda: 0/1 = D⁺ z c/b (8 460; 4 627 / 2 602 / 1 231 w binach pT), 2–7 = częściowo dopasowany charm (17 195; podkody nieudokumentowane lokalnie),
  8–15 = kombinatoryka (podkody nieudokumentowane; 12 = 90 %). mb nie brało udziału w treningu (sygnał = fd, tło = pasma danych) → score niezależny.
- Po WP70: N_true MC 3 274 / 1 843 / 849 (H) wobec N_fit danych 3 032 / 1 603 / 631 (427929) i 4 142 / 3 180 / 6 205 (428770). ε sygnału mb przy progach z fd: 0,69–0,71 (zgodne z 0,70).
- S/B prawdy MC vs fit danych 427929: stosunek 0,82–1,30 (wszystkie biny „ważne”; H+DLZ 8–30 GeV 1,297 — na granicy). Wobec 428770: 2,2–3,0× (MC wprost nieważne dla tego runu).

## 2. Konstrukcja
(a) MC wprost; (b1) tło MC (prawda ≥ 2) ważone exp(θ·x), regresja logistyczna dane vs MC tylko w pasmach |Δm| > 45 MeV (masa do t², logit score H i H+DLZ, ln pT, wyrazy mieszane), sygnał × k tak, by N_H(WP70) = dane;
(b0) to samo bez masy. Fit i kształty jak K1c (`fastfit.py`, 3–5 GeV kształt ustalony, wyżej swobodny). 200 replik bootstrapu po zdarzeniach MC (sparowane), 200 pseudo-eksperymentów Poissona na komórkę, pol2/3/4, WP 50/70/80.
Domknięcie przeważenia: pasma po cięciu MC/dane 0,99–1,01, χ²/n 0,9–1,6; **liczba kandydatów w oknie piku (nieużywana w dopasowaniu wag) MC/dane 0,95–1,04** (błąd stat. MC 1–4 %) — przy WP50/70/80, dla obu klasyfikatorów i obu runów.
S/B H+DLZ po przeważeniu (przewidywanie): 0,82–1,08 × dane.

## 3. Obciążenie b = N_fit/N_true, WP70 (zadanie 3)
| wariant | pT | pol2: H / H+DLZ | **pol3: H / H+DLZ** | pol4: H / H+DLZ | Δ(DLZ − H), pol3 |
|---|---|---|---|---|---|
| (a) MC wprost | 3–5 | 1,00 ± 0,05 / 1,10 ± 0,04 | **1,04 ± 0,05 / 1,09 ± 0,04** | 1,05 ± 0,05 / 1,09 ± 0,04 | +0,05 ± 0,04 |
| | 5–8 | 1,05 ± 0,05 / 1,13 ± 0,05 | **1,11 ± 0,07 / 1,12 ± 0,06** | 1,12 ± 0,08 / 1,13 ± 0,06 | +0,01 ± 0,05 |
| | 8–30 | 0,97 ± 0,05 / 0,98 ± 0,04 | **0,98 ± 0,07 / 1,01 ± 0,05** | 0,98 ± 0,07 / 1,01 ± 0,06 | +0,03 ± 0,04 |
| (b1) → 427929 | 3–5 | 1,01 ± 0,06 / 1,14 ± 0,05 | **1,06 ± 0,06 / 1,12 ± 0,05** | 1,06 ± 0,06 / 1,12 ± 0,05 | +0,06 ± 0,05 |
| | 5–8 | 1,05 ± 0,05 / 1,14 ± 0,05 | **1,11 ± 0,07 / 1,12 ± 0,06** | 1,12 ± 0,08 / 1,13 ± 0,07 | +0,02 ± 0,05 |
| | 8–30 | 0,98 ± 0,06 / 0,99 ± 0,05 | **0,98 ± 0,08 / 1,02 ± 0,06** | 0,98 ± 0,08 / 1,02 ± 0,07 | +0,04 ± 0,05 |
| (b1) → 428770 | 3–5 | 1,03 ± 0,13 / 1,29 ± 0,10 | **1,13 ± 0,14 / 1,26 ± 0,11** | 1,14 ± 0,14 / 1,26 ± 0,11 | +0,13 ± 0,11 |
| | 5–8 | 1,24 ± 0,16 / 1,54 ± 0,16 | **1,51 ± 0,22 / 1,52 ± 0,18** (57 % udanych replik) | 1,64 ± 0,21 / 1,56 ± 0,18 | 0,00 ± 0,19 |
| | 8–30 | 0,84 ± 0,10 / 0,98 ± 0,09 | **0,86 ± 0,12 / 0,99 ± 0,11** | 0,84 ± 0,12 / 0,99 ± 0,11 | +0,13 ± 0,08 |

- b0 ≡ b1: różnica b ≤ 0,003 (427929) i ≤ 0,013 (428770) → wynik nie zależy od masowej części wagi (brak kołowości przez interpolację).
- Pull width: (a) bootstrap 0,91–1,08; pseudo-eksperymenty 0,89–1,08 (427929) i 0,80–1,10 (428770; 0,80 w 5–8 GeV) → błąd fitu wiarygodny, jak w K1c.
- (b1) → 427929 to praktycznie ta sama próbka co (a) (wagi 0,86–1,16) — nie jest niezależnym potwierdzeniem. (b1) → 428770 to ta sama fluktuacja/struktura tła MC przemnożona przez 2–3× większe B/S
  (fałszywy sygnał w jednostkach B w oknie jest identyczny jak w (a)); wagi 3,5–16 (p99 do 62) zwiększają błąd do 11–22 %; χ²/ndf ważonych histogramów 3–21 odzwierciedla statystykę MC, nie model.

## 4. Co fit bierze za sygnał (leave-one-out wg prawdy, WP70, pol3, wariant (a); % N_true)
| pT | klas. | b − 1 | sam kształt sygnału | charm częściowy 2–7 (razem; największe kody) | kombinatoryka: kod 12 / 8 / 9 / 11 | N_fit na samym tle |
|---|---|---|---|---|---|---|
| 3–5 | H | +4,2 | −1,3 | −0,1 (kod 2 +0,3; 4 −0,4; 6 −0,4) | +5,1 / +1,8 / +0,1 / −0,4 | +5,6 (1,0 % B okna) |
| 3–5 | H+DLZ | +9,3 | −0,8 | +0,5 (2 +0,9; 3 +0,7; 4 −0,6) | +8,5 / +2,4 / −1,2 / +0,8 | +10,1 (2,9 % B) |
| 5–8 | H | +10,5 | −0,7 | +0,7 (4 +1,6; 6 −0,9) | +3,3 / +1,6 / +2,6 / +6,9 | +6,4 (4,7 % B) |
| 5–8 | H+DLZ | +11,7 | +0,5 | −0,5 (4 +1,6; 6 −1,7) | +4,8 / +3,5 / +2,6 / +4,4 | +7,4 (7,1 % B) |
| 8–30 | H | −2,3 | −3,8 | −4,0 (2 −1,8; 4 +1,7; 6 −3,5) | 0,0 / +1,2 / −0,2 / −0,7 | +1,0 |
| 8–30 | H+DLZ | +1,1 | −3,7 | 0,0 (4 +2,4; 6 −2,8) | −0,1 / +2,2 / −0,2 / −0,3 | +5,3 |
- **Fałszywy sygnał to kombinatoryka (kody 12, 8, 11), nie składowe piko-podobne z charmu**: kody 2–7 stanowią 3–26 % tła w oknie, ale ich wkład netto do N_fit to |≤ 1 %| (−4 % dla H w 8–30 GeV); kod 4 daje +1,6…+2,4 %, kod 6 −1…−3,5 % (znoszą się).
- Fit odzyskuje prawdziwy sygnał z dokładnością −1 % (3–5, 5–8) i −4 % (8–30 GeV; kształt fd vs mb / ogony).

## 5. Zysk „prawdziwy” na MC (zadanie 4), pol3/pol3, WP70, g = 1 − [σ/N_true]_DLZ/[σ/N_true]_H [%]
| | 3–5 | 5–8 | 8–30 | średnio | dane K1c |
|---|---|---|---|---|---|
| (a) | +14,2 ± 0,5 | +8,3 ± 4,6 | +14,0 ± 4,3 | **+12,2 ± 2,1** | 427929: 19,2 / 8,3 / 4,4 / **10,7 ± 2,0** |
| (b1) → 427929 | +14,9 ± 0,6 | +8,3 ± 4,8 | +14,3 ± 4,8 | **+12,5 ± 2,3** | jw. |
| (b1) → 428770 | +17,1 ± 0,6 | +20,4 ± 15,2 | +7,4 ± 10,0 | **+15,0 ± 6,1** | 428770: 29,7 / 9,3 / 5,8 / **14,9 ± 1,9** |
Zysk liczony względem N_true (bez premii za zawyżone N_fit) zgadza się ze średnim zyskiem z danych w obu runach. W definicji danych (σ/N_fit) MC daje +14,8 ± 1,9 % (a) i +21,8 ± 5,1 % (428770), czyli
~3–7 pkt zysku w tej definicji to większe N_fit H+DLZ; w 3–5 GeV dane (19 % i 30 %) są powyżej zysku prawdziwego MC (14–17 %). Przy WP80 zysk prawdziwy średnio +9,4 ± 2,4 % (dane 15,6 %), przy WP50 +8,6 ± 1,7 % (6,3 %).

## 6. POST HOC (po obejrzeniu wyników; nie zmienia werdyktu) — `posthoc.md`
- Δ(DLZ − H) w 3–5 GeV vs WP, pol3, (a): WP50 +0,069 ± 0,027 (2,5σ), WP70 +0,050 ± 0,041 (1,2σ), **WP80 +0,159 ± 0,047 (3,4σ)**; w 5–8 i 8–30 GeV ≤ 1,5σ. Kryterium C2 było zapisane tylko dla WP70 i tam przechodzi,
  ale obraz po WP (cięcia zagnieżdżone, więc skorelowane) wskazuje, że w 3–5 GeV H+DLZ ma realnie większe obciążenie dodatnie niż H.
- Fałszywy sygnał na samym tle MC w 3–5 GeV: H+DLZ +3,8 / +2,9 / +2,0 % B okna (WP50/70/80; 1,5–2,4σ), H −1,4 / +1,0 / −0,1 %. Stały ułamek B (nie maleje jak 1/√B) pasuje do rzeźbienia kombinatoryki pod pikiem
  przez H+DLZ (por. K2: R_sculpt 0,93 → 1,03), nie do fluktuacji. Przy S/B danych: +9 % N_S w 427929 i ~+26 % w 428770 (3–5 GeV, WP70), jeśli efekt przenosi się na dane.
- 5–8 GeV: nadmiar wspólny dla obu klasyfikatorów (H +10,5 ± 6,8 %, H+DLZ +11,7 ± 5,9 %; Δ ≈ 0), głównie kod 11 i 12 — 1,5–2σ, może być fluktuacją tej jednej próbki MC.

## 7. Ograniczenia
1. **Statystyka MC ≈ statystyka 427929**: σ_b 4–8 % na komórkę, więc „|b − 1| < 5 %” da się co najwyżej nie obalić, nie wykazać. Bootstrap opisuje rozrzut wokół tej jednej realizacji MC — b zawiera jej fluktuację.
2. **428770**: MC minbias nie opisuje próbki jetowej (S/B 2–3× niższe, wagi do 16 w 8–30 GeV, n_eff/n 0,73); wynik (b1) → 428770 to ekstrapolacja struktury tła MC na inne S/B, z błędem 11–22 % i niestabilnym fitem w 5–8 GeV. Nie jest pomiarem obciążenia w tym runie.
3. **Kołowość**: waga dopasowana tylko w pasmach, gładka (rzędu 2) w masie; b0 (bez masy) daje to samo, a liczba kandydatów w oknie zgadza się z danymi na 1–5 % — interpolacja nie steruje wynikiem. Zostaje: k sygnału pochodzi z N_H fitu danych
   (efekt drugiego rzędu) oraz założenie, że udział składowych piko-podobnych w tle danych jest jak w MC (przeważenie ich nie koryguje osobno).
4. Tło pod pikiem jest „prawdziwe” tylko w sensie MC: jeśli rzeźbienie kombinatoryki w danych różni się od MC (ε_B dane/MC = 0,7–0,9), obciążenie w danych może być inne.
5. Podkody prawdy 2–7 i 8–15 opisane tylko numerami (brak lokalnej dokumentacji). Wagi przeważenia nie były dopasowywane ponownie w replikach. Kształt swobodny w 5–8 i 8–30 GeV (jak K1c).
6. Kryteria: brak odstępstw; §6 jest post hoc. Próg C1 przekroczony minimalnie — wynik formalny „POZA NORMĄ” jest statystycznie słaby (2σ), a „w normie” niewykazalne tą próbką.

## 8. Wniosek
Niezależna od pasm metoda działa (domknięcie w oknie 1–5 %, brak zależności od przeważenia) i daje: (i) zysk prawdziwy H+DLZ ≈ +12 % (427929) / +15 % (428770), zgodny z danymi; (ii) fałszywy sygnał pochodzi z kombinatoryki, nie z refleksji;
(iii) obciążenie H+DLZ w 3–5 GeV +9 ± 4 % (H +4 ± 5 %), z poszlaką post hoc (3,4σ przy WP80), że różnica DLZ − H jest realna. Warunek „obciążenie w normie” nie jest spełniony → **H zostaje nominałem dla D⁺**.
Do rozstrzygnięcia potrzeba ≥ 4× większego MC minbias (σ_b ≈ 2–3 %) albo kontroli na D⁰ z tagiem D*.


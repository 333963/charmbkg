# FEAT_charmnet — czy dodatkowe wejścia CharmNet poprawiają XGB H? (raport końcowy, 16.11.2026)

**Pytanie.** Czy cechy, których CharmNet używa, a XGB H (18 cech) nie — izolacja (ptcone40/60/80/100), absImpact, |z0| córek,
stosunki pT córek, α Armenteros (D⁰); maxdR, min pT córki, mKπ/m_D (Dalitz), dE/dx córek (D⁺) — podnoszą Zproxy = ε_S/√ε_B
i czy są dobrze modelowane w danych.

**Krótko.**
- Izolacja daje na MC +17–21 % Zproxy, ale jest źle modelowana (C2ST H 0,578 → 0,774 dla D⁰, 0,589 → 0,784 dla D⁺) → **NO-GO** bez kalibracji na danych.
- dE/dx (D⁺) nic nie daje (+0,7 %) i jest bardzo źle modelowane (C2ST 0,915) → **twarde NO-GO**.
- **D⁺ Dalitz (mKπ/m_D): GO** — +9–11 % na MC, modelowanie bez zmian (C2ST +0,001), bez rzeźbienia, po przeniesieniu na tło z danych +17–20 % (§3).
- Pozostałe pojedyncze grupy mają zysk poniżej progu szumu (≤ 4 %). Suma cech dobrze modelowanych dla D⁰ (H+ALL_noiso):
  +7 / +5 / +5 % na MC, +6 / +4 / +3 % po przeniesieniu na dane → **GO tylko warunkowo (granicznie)**.

## 0. Ustawienia

- Skrypty: `scripts/study.py` (trening i bootstrap), `modelling.py` (C2ST), `dstar_iso*.py` (izolacja sygnału z tagu D*),
  `sculpt.py`, `transfer.py`. Logi: `logs/`, wyniki: `out/`.
- XGB z hiperparametrami H, 5 foldów po zdarzeniach, cap 300k, bootstrap Poissona po zdarzeniach (200 replik, sparowany
  między wariantami → błąd ΔZ/Z_H jest sparowany).
- Test **mc**: sygnał fd (truth 0/1), tło mb truth ≥ 8 w oknie |Δm| < 55/45 MeV.
- Test **data**: tło treningowe i ε_B = dalekie SB danych 428770 (> 100/110 MeV, z wetami).
- Progi WP z sygnału MC per bin pT.
- **Szum treningu** (H_seed2 vs H): 0,2–1,3 % (D⁰), 0,1–2,5 % (D⁺; 8–30 do 2,5 %). Przyjęty szum: 1,5 %. **Próg zysku: > 4,5 %**.
- Kryterium GO (wszystkie trzy warunki):
  1. zysk > 3× szum;
  2. C2ST(dane SB vs MC tło, zestaw) ≤ C2ST(H) + 0,01;
  3. R_sculpt zgodne z H (w granicach ~2σ).

## 1. Zproxy na MC (test mc), ΔZ/Z_H [%] ± bootstrap (sparowany)

### D⁰

| wariant | WP70 3–5 | WP70 5–8 | WP70 8–30 | WP50 3–5 | WP80 3–5 |
|---|---|---|---|---|---|
| H (Z) | 3,52 ± 0,04 | 3,88 ± 0,07 | 3,87 ± 0,18 | 5,78 ± 0,08 | 2,68 ± 0,03 |
| H_seed2 | −0,2 ± 0,4 | +1,2 ± 0,9 | −1,3 ± 2,4 | +0,1 | +0,2 |
| +iso2 (ptcone40,100) | **+17,0 ± 0,8** | **+18,5 ± 2,1** | **+14,1 ± 4,4** | +16,8 | +17,3 |
| +iso4 | +16,8 ± 0,8 | +20,6 ± 2,1 | +13,7 ± 4,4 | +18,0 | +17,0 |
| +absimp | +1,3 ± 0,4 | +1,9 ± 0,8 | +3,7 ± 2,4 | +3,4 | +1,2 |
| +trkz0 | +1,5 ± 0,4 | +3,2 ± 0,9 | +1,1 ± 2,4 | +1,9 | +1,2 |
| +ptratio | +3,1 ± 0,5 | +1,6 ± 0,8 | −1,3 ± 2,2 | +2,5 | +2,9 |
| +apalpha | −0,5 ± 0,3 | +0,3 ± 0,7 | −1,8 ± 2,1 | 0,0 | −0,4 |
| +ALL_noiso | **+7,3 ± 0,6** | **+5,4 ± 1,1** | +4,6 ± 2,7 | +8,5 | +6,6 |
| +ALL | +25,6 ± 1,1 | +23,6 ± 2,1 | +16,2 ± 4,6 | +27,6 | +24,0 |

### D⁺

| wariant | WP70 3–5 | WP70 5–8 | WP70 8–30 | WP50 3–5 | WP80 3–5 |
|---|---|---|---|---|---|
| H (Z) | 6,68 ± 0,09 | 5,57 ± 0,13 | 5,16 ± 0,25 | 9,97 ± 0,18 | 4,98 ± 0,05 |
| H_seed2 | −0,4 ± 0,6 | −0,7 ± 1,2 | +2,5 ± 3,2 | +0,8 | −0,1 |
| +iso2 | **+20,1 ± 1,2** | **+21,4 ± 2,3** | **+26,1 ± 6,2** | +25,1 | +21,7 |
| +iso4 | +21,1 ± 1,2 | +21,5 ± 2,4 | +29,4 ± 6,7 | +25,2 | +21,4 |
| +absimp | +3,9 ± 0,7 | +1,7 ± 1,2 | +4,8 ± 3,0 | +5,8 | +4,1 |
| +maxdR | +4,2 ± 0,6 | +0,8 ± 1,1 | −0,4 ± 2,9 | +1,8 | +5,1 |
| +ptratio (dMinPt + 3 stosunki) | +4,4 ± 0,8 | +1,7 ± 1,3 | +0,8 ± 3,4 | +5,6 | +5,9 |
| +dalitz (mKπ₁,₂/m_D) | **+10,3 ± 0,9** | **+11,4 ± 1,7** | +9,2 ± 4,7 | +12,8 | +11,0 |
| +pid (dE/dx 3 córek) | +0,7 ± 0,6 | −0,4 ± 1,1 | +0,4 ± 2,8 | +1,0 | +0,8 |
| +ALL_noiso (z pid!) | +25,5 ± 1,2 | +16,4 ± 2,2 | +23,0 ± 5,8 | +27,0 | +24,9 |
| +ALL | +53,5 ± 2,1 | +45,5 ± 4,0 | +57,5 ± 10,9 | +59,7 | +49,4 |

**Uwagi.**
- ALL > suma pojedynczych grup — cechy wzmacniają się nawzajem, zwłaszcza w D⁺: dalitz + ptratio + absimp + maxdR ≈ +23 %,
  a ALL_noiso daje +25 %.
- pT 8–30: błędy statystyczne są 2–5 %, więc nic tam nie jest rozstrzygające.

## 2. Test na danych (sygnał MC vs tło z dalekich SB danych), WP70, ΔZ/Z_H [%]

| gatunek | wariant | 3–5 | 5–8 | 8–30 |
|---|---|---|---|---|
| D⁰ | H_seed2 | +1,1 ± 0,6 | +0,1 ± 1,1 | +0,9 ± 1,6 |
| D⁰ | +ALL_noiso | +10,3 ± 0,9 | +4,5 ± 1,5 | +0,2 ± 1,8 |
| D⁰ | +iso2 | +83 ± 3 | +78 ± 6 | +56 ± 7 |
| D⁺ | H_seed2 | 0,0 ± 0,6 | −0,2 ± 0,9 | +0,2 ± 1,2 |
| D⁺ | +ALL_noiso | +153 ± 4 | +150 ± 8 | +159 ± 11 |
| D⁺ | +iso2 | +99 ± 3 | +99 ± 6 | +88 ± 8 |

**Ten test jest zawyżony dla każdej źle modelowanej cechy.** Klasyfikator uczy się wtedy „MC vs dane” zamiast „sygnał vs tło”.
- Izolacja: sygnał MC jest bardziej izolowany niż tło danych, ale także niż prawdziwy sygnał w danych (§3).
- D⁺ ALL_noiso: +150 % to artefakt dE/dx (C2ST 0,915; §5).
- Jedyny wiarygodny wiersz to D⁰ ALL_noiso: +10 / +4,5 / 0 %, zgodny z MC.

## 3. Przeniesienie na dane (`transfer.py`, model z testu mc), WP70

Definicje:
- ε_B = dalekie SB danych, ocenione modelem MC (OOF).
- ε_S^dane = ε_S MC przeważone w 5 klasach ptcone40 stosunkiem f_dane/f_MC, mierzonym na tagu D*
  (`out/dstar_iso_fit.json`; S ≈ 200; klasa ptcone40 < 1 GeV: dane 0,19 ± 0,07 vs MC 0,50).
- Niepewność ε_S^dane od frakcji D* (toy, 2000 prób, przedział 16–84 %) policzona osobno: `scripts/toy_rw.py`.

### D⁰

| wariant | bin | ε_S MC | ε_S^dane (toy 68 %) | ε_B (SB danych) | Z_transfer/Z_H |
|---|---|---|---|---|---|
| H | 3–5 / 5–8 / 8–30 | 0,70 | 0,720 / 0,723 / 0,715 | 0,032 / 0,024 / 0,036 | 1 |
| +iso2 | 3–5 | 0,70 | 0,648 [0,631; 0,666] | 0,012 | **1,44** |
| +iso2 | 5–8 | 0,70 | 0,642 [0,625; 0,665] | 0,009 | 1,43 |
| +iso2 | 8–30 | 0,70 | 0,659 [0,648; 0,673] | 0,020 | 1,24 |
| +ALL_noiso | 3–5 / 5–8 / 8–30 | 0,70 | ≈ H | 0,028 / 0,022 / 0,034 | **1,06 / 1,03 / 1,03** |
| +absimp / trkz0 / ptratio / apalpha | — | — | ≈ H | — | 0,999–1,04 |

### D⁺

Błąd bootstrapu Z_transfer: 0,6–2,9 % (względny).

| wariant | ε_S^dane 3–5 / 5–8 / 8–30 | ε_B (SB danych) 3–5 / 5–8 / 8–30 | Z_transfer/Z_H 3–5 / 5–8 / 8–30 |
|---|---|---|---|
| H | 0,707 / 0,705 / 0,704 | 0,0069 / 0,0072 / 0,0153 | 1 |
| H_seed2 | ≈ H | ≈ H | 0,991 / 1,001 / 1,026 |
| +iso2 | 0,615 [0,58; 0,65] / 0,578 [0,53; 0,63] / 0,607 [0,58; 0,64] | 0,0022 / 0,0024 / 0,0058 | 1,55 / 1,44 / 1,40 (na samym MC-S: 1,79 / 1,75 / 1,63) |
| +dalitz | ≈ H | 0,0048 / 0,0053 / 0,0111 | **1,20 / 1,17 / 1,18** |
| +ptratio | ≈ H | 0,0060 / 0,0068 / 0,0137 | 1,07 / 1,03 / 1,05 |
| +maxdR | ≈ H | 0,0058 / 0,0069 / 0,0145 | 1,09 / 1,03 / 1,03 |
| +absimp | ≈ H | 0,0064 / 0,0069 / 0,0156 | 1,04 / 1,03 / 0,99 |
| +pid | ≈ H | 0,0068 / 0,0073 / 0,0150 | 1,01 / 1,00 / 1,01 |
| +ALL_noiso | 0,709 / 0,703 / 0,704 | 0,0035 / 0,0048 / 0,0102 | **1,41 / 1,23 / 1,22** |
| +ALL | 0,617 / 0,580 / 0,610 | 0,0012 / 0,0017 / 0,0039 | 2,12 / 1,72 / 1,71 |

**Wnioski dla D⁺.**
- Model MC zastosowany do tła z danych daje dla Dalitz **większy** zysk niż na tle MC: +17–20 % (MC: +9–11 %).
- ALL_noiso: +22–41 %. Model MC prawie nie korzysta z dE/dx (pid sam: 1,00–1,01), więc ten zysk pochodzi z cech
  dobrze modelowanych. Zastrzeżenie: ε_S^dane nie jest korygowane na dE/dx sygnału.
- ε_S^dane dla cech nieizolacyjnych ≈ ε_S MC — korekta izolacji praktycznie nie działa na te modele.

### Interpretacja izolacji

Po przeważeniu ptcone40 izolacja nadal daje formalnie +24–44 % (D⁰). **Tej liczbie nie należy ufać**, z trzech powodów:
- przeważenie jest jednowymiarowe (tylko ptcone40, 5 klas) i opiera się na S ≈ 200 z dwóch plików;
- ptcone100 — gorzej modelowane (C2ST 1D 0,76) — nie jest korygowane;
- korelacje izolacji z innymi cechami nie są korygowane.

Wynik pokazuje jednak, że **izolacja ma potencjał** także w danych: tło w danych jest jeszcze mniej izolowane niż w MC
(ε_B spada 2,5×). Warunkiem jest kalibracja izolacji sygnału na danych (D* tag / sWeights, większa statystyka v2.1).

## 4. Rzeźbienie tła (`sculpt.py`)

**Metoda.** R = ε_B(okno) obserwowane / ε_B(okno) z liniowego dopasowania ε_B(m) w SB.
- „mc” = mb truth ≥ 8, okno sygnału.
- „data” = pasmo [okno, 100/110 MeV) danych, z dopasowaniem w dalekich SB.

Wyniki: WP70, modele z testu mc (`out/*_mc/sculpt.json`).

| gatunek | wariant | R (MC truth8, SR) | R (dane, pasmo) | nachylenie ε_B [/100 MeV] |
|---|---|---|---|---|
| D⁰ | H | 0,947 ± 0,011 | 0,912 ± 0,020 | −0,25 |
| D⁰ | +iso2 | 0,955 ± 0,013 | 0,894 ± 0,031 | −0,25 |
| D⁰ | +absimp / trkz0 / ptratio / apalpha | 0,941–0,958 | 0,913–0,929 | −0,22…−0,25 |
| D⁰ | +ALL_noiso | 0,966 ± 0,012 | 0,912 ± 0,021 | −0,20 |
| D⁰ | +ALL | 0,963 ± 0,014 | 0,914 ± 0,033 | −0,21 |
| D⁺ | H | 0,890 ± 0,013 | 0,946 ± 0,015 | −0,23 |
| D⁺ | +dalitz | 0,953 ± 0,015 | 0,970 ± 0,018 | −0,26 |
| D⁺ | +iso2 | 0,909 ± 0,016 | 0,928 ± 0,027 | −0,22 |
| D⁺ | +absimp / maxdR / ptratio / pid | 0,897–0,922 | 0,929–0,956 | −0,22…−0,26 |
| D⁺ | +ALL_noiso | 0,958 ± 0,018 | 0,987 ± 0,020 | −0,28 |
| D⁺ | +ALL | 0,924 ± 0,021 | 1,001 ± 0,036 | −0,27 |

**Wniosek.** Żadna grupa nie rzeźbi bardziej niż H. Wszystkie R są w granicach ~2σ od H albo bliżej 1 niż H.
Lekkie niedopasowanie R < 1 występuje już w H — to krzywizna ε_B(m), nie efekt nowych cech.

Dalitz (mKπ/m_D zawiera m_D) nie wprowadza piku: R = 0,953 przy 0,890 dla H, nachylenie bez zmian.
Warianty z modelami testu data (`out/*_data/sculpt.json`) dają te same wnioski.

## 5. Modelowanie (C2ST tło: dane SB vs MC mb, `out/modelling_*.json`)

| zestaw | D⁰ C2ST | Δ vs H | D⁺ C2ST | Δ vs H |
|---|---|---|---|---|
| H | 0,578 | — | 0,589 | — |
| +iso2 / iso4 | 0,774 / 0,775 | **+0,196** | 0,784 / 0,785 | **+0,195** |
| +absimp | 0,581 | +0,003 | 0,591 | +0,001 |
| +trkz0 | 0,580 | +0,002 | — | — |
| +ptratio | 0,580 | +0,001 | 0,594 | +0,004 |
| +apalpha | 0,579 | +0,000 | — | — |
| +maxdR | — | — | 0,590 | +0,001 |
| +dalitz | — | — | 0,590 | +0,001 |
| +pid (dE/dx) | — | — | **0,915** | **+0,326** |
| +ALL_noiso | 0,584 | +0,006 | 0,911 (pid) | +0,32 |
| +ALL | 0,775 | +0,196 | 0,919 | +0,33 |

### Per cecha (C2ST 1D / KS)

**Izolacja.** ptcone40 0,69 / 0,30; ptcone60 0,72; ptcone80 0,74; ptcone100 0,76 / 0,40 — oba gatunki.
- Mediana ptcone40 tła: dane 6,3 vs MC 2,4 GeV (D⁰); 5,8 vs 2,3 GeV (D⁺).
- Sygnał z tagu D*: f(ptcone40 < 1 GeV) dane 0,19 ± 0,07 vs MC 0,49.

**dE/dx D⁺.**
- C2ST 1D jest niskie (0,50–0,52), ale KS 0,29–0,71: rozkłady różnią się wartościami dyskretnymi / domyślnymi
  (DistPerp = 0, flagi Valid).
- Łącznie z resztą zestawu: 0,915. Wygląda to na **różnicę w kodowaniu wartości domyślnych / kalibracji dE/dx między
  MC a danymi**, nie fizykę.

**Pozostałe nowe cechy:** C2ST 1D 0,50–0,54 — tyle samo, co cechy H (0,50–0,55).

## 6. Werdykt

Kryterium: zysk WP70 > 4,5 % (3–5 lub 5–8) **i** ΔC2ST ≤ 0,01 **i** brak rzeźbienia.

| gatunek | grupa | zysk MC (3–5 / 5–8) | ΔC2ST | rzeźbienie | werdykt |
|---|---|---|---|---|---|
| D⁰ | izolacja | +17 / +19 % | +0,196 ✗ | nie | **NO-GO** (potencjał po kalibracji na danych) |
| D⁰ | absImpact | +1,3 / +1,9 % ✗ | +0,003 | nie | NO-GO (brak zysku) |
| D⁰ | \|z0\| córek | +1,5 / +3,2 % ✗ | +0,002 | nie | NO-GO (brak zysku) |
| D⁰ | stosunki pT | +3,1 / +1,6 % ✗ | +0,001 | nie | NO-GO (brak zysku) |
| D⁰ | α Armenteros | ~0 ✗ | 0,000 | nie | NO-GO |
| D⁰ | razem bez izolacji | +7,3 / +5,4 % ✓ (dane: +6 / +3 %) | +0,006 ✓ | nie ✓ | **GO warunkowe** — formalnie spełnia kryterium, ale zysk ≈ 4–5× szum, a na danych +3–6 % |
| D⁺ | izolacja | +20 / +21 % | +0,195 ✗ | nie | **NO-GO** (potencjał po kalibracji) |
| D⁺ | dE/dx | +0,7 / −0,4 % ✗ | +0,326 ✗ | nie | **NO-GO (twarde)** |
| D⁺ | Dalitz mKπ/m_D | **+10,3 / +11,4 %** ✓ (dane: +20 / +17 %) | +0,001 ✓ | nie ✓ (R lepsze niż H) | **GO** |
| D⁺ | stosunki pT + min pT | +4,4 / +1,7 % (dane: +7 / +3 %) | +0,004 | nie | NO-GO samodzielnie; wzmacnia Dalitz (ALL_noiso) |
| D⁺ | maxdR | +4,2 / +0,8 % ✗ (dane: +9 / +3 %) | +0,001 | nie | NO-GO samodzielnie |
| D⁺ | absImpact | +3,9 / +1,7 % ✗ | +0,001 | nie | NO-GO samodzielnie |

**Brakujący wariant.** Nie policzono **D⁺ „ALL_noiso_nopid”** (dalitz + ptratio + absimp + maxdR), bo ALL_noiso zawiera
dE/dx. Pid sam daje +0,7 % (MC) i +1 % (transfer), a ALL_noiso +16–26 % (MC) i +22–41 % (transfer). To sugeruje podobny zysk bez pid, ale wymaga potwierdzenia C2ST (dla ALL_noiso 0,911 przez pid). To najważniejszy następny bieg:
```
study.py --species Dplus --test mc,data
```
z nowym wariantem w `GROUPS` oraz `modelling.py`, `sculpt.py` i `transfer.py`.

## 7. Rekomendacje

### Projekt A (XGB H, analiza po odsłonięciu)

Zmiana klasyfikatora po odsłonięciu **wymaga osobnej decyzji**. Proponowana kolejność, jeśli zapadnie decyzja „tak”:
1. D⁺: dodać Dalitz (mKπ₁,₂/m_D), ewentualnie ze stosunkami pT. To realny zysk: +10 % Zproxy na MC i +17–20 % na tle z danych
   (≈ +20–40 % efektywnej statystyki), dobrze modelowany, bez rzeźbienia. Komplet bez izolacji i dE/dx
   (dalitz + ptratio + maxdR + absimp) może dać +20–40 % — po biegu kontrolnym (§6).
2. D⁰: dodatkowe cechy dają +3–7 % — nie warto zmieniać klasyfikatora tylko dla nich.
3. Izolacji i dE/dx nie dodawać. Izolacja może wrócić dopiero po pomiarze rozkładu izolacji sygnału w danych
   (tag D*, sWeights, v2.1), z niepewnością systematyczną ε_S.

## Pliki
- Zproxy: `out/{D0,Dplus}_{mc,data}/zproxy.csv`
- Rzeźbienie: `out/*/sculpt.json`, logi `logs/sculpt_*.log`
- Transfer: `out/transfer_{D0,Dplus}.csv`, `logs/transfer_*.log`
- Modelowanie: `out/modelling_*.json`
- Izolacja z tagu D*: `out/dstar_iso_fit.json`

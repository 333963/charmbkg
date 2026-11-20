# FITBIAS_Dplus — nominalny model fitu masy D⁺→K⁻π⁺π⁺ (projekt B)

20.11.2026, ~02:15 CET. **PRELIMINARY, w całości POST HOC** (uzyski M0, p fitu 3–5 GeV, domknięcia FITBIAS i sygnały pozorne BKGCAT znane przed zapisem kryteriów).
Kryteria: `KRYTERIA.md` — §1–7 zapisane przed jakimkolwiek fitem, §8 (wynik rankingu) przed fitami okna sygnału danych. Tylko piki inkluzywne, WP70 drzewa H, bez AFP. Praca lokalna; z lxplus tylko odczyt (`bin/lx -pull` 78 plików `bkgcat/parts`), bez condora.
Kod fitu: `FITBIAS_2026-11-19/skrypty/fb_fit.py` bez zmian (import); szablony wprost z drzewa H (`BKGCAT/in/bkgcat.npz`). Kontrola mechaniki: A0 odtwarza uzyski M0 z FITBIAS/XSEC co do cyfry (3547 / 1651 / 687; 4118 / 3451 / 6458).
Liczba porównań: 9 ustawień × 3 biny (MC: domknięcie pełne, połówki, Asimov, rozkład) + 6 komórek zabawek × 200 + 6 komórek danych; jedna decyzja; bez poprawki na wielokrotność.

## Odpowiedź

1. **Obecny model (A0 = M0) na MC, WP70:** b = N_fit/N_prawda − 1 = **+4,8 ± 8,0 / +11,4 ± 6,4 / −0,8 ± 7,2 %** (3–5 / 5–8 / 8–30 GeV; średnio +5,7 ± 4,1 %). Kryterium |b| ≤ 5 % albo ≤ 1σ_MC nie jest spełnione w 5–8 GeV (1,8σ).
2. **Żadne z 9 ustawień nie spełnia kryterium obciążenia (K1)** — bin 5–8 GeV ma +9…+16 % we wszystkich ustawieniach bez szablonu PartComb, a ustawienie z szablonem PartComb (A5) jest poprawne tylko w próbce (poza próbką −14…−38 %).
   Zgodnie z zapisanym przypadkiem „nikt nie spełnia K1” wybór wg Q: **A2 = szablon sygnału z ograniczeniem Gaussa (dmu 0,5 MeV; sc 0,10) + szablon DMM + CM wprost z drzewa H przy WP70 (κ_MC × 1, związany z N_S) + exp(Czebyszew 3), 1,70–2,20 GeV** (odpowiednik M4 z FITBIAS z szablonem wprost).
   Q = 0,074 wobec 0,090 (A0); obciążenie A2: +4,6 ± 7,7 / +9,2 ± 6,2 / −0,9 ± 7,3 %.
3. **Dane 427929 (wynik główny):** 3547 → **3404 ± 292 (stat) ± 142 (model)**; 1651 → **1630 ± 102 ± 21**; 687 → **651 ± 53 ± 36**. Zmiana −4,0 / −1,3 / −5,2 %; p(χ²) = 0,034 / 0,89 / 0,23.
4. **D⁺ 3–5 GeV, p = 0,014–0,04:** niezgodność **nie siedzi w oknie sygnału** (χ² 20/18), tylko w pasmach: dołek −3,0σ w 1,800–1,825 GeV (tuż pod pikiem) oraz −2,1σ / +2,1σ / −2,7σ w 2,00–2,15 GeV. Nie naprawia jej rząd tła ani szablon; zawężenie zakresu daje p = 0,10–0,12.
   Ta sama komórka w 428770 (3× więcej kandydatów) ma p = 0,85, a mb_nd ma p = 0,05–0,08 przy resztach rozproszonych → fluktuacja pasm, nie brakująca składowa. Skutek uboczny: dołek pod pikiem ciągnie szerokość (sc = 1,16 ± 0,07) i uzysk zależy od niej na 13 % (kształt stały 3136, wolny 3547).
5. **Wpływ:** σ(D⁺) i D⁺/D⁰ (D⁰ bez zmian) **−4,0 / −1,3 / −5,2 %** względem obecnych uzysków; wszystko w granicach σ_stat. Systematyka modelu 4,2 / 1,3 / 5,5 %; osobno niepewność domknięcia (fit_bias) 9,0 / 11,1 / 7,3 %, ograniczona statystyką jednej realizacji mb_nd.

## 1. Test na MC (mb_nd, pełne tło, WP70) — `tabele/mc.md`, `out/mc.json`

Prawda sygnału = tryb ekskluzywny K⁻π⁺π⁺ (+ radiacyjny): N_prawda = 3245 / 1827 / 843 (Σ kluczy `Sig*` = `mb|sig` z WPSCAN). κ(DMM + CM) = 0,527 / 0,534 / 0,507; κ(PartComb) = 10,2 / 4,0 / 2,0.

| ustawienie | b pełne [%]: 3–5 | 5–8 | 8–30 | połówki (poza próbką), średnia: 3–5 | 5–8 | 8–30 | Asimov poza próbką b̄ [%] |
|---|---|---|---|---|---|---|---|
| **A0** wolny, bez szablonu (nominał) | +4,8 ± 8,0 | **+11,4 ± 6,4** | −0,8 ± 7,2 | +4,5 ± 7,2 | +13,5 ± 6,3 | −1,9 ± 7,1 | 4,4 |
| A1 Gauss | +4,9 ± 7,8 | +12,6 ± 6,9 | −0,8 ± 7,4 | +6,0 | +10,9 | −1,0 | 4,1 |
| **A2** Gauss + DMM + CM | +4,6 ± 7,7 | **+9,2 ± 6,2** | −0,9 ± 7,3 | +3,5 ± 7,1 | +6,9 ± 6,1 | −5,7 ± 7,0 | **2,4** |
| A3 = A2, 1,75–2,10 | +6,3 ± 9,1 | +9,2 ± 7,4 | −3,9 ± 8,3 | +7,6 | +4,1 | −8,0 | 2,7 |
| A4 wolny + DMM + CM | +4,5 ± 9,3 (!) | +10,6 ± 6,3 | +5,9 ± 7,5 (!) | +2,9 | +13,4 | +2,8 | 4,2 |
| A5 Gauss + DMM + CM + PartComb | +2,6 ± 7,0 | +2,6 ± 6,0 | −0,0 ± 7,0 | **−17,7** | **−18,6** | **−33,7** | 23,0 |
| A6 Gauss, 1,75–2,10 | +6,3 ± 9,0 | +16,0 ± 7,5 | −1,5 ± 8,4 | +8,5 | +9,0 | −2,2 | 3,8 |
| A7 Gauss + DMM + CM, 1,78–2,05, Czeb. 2 | +7,4 ± 9,8 | +15,9 ± 6,9 | −4,2 ± 7,7 | +7,9 | +10,9 | −7,5 | 4,7 |
| A8 = A2 z Czeb. 4 | +4,4 ± 9,2 | +11,4 ± 6,3 | −1,7 ± 7,4 | +3,3 | +5,5 | −6,5 | 2,9 |

(!) = fit bez dodatnio określonej macierzy (wolny kształt + szablon). σ_b = max(1,3 × HESSE, rms 40 replik Poissona).

**Która składowa (A0, kształt stały, reżim addytywny; N_spur/N_prawda w %):**

| składowa | 3–5 | 5–8 | 8–30 | średnia ważona |
|---|---|---|---|---|
| DMM (D⁺→KKπ, πππ(π⁰), Kπππ⁰ przy dolnej krawędzi) | −1,2 ± 0,4 | −1,5 ± 0,7 | −3,0 ± 1,1 | **−1,4 ± 0,3** |
| CM: D_s (KKπ z błędną masą, πππ) | +1,8 ± 0,5 | +1,4 ± 0,8 | −0,6 ± 1,2 | **+1,4 ± 0,4** |
| CM: D⁰, D*⁰, D*⁺ (K3π …) | +0,1 ± 0,5 | −0,1 ± 0,8 | −0,6 ± 1,2 | 0,0 ± 0,4 |
| PartComb | +3,2 ± 3,4 | **+5,6 ± 2,8** | +3,1 ± 3,3 | **+4,2 ± 1,8** |
| kombinatoryka bez powabu | +3,0 ± 3,3 | +0,9 ± 1,5 | −1,1 ± 1,3 | 0,0 ± 1,0 |
| kombinatoryka z b | −2,4 ± 1,3 | +0,1 ± 1,0 | +1,7 ± 1,8 | −0,4 ± 0,7 |
| suma / b przy kształcie stałym | +4,4 / +5,3 | +6,4 / +4,9 | −0,4 / −2,5 | |
| swoboda kształtu (b własne − b stałe) | −0,5 | **+6,5** (sc = 1,08) | +1,7 | |

- **Obciążenie obecnego modelu to PartComb (+4,2 ± 1,8 %, 2,3σ, spójny znak w 3 binach) plus w 5–8 GeV wolna szerokość (+6,5 %).** DMM i D_s są rzeczywiste (4–5σ), ale przeciwnego znaku i znoszą się (suma 0,0 ± 0,5 %): szablon DMM + CM (A2) niczego dużego nie naprawia — porządkuje model i zmniejsza b̄ Asimova z 4,4 do 2,4 %.
- PartComb per rodzic: 3–5 GeV rodzic D⁺ +3,3 ± 1,5 %; 5–8 GeV D*⁰ +4,4 ± 1,5 i D⁰ +2,7 ± 1,2 %; 8–30 GeV D⁺ +2,5 ± 1,9, D_s +2,3 ± 1,3 % — **żaden rodzic się nie powtarza**; kształty w binach 25 MeV są gładkie, bez piku D* przy 2,01 GeV i bez progu.
  „Próg poniżej piku” z D⁺→Kππ + π⁰ jest w DMM (36 / 14 / 8 kandydatów, < 1,735 GeV), nie w PartComb.
- Połówki (5–8 GeV, A0): +22,9 ± 9,1 % w jednej i +4,2 ± 8,6 % w drugiej — nadwyżka siedzi w połowie plików (różnica 1,5σ). **Struktury PartComb nie da się tą statystyką ani potwierdzić, ani wykluczyć** (jak w BKGCAT).
- **Szablon PartComb (A5) nie działa:** w próbce domyka (+2,6 / +2,6 / 0,0 %, p = 0,99–1,00 = przeuczenie na własnym szumie), poza próbką −14…−38 % — szum szablonu o normalizacji 2–10 × N_S ściąga N_S w dół. Wygładzony szablon PartComb byłby nieodróżnialny od wielomianu tła.
- Szum szablonu DMM + CM (A2) poza próbką, szablon z połowy statystyki: −0,8 ± 0,7 / −2,2 ± 1,0 / −5,0 ± 1,3 % (`out/mc_oos.json`); przy pełnej statystyce szablonu oczekiwane ≈ połowa (−0,4 / −1,1 / −2,5 %) — mieści się w σ_model, zapisane w json jako `template_stat_bias_rel`.

## 2. Zabawki przy czystości danych — `tabele/select.md`, `out/toys.json`

Tło z pasm danych WP70 (okno ±45 MeV zaślepione), N_inj = znany uzysk M0, wstrzyknięty sygnał fd(dmu_t, sc_t) + κ(DMM + CM); sc_t ~ N(1,05; 0,10) (priorytet fitu w 1,00 — celowe niedopasowanie). S/B w oknie: 427929 0,17 / 0,74 / 1,29; 428770 0,065 / 0,24 / 0,68.

| ustawienie | średnie obciążenie % | mediana σ_toy % | rms pull (błąd 1,3 × HESSE) | rms(N)/HESSE | K2 (komórek) | Q |
|---|---|---|---|---|---|---|
| A0 | −1,6 (8–30: −4,3 / −3,6) | 7,9 | 1,01–1,17 | 1,39 | 3/6 | 0,090 |
| A1 | −2,7 | 7,0 | 0,93–1,22 | 1,31 | 2/6 | 0,082 |
| **A2** | −2,1 (−0,5…−3,6) | 7,0 | 0,92–1,05 | 1,32 | 2/6 | **0,074** |
| A4 | −0,4 | 7,9 | 0,98–1,11 | 1,42 | **6/6** | 0,089 |
| A8 | −2,1 | 7,2 | 0,87–1,18 | 1,33 | 2/6 | 0,078 |
| A5 | −16,5 | 6,0 | — | — | 0/6 | 0,238 |

- **Błąd HESSE trzeba skalować × 1,3** (rms(N)/HESSE = 1,3–1,4 we wszystkich ustawieniach; po skalowaniu rms pull ≈ 1) — to samo co dla M8.
- Ustawienia Gaussa mają w zabawkach −2 % z założonego przesunięcia prawdziwej szerokości (1,05 wobec priorytetu 1,00); A0/A1 dodatkowo −4 % przy 8–30 GeV (brak szablonu DMM + CM). K2 spełnia tylko A4 (wolny kształt + szablon), ale ma większy rozrzut i niestabilną macierz na MC.
- Ranking Q (zapisany przed danymi): A2 < A8 < A1 < A3 < A4 < A0 < A6 < A7 ≪ A5. Prostsze od A2 w 5 %: brak → **A2**.

## 3. Dane — `tabele/data.md`, `out/data.json`, **`out/dplus_yields.json`**

σ_stat = max(1,3 × HESSE, bootstrap po zdarzeniach); σ_model = max(|N − N_A0|, |N − N_A4|); fit_bias = |b| ⊕ σ_b z domknięcia mb_nd (nie w σ_model; uzysk nie jest poprawiany).

| run | bin | N obecny (A0) | p | **N nowy ± stat ± model** | ΔN/N | p | sc | fit_bias (MC) |
|---|---|---|---|---|---|---|---|---|
| 427929 | 3–5 | 3547 ± 276 | 0,041 | **3404 ± 292 ± 142** | −4,0 % | 0,034 | 1,11 ± 0,04 | 9,0 % |
| 427929 | 5–8 | 1651 ± 76 | 0,92 | **1630 ± 102 ± 21** | −1,3 % | 0,89 | 0,98 ± 0,03 | 11,1 % |
| 427929 | 8–30 | 687 ± 59 | 0,35 | **651 ± 53 ± 36** | −5,2 % | 0,23 | 0,86 ± 0,02 | 7,3 % |
| 428770 | 3–5 | 4118 ± 315 | 0,85 | **4098 ± 424 ± 20** | −0,5 % | 0,85 | 0,98 ± 0,05 | 9,0 % |
| 428770 | 5–8 | 3451 ± 199 | 0,31 | **3313 ± 206 ± 139** | −4,0 % | 0,15 | 0,96 ± 0,02 | 11,1 % |
| 428770 | 8–30 | 6458 ± 216 | 0,051 | **6508 ± 254 ± 91** (A8) | +0,8 % | 0,016 | 1,01 ± 0,01 | 7,6 % |

- **Bramka K3 (p ≥ 0,01):** A2 przechodzi w 5 z 6 komórek. **428770 8–30 GeV odrzuca A2 (p = 0,004)**; następne w rankingu A8 (ten sam model z Czeb. 4) ma p = 0,016 → użyte w tej jednej komórce (σ_model zawiera |A8 − A2|). To ten sam obraz co dla D⁰ M8: skład/kształt z minimum-bias MC słabo przenosi się na próbkę dżetową przy wysokim pT. Dla 427929 nominał A2 bez wyjątków.
- Kolumna „obecny ± σ” = max(HESSE, bootstrap), jak w XSEC. Nowy błąd statystyczny jest większy o czynnik 1,3 z zabawek (nie dlatego, że fit jest gorszy: HESSE A2 jest mniejsze niż A0 w 3–5 GeV: 224 vs 257).
- **427929 8–30:** szerokość 0,86 ± 0,02 (w 428770: 1,01 ± 0,01); fit wolny A0 stał w minimum lokalnym (sc = 0,92, χ² 98; A1/A4 znajdują 0,85, χ² 95). Bootstrap (53) > HESSE (38) — komórka niestabilna (stąd też „skok +5,3 %” z BKGCAT). Kształt stały: 701 (+7,8 %).
- Pozostałe ustawienia (opisowo): rozrzut A0–A4, A6–A8 w 427929: 3362–3547 / 1596–1654 / 627–687.

### D⁺ 3–5 GeV: gdzie jest zły χ² (427929; A0 118/93, A2 119/93)

| obszar [GeV] | 1,70–1,78 | 1,78–1,825 | okno 1,825–1,915 | 1,915–2,00 | 2,00–2,05 | 2,05–2,20 |
|---|---|---|---|---|---|---|
| χ²/biny (A2), Σ(dane − model) | 17,5/16 (+50) | **15,6/9 (−165)** | 20,0/18 (+136) | 15,4/17 (+23) | 10,8/10 (−120) | **40,2/30** (+76) |
| to samo, 428770 | 16,5/16 | 7,1/9 | 10,1/18 | 12,5/17 | 9,7/10 | 22,9/30 |

- Reszty w binach 25 MeV (A2): **−3,1σ w 1,800–1,825**, +1,8σ w 1,825–1,850, −2,1σ w 2,000–2,025, +2,2σ w 2,075–2,100, −2,6σ w 2,125–2,150; reszta |·| ≤ 1,2.
- Nie pomaga: Czeb. 4 (p = 0,033), szablon DMM + CM (0,034), kształt stały (0,028). Zawężenie zakresu: 1,75–2,10 → p = 0,12; 1,78–2,05 → 0,10 (odcina fluktuacje > 2,07 GeV, dołek pod pikiem zostaje).
- **Brak piku D*⁺→D⁰π⁺ przy 2,010 GeV:** dodatkowy Gauss (σ = 8 MeV) woli normalizację ujemną w 427929 3–5 (−166 kandydatów, −2Δln L = 6,0), 428770 5–8 (−163; 7,1) i 428770 8–30 (−128; 10,4); w pozostałych 3 komórkach ≈ 0. mb_nd tego nie ma (3–5: +137, 4,6; reszta 0).
  Trzy lokalne niedobory 2,4–3,2σ w jednym z góry wskazanym miejscu to więcej niż przypadek pojedynczej komórki — **otwarte** (możliwe rzeźbienie masy przez score w danych inaczej niż w MC); wpływ na N_S ≤ 2,8 %, w σ_model.
- **Wniosek:** p = 0,014–0,04 pochodzi z pasm (dołek tuż pod pikiem + rozrzut powyżej 2 GeV), nie z kształtu sygnału i nie z brakującej składowej pikującej; w 428770 i w MC tej struktury nie ma. Ostrzeżenie o uzysku jest jednak zasadne z innego powodu: dołek pod pikiem poszerza dopasowany sygnał (sc 1,16 wolny, 1,11 z ograniczeniem) i **uzysk 3–5 GeV zależy od szerokości na 13 %**. A2 (ograniczenie Gaussa) daje 3404, tj. −4 %; R = σ_cięcia/σ_WP70 z XSEC_atlas rośnie z 0,77 do ≈ 0,80 ± 0,12.

## 4. Werdykt

- **Nominał D⁺: A2** — sygnał: szablon fd_nd z ograniczeniem Gaussa (dmu 0 ± 0,5 MeV; sc 1 ± 0,10); składowa związana: DMM + CM z mb_nd wprost z drzewa H przy WP70 (`BKGCAT/out/peak_templates_H.npz`, klucze `Dplus|<bin>|WP70|{h5,kap}`; κ × 1); tło exp(Czebyszew 3); 1,70–2,20 GeV; biny 5 MeV. Błąd = max(1,3 × HESSE, bootstrap).
- Końcowe uzyski 427929: **3404 ± 292 ± 142; 1630 ± 102 ± 21; 651 ± 53 ± 36** (stat ± model). Systematyka modelu **4,2 / 1,3 / 5,5 %**. Do budżetu przekroju dodatkowo niepewność domknięcia fit_bias 9,0 / 11,1 / 7,3 % (decyzja o złożeniu: XSEC_final; nie jest to pomiar obciążenia, tylko granica z jednej realizacji MC).
- Wpływ na σ(D⁺) i na D⁺/D⁰ względem obecnych uzysków: **−4,0 / −1,3 / −5,2 %** (D⁰ bez zmian). Przekroju nie przeliczałem; plików `XSEC_final_2026-11-20/` nie dotykałem.
- **Czy σ(D⁺) przestaje być „scenariuszem”:** od strony fitu — model jest ustalony procedurą jak dla D⁰ i ma jawną systematykę; ale **kryterium domknięcia nie zostało spełnione** (5–8 GeV +9,2 ± 6,2 %). Scenariusz SF = 1 ± 10 % to osobna sprawa (SFABS / XSEC_atlas), której to zadanie nie zmienia.

## Zastrzeżenia
- Post hoc; jedna realizacja mb_nd (statystyka ≈ run 427929): domknięcie ma błąd 6–8 % i nie rozróżnia ustawień A0–A4, A6–A8; ranking rozstrzyga Asimov poza próbką (b̄ 2,4 vs 4,4 %) przy σ_toy 7–8 % — przewaga A2 nad A0 w Q to 18 %, nad A1 10 %.
- K1 niespełnione przez nic, K2 tylko przez A4; wybór wg zapisanego przypadku awaryjnego. Zbiór „spełniających K1” do σ_model jest pusty → σ_model z A0 i A4 (doprecyzowane w KRYTERIA §8 przed danymi).
- Szablon DMM + CM i κ z tego samego mb_nd co domknięcie pełne (w próbce); poza próbką szum szablonu daje −0,4…−2,5 % (oszacowanie). BR w MC ≠ PDG dla trybów szablonu (D⁺→K⁺K⁻π⁺π⁰ × 4,5; π⁺π⁺π⁻ × 1,47) — nie przeważałem (BKGCAT: wpływ na uzyski ≤ 0,7 %, poza niestabilną komórką 8–30); D_s nie ma w TruthDmesons.
- Zabawki: tło z tej samej rodziny co model (nie testują tła pod pikiem); czystość ustawiona znanym uzyskiem M0; środek sc_t = 1,05 to założenie — przy sc_t = 1,00 ustawienia Gaussa nie miałyby obciążenia −2 %.
- Ograniczenie Gaussa (0,5 MeV; 10 %) pochodzi z PEAK_pos (ten sam pik inkluzywny). W 427929 dane ciągną sc do 1,11 (3–5) i 0,86 (8–30), w 428770 jest 0,96–1,01.
- 428770 8–30: A2 odrzucone (p = 0,004), A8 przechodzi ledwo (0,016); komórka ma słabe p także w A0 (0,051). Dla 428770 liczby są drugorzędne (próbka dżetowa).
- Pierwsza wersja rozkładu „poza próbką” dla składowych modelowanych (`mc.json`: `dec.AtoB/BtoA`) była błędna (fit bez tła, N_B na granicy 0, sztuczne −2…−5 %); zastąpiona przez `fd_mc_oos.py` / `mc_oos.json`. Na wybór nie wpływała (nie wchodzi do Q).
- 30 replik bootstrapu (błąd rms ±13 %); 40 replik Poissona w σ_b.

## Czego nie dało się ustalić
- Czy nadwyżka PartComb (+4,2 ± 1,8 %) jest strukturą — potrzeba ≥ 4× większego mb_nd (albo niezależnej realizacji); do tego czasu jest w fit_bias.
- Źródła niedoboru przy 2,01 GeV w trzech komórkach danych (brak w MC).
- Rzeczywistej skali szerokości D⁺ w danych z niezależnego źródła (D* nie daje D⁺).

## Pozostałości
- `in/parts/` (78 plików, 8,6 MB — kopia `/eos/user/k/klyczak/charm_q/bkgcat/parts`), `in/half_A`, `in/half_B` (dowiązania), `in/bkgcat_half{A,B}.npz`; `skrypty/__pycache__`. Na lxplus i EOS nic nie zapisano.

## Pliki
`KRYTERIA.md`, `RAPORT.md`, `NOTA_tekst_EN.md`; `skrypty/`: `fd_lib.py` (ustawienia, szablony, fit), `fd_mc.py`, `fd_mc_oos.py`, `fd_toys.py`, `fd_select.py`, `fd_data.py`, `fd_report.py`;
`tabele/`: `mc.md`, `select.md`, `data.md`; `out/`: `mc.json`, `mc_oos.json`, `toys.json`, `select.json`, `data.json`, **`dplus_yields.json`** (`cells[]`: run, bin, model, N, sigma_stat, sigma_model, p + pola pomocnicze: N_A0, dN_rel_vs_A0, fit_bias_rel, closure_b, sc, N_fixshape …).

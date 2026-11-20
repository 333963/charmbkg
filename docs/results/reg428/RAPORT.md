# REG428 — run dżetowy 428770 per region triggera: uzyski, SF, wydajność triggera dżetowego dla zdarzeń z D, przekrój z R3 (projekt B)

20.11.2026, ~00:30 CET. **Całość post hoc** (uzyski WP70, model fitu i wynik 427929 były znane; regiony zdefiniowane w nocie §6.11 przed tą pracą). Tylko bity triggera — bez selekcji protonu offline i bez ξ; wyników `TAGGED_2026-11-19` nie czytano.
Pojedyncze D, drzewo H, WP70. Liczby: `tabele/tabele.md` (T0–T5, generowane), `out/summary.json`.

## 0. Odpowiedź

1. **Uzyski per region (T1).** R3x (`noalg_L1J12 ∧ ¬R1 ∧ ¬R2`) niesie 70–87 % uzysku runu, „reszta” (żaden z trzech łańcuchów; bez świetlności) 12–24 %, R2 2–7 %, R1 < 1 % (zgodnie z XSEC §8.4: ~80 D przy błędzie ~40). Suma regionów rozłącznych = całość w 0–4 %.
   D⁰ 8–30 GeV: szablon DMM + CM (M8) odrzucony w całości (p = 1·10⁻⁷), w R3x (3·10⁻⁴) i w reszcie (3·10⁻⁴), nie w R2 (0,53; mała statystyka) → w tym binie M0 dla całego runu; przy 5–8 GeV M8 przechodzi wszędzie poza resztą (p = 0,009).
2. **SF drzewa z tagu D* (T2).** 8–30 GeV: R3x 0,990 ± 0,050, R2 0,90 ± 0,15, reszta 0,830 ± 0,076; zgodność trzech regionów χ²/ndf = 3,2/2 (p = 0,20); łańcuch J12 minus ¬J12 = +0,150 ± 0,086 (1,7σ). 5–8 GeV: 0,87 ± 0,15 / 0,60 ± 0,21 / 0,89 ± 0,17, p = 0,50. 3–5 GeV i R1: nieokreślone. Różnicy SF drzewa między regionami nie widać w tej precyzji (5–15 %); kierunek (reszta niżej) ten sam co dla sieci w §5.11.
3. **Wydajność triggera dżetowego dla zdarzeń z D (T3).** Metoda z planu (R1 jako tryger ortogonalny) **nie działa**: w R1 jest 20–1126 kandydatów WP70 na komórkę (cała masa), z bitem J12 2–21, sygnału ∧ J12 zgodne z zerem. Zamiast tego odniesienie MBTS w runie 427929 (ten sam `HLT_noalg_L1J12`, preskala 1, ta sama świetlność co R1, ale ~100× więcej D):
   **ε_J12(D) = 2,5 ± 0,8 % (3–5 GeV), 2,7 ± 0,9 % (5–8), 12,0 ± 1,7 % (8–30)** (średnia D⁰ ⊕ D⁺; per gatunek błąd 16–71 %). Rośnie z pT: ~25–38 % przy 12–30 GeV. Na poziomie kandydatów (tło) ułamek jest ~2× większy (3 / 6 / 21 %) — proxy z kandydatów zawyża.
   ε_j20(D | L1_AFP) dla R2: z R1 błąd 65–120 % (T5) — niemierzalna.
4. **Przekrój z R3 (T4).** σ_R3 = N_J12/(C · SF · B · L · ε_J12): D⁰ 5–8: 504 ± 354 µb (427929: 376 ± 31 stat.), D⁰ 8–30: 155 ± 47 (148 tym samym modelem M0; nominał 132 ± 11), D⁺ 8–30: 59 ± 10 (55 ± 3; scenariusz SF). Błąd 17–72 %, **≤ 20 % tylko dla D⁺ 8–30 GeV**.
   To nie jest niezależny pomiar: σ_R3/σ_427929 ≡ ρ = tempo (D ∧ J12)/L w 428770 do tego samego w 427929. ρ = 1,07 ± 0,15 (pT ≥ 5 GeV, 4 komórki, p = 0,93); 0,85 ± 0,12 ze wszystkich sześciu (D⁺ 3–5: 0,49 ± 0,20). Tempo wszystkich zdarzeń J12 na jednostkę świetlności jest w 428770 o 16 % wyższe (T0), więc przeniesienie ε między runami ma nieusuwalną systematykę tego rzędu, a ρ nie rozróżnia 1,00 od 1,16.
5. **Werdykt** — tabela §5.

## 1. Dane i metoda

- **Bity per zdarzenie.** Tabele score etapu 4 (`dnet4/app/d150`, `cls`) mają tylko GRL, MBTS, `noalg_L1J12`; h5 v2.0 nie ma `HLT_j20_L1AFP_A_OR_C`. Bity wzięte z ROOT v2.0 (`charm_raw_v20/data_428770`, 97 plików, 7 989 496 zdarzeń, wszystkie w GRL) i dołączone po `event_uid` = task << 36 | plik << 22 | wpis.
  Walidacja złączenia: 0 niezgodności bitu GRL i J12 względem `cls` na 77 M wierszy kandydatów; w tabeli D* 0 niezgodności numeru zdarzenia, LB, bitów J12 i AFP na 905 844 wierszy. Klasa „wszystkie” i MB427929 odtwarzają histogramy WPSCAN co do wpisu (także 30 replik), uzyski FITBIAS co do cyfry, SF z tagu D* etapu 4 (0,876 ± 0,119; 0,968 ± 0,043).
- **Regiony (rozłączne):** R1 = `HLT_noalg_L1AFP_A_OR_C`; R2 = `HLT_j20_L1AFP_A_OR_C ∧ ¬R1`; R3x = `HLT_noalg_L1J12 ∧ ¬R1 ∧ ¬R2`; reszta. Osobno całe łańcuchy J12 (29,41 nb⁻¹) i j20_L1AFP (28,12 nb⁻¹) oraz warianty bez LB 4–16. Zdarzenia: T0.
  `HLT_j20_L1J12` jest ustawiony w 99,96 % zdarzeń J12; P(L1_AFP przed preskalą | J12) = 0,58 %; P(j20_L1AFP | J12 ∧ L1_AFP) = 95,8 %; P(J12 | R1) = 0,064 %, P(j20 | R1) = 0,93 %.
- **Fity:** maszyneria FITBIAS (`fb_fit.py`, `fb_var.py`, bez zmian): D⁰ M0 i M8 (κ × 1,364), D⁺ M0; σ = max(HESSE, rms 30 replik bootstrapu po zdarzeniach). Test odrzucenia szablonu: p(χ²) M8 < 0,01 → M0 (**próg mój, post hoc**; FITBIAS progu nie zapisał: przyjął 0,038, odrzucił 1·10⁻⁷). Jeden model na bin dla całego runu, wybrany na klasie „wszystkie”.
  Klasy 427929 (MB ∧ J12, MB ∧ ¬J12) nie mają replik — błąd HESSE, który wg FITBIAS zaniża rozrzut o ~25–30 %.
- **SF z tagu D*:** `btrig_common.fit_eps` → `charmtag.dstar_refit.fit_passfail` (tło progowe, przedział z profilu), progi i ε_MC z `dnet4_results.json`; bez odejmowania pasm.
- **Wydajność triggera:** ε = N_pass/(N_pass + N_fail) z fitów piku osobno w zdarzeniach z bitem i bez; sprawdzian z odjęcia pasm (okno ± win, pasma symetryczne za przerwą 10 MeV) w binach analizy i w 8 drobnych binach pT.

## 2. Uzyski i SF per region — T1, T2

Udział R3x w uzysku: D⁰ 73 / 74 / 79 %, D⁺ 70 / 75 / 87 %; reszta 23 / 22 / 14 % i 24 / 20 / 12 %. Uzysk na jednostkę świetlności w łańcuchu J12 jest taki sam z LB 4–16 i bez nich (T4, ostatnia kolumna) — weto ¬R2 nie ma znaczenia liczbowego (R2 ∩ J12 = 0,6 % zdarzeń J12).
D⁺ nie ma próbki tagowanej; SF_rel z inkluzywnego piku per region nie powtarzałem (w R2 i R1 statystyka kategorii luźnej nie wystarcza; w R3x byłby to wynik WPTIGHT dla „J12only”).

## 3. Wydajność triggera dżetowego — T3, T5

- R1 jest regionem L1_AFP: zawiera 0,5 % zdarzeń z D (XSEC §8.4) przy preskali 60, czyli ~10⁻⁴ wszystkich D runu. Do ε_J12 z błędem 20 % trzeba ≥ ~25–50 sygnałowych D z bitem J12 na komórkę przy S/B jak w R1; jest ich ≲ 1.
  Oczekiwane dziś: ~15 D (oba gatunki) × 12 % ≈ 2 przy 8–30 GeV, ~10 × 2,7 % ≈ 0,3 przy 5–8 GeV. **Potrzebna statystyka R1: ≥ ~30× (8–30 GeV) i ≥ ~150× (5–8 GeV)** (szacunek rzędu wielkości) — tyle co zdjęcie preskali 60 albo więcej. I nawet wtedy byłaby to ε dla zdarzeń z D i sygnałem w AFP, nie dla wszystkich D.
- Odniesienie MBTS (427929): ε_J12(D) ma błąd 14 % (8–30, oba gatunki razem), 33–34 % (3–8 GeV). Do 20 % per gatunek: D⁰ 8–30 ×2,1 świetlności minimum-bias (~1 nb⁻¹), D⁰ 5–8 ×12 (~6 nb⁻¹), D⁺ 5–8 ×3,5.
- J12 wybiera 2–3 % zdarzeń z D poniżej 8 GeV i 12 % powyżej, a mimo to niesie ¾ uzysku runu — próbka J12 to wąski, twardy wycinek produkcji D (ekstrapolacja ×8–40 do przekroju inkluzywnego).
- R2: ε_j20(D | AFP) wymaga N_R1 w mianowniku → 65–120 % błędu. Na poziomie kandydatów w R1 bit j20 ma 19 / 42 / 35 % (D⁰), 22 / 43 / 59 % (D⁺) wierszy — to tło, górne oszacowanie (dla J12 kandydaci zawyżają ×2).

## 4. Przekrój z R3 i spójność z 427929 — T4

Wzór i składniki jak w XSEC (`DEFINICJA.md`): C z mb_nd bez wag, SF z tagu D* runu 427929 (D⁰ ≥ 5 GeV; reszta scenariusz SF = 1), B PDG, L = 29,41 ± 0,68 nb⁻¹, ε_J12 z T3. Założenie: C · SF dla zdarzeń z D ∧ J12 jest takie samo w obu runach (wtedy skraca się w ρ).
Nieujęte: przeniesienie ε_J12 między runami (tempo J12: 152,5 vs 131,4 zdarzeń/µb⁻¹, ×1,16; XSEC §8.1: dryf 6 % w 428770, czułość na pociągi paczek), różnica ⟨μ⟩ (0,045 vs 0,005), zależność C od krotności w zdarzeniach dżetowych.
Wynik jest zgodny z 427929 w każdej komórce (największe odchylenie D⁺ 3–5 GeV: −2,6σ; fit MB ∧ J12 na granicy skali szerokości), ale jego moc to ±15 % w najlepszym razie (pT ≥ 5 GeV łącznie) — tyle, ile wynosi efekt, który miałby wykryć.

## 5. Werdykt

| region | L [nb⁻¹] | D w regionie (D⁰ / D⁺, 3–30 GeV) | ilościowo użyteczny? | dlaczego |
|---|---|---|---|---|
| R1 `noalg_L1AFP_A_OR_C` | 0,469 | ~60 ± 30 / ~30 ± 20 | **nie** jako punkt inkluzywny; nie jako odniesienie triggera | mierzy σ × P(L1_AFP \| D) ≈ 0,4 % σ (XSEC §8.4); za mało D o czynnik ≥ 30–150 |
| R2 `j20_L1AFP ∧ ¬R1` | 27,65 | ~610 / ~530 | **nie** dla przekroju (ani inkluzywnego, ani „z AFP”) | potrzebne P(L1_AFP \| D) (±50 %) i ε_j20(D \| AFP) (niemierzalna: 65–120 %); zostaje jako próbka kształtów / SF (0,60–0,90 ± 0,15–0,21) |
| R3 łańcuch `noalg_L1J12` | 29,41 | ~14 500 / ~11 100 | **tylko jako test spójności** przy pT ≥ 8 GeV (D⁺ 17 %, D⁰ 30 %); **nie** jako niezależny przekrój | ε_J12(D) = 2,5–12 % znane tylko z 427929 (14–34 %) i nieprzenoszalne lepiej niż ~16 %; σ_R3/σ_427929 ≡ stosunek temp |
| reszta (inne łańcuchy skimu) | brak | ~3600 / ~2400 | **nie** | brak definicji świetlności |
| cały run | — | 19 000 / 14 000 | **tak, względnie:** kształty, SF z tagu D* (R3x 8–30 GeV: 0,990 ± 0,050), testy modelu fitu, stosunki w obrębie jednej klasy triggera | statystyka 3–10× większa niż 427929; pomiar bezwzględny i D⁺/D⁰ zostają w 427929 |

## 6. Zastrzeżenia

- Post hoc; 84 fity uzysku + 72 fity Δm + ułamki z pasm; bez poprawki na wielokrotność. Próg p < 0,01 dla M8 ustalony po obejrzeniu wyników FITBIAS.
- ε_J12(D) pochodzi z innego runu (inne paczki, ⟨μ⟩ 10× mniejsze); zgodność tempa D ∧ J12 (ρ) jest jedynym sprawdzianem i ma błąd 15–35 %.
- Fity MB ∧ J12 (43–98 zdarzeń sygnału): trzy z sześciu na granicy skali szerokości; błąd HESSE bez bootstrapu. Odjęcie pasm dla D⁰ poniżej 5 GeV jest nieużyteczne (odbicie K↔π).
- σ(427929) w T4 liczone tu z tych samych składników co XSEC, ale bez jego systematyk; dla D⁰ 8–30 modelem M0 (148 µb), nie nominalnym M8 (132 µb).
- D⁰ 3–5 GeV i D⁺: SF = 1 (scenariusz) — liczby σ w tych komórkach nie są pomiarem.
- R1 nie był tu fitowany jako wynik (T1 podaje go tylko dla domknięcia sumy); obowiązują liczby XSEC §8.4 (zgodne: D⁰ 45 ± 30 / 9 ± 10; w 8–30 GeV tu M0: 9 ± 6, tam M8: 5 ± 4).

## 7. Pliki

- `skrypty/`: `r0_probe.py`, `r1_extract.py` (lxplus: bits → hist → merge → dstar → bsum), `r2_fit.py` (lokalnie: yields, dstar, sbeff), `r3_report.py` (tabele).
- `runs/lxplus/`: `reg_hist.npz` (masa 1 MeV × przedział score × klasa; repliki WP70), `dstar.npz` (wiersze D* w oknie masy z bitami), `bits_summary.json`.
- `out/`: `yields.json`, `dstar.json`, `sbeff.json`, `summary.json`; `tabele/tabele.md`; `AKAPIT_DO_NOTY.md`.
- Przebiegi interaktywne na lxplus9103 (łącznie ~7 min, 2 procesy), bez condora, bez zmian w ORCH.
- **Pozostałości:** lxplus9103 `/tmp/klyczak_reg428` (130 MB: `bits/` 98 plików, `hist/` 40 plików, `reg_hist.npz`, `dstar.npz`, skrypty) — nieusunięte; lokalnie `out/dstar.log` z pierwszego, nieudanego startu (błąd importu `charmtag`).

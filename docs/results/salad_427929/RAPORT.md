# SALAD_427929 — SALAD v1.0 na runie 427929 (MBTS) i w klasach trygera 428770

Projekt B (tło podobne do sygnału), pojedyncze D⁰ i D⁺. Stan: **paczka gotowa i przetestowana, joby niewysłane, wyników fizycznych nie ma.**
Ten plik zawiera procedurę i kryteria zapisane 17.11.2026 przed obejrzeniem jakichkolwiek liczb z pełnej próbki.
Liczby z testu (§6) pochodzą z 12 + 1 + 1 + 1 plików w trybie skróconym i służą tylko do sprawdzenia kodu.
Tabele wyników wygeneruje `analyze.py` do `WYNIKI.md` po pobraniu json; wnioski dopisać w §7.

## 1. Pytanie

SALAD v1.0 (MC tła przeważone do danych klasyfikatorem uczonym w pasmach bocznych, 18 cech H, masa jako wejście) był dotąd liczony tylko na runie 428770, gdzie 82 % zdarzeń zapisał tryger dżetowy L1_J12. `B_trig/BKG_427929/RAPORT.md` pokazał na próbce inkluzywnej bez 18 cech, że w 427929 przeważenie zawierające score domyka ε_B pod pikiem (χ² = 13,2/15), ale samego SALAD v1.0 nie uruchomiono. Tu:

1. Czy SALAD v1.0 w konfiguracji produkcyjnej domyka ε_B w czystym minbias (427929, `HLT_mb_mbts_L1MBTS_1`, GRL) w granicach budżetu?
2. Czy wagi wyuczone w jednej klasie trygera przenoszą się na inną i ile różnicy dane/MC pochłania SALAD w próbce dżetowej?
3. Czy czułość na skład tła (B1) jest taka sama przy statystyce 427929?

Zakres: tylko pasma boczne i tło. W obszarze sygnału danych liczona jest wyłącznie interpolacja z pasm bocznych; kandydaci danych z SR nie wchodzą do żadnego histogramu, treningu ani C2ST. Żadnych fitów piku. Pary D–D poza zakresem.

## 2. Procedura (zamrożona w kodzie paczki `Q_lxplus/salad427`, wersja `salad427-1.0`)

**Dane.** h5 v2.0 `/eos/user/k/klyczak/charm_h5_v20/parts/{data_427929, data_428770, mb_nd, fd_nd}`. Selekcja jak B1: pT 3–30 GeV, L_xy ≥ 0,25 mm, okno masy D⁰ 1,75–2,00, D⁺ 1,70–2,20 GeV; D⁺: weto D*; D⁰: weto refleksji |m_swap − m₀| > 25 MeV wszędzie. Klucz zdarzenia: `event_uid`.

**Regiony** (z `salad_v1.yaml`, bez zmian): D⁰ SR < 55 MeV, pasmo walidacyjne 55–75, SB > 75; D⁺ SR < 45, pasmo 45–70, SB > 70 MeV.

**SALAD.** `charmbkg.bkggen.salad_v1.fit_salad_v1` z zamrożonego yaml (plik w wersji 1.1): `hgb`, 18 cech, t = (m − m₀)/0,1, obcięcie wag 0,05–20, `max_n` 400 tys. MC tła: D⁰ truth ≥ 2 bez sygnału 0/1 (`mc_truth_min_species`), D⁺ truth ≥ 8. Wariant kontrolny D⁰: truth ≥ 8.

**Klasy** (bity `trigPassed_*`, definicje jak `B_trig`):

| klasa | próbka | definicja |
|---|---|---|
| MB427929 | 427929 | GRL i `HLT_mb_mbts_L1MBTS_1` |
| J12only | 428770 | `HLT_noalg_L1J12` i żaden z 7 łańcuchów L1MBTS (25 % zdarzeń wg hasza) |
| MBTSbased | 428770 | dowolny łańcuch L1MBTS |
| AFPonly | 428770 | bit AFP (`noalg_L1AFP_A_OR_C` lub `mb_sptrk_pt4_L1AFP_A_OR_C`), bez J12 i bez MBTS |

**Tagger H427.** Jak `b1_dataeff_d0.py`: HGB (200 drzew) fd truth 0/1 vs SB danych, progi WP z fd per bin pT. Różnice wobec B1: 5 foldów po zdarzeniach zamiast parzystości (score OOF), tło treningowe zawsze z SB klasy MB427929, jeden tagger dla wszystkich klas. WP 95/80/70/50, biny pT 3–5/5–8/8–30.

**Per klasa:**
- bieg „prod”: jeden fit na całym SB, bramki v1.0 z yaml (n_eff/N(SR) ≥ 0,5; N_SB ≥ 20 tys.; C2ST pasma ≤ 0,53);
- k-fold: 5 foldów po zdarzeniach, model foldu k uczony na danych SB i MC spoza foldu, wagi dla MC z foldu k; ziarna 0 i 1;
- ε_B: w SB (OOF) i w paśmie walidacyjnym dane zliczone wprost; pod pikiem ε_B(MC ważone, SR) wobec interpolacji danych N_pass/N_all z dopasowań ML exp(pol2) do binów SB (5 MeV);
- błędy: dane — bootstrap Poissona histogramu (400 replik); MC — ważony dwumianowy; systematyka interpolacji — największe odchylenie wariantów {exp(pol1), SB symetryczne, liniowa lewe/prawe}; trening — |R(ziarno 0) − R(ziarno 1)|/√2;
- C2ST (HGB, 3 foldy, 120 tys. na stronę) w SB i paśmie, przed i po, plus podłoga (połówki danych SB);
- n_eff/N w SR, paśmie i SB per bin pT; udział wag obciętych;
- rzeźbienie masy: corr(w, m) w SB, średnia waga SR/pasmo/SB, χ² kształtu masy w binach SB dane vs MC przed i po, N_SR/N_SB danych (pol2) wobec MC;
- obciążenie estymatora: ta sama interpolacja na MC ważonym wobec prawdy w SR.

**Transfer.** Wagi OOF klasy a porównane z danymi klasy b (12 par): R pod pikiem, w paśmie i w SB, C2ST SB i pasma. Dodatkowo C2ST dane–dane w SB (każda klasa 428770 vs MB427929).

**Adwersarialny B1.** Na MC v2.0: pseudo-dane = MC parzystość 1 zniekształcone (tryby `smooth`, `cat_global`, `cat_local`, `cat_step`, `cuts_local`; współczynniki i σ = 35 MeV jak `b1_study.py`), zmniejszone do liczby kandydatów SB klasy MB427929; „MC” = parzystość 0; tagger H427. To nie jest test na danych 427929: sprawdza czułość B1 przy tej statystyce i na h5 v2.0.

## 3. Kryteria (zapisane przed wynikami)

R = ε_B(MC ważone SALAD, SR) / ε_B(dane, interpolacja pol2); σ = stat danych ⊕ stat MC ⊕ syst interpolacji ⊕ rozrzut ziaren.

| nr | pytanie | kryterium |
|---|---|---|
| K1 | bramki v1.0 w klasie | n_eff/N(SR) ≥ 0,5 i N_SB ≥ 20 tys. i C2ST pasma ≤ 0,53 (bieg „prod”). Niespełnienie którejkolwiek = SALAD v1.0 nieakceptowany w tej klasie, niezależnie od R |
| K2 | domknięcie pod pikiem wobec budżetu | budżet B: D⁺ ±5 / ±10 / ±25 %, D⁰ ±5 / ±10 / ±15 % przy WP80/70/50 (D⁰ z truth ≥ 8: ±20 % przy WP50). **ZGODNE**: \|R − 1\| ≤ B i σ ≤ B. **zgodne, niemiarodajne**: \|R − 1\| ≤ B, ale σ > B. **nierozstrzygnięte**: \|R − 1\| > B i \|R − 1\| − B < 2σ. **PRZEKROCZONE**: \|R − 1\| − B ≥ 2σ. WP95 nie ma budżetu: liczba tylko raportowana |
| K3 | warunek użycia v1.0 (B1 §4) | \|R − 1\| ≤ 5 % przy WP80 i ≤ 10 % przy WP70 w binie „all” i w 3–5 GeV, klasa MB427929 |
| K4 | pasmo walidacyjne | \|R_pasmo − 1\| ≤ B przy WP80/70/50 i odchylenie < 2σ; w paśmie jest ogon sygnału (szacunek < 0,5 % kandydatów po WP50), więc R_pasmo < 1 o tyle nie jest niedomknięciem |
| K5 | C2ST | po przeważeniu ≤ 0,53 w SB (OOF) i paśmie. Warunek konieczny, nie test (B1: ślepy na skład) |
| K6 | rzeźbienie masy | \|corr(w, m)\| < 0,05 w SB i \|średnia w(SR)/średnia w(SB) − 1\| < 0,03; χ²/ndf kształtu masy w SB po przeważeniu nie większe niż przed |
| K7 | transfer a → b | **przenosi się**: K2 = ZGODNE przy WP80 i WP70 (bin „all”) i C2ST SB po ≤ 0,53. **przenosi ε_B, nie rozkłady**: K2 zgodne, C2ST > 0,53. **nie przenosi się**: K2 PRZEKROCZONE przy WP80 lub WP70. Reszta: nierozstrzygnięte |
| K8 | ile różnicy pochłania SALAD w próbce dżetowej | miara: n_eff/N(SR) i AUC_SB klasy J12only wobec MB427929. „Pochłania tryger” gdy n_eff/N(J12only) < 0,5 lub AUC_SB(J12only) − AUC_SB(MB427929) > 0,05 |
| K9 | B1 przy statystyce 427929 | tryb przechodzi, gdy \|pred/prawda − 1\| ≤ 10 % przy WP80 i WP70; porównanie z B1 (D⁺ `cat_local` 0,845 / 0,809 / 0,70; D⁰ 0,97 / 0,95 / 0,85). Zmiana czułości = różnica > 2σ_stat prawdy |

Odczyt łączny:
- pytanie 1: SALAD v1.0 potwierdzony w minbias, jeśli MB427929 spełnia K1, K3 i K5, a K2 nie ma werdyktu PRZEKROCZONE w żadnym binie pT przy WP80/70/50;
- wynik K2 „niemiarodajne” lub „nierozstrzygnięte” w binie 8–30 GeV jest spodziewany (mało MC i danych) i nie jest ani potwierdzeniem, ani odrzuceniem;
- przy wielu komórkach (4 klasy × 4 biny × 3 WP) pojedyncze „nierozstrzygnięte” przy 2σ nie zmienia odczytu; liczy się wzór po WP i pT.

Hipotezy, które te liczby mogą obalić: (a) pełny SALAD na 18 cechach domyka 427929 tak jak przeważenie V2 ze score (BKG_427929); (b) wagi z 427929 nie przenoszą się na J12only, bo różnica siedzi w aktywności zdarzenia, której 18 cech nie zawiera.

## 4. Pliki

| plik | rola |
|---|---|
| `Q_lxplus/salad427/s427_lib.py` | stałe, klasy trygera, hash/foldy, interpolacja SB (ML exp(pol-k)), błędy ważone |
| `Q_lxplus/salad427/s427_extract.py` | h5 v2.0 → `parts/s427_<próbka>_<g>.npz` (D⁰ i D⁺, 18 cech, bity trygera) |
| `Q_lxplus/salad427/s427_study.py` | tagger H427, SALAD prod + k-fold per klasa, transfer, B1 → `s427_<gatunek>.json`, `_adv.json`, `_w.npz` |
| `Q_lxplus/salad427/{s427_common.sh, s427_job.sh, s427.sub, make_dag_s427.sh, submit_s427.sh}` | DAG: 39 × ext → 2 × study → fin (batch `charm_salad427`) |
| `Q_lxplus/salad427/src/charmbkg` | kopia pakietu z `~/Projects/charm/charm-bkg` (stan roboczy, `SNAPSHOT.txt`) |
| `ORCH/stages/salad427.sh` | etap ORCH (sync, check, submit, done, fetch → `SALAD_427929/runs/lxplus/`) |
| `SALAD_427929/analyze.py` | lokalnie, < 5 s: json → `WYNIKI.md` z werdyktami K2 |

## 5. Uruchomienie

Wpis do `ORCH/pipeline.yaml`:

```yaml
  - name: salad427    # SALAD427: SALAD v1.0 na 427929 (MBTS) i w klasach trygera 428770, transfer wag, B1 (projekt B) → SALAD_427929/runs/lxplus; analiza SALAD_427929/analyze.py
    needs: [h5_v20]
    max_hours: 24
```

Szacunek: ekstrakcja 39 węzłów po 5–20 min; study D⁰ 2–3 h, D⁺ 4–7 h (8 CPU, 24 GB, flavour `tomorrow`); plus kolejka. Po pobraniu:
`cd ~/Documents/charm && .venv_h5/bin/python charm_work/SALAD_427929/analyze.py` → `charm_work/SALAD_427929/WYNIKI.md`.

## 6. Test (17.11.2026, lxplus993, `/tmp/klyczak_s427test`, bez condora)

- Ekstrakcja: 12 plików 427929 (10 s), 1,2 mln wierszy jednego pliku 428770 (25 s), 1 plik mb_nd (6 s), 1 plik fd_nd (7 s). Kandydaci D⁰ / D⁺: 85 688 / 311 377; 64 050 / 243 116; 39 300 / 135 707; 16 168 / 7 521.
- Study D⁰ w trybie skróconym (`S427_FAST=1`: 2 foldy, 30 drzew, 20 replik, C2ST na 5 tys.): całość 10 min na 2 wątkach, kod wyjścia 0, oba json zapisane; `analyze.py` składa z nich wszystkie tabele.
- Liczby z testu nie są wynikiem: MC tła 38 tys. kandydatów, σ(R) przy WP80 = 0,07, bramka N_SB wyłączona w trybie skróconym.
- Study D⁺ w trybie skróconym: 12 min na 2 wątkach, kod wyjścia 0, oba json zapisane.
- `make_dag_s427.sh` (DAG testowy w `/tmp`): 39 węzłów ext + 2 study + fin; `submit_s427.sh check`: OK (schemat h5 v2.0, środowisko, kwota AFS 56 %); `DRYRUN=1 s427_job.sh study Dplus`: OK.
- Katalog testowy `/tmp/klyczak_s427test` na lxplus993 (ok. 60 MB) został na miejscu.

## 7. Wyniki i wnioski

Do uzupełnienia po biegu, wg §3. Tabele: `WYNIKI.md`.

Stan 17.11.2026 14:10 — pełne runy, tabele w `WYNIKI.md`.

1. **Pytanie 1 (SALAD v1.0 w czystym minbias):** D⁰ — potwierdzony. Bramki K1 spełnione (n_eff/N 0,93, C2ST pasma 0,508). Pod pikiem R − 1 = +3,6 ± 4,7 % (WP80), +4,1 ± 6,6 % (WP70), +6,1 ± 11,1 % (WP50) w binie „all” i +3,2 / +3,4 / +8,0 % w 3–5 GeV: K2 ZGODNE, K3 spełnione. Biny 5–8 i 8–30 GeV nierozstrzygnięte (+7,8 ± 5,7 % i +8,6 ± 10,0 % przy WP80). Werdyktu PRZEKROCZONE nie ma w żadnej komórce. D⁺ — niemiarodajne pod pikiem: R − 1 = +1,7 % przy σ = 17 % (WP80), błąd w całości z interpolacji pasm bocznych. W paśmie walidacyjnym D⁺ SALAD domyka 3–5 GeV (1,186 → 1,008 ± 0,009 przy WP80), 5–8 GeV częściowo (1,045 ± 0,029), 8–30 GeV nie (0,873 ± 0,053).
2. **Pytanie 2 (transfer między klasami trygera):** D⁰, wagi z MB427929 na danych J12only: R = 1,052 / 1,018 / 0,993 ± 0,034 (WP80/70/50, „all”), ale C2ST 0,637 → „przenosi ε_B, nie rozkłady”; w binie 5–8 GeV nie przenosi (1,171 ± 0,060). D⁺: każdy transfer ma σ ≈ 16–20 %, nierozstrzygnięte; wagi z J12only na MB427929 dają 1,20 / 1,27 / 1,41. K8: AUC_SB J12only − MB427929 = 0,084 (D⁰) i 0,065 (D⁺) > 0,05, n_eff/N 0,57–0,61 → w próbce dżetowej SALAD pochłania różnicę trygerową.
3. **Pytanie 3 (czułość na skład, B1):** bez zmian względem 428770. `cat_local`: D⁰ 0,974 / 0,941 / 0,873 (B1: 0,97 / 0,95 / 0,85), D⁺ 0,857 / 0,797 / 0,727 (B1: 0,845 / 0,809 / 0,70); SALAD nie poprawia względem MC surowego, C2ST 0,50. Tryb `smooth` domyka (D⁺ 0,995 / 0,967 / 0,940 wobec 0,867 / 0,787 / 0,666 bez wag).

Hipoteza (a) potwierdzona dla D⁰; (b) częściowo: ε_B przenosi się w binie „all”, rozkłady nie.

## 8. Zastrzeżenia znane z góry

- Tagger H427 to HGB uczony na SB 427929, nie nominalny XGB `D0_H`; progi WP są jego własne. ε_B w liczbach bezwzględnych nie jest porównywalne z `BKG_427929`, stosunki R tak.
- Jeden tagger dla wszystkich klas: w klasach 428770 tło treningowe taggera pochodzi z innego runu. To zamierzone (stały tagger, zmienia się tylko model tła), ale różni się od B1, gdzie tagger uczono na SB tej samej próbki.
- J12only to 25 % zdarzeń klasy; „all 428770” nie jest liczone.
- Interpolacja exp(pol2) ma własne obciążenie (BKG_427929: 0,965–0,986 przy 3–5 GeV na MC); raportowane obok, nie odejmowane.
- Lewe pasmo D⁰ ma 40 MeV, prawe 60 MeV.
- Wariancja treningu z dwóch ziaren to oszacowanie zgrubne.
- MC to tylko mb_nd (bez sd/dd/cd), wagi zdarzeń 1.
- Adwersarialny B1 jest testem na MC; składu tła w danych 427929 nie mierzy.
- Kopia `charmbkg` w paczce zawiera stan roboczy pakietu (w tym pliki `afp_mix*`, nieużywane tutaj).

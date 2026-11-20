# SD_var: niepewność modelowa R = σ(D + proton w oknie AFP)/σ_incl(D) — wynik etapu `sdvar` (19.11.2026, projekt B)

**Poziom cząstek, Pythia 8.317, bez danych ATLAS.** Definicje jak `SD_pred` §2: D z antycząstkami, 2 ≤ pT < 30 GeV, |η| < 2,2; proton 0,030 ≤ ξ ≤ 0,100 po stronie A lub C (strony zsumowane); błędy „±” tylko stat. MC.
Kryterium ustalone przed wynikami (`PLAN.md` §4, `Q_lxplus/sd_var/scenarios.txt`) i wykonane bez zmian przez `analyze_var.py`. Pełne tabele: `out/tables.md`, `out/results.json`.

## 1. Wynik główny

| | R(D⁰) [%] | zakres [%] | uzysk D⁰ w 427929, pT 3–30 (scenariusz M z TAGREG §3) |
|---|---|---|---|
| dotąd (kształt SaS, zakres „na oko”) | 0,75 | 0,07 – 1,6 | 39 (7 – 152) |
| **teraz: nominał = oficjalne A3, obwiednia zbioru A** | **0,384 ± 0,005** | **0,154 – 1,53** (× 0,40 … × 3,97) | **18,6 (7,4 – 74)** |
| scenariusz twardy (typ 4), osobno | 1,061 ± 0,007 | 0,95 – 2,27 (× 0,89 … × 2,14) | 35 (31 – 75) |

- Centralna wartość spada dwukrotnie (0,75 → 0,38 %), bo oficjalne A3 ma strumień DL ε = 0,07, a nie kształt SaS. Zakres zwęża się z ×23 do ×10 (max/min), ale pozostaje asymetryczny: w górę ×4, w dół ×2,5.
- D⁺ i D*⁺: R = 0,389 ± 0,007 i 0,390 ± 0,009 %, te same czynniki (× 0,38 … × 3,80; × 0,39 … × 3,99). **R(D⁺)/R(D⁰) = 0,98–1,05 we wszystkich 44 scenariuszach** — stosunek gatunków jest niewrażliwy na model.
- Zależność od pT (nominał): 0,433 / 0,375 / 0,324 / 0,268 % w binach 2–3 / 3–5 / 5–8 / 8–30 GeV; spadek 2–3 → 8–30 o czynnik 1,5–1,9 w scenariuszach miękkich zbioru A, w twardym ×4,9.
- Przewidywany stosunek tagowane/inkluzywne po detektorze: P = R × A·ε_AFP × ε_ξcal = 0,384 % × 0,55 × 0,99 = **0,21 % (0,08 – 0,83 %)** wobec dotychczasowych 0,41 % (0,07 – 1,6 %); A·ε_AFP (0,6–0,9 w TAGREG) nie jest tu wariowane.
- Scenariusz S (kotwica z pasm bocznych, 35 (18–71) D⁰) nie zależy od generatora; leży w górnej połowie nowego zakresu M i jest zgodny ze scenariuszem twardym.

## 2. Wkłady grup (D⁰; czynniki względem odniesienia grupy)

| grupa | zakres | co wyznacza kraniec |
|---|---|---|
| M pełne modele (wobec A3off) | × 1,00 … × 3,17 | ABMST 1,217 %; SaS 0,750 %; A3 z ε = 0,037: 0,639 %; MBR 0,431 i 0,442 % — żaden model nie daje mniej niż A3 |
| P parametry kształtu (wobec DL ε = 0,07, σ SaS) | × 0,81 … × 1,54 | ε = 0,085 … ε = 0,04; α′ 0,10/0,40: × 0,84 / × 1,09; Berger–Streng × 0,91 |
| X układ X | × 0,47 … × 1,83 | σ_ref(Pom–p) = 25 mb … 5 mb; PDF Pomeronu × 0,88 – × 1,46; pT0Ref × 0,94 / × 1,06; CRqcd × 1,44 |
| L proton wiodący ND | × 0,92 … × 1,57 | ND Monash … ND CRqcd |
| normalizacja σ | × 0,63 … × 1,07 | σ_SD × 0,64 … σ_DD × 2 (CD w A3 = 0) |
| **razem** (większa z M, P; potem kwadratura w ln R) | **× 0,40 … × 3,97** | w górę dominuje M (ABMST), w dół X (σ_ref = 25 mb) i normalizacja σ_SD |

Obwiednia A + B (tylko informacyjnie): M × 0,96 … × 4,80 (ABMST bez przerwy), P × 0,31 … × 1,98 (strumienie H1; ε = 0,02), X × 0,47 … × 2,49 (PDF jak π⁰).

Obserwacje:
1. **Największa pojedyncza dźwignia to model σ i kształtu (ABMST ×3,2), druga — MPI w układzie Pomeron–proton** (`Diffraction:sigmaRefPomP` 5 → 25 mb zmienia R o czynnik 3,9). Tego parametru nie wiąże żaden pomiar.
2. **CRqcd podnosi składową ND z protonem wiodącym 3,4 razy** (1,38 → 4,76 µb D⁰ + p; udział nie-SD 17 → 36–39 %). Tło fizyczne „proton z resztek wiązki” jest więc niepewne co najmniej o czynnik 3.
3. Kontrole: σ_ND z `ref_A3` = 55,587 mb (zgodne z `nd_A3`); każdy wariant kształtu zmienia P(proton w oknie | SD); SaS mode 0 wobec mode 1 nadal 0,816 wobec 0,750 % (9 %, niewyjaśnione).
4. Udział zdarzeń z |t| < 1 GeV²: R_t1/R = 0,93 (nominał), 0,80 (z ND CRqcd) — okno t AFP zmniejszy głównie składową ND.

## 3. Twarda dyfrakcja (osobny mechanizm, nie wchodzi do sumy)

R(D⁰) = 1,061 % (DL-A3); zbiór A: pT0Ref 2,30 / 2,60 → 0,945 / 1,183 % (× 0,89 / × 1,12), bSelHard 2 → × 1,02, PDF GKG18 → × 1,08, strumień + PDF H1 Fit A → × 1,45, strumień SaS → × 1,41, H1 Fit B → × 2,14.
Zbiór B: pT0Ref = 3,00 → × 1,48, bSelHard 3 → × 1,13. Wrażliwość na pT0Ref (R7) jest tu umiarkowana (±11 % dla 2,30–2,60; +48 % dla 3,00); dominuje strumień.

## 4. Braki i odstępstwa (jawnie)

1. **`dv_CRqcd`: 9 z 10 części** — ziarno 1230008 kończy się naruszeniem ochrony pamięci Pythii (3 próby, deterministyczne). **`nv_CRqcd`: 16 z 20 części** — ziarna 1260008, 1260012, 1260013, 1260017 zawisły (limit 8 h; pozostałe 16 jobów trwało 1,0–2,0 h, więc to pętla w Pythii z rekoneksją QCD, nie tempo).
   Brakujących części nie dosyłano: statystyka wystarcza (błąd stat. R scenariuszy CRqcd i NDCRqcd = 1,5 %, próg kryterium 5 % → oba wchodzą do obwiedni), a priorytet konta jest dziś potrzebny etapom sieci.
   Utrata całych jobów według ziarna nie obciąża R, o ile zawieszenia są rzadkie (4 na 2·10⁷ zdarzeń); nie sprawdzono, czy zawieszające się zdarzenia mają szczególną topologię.
2. `merge` uruchomiony interaktywnie na lxplus (ten sam `sd_job.sh merge`, ~40 s), bo węzeł `merge` DAG-a był zablokowany przez nieudany `g_dv_CRqcd_8`. `DONE` na EOS → pętla ORCH sama pobrała wyniki; `state/sdvar` = `done` (36/36 modeli + `xsec.json`), bez ręcznej edycji stanu.
3. **`dv_A3alt` nie jest dosłownie oficjalnym A3 ALT.** Plik `Pythia8_A3_ALT_NNPDF23LO_Common.py` nie ustawia `SigmaDiffractive:mode = 0`, więc w Pythii 8.3 dawałby kształt SaS. Nasz wariant to A3 z ε = 0,037 i mode = 0; w zbiorze A pozostaje z reguły (ii) (parametr w zakresie pomiaru), nie (i). Na obwiednię nie wpływa (krańce wyznaczają ABMST i σ_ref).
4. **Podwójne liczenie CRqcd:** ten sam efekt (proton wiodący w ND) wchodzi przez grupę X (CRqcd × 1,44) i L (NDCRqcd × 1,57). Reguły nie zmieniono; dla informacji: bez grupy L górny czynnik wynosi × 3,68 zamiast × 3,97.
5. Kryteriów, zbiorów A/B i stałych (`STAT_MAX`, skale σ) nie zmieniano po obejrzeniu wyników.

## 5. Ustawienia A3 w repozytorium Athena — potwierdzone

Odczyt surowych plików z `gitlab.cern.ch/atlas/athena`, gałąź `main`, `Generators/Pythia8_i/share/common/` (19.11.2026, bez streszczania):
- `Pythia8_A3_NNPDF23LO_Common.py`: `SigmaDiffractive:PomFlux = 4`, `PomFluxEpsilon = 0.07`, `PomFluxAlphaPrime = 0.25`, `SigmaElastic:rho = 0.0`, **`SigmaDiffractive:mode = 0`**, `SigmaTotal:sigmaAXB =0`, `SigmaDiffractive:OwnbMinDD = 5`, `ColourReconnection:range = 1.8` (+ MPI jak w `base_A3.cmnd`). Zgodne z `models/a3off.cmnd`.
- `Pythia8_A3_ALT_NNPDF23LO_Common.py`: ε = 0,037, α′ = 0,25, **bez** `mode = 0`, `sigmaAXB`, `OwnbMinDD`, `rho`.
Niezweryfikowane: która wersja pliku (i która wersja Pythii) posłużyła do próbek mc23 używanych w analizie; liczby literaturowe oznaczone „z pamięci” w `PLAN.md` §3.

## 6. Zastrzeżenia

1. Obwiednia jest zdominowana przez dwa warianty (ABMST w górę, σ_ref(Pom–p) = 25 mb w dół); to rozrzut modeli, nie przedział ufności.
2. Skala σ_SD × 0,64 pochodzi z obszaru ALFA (ξ < 0,025); przeniesienie do okna 0,03–0,10 jest założeniem.
3. Brak Reggeonu, POMWIG/innego generatora, akceptancji AFP w t i φ. Proton wiodący w ND: tylko Monash i CRqcd.
4. R miękkie i twarde mają różne mianowniki (`SD_pred` §7.1) — porównywać R, nie σ.

## 7. Do decyzji / do zrobienia

1. **W kolejce nadal są 4 ponowienia `g_nv_CRqcd_{8,12,13,17}` (batch `charm_sdvar`, DAGMan 12789575)** — zawisną ponownie i zużyją do ~32–64 CPU·h. Trzeba je usunąć ręcznie: `condor_rm -constraint 'JobBatchName=="charm_sdvar"'`.
2. Trzy puste pliki `core.1255488`, `core.3545897`, `core.4069290` (0 B) w `~/private/charm_q/Q_lxplus/sd_var` na AFS oraz `lxplus932:/tmp/klyczak_sdvar` — nieusunięte.
3. Czy przyjąć R = 0,38 % (0,15–1,5 %) jako nowy scenariusz M w TAGREG (uzysk D⁰ 19 (7–74) zamiast 39 (7–152)).
4. Czy liczyć podwójnie CRqcd (pkt 4.4) i czy dosłać brakujące części CRqcd przy wolnym priorytecie (~7 CPU·h, nowe ziarna, z limitem czasu na job).

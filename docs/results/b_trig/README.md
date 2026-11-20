# B_trig — odporność wyników projektu B na trigger i niezależny run (17.11.2026)

Projekt B (tło podobne do sygnału; repo `~/Projects/charm/charm-bkg`, pakiet `charmbkg`). Nie mieszać z A (pary D–D).
Procedura z bramkami: `PROCEDURA.md`. Stan: kod gotowy i przetestowany na MC/toy/lokalnych h5 v1.2. Na danych v2.0 nic nie uruchomiono
(poza odczytowym licznikiem D* z ROOT, `lx/dshist_root.py` → `out/dshist_full.json`, `out/btrig_sens.md`). Nic nie wysłano.

| plik | rola |
|---|---|
| `btrig_common.py` | klasy triggera, nazwy kolumn (tp_/ps_, trigPassed_/trigPS_, EventInfo_*), złączenie (run, event), fit pass/fail, χ²/pulle, ε_B ML |
| `trg_root_extract.py` | ścieżka A: tabela zdarzeń z ROOT v2.0 (bity, PS, nSelBL, lead_pt, GRL, μ, liczniki D*) |
| `btrig_extract.py` | ścieżka B: h5 (v1.2/v2.0) → pary D0←D* z score H (jak EV) + ptcone + e_* (trigger/nSelBL/lead_pt na wierszu) + all_* z kluczem + tabela ev_* ; `--lxy-min 0.25` (zgodność v2.0 ↔ v1.2) |
| `btrig_emul.py` | emulacja triggera w MC: mapa ε_k(nSelBL, lead_pt) z 427929 + prescale’e 428770 (A) albo gęstość (B) |
| `btrig_analyze.py` | `effmap`, `sf`, `epsb`, `iso` (lokalnie) |
| `btrig_sens.py`, `btrig_toy.py` | uzysk D* per klasa/run i czułość; toy pass/fail przy statystyce klasy |
| `lx/` | paczka lxplus (wzór `Q_lxplus/dp_extract`): `bt_common.sh`, `bt_job.sh`, `bt.sub`, `make_dag_bt.sh`, `submit_bt.sh`, `bt_merge.py`, `dshist_root.py` |
| `orch/b_trig.sh` | etap ORCH — **niezainstalowany** (kopiować do `ORCH/stages/` + wpis w `pipeline.yaml` po zgodzie) |
| `tests/test_btrig.py` | 11 testów (konwencje nazw, klasy z −1, złączenie, χ², mapy emulacji, toy fit pass/fail) |

## Uruchomienie, gdy h5_v20 będzie na EOS (`/eos/user/k/klyczak/charm_h5_v20/parts/DONE`)
1. ORCH (zalecane): `cp charm_work/B_trig/orch/b_trig.sh charm_work/ORCH/stages/` i dopisać do `pipeline.yaml`
   (`- name: b_trig / needs: [h5_v20] / max_hours: 8`). Etap synchronizuje `B_trig/lx/` + 3 skrypty .py do `~/private/charm_q/Q_lxplus/b_trig`,
   `submit_bt.sh check` (schemat h5_v20, model D0_H, kolejka, AFS), submit DAG (`charm_btrig`), wynik → `charm_work/B_trig/runs/lxplus/`.
   Ścieżka A przed h5_v20 (tylko tabele zdarzeń z ROOT, ~30 min): `needs: [rucio_v20]` i `BT_PATH=A` w etapie.
2. Ręcznie na lxplus (sesja ssh): `cd ~/private/charm_q/Q_lxplus/b_trig && ./submit_bt.sh check && ./submit_bt.sh submit`.
3. Lokalnie (Python `.venv_h5`, `cd charm_work/B_trig`, `R=runs/lxplus`, `O=out`):
```
P=../../.venv_h5/bin/python
$P btrig_analyze.py effmap --ref-table $R/trg_data_427929.npz --closure-table $R/trg_data_428770.npz --out $O          # G5
$P btrig_analyze.py sf   --data $R/bt_data_428770.npz@428770 --data $R/bt_data_427929.npz@427929 --mc $R/bt_fd_nd.npz \
                         --emul-map $O/effmap_427929.json --nproc 8 --out $O                                         # G2, G6
$P btrig_analyze.py epsb --data $R/bt_data_428770.npz@428770 --data $R/bt_data_427929.npz@427929 --mc $R/bt_fd_nd.npz \
                         --mc-bkg $R/bt_mb_nd.npz --emul-map $O/effmap_427929.json --wps 0.80,0.70,0.50 --out $O         # G7
$P btrig_analyze.py iso  --data $R/bt_data_428770.npz@428770 --data $R/bt_data_427929.npz@427929 --mc $R/bt_fd_nd.npz \
                         --emul-map $O/effmap_427929.json --nproc 8 --out $O                                         # G8
```
   Szybka ścieżka A (bez h5_v20, tylko 428770): `sf --data ../EV_dstar/lxplus/ev_data.npz@428770 --event-table $R/trg_data_428770.npz
   --mc ../EV_dstar/lxplus/ev_fd.npz` (EV ma pary z pełnej lumi v1.1; złączenie po (run, event); bramka G1).
   Alternatywa bez ORCH: pobrać slim (`bin/lx -pull /eos/user/k/klyczak/charm_q/b_trig/ runs/lxplus/`).

## Testy (wykonane 17.11)
`../../.venv_h5/bin/python -m pytest -q tests` (11/11). Łańcuch na 2 lokalnych plikach h5 v1.2 danych + 8 fd z bitami toy:
`btrig_extract` (62 263 pary, 225 k D0 all), `sf` (36 fitów), `epsb` (G7 = „nie” we wszystkich binach dla losowych bitów — domknięcie),
`iso`, `bt_merge` (wykrywa duplikaty (run, event)), `make_dag_bt.sh` (atrapa EOS).

## Stan 17.11 ~02:30 (etap 2: dokończenie przed biegiem na danych)
- Sprawdzone ODCZYTEM na faktycznym h5_v20 (lxplus, /tmp, posprzątane): `bt_job.sh ext` na data_427929 _000003 (546 par, 24 s), fd_nd _000002 (9080 par, 14 s),
  data_428770 _000003 (12 953 pary, 18 s); wszystkie kolumny trigger/PS/nSelBL/lead_pt/GRL/ptcone znalezione; `merge` (dup = 0), `check` (schemat OK), DAG 70 węzłów
  (67 ext po 10 plików + 2 trg + merge); `trg_root_extract` 2 pliki ROOT 9 s.
- Zmiany paczki: `merge` czeka też na `trg_*` (DONE = komplet), listy w `lists/g<BT_GROUP>/`, ext `longlunch`, trg 8 GB, `BT_NODONE=1` (pominięcie testu parts/DONE).
- Etap zainstalowany: `ORCH/stages/b_trig.sh` (= `orch/b_trig.sh`, + `runs/lxplus/diag.txt`). Wpis do `pipeline.yaml` i start — do zrobienia.
- Format wyjścia: npz (nie parquet) — `btrig_analyze.py` czyta npz; score H liczony w jobie modelem `runs/lxplus_v11/models/D0_H` + `ev_dstar/src` (jak EV).
- **Analiza jedną komendą:** `bash charm_work/B_trig/run_analysis.sh` → `runs/lxplus/out/RAPORT_DANE.md` (+ `gates.json`, `btrig_*.md`); raport przez `btrig_report.py`,
  odświeżany po każdym kroku (effmap → epsb → sf → iso). Fit pass/fail jest wolny (~1–3 min CPU/fit): całość ~1,5–3 h na 8 rdzeniach; `FAST=1` = bez skanu profilu.
- Test komendy na 1 pliku/próbkę v2.0: łańcuch effmap/sf/epsb + raport przechodzi (iso przerwane dla czasu), fity nieokreślone (1/97 statystyki) — bramki G2/G6/G8
  sprawdzone tylko mechanicznie (syntetyczne wiersze). G5 na tabelach „tylko zdarzenia z kandydatem” nie domyka się (χ²/ndf ≫ 2) — ocenić na pełnych `trg_*.npz`; bramka miękka.

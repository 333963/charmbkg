#!/bin/bash
# B_trig — lokalna analiza po pobraniu wyników etapu ORCH b_trig (projekt B; pojedyncze D, nie ślepe).
#   bash charm_work/B_trig/run_analysis.sh [KATALOG_WEJŚCIA=runs/lxplus] [KATALOG_WYJŚCIA=<wejście>/out]
# Wejście: bt_{data_428770,data_427929,fd_nd,mb_nd}.npz (+ trg_data_{427929,428770}.npz, bt_summary.json). Wyjście: btrig_{effmap,sf,epsb,iso}.{md,json}
#   + RAPORT_DANE.md (bramki G1–G8 z PROCEDURA.md), odświeżany po każdym kroku (kolejność: effmap, epsb [szybkie], sf, iso).
#   Zmienne: NPROC (8), FAST=1 (bez skanu profilu → błędy HESSE), WPS (0.95,0.80,0.50; pełne: 0.95,0.90,0.80,0.70,0.50), EPSB_WPS, STEPS.
#   Czas: fit pass/fail (odporny, wielostartowy) ~1–3 min CPU/fit; sf = 9 klas × 3 biny × |WPS| fitów, iso = 180 fitów → łącznie ~1,5–3 h na 8 rdzeniach.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
P=${PYTHON:-$HERE/../../.venv_h5/bin/python}
R=${1:-$HERE/runs/lxplus}; O=${2:-$R/out}
NPROC=${NPROC:-8}; WPS=${WPS:-0.95,0.80,0.50}; EPSB_WPS=${EPSB_WPS:-0.80,0.70,0.50}
FAST=(); [ "${FAST_FIT:-${FAST:-0}}" = 1 ] && FAST=(--fast)
mkdir -p "$O"; LOG=$O/run_analysis.log; : > "$LOG"
say() { echo "[$(date +%T)] $*" | tee -a "$LOG"; }
DATA=(); for r in 428770 427929; do [ -s "$R/bt_data_$r.npz" ] && DATA+=(--data "$R/bt_data_$r.npz@$r") || say "UWAGA: brak $R/bt_data_$r.npz"; done
[ ${#DATA[@]} -gt 0 ] || { say "STOP: brak bt_data_*.npz w $R"; exit 1; }
[ -s "$R/bt_fd_nd.npz" ] || { say "STOP: brak $R/bt_fd_nd.npz"; exit 1; }
MC=(--mc "$R/bt_fd_nd.npz"); MB=(); [ -s "$R/bt_mb_nd.npz" ] && MB=(--mc-bkg "$R/bt_mb_nd.npz")
# tabele zdarzeń do mapy emulacji: trg_* (ROOT, wszystkie zdarzenia) albo ev_* z bt_* (tylko zdarzenia z kandydatem)
T9=$R/trg_data_427929.npz; [ -s "$T9" ] || T9=$R/bt_data_427929.npz
T8=$R/trg_data_428770.npz; [ -s "$T8" ] || T8=$R/bt_data_428770.npz
EMUL=(); rc=0
if [ -s "$T9" ]; then
  CL=(); [ -s "$T8" ] && CL=(--closure-table "$T8")
  say "effmap: ref $T9 closure $T8"
  if "$P" "$HERE/btrig_analyze.py" effmap --ref-table "$T9" "${CL[@]}" --out "$O" >>"$LOG" 2>&1 && [ -s "$O/effmap_427929.json" ]; then
    EMUL=(--emul-map "$O/effmap_427929.json"); else say "UWAGA: effmap nieudane → wagi MC = 1 (G5 = brak)"; fi
fi
for step in ${STEPS:-epsb sf iso}; do
  extra=(); wps=$WPS
  [ $step = epsb ] && { extra=("${MB[@]}"); wps=$EPSB_WPS; }
  say "$step …"
  "$P" "$HERE/btrig_analyze.py" $step "${DATA[@]}" "${MC[@]}" "${extra[@]}" "${EMUL[@]}" --wps "$wps" --nproc "$NPROC" "${FAST[@]}" --out "$O" >>"$LOG" 2>&1 \
    || { say "BŁĄD: $step (log: $LOG)"; rc=1; }
  "$P" "$HERE/btrig_report.py" --in "$R" --out "$O" >>"$LOG" 2>&1 || rc=1      # raport częściowy po każdym kroku
done
grep -A12 '^zapisano .*RAPORT_DANE' "$LOG" | tail -13
say "koniec rc=$rc → $O/RAPORT_DANE.md"
exit $rc

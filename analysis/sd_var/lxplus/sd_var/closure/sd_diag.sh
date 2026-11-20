#!/bin/bash
# sd_diag.sh diag <model> <seed> <nev> <cmnd1,...>  — diagnostyka awarii CRqcd (19.11.2026): ten sam bin/pythia_sd pod gdb.
#  Awaria (SIGSEGV) → ślad stosu. Brak końca po DIAG_DEADLINE s (domyślnie 10800) → SIGINT do Pythii → ślad stosu miejsca pętli.
#  Nic nie zapisuje na EOS; wynik tylko w stdout (log na AFS, mały). Scratch mktemp + trap.
set -eo pipefail
SD=$(cd "$(dirname "$0")/.." && pwd); export SD
source "$SD/sd_common.sh"
step=$1; model=$2; seed=$3; nev=$4; cm=$5
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/sddiag_${model}_${seed}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"
echo "[$(date -u +%FT%TZ)] sd_diag $model seed $seed nev $nev host $(hostname) deadline ${DIAG_DEADLINE:-10800}s"
q_env
files=(); for c in ${cm//,/ }; do files+=("$SD/models/$c.cmnd"); done
LCG_LD=$LD_LIBRARY_PATH
ulimit -c 0; echo "ulimit -s: $(ulimit -s)"
env -u LD_LIBRARY_PATH -u PYTHONPATH -u PYTHONHOME /usr/bin/gdb -q -batch \
  -ex "set environment LD_LIBRARY_PATH $LCG_LD" -ex "set pagination off" -ex "set print frame-arguments none" -ex run \
  -ex "echo \n### SLAD STOSU (pierwsze 40 ramek)\n" -ex "bt 40" -ex "echo \n### GLEBOKOSC STOSU (ostatnie 5 ramek)\n" -ex "bt -5" \
  -ex "echo \n### rejestry\n" -ex "info registers rip rsp" -ex "info proc stat" \
  --args "$SD/bin/pythia_sd" "$nev" "$seed" part.tsv "${files[@]}" > gdb.log 2>&1 &
G=$!
t0=$(date +%s); DL=${DIAG_DEADLINE:-10800}; sent=0
while kill -0 $G 2>/dev/null; do
  sleep 20
  if [ $sent = 0 ] && [ $(( $(date +%s) - t0 )) -ge "$DL" ]; then
    echo "[$(date -u +%FT%TZ)] termin $DL s — SIGINT do pythia_sd (CPU dotąd: $(ps -o cputime= -C pythia_sd | head -1))"
    pkill -INT -x pythia_sd -u "$(id -u)" || true; sent=1
  fi
  if [ $sent = 1 ] && [ $(( $(date +%s) - t0 )) -ge $((DL + 300)) ]; then kill -9 $G 2>/dev/null || true; pkill -9 -x pythia_sd -u "$(id -u)" || true; fi
done
wait $G || true
echo "[$(date -u +%FT%TZ)] koniec po $(( $(date +%s) - t0 )) s; sent=$sent; part.tsv: $(wc -l < part.tsv 2>/dev/null || echo brak) wierszy, ostatni: $(tail -n 1 part.tsv 2>/dev/null | cut -c1-40)"
echo "### gdb.log (bez tabeli ustawień)"; grep -vE '^ \|| \*---' gdb.log | tail -n 150
echo "[$(date -u +%FT%TZ)] OK sd_diag $model $seed"

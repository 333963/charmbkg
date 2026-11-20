#!/bin/bash
# smoke.sh — test dymny BEZ condora (interaktywnie, ≤ 10 min): każdy model × 1 job × mało zdarzeń + xsec + merge do katalogu lokalnego.
# Użycie: SD_OUT=/tmp/klyczak_sdvar/out ./smoke.sh [nev_gen=3000] [nev_xsec=5000] [równolegle=4]
set -uo pipefail
SD=$(cd "$(dirname "$0")" && pwd); export SD
: "${SD_OUT:?ustaw SD_OUT na katalog lokalny (nie EOS)}"; case $SD_OUT in /eos/*) echo "STOP: smoke nie pisze na EOS"; exit 1 ;; esac
export SD_OUT; NG=${1:-3000}; NX=${2:-5000}; P=${3:-4}
source "$SD/sd_common.sh"
sd_build || exit 1
mkdir -p "$SD_OUT/log"; t0=$(date +%s)
{ while read -r idx model nj nev cm; do
    [ "$nj" -gt 0 ] && echo "gen $model $((1000000 + idx * 10000 + 1)) $NG $cm"
    sd_needs_xsec "$model" && echo "xsec $model $((1000000 + idx * 10000)) $NX $cm"
  done < <(sd_models); } > "$SD_OUT/log/tasks.txt"
xargs -P "$P" -L 1 bash -c 's=$(date +%s); "$SD/sd_job.sh" "$0" "$1" "$2" "$3" "$4" > "$SD_OUT/log/$0_$1.log" 2>&1; echo "$0 $1 rc=$? $(( $(date +%s) - s )) s"' < "$SD_OUT/log/tasks.txt" | tee "$SD_OUT/log/summary.txt"
"$SD/sd_job.sh" merge > "$SD_OUT/log/merge.log" 2>&1; echo "merge rc=$?" | tee -a "$SD_OUT/log/summary.txt"
echo "razem $(( $(date +%s) - t0 )) s; błędy: $(grep -vc 'rc=0' "$SD_OUT/log/summary.txt")"

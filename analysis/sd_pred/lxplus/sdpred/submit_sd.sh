#!/bin/bash
# submit_sd.sh check|submit|status|fetch   (lxplus, ~/private/charm_q/Q_lxplus/sdpred)
#  check  — pliki paczki, cmnd wszystkich modeli, kompilacja bin/pythia_sd + test 200 zdarzeń, nic w kolejce (charm_sd), brak $SD_OUT/DONE (FORCE=1), AFS < 85 %
#  submit — check + make_dag_sd.sh → condor_submit_dag -f (batch charm_sd, -maxjobs $MAXJOBS) w CZYSTYM env; DRYRUN=1 → tylko DAG
# Zmienne: SD_OUT, SD_ONLY=model1,model2, SD_NEV, SD_NJOBS (test), FORCE=1, MAXJOBS=150
set -uo pipefail
SD=$(cd "$(dirname "$0")" && pwd); export SD
source "$SD/sd_common.sh"
cmd=${1:-check}
check() {
  local ok=0 f q n
  for f in pythia_sd.cc sd_reduce.py sd_merge.py sd_job.sh sd.sub models.txt make_dag_sd.sh; do [ -f "$SD/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  while read -r _ model nj nev cm; do for f in ${cm//,/ }; do [ -f "$SD/models/$f.cmnd" ] || { echo "STOP: brak models/$f.cmnd ($model)"; ok=1; }; done; done < <(sd_models)
  n=$(sd_models | awk '{j+=$3; e+=$3*$4} END{printf "%d modeli-jobów, %.3g zdarzeń", j, e}'); echo "  plan: $(sd_models | wc -l) modeli, $n"
  sd_build || { echo "STOP: kompilacja pythia_sd"; ok=1; }
  if [ $ok = 0 ]; then
    ( q_env >/dev/null 2>&1; t=$(mktemp -d); cd "$t" && "$SD/bin/pythia_sd" 200 1 t.tsv "$SD/models/base_A3.cmnd" "$SD/models/proc_hq.cmnd" "$SD/models/hard_t3.cmnd" > t.log 2>&1 \
        && python3 "$SD/sd_reduce.py" t.tsv t.tsv.json t.npz | sed 's/^/  test: /'; rc=$?; cd /; rm -rf "$t"; exit $rc ) || { echo "STOP: test pythia_sd/sd_reduce (doHard)"; ok=1; }
  fi
  condor_q -batch 2>/dev/null | grep -qE 'charm_sd|sd\.dag' && { echo "STOP: sdpred w kolejce"; ok=1; }
  [ -s "$SD_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: już policzone ($(cat "$SD_OUT/DONE")); FORCE=1 → nowy bieg (części istniejące są pomijane)"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  SD_OUT=$SD_OUT SD_ONLY=$SD_ONLY"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_sd.sh submit)" ;;
  submit) check || exit 1
          chmod +x "$SD"/*.sh "$SD"/*.py
          "$SD/make_dag_sd.sh" || exit 1
          if [ "${DRYRUN:-0}" = 1 ]; then grep -c '^JOB' "$SD/sd.dag"; grep -h '^VARS' "$SD/sd.dag" | head -3; exit 0; fi
          [ "${FORCE:-0}" = 1 ] && [ -f "$SD_OUT/DONE" ] && mv -f "$SD_OUT/DONE" "$SD_OUT/DONE.old_$(date +%s)"
          cd "$SD" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_sd -maxjobs "$MAXJOBS" sd.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_sd|Total for klyczak'
          for d in "$SD_OUT"/*/; do [ "$(basename "$d")" = merged ] || echo "  $(basename "$d"): $(ls "$d" 2>/dev/null | grep -c '\.npz$') npz"; done; ls "$SD_OUT/merged" 2>/dev/null | tr '\n' ' '; echo ;;
  fetch)  find "$SD_OUT/merged" -maxdepth 1 -type f 2>/dev/null; ls "$SD_OUT/DONE" 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac

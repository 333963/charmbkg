#!/bin/bash
# submit_emb.sh check|submit|status|fetch   (na lxplus, w ~/private/charm_q/Q_lxplus/sigemb)
#  check  — wejścia: ROOT v2.0 (dpt3to150) per próbka, emb_donors.npz, ../trkjet2/tj2_kernel.py; nic w kolejce (charm_sigemb); brak $EMB_OUT/DONE; AFS < 85 %
#  submit — check + make_dag_emb.sh → condor_submit_dag emb.dag (batch charm_sigemb) w CZYSTYM env; DRYRUN=1 → tylko DAG + komendy
set -uo pipefail
EMB=$(cd "$(dirname "$0")" && pwd); export EMB
source "$EMB/emb_common.sh"
cmd=${1:-check}
check() {
  local ok=0 s q r
  for f in emb_extract.py emb_merge.py emb_donors.npz; do [ -s "$EMB/$f" ] || { echo "STOP: brak $f"; ok=1; }; done
  [ -s "$Q/trkjet2/tj2_kernel.py" ] || { echo "STOP: brak $Q/trkjet2/tj2_kernel.py (jądro TRKJET2)"; ok=1; }
  for s in $EMB_SAMPLES; do
    r=$(ls "$RAWV20/$s" 2>/dev/null | grep -c 'dpt3to150.*\.root$')
    echo "  $s: ROOT $r, węzłów $(emb_np "$s"), prescale $(emb_ps "$s")"
    [ "$r" -gt 0 ] || { echo "STOP: $s — brak ROOT v2.0"; ok=1; }
  done
  condor_q -batch 2>/dev/null | grep -qE 'charm_sigemb|emb\.dag' && { echo "STOP: SIGEMB w kolejce"; ok=1; }
  [ -s "$EMB_OUT/DONE" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: SIGEMB już policzone ($(cat "$EMB_OUT/DONE")); FORCE=1 nadpisuje albo inny EMB_OUT"; ok=1; }
  q=$(fs lq ~/private 2>/dev/null | awk 'NR==2{gsub("%","",$4); print $4}')
  [ -n "$q" ] && [ "$q" -ge 85 ] && { echo "STOP: kwota AFS ${q}%"; ok=1; }
  echo "kwota AFS: ${q:-?}%  EMB_OUT=$EMB_OUT EMB_NPARTS=$EMB_NPARTS EMB_PRESCALE=$EMB_PRESCALE"; return $ok
}
case $cmd in
  check)  check && echo "OK — gotowe (./submit_emb.sh submit)" ;;
  submit) [ -z "${EMB_DAG:-}" ] || { echo "STOP: EMB_DAG tylko do testów make_dag_emb.sh"; exit 1; }
          check || exit 1
          "$EMB/make_dag_emb.sh" || exit 1
          chmod +x "$EMB"/*.sh "$EMB"/*.py
          if [ "${DRYRUN:-0}" = 1 ]; then for s in $EMB_SAMPLES; do DRYRUN=1 "$EMB/emb_job.sh" ext "$s" 0; done; DRYRUN=1 "$EMB/emb_job.sh" merge; exit 0; fi
          ( emb_env ) || exit 1   # podpowłoka: LCG psuje condor_dagman (GLIBCXX_3.4.31) — submit w czystym env
          cd "$EMB" && env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -f -batch-name charm_sigemb emb.dag ;;
  status) condor_q -batch 2>/dev/null | grep -E 'charm_sigemb|Total for klyczak'; ls -la "$EMB_OUT" 2>/dev/null | head -20; ls "$EMB_OUT/parts" 2>/dev/null | wc -l ;;
  fetch)  find "$EMB_OUT" -maxdepth 1 \( -name 'emb_*.npz' -o -name '*.json' -o -name DONE \) 2>/dev/null ;;
  *) echo "użycie: $0 check|submit|status|fetch"; exit 2 ;;
esac

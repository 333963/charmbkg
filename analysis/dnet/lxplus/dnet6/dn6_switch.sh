#!/bin/bash
# dn6_switch.sh cpu   (na lxplus, w ~/private/charm_q/Q_lxplus/dnet6) — przejście treningu z GPU na CPU po foldach BEZ podwójnych wyników:
#   1) joby GPU (gtrain) czekające w kolejce (JobStatus 1/5) są wycofywane (to samo zrobiłby ich periodic_remove po DN6_GPU_WAIT s);
#      job GPU, który JUŻ LICZY, zostaje (kończy swój tag w ~15 min) — chyba że FORCE=1 (wtedy też wycofany; jego gotowe foldy zostają na EOS i są użyte);
#   2) czeka (≤ DN6_SWITCH_WAIT s, domyślnie 1500), aż DAGMan charm_dnet6 zniknie z kolejki (pada, bo węzeł gtrain nie powiódł się);
#   3) DN6_TRAIN=cpu ./submit_dn6.sh submit — nowy DAG tylko dla tagów bez meta na EOS (fold × 5 + merge), potem ana (+ app).
# Brak DAG-a w kolejce → od razu krok 3. Niczego nie kasuje z dysku. Stan: ./submit_dn6.sh status
set -uo pipefail
cd "$(dirname "$0")" || exit 1
[ "${1:-}" = cpu ] || { echo "użycie: $0 cpu"; exit 2; }
source ./dn6_common.sh
[ -s "$DN6_OUT/DONE" ] && { echo "STOP: DONE już jest ($(cat "$DN6_OUT/DONE")) — nic do przełączania"; exit 0; }
C='JobBatchName == "charm_dnet6" && JobUniverse != 7 && RequestGPUs > 0'
if condor_q -batch 2>/dev/null | grep -qE 'charm_dnet6( |$)'; then
  nrun=$(condor_q -nobatch -constraint "$C && JobStatus == 2" -af ClusterId 2>/dev/null | wc -l)
  nidle=$(condor_q -nobatch -constraint "$C && JobStatus != 2" -af ClusterId 2>/dev/null | wc -l)
  ncpu=$(condor_q -nobatch -constraint 'JobBatchName == "charm_dnet6" && JobUniverse != 7 && RequestGPUs =!= 1' -af ClusterId 2>/dev/null | wc -l)
  echo "w kolejce: GPU liczące $nrun, GPU czekające $nidle, CPU $ncpu"
  if [ "$nrun" -eq 0 ] && [ "$nidle" -eq 0 ]; then echo "STOP: DAG charm_dnet6 żyje i nie ma jobów GPU (tryb CPU już biegnie albo ana/app) — nic nie robię"; exit 1; fi
  if [ "${FORCE:-0}" = 1 ]; then condor_rm -constraint "$C" || exit 1; else [ "$nidle" -gt 0 ] && { condor_rm -constraint "$C && JobStatus != 2" || exit 1; }; fi
  t0=$(date +%s)
  while condor_q -batch 2>/dev/null | grep -qE 'charm_dnet6( |$)'; do
    [ $(( $(date +%s) - t0 )) -gt "${DN6_SWITCH_WAIT:-1500}" ] && { echo "STOP: DAG nadal w kolejce po ${DN6_SWITCH_WAIT:-1500}s (liczy job GPU?) — ponów później albo FORCE=1"; exit 1; }
    sleep 20
  done
fi
DN6_TRAIN=cpu ./submit_dn6.sh submit

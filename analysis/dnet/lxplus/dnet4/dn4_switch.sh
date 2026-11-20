#!/bin/bash
# dn4_switch.sh gpu|cpu   (na lxplus, w ~/private/charm_q/Q_lxplus/dnet4) — przełącza NIEUKOŃCZONE węzły treningu sieci (train_d0s*, train_dps*, train_d0sb)
# biegnącego DAG-a dnet4 między GPU (dn4_gpu.sub) a CPU (dn4.sub): condor_rm DAGMana → rescue DAG (węzły DONE zachowane) → podmiana pliku .sub w dn4.dag
# (kopia dn4.dag.bak_<epoch>) → ponowny condor_submit_dag BEZ -f (używa rescue). Niczego nie kasuje. Odmawia, gdy któryś trening sieci właśnie biegnie
# (FORCE=1 pomija; wtedy sprawdź models/<tag> na EOS — tag bez meta/<tag>.json jest niegotowy i zostanie policzony od nowa).
set -uo pipefail
cd "$(dirname "$0")" || exit 1
mode=${1:?gpu|cpu}; case $mode in gpu) NEW=dn4_gpu.sub ;; cpu) NEW=dn4.sub ;; *) echo "użycie: $0 gpu|cpu"; exit 2 ;; esac
D=$PWD
run=$(condor_q -nobatch -constraint 'JobStatus == 2' -af Args 2>/dev/null | grep -cE 'dn4_job.sh train (d0s|dps|d0sb)')
[ "$run" -gt 0 ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: $run treningów sieci biegnie (FORCE=1 pomija)"; exit 1; }
id=$(condor_q -nobatch -constraint 'JobBatchName == "charm_dnet4" && JobUniverse == 7' -af ClusterId 2>/dev/null | head -1)
if [ -n "$id" ]; then
  n0=$(ls dn4.dag.rescue* 2>/dev/null | wc -l)
  condor_rm "$id" || exit 1
  for i in $(seq 1 60); do [ "$(ls dn4.dag.rescue* 2>/dev/null | wc -l)" -gt "$n0" ] && ! condor_q -batch 2>/dev/null | grep -q charm_dnet4 && break; sleep 2; done
  ls dn4.dag.rescue* >/dev/null 2>&1 || { echo "STOP: brak rescue DAG po condor_rm"; exit 1; }
  condor_q -batch 2>/dev/null | grep -q charm_dnet4 && { echo "STOP: joby charm_dnet4 nadal w kolejce"; exit 1; }
fi
cp -p dn4.dag "dn4.dag.bak_$(date +%s)" || exit 1
sed -i -E "s#^(JOB train_(d0s[0-9]+|dps[0-9]+|d0sb)) .*#\1 $D/$NEW#" dn4.dag
grep -E '^JOB train_' dn4.dag
env -u LD_LIBRARY_PATH -u PYTHONPATH PATH=/usr/bin:/bin:/usr/sbin:/sbin condor_submit_dag -batch-name charm_dnet4 dn4.dag

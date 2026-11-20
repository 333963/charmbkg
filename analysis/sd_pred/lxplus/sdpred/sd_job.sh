#!/bin/bash
# Węzeł sdpred. Użycie:
#   sd_job.sh gen <model> <seed> <nev> <cmnd1,cmnd2,...>  → pythia_sd + sd_reduce → $SD_OUT/<model>/part_<seed>.{tsv.gz,tsv.json,npz} (npz ostatni = znacznik; istnieje → pomija, FORCE=1 liczy)
#   sd_job.sh merge                                        → sd_merge.py → $SD_OUT/merged/<model>.json + summary.json + DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi.
set -eo pipefail
SD=$(cd "$(dirname "$0")" && pwd); export SD
source "$SD/sd_common.sh"
step=$1; model=${2:--}; seed=${3:-0}; nev=${4:-0}; cm=${5:--}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/sd_${step}_${model}_${seed}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"
echo "[$(date -u +%FT%TZ)] sdpred $step $model seed $seed nev $nev host $(hostname) SD_OUT=$SD_OUT"
q_env
case $step in
  gen)
    o=$SD_OUT/$model/part_$seed
    if [ "${FORCE:-0}" != 1 ] && eos_exists "$o.npz"; then echo "pomijam (istnieje) $o.npz"; exit 0; fi
    files=(); for c in ${cm//,/ }; do files+=("$SD/models/$c.cmnd"); done
    "$SD/bin/pythia_sd" "$nev" "$seed" part.tsv "${files[@]}" > gen.log 2>&1 || { tail -40 gen.log; exit 1; }
    grep -E 'PYTHIA (Error|Abort)' gen.log | sort | uniq -c | sort -rn | head -5 || true
    gzip -f part.tsv
    python3 -u "$SD/sd_reduce.py" part.tsv.gz part.tsv.json part.npz
    put_output part.tsv.gz "$o.tsv.gz"; put_output part.tsv.json "$o.tsv.json"; put_output part.npz "$o.npz" ;;
  merge)
    python3 -u "$SD/sd_merge.py" "$SD_OUT" "$WD/merged"
    for f in "$WD"/merged/*.json; do put_output "$f" "$SD_OUT/merged/$(basename "$f")"; done
    cp "$SD/models.txt" "$WD/models.txt"; put_output "$WD/models.txt" "$SD_OUT/merged/models.txt"
    date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$SD_OUT/DONE" ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK sdpred $step $model $seed"

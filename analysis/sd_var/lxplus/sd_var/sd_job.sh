#!/bin/bash
# Węzeł sd_var. Użycie:
#   sd_job.sh gen <model> <seed> <nev> <cmnd1,cmnd2,...>   → pythia_sd + sd_reduce → $SD_OUT/<model>/part_<seed>.{tsv.gz,tsv.json,npz} (npz = znacznik; istnieje → pomija, FORCE=1 liczy)
#   sd_job.sh xsec <model> <seed> <nev> <cmnd1,...>        → ta sama konfiguracja z proc_all zamiast proc_dif (SoftQCD:all) → $SD_OUT/xsec/<model>.json (σ per proces)
#   sd_job.sh merge                                         → sd_merge.py → $SD_OUT/merged/<model>.json + summary.json + xsec.json + DONE
# Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi.
set -eo pipefail
SD=$(cd "$(dirname "$0")" && pwd); export SD
source "$SD/sd_common.sh"
step=$1; model=${2:--}; seed=${3:-0}; nev=${4:-0}; cm=${5:--}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/sdvar_${step}_${model}_${seed}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"
echo "[$(date -u +%FT%TZ)] sd_var $step $model seed $seed nev $nev host $(hostname) SD_OUT=$SD_OUT"
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
  xsec)
    o=$SD_OUT/xsec/$model.json
    if [ "${FORCE:-0}" != 1 ] && eos_exists "$o"; then echo "pomijam (istnieje) $o"; exit 0; fi
    files=(); for c in ${cm//,/ }; do [ "$c" = proc_dif ] && c=proc_all; files+=("$SD/models/$c.cmnd"); done
    "$SD/bin/pythia_sd" "$nev" "$seed" x.tsv "${files[@]}" > gen.log 2>&1 || { tail -40 gen.log; exit 1; }
    put_output x.tsv.json "$o" ;;
  merge)
    python3 -u "$SD/sd_merge.py" "$SD_OUT" "$WD/merged" $(ls "$SD_OUT" | grep -E '^(dv|nv|hv)_' || true)
    python3 - "$SD_OUT/xsec" "$WD/merged/xsec.json" <<'PY'
import glob, json, os, sys
N = {"101": "ND", "103": "SD", "104": "SD", "105": "DD", "106": "CD"}; out = {}
for p in sorted(glob.glob(os.path.join(sys.argv[1], "*.json"))):
    j = json.load(open(p)); s = {}
    for k, v in j["proc"].items(): s[N.get(k.split("_")[0], k)] = s.get(N.get(k.split("_")[0], k), 0.0) + v["sigma_mb"]
    s["n_ok"] = j["n_ok"]; s["total"] = j["sigma_gen_mb_total"]; s["cmnd"] = ",".join(os.path.basename(c)[:-5] for c in j["cmnd"].split(","))
    out[os.path.basename(p)[:-5]] = s
json.dump(out, open(sys.argv[2], "w"), indent=1); print("xsec:", len(out), "modeli")
PY
    for f in "$WD"/merged/*.json; do put_output "$f" "$SD_OUT/merged/$(basename "$f")"; done
    for f in models.txt scenarios.txt; do cp "$SD/$f" "$WD/$f"; put_output "$WD/$f" "$SD_OUT/merged/$f"; done
    date -u +%FT%TZ > "$WD/DONE"; put_output "$WD/DONE" "$SD_OUT/DONE" ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK sd_var $step $model $seed"

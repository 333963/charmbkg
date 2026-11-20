#!/bin/bash
# Węzeł TAGTRIG. Użycie: tt_job.sh <data_RUN> <część> | tt_job.sh lumi | tt_job.sh finish
#  data_*: tt_bits.py (ROOT v2.0 kopiowany xrdcp pojedynczo do scratch) → $TT_OUT/parts/<próbka>.p<k>.{events,ps,lb}.parquet + .meta.json
#  lumi:   tt_lumi.sh w osobnym, czystym środowisku (Athena, COOL tylko odczyt); ZAWSZE kod 0 — pisze $TT_OUT/lumi_STATUS.txt
#  finish: tt_merge.py → pliki końcowe + summary.json + DONE
# Scratch mktemp + trap (kasowany); wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
TT=$(cd "$(dirname "$0")" && pwd); export TT
[ -n "${TT_SAMPLES:-}" ] && TT_SAMPLES=${TT_SAMPLES//,/ }   # z penv (przecinki) → spacje
[ -n "${TT_NFILES:-}" ] && TT_NFILES=${TT_NFILES//,/ }
source "$TT/tt_common.sh"
step=$1; part=${2:-0}
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/tt_${step}_${part}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] TAGTRIG $step/$part host $(hostname) WD=$WD TT_OUT=$TT_OUT RAW20=$RAW20"
case $step in
  lumi)
    # niekrytyczny: żadna porażka nie może dać kodu ≠ 0 (finish czeka na ten węzeł; dodatkowo SCRIPT POST /bin/true w DAG)
    set +e
    mkdir -p "$WD/lumi"
    run timeout "$TT_LUMI_TIMEOUT" env -i HOME="${HOME:-/tmp}" USER="${USER:-klyczak}" PATH=/usr/bin:/bin:/usr/sbin:/sbin TMPDIR="$WD" \
        TT_ASETUP="$TT_ASETUP" TT_LUMITAG="$TT_LUMITAG" TT_GRL_XML="$TT_GRL_XML" TT_SAMPLES="$TT_SAMPLES" \
        /bin/bash --noprofile --norc "$TT/tt_lumi.sh" "$WD/lumi"
    rc=$?
    if [ "${DRYRUN:-0}" != 1 ]; then
      [ -s "$WD/lumi/lumi_STATUS.txt" ] || echo "tt_lumi.sh nie zostawił statusu" > "$WD/lumi/lumi_STATUS.txt"
      echo "tt_lumi.sh kod $rc ($([ $rc = 124 ] && echo "TIMEOUT ${TT_LUMI_TIMEOUT}s" || echo koniec)) $(date -u +%FT%TZ)" >> "$WD/lumi/lumi_STATUS.txt"
      rm -f "$WD"/lumi/.asetup.save
      mkdir -p "$TT_OUT" 2>/dev/null
      for f in "$WD"/lumi/*; do [ -f "$f" ] && { tt_put "$f" "$TT_OUT/$(basename "$f")" || echo "UWAGA: nie zapisano $(basename "$f")"; }; done
      cat "$WD/lumi/lumi_STATUS.txt"
    fi
    echo "[$(date -u +%FT%TZ)] OK TAGTRIG lumi (niekrytyczny, kod 0)"
    exit 0 ;;
esac
if [ "${DRYRUN:-0}" != 1 ]; then tt_env; mkdir -p "$TT_OUT/parts" 2>/dev/null || true; fi
case $step in
  finish)
    run python3 -u "$TT/tt_merge.py" --out "$TT_OUT" --wd "$WD/fin" --samples "$TT_SAMPLES" --nparts "$TT_NPARTS" --nfiles "$TT_NFILES" --max-mb "$TT_MAX_MB"
    if [ "${DRYRUN:-0}" != 1 ]; then
      for f in "$WD"/fin/bits_*.parquet "$WD"/fin/*.csv "$WD"/fin/chains_*.json "$WD"/fin/summary.json; do tt_put "$f" "$TT_OUT/$(basename "$f")"; done
      date -u +%FT%TZ > "$WD/DONE"; tt_put "$WD/DONE" "$TT_OUT/DONE"
    fi ;;
  data_*)
    o=$step.p$part
    run python3 -u "$TT/tt_bits.py" --sample "$step" --rootdir "$RAW20/$step" --ptprod "$PTPROD" --part "$part" --nparts "$TT_NPARTS" \
        --out-prefix "$WD/$o" --wd "$WD" --step "$TT_STEP" ${TT_MAX_FILES:+--max-files "$TT_MAX_FILES"} ${TT_MAX_EVENTS:+--max-events "$TT_MAX_EVENTS"}
    if [ "${DRYRUN:-0}" != 1 ]; then
      for x in events.parquet ps.parquet lb.parquet meta.json; do tt_put "$WD/$o.$x" "$TT_OUT/parts/$o.$x"; done
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK TAGTRIG $step/$part"

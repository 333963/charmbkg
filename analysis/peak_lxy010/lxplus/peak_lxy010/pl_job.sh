#!/bin/bash
# Węzeł PEAK_LXY010. Użycie: pl_job.sh <próbka|finish>
# próbka: python3 -m charmtag.hyp_extract (h5 v2.0 → D0 z tagiem D*, obie hipotezy, cechy H, score D0_H z modeli foldów) → $PL_OUT/<próbka>.parquet
# finish: extract_summary.json (metadane + kontrola L_xy: min i udział < 0,25 mm per próbka), thresholds.json (progi WP nominalne), DONE.
# Scratch mktemp + trap (kasowany); wyjścia tylko EOS; AFS tylko logi. DRYRUN=1 → tylko komendy.
set -eo pipefail
PL=$(cd "$(dirname "$0")" && pwd); export PL
source "$PL/pl_common.sh"
step=$1
set -u
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/pl_${step}_XXXXXX")
cleanup() { cd /; rm -rf "$WD"; }
trap cleanup EXIT
cd "$WD"; export MPLCONFIGDIR=$WD/mpl TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
echo "[$(date -u +%FT%TZ)] PEAK_LXY010 $step host $(hostname) WD=$WD PL_OUT=$PL_OUT H5V20=$H5V20 NOMINAL=$NOMINAL"
if [ "${DRYRUN:-0}" != 1 ]; then pl_env; mkdir -p "$PL_OUT"; fi
case $step in
  finish)
    run python3 - "$PL_OUT" "$WD" <<'PY'
import glob, json, os, sys
import pyarrow.parquet as pq, pyarrow.compute as pc
out, wd = sys.argv[1], sys.argv[2]
S, thr = {}, None
for m in sorted(glob.glob(os.path.join(out, "*.meta.json"))):
    j = json.load(open(m))
    thr = thr or j.get("thresholds")
    s = j["sample"]
    S[s] = {k: j.get(k) for k in ("n_rows_out", "bytes", "frac", "score_mode", "hypo_swap_counts", "v12", "rows_total")}
    S[s]["dstar"] = (j.get("dstar") or {}).get("counts")
    t = pq.read_table(os.path.join(out, s + ".parquet"), columns=["Lxy", "pt", "score"])
    lx = t["Lxy"]
    S[s]["Lxy_min"] = pc.min(lx).as_py(); S[s]["frac_Lxy_lt_025"] = pc.mean(pc.cast(pc.less(lx, 0.25), "float64")).as_py()
    S[s]["pt_min"] = pc.min(t["pt"]).as_py(); S[s]["n_score_nan"] = int(pc.sum(pc.cast(pc.is_nan(t["score"]), "int64")).as_py() or 0)
need = os.environ["SAMPLES"].split()
miss = [s for s in need if s not in S]
if miss:
    sys.exit(f"brak wyników: {miss}")
bad = [s for s in need if not (S[s]["Lxy_min"] is not None and S[s]["Lxy_min"] < 0.2 and S[s]["frac_Lxy_lt_025"] > 0.05)]
if bad:
    sys.exit(f"STOP: próbki bez kandydatów L_xy < 0,25 mm (wejście nie jest h5 v2.0?): {bad} {[(s, S[s]['Lxy_min']) for s in bad]}")
nan = [s for s in need if S[s]["n_score_nan"]]
if nan:
    sys.exit(f"STOP: score NaN (brak modeli D0_H?): {nan}")
tot = sum(v["bytes"] or 0 for v in S.values())
json.dump({"samples": S, "total_bytes": tot, "input": os.environ.get("H5V20"), "score_mode": os.environ.get("SCORE_MODE"),
           "extra": os.environ.get("PL_EXTRA", "")}, open(os.path.join(wd, "extract_summary.json"), "w"), indent=1)
if thr:
    json.dump({"thresholds": thr, "source": "models/D0_H/meta.json"}, open(os.path.join(wd, "thresholds.json"), "w"), indent=1)
print(f"razem {tot/1e9:.2f} GB", {k: (v["n_rows_out"], v["score_mode"], round(v["Lxy_min"], 3), round(v["frac_Lxy_lt_025"], 3)) for k, v in S.items()})
PY
    if [ "${DRYRUN:-0}" != 1 ]; then
      pl_put "$WD/extract_summary.json" "$PL_OUT/extract_summary.json"
      [ -f "$WD/thresholds.json" ] && pl_put "$WD/thresholds.json" "$PL_OUT/thresholds.json"
      date -u +%FT%TZ > "$WD/DONE"; pl_put "$WD/DONE" "$PL_OUT/DONE"
    fi ;;
  fd_nd|mb_nd|data_428770|data_427929)
    set -f   # glob `*_none.h5` przekazywany dosłownie (rozwija go Python: sortowanie jak config.sample_paths)
    # shellcheck disable=SC2046
    run python3 -u -m charmtag.hyp_extract $(pl_args "$step") --out "$WD/$step.parquet"
    set +f
    if [ "${DRYRUN:-0}" != 1 ]; then
      pl_put "$WD/$step.parquet" "$PL_OUT/$step.parquet"
      pl_put "$WD/$step.meta.json" "$PL_OUT/$step.meta.json"
    fi ;;
  *) echo "nieznany krok $step"; exit 2 ;;
esac
echo "[$(date -u +%FT%TZ)] OK PEAK_LXY010 $step"

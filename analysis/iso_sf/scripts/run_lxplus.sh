#!/bin/bash
# ISO-FULL (16.11): analiza SF izolacji na pełnej próbce z lxplus (etap ORCH iso_full → ISO_sf/runs/lxplus).
# Nie rusza ISO_sf/out (wyniki 2-plikowe). Wyniki: runs/lxplus/out/{sf_fits.json,sf_table.md,sf_summary.json}, logi tamże.
# Użycie: scripts/run_lxplus.sh [NPROC]      zmienne: RUN_DIR (domyślnie ISO_sf/runs/lxplus), FORCE=1 (nadpisz runs/.../out)
set -euo pipefail
B=$(cd "$(dirname "$0")/.." && pwd); PY=${PY:-$HOME/Documents/charm/.venv_h5/bin/python}
RUN=${RUN_DIR:-$B/runs/lxplus}; OUT=$RUN/out
[ -s "$RUN/DONE" ] || { echo "STOP: brak $RUN/DONE (etap ORCH iso_full niepobrany)"; exit 1; }
for f in data_428770_tag.parquet fd_nd_tag.parquet; do [ -s "$RUN/$f" ] || { echo "STOP: brak $RUN/$f"; exit 1; }; done
[ -s "$OUT/sf_fits.json" ] && [ "${FORCE:-0}" != 1 ] && { echo "STOP: $OUT/sf_fits.json istnieje (FORCE=1 nadpisze)"; exit 1; }
mkdir -p "$OUT"
"$PY" - "$RUN" <<'PY'
import json, sys, pyarrow.parquet as pq
r = sys.argv[1]
for s in ("data_428770", "fd_nd", "mb_nd"):
    try:
        m = json.load(open(f"{r}/{s}_tag.meta.json")); n = pq.ParquetFile(f"{r}/{s}_tag.parquet").metadata.num_rows
        t4 = m.get("n_rows_t4"); print(f"{s}: {n} wierszy, {m.get('n_files', len(m.get('files', [])))} plików, T4 {t4}" + ("" if t4 in (None, n) else "  ⚠ NIEZGODNE"))
    except FileNotFoundError:
        print(f"{s}: brak")
PY
cd "$B/scripts"
ISO_IN=$RUN ISO_OUT=$OUT "$PY" sf_fit.py "${1:-6}" > "$OUT/sf_fit.log" 2>&1
ISO_OUT=$OUT "$PY" summarize.py > "$OUT/summarize.log" 2>&1
echo "OK → $OUT/sf_table.md"; head -25 "$OUT/sf_table.md"

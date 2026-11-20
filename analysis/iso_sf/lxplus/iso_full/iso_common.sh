# iso_common.sh — paczka ISO-FULL (16.11; izolacja D0 z tagu D*, wszystkie pliki h5 v1.2; source, nie uruchamiać). Wymaga ../common.sh.
IX=${IX:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$IX/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                    # EOS_OUT, q_env, put_output
: "${H5V12:=/eos/user/k/klyczak/charm_h5_v12/parts}"
: "${TAG_DIR:=$EOS_OUT/charm_q/hyp_d0_dstar}"            # T4 (hyp_dstar): <próbka>.parquet + .meta.json (tylko odczyt)
: "${ISO_OUT:=$EOS_OUT/charm_q/iso_full}"                # wyjścia (EOS): <próbka>_tag.parquet + .meta.json, iso_summary.json, DONE; części w parts/
: "${ISO_SAMPLES:=data_428770 fd_nd mb_nd}"
: "${ISO_NPARTS_DATA:=3}"                                # węzły na dane (39 plików → po 13); MC po 1
: "${ISO_DISK:=3GB}"
export IX Q H5V12 TAG_DIR ISO_OUT ISO_SAMPLES ISO_NPARTS_DATA ISO_DISK
ix_nparts() { case $1 in data_*) echo "$ISO_NPARTS_DATA" ;; *) echo 1 ;; esac; }
ix_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else cp -f "$1" "$2" && rm -f "$1"; fi; }
ix_env() {
  q_env
  export PYTHONDONTWRITEBYTECODE=1
  python3 -c 'import numpy, pandas, h5py, pyarrow.parquet' || { echo "BŁĄD: brak numpy/pandas/h5py/pyarrow w q_env"; return 1; }
}

# tt_common.sh — paczka TAGTRIG (18.11.2026, projekt B; uzupełnienie części A raportu TAGREG_2026-11-18). source, nie uruchamiać. Wymaga ../common.sh.
#   bits:   ROOT v2.0 (CharmAnalysis) → bity WSZYSTKICH łańcuchów w regionie tagowanym AFP + preskale/liczniki per LB × łańcuch (zdarzenia GRL)
#   lumi:   best-effort, NIEKRYTYCZNY: klucze i preskale L1/HLT per LB z COOL/TriggerDB (tylko odczyt) w środowisku Athena; świetlność tylko gdy podano tag
#   finish: scalenie → bits_<run>.parquet, ps_perLB_<run>.csv, lb_summary_<run>.csv, chains_<run>.json, summary.json, DONE
# W regionie tagowanym tylko liczniki i bity — żadnych mas kandydatów, histogramów ani fitów.
TT=${TT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$TT/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"                                    # EOS_OUT, q_env, put_output
: "${RAW20:=$EOS_OUT/charm_raw_v20}"                     # wejście (tylko odczyt): <próbka>/**/*dpt3to150*.root (etap rucio_v20)
: "${PTPROD:=dpt3to150}"                                 # NIE mieszać z dpt2to3
: "${TT_OUT:=$EOS_OUT/charm_q/tagtrig}"                  # wyjście (EOS); części w parts/
: "${TT_SAMPLES:=data_428770 data_427929}"
: "${TT_NPARTS:=4}"                                      # węzłów bits na próbkę
: "${TT_NFILES:=data_428770:97 data_427929:458}"         # oczekiwane pliki ROOT (h5_v20/samples_v20.txt, 17.11)
: "${TT_DISK:=8GB}"
: "${TT_STEP:=100MB}"                                    # porcja uproot (bez spacji — trafia do environment w DAG)
: "${TT_MAX_MB:=300}"                                    # limit łącznego rozmiaru wyjść końcowych
: "${TT_LUMI:=1}"                                        # 0 = bez węzła lumi
: "${TT_ASETUP:=Athena,25.0,latest;Athena,24.0,latest}"  # kandydaci asetup (rozdzielone ';'), pierwszy z AtlCoolTrigger.py wygrywa
: "${TT_LUMITAG:=}"                                      # tag świetlności (OflLumi…) — NIEZNANY lokalnie, nie zgadywać; pusty = bez iLumiCalc
: "${TT_GRL_XML:=}"                                      # ścieżka GRL_MinBiasAFP_lowmu2022.xml na lxplus — brak lokalnie; pusty = iLumiCalc -r <run> (gdy jest tag)
: "${TT_LUMI_TIMEOUT:=2400}"                             # s; po tym czasie węzeł lumi kończy z kodem 0 i statusem
: "${TT_CHECK_EVENTS:=20000}"                            # suchy przebieg w check: zdarzeń z 1. pliku każdej próbki (0 = pomiń)
export TT Q RAW20 PTPROD TT_OUT TT_SAMPLES TT_NPARTS TT_NFILES TT_DISK TT_STEP TT_MAX_MB TT_LUMI TT_ASETUP TT_LUMITAG TT_GRL_XML TT_LUMI_TIMEOUT TT_CHECK_EVENTS
tt_nexp() { local kv; for kv in $TT_NFILES; do [ "${kv%%:*}" = "$1" ] && { echo "${kv##*:}"; return; }; done; echo 0; }
tt_nroot() { find "$RAW20/$1" -name "*${PTPROD}*.root" 2>/dev/null | wc -l; }
tt_put() { if declare -F put_output >/dev/null; then put_output "$1" "$2"; else mkdir -p "$(dirname "$2")" && cp -f "$1" "$2"; fi; }
tt_env() {   # LCG + venv (q_env): uproot/awkward/pyarrow/pandas. NIE wołać przed condor_submit_dag ani w węźle lumi.
  q_env
  export PYTHONDONTWRITEBYTECODE=1
  python3 -c 'import numpy, awkward, uproot, pandas, pyarrow.parquet' || { echo "BŁĄD: brak numpy/awkward/uproot/pandas/pyarrow w q_env"; return 1; }
}

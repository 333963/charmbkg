#!/bin/bash
# TAGTRIG / węzeł lumi — best-effort, NIEKRYTYCZNY, COOL/TriggerDB TYLKO ODCZYT. Użycie: tt_lumi.sh <katalog_wyjścia>
# Uruchamiany przez tt_job.sh w czystym środowisku (env -i; bez LCG). Zawsze kończy kodem 0 i zostawia <wyjście>/lumi_STATUS.txt.
#  1. setupATLAS + asetup (pierwszy kandydat z TT_ASETUP, w którym jest AtlCoolTrigger.py)
#  2. per run: AtlCoolTrigger.py -r <run> -v → cool_keys_<run>.txt (klucze SMK / L1PSK / HLTPSK per zakres LB)
#              tt_cool.py → coolps_<run>.csv (preskale L1 i HLT osobno per zakres LB, wybrane łańcuchy AFP/J12/MBTS) + cool_raw_<run>.json
#  3. świetlność per LB (iLumiCalc) TYLKO gdy podano TT_LUMITAG (tag nieznany w materiałach lokalnych — nie zgadywany);
#     GRL: TT_GRL_XML (xml z repo CharmAnalysis data/GRL_MinBiasAFP_lowmu2022.xml; brak lokalnie) albo cały run (-r).
# NIETESTOWANE lokalnie (brak Atheny) — składnia narzędzi z pamięci; każdy krok ma własny wpis w statusie.
OUT=${1:?katalog wyjścia}; mkdir -p "$OUT"; cd "$OUT" || exit 0
TT=$(cd "$(dirname "$0")" && pwd)
ST=$OUT/lumi_STATUS.txt
say() { echo "$*" | tee -a "$ST"; }
RUNS=$(for s in ${TT_SAMPLES:-data_428770 data_427929}; do echo "${s#data_}"; done | tr ',\n' '  ')
HLT="HLT_noalg_L1AFP_A_OR_C HLT_mb_sptrk_pt4_L1AFP_A_OR_C HLT_noalg_L1J12 HLT_mb_mbts_L1MBTS_1 HLT_mb_sptrk_L1MBTS_1 HLT_mb_sptrk_pt2_L1MBTS_2 HLT_mb_sptrk_pt4_L1MBTS_2 HLT_j20_L1MBTS_2 HLT_noalg_mb_L1MBTS_1 HLT_noalg_mb_L1MBTS_2 HLT_noalg_L1RD0_FILLED"
L1="L1_AFP_A_OR_C L1_J12 L1_MBTS_1 L1_MBTS_2 L1_RD0_FILLED"
say "TAGTRIG lumi $(date -u +%FT%TZ) host $(hostname) runy: $RUNS"
say "ZAKRES: tylko odczyt COOL/TriggerDB. Świetlność: $([ -n "${TT_LUMITAG:-}" ] && echo "iLumiCalc z tagiem ${TT_LUMITAG}" || echo 'NIE liczona — brak tagu (TT_LUMITAG)')"
[ -n "${TT_LUMITAG:-}" ] || say "BRAK: tag świetlności 2022 low-mu (OflLumi…; strona LuminosityForPhysics) — nie ma go w materiałach lokalnych, nie zgadywano"
[ -n "${TT_GRL_XML:-}" ] || say "BRAK: GRL xml (CharmAnalysis data/GRL_MinBiasAFP_lowmu2022.xml) — nie ma go lokalnie; znane zakresy: 428770 LB 4–25 i 210–324, 427929 LB 119–324 z przerwami"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
if [ ! -f "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" ]; then
  say "BRAK: $ATLAS_LOCAL_ROOT_BASE (cvmfs atlas niedostępny na węźle) — nic nie odczytano"; say "WYNIK: brak danych COOL"; exit 0
fi
# shellcheck disable=SC1091
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet >/dev/null 2>&1
if [ -z "${AtlasSetup:-}" ] || [ ! -f "$AtlasSetup/scripts/asetup.sh" ]; then
  say "BRAK: AtlasSetup po atlasLocalSetup.sh (AtlasSetup='${AtlasSetup:-}') — nic nie odczytano"; say "WYNIK: brak danych COOL"; exit 0
fi
REL=""
IFS=';' read -ra CAND <<< "${TT_ASETUP:-Athena,25.0,latest}"
for r in "${CAND[@]}"; do
  # shellcheck disable=SC1091
  if ( source "$AtlasSetup/scripts/asetup.sh" "$r" >/dev/null 2>&1 && command -v AtlCoolTrigger.py >/dev/null 2>&1 ); then REL=$r; break; fi
  say "asetup $r: brak wydania albo AtlCoolTrigger.py"
done
if [ -z "$REL" ]; then say "BRAK: żaden kandydat asetup (${TT_ASETUP:-}) nie dał AtlCoolTrigger.py"; say "WYNIK: brak danych COOL"; exit 0; fi
# shellcheck disable=SC1091
source "$AtlasSetup/scripts/asetup.sh" "$REL" > "$OUT/lumi_asetup.log" 2>&1
say "asetup $REL: OK ($(command -v AtlCoolTrigger.py))"
AtlCoolTrigger.py -h > "$OUT/cool_help.txt" 2>&1
got=0
for run in $RUNS; do
  timeout 600 AtlCoolTrigger.py -r "$run" -v > "$OUT/cool_keys_$run.txt" 2> "$OUT/cool_keys_$run.err"; rc=$?
  say "AtlCoolTrigger.py -r $run -v: kod $rc, $(wc -l < "$OUT/cool_keys_$run.txt") linii → cool_keys_$run.txt"
  [ $rc = 0 ] && [ -s "$OUT/cool_keys_$run.txt" ] && got=1
  [ -s "$OUT/cool_keys_$run.err" ] || rm -f "$OUT/cool_keys_$run.err"
  # shellcheck disable=SC2086
  timeout 900 python "$TT/tt_cool.py" --run "$run" --out "$OUT/coolps_$run.csv" --raw "$OUT/cool_raw_$run.json" --hlt $HLT --l1 $L1 > "$OUT/cool_py_$run.log" 2>&1; rc=$?
  grep '^COOL:' "$OUT/cool_py_$run.log" | while read -r l; do say "  $l"; done
  say "tt_cool.py run $run: kod $rc, $([ -s "$OUT/coolps_$run.csv" ] && echo "$(($(wc -l < "$OUT/coolps_$run.csv") - 1)) wierszy → coolps_$run.csv" || echo 'bez coolps (patrz cool_py log)')"
  [ -s "$OUT/coolps_$run.csv" ] && got=1
done
if [ -n "${TT_LUMITAG:-}" ]; then
  if command -v iLumiCalc >/dev/null 2>&1; then
    trg=""; for c in $HLT; do trg+=" --trigger=$c"; done
    for run in $RUNS; do
      mkdir -p "$OUT/ilc_$run"
      if [ -n "${TT_GRL_XML:-}" ] && [ -f "$TT_GRL_XML" ]; then src="--xml=$TT_GRL_XML"; else src=""; say "iLumiCalc run $run: bez GRL xml → cały run (NIE świetlność GRL)"; fi
      # shellcheck disable=SC2086
      ( cd "$OUT/ilc_$run" && timeout 900 iLumiCalc --lumitag="$TT_LUMITAG" --livetrigger=L1_MBTS_1 $trg $src -r "$run" --plots ) > "$OUT/ilumicalc_$run.log" 2>&1; rc=$?
      say "iLumiCalc run $run (lumitag $TT_LUMITAG): kod $rc → ilumicalc_$run.log"
      for f in "$OUT/ilc_$run"/*.root; do [ -f "$f" ] && [ "$(stat -c %s "$f")" -lt 20000000 ] && mv "$f" "$OUT/ilumicalc_${run}_$(basename "$f")"; done
      rm -rf "$OUT/ilc_$run"
    done
  else
    say "BRAK: iLumiCalc w $REL"
  fi
fi
[ $got = 1 ] && say "WYNIK: są dane COOL (klucze/preskale) — sprawdzić treść; świetlność per LB: $([ -n "${TT_LUMITAG:-}" ] && echo 'patrz ilumicalc_*.log' || echo NIE MA)" || say "WYNIK: brak danych COOL"
exit 0

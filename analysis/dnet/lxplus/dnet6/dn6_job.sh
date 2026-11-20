#!/bin/bash
# Węzeł DAG DNET etap 6. Użycie:
#   dn6_job.sh gtrain <blok>           GPU: 5 foldów po kolei + łączenie w jednym jobie → $DN6_OUT/{models,folds,sc,meta}/<tag>
#   dn6_job.sh fold <blok> <k>         CPU: trening JEDNEGO foldu → models/<tag>/fold<k>.pt + folds/<tag>/fold<k>.json (json na końcu = „gotowe”)
#   dn6_job.sh merge <blok>            modele foldów z EOS → score OOF wszystkich tabel etapu 1 → sc/<tag>.npz, models/<tag>/meta.json, meta/<tag>.json
#   dn6_job.sh ana                     dn6_ana.py (P1 ze score) → dnet6_p1.json; dn6_study.py (maska: liczniki) → dnet6_mask.json;
#                                      analyze.py compute --sf-all (G na pasmach, rzeźbienie, SF z tagu D*) → dnet6_results.json; DONE
#   dn6_job.sh app <próbka> <k> <n>    (miękki, po ana; bramka: dnet6_p1.json p1_gate_seed1 = true) h5 → cechy → tory → dn6_apply (drzewo H + sieć z maską,
#                                      ziarno 1; 100 % zdarzeń) → $DN6_OUT/app/d150/app_<próbka>_<k>_<NNNNNN>.npz
#   dn6_job.sh anar                    (miękki) dn4_ana.py r --species D0 --force-data → dnet6_R.json (G_fit z fitu piku; R opisowo); DONE_R
# Każdy krok pomija to, co już jest na EOS (tag z meta = gotowy; DN_REDO=1 liczy od nowa) → GPU i CPU nie dublują wyników tego samego tagu.
# DN6_SOFTNODE=1: niepowodzenie → znacznik $DN6_OUT/fail/<węzeł> i kod 0. Scratch mktemp + trap; wyjścia tylko EOS; AFS tylko logi.
set -eo pipefail
DN6=$(cd "$(dirname "$0")" && pwd); export DN6
for v in DN_SAMPLES DN_GROUPS DN_BLOCKS DN2_NPARTS DN2_FRAC DN2_BLOCKS DN4_BLOCKS DN4_SOFT DN4_EXTDP DN4_APP150 DN4_APP23 DN5_BLOCKS DN5_GPU_BLOCKS DN6_BLOCKS DN6_APP150; do
  [ -n "${!v:-}" ] && printf -v "$v" '%s' "${!v//+/ }"
done
[ -n "${DN6_KEEP_B64:-}" ] && DN6_KEEP=$(echo "$DN6_KEEP_B64" | base64 -d) && export DN6_KEEP
[ -n "${DN6_APP_TAG_B64:-}" ] && DN6_APP_TAG=$(echo "$DN6_APP_TAG_B64" | base64 -d) && export DN6_APP_TAG
source "$DN6/dn6_common.sh"
step=$1; shift
set -u
NODE="${step}_$(echo "$*" | tr ' ' '_')"; NODE=${NODE%_}
T_START=$(date +%s)
WD=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/dn6_${step}_XXXXXX")
cleanup() { cd /; case "$WD" in */dn6_"${step}"_??????) rm -rf "${WD:?}" ;; esac; }   # tylko własny katalog mktemp
trap cleanup EXIT
cd "$WD"; export TMPDIR=$WD
run() { echo "+ $*"; [ "${DRYRUN:-0}" = 1 ] || "$@"; }
put() { [ "${DRYRUN:-0}" = 1 ] || dn_put "$1" "$2"; }
echo "[$(date -u +%FT%TZ)] DNET6 $step $* host $(hostname) WD=$WD DN_OUT=$DN_OUT DN2_OUT=$DN2_OUT DN6_OUT=$DN6_OUT maska dm<$DN6_DM_MAX znak $DN6_SIGN soft=${DN6_SOFTNODE:-0}"
COM=(--data "$DN_OUT/data" --env "$DN2_OUT/parts" --env-dp "$DN4_OUT/parts_dp" --hsc "$DN_OUT" --epochs "$DN6_EPOCHS" --cap-sig "$DN6_CAP_SIG" --cap-bkg "$DN6_CAP_BKG")
body() {
  local b k tag f s n nf i id o mode gate
  if [ "${DRYRUN:-0}" != 1 ]; then
    case $step in
      gtrain) DN2_USE_GPU=1; DN2_CPUS=$DN6_CPUS; export DN2_USE_GPU DN2_CPUS; dn2_env ;;
      *) dn4_env "$DN6_CPUS" ;;
    esac
    mkdir -p "$DN6_OUT"/{models,folds,sc,meta,app/d150,fail} 2>/dev/null || true
  fi
  case $step in
    fold)
      b=$1; k=$2; tag=$(dn6_tag "$b"); [ -n "$tag" ] || { echo "nieznany blok $b"; return 2; }
      if dn6_ready "$tag" && [ "${DN_REDO:-0}" != 1 ]; then echo "$tag już gotowy (meta na EOS) — fold $k pomijam"; return 0; fi
      if [ -s "$DN6_OUT/models/$tag/fold$k.pt" ] && [ -s "$DN6_OUT/folds/$tag/fold$k.json" ] && [ "${DN_REDO:-0}" != 1 ]; then echo "$tag fold $k już jest — pomijam"; return 0; fi
      run python3 -u "$DN6/dn6_train.py" fold "${COM[@]}" --out "$WD/out" --block "$b" --fold "$k" --threads "$DN6_CPUS" --device cpu
      put "$WD/out/models/$tag/fold$k.pt" "$DN6_OUT/models/$tag/fold$k.pt"
      put "$WD/out/folds/$tag/fold$k.json" "$DN6_OUT/folds/$tag/fold$k.json" ;;
    merge)
      b=$1; tag=$(dn6_tag "$b"); [ -n "$tag" ] || { echo "nieznany blok $b"; return 2; }
      if dn6_ready "$tag" && [ "${DN_REDO:-0}" != 1 ]; then echo "$tag już złączony — pomijam"; return 0; fi
      for ((k = 0; k < DN_FOLDS; k++)); do [ -s "$DN6_OUT/models/$tag/fold$k.pt" ] && [ -s "$DN6_OUT/folds/$tag/fold$k.json" ] || { echo "brak foldu $k dla $tag"; return 3; }; done
      run python3 -u "$DN6/dn6_train.py" merge "${COM[@]}" --out "$WD/out" --block "$b" --models "$DN6_OUT/models" --folds-dir "$DN6_OUT/folds" --threads "$DN6_CPUS"
      put "$WD/out/sc/$tag.npz" "$DN6_OUT/sc/$tag.npz"
      put "$WD/out/models/$tag/meta.json" "$DN6_OUT/models/$tag/meta.json"
      put "$WD/out/meta/$tag.json" "$DN6_OUT/meta/$tag.json" ;;                      # meta na końcu = „gotowe”
    gtrain)
      b=$1; tag=$(dn6_tag "$b"); [ -n "$tag" ] || { echo "nieznany blok $b"; return 2; }
      if dn6_ready "$tag" && [ "${DN_REDO:-0}" != 1 ]; then echo "$tag już złączony — pomijam"; return 0; fi
      python3 -c 'import torch; t = torch.randn(64, 8, device="cuda"); print("CUDA OK", torch.cuda.get_device_name(0), float((t @ t.T).sum()))' \
        || { echo "BŁĄD: GPU nie działa na $(hostname) — węzeł gtrain nie liczy na CPU"; return 5; }
      mkdir -p "$WD/out/models/$tag" "$WD/out/folds/$tag"
      for ((k = 0; k < DN_FOLDS; k++)); do
        if [ -s "$DN6_OUT/models/$tag/fold$k.pt" ] && [ -s "$DN6_OUT/folds/$tag/fold$k.json" ] && [ "${DN_REDO:-0}" != 1 ]; then   # fold z wcześniejszej próby (urządzenie zapisane w meta.fit)
          cp "$DN6_OUT/models/$tag/fold$k.pt" "$WD/out/models/$tag/"; cp "$DN6_OUT/folds/$tag/fold$k.json" "$WD/out/folds/$tag/"; echo "$tag fold $k już jest — biorę z EOS"; continue
        fi
        run python3 -u "$DN6/dn6_train.py" fold "${COM[@]}" --out "$WD/out" --block "$b" --fold "$k" --threads "$DN6_CPUS" --device cuda
        grep -q '"device": "cuda"' "$WD/out/folds/$tag/fold$k.json" || { echo "BŁĄD: fold $k policzony nie na GPU"; return 5; }
        put "$WD/out/models/$tag/fold$k.pt" "$DN6_OUT/models/$tag/fold$k.pt"
        put "$WD/out/folds/$tag/fold$k.json" "$DN6_OUT/folds/$tag/fold$k.json"
      done
      if dn6_ready "$tag" && [ "${DN_REDO:-0}" != 1 ]; then echo "$tag złączony w międzyczasie przez inny job — nie nadpisuję"; return 0; fi
      run python3 -u "$DN6/dn6_train.py" merge "${COM[@]}" --out "$WD/out" --block "$b" --models "$WD/out/models" --folds-dir "$WD/out/folds" --threads "$DN6_CPUS" --device cuda
      put "$WD/out/sc/$tag.npz" "$DN6_OUT/sc/$tag.npz"
      put "$WD/out/models/$tag/meta.json" "$DN6_OUT/models/$tag/meta.json"
      put "$WD/out/meta/$tag.json" "$DN6_OUT/meta/$tag.json" ;;
    ana)
      fast=(); [ "$DN_FAST" = 1 ] && fast=(--fast)
      n=0; for b in $DN6_BLOCKS; do if dn6_ready "$(dn6_tag "$b")"; then n=$((n + 1)); fi; done
      echo "gotowych tagów: $n z $(echo "$DN6_BLOCKS" | wc -w)"
      [ "$n" -gt 0 ] || [ "${DRYRUN:-0}" = 1 ] || { echo "brak jakiegokolwiek gotowego tagu"; return 3; }
      run python3 -u "$DN6/dn6_ana.py" --data "$DN_OUT/data" --env "$DN2_OUT/parts" --sc "$DN_OUT,$DN2_OUT,$DN6_OUT" --tags "$DN6_KEEP" --out "$WD/dnet6_p1.json"
      put "$WD/dnet6_p1.json" "$DN6_OUT/dnet6_p1.json"
      run python3 -u "$DN6/dn6_study.py" --data "$DN_OUT/data" --env "$DN2_OUT/parts" --out "$WD/dnet6_mask.json" || echo "UWAGA: dn6_study nieudane (opisowe; nie blokuje)"
      [ -s "$WD/dnet6_mask.json" ] && put "$WD/dnet6_mask.json" "$DN6_OUT/dnet6_mask.json"
      run python3 -u "$DN1/analyze.py" compute --data "$DN_OUT/data" --sc "$DN_OUT,$DN2_OUT,$DN6_OUT" --keep "$DN6_KEEP" --sf-all \
          --out "$WD/dnet6_results.json" --nproc "$DN6_CPUS" --nboot "$DN_NBOOT" "${fast[@]}"
      put "$WD/dnet6_results.json" "$DN6_OUT/dnet6_results.json"
      if [ "${DRYRUN:-0}" != 1 ]; then date -u +%FT%TZ > "$WD/DONE"; put "$WD/DONE" "$DN6_OUT/DONE"; fi ;;
    app)
      s=$1; k=$2; n=$3; mode=d150
      if [ "${DRYRUN:-0}" != 1 ] && [ "${DN6_APP_NOGATE:-0}" != 1 ]; then
        gate=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("p1_gate_seed1"))' "$DN6_OUT/dnet6_p1.json" 2>/dev/null || echo "?")
        [ "$gate" = True ] || { echo "bramka P1 (ziarno 1) = $gate → app nie jest liczone (KRYTERIA_etap6.md §5)"; return 0; }
      fi
      dn6_ready "$DN6_APP_TAG" || [ "${DRYRUN:-0}" = 1 ] || { echo "brak modelu $DN6_APP_TAG"; return 3; }
      [ -s "$DN4_OUT/models/$DN6_APP_REF/meta.json" ] || [ "${DRYRUN:-0}" = 1 ] || { echo "brak drzewa odniesienia $DN4_OUT/models/$DN6_APP_REF"; return 3; }
      mkdir -p "$WD/models"; ln -s "$DN6_OUT/models/$DN6_APP_TAG" "$WD/models/$DN6_APP_TAG"; ln -s "$DN4_OUT/models/$DN6_APP_REF" "$WD/models/$DN6_APP_REF"
      mapfile -t ALL < <(dn4_files "$mode" "$s")
      nf=${#ALL[@]}; [ "$nf" -gt 0 ] || { echo "brak h5 dla $mode $s"; return 2; }
      FILES=(); for ((i = k; i < nf; i += n)); do FILES+=("${ALL[$i]}"); done
      [ ${#FILES[@]} -gt 0 ] || { echo "grupa $k/$n pusta ($nf plików)"; return 0; }
      echo "pliki h5 ($s $k/$n): ${#FILES[@]} z $nf"
      case $s in data_*) tag=data ;; fd_*) tag=fd ;; *) tag=mb ;; esac
      export DN4_PT=3,30 DN4_MR_D0
      for f in "${FILES[@]}"; do
        id=$(basename "$f" | sed -E 's/.*_([0-9]{6})_none\.h5$/\1/')
        case $id in [0-9][0-9][0-9][0-9][0-9][0-9]) ;; *) echo "nieoczekiwana nazwa h5: $f"; return 2 ;; esac
        o="app_${s}_${k}_${id}.npz"
        if [ -s "$DN6_OUT/app/$mode/$o" ] && [ "${DN_REDO:-0}" != 1 ]; then echo "$o już jest — pomijam"; continue; fi
        # katalog pliku: nowy podkatalog scratchu na każdy plik h5 (bez kasowania w pętli; całość sprząta trap przy wyjściu)
        P="$WD/f_$id"; mkdir -p "$P/h5/$s" "$P/data" "$P/parts"; ln -s "$f" "$P/h5/$s/"
        run python3 -u "$DN4/dn4_extract.py" --tag "$tag" --frac-d0 1.0 --frac-dp 1.0 --lxy-min "$DN_LXY" --species D0 --out "$P/data/dn_${s}.npz" "$f"
        run python3 -u "$DN2/dn2_extract.py" --sample "$s" --data "$P/data" --rawdir "$(dn4_rawdir "$mode")/$s" --h5dir "$P/h5" --part 0 --nparts 1 \
            --species D0 --frac-a 1.0 --out "$P/parts/dn2_${s}_0.npz"
        run python3 -u "$DN6/dn6_apply.py" --sample "$s" --data "$P/data" --env "$P/parts" --models "$WD/models" --species D0 --threads "$DN6_CPUS" --out "$P/$o"
        put "$P/$o" "$DN6_OUT/app/$mode/$o"
        : > "$P/data/dn_${s}.npz"; : > "$P/parts/dn2_${s}_0.npz"; : > "$P/$o"          # zwolnienie dysku scratch (pliki obcięte do 0 B)
      done ;;
    anar)
      # 20.11: liczenie odporne na EOS FUSE — dawne `ls …/app_*.npz | wc -l` pod `set -eo pipefail` kończyło węzeł kodem ls (rc = 2), gdy świeżo zapisane
      # pliki nie były jeszcze widoczne/stat-owalne na węźle (job 12794123); teraz: listing katalogu, bez statów, do 6 prób co 30 s
      nf=0
      for i in 1 2 3 4 5 6; do
        nf=$( { ls "$DN6_OUT/app/d150" 2>/dev/null || true; } | grep -c '^app_.*\.npz$' || true )
        [ "${nf:-0}" -gt 0 ] && break
        echo "app/d150: brak widocznych plików (próba $i/6) — czekam 30 s"; sleep 30
      done
      echo "app/d150: $nf plików app"
      [ "$nf" -gt 0 ] || { echo "brak plików app (bramka P1 zamknięta albo app nieudane) — anar nic nie liczy"; return 0; }
      run python3 -u "$DN4/dn4_ana.py" r --app "$DN6_OUT/app/d150" --species D0 --force-data --out "$WD/dnet6_R.json" --nboot "$DN4_NBOOT"
      put "$WD/dnet6_R.json" "$DN6_OUT/dnet6_R.json"
      if [ "${DRYRUN:-0}" != 1 ]; then date -u +%FT%TZ > "$WD/DONE_R"; put "$WD/DONE_R" "$DN6_OUT/DONE_R"; fi ;;
    *) echo "nieznany krok $step"; return 2 ;;
  esac
  return 0
}
set +e
( set -eo pipefail; body "$@" )        # podpowłoka POZA listą ||/&& — inaczej bash ignoruje set -e w jej wnętrzu
rc=$?
if [ $rc -ne 0 ] && [ "${DN6_SOFTNODE:-0}" = 1 ]; then
  echo "[$(date -u +%FT%TZ)] MIĘKKI BŁĄD DNET6 $NODE rc=$rc — znacznik fail/$NODE, DAG idzie dalej"
  echo "$(date -u +%FT%TZ) rc=$rc host=$(hostname)" > "$WD/fail_$NODE"; dn_put "$WD/fail_$NODE" "$DN6_OUT/fail/$NODE" || true
  exit 0
fi
[ $rc -eq 0 ] && echo "[$(date -u +%FT%TZ)] OK DNET6 $NODE ($(( $(date +%s) - T_START )) s)"
exit $rc

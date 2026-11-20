#!/bin/bash
# TAGXSEC (20.11.2026): score drzewa H dla mb_sd — ten sam łańcuch co węzeł `dn4_job.sh app d150` (dn4_extract → dn2_extract → dn4_apply),
# tylko modele .xgb., wyjście do /tmp/klyczak_tagxsec/app (nic nie trafia do charm_q/dnet4). Interaktywnie, 2 wątki.
set -eo pipefail
S=${1:-mb_sd}
DN4=$HOME/private/charm_q/Q_lxplus/dnet4; export DN4
source "$DN4/dn4_common.sh"
W=/tmp/klyczak_tagxsec; mkdir -p "$W/app" "$W/wd"
DN4_APP_CPUS=2; dn4_env 2
export DN4_PT=3,30 DN4_MR_D0
for f in $(dn4_files d150 "$S"); do
  id=$(basename "$f" | sed -E 's/.*_([0-9]{6})_none\.h5$/\1/'); o="$W/app/app_${S}_0_${id}.npz"
  [ -s "$o" ] && { echo "$o jest"; continue; }
  WD=$W/wd/$id; mkdir -p "$WD/h5/$S" "$WD/data" "$WD/parts"; ln -sf "$f" "$WD/h5/$S/"
  python3 -u "$DN4/dn4_extract.py" --tag mb --frac-d0 1.0 --frac-dp 1.0 --lxy-min "$DN_LXY" --out "$WD/data/dn_${S}.npz" "$f"
  python3 -u "$DN2/dn2_extract.py" --sample "$S" --data "$WD/data" --rawdir "$RAWV20/$S" --h5dir "$WD/h5" --part 0 --nparts 1 \
      --species D0,Dplus --frac-a 1.0 --out "$WD/parts/dn2_${S}_0.npz"
  python3 -u "$DN4/dn4_apply.py" --sample "$S" --data "$WD/data" --env "$WD/parts" --models "$DN4_OUT/models" --tags '\.xgb\.' --threads 2 --out "$o"
done
echo "OK $(ls $W/app | wc -l) plików"

#!/bin/bash
# l1: lista tagów świetlności w COOL (odczyt przez Frontier) + pliki GRL/ilumicalc 2022 na cvmfs.
cd /tmp/klyczak_lumi || exit 1
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh -q
asetup Athena,24.0,latest > asetup.log 2>&1
G=/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists
echo "== GRL data22"; for d in $G/data22_13p6TeV/*; do echo "-- $d"; ls $d; done
echo "== grep runs in GRL xml"; grep -l -E "<Run>(427929|428770)</Run>" $G/data22_13p6TeV/*/*.xml $G/*/*/*.xml 2>/dev/null | sort -u
python - <<'PY'
from PyCool import cool
from CoolConvUtilities.AtlCoolLib import indirectOpen
for dbn in ("COOLOFL_TRIGGER/CONDBR2",):
    db = indirectOpen(dbn, True, True)
    for fn in ("/TRIGGER/OFLLUMI/OflPrefLumi", "/TRIGGER/OFLLUMI/LumiAccounting"):
        try:
            f = db.getFolder(fn); tags = sorted(str(t) for t in f.listTags())
            print("FOLDER", fn, "ntags", len(tags))
            for t in tags:
                try: print("  ", t, "| inserted", f.tagInsertionTime(t), "|", f.tagDescription(t))
                except Exception as e: print("  ", t)
        except Exception as e: print("ERR", fn, e)
PY

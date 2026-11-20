#!/bin/bash
# l2: iLumiCalc (odczyt COOL przez Frontier) dla jednego runu/GRL, listy łańcuchów i jednego livetrigger.
# Użycie: l2_ilumicalc.sh <etykieta> <tag> <livetrigger> <xml> <run> <trig1> [trig2 ...]   (limit 9 min)
L=$1; TAG=$2; LIVE=$3; XML=$4; RUN=$5; shift 5
W=/tmp/klyczak_lumi/$L; mkdir -p $W; cd $W || exit 1
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh -q
asetup Athena,24.0,latest > asetup.log 2>&1
# GRL ograniczony do jednego runu
python - "$XML" "$RUN" > grl_$RUN.xml <<'PY'
import sys, re
x = open(sys.argv[1]).read(); run = sys.argv[2]
blocks = re.findall(r"<LumiBlockCollection>.*?</LumiBlockCollection>", x, re.S)
keep = [b for b in blocks if f"<Run>{run}</Run>" in b]
head = x.split("<LumiBlockCollection>")[0]
print(head + "\n".join(keep) + "\n</NamedLumiRange>\n</LumiRangeCollection>")
PY
T=""; for t in "$@"; do T="$T --trigger=$t"; done
echo "CMD: iLumiCalc --lumitag=$TAG --livetrigger=$LIVE $T --xml=grl_$RUN.xml --plots --verbose" > cmd.txt; date -u >> cmd.txt
timeout 540 iLumiCalc --lumitag=$TAG --livetrigger=$LIVE $T --xml=grl_$RUN.xml --plots --verbose > ilumicalc.log 2>&1
echo "rc=$?" >> cmd.txt; date -u >> cmd.txt
cat cmd.txt; ls -la; grep -n -i -E "error|fatal|not found|Total|IntL|Recorded|Delivered|L1 ?Acc|L3 ?Acc|Livetime|Prescale|Good LBs|Bad LBs|Lumi ?tag|channel" ilumicalc.log | tail -60

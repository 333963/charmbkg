#!/bin/bash
# l0: środowisko + pomoc iLumiCalc (tylko odczyt). Uruchamiane na lxplus w /tmp/klyczak_lumi.
cd /tmp/klyczak_lumi || exit 1
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh -q
asetup Athena,24.0,latest > asetup.log 2>&1
tail -2 asetup.log
which iLumiCalc iLumiCalc.exe 2>&1
iLumiCalc --help > ilumicalc_help.txt 2>&1; echo "help rc=$?"
cat ilumicalc_help.txt

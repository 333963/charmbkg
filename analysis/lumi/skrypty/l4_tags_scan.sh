#!/bin/bash
# l4: iLumiCalc bez triggera dla GRL analizy, kolejne tagi świetlności (sprawdzenie, które istnieją i czy L_del się różni). Limit 8 min.
W=/tmp/klyczak_lumi/tagscan; mkdir -p $W; cd $W || exit 1
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh -q
asetup Athena,24.0,latest > asetup.log 2>&1
for TAG in "$@"; do
  mkdir -p $TAG; ( cd $TAG; timeout 110 iLumiCalc --lumitag=$TAG --livetrigger=L1_EM22VHI --xml=/tmp/klyczak_lumi/GRL.xml > ilumicalc.log 2>&1; echo "TAG $TAG rc=$?"
    sed 's/\x1b\[[0-9;]*m//g' ilumicalc.log | grep -E "ERROR|IntL delivered|Total IntL|Good LBs|Bad LBs|^--- LumiCalculator      : +4[0-9]{5} " | cut -c1-200 )
done

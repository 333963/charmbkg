#!/bin/bash
# l3_run.sh <run> <items>  (limit 8 min)
mkdir -p /tmp/klyczak_lumi/cool; cd /tmp/klyczak_lumi/cool || exit 1
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh -q
asetup Athena,24.0,latest > asetup.log 2>&1
timeout 480 python /tmp/klyczak_lumi/l3_cool_counters.py "$@" 2>&1 | grep -v "^Warning in\|TClass::Init" ; echo "rc=${PIPESTATUS[0]}"

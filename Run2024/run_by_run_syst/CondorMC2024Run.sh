#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.460
export HOME=/sphenix/u/egm2153
export MYINSTALL=$HOME/install
source $OPT_SPHENIX/bin/setup_local.sh $MYINSTALL
export SPHENIX=$MYINSTALL
root.exe -q -b 'dETdeta_multi_run_analysis.C("'$1'",-10,10,1,1,1,'$2','$3',0,0,0,"","","","'$4'","")'
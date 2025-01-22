#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.450
export HOME=/sphenix/u/egm2153
export MYINSTALL=$HOME/install
source $OPT_SPHENIX/bin/setup_local.sh $MYINSTALL
export SPHENIX=$MYINSTALL
root.exe -q -b 'dETdeta_multi_run_analysis.C("'$1'",'$2','$3','$4','$5',1,0,5,0,0,0,"","","","'$6'","'$7'")'
#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n
export HOME=/sphenix/u/egm2153
export MYINSTALL=$HOME/install
source $OPT_SPHENIX/bin/setup_local.sh $MYINSTALL
export SPHENIX=$MYINSTALL

echo "------------------setting up environment--------------------"
export Cur_dir=$(pwd)
echo "running area:" ${Cur_dir}
echo "-------------------------running----------------------------"
cd ${Cur_dir}
ls
date
# $1 = runnumber, $2 = generator, $3 = minus_z, $4 = plus_z, $5 = dataormc, $6 = reweighting, $7 = central
root.exe -q -b "dETdeta_vertex_reweighting.C(\"$1\", \"$2\")"
root.exe -q -b "dETdeta_hot_dead_map.C(\"$1\", \"$3\", \"$4\")"
root.exe -q -b "dETdeta_analysis.C(\"$1\", \"$2\", \"$3\", \"$4\", \"$5\", \"$6\", \"$7\")"
date
echo "JOB COMPLETE!"
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

root.exe -q -b 'dETdeta_analysis.C('$1',"",-20,20,0,0,1,0,5,1,30,0,"4_4_24")' > "run_by_run_studies/run$1_log.txt"
root.exe -q -b 'dETdeta_analysis.C('$1',"",-20,20,0,0,1,5,10,1,30,0,"4_4_24")' > "run_by_run_studies/run$1_log.txt"
root.exe -q -b 'dETdeta_analysis.C('$1',"",-20,20,0,0,1,10,20,1,30,0,"4_4_24")' > "run_by_run_studies/run$1_log.txt"
root.exe -q -b 'dETdeta_analysis.C('$1',"",-20,20,0,0,1,20,30,1,30,0,"4_4_24")' > "run_by_run_studies/run$1_log.txt"
root.exe -q -b 'dETdeta_analysis.C('$1',"",-20,20,0,0,1,30,40,1,30,0,"4_4_24")' > "run_by_run_studies/run$1_log.txt"
root.exe -q -b 'dETdeta_analysis.C('$1',"",-20,20,0,0,1,40,50,1,30,0,"4_4_24")' > "run_by_run_studies/run$1_log.txt"
root.exe -q -b 'dETdeta_analysis.C('$1',"",-20,20,0,0,1,50,60,1,30,0,"4_4_24")' > "run_by_run_studies/run$1_log.txt"
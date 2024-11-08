#!/bin/bash

root.exe -q -b 'dETdeta_multi_run_analysis.C("'$1'",-20,20,1,1,1,30,40,0,0,0,"","p015","MC")' >> "MC/$1_ext_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("'$1'",-20,20,1,1,1,40,50,0,0,0,"","p015","MC")' >> "MC/$1_ext_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("'$1'",-20,20,1,1,1,50,60,0,0,0,"","p015","MC")' >> "MC/$1_ext_log.txt"
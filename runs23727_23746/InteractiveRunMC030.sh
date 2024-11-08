#!/bin/bash

#root.exe -q -b 'dETdeta_multi_run_analysis.C("'$1'",-20,20,1,1,1,0,5,0,0,0,"","p015","MC")' >> "MC/$1_log.txt"
#root.exe -q -b 'dETdeta_multi_run_analysis.C("'$1'",-20,20,1,1,1,5,10,0,0,0,"","p015","MC")' >> "MC/$1_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("'$1'",-20,20,1,1,1,10,20,0,0,0,"","p015","MC")' >> "MC/$1_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("'$1'",-20,20,1,1,1,20,30,0,0,0,"","p015","MC")' >> "MC/$1_log.txt"

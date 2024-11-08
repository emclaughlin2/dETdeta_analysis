#!/bin/bash

root.exe -q -b 'dETdeta_analysis.C('$1',"reweight_hijing",-20,20,1,1,1,0,5,0,0,0,"4_4_24")' >> "run_by_run_studies/hijing_run$1_log.txt"
root.exe -q -b 'dETdeta_analysis.C('$1',"reweight_hijing",-20,20,1,1,1,5,10,0,0,0,"4_4_24")' >> "run_by_run_studies/hijing_run$1_log.txt"
root.exe -q -b 'dETdeta_analysis.C('$1',"reweight_hijing",-20,20,1,1,1,10,20,0,0,0,"4_4_24")' >> "run_by_run_studies/hijing_run$1_log.txt"
root.exe -q -b 'dETdeta_analysis.C('$1',"reweight_hijing",-20,20,1,1,1,20,30,0,0,0,"4_4_24")' >> "run_by_run_studies/hijing_run$1_log.txt"
#root.exe -q -b 'dETdeta_analysis.C('$1',"reweight_hijing",-20,20,1,1,1,30,40,0,0,0,"4_4_24")' >> "run_by_run_studies/hijing_run$1_log.txt"
#root.exe -q -b 'dETdeta_analysis.C('$1',"reweight_hijing",-20,20,1,1,1,40,50,0,0,0,"4_4_24")' >> "run_by_run_studies/hijing_run$1_log.txt"
#root.exe -q -b 'dETdeta_analysis.C('$1',"reweight_hijing",-20,20,1,1,1,50,60,0,0,0,"4_4_24")' >> "run_by_run_studies/hijing_run$1_log.txt"
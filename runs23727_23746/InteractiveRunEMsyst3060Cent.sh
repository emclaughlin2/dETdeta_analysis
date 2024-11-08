#!/bin/bash

root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,30,40,1,30,0,"'$1'","p015","emcal_syst")' >> "emcal_syst/syst0$1_ext_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,40,50,1,30,0,"'$1'","p015","emcal_syst")' >> "emcal_syst/syst0$1_ext_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,50,60,1,30,0,"'$1'","p015","emcal_syst")' >> "emcal_syst/syst0$1_ext_log.txt"

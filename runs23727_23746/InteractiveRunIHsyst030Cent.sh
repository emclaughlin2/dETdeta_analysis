#!/bin/bash

root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,0,5,1,30,0,"","'$1'","","p015","hcal_syst")' >> "hcal_syst/ihsyst0$1_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,5,10,1,30,0,"","'$1'","","p015","hcal_syst")' >> "hcal_syst/ihsyst0$1_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,10,20,1,30,0,"","'$1'","","p015","hcal_syst")' >> "hcal_syst/ihsyst0$1_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,20,30,1,30,0,"","'$1'","","p015","hcal_syst")' >> "hcal_syst/ihsyst0$1_log.txt"
#!/bin/bash

root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,30,40,1,'$1',0,"","p015","zs_syst")' >> "zs_syst/syst$1_ext_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,40,50,1,'$1',0,"","p015","zs_syst")' >> "zs_syst/syst$1_ext_log.txt"
root.exe -q -b 'dETdeta_multi_run_analysis.C("",-20,20,0,0,1,50,60,1,'$1',0,"","p015","zs_syst")' >> "zs_syst/syst$1_ext_log.txt"

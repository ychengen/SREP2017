###################################################################
# Script:	step3-GSWP3.sh
# Purpose:	Compute climate means from GSWP3 data
# 
# Created by:	Cheng-En Yang @ University of Tennessee
# Last update:	March 26, 2018
###################################################################
#!/bin/bash

###################################################################
## Regrid GSWP3 data
###################################################################
cd step3-GSWP3
./step1-run_regrid.sh
wait

###################################################################
## Compute climate mean
###################################################################
./step2-run_climate_mean.sh
wait
cd ../


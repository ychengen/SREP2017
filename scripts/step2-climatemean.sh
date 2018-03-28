###################################################################
# Script:	step2-CMIP5_master.sh
# Purpose:	Compute annul and climate means from CMIP5 outputs
# 
# Created by:	Cheng-En Yang @ University of Tennessee
# Last update:	March 26, 2018
###################################################################
#!/bin/bash

module load nco

###################################################################
## Annual mean
###################################################################
## declare an array storing ESM names
declare -a arr=("BNU-ESM" "HadGEM2-CC" "HadGEM2-ES" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM" "MIROC-ESM-CHEM")

## output directory for the contemporary period (1982-2005)
if [ ! -d CMIP5/annual_mean/modern ]; then
   mkdir -p CMIP5/annual_mean/modern
fi

cd step2-climatemean
./wrapper_run_modern.sh
cd ..

## output directory for the pre-industrial period (1861-1885)
if [ ! -d CMIP5/annual_mean/PI ]; then
  mkdir -p CMIP5/annual_mean/PI
fi

cd step2-climatemean
./wrapper_run_PI.sh
cd ..

wait

###################################################################
## Climate mean
###################################################################

cd CMIP5/annual_mean/

### modern period (1982-2005)
cd modern

for index in "${arr[@]}"
   do
      echo "Computing MODERN climate mean of $index"
      ncra -F -d time,1,24 $index.nc temp1_$index.nc
      ncwa -O -a time temp1_${index}.nc temp2_${index}.nc
      ncks -a -x -v height,time,time_bnds temp2_${index}.nc climate_${index}.nc
done

wait

rm -f temp1_*.nc temp2_*.nc

if [ ! -d ../../climate_mean/modern ]; then
  mkdir -p ../../climate_mean/modern
fi
mv climate_*.nc ../../climate_mean/modern
cd ..

## pre-industrial period (1861-1885)
cd PI

for index in "${arr[@]}"
   do
      echo "Computing PI climate mean of $index"
      ncra -F -d time,1,25 $index.nc temp1_$index.nc
      ncwa -O -a time temp1_${index}.nc temp2_${index}.nc
      ncks -a -x -v height,time,time_bnds temp2_${index}.nc climate_${index}.nc
done

wait

rm -f temp1_*.nc temp2_*.nc
if [ ! -d ../../climate_mean/PI ]; then
  mkdir -p ../../climate_mean/PI
fi
mv climate_*.nc ../../climate_mean/PI


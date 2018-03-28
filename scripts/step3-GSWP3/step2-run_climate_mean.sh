##########################
###  For climate mean  ###
##########################
#!/bin/bash

module load nco

cd ../../Observation/GSWP3

if [ ! -d climate_mean/modern ]; then
   mkdir -p climate_mean/modern
fi
if [ ! -d climate_mean/PI ]; then
   mkdir -p climate_mean/PI
fi

## declare an array variable
declare -a arr=("BNU-ESM" "HadGEM2-CC" "HadGEM2-ES" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM" "MIROC-ESM-CHEM")

## now loop through the above array to extract data during 1982-2005
## modern time period (1982-2005)
cd annual_mean
for index in "${arr[@]}"
   do
      ncks -a -F -d time,133,156 GSWP3_$index.nc temp_1_$index.nc
      ncwa -O -a time temp_1_$index.nc temp_2_$index.nc
      ncks -v time -x temp_2_$index.nc climate_$index.nc
done

mv climate_*.nc ../climate_mean/modern/

rm -f out_*.nc temp_*.nc

## pre-industrial time period (1861-1885)
for index in "${arr[@]}"
   do
      ncks -a -F -d time,12,36 GSWP3_$index.nc temp_1_$index.nc
      ncwa -O -a time temp_1_$index.nc temp_2_$index.nc
      ncks -v time -x temp_2_$index.nc climate_$index.nc
done

mv climate_*.nc ../climate_mean/PI/

rm -f out_*.nc temp_*.nc


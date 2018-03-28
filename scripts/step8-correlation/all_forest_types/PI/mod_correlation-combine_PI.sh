#!/bin/bash

module load nco

cd ../../../../

if [ ! -d postanalysis ]; then
   mkdir -p postanalysis
fi
cd postanalysis
if [ ! -d correlation/all_forest_types ]
   mkdir -p correlation/all_forest_types
fi
cd correlation

#declare -a arr1=("0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9")
declare -a arr1=("0" "0.2" "0.5" "0.8")
declare -a arr2=("BNU-ESM" "HadGEM2-CC" "HadGEM2-ES" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM" "MIROC-ESM-CHEM")

for index1 in "${arr1[@]}"
   do
      mkdir -p all_forest_types/PI/frac$index1
done

cd ../../scripts/step8-correlation/all_forest_types/PI
ncl mod_correlation_PI-1.ncl
ncl mod_correlation_PI-2.ncl

cd ../../../../postanalysis/correlation/all_forest_types/PI/

for index1 in "${arr1[@]}"
   do
      cd frac$index1

      for index2 in "${arr2[@]}"
         do
            echo "$index2"
            ncks -A correlation-${index2}-MOD-2.nc correlation-${index2}-MOD-1.nc
            mv correlation-${index2}-MOD-1.nc correlation-${index2}-MODEL.nc
      done
      rm -f correlation*MOD-2.nc
      cd ..
done

echo "Done..."

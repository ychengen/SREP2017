#!/bin/bash

module load ncl

declare -a arr=("BNU-ESM" "HadGEM2-CC" "HadGEM2-ES" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM" "MIROC-ESM-CHEM")

cd ../../Observation/GSWP3

if [ ! -d annual_mean ]; then
   mkdir -p annual_mean
fi
if [ ! -d climate_mean ]; then
   mkdir -p climate_mean
fi

cd ../../scripts/step3-GSWP3

for index in "${arr[@]}"
   do 
      echo 'Processing ' ${index}
      sed -e "s/modelname/${index}/g" regrid_template.ncl > create_regrid_${index}.ncl
      ncl create_regrid_${index}.ncl
      rm -f create_regrid_${index}.ncl
done

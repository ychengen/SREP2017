#####################################################################
# Script:	step1-timeslice.sh
# Purpose:	Extracting contemporary and pre-industrial data from
#		CMIP5 outputs
# 
# Created by:	Cheng-En Yang @ University of Tennessee
# Last update:	March 26, 2018
#####################################################################
#!/bin/bash

module load nco

## declare an array storing ESM names
declare -a arr=("BNU-ESM" "HadGEM2-CC" "HadGEM2-ES" "IPSL-CM5A-LR" "IPSL-CM5A-MR" "IPSL-CM5B-LR" "MIROC-ESM" "MIROC-ESM-CHEM")

## check input files from CMIP5
for i in "${arr[@]}"
   do
      if [ ! -d CMIP5/${i} ]; then
         echo "Error! Input folder \"${i}\" does not exist under /CMIP5..."
         echo ""
         exit 1
      fi
done

## output directory for the contemporary period (1982-2005)

for i in "${arr[@]}"
   do
      echo "Processing ${i} - modern..."
      if [ ! -d CMIP5/$i/1982-2005 ]; then
         mkdir -p CMIP5/$i/1982-2005
      fi

      cd step1-timeslice
      ./nco_extract-${i}.sh
      cd ../
done

## output directory for the pre-industrial period (1861-1885)

for i in "${arr[@]}"
   do
      echo "Processing ${i} - PI..."
      if [ ! -d CMIP5/$i/1861-1885 ]; then
         mkdir -p CMIP5/$i/1861-1885
      fi

      cd step1-timeslice
      ./nco_extract-${i}-pi.sh
      cd ../
done


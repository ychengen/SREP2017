#######################################################################
# Script:       step6-upscale.sh
# Purpose:      Upscale the remapped GLC2000-MPIBGI observational data
#		to the grid resolution of each ESM
#
# Note:		Important!! Modify the input paths before running the
#		programs!
# 
# Created by:   Cheng-En Yang @ University of Tennessee
# Last update:  March 26, 2018
#######################################################################
#!/bin/bash

## Create upscaled land cover types according to the grid resolution
##  of each ESM

if [ ! -d ../CMIP5/landcover_GLC2000_for_CMIP5 ]; then
   mkdir -p ../CMIP5/landcover_GLC2000_for_CMIP5
fi

## Important: modify the input and output paths accordingly in the
##            Fortran files residing in each compartment folder first!!

cd step6-upscale
./run.csh

## compute min and max range based on the uncertainty 
if [ ! -d ../CMIP5/landcover_GLC2000_for_CMIP5/uncertainty ]; then
   mkdir -p ../CMIP5/landcover_GLC2000_for_CMIP5/uncertainty/max
   mkdir -p ../CMIP5/landcover_GLC2000_for_CMIP5/uncertainty/min
fi

./run_uncertainty_max.csh
./run_uncertainty_min.csh


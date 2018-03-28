#####################################################################
# Script:	run_whittaker.sh
# Purpose:	scatter plots on a temperature-precipitation
#		space
# Note:		Must use NCAR NCL verion 6.4.0 or higher
# Created by:	Cheng-En Yang @ University of Tennessee at Knoxville
# Last modified:03/27/2018
#####################################################################
#!/bin/bash

module load ncl/6.4.0

## Create observation and model plots
ncl whittaker_pft_modify_threshold.ncl

## Create obsevaional biomass on the modeled climate space
ncl whittaker_pft_modify_threshold_ObsClim_ModBio.ncl

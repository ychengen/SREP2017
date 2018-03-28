#!/bin/bash

module load nco

cd CMIP5/BNU-ESM

ncea -F -d time,133,432 -v cLeaf cLeaf_Lmon_BNU-ESM_historical_r1i1p1_185001-200512.nc  cLeaf_Lmon_BNU-ESM_historical_r1i1p1_186101-188512.nc
ncea -F -d time,133,432 -v cWood cWood_Lmon_BNU-ESM_historical_r1i1p1_185001-200512.nc  cWood_Lmon_BNU-ESM_historical_r1i1p1_186101-188512.nc
ncea -F -d time,133,432 -v cVeg  cVeg_Lmon_BNU-ESM_historical_r1i1p1_185001-200512.nc   cVeg_Lmon_BNU-ESM_historical_r1i1p1_186101-188512.nc
ncea -F -d time,133,432 -v cRoot cRoot_Lmon_BNU-ESM_historical_r1i1p1_185001-200512.nc  cRoot_Lmon_BNU-ESM_historical_r1i1p1_186101-188512.nc
ncea -F -d time,133,432 -v landCoverFrac landCoverFrac_Lmon_BNU-ESM_historical_r1i1p1_185001-200512.nc  landCoverFrac_Lmon_BNU-ESM_historical_r1i1p1_186101-188512.nc
ncea -F -d time,133,432 -v treeFrac treeFrac_Lmon_BNU-ESM_historical_r1i1p1_185001-200512.nc  treeFrac_Lmon_BNU-ESM_historical_r1i1p1_186101-188512.nc
ncea -F -d time,133,432 -v pr    pr_Amon_BNU-ESM_historical_r1i1p1_185001-200512.nc     pr_Amon_BNU-ESM_historical_r1i1p1_186101-188512.nc
ncea -F -d time,133,432 -v rsds  rsds_Amon_BNU-ESM_historical_r1i1p1_185001-200512.nc   rsds_Amon_BNU-ESM_historical_r1i1p1_186101-188512.nc
ncea -F -d time,133,432 -v tas   tas_Amon_BNU-ESM_historical_r1i1p1_185001-200512.nc    tas_Amon_BNU-ESM_historical_r1i1p1_186101-188512.nc

mv *186101-188512.nc ./1861-1885/

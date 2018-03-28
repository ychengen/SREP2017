#!/bin/bash

module load nco

cd CMIP5/MIROC-ESM

ncea -F -d time,1585,1872 -v cLeaf cLeaf_Lmon_MIROC-ESM_historical_r1i1p1_185001-200512.nc  cLeaf_Lmon_MIROC-ESM_historical_r1i1p1_198201-200512.nc
ncea -F -d time,1585,1872 -v cWood cWood_Lmon_MIROC-ESM_historical_r1i1p1_185001-200512.nc  cWood_Lmon_MIROC-ESM_historical_r1i1p1_198201-200512.nc
ncea -F -d time,1585,1872 -v cRoot cRoot_Lmon_MIROC-ESM_historical_r1i1p1_185001-200512.nc  cRoot_Lmon_MIROC-ESM_historical_r1i1p1_198201-200512.nc
ncea -F -d time,1585,1872 -v cVeg  cVeg_Lmon_MIROC-ESM_historical_r1i1p1_185001-200512.nc   cVeg_Lmon_MIROC-ESM_historical_r1i1p1_198201-200512.nc
ncea -F -d time,1585,1872 -v landCoverFrac  landCoverFrac_Lmon_MIROC-ESM_historical_r1i1p1_185001-200512.nc   landCoverFrac_Lmon_MIROC-ESM_historical_r1i1p1_198201-200512.nc
ncea -F -d time,1585,1872 -v treeFrac  treeFrac_Lmon_MIROC-ESM_historical_r1i1p1_185001-200512.nc   treeFrac_Lmon_MIROC-ESM_historical_r1i1p1_198201-200512.nc
ncea -F -d time,1585,1872 -v pr    pr_Amon_MIROC-ESM_historical_r1i1p1_185001-200512.nc     pr_Amon_MIROC-ESM_historical_r1i1p1_198201-200512.nc
ncea -F -d time,1585,1872 -v rsds  rsds_Amon_MIROC-ESM_historical_r1i1p1_185001-200512.nc   rsds_Amon_MIROC-ESM_historical_r1i1p1_198201-200512.nc
ncea -F -d time,1585,1872 -v tas   tas_Amon_MIROC-ESM_historical_r1i1p1_185001-200512.nc    tas_Amon_MIROC-ESM_historical_r1i1p1_198201-200512.nc

mv *1982*.nc ./1982-2005/

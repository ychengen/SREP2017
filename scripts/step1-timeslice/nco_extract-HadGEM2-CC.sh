#!/bin/bash

module load nco

cd CMIP5/HadGEM2-CC

ncea -F -d time,265,300 -v cLeaf cLeaf_Lmon_HadGEM2-CC_historical_r1i1p1_195912-198411.nc temp1.nc
ncrcat temp1.nc cLeaf_Lmon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc cLeaf_Lmon_HadGEM2-CC_historical_r1i1p1_198112-200511.nc

ncea -F -d time,265,300 -v cWood cWood_Lmon_HadGEM2-CC_historical_r1i1p1_195912-198411.nc temp2.nc
ncrcat temp2.nc cWood_Lmon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc cWood_Lmon_HadGEM2-CC_historical_r1i1p1_198112-200511.nc

ncea -F -d time,265,300 -v cRoot cRoot_Lmon_HadGEM2-CC_historical_r1i1p1_195912-198411.nc temp3.nc
ncrcat temp3.nc cRoot_Lmon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc cRoot_Lmon_HadGEM2-CC_historical_r1i1p1_198112-200511.nc

ncea -F -d time,265,300 -v cVeg cVeg_Lmon_HadGEM2-CC_historical_r1i1p1_195912-198411.nc temp4.nc
ncrcat temp4.nc cVeg_Lmon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc cVeg_Lmon_HadGEM2-CC_historical_r1i1p1_198112-200511.nc

ncea -F -d time,265,300 -v landCoverFrac landCoverFrac_Lmon_HadGEM2-CC_historical_r1i1p1_195912-198411.nc temp5.nc
ncrcat temp5.nc landCoverFrac_Lmon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc landCoverFrac_Lmon_HadGEM2-CC_historical_r1i1p1_198112-200511.nc

ncea -F -d time,265,300 -v treeFrac treeFrac_Lmon_HadGEM2-CC_historical_r1i1p1_195912-198411.nc temp6.nc
ncrcat temp6.nc treeFrac_Lmon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc treeFrac_Lmon_HadGEM2-CC_historical_r1i1p1_198112-200511.nc

ncea -F -d time,265,300 -v pr pr_Amon_HadGEM2-CC_historical_r1i1p1_195912-198411.nc temp7.nc
ncrcat temp7.nc pr_Amon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc pr_Amon_HadGEM2-CC_historical_r1i1p1_198112-200511.nc

ncea -F -d time,265,300 -v rsds rsds_Amon_HadGEM2-CC_historical_r1i1p1_195912-198411.nc temp8.nc
ncrcat temp8.nc rsds_Amon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc rsds_Amon_HadGEM2-CC_historical_r1i1p1_198112-200511.nc

ncea -F -d time,265,300 -v tas tas_Amon_HadGEM2-CC_historical_r1i1p1_195912-198411.nc temp9.nc
ncrcat temp9.nc tas_Amon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc tas_Amon_HadGEM2-CC_historical_r1i1p1_198112-200511.nc

mv *198112*.nc ./1982-2005/
rm temp1.nc temp2.nc temp3.nc temp4.nc temp5.nc temp6.nc temp7.nc temp8.nc temp9.nc


#!/bin/bash

module load nco

cd CMIP5/HadGEM2-CC

ncea -F -d time,13,300 -v cLeaf cLeaf_Lmon_HadGEM2-CC_historical_r1i1p1_185912-188411.nc temp11.nc
ncea -F -d time,1,12 -v cLeaf cLeaf_Lmon_HadGEM2-CC_historical_r1i1p1_188412-190911.nc temp12.nc
ncrcat temp11.nc temp12.nc cLeaf_Lmon_HadGEM2-CC_historical_r1i1p1_186012-188511.nc

ncea -F -d time,13,300 -v cWood cWood_Lmon_HadGEM2-CC_historical_r1i1p1_185912-188411.nc temp21.nc
ncea -F -d time,1,12 -v cWood cWood_Lmon_HadGEM2-CC_historical_r1i1p1_188412-190911.nc temp22.nc
ncrcat temp21.nc temp22.nc cWood_Lmon_HadGEM2-CC_historical_r1i1p1_186012-188511.nc

ncea -F -d time,13,300 -v cRoot cRoot_Lmon_HadGEM2-CC_historical_r1i1p1_185912-188411.nc temp31.nc
ncea -F -d time,1,12 -v cRoot cRoot_Lmon_HadGEM2-CC_historical_r1i1p1_188412-190911.nc temp32.nc
ncrcat temp31.nc temp32.nc cRoot_Lmon_HadGEM2-CC_historical_r1i1p1_186012-188511.nc

ncea -F -d time,13,300 -v cVeg cVeg_Lmon_HadGEM2-CC_historical_r1i1p1_185912-188411.nc temp41.nc
ncea -F -d time,1,12 -v cVeg cVeg_Lmon_HadGEM2-CC_historical_r1i1p1_188412-190911.nc temp42.nc
ncrcat temp41.nc temp42.nc cVeg_Lmon_HadGEM2-CC_historical_r1i1p1_186012-188511.nc

ncea -F -d time,13,300 -v landCoverFrac landCoverFrac_Lmon_HadGEM2-CC_historical_r1i1p1_185912-188411.nc temp51.nc
ncea -F -d time,1,12 -v landCoverFrac landCoverFrac_Lmon_HadGEM2-CC_historical_r1i1p1_188412-190911.nc temp52.nc
ncrcat temp51.nc temp52.nc landCoverFrac_Lmon_HadGEM2-CC_historical_r1i1p1_186012-188511.nc

ncea -F -d time,13,300 -v treeFrac treeFrac_Lmon_HadGEM2-CC_historical_r1i1p1_185912-188411.nc temp61.nc
ncea -F -d time,1,12 -v treeFrac treeFrac_Lmon_HadGEM2-CC_historical_r1i1p1_188412-190911.nc temp62.nc
ncrcat temp61.nc temp62.nc treeFrac_Lmon_HadGEM2-CC_historical_r1i1p1_186012-188511.nc

ncea -F -d time,13,300 -v pr pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-188411.nc temp71.nc
ncea -F -d time,1,12 -v pr pr_Amon_HadGEM2-CC_historical_r1i1p1_188412-190911.nc temp72.nc
ncrcat temp71.nc temp72.nc pr_Amon_HadGEM2-CC_historical_r1i1p1_186012-188511.nc

ncea -F -d time,13,300 -v rsds rsds_Amon_HadGEM2-CC_historical_r1i1p1_185912-188411.nc temp81.nc
ncea -F -d time,1,12 -v rsds rsds_Amon_HadGEM2-CC_historical_r1i1p1_188412-190911.nc temp82.nc
ncrcat temp81.nc temp82.nc rsds_Amon_HadGEM2-CC_historical_r1i1p1_186012-188511.nc

ncea -F -d time,13,300 -v tas tas_Amon_HadGEM2-CC_historical_r1i1p1_185912-188411.nc temp91.nc
ncea -F -d time,1,12 -v tas tas_Amon_HadGEM2-CC_historical_r1i1p1_188412-190911.nc temp92.nc
ncrcat temp91.nc temp92.nc tas_Amon_HadGEM2-CC_historical_r1i1p1_186012-188511.nc

mv *186012-188511.nc ./1861-1885/
rm temp11.nc temp12.nc temp31.nc temp41.nc temp51.nc temp61.nc temp71.nc temp81.nc temp91.nc
rm temp21.nc temp22.nc temp32.nc temp42.nc temp52.nc temp62.nc temp72.nc temp82.nc temp92.nc


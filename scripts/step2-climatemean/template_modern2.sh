#!/bin/bash

#### Only for HadGEM2 models

module load nco

cd ../../CMIP5/modelnm/1982-2005

for i in {1..24}
  do
    let "index = $i+1981"
    let "a = $i-1"
    let "start = $a * 12 + 1"
    end=$(($start+11))

    ncra -F -d time,$start,$end -v cLeaf cLeaf_Lmon_modelnm_historical_r1i1p1_198112-200511.nc  cLeaf_Lmon_modelnm_historical_r1i1p1_$index.nc
    ncra -F -d time,$start,$end -v cRoot cRoot_Lmon_modelnm_historical_r1i1p1_198112-200511.nc  cRoot_Lmon_modelnm_historical_r1i1p1_$index.nc
    ncra -F -d time,$start,$end -v cVeg cVeg_Lmon_modelnm_historical_r1i1p1_198112-200511.nc  cVeg_Lmon_modelnm_historical_r1i1p1_$index.nc
    ncra -F -d time,$start,$end -v cWood cWood_Lmon_modelnm_historical_r1i1p1_198112-200511.nc  cWood_Lmon_modelnm_historical_r1i1p1_$index.nc
    ncra -F -d time,$start,$end -v landCoverFrac landCoverFrac_Lmon_modelnm_historical_r1i1p1_198112-200511.nc  landCoverFrac_Lmon_modelnm_historical_r1i1p1_$index.nc
    ncra -F -d time,$start,$end -v treeFrac treeFrac_Lmon_modelnm_historical_r1i1p1_198112-200511.nc  treeFrac_Lmon_modelnm_historical_r1i1p1_$index.nc
    ncra -F -d time,$start,$end -v pr pr_Amon_modelnm_historical_r1i1p1_198112-200511.nc  pr_Amon_modelnm_historical_r1i1p1_$index.nc
    ncra -F -d time,$start,$end -v tas tas_Amon_modelnm_historical_r1i1p1_198112-200511.nc  tas_Amon_modelnm_historical_r1i1p1_$index.nc
    ncra -F -d time,$start,$end -v rsds rsds_Amon_modelnm_historical_r1i1p1_198112-200511.nc  rsds_Amon_modelnm_historical_r1i1p1_$index.nc

done

mv *_198112-200511.nc ../

ncrcat -v cLeaf cLeaf_Lmon_modelnm_historical_r1i1p1_198[23456789].nc cLeaf_Lmon_modelnm_historical_r1i1p1_199[0123456789].nc cLeaf_Lmon_modelnm_historical_r1i1p1_200[012345].nc test_cLeaf.nc
ncrcat -v cRoot cRoot_Lmon_modelnm_historical_r1i1p1_198[23456789].nc cRoot_Lmon_modelnm_historical_r1i1p1_199[0123456789].nc cRoot_Lmon_modelnm_historical_r1i1p1_200[012345].nc test_cRoot.nc
ncrcat -v cVeg cVeg_Lmon_modelnm_historical_r1i1p1_198[23456789].nc cVeg_Lmon_modelnm_historical_r1i1p1_199[0123456789].nc cVeg_Lmon_modelnm_historical_r1i1p1_200[012345].nc test_cVeg.nc
ncrcat -v cWood cWood_Lmon_modelnm_historical_r1i1p1_198[23456789].nc cWood_Lmon_modelnm_historical_r1i1p1_199[0123456789].nc cWood_Lmon_modelnm_historical_r1i1p1_200[012345].nc test_cWood.nc
ncrcat -v landCoverFrac landCoverFrac_Lmon_modelnm_historical_r1i1p1_198[23456789].nc landCoverFrac_Lmon_modelnm_historical_r1i1p1_199[0123456789].nc landCoverFrac_Lmon_modelnm_historical_r1i1p1_200[012345].nc test_landCoverFrac.nc
ncrcat -v treeFrac treeFrac_Lmon_modelnm_historical_r1i1p1_198[23456789].nc treeFrac_Lmon_modelnm_historical_r1i1p1_199[0123456789].nc treeFrac_Lmon_modelnm_historical_r1i1p1_200[012345].nc test_treeFrac.nc
ncrcat -v pr pr_Amon_modelnm_historical_r1i1p1_198[23456789].nc pr_Amon_modelnm_historical_r1i1p1_199[0123456789].nc pr_Amon_modelnm_historical_r1i1p1_200[012345].nc test_pr.nc
ncrcat -v rsds rsds_Amon_modelnm_historical_r1i1p1_198[23456789].nc rsds_Amon_modelnm_historical_r1i1p1_199[0123456789].nc rsds_Amon_modelnm_historical_r1i1p1_200[012345].nc test_rsds.nc
ncrcat -v tas tas_Amon_modelnm_historical_r1i1p1_198[23456789].nc tas_Amon_modelnm_historical_r1i1p1_199[0123456789].nc tas_Amon_modelnm_historical_r1i1p1_200[012345].nc test_tas.nc

ncks -A test_cRoot.nc test_cVeg.nc
ncks -A test_cVeg.nc test_cWood.nc
ncks -A test_cWood.nc test_landCoverFrac.nc
ncks -A test_landCoverFrac.nc test_treeFrac.nc
ncks -A test_treeFrac.nc test_pr.nc
ncks -A test_pr.nc test_rsds.nc
ncks -A test_rsds.nc test_tas.nc
ncks -A test_tas.nc test_cLeaf.nc

wait

mv test_cLeaf.nc ../../annual_mean/modern/modelnm.nc

rm -f *Lmon*.nc *Amon*.nc test_*.nc
mv ../*_198112-200511.nc .


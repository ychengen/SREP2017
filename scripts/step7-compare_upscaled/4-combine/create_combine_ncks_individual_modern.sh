##############################################
#!/bin/bash
#
# Usage: 1. Change Modern or PI cases
#        2. Substitute fraction number from 0
#            to 0.9 (note 0, not 0.0) 
##############################################
module load nco
set echo

i="0"

### For case frac0
cd ../../../
if [ ! -d postanalysis ]; then
   mkdir postanalysis
fi
cd postanalysis
ln -s ../Observation/GSWP3/regrid_for_CMIP5/climate_mean/modern/masked/individual/frac0/*.nc .
ln -s ../CMIP5/climate_mean/modern/masked/individual/frac0/mask_climate_*.nc .
mkdir -p input_data/modern/individual/frac0

# BNU-ESM   
ncks -A mask_climate_BNU-ESM_cWood.nc mask_climate_BNU-ESM_cLeaf.nc
ncks -A mask_climate_BNU-ESM_cLeaf.nc mask_climate_BNU-ESM_Bbio.nc
ncks -A mask_climate_BNU-ESM_Bbio.nc mask_climate_BNU-ESM_Tbio.nc
ncks -A mask_climate_BNU-ESM_Tbio.nc mask_met_climate_BNU-ESM_needleleaf.nc
ncks -A mask_met_climate_BNU-ESM_needleleaf.nc mask_met_climate_BNU-ESM_evergreen.nc
ncks -A mask_met_climate_BNU-ESM_evergreen.nc mask_met_climate_BNU-ESM_deciduous.nc
ncks -A mask_met_climate_BNU-ESM_deciduous.nc mask_met_climate_BNU-ESM_broadleaf.nc
   
mv mask_met_climate_BNU-ESM_broadleaf.nc mask_climate_all_BNU-ESM.nc
rm -f mask_climate_BNU-ESM_*.nc mask_met_climate_BNU-ESM_*.nc
   
# HadGEM2-CC   
ncks -A mask_climate_HadGEM2-CC_cWood.nc mask_climate_HadGEM2-CC_cLeaf.nc
ncks -A mask_climate_HadGEM2-CC_cLeaf.nc mask_climate_HadGEM2-CC_Bbio.nc
ncks -A mask_climate_HadGEM2-CC_Bbio.nc mask_climate_HadGEM2-CC_Tbio.nc
ncks -A mask_climate_HadGEM2-CC_Tbio.nc mask_met_climate_HadGEM2-CC_needleleaf.nc
ncks -A mask_met_climate_HadGEM2-CC_needleleaf.nc mask_met_climate_HadGEM2-CC_evergreen.nc
ncks -A mask_met_climate_HadGEM2-CC_evergreen.nc mask_met_climate_HadGEM2-CC_deciduous.nc
ncks -A mask_met_climate_HadGEM2-CC_deciduous.nc mask_met_climate_HadGEM2-CC_broadleaf.nc

mv mask_met_climate_HadGEM2-CC_broadleaf.nc mask_climate_all_HadGEM2-CC.nc
rm -f mask_climate_HadGEM2-CC_*.nc mask_met_climate_HadGEM2-CC_*.nc
   
# HadGEM2-ES   
ncks -A mask_climate_HadGEM2-ES_cWood.nc mask_climate_HadGEM2-ES_cLeaf.nc
ncks -A mask_climate_HadGEM2-ES_cLeaf.nc mask_climate_HadGEM2-ES_Bbio.nc
ncks -A mask_climate_HadGEM2-ES_Bbio.nc mask_climate_HadGEM2-ES_Tbio.nc
ncks -A mask_climate_HadGEM2-ES_Tbio.nc mask_met_climate_HadGEM2-ES_needleleaf.nc
ncks -A mask_met_climate_HadGEM2-ES_needleleaf.nc mask_met_climate_HadGEM2-ES_evergreen.nc
ncks -A mask_met_climate_HadGEM2-ES_evergreen.nc mask_met_climate_HadGEM2-ES_deciduous.nc
ncks -A mask_met_climate_HadGEM2-ES_deciduous.nc mask_met_climate_HadGEM2-ES_broadleaf.nc

mv mask_met_climate_HadGEM2-ES_broadleaf.nc mask_climate_all_HadGEM2-ES.nc
rm -f mask_climate_HadGEM2-ES_*.nc mask_met_climate_HadGEM2-ES_*.nc
   
# IPSL-CM5A-LR   
ncks -A mask_climate_IPSL-CM5A-LR_cWood.nc mask_climate_IPSL-CM5A-LR_cLeaf.nc
ncks -A mask_climate_IPSL-CM5A-LR_cLeaf.nc mask_climate_IPSL-CM5A-LR_Bbio.nc
ncks -A mask_climate_IPSL-CM5A-LR_Bbio.nc mask_climate_IPSL-CM5A-LR_Tbio.nc
ncks -A mask_climate_IPSL-CM5A-LR_Tbio.nc mask_met_climate_IPSL-CM5A-LR_needleleaf.nc
ncks -A mask_met_climate_IPSL-CM5A-LR_needleleaf.nc mask_met_climate_IPSL-CM5A-LR_evergreen.nc
ncks -A mask_met_climate_IPSL-CM5A-LR_evergreen.nc mask_met_climate_IPSL-CM5A-LR_deciduous.nc
ncks -A mask_met_climate_IPSL-CM5A-LR_deciduous.nc mask_met_climate_IPSL-CM5A-LR_broadleaf.nc
   
mv mask_met_climate_IPSL-CM5A-LR_broadleaf.nc mask_climate_all_IPSL-CM5A-LR.nc
rm -f mask_climate_IPSL-CM5A-LR_*.nc mask_met_climate_IPSL-CM5A-LR_*.nc
   
# IPSL-CM5A-MR   
ncks -A mask_climate_IPSL-CM5A-MR_cWood.nc mask_climate_IPSL-CM5A-MR_cLeaf.nc
ncks -A mask_climate_IPSL-CM5A-MR_cLeaf.nc mask_climate_IPSL-CM5A-MR_Bbio.nc
ncks -A mask_climate_IPSL-CM5A-MR_Bbio.nc mask_climate_IPSL-CM5A-MR_Tbio.nc
ncks -A mask_climate_IPSL-CM5A-MR_Tbio.nc mask_met_climate_IPSL-CM5A-MR_needleleaf.nc
ncks -A mask_met_climate_IPSL-CM5A-MR_needleleaf.nc mask_met_climate_IPSL-CM5A-MR_evergreen.nc
ncks -A mask_met_climate_IPSL-CM5A-MR_evergreen.nc mask_met_climate_IPSL-CM5A-MR_deciduous.nc
ncks -A mask_met_climate_IPSL-CM5A-MR_deciduous.nc mask_met_climate_IPSL-CM5A-MR_broadleaf.nc

mv mask_met_climate_IPSL-CM5A-MR_broadleaf.nc mask_climate_all_IPSL-CM5A-MR.nc
rm -f mask_climate_IPSL-CM5A-MR_*.nc mask_met_climate_IPSL-CM5A-MR_*.nc
   
# IPSL-CM5B-LR   
ncks -A mask_climate_IPSL-CM5B-LR_cWood.nc mask_climate_IPSL-CM5B-LR_cLeaf.nc
ncks -A mask_climate_IPSL-CM5B-LR_cLeaf.nc mask_climate_IPSL-CM5B-LR_Bbio.nc
ncks -A mask_climate_IPSL-CM5B-LR_Bbio.nc mask_climate_IPSL-CM5B-LR_Tbio.nc
ncks -A mask_climate_IPSL-CM5B-LR_Tbio.nc mask_met_climate_IPSL-CM5B-LR_needleleaf.nc
ncks -A mask_met_climate_IPSL-CM5B-LR_needleleaf.nc mask_met_climate_IPSL-CM5B-LR_evergreen.nc
ncks -A mask_met_climate_IPSL-CM5B-LR_evergreen.nc mask_met_climate_IPSL-CM5B-LR_deciduous.nc
ncks -A mask_met_climate_IPSL-CM5B-LR_deciduous.nc mask_met_climate_IPSL-CM5B-LR_broadleaf.nc
   
mv mask_met_climate_IPSL-CM5B-LR_broadleaf.nc mask_climate_all_IPSL-CM5B-LR.nc
rm -f mask_climate_IPSL-CM5B-LR_*.nc mask_met_climate_IPSL-CM5B-LR_*.nc
   
# MIROC-ESM-CHEM   
ncks -A mask_climate_MIROC-ESM-CHEM_cWood.nc mask_climate_MIROC-ESM-CHEM_cLeaf.nc
ncks -A mask_climate_MIROC-ESM-CHEM_cLeaf.nc mask_climate_MIROC-ESM-CHEM_Bbio.nc
ncks -A mask_climate_MIROC-ESM-CHEM_Bbio.nc mask_climate_MIROC-ESM-CHEM_Tbio.nc
ncks -A mask_climate_MIROC-ESM-CHEM_Tbio.nc mask_met_climate_MIROC-ESM-CHEM_needleleaf.nc
ncks -A mask_met_climate_MIROC-ESM-CHEM_needleleaf.nc mask_met_climate_MIROC-ESM-CHEM_evergreen.nc
ncks -A mask_met_climate_MIROC-ESM-CHEM_evergreen.nc mask_met_climate_MIROC-ESM-CHEM_deciduous.nc
ncks -A mask_met_climate_MIROC-ESM-CHEM_deciduous.nc mask_met_climate_MIROC-ESM-CHEM_broadleaf.nc
   
mv mask_met_climate_MIROC-ESM-CHEM_broadleaf.nc mask_climate_all_MIROC-ESM-CHEM.nc
rm -f mask_climate_MIROC-ESM-CHEM_*.nc mask_met_climate_MIROC-ESM-CHEM_*.nc
   
# MIROC-ESM   
ncks -A mask_climate_MIROC-ESM_cWood.nc mask_climate_MIROC-ESM_cLeaf.nc
ncks -A mask_climate_MIROC-ESM_cLeaf.nc mask_climate_MIROC-ESM_Bbio.nc
ncks -A mask_climate_MIROC-ESM_Bbio.nc mask_climate_MIROC-ESM_Tbio.nc
ncks -A mask_climate_MIROC-ESM_Tbio.nc mask_met_climate_MIROC-ESM_needleleaf.nc
ncks -A mask_met_climate_MIROC-ESM_needleleaf.nc mask_met_climate_MIROC-ESM_evergreen.nc
ncks -A mask_met_climate_MIROC-ESM_evergreen.nc mask_met_climate_MIROC-ESM_deciduous.nc
ncks -A mask_met_climate_MIROC-ESM_deciduous.nc mask_met_climate_MIROC-ESM_broadleaf.nc
   
mv mask_met_climate_MIROC-ESM_broadleaf.nc mask_climate_all_MIROC-ESM.nc
rm -f mask_climate_MIROC-ESM_*.nc mask_met_climate_MIROC-ESM_*.nc
   
mv *.nc input_data/modern/individual/frac0/

###################################
### For cases frac0.1 - frac0.9 ###
###################################
while [ $i -lt 9 ]
do
   i=$[$i+1]

   ln -s ../Observation/GSWP3/regrid_for_CMIP5/climate_mean/modern/masked/individual/frac0.$i/*.nc .
   ln -s ../CMIP5/climate_mean/modern/masked/individual/frac0.$i/mask_climate_*.nc .
   mkdir -p input_data/modern/individual/frac0.$i

# BNU-ESM   
   ncks -A mask_climate_BNU-ESM_cWood.nc mask_climate_BNU-ESM_cLeaf.nc
   ncks -A mask_climate_BNU-ESM_cLeaf.nc mask_climate_BNU-ESM_Bbio.nc
   ncks -A mask_climate_BNU-ESM_Bbio.nc mask_climate_BNU-ESM_Tbio.nc
   ncks -A mask_climate_BNU-ESM_Tbio.nc mask_met_climate_BNU-ESM_needleleaf.nc
   ncks -A mask_met_climate_BNU-ESM_needleleaf.nc mask_met_climate_BNU-ESM_evergreen.nc
   ncks -A mask_met_climate_BNU-ESM_evergreen.nc mask_met_climate_BNU-ESM_deciduous.nc
   ncks -A mask_met_climate_BNU-ESM_deciduous.nc mask_met_climate_BNU-ESM_broadleaf.nc
   
   mv mask_met_climate_BNU-ESM_broadleaf.nc mask_climate_all_BNU-ESM.nc
   rm -f mask_climate_BNU-ESM_*.nc mask_met_climate_BNU-ESM_*.nc
   
# HadGEM2-CC   
   ncks -A mask_climate_HadGEM2-CC_cWood.nc mask_climate_HadGEM2-CC_cLeaf.nc
   ncks -A mask_climate_HadGEM2-CC_cLeaf.nc mask_climate_HadGEM2-CC_Bbio.nc
   ncks -A mask_climate_HadGEM2-CC_Bbio.nc mask_climate_HadGEM2-CC_Tbio.nc
   ncks -A mask_climate_HadGEM2-CC_Tbio.nc mask_met_climate_HadGEM2-CC_needleleaf.nc
   ncks -A mask_met_climate_HadGEM2-CC_needleleaf.nc mask_met_climate_HadGEM2-CC_evergreen.nc
   ncks -A mask_met_climate_HadGEM2-CC_evergreen.nc mask_met_climate_HadGEM2-CC_deciduous.nc
   ncks -A mask_met_climate_HadGEM2-CC_deciduous.nc mask_met_climate_HadGEM2-CC_broadleaf.nc

   mv mask_met_climate_HadGEM2-CC_broadleaf.nc mask_climate_all_HadGEM2-CC.nc
   rm -f mask_climate_HadGEM2-CC_*.nc mask_met_climate_HadGEM2-CC_*.nc
   
# HadGEM2-ES   
   ncks -A mask_climate_HadGEM2-ES_cWood.nc mask_climate_HadGEM2-ES_cLeaf.nc
   ncks -A mask_climate_HadGEM2-ES_cLeaf.nc mask_climate_HadGEM2-ES_Bbio.nc
   ncks -A mask_climate_HadGEM2-ES_Bbio.nc mask_climate_HadGEM2-ES_Tbio.nc
   ncks -A mask_climate_HadGEM2-ES_Tbio.nc mask_met_climate_HadGEM2-ES_needleleaf.nc
   ncks -A mask_met_climate_HadGEM2-ES_needleleaf.nc mask_met_climate_HadGEM2-ES_evergreen.nc
   ncks -A mask_met_climate_HadGEM2-ES_evergreen.nc mask_met_climate_HadGEM2-ES_deciduous.nc
   ncks -A mask_met_climate_HadGEM2-ES_deciduous.nc mask_met_climate_HadGEM2-ES_broadleaf.nc

   mv mask_met_climate_HadGEM2-ES_broadleaf.nc mask_climate_all_HadGEM2-ES.nc
   rm -f mask_climate_HadGEM2-ES_*.nc mask_met_climate_HadGEM2-ES_*.nc
   
# IPSL-CM5A-LR   
   ncks -A mask_climate_IPSL-CM5A-LR_cWood.nc mask_climate_IPSL-CM5A-LR_cLeaf.nc
   ncks -A mask_climate_IPSL-CM5A-LR_cLeaf.nc mask_climate_IPSL-CM5A-LR_Bbio.nc
   ncks -A mask_climate_IPSL-CM5A-LR_Bbio.nc mask_climate_IPSL-CM5A-LR_Tbio.nc
   ncks -A mask_climate_IPSL-CM5A-LR_Tbio.nc mask_met_climate_IPSL-CM5A-LR_needleleaf.nc
   ncks -A mask_met_climate_IPSL-CM5A-LR_needleleaf.nc mask_met_climate_IPSL-CM5A-LR_evergreen.nc
   ncks -A mask_met_climate_IPSL-CM5A-LR_evergreen.nc mask_met_climate_IPSL-CM5A-LR_deciduous.nc
   ncks -A mask_met_climate_IPSL-CM5A-LR_deciduous.nc mask_met_climate_IPSL-CM5A-LR_broadleaf.nc
   
   mv mask_met_climate_IPSL-CM5A-LR_broadleaf.nc mask_climate_all_IPSL-CM5A-LR.nc
   rm -f mask_climate_IPSL-CM5A-LR_*.nc mask_met_climate_IPSL-CM5A-LR_*.nc
   
# IPSL-CM5A-MR   
   ncks -A mask_climate_IPSL-CM5A-MR_cWood.nc mask_climate_IPSL-CM5A-MR_cLeaf.nc
   ncks -A mask_climate_IPSL-CM5A-MR_cLeaf.nc mask_climate_IPSL-CM5A-MR_Bbio.nc
   ncks -A mask_climate_IPSL-CM5A-MR_Bbio.nc mask_climate_IPSL-CM5A-MR_Tbio.nc
   ncks -A mask_climate_IPSL-CM5A-MR_Tbio.nc mask_met_climate_IPSL-CM5A-MR_needleleaf.nc
   ncks -A mask_met_climate_IPSL-CM5A-MR_needleleaf.nc mask_met_climate_IPSL-CM5A-MR_evergreen.nc
   ncks -A mask_met_climate_IPSL-CM5A-MR_evergreen.nc mask_met_climate_IPSL-CM5A-MR_deciduous.nc
   ncks -A mask_met_climate_IPSL-CM5A-MR_deciduous.nc mask_met_climate_IPSL-CM5A-MR_broadleaf.nc

   mv mask_met_climate_IPSL-CM5A-MR_broadleaf.nc mask_climate_all_IPSL-CM5A-MR.nc
   rm -f mask_climate_IPSL-CM5A-MR_*.nc mask_met_climate_IPSL-CM5A-MR_*.nc
   
# IPSL-CM5B-LR   
   ncks -A mask_climate_IPSL-CM5B-LR_cWood.nc mask_climate_IPSL-CM5B-LR_cLeaf.nc
   ncks -A mask_climate_IPSL-CM5B-LR_cLeaf.nc mask_climate_IPSL-CM5B-LR_Bbio.nc
   ncks -A mask_climate_IPSL-CM5B-LR_Bbio.nc mask_climate_IPSL-CM5B-LR_Tbio.nc
   ncks -A mask_climate_IPSL-CM5B-LR_Tbio.nc mask_met_climate_IPSL-CM5B-LR_needleleaf.nc
   ncks -A mask_met_climate_IPSL-CM5B-LR_needleleaf.nc mask_met_climate_IPSL-CM5B-LR_evergreen.nc
   ncks -A mask_met_climate_IPSL-CM5B-LR_evergreen.nc mask_met_climate_IPSL-CM5B-LR_deciduous.nc
   ncks -A mask_met_climate_IPSL-CM5B-LR_deciduous.nc mask_met_climate_IPSL-CM5B-LR_broadleaf.nc
   
   mv mask_met_climate_IPSL-CM5B-LR_broadleaf.nc mask_climate_all_IPSL-CM5B-LR.nc
   rm -f mask_climate_IPSL-CM5B-LR_*.nc mask_met_climate_IPSL-CM5B-LR_*.nc
   
# MIROC-ESM-CHEM   
   ncks -A mask_climate_MIROC-ESM-CHEM_cWood.nc mask_climate_MIROC-ESM-CHEM_cLeaf.nc
   ncks -A mask_climate_MIROC-ESM-CHEM_cLeaf.nc mask_climate_MIROC-ESM-CHEM_Bbio.nc
   ncks -A mask_climate_MIROC-ESM-CHEM_Bbio.nc mask_climate_MIROC-ESM-CHEM_Tbio.nc
   ncks -A mask_climate_MIROC-ESM-CHEM_Tbio.nc mask_met_climate_MIROC-ESM-CHEM_needleleaf.nc
   ncks -A mask_met_climate_MIROC-ESM-CHEM_needleleaf.nc mask_met_climate_MIROC-ESM-CHEM_evergreen.nc
   ncks -A mask_met_climate_MIROC-ESM-CHEM_evergreen.nc mask_met_climate_MIROC-ESM-CHEM_deciduous.nc
   ncks -A mask_met_climate_MIROC-ESM-CHEM_deciduous.nc mask_met_climate_MIROC-ESM-CHEM_broadleaf.nc
   
   mv mask_met_climate_MIROC-ESM-CHEM_broadleaf.nc mask_climate_all_MIROC-ESM-CHEM.nc
   rm -f mask_climate_MIROC-ESM-CHEM_*.nc mask_met_climate_MIROC-ESM-CHEM_*.nc
   
# MIROC-ESM   
   ncks -A mask_climate_MIROC-ESM_cWood.nc mask_climate_MIROC-ESM_cLeaf.nc
   ncks -A mask_climate_MIROC-ESM_cLeaf.nc mask_climate_MIROC-ESM_Bbio.nc
   ncks -A mask_climate_MIROC-ESM_Bbio.nc mask_climate_MIROC-ESM_Tbio.nc
   ncks -A mask_climate_MIROC-ESM_Tbio.nc mask_met_climate_MIROC-ESM_needleleaf.nc
   ncks -A mask_met_climate_MIROC-ESM_needleleaf.nc mask_met_climate_MIROC-ESM_evergreen.nc
   ncks -A mask_met_climate_MIROC-ESM_evergreen.nc mask_met_climate_MIROC-ESM_deciduous.nc
   ncks -A mask_met_climate_MIROC-ESM_deciduous.nc mask_met_climate_MIROC-ESM_broadleaf.nc
   
   mv mask_met_climate_MIROC-ESM_broadleaf.nc mask_climate_all_MIROC-ESM.nc
   rm -f mask_climate_MIROC-ESM_*.nc mask_met_climate_MIROC-ESM_*.nc

   mv *.nc input_data/modern/individual/frac0.$i/
done

##############################################
#!/bin/bash
#
# Usage: 1. Change Modern or PI cases
#        2. Substitute fraction number from 0
#            to 0.9 (note 0, not 0.0) 
#        3. Only for HadGEM2-ES wih detail PFT
##############################################
module load nco
set echo

i="0"

### For case frac0
cd ../../postanalysis
ln -s ../../Observation/GSWP3/regrid_for_CMIP5/climate_mean/modern/masked/individual/detail_pft/frac0/*.nc .
ln -s ../../CMIP5/climate_mean/modern/masked/individual/detail_pft/mask_climate_*.nc .
mkdir -p input_data/modern/individual/detail_pft/frac0

# HadGEM2-ES   
ncks -A mask_climate_HadGEM2-ES_cWood.nc mask_climate_HadGEM2-ES_cLeaf.nc
ncks -A mask_climate_HadGEM2-ES_cLeaf.nc mask_climate_HadGEM2-ES_Bbio.nc
ncks -A mask_climate_HadGEM2-ES_Bbio.nc  mask_climate_HadGEM2-ES_Tbio.nc
ncks -A mask_climate_HadGEM2-ES_Tbio.nc  mask_met_climate_HadGEM2-ES_lumped.nc
ncks -A mask_met_climate_HadGEM2-ES_lumped.nc mask_met_climate_HadGEM2-ES_needleleaf.nc
ncks -A mask_met_climate_HadGEM2-ES_needleleaf.nc mask_met_climate_HadGEM2-ES_evergreen.nc
ncks -A mask_met_climate_HadGEM2-ES_evergreen.nc mask_met_climate_HadGEM2-ES_deciduous.nc
ncks -A mask_met_climate_HadGEM2-ES_deciduous.nc mask_met_climate_HadGEM2-ES_broadleaf.nc

mv mask_met_climate_HadGEM2-ES_broadleaf.nc mask_climate_all_HadGEM2-ES.nc

mv mask_climate_all*.nc input_data/modern/individual/detail_pft/frac0/
rm -f mask_*.nc

###################################
### For cases frac0.1 - frac0.9 ###
###################################
while [ $i -lt 9 ]
do
   i=$[$i+1]

   ln -s ../../Observation/GSWP3/regrid_for_CMIP5/climate_mean/modern/masked/individual/detail_pft/frac0.$i/*.nc .
   ln -s ../../CMIP5/climate_mean/modern/masked/individual/.$i/l_pft/mask_climate_*.nc .
   mkdir -p input_data/modern/individual/detail_pft/frac0.$i
   
# HadGEM2-ES   
   ncks -A mask_climate_HadGEM2-ES_cWood.nc mask_climate_HadGEM2-ES_cLeaf.nc
   ncks -A mask_climate_HadGEM2-ES_cLeaf.nc mask_climate_HadGEM2-ES_Bbio.nc
   ncks -A mask_climate_HadGEM2-ES_Bbio.nc mask_climate_HadGEM2-ES_Tbio.nc
   ncks -A mask_climate_HadGEM2-ES_Tbio.nc mask_met_climate_HadGEM2-ES_lumped.nc
   ncks -A mask_met_climate_HadGEM2-ES_lumped.nc mask_met_climate_HadGEM2-ES_needleleaf.nc
   ncks -A mask_met_climate_HadGEM2-ES_needleleaf.nc mask_met_climate_HadGEM2-ES_evergreen.nc
   ncks -A mask_met_climate_HadGEM2-ES_evergreen.nc mask_met_climate_HadGEM2-ES_deciduous.nc
   ncks -A mask_met_climate_HadGEM2-ES_deciduous.nc mask_met_climate_HadGEM2-ES_broadleaf.nc

   mv mask_met_climate_HadGEM2-ES_broadleaf.nc mask_climate_all_HadGEM2-ES.nc
   
   mv mask_climate_all*.nc input_data/modern/individual/detail_pft/frac0.$i/
   rm -f mask_*.nc
done

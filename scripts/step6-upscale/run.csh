#!/bin/csh -f

module swap PE-intel PE-pgi
module load netcdf/4.3.3.1

set echo
foreach casenm ('Bbio' 'Tbio' 'cLeaf' 'cWood')
   foreach modelnm ('BNU-ESM' 'HadGEM2-CC' 'HadGEM2-ES' 'IPSL-CM5A-LR' 'IPSL-CM5A-MR' 'IPSL-CM5B-LR' 'MIROC-ESM' 'MIROC-ESM-CHEM')

      sed -s '1s/modelname/'$modelnm'/' template.Makefile > temp.Makefile
      sed -s '1s/biome/'$casenm'/' temp.Makefile > Makefile
      mv Makefile $casenm
      cd $casenm
      make
      ./regrid-BGI_GLC_${modelnm}_${casenm}.exe
      rm -f Makefile regrid-BGI_GLC_${modelnm}_${casenm}.exe
      cd ..
      rm -f temp.Makefile 
   end

   mv ${casenm}/*.nc ../../CMIP5/landcover_GLC2000_for_CMIP5/

end


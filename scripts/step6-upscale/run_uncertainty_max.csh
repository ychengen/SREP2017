#!/bin/csh -f

set echo

foreach casenm ('Bbio' 'Tbio' 'cLeaf' 'cWood')
   foreach modelnm ('BNU-ESM' 'HadGEM2-CC' 'HadGEM2-ES' 'IPSL-CM5A-LR' 'IPSL-CM5A-MR' 'IPSL-CM5B-LR' 'MIROC-ESM' 'MIROC-ESM-CHEM')
      sed -s '1s/modelname/'$modelnm'/' template_max.Makefile > temp_max.Makefile
      sed -s '1s/biome/'$casenm'/' temp_max.Makefile > Makefile
      mv Makefile $casenm/max
      cd $casenm/max
      make
      ./regrid-BGI_GLC_${modelnm}_${casenm}_max.exe
      rm -f regrid-BGI_GLC_${modelnm}_${casenm}_max.exe
      cd ../../
      rm -f temp_max.Makefile 
   end

   mv ${casenm}/max/*.nc ../../CMIP5/landcover_GLC2000_for_CMIP5/uncertainty/max/

end


#!/bin/csh -f

set echo

foreach casenm ('Bbio' 'Tbio' 'cLeaf' 'cWood')
   foreach modelnm ('BNU-ESM' 'HadGEM2-CC' 'HadGEM2-ES' 'IPSL-CM5A-LR' 'IPSL-CM5A-MR' 'IPSL-CM5B-LR' 'MIROC-ESM' 'MIROC-ESM-CHEM')
      sed -s '1s/modelname/'$modelnm'/' template_min.Makefile > temp_min.Makefile
      sed -s '1s/biome/'$casenm'/' temp_min.Makefile > Makefile
      mv Makefile $casenm/min
      cd $casenm/min
      make
      ./regrid-BGI_GLC_${modelnm}_${casenm}_min.exe
      rm -f regrid-BGI_GLC_${modelnm}_${casenm}_min.exe
      cd ../../
      rm -f temp_min.Makefile 
   end

   mv ${casenm}/min/*.nc ../../CMIP5/landcover_GLC2000_for_CMIP5/uncertainty/min/

end


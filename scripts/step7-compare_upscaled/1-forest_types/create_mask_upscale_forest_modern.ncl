;==========================================================
; Plot the upscaled forest fraction between BGI and CMIP5
;  Note1: only boreal and temperate forests are used.
;  Note2: BGI data is already forest type only, no fraction
;         needs to be involved; models' forest fractions are
;         required.
;  Note3: Only keep grids with both observation and model
;          are available.
;  Note4: Modern time (1982-2005)
;==========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  casenm  = (/"Bbio","Tbio","cLeaf","cWood"/)
  orignm  = (/"cRoot","cVeg","cLeaf","cWood"/)
  modelnm = (/"BNU-ESM","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR",\
              "IPSL-CM5B-LR","MIROC-ESM-CHEM","MIROC-ESM"/)

;;; For forest grid definition, by determining the forest fraction (0.0 to 1.0)
  forest_threshold = (/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)

;----------------------------------------------------------------------
;;; Start looping

do k=0,dimsizes(forest_threshold)-1
do i=0,dimsizes(casenm)-1
do j=0,dimsizes(modelnm)-1
  print (casenm(i)+", "+modelnm(j)+", frac "+forest_threshold(k))

;;; Observation biomes only
  indir11= "../../CMIP5/landcover_GLC2000_for_CMIP5/"
  fname11 = "landcover_BGI_"+casenm(i)+"_"+modelnm(j)+".nc"

;;; Observation met only
  indir12= "../../Observation/GSWP3/regrid_for_CMIP5/climate_mean/modern/"
  fname12 = "climate_"+modelnm(j)+".nc"

;;; CMIP5 models all variables
  indir2 = "../../CMIP5/climate_mean/modern/"
  fname2  = "climate_"+modelnm(j)

;;; Output path
  outdir = "../../CMIP5/climate_mean/modern/masked/frac"+forest_threshold(k)+"/"

  if (isfilepresent(outdir+"mask_"+fname2+"_"+casenm(i)+".nc"))
    print ("File exist. Delete it.")
    system("/bin/rm -f "+outdir+"mask_"+fname2+"_"+casenm(i)+".nc")
  end if
  outfile = addfile(outdir+"mask_"+fname2+"_"+casenm(i)+".nc","c")

;;; Observation biomes only
  infile11= addfile (indir11+fname11,"r")

;;; Observation met only
  infile12= addfile (indir12+fname12,"r")

;;; CMIP5 models all variables
  infile2 = addfile (indir2+fname2+".nc","r")

;;; Get variables
  lat = infile11->lat
  lon = infile11->lon

;;; observation
  temp1_forestFrac  = infile11->veg_frac
  temp1_total_biome = infile11->veg_total
  cell_area         = infile11->cell_area

  obs_pr   = infile12->pr
  obs_rsds = infile12->rsds
  obs_tas  = infile12->tas

;;; model
  temp2_forestFrac  = infile2->landCoverFrac
  temp2_biome_dens  = infile2->$orignm(i)$

  model_pr   = infile2->pr
  model_rsds = infile2->rsds
  model_tas  = infile2->tas

;;; Mask model output to the range 30N-80N
  do zzz = 0,dimsizes(lat)-1
    if (lat(zzz).lt.30.0 .or. lat(zzz).gt.80.0) then
      temp2_forestFrac(:,zzz,:) = 9.96921e36
      temp2_biome_dens(zzz,:) = 9.96921e36
    end if
  end do

;---------------------------------------------------------------------
;;; Land fraction in a grid ****
  indir0 = "../../CMIP5/cellarea/"
;  fname01 = "sftlf_fx_"+modelnm(j)+"_historical_r0i0p0.nc"
  fname02 = "areacella_fx_"+modelnm(j)+"_historical_r0i0p0.nc"
;  infile01= addfile (indir0+fname01,"r")
  infile02= addfile (indir0+fname02,"r")
;  sftlf  = infile01->sftlf      ;; %
  areacella  = infile02->areacella

;---------------------------------------------------------------------
  temp1_forestFrac&lon = lon
  temp1_forestFrac@_FillValue  = 9.96921e36
  temp1_total_biome&lon = lon
  temp1_total_biome@_FillValue = 9.96921e36
  temp2_forestFrac&lon = lon
  temp2_forestFrac@_FillValue  = 9.96921e36
  temp2_biome_dens&lon = lon
  temp2_biome_dens@_FillValue  = 9.96921e36

  temp_obs_forestFrac   = dim_sum_n( temp1_forestFrac(0:9,:,:),0 )
  obs_forestFrac = mask(temp_obs_forestFrac, temp_obs_forestFrac.le. \
                        forest_threshold(k), False)
  obs_forestFrac!0 = "lat" 
  obs_forestFrac!1 = "lon" 
  obs_forestFrac&lat = lat
  obs_forestFrac&lon = lon
  temp_obs_total_biome  = dim_sum_n( temp1_total_biome(0:9,:,:),0 )
  obs_total_biome= mask(temp_obs_total_biome, obs_forestFrac.le. \
                        forest_threshold(k), False)
  obs_total_biome!0 = "lat"
  obs_total_biome!1 = "lon"
  obs_total_biome&lat = lat
  obs_total_biome&lon = lon

  if (j.eq.0) then
    model_forestFrac = temp2_forestFrac(0,:,:) + temp2_forestFrac(1,:,:) + \
                       temp2_forestFrac(2,:,:) + temp2_forestFrac(4,:,:) + \
                       temp2_forestFrac(6,:,:) + temp2_forestFrac(7,:,:)
  else if (j.eq.1 .or. j.eq.2) then
    model_forestFrac = temp2_forestFrac(0,:,:) + temp2_forestFrac(1,:,:)
  else if (j.eq.6 .or. j.eq.7) then
    model_forestFrac = temp2_forestFrac(1,:,:) + temp2_forestFrac(2,:,:) + \
                       temp2_forestFrac(3,:,:) + temp2_forestFrac(4,:,:)
  else
    model_forestFrac = temp2_forestFrac(3,:,:) + temp2_forestFrac(4,:,:) + \
                       temp2_forestFrac(5,:,:) + temp2_forestFrac(6,:,:) + \
                       temp2_forestFrac(7,:,:) + temp2_forestFrac(8,:,:)
  end if
  end if
  end if

  temp_model_forestFrac = 0.01*model_forestFrac
  model_forestFrac = mask (temp_model_forestFrac, temp_model_forestFrac.le. \
                           forest_threshold(k), False)
  model_forestFrac!0 = "lat"
  model_forestFrac!1 = "lon"
  model_forestFrac&lat = lat
  model_forestFrac&lon = lon

;  model_total_biome_forest= temp2_biome_dens*cell_area*model_forestFrac
  model_total_biome_forest= temp2_biome_dens*areacella*model_forestFrac
  model_total_biome_forest!1 = "lon"
  model_total_biome_forest&lon = lon

;;; Masking (common grids, both OBS and MODEL are non-missing values)

;; Forest fraction
  mask_obs_forestFrac = \
    mask( obs_forestFrac, ismissing(obs_forestFrac) .or. \
          ismissing(model_forestFrac) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)
  mask_model_forestFrac = \
    mask( model_forestFrac, ismissing(obs_forestFrac) .or. \
          ismissing(model_forestFrac) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)

  mask_obs_forestFrac!0 = "lat" 
  mask_obs_forestFrac!1 = "lon" 
  mask_obs_forestFrac&lat = lat
  mask_obs_forestFrac&lon = lon
  mask_model_forestFrac!0 = "lat" 
  mask_model_forestFrac!1 = "lon" 
  mask_model_forestFrac&lat = lat
  mask_model_forestFrac&lon = lon

;; Forest density
  mask_obs_biome_dens = mask(obs_total_biome/cell_area, \
    ismissing(mask_obs_forestFrac) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)
;  mask_model_forestDen = mask(model_total_biome_forest/cell_area, \
  mask_model_forestDen = mask(model_total_biome_forest/areacella, \
    ismissing(mask_model_forestFrac) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)

  mask_model_forestDen@_FillValue = temp2_biome_dens@_FillValue

;; Forest biomass
  mask_obs_total_biome = mask(obs_total_biome, \
    ismissing(mask_obs_forestFrac) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)
  mask_model_total_biome_forest = mask(model_total_biome_forest, \
    ismissing(mask_model_forestFrac) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)

  copyatt(mask_obs_total_biome, obs_total_biome)
  copyatt(mask_model_total_biome_forest, model_total_biome_forest)

;; Differences
  diff_forestFrac = mask_model_forestFrac - mask_obs_forestFrac
  diff_forestDen  = mask_model_forestDen - mask_obs_biome_dens
  diff_total_biome_forest = mask_model_total_biome_forest - mask_obs_total_biome

;-----------------------------------------------------------------------
;;; Assign attributes
  mask_obs_forestFrac!0 = "lat"
  mask_obs_forestFrac!1 = "lon"
  mask_obs_forestFrac&lat = lat
  mask_obs_forestFrac&lon = lon
  mask_obs_forestFrac@comment = "observation"
  mask_obs_forestFrac@unit = ""
  mask_obs_forestFrac@long_name = "forest fraction"

  mask_model_forestFrac!0 = "lat"
  mask_model_forestFrac!1 = "lon"
  mask_model_forestFrac&lat = lat
  mask_model_forestFrac&lon = lon
  mask_model_forestFrac@comment = "model"
  mask_model_forestFrac@unit= ""
  mask_model_forestFrac@long_name = "Forest fraction"

  diff_forestFrac!0 = "lat"
  diff_forestFrac!1 = "lon"
  diff_forestFrac&lat = lat
  diff_forestFrac&lon = lon
  diff_forestFrac@comment = "model - obs"
  diff_forestFrac@unit= ""
  diff_forestFrac@long_name="Forest fraction"

  mask_obs_total_biome!0 = "lat"
  mask_obs_total_biome!1 = "lon"
  mask_obs_total_biome&lat = lat
  mask_obs_total_biome&lon = lon
  mask_obs_total_biome@comment = "observation"
  mask_obs_total_biome@unit = "kg"
  mask_obs_total_biome@long_name = "Forest biomass"

  mask_model_total_biome_forest!0 = "lat"
  mask_model_total_biome_forest!1 = "lon"
  mask_model_total_biome_forest&lat = lat
  mask_model_total_biome_forest&lon = lon
  mask_model_total_biome_forest@comment = "model"
  mask_model_total_biome_forest@unit = "kg"
  mask_model_total_biome_forest@long_name = "Forest biomass"

  diff_total_biome_forest!0 = "lat"
  diff_total_biome_forest!1 = "lon"
  diff_total_biome_forest&lat = lat
  diff_total_biome_forest&lon = lon
  diff_total_biome_forest@comment = "model - obs"
  diff_total_biome_forest@unit = "kg"
  diff_total_biome_forest@long_name = "Forest biomass"

  mask_obs_biome_dens!0 = "lat"
  mask_obs_biome_dens!1 = "lon"
  mask_obs_biome_dens&lat = lat
  mask_obs_biome_dens&lon = lon
  mask_obs_biome_dens@comment = "observation"
  mask_obs_biome_dens@unit = "kg/m2"
  mask_obs_biome_dens@long_name = "Forest biomass density"

  mask_model_forestDen!0 = "lat"
  mask_model_forestDen!1 = "lon"
  mask_model_forestDen&lat = lat
  mask_model_forestDen&lon = lon
  mask_model_forestDen@comment = "model"
  mask_model_forestDen@unit = "kg/m2"
  mask_model_forestDen@long_name = "Forest biomass density"

  diff_forestDen!0 = "lat"
  diff_forestDen!1 = "lon"
  diff_forestDen&lat = lat
  diff_forestDen&lon = lon
  diff_forestDen@comment = "model - obs"
  diff_forestDen@unit = "kg/m2"
  diff_forestDen@long_name = "Forest biomass density"

;-----------------------------------------------------------------------
;;; Write out data

  newvar01 = "mask_obs_forestFrac_"+casenm(i)
  newvar02 = "mask_obs_total_biome_"+casenm(i)
  newvar03 = "mask_obs_biome_dens_"+casenm(i)
  newvar04 = "mask_model_forestFrac_"+casenm(i)
  newvar05 = "mask_model_total_biome_"+casenm(i)
  newvar06 = "mask_model_biome_dens_"+casenm(i)
  newvar07 = "mask_diff_forestFrac_"+casenm(i)
  newvar08 = "mask_diff_total_biome_"+casenm(i)
  newvar09 = "mask_diff_biome_dens_"+casenm(i)

  outfile->lat = lat
  outfile->lon = lon
  outfile->$newvar01$ = mask_obs_forestFrac
  outfile->$newvar02$ = mask_obs_total_biome   ;(kg)
  outfile->$newvar03$ = mask_obs_biome_dens
  outfile->$newvar04$ = mask_model_forestFrac
  outfile->$newvar05$ = mask_model_total_biome_forest   ;(kg)
  outfile->$newvar06$ = mask_model_forestDen
  outfile->$newvar07$ = diff_forestFrac
  outfile->$newvar08$ = diff_total_biome_forest   ;(kg)
  outfile->$newvar09$ = diff_forestDen

  delete(lon)
  delete(lat)
  delete(obs_pr)
  delete(obs_rsds)
  delete(obs_tas)
  delete(model_pr)
  delete(model_rsds)
  delete(model_tas)
  delete(obs_total_biome)
  delete(temp_obs_total_biome)
  delete(temp1_forestFrac)
  delete(temp1_total_biome)
  delete(cell_area)
  delete(temp2_forestFrac)
  delete(temp2_biome_dens)
  delete(obs_forestFrac)
  delete(temp_obs_forestFrac)
  delete(model_forestFrac)
  delete(temp_model_forestFrac)
  delete(model_total_biome_forest)
  delete(mask_model_forestFrac)
  delete(mask_obs_forestFrac)
  delete(mask_obs_biome_dens)
  delete(mask_model_forestDen)
  delete(mask_obs_total_biome)
  delete(mask_model_total_biome_forest)
  delete(diff_forestFrac)
  delete(diff_forestDen)
  delete(diff_total_biome_forest)
;  delete(sftlf)
  delete(areacella)
  delete(newvar01)
  delete(newvar02)
  delete(newvar03)
  delete(newvar04)
  delete(newvar05)
  delete(newvar06)
  delete(newvar07)
  delete(newvar08)
  delete(newvar09)
  delete(outfile)
  delete(infile11)
  delete(infile12)
  delete(infile2)
  delete(fname11)
  delete(fname12)
  delete(fname2)
  delete(fname01)
  delete(fname02)

end do
end do
end do

end
exit

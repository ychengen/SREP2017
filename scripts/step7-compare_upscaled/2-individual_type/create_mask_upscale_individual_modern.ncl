;==========================================================
; Plot the upscaled forest fraction between BGI and CMIP5
;  Note1: only boreal and temperate forests are used.
;  Note2: BGI data is already forest type only, no fraction
;          needs to be involved; models' forest fractions
;          are required.
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

;;; For forest grid definition, by determining the forest fraction (0.0 to 0.0)
  forest_threshold = (/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)

;----------------------------------------------------------------------
;;; Start looping

do k=0,dimsizes(forest_threshold)-1
do i=0,dimsizes(casenm)-1
do j=0,dimsizes(modelnm)-1
  print (casenm(i)+", "+modelnm(j)+", frac"+forest_threshold(k))

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
  outdir = "../../CMIP5/climate_mean/modern/masked/individual/frac"+forest_threshold(k)+"/"

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
;;; Land fraction in a grid
  indir0 = "../../CMIP5/cellarea/"
;  fname01 = "sftlf_fx_"+modelnm(j)+"_historical_r0i0p0.nc"
  fname02 = "areacella_fx_"+modelnm(j)+"_historical_r0i0p0.nc"
;  infile01= addfile (indir0+fname01,"r")
  infile02= addfile (indir0+fname02,"r")
;  sftlf  = infile01->sftlf
  areacella = infile02->areacella

;---------------------------------------------------------------------
  temp1_forestFrac&lon =lon
  temp1_forestFrac@_FillValue  = 9.96921e36
  temp1_total_biome&lon =lon
  temp1_total_biome@_FillValue = 9.96921e36
  temp2_forestFrac&lon =lon
  temp2_forestFrac@_FillValue  = 9.96921e36
  temp2_biome_dens&lon =lon
  temp2_biome_dens@_FillValue  = 9.96921e36

;;; Forest fraction
;;;;; Observation
  do c=0,21
    do p=0,dimsizes(lat)-1
      do q=0,dimsizes(lon)-1
        if (ismissing(temp1_forestFrac(c,p,q))) then
          temp1_forestFrac(c,p,q) = 0.
        end if
      end do
    end do
  end do

  temp_obs_forestFrac_evergreen  = temp1_forestFrac(0,:,:)+temp1_forestFrac(3,:,:)
  temp_obs_forestFrac_deciduous  = temp1_forestFrac(1,:,:)+temp1_forestFrac(2,:,:)  \
                                  +temp1_forestFrac(4,:,:)
  temp_obs_forestFrac_broadleaf  = temp1_forestFrac(0,:,:)+temp1_forestFrac(1,:,:)  \
                                  +temp1_forestFrac(2,:,:)
  temp_obs_forestFrac_needleleaf = temp1_forestFrac(3,:,:)+temp1_forestFrac(4,:,:)
  obs_forestFrac_evergreen = mask( temp_obs_forestFrac_evergreen,  \
                                   temp_obs_forestFrac_evergreen.le. \
                                   forest_threshold(k), False )
  obs_forestFrac_deciduous = mask( temp_obs_forestFrac_deciduous,  \
                                   temp_obs_forestFrac_deciduous.le. \
                                   forest_threshold(k), False )
  obs_forestFrac_broadleaf = mask( temp_obs_forestFrac_broadleaf,  \
                                   temp_obs_forestFrac_broadleaf.le. \
                                   forest_threshold(k), False )
  obs_forestFrac_needleleaf= mask( temp_obs_forestFrac_needleleaf,  \
                                   temp_obs_forestFrac_needleleaf.le. \
                                   forest_threshold(k), False )
  obs_forestFrac_evergreen!0 = "lat"
  obs_forestFrac_evergreen!1 = "lon"
  obs_forestFrac_evergreen&lat = lat
  obs_forestFrac_evergreen&lon = lon
  obs_forestFrac_deciduous!0 = "lat"
  obs_forestFrac_deciduous!1 = "lon"
  obs_forestFrac_deciduous&lat = lat
  obs_forestFrac_deciduous&lon = lon
  obs_forestFrac_broadleaf!0 = "lat"
  obs_forestFrac_broadleaf!1 = "lon"
  obs_forestFrac_broadleaf&lat = lat
  obs_forestFrac_broadleaf&lon = lon
  obs_forestFrac_needleleaf!0 = "lat"
  obs_forestFrac_needleleaf!1 = "lon"
  obs_forestFrac_needleleaf&lat = lat
  obs_forestFrac_needleleaf&lon = lon

;;;;; Model
  if (j.eq.0) then
    temp_model_forestFrac_evergreen  = temp2_forestFrac(0,:,:)+temp2_forestFrac(1,:,:)  \
                                      +temp2_forestFrac(4,:,:)
    temp_model_forestFrac_deciduous  = temp2_forestFrac(2,:,:)+temp2_forestFrac(6,:,:)  \
                                      +temp2_forestFrac(7,:,:)
    temp_model_forestFrac_broadleaf  = temp2_forestFrac(4,:,:)+temp2_forestFrac(6,:,:)  \
                                      +temp2_forestFrac(7,:,:)
    temp_model_forestFrac_needleleaf = temp2_forestFrac(0,:,:)+temp2_forestFrac(1,:,:)  \
                                      +temp2_forestFrac(2,:,:)

    model_forestFrac_evergreen = mask( temp_model_forestFrac_evergreen,  \
                                       temp_model_forestFrac_evergreen.le. \
                                       forest_threshold(k)*100., False )
    model_forestFrac_deciduous = mask( temp_model_forestFrac_deciduous,  \
                                       temp_model_forestFrac_deciduous.le. \
                                       forest_threshold(k)*100., False )
    model_forestFrac_broadleaf = mask( temp_model_forestFrac_broadleaf,  \
                                       temp_model_forestFrac_broadleaf.le. \
                                       forest_threshold(k)*100., False )
    model_forestFrac_needleleaf= mask( temp_model_forestFrac_needleleaf,  \
                                       temp_model_forestFrac_needleleaf.le. \
                                       forest_threshold(k)*100., False )

  else if (j.eq.1 .or. j.eq.2) then
    temp_model_forestFrac_broadleaf = temp2_forestFrac(0,:,:)
    temp_model_forestFrac_needleleaf= temp2_forestFrac(1,:,:)
    model_forestFrac_broadleaf = mask( temp_model_forestFrac_broadleaf,  \
                                       temp_model_forestFrac_broadleaf.le. \
                                       forest_threshold(k)*100., False )
    model_forestFrac_needleleaf= mask( temp_model_forestFrac_needleleaf,  \
                                       temp_model_forestFrac_needleleaf.le. \
                                       forest_threshold(k)*100., False )

;; no class for evergreen and deciduous
    model_forestFrac_evergreen = model_forestFrac_broadleaf 
    model_forestFrac_deciduous = model_forestFrac_broadleaf
    model_forestFrac_evergreen = 9.96921e36
    model_forestFrac_deciduous = 9.96921e36

  else if (j.eq.6 .or. j.eq.7) then
    temp_model_forestFrac_evergreen = temp2_forestFrac(1,:,:) + temp2_forestFrac(3,:,:)
    temp_model_forestFrac_deciduous = temp2_forestFrac(2,:,:) + temp2_forestFrac(4,:,:)
    model_forestFrac_evergreen = mask( temp_model_forestFrac_evergreen,  \
                                       temp_model_forestFrac_evergreen.le. \
                                       forest_threshold(k)*100., False )
    model_forestFrac_deciduous = mask( temp_model_forestFrac_deciduous,  \
                                       temp_model_forestFrac_deciduous.le. \
                                       forest_threshold(k)*100., False )

;; no class for broadleaf and needleleaf
    model_forestFrac_broadleaf = model_forestFrac_evergreen 
    model_forestFrac_needleleaf= model_forestFrac_evergreen
    model_forestFrac_broadleaf = 9.96921e36
    model_forestFrac_needleleaf= 9.96921e36

  else
    temp_model_forestFrac_evergreen  = temp2_forestFrac(3,:,:)+temp2_forestFrac(4,:,:) \
                                      +temp2_forestFrac(6,:,:)
    temp_model_forestFrac_deciduous  = temp2_forestFrac(5,:,:)+temp2_forestFrac(7,:,:)  \
                                      +temp2_forestFrac(8,:,:)
    temp_model_forestFrac_broadleaf  = temp2_forestFrac(4,:,:)+temp2_forestFrac(5,:,:)  \
                                      +temp2_forestFrac(7,:,:)
    temp_model_forestFrac_needleleaf = temp2_forestFrac(3,:,:)+temp2_forestFrac(6,:,:) \
                                      +temp2_forestFrac(8,:,:)

    model_forestFrac_evergreen = mask( temp_model_forestFrac_evergreen,  \
                                       temp_model_forestFrac_evergreen.le. \
                                       forest_threshold(k)*100., False )
    model_forestFrac_deciduous = mask( temp_model_forestFrac_deciduous,  \
                                       temp_model_forestFrac_deciduous.le. \
                                       forest_threshold(k)*100., False )
    model_forestFrac_broadleaf = mask( temp_model_forestFrac_broadleaf,  \
                                       temp_model_forestFrac_broadleaf.le. \
                                       forest_threshold(k)*100., False )
    model_forestFrac_needleleaf= mask( temp_model_forestFrac_needleleaf,  \
                                       temp_model_forestFrac_needleleaf.le. \
                                       forest_threshold(k)*100., False )
  end if
  end if
  end if

  model_forestFrac_evergreen!0 = "lat"
  model_forestFrac_evergreen!1 = "lon"
  model_forestFrac_evergreen&lat = lat
  model_forestFrac_evergreen&lon = lon
  model_forestFrac_deciduous!0 = "lat"
  model_forestFrac_deciduous!1 = "lon"
  model_forestFrac_deciduous&lat = lat
  model_forestFrac_deciduous&lon = lon
  model_forestFrac_broadleaf!0 = "lat"
  model_forestFrac_broadleaf!1 = "lon"
  model_forestFrac_broadleaf&lat = lat
  model_forestFrac_broadleaf&lon = lon
  model_forestFrac_needleleaf!0 = "lat"
  model_forestFrac_needleleaf!1 = "lon"
  model_forestFrac_needleleaf&lat = lat
  model_forestFrac_needleleaf&lon = lon

  temp_model_forestFrac_evergreen  = 0.01*model_forestFrac_evergreen
  temp_model_forestFrac_deciduous  = 0.01*model_forestFrac_deciduous
  temp_model_forestFrac_broadleaf  = 0.01*model_forestFrac_broadleaf
  temp_model_forestFrac_needleleaf = 0.01*model_forestFrac_needleleaf
  model_forestFrac_evergreen = mask( temp_model_forestFrac_evergreen, \
                                     temp_model_forestFrac_evergreen.le. \
                                     forest_threshold(k), False )
  model_forestFrac_deciduous = mask( temp_model_forestFrac_deciduous, \
                                     temp_model_forestFrac_deciduous.le. \
                                     forest_threshold(k), False )
  model_forestFrac_broadleaf = mask( temp_model_forestFrac_broadleaf, \
                                     temp_model_forestFrac_broadleaf.le. \
                                     forest_threshold(k), False )
  model_forestFrac_needleleaf = mask( temp_model_forestFrac_needleleaf, \
                                      temp_model_forestFrac_needleleaf.le. \
                                      forest_threshold(k), False )
  model_forestFrac_evergreen!0 = "lat"
  model_forestFrac_evergreen!1 = "lon"
  model_forestFrac_evergreen&lat = lat
  model_forestFrac_evergreen&lon = lon
  model_forestFrac_deciduous!0 = "lat"
  model_forestFrac_deciduous!1 = "lon"
  model_forestFrac_deciduous&lat = lat
  model_forestFrac_deciduous&lon = lon
  model_forestFrac_broadleaf!0 = "lat"
  model_forestFrac_broadleaf!1 = "lon"
  model_forestFrac_broadleaf&lat = lat
  model_forestFrac_broadleaf&lon = lon
  model_forestFrac_needleleaf!0 = "lat"
  model_forestFrac_needleleaf!1 = "lon"
  model_forestFrac_needleleaf&lat = lat
  model_forestFrac_needleleaf&lon = lon

;;; Forest total biomass
;;;;; Observation
  do c=0,21
    do p=0,dimsizes(lat)-1
      do q=0,dimsizes(lon)-1
        if (ismissing(temp1_total_biome(c,p,q))) then
          temp1_total_biome(c,p,q) = 0.
        end if
      end do
    end do
  end do

  temp_obs_total_biome_evergreen  = temp1_total_biome(0,:,:)+ temp1_total_biome(3,:,:)
  temp_obs_total_biome_deciduous  = temp1_total_biome(1,:,:)+ temp1_total_biome(2,:,:) \
                                   +temp1_total_biome(4,:,:)
  temp_obs_total_biome_broadleaf  = temp1_total_biome(0,:,:)+ temp1_total_biome(1,:,:) \
                                   +temp1_total_biome(2,:,:)
  temp_obs_total_biome_needleleaf = temp1_total_biome(3,:,:)+ temp1_total_biome(4,:,:)
  obs_total_biome_evergreen  = mask( temp_obs_total_biome_evergreen, \
                                     obs_forestFrac_evergreen.le.0., False)
  obs_total_biome_deciduous  = mask( temp_obs_total_biome_deciduous, \
                                     obs_forestFrac_deciduous.le.0., False)
  obs_total_biome_broadleaf  = mask( temp_obs_total_biome_broadleaf, \
                                     obs_forestFrac_broadleaf.le.0., False)
  obs_total_biome_needleleaf = mask( temp_obs_total_biome_needleleaf, \
                                     obs_forestFrac_needleleaf.le.0., False)
  obs_total_biome_evergreen!0 = "lat"
  obs_total_biome_evergreen!1 = "lon"
  obs_total_biome_evergreen&lat = lat
  obs_total_biome_evergreen&lon = lon
  obs_total_biome_deciduous!0 = "lat"
  obs_total_biome_deciduous!1 = "lon"
  obs_total_biome_deciduous&lat = lat
  obs_total_biome_deciduous&lon = lon
  obs_total_biome_broadleaf!0 = "lat"
  obs_total_biome_broadleaf!1 = "lon"
  obs_total_biome_broadleaf&lat = lat
  obs_total_biome_broadleaf&lon = lon
  obs_total_biome_needleleaf!0 = "lat"
  obs_total_biome_needleleaf!1 = "lon"
  obs_total_biome_needleleaf&lat = lat
  obs_total_biome_needleleaf&lon = lon

;;;;; Model
  model_total_biome_evergreen = temp2_biome_dens*areacella*model_forestFrac_evergreen
  model_total_biome_deciduous = temp2_biome_dens*areacella*model_forestFrac_deciduous
  model_total_biome_broadleaf = temp2_biome_dens*areacella*model_forestFrac_broadleaf
  model_total_biome_needleleaf = temp2_biome_dens*areacella*model_forestFrac_needleleaf
  model_total_biome_evergreen!0 = "lat"
  model_total_biome_evergreen!1 = "lon"
  model_total_biome_evergreen&lat = lat
  model_total_biome_evergreen&lon = lon
  model_total_biome_deciduous!0 = "lat"
  model_total_biome_deciduous!1 = "lon"
  model_total_biome_deciduous&lat = lat
  model_total_biome_deciduous&lon = lon
  model_total_biome_broadleaf!0 = "lat"
  model_total_biome_broadleaf!1 = "lon"
  model_total_biome_broadleaf&lat = lat
  model_total_biome_broadleaf&lon = lon
  model_total_biome_needleleaf!0 = "lat"
  model_total_biome_needleleaf!1 = "lon"
  model_total_biome_needleleaf&lat = lat
  model_total_biome_needleleaf&lon = lon

;;;;======================

;;; masking

;; Forest fraction
  mask_obs_forestFrac_evergreen = \
    mask( obs_forestFrac_evergreen, ismissing(obs_forestFrac_evergreen) .or. \
          ismissing(model_forestFrac_evergreen) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)
  mask_obs_forestFrac_deciduous = \
    mask( obs_forestFrac_deciduous, ismissing(obs_forestFrac_deciduous) .or. \
          ismissing(model_forestFrac_deciduous) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)
  mask_obs_forestFrac_broadleaf = \
    mask( obs_forestFrac_broadleaf, ismissing(obs_forestFrac_broadleaf) .or. \
          ismissing(model_forestFrac_broadleaf) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)
  mask_obs_forestFrac_needleleaf = \
    mask( obs_forestFrac_needleleaf, ismissing(obs_forestFrac_needleleaf) .or. \
          ismissing(model_forestFrac_needleleaf) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)

  mask_model_forestFrac_evergreen = \
    mask( model_forestFrac_evergreen, ismissing(obs_forestFrac_evergreen) .or. \
          ismissing(model_forestFrac_evergreen) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)
  mask_model_forestFrac_deciduous = \
    mask( model_forestFrac_deciduous, ismissing(obs_forestFrac_deciduous) .or. \
          ismissing(model_forestFrac_deciduous) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)
  mask_model_forestFrac_broadleaf = \
    mask( model_forestFrac_broadleaf, ismissing(obs_forestFrac_broadleaf) .or. \
          ismissing(model_forestFrac_broadleaf) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)
  mask_model_forestFrac_needleleaf = \
    mask( model_forestFrac_needleleaf, ismissing(obs_forestFrac_needleleaf) .or. \
          ismissing(model_forestFrac_needleleaf) .or. ismissing(obs_pr) .or. \
          ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
          ismissing(model_pr) .or. ismissing(model_rsds) .or. \
          ismissing(model_tas), False)

  mask_obs_forestFrac_evergreen!0 = "lat" 
  mask_obs_forestFrac_evergreen!1 = "lon" 
  mask_obs_forestFrac_evergreen&lat = lat
  mask_obs_forestFrac_evergreen&lon = lon
  mask_obs_forestFrac_deciduous!0 = "lat" 
  mask_obs_forestFrac_deciduous!1 = "lon" 
  mask_obs_forestFrac_deciduous&lat = lat
  mask_obs_forestFrac_deciduous&lon = lon
  mask_obs_forestFrac_broadleaf!0 = "lat" 
  mask_obs_forestFrac_broadleaf!1 = "lon" 
  mask_obs_forestFrac_broadleaf&lat = lat
  mask_obs_forestFrac_broadleaf&lon = lon
  mask_obs_forestFrac_needleleaf!0 = "lat" 
  mask_obs_forestFrac_needleleaf!1 = "lon" 
  mask_obs_forestFrac_needleleaf&lat = lat
  mask_obs_forestFrac_needleleaf&lon = lon
  mask_model_forestFrac_evergreen!0 = "lat" 
  mask_model_forestFrac_evergreen!1 = "lon" 
  mask_model_forestFrac_evergreen&lat = lat
  mask_model_forestFrac_evergreen&lon = lon
  mask_model_forestFrac_deciduous!0 = "lat" 
  mask_model_forestFrac_deciduous!1 = "lon" 
  mask_model_forestFrac_deciduous&lat = lat
  mask_model_forestFrac_deciduous&lon = lon
  mask_model_forestFrac_broadleaf!0 = "lat" 
  mask_model_forestFrac_broadleaf!1 = "lon" 
  mask_model_forestFrac_broadleaf&lat = lat
  mask_model_forestFrac_broadleaf&lon = lon
  mask_model_forestFrac_needleleaf!0 = "lat" 
  mask_model_forestFrac_needleleaf!1 = "lon" 
  mask_model_forestFrac_needleleaf&lat = lat
  mask_model_forestFrac_needleleaf&lon = lon

;; Forest density
  mask_obs_biome_dens_evergreen = mask(obs_total_biome_evergreen/cell_area, \
    ismissing(mask_obs_forestFrac_evergreen) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)
  mask_obs_biome_dens_deciduous = mask(obs_total_biome_deciduous/cell_area, \
    ismissing(mask_obs_forestFrac_deciduous) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)
  mask_obs_biome_dens_broadleaf = mask(obs_total_biome_broadleaf/cell_area, \
    ismissing(mask_obs_forestFrac_broadleaf) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)
  mask_obs_biome_dens_needleleaf = mask(obs_total_biome_needleleaf/cell_area, \
    ismissing(mask_obs_forestFrac_needleleaf) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)

  mask_model_forestDen_evergreen = mask(model_total_biome_evergreen/areacella, \
    ismissing(mask_model_forestFrac_evergreen) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)
  mask_model_forestDen_deciduous = mask(model_total_biome_deciduous/areacella, \
    ismissing(mask_model_forestFrac_deciduous) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)
  mask_model_forestDen_broadleaf = mask(model_total_biome_broadleaf/areacella, \
    ismissing(mask_model_forestFrac_broadleaf) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)
  mask_model_forestDen_needleleaf = mask(model_total_biome_needleleaf/areacella,\
    ismissing(mask_model_forestFrac_needleleaf) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or.\
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)

  copyatt(mask_model_forestDen_evergreen, temp2_biome_dens)
  copyatt(mask_model_forestDen_deciduous, temp2_biome_dens)
  copyatt(mask_model_forestDen_broadleaf, temp2_biome_dens)
  copyatt(mask_model_forestDen_needleleaf, temp2_biome_dens)

;; Forest biomass
  mask_obs_total_biome_evergreen = mask(obs_total_biome_evergreen, \
    ismissing(mask_obs_forestFrac_evergreen) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)
  mask_obs_total_biome_deciduous = mask(obs_total_biome_deciduous, \
    ismissing(mask_obs_forestFrac_deciduous) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)
  mask_obs_total_biome_broadleaf = mask(obs_total_biome_broadleaf, \
    ismissing(mask_obs_forestFrac_broadleaf) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)
  mask_obs_total_biome_needleleaf = mask(obs_total_biome_needleleaf, \
    ismissing(mask_obs_forestFrac_needleleaf) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas), False)

  mask_model_total_biome_evergreen = mask(model_total_biome_evergreen, \
    ismissing(mask_model_forestFrac_evergreen) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)
  mask_model_total_biome_deciduous = mask(model_total_biome_deciduous, \
    ismissing(mask_model_forestFrac_deciduous) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)
  mask_model_total_biome_broadleaf = mask(model_total_biome_broadleaf, \
    ismissing(mask_model_forestFrac_broadleaf) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)
  mask_model_total_biome_needleleaf = mask(model_total_biome_needleleaf, \
    ismissing(mask_model_forestFrac_needleleaf) .or. ismissing(obs_pr) .or. \
    ismissing(obs_rsds) .or. ismissing(obs_tas) .or. \
    ismissing(model_pr) .or. ismissing(model_rsds) .or. \
    ismissing(model_tas),False)

  copyatt(mask_obs_total_biome_evergreen, obs_total_biome_evergreen)
  copyatt(mask_obs_total_biome_deciduous, obs_total_biome_deciduous)
  copyatt(mask_obs_total_biome_broadleaf, obs_total_biome_broadleaf)
  copyatt(mask_obs_total_biome_needleleaf, obs_total_biome_needleleaf)
  copyatt(mask_model_total_biome_evergreen, model_total_biome_evergreen)
  copyatt(mask_model_total_biome_deciduous, model_total_biome_deciduous)
  copyatt(mask_model_total_biome_broadleaf, model_total_biome_broadleaf)
  copyatt(mask_model_total_biome_needleleaf, model_total_biome_needleleaf)

;; Differences

  if (j.eq.1 .or. j.eq.2) then    ;;; no specific tree types
    diff_forestFrac_evergreen = mask_model_forestFrac_evergreen - mask_obs_forestFrac_evergreen
    diff_forestFrac_deciduous = mask_model_forestFrac_deciduous - mask_obs_forestFrac_deciduous
    diff_forestDen_evergreen = mask_model_forestDen_evergreen - mask_obs_biome_dens_evergreen
    diff_forestDen_deciduous = mask_model_forestDen_deciduous - mask_obs_biome_dens_deciduous
    diff_total_biome_evergreen = mask_model_total_biome_evergreen - mask_obs_total_biome_evergreen
    diff_total_biome_deciduous = mask_model_total_biome_deciduous - mask_obs_total_biome_deciduous
    diff_forestFrac_evergreen = 9.96921e36
    diff_forestFrac_deciduous = 9.96921e36
    diff_forestDen_evergreen = 9.96921e36
    diff_forestDen_deciduous = 9.96921e36
    diff_total_biome_evergreen = 9.96921e36
    diff_total_biome_deciduous = 9.96921e36
  else
    diff_forestFrac_evergreen = mask_model_forestFrac_evergreen - mask_obs_forestFrac_evergreen
    diff_forestFrac_deciduous = mask_model_forestFrac_deciduous - mask_obs_forestFrac_deciduous
    diff_forestDen_evergreen = mask_model_forestDen_evergreen - mask_obs_biome_dens_evergreen
    diff_forestDen_deciduous = mask_model_forestDen_deciduous - mask_obs_biome_dens_deciduous
    diff_total_biome_evergreen = mask_model_total_biome_evergreen - mask_obs_total_biome_evergreen
    diff_total_biome_deciduous = mask_model_total_biome_deciduous - mask_obs_total_biome_deciduous
  end if

  if (j.eq.6 .or. j.eq.7) then   ;;; no specific tree types
    diff_forestFrac_broadleaf = mask_model_forestFrac_broadleaf - mask_obs_forestFrac_broadleaf
    diff_forestFrac_needleleaf = mask_model_forestFrac_needleleaf - mask_obs_forestFrac_needleleaf
    diff_forestDen_broadleaf = mask_model_forestDen_broadleaf - mask_obs_biome_dens_broadleaf
    diff_forestDen_needleleaf = mask_model_forestDen_needleleaf - mask_obs_biome_dens_needleleaf
    diff_total_biome_broadleaf = mask_model_total_biome_broadleaf - mask_obs_total_biome_broadleaf
    diff_total_biome_needleleaf = mask_model_total_biome_needleleaf - mask_obs_total_biome_needleleaf
    diff_forestFrac_broadleaf = 9.96921e36
    diff_forestFrac_needleleaf = 9.96921e36
    diff_forestDen_broadleaf = 9.96921e36
    diff_forestDen_needleleaf = 9.96921e36
    diff_total_biome_broadleaf = 9.96921e36
    diff_total_biome_needleleaf = 9.96921e36
  else
    diff_forestFrac_broadleaf = mask_model_forestFrac_broadleaf - mask_obs_forestFrac_broadleaf
    diff_forestFrac_needleleaf = mask_model_forestFrac_needleleaf - mask_obs_forestFrac_needleleaf
    diff_forestDen_broadleaf = mask_model_forestDen_broadleaf - mask_obs_biome_dens_broadleaf
    diff_forestDen_needleleaf = mask_model_forestDen_needleleaf - mask_obs_biome_dens_needleleaf
    diff_total_biome_broadleaf = mask_model_total_biome_broadleaf - mask_obs_total_biome_broadleaf
    diff_total_biome_needleleaf = mask_model_total_biome_needleleaf - mask_obs_total_biome_needleleaf
  end if

  copyatt(diff_forestFrac_evergreen, mask_model_forestFrac_evergreen)
  copyatt(diff_forestFrac_deciduous, mask_model_forestFrac_deciduous)
  copyatt(diff_forestDen_evergreen, mask_model_forestDen_evergreen)
  copyatt(diff_forestDen_deciduous, mask_model_forestDen_deciduous)
  copyatt(diff_total_biome_evergreen, mask_model_total_biome_evergreen)
  copyatt(diff_total_biome_deciduous, mask_model_total_biome_deciduous)

;-----------------------------------------------------------------------
;;; Assign attributes
  mask_obs_forestFrac_evergreen!0 = "lat"
  mask_obs_forestFrac_evergreen!1 = "lon"
  mask_obs_forestFrac_evergreen&lat = lat
  mask_obs_forestFrac_evergreen&lon = lon
  mask_obs_forestFrac_evergreen@comment = "observation"
  mask_obs_forestFrac_evergreen@unit = ""
  mask_obs_forestFrac_evergreen@long_name = "evergreen forest fraction"
  mask_obs_forestFrac_deciduous!0 = "lat"
  mask_obs_forestFrac_deciduous!1 = "lon"
  mask_obs_forestFrac_deciduous&lat = lat
  mask_obs_forestFrac_deciduous&lon = lon
  mask_obs_forestFrac_deciduous@comment = "observation"
  mask_obs_forestFrac_deciduous@unit = ""
  mask_obs_forestFrac_deciduous@long_name = "deciduous forest fraction"
  mask_obs_forestFrac_broadleaf!0 = "lat"
  mask_obs_forestFrac_broadleaf!1 = "lon"
  mask_obs_forestFrac_broadleaf&lat = lat
  mask_obs_forestFrac_broadleaf&lon = lon
  mask_obs_forestFrac_broadleaf@comment = "observation"
  mask_obs_forestFrac_broadleaf@unit = ""
  mask_obs_forestFrac_broadleaf@long_name = "broadleaf forest fraction"
  mask_obs_forestFrac_needleleaf!0 = "lat"
  mask_obs_forestFrac_needleleaf!1 = "lon"
  mask_obs_forestFrac_needleleaf&lat = lat
  mask_obs_forestFrac_needleleaf&lon = lon
  mask_obs_forestFrac_needleleaf@comment = "observation"
  mask_obs_forestFrac_needleleaf@unit = ""
  mask_obs_forestFrac_needleleaf@long_name = "needleleaf forest fraction"

  mask_model_forestFrac_evergreen!0 = "lat"
  mask_model_forestFrac_evergreen!1 = "lon"
  mask_model_forestFrac_evergreen&lat = lat
  mask_model_forestFrac_evergreen&lon = lon
  mask_model_forestFrac_evergreen@comment = "model"
  mask_model_forestFrac_evergreen@unit = ""
  mask_model_forestFrac_evergreen@long_name = "evergreen forest fraction"
  mask_model_forestFrac_deciduous!0 = "lat"
  mask_model_forestFrac_deciduous!1 = "lon"
  mask_model_forestFrac_deciduous&lat = lat
  mask_model_forestFrac_deciduous&lon = lon
  mask_model_forestFrac_deciduous@comment = "model"
  mask_model_forestFrac_deciduous@unit = ""
  mask_model_forestFrac_deciduous@long_name = "deciduous forest fraction"
  mask_model_forestFrac_broadleaf!0 = "lat"
  mask_model_forestFrac_broadleaf!1 = "lon"
  mask_model_forestFrac_broadleaf&lat = lat
  mask_model_forestFrac_broadleaf&lon = lon
  mask_model_forestFrac_broadleaf@comment = "model"
  mask_model_forestFrac_broadleaf@unit = ""
  mask_model_forestFrac_broadleaf@long_name = "broadleaf forest fraction"
  mask_model_forestFrac_needleleaf!0 = "lat"
  mask_model_forestFrac_needleleaf!1 = "lon"
  mask_model_forestFrac_needleleaf&lat = lat
  mask_model_forestFrac_needleleaf&lon = lon
  mask_model_forestFrac_needleleaf@comment = "model"
  mask_model_forestFrac_needleleaf@unit = ""
  mask_model_forestFrac_needleleaf@long_name = "needleleaf forest fraction"

  diff_forestFrac_evergreen!0 = "lat"
  diff_forestFrac_evergreen!1 = "lon"
  diff_forestFrac_evergreen&lat = lat
  diff_forestFrac_evergreen&lon = lon
  diff_forestFrac_evergreen@comment = "model - obs"
  diff_forestFrac_evergreen@unit= ""
  diff_forestFrac_evergreen@long_name="evergreen orest fraction"
  diff_forestFrac_deciduous!0 = "lat"
  diff_forestFrac_deciduous!1 = "lon"
  diff_forestFrac_deciduous&lat = lat
  diff_forestFrac_deciduous&lon = lon
  diff_forestFrac_deciduous@comment = "model - obs"
  diff_forestFrac_deciduous@unit= ""
  diff_forestFrac_deciduous@long_name="deciduous forest fraction"
  diff_forestFrac_broadleaf!0 = "lat"
  diff_forestFrac_broadleaf!1 = "lon"
  diff_forestFrac_broadleaf&lat = lat
  diff_forestFrac_broadleaf&lon = lon
  diff_forestFrac_broadleaf@comment = "model - obs"
  diff_forestFrac_broadleaf@unit= ""
  diff_forestFrac_broadleaf@long_name="broadleaf forest fraction"
  diff_forestFrac_needleleaf!0 = "lat"
  diff_forestFrac_needleleaf!1 = "lon"
  diff_forestFrac_needleleaf&lat = lat
  diff_forestFrac_needleleaf&lon = lon
  diff_forestFrac_needleleaf@comment = "model - obs"
  diff_forestFrac_needleleaf@unit= ""
  diff_forestFrac_needleleaf@long_name="needleleaf forest fraction"

  mask_obs_total_biome_evergreen!0 = "lat"
  mask_obs_total_biome_evergreen!1 = "lon"
  mask_obs_total_biome_evergreen&lat = lat
  mask_obs_total_biome_evergreen&lon = lon
  mask_obs_total_biome_evergreen@comment = "observation"
  mask_obs_total_biome_evergreen@unit = "kg"
  mask_obs_total_biome_evergreen@long_name = "evergreen forest biomass"
  mask_obs_total_biome_deciduous!0 = "lat"
  mask_obs_total_biome_deciduous!1 = "lon"
  mask_obs_total_biome_deciduous&lat = lat
  mask_obs_total_biome_deciduous&lon = lon
  mask_obs_total_biome_deciduous@comment = "observation"
  mask_obs_total_biome_deciduous@unit = "kg"
  mask_obs_total_biome_deciduous@long_name = "deciduous forest biomass"
  mask_obs_total_biome_broadleaf!0 = "lat"
  mask_obs_total_biome_broadleaf!1 = "lon"
  mask_obs_total_biome_broadleaf&lat = lat
  mask_obs_total_biome_broadleaf&lon = lon
  mask_obs_total_biome_broadleaf@comment = "observation"
  mask_obs_total_biome_broadleaf@unit = "kg"
  mask_obs_total_biome_broadleaf@long_name = "broadleaf forest biomass"
  mask_obs_total_biome_needleleaf!0 = "lat"
  mask_obs_total_biome_needleleaf!1 = "lon"
  mask_obs_total_biome_needleleaf&lat = lat
  mask_obs_total_biome_needleleaf&lon = lon
  mask_obs_total_biome_needleleaf@comment = "observation"
  mask_obs_total_biome_needleleaf@unit = "kg"
  mask_obs_total_biome_needleleaf@long_name = "needleleaf forest biomass"

  mask_model_total_biome_evergreen!0 = "lat"
  mask_model_total_biome_evergreen!1 = "lon"
  mask_model_total_biome_evergreen&lat = lat
  mask_model_total_biome_evergreen&lon = lon
  mask_model_total_biome_evergreen@comment = "model"
  mask_model_total_biome_evergreen@unit = "kg"
  mask_model_total_biome_evergreen@long_name = "evergreen forest biomass"
  mask_model_total_biome_deciduous!0 = "lat"
  mask_model_total_biome_deciduous!1 = "lon"
  mask_model_total_biome_deciduous&lat = lat
  mask_model_total_biome_deciduous&lon = lon
  mask_model_total_biome_deciduous@comment = "model"
  mask_model_total_biome_deciduous@unit = "kg"
  mask_model_total_biome_deciduous@long_name = "deciduous forest biomass"
  mask_model_total_biome_broadleaf!0 = "lat"
  mask_model_total_biome_broadleaf!1 = "lon"
  mask_model_total_biome_broadleaf&lat = lat
  mask_model_total_biome_broadleaf&lon = lon
  mask_model_total_biome_broadleaf@comment = "model"
  mask_model_total_biome_broadleaf@unit = "kg"
  mask_model_total_biome_broadleaf@long_name = "broadleaf forest biomass"
  mask_model_total_biome_needleleaf!0 = "lat"
  mask_model_total_biome_needleleaf!1 = "lon"
  mask_model_total_biome_needleleaf&lat = lat
  mask_model_total_biome_needleleaf&lon = lon
  mask_model_total_biome_needleleaf@comment = "model"
  mask_model_total_biome_needleleaf@unit = "kg"
  mask_model_total_biome_needleleaf@long_name = "needleleaf forest biomass"

  diff_total_biome_evergreen!0 = "lat"
  diff_total_biome_evergreen!1 = "lon"
  diff_total_biome_evergreen&lat = lat
  diff_total_biome_evergreen&lon = lon
  diff_total_biome_evergreen@comment = "model - obs"
  diff_total_biome_evergreen@unit = "kg"
  diff_total_biome_evergreen@long_name = "evergreen forest biomass"
  diff_total_biome_deciduous!0 = "lat"
  diff_total_biome_deciduous!1 = "lon"
  diff_total_biome_deciduous&lat = lat
  diff_total_biome_deciduous&lon = lon
  diff_total_biome_deciduous@comment = "model - obs"
  diff_total_biome_deciduous@unit = "kg"
  diff_total_biome_deciduous@long_name = "deciduous forest biomass"
  diff_total_biome_broadleaf!0 = "lat"
  diff_total_biome_broadleaf!1 = "lon"
  diff_total_biome_broadleaf&lat = lat
  diff_total_biome_broadleaf&lon = lon
  diff_total_biome_broadleaf@comment = "model - obs"
  diff_total_biome_broadleaf@unit = "kg"
  diff_total_biome_broadleaf@long_name = "broadleaf forest biomass"
  diff_total_biome_needleleaf!0 = "lat"
  diff_total_biome_needleleaf!1 = "lon"
  diff_total_biome_needleleaf&lat = lat
  diff_total_biome_needleleaf&lon = lon
  diff_total_biome_needleleaf@comment = "model - obs"
  diff_total_biome_needleleaf@unit = "kg"
  diff_total_biome_needleleaf@long_name = "needleleaf forest biomass"

  mask_obs_biome_dens_evergreen!0 = "lat"
  mask_obs_biome_dens_evergreen!1 = "lon"
  mask_obs_biome_dens_evergreen&lat = lat
  mask_obs_biome_dens_evergreen&lon = lon
  mask_obs_biome_dens_evergreen@comment = "observation"
  mask_obs_biome_dens_evergreen@unit = "kg/m2"
  mask_obs_biome_dens_evergreen@long_name = "evergreen forest density"
  mask_obs_biome_dens_deciduous!0 = "lat"
  mask_obs_biome_dens_deciduous!1 = "lon"
  mask_obs_biome_dens_deciduous&lat = lat
  mask_obs_biome_dens_deciduous&lon = lon
  mask_obs_biome_dens_deciduous@comment = "observation"
  mask_obs_biome_dens_deciduous@unit = "kg/m2"
  mask_obs_biome_dens_deciduous@long_name = "deciduous forest density"
  mask_obs_biome_dens_broadleaf!0 = "lat"
  mask_obs_biome_dens_broadleaf!1 = "lon"
  mask_obs_biome_dens_broadleaf&lat = lat
  mask_obs_biome_dens_broadleaf&lon = lon
  mask_obs_biome_dens_broadleaf@comment = "observation"
  mask_obs_biome_dens_broadleaf@unit = "kg/m2"
  mask_obs_biome_dens_broadleaf@long_name = "broadleaf forest density"
  mask_obs_biome_dens_needleleaf!0 = "lat"
  mask_obs_biome_dens_needleleaf!1 = "lon"
  mask_obs_biome_dens_needleleaf&lat = lat
  mask_obs_biome_dens_needleleaf&lon = lon
  mask_obs_biome_dens_needleleaf@comment = "observation"
  mask_obs_biome_dens_needleleaf@unit = "kg/m2"
  mask_obs_biome_dens_needleleaf@long_name = "needleleaf forest density"

  mask_model_forestDen_evergreen!0 = "lat"
  mask_model_forestDen_evergreen!1 = "lon"
  mask_model_forestDen_evergreen&lat = lat
  mask_model_forestDen_evergreen&lon = lon
  mask_model_forestDen_evergreen@comment = "model"
  mask_model_forestDen_evergreen@unit = "kg/m2"
  mask_model_forestDen_evergreen@long_name = "evergreen forest density"
  mask_model_forestDen_deciduous!0 = "lat"
  mask_model_forestDen_deciduous!1 = "lon"
  mask_model_forestDen_deciduous&lat = lat
  mask_model_forestDen_deciduous&lon = lon
  mask_model_forestDen_deciduous@comment = "model"
  mask_model_forestDen_deciduous@unit = "kg/m2"
  mask_model_forestDen_deciduous@long_name = "deciduous forest density"
  mask_model_forestDen_broadleaf!0 = "lat"
  mask_model_forestDen_broadleaf!1 = "lon"
  mask_model_forestDen_broadleaf&lat = lat
  mask_model_forestDen_broadleaf&lon = lon
  mask_model_forestDen_broadleaf@comment = "model"
  mask_model_forestDen_broadleaf@unit = "kg/m2"
  mask_model_forestDen_broadleaf@long_name = "broadleaf forest density"
  mask_model_forestDen_needleleaf!0 = "lat"
  mask_model_forestDen_needleleaf!1 = "lon"
  mask_model_forestDen_needleleaf&lat = lat
  mask_model_forestDen_needleleaf&lon = lon
  mask_model_forestDen_needleleaf@comment = "model"
  mask_model_forestDen_needleleaf@unit = "kg/m2"
  mask_model_forestDen_needleleaf@long_name = "needleleaf forest density"

  diff_forestDen_evergreen!0 = "lat"
  diff_forestDen_evergreen!1 = "lon"
  diff_forestDen_evergreen&lat = lat
  diff_forestDen_evergreen&lon = lon
  diff_forestDen_evergreen@comment = "model - obs"
  diff_forestDen_evergreen@unit = "kg/m2"
  diff_forestDen_evergreen@long_name = "evergreen forest density"
  diff_forestDen_deciduous!0 = "lat"
  diff_forestDen_deciduous!1 = "lon"
  diff_forestDen_deciduous&lat = lat
  diff_forestDen_deciduous&lon = lon
  diff_forestDen_deciduous@comment = "model - obs"
  diff_forestDen_deciduous@unit = "kg/m2"
  diff_forestDen_deciduous@long_name = "deciduous forest density"
  diff_forestDen_broadleaf!0 = "lat"
  diff_forestDen_broadleaf!1 = "lon"
  diff_forestDen_broadleaf&lat = lat
  diff_forestDen_broadleaf&lon = lon
  diff_forestDen_broadleaf@comment = "model - obs"
  diff_forestDen_broadleaf@unit = "kg/m2"
  diff_forestDen_broadleaf@long_name = "broadleaf forest density"
  diff_forestDen_needleleaf!0 = "lat"
  diff_forestDen_needleleaf!1 = "lon"
  diff_forestDen_needleleaf&lat = lat
  diff_forestDen_needleleaf&lon = lon
  diff_forestDen_needleleaf@comment = "model - obs"
  diff_forestDen_needleleaf@unit = "kg/m2"
  diff_forestDen_needleleaf@long_name = "needleleaf forest density"

;-----------------------------------------------------------------------
;;; Write out data

  newvar01 = "mask_obs_forestFrac_evergreen_"+casenm(i)
  newvar11 = "mask_obs_forestFrac_deciduous_"+casenm(i)
  newvar21 = "mask_obs_forestFrac_broadleaf_"+casenm(i)
  newvar31 = "mask_obs_forestFrac_needleleaf_"+casenm(i)
  newvar02 = "mask_obs_total_biome_evergreen_"+casenm(i)
  newvar12 = "mask_obs_total_biome_deciduous_"+casenm(i)
  newvar22 = "mask_obs_total_biome_broadleaf_"+casenm(i)
  newvar32 = "mask_obs_total_biome_needleleaf_"+casenm(i)
  newvar03 = "mask_obs_biome_dens_evergreen_"+casenm(i)
  newvar13 = "mask_obs_biome_dens_deciduous_"+casenm(i)
  newvar23 = "mask_obs_biome_dens_broadleaf_"+casenm(i)
  newvar33 = "mask_obs_biome_dens_needleleaf_"+casenm(i)
  newvar04 = "mask_model_forestFrac_evergreen_"+casenm(i)
  newvar14 = "mask_model_forestFrac_deciduous_"+casenm(i)
  newvar24 = "mask_model_forestFrac_broadleaf_"+casenm(i)
  newvar34 = "mask_model_forestFrac_needleleaf_"+casenm(i)
  newvar05 = "mask_model_total_biome_evergreen_"+casenm(i)
  newvar15 = "mask_model_total_biome_deciduous_"+casenm(i)
  newvar25 = "mask_model_total_biome_broadleaf_"+casenm(i)
  newvar35 = "mask_model_total_biome_needleleaf_"+casenm(i)
  newvar06 = "mask_model_biome_dens_evergreen_"+casenm(i)
  newvar16 = "mask_model_biome_dens_deciduous_"+casenm(i)
  newvar26 = "mask_model_biome_dens_broadleaf_"+casenm(i)
  newvar36 = "mask_model_biome_dens_needleleaf_"+casenm(i)
  newvar07 = "mask_diff_forestFrac_evergreen_"+casenm(i)
  newvar17 = "mask_diff_forestFrac_deciduous_"+casenm(i)
  newvar27 = "mask_diff_forestFrac_broadleaf_"+casenm(i)
  newvar37 = "mask_diff_forestFrac_needleleaf_"+casenm(i)
  newvar08 = "mask_diff_total_biome_evergreen_"+casenm(i)
  newvar18 = "mask_diff_total_biome_deciduous_"+casenm(i)
  newvar28 = "mask_diff_total_biome_broadleaf_"+casenm(i)
  newvar38 = "mask_diff_total_biome_needleleaf_"+casenm(i)
  newvar09 = "mask_diff_biome_dens_evergreen_"+casenm(i)
  newvar19 = "mask_diff_biome_dens_deciduous_"+casenm(i)
  newvar29 = "mask_diff_biome_dens_broadleaf_"+casenm(i)
  newvar39 = "mask_diff_biome_dens_needleleaf_"+casenm(i)

  outfile->lat = lat
  outfile->lon = lon
  outfile->$newvar01$ = mask_obs_forestFrac_evergreen
  outfile->$newvar11$ = mask_obs_forestFrac_deciduous
  outfile->$newvar21$ = mask_obs_forestFrac_broadleaf
  outfile->$newvar31$ = mask_obs_forestFrac_needleleaf
  outfile->$newvar02$ = mask_obs_total_biome_evergreen   ;(kg)
  outfile->$newvar12$ = mask_obs_total_biome_deciduous   ;(kg)
  outfile->$newvar22$ = mask_obs_total_biome_broadleaf   ;(kg)
  outfile->$newvar32$ = mask_obs_total_biome_needleleaf  ;(kg)
  outfile->$newvar03$ = mask_obs_biome_dens_evergreen
  outfile->$newvar13$ = mask_obs_biome_dens_deciduous
  outfile->$newvar23$ = mask_obs_biome_dens_broadleaf
  outfile->$newvar33$ = mask_obs_biome_dens_needleleaf
  outfile->$newvar04$ = mask_model_forestFrac_evergreen
  outfile->$newvar14$ = mask_model_forestFrac_deciduous
  outfile->$newvar24$ = mask_model_forestFrac_broadleaf
  outfile->$newvar34$ = mask_model_forestFrac_needleleaf
  outfile->$newvar05$ = mask_model_total_biome_evergreen   ;(kg)
  outfile->$newvar15$ = mask_model_total_biome_deciduous   ;(kg)
  outfile->$newvar25$ = mask_model_total_biome_broadleaf   ;(kg)
  outfile->$newvar35$ = mask_model_total_biome_needleleaf  ;(kg)
  outfile->$newvar06$ = mask_model_forestDen_evergreen
  outfile->$newvar16$ = mask_model_forestDen_deciduous
  outfile->$newvar26$ = mask_model_forestDen_broadleaf
  outfile->$newvar36$ = mask_model_forestDen_needleleaf
  outfile->$newvar07$ = diff_forestFrac_evergreen
  outfile->$newvar17$ = diff_forestFrac_deciduous
  outfile->$newvar27$ = diff_forestFrac_broadleaf
  outfile->$newvar37$ = diff_forestFrac_needleleaf
  outfile->$newvar08$ = diff_total_biome_evergreen   ;(kg)
  outfile->$newvar18$ = diff_total_biome_deciduous   ;(kg)
  outfile->$newvar28$ = diff_total_biome_broadleaf   ;(kg)
  outfile->$newvar38$ = diff_total_biome_needleleaf  ;(kg)
  outfile->$newvar09$ = diff_forestDen_evergreen
  outfile->$newvar19$ = diff_forestDen_deciduous
  outfile->$newvar29$ = diff_forestDen_broadleaf
  outfile->$newvar39$ = diff_forestDen_needleleaf

  delete(lon)
  delete(lat)
  delete(obs_pr)
  delete(obs_rsds)
  delete(obs_tas)
  delete(model_pr)
  delete(model_rsds)
  delete(model_tas)
  delete(obs_total_biome_evergreen)
  delete(obs_total_biome_deciduous)
  delete(obs_total_biome_broadleaf)
  delete(obs_total_biome_needleleaf)
  delete(temp_obs_total_biome_evergreen)
  delete(temp_obs_total_biome_deciduous)
  delete(temp_obs_total_biome_broadleaf)
  delete(temp_obs_total_biome_needleleaf)
  delete(temp1_forestFrac)
  delete(temp1_total_biome)
  delete(cell_area)
  delete(areacella)
;  delete(sftlf)
  delete(temp2_forestFrac)
  delete(temp2_biome_dens)
  delete(obs_forestFrac_evergreen)
  delete(obs_forestFrac_deciduous)
  delete(obs_forestFrac_broadleaf)
  delete(obs_forestFrac_needleleaf)
  delete(temp_obs_forestFrac_evergreen)
  delete(temp_obs_forestFrac_deciduous)
  delete(temp_obs_forestFrac_broadleaf)
  delete(temp_obs_forestFrac_needleleaf)
  delete(model_forestFrac_evergreen)
  delete(model_forestFrac_deciduous)
  delete(model_forestFrac_broadleaf)
  delete(model_forestFrac_needleleaf)
  delete(temp_model_forestFrac_evergreen)
  delete(temp_model_forestFrac_deciduous)
  delete(temp_model_forestFrac_broadleaf)
  delete(temp_model_forestFrac_needleleaf)
  delete(model_total_biome_evergreen)
  delete(model_total_biome_deciduous)
  delete(model_total_biome_broadleaf)
  delete(model_total_biome_needleleaf)
  delete(mask_model_forestFrac_evergreen)
  delete(mask_model_forestFrac_deciduous)
  delete(mask_model_forestFrac_broadleaf)
  delete(mask_model_forestFrac_needleleaf)
  delete(mask_obs_forestFrac_evergreen)
  delete(mask_obs_forestFrac_deciduous)
  delete(mask_obs_forestFrac_broadleaf)
  delete(mask_obs_forestFrac_needleleaf)
  delete(mask_obs_biome_dens_evergreen)
  delete(mask_obs_biome_dens_deciduous)
  delete(mask_obs_biome_dens_broadleaf)
  delete(mask_obs_biome_dens_needleleaf)
  delete(mask_model_forestDen_evergreen)
  delete(mask_model_forestDen_deciduous)
  delete(mask_model_forestDen_broadleaf)
  delete(mask_model_forestDen_needleleaf)
  delete(mask_obs_total_biome_evergreen)
  delete(mask_obs_total_biome_deciduous)
  delete(mask_obs_total_biome_broadleaf)
  delete(mask_obs_total_biome_needleleaf)
  delete(mask_model_total_biome_evergreen)
  delete(mask_model_total_biome_deciduous)
  delete(mask_model_total_biome_broadleaf)
  delete(mask_model_total_biome_needleleaf)
  delete(diff_forestFrac_evergreen)
  delete(diff_forestFrac_deciduous)
  delete(diff_forestFrac_broadleaf)
  delete(diff_forestFrac_needleleaf)
  delete(diff_forestDen_evergreen)
  delete(diff_forestDen_deciduous)
  delete(diff_forestDen_broadleaf)
  delete(diff_forestDen_needleleaf)
  delete(diff_total_biome_evergreen)
  delete(diff_total_biome_deciduous)
  delete(diff_total_biome_broadleaf)
  delete(diff_total_biome_needleleaf)
  delete(newvar01)
  delete(newvar02)
  delete(newvar03)
  delete(newvar04)
  delete(newvar05)
  delete(newvar06)
  delete(newvar07)
  delete(newvar08)
  delete(newvar09)
  delete(newvar11)
  delete(newvar12)
  delete(newvar13)
  delete(newvar14)
  delete(newvar15)
  delete(newvar16)
  delete(newvar17)
  delete(newvar18)
  delete(newvar19)
  delete(newvar21)
  delete(newvar22)
  delete(newvar23)
  delete(newvar24)
  delete(newvar25)
  delete(newvar26)
  delete(newvar27)
  delete(newvar28)
  delete(newvar29)
  delete(newvar31)
  delete(newvar32)
  delete(newvar33)
  delete(newvar34)
  delete(newvar35)
  delete(newvar36)
  delete(newvar37)
  delete(newvar38)
  delete(newvar39)
  delete(outfile)
  delete(infile11)
  delete(infile12)
  delete(infile2)
  delete(fname11)
  delete(fname12)
  delete(fname2)

end do
end do
end do

end
exit

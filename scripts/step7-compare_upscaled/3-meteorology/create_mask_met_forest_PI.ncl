;==========================================================
; Plot the upscaled forest fraction between BGI and CMIP5
;  Note1: only boreal and temperate forests are used.
;  Note2: BGI data is already forest type only, no fraction
;         needs to be involved; models' forest fractions are
;         required.
;  Note3: Pre-Industrial time (1861-1885)
;==========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

 modelnm = (/"BNU-ESM","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR",\
              "IPSL-CM5B-LR","MIROC-ESM-CHEM","MIROC-ESM"/)

;;; For forest grid definition, by determining the forest fraction (0.0 to 1.0)
  forest_threshold = (/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)

;----------------------------------------------------------------------
;;; Start looping

do k=0,dimsizes(forest_threshold)-1
do i=0,dimsizes(modelnm)-1
  print ("met "+modelnm(i))

;;; Observation met
  indir11= "../../Observation/GSWP3/regrid_for_CMIP5/climate_mean/PI/"
  fname11 = "climate_"+modelnm(i)+".nc"
  infile11= addfile (indir11+fname11,"r")

;;; CMIP5 models all variables
  indir2 = "../../CMIP5/climate_mean/PI/"
  fname2  = "climate_"+modelnm(i)
  infile2 = addfile (indir2+fname2+".nc","r")

;;; Mask file
  indir3 = indir2+"masked/frac"+forest_threshold(k)+"/"
  fname3 = "mask_climate_"+modelnm(i)+"_Tbio.nc"
  infile3 = addfile (indir3+fname3,"r")

;;; Output path
  outdir = indir11+"masked/frac"+forest_threshold(k)+"/"

  if (isfilepresent(outdir+"mask_met_"+fname2+".nc"))
    print ("File exist. Delete it.")
    system("/bin/rm -f "+outdir+"mask_met_"+fname2+".nc")
  end if

  outfile = addfile(outdir+"mask_met_"+fname2+".nc","c")

;;---------------------------------------------------------------------
;;; Get variables
  lat = infile11->lat
  lon = infile11->lon

;;; observation
  obs_pr   = infile11->pr
  obs_rsds = infile11->rsds
  obs_tas  = infile11->tas

;;; model
  temp_model_pr = infile2->pr
  model_rsds    = infile2->rsds
  model_tas     = infile2->tas

;;; mask file   ;; use forest fraction as filter
   mask_var = "mask_obs_forestFrac_Tbio"
   mask_obs_forestFrac = infile3->$mask_var$

;---------------------------------------------------------------------
  model_pr= temp_model_pr*86400.
  copyatt(model_pr,temp_model_pr)
  model_pr@_FillValue = 9.96921e36
  model_pr@missing_value = 9.96921e36
  model_pr@unit = "mm/day"

;;; masking

;; Pricipitation
  mask_obs_pr = mask( obs_pr, ismissing(mask_obs_forestFrac) .or. \
    ismissing(obs_pr) .or. ismissing(obs_rsds) .or. \
    ismissing(obs_tas) .or. ismissing(model_pr) .or. \
    ismissing(model_rsds) .or. ismissing(model_tas) ,False)
  mask_model_pr = mask( model_pr, ismissing(mask_obs_forestFrac) .or. \
    ismissing(obs_pr) .or. ismissing(obs_rsds) .or. \
    ismissing(obs_tas) .or. ismissing(model_pr) .or. \
    ismissing(model_rsds) .or. ismissing(model_tas) ,False)

  mask_obs_pr!0 = "lat" 
  mask_obs_pr!1 = "lon" 
  mask_obs_pr&lat = lat
  mask_obs_pr&lon = lon
  mask_model_pr!0 = "lat" 
  mask_model_pr!1 = "lon" 
  mask_model_pr&lat = lat
  mask_model_pr&lon = lon

;; Radiation
  mask_obs_rsds = mask( obs_rsds, ismissing(mask_obs_forestFrac) .or. \
    ismissing(obs_pr) .or. ismissing(obs_rsds) .or. \
    ismissing(obs_tas) .or. ismissing(model_pr) .or. \
    ismissing(model_rsds) .or. ismissing(model_tas) ,False)
  mask_model_rsds = mask( model_rsds, ismissing(mask_obs_forestFrac) .or. \
    ismissing(obs_pr) .or. ismissing(obs_rsds) .or. \
    ismissing(obs_tas) .or. ismissing(model_pr) .or. \
    ismissing(model_rsds) .or. ismissing(model_tas) ,False)

  mask_obs_rsds!0 = "lat" 
  mask_obs_rsds!1 = "lon" 
  mask_obs_rsds&lat = lat
  mask_obs_rsds&lon = lon
  mask_model_rsds!0 = "lat" 
  mask_model_rsds!1 = "lon" 
  mask_model_rsds&lat = lat
  mask_model_rsds&lon = lon

;; Temperature
  mask_obs_tas = mask( obs_tas, ismissing(mask_obs_forestFrac) .or. \
    ismissing(obs_pr) .or. ismissing(obs_rsds) .or. \
    ismissing(obs_tas) .or. ismissing(model_pr) .or. \
    ismissing(model_rsds) .or. ismissing(model_tas) ,False)
  mask_model_tas = mask( model_tas, ismissing(mask_obs_forestFrac) .or. \
    ismissing(obs_pr) .or. ismissing(obs_rsds) .or. \
    ismissing(obs_tas) .or. ismissing(model_pr) .or. \
    ismissing(model_rsds) .or. ismissing(model_tas) ,False)

  mask_obs_tas!0 = "lat" 
  mask_obs_tas!1 = "lon" 
  mask_obs_tas&lat = lat
  mask_obs_tas&lon = lon
  mask_model_tas!0 = "lat" 
  mask_model_tas!1 = "lon" 
  mask_model_tas&lat = lat
  mask_model_tas&lon = lon

;; Differences
  diff_pr = mask_model_pr - mask_obs_pr
  diff_rsds  = mask_model_rsds - mask_obs_rsds
  diff_tas = mask_model_tas - mask_obs_tas

;-----------------------------------------------------------------------
;;; Assign attributes
  mask_obs_pr!0 = "lat"
  mask_obs_pr!1 = "lon"
  mask_obs_pr&lat = lat
  mask_obs_pr&lon = lon
  mask_obs_pr@comment = "observation"
  mask_obs_pr@unit = "mm/day"
  mask_obs_pr@long_name = "precipitation rate"

  mask_model_pr!0 = "lat"
  mask_model_pr!1 = "lon"
  mask_model_pr&lat = lat
  mask_model_pr&lon = lon
  mask_model_pr@comment = "model"
  mask_model_pr@unit= "mm/day"
  mask_model_pr@long_name = "precipitation rate"

  diff_pr!0 = "lat"
  diff_pr!1 = "lon"
  diff_pr&lat = lat
  diff_pr&lon = lon
  diff_pr@comment = "model - obs"
  diff_pr@unit= "mm/day"
  diff_pr@long_name="precipitation rate"

  mask_obs_rsds!0 = "lat"
  mask_obs_rsds!1 = "lon"
  mask_obs_rsds&lat = lat
  mask_obs_rsds&lon = lon
  mask_obs_rsds@comment = "observation"
  mask_obs_rsds@unit = "W/m2"
  mask_obs_rsds@long_name = "downward shortwave radiation"

  mask_model_rsds!0 = "lat"
  mask_model_rsds!1 = "lon"
  mask_model_rsds&lat = lat
  mask_model_rsds&lon = lon
  mask_model_rsds@comment = "model"
  mask_model_rsds@unit = "W/m2"
  mask_model_rsds@long_name = "downward shortwave radiation"

  diff_rsds!0 = "lat"
  diff_rsds!1 = "lon"
  diff_rsds&lat = lat
  diff_rsds&lon = lon
  diff_rsds@comment = "model - obs"
  diff_rsds@unit = "W/m2"
  diff_rsds@long_name = "downward shortwave radiation"

  mask_obs_tas!0 = "lat"
  mask_obs_tas!1 = "lon"
  mask_obs_tas&lat = lat
  mask_obs_tas&lon = lon
  mask_obs_tas@comment = "observation"
  mask_obs_tas@unit = "K"
  mask_obs_tas@long_name = "sufrace temperature"

  mask_model_tas!0 = "lat"
  mask_model_tas!1 = "lon"
  mask_model_tas&lat = lat
  mask_model_tas&lon = lon
  mask_model_tas@comment = "model"
  mask_model_tas@unit = "K"
  mask_model_tas@long_name = "surface temperature"

  diff_tas!0 = "lat"
  diff_tas!1 = "lon"
  diff_tas&lat = lat
  diff_tas&lon = lon
  diff_tas@comment = "model - obs"
  diff_tas@unit = "K"
  diff_tas@long_name = "surface temperature"

;-----------------------------------------------------------------------
;;; Write out data
  outfile->lat = lat
  outfile->lon = lon
  outfile->mask_obs_pr    = mask_obs_pr
  outfile->mask_obs_rsds  = mask_obs_rsds
  outfile->mask_obs_tas   = mask_obs_tas
  outfile->mask_model_pr  = mask_model_pr
  outfile->mask_model_rsds= mask_model_rsds
  outfile->mask_model_tas = mask_model_tas
  outfile->diff_pr        = diff_pr
  outfile->diff_rsds      = diff_rsds
  outfile->diff_tas       = diff_tas

  delete(lon)
  delete(lat)
  delete(obs_pr)
  delete(obs_tas)
  delete(obs_rsds)
  delete(temp_model_pr)
  delete(model_pr)
  delete(model_tas)
  delete(model_rsds)
  delete(mask_obs_pr)
  delete(mask_obs_tas)
  delete(mask_obs_rsds)
  delete(mask_model_pr)
  delete(mask_model_tas)
  delete(mask_model_rsds)
  delete(diff_pr)
  delete(diff_tas)
  delete(diff_rsds)
  delete(mask_obs_forestFrac)
  delete(outfile)
  delete(infile11)
  delete(infile2)
  delete(infile3)
  delete(fname11)
  delete(fname2)
  delete(fname3)

end do
end do

end
exit

;===========================================================
; Plot the upscaled forest fraction between BGI and CMIP5
;  Note1: only boreal and temperate forests are used.
;  Note2: BGI data is already forest type only, no fraction
;         needs to be involved; models' forest fractions are
;         required.
;  Note3: for PI time (1861-1885)
;===========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  biome   = (/"broadleaf","needleleaf","evergreen","deciduous"/)
  casenm  = (/"Tbio", "cLeaf", "cWood", "Bbio"/)
  modelnm = (/"BNU-ESM", "HadGEM2-CC", "HadGEM2-ES", "IPSL-CM5A-LR", "IPSL-CM5A-MR", \
              "IPSL-CM5B-LR", "MIROC-ESM", "MIROC-ESM-CHEM"/)
  frac_num= (/0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)

  print ("Model, Biome, Component, Fraction, MassMod(Pg), MeanMOD(Pg), StdMOD(Pg), MassOBS(Pg), MeanOBS(Pg), StdOBS(Pg)")

;----------------------------------------------------------------------
;;; Start looping

do b=0, dimsizes(biome)-1
do z=0, dimsizes(frac_num)-1
do y=0, dimsizes(modelnm)-1
do x=0, dimsizes(casenm)-1

  indir20 = "../../../postanalysis/input_data/PI/individual/frac"+frac_num(z)+"/"

  fname2  = "mask_climate_all_"+modelnm(y)

  infile20= addfile (indir20+fname2+".nc","r")   ; biomass, forest frac

;;; Get variables
  lat = infile20->lat
  lon = infile20->lon

  obs_str = "mask_obs_biome_dens_"+biome(b)+"_"+casenm(x)
  obs_str2= "mask_obs_total_biome_"+biome(b)+"_"+casenm(x)
  model_str = "mask_model_biome_dens_"+biome(b)+"_"+casenm(x)
  model_str2= "mask_model_total_biome_"+biome(b)+"_"+casenm(x)

;;;;; forest frac > x%
  obs_biome_dens1 = tofloat(infile20->$obs_str$)
  obs_total_biome1 = tofloat(infile20->$obs_str2$)
  model_biome_dens1 = tofloat(infile20->$model_str$)
  model_total_biome1 = tofloat(infile20->$model_str2$)

;;; Mask for scatter data
  do p=0,dimsizes(lat)-1
    do q=0,dimsizes(lon)-1
      if ( ismissing(obs_total_biome1(p,q)) .or. ismissing(obs_biome_dens1(p,q)) ) then
        obs_total_biome1(p,q)= obs_total_biome1@_FillValue
        obs_biome_dens1(p,q) = obs_biome_dens1@_FillValue
      else
        obs_total_biome1(p,q)= obs_total_biome1(p,q)
        obs_biome_dens1(p,q) = obs_biome_dens1(p,q)
      end if

      if ( ismissing(model_total_biome1(p,q)) .or. ismissing(model_biome_dens1(p,q)) ) then
        model_total_biome1(p,q)= model_total_biome1@_FillValue
        model_biome_dens1(p,q) = model_biome_dens1@_FillValue
      else
        model_total_biome1(p,q)= model_total_biome1(p,q)
        model_biome_dens1(p,q) = model_biome_dens1(p,q)
      end if
    end do
  end do

;------------------------------------------------------------------------------
;;; Print out results
;  do p=0,dimsizes(lat)-1
;    do q=0,dimsizes(lon)-1
;      if (.not.ismissing(model_total_biome1(p,q)) .and. \
;          .not.ismissing(model_biome_dens1(p,q)) ) then
;        print (model_total_biome1(p,q)+", "+model_biome_dens1(p,q))
;      end if
;    end do
;  end do


   print ( modelnm(y)+", "+biome(b)+","+casenm(x)+", "+frac_num(z)+", " \
          +sum(model_total_biome1)/1.e12+", "+avg(model_total_biome1)/1.e12+", " \
          +stddev(model_total_biome1)/1.e12+", " \
          +sum(obs_total_biome1)/1.e12+", "+avg(obs_total_biome1)/1.e12+", " \
          +stddev(obs_total_biome1)/1.e12 )

;;; clean out pointers
  delete(obs_total_biome1)
  delete(obs_biome_dens1)
  delete(model_total_biome1)
  delete(model_biome_dens1)
  delete(lon)
  delete(lat)

end do
end do
end do
end do

end
exit

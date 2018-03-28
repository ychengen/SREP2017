;===========================================================
; Plot the upscaled forest fraction between BGI and CMIP5
;  Note1: only boreal and temperate forests are used.
;  Note2: BGI data is already forest type only, no fraction
;         needs to be involved; models' forest fractions are
;         required.
;  Note3: for modern time (1982-2005)
;===========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  biome   = (/"lump"/)
  casenm1  = (/"cVeg", "cLeaf", "cWood", "cRoot"/)
  casenm2  = (/"Tbio", "cLeaf", "cWood", "Bbio"/)
  modelnm = (/"BNU-ESM", "HadGEM2-CC", "HadGEM2-ES", "IPSL-CM5A-LR", "IPSL-CM5A-MR", \
              "IPSL-CM5B-LR", "MIROC-ESM", "MIROC-ESM-CHEM"/)
  frac_num= (/0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)

  print ("Model, Biome, Component, Fraction, MassMod(Pg), MeanMOD(Pg), StdMOD(Pg), MassOBS(Pg), MeanOBS(Pg), StdOBS(Pg)")

;----------------------------------------------------------------------
;;; Start looping

do b=0, dimsizes(biome)-1
do z=0, dimsizes(frac_num)-1
do y=0, dimsizes(modelnm)-1
do x=0, dimsizes(casenm2)-1

  indir10 = "../../CMIP5/climate_mean/modern/masked/frac"+frac_num(z)+"/"
  fname1  = "mask_climate_"+modelnm(y)+"_"+casenm1(x)+"_max"
  infile10= addfile (indir10+fname1+".nc","r")   ; biomass, forest frac
  
  indir20 = "../../../postanalysis/input_data/modern/frac"+frac_num(z)+"/"
  fname2  = "mask_climate_all_"+modelnm(y)
  infile20= addfile (indir20+fname2+".nc","r")   ; biomass, forest frac

;;; Get variables
  lat = infile20->lat
  lon = infile20->lon

  obs_str1= "mask_obs_total_biome_"+casenm1(x)
  obs_str2= "mask_obs_total_biome_"+casenm2(x)
  model_str2= "mask_model_total_biome_"+casenm2(x)

;;;;; forest frac > x%
  obs_total_biome1 = tofloat(infile10->$obs_str1$)
  obs_total_biome2 = tofloat(infile20->$obs_str2$)
  obs_total_biome1 = mask (obs_total_biome1,ismissing(obs_total_biome2),False)
  model_total_biome1 = tofloat(infile20->$model_str2$)

;;; Mask for scatter data
  do p=0,dimsizes(lat)-1
    do q=0,dimsizes(lon)-1
      if ( ismissing(obs_total_biome1(p,q)) ) then
        obs_total_biome1(p,q)= obs_total_biome1@_FillValue
      else
        obs_total_biome1(p,q)= obs_total_biome1(p,q)
      end if

      if ( ismissing(model_total_biome1(p,q)) ) then
        model_total_biome1(p,q)= model_total_biome1@_FillValue
      else
        model_total_biome1(p,q)= model_total_biome1(p,q)
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


   print ( modelnm(y)+", "+biome(b)+","+casenm2(x)+", "+frac_num(z)+", " \
          +sum(model_total_biome1)/1.e12+", "+avg(model_total_biome1)/1.e12+", " \
          +stddev(model_total_biome1)/1.e12+", " \
          +sum(obs_total_biome1)/1.e12+", "+avg(obs_total_biome1)/1.e12+", " \
          +stddev(obs_total_biome1)/1.e12 )

;;; clean out pointers
  delete(obs_total_biome1)
  delete(obs_total_biome2)
  delete(model_total_biome1)
  delete(lon)
  delete(lat)

end do
end do
end do
end do

end
exit

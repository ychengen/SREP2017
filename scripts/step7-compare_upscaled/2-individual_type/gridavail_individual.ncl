;==========================================================
; Plot the upscaled forest fraction between BGI and CMIP5
;  Note1: only boreal and temperate forests are used.
;  Note2: BGI data is already forest type only, no fraction
;         needs to be involved; models' forest fractions are
;         required.
;==========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  areamethod = 1   ;;; cell area: 0 from model output, 1 from upscaled observation
  comp    = (/"broadleaf","needleleaf","evergreen","deciduous"/)
  casenm  = (/"Tbio","Bbio","cLeaf","cWood"/)
  orignm  = (/"cVeg","cRoot","cLeaf","cWood"/)
  modelnm = (/"BNU-ESM","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR",\
              "IPSL-CM5B-LR","MIROC-ESM","MIROC-ESM-CHEM"/)

;;; Forest fraction threshold (deternmine if a forest grid)
  forest_frac = (/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)
;----------------------------------------------------------------------
;;; Start looping

print ("Model "+"PFT "+"Comp "+"Threshold "+"ForestArea(10^9ha) "+"TotalBio_OBS(Pg) "+"TotalBio_MOD(Pg) "+"GridAvail")

do h=0,dimsizes(comp)-1
do k=0,dimsizes(forest_frac)-1
do i=0,dimsizes(casenm)-1
do j=0,dimsizes(modelnm)-1

;;; land area info (from models)
  indir1 = "../../CMIP5/cellarea/"
  fname11 = "areacella_fx_"+modelnm(j)+"_historical_r0i0p0.nc"
  fname12 = "sftlf_fx_"+modelnm(j)+"_historical_r0i0p0.nc"
  infile11= addfile (indir1+fname11+".nc","r")
  infile12= addfile (indir1+fname12+".nc","r")

;;; land area info (from upscaled observations)
  if (areamethod.ne.0) then
    indir0 = "../../CMIP5/landcover_GLC2000_for_CMIP5/"
    fname0 = "landcover_BGI_"+casenm(i)+"_"+modelnm(j)+".nc"
    infile0= addfile (indir0+fname0+".nc","r")
  end if

;;; biomass info
  indir2 = "../../CMIP5/climate_mean/modern/masked/individual/frac"+forest_frac(k)+"/"
  fname2  = "mask_climate_"+modelnm(j)+"_"+casenm(i)
  infile2 = addfile (indir2+fname2+".nc","r")

;;; Get variables
  lat = infile11->lat
  lon = infile11->lon

  areacella  = infile11->areacella
  sftlf      = 0.01*(infile12->sftlf)
;  cell_area = infile1->cell_area

  str1 = "mask_obs_forestFrac_"+comp(h)+"_"+casenm(i)
  str2 = "mask_obs_total_biome_"+comp(h)+"_"+casenm(i)
  str3 = "mask_obs_biome_dens_"+comp(h)+"_"+casenm(i)
  str4 = "mask_model_forestFrac_"+comp(h)+"_"+casenm(i)
  str5 = "mask_model_total_biome_"+comp(h)+"_"+casenm(i)
  str6 = "mask_model_biome_dens_"+comp(h)+"_"+casenm(i)
  str7 = "mask_diff_forestFrac_"+comp(h)+"_"+casenm(i)
  str8 = "mask_diff_total_biome_"+comp(h)+"_"+casenm(i)
  str9 = "mask_diff_biome_dens_"+comp(h)+"_"+casenm(i)

  obs_forestFrac = infile2->$str1$
  temp_obs_total_biome= infile2->$str2$
  obs_biome_dens = infile2->$str3$
  model_forestFrac = infile2->$str4$
  temp_model_total_biome= infile2->$str5$
  model_biome_dens = infile2->$str6$
  diff_forestFrac = infile2->$str7$
  temp_diff_total_biome = infile2->$str8$
  diff_biome_dens = infile2->$str9$

  obs_total_biome = temp_obs_total_biome*1.e-12
  model_total_biome = temp_model_total_biome*1.e-12
  diff_total_biome = temp_diff_total_biome*1.e-12
  copyatt (obs_total_biome, temp_obs_total_biome) 
  copyatt (model_total_biome, temp_model_total_biome) 
  copyatt (diff_total_biome, temp_diff_total_biome) 
  obs_total_biome@unit = "Pg"
  model_total_biome@unit = "Pg"
  diff_total_biome@unit = "Pg"

;;; Choice of using cell area
  if (areamethod .eq. 0) then     ;; model
    cell_area = areacella*sftlf*model_forestFrac
  else                            ;; upscaled observation
    cell_area = infile0->cell_area
    cell_area = cell_area*doubletofloat(obs_forestFrac)
  end if

;  obs_forestFrac = lonFlip(obs_forestFrac)
;  obs_total_biome= lonFlip(obs_total_biome)
;  obs_biome_dens = lonFlip(obs_biome_dens)
;  model_forestFrac = lonFlip(model_forestFrac)
;  model_total_biome= lonFlip(model_total_biome)
;  model_biome_dens = lonFlip(model_biome_dens)
;  diff_forestFrac = lonFlip(diff_forestFrac)
;  diff_total_biome= lonFlip(diff_total_biome)
;  diff_biome_dens = lonFlip(diff_biome_dens)
;  diff_forestFrac = diff_forestFrac
;  diff_total_biome= diff_total_biome
;  diff_biome_dens = diff_biome_dens
;  gridavail = num(diff_forestFrac)
;  print (modelnm(j)+" "+casenm(i)+" "+forest_frac(k)+" "+gridavail)
  totbiomass_obs = sum(mask (obs_total_biome, ismissing(diff_forestFrac), False) )
  totbiomass_mod = sum(mask (model_total_biome, ismissing(diff_forestFrac), False) )
  totalarea = sum(mask (cell_area, ismissing(diff_forestFrac), False) )
  totalarea = totalarea/1.e4/1.e9
  gridavail1= num(diff_forestFrac)
  gridavail2= num(mask (cell_area, ismissing(diff_forestFrac), False))
  print (modelnm(j)+" "+comp(h)+" "+casenm(i)+" "+forest_frac(k)+" "+totalarea+" "+totbiomass_obs+" "+totbiomass_mod+" "+gridavail2)

;;; clean out pointers

  delete(obs_forestFrac)
  delete(temp_obs_total_biome)
  delete(obs_total_biome)
  delete(obs_biome_dens)
  delete(model_forestFrac)
  delete(temp_model_total_biome)
  delete(model_total_biome)
  delete(model_biome_dens)
  delete(diff_forestFrac)
  delete(temp_diff_total_biome)
  delete(diff_total_biome)
  delete(diff_biome_dens)
  delete(sftlf)
  delete(areacella)
  delete(cell_area)
  delete(lon)
  delete(lat)
  delete(infile2)
  delete(fname2)
  delete(fname12)
  delete(fname11)

end do
end do
end do
end do

end
exit

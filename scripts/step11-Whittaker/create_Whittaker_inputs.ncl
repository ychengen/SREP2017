;==========================================================
;  Create a new dataset as inputs of Whittaker diagram
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

;;; For forest grid definition, by determining the forest fraction (0.0 to 0.9)
;  forest_threshold = (/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)
  forest_threshold = (/0.0/)
  class = (/1,2,3,4,5,6,7,8,9/) 
  FillValue  = 9.96921e36

;----------------------------------------------------------------------
;;; Start looping

do k=0,dimsizes(forest_threshold)-1
do i=1,1 ; for Tbio only;  0,dimsizes(casenm)-1
do j=0, 3; dimsizes(modelnm)-1
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
  outdir = "../../CMIP5/climate_mean/modern/masked/individual/frac"+forest_threshold(k)+"/Whittaker/"

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
  temp_obs_pr  = obs_pr
  temp_obs_tas = obs_tas
  temp_obs_rsds =obs_rsds

;;; model
  temp2_forestFrac  = infile2->landCoverFrac
  temp2_biome_dens  = infile2->$orignm(i)$

  model_pr   = infile2->pr
  model_rsds = infile2->rsds
  model_tas  = infile2->tas
  temp_model_pr  = model_pr
  temp_model_tas = model_tas
  temp_model_rsds= model_rsds

;;; Mask model output to the range 30N-80N
  do zzz = 0,dimsizes(lat)-1
    if (lat(zzz).lt.30.0 .or. lat(zzz).gt.80.0) then
      temp2_forestFrac(:,zzz,:) = FillValue
      temp2_biome_dens(zzz,:) = FillValue
      obs_pr(zzz,:) = FillValue
      obs_tas(zzz,:) = FillValue
      obs_rsds(zzz,:) = FillValue
      model_pr(zzz,:) = FillValue
      model_tas(zzz,:) = FillValue
      model_rsds(zzz,:) = FillValue
    end if
  end do

  copyatt(obs_pr, temp_obs_pr)
  copyatt(obs_tas, temp_obs_tas)
  copyatt(obs_rsds, temp_obs_rsds)
  copyatt(model_pr, temp_model_pr)
  copyatt(model_tas, temp_model_tas)
  copyatt(model_rsds, temp_model_rsds)

  delete([/temp_model_pr, temp_model_tas, temp_model_rsds/])
  delete([/temp_obs_pr, temp_obs_tas, temp_obs_rsds/])
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
  temp1_forestFrac@_FillValue  = FillValue
  temp1_total_biome&lon =lon
  temp1_total_biome@_FillValue = FillValue
  temp2_forestFrac&lon =lon
  temp2_forestFrac@_FillValue  = FillValue
  temp2_biome_dens&lon =lon
  temp2_biome_dens@_FillValue  = FillValue

;;; Forest fraction
;;;;; Observation
  do c=0,21
    do p=0,dimsizes(lat)-1
      do q=0,dimsizes(lon)-1
        if (ismissing(temp1_forestFrac(c,p,q))) then
          temp1_forestFrac(c,p,q) = 0.
        else if (c .ge. 5) then
          temp1_forestFrac(c,p,q) = 0.
        end if
        end if
      end do
    end do
  end do

;;; new 9 categories for forest in GLC
  temp_obs_forestFrac = temp1_forestFrac(0:dimsizes(class)-1,:,:)
  obs_forestFrac      = mask( temp_obs_forestFrac,  \
                             temp_obs_forestFrac .le. forest_threshold(k), False)

  obs_forestFrac!0   = "class"
  obs_forestFrac!1   = "lat"
  obs_forestFrac!2   = "lon"
  obs_forestFrac&class = class
  obs_forestFrac&lat = lat
  obs_forestFrac&lon = lon

;;;;; Model
  temp_model_forestFrac = temp2_forestFrac(0:dimsizes(class)-1,:,:)

  model_forestFrac      = mask( temp_model_forestFrac,  \
                                temp_model_forestFrac .le. forest_threshold(k)*100., False )
  model_forestFrac      = 0.01*model_forestFrac
  model_forestFrac!0 = "class"
  model_forestFrac!1 = "lat"
  model_forestFrac!2 = "lon"
  model_forestFrac&class = class
  model_forestFrac&lat = lat
  model_forestFrac&lon = lon

;;; Forest total biomass
;;;;; Observation
  do c=0,21
    do p=0,dimsizes(lat)-1
      do q=0,dimsizes(lon)-1
        if (ismissing(temp1_total_biome(c,p,q))) then
          temp1_total_biome(c,p,q) = 0.
        else if (c .ge. 5) then
          temp1_total_biome(c,p,q) = 0.
        end if
        end if
      end do
    end do
  end do

  temp_obs_total_biome  = temp1_total_biome(0:dimsizes(class)-1,:,:)
  obs_total_biome       = mask( temp_obs_total_biome, obs_forestFrac.le.0., False)
  obs_total_biome!0     = "class"
  obs_total_biome!1     = "lat"
  obs_total_biome!2     = "lon"
  obs_total_biome&class = class
  obs_total_biome&lat   = lat
  obs_total_biome&lon   = lon

;;;;; Model
  model_total_biome     = model_forestFrac
  model_total_biome     = model_forestFrac@_FillValue

  do zzz=0, dimsizes(model_forestFrac(:,0,0))-1
    model_total_biome(zzz,:,:) = temp2_biome_dens(:,:)*areacella(:,:)*model_forestFrac(zzz,:,:)
  end do

  model_total_biome!0   = "class"
  model_total_biome!1   = "lat"
  model_total_biome!2   = "lon"
  model_total_biome&class = class
  model_total_biome&lat = lat
  model_total_biome&lon = lon

;-----------------------------------------------------------------------
;;; Assign attributes
  obs_forestFrac!0 = "class"
  obs_forestFrac!1 = "lat"
  obs_forestFrac!2 = "lon"
  obs_forestFrac&class = class
  obs_forestFrac&lat = lat
  obs_forestFrac&lon = lon
  obs_forestFrac@comment = "observation"
  obs_forestFrac@unit = ""
  obs_forestFrac@long_name = "forest fraction"

  model_forestFrac!0 = "class"
  model_forestFrac!1 = "lat"
  model_forestFrac!2 = "lon"
  model_forestFrac&class = class
  model_forestFrac&lat = lat
  model_forestFrac&lon = lon
  model_forestFrac@comment = "model"
  model_forestFrac@unit = ""
  model_forestFrac@long_name = "forest fraction"

  obs_total_biome!0 = "class"
  obs_total_biome!1 = "lat"
  obs_total_biome!2 = "lon"
  obs_total_biome&class = class
  obs_total_biome!1 = "lat"
  obs_total_biome&lat = lat
  obs_total_biome&lon = lon
  obs_total_biome@comment = "observation"
  obs_total_biome@unit = "kg"
  obs_total_biome@long_name = "forest biomass"

  model_total_biome!0 = "class"
  model_total_biome!1 = "lat"
  model_total_biome!2 = "lon"
  model_total_biome&class = class
  model_total_biome&lat = lat
  model_total_biome&lon = lon
  model_total_biome@comment = "model"
  model_total_biome@unit = "kg"
  model_total_biome@long_name = "forest biomass"

;-----------------------------------------------------------------------
;;; Write out data

  newvar11 = "obs_forestFrac_"+casenm(i)
  newvar12 = "obs_total_biome_"+casenm(i)
  newvar13 = "obs_pr_"+casenm(i)
  newvar14 = "obs_tas_"+casenm(i)
  newvar15 = "obs_rsds_"+casenm(i)
  newvar21 = "model_forestFrac_"+casenm(i)
  newvar22 = "model_total_biome_"+casenm(i)
  newvar23 = "model_pr_"+casenm(i)
  newvar24 = "model_tas_"+casenm(i)
  newvar25 = "model_rsds_"+casenm(i)

  outfile->lat = lat
  outfile->lon = lon
  outfile->$newvar11$ = obs_forestFrac
  outfile->$newvar12$ = obs_total_biome   ;(kg)
  outfile->$newvar13$ = obs_pr   ;(kg)
  outfile->$newvar14$ = obs_tas  ;(kg)
  outfile->$newvar15$ = obs_rsds ;(kg)
  outfile->$newvar21$ = model_forestFrac
  outfile->$newvar22$ = model_total_biome   ;(kg)
  outfile->$newvar23$ = model_pr   ;(kg)
  outfile->$newvar24$ = model_tas  ;(kg)
  outfile->$newvar25$ = model_rsds ;(kg)

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
  delete(areacella)
;  delete(sftlf)
  delete(temp2_forestFrac)
  delete(temp2_biome_dens)
  delete(obs_forestFrac)
  delete(temp_obs_forestFrac)
  delete(model_forestFrac)
  delete(temp_model_forestFrac)
  delete(model_total_biome)
  delete(newvar11)
  delete(newvar12)
  delete(newvar13)
  delete(newvar14)
  delete(newvar15)
  delete(newvar21)
  delete(newvar22)
  delete(newvar23)
  delete(newvar24)
  delete(newvar25)
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

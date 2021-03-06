;==========================================================
; Print the upscaled forest fraction between BGI and CMIP5
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

  casenm  = (/"Bbio","Tbio","cLeaf","cWood"/)
  orignm  = (/"cRoot","cVeg","cLeaf","cWood"/)
  modelnm = (/"BNU-ESM","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR",\
              "IPSL-CM5B-LR","MIROC-ESM-CHEM","MIROC-ESM"/)
;----------------------------------------------------------------------
;;; Start looping

do i=0,dimsizes(casenm)-1
do j=0,dimsizes(modelnm)-1
  print (casenm(i)+" "+modelnm(j))

  indir2 = "../../CMIP5/climate_mean/masked/"
  fname2  = "mask_climate_"+modelnm(j)+"_"+casenm(i)

;;; CMIP5 models all variables
  infile2 = addfile (indir2+fname2+".nc","r")

;;; Get variables
  lat = infile2->lat
  temp_lon = infile2->lon
  lon = lonFlip(temp_lon)

  obs_forestFrac = infile2->mask_obs_forestFrac 
  temp_obs_total_biome= infile2->mask_obs_total_biome
  obs_biome_dens = infile2->mask_obs_biome_dens 
  model_forestFrac = infile2->mask_model_forestFrac 
  temp_model_total_biome= infile2->mask_model_total_biome
  model_biome_dens = infile2->mask_model_biome_dens 
  diff_forestFrac = infile2->diff_forestFrac 
  temp_diff_total_biome = infile2->diff_total_biome
  diff_biome_dens = infile2->diff_biome_dens 

;  diff_forestFrac = model_forestFrac - obs_forestFrac
;  temp_diff_total_biome = temp_model_total_biome - temp_obs_total_biome
;  diff_biome_dens = model_biome_dens - obs_biome_dens
;  copyatt (diff_forestFrac, model_forestFrac) 
;  copyatt (temp_diff_total_biome, temp_model_total_biome) 
;  copyatt (diff_biome_dens, model_biome_dens) 

  obs_total_biome = temp_obs_total_biome*1.e-10
  model_total_biome = temp_model_total_biome*1.e-10
  diff_total_biome = temp_diff_total_biome*1.e-10
  copyatt (obs_total_biome, temp_obs_total_biome) 
  copyatt (model_total_biome, temp_model_total_biome) 
  copyatt (diff_total_biome, temp_diff_total_biome) 
  obs_total_biome@unit = "Pg"
  model_total_biome@unit = "Pg"
  diff_total_biome@unit = "Pg"

;  temp = tofloat(obs_forestFrac)
;  temp = tofloat(obs_total_biome)
;  temp = tofloat(obs_biome_dens)
  temp = tofloat(model_forestFrac)
;  temp = tofloat(model_total_biome)
;  temp = tofloat(model_biome_dens)
;  temp = tofloat(diff_forestFrac)
;  temp = tofloat(diff_total_biome)
;  temp = tofloat(diff_biome_dens)

  zz = 0
  do p=0,dimsizes(lat)-1
    do q=0,dimsizes(lon)-1
      if (.not.ismissing(temp(p,q))) then
        print(zz+" "+temp(p,q))
        zz = zz+1
      end if
    end do
  end do

;----------------------------------------------------------------------
;; Statistics
;
;  opt        = False
;  stat_obs_forestFrac   = stat_dispersion(obs_forestFrac, opt )
;  stat_model_forestFrac = stat_dispersion(model_forestFrac, opt )
;  stat_obs_total_biome  = stat_dispersion(obs_total_biome, opt )
;  stat_model_total_biome= stat_dispersion(model_total_biome, opt )
;  stat_obs_biome_dens   = stat_dispersion(obs_biome_dens, opt )
;  stat_model_biome_dens = stat_dispersion(model_biome_dens, opt )
;
;  nStat = num(.not.ismissing(stat_obs_forestFrac))
;  statLabel = new ( nStat, "string")
;
;  statLabel = (/ "    mean", "  stddev", "     min" \
;               , "  lowDec", "  LowOct", "  LowSex" \
;               , "LowQuart", "  LowTri" \
;               , "  Median" \
;               , "   HiTri", " HiQuart","   HiSex","   HiOct", "   HiDec" \
;               , "     Max", "   Range","  Disper","     RMS" \
;               , "    nTot", "    nUse","    nMsg","    %Msg" /)
;
;  print(statLabel+"   " \
;       +sprintf("%9.5f", stata(0:nStat-1))+"   " \
;       +sprintf("%9.5f", statb(0:nStat-1))+"   " \
;       +sprintf("%9.5f", statc(0:nStat-1)) )

;-----------------------------------------------------------------------
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
  delete(temp_lon)
  delete(lon)
  delete(lat)
  delete(infile2)
  delete(fname2)
  delete(temp)

end do
end do

end
exit

;***********************************************************************
; Purpose: compute the relative errors between modern time data
;      For modern model - modern observation!
;      Forest only
;
;***********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
                                  
begin

  biotype = (/"lumped"/)
  ffrac   = (/"0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"/)
  modelnm = (/"BNU-ESM","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR", \
              "IPSL-CM5B-LR","MIROC-ESM","MIROC-ESM-CHEM"/)
  comp    = (/"Tbio", "cLeaf", "cWood", "Bbio"/)
  missingval = 9.96921e36

  print ("Model, "+"Biome, "+"Comp, "+"Threshold, "+"meanRE_mass_grid, "+"meanRE_dens_grid, "+"meanRE_mass_avg,"+"meanRE_dens_avg")

  do btype=0,dimsizes(biotype)-1
   do i=0,dimsizes(ffrac)-1
    do k=0,dimsizes(comp)-1
     do j=0,dimsizes(modelnm)-1

    indir0 = "../../../postanalysis/input_data/modern/frac"+ffrac(i)+"/"
    indir1 = "../../../postanalysis/input_data/modern/frac"+ffrac(i)+"/"
    fname  = "mask_climate_all_"+modelnm(j)+".nc"
    infile0= addfile(indir0+fname,"r")
    infile1= addfile(indir1+fname,"r")

    string11 = "mask_obs_total_biome_"+comp(k)
    string21 = "mask_model_total_biome_"+comp(k)
    string12 = "mask_obs_biome_dens_"+comp(k)
    string22 = "mask_model_biome_dens_"+comp(k)
    obs_var1 = infile0->$string11$
    modern_var1  = infile1->$string21$
    obs_var2 = infile0->$string12$
    modern_var2  = infile1->$string22$

;;; Compute relative error, modern-MODERN OBS / MODERN OBS
    re_var11 = avg( (modern_var1 - obs_var1)/obs_var1 )
    re_var21 = avg( (modern_var2 - obs_var2)/obs_var2 )
    re_var12 = ( avg(modern_var1)-avg(obs_var1) ) / avg(obs_var1)
    re_var22 = ( avg(modern_var2)-avg(obs_var2) ) / avg(obs_var2)

    print (modelnm(j)+", "+biotype+", "+comp(k)+", "+ffrac(i)+", "+re_var11+", "+re_var21+", "+re_var12+", "+re_var22)

    delete(obs_var1)
    delete(modern_var1)
    delete(obs_var2)
    delete(modern_var2)
    delete(re_var11)
    delete(re_var12)
    delete(re_var21)
    delete(re_var22)

     end do
    end do
   end do
  end do

end

;***********************************************************************
; Purpose: compute the relative errors between PI and modern time data
;      For preindustrial model - modern observation!
;      Individual only
;
;***********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
                                  
begin

  biotype = (/"broadleaf", "needleleaf", "evergreen", "deciduous"/)
  ffrac   = (/"0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"/)
  modelnm = (/"BNU-ESM","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR", \
              "IPSL-CM5B-LR","MIROC-ESM","MIROC-ESM-CHEM"/)
  comp    = (/"Tbio", "cLeaf", "cWood", "Bbio"/)
  missingval = 9.96921e36

  print ("Model, "+"Biome, "+"Comp, "+"Threshold, "+"meanRE_mass_grid, "+"meanRE_dens_grid, "+"meanRE_mass_avg, "+"meanRE_dens_avg")

  do btype=0,dimsizes(biotype)-1
   do i=0,dimsizes(ffrac)-1
    do k=0,dimsizes(comp)-1
     do j=0,dimsizes(modelnm)-1

    indir0 = "../../../postanalysis/input_data/modern/individual/frac"+ffrac(i)+"/"
    indir1 = "../../../postanalysis/input_data/PI/individual/frac"+ffrac(i)+"/"
    fname  = "mask_climate_all_"+modelnm(j)+".nc"
    infile0= addfile(indir0+fname,"r")
    infile1= addfile(indir1+fname,"r")

    string11 = "mask_obs_total_biome_"+biotype(btype)+"_"+comp(k)
    string21 = "mask_model_total_biome_"+biotype(btype)+"_"+comp(k)
    string12 = "mask_obs_biome_dens_"+biotype(btype)+"_"+comp(k)
    string22 = "mask_model_biome_dens_"+biotype(btype)+"_"+comp(k)
    obs_var1 = infile0->$string11$
    PI_var1  = infile1->$string21$
    obs_var2 = infile0->$string12$
    PI_var2  = infile1->$string22$

;;; Compute relative error, PI-MODERN / MODERN
    re_var11 = avg( (PI_var1 - obs_var1)/obs_var1 )
    re_var21 = avg( (PI_var2 - obs_var2)/obs_var2 )
    re_var12 = ( avg(PI_var1)-avg(obs_var1) ) / avg(obs_var1)
    re_var22 = ( avg(PI_var2)-avg(obs_var2) ) / avg(obs_var2)

    print (modelnm(j)+", "+biotype(btype)+", "+comp(k)+", "+ffrac(i)+", "+re_var11+", "+re_var21+", "+re_var12+", "+re_var22)

    delete(obs_var1)
    delete(PI_var1)
    delete(obs_var2)
    delete(PI_var2)
    delete(re_var11)
    delete(re_var12)
    delete(re_var21)
    delete(re_var22)

     end do
    end do
   end do
  end do

end

;***********************************************************************
; Purpose: Find the ratio and cross correlation coeff for Taylor diagram
;    ratio are Case_Variance/Reference_Variance
;    cc    are the cross correlation coef of Case to Reference
;
;***********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
                                  
begin

  print ("Model "+" Threshold "+" Var_OBS"+" Var_MOD"+" Ratio"+" Corr"+" AvailGrid")
  ffrac   = (/"0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"/)
  modelnm = (/"BNU-ESM","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR", \
              "IPSL-CM5B-LR","MIROC-ESM","MIROC-ESM-CHEM"/)
  comp    = (/"Tbio"/)
  missingval = 9.96921e36

  do i=0,dimsizes(ffrac)-1
    do j=0,dimsizes(modelnm)-1
      do k=0,dimsizes(comp)-1

    indir0 = "../../CMIP5/landcover_GLC2000_for_CMIP5/"
    fname0 = "landcover_BGI_"+comp(k)+"_"+modelnm(j)+".nc"
    infile0= addfile(indir0+fname0,"r")

    indir1 = "../../postanalysis/input_data/modern/frac"+ffrac(i)+"/"
    fname1 = "mask_climate_all_"+modelnm(j)+".nc"
    infile1= addfile(indir1+fname1,"r")

    indir99  = "../../CMIP5/cellarea/"
    fname991 = "areacella_fx_"+modelnm(j)+"_historical_r0i0p0.nc"
    fname992 = "sftlf_fx_"+modelnm(j)+"_historical_r0i0p0.nc"
    infile991= addfile (indir99+fname991,"r")
    infile992= addfile (indir99+fname992,"r")

    var_obs   = "mask_obs_biome_dens_"+comp(k)
    var_model   = "mask_model_biome_dens_"+comp(k)
;    var_obs   = "mask_obs_total_biome_"+comp(k)
;    var_model   = "mask_model_total_biome_"+comp(k)

    cellarea = infile0->cell_area     ;;; from unpscaled observation data
    areacella= infile991->areacella
    sftlf    = 0.01*(infile992->sftlf)
    var1 = infile1->$var_obs$
    var2 = infile1->$var_model$

;    print(num(var1)+" "+num(var2))
;;; Compute areal weights
    tempweight = mask (cellarea, ismissing(var1), False)
;    tempweight = mask (areacella*sftlf, ismissing(var1), False)
    if (num(tempweight).gt.0) then
      cellweight = tempweight/sum(tempweight)
;      cellweight = tempweigh/sum(tempweight)
    else
      cellweight = missingval
    end if
    availgrid = num(cellweight)

    var11 = var1*cellweight
    var22 = var2*cellweight

;;; Compute statistics
    if (all(ismissing(var11)) .or. all(ismissing(var22))) then
      ratio = missingval
      cc    = missingval
    else
      variance_obs   = variance(var11)
      variance_model = variance(var22)

      if (variance_obs.eq.0. .or. variance_model.eq.0.) then
        ratio = missingval
        cc    = missingval
      else
        ratio = variance_model / variance_obs
;        cc    = pattern_cor (var1, var2, 1.0, 1)   ;; no weight, uncentered
;        cc    = pattern_cor (var1, var2, 1.0, 0)   ;; no weight, centered
        cc    = pattern_cor (var1, var2, cellweight, 1)   ;; weight, uncentered
      end if
    end if

    print (modelnm(j)+" "+ffrac(i)+" "+variance_obs+" "+variance_model+" "+ratio+" "+cc+" "+availgrid)

    delete(var1)
    delete(var2)
    delete(var11)
    delete(var22)
    delete(areacella)
    delete(sftlf)
    delete(cellarea)
    delete(tempweight)
    delete(cellweight)

      end do
    end do
  end do

end

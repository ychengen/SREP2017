;;;; ONLY WORKS with NCL 6.4.0 or higher version!!!
;;;; Modified to 4 categories only: broadleaf evergreen, needleleaf evergreen,
;;;;                                broadleaf deciduous, needleleaf decidious
;;;; Plot observed biomass on modeled climate space
;----------------------------------------------------------------------
; Procedure for adding a labelbar at a given NDC location, given
; the levels and colors to use.
;----------------------------------------------------------------------
undef("add_labelbar")
procedure add_labelbar(wks,plot,levels,colors,transpar)
local lbres, labels

begin

  nlevels = dimsizes(levels)
;----------------------------------------------------------------------
; Draw a labelbar
;----------------------------------------------------------------------
  lbres                    = True
  lbres@vpWidthF           = 0.12             ; width
  lbres@vpHeightF          = 0.40             ; height
  lbres@lbPerimOn          = False            ; Turn off perimeter.
  lbres@lbOrientation      = "Vertical" ; "Horizontal"     ; Default is vertical.
  lbres@lbLabelAlignment   = "BoxCenters"; "InteriorEdges"  ; Default is "BoxCenters"
  lbres@lbFillColors       = colors           ; Colors for boxes.
  lbres@lbMonoFillPattern  = True             ; Fill them all solid.
  lbres@lbLabelFontHeightF = 0.018            ; label font height
  lbres@lbFillOpacityF     = transpar

;  labels = sprintf("%2g",levels)
  labels = levels

;  lbid = gsn_create_labelbar(wks,nlevels+1,labels,lbres)
  lbid = gsn_create_labelbar(wks,nlevels,labels,lbres)

; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot. Here, we are using the top right
; corner of the labelbar as the point which we are going to position
; it, and then we use amParallelPosF and amOrthogonalPosF to indicate
; where we want to place it.
;
; amParallelPosF/amOrthogonalPosF                                       ;
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;
  amres                  = True
  amres@amJust           = "TopCenter"
  amres@amParallelPosF   =  0.65   ; keep labelbar centered
  amres@amOrthogonalPosF =  -0.30  ; move down and outside of plot

;
; Give both annotation id and labelbar id unique names.
;
; Attaching them to plot with unique names ensures that
; labelbar "lives" outside this procedure.
;
  tmpid1 = "anno"+unique_string("id")
  tmpid2 = "lbar"+unique_string("id")

  plot@$tmpid1$ = gsn_add_annotation(plot,lbid,amres)
  plot@$tmpid2$ = lbid
end

; ====================================================================================
;;; Main codes

begin

  transpar = 0.7
  threshold = 0.2
  modelnm = (/"BNU-ESM","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR",\
              "IPSL-CM5B-LR", "MIROC-ESM","MIROC-ESM-CHEM"/)
  casenm  = "ObsBiom_ModClim"
  casenm1 = "obs"  ;; obs or model
  casenm2 = "model"  ;; obs or model
  comp   = "Tbio"

do mmm = 0, dimsizes(modelnm)-1

;  indir  = "/lustre/atlas2/cli017/scratch/ceyang1/climate_mean/modern/masked/individual/frac0/whittaker/"
  indir  = "./"
  infile = addfile(indir+"mask_climate_"+modelnm(mmm)+"_"+comp+".nc","r")

  var1= casenm2+"_pr_"+comp
  var2= casenm2+"_tas_"+comp
;  var3= casenm+"_rsds_"+comp
  var4= casenm1+"_forestFrac_"+comp
  var5= casenm1+"_total_biome_"+comp

  lat  = infile->lat
  lon  = infile->lon
  pr   = tofloat(infile->$var1$)     ;; (lat, lon)
  tas  = tofloat(infile->$var2$)     ;; (lat, lon)
;  rsds = tofloat(infile->$var3$)     ;; (lat, lon)
  frac = tofloat(infile->$var4$)     ;; (class, lat, lon)
  bmass= tofloat(infile->$var5$)     ;; (class, lat, lon)

  pr   = pr*86400.*365./10.   ;;; modeled climate, convert kg m-2 s-1 to cm/yr
  tas  = tas-273.       ;;; convert K to Celcius
  bmass= bmass/1.e12    ;;; convert kg C to Pg C

;-------Options--------
;;; new 4 categories: 1-BET, 2-NET, 3-BDT, 4-NDT
;;; new 2 categories for HadGEM: 1-broadleaf. 2-needleleaf
;;; OBS1: 1->1, 4->2, 2+3->3, 5->4
;;; OBS2: 1+2+3->1, 4+5->2
;;; MOHC: 1->1, 2->2
;;; BNU: 4+5->1, 1+2->2, 6+7+8->3, 3->4
;;; IPSL: 2+5->1, 4+7->2, 3+6+8->3, 9->4
;;; MIROC: 1->1, 2+4->2, 3->3, 5->4

  if (mmm.eq.1 .or. mmm.eq.2) then
    levels = (/"BT","NT"/)
  else
    levels = (/"BET","NET","BDT","NDT"/)
  end if
  nlevels = dimsizes(levels)

  temp_pr   = new ((/nlevels,dimsizes(lat),dimsizes(lon)/),float)
  temp_tas  = new ((/nlevels,dimsizes(lat),dimsizes(lon)/),float)
  temp_frac = new ((/nlevels,dimsizes(lat),dimsizes(lon)/),float)
  temp_bmass= new ((/nlevels,dimsizes(lat),dimsizes(lon)/),float)

  if (mmm.eq.1 .or. mmm.eq.2) then
    temp_frac(0,:,:) = dim_sum_n(frac(0:2,:,:),0)
    temp_frac(1,:,:) = dim_sum_n(frac(3:4,:,:),0)
    temp_bmass(0,:,:) = dim_sum_n(bmass(0:2,:,:),0)
    temp_bmass(1,:,:) = dim_sum_n(bmass(3:4,:,:),0)
  else
    temp_frac(0,:,:) = frac(0,:,:)
    temp_frac(1,:,:) = frac(3,:,:)
    temp_frac(2,:,:) = dim_sum_n(frac(1:2,:,:),0)
    temp_frac(3,:,:) = frac(4,:,:)
    temp_bmass(0,:,:) = bmass(0,:,:)
    temp_bmass(1,:,:) = bmass(3,:,:)
    temp_bmass(2,:,:) = dim_sum_n(bmass(1:2,:,:),0)
    temp_bmass(3,:,:) = bmass(4,:,:)
  end if

  do nnn = 0,nlevels-1
    temp_pr(nnn,:,:)   = mask( pr(:,:), (ismissing(pr(:,:)) .or. ismissing(tas(:,:)) .or. \
                               ismissing(temp_bmass(nnn,:,:)) .or. temp_bmass(nnn,:,:).le.0. \
                               .or. ismissing(temp_frac(nnn,:,:)) .or. \
                               temp_frac(nnn,:,:).le.threshold), False )
    temp_tas(nnn,:,:)  = mask( tas(:,:), (ismissing(pr(:,:)) .or. ismissing(tas(:,:)) .or. \
                               ismissing(temp_bmass(nnn,:,:)) .or. temp_bmass(nnn,:,:).le.0. \
                               .or. ismissing(temp_frac(nnn,:,:)) .or. \
                               temp_frac(nnn,:,:).le.threshold), False )
  end do

;====================================================
  if (mmm.eq.1 .or. mmm.eq.2) then
    p1 = ndtooned (temp_pr(0,:,:))
    p2 = ndtooned (temp_pr(1,:,:))
    t1 = ndtooned (temp_tas(0,:,:))
    t2 = ndtooned (temp_tas(1,:,:))
    b1 = ndtooned (temp_bmass(0,:,:))
    b2 = ndtooned (temp_bmass(1,:,:))
  else
    p1 = ndtooned (temp_pr(0,:,:))
    p2 = ndtooned (temp_pr(1,:,:))
    p3 = ndtooned (temp_pr(2,:,:))
    p4 = ndtooned (temp_pr(3,:,:))
    t1 = ndtooned (temp_tas(0,:,:))
    t2 = ndtooned (temp_tas(1,:,:))
    t3 = ndtooned (temp_tas(2,:,:))
    t4 = ndtooned (temp_tas(3,:,:))
    b1 = ndtooned (temp_bmass(0,:,:))
    b2 = ndtooned (temp_bmass(1,:,:))
    b3 = ndtooned (temp_bmass(2,:,:))
    b4 = ndtooned (temp_bmass(3,:,:))
  end if
  delete ([/temp_pr,temp_tas,temp_frac,temp_bmass,pr,tas,frac,bmass/])

;========================================================================
  wks = gsn_open_wks("ps","Whittaker_"+casenm+"_"+modelnm(mmm)+"_"+comp+"_"+threshold)
;  colors = (/"red","blue","orange","darkgreen","cyan", \
;             "magenta","gold","brown","mediumpurple2"/)
;  marksty= (/16, 7, 6, 8, 15, 3, 0, 9/)
  colors = (/"red","blue","orange","darkgreen"/)
  marksty= (/4, 7, 6, 9/)

  res                   = True
  res@gsnDraw           = False
  res@gsnFrame          = False
;  res@gsnMaximize       = True           ; maximize plot
  res@xyMarkLineMode     = "Markers"
  res@xyMarkerThicknessF = 2.5
;  res@xyMarker          = marksty
  res@xyMarkerOpacityF  = transpar
;  res@xyMarkerColors    = "darkgreen"
  res@trXReverse        = True
  res@trXMinF           = -30.
  res@trXMaxF           = 30.
  res@trYMinF           = 0.
  res@trYMaxF           = 500.

;  res@gsnStringFont    = "helvatica"
  res@tiYAxisString    = "Annual precipitation (cm)"
  res@tiXAxisString    = "Temperature (~S~o~N~C)"
  res@tiMainString     = modelnm(mmm); + " - " + comp
  res@gsnLeftStringFontHeightF = 0.018
  res@gsnRightStringFontHeightF = 0.018
  res@gsnLeftString    = casenm
  res@gsnRightString   = "F~B~f~N~ = "+threshold

  do gg = 1,nlevels
    if (gg.eq.1) then
      res@xyMarkerColors = colors(0)
      res@xyMarker       = marksty(0)
      plot1 = gsn_csm_xy(wks,t1,p1,res)
    else if (gg.eq.2) then
      res@xyMarkerColors = colors(1)
      res@xyMarker       = marksty(1)
      plot2 = gsn_csm_xy(wks,t2,p2,res)
    else if (gg.eq.3) then
      res@xyMarkerColors = colors(2)
      res@xyMarker       = marksty(2)
      plot3 = gsn_csm_xy(wks,t3,p3,res)
    else if (gg.eq.4) then
      res@xyMarkerColors = colors(3)
      res@xyMarker       = marksty(3)
      plot4 = gsn_csm_xy(wks,t4,p4,res)
    end if
    end if
    end if
    end if
  end do

  if (mmm.eq.1 .or. mmm.eq.2) then
    overlay(plot1,plot2)
  else
    overlay(plot1,plot2)
    overlay(plot1,plot3)
    overlay(plot1,plot4)
  end if

  add_labelbar(wks,plot1,levels,colors,transpar)
  draw(plot1)
  frame(wks)

  if (mmm.eq.1 .or. mmm.eq.2) then
    delete ([/p1,t1,p2,t2/])
    delete ([/b1,b2/])
  else
    delete ([/p1,t1,p2,t2,p3,t3,p4,t4/])
    delete ([/b1,b2,b3,b4/])
  end if
  delete ([/lat, lon, levels, nlevels/])
end do

end

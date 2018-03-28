;************************************************
; taylor.ncl: requires NCL 4.2.0 or higher
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "./2-taylor_diagram_select.ncl"    ;; only upto ratio=2.5 and Ff=0, 0.1, 0.5, 0.9
;************************************************

begin
; variables compared
  var       = (/ "BNU-ESM","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR",\
                  "IPSL-CM5B-LR", "MIROC-ESM", "MIROC-ESM-CHEM" /) 
  nVar      = dimsizes(var)               ; # of Variables

  cases     = (/"F~B~f~N~=0","F~B~f~N~=0.1","F~B~f~N~=0.5","F~B~f~N~=0.8","F~B~f~N~=0.9"/)
  nCase     = dimsizes( cases  )

; *******************************************************************
; Open wks and set plot options
; *******************************************************************
  wks   = gsn_open_wks("pdf","taylor_diagram_select") 
  Colors  = (/ "red","blue4","magenta","brown","black"/)  
  plot  = new(1,graphic)        

  opt            = True                                  ; taylor diagram with options
  opt@Markers    = (/16,16,16,16,16/)     ; make all "cases" solid fill
  opt@Colors     = Colors ; (/"red","blue","green","magenta","orange" /)          
  opt@varLabels  = var
  opt@caseLabels = cases ; (/"Case A", "Case B","Case C","Case D","Case E"/)
  
  opt@taylorDraw  = True                           ; don't draw
  opt@taylorFrame = True                           ; don't advance frame

; *******************************************************************
; Generate plot.
;   ratio are the ratio:  Case_Variance/Reference_Variance
;   cc    are the cross correlation coef of Case to Reference
; *******************************************************************
  ratio     = asciiread ("ratio_select.txt",(/nCase, nVar/),"float")
  cc        = asciiread ("corr_select.txt", (/nCase, nVar/),"float")

  mask_ratio= mask (ratio, ratio.gt.2.5 .or. ratio.le.0.0, False)
  mask_cc   = mask (cc   , ismissing(mask_ratio), False)
  plot      = taylor_diagram(wks,mask_ratio,mask_cc,opt)

exit
end 

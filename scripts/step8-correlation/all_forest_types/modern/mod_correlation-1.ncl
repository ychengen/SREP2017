;==========================================================
; Compute statistics of the BGI dataset
; For modern time (1982-2005)
;==========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  modelnm = (/"BNU-ESM", "HadGEM2-CC", "HadGEM2-ES", "IPSL-CM5A-LR", \
              "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC-ESM-CHEM", "MIROC-ESM"/)
;  forest_frac = (/0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)
  forest_frac = (/0, 0.2, 0.5, 0.8/)

;;; Start looping
do fc=0,dimsizes(forest_frac)-1
  indir = "../../../../postanalysis/"
  indir1= indir+"input_data/modern/frac"+forest_frac(fc)+"/"

  do z=0,dimsizes(modelnm)-1
    print (modelnm(z)+" ")

    fname = "mask_climate_all_"+modelnm(z)+".nc"
    infile= addfile (indir1+fname,"r")

    temp_lat  = infile->lat
    temp_lon  = infile->lon

    lat = tofloat(temp_lat)
    lon = tofloat(temp_lon)
    lat!0 = "lat"
    lat@units = "degree_north"
    lat@long_name = "latitude"
    lon!0 = "lon"
    lon@units = "degree_east"
    lon@long_name = "longitude"

;;; Meteorology
      model_pr = tofloat(infile->mask_model_pr)
      model_rsds = tofloat(infile->mask_model_rsds)
      model_tas = tofloat(infile->mask_model_tas)

;;; Biomass
      model_Bbio = tofloat(infile->mask_model_biome_dens_Bbio)
      model_Tbio = tofloat(infile->mask_model_biome_dens_Tbio)
      model_cLeaf= tofloat(infile->mask_model_biome_dens_cLeaf)
      model_cWood= tofloat(infile->mask_model_biome_dens_cWood)
      model_Abio = model_Tbio - model_Bbio

      temp_model_Bbio = mask(model_Bbio, model_Bbio.le.0. , False)
      temp_model_cWood = mask(model_cWood, model_cWood.le.0. , False)

      model_ratio0 = model_Abio/temp_model_Bbio
      model_ratio1 = model_cLeaf/temp_model_cWood
      model_ratio2 = model_cLeaf/temp_model_Bbio
      model_ratio3 = model_Bbio/temp_model_cWood

;;; Output file
      outdir  = indir+"correlation/all_forest_types/modern/frac"+forest_frac(fc)+"/"
      fout    = outdir+"correlation-"+modelnm(z)+"-MOD-1.nc"
      if (isfilepresent(fout))
        print("Warning: '"+fout+"' exists.")
        print("Will remove it.")
        system ("/bin/rm -f "+fout)
      end if
      outfile = addfile (fout,"c")
      outfile->lat = lat
      outfile->lon = lon

;----------------------------------------------------------
;;; Using 11-grid moving window, note at north-south bounday
;;; grids are less than 5
;;; Boundary at Southern hemisphere is not modified!!

;;;;; PR vs TAS

  corr_model_prtas = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prtas = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_tas(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_tas(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prtas(i,j) = 9.96921e+36
        else
          temp_corr_model_prtas(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prtas(i,j) = 9.96921e+36
        else
          temp_corr_model_prtas(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prtas(i,j) = temp_corr_model_prtas(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prtas(i,j) = temp_corr_model_prtas(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prtas(i,j) = 9.96921e+36
        else
          corr_model_prtas(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_tas(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prtas(i,j) = 9.96921e+36
        else
          corr_model_prtas(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_tas(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prtas = mask(corr_model_prtas, ismissing(model_pr), False)
  corr_model_prtas!0   = "lat"
  corr_model_prtas!1   = "lon"
  corr_model_prtas&lat = lat
  corr_model_prtas&lon = lon
  corr_model_prtas@comment = "Model"
  corr_model_prtas@long_name = "Correlation of precipitation and surface temperature"
  outfile->corr_model_prtas = corr_model_prtas
  delete(corr_model_prtas)
  delete(temp_corr_model_prtas)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs RSDS

  corr_model_prrsds = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prrsds = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_rsds(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_rsds(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prrsds(i,j) = 9.96921e+36
        else
          temp_corr_model_prrsds(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prrsds(i,j) = 9.96921e+36
        else
          temp_corr_model_prrsds(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prrsds(i,j) = temp_corr_model_prrsds(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prrsds(i,j) = temp_corr_model_prrsds(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prrsds(i,j) = 9.96921e+36
        else
          corr_model_prrsds(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_rsds(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prrsds(i,j) = 9.96921e+36
        else
          corr_model_prrsds(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_rsds(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prrsds = mask(corr_model_prrsds, ismissing(model_pr), False)
  corr_model_prrsds!0   = "lat"
  corr_model_prrsds!1   = "lon"
  corr_model_prrsds&lat = lat
  corr_model_prrsds&lon = lon
  corr_model_prrsds@comment = "Model"
  corr_model_prrsds@long_name = "Correlation of precipitation and shortwave radiation"
  outfile->corr_model_prrsds = corr_model_prrsds
  delete(corr_model_prrsds)
  delete(temp_corr_model_prrsds)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs Abio

  corr_model_prAbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prAbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_Abio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_Abio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prAbio(i,j) = 9.96921e+36
        else
          temp_corr_model_prAbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prAbio(i,j) = 9.96921e+36
        else
          temp_corr_model_prAbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prAbio(i,j) = temp_corr_model_prAbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prAbio(i,j) = temp_corr_model_prAbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prAbio(i,j) = 9.96921e+36
        else
          corr_model_prAbio(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_Abio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prAbio(i,j) = 9.96921e+36
        else
          corr_model_prAbio(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_Abio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prAbio = mask(corr_model_prAbio, ismissing(model_pr), False)
  corr_model_prAbio!0   = "lat"
  corr_model_prAbio!1   = "lon"
  corr_model_prAbio&lat = lat
  corr_model_prAbio&lon = lon
  corr_model_prAbio@comment = "Model"
  corr_model_prAbio@long_name = "Correlation of precipitation and aboveground biomass"
  outfile->corr_model_prAbio = corr_model_prAbio
  delete(corr_model_prAbio)
  delete(temp_corr_model_prAbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs Bbio

  corr_model_prBbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prBbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_Bbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_Bbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prBbio(i,j) = 9.96921e+36
        else
          temp_corr_model_prBbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prBbio(i,j) = 9.96921e+36
        else
          temp_corr_model_prBbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prBbio(i,j) = temp_corr_model_prBbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prBbio(i,j) = temp_corr_model_prBbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prBbio(i,j) = 9.96921e+36
        else
          corr_model_prBbio(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_Bbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prBbio(i,j) = 9.96921e+36
        else
          corr_model_prBbio(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_Bbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prBbio = mask(corr_model_prBbio, ismissing(model_pr), False)
  corr_model_prBbio!0   = "lat"
  corr_model_prBbio!1   = "lon"
  corr_model_prBbio&lat = lat
  corr_model_prBbio&lon = lon
  corr_model_prBbio@comment = "Model"
  corr_model_prBbio@long_name = "Correlation of precipitation and belowground biomass"
  outfile->corr_model_prBbio = corr_model_prBbio
  delete(corr_model_prBbio)
  delete(temp_corr_model_prBbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs Tbio

  corr_model_prTbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prTbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_Tbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_Tbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prTbio(i,j) = 9.96921e+36
        else
          temp_corr_model_prTbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prTbio(i,j) = 9.96921e+36
        else
          temp_corr_model_prTbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prTbio(i,j) = temp_corr_model_prTbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prTbio(i,j) = temp_corr_model_prTbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prTbio(i,j) = 9.96921e+36
        else
          corr_model_prTbio(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_Tbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prTbio(i,j) = 9.96921e+36
        else
          corr_model_prTbio(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_Tbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prTbio = mask(corr_model_prTbio, ismissing(model_pr), False)
  corr_model_prTbio!0   = "lat"
  corr_model_prTbio!1   = "lon"
  corr_model_prTbio&lat = lat
  corr_model_prTbio&lon = lon
  corr_model_prTbio@comment = "Model"
  corr_model_prTbio@long_name = "Correlation of precipitation and total biomass"
  outfile->corr_model_prTbio = corr_model_prTbio
  delete(corr_model_prTbio)
  delete(temp_corr_model_prTbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs cLeaf

  corr_model_prcLeaf = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prcLeaf = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_cLeaf(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_cLeaf(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prcLeaf(i,j) = 9.96921e+36
        else
          temp_corr_model_prcLeaf(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prcLeaf(i,j) = 9.96921e+36
        else
          temp_corr_model_prcLeaf(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prcLeaf(i,j) = temp_corr_model_prcLeaf(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prcLeaf(i,j) = temp_corr_model_prcLeaf(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prcLeaf(i,j) = 9.96921e+36
        else
          corr_model_prcLeaf(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prcLeaf(i,j) = 9.96921e+36
        else
          corr_model_prcLeaf(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_cLeaf(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prcLeaf = mask(corr_model_prcLeaf, ismissing(model_pr), False)
  corr_model_prcLeaf!0   = "lat"
  corr_model_prcLeaf!1   = "lon"
  corr_model_prcLeaf&lat = lat
  corr_model_prcLeaf&lon = lon
  corr_model_prcLeaf@comment = "Model"
  corr_model_prcLeaf@long_name = "Correlation of precipitation and leaf biomass"
  outfile->corr_model_prcLeaf = corr_model_prcLeaf
  delete(corr_model_prcLeaf)
  delete(temp_corr_model_prcLeaf)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs cWood

  corr_model_prcWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prcWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prcWood(i,j) = 9.96921e+36
        else
          temp_corr_model_prcWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prcWood(i,j) = 9.96921e+36
        else
          temp_corr_model_prcWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prcWood(i,j) = temp_corr_model_prcWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prcWood(i,j) = temp_corr_model_prcWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prcWood(i,j) = 9.96921e+36
        else
          corr_model_prcWood(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prcWood(i,j) = 9.96921e+36
        else
          corr_model_prcWood(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prcWood = mask(corr_model_prcWood, ismissing(model_pr), False)
  corr_model_prcWood!0   = "lat"
  corr_model_prcWood!1   = "lon"
  corr_model_prcWood&lat = lat
  corr_model_prcWood&lon = lon
  corr_model_prcWood@comment = "Model"
  corr_model_prcWood@long_name = "Correlation of precipitation and wood biomass"
  outfile->corr_model_prcWood = corr_model_prcWood
  delete(corr_model_prcWood)
  delete(temp_corr_model_prcWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; PR vs ratio0 (Abio/Bbio)

  corr_model_prratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prratio0(i,j) = 9.96921e+36
        else
          temp_corr_model_prratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prratio0(i,j) = 9.96921e+36
        else
          temp_corr_model_prratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prratio0(i,j) = temp_corr_model_prratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prratio0(i,j) = temp_corr_model_prratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prratio0(i,j) = 9.96921e+36
        else
          corr_model_prratio0(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prratio0(i,j) = 9.96921e+36
        else
          corr_model_prratio0(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prratio0 = mask(corr_model_prratio0, ismissing(model_pr), False)
  corr_model_prratio0!0   = "lat"
  corr_model_prratio0!1   = "lon"
  corr_model_prratio0&lat = lat
  corr_model_prratio0&lon = lon
  corr_model_prratio0@comment = "Model, ratio=Abio/Bbio"
  corr_model_prratio0@long_name = "Correlation of precipitation and Abio/Bbio"
  outfile->corr_model_prratio0 = corr_model_prratio0
  delete(corr_model_prratio0)
  delete(temp_corr_model_prratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; PR vs ratio1 (cLeaf/Bbio)

  corr_model_prratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prratio1(i,j) = 9.96921e+36
        else
          temp_corr_model_prratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prratio1(i,j) = 9.96921e+36
        else
          temp_corr_model_prratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prratio1(i,j) = temp_corr_model_prratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prratio1(i,j) = temp_corr_model_prratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prratio1(i,j) = 9.96921e+36
        else
          corr_model_prratio1(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prratio1(i,j) = 9.96921e+36
        else
          corr_model_prratio1(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prratio1 = mask(corr_model_prratio1, ismissing(model_pr), False)
  corr_model_prratio1!0   = "lat"
  corr_model_prratio1!1   = "lon"
  corr_model_prratio1&lat = lat
  corr_model_prratio1&lon = lon
  corr_model_prratio1@comment = "Model, ratio=cLeaf/Bbio"
  corr_model_prratio1@long_name = "Correlation of precipitation and cLeaf/Bbio"
  outfile->corr_model_prratio1 = corr_model_prratio1
  delete(corr_model_prratio1)
  delete(temp_corr_model_prratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; PR vs ratio2 (cLeaf/cWood)

  corr_model_prratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prratio2(i,j) = 9.96921e+36
        else
          temp_corr_model_prratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prratio2(i,j) = 9.96921e+36
        else
          temp_corr_model_prratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prratio2(i,j) = temp_corr_model_prratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prratio2(i,j) = temp_corr_model_prratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prratio2(i,j) = 9.96921e+36
        else
          corr_model_prratio2(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prratio2(i,j) = 9.96921e+36
        else
          corr_model_prratio2(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prratio2 = mask(corr_model_prratio2, ismissing(model_pr), False)
  corr_model_prratio2!0   = "lat"
  corr_model_prratio2!1   = "lon"
  corr_model_prratio2&lat = lat
  corr_model_prratio2&lon = lon
  corr_model_prratio2@comment = "Model, ratio=cLeaf/cWood"
  corr_model_prratio2@long_name = "Correlation of precipitation and cLeaf/cWood"
  outfile->corr_model_prratio2 = corr_model_prratio2
  delete(corr_model_prratio2)
  delete(temp_corr_model_prratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; PR vs ratio3 (Bbio/cWood)

  corr_model_prratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_prratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_pr(i,j-10)
      temp_array2(i,j) = model_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_prratio3(i,j) = 9.96921e+36
        else
          temp_corr_model_prratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_prratio3(i,j) = 9.96921e+36
        else
          temp_corr_model_prratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_prratio3(i,j) = temp_corr_model_prratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_prratio3(i,j) = temp_corr_model_prratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_prratio3(i,j) = 9.96921e+36
        else
          corr_model_prratio3(i,j) = pattern_cor( model_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_pr(i-5:i+5,j-5:j+5) ))) then
          corr_model_prratio3(i,j) = 9.96921e+36
        else
          corr_model_prratio3(i,j) = pattern_cor( model_pr(i-5:i+5,j-5:j+5),model_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_prratio3 = mask(corr_model_prratio3, ismissing(model_pr), False)
  corr_model_prratio3!0   = "lat"
  corr_model_prratio3!1   = "lon"
  corr_model_prratio3&lat = lat
  corr_model_prratio3&lon = lon
  corr_model_prratio3@comment = "Model, ratio=Bbio/cWood"
  corr_model_prratio3@long_name = "Correlation of precipitation and Bbio/cWood"
  outfile->corr_model_prratio3 = corr_model_prratio3
  delete(corr_model_prratio3)
  delete(temp_corr_model_prratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; TAS vs RSDS

  corr_model_tasrsds = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tasrsds = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_rsds(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_rsds(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tasrsds(i,j) = 9.96921e+36
        else
          temp_corr_model_tasrsds(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tasrsds(i,j) = 9.96921e+36
        else
          temp_corr_model_tasrsds(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tasrsds(i,j) = temp_corr_model_tasrsds(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tasrsds(i,j) = temp_corr_model_tasrsds(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tasrsds(i,j) = 9.96921e+36
        else
          corr_model_tasrsds(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_rsds(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tasrsds(i,j) = 9.96921e+36
        else
          corr_model_tasrsds(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_rsds(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tasrsds = mask(corr_model_tasrsds, ismissing(model_tas), False)
  corr_model_tasrsds!0   = "lat"
  corr_model_tasrsds!1   = "lon"
  corr_model_tasrsds&lat = lat
  corr_model_tasrsds&lon = lon
  corr_model_tasrsds@comment = "Model"
  corr_model_tasrsds@long_name = "Correlation of surface temperature and shortwave radiation"
  outfile->corr_model_tasrsds = corr_model_tasrsds
  delete(corr_model_tasrsds)
  delete(temp_corr_model_tasrsds)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs Abio

  corr_model_tasAbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tasAbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_Abio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_Abio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tasAbio(i,j) = 9.96921e+36
        else
          temp_corr_model_tasAbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tasAbio(i,j) = 9.96921e+36
        else
          temp_corr_model_tasAbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tasAbio(i,j) = temp_corr_model_tasAbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tasAbio(i,j) = temp_corr_model_tasAbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tasAbio(i,j) = 9.96921e+36
        else
          corr_model_tasAbio(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_Abio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tasAbio(i,j) = 9.96921e+36
        else
          corr_model_tasAbio(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_Abio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tasAbio = mask(corr_model_tasAbio, ismissing(model_tas), False)
  corr_model_tasAbio!0   = "lat"
  corr_model_tasAbio!1   = "lon"
  corr_model_tasAbio&lat = lat
  corr_model_tasAbio&lon = lon
  corr_model_tasAbio@comment = "Model, Tbio-Bbio"
  corr_model_tasAbio@long_name = "Correlation of surface temperature and aboveground biomass"
  outfile->corr_model_tasAbio = corr_model_tasAbio
  delete(corr_model_tasAbio)
  delete(temp_corr_model_tasAbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs Bbio

  corr_model_tasBbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tasBbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_Bbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_Bbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tasBbio(i,j) = 9.96921e+36
        else
          temp_corr_model_tasBbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tasBbio(i,j) = 9.96921e+36
        else
          temp_corr_model_tasBbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tasBbio(i,j) = temp_corr_model_tasBbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tasBbio(i,j) = temp_corr_model_tasBbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tasBbio(i,j) = 9.96921e+36
        else
          corr_model_tasBbio(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_Bbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tasBbio(i,j) = 9.96921e+36
        else
          corr_model_tasBbio(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_Bbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tasBbio = mask(corr_model_tasBbio, ismissing(model_tas), False)
  corr_model_tasBbio!0   = "lat"
  corr_model_tasBbio!1   = "lon"
  corr_model_tasBbio&lat = lat
  corr_model_tasBbio&lon = lon
  corr_model_tasBbio@comment = "Model"
  corr_model_tasBbio@long_name = "Correlation of surface temperature and belowground biomass"
  outfile->corr_model_tasBbio = corr_model_tasBbio
  delete(corr_model_tasBbio)
  delete(temp_corr_model_tasBbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs Tbio

  corr_model_tasTbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tasTbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_Tbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_Tbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tasTbio(i,j) = 9.96921e+36
        else
          temp_corr_model_tasTbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tasTbio(i,j) = 9.96921e+36
        else
          temp_corr_model_tasTbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tasTbio(i,j) = temp_corr_model_tasTbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tasTbio(i,j) = temp_corr_model_tasTbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tasTbio(i,j) = 9.96921e+36
        else
          corr_model_tasTbio(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_Tbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tasTbio(i,j) = 9.96921e+36
        else
          corr_model_tasTbio(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_Tbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tasTbio = mask(corr_model_tasTbio, ismissing(model_tas), False)
  corr_model_tasTbio!0   = "lat"
  corr_model_tasTbio!1   = "lon"
  corr_model_tasTbio&lat = lat
  corr_model_tasTbio&lon = lon
  corr_model_tasTbio@comment = "Model"
  corr_model_tasTbio@long_name = "Correlation of surface temperature and total biomass"
  outfile->corr_model_tasTbio = corr_model_tasTbio
  delete(corr_model_tasTbio)
  delete(temp_corr_model_tasTbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs cLeaf

  corr_model_tascLeaf = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tascLeaf = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_cLeaf(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_cLeaf(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tascLeaf(i,j) = 9.96921e+36
        else
          temp_corr_model_tascLeaf(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tascLeaf(i,j) = 9.96921e+36
        else
          temp_corr_model_tascLeaf(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tascLeaf(i,j) = temp_corr_model_tascLeaf(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tascLeaf(i,j) = temp_corr_model_tascLeaf(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tascLeaf(i,j) = 9.96921e+36
        else
          corr_model_tascLeaf(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tascLeaf(i,j) = 9.96921e+36
        else
          corr_model_tascLeaf(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_cLeaf(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tascLeaf = mask(corr_model_tascLeaf, ismissing(model_tas), False)
  corr_model_tascLeaf!0   = "lat"
  corr_model_tascLeaf!1   = "lon"
  corr_model_tascLeaf&lat = lat
  corr_model_tascLeaf&lon = lon
  corr_model_tascLeaf@comment = "Model"
  corr_model_tascLeaf@long_name = "Correlation of surface temperature and leaf biomass"
  outfile->corr_model_tascLeaf = corr_model_tascLeaf
  delete(corr_model_tascLeaf)
  delete(temp_corr_model_tascLeaf)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs cWood

  corr_model_tascWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tascWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tascWood(i,j) = 9.96921e+36
        else
          temp_corr_model_tascWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tascWood(i,j) = 9.96921e+36
        else
          temp_corr_model_tascWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tascWood(i,j) = temp_corr_model_tascWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tascWood(i,j) = temp_corr_model_tascWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tascWood(i,j) = 9.96921e+36
        else
          corr_model_tascWood(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tascWood(i,j) = 9.96921e+36
        else
          corr_model_tascWood(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tascWood = mask(corr_model_tascWood, ismissing(model_tas), False)
  corr_model_tascWood!0   = "lat"
  corr_model_tascWood!1   = "lon"
  corr_model_tascWood&lat = lat
  corr_model_tascWood&lon = lon
  corr_model_tascWood@comment = "Model"
  corr_model_tascWood@long_name = "Correlation of surface temperature and wood biomass"
  outfile->corr_model_tascWood = corr_model_tascWood
  delete(corr_model_tascWood)
  delete(temp_corr_model_tascWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs ratio0 (Abio/Bbio)

  corr_model_tasratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tasratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tasratio0(i,j) = 9.96921e+36
        else
          temp_corr_model_tasratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tasratio0(i,j) = 9.96921e+36
        else
          temp_corr_model_tasratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tasratio0(i,j) = temp_corr_model_tasratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tasratio0(i,j) = temp_corr_model_tasratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tasratio0(i,j) = 9.96921e+36
        else
          corr_model_tasratio0(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tasratio0(i,j) = 9.96921e+36
        else
          corr_model_tasratio0(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tasratio0 = mask(corr_model_tasratio0, ismissing(model_tas), False)
  corr_model_tasratio0!0   = "lat"
  corr_model_tasratio0!1   = "lon"
  corr_model_tasratio0&lat = lat
  corr_model_tasratio0&lon = lon
  corr_model_tasratio0@comment = "Model, ratio=Abio/Bbio"
  corr_model_tasratio0@long_name = "Correlation of surface temperature and Abio/Bbio"
  outfile->corr_model_tasratio0 = corr_model_tasratio0
  delete(corr_model_tasratio0)
  delete(temp_corr_model_tasratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs ratio1 (cLeaf/Bbio)

  corr_model_tasratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tasratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tasratio1(i,j) = 9.96921e+36
        else
          temp_corr_model_tasratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tasratio1(i,j) = 9.96921e+36
        else
          temp_corr_model_tasratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tasratio1(i,j) = temp_corr_model_tasratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tasratio1(i,j) = temp_corr_model_tasratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tasratio1(i,j) = 9.96921e+36
        else
          corr_model_tasratio1(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tasratio1(i,j) = 9.96921e+36
        else
          corr_model_tasratio1(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tasratio1 = mask(corr_model_tasratio1, ismissing(model_tas), False)
  corr_model_tasratio1!0   = "lat"
  corr_model_tasratio1!1   = "lon"
  corr_model_tasratio1&lat = lat
  corr_model_tasratio1&lon = lon
  corr_model_tasratio1@comment = "Model, ratio=cLeaf/Bbio"
  corr_model_tasratio1@long_name = "Correlation of surface temperature and cLeaf/Bbio"
  outfile->corr_model_tasratio1 = corr_model_tasratio1
  delete(corr_model_tasratio1)
  delete(temp_corr_model_tasratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs ratio2 (cLeaf/cWood)

  corr_model_tasratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tasratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tasratio2(i,j) = 9.96921e+36
        else
          temp_corr_model_tasratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tasratio2(i,j) = 9.96921e+36
        else
          temp_corr_model_tasratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tasratio2(i,j) = temp_corr_model_tasratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tasratio2(i,j) = temp_corr_model_tasratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tasratio2(i,j) = 9.96921e+36
        else
          corr_model_tasratio2(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tasratio2(i,j) = 9.96921e+36
        else
          corr_model_tasratio2(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tasratio2 = mask(corr_model_tasratio2, ismissing(model_tas), False)
  corr_model_tasratio2!0   = "lat"
  corr_model_tasratio2!1   = "lon"
  corr_model_tasratio2&lat = lat
  corr_model_tasratio2&lon = lon
  corr_model_tasratio2@comment = "Model, ratio=cLeaf/cWood"
  corr_model_tasratio2@long_name = "Correlation of surface temperature and cLeaf/cWood"
  outfile->corr_model_tasratio2 = corr_model_tasratio2
  delete(corr_model_tasratio2)
  delete(temp_corr_model_tasratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs ratio3 (Bbio/cWood)

  corr_model_tasratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_tasratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_tas(i,j-10)
      temp_array2(i,j) = model_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_tasratio3(i,j) = 9.96921e+36
        else
          temp_corr_model_tasratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_tasratio3(i,j) = 9.96921e+36
        else
          temp_corr_model_tasratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_tasratio3(i,j) = temp_corr_model_tasratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_tasratio3(i,j) = temp_corr_model_tasratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_tasratio3(i,j) = 9.96921e+36
        else
          corr_model_tasratio3(i,j) = pattern_cor( model_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_tas(i-5:i+5,j-5:j+5) ))) then
          corr_model_tasratio3(i,j) = 9.96921e+36
        else
          corr_model_tasratio3(i,j) = pattern_cor( model_tas(i-5:i+5,j-5:j+5),model_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_tasratio3 = mask(corr_model_tasratio3, ismissing(model_tas), False)
  corr_model_tasratio3!0   = "lat"
  corr_model_tasratio3!1   = "lon"
  corr_model_tasratio3&lat = lat
  corr_model_tasratio3&lon = lon
  corr_model_tasratio3@comment = "Model, ratio=Bbio/cWood"
  corr_model_tasratio3@long_name = "Correlation of surface temperature and Bbio/cWood"
  outfile->corr_model_tasratio3 = corr_model_tasratio3
  delete(corr_model_tasratio3)
  delete(temp_corr_model_tasratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs Abio

  corr_model_rsdsAbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_rsdsAbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_Abio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,j-10)
      temp_array2(i,j) = model_Abio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_rsdsAbio(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsAbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_rsdsAbio(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsAbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_rsdsAbio(i,j) = temp_corr_model_rsdsAbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_rsdsAbio(i,j) = temp_corr_model_rsdsAbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_rsdsAbio(i,j) = 9.96921e+36
        else
          corr_model_rsdsAbio(i,j) = pattern_cor( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_Abio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_model_rsdsAbio(i,j) = 9.96921e+36
        else
          corr_model_rsdsAbio(i,j) = pattern_cor( model_rsds(i-5:i+5,j-5:j+5),model_Abio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_rsdsAbio = mask(corr_model_rsdsAbio, ismissing(model_rsds), False)
  corr_model_rsdsAbio!0   = "lat"
  corr_model_rsdsAbio!1   = "lon"
  corr_model_rsdsAbio&lat = lat
  corr_model_rsdsAbio&lon = lon
  corr_model_rsdsAbio@comment = "Model, Tbio-Bbio"
  corr_model_rsdsAbio@long_name = "Correlation of shortwave radiation and aboveground biomass"
  outfile->corr_model_rsdsAbio = corr_model_rsdsAbio
  delete(corr_model_rsdsAbio)
  delete(temp_corr_model_rsdsAbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs Bbio

  corr_model_rsdsBbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_rsdsBbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_Bbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,j-10)
      temp_array2(i,j) = model_Bbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_rsdsBbio(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsBbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_rsdsBbio(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsBbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_rsdsBbio(i,j) = temp_corr_model_rsdsBbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_rsdsBbio(i,j) = temp_corr_model_rsdsBbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_rsdsBbio(i,j) = 9.96921e+36
        else
          corr_model_rsdsBbio(i,j) = pattern_cor( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_Bbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_model_rsdsBbio(i,j) = 9.96921e+36
        else
          corr_model_rsdsBbio(i,j) = pattern_cor( model_rsds(i-5:i+5,j-5:j+5),model_Bbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_rsdsBbio = mask(corr_model_rsdsBbio, ismissing(model_rsds), False)
  corr_model_rsdsBbio!0   = "lat"
  corr_model_rsdsBbio!1   = "lon"
  corr_model_rsdsBbio&lat = lat
  corr_model_rsdsBbio&lon = lon
  corr_model_rsdsBbio@comment = "Model"
  corr_model_rsdsBbio@long_name = "Correlation of shortwave radiation and belowground biomass"
  outfile->corr_model_rsdsBbio = corr_model_rsdsBbio
  delete(corr_model_rsdsBbio)
  delete(temp_corr_model_rsdsBbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs Tbio

  corr_model_rsdsTbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_rsdsTbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_Tbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,j-10)
      temp_array2(i,j) = model_Tbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_rsdsTbio(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsTbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_rsdsTbio(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsTbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_rsdsTbio(i,j) = temp_corr_model_rsdsTbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_rsdsTbio(i,j) = temp_corr_model_rsdsTbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_rsdsTbio(i,j) = 9.96921e+36
        else
          corr_model_rsdsTbio(i,j) = pattern_cor( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_Tbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_model_rsdsTbio(i,j) = 9.96921e+36
        else
          corr_model_rsdsTbio(i,j) = pattern_cor( model_rsds(i-5:i+5,j-5:j+5),model_Tbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_rsdsTbio = mask(corr_model_rsdsTbio, ismissing(model_rsds), False)
  corr_model_rsdsTbio!0   = "lat"
  corr_model_rsdsTbio!1   = "lon"
  corr_model_rsdsTbio&lat = lat
  corr_model_rsdsTbio&lon = lon
  corr_model_rsdsTbio@comment = "Model"
  corr_model_rsdsTbio@long_name = "Correlation of shortwave radiation and total biomass"
  outfile->corr_model_rsdsTbio = corr_model_rsdsTbio
  delete(corr_model_rsdsTbio)
  delete(temp_corr_model_rsdsTbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs cLeaf

  corr_model_rsdscLeaf = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_rsdscLeaf = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_cLeaf(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,j-10)
      temp_array2(i,j) = model_cLeaf(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_rsdscLeaf(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdscLeaf(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_rsdscLeaf(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdscLeaf(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_rsdscLeaf(i,j) = temp_corr_model_rsdscLeaf(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_rsdscLeaf(i,j) = temp_corr_model_rsdscLeaf(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_rsdscLeaf(i,j) = 9.96921e+36
        else
          corr_model_rsdscLeaf(i,j) = pattern_cor( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_model_rsdscLeaf(i,j) = 9.96921e+36
        else
          corr_model_rsdscLeaf(i,j) = pattern_cor( model_rsds(i-5:i+5,j-5:j+5),model_cLeaf(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_rsdscLeaf = mask(corr_model_rsdscLeaf, ismissing(model_rsds), False)
  corr_model_rsdscLeaf!0   = "lat"
  corr_model_rsdscLeaf!1   = "lon"
  corr_model_rsdscLeaf&lat = lat
  corr_model_rsdscLeaf&lon = lon
  corr_model_rsdscLeaf@comment = "Model"
  corr_model_rsdscLeaf@long_name = "Correlation of shortwave radiation and leaf biomass"
  outfile->corr_model_rsdscLeaf = corr_model_rsdscLeaf
  delete(corr_model_rsdscLeaf)
  delete(temp_corr_model_rsdscLeaf)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs cWood

  corr_model_rsdscWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_rsdscWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,j-10)
      temp_array2(i,j) = model_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_rsdscWood(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdscWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_rsdscWood(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdscWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_rsdscWood(i,j) = temp_corr_model_rsdscWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_rsdscWood(i,j) = temp_corr_model_rsdscWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_rsdscWood(i,j) = 9.96921e+36
        else
          corr_model_rsdscWood(i,j) = pattern_cor( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_model_rsdscWood(i,j) = 9.96921e+36
        else
          corr_model_rsdscWood(i,j) = pattern_cor( model_rsds(i-5:i+5,j-5:j+5),model_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_rsdscWood = mask(corr_model_rsdscWood, ismissing(model_rsds), False)
  corr_model_rsdscWood!0   = "lat"
  corr_model_rsdscWood!1   = "lon"
  corr_model_rsdscWood&lat = lat
  corr_model_rsdscWood&lon = lon
  corr_model_rsdscWood@comment = "Model"
  corr_model_rsdscWood@long_name = "Correlation of shortwave radiation and wood biomass"
  outfile->corr_model_rsdscWood = corr_model_rsdscWood
  delete(corr_model_rsdscWood)
  delete(temp_corr_model_rsdscWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs ratio0 (Abio/Bbio)

  corr_model_rsdsratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_rsdsratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,j-10)
      temp_array2(i,j) = model_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_rsdsratio0(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_rsdsratio0(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_rsdsratio0(i,j) = temp_corr_model_rsdsratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_rsdsratio0(i,j) = temp_corr_model_rsdsratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_rsdsratio0(i,j) = 9.96921e+36
        else
          corr_model_rsdsratio0(i,j) = pattern_cor( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_model_rsdsratio0(i,j) = 9.96921e+36
        else
          corr_model_rsdsratio0(i,j) = pattern_cor( model_rsds(i-5:i+5,j-5:j+5),model_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_rsdsratio0 = mask(corr_model_rsdsratio0, ismissing(model_rsds), False)
  corr_model_rsdsratio0!0   = "lat"
  corr_model_rsdsratio0!1   = "lon"
  corr_model_rsdsratio0&lat = lat
  corr_model_rsdsratio0&lon = lon
  corr_model_rsdsratio0@comment = "Model, ratio=Abio/Bbio"
  corr_model_rsdsratio0@long_name = "Correlation of shortwave radiation and Abio/Bbio"
  outfile->corr_model_rsdsratio0 = corr_model_rsdsratio0
  delete(corr_model_rsdsratio0)
  delete(temp_corr_model_rsdsratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs ratio1 (cLeaf/Bbio)

  corr_model_rsdsratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_rsdsratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,j-10)
      temp_array2(i,j) = model_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_rsdsratio1(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_rsdsratio1(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_rsdsratio1(i,j) = temp_corr_model_rsdsratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_rsdsratio1(i,j) = temp_corr_model_rsdsratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_rsdsratio1(i,j) = 9.96921e+36
        else
          corr_model_rsdsratio1(i,j) = pattern_cor( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_model_rsdsratio1(i,j) = 9.96921e+36
        else
          corr_model_rsdsratio1(i,j) = pattern_cor( model_rsds(i-5:i+5,j-5:j+5),model_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_rsdsratio1 = mask(corr_model_rsdsratio1, ismissing(model_rsds), False)
  corr_model_rsdsratio1!0   = "lat"
  corr_model_rsdsratio1!1   = "lon"
  corr_model_rsdsratio1&lat = lat
  corr_model_rsdsratio1&lon = lon
  corr_model_rsdsratio1@comment = "Model, ratio=cLeaf/Bbio"
  corr_model_rsdsratio1@long_name = "Correlation of shortwave radiation and cLeaf/Bbio"
  outfile->corr_model_rsdsratio1 = corr_model_rsdsratio1
  delete(corr_model_rsdsratio1)
  delete(temp_corr_model_rsdsratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs ratio2 (cLeaf/cWood)

  corr_model_rsdsratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_rsdsratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,j-10)
      temp_array2(i,j) = model_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_rsdsratio2(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_rsdsratio2(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_rsdsratio2(i,j) = temp_corr_model_rsdsratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_rsdsratio2(i,j) = temp_corr_model_rsdsratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_rsdsratio2(i,j) = 9.96921e+36
        else
          corr_model_rsdsratio2(i,j) = pattern_cor( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_model_rsdsratio2(i,j) = 9.96921e+36
        else
          corr_model_rsdsratio2(i,j) = pattern_cor( model_rsds(i-5:i+5,j-5:j+5),model_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_rsdsratio2 = mask(corr_model_rsdsratio2, ismissing(model_rsds), False)
  corr_model_rsdsratio2!0   = "lat"
  corr_model_rsdsratio2!1   = "lon"
  corr_model_rsdsratio2&lat = lat
  corr_model_rsdsratio2&lon = lon
  corr_model_rsdsratio2@comment = "Model, ratio=cLeaf/cWood"
  corr_model_rsdsratio2@long_name = "Correlation of shortwave radiation and cLeaf/cWood"
  outfile->corr_model_rsdsratio2 = corr_model_rsdsratio2
  delete(corr_model_rsdsratio2)
  delete(temp_corr_model_rsdsratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs ratio3 (Bbio/cWood)

  corr_model_rsdsratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_model_rsdsratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = model_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = model_rsds(i,j-10)
      temp_array2(i,j) = model_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_model_rsdsratio3(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_model_rsdsratio3(i,j) = 9.96921e+36
        else
          temp_corr_model_rsdsratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_model_rsdsratio3(i,j) = temp_corr_model_rsdsratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_model_rsdsratio3(i,j) = temp_corr_model_rsdsratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_model_rsdsratio3(i,j) = 9.96921e+36
        else
          corr_model_rsdsratio3(i,j) = pattern_cor( model_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               model_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( model_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_model_rsdsratio3(i,j) = 9.96921e+36
        else
          corr_model_rsdsratio3(i,j) = pattern_cor( model_rsds(i-5:i+5,j-5:j+5),model_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_model_rsdsratio3 = mask(corr_model_rsdsratio3, ismissing(model_rsds), False)
  corr_model_rsdsratio3!0   = "lat"
  corr_model_rsdsratio3!1   = "lon"
  corr_model_rsdsratio3&lat = lat
  corr_model_rsdsratio3&lon = lon
  corr_model_rsdsratio3@comment = "Model, ratio=Bbio/cWood"
  corr_model_rsdsratio3@long_name = "Correlation of shortwave radiation and Bbio/cWood"
  outfile->corr_model_rsdsratio3 = corr_model_rsdsratio3
  delete(corr_model_rsdsratio3)
  delete(temp_corr_model_rsdsratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------

;;; Delete variables

  delete(model_pr)
  delete(model_rsds)
  delete(model_tas)
  delete(temp_model_Bbio)
  delete(model_Bbio)
  delete(model_Tbio)
  delete(model_cLeaf)
  delete(temp_model_cWood)
  delete(model_cWood)
  delete(model_Abio)
  delete(model_ratio0)
  delete(model_ratio1)
  delete(model_ratio2)
  delete(model_ratio3)

  delete(temp_lat)
  delete(lat)
  delete(temp_lon)
  delete(lon)
  delete(outfile)
  delete(infile)
end do
end do

end
exit

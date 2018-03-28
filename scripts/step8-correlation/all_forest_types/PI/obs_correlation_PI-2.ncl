;==========================================================
; Compute statistics of the BGI dataset
; For preidustrial time (1861-1885)
;==========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  obsnm = (/"BNU-ESM", "HadGEM2-CC", "HadGEM2-ES", "IPSL-CM5A-LR", \
            "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC-ESM-CHEM", "MIROC-ESM"/)
;  forest_frac = (/0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)
  forest_frac = (/0, 0.2, 0.5, 0.8/)

;;; Start looping
do fc=0,dimsizes(forest_frac)-1
  indir = "../../../../postanalysis/"
  indir1= indir+"input_data/PI/frac"+forest_frac(fc)+"/"
  do z=0,dimsizes(obsnm)-1
    print (obsnm(z)+" ")

    fname = "mask_climate_all_"+obsnm(z)+".nc"
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

;;; Biomass
      obs_Bbio = tofloat(infile->mask_obs_biome_dens_Bbio)
      obs_Tbio = tofloat(infile->mask_obs_biome_dens_Tbio)
      obs_cLeaf= tofloat(infile->mask_obs_biome_dens_cLeaf)
      obs_cWood= tofloat(infile->mask_obs_biome_dens_cWood)
      obs_Abio = obs_Tbio - obs_Bbio

      temp_obs_Bbio = mask(obs_Bbio, obs_Bbio.le.0. , False)
      temp_obs_cWood = mask(obs_cWood, obs_cWood.le.0. , False)

      obs_ratio0 = obs_Abio/temp_obs_Bbio
      obs_ratio1 = obs_cLeaf/temp_obs_cWood
      obs_ratio2 = obs_cLeaf/temp_obs_Bbio
      obs_ratio3 = obs_Bbio/temp_obs_cWood

;;; Output file
      outdir  = indir+"correlation/all_forest_types/PI/frac"+forest_frac(fc)+"/"
      fout    = outdir+"correlation-"+obsnm(z)+"-OBS-2.nc"
      if (isfilepresent(fout))
        print("Warning: '"+fout+"' exists.")
        print("Will remove it.")
        system ("/bin/rm -f "+fout)
      end if
      outfile = addfile (fout,"c")
      outfile->lat = lat
      outfile->lon = lon

;----------------------------------------------------------
;;; Using 11-grid moving window, note at north-south bounday,
;;; grids are less than 5
;;; Boundary at Southern hemisphere is not modified!!

;;; Tbio vs Abio

  corr_obs_TbioAbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_TbioAbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,j-10)
      temp_array2(i,j) = obs_Abio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_TbioAbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_TbioAbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_TbioAbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_TbioAbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_TbioAbio(i,j) = temp_corr_obs_TbioAbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_TbioAbio(i,j) = temp_corr_obs_TbioAbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_TbioAbio(i,j) = 9.96921e+36
        else
          corr_obs_TbioAbio(i,j) = pattern_cor( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_TbioAbio(i,j) = 9.96921e+36
        else
          corr_obs_TbioAbio(i,j) = pattern_cor( obs_Tbio(i-5:i+5,j-5:j+5),obs_Abio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_TbioAbio = mask(corr_obs_TbioAbio, ismissing(obs_Tbio), False)
  corr_obs_TbioAbio!0   = "lat"
  corr_obs_TbioAbio!1   = "lon"
  corr_obs_TbioAbio&lat = lat
  corr_obs_TbioAbio&lon = lon
  corr_obs_TbioAbio@comment = "Observation, Tbio-Bbio, preindustrial"
  corr_obs_TbioAbio@long_name = "Correlation of total biomass and aboveground biomass"
  outfile->corr_obs_TbioAbio = corr_obs_TbioAbio
  delete(corr_obs_TbioAbio)
  delete(temp_corr_obs_TbioAbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Tbio vs Bbio

  corr_obs_TbioBbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_TbioBbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,j-10)
      temp_array2(i,j) = obs_Bbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_TbioBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_TbioBbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_TbioBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_TbioBbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_TbioBbio(i,j) = temp_corr_obs_TbioBbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_TbioBbio(i,j) = temp_corr_obs_TbioBbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_TbioBbio(i,j) = 9.96921e+36
        else
          corr_obs_TbioBbio(i,j) = pattern_cor( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_TbioBbio(i,j) = 9.96921e+36
        else
          corr_obs_TbioBbio(i,j) = pattern_cor( obs_Tbio(i-5:i+5,j-5:j+5),obs_Bbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_TbioBbio = mask(corr_obs_TbioBbio, ismissing(obs_Tbio), False)
  corr_obs_TbioBbio!0   = "lat"
  corr_obs_TbioBbio!1   = "lon"
  corr_obs_TbioBbio&lat = lat
  corr_obs_TbioBbio&lon = lon
  corr_obs_TbioBbio@comment = "Observation, preindustrial"
  corr_obs_TbioBbio@long_name = "Correlation of total biomass and belowground biomass"
  outfile->corr_obs_TbioBbio = corr_obs_TbioBbio
  delete(corr_obs_TbioBbio)
  delete(temp_corr_obs_TbioBbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Tbio vs cLeaf

  corr_obs_TbiocLeaf = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_TbiocLeaf = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,j-10)
      temp_array2(i,j) = obs_cLeaf(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_TbiocLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_TbiocLeaf(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_TbiocLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_TbiocLeaf(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_TbiocLeaf(i,j) = temp_corr_obs_TbiocLeaf(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_TbiocLeaf(i,j) = temp_corr_obs_TbiocLeaf(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_TbiocLeaf(i,j) = 9.96921e+36
        else
          corr_obs_TbiocLeaf(i,j) = pattern_cor( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_TbiocLeaf(i,j) = 9.96921e+36
        else
          corr_obs_TbiocLeaf(i,j) = pattern_cor( obs_Tbio(i-5:i+5,j-5:j+5),obs_cLeaf(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_TbiocLeaf = mask(corr_obs_TbiocLeaf, ismissing(obs_Tbio), False)
  corr_obs_TbiocLeaf!0   = "lat"
  corr_obs_TbiocLeaf!1   = "lon"
  corr_obs_TbiocLeaf&lat = lat
  corr_obs_TbiocLeaf&lon = lon
  corr_obs_TbiocLeaf@comment = "Observation, preindustrial"
  corr_obs_TbiocLeaf@long_name = "Correlation of total biomass and leaf biomass"
  outfile->corr_obs_TbiocLeaf = corr_obs_TbiocLeaf
  delete(corr_obs_TbiocLeaf)
  delete(temp_corr_obs_TbiocLeaf)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Tbio vs cWood

  corr_obs_TbiocWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_TbiocWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,j-10)
      temp_array2(i,j) = obs_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_TbiocWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_TbiocWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_TbiocWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_TbiocWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_TbiocWood(i,j) = temp_corr_obs_TbiocWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_TbiocWood(i,j) = temp_corr_obs_TbiocWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_TbiocWood(i,j) = 9.96921e+36
        else
          corr_obs_TbiocWood(i,j) = pattern_cor( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_TbiocWood(i,j) = 9.96921e+36
        else
          corr_obs_TbiocWood(i,j) = pattern_cor( obs_Tbio(i-5:i+5,j-5:j+5),obs_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_TbiocWood = mask(corr_obs_TbiocWood, ismissing(obs_Tbio), False)
  corr_obs_TbiocWood!0   = "lat"
  corr_obs_TbiocWood!1   = "lon"
  corr_obs_TbiocWood&lat = lat
  corr_obs_TbiocWood&lon = lon
  corr_obs_TbiocWood@comment = "Observation, preindustrial"
  corr_obs_TbiocWood@long_name = "Correlation of total biomass and wood biomass"
  outfile->corr_obs_TbiocWood = corr_obs_TbiocWood
  delete(corr_obs_TbiocWood)
  delete(temp_corr_obs_TbiocWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Tbio vs ratio0 (Abio/Bbio)

  corr_obs_Tbioratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Tbioratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,j-10)
      temp_array2(i,j) = obs_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Tbioratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_Tbioratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Tbioratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_Tbioratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Tbioratio0(i,j) = temp_corr_obs_Tbioratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Tbioratio0(i,j) = temp_corr_obs_Tbioratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Tbioratio0(i,j) = 9.96921e+36
        else
          corr_obs_Tbioratio0(i,j) = pattern_cor( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Tbioratio0(i,j) = 9.96921e+36
        else
          corr_obs_Tbioratio0(i,j) = pattern_cor( obs_Tbio(i-5:i+5,j-5:j+5),obs_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Tbioratio0 = mask(corr_obs_Tbioratio0, ismissing(obs_Tbio), False)
  corr_obs_Tbioratio0!0   = "lat"
  corr_obs_Tbioratio0!1   = "lon"
  corr_obs_Tbioratio0&lat = lat
  corr_obs_Tbioratio0&lon = lon
  corr_obs_Tbioratio0@comment = "Observation, ratio=Abio/Bbio, preindustrial"
  corr_obs_Tbioratio0@long_name = "Correlation of total biomass and Abio/Bbio"
  outfile->corr_obs_Tbioratio0 = corr_obs_Tbioratio0
  delete(corr_obs_Tbioratio0)
  delete(temp_corr_obs_Tbioratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Tbio vs ratio1 (cLeaf/Bbio)

  corr_obs_Tbioratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Tbioratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,j-10)
      temp_array2(i,j) = obs_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Tbioratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_Tbioratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Tbioratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_Tbioratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Tbioratio1(i,j) = temp_corr_obs_Tbioratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Tbioratio1(i,j) = temp_corr_obs_Tbioratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Tbioratio1(i,j) = 9.96921e+36
        else
          corr_obs_Tbioratio1(i,j) = pattern_cor( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Tbioratio1(i,j) = 9.96921e+36
        else
          corr_obs_Tbioratio1(i,j) = pattern_cor( obs_Tbio(i-5:i+5,j-5:j+5),obs_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Tbioratio1 = mask(corr_obs_Tbioratio1, ismissing(obs_Tbio), False)
  corr_obs_Tbioratio1!0   = "lat"
  corr_obs_Tbioratio1!1   = "lon"
  corr_obs_Tbioratio1&lat = lat
  corr_obs_Tbioratio1&lon = lon
  corr_obs_Tbioratio1@comment = "Observation, ratio=cLeaf/Bbio, preindustrial"
  corr_obs_Tbioratio1@long_name = "Correlation of total biomass and cLeaf/Bbio"
  outfile->corr_obs_Tbioratio1 = corr_obs_Tbioratio1
  delete(corr_obs_Tbioratio1)
  delete(temp_corr_obs_Tbioratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Tbio vs ratio2 (cLeaf/cWood)

  corr_obs_Tbioratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Tbioratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,j-10)
      temp_array2(i,j) = obs_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Tbioratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_Tbioratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Tbioratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_Tbioratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Tbioratio2(i,j) = temp_corr_obs_Tbioratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Tbioratio2(i,j) = temp_corr_obs_Tbioratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Tbioratio2(i,j) = 9.96921e+36
        else
          corr_obs_Tbioratio2(i,j) = pattern_cor( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Tbioratio2(i,j) = 9.96921e+36
        else
          corr_obs_Tbioratio2(i,j) = pattern_cor( obs_Tbio(i-5:i+5,j-5:j+5),obs_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Tbioratio2 = mask(corr_obs_Tbioratio2, ismissing(obs_Tbio), False)
  corr_obs_Tbioratio2!0   = "lat"
  corr_obs_Tbioratio2!1   = "lon"
  corr_obs_Tbioratio2&lat = lat
  corr_obs_Tbioratio2&lon = lon
  corr_obs_Tbioratio2@comment = "Observation, ratio=cLeaf/cWood, preindustrial"
  corr_obs_Tbioratio2@long_name = "Correlation of total biomass and cLeaf/cWood"
  outfile->corr_obs_Tbioratio2 = corr_obs_Tbioratio2
  delete(corr_obs_Tbioratio2)
  delete(temp_corr_obs_Tbioratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Tbio vs ratio3 (Bbio/cWood)

  corr_obs_Tbioratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Tbioratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Tbio(i,j-10)
      temp_array2(i,j) = obs_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Tbioratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_Tbioratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Tbioratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_Tbioratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Tbioratio3(i,j) = temp_corr_obs_Tbioratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Tbioratio3(i,j) = temp_corr_obs_Tbioratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Tbioratio3(i,j) = 9.96921e+36
        else
          corr_obs_Tbioratio3(i,j) = pattern_cor( obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Tbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Tbioratio3(i,j) = 9.96921e+36
        else
          corr_obs_Tbioratio3(i,j) = pattern_cor( obs_Tbio(i-5:i+5,j-5:j+5),obs_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Tbioratio3 = mask(corr_obs_Tbioratio3, ismissing(obs_Tbio), False)
  corr_obs_Tbioratio3!0   = "lat"
  corr_obs_Tbioratio3!1   = "lon"
  corr_obs_Tbioratio3&lat = lat
  corr_obs_Tbioratio3&lon = lon
  corr_obs_Tbioratio3@comment = "Observation, ratio=Bbio/cWood, preindustrial"
  corr_obs_Tbioratio3@long_name = "Correlation of total biomass and Bbio/cWood"
  outfile->corr_obs_Tbioratio3 = corr_obs_Tbioratio3
  delete(corr_obs_Tbioratio3)
  delete(temp_corr_obs_Tbioratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Abio vs Bbio

  corr_obs_AbioBbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_AbioBbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,j-10)
      temp_array2(i,j) = obs_Bbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_AbioBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_AbioBbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_AbioBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_AbioBbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_AbioBbio(i,j) = temp_corr_obs_AbioBbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_AbioBbio(i,j) = temp_corr_obs_AbioBbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_AbioBbio(i,j) = 9.96921e+36
        else
          corr_obs_AbioBbio(i,j) = pattern_cor( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_AbioBbio(i,j) = 9.96921e+36
        else
          corr_obs_AbioBbio(i,j) = pattern_cor( obs_Abio(i-5:i+5,j-5:j+5),obs_Bbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_AbioBbio = mask(corr_obs_AbioBbio, ismissing(obs_Abio), False)
  corr_obs_AbioBbio!0   = "lat"
  corr_obs_AbioBbio!1   = "lon"
  corr_obs_AbioBbio&lat = lat
  corr_obs_AbioBbio&lon = lon
  corr_obs_AbioBbio@comment = "Observation, preindustrial"
  corr_obs_AbioBbio@long_name = "Correlation of aboveground biomass and belowground biomass"
  outfile->corr_obs_AbioBbio = corr_obs_AbioBbio
  delete(corr_obs_AbioBbio)
  delete(temp_corr_obs_AbioBbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Abio vs cLeaf

  corr_obs_AbiocLeaf = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_AbiocLeaf = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,j-10)
      temp_array2(i,j) = obs_cLeaf(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_AbiocLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_AbiocLeaf(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_AbiocLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_AbiocLeaf(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_AbiocLeaf(i,j) = temp_corr_obs_AbiocLeaf(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_AbiocLeaf(i,j) = temp_corr_obs_AbiocLeaf(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_AbiocLeaf(i,j) = 9.96921e+36
        else
          corr_obs_AbiocLeaf(i,j) = pattern_cor( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_AbiocLeaf(i,j) = 9.96921e+36
        else
          corr_obs_AbiocLeaf(i,j) = pattern_cor( obs_Abio(i-5:i+5,j-5:j+5),obs_cLeaf(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_AbiocLeaf = mask(corr_obs_AbiocLeaf, ismissing(obs_Abio), False)
  corr_obs_AbiocLeaf!0   = "lat"
  corr_obs_AbiocLeaf!1   = "lon"
  corr_obs_AbiocLeaf&lat = lat
  corr_obs_AbiocLeaf&lon = lon
  corr_obs_AbiocLeaf@comment = "Observation, preindustrial"
  corr_obs_AbiocLeaf@long_name = "Correlation of aboveground biomass and leaf biomass"
  outfile->corr_obs_AbiocLeaf = corr_obs_AbiocLeaf
  delete(corr_obs_AbiocLeaf)
  delete(temp_corr_obs_AbiocLeaf)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Abio vs cWood

  corr_obs_AbiocWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_AbiocWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,j-10)
      temp_array2(i,j) = obs_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_AbiocWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_AbiocWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_AbiocWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_AbiocWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_AbiocWood(i,j) = temp_corr_obs_AbiocWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_AbiocWood(i,j) = temp_corr_obs_AbiocWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_AbiocWood(i,j) = 9.96921e+36
        else
          corr_obs_AbiocWood(i,j) = pattern_cor( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_AbiocWood(i,j) = 9.96921e+36
        else
          corr_obs_AbiocWood(i,j) = pattern_cor( obs_Abio(i-5:i+5,j-5:j+5),obs_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_AbiocWood = mask(corr_obs_AbiocWood, ismissing(obs_Abio), False)
  corr_obs_AbiocWood!0   = "lat"
  corr_obs_AbiocWood!1   = "lon"
  corr_obs_AbiocWood&lat = lat
  corr_obs_AbiocWood&lon = lon
  corr_obs_AbiocWood@comment = "Observation, preindustrial"
  corr_obs_AbiocWood@long_name = "Correlation of aboveground biomass and wood biomass"
  outfile->corr_obs_AbiocWood = corr_obs_AbiocWood
  delete(corr_obs_AbiocWood)
  delete(temp_corr_obs_AbiocWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Abio vs ratio0 (Abio/Bbio)

  corr_obs_Abioratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Abioratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,j-10)
      temp_array2(i,j) = obs_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Abioratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_Abioratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Abioratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_Abioratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Abioratio0(i,j) = temp_corr_obs_Abioratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Abioratio0(i,j) = temp_corr_obs_Abioratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Abioratio0(i,j) = 9.96921e+36
        else
          corr_obs_Abioratio0(i,j) = pattern_cor( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Abioratio0(i,j) = 9.96921e+36
        else
          corr_obs_Abioratio0(i,j) = pattern_cor( obs_Abio(i-5:i+5,j-5:j+5),obs_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Abioratio0 = mask(corr_obs_Abioratio0, ismissing(obs_Abio), False)
  corr_obs_Abioratio0!0   = "lat"
  corr_obs_Abioratio0!1   = "lon"
  corr_obs_Abioratio0&lat = lat
  corr_obs_Abioratio0&lon = lon
  corr_obs_Abioratio0@comment = "Observation, ratio=Abio/Bbio, preindustrial"
  corr_obs_Abioratio0@long_name = "Correlation of aboveground biomass and Abio/Bbio"
  outfile->corr_obs_Abioratio0 = corr_obs_Abioratio0
  delete(corr_obs_Abioratio0)
  delete(temp_corr_obs_Abioratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Abio vs ratio1 (cLeaf/Bbio)

  corr_obs_Abioratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Abioratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,j-10)
      temp_array2(i,j) = obs_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Abioratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_Abioratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Abioratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_Abioratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Abioratio1(i,j) = temp_corr_obs_Abioratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Abioratio1(i,j) = temp_corr_obs_Abioratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Abioratio1(i,j) = 9.96921e+36
        else
          corr_obs_Abioratio1(i,j) = pattern_cor( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Abioratio1(i,j) = 9.96921e+36
        else
          corr_obs_Abioratio1(i,j) = pattern_cor( obs_Abio(i-5:i+5,j-5:j+5),obs_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Abioratio1 = mask(corr_obs_Abioratio1, ismissing(obs_Abio), False)
  corr_obs_Abioratio1!0   = "lat"
  corr_obs_Abioratio1!1   = "lon"
  corr_obs_Abioratio1&lat = lat
  corr_obs_Abioratio1&lon = lon
  corr_obs_Abioratio1@comment = "Observation, ratio=cLeaf/Bbio, preindustrial"
  corr_obs_Abioratio1@long_name = "Correlation of aboveground biomass and cLeaf/Bbio"
  outfile->corr_obs_Abioratio1 = corr_obs_Abioratio1
  delete(corr_obs_Abioratio1)
  delete(temp_corr_obs_Abioratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Abio vs ratio2 (cLeaf/cWood)

  corr_obs_Abioratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Abioratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,j-10)
      temp_array2(i,j) = obs_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Abioratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_Abioratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Abioratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_Abioratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Abioratio2(i,j) = temp_corr_obs_Abioratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Abioratio2(i,j) = temp_corr_obs_Abioratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Abioratio2(i,j) = 9.96921e+36
        else
          corr_obs_Abioratio2(i,j) = pattern_cor( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Abioratio2(i,j) = 9.96921e+36
        else
          corr_obs_Abioratio2(i,j) = pattern_cor( obs_Abio(i-5:i+5,j-5:j+5),obs_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Abioratio2 = mask(corr_obs_Abioratio2, ismissing(obs_Abio), False)
  corr_obs_Abioratio2!0   = "lat"
  corr_obs_Abioratio2!1   = "lon"
  corr_obs_Abioratio2&lat = lat
  corr_obs_Abioratio2&lon = lon
  corr_obs_Abioratio2@comment = "Observation, ratio=cLeaf/cWood, preindustrial"
  corr_obs_Abioratio2@long_name = "Correlation of aboveground biomass and cLeaf/cWood"
  outfile->corr_obs_Abioratio2 = corr_obs_Abioratio2
  delete(corr_obs_Abioratio2)
  delete(temp_corr_obs_Abioratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Abio vs ratio3 (Bbio/cWood)

  corr_obs_Abioratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Abioratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Abio(i,j-10)
      temp_array2(i,j) = obs_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Abioratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_Abioratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Abioratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_Abioratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Abioratio3(i,j) = temp_corr_obs_Abioratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Abioratio3(i,j) = temp_corr_obs_Abioratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Abioratio3(i,j) = 9.96921e+36
        else
          corr_obs_Abioratio3(i,j) = pattern_cor( obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Abio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Abioratio3(i,j) = 9.96921e+36
        else
          corr_obs_Abioratio3(i,j) = pattern_cor( obs_Abio(i-5:i+5,j-5:j+5),obs_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Abioratio3 = mask(corr_obs_Abioratio3, ismissing(obs_Abio), False)
  corr_obs_Abioratio3!0   = "lat"
  corr_obs_Abioratio3!1   = "lon"
  corr_obs_Abioratio3&lat = lat
  corr_obs_Abioratio3&lon = lon
  corr_obs_Abioratio3@comment = "Observation, ratio=Bbio/cWood, preindustrial"
  corr_obs_Abioratio3@long_name = "Correlation of aboveground biomass and Bbio/cWood"
  outfile->corr_obs_Abioratio3 = corr_obs_Abioratio3
  delete(corr_obs_Abioratio3)
  delete(temp_corr_obs_Abioratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Bbio vs cLeaf

  corr_obs_BbiocLeaf = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_BbiocLeaf = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,j-10)
      temp_array2(i,j) = obs_cLeaf(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_BbiocLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_BbiocLeaf(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_BbiocLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_BbiocLeaf(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_BbiocLeaf(i,j) = temp_corr_obs_BbiocLeaf(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_BbiocLeaf(i,j) = temp_corr_obs_BbiocLeaf(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_BbiocLeaf(i,j) = 9.96921e+36
        else
          corr_obs_BbiocLeaf(i,j) = pattern_cor( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_BbiocLeaf(i,j) = 9.96921e+36
        else
          corr_obs_BbiocLeaf(i,j) = pattern_cor( obs_Bbio(i-5:i+5,j-5:j+5),obs_cLeaf(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_BbiocLeaf = mask(corr_obs_BbiocLeaf, ismissing(obs_Bbio), False)
  corr_obs_BbiocLeaf!0   = "lat"
  corr_obs_BbiocLeaf!1   = "lon"
  corr_obs_BbiocLeaf&lat = lat
  corr_obs_BbiocLeaf&lon = lon
  corr_obs_BbiocLeaf@comment = "Observation, preindustrial"
  corr_obs_BbiocLeaf@long_name = "Correlation of belowground biomass and leaf biomass"
  outfile->corr_obs_BbiocLeaf = corr_obs_BbiocLeaf
  delete(corr_obs_BbiocLeaf)
  delete(temp_corr_obs_BbiocLeaf)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Bbio vs cWood

  corr_obs_BbiocWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_BbiocWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,j-10)
      temp_array2(i,j) = obs_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_BbiocWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_BbiocWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_BbiocWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_BbiocWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_BbiocWood(i,j) = temp_corr_obs_BbiocWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_BbiocWood(i,j) = temp_corr_obs_BbiocWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_BbiocWood(i,j) = 9.96921e+36
        else
          corr_obs_BbiocWood(i,j) = pattern_cor( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_BbiocWood(i,j) = 9.96921e+36
        else
          corr_obs_BbiocWood(i,j) = pattern_cor( obs_Bbio(i-5:i+5,j-5:j+5),obs_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_BbiocWood = mask(corr_obs_BbiocWood, ismissing(obs_Bbio), False)
  corr_obs_BbiocWood!0   = "lat"
  corr_obs_BbiocWood!1   = "lon"
  corr_obs_BbiocWood&lat = lat
  corr_obs_BbiocWood&lon = lon
  corr_obs_BbiocWood@comment = "Observation, preindustrial"
  corr_obs_BbiocWood@long_name = "Correlation of belowground biomass and wood biomass"
  outfile->corr_obs_BbiocWood = corr_obs_BbiocWood
  delete(corr_obs_BbiocWood)
  delete(temp_corr_obs_BbiocWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Bbio vs ratio0 (Abio/Bbio)

  corr_obs_Bbioratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Bbioratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,j-10)
      temp_array2(i,j) = obs_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Bbioratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_Bbioratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Bbioratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_Bbioratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Bbioratio0(i,j) = temp_corr_obs_Bbioratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Bbioratio0(i,j) = temp_corr_obs_Bbioratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Bbioratio0(i,j) = 9.96921e+36
        else
          corr_obs_Bbioratio0(i,j) = pattern_cor( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Bbioratio0(i,j) = 9.96921e+36
        else
          corr_obs_Bbioratio0(i,j) = pattern_cor( obs_Bbio(i-5:i+5,j-5:j+5),obs_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Bbioratio0 = mask(corr_obs_Bbioratio0, ismissing(obs_Bbio), False)
  corr_obs_Bbioratio0!0   = "lat"
  corr_obs_Bbioratio0!1   = "lon"
  corr_obs_Bbioratio0&lat = lat
  corr_obs_Bbioratio0&lon = lon
  corr_obs_Bbioratio0@comment = "Observation, ratio=Abio/Bbio, preindustrial"
  corr_obs_Bbioratio0@long_name = "Correlation of belowground biomass and Abio/Bbio"
  outfile->corr_obs_Bbioratio0 = corr_obs_Bbioratio0
  delete(corr_obs_Bbioratio0)
  delete(temp_corr_obs_Bbioratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Bbio vs ratio1 (cLeaf/Bbio)

  corr_obs_Bbioratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Bbioratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,j-10)
      temp_array2(i,j) = obs_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Bbioratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_Bbioratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Bbioratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_Bbioratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Bbioratio1(i,j) = temp_corr_obs_Bbioratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Bbioratio1(i,j) = temp_corr_obs_Bbioratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Bbioratio1(i,j) = 9.96921e+36
        else
          corr_obs_Bbioratio1(i,j) = pattern_cor( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Bbioratio1(i,j) = 9.96921e+36
        else
          corr_obs_Bbioratio1(i,j) = pattern_cor( obs_Bbio(i-5:i+5,j-5:j+5),obs_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Bbioratio1 = mask(corr_obs_Bbioratio1, ismissing(obs_Bbio), False)
  corr_obs_Bbioratio1!0   = "lat"
  corr_obs_Bbioratio1!1   = "lon"
  corr_obs_Bbioratio1&lat = lat
  corr_obs_Bbioratio1&lon = lon
  corr_obs_Bbioratio1@comment = "Observation, ratio=cLeaf/Bbio, preindustrial"
  corr_obs_Bbioratio1@long_name = "Correlation of belowground biomass and cLeaf/Bbio"
  outfile->corr_obs_Bbioratio1 = corr_obs_Bbioratio1
  delete(corr_obs_Bbioratio1)
  delete(temp_corr_obs_Bbioratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Bbio vs ratio2 (cLeaf/cWood)

  corr_obs_Bbioratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Bbioratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,j-10)
      temp_array2(i,j) = obs_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Bbioratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_Bbioratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Bbioratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_Bbioratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Bbioratio2(i,j) = temp_corr_obs_Bbioratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Bbioratio2(i,j) = temp_corr_obs_Bbioratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Bbioratio2(i,j) = 9.96921e+36
        else
          corr_obs_Bbioratio2(i,j) = pattern_cor( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Bbioratio2(i,j) = 9.96921e+36
        else
          corr_obs_Bbioratio2(i,j) = pattern_cor( obs_Bbio(i-5:i+5,j-5:j+5),obs_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Bbioratio2 = mask(corr_obs_Bbioratio2, ismissing(obs_Bbio), False)
  corr_obs_Bbioratio2!0   = "lat"
  corr_obs_Bbioratio2!1   = "lon"
  corr_obs_Bbioratio2&lat = lat
  corr_obs_Bbioratio2&lon = lon
  corr_obs_Bbioratio2@comment = "Observation, ratio=cLeaf/cWood, preindustrial"
  corr_obs_Bbioratio2@long_name = "Correlation of belowground biomass and cLeaf/cWood"
  outfile->corr_obs_Bbioratio2 = corr_obs_Bbioratio2
  delete(corr_obs_Bbioratio2)
  delete(temp_corr_obs_Bbioratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; Bbio vs ratio3 (Bbio/cWood)

  corr_obs_Bbioratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_Bbioratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_Bbio(i,j-10)
      temp_array2(i,j) = obs_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_Bbioratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_Bbioratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_Bbioratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_Bbioratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_Bbioratio3(i,j) = temp_corr_obs_Bbioratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_Bbioratio3(i,j) = temp_corr_obs_Bbioratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_Bbioratio3(i,j) = 9.96921e+36
        else
          corr_obs_Bbioratio3(i,j) = pattern_cor( obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_Bbio(i-5:i+5,j-5:j+5) ))) then
          corr_obs_Bbioratio3(i,j) = 9.96921e+36
        else
          corr_obs_Bbioratio3(i,j) = pattern_cor( obs_Bbio(i-5:i+5,j-5:j+5),obs_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_Bbioratio3 = mask(corr_obs_Bbioratio3, ismissing(obs_Bbio), False)
  corr_obs_Bbioratio3!0   = "lat"
  corr_obs_Bbioratio3!1   = "lon"
  corr_obs_Bbioratio3&lat = lat
  corr_obs_Bbioratio3&lon = lon
  corr_obs_Bbioratio3@comment = "Observation, ratio=Bbio/cWood, preindustrial"
  corr_obs_Bbioratio3@long_name = "Correlation of belowground biomass and Bbio/cWood"
  outfile->corr_obs_Bbioratio3 = corr_obs_Bbioratio3
  delete(corr_obs_Bbioratio3)
  delete(temp_corr_obs_Bbioratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; cLeaf vs cWood

  corr_obs_cLeafcWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_cLeafcWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,j-10)
      temp_array2(i,j) = obs_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_cLeafcWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafcWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_cLeafcWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafcWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_cLeafcWood(i,j) = temp_corr_obs_cLeafcWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_cLeafcWood(i,j) = temp_corr_obs_cLeafcWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_cLeafcWood(i,j) = 9.96921e+36
        else
          corr_obs_cLeafcWood(i,j) = pattern_cor( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:i+5,j-5:j+5) ))) then
          corr_obs_cLeafcWood(i,j) = 9.96921e+36
        else
          corr_obs_cLeafcWood(i,j) = pattern_cor( obs_cLeaf(i-5:i+5,j-5:j+5),obs_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_cLeafcWood = mask(corr_obs_cLeafcWood, ismissing(obs_cLeaf), False)
  corr_obs_cLeafcWood!0   = "lat"
  corr_obs_cLeafcWood!1   = "lon"
  corr_obs_cLeafcWood&lat = lat
  corr_obs_cLeafcWood&lon = lon
  corr_obs_cLeafcWood@comment = "Observation, preindustrial"
  corr_obs_cLeafcWood@long_name = "Correlation of leaf biomass and wood biomass"
  outfile->corr_obs_cLeafcWood = corr_obs_cLeafcWood
  delete(corr_obs_cLeafcWood)
  delete(temp_corr_obs_cLeafcWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; cLeaf vs ratio0 (Abio/Bbio)

  corr_obs_cLeafratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_cLeafratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,j-10)
      temp_array2(i,j) = obs_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_cLeafratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_cLeafratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_cLeafratio0(i,j) = temp_corr_obs_cLeafratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_cLeafratio0(i,j) = temp_corr_obs_cLeafratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_cLeafratio0(i,j) = 9.96921e+36
        else
          corr_obs_cLeafratio0(i,j) = pattern_cor( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:i+5,j-5:j+5) ))) then
          corr_obs_cLeafratio0(i,j) = 9.96921e+36
        else
          corr_obs_cLeafratio0(i,j) = pattern_cor( obs_cLeaf(i-5:i+5,j-5:j+5),obs_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_cLeafratio0 = mask(corr_obs_cLeafratio0, ismissing(obs_cLeaf), False)
  corr_obs_cLeafratio0!0   = "lat"
  corr_obs_cLeafratio0!1   = "lon"
  corr_obs_cLeafratio0&lat = lat
  corr_obs_cLeafratio0&lon = lon
  corr_obs_cLeafratio0@comment = "Observation, ratio=Abio/Bbio, preindustrial"
  corr_obs_cLeafratio0@long_name = "Correlation of leaf biomass and Abio/Bbio"
  outfile->corr_obs_cLeafratio0 = corr_obs_cLeafratio0
  delete(corr_obs_cLeafratio0)
  delete(temp_corr_obs_cLeafratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; cLeaf vs ratio1 (cLeaf/Bbio)

  corr_obs_cLeafratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_cLeafratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,j-10)
      temp_array2(i,j) = obs_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_cLeafratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_cLeafratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_cLeafratio1(i,j) = temp_corr_obs_cLeafratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_cLeafratio1(i,j) = temp_corr_obs_cLeafratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_cLeafratio1(i,j) = 9.96921e+36
        else
          corr_obs_cLeafratio1(i,j) = pattern_cor( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:i+5,j-5:j+5) ))) then
          corr_obs_cLeafratio1(i,j) = 9.96921e+36
        else
          corr_obs_cLeafratio1(i,j) = pattern_cor( obs_cLeaf(i-5:i+5,j-5:j+5),obs_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_cLeafratio1 = mask(corr_obs_cLeafratio1, ismissing(obs_cLeaf), False)
  corr_obs_cLeafratio1!0   = "lat"
  corr_obs_cLeafratio1!1   = "lon"
  corr_obs_cLeafratio1&lat = lat
  corr_obs_cLeafratio1&lon = lon
  corr_obs_cLeafratio1@comment = "Observation, ratio=cLeaf/Bbio, preindustrial"
  corr_obs_cLeafratio1@long_name = "Correlation of leaf biomass and cLeaf/Bbio"
  outfile->corr_obs_cLeafratio1 = corr_obs_cLeafratio1
  delete(corr_obs_cLeafratio1)
  delete(temp_corr_obs_cLeafratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; cLeaf vs ratio2 (cLeaf/cWood)

  corr_obs_cLeafratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_cLeafratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,j-10)
      temp_array2(i,j) = obs_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_cLeafratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_cLeafratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_cLeafratio2(i,j) = temp_corr_obs_cLeafratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_cLeafratio2(i,j) = temp_corr_obs_cLeafratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_cLeafratio2(i,j) = 9.96921e+36
        else
          corr_obs_cLeafratio2(i,j) = pattern_cor( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:i+5,j-5:j+5) ))) then
          corr_obs_cLeafratio2(i,j) = 9.96921e+36
        else
          corr_obs_cLeafratio2(i,j) = pattern_cor( obs_cLeaf(i-5:i+5,j-5:j+5),obs_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_cLeafratio2 = mask(corr_obs_cLeafratio2, ismissing(obs_cLeaf), False)
  corr_obs_cLeafratio2!0   = "lat"
  corr_obs_cLeafratio2!1   = "lon"
  corr_obs_cLeafratio2&lat = lat
  corr_obs_cLeafratio2&lon = lon
  corr_obs_cLeafratio2@comment = "Observation, ratio=cLeaf/cWood, preindustrial"
  corr_obs_cLeafratio2@long_name = "Correlation of leaf biomass and cLeaf/cWood"
  outfile->corr_obs_cLeafratio2 = corr_obs_cLeafratio2
  delete(corr_obs_cLeafratio2)
  delete(temp_corr_obs_cLeafratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; cLeaf vs ratio3 (Bbio/cWood)

  corr_obs_cLeafratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_cLeafratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_cLeaf(i,j-10)
      temp_array2(i,j) = obs_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_cLeafratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_cLeafratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_cLeafratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_cLeafratio3(i,j) = temp_corr_obs_cLeafratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_cLeafratio3(i,j) = temp_corr_obs_cLeafratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_cLeafratio3(i,j) = 9.96921e+36
        else
          corr_obs_cLeafratio3(i,j) = pattern_cor( obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cLeaf(i-5:i+5,j-5:j+5) ))) then
          corr_obs_cLeafratio3(i,j) = 9.96921e+36
        else
          corr_obs_cLeafratio3(i,j) = pattern_cor( obs_cLeaf(i-5:i+5,j-5:j+5),obs_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_cLeafratio3 = mask(corr_obs_cLeafratio3, ismissing(obs_cLeaf), False)
  corr_obs_cLeafratio3!0   = "lat"
  corr_obs_cLeafratio3!1   = "lon"
  corr_obs_cLeafratio3&lat = lat
  corr_obs_cLeafratio3&lon = lon
  corr_obs_cLeafratio3@comment = "Observation, ratio=Bbio/cWood, preindustrial"
  corr_obs_cLeafratio3@long_name = "Correlation of leaf biomass and Bbio/cWood"
  outfile->corr_obs_cLeafratio3 = corr_obs_cLeafratio3
  delete(corr_obs_cLeafratio3)
  delete(temp_corr_obs_cLeafratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; cWood vs ratio0 (Abio/Bbio)

  corr_obs_cWoodratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_cWoodratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_cWood(i,j-10)
      temp_array2(i,j) = obs_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_cWoodratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_cWoodratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_cWoodratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_cWoodratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_cWoodratio0(i,j) = temp_corr_obs_cWoodratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_cWoodratio0(i,j) = temp_corr_obs_cWoodratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_cWoodratio0(i,j) = 9.96921e+36
        else
          corr_obs_cWoodratio0(i,j) = pattern_cor( obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cWood(i-5:i+5,j-5:j+5) ))) then
          corr_obs_cWoodratio0(i,j) = 9.96921e+36
        else
          corr_obs_cWoodratio0(i,j) = pattern_cor( obs_cWood(i-5:i+5,j-5:j+5),obs_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_cWoodratio0 = mask(corr_obs_cWoodratio0, ismissing(obs_cWood), False)
  corr_obs_cWoodratio0!0   = "lat"
  corr_obs_cWoodratio0!1   = "lon"
  corr_obs_cWoodratio0&lat = lat
  corr_obs_cWoodratio0&lon = lon
  corr_obs_cWoodratio0@comment = "Observation, ratio=Abio/Bbio, preindustrial"
  corr_obs_cWoodratio0@long_name = "Correlation of wood biomass and Abio/Bbio"
  outfile->corr_obs_cWoodratio0 = corr_obs_cWoodratio0
  delete(corr_obs_cWoodratio0)
  delete(temp_corr_obs_cWoodratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; cWood vs ratio1 (cLeaf/Bbio)

  corr_obs_cWoodratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_cWoodratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_cWood(i,j-10)
      temp_array2(i,j) = obs_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_cWoodratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_cWoodratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_cWoodratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_cWoodratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_cWoodratio1(i,j) = temp_corr_obs_cWoodratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_cWoodratio1(i,j) = temp_corr_obs_cWoodratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_cWoodratio1(i,j) = 9.96921e+36
        else
          corr_obs_cWoodratio1(i,j) = pattern_cor( obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cWood(i-5:i+5,j-5:j+5) ))) then
          corr_obs_cWoodratio1(i,j) = 9.96921e+36
        else
          corr_obs_cWoodratio1(i,j) = pattern_cor( obs_cWood(i-5:i+5,j-5:j+5),obs_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_cWoodratio1 = mask(corr_obs_cWoodratio1, ismissing(obs_cWood), False)
  corr_obs_cWoodratio1!0   = "lat"
  corr_obs_cWoodratio1!1   = "lon"
  corr_obs_cWoodratio1&lat = lat
  corr_obs_cWoodratio1&lon = lon
  corr_obs_cWoodratio1@comment = "Observation, ratio=cLeaf/Bbio, preindustrial"
  corr_obs_cWoodratio1@long_name = "Correlation of wood biomass and cLeaf/Bbio"
  outfile->corr_obs_cWoodratio1 = corr_obs_cWoodratio1
  delete(corr_obs_cWoodratio1)
  delete(temp_corr_obs_cWoodratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; cWood vs ratio2 (cLeaf/cWood)

  corr_obs_cWoodratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_cWoodratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_cWood(i,j-10)
      temp_array2(i,j) = obs_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_cWoodratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_cWoodratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_cWoodratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_cWoodratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_cWoodratio2(i,j) = temp_corr_obs_cWoodratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_cWoodratio2(i,j) = temp_corr_obs_cWoodratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_cWoodratio2(i,j) = 9.96921e+36
        else
          corr_obs_cWoodratio2(i,j) = pattern_cor( obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cWood(i-5:i+5,j-5:j+5) ))) then
          corr_obs_cWoodratio2(i,j) = 9.96921e+36
        else
          corr_obs_cWoodratio2(i,j) = pattern_cor( obs_cWood(i-5:i+5,j-5:j+5),obs_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_cWoodratio2 = mask(corr_obs_cWoodratio2, ismissing(obs_cWood), False)
  corr_obs_cWoodratio2!0   = "lat"
  corr_obs_cWoodratio2!1   = "lon"
  corr_obs_cWoodratio2&lat = lat
  corr_obs_cWoodratio2&lon = lon
  corr_obs_cWoodratio2@comment = "Observation, ratio=cLeaf/cWood, preindustrial"
  corr_obs_cWoodratio2@long_name = "Correlation of wood biomass and cLeaf/cWood"
  outfile->corr_obs_cWoodratio2 = corr_obs_cWoodratio2
  delete(corr_obs_cWoodratio2)
  delete(temp_corr_obs_cWoodratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; cWood vs ratio3 (Bbio/cWood)

  corr_obs_cWoodratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_cWoodratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_cWood(i,j-10)
      temp_array2(i,j) = obs_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_cWoodratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_cWoodratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_cWoodratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_cWoodratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_cWoodratio3(i,j) = temp_corr_obs_cWoodratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_cWoodratio3(i,j) = temp_corr_obs_cWoodratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_cWoodratio3(i,j) = 9.96921e+36
        else
          corr_obs_cWoodratio3(i,j) = pattern_cor( obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_cWood(i-5:i+5,j-5:j+5) ))) then
          corr_obs_cWoodratio3(i,j) = 9.96921e+36
        else
          corr_obs_cWoodratio3(i,j) = pattern_cor( obs_cWood(i-5:i+5,j-5:j+5),obs_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_cWoodratio3 = mask(corr_obs_cWoodratio3, ismissing(obs_cWood), False)
  corr_obs_cWoodratio3!0   = "lat"
  corr_obs_cWoodratio3!1   = "lon"
  corr_obs_cWoodratio3&lat = lat
  corr_obs_cWoodratio3&lon = lon
  corr_obs_cWoodratio3@comment = "Observation, ratio=Bbio/cWood, preindustrial"
  corr_obs_cWoodratio3@long_name = "Correlation of wood biomass and Bbio/cWood"
  outfile->corr_obs_cWoodratio3 = corr_obs_cWoodratio3
  delete(corr_obs_cWoodratio3)
  delete(temp_corr_obs_cWoodratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------

;;; Delete variables

  delete(temp_obs_Bbio)
  delete(obs_Bbio)
  delete(obs_Tbio)
  delete(obs_cLeaf)
  delete(temp_obs_cWood)
  delete(obs_cWood)
  delete(obs_Abio)
  delete(obs_ratio0)
  delete(obs_ratio1)
  delete(obs_ratio2)
  delete(obs_ratio3)

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

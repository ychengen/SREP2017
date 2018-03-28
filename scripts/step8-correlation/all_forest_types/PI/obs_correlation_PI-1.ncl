;==========================================================
; Compute statistics of the BGI dataset
; For preindustrial time (1861-1885)
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

;;; Meteorology
      obs_pr = tofloat(infile->mask_obs_pr)
      obs_rsds = tofloat(infile->mask_obs_rsds)
      obs_tas = tofloat(infile->mask_obs_tas)

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
      fout    = outdir+"correlation-"+obsnm(z)+"-OBS-1.nc"
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

;;;;; PR vs TAS

  corr_obs_prtas = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prtas = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_tas(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_tas(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prtas(i,j) = 9.96921e+36
        else
          temp_corr_obs_prtas(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prtas(i,j) = 9.96921e+36
        else
          temp_corr_obs_prtas(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prtas(i,j) = temp_corr_obs_prtas(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prtas(i,j) = temp_corr_obs_prtas(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prtas(i,j) = 9.96921e+36
        else
          corr_obs_prtas(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_tas(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prtas(i,j) = 9.96921e+36
        else
          corr_obs_prtas(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_tas(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prtas = mask(corr_obs_prtas, ismissing(obs_pr), False)
  corr_obs_prtas!0   = "lat"
  corr_obs_prtas!1   = "lon"
  corr_obs_prtas&lat = lat
  corr_obs_prtas&lon = lon
  corr_obs_prtas@comment = "Observation, preindustrial"
  corr_obs_prtas@long_name = "Correlation of precipitation and surface temperature"
  outfile->corr_obs_prtas = corr_obs_prtas
  delete(corr_obs_prtas)
  delete(temp_corr_obs_prtas)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs RSDS

  corr_obs_prrsds = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prrsds = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_rsds(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prrsds(i,j) = 9.96921e+36
        else
          temp_corr_obs_prrsds(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prrsds(i,j) = 9.96921e+36
        else
          temp_corr_obs_prrsds(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prrsds(i,j) = temp_corr_obs_prrsds(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prrsds(i,j) = temp_corr_obs_prrsds(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prrsds(i,j) = 9.96921e+36
        else
          corr_obs_prrsds(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prrsds(i,j) = 9.96921e+36
        else
          corr_obs_prrsds(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_rsds(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prrsds = mask(corr_obs_prrsds, ismissing(obs_pr), False)
  corr_obs_prrsds!0   = "lat"
  corr_obs_prrsds!1   = "lon"
  corr_obs_prrsds&lat = lat
  corr_obs_prrsds&lon = lon
  corr_obs_prrsds@comment = "Observation, preindustrial"
  corr_obs_prrsds@long_name = "Correlation of precipitation and shortwave radiation"
  outfile->corr_obs_prrsds = corr_obs_prrsds
  delete(corr_obs_prrsds)
  delete(temp_corr_obs_prrsds)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs Abio

  corr_obs_prAbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prAbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_Abio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prAbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_prAbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prAbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_prAbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prAbio(i,j) = temp_corr_obs_prAbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prAbio(i,j) = temp_corr_obs_prAbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prAbio(i,j) = 9.96921e+36
        else
          corr_obs_prAbio(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prAbio(i,j) = 9.96921e+36
        else
          corr_obs_prAbio(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_Abio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prAbio = mask(corr_obs_prAbio, ismissing(obs_pr), False)
  corr_obs_prAbio!0   = "lat"
  corr_obs_prAbio!1   = "lon"
  corr_obs_prAbio&lat = lat
  corr_obs_prAbio&lon = lon
  corr_obs_prAbio@comment = "Observation, preindustrial"
  corr_obs_prAbio@long_name = "Correlation of precipitation and aboveground biomass"
  outfile->corr_obs_prAbio = corr_obs_prAbio
  delete(corr_obs_prAbio)
  delete(temp_corr_obs_prAbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs Bbio

  corr_obs_prBbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prBbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_Bbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_prBbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_prBbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prBbio(i,j) = temp_corr_obs_prBbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prBbio(i,j) = temp_corr_obs_prBbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prBbio(i,j) = 9.96921e+36
        else
          corr_obs_prBbio(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prBbio(i,j) = 9.96921e+36
        else
          corr_obs_prBbio(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_Bbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prBbio = mask(corr_obs_prBbio, ismissing(obs_pr), False)
  corr_obs_prBbio!0   = "lat"
  corr_obs_prBbio!1   = "lon"
  corr_obs_prBbio&lat = lat
  corr_obs_prBbio&lon = lon
  corr_obs_prBbio@comment = "Observation, preindustrial"
  corr_obs_prBbio@long_name = "Correlation of precipitation and belowground biomass"
  outfile->corr_obs_prBbio = corr_obs_prBbio
  delete(corr_obs_prBbio)
  delete(temp_corr_obs_prBbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs Tbio

  corr_obs_prTbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prTbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_Tbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prTbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_prTbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prTbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_prTbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prTbio(i,j) = temp_corr_obs_prTbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prTbio(i,j) = temp_corr_obs_prTbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prTbio(i,j) = 9.96921e+36
        else
          corr_obs_prTbio(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prTbio(i,j) = 9.96921e+36
        else
          corr_obs_prTbio(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_Tbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prTbio = mask(corr_obs_prTbio, ismissing(obs_pr), False)
  corr_obs_prTbio!0   = "lat"
  corr_obs_prTbio!1   = "lon"
  corr_obs_prTbio&lat = lat
  corr_obs_prTbio&lon = lon
  corr_obs_prTbio@comment = "Observation, preindustrial"
  corr_obs_prTbio@long_name = "Correlation of precipitation and total biomass"
  outfile->corr_obs_prTbio = corr_obs_prTbio
  delete(corr_obs_prTbio)
  delete(temp_corr_obs_prTbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs cLeaf

  corr_obs_prcLeaf = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prcLeaf = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_cLeaf(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prcLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_prcLeaf(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prcLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_prcLeaf(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prcLeaf(i,j) = temp_corr_obs_prcLeaf(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prcLeaf(i,j) = temp_corr_obs_prcLeaf(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prcLeaf(i,j) = 9.96921e+36
        else
          corr_obs_prcLeaf(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prcLeaf(i,j) = 9.96921e+36
        else
          corr_obs_prcLeaf(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_cLeaf(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prcLeaf = mask(corr_obs_prcLeaf, ismissing(obs_pr), False)
  corr_obs_prcLeaf!0   = "lat"
  corr_obs_prcLeaf!1   = "lon"
  corr_obs_prcLeaf&lat = lat
  corr_obs_prcLeaf&lon = lon
  corr_obs_prcLeaf@comment = "Observation, preindustrial"
  corr_obs_prcLeaf@long_name = "Correlation of precipitation and leaf biomass"
  outfile->corr_obs_prcLeaf = corr_obs_prcLeaf
  delete(corr_obs_prcLeaf)
  delete(temp_corr_obs_prcLeaf)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; PR vs cWood

  corr_obs_prcWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prcWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prcWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_prcWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prcWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_prcWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prcWood(i,j) = temp_corr_obs_prcWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prcWood(i,j) = temp_corr_obs_prcWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prcWood(i,j) = 9.96921e+36
        else
          corr_obs_prcWood(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prcWood(i,j) = 9.96921e+36
        else
          corr_obs_prcWood(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prcWood = mask(corr_obs_prcWood, ismissing(obs_pr), False)
  corr_obs_prcWood!0   = "lat"
  corr_obs_prcWood!1   = "lon"
  corr_obs_prcWood&lat = lat
  corr_obs_prcWood&lon = lon
  corr_obs_prcWood@comment = "Observation, preindustrial"
  corr_obs_prcWood@long_name = "Correlation of precipitation and wood biomass"
  outfile->corr_obs_prcWood = corr_obs_prcWood
  delete(corr_obs_prcWood)
  delete(temp_corr_obs_prcWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; PR vs ratio0 (Abio/Bbio)

  corr_obs_prratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_prratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_prratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prratio0(i,j) = temp_corr_obs_prratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prratio0(i,j) = temp_corr_obs_prratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prratio0(i,j) = 9.96921e+36
        else
          corr_obs_prratio0(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prratio0(i,j) = 9.96921e+36
        else
          corr_obs_prratio0(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prratio0 = mask(corr_obs_prratio0, ismissing(obs_pr), False)
  corr_obs_prratio0!0   = "lat"
  corr_obs_prratio0!1   = "lon"
  corr_obs_prratio0&lat = lat
  corr_obs_prratio0&lon = lon
  corr_obs_prratio0@comment = "Observation, ratio=Abio/Bbio, preindustrial"
  corr_obs_prratio0@long_name = "Correlation of precipitation and Abio/Bbio"
  outfile->corr_obs_prratio0 = corr_obs_prratio0
  delete(corr_obs_prratio0)
  delete(temp_corr_obs_prratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; PR vs ratio1 (cLeaf/Bbio)

  corr_obs_prratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_prratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_prratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prratio1(i,j) = temp_corr_obs_prratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prratio1(i,j) = temp_corr_obs_prratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prratio1(i,j) = 9.96921e+36
        else
          corr_obs_prratio1(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prratio1(i,j) = 9.96921e+36
        else
          corr_obs_prratio1(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prratio1 = mask(corr_obs_prratio1, ismissing(obs_pr), False)
  corr_obs_prratio1!0   = "lat"
  corr_obs_prratio1!1   = "lon"
  corr_obs_prratio1&lat = lat
  corr_obs_prratio1&lon = lon
  corr_obs_prratio1@comment = "Observation, ratio=cLeaf/Bbio, preindustrial"
  corr_obs_prratio1@long_name = "Correlation of precipitation and cLeaf/Bbio"
  outfile->corr_obs_prratio1 = corr_obs_prratio1
  delete(corr_obs_prratio1)
  delete(temp_corr_obs_prratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; PR vs ratio2 (cLeaf/cWood)

  corr_obs_prratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_prratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_prratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prratio2(i,j) = temp_corr_obs_prratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prratio2(i,j) = temp_corr_obs_prratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prratio2(i,j) = 9.96921e+36
        else
          corr_obs_prratio2(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prratio2(i,j) = 9.96921e+36
        else
          corr_obs_prratio2(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prratio2 = mask(corr_obs_prratio2, ismissing(obs_pr), False)
  corr_obs_prratio2!0   = "lat"
  corr_obs_prratio2!1   = "lon"
  corr_obs_prratio2&lat = lat
  corr_obs_prratio2&lon = lon
  corr_obs_prratio2@comment = "Observation, ratio=cLeaf/cWood, preindustrial"
  corr_obs_prratio2@long_name = "Correlation of precipitation and cLeaf/cWood"
  outfile->corr_obs_prratio2 = corr_obs_prratio2
  delete(corr_obs_prratio2)
  delete(temp_corr_obs_prratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; PR vs ratio3 (Bbio/cWood)

  corr_obs_prratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_prratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_pr(i,j-10)
      temp_array2(i,j) = obs_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_prratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_prratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_prratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_prratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_prratio3(i,j) = temp_corr_obs_prratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_prratio3(i,j) = temp_corr_obs_prratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_prratio3(i,j) = 9.96921e+36
        else
          corr_obs_prratio3(i,j) = pattern_cor( obs_pr(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_pr(i-5:i+5,j-5:j+5) ))) then
          corr_obs_prratio3(i,j) = 9.96921e+36
        else
          corr_obs_prratio3(i,j) = pattern_cor( obs_pr(i-5:i+5,j-5:j+5),obs_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_prratio3 = mask(corr_obs_prratio3, ismissing(obs_pr), False)
  corr_obs_prratio3!0   = "lat"
  corr_obs_prratio3!1   = "lon"
  corr_obs_prratio3&lat = lat
  corr_obs_prratio3&lon = lon
  corr_obs_prratio3@comment = "Observation, ratio=Bbio/cWood, preindustrial"
  corr_obs_prratio3@long_name = "Correlation of precipitation and Bbio/cWood"
  outfile->corr_obs_prratio3 = corr_obs_prratio3
  delete(corr_obs_prratio3)
  delete(temp_corr_obs_prratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;;;; TAS vs RSDS

  corr_obs_tasrsds = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tasrsds = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_rsds(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tasrsds(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasrsds(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tasrsds(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasrsds(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tasrsds(i,j) = temp_corr_obs_tasrsds(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tasrsds(i,j) = temp_corr_obs_tasrsds(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tasrsds(i,j) = 9.96921e+36
        else
          corr_obs_tasrsds(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tasrsds(i,j) = 9.96921e+36
        else
          corr_obs_tasrsds(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_rsds(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tasrsds = mask(corr_obs_tasrsds, ismissing(obs_tas), False)
  corr_obs_tasrsds!0   = "lat"
  corr_obs_tasrsds!1   = "lon"
  corr_obs_tasrsds&lat = lat
  corr_obs_tasrsds&lon = lon
  corr_obs_tasrsds@comment = "Observation, preindustrial"
  corr_obs_tasrsds@long_name = "Correlation of surface temperature and shortwave radiation"
  outfile->corr_obs_tasrsds = corr_obs_tasrsds
  delete(corr_obs_tasrsds)
  delete(temp_corr_obs_tasrsds)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs Abio

  corr_obs_tasAbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tasAbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_Abio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tasAbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasAbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tasAbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasAbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tasAbio(i,j) = temp_corr_obs_tasAbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tasAbio(i,j) = temp_corr_obs_tasAbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tasAbio(i,j) = 9.96921e+36
        else
          corr_obs_tasAbio(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tasAbio(i,j) = 9.96921e+36
        else
          corr_obs_tasAbio(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_Abio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tasAbio = mask(corr_obs_tasAbio, ismissing(obs_tas), False)
  corr_obs_tasAbio!0   = "lat"
  corr_obs_tasAbio!1   = "lon"
  corr_obs_tasAbio&lat = lat
  corr_obs_tasAbio&lon = lon
  corr_obs_tasAbio@comment = "Observation, Tbio-Bbio, preindustrial"
  corr_obs_tasAbio@long_name = "Correlation of surface temperature and aboveground biomass"
  outfile->corr_obs_tasAbio = corr_obs_tasAbio
  delete(corr_obs_tasAbio)
  delete(temp_corr_obs_tasAbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs Bbio

  corr_obs_tasBbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tasBbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_Bbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tasBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasBbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tasBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasBbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tasBbio(i,j) = temp_corr_obs_tasBbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tasBbio(i,j) = temp_corr_obs_tasBbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tasBbio(i,j) = 9.96921e+36
        else
          corr_obs_tasBbio(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tasBbio(i,j) = 9.96921e+36
        else
          corr_obs_tasBbio(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_Bbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tasBbio = mask(corr_obs_tasBbio, ismissing(obs_tas), False)
  corr_obs_tasBbio!0   = "lat"
  corr_obs_tasBbio!1   = "lon"
  corr_obs_tasBbio&lat = lat
  corr_obs_tasBbio&lon = lon
  corr_obs_tasBbio@comment = "Observation, preindustrial"
  corr_obs_tasBbio@long_name = "Correlation of surface temperature and belowground biomass"
  outfile->corr_obs_tasBbio = corr_obs_tasBbio
  delete(corr_obs_tasBbio)
  delete(temp_corr_obs_tasBbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs Tbio

  corr_obs_tasTbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tasTbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_Tbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tasTbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasTbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tasTbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasTbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tasTbio(i,j) = temp_corr_obs_tasTbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tasTbio(i,j) = temp_corr_obs_tasTbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tasTbio(i,j) = 9.96921e+36
        else
          corr_obs_tasTbio(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tasTbio(i,j) = 9.96921e+36
        else
          corr_obs_tasTbio(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_Tbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tasTbio = mask(corr_obs_tasTbio, ismissing(obs_tas), False)
  corr_obs_tasTbio!0   = "lat"
  corr_obs_tasTbio!1   = "lon"
  corr_obs_tasTbio&lat = lat
  corr_obs_tasTbio&lon = lon
  corr_obs_tasTbio@comment = "Observation, preindustrial"
  corr_obs_tasTbio@long_name = "Correlation of surface temperature and total biomass"
  outfile->corr_obs_tasTbio = corr_obs_tasTbio
  delete(corr_obs_tasTbio)
  delete(temp_corr_obs_tasTbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs cLeaf

  corr_obs_tascLeaf = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tascLeaf = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_cLeaf(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tascLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_tascLeaf(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tascLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_tascLeaf(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tascLeaf(i,j) = temp_corr_obs_tascLeaf(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tascLeaf(i,j) = temp_corr_obs_tascLeaf(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tascLeaf(i,j) = 9.96921e+36
        else
          corr_obs_tascLeaf(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tascLeaf(i,j) = 9.96921e+36
        else
          corr_obs_tascLeaf(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_cLeaf(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tascLeaf = mask(corr_obs_tascLeaf, ismissing(obs_tas), False)
  corr_obs_tascLeaf!0   = "lat"
  corr_obs_tascLeaf!1   = "lon"
  corr_obs_tascLeaf&lat = lat
  corr_obs_tascLeaf&lon = lon
  corr_obs_tascLeaf@comment = "Observation, preindustrial"
  corr_obs_tascLeaf@long_name = "Correlation of surface temperature and leaf biomass"
  outfile->corr_obs_tascLeaf = corr_obs_tascLeaf
  delete(corr_obs_tascLeaf)
  delete(temp_corr_obs_tascLeaf)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs cWood

  corr_obs_tascWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tascWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tascWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_tascWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tascWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_tascWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tascWood(i,j) = temp_corr_obs_tascWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tascWood(i,j) = temp_corr_obs_tascWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tascWood(i,j) = 9.96921e+36
        else
          corr_obs_tascWood(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tascWood(i,j) = 9.96921e+36
        else
          corr_obs_tascWood(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tascWood = mask(corr_obs_tascWood, ismissing(obs_tas), False)
  corr_obs_tascWood!0   = "lat"
  corr_obs_tascWood!1   = "lon"
  corr_obs_tascWood&lat = lat
  corr_obs_tascWood&lon = lon
  corr_obs_tascWood@comment = "Observation, preindustrial"
  corr_obs_tascWood@long_name = "Correlation of surface temperature and wood biomass"
  outfile->corr_obs_tascWood = corr_obs_tascWood
  delete(corr_obs_tascWood)
  delete(temp_corr_obs_tascWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs ratio0 (Abio/Bbio)

  corr_obs_tasratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tasratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tasratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tasratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tasratio0(i,j) = temp_corr_obs_tasratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tasratio0(i,j) = temp_corr_obs_tasratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tasratio0(i,j) = 9.96921e+36
        else
          corr_obs_tasratio0(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tasratio0(i,j) = 9.96921e+36
        else
          corr_obs_tasratio0(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tasratio0 = mask(corr_obs_tasratio0, ismissing(obs_tas), False)
  corr_obs_tasratio0!0   = "lat"
  corr_obs_tasratio0!1   = "lon"
  corr_obs_tasratio0&lat = lat
  corr_obs_tasratio0&lon = lon
  corr_obs_tasratio0@comment = "Observation, ratio=Abio/Bbio, preindustrial"
  corr_obs_tasratio0@long_name = "Correlation of surface temperature and Abio/Bbio"
  outfile->corr_obs_tasratio0 = corr_obs_tasratio0
  delete(corr_obs_tasratio0)
  delete(temp_corr_obs_tasratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs ratio1 (cLeaf/Bbio)

  corr_obs_tasratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tasratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tasratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tasratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tasratio1(i,j) = temp_corr_obs_tasratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tasratio1(i,j) = temp_corr_obs_tasratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tasratio1(i,j) = 9.96921e+36
        else
          corr_obs_tasratio1(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tasratio1(i,j) = 9.96921e+36
        else
          corr_obs_tasratio1(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tasratio1 = mask(corr_obs_tasratio1, ismissing(obs_tas), False)
  corr_obs_tasratio1!0   = "lat"
  corr_obs_tasratio1!1   = "lon"
  corr_obs_tasratio1&lat = lat
  corr_obs_tasratio1&lon = lon
  corr_obs_tasratio1@comment = "Observation, ratio=cLeaf/Bbio, preindustrial"
  corr_obs_tasratio1@long_name = "Correlation of surface temperature and cLeaf/Bbio"
  outfile->corr_obs_tasratio1 = corr_obs_tasratio1
  delete(corr_obs_tasratio1)
  delete(temp_corr_obs_tasratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs ratio2 (cLeaf/cWood)

  corr_obs_tasratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tasratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tasratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tasratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tasratio2(i,j) = temp_corr_obs_tasratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tasratio2(i,j) = temp_corr_obs_tasratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tasratio2(i,j) = 9.96921e+36
        else
          corr_obs_tasratio2(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tasratio2(i,j) = 9.96921e+36
        else
          corr_obs_tasratio2(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tasratio2 = mask(corr_obs_tasratio2, ismissing(obs_tas), False)
  corr_obs_tasratio2!0   = "lat"
  corr_obs_tasratio2!1   = "lon"
  corr_obs_tasratio2&lat = lat
  corr_obs_tasratio2&lon = lon
  corr_obs_tasratio2@comment = "Observation, ratio=cLeaf/cWood, preindustrial"
  corr_obs_tasratio2@long_name = "Correlation of surface temperature and cLeaf/cWood"
  outfile->corr_obs_tasratio2 = corr_obs_tasratio2
  delete(corr_obs_tasratio2)
  delete(temp_corr_obs_tasratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; TAS vs ratio3 (Bbio/cWood)

  corr_obs_tasratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_tasratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_tas(i,j-10)
      temp_array2(i,j) = obs_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_tasratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_tasratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_tasratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_tasratio3(i,j) = temp_corr_obs_tasratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_tasratio3(i,j) = temp_corr_obs_tasratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_tasratio3(i,j) = 9.96921e+36
        else
          corr_obs_tasratio3(i,j) = pattern_cor( obs_tas(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_tas(i-5:i+5,j-5:j+5) ))) then
          corr_obs_tasratio3(i,j) = 9.96921e+36
        else
          corr_obs_tasratio3(i,j) = pattern_cor( obs_tas(i-5:i+5,j-5:j+5),obs_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_tasratio3 = mask(corr_obs_tasratio3, ismissing(obs_tas), False)
  corr_obs_tasratio3!0   = "lat"
  corr_obs_tasratio3!1   = "lon"
  corr_obs_tasratio3&lat = lat
  corr_obs_tasratio3&lon = lon
  corr_obs_tasratio3@comment = "Observation, ratio=Bbio/cWood, preindustrial"
  corr_obs_tasratio3@long_name = "Correlation of surface temperature and Bbio/cWood"
  outfile->corr_obs_tasratio3 = corr_obs_tasratio3
  delete(corr_obs_tasratio3)
  delete(temp_corr_obs_tasratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs Abio

  corr_obs_rsdsAbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_rsdsAbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Abio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,j-10)
      temp_array2(i,j) = obs_Abio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_rsdsAbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsAbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_rsdsAbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsAbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_rsdsAbio(i,j) = temp_corr_obs_rsdsAbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_rsdsAbio(i,j) = temp_corr_obs_rsdsAbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_rsdsAbio(i,j) = 9.96921e+36
        else
          corr_obs_rsdsAbio(i,j) = pattern_cor( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Abio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_obs_rsdsAbio(i,j) = 9.96921e+36
        else
          corr_obs_rsdsAbio(i,j) = pattern_cor( obs_rsds(i-5:i+5,j-5:j+5),obs_Abio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_rsdsAbio = mask(corr_obs_rsdsAbio, ismissing(obs_rsds), False)
  corr_obs_rsdsAbio!0   = "lat"
  corr_obs_rsdsAbio!1   = "lon"
  corr_obs_rsdsAbio&lat = lat
  corr_obs_rsdsAbio&lon = lon
  corr_obs_rsdsAbio@comment = "Observation, Tbio-Bbio, preindustrial"
  corr_obs_rsdsAbio@long_name = "Correlation of shortwave radiation and aboveground biomass"
  outfile->corr_obs_rsdsAbio = corr_obs_rsdsAbio
  delete(corr_obs_rsdsAbio)
  delete(temp_corr_obs_rsdsAbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs Bbio

  corr_obs_rsdsBbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_rsdsBbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Bbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,j-10)
      temp_array2(i,j) = obs_Bbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_rsdsBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsBbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_rsdsBbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsBbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_rsdsBbio(i,j) = temp_corr_obs_rsdsBbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_rsdsBbio(i,j) = temp_corr_obs_rsdsBbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_rsdsBbio(i,j) = 9.96921e+36
        else
          corr_obs_rsdsBbio(i,j) = pattern_cor( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Bbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_obs_rsdsBbio(i,j) = 9.96921e+36
        else
          corr_obs_rsdsBbio(i,j) = pattern_cor( obs_rsds(i-5:i+5,j-5:j+5),obs_Bbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_rsdsBbio = mask(corr_obs_rsdsBbio, ismissing(obs_rsds), False)
  corr_obs_rsdsBbio!0   = "lat"
  corr_obs_rsdsBbio!1   = "lon"
  corr_obs_rsdsBbio&lat = lat
  corr_obs_rsdsBbio&lon = lon
  corr_obs_rsdsBbio@comment = "Observation, preindustrial"
  corr_obs_rsdsBbio@long_name = "Correlation of shortwave radiation and belowground biomass"
  outfile->corr_obs_rsdsBbio = corr_obs_rsdsBbio
  delete(corr_obs_rsdsBbio)
  delete(temp_corr_obs_rsdsBbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs Tbio

  corr_obs_rsdsTbio = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_rsdsTbio = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_Tbio(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,j-10)
      temp_array2(i,j) = obs_Tbio(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_rsdsTbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsTbio(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_rsdsTbio(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsTbio(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_rsdsTbio(i,j) = temp_corr_obs_rsdsTbio(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_rsdsTbio(i,j) = temp_corr_obs_rsdsTbio(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_rsdsTbio(i,j) = 9.96921e+36
        else
          corr_obs_rsdsTbio(i,j) = pattern_cor( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_Tbio(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_obs_rsdsTbio(i,j) = 9.96921e+36
        else
          corr_obs_rsdsTbio(i,j) = pattern_cor( obs_rsds(i-5:i+5,j-5:j+5),obs_Tbio(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_rsdsTbio = mask(corr_obs_rsdsTbio, ismissing(obs_rsds), False)
  corr_obs_rsdsTbio!0   = "lat"
  corr_obs_rsdsTbio!1   = "lon"
  corr_obs_rsdsTbio&lat = lat
  corr_obs_rsdsTbio&lon = lon
  corr_obs_rsdsTbio@comment = "Observation, preindustrial"
  corr_obs_rsdsTbio@long_name = "Correlation of shortwave radiation and total biomass"
  outfile->corr_obs_rsdsTbio = corr_obs_rsdsTbio
  delete(corr_obs_rsdsTbio)
  delete(temp_corr_obs_rsdsTbio)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs cLeaf

  corr_obs_rsdscLeaf = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_rsdscLeaf = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cLeaf(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,j-10)
      temp_array2(i,j) = obs_cLeaf(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_rsdscLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdscLeaf(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_rsdscLeaf(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdscLeaf(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_rsdscLeaf(i,j) = temp_corr_obs_rsdscLeaf(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_rsdscLeaf(i,j) = temp_corr_obs_rsdscLeaf(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_rsdscLeaf(i,j) = 9.96921e+36
        else
          corr_obs_rsdscLeaf(i,j) = pattern_cor( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cLeaf(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_obs_rsdscLeaf(i,j) = 9.96921e+36
        else
          corr_obs_rsdscLeaf(i,j) = pattern_cor( obs_rsds(i-5:i+5,j-5:j+5),obs_cLeaf(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_rsdscLeaf = mask(corr_obs_rsdscLeaf, ismissing(obs_rsds), False)
  corr_obs_rsdscLeaf!0   = "lat"
  corr_obs_rsdscLeaf!1   = "lon"
  corr_obs_rsdscLeaf&lat = lat
  corr_obs_rsdscLeaf&lon = lon
  corr_obs_rsdscLeaf@comment = "Observation, preindustrial"
  corr_obs_rsdscLeaf@long_name = "Correlation of shortwave radiation and leaf biomass"
  outfile->corr_obs_rsdscLeaf = corr_obs_rsdscLeaf
  delete(corr_obs_rsdscLeaf)
  delete(temp_corr_obs_rsdscLeaf)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs cWood

  corr_obs_rsdscWood = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_rsdscWood = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_cWood(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,j-10)
      temp_array2(i,j) = obs_cWood(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_rsdscWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdscWood(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_rsdscWood(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdscWood(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_rsdscWood(i,j) = temp_corr_obs_rsdscWood(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_rsdscWood(i,j) = temp_corr_obs_rsdscWood(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_rsdscWood(i,j) = 9.96921e+36
        else
          corr_obs_rsdscWood(i,j) = pattern_cor( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_cWood(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_obs_rsdscWood(i,j) = 9.96921e+36
        else
          corr_obs_rsdscWood(i,j) = pattern_cor( obs_rsds(i-5:i+5,j-5:j+5),obs_cWood(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_rsdscWood = mask(corr_obs_rsdscWood, ismissing(obs_rsds), False)
  corr_obs_rsdscWood!0   = "lat"
  corr_obs_rsdscWood!1   = "lon"
  corr_obs_rsdscWood&lat = lat
  corr_obs_rsdscWood&lon = lon
  corr_obs_rsdscWood@comment = "Observation, preindustrial"
  corr_obs_rsdscWood@long_name = "Correlation of shortwave radiation and wood biomass"
  outfile->corr_obs_rsdscWood = corr_obs_rsdscWood
  delete(corr_obs_rsdscWood)
  delete(temp_corr_obs_rsdscWood)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs ratio0 (Abio/Bbio)

  corr_obs_rsdsratio0 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_rsdsratio0 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio0(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,j-10)
      temp_array2(i,j) = obs_ratio0(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_rsdsratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsratio0(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_rsdsratio0(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsratio0(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_rsdsratio0(i,j) = temp_corr_obs_rsdsratio0(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_rsdsratio0(i,j) = temp_corr_obs_rsdsratio0(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_rsdsratio0(i,j) = 9.96921e+36
        else
          corr_obs_rsdsratio0(i,j) = pattern_cor( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio0(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_obs_rsdsratio0(i,j) = 9.96921e+36
        else
          corr_obs_rsdsratio0(i,j) = pattern_cor( obs_rsds(i-5:i+5,j-5:j+5),obs_ratio0(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_rsdsratio0 = mask(corr_obs_rsdsratio0, ismissing(obs_rsds), False)
  corr_obs_rsdsratio0!0   = "lat"
  corr_obs_rsdsratio0!1   = "lon"
  corr_obs_rsdsratio0&lat = lat
  corr_obs_rsdsratio0&lon = lon
  corr_obs_rsdsratio0@comment = "Observation, ratio=Abio/Bbio, preindustrial"
  corr_obs_rsdsratio0@long_name = "Correlation of shortwave radiation and Abio/Bbio"
  outfile->corr_obs_rsdsratio0 = corr_obs_rsdsratio0
  delete(corr_obs_rsdsratio0)
  delete(temp_corr_obs_rsdsratio0)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs ratio1 (cLeaf/Bbio)

  corr_obs_rsdsratio1 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_rsdsratio1 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio1(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,j-10)
      temp_array2(i,j) = obs_ratio1(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_rsdsratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsratio1(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_rsdsratio1(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsratio1(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_rsdsratio1(i,j) = temp_corr_obs_rsdsratio1(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_rsdsratio1(i,j) = temp_corr_obs_rsdsratio1(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_rsdsratio1(i,j) = 9.96921e+36
        else
          corr_obs_rsdsratio1(i,j) = pattern_cor( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio1(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_obs_rsdsratio1(i,j) = 9.96921e+36
        else
          corr_obs_rsdsratio1(i,j) = pattern_cor( obs_rsds(i-5:i+5,j-5:j+5),obs_ratio1(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_rsdsratio1 = mask(corr_obs_rsdsratio1, ismissing(obs_rsds), False)
  corr_obs_rsdsratio1!0   = "lat"
  corr_obs_rsdsratio1!1   = "lon"
  corr_obs_rsdsratio1&lat = lat
  corr_obs_rsdsratio1&lon = lon
  corr_obs_rsdsratio1@comment = "Observation, ratio=cLeaf/Bbio, preindustrial"
  corr_obs_rsdsratio1@long_name = "Correlation of shortwave radiation and cLeaf/Bbio"
  outfile->corr_obs_rsdsratio1 = corr_obs_rsdsratio1
  delete(corr_obs_rsdsratio1)
  delete(temp_corr_obs_rsdsratio1)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs ratio2 (cLeaf/cWood)

  corr_obs_rsdsratio2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_rsdsratio2 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio2(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,j-10)
      temp_array2(i,j) = obs_ratio2(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_rsdsratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsratio2(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_rsdsratio2(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsratio2(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_rsdsratio2(i,j) = temp_corr_obs_rsdsratio2(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_rsdsratio2(i,j) = temp_corr_obs_rsdsratio2(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_rsdsratio2(i,j) = 9.96921e+36
        else
          corr_obs_rsdsratio2(i,j) = pattern_cor( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_obs_rsdsratio2(i,j) = 9.96921e+36
        else
          corr_obs_rsdsratio2(i,j) = pattern_cor( obs_rsds(i-5:i+5,j-5:j+5),obs_ratio2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_rsdsratio2 = mask(corr_obs_rsdsratio2, ismissing(obs_rsds), False)
  corr_obs_rsdsratio2!0   = "lat"
  corr_obs_rsdsratio2!1   = "lon"
  corr_obs_rsdsratio2&lat = lat
  corr_obs_rsdsratio2&lon = lon
  corr_obs_rsdsratio2@comment = "Observation, ratio=cLeaf/cWood, preindustrial"
  corr_obs_rsdsratio2@long_name = "Correlation of shortwave radiation and cLeaf/cWood"
  outfile->corr_obs_rsdsratio2 = corr_obs_rsdsratio2
  delete(corr_obs_rsdsratio2)
  delete(temp_corr_obs_rsdsratio2)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------
;;; RSDS vs ratio3 (Bbio/cWood)

  corr_obs_rsdsratio3 = new((/dimsizes(lat),dimsizes(lon)/),float)
  temp_corr_obs_rsdsratio3 = new((/dimsizes(lat),20/),float)

  temp_array1 = new((/dimsizes(lat),20/),float)
  temp_array2 = new((/dimsizes(lat),20/),float)

  do i = 0, dimsizes(lat)-1
    do j = 0,9                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,dimsizes(lon)-10+j)
      temp_array2(i,j) = obs_ratio3(i,dimsizes(lon)-10+j)
    end do
    do j = 10,19                ;; longitude boundary points
      temp_array1(i,j) = obs_rsds(i,j-10)
      temp_array2(i,j) = obs_ratio3(i,j-10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:dimsizes(lat)-1,j-5:j+5))) .or. \
            all(ismissing(temp_array2(i-5:dimsizes(lat)-1,j-5:j+5))) ) then
          temp_corr_obs_rsdsratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsratio3(i,j) = pattern_cor( temp_array1(i-5:dimsizes(lat)-1,j-5:j+5), \
                                                    temp_array2(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,14  ;; longitude
        if ( all(ismissing(temp_array1(i-5:i+5,j-5:j+5))) .or. \
             all(ismissing(temp_array2(i-5:i+5,j-5:j+5))) ) then
          temp_corr_obs_rsdsratio3(i,j) = 9.96921e+36
        else
          temp_corr_obs_rsdsratio3(i,j) = pattern_cor( temp_array1(i-5:i+5,j-5:j+5), \
                                                    temp_array2(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  do i=5, dimsizes(lat)-1
    do j=0, 4
      corr_obs_rsdsratio3(i,j) = temp_corr_obs_rsdsratio3(i,10+j)
    end do
    do j=dimsizes(lon)-5,dimsizes(lon)-1
      corr_obs_rsdsratio3(i,j) = temp_corr_obs_rsdsratio3(i,j-dimsizes(lon)+10)
    end do
  end do

  do i = 5, dimsizes(lat)-1
    if (i.gt.dimsizes(lat)-6) then
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5) ))) then
          corr_obs_rsdsratio3(i,j) = 9.96921e+36
        else
          corr_obs_rsdsratio3(i,j) = pattern_cor( obs_rsds(i-5:dimsizes(lat)-1,j-5:j+5), \
                                               obs_ratio3(i-5:dimsizes(lat)-1,j-5:j+5),1.0,0 )
        end if
      end do
    else
      do j = 5,dimsizes(lon)-6  ;; longitude
        if ( all(ismissing( obs_rsds(i-5:i+5,j-5:j+5) ))) then
          corr_obs_rsdsratio3(i,j) = 9.96921e+36
        else
          corr_obs_rsdsratio3(i,j) = pattern_cor( obs_rsds(i-5:i+5,j-5:j+5),obs_ratio3(i-5:i+5,j-5:j+5),1.0,0 )
        end if
      end do
    end if
  end do

  corr_obs_rsdsratio3 = mask(corr_obs_rsdsratio3, ismissing(obs_rsds), False)
  corr_obs_rsdsratio3!0   = "lat"
  corr_obs_rsdsratio3!1   = "lon"
  corr_obs_rsdsratio3&lat = lat
  corr_obs_rsdsratio3&lon = lon
  corr_obs_rsdsratio3@comment = "Observation, ratio=Bbio/cWood, preindustrial"
  corr_obs_rsdsratio3@long_name = "Correlation of shortwave radiation and Bbio/cWood"
  outfile->corr_obs_rsdsratio3 = corr_obs_rsdsratio3
  delete(corr_obs_rsdsratio3)
  delete(temp_corr_obs_rsdsratio3)
  delete(temp_array1)
  delete(temp_array2)

;----------------------------------------------------------

;;; Delete variables

  delete(obs_pr)
  delete(obs_rsds)
  delete(obs_tas)
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

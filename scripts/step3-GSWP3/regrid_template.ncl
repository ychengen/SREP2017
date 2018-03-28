;*************************************************
; regrid.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

  modelnm= "modelname"
  indir1 = "../../CMIP5/climate_mean/modern/"
  fname1 = indir1+"climate_"+modelnm+".nc"
  infile1= addfile(fname1,"r")

  indir2 = "../../Observation/GSWP3/"
  fname21= indir2+"regrid-GSWP3_pr.nc"
  fname22= indir2+"regrid-GSWP3_rsds.nc"
  fname23= indir2+"regrid-GSWP3_tas.nc"
  infile21= addfile(fname21,"r")
  infile22= addfile(fname22,"r")
  infile23= addfile(fname23,"r")

  outdir = "../../Observation/GSWP3/annual_mean/"
  fname3 = outdir+"GSWP3_"+modelnm+".nc"
  system("rm -f " + fname3)
  outfile= addfile(fname3,"c")

  time  = infile21->time
  lat   = infile21->lat
  temp_lon = infile21->lon

  temp_pr = infile21->pr
  temp_rsds = infile22->rsds
  temp_tas = infile23->tas

  lon     = temp_lon+180.0
  lon!0   = "lon"
  lon&lon = lon

  pr   = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)
  rsds = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)
  tas  = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)

  do i=dimsizes(lon)/2,dimsizes(lon)-1
    pr(:,:,i-360) = temp_pr(:,:,i)
    rsds(:,:,i-360) = temp_rsds(:,:,i)
    tas(:,:,i-360) = temp_tas(:,:,i)
  end do
  do i=0,dimsizes(lon)/2-1
    pr(:,:,360+i) = temp_pr(:,:,i)
    rsds(:,:,360+i) = temp_rsds(:,:,i)
    tas(:,:,360+i) = temp_tas(:,:,i)
  end do

  pr!0 = "time"
  pr!1 = "lat"
  pr!2 = "lon"
  pr&time = time
  pr&lat  = lat
  pr&lon  = lon
  rsds!0 = "time"
  rsds!1 = "lat"
  rsds!2 = "lon"
  rsds&time = time
  rsds&lat  = lat
  rsds&lon  = lon
  tas!0 = "time"
  tas!1 = "lat"
  tas!2 = "lon"
  tas&time = time
  tas&lat  = lat
  tas&lon  = lon

;-------------------------------------------------------------
; interpolate to new grid

  lat2 = infile1->lat
  lon2 = infile1->lon

  newlat = doubletofloat(lat2)
  newlon = doubletofloat(lon2)

  newlat!0 ="lat"
  newlat&lat = newlat
  newlat@units = "degree_north"
  newlat@standard_name = "latitude"
  newlat@long_name = "latitude"

  newlon!0 = "lon"
  newlon&lon = newlon
  newlon@units = "degree_east"
  newlon@standard_name = "longitude"
  newlon@long_name = "longitude"

  newpr   = linint2(pr&lon,pr&lat,pr,True,newlon,newlat,0)
  newrsds = linint2(rsds&lon,rsds&lat,rsds,True,newlon,newlat,0)
  newtas  = linint2(tas&lon,tas&lat,tas,True,newlon,newlat,0)

  newpr!0   = "time"
  newpr!1   = "lat"
  newpr!2   = "lon" 
  newpr&time= time
  newpr&lat = newlat
  newpr&lon = newlon
  newpr@_FillValue = 9.96921e+36
  newpr@missing_value = 9.96921e+36 
  newpr@cell_methods = "time: annual mean"
  newpr@comment = "Rain + Snow"
  newpr@units = "mm/day"
  newpr@standard_name = "precipitation_flux"
  newpr@long_name = "atmospheric rain"

  newrsds!0   = "time"
  newrsds!1   = "lat"
  newrsds!2   = "lon" 
  newrsds&time= time
  newrsds&lat = newlat
  newrsds&lon = newlon
  newrsds@_FillValue = 9.96921e+36
  newrsds@missing_value = 9.96921e+36 
  newrsds@cell_methods = "time: annual mean"
  newrsds@comment = "FSDS unchanged"
  newrsds@units = "W m-2"
  newrsds@original_name = "FSDS"
  newrsds@standard_name = "surface_downwelling_shortwave_flux_in_air"
  newrsds@long_name = "Atmospheric Incident Solar Radition"

  newtas!0   = "time"
  newtas!1   = "lat"
  newtas!2   = "lon" 
  newtas&time= time
  newtas&lat = newlat
  newtas&lon = newlon
  newtas@_FillValue = 9.96921e+36
  newtas@missing_value = 9.96921e+36 
  newtas@cell_methods = "time: annual mean"
  newtas@comment = "TAS unchanged"
  newtas@units = "K"
  newtas@original_name = "TSA"
  newtas@standard_name = "air_temperature"
  newtas@long_name = "Air Temperature"

;;; Output to new NetCDF file

;---Create variable to hold global file attributes
      global = True
      copy_VarAtts(infile21, global)

      if (isatt(infile21,"title")) then
        global@TITLE = "REMAPPED: " + infile21@title
      end if
      global@remap         = "NCL: ESMF_regrid_with_weights (NCL version '" + \
                             get_ncl_version() + "')"
      global@remap_method  = "Bilinear"
      global@creation_date = systemfunc("date")
    
      fileattdef( outfile, global )        ; copy global file attributes

  outfile->time = time
  outfile->lat  = newlat
  outfile->lon  = newlon
  outfile->pr   = newpr
  outfile->rsds = newrsds
  outfile->tas  = newtas

end
exit

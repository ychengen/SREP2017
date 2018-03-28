;;;========================================================
;;; Compute the distance/area for each lat/lon grid
;;;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  case = (/"Tbio"/)
  indir0 = "../Observation/BGI2014/1km/"

;; delta1 and delta2 are based on the given GEOTIFF info from BGI data

  delta1= 0.00999999977639732
  delta2= 0.00499999988814866
  npts  = 0      ;; number of interpolations between two points
  lon0  = new(1, double)
  lon1  = new(1, double)
  lon2  = new(1, double)
  lat0  = new(1, double)
  lat1  = new(1, double)
  lat2  = new(1, double)
  dist1 = new(1, double)
  dist2 = new(1, double)
  
  z = 0

  fname1 = indir0+"BGI_"+case(z)+"_1km"
  infile1= addfile (fname1+".nc","r")

  fname2 = indir0+"BGI_1km_gridarea"
  if(isfilepresent(fname2))
    print("Warning: '" + fname2 + "' exists.")
    print("Will remove it.")
    system("/bin/rm " + fname2)
  end if
  outfile= addfile (fname2+".nc","c")

  vname1 = getfilevarnames(infile1)

  lat = infile1->$vname1(2)$
  lon = infile1->$vname1(3)$
  lat0 = lat(0)-delta2
  lon0 = lon(0)-delta2

  area  = new((/dimsizes(lat),dimsizes(lon)/), double)
  area!0= "lat"
  area!1= "lon"
  area&lat = lat
  area&lon = lon
  area@units    = "square meter"
  area@longname = "cell area"

;;;---------------------------------------------------------------------

  do x=0,dimsizes(lon)-1
    lon1 = lon(0)+x*delta1
    lon2 = lon1+delta1

    do y=0,dimsizes(lat)-1
      lat1 = lat(0)+y*delta1
      lat2 = lat1+delta1

;;;  abs(iu) = 1   ; return the distance in radians
;;;          = 2   ; return the distance in degrees
;;;          = 3   ; return the distance in meters
;;;          = 4   ; return the distance in kilometers

        dist1 = gc_latlon(lat1,lon1,lat2,lon1,npts,3)
        dist2 = gc_latlon(lat1,lon1,lat1,lon2,npts,3)
        area(y,x) = (dist1@spacing) * (dist2@spacing)
;;;        print ("(lat, lon, dy, dx, area) = "+"("+lat(y)+","+lon(x)+
;;;               ","+dy+","+dx+","+dy*dx+")")
    end do

  end do

;;;---------------------------------------------------------------------
;;; Write out data

  outfile->lat = lat
  outfile->lon = lon
  outfile->gridarea = area 

end
exit

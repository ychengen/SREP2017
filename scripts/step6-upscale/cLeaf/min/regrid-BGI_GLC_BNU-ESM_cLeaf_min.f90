!!!=================================================================
!!! Regrid 1km to BNU-ESM resolution for 22 vegetation types - cLeaf
!!!=================================================================

PROGRAM REGRID_TO_ESM

use netcdf
implicit none

integer :: nc_id0, lat_id0, lon_id0
integer :: nc_id1, lat_id1, lon_id1
integer :: nc_id2, lat_id2, lon_id2
integer :: nc_id3, lat_id3, lon_id3
integer :: nc_id4, lat_id4, lon_id4
integer :: latbnd_id4, lonbnd_id4
integer :: var_id0, var_id1, var_id2,veg_id3
integer :: var_id3, var_id4,var_id5,var_id6,var_id7,var_id8, &
           var_id9, var_id10, var_id11, var_id12, var_id13,  &
           var_id14, var_id15
integer :: status, lat_dimid, lon_dimid, veg_dimid

character(len=*), parameter :: FNAME0 =     &
  "/lustre/atlas/scratch/ceyang1/cli017/BGC-Biomass/Observation/BGI2014/1km/BGI_1km_gridarea.nc"
character(len=*), parameter :: FNAME1 =     &
  "/lustre/atlas/scratch/ceyang1/cli017/BGC-Biomass/Observation/BGI2014/1km/BGI_cLeaf_min.nc"
character(len=*), parameter :: FNAME2 =     &
  "/lustre/atlas/scratch/ceyang1/cli017/BGC-Biomass/Observation/GLC2000/GLC2000_BGI_un.nc"
character(len=*), parameter :: FNAME4 =     &
  "/lustre/atlas/scratch/ceyang1/cli017/BGC-Biomass/CMIP5/BNU-ESM/1982-2005/cLeaf_Lmon_BNU-ESM_historical_r1i1p1_198201-200512.nc"
character(len=*), parameter :: FOUTNAME =   &
  "/lustre/atlas/scratch/ceyang1/cli017/BGC-Biomass/CMIP5/landcover_GLC2000_for_CMIP5/cLeaf/min/landcover_BGI_cLeaf_min_BNU-ESM.nc"

integer, parameter :: NDIMS=3, NTYPES=22, NLATS1=5000, NLONS1=36000, &
                      NLATS2=64, NLONS2=128, NTIME2=288
integer :: dimid1(NDIMS-1), dimid2(NDIMS), start(NDIMS), count(NDIMS)
double precision :: lat_in(NLATS1), lon_in(NLONS1), lat_in4(NLATS2), lon_in4(NLONS2), &
                    latbnd_in(2,NLATS2), lonbnd_in(2,NLONS2)
real :: lat_out(NLATS2), lon_out(NLONS2)
double precision, dimension(:,:), allocatable :: varin0, varin1
integer, dimension(:,:), allocatable :: varin2

character(len=*), parameter ::                           &
  LONGNAME="longname", UNITS="units", VEG_NAME="class",  &
  LAT_NAME="lat", LON_NAME="lon", COMMENT1="comment1",   &
  COMMENT2="comment2", COMMENT3="comment3",              &
  COMMENT4="comment4", COMMENT5="comment5",              &
  COMMENT="comment", VEG_UNITS="",                       &
  LAT_LONGNAME="latitude", LON_LONGNAME="longitude",     &
  VEG_LONGNAME="land type", FILLNAME="_FillValue",       &
  LAT_UNITS="degrees_north", LON_UNITS="degrees_east"   

character(len=*), parameter :: &
VEG_COMMENT1=                  &
"1-Tree Cover, broadleaved, evergreen; 2-Tree Cover, broadleaved, deciduous, closed; 3-Tree Cover, broadleaved, deciduous, open; 4-Tree Cover, needle-leaved, evergreen; 5-Tree Cover, needle-leaved, deciduous;"

character(len=*), parameter :: &
VEG_COMMENT2=                  &
"6-Tree Cover, mixed leaf type; 7-Tree Cover, regularly flooded, fresh water; 8-Tree Cover, regularly flooded, saline water; 9-Mosaic: Tree Cover/Other natural vegetation, 10-Tree Cover, burnt;"

character(len=*), parameter :: &
VEG_COMMENT3=                  &
"11-Shrub Cover, closed/open, evergreen; 12-Shrub Cover, closed/open, deciduous; 13-Herbaceous Cover, closed/open; 14-Sparse herbaceous or sparse shrub cover; 15-Regularly flooded shrub and/or herbaceous cover;"

character(len=*), parameter :: &
VEG_COMMENT4=                  &
"16-Cultivated and managed areas; 17-Mosaic: Cropland/Tree Cover/Other natural vege; 18-Mosaic: Cropland/Shrub and/or grass cover; 19-Bare Areas; 20-Water Bodies;"

character(len=*), parameter :: &
VEG_COMMENT5=                  &
"21-Snow and Ice; 22-Artificial surfaces and associated areas;"

character(len=*), parameter :: CELL_COMMENT = "per cell area"
character(len=*), parameter :: FOREST_COMMENT = "per forest area"

character(len=*), parameter :: LATBND_NAME="lat_bnds",        &
  LONBND_NAME="lon_bnds", VARIN0_NAME="gridarea",             &
  VARIN111_NAME="cLeaf_min",                                  &
  VARIN1_NAME="Band1", VOUT1_NAME="cell_total",               &
  VOUT2_NAME="cell_area", VOUT3_NAME="cell_mean",             &
  VOUT4_NAME="veg_total", VOUT5_NAME="veg_area",              &
  VOUT6_NAME="veg_mean", VOUT7_NAME="veg_frac",               &
  VOUT8_NAME="forest_total", VOUT9_NAME="forest_area",        &
  VOUT10_NAME="forest_mean", VOUT11_NAME="forest_frac",       &
  VOUT12_NAME="veg_grid_mean", VOUT13_NAME="forest_grid_mean",&
  VOUT1_LONGNAME="biomass", VOUT2_LONGNAME="cell area",       &
  VOUT3_LONGNAME="cell mean biomass flux",                    &
  VOUT4_LONGNAME="biomass per land type",                     &
  VOUT5_LONGNAME="land area per land type",                   &
  VOUT6_LONGNAME="mean land type biomass flux",               &
  VOUT12_LONGNAME="mean land type biomass flux per cell area",&
  VOUT7_LONGNAME="fraction per land type",                    &
  VOUT8_LONGNAME="forest type only biomass",                  &
  VOUT9_LONGNAME="forest type only area",                     &
  VOUT10_LONGNAME="mean forest biomass flux" ,                &
  VOUT13_LONGNAME="mean forest biomass flux per cell area",   &
  VOUT11_LONGNAME="fraction per forest type" ,                &
  VOUT1_UNITS="kilogram", VOUT2_UNITS="square meter",         &
  VOUT3_UNITS="kilogram per square meter", VOUT7_UNITS="",    &
  VOUT4_UNITS="kilogram", VOUT5_UNITS="square meter",         &
  VOUT6_UNITS="kilogram per square meter", VOUT11_UNITS="",   &
  VOUT8_UNITS="kilogram",  VOUT9_UNITS="square meter",        &
  VOUT10_UNITS="kilogram per square meter",                   &
  VOUT12_UNITS="kilogram per square meter",                   &
  VOUT13_UNITS="kilogram per square meter"

double precision, parameter :: FILLDVAL=9.969209968386869e+36

integer :: i, j, p, x, y
integer :: veg_type(NTYPES)
double precision :: cell_area(NLONS2,NLATS2),                 &
  forest_area(NLONS2,NLATS2), veg_area(NLONS2,NLATS2,NTYPES), &
  cell_total(NLONS2,NLATS2), veg_total(NLONS2,NLATS2,NTYPES), &
  cell_mean(NLONS2,NLATS2), veg_mean(NLONS2,NLATS2,NTYPES),   &
  veg_frac(NLONS2,NLATS2,NTYPES), forest_total(NLONS2,NLATS2),&
  forest_mean(NLONS2,NLATS2), forest_frac(NLONS2,NLATS2),     &
  grid_mean(NLONS2,NLATS2), veg_grid_mean(NLONS2,NLATS2,NTYPES)

!!!---------------------------------------------------------------------

!!! Allocate memory for input data (note: column, row)
allocate(varin0(NLONS1, NLATS1))        !! grid area
allocate(varin1(NLONS1, NLATS1))        !! observational biomass
allocate(varin2(NLONS1, NLATS1))        !! landcover type

!!! Check input file and variable names
status = nf90_open(FNAME0, nf90_NoWrite, nc_id0)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_open(FNAME1, nf90_NoWrite, nc_id1)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_open(FNAME2, nf90_NoWrite, nc_id2)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_open(FNAME4, nf90_NoWrite, nc_id4)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_inq_varid(nc_id0, LAT_NAME, lat_id0)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id0, LON_NAME, lon_id0)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id0, VARIN0_NAME, var_id0)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id1, LAT_NAME, lat_id1)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id1, LON_NAME, lon_id1)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id1, VARIN111_NAME, var_id1)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id2, LAT_NAME, lat_id2)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id2, LON_NAME, lon_id2)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id2, VARIN1_NAME, var_id2)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id4, LAT_NAME, lat_id4)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id4, LON_NAME, lon_id4)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id4, LATBND_NAME, latbnd_id4)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_id4, LONBND_NAME, lonbnd_id4)
if(status /= nf90_NoErr) call handle_err(status)

!!! Get the latitude, longitude, and variable data from inputs
status = nf90_get_var(nc_id1, lat_id1, lat_in)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(nc_id1, lon_id1, lon_in)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(nc_id4, lat_id4, lat_in4)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(nc_id4, lon_id4, lon_in4)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(nc_id4, latbnd_id4, latbnd_in)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(nc_id4, lonbnd_id4, lonbnd_in)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_get_var(nc_id0, var_id0, varin0)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(nc_id1, var_id1, varin1)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_var(nc_id2, var_id2, varin2)
if(status /= nf90_NoErr) call handle_err(status)

print *,"... Stage 1 done ..."

!!!---------------------------------------------------------------------
!!! Begin regridding

!!! Define new NLATS2xNLONS2 lat/lon and 22 vegetation types

DO i=1,NLONS2
  lon_out(i) = lon_in4(i)
ENDDO
DO j=1,NLATS2
  lat_out(j) = lat_in4(j)
ENDDO
DO p=1,22
  veg_type(p) = p
ENDDO

!!! Initialize variables with missing values
DO j=1,NLATS2
  DO i=1,NLONS2
    cell_total(i,j)  = FILLDVAL
    cell_area(i,j)   = FILLDVAL
    cell_mean(i,j)   = FILLDVAL    !! for all land types 
    forest_total(i,j)= FILLDVAL
    forest_area(i,j) = FILLDVAL
    forest_mean(i,j) = FILLDVAL
    forest_frac(i,j) = FILLDVAL    !! for forest only and per forest area
    grid_mean(i,j)   = FILLDVAL    !! for forest only and per cell area

    DO p=1,NTYPES
      veg_total(i,j,p)= FILLDVAL
      veg_area(i,j,p) = FILLDVAL
      veg_mean(i,j,p) = FILLDVAL        !! per forest area
      veg_frac(i,j,p) = FILLDVAL
      veg_grid_mean(i,j,p) = FILLDVAL   !! per cell area
    ENDDO
  ENDDO
ENDDO

print *,"... Stage 2 done ..."

!!!---------------------------------------------------------------------
!!! Find the fitted data and assign its data to the specific 0.5 degree

!DO j=241,340              ! NLATS for 0.5 degree
DO j=1,NLATS2             ! NLATS for model resolution

 IF ( lat_in4(j).ge.30.0 .and. lat_in4(j).lt.80.0) THEN   ! Check if the latitude resides in the larger grid range

  DO i=1,NLONS2           ! NLONS for model resolution

!    print *, 'test', i, j, lat_out(j), lon_out(i)

    cell_area(i,j) = 0.0
    cell_total(i,j) = 0.0
    forest_area(i,j) = 0.0
    forest_total(i,j) = 0.0

    DO p=1,22
      veg_area(i,j,p) = 0.0
      veg_total(i,j,p) = 0.0
    ENDDO

!!!! Original 0.01 degree check, aggregate to model grid size

!    DO y=50*(j-241)+1,50*(j-241)+50
    DO y=1,NLATS1

     IF (lat_in(y).ge.latbnd_in(1,j) .and. lat_in(y).lt.latbnd_in(2,j)) THEN

!      DO x=50*i+1,50*i+50
      DO x=1,NLONS1
        IF (lon_in(x).lt.0.) THEN
          lon_in(x) = 360.0+lon_in(x)    !! convert -180-0-180 to 0-360
        ENDIF
!write (*,*) i,"=",lon_in4(i),",",j,"=",lat_in4(j),",",x,"=",lon_in(x),",",y,"=",lat_in(y)

       IF (i.eq.1 .and. lon_in(x).ge.lonbnd_in(2,NLONS2)) THEN   !! add toward the first lon grid
        cell_area(i,j) = cell_area(i,j)+varin0(x,y)            ! m^2

!!! modified by CYZ
!        IF ( varin1(x,y)>0.0 .and. varin1(x,y)<1.e+36 ) THEN
        IF ( varin1(x,y)> -3.4e+38 .and. varin1(x,y)<1.e+36 ) THEN

!!!! For a whole half degree grid biomass flux, all land types, kg
          cell_total(i,j)= cell_total(i,j)+varin1(x,y)*varin0(x,y)

!!!! For 22 land types
          IF ( int(varin2(x,y))>=1 .and. int(varin2(x,y))<=22 ) THEN
            veg_area(i,j,int(varin2(x,y))) =  &
                              veg_area(i,j,int(varin2(x,y)))+varin0(x,y)
            veg_total(i,j,int(varin2(x,y)))=  &
                 veg_total(i,j,int(varin2(x,y)))+varin1(x,y)*varin0(x,y)
          ELSE
            CONTINUE
          ENDIF

!!!! For forest types (1-10) only
          IF ( int(varin2(x,y))>=1 .and. int(varin2(x,y))<=10 ) THEN
            forest_area(i,j) = forest_area(i,j)+varin0(x,y)
            forest_total(i,j)= forest_total(i,j)   &
                             + varin1(x,y)*varin0(x,y)
          ELSE
            CONTINUE
          ENDIF

        ELSE
          CONTINUE
        ENDIF

       ELSEIF (lon_in(x).ge.lonbnd_in(1,i) .and. lon_in(x).lt.lonbnd_in(2,i) ) THEN

        cell_area(i,j) = cell_area(i,j)+varin0(x,y)            ! m^2

!!! modified by CYZ
!        IF ( varin1(x,y)>0.0 .and. varin1(x,y)<1.e+36 ) THEN
        IF ( varin1(x,y)> -3.4e+38 .and. varin1(x,y)<1.e+36 ) THEN

!!!! For a whole half degree grid biomass flux, all land types, kg
          cell_total(i,j)= cell_total(i,j)+varin1(x,y)*varin0(x,y)

!!!! For 22 land types
          IF ( int(varin2(x,y))>=1 .and. int(varin2(x,y))<=22 ) THEN
            veg_area(i,j,int(varin2(x,y))) =  &
                              veg_area(i,j,int(varin2(x,y)))+varin0(x,y)
            veg_total(i,j,int(varin2(x,y)))=  &
                 veg_total(i,j,int(varin2(x,y)))+varin1(x,y)*varin0(x,y)
          ELSE
            CONTINUE
          ENDIF

!!!! For forest types (1-10) only
          IF ( int(varin2(x,y))>=1 .and. int(varin2(x,y))<=10 ) THEN
            forest_area(i,j) = forest_area(i,j)+varin0(x,y)
            forest_total(i,j)= forest_total(i,j)   &
                             + varin1(x,y)*varin0(x,y)
          ELSE
            CONTINUE
          ENDIF

        ELSE
          CONTINUE
        ENDIF

       ELSE
         CONTINUE
       ENDIF

      ENDDO

     ENDIF

    ENDDO
  ENDDO

 ENDIF

ENDDO                  ! NLATS for model resolution

print *,"... Stage 3 done ..."

!write (*,999)

!DO j=241,340           ! NLATS for 0.5 degree, 360
DO j=1,NLATS2          ! NLATS for model resolution
  DO i=1,NLONS2      ! NLONS for 0.5 degree, 720

!! For whole grid mean flux, all land types
!!! modified by CYZ
!    IF ( cell_area(i,j)>0.0 .and. cell_total(i,j)>0.0 .and. &
    IF ( cell_area(i,j)>0.0 .and. cell_total(i,j)> -3.4e+38 .and. &
         cell_area(i,j).ne.FILLDVAL .and. cell_total(i,j).ne.FILLDVAL) THEN
      cell_mean(i,j) = cell_total(i,j)/cell_area(i,j)
    ELSE
      cell_total(i,j)= FILLDVAL
      cell_area(i,j) = FILLDVAL
      cell_mean(i,j) = FILLDVAL
    ENDIF

!! Forest only mean flux
!!! Per cell area
!!! modified by CYZ
!    IF ( cell_area(i,j)>0.0 .and. forest_total(i,j)>0.0 .and. &
    IF ( cell_area(i,j)>0.0 .and. forest_total(i,j)> -3.4e+38 .and. &
         cell_area(i,j).ne.FILLDVAL .and. forest_total(i,j).ne.FILLDVAL) THEN
      grid_mean(i,j) = forest_total(i,j)/cell_area(i,j)
    ELSE
      grid_mean(i,j) = FILLDVAL
    ENDIF

!!! Per forest area
    IF ( cell_area(i,j)>0.0 .and. forest_area(i,j)>0.0 .and. &
!!! modified by CYZ
!         forest_total(i,j)>0.0 .and. cell_area(i,j).ne.FILLDVAL .and. &
         forest_total(i,j)> -3.4e+38 .and. cell_area(i,j).ne.FILLDVAL .and. &
         forest_total(i,j).ne.FILLDVAL .and. forest_area(i,j).ne.FILLDVAL) THEN
      forest_mean(i,j) = forest_total(i,j)/forest_area(i,j)
      forest_frac(i,j) = forest_area(i,j)/cell_area(i,j)
    ELSE
      forest_total(i,j)= FILLDVAL
      forest_area(i,j) = FILLDVAL
      forest_mean(i,j) = FILLDVAL
      forest_frac(i,j) = FILLDVAL
    ENDIF

    DO p=1,NTYPES
      IF ( cell_area(i,j)>0.0 .and. veg_area(i,j,p)>0.0 .and.  &
!!! modified by CYZ
!           veg_total(i,j,p)>0.0 .and. cell_area(i,j).ne.FILLDVAL .and. &
           veg_total(i,j,p)> -3.4e+38 .and. cell_area(i,j).ne.FILLDVAL .and. &
           veg_area(i,j,p).ne.FILLDVAL .and.veg_total(i,j,p).ne.FILLDVAL) THEN

!! For specific land type
!!! Per cell area
        veg_grid_mean(i,j,p) = veg_total(i,j,p)/cell_area(i,j)

!!! Per forest area
        veg_mean(i,j,p) = veg_total(i,j,p)/veg_area(i,j,p)
        veg_frac(i,j,p) = veg_area(i,j,p)/cell_area(i,j)
!          write (*,1000) i, j, p, cell_mean(i,j), cell_total(i,j), &
!                         cell_area(i,j), cell_size(i,j), &
!                         veg_mean(i,j,p), veg_total(i,j,p), &
!                         veg_area(i,j,p),veg_frac(i,j,p)
      ELSE
        veg_total(i,j,p)= FILLDVAL
        veg_area(i,j,p) = FILLDVAL
        veg_mean(i,j,p) = FILLDVAL
        veg_grid_mean(i,j,p) = FILLDVAL
        veg_frac(i,j,p) = FILLDVAL
      ENDIF
   ENDDO
  ENDDO
ENDDO

!999  FORMAT('LON',',','LAT',',','TYPE',',',4X,'CELL_AVE',',',   &
!            4X,'CELL_SUM',',','CELL_NUM',',','GRID_NUM',',',    &
!            5X,'VEG_AVE',',',5X,'VEG_SUM',',','VEG_NUM',',',    &
!            4X,'VEG_FRAC')
!
!1000 FORMAT(I3,',',I3,',',I4,',',D12.7,',',D12.7,',',I8,',',I8, &
!            ',',D12.7,',',D12.7,',',I7,',',D12.7)

print *,"... Stage 4 done ..."

!!! Close files to free up any internal netCDF resources associated

status = nf90_close(nc_id0)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_close(nc_id1)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_close(nc_id2)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_close(nc_id4)
if(status /= nf90_NoErr) call handle_err(status)

!!!---------------------------------------------------------------------

!! Write out results
!! Create output file
status = nf90_create(FOUTNAME, NF90_CLOBBER, nc_id3)
if(status /= nf90_NoErr) call handle_err(status)

!! Create dimensions
status = nf90_def_dim(nc_id3, LAT_NAME, NLATS2, lat_dimid)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_def_dim(nc_id3, LON_NAME, NLONS2, lon_dimid)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_def_dim(nc_id3, VEG_NAME, NTYPES, veg_dimid)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_def_var(nc_id3, LAT_NAME, NF90_REAL, lat_dimid, lat_id3)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_def_var(nc_id3, LON_NAME, NF90_REAL, lon_dimid, lon_id3)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_def_var(nc_id3, VEG_NAME, NF90_INT, veg_dimid, veg_id3)
if(status /= nf90_NoErr) call handle_err(status)

dimid1 = (/ lon_dimid, lat_dimid /)
dimid2 = (/ lon_dimid, lat_dimid, veg_dimid /)

!! For cell sum cLeaf = cell_total(i,j)
status = nf90_def_var(nc_id3, VOUT1_NAME, NF90_DOUBLE, dimid1, var_id3)
if(status /= nf90_NoErr) call handle_err(status)

!! For cell area cLeaf = cell_area(i,j)
status = nf90_def_var(nc_id3, VOUT2_NAME, NF90_DOUBLE, dimid1, var_id4)
if(status /= nf90_NoErr) call handle_err(status)

!! For cell mean cLeaf = cell_mean(i,j)
status = nf90_def_var(nc_id3, VOUT3_NAME, NF90_DOUBLE, dimid1, var_id5)
if(status /= nf90_NoErr) call handle_err(status)

!! For each land type sum cLeaf = veg_total(i,j,p)
status = nf90_def_var(nc_id3, VOUT4_NAME, NF90_DOUBLE, dimid2, var_id6)
if(status /= nf90_NoErr) call handle_err(status)

!! For each land type area cLeaf = veg_area(i,j,p)
status = nf90_def_var(nc_id3, VOUT5_NAME, NF90_DOUBLE, dimid2, var_id7)
if(status /= nf90_NoErr) call handle_err(status)

!! For each land type mean cLeaf, per forest area = veg_mean(i,j,p)
status = nf90_def_var(nc_id3, VOUT6_NAME, NF90_DOUBLE, dimid2, var_id8)
if(status /= nf90_NoErr) call handle_err(status)

!! For each land type mean cLeaf, per cell area = veg_grid_mean(i,j,p)
status = nf90_def_var(nc_id3, VOUT12_NAME, NF90_DOUBLE, dimid2, var_id14)
if(status /= nf90_NoErr) call handle_err(status)

!! Grid fraction for each vegetation type = veg_frac(i,j,p)
status = nf90_def_var(nc_id3, VOUT7_NAME, NF90_DOUBLE, dimid2, var_id9)
if(status /= nf90_NoErr) call handle_err(status)

!! For forest type only sum cLeaf = forest_total(i,j)
status = nf90_def_var(nc_id3, VOUT8_NAME, NF90_DOUBLE, dimid1, var_id10)
if(status /= nf90_NoErr) call handle_err(status)

!! For forest type only area cLeaf = forest_area(i,j)
status = nf90_def_var(nc_id3, VOUT9_NAME, NF90_DOUBLE, dimid1, var_id11)
if(status /= nf90_NoErr) call handle_err(status)

!! For forest type only mean cLeaf, per forest area = forest_mean(i,j)
status = nf90_def_var(nc_id3, VOUT10_NAME, NF90_DOUBLE, dimid1, var_id12)
if(status /= nf90_NoErr) call handle_err(status)

!! For forest type only mean cLeaf, per cell area = grid_mean(i,j)
status = nf90_def_var(nc_id3, VOUT13_NAME, NF90_DOUBLE, dimid1, var_id15)
if(status /= nf90_NoErr) call handle_err(status)

!! Grid fraction of forest only type = forest_frac(i,j)
status = nf90_def_var(nc_id3, VOUT11_NAME, NF90_DOUBLE, dimid1, var_id13)
if(status /= nf90_NoErr) call handle_err(status)

!! Put attributes
status = nf90_put_att(nc_id3, lat_id3, LONGNAME, LAT_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, lon_id3, LONGNAME, LON_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, veg_id3, LONGNAME, VEG_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, lat_id3, UNITS, LAT_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, lon_id3, UNITS, LON_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, veg_id3, UNITS, VEG_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, veg_id3, COMMENT1, VEG_COMMENT1)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, veg_id3, COMMENT2, VEG_COMMENT2)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, veg_id3, COMMENT3, VEG_COMMENT3)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, veg_id3, COMMENT4, VEG_COMMENT4)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, veg_id3, COMMENT5, VEG_COMMENT5)
if(status /= nf90_NoErr) call handle_err(status)

status =nf90_put_att(nc_id3, var_id3, LONGNAME, VOUT1_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id4, LONGNAME, VOUT2_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id5, LONGNAME, VOUT3_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id6, LONGNAME, VOUT4_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id7, LONGNAME, VOUT5_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id8, LONGNAME, VOUT6_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id14, LONGNAME, VOUT12_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id9, LONGNAME, VOUT7_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id10, LONGNAME, VOUT8_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id11, LONGNAME, VOUT9_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id12, LONGNAME, VOUT10_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id15, LONGNAME, VOUT13_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id13, LONGNAME, VOUT11_LONGNAME)
if(status /= nf90_NoErr) call handle_err(status)

status =nf90_put_att(nc_id3, var_id3, UNITS, VOUT1_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id4, UNITS, VOUT2_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id5, UNITS, VOUT3_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id6, UNITS, VOUT4_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id7, UNITS, VOUT5_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id8, UNITS, VOUT6_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id14, UNITS, VOUT12_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id9, UNITS, VOUT7_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id10, UNITS, VOUT8_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id11, UNITS, VOUT9_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id12, UNITS, VOUT10_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id15, UNITS, VOUT13_UNITS)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id13, UNITS, VOUT11_UNITS)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_put_att(nc_id3, var_id8, COMMENT, FOREST_COMMENT)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, var_id12, COMMENT, FOREST_COMMENT)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, var_id14, COMMENT, CELL_COMMENT)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_att(nc_id3, var_id15, COMMENT, CELL_COMMENT)
if(status /= nf90_NoErr) call handle_err(status)

status =nf90_put_att(nc_id3, var_id3, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id4, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id5, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id6, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id7, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id8, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id14, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id9, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id10, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id11, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id12, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id15, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)
status =nf90_put_att(nc_id3, var_id13, FILLNAME, FILLDVAL)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_enddef(nc_id3)
if(status /= nf90_NoErr) call handle_err(status)

print *,"... Stage 5 done ..."

!!!---------------------------------------------------------------------

!! Write out variables values
status = nf90_put_var(nc_id3, lat_id3, lat_out)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, lon_id3, lon_out)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, veg_id3, veg_type)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, var_id3, cell_total)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, var_id4, cell_area)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, var_id5, cell_mean)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, var_id10, forest_total)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, var_id11, forest_area)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, var_id12, forest_mean)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, var_id15, grid_mean)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3, var_id13, forest_frac)
if(status /= nf90_NoErr) call handle_err(status)

start = (/ 1, 1, 1 /)
count = (/ NLONS2, NLATS2, NTYPES /)
status = nf90_put_var(nc_id3,var_id6,veg_total,start=start,count=count)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3,var_id7,veg_area,start=start,count=count)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3,var_id8,veg_mean,start=start,count=count)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3,var_id14,veg_grid_mean,start=start,count=count)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_put_var(nc_id3,var_id9,veg_frac,start=start,count=count)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_close(nc_id3)
if(status /= nf90_NoErr) call handle_err(status)

print *,"... Stage 6 done ...Finish."

END

SUBROUTINE handle_err(status)
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  include 'netcdf.inc'

  INTEGER status
  CHARACTER*80 err_message

  err_message=NF_STRERROR(status)
  WRITE (nuerr,*) 'Error:', err_message

  STOP
END

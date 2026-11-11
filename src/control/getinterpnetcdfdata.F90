module getinterpnetcdfdata

! Description:
!   Routines for extracting a column from a netcdf file
!
! Author:
!
! Modules Used:
!
  use cam_abortutils, only: endrun
  use pmgrid,         only: plev
  use cam_logfile,    only: iulog
  use shr_kind_mod,   only: r8_mod => shr_kind_r8

  implicit none
  private
!
! Public Methods:
!
  public getinterpncdata

  ! The column for each field read, found once per field rather than at every
  ! read (the search reads the coordinate arrays).  Keyed by file, field and
  ! position; the column depends on the field's own horizontal dimensions
  ! (e.g. ncol_d or ncol on a CAMIOP file; lon/lat on an observed IOP).
  integer, parameter :: ncache = 256
  integer,           save :: ncached = 0
  integer,           save :: cache_ncid(ncache), cache_latidx(ncache), cache_lonidx(ncache)
  integer,           save :: cache_nh(ncache)
  real(r8_mod),      save :: cache_lat(ncache), cache_lon(ncache)
  character(len=64), save :: cache_var(ncache), cache_latdim(ncache), cache_hdims(2,ncache)

contains

subroutine getinterpncdata( NCID, camlat, camlon, TimeIdx, &
   varName, have_surfdat, surfdat, fill_ends, scm_crm_mode, &
   press, npress, ps, hyam, hybm, outData, STATUS )

!     getinterpncdata: extracts the entire level dimension for a
!     particular lat,lon,time from a netCDF file
!     and interpolates it onto the input pressure levels, placing
!     result in outData, and the error status inx STATUS

   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use scam_column_mod, only: scam_find_column, scam_horiz_dims
   use netcdf
   implicit none
!-----------------------------------------------------------------------


!     ---------- inputs ------------

   integer, intent(in)  :: NCID          ! NetCDF ID
   integer, intent(in)  :: TimeIdx       ! time index
   real(r8), intent(in) :: camlat,camlon ! target lat and lon to be extracted
   logical, intent(in)  :: have_surfdat  ! is surfdat provided
   logical, intent(in)  :: fill_ends     ! extrapolate the end values
   logical, intent(in)  :: scm_crm_mode  ! scam column radiation mode
   integer, intent(in)  :: npress        ! number of dataset pressure levels
   real(r8), intent(in) :: press(npress) ! dataset pressure levels
   real(r8), intent(in) :: ps            ! surface pressure
   real(r8), intent(in) :: hyam(:)       ! dataset hybrid midpoint pressure levels
   real(r8), intent(in) :: hybm(:)       ! dataset hybrid midpoint pressure levels

!     ---------- outputs ----------

   real(r8), intent(inout)  :: outData(:)    ! interpolated output data
   integer, intent(out)   :: STATUS        ! return status of netcdf calls

!     -------  locals ---------

   real(r8)  surfdat       ! surface value to be added before interpolation
   integer nlev          ! number of levels in dataset
   integer     latIdx        ! latitude index
   integer     lonIdx        ! longitude index
   real(r8),allocatable :: tmp(:)
   real(r8)        closelat,closelon
   real(r8)        dx, dy, m             ! slope for interpolation of surfdat
   integer     varID
   integer     var_ndims
   integer     dims_set
   integer     i
   integer     var_dimIDs( NF90_MAX_VAR_DIMS )
   integer     start( NF90_MAX_VAR_DIMS )
   integer     count( NF90_MAX_VAR_DIMS )

   character   varName*(*)
   character   dim_name*( 256 )
   real(r8)        missing_val
   logical     usable_var
   integer     ic, nh
   logical     hfound
   character(len=64) :: hdims(2), latdim

!     -------  code ---------

   ! the column for this field (see the cache above)
   ic = 0
   do i = 1, ncached
      if (cache_ncid(i) == ncid .and. cache_lat(i) == camlat .and. cache_lon(i) == camlon &
           .and. trim(cache_var(i)) == trim(varName)) then
         ic = i
         exit
      end if
   end do
   if (ic == 0) then
      latidx = 1; lonidx = 1; latdim = ' '; nh = 0
      call scam_horiz_dims(ncid, varName, hdims, nh, hfound)
      if (hfound) then
         call scam_find_column(ncid, hdims(1:nh), camlat, camlon, closelat, closelon, &
              latidx, lonidx, lat_dimname=latdim)
      end if
      ncached = min(ncached + 1, ncache)
      ic = ncached
      cache_ncid(ic) = ncid; cache_lat(ic) = camlat; cache_lon(ic) = camlon
      cache_var(ic) = varName; cache_latidx(ic) = latidx; cache_lonidx(ic) = lonidx
      cache_latdim(ic) = latdim; cache_nh(ic) = nh; cache_hdims(:,ic) = hdims
   end if
   latidx = cache_latidx(ic)
   lonidx = cache_lonidx(ic)
   nh     = cache_nh(ic)
   hdims  = cache_hdims(:,ic)
   latdim = cache_latdim(ic)

!
! Check mode: double or single precision
!

!
! Get var ID.  Check to make sure it's there.
!
   STATUS = NF90_INQ_VARID( NCID, varName, varID )

   if ( STATUS .NE. NF90_NOERR )  return

!
! Check the var variable's information with what we are expecting
! it to be.
!

   STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, ndims=var_ndims )
   if ( var_ndims .GT. 4 ) then
      write(iulog,* ) 'ERROR - extractdata.F: The input var',varName, &
         'has', var_ndims, 'dimensions'
      STATUS = -1
   endif

!
!     surface variables
!
   if ( var_ndims .EQ. 0 ) then
      STATUS = NF90_GET_VAR( NCID, varID, outData )
      return
   endif

   STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, dimids=var_dimIDs )
   if ( STATUS .NE. NF90_NOERR ) then
      write(iulog,* ) 'ERROR - extractdata.F:Cant get dimension IDs for', varName
      return
   endif
!
!     Initialize the start and count arrays
!
   dims_set = 0
   nlev = 1
   do i =  var_ndims, 1, -1

      usable_var = .false.
      STATUS = NF90_INQUIRE_DIMENSION( NCID, var_dimIDs( i ), dim_name )

      ! horizontal dimensions of this field: the latitude dimension of a structured
      ! grid gets latIdx, any other (longitude, or the column of an unstructured grid)
      ! gets lonIdx
      if ( nh == 2 .and. trim(dim_name) == trim(latdim) ) then
         start( i ) =  latIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      else if ( nh >= 1 .and. any(hdims(1:nh) == dim_name) ) then
         start( i ) = lonIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lev' ) then
         STATUS = NF90_INQUIRE_DIMENSION( NCID, var_dimIDs( i ), len=nlev )
         start( i ) = 1
         count( i ) = nlev       ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'ilev' ) then
         STATUS = NF90_INQUIRE_DIMENSION( NCID, var_dimIDs( i ), len=nlev )
         start( i ) = 1
         count( i ) = nlev        ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'time' .OR. dim_name .EQ. 'tsec' ) then
         start( i ) = TimeIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( usable_var .EQV. .false. ) then
         write(iulog,* )'ERROR - extractdata.F: The input var ',varName, &
            ' has an unusable dimension ', dim_name
         STATUS = 1
      endif
   end do

   if ( dims_set .NE. var_ndims ) then
      write(iulog,* )'ERROR - extractdata.F: Could not find all the', &
         ' dimensions for input var ', varName
      write(iulog,* )'Found ',dims_set, ' of ',var_ndims
      STATUS = 1
   endif

   allocate(tmp(nlev+1))

   STATUS = NF90_GET_VAR( NCID, varID, tmp, start, count )

   if ( STATUS .NE. NF90_NOERR ) then
      write(iulog,* )'ERROR - extractdata.F: Could not get data for input var ', varName
      return
   endif

   if ( nlev .eq. 1 ) then
      outdata(1) = tmp(1)
      return                 ! no need to do interpolation
   endif
!   if ( use_camiop .and. nlev.eq.plev) then
   if ( nlev.eq.plev .or. nlev.eq.plev+1) then
      outData(:nlev)= tmp(:nlev)! no need to do interpolation
   else
!
!     add the surface data if available, else
!     fill in missing surface data by extrapolation
!
      if(.not.scm_crm_mode) then
         if ( have_surfdat ) then
            tmp(npress) = surfdat
         else
            dy = press(npress-1) - press(npress-2)
            dx = tmp(npress-1) - tmp(npress-2)
            if ( dx .ne. 0.0_r8 ) then
               m = dy/dx
               tmp(npress) = ((press(npress) - press(npress-1)) / m ) + tmp(npress-1)
            else
               tmp(npress) = tmp(npress-1)
            endif
            surfdat = tmp(npress)
         endif
      endif

#if DEBUG > 1
!
!     check data for missing values
!

      STATUS = NF90_GET_ATT( NCID, varID, 'missing_value', missing_val )
      if ( STATUS .NE. NF90_NOERR ) then
         missing_val = -9999999.0_r8
      endif
!
! reset status to zero
!
      STATUS = 0
!
      do i=1, npress
         if ( tmp(i) .eq. missing_val ) then
            write(iulog,*) 'ERROR - missing value found in ', varname
            write(iulog,*) 'time,lat,lon,lev = ' ,timeidx, latidx, lonidx, i
            stop
         endif
      enddo
#endif
!
      call interplevs( tmp(:npress), press, npress, ps, fill_ends, hyam, hybm, outdata )

   endif

   deallocate(tmp)
   return
 end subroutine getinterpncdata

subroutine interplevs( inputdata,   dplevs,   nlev, &
                       ps, fill_ends, hyam, hybm, outdata)

   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use interpolate_data, only: lininterp
   implicit none

!
!     WARNING: ps, siga and sigb must be initialized before calling this routine
!

!------------------------------Commons----------------------------------

!-----------------------------------------------------------------------


!     ------- inputs -----------
   integer, intent(in) :: nlev                 ! num press levels in dataset

   real(r8), intent(in) :: ps                  ! surface pressure
   real(r8), intent(in) :: hyam(:)             ! a midpoint pressure
   real(r8), intent(in) :: hybm(:)             ! b midpoint pressure
   real(r8), intent(in) :: inputdata(nlev)     ! data from netcdf dataset
   real(r8), intent(in) :: dplevs(nlev)        ! input data pressure levels

   logical, intent(in) :: fill_ends            ! fill in missing end values(used for
                                               ! global model datasets)


! ------- outputs ----------
   real(r8), intent(inout) :: outdata(:)      ! interpolated column data

! ------- locals -----------

   real(r8) mplevs( PLEV )
   real(r8) interpdata( PLEV )


   integer dstart_lev, dend_lev
   integer mstart_lev, mend_lev
   integer data_nlevs, model_nlevs, i
   integer STATUS

!
!     Initialize  model_pressure_levels.  ps should be set in the calling
!     routine to the value in the dataset
!
   do i = 1, plev
      mplevs( i ) = 1000.0_r8 * hyam( i ) + ps * hybm( i ) / 100.0_r8
   end do
!
!     the following algorithm assumes that pressures are increasing in the
!     arrays
!
!
!     Find the data pressure levels that are just outside the range
!     of the model pressure levels, and that contain valid values
!
   dstart_lev = 1
   do i= 1, nlev
      if ( dplevs(i) .LE. mplevs(1) ) dstart_lev  = i
   end do

   dend_lev = nlev
   do i= nlev, 1, -1
      if ( dplevs(i) .GE. mplevs(plev) ) then
         dend_lev  = i
      endif
   end do
!
!     Find the model pressure levels that are just inside the range
!     of the data pressure levels
!
   mstart_lev = 1
   do i=plev, 1, -1
      if ( mplevs( i ) .GE. dplevs( dstart_lev ) )  mstart_lev = i
   end do

   mend_lev = plev
   do i=1,plev
      if ( mplevs( i ) .LE. dplevs( dend_lev ) ) mend_lev = i
   end do

   data_nlevs = dend_lev - dstart_lev +1
   model_nlevs = mend_lev - mstart_lev +1

   call lininterp (inputdata(dstart_lev:dend_lev),dplevs(dstart_lev:dend_lev),data_nlevs, &
		   interpdata,mplevs(mstart_lev:mend_lev),model_nlevs)
!
!     interpolate data onto the model pressure levels
!
!!$   call lininterp (inputdata,dplevs,nlev, &
!!$		   outdata(:plev),mplevs,plev)
   do i=1 , model_nlevs
      outdata( i+mstart_lev-1 ) = interpdata( i )
   end do
!
!     fill in the missing end values
!           (usually  done if this is global model dataset)
!
   if ( fill_ends ) then
      do i=1, mstart_lev
         outdata(i) = inputdata(1)
      end do
      do i= mend_lev, plev
         outdata(i) = inputdata(nlev)
      end do
   end if

   return
end subroutine interplevs
end module getinterpnetcdfdata

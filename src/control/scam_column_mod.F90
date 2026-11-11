module scam_column_mod

!-------------------------------------------------------------------------------
! Single-column (SCAM) nearest-column search.
!
! Given the horizontal dimension(s) of the field a caller is about to read from a
! file, find the point of that field nearest a target latitude/longitude:
!
!   * one dimension (unstructured: ncol, ncol_d, nCells, ...): the latitude and
!     longitude coordinate variables on that dimension; nearest by great-circle
!     distance.  Returns the column as lon_idx, with lat_idx = 1.
!   * two dimensions (structured lat-lon: lon, lat in either order): the coordinate
!     variable on each dimension; nearest latitude and nearest longitude
!     (longitude wrapped across 0/360).  Returns lon_idx and lat_idx.
!
! Coordinate variables are the 1-D variables on a dimension identified as latitude
! or longitude by their standard_name or units attribute (CF), or else by name.
! Values in radians are converted to degrees.  Since the coordinate is taken from
! the field's own dimension, files with several coordinate sets (lat on ncol and
! lat_d on ncol_d, for example) are handled correctly.
!
! Replaces CAM's use of shr_scam_getCloseLatLon (CESM share), which did not know
! which coordinate set belonged to the field being read.
!-------------------------------------------------------------------------------

   use shr_kind_mod,   only: r8 => shr_kind_r8
   use shr_const_mod,  only: pi => shr_const_pi
   use cam_abortutils, only: endrun

   implicit none
   private
   save

   public :: scam_find_column   ! nearest point for a field with the given horizontal dims
   public :: scam_horiz_dims    ! horizontal dimension names of a variable on a file

   interface scam_find_column
      module procedure find_column_nc
      module procedure find_column_pio
   end interface scam_find_column

   interface scam_horiz_dims
      module procedure horiz_dims_nc
      module procedure horiz_dims_pio
   end interface scam_horiz_dims

   integer, parameter :: namelen = 256
   integer, parameter :: KIND_NONE = 0, KIND_LAT = 1, KIND_LON = 2

   ! metadata of the file's 1-D variables (candidates for coordinates)
   type coord_meta_t
      integer :: nvars = 0
      character(len=namelen), allocatable :: name(:)     ! variable name
      character(len=namelen), allocatable :: dimname(:)  ! its dimension (1-D variables only)
      integer,                allocatable :: kind(:)     ! KIND_NONE/LAT/LON
      logical,                allocatable :: radians(:)  ! values in radians
   end type coord_meta_t

!===============================================================================
contains
!===============================================================================

   !----------------------------------------------------------------------------
   ! Classification of a variable as latitude or longitude
   !----------------------------------------------------------------------------
   integer function coord_kind(name, units, stdname)
      character(len=*), intent(in) :: name, units, stdname
      character(len=namelen) :: n, u, s
      n = to_lower(name); u = to_lower(units); s = to_lower(stdname)
      coord_kind = KIND_NONE
      if (trim(s) == 'latitude')  then; coord_kind = KIND_LAT; return; end if
      if (trim(s) == 'longitude') then; coord_kind = KIND_LON; return; end if
      select case (trim(u))
      case ('degrees_north', 'degree_north', 'degrees_n', 'degree_n')
         coord_kind = KIND_LAT; return
      case ('degrees_east', 'degree_east', 'degrees_e', 'degree_e')
         coord_kind = KIND_LON; return
      end select
      select case (trim(n))
      case ('lat', 'lat_d', 'latitude', 'latcell', 'lat_cell', 'yc', 'lsmlat', 'grid1d_lat')
         coord_kind = KIND_LAT
      case ('lon', 'lon_d', 'longitude', 'loncell', 'lon_cell', 'xc', 'lsmlon', 'grid1d_lon')
         coord_kind = KIND_LON
      end select
   end function coord_kind

   logical function units_radians(units)
      character(len=*), intent(in) :: units
      character(len=namelen) :: u
      u = to_lower(units)
      units_radians = (trim(u) == 'radians' .or. trim(u) == 'radian' .or. trim(u) == 'rad')
   end function units_radians

   function to_lower(str) result(low)
      character(len=*), intent(in) :: str
      character(len=namelen) :: low
      integer :: i, c
      low = str
      do i = 1, len_trim(low)
         c = iachar(low(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) low(i:i) = achar(c + 32)
      end do
   end function to_lower

   !----------------------------------------------------------------------------
   ! Choice of the coordinate variable on a dimension (0: none)
   !----------------------------------------------------------------------------
   integer function select_coord(meta, dimname, kind)
      type(coord_meta_t), intent(in) :: meta
      character(len=*),   intent(in) :: dimname
      integer,            intent(in) :: kind
      integer :: n
      select_coord = 0
      do n = 1, meta%nvars
         if (meta%kind(n) == kind .and. trim(meta%dimname(n)) == trim(dimname)) then
            select_coord = n
            return
         end if
      end do
   end function select_coord

   !----------------------------------------------------------------------------
   ! Nearest points
   !----------------------------------------------------------------------------
   integer function nearest_column(lats, lons, tlat, tlon)
      ! great-circle distance (haversine); the first of equal minima wins
      real(r8), intent(in) :: lats(:), lons(:), tlat, tlon
      real(r8) :: d2r, a, amin
      integer  :: i
      d2r = pi/180._r8
      amin = huge(1._r8)
      nearest_column = 0
      do i = 1, size(lats)
         a = sin(0.5_r8*d2r*(lats(i) - tlat))**2 + &
             cos(d2r*tlat)*cos(d2r*lats(i))*sin(0.5_r8*d2r*(lons(i) - tlon))**2
         if (a < amin) then
            amin = a
            nearest_column = i
         end if
      end do
   end function nearest_column

   integer function nearest_lat(lats, tlat)
      real(r8), intent(in) :: lats(:), tlat
      nearest_lat = minloc(abs(lats - tlat), dim=1)
   end function nearest_lat

   integer function nearest_lon(lons, tlon)
      real(r8), intent(in) :: lons(:), tlon
      real(r8) :: d(size(lons))
      d = abs(modulo(lons - tlon + 180._r8, 360._r8) - 180._r8)   ! wrapped |dlon|
      nearest_lon = minloc(d, dim=1)
   end function nearest_lon

   !----------------------------------------------------------------------------
   ! The search, common to both file interfaces
   !----------------------------------------------------------------------------
   subroutine search(dimnames, lat_vals, lon_vals, tlat, tlon, close_lat, close_lon, &
        lat_idx, lon_idx, ok)
      ! lat_vals/lon_vals: the latitude and longitude coordinate values (degrees)
      character(len=*),   intent(in)  :: dimnames(:)
      real(r8),           intent(in)  :: lat_vals(:), lon_vals(:)
      real(r8),           intent(in)  :: tlat, tlon
      real(r8),           intent(out) :: close_lat, close_lon
      integer,            intent(out) :: lat_idx, lon_idx
      logical,            intent(out) :: ok
      ok = .false.
      if (size(dimnames) == 1) then
         if (size(lat_vals) /= size(lon_vals)) return
         lon_idx = nearest_column(lat_vals, lon_vals, tlat, tlon)
         lat_idx = 1
         close_lat = lat_vals(lon_idx)
         close_lon = lon_vals(lon_idx)
      else
         lat_idx = nearest_lat(lat_vals, tlat)
         lon_idx = nearest_lon(lon_vals, tlon)
         close_lat = lat_vals(lat_idx)
         close_lon = lon_vals(lon_idx)
      end if
      ok = .true.
   end subroutine search

   subroutine coord_vars_for(meta, dimnames, ilat_var, ilon_var, latdim, londim)
      ! Latitude and longitude coordinate variables for the field's horizontal dims
      type(coord_meta_t), intent(in)  :: meta
      character(len=*),   intent(in)  :: dimnames(:)
      integer,            intent(out) :: ilat_var, ilon_var
      character(len=*),   intent(out) :: latdim, londim
      integer :: k
      ilat_var = 0; ilon_var = 0; latdim = ' '; londim = ' '
      if (size(dimnames) == 1) then
         ilat_var = select_coord(meta, dimnames(1), KIND_LAT)
         ilon_var = select_coord(meta, dimnames(1), KIND_LON)
         latdim = dimnames(1); londim = dimnames(1)
      else
         do k = 1, size(dimnames)
            if (ilat_var == 0) then
               ilat_var = select_coord(meta, dimnames(k), KIND_LAT)
               if (ilat_var > 0) latdim = dimnames(k)
            end if
            if (ilon_var == 0) then
               ilon_var = select_coord(meta, dimnames(k), KIND_LON)
               if (ilon_var > 0) londim = dimnames(k)
            end if
         end do
      end if
   end subroutine coord_vars_for

   subroutine finish(ok, found, caller, dimnames)
      logical,           intent(in)  :: ok
      logical, optional, intent(out) :: found
      character(len=*),  intent(in)  :: caller
      character(len=*),  intent(in)  :: dimnames(:)
      integer :: k
      character(len=512) :: dims
      if (present(found)) then
         found = ok
      else if (.not. ok) then
         dims = ' '
         do k = 1, size(dimnames)
            dims = trim(dims)//' '//trim(dimnames(k))
         end do
         call endrun(caller//': no latitude/longitude coordinates found on dimension(s)'//trim(dims))
      end if
   end subroutine finish

   !============================================================================
   ! netCDF interface
   !============================================================================
   subroutine meta_nc(ncid, meta)
      use netcdf, only: nf90_inquire, nf90_inquire_variable, nf90_inquire_dimension, &
                        nf90_get_att, nf90_noerr, nf90_max_var_dims
      integer,            intent(in)  :: ncid
      type(coord_meta_t), intent(out) :: meta
      integer :: n, nd, st, dimids(nf90_max_var_dims)
      character(len=namelen) :: units, stdname
      st = nf90_inquire(ncid, nVariables=meta%nvars)
      allocate(meta%name(meta%nvars), meta%dimname(meta%nvars), meta%kind(meta%nvars), &
               meta%radians(meta%nvars))
      do n = 1, meta%nvars
         meta%dimname(n) = ' '; meta%kind(n) = KIND_NONE; meta%radians(n) = .false.
         st = nf90_inquire_variable(ncid, n, name=meta%name(n), ndims=nd, dimids=dimids)
         if (nd /= 1) cycle
         st = nf90_inquire_dimension(ncid, dimids(1), name=meta%dimname(n))
         units = ' '; stdname = ' '
         if (nf90_get_att(ncid, n, 'units', units) /= nf90_noerr) units = ' '
         if (nf90_get_att(ncid, n, 'standard_name', stdname) /= nf90_noerr) stdname = ' '
         meta%kind(n) = coord_kind(meta%name(n), units, stdname)
         meta%radians(n) = units_radians(units)
      end do
   end subroutine meta_nc

   subroutine read_nc(ncid, meta, ivar, vals)
      use netcdf, only: nf90_inquire_variable, nf90_inquire_dimension, nf90_get_var, &
                        nf90_noerr, nf90_max_var_dims
      integer,               intent(in)  :: ncid
      type(coord_meta_t),    intent(in)  :: meta
      integer,               intent(in)  :: ivar
      real(r8), allocatable, intent(out) :: vals(:)
      integer :: st, n, dimids(nf90_max_var_dims)
      st = nf90_inquire_variable(ncid, ivar, dimids=dimids)
      st = nf90_inquire_dimension(ncid, dimids(1), len=n)
      allocate(vals(n))
      st = nf90_get_var(ncid, ivar, vals)
      if (st /= nf90_noerr) call endrun('scam_column_mod: cannot read '//trim(meta%name(ivar)))
      if (meta%radians(ivar)) vals = vals*180._r8/pi
   end subroutine read_nc

   subroutine find_column_nc(ncid, dimnames, target_lat, target_lon, close_lat, close_lon, &
        lat_idx, lon_idx, found, lat_dimname, lon_dimname)
      integer,                    intent(in)  :: ncid         ! netCDF file id
      character(len=*),           intent(in)  :: dimnames(:)  ! horizontal dims of the field
      real(r8),                   intent(in)  :: target_lat, target_lon
      real(r8),                   intent(out) :: close_lat, close_lon
      integer,                    intent(out) :: lat_idx, lon_idx
      logical,          optional, intent(out) :: found
      character(len=*), optional, intent(out) :: lat_dimname, lon_dimname
      type(coord_meta_t)    :: meta
      real(r8), allocatable :: lat_vals(:), lon_vals(:)
      integer :: ilat_var, ilon_var
      logical :: ok
      character(len=namelen) :: latdim, londim
      call meta_nc(ncid, meta)
      call coord_vars_for(meta, dimnames, ilat_var, ilon_var, latdim, londim)
      ok = .false.
      if (ilat_var > 0 .and. ilon_var > 0) then
         call read_nc(ncid, meta, ilat_var, lat_vals)
         call read_nc(ncid, meta, ilon_var, lon_vals)
         call search(dimnames, lat_vals, lon_vals, target_lat, target_lon, close_lat, &
              close_lon, lat_idx, lon_idx, ok)
      end if
      if (present(lat_dimname)) lat_dimname = latdim
      if (present(lon_dimname)) lon_dimname = londim
      call finish(ok, found, 'scam_find_column', dimnames)
   end subroutine find_column_nc

   subroutine horiz_dims_nc(ncid, varname, dimnames, ndims, found)
      ! The dimensions of varname that carry a latitude or longitude coordinate
      use netcdf, only: nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension, &
                        nf90_noerr, nf90_max_var_dims
      integer,           intent(in)  :: ncid
      character(len=*),  intent(in)  :: varname
      character(len=*),  intent(out) :: dimnames(:)
      integer,           intent(out) :: ndims
      logical, optional, intent(out) :: found
      type(coord_meta_t) :: meta
      integer :: vid, nd, k, dimids(nf90_max_var_dims)
      character(len=namelen) :: dname
      ndims = 0
      if (present(found)) found = .false.
      if (nf90_inq_varid(ncid, varname, vid) /= nf90_noerr) return
      call meta_nc(ncid, meta)
      if (nf90_inquire_variable(ncid, vid, ndims=nd, dimids=dimids) /= nf90_noerr) return
      do k = 1, nd
         if (nf90_inquire_dimension(ncid, dimids(k), name=dname) /= nf90_noerr) cycle
         if (select_coord(meta, dname, KIND_LAT) > 0 .or. select_coord(meta, dname, KIND_LON) > 0) then
            ndims = ndims + 1
            if (ndims <= size(dimnames)) dimnames(ndims) = dname
         end if
      end do
      if (present(found)) found = (ndims >= 1 .and. ndims <= size(dimnames))
   end subroutine horiz_dims_nc

   !============================================================================
   ! PIO interface
   !============================================================================
   subroutine meta_pio(file, meta)
      use pio, only: file_desc_t, pio_inquire, pio_inquire_variable, pio_inquire_dimension, &
                     pio_get_att, pio_noerr, pio_seterrorhandling, PIO_BCAST_ERROR
      type(file_desc_t),  intent(inout) :: file
      type(coord_meta_t), intent(out)   :: meta
      integer :: n, nd, st, err_handling, dimids(8), dlen
      character(len=namelen) :: units, stdname
      call pio_seterrorhandling(file, PIO_BCAST_ERROR, err_handling)
      st = pio_inquire(file, nVariables=meta%nvars)
      allocate(meta%name(meta%nvars), meta%dimname(meta%nvars), meta%kind(meta%nvars), &
               meta%radians(meta%nvars))
      do n = 1, meta%nvars
         meta%dimname(n) = ' '; meta%kind(n) = KIND_NONE; meta%radians(n) = .false.
         st = pio_inquire_variable(file, n, name=meta%name(n), ndims=nd)
         if (nd /= 1) cycle
         st = pio_inquire_variable(file, n, dimids=dimids(1:1))
         st = pio_inquire_dimension(file, dimids(1), name=meta%dimname(n), len=dlen)
         units = ' '; stdname = ' '
         if (pio_get_att(file, n, 'units', units) /= pio_noerr) units = ' '
         if (pio_get_att(file, n, 'standard_name', stdname) /= pio_noerr) stdname = ' '
         meta%kind(n) = coord_kind(meta%name(n), units, stdname)
         meta%radians(n) = units_radians(units)
      end do
      call pio_seterrorhandling(file, err_handling)
   end subroutine meta_pio

   subroutine read_pio(file, meta, ivar, vals)
      use pio, only: file_desc_t, pio_inquire_variable, pio_inquire_dimension, pio_get_var, &
                     pio_noerr
      type(file_desc_t),     intent(inout) :: file
      type(coord_meta_t),    intent(in)    :: meta
      integer,               intent(in)    :: ivar
      real(r8), allocatable, intent(out)   :: vals(:)
      integer :: st, n, dimids(1)
      st = pio_inquire_variable(file, ivar, dimids=dimids)
      st = pio_inquire_dimension(file, dimids(1), len=n)
      allocate(vals(n))
      st = pio_get_var(file, ivar, vals)
      if (st /= pio_noerr) call endrun('scam_column_mod: cannot read '//trim(meta%name(ivar)))
      if (meta%radians(ivar)) vals = vals*180._r8/pi
   end subroutine read_pio

   subroutine find_column_pio(file, dimnames, target_lat, target_lon, close_lat, close_lon, &
        lat_idx, lon_idx, found, lat_dimname, lon_dimname)
      use pio, only: file_desc_t
      type(file_desc_t),          intent(inout) :: file
      character(len=*),           intent(in)    :: dimnames(:)
      real(r8),                   intent(in)    :: target_lat, target_lon
      real(r8),                   intent(out)   :: close_lat, close_lon
      integer,                    intent(out)   :: lat_idx, lon_idx
      logical,          optional, intent(out)   :: found
      character(len=*), optional, intent(out)   :: lat_dimname, lon_dimname
      type(coord_meta_t)    :: meta
      real(r8), allocatable :: lat_vals(:), lon_vals(:)
      integer :: ilat_var, ilon_var
      logical :: ok
      character(len=namelen) :: latdim, londim
      call meta_pio(file, meta)
      call coord_vars_for(meta, dimnames, ilat_var, ilon_var, latdim, londim)
      ok = .false.
      if (ilat_var > 0 .and. ilon_var > 0) then
         call read_pio(file, meta, ilat_var, lat_vals)
         call read_pio(file, meta, ilon_var, lon_vals)
         call search(dimnames, lat_vals, lon_vals, target_lat, target_lon, close_lat, &
              close_lon, lat_idx, lon_idx, ok)
      end if
      if (present(lat_dimname)) lat_dimname = latdim
      if (present(lon_dimname)) lon_dimname = londim
      call finish(ok, found, 'scam_find_column', dimnames)
   end subroutine find_column_pio

   subroutine horiz_dims_pio(file, varname, dimnames, ndims, found)
      use pio, only: file_desc_t, pio_inq_varid, pio_inquire_variable, pio_inquire_dimension, &
                     pio_noerr, pio_seterrorhandling, PIO_BCAST_ERROR
      type(file_desc_t), intent(inout) :: file
      character(len=*),  intent(in)    :: varname
      character(len=*),  intent(out)   :: dimnames(:)
      integer,           intent(out)   :: ndims
      logical, optional, intent(out)   :: found
      type(coord_meta_t) :: meta
      integer :: vid, nd, k, st, err_handling, dimids(8)
      character(len=namelen) :: dname
      ndims = 0
      if (present(found)) found = .false.
      call pio_seterrorhandling(file, PIO_BCAST_ERROR, err_handling)
      st = pio_inq_varid(file, varname, vid)
      call pio_seterrorhandling(file, err_handling)
      if (st /= pio_noerr) return
      call meta_pio(file, meta)
      st = pio_inquire_variable(file, vid, ndims=nd)
      st = pio_inquire_variable(file, vid, dimids=dimids(1:nd))
      do k = 1, nd
         st = pio_inquire_dimension(file, dimids(k), name=dname)
         if (select_coord(meta, dname, KIND_LAT) > 0 .or. select_coord(meta, dname, KIND_LON) > 0) then
            ndims = ndims + 1
            if (ndims <= size(dimnames)) dimnames(ndims) = dname
         end if
      end do
      if (present(found)) found = (ndims >= 1 .and. ndims <= size(dimnames))
   end subroutine horiz_dims_pio

end module scam_column_mod

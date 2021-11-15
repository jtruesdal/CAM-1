module scamMod
  !----------------------------------------------------------------------
  !
  ! This module provides control variables and namelist functionality
  ! for SCAM.
  !
  ! As with global CAM, SCAM is initialized with state information
  ! from the initial and boundary files. For each succeeding timestep
  ! SCAM calculates the physics tendencies and combines these with
  ! the advective tendencies provided by the IOP file to produce
  ! a forcast. Generally, the control variables in this module
  ! determine what data and parameterizations will be used to make
  ! the forecast. For instance, many of the control variables in
  ! this module provide flexibility to affect the forecast by overriding
  ! parameterization prognosed tendencies with observed tendencies
  ! of a particular field program recorded on the IOP file.
  ! 
  ! Public functions/subroutines:
  !   scam_readnl
  !-----------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8, r4 => shr_kind_r4
use pmgrid,         only: plon, plat, plev, plevp
use constituents,   only: pcnst
use shr_scam_mod,   only: shr_scam_getCloseLatLon
use dycore,         only: dycore_is
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun

use pio,                 only: file_desc_t, var_desc_t, pio_int, &
                               pio_double, pio_char, pio_inq_dimid, pio_inq_varid, &
                               pio_def_dim, pio_def_var, pio_put_att, pio_noerr, &
                               pio_seterrorhandling, pio_internal_error, pio_bcast_error, &
                               pio_inq_dimlen, pio_put_var, pio_get_var
use spmd_utils,          only: masterproc,npes

implicit none
private

! PUBLIC INTERFACES:

public :: scam_readnl               ! read SCAM namelist options 
public :: scam_restart_init         ! initialize scam restart history data
public :: scam_restart_write        ! write SCAM restart history data
public :: scam_restart_read         ! read SCAM restart history data

integer, parameter :: nlen=26
integer, parameter :: scam_fieldname_len=26
integer, parameter :: errmsglen=132

character(len=errmsglen) :: errmsg

type restart_scamvar_t

   real(r8), pointer              :: v0d => null()
   real(r8), pointer              :: v1d(:) => null()
   real(r8), pointer              :: v2d(:,:) => null()
   real(r8), pointer              :: v3d(:, :, :) => null()
   logical, pointer               :: l0d => null()
   character(len=scam_fieldname_len),   pointer   :: c1d(:) => null()
   character(len=scam_fieldname_len),   pointer   :: c0d => null()

   type(var_desc_t), pointer        :: vdesc => null()
   integer                          :: type
   integer                          :: ndims
   integer                          :: dims(4)
   character(len=nlen)              :: name
end type restart_scamvar_t

type restart_dim_t
   integer                          :: len
   integer                          :: dimid
   character(len=nlen) :: name
   logical                          :: define
end type restart_dim_t

!   The size of these parameters should match the assignments in scam_vars_setnames and scam_dims_setnames below
!
integer, parameter           :: scamvarcnt              = 73
integer, parameter           :: scamdimcnt              = 8
type(restart_scamvar_t)      :: scamvars(scamvarcnt)
type(restart_dim_t)          :: scamdims(scamdimcnt)
integer, parameter           :: latiopbnd_dim_ind      =  1
integer, parameter           :: loniopbnd_dim_ind      =  2
integer, parameter           :: hdimid1_dim_ind        =  3
integer, parameter           :: hdimid2_dim_ind        =  4
integer, parameter           :: lev_dim_ind           =  5
integer, parameter           :: ilev_dim_ind          =  6
integer, parameter           :: pcnst_dim_ind          =  7
integer, parameter           :: scam_fieldname_len_dim_ind  =  8

! PUBLIC MODULE DATA:

real(r8), public ::  pressure_levels(plev)
real(r8), public ::  scmlat   ! input namelist latitude for scam
real(r8), public ::  scmlon   ! input namelist longitude for scam


integer, parameter :: num_switches = 20
integer, parameter :: max_path_len = 128

logical, public ::  single_column         ! Using IOP file or not
logical, public ::  use_iop               ! Using IOP file or not
logical, public ::  use_pert_init         ! perturb initial values
logical, public ::  use_pert_frc          ! perturb forcing 
logical, public ::  switch(num_switches)  ! Logical flag settings from GUI
logical, public ::  l_uvphys              ! If true, update u/v after TPHYS
logical, public ::  l_uvadvect            ! If true, T, U & V will be passed to SLT
logical, public ::  l_conv                ! use flux divergence terms for T and q?     
logical, public ::  l_divtr               ! use flux divergence terms for constituents?
logical, public ::  l_diag                ! do we want available diagnostics?

integer, public ::  error_code            ! Error code from netCDF reads
integer, public ::  initTimeIdx
integer, public ::  seedval

character*(max_path_len), public ::  modelfile
character*(max_path_len), public ::  analysisfile
character*(max_path_len), public ::  sicfile
character*(max_path_len), public ::  userfile
character*(max_path_len), public ::  sstfile
character*(max_path_len), public ::  lsmpftfile
character*(max_path_len), public ::  pressfile
character*(max_path_len), public ::  topofile
character*(max_path_len), public ::  ozonefile
character*(max_path_len), public ::  iopfile
character*(max_path_len), public ::  absemsfile
character*(max_path_len), public ::  aermassfile
character*(max_path_len), public ::  aeropticsfile
character*(max_path_len), public ::  timeinvfile
character*(max_path_len), public ::  lsmsurffile
character*(max_path_len), public ::  lsminifile

! note that scm_zadv_q is set to slt to be consistent with CAM BFB testing


character(len=scam_fieldname_len), public    :: scm_zadv_T  = 'eulc            '
character(len=scam_fieldname_len), public    :: scm_zadv_q  = 'slt             '
character(len=scam_fieldname_len), public    :: scm_zadv_uv = 'eulc            '

real(r8), public ::  fixmascam
real(r8), public ::  betacam
real(r8), public ::  alphacam(pcnst)
real(r8), public ::  dqfxcam(plon,plev,pcnst)

real(r8), public ::      divq3d(plev,pcnst)  ! 3D q advection
real(r8), public ::      divt3d(plev)        ! 3D T advection
real(r8), public ::      divu3d(plev)        ! 3D U advection
real(r8), public ::      divv3d(plev)        ! 3D V advection
real(r8), public ::      vertdivq(plev,pcnst)! vertical q advection
real(r8), public ::      vertdivt(plev)      ! vertical T advection
real(r8), public ::      vertdivu(plev)      ! vertical T advection
real(r8), public ::      vertdivv(plev)      ! vertical T advection
real(r8), public ::      ptend               ! surface pressure tendency
real(r8), public ::      qdiff(plev)         ! model minus observed humidity
real(r8), public ::      qobs(plev)          ! actual W.V. Mixing ratio
real(r8), public ::      qinitobs(plev,pcnst)! initial tracer field
real(r8), public ::      cldliqobs(plev)     ! actual W.V. Mixing ratio
real(r8), public ::      cldiceobs(plev)     ! actual W.V. Mixing ratio
real(r8), public ::      numliqobs(plev)     ! actual 
real(r8), public ::      numiceobs(plev)     ! actual 
real(r8), public ::      precobs(1)          ! observed precipitation 
real(r8), public ::      lhflxobs(1)         ! observed surface latent heat flux 
real(r8), public ::      shflxobs(1)         ! observed surface sensible heat flux
real(r8), public ::      q1obs(plev)         ! observed apparent heat source
real(r8), public ::      q2obs(plev)         ! observed apparent heat sink
real(r8), public ::      tdiff(plev)         ! model minus observed temp 
real(r8), public ::      tground(1)          ! ground temperature
real(r8), public ::      tobs(plev)          ! actual temperature
real(r8), public ::      tsair(1)            ! air temperature at the surface
real(r8), public ::      udiff(plev)         ! model minus observed uwind
real(r8), public ::      uobs(plev)          ! actual u wind
real(r8), public ::      vdiff(plev)         ! model minus observed vwind
real(r8), public ::      vobs(plev)          ! actual v wind
real(r8), public ::      cldobs(plev)        ! observed cld
real(r8), public ::      clwpobs(plev)       ! observed clwp
real(r8), public ::      aldirobs(1)         ! observed aldir
real(r8), public ::      aldifobs(1)         ! observed aldif
real(r8), public ::      asdirobs(1)         ! observed asdir
real(r8), public ::      asdifobs(1)         ! observed asdif

real(r8), public ::      wfld(plev)          ! Vertical motion (slt)
real(r8), public ::      wfldh(plevp)        ! Vertical motion (slt)
real(r8), public ::      divq(plev,pcnst)    ! Divergence of moisture
real(r8), public ::      divt(plev)          ! Divergence of temperature
real(r8), public ::      divu(plev)          ! Horiz Divergence of E/W
real(r8), public ::      divv(plev)          ! Horiz Divergence of N/S
                                             ! mo_drydep algorithm
real(r8), public, pointer :: loniopbnd(:)
real(r8), public, pointer :: latiopbnd(:)
real(r8), public, allocatable :: psinit(:,:) ! initial ps to set etamid/etaint for exact restarts

integer, public ::     iopTimeIdx            ! index into iop dataset
integer, public ::     steplength            ! Length of time-step
integer, public ::     base_date             ! Date in (yyyymmdd) of start time
integer, public ::     base_secs             ! Time of day of start time (sec)

! SCAM public data defaults

logical, public ::  doiopupdate            = .false. ! do we need to read next iop timepoint
logical, public ::  have_lhflx             = .false. ! dataset contains lhflx 
logical, public ::  have_shflx             = .false. ! dataset contains shflx
logical, public ::  have_tg                = .false. ! dataset contains tg
logical, public ::  have_tsair             = .false. ! dataset contains tsair
logical, public ::  have_divq              = .false. ! dataset contains divq 
logical, public ::  have_divt              = .false. ! dataset contains divt
logical, public ::  have_divq3d            = .false. ! dataset contains divq3d 
logical, public ::  have_vertdivu          = .false. ! dataset contains vertdivu
logical, public ::  have_vertdivv          = .false. ! dataset contains vertdivv
logical, public ::  have_vertdivt          = .false. ! dataset contains vertdivt
logical, public ::  have_vertdivq          = .false. ! dataset contains vertdivq 
logical, public ::  have_divt3d            = .false. ! dataset contains divt3d
logical, public ::  have_divu3d            = .false. ! dataset contains divu3d
logical, public ::  have_divv3d            = .false. ! dataset contains divv3d
logical, public ::  have_divu              = .false. ! dataset contains divu
logical, public ::  have_divv              = .false. ! dataset contains divv 
logical, public ::  have_omega             = .false. ! dataset contains omega
logical, public ::  have_phis              = .false. ! dataset contains phis
logical, public ::  have_ptend             = .false. ! dataset contains ptend
logical, public ::  have_ps                = .false. ! dataset contains ps
logical, public ::  have_q                 = .false. ! dataset contains q
logical, public ::  have_q1                = .false. ! dataset contains Q1
logical, public ::  have_q2                = .false. ! dataset contains Q2
logical, public ::  have_prec              = .false. ! dataset contains prec 
logical, public ::  have_t                 = .false. ! dataset contains t
logical, public ::  have_u                 = .false. ! dataset contains u 
logical, public ::  have_v                 = .false. ! dataset contains v 
logical, public ::  have_cld               = .false. ! dataset contains cld
logical, public ::  have_cldliq            = .false. ! dataset contains cldliq
logical, public ::  have_cldice            = .false. ! dataset contains cldice
logical, public ::  have_numliq            = .false. ! dataset contains numliq
logical, public ::  have_numice            = .false. ! dataset contains numice
logical, public ::  have_clwp              = .false. ! dataset contains clwp
logical, public ::  have_aldir             = .false. ! dataset contains aldir
logical, public ::  have_aldif             = .false. ! dataset contains aldif
logical, public ::  have_asdir             = .false. ! dataset contains asdir
logical, public ::  have_asdif             = .false. ! dataset contains asdif
logical, public ::  use_camiop             = .false. ! use cam generated forcing 
logical, public ::  use_3dfrc              = .false. ! use 3d forcing
logical, public ::  isrestart              = .false. ! If this is a restart step or not
  
! SCAM namelist defaults

logical, public ::  scm_backfill_iop_w_init = .false. ! Backfill missing IOP data from initial file
logical, public ::  scm_relaxation         = .false. ! Use relaxation
logical, public ::  scm_crm_mode           = .false. ! Use column radiation mode
logical, public ::  scm_cambfb_mode        = .false. ! Use extra CAM IOP fields to assure bit for bit match with CAM run
logical, public ::  scm_use_obs_T          = .false. ! Use the SCAM-IOP specified observed T   at each time step instead of forecasting.
logical, public ::  scm_force_latlon       = .false. ! force scam to use the lat lon fields specified in the scam namelist not what is closest to iop avail lat lon
real*8, public              ::  scm_relax_top_p       = 1.e36_r8 ! upper bound for scm relaxation
real*8, public              ::  scm_relax_bot_p       = -1.e36_r8 !  lower bound for scm relaxation
real*8, public              ::  scm_relax_tau_sec       = 10800._r8  ! relaxation time constant (sec)

! +++BPM:
! modification... allow a linear ramp in relaxation time scale:
logical, public :: scm_relax_linear = .false.
real*8, public    :: scm_relax_tau_bot_sec = 10800._r8
real*8, public    :: scm_relax_tau_top_sec = 10800._r8
character(len=scam_fieldname_len), public  :: scm_relax_fincl(pcnst)

!
! note that scm_use_obs_uv is set to true to be consistent with CAM BFB testing
!

logical, public ::  scm_use_obs_uv         = .true. ! Use the SCAM-IOP specified observed u,v at each time step instead of forecasting.

logical, public ::  scm_use_obs_qv         = .false. ! Use the SCAM-IOP specified observed qv  at each time step instead of forecasting.
logical, public ::  scm_iop_lhflxshflxTg   = .false. !turn off LW rad
logical, public ::  scm_iop_Tg             = .false. !turn off LW rad

character(len=200), public ::  scm_clubb_iop_name   ! IOP name for CLUBB

!=======================================================================
contains
!=======================================================================

subroutine scam_readnl(nlfile,single_column_in,scmlat_in,scmlon_in)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use dycore,          only: dycore_is
  use wrap_nf,         only: wrap_open
  use spmd_utils,      only : masterproc,npes
  use netcdf,          only : nf90_inquire_attribute,NF90_NOERR,NF90_GLOBAL,NF90_NOWRITE


!---------------------------Arguments-----------------------------------

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input  (nlfile=atm_in)
  logical,          intent(in) :: single_column_in
  real(r8),         intent(in) :: scmlat_in
  real(r8),         intent(in) :: scmlon_in

  ! Local variables
  character(len=*), parameter :: sub = 'scam_readnl'
  integer :: unitn, ierr, i
  integer  :: ncid
  integer  :: iatt
  integer  :: latidx, lonidx
  logical  :: adv
  real(r8) :: ioplat,ioplon

! this list should include any variable that you might want to include in the namelist
  namelist /scam_nl/ iopfile, scm_iop_lhflxshflxTg, scm_iop_Tg, scm_relaxation, &
       scm_relax_top_p,scm_relax_bot_p,scm_relax_tau_sec, &
       scm_cambfb_mode,scm_crm_mode,scm_zadv_uv,scm_zadv_T,scm_zadv_q,&
       scm_use_obs_T, scm_use_obs_uv, scm_use_obs_qv, &
       scm_relax_linear, scm_relax_tau_top_sec, &
       scm_relax_tau_bot_sec, scm_force_latlon, scm_relax_fincl, scm_backfill_iop_w_init

  single_column=single_column_in

  iopfile            = ' '
  scm_clubb_iop_name = ' '
  scm_relax_fincl(:) = ' '
  
  if( single_column ) then
     if( npes.gt.1) call endrun('SCAM_READNL: SCAM doesnt support using more than 1 pe.')

     if (.not. dycore_is('EUL') .or. plon /= 1 .or. plat /=1 ) then 
        call endrun('SCAM_SETOPTS: must compile model for SCAM mode when namelist parameter single_column is .true.')
     endif

     scmlat=scmlat_in
     scmlon=scmlon_in
     
     if( scmlat .lt. -90._r8 .or. scmlat .gt. 90._r8 ) then
        call endrun('SCAM_READNL: SCMLAT must be between -90. and 90. degrees.')
     elseif( scmlon .lt. 0._r8 .or. scmlon .gt. 360._r8 ) then
        call endrun('SCAM_READNL: SCMLON must be between 0. and 360. degrees.')
     end if
     
     ! Read namelist
     if (masterproc) then
        unitn = getunit()
        open( unitn, file=trim(nlfile), status='old' )
        call find_group_name(unitn, 'scam_nl', status=ierr)
        if (ierr == 0) then
           read(unitn, scam_nl, iostat=ierr)
           if (ierr /= 0) then
              call endrun(sub // ':: ERROR reading namelist')
           end if
        end if
        close(unitn)
        call freeunit(unitn)
     end if
     
     ! Error checking:
     
     iopfile = trim(iopfile)
     if( iopfile .ne. "" ) then 
        use_iop = .true.
     else
        call endrun('SCAM_READNL: must specify IOP file for single column mode')
     endif

     call wrap_open( iopfile, NF90_NOWRITE, ncid )

     if( nf90_inquire_attribute( ncid, NF90_GLOBAL, 'CAM_GENERATED_FORCING', iatt ) .EQ. NF90_NOERR ) then
        use_camiop = .true.
     else
        use_camiop = .false.
     endif
     
     ! If we are not forcing the lat and lon from the namelist use the closest lat and lon that is found in the IOP file.
     if (.not.scm_force_latlon) then
        call shr_scam_GetCloseLatLon( ncid, scmlat, scmlon, ioplat, ioplon, latidx, lonidx )
        write(iulog,*) 'SCAM_READNL: using closest IOP column to lat/lon specified in drv_in'
        write(iulog,*) '   requested lat,lon    =',scmlat,', ',scmlon
        write(iulog,*) '   closest IOP lat,lon  =',ioplat,', ',ioplon
     
        scmlat = ioplat
        scmlon = ioplon
     end if
     
     if (masterproc) then
        write (iulog,*) 'Single Column Model Options: '
        write (iulog,*) '============================='
        write (iulog,*) '  iopfile                     = ',trim(iopfile)
        write (iulog,*) '  scm_backfill_iop_w_init     = ',scm_backfill_iop_w_init
        write (iulog,*) '  scm_cambfb_mode             = ',scm_cambfb_mode
        write (iulog,*) '  scm_crm_mode                = ',scm_crm_mode
        write (iulog,*) '  scm_force_latlon            = ',scm_force_latlon
        write (iulog,*) '  scm_iop_Tg                  = ',scm_iop_Tg
        write (iulog,*) '  scm_iop_lhflxshflxTg        = ',scm_iop_lhflxshflxTg
        write (iulog,*) '  scm_relaxation              = ',scm_relaxation
        write (iulog,*) '  scm_relax_bot_p             = ',scm_relax_bot_p
        write (iulog,*) '  scm_relax_linear            = ',scm_relax_linear
        write (iulog,*) '  scm_relax_tau_bot_sec       = ',scm_relax_tau_bot_sec
        write (iulog,*) '  scm_relax_tau_sec           = ',scm_relax_tau_sec
        write (iulog,*) '  scm_relax_tau_top_sec       = ',scm_relax_tau_top_sec
        write (iulog,*) '  scm_relax_top_p             = ',scm_relax_top_p
        write (iulog,*) '  scm_use_obs_T               = ',scm_use_obs_T
        write (iulog,*) '  scm_use_obs_qv              = ',scm_use_obs_qv
        write (iulog,*) '  scm_use_obs_uv              = ',scm_use_obs_uv
        write (iulog,*) '  scm_zadv_T                  = ',trim(scm_zadv_T)
        write (iulog,*) '  scm_zadv_q                  = ',trim(scm_zadv_q)
        write (iulog,*) '  scm_zadv_uv                 = ',trim(scm_zadv_uv)
        write (iulog,*) '  scm_relax_finc: '
        ! output scm_relax_fincl character array
        do i=1,pcnst
           if (scm_relax_fincl(i) .ne. '') then
              adv = mod(i,4)==0
              if (adv) then
                 write (iulog, "(A18)") "'"//trim(scm_relax_fincl(i))//"',"
              else
                 write (iulog, "(A18)", ADVANCE="NO") "'"//trim(scm_relax_fincl(i))//"',"
              end if
           else
              exit
           end if
        end do
        print *
     end if
  end if
     
end subroutine scam_readnl

function scamvar_getdesc(name) result(vdesc)
  character(len=*), intent(in) :: name
  type(var_desc_t), pointer    :: vdesc
  character(len=errmsglen)     :: errmsg
  integer :: i
  
  nullify(vdesc)
  do i=1,scamvarcnt
     if(name .eq. scamvars(i)%name) then
        vdesc=>scamvars(i)%vdesc
        exit
     end if
  end do
  if(.not.associated(vdesc)) then
     errmsg = 'Could not find scam variable '//name
     call endrun(errmsg)
  end if
end function scamvar_getdesc

subroutine set_scam_var(name, index, v0, v1, v2, v3, l0, c0, c1 )
  use cam_abortutils,      only: endrun
  
  character(len=*), intent(in) :: name
  integer, intent(in) :: index
  real(r8), target, optional :: v0, v1(:), v2(:,:), v3(:,:,:)
  logical, target, optional  :: l0  
  character(len=scam_fieldname_len), target, optional :: c1(:)
  character(len=scam_fieldname_len), target, optional :: c0

  scamvars(index)%name=name
  scamvars(index)%type=pio_double
  if(present(v0)) then
     scamvars(index)%ndims = 0
     scamvars(index)%v0d => v0
  else if(present(v1)) then
     scamvars(index)%ndims = 1
     scamvars(index)%v1d => v1
  else if(present(v2)) then
     scamvars(index)%ndims = 2
     scamvars(index)%v2d => v2
  else if(present(v3)) then
     scamvars(index)%ndims = 3
     scamvars(index)%v3d => v3
  else if(present(l0)) then
     scamvars(index)%ndims = 0
     scamvars(index)%l0d => l0
     scamvars(index)%type=pio_int
  else if(present(c0)) then
     scamvars(index)%c0d => c0
     scamvars(index)%type=pio_char
     scamvars(index)%ndims=0
  else if(present(c1)) then
     scamvars(index)%c1d => c1
     scamvars(index)%type=pio_char
     scamvars(index)%ndims=1
  else
     call endrun('bad ndims in call to set_scam_var')
  end if
  allocate(scamvars(index)%vdesc)
  
end subroutine set_scam_var


subroutine scam_vars_setnames()

  ! Local variable
  integer :: rvindex


  rvindex = 1
  call set_scam_var('latiopbnd', rvindex, v1=latiopbnd)
  scamvars(rvindex)%dims(1) = latiopbnd_dim_ind
  
  rvindex = rvindex + 1
  call set_scam_var('loniopbnd', rvindex, v1=loniopbnd)
  scamvars(rvindex)%dims(1) = loniopbnd_dim_ind
  
  rvindex = rvindex + 1
  call set_scam_var('psinit', rvindex, v2=psinit)
  scamvars(rvindex)%dims(1) = hdimid1_dim_ind
  scamvars(rvindex)%dims(2) = hdimid2_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('divq', rvindex, v2=divq)
  scamvars(rvindex)%dims(1) = lev_dim_ind
  scamvars(rvindex)%dims(2) = pcnst_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('divq3d', rvindex, v2=divq3d)
  scamvars(rvindex)%dims(1) = lev_dim_ind
  scamvars(rvindex)%dims(2) = pcnst_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('divt', rvindex, v1=divt)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('divu', rvindex, v1=divu)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('divt3d', rvindex, v1=divt3d)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('divu3d', rvindex, v1=divu3d)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('have_divv', rvindex, l0=have_divv)

  rvindex = rvindex + 1
  call set_scam_var('divv', rvindex, v1=divv)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('divv3d', rvindex, v1=divv3d)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('have_aldif', rvindex, l0=have_aldif)

  rvindex = rvindex + 1
  call set_scam_var('have_aldir', rvindex, l0=have_aldir)

  rvindex = rvindex + 1
  call set_scam_var('have_asdif', rvindex, l0=have_asdif)

  rvindex = rvindex + 1
  call set_scam_var('have_asdir', rvindex, l0=have_asdir)

  rvindex = rvindex + 1
  call set_scam_var('have_cld', rvindex, l0=have_cld)

  rvindex = rvindex + 1
  call set_scam_var('have_cldice', rvindex, l0=have_cldice)

  rvindex = rvindex + 1
  call set_scam_var('have_cldliq', rvindex, l0=have_cldliq)

  rvindex = rvindex + 1
  call set_scam_var('have_clwp', rvindex, l0=have_clwp)

  rvindex = rvindex + 1
  call set_scam_var('have_divq', rvindex, l0=have_divq)

  rvindex = rvindex + 1
  call set_scam_var('have_divq3d', rvindex, l0=have_divq3d)

  rvindex = rvindex + 1
  call set_scam_var('have_divt', rvindex, l0=have_divt)

  rvindex = rvindex + 1
  call set_scam_var('have_divt3d', rvindex, l0=have_divt3d)

  rvindex = rvindex + 1
  call set_scam_var('have_divu', rvindex, l0=have_divu)

  rvindex = rvindex + 1
  call set_scam_var('have_divu3d', rvindex, l0=have_divu3d)

  rvindex = rvindex + 1
  call set_scam_var('have_divv3d', rvindex, l0=have_divv3d)

  rvindex = rvindex + 1
  call set_scam_var('have_numice', rvindex, l0=have_numice)

  rvindex = rvindex + 1
  call set_scam_var('have_numliq', rvindex, l0=have_numliq)

  rvindex = rvindex + 1
  call set_scam_var('have_omega', rvindex, l0=have_omega)

  rvindex = rvindex + 1
  call set_scam_var('have_phis', rvindex, l0=have_phis)

  rvindex = rvindex + 1
  call set_scam_var('have_prec', rvindex, l0=have_prec)

  rvindex = rvindex + 1
  call set_scam_var('have_ps', rvindex, l0=have_ps)

  rvindex = rvindex + 1
  call set_scam_var('have_ptend', rvindex, l0=have_ptend)

  rvindex = rvindex + 1
  call set_scam_var('have_q', rvindex, l0=have_q)

  rvindex = rvindex + 1
  call set_scam_var('have_q1', rvindex, l0=have_q1)

  rvindex = rvindex + 1
  call set_scam_var('have_q2', rvindex, l0=have_q2)

  rvindex = rvindex + 1
  call set_scam_var('have_t', rvindex, l0=have_t)

  rvindex = rvindex + 1
  call set_scam_var('have_u', rvindex, l0=have_u)

  rvindex = rvindex + 1
  call set_scam_var('have_v', rvindex, l0=have_v)

  rvindex = rvindex + 1
  call set_scam_var('have_vertdivq', rvindex, l0=have_vertdivq)

  rvindex = rvindex + 1
  call set_scam_var('have_vertdivt', rvindex, l0=have_vertdivt)

  rvindex = rvindex + 1
  call set_scam_var('have_vertdivu', rvindex, l0=have_vertdivu)

  rvindex = rvindex + 1
  call set_scam_var('have_vertdivv', rvindex, l0=have_vertdivv)

  rvindex = rvindex + 1
  call set_scam_var('qdiff', rvindex, v1=qdiff)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('qobs', rvindex, v1=qobs)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('scm_relax_bot_p', rvindex, v0=scm_relax_bot_p)

  rvindex = rvindex + 1
  call set_scam_var('scm_relax_linear', rvindex, l0=scm_relax_linear)

  rvindex = rvindex + 1
  call set_scam_var('scm_relax_tau_bot_sec', rvindex, v0=scm_relax_tau_bot_sec)

  rvindex = rvindex + 1
  call set_scam_var('scm_relax_tau_sec', rvindex, v0=scm_relax_tau_sec)

  rvindex = rvindex + 1
  call set_scam_var('scm_relax_tau_top_sec', rvindex, v0=scm_relax_tau_top_sec)

  rvindex = rvindex + 1
  call set_scam_var('scm_relax_top_p', rvindex, v0=scm_relax_top_p)

  rvindex = rvindex + 1
  call set_scam_var('scm_relaxation', rvindex, l0=scm_relaxation)

  rvindex = rvindex + 1
  call set_scam_var('scm_use_obs_qv', rvindex, l0=scm_use_obs_qv)

  rvindex = rvindex + 1
  call set_scam_var('scm_use_obs_t', rvindex, l0=scm_use_obs_t)

  rvindex = rvindex + 1
  call set_scam_var('scm_use_obs_uv', rvindex, l0=scm_use_obs_uv)

  rvindex = rvindex + 1
  call set_scam_var('scm_zadv_q', rvindex, c0=scm_zadv_q)
  scamvars(rvindex)%dims(1) = scam_fieldname_len_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('scm_zadv_t', rvindex, c0=scm_zadv_t)
  scamvars(rvindex)%dims(1) = scam_fieldname_len_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('scm_zadv_uv', rvindex, c0=scm_zadv_uv)
  scamvars(rvindex)%dims(1) = scam_fieldname_len_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('tdiff', rvindex, v1=tdiff)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('tobs', rvindex, v1=tobs)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('uobs', rvindex, v1=uobs)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('use_3dfrc', rvindex, l0=use_3dfrc)

  rvindex = rvindex + 1
  call set_scam_var('use_camiop', rvindex, l0=use_camiop)

  rvindex = rvindex + 1
  call set_scam_var('vertdivq', rvindex, v2=vertdivq)
  scamvars(rvindex)%dims(1) = lev_dim_ind
  scamvars(rvindex)%dims(2) = pcnst_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('vertdivt', rvindex, v1=vertdivt)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('vertdivu', rvindex, v1=vertdivu)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('vertdivv', rvindex, v1=vertdivv)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('vobs', rvindex, v1=vobs)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('wfld', rvindex, v1=wfld)
  scamvars(rvindex)%dims(1) = lev_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('qinitobs', rvindex, v2=qinitobs)
  scamvars(rvindex)%dims(1) = lev_dim_ind
  scamvars(rvindex)%dims(2) = pcnst_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('scm_relax_fincl', rvindex, c1=scm_relax_fincl)
  scamvars(rvindex)%dims(1) = scam_fieldname_len_dim_ind
  scamvars(rvindex)%dims(2) = pcnst_dim_ind

  rvindex = rvindex + 1
  call set_scam_var('wfldh', rvindex, v1=wfldh)
  scamvars(rvindex)%dims(1) = ilev_dim_ind

end subroutine scam_vars_setnames

subroutine scam_dims_setnames()

  scamdims(latiopbnd_dim_ind)%name = 'latiopbnd'
  scamdims(latiopbnd_dim_ind)%len  = 2
  scamdims(latiopbnd_dim_ind)%define  = .true.

  scamdims(loniopbnd_dim_ind)%name = 'loniopbnd'
  scamdims(loniopbnd_dim_ind)%len  = 2
  scamdims(loniopbnd_dim_ind)%define  = .true.

  scamdims(hdimid1_dim_ind)%name = 'hdimid1'
  scamdims(hdimid1_dim_ind)%define  = .false.

  scamdims(hdimid2_dim_ind)%name = 'hdimid2'
  scamdims(hdimid2_dim_ind)%define  = .false.
  
  scamdims(lev_dim_ind)%name = 'lev'
  scamdims(lev_dim_ind)%define  = .false.
  
  scamdims(ilev_dim_ind)%name = 'ilev'
  scamdims(ilev_dim_ind)%define  = .false.
  
  scamdims(pcnst_dim_ind)%name = 'pcnst'
  scamdims(pcnst_dim_ind)%define  = .false.

  scamdims(scam_fieldname_len_dim_ind)%name    = 'scam_fieldname_len'
  scamdims(scam_fieldname_len_dim_ind)%len  = scam_fieldname_len
  scamdims(scam_fieldname_len_dim_ind)%define  = .true.
  
end subroutine scam_dims_setnames

subroutine scam_restart_init (File,hdimids,vdimids)
  
  !---------------------------------------------------------------------------
  !
  ! Arguments
  !
  type(file_desc_t), intent(inout) :: File                 ! Pio file Handle
  integer, optional, intent(in)    :: hdimids(:),vdimids(:)! lat/lon/lev dimids
  !
  ! Local
  !
  integer :: dimids(4), ndims,mcdimid,pcnstdimid,latdimid,londimid,levdimid,ilevdimid
  integer :: ierr, i, k, err_handling

  call scam_dims_setnames()
  call scam_vars_setnames()

!jt  call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
  call pio_seterrorhandling(File, PIO_INTERNAL_ERROR, err_handling)
  if (present(hdimids)) then
     londimid=hdimids(1)
     latdimid=hdimids(2)
     ilevdimid=vdimids(1)
     levdimid=vdimids(2)
  else
     ierr = pio_inq_dimid(File,'lon', londimid)
     ierr = pio_inq_dimid(File,'lat', latdimid)
     ierr = pio_inq_dimid(File,'lev', levdimid)
     ierr = pio_inq_dimid(File,'ilev', ilevdimid)
  end if

  do i=1,scamdimcnt
     if (scamdims(i)%define) &
          ierr = pio_def_dim(File,  scamdims(i)%name, scamdims(i)%len, scamdims(i)%dimid)
  end do
  ! Set the dimensions already contained in the file
  scamdims(hdimid1_dim_ind)%dimid=londimid
  scamdims(hdimid2_dim_ind)%dimid=latdimid
  scamdims(lev_dim_ind)%dimid=levdimid
  scamdims(ilev_dim_ind)%dimid=ilevdimid

  ierr = pio_inq_dimid(File,'pcnst', pcnstdimid)
  scamdims(pcnst_dim_ind)%dimid=pcnstdimid

  do i = 1, scamvarcnt
     ndims = scamvars(i)%ndims
     dimids(:)=0
     do k = 1 ,ndims
        dimids(k) = scamdims(scamvars(i)%dims(k))%dimid
     end do
     if (scamvars(i)%ndims==0) then
        ierr = pio_def_var(File, scamvars(i)%name, scamvars(i)%type, scamvars(i)%vdesc)
     else
        ierr = pio_def_var(File, scamvars(i)%name, scamvars(i)%type, dimids(1:scamvars(i)%ndims), scamvars(i)%vdesc)
     end if
  end do

  call pio_seterrorhandling(File, err_handling)

end subroutine scam_restart_init

subroutine scam_restart_write(File)
!---------------------------Arguments-----------------------------------

  type(file_desc_t), intent(inout) :: File         ! PIO restart file pointer
  !
  ! Local workspace
  !
  integer :: ierr, i, err_handling, itmp
  type(var_desc_t), pointer :: vdesc
  call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
  
  !-----------------------------------------------------------------------
  ! Write the scam restart data
  !-----------------------------------------------------------------------
  
  nullify(vdesc)
  do i=1,scamvarcnt
     vdesc=>scamvars(i)%vdesc
     select case (scamvars(i)%type)
     case(pio_double)
        select case (scamvars(i)%ndims)
        case(0)
           ierr= pio_put_var(File, vdesc, scamvars(i)%v0d)
        case(1)
           ierr= pio_put_var(File, vdesc, scamvars(i)%v1d)
        case(2)
           ierr= pio_put_var(File, vdesc, scamvars(i)%v2d)
        case(3)
           ierr= pio_put_var(File, vdesc, scamvars(i)%v3d)
        case default
           errmsg='scam_restart_write: Error writing restart variable: '//scamvars(i)%name
           call endrun(trim(errmsg))
        end select
     case(pio_int)
           select case (scamvars(i)%ndims)
              case(0)
                 itmp=scamvars(i)%l0d
                 ierr= pio_put_var(File, vdesc, itmp)
              case default
                 errmsg='scam_restart_write: Error writing restart variable: '//scamvars(i)%name
                 call endrun(trim(errmsg))
              end select
     case(pio_char)
        if (scamvars(i)%ndims==0) then
           ierr= pio_put_var(File, vdesc, scamvars(i)%c0d)
        else
           ierr= pio_put_var(File, vdesc, scamvars(i)%c1d)
        end if
     case default
        errmsg='scam_restart_write: Error writing restart variable: '//scamvars(i)%name
        call endrun(trim(errmsg))
     end select
  end do

  call pio_seterrorhandling(File, err_handling)
  
  return
end subroutine scam_restart_write
subroutine scam_restart_read(File)
!---------------------------Arguments-----------------------------------

  type(file_desc_t), intent(inout) :: File         ! PIO restart file pointer
  !
  ! Local workspace
  !
  integer :: ierr, i, err_handling, itmp
  type(var_desc_t), pointer :: vdesc
  integer :: latiopbnd_dimid
  integer :: loniopbnd_dimid
  integer :: lat_dimid
  integer :: lon_dimid
  
  integer                          :: nlatiopbnd
  integer                          :: nloniopbnd
  integer                          :: nlats,nlons
  
  call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)

  ierr = pio_inq_dimid(File, 'latiopbnd', latiopbnd_dimid)
  ierr = pio_inq_dimid(File, 'loniopbnd', loniopbnd_dimid)
  ierr = pio_inq_dimid(File, 'lat', lat_dimid)
  ierr = pio_inq_dimid(File, 'lon', lon_dimid)
  
  ierr = pio_inq_dimlen(File, latiopbnd_dimid, nlatiopbnd)
  ierr = pio_inq_dimlen(File, loniopbnd_dimid, nloniopbnd)
  ierr = pio_inq_dimlen(File, lat_dimid, nlats)
  ierr = pio_inq_dimlen(File, lon_dimid, nlons)
  
  allocate(latiopbnd(nlatiopbnd))
  allocate(loniopbnd(nloniopbnd))
  allocate(psinit(nlons,nlats))

  call scam_dims_setnames()
  call scam_vars_setnames()

  
  !-----------------------------------------------------------------------
  ! Read the scam restart data
  !-----------------------------------------------------------------------
  
!jt  nullify(vdesc)
  do i=1,scamvarcnt
!jt     vdesc=>scamvars(i)%vdesc
     ierr = PIO_Inq_varid(File, scamvars(i)%name, scamvars(i)%vdesc)
     select case (scamvars(i)%type)
     case(pio_double)
        select case (scamvars(i)%ndims)
        case(0)
           ierr= pio_get_var(File, scamvars(i)%vdesc, scamvars(i)%v0d)
        case(1)
           ierr= pio_get_var(File, scamvars(i)%vdesc, scamvars(i)%v1d)
        case(2)
           ierr= pio_get_var(File, scamvars(i)%vdesc, scamvars(i)%v2d)
        case(3)
           ierr= pio_get_var(File, scamvars(i)%vdesc, scamvars(i)%v3d)
        case default
           errmsg='scam_restart_read: Error writing restart variable: '//scamvars(i)%name
           call endrun(trim(errmsg))
        end select
     case(pio_int)
           select case (scamvars(i)%ndims)
              case(0)
                 ierr= pio_get_var(File, scamvars(i)%vdesc, itmp)
                 scamvars(i)%l0d=itmp
              case default
                 errmsg='scam_restart_read: Error writing restart variable: '//scamvars(i)%name
                 call endrun(trim(errmsg))
              end select
     case(pio_char)
           select case (scamvars(i)%ndims)
              case(0)
                 ierr= pio_get_var(File, scamvars(i)%vdesc, scamvars(i)%c0d)
              case(1)
                 ierr= pio_get_var(File, scamvars(i)%vdesc, scamvars(i)%c1d)
              case default
                 errmsg='scam_restart_read: Error reading restart variable: '//scamvars(i)%name
                 call endrun(trim(errmsg))
              end select
     case default
        errmsg='scam_restart_read: Error reading restart variable: '//scamvars(i)%name
        call endrun(trim(errmsg))
     end select
  end do

  call pio_seterrorhandling(File, err_handling)
  
  return
end subroutine scam_restart_read

subroutine scam_restart_read1(File)
  !-----------------------------------------------------------------------
  !
  ! Arguments
  !
  type(file_desc_t), intent(inout) :: File            ! unit number
  !
  ! Local workspace
  !
  integer :: ierr, i
  type(var_desc_t), pointer :: vdesc
  integer :: latiopbnd_dimid
  integer :: loniopbnd_dimid
  integer :: lat_dimid
  integer :: lon_dimid
  
  integer                          :: nlatiopbnd
  integer                          :: nloniopbnd
  integer                          :: nlats,nlons
  integer                          :: err_handling

  call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
  
  ierr = pio_inq_dimid(File, 'latiopbnd', latiopbnd_dimid)
  ierr = pio_inq_dimid(File, 'loniopbnd', loniopbnd_dimid)
  ierr = pio_inq_dimid(File, 'lat', lat_dimid)
  ierr = pio_inq_dimid(File, 'lon', lon_dimid)
  
  ierr = pio_inq_dimlen(File, latiopbnd_dimid, nlatiopbnd)
  ierr = pio_inq_dimlen(File, loniopbnd_dimid, nloniopbnd)
  ierr = pio_inq_dimlen(File, lat_dimid, nlats)
  ierr = pio_inq_dimlen(File, lon_dimid, nlons)
  
  allocate(latiopbnd(nlatiopbnd))
  allocate(loniopbnd(nloniopbnd))
  allocate(psinit(nlons,nlats))
   
  !-----------------------------------------------------------------------
  ! Read the scam restart data
  !-----------------------------------------------------------------------
  
  do i=1,scamvarcnt
     nullify(vdesc)
     vdesc=>scamvars(i)%vdesc
!     ndims=scamvars(i)%ndims
!     ierr = pio_get_var(File, vdesc, scamvars(i)%data(ndims))
  end do
  
  call pio_seterrorhandling(File, err_handling)
  return
end subroutine scam_restart_read1
end module scamMod

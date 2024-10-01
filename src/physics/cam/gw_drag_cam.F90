module gw_drag_cam

!--------------------------------------------------------------------------
! CAM and WACCM gravity wave parameterizations were merged by Sean Patrick
! Santos in Summer 2013, and at the same time, gw_drag was split into
! various modules. This is the CAM interface and driver module. The below
! notes are for the old CAM and WACCM versions of gw_drag.
!--------------------------------------------------------------------------
! This file came from wa17 and was modified by Fabrizio: 07-02-2004
! Standard gw_drag with modification (6) of latitude profile of gw spectrum
!--------------------------------------------------------------------------
! Purpose:
!
! Module to compute the forcing due to parameterized gravity waves. Both an
! orographic and an internal source spectrum are considered.
!
! Author: Byron Boville
!
!--------------------------------------------------------------------------
  use shr_kind_mod,   only: r8=>shr_kind_r8, cl=>shr_kind_cl
  use shr_log_mod,    only: errMsg => shr_log_errMsg
  use shr_assert_mod, only: shr_assert

  use ppgrid,         only: pcols, pver, begchunk, endchunk
  use constituents,   only: pcnst
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use spmd_utils,     only: masterproc
  use cam_history,    only: outfld
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use error_messages, only: alloc_err


  use ref_pres,       only: do_molec_diff, nbot_molec
  use physconst,      only: cpair

  ! These are the actual switches for different gravity wave sources.
  use phys_control,   only: use_gw_oro, use_gw_front, use_gw_front_igw, &
                            use_gw_convect_dp, use_gw_convect_sh,       &
                            use_simple_phys, use_gw_movmtn_pbl

  use gw_drag,        only: gw_drag_init
  use gw_common,      only: GWBand
  use gw_convect,     only: BeresSourceDesc
  use gw_movmtn,      only: MovMtnSourceDesc
  use gw_front,       only: CMSourceDesc

! Typical module header
  implicit none
  private
  save

!
! PUBLIC: interfaces
!
  public :: gw_drag_cam_readnl           ! Read namelist
  public :: gw_drag_cam_init             ! Initialization
  public :: gw_drag_cam_tend             ! interface to actual parameterization

!
! PRIVATE: Rest of the data and interfaces are private to this module
!
  real(r8), parameter :: unset_r8 = huge(1._r8)

  ! A mid-scale "band" with only stationary waves (l = 0).
  type(GWBand) :: band_oro
  ! Medium scale waves.
  type(GWBand) :: band_mid
  ! Long scale waves for IGWs.
  type(GWBand) :: band_long
  ! Medium scale waves for moving mountain
  type(GWBand) :: band_movmtn

  ! Top level for gravity waves.
  integer, parameter :: ktop = 1
  ! Bottom level for frontal waves.
  integer :: kbot_front

  ! Factor for SH orographic waves.
  real(r8) :: gw_oro_south_fac = 1._r8

  ! Frontogenesis function critical threshold.
  real(r8) :: frontgfc = unset_r8

  ! Tendency efficiencies.

  ! Ridge scheme.
  logical  :: use_gw_rdg_beta      = .false.
  integer  :: n_rdg_beta           = 1
  real(r8) :: effgw_rdg_beta       = unset_r8
  real(r8) :: effgw_rdg_beta_max   = unset_r8
  real(r8) :: rdg_beta_cd_llb      = unset_r8  ! Low-level obstacle drag coefficient Ridge scheme.
  logical  :: trpd_leewv_rdg_beta  = .false.

  logical  :: use_gw_rdg_gamma     = .false.
  integer  :: n_rdg_gamma          = -1
  real(r8) :: effgw_rdg_gamma      = unset_r8
  real(r8) :: effgw_rdg_gamma_max  = unset_r8
  real(r8) :: rdg_gamma_cd_llb     = unset_r8
  logical  :: trpd_leewv_rdg_gamma = .false.
  character(len=cl) :: bnd_rdggm   = 'bnd_rdggm' ! full pathname for meso-Gamma ridge dataset

  ! Orography.
  real(r8) :: effgw_oro = unset_r8
  ! C&M scheme.
  real(r8) :: effgw_cm = unset_r8
  ! C&M scheme (inertial waves).
  real(r8) :: effgw_cm_igw = unset_r8
  ! Beres (deep convection).
  real(r8) :: effgw_beres_dp = unset_r8
  ! Beres (shallow convection).
  real(r8) :: effgw_beres_sh = unset_r8

  ! Horzontal wavelengths [m].
  real(r8), parameter :: wavelength_mid = 1.e5_r8
  real(r8), parameter :: wavelength_long = 3.e5_r8

  ! Background stress source strengths.
  real(r8) :: taubgnd = unset_r8
  real(r8) :: taubgnd_igw = unset_r8

  ! Whether or not to use a polar taper for frontally generated waves.
  logical :: gw_polar_taper = .false.

  ! Whether or not to enforce an upper boundary condition of tau = 0.
  ! (Like many variables, this is only here to hold the value between
  ! the readnl phase and the init phase of the CAM physics; only gw_common
  ! should actually use it.)
  logical :: tau_0_ubc = .false.

  ! Whether or not to limit tau *before* applying any efficiency factors.
  logical :: gw_limit_tau_without_eff = .false.

  ! Whether or not to apply tendency max
  logical :: gw_apply_tndmax = .true.

  ! Files to read Beres source spectra from.
  character(len=cl) :: gw_drag_file = ""
  character(len=cl) :: gw_drag_file_sh = ""
  character(len=cl) :: gw_drag_file_mm = ""

  ! Indices into pbuf
  integer :: kvt_idx      = -1
  integer :: ttend_dp_idx = -1
  integer :: ttend_sh_idx = -1
  integer :: frontgf_idx  = -1
  integer :: frontga_idx  = -1
  integer :: sgh_idx      = -1

  ! From CLUBB
  integer :: ttend_clubb_idx  = -1
  integer :: upwp_clubb_gw_idx   = -1
  integer :: vpwp_clubb_gw_idx   = -1
  integer :: thlp2_clubb_gw_idx  = -1
  integer :: wpthlp_clubb_gw_idx  = -1

  ! anisotropic ridge fields
  integer, parameter :: prdg = 16

     ! Meso Gamma
  real(r8), allocatable, dimension(:,:), target :: &
     rdg_gbxar
     ! Meso Beta
  real(r8), allocatable, dimension(:,:,:), target :: &
     rdg_hwdth,  &
     rdg_clngt,  &
     rdg_mxdis,  &
     rdg_anixy,  &
     rdg_angll

     ! Meso Gamma
  real(r8), allocatable, dimension(:,:), target :: &
     rdg_gbxarg

  real(r8), allocatable, dimension(:,:,:), target :: &
     rdg_hwdthg, &
     rdg_clngtg, &
     rdg_mxdisg, &
     rdg_anixyg, &
     rdg_angllg

  ! Water constituent indices for budget
  integer :: ixcldliq = -1
  integer :: ixcldice = -1

  ! Prefixes for history field names
  character(len=1), parameter :: cm_pf = " "
  character(len=1), parameter :: cm_igw_pf = "I"
  character(len=1), parameter :: beres_dp_pf = "B"
  character(len=1), parameter :: beres_sh_pf = "S"

  ! namelist
  logical          :: history_amwg                   ! output the variables used by the AMWG diag package
  logical  :: gw_lndscl_sgh = .true. ! scale SGH by land frac
  real(r8) :: gw_prndl = 0.25_r8
  real(r8) :: gw_qbo_hdepth_scaling = 1._r8 ! heating depth scaling factor

  ! Width of gaussian used to create frontogenesis tau profile [m s-1].
  real(r8) :: front_gaussian_width = -huge(1._r8)

  real(r8) :: alpha_gw_movmtn
  real(r8) :: gw_bot_taper_pres

  logical :: gw_top_taper=.false.

  ! Maximum wave number and width of spectrum bins.
  integer :: pgwv = -1
  real(r8) :: gw_dc = unset_r8
  integer :: pgwv_long = -1
  real(r8) :: gw_dc_long = unset_r8

  ! fcrit2 for the mid-scale waves has been made a namelist variable to
  ! facilitate backwards compatibility with the CAM3 version of this
  ! parameterization.  In CAM3, fcrit2=0.5.
  real(r8) :: fcrit2 = unset_r8   ! critical froude number squared


!==========================================================================
contains
!==========================================================================

subroutine gw_drag_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8, &
                             mpi_character, mpi_logical, mpi_integer
  use gw_rdg_cam,      only: gw_rdg_cam_readnl

  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'gw_drag_cam_readnl'

  namelist /gw_drag_nl/ pgwv, gw_dc, pgwv_long, gw_dc_long, tau_0_ubc, &
       effgw_beres_dp, effgw_beres_sh, effgw_cm, effgw_cm_igw, effgw_oro, &
       fcrit2, frontgfc, gw_drag_file, gw_drag_file_sh, gw_drag_file_mm, taubgnd, &
       taubgnd_igw, gw_polar_taper, &
       use_gw_rdg_beta, n_rdg_beta, effgw_rdg_beta, effgw_rdg_beta_max, &
       rdg_beta_cd_llb, trpd_leewv_rdg_beta, &
       use_gw_rdg_gamma, n_rdg_gamma, effgw_rdg_gamma, effgw_rdg_gamma_max, &
       rdg_gamma_cd_llb, trpd_leewv_rdg_gamma, bnd_rdggm, &
       gw_oro_south_fac, gw_limit_tau_without_eff, &
       gw_lndscl_sgh, gw_prndl, gw_apply_tndmax, gw_qbo_hdepth_scaling, &
       gw_top_taper, front_gaussian_width, alpha_gw_movmtn, gw_bot_taper_pres

!----------------------------------------------------------------------

  if (use_simple_phys) return

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'gw_drag_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, gw_drag_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(sub // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)
  end if

  call mpi_bcast(pgwv, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: pgwv")
  call mpi_bcast(gw_dc, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_dc")
  call mpi_bcast(pgwv_long, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: pgwv_long")
  call mpi_bcast(gw_dc_long, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_dc_long")
  call mpi_bcast(tau_0_ubc, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: tau_0_ubc")
  call mpi_bcast(effgw_beres_dp, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_beres_dp")
  call mpi_bcast(effgw_beres_sh, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_beres_sh")
  call mpi_bcast(effgw_cm, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_cm")
  call mpi_bcast(effgw_cm_igw, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_cm_igw")
  call mpi_bcast(effgw_oro, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_oro")

  call mpi_bcast(use_gw_rdg_beta, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: use_gw_rdg_beta")
  call mpi_bcast(n_rdg_beta, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: n_rdg_beta")
  call mpi_bcast(effgw_rdg_beta, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_rdg_beta")
  call mpi_bcast(effgw_rdg_beta_max, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_rdg_beta_max")
  call mpi_bcast(rdg_beta_cd_llb, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rdg_beta_cd_llb")
  call mpi_bcast(trpd_leewv_rdg_beta, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: trpd_leewv_rdg_beta")

  call mpi_bcast(use_gw_rdg_gamma, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: use_gw_rdg_gamma")
  call mpi_bcast(n_rdg_gamma, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: n_rdg_gamma")
  call mpi_bcast(effgw_rdg_gamma, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_rdg_gamma")
  call mpi_bcast(effgw_rdg_gamma_max, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_rdg_gamma_max")
  call mpi_bcast(rdg_gamma_cd_llb, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rdg_gamma_cd_llb")
  call mpi_bcast(trpd_leewv_rdg_gamma, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: trpd_leewv_rdg_gamma")
  call mpi_bcast(bnd_rdggm, len(bnd_rdggm), mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: bnd_rdggm")

  call mpi_bcast(gw_oro_south_fac, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_oro_south_fac")
  call mpi_bcast(fcrit2, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fcrit2")
  call mpi_bcast(frontgfc, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frontgfc")
  call mpi_bcast(taubgnd, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: taubgnd")
  call mpi_bcast(taubgnd_igw, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: taubgnd_igw")

  call mpi_bcast(gw_polar_taper, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_polar_taper")
  call mpi_bcast(gw_limit_tau_without_eff, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_limit_tau_without_eff")
  call mpi_bcast(gw_apply_tndmax, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_apply_tndmax")

  call mpi_bcast(gw_top_taper, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_top_taper")

  call mpi_bcast(gw_bot_taper_pres, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_bot_taper_pres")

  call mpi_bcast(gw_lndscl_sgh, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_lndscl_sgh")
  call mpi_bcast(gw_prndl, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_prndl")
  call mpi_bcast(gw_qbo_hdepth_scaling, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_qbo_hdepth_scaling")

  call mpi_bcast(gw_drag_file, len(gw_drag_file), mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_drag_file")
  call mpi_bcast(gw_drag_file_sh, len(gw_drag_file_sh), mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_drag_file_sh")
  call mpi_bcast(gw_drag_file_mm, len(gw_drag_file_mm), mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_drag_file_mm")

  call mpi_bcast(front_gaussian_width, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: front_gaussian_width")

  call mpi_bcast(alpha_gw_movmtn, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: alpha_gw_movmtn")

  ! Check if fcrit2 was set.
  call shr_assert(fcrit2 /= unset_r8, &
       "gw_drag_cam_readnl: fcrit2 must be set via the namelist."// &
       errMsg(__FILE__, __LINE__))

  ! Check if pgwv was set.
  call shr_assert(pgwv >= 0, &
       "gw_drag_cam_readnl: pgwv must be set via the namelist and &
       &non-negative."// &
       errMsg(__FILE__, __LINE__))

  ! Check if gw_dc was set.
  call shr_assert(gw_dc /= unset_r8, &
       "gw_drag_cam_readnl: gw_dc must be set via the namelist."// &
       errMsg(__FILE__, __LINE__))

  band_oro = GWBand(0, gw_dc, fcrit2, wavelength_mid)
  band_mid = GWBand(pgwv, gw_dc, 1.0_r8, wavelength_mid)
  band_long = GWBand(pgwv_long, gw_dc_long, 1.0_r8, wavelength_long)
  band_movmtn = GWBand(0, gw_dc, 1.0_r8, wavelength_mid)

  if (use_gw_rdg_gamma .or. use_gw_rdg_beta) then
     call gw_rdg_cam_readnl(nlfile)
  end if

end subroutine gw_drag_cam_readnl

!==========================================================================

subroutine gw_drag_cam_init()

  !-----------------------------------------------------------------------
  ! Time independent initialization for multiple gravity wave
  ! parameterization.
  !-----------------------------------------------------------------------

  use cam_history,      only: addfld, add_default, horiz_only
  use cam_history,      only: register_vector_field
  use air_composition,  only: cpairv
  use interpolate_data, only: lininterp
  use phys_control,     only: phys_getopts
  use physics_buffer,   only: pbuf_get_index
  use constituents,     only: cnst_get_ind

  use cam_initfiles,    only: topo_file_get_id

  ! temporary for restart with ridge scheme
  use cam_initfiles,    only: bnd_topo

  use cam_pio_utils,    only: cam_pio_openfile
  use cam_grid_support, only: cam_grid_check, cam_grid_id
  use cam_grid_support, only: cam_grid_get_dim_names
  use pio,              only: file_desc_t, pio_nowrite, pio_closefile
  use ncdio_atm,        only: infld
  use ioFileMod,        only: getfil

  use ref_pres,   only: pref_edge, pref_mid
  use physconst,  only: gravit, rair, rearth, pi

  !---------------------------Local storage-------------------------------

  integer          :: i, l, k
  character(len=1) :: cn

  ! Index for levels at specific pressures.
  integer :: kfront

  ! output tendencies and state variables for CAM4 temperature,
  ! water vapor, cloud ice and cloud liquid budgets.
  logical :: history_budget
  ! output history file number for budget fields
  integer :: history_budget_histfile_num
  ! output variables of interest in WACCM runs
  logical :: history_waccm

  ! Read data from file
  type(file_desc_t), pointer :: fh_topo
  type(file_desc_t)  :: fh_rdggm
  integer            :: grid_id
  character(len=8)   :: dim1name, dim2name
  logical            :: found
  character(len=cl) :: bnd_rdggm_loc   ! filepath of topo file on local disk

  ! Allow reporting of error messages.
  character(len=128) :: errstring
  character(len=*), parameter :: sub = 'gw_init'

  ! temporary workaround for restart w/ ridge scheme
  character(len=cl) :: bnd_topo_loc   ! filepath of topo file on local disk

  ! Beres settings and table.
  type(BeresSourceDesc) :: beres_dp_desc
  type(BeresSourceDesc) :: beres_sh_desc

  ! Moving mountain settings and table.
  type(MovMtnSourceDesc) :: movmtn_desc

  ! Frontogenesis wave settings.
  type(CMSourceDesc) :: cm_desc
  type(CMSourceDesc) :: cm_igw_desc
  character(len=512)              :: errormsg
  integer                         :: errorflg


  !-----------------------------------------------------------------------

  if (do_molec_diff) then
     kvt_idx     = pbuf_get_index('kvt')
  end if

  if (masterproc) then
     write(iulog,*) ' '
     write(iulog,*) "GW_DRAG: band_mid%ngwv = ", band_mid%ngwv
     do l = -band_mid%ngwv, band_mid%ngwv
        write (iulog,'(A,I0,A,F7.2)') &
             "GW_DRAG: band_mid%cref(",l,") = ",band_mid%cref(l)
     enddo
     write(iulog,*) 'GW_DRAG: band_mid%kwv = ', band_mid%kwv
     write(iulog,*) 'GW_DRAG: band_mid%fcrit2 = ', band_mid%fcrit2
     write(iulog,*) ' '
     write(iulog,*) "GW_DRAG: band_long%ngwv = ", band_long%ngwv
     do l = -band_long%ngwv, band_long%ngwv
        write (iulog,'(A,I2,A,F7.2)') &
             "GW_DRAG: band_long%cref(",l,") = ",band_long%cref(l)
     enddo
     write(iulog,*) 'GW_DRAG: band_long%kwv = ', band_long%kwv
     write(iulog,*) 'GW_DRAG: band_long%fcrit2 = ', band_long%fcrit2
     write(iulog,*) ' '
  end if

  ! Used to decide whether temperature tendencies should be output.
  call phys_getopts( history_budget_out = history_budget, &
       history_budget_histfile_num_out = history_budget_histfile_num, &
       history_waccm_out = history_waccm, &
       history_amwg_out   = history_amwg  )

  if ( use_gw_oro ) then

     if (effgw_oro == unset_r8) then
        call endrun("gw_init: Orographic gravity waves enabled, &
             &but effgw_oro was not set.")
     end if
  end if

  if (use_gw_oro .or. use_gw_rdg_beta .or. use_gw_rdg_gamma) then

     sgh_idx = pbuf_get_index('SGH')

     ! Declare history variables for orographic term
     call addfld ('TAUAORO',    (/ 'ilev' /), 'I','N m-2',  &
          'Total stress from original OGW scheme')
     call addfld ('TTGWORO',    (/ 'lev' /), 'A','K s-1',  &
          'T tendency - orographic gravity wave drag')
     call addfld ('TTGWSDFORO', (/ 'lev' /), 'A','K s-1',  &
          'T tendency - orographic gravity wave, diffusion.')
     call addfld ('TTGWSKEORO', (/ 'lev' /), 'A','K s-1',  &
          'T tendency - orographic gravity wave, breaking KE.')
     call addfld ('UTGWORO',    (/ 'lev' /), 'A','m s-2', &
          'U tendency - orographic gravity wave drag')
     call addfld ('VTGWORO',    (/ 'lev' /), 'A','m s-2', &
          'V tendency - orographic gravity wave drag')
     call register_vector_field('UTGWORO', 'VTGWORO')
     call addfld ('TAUGWX',     horiz_only,  'A','N m-2', &
          'Zonal gravity wave surface stress')
     call addfld ('TAUGWY',     horiz_only,  'A','N m-2', &
          'Meridional gravity wave surface stress')
     call register_vector_field('TAUGWX', 'TAUGWY')

     if (history_amwg) then
        call add_default('TAUGWX  ', 1, ' ')
        call add_default('TAUGWY  ', 1, ' ')
     end if

     if (history_waccm) then
        call add_default('UTGWORO ', 1, ' ')
        call add_default('VTGWORO ', 1, ' ')
        call add_default('TAUGWX  ', 1, ' ')
        call add_default('TAUGWY  ', 1, ' ')
     end if

  end if

  if (use_gw_rdg_beta .or. use_gw_rdg_gamma) then
     grid_id = cam_grid_id('physgrid')
     if (.not. cam_grid_check(grid_id)) then
        call endrun(sub//': ERROR: no "physgrid" grid')
     end if
     call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
  end if

  if (use_gw_rdg_beta) then

     if (effgw_rdg_beta == unset_r8) then
        call endrun(sub//": ERROR: Anisotropic OGW enabled, &
                     &but effgw_rdg_beta was not set.")
     end if

     fh_topo => topo_file_get_id()
     bnd_topo_loc = ' '
     if (.not. associated(fh_topo)) then

        ! Try to open topo file here.  This workaround will not be needed
        ! once the refactored initialization sequence is on trunk.

        allocate(fh_topo)
        ! Error exit is from getfil if file not found.
        call getfil(bnd_topo, bnd_topo_loc)
        call cam_pio_openfile(fh_topo, bnd_topo_loc, PIO_NOWRITE)

     end if

     ! Get beta ridge data
     allocate( &
        rdg_gbxar(pcols,begchunk:endchunk),      &
        rdg_hwdth(pcols,prdg,begchunk:endchunk), &
        rdg_clngt(pcols,prdg,begchunk:endchunk), &
        rdg_mxdis(pcols,prdg,begchunk:endchunk), &
        rdg_anixy(pcols,prdg,begchunk:endchunk), &
        rdg_angll(pcols,prdg,begchunk:endchunk)  )

     call infld('GBXAR', fh_topo, dim1name, dim2name, 1, pcols, &
                         begchunk, endchunk, rdg_gbxar, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: GBXAR not found on topo file')

     rdg_gbxar = rdg_gbxar * (rearth/1000._r8)*(rearth/1000._r8) ! transform to km^2

     call infld('HWDTH', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_hwdth, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: HWDTH not found on topo file')

     call infld('CLNGT', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_clngt, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: CLNGT not found on topo file')

     call infld('MXDIS', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_mxdis, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: MXDIS not found on topo file')

     call infld('ANIXY', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_anixy, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: ANIXY not found on topo file')

     call infld('ANGLL', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_angll, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: ANGLL not found on topo file')

     ! close topo file only if it was opened here
     if (len_trim(bnd_topo_loc) > 0) then
        call pio_closefile(fh_topo)
     end if

     call addfld ('WBR_HT1',     horiz_only,  'I','m', &
          'Wave breaking height for DSW')
     call addfld ('TLB_HT1',     horiz_only,  'I','m', &
          'Form drag layer height')
     call addfld ('BWV_HT1',     horiz_only,  'I','m', &
          'Bottom of freely-propagating OGW regime')
     call addfld ('TAUDSW1',     horiz_only,  'I','Nm-2', &
          'DSW enhanced drag')
     call addfld ('TAUORO1',     horiz_only,  'I','Nm-2', &
          'lower BC on propagating wave stress')
     call addfld ('UBMSRC1',     horiz_only,  'I','ms-1', &
          'below-peak-level on-ridge wind')
     call addfld ('USRC1',     horiz_only,  'I','ms-1', &
          'below-peak-level Zonal wind')
     call addfld ('VSRC1',     horiz_only,  'I','ms-1', &
          'below-peak-level Meridional wind')
     call addfld ('NSRC1',     horiz_only,  'I','s-1', &
          'below-peak-level stratification')
     call addfld ('MXDIS1',     horiz_only,  'I','m', &
          'Ridge/obstacle height')
     call addfld ('ANGLL1',     horiz_only,  'I','degrees', &
          'orientation clockwise w/resp north-south')
     call addfld ('ANIXY1',     horiz_only,  'I','1', &
          'Ridge quality')
     call addfld ('HWDTH1',     horiz_only,  'I','km', &
          'Ridge width')
     call addfld ('CLNGT1',     horiz_only,  'I','km', &
          'Ridge length')
     call addfld ('GBXAR1',     horiz_only,  'I','km+2', &
          'grid box area')

     call addfld ('Fr1_DIAG',     horiz_only,  'I','1', &
          'Critical Froude number for linear waves')
     call addfld ('Fr2_DIAG',     horiz_only,  'I','1', &
          'Critical Froude number for blocked flow')
     call addfld ('Frx_DIAG',     horiz_only,  'I','1', &
          'Obstacle Froude Number')

     call addfld('UEGW',  (/ 'lev' /) , 'A'  ,'s-1' ,  &
          'Zonal wind profile-entry to GW ' )
     call addfld('VEGW',  (/ 'lev' /) , 'A'  ,'s-1' ,  &
          'Merdional wind profile-entry to GW ' )
     call register_vector_field('UEGW','VEGW')
     call addfld('TEGW',  (/ 'lev' /) , 'A'  ,'K' ,  &
          'Temperature profile-entry to GW ' )
     call addfld('ZEGW',  (/ 'ilev' /) , 'A'  ,'m' ,  &
          'interface geopotential heights in GW code ' )
     call addfld('ZMGW',  (/ 'lev' /) , 'A'  ,'m' ,  &
          'midlayer geopotential heights in GW code ' )

     call addfld('TAUM1_DIAG' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld('TAU1RDGBETAM' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld('UBM1BETA',  (/ 'lev' /) , 'A'  ,'s-1' ,  &
          'On-ridge wind profile           ' )
     call addfld('UBT1RDGBETA' , (/ 'lev' /) , 'I'  ,'m s-1' , &
          'On-ridge wind tendency from ridge 1     ')

     do i = 1, 6
        write(cn, '(i1)') i
        call addfld('TAU'//cn//'RDGBETAY' , (/ 'ilev' /), 'I', 'N m-2', &
          'Ridge based momentum flux profile')
        call addfld('TAU'//cn//'RDGBETAX' , (/ 'ilev' /), 'I', 'N m-2', &
          'Ridge based momentum flux profile')
        call register_vector_field('TAU'//cn//'RDGBETAX','TAU'//cn//'RDGBETAY')
        call addfld('UT'//cn//'RDGBETA',    (/ 'lev' /),  'I', 'm s-1', &
          'U wind tendency from ridge '//cn)
        call addfld('VT'//cn//'RDGBETA',    (/ 'lev' /),  'I', 'm s-1', &
          'V wind tendency from ridge '//cn)
        call register_vector_field('UT'//cn//'RDGBETA','VT'//cn//'RDGBETA')
     end do

     call addfld('TAUARDGBETAY' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld('TAUARDGBETAX' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call register_vector_field('TAUARDGBETAX','TAUARDGBETAY')

     if (history_waccm) then
        call add_default('TAUARDGBETAX', 1, ' ')
        call add_default('TAUARDGBETAY  ', 1, ' ')
     end if

  end if

  if (use_gw_rdg_gamma) then

     if (effgw_rdg_gamma == unset_r8) then
        call endrun(sub//": ERROR: Anisotropic OGW enabled, but effgw_rdg_gamma was not set.")
     end if

     call getfil(bnd_rdggm, bnd_rdggm_loc, iflag=1, lexist=found)
     if (found) then
        call cam_pio_openfile(fh_rdggm, bnd_rdggm_loc, PIO_NOWRITE)
     else
        call endrun(sub//': ERROR: file for gamma ridges not found: bnd_rdggm='// &
                   trim(bnd_rdggm))
     end if

     if (.not. allocated(rdg_gbxar)) then
        allocate(rdg_gbxar(pcols,begchunk:endchunk))
        call infld('GBXAR', fh_rdggm, dim1name, dim2name, 1, pcols, &
                            begchunk, endchunk, rdg_gbxar, found, gridname='physgrid')
        if (.not. found) call endrun(sub//': ERROR: GBXAR not found on bnd_rdggm')
        rdg_gbxar = rdg_gbxar * (rearth/1000._r8)*(rearth/1000._r8) ! transform to km^2
     end if

     ! Get meso-gamma ridge data
     allocate( &
        rdg_hwdthg(pcols,prdg,begchunk:endchunk), &
        rdg_clngtg(pcols,prdg,begchunk:endchunk), &
        rdg_mxdisg(pcols,prdg,begchunk:endchunk), &
        rdg_anixyg(pcols,prdg,begchunk:endchunk), &
        rdg_angllg(pcols,prdg,begchunk:endchunk)  )

     call infld('HWDTH', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_hwdthg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: HWDTH not found on bnd_rdggm')

     call infld('CLNGT', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_clngtg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: CLNGT not found on bnd_rdggm')

     call infld('MXDIS', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_mxdisg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: MXDIS not found on bnd_rdggm')

     call infld('ANIXY', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_anixyg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: ANIXY not found on bnd_rdggm')

     call infld('ANGLL', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_angllg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: ANGLL not found on bnd_rdggm')

     call pio_closefile(fh_rdggm)

     call addfld ('TAU1RDGGAMMAM' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld ('UBM1GAMMA',  (/ 'lev' /) , 'A'  ,'s-1' ,  &
          'On-ridge wind profile           ' )
     call addfld ('UBT1RDGGAMMA' , (/ 'lev' /) , 'I'  ,'m s-1' , &
          'On-ridge wind tendency from ridge 1     ')

     do i = 1, 6
        write(cn, '(i1)') i
        call addfld('TAU'//cn//'RDGGAMMAY', (/ 'ilev' /), 'I', 'N m-2', &
          'Ridge based momentum flux profile')
        call addfld('TAU'//cn//'RDGGAMMAX', (/ 'ilev' /), 'I', 'N m-2', &
          'Ridge based momentum flux profile')
        call addfld('UT'//cn//'RDGGAMMA' , (/ 'lev' /),  'I', 'm s-1', &
          'U wind tendency from ridge '//cn)
        call addfld('VT'//cn//'RDGGAMMA' , (/ 'lev' /),  'I', 'm s-1', &
          'V wind tendency from ridge '//cn)
        call register_vector_field('UT'//cn//'RDGGAMMA','VT'//cn//'RDGGAMMA')
     end do

     call addfld ('TAUARDGGAMMAY' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld ('TAUARDGGAMMAX' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call register_vector_field('TAUARDGGAMMAX','TAUARDGGAMMAY')
     call addfld ('TAURDGGMX',     horiz_only,  'A','N m-2', &
          'Zonal gravity wave surface stress')
     call addfld ('TAURDGGMY',     horiz_only,  'A','N m-2', &
          'Meridional gravity wave surface stress')
     call register_vector_field('TAURDGGMX','TAURDGGMY')
     call addfld ('UTRDGGM' , (/ 'lev' /) , 'I'  ,'m s-1' , &
          'U wind tendency from ridge 6     ')
     call addfld ('VTRDGGM' , (/ 'lev' /) , 'I'  ,'m s-1' , &
          'V wind tendency from ridge 6     ')
     call register_vector_field('UTRDGGM','VTRDGGM')
  end if

  if (use_gw_front .or. use_gw_front_igw) then

     frontgf_idx = pbuf_get_index('FRONTGF')
     frontga_idx = pbuf_get_index('FRONTGA')

     call shr_assert(unset_r8 /= frontgfc, &
          "gw_init: Frontogenesis enabled, but frontgfc was &
          & not set!"// &
          errMsg(__FILE__, __LINE__))

     call addfld ('FRONTGF', (/ 'lev' /), 'A', 'K^2/M^2/S', &
          'Frontogenesis function at gws src level')
     call addfld ('FRONTGFA', (/ 'lev' /), 'A', 'K^2/M^2/S', &
          'Frontogenesis function at gws src level')

     if (history_waccm) then
        call add_default('FRONTGF', 1, ' ')
        call add_default('FRONTGFA', 1, ' ')
     end if

  end if

  if (use_gw_front) then

     call shr_assert(all(unset_r8 /= [ effgw_cm, taubgnd ]), &
          "gw_init: Frontogenesis mid-scale waves enabled, but not &
          &all required namelist variables were set!"// &
          errMsg(__FILE__, __LINE__))

     if (masterproc) then
        write(iulog,*) 'gw_init: gw spectrum taubgnd, ', &
             'effgw_cm = ',taubgnd, effgw_cm
        write(iulog,*) ' '
     end if


!!$     ! Output for gravity waves from frontogenesis.
!!$     call gw_spec_addflds(prefix=cm_pf, scheme="C&M", band=band_mid, &
!!$          history_defaults=history_waccm)

  end if

  if (use_gw_front_igw) then

     call shr_assert(all(unset_r8 /= [ effgw_cm_igw, taubgnd_igw ]), &
          "gw_init: Frontogenesis inertial waves enabled, but not &
          &all required namelist variables were set!"// &
          errMsg(__FILE__, __LINE__))

     if (masterproc) then
        write(iulog,*) 'gw_init: gw spectrum taubgnd_igw, ', &
             'effgw_cm_igw = ',taubgnd_igw, effgw_cm_igw
        write(iulog,*) ' '
     end if

!!$     ! Output for gravity waves from frontogenesis.
!!$     call gw_spec_addflds(prefix=cm_igw_pf, scheme="C&M IGW", &
!!$          band=band_long, history_defaults=history_waccm)

  end if

  ! ========= Moving Mountain initialization! ==========================
  if (use_gw_movmtn_pbl) then

     ! get pbuf indices for CLUBB couplings
     ttend_clubb_idx     = pbuf_get_index('TTEND_CLUBB')
     thlp2_clubb_gw_idx  = pbuf_get_index('THLP2_CLUBB_GW')
     upwp_clubb_gw_idx   = pbuf_get_index('UPWP_CLUBB_GW')
     vpwp_clubb_gw_idx   = pbuf_get_index('VPWP_CLUBB_GW')
     wpthlp_clubb_gw_idx  = pbuf_get_index('WPTHLP_CLUBB_GW')

     if (masterproc) then
        write (iulog,*) 'Moving Mountain development code call init_movmtn'
     end if


     ! Confirm moving mountain file is enabled
     call shr_assert(trim(gw_drag_file_mm) /= "", &
          "gw_init: No gw_drag_file provided for DP GW moving mountain lookup &
          &table. Set this via namelist."// &
          errMsg(__FILE__, __LINE__))

     call gw_drag_cam_movmtn_init(gw_drag_file_mm, band_movmtn, movmtn_desc)

     do k = 0, pver
        ! 950 hPa index
        if (pref_edge(k+1) < 95000._r8) movmtn_desc%k = k+1
     end do

    ! Don't use deep convection heating depths below this limit.
     movmtn_desc%min_hdepth = 1._r8
     if (masterproc) then
        write (iulog,*) 'Moving mountain deep level =',movmtn_desc%k
     end if

     call addfld ('GWUT_MOVMTN',(/ 'lev' /), 'I','m s-2', &
          'Mov Mtn dragforce - ubm component')
     call addfld ('UTGW_MOVMTN',(/ 'lev' /), 'I','m s-2', &
          'Mov Mtn dragforce - u component')
     call addfld ('VTGW_MOVMTN',(/ 'lev' /), 'I','m s-2', &
          'Mov Mtn dragforce - v component')
     call addfld('TAU_MOVMTN', (/ 'ilev' /), 'I', 'N m-2', &
          'Moving Mountain momentum flux profile')
     call addfld('U_MOVMTN_IN', (/ 'lev' /), 'I', 'm s-1', &
          'Moving Mountain - midpoint zonal input wind')
     call addfld('V_MOVMTN_IN', (/ 'lev' /), 'I', 'm s-1', &
          'Moving Mountain - midpoint meridional input wind')
     call addfld('UBI_MOVMTN', (/ 'ilev' /), 'I', 'm s-1', &
          'Moving Mountain - interface wind in direction of wave')
     call addfld('UBM_MOVMTN', (/ 'lev' /), 'I', 'm s-1', &
          'Moving Mountain - midpoint wind in direction of wave')
     call addfld ('HDEPTH_MOVMTN',horiz_only,'I','km', &
          'Heating Depth')
     call addfld ('UCELL_MOVMTN',horiz_only,'I','m s-1', &
          'Gravity Wave Moving Mountain - Source-level X-wind')
     call addfld ('VCELL_MOVMTN',horiz_only,'I','m s-1', &
          'Gravity Wave Moving Mountain - Source-level Y-wind')
     call addfld ('CS_MOVMTN',horiz_only,'I','m s-1', &
          'Gravity Wave Moving Mountain - phase speed in direction of wave')
     call addfld ('STEER_LEVEL_MOVMTN',horiz_only,'I','1', &
          'Gravity Wave Moving Mountain - steering level for movmtn GW')
     call addfld ('SRC_LEVEL_MOVMTN',horiz_only,'I','1', &
          'Gravity Wave Moving Mountain - launch level for movmtn GW')
     call addfld ('TND_LEVEL_MOVMTN',horiz_only,'I','1', &
          'Gravity Wave Moving Mountain - tendency lowest level for movmtn GW')
     call addfld ('NETDT_MOVMTN',(/ 'lev' /),'I','K s-1', &
          'Gravity Wave Moving Mountain - Net heating rate')
     call addfld ('TTEND_CLUBB',(/ 'lev' /),'A','K s-1', &
          'Gravity Wave Moving Mountain - CLUBB Net heating rate')
     call addfld ('THLP2_CLUBB_GW',(/ 'ilev' /),'A','K+2', &
          'Gravity Wave Moving Mountain - THLP variance from CLUBB to GW')
     call addfld ('WPTHLP_CLUBB_GW',(/ 'ilev' /),'A','Km s-2', &
          'Gravity Wave Moving Mountain - WPTHLP from CLUBB to GW')
     call addfld ('UPWP_CLUBB_GW',(/ 'ilev' /),'A','m+2 s-2', &
          'Gravity Wave Moving Mountain - X-momflux from CLUBB to GW')
     call addfld ('VPWP_CLUBB_GW',(/ 'ilev' /),'A','m+2 s-2', &
          'Gravity Wave Moving Mountain - Y-momflux from CLUBB to GW')
     call addfld ('XPWP_SRC_MOVMTN',horiz_only,'I','m+2 s-2', &
          'Gravity Wave Moving Mountain - flux source for moving mtn')

  end if

  if (use_gw_convect_dp) then

     ttend_dp_idx    = pbuf_get_index('TTEND_DP')

     ! Set the deep scheme specification components.
     beres_dp_desc%storm_shift = .true.

     do k = 0, pver
        ! 700 hPa index
        if (pref_edge(k+1) < 70000._r8) beres_dp_desc%k = k+1
     end do

     if (masterproc) then
        write (iulog,*) 'Beres deep level =',beres_dp_desc%k
     end if

     ! Don't use deep convection heating depths below this limit.
     ! This is important for QBO. Bad result if left at 2.5 km.
     beres_dp_desc%min_hdepth = 1000._r8

     ! Read Beres file.

     call shr_assert(trim(gw_drag_file) /= "", &
          "gw_init: No gw_drag_file provided for Beres deep &
          &scheme. Set this via namelist."// &
          errMsg(__FILE__, __LINE__))

     call gw_drag_cam_beres_init(gw_drag_file, band_mid, beres_dp_desc)

!!$     ! Output for gravity waves from the Beres scheme (deep).
!!$     call gw_spec_addflds(prefix=beres_dp_pf, scheme="Beres (deep)", &
!!$          band=band_mid, history_defaults=history_waccm)

     call addfld ('NETDT',(/ 'lev' /), 'A','K s-1', &
          'Net heating rate')
     call addfld ('MAXQ0',horiz_only  ,  'A','K day-1', &
          'Max column heating rate')
     call addfld ('HDEPTH',horiz_only,    'A','km', &
          'Heating Depth')

     if (history_waccm) then
        call add_default('NETDT    ', 1, ' ')
        call add_default('HDEPTH   ', 1, ' ')
        call add_default('MAXQ0    ', 1, ' ')
     end if

  end if

  if (use_gw_convect_sh) then

     ttend_sh_idx    = pbuf_get_index('TTEND_SH')

     ! Set the shallow scheme specification components.
     beres_sh_desc%storm_shift = .false.

     do k = 0, pver
        ! 900 hPa index
        if (pref_edge(k+1) < 90000._r8) beres_sh_desc%k = k+1
     end do

     if (masterproc) then
        write (iulog,*) 'Beres shallow level =',beres_sh_desc%k
     end if

     ! Use all heating depths for shallow convection.
     beres_sh_desc%min_hdepth = 0._r8

     ! Read Beres file.

     call shr_assert(trim(gw_drag_file_sh) /= "", &
          "gw_init: No gw_drag_file_sh provided for Beres shallow &
          &scheme. Set this via namelist."// &
          errMsg(__FILE__, __LINE__))

     call gw_drag_cam_beres_init(gw_drag_file_sh, band_mid, beres_sh_desc)

!!$     ! Output for gravity waves from the Beres scheme (shallow).
!!$     call gw_spec_addflds(prefix=beres_sh_pf, scheme="Beres (shallow)", &
!!$          band=band_mid, history_defaults=history_waccm)
!!$
     call addfld ('SNETDT',(/ 'lev' /), 'A','K s-1', &
          'Net heating rate')
     call addfld ('SMAXQ0',horiz_only  ,  'A','K day-1', &
          'Max column heating rate')
     call addfld ('SHDEPTH',horiz_only,    'A','km', &
          'Heating Depth')

     if (history_waccm) then
        call add_default('SNETDT   ', 1, ' ')
        call add_default('SHDEPTH  ', 1, ' ')
        call add_default('SMAXQ0   ', 1, ' ')
     end if

  end if

  call addfld ('EKGW' ,(/ 'ilev' /), 'A','M2/S', &
       'Effective Kzz due to diffusion by gravity waves')

  if (history_waccm) then
     call add_default('EKGW', 1, ' ')
  end if

  call addfld ('UTGW_TOTAL',    (/ 'lev' /), 'A','m s-2', &
       'Total U tendency due to gravity wave drag')
  call addfld ('VTGW_TOTAL',    (/ 'lev' /), 'A','m s-2', &
       'Total V tendency due to gravity wave drag')
  call register_vector_field('UTGW_TOTAL', 'VTGW_TOTAL')

  ! Total temperature tendency output.
  call addfld ('TTGW', (/ 'lev' /), 'A', 'K s-1',  &
       'T tendency - gravity wave drag')

  ! Water budget terms.
  call addfld('QTGW',(/ 'lev' /), 'A','kg/kg/s', &
       'Q tendency - gravity wave drag')
  call addfld('CLDLIQTGW',(/ 'lev' /), 'A','kg/kg/s', &
       'CLDLIQ tendency - gravity wave drag')
  call addfld('CLDICETGW',(/ 'lev' /), 'A','kg/kg/s', &
       'CLDICE tendency - gravity wave drag')

  if ( history_budget ) then
     call add_default('TTGW', history_budget_histfile_num, ' ')
     call add_default('QTGW', history_budget_histfile_num, ' ')
     call add_default('CLDLIQTGW', history_budget_histfile_num, ' ')
     call add_default('CLDICETGW', history_budget_histfile_num, ' ')
  end if

  ! Get indices to actually output the above.
  call cnst_get_ind("CLDLIQ", ixcldliq)
  call cnst_get_ind("CLDICE", ixcldice)

  call gw_drag_init( &
       pver,  &
  gravit,  &
  rair,  &
  pi,  &
  pgwv,  &
  gw_dc,  &
  pgwv_long,  &
  gw_dc_long,  &
  tau_0_ubc,  &
  pref_edge,  &
  pref_mid,  &
  gw_bot_taper_pres,  &
  effgw_beres_dp,  &
  effgw_beres_sh,  &
  effgw_cm,  &
  effgw_cm_igw,  &
  effgw_oro,  &
  fcrit2,  &
  frontgfc,  &
  gw_drag_file,  &
  gw_drag_file_sh,  &
  gw_drag_file_mm,  &
  taubgnd,  &
  taubgnd_igw,  &
  gw_polar_taper,  &
  use_gw_rdg_beta,  &
  n_rdg_beta,  &
  effgw_rdg_beta,  &
  effgw_rdg_beta_max,  &
  rdg_beta_cd_llb,  &
  trpd_leewv_rdg_beta,  &
  use_gw_rdg_gamma,  &
  n_rdg_gamma,  &
  effgw_rdg_gamma,  &
  effgw_rdg_gamma_max,  &
  rdg_gamma_cd_llb,  &
  trpd_leewv_rdg_gamma,  &
  bnd_rdggm,  &
  gw_oro_south_fac,  &
  gw_limit_tau_without_eff,  &
  gw_lndscl_sgh,  &
  gw_prndl,  &
  gw_apply_tndmax,  &
  gw_qbo_hdepth_scaling,  &
  gw_top_taper,  &
  front_gaussian_width,  &
  alpha_gw_movmtn,  &
  use_gw_front,  &
  use_gw_oro,  &
  use_gw_front_igw,  &
  use_gw_convect_dp,  &
  use_gw_convect_sh,  &
  use_simple_phys,  &
  use_gw_movmtn_pbl, &
  iulog,  &
  ktop,  &
  masterproc, &
  do_molec_diff,  &
  wavelength_mid,  &
  wavelength_long,  &
!!$  gw_rdg_do_divstream,  &
!!$  gw_rdg_C_BetaMax_DS,  &
!!$  gw_rdg_C_GammaMax,  &
!!$  gw_rdg_Frx0,  &
!!$  gw_rdg_Frx1,  &
!!$  gw_rdg_C_BetaMax_SM,  &
!!$  gw_rdg_Fr_c,  &
!!$  gw_rdg_do_smooth_regimes,  &
!!$  gw_rdg_do_adjust_tauoro,  &
!!$  gw_rdg_do_backward_compat,  &
!!$  gw_rdg_orohmin,  &
!!$  gw_rdg_orovmin,  &
!!$  gw_rdg_orostratmin,  &
!!$  gw_rdg_orom2min,  &
!!$  gw_rdg_do_vdiff,  &
  rdg_gbxar,  &
  rdg_hwdth,  &
  rdg_clngt,  &
  rdg_mxdis,  &
  rdg_anixy,  &
  rdg_angll,  &
  rdg_gbxarg,  &
  rdg_hwdthg,  &
  rdg_clngtg,  &
  rdg_mxdisg,  &
  rdg_anixyg,  &
  rdg_angllg,  &
    beres_dp_desc%storm_shift, &
    beres_dp_desc%k, &
    beres_dp_desc%min_hdepth, &
    beres_dp_desc%maxh, &
    beres_dp_desc%maxuh, &
    beres_dp_desc%hd, &
    beres_dp_desc%mfcc, &
    beres_sh_desc%storm_shift, &
    beres_sh_desc%k, &
    beres_sh_desc%min_hdepth, &
    beres_sh_desc%maxh, &
    beres_sh_desc%maxuh, &
    beres_sh_desc%hd, &
    beres_sh_desc%mfcc, &
    movmtn_desc%storm_shift, &
    movmtn_desc%k, &
    movmtn_desc%min_hdepth, &
    movmtn_desc%maxh, &
    movmtn_desc%maxuh, &
    movmtn_desc%hd, &
    movmtn_desc%uh, &
    movmtn_desc%mfcc, &
    cm_desc%ksrc, &
    cm_desc%kfront, &
    cm_desc%frontgfc, &
    cm_desc%src_tau, &
    cm_igw_desc%ksrc, &
    cm_igw_desc%kfront, &
    cm_igw_desc%frontgfc, &
    cm_igw_desc%src_tau, &
  errormsg,  &
  errorflg )


end subroutine gw_drag_cam_init

!==========================================================================

subroutine gw_drag_cam_beres_init(file_name, band, desc)

  use ioFileMod, only: getfil
  use pio, only: file_desc_t, pio_nowrite, pio_inq_varid, pio_get_var, &
       pio_closefile
  use cam_pio_utils, only: cam_pio_openfile

  character(len=*), intent(in) :: file_name
  type(GWBand), intent(in) :: band

  type(BeresSourceDesc), intent(inout) :: desc

  type(file_desc_t) :: gw_file_desc

  ! PIO variable ids and error code.
  integer :: mfccid, hdid, stat

  ! Number of wavenumbers in the input file.
  integer :: ngwv_file

  ! Full path to gw_drag_file.
  character(len=cl) :: file_path

  character(len=cl) :: msg

  !----------------------------------------------------------------------
  ! read in look-up table for source spectra
  !-----------------------------------------------------------------------

  call getfil(file_name, file_path)
  call cam_pio_openfile(gw_file_desc, file_path, pio_nowrite)

  ! Get HD (heating depth) dimension.

  desc%maxh = get_pio_dimlen(gw_file_desc, "HD", file_path)

  ! Get MW (mean wind) dimension.

  desc%maxuh = get_pio_dimlen(gw_file_desc, "MW", file_path)

  ! Get PS (phase speed) dimension.

  ngwv_file = get_pio_dimlen(gw_file_desc, "PS", file_path)

  ! Number in each direction is half of total (and minus phase speed of 0).
  desc%maxuh = (desc%maxuh-1)/2
  ngwv_file = (ngwv_file-1)/2

  call shr_assert(ngwv_file >= band%ngwv, &
       "gw_drag_cam_beres_init: PhaseSpeed in lookup table file does not cover the whole &
       &spectrum implied by the model's ngwv. ")

  ! Allocate hd and get data.

  allocate(desc%hd(desc%maxh), stat=stat, errmsg=msg)

  call shr_assert(stat == 0, &
       "gw_drag_cam_beres_init: Allocation error (hd): "//msg// &
       errMsg(__FILE__, __LINE__))

  stat = pio_inq_varid(gw_file_desc,'HD',hdid)

  call handle_pio_error(stat, &
       'Error finding HD in: '//trim(file_path))

  stat = pio_get_var(gw_file_desc, hdid, start=[1], count=[desc%maxh], &
       ival=desc%hd)

  call handle_pio_error(stat, &
       'Error reading HD from: '//trim(file_path))

  ! While not currently documented in the file, it uses kilometers. Convert
  ! to meters.
  desc%hd = desc%hd*1000._r8

  ! Allocate mfcc. "desc%maxh" and "desc%maxuh" are from the file, but the
  ! model determines wavenumber dimension.

  allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,&
       -band%ngwv:band%ngwv), stat=stat, errmsg=msg)

  call shr_assert(stat == 0, &
       "gw_drag_cam_beres_init: Allocation error (mfcc): "//msg// &
       errMsg(__FILE__, __LINE__))

  ! Get mfcc data.

  stat = pio_inq_varid(gw_file_desc,'mfcc',mfccid)

  call handle_pio_error(stat, &
       'Error finding mfcc in: '//trim(file_path))

  stat = pio_get_var(gw_file_desc, mfccid, &
       start=[1,1,ngwv_file-band%ngwv+1], count=shape(desc%mfcc), &
       ival=desc%mfcc)

  call handle_pio_error(stat, &
       'Error reading mfcc from: '//trim(file_path))

  call pio_closefile(gw_file_desc)

  if (masterproc) then

     write(iulog,*) "Read in source spectra from file."
     write(iulog,*) "mfcc max, min = ", &
          maxval(desc%mfcc),", ",minval(desc%mfcc)

  endif

end subroutine gw_drag_cam_beres_init

!==============================================================
subroutine gw_drag_cam_movmtn_init(file_name, band, desc)

  use ioFileMod, only: getfil
  use pio, only: file_desc_t, pio_nowrite, pio_inq_varid, pio_get_var, &
       pio_closefile
  use cam_pio_utils, only: cam_pio_openfile

  character(len=*), intent(in) :: file_name
  type(GWBand), intent(in) :: band

  type(MovMtnSourceDesc), intent(inout) :: desc

  type(file_desc_t) :: gw_file_desc

  ! PIO variable ids and error code.
  integer :: mfccid, uhid, hdid, stat

  ! Number of wavenumbers in the input file.
  integer :: ngwv_file

  ! Full path to gw_drag_file.
  character(len=cl) :: file_path

  character(len=cl) :: msg

  !----------------------------------------------------------------------
  ! read in look-up table for source spectra
  !-----------------------------------------------------------------------

  call getfil(file_name, file_path)

  call cam_pio_openfile(gw_file_desc, file_path, pio_nowrite)

  ! Get HD (heating depth) dimension.

  desc%maxh = 15 !get_pio_dimlen(gw_file_desc, "HD", file_path)

  ! Get MW (mean wind) dimension.

 desc%maxuh = 241 ! get_pio_dimlen(gw_file_desc, "MW", file_path)

  ! Get PS (phase speed) dimension.

  ngwv_file = 0 !get_pio_dimlen(gw_file_desc, "PS", file_path)

  ! Number in each direction is half of total (and minus phase speed of 0).
  desc%maxuh = (desc%maxuh-1)/2
  ngwv_file = (ngwv_file-1)/2

  call shr_assert(ngwv_file >= band%ngwv, &
       "gw_drag_cam_movmtn_init: PhaseSpeed in lookup table inconsistent with moving mountain")

  ! Allocate hd and get data.

  allocate(desc%hd(desc%maxh), stat=stat, errmsg=msg)

  call shr_assert(stat == 0, &
       "gw_drag_cam_movmtn_init: Allocation error (hd): "//msg// &
       errMsg(__FILE__, __LINE__))

  stat = pio_inq_varid(gw_file_desc,'HDEPTH',hdid)

  call handle_pio_error(stat, &
       'Error finding HD in: '//trim(file_path))

  stat = pio_get_var(gw_file_desc, hdid, start=[1], count=[desc%maxh], &
       ival=desc%hd)

  call handle_pio_error(stat, &
       'Error reading HD from: '//trim(file_path))

  ! While not currently documented in the file, it uses kilometers. Convert
  ! to meters.
  desc%hd = desc%hd*1000._r8

 ! Allocate wind and get data.

  allocate(desc%uh(desc%maxuh), stat=stat, errmsg=msg)

  call shr_assert(stat == 0, &
       "gw_drag_cam_movmtn_init: Allocation error (uh): "//msg// &
       errMsg(__FILE__, __LINE__))

  stat = pio_inq_varid(gw_file_desc,'UARR',uhid)

  call handle_pio_error(stat, &
       'Error finding UH in: '//trim(file_path))

  stat = pio_get_var(gw_file_desc, uhid, start=[1], count=[desc%maxuh], &
       ival=desc%uh)

  call handle_pio_error(stat, &
       'Error reading UH from: '//trim(file_path))

  ! Allocate mfcc. "desc%maxh" and "desc%maxuh" are from the file, but the
  ! model determines wavenumber dimension.

  allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,&
       -band%ngwv:band%ngwv), stat=stat, errmsg=msg)

  call shr_assert(stat == 0, &
       "gw_drag_cam_movmtn_init: Allocation error (mfcc): "//msg// &
       errMsg(__FILE__, __LINE__))

  ! Get mfcc data.

  stat = pio_inq_varid(gw_file_desc,'NEWMF',mfccid)

  call handle_pio_error(stat, &
       'Error finding mfcc in: '//trim(file_path))

  stat = pio_get_var(gw_file_desc, mfccid, &
       start=[1,1], count=shape(desc%mfcc), &
       ival=desc%mfcc)

  call handle_pio_error(stat, &
       'Error reading mfcc from: '//trim(file_path))

  call pio_closefile(gw_file_desc)

  if (masterproc) then

     write(iulog,*) "Read in Mov Mountain source file."

  endif

end subroutine gw_drag_cam_movmtn_init
!==========================================================================

! Utility to reduce the repetitiveness of reads during initialization.
function get_pio_dimlen(file_desc, dim_name, file_path) result(dimlen)

  use pio, only: file_desc_t, pio_inq_dimid, pio_inq_dimlen

  type(file_desc_t), intent(in) :: file_desc
  character(len=*), intent(in) :: dim_name

  ! File path, for use in error messages only.
  character(len=*), intent(in) :: file_path

  integer :: dimlen

  integer :: dimid, stat

  stat = pio_inq_dimid(file_desc, dim_name, dimid)

  call handle_pio_error(stat, &
       "Error finding dimension "//dim_name//" in: "//file_path)

  stat = pio_inq_dimlen(file_desc, dimid, dimlen)

  call handle_pio_error(stat, &
       "Error reading dimension "//dim_name//" from: "//file_path)

end function get_pio_dimlen

!==========================================================================

! In fact, we'd usually expect PIO errors to abort the run before you can
! even check the error code. But just in case, use this little assert.
subroutine handle_pio_error(stat, message)
  use pio, only: pio_noerr
  integer, intent(in) :: stat
  character(len=*) :: message

  call shr_assert(stat == pio_noerr, &
       "PIO error:"//trim(message)// &
       errMsg(__FILE__, __LINE__))

end subroutine handle_pio_error

!==========================================================================

subroutine gw_drag_cam_tend(state, pbuf, dt, ptend, cam_in, flx_heat)
  !-----------------------------------------------------------------------
  ! Interface for multiple gravity wave drag parameterization.
  !-----------------------------------------------------------------------

  use physics_types,   only: physics_state_copy, set_dry_to_wet
  use constituents,    only: cnst_type
  use physics_buffer,  only: physics_buffer_desc, pbuf_get_field
  use camsrfexch,      only: cam_in_t
  ! Location-dependent cpair
  use air_composition, only: cpairv
  use physconst,       only: pi
  use time_manager,    only: get_step_size
  !------------------------------Arguments--------------------------------
  type(physics_state), intent(in) :: state   ! physics state structure
  type(physics_buffer_desc), pointer :: pbuf(:) ! Physics buffer
  real(r8), intent(in) :: dt                    ! time step
  ! Parameterization net tendencies.
  type(physics_ptend), intent(out):: ptend
  type(cam_in_t), intent(in) :: cam_in
  real(r8), intent(out) :: flx_heat(pcols)

  !---------------------------Local storage-------------------------------

  type(physics_state) :: state1     ! Local copy of state variable

  integer :: lchnk                  ! chunk identifier
  integer :: ncol                   ! number of atmospheric columns
  integer :: istat

  integer :: i, k                   ! loop indices


  real(r8) :: ttgw(state%ncol,pver) ! temperature tendency
  real(r8) :: utgw(state%ncol,pver) ! zonal wind tendency
  real(r8) :: vtgw(state%ncol,pver) ! meridional wind tendency

  real(r8) :: ni(state%ncol,pver+1) ! interface Brunt-Vaisalla frequency
  real(r8) :: nm(state%ncol,pver)   ! midpoint Brunt-Vaisalla frequency
  real(r8) :: rhoi(state%ncol,pver+1)     ! interface density
  real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
  real(r8) :: tau0x(state%ncol)     ! c=0 sfc. stress (zonal)
  real(r8) :: tau0y(state%ncol)     ! c=0 sfc. stress (meridional)
  real(r8) :: ubi(state%ncol,pver+1)! projection of wind at interfaces
  real(r8) :: ubm(state%ncol,pver)  ! projection of wind at midpoints
  real(r8) :: xv(state%ncol)        ! unit vector of source wind (x)
  real(r8) :: yv(state%ncol)        ! unit vector of source wind (y)

  integer :: m                      ! dummy integers
  real(r8) :: qtgw(state%ncol,pver,pcnst) ! constituents tendencies

  ! Reynolds stress for waves propagating in each cardinal direction.
  real(r8) :: taucd(state%ncol,pver+1,4)

  ! gravity wave wind tendency for each wave
  real(r8), allocatable :: gwut(:,:,:)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8) :: dttdf(state%ncol,pver)
  real(r8) :: dttke(state%ncol,pver)

  ! Wave phase speeds for each column
  real(r8), allocatable :: phase_speeds(:,:)

  ! Efficiency for a gravity wave source.
  real(r8) :: effgw(state%ncol)

  ! pbuf fields
  ! Molecular diffusivity
  real(r8), pointer :: kvt(:,:)
  real(r8) :: kvtt(state%ncol,pver+1)

  ! Frontogenesis
  real(r8), pointer :: frontgf(:,:)
  real(r8), pointer :: frontga(:,:)

  ! Temperature change due to deep convection.
  real(r8), pointer :: ttend_dp(:,:)
  ! Temperature change due to shallow convection.
  real(r8), pointer :: ttend_sh(:,:)

  !  New couplings from CLUBB
  real(r8), pointer :: ttend_clubb(:,:)
  real(r8), pointer :: thlp2_clubb_gw(:,:)
  real(r8), pointer :: wpthlp_clubb_gw(:,:)
  real(r8), pointer :: upwp_clubb_gw(:,:)
  real(r8), pointer :: vpwp_clubb_gw(:,:)
  real(r8) :: xpwp_clubb(state%ncol,pver+1)


  ! Standard deviation of orography.
  real(r8), pointer :: sgh(:)

  ! gridbox area
  real(r8), pointer :: gbxar(:)

     ! Beta ridges
  ! width of ridges.
  real(r8), pointer :: hwdth(:,:)
  ! length of ridges.
  real(r8), pointer :: clngt(:,:)
  ! Maximum deviations of ridges.
  real(r8), pointer :: mxdis(:,:)
  ! orientation of ridges.
  real(r8), pointer :: angll(:,:)
  ! anisotropy of ridges.
  real(r8), pointer :: anixy(:,:)

     ! Gamma ridges
  ! width of ridges.
  real(r8), pointer :: hwdthg(:,:)
  ! length of ridges.
  real(r8), pointer :: clngtg(:,:)
  ! Maximum deviations of ridges.
  real(r8), pointer :: mxdisg(:,:)
  ! orientation of ridges.
  real(r8), pointer :: angllg(:,:)
  ! anisotropy of ridges.
  real(r8), pointer :: anixyg(:,:)

  ! Indices of gravity wave source and lowest level where wind tendencies
  ! are allowed.
  integer :: src_level(state%ncol)
  integer :: tend_level(state%ncol)

  ! Convective source heating depth.
  ! heating depth
  real(r8) :: hdepth(state%ncol)
  ! maximum heating rate
  real(r8) :: maxq0(state%ncol)

  ! Scale sgh to account for landfrac.
  real(r8) :: sgh_scaled(state%ncol)

  ! Parameters for the IGW polar taper.
  real(r8), parameter :: degree2radian = pi/180._r8
  real(r8), parameter :: al0 = 82.5_r8 * degree2radian
  real(r8), parameter :: dlat0 = 5.0_r8 * degree2radian

  ! effective gw diffusivity at interfaces needed for output
  real(r8) :: egwdffi(state%ncol,pver+1)
  ! sum from the two types of spectral GW
  real(r8) :: egwdffi_tot(state%ncol,pver+1)

  ! Momentum fluxes used by fixer.
  real(r8) :: um_flux(state%ncol), vm_flux(state%ncol)
  ! Energy change used by fixer.
  real(r8) :: de(state%ncol)

  ! Which constituents are being affected by diffusion.
  logical  :: lq(pcnst)

  ! Contiguous copies of state arrays.
  real(r8) :: dse(state%ncol,pver)
  real(r8) :: t(state%ncol,pver)
  real(r8) :: u(state%ncol,pver)
  real(r8) :: v(state%ncol,pver)
  real(r8) :: q(state%ncol,pver,pcnst)
  real(r8) :: piln(state%ncol,pver+1)
  real(r8) :: zm(state%ncol,pver)
  real(r8) :: zi(state%ncol,pver+1)

  !------------------------------------------------------------------------
  dt = get_step_size()
  ! Make local copy of input state.
  call physics_state_copy(state, state1)

  ! constituents are all treated as wet mmr
  call set_dry_to_wet(state1, convert_cnst_type='dry')

  dse = state1%s(:ncol,:)
  t = state1%t(:ncol,:)
  u = state1%u(:ncol,:)
  v = state1%v(:ncol,:)
  q = state1%q(:ncol,:,:)
  piln = state1%lnpint(:ncol,:)
  zm = state1%zm(:ncol,:)
  zi = state1%zi(:ncol,:)

  lq = .true.
  call physics_ptend_init(ptend, state1%psetcols, "Gravity wave drag", &
       ls=.true., lu=.true., lv=.true., lq=lq)

  if (do_molec_diff) then
     call pbuf_get_field(pbuf, kvt_idx, kvt)  ! kvt_in(1:pcols,1:pver+1)
  end if

  if (use_gw_front .or. use_gw_front_igw) then
     ! Get frontogenesis physics buffer fields set by dynamics.
     call pbuf_get_field(pbuf, frontgf_idx, frontgf)
     call pbuf_get_field(pbuf, frontga_idx, frontga)
  end if

  if (use_gw_movmtn_pbl) then
     ! Set up heating
     call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)

     !   New couplings from CLUBB
     call pbuf_get_field(pbuf, ttend_clubb_idx, ttend_clubb)
     call pbuf_get_field(pbuf, thlp2_clubb_gw_idx, thlp2_clubb_gw)
     call pbuf_get_field(pbuf, wpthlp_clubb_gw_idx, wpthlp_clubb_gw)
     call pbuf_get_field(pbuf, upwp_clubb_gw_idx, upwp_clubb_gw)
     call pbuf_get_field(pbuf, vpwp_clubb_gw_idx, vpwp_clubb_gw)
  end if
  if (use_gw_convect_dp) then
     ! Set up heating
     call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)
  end if
  if (use_gw_convect_sh) then
     ! Set up heating
     call pbuf_get_field(pbuf, ttend_sh_idx, ttend_sh)
  end if
  if (use_gw_oro) then
     call pbuf_get_field(pbuf, sgh_idx, sgh)
  end if

  call gw_drag_run( &
     ncol, pver, dt, cpair, cpairv, pi, frontgf, frontga, &
     degree2radian,al0,dlat0, &
     pint, piln, pdel, pdeldry, zm, zi, lat, cam_in%landfrac, &
     state_s, state_t, state_u, state_v, state_q, &
     landfrac, sgh, kvt, ttend_dp, ttend_sh, ttend_clubb, &
     thlp2_clubb_gw,wpthlp_clubb_gw,upwp_clubb_gw, vpwp_clubb_gw, &
     ptend%s, ptend%q, ptend%u, ptend%v, scheme_name, nbot_molec, flx_heat, &
     errormsg, errorflg)

  ! Convert the tendencies for the dry constituents to dry air basis.
  do m = 1, pcnst
     if (cnst_type(m).eq.'dry') then
        do k = 1, pver
           do i = 1, ncol
              ptend%q(i,k,m) = ptend%q(i,k,m)*state1%pdel(i,k)/state1%pdeldry(i,k)
           end do
        end do
     end if
  end do

  ! Write totals to history file.
  call outfld('EKGW', egwdffi_tot , ncol, lchnk)
  call outfld('TTGW', ptend%s/cpairv(:,:,lchnk),  pcols, lchnk)

  call outfld('UTGW_TOTAL', ptend%u, pcols, lchnk)
  call outfld('VTGW_TOTAL', ptend%v, pcols, lchnk)

  call outfld('QTGW', ptend%q(:,:,1), pcols, lchnk)
  call outfld('CLDLIQTGW', ptend%q(:,:,ixcldliq), pcols, lchnk)
  call outfld('CLDICETGW', ptend%q(:,:,ixcldice), pcols, lchnk)


end subroutine gw_drag_cam_tend


!==========================================================================

end module gw_drag_cam

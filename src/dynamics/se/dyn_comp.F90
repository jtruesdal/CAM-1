Module dyn_comp

  use bndry_mod,               only: bndry_exchangev
  use cam_abortutils,          only: endrun
  use cam_control_mod,         only: initial_run
  use cam_grid_support,        only: cam_grid_id, cam_grid_get_gcid, &
                                  cam_grid_dimensions, cam_grid_get_dim_names, &
                                  cam_grid_get_latvals, cam_grid_get_lonvals,  &
                                  max_hcoordname_len
  use cam_history_support,     only: max_fieldname_len
  use cam_initfiles,           only: initial_file_get_id, topo_file_get_id, pertlim
  use cam_logfile,             only: iulog
  use cam_map_utils,           only: iMap
  use dimensions_mod,          only: nelemd, nlev, np, npsq
  use dyn_grid,                only: timelevel, dom_mt, hvcoord, ini_grid_hdim_name
  use dyn_grid,                only: get_horiz_grid_dim_d, dyn_decomp, fv_nphys, ini_grid_name
  use edge_mod,                only: edgevpack_nlyr, edgevunpack_nlyr, edge_g
  use element_mod,             only: element_t
  use element_state,           only: elem_state_t
  use hybrid_mod,              only: hybrid_t
  use inic_analytic,           only: analytic_ic_active, analytic_ic_set_ic
  use ncdio_atm,               only: infld
  use parallel_mod,            only: par, initmp
  use perf_mod,                only: t_startf, t_stopf
  use physconst,               only: pi
  use pio,                     only: file_desc_t, io_desc_t, pio_double, PIO_BCAST_ERROR, &
                                     pio_get_local_array_size, pio_freedecomp, PIO_NOERR, &
                                     var_desc_t, PIO_inq_varid, pio_inq_dimid, pio_inq_dimlen, &
                                     pio_seterrorhandling
  use shr_kind_mod,            only: r8 => shr_kind_r8, shr_kind_cl
  use shr_const_mod,           only: SHR_CONST_PI
  use shr_sys_mod,             only: shr_sys_flush
  use spmd_utils,              only: iam, npes_cam => npes, masterproc
  use thread_mod,              only: nthreads, hthreads, vthreads, omp_get_max_threads, omp_get_thread_num
  use time_mod,                only: nsplit
  use time_manager,            only: is_first_step

  implicit none
  private
  save


public ::          &
     dyn_import_t, &
     dyn_export_t, &
     dyn_readnl,   &
     dyn_register, &
     dyn_init,     &
     dyn_run,      &
     dyn_final


  type dyn_import_t
     type (element_t), pointer :: elem(:) => null()
  end type dyn_import_t

  type dyn_export_t
     type (element_t), pointer :: elem(:) => null()
  end type dyn_export_t

  integer, parameter  ::  DYN_RUN_SUCCESS           = 0
  integer, parameter  ::  DYN_RUN_FAILURE           = -1

  ! !DESCRIPTION: This module implements the SE Dynamical Core as
  !               an ESMF gridded component.  It is specific to SE
  !               and does not use ESMF.
  !
  ! \paragraph{Overview}
  !
  !   This module contains an ESMF wrapper for the SE
  !   Dynamical Core used in the Community Atmospheric Model.
  !
  ! !REVISION HISTORY:
  !
  !  JPE  06.05.31:  created
  !  Aaron Donahue 17.04.11: Fixed bug in write_grid_mapping which caused
  !       a segmentation fault when dyn_npes<npes
  !  MT 2020.06.30: remove write_grid_mapping - moved to homme/src/tool
  !
  !----------------------------------------------------------------------

  ! Enumeration of DYNAMICS_IN_COUPLINGS


!jt  logical, parameter         :: DEBUG = .true.

  real(r8), parameter        :: ONE    = 1.0_r8

  character(*), parameter, public :: MODULE_NAME = "dyn_comp"
  character(*), parameter, public :: VERSION     = "$Id$"
!jt  type (domain1d_t), pointer, public :: dom_mt(:) => null()

  ! Frontogenesis indices
  integer, public :: frontgf_idx = -1
  integer, public :: frontga_idx = -1

  interface read_dyn_var
  module procedure read_dyn_field_2d
  module procedure read_dyn_field_3d
end interface read_dyn_var

real(r8), parameter :: rad2deg = 180.0_r8 / pi
real(r8), parameter :: deg2rad = pi / 180.0_r8

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dyn_readnl(NLFileName)
  use namelist_utils, only: find_group_name
  use namelist_mod,   only: homme_set_defaults, homme_postprocess_namelist
  use units,          only: getunit, freeunit
  use spmd_utils,     only: masterproc, masterprocid, mpicom, npes
  use spmd_utils,     only: mpi_real8, mpi_integer, mpi_character, mpi_logical
  use control_mod,    only: hypervis_order, hypervis_subcycle
  use control_mod,    only: hypervis_subcycle_q, integration, statefreq, runtype
  use control_mod,    only: nu, nu_div, nu_p, nu_q, nu_top, qsplit, rsplit
  use control_mod,    only: vert_remap_q_alg, tstep_type, rk_stage_user
  use control_mod,    only: ftype, limiter_option, partmethod
  use control_mod,    only: topology
!jt  use control_mod,    only: fine_ne, hypervis_power, hypervis_scaling
  use control_mod,    only: hypervis_scaling
!jt  use control_mod,    only: max_hypervis_courant
!jt  use fvm_mod,        only: fvm_ideal_test, fvm_test_type
!jt  use fvm_mod,        only: fvm_get_test_type
  use dimensions_mod, only: qsize, qsize_d, ntrac, ntrac_d, npsq, ne, npart, lcp_moist
  use constituents,   only: pcnst
  use params_mod,     only: SFCURVE
!jt  use native_mapping, only: native_mapping_readnl
!!XXgoldyXX: v For future CSLAM / physgrid commit
!    use dp_grids,       only: fv_nphys, fv_nphys2, nphys_pts, write_phys_grid, phys_grid_file
!!XXgoldyXX: ^ For future CSLAM / physgrid commit

  ! Dummy argument
  character(len=*), intent(in) :: NLFileName

  ! Local variables
  integer                      :: unitn, ierr

  ! SE Namelist variables
!jt  integer                      :: se_fine_ne
  integer                      :: se_ftype
  integer                      :: se_hypervis_order
!jt  real(r8)                     :: se_hypervis_power
  real(r8)                     :: se_hypervis_scaling
  integer                      :: se_hypervis_subcycle
  integer                      :: se_hypervis_subcycle_q
  logical                      :: se_lcp_moist
  integer                      :: se_limiter_option
!jt  real(r8)                     :: se_max_hypervis_courant
  character(len=SHR_KIND_CL)   :: se_mesh_file
  integer                      :: se_ne
  integer                      :: se_npes
  integer                      :: se_nsplit
  real(r8)                     :: se_nu
  real(r8)                     :: se_nu_div
  real(r8)                     :: se_nu_p
  real(r8)                     :: se_nu_s
  real(r8)                     :: se_nu_q
  real(r8)                     :: se_nu_top
  integer                      :: se_qsplit
!jt  logical                      :: se_refined_mesh
  integer                      :: se_rsplit
  integer                      :: se_statefreq
  integer                      :: se_tstep_type
  integer                      :: se_vert_remap_q_alg
!!XXgoldyXX: v For future CSLAM / physgrid commit
!    character(len=METHOD_LEN)     :: se_tracer_transport_method
!    character(len=METHOD_LEN)     :: se_cslam_ideal_test
!    character(len=METHOD_LEN)     :: se_cslam_test_type
!    character(len=METHOD_LEN)     :: se_write_phys_grid
!    character(len=shr_kind_cl)    :: se_phys_grid_file
!    integer                       :: se_fv_nphys = 0
!!XXgoldyXX: ^ For future CSLAM / physgrid commit

  namelist /dyn_se_inparm/      &
!jt       se_fine_ne,              & ! For refined meshes
       se_ftype,                & ! forcing type
       se_hypervis_order,       &
!jt       se_hypervis_power,       &
       se_hypervis_scaling,     &
       se_hypervis_subcycle,    &
       se_hypervis_subcycle_q,  &
       se_lcp_moist,       &
       se_limiter_option,       &
!jt       se_max_hypervis_courant, &
       se_mesh_file,            & ! Refined mesh definition file
       se_ne,                   &
       se_npes,                 &
       se_nsplit,               & ! # of dynamics steps per physics timestep
       se_nu,                   &
       se_nu_div,               &
       se_nu_p,                 &
       se_nu_s,                 &
       se_nu_q,                 &
       se_nu_top,               &
       se_qsplit,               &
!jt       se_refined_mesh,         &
       se_rsplit,               &
       se_statefreq,            & ! number of steps per printstate call
       se_tstep_type,           &
       se_vert_remap_q_alg
!!XXgoldyXX: v For future physgrid commit
!         se_fv_nphys,          & ! Linear size of FV physics grid
!         se_write_phys_grid,   &
!         se_phys_grid_file,    &
!!XXgoldyXX: ^ For future physgrid commit

!!XXgoldyXX: v For future CSLAM / physgrid commit
!    namelist /cslam_nl/ se_tracer_transport_method, se_cslam_ideal_test, se_cslam_test_type
!!XXgoldyXX: ^ For future CSLAM / physgrid commit

  !--------------------------------------------------------------------------

 ! namelist default values should be in namelist (but you know users . . .)
  ! NB: Of course, these should keep up with what is in namelist_defaults ...
!jt  se_fine_ne              = -1
  se_ftype                = 0
  se_hypervis_order       = 2
!jt  se_hypervis_power       = 0
  se_hypervis_scaling     = 0
  se_hypervis_subcycle    = 3
  se_hypervis_subcycle_q  = 1
  se_limiter_option       = 8
!jt  se_max_hypervis_courant = 1.0e99_r8
  se_mesh_file            = ''
  se_ne                   = -1
  se_npes                 = npes
  se_nsplit               = 2
  se_nu                   = 1.0e15_r8
  se_nu_div               = 2.5e15_r8
  se_nu_p                 = 1.0e15_r8
  se_nu_q                 = -1.0_r8
  se_nu_top               = 2.5e5_r8
  se_qsplit               = 1
!jt  se_refined_mesh         = .false.
  se_rsplit               = 3
  se_statefreq            = 480
  se_tstep_type           = 5
  se_vert_remap_q_alg     = 1
!!XXgoldyXX: v For future CSLAM / physgrid commit
!    character(len=METHOD_LEN)     :: se_tracer_transport_method
!    character(len=METHOD_LEN)     :: se_cslam_ideal_test
!    character(len=METHOD_LEN)     :: se_cslam_test_type
!    character(len=METHOD_LEN)     :: se_write_phys_grid
!    character(len=shr_kind_cl)    :: se_phys_grid_file
!    integer                       :: se_fv_nphys = 0
!!XXgoldyXX: ^ For future CSLAM / physgrid commit

  ! Read the namelist (dyn_se_inparm)
  call MPI_barrier(mpicom, ierr)
  if (masterproc) then
    write(iulog, *) "dyn_readnl: reading dyn_se_inparm namelist..."
    unitn = getunit()
    open( unitn, file=trim(NLFileName), status='old' )
    call find_group_name(unitn, 'dyn_se_inparm', status=ierr)
    if (ierr == 0) then
      read(unitn, dyn_se_inparm, iostat=ierr)
      if (ierr /= 0) then
        call endrun('dyn_readnl: ERROR reading dyn_se_inparm namelist')
      end if
    end if
    close(unitn)
    call freeunit(unitn)
#ifndef _USEMETIS
      ! override METIS options to SFCURVE
      if (partmethod>=0 .and. partmethod<=3) partmethod=SFCURVE
#endif
       ! ========================
       ! if this is a restart run
       ! ========================
       if(runtype .eq. 1) then
          write(iulog,*)"readnl: restartfile = ",restartfile
       else if(runtype < 0) then
          write(iulog,*)'readnl: runtype=', runtype,' interpolation mode '
       endif


       if((integration .ne. "explicit").and.(integration .ne. "runge_kutta").and. &
                    (integration .ne. "full_imp")) then
          call abortmp('integration must be explicit, full_imp, or runge_kutta')
       end if

       if (integration == "full_imp") then
          if (tstep_type<10) then
             ! namelist did not set a valid tstep_type. pick one:
             tstep_type=11   ! backward euler
             !tstep_type=12  ! BDF2 with BE bootstrap
          endif
       endif

       ierr = timestep_make_subcycle_parameters_consistent(par, rsplit, qsplit, &
            dt_remap_factor, dt_tracer_factor)

       limiter_option=se_limiter_option
       partmethod = se_partmethod
       ne         = se_ne
       npes       = se_npes
       ne_x       = se_ne_x
       ne_y       = se_ne_y
       Lx         = se_lx
       Ly         = se_ly
       topology   = se_topology
       geometry   = se_geometry
       qsize      = qsize_d
       nsplit     = se_nsplit
       tstep      = se_tstep
       if (tstep > 0) then
          if (par%masterproc .and. nsplit > 0) then
             write(iulog,'(a,i3,a)') &
                  'se_tstep and se_nsplit were specified; changing se_nsplit from ', &
                  nsplit, ' to -1.'
          end if
          nsplit = -1
       end if
  end if

  call MPI_barrier(par%comm,ierr)

  npart  = par%nprocs

  ! Broadcast namelist values to all PEs
!jt  call MPI_bcast(se_fine_ne, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_ftype, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_hypervis_order, 1, mpi_integer, masterprocid, mpicom, ierr)
!jt  call MPI_bcast(se_hypervis_power, 1, mpi_real8, masterprocid, mpicom, ierr)
  call MPI_bcast(se_hypervis_scaling, 1, mpi_real8, masterprocid, mpicom, ierr)
  call MPI_bcast(se_hypervis_subcycle, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_hypervis_subcycle_q, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_limiter_option, 1, mpi_integer, masterprocid, mpicom, ierr)
!jt  call MPI_bcast(se_max_hypervis_courant, 1, mpi_real8, masterprocid, mpicom, ierr)
  call MPI_bcast(se_mesh_file, SHR_KIND_CL,  mpi_character, masterprocid, mpicom, ierr)
  call MPI_bcast(se_ne, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_npes, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_nsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_nu, 1, mpi_real8, masterprocid, mpicom, ierr)
  call MPI_bcast(se_nu_div, 1, mpi_real8, masterprocid, mpicom, ierr)
  call MPI_bcast(se_nu_p, 1, mpi_real8, masterprocid, mpicom, ierr)
  call MPI_bcast(se_nu_s, 1, mpi_real8, masterprocid, mpicom, ierr)
  call MPI_bcast(se_nu_q, 1, mpi_real8, masterprocid, mpicom, ierr)
  call MPI_bcast(se_nu_top, 1, mpi_real8, masterprocid, mpicom, ierr)
  call MPI_bcast(se_qsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
!jt  call MPI_bcast(se_refined_mesh, 1, mpi_logical, masterprocid, mpicom, ierr)
  call MPI_bcast(se_rsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_statefreq, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_tstep_type, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(se_vert_remap_q_alg, 1, mpi_integer, masterprocid, mpicom, ierr)
!!XXgoldyXX: v For future physgrid commit
!    call MPI_bcast(fv_nphys, 1, mpi_integer, masterprocid, mpicom, ierr)
!    call MPI_bcast(write_phys_grid, 80,  mpi_character, masterprocid, mpicom, ierr)
!    call MPI_bcast(phys_grid_file,  256, mpi_character, masterprocid, mpicom, ierr)
!    fv_nphys2 = fv_nphys * fv_nphys
!    if (fv_nphys > 0) then
!      nphys_pts = fv_nphys2
!    else
!      nphys_pts = npsq
!    end if
!!XXgoldyXX: ^ For future physgrid commit

  ! Initialize the SE structure that holds the MPI decomposition information
  if (se_npes <= 0) then
    se_npes = npes
  end if

!jt initialize nh dycore
  par = initmp(se_npes)


  ! Fix up unresolved default values
  ! default diffusion coefficiets
  if (se_nu_q < 0) then
    se_nu_q = se_nu
  end if
  if (se_nu_div < 0) then
    se_nu_div = se_nu
  end if
!!$  ! Go ahead and enforce ne = 0 for refined mesh runs
!!$  if (se_refined_mesh) then
!!$    se_ne = 0
!!$  end if

   ! Set HOMME defaults
   call homme_set_defaults()
   ! Set HOMME variables not in CAM's namelist but with different CAM defaults
   partmethod               = SFCURVE
   npart                    = se_npes
   ! CAM requires forward-in-time, subcycled dynamics
   ! RK2 3 stage tracers, sign-preserving conservative
   rk_stage_user            = 3
   topology                 = "cube"
   ! Finally, set the HOMME variables which have different names
!jt   fine_ne                  = se_fine_ne
   ftype                    = se_ftype
!jt   statediag_numtrac        = MIN(se_statediag_numtrac,pcnst)
!jt   hypervis_power           = se_hypervis_power
   hypervis_scaling         = se_hypervis_scaling
   hypervis_subcycle        = se_hypervis_subcycle
!!$   if (hypervis_subcycle_sponge<0) then
!!$     hypervis_subcycle_sponge = hypervis_subcycle
!!$   else
!!$     hypervis_subcycle_sponge = se_hypervis_subcycle_sponge
!!$   end if
   hypervis_subcycle_q      = se_hypervis_subcycle_q
   limiter_option           = se_limiter_option
!jt   max_hypervis_courant     = se_max_hypervis_courant
!jt   refined_mesh             = se_refined_mesh
   ne                       = se_ne
   nsplit                   = se_nsplit
   nu                       = se_nu
   nu_div                   = se_nu_div
   nu_p                     = se_nu_p
   nu_s                     = se_nu_s
   nu_q                     = se_nu_q !for tracer-wind consistency nu_q must me equal to nu_p
   nu_top                   = se_nu_top
!!$   sponge_del4_nu_fac       = se_sponge_del4_nu_fac
!!$   sponge_del4_nu_div_fac   = se_sponge_del4_nu_div_fac
!!$   sponge_del4_lev          = se_sponge_del4_lev
   qsplit                   = se_qsplit
   rsplit                   = se_rsplit
   statefreq                = se_statefreq
   tstep_type               = se_tstep_type
!!$   vert_remap_uvTq_alg      = set_vert_remap(se_vert_remap_T, se_vert_remap_uvTq_alg)
!!$   vert_remap_tracer_alg    = set_vert_remap(se_vert_remap_T, se_vert_remap_tracer_alg)
!jt   fv_nphys                 = se_fv_nphys
   lcp_moist                = se_lcp_moist
!!$   large_Courant_incr       = se_large_Courant_incr
!!$   fvm_supercycling         = se_fvm_supercycling
!!$   fvm_supercycling_jet     = se_fvm_supercycling_jet
!!$   kmin_jet                 = se_kmin_jet
!!$   kmax_jet                 = se_kmax_jet
!jt   variable_nsplit          = .false.
!jt   phys_dyn_cp              = se_phys_dyn_cp
!jt   molecular_diff           = se_molecular_diff

!!$   if (fv_nphys > 0) then
!!$      ! Use finite volume physics grid and CSLAM for tracer advection
!!$      nphys_pts = fv_nphys*fv_nphys
!!$      qsize = thermodynamic_active_species_num ! number tracers advected by GLL
!!$      ntrac = pcnst                    ! number tracers advected by CSLAM
!!$   else
      ! Use GLL grid for physics and tracer advection
!jt      nphys_pts = npsq
      qsize = pcnst
      ntrac = 0
!!$   end if

   if (rsplit < 1) then
      call endrun('dyn_readnl: rsplit must be > 0')
   end if

   ! if restart or branch run
   if (.not. initial_run) then
      runtype = 1
   end if

   ! HOMME wants 'none' to indicate no mesh file
   if (len_trim(se_mesh_file) == 0) then
      se_mesh_file = 'none'
!!$      if (se_refined_mesh) then
!!$         call endrun('dyn_readnl ERROR: se_refined_mesh=.true. but no se_mesh_file')
!!$      end if
   end if
   call homme_postprocess_namelist(se_mesh_file, par)

  if (masterproc) then
    write(iulog, '(a,i0)') 'dyn_readnl: se_ftype = ',se_ftype
    write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_order = ',se_hypervis_order
    write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_subcycle = ',se_hypervis_subcycle
    write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_subcycle_q = ',se_hypervis_subcycle_q
    write(iulog, '(a,i0)') 'dyn_readnl: se_limiter_option = ',se_limiter_option
!!$    if (.not. se_refined_mesh) then
      write(iulog, '(a,i0)') 'dyn_readnl: se_ne = ',se_ne
!!$    end if
    write(iulog, '(a,i0)') 'dyn_readnl: se_npes = ',se_npes
    write(iulog, '(a,i0)') 'dyn_readnl: se_nsplit = ',se_nsplit
    write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu = ',se_nu
    write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_div = ',se_nu_div
    write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_p = ',se_nu_p
    write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_s = ',se_nu_s
    write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_q = ',se_nu_q
    write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_top = ',se_nu_top
    write(iulog, '(a,i0)') 'dyn_readnl: se_qsplit = ',se_qsplit
    write(iulog, '(a,i0)') 'dyn_readnl: se_rsplit = ',se_rsplit
    write(iulog, '(a,i0)') 'dyn_readnl: se_statefreq = ',se_statefreq
    write(iulog, '(a,i0)') 'dyn_readnl: se_tstep_type = ',se_tstep_type
    write(iulog, '(a,i0)') 'dyn_readnl: se_vert_remap_q_alg = ',se_vert_remap_q_alg
!!$    if (se_refined_mesh) then
!!$      write(iulog, *) 'dyn_readnl: Refined mesh simulation'
!!$      write(iulog, *) 'dyn_readnl: se_mesh_file = ',trim(se_mesh_file)
!jt      if (abs(se_hypervis_power) < 1.0e-12_r8) then
!jt        write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_power = ',se_hypervis_power, ', (tensor hyperviscosity)'
!jt        write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_scaling = ',se_hypervis_scaling
!jt      else if (abs(se_hypervis_power - 3.322_r8) < 1.0e-12_r8) then
!jt        write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_power = ',se_hypervis_power, ', (scalar hyperviscosity)'
!jt        write(iulog, '(a,i0)') 'dyn_readnl: se_fine_ne = ',se_fine_ne
!jt      else
!jt        write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_power = ',se_hypervis_power
!jt        write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_scaling = ',se_hypervis_scaling
!jt        write(iulog, '(a,e11.4)') 'dyn_readnl: se_fine_ne = ',se_fine_ne
!jt      end if
!jt      write(iulog, '(a,e11.4)') 'dyn_readnl: se_max_hypervis_courant = ',se_max_hypervis_courant
!!$    end if


!!XXgoldyXX: v For future physgrid commit
!      write(iulog,*) 'dyn_readnl: fv_nphys = ', fv_nphys, ', nphys_pts = ', nphys_pts
!      if (fv_nphys > 0) then
!        if (trim(write_phys_grid) == 'grid') then
!          write(iulog,*) "dyn_readnl: write physics grid file = ", trim(phys_grid_file)
!        else if (trim(write_phys_grid) == 'interp') then
!          write(iulog,*) "dyn_readnl: write physics interp file = ", trim(phys_grid_file)
!        else
!          write(iulog,*) "dyn_readnl: do not write physics grid or interp file"
!        end if
!      end if
!!XXgoldyXX: ^ For future physgrid commit
 end if

!jt! Create mapping files using SE basis functions if requested
!jt call native_mapping_readnl(NLFileName)

end subroutine dyn_readnl

!=============================================================================================

subroutine dyn_register()

   use physics_buffer,  only: pbuf_add_field, dtype_r8
   use ppgrid,          only: pcols, pver
   use phys_control,    only: use_gw_front, use_gw_front_igw

   ! These fields are computed by the dycore and passed to the physics via the
   ! physics buffer.

   if (use_gw_front .or. use_gw_front_igw) then
      call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), &
         frontgf_idx)
      call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), &
         frontga_idx)
   end if

end subroutine dyn_register

!=============================================================================================

subroutine dyn_init(dyn_in, dyn_out)

    use dyn_grid,         only: elem
    use cam_control_mod,  only: aqua_planet, ideal_phys, adiabatic
    use cam_instance,     only: inst_index
!jt    use native_mapping,   only: create_native_mapping_files
    use cam_pio_utils,    only: clean_iodesc_list

    use prim_driver_mod,  only: prim_init2
    use prim_si_mod,      only: prim_set_mass
    use hybrid_mod,       only: hybrid_create
    use parallel_mod,     only: par
    use control_mod,      only: runtype
!jt    use comsrf,           only: sgh, sgh30
    use element_ops,      only: set_thermostate
!jt    use nctopo_util_mod,  only: nctopo_util_driver

    type (dyn_import_t), intent(out) :: dyn_in
    type (dyn_export_t), intent(out) :: dyn_out

    integer :: ithr, nets, nete, ie, k, tlev
    real(r8), parameter :: Tinit=300.0_r8
    type(hybrid_t) :: hybrid
    real(r8) :: temperature(np,np,nlev),ps(np,np)
   !----------------------------------------------------------------------------

   ! Initialize the import/export objects
   dyn_in%elem  => elem
!jt   dyn_in%fvm   => fvm

   dyn_out%elem => elem
!jt   dyn_out%fvm  => fvm

!jt   ! Create mapping files using SE basis functions if requested
!jt   call create_native_mapping_files(par, elem, 'native')
!jt   call create_native_mapping_files(par, elem, 'bilin')

   call set_phis(dyn_in)

   if (initial_run) then
      call read_inidat(dyn_in)
      call clean_iodesc_list()
   end if

    if(par%dynproc) then

#ifdef HORIZ_OPENMP
       if (iam==0) write (iulog,*) "dyn_init: hthreads=",hthreads,&
                                   "max_threads=",omp_get_max_threads()
       !$OMP PARALLEL NUM_THREADS(hthreads), DEFAULT(SHARED), PRIVATE(ie,ithr,nets,nete,hybrid)
#endif
#ifdef COLUMN_OPENMP
       call omp_set_num_threads(vthreads)
#endif
       ithr=omp_get_thread_num()
       nets=dom_mt(ithr)%start
       nete=dom_mt(ithr)%end
       hybrid = hybrid_create(par,ithr,hthreads)

!jt    moisture='moist'

       if(adiabatic) then
!jt         moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%q(:,:,:,:)=0.0_r8
                elem(ie)%derived%fq(:,:,:,:)=0.0_r8
             end do
          end if
       else if(ideal_phys) then
!jt         moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%ps_v(:,:,:) =hvcoord%ps0

                elem(ie)%state%phis(:,:)=0.0_r8

                elem(ie)%state%v(:,:,:,:,:) =0.0_r8

                elem(ie)%state%q(:,:,:,:)=0.0_r8

                temperature(:,:,:)=300.0_r8
                ps=hvcoord%ps0
                call set_thermostate(elem(ie),ps,temperature,hvcoord)
             end do
          end if
       else if(aqua_planet .and. runtype==0)  then
          do ie=nets,nete
             elem(ie)%state%phis(:,:)=0.0_r8
          end do
!jt          if(allocated(sgh)) sgh=0.0_r8
!jt          if(allocated(sgh30)) sgh30=0.0_r8
       end if

       do ie=nets,nete
          elem(ie)%derived%FM=0.0_r8
          elem(ie)%derived%FT=0.0_r8
          elem(ie)%derived%FQ=0.0_r8
#ifdef MODEL_THETA_L
          elem(ie)%derived%FPHI=0.0_r8
          elem(ie)%derived%FVTheta=0.0_r8
#endif
       end do

       ! scale PS to achieve prescribed dry mass
       if (runtype == 0) then
          ! new run, scale mass to value given in namelist, if needed
          call prim_set_mass(elem, TimeLevel,hybrid,hvcoord,nets,nete)
       endif

       call t_startf('prim_init2')
       call prim_init2(elem,hybrid,nets,nete, TimeLevel, hvcoord)
       call t_stopf('prim_init2')

!jt      ! This subroutine is used to create nc_topo files, if requested
!jt       call nctopo_util_driver(elem,hybrid,nets,nete)

#ifdef HORIZ_OPENMP
       !$OMP END PARALLEL
#endif
    end if

!jt    if (inst_index == 1) then
!jt      call write_grid_mapping(par, elem)
!jt   end if


end subroutine dyn_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  RUN --- Driver for the
  !
  ! !INTERFACE:
  subroutine dyn_run( dyn_state, rc )

    ! !USES:
!jt    use iop_data_mod,     only: single_column, dp_crm, use_3dfrc
!jt    use se_iop_intr_mod,  only: apply_iop_forcing
    use parallel_mod,     only : par
    use prim_driver_mod,  only: prim_run_subcycle
    use dimensions_mod,   only : nlev
    use time_mod,         only: tstep
    use hybrid_mod,       only: hybrid_create
    implicit none


    type (dyn_export_t), intent(inout)       :: dyn_state   !  container
    type(hybrid_t) :: hybrid

    integer, intent(out)               :: rc      ! Return code
    integer ::  n
    integer :: nets, nete, ithr
    integer :: ie
    logical :: single_column_in, do_prim_run

    ! !DESCRIPTION:
    !
    if(par%dynproc) then
#ifdef HORIZ_OPENMP
       !if (iam==0) write (iulog,*) "dyn_run: hthreads=",hthreads,&
       !                            "max_threads=",omp_get_max_threads()
       !$OMP PARALLEL NUM_THREADS(hthreads), DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid,n)
#endif
#ifdef COLUMN_OPENMP
       ! nested threads
       call omp_set_num_threads(vthreads)
#endif
       ithr=omp_get_thread_num()
       nets=dom_mt(ithr)%start
       nete=dom_mt(ithr)%end
       hybrid = hybrid_create(par,ithr,hthreads)

       do_prim_run = .true. ! Always do prim_run_subcycle
                            ! Unless turned off by SCM for specific cases

!jt       single_column_in = single_column

       ! if doubly period CRM mode we want dycore to operate in non-SCM mode,
       !   (typically the single_column flag is true in this case
       !   because we need to take advantage of SCM infrastructure and forcing)
       !   thus turn this switch to false for dycore input.  NOTE that
       !   dycore in SCM mode means that only the large scale vertical
       !   advection is computed (i.e. no horizontal communication)
!jt       if (dp_crm) then
!jt         single_column_in = .false.
!jt       endif

       ! if true SCM mode (NOT DP-CRM mode) do not call
       !   dynamical core if 3D forcing is prescribed
       !   (since large scale vertical advection is accounted for
       !   in that forcing)
!jt       if (single_column .and. .not. dp_crm) then
!jt         if (use_3dfrc) do_prim_run = .false.
!jt       endif

       if (do_prim_run) then
         do n=1,nsplit
           ! forward-in-time RK, with subcycling
           call t_startf('prim_run_subcycle')
           call prim_run_subcycle(dyn_state%elem,hybrid,nets,nete,&
               tstep, single_column_in, TimeLevel, hvcoord, n)
           call t_stopf('prim_run_subcycle')
         end do
       endif

!jt       if (single_column) then
!jt         call apply_iop_forcing(dyn_state%elem,hvcoord,hybrid,TimeLevel,3,.false.,nets,nete)
!jt       endif

#ifdef HORIZ_OPENMP
       !$OMP END PARALLEL
#endif
    end if
    rc = DYN_RUN_SUCCESS

    !EOC
end subroutine dyn_run

!=============================================================================================

subroutine dyn_final(DYN_STATE, RESTART_FILE)

  type (elem_state_t), target     :: DYN_STATE
  character(LEN=*)   , intent(IN) :: RESTART_FILE

end subroutine dyn_final

!=============================================================================================


subroutine read_inidat(dyn_in)

   use shr_vmath_mod,       only: shr_vmath_log
   use hycoef,              only: ps0
   use constituents,            only: cnst_name, cnst_read_iv, qmin,cnst_is_a_water_species
   use cam_control_mod,     only: ideal_phys, aqua_planet
   use cam_initfiles,       only: initial_file_get_id, topo_file_get_id, pertlim
   use cam_history_support, only: max_fieldname_len
   use cam_grid_support,    only: cam_grid_get_local_size, cam_grid_get_gcid
   use const_init,          only: cnst_init_default

    use chemistry,               only: chem_implements_cnst, chem_init_cnst
   use carma_intr,          only: carma_implements_cnst, carma_init_cnst
    use tracers,                 only: tracers_implements_cnst, tracers_init_cnst
    use aoa_tracers,             only: aoa_tracers_implements_cnst, aoa_tracers_init_cnst
    use clubb_intr,              only: clubb_implements_cnst, clubb_init_cnst
   use rk_stratiform,       only: rk_stratiform_implements_cnst, rk_stratiform_init_cnst
    use microp_driver,           only: microp_driver_implements_cnst, microp_driver_init_cnst
    use phys_control,            only: phys_getopts
    use co2_cycle,               only: co2_implements_cnst, co2_init_cnst

   use dof_mod,             only: putUniquePoints
   use edge_mod,            only : edgevpack_nlyr, edgevunpack_nlyr, edge_g
!jt   use nctopo_util_mod,     only: nctopo_util_inidat

!jt   use iop_data_mod,            only: setiopupdate, setiopupdate_init, readiopdata
!jt   use se_iop_intr_mod,         only: iop_setinitial, iop_broadcast
   use element_ops,             only: set_thermostate
   use gllfvremap_mod,          only: gfr_fv_phys_to_dyn_topo

   type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

!jt   type(file_desc_t), pointer :: fh_ini, fh_topo
    type(file_desc_t), pointer :: fh_ini
    real(r8), parameter :: rad2deg = 180.0 / SHR_CONST_PI
    type(element_t), pointer :: elem(:)
    real(r8), allocatable :: tmp(:,:,:)    ! (npsq,nlev,nelemd)
    real(r8), allocatable :: tmp_point(:,:)! (npsq,nlev)
    real(r8), allocatable :: qtmp(:,:,:,:,:)    ! (np,np,nlev,nelemd,n)
    real(r8), allocatable :: dbuf2(:,:)    ! (npsq,nelemd)
    real(r8), allocatable :: dbuf3(:,:,:)  ! (npsq,nlev,nelemd)
    real(r8) :: ps(np,np)
!jt    logical,  allocatable :: tmpmask(:,:)  ! (npsq,nlev,nelemd) unique grid val
    real(r8), allocatable :: phis_tmp(:,:) ! (nphys_sq,nelemd)
    logical,  allocatable :: pmask(:)           ! (npsq*nelemd) unique grid vals
    integer :: nphys_sq                    ! # of fv physics columns per element
    integer :: ie, k, t
!jt    integer :: indx_scm, ie_scm, i_scm, j_scm
    character(len=max_fieldname_len) :: fieldname, fieldname2
    logical                          :: inic_wet           ! true if initial condition is based on
                                                           ! wet pressure and water species
    logical :: found
    integer :: kptr, m_cnst
    real(r8) :: p_ref(nlev)

    integer,parameter :: pcnst = PCNST
    integer(iMap), pointer :: ldof(:) => NULL() ! Basic (2D) grid dof
    integer,       pointer :: gcid(:) => NULL() ! ID based on ldof with no holes

    character(len=max_fieldname_len) :: ncol_name
    character(len=max_fieldname_len) :: grid_name
    real(r8), allocatable            :: latvals(:),latvals_phys(:)
    real(r8), allocatable            :: lonvals(:),lonvals_phys(:)
    real(r8), pointer                :: latvals_deg(:)
    real(r8), pointer                :: lonvals_deg(:)
    logical :: read_pg_grid
    integer :: rndm_seed_sz
    integer :: pio_errtype
    integer, allocatable :: rndm_seed(:)
    real(r8) :: pertval
    integer :: sysclk
    integer :: i, j, indx, tl
    real(r8), parameter :: D0_0 = 0.0_r8
    real(r8), parameter :: D0_5 = 0.5_r8
    real(r8), parameter :: D1_0 = 1.0_r8
    real(r8), parameter :: D2_0 = 2.0_r8
!jt    real(r8) :: scmposlon, minpoint, testlat, testlon, testval
    character*16 :: subname='READ_INIDAT'
    integer :: nlev_tot

    logical :: iop_update_surface

    tl = 1

   fh_ini  => initial_file_get_id()
!jt   fh_topo => topo_file_get_id()

   if(iam < par%nprocs) then
       elem=> dyn_in%elem
    else
       nullify(elem)
    end if

    allocate(tmp(npsq,nlev,nelemd))
    allocate(tmp_point(1,nlev)) ! To find input at a single location
    allocate(qtmp(np,np,nlev,nelemd,pcnst))
!jt    allocate(qtmp(npsq*nelemd,nlev))
    tmp = 0.0_r8
    qtmp = 0.0_r8

    if (fv_nphys>0) then
      nphys_sq = fv_nphys*fv_nphys
      allocate(phis_tmp(nphys_sq,nelemd))
    end if

    if (par%dynproc) then
      if(elem(1)%idxP%NumUniquePts <=0 .or. elem(1)%idxP%NumUniquePts > np*np) then
         write(iulog,*)  elem(1)%idxP%NumUniquePts
         call endrun(trim(subname)//': invalid idxP%NumUniquePts')
      end if
    end if

    ! Set mask to indicate which columns are active
    nullify(ldof)
    call cam_grid_get_gcid(cam_grid_id(ini_grid_name), ldof)
    allocate(pmask(npsq*nelemd))
    pmask(:) = (ldof /= 0)
    ! lat/lon needed in radians
    latvals_deg => cam_grid_get_latvals(cam_grid_id(ini_grid_name))
    lonvals_deg => cam_grid_get_lonvals(cam_grid_id(ini_grid_name))
    allocate(latvals(np*np*nelemd))
    allocate(lonvals(np*np*nelemd))
    latvals(:) = latvals_deg(:)*deg2rad
    lonvals(:) = lonvals_deg(:)*deg2rad

    ! Set PIO to return error codes when reading data from IC file.
    call pio_seterrorhandling(fh_ini, PIO_BCAST_ERROR, pio_errtype)

!jt!   Determine column closest to SCM point
!jt    if (single_column .and. .not. scm_multcols .and. par%dynproc) then
!jt      if (scmlon .lt. 0._r8) then
!jt        scmposlon=scmlon+360._r8
!jt      else
!jt        scmposlon=scmlon
!jt      endif
!jt      minpoint=10000.0_r8
!jt      ie_scm=0
!jt      i_scm=0
!jt      j_scm=0
!jt      indx_scm=0
!jt      do ie=1, nelemd
!jt        indx=1
!jt        do j=1, np
!jt          do i=1, np
!jt            testlat=elem(ie)%spherep(i,j)%lat * rad2deg
!jt            testlon=elem(ie)%spherep(i,j)%lon * rad2deg
!jt            if (testlon .lt. 0._r8) testlon=testlon+360._r8
!jt            testval=abs(scmlat-testlat)+abs(scmposlon-testlon)
!jt            if (testval .lt. minpoint) then
!jt              ie_scm=ie
!jt              indx_scm=indx
!jt              i_scm=i
!jt              j_scm=j
!jt              minpoint=testval
!jt            endif
!jt            indx=indx+1
!jt          enddo
!jt        enddo
!jt      enddo
!jt
!jt      if (ie_scm == 0 .or. i_scm == 0 .or. j_scm == 0 .or. indx_scm == 0) then
!jt        call endrun('Could not find closest SCM point on input datafile')
!jt      endif
!jt
!jt    endif ! single_column

!jt    if (scm_multcols) then
!jt      indx_scm = 1
!jt    endif

      ! Read in 3-D fields

      if (dyn_field_exists(fh_ini, 'U')) then
         call read_dyn_var('U', fh_ini, ini_grid_hdim_name, dbuf3)
      else
         call endrun(trim(subname)//': U not found')
      end if
      do ie = 1, nelemd
         elem(ie)%state%v = 0.0_r8
         indx = 1
         do j = 1, np
            do i = 1, np
               elem(ie)%state%v(i,j,1,:,1) = dbuf3(indx,:,ie)
               indx = indx + 1
            end do
         end do
      end do

      if (dyn_field_exists(fh_ini, 'V')) then
         call read_dyn_var('V', fh_ini, ini_grid_hdim_name, dbuf3)
      else
         call endrun(trim(subname)//': V not found')
      end if
      do ie = 1, nelemd
         indx = 1
         do j = 1, np
            do i = 1, np
               elem(ie)%state%v(i,j,2,:,1) = dbuf3(indx,:,ie)
               indx = indx + 1
            end do
         end do
      end do

      if (dyn_field_exists(fh_ini, 'T')) then
         call read_dyn_var('T', fh_ini, ini_grid_hdim_name, dbuf3)
      else
         call endrun(trim(subname)//': T not found')
      end if
      do ie=1,nelemd
#ifdef MODEL_THETA_L
         elem(ie)%derived%FT=0.0_r8
#else
         elem(ie)%state%T=0.0_r8
#endif
         indx = 1
         do j = 1, np
            do i = 1, np
#ifdef MODEL_THETA_L
               elem(ie)%derived%FT(i,j,:) = dbuf3(indx,:,ie)
#else
               elem(ie)%state%T(i,j,:,1) = dbuf3(indx,:,ie)
#endif
               indx = indx + 1
            end do
         end do
      end do

      if (pertlim .ne. 0.0_r8) then
         if (masterproc) then
            write(iulog,*) trim(subname), ': Adding random perturbation bounded', &
               'by +/- ', pertlim, ' to initial temperature field'
         end if

         call random_seed(size=rndm_seed_sz)
         allocate(rndm_seed(rndm_seed_sz))

         do ie = 1, nelemd
            ! seed random number generator based on element ID
            ! (possibly include a flag to allow clock-based random seeding)
            rndm_seed = elem(ie)%GlobalId
            call random_seed(put=rndm_seed)
            do i = 1, np
               do j = 1, np
                  do k = 1, nlev
                     call random_number(pertval)
                     pertval = 2.0_r8*pertlim*(0.5_r8 - pertval)
#ifdef MODEL_THETA_L
                     elem(ie)%derived%FT(i,j,k) = elem(ie)%derived%FT(i,j,k)*(1.0_r8 + pertval)
#else
                     elem(ie)%state%T(i,j,k,1) = elem(ie)%state%T(i,j,k,1)*(1.0_r8 + pertval)
#endif
                  end do
               end do
            end do
         end do

         deallocate(rndm_seed)
      end if

      ! Cleanup
      deallocate(dbuf2)
      deallocate(dbuf3)

   allocate(dbuf3(npsq,nlev,nelemd))

   do m_cnst = 1, pcnst

      if (analytic_ic_active() .and. cnst_is_a_water_species(cnst_name(m_cnst))) cycle

      found = .false.
      if (cnst_read_iv(m_cnst)) then
         found = dyn_field_exists(fh_ini, trim(cnst_name(m_cnst)), required=.false.)
      end if

      if (found) then
         call read_dyn_var(trim(cnst_name(m_cnst)), fh_ini, ini_grid_hdim_name, dbuf3)
      else
         call cnst_init_default(m_cnst, latvals, lonvals, dbuf3, pmask)
      end if

      do ie = 1, nelemd
         ! Copy tracers defined on GLL grid into Eulerian array
         ! Make sure tracers have at least minimum value
         do k=1, nlev
            indx = 1
            do j = 1, np
               do i = 1, np
                  ! Set qtmp at the unique columns only: zero non-unique columns
                  if (pmask(((ie - 1) * npsq) + indx)) then
                     qtmp(i,j, k, ie, m_cnst) = max(qmin(m_cnst),dbuf3(indx,k,ie))
                  else
                     qtmp(i,j, k, ie, m_cnst) = 0.0_r8
                  end if
                  indx = indx + 1
               end do
            end do
         end do
      end do

   end do ! pcnst

   ! Cleanup
   deallocate(dbuf3)
   ! Put the error handling back the way it was
   call pio_seterrorhandling(fh_ini, pio_errtype)

   ! Cleanup
   deallocate(pmask)
   deallocate(latvals)
   deallocate(lonvals)

   if (associated(ldof)) then
      deallocate(ldof)
      nullify(ldof)
   end if

   ! Read ICs from file.  Assume all fields in the initial file are on the GLL grid.

   allocate(dbuf2(npsq,nelemd))
   allocate(dbuf3(npsq,nlev,nelemd))

   ! Read 2-D field

   fieldname  = 'PS'
   fieldname2 = 'PSDRY'
   if (dyn_field_exists(fh_ini, trim(fieldname), required=.false.)) then
      inic_wet = .true.
      call read_dyn_var(trim(fieldname), fh_ini, ini_grid_hdim_name, dbuf2)
   elseif (dyn_field_exists(fh_ini, trim(fieldname2), required=.false.)) then
      inic_wet = .false.
      call read_dyn_var(trim(fieldname2), fh_ini, ini_grid_hdim_name, dbuf2)
   else
      call endrun(trim(subname)//': PS or PSDRY must be on GLL grid')
   end if

!!$    fieldname = 'PS'
!!$    tmp(:,1,:) = 0.0_r8
!!$    call t_startf('read_inidat_infld')
!!$!jt    if (.not. scm_multcols) then
!!$      call infld(fieldname, fh_ini, ncol_name,      &
!!$           1, npsq, 1, nelemd, tmp(:,1,:), found, gridname=grid_name)
!!$!jt    else
!!$!jt      call infld(fieldname, fh_ini, ncol_name,      &
!!$!jt           1, 1, 1, 1, tmp(:,1,:), found, gridname=grid_name)
!!$!jt    endif
!!$    call t_stopf('read_inidat_infld')
!!$    if(.not. found) then
!!$       call endrun('Could not find PS field on input datafile')
!!$    end if
!!$
!!$    ! Check read-in data to make sure it is in the appropriate units
!!$    allocate(tmpmask(npsq,nelemd))
!!$    tmpmask = (reshape(ldof, (/npsq,nelemd/)) /= 0)
!!$
!!$!jt    if(minval(tmp(:,1,:), mask=tmpmask) < 10000._r8 .and. .not. scm_multcols) then
!!$    if(minval(tmp(:,1,:), mask=tmpmask) < 10000._r8) then
!!$       call endrun('Problem reading ps field')
!!$    end if

!jt    if (scm_multcols) then
!jt      if (tmp(1,1,1) < 10000._r8) then
!jt        call endrun('Problem reading ps field')
!jt      endif
!jt    endif

!jt    deallocate(tmpmask)

#ifndef planet_mars
      if (iam < par%nprocs) then
         if (minval(dbuf2, mask=reshape(pmask, (/npsq,nelemd/))) < 10000._r8) then
            call endrun(trim(subname)//': Problem reading ps or psdry field -- bad values')
         end if
      end if
#endif
      do ie = 1, nelemd
         indx = 1
         do j = 1, np
            do i = 1, np
               elem(ie)%state%ps_v(i,j,tl) = dbuf2(indx,ie) ! can be either wet or dry ps
               indx = indx + 1
            end do
         end do
      end do


!!$      do ie=1,nelemd
!!$       elem(ie)%state%ps_v=0.0_r8
!!$          indx = 1
!!$          do j = 1, np
!!$             do i = 1, np
!!$                elem(ie)%state%ps_v(i,j,tl) = tmp(indx,1,ie)
!!$!jt                if (single_column .and. .not. scm_multcols) elem(ie)%state%ps_v(i,j,tl) = tmp(indx_scm,1,ie_scm)
!!$!jt                if (scm_multcols) elem(ie)%state%ps_v(i,j,tl) = tmp(1,1,1)
!!$                indx = indx + 1
!!$             end do
!!$          end do
!!$    end do

!!$    read_pg_grid = .false.
!!$    if ( (ideal_phys .or. aqua_planet)) then
!!$       tmp(:,1,:) = 0._r8
!!$       if (fv_nphys > 0) phis_tmp(:,:) = 0._r8
!!$    else
!!$      fieldname = 'PHIS'
!!$      tmp(:,1,:) = 0.0_r8
!!$      if (fv_nphys == 0) then
!!$         call t_startf('read_inidat_infld')
!!$!jt         if (.not. scm_multcols) then
!!$           call infld(fieldname, ncid_topo, ncol_name,      &
!!$              1, npsq, 1, nelemd, tmp(:,1,:), found, gridname=grid_name)
!!$!jt         else
!!$!jt           call infld(fieldname, ncid_topo, ncol_name,      &
!!$!jt              1, 1, 1, 1, tmp(:,1,:), found, gridname=grid_name)
!!$!jt         endif
!!$         call t_stopf('read_inidat_infld')
!!$
!!$      else
!!$         ! Attempt to read a mixed GLL-FV topo file, which contains PHIS_d in
!!$         ! addition to PHIS.
!!$
!!$         call t_startf('read_inidat_infld')
!!$         call infld(trim(fieldname) // '_d', ncid_topo, ncol_name, &
!!$              1, npsq, 1, nelemd, tmp(:,1,:), found, gridname=grid_name)
!!$         call t_stopf('read_inidat_infld')
!!$
!!$         if (found) then
!!$            if (masterproc) then
!!$               write(iulog,*) 'reading GLL ', trim(fieldname) // '_d', &
!!$                    ' on gridname ', trim(grid_name)
!!$            end if
!!$         else
!!$            ! Pure-FV topo file, so read FV PHIS and map it to GLL.
!!$            if (masterproc) then
!!$               write(iulog,*) 'reading FV ', trim(fieldname), &
!!$                    ' on gridname physgrid_d'
!!$            end if
!!$            read_pg_grid = .true.
!!$
!!$            call t_startf('read_inidat_infld')
!!$            call infld(fieldname, ncid_topo, ncol_name, 1, nphys_sq, &
!!$                 1, nelemd, phis_tmp, found, gridname='physgrid_d')
!!$            call t_stopf('read_inidat_infld')
!!$
!!$            call gfr_fv_phys_to_dyn_topo(par, dom_mt, elem, phis_tmp)
!!$         end if
!!$      end if
!!$      if(.not. found) then
!!$         call endrun('Could not find PHIS field on input datafile')
!!$      end if
!!$    end if
!!$
!!$    if (.not. read_pg_grid) then
!!$      do ie=1,nelemd
!!$         elem(ie)%state%phis=0.0_r8
!!$         indx = 1
!!$         do j = 1, np
!!$            do i = 1, np
!!$               elem(ie)%state%phis(i,j) = tmp(indx,1,ie)
!!$!jt               if (single_column .and. .not. scm_multcols) elem(ie)%state%phis(i,j) = tmp(indx_scm,1,ie_scm)
!!$!jt               if (scm_multcols) elem(ie)%state%phis(i,j) = tmp(1,1,1)
!!$               indx = indx + 1
!!$            end do
!!$         end do
!!$      end do
!!$    end if ! not read_pg_grid

!jt    if (single_column) then
!jt      iop_update_surface = .false.
!jt      if (masterproc) call setiopupdate_init()
!jt      if (masterproc) call readiopdata(iop_update_surface,hyam,hybm)
!jt      if (scm_multcols) call iop_broadcast()
!jt      call iop_setinitial(elem)
!jt    endif

!!$      if (dp_crm) then
!!$        ! Define reference pressure, to potentially restrict initial perturbations
!!$        !  to a certain height if requested
!!$        do k=1,nlev
!!$          p_ref(k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*hvcoord%ps0
!!$        enddo
!!$      endif

!jt    if (.not. single_column) then
    if (.true.) then

      ! once we've read all the fields we do a boundary exchange to
      ! update the redundent columns in the dynamics
!jt      nlev_tot=(3+pcnst)*nlev+2
      nlev_tot=(3+pcnst)*nlev+1

#ifdef MODEL_THETA_L
      do ie=1,nelemd
        kptr=0
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%ps_v(:,:,tl),1,kptr,nlev_tot)
!!$        kptr=kptr+1
!!$        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%v(:,:,:,:,tl),2*nlev,kptr,nlev_tot)
        kptr=kptr+2*nlev
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FT(:,:,:),nlev,kptr,nlev_tot)
        kptr=kptr+nlev
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
      end do
#else
      do ie=1,nelemd
        kptr=0
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%ps_v(:,:,tl),1,kptr,nlev_tot)
!!$        kptr=kptr+1
!!$        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%v(:,:,:,:,tl),2*nlev,kptr,nlev_tot)
        kptr=kptr+2*nlev
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%T(:,:,:,tl),nlev,kptr,nlev_tot)
        kptr=kptr+nlev
        call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
      end do
#endif
      if(par%dynproc) then
        call bndry_exchangeV(par,edge_g)
      end if
#ifdef MODEL_THETA_L
      do ie=1,nelemd
        kptr=0
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%ps_v(:,:,tl),1,kptr,nlev_tot)
!!$        kptr=kptr+1
!!$        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%v(:,:,:,:,tl),2*nlev,kptr,nlev_tot)
        kptr=kptr+2*nlev
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FT(:,:,:),nlev,kptr,nlev_tot)
        kptr=kptr+nlev
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
      end do
#else
      do ie=1,nelemd
        kptr=0
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%ps_v(:,:,tl),1,kptr,nlev_tot)
!!$        kptr=kptr+1
!!$        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
        kptr=kptr+1
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%v(:,:,:,:,tl),2*nlev,kptr,nlev_tot)
        kptr=kptr+2*nlev
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%T(:,:,:,tl),nlev,kptr,nlev_tot)
        kptr=kptr+nlev
        call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
      end do
#endif
    endif

!$omp parallel do private(ie, ps, t, m_cnst)
    do ie=1,nelemd
       ps=elem(ie)%state%ps_v(:,:,tl)
#ifdef MODEL_THETA_L
       elem(ie)%state%w_i = 0.0
       call set_thermostate(elem(ie),ps,elem(ie)%derived%FT,hvcoord,elem(ie)%state%Q(:,:,:,1))
       !FT used as tmp array - reset
       elem(ie)%derived%FT = 0.0
#else
       call set_thermostate(elem(ie),ps,elem(ie)%state%T(:,:,:,tl),hvcoord)
#endif
    end do

   ! Cleanup
    deallocate(tmp)
    deallocate(tmp_point)
    deallocate(qtmp)
    if (fv_nphys>0) then
      deallocate(phis_tmp)
    end if

  end subroutine read_inidat

!=============================================================================================
!========================================================================================

subroutine set_phis(dyn_in)

   ! Set PHIS according to the following rules.
   !
   ! 1) If a topo file is specified use it.  This option has highest precedence.
   ! 2) If not using topo file, but analytic_ic option is on, use analytic phis.
   ! 3) Set phis = 0.0.
   !
   ! If using the physics grid then the topo file will be on that grid since its
   ! contents are primarily for the physics parameterizations, and the values of
   ! PHIS should be consistent with the values of sub-grid variability (e.g., SGH)
   ! which are computed on the physics grid.  In this case phis on the physics grid
   ! will be interpolated to the GLL grid.

  use dyn_tests_utils,        only: vcoord=>vc_dry_pressure

   ! Arguments
   type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

   ! local variables
   type(file_desc_t), pointer       :: fh_topo

   type(element_t), pointer         :: elem(:)

   real(r8), allocatable            :: phis_tmp(:,:)      ! (npsp,nelemd)
   real(r8), allocatable            :: phis_phys_tmp(:,:) ! (fv_nphys**2,nelemd)

   integer                          :: i, ie, indx, j, kptr
   integer                          :: ierr, pio_errtype

   character(len=max_fieldname_len) :: fieldname
   character(len=max_hcoordname_len):: grid_name
   integer                          :: dims(2)
   integer                          :: dyn_cols
   integer                          :: ncol_did
   integer                          :: ncol_size

   integer(iMap), pointer           :: ldof(:)            ! Basic (2D) grid dof
   logical,  allocatable            :: pmask(:)           ! (npsq*nelemd) unique columns

   ! Variables for analytic initial conditions
   integer,  allocatable            :: glob_ind(:)
   logical,  allocatable            :: pmask_phys(:)
   real(r8), pointer                :: latvals_deg(:)
   real(r8), pointer                :: lonvals_deg(:)
   real(r8), allocatable            :: latvals(:)
   real(r8), allocatable            :: lonvals(:)
   real(r8), allocatable            :: latvals_phys(:)
   real(r8), allocatable            :: lonvals_phys(:)
   integer                          :: nlev_tot

   character(len=*), parameter      :: sub='set_phis'
   !----------------------------------------------------------------------------

   fh_topo => topo_file_get_id()

   if (iam < par%nprocs) then
      elem => dyn_in%elem
   else
      nullify(elem)
   end if

   allocate(phis_tmp(npsq,nelemd))
   phis_tmp = 0.0_r8

!!$   if (fv_nphys > 0) then
!!$      allocate(phis_phys_tmp(fv_nphys**2,nelemd))
!!$      phis_phys_tmp = 0.0_r8
!!$      do ie=1,nelemd
!!$        elem(ie)%sub_elem_mass_flux=0.0_r8
!!$#ifdef waccm_debug
!!$        dyn_in%fvm(ie)%CSLAM_gamma = 0.0_r8
!!$#endif
!!$      end do
!!$   end if

   ! Set mask to indicate which columns are active in GLL grid.
   nullify(ldof)
   call cam_grid_get_gcid(cam_grid_id('GLL'), ldof)
   allocate(pmask(npsq*nelemd))
   pmask(:) = (ldof /= 0)
   deallocate(ldof)

   if (associated(fh_topo)) then

      ! Set PIO to return error flags.
      call pio_seterrorhandling(fh_topo, PIO_BCAST_ERROR, pio_errtype)

      ! Set name of grid object which will be used to read data from file
      ! into internal data structure via PIO.
      if (fv_nphys == 0) then
         grid_name = 'GLL'
      else
         grid_name = 'physgrid_d'
      end if

      ! Get number of global columns from the grid object and check that
      ! it matches the file data.
      call cam_grid_dimensions(grid_name, dims)
      dyn_cols = dims(1)

      ! The dimension of the unstructured grid in the TOPO file is 'ncol'.
      ierr = pio_inq_dimid(fh_topo, 'ncol', ncol_did)
      if (ierr /= PIO_NOERR) then
         call endrun(sub//': dimension ncol not found in bnd_topo file')
      end if
      ierr = pio_inq_dimlen(fh_topo, ncol_did, ncol_size)
      if (ncol_size /= dyn_cols) then
         if (masterproc) then
            write(iulog,*) sub//': ncol_size=', ncol_size, ' : dyn_cols=', dyn_cols
         end if
         call endrun(sub//': ncol size in bnd_topo file does not match grid definition')
      end if

      fieldname = 'PHIS'
      if (dyn_field_exists(fh_topo, trim(fieldname))) then
         if (fv_nphys == 0) then
           call read_dyn_var(fieldname, fh_topo, 'ncol', phis_tmp)
         else
            call endrun(sub//': fv_nphys > 0 not implemented for nh dycore')
!!$           call read_phys_field_2d(fieldname, fh_topo, 'ncol', phis_phys_tmp)
!!$           call map_phis_from_physgrid_to_gll(dyn_in%fvm, elem, phis_phys_tmp, &
!!$                phis_tmp, pmask)
         end if
      else
         call endrun(sub//': Could not find PHIS field on input datafile')
      end if

      ! Put the error handling back the way it was
      call pio_seterrorhandling(fh_topo, pio_errtype)

   else if (analytic_ic_active() .and. (iam < par%nprocs)) then

      ! lat/lon needed in radians
      latvals_deg => cam_grid_get_latvals(cam_grid_id('GLL'))
      lonvals_deg => cam_grid_get_lonvals(cam_grid_id('GLL'))
      allocate(latvals(np*np*nelemd))
      allocate(lonvals(np*np*nelemd))
      latvals(:) = latvals_deg(:)*deg2rad
      lonvals(:) = lonvals_deg(:)*deg2rad

      allocate(glob_ind(npsq*nelemd))
      j = 1
      do ie = 1, nelemd
         do i = 1, npsq
            ! Create a global(ish) column index
            glob_ind(j) = elem(ie)%GlobalId
            j = j + 1
         end do
      end do
      call analytic_ic_set_ic(vcoord, latvals, lonvals, glob_ind, &
                              PHIS_OUT=phis_tmp, mask=pmask(:))
      deallocate(glob_ind)

!!$      if (fv_nphys > 0) then
!!$
!!$         ! initialize PHIS on physgrid
!!$         allocate(latvals_phys(fv_nphys*fv_nphys*nelemd))
!!$         allocate(lonvals_phys(fv_nphys*fv_nphys*nelemd))
!!$         indx = 1
!!$         do ie = 1, nelemd
!!$            do j = 1, fv_nphys
!!$               do i = 1, fv_nphys
!!$                  latvals_phys(indx) = dyn_in%fvm(ie)%center_cart_physgrid(i,j)%lat
!!$                  lonvals_phys(indx) = dyn_in%fvm(ie)%center_cart_physgrid(i,j)%lon
!!$                  indx = indx + 1
!!$               end do
!!$            end do
!!$         end do
!!$
!!$         allocate(pmask_phys(fv_nphys*fv_nphys*nelemd))
!!$         pmask_phys(:) = .true.
!!$         allocate(glob_ind(fv_nphys*fv_nphys*nelemd))
!!$
!!$         j = 1
!!$         do ie = 1, nelemd
!!$            do i = 1, fv_nphys*fv_nphys
!!$               ! Create a global(ish) column index
!!$               glob_ind(j) = elem(ie)%GlobalId
!!$               j = j + 1
!!$            end do
!!$         end do
!!$
!!$         call analytic_ic_set_ic(vcoord, latvals_phys, lonvals_phys, glob_ind, &
!!$                                 PHIS_OUT=phis_phys_tmp, mask=pmask_phys)
!!$
!!$         deallocate(latvals_phys)
!!$         deallocate(lonvals_phys)
!!$         deallocate(pmask_phys)
!!$         deallocate(glob_ind)
!!$      end if
!!$
   end if

   deallocate(pmask)

   ! Set PHIS in element objects
   do ie = 1, nelemd
      elem(ie)%state%phis = 0.0_r8
      indx = 1
      do j = 1, np
         do i = 1, np
            elem(ie)%state%phis(i,j) = phis_tmp(indx, ie)
            indx = indx + 1
         end do
      end do
   end do
!!$   if (fv_nphys > 0) then
!!$      do ie = 1, nelemd
!!$         dyn_in%fvm(ie)%phis_physgrid = RESHAPE(phis_phys_tmp(:,ie),(/fv_nphys,fv_nphys/))
!!$      end do
!!$   end if

   deallocate(phis_tmp)
!!$   if (fv_nphys > 0) then
!!$      deallocate(phis_phys_tmp)
!!$   end if
   nlev_tot=1
   ! boundary exchange to update the redundent columns in the element objects
   do ie = 1, nelemd
      kptr = 0
      call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
   end do
   if(iam < par%nprocs) then
      call bndry_exchangeV(par,edge_g)
   end if
   do ie = 1, nelemd
      kptr = 0
      call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis,1,kptr,nlev_tot)
   end do

end subroutine set_phis

!========================================================================================

subroutine check_file_layout(file, elem, dyn_cols, file_desc, dyn_ok)

   ! This routine is only called when data will be read from the initial file.  It is not
   ! called when the initial file is only supplying vertical coordinate info.

   type(file_desc_t), pointer       :: file
   type(element_t),   pointer       :: elem(:)
   integer,           intent(in)    :: dyn_cols
   character(len=*),  intent(in)    :: file_desc
   logical,           intent(in)    :: dyn_ok ! .true. iff ncol_d is okay

   integer                          :: ncol_did, ncol_size
   integer                          :: ierr
   integer                          :: ie, i, j
   integer                          :: grid_id
   integer                          :: indx
   real(r8)                         :: dbuf2(npsq, nelemd)
   logical                          :: found
   character(len=max_fieldname_len) :: coordname

   character(len=*), parameter      :: sub = 'check_file_layout'
   !----------------------------------------------------------------------------

   ! Check that number of columns in IC file matches grid definition.
   if (trim(ini_grid_hdim_name) == 'none') then
      call endrun(sub//': ERROR: no horizontal dimension in initial data file. &
         &Cannot read data from file')
   end if

   ierr = pio_inq_dimid(file, trim(ini_grid_hdim_name), ncol_did)
   if (ierr /= PIO_NOERR) then
      call endrun(sub//': ERROR: '//trim(ini_grid_hdim_name)//' dimension not found in ' &
         //trim(file_desc)//' file')
   end if

   ierr = pio_inq_dimlen(file, ncol_did, ncol_size)
   if (ncol_size /= dyn_cols) then
      if (masterproc) then
         write(iulog, '(a,2(a,i0))') trim(sub), ': ncol_size=', ncol_size, &
             ' : dyn_cols=', dyn_cols
      end if
      call endrun(sub//': ERROR: dimension '//trim(ini_grid_hdim_name)//' size not same as in ncdata file')
   end if

   ! Set coordinate name associated with ini_grid_hdim_name.
   if (trim(ini_grid_hdim_name) == 'ncol') then
      coordname = 'lat'
   else
      coordname = 'lat_d'
   end if

   !! Check to make sure file is in correct order
   call read_dyn_var(coordname, file, ini_grid_hdim_name, dbuf2)
   found = .true.
   do ie = 1, nelemd
      indx = 1
      do j = 1, np
         do i = 1, np
            if (abs(dbuf2(indx,ie)) > 1.e-12_r8) then
               if (abs((elem(ie)%spherep(i,j)%lat*rad2deg - dbuf2(indx,ie)) / &
                    dbuf2(indx,ie)) > 1.0e-10_r8) then
                  write(iulog, '(2a,4(i0,a),f11.5,a,f11.5)')                  &
                       "ncdata file latitudes not in correct column order",   &
                       ' on task ', iam, ': elem(', ie, ')%spherep(', i,      &
                       ', ', j, ')%lat = ', elem(ie)%spherep(i,j)%lat,        &
                       ' /= ', dbuf2(indx, ie)*deg2rad
                  call shr_sys_flush(iulog)
                  found = .false.
               end if
            end if
            indx = indx + 1
         end do
      end do
   end do
   if (.not. found) then
      call endrun("ncdata file latitudes not in correct column order")
   end if

   if (trim(ini_grid_hdim_name) == 'ncol') then
      coordname = 'lon'
   else
      coordname = 'lon_d'
   end if

   call read_dyn_var(coordname, file, ini_grid_hdim_name, dbuf2)
   do ie = 1, nelemd
      indx = 1
      do j = 1, np
         do i = 1, np
            if (abs(dbuf2(indx,ie)) > 1.e-12_r8) then
               if (abs((elem(ie)%spherep(i,j)%lon*rad2deg - dbuf2(indx,ie)) / &
                    dbuf2(indx,ie)) > 1.0e-10_r8) then
                  write(iulog, '(2a,4(i0,a),f11.5,a,f11.5)')                  &
                       "ncdata file longitudes not in correct column order",  &
                       ' on task ', iam, ': elem(', ie, ')%spherep(', i,      &
                       ', ', j, ')%lon = ', elem(ie)%spherep(i,j)%lon,        &
                       ' /= ', dbuf2(indx, ie)*deg2rad
                  call shr_sys_flush(iulog)
                  found = .false.
               end if
            end if
            indx = indx + 1
         end do
      end do
   end do
   if (.not. found) then
      call endrun("ncdata file longitudes not in correct column order")
   end if
end subroutine check_file_layout

!========================================================================================

logical function dyn_field_exists(fh, fieldname, required)


   type(file_desc_t), intent(in) :: fh
   character(len=*),  intent(in) :: fieldname
   logical, optional, intent(in) :: required

   ! Local variables
   logical                  :: found
   logical                  :: field_required
   integer                  :: ret
   type(var_desc_t)         :: varid
   character(len=128)       :: errormsg
   !--------------------------------------------------------------------------

   if (present(required)) then
      field_required = required
   else
      field_required = .true.
   end if

   ret = PIO_inq_varid(fh, trim(fieldname), varid)
   found = (ret == PIO_NOERR)
   if (.not. found) then
      if (field_required) then
         write(errormsg, *) trim(fieldname),' was not present in the input file.'
         call endrun('DYN_FIELD_EXISTS: '//errormsg)
      end if
   end if

   dyn_field_exists = found

end function dyn_field_exists

!========================================================================================

subroutine read_dyn_field_2d(fieldname, fh, dimname, buffer)

   ! Dummy arguments
   character(len=*),  intent(in)    :: fieldname
   type(file_desc_t), intent(inout) :: fh
   character(len=*),  intent(in)    :: dimname
   real(r8),          intent(inout) :: buffer(:, :)

   ! Local variables
   logical                  :: found
   real(r8)                 :: fillvalue
   !----------------------------------------------------------------------------

   buffer = 0.0_r8
   call infld(trim(fieldname), fh, dimname, 1, npsq, 1, nelemd, buffer,    &
        found, gridname=ini_grid_name, fillvalue=fillvalue)
   if(.not. found) then
      call endrun('READ_DYN_FIELD_2D: Could not find '//trim(fieldname)//' field on input datafile')
   end if

   ! This code allows use of compiler option to set uninitialized values
   ! to NaN.  In that case infld can return NaNs where the element GLL points
   ! are not "unique columns"
   ! Set NaNs or fillvalue points to zero
   where (isnan(buffer) .or. (buffer==fillvalue)) buffer = 0.0_r8

end subroutine read_dyn_field_2d

!========================================================================================

subroutine read_dyn_field_3d(fieldname, fh, dimname, buffer)

   ! Dummy arguments
   character(len=*),  intent(in)    :: fieldname
   type(file_desc_t), intent(inout) :: fh
   character(len=*),  intent(in)    :: dimname
   real(r8),          intent(inout) :: buffer(:,:,:)

   ! Local variables
   logical                  :: found
   real(r8)                 :: fillvalue
   !----------------------------------------------------------------------------

   buffer = 0.0_r8
   call infld(trim(fieldname), fh, dimname, 'lev',  1, npsq, 1, nlev,         &
        1, nelemd, buffer, found, gridname=ini_grid_name, fillvalue=fillvalue)
   if(.not. found) then
      call endrun('READ_DYN_FIELD_3D: Could not find '//trim(fieldname)//' field on input datafile')
   end if

   ! This code allows use of compiler option to set uninitialized values
   ! to NaN.  In that case infld can return NaNs where the element GLL points
   ! are not "unique columns"
   ! Set NaNs or fillvalue points to zero
   where (isnan(buffer) .or. (buffer == fillvalue)) buffer = 0.0_r8

end subroutine read_dyn_field_3d

!========================================================================================

subroutine read_phys_field_2d(fieldname, fh, dimname, buffer)

   ! Dummy arguments
   character(len=*),  intent(in)    :: fieldname
   type(file_desc_t), intent(inout) :: fh
   character(len=*),  intent(in)    :: dimname
   real(r8),          intent(inout) :: buffer(:, :)

   ! Local variables
   logical                  :: found
   !----------------------------------------------------------------------------

   call infld(trim(fieldname), fh, dimname, 1, fv_nphys**2, 1, nelemd, buffer,    &
      found, gridname='physgrid_d')
   if(.not. found) then
      call endrun('READ_PHYS_FIELD_2D: Could not find '//trim(fieldname)//' field on input datafile')
   end if

end subroutine read_phys_field_2d

!========================================================================================


end module dyn_comp

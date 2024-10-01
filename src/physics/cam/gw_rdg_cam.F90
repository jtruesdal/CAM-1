module gw_rdg_cam

!
! This module handles gravity waves from orographic sources, and was
! extracted from gw_drag in May 2013.
!
use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,only: masterproc
use cam_abortutils, only: endrun


implicit none
private
save

! Public interface
public :: gw_rdg_cam_readnl

!==========================================================================
contains
!==========================================================================

subroutine gw_rdg_cam_readnl(nlfile)
  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8, mpi_logical

  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'gw_rdg_readnl'

  logical  :: gw_rdg_do_divstream, gw_rdg_do_smooth_regimes, gw_rdg_do_adjust_tauoro, &
              gw_rdg_do_backward_compat

  logical  :: gw_rdg_do_vdiff=.true.

  real(r8) :: gw_rdg_C_BetaMax_DS, gw_rdg_C_GammaMax, &
              gw_rdg_Frx0, gw_rdg_Frx1, gw_rdg_C_BetaMax_SM, gw_rdg_Fr_c, &
              gw_rdg_orohmin, gw_rdg_orovmin, gw_rdg_orostratmin, gw_rdg_orom2min

  namelist /gw_rdg_nl/ gw_rdg_do_divstream, gw_rdg_C_BetaMax_DS, gw_rdg_C_GammaMax, &
                       gw_rdg_Frx0, gw_rdg_Frx1, gw_rdg_C_BetaMax_SM, gw_rdg_Fr_c, &
                       gw_rdg_do_smooth_regimes, gw_rdg_do_adjust_tauoro, &
                       gw_rdg_do_backward_compat, gw_rdg_orohmin, gw_rdg_orovmin, &
                       gw_rdg_orostratmin, gw_rdg_orom2min, gw_rdg_do_vdiff

  !----------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'gw_rdg_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, gw_rdg_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(sub // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)
  end if

  ! Broadcast the local variables

  call mpi_bcast(gw_rdg_do_divstream, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_divstream")
  call mpi_bcast(gw_rdg_do_smooth_regimes, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_smooth_regimes")
  call mpi_bcast(gw_rdg_do_adjust_tauoro, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_adjust_tauoro")
  call mpi_bcast(gw_rdg_do_backward_compat, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_backward_compat")

  call mpi_bcast(gw_rdg_C_BetaMax_DS, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_C_BetaMax_DS")
  call mpi_bcast(gw_rdg_C_GammaMax, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_C_GammaMax")
  call mpi_bcast(gw_rdg_Frx0, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_Frx0")
  call mpi_bcast(gw_rdg_Frx1, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_Frx1")
  call mpi_bcast(gw_rdg_C_BetaMax_SM, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_C_BetaMax_SM")
  call mpi_bcast(gw_rdg_Fr_c, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_Fr_c")
  call mpi_bcast(gw_rdg_orohmin, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_orohmin")
  call mpi_bcast(gw_rdg_orovmin, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_orovmin")
  call mpi_bcast(gw_rdg_orostratmin, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_orostratmin")
  call mpi_bcast(gw_rdg_orom2min, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_orom2min")

  call mpi_bcast(gw_rdg_do_vdiff, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_vdiff")

  if (gw_rdg_Fr_c > 1.0_r8) call endrun(sub//": FATAL: gw_rdg_Fr_c must be <= 1")

  call gw_rdg_init(gw_rdg_do_divstream, gw_rdg_C_BetaMax_DS, gw_rdg_C_GammaMax, gw_rdg_Frx0, &
       gw_rdg_Frx1, gw_rdg_C_BetaMax_SM, gw_rdg_Fr_c, gw_rdg_do_smooth_regimes, gw_rdg_do_adjust_tauoro, &
       gw_rdg_do_backward_compat, gw_rdg_orohmin, gw_rdg_orovmin, gw_rdg_orostratmin, gw_rdg_orom2min)

end subroutine gw_rdg_cam_readnl

end module gw_rdg_cam

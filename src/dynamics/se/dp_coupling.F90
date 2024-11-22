!
!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------
module dp_coupling
  use air_composition, only:  rairv
  use bndry_mod,      only: bndry_exchangeV
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use constituents,   only: pcnst, cnst_name
  use cam_history,    only: outfld, write_inithist, hist_fld_active
  use dimensions_mod, only: np, npsq, nelemd, nlev, fv_nphys
  use dof_mod,        only: UniquePoints, PutUniquePoints
  use dyn_comp,       only: dyn_export_t, dyn_import_t
  use dyn_grid,       only: TimeLevel, hvcoord, dom_mt
  use element_ops,    only: get_temperature
  use element_mod,    only: element_t
  use kinds,          only: real_kind, int_kind
  use physics_types,  only: physics_state, physics_tend
  use phys_grid,      only: get_ncols_p, get_gcol_all_p
  use phys_grid,      only: get_dyn_col_p, columns_on_task, get_chunk_info_p
  use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
  use perf_mod,       only: t_startf, t_stopf, t_barrierf
  use parallel_mod,   only: par
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
  use spmd_utils,     only: mpicom, iam
  use qneg_module,    only: qneg3
  use thread_mod,     only: max_num_threads
!jt  use iop_data_mod,   only: single_column
  private
  public :: d_p_coupling, p_d_coupling
!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
  subroutine d_p_coupling(phys_state, phys_tend,  pbuf2d, dyn_out)
    use dyn_comp,              only: frontgf_idx, frontga_idx
    use gllfvremap_mod,        only: gfr_dyn_to_fv_phys
    use gravity_waves_sources, only: gws_src_fnct
    use phys_control,          only: use_gw_front, use_gw_front_igw
    use physics_buffer,        only: physics_buffer_desc, pbuf_get_chunk, &
                                     pbuf_get_field
    use shr_vmath_mod,         only: shr_vmath_exp
    use time_manager,          only: is_first_step
    implicit none
!-----------------------------------------------------------------------
! !INPUT PARAMETERS:
!
    type(dyn_export_t), intent(inout) :: dyn_out    ! dynamics export
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
! !OUTPUT PARAMETERS:

    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend


! LOCAL VARIABLES
    type(element_t), pointer  :: elem(:)                       ! pointer to dyn_out element array
    integer                   :: ie                            ! indices over elements
    integer                   :: lchnk, icol, ilyr             ! indices over chunks, columns, layers
    real(r8)                  :: ps_tmp(npsq,nelemd)           ! temporary array to hold ps
    real(r8)                  :: phis_tmp(npsq,nelemd)         ! temporary array to hold phis
    real(r8)                  :: T_tmp(npsq,pver,nelemd)       ! temporary array to hold T
    real(r8)                  :: uv_tmp(npsq,2,pver,nelemd)    ! temporary array to hold u and v
    real(r8)                  :: q_tmp(npsq,pver,pcnst,nelemd) ! temporary to hold advected constituents
    real(r8)                  :: omega_tmp(npsq,pver,nelemd)   ! temporary array to hold omega

    ! Frontogenesis
    real(r8), allocatable     :: frontgf(:,:,:) ! temporary arrays to hold frontogenesis
    real(r8), allocatable     :: frontga(:,:,:) !   function (frontgf) and angle (frontga)
    ! Pointers to pbuf
    real(r8), pointer         :: pbuf_frontgf(:,:)
    real(r8), pointer         :: pbuf_frontga(:,:)

    integer                   :: ncols,i,j,ierr
    integer                   :: col_ind        ! index over columns
    integer                   :: blk_ind(1)     ! element offset

    integer                   ::  m
    integer                   :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)
    integer                   :: nphys, nphys_sq               ! physics grid parameters
    real (r8)                 :: temperature(np,np,nlev)       ! Temperature from dynamics

    integer :: tl_f, ncol_d

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    !----------------------------------------------------------------------

    if (.not. local_dp_map) then
       call endrun('d_p_coupling: Weak scaling does not support load balancing')
    end if

    nullify(pbuf_chnk)
    nullify(pbuf_frontgf)
    nullify(pbuf_frontga)

    if (fv_nphys > 0) then
      nphys = fv_nphys
    else
      nphys = np
    end if
    nphys_sq = nphys*nphys

    if (use_gw_front .or. use_gw_front_igw) then

       allocate(frontgf(nphys_sq,pver,nelemd), stat=ierr)
       if (ierr /= 0) call endrun("dp_coupling: Allocate of frontgf failed.")

       allocate(frontga(nphys_sq,pver,nelemd), stat=ierr)
       if (ierr /= 0) call endrun("dp_coupling: Allocate of frontga failed.")

    end if

    if( par%dynproc) then

       elem => dyn_out%elem

       tl_f = TimeLevel%n0  ! time split physics (with forward-in-time RK)

      if (use_gw_front .or. use_gw_front_igw) call gws_src_fnct(elem, tl_f, nphys, frontgf, frontga)

      if (fv_nphys > 0) then
        !-----------------------------------------------------------------------
        ! Map dynamics state to FV physics grid
        !-----------------------------------------------------------------------
        call t_startf('dyn_to_fv_phys')
        call gfr_dyn_to_fv_phys(par, dom_mt, tl_f, hvcoord, elem, ps_tmp, phis_tmp, &
             T_tmp, uv_tmp, omega_tmp, q_tmp)
        call t_stopf('dyn_to_fv_phys')

        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
      else ! fv_nphys > 0
        !-----------------------------------------------------------------------
        ! Physics on GLL grid: collect unique points before copying
        !-----------------------------------------------------------------------
       call t_startf('UniquePoints')
       do ie=1,nelemd
          ncols = elem(ie)%idxP%NumUniquePts
          call get_temperature(elem(ie),temperature,hvcoord,tl_f)
          call UniquePoints(elem(ie)%idxP, elem(ie)%state%ps_v(:,:,tl_f), ps_tmp(1:ncols,ie))
          call UniquePoints(elem(ie)%idxP, elem(ie)%state%phis, phis_tmp(1:ncols,ie))
          call UniquePoints(elem(ie)%idxP,  nlev, temperature,                  T_tmp(1:ncols,:,ie))
          call UniquePoints(elem(ie)%idxP,  nlev,elem(ie)%derived%omega_p,      omega_tmp(1:ncols,:,ie))
          call UniquePoints(elem(ie)%idxP,2,nlev,elem(ie)%state%V(:,:,:,:,tl_f),uv_tmp(1:ncols,:,:,ie))
          call UniquePoints(elem(ie)%idxP, nlev,pcnst, elem(ie)%state%Q(:,:,:,:), Q_tmp(1:ncols,:,:,ie))
       end do
       call t_stopf('UniquePoints')
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
      end if ! fv_nphys > 0

    else ! par%dynproc

       ps_tmp(:,:)        = 0._r8
       T_tmp(:,:,:)       = 0._r8
       uv_tmp(:,:,:,:)    = 0._r8
       omega_tmp(:,:,:)   = 0._r8
       phis_tmp(:,:)      = 0._r8
       q_tmp(:,:,:,:)     = 0._r8
       if (use_gw_front .or. use_gw_front_igw) then
          frontgf(:,:,:) = 0._r8
          frontga(:,:,:) = 0._r8
       end if

    end if ! par%dynproc

    if (use_gw_front .or. use_gw_front_igw) then
       call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
       call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
    end if

    !$omp parallel do num_threads(max_num_threads) private (col_ind, lchnk, icol, ie, blk_ind, ilyr, m)
    do col_ind = 1, columns_on_task
       call get_dyn_col_p(col_ind, ie, blk_ind)
       call get_chunk_info_p(col_ind, lchnk, icol)
       phys_state(lchnk)%ps(icol)=ps_tmp(blk_ind(1),ie)
       phys_state(lchnk)%phis(icol)=phis_tmp(blk_ind(1),ie)
       do ilyr=1,pver
          phys_state(lchnk)%t(icol,ilyr)=T_tmp(blk_ind(1),ilyr,ie)
          phys_state(lchnk)%u(icol,ilyr)=uv_tmp(blk_ind(1),1,ilyr,ie)
          phys_state(lchnk)%v(icol,ilyr)=uv_tmp(blk_ind(1),2,ilyr,ie)
          phys_state(lchnk)%omega(icol,ilyr)=omega_tmp(blk_ind(1),ilyr,ie)

          if (use_gw_front .or. use_gw_front_igw) then
             pbuf_frontgf(icol,ilyr) = frontgf(blk_ind(1),ilyr,ie)
             pbuf_frontga(icol,ilyr) = frontga(blk_ind(1),ilyr,ie)
          endif
       end do ! ilyr

       do m=1,pcnst
          do ilyr=1,pver
             phys_state(lchnk)%q(icol,ilyr,m)=Q_tmp(blk_ind(1),ilyr,m,ie)
          end do ! ilyr
       end do ! m
    end do ! icol

    call t_startf('derived_phys')
    call derived_phys(phys_state,phys_tend,pbuf2d)
    call t_stopf('derived_phys')

!for theta there is no need to multiply omega_p by p
#ifndef MODEL_THETA_L
!$omp parallel do private (lchnk, ncols, ilyr, icol)
    do lchnk=begchunk,endchunk
       ncols=get_ncols_p(lchnk)
       do ilyr=1,pver
          do icol=1,ncols
!!$          if (.not.single_column) then
            phys_state(lchnk)%omega(icol,ilyr) = phys_state(lchnk)%omega(icol,ilyr) &
                                                *phys_state(lchnk)%pmid(icol,ilyr)
!!$          end if
        end do ! icol
      end do ! ilyr
    end do ! lchnk
#endif

   if (write_inithist() ) then
      if (fv_nphys > 0) then

        ncol_d = np*np
        do ie = 1,nelemd
          ncols = elem(ie)%idxP%NumUniquePts
          call outfld('PS&IC',elem(ie)%state%ps_v(:,:,tl_f),  ncol_d,ie)
          call outfld('U&IC', elem(ie)%state%V(:,:,1,:,tl_f), ncol_d,ie)
          call outfld('V&IC', elem(ie)%state%V(:,:,2,:,tl_f), ncol_d,ie)
          call get_temperature(elem(ie),temperature,hvcoord,tl_f)
          call outfld('T&IC',temperature,ncol_d,ie)
          do m = 1,pcnst
            call outfld(trim(cnst_name(m))//'&IC',elem(ie)%state%Q(:,:,:,m), ncol_d,ie)
          end do ! m
        end do ! ie

      else

         do lchnk=begchunk,endchunk
            call outfld('T&IC',phys_state(lchnk)%t,pcols,lchnk)
            call outfld('U&IC',phys_state(lchnk)%u,pcols,lchnk)
            call outfld('V&IC',phys_state(lchnk)%v,pcols,lchnk)
            call outfld('PS&IC',phys_state(lchnk)%ps,pcols,lchnk)
            do m=1,pcnst
               call outfld(trim(cnst_name(m))//'&IC',phys_state(lchnk)%q(1,1,m), pcols,lchnk)
            end do ! m
         end do ! lchnk

      end if ! fv_nphys > 0
    end if ! write_inithist

  end subroutine d_p_coupling
  !=================================================================================================
  !=================================================================================================
  subroutine p_d_coupling(phys_state, phys_tend,  dyn_in)
    use shr_vmath_mod, only: shr_vmath_log
    use cam_control_mod, only : adiabatic
    use control_mod,             only: ftype
    use edge_mod,                only: edge_g, edgeVpack_nlyr, edgeVunpack_nlyr
    use gllfvremap_mod,          only: gfr_fv_phys_to_dyn
    use time_manager,            only: get_step_size
    implicit none
    ! INPUT PARAMETERS:
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend
    ! OUTPUT PARAMETERS:
    type(dyn_import_t),  intent(inout)   :: dyn_in

    ! LOCAL VARIABLES
    integer :: ic , ncols                        ! index
    type(element_t), pointer :: elem(:)          ! pointer to dyn_in element array
    integer                  :: ie, iep          ! indices over elements
    integer                  :: lchnk, icol, ilyr! indices over chunks, columns, layers

    real(r8)    :: T_tmp(npsq,pver,nelemd)       ! temporary array to hold T
    real(r8)    :: uv_tmp(npsq,2,pver,nelemd)    ! temporary array to hold uv
    real(r8)    :: q_tmp(npsq,pver,pcnst,nelemd) ! temporary array to hold q
    real(r8)    :: fmtmp(np,np,nlev)
    integer     :: blk_ind(1), m, i, j, k
    integer     :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)
    integer     :: nphys, nphys_sq           ! physics grid parameters
    integer     :: col_ind                   ! index over columns
    integer     :: nlev_tot, kptr            ! level indices

    if (.not. local_dp_map) then
       call endrun('p_d_coupling: Weak scaling does not support load balancing')
    end if

    if (iam .lt. par%nprocs) then
       elem => dyn_in%elem
    else
       nullify(elem)
    end if

    if (fv_nphys > 0) then
      nphys = fv_nphys
    else
      nphys = np
    end if
    nphys_sq = nphys*nphys

    T_tmp=0.0_r8
    uv_tmp=0.0_r8
    q_tmp=0.0_r8

    if(adiabatic) return

    call t_startf('pd_copy')

    !$omp parallel do num_threads(max_num_threads) private (col_ind, lchnk, icol, ie, blk_ind, ilyr, m)
    do col_ind = 1, columns_on_task
       call get_dyn_col_p(col_ind, ie, blk_ind)
       call get_chunk_info_p(col_ind, lchnk, icol)
       do ilyr=1,pver
          T_tmp(blk_ind(1),ilyr,ie)      = phys_tend(lchnk)%dtdt(icol,ilyr)
          uv_tmp(blk_ind(1),1,ilyr,ie)   = phys_tend(lchnk)%dudt(icol,ilyr)
          uv_tmp(blk_ind(1),2,ilyr,ie)   = phys_tend(lchnk)%dvdt(icol,ilyr)

          do m=1,pcnst
             q_tmp(blk_ind(1),ilyr,m,ie) = phys_state(lchnk)%q(icol,ilyr,m)
          end do
       end do
    end do

    call t_stopf('pd_copy')

    if (par%dynproc) then
       if (fv_nphys > 0) then
          call t_startf('fv_phys_to_dyn')
          ! Map FV physics state to dynamics grid
          call gfr_fv_phys_to_dyn(par, dom_mt, TimeLevel%n0, hvcoord, elem, T_tmp, &
               uv_tmp, q_tmp)
          call t_stopf('fv_phys_to_dyn')

       else ! physics is on GLL nodes

          call t_startf('putUniquePoints')
          do ie=1,nelemd
             ncols = elem(ie)%idxP%NumUniquePts
             call putUniquePoints(elem(ie)%idxP,    nlev,       T_tmp(1:ncols,:,ie),   &
                  elem(ie)%derived%fT(:,:,:))
             call putUniquePoints(elem(ie)%idxP, 2, nlev,       uv_tmp(1:ncols,:,:,ie), &
                  elem(ie)%derived%fM(:,:,:,:))
             call putUniquePoints(elem(ie)%idxP,    nlev,pcnst, q_tmp(1:ncols,:,:,ie),   &
                  elem(ie)%derived%fQ(:,:,:,:))
          end do ! ie
          call t_stopf('putUniquePoints')

       end if ! fv_nphys > 0
    end if ! par%dynproc

    ! Boundary exchange for physics forcing terms.
   ! For physics on GLL grid, for points with duplicate degrees of freedom,
   ! putuniquepoints() set one of the element values and set the others to zero,
   ! so do a simple sum (boundary exchange with no weights).
   ! For physics grid, we interpolated into all points, so do weighted average.

   call t_startf('p_d_coupling:bndry_exchange')

   nlev_tot=(3+pcnst)*nlev
   do ie = 1, nelemd
      if (fv_nphys > 0) then
         do k = 1, nlev
            dyn_in%elem(ie)%derived%FM(:,:,1,k) =                          &
                 dyn_in%elem(ie)%derived%FM(:,:,1,k) *                     &
                 dyn_in%elem(ie)%spheremp(:,:)
            dyn_in%elem(ie)%derived%FM(:,:,2,k) =                          &
                 dyn_in%elem(ie)%derived%FM(:,:,2,k) *                     &
                 dyn_in%elem(ie)%spheremp(:,:)
            dyn_in%elem(ie)%derived%FT(:,:,k) =                            &
                 dyn_in%elem(ie)%derived%FT(:,:,k) *                       &
                 dyn_in%elem(ie)%spheremp(:,:)
            do m = 1, pcnst
               dyn_in%elem(ie)%derived%FQ(:,:,k,m) =                       &
                    dyn_in%elem(ie)%derived%FQ(:,:,k,m) *                  &
                    dyn_in%elem(ie)%spheremp(:,:)
            end do
         end do
      end if
      kptr = 0
      ! fmtmp can be removed if theta and preqx model had the same size FM array
      fmtmp=dyn_in%elem(ie)%derived%FM(:,:,1,:)
      call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,fmtmp,nlev,kptr,nlev_tot)
      kptr=kptr+nlev

      fmtmp=dyn_in%elem(ie)%derived%FM(:,:,2,:)
      call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,fmtmp,nlev,kptr,nlev_tot)
      kptr=kptr+nlev


      call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,dyn_in%elem(ie)%derived%FT(:,:,:),nlev,kptr,nlev_tot)
      kptr=kptr+nlev
      call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,dyn_in%elem(ie)%derived%FQ(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
   end do

   if (iam < par%nprocs) then
     call bndry_exchangeV(par, edge_g)
   end if

   do ie = 1, nelemd
      kptr = 0
      ! fmtmp can be removed if theta and preqx model had the same size FM array
      fmtmp=dyn_in%elem(ie)%derived%FM(:,:,1,:)
      call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,fmtmp,nlev,kptr,nlev_tot)
      kptr=kptr+nlev

      fmtmp=dyn_in%elem(ie)%derived%FM(:,:,2,:)
      call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,fmtmp,nlev,kptr,nlev_tot)
      kptr=kptr+nlev


      call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,dyn_in%elem(ie)%derived%FT(:,:,:),nlev,kptr,nlev_tot)
      kptr=kptr+nlev
!jt      call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,dyn_in%elem(ie)%derived%FQ(:,:,:,:),nlev*qsize,kptr,nlev_tot)
      call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,dyn_in%elem(ie)%derived%FQ(:,:,:,:),nlev*pcnst,kptr,nlev_tot)
      if (fv_nphys > 0) then
         do k = 1, nlev
            dyn_in%elem(ie)%derived%FM(:,:,1,k) =                             &
                 dyn_in%elem(ie)%derived%FM(:,:,1,k) *                        &
                 dyn_in%elem(ie)%rspheremp(:,:)
            dyn_in%elem(ie)%derived%FM(:,:,2,k) =                             &
                 dyn_in%elem(ie)%derived%FM(:,:,2,k) *                        &
                 dyn_in%elem(ie)%rspheremp(:,:)
            dyn_in%elem(ie)%derived%FT(:,:,k) =                               &
                 dyn_in%elem(ie)%derived%FT(:,:,k) *                          &
                 dyn_in%elem(ie)%rspheremp(:,:)
!jtn            do m = 1, qsize
            do m = 1, pcnst
               dyn_in%elem(ie)%derived%FQ(:,:,k,m) =                          &
                    dyn_in%elem(ie)%derived%FQ(:,:,k,m) *                     &
                    dyn_in%elem(ie)%rspheremp(:,:)
            end do
         end do
      end if
   end do
   call t_stopf('p_d_coupling:bndry_exchange')


  end subroutine p_d_coupling
  !=================================================================================================
  !=================================================================================================
  subroutine derived_phys(phys_state, phys_tend, pbuf2d)
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk
    use constituents,  only: qmin
    use physconst,     only: cpair, gravit, zvir, cappa
    use spmd_utils,    only: masterproc
    use ppgrid,        only: pver
    use geopotential,  only: geopotential_t
    use physics_types, only: set_state_pdry, set_wet_to_dry
    use check_energy,  only: check_energy_timestep_init
    use hycoef,   only : hyam, hybm, hyai, hybi, ps0
    use shr_vmath_mod, only: shr_vmath_log
    use gmean_mod,     only: gmean


    implicit none
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc),      pointer     :: pbuf2d(:,:)

    integer :: lchnk
    real(r8) :: qbot                 ! bottom level q before change
    real(r8) :: qbotm1               ! bottom-1 level q before change
    real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
    real(r8) :: qmavl                ! available q at level pver-1

    real(r8) :: ke(pcols,begchunk:endchunk)
    real(r8) :: se(pcols,begchunk:endchunk)
    real(r8) :: ke_glob(1),se_glob(1)
    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer

    integer :: m, i, k, ncol

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

!
! Evaluate derived quantities
!
!$omp parallel do private (lchnk, ncol, k, i, zvirv, pbuf_chnk)
    do lchnk = begchunk,endchunk
       ncol = get_ncols_p(lchnk)
       do k=1,nlev
          do i=1,ncol
             phys_state(lchnk)%pint(i,k)=hyai(k)*ps0+hybi(k)*phys_state(lchnk)%ps(i)
             phys_state(lchnk)%pmid(i,k)=hyam(k)*ps0+hybm(k)*phys_state(lchnk)%ps(i)
          end do
        call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k), &
                           phys_state(lchnk)%lnpint(1:ncol,k),ncol)
        call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k), &
                           phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
       end do
       do i=1,ncol
          phys_state(lchnk)%pint(i,pverp)=hyai(pverp)*ps0+hybi(pverp)*phys_state(lchnk)%ps(i)
       end do
      call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp), &
                         phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)

       do k=1,nlev
          do i=1,ncol
          phys_state(lchnk)%pdel (i,k)  = phys_state(lchnk)%pint(i,k+1) &
                                         -phys_state(lchnk)%pint(i,k)
             phys_state(lchnk)%rpdel(i,k) = 1._r8/phys_state(lchnk)%pdel(i,k)
             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                             / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do

!-----------------------------------------------------------------------------------
!  Need to fill zvirv 2D variables to be compatible with geopotential_t interface
!-----------------------------------------------------------------------------------
       zvirv(:,:) = zvir
!
! Compute initial geopotential heights
!

      ! Compute initial geopotential heights
      call geopotential_t(phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  ,&
                          phys_state(lchnk)%pint  , phys_state(lchnk)%pmid    ,&
                          phys_state(lchnk)%pdel  , phys_state(lchnk)%rpdel   ,&
                          phys_state(lchnk)%t     , phys_state(lchnk)%q(:,:,1),&
                          rairv(:,:,lchnk)        , gravit, zvirv             ,&
                          phys_state(lchnk)%zi    , phys_state(lchnk)%zm      ,&
                          ncol)

       ! Compute initial dry static energy s = g*z + c_p*T
       do k = 1, pver
          do i=1,ncol
             phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                                      + gravit*phys_state(lchnk)%zm(i,k)  &
                                      + phys_state(lchnk)%phis(i)
          end do
       end do

! NOTE:  if a tracer is marked "dry", that means physics wants it dry
!        if dycore advects it wet, it should be converted here
!        FV dycore does this, and in physics/cam/tphysac.F90 it will
!        be converted back to wet, BUT ONLY FOR FV dycore
!
!        EUL: advects all tracers (except q1) as dry.  so it never
!        calls this.
!
!        SE:  we should follow FV and advect all tracers wet (especially
!        since we will be switching to conservation form of advection).
!        So this is broken since dry tracers will never get converted
!        back to wet. (in APE, all tracers are wet, so it is ok for now)
!
! Convert dry type constituents from moist to dry mixing ratio
       call set_state_pdry(phys_state(lchnk))	 ! First get dry pressure to use for this timestep
       call set_wet_to_dry(phys_state(lchnk))    ! Dynamics had moist, physics wants dry.
!
! Ensure tracers are all positive
!
!jt       Could cause differences if checking BFB
       call qneg3('D_P_COUPLING',lchnk  ,ncol    ,pcols   ,pver    , &
                  1, pcnst, qmin  ,phys_state(lchnk)%q)

! Compute energy and water integrals of input state
       pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
       call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)


#if 0
       ke(:,lchnk) = 0._r8
       se(:,lchnk) = 0._r8
!       wv = 0._r8
!       wl = 0._r8
!       wi = 0._r8
       do k = 1, pver
          do i = 1, ncol
             ke(i,lchnk) = ke(i,lchnk) + ( 0.5_r8*(phys_state(lchnk)%u(i,k)**2 + &
                  phys_state(lchnk)%v(i,k)**2)*phys_state(lchnk)%pdel(i,k) )/gravit
             se(i,lchnk) = se(i,lchnk) + phys_state(lchnk)%s(i,k         )*phys_state(lchnk)%pdel(i,k)/gravit
!             wv = wv + phys_state(lchnk)%q(i,k,1       )*phys_state(lchnk)%pdel(i,k)
!             wl = wl + phys_state(lchnk)%q(i,k,ixcldliq)*phys_state(lchnk)%pdel(i,k)
!             wi = wi + phys_state(lchnk)%q(i,k,ixcldice)*phys_state(lchnk)%pdel(i,k)
          end do
       end do
#endif
    end do

#if 0
! This wont match SE exactly.  SE computes KE at half levels
! SE includes cp_star (physics SE uses cp )
! CAM stdout of total energy also includes latent energy of Q,Q1,Q2
! making it a little larger
    call gmean(ke,ke_glob,1)
    call gmean(se,se_glob,1)
    if (masterproc) then
       write(iulog,'(a,e20.8,a,e20.8)') 'KE = ',ke_glob(1),' SE = ',se_glob(1)
       write(iulog,'(a,e20.8)') 'TOTE = ',ke_glob(1)+se_glob(1)
    endif
#endif

  end subroutine derived_phys


end module dp_coupling

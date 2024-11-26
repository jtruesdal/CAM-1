module namelist_mod
  !-----------------
  use cam_logfile,    only: iulog
  !-----------------
  use params_mod,     only: recursive, sfcurve
  !-----------------
  use shr_string_mod, only: shr_string_toUpper
  use shr_kind_mod,   only: r8=>shr_kind_r8
  !-----------------
  use control_mod,    only:   &
       partmethod,            & ! Mesh partitioning method (METIS)
       multilevel,            &
       numnodes,              &
       tasknum,               & ! used dg model in AIX machine
       remapfreq,             & ! number of steps per remapping call
       statefreq,             & ! number of steps per printstate call
       runtype,               &
       cubed_sphere_map,      &
       limiter_option,        &
       nu_top,                &
!jt       hypervis_scaling,      & ! use tensor HV instead of scalar coefficient
       hypervis_scaling         ! use tensor HV instead of scalar coefficient
!jt       hypervis_power,        &
!jt       columnpackage

  !-----------------
  use thread_mod, only : omp_get_max_threads, max_num_threads, horz_num_threads, vert_num_threads, tracer_num_threads
  !-----------------
  use dimensions_mod, only : ne, np, npdg, nnodes, nmpi_per_node, npart, qsize, qsize_d, set_mesh_dimensions
  !-----------------
  !-----------------
  use cam_abortutils, only: endrun
!jt  use parallel_mod,   only: parallel_t, partitionfornodes, useframes
  use parallel_mod,   only: parallel_t, partitionfornodes
  !-----------------


  use interpolate_mod, only : set_interp_parameter, get_interp_parameter

!=============================================================================!
  implicit none
  private
!
! This module should contain no global data and should only be 'use'd to
!    call one of the public interfaces below
!
  public :: homme_set_defaults
  public :: homme_postprocess_namelist

 contains

  ! ============================================
  ! homme_set_defaults:
  !
  !  Set default values for namelist variables
  !
  ! ============================================
  subroutine homme_set_defaults()
    npart               = 1
!jt    useframes           = 0
    multilevel          = 1
    numnodes            = -1
    runtype             = 0
    statefreq           = 1
    remapfreq           = 240
    tasknum             =-1
!jt    columnpackage       = "none"
    nu_top              = 0
    ne                  = 0

  end subroutine homme_set_defaults

  subroutine homme_postprocess_namelist(mesh_file, par)
    use mesh_mod,        only: MeshOpen
    use dimensions_mod,  only: ntrac, ne, ne_x, ne_y
    use control_mod,     only: nu, nu_div, nu_p, nu_s, nu_q, rsplit, &
                               vert_remap_q_alg, vert_remap_u_alg
    use physical_constants, only : scale_factor, scale_factor_inv, domain_size, laplacian_rigid_factor

!!$    use control_mod,     only: dcmip16_mu, dcmip16_mu_s, dcmip16_mu_q

    ! Dummy arguments
    character(len=*),  intent(in) :: mesh_file
    type (parallel_t), intent(in) :: par

    ! Local variable
    character(len=*), parameter :: subname = 'HOMME_POSTPROCESS_NAMELIST: '

    ! set defautl for dynamics remap
    if (vert_remap_u_alg == -2) vert_remap_u_alg = vert_remap_q_alg
    
    ! more thread error checks:  
#ifdef HORIZ_OPENMP
    if(par%masterproc) write(iulog,*)'-DHORIZ_OPENMP enabled'
#else
    if(par%masterproc) write(iulog,*)'-DHORIZ_OPENMP disabled'
#endif
#ifdef COLUMN_OPENMP
    if(par%masterproc) write(iulog,*)'-DCOLUMN_OPENMP enabled'
#else
    if(par%masterproc) write(iulog,*)'-DCOLUMN_OPENMP disabled'
#endif

    if (ne /=0 .or. ne_x /=0 .or. ne_y /=0) then
       if (mesh_file /= "none" .and. mesh_file /= "/dev/null") then
          write (*,*) "namelist_mod: mesh_file:",trim(mesh_file), &
               " and ne/ne_x/ne_y:",ne,ne_x,ne_y," are both specified in the input file."
          write (*,*) "Specify one or the other, but not both."
          call abortmp("Do not specify ne (or ne_x, ne_y) if using a mesh file input.")
       end if
    end if
    if (par%masterproc) write (iulog,*) "Mesh File:", trim(mesh_file)
    if (ne.eq.0 .and. ne_x .eq. 0 .and. ne_y .eq. 0) then
#ifndef HOMME_WITHOUT_PIOLIBRARY
       call set_mesh_dimensions()
       if (par%masterproc) write (iulog,*) "Opening Mesh File:", trim(mesh_file)
       call MeshOpen(mesh_file, par)
#else
       call abortmp("Build is without PIO library, mesh runs (ne=0) are not supported.")
#endif
    end if
    ! set map
    if (cubed_sphere_map<0) then
#if ( defined MODEL_THETA_C || defined MODEL_THETA_L ) 
       cubed_sphere_map=2  ! theta model default = element local
#else
       cubed_sphere_map=0  ! default is equi-angle gnomonic
#endif
    endif
    if (ne.eq.0 .and. ne_x .eq. 0 .and. ne_y .eq. 0) cubed_sphere_map=2  ! must use element_local for var-res grids
    if (par%masterproc) write (iulog,*) "Reference element projection: cubed_sphere_map=",cubed_sphere_map

    scale_factor = rearth
    scale_factor_inv = rrearth
    domain_size = 4.0D0*DD_PI
    laplacian_rigid_factor = rrearth

#ifdef _PRIM
    if (limiter_option==8 .or. limiter_option==84 .or. limiter_option == 9) then
       if (hypervis_subcycle_q/=1 .and. transport_alg == 0) then
          call abortmp('limiter 8,84,9 require hypervis_subcycle_q=1')
       endif
    endif
    if (transport_alg == 0 .and. dt_remap_factor > 0 .and. dt_remap_factor < dt_tracer_factor) then
       call abortmp('Only SL transport supports vertical remap time step < tracer time step.')
    end if
#endif

    if((prescribed_wind/=0).and.(prescribed_wind/=1))then
          call abortmp('prescribed_wind should be either 0 or 1')
    endif

    nmpi_per_node=1

    ! some default diffusion coefficiets
    if(nu_s<0)    nu_s  = nu
    if(nu_q<0)    nu_q  = nu
    if(nu_div<0)  nu_div= nu
    if(nu_p<0) then                                                                           
       if (rsplit==0) then                                                                    
          nu_p=0  ! eulerian code traditionally run with nu_p=0                               
       else                                                                                   
          nu_p=nu                                                                             
       endif                                                                                  
    endif 
!!$    if(dcmip16_mu_s<0)    dcmip16_mu_s  = dcmip16_mu
!!$    if(dcmip16_mu_q<0)    dcmip16_mu_q  = dcmip16_mu_s

    nnodes = npart/nmpi_per_node
    if(numnodes > 0 ) then
        nnodes = numnodes
        nmpi_per_node = npart/nnodes
    endif


  end subroutine homme_postprocess_namelist
end module namelist_mod

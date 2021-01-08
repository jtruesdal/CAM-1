!-----------------------------------------------------------------------
!  module memstats
!
!> \brief  Provides memory usage statistics, written to the CAM log file
!> \author Michael Duda
!> \date   3 December 2019
!> \details
!>  This module provides routines for measuring heap memory usage statistics
!>  and writing those statistics to CAM's log file.
!>
!>  Before any other routine in this module may be used, the memstats_init
!>  routine must first be called to initialize the module. It is also
!>  recommended to call the memstats_finalize routine after all calls to
!>  routines in the module have been made.
!
!-----------------------------------------------------------------------
module memstats

   use iso_c_binding, only : c_long

   public :: memstats_init, &
             memstats_finalize, &
             memstats_log

   private

   integer :: mpi_comm
   integer :: logunit
   logical :: initialized = .false.
   integer :: comm_size
   integer :: comm_rank


contains


    !-----------------------------------------------------------------------
    !  routine memstats_init
    !
    !> \brief  Pre-computes thread loop bounds for cell, edge, and vertex elements
    !> \author Michael Duda
    !> \date   3 December 2019
    !> \details
    !>  This routine is responsible for computing thread loop bounds for cell,
    !>  edge, and vertex elements in each block of the input blocklist argument.
    !>  Starting and ending loop bounds are computed for these three element
    !>  types for all elements (e.g., nCells) as well as owned elements (e.g.,
    !>  nCellsSolve).
    !
    !-----------------------------------------------------------------------
    subroutine memstats_init(mpi_comm_, logunit_)

        use mpi, only : MPI_Comm_size, MPI_Comm_rank, MPI_SUCCESS

        implicit none

        integer, intent(in) :: mpi_comm_
        integer, intent(in) :: logunit_

        integer :: mpi_ierr

        mpi_comm = mpi_comm_
        logunit = logunit_

        call MPI_Comm_size(mpi_comm, comm_size, mpi_ierr)
        if (mpi_ierr /= MPI_SUCCESS) then
            write(0,'(a)') 'MEM memstats_init: MPI_Comm_size failed'
            return
        end if

        call MPI_Comm_rank(mpi_comm, comm_rank, mpi_ierr)
        if (mpi_ierr /= MPI_SUCCESS) then
            write(0,'(a)') 'MEM memstats_init: MPI_Comm_rank failed'
            return
        end if

        initialized = .true.

        if (comm_rank == 0 .and. mpi_ierr == MPI_SUCCESS) then
            write(logunit,'(a)') 'MEM initializing memory statistics'
        end if

    end subroutine memstats_init


    !-----------------------------------------------------------------------
    !  routine memstats_finalize
    !
    !> \brief  Pre-computes thread loop bounds for cell, edge, and vertex elements
    !> \author Michael Duda
    !> \date   3 December 2019
    !> \details
    !>  This routine is responsible for computing thread loop bounds for cell,
    !>  edge, and vertex elements in each block of the input blocklist argument.
    !>  Starting and ending loop bounds are computed for these three element
    !>  types for all elements (e.g., nCells) as well as owned elements (e.g.,
    !>  nCellsSolve).
    !
    !-----------------------------------------------------------------------
    subroutine memstats_finalize()

        implicit none

        !
        ! At present, this routine is a no-op
        !

    end subroutine memstats_finalize


    !-----------------------------------------------------------------------
    !  routine memstats_log
    !
    !> \brief  Pre-computes thread loop bounds for cell, edge, and vertex elements
    !> \author Michael Duda
    !> \date   3 December 2019
    !> \details
    !>  This routine is responsible for computing thread loop bounds for cell,
    !>  edge, and vertex elements in each block of the input blocklist argument.
    !>  Starting and ending loop bounds are computed for these three element
    !>  types for all elements (e.g., nCells) as well as owned elements (e.g.,
    !>  nCellsSolve).
    !
    !-----------------------------------------------------------------------
    subroutine memstats_log(mesg)

        use mpi, only : MPI_Gather, MPI_INTEGER, MPI_SUCCESS

        implicit none

        character(len=*), intent(in) :: mesg

        character(len=128) :: mesg_log
        integer(c_long) :: maxrss
        integer :: maxrss_i4
        integer, dimension(:), allocatable :: maxrss_array
        integer :: mpi_ierr

        integer :: nsp
        character(len=64), parameter :: spaces = '                                                                '

        interface
            subroutine get_maxrss(maxrss) bind(c)
                use iso_c_binding, only : c_long
                integer(c_long) :: maxrss
            end subroutine get_maxrss
        end interface

        if (initialized) then
            maxrss = mem_usage()
            maxrss_i4 = maxrss
            allocate(maxrss_array(comm_size))
            call MPI_Gather(maxrss_i4, 1, MPI_INTEGER, maxrss_array, 1, MPI_INTEGER, 0, mpi_comm, mpi_ierr)
            if (comm_rank == 0 .and. mpi_ierr == MPI_SUCCESS) then
                nsp = max(0,64-len_trim(mesg))
                write(mesg_log,'(a,a,i9,a1,i9,a1,i9,a)') 'MEM [min/max/sum] '//trim(mesg)//' -> ', spaces(1:nsp), &
                                                      minval(maxrss_array), '/', maxval(maxrss_array), '/', sum(maxrss_array), ' kB'
                write(logunit,'(a)') trim(mesg_log)
            else if (comm_rank == 0) then
                write(logunit,'(a)') 'MEM memstats_log: MPI_Gather failed'
            end if
            deallocate(maxrss_array)
        else
            write(0,'(a)') 'MEM memstats_log: memstats module has not been initialized'
        end if

    end subroutine memstats_log


    !-----------------------------------------------------------------------
    !  routine mem_usage
    !
    !> \brief  Returns the maximum resident set size, in kB
    !> \author Michael Duda
    !> \date   3 December 2019
    !> \details
    !>  This routine returns the maximum resident set size, in kB, by calling
    !>  a C routine.
    !
    !-----------------------------------------------------------------------
    function mem_usage() result(maxrss)

        implicit none

        integer(c_long) :: maxrss

        interface
            subroutine get_maxrss(maxrss) bind(c)
                use iso_c_binding, only : c_long
                integer(c_long) :: maxrss
            end subroutine get_maxrss
        end interface

        call get_maxrss(maxrss)

    end function mem_usage

end module memstats

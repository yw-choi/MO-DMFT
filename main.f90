program MO_DMFT_ED

    use precision
    use parallel_params
    use fdf

    ! DMFT-ED realted modules
    use ed_config
    ! use ed_hamiltonian
    ! use ed_solver
    ! use ed_green_ftn
    ! use ed_minimize
    ! use ed_converged

    implicit none

    ! local variables
    integer :: iloop
    logical :: converged

    ! initializations that is also done in SIESTA part.
    call MPI_Init(ierr)
    call MPI_Comm_size(comm,Nodes,ierr)
    call MPI_Comm_rank(comm,Node,ierr)
    call fdf_init('input.fdf', 'fdf.out')

    ! /////////////////////////////////////////////////////////////////////// !
    ! DMFT PART START
    ! /////////////////////////////////////////////////////////////////////// !

    ! read input parameters
    call ed_config_init

    ! Set band structure 
    ! call ed_set_band_structure

    ! DMFT Loop
    iloop = 0
    converged = .false.
    ! do while(.not.converged.and.iloop<nloop)
        ! call ed_solve

    !     call ed_green_ftn

    !     call ed_delta_new

    !     call ed_minimize 

    !     call ed_converged(converged)

    !     iloop = iloop + 1
    ! enddo
    
    ! @TODO post-processing

    ! /////////////////////////////////////////////////////////////////////// !
    ! DMFT PART END
    ! /////////////////////////////////////////////////////////////////////// !

    call MPI_Finalize(ierr)

contains

end program MO_DMFT_ED


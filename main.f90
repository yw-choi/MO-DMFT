program MO_DMFT_ED

    use precision
    use parallel_params
    use fdf
    use ionew

    ! DMFT-ED realted modules
    use ed_config
    use ed_hamiltonian
    use ed_solver
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
    call io_setup

    if (node.eq.0) then
        write(6,*) "--------------------------------------------------"
        write(6,*) "Number of processors = ", nodes
        write(6,*) "--------------------------------------------------"
    endif

    call timer('DMFT',0)

    if (node.eq.0) then
        write(6,*) "--------------------------------------------------"
        write(6,*) "DMFT PART"
        write(6,*) "--------------------------------------------------"
    endif
    call timestamp2("DMFT PART START")

    ! read input parameters
    call ed_config_init

    if (node.eq.0) then
        write(6,*) "--------------------------------------------------"
        write(6,*) "Seting up band structure"
        write(6,*) "--------------------------------------------------"
    endif
    ! Set band structure 
    call ed_set_band_structure

    if (node.eq.0) then
        write(6,*) "--------------------------------------------------"
        write(6,*) "DMFT loop start"
        write(6,*) "--------------------------------------------------"
    endif
    iloop = 0
    converged = .false.
    do while(.not.converged.and.iloop<nloop)
        call ed_solve

        ! call ed_green_ftn

        ! call ed_delta_new

        ! call ed_minimize 

        ! call ed_converged(converged)

        iloop = iloop + 1
    enddo
    
    ! @TODO post-processing

    call timestamp2("DMFT PART END")
    call timer('DMFT',2)
    call timer('DMFT',3)

    call MPI_Finalize(ierr)

end program MO_DMFT_ED


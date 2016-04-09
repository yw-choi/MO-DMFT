program MO_DMFT_ED

    use precision
    use parallel_params
    use fdf
    use ionew

    ! DMFT-ED realted modules
    use ed_config
    use ed_basis
    use ed_hamiltonian
    use ed_solver
    use ed_green
    ! use ed_minimize
    ! use ed_converged

    implicit none

    ! local variables
    integer :: iloop
    logical :: converged

    integer :: nev_calc

    ! initializations that is also done in SIESTA part.
    call MPI_Init(ierr)
    call MPI_Comm_size(comm,Nodes,ierr)
    call MPI_Comm_rank(comm,Node,ierr)
    call fdf_init('input.fdf', 'fdf.out')
    call io_setup
    call alloc_report( level=4, file='DMFT.alloc', threshold=0.0_dp, printNow=.false. )

    if (node.eq.0) then
        write(6,*) "Number of processors = ", nodes
    endif

    call timer('DMFT',0)
    call timer('DMFT',1)

    call timestamp2("DMFT INITIALIZATION START")

    call ed_read_options
    call ed_hamiltonian_init
    call ed_basis_init
    call ed_set_band_structure

    call timestamp2("DMFT LOOP START")

    iloop = 0
    converged = .false.
    do while(.not.converged.and.iloop<nloop)
        call ed_solve(iloop,nev_calc)

        call ed_calc_green_ftn(nev_calc)

        ! call ed_delta_new

        ! call ed_minimize 

        ! call ed_converged(converged)

        iloop = iloop + 1
    enddo

    call timestamp2("DMFT LOOP END")
    
    ! @TODO post-processing

    call alloc_report( printNow=.true. )
    call timestamp2("DMFT PART END")
    call timer('DMFT',2)
    call timer('all',3)

    call MPI_Finalize(ierr)

end program MO_DMFT_ED


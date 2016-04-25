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
    integer :: iloop, i
    logical :: converged

    integer :: nev_calc

    ! @TODO the follownig variables should be refactored to a separate module
    integer :: nxsize,korb
    real(dp) :: tol, xmin
    real(dp), allocatable :: x(:)

    ! initializations that is also done in SIESTA part.
    call MPI_Init(ierr)
    call MPI_Comm_size(comm,Nodes,ierr)
    call MPI_Comm_rank(comm,Node,ierr)
    call fdf_init('input.fdf', 'fdf.out')
    call io_setup

    if (node.eq.0) then
        write(6,"(a)") repeat("=",80)
        write(6,*) "Multi-orbital DMFT calculation with ED solver"
        write(6,*) "Number of processors = ", nodes
        write(6,*) 
    endif

    call timer('DMFT',0)
    call timer('DMFT',1)

    call timestamp2("DMFT INITIALIZATION START")

    call ed_read_options
    call ed_hamiltonian_init
    call ed_set_band_structure
    call ed_green_init

    nxsize = nbath + norb + nbath*norb
    allocate(x(nxsize))

    call timestamp2("DMFT LOOP START")

    iloop = 0
    converged = .false.
    dmft_loop: do while(.not.converged.and.iloop<nloop)
        call ed_solve(iloop,nev_calc)

        call ed_calc_green_ftn(nev_calc)

        ! open(unit=137,file="green.dump",status="replace",form="formatted")
        ! do i=1,nwloc
        !     write(137,"(4F20.16)") omega(i), real(Gr(1,i)), aimag(Gr(1,i))
        ! enddo
        ! close(137)
        ! stop
        call ed_delta_new

        call n_from_gksum

        ! *********************************************************************
        ! Update ek, vk
        ! *********************************************************************
        call ev_to_x(ek,vk,ef,Nsite,Nbath,nxsize,x)

        call minimization(x,nwloc,Nsite,nxsize,omega,D_ev,Nbath,Norb,comm,xmin)
        if(node.eq.0) then 
            call x_to_ev(x,nxsize,Nsite,Nbath,Norb,ef,ek,vk)
            tol = 0.D0
            do korb = 1, Norb
                tol = tol + &
                    sum(abs(Gr_prev(korb,:)-Gr(korb,:)))/float(Nw)
            enddo
            tol = tol/float(Norb)
            if(node.eq.0) &
                write(6,'(i5,3x,a,x,e)') iloop, "d_max=", tol
            Gr_prev = Gr
        endif
        call mpi_barrier(comm,ierr)
        call mpi_bcast(ef(1),norb,mpi_double_precision,0,comm,ierr)
        call mpi_bcast(ek(1),Nsite,mpi_double_precision,0,comm,ierr)
        call mpi_bcast(vk(1,1),Nbath*Norb,mpi_double_precision,0,comm,ierr)
        call mpi_bcast(tol,1,mpi_double_precision,0,comm,ierr)

        converged = tol.lt.scf_tol

        call calc_dev
        iloop = iloop + 1
    enddo dmft_loop

    call timestamp2("DMFT LOOP END")
    
    ! @TODO post-processing

    call timestamp2("DMFT PART END")
    call timer('DMFT',2)
    call timer('all',3)

    call MPI_Finalize(ierr)

end program MO_DMFT_ED


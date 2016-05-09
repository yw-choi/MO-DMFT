program MO_DMFT_ED

    
    use parallel_params
    use fdf
    use ionew

    ! DMFT-ED realted modules
    use ed_config
    use ed_basis
    use ed_hamiltonian
    use ed_solver
    use ed_green

    implicit none

    ! local variables
    integer :: iloop
    logical :: converged

    ! @TODO the follownig variables should be refactored to a separate module
    integer :: nxsize,korb, i, k
    double precision :: tol, xmin, dw
    double precision, allocatable :: x(:), Zq(:), Aff(:,:)
    double complex, allocatable :: Sigma(:,:)
    character(len=200) :: fn

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

    nxsize = nbath + nbath*norb + norb 
    allocate(x(nxsize))

    call timestamp2("DMFT LOOP START")

    iloop = 0
    converged = .false.
    dmft_loop: do while(.not.converged.and.iloop<nloop)
        if (node.eq.0) then
            write(6,*)
            write(6,'(a)') repeat("=",80)
            write(6,*) "DMFT iteration ", iloop
            write(6,'(a)') repeat("=",80)
        endif

        call timer("solver",1)
        call ed_solve
        call timer("solver",2)

        call timer("greenftn",1)
        call ed_calc_green_ftn(nev_calc)
        call timer("greenftn",2)

        call timer("delta_new",1)
        call ed_delta_new
        call timer("delta_new",2)

        call n_from_gksum

        ! *********************************************************************
        ! Update ek, vk
        ! *********************************************************************
        call timer("minimize",1)
        if (node.eq.0) then
            write(6,*) "minimization: Finding new ek, vk..."
        endif
        call ev_to_x(ek,vk,ef,Nsite,Nbath,nxsize,x)

        call mpi_barrier(comm,ierr)
        call minimization(x,nwloc,Nsite,nxsize,omega,D_ev,Nbath,Norb,comm,xmin)
        if(node.eq.0) then 
            call x_to_ev(x,nxsize,Nsite,Nbath,Norb,ef,ek,vk)
            tol = 0.D0
            do korb = 1, Norb
                tol = tol + &
                    sum(abs(Gr_prev(korb,:)-Gr(korb,:)))/float(Nw)
            enddo
            tol = tol/float(Norb)
            Gr_prev = Gr
        endif
        call mpi_barrier(comm,ierr)
        call mpi_bcast(ef(1),norb,mpi_double_precision,0,comm,ierr)
        call mpi_bcast(ek(1),Nsite,mpi_double_precision,0,comm,ierr)
        call mpi_bcast(vk(1,norb+1),Nbath*Norb,mpi_double_precision,0,comm,ierr)
        call mpi_bcast(tol,1,mpi_double_precision,0,comm,ierr)

        if (node.eq.0) then
            write(6,*) "minimization: Updated ek and vk."
        endif
        call timer("minimize",2)

        call timer("calc_dev",1)
        if (node.eq.0) then
            write(6,*) "calc_dev: Calculating new D_ev"
        endif
        call calc_dev
        call timer("calc_dev",2)

        converged = tol.lt.scf_tol
        if (node.eq.0) then
            write(6,*)
            write(6,"(x,a,ES12.5)") "scf_diff = ", tol
            write(6,*)
        endif
        iloop = iloop + 1
    enddo dmft_loop

    if (node.eq.0) then
        write(6,'(a)') repeat("=",80)
        write(6,*)
        if (converged) then
            write(6,*) "DMFT SCF Convergence has reached at loop ", iloop
            write(6,*)
        else
            write(6,"(a,i,a)") "[FAILURE] DMFT SCF loop did not converge within ",&
                               nloop," iterations."
            call die
        endif
    endif

    call timestamp2("DMFT LOOP END")
    
    if (node.eq.0) then
        write(6,*)
        write(6,*) "Calculating the self-energy..."
    endif

    ! Self-energy
    allocate(Sigma(norb,nwloc))
    do i = 1, nwloc
        do k = 1, Norb
            Sigma(k,i) = cmplx(0.0D0,omega(i)) - D_ev(k,i) - 1.0D0/Gr(k,i)
        enddo
    enddo

    if(node.eq.0) then
        open(unit=100,file="Sigma_dia.dat",form="formatted")
        do i = 1, nwloc
            write(100,'(100f10.5)') omega(i), (Sigma(k,i),k=1,norb)  ! cautious
        enddo
        close(100)

        do k=1,norb
            write(fn,"(A5,I2.2,A4)") "green",k,".dat"
            open(unit=177,file=fn,status="replace",form="formatted")
            do i=1,nwloc
                write(177,*) omega(i), real(Gr(k,i)), aimag(Gr(k,i))
            enddo
            close(177)
        enddo
    endif

    allocate(Zq(norb))
    if(node.eq.0) then
        do i = 1, norb
            Zq(i) = (aimag(Sigma(i,1)))/(omega(1))
            Zq(i) = 1.D0/(1.D0-Zq(i))
        enddo
        write(6,'(a,2x,2e)') "Quasiparticle weight = ",(Zq(i),i=1,norb)
    endif
    deallocate(Sigma,Zq)
    deallocate(omega)

    allocate(omega(Nw))
    allocate(Aff(norb,nw))
    dw = 0.01D0 
    if(node.eq.0) then
        write(6,*)
        write(6,*) "Calculating the spectral function..."
        do i = 1, Nw
            omega(i) = (-Nw/2+i)*dw
        enddo
        call dos(omega,Aff)
        open(unit=101,file="SpectralFtn.dat",form="formatted",status="unknown")
        do i = -Nw/2, Nw/2-1
            write(101,'(6e)') i*dw,(Aff(k,i+Nw/2+1),k=1,norb)
        enddo
        close(101)
        write(6,*) ".... DOS CALCULATION COMPLETED"
    endif
    call mpi_barrier(comm,ierr)

    if(node.eq.0) then
        write(6,'(a,1x,i4,a)') &
            "By ",iloop," iterations, convergence archived" 
        write(6,*) "Converged ek(i), vk(i)"
        do korb = 1, Norb
        write(6,'(a,2x,i2)') "Orbital", korb
        do i = Norb+1, Nsite
        write(6,'(f10.5,4x,f10.5)') ek(i),vk(korb,i)
        enddo
        write(6,*)
        enddo
    endif
    call mpi_barrier(comm,ierr)

    call timestamp2("DMFT PART END")
    call timer('DMFT',2)
    call timer('all',3)

    call MPI_Finalize(ierr)
end program MO_DMFT_ED


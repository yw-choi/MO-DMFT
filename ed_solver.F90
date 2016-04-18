module ed_solver

    use precision
    use ed_config
    use ed_hamiltonian
    use ed_basis
    use parallel_params
    use sys
    use ed_utils
    use ed_io
    implicit none

    public
contains

    subroutine ed_solve(iloop,nev_calc)
        integer, intent(in) :: iloop
        integer, intent(out) :: nev_calc

        integer :: isector, i, ne_up, ne_down
        real(dp), allocatable :: eigvec(:,:)

        real(dp) :: eigval(nev), eigval_all(nev*nsector)
        real(dp) :: pev(nev*nsector)
        integer :: ind(nsector*nev)
        type(basis_t) :: basis


        call timer('ed_solve',1)
        do isector=1,nsector
            ne_up = nup(isector)
            ne_down = nelec(isector) - nup(isector)

            if (node.eq.0) then
                write(*,*) "ed_solver: iloop = ",iloop," , isector = ",isector
                write(*,*), "ed_solver: ne_up=", ne_up, " ne_down=",ne_down
            endif

            basis = generate_basis( ne_up, ne_down )

            allocate(eigvec(basis%nloc,nev))

            call diag(basis,eigval,eigvec)

            do i=1,nev
                eigval_all((isector-1)*nev+i) = eigval(i)
                call export_eigvec(isector,i,node,basis%nloc,nev,eigvec(:,i))
            enddo

            deallocate(eigvec)
        enddo

        call sort(nev*nsector,eigval_all,ind)
        call boltzmann_factor(eigval_all,nev*nsector,beta,pev)

        nev_calc = 0

        do i = 1, nev
            if(pev(ind(i)).gt.0.001D0) nev_calc = nev_calc + 1
        enddo
     
        if(Node.eq.0) then
            write(6,'(i3,2x,a)') nev_calc, &
                "levels will be considered( p > 1/1000 )"
            write(6,*) "       Energy               probability" 
            do i = 1, nev_calc
                write(6,*) eigval_all(ind(i)), pev(ind(i))
            enddo
        endif

        if(node.eq.0)  then
            call export_eigval(nev_calc,eigval_all,pev,ind)
        endif

        call mpi_barrier(comm,ierr)
        call timer('ed_solve',2)
        return
    end subroutine ed_solve

    subroutine diag(basis, eigval, eigvec)
#ifdef DEBUG
        include 'debug.h'
#endif
        include 'stat.h'
        type(basis_t), intent(in) :: basis
        real(dp), intent(out) :: eigvec(basis%nloc,nev)
        real(dp), intent(out) :: eigval(nev)

        integer  maxnloc,maxnev,maxncv,ldv
        parameter  (maxnloc=3500000,maxnev=40,maxncv=60,ldv=maxnloc)
        real(dp) :: v(ldv,maxncv), workl(maxncv*(maxncv+8)),        &
                    workd(3*maxnloc), d(maxncv,2), resid(maxnloc),  &
                    ax(maxnloc)
        logical :: select(maxncv)
        character, parameter :: bmat = 'I'
        character(len=2), parameter :: which = 'SA'
        integer :: iparam(11), ipntr(11), lworkl, info, ido, nconv, i, j, ncv
        real(dp) :: sigma, tol
        real(dp) :: pdnorm2
        external :: pdnorm2, daxpy

#ifdef DEBUG
        ndigit = -3
        logfil = 6
        msaupd = 1
#endif 
        ncv = nev*2

        lworkl = ncv*(ncv+8)
        tol = 0.0

        iparam(1) = 1   ! ishfts
        iparam(3) = 500 ! maxitr
        iparam(7) = 1   ! mode

        info = 0
        ido = 0

        do
            call pdsaupd( comm, ido, bmat, basis%nloc, which, nev, tol, resid, &
                ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
            if (ido .eq. -1 .or. ido .eq. 1) then
                call multiply_H(basis,workd(ipntr(1)),workd(ipntr(2)))
            else
                exit
            endif
        enddo
        if ( info .lt. 0 ) then
            if ( node.eq. 0 ) then
                print *, ' Error with pdsaupd, info = ', info
                print *, iparam(5)
                call die
            endif
        else
            call pdseupd(comm,.true.,'All', select,d,v,ldv,sigma, &
                bmat, basis%nloc, which, nev, tol, resid, ncv, v, ldv, &
                iparam, ipntr, workd, workl, lworkl, ierr )
            if ( ierr .ne. 0) then
                if ( node .eq. 0 ) then
                    print *, ' Error with pdseupd, info = ', ierr
                    call die
                endif
            else
                nconv =  iparam(5)
                do j=1, nconv
                    call multiply_H(basis,v(1:ldv,j),ax)
                    call daxpy(basis%nloc, -d(j,1), v(1,j), 1, ax, 1)
                    d(j,2) = pdnorm2( comm, basis%nloc, ax, 1 )
                enddo
                call pdmout(comm, 6, nconv, 2, d, ncv, -6, &
                    'Ritz values and direct residuals')
            end if
        endif

        do i=1,nev
            eigvec(1:basis%nloc,i) = v(:,i)
            eigval(i) = d(i,1)
        enddo
    end subroutine diag

end module ed_solver

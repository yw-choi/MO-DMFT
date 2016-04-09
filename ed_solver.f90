module ed_solver

    use precision
    use ed_config
    use ed_hamiltonian
    use ed_basis
    use parallel_params
    use sys
    use ed_utils
    use alloc

    implicit none

    real(dp), pointer :: eigval(:), eigval_all(:)
    real(dp), pointer :: eigvec(:,:)

    public
contains

    subroutine ed_solver_init

        call re_alloc(eigval,1,nev,name="eigval",routine="ed_solve")
        call re_alloc(eigval_all,1,nev*nsector,name="eigval_all",routine="ed_solve")

    end subroutine ed_solver_init

    subroutine ed_solve(iloop)
        integer, intent(in) :: iloop
        integer :: isector,i
        call timer('ed_solve',1)

        do isector=1,nsector
            call prepare_basis_for_sector(isector)
            call re_alloc(eigvec,1,nbasis_loc,1,nev,"eigvec","ed_solve", &
                         copy=.false.,shrink=.true.)

            call diag(isector)

            do i=1,nev
                eigval_all((isector-1)*nev+i) = eigval(i)
            enddo
        enddo

        call timer('ed_solve',2)
    end subroutine ed_solve

    subroutine diag(isector)
#ifdef DEBUG
        include 'debug.h'
#endif
        include 'stat.h'
        integer, intent(in) :: isector

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
#endif DEBUG
        ncv = nev*2

        lworkl = ncv*(ncv+8)
        tol = 0.0

        iparam(1) = 1   ! ishfts
        iparam(3) = 500 ! maxitr
        iparam(7) = 1   ! mode

        info = 0
        ido = 0

        do
            call pdsaupd( comm, ido, bmat, nbasis_loc, which, nev, tol, resid, &
                ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
            if (ido .eq. -1 .or. ido .eq. 1) then
                call multiply_H(nbasis(isector),nbasis_loc,workd(ipntr(1)),workd(ipntr(2)))
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
                bmat, nbasis_loc, which, nev, tol, resid, ncv, v, ldv, &
                iparam, ipntr, workd, workl, lworkl, ierr )
            if ( ierr .ne. 0) then
                if ( node .eq. 0 ) then
                    print *, ' Error with pdseupd, info = ', ierr
                    call die
                endif
            else
                nconv =  iparam(5)
                do j=1, nconv
                    call multiply_H(nbasis(isector),nbasis_loc,v(1:ldv,j),ax)
                    call daxpy(nbasis_loc, -d(j,1), v(1,j), 1, ax, 1)
                    d(j,2) = pdnorm2( comm, nbasis_loc, ax, 1 )
                enddo
                call pdmout(comm, 6, nconv, 2, d, ncv, -6, &
                    'Ritz values and direct residuals')
            end if
        endif

        do i=1,nev
            eigvec(:,i) = v(:,i)
            eigval(i) = d(i,1)
        enddo
    end subroutine diag

end module ed_solver

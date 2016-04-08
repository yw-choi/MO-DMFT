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
        integer :: isector,i,iunit
        character(len=100) :: fname

        call timer('ed_solve',1)

        do isector=1,nsector
            call prepare_basis_for_sector(isector)
            call re_alloc(eigvec,1,nbasis_loc,1,nev,"eigvec","ed_solve", &
                         copy=.false.,shrink=.true.)
#ifdef DEBUG
            iunit = 11+node
            write(fname,"(A,I1)") "basis.node-",node
            open(unit=iunit,file=fname,status="replace")
            do i=1,nbasis_loc
                write(iunit,"(I5,B)") i,ed_basis_get(offsets(node)+i)
            enddo
            close(iunit)
#endif
        enddo

        call timer('ed_solve',2)
    end subroutine ed_solve

end module ed_solver

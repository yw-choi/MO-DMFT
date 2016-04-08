module ed_solver

    use precision
    use ed_config
    use ed_hamiltonian
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
        integer :: isector
        call timer('ed_solve',1)

        do isector=1,nsector
            call prepare_hamiltonian_sector(isector)
            call re_alloc(eigvec,1,nbasis_loc,1,nev,"eigvec","ed_solve", &
                         copy=.false.,shrink=.true.)

            call end_hamiltonian_sector
        enddo

        call timer('ed_solve',2)
    end subroutine ed_solve

end module ed_solver

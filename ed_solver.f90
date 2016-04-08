module ed_solver

    use precision
    use ed_config
    use ed_hamiltonian
    use parallel_params
    use sys
    use ed_utils

    implicit none

    real(dp), allocatable :: eigval(:), eigval_all(:)
    real(dp), allocatable :: eigvec(:,:)

contains

    subroutine ed_solve
        integer :: isector

        if (.not.allocated(eigval)) allocate(eigval(nev))
        if (.not.allocated(eigval_all)) allocate(eigval_all(nev*nsector))

        do isector=1,nsector
            call prepare_hamiltonian_sector
            allocate(eigvec(nbasis_loc,nev))

            

            call end_hamiltonian_sector
            deallocate(eigvec)
        enddo

    end subroutine ed_solve

end module ed_solver

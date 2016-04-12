module ed_lanczos
    use precision
    use ed_config
    use ed_basis
    use ed_utils
    use ed_hamiltonian
    use parallel_params
    
    implicit none

contains

    subroutine lanczos_iteration(basis,vec,nstep_calc,aout,bout)

        type(basis_t), intent(in) :: basis
        real(dp), intent(in) :: vec(basis%nloc) ! starting vector for lanczos iteration

        integer, intent(out) :: nstep_calc
        real(dp), allocatable, intent(out) :: aout(:), bout(:)

        ! lanczos vectors
        real(dp), allocatable :: a(:), b(:)
        real(dp), allocatable :: v(:,:), w(:)
        real(dp) :: norm_v
        
        integer :: j
        
        allocate(v(basis%nloc,2),w(basis%nloc))
        allocate(a(nstep),b(nstep))

        v(1:basis%nloc,2) = 0.0_dp
        w(1:basis%nloc) = 0.0_dp

        a(1:nstep) = 0.0_dp
        b(1:nstep) = 0.0_dp

        v(1:basis%nloc,2) = vec

        ! Lanczos steps
        ! ref: https://en.wikipedia.org/wiki/Lanczos_algorithm#Iteration 

        ! normalize the initial vector
        norm_v = mpi_norm( v(1:basis%nloc,2), basis%nloc )
        v(1:basis%nloc,2) = v(1:basis%nloc,2)/norm_v
        
        ! v(:,1) = v_(j-1)
        ! v(:,2) = v_j
        ! w(:)   = w_j
        lanczos_loop: do j=1,nstep-1
            nstep_calc = j

            ! w_j = H*v_j
            call multiply_H( basis, v(:,2), w(:) )

            ! a_j = dot(w_j,v_j)
            a(j) = mpi_dot_product(w(:), v(:,2), basis%nloc)

            ! w_j = w_j - a_j * v_j - b_j * v_(j-1)
            w(:) = w(:) - a(j)*v(:,2) - b(j)*v(:,1)

            ! b_(j+1) = norm(w_j)
            b(j+1) = mpi_norm(w(:), basis%nloc)

            if (b(j+1).lt.1e-12) then
                exit lanczos_loop
            endif

            ! v_(j-1) = v_j
            v(:,1) = v(:,2)

            ! v_(j+1) = w_j/b(j+1)
            v(:,2) = w(:)/b(j+1)
        enddo lanczos_loop

        if (nstep_calc.eq.nstep-1) then
            call multiply_H( basis, v(:,2), w(:) )
            a(nstep) = mpi_dot_product(w(:),v(:,2),basis%nloc)
            nstep_calc = nstep
        endif

        if (allocated(aout)) deallocate(aout)
        if (allocated(bout)) deallocate(bout)
        allocate(aout(nstep_calc),bout(nstep_calc))

        aout(1:nstep_calc) = a(1:nstep_calc)
        bout(1:nstep_calc) = b(1:nstep_calc)

    end subroutine lanczos_iteration

end module ed_lanczos
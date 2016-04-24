module ed_lanczos
    use precision
    use ed_config
    use ed_basis
    use ed_utils
    use ed_hamiltonian
    use parallel_params
    use sys
    
    implicit none

contains

    subroutine lanczos_iteration(basis,vec,nstep_calc,aout,bout)
        include 'mpif.h'

        type(basis_t), intent(in) :: basis
        real(dp), intent(in) :: vec(basis%nloc) ! starting vector for lanczos iteration

        integer, intent(out) :: nstep_calc
        real(dp), allocatable, intent(out) :: aout(:), bout(:)

        ! lanczos vectors
        real(dp), allocatable :: a(:), b(:)
        real(dp), allocatable :: v(:,:), w(:)
        real(dp) :: norm_v
        
        integer :: j, ierr
        
        allocate(v(basis%nloc,3),w(basis%nloc))
        allocate(a(0:nstep),b(0:nstep))

        v(1:basis%nloc,2) = 0.0_dp
        w(1:basis%nloc) = 0.0_dp

        a(0:nstep) = 0.0_dp
        b(0:nstep) = 0.0_dp

        v(:,1) = vec
        ! normalize the initial vector
        norm_v = mpi_dot_product( v(:,1), v(:,1), basis%nloc)
        v(:,1) = v(:,1)/sqrt(norm_v)
        call multiply_H(basis,v(:,1),w(:))
        a(0) = mpi_dot_product(v(:,1),w(:),basis%nloc)
        b(0) = 0.0_dp

        v(:,2) = w(:) - a(0)*v(:,1)
        b(1) = mpi_dot_product(v(:,2),v(:,2),basis%nloc)
        call multiply_H(basis,v(:,2),w(:))
        a(1) = mpi_dot_product( w(:), v(:,2), basis%nloc)
        a(1) = a(1)/b(1)
        v(:,2) = v(:,2)/sqrt(b(1))

        nstep_calc = nstep
        lanczos_loop: do j = 2, nstep
            v(:,3) = w(:)/sqrt(b(j-1)) - a(j-1)*v(:,2) - sqrt(b(j-1))*v(:,1)

            b(j) = mpi_dot_product(v(:,3),v(:,3),basis%nloc)
            call multiply_H(basis,v(:,3),w(:))

            a(j) = mpi_dot_product(v(:,3),w(:),basis%nloc)
            a(j) = a(j)/b(j)
            v(:,3) = v(:,3)/sqrt(b(j))
            if(b(j).lt.1.e-12) then
                nstep_calc = j
                if (node.eq.0) then
                    write(*,*) "lanczos: b_(j+1) = 0 at j=", j
                endif
                exit lanczos_loop
            endif
            v(:,1) = v(:,2)
            v(:,2) = v(:,3)
        enddo lanczos_loop

        ! ! Lanczos steps
        ! ! ref: https://en.wikipedia.org/wiki/Lanczos_algorithm#Iteration 

        ! ! normalize the initial vector
        ! norm_v = mpi_norm( v(1:basis%nloc,2), basis%nloc)
        ! v(1:basis%nloc,2) = v(1:basis%nloc,2)/norm_v

        ! ! v(:,1) = v_(j-1)
        ! ! v(:,2) = v_j
        ! ! w(:)   = w_j
        ! lanczos_loop: do j=1,nstep-1
        !     nstep_calc = j

        !     ! w_j = H*v_j
        !     call multiply_H( basis, v(1:basis%nloc,2), w(1:basis%nloc) )

        !     ! a_j = dot(w_j,v_j)
        !     a(j) = mpi_dot_product(w(1:basis%nloc), v(1:basis%nloc,2), basis%nloc)

        !     ! w_j = w_j - a_j * v_j - b_j * v_(j-1)
        !     w(:) = w(:) - a(j)*v(:,2) - b(j)*v(:,1)

        !     ! b_(j+1) = norm(w_j)
        !     b(j+1) = mpi_norm(w(:), basis%nloc)

        !     if (b(j+1).lt.1e-8) then
        !         if (node.eq.0) then
        !             write(*,*) "lanczos: b_(j+1) = 0 at j=", j
        !         endif
        !         exit lanczos_loop
        !     endif

        !     ! v_(j-1) = v_j
        !     v(:,1) = v(:,2)

        !     ! v_(j+1) = w_j/b(j+1)
        !     v(:,2) = w(:)/b(j+1)
        !     print "(I,F20.16)",j, mpi_norm(v(:,2),basis%nloc)
        ! enddo lanczos_loop

        ! if (nstep_calc.eq.nstep-1) then
        !     call multiply_H( basis, v(:,2), w(:) )
        !     a(nstep) = mpi_dot_product(w(:),v(:,2),basis%nloc)
        !     nstep_calc = nstep
        ! endif

        if (allocated(aout)) deallocate(aout)
        if (allocated(bout)) deallocate(bout)
        allocate(aout(0:nstep_calc),bout(0:nstep_calc))

        aout(0:nstep_calc) = a(0:nstep_calc)
        bout(0:nstep_calc) = b(0:nstep_calc)

        deallocate(v,w,a,b)
    end subroutine lanczos_iteration

end module ed_lanczos

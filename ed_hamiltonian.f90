module ed_hamiltonian
    use precision
    use ed_operators
    use ed_config
    use parallel_params
    use ed_utils
    use ed_basis

    implicit none
    include 'mpif.h'

    complex(dp), allocatable :: Hk(:,:,:)
    real(dp), allocatable :: ef(:)

    ! current sector of the hamiltonain.
    integer, private :: isector

    public
contains

    ! Hamiltonian related initialization before entering DMFT loop
    subroutine ed_hamiltonian_init

        allocate(Hk(norb,norb,nq),ef(norb))

    end subroutine ed_hamiltonian_init

    subroutine ed_set_band_structure
        integer nw
        parameter(nw=1000) 
        real(dp):: kx,ky,pi,de,energy
        parameter(pi=acos(-1.D0))
        real(dp):: dq !gb gaussian broadening
        integer:: i,j,ik,nk,iorb

        de = 0.01D0

        Hk(:,:,:) = dcmplx(0.D0,0.D0)
        nk = int(dsqrt(dfloat(Nq)))
        dq = 2.D0*pi/float(nk)

        ik = 0
        do i = 1, nk
            do j = 1, nk
                ik = ik + 1
                kx = dq*float(i-1)
                ky = dq*float(j-1)
                call find_hk(kx,ky,Hk(:,:,ik))
            enddo
        enddo
        if(ik.ne.Nq) stop "ik =\= Nq"

        do i = 1, norb
            ef(i) = sum(real(Hk(i,i,:)))/float(Nq)
            ek(i) = ef(i)
        enddo

        if (Node.eq.0) then

            write(6,*) "ed_hamiltonian: initial impurity levels"
            do i = 1, Norb
                write(6,'(i2,3x,e)') i, ef(i)
            enddo

        endif

    end subroutine ed_set_band_structure

    subroutine find_hk(kx,ky,Hik)
        !      tight-binding parameters from PRL vol. 84, 1591 (2007)
        !      A. Liebsch, A. Lichtenstein 
        implicit none
        integer i,j,k
        real(dp) kx, ky, coskx, cosky
        complex(dp) Hik(Norb,Norb)

        Hik(:,:) = dcmplx(0.0_dp,0.0_dp)

        coskx = cos(kx)
        cosky = cos(ky)

        do i = 1, norb
            Hik(i,i) = -0.5_dp*(coskx + cosky)
        enddo
    end subroutine find_hk

    subroutine multiply_H(basis,A,B)
        type(basis_t), intent(in) :: basis
        real(dp), intent(in) :: A(basis%nloc)
        real(dp), intent(out) :: B(basis%nloc)

        real(dp), allocatable :: A_all(:)
        integer :: i,j,k
        integer :: ispin,jspin,iorb,jorb,ibath
        integer(kind=kind_basis) :: basis_i, basis_j

        real(dp) :: coeff_sum, coeff

        allocate(A_all(basis%ntot))
        call mpi_allgatherv(A,basis%nloc,mpi_double_precision,A_all,&
            basis%nlocals,basis%offsets,mpi_double_precision,comm,ierr)

        B(1:basis%nloc) = 0.0_dp
        iloop: do i=1,basis%nloc
            ! i-th basis state
            basis_i = ed_basis_get(basis,basis%offsets(node)+i)
            
            ! [diagonal elements]
            ! B(i) = H(i,i)*A(i)
            coeff_sum = 0.0_dp
            ! 1. onsite 
            do ispin=1,2
                do iorb=1,Norb
                    call onsite(basis_i,iorb,ispin,basis_j,coeff)
                    coeff_sum = coeff_sum+coeff 
                enddo
                do ibath=1,Nbath
                    call onsite(basis_i,Norb+ibath,ispin,basis_j,coeff)
                    coeff_sum = coeff_sum+coeff
                enddo
            enddo
            ! 2. density-density : intra-orbital
            do iorb=1,Norb
                call dens_dens(basis_i,iorb,iorb,1,2,basis_j,coeff)
                coeff_sum = coeff_sum + coeff 
            enddo
            ! 3. density-density : inter-orbital
            do iorb=1,Norb
                do jorb=iorb+1,Norb
                    do ispin=1,2
                        do jspin=1,2
                            call dens_dens(basis_i,iorb,jorb,ispin,jspin,basis_j,coeff)
                            coeff_sum = coeff_sum + coeff
                        enddo
                    enddo
                enddo
            enddo

            B(i) = coeff_sum*A_all(basis%offsets(node)+i)

            ! [off-diagonal elements]
            ! B(i) = sum_i H(i,j)*A(j)
            ! 1. hybridization
            do ispin=1,2
                do iorb=1,Norb
                    do ibath=1,Nbath
                        call hybridization1(basis_i,iorb,ibath,ispin,basis_j,coeff)
                        if (basis_j.ne.0) then
                            j = ed_basis_idx(basis,basis_j) 
                            B(i) = B(i) + coeff*A_all(j)
                        endif
                        call hybridization2(basis_i,iorb,ibath,ispin,basis_j,coeff)
                        if (basis_j.ne.0) then
                            j = ed_basis_idx(basis,basis_j) 
                            B(i) = B(i) + coeff*A_all(j)
                        endif
                    enddo
                enddo
            enddo
            ! 2. spin-flip
            do iorb=1,Norb
                do jorb=1,Norb
                    if (iorb.eq.jorb) then
                        cycle
                    endif
                    call spin_flip(basis_i,iorb,jorb,basis_j,coeff)
                    if (basis_j.ne.0) then
                        j = ed_basis_idx(basis,basis_j) 
                        B(i) = B(i) + coeff*A_all(j)
                    endif
                enddo
            enddo
            ! 3. pair-exchange
            do iorb=1,Norb
                do jorb=1,Norb
                    if (iorb.eq.jorb) then
                        cycle
                    endif
                    call pair_exchange(basis_i,iorb,jorb,basis_j,coeff)
                    if (basis_j.ne.0) then
                        j = ed_basis_idx(basis,basis_j) 
                        B(i) = B(i) + coeff*A_all(j)
                    endif
                enddo
            enddo
        enddo iloop

        deallocate(A_all)
    end subroutine multiply_H

end module ed_hamiltonian

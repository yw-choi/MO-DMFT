module ed_hamiltonian
    use precision
    use ed_config
    use parallel_params
    use ed_utils
    use alloc
    use ed_basis

    implicit none
    include 'mpif.h'

    complex(dp), pointer :: Hk(:,:,:)
    real(dp), pointer :: ef(:)

    ! current sector of the hamiltonain.
    integer, private :: isector

    public
contains

    ! Hamiltonian related initialization before entering DMFT loop
    subroutine ed_hamiltonian_init
        call re_alloc(Hk,1,norb,1,norb,1,nq,name="Hk",routine="ed_hamiltonian_init")
        call re_alloc(ef,1,norb,name="ef",routine="ed_hamiltonian_init")
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
        enddo

        if (Node.eq.0) then
            write(6,*) " ************************ "
            write(6,*) "     IMPURITY LEVELS      "
            write(6,*) " ************************ "

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
end module ed_hamiltonian

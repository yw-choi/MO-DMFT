module ed_hamiltonian
    use precision
    use ed_config
    use parallel_params
    use ed_utils

    implicit none

    complex(dp), allocatable :: Hk(:,:,:)
    real(dp), allocatable :: ef(:)

    integer :: nbasis_up 
    integer :: nbasis_down
    integer :: nbasis       ! number of basis in the sector.
    integer :: nbasis_loc   ! number of basis in the sector local to the node.

    integer :: isector
    integer, allocatable :: nlocals(:), offsets(:)
contains

    subroutine prepare_sector(isec)
        integer, intent(in) :: isec

        integer :: nd, nam, i

        if (.not.allocated(nlocals)) allocate(nlocals(0:Nodes-1))
        if (.not.allocated(offsets)) allocate(offsets(0:Nodes-1))

        isector = isec

        nd = nelec(isector) - nup(isector)

        nbasis_up = icom(nelec(isector),nup(isector))
        nbasis_down = icom(nelec(isector),nd)
        nbasis = nbasis_up*nbasis_down
        
        nbasis_loc = nbasis/nodes
        nam = mod(nbasis,nodes)
        if (node.lt.nam) nbasis_loc = nbasis_loc + 1

        call mpi_allgather(nbasis_loc,1,mpi_integer,nlocals(0),1,mpi_integer,comm,ierr)
        offsets(0) = 0                ! note : starting from 0
        do i = 1, nodes-1
            offsets(i) = offsets(i-1) + nlocals(i-1)
        enddo

        if (node.eq.0) then
            write(6,*) "[Dimension of the sector]"
        endif
    end subroutine prepare_sector

    subroutine end_sector

    end subroutine end_sector

    subroutine ed_set_band_structure
        integer nw
        parameter(nw=1000) 
        real(dp):: kx,ky,pi,de,energy
        parameter(pi=acos(-1.D0))
        real(dp):: dq !gb gaussian broadening
        real(dp), allocatable::dos(:,:)
        integer:: i,j,ik,nk,iorb

        allocate(Hk(norb,norb,nq),ef(norb))
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

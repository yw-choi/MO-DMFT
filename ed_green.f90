module ed_green
    use precision
    use parallel_params
    use ed_io
    use ed_config
    use ed_basis
    use alloc

    implicit none

    include 'mpif.h'

    integer, allocatable :: nlocals_w(:), offsets_w(:)
    real(dp), allocatable :: omega(:)
    integer :: nwloc

    complex(dp), allocatable :: Gr(:,:)

    integer, parameter :: pi = acos(-1.0_dp)
contains

    subroutine ed_green_init

        integer :: namw,i,ishift

        nwloc = nw/nodes
        namw = mod(nw,nodes)
        if(node.lt.namw) nwloc=nwloc+1

        allocate(nlocals_w(0:nodes-1),offsets_w(0:nodes-1))
        call mpi_allgather(nwloc,1,mpi_integer,nlocals_w(0),1,mpi_integer,comm,ierr)

        offsets_w(0) = 0
        do i = 1, nodes-1
            offsets_w(i) = offsets_w(i-1) + nlocals_w(i-1)
        enddo

        allocate(omega(1:nwloc))
        if(node.lt.namw) then
            ishift = node*nwloc
            do i = 1, nwloc  ! parallelization over frequecies
                omega(i) = (2.D0*float(ishift+i-1)+1)*pi/beta
            enddo
        else
            ishift = namw+nodes*nwloc
            do i = 1, nwloc  ! parallelization over frequecies
                omega(i) = (2.D0*float(ishift+i-1)+1)*pi/beta
            enddo
        endif
        
        allocate(Gr(1:Norb,1:nwloc))
    end subroutine ed_green_init

    subroutine ed_calc_green_ftn(nev_calc)
        integer, intent(in) :: nev_calc
        
        real(dp), allocatable :: eigval(:), pev(:)
        integer, allocatable :: ind(:)

        real(dp) :: Z
        double precision:: no_up(Norb),no_dn(Norb),no_up_tmp(Norb),&
                           nou,nod,nocc(Norb),ap(0:Nstep),bp(0:Nstep),&
                           an(0:Nstep),bn(0:Nstep),factor

        integer :: i, isector, ilevel, nloc, iorb, iev, ispin

        allocate(eigval(nev_calc),pev(nev_calc),ind(nev_calc))
        call import_eigval(nev_calc,eigval,pev,ind)

        Z = 0.D0
        do i = 1, nev_calc
           Z = Z + pev(i)
        enddo

        do i = 1, nev_calc 
           pev(i) = pev(i)/Z
        enddo

        Gr(:,:) = 0.D0 

        no_up(:) = 0.D0
        no_dn(:) = 0.D0

        nevloop: do iev=1,nev_calc
            isector = (ind(iev)-1)/nev+1
            ilevel = mod(ind(iev)-1,nev)+1

            !call prepare_basis_for_sector(isector,nloc)    

            ! @TODO import eigenvector

            do ispin=1,2
                do iorb = 1,norb
                    ! call calc_green_diag(nloc,eigvec,eigval(iev),pev(iev), &
                    !                 iorb, ispin)
                enddo
            enddo

        enddo nevloop

    end subroutine ed_calc_green_ftn

    ! subroutine calc_green_diag(nloc,eigvec,eigval,pev,iorb,ispin)
    !     integer, intent(in) :: nloc, iorb, ispin
    !     real(dp), intent(in) :: eigvec(nloc), eigval, pev


    ! end subroutine green_diag

    ! subroutine lanczos(
    !     real(dp), allocatable, intent(out) :: 
    !    integer :: istep

    !    do istep=1,nstep
            
    !    enddo

    ! end subroutine lanczos
end module ed_green

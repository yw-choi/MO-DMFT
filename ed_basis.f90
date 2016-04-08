module ed_basis
    use ed_config
    use parallel_params
    use ed_utils
    use alloc
    use precision
    use ionew
    
    implicit none

    integer, pointer :: nlocals(:), offsets(:)
    integer, pointer :: nbasis_up(:)
    integer, pointer :: nbasis_down(:)
    integer, pointer :: nbasis(:)       
    integer :: nbasis_loc ! number of basis in the sector local to the node.

    ! For up,down basis, 4-bit integer is sufficient,
    ! because the maximum number of sites will not be more than 31.
    integer, pointer :: basis_up(:)
    integer, pointer :: basis_down(:)

    integer, pointer :: basis_up_idx(:)
    integer, pointer :: basis_down_idx(:)

    integer, private :: isector

    public
contains

    subroutine ed_basis_init
        integer :: i, nd

        call re_alloc(nlocals,0,Nodes-1,name="nlocals",routine="ed_hamiltonian_init")
        call re_alloc(offsets,0,Nodes-1,name="offsets",routine="ed_hamiltonian_init")

        call re_alloc(nbasis_up,1,nsector,name="nbasis_up",routine="ed_hamiltonian_init")
        call re_alloc(nbasis_down,1,nsector,name="nbasis_down",routine="ed_hamiltonian_init")
        call re_alloc(nbasis,1,nsector,name="nbasis",routine="ed_hamiltonian_init")

        do i=1,nsector
            nd = nelec(i) - nup(i)
            nbasis_up(i) = icom(nelec(i),nup(i))
            nbasis_down(i) = icom(nelec(i),nd)
            nbasis(i) = nbasis_up(i)*nbasis_down(i)
        enddo

        if (node.eq.0) then
            write(6,"(A)") "|------------------------------------------------------------------------|"
            write(6,"(A)") "  Dimension of each sector" 
            write(6,"(A)") "|------------------------------------------------------------------------|"
            write(6,"(A)") "| isector | Ne      | Nup     | Ndown   | nbasis   | nbasis_u | nbasis_d |"
            write(6,"(A)") "|------------------------------------------------------------------------|"
            do i=1,nsector
                write(6,"(7(I9,x))") i,Nelec(i),Nup(i),(Nelec(i)-nup(i)),nbasis(i),nbasis_up(i),nbasis_down(i)
            enddo
            write(6,"(A)") "|------------------------------------------------------------------------|"
        endif
    end subroutine ed_basis_init

    subroutine prepare_basis_for_sector(isec)
        include 'mpif.h'
        integer, intent(in) :: isec
        integer :: nam, i
        isector = isec

        nbasis_loc = nbasis(isector)/nodes
        nam = mod(nbasis(isector),nodes)
        if (node.lt.nam) nbasis_loc = nbasis_loc + 1

        call mpi_allgather(nbasis_loc,1,mpi_integer,nlocals(0),1,mpi_integer,comm,ierr)

        offsets(0) = 0 
        do i = 1, nodes-1
            offsets(i) = offsets(i-1) + nlocals(i-1)
        enddo

        call generate_basis(isector)
    end subroutine prepare_basis_for_sector

    subroutine generate_basis(isector)
        integer, intent(in) :: isector
        integer :: ispin, i, j
        integer :: minrange, maxrange, counts, nbit
        integer :: nud(2)

        call re_alloc(basis_up,1,nbasis_up(isector),'basis_up','prepare_basis_for_sector')
        call re_alloc(basis_down,1,nbasis_down(isector),'basis_down','prepare_basis_for_sector')

        nud(1) = nup(isector)
        nud(2) = nelec(isector) - nup(isector)

        do ispin=1,2
            minrange = 0
            maxrange = 0

            do i=1,nud(ispin)
                minrange = minrange + 2**(i-1)
                maxrange = maxrange + 2**(Nsite-i)
            enddo
            
            if (ispin.eq.1) then
                call re_alloc(basis_up_idx,minrange,maxrange,'basis_up_idx','generate_basis')
            else
                call re_alloc(basis_down_idx,minrange,maxrange,'basis_down_idx', &
                              'generate_basis')
            endif
            
            counts = 0

            do i=minrange,maxrange
                nbit = 0
                do j=0,Nsite-1
                    if (BTEST(i,j)) then
                        nbit = nbit + 1
                    endif
                enddo

                if (nbit.eq.nud(ispin)) then
                    counts = counts + 1
                    if (ispin.eq.1) then
                        basis_up(counts) = i
                        basis_up_idx(i) = counts
                    else
                        basis_down(counts) = i
                        basis_down_idx(i) = counts
                    endif
                endif
            enddo
        enddo

    end subroutine generate_basis

    integer(kind=kind_basis) function ed_basis_get(ib_g)
        integer :: ib_g
        integer :: iup, idown

        iup = mod(ib_g-1,nbasis_up(isector))+1
        idown = (ib_g-1)/nbasis_up(isector)+1

        ed_basis_get=basis_up(iup)+2**(Nsite)*basis_down(idown)
    end function ed_basis_get

    integer function get_basis_idx(basis)
        integer(kind=kind_basis) :: basis
        integer :: bup, bdown

        bup = mod(basis,2**Nsite)
        bdown = basis/(2**Nsite)

        get_basis_idx = basis_up_idx(bup) + &
                        (basis_down_idx(bdown)-1)*nbasis_up(isector)

    end function get_basis_idx

#ifdef DEBUG
    subroutine dump_basis
        integer :: i,iunit
        character(len=100) :: fname

        call io_assign(iunit)
        write(fname,"(A,I1)") "basis.node.",node
        open(unit=iunit,file=fname,status="replace")
        do i=1,nbasis_loc
            write(iunit,"(I5,B)") i,ed_basis_get(offsets(node)+i)
        enddo
        close(iunit)
    end subroutine dump_basis
#endif

end module ed_basis

module ed_basis
    use ed_config
    use parallel_params
    use ed_utils
    use alloc
    use precision
    use ionew
    
    implicit none

    type basis_t
        integer :: nloc

        integer :: ntot
        integer :: nup
        integer :: ndown

        ! For up,down basis, 4-bit integer is sufficient,
        ! because the maximum number of sites will not be more than 31.
        integer, pointer :: up(:)
        integer, pointer :: down(:)

        integer, pointer :: idx_up(:)
        integer, pointer :: idx_down(:)

        integer, allocatable :: nlocals(:)
        integer, allocatable :: offsets(:)
    end type basis_t

    integer, private :: isector

    public
contains

    subroutine prepare_basis( ne_up, ne_down, basis )
        include 'mpif.h'
        integer, intent(in) :: ne_up, ne_down
        type(basis_t), intent(out) :: basis

        integer :: nam, i

        basis%nbasis_up   = icom(Nsite,ne_up)
        basis%nbasis_down = icom(Nsite,ne_down)
        basis%nbasis      = basis%nbasis_up * basis%nbasis_down

        basis%nbasis_loc  = basis%nbasis/nodes
        nam = mod(basis%nbasis,nodes)
        if (node.lt.nam) basis%nbasis_loc = basis%nbasis_loc + 1

        call mpi_allgather(nbasis%nbasis_loc,1,mpi_integer,basis(0),1,mpi_integer,comm,ierr)

        offsets(0) = 0 
        do i = 1, nodes-1
            offsets(i) = offsets(i-1) + nlocals(i-1)
        enddo

        call re_alloc(basis_up,1,nbasis_up(isector),'basis_up','prepare_basis_for_sector')
        call re_alloc(basis_down,1,nbasis_down(isector),'basis_down','prepare_basis_for_sector')
        call generate_basis(nup(isector),nelec(isector)-nup(isector),nbasis_up(isector),&
                            nbasis_down(isector),basis_up,basis_down,basis_up_idx,basis_down_idx)
    end subroutine prepare_basis_for_sector

    subroutine generate_basis(neup,nedown,nbup,nbdown,b_up,b_down,bidx_up,bidx_down)
        integer, intent(in) :: neup, nedown, nbup, nbdown
        integer, pointer, intent(out) :: b_up(:), b_down(:), bidx_up(:), bidx_down(:)

        integer :: ispin, i, j
        integer :: minrange, maxrange, counts, nbit
        integer :: nud(2)

        nud(1) = neup
        nud(2) = nedown

        do ispin=1,2
            minrange = 0
            maxrange = 0

            do i=1,nud(ispin)
                minrange = minrange + 2**(i-1)
                maxrange = maxrange + 2**(Nsite-i)
            enddo
            
            if (ispin.eq.1) then
                call re_alloc(bidx_up,minrange,maxrange,routine='generate_basis')
            else
                call re_alloc(bidx_down,minrange,maxrange,routine='generate_basis')
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
                        b_up(counts) = i
                        bidx_up(i) = counts
                    else
                        b_down(counts) = i
                        bidx_down(i) = counts
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
    subroutine dump_basis(nbasis_loc)
        integer :: nbasis_loc
        integer :: i,iunit,j
        integer(kind=kind_basis) :: basis
        character(len=100) :: fname

        call io_assign(iunit)
        write(fname,"(A,I1)") "basis.node.",node
        open(unit=iunit,file=fname,status="replace")
        do i=1,nbasis_loc
            write(iunit,"(I10,3x)",advance="no") i
            basis = ed_basis_get(offsets(node)+i)
            do j=1,2*Nsite
                if (BTEST(basis,j-1)) then
                    write(iunit,"(I1)",advance="no") 1
                else
                    write(iunit,"(I1)",advance="no") 0
                endif
            enddo
            write(iunit,*)
        enddo
        close(iunit)
    end subroutine dump_basis
#endif

end module ed_basis

module ed_basis
    use ed_config, only: nsite, kind_basis
    use parallel_params, only: nodes, node,comm,ierr
    use sys
    use ed_utils, only: icom
    use precision
    
    implicit none

    type basis_t
        integer :: nloc

        integer :: ntot
        integer :: nup
        integer :: ndown

        integer :: ne_up
        integer :: ne_down

        ! For up,down basis, 4-bit integer is sufficient,
        ! because the maximum number of sites will not be more than 31.
        integer, allocatable :: up(:)
        integer, allocatable :: down(:)

        integer, allocatable :: idx_up(:)
        integer, allocatable :: idx_down(:)

        integer, allocatable :: nlocals(:)
        integer, allocatable :: offsets(:)
    end type basis_t

    public
contains

    type(basis_t) function generate_basis( ne_up, ne_down ) result(basis)
        include 'mpif.h'

        integer, intent(in) :: ne_up, ne_down

        ! local variables
        integer :: ispin, i, j, nam
        integer :: minrange, maxrange, counts, nbit
        integer :: nud(2)

        basis%ne_up = ne_up
        basis%ne_down = ne_down

        basis%nup   = icom(Nsite,ne_up)
        basis%ndown = icom(Nsite,ne_down)
        basis%ntot  = basis%nup * basis%ndown

        basis%nloc  = basis%ntot/nodes
        nam = mod(basis%ntot,nodes)
        if (node.lt.nam) basis%nloc = basis%nloc + 1

        allocate(basis%nlocals(0:nodes-1),basis%offsets(0:nodes-1))
        call mpi_allgather(basis%nloc,1,mpi_integer,basis%nlocals(0),1,mpi_integer,comm,ierr)
        if (ierr.ne.0) then
            call die("MPIError in generate_basis")
            return
        endif
        call mpi_barrier(mpi_comm_world,ierr)

        basis%offsets(0) = 0 
        do i = 1, nodes-1
            basis%offsets(i) = basis%offsets(i-1) + basis%nlocals(i-1)
        enddo

        allocate(basis%up(basis%nup), basis%down(basis%ndown))

        nud(1) = ne_up
        nud(2) = ne_down

        do ispin=1,2
            ! ref : arXiv:1307.7542, Appendix A
            minrange = 0
            maxrange = 0

            do i=1,nud(ispin)
                minrange = minrange + 2**(i-1)
                maxrange = maxrange + 2**(Nsite-i)
            enddo
            
            if (ispin.eq.1) then
                allocate(basis%idx_up(minrange:maxrange))
            else
                allocate(basis%idx_down(minrange:maxrange))
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
                        basis%up(counts) = i
                        basis%idx_up(i) = counts
                    else
                        basis%down(counts) = i
                        basis%idx_down(i) = counts
                    endif
                endif
            enddo
        enddo
    end function generate_basis

    ! ref : arXiv:1307.7542 eq (6)
    integer(kind=kind_basis) function ed_basis_get(basis,idx_loc) 
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: idx_loc

        ! local variables
        integer :: iup, idown, idx

        ! local idx to global idx
        idx = idx_loc + basis%offsets(node)

        iup = mod(idx-1,basis%nup)+1
        idown = (idx-1)/basis%nup+1

        ed_basis_get = basis%up(iup)+2**(Nsite)*basis%down(idown)
    end function ed_basis_get

    ! ref : arXiv:1307.7542 eq (8)
    integer function ed_basis_idx(basis, basis_i)
        type(basis_t), intent(in) :: basis
        integer(kind=kind_basis) :: basis_i

        ! local variables
        integer :: basis_i_up, basis_i_down
        
        basis_i_up = mod(basis_i,2**Nsite)
        basis_i_down = basis_i/(2**Nsite)

        ed_basis_idx = (basis%idx_down(basis_i_down)-1)*basis%nup + &
                             basis%idx_up(basis_i_up)
    end function ed_basis_idx

end module ed_basis

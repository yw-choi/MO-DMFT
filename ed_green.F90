module ed_green
    use precision
    use parallel_params
    use ed_io
    use ed_config
    use ed_basis
    use ed_operators
    use ed_lanczos
    use ed_utils
    use sys

    implicit none

    integer, allocatable, save :: nlocals_w(:), offsets_w(:)
    real(dp), allocatable, save :: omega(:)
    integer, save :: nwloc

    complex(dp), allocatable, save :: Gr(:,:), D_ev(:,:), Gksum(:,:)
    real(dp), allocatable, save :: nocc(:)

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
        
        allocate(Gr(Norb,nwloc),D_ev(Norb,nwloc),Gksum(norb,nwloc))
        allocate(nocc(norb))
        ! allocate(nocc(norb,2))
    end subroutine ed_green_init

    subroutine ed_calc_green_ftn(nev_calc)
        integer, intent(in) :: nev_calc
        
        real(dp), allocatable :: eigval(:), pev(:)
        real(dP), allocatable :: eigvec(:)
        integer, allocatable :: ind(:)
        real(dp) :: Z, factor, nocc_up, nocc_down
        integer :: i, isector, ilevel, nloc, iorb, iev, ispin, nd
        type(basis_t) :: basis

        allocate(eigval(nev_calc),pev(nev_calc),ind(nev_calc))

        if (node.eq.0) then
            write(6,*) "ed_green: import eigvalues"
        endif
        call import_eigval(nev_calc,eigval,pev,ind)

        Z = 0.D0
        do i = 1, nev_calc
           Z = Z + pev(i)
        enddo

        do i = 1, nev_calc 
           pev(i) = pev(i)/Z
        enddo

        Gr(:,:) = 0.D0 

        ! nocc(:,:) = 0.D0

        nevloop: do iev=1,nev_calc
            isector = (ind(iev)-1)/nev+1
            ilevel = mod(ind(iev)-1,nev)+1

            if (node.eq.0) then
                write(6,*) "ed_green: iev=",iev,", isector=",isector,",ilevel=",ilevel
            endif

            nd = nelec(isector)-nup(isector)
            if (nd.eq.nup(isector)) then
                factor=1.0_dp
            else
                factor=0.5_dp
            endif

            call import_eigvec(isector,ilevel,node,nloc,eigvec)
            basis = generate_basis( nup(isector), nelec(isector)-nup(isector))

            if (nloc.ne.basis%nloc) then
                write(6,*) "ed_calc_green_ftn: basis dimension mismatch"
                write(6,*) "node=",node," nloc_read=",nloc," nloc=",basis%nloc
                call die
                return 
            endif

            ! spin up part 
            do iorb = 1,norb
                call green_diag(basis,eigvec,eigval(iev),pev(iev),factor,iorb, 1)
            enddo

            if (nd.ne.nup(isector)) then
                ! spin down part
                do iorb = 1,norb
                    call green_diag(basis,eigvec,eigval(iev),pev(iev),factor,iorb, 2)
                enddo

            endif
        enddo nevloop
    end subroutine ed_calc_green_ftn

    subroutine green_diag(basis,eigvec,eigval,pev,factor,iorb,ispin)
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: iorb, ispin
        real(dp), intent(in) :: eigvec(basis%nloc), eigval, pev, factor

        ! local variables
        integer :: i, j
        type(basis_t) :: basis_out
        real(dp), allocatable :: v(:)
        real(dp) :: nocc_i
        integer :: nstep_calc
        complex(dp) :: cmplx_omega, grx

        ! lanczos matrix elements
        real(dp), allocatable :: a(:), b(:) 

        ! One more particle at iorb,ispin
        call apply_c(basis, eigvec, 1, iorb, ispin, basis_out, v)
        call lanczos_iteration(basis_out, v, nstep_calc, a, b)

        nocc_i = mpi_dot_product(v,v,basis_out%nloc)
        b(0) = 1-nocc_i
        do i=1,nwloc
            cmplx_omega = cmplx(eigval,omega(i))
            grx = b(nstep_calc)/(cmplx_omega-a(nstep_calc))
            do j = nstep_calc-1, 0, -1
                grx = b(j)/(cmplx_omega-a(j)-grx)
            enddo 

            Gr(iorb,i) = Gr(iorb,i) + pev*factor*grx
        enddo

        ! One less particle at iorb,ispin
        call apply_c(basis, eigvec, 2, iorb, ispin, basis_out, v)
        call lanczos_iteration(basis_out, v, nstep_calc, a, b)

        b(0) = nocc_i
        do i=1,nwloc
            cmplx_omega = cmplx(-eigval,omega(i))
            grx = b(nstep_calc)/(cmplx_omega+a(nstep_calc))
            do j = nstep_calc-1, 0, -1
                grx = b(j)/(cmplx_omega+a(j)-grx)
            enddo 

            Gr(iorb,i) = Gr(iorb,i) + pev*factor*grx
        enddo
    end subroutine green_diag

    subroutine find_nocc(basis,vec,iorb,no_up,no_dn)
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: iorb
        real(dp), intent(in) :: vec(basis%nloc)

        real(dp), intent(out) :: no_up, no_dn

        integer :: i
        real(dp) :: ind_up(basis%nloc), ind_dn(basis%nloc)
        real(dp) :: no_up_tmp, no_dn_tmp
        integer(kind=kind_basis) :: basis_i

        ind_up(:) = 0.D0
        ind_dn(:) = 0.D0
        do i = 1, basis%nloc
            basis_i = ed_basis_get(basis,i)
            if(BTEST(basis_i,iorb)) ind_up(i) = 1.D0
            if(BTEST(basis_i,Nsite+iorb)) ind_dn(i) = 1.D0
        enddo

        no_up_tmp = sum(ind_up(:)*vec(:)*vec(:))
        no_dn_tmp = sum(ind_dn(:)*vec(:)*vec(:))

        call mpi_allreduce(no_up_tmp,no_up,1,mpi_double_precision,mpi_sum,&
            comm,ierr)
        call mpi_allreduce(no_dn_tmp,no_dn,1,mpi_double_precision,mpi_sum,&
            comm,ierr)

        return
    end subroutine find_nocc

    subroutine ed_delta_new

        complex(dp) :: Atmp(Norb,Norb),Btmp(Norb,Norb),Hk(Norb,Norb,Nq)
        complex(dp) :: Ptmp(Norb,Norb)
        integer :: iw, jq, iorb, korb

        do iw = 1,nwloc
            Gksum(:,iw) = cmplx(0.D0,0.D0)
            do jq = 1, Nq
                do iorb = 1, Norb
                    do korb = 1, Norb
                        Ptmp(iorb,korb)=(1.D0/Gr(iorb,iw)+D_ev(iorb,iw) &
                            +rMu)*kdel(iorb,korb)-Hk(iorb,korb,jq)
                    enddo
                enddo
                Atmp(:,:) = Ptmp(:,:)
                call cinv(Atmp,Norb,Norb,Btmp)
                do iorb = 1, Norb
                    Gksum(iorb,iw) = Gksum(iorb,iw)+Btmp(iorb,iorb)
                enddo
            enddo
            Gksum(:,iw) = Gksum(:,iw)/float(Nq)
            D_ev(:,iw) = 1.D0/Gr(:,iw)-1.D0/Gksum(:,iw) + D_ev(:,iw)
        enddo
    end subroutine ed_delta_new

    real(dp) function kdel(iorb,korb)
        integer:: iorb, korb

        if(iorb.eq.korb) kdel = 1.D0
        if(iorb.ne.korb) kdel = 0.D0

        return
    end function kdel

    subroutine n_from_gksum
        integer i,j,k,iw,iorb,nw,master
        complex(dp), allocatable:: Gksum_tot(:,:)

        call mpi_allreduce(nwloc,nw,1,mpi_integer,mpi_sum,comm,ierr)

        if (node.eq.0) then
            allocate(Gksum_tot(norb,nw))
        endif

        do iorb = 1, norb
            call mpi_gatherv(gksum(iorb,:),nwloc,mpi_double_complex,&
                             Gksum_tot(iorb,:),nlocals_w,offsets_w,mpi_double_complex,&
                             0,comm,ierr)
        enddo

        if(node.eq.0) then
            write(6,*) "Occupations from local Green function"
            do iorb = 1, norb
                nocc(iorb) = 2.D0*sum(real(gksum_tot(iorb,:)))/beta+0.5D0
                write(6,'(a,i2,3x,f10.7)') "Orbtal",iorb,nocc(iorb)
            enddo
            write(6,'(a,5x,f10.7)') "sum =", sum(nocc(:))
            deallocate(gksum_tot)
        endif
    end subroutine n_from_gksum
end module ed_green

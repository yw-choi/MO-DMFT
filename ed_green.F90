module ed_green

    use parallel_params
    use ed_io
    use ed_config
    use ed_basis
    use ed_operators
    use ed_lanczos
    use ed_utils
    use ed_solver
    use sys

    implicit none

    integer, allocatable, save :: nlocals_w(:), offsets_w(:)
    double precision, allocatable, save :: omega(:)
    integer, save :: nwloc

    double complex, allocatable, save :: Gr(:,:), Gr_prev(:,:), D_ev(:,:), Gksum(:,:)
    double precision, allocatable, save :: nocc(:)

    integer, parameter :: IO_AB = 89
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
                omega(i) = (2.0D0*float(ishift+i-1)+1)*pi/beta
            enddo
        else
            ishift = namw+node*nwloc
            do i = 1, nwloc  ! parallelization over frequecies
                omega(i) = (2.0D0*float(ishift+i-1)+1)*pi/beta
            enddo
        endif

        allocate(Gr(Norb,nwloc),Gr_prev(Norb,nwloc),D_ev(Norb,nwloc),Gksum(norb,nwloc))
        allocate(nocc(norb))

        call calc_dev
    end subroutine ed_green_init

    subroutine ed_calc_green_ftn(nev_calc)
        integer, intent(in) :: nev_calc

        double precision, allocatable :: eigvec(:)
        double precision :: Z, factor, ap(0:nstep), bp(0:nstep), an(0:nstep), bn(0:nstep)
        integer :: i, nloc, iorb, iev, ispin, nd, isector, ilevel, nstep_calcp, nstep_calcn
        logical even, ex
        type(basis_t) :: basis

        if(node.eq.0) then
            inquire(file="apbpanbn.dat",exist=ex)
            if (ex) then
                open(IO_AB,file="apbpanbn.dat",form="unformatted")
                close(IO_AB,status="delete")
            endif
            open(unit=IO_AB,file="apbpanbn.dat",form="unformatted",status="new")
        endif

        Gr(:,:) = 0.0D0

        nevloop: do iev=1,nev_calc
            if (node.eq.0) then
                write(6,"(x,a,I3)") "ed_green: calculating Green's function of iev = ",iev
            endif
            isector = eigval(iev)%sector
            ilevel  = eigval(iev)%level

            nd = nelec(isector)-nup(isector)
            if (nd.eq.nup(isector)) then
                factor=1.0D0
                even = .true.
            else
                factor=0.5D0
                even = .false.
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
                call green_diag(basis,eigvec,eigval(iev)%val,eigval(iev)%prob,factor,&
                                iorb,1,nstep_calcp,nstep_calcn,ap,bp,an,bn)

                if (node.eq.0) then
                    write(IO_AB) nstep_calcp, nstep_calcn
                    write(IO_AB) iev,iorb,eigval(iev)%val,eigval(iev)%prob,even,&
                                 (ap(i),i=0,nstep_calcp),(bp(i),i=0,nstep_calcp),&
                                 (an(i),i=0,nstep_calcn),(bn(i),i=0,nstep_calcn)
                endif
            enddo

            if (nd.ne.nup(isector)) then
                ! spin down part
                do iorb = 1,norb
                    call green_diag(basis,eigvec,eigval(iev)%val,eigval(iev)%prob,factor,&
                                   iorb,2,nstep_calcp,nstep_calcn,ap,bp,an,bn)
                    if (node.eq.0) then
                        write(IO_AB) nstep_calcp, nstep_calcn
                        write(IO_AB) iev,iorb,eigval(iev)%val,eigval(iev)%prob,even,&
                                     (ap(i),i=0,nstep_calcp),(bp(i),i=0,nstep_calcp),&
                                     (an(i),i=0,nstep_calcn),(bn(i),i=0,nstep_calcn)
                    endif
                enddo

            endif
        enddo nevloop
    end subroutine ed_calc_green_ftn

    subroutine green_diag(basis,eigvec,eigval,pev,factor,iorb,ispin,nstep_calcp,nstep_calcn,ap,bp,an,bn)
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: iorb, ispin
        double precision, intent(in) :: eigvec(basis%nloc), eigval, pev, factor
        double precision, intent(out) :: ap(0:nstep), bp(0:nstep), an(0:nstep), bn(0:nstep)
        integer, intent(out) :: nstep_calcp, nstep_calcn

        ! local variables
        integer :: i, j
        type(basis_t) :: basis_out
        double precision, allocatable :: v(:)
        double precision :: normsq
        double complex :: cmplx_omega, grx

        ! lanczos matrix elements
        double precision, allocatable :: a(:), b(:)

        ! One more particle at iorb,ispin
        call apply_c(basis, eigvec, 1, iorb, ispin, basis_out, v)
        call lanczos_iteration(basis_out, v, nstep_calcp, ap, bp)

        normsq = mpi_dot_product(v,v,basis_out%nloc)
        bp(0) = normsq

        do i=1,nwloc
            cmplx_omega = cmplx(eigval,omega(i))
            grx = bp(nstep_calcp)/(cmplx_omega-ap(nstep_calcp))
            do j = nstep_calcp-1, 0, -1
                grx = bp(j)/(cmplx_omega-ap(j)-grx)
            enddo

            Gr(iorb,i) = Gr(iorb,i) + pev*factor*grx
        enddo

        ! One less particle at iorb,ispin
        call apply_c(basis, eigvec, 2, iorb, ispin, basis_out, v)
        call lanczos_iteration(basis_out, v, nstep_calcn, an, bn)

        bn(0) = 1.0D0-normsq

        do i=1,nwloc
            cmplx_omega = cmplx(-eigval,omega(i))
            grx = bn(nstep_calcn)/(cmplx_omega+an(nstep_calcn))
            do j = nstep_calcn-1, 0, -1
                grx = bn(j)/(cmplx_omega+an(j)-grx)
            enddo

            Gr(iorb,i) = Gr(iorb,i) + pev*factor*grx
        enddo
    end subroutine green_diag

    subroutine ed_delta_new
        ! double complex :: Atmp(Norb,Norb),Btmp(Norb,Norb)
        double complex :: tmp
        integer :: iw, jq, iorb, korb

        if (node.eq.0) then
            write(6,*) "ed_delta_new: Calculating delta_new..."
        endif

        ! @TODO H(k) is assumed to be diagonal
        do iw = 1,nwloc
            Gksum(:,iw) = cmplx(0.0D0,0.0D0)
            do jq = 1, Nq
                ! do iorb = 1, Norb
                !     do korb = 1, Norb
                !         Atmp(iorb,korb)=(1.0D0/Gr(iorb,iw)+D_ev(iorb,iw) &
                !             +rMu)*kdel(iorb,korb)-Hk(iorb,korb,jq)
                !     enddo
                ! enddo
                ! call cinv(Atmp,Norb,Norb,Btmp)
                do iorb = 1, Norb
                    tmp = (1.0D0/Gr(iorb,iw)+D_ev(iorb,iw)+rMu)-Hk(iorb,iorb,jq)
                    Gksum(iorb,iw) = Gksum(iorb,iw)+ 1.0D0/tmp
                enddo
            enddo
            Gksum(:,iw) = Gksum(:,iw)/float(Nq)
            D_ev(:,iw) = 1.0D0/Gr(:,iw)-1.0D0/Gksum(:,iw) + D_ev(:,iw)
        enddo
    end subroutine ed_delta_new

    double precision function kdel(iorb,korb)
        integer:: iorb, korb

        if(iorb.eq.korb) kdel = 1.0D0
        if(iorb.ne.korb) kdel = 0.0D0

        return
    end function kdel

    subroutine n_from_gksum
        integer i,j,k,iw,iorb,nw
        double complex, allocatable:: Gksum_tot(:,:)

        call mpi_allreduce(nwloc,nw,1,mpi_integer,mpi_sum,comm,ierr)

        allocate(Gksum_tot(norb,nw))

        do iorb = 1, norb
            call mpi_gatherv(gksum(iorb,:),nwloc,mpi_double_complex,&
                Gksum_tot(iorb,1:nw),nlocals_w,offsets_w,mpi_double_complex,&
                0,comm,ierr)
        enddo

        if (node.eq.0) then
            write(6,*)
            write(6,*) "Occupations from local Green function"
            write(6,*)
            write(6,"(a)") " Orbital    Occupancy"
            do iorb = 1, norb
                nocc(iorb) = 2.0D0*sum(real(gksum_tot(iorb,:)))/beta+0.50D0
                write(6,'(x,i7,4x,f9.6)') iorb,nocc(iorb)
            enddo
            write(6,'(a,f9.6)') "   Total    ", sum(nocc(:))
            write(6,*)
        endif
        deallocate(gksum_tot)
    end subroutine n_from_gksum

    !############################# DEFINITION OF D_EV ###############################
    !
    !         D_ev(io,iw) = \sum_k vk(ko)*vk(ko)/(iw-ek) + ef(ko)
    !
    !################################################################################
    subroutine calc_dev
        integer ::i,korb,iorb,k

        do i = 1,nwloc
            do korb = 1, Norb
                D_ev(korb,i) = cmplx(0.0D0,0.0D0)
                do k = norb+1,nsite
                    D_ev(korb,i) = D_ev(korb,i) + &
                        vk(korb,k)*vk(korb,k)/(cmplx(0.0D0,omega(i))-ek(k))
                enddo
                D_ev(korb,i) = D_ev(korb,i) + ef(korb)
            enddo
        enddo
    end subroutine calc_dev
end module ed_green

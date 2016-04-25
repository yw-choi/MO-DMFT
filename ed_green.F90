module ed_green
    use precision
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
    real(dp), allocatable, save :: omega(:)
    integer, save :: nwloc

    complex(dp), allocatable, save :: Gr(:,:), Gr_prev(:,:), D_ev(:,:), Gksum(:,:)
    real(dp), allocatable, save :: nocc(:)

    integer, parameter :: IO_AB = 109
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
                omega(i) = (2.0_dp*float(ishift+i-1)+1)*pi/beta
            enddo
        else
            ishift = namw+node*nwloc
            do i = 1, nwloc  ! parallelization over frequecies
                omega(i) = (2.0_dp*float(ishift+i-1)+1)*pi/beta
            enddo
        endif
        
        allocate(Gr(Norb,nwloc),Gr_prev(Norb,nwloc),D_ev(Norb,nwloc),Gksum(norb,nwloc))
        allocate(nocc(norb))

        call calc_dev
    end subroutine ed_green_init

    subroutine ed_calc_green_ftn(nev_calc)
        integer, intent(in) :: nev_calc
        
        real(dP), allocatable :: eigvec(:)
        real(dp) :: Z, factor, ap(0:nstep), bp(0:nstep), an(0:nstep), bn(0:nstep)
        integer :: i, nloc, iorb, iev, ispin, nd, isector, ilevel
        logical even
        type(basis_t) :: basis

        if(node.eq.0) open(unit=IO_AB,file="apbpanbn.dat",form="unformatted")

        Gr(:,:) = 0.0_dp 

        nevloop: do iev=1,nev_calc
            if (node.eq.0) then
                write(6,"(x,a,I3)") "ed_green: calculating Green's function of iev = ",iev
            endif
            isector = eigval(iev)%sector
            ilevel  = eigval(iev)%level

            nd = nelec(isector)-nup(isector)
            if (nd.eq.nup(isector)) then
                factor=1.0_dp
                even = .true.
            else
                factor=0.5_dp
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
                                iorb,1,ap,bp,an,bn)
                 write(IO_AB) iev,iorb,eigval(iev)%val,eigval(iev)%prob,even,&
                           (ap(i),i=0,Nstep),(bp(i),i=0,Nstep),&
                           (an(i),i=0,Nstep),(bn(i),i=0,Nstep)
            enddo

            if (nd.ne.nup(isector)) then
                ! spin down part
                do iorb = 1,norb
                    call green_diag(basis,eigvec,eigval(iev)%val,eigval(iev)%prob,factor&
                                    ,iorb,2,ap,bp,an,bn)
                    write(IO_AB) iev,iorb,eigval(iev)%val,eigval(iev)%prob,even,&
                                (ap(i),i=0,Nstep),(bp(i),i=0,Nstep),&
                                (an(i),i=0,Nstep),(bn(i),i=0,Nstep)
                enddo

            endif
        enddo nevloop
    end subroutine ed_calc_green_ftn

    subroutine green_diag(basis,eigvec,eigval,pev,factor,iorb,ispin,ap,bp,an,bn)
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: iorb, ispin
        real(dp), intent(in) :: eigvec(basis%nloc), eigval, pev, factor
        real(dp), intent(out) :: ap(0:nstep), bp(0:nstep), an(0:nstep), bn(0:nstep)

        ! local variables
        integer :: i, j
        type(basis_t) :: basis_out
        real(dp), allocatable :: v(:)
        real(dp) :: normsq
        integer :: nstep_calc
        complex(dp) :: cmplx_omega, grx

        ! lanczos matrix elements
        real(dp), allocatable :: a(:), b(:) 

        ! One more particle at iorb,ispin
        call apply_c(basis, eigvec, 1, iorb, ispin, basis_out, v)
        call lanczos_iteration(basis_out, v, nstep_calc, ap, bp)

        normsq = mpi_dot_product(v,v,basis_out%nloc)
        bp(0) = normsq

        do i=1,nwloc
            cmplx_omega = cmplx(eigval,omega(i))
            grx = bp(nstep_calc)/(cmplx_omega-ap(nstep_calc))
            do j = nstep_calc-1, 0, -1
                grx = bp(j)/(cmplx_omega-ap(j)-grx)
            enddo 

            Gr(iorb,i) = Gr(iorb,i) + pev*factor*grx
        enddo

        ! One less particle at iorb,ispin
        call apply_c(basis, eigvec, 2, iorb, ispin, basis_out, v)
        call lanczos_iteration(basis_out, v, nstep_calc, an, bn)

        bn(0) = 1.0_dp-normsq

        do i=1,nwloc
            cmplx_omega = cmplx(-eigval,omega(i))
            grx = bn(nstep_calc)/(cmplx_omega+an(nstep_calc))
            do j = nstep_calc-1, 0, -1
                grx = bn(j)/(cmplx_omega+an(j)-grx)
            enddo 

            Gr(iorb,i) = Gr(iorb,i) + pev*factor*grx
        enddo
    end subroutine green_diag

    subroutine ed_delta_new
        ! complex(dp) :: Atmp(Norb,Norb),Btmp(Norb,Norb)
        complex(dp) :: tmp
        integer :: iw, jq, iorb, korb

        if (node.eq.0) then
            write(6,*) "ed_delta_new: Calculating delta_new..."
        endif

        ! @TODO H(k) is assumed to be diagonal
        do iw = 1,nwloc
            Gksum(:,iw) = cmplx(0.0_dp,0.0_dp)
            do jq = 1, Nq
                ! do iorb = 1, Norb
                !     do korb = 1, Norb
                !         Atmp(iorb,korb)=(1.0_dp/Gr(iorb,iw)+D_ev(iorb,iw) &
                !             +rMu)*kdel(iorb,korb)-Hk(iorb,korb,jq)
                !     enddo
                ! enddo
                ! call cinv(Atmp,Norb,Norb,Btmp)
                do iorb = 1, Norb
                    tmp = (1.0_dp/Gr(iorb,iw)+D_ev(iorb,iw)+rMu)-Hk(iorb,iorb,jq)
                    Gksum(iorb,iw) = Gksum(iorb,iw)+ 1.0_dp/tmp
                enddo
            enddo
            Gksum(:,iw) = Gksum(:,iw)/float(Nq)
            D_ev(:,iw) = 1.0_dp/Gr(:,iw)-1.0_dp/Gksum(:,iw) + D_ev(:,iw)
        enddo
    end subroutine ed_delta_new

    real(dp) function kdel(iorb,korb)
        integer:: iorb, korb

        if(iorb.eq.korb) kdel = 1.0_dp
        if(iorb.ne.korb) kdel = 0.0_dp

        return
    end function kdel

    subroutine n_from_gksum
        integer i,j,k,iw,iorb,nw
        complex(dp), allocatable:: Gksum_tot(:,:)

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
                nocc(iorb) = 2.0_dp*sum(real(gksum_tot(iorb,:)))/beta+0.50_dp
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
        complex(dp) ::D_ev(Norb,nwloc)

        do i = 1,nwloc
            do korb = 1, Norb
                D_ev(korb,i) = cmplx(0.0_dp,0.0_dp)
                do k = 1,Nbath
                    D_ev(korb,i) = D_ev(korb,i) + &
                        vk(korb,k)*vk(korb,k)/(cmplx(0.0_dp,omega(i))-ek(Norb+k))
                enddo
                D_ev(korb,i) = D_ev(korb,i) + ef(korb)
            enddo
        enddo
    end subroutine calc_dev
end module ed_green

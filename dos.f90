subroutine dos(omega_real,Aff)
    use ed_config
    use ed_solver
    implicit none

    integer :: iev,iev_read,iorb,iorb_read,istep,nstep_calcp,nstep_calcn
    double precision :: dw,Z,ev_read,prob_read
    double precision :: omega_real(nw)
    double precision :: Aff(Norb,nw),factor
    double precision :: ap(0:Nstep),bp(0:Nstep),an(0:Nstep),bn(0:Nstep)
    logical even

    Aff(:,:) = 0.D0

    open(unit=119,file="apbpanbn.dat",form="unformatted")

    do iev = 1, nev_calc
        do iorb = 1, Norb
            read(119) nstep_calcp, nstep_calcn
            read(119) iev_read,iorb_read,ev_read,prob_read,even,&
                (ap(istep),istep=0,nstep_calcp),(bp(istep),istep=0,nstep_calcp), &
                (an(istep),istep=0,nstep_calcn),(bn(istep),istep=0,nstep_calcn)

            if (iev.ne.iev_read.or.iorb.ne.iorb_read) then
                print *, "Wrong lanczos coefficient : ",iev,iorb,iev_read,iorb_read
                stop
            elseif (abs(ev_read-eigval(iev)%val).gt.1d-10.or. &
                    abs(prob_read-eigval(iev)%prob).gt.1d-10) then
                print *, "Eigenvalue mismatch : ",iev,iorb,ev_read,prob_read, &
                         eigval(iev)%val, eigval(iev)%prob
                stop
            endif

            call dos_diag(iorb,aff,omega_real,nev_calc,eigval(iev)%val,eigval(iev)%prob,&
                          even,nstep_calcp,nstep_calcn,ap,bp,an,bn)
        enddo

        if(even)  goto 1000

        do iorb  = 1, Norb
            read(119) nstep_calcp, nstep_calcn
            read(119) iev_read,iorb_read,ev_read,prob_read,even,&
                (ap(istep),istep=0,nstep_calcp),(bp(istep),istep=0,nstep_calcp), &
                (an(istep),istep=0,nstep_calcn),(bn(istep),istep=0,nstep_calcn)
            if (iev.ne.iev_read.or.iorb.ne.iorb_read) then
                print *, "Wrong lanczos coefficient : ",iev,iorb,iev_read,iorb_read
                stop
            elseif (abs(ev_read-eigval(iev)%val).gt.1d-10.or. &
                    abs(prob_read-eigval(iev)%prob).gt.1d-10) then
                print *, "Eigenvalue mismatch : ",iev,iorb,ev_read,prob_read, &
                         eigval(iev)%val, eigval(iev)%prob
                stop
            endif

            call dos_diag(iorb,aff,omega_real,nev_calc,eigval(iev)%val,eigval(iev)%prob,&
                          even,nstep_calcp,nstep_calcn,ap,bp,an,bn)
        enddo
1000    continue
    enddo
    close(119)

    return
end subroutine dos

subroutine dos_diag(iorb,aff,omega_real,nev_calc,ev,prob,even,nstep_calcp,nstep_calcn,&
                    ap,bp,an,bn)

    use ed_config
    implicit none

    integer :: iorb, nev_calc, nstep_calcp, nstep_calcn
    logical :: even
    double precision :: omega_real(nw),ev,prob
    double precision :: an(0:nstep_calcn),bn(0:nstep_calcn),ap(0:nstep_calcp),bp(0:nstep_calcp)
    double precision:: Aff(Norb,Nw)
    
    integer :: iw, istep
    double precision :: factor
    double complex:: gr_tmp, cpx_omega, grx

    if (even) then
        factor = 1.0D0
    else
        factor = 0.5D0
    endif

    do iw = 1, Nw
        cpx_omega = cmplx(ev+omega_real(iw),small)
        gr_tmp = bp(nstep_calcp)/(cpx_omega-ap(nstep_calcp))
        do istep = nstep_calcp-1, 0, -1
            gr_tmp = bp(istep)/(cpx_omega-ap(istep)-gr_tmp)
        enddo

        grx = gr_tmp

        cpx_omega = cmplx(omega_real(iw)-ev,small)
        gr_tmp = bn(nstep_calcn)/(cpx_omega+an(nstep_calcn))
        do istep = nstep_calcn-1, 0, -1
            gr_tmp = bn(istep)/(cpx_omega+an(istep)-gr_tmp)
        enddo
        grx = grx + gr_tmp
        Aff(iorb,iw) = -aimag(grx)/pi
    enddo

    return
end subroutine dos_diag

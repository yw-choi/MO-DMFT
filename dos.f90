        subroutine dos(Nstep,Norb,small,Nw,omega_real,nev_calc,&
                       Aff)
        implicit none
    

        double precision small,dw,Z
        integer:: Norb,Nstep,nev_calc,nw,ishift,io,itmp,k,i
        double precision::omega_real(nw),pev(nev_calc),E0
        double precision::Aff(Norb,nw),Aff_tmp(Norb,nw),factor
        double precision:: ap(0:Nstep),bp(0:Nstep),an(0:Nstep),bn(0:Nstep)
        logical even

        write(6,*) 
        write(6,*) "************* SPECTRAL FUNCTION CALCULATION  &
                   **************"
        write(6,*) 

        Aff(:,:) = 0.D0 
        Aff_tmp(:,:) = 0.D0 

        open(unit=11,file="apbpanbn.dat",form="unformatted")

        do i = 1, nev_calc
!          Diagonal component
           do io = 1, Norb
              read(11) itmp,itmp,E0,pev(i),even,&
                                  (ap(k),k=0,Nstep),(bp(k),k=0,Nstep), &
                                  (an(k),k=0,Nstep),(bn(k),k=0,Nstep)

              call dos_diag(Nstep,Norb,Nw,io,Aff_tmp,omega_real,E0,&
                            small,ap,bp,an,bn)
           enddo
           factor = 0.5D0
           if(even)  factor = 1.D0
           Aff(:,:) = Aff(:,:) + pev(i)*Aff_tmp(:,:)*factor
           if(even)  goto 1000

           do io  = 1, Norb
              read(11) itmp,itmp,E0,pev(i),even,&
                               (ap(k),k=0,Nstep),(bp(k),k=0,Nstep), &
                               (an(k),k=0,Nstep),(bn(k),k=0,Nstep)

              call dos_diag(Nstep,Norb,Nw,io,Aff_tmp,omega_real,E0,&
                         small,ap,bp,an,bn)
           enddo
           Aff(:,:) = Aff(:,:) + pev(i)*Aff_tmp(:,:)*factor
 1000      continue
        enddo
        close(11)

        return
        end

        SUBROUTINE dos_DIAG(Nstep,Norb,Nw,io,Aff,omega_real,E0,&
                            small,ap,bp,an,bn)

        implicit none

        integer i,j,k,io,Nstep,nw,Norb
        double precision:: omega_real(Nw),small
        double precision:: an(0:Nstep),bn(0:Nstep),ap(0:Nstep),bp(0:Nstep)

        double complex:: ztmp, cpx_omega,grx
        double precision:: Aff(Norb,Nw),pi,E0
        parameter(pi=acos(-1.0D0))

        do j = 1, Nw

           cpx_omega=dcmplx(E0+omega_real(j),small)
           grx = bp(Nstep)/(cpx_omega-ap(Nstep))
           do k = Nstep-1, 0, -1
              grx = bp(k)/(cpx_omega-ap(k)-grx)
           enddo 

           ztmp = grx

           cpx_omega = dcmplx(omega_real(j)-E0,small)
           grx = bn(Nstep)/(cpx_omega+an(Nstep))
           do k = Nstep-1, 0, -1
              grx = bn(k)/(cpx_omega+an(k)-grx)
           enddo
           grx = ztmp+grx
           Aff(io,j) = -aimag(grx)/pi
       enddo

       return
       end


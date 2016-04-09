      FUNCTION FUNC(x,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,comm)
       
      implicit none

      include 'mpif.h'
      integer comm,nprocs,ierr,taskid,master
      parameter(master=0)
      
      integer::iw,ksite,iorb,korb,Ns,nwloc,nxsize,Nbath,Norb,nw,i
      double complex::cpx_omega,delta(Norb,nwloc),cpx_temp,&
                      D_ev(Norb,nwloc)
      double precision::x(nxsize),func,func_loc,ek(Ns),vk(Norb,Norb+1:Ns),&
                      omega(nwloc),ef(Norb),kdel

      call mpi_comm_size(comm,nprocs,ierr)
      call mpi_comm_rank(comm,taskid,ierr)

      call x_to_ev(x,nxsize,Ns,Nbath,Norb,ef,ek,vk)

      do iw = 1, nwloc
         do korb = 1, Norb
            delta(korb,iw) = dcmplx(0.D0,0.D0)
            do ksite = Norb+1, Ns
               delta(korb,iw) = delta(korb,iw) + &
                                vk(korb,ksite)*vk(korb,ksite) &
                                /(dcmplx(0.D0,omega(iw))-ek(ksite))
            enddo
            delta(korb,iw) = delta(korb,iw) + ef(korb)
         enddo
      enddo

      func_loc = 0.D0
      do iw = 1, nwloc
         do korb = 1, Norb
            cpx_temp = delta(korb,iw) - D_ev(korb,iw)
            func_loc = func_loc + abs(cpx_temp)*abs(cpx_temp)/omega(iw)
         enddo
      enddo
      
      call mpi_allreduce(nwloc,nw,1,mpi_integer,mpi_sum,comm,ierr)
      call mpi_allreduce(func_loc,func,1,mpi_double_precision,mpi_sum,comm,ierr)

      func = func/float(Norb*Norb*Nw)

      return
      end

      subroutine dfunc(x,p,nwloc,Ns,nxsize,omega,D_ev,Nbath,Norb,comm)

      implicit none

      include 'mpif.h'
      integer comm,nprocs,taskid,ierr,master
      parameter(master=0)

      integer::i,j,k,iorb,korb,Ns,nwloc,nxsize,Nbath,Norb,iw,Nek,&
               ibath,ksite,iadd,nw
      double complex:: cpx_omega,delta(Norb,nwloc), cpx_temp2,    &
                      D_ev(Norb,nwloc)
      double precision::x(nxsize),p(nxsize),p_loc(nxsize),ek(Ns),&
                        vk(Norb,Norb+1:Ns),omega(nwloc),ef(Norb),&
                        fp1,fp2,func,ptmp(nxsize),kdel
   
!     do i = 1, nxsize
!        x(i) = x(i) + 0.00001
!        fp1 = func(x,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,comm)
!        x(i) = x(i) - 0.00002
!        fp2 = func(x,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,comm)
!        ptmp(i)= (fp1-fp2)/0.00002 
!        x(i) = x(i) + 0.00001
!     enddo

      call x_to_ev(x,nxsize,Ns,Nbath,Norb,ef,ek,vk)

      do iw = 1,nwloc
         do korb = 1, Norb
            delta(korb,iw) = dcmplx(0.D0,0.D0)
            do ksite = Norb+1, Ns
               delta(korb,iw) = delta(korb,iw) + &
                                 vk(korb,ksite)*vk(korb,ksite) &
                                   /(dcmplx(0.D0,omega(iw))-ek(ksite))
            enddo
            delta(korb,iw) = delta(korb,iw) + ef(korb)
         enddo
      enddo

      iadd = 0
      Nek = Nbath

      do ibath = 1, Nek
         p_loc(ibath) = 0.D0
         do korb = 1, Norb
            do iw = 1, nwloc
                p_loc(ibath) = p_loc(ibath) + 2.D0*real(             &
                    dconjg(delta(korb,iw)-D_ev(korb,iw))             &
                   *vk(korb,Norb+ibath)*vk(korb,Norb+ibath)          &
                   /(dcmplx(0.D0,omega(iw))-ek(Norb+ibath))          &
                   /(dcmplx(0.D0,omega(iw))-ek(Norb+ibath))          &
                   )/omega(iw)
            enddo
         enddo
      enddo

      do iorb = 1, Norb
         do ksite = Norb+1, Ns
            i = Nek+ksite-Norb+(iorb-1)*Nbath
            p_loc(i) = 0.D0
            do iw = 1, nwloc
               p_loc(i) = p_loc(i)+2.D0*real(                       &
                          dconjg(delta(iorb,iw)-D_ev(iorb,iw))      & 
                         *vk(iorb,ksite)/(dcmplx(0.D0,omega(iw))-ek(ksite))   &
                         + dconjg(delta(iorb,iw)-D_ev(iorb,iw)) &
                         *vk(iorb,ksite)/(dcmplx(0.D0,omega(iw))-ek(ksite))   &
                         )/omega(iw)
            enddo
          enddo
      enddo

      do iorb = 1, Norb
         i = nxsize-Norb+iorb
         p_loc(i) = 0.D0
         do iw = 1,nwloc
            p_loc(i) = p_loc(i) + 2.D0*real(delta(iorb,iw)-D_ev(iorb,iw))/omega(iw)
         enddo
      enddo

      do i = 1, nxsize
         call mpi_allreduce(p_loc(i),p(i),1,mpi_double_precision,mpi_sum,comm,ierr)
      enddo
      call mpi_allreduce(nwloc,nw,1,mpi_integer,mpi_sum,comm,ierr)

      p(:) = p(:)/float(norb*norb)/float(nw)

!     if(taskid.eq.master) then
!        do i = 1, nxsize
!           write(6,*) i, p(i), ptmp(i)
!        enddo
!     endif
!     if(taskid.eq.master) write(6,*)  
!     stop

      return
      end
 
      subroutine minimization(x,nwloc,Ns,nxsize,omega,D_ev,Nbath,Norb,comm,dab)
      
      implicit none

      include 'mpif.h' 
      integer comm,nprocs,taskid,ierr

      double precision, external:: func
      integer:: iter, iorb,k,Ns,nwloc,nxsize,Nbath,Norb
      double precision:: dab, x(nxsize),FTOL,fret,omega(nwloc)
      double complex:: D_ev(Norb,nwloc)
      PARAMETER(FTOL=1.0D-7)

      call mpi_comm_size(comm,nprocs,ierr)
      call mpi_comm_rank(comm,taskid,ierr)

      call frprmn(x,nxsize,FTOL,iter,nxsize,nwloc,Ns,omega,D_ev,Nbath,&
                  Norb,dab,comm)
      return
      end
 
      subroutine ev_to_x(ek,vk,ef,Ns,Nbath,nxsize,x)

      implicit none
      integer:: i,k,nxsize,Ns,Nbath,iorb,iinit,ifina,Norb
      double precision x(nxsize)
      double precision ek(Ns),vk(Ns-Nbath,Ns-Nbath+1:Ns),ef(Ns-Nbath)

      Norb = Ns - Nbath

      do i = 1, Nbath
         x(i) = ek(Norb+i)
      enddo

      do iorb = 1, Norb
         iinit = Nbath+1+(iorb-1)*Nbath
         ifina = Nbath+iorb*Nbath
         do k = iinit, ifina
            x(k) = vk(iorb,Norb+k-iinit+1)
         enddo
      enddo

      do i = 1, Norb
         x(nxsize-Norb+i) = ef(i)
      enddo

      return
      end

      subroutine x_to_ev(x,nxsize,Ns,Nbath,Norb,ef,ek,vk)

      implicit none
      integer:: i,nxsize,Ns,Nbath,Norb,iorb,ia,ib
      double precision x(nxsize),ek(Ns),vk(Norb,Norb+1:Ns),ef(Norb)

      do i =1, Norb
         ek(i) = x(nxsize-Norb+i)
         ef(i) = ek(i)
      enddo

      do i = Norb+1, Norb+Nbath
           ek(i) = x(i-Norb)
      enddo

      do iorb = 1, Norb
         ia = (iorb-1)*Nbath
         do i = Norb+1, Ns
            vk(iorb,i) = x(i+Nbath+ia-Norb)
         enddo
      enddo

      return
      end


      subroutine proj(u,v,x,n)

      implicit none 

      integer n      
      double precision u(n),v(n),x(n)

      x(:) = dot_product(u,v)/dot_product(u,u)*u(:)
      
      return
      end

      subroutine Grahm_Schmidt(nvec,ndim,vvec,uvec)

      implicit none
      integer nvec,ndim,i,j,k
      double precision uvec(ndim,nvec),vvec(ndim,nvec),x(ndim),y(ndim)

      uvec(:,1) = vvec(:,1) 
      do i = 2, nvec
         y(:) = 0.D0
         do j = 1, i-1
            call proj(uvec(:,j),vvec(:,i),x,ndim)
            y(:) = y(:) + x(:)
         enddo      
         uvec(:,i) = vvec(:,i) - y(:) 
         uvec(:,i) = uvec(:,i)/dsqrt(dot_product(uvec(:,i),uvec(:,i)))
      enddo


!/////////////////////////////////
!
!    Check orthonormalization
!
!/////////////////////////////////

!    do i = 1, nvec
!       do j = 1, nvec
!          write(6,*) i,j, dot_product(uvec(:,i),uvec(:,j))
!       enddo
!    enddo; stop
 
      return
      end


      subroutine msum(N,A,dab)

      implicit none

      integer N,i,j
      double complex A(N,N), dab

      do i = 1, N
         do j = 1, N
            dab = dab + A(i,j)
         enddo
      enddo

      return
      end
 
      subroutine mpi_dot_product(A,B,n,comm,dab)

      implicit none
      include 'mpif.h'
  
      integer n,taskid,ierr,comm,nprocs,MASTER
      parameter(MASTER=0)
      double precision A(n), B(n), dab, dab_tmp

      call mpi_comm_rank(comm,taskid,ierr)
      call mpi_comm_size(comm,nprocs,ierr) 

      dab_tmp = dot_product(A(:),B(:)) 
      call mpi_allreduce(dab_tmp,dab,1,mpi_double_precision,&
                         mpi_sum,comm,ierr)

      return    
      end


       subroutine spin_up_down_exchange(ns,nloc,comm,displ,nlocal,ne,&
                                        nu,ncpu,p0,p0p)

       use ed_utils
       implicit none  
       include 'mpif.h'
 
        integer ierr,taskid,nprocs,comm,req,master,nloc,nglob
        parameter(master=0) 
        integer i,k,ncpu,nd
        integer ns,nev,ipos,nu,ne,mypos,displ(0:ncpu-1),&
               nlocal(0:ncpu-1)
        integer imok, inam, nsCnu, nsCnd
        double precision p0(nloc),p0p(nloc)

        double precision, allocatable:: all_tmp(:)

       call mpi_comm_size(comm,nprocs,ierr)
       call mpi_comm_rank(comm,taskid,ierr)
      
       call mpi_allreduce(nloc,nglob,1,mpi_integer,mpi_sum,comm,ierr)

       nd = ne - nu
       nsCnu = icom(Ns,nu)
       nsCnd = icom(Ns,nd)

       allocate(all_tmp(nglob))
!      allocate(p0p(nloc))

       call mpi_allgatherv(p0,nloc,mpi_double_precision,all_tmp,&
                           nlocal,displ,mpi_double_precision,comm,ierr)

       do i = 1, nloc
          mypos = i + displ(taskid)
          imok = (mypos-1)/nsCnu
          ipos = imok + 1 + nsCnd * mod(mypos-1,nsCnu)
          p0p(i) = all_tmp(ipos)
       enddo
!      p0(:) = p0p(:)

       deallocate(all_tmp)

       return
       end

       subroutine boltzmann_factor(elv,nev,beta,p)
    
       implicit none

       integer:: i, k, nev
       double precision:: elv(nev), beta, p(nev), Z, en,betatmp

!      write(6,*) "This calculation assumes zero temperature!!!!"
!      betatmp = beta 
!      beta = 10000*beta
       do i = 1, nev
          Z = 0.D0
          en = elv(i)
          do k = 1, nev
             if(-beta*(elv(k)-en).gt.10.D0) then
               p(i) = 0.D0
               goto 100
             else
               Z = Z + exp(-beta*(elv(k)-en)) 
             endif
          enddo
          p(i) = 1.D0/Z
  100     continue 
       enddo
!      beta = betatmp
       return
       end
    
      SUBROUTINE SORT(N,ARRIN,INDX)

!     SORTS AN ARRAY BY THE HEAPSORT METHOD
!     W. H. PREUSS ET AL. NUMERICAL RECIPES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     this routine is not useful for samll size of n.
!     in particular n=1 lead to operation error.(arrin(0) is not defined.)
!     if n is greater than ~ 20 this routine is more fast.
!
      DIMENSION ARRIN(N),INDX(N)
!     DIMENSION ARRIN(1),INDX(1)
      IF (N .EQ. 1) THEN
        INDX(1) = 1
        RETURN
      ENDIF
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
 
      SUBROUTINE write_wffile_name(isector,ieigen,taskid,wffile) 

      implicit none
      integer isector, ieigen,taskid
      character*20 wffile

      if(taskid.lt.10) then
        if(isector.lt.10.and.ieigen.lt.10) then
           write(wffile,'(a,i1,a,i1,a,i1)') "wfsector",isector,"level"&
                                             ,ieigen,"id",taskid
        elseif(isector.ge.10.and.ieigen.lt.10) then
           write(wffile,'(a,i2,a,i1,a,i1)') "wfsector",isector,"level"&
                                             ,ieigen,"id",taskid
        elseif(isector.lt.10.and.ieigen.ge.10) then
           write(wffile,'(a,i1,a,i2,a,i1)') "wfsector",isector,"level"&
                                             ,ieigen,"id",taskid
        elseif(isector.ge.10.and.ieigen.ge.10) then
           write(wffile,'(a,i2,a,i2,a,i1)') "wfsector",isector,"level"&
                                             ,ieigen,"id",taskid
        endif
      else
        if(isector.lt.10.and.ieigen.lt.10) then
           write(wffile,'(a,i1,a,i1,a,i2)') "wfsector",isector,"level"&
                                             ,ieigen,"id",taskid
        elseif(isector.ge.10.and.ieigen.lt.10) then
           write(wffile,'(a,i2,a,i1,a,i2)') "wfsector",isector,"level"&
                                             ,ieigen,"id",taskid
        elseif(isector.lt.10.and.ieigen.ge.10) then
           write(wffile,'(a,i1,a,i2,a,i2)') "wfsector",isector,"level"&
                                             ,ieigen,"id",taskid
        elseif(isector.ge.10.and.ieigen.ge.10) then
           write(wffile,'(a,i2,a,i2,a,i2)') "wfsector",isector,"level"&
                                             ,ieigen,"id",taskid
        endif
      endif

      return
      end
   
      subroutine n_from_gksum(norb,nwloc,Gksum,displ_w,nloc_w,ncpu,beta,comm)

      implicit none
      include "mpif.h"

      integer i,j,k,iw,iorb,nw,master
      parameter(master=0)
      integer norb, nwloc, ncpu, comm, taskid, nprocs, ierr
      integer displ_w(0:ncpu-1), nloc_w(0:ncpu-1)
      double precision nocc(norb),beta
      double complex Gksum(norb,nwloc)
      double complex, allocatable:: Gksum_tot(:,:)
     
      call mpi_comm_size(comm,nprocs,ierr)
      call mpi_comm_rank(comm,taskid,ierr)

      call mpi_allreduce(nwloc,nw,1,mpi_integer,mpi_sum,comm,ierr)

      if(taskid.eq.master) write(6,'(a,2x,i8)') "nw =", nw

      allocate(Gksum_tot(norb,nw))
      do iorb = 1, norb
         call mpi_gatherv(gksum(iorb,:),nwloc,mpi_double_complex,&
                       Gksum_tot(iorb,:),&
                       nloc_w,displ_w,mpi_double_complex,&
                       master,comm,ierr)
      enddo
      
      if(taskid.eq.master) then
        write(6,*) "Occupations from local Green function"
!       do i = 1, nw
!          write(100,'(4e)') (gksum_tot(iorb,i),iorb=1,2)
!       enddo
        do iorb = 1, norb
           nocc(iorb) = 2.D0*sum(real(gksum_tot(iorb,:)))/beta+0.5D0
           write(6,'(a,i2,3x,f10.7)') "Orbtal",iorb,nocc(iorb)
        enddo
        write(6,'(a,5x,f10.7)') "sum =", sum(nocc(:))
      endif
      
      deallocate(gksum_tot)
      return
      end

      subroutine find_gtau(Gr,norb,nw,displ_w,nlocal_w,beta,ncpu,comm)     

      implicit none

      include 'mpif.h'

      integer comm, ierr,norb,nw,ncpu,nwtot,i,j,k,nprocs,taskid
      integer displ_w(0:ncpu-1),nlocal_w(0:ncpu-1) 
      double precision beta, pi, dtau
      double complex Gr(norb,nw)
      double complex, allocatable:: Gr_tot(:),cin_fft(:),cout_fft(:),chi_w(:,:)
      double precision, allocatable:: Gr_tau(:,:),in_fft(:)
 
      pi = acos(-1.D0)
      call mpi_comm_size(comm,nprocs,ierr)
      call mpi_comm_rank(comm,taskid,ierr)

      call mpi_allreduce(nw,nwtot,1,mpi_integer,mpi_sum,comm,ierr) 
      dtau = beta/nwtot/2.D0

      if(taskid.eq.0) write(6,*) "nwtot=",nwtot
       
      allocate(Gr_tot(nwtot),Gr_tau(norb,4*nwtot),chi_w(norb,nwtot))
      do i = 1, norb  
         call mpi_allgatherv(Gr(i,:),nw,mpi_double_complex,Gr_tot(:),&
                             nlocal_w,displ_w,mpi_double_complex,comm,ierr)
!        if(taskid.eq.0) then
!           do j = 1, nwtot
!              write(700,'(2e)') Gr_tot(j)
!           enddo
!        endif 
   
         allocate(cin_fft(nwtot*4),cout_fft(nwtot*4))
         cin_fft(:) = dcmplx(0.D0,0.D0)
    
         do j = 1, nwtot
             cin_fft(2*j) = Gr_tot(j) - 1.D0/dcmplx(0.D0,(2*j-1)*pi/beta)
         enddo 
         do j = 2*nwtot+2,4*nwtot
             cin_fft(j) = dconjg(cin_fft(4*nwtot-j+2))
         enddo
         call cfft_1d_bk(4*nwtot,cin_fft,cout_fft)  
         cout_fft(:) = cout_fft(:)/beta
         Gr_tau(i,:) = real(cout_fft(:))+0.5D0
!        if(taskid.eq.0) then
!           do j = 1, 4*nwtot
!              write(900,'(1e)') Gr_tau(i,j)
!           enddo
!        endif
         deallocate(cin_fft,cout_fft) 
         allocate(in_fft(2*nwtot),cout_fft(nwtot+1))
         in_fft(1) = -Gr_tau(i,1)*Gr_tau(i,2*nwtot+1)
         do j = 2, 2*nwtot
             in_fft(j) = Gr_tau(i,j)*Gr_tau(i,2*nwtot-j+2)
         enddo
         call fft_1d(2*nwtot,in_fft,cout_fft) 
         do j = 1, nwtot
            chi_w(i,j) = cout_fft(j)*dtau
         enddo
         if(taskid.eq.0) then
            do j = 1, nwtot
               write(900,'(2e)') chi_w(i,j)
            enddo
         endif
         deallocate(in_fft,cout_fft)
      enddo
      deallocate(Gr_tot,chi_w)
      return
      end


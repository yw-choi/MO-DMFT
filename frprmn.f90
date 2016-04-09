       subroutine frprmn(p,n,ftol,iter,nxsize,nwloc,Ns,omega,D_ev,Nbath,&
                         Norb,fret,comm)

       implicit none

       include 'mpif.h'
       integer comm,nprocs,taskid,ierr,master
       parameter(master=0)

       integer:: iter, n, nmax, itmax,i
       double precision:: fret, ftol, p(n), eps, func
       parameter(nmax=60, itmax=20000, eps=1.d-10)
      
!      USES dfunc, func, linmin
!      Given a starting point p that is a vector of length n, Fletcher-Reeves-Polak-Ribiere minimization
!      is performed on a function func, using its gradient as calculated by a routine dfunc.
!      The convergence tolerance on the function value is input as ftol. Returned quantities are p 
!      (the location of the minimum), iter ( the # of iterations that were performed), and fret (the minimum value of function)
!      The routine linmin is called to perform line minimizations
!      parameter: nmax is the maximum anticipated value of n; itmax is the maximum allowd number of iterations;
!      eps is a small number of rectify special case of converging to exactly zero function value.

       integer its, j 
       double precision dgg, fp, gam, gg, g(nmax), h(nmax), xi(nmax),zero
       integer nwloc,Ns,Nbath,Norb,nxsize
       double precision omega(nwloc)
       double complex D_ev(Norb,nwloc)      
 
       call mpi_comm_size(comm,nprocs,ierr)
       call mpi_comm_rank(comm,taskid,ierr)

       zero=0.D0

       fp = func(p,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,comm)     ! initialization
       if(taskid.eq.master)  write(6,'("Initial difference :",e)') fp
       call dfunc(p,xi,nwloc,Ns,nxsize,omega,D_ev,Nbath,Norb,comm)

       if(taskid.eq.0) then
          write(6,*)   "               x                    df/dx      "
          do i = 1, n
             write(6,*) p(i),xi(i)
          enddo
          write(6,*) 
       endif

       do j = 1, n
          g(j) = -xi(j)
          h(j) = g(j)
          xi(j) = h(j)
       enddo
       
       do its = 1, itmax
          iter = its
          call linmin(p,xi,n,fret,Ns,nwloc,nxsize,omega,&
                      D_ev,Nbath,Norb,comm)
          if(2.D0*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+eps)) then
            if(taskid.eq.master) &
               write(6,*) "After",its,"iteration converged." 
             return
          endif
          fp=fret
          call dfunc(p,xi,nwloc,Ns,nxsize,omega,D_ev,Nbath,Norb,comm)
          gg = 0.D0
          dgg = 0.D0
      
          do j = 1, n 
             gg = gg + g(j)*g(j)
!            dgg = dgg + xi(j)*x(j)
             dgg = dgg+(xi(j)+g(j))*xi(j)
          enddo
  
          if(gg.eq.zero) return

          gam = dgg/gg
          do j =1, n
             g(j) = -xi(j)
             h(j) = g(j)+gam*h(j)
             xi(j)=h(j)
          enddo
      enddo

      stop "frprmn maximum iterations exceeded!"
      return
      end

      SUBROUTINE linmin(p,xi,n,fret,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,comm)

      implicit none
      INTEGER n,NMAX  
      double precision fret,p(n),xi(n),TOL
      PARAMETER (NMAX=60,TOL=1.e-8)  
!U    USES brent,f1dim,mnbrak  
      INTEGER j,ncom  
      double precision ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brent  
!     COMMON /f1com/ pcom,xicom,ncom  
      EXTERNAL f1dim  

      integer Ns,nwloc,nxsize,Nbath,Norb,comm
      double complex D_ev(Norb,nwloc)
      double precision omega(nwloc) !, ef(Norb)

      do j=1,n  
        pcom(j)=p(j)  
        xicom(j)=xi(j)  
      enddo
      ax=0.  
      xx=1.  
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim,Ns,nwloc,nxsize,omega,D_ev,&
                  Nbath,Norb,pcom,xicom,comm)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin,Ns,nwloc,nxsize,omega,D_ev,&
                  Nbath,Norb,pcom,xicom,comm)
      do j=1,n  
        xi(j)=xmin*xi(j)  
        p(j)=p(j)+xi(j)  
      enddo
      return  
      END  

      FUNCTION brent(ax,bx,cx,f,tol,xmin,Ns,nwloc,nxsize,omega,D_ev,&
                     Nbath,Norb,pcom,xicom,comm)

      implicit none
      INTEGER ITMAX  
      double precision brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS  
      EXTERNAL f  
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)  
      INTEGER iter  
      double precision a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm  
      double precision pcom, xicom ! cklee

      integer Ns,nwloc,nxsize,Nbath,Norb,comm
      double precision omega(nwloc)!,ef(Norb)
      double complex D_ev(Norb,nwloc)

      a=min(ax,cx)  
      b=max(ax,cx)  
      v=bx  
      w=v  
      x=v  
      e=0.  
      fx=f(x,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)
      fv=fx  
      fw=fx  
      do 11 iter=1,ITMAX  
        xm=0.5*(a+b)  
        tol1=tol*abs(x)+ZEPS  
        tol2=2.*tol1  
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3  
        if(abs(e).gt.tol1) then  
          r=(x-w)*(fx-fv)  
          q=(x-v)*(fx-fw)  
          p=(x-v)*q-(x-w)*r  
          q=2.*(q-r)  
          if(q.gt.0.) p=-p  
          q=abs(q)  
          etemp=e  
          e=d  
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))  goto 1    
          d=p/q  
          u=x+d  
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)  
          goto 2  
        endif  
1       if(x.ge.xm) then  
          e=a-x  
        else  
          e=b-x  
        endif  
        d=CGOLD*e  
2       if(abs(d).ge.tol1) then  
          u=x+d  
        else  
          u=x+sign(tol1,d)  
        endif  
        fu=f(u,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)
        if(fu.le.fx) then  
          if(u.ge.x) then  
            a=x  
          else  
            b=x  
          endif  
          v=w  
          fv=fw  
          w=x  
          fw=fx  
          x=u  
          fx=fu  
        else  
          if(u.lt.x) then  
            a=u  
          else  
            b=u  
          endif  
          if(fu.le.fw .or. w.eq.x) then  
            v=w  
            fv=fw  
            w=u  
            fw=fu  
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then  
            v=u  
            fv=fu  
          endif  
        endif  
11    continue  
      stop 'brent exceed maximum iterations'  

3     xmin=x  
      brent=fx  
      return  
      END

      FUNCTION f1dim(x,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)

      implicit none
      INTEGER NMAX

      integer Ns,nwloc,nxsize,Nbath,Norb,comm  
      double precision omega(nwloc)!,ef(Norb)
      double complex D_ev(Norb,nwloc)

      double precision f1dim,func,x  
!     PARAMETER (NMAX=90)  
!U    USES func  
      INTEGER j,ncom  
      double precision pcom(nxsize),xicom(nxsize),xt(nxsize)
!     COMMON /f1com/ pcom,xicom,ncom  
      
      do j=1,nxsize
        xt(j)=pcom(j)+x*xicom(j)  
      enddo

      f1dim=func(xt,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,comm)
      return  
      END  

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func,Ns,nwloc,nxsize,&
                        omega,D_ev,Nbath,Norb,pcom,xicom,comm)

      implicit none

      integer Ns,nwloc,nxsize,Nbath,Norb,comm  
      double precision omega(nwloc)
      double complex D_ev(Norb,nwloc)

      double precision ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY  
      EXTERNAL func  
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)  
      double precision dum,fu,q,r,u,ulim  
      double precision pcom, xicom
      
      fa=func(ax,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)
      fb=func(bx,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)  

      if(fb.gt.fa)then  
        dum=ax  
        ax=bx  
        bx=dum  
        dum=fb  
        fb=fa  
        fa=dum  
      endif  
      cx=bx+GOLD*(bx-ax)  
      fc=func(cx,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)  
1     if(fb.ge.fc)then  
        r=(bx-ax)*(fb-fc)  
        q=(bx-cx)*(fb-fa)  
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))  
        ulim=bx+GLIMIT*(cx-bx)  
        if((bx-u)*(u-cx).gt.0.)then  
          fu=func(u,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)  
          if(fu.lt.fc)then  
            ax=bx  
            fa=fb  
            bx=u  
            fb=fu  
            return  
          else if(fu.gt.fb)then  
            cx=u  
            fc=fu  
            return  
          endif  
          u=cx+GOLD*(cx-bx)  
          fu=func(u,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)  
        else if((cx-u)*(u-ulim).gt.0.)then  
          fu=func(u,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)  
          if(fu.lt.fc)then  
            bx=cx  
            cx=u  
            u=cx+GOLD*(cx-bx)  
            fb=fc  
            fc=fu  
            fu=func(u,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)  
          endif  
        else if((u-ulim)*(ulim-cx).ge.0.)then  
          u=ulim  
          fu=func(u,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)  
        else  
          u=cx+GOLD*(cx-bx)  
          fu=func(u,Ns,nwloc,nxsize,omega,D_ev,Nbath,Norb,pcom,xicom,comm)  
        endif  
        ax=bx  
        bx=cx  
        cx=u  
        fa=fb  
        fb=fc  
        fc=fu  
        goto 1  
      endif  
      return  
      END  

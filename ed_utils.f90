module ed_utils
    use precision
    use parallel_params
    implicit none

contains

    real(dp) function mpi_dot_product(A,B,n) result(dab)
        include 'mpif.h'

        integer n
        real(dp) :: A(n), B(n), dab_tmp

        dab_tmp = sum(A(1:n)*B(1:n)) 
        call mpi_allreduce(dab_tmp,dab,1,mpi_double_precision,&
            mpi_sum,comm,ierr)

        return    
    end function mpi_dot_product

    real(dp) function mpi_norm(A,n) result(dab)
        include 'mpif.h'
        integer n
        real(dp) :: A(n)
        dab = sqrt(mpi_dot_product(A,A,n))
        return    
    end function mpi_norm

    integer function ifact(n)

        integer:: n,i

        ifact = 1
        do i = 1, n
        ifact = ifact*i
        enddo

        return
    end function ifact

    integer function iP(n,m)

        integer:: n,m,i

        iP = 1

        do i = n-m+1, n
        iP = iP*i
        enddo 

        return
    end function iP


    integer function iCom(n,m) 
        integer:: n, m, k

        if((n.eq.m).or.(m.eq.0)) then
            iCom = 1
        else
            if(m.gt.n/2) then
                k = n-m
                iCom = iP(n,k)/ifact(k)
            elseif(m.le.n/2) then
                iCom = iP(n,m)/ifact(m)
            endif 
        endif

        return
    end function icom

! Compute inverse of complex matrix

      subroutine cinv( A, N, MAXN, Ainv )

      integer*4 N, MAXN
      double complex A(MAXN,MAXN), Ainv(MAXN,MAXN)

! Inputs
!   A       Matrix A to be inverted
!   N       Elements used in matrix A (N by N)
!  MAXN     Matrix dimenstions as A(MAXN,MAXN)
! Outputs
!  Ainv     Inverse of matrix A

      integer*4 MAXMAXN
      parameter( MAXMAXN = 200 )
      integer*4 i, j, k, index(MAXMAXN), jPivot, indexJ
      double precision scale(MAXMAXN), scaleMax, ratio, ratioMax
      double complex AA(MAXMAXN,MAXMAXN), B(MAXMAXN,MAXMAXN), coeff, sum

      if( MAXN .gt. MAXMAXN ) then
        write(*,*) 'ERROR in cinv: Matrix too large'
        stop
      endif

      !* Matrix B is initialized to the identity matrix
      do i=1,N
       do j=1,N
         AA(i,j) = A(i,j)  ! Copy matrix so as not to overwrite
         B(i,j) = 0.0
       enddo
       B(i,i) = 1.0
      enddo

      !* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
      do i=1,N
        index(i) = i     ! Initialize row index list
        scaleMax = 0.0
        do j=1,N
          if( abs(AA(i,j)) .gt. scaleMax ) then
            scaleMax = abs(AA(i,j))
          endif
        enddo
        scale(i) = scaleMax
      enddo

      !* Loop over rows k = 1, ..., (N-1)
      do k=1,(N-1)
        !* Select pivot row from max( |a(j,k)/s(j)| )
        ratiomax = 0.0
        jPivot = k
        do i=k,N
          ratio = abs(AA(index(i),k))/scale(index(i))
          if( ratio .gt. ratiomax ) then
            jPivot=i
            ratiomax = ratio
          endif
        enddo
        !* Perform pivoting using row index list
        indexJ = index(k)
        if( jPivot .ne. k ) then     ! Pivot
          indexJ = index(jPivot)
          index(jPivot) = index(k)   ! Swap index jPivot and k
          index(k) = indexJ
        endif
        !* Perform forward elimination
        do i=k+1,N
          coeff = AA(index(i),k)/AA(indexJ,k)
          do j=k+1,N
            AA(index(i),j) = AA(index(i),j) - coeff*AA(indexJ,j)
          enddo
          AA(index(i),k) = coeff
          do j=1,N
            B(index(i),j) = B(index(i),j) - AA(index(i),k)*B(indexJ,j)
          enddo
        enddo
      enddo

      !* Perform backsubstitution
      do k=1,N
        Ainv(N,k) = B(index(N),k)/AA(index(N),N)
        do i=N-1,1,-1
          sum = B(index(i),k)
          do j=i+1,N
            sum = sum - AA(index(i),j)*Ainv(j,k)
          enddo
          Ainv(i,k) = sum/AA(index(i),i)
        enddo
      enddo

      return
      end subroutine cinv
end module ed_utils 

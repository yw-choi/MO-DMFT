module ed_utils
    use precision
    use parallel_params
    implicit none

contains

    real(dp) function mpi_dot_product(A,B,n) result(dab)
        include 'mpif.h'

        integer n
        real(dp) :: A(n), B(n), dab_tmp

        dab_tmp = dot_product(A(1:n),B(1:n)) 
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

end module ed_utils 

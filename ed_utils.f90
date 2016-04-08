module ed_utils

    implicit none

contains

    
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

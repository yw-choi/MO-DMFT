module ed_green
    use precision
    use parallel_params
    use ed_io
    use ed_config

    implicit none

contains

    subroutine ed_calc_green_ftn(nev_calc)
        integer, intent(in) :: nev_calc
        
        real(dp), allocatable :: eigval(:), pev(:)

        integer :: i

        call import_eigval(nev_calc,eigval,pev)

        if (node.eq.0) then
            do i=1,nev_calc
                print *, i,eigval(i),pev(i)
            enddo
        endif

    end subroutine ed_calc_green_ftn

end module ed_green

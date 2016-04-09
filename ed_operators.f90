module ed_operators
    use ed_config

    implicit none

contains
    integer function get_bitidx(isite,ispin)
        integer :: isite, ispin
        get_bitidx = (ispin-1)*Nsite + isite-1
    end function get_bitidx

    subroutine onsite(basis_in,isite,ispin,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: isite,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        real(dp), intent(out) :: coeff

        integer :: bitidx,sgn
        bitidx = get_bitidx(isite,ispin)

        call number_op(basis_in,isite,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0_dp
            return
        endif
        
        if (isite.le.Norb) then
            ! impurity orbitals
            coeff = ek(isite) - rMu
        else
            ! bath orbitals
            coeff = ek(isite)
        endif
    end subroutine onsite

    subroutine hybridization1(basis_in,iorb,ibath,ispin,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,ibath,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        real(dp), intent(out) :: coeff

        integer :: sgn
        coeff = 1.0_dp
        call destruction_op(basis_in,Norb+ibath,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0_dp
            return
        endif
        coeff = coeff*sgn
        call creation_op(basis_out,iorb,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0_dp
            return
        endif
        coeff = coeff*sgn
        coeff = coeff*vk(iorb,ibath)
    end subroutine hybridization1

    subroutine hybridization2(basis_in,iorb,ibath,ispin,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,ibath,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        real(dp), intent(out) :: coeff

        integer :: sgn
        coeff = 1.0_dp
        call destruction_op(basis_in,iorb,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0_dp
            return
        endif
        coeff = coeff*sgn
        call creation_op(basis_out,Norb+ibath,ispin,basis_out,sgn)
        if (sgn.eq.0) then
            coeff = 0.0_dp
            return
        endif
        coeff = coeff*sgn
        coeff = coeff*vk(iorb,ibath)
    end subroutine hybridization2

    subroutine dens_dens(basis_in,iorb,jorb,ispin,jspin,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,jorb,ispin,jspin
        integer(kind=kind_basis), intent(out) :: basis_out
        real(dp), intent(out) :: coeff

        integer :: n1,n2

        call number_op(basis_in,iorb,ispin,basis_out,n1)
        call number_op(basis_in,jorb,jspin,basis_out,n2)

        if (n1*n2.eq.0) then
            coeff = 0
            return
        endif

        if (iorb.eq.jorb) then
            if (ispin.eq.jspin) then
                coeff = 0.0_dp
                basis_out = 0
                return
            endif

            ! intra orbital coupling
            coeff = U
        else if (ispin.eq.jspin) then
            ! inter orbital coupling, same spin
            coeff = U - 3*Jex ! Uprime - J
        else
            ! inter orbital coupling, different spin
            coeff = U - 2*Jex ! Uprime 
        endif

    end subroutine dens_dens

    subroutine spin_flip(basis_in,iorb,jorb,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,jorb
        integer(kind=kind_basis), intent(out) :: basis_out
        real(dp), intent(out) :: coeff

        integer :: sgn_tot,sgn

        sgn_tot=1
        coeff = 0.0_dp

        if (iorb.eq.jorb) then
            basis_out = 0
            return
        endif

        ! c_j,u
        call destruction_op(basis_in,jorb,1,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c^+_j,d
        call creation_op(basis_out,jorb,2,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c_i,d
        call destruction_op(basis_out,iorb,2,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c^+_i,u
        call creation_op(basis_out,iorb,1,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        coeff = -Jex*sgn_tot ! Jprime

        return
    end subroutine spin_flip

    subroutine pair_exchange(basis_in,iorb,jorb,basis_out,coeff)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: iorb,jorb
        integer(kind=kind_basis), intent(out) :: basis_out
        real(dp), intent(out) :: coeff

        integer :: sgn_tot,sgn

        sgn_tot=1
        coeff = 0.0_dp

        if (iorb.eq.jorb) then
            basis_out = 0
            return
        endif

        ! c_j,d
        call destruction_op(basis_in,jorb,2,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c_j,u
        call destruction_op(basis_out,jorb,1,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c^+_i,d
        call creation_op(basis_out,iorb,2,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        ! c^+_i,u
        call creation_op(basis_out,iorb,1,basis_out,sgn)
        if (sgn.eq.0) then
            return 
        else
            sgn_tot = sgn_tot*sgn
        endif

        coeff = -Jex*sgn_tot ! Jprime

        return
    end subroutine pair_exchange

    subroutine number_op(basis_in,isite,ispin,basis_out,sgn)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: isite,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        integer, intent(out) :: sgn

        integer :: bitidx

        bitidx = get_bitidx(isite,ispin)
        if (BTEST(basis_in,bitidx)) then
            sgn=1
            basis_out = basis_in
        else
            sgn=0
            basis_out = 0
        endif

    end subroutine number_op

    subroutine destruction_op(basis_in,isite,ispin,basis_out,sgn)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: isite,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        integer, intent(out) :: sgn

        ! internal variables
        integer :: i, bitidx, n

        bitidx = get_bitidx(isite,ispin)
        if (BTEST(basis_in,bitidx)) then
            n = 0
            do i=0,bitidx-1
                if (BTEST(basis_in,i)) then
                    n = n + 1
                endif
            enddo
            basis_out = IBCLR(basis_in,bitidx)
            if (mod(n,2).eq.1) then
                sgn = -1
            else
                sgn = +1
            endif
        else
            basis_out = 0
            sgn = 0
        end if
    end subroutine destruction_op

    subroutine creation_op(basis_in,isite,ispin,basis_out,sgn)
        integer(kind=kind_basis), intent(in) :: basis_in
        integer, intent(in) :: isite,ispin
        integer(kind=kind_basis), intent(out) :: basis_out
        integer, intent(out) :: sgn

        ! internal variables
        integer :: i, bitidx, n

        bitidx = get_bitidx(isite,ispin)
        if (.not.BTEST(basis_in,bitidx)) then
            n = 0
            do i=0,bitidx-1
                if (BTEST(basis_in,i)) then
                    n = n + 1
                endif
            enddo
            basis_out = IBSET(basis_in,bitidx)
            if (mod(n,2).eq.1) then
                sgn = -1
            else
                sgn = +1
            endif
        else
            basis_out = 0
            sgn = 0
        end if
    end subroutine creation_op
end module ed_operators

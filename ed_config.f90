module ed_config
    use precision
    use fdf

    implicit none
    integer, parameter :: kind_basis = 4

    ! dimensions of the problem
    integer :: Norb
    integer :: Nbath
    integer :: Nsite

    integer :: Nsector
    integer, allocatable :: Nelec(:)
    integer, allocatable :: Nup(:)

    ! physical parameters
    real(dp) :: U
    real(dp) :: Jex
    real(dp) :: rMu
    real(dp) :: beta

    real(dp), allocatable :: ek(:)
    real(dp), allocatable :: vk(:,:)
    
    ! calculation parameters
    real(dp) :: small
    real(dp) :: scf_tol
    integer :: Nstep
    integer :: Niter
    integer :: Nev

    integer :: Nq
    integer :: Nw

!-- Local Variables
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
contains
    
    subroutine ed_config_init
        integer :: i,j

        Norb = fdf_integer("DMFT.Norbs",2)
        Nbath = fdf_integer("DMFT.Nbaths",2)
        Nsite = Norb+Nbath

        Nsector = fdf_integer("DMFT.Nsectors",2)
        allocate(nelec(nsector),nup(nsector))
        if (fdf_block('DMFT.Sectors', bfdf)) then
            i = 1
            do while((fdf_bline(bfdf, pline)) .and. (i .le. nsector))
                nelec(i) = fdf_bintegers(pline,1)
                nup(i) = fdf_bintegers(pline,2)
                i = i + 1
            enddo
        endif

        U = fdf_double("DMFT.U", 1.0_dp)
        Jex = fdf_double("DMFT.J",0.3_dp) 
        rMu = fdf_double("DMFT.Mu",1.0_dp)
        beta = fdf_double("DMFT.beta", 100.0_dp)

        allocate(ek(nsite),vk(norb,nbath))
        if (fdf_block('DMFT.Baths', bfdf)) then
            i = 1
            do while((fdf_bline(bfdf, pline)) .and. (i .le. nbath))
                ek(Norb+i) = fdf_breals(pline,1)
                i = i + 1
            enddo
        endif

        if (fdf_block('DMFT.BathCouplings', bfdf)) then
            i = 1
            do while((fdf_bline(bfdf, pline)) .and. (i .le. nbath))
                do j=1,norb
                    vk(j,i) = fdf_breals(pline,j)
                enddo
                i = i + 1
            enddo
        endif

        small = fdf_double("DMFT.small", 0.05_dp)

        Nstep = fdf_integer("DMFT.ContinuedFractionSteps", 50)

        Niter = fdf_integer("DMFT.MaxIterations", 100)
        
        scf_tol = fdf_double("DMFT.ScfTolerance", 1d-4)

        nev = fdf_integer("DMFT.Nev",30)
        Nq = fdf_integer("DMFT.Nq",10000)
        Nw = fdf_integer("DMFT.Nw",2000)

    end subroutine ed_config_init



end module ed_config

module ed_config
    use precision
    use fdf
    use parallel_params

    implicit none
    real(dp), parameter :: pi = 4.0_dp*ATAN(1.0_dp)
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
    integer :: Nloop
    integer :: Nev

    integer :: Nq
    integer :: Nw

!-- Local Variables
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
contains

    subroutine ed_read_options
        integer :: i,j
        character(len=100) :: text
        if (node.eq.0) then
            write(6,*)
            write(6,'(a)') repeat("=",80)
            write(6,'(4x,a)') "Input Parameters"
            write(6,'(a)') repeat("=",80)
        endif

        Norb = fdf_integer("DMFT.Norbs",2)
        Nbath = fdf_integer("DMFT.Nbaths",2)
        Nsite = Norb+Nbath

        if (node.eq.0) then
            write(text,*) "Number of impurity orbitals"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Norb
            write(text,*) "Number of bath orbitals"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nbath
            write(text,*) "Number of sites"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nsite
        endif

        Nsector = fdf_integer("DMFT.Nsectors",2)
        if (node.eq.0) then
            write(text,*) "Number of sectors"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nsector
        endif
        allocate(nelec(1:nsector),nup(1:nsector))

        if (fdf_block('DMFT.Sectors', bfdf)) then
            i = 1
            do while((fdf_bline(bfdf, pline)) .and. (i .le. nsector))
                nelec(i) = fdf_bintegers(pline,1)
                nup(i) = fdf_bintegers(pline,2)
                i = i + 1
            enddo
        endif

        U = fdf_double("DMFT.U", 1.0D0)
        Jex = fdf_double("DMFT.J",0.3D0) 
        rMu = fdf_double("DMFT.Mu",1.0D0)
        beta = fdf_double("DMFT.beta", 100.0D0)
        if (node.eq.0) then
            write(text,*) "U"
            write(6,'(3x,a40,2x,a,2x,F8.3)') text, '=', U
            write(text,*) "J"
            write(6,'(3x,a40,2x,a,2x,F8.3)') text, '=', Jex
            write(text,*) "chemical potential"
            write(6,'(3x,a40,2x,a,2x,F8.3)') text, '=', rMu
            write(text,*) "beta (inverse temperature)"
            write(6,'(3x,a40,2x,a,2x,F8.3)') text, '=', beta
        endif

        allocate(ek(1:nsite),vk(1:norb,1:nbath))
        if (fdf_block('DMFT.Baths', bfdf)) then
            i = 1
            do while( (i .le. nbath) .and. (fdf_bline(bfdf, pline)))
                ek(Norb+i) = fdf_breals(pline,1)
                i = i + 1
            enddo
        endif

        if (fdf_block('DMFT.BathCouplings', bfdf)) then
            i = 1
            do while( i.le.nbath .and. (fdf_bline(bfdf, pline)))
                do j=1,norb
                    vk(j,i) = fdf_breals(pline,j)
                enddo
                i = i + 1
            enddo
        endif

        small = fdf_double("DMFT.small", 0.05D0)

        Nstep = fdf_integer("DMFT.ContinuedFractionSteps", 50)

        Nloop = fdf_integer("DMFT.MaxIterations", 100)
        
        scf_tol = fdf_double("DMFT.ScfTolerance", 1d-4)

        nev = fdf_integer("DMFT.Nev",30)
        Nq = fdf_integer("DMFT.Nq",10000)
        Nw = fdf_integer("DMFT.Nw",2000)

        if (node.eq.0) then
            write(text,*) "Broadening parameter"
            write(6,'(3x,a40,2x,a,2x,F8.3)') text, '=', small
            write(text,*) "Number of steps in continued fraction"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nstep
            write(text,*) "Maximum DMFT iterations"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', Nloop
            write(text,*) "DMFT SCF tolerance"
            write(6,'(3x,a40,2x,a,2x,ES8.1)') text, '=', scf_tol
            write(text,*) "Number of eigenvalues to be computed"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', nev
            write(text,*) "Number of k-points in band structure"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', nq
            write(text,*) "Number of frequency points"
            write(6,'(3x,a40,2x,a,2x,I8)') text, '=', nw
        endif

        if (node.eq.0) then
            write(6,'(a)') repeat("=",80)
            write(6,*)
        endif
    end subroutine ed_read_options

end module ed_config

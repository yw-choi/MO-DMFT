module ed_config
    use precision
    use fdf

    implicit none

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

    integer :: kind_basis = 4
contains
    
    subroutine ed_config_init

        Norb = fdf_integer("DMFT.Norbs",2)
        Nbath = fdf_integer("DMFT.Nbaths",2)
        Nsite = Norb+Nbath


    end subroutine ed_config_init



end module ed_config

module ed_io
    use precision
    use parallel_params
    use ed_config
    use sys

    implicit none

    integer, parameter, private :: EIGVEC_IUNIT_BASE=100
    integer, parameter, private :: EIGVAL_IUNIT_BASE=101
contains

    subroutine export_eigvec(isector,iev,node,nloc,nev,eigvec)
        integer, intent(in) :: isector, iev, node, nloc, nev
        real(dp), intent(in) :: eigvec(nloc)

        integer :: i,iunit
        character(len=150) :: fname

        call eigvec_filename(isector,iev,node,fname)

        iunit = EIGVEC_IUNIT_BASE + node
#ifdef DATA_PLAINTEXT
        open(unit=iunit,file=fname,status="replace")
        write(iunit,"(F10.3)") (eigvec(i),i=1,nloc)
        close(iunit)
#else
        open(unit=iunit,file=fname,status="replace",form="unformatted")
        write(iunit) (eigvec(i),i=1,nloc)
        close(iunit)
#endif
    end subroutine

    subroutine export_eigval(nev_calc,eigval_all,pev,ind)
        integer, intent(in) :: nev_calc, ind(nsector*nev)
        real(dp), intent(in) :: eigval_all(nev*nsector)
        real(dp), intent(in) :: pev(nev*nsector)

        integer :: k

#ifdef DATA_PLAINTEXT
        open(unit=EIGVAL_IUNIT_BASE,file="eigenvalues.dat",form="formatted",&
            status="replace")
        write(EIGVAL_IUNIT_BASE,"(I4)") nev_calc
        write(EIGVAL_IUNIT_BASE,"(I4,2F10.3)") (ind(k),eigval_all(ind(k)),pev(ind(k)),k=1,nev_calc)
#else
        open(unit=EIGVAL_IUNIT_BASE,file="eigenvalues.dat",form="unformatted",&
            status="replace")
        write(EIGVAL_IUNIT_BASE) nev_calc
        write(EIGVAL_IUNIT_BASE) (ind(k),eigval_all(ind(k)),pev(ind(k)),k=1,nev_calc)
#endif
        close(EIGVAL_IUNIT_BASE)

    end subroutine export_eigval

    subroutine import_eigval(nev_calc,eigval,pev)
        integer, intent(in) :: nev_calc
        real(dp), allocatable, intent(out) :: eigval(:), pev(:)

        integer :: i, iunit, idx, nev_calc_read
        iunit = EIGVAL_IUNIT_BASE
#ifdef DATA_PLAINTEXT
        open(unit=iunit,file="eigenvalues.dat",status="old",form="formatted")
#else
        open(unit=iunit,file="eigenvalues.dat",status="old",form="unformatted")
#endif

#ifdef DATA_PLAINTEXT
        read(iunit,*) nev_calc_read
#else
        read(iunit) nev_calc_read
#endif
        if (nev_calc_read.ne.nev_calc) then
            if (node.eq.0) then
                write(6,*) "nev_calc mismatch."
            endif
            call die
            return
        endif

        allocate(eigval(nev_calc),pev(nev_calc))
        do i=1,nev_calc
#ifdef DATA_PLAINTEXT
            read(iunit,*) idx, eigval(i), pev(i)
#else
            read(iunit) idx, eigval(i), pev(i)
#endif
        enddo
        close(iunit)

    end subroutine import_eigval

    subroutine eigvec_filename(isector,iev,inode,fname)
        integer, intent(in) :: isector, iev, inode
        character(len=150), intent(out) :: fname

        write(fname,"(A7,I2.2,A1,I3.3,A1,I4.4)") "eigvec-",isector,"-",iev,"-",inode

    end subroutine eigvec_filename

end module ed_io

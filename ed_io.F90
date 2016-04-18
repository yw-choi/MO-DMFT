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
        write(iunit,"(I)") nloc
        write(iunit,"(F24.16)") (eigvec(i),i=1,nloc)
#else
        open(unit=iunit,file=fname,status="replace",form="unformatted")
        write(iunit) nloc
        write(iunit) (eigvec(i),i=1,nloc)
#endif
        close(iunit)
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
        write(EIGVAL_IUNIT_BASE,"(I4,2F24.16)") (ind(k),eigval_all(ind(k)),pev(ind(k)),k=1,nev_calc)
#else
        open(unit=EIGVAL_IUNIT_BASE,file="eigenvalues.dat",form="unformatted",&
            status="replace")
        write(EIGVAL_IUNIT_BASE) nev_calc
        write(EIGVAL_IUNIT_BASE) (ind(k),eigval_all(ind(k)),pev(ind(k)),k=1,nev_calc)
#endif
        close(EIGVAL_IUNIT_BASE)

    end subroutine export_eigval

    subroutine import_eigval(nev_calc,eigval,pev,ind)
        include 'mpif.h'
        integer, intent(in) :: nev_calc
        real(dp), intent(out) :: eigval(nev_calc), pev(nev_calc)
        integer, intent(out) :: ind(nev_calc)

        integer :: i, iunit, nev_calc_read
        if (node.eq.0) then
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
                write(6,*) "import_eigval: nev_calc mismatch."
            endif
            call die
            return
        endif

        do i=1,nev_calc
#ifdef DATA_PLAINTEXT
            read(iunit,*) ind(i), eigval(i), pev(i)
#else
            read(iunit) ind(i), eigval(i), pev(i)
#endif
        enddo
        close(iunit)
        endif

        call MPI_Bcast(nev_calc_read,1,mpi_integer,0,comm,ierr)
        call MPI_Bcast(ind,nev_calc,mpi_integer,0,comm,ierr)
        call MPI_Bcast(eigval,nev_calc,mpi_double_precision,0,comm,ierr)
        call MPI_Bcast(pev,nev_calc,mpi_double_precision,0,comm,ierr)

    end subroutine import_eigval

    subroutine import_eigvec( isector, iev, node, nloc, eigvec )
        integer, intent(in) :: isector, iev, node
        integer, intent(out) :: nloc
        real(dp), allocatable, intent(out) :: eigvec(:)

        ! local variables
        logical :: file_exist
        integer :: iunit, i
        character(len=150) :: fname

        call eigvec_filename(isector,iev,node,fname)
        iunit = EIGVEC_IUNIT_BASE + node

        inquire(file=fname,exist=file_exist)

        if (.not.file_exist) then
            write(6,*) "import_eigvec: File not found, ", fname
            call die
            return
        endif
        
#ifdef DATA_PLAINTEXT
        open(unit=iunit,file=fname,status="old",form="formatted")
        read(iunit,"(I)") nloc
        allocate(eigvec(nloc))
        read(iunit,"(F24.16)") (eigvec(i),i=1,nloc)
#else
        open(unit=iunit,file=fname,status="old",form="unformatted")
        read(iunit) nloc
        allocate(eigvec(nloc))
        read(iunit) (eigvec(i),i=1,nloc)
#endif
        close(iunit)

    end subroutine import_eigvec

    subroutine eigvec_filename(isector,iev,inode,fname)
        integer, intent(in) :: isector, iev, inode
        character(len=150), intent(out) :: fname

        write(fname,"(A7,I2.2,A1,I3.3,A1,I4.4)") "eigvec-",isector,"-",iev,"-",inode

    end subroutine eigvec_filename

end module ed_io

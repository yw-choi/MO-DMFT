# 
SIESTA_ARCH     = ifort-mpich
FC              = mpif90
FFLAGS          = -traceback -O3 -fast -no-ipo
LIBS            = -mkl=sequential ./lib/ARPACK/libarpack_OSX.a \
                ./lib/ARPACK/parpack_MPI-OSX.a

################################################################

RANLIB          = ranlib
DEFS            = -DMPI -DDEBUG
#
.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#

# 
SIESTA_ARCH     = ifort-mpich
FC              = mpif90
FFLAGS          = -traceback
LIBS            = -mkl=sequential

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

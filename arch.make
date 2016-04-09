# 
SIESTA_ARCH   = ifort-mpich
FC            = mpif90
FFLAGS        = -traceback -fast -no-ipo -xSSE4.2 \
                -I$(MKLROOT)/include/fftw

MKL_FFTW_PATH = $(MKLROOT)/interfaces/fftw3xf
FFTW          = $(MKL_FFTW_PATH)/libfftw3xf_intel.a
ARPACK        = ./lib/ARPACK/libarpack_OSX.a
PARPACK       = ./lib/ARPACK/parpack_MPI-OSX.a
LIBS          = -mkl=sequential $(PARPACK) $(ARPACK) $(FFTW)

################################################################

RANLIB          = ranlib
DEFS            = -DMPI -DDEBUG -DDATA_PLAINTEXT
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

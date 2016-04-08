# Makefile for siesta
#
.SUFFIXES: .f .F .o .a  .f90 .F90

VPATH=.
#
default: main
#
include arch.make

# DMFT_OBJS = ed_operators.o
DMFT_OBJS = ed_solver.o ed_config.o ed_hamiltonian.o ed_utils.o ed_basis.o

MOD_OBJS = alloc.o memory.o memoryinfo.o sys.o parallel_params.o precision.o timer.o timestamp2.o ionew.o
OBJS = bsd.o main.o $(DMFT_OBJS)

COM_OBJS=$(OBJS) $(SYSOBJ)
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)

##################
### FDF Module ###
##################
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "VPATH=$(VPATH)/fdf" \
		    "FPPFLAGS=$(FPPFLAGS)" module )

# Dependencies
main.o ed_config.o: $(FDF)

main.o: parallel_params.o ed_config.o ed_hamiltonian.o ed_solver.o ionew.o
alloc.o: parallel_params.o precision.o sys.o ionew.o
ed_solver.o: ed_utils.o ed_hamiltonian.o
ed_hamiltonian.o: ed_utils.o ed_basis.o
ed_basis.o: alloc.o
memory.o: memoryinfo.o 
main: $(FDF) $(ALL_OBJS)
	$(FC) -o main.x $(LDFLAGS) $(ALL_OBJS) $(FDF) $(LIBS) 
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f main *.o  *.a
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)
#
%.o:%.mod



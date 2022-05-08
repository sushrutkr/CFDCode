FC=gfortran
# FC=ifort
# FFLAGS=-O3
# $(FC) $(FFLAGS) -o $(FNAME) $(OBJ)

SRC=modules.f90 UTIL_PRE_SIM.f90 UTIL_POST_PROC.f90 TDMA.f90 UTIL_SOLVE.f90 UTIL_IMMERSED_BNDRY.f90 UTIL_BC.f90 UTIL_STATS.f90 MAIN.f90

OBJ=${SRC:.f90=.o}
FNAME=CFDCode

%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

default : $(OBJ)
	$(FC) -o $(FNAME) $(OBJ)

clean:
	@rm *.mod *.o 
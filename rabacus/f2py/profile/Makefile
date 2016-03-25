#!/bin/sh

OPT=

FC= gfortran
#FFLAGS= -g -p -cpp -fbounds-check -ffpe-trap=invalid,zero,overflow -fbacktrace
FFLAGS= -g -p -fpp -cpp -O3 -mtune=native -ftree-vectorize 
FFLAGS= -g -p -cpp -O3 -mtune=native # -fopenmp 
#OPT= -DUseOMP


FC= ifort
#FFLAGS= -cpp -g -p -traceback
#FFLAGS= -cpp -g -p -O3 -openmp 
OPT= -DUseOMP
FFLAGS= -cpp -g -p -xHost -O3 -ipo -funroll-loops -heap-arrays -mcmodel=medium -openmp




OBJS=   types.o \
        zhang_jin.o \
        zhang_jin_f.o \
        slatec.o \
	special_functions.o \
	physical_constants.o \
	utils.o \
	hui_gnedin_97.o \
	chem_cool_rates.o \
	ion_solver.o \
	verner_96.o \
	photo_xsections.o \
	source_plane.o \
	source_background.o \
	slab_base.o \
	slab_bgnd.o


# Targets
#=============================================================================
#
# Implicit Rules
#----------------
#
# $@ = name of target 
# $< = name of first dependency
# $^ = name of all dependencies with duplicates removed
# $? = name of all dependencies newer than the target
# $+ = name of all dependencies w/o duplicates removed
#


APPS= screen test_slab_bgnd


all:$(APPS)


# Just Dummy Reporting
#=============================================================================
screen: Makefile
	@echo
	@echo "FC=     " $(FC)
	@echo "FFLAGS= " $(FFLAGS)
	@echo "OPT=    " $(OPT)
	@echo

#
# Main Urchin applications
#=============================================================================
test_slab_bgnd: $(OBJS) test_slab_bgnd.o
	$(FC) $(FFLAGS) $^ $(OPT)  -o $@


# Implicit Rules
#=============================================================================

%.o: %.F90 
	$(FC) $(FFLAGS) $(OPT) -c $< -o $@

%.o: %.f90 
	$(FC) $(FFLAGS) $(OPT) -c $< -o $@


# Standard Cleaning Targets
#=============================================================================

clean :
	rm -f *.o *.mod 

cleanall :
	rm -f *.o *.mod $(APPS) 

tidy :
	rm -f *~ 
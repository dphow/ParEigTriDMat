# ============================================================================
# Name        : Makefile
# Author      : Daniel Howard
# Version     :
# Copyright   : GBU General Public License 2018
# Description : Makefile for Hello MPI World in Fortran
# ============================================================================

.PHONY: all clean

all: src/ParEigTriDMatScaLA.f90
	mpif90 -O2 -g -o bin/ParEigTriDMatScaLA \
		src/ParEigTriDMatScaLA.f90 \
		src/SolveSecularEq.f90 \
		src/TDQREigensolver.f90

clean:
	rm -f bin/ParEigTriDMatScaLA *.mod

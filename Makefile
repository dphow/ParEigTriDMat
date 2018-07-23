# ============================================================================
# Name        : Makefile
# Author      : Daniel Howard
# Version     :
# Copyright   : GBU General Public License 2018
# Description : Makefile for Hello MPI World in Fortran
# ============================================================================

.PHONY: all clean

all: src/ParEigTriDMatScaLA.f90
	mpifort -O3 -g -o bin/ParEigTriDMatScaLA \
		src/SecularSolver.f90 src/EigenSolver.f90 \
		src/ParEigTriDMatScaLA.f90

src: src/SolveSecularEq.f90 src/TDQREigensolver.f90
	mpifort -O3 -c src/SecularSolver.f90 src/EigenSolver.f90
	
clean:
	rm -f bin/ParEigTriDMatScaLA *.mod
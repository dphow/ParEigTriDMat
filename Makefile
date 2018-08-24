# ============================================================================
# Name        : Makefile
# Author      : Daniel Howard
# Version     :
# Copyright   : GBU General Public License 2018
# Description : Makefile for MPI based Cuppen's Algorithm in Fortran
# ============================================================================

.PHONY: all clean

all: src/ParEigTriDMat.f90
	mpiifort -O3 -g -o bin/ParEigTriDMat \
		src/SecularSolver.f90 src/EigenSolver.f90 \
		src/ParEigTriDMat.f90

src: src/SolveSecularEq.f90 src/TDQREigensolver.f90
	mpiifort -O3 -c -g src/SecularSolver.f90 src/EigenSolver.f90
	
clean:
	rm -f bin/ParEigTriDMat *.mod
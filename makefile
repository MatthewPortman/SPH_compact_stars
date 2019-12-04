# Author: Matthew Portman
# Date: 5/1/18
# Makefile for executable sph_test

all: sph_toy_cuda sph_test

sph_test: sph_test.F90 parameters_test.o
	gfortran -o sph_test sph_test.F90 parameters_test.o #-Wall # Generate the executable 'sph_test'.

sph_toy_cuda: sph_toy_cuda.F90 parameters.o
	pgf90 -o sph_toy_cuda sph_toy_cuda.F90 parameters.o -Mcuda -pgf90libs #-Wall # Generate the executable 'sph_test'.

parameters.o: parameters.F90
	pgf90 -c parameters.F90 -Mcuda -pgf90libs

parameters_test.o: parameters_test.F90
	gfortran -c parameters_test.F90

clean: 
	rm sph_toy_cuda sph_test parameters.mod parameters.o #Clean out the target file when 'make clean' is called.

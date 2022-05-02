#!/bin/sh
#
# This script compiles and executes DISORT
#

rm -f a.out
cat disort_py.py BDREF.f DISOBRDF.f ERRPACK.f LINPAK.f LAPACK.f RDI1MACH.f > code.f
#gfortran -O0 -g -fcheck=all -fdump-core -ffpe-trap=invalid,zero,overflow,underflow,denormal -Wall rte_driver.f90 code.f 
#gfortran -O0 -g -fcheck=all -fdump-core -ffpe-trap=invalid,zero,overflow -Wall rte_driver.f90 code.f 
#gfortran -O3 -g -fcheck=all -fdump-core -fbounds-check -Wall disort_program.f90 code.f -o disort_program.exe
gfortran -O3 -g -fcheck=all -fdump-core -fbounds-check -Wall mydisotest.f90 code.f -o disort_program.exe
chmod u+x ./disort_program.exe
time ./disort_program.exe
rm -f code.f
# disort_driver disort_py.f 


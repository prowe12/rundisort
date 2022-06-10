#!/bin/sh
#
# This script compiles and executes DISORT
#

rm -f a.out
cat disort_driver.f disort_py.f BDREF.f DISOBRDF.f ERRPACK.f LINPAK.f LAPACK.f RDI1MACH.f > code.f
#gfortran -O0 -g -fcheck=all -fdump-core -ffpe-trap=invalid,zero,overflow,underflow,denormal -Wall rte_driver.f90 code.f 
#gfortran -O0 -g -fcheck=all -fdump-core -ffpe-trap=invalid,zero,overflow -Wall rte_driver.f90 code.f 
gfortran -O3 -g -fcheck=all -fdump-core -fbounds-check -Wall mydisotest.f90 code.f -o mydisotest.exe
chmod u+x ./mydisotest.exe
time ./mydisotest.exe
rm -f code.f
# disort_driver disort_py.f 


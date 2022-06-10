#!/bin/sh
#
# This script compiles and executes DISORT
#

rm -f a.out
cat disort4_driver.f disort_py.f BDREF.f DISOBRDF.f ERRPACK.f LINPAK.f LAPACK.f RDI1MACH.f > code.f
#gfortran -O0 -g -fcheck=all -fdump-core -ffpe-trap=invalid,zero,overflow,underflow,denormal -Wall rte_driver.f90 code.f 
#gfortran -O0 -g -fcheck=all -fdump-core -ffpe-trap=invalid,zero,overflow -Wall rte_driver.f90 code.f 
#gfortran -O3 -g -fcheck=all -fdump-core -fbounds-check -Wall disort_program.f90 code.f -o disort_program.exe
gfortran -O3 -g -fcheck=all -fdump-core -fbounds-check -Wall mydisortdrivertest.f90 code.f -o mydisortdrivertest.exe
chmod u+x ./mydisortdrivertest.exe
time ./mydisortdrivertest.exe
rm -f code.f
# disort_driver disort_py.f 


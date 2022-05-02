#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 13:04:15 2019

@author: prowe

Purpose: Run DISORT using inputs specified in a namelist file

Instructions:
    1) You will need to install f90nml:
       $ pip install f90nml
    2) Install rundisort_py
    3) Change the path below as needed for your setup.
    4) Make sure you have the file "sample.nml" in your path.
    5) Results should agree with those from sampleRun.

Copyright 2019 by Penny Rowe and NorthWest Research Associates
"""

# .. Third-party modules
import numpy as np
import f90nml

# .. Provided module
from disort_driver_py import disort_driver


# .. Enter the namelist file below:
NAMELIST = 'sample.nml'


# .. Load in the namelist
nml = f90nml.read(NAMELIST)
ip = nml['disort_nml']


# .. Rearrange PMOM into matrix based on NCLDLYR and NMOM
#    First convert PMOM into a 1x(NMOMxNCLDYR) matrix
PMOM = np.array(ip['PMOM']).reshape(-1, 1)
PMOM = PMOM.reshape(ip['NMOM']+1, ip['NCLDLYR'])

# .. Get lengths of variables
MAXMOM = PMOM.shape[0]-1
MAXULV = len(ip['UTAU'])
MAXCLY = len(ip['DTAUC'])
MAXUMU = len(ip['UMU'])
MAXPHI = len(ip['PHI'])

# .. Call DISORT
RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED, errmsg, errflag = \
  disort_driver(ip['NLYR'], ip['DTAUC'], ip['SSALB'], ip['NMOM'],
                ip['NCLDLYR'], ip['CLDLYR'], PMOM, ip['TEMPER'],
                ip['WVNMLO'], ip['WVNMHI'], ip['USRTAU'], ip['NTAU'],
                ip['UTAU'], ip['NSTR'], ip['USRANG'], ip['NUMU'], ip['UMU'],
                ip['NPHI'], ip['PHI'], ip['IBCND'], ip['FBEAM'], ip['UMU0'],
                ip['PHI0'], ip['FISOT'], ip['LAMBER'], ip['ALBEDO'],
                ip['BTEMP'], ip['TTEMP'], ip['TEMIS'], ip['PLANK'],
                ip['ONLYFL'], ip['HEADER'],
                MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM)

# .. Clean up the string (which is bytes)
(errmsg.decode("utf-8")).rstrip()


# .. Print out the results

print("RFLDN: Diffuse down flux, sfc: ")
print(RFLDN[0])

print(" ")
print("FLUP: Diffuse upward flux, toa")
print(FLUP[1])

print(" ")
print("UU: down radiance(z. ang), sfc")
print(UU[:, 0])

print(" ")
print("UU: up radiance(z. ang), toa")
print(UU[:, 1])

print('Error message is:', errmsg)
print('Error flag is:', errflag)

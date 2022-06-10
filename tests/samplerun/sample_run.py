#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr 2, 2021
Last modified Apr 12, 2022s

@author: prowe
"""


# DISORT test runs

# Modules
import numpy as np

# DISORT modules
import disort_driver_py
import disort4_driver_py


def print_results(RFLDN, FLUP, UU):
    """
    Print the results
    @param RFLDN Diffuse down flux, sfc
    @param FLUP Diffuse upward flux, to
    @param UU Upwelling radiance(zenith ang), at toa
    """

    print("RFLDN: Diffuse down flux, sfc: ")
    print(RFLDN[0])

    print(" ")
    print("FLUP: Diffuse upward flux, toa")
    print(FLUP[1])

    print(" ")
    print("UU: down radiance(z. ang), sfc")
    print(UU[:,0])


    print(" ")
    print("UU: up radiance(z. ang), toa")
    print(UU[:,1])



def test1():
    """Test 1: All DTAUC > 3E-1"""
    
    #&disortinput
    NSTR = 16
    NLYR = 3
    NTAU = 2
    UTAU =  [2.59751963615417E+00,  0.00000000000000E+00]
    NPHI = 1
    PHI =  [0.00000000000000E+00]
    IBCND = 0
    UMU0 =  0.00000000000000E+00
    PHI0 =  5.09752949553775E+01
    ALBEDO =  2.00000000000000E-02
    SSALB =  [3.97345705305091E-01,  0.00000000000000E+00,  0.00000000000000E+00]
    FBEAM =  0.00000000000000E+00
    FISOT =  0.00000000000000E+00
    DTAUC =  [1.79751963461538E+00,  3.00000000000000E-01,  5.00000000000000E-01]
    WVNMLO =  9.49500000000000E+02
    WVNMHI =  9.50500000000000E+02
    USRTAU = 1
    NUMU = 6
    UMU = np.array([-1.00E+00, -5.0E-01, -1.0E-02, 1.0E-02,  5.0E-01,  1.0E+00])
    USRANG = 1
    LAMBER = 1
    TEMIS =  0.00000000000000E+00
    PLANK = 1
    ONLYFL = 0
    TEMPER =  [2.730E+02,  2.6300E+02,  2.530000E+02, 2.4300E+02]
    TTEMP =  2.73000000000000E+02
    BTEMP =  2.43000000000000E+02
    NMOM = 154
    NCLDLYR = 1
    HEADER = ' '
    
    #&PMOMINPUT
    CLDLYR = 1;
    NCLDLYR = 1
    PMOM_CLD =  np.array([[\
      1.00000000000000E+00,  9.48125750204321E-01,  8.77985610583193E-01,
      7.98549223777405E-01,  7.18558236968330E-01,  6.42693810867911E-01,
      5.73678519005854E-01,  5.13010454531465E-01,  4.61163809320369E-01,
      4.17835278053646E-01,  3.82217091815565E-01,  3.53202812514318E-01,
      3.29630262810088E-01,  3.10378127243494E-01,  2.94476419822557E-01,
      2.81109058899849E-01,  2.69634149142829E-01,  2.59550738553612E-01,
      2.50486219276792E-01,  2.42162898636502E-01,  2.34381589710650E-01,
      2.26998358208706E-01,  2.19912700250174E-01,  2.13053710650324E-01,
      2.06373137171499E-01,  1.99837921750955E-01,  1.93426277924288E-01,
      1.87124089445595E-01,  1.80922805867113E-01,  1.74817473501186E-01,
      1.68805945029549E-01,  1.62887921212190E-01,  1.57064197780938E-01,
      1.51336163716527E-01,  1.45705774258871E-01,  1.40175515510757E-01,
      1.34748282156670E-01,  1.29427165008278E-01,  1.24214867003941E-01,
      1.19115812402081E-01,  1.14131678100628E-01,  1.09265088145433E-01,
      1.04518426179420E-01,  9.98938211440475E-02,  9.53931322684937E-02,
      9.10179378236864E-02,  8.67695264463031E-02,  8.26488900299741E-02,
      7.86567217887320E-02,  7.47934137676595E-02,  7.10590600516277E-02,
      6.74534590639774E-02,  6.39761219368410E-02,  6.06262790034489E-02,
      5.74028919069218E-02,  5.43046635622176E-02,  5.13300537839376E-02,
      4.84772920089445E-02,  4.57443953448331E-02,  4.31291827850950E-02,
      4.06292944321368E-02,  3.82422073975969E-02,  3.59652549501067E-02,
      3.37956430388841E-02,  3.17304688120563E-02,  2.97667367246471E-02,
      2.79013765472455E-02,  2.61312583318870E-02,  2.44532086757577E-02,
      2.28640244689388E-02,  2.13604877336703E-02,  1.99393772957090E-02,
      1.85974819796405E-02,  1.73316108250826E-02,  1.61386030715892E-02,
      1.50153370081645E-02,  1.39587384034318E-02,  1.29657877683717E-02,
      1.20335273102900E-02,  1.11590656592144E-02,  1.03395833355393E-02,
      9.57233608480240E-03,  8.85465862063461E-03,  8.18396736635522E-03,
      7.55776223760124E-03,  6.97362891687143E-03,  6.42923925208936E-03,
      5.92235238976360E-03,  5.45081514875977E-03,  5.01256183446597E-03,
      4.60561382945866E-03,  4.22807874112366E-03,  3.87814957034229E-03,
      3.55410328428075E-03,  3.25429969685588E-03,  2.97717974730028E-03,
      2.72126389294159E-03,  2.48514989841363E-03,  2.26751060581564E-03,
      2.06709171683086E-03,  1.88270940171437E-03,  1.71324814312041E-03,
      1.55765791260256E-03,  1.41495165953736E-03,  1.28420258937465E-03,
      1.16454133942007E-03,  1.05515343639590E-03,  9.55276403275152E-04,
      8.64197063736000E-04,  7.81248728693611E-04,  7.05808681537369E-04,
      6.37295338824095E-04,  5.75166007950539E-04,  5.18914184753112E-04,
      4.68067120626412E-04,  4.22183839118623E-04,  3.80852631574527E-04,
      3.43689174437717E-04,  3.10334936817090E-04,  2.80455032813744E-04,
      2.53736927707329E-04,  2.29888946797725E-04,  2.08638550292758E-04,
      1.89731552812827E-04,  1.72931164311424E-04,  1.58016805694613E-04,
      1.44783399392560E-04,  1.33040544441120E-04,  1.22611918635706E-04,
      1.13334614532940E-04,  1.05058902705752E-04,  9.76477054657584E-05,
      9.09761810019623E-05,  8.49312807437148E-05,  7.94111393386365E-05,
      7.43248781145902E-05,  6.95921945644890E-05,  6.51430113324402E-05,
      6.09169142411044E-05,  5.68629435342834E-05,  5.29387527304218E-05,
      4.91102115330458E-05,  4.53509002765885E-05,  4.16412788029823E-05,
      3.79681930158634E-05,  3.43242766269775E-05,  3.07072106394779E-05,
      2.71188751747762E-05,  2.35647473451367E-05,  2.00532694158332E-05,
      1.65953612160422E-05,  1.32037184522940E-05,  9.89202934043430E-06,
      6.67410987669309E-06,  3.56329868017677E-06]]).T
    #/
    
    
    NCLDLYR = PMOM_CLD.shape[1]
    MAXMOM = PMOM_CLD.shape[0]-1
    MAXULV = len(UTAU)
    MAXCLY = len(DTAUC)
    MAXUMU = len(UMU)
    MAXPHI = len(PHI)
    
    
    RFLDIR, \
    RFLDN, \
    FLUP, \
    DFDT, \
    UAVG, \
    UU, \
    ALBMED, \
    TRNMED, \
    ERRMSG, \
    ERRFLAG = disort_driver_py.disort_driver(NLYR, DTAUC, SSALB, NMOM, NCLDLYR,
                                             CLDLYR, PMOM_CLD, TEMPER, WVNMLO,
                                             WVNMHI, USRTAU, NTAU, UTAU, NSTR,
                                             USRANG, NUMU, UMU, NPHI, PHI, IBCND,
                                             FBEAM, UMU0, PHI0, FISOT, LAMBER,
                                             ALBEDO, BTEMP, TTEMP, TEMIS, PLANK,
                                             ONLYFL, HEADER, MAXCLY, MAXULV,
                                             MAXUMU, MAXPHI, MAXMOM)
    print_results(RFLDN, FLUP, UU)
    
    
    RFLDIR, \
    RFLDN, \
    FLUP, \
    DFDT, \
    UAVG, \
    UU, \
    ALBMED, \
    TRNMED, \
    ERRMSG, \
    ERRFLAG = disort4_driver_py.disort_driver(NLYR, NMOM, NSTR,
                                    NUMU, NPHI, NTAU,
                                    USRANG, USRTAU, IBCND, ONLYFL,
                                    PLANK, LAMBER, DTAUC, SSALB,
                                    NCLDLYR, CLDLYR, PMOM_CLD, TEMPER,
                                    WVNMLO, WVNMHI,
                                    UTAU, UMU0, PHI0, UMU, PHI, FBEAM,
                                    FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,
                                    HEADER, MAXCLY, MAXULV,
                                    MAXUMU, MAXPHI, MAXMOM )
    print_results(RFLDN, FLUP, UU)
    
    
    # Expected results:
    
    # RFLDN: Diffuse down flux, sfc:
    # 0.13878909
    
    # FLUP: Diffuse upward flux, toa
    # 0.1894969
    
    # UU: down radiance(z. ang), sfc
    # [[0.04285485]
    #  [0.04521359]
    #  [0.037155  ]
    #  [0.03711424]
    #  [0.03711424]
    #  [0.03711424]]
    
    # UU: up radiance(z. ang), toa
    # [[0.        ]
    #  [0.        ]
    #  [0.        ]
    #  [0.05801807]
    #  [0.0624982 ]
    #  [0.0567358 ]]



# .. Test 2: Some DTAUC < 1E-5
DTAUC =  [1.8E+00,  3.0E-8,  5.00000000000000E-01]
UTAU =  [sum(DTAUC),  0.0000E+00]

RFLDIR, \
RFLDN, \
FLUP, \
DFDT, \
UAVG, \
UU, \
ALBMED, \
TRNMED, \
ERRMSG, \
ERRFLAG = disort_driver_py.disort_driver(NLYR, DTAUC, SSALB, NMOM, NCLDLYR,
                                          CLDLYR, PMOM_CLD, TEMPER, WVNMLO,
                                          WVNMHI, USRTAU, NTAU, UTAU, NSTR,
                                          USRANG, NUMU, UMU, NPHI, PHI, IBCND,
                                          FBEAM, UMU0, PHI0, FISOT, LAMBER,
                                          ALBEDO, BTEMP, TTEMP, TEMIS, PLANK,
                                          ONLYFL, HEADER, MAXCLY, MAXULV,
                                          MAXUMU, MAXPHI, MAXMOM)
print_results(RFLDN, FLUP, UU)

del RFLDN, FLUP, UU

UTAU =  [sum(DTAUC),  0.0000E+00]

RFLDIR, \
RFLDN, \
FLUP, \
DFDT, \
UAVG, \
UU, \
ALBMED, \
TRNMED, \
ERRMSG, \
ERRFLAG = disort4_driver_py.disort_driver(NLYR, NMOM, NSTR,
                                NUMU, NPHI, NTAU,
                                USRANG, USRTAU, IBCND, ONLYFL,
                                PLANK, LAMBER, DTAUC, SSALB,
                                NCLDLYR, CLDLYR, PMOM_CLD, TEMPER,
                                WVNMLO, WVNMHI,
                                UTAU, UMU0, PHI0, UMU, PHI, FBEAM,
                                FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,
                                HEADER, MAXCLY, MAXULV,
                                MAXUMU, MAXPHI, MAXMOM )
print_results(RFLDN, FLUP, UU)


# UU: down radiance(z. ang), sfc
#     3e-4
# [[0.04131474]
#  [0.0457374 ]
#  [0.03715499]
#  [0.03710805]
#  [0.03710718]
#  [0.03710717]]

# 3e-5
# UU: down radiance(z. ang), sfc
# [[0.04132137]
#  [0.0457426 ]
#  [0.037155  ]
#  [0.03710745]
#  [0.03710745]
#  [0.03710745]]

# 3e-6
# UU: down radiance(z. ang), sfc
# [[0.04125511]
#  [0.0458356 ]
#  [0.03715499]
#  [0.03710451]
#  [0.03710365]
#  [0.03710364]]

# 3e-7
# UU: down radiance(z. ang), sfc
# [[0.04221911]
#  [0.04518124]
#  [0.037155  ]
#  [0.03706822]
#  [0.03706822]
#  [0.03706822]]

# 3e-8
# UU: down radiance(z. ang), sfc
# [[0.04126661]
#  [0.04560908]
#  [0.037155  ]
#  [0.03787563]
#  [0.03787563]
#  [0.03787563]]

# 3e-9
# UU: down radiance(z. ang), sfc
# [[0.04143301]
#  [0.04582589]
#  [0.037155  ]
#  [0.03458402]
#  [0.03458402]
#  [0.03458402]]

# 3e-11
# UU: down radiance(z. ang), sfc
# [[0.03379185]
#  [0.04718457]
#  [0.037155  ]
#  [0.11340915]
#  [0.11340915]
#  [0.11340915]]
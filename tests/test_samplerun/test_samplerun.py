#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr 2, 2021
Last modified Apr 12, 2022s

@author: prowe

Copyright Penny M. Rowe and NorthWest Research Associates
"""


# DISORT test runs

# Modules
import numpy as np
import pytest

# DISORT modules
from rundisort import disort_driver_py

# Setup and teardown: the fixtures
@pytest.fixture
def get_disortinput():
    """
    Get inputs to DISORT as a dictionary
    @yield disortinput
    """

    disortinput = {
    'NSTR': 16,
    'NLYR': 3,
    'NTAU': 2,
    'NPHI': 1,
    'PHI': [0.00000000000000E+00],
    'IBCND': 0,
    'UMU0': 0.00000000000000E+00,
    'PHI0': 5.09752949553775E+01,
    'ALBEDO': 2.00000000000000E-02,
    'SSALB': [3.97345705305091E-01, 0.00000000000000E+00, 0.00000000000000E+00],
    'FBEAM': 0.00000000000000E+00,
    'FISOT': 0.00000000000000E+00,
    'WVNMLO': 949.50,
    'WVNMHI': 950.50,
    'USRTAU': 1,
    'NUMU': 6,
    'UMU': np.array([-1.00E+00, -5.0E-01, -1.0E-02, 1.0E-02, 5.0E-01, 1.0E+00]),
    'USRANG': 1,
    'LAMBER': 1,
    'TEMIS': 0.00000000000000E+00,
    'PLANK': 1,
    'ONLYFL': 0,
    'TEMPER': [2.730E+02, 2.6300E+02, 2.530000E+02, 2.4300E+02],
    'TTEMP': 2.73000000000000E+02,
    'BTEMP': 2.43000000000000E+02,
    'NMOM': 154,
    'NCLDLYR': 1,
    'HEADER': ' ',
    'CLDLYR': 1,
    'NCLDLYR': 1,
    'PMOM_CLD': np.array([[\
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
      6.67410987669309E-06,  3.56329868017677E-06]]).T,
    }

    yield disortinput


def test1(get_disortinput):
    """Test 1: All DTAUC > 3E-1"""

    disortinput = get_disortinput
    DTAUC = [1.79751963461538E+00, 3.00000000000000E-01, 5.00000000000000E-01]
    UTAU = [sum(DTAUC),  0.0000E+00]
    #UTAU = [2.59751963615417E+00, 0.00000000000000E+00]

    NCLDLYR = disortinput['PMOM_CLD'].shape[1]
    MAXMOM = disortinput['PMOM_CLD'].shape[0]-1
    MAXULV = len(UTAU)
    MAXCLY = len(DTAUC)
    MAXUMU = len(disortinput['UMU'])
    MAXPHI = len(disortinput['PHI'])

    RFLDIR, \
    RFLDN, \
    FLUP, \
    DFDT, \
    UAVG, \
    UU, \
    ALBMED, \
    TRNMED, \
    ERRMSG, \
    ERRFLAG = disort_driver_py.disort_driver(disortinput['NLYR'],
                                             DTAUC,
                                             disortinput['SSALB'],
                                             disortinput['NMOM'],
                                             NCLDLYR,
                                             disortinput['CLDLYR'],
                                             disortinput['PMOM_CLD'],
                                             disortinput['TEMPER'],
                                             disortinput['WVNMLO'],
                                             disortinput['WVNMHI'],
                                             disortinput['USRTAU'],
                                             disortinput['NTAU'],
                                             UTAU,
                                             disortinput['NSTR'],
                                             disortinput['USRANG'],
                                             disortinput['NUMU'],
                                             disortinput['UMU'],
                                             disortinput['NPHI'],
                                             disortinput['PHI'],
                                             disortinput['IBCND'],
                                             disortinput['FBEAM'],
                                             disortinput['UMU0'],
                                             disortinput['PHI0'],
                                             disortinput['FISOT'],
                                             disortinput['LAMBER'],
                                             disortinput['ALBEDO'],
                                             disortinput['BTEMP'],
                                             disortinput['TTEMP'],
                                             disortinput['TEMIS'],
                                             disortinput['PLANK'],
                                             disortinput['ONLYFL'],
                                             disortinput['HEADER'],
                                             MAXCLY, MAXULV, MAXUMU,
                                             MAXPHI, MAXMOM)

    # RFLDN: Diffuse down flux, sfc:
    assert np.isclose(RFLDN[0], 0.13878909)

    # FLUP: Diffuse upward flux, toa
    assert np.isclose(FLUP[1], 0.1894969)

    # UU[:,0]: down radiance(z. ang), sfc
    assert np.allclose(UU[:,0], [[0.04285485],
                                 [0.04521359],
                                 [0.037155  ],
                                 [0.03711424],
                                 [0.03711424],
                                 [0.03711424]])

    # UU[:,1]: up radiance(z. ang), toa
    assert np.allclose(UU[:,1], [[0.],
                                 [0.],
                                 [0.],
                                 [0.05801807],
                                 [0.0624982 ],
                                 [0.0567358 ]])


def test_3em3(get_disortinput):
    """
    Second DTAUC = 3E-3
    This optical depth is small compared to the other ods (1.8 and 0.5),
    Here we test that the answers are as expected, and also show that the
    absolute differences for the DISORT output for this case compared to
    the case below with the 2nd DTAUC = 3E-4 are within
    1E-4 W/m2 for fluxes,
    and within 2E-5 W/(m2 str-1 cm-1) for radiances.
    e.g. 0.1 mW/m2 for fluxes, and 0.2 mW/(m2 str-1 cm-1) for radiances.
    Because we should not be running in to precision errors yet, 
    this value, 0.2 mW/(m2 str-1 cm-1), is  expected to be the difference 
    due to increasing the 2nd DTAUC from 0.0003 to 0.003. Therefore,
    to "zero-out" the optical depth in a level, use an optical depth <= 0.0003
    """

    DTAUC = [1.8E+00, 3.0E-3, 5.0E-01]

    UTAU = [sum(DTAUC), 0.0000E+00]
    MAXULV = len(UTAU)
    MAXCLY = len(DTAUC)

    disortinput = get_disortinput
    NCLDLYR = disortinput['PMOM_CLD'].shape[1]
    MAXMOM = disortinput['PMOM_CLD'].shape[0]-1
    MAXUMU = len(disortinput['UMU'])
    MAXPHI = len(disortinput['PHI'])

    RFLDIR, \
    RFLDN, \
    FLUP, \
    DFDT, \
    UAVG, \
    UU, \
    ALBMED, \
    TRNMED, \
    ERRMSG, \
    ERRFLAG = disort_driver_py.disort_driver(disortinput['NLYR'],
                                             DTAUC,
                                             disortinput['SSALB'],
                                             disortinput['NMOM'],
                                             NCLDLYR,
                                             disortinput['CLDLYR'],
                                             disortinput['PMOM_CLD'],
                                             disortinput['TEMPER'],
                                             disortinput['WVNMLO'],
                                             disortinput['WVNMHI'],
                                             disortinput['USRTAU'],
                                             disortinput['NTAU'],
                                             UTAU,
                                             disortinput['NSTR'],
                                             disortinput['USRANG'],
                                             disortinput['NUMU'],
                                             disortinput['UMU'],
                                             disortinput['NPHI'],
                                             disortinput['PHI'],
                                             disortinput['IBCND'],
                                             disortinput['FBEAM'],
                                             disortinput['UMU0'],
                                             disortinput['PHI0'],
                                             disortinput['FISOT'],
                                             disortinput['LAMBER'],
                                             disortinput['ALBEDO'],
                                             disortinput['BTEMP'],
                                             disortinput['TTEMP'],
                                             disortinput['TEMIS'],
                                             disortinput['PLANK'],
                                             disortinput['ONLYFL'],
                                             disortinput['HEADER'],
                                             MAXCLY, MAXULV, MAXUMU,
                                             MAXPHI, MAXMOM)

    # RFLDN: Diffuse down flux, sfc:
    assert np.isclose(RFLDN[0], 0.13770026)

    # FLUP[1]: Diffuse upward flux, toa
    assert np.isclose(FLUP[1], 0.18716988)

    # UU[:,0]: down radiance(z. ang), sfc
    assert np.allclose(UU[:,0], [[0.04133346],
                                 [0.04573157],
                                 [0.03715499],
                                 [0.03710731],
                                 [0.03710731],
                                 [0.03710731]])

    # UU[:,1]: up radiance(z. ang), toa
    assert np.allclose(UU[:,1], [[0.        ],
                                 [0.        ],
                                 [0.        ],
                                 [0.05799797],
                                 [0.06197565],
                                 [0.05564552]])

    # RFLDN: Diffuse down flux, sfc compared to result with 2nd DTAUC=3e-4
    assert np.isclose(RFLDN[0]-.13767686, 0, atol=1e-4)

    # FLUP: Diffuse upward flux, toa compared to result with 2nd DTAUC=3e-4
    assert np.isclose(FLUP[1]- 0.18715, 0, atol=1e-4)

    # UU[:,0]: down rad(z. ang), sfc compared to result with 2nd DTAUC=3e-4
    assert np.allclose(UU[0,0] - 0.04131474, 0, atol=2e-5)
    assert np.allclose(UU[1,0] - 0.04573740, 0, atol=1e-5)
    assert np.allclose(UU[2,0] - 0.03715499, 0, atol=1e-6)
    assert np.allclose(UU[3,0] - 0.03710805, 0, atol=1e-6)
    assert np.allclose(UU[4,0] - 0.03710718, 0, atol=1e-6)
    assert np.allclose(UU[5,0] - 0.03710717, 0, atol=1e-6)

    # UU[:,1]: up radiance(z. ang), toa compared to result with 2nd DTAUC=3e-4
    assert np.allclose(UU[0,1], 0, atol=1e-5)
    assert np.allclose(UU[1,1], 0, atol=1e-5)
    assert np.allclose(UU[2,1], 0, atol=1e-6)
    assert np.allclose(UU[3,1] - 0.05799777, 0, atol=1e-6)
    assert np.allclose(UU[4,1] - 0.06197017, 0, atol=1e-5)
    assert np.allclose(UU[5,1] - 0.05563463, 0, atol=1e-4)


def test_3em4(get_disortinput):
    """
    Second DTAUC = 3E-4
    This optical depth is very small compared to the other ods (1.8 and 0.5),
    and so we assume it is not very different from zero.
    Thus decreasing it further should make a negligible difference.
    However, as the next tests show, instead the difference grows as
    DTAUC shrinks, due to the single precision.
    """

    DTAUC = [1.8E+00, 3.0E-4, 5.0E-01]

    UTAU = [sum(DTAUC), 0.0000E+00]
    MAXULV = len(UTAU)
    MAXCLY = len(DTAUC)

    disortinput = get_disortinput
    NCLDLYR = disortinput['PMOM_CLD'].shape[1]
    MAXMOM = disortinput['PMOM_CLD'].shape[0]-1
    MAXUMU = len(disortinput['UMU'])
    MAXPHI = len(disortinput['PHI'])

    RFLDIR, \
    RFLDN, \
    FLUP, \
    DFDT, \
    UAVG, \
    UU, \
    ALBMED, \
    TRNMED, \
    ERRMSG, \
    ERRFLAG = disort_driver_py.disort_driver(disortinput['NLYR'],
                                             DTAUC,
                                             disortinput['SSALB'],
                                             disortinput['NMOM'],
                                             NCLDLYR,
                                             disortinput['CLDLYR'],
                                             disortinput['PMOM_CLD'],
                                             disortinput['TEMPER'],
                                             disortinput['WVNMLO'],
                                             disortinput['WVNMHI'],
                                             disortinput['USRTAU'],
                                             disortinput['NTAU'],
                                             UTAU,
                                             disortinput['NSTR'],
                                             disortinput['USRANG'],
                                             disortinput['NUMU'],
                                             disortinput['UMU'],
                                             disortinput['NPHI'],
                                             disortinput['PHI'],
                                             disortinput['IBCND'],
                                             disortinput['FBEAM'],
                                             disortinput['UMU0'],
                                             disortinput['PHI0'],
                                             disortinput['FISOT'],
                                             disortinput['LAMBER'],
                                             disortinput['ALBEDO'],
                                             disortinput['BTEMP'],
                                             disortinput['TTEMP'],
                                             disortinput['TEMIS'],
                                             disortinput['PLANK'],
                                             disortinput['ONLYFL'],
                                             disortinput['HEADER'],
                                             MAXCLY, MAXULV, MAXUMU,
                                             MAXPHI, MAXMOM)

    # RFLDN: Diffuse down flux, sfc:
    assert np.isclose(RFLDN[0], 0.13767686)

    # FLUP: Diffuse upward flux, toa
    assert np.isclose(FLUP[1], 0.18715039)

    # UU[:,0]: down radiance(z. ang), sfc
    assert np.allclose(UU[:,0], [[0.04131474],
                                 [0.0457374 ],
                                 [0.03715499],
                                 [0.03710805],
                                 [0.03710718],
                                 [0.03710717]])

    # UU[:,1]: up radiance(z. ang), toa
    assert np.allclose(UU[:,1], [[0.        ],
                                 [0.        ],
                                 [0.        ],
                                 [0.05799777],
                                 [0.06197017],
                                 [0.05563463]])


def test_3em5(get_disortinput):
    """
    Second DTAUC = 3E-5
    When compared to the results with the second dTAUC = 3E-4,
    here we see that all absolute differences for the DISORT output
    are within 1E-4 W/m2 for fluxes, and 1E-5 W/(m2 str-1 cm-1) for radiances.
    e.g. 0.1 mW/m2 for fluxes, and 0.01 mW/(m2 str-1 cm-1) for radiances.
    Since 0.01 mW/(m2 str-1 cm-1) is below the typical error level, 
    3E-5 seems like a reasonable value to use to zero-out optical depths
    """

    DTAUC = [1.8E+00,  3.0E-5,  5.0E-01]

    UTAU = [sum(DTAUC),  0.0000E+00]
    MAXULV = len(UTAU)
    MAXCLY = len(DTAUC)

    disortinput = get_disortinput
    NCLDLYR = disortinput['PMOM_CLD'].shape[1]
    MAXMOM = disortinput['PMOM_CLD'].shape[0]-1
    MAXUMU = len(disortinput['UMU'])
    MAXPHI = len(disortinput['PHI'])

    RFLDIR, \
    RFLDN, \
    FLUP, \
    DFDT, \
    UAVG, \
    UU, \
    ALBMED, \
    TRNMED, \
    ERRMSG, \
    ERRFLAG = disort_driver_py.disort_driver(disortinput['NLYR'],
                                             DTAUC,
                                             disortinput['SSALB'],
                                             disortinput['NMOM'],
                                             NCLDLYR,
                                             disortinput['CLDLYR'],
                                             disortinput['PMOM_CLD'],
                                             disortinput['TEMPER'],
                                             disortinput['WVNMLO'],
                                             disortinput['WVNMHI'],
                                             disortinput['USRTAU'],
                                             disortinput['NTAU'],
                                             UTAU,
                                             disortinput['NSTR'],
                                             disortinput['USRANG'],
                                             disortinput['NUMU'],
                                             disortinput['UMU'],
                                             disortinput['NPHI'],
                                             disortinput['PHI'],
                                             disortinput['IBCND'],
                                             disortinput['FBEAM'],
                                             disortinput['UMU0'],
                                             disortinput['PHI0'],
                                             disortinput['FISOT'],
                                             disortinput['LAMBER'],
                                             disortinput['ALBEDO'],
                                             disortinput['BTEMP'],
                                             disortinput['TTEMP'],
                                             disortinput['TEMIS'],
                                             disortinput['PLANK'],
                                             disortinput['ONLYFL'],
                                             disortinput['HEADER'],
                                             MAXCLY, MAXULV, MAXUMU,
                                             MAXPHI, MAXMOM)

    # RFLDN: Diffuse down flux, sfc:
    assert np.isclose(RFLDN[0]-.13767686, 0, atol=1e-4)

    # FLUP: Diffuse upward flux, toa
    assert np.isclose(FLUP[1]- 0.18715039, 0, atol=1e-4)

    # UU[:,0]: down radiance(z. ang), sfc
    assert np.allclose(UU[0,0] - 0.04131474, 0, atol=1e-5)
    assert np.allclose(UU[1,0] - 0.04573740, 0, atol=1e-5)
    assert np.allclose(UU[2,0] - 0.03715499, 0, atol=1e-6)
    assert np.allclose(UU[3,0] - 0.03710805, 0, atol=1e-6)
    assert np.allclose(UU[4,0] - 0.03710718, 0, atol=1e-6)
    assert np.allclose(UU[5,0] - 0.03710717, 0, atol=1e-6)

    # UU[:,1]: up radiance(z. ang), toa
    assert np.allclose(UU[0,1], 0, atol=1e-5)
    assert np.allclose(UU[1,1], 0, atol=1e-5)
    assert np.allclose(UU[2,1], 0, atol=1e-6)
    assert np.allclose(UU[3,1] - 0.05799777, 0, atol=1e-6)
    assert np.allclose(UU[4,1] - 0.06197017, 0, atol=1e-5)
    assert np.allclose(UU[5,1] - 0.05563463, 0, atol=1e-4)



def test_3em6(get_disortinput):
    """
    Second DTAUC = 3E-6
    When compared to the results with the second dTAUC = 3E-5,
    here we see that all absolute differences for the DISORT output
    are within 1E-3 W/m2 for fluxes, and 1E-4 W/(m2 str-1 cm-1) for radiances.
    e.g. 1 mW/m2 for fluxes, and 0.1 mW/(m2 str-1 cm-1) for radiances.
    Given how small the optical depths are, these differences seem likely to
    be due to increasing round-off error rather than the decreasing optical
    depths. Thus we keep the optical depths above 1E-5
    """

    DTAUC = [1.8E+00,  3.0E-6,  5.0E-01]

    UTAU = [sum(DTAUC),  0.0000E+00]
    MAXULV = len(UTAU)
    MAXCLY = len(DTAUC)

    disortinput = get_disortinput
    NCLDLYR = disortinput['PMOM_CLD'].shape[1]
    MAXMOM = disortinput['PMOM_CLD'].shape[0]-1
    MAXUMU = len(disortinput['UMU'])
    MAXPHI = len(disortinput['PHI'])

    RFLDIR, \
    RFLDN, \
    FLUP, \
    DFDT, \
    UAVG, \
    UU, \
    ALBMED, \
    TRNMED, \
    ERRMSG, \
    ERRFLAG = disort_driver_py.disort_driver(disortinput['NLYR'],
                                             DTAUC,
                                             disortinput['SSALB'],
                                             disortinput['NMOM'],
                                             NCLDLYR,
                                             disortinput['CLDLYR'],
                                             disortinput['PMOM_CLD'],
                                             disortinput['TEMPER'],
                                             disortinput['WVNMLO'],
                                             disortinput['WVNMHI'],
                                             disortinput['USRTAU'],
                                             disortinput['NTAU'],
                                             UTAU,
                                             disortinput['NSTR'],
                                             disortinput['USRANG'],
                                             disortinput['NUMU'],
                                             disortinput['UMU'],
                                             disortinput['NPHI'],
                                             disortinput['PHI'],
                                             disortinput['IBCND'],
                                             disortinput['FBEAM'],
                                             disortinput['UMU0'],
                                             disortinput['PHI0'],
                                             disortinput['FISOT'],
                                             disortinput['LAMBER'],
                                             disortinput['ALBEDO'],
                                             disortinput['BTEMP'],
                                             disortinput['TTEMP'],
                                             disortinput['TEMIS'],
                                             disortinput['PLANK'],
                                             disortinput['ONLYFL'],
                                             disortinput['HEADER'],
                                             MAXCLY, MAXULV, MAXUMU,
                                             MAXPHI, MAXMOM)

    # RFLDN: Diffuse down flux, sfc:
    assert np.isclose(RFLDN[0]-0.1377229, 0, atol=1e-3)

    # FLUP: Diffuse upward flux, toa
    assert np.isclose(FLUP[1]- 0.18712741, 0, atol=1e-3)

    # UU[:,0]: down radiance(z. ang), sfc
    assert np.allclose(UU[0,0] - 0.04132137, 0, atol=1e-4)
    assert np.allclose(UU[1,0] - 0.0457426, 0, atol=1e-4)
    assert np.allclose(UU[2,0] - 0.037155, 0, atol=1e-6)
    assert np.allclose(UU[3,0] - 0.03710745, 0, atol=1e-5)
    assert np.allclose(UU[4,0] - 0.03710745, 0, atol=1e-5)
    assert np.allclose(UU[5,0] - 0.03710745, 0, atol=1e-5)

    # UU[:,1]: up radiance(z. ang), toa
    assert np.allclose(UU[0,1], 0, atol=1e-5)
    assert np.allclose(UU[1,1], 0, atol=1e-5)
    assert np.allclose(UU[2,1], 0, atol=1e-6)
    assert np.allclose(UU[3,1] - 0.05799777, 0, atol=1e-5)
    assert np.allclose(UU[4,1] - 0.06197017, 0, atol=1e-4)
    assert np.allclose(UU[5,1] - 0.05563463, 0, atol=1e-4)


def test_3em8(get_disortinput):
    """
    2nd DTAUC = 3E-8
    This test demonstrates how errors explode for DTAUC << 1e-5
    """

    DTAUC = [1.8E+00,  3.0E-8,  5.0E-01]

    UTAU = [sum(DTAUC),  0.0000E+00]
    MAXULV = len(UTAU)
    MAXCLY = len(DTAUC)

    disortinput = get_disortinput
    NCLDLYR = disortinput['PMOM_CLD'].shape[1]
    MAXMOM = disortinput['PMOM_CLD'].shape[0]-1
    MAXUMU = len(disortinput['UMU'])
    MAXPHI = len(disortinput['PHI'])

    RFLDIR, \
    RFLDN, \
    FLUP, \
    DFDT, \
    UAVG, \
    UU, \
    ALBMED, \
    TRNMED, \
    ERRMSG, \
    ERRFLAG = disort_driver_py.disort_driver(disortinput['NLYR'],
                                             DTAUC,
                                             disortinput['SSALB'],
                                             disortinput['NMOM'],
                                             NCLDLYR,
                                             disortinput['CLDLYR'],
                                             disortinput['PMOM_CLD'],
                                             disortinput['TEMPER'],
                                             disortinput['WVNMLO'],
                                             disortinput['WVNMHI'],
                                             disortinput['USRTAU'],
                                             disortinput['NTAU'],
                                             UTAU,
                                             disortinput['NSTR'],
                                             disortinput['USRANG'],
                                             disortinput['NUMU'],
                                             disortinput['UMU'],
                                             disortinput['NPHI'],
                                             disortinput['PHI'],
                                             disortinput['IBCND'],
                                             disortinput['FBEAM'],
                                             disortinput['UMU0'],
                                             disortinput['PHI0'],
                                             disortinput['FISOT'],
                                             disortinput['LAMBER'],
                                             disortinput['ALBEDO'],
                                             disortinput['BTEMP'],
                                             disortinput['TTEMP'],
                                             disortinput['TEMIS'],
                                             disortinput['PLANK'],
                                             disortinput['ONLYFL'],
                                             disortinput['HEADER'],
                                             MAXCLY, MAXULV, MAXUMU,
                                             MAXPHI, MAXMOM)


    # THE ACTUAL VALUES - note that the down flux no longer makes sense
    # RFLDN: Diffuse down flux, sfc:
    assert np.isclose(RFLDN[0], .25838745)       # compare to 0.13767686

    # FLUP: Diffuse upward flux, toa
    assert np.isclose(FLUP[1], 0.14003532)       # compare to 0.18715039

    # UU[:,0]: down radiance(z. ang), sfc        # compare to
    assert np.allclose(UU[:,0], [[0.04126661],   # [0.04131474]
                                 [0.04560908],   # [0.0457374 ]
                                 [0.037155  ],   # [0.03715499]
                                 [0.03787563],   # [0.03710805]
                                 [0.03787563],   # [0.03710718]
                                 [0.03787563]])  # [0.03710717]

    # UU[:,1]: up radiance(z. ang), toa          # compare to
    assert np.allclose(UU[:,1], [[0.        ],   # 0
                                 [0.        ],   # 0
                                 [0.        ],   # 0
                                 [0.05766603],   # [0.05799777]
                                 [0.05684418],   # [0.06197017]
                                 [0.04508059]])  # [0.05563463]


    # COMPARISON TO EXAMPLE WITH SECOND DTAUC = 3E-4
    # RFLDN: Diffuse down flux, sfc:
    assert np.isclose(RFLDN[0] -.13767686, 0, atol=0.2)

    # FLUP: Diffuse upward flux, toa
    assert np.isclose(FLUP[1] - 0.18715039, 0, atol=0.1)

    # UU[:,0]: down radiance(z. ang), sfc
    assert np.allclose(UU[0,0] - 0.04131474, 0, atol=1e-4)
    assert np.allclose(UU[1,0] - 0.04573740, 0, atol=1e-3)
    assert np.allclose(UU[2,0] - 0.03715499, 0, atol=1e-6)
    assert np.allclose(UU[3,0] - 0.03710805, 0, atol=1e-3)
    assert np.allclose(UU[4,0] - 0.03710718, 0, atol=1e-3)
    assert np.allclose(UU[5,0] - 0.03710717, 0, atol=1e-3)

    # UU[:,1]: up radiance(z. ang), toa
    assert np.allclose(UU[0,1], 0, atol=1e-5)
    assert np.allclose(UU[1,1], 0, atol=1e-5)
    assert np.allclose(UU[2,1], 0, atol=1e-6)
    assert np.allclose(UU[3,1] - 0.05799777, 0, atol=1e-3)
    assert np.allclose(UU[4,1] - 0.06197017, 0, atol=1e-2)
    assert np.allclose(UU[5,1] - 0.05563463, 0, atol=0.1)

c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: DISOTEST.f,v 2.1 2000/04/03 21:21:55 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      PROGRAM  TESTDO

c    Runs test problems for DISORT and checks answers. These
c    problems test almost all logical branches in DISORT.

c    As distributed, runs all problems.  To run just a subset,
c    change DATA statement for DOPROB;  e.g. to run just problem 4,
c    set DOPROB(4) TRUE and all others FALSE.

c    It is HIGHLY recommended that you use the code below as a template
c    for creating your own CALLs to DISORT, rather than starting from
c    scratch.  This will prevent a lot of mistakes and ensure that every
c    input argument gets a value.  Note in particular how GETMOM is
c    sometimes called to fill an array section of PMOM (for one layer);
c    several people have done this incorrectly in attempting to write it
c    ab initio (passing array sections for arrays that do not start at
c    element 1 is tricky).

c    Note that the ratio to the 'correct answer' may occasionally be
c    significantly different from unity -- even so different that
c    the ratio just prints as ****** rather than a number.  However,
c    this mostly occurs for values of flux or intensity that are very
c    small compared to the forcing functions (that is, small compared
c    to internal thermal emission and/or radiation incident at the
c    boundaries).  The printed number 'SERIOUSLY NON-UNIT RATIOS'
c    attempts to count just the cases where there is a real disagreement
c    and not those where quantitites are down at their noise level
c    (defined as 10^(-6) times their maximum value).

c    Further documentation can be found in the file DISOTEST.doc.


c  Routines called :

c    DISORT:   The discrete ordinates radiative transfer program

c    BDREF:    Sets bidirectional reflectance of lower boundary

c    GETMOM:   Sets phase function Legendre coefficients

c    PRTFIN:   Prints fluxes and intensities and their ratios to
c              the correct values

c    CHEKDO:   Data block containing correct fluxes and intensities

c    RATIO :   Ratio of calculated to correct value with underflow
c              and overflow protection (kept in file DISORT.f)

c+---------------------------------------------------------------------+
c
c                 ** DISORT I/O specifications **

      INTEGER  MAXCLY, MAXMOM, MAXPHI, MAXULV, MAXUMU
      PARAMETER ( MAXCLY = 6, MAXMOM = 299, MAXPHI = 3, MAXULV = 5,
     &            MAXUMU = 10 )
      CHARACTER  HEADER*127
      LOGICAL  LAMBER, PLANK, ONLYFL, PRNT(5), USRANG, USRTAU
      INTEGER  IBCND, NMOM, NLYR, NUMU, NSTR, NPHI, NTAU
      REAL     ACCUR, ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,
     &         PHI( MAXPHI ), PMOM( 0:MAXMOM, MAXCLY ),
     &         PHI0, SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,
     &         WVNMLO, WVNMHI, UMU( MAXUMU ), UMU0, UTAU( MAXULV )

      REAL     RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),
     &         DFDT( MAXULV ), UAVG( MAXULV ),
     &         UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),
     &         TRNMED( MAXUMU )

c+---------------------------------------------------------------------+

c                     ** Correct answers **

      INTEGER  MXPROB, MXCASE, MAXTAU, MAXMU, MAXAZ
      PARAMETER  ( MXPROB = 9, MXCASE = 8, MAXTAU = 5, MAXMU = 10,
     &             MAXAZ = 3 )
      REAL  TSTFIR( MAXTAU, MXCASE, MXPROB ),
     &      TSTFDN( MAXTAU, MXCASE, MXPROB ),
     &      TSTFUP( MAXTAU, MXCASE, MXPROB ),
     &      TSTDFD( MAXTAU, MXCASE, MXPROB ),
     &      TSTUU ( MAXTAU, MAXMU, MAXAZ, MXCASE, MXPROB )
      COMMON / DOCHEK / TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU

c+---------------------------------------------------------------------+

      INTEGER  MXTAU, MXMU, MXPHI
      PARAMETER     ( MXTAU = 5, MXMU = 10, MXPHI = 3 )
      CHARACTER*1   ABC(18)*1, TITLE*100, BLANKS*3
      LOGICAL       DOPROB( 13 )
      INTEGER       ICAS, IOD, ISS, IU, J, K, LC, LENTIT, LU, NPROB
      REAL          CMPFIR( MXTAU ), CMPFDN( MXTAU ), CMPFUP( MXTAU ),
     &              CMPDFD( MXTAU ), CMPUU ( MXTAU, MXMU, MXPHI ),
     &              PI

c     .. External Subroutines ..

      EXTERNAL  DISORT, ERRMSG, GETMOM, PRTFIN
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ASIN, FLOAT, INDEX
c     ..
      DATA  ABC / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
     &            'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r' /,
     &      PRNT / .TRUE., 3*.FALSE., .TRUE. /, ACCUR / 0.0 /,
     &      BLANKS / '   ' /,  DOPROB / 13*.TRUE. /


      PI = 2.* ASIN( 1.0 )


      IF( DOPROB(1) )  THEN

c **********************************************************************
c ****  Test Problem 1:  Isotropic Scattering                       ****
c ****  (Compare to Ref. VH1, Table 12)                             ****
c **********************************************************************

      NSTR = 16
      NLYR = 1
      NMOM = NSTR
      CALL  GETMOM( 1, 0.0, NMOM, PMOM )
      USRTAU    = .TRUE.
      NTAU      = 2
      UTAU( 1 ) = 0.0
      USRANG    = .TRUE.
      NUMU      = 6
      UMU( 1 )  = -1.0
      UMU( 2 )  = -0.5
      UMU( 3 )  = -0.1
      UMU( 4 )  =  0.1
      UMU( 5 )  =  0.5
      UMU( 6 )  =  1.0
      NPHI      = 1
      PHI( 1 )  = 0.0
      IBCND     = 0
      UMU0      = 0.1
      PHI0      = 0.0
      LAMBER    = .TRUE.
      ALBEDO    = 0.0
      PLANK     = .FALSE.
      ONLYFL    = .FALSE.

      DO 10  ICAS = 1, 6

         IF ( ICAS.EQ.1 ) THEN

            UTAU( 2 )  = 0.03125
            SSALB( 1 ) = 0.2
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.2 ) THEN

            UTAU( 2 )  = 0.03125
            SSALB( 1 ) = 1.0
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.3 ) THEN

            UTAU( 2 )  = 0.03125
            SSALB( 1 ) = 0.99
            FBEAM      = 0.0
            FISOT      = 1.0

         ELSE IF ( ICAS.EQ.4 ) THEN

            UTAU( 2 )  = 32.0
            SSALB( 1 ) = 0.2
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.5 ) THEN

            UTAU( 2 )  = 32.0
            SSALB( 1 ) = 1.0
            FBEAM      = PI / UMU0
            FISOT      = 0.0

         ELSE IF ( ICAS.EQ.6 ) THEN

            UTAU( 2 )  = 32.0
            SSALB( 1 ) = 0.99
            FBEAM      = 0.0
            FISOT      = 1.0

         END IF

         DTAUC( 1 ) = UTAU( 2 )

         WRITE( HEADER,'(3A,F9.5,A,F5.2)') 'Test Case No. 1',ABC(ICAS),
     &          ':  Isotropic Scattering, Ref. VH1, Table 12:  b =',
     &          UTAU(2), ', a =', SSALB(1)

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         NPROB = 1
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &       ERRMSG( 'Out of bounds in exact-answer arrays', .FALSE. )

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,
     &                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),
     &                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),
     &                 TSTUU(1,1,1,ICAS,NPROB), MAXTAU, MAXMU, MAXAZ )
  10  CONTINUE

      ENDIF


      IF( DOPROB(2) )  THEN

c **********************************************************************
c ****  Test Problem 2:  Rayleigh Scattering, Beam Source           ****
c ****  (Compare To Ref. SW, Table 1)                               ****
c **********************************************************************

      NSTR = 16
      NLYR = 1
      NMOM = NSTR
      CALL  GETMOM( 2, 0.0, NMOM, PMOM )
      USRTAU    = .TRUE.
      NTAU      = 2
      UTAU( 1 ) = 0.0
      USRANG    = .TRUE.
      NUMU      = 6
      UMU( 1 )  = -0.981986
      UMU( 2 )  = -0.538263
      UMU( 3 )  = -0.018014
      UMU( 4 )  =  0.018014
      UMU( 5 )  =  0.538263
      UMU( 6 )  =  0.981986
      NPHI      = 1
      PHI( 1 )  = 0.0
      IBCND     = 0
      FBEAM     = PI
      UMU0      = 0.080442
      PHI0      = 0.0
      FISOT     = 0.0
      LAMBER    = .TRUE.
      ALBEDO    = 0.0
      PLANK     = .FALSE.
      ONLYFL    = .FALSE.

      ICAS = 0
      DO 20  IOD = 1, 2

         IF ( IOD.EQ.1 )  UTAU( 2 ) = 0.2
         IF ( IOD.EQ.2 )  UTAU( 2 ) = 5.0

         DTAUC( 1 ) = UTAU( 2 )

         DO 20  ISS = 1, 2

            IF( ISS.EQ.1 )  SSALB( 1 ) = 0.5
            IF( ISS.EQ.2 )  SSALB( 1 ) = 1.0
            ICAS = ICAS + 1

            WRITE( HEADER, '(3A,F5.2,A,F9.6,A,F4.2)')
     &             'Test Case No. 2', ABC(ICAS),
     &             ', Rayleigh Scattering, Ref. SW, Table 1:  tau =',
     &             UTAU(2), ', mu0 =', UMU0, ', ss-albedo =', SSALB(1)

            CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                    WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                    USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                    UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                    TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                    HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                    MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                    ALBMED, TRNMED )

            NPROB = 2
            IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

            CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,
     &                    TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),
     &                    TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),
     &                    TSTUU(1,1,1,ICAS,NPROB), MAXTAU, MAXMU,
     &                    MAXAZ )
   20 CONTINUE

      ENDIF


      IF( DOPROB(3) )  THEN

c **********************************************************************
c ****  Test Problem 3:  Henyey-Greenstein Scattering               ****
c ****  (Compare To Ref. VH2, Table 37)                             ****
c **********************************************************************

      NSTR = 16
      NLYR = 1
      SSALB( 1 ) = 1.0
      NMOM = 32
      CALL  GETMOM( 3, 0.75, NMOM, PMOM )
      USRTAU    = .TRUE.
      NTAU      = 2
      UTAU( 1 ) = 0.0
      USRANG    = .TRUE.
      NUMU      = 6
      UMU( 1 )  = -1.0
      UMU( 2 )  = -0.5
      UMU( 3 )  = -0.1
      UMU( 4 )  =  0.1
      UMU( 5 )  =  0.5
      UMU( 6 )  =  1.0
      NPHI      = 1
      PHI( 1 )  = 0.0
      IBCND     = 0
      UMU0      = 1.0
      PHI0      = 0.0
      FBEAM      = PI / UMU0
      FISOT      = 0.0
      LAMBER    = .TRUE.
      ONLYFL    = .FALSE.
      ALBEDO    = 0.0
      PLANK     = .FALSE.

      DO 30  ICAS = 1, 2

         IF ( ICAS.EQ.1 ) THEN

            UTAU( 2 )  = 1.0

         ELSE IF ( ICAS.EQ.2 ) THEN

            UTAU( 2 )  = 8.0

         END IF

         DTAUC( 1 ) = UTAU( 2 )

         WRITE( HEADER, '(3A,F9.5,A,F5.2)') 'Test Case No. 3',ABC(ICAS),
     &         ', Henyey-Greenstein Scattering, Ref. VH2, Table 37,'
     &         //' g = 0.75, b =', UTAU(2), ', a =', SSALB(1)

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         NPROB = 3
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,
     &                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),
     &                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),
     &                 TSTUU(1,1,1,ICAS,NPROB), MAXTAU, MAXMU, MAXAZ )
30    CONTINUE

      ENDIF


      IF( DOPROB(4) )  THEN

c **********************************************************************
c ****  Test Problem 4:  Haze-L Scattering, Beam Source             ****
c ****  (Compare to Ref. GS, Tables 12-16)                          ****
c **********************************************************************

      NSTR = 32
      NLYR = 1
      NMOM = NSTR
      CALL  GETMOM( 4, 0.0, NMOM, PMOM )
      DTAUC( 1 ) = 1.0
      USRTAU     = .TRUE.
      NTAU       = 3
      UTAU( 1 )  = 0.0
      UTAU( 2 )  = 0.5
      UTAU( 3 )  = 1.0
      USRANG     = .TRUE.
      NUMU       = 6
      UMU( 1 )   = -1.0
      UMU( 2 )   = -0.5
      UMU( 3 )   = -0.1
      UMU( 4 )   =  0.1
      UMU( 5 )   =  0.5
      UMU( 6 )   =  1.0
      IBCND      = 0
      FBEAM      = PI
      PHI0       = 0.0
      FISOT      = 0.0
      LAMBER     = .TRUE.
      ALBEDO     = 0.0
      PLANK      = .FALSE.
      ONLYFL     = .FALSE.

      DO 40 ICAS = 1, 3

         WRITE( TITLE, '(3A)' ) 'Test Case No. 4', ABC(ICAS),
     &          ', Haze-L Scattering, Ref. GS, Table '
         LENTIT = INDEX( TITLE,BLANKS )

         IF ( ICAS.EQ.1 ) THEN

            SSALB( 1 ) = 1.0
            NPHI       = 1
            PHI( 1 )   = 0.0
            UMU0       = 1.0
            HEADER = TITLE(1:LENTIT) // ' 12'

         ELSE IF ( ICAS.EQ.2 ) THEN

            SSALB( 1 ) = 0.9
            NPHI       = 1
            PHI( 1 )   = 0.0
            UMU0       = 1.0
            HEADER = TITLE(1:LENTIT) // ' 13'

         ELSE IF ( ICAS.EQ.3 ) THEN

            SSALB( 1 ) = 0.9
            NPHI       = 3
            PHI( 1 )   = 0.0
            PHI( 2 )   = 90.0
            PHI( 3 )   = 180.0
            UMU0       = 0.5
            HEADER = TITLE(1:LENTIT) // ' 14-16'

         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         NPROB = 4
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,
     &                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),
     &                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),
     &                 TSTUU(1,1,1,ICAS,NPROB), MAXTAU, MAXMU, MAXAZ )
40    CONTINUE

      ENDIF


      IF( DOPROB(5) )  THEN

c **********************************************************************
c ****  Test Problem 5:  Cloud C.1 Scattering, Beam Source          ****
c ****  (Compare to Ref. GS, Tables 19-20)                          ****
c **********************************************************************

      NSTR = 48
      NLYR = 1
      NMOM = 299
      CALL  GETMOM( 5, 0.0, NMOM, PMOM )
      DTAUC( 1 ) = 64.0
      USRTAU     = .TRUE.
      NTAU       = 3
      USRANG     = .TRUE.
      NUMU       = 6
      UMU( 1 )   = -1.0
      UMU( 2 )   = -0.5
      UMU( 3 )   = -0.1
      UMU( 4 )   =  0.1
      UMU( 5 )   =  0.5
      UMU( 6 )   =  1.0
      NPHI       = 1
      PHI( 1 )   = 0.0
      IBCND      = 0
      FBEAM      = PI
      UMU0       = 1.0
      PHI0       = 0.0
      FISOT      = 0.0
      LAMBER     = .TRUE.
      ALBEDO     = 0.0
      PLANK      = .FALSE.
      ONLYFL     = .FALSE.

      DO 50 ICAS = 1, 2

         WRITE( TITLE, '(3A)' ) 'Test Case No. 5', ABC(ICAS),
     &          ', Cloud C.1 Scattering, Ref. GS, Table '
         LENTIT = INDEX( TITLE,BLANKS )

         IF ( ICAS.EQ.1 ) THEN

            UTAU( 1 )  = 0.0
            UTAU( 2 )  = 32.0
            UTAU( 3 )  = 64.0
            SSALB( 1 ) = 1.0
            HEADER = TITLE(1:LENTIT) // ' 19'

         END IF

         IF ( ICAS.EQ.2 ) THEN

            UTAU( 1 )  =  3.2
            UTAU( 2 )  = 12.8
            UTAU( 3 )  = 48.0
            SSALB( 1 ) = 0.9
            HEADER = TITLE(1:LENTIT) // ' 20'

         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         NPROB = 5
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,
     &                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),
     &                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),
     &                 TSTUU(1,1,1,ICAS,NPROB), MAXTAU, MAXMU, MAXAZ )
50    CONTINUE

      ENDIF


      IF( DOPROB(6) )  THEN

c **********************************************************************
c ****  Test Problem 6:  No Scattering, Increasingly Complex Sources****
c **********************************************************************

      NSTR       = 16
      NLYR       = 1
      SSALB( 1 ) = 0.0
      WVNMLO     = 0.0
      WVNMHI     = 50000.
      USRTAU     = .TRUE.
      USRANG     = .TRUE.
      NUMU       = 4
      UMU( 1 )   = -1.0
      UMU( 2 )   = -0.1
      UMU( 3 )   =  0.1
      UMU( 4 )   =  1.0
      NPHI       = 1
      PHI( 1 )   = 90.0
      IBCND      = 0
      FBEAM      = 200.0
      UMU0       = 0.5
      PHI0       = 0.0
      FISOT      = 0.0
      TEMIS      = 1.0
      ONLYFL     = .FALSE.

      DO 60  ICAS = 1, 8

         WRITE( TITLE, '(3A)' ) 'Test Case No. 6', ABC(ICAS),
     &          ': No Scattering; Source = Beam'
         LENTIT = INDEX( TITLE, BLANKS )

         IF ( ICAS.EQ.1 ) THEN

            NTAU = 2
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 0.0

         ELSE IF ( ICAS.GT.1 ) THEN

            NTAU = 3
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 0.5
            UTAU( 3 ) = 1.0

         END IF

         IF ( ICAS.EQ.1 ) THEN
c                                    ** Transparent medium, beam source
            DTAUC( 1 ) = 0.0
            LAMBER = .TRUE.
            ALBEDO = 0.0
            PLANK = .FALSE.
            HEADER = TITLE(1:LENTIT) // '; Bottom Albedo = 0'

         ELSE IF ( ICAS.EQ.2 ) THEN
c                                    ** Add some optical depth
            DTAUC( 1 ) = 1.0
            HEADER = TITLE(1:LENTIT) // '; Bottom Albedo = 0'

         ELSE IF ( ICAS.EQ.3 ) THEN
c                                   ** Add some isotropic reflection
            LAMBER = .TRUE.
            ALBEDO = 0.50
            PLANK = .FALSE.
            HEADER = TITLE(1:LENTIT) // '; Bottom Albedo=0.5 Lambert'

         ELSE IF ( ICAS.EQ.4 ) THEN
c                                   ** Use non-isotropic reflection
            DTAUC( 1 ) = 1.0
            LAMBER = .FALSE.
            PLANK = .FALSE.
            HEADER = TITLE(1:LENTIT) // '; Bottom Albedo = Non-Lambert'

         ELSE IF ( ICAS.EQ.5 ) THEN
c                                   ** Add some bottom-boundary emission
            DTAUC( 1 ) = 1.0
            TEMPER( 0 ) = 0.0
            TEMPER( 1 ) = 0.0
            LAMBER = .FALSE.
            BTEMP = 300.0
            TTEMP = 0.0
            PLANK = .TRUE.
            HEADER = TITLE(1:LENTIT) //
     &               ', Bottom Emission; Bott Alb = Non-Lambert'

         ELSE IF ( ICAS.EQ.6 ) THEN
c                                   ** Add some top-boundary diffuse
c                                      incidence (prescribed + emitted)
            DTAUC( 1 ) = 1.0
            TEMPER( 0 ) = 0.0
            TEMPER( 1 ) = 0.0
            FISOT  = 100.0 / PI
            LAMBER = .FALSE.
            BTEMP = 300.0
            TTEMP = 250.0
            PLANK = .TRUE.
            HEADER = TITLE(1:LENTIT) //
     &               ', Bottom+Top Emission; Bott Alb = Non-Lambert'

         ELSE IF ( ICAS.EQ.7 ) THEN
c                                   ** Add some internal emission
            DTAUC( 1 ) = 1.0
            TEMPER( 0 ) = 250.0
            TEMPER( 1 ) = 300.0
            LAMBER = .FALSE.
            BTEMP = 300.0
            TTEMP = 250.0
            PLANK = .TRUE.
            HEADER = TITLE(1:LENTIT) //
     &         ', Bottom+Top+Internal Emission; Bott Alb = Non-Lambert'

         ELSE IF ( ICAS.EQ.8 ) THEN
c                                   ** Increase the optical depth
            DTAUC( 1 ) = 10.0
            TEMPER( 0 ) = 250.0
            TEMPER( 1 ) = 300.0
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 1.0
            UTAU( 3 ) = 10.0
            LAMBER = .FALSE.
            BTEMP = 300.0
            TTEMP = 250.0
            PLANK = .TRUE.
            HEADER = TITLE(1:LENTIT) //
     &         ', Bottom+Top+Internal Emission; Bott Alb = Non-Lambert'

         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         NPROB = 6
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,
     &                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),
     &                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),
     &                 TSTUU(1,1,1,ICAS,NPROB), MAXTAU, MAXMU, MAXAZ )
60    CONTINUE

      ENDIF


      IF( DOPROB(7) )  THEN

c **********************************************************************
c ****  Test Problem 7:  Absorption + Scattering + All Possible     ****
c ****  Sources, Various Surface Reflectivities ( One Layer )       ****
c **** (Compare 7a,f Fluxes and Intensities to Ref. KS, Tables I-II ****
c **********************************************************************

      NLYR = 1
      USRTAU = .TRUE.
      USRANG = .TRUE.

      DO 70  ICAS = 1, 5

         WRITE( TITLE, '(2A)' ) 'Test Case No. 7', ABC(ICAS)
         LENTIT = INDEX( TITLE, BLANKS )

         IF ( ICAS.EQ.1 ) THEN

            NSTR = 16
            NMOM = 16
            DTAUC( 1 ) = 1.0
            SSALB( 1 ) = 0.1
            CALL  GETMOM( 3, 0.05, NMOM, PMOM )
            TEMPER(0) = 200.0
            TEMPER(1) = 300.0
            WVNMLO = 300.0
            WVNMHI = 800.0
            NTAU = 2
            UTAU(1) = 0.0
            UTAU(2) = 1.0
            NUMU   = 2
            UMU(1) = -1.
            UMU(2) =  1.
            NPHI = 1
            PHI(1) = 0.0
            IBCND = 0
            FBEAM  = 0.0
            FISOT = 0.0
            LAMBER = .TRUE.
            ALBEDO = 0.0
            PLANK = .TRUE.
            BTEMP  = 0.0
            TTEMP  = 0.0
            TEMIS = 1.0
            ONLYFL = .TRUE.
            HEADER = TITLE(1:LENTIT) // ': Absorption + Scattering, '//
     &        'Internal Thermal Sources; Ref. KS, Table I, '//
     &        'tau = 1.0, a = 0.1, g = 0.05'

         ELSE IF ( ICAS.EQ.2 ) THEN

            DTAUC( 1 ) = 100.0
            SSALB( 1 ) = 0.95
            CALL  GETMOM( 3, 0.75, NMOM, PMOM )
            TEMPER(0) = 200.0
            TEMPER(1) = 300.0
            WVNMLO = 2702.99
            WVNMHI = 2703.01
            NTAU = 2
            UTAU(1) = 0.0
            UTAU(2) = 100.0
            USRANG = .TRUE.
            NUMU   = 2
            UMU(1) = -1.
            UMU(2) =  1.
            NPHI = 1
            PHI(1) = 0.0
            IBCND = 0
            FBEAM  = 0.0
            FISOT = 0.0
            LAMBER = .TRUE.
            ALBEDO = 0.0
            PLANK = .TRUE.
            BTEMP  = 0.0
            TTEMP  = 0.0
            TEMIS = 1.0
            ONLYFL = .FALSE.
            HEADER = TITLE(1:LENTIT) // ': Absorption + Scattering, '//
     &        'Internal Thermal Sources; Ref. KS, Table II, '//
     &        'tau = 100.0, a = 0.95, g = 0.75'

         ELSE IF ( ICAS.EQ.3 ) THEN

            NSTR = 12
            NMOM = 12
            TEMPER( 0 ) = 300.0
            TEMPER( 1 ) = 200.0
            CALL  GETMOM( 3, 0.8, NMOM, PMOM )
            DTAUC( 1 )  = 1.0
            SSALB( 1 )  = 0.5
            WVNMLO      = 0.0
            WVNMHI      = 50000.
            NTAU        = 3
            UTAU( 1 )   = 0.0
            UTAU( 2 )   = 0.5
            UTAU( 3 )   = 1.0
            USRANG      = .TRUE.
            NUMU        = 4
            UMU( 1 )    = -1.0
            UMU( 2 )    = -0.1
            UMU( 3 )    =  0.1
            UMU( 4 )    =  1.0
            NPHI        = 2
            PHI( 1 )    = 0.0
            PHI( 2 )    = 90.0
            IBCND       = 0
            FBEAM       = 200.0
            UMU0        = 0.5
            PHI0        = 0.0
            FISOT       = 100.0
            BTEMP       = 320.0
            TTEMP       = 100.0
            TEMIS       = 1.0
            PLANK       = .TRUE.
            ONLYFL      = .FALSE.
            LAMBER = .TRUE.
            ALBEDO = 0.0
            HEADER = TITLE(1:LENTIT) // ': Absorption + '//
     &        'Henyey-Greenstein Scattering, All Sources, '//
     &        'Bottom Albedo = 0'

         ELSE IF ( ICAS.EQ.4 ) THEN

            LAMBER = .TRUE.
            ALBEDO = 1.0
            HEADER = TITLE(1:LENTIT) // ': Absorption + '//
     &        'Henyey-Greenstein Scattering, All Sources, '//
     &        'Bottom Albedo = 1'

         ELSE IF ( ICAS.EQ.5 ) THEN

            LAMBER = .FALSE.
            HEADER = TITLE(1:LENTIT) // ': Absorption + '//
     &        'Henyey-Greenstein Scattering, All Sources, '//
     &        'Bottom Albedo = BDR Function'

         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         NPROB = 7
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,
     &                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),
     &                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),
     &                 TSTUU(1,1,1,ICAS,NPROB), MAXTAU, MAXMU, MAXAZ )
70    CONTINUE

      ENDIF


      IF( DOPROB(8) )  THEN

c **********************************************************************
c ****  Test Problem 8:  Absorbing/Isotropic-Scattering Medium      ****
c ****  With Two Computational Layers                               ****
c **** (Compare Fluxes To Ref. OS, Table 1)                         ****
c **********************************************************************

      NSTR = 8
      NLYR = 2
      NMOM = NSTR
      CALL  GETMOM( 1, 0.0, NMOM, PMOM(0,1) )
      CALL  GETMOM( 1, 0.0, NMOM, PMOM(0,2) )
      USRTAU   = .TRUE.
      USRANG   = .TRUE.
      NUMU     = 4
      UMU( 1 ) = -1.0
      UMU( 2 ) = -0.2
      UMU( 3 ) =  0.2
      UMU( 4 ) =  1.0
      NPHI     = 1
      PHI( 1 ) = 60.0
      IBCND    = 0
      FBEAM    = 0.0
      FISOT    = 1.0 / PI
      LAMBER   = .TRUE.
      ALBEDO   = 0.0
      PLANK    = .FALSE.
      ONLYFL   = .FALSE.

      DO 80  ICAS = 1, 3

         IF ( ICAS.EQ.1 ) THEN

            DTAUC( 1 ) = 0.25
            DTAUC( 2 ) = 0.25
            SSALB( 1 ) = 0.5
            SSALB( 2 ) = 0.3
            NTAU       = 3
            UTAU( 1 )  = 0.0
            UTAU( 2 )  = 0.25
            UTAU( 3 )  = 0.5
            HEADER = 'Test Case No. 8a:  Ref. OS, Table 1,'
     &               // ' Line 4 (Two Inhomogeneous Layers)'

         ELSE IF ( ICAS.EQ.2 ) THEN

            DTAUC( 1 ) = 0.25
            DTAUC( 2 ) = 0.25
            SSALB( 1 ) = 0.8
            SSALB( 2 ) = 0.95
            NTAU       = 3
            UTAU( 1 )  = 0.0
            UTAU( 2 )  = 0.25
            UTAU( 3 )  = 0.5
            HEADER = 'Test Case No. 8b:  Ref. OS, Table 1,'
     &               // ' Line 1 (Two Inhomogeneous Layers)'

         ELSE IF ( ICAS.EQ.3 ) THEN

            DTAUC( 1 ) = 1.0
            DTAUC( 2 ) = 2.0
            SSALB( 1 ) = 0.8
            SSALB( 2 ) = 0.95
            NTAU       = 3
            UTAU( 1 )  = 0.0
            UTAU( 2 )  = 1.0
            UTAU( 3 )  = 3.0
            HEADER = 'Test Case No. 8c:  Ref. OS, Table 1,'
     &               // ' Line 13 (Two Inhomogeneous Layers)'
         ENDIF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         NPROB = 8
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,
     &                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),
     &                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),
     &                 TSTUU(1,1,1,ICAS,NPROB), MAXTAU, MAXMU, MAXAZ )
  80  CONTINUE

      ENDIF


      IF( DOPROB(9) )  THEN

c **********************************************************************
c ****  Test Problem 9:  General Emitting/Absorbing/Scattering      ****
c ****  Medium with Every Computational Layer Different.            ****
c **** (Compare 9a,b Fluxes to Ref. DGIS, Tables VI-VII, beta = 0)  ****
c **********************************************************************

      NSTR = 8
      NLYR = 6
      NMOM = 8
      DO 86  LC = 1, NLYR
         DTAUC( LC ) = LC
         SSALB( LC ) = 0.6 + LC*0.05
 86   CONTINUE
      USRTAU    = .TRUE.
      NTAU      = 5
      UTAU( 1 ) = 0.0
      UTAU( 2 ) = 1.05
      UTAU( 3 ) = 2.1
      UTAU( 4 ) = 6.0
      UTAU( 5 ) = 21.0
      USRANG    = .TRUE.
      NUMU      = 4
      UMU( 1 )  = -1.0
      UMU( 2 )  = -0.2
      UMU( 3 )  =  0.2
      UMU( 4 )  =  1.0
      NPHI      = 1
      PHI( 1 )  = 60.0
      IBCND     = 0
      FBEAM     = 0.0
      FISOT     = 1.0 / PI
      LAMBER    = .TRUE.
      ONLYFL    = .FALSE.

      DO 90  ICAS = 1, 3

         IF ( ICAS.EQ.1 ) THEN

            DO 87  LC = 1, NLYR
               CALL  GETMOM( 1, 0.0, NMOM, PMOM(0,LC) )
 87         CONTINUE
            ALBEDO = 0.0
            PLANK  = .FALSE.
            HEADER = 'Test Case No. 9a:  Ref. DGIS, Tables VI-VII,'
     &               // ' beta=l=0 (multiple inhomogeneous layers)'

         ELSE IF ( ICAS.EQ.2 ) THEN

            PMOM(0,1) = 1.0
            PMOM(1,1) = 2.00916/3.
            PMOM(2,1) = 1.56339/5.
            PMOM(3,1) = 0.67407/7.
            PMOM(4,1) = 0.22215/9.
            PMOM(5,1) = 0.04725/11.
            PMOM(6,1) = 0.00671/13.
            PMOM(7,1) = 0.00068/15.
            PMOM(8,1) = 0.00005/17.
            DO 88  LC = 2, NLYR
               DO 88  K = 0, 8
                  PMOM(K,LC) = PMOM(K,1)
 88         CONTINUE
            HEADER = 'Test Case No. 9b:  Ref. DGIS, Tables VI-VII,'
     &               // ' beta=0,l=8 (multiple inhomogeneous layers)'

         ELSE IF ( ICAS.EQ.3 ) THEN

            TEMPER( 0 ) = 600.0
            DO 89  LC = 1, NLYR
               CALL  GETMOM( 3, FLOAT(LC)/7.0, NMOM, PMOM(0,LC) )
               TEMPER( LC ) = 600.0 + LC*10.0
 89         CONTINUE
            NPHI = 3
            PHI( 1 ) = 60.0
            PHI( 2 ) = 120.0
            PHI( 3 ) = 180.0
            FBEAM    = PI
            UMU0     = 0.5
            PHI0     = 0.0
            FISOT    = 1.0
            ALBEDO   =  0.5
            PLANK    = .TRUE.
            WVNMLO   =  999.0
            WVNMHI   = 1000.0
            BTEMP    = 700.0
            TTEMP    = 550.0
            TEMIS    = 1.0
            HEADER = 'Test Case No. 9c:  Generalization of 9A '//
     &               'to include all possible complexity'

         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         NPROB = 9
         IF( NPROB.GT.MXPROB .OR. ICAS.GT.MXCASE )  CALL
     &         ERRMSG( 'Out of bounds in exact-answer arrays',.FALSE.)

         CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                 MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT, UU,
     &                 TSTFIR(1,ICAS,NPROB), TSTFDN(1,ICAS,NPROB),
     &                 TSTFUP(1,ICAS,NPROB), TSTDFD(1,ICAS,NPROB),
     &                 TSTUU(1,1,1,ICAS,NPROB), MAXTAU, MAXMU, MAXAZ )
  90  CONTINUE

      ENDIF


      IF( DOPROB(10) )  THEN

c **********************************************************************
c ****  Test Problem 10: Compare USRANG = True With USRANG = False  ****
c ****  take Problem 9c (our most general case) but only 4 Streams  ****
c **********************************************************************

      NSTR = 4
      NLYR = 6
      NMOM = NSTR
      TEMPER( 0 ) = 600.0
      DO 97  LC = 1, NLYR
         DTAUC( LC ) = LC
         SSALB( LC ) = 0.6 + LC*0.05
         CALL  GETMOM( 3, FLOAT(LC)/(NLYR+1), NMOM, PMOM(0,LC) )
         TEMPER( LC ) = 600.0 + LC*10.0
  97  CONTINUE
      USRTAU    = .TRUE.
      NTAU      = 3
      UTAU( 1 ) = 0.0
      UTAU( 2 ) = 2.1
      UTAU( 3 ) = 21.0
      NPHI      = 2
      PHI( 1 )  = 60.0
      PHI( 2 )  = 120.0
      IBCND     = 0
      FBEAM     = PI
      UMU0      = 0.5
      PHI0      = 0.0
      FISOT     = 1.0
      LAMBER    = .TRUE.
      ALBEDO    =  0.5
      PLANK     = .TRUE.
      WVNMLO    =  999.0
      WVNMHI    = 1000.0
      BTEMP     = 700.0
      TTEMP     = 550.0
      TEMIS     = 1.0
      ONLYFL    = .FALSE.

      DO 100  ICAS = 1, 2

         IF ( ICAS.EQ.1 ) THEN

            USRANG = .TRUE.
            NUMU      = 4
            UMU( 1 )  = - 0.788675129
            UMU( 2 )  = - 0.211324871
            UMU( 3 )  =   0.211324871
            UMU( 4 )  =   0.788675129
            PRNT( 2 ) = .TRUE.
            PRNT( 3 ) = .TRUE.
            HEADER = 'Test Case No. 10a:  like 9c, USRANG = True'

         ELSE IF ( ICAS.EQ.2 ) THEN

            USRANG    = .FALSE.
            NUMU      = 0
            PRNT( 2 ) = .FALSE.
            PRNT( 3 ) = .FALSE.
            HEADER = 'Test Case No. 10b:  like 9C, USRANG = False'

         ENDIF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         IF ( ICAS.EQ.1 ) THEN
c                               ** Save results to compare to case 2
            DO 98  LU = 1, NTAU
               CMPFIR( LU ) = RFLDIR( LU )
               CMPFDN( LU ) = RFLDN( LU )
               CMPFUP( LU ) = FLUP( LU )
               CMPDFD( LU ) = DFDT( LU )
               DO 98  IU = 1, NUMU
                  DO 98  J = 1, NPHI
                     CMPUU( LU, IU, J ) = UU( IU, LU, J )
   98       CONTINUE

         ELSE IF ( ICAS.EQ.2 ) THEN

            CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT,
     &                    UU, CMPFIR, CMPFDN, CMPFUP, CMPDFD, CMPUU,
     &                    MXTAU, MXMU, MXPHI )
         END IF

 100  CONTINUE

      ENDIF


      IF( DOPROB(11) )  THEN

c **********************************************************************
c ****  Test Problem 11: Single-Layer vs. Multiple Layers           ****
c ****  11a: Results at user levels for one computational layer     ****
c ****  11b: Single layer of 11a subdivided into multiple           ****
c ****       computational layers at the 11a user levels            ****
c **********************************************************************

      NSTR      = 16
      NMOM      = NSTR
      USRANG    = .TRUE.
      NUMU      = 4
      UMU( 1 )  = -1.0
      UMU( 2 )  = -0.1
      UMU( 3 )  =  0.1
      UMU( 4 )  =  1.0
      NPHI      = 2
      PHI( 1 )  = 0.0
      PHI( 2 )  = 90.0
      IBCND     = 0
      FBEAM     = 1.0
      UMU0      = 0.5
      PHI0      = 0.0
      FISOT     = 0.5 / PI
      LAMBER    = .TRUE.
      ALBEDO    = 0.5
      PLANK     = .FALSE.
      ONLYFL    = .FALSE.

      DO 110  ICAS = 1, 2

         IF ( ICAS.EQ.1 ) THEN

            NLYR = 1
            DTAUC( 1 ) = 1.0
            SSALB( 1 ) = 0.9
            CALL  GETMOM( 1, 0.0, NMOM, PMOM )
            USRTAU    = .TRUE.
            NTAU      = 4
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 0.05
            UTAU( 3 ) = 0.5
            UTAU( 4 ) = 1.0
            PRNT( 2 ) = .TRUE.
            PRNT( 3 ) = .TRUE.
            HEADER = 'Test Case No. 11a: One Isotropic-Scattering Layer'

         ELSE IF ( ICAS.EQ.2 ) THEN

            NLYR = NTAU - 1
            DO 107 LC = 1, NLYR
               DTAUC( LC ) = UTAU(LC+1) - UTAU(LC)
               SSALB( LC ) = 0.9
               CALL  GETMOM( 1, 0.0, NMOM, PMOM(0,LC) )
107         CONTINUE
            USRTAU    = .FALSE.
            PRNT( 2 ) = .FALSE.
            PRNT( 3 ) = .FALSE.
            HEADER = 'Test Case No. 11b: Same as 11a but treated as' //
     &               ' multiple layers'
         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         IF ( ICAS.EQ.1 ) THEN
c                               ** Save results to compare to case 2
            DO 108  LU = 1, NTAU
               CMPFIR( LU ) = RFLDIR( LU )
               CMPFDN( LU ) = RFLDN( LU )
               CMPFUP( LU ) = FLUP( LU )
               CMPDFD( LU ) = DFDT( LU )
               DO 108  IU = 1, NUMU
                  DO 108  J = 1, NPHI
                     CMPUU( LU, IU, J ) = UU( IU, LU, J )
  108       CONTINUE

         ELSE IF ( ICAS.EQ.2 ) THEN

            CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, RFLDIR, RFLDN,FLUP, DFDT,
     &                    UU, CMPFIR, CMPFDN, CMPFUP, CMPDFD, CMPUU,
     &                    MXTAU, MXMU, MXPHI )
         END IF

110   CONTINUE

      ENDIF


      IF( DOPROB(12) )  THEN

c **********************************************************************
c ****  Test Problem 12: Test Absorption-Optical-Depth Shortcut     ****
c ****  compares cases where the DISORT shortcut for absorption     ****
c ****  optical depth .GT. 10 is not used (12a), then is used (12b) ****
c ****  (this shortcut is only employed when  PLANK = False.)       ****
c **********************************************************************

      NSTR    = 20
      NMOM    = NSTR
      USRANG  = .TRUE.
      NUMU     = 4
      UMU( 1 ) = -1.0
      UMU( 2 ) = -0.1
      UMU( 3 ) =  0.1
      UMU( 4 ) =  1.0
      NPHI     = 1
      PHI( 1 ) = 0.0
      IBCND    = 0
      FBEAM    = 1.0
      UMU0     = 1.0
      PHI0     = 0.0
      FISOT    = 0.0
      LAMBER   = .TRUE.
      ALBEDO   = 1.0
      PLANK    = .FALSE.
      ONLYFL   = .FALSE.

      DO 120  ICAS = 1, 2

         IF ( ICAS.EQ.1 ) THEN

            NLYR = 1
            DTAUC( 1 ) = 20.1
            SSALB( 1 ) = 0.5
            CALL  GETMOM( 3, 0.9, NMOM, PMOM )
            USRTAU    = .TRUE.
            NTAU      = 4
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 10.0
            UTAU( 3 ) = 19.9
            UTAU( 4 ) = 20.1
            PRNT( 2 ) = .TRUE.
            PRNT( 3 ) = .TRUE.
            HEADER = 'Test Case No. 12a:  Overhead Beam Striking '//
     &               'Absorbing/Scattering Medium'

         ELSE IF ( ICAS.EQ.2 ) THEN

            NLYR = NTAU - 1
            DO 117 LC = 1, NLYR
               DTAUC( LC ) = UTAU(LC+1) - UTAU(LC)
               SSALB( LC ) = 0.5
               CALL  GETMOM( 3, 0.9, NMOM, PMOM(0,LC) )
117         CONTINUE
            USRTAU    = .FALSE.
            PRNT( 2 ) = .FALSE.
            PRNT( 3 ) = .FALSE.
            HEADER = 'Test Case No. 12b: Same as 12a but uses shortcut'
     &               // ' for absorption optical depth .GT. 10'
         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

         IF ( ICAS.EQ.1 ) THEN
c                               ** Save results to compare to case 2
            DO 118  LU = 1, NTAU
               CMPFIR( LU ) = RFLDIR( LU )
               CMPFDN( LU ) = RFLDN( LU )
               CMPFUP( LU ) = FLUP( LU )
               CMPDFD( LU ) = DFDT( LU )
               DO 118  IU = 1, NUMU
                  DO 118  J = 1, NPHI
                     CMPUU( LU, IU, J ) = UU( IU, LU, J )
  118       CONTINUE

         ELSE IF ( ICAS.EQ.2 ) THEN

            CALL  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT,
     &                    UU, CMPFIR, CMPFDN, CMPFUP, CMPDFD, CMPUU,
     &                    MXTAU, MXMU, MXPHI )
         END IF

120   CONTINUE

      ENDIF


      IF( DOPROB(13) )  THEN

c **********************************************************************
c ****  Test Problem 13: Test shortcut for flux albedo, transmission ***
c **** ( shortcut gives flux albedo, transmission of entire medium  ****
c ****   as a function of sun angle )                               ****
c ****  13a,c = Shortcut;  13b,d = Brute Force Method               ****
c **********************************************************************

      NSTR   = 16
      NMOM   = NSTR
      NPHI   = 0
      PHI0   = 0.0
      ALBEDO = 0.5

      DO 130  ICAS = 1, 4

         IF ( ICAS.EQ.1 ) THEN

            IBCND      = 1
            NLYR       = 1
            DTAUC( 1 ) = 1.0
            SSALB( 1 ) = 0.99
            CALL  GETMOM( 3, 0.8, NMOM, PMOM )
            PRNT( 4 )  = .TRUE.
            PRNT( 2 )  = .FALSE.
            USRANG     = .TRUE.
            NUMU       = 1
            UMU( 1 )   =  0.5
            HEADER = 'Test Case No. 13a:  Albedo and Transmissivity'//
     &               ' from Shortcut, Single Layer'

         ELSE IF ( ICAS.EQ.2 ) THEN

            IBCND     = 0
            USRTAU    = .TRUE.
            NTAU      = 2
            UTAU( 1 ) = 0.0
            UTAU( 2 ) = 1.0
            UMU0      = 0.5
            FBEAM     = 1.0 / UMU0
            FISOT     = 0.0
            LAMBER    = .TRUE.
            PLANK     = .FALSE.
            ONLYFL    = .TRUE.
            PRNT( 4 ) = .FALSE.
            PRNT( 2 ) = .TRUE.
            HEADER = 'Test Case No. 13b:  Albedo and Transmissivity'//
     &               ' by Regular Method, Single Layer'

         ELSE IF ( ICAS.EQ.3 ) THEN

            IBCND     = 1
            PRNT( 4 ) = .TRUE.
            PRNT( 2 ) = .FALSE.
            NLYR      = 2
            DO 125  LC = 1, NLYR
               DTAUC( LC ) = 1.0 / NLYR
               CALL  GETMOM( 3, 0.8, NMOM, PMOM(0,LC) )
  125       CONTINUE
            SSALB( 1 ) = 0.99
            SSALB( 2 ) = 0.50
            HEADER = 'Test Case No. 13c:  Albedo and Transmissivity'//
     &               ' from Shortcut, Multiple Layer'

         ELSE IF ( ICAS.EQ.4 ) THEN

            IBCND = 0
            PRNT( 4 ) = .FALSE.
            PRNT( 2 ) = .TRUE.
            HEADER = 'Test Case No. 13d:  Albedo and Transmissivity'//
     &               ' by Regular Method, Multiple Layer'
         END IF

         CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )
130   CONTINUE

      ENDIF

      STOP
      END

      REAL FUNCTION  BDREF( WVNMLO, WVNMHI, MU, MUP, DPHI )

c      Supplies surface bi-directional reflectivity.
c
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c              angles) are positive.
c
c  INPUT:
c
c    WVNMLO : Lower wavenumber (inv cm) of spectral interval
c
c    WVNMHI : Upper wavenumber (inv cm) of spectral interval
c
c    MU     : Cosine of angle of reflection (positive)
c
c    MUP    : Cosine of angle of incidence (positive)
c
c    DPHI   : Difference of azimuth angles of incidence and reflection
c                (radians)
c
c  LOCAL VARIABLES:
c
c    IREF   : bidirectional reflectance options
c             1 - Hapke's BDR model
c
c    B0     : empirical factor to account for the finite size of
c             particles in Hapke's BDR model
c
c    B      : term that accounts for the opposition effect
c             (retroreflectance, hot spot) in Hapke's BDR model
c
c    CTHETA : cosine of phase angle in Hapke's BDR model
c
c    GAMMA  : albedo factor in Hapke's BDR model
c
c    H0     : H( mu0 ) in Hapke's BDR model
c
c    H      : H( mu ) in Hapke's BDR model
c
c    HH     : angular width parameter of opposition effect in Hapke's
c             BDR model
c
c    P      : scattering phase function in Hapke's BDR model
c
c    THETA  : phase angle (radians); the angle between incidence and
c             reflection directions in Hapke's BDR model
c
c    W      : single scattering albedo in Hapke's BDR model
c
c
c   Called by- DREF, SURFAC
c +-------------------------------------------------------------------+
c     .. Scalar Arguments ..

      REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
c     ..
c     .. Local Scalars ..

      INTEGER   IREF
      REAL      B0, B, CTHETA, GAMMA, H0, H, HH, P, THETA, W
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC COS, SQRT
c     ..

      IREF = 1

      IF ( IREF.EQ.1 ) THEN

c                              ** Hapke's BRDF model (times Pi/Mu0)
c                              ** (Hapke, B., Theory of reflectance
c                              ** and emittance spectroscopy, Cambridge
c                              ** University Press, 1993, Eq. 8.89 on
c                              ** page 233. Parameters are from
c                              ** Fig. 8.15 on page 231, expect for w.)

         CTHETA = MU * MUP + (1.-MU**2)**.5 * (1.-MUP**2)**.5
     &            * COS( DPHI )
         THETA = ACOS( CTHETA )

         P    = 1. + 0.5 * CTHETA

         HH   = 0.06
         B0   = 1.0
         B    = B0 * HH / ( HH + TAN( THETA/2.) )

         W = 0.6
         GAMMA = SQRT( 1. - W )
         H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
         H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )

         BDREF = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

      ELSE

         CALL ERRMSG( 'BDREF--Need to supply surface BDRF model',
     &                 .TRUE.)

      END IF

      RETURN
      END

      SUBROUTINE  GETMOM( IPHAS, GG, NMOM, PMOM )

c        Calculate phase function Legendre expansion coefficients
c        in various special cases


c       INPUT: IPHAS   Phase function options
c                      1 : Isotropic
c                      2 : Rayleigh
c                      3 : Henyey-Greenstein with asymmetry factor GG
c                      4 : Haze L as specified by Garcia/Siewert
c                      5 : Cloud C.1 as specified by Garcia/Siewert

c              GG      Asymmetry factor for Henyey-Greenstein case

c              NMOM    Index of highest Legendre coefficient needed
c                        ( = number of streams 'NSTR'  chosen
c                         for the discrete ordinate method)

c      OUTPUT: PMOM(K)  Legendre expansion coefficients (K=0 to NMOM)
c                         (be sure to dimension '0:maxval' in calling
c                          program)

c      Reference:  Garcia, R. and C. Siewert, 1985: Benchmark Results
c                     in Radiative Transfer, Transp. Theory and Stat.
c                     Physics 14, 437-484, Tables 10 And 17
c ------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   IPHAS, NMOM
      REAL      GG
c     ..
c     .. Array Arguments ..

      REAL      PMOM( 0:NMOM )
c     ..
c     .. Local Scalars ..

      INTEGER   K
c     ..
c     .. Local Arrays ..

      REAL      CLDMOM( 299 ), HAZELM( 82 )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..

      DATA HAZELM /  2.41260, 3.23047, 3.37296, 3.23150, 2.89350,
     A               2.49594, 2.11361, 1.74812, 1.44692, 1.17714,
     B               0.96643, 0.78237, 0.64114, 0.51966, 0.42563,
     C               0.34688, 0.28351, 0.23317, 0.18963, 0.15788,
     D               0.12739, 0.10762, 0.08597, 0.07381, 0.05828,
     E               0.05089, 0.03971, 0.03524, 0.02720, 0.02451,
     F               0.01874, 0.01711, 0.01298, 0.01198, 0.00904,
     G               0.00841, 0.00634, 0.00592, 0.00446, 0.00418,
     H               0.00316, 0.00296, 0.00225, 0.00210, 0.00160,
     I               0.00150, 0.00115, 0.00107, 0.00082, 0.00077,
     J               0.00059, 0.00055, 0.00043, 0.00040, 0.00031,
     K               0.00029, 0.00023, 0.00021, 0.00017, 0.00015,
     L               0.00012, 0.00011, 0.00009, 0.00008, 0.00006,
     M               0.00006, 0.00005, 0.00004, 0.00004, 0.00003,
     N               0.00003, 3*0.00002, 8*0.00001 /

      DATA  ( CLDMOM(K), K = 1, 159 ) /
     A  2.544,  3.883,  4.568,  5.235,  5.887,  6.457,  7.177,  7.859,
     B  8.494,  9.286,  9.856, 10.615, 11.229, 11.851, 12.503, 13.058,
     C 13.626, 14.209, 14.660, 15.231, 15.641, 16.126, 16.539, 16.934,
     D 17.325, 17.673, 17.999, 18.329, 18.588, 18.885, 19.103, 19.345,
     E 19.537, 19.721, 19.884, 20.024, 20.145, 20.251, 20.330, 20.401,
     F 20.444, 20.477, 20.489, 20.483, 20.467, 20.427, 20.382, 20.310,
     G 20.236, 20.136, 20.036, 19.909, 19.785, 19.632, 19.486, 19.311,
     H 19.145, 18.949, 18.764, 18.551, 18.348, 18.119, 17.901, 17.659,
     I 17.428, 17.174, 16.931, 16.668, 16.415, 16.144, 15.883, 15.606,
     J 15.338, 15.058, 14.784, 14.501, 14.225, 13.941, 13.662, 13.378,
     K 13.098, 12.816, 12.536, 12.257, 11.978, 11.703, 11.427, 11.156,
     L 10.884, 10.618, 10.350, 10.090,  9.827,  9.574,  9.318,  9.072,
     M  8.822, 8.584, 8.340, 8.110, 7.874, 7.652, 7.424, 7.211, 6.990,
     N  6.785, 6.573, 6.377, 6.173, 5.986, 5.790, 5.612, 5.424, 5.255,
     O  5.075, 4.915, 4.744, 4.592, 4.429, 4.285, 4.130, 3.994, 3.847,
     P  3.719, 3.580, 3.459, 3.327, 3.214, 3.090, 2.983, 2.866, 2.766,
     Q  2.656, 2.562, 2.459, 2.372, 2.274, 2.193, 2.102, 2.025, 1.940,
     R  1.869, 1.790, 1.723, 1.649, 1.588, 1.518, 1.461, 1.397, 1.344,
     S  1.284, 1.235, 1.179, 1.134, 1.082, 1.040, 0.992, 0.954, 0.909 /
      DATA  ( CLDMOM(K), K = 160, 299 ) /
     T  0.873, 0.832, 0.799, 0.762, 0.731, 0.696, 0.668, 0.636, 0.610,
     U  0.581, 0.557, 0.530, 0.508, 0.483, 0.463, 0.440, 0.422, 0.401,
     V  0.384, 0.364, 0.349, 0.331, 0.317, 0.301, 0.288, 0.273, 0.262,
     W  0.248, 0.238, 0.225, 0.215, 0.204, 0.195, 0.185, 0.177, 0.167,
     X  0.160, 0.151, 0.145, 0.137, 0.131, 0.124, 0.118, 0.112, 0.107,
     Y  0.101, 0.097, 0.091, 0.087, 0.082, 0.079, 0.074, 0.071, 0.067,
     Z  0.064, 0.060, 0.057, 0.054, 0.052, 0.049, 0.047, 0.044, 0.042,
     A  0.039, 0.038, 0.035, 0.034, 0.032, 0.030, 0.029, 0.027, 0.026,
     B  0.024, 0.023, 0.022, 0.021, 0.020, 0.018, 0.018, 0.017, 0.016,
     C  0.015, 0.014, 0.013, 0.013, 0.012, 0.011, 0.011, 0.010, 0.009,
     D  0.009, 3*0.008, 2*0.007, 3*0.006, 4*0.005, 4*0.004, 6*0.003,
     E  9*0.002, 18*0.001 /


      IF ( IPHAS.LT.1 .OR. IPHAS.GT.5 )
     &     CALL ERRMSG( 'GETMOM--bad input variable IPHAS',.TRUE.)

      IF ( IPHAS.EQ.3 .AND. (GG.LE.-1.0 .OR. GG.GE.1.0) )
     &     CALL ERRMSG( 'GETMOM--bad input variable GG',.TRUE.)

      IF ( NMOM.LT.2 )
     &     CALL ERRMSG( 'GETMOM--bad input variable NMOM',.TRUE.)


      PMOM(0) = 1.0
      DO  10  K = 1, NMOM
         PMOM(K) = 0.0
   10 CONTINUE


      IF ( IPHAS.EQ.2 )  THEN
c                                       ** Rayleigh phase function
         PMOM(2) = 0.1

      ELSE IF ( IPHAS.EQ.3 ) THEN
c                                       ** Henyey-Greenstein phase fcn
         DO  20  K = 1, NMOM
            PMOM(K) = GG**K
   20    CONTINUE

      ELSE IF ( IPHAS.EQ.4 ) THEN
c                                        ** Haze-L phase function
         DO  30  K = 1, MIN(82,NMOM)
            PMOM(K) = HAZELM(K) / ( 2*K+1 )
   30    CONTINUE

      ELSE IF ( IPHAS.EQ.5 ) THEN
c                                        ** Cloud C.1 phase function
         DO  40  K = 1, MIN(298,NMOM)
            PMOM(K) = CLDMOM(K) / ( 2*K+1 )
40       CONTINUE

      END IF

      END

      SUBROUTINE  PRTFIN( UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                    MAXUMU, ONLYFL, RFLDIR, RFLDN, FLUP, DFDT,
     &                    UU, TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU,
     &                    MAXTAU, MAXMU, MAXAZ )

c        Print DISORT results and, directly beneath them, their
c        ratios to the correct answers;  print number of non-unit
c        ratios that occur but try to count just the cases where
c        there is a real disagreement and not those where flux or
c        intensity are down at their noise level (defined as 10^(-6)
c        times their maximum value).  d(flux)/d(tau) is treated the
c        same as fluxes in this noise estimation even though it
c        is a different type of quantity (although with flux units).

c     INPUT :   TSTFIR  correct direct flux
c               TSTFDN  correct diffuse down flux
c               TSTFUP  correct diffuse up flux
c               TSTDFD  correct d(flux)/d(optical depth)
c               TSTUU   correct intensity
c               (remaining input = DISORT I/O variables)

c --------------------------------------------------------------------

c     .. Parameters ..

      INTEGER   MAXRAT
      PARAMETER ( MAXRAT = 100 )
c     ..
c     .. Scalar Arguments ..

      LOGICAL   ONLYFL
      INTEGER   MAXAZ, MAXMU, MAXTAU, MAXULV, MAXUMU, NPHI, NTAU, NUMU
c     ..
c     .. Array Arguments ..

      REAL      DFDT( * ), FLUP( * ), PHI( * ), RFLDIR( * ), RFLDN( * ),
     &          TSTDFD( * ), TSTFDN( * ), TSTFIR( * ), TSTFUP( * ),
     &          TSTUU( MAXTAU, MAXMU, MAXAZ ), UMU( * ), UTAU( * ),
     &          UU( MAXUMU, MAXULV, * )
c     ..
c     .. Local Scalars ..

      INTEGER  IU, J, LU, NUMBAD
      REAL     FLXMAX, FNOISE, RAT, RAT1, RAT2, RAT3, RAT4, UMAX, UNOISE
c     ..
c     .. Local Arrays ..

      REAL      RATV( MAXRAT )
c     ..
c     .. External Functions ..

      REAL      RATIO
      EXTERNAL  RATIO
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Statement Functions ..

      LOGICAL   BADRAT
c     ..
c     .. Statement Function definitions ..

      BADRAT( RAT ) = (RAT.LT.0.99) .OR. (RAT.GT.1.01)
c     ..


      IF ( NTAU.GT.MAXTAU .OR. NUMU.GT.MAXMU .OR. NPHI.GT.MAXAZ )  CALL
     &   ERRMSG( 'PRTFIN--out of bounds in comparator arrays', .TRUE.)

      FLXMAX = 0.0
      DO 5  LU = 1, NTAU
         FLXMAX = MAX( FLXMAX, TSTFIR(LU), TSTFDN(LU), TSTFUP(LU) )
 5    CONTINUE
      FNOISE = 1.E-6 * FLXMAX
      IF( FLXMAX.LE.0.0 )
     &    CALL ERRMSG( 'PRTFIN--all fluxes zero or negative', .FALSE.)
      IF( FNOISE.LE.0.0 )
     &    CALL ERRMSG( 'PRTFIN--all fluxes near underflowing', .FALSE.)

      NUMBAD = 0

      WRITE(*,'(//,A,/,A,/,A)')
     &  '                  <-------------- FLUXES -------------->',
     &  '    Optical       Downward       Downward         Upward'//
     &  '    d(Net Flux)',
     &  '      Depth         Direct        Diffuse        Diffuse'//
     &  '    / d(Op Dep)'

      DO 10  LU = 1, NTAU

         WRITE( *,'(0P,F11.4,1P,4E15.4)')  UTAU(LU), RFLDIR(LU),
     &          RFLDN(LU), FLUP(LU), DFDT(LU)
         RAT1 = RATIO( RFLDIR(LU), TSTFIR(LU) )
         RAT2 = RATIO( RFLDN(LU),  TSTFDN(LU) )
         RAT3 = RATIO(  FLUP(LU),  TSTFUP(LU) )
         RAT4 = RATIO(  DFDT(LU),  TSTDFD(LU) )
         WRITE( *,'(11X,4( ''    ('',F9.4,'')''))')
     &          RAT1, RAT2, RAT3, RAT4

         IF( BADRAT(RAT1) .AND. ABS(RFLDIR(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT2) .AND. ABS(RFLDN(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT3) .AND. ABS(FLUP(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1
         IF( BADRAT(RAT4) .AND. ABS(DFDT(LU)).GT.FNOISE )
     &        NUMBAD = NUMBAD+1

10    CONTINUE


      IF ( ONLYFL )  GO TO 100

      IF ( NUMU.GT.MAXRAT .OR. NPHI.GT.MAXRAT )
     &     CALL ERRMSG( 'PRTFIN--increase parameter MAXRAT', .TRUE.)


c                                       ** Print intensities

      IF ( NPHI.GT.8 ) CALL ERRMSG
     &      ( 'PRTFIN--intensity FORMATs inadequate',.FALSE.)

      UMAX = 0.0
      DO 36  LU = 1, NTAU
         DO 35  IU = 1, NUMU
            DO 34  J = 1, NPHI
               UMAX = MAX( UMAX, TSTUU(LU,IU,J) )
 34         CONTINUE
 35      CONTINUE
 36   CONTINUE
      UNOISE = 1.E-6 * UMAX
      IF( UMAX.LE.0.0 )  CALL ERRMSG
     &     ( 'PRTFIN--all intensities zero or negative',.FALSE.)
      IF( UNOISE.LE.0.0 ) CALL ERRMSG
     &     ( 'PRTFIN--all intensities near underflowing',.FALSE.)

      WRITE( *,'(//,A,//,A,/,A,/,A,8(F10.1,4X))' )
     &            ' ********  I N T E N S I T I E S  *********',
     &            '             Polar   Azimuthal Angles (Degrees)',
     &            '   Optical   Angle',
     &            '     Depth  Cosine', ( PHI(J), J = 1, NPHI )

      DO 60  LU = 1, NTAU

         DO 50  IU = 1, NUMU

            IF( IU.EQ.1 ) WRITE( *,'(/,0P,F10.3,F8.3,1P,8E14.4)')
     &          UTAU(LU), UMU(IU), ( UU( IU,LU,J ), J = 1, NPHI )

            IF( IU.GT.1 ) WRITE( *,'(10X,0P,F8.3, 1P,8E14.4)')
     &                    UMU(IU), ( UU( IU,LU,J ), J = 1, NPHI )

            DO 40  J = 1, NPHI
               RATV(J) = RATIO( UU(IU,LU,J), TSTUU(LU,IU,J) )
               IF( BADRAT(RATV(J)) .AND. ABS(UU(IU,LU,J)).GT.UNOISE )
     &                  NUMBAD = NUMBAD + 1
   40       CONTINUE

            WRITE( *,'(18X, 8(:,''   ('',F9.4,'')''))')
     &           ( RATV(J), J = 1, NPHI )

   50    CONTINUE
   60 CONTINUE

  100 CONTINUE
      IF( NUMBAD.GT.0 )  WRITE( *,300)  ' ====  ', NUMBAD,
     &    '  SERIOUSLY NON-UNIT RATIOS    ===='

      RETURN

300   FORMAT( //,1X,45('='),/,A,I4,A,/,1X,45('=') )
      END

      BLOCK DATA  CHEKDO

c       Correct answers to test problems ( as produced by DISORT
c       running entirely in double precision (14 significant digits)
c       on a Digital Alpha Workstation computer ).

c     .. Parameters ..

      INTEGER   MXPROB, MXCASE, MAXTAU, MAXMU, MAXAZ
      PARAMETER ( MXPROB = 9, MXCASE = 8, MAXTAU = 5, MAXMU = 10,
     &            MAXAZ = 3 )
c     ..
c     .. Local Scalars ..

      INTEGER   I, J, K
c     ..
c     .. Common blocks ..

      COMMON    / DOCHEK / TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU

      REAL      TSTDFD( MAXTAU, MXCASE, MXPROB ),
     &          TSTFDN( MAXTAU, MXCASE, MXPROB ),
     &          TSTFIR( MAXTAU, MXCASE, MXPROB ),
     &          TSTFUP( MAXTAU, MXCASE, MXPROB ),
     &          TSTUU( MAXTAU, MAXMU, MAXAZ, MXCASE, MXPROB )
c     ..

c ********************* Test Case 1A *********************************

      DATA (TSTFIR(I,1,1), I = 1, 2) / 3.14159E+00, 2.29844E+00 /
      DATA (TSTFDN(I,1,1), I = 1, 2) / 0.0, 7.94108E-02 /
      DATA (TSTFUP(I,1,1), I = 1, 2) / 7.99451E-02, 0.0 /
      DATA (TSTDFD(I,1,1), I = 1, 2) / 2.54067E+01, 1.86531E+01 /
      DATA ((TSTUU(I,J,1,1,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.17771E-01, 2.64170E-02, 1.34041E-02,
     &    1.33826E-02, 2.63324E-02, 1.15898E-01, 3*0.0 /

c ********************* Test Case 1B *********************************

      DATA (TSTFIR(I,2,1), I = 1, 2) / 3.14159E+00, 2.29844E+00 /
      DATA (TSTFDN(I,2,1), I = 1, 2) / 0.0, 4.20233E-01 /
      DATA (TSTFUP(I,2,1), I = 1, 2) / 4.22922E-01, 0.0 /
      DATA (TSTDFD(I,2,1), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 6.22884E-01, 1.39763E-01, 7.09192E-02,
     &    7.08109E-02, 1.39337E-01, 6.13458E-01, 3*0.0 /

c ********************* Test Case 1C *********************************

      DATA (TSTFIR(I,3,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,3,1), I = 1, 2) / 3.14159E+00, 3.04897E+00 /
      DATA (TSTFUP(I,3,1), I = 1, 2) / 9.06556E-02, 0.0 /
      DATA (TSTDFD(I,3,1), I = 1, 2) / 6.66870E-02, 5.88936E-02 /
      DATA ((TSTUU(I,J,1,3,1), J = 1, 6), I = 1, 2)
     &  / 3*1.0, 1.33177E-01, 2.99879E-02, 1.52233E-02,
     &    9.84447E-01, 9.69363E-01, 8.63946E-01, 3*0.0 /

c ********************* Test Case 1D *********************************

      DATA (TSTFIR(I,4,1), I = 1, 2) / 3.14159E+00, 0.00000E+00 /
      DATA (TSTFDN(I,4,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFUP(I,4,1), I = 1, 2) / 2.59686E-01, 0.0 /
      DATA (TSTDFD(I,4,1), I = 1, 2) / 2.57766E+01, 0.0 /
      DATA ((TSTUU(I,J,1,4,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 2.62972E-01, 9.06967E-02, 5.02853E-02,
     &    1.22980E-15, 1.30698E-17, 6.88840E-18, 3*0.0 /

c ********************* Test Case 1E *********************************

      DATA (TSTFIR(I,5,1), I = 1, 2) / 3.14159E+00, 0.00000E+00 /
      DATA (TSTFDN(I,5,1), I = 1, 2) / 0.0, 6.76954E-02 /
      DATA (TSTFUP(I,5,1), I = 1, 2) / 3.07390E+00, 0.0 /
      DATA (TSTDFD(I,5,1), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,5,1), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.93321E+00, 1.02732E+00, 7.97199E-01,
     &    2.71316E-02, 1.87805E-02, 1.16385E-02, 3*0.0 /

c ********************* Test Case 1F *********************************

      DATA (TSTFIR(I,6,1), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,6,1), I = 1, 2) / 3.14159E+00, 4.60048E-03 /
      DATA (TSTFUP(I,6,1), I = 1, 2) / 2.49618E+00, 0.0 /
      DATA (TSTDFD(I,6,1), I = 1, 2) / 1.14239E-01, 7.93633E-05 /
      DATA ((TSTUU(I,J,1,6,1), J = 1, 6), I = 1, 2)
     &  / 3*1.0, 8.77510E-01, 8.15136E-01, 7.52715E-01,
     &    1.86840E-03, 1.26492E-03, 7.79280E-04, 3*0.0 /


c ********************* Test Case 2A *********************************

      DATA (TSTFIR(I,1,2), I = 1, 2) / 2.52716E-01, 2.10311E-02 /
      DATA (TSTFDN(I,1,2), I = 1, 2) / 0.0, 4.41791E-02 /
      DATA (TSTFUP(I,1,2), I = 1, 2) / 5.35063E-02, 0.0 /
      DATA (TSTDFD(I,1,2), I = 1, 2) / 1.66570E+00, 1.89848E-01 /
      DATA ((TSTUU(I,J,1,1,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.61796E-01, 2.11501E-02, 7.86713E-03,
     &    7.71897E-03, 2.00778E-02, 2.57685E-02, 3*0.0 /

c ********************* Test Case 2B *********************************

      DATA (TSTFIR(I,2,2), I = 1, 2) / 2.52716E-01, 2.10311E-02 /
      DATA (TSTFDN(I,2,2), I = 1, 2) / 0.0, 1.06123E-01 /
      DATA (TSTFUP(I,2,2), I = 1, 2) / 1.25561E-01, 0.0 /
      DATA (TSTDFD(I,2,2), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 3.47678E-01, 4.87120E-02, 1.89387E-02,
     &    1.86027E-02, 4.64061E-02, 6.77603E-02, 3*0.0 /

c ********************* Test Case 2C *********************************

      DATA (TSTFIR(I,3,2), I = 1, 2) / 2.52716E-01, 2.56077E-28 /
      DATA (TSTFDN(I,3,2), I = 1, 2) / 0.0, 2.51683E-04 /
      DATA (TSTFUP(I,3,2), I = 1, 2) / 6.24730E-02, 0.0 /
      DATA (TSTDFD(I,3,2), I = 1, 2) / 1.67462E+00, 1.75464E-04 /
      DATA ((TSTUU(I,J,1,3,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 1.62566E-01, 2.45786E-02, 1.01498E-02,
     &    1.70004E-04, 3.97168E-05, 1.32472E-05, 3*0.0 /

c ********************* Test Case 2D *********************************

      DATA (TSTFIR(I,4,2), I = 1, 2) / 2.52716E-01, 0.0 /
      DATA (TSTFDN(I,4,2), I = 1, 2) / 0.0, 2.68008E-02 /
      DATA (TSTFUP(I,4,2), I = 1, 2) / 2.25915E-01, 0.0 /
      DATA (TSTDFD(I,4,2), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,4,2), J = 1, 6), I = 1, 2)
     &  / 3*0.0, 3.64010E-01, 8.26993E-02, 4.92370E-02,
     &    1.05950E-02, 7.69337E-03, 3.79276E-03, 3*0.0 /


c ********************* Test Case 3A *********************************

      DATA (TSTFIR(I,1,3), I = 1, 2) / 3.14159E+00, 1.15573E+00 /
      DATA (TSTFDN(I,1,3), I = 1, 2) / 0.0, 1.73849E+00 /
      DATA (TSTFUP(I,1,3), I = 1, 2) / 2.47374E-01, 0.0 /
      DATA (TSTDFD(I,1,3), I = 1, 2) / 0.0, 0.0 /
      DATA ((TSTUU(I,J,1,1,3), J = 1, 6), I = 1, 2) /
     &    3*0.0, 1.51159E-01, 1.01103E-01, 3.95460E-02,
     &    3.05855E+00, 2.66648E-01, 2.13750E-01, 3*0.0 /

c ********************* Test Case 3B *********************************

      DATA (TSTFIR(I,2,3), I = 1, 2) / 3.14159E+00, 1.05389E-03 /
      DATA (TSTFDN(I,2,3), I = 1, 2) / 0.0, 1.54958E+00 /
      DATA (TSTFUP(I,2,3), I = 1, 2) / 1.59096E+00, 0.0 /
      DATA (TSTDFD(I,2,3), I = 1, 2) / 2*0.0 /
      DATA ((TSTUU(I,J,1,2,3), J = 1, 6), I = 1, 2) /
     &    3*0.0, 3.79740E-01, 5.19598E-01, 4.93302E-01,
     &    6.69581E-01, 4.22350E-01, 2.36362E-01, 3*0.0 /


c ********************* Test Case 4A *********************************

      DATA (TSTFIR(I,1,4), I = 1, 3)
     &  / 3.14159E+00, 1.90547E+00, 1.15573E+00 /
      DATA (TSTFDN(I,1,4), I = 1, 3)
     &  / 0.0, 1.17401E+00, 1.81264E+00 /
      DATA (TSTFUP(I,1,4), I = 1, 3)
     &  / 1.73223E-01, 1.11113E-01, 0.0 /
      DATA (TSTDFD(I,1,4), I = 1, 3) / 3*0.0 /
      DATA ((TSTUU(I,J,1,1,4), J = 1, 6), I = 1, 3)
     &  / 3*0.0, 9.26837E-02,
     &    6.59569E-02, 3.64755E-02, 2.51608E+00, 1.19287E-01,
     &    1.34962E-01, 1.23887E-01, 4.02058E-02, 1.77746E-02,
     &    3.37302E+00, 2.19835E-01, 1.56893E-01, 3*0.0 /

c ********************* Test Case 4B *********************************

      DATA (TSTFIR(I,2,4), I = 1, 3)
     &  / 3.14159E+00, 1.90547E+00, 1.15573E+00 /
      DATA (TSTFDN(I,2,4), I = 1, 3)
     &  / 0.0, 1.01517E+00, 1.51554E+00 /
      DATA (TSTFUP(I,2,4), I = 1, 3)
     &  / 1.23665E-01, 7.88690E-02, 0.0 /
      DATA (TSTDFD(I,2,4), I = 1, 3)
     &  / 3.43724E-01, 3.52390E-01, 3.19450E-01 /
      DATA ((TSTUU(I,J,1,2,4), J = 1, 6), I = 1, 3)
     &  / 3*0.0, 6.53056E-02,
     &    4.55144E-02, 2.82693E-02, 2.24258E+00, 9.66049E-02,
     &    9.61335E-02, 8.43278E-02, 2.79473E-02, 1.38835E-02,
     &    2.97057E+00, 1.67698E-01, 1.08115E-01, 3*0.0 /

c ********************* Test Case 4C *********************************

      DATA (TSTFIR(I,3,4), I = 1, 3)
     &  / 1.57080E+00, 5.77864E-01, 2.12584E-01 /
      DATA (TSTFDN(I,3,4), I = 1, 3)
     &  / 0.0, 7.02764E-01, 8.03294E-01 /
      DATA (TSTFUP(I,3,4), I = 1, 3)
     &  / 2.25487E-01, 1.23848E-01, 0.0 /
      DATA (TSTDFD(I,3,4), I = 1, 3)
     &  / 3.85003E-01, 3.37317E-01, 2.16403E-01 /
      DATA (((TSTUU(I,J,K,3,4), J = 1, 6), I = 1, 3), K = 1, 3)
     &  / 3*0.0, 8.70812E-01,
     &    2.24960E-01, 2.27572E-02, 4.77016E-02, 3.02631E+00,
     &    1.41195E+00, 6.97692E-01, 1.09130E-01, 9.32861E-03,
     &    8.38488E-02, 2.70538E+00, 8.76523E-01, 6*0.0,
     &    8.88117E-02, 5.77411E-02, 2.27572E-02,
     &    4.77016E-02, 5.80971E-02, 1.04502E-01, 9.16071E-02,
     &    2.95842E-02, 9.32861E-03, 8.38488E-02, 9.42187E-02,
     &    8.95457E-02, 6*0.0, 6.98247E-02,
     &    5.02877E-02, 2.27572E-02, 4.77016E-02, 2.58544E-02,
     &    6.25954E-02, 5.91273E-02, 2.47702E-02, 9.32861E-03,
     &    8.38488E-02, 3.99383E-02, 4.67155E-02, 3*0.0 /


c ********************* Test Case 5A *********************************

      DATA (TSTFIR(I,1,5), I = 1, 3)
     &  / 3.14159E+00, 3.97856E-14, 5.03852E-28 /
      DATA (TSTFDN(I,1,5), I = 1, 3)
     &  / 0.0, 2.24768E+00, 4.79851E-01 /
      DATA (TSTFUP(I,1,5), I = 1, 3)
     &  / 2.66174E+00, 1.76783E+00, 0.0 /
      DATA (TSTDFD(I,1,5), I = 1, 3) / 3*0.0 /
      DATA ((TSTUU(I,J,1,1,5), J = 1, 6), I = 1, 3)
     &  / 3*0.0, 4.58927E-01,
     &    7.72983E-01, 1.07196E+00, 7.53662E-01, 6.96362E-01,
     &    6.50541E-01, 6.27631E-01, 5.81809E-01, 5.24532E-01,
     &    1.95230E-01, 1.31990E-01, 7.20655E-02, 3*0.0 /

c ********************* Test Case 5B *********************************

      DATA (TSTFIR(I,2,5), I = 1, 3)
     &  / 1.28058E-01, 8.67322E-06, 4.47729E-21 /
      DATA (TSTFDN(I,2,5), I = 1, 3)
     &  / 1.74767E+00, 2.33975E-01, 6.38345E-05 /
      DATA (TSTFUP(I,2,5), I = 1, 3)
     &  / 2.70485E-01, 3.74252E-02, 1.02904E-05 /
      DATA (TSTDFD(I,2,5), I = 1, 3)
     &  / 3.10129E-01, 4.52671E-02, 1.25021E-05 /
      DATA ((TSTUU(I,J,1,2,5), J = 1, 6), I = 1, 3)
     &  / 6.79623E+01, 2.21027E-01, 1.36619E-01, 1.14084E-01,
     &    8.73870E-02, 8.81626E-02, 2.05706E-01, 4.92736E-02,
     &    2.65449E-02, 2.02154E-02, 1.29661E-02, 9.51334E-03,
     &    3.41286E-05, 1.39916E-05, 7.47039E-06, 5.65602E-06,
     &    3.58245E-06, 2.57858E-06 /

c ********************* Test Case 6A *********************************

      DATA (TSTFIR(I,1,6), I = 1, 2) / 2*100.0 /
      DATA (TSTFDN(I,1,6), I = 1, 2) / 2*0.0 /
      DATA (TSTFUP(I,1,6), I = 1, 2) / 2*0.0 /
      DATA (TSTDFD(I,1,6), I = 1, 2) / 2*200.0 /
      DATA ((TSTUU(I,J,1,1,6), J = 1, 4), I = 1, 2) / 8*0.0 /

c ********************* Test Case 6B *********************************

      DATA (TSTFIR(I,2,6), I = 1, 3)
     &  / 1.000000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,2,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,2,6), I = 1, 3) / 3*0.0 /
      DATA (TSTDFD(I,2,6), I = 1, 3)
     &  / 2.00000E+02, 7.35759E+01, 2.70671E+01 /
      DATA ((TSTUU(I,J,1,2,6), J = 1, 4), I = 1, 3) / 12*0.0 /

c ********************* Test Case 6C *********************************

      DATA (TSTFIR(I,3,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,3,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,3,6), I = 1, 3)
     &  / 1.48450E+00, 2.99914E+00, 6.76676E+00 /
      DATA (TSTDFD(I,3,6), I = 1, 3)
     &  / 2.02010E+02, 7.79962E+01, 4.06006E+01 /
      DATA ((TSTUU(I,J,1,3,6), J = 1, 4), I = 1, 3)
     &  / 2*0.0, 9.77882E-05, 7.92386E-01,
     &    2*0.0, 1.45131E-02, 1.30642E+00,
     &    2*0.0, 2.15393E+00, 2.15393E+00 /

c ********************* Test Case 6D *********************************

      DATA (TSTFIR(I,4,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,4,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,4,6), I = 1, 3)
     &  / 6.70783E-01, 1.39084E+00, 3.31655E+00 /
      DATA (TSTDFD(I,4,6), I = 1, 3)
     &  / 2.00936E+02, 7.57187E+01, 3.45317E+01 /
      DATA ((TSTUU(I,J,1,4,6), J = 1, 4), I = 1, 3)
     &  / 2*0.0, 6.80068E-05, 3.15441E-01,
     &    2*0.0, 1.00931E-02, 5.20074E-01,
     &    2*0.0, 1.49795E+00, 8.57458E-01 /

c ********************* Test Case 6E *********************************

      DATA (TSTFIR(I,5,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,5,6), I = 1, 3) / 3*0.0 /
      DATA (TSTFUP(I,5,6), I = 1, 3)
     &  / 7.95458E+01, 1.59902E+02, 3.56410E+02 /
      DATA (TSTDFD(I,5,6), I = 1, 3)
     &  / 3.07079E+02, 3.07108E+02, 7.17467E+02 /
      DATA ((TSTUU(I,J,1,5,6), J = 1, 4), I = 1, 3)
     &  / 2*0.0, 4.53789E-03, 4.33773E+01,
     &    2*0.0, 6.73483E-01, 7.15170E+01,
     &    2*0.0, 9.99537E+01, 1.17912E+02 /

c ********************* Test Case 6F *********************************

      DATA (TSTFIR(I,6,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,6,6), I = 1, 3)
     &  / 3.21497E+02, 1.42493E+02, 7.05305E+01 /
      DATA (TSTFUP(I,6,6), I = 1, 3)
     &  / 8.27917E+01, 1.66532E+02, 3.71743E+02 /
      DATA (TSTDFD(I,6,6), I = 1, 3)
     &  / 9.54523E+02, 5.27085E+02, 8.45341E+02 /
      DATA ((TSTUU(I,J,1,6,6), J = 1, 4), I = 1, 3)
     &  / 1.02336E+02, 1.02336E+02, 4.80531E-03, 4.50168E+01,
     &    6.20697E+01, 6.89532E-01, 7.13172E-01, 7.42191E+01,
     &    3.76472E+01, 4.64603E-03, 1.05844E+02, 1.22368E+02 /

c ********************* Test Case 6G *********************************

      DATA (TSTFIR(I,7,6), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,7,6), I = 1, 3)
     &  / 3.21497E+02, 3.04775E+02, 3.63632E+02 /
      DATA (TSTFUP(I,7,6), I = 1, 3)
     &  / 3.35292E+02, 4.12540E+02, 4.41125E+02 /
      DATA (TSTDFD(I,7,6), I = 1, 3)
     &  / 5.80394E+02, 1.27117E+02, -1.68003E+02 /
      DATA ((TSTUU(I,J,1,7,6), J = 1, 4), I = 1, 3)
     &  / 1.02336E+02, 1.02336E+02, 7.80733E+01, 1.16430E+02,
     &    9.78748E+01, 1.01048E+02, 1.15819E+02, 1.34966E+02,
     &    1.10061E+02, 1.38631E+02, 1.38695E+02, 1.40974E+02 /

c ********************* Test Case 6H *********************************

      DATA (TSTFIR(I,8,6), I = 1, 3)
     &  / 1.00000E+02, 1.35335E+01, 2.06115E-07 /
      DATA ( TSTFDN(I,8,6), I = 1, 3)
     &  / 3.21497E+02, 2.55455E+02, 4.43444E+02 /
      DATA (TSTFUP(I,8,6), I = 1, 3)
     &  / 2.37350E+02, 2.61130E+02, 4.55861E+02 /
      DATA (TSTDFD(I,8,6), I = 1, 3)
     &  / 4.23780E+02, 6.19828E+01, -3.11658E+01 /
      DATA ((TSTUU(I,J,1,8,6), J = 1, 4), I = 1, 3)
     &  / 1.02336E+02, 1.02336E+02, 7.12616E+01, 7.80736E+01,
     &    8.49992E+01, 7.73186E+01, 7.88310E+01, 8.56423E+01,
     &    1.38631E+02, 1.45441E+02, 1.44792E+02, 1.45163E+02 /


c ********************* Test Case 7A *********************************

      DATA (TSTFIR(I,1,7), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,1,7), I = 1, 2) / 0.0, 1.21204E+02 /
      DATA (TSTFUP(I,1,7), I = 1, 2) / 8.62936E+01,  0.0 /
      DATA (TSTDFD(I,1,7), I = 1, 2) /-5.13731E+01,-5.41036E+02 /

c ********************* Test Case 7B *********************************

      DATA (TSTFIR(I,2,7), I = 1, 2) / 2*0.0 /
      DATA (TSTFDN(I,2,7), I = 1, 2) / 0.0, 2.07786E-05 /
      DATA (TSTFUP(I,2,7), I = 1, 2) / 1.10949E-06,  0.0 /
      DATA (TSTDFD(I,2,7), I = 1, 2) / 8.23219E-08, -5.06461E-06 /
      DATA ((TSTUU(I,J,1,2,7), J = 1, 2), I = 1, 2) /
     &    0.00000E+00, 4.65744E-07, 7.52311E-06, 0.00000E+00 /

c ********************* Test Case 7C *********************************

      DATA (TSTFIR(I,3,7), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,3,7), I = 1, 3)
     &  / 3.19830E+02, 3.54099E+02, 3.01334E+02 /
      DATA (TSTFUP(I,3,7), I = 1, 3)
     &  / 4.29572E+02, 4.47018E+02, 5.94576E+02 /
      DATA (TSTDFD(I,3,7), I = 1, 3)
     &  /-8.04270E+01, 2.51589E+02, 7.15964E+02 /
      DATA (((TSTUU(I,J,K,3,7), J = 1, 4), I = 1, 3), K = 1, 2)
     &  / 1.01805E+02, 1.01805E+02, 1.46775E+02, 1.49033E+02,
     &    1.06583E+02, 1.28565E+02, 1.04464E+02, 1.59054E+02,
     &    9.66519E+01, 8.65854E+01, 1.89259E+02, 1.89259E+02,
     &    1.01805E+02, 1.01805E+02, 1.29641E+02, 1.49033E+02,
     &    1.06583E+02, 1.06408E+02, 9.48418E+01, 1.59054E+02,
     &    9.66519E+01, 7.49310E+01, 1.89259E+02, 1.89259E+02 /

c ********************* Test Case 7D *********************************

      DATA (TSTFIR(I,4,7), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,4,7), I = 1, 3)
     &  / 3.19830E+02, 3.50555E+02, 2.92063E+02 /
      DATA (TSTFUP(I,4,7), I = 1, 3)
     &  / 3.12563E+02, 2.68126E+02, 3.05596E+02 /
      DATA (TSTDFD(I,4,7), I = 1, 3)
     &  /-1.68356E+02, 1.01251E+02, 4.09326E+02 /
      DATA (((TSTUU(I,J,K,4,7), J = 1, 4), I = 1, 3), K = 1, 2)
     &  / 1.01805E+02, 1.01805E+02, 1.40977E+02, 9.62764E+01,
     &    1.06203E+02, 1.23126E+02, 9.19545E+01, 8.89528E+01,
     &    9.56010E+01, 7.25576E+01, 9.72743E+01, 9.72743E+01,
     &    1.01805E+02, 1.01805E+02, 1.23843E+02, 9.62764E+01,
     &    1.06203E+02, 1.00969E+02, 8.23318E+01, 8.89528E+01,
     &    9.56010E+01, 6.09031E+01, 9.72743E+01, 9.72743E+01 /

c ********************* Test Case 7E *********************************

      DATA (TSTFIR(I,5,7), I = 1, 3)
     &  / 1.00000E+02, 3.67879E+01, 1.35335E+01 /
      DATA (TSTFDN(I,5,7), I = 1, 3)
     &  / 3.19830E+02, 3.53275E+02, 2.99002E+02 /
      DATA (TSTFUP(I,5,7), I = 1, 3)
     &  / 4.04300E+02, 4.07843E+02, 5.29248E+02 /
      DATA (TSTDFD(I,5,7), I = 1, 3)
     &  /-9.98568E+01, 2.17387E+02, 6.38461E+02 /
      DATA (((TSTUU(I,J,K,5,7), J = 1, 4), I = 1, 3), K = 1, 2)
     &  / 1.01805E+02, 1.01805E+02, 1.45448E+02, 1.38554E+02,
     &    1.06496E+02, 1.27296E+02, 1.01395E+02, 1.45229E+02,
     &    9.63993E+01, 8.29009E+01, 1.60734E+02, 1.71307E+02,
     &    1.01805E+02, 1.01805E+02, 1.28281E+02, 1.38554E+02,
     &    1.06496E+02, 1.05111E+02, 9.16726E+01, 1.45229E+02,
     &    9.63993E+01, 7.11248E+01, 1.59286E+02, 1.71307E+02 /

c ********************* Test Case 8A *********************************

      DATA (TSTFIR(I,1,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,1,8), I = 1, 3) /
     &  1.00000E+00, 7.22235E-01, 5.13132E-01 /
      DATA (TSTFUP(I,1,8), I = 1, 3) / 9.29633E-02, 2.78952E-02, 0.0 /
      DATA (TSTDFD(I,1,8), I = 1, 3) /
     &  1.12474E+00, 6.51821E-01, 5.63361E-01 /
      DATA ((TSTUU(I,J,1,1,8), J = 1, 4), I = 1, 3) /
     &  2*3.18310E-01, 5.62566E-02, 1.94423E-02,
     &  2.62711E-01, 1.36952E-01, 1.84909E-02, 5.52188E-03,
     &  2.10014E-01, 5.60376E-02, 2*0.0 /

c ********************* Test Case 8B *********************************

      DATA (TSTFIR(I,2,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,2,8), I = 1, 3) /
     &  1.00000E+00, 7.95332E-01, 6.50417E-01 /
      DATA (TSTFUP(I,2,8), I = 1, 3) / 2.25136E-01, 1.26349E-01, 0.0 /
      DATA (TSTDFD(I,2,8), I = 1, 3) /
     &  5.12692E-01, 3.56655E-01, 5.68095E-02 /
      DATA ((TSTUU(I,J,1,2,8), J = 1, 4), I = 1, 3) /
     &  2*3.18310E-01, 1.23687E-01, 4.95581E-02,
     &  2.77499E-01, 1.83950E-01, 8.35695E-02, 2.50575E-02,
     &  2.40731E-01, 1.29291E-01, 2*0.0 /

c ********************* Test Case 8C *********************************

      DATA (TSTFIR(I,3,8), I = 1, 3) / 3*0.0 /
      DATA (TSTFDN(I,3,8), I = 1, 3) /
     &  1.00000E+00, 4.86157E-01, 1.59984E-01 /
      DATA (TSTFUP(I,3,8), I = 1, 3) / 3.78578E-01, 2.43397E-01, 0.0 /
      DATA (TSTDFD(I,3,8), I = 1, 3) /
     &  5.65095E-01, 2.76697E-01, 1.35679E-02 /
      DATA ((TSTUU(I,J,1,3,8), J = 1, 4), I = 1, 3) /
     &  2*3.18310E-01, 1.49335E-01, 1.04766E-01,
     &  1.89020E-01, 9.88158E-02, 9.65192E-02, 6.54445E-02,
     &  6.84762E-02, 2.96698E-02, 2*0.0 /


c ********************* Test Case 9A *********************************

      DATA (TSTFIR(I,1,9), I = 1, 5) / 5*0.0 /
      DATA (TSTFDN(I,1,9), I = 1, 5) /
     &  1.00000E+00, 3.55151E-01, 1.44265E-01, 6.71445E-03, 6.16968E-07/
      DATA (TSTFUP(I,1,9), I = 1, 5) /
     &  2.27973E-01, 8.75098E-02, 3.61819E-02, 2.19291E-03, 0.0 /
      DATA (TSTDFD(I,1,9), I = 1, 5) /
     &  8.82116E-01, 2.32366E-01, 9.33443E-02, 3.92782E-03, 1.02500E-07/
      DATA ((TSTUU(I,J,1,1,9), J = 1, 4), I = 1, 5) /
     &  2*3.18310E-01, 9.98915E-02, 5.91345E-02,
     &  1.53507E-01, 5.09531E-02, 3.67006E-02, 2.31903E-02,
     &  7.06614E-02, 2.09119E-02, 1.48545E-02, 9.72307E-03,
     &  3.72784E-03, 1.08815E-03, 8.83316E-04, 5.94743E-04,
     &  2.87656E-07, 1.05921E-07, 2*0.0 /

c ********************* Test Case 9B *********************************

      DATA (TSTFIR(I,2,9), I = 1, 5) / 5*0.0 /
      DATA (TSTFDN(I,2,9), I = 1, 5) /
     &  1.00000E+00, 4.52357E-01, 2.36473E-01, 2.76475E-02, 7.41853E-05/
      DATA (TSTFUP(I,2,9), I = 1, 5) /
     &  1.00079E-01, 4.52014E-02, 2.41941E-02, 4.16016E-03, 0.0 /
      DATA (TSTDFD(I,2,9), I = 1, 5) /
     &  8.04577E-01, 2.55330E-01, 1.30976E-01, 1.36227E-02, 1.22022E-05/
      DATA ((TSTUU(I,J,1,2,9), J = 1, 4), I = 1, 5) /
     &  2*3.18310E-01, 7.39198E-02, 1.32768E-02,
     &  1.96609E-01, 5.92369E-02, 3.00230E-02, 7.05566E-03,
     &  1.15478E-01, 3.01809E-02, 1.52672E-02, 4.06932E-03,
     &  1.46177E-02, 3.85590E-03, 2.38301E-03, 7.77890E-04,
     &  3.37742E-05, 1.20858E-05, 2*0.0 /

c ********************* Test Case 9C *********************************

      DATA (TSTFIR(I,3,9), I = 1, 5) /
     &  1.57080E+00, 1.92354E-01, 2.35550E-02, 9.65131E-06, 9.03133E-19/
      DATA (TSTFDN(I,3,9), I = 1, 5 ) /
     &  6.09217E+00, 4.97279E+00, 4.46616E+00, 4.22731E+00, 4.73767E+00/
      DATA (TSTFUP(I,3,9), I = 1, 5) /
     &  4.68414E+00, 4.24381E+00, 4.16941E+00, 4.30667E+00, 5.11524E+00/
      DATA (TSTDFD(I,3,9), I = 1, 5) /
     &  3.49563E+00, 8.81206E-01, 3.50053E-01, 1.93471E-02, 7.15349E-02/
      DATA (((TSTUU(I,J,K,3,9), J = 1, 4), I = 1, 5), K = 1, 3)
     & / 1.93920E+00, 1.93920E+00, 1.61855E+00, 1.43872E+00,
     &   1.66764E+00, 1.44453E+00, 1.38339E+00, 1.33890E+00,
     &   1.48511E+00, 1.35009E+00, 1.33079E+00, 1.32794E+00,
     &   1.34514E+00, 1.35131E+00, 1.35980E+00, 1.37918E+00,
     &   1.48927E+00, 1.54270E+00, 1.62823E+00, 1.62823E+00,
     &   1.93920E+00, 1.93920E+00, 1.57895E+00, 1.43872E+00,
     &   1.66764E+00, 1.42925E+00, 1.37317E+00, 1.33890E+00,
     &   1.48511E+00, 1.34587E+00, 1.32921E+00, 1.32794E+00,
     &   1.34514E+00, 1.35129E+00, 1.35979E+00, 1.37918E+00,
     &   1.48927E+00, 1.54270E+00, 1.62823E+00, 1.62823E+00,
     &   1.93920E+00, 1.93920E+00, 1.56559E+00, 1.43872E+00,
     &   1.66764E+00, 1.42444E+00, 1.37034E+00, 1.33890E+00,
     &   1.48511E+00, 1.34469E+00, 1.32873E+00, 1.32794E+00,
     &   1.34514E+00, 1.35128E+00, 1.35979E+00, 1.37918E+00,
     &   1.48927E+00, 1.54270E+00, 2*1.62823E+00 /

      END
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: ErrPack.f,v 2.1 2000/03/27 21:40:49 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE  ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error

      LOGICAL       FATAL, MsgLim
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /


      IF ( FATAL )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         WRITE ( *,99 )
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     &   'They will no longer be printed  <<<<<<<', // )
      END

      LOGICAL FUNCTION  WrtBad ( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

      CHARACTER*(*)  VarNam
      INTEGER        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,
     &                     '  in error  ****'
      IF ( NumMsg.EQ.MaxMsg )
     &   CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )

      RETURN
      END

      LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

      CHARACTER*(*)  DimNam
      INTEGER        MinVal


      WRITE ( *, '(/,3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,
     &                     '  should be increased to at least ', MinVal
      WrtDim = .TRUE.

      RETURN
      END

      LOGICAL FUNCTION  TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

      CHARACTER*(*)  VarNam
      REAL           RelErr


      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' )
     &       ' Output variable ', VarNam,' differed by ', 100.*RelErr,
     &       ' per cent from correct value.  Self-test failed.'

      RETURN
      END

c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: LINPAK.f,v 2.1 2000/03/27 21:40:49 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c Call tree:
c
c    SGBCO
c       SASUM
c       SDOT
c       SAXPY
c       SGBFA
c           ISAMAX
c           SAXPY
c           SSCAL
c       SSCAL
c   SGBSL
c       SDOT
c       SAXPY
c   SGECO
c       SASUM
c       SDOT
c       SAXPY
c       SGEFA
c           ISAMAX
c           SAXPY
c           SSCAL
c       SSCAL
c   SGESL
c       SDOT
c       SAXPY
c   SSWAP
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      SUBROUTINE SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )

c         Factors a real band matrix by Gaussian elimination
c         and estimates the condition of the matrix.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     If  RCOND  is not needed, SGBFA is slightly faster.
c     To solve  A*X = B , follow SBGCO by SGBSL.

c     input:

C        ABD     REAL(LDA, N)
c                contains the matrix in band storage.  The columns
c                of the matrix are stored in the columns of  ABD  and
c                the diagonals of the matrix are stored in rows
c                ML+1 through 2*ML+MU+1 of  ABD .
c                See the comments below for details.

C        LDA     INTEGER
c                the leading dimension of the array  ABD .
c                LDA must be .GE. 2*ML + MU + 1 .

C        N       INTEGER
c                the order of the original matrix.

C        ML      INTEGER
c                number of diagonals below the main diagonal.
c                0 .LE. ML .LT. N .

C        MU      INTEGER
c                number of diagonals above the main diagonal.
c                0 .LE. MU .LT. N .
c                more efficient if  ML .LE. MU .

c     on return

c        ABD     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                The factorization can be written  A = L*U  where
c                L  is a product of permutation and unit lower
c                triangular matrices and  U  is upper triangular.

C        IPVT    INTEGER(N)
c                an integer vector of pivot indices.

C        RCOND   REAL
c                an estimate of the reciprocal condition of  A .
c                For the system  A*X = B , relative perturbations
c                in  A  and  B  of size  epsilon  may cause
c                relative perturbations in  X  of size  epsilon/RCOND .
c                If  RCOND  is so small that the logical expression
c                           1.0 + RCOND .EQ. 1.0
c                is true, then  A  may be singular to working
c                precision.  In particular,  RCOND  is zero  if
c                exact singularity is detected or the estimate
c                underflows.

C        Z       REAL(N)
c                a work vector whose contents are usually unimportant.
c                If  A  is close to a singular matrix, then  Z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .

c     Band storage

c           If  A  is a band matrix, the following program segment
c           will set up the input.

c                   ML = (band width below the diagonal)
c                   MU = (band width above the diagonal)
c                   M = ML + MU + 1
c                   DO 20 J = 1, N
c                      I1 = MAX(1, J-MU)
c                      I2 = MIN(N, J+ML)
c                      DO 10 I = I1, I2
c                         K = I - J + M
c                         ABD(K,J) = A(I,J)
c                10    CONTINUE
c                20 CONTINUE

c           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
c           In addition, the first  ML  rows in  ABD  are used for
c           elements generated during the triangularization.
c           The total number of rows needed in  ABD  is  2*ML+MU+1 .
c           The  ML+MU by ML+MU  upper left triangle and the
c           ML by ML  lower right triangle are not referenced.

c     Example:  if the original matrix is

c           11 12 13  0  0  0
c           21 22 23 24  0  0
c            0 32 33 34 35  0
c            0  0 43 44 45 46
c            0  0  0 54 55 56
c            0  0  0  0 65 66

c      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain

c            *  *  *  +  +  +  , * = not used
c            *  * 13 24 35 46  , + = used for pivoting
c            * 12 23 34 45 56
c           11 22 33 44 55 66
c           21 32 43 54 65  *

c --------------------------------------------------------------------


c     .. Scalar Arguments ..

      INTEGER   LDA, ML, MU, N
      REAL      RCOND
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ABD( LDA, * ), Z( * )
c     ..
c     .. Local Scalars ..

      INTEGER   INFO, IS, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM
      REAL      ANORM, EK, S, SM, T, WK, WKM, YNORM
c     ..
c     .. External Functions ..

      REAL      SASUM, SDOT
      EXTERNAL  SASUM, SDOT
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY, SGBFA, SSCAL
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MAX, MIN, SIGN
c     ..


c                       ** compute 1-norm of A
      ANORM  = 0.0E0
      L  = ML + 1
      IS = L + MU

      DO 10 J = 1, N

         ANORM  = MAX( ANORM, SASUM( L,ABD( IS,J ),1 ) )

         IF( IS.GT.ML + 1 ) IS = IS - 1

         IF( J.LE.MU ) L  = L + 1

         IF( J.GE.N - ML ) L  = L - 1

   10 CONTINUE
c                                               ** factor

      CALL SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

c     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
c     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E.
c     trans(A) is the transpose of A.  The components of E  are
c     chosen to cause maximum local growth in the elements of W  where
c     trans(U)*W = E.  The vectors are frequently rescaled to avoid
c     overflow.

c                     ** solve trans(U)*W = E
      EK = 1.0E0

      DO 20 J = 1, N
         Z( J ) = 0.0E0
   20 CONTINUE


      M  = ML + MU + 1
      JU = 0

      DO 50 K = 1, N

         IF( Z( K ).NE.0.0E0 ) EK = SIGN( EK, -Z( K ) )

         IF( ABS( EK - Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN

            S  = ABS( ABD( M,K ) ) / ABS( EK - Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            EK = S*EK

         END IF

         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )

         IF( ABD( M,K ).NE.0.0E0 ) THEN

            WK   = WK / ABD( M, K )
            WKM  = WKM / ABD( M, K )

         ELSE

            WK   = 1.0E0
            WKM  = 1.0E0

         END IF

         KP1  = K + 1
         JU   = MIN( MAX( JU,MU + IPVT( K ) ), N )
         MM   = M

         IF( KP1.LE.JU ) THEN

            DO 30 J = KP1, JU
               MM     = MM - 1
               SM     = SM + ABS( Z( J ) + WKM*ABD( MM,J ) )
               Z( J ) = Z( J ) + WK*ABD( MM, J )
               S      = S + ABS( Z( J ) )
   30       CONTINUE

            IF( S.LT.SM ) THEN

               T  = WKM - WK
               WK = WKM
               MM = M

               DO 40 J = KP1, JU
                  MM = MM - 1
                  Z( J ) = Z( J ) + T*ABD( MM, J )
   40          CONTINUE

            END IF

         END IF

         Z( K ) = WK

   50 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

c                         ** solve trans(L)*Y = W
      DO 60 KB = 1, N
         K  = N + 1 - KB
         LM = MIN( ML, N - K )

         IF( K.LT.N )
     &       Z( K ) = Z( K ) + SDOT( LM, ABD( M+1, K ), 1, Z( K+1 ), 1 )

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

         END IF

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T

   60 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = 1.0E0
c                         ** solve L*V = Y
      DO 70 K = 1, N

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
         LM     = MIN( ML, N - K )

         IF( K.LT.N )
     &       CALL SAXPY( LM, T, ABD( M+1, K ), 1, Z( K+1 ), 1 )

         IF( ABS( Z(K) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z(K) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

   70 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = S*YNORM

c                           ** solve  U*Z = W
      DO 80 KB = 1, N

         K  = N + 1 - KB

         IF( ABS( Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN

            S  = ABS( ABD( M,K ) ) / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

         IF( ABD( M,K ).NE.0.0E0 ) Z( K ) = Z( K ) / ABD( M, K )
         IF( ABD( M,K ).EQ.0.0E0 ) Z( K ) = 1.0E0

         LM = MIN( K, M ) - 1
         LA = M - LM
         LZ = K - LM
         T  = -Z( K )

         CALL SAXPY( LM, T, ABD( LA,K ), 1, Z( LZ ), 1 )

   80 CONTINUE
c                              ** make znorm = 1.0

      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = S*YNORM
      IF( ANORM.NE.0.0E0 ) RCOND  = YNORM / ANORM
      IF( ANORM.EQ.0.0E0 ) RCOND  = 0.0E0

      END

      SUBROUTINE SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

c         Factors a real band matrix by elimination.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     SGBFA is usually called by SBGCO, but it can be called
c     directly with a saving in time if  RCOND  is not needed.

c     Input:  same as SGBCO

c     On return:

c        ABD,IPVT    same as SGBCO

c        INFO    INTEGER
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  This is not an error
c                     condition for this subroutine, but it does
c                     indicate that SGBSL will divide by zero if
c                     called.  Use  RCOND  in SBGCO for a reliable
c                     indication of singularity.

c     (see SGBCO for description of band storage mode)

c ----------------------------------------------------------------


c     .. Scalar Arguments ..

      INTEGER   INFO, LDA, ML, MU, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ABD( LDA, * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, I0, J, J0, J1, JU, JZ, K, KP1, L, LM, M, MM, NM1
      REAL      T
c     ..
c     .. External Functions ..

      INTEGER   ISAMAX
      EXTERNAL  ISAMAX
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY, SSCAL
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MAX, MIN
c     ..


      M    = ML + MU + 1
      INFO = 0
c                        ** zero initial fill-in columns
      J0 = MU + 2
      J1 = MIN( N, M ) - 1

      DO 20 JZ = J0, J1

         I0 = M + 1 - JZ

         DO 10 I = I0, ML
            ABD( I, JZ ) = 0.0E0
   10    CONTINUE

   20 CONTINUE

      JZ = J1
      JU = 0
c                       ** Gaussian elimination with partial pivoting
      NM1  = N - 1

      DO 50 K = 1, NM1

         KP1 = K + 1
c                                  ** zero next fill-in column
         JZ = JZ + 1

         IF( JZ.LE.N ) THEN

            DO 30 I = 1, ML
               ABD( I, JZ ) = 0.0E0
   30       CONTINUE

         END IF
c                                  ** find L = pivot index
         LM  = MIN( ML, N - K )
         L   = ISAMAX( LM + 1, ABD( M, K ), 1 ) + M - 1
         IPVT( K ) = L + K - M

         IF( ABD( L,K ).EQ.0.0E0 ) THEN
c                                      ** zero pivot implies this column
c                                      ** already triangularized
            INFO = K

         ELSE
c                                ** interchange if necessary
            IF( L.NE.M ) THEN

               T           = ABD( L, K )
               ABD( L, K ) = ABD( M, K )
               ABD( M, K ) = T
            END IF
c                                      ** compute multipliers
            T  = - 1.0E0 / ABD( M, K )

            CALL SSCAL( LM, T, ABD( M + 1,K ), 1 )

c                               ** row elimination with column indexing

            JU = MIN( MAX( JU,MU + IPVT( K ) ), N )
            MM = M

            DO 40 J = KP1, JU

               L  = L - 1
               MM = MM - 1
               T  = ABD( L, J )

               IF( L.NE.MM ) THEN

                  ABD( L, J ) = ABD( MM, J )
                  ABD( MM, J ) = T

               END IF

               CALL SAXPY( LM, T, ABD( M+1, K ), 1, ABD( MM+1, J ), 1)

   40       CONTINUE

         END IF

   50 CONTINUE


      IPVT( N ) = N
      IF( ABD( M,N ).EQ.0.0E0 ) INFO = N

      END

      SUBROUTINE SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )

c         Solves the real band system
c            A * X = B  or  transpose(A) * X = B
c         using the factors computed by SBGCO or SGBFA.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     Input:

C        ABD     REAL(LDA, N)
c                the output from SBGCO or SGBFA.

C        LDA     INTEGER
c                the leading dimension of the array  ABD .

C        N       INTEGER
c                the order of the original matrix.

C        ML      INTEGER
c                number of diagonals below the main diagonal.

C        MU      INTEGER
c                number of diagonals above the main diagonal.

C        IPVT    INTEGER(N)
c                the pivot vector from SBGCO or SGBFA.

C        B       REAL(N)
c                the right hand side vector.

C        JOB     INTEGER
c                = 0         to solve  A*X = B ,
c                = nonzero   to solve  transpose(A)*X = B

c     On return

c        B       the solution vector  X

c     Error condition

c        A division by zero will occur if the input factor contains a
c        zero on the diagonal.  Technically, this indicates singularity,
c        but it is often caused by improper arguments or improper
c        setting of LDA .  It will not occur if the subroutines are
c        called correctly and if SBGCO has set RCOND .GT. 0.0
c        or SGBFA has set INFO .EQ. 0 .

c     To compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue

c --------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   JOB, LDA, ML, MU, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ABD( LDA, * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   K, KB, L, LA, LB, LM, M, NM1
      REAL      T
c     ..
c     .. External Functions ..

      REAL      SDOT
      EXTERNAL  SDOT
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..


      M   = MU + ML + 1
      NM1 = N - 1

      IF( JOB.EQ.0 ) THEN
c                           ** solve  A * X = B

c                               ** first solve L*Y = B
         IF( ML.NE.0 ) THEN

            DO 10 K = 1, NM1

               LM = MIN( ML, N - K )
               L  = IPVT( K )
               T  = B( L )

               IF( L.NE.K ) THEN

                  B( L ) = B( K )
                  B( K ) = T

               END IF

               CALL SAXPY( LM, T, ABD( M + 1,K ), 1, B( K + 1 ), 1 )

   10       CONTINUE

         END IF

c                           ** now solve  U*X = Y
         DO 20 KB = 1, N

            K      = N + 1 - KB
            B( K ) = B( K ) / ABD( M, K )
            LM     = MIN( K, M ) - 1
            LA     = M - LM
            LB     = K - LM
            T      = -B( K )

            CALL SAXPY( LM, T, ABD( LA,K ), 1, B( LB ), 1 )

   20    CONTINUE


      ELSE
c                          ** solve  trans(A) * X = B

c                                  ** first solve  trans(U)*Y = B
         DO 30 K = 1, N

            LM     = MIN( K, M ) - 1
            LA     = M - LM
            LB     = K - LM
            T      = SDOT( LM, ABD( LA,K ), 1, B( LB ), 1 )
            B( K ) = ( B( K ) - T ) / ABD( M, K )

   30    CONTINUE

c                                  ** now solve trans(L)*X = Y
         IF( ML.NE.0 ) THEN

            DO 40 KB = 1, NM1

               K      = N - KB
               LM     = MIN( ML, N - K )
               B( K ) = B( K ) + SDOT( LM, ABD( M+1, K ), 1,
     &                                 B( K+1 ), 1 )
               L      = IPVT( K )

               IF( L.NE.K ) THEN

                  T    = B( L )
                  B( L ) = B( K )
                  B( K ) = T

               END IF

   40       CONTINUE

         END IF

      END IF

      END

      SUBROUTINE SGECO( A, LDA, N, IPVT, RCOND, Z )

c         Factors a real matrix by Gaussian elimination
c         and estimates the condition of the matrix.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c         If  RCOND  is not needed, SGEFA is slightly faster.
c         To solve  A*X = B , follow SGECO by SGESL.

c     On entry

c        A       REAL(LDA, N)
c                the matrix to be factored.

c        LDA     INTEGER
c                the leading dimension of the array  A .

c        N       INTEGER
c                the order of the matrix  A .

c     On return

c        A       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                The factorization can be written  A = L*U , where
c                L  is a product of permutation and unit lower
c                triangular matrices and  U  is upper triangular.

c        IPVT    INTEGER(N)
c                an integer vector of pivot indices.

c        RCOND   REAL
c                an estimate of the reciprocal condition of  A .
c                For the system  A*X = B , relative perturbations
c                in  A  and  B  of size  epsilon  may cause
c                relative perturbations in  X  of size  epsilon/RCOND .
c                If  RCOND  is so small that the logical expression
c                           1.0 + RCOND .EQ. 1.0
c                is true, then  A  may be singular to working
c                precision.  In particular,  RCOND  is zero  if
c                exact singularity is detected or the estimate
c                underflows.

C        Z       REAL(N)
c                a work vector whose contents are usually unimportant.
c                If  A  is close to a singular matrix, then  Z  is
c                an approximate null vector in the sense that
c                norm(A*Z) = RCOND*norm(A)*norm(Z) .

c ------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   LDA, N
      REAL      RCOND
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      A( LDA, * ), Z( * )
c     ..
c     .. Local Scalars ..

      INTEGER   INFO, J, K, KB, KP1, L
      REAL      ANORM, EK, S, SM, T, WK, WKM, YNORM
c     ..
c     .. External Functions ..

      REAL      SASUM, SDOT
      EXTERNAL  SASUM, SDOT
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY, SGEFA, SSCAL
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MAX, SIGN
c     ..


c                        ** compute 1-norm of A
      ANORM  = 0.0E0
      DO 10 J = 1, N
         ANORM  = MAX( ANORM, SASUM( N,A( 1,J ),1 ) )
   10 CONTINUE
c                                      ** factor

      CALL SGEFA( A, LDA, N, IPVT, INFO )

c     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
c     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E .
c     trans(A) is the transpose of A.  The components of E  are
c     chosen to cause maximum local growth in the elements of W  where
c     trans(U)*W = E.  The vectors are frequently rescaled to avoid
c     overflow.

c                        ** solve trans(U)*W = E
      EK = 1.0E0

      DO 20 J = 1, N
         Z( J ) = 0.0E0
   20 CONTINUE


      DO 50 K = 1, N

         IF( Z( K ).NE.0.0E0 ) EK = SIGN( EK, -Z( K ) )

         IF( ABS( EK - Z( K ) ).GT.ABS( A( K,K ) ) ) THEN

            S  = ABS( A( K,K ) ) / ABS( EK - Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            EK = S*EK

         END IF

         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )

         IF( A( K,K ).NE.0.0E0 ) THEN

            WK   = WK / A( K, K )
            WKM  = WKM / A( K, K )

         ELSE

            WK   = 1.0E0
            WKM  = 1.0E0

         END IF

         KP1  = K + 1

         IF( KP1.LE.N ) THEN

            DO 30 J = KP1, N
               SM     = SM + ABS( Z( J ) + WKM*A( K,J ) )
               Z( J ) = Z( J ) + WK*A( K, J )
               S      = S + ABS( Z( J ) )
   30       CONTINUE

            IF( S.LT.SM ) THEN

               T  = WKM - WK
               WK = WKM

               DO 40 J = KP1, N
                  Z( J ) = Z( J ) + T*A( K, J )
   40          CONTINUE

            END IF

         END IF

         Z( K ) = WK

   50 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )
c                                ** solve trans(L)*Y = W
      DO 60 KB = 1, N
         K  = N + 1 - KB

         IF( K.LT.N )
     &       Z( K ) = Z( K ) + SDOT( N - K, A( K+1, K ), 1, Z( K+1 ), 1)

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

         END IF

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
   60 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )
c                                 ** solve L*V = Y
      YNORM  = 1.0E0

      DO 70 K = 1, N
         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T

         IF( K.LT.N ) CALL SAXPY( N - K, T, A( K + 1,K ), 1, Z( K + 1 ),
     &                            1 )

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM
         END IF

   70 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )
c                                  ** solve  U*Z = V
      YNORM  = S*YNORM

      DO 80 KB = 1, N

         K  = N + 1 - KB

         IF( ABS( Z( K ) ).GT.ABS( A( K,K ) ) ) THEN

            S  = ABS( A( K,K ) ) / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

         IF( A( K,K ).NE.0.0E0 ) Z( K ) = Z( K ) / A( K, K )

         IF( A( K,K ).EQ.0.0E0 ) Z( K ) = 1.0E0

         T  = -Z( K )

         CALL SAXPY( K - 1, T, A( 1,K ), 1, Z( 1 ), 1 )

   80 CONTINUE
c                                   ** make znorm = 1.0
      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = S*YNORM

      IF( ANORM.NE.0.0E0 ) RCOND = YNORM / ANORM
      IF( ANORM.EQ.0.0E0 ) RCOND = 0.0E0

      END

      SUBROUTINE SGEFA( A, LDA, N, IPVT, INFO )

c         Factors a real matrix by Gaussian elimination.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     SGEFA is usually called by SGECO, but it can be called
c     directly with a saving in time if  RCOND  is not needed.
c     (time for SGECO) = (1 + 9/N) * (time for SGEFA) .

c     Input:  same as SGECO

c     On return:

c        A,IPVT  same as SGECO

c        INFO    INTEGER
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  This is not an error
c                     condition for this subroutine, but it does
c                     indicate that SGESL or SGEDI will divide by zero
c                     if called.  Use  RCOND  in SGECO for a reliable
c                     indication of singularity.

c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INFO, LDA, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      A( LDA, * )
c     ..
c     .. Local Scalars ..

      INTEGER   J, K, KP1, L, NM1
      REAL      T
c     ..
c     .. External Functions ..

      INTEGER   ISAMAX
      EXTERNAL  ISAMAX
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY, SSCAL
c     ..


c                      ** Gaussian elimination with partial pivoting
      INFO = 0
      NM1  = N - 1

      DO 20 K = 1, NM1

         KP1  = K + 1
c                                            ** find L = pivot index

         L  = ISAMAX( N - K + 1, A( K,K ), 1 ) + K - 1
         IPVT( K ) = L

         IF( A( L,K ).EQ.0.0E0 ) THEN
c                                     ** zero pivot implies this column
c                                     ** already triangularized
            INFO = K

         ELSE
c                                     ** interchange if necessary
            IF( L.NE.K ) THEN

               T         = A( L, K )
               A( L, K ) = A( K, K )
               A( K, K ) = T

            END IF
c                                     ** compute multipliers
            T  = -1.0E0 / A( K, K )

            CALL SSCAL( N - K, T, A( K + 1,K ), 1 )

c                              ** row elimination with column indexing
            DO 10 J = KP1, N

               T  = A( L, J )

               IF( L.NE.K ) THEN

                  A( L, J ) = A( K, J )
                  A( K, J ) = T

               END IF

               CALL SAXPY( N-K, T, A( K+1, K ), 1, A( K+1, J ), 1 )

   10       CONTINUE

         END IF

   20 CONTINUE


      IPVT( N ) = N
      IF( A( N,N ) .EQ. 0.0E0 ) INFO = N

      END

      SUBROUTINE SGESL( A, LDA, N, IPVT, B, JOB )

c         Solves the real system
c            A * X = B  or  transpose(A) * X = B
c         using the factors computed by SGECO or SGEFA.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     On entry

c        A       REAL(LDA, N)
c                the output from SGECO or SGEFA.

c        LDA     INTEGER
c                the leading dimension of the array  A

c        N       INTEGER
c                the order of the matrix  A

c        IPVT    INTEGER(N)
c                the pivot vector from SGECO or SGEFA.

c        B       REAL(N)
c                the right hand side vector.

c        JOB     INTEGER
c                = 0         to solve  A*X = B ,
c                = nonzero   to solve  transpose(A)*X = B

c     On return

c        B       the solution vector  X

c     Error condition

c        A division by zero will occur if the input factor contains a
c        zero on the diagonal.  Technically, this indicates singularity,
c        but it is often caused by improper arguments or improper
c        setting of LDA.  It will not occur if the subroutines are
c        called correctly and if SGECO has set RCOND .GT. 0.0
c        or SGEFA has set INFO .EQ. 0 .

c     To compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue

c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   JOB, LDA, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      A( LDA, * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   K, KB, L, NM1
      REAL      T
c     ..
c     .. External Functions ..

      REAL      SDOT
      EXTERNAL  SDOT
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY
c     ..


      NM1  = N - 1

      IF( JOB.EQ.0 ) THEN
c                                 ** solve  A * X = B

c                                     ** first solve  L*Y = B
         DO 10 K = 1, NM1

            L  = IPVT( K )
            T  = B( L )

            IF( L.NE.K ) THEN

               B( L ) = B( K )
               B( K ) = T

            END IF

            CALL SAXPY( N - K, T, A( K+1, K ), 1, B( K+1 ), 1 )

   10    CONTINUE
c                                    ** now solve  U*X = Y
         DO 20 KB = 1, N

            K      = N + 1 - KB
            B( K ) = B( K ) / A( K, K )
            T      = - B( K )

            CALL SAXPY( K-1, T, A( 1, K ), 1, B(1), 1 )

   20    CONTINUE


      ELSE
c                         ** solve  trans(A) * X = B

c                                    ** first solve  trans(U)*Y = B
         DO 30 K = 1, N

            T      = SDOT( K - 1, A( 1,K ), 1, B( 1 ), 1 )
            B( K ) = ( B( K ) - T ) / A( K, K )

   30    CONTINUE

c                                    ** now solve  trans(l)*x = y
         DO 40 KB = 1, NM1

            K      = N - KB
            B( K ) = B( K ) + SDOT( N - K, A( K+1, K ), 1, B( K+1 ), 1)
            L      = IPVT( K )

            IF( L.NE.K ) THEN

               T      = B( L )
               B( L ) = B( K )
               B( K ) = T

            END IF

   40    CONTINUE

      END IF

      END

      REAL FUNCTION SASUM( N, SX, INCX )

c  INPUT--    N  Number of elements in vector to be summed
c            SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
c          INCX  Spacing of vector elements in SX

c  OUTPUT-- SASUM   Sum from 0 to N-1 of  ABS(SX(1+I*INCX))
c ----------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, N
c     ..
c     .. Array Arguments ..

      REAL      SX( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MOD
c     ..

      SASUM  = 0.0

      IF( N.LE.0 ) RETURN

      IF( INCX.NE.1 ) THEN
c                                          ** non-unit increments
         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SASUM  = SASUM + ABS( SX( I ) )
   10    CONTINUE

      ELSE
c                                          ** unit increments
         M  = MOD( N, 6 )

         IF( M.NE.0 ) THEN
c                             ** clean-up loop so remaining vector
c                             ** length is a multiple of 6.
            DO 20 I = 1, M
               SASUM  = SASUM + ABS( SX( I ) )
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 6
            SASUM  = SASUM + ABS( SX( I ) ) + ABS( SX( I + 1 ) ) +
     &               ABS( SX( I + 2 ) ) + ABS( SX( I + 3 ) ) +
     &               ABS( SX( I + 4 ) ) + ABS( SX( I + 5 ) )
   30    CONTINUE

      END IF

      END

      SUBROUTINE SAXPY( N, SA, SX, INCX, SY, INCY )

c          Y = A*X + Y  (X, Y = vectors, A = scalar)

c  INPUT--
c        N  Number of elements in input vectors X and Y
c       SA  Single precision scalar multiplier A
c       SX  Sing-prec array containing vector X
c     INCX  Spacing of elements of vector X in SX
c       SY  Sing-prec array containing vector Y
c     INCY  Spacing of elements of vector Y in SY

c OUTPUT--
c       SY   For I = 0 to N-1, overwrite  SY(LY+I*INCY) with
c                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
c            where LX = 1          if INCX .GE. 0,
c                     = (-INCX)*N  if INCX .LT. 0
c            and LY is defined analogously using INCY.
c ------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
      REAL      SA
c     ..
c     .. Array Arguments ..

      REAL      SX( * ), SY( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SY( I ) = SY( I ) + SA*SX( I )
   10    CONTINUE

      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

c                                        ** equal, unit increments
         M  = MOD( N, 4 )

         IF( M.NE.0 ) THEN
c                            ** clean-up loop so remaining vector length
c                            ** is a multiple of 4.
            DO 20 I = 1, M
               SY( I ) = SY( I ) + SA*SX( I )
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 4
            SY( I ) = SY( I ) + SA*SX( I )
            SY( I + 1 ) = SY( I + 1 ) + SA*SX( I + 1 )
            SY( I + 2 ) = SY( I + 2 ) + SA*SX( I + 2 )
            SY( I + 3 ) = SY( I + 3 ) + SA*SX( I + 3 )
   30    CONTINUE


      ELSE
c               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1
         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            SY( IY ) = SY( IY ) + SA*SX( IX )
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE

      END IF

      END

      REAL FUNCTION SDOT( N, SX, INCX, SY, INCY )

c        Single-prec dot product of vectors  X  and  Y

c  INPUT--
c        N  Number of elements in input vectors X and Y
c       SX  Sing-prec array containing vector X
c     INCX  Spacing of elements of vector X in SX
c       SY  Sing-prec array containing vector Y
c     INCY  Spacing of elements of vector Y in SY

c OUTPUT--
c     SDOT   Sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
c            where  LX = 1          if INCX .GE. 0,
c                      = (-INCX)*N  if INCX .LT. 0,
c            and LY is defined analogously using INCY.
c ------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
c     ..
c     .. Array Arguments ..

      REAL      SX( * ), SY( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      SDOT = 0.0

      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SDOT = SDOT + SX( I )*SY( I )
   10    CONTINUE


      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

c                                        ** equal, unit increments
         M  = MOD( N, 5 )

         IF( M.NE.0 ) THEN
c                            ** clean-up loop so remaining vector length
c                            ** is a multiple of 4.
            DO 20 I = 1, M
               SDOT = SDOT + SX( I )*SY( I )
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 5
            SDOT = SDOT + SX( I )*SY( I ) + SX( I + 1 )*SY( I + 1 ) +
     &               SX( I + 2 )*SY( I + 2 ) + SX( I + 3 )*SY( I + 3 ) +
     &               SX( I + 4 )*SY( I + 4 )
   30    CONTINUE

      ELSE
c               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1

         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            SDOT = SDOT + SX( IX )*SY( IY )
            IX   = IX + INCX
            IY   = IY + INCY
   40    CONTINUE

      END IF

      END

      SUBROUTINE SSCAL( N, SA, SX, INCX )

c         Multiply vector SX by scalar SA

c  INPUT--  N  Number of elements in vector
c          SA  Single precision scale factor
c          SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
c        INCX  Spacing of vector elements in SX

c OUTPUT-- SX  Replace  SX(1+I*INCX)  with  SA * SX(1+I*INCX)
c                for I = 0 to N-1
c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, N
      REAL      SA
c     ..
c     .. Array Arguments ..

      REAL      SX( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      IF( N.LE.0 ) RETURN

      IF( INCX.NE.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SX( I ) = SA*SX( I )
   10    CONTINUE


      ELSE

         M  = MOD( N, 5 )

         IF( M.NE.0 ) THEN
c                           ** clean-up loop so remaining vector length
c                           ** is a multiple of 5.
            DO 20 I = 1, M
               SX( I ) = SA*SX( I )
   20       CONTINUE

         END IF
c                             ** unroll loop for speed
         DO 30 I = M + 1, N, 5
            SX( I ) = SA*SX( I )
            SX( I + 1 ) = SA*SX( I + 1 )
            SX( I + 2 ) = SA*SX( I + 2 )
            SX( I + 3 ) = SA*SX( I + 3 )
            SX( I + 4 ) = SA*SX( I + 4 )
   30    CONTINUE

      END IF

      END

      SUBROUTINE SSWAP( N, SX, INCX, SY, INCY )

c          Interchange s.p vectors  X  and  Y, as follows:

c     For I = 0 to N-1, interchange  SX(LX+I*INCX) and SY(LY+I*INCY),
c     where LX = 1          if INCX .GE. 0,
c              = (-INCX)*N  if INCX .LT. 0
c     and LY is defined analogously using INCY.


c  INPUT--
c        N  Number of elements in input vectors X and Y
c       SX  Sing-prec array containing vector X
c     INCX  Spacing of elements of vector X in SX
c       SY  Sing-prec array containing vector Y
c     INCY  Spacing of elements of vector Y in SY

c OUTPUT--
c       SX  Input vector SY (unchanged if N .LE. 0)
c       SY  Input vector SX (unchanged IF N .LE. 0)
c --------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
c     ..
c     .. Array Arguments ..

      REAL      SX( * ), SY( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M
      REAL      STEMP1, STEMP2, STEMP3
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N-1 )*INCX, INCX
            STEMP1 = SX( I )
            SX( I ) = SY( I )
            SY( I ) = STEMP1
   10    CONTINUE


      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

c                                        ** equal, unit increments
         M  = MOD( N, 3 )

         IF( M.NE.0 ) THEN
c                            ** clean-up loop so remaining vector length
c                            ** is a multiple of 3.
            DO 20 I = 1, M
               STEMP1 = SX( I )
               SX( I ) = SY( I )
               SY( I ) = STEMP1
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 3
            STEMP1 = SX( I )
            STEMP2 = SX( I + 1 )
            STEMP3 = SX( I + 2 )
            SX( I ) = SY( I )
            SX( I + 1 ) = SY( I + 1 )
            SX( I + 2 ) = SY( I + 2 )
            SY( I ) = STEMP1
            SY( I + 1 ) = STEMP2
            SY( I + 2 ) = STEMP3
   30    CONTINUE


      ELSE
c               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1

         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            STEMP1 = SX( IX )
            SX( IX ) = SY( IY )
            SY( IY ) = STEMP1
            IX   = IX + INCX
            IY   = IY + INCY
   40    CONTINUE

      END IF

      END

      INTEGER FUNCTION ISAMAX( N, SX, INCX )

c INPUT--  N     Number of elements in vector of interest
c          SX    Sing-prec array, length 1+(N-1)*INCX, containing vector
c          INCX  Spacing of vector elements in SX

c OUTPUT-- ISAMAX   First I, I = 1 to N, to maximize
c                         ABS(SX(1+(I-1)*INCX))
c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   INCX, N
c     ..
c     .. Array Arguments ..

      REAL      SX( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, II
      REAL      SMAX, XMAG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS
c     ..


      IF( N.LE.0 ) THEN

         ISAMAX = 0

      ELSE IF( N.EQ.1 ) THEN

         ISAMAX = 1

      ELSE

         SMAX = 0.0
         II   = 1

         DO 10 I = 1, 1 + ( N-1 )*INCX, INCX

            XMAG = ABS( SX( I ) )

            IF( SMAX.LT.XMAG ) THEN

               SMAX   = XMAG
               ISAMAX = II

            END IF

            II = II + 1

   10    CONTINUE

      END IF

      END

c ---------------------------------------------------------------------
c  Fortran-90 versions of machine-constant routines R1MACH, D1MACH, I1MACH
c
c  {R,D,I}1MACH revisited: no more uncommenting DATA statements
c  
c  Presented at the IFIP WG 2.5 International Workshop on 
c  "Current Directions in Numerical Software and High Performance 
c  Computing", 19 - 20 October 1995, Kyoto, Japan. 
c  
c  The widely-used original routines were modified to use Fortran-90 
c  intrinsic functions.  This was not completely possible with I1MACH, 
c  which returns some parameters (logical unit numbers of standard
c  input, standard output, and standard error) that may require
c  user customization. 
c  
c  David Gay (dmg@bell-labs.com)
c  Eric Grosse (ehg@bell-labs.com)
c  Bell Laboratories
c  700 Mountain Avenue
c  Murray Hill, New Jersey 07974-0636
c  USA 
c  
c  References:
c  
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996.
c
c http://www.nsc.liu.se/~boein/ifip/kyoto/workshop-info/proceedings/einarsson
c    /d1mach.html  (THIS WEB SITE WORKED AS OF APR 2000)
c -------------------------------------------------------------------------


      REAL FUNCTION R1MACH (I)
c
c   R1MACH can be used to obtain machine-dependent parameters for
c   single precision numbers.  The results for various values of I are:
c
c   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c   R1MACH(3) = B**(-T), the smallest relative spacing.
c   R1MACH(4) = B**(1-T), the largest relative spacing.
c   R1MACH(5) = LOG10(B)
c
c   Assume single precision numbers are represented in the T-digit,
c   base-B form
c
c              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
c   I1MACH(10) = B, the base.
c   I1MACH(11) = T, the number of base-B digits.
c   I1MACH(12) = EMIN, the smallest exponent E.
c   I1MACH(13) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c     Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c     August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: I
      REAL :: B, X = 1.0

      B = RADIX(X)

      SELECT CASE (I)
        CASE (1)
          R1MACH = TINY(X)            ! smallest positive magnitude.
        CASE (2)
          R1MACH = HUGE(X)            ! largest magnitude.
        CASE (3)
          R1MACH = B**(-DIGITS(X))    ! smallest relative spacing.
        CASE (4)
          R1MACH = B**(1-DIGITS(X))   ! largest relative spacing.
        CASE (5)
          R1MACH = LOG10(B)
        CASE DEFAULT
          STOP 'R1MACH -- input argument out of bounds'
      END SELECT

      RETURN
      END FUNCTION R1MACH


      DOUBLE PRECISION FUNCTION D1MACH (I)
c
c   D1MACH can be used to obtain machine-dependent parameters for
c   double precision numbers.  The results for various values of I are:
c
c   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c   D1MACH(3) = B**(-T), the smallest relative spacing.
c   D1MACH(4) = B**(1-T), the largest relative spacing.
c   D1MACH(5) = LOG10(B)
c
c   Assume double precision numbers are represented in the T-digit,
c   base-B form
c
c        sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
c   I1MACH(10) = B, the base.
c   I1MACH(11) = T, the number of base-B digits.
c   I1MACH(12) = EMIN, the smallest exponent E.
c   I1MACH(13) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: I
      DOUBLE PRECISION :: B, X = 1.D0

      B = RADIX(X)

      SELECT CASE (I)
        CASE (1)
          D1MACH = TINY(X)            ! smallest positive magnitude.
        CASE (2)
          D1MACH = HUGE(X)            ! largest magnitude.
        CASE (3)
          D1MACH = B**(-DIGITS(X))    ! smallest relative spacing.
        CASE (4)
          D1MACH = B**(1-DIGITS(X))   ! largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          STOP 'D1MACH -- input arg out of bounds'
      END SELECT

      RETURN
      END FUNCTION D1MACH


      INTEGER FUNCTION I1MACH (I)
c
c   I1MACH can be used to obtain machine-dependent parameters for the
c   local machine environment.  The results for various values of I are:
c
c   I/O unit numbers (**MAY REQUIRE USER CUSTOMIZATION**):
c     I1MACH( 1) = the standard input unit.
c     I1MACH( 2) = the standard output unit.
c     I1MACH( 3) = the standard punch unit (obsolete, will cause error)
c     I1MACH( 4) = the standard error message unit.
c                  (the error message unit is usually 0 in UNIX systems)
c
c   Words:
c     I1MACH( 5) = the number of bits per integer storage unit.
c     I1MACH( 6) = the number of characters per integer storage unit.
c                  (obsolete, will cause an error)
c
c   Integers:
c     assume integers are represented in the S-digit, base-A form
c
c          sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
c
c     where 0 <= X(I) < A for I=0,...,S-1.
c
c     I1MACH( 7) = A, the base.
c     I1MACH( 8) = S, the number of base-A digits.
c     I1MACH( 9) = A**S - 1, the largest magnitude.
c
c   Floating-Point Numbers:
c     Assume floating-point numbers are represented in the T-digit,
c     base-B form
c                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c     where 0 <= X(I) .LT. B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c     I1MACH(10) = B, the base.
c
c   Single-Precision:
c     I1MACH(11) = T, the number of base-B digits.
c     I1MACH(12) = EMIN, the smallest exponent E.
c     I1MACH(13) = EMAX, the largest exponent E.
c
c   Double-Precision:
c     I1MACH(14) = T, the number of base-B digits.
c     I1MACH(15) = EMIN, the smallest exponent E.
c     I1MACH(16) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   750101  DATE WRITTEN
c   960411  Modified for Fortran 90 (BE after suggestions by Eric Grosse)    
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: I
      REAL :: X_single  = 1.0
      DOUBLE PRECISION :: X_double = 1.D0

      SELECT CASE (I)
        CASE (1)
          I1MACH = 5 ! Input unit
        CASE (2)
          I1MACH = 6 ! Output unit
        CASE (3)
          STOP 'I1MACH: input arg = 3 is obsolete'
        CASE (4)
          I1MACH = 0 ! Error message unit
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          STOP 'I1MACH: input arg = 6 is obsolete'
        CASE (7)
          I1MACH = RADIX(1)
        CASE (8)
          I1MACH = BIT_SIZE(I) - 1
        CASE (9)
          I1MACH = HUGE(1)
        CASE (10)
          I1MACH = RADIX(X_single)
        CASE (11)
          I1MACH = DIGITS(X_single)
        CASE (12)
          I1MACH = MINEXPONENT(X_single)
        CASE (13)
          I1MACH = MAXEXPONENT(X_single)
        CASE (14)
          I1MACH = DIGITS(X_double)
        CASE (15)
          I1MACH = MINEXPONENT(X_double)
        CASE (16)
          I1MACH = MAXEXPONENT(X_double) 
        CASE DEFAULT
          STOP 'I1MACH: input argument out of bounds'
      END SELECT

      RETURN
      END FUNCTION I1MACH


c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: DISOTEST.f,v 2.1 2000/04/03 21:21:55 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE  disort_driver(NLYR, DTAUC, SSALB, NMOM,
     &                   NCLDLYR, CLDLYR, PMOM_CLD, TEMPER,
     &                   WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                   USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                   UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                   TTEMP, TEMIS, PLANK, ONLYFL, HEADER,
     &                   MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM,
     &                   RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, 
     &                   TRNMED, msg, errflag)


c    Runs DISORT using inputs from Python.
c
c    Based on DISOTEST codes but with modifications by
c    Penny M. Rowe and Steven Neshyba

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
c
c


      INTEGER  MAXCLY, MAXMOM, MAXPHI, MAXULV, MAXUMU, CLDLYR(*),
     &         NCLDLYR, I, J, ICLDLYR
      CHARACTER  HEADER*127, msg*127
      LOGICAL  LAMBER, PLANK, ONLYFL, PRNT(5), USRANG, USRTAU,
     &         errflag
      INTEGER  IBCND, NMOM, NLYR, NUMU, NSTR, NPHI, NTAU
      REAL     ACCUR, ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,
     &         PHI( MAXPHI ), PMOM( 0:MAXMOM, MAXCLY ), 
     &         PMOM_CLD(0:MAXMOM,*),
     &         PHI0, SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,
     &         WVNMLO, WVNMHI, UMU( MAXUMU ), UMU0, UTAU( MAXULV )

      REAL     RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),
     &         DFDT( MAXULV ), UAVG( MAXULV ),
     &         UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),
     &         TRNMED( MAXUMU )


c     .. External Subroutines ..

      EXTERNAL  DISORT, ERRMSG, GETMOM, PRTFIN
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ASIN, FLOAT, INDEX

      DATA  PRNT / .FALSE., 3*.FALSE., .FALSE. /, ACCUR / 0.0 /



c     Lines added for f2py:
cf2py intent(in)  :: NLYR, DTAUC, SSALB, NMOM,
cf2py intent(in)  :: NCLDLYR,CLDLYR, PMOM_CLD,
cf2py intent(in)  :: TEMPER,WVNMLO,WVNMHI, USRTAU, NTAU, UTAU, NSTR,
cf2py intent(in)  :: USRANG,NUMU,UMU, NPHI, PHI, IBCND, FBEAM,
cf2py intent(in)  :: UMU0, PHI0,FISOT, LAMBER, ALBEDO, BTEMP,
cf2py intent(in)  :: TTEMP, TEMIS,PLANK, ONLYFL,HEADER,
cf2py intent(in)  :: MAXCLY,MAXULV, MAXUMU, MAXPHI, MAXMOM

cf2py intent(out) :: RFLDIR,RFLDN,FLUP,DFDT,UAVG,UU,ALBMED,TRNMED
cf2py intent(out) :: msg, errflag




	 
c    .. Initialize PMOM matrix to zeros
      DO I=0,MAXMOM
        DO J=1,MAXCLY
          PMOM(I,J) = 0.0
        ENDDO
      ENDDO


c   ..Fill in the moments from the namelist
      DO ICLDLYR=1,NCLDLYR
        DO I=0,NMOM
           PMOM(I, CLDLYR(ICLDLYR)) = PMOM_CLD(I,ICLDLYR)
        ENDDO
      ENDDO

c   ..** WARNING **  If UMU0 equals one of the
c             computational polar angle cosines, a singularity
c             occurs;  hence this is treated as a fatal
c             error. The problem is most likely to
c             occur when NSTR/2 is odd and UMU0 = 0.5;
c             otherwise, it is almost impossible to hit a
c             computational angle by chance.  The problem can
c             easily be corrected by changing NSTR.
c
c   ..Here we check for cases where this occurred and exclude them
c       IF (round(UMU0,5)==0.23724) and (NSTR==16) and (NMOM>=19): NSTR = 18




c   ..Don't let number streams exceed number of legendre moments
c     (this has to be handled BEFORE calling this code)
c      IF (NSTR .GT. NMOM) THEN
c		      NSTR = NMOM
c      ENDIF

c    .. UTAU specifies the optical depth (layer height) at which
c    .. we want values reported. This is often either 0 (TOA)
c    .. or ~SUM(DTAUC) (the surface). However,
c    .. sometimes UTAU is bigger than sum(DTAUC) due to
c    .. differences in how single precision numbers are computed.
c    .. Don't allow this, instead reset to sum(DTAUC).
c    .. P. Rowe, 2015/08/18
      DO I=1,NTAU
        IF (UTAU(I) .GT. SUM(DTAUC)) THEN
          UTAU(I) = SUM(DTAUC)
        ENDIF
      ENDDO



c
c   ..Call disort..
      CALL  DISORT(
     &                 NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED, msg, errflag )


      RETURN
      END





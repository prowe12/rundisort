c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: DISOTEST.f,v 2.1 2000/04/03 21:21:55 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE  disort_driver(NLYR, NMOM, NSTR, NUMU, NPHI, NTAU, 
     &                   USRANG, USRTAU, IBCND, ONLYFL, 
     &                   PLANK, LAMBER,  
     &                   DTAUC, SSALB, PMOM_CLD, TEMPER, WV NMLO, WVNMHI,
     &                   UTAU, UMU0, PHI0, UMU, PHI, FBEAM,
     &                   FISOT, ALBEDO, BTEMP, TTEMP, TEMIS, 
     &                   HEADER,
     &                   RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, 
     &                   ALBMED,TRNMED,     
     &                   NCLDLYR, CLDLYR, msg, errflag)


c    Based on DISOTEST codes but with modifications by
c    Penny M. Rowe and Steven Neshyba
c
C  DISORT variables not in call, which are set below:
c     PRNT
c     DELTAMPLUS, DO_PSEUDO_SPHERE,
c     PMOM (instead PMOM_CLD)
c     EARTH_RADIUS, H_LYR
c     RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,
c     ACCUR
c
c  Extra variables in call (not in DISORT call)
c     &                   NCLDLYR, CLDLYR, 
c     &                   msg, errflag)
c
c  Routines called :

c    DISORT:   The discrete ordinates radiative transfer program

c+---------------------------------------------------------------------+
c
c                 ** DISORT I/O specifications **
c
c
c
c Note: the following assignments are now made in DISORT.f:
C     NLYR = MAXCLY
C     NMOM = MAXMOM
C     NSTR = MAXCMU   (not used here)
C     NUMU = MAXUMU
C     NPHI = MAXPHI
C     NTAU = MAXULV



      INTEGER  NCLDLYR, ICLDLYR
      INTEGER  NLYR, NMOM, NPHI, NTAU, NUMU, CLDLYR(NCLDLYR),
     &         I, J, 
      CHARACTER  HEADER*127, msg*127
      LOGICAL  LAMBER, PLANK, ONLYFL, PRNT(5), USRANG, USRTAU,
     &         errflag,
     &         DO_PSEUDO_SPHERE, DELTAMPLUS
      INTEGER  IBCND, NSTR
      REAL     ACCUR, ALBEDO, BTEMP, DTAUC( NLYR ), FBEAM, FISOT,
     &         PHI( NPHI ), PMOM( 0:NMOM, NLYR ), 
     &         PMOM_CLD(0:NMOM,*),
     &         PHI0, SSALB( NLYR ), TEMPER( 0:NLYR ), TEMIS, TTEMP,
     &         WVNMLO, WVNMHI, UMU( NUMU ), UMU0, UTAU( NTAU )

      REAL     RFLDIR( NTAU ), RFLDN( NTAU ), FLUP( NTAU ),
     &         DFDT( NTAU ), UAVG( NTAU ),
     &         UU( NUMU, NTAU, NPHI ), ALBMED( NUMU ),
     &         TRNMED( NUMU )

c     .. New stuff, might not work
      REAL(kind=4),parameter :: EARTH_RADIUS = 6371.0 
      
      REAL(kind=4) :: RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), 
     &         RHOU(NUMU, 0:NSTR/2, 0:(NSTR-1)), &
     &         EMUST(NUMU), BEMST(NSTR/2), RHO_ACCURATE(NUMU, NPHI)          

      REAL(kind=4) :: DFDT( NTAU )




c     .. External Subroutines ..

      EXTERNAL  DISORT, ERRMSG, GETMOM, PRTFIN
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ASIN, FLOAT, INDEX

      DATA  PRNT / 5*.FALSE. /, ACCUR / 0.0 /

c     .. New unfamiliar variables
      DO_PSEUDO_SPHERE = .FALSE.
      DELTAMPLUS = .FALSE.
      H_LYR = 0.0; RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; 
      BEMST = 0.0; RHO_ACCURATE = 0.0; DFDT = 0.0; 


c     Lines added for f2py:
cf2py intent(in)  :: NLYR, NMOM, NSTR, NUMU, NPHI, NTAU, 
cf2py intent(in)  :: DTAUC, SSALB, 
cf2py intent(in)  :: NCLDLYR,CLDLYR, PMOM_CLD,
cf2py intent(in)  :: TEMPER, WVNMLO, WVNMHI, USRTAU, UTAU,
cf2py intent(in)  :: USRANG, UMU, PHI, IBCND, FBEAM,
cf2py intent(in)  :: UMU0, PHI0,FISOT, LAMBER, ALBEDO, BTEMP,
cf2py intent(in)  :: TTEMP, TEMIS,PLANK, ONLYFL,HEADER,

cf2py intent(out) :: RFLDIR,RFLDN,FLUP,DFDT,UAVG,UU,ALBMED,TRNMED
cf2py intent(out) :: msg, errflag




	 
c    .. Initialize PMOM matrix to zeros
      DO I=0,NMOM
        DO J=1,NLYR
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
c      CALL  DISORT(
c     &                 NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
c     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
c     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
c     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
c     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
c     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
c     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
c     &                 ALBMED, TRNMED, msg, errflag )


C     CALL DISORT( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU,           &
C    &             USRANG, USRTAU, IBCND, ONLYFL, PRNT,          &
C    &             PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,  &          
C    &             DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,   & 
C    &             UTAU, UMU0, PHI0, UMU, PHI, FBEAM,            &                        
C    &             FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,           &
C    &             EARTH_RADIUS, H_LYR,                          &
C    &             RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,       &
C    &             ACCUR,  HEADER,                               &
C    &             RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,          &
C    &             ALBMED, TRNMED )


      RETURN
      END





c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c $Rev: 90 $ $Date: 2017-11-30 20:01:24 -0500 (Thu, 30 Nov 2017) $
c FORTRAN 77
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE disort_driver( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU,
     &                          USRANG, USRTAU, IBCND, ONLYFL,
     &                          PLANK, LAMBER, DTAUC, SSALB, 
     &                          NCLDLYR, CLDLYR, PMOM_CLD, TEMPER, 
     &                          WVNMLO, WVNMHI,
     &                          UTAU, UMU0, PHI0, UMU, PHI, FBEAM,
     &                          FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,
     &                          HEADER, 
     &                          MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM,
     &                          RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                          ALBMED, TRNMED, 
     &                          msg, errflag) 

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


      INTEGER  NLYR, NMOM, NSTR, NUMU, NPHI, NTAU
      INTEGER  MAXCLY, MAXMOM, MAXPHI, MAXULV, MAXUMU

c     ..
c     .. Scalar Arguments ..
      CHARACTER HEADER*127
      LOGICAL   LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO

c     ..
c     .. Array Arguments ..
      LOGICAL   PRNT( 5 )
c      REAL      ALBMED( NUMU ), DFDT( NTAU ), DTAUC( NLYR ),
c     &          FLUP( NTAU ), PHI( NPHI ), 
c     &          PMOM_CLD( 0:NMOM, NCLDLYR ), PMOM( 0:NMOM, NLYR ),
c     &          RFLDIR( NTAU ), RFLDN( NTAU ), SSALB( NLYR ),
c     &          TEMPER( 0:NLYR ), TRNMED( NUMU ), UAVG( NTAU ),
c     &          UMU( NUMU ), UTAU( NTAU ),
c     &          UU( NUMU, NTAU, NPHI )

      REAL     DTAUC( MAXCLY ), 
     &         PHI( MAXPHI ), PMOM( 0:MAXMOM, MAXCLY ), 
     &         PMOM_CLD(0:MAXMOM,*),
     &         SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), 
     &         UMU( MAXUMU ), UTAU( MAXULV )

      REAL     RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),
     &         DFDT( MAXULV ), UAVG( MAXULV ),
     &         UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),
     &         TRNMED( MAXUMU )



c     ..
c     .. Version 3 .. 
      REAL      RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), 
     &          RHOU(NSTR,   0:NSTR/2, 0:(NSTR-1)),
     &          EMUST(MAXUMU), BEMST(NSTR/2)
      REAL      RHO_ACCURATE(MAXUMU, MAXPHI)  

c     .. Version 3: spherical correction ..
      REAL(kind=4),parameter :: EARTH_RADIUS = 6371.0     
      REAL      H_LYR(0:MAXCLY)
      LOGICAL   DO_PSEUDO_SPHERE

c     .. Version 3: deltam plus
      LOGICAL   DELTAMPLUS

c     .. Extra variables for running from Python
      CHARACTER msg*127
      INTEGER   NCLDLYR, ICLDLYR, CLDLYR( * )
      LOGICAL   errflag
c      INTEGER   NCLDLYR, ICLDLYR, CLDLYR(NCLDLYR)
      

      DATA      PRNT  / .FALSE., 3*.FALSE., .FALSE. /



c     Lines added for f2py:
cf2py intent(in)  :: NLYR, NMOM, NSTR, NUMU, NPHI, NTAU
cf2py intent(in)  :: USRANG, USRTAU, IBCND, ONLYFL
cf2py intent(in)  :: PLANK, LAMBER, DTAUC, SSALB
cf2py intent(in)  :: NCLDLYR, CLDLYR, PMOM_CLD, TEMPER
cf2py intent(in)  :: WVNMLO, WVNMHI
cf2py intent(in)  :: UTAU, UMU0, PHI0, UMU, PHI, FBEAM
cf2py intent(in)  :: FISOT, ALBEDO, BTEMP, TTEMP, TEMIS
cf2py intent(in)  :: HEADER
cf2py intent(in)  :: MAXCLY,MAXULV, MAXUMU, MAXPHI, MAXMOM

cf2py intent(out) :: RFLDIR,RFLDN,FLUP,DFDT,UAVG,UU
cf2py intent(out) :: ALBMED,TRNMED
cf2py intent(out) :: msg, errflag

c     I am not sure what these do. Set to false and 0.0 as in disotest      
      DO_PSEUDO_SPHERE = .FALSE.
      DELTAMPLUS = .FALSE.
      H_LYR = 0.0
      RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; BEMST = 0.0; 
      RHO_ACCURATE = 0.0; ACCUR = 0.0

c    .. Initialize H_LYR
      DO I=0,NLYR
        H_LYR(I) = 0.0
      ENDDO

c    .. Initialize PMOM matrix to zeros
      DO I=0,NMOM
        DO J=1,NLYR
          PMOM(I,J) = 0.0
        ENDDO
      ENDDO

c   ..Fill in the moments from PMOM_CLD
      DO ICLDLYR=1,NCLDLYR
        DO I=0,NMOM
           PMOM(I, CLDLYR(ICLDLYR)) = PMOM_CLD(I,ICLDLYR)
        ENDDO
      ENDDO

c    .. UTAU specifies the optical depth (layer height) at which
c    .. we want values reported. This is often either 0 (TOA)
c    .. or ~SUM(DTAUC) (the surface). However,
c    .. sometimes UTAU is bigger than sum(DTAUC) due to
c    .. differences in how single precision numbers are computed.
c    .. Don't allow this, instead reset to sum(DTAUC).
c    .. P. Rowe, 2015/08/18
c      DO I=1,NTAU
c        IF (UTAU(I) .GT. SUM(DTAUC)) THEN
c          UTAU(I) = SUM(DTAUC)
c        ENDIF
c      ENDDO
      
      errflag = .FALSE.
      msg = ' '


      CALL DISORT(
     &             NLYR, NMOM, NSTR, NUMU, NPHI, NTAU,           
     &             USRANG, USRTAU, IBCND, ONLYFL, PRNT,          
     &             PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,            
     &             DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,    
     &             UTAU, UMU0, PHI0, UMU, PHI, FBEAM,                                    
     &             FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,          
     &             EARTH_RADIUS, H_LYR,                         
     &             RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,       
     &             ACCUR,  HEADER,                               
     &             RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,          
     &             ALBMED, TRNMED )


      RETURN
      END
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -





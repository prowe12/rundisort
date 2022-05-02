c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c $Rev: 90 $ $Date: 2017-11-30 20:01:24 -0500 (Thu, 30 Nov 2017) $
c FORTRAN 77
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE disort_driver( NLYR, NMOM, NSTR, 
     &                   NUMU, NPHI, NTAU,
     &                   USRANG, USRTAU, IBCND, ONLYFL,
     &                   PLANK, LAMBER, 
     &                   DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,
     &                   UTAU, UMU0, PHI0, UMU, PHI, FBEAM,
     &                   FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,
     &                   HEADER,
     &                   RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                   ALBMED, TRNMED, 
     &                   NCLDLYR, CLDLYR, msg, errflag) 

c      USE PARAMETERS 
      INTEGER   NLYR, NMOM, NPHI, NTAU, NUMU, NSTR

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
      REAL      ALBMED( NUMU ), DFDT( NTAU ), DTAUC( NLYR ),
     &          FLUP( NTAU ), PHI( NPHI ), 
     &          PMOM( 0:NMOM, NLYR ),
     &          RFLDIR( NTAU ), RFLDN( NTAU ), SSALB( NLYR ),
     &          TEMPER( 0:NLYR ), TRNMED( NUMU ), UAVG( NTAU ),
     &          UMU( NUMU ), UTAU( NTAU ),
     &          UU( NUMU, NTAU, NPHI )

c     ..
c     .. Version 3 .. 
      REAL      RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), 
     &          RHOU(NUMU,   0:NSTR/2, 0:(NSTR-1)),
     &          EMUST(NUMU), BEMST(NSTR/2)
      REAL      UMU0DI, UMU0SQ, DENOM
      REAL      RHO_ACCURATE(NUMU,NPHI)  
      REAL      EIGEN_MAT(NSTR/2, NSTR/2 )  

c     .. Version 3: spherical correction ..
      REAL(kind=4),parameter :: EARTH_RADIUS = 6371.0     
      REAL      H_LYR(0:NLYR)
      REAL      UMU0L(NLYR)
      LOGICAL   DO_PSEUDO_SPHERE

c     .. Version 3: deltam plus
      LOGICAL   DELTAMPLUS

c     .. Extra variables for running from Python
      CHARACTER msg*127
      INTEGER   NCLDLYR, ICLDLYR, CLDLYR(NCLDLYR)
      LOGICAL   errflag
      

      DATA      PRNT  / .FALSE., 3*.FALSE., .FALSE. /

c Note: the following assignments are now made in DISORT.f:
C     NLYR = MAXCLY
C     NMOM = MAXMOM
C     NSTR = MAXCMU
C     NUMU = MAXUMU
C     NPHI = MAXPHI
C     NTAU = MAXULV


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

c     I am not sure what these do. Set to false and 0.0 as in disotest      
      DO_PSEUDO_SPHERE = .FALSE.
      DELTAMPLUS = .FALSE.
      H_LYR(:) = 0.0
      RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; BEMST = 0.0; RHO_ACCURATE = 0.0;
      ACCUR = 0.0


      CALL DISORT( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU,           
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





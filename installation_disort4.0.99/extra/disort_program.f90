

program DISORT_PROGRAM
implicit none

real(kind=4),parameter :: PI = 2.*ASIN(1.0)
real(kind=4),parameter :: EARTH_RADIUS = 6371.0 
integer, parameter :: NLYR = 6, NTAU = 5, NUMU = 4
integer, parameter :: NSTR = 8, NMOM = 8, NPHI = 1
integer, parameter :: NCLDLYR = 2

INTEGER ::  K, LC
INTEGER :: IBCND
INTEGER :: CLDLYR(NCLDLYR)
CHARACTER  HEADER*127, msg*127

LOGICAL  PRNT(5), errflag
DATA PRNT  / .TRUE., 3*.FALSE., .TRUE. /

LOGICAL  USRANG, USRTAU, ONLYFL, &
         PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE
         
real(kind=4) :: ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, &
                PHI0, TEMIS, TTEMP, WVNMLO, WVNMHI, UMU0 

real(kind=4) :: ALBMED(NUMU)
real(kind=4) :: SSALB( NLYR ), PMOM_CLD( 0:NMOM, NCLDLYR ), &
                TEMPER( 0:NLYR ), UTAU( NTAU ), UMU( NUMU ), &
                PHI( NPHI ), H_LYR( 0:NLYR ) 
real(kind=4) :: RFLDIR( NTAU ), RFLDN( NTAU ), FLUP( NTAU ), UAVG( NTAU ), &
                TRNMED( NUMU ), UU( NUMU, NTAU, NPHI )   
real(kind=4) :: RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), RHOU(NUMU, 0:NSTR/2, 0:(NSTR-1)), &
                EMUST(NUMU), BEMST(NSTR/2), RHO_ACCURATE(NUMU, NPHI)          
real(kind=4) :: DTAUC( NLYR )
real(kind=4) :: DFDT( NTAU )


ACCUR = 0.0

USRTAU    = .TRUE.
USRANG    = .TRUE.
LAMBER    = .TRUE.
PLANK     = .FALSE.
ONLYFL    = .FALSE.
DO_PSEUDO_SPHERE = .FALSE.
DELTAMPLUS = .FALSE.
         
IBCND  = 0     
TEMPER = 0.0; 
RFLDIR = 0.0; RFLDN = 0.0; FLUP = 0.0; UAVG = 0.0; UU = 0.0;
TRNMED = 0.0; 
ALBMED = 0.0; 

CLDLYR (1) = 4; CLDLYR (2) = 5
UMU( 1 )  = -1.0; UMU( 2 )  = -0.2; UMU( 3 ) = 0.2; UMU( 4 ) = 1.0;
PHI( 1 )  = 60.0; 
UTAU( 1 ) = 0.0;  UTAU( 2 ) = 1.05; UTAU( 3 ) = 2.1; UTAU( 4 ) = 6.0; UTAU( 5 ) = 21.0; 

DO LC = 1, NLYR
  DTAUC( LC ) = LC/2
  SSALB( LC ) = 0
ENDDO

DO LC = 1, NCLDLYR
  PMOM_CLD(0,LC) = 1.0
  DO  K = 1, NMOM
    PMOM_CLD(K,LC) = 0.0
  ENDDO
  DTAUC( LC ) = LC
  SSALB( LC ) = 0.6 + LC*0.05
ENDDO

! New unfamiliar variables
H_LYR = 0.0; RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; 
BEMST = 0.0; RHO_ACCURATE = 0.0; DFDT = 0.0; 


FBEAM      = 0.0;   
FISOT      = 1.0/PI
PLANK      = .FALSE.;     
LAMBER     = .TRUE.; 
ALBEDO   = 0.0; 
HEADER = ''

! CALL DISORT( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU,           &
!                  USRANG, USRTAU, IBCND, ONLYFL, PRNT,          &
!                  PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,  &          
!                  DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,   & 
!                  UTAU, UMU0, PHI0, UMU, PHI, FBEAM,            &                        
!                  FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,           &
!                  EARTH_RADIUS, H_LYR,                          &
!                  RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,       &
!                  ACCUR,  HEADER,                               &
!                  RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,          &
!                  ALBMED, TRNMED )


! CALL DISORT( MAXCLY, MAXMOM, MAXCMU, 
!      &                   MAXUMU, MAXPHI, MAXULV,
!      &                   USRANG, USRTAU, IBCND, ONLYFL, PRNT,
!      &                   PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,
!      &                   DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,
!      &                   UTAU, UMU0, PHI0, UMU, PHI, FBEAM,
!      &                   FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,
!      &                   EARTH_RADIUS, H_LYR, 
!      &                   RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,
!      &                   ACCUR,  HEADER,
!      &                   RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
!      &                   ALBMED, TRNMED )    
     
! Note: the following assignments are now made in DISORT.f:
!     NLYR = MAXCLY
!     NMOM = MAXMOM
!     NSTR = MAXCMU
!     NUMU = MAXUMU
!     NPHI = MAXPHI
!     NTAU = MAXULV
     
! CALL disort_driver(NLYR, NMOM, NSTR, NUMU, NPHI, NTAU, &
!                     USRANG, USRTAU, IBCND, ONLYFL, &
!                     PLANK, LAMBER,  &
!                     DTAUC, SSALB, PMOM_CLD, TEMPER, WVNMLO, WVNMHI,&
!                     UTAU, UMU0, PHI0, UMU, PHI, FBEAM,&
!                     FISOT, ALBEDO, BTEMP, TTEMP, TEMIS, &
!                     HEADER,&
!                     RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, &
!                     ALBMED,TRNMED,     &
!                     NCLDLYR, CLDLYR, msg, errflag)                       

WRITE(*,*) "Complete again"

end program DISORT_PROGRAM


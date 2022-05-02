! c **********************************************************************
module disort_vars
implicit none 

LOGICAL DOPROB( 17 )
DATA DOPROB / 17*.TRUE. /
LOGICAL  USRANG, USRTAU, ONLYFL, PRNT(5), &
         PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE
INTEGER  IBCND, NMOM, NLYR, NUMU, NSTR, NPHI, NTAU
INTEGER  NUMU_O
LOGICAL  DEBUG
real(kind=4) :: ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, &
                PHI0, TEMIS, TTEMP, WVNMLO, WVNMHI, UMU0 
real(kind=4),parameter :: EARTH_RADIUS = 6371.0     
DATA PRNT  / .TRUE., 3*.FALSE., .TRUE. /

INTEGER       BRDF_TYPE
REAL          BRDF_ARG(4)
LOGICAL       DO_SHADOW
REAL          WIND_SPD, REFRAC_INDEX
REAL          B0, HH, W
REAL          K_VOL, K_ISO, K_GEO
REAL          RHO_0, KAPPA, G, H0
REAL          FLUX_UP, DFDTAU
INTEGER       NMUG
REAL          BDREF
EXTERNAL      BDREF

real(kind=4),dimension(:),allocatable     :: DTAUC, PHI, SSALB, TEMPER, UMU, UTAU                             
real(kind=4),dimension(:,:),allocatable   :: PMOM_cld          
real(kind=4),dimension(:,:,:),allocatable :: RHOQ, RHOU 
real(kind=4),dimension(:),allocatable     :: EMUST, BEMST   
real(kind=4),dimension(:,:),allocatable   :: RHO_ACCURATE                             
real(kind=4),dimension(:),allocatable     :: RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED
real(kind=4),dimension(:,:,:),allocatable :: UU
real(kind=4),dimension(:),allocatable     :: H_LYR

contains 

subroutine allocate_disort_allocatable_arrays(NLYR, NMOM, NSTR, NUMU, NPHI, NTAU)
implicit none
integer,intent(in) :: NLYR, NMOM, NSTR, NUMU, NPHI, NTAU 

allocate( DTAUC( NLYR ), SSALB( NLYR ), PMOM_cld( 0:NMOM, NLYR ), &
          TEMPER( 0:NLYR ), UTAU( NTAU ), UMU( NUMU ), PHI( NPHI ), H_LYR( 0:NLYR ) )  
allocate( RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), RHOU(NUMU, 0:NSTR/2, 0:(NSTR-1)), &
          EMUST(NUMU), BEMST(NSTR/2), RHO_ACCURATE(NUMU, NPHI) )                
allocate( RFLDIR( NTAU ), RFLDN( NTAU ), FLUP( NTAU ), DFDT( NTAU ), UAVG( NTAU ),&
          ALBMED( NUMU ), TRNMED( NUMU ), UU( NUMU, NTAU, NPHI ) )   
DTAUC = 0.0; SSALB = 0.0; PMOM_cld = 0.0; TEMPER = 0.0; UTAU = 0.0; UMU = 0.0; PHI = 0.0;
H_LYR = 0.0; RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; BEMST = 0.0; RHO_ACCURATE = 0.0;
RFLDIR = 0.0; RFLDN = 0.0; FLUP = 0.0; DFDT = 0.0; UAVG = 0.0; UU = 0.0;
ALBMED = 0.0; TRNMED = 0.0; 

end subroutine allocate_disort_allocatable_arrays

subroutine deallocate_disort_allocatable_arrays()

deallocate( DTAUC, SSALB, PMOM_cld, TEMPER, UTAU, UMU, PHI, H_LYR )  
deallocate( RHOQ, RHOU, EMUST, BEMST, RHO_ACCURATE )                
deallocate( RFLDIR, RFLDN, FLUP, DFDT, UAVG, ALBMED, TRNMED, UU )  

end subroutine deallocate_disort_allocatable_arrays

end module disort_vars

! c **********************************************************************






program DISORT_PROGRAM
use disort_vars
implicit none


real(kind=4),parameter :: PI = 2.*ASIN(1.0)
INTEGER :: IOD, IU, I, J, K, LC, LENTIT, LU, NPROB
CHARACTER  HEADER*127

! Extra variables for running from Python
CHARACTER msg*127
INTEGER,parameter :: NCLDLYR = 2
INTEGER CLDLYR(NCLDLYR)
LOGICAL errflag

ACCUR = 0.0

USRTAU    = .TRUE.
USRANG    = .TRUE.
LAMBER    = .TRUE.
PLANK     = .FALSE.
ONLYFL    = .FALSE.
DO_PSEUDO_SPHERE = .FALSE.
DELTAMPLUS = .FALSE.

NSTR = 8; IF(MOD(NSTR,2).NE.0) NSTR = NSTR+1;
NLYR = 6; 
NMOM = NSTR 
NTAU      = 5; IF(.NOT.USRTAU) NTAU = NLYR + 1
NUMU      = 4; IF(.NOT.USRANG) NUMU = NSTR
NPHI      = 1; 
IBCND     = 0     
         
call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )

!UMU0      = 0.5;  PHI0      = 0.0
UMU( 1 )  = -1.0; UMU( 2 )  = -0.2; UMU( 3 ) = 0.2; UMU( 4 ) = 1.0;
PHI( 1 )  = 60.0; 
UTAU( 1 ) = 0.0;  UTAU( 2 ) = 1.05; UTAU( 3 ) = 2.1; UTAU( 4 ) = 6.0; UTAU( 5 ) = 21.0; 

CLDLYR(1)=4; CLDLYR(2)=5
DO LC = 1, NLYR
  DTAUC( LC ) = LC
  SSALB( LC ) = 0
ENDDO

DO LC = 1, NCLDLYR
  PMOM_cld(0,LC) = 1.0
  DO  K = 1, NMOM
    PMOM_cld(K,LC) = 0.0
  SSALB( CLDLYR(LC) ) = 0.6 + LC*0.05
  ENDDO
ENDDO


FBEAM      = 0.0;   FISOT      = 1.0/PI
PLANK      = .FALSE.;     
LAMBER     = .TRUE.; ALBEDO   = 0.0; 

HEADER = ''

CALL disort_driver( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU,           &
                 USRANG, USRTAU, IBCND, ONLYFL,                &
                 PLANK, LAMBER,  DTAUC, SSALB, &          
                 NCLDLYR, CLDLYR, PMOM_cld, &
                 TEMPER, WVNMLO, WVNMHI,   & 
                 UTAU, UMU0, PHI0, UMU, PHI, FBEAM,              &                        
                 FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,             &
                 HEADER,                         &
                 RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,            &
                 ALBMED, TRNMED,          &
                 msg, errflag)                       
                 
!                     DTAUC, SSALB, PMOM_CLD, TEMPER, WVNMLO, WVNMHI,&
!                     RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, &
!                     ALBMED,TRNMED,     &
!                     NCLDLYR, CLDLYR, msg, errflag)                       

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
                           
                    
                   
call deallocate_disort_allocatable_arrays()



WRITE(*,*) "Still Success"

end program DISORT_PROGRAM


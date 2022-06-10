! c **********************************************************************
module disort_vars
implicit none 

real(kind=4),dimension(:,:,:),allocatable :: RHOQ, RHOU 
real(kind=4),dimension(:),allocatable     :: EMUST, BEMST   
real(kind=4),dimension(:,:),allocatable   :: RHO_ACCURATE                             
real(kind=4),dimension(:),allocatable     :: PHI, SSALB, TEMPER, UMU, UTAU                             
real(kind=4),dimension(:,:),allocatable   :: PMOM          
real(kind=4),dimension(:),allocatable     :: RFLDIR, RFLDN, FLUP, UAVG, TRNMED
real(kind=4),dimension(:,:,:),allocatable :: UU
real(kind=4),dimension(:),allocatable     :: H_LYR
!real(kind=4),dimension(:),allocatable     :: ALBMED
real(kind=4),dimension(:),allocatable     :: DTAUC                             
real(kind=4),dimension(:),allocatable     :: DFDT

contains 

subroutine allocate_disort_allocatable_arrays(NLYR, NMOM, NSTR, NUMU, NPHI, NTAU)
implicit none
integer,intent(in) :: NLYR, NMOM, NSTR, NUMU, NPHI, NTAU 

allocate( SSALB( NLYR ), PMOM( 0:NMOM, NLYR ), &
          TEMPER( 0:NLYR ), UTAU( NTAU ), UMU( NUMU ), PHI( NPHI ), H_LYR( 0:NLYR ) )  
allocate( RFLDIR( NTAU ), RFLDN( NTAU ), FLUP( NTAU ), UAVG( NTAU ),&
          TRNMED( NUMU ), UU( NUMU, NTAU, NPHI ) )   
allocate( RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), RHOU(NUMU, 0:NSTR/2, 0:(NSTR-1)), &
          EMUST(NUMU), BEMST(NSTR/2), RHO_ACCURATE(NUMU, NPHI) )                
allocate( DTAUC( NLYR ))
!allocate( ALBMED( NUMU ))
allocate( DFDT( NTAU ))

SSALB = 0.0; PMOM = 0.0; TEMPER = 0.0; UTAU = 0.0; UMU = 0.0; PHI = 0.0;
RFLDIR = 0.0; RFLDN = 0.0; FLUP = 0.0; UAVG = 0.0; UU = 0.0;
TRNMED = 0.0; 
H_LYR = 0.0; RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; BEMST = 0.0; RHO_ACCURATE = 0.0;
DTAUC = 0.0; DFDT = 0.0; 
!ALBMED = 0.0; 

end subroutine allocate_disort_allocatable_arrays

subroutine deallocate_disort_allocatable_arrays()

deallocate( DTAUC, SSALB, PMOM, TEMPER, UTAU, UMU, PHI, H_LYR )  
deallocate( RHOQ, RHOU, EMUST, BEMST, RHO_ACCURATE )                
deallocate( RFLDIR, RFLDN, FLUP, DFDT, UAVG, TRNMED, UU )  
!deallocate( ALBMED)  

end subroutine deallocate_disort_allocatable_arrays

end module disort_vars

! c **********************************************************************






program DISORT_PROGRAM
use disort_vars
implicit none



real(kind=4),parameter :: PI = 2.*ASIN(1.0)
real(kind=4),parameter :: EARTH_RADIUS = 6371.0     
INTEGER :: IOD, IU, I, J, K, LC, LENTIT, LU, NPROB
INTEGER :: IBCND, NMOM, NLYR, NUMU, NSTR, NPHI, NTAU
CHARACTER  HEADER*127

LOGICAL  PRNT(5)
DATA PRNT  / .TRUE., 3*.FALSE., .TRUE. /

LOGICAL  USRANG, USRTAU, ONLYFL, &
         PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE
         
real(kind=4) :: ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, &
                PHI0, TEMIS, TTEMP, WVNMLO, WVNMHI, UMU0 

real, allocatable :: ALBMED

!real, allocatable :: ALBMED, DFDT, DTAUC, 
!      &          FLUP( MAXULV ), PHI( MAXPHI ), 
!      &          PMOM( 0:MAXMOM, MAXCLY ),
!      &          RFLDIR( MAXULV ), RFLDN( MAXULV ), SSALB( MAXCLY ),
!      &          TEMPER( 0:MAXCLY ), TRNMED( MAXUMU ), UAVG( MAXULV ),
!      &          UMU( MAXUMU ), UTAU( MAXULV ),
!      &          UU( MAXUMU, MAXULV, MAXPHI )


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
WRITE(*,*) NUMU
NPHI      = 1; 
IBCND     = 0     
         
!allocate( ALBMED( NUMU ))
call allocate_disort_allocatable_arrays( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU )


!allocate(DFDT( MAXULV ))
!allocate(DTAUC( MAXCLY ))
!      &          FLUP( MAXULV ), PHI( MAXPHI ), 
!      &          PMOM( 0:MAXMOM, MAXCLY ),
!      &          RFLDIR( MAXULV ), RFLDN( MAXULV ), SSALB( MAXCLY ),
!      &          TEMPER( 0:MAXCLY ), TRNMED( MAXUMU ), UAVG( MAXULV ),
!      &          UMU( MAXUMU ), UTAU( MAXULV ),
!      &          UU( MAXUMU, MAXULV, MAXPHI ))
     
! REAL(kind=4) :: ALBMED( MAXUMU ), DFDT( MAXULV ), DTAUC( MAXCLY ),
!      &          FLUP( MAXULV ), PHI( MAXPHI ), 
!      &          PMOM( 0:MAXMOM, MAXCLY ),
!      &          RFLDIR( MAXULV ), RFLDN( MAXULV ), SSALB( MAXCLY ),
!      &          TEMPER( 0:MAXCLY ), TRNMED( MAXUMU ), UAVG( MAXULV ),
!      &          UMU( MAXUMU ), UTAU( MAXULV ),
!      &          UU( MAXUMU, MAXULV, MAXPHI )
         
!RHOQ, RHOU, EMUST, BEMST RHO_ACCURATE, H_LYR

ALBMED = 0.0; 

!UMU0      = 0.5;  PHI0      = 0.0
UMU( 1 )  = -1.0; UMU( 2 )  = -0.2; UMU( 3 ) = 0.2; UMU( 4 ) = 1.0;
PHI( 1 )  = 60.0; 
UTAU( 1 ) = 0.0;  UTAU( 2 ) = 1.05; UTAU( 3 ) = 2.1; UTAU( 4 ) = 6.0; UTAU( 5 ) = 21.0; 

DO LC = 1, NLYR
  PMOM(0,LC) = 1.0
  DO  K = 1, NMOM
    PMOM(K,LC) = 0.0
  ENDDO
  DTAUC( LC ) = LC
  SSALB( LC ) = 0.6 + LC*0.05
ENDDO


FBEAM      = 0.0;   FISOT      = 1.0/PI
PLANK      = .FALSE.;     
LAMBER     = .TRUE.; ALBEDO   = 0.0; 

HEADER = ''

CALL DISORT( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU,           &
                 USRANG, USRTAU, IBCND, ONLYFL, PRNT,          &
                 PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,  &          
                 DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,   & 
                 UTAU, UMU0, PHI0, UMU, PHI, FBEAM,            &                        
                 FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,           &
                 EARTH_RADIUS, H_LYR,                          &
                 RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,       &
                 ACCUR,  HEADER,                               &
                 RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,          &
                 ALBMED, TRNMED )
                           
                    
                   
call deallocate_disort_allocatable_arrays()



WRITE(*,*) "Com-plete"

end program DISORT_PROGRAM


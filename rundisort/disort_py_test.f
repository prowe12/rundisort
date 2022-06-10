c   Modifications made by P. Rowe:
c       Parameter MXCLY (max layers) changed from 6 to 120
c       Commented out the following lines which display header
c          IF( .NOT.PASS1 .AND. LEN( HEADER ).NE.0 )
c        &    WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' )
c        &    ' DISORT: '//HEADER

c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: DISORT.f,v 2.1 2000/04/04 18:21:55 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,
     &                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
     &                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
     &                   FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,
     &                   PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY,
     &                   MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,
     &                   FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

c *******************************************************************
c       Plane-parallel discrete ordinates radiative transfer program
c             ( see DISORT.doc for complete documentation )
c *******************************************************************
c
c +------------------------------------------------------------------+
c  Calling Tree (omitting calls to ERRMSG):
c  (routines in parentheses are not in this file)
c
c  DISORT-+-(R1MACH)
c         +-SLFTST-+-(TSTBAD)
c         +-ZEROIT
c         +-CHEKIN-+-(WRTBAD)
c         |        +-(WRTDIM)
c         |        +-DREF
c         +-ZEROAL
c         +-SETDIS-+-QGAUSN-+-(D1MACH)
c         +-PRTINP
c         +-ALBTRN-+-LEPOLY
c         |        +-ZEROIT
c         |        +-SOLEIG-+-ASYMTX-+-(D1MACH)
c         |        +-TERPEV
c         |        +-SETMTX-+-ZEROIT
c         |        +-(SGBCO)
c         |        +-SOLVE1-+-ZEROIT
c         |        |        +-(SGBSL)
c         |        +-ALTRIN
c         |        +-SPALTR
c         |        +-PRALTR
c         +-PLKAVG-+-(R1MACH)
c         +-LEPOLY
c         +-SURFAC-+-QGAUSN-+-(D1MACH)
c         |        +-BDREF
c         |        +-ZEROIT
c         +-SOLEIG-+-ASYMTX-+-(D1MACH)
c         +-UPBEAM-+-(SGECO)
c         |        +-(SGESL)
c         +-UPISOT-+-(SGECO)
c         |        +-(SGESL)
c         +-TERPEV
c         +-TERPSO
c         +-SETMTX-+-ZEROIT
c         +-SOLVE0-+-ZEROIT
c         |        +-(SGBCO)
c         |        +-(SGBSL)
c         +-FLUXES--ZEROIT
c         +-ZEROIT
c         +-USRINT
c         +-CMPINT
c         +-PRAVIN
c         +-ZEROIT
c         +-RATIO--(R1MACH)
c         +-INTCOR-+-SINSCA
c         |        +-SECSCA-+-XIFUNC
c         +-PRTINT
c
c *** Intrinsic Functions used in DISORT package which take
c     non-negligible amount of time:
c
c    EXP :  Called by- ALBTRN, ALTRIN, CMPINT, FLUXES, SETDIS,
c                      SETMTX, SPALTR, USRINT, PLKAVG
c
c    SQRT : Called by- ASYMTX, SOLEIG
c
c +-------------------------------------------------------------------+
c
c  Index conventions (for all DO-loops and all variable descriptions):
c
c     IU     :  for user polar angles
c
c  IQ,JQ,KQ  :  for computational polar angles ('quadrature angles')
c
c   IQ/2     :  for half the computational polar angles (just the ones
c               in either 0-90 degrees, or 90-180 degrees)
c
c     J      :  for user azimuthal angles
c
c     K,L    :  for Legendre expansion coefficients or, alternatively,
c               subscripts of associated Legendre polynomials
c
c     LU     :  for user levels
c
c     LC     :  for computational layers (each having a different
c               single-scatter albedo and/or phase function)
c
c    LEV     :  for computational levels
c
c    MAZIM   :  for azimuthal components in Fourier cosine expansion
c               of intensity and phase function
c
c +------------------------------------------------------------------+
c
c               I N T E R N A L    V A R I A B L E S
c
c   AMB(IQ/2,IQ/2)    First matrix factor in reduced eigenvalue problem
c                     of Eqs. SS(12), STWJ(8E), STWL(23f)
c                     (used only in SOLEIG)
c
c   APB(IQ/2,IQ/2)    Second matrix factor in reduced eigenvalue problem
c                     of Eqs. SS(12), STWJ(8E), STWL(23f)
c                     (used only in SOLEIG)
c
c   ARRAY(IQ,IQ)      Scratch matrix for SOLEIG, UPBEAM and UPISOT
c                     (see each subroutine for definition)
c
c   B()               Right-hand side vector of Eq. SC(5) going into
c                     SOLVE0,1;  returns as solution vector
c                     vector  L, the constants of integration
c
c   BDR(IQ/2,0:IQ/2)  Bottom-boundary bidirectional reflectivity for a
c                     given azimuthal component.  First index always
c                     refers to a computational angle.  Second index:
c                     if zero, refers to incident beam angle UMU0;
c                     if non-zero, refers to a computational angle.
c
c   BEM(IQ/2)         Bottom-boundary directional emissivity at compu-
c                     tational angles.
c
c   BPLANK            Intensity emitted from bottom boundary
c
c   CBAND()           Matrix of left-hand side of the linear system
c                     Eq. SC(5), scaled by Eq. SC(12);  in banded
c                     form required by LINPACK solution routines
c
c   CC(IQ,IQ)         C-sub-IJ in Eq. SS(5)
c
c   CMU(IQ)           Computational polar angles (Gaussian)
c
c   CWT(IQ)           Quadrature weights corresponding to CMU
c
c   CORINT            When set TRUE, correct intensities for
c                     delta-scaling effects (see Nakajima and Tanaka,
c                     1988). When FALSE, intensities are not corrected.
c                     In general, CORINT should be set true when beam
c                     source is present (FBEAM is not zero) and DELTAM
c                     is TRUE in a problem including scattering.
c                     However, execution is faster when CORINT is FALSE,
c                     and intensities outside the aureole may still be
c                     accurate enough.  When CORINT is TRUE, it is
c                     important to have a sufficiently high order of
c                     Legendre approximation of the phase function. This
c                     is because the intensities are corrected by
c                     calculating the single-scattered radiation, for
c                     which an adequate representation of the phase
c                     function is crucial.  In case of a low order
c                     Legendre approximation of an otherwise highly
c                     anisotropic phase function, the intensities might
c                     actually be more accurate when CORINT is FALSE.
c                     When only fluxes are calculated (ONLYFL is TRUE),
c                     or there is no beam source (FBEAM=0.0), or there
c                     is no scattering (SSALB=0.0 for all layers) CORINT
c                     is set FALSE by the code.
c
c   DELM0             Kronecker delta, delta-sub-M0, where M = MAZIM
c                     is the number of the Fourier component in the
c                     azimuth cosine expansion
c
c   DELTAM            TRUE,  use delta-M method ( see Wiscombe, 1977 );
c                     FALSE, do not use delta-M method. In general, for
c                     a given number of streams, intensities and
c                     fluxes will be more accurate for phase functions
c                     with a large forward peak if DELTAM is set true.
c                     Intensities close to the forward scattering
c                     direction are often less accurate, however, when
c                     the delta-M method is applied. The intensity
c                     correction of Nakajima and Tanaka is used to
c                     improve the accuracy of the intensities.
c
c   DITHER            Small quantity subtracted from single-scattering
c                     albedos of unity, in order to avoid using special
c                     case formulas;  prevents an eigenvalue of exactly
c                     zero from occurring, which would cause an
c                     immediate overflow
c
c   DTAUCP(LC)        Computational-layer optical depths (delta-M-scaled
c                     if DELTAM = TRUE, otherwise equal to DTAUC)
c
c   EMU(IU)           Bottom-boundary directional emissivity at user
c                     angles.
c
c   EVAL(IQ)          Temporary storage for eigenvalues of Eq. SS(12)
c
c   EVECC(IQ,IQ)      Complete eigenvectors of SS(7) on return from
c                     SOLEIG; stored permanently in  GC
c
c   EXPBEA(LC)        Transmission of direct beam in delta-M optical
c                     depth coordinates
c
c   FLYR(LC)          Separated fraction in delta-M method
c
c   GL(K,LC)          Phase function Legendre polynomial expansion
c                     coefficients, calculated from PMOM by
c                     including single-scattering albedo, factor
c                     2K+1, and (if DELTAM=TRUE) the delta-M
c                     scaling
c
c   GC(IQ,IQ,LC)      Eigenvectors at polar quadrature angles,
c                     g  in Eq. SC(1)
c
c   GU(IU,IQ,LC)      Eigenvectors interpolated to user polar angles
c                     ( g  in Eqs. SC(3) and S1(8-9), i.e.
c                       G without the L factor )
c
c   IPVT(LC*IQ)       Integer vector of pivot indices for LINPACK
c                     routines
c
c   KK(IQ,LC)         Eigenvalues of coeff. matrix in Eq. SS(7)
c
c   KCONV             Counter in azimuth convergence test
c
c   LAYRU(LU)         Computational layer in which user output level
c                     UTAU(LU) is located
c
c   LL(IQ,LC)         Constants of integration L in Eq. SC(1),
c                     obtained by solving scaled version of Eq. SC(5)
c
c   LYRCUT            TRUE, radiation is assumed zero below layer
c                     NCUT because of almost complete absorption
c
c   NAZ               Number of azimuthal components considered
c
c   NCUT              Computational layer number in which absorption
c                     optical depth first exceeds ABSCUT
c
c   OPRIM(LC)         Single scattering albedo after delta-M scaling
c
c   PASS1             TRUE on first entry, FALSE thereafter
c
c   PKAG(0:LC)        Integrated Planck function for internal emission
c
c   PRNTU0(L)         logical flag to trigger printing of azimuthally-
c                     averaged intensities:
c                       L    quantities printed
c                      --    ------------------
c                       1    azimuthally-averaged intensities at user
c                               levels and computational polar angles
c                       2    azimuthally-averaged intensities at user
c                               levels and user polar angles
c
c   PSI0(IQ)          Sum just after square bracket in  Eq. SD(9)
c
c   PSI1(IQ)          Sum in  Eq. STWL(31d)
c
c   RMU(IU,0:IQ)      Bottom-boundary bidirectional reflectivity for a
c                     given azimuthal component.  First index always
c                     refers to a user angle.  Second index:
c                     if zero, refers to incident beam angle UMU0;
c                     if non-zero, refers to a computational angle.
c
c   SQT(k)            Square root of k (used only in LEPOLY for
c                     computing associated Legendre polynomials)
c
c   TAUC(0:LC)        Cumulative optical depth (un-delta-M-scaled)
c
c   TAUCPR(0:LC)      Cumulative optical depth (delta-M-scaled if
c                     DELTAM = TRUE, otherwise equal to TAUC)
c
c   TPLANK            Intensity emitted from top boundary
c
c   UUM(IU,LU)        Expansion coefficients when the intensity
c                     (u-super-M) is expanded in Fourier cosine series
c                     in azimuth angle
c
c   U0C(IQ,LU)        Azimuthally-averaged intensity at quadrature
c                     angle
c
c   U0U(IU,LU)        If ONLYFL = FALSE, azimuthally-averaged intensity
c                     at user angles and user levels
c
c                     If ONLYFL = TRUE and MAXUMU.GE.NSTR,
c                     azimuthally-averaged intensity at computational
c                     (Gaussian quadrature) angles and user levels;
c                     the corresponding quadrature angle cosines are
c                     returned in UMU.  If MAXUMU.LT.NSTR, U0U will be
c                     zeroed, and UMU, NUMU will not be set.
c
c   UTAUPR(LU)        Optical depths of user output levels in delta-M
c                     coordinates;  equal to  UTAU(LU) if no delta-M
c
c   WK()              scratch array
c
c   XR0(LC)           X-sub-zero in expansion of thermal source func-
c                     tion preceding Eq. SS(14)(has no mu-dependence);
c                     b-sub-zero in Eq. STWL(24d)
c
c   XR1(LC)           X-sub-one in expansion of thermal source func-
c                     tion; see  Eqs. SS(14-16); b-sub-one in STWL(24d)
c
c   YLM0(L)           Normalized associated Legendre polynomial
c                     of subscript L at the beam angle (not saved
c                     as function of superscipt M)
c
c   YLMC(L,IQ)        Normalized associated Legendre polynomial
c                     of subscript L at the computational angles
c                     (not saved as function of superscipt M)
c
c   YLMU(L,IU)        Normalized associated Legendre polynomial
c                     of subscript L at the user angles
c                     (not saved as function of superscipt M)
c
c   Z()               scratch array used in SOLVE0, ALBTRN to solve
c                     a linear system for the constants of integration
c
c   Z0(IQ)            Solution vectors Z-sub-zero of Eq. SS(16)
c
c   Z0U(IU,LC)        Z-sub-zero in Eq. SS(16) interpolated to user
c                     angles from an equation derived from SS(16)
c
c   Z1(IQ)            Solution vectors Z-sub-one  of Eq. SS(16)
c
c   Z1U(IU,LC)        Z-sub-one in Eq. SS(16) interpolated to user
c                     angles from an equation derived from SS(16)
c
c   ZBEAM(IU,LC)      Particular solution for beam source
c
c   ZJ(IQ)            Right-hand side vector  X-sub-zero in
c                     Eq. SS(19), also the solution vector
c                     Z-sub-zero after solving that system
c
c   ZZ(IQ,LC)         Permanent storage for the beam source vectors ZJ
c
c   ZPLK0(IQ,LC)      Permanent storage for the thermal source
c                     vectors  Z0  obtained by solving  Eq. SS(16)
c
c   ZPLK1(IQ,LC)      Permanent storage for the thermal source
c                     vectors  Z1  obtained by solving  Eq. SS(16)
c
c +-------------------------------------------------------------------+
c
c  LOCAL SYMBOLIC DIMENSIONS (have big effect on storage requirements):
c
c       MXCLY  = Max no. of computational layers
c       MXULV  = Max no. of output levels
c       MXCMU  = Max no. of computation polar angles
c       MXUMU  = Max no. of output polar angles
c       MXPHI  = Max no. of output azimuthal angles
c       MXSQT  = Max no. of square roots of integers (for LEPOLY)
c +-------------------------------------------------------------------+


c     FUNCTION INPUTS
c
c      SUBROUTINE DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,
c     &                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
c     &                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
c     &                   FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,
c     &                   PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY,
c     &                   MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,
c     &                   FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

c     Lines added for f2py:
cf2py intent(in)  :: NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,
cf2py intent(in)  :: WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
cf2py intent(in)  :: UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
cf2py intent(in)  :: FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,
cf2py intent(in)  :: PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY,
cf2py intent(in)  :: MAXULV, MAXUMU, MAXPHI, MAXMOM

cf2py intent(out) :: RFLDIR,RFLDN,FLUP,DFDT,UAVG,UU,ALBMED,TRNMED

c     .. Parameters ..

      INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI,
     &          MXSQT
      PARAMETER ( MXCLY = 120, MXULV = 5, MXCMU = 48, MXUMU = 10,
     &          MXPHI = 3, MI = MXCMU / 2, MI9M2 = 9*MI - 2,
     &          NNLYRI = MXCMU*MXCLY, MXSQT = 1000 )
c     ..
c     .. Scalar Arguments ..

      CHARACTER HEADER*127
      LOGICAL   LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCLY, MAXMOM, MAXPHI, MAXULV, MAXUMU, NLYR,
     &          NMOM, NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO
c     ..
c     .. Array Arguments ..

      LOGICAL   PRNT( 5 )
      REAL      ALBMED( MAXUMU ), DFDT( MAXULV ), DTAUC( MAXCLY ),
     &          FLUP( MAXULV ), PHI( MAXPHI ), PMOM( 0:MAXMOM, MAXCLY ),
     &          RFLDIR( MAXULV ), RFLDN( MAXULV ), SSALB( MAXCLY ),
     &          TEMPER( 0:MAXCLY ), TRNMED( MAXUMU ), UAVG( MAXULV ),
     &          UMU( MAXUMU ), UTAU( MAXULV ),
     &          UU( MAXUMU, MAXULV, MAXPHI )
c     ..
c     .. Local Scalars ..

      LOGICAL   COMPAR, CORINT, DELTAM, LYRCUT, PASS1
      INTEGER   IQ, IU, J, KCONV, L, LC, LEV, LU, MAZIM, NAZ, NCOL,
     &          NCOS, NCUT, NN, NS
      REAL      ANGCOS, AZERR, AZTERM, BPLANK, COSPHI, DELM0, DITHER,
     &          DUM, PI, RPD, SGN, TPLANK
c     ..
c     .. Local Arrays ..

      LOGICAL   PRNTU0( 2 )
      INTEGER   IPVT( NNLYRI ), LAYRU( MXULV )

      REAL      AMB( MI, MI ), APB( MI, MI ), ARRAY( MXCMU, MXCMU ),
     &          B( NNLYRI ), BDR( MI, 0:MI ), BEM( MI ),
     &          CBAND( MI9M2, NNLYRI ), CC( MXCMU, MXCMU ),
     &          CMU( MXCMU ), CWT( MXCMU ), DTAUCP( MXCLY ),
     &          EMU( MXUMU ), EVAL( MI ), EVECC( MXCMU, MXCMU ),
     &          EXPBEA( 0:MXCLY ), FLDIR( MXULV ), FLDN( MXULV ),
     &          FLYR( MXCLY ), GC( MXCMU, MXCMU, MXCLY ),
     &          GL( 0:MXCMU, MXCLY ), GU( MXUMU, MXCMU, MXCLY ),
     &          KK( MXCMU, MXCLY ), LL( MXCMU, MXCLY ), OPRIM( MXCLY ),
     &          PHASA( MXCLY ), PHAST( MXCLY ), PHASM( MXCLY ),
     &          PHIRAD( MXPHI ), PKAG( 0:MXCLY ), PSI0( MXCMU ),
     &          PSI1( MXCMU ), RMU( MXUMU, 0:MI ), SQT( MXSQT ),
     &          TAUC( 0:MXCLY ), TAUCPR( 0:MXCLY ), U0C( MXCMU, MXULV ),
     &          U0U( MXUMU, MXULV ), UTAUPR( MXULV ),
     &          UUM( MXUMU, MXULV ), WK( MXCMU ), XR0( MXCLY ),
     &          XR1( MXCLY ), YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ),
     &          YLMU( 0:MXCMU, MXUMU ), Z( NNLYRI ), Z0( MXCMU ),
     &          Z0U( MXUMU, MXCLY ), Z1( MXCMU ), Z1U( MXUMU, MXCLY ),
     &          ZBEAM( MXUMU, MXCLY )
      REAL      ZJ( MXCMU ), ZPLK0( MXCMU, MXCLY ),
     &          ZPLK1( MXCMU, MXCLY ), ZZ( MXCMU, MXCLY )

      DOUBLE PRECISION AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ),
     &                 WKD( MXCMU )
c     ..
c     .. External Functions ..

      REAL      PLKAVG, R1MACH, RATIO
      EXTERNAL  PLKAVG, R1MACH, RATIO
c     ..
c     .. External Subroutines ..

      EXTERNAL  ALBTRN, CHEKIN, CMPINT, FLUXES, INTCOR, LEPOLY, PRAVIN,
     &          PRTINP, PRTINT, SETDIS, SETMTX, SLFTST, SOLEIG, SOLVE0,
     &          SURFAC, TERPEV, TERPSO, UPBEAM, UPISOT, USRINT, ZEROAL,
     &          ZEROIT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, ASIN, COS, FLOAT, LEN, MAX, SQRT
c     ..
      SAVE      DITHER, PASS1, PI, RPD, SQT
      DATA      PASS1 / .TRUE. /, PRNTU0 / 2*.FALSE. /




      OPEN(8,FILE='dumfile',STATUS='unknown',FORM='FORMATTED')
      WRITE (8,*) 'hi Penny'
	  CLOSE(8)


      END
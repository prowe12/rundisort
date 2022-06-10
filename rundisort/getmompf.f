C $Id: getmompf.f,v 1.4 2007/10/11 21:03:27 dturner Release_2_1 $
C
      SUBROUTINE GETMOMPF(PHASE, ANGLE, NPHASE, MAXTABLEG, 
     :     SIZEPARAM, THETADEG, LEGEN, NLEG, VERBOSE, NEWPHASE, DONOR)

C       1/13/03
C     New version of cloudprp with increased flexibility.  Cloudprp.f is
C     now broken up into two programs.  The first one, create_tables.f 
C     is used to create tables of scattering properties of liquid, ice
C     or aerosol particles.  Liquid and aerosol properties are calculated
C     using Mie theory.  Ice properties are calculated using Ping Yang's
C     database of scattering properties for individual ice crystals.

C     The second program combines the optical properties from the
C     scattering tables with input information on mass content, particle
C     size and particle type at each grid point to create a .prp file
C     to be used by SHDOM.

C     Original author: Sally McFarlane (SAM)
C	2/26/03 Subroutine make_legendre_table modified by SAM

C     9/21/04 JMC: modified from create_legentable.f from SAM to interface with 
C     lbldis by D. Turner; creates legendre coefficients when provided the phase
C     function.  JMC -> Jennifer Comstock

C     2/11/2004 DDT: modified to allow this routine to normalize arbitrary
C	phase functions without fitting the truncation peak.  Therefore, the 
C 	most common use of this code is to call it once to normalize the phase
C	function (with the "DONOR" variable set to true), and then to call it
C	again to fit the forward peak and compute the legendre coefficients
C	("DONOR" set to false).

C   Inputs:
C     PHASE:     phase function
C     ANGLE: 	 the angles [deg] of the phase function
C     NPHASE:    the number of angles in the phase function
C     MAXTABLEG: maximum number of legendre coefficients
C     SIZEPARAM: the size parameter (2*Pi*r / lambda)
C     THETADEG:  width of broadened peak (<1.0)
C     VERBOSE:   If 1, then write output to screen; otherwise keep quiet.

C   Outputs:
C     LEGEN:     legendre coefficients
C     NLEG:      minimum number of legendre coefficients that should be used.
C     NEWPHASE:	 the new normalized phase function
C     DONOR:	 if set to 1, the phase function is normalized, but the fitting of
C			the forward peak and the computation of the Legendre
C			coefficients is not performed

      IMPLICIT NONE
      INTEGER  J, L

C Parameter variables 
      INTEGER   MAXLEG, MAXNQUAD
      PARAMETER (MAXLEG=5000, MAXNQUAD=3000 )
      REAL*8	  EPS1, EPS2
      PARAMETER (EPS1=0.01D0, EPS2=1.D-8)  
      REAL	  TRUNC
      PARAMETER   (TRUNC = 1.0E-1)

C Variables for input/output
      REAL*8      PHASE(*), NEWPHASE(*)
      REAL*8      LEGEN(0:MAXLEG), ANGLE(*), SIZEPARAM
      INTEGER     NLEG, MAXTABLEG, NPHASE, VERBOSE, DONOR
      REAL*8      THETADEG

C Local variables
      LOGICAL	PEAKFLAG
      REAL*8  	PI
      REAL*8	INTEG, SCALEFACTOR
      REAL*8    GAUSSMU(MAXNQUAD), GAUSSWTS(MAXNQUAD)
      INTEGER   NQUAD, NLEGEN
c---

      PI = ACOS(-1.D0)
      if(DONOR .eq. 1) then
        PEAKFLAG=.FALSE.      ! DONOR == 0 implies only normalize phase function
      else
        PEAKFLAG=.TRUE.       ! DONOR == 1 implies fit forward peak and generate L coefs
      endif

C     Initialize variables
      DO L=0,MAXLEG
         LEGEN(L)=0.0
      ENDDO
      DO L=1,NPHASE
         NEWPHASE(L) = 0.D0
      ENDDO

C     Estimate the number of Legendre terms needed, but limit NQUAD and 
C     NLEGEN to MAXNQUAD and MAXLEG, respectively.  This is just a
C     rough estimate.  You will need more Legendre/quadrature terms 
C     if you do not broaden the forward peak of the phase function.

      NLEGEN = MIN(MAXLEG,INT(30.0*SIZEPARAM))
      NQUAD = MIN(MAXNQUAD,NLEGEN)

      IF (NLEGEN .GT. MAXLEG .OR. NQUAD .GT. MAXNQUAD) THEN
         IF (VERBOSE .GE. 1) THEN
           WRITE(*,*) 'GETMOMPF: NLEGEN or NQUAD too large; ',
     .        'setting to MAXLEG and MAXNQUAD respectively.'
         ENDIF
         NLEGEN = MAXLEG
         NQUAD = MAXNQUAD
      ENDIF

C     Check whether the combined phase function is normalized.
C     Get the Gauss-Legendre quadrature abscissas and weights
      CALL GAUSQUAD (NQUAD, GAUSSMU, GAUSSWTS)

C     Calculate the integral of the phase function; normalized to 1.
C     \int_0^{2\pi} d\phi \int_0^{\pi} 
C     d\theta \sin(\theta) p(\theta,\phi)/(4*\pi) = 1.
C     \int_0^{\pi} d\theta \sin(\theta) p(\theta)/2 = 1.

      CALL CALC_INTEG(NPHASE, NQUAD, ANGLE, PHASE,
     .     GAUSSMU, GAUSSWTS, EPS1, INTEG)

      IF ( ABS(INTEG-1.d0) .GT. EPS1) THEN
         IF (VERBOSE .GE. 2) THEN
          WRITE(*,*) 'Combined tabulated phase function ',
     .        ' not normalized. Integrated value = ', INTEG
         ENDIF
         
         IF (PEAKFLAG) THEN
C     Fit the forward peak with a Gaussian function,
            IF (VERBOSE .GE. 2) THEN
               WRITE(*,*) 'Fitting Gaussian to forward peak'
            ENDIF

            CALL FIT_FORWARDPEAK(NPHASE,THETADEG,ANGLE,
     .           PHASE,EPS1,EPS2,NEWPHASE)
c ddt: If nphase is negative, then there was an error in the routine called above
            if(nphase .lt. 0) then
              return
            endif

            CALL CALC_INTEG(NPHASE, NQUAD, ANGLE, NEWPHASE,
     .           GAUSSMU, GAUSSWTS, EPS1, INTEG)
            IF (VERBOSE .GE. 2) THEN
               WRITE(*,*) 'Integration after forward peak fit: ', INTEG
            ENDIF

         ELSE
C     Renormalize by multiplicative scaling.
            IF (VERBOSE .GE. 2) THEN
               WRITE(*,*) 'Renormalizing by multiplicative scaling'
            ENDIF
            SCALEFACTOR = 1.D0/INTEG
            DO J = 1, NPHASE
               NEWPHASE(J) = PHASE(J)*SCALEFACTOR
            ENDDO
         ENDIF
         
      ELSE
         IF (VERBOSE .GE. 2) THEN
            WRITE(*,*) 'Phase function normalized; ', 
     .        'not fitting forward peak ', INTEG
         ENDIF
         DO J = 1, NPHASE
            NEWPHASE(J) = PHASE(J)
         ENDDO
      ENDIF

C     Calculate the Legendre series using Gaussian quadrature angles

      CALL LEGCALC(MAXLEG, NLEGEN,NQUAD,NPHASE,ANGLE,
     .     NEWPHASE, GAUSSMU, GAUSSWTS, NLEG, LEGEN)

C     If the Legendre series is being truncated, then fit the
C     forward peak, even if the normalization is okay.
      IF (LEGEN(MAXLEG) .GT. TRUNC) THEN
         IF (VERBOSE .GE. 1) THEN
           WRITE(*,*) 'WARNING: Legendre series is being truncated.',
     .        ' Thefore, the forward peak is being fit even though ',
     .        ' the phase function is normalized.  Another option ',
     .        ' is to increase the value of MAXLEG.'
         ENDIF
         IF (.NOT. PEAKFLAG) THEN
            THETADEG = 0.25
         ENDIF
         CALL FIT_FORWARDPEAK(NPHASE,THETADEG,ANGLE,
     .        PHASE,EPS1,EPS2,NEWPHASE)
c ddt: If nphase is negative, then there was an error in the routine called above
            if(nphase .lt. 0) then
              return
            endif
         CALL LEGCALC(MAXLEG, NLEGEN,NQUAD,NPHASE,ANGLE,
     .        NEWPHASE, GAUSSMU, GAUSSWTS, NLEG, LEGEN)

         CALL CALC_INTEG(NPHASE, NQUAD, ANGLE, NEWPHASE,
     .        GAUSSMU, GAUSSWTS, EPS1, INTEG)
         IF (VERBOSE .GE. 2) THEN
            WRITE(*,*) 'Integration of phase function after ',
     .        'forward peak fit: ', INTEG
         ENDIF

      ENDIF

C     Make sure the Legendre polynomials are correctly normalized,
C     If not, use simple multiplicative factor to rescale them
      
      IF ( ABS(LEGEN(0) - 1.d0) .GT. EPS1) THEN
         IF (VERBOSE .GE. 1) THEN
            WRITE(*,*)  
     .        ' Legendre polynomials still not correctly normalized; ',
     .        'using multiplicative scaling to renormalize - but ',
     .        'you might want to increase MAXLEG'
         ENDIF
      ENDIF
      SCALEFACTOR = 1.0/LEGEN(0)
      DO L = 0, NLEGEN
         LEGEN(L) = LEGEN(L) * SCALEFACTOR
         LEGEN(L) = LEGEN(L) / (2.0 * L + 1.0) 
      ENDDO 

      RETURN
      END


C ------------------------------------------------------

      SUBROUTINE FIT_FORWARDPEAK(NPHASE,THETADEG,ANGLE,
     .                PHASE,EPS1,EPS2,NEWPHASE)

C Subroutine to fit the forward peak of the phase function
C  with a Gaussian such
C       that the entire phase function will be normalized and the
C       Gaussian peak will match the remainder of the phase function
C       at theta_0.
          
C These conditions can be rewritten as:
C       1) \int_{\cos(\theta_0)}^1 \frac{A}{\sigma \sqrt(2.\pi)}
C               \exp(-0.5\theta^2/\sigma^2)
C               d\cos(\theta) = 2. - truncinteg
C       2) \frac{A}{\sigma \sqrt(2.*\pi)}
C               \exp(\exp(-0.5\theta_0^2/\sigma^2) = P(\theta_0)
C where A is the amplitude and \sigma the width of the Gaussian
C distribution, and P(\theta_0) is the value of the original
C phase function at \theta_0.
C
C The first constraint can be rewritten as:
C       \int_0^{\theta_0} \frac{A}{\sigma \sqrt(2 \pi)}
C       \exp(-0.5\theta^2/\sigma^2) \sin(\theta) d(\theta) = 2. - TRUNCINTEG
C  Assuming that \theta << 1, you can use the small angle approximation,
C  \sin(\theta) ~ \theta.  Therefore the constraint becomes
C       \int_0^{\theta_0} \frac{A}{\sigma \sqrt(2 \pi)}
C       \exp(-0.5\theta^2/\sigma^2) \theta d(\theta) = 2. - TRUNCINTEG
C  This equation can be solved analytically to give:
C       \frac{A \sigma}{\sqrt(2 \pi)}(1-\exp(-0.5\theta_0^2/\sigma^2) = 2 - TRUNCINTEG
C
C  To solve the equations, find A in terms of P(\theta_0):
C       A = P(\theta_0) \sigma \sqrt(2\pi) \exp(-0.5\theta_0^2/\sigma^2)
C  Now substitute this value for A into the analytic form of the
C  integral constraint, which gives:
C       P(\theta_0)\sigma^2(\exp(0.5\theta_0^2/\sigma^2) - 2 + truncinteg = 0
C  The problem has then become simply to find the roots of this equation.


      IMPLICIT 	NONE
      INTEGER	I, J, K, L

C Input variables
      INTEGER	NPHASE
      REAL*8	THETADEG, ANGLE(*)
      REAL*8    PHASE(*)  

C Output variables
      REAL*8 	NEWPHASE(*)

C Local variables
      REAL	THETA, THETACOS
      INTEGER	THETAINDEX
      REAL*8	PI,P1, MU(NPHASE), TRUNCINTEG, INTEG
      REAL	RADANGLE(NPHASE)
      DOUBLE PRECISION GAUSSFUNC        
      INTEGER   NINTERV, NMAX, NROOTS   
      PARAMETER (NINTERV=20,NMAX=5)
      REAL*8 RTBIS
      REAL*8 EPS1, EPS2
      REAL*8  XB1(NMAX), XB2(NMAX)
      REAL*8  AMP(NMAX), SIGMA(NMAX) 

C----------------  Begin processing ---------------------------

      PI = ACOS(-1.D0)

C Set the width of the forward peak in radians
      THETA = THETADEG*PI/180.
      THETACOS = COS(THETA)
C Convert the tabulated angles into radians
      DO J = 1, NPHASE
        RADANGLE(J) = ANGLE(J)*PI/180.
      ENDDO
          
C Get the value of the phase function at theta and get thetaindex
      CALL PHASE_ANGLE(THETACOS,PHASE,ANGLE,NPHASE, P1,
     .          THETAINDEX)
      THETAINDEX = THETAINDEX - 1

C Do the integration of the truncated part of the phase function
C (use trapezoid rule for speed)
      DO J = 1, NPHASE
        MU(J) = COS(ANGLE(J)*PI/180.)
      ENDDO
      TRUNCINTEG = 0.0
      DO J = THETAINDEX, NPHASE-1
        TRUNCINTEG = TRUNCINTEG +
     .  	ABS(MU(J)-MU(J+1))*
     .          (PHASE(J)+PHASE(J+1))/2.
      ENDDO

C Now begin to find the roots of the equation:
C P(\theta_0)\sigma^2(\exp(0.5\theta_0^2/\sigma^2) - 2 + TRUNCINTEG = 0

C The first step is to bracket the roots with brak -- this routine
C splits the given interval up into ninterv subintervals, and
C returns a maximum of nmax roots.  Note that the intervals  
C are in log-space, so the first endpoint cannot be zero.    


      CALL BRAK(.001d0, 100.d0, NINTERV, XB1, XB2,
     .          NMAX,NROOTS,THETA,P1,TRUNCINTEG)

c      NROOTS = 0
c      DO J = 1, NMAX
c        IF (XB1(J) .GT. 0 .AND. XB2(J) .GT. 0) THEN
c          NROOTS = J
c        ENDIF
c      ENDDO

c ddt: modified to return a negative value if roots could not be found;
c ddt:   this allows the MIXCRA code to exit in a 'nicer' way 
      IF (NROOTS .LT. 1) THEN
        WRITE(*,*) 'No roots found; check theta'
        nphase = -1
        return
      ELSE IF (NROOTS .GT. 1) THEN
        WRITE(*,*) 'WARNING: more than one root found, ',
     .          'check bracketing interval'
      ENDIF

C Now take the bracketed intervals and try to actually find
C the roots -- roots must be within EPS of 0.0.  Note
C that there may be more than one root - in that case,
C notify user, but just use first root.

      DO J = 1, NROOTS
        SIGMA(J) =  RTBIS(XB1(J), XB2(J), EPS2, THETA,
     .          P1, TRUNCINTEG)
c Given the width, find the amplitude of the Gaussian
        AMP(J) = P1*SIGMA(J)*SQRT(2.*PI) *
     .          EXP((THETA**2)/(2.*SIGMA(J)**2))
      ENDDO

C Do a couple of tests to make sure the roots actually
C       fit the constraints correctly.  Check that
C       constraints are matched within EPS1.
          
      IF ( ABS(P1*SIGMA(1)**2*(EXP(THETA**2/(2.*SIGMA(1)**2))-1)
     .          + TRUNCINTEG - 2.0) .GT. EPS1) THEN
        STOP 'Normalization constraint does not agree within EPS1'
      ENDIF
      IF ( ABS( (AMP(1) / (SIGMA(1)*SQRT(2.*PI)) ) *
     .          EXP(-(THETA**2)/(2.*SIGMA(1)**2))
     .          - p1) .GT. EPS1) THEN
        STOP 'Phase function does not agree w/in EPS1 at theta_0'
      ENDIF

C Now add this gaussian onto the phase function
     
      DO J = 1, thetaindex-1
        NEWPHASE(J) = GAUSSFUNC(RADANGLE(J),
     .          SIGMA(1), AMP(1))
      ENDDO
      DO J = THETAINDEX, NPHASE
        NEWPHASE(J) = PHASE(J)
      ENDDO
          
C Check whether entire phase function is normalized within eps. 
      INTEG = 0.0
      DO j = 1, nphase-1
        INTEG = INTEG +
     .      abs(MU(J)-MU(J+1))*
     .      (NEWPHASE(J)+NEWPHASE(J+1))/2.
      ENDDO
      IF (ABS(integ-2.d0) .GT. EPS1) THEN
        STOP 'New phase function is not normalized within EPS1'
      ENDIF

      RETURN
      END


C ------------------------------------------------------------

      SUBROUTINE LEGCALC(MAXLEG,NLEGEN,NQUAD,NPHASE,TABANGLE,
     .		PHASE, QUADMU, QUADWTS, NLEG, LEGEN)


C       Computes the Legendre series expansion of a phase
C     function PHASE tabulated at a set of angles TABANGLE
C     given a set of quadrature angles with cosines QUADMU 
C     and weights QUADWTS.

      IMPLICIT NONE

C Input variables
      INTEGER	MAXLEG, NQUAD, NLEGEN, NPHASE
      REAL*8	PHASE(NPHASE),QUADMU(*),QUADWTS(*)
      REAL*8	TABANGLE(NPHASE)

C Output variables
      INTEGER	NLEG
      REAL*8  	LEGEN(0:MAXLEG)

C Local variables
      INTEGER	J, L
      REAL	QUADFACTOR(NQUAD), QUADANGLE
      INTEGER	QUADINDEX(NQUAD)
      REAL*8	P1,PI, INTEG, PL1, PL, PL2, COEF1(0:NLEGEN)


      PI = ACOS(-1.D0)
      DO L = 0, NLEGEN
        COEF1(L) = 0.0
      ENDDO

C Find the index and interpolation factor
      DO J = 1, NQUAD
        QUADANGLE = ACOS( QUADMU(J) ) * 180.0/PI
        CALL BACKINTERPOLATE(NPHASE, TABANGLE, QUADANGLE,
     .                QUADFACTOR(J),QUADINDEX(J))
      ENDDO

      INTEG = 0.0
      DO J = 1, NQUAD
        P1 = PHASE(QUADINDEX(J)) +
     .          ( PHASE(QUADINDEX(J)+1) -
     .            PHASE(QUADINDEX(J)) )*QUADFACTOR(J)

C     Use upward recurrence to find Legendre polynomials
          
        PL1 = 1.0
        PL = 1.0
        DO L = 0, NLEGEN
          IF (L .GT. 0) PL = (2*L-1)*QUADMU(J)*PL1/L - (L-1)*PL2/L
          COEF1(L) = COEF1(L) + PL*QUADWTS(J)*P1
          PL2 = PL1
	  PL1 = PL
        ENDDO
        INTEG = INTEG + P1*QUADWTS(J)
      ENDDO

      DO L = 0, NLEGEN
        LEGEN(L) = (2*L+1)/2.0 * COEF1(L)
        IF (LEGEN(L) .GT. 1.0E-5) THEN
          NLEG = L
        ELSE
          GOTO 50
        ENDIF
      ENDDO
50    CONTINUE

      RETURN
      END



C ----------------------------------------------------
      SUBROUTINE CALC_INTEG(NPHASE,NQUAD,TABANGLE,PHASE,
     .		QUADMU,QUADWTS,EPS1,INTEG)

C	Subroutine calculates the integral of a phase function 
c	PHASE, tabulated at angles TABANG given a set of quadrature 
C       angle cosines, QUADMU, and weights, QUADWTS

      IMPLICIT	NONE

C Input variables
      INTEGER	NPHASE, NQUAD
      REAL*8	PHASE(*),QUADMU(*),QUADWTS(*), EPS1
      REAL*8	TABANGLE(*)
C Output variables
      REAL*8	INTEG
C Local variables
      INTEGER	J
      REAL	QUADANGLE, QUADFACTOR(NQUAD)
      REAL*8	P1, PI
      INTEGER	QUADINDEX(NQUAD)

      PI = ACOS(-1.D0)

C Get the index and factor to interpolate the tabulated
C phase function to the quadrature angles.
      DO J = 1, NQUAD
        QUADANGLE = ACOS( QUADMU(J) ) * 180.0/PI
        CALL BACKINTERPOLATE(NPHASE,TABANGLE, QUADANGLE,
     .                QUADFACTOR(J),QUADINDEX(J))
      ENDDO
            
      INTEG = 0.D0
      DO J = 1, NQUAD
        P1 = PHASE(QUADINDEX(J)) +
     .      (PHASE(QUADINDEX(J)+1) -
     .       PHASE(QUADINDEX(J)) )*QUADFACTOR(J)  
        INTEG = INTEG + P1*QUADWTS(J)
      ENDDO
      INTEG = INTEG/2.D0

      RETURN
      END


c------------------------------------------------------------------
      SUBROUTINE GAUSQUAD (N, XA, WT)
C        Generates the abscissas (X) and weights (W) for an N point
C      Gauss-Legendre quadrature.  
      IMPLICIT NONE
      INTEGER  N
      REAL*8   XA(*), WT(*)
      INTEGER  K, I, J, L
      REAL*8   X, XP, PL, PL1, PL2, DPL, TINY
      PARAMETER (TINY=3.0D-13)

      K = (N+1)/2
      DO 130 J = 1, K
        X = COS(3.141592654*(J-.25)/(N+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO 120 L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
120       CONTINUE
          DPL = N*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        XA(J)     = -X
        XA(N-J+1) = X
        WT(J  )   = 2.0D0/((1.0D0-X*X)*DPL*DPL)
        WT(N-J+1) = WT(J)
130   CONTINUE

      RETURN
      END
c***************************************************************************
      

C     obtain the phase function value P11 that corresponds to MU

      SUBROUTINE PHASE_ANGLE(MU,PHASE,ANGLE,NPHASE,P11,PINDEX)
      IMPLICIT NONE
      INTEGER   NPHASE, I, PINDEX
      REAL*8    ANGLE(NPHASE)
      REAL      MU, MU_ANG
      REAL*8    PI, P11, PHASE(NPHASE)
      LOGICAL   FOUND

      
      PI = ACOS(-1.0D0)
      MU_ANG = ACOS( MU ) * 180.0/PI
      FOUND = .FALSE.
      I=1
         
      DO WHILE (.NOT. FOUND)
         IF (ANGLE(I) .GE. MU_ANG) THEN
            P11 = PHASE(I)
            PINDEX = I
            FOUND = .TRUE.
         ELSE IF (I .GT. NPHASE) THEN
            PRINT *, 'PHASE_ANGLE: ANGLE NOT FOUND.... ABORT!'
            STOP
         ENDIF
         I=I+1
      ENDDO
         
      RETURN
      END
      
         
C ----------------------------------------------------------------

      SUBROUTINE brak(x1,x2,n,xb1,xb2,nb,nbb, 
     .          theta,p1,trunc)
      IMPLICIT NONE
      INTEGER n,nb
      DOUBLE PRECISION x1,x2,xb1(nb),xb2(nb)
      INTEGER i,nbb
      DOUBLE PRECISION dx,fc,fp,x
      DOUBLE PRECISION p1,trunc
      REAL  THETA

      nbb=0
      x=x1
      dx=(LOG10(x2)-LOG10(x1))/n
      CALL funcs(THETA,X,p1,trunc,fp)
      do 11 i=1,n
         x=10.**(LOG10(x)+dx)
         CALL funcs(THETA,X,p1,trunc,fc)
         if(fc*fp.lt.0.) then
            nbb=nbb+1
            xb1(nbb)=10.**(LOG10(x)-dx)
            xb2(nbb)=x
            if(nbb.eq.nb) goto 1
         endif
         fp=fc
 11   continue
 1    continue
c      nb=nbb

      return
      END   
      
c----------------------------------------------------------------------
      FUNCTION rtbis(x1,x2,xacc,theta,P1,trunc)
      implicit none

      INTEGER JMAX
      DOUBLE PRECISION rtbis,x1,x2,xacc
      PARAMETER (JMAX=400)
      INTEGER j
      DOUBLE PRECISION dx,f,fmid,xmid
      REAL theta
      DOUBLE PRECISION P1, trunc

      CALL funcs(theta,x2,P1,trunc,fmid)
      CALL funcs(theta,x1,p1,trunc,f)
      if(f*fmid.ge.0.) pause 'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1 
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        CALL funcs(theta,xmid,P1,trunc,fmid)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      pause 'too many bisections in rtbis'
      END

            
         

C ***********************************************

      SUBROUTINE BACKINTERPOLATE(nlev,array, z, x, j)

      implicit none
      integer nlev
      real    z, x
      real*8  array(nlev)
      integer j, jl, ju, jm, k


        jl = 0
        ju = nlev + 1
10      if (ju-jl .gt. 1) then
          jm = (ju+jl)/2
           if  ( (array(nlev) .ge. array(1)) .eqv.
     .           (z .ge. array(jm))) then
            jl=jm
           else
            ju=jm
          endif
          goto 10  
        endif
        j = jl
        if (j .le. 0 ) then
          j = 1
          x = 0.0
        else if (j .ge. nlev) then
          j = nlev-1
          x = 1.0
        else   
          x =  (z - array(j))/(array(j+1)-array(j))
        endif
         
      RETURN
      END

c***************************************************************************

      SUBROUTINE FUNCS(THETA, SIG, P1, TRUNCINTEG,   
     .          OUTPUTFUNC)


      IMPLICIT NONE
      DOUBLE PRECISION OUTPUTFUNC, SIG, PI,
     .  P1, TRUNCINTEG, TEMPTHETA
      REAL THETA
      

C Angle must be in radians

      PI = ACOS(-1.D0) 
      TEMPTHETA = THETA*1.d0

      OUTPUTFUNC = (P1*sig**2) * 
     .          (EXP( (TEMPTHETA**2)/(2.*SIG**2) )-1.) +
     .          TRUNCINTEG - 2.D0

      RETURN
      END

C -------------------------------------------------------

      FUNCTION GAUSSFUNC(THETA, SIG, AMP)
     
      IMPLICIT NONE
      DOUBLE PRECISION GAUSSFUNC, SIG, PI, AMP
      REAL THETA

C Angle must be in radians
      
      PI = ACOS(-1.D0)

      GAUSSFUNC = ( AMP / (SIG*SQRT(2.*PI)) ) *
     .          EXP( -( (THETA)**2 ) /
     .          (2.*SIG**2))

      RETURN
      END


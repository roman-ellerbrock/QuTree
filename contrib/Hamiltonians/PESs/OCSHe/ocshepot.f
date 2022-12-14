      subroutine init_pot()
      implicit none

c     initialisation of potential
      open(5,file='ocsheparas.dat',status='old',action='read')
	stop "did not compile"
c      call EXTINT
      close(5)
      
      end

c     calculates the energy at one point He-OCS potential energy surface 
c     of joanna m. m. howson and jeremy d. hutson from 2001
      subroutine get_pot_en(x,y,z,en)
      implicit none
      
c     coordinates
      real*8 z   !coordinate along the molecular axis
      real*8 x,y !coordinates perpendicular to 
                 !the molecular axis
      
      double precision en !return value (energy)

      double precision r     !distance to center of mass
      double precision costh !angle to molecular achses

c     functions
      double precision extpot
      
      r=sqrt(x**2+y**2+z**2)
      costh=z/r

      en = extpot(r,z/r)
      if((en.ge.10000.).or.(en.ne.en)) en = 10000.

      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% FUNCTION EXTPOT(R,COSTH)
C% Author: J. M. M. Howson and J. D. Hutson 2001
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION EXTPOT(R,COSTH)
c
c     Program to project out in legendre polys
c     in the angular co-ordinate.
c
      implicit double precision (a-h,o-z)
      integer mxparm, mxdata, mxsite
      save
      
      PARAMETER (MXSITE=10)
      parameter (mxparm=500,mxdata=500)
      
      double precision parm0(mxparm), rdata(mxdata),RLAM(50),BLAM(50)
      double precision p0(mxparm),vx(10),ax(10),sumleg
      double precision alpha(19,13),x(19),rad_kernel,rnw,cthnw
      double precision x_k,r,costh,pi,plm,pot,vpp,ra,ro,rs,rc,rb
      double precision eps_scal, r_scal, extpot, extint,ang,aur,aue
      
      integer n_r,n_lam, np,k
      integer i,j,n,npts
c      common npts
C     
      DATA RLAM/50*0.D0/, BLAM/50*0.D0/
      DATA AUE/219474.6354D0/
      DATA AUR/0.52917706D0/
C     
      ro = 1.68029d0
      rs = 1.03731d0

      ra = (r*r+rs*rs+2.0d0*r*rs*costh)**0.5d0
      rb = (r*r+ro*ro-2.0d0*r*ro*costh)**0.5d0

      rnw = (ra + rb) /(ro+rs)
      cthnw = (ra - rb) /(ro+rs)

      eps_scal = sumleg(vx,np,cthnw)
      
      r_scal = rnw*sumleg(ax,np,cthnw)
c     write(*,*)'r=',r_scal,'costh=',costh,'eps_scal=',eps_scal
      
      pi=dacos(-1.d0)
      
      pot=0.0d0

      do j = 1,n_lam
         
         vpp = 0.0d0
         do i=1,npts 
            x_k=x(i)
            vpp=vpp+alpha(i,j)*rad_kernel(x_k,r_scal)
         enddo

         pot = pot+ vpp*plm(j-1,0,cthnw)
c         pot = pot+ vpp*legp(j-1,cthnw)
      enddo
*
      extpot=eps_scal*pot
*
      return

	stop "did not compile"
c      ENTRY EXTINT
c     read number of scaling paras
      read(5,*) np
c     read scaling paras
      do i=1,np
	 read(5,*) vx(i), ax(i)
      enddo
c     read number of alphas(=no. of rvalues and theta values)
      read(5,*) npts, n_lam
c     NB:x(i)=radial pt at which have abinitio pt
c     NB:y(i)=(1-cos(theta))/2.0
c     where theta is angular pt at which have abinitio pt
c     read x,y and alphas for  RK
      
      do j = 1,n_lam
         do i=1, npts
            read(5,*) x(i),alpha(i,j)
         enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% FUNCTION RAD_KERNEL(X,R)
C% Author: J. M. M. Howson and J. D. Hutson 2001
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function rad_kernel(x,r)
      implicit none
      real*8 rad_kernel,x,r,x_l,x_s,R0
      
      x_l=r
      x_s=x
      if(r.lt.x) then
         x_s=r
         x_l=x
      endif
c     kernel for protonated systems (long range prop to R-4)
c     rad_kernel=2.d0/(15.d0*x_l**5)*(1.d0-5.d0*x_s/(7.d0*x_l))
c     kernel for protonated systems (long range prop to R-6)m=5,n=2
c     rad_kernel=2.d0/(21.d0*x_l**6.0d0)*(1.d0-3.d0*x_s/(4.d0*x_l))
c     m=6, n=2
      rad_kernel=1.d0/(14.d0*x_l**7.0d0)*(1.d0-7.d0*x_s/(9.d0*x_l))
c     Increase order of smoothness.n=3
c     rad_kernel=1.0d0/(28.0d0*x_l**7.0d0)*(1.d0-7.d0*x_s/(5.d0*x_l)
c     . + 28.d0*x_s*x_s/(55.d0*x_l*x_l))
c     n=3 m=5
c     rad_kernel=3.d0/(56.0d0*x_l**6)*(1.d0-4.d0*x_s/(3.d0*x_l)
c     . + 7.d0*x_s*x_s/(15.d0*x_l*x_l))
      
c     rad_kernel=1.d0/(14.0d0*(R0+x_l)**7)*
c     .    (1.d0-7.d0*(x_s+R0)/(9.d0*(x_l+R0)))
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% FUNCTION PLM(L,MM,X)
C% Author: J. M. M. Howson and J. D. Hutson 2001
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION PLM(L,MM,X)
C
C     COMPUTES NORMALIZED ASSOC. LEGENDRE POLYNOMIALS BY RECURSION.
C     THE VALUES RETURNED ARE NORMALIZED FOR INTEGRATION OVER X
C     (I.E. INTEGRATION OVER COS THETA BUT NOT PHI).
C     NOTE THAT THE NORMALIZATION GIVES 
C           PLM(L,0,1)=SQRT(L+0.5)
C           PLM(L,0,X)=SQRT(L+0.5) P(L,X)
C     FOR M.NE.0, THE VALUE RETURNED DIFFERS FROM THE USUAL
C           DEFINITION OF THE ASSOCIATED LEGENDRE POLYNOMIAL
C           (E.G. EDMONDS PAGES 23-24) 
C           BY A FACTOR OF (-1)**M*SQRT(L+0.5)*SQRT((L-M)!/(L+M)!)
C     THUS THE SPHERICAL HARMONICS ARE 
C          CLM = PLM * EXP(I*M*PHI) / SQRT(L+0.5)
C          YLM = PLM * EXP(I*M*PHI) / SQRT(2*PI)
C
C     PROGRAM OF R. NERF MODIFED BY S. GREEN.
C       MOD FEB. 82 BY S.G. ACCORDING TO R.T PACK'S SUGGESTION
C       FOR IMPROVED ACCURACY LARGE ABS(X)
C     MODIFIED AUG 88 BY S.G. TO KEEP RAT FROM BLOWING UP AT LARGE L
C       BY INCLUDING FACTOR (-Z)**L IN PRECEDING LOOP
C     ERROR IN STMT. NO. 211 CORRECTED MAY 93 BY SG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      M=IABS(MM)
      IF ( M.GT.L .OR. L.LT.0)  GO TO 9999
      XL=L
      XM=M
      P1=1.D0
      P2=X
      IF(L.EQ.0) GO TO 20
C *** USE ALTERNATE RECURSION FOR LARGE M OR LARGE ABS(X)
      IF (M.EQ.0)  GO TO 10
      XTEST=0.49999D0*(1.D0+1.D0/XM)
      MTEST=L/3
      IF (M.GT.MTEST .OR. ABS(X).GT.XTEST)  GO TO 210
C ***
   10 DO 100 I=1,L
      XI=I
      P3=((2.D0*XI+1.D0)*X*P2-XI*P1)/(XI+1.D0)
      P1=P2
  100 P2=P3
      IF (M.EQ.0) GO TO 20
C AT END OF LOOP P1=P(L,0,X)
      IF (ABS(X).GT.1.D0)  GO TO 9999
      Z=SQRT(1.D0-X*X)
      IF (Z.LE.1.D-10) GO TO 999
  201 P2=(XL+1.D0)*(P2-X*P1)/Z
      DO 200 I=1,M
      XI=I
      P3=-2.D0*X/Z*P2*XI-(XL+XI)*(XL-XI+1.D0)*P1
      P1=P2
  200 P2=P3
      GO TO 20
C ***
C     BELOW RECURS DOWN IN M, AS SUGGESTED BY R. T PACK
  210 IF (ABS(X).GT.1.D0)  GO TO 9999
      Z=SQRT(1.D0-X*X)
      IF (Z.LE.1.D-10)  GO TO 999
C     CALCULATE RATIO OF FACTORIALS FOR PLL
      RAT=1.D0
      XI=0.D0
      DO 211 I=1,L
      XI=XI+1.D0
  211 RAT=-0.5D0*Z*(XL+XI)*RAT
C     CALCULATE PLL   ----  N.B. ABOVE INCL (-Z)**L COMPUTATION
C     P1=RAT*(-Z)**L
      P1=RAT
      IF (M.EQ.L)  GO TO 20
      P2=P1
      P1=-X*P2/Z
      IF (M.EQ.L-1)  GO TO 20
      LM1=L-M-1
C     RECUR DOWNWARD IN M
      DO 213 I=1,LM1
      MU=L-I-1
      XMU=MU
      P3=P2
      P2=P1
      P1=2.D0*(XMU+1.D0)*X*P2/Z+P3
  213 P1=-P1/((XL-XMU)*(XL+XMU+1.D0))
      GO TO 20
C ***
C     NORMALIZATION . . .
   20 XNORM=(2.D0*XL+1.D0)/2.D0
      IF (M.LE.0)  GO TO 1000
      XLM=XL+1.D0
      XLP=XL
      DO 1100 I=1,M
      XLM=XLM-1.D0
      XLP=XLP+1.D0
 1100 XNORM=XNORM/(XLM*XLP)
 1000 PLM=P1*SQRT(XNORM)
      RETURN
 9999 WRITE(6,699)  L,MM,X
  699 FORMAT('0 * * * ERROR.  ARGUMENT OUT OF RANGE FOR PLM(',2I6,D16.8,
     1  ' ).')
C     IF Z=0, THEN X=1 AND PLM(1.)=0 FOR M.GT.0.
C       IN THAT CASE, WE HAVE BRANCHED TO 999
  999 PLM=0.D0
      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% FUNCTION SUMLEG(COEFF,NP,X)
C% Author: J. M. M. Howson and J. D. Hutson 2001
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION SUMLEG(COEFF,NP,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C   
C     SUMLEG EVALUATES A LEGENDRE SERIES AT A GIVEN ANGLE THETA, USING
C     THE RECURSION RELATIONSHIP FOR LEGENDRE POLYNOMIALS.
C
C     ON INPUT, COEFF IS THE LEGENDRE SERIES, STARTING AT P0
C               NP    IS THE ORDER OF THE LEGENDRE SERIES
C               X     IS COS(THETA)
C
      DIMENSION COEFF(NP)
      IF(NP.GE.1) GOTO 1
      WRITE(6,601)NP
  601 FORMAT(/' **** ERROR IN SUMLEG: NP =',I5)
      STOP
    1 SUMLEG=COEFF(1)
      IF(NP.EQ.1) RETURN
      SUMLEG=SUMLEG+X*COEFF(2)
      IF(NP.EQ.2) RETURN
      P0=1.D0
      P1=X
      DO 10 K=3,NP
      TEMP=(DBLE(K+K-3)*X*P1 - DBLE(K-2)*P0) / DBLE(K-1)
      P0=P1
      P1=TEMP
   10 SUMLEG=SUMLEG+P1*COEFF(K)
      RETURN
      END


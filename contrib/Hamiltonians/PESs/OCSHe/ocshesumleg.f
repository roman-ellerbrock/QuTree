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


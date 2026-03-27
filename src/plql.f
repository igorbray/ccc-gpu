C  FILE 'plql.f'
C
      SUBROUTINE FAKRED
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'par.f'
      parameter (lll=lmax*5+20)
      COMMON/REDFAK/YNN,FAK(lll)
C
C ****  COMPUTE AND STORE,  FAK(N) = (N-1)!/(YNN**(N-1))
C
      YNN = 1e4 !500.0D0
      FAK(1) = 1.0D0
      DO 10 I = 1,lll-1
      FAK(I+1) = FAK(I)/YNN*DBLE(I)
   10 CONTINUE
      RETURN
      END

      function fakdiv(n,m)      !(n-1)!/(m-1)! for m >= n
      real*8 fakdiv
      fakdiv = 1d0
      do i = n, m-1
         fakdiv = fakdiv / dfloat(i)
      enddo
      return
      end
      
      DOUBLE PRECISION FUNCTION PL(L,X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C ****  CALCULATES LEGENDRE POLYNOMIAL, ORDER L, ARGUMENT X
C
      PL = 0.0D0
        IF((ABS(X).LE.1.0D0).AND.(L.GE.0)) THEN
        IF(L.GT.5) GO TO 200
          IF(L.LE.1) THEN
            IF(L.EQ.0) THEN
            PL = 1.0D0
            ELSE
            PL = X
            END IF
          ELSE
          X2 = X*X
            IF(L.LE.3) THEN
              IF(L.EQ.2) THEN
              PL = 1.5D0*(X2-0.33333333333333D0)
              ELSE
              PL = 2.5D0*X*(X2-0.6D0)
              END IF
            ELSE
           PL=4.375D0*((X2-0.857142857142857D0)*X2+0.857142857142857D-1)
           IF(L.EQ.4) RETURN
           PL=7.875D0*X*((X2-1.11111111111111D0)*X2+0.238095238095238D0)
            END IF
          END IF
          RETURN
C
C ****  EVALUATE THE LEGENDRE POLYNOMIAL USING A RECURSION FORMULA
C ****  IF ITS ORDER IS TOO LARGE FOR DIRECT EVALUATION.
C
  200   CONTINUE
        X2 = X*X
        P0=4.375D0*((X2-0.857142857142857D0)*X2+0.857142857142857D-1)
        P1=7.875D0*X*((X2-1.11111111111111D0)*X2+0.238095238095238D0)
        DO 100 I = 6,L
        YL = DBLE(I-1)/DBLE(I)
        P2 = X*P1*(1.0D0+YL) - YL*P0
        P0 = P1
        P1 = P2
  100   CONTINUE
        PL = P2
        ELSE
        WRITE(6,1000)L,X
        END IF
      RETURN
 1000 FORMAT(1H0,'ATTEMPT TO FIND P',I2,'(',F5.2,')')
      END
      DOUBLE PRECISION FUNCTION QL(L,X)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C ****  COMPUTE LEGENDRE FUNCTION OF THE SECOND KIND
C ****  FOR ORDER L AND ARGUMENT X.
C
      QL = 0.0D0
        IF((L.LT.0).OR.(ABS(X).EQ.1.0D0)) THEN
        WRITE(6,1000) L,X
        RETURN
        END IF
C
C ****  IF ABS(X).GT.1.07 USE HYPERGEOMETRIC REPRESENTATION
C ****  TO AVOID ROUNDOFF. ABRAMOWITZ AND STEGUN PAGE 332.
C
      IF(ABS(X).GT.1.07D0) GO TO 300
C
C ****  EVALUATE QL DIRECTLY IF L IS SUFFICIENTLY SMALL
C
      A = 0.5D0*LOG((1.0D0+X)/ABS(1.0D0-X))
      X2 = X*X
      IF(L.GE.4) GO TO 100
        IF(L.LE.1) THEN
        QL = A
        IF(L.EQ.1) QL = X*A - 1.0D0
        ELSE
        IF(L.EQ.2) QL = 1.5D0*((X2-0.333333333333333D0)*A - X)
        IF(L.EQ.3) QL = 2.5D0*(X*(X2-0.6D0)*A-X2)+2.0D0/3.0D0
        END IF
      RETURN
C
C ****  USE THE RECURSION RELATION WHEN L >= 4.
C
  100 CONTINUE
      Q0 = 1.5D0*((3.0D0*X2-0.333333333333333D0)*A - X)
      Q1 = 2.5D0*(X*(X2-0.6D0)*A-X2)+2.0D0/3.0D0
      DO 200 LL = 4,L
      DLL = DBLE(LL)
      YL = -1.0D0/DLL
      QL = (2.0D0+YL)*X*Q1 - (1.0D0+YL)*Q0
      Q0 = Q1
      Q1 = QL
  200 CONTINUE
      RETURN
C
C ****  HYPERGEOMETRIC REPRESENTATION
C
  300 CONTINUE
      FAC = 1.0D0/X
        IF(L.GT.0) THEN
        DO 360 LM = 1,L
        FAC = FAC*DBLE(LM)/(DBLE(LM)+0.5D0)
  360   CONTINUE
        FAC = FAC/((2.0D0*X)**L)
        END IF
C
      AP = 0.5D0*DBLE(L+2)
      BP = 0.5D0*DBLE(L+1)
      CP = 0.5D0*DBLE(2*L+3)
      Z = 1.0D0/(X*X)
      QHYPE = 1.0D0
      T = 1.0D0
      DP = 1.0D0
      DO 500 I = 1,300
      T = T*AP*BP*Z/(CP*DP)
      RAT = ABS(T/QHYPE)
      IF(RAT.LE.1D-12.OR.ABS(T).LE.1D-12) GO TO 17
      QHYPE = QHYPE + T
      AP = AP + 1.0D0
      BP = BP + 1.0D0
      CP = CP + 1.0D0
      DP = DP + 1.0D0
  500 CONTINUE
C
C ****  ERROR EXITS
C
      WRITE(6,14) RAT,T
   17 QL = FAC*QHYPE
      RETURN
 1000 FORMAT(1H0,'QL(X) NOT DEFINED FOR L =',I3,' AND X=',F10.3)
   14 FORMAT(' WARNING: QL HAS NOT CONVERGED TO REQUIRED ACCURACY'
     1,'RATIO = ',D16.4,' LAST TERM=',D16.4)
      END
      DOUBLE PRECISION FUNCTION PLM(L,M,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'par.f'
      COMMON/DFCM/DF1(lmax),DF2(lmax)
C
C ****  CALCULATE ASSOCIATED LEGENDRE POLY FROM A RECURSION RELATION
C ****  IN MESSIAH. NOTE: M CAN'T BE NEGATIVE
C
      PLM = 0.0D0
C
C ***  CHECK THAT PARAMETERS ARE O.K.
C 
        IF((ABS(X).GT.1D0).OR.(L.LT.0)) THEN
        WRITE(6,1000)X,L
        STOP
        END IF
        IF((M.LT.0).OR.(M.GT.L)) THEN
        WRITE(6,1100)L,M
        STOP
        END IF
C
C ****  DIRECT EVALUATION WHEN L = 0 OR 1.
C
      IF(L.GT.1) GO TO 200
      PLM = 1.0D0
      IF(L.EQ.0) RETURN
        IF(M.EQ.0) THEN
        PLM = X
        ELSE
        PLM = DSQRT(1.0D0-X*X)
        END IF
      RETURN
C
C ****  HERE L.GE.2
C ****  FIRST, CALCULATE PLM(M,M,X)
C
  200 CONTINUE
C
      if (m.eq.0) then
         omx = 1.0
      else 
         OMX = DSQRT(1.0D0-X*X)**M
      endif 
      A = OMX*DF1(M+1)
        IF(L.LE.M+1) THEN
        PLM = A
        IF(L.EQ.M) RETURN
        B = A*DBLE(2*M+1)*X
        PLM = B
        ELSE
C
        B = A*DBLE(2*M+1)*X
        DM = DBLE(M)
        MPTWO = M + 2
        DO 600 LL = MPTWO,L
        DLL = DBLE(LL)
        PLM = ((2.0D0*DLL-1.0D0)*X*B-(DLL+DM-1.0D0)*A)/(DLL-DM)
        A = B
        B = PLM
  600   CONTINUE
        END IF
      RETURN
 1000 FORMAT(1H0,'ATTEMPT TO FIND PLM(',F5.2,') WITH L =',I3)
 1100 FORMAT(1H0,'ATTEMPT TO FIND PLM WITH L=',I3,' AND M=',I3)
      END

      FUNCTION YLM(L,M,THETA,PHI)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMPLEX*16 EIPHI,ZLM,YLM
      PARAMETER( PI4=12.5663706143591D0 )
      include 'par.f'
      parameter (lll=lmax*5+20)
      COMMON/WYGNER/XXXX,FL(lll)
      COMMON/REDFAK/YNN,FAK(lll)
C
C ****  SPHERICAL HARMONIC FUNCTION OF ORDER L, MAGNETIC PROJECTION M
C ****  WITH ARGUMENTS THETA AND PHI (IN RADIANS)
C
      ZLM = DCMPLX(0.0D0,0.0D0)
      X = COS(THETA)
      MABS = IABS(M)
      if (l+mabs+1.gt.lll) stop 'L too large in YLM'
      PHIM = DBLE(MABS)*PHI
      APHI = ABS(PHIM)
C
        IF((MABS.GT.0).AND.(APHI.GT.1.0D-10)) THEN
        CS = COS(PHIM)
        SS = SIN(PHIM)
        EIPHI = DCMPLX(CS,SS)
        ELSE
        EIPHI = DCMPLX(1.0D0,0.0D0)
        END IF
C
        IF(L.LE.1.) THEN
          IF(L.EQ.0) THEN
          SP = 1.0D0
          ELSE
            IF(M.EQ.0) THEN
            SP = X
            ELSE
            SP = DSQRT(1.0D0-X*X)
            END IF
          END IF
        ELSE
        SP = PLM(L,MABS,X)
        END IF
C
      RL = DBLE(2*L+1)/PI4
      SP = SP*SQRT(RL*FAK(L-MABS+1)/FAK(L+MABS+1))/(YNN**MABS)
C
        IF((M.GE.0).AND.(MOD(MABS,2).NE.0)) THEN
        ZLM = DCMPLX(-SP,0.0D0)
        ELSE
        ZLM = DCMPLX(SP,0.0D0)
        END IF
C
      YLM = ZLM*EIPHI
      RETURN
      END
      SUBROUTINE DFSET
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'par.f'
      COMMON/DFCM/DF1(lmax),DF2(lmax)
C
C ****  INITIALIZATION OF DOUBLE FACTORIAL FUNCTION
C ****  DF1(I+1) = (2*I+1)!!
C ****  DF2(I+1) = (2*I)!!
C
      DF1(1) = 1.0D0
      DF2(1) = 1.0D0
C This patch remains until we have completed
C transferring all programmes onto the Encore.      
Ctemporary-uncomment when everything works    DO 100 I = 2,lmax
C This routine crashes without this patch
C if lmax > 147.      
      DO 100 I = 2,min(140,lmax)
c$$$      DO 100 I = 2,lmax
      DF1(I) = DBLE(2*I-3)*DF1(I-1)
      DF2(I) = DBLE(2*I-2)*DF2(I-1)
  100 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION RYLM(L,M,THETA)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      include 'par.f'
      parameter (lll=lmax*5+20)
      PARAMETER( PI4=12.5663706143591D0 )
      COMMON/REDFAK/YNN,FAK(lll)
C
C ****  SPHERICAL HARMONIC FUNCTION OF ORDER L, MAGNETIC PROJECTION M
C ****  WITH ARGUMENTS THETA (IN RADIANS). THE AZIMUTHAL ANGLE, PHI
C ****  IS ZERO. THE SPHERICAL HARMONIC IS A REAL FUNCTION IN THIS CASE.
C
      SP = 0.0D0
      X = COS(THETA)
      MABS = IABS(M)
      SP = PLM(L,MABS,X)
      RL = DBLE(2*L+1)/PI4
c$$$      if (mabs.gt.1)
c$$$     >   print*,'n,m,fakred:',FAK(L-MABS+1)/FAK(L+MABS+1)/YNN**(2*MABS),
c$$$     >   fakdiv(L-MABS+1,L+MABS+1),L-MABS+1,L+MABS+1
      SP = SP*SQRT(RL*fakdiv(L-MABS+1,L+MABS+1))
c$$$      SP = SP*SQRT(RL*FAK(L-MABS+1)/FAK(L+MABS+1))/(YNN**MABS)
      IF(M.GE.0.AND.MOD(MABS,2).NE.0) SP = -SP
      if (theta.lt.0d0) sp = sp * (-1) ** m
      RYLM = SP
      RETURN
      END
      SUBROUTINE YLMVEC(L,X,Y,Z,YL)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMPLEX*16 EIPHI,EIM
      include 'par.f'
      PARAMETER ( LU=2 , LMAX2=2*LU+1 )
      PARAMETER ( PI4=12.5663706143591D0 )
      parameter (lll=lmax*5+20)
      COMPLEX*16 YL(-LMAX2:LMAX2)
      COMMON/REDFAK/YNN,FAK(lll)
C
C ****  THIS SUBROUTINE CALCULATES THE SPHERICAL HARMONIC VECTOR
C ****  FOR CARTESIAN COORDINATES.
C
      R = SQRT(X*X+Y*Y+Z*Z)
        IF(R.LT.1.0D-14) THEN
        X1 = X*1.0D10
        Y1 = Y*1.0D10
        Z1 = Z*1.0D10
        R = SQRT(X1*X1+Y1*Y1+Z1*Z1)
          IF(R.LT.1.0D-14) THEN
          WRITE(6,1000)
          RETURN
          END IF
        ELSE
        X1 = X
        Y1 = Y
        Z1 = Z
        END IF
C
      C = Z1/R
      RS = R*SQRT(1.0D0-C*C)
C
        IF(RS.GT.1.0D-36) THEN
        CP = X/RS
        SP = Y/RS
        ELSE
        CP = X/RS
        SP = Y/RS
        END IF
C
      FL4 = DBLE(2*L+1)/PI4
      EIPHI = DCMPLX(CP,SP)
      EIM = DCMPLX(1.0D0,0.0D0)
      RPH = 1.0D0
      DO 100 M = 0,L
      F = SQRT(FL4*FAK(L-M+1)/FAK(L+M+1))*PLM(L,M,C)
      YL(M) = F*RPH*EIPHI
      RPH = -RPH
      EIM = EIM*EIPHI
  100 CONTINUE
C
        IF(MOD(L,2).NE.0) THEN
        RPH = -1.0D0
        ELSE
        RPH = 1.0D0
        END IF
      DO 200 M = -L,-1
      YL(M) = RPH*CONJG(YL(-M))
      RPH = -RPH
  200 CONTINUE
      RETURN
 1000 FORMAT(2X,' ZERO RETURNED IN YLMVEC, TOO CLOSE TO ORIGIN')
      END
      DOUBLE PRECISION FUNCTION XPLM(L,M,X)
C
C     CALCULATE ASSOCIATED LEGENDRE POLY FROM A RECURSION RELATION
C     IN MESSIAH. NOTE: M CAN'T BE NEGATIVE
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(ABS(X).GT.1D0) THEN
                             WRITE(6,1000) X
1000   FORMAT(' WARNING : ERROR IN XPLM   X MUST BE BETWEEN [-1,1]',
     *        '   HERE X=',G22.16)
                            STOP
                        END IF
      IF(L.LT.0) WRITE(6,1) L
1     FORMAT(/,' ? L IN XPLM = ',I6,/)
      IF(M.LT.0) WRITE(6,2) M
2     FORMAT(/,' ? M IN XPLM = ',I6,/)
      IF(M.GT.L) WRITE(6,3) L,M
3     FORMAT(/,' ? L = ',I6,' & M = ',I6,' IN XPLM'/)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(L.NE.0) GO TO 4
      XPLM=1.D0
      RETURN
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
4     IF(L.NE.1) GO TO 6
      IF(M.NE.0) GO TO 5
      XPLM=X
      RETURN
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
5     XPLM=DSQRT(1.D0-X*X)
      RETURN
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     HERE  L.GE.2
C
6     IF(M.NE.0) GO TO 7
      A=1.D0
      B=X
      DM=DBLE(M)
      GO TO 10
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     CALCULATE XPLM(M,M,X)
C
7     DF=1.D0
      ITMONE=2*M-1
      DO 8 LL=1,ITMONE,2
      DLL=DBLE(LL)
8     DF=DF*DLL
      if (m.eq.0) then
         omx = 1.0
      else 
         OMX=DSQRT(1.D0-X*X)**M
      endif 
      A=OMX*DF
      IF(L.NE.M) GO TO 9
      XPLM=A
      RETURN
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
9     DM=DBLE(M)
      B=A*(2.0D0*DM+1.D0)*X
      IF(L.NE.M+1) GO TO 10
      XPLM=B
      RETURN
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
10      MPTWO=M+2
       DO 11 LL=MPTWO,L
      DLL=DBLE(LL)
      XPLM=((2.0D0*DLL-1.D0)*X*B-(DLL+DM-1.D0)*A)/(DLL-DM)
      A=B
      B=XPLM
11    CONTINUE
      RETURN
      END

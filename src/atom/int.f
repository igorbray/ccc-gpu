      SUBROUTINE  INTZZZ(R5,R8,RO1,BET,
     *   ALFA,R7,H,R2,IS,NE,NE5,NE5IS,KT,IN,K,r1)
      include 'paratom.f'
      real*8 x(nsizek7),y(nsizek7)
      REAL*8ALFA,BET,H,Z,RO1,RO2,EPS,R2(KT ),
     *   R5(NE5 ),R7( KT),R8(NE5IS),ROT,
     *   AA1,T2,AKSI,UU,VV ,R,r1(ne5)
      INTEGER IN(IS)

      DO 1 JJ=1,IS
         IAA=JJ
c!!!
c
c     write wave functions to file
c
c$$$         open(10,file='wf')
c$$$         open(20,file='en')
         IF(K.EQ.3.AND.JJ.NE.IS)  GOTO 1
         MSK= NE5 *  (IAA-1)
         DO 10 IW=1,NE5
 10      R5(IW)=R8(IW+MSK)
         
         write(20, *)  R5(1)
         write(66,   2) R5(1),R5(NE+3)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         IF(R5(1).LT.0.0)  write(66, 201) IN(JJ)
 201     FORMAT(' Principal quantum number=',6X,I5)
         write(66, 202) R5(2)
 202     FORMAT(' Orbital moment=     ',
     *      9X,     E15.8// )
         write(66, 5) JJ
 5       FORMAT( /  '   Wave function',I3,' shell'
     *      /3('  R',11X,'P(R)',18X))
c$$$         DO 3 IWW=1,KT
c$$$            ROT=ALFA*R7(IWW)+BET* DLOG(R7(IWW))
c$$$            AA1=(ROT-RO1)/H
c$$$            IT1=AA1
c$$$            if (it1.gt.ne-3) then
c$$$               print*,'Warning: IT1 reset to NE - 3, may need to'
c$$$     >            //' increase NE'
c$$$               it1 = ne - 3
c$$$            end if 
c$$$            T2=AA1-IT1
c$$$            IF(IT1.EQ.0) GOTO 3
c$$$            IF(IT1.GT.NE) GOTO 1
c$$$            AKSI=1.-T2
c$$$            UU=T2*(1.-T2**2)/6.
c$$$            VV=AKSI*(1.-AKSI**2)/6.
c$$$            R2(IWW)=DSQRT(R7(IWW)/(ALFA*
c$$$     *         R7(IWW)+BET))*(R5(IT1+2)*AKSI
c$$$     *         +R5(IT1+3)*T2-VV*(R5(IT1+3)-2.
c$$$     *         *R5(IT1+2)+R5(IT1+1))-
c$$$     *         UU*(R5(IT1
c$$$     *         +4)-2.*R5(IT1+3)+R5(IT1+2)))
c$$$ 3       CONTINUE

         do i = 1, ne
            x(i) = r1(2+i)
            y(i) = R8(2+I+MSK) * sqrt(x(i)/(alfa * x(i) + bet))
         end do 
         call INTRPL(ne,x,y,kt,r7,r2)
c$$$         write(66, 4) (R7(IW),R2(IW),IW=1,KT )
         write(10) (R7(IW),R2(IW),IW=1,KT )
c$$$         write(10, 40) (R7(IW),R2(IW),IW=1,KT )
 4       FORMAT(3(2X,F8.5,2X,E20.13,4X))
         
 40      FORMAT(2X,F8.5,2X,F20.13,4X)
         write(66, 14)
 14      FORMAT(/)
 1    CONTINUE
c
c$$$      close(10)
c$$$      close(20)
c
      RETURN
      END

C  FILE 'INTRPL.FTN77'
C
      SUBROUTINE INTRPLOLD(L,X,Y,N,U,V)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      include 'paratom.f'
      PARAMETER ( NR=NSIZEK7 )
C
C     DOUBLE PRECISION INTERPOLATION OF A SINGLE VALUED FUNCTION
C     THIS SUBROUTINE INTERPOLATES, FROM VALUES OF THE FUNCTION
C     GIVEN  AS ORDINATES OF INPUT DATA POINTS IN AN X-Y PLANE
C     AND FOR A GIVEN SET OF X VALUES(ABSCISSAE),THE VALUES OF
C     A SINGLE VALUED FUNCTION Y=Y(X).
C
C     THE SUBROUTINE ALSO CALCULATES FIRST DERIVATIVES DV(X)
C     AND THE INTEGRAL OVER THE WHOLE RANGE OF X FINT.
C     DV(X) IS IN COMMON/DV/, CALCULATED IF IDERIV.NE.0
C     FINT IS IN COMMON/INTEG/, CALCULATED IF INT.NE.0
C     QQ(4,L) ARE THE SPLINE COEFFICIENTS. THEY ARE STORED IN COMMON/QQ/
C
C     THE INPUT PARAMETERS ARE;
C
C     L=NUMBER OF DATA POINTS
C     (MUST BE TWO OR GREATER)
C     X=ARRAY OF DIMENSION L STORING THE X VALUES
C     OF INPUT DATA POINTS (IN ASCENDING ORDER)
C     Y=ARRAY OF DIMENSION L STORING THE Y VALUES OF INPUT DATA POINTS
C     N=NUMBER OF POINTS AT WHICH INTERPOLATION OF THE Y-VALUES
C     IS REQUIRED (MUST BE 1 OR GREATER)
C     U=ARRAY OF DIMENSION N STORING THE X VALUES
C     OF THE DESIRED POINTS
C
C     THE OUTPUT PARAMETER IS V=ARRAY OF DIMENSION N WHERE THE
C     INTERPOLATED Y VALUES ARE TO BE DISPLAYED
C
      COMMON/DV/ IDERIV,DV(NR)
      COMMON/INTEG/ INT,FINT
      COMMON/QQ/ QQ(4,NR)
      DIMENSION X(L),Y(L),U(N),V(N)
      EQUIVALENCE (P0,X3),(Q0,Y3),(Q1,T3)
      REAL*8 M1,M2,M3,M4,M5
      EQUIVALENCE (UK,DX),(IMN,X2,A1,M1),(IMX,X5,A5,M5),
     1         (J,SW,SA),(Y2,W2,W4,Q2),(Y5,W3,Q3)
C
C     PRELIMINARY PROCESSING
C
      L0 = L
      LM1 = L0 - 1
      LM2 = LM1 - 1
      LP1 = L0 + 1
      N0 = N
      NQQ = NR
        IF(N0.GT.NQQ) THEN
        NQQV = NQQ
        WRITE(6,2089) NQQV,N0
        STOP
        END IF
      IF(LM2.LT.0) GO TO 90
      IF(N0.LE.0) GO TO 91
      DO 11 I = 2,L0
      IF(X(I-1)-X(I))11,95,96
   11 CONTINUE
      IPV = 0
C
C     MAIN LOOP
C
      FINT = 0.0D0
      DO 80 K = 1,N0
      UK = U(K)
C
C     BINARY SEARCH TO LOCATE THE DESIRED POINT
C
      IF(LM2.EQ.0) GO TO 27
      IF(UK.GE.X(L0)) GO TO 26
      IF(UK.LT.X(1)) GO TO 25
      IMN = 2
      IMX = L0
   21 I = (IMN+IMX)/2
      IF(UK.GE.X(I)) GO TO 23
      IMX = I
      GO TO 24
   23 IMN = I + 1
   24 IF(IMX.GT.IMN) GO TO 21
      I = IMX
      GO TO 30
   25 I = 1
      GO TO 30
   26 I = LP1
      GO TO 30
   27 I = 2
C
C     CHECK IF I = IPV
C
   30 IF(I.EQ.IPV) GO TO 70
      IPV = I
C
C     ROUTINES TO PICK UP NECESSARY X AND Y VALUES AND TO
C     ESTIMATE THEM IF NECESSARY
C
      J = I
      IF(J.EQ.1) J = 2
      IF(J.EQ.LP1) J = L0
      X3 = X(J-1)
      Y3 = Y(J-1)
      X4 = X(J)
      Y4 = Y(J)
      A3 = X4 - X3
      M3 = (Y4 - Y3)/A3
      IF(LM2.EQ.0) GO TO 43
      IF(J.EQ.2) GO TO 41
      X2 = X(J-2)
      Y2 = Y(J-2)
      A2 = X3 - X2
      M2 = (Y3 - Y2)/A2
      IF(J.EQ.L0) GO TO 42
   41 X5 = X(J+1)
      Y5 = Y(J+1)
      A4 = X5 - X4
      M4 = (Y5 - Y4)/A4
      IF(J.EQ.2) M2 = M3 + M3 - M4
      GO TO 45
   42 M4 = M3 + M3 - M2
      GO TO 45
   43 M2 = M3
   45 IF(J.LE.3) GO TO 46
      A1 = X2 - X(J-3)
      M1 = (Y2 - Y(J-3))/A1
      GO TO 47
   46 M1 = M2 + M2 - M3
   47 IF(J.GE.LM1) GO TO 48
      A5 = X(J+2) - X5
      M5 = (Y(J+2)-Y5)/A5
      GO TO 50
   48 M5 = M4 + M4 - M3
C
C     NUMERICAL DIFFERENTIATION
C
   50 IF(I.EQ.LP1) GO TO 52
      W2 = ABS(M4-M3)
      W3 = ABS(M2-M1)
      SW = W2 + W3
      IF(SW.NE.0.0D0) GO TO 51
      W2 = 0.5D0
      W3 = 0.5D0
      SW = 1.0D0
   51 T3 = (W2*M2 + W3*M3)/SW
      IF(I.EQ.1) GO TO 54
   52 W3 = ABS(M5-M4)
      W4 = ABS(M3-M2)
      SW = W3 + W4
      IF(SW.NE.0.0D0) GO TO 53
      W3 = 0.5D0
      W4 = 0.5D0
      SW = 1.0D0
   53 T4 = (W3*M3 + W4*M4)/SW
      IF(I.NE.LP1) GO TO 60
      T3 = T4
      SA = A2 + A3
      T4 = 0.5D0*(M4 + M5 - A2*(A2-A3)*(M2-M3)/(SA*SA))
      X3 = X4
      Y3 = Y4
      A3 = A2
      M3 = M4
      GO TO 60
   54 T4 = T3
      SA = A3 + A4
      T3 = 0.5D0*(M1 + M2 - A4*(A3-A4)*(M3-M4)/(SA*SA))
      X3 = X3 - A4
      Y3 = Y3 - M2*A4
      A3 = A4
      M3 = M2
C
C     DETERMINATION OF THE COEFFICIENTS
C
   60 Q2 = (2.0D0*(M3-T3) + M3 - T4)/A3
      Q3 = (-M3 - M3 + T3 + T4)/(A3*A3)
C
C     COMPUTATION OF THE POLYNOMIAL
C
   70 DX = UK - P0
      V(K) = Q0 + DX*(Q1 + DX*(Q2 + DX*Q3))
      IF(IDERIV.EQ.0) GO TO 72
      DV(K) = Q1 + DX*(2.0D0*Q2 + DX*3.0D0*Q3)
      QQ(1,K) = Q0
      QQ(2,K) = Q1
      QQ(3,K) = Q2
      QQ(4,K) = Q3
   72 IF(INT.EQ.0) GO TO 80
      IF(K.LT.N0) DD = U(K+1) - U(K)
      IF(K.LT.N0) FINT = FINT
     1 + DD*(Q0 + DD*(Q1/2.0D0+DD*(Q2/3.0D0+DD*Q3/4.0D0)))
   80 CONTINUE
      RETURN
C
C     ERROR EXITS
C
   90 WRITE(6,2090)
      GO TO 99
   91 WRITE(6,2091)
      GO TO 99
   95 WRITE(6,2095)
      GO TO 97
   96 WRITE(6,2096)
   97 WRITE(6,2097)I,X(I)
   99 WRITE(6,2099) L0,N0
      RETURN
C
C     FORMAT STATEMENTS
C
 2089 FORMAT( 'WARNING ERROR IN INTRPL. MAX ALLOWED VALUE OF N0 IS',
     *I5,' HERE N0 IS',I5)
 2090 FORMAT(1X/' N  =  1 OR LESS.'/)
 2091 FORMAT(1X/' N  =  0 OR LESS.'/)
 2095 FORMAT(1X/' IDENTICAL X VALUES.'/)
 2096 FORMAT(1X/' X VALUES OUT OF SEQUENCE.'/)
 2097 FORMAT(6H   I =,I7,10X,6HX(I) =,E12.3)
 2099 FORMAT(6H   L =,I7,10X,3HN =,I7/
     1 ' ERRROR DETECTED IN ROUTINE INTRPL')
      END
cx      BLOCK DATA
cx      IMPLICIT REAL*8(A-H,O-Z)
cx      IMPLICIT INTEGER*4(I-N)
cx      include 'par.for'
c          PARAMETER ( NR=1000 )
cx      COMMON/DV/ IDERIV,DV(NR)
cx      COMMON/INTEG/ INT,FINT
cx      DATA IDERIV,INT/1,0/
cx      END
C     ALGORITHM 474 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN COMM. ACM, VOL. VV, NO. NN,
C     P. 000.
      SUBROUTINE ITPLBV(IU, LX, LY, X, Y, Z, N, U, V, W)                   A  10
C BIVARIATE INTERPOLATION
C THIS SUBROUTINE INTERPOLATES, FROM VALUES OF THE FUNCTION
C GIVEN AT INPUT GRID POINTS IN AN X-Y PLANE AND FOR A GIVEN
C SET OF POINTS IN THE PLANE, THE VALUES OF A SINGLE-VALUED
C BIVARIATE FUNCTION Z = Z(X,Y).
C THE METHOD IS BASED ON A PIECE-WISE FUNCTION COMPOSED OF
C A SET OF BICUBIC POLYNOMIALS IN X AND Y.  EACH POLYNOMIAL
C IS APPLICABLE TO A RECTANGLE OF THE INPUT GRID IN THE X-Y
C PLANE.  EACH POLYNOMIAL IS DETERMINED LOCALLY.
C THE INPUT PARAMETERS ARE
C IU  = LOGICAL UNIT NUMBER OF STANDARD OUTPUT UNIT
C LX  = NUMBER OF INPUT GRID POINTS IN THE X COORDINATE
C       (MUST BE 2 OR GREATER)
C LY  = NUMBER OF INPUT GRID POINTS IN THE Y COORDINATE
C       (MUST BE 2 OR GREATER)
C X   = ARRAY OF DIMENSION LX STORING THE X COORDINATES
C       OF INPUT GRID POINTS (IN ASCENDING ORDER)
C Y   = ARRAY OF DIMENSION LY STORING THE Y COORDINATES
C       OF INPUT GRID POINTS (IN ASCENDING ORDER)
C Z   = DOUBLY-DIMENSIONED ARRAY OF DIMENSION (LX,LY)
C       STORING THE VALUES OF THE FUNCTION (Z VALUES)
C       AT INPUT GRID POINTS
C N   = NUMBER OF POINTS AT WHICH INTERPOLATION OF THE
C       Z VALUE IS DESIRED (MUST BE 1 OR GREATER)
C U   = ARRAY OF DIMENSION N STORING THE X COORDINATES
C       OF DESIRED POINTS
C V   = ARRAY OF DIMENSION N STORING THE Y COORDINATES
C       OF DESIRED POINTS
C THE OUTPUT PARAMETER IS
C W   = ARRAY OF DIMENSION N WHERE THE INTERPOLATED Z
C       VALUES AT DESIRED POINTS ARE TO BE DISPLAYED
C       SOME VARIABLES INTERNALLY USED ARE
C ZA  = DIVIDED DIFFERENCE OF Z WITH RESPECT TO X
C ZB  = DIVIDED DIFFERENCE OF Z WITH RESPECT TO Y
C ZAB = SECOND ORDER DIVIDED DIFFERENCE OF Z WITH
C       RESPECT TO X AND Y
C ZX  = PARTIAL DERIVATIVE OF Z WITH RESPECT TO X
C ZY  = PARTIAL DERIVATIVE OF Z WITH RESPECT TO Y
C ZXY = SECOND ORDER PARTIAL DERIVATIVE OF Z WITH
C       RESPECT TO X AND Y
C DECLARATION STATEMENTS
      DIMENSION X(LX), Y(LY), Z(LX,LY), U(N), V(N), W(N)
      DIMENSION ZA(5,2), ZB(2,5), ZAB(3,3), ZX(4,4), ZY(4,4),
     * ZXY(4,4)
      EQUIVALENCE (Z3A1,ZA(1,1)), (Z3A2,ZA(2,1)), (Z3A3,ZA(3,1)),
     * (Z3A4,ZA(4,1)), (Z3A5,ZA(5,1)), (Z4A1,ZA(1,2)), (Z4A2,ZA(2,2)),
     * (Z4A3,ZA(3,2)), (Z4A4,ZA(4,2)), (Z4A5,ZA(5,2)), (Z3B1,ZB(1,1)),
     * (Z3B2,ZB(1,2)), (Z3B3,ZB(1,3)), (Z3B4,ZB(1,4)), (Z3B5,ZB(1,5)),
     * (Z4B1,ZB(2,1)), (Z4B2,ZB(2,2)), (Z4B3,ZB(2,3)), (Z4B4,ZB(2,4)),
     * (Z4B5,ZB(2,5)), (ZA2B2,ZAB(1,1)), (ZA3B2,ZAB(2,1)),
     * (ZA4B2,ZAB(3,1)), (ZA2B3,ZAB(1,2)), (ZA3B3,ZAB(2,2)),
     * (ZA4B3,ZAB(3,2)), (ZA2B4,ZAB(1,3)), (ZA3B4,ZAB(2,3)),
     * (ZA4B4,ZAB(3,3)), (ZX33,ZX(2,2)), (ZX43,ZX(3,2)),
     * (ZX34,ZX(2,3)), (ZX44,ZX(3,3)), (ZY33,ZY(2,2)),
     * (ZY43,ZY(3,2)), (ZY34,ZY(2,3)), (ZY44,ZY(3,3)),
     * (ZXY33,ZXY(2,2)), (ZXY43,ZXY(3,2)), (ZXY34,ZXY(3,2)),
     * (ZXY44,ZXY(3,3)), (P00,Z33), (P01,ZY33), (P10,ZX33),
     * (P11,ZXY33)
      EQUIVALENCE (LX0,ZX(1,1)), (LXM1,ZX(4,1)), (LXM2,ZX(1,4)),
     * (LXP1,ZX(4,4)), (LY0,ZY(1,1)), (LYM1,ZY(4,1)), (LYM2,ZY(1,4)),
     * (LYP1,ZY(4,4)), (IX,ZXY(1,1)), (IY,ZXY(4,1)), (IXPV,ZXY(1,4)),
     * (IYPV,ZXY(4,4)), (IMN,JX), (IMX,JY), (JXM2,JX1),
     * (JYM2,JY1), (UK,DX), (VK,DY), (A1,A5,B1,B5,ZX(2,1),A,Q0),
     * (A2,ZX(1,2),B,Q1), (A4,ZX(4,2),C,Q2), (B2,ZY(2,1),D,Q3),
     * (B4,ZY(2,4),E), (X2,ZX(3,1),A3SQ), (X4,ZX(1,3)), (X5,ZX(4,3)),
     * (Y2,ZX(2,4)), (Y4,ZY(3,1),B3SQ), (Y5,ZX(3,4),P02),
     * (Z23,ZY(1,2),P03), (Z24,ZY(4,2),P12), (Z32,ZY(1,3),P13),
     * (Z34,ZY(4,3),P20), (Z35,ZY(3,4),P21), (Z42,ZXY(2,1),P22),
     * (Z43,ZXY(1,2),P23), (Z44,ZXY(3,1),P30), (Z45,ZXY(4,2),P31),
     * (Z53,ZXY(1,3),P32), (Z54,ZXY(4,3),P33), (W2,WY2,W4),
     * (W3,WY3,W1,W5), (WX2,ZXY(2,4)), (WX3,ZXY(3,4))
C PRELIMINARY PROCESSING
C SETTING OF SOME INPUT PARAMETERS TO LOCAL VARIABLES
      IU0 = IU
      LX0 = LX
      LXM1 = LX0 - 1
      LXM2 = LXM1 - 1
      LXP1 = LX0 + 1
      LY0 = LY
      LYM1 = LY0 - 1
      LYM2 = LYM1 - 1
      LYP1 = LY0 + 1
      N0 = N
C ERROR CHECK
      IF (LXM2.LT.0) GO TO 710
      IF (LYM2.LT.0) GO TO 720
      IF (N0.LT.1) GO TO 730
      DO 10 IX=2,LX0
        IF (X(IX-1)-X(IX)) 10, 740, 750
   10 CONTINUE
      DO 20 IY=2,LY0
        IF (Y(IY-1)-Y(IY)) 20, 770, 780
   20 CONTINUE
C INITIAL SETTING OF PREVIOUS VALUES OF IX AND IY
      IXPV = 0
      IYPV = 0
C MAIN DO-LOOP
      DO 700 K=1,N0
        UK = U(K)
        VK = V(K)
C ROUTINES TO LOCATE THE DESIRED POINT
C TO FIND OUT THE IX VALUE FOR WHICH
C (U(K).GE.X(IX-1)).AND.(U(K).LT.X(IX))
        IF (LXM2.EQ.0) GO TO 80
        IF (UK.GE.X(LX0)) GO TO 70
        IF (UK.LT.X(1)) GO TO 60
        IMN = 2
        IMX = LX0
   30   IX = (IMN+IMX)/2
        IF (UK.GE.X(IX)) GO TO 40
        IMX = IX
        GO TO 50
   40   IMN = IX + 1
   50   IF (IMX.GT.IMN) GO TO 30
        IX = IMX
        GO TO 90
   60   IX = 1
        GO TO 90
   70   IX = LXP1
        GO TO 90
   80   IX = 2
C TO FIND OUT THE IY VALUE FOR WHICH
C (V(K).GE.Y(IY-1)).AND.(V(K).LT.Y(IY))
   90   IF (LYM2.EQ.0) GO TO 150
        IF (VK.GE.Y(LY0)) GO TO 140
        IF (VK.LT.Y(1)) GO TO 130
        IMN = 2
        IMX = LY0
  100   IY = (IMN+IMX)/2
        IF (VK.GE.Y(IY)) GO TO 110
        IMX = IY
        GO TO 120
  110   IMN = IY + 1
  120   IF (IMX.GT.IMN) GO TO 100
        IY = IMX
        GO TO 160
  130   IY = 1
        GO TO 160
  140   IY = LYP1
        GO TO 160
  150   IY = 2
C TO CHECK IF THE DESIRED POINT IS IN THE SAME RECTANGLE
C AS THE PREVIOUS POINT.  IF YES, SKIP TO THE COMPUTATION
C OF THE POLYNOMIAL
  160   IF (IX.EQ.IXPV .AND. IY.EQ.IYPV) GO TO 690
        IXPV = IX
        IYPV = IY
C ROUTINES TO PICK UP NECESSARY X, Y, AND Z VALUES, TO
C COMPUTE THE ZA, ZB, AND ZAB VALUES, AND TO ESTIMATE THEM
C WHEN NECESSARY
        JX = IX
        IF (JX.EQ.1) JX = 2
        IF (JX.EQ.LXP1) JX = LX0
        JY = IY
        IF (JY.EQ.1) JY = 2
        IF (JY.EQ.LYP1) JY = LY0
        JXM2 = JX - 2
        JXML = JX - LX0
        JYM2 = JY - 2
        JYML = JY - LY0
C IN THE CORE AREA, I.E., IN THE RECTANGLE THAT CONTAINS
C THE DESIRED POINT
        X3 = X(JX-1)
        X4 = X(JX)
        A3 = 1.0/(X4-X3)
        Y3 = Y(JY-1)
        Y4 = Y(JY)
        B3 = 1.0/(Y4-Y3)
        Z33 = Z(JX-1,JY-1)
        Z43 = Z(JX,JY-1)
        Z34 = Z(JX-1,JY)
        Z44 = Z(JX,JY)
        Z3A3 = (Z43-Z33)*A3
        Z4A3 = (Z44-Z34)*A3
        Z3B3 = (Z34-Z33)*B3
        Z4B3 = (Z44-Z43)*B3
        ZA3B3 = (Z4B3-Z3B3)*A3
C IN THE X DIRECTION
        IF (LXM2.EQ.0) GO TO 230
        IF (JXM2.EQ.0) GO TO 170
        X2 = X(JX-2)
        A2 = 1.0/(X3-X2)
        Z23 = Z(JX-2,JY-1)
        Z24 = Z(JX-2,JY)
        Z3A2 = (Z33-Z23)*A2
        Z4A2 = (Z34-Z24)*A2
        IF (JXML.EQ.0) GO TO 180
  170   X5 = X(JX+1)
        A4 = 1.0/(X5-X4)
        Z53 = Z(JX+1,JY-1)
        Z54 = Z(JX+1,JY)
        Z3A4 = (Z53-Z43)*A4
        Z4A4 = (Z54-Z44)*A4
        IF (JXM2.NE.0) GO TO 190
        Z3A2 = Z3A3 + Z3A3 - Z3A4
        Z4A2 = Z4A3 + Z4A3 - Z4A4
        GO TO 190
  180   Z3A4 = Z3A3 + Z3A3 - Z3A2
        Z4A4 = Z4A3 + Z4A3 - Z4A2
  190   ZA2B3 = (Z4A2-Z3A2)*B3
        ZA4B3 = (Z4A4-Z3A4)*B3
        IF (JX.LE.3) GO TO 200
        A1 = 1.0/(X2-X(JX-3))
        Z3A1 = (Z23-Z(JX-3,JY-1))*A1
        Z4A1 = (Z24-Z(JX-3,JY))*A1
        GO TO 210
  200   Z3A1 = Z3A2 + Z3A2 - Z3A3
        Z4A1 = Z4A2 + Z4A2 - Z4A3
  210   IF (JX.GE.LXM1) GO TO 220
        A5 = 1.0/(X(JX+2)-X5)
        Z3A5 = (Z(JX+2,JY-1)-Z53)*A5
        Z4A5 = (Z(JX+2,JY)-Z54)*A5
        GO TO 240
  220   Z3A5 = Z3A4 + Z3A4 - Z3A3
        Z4A5 = Z4A4 + Z4A4 - Z4A3
        GO TO 240
  230   Z3A2 = Z3A3
        Z4A2 = Z4A3
        GO TO 180
C IN THE Y DIRECTION
  240   IF (LYM2.EQ.0) GO TO 310
        IF (JYM2.EQ.0) GO TO 250
        Y2 = Y(JY-2)
        B2 = 1.0/(Y3-Y2)
        Z32 = Z(JX-1,JY-2)
        Z42 = Z(JX,JY-2)
        Z3B2 = (Z33-Z32)*B2
        Z4B2 = (Z43-Z42)*B2
        IF (JYML.EQ.0) GO TO 260
  250   Y5 = Y(JY+1)
        B4 = 1.0/(Y5-Y4)
        Z35 = Z(JX-1,JY+1)
        Z45 = Z(JX,JY+1)
        Z3B4 = (Z35-Z34)*B4
        Z4B4 = (Z45-Z44)*B4
        IF (JYM2.NE.0) GO TO 270
        Z3B2 = Z3B3 + Z3B3 - Z3B4
        Z4B2 = Z4B3 + Z4B3 - Z4B4
        GO TO 270
  260   Z3B4 = Z3B3 + Z3B3 - Z3B2
        Z4B4 = Z4B3 + Z4B3 - Z4B2
  270   ZA3B2 = (Z4B2-Z3B2)*A3
        ZA3B4 = (Z4B4-Z3B4)*A3
        IF (JY.LE.3) GO TO 280
        B1 = 1.0/(Y2-Y(JY-3))
        Z3B1 = (Z32-Z(JX-1,JY-3))*B1
        Z4B1 = (Z42-Z(JX,JY-3))*B1
        GO TO 290
  280   Z3B1 = Z3B2 + Z3B2 - Z3B3
        Z4B1 = Z4B2 + Z4B2 - Z4B3
  290   IF (JY.GE.LYM1) GO TO 300
        B5 = 1.0/(Y(JY+2)-Y5)
        Z3B5 = (Z(JX-1,JY+2)-Z35)*B5
        Z4B5 = (Z(JX,JY+2)-Z45)*B5
        GO TO 320
  300   Z3B5 = Z3B4 + Z3B4 - Z3B3
        Z4B5 = Z4B4 + Z4B4 - Z4B3
        GO TO 320
  310   Z3B2 = Z3B3
        Z4B2 = Z4B3
        GO TO 260
C IN THE DIAGONAL DIRECTIONS
  320   IF (LXM2.EQ.0) GO TO 400
        IF (LYM2.EQ.0) GO TO 410
        IF (JXML.EQ.0) GO TO 350
        IF (JYM2.EQ.0) GO TO 330
        ZA4B2 = ((Z53-Z(JX+1,JY-2))*B2-Z4B2)*A4
        IF (JYML.EQ.0) GO TO 340
  330   ZA4B4 = ((Z(JX+1,JY+1)-Z54)*B4-Z4B4)*A4
        IF (JYM2.NE.0) GO TO 380
        ZA4B2 = ZA4B3 + ZA4B3 - ZA4B4
        GO TO 380
  340   ZA4B4 = ZA4B3 + ZA4B3 - ZA4B2
        GO TO 380
  350   IF (JYM2.EQ.0) GO TO 360
        ZA2B2 = (Z3B2-(Z23-Z(JX-2,JY-2))*B2)*A2
        IF (JYML.EQ.0) GO TO 370
  360   ZA2B4 = (Z3B4-(Z(JX-2,JY+1)-Z24)*B4)*A2
        IF (JYM2.NE.0) GO TO 390
        ZA2B2 = ZA2B3 + ZA2B3 - ZA2B4
        GO TO 390
  370   ZA2B4 = ZA2B3 + ZA2B3 - ZA2B2
        GO TO 390
  380   IF (JXM2.NE.0) GO TO 350
        ZA2B2 = ZA3B2 + ZA3B2 - ZA4B2
        ZA2B4 = ZA3B4 + ZA3B4 - ZA4B4
        GO TO 420
  390   IF (JXML.NE.0) GO TO 420
        ZA4B2 = ZA3B2 + ZA3B2 - ZA2B2
        ZA4B4 = ZA3B4 + ZA3B4 - ZA2B4
        GO TO 420
  400   ZA2B2 = ZA3B2
        ZA4B2 = ZA3B2
        ZA2B4 = ZA3B4
        ZA4B4 = ZA3B4
        GO TO 420
  410   ZA2B2 = ZA2B3
        ZA2B4 = ZA2B3
        ZA4B2 = ZA4B3
        ZA4B4 = ZA4B3
C NUMERICAL DIFFERENTIATION   ---   TO DETERMINE PARTIAL
C DERIVATIVES ZX, ZY, AND ZXY AS WEIGHTED MEANS OF DIVIDED
C DIFFERENCES ZA, ZB, AND ZAB, RESPECTIVELY
  420   DO 480 JY=2,3
          DO 470 JX=2,3
            W2 = ABS(ZA(JX+2,JY-1)-ZA(JX+1,JY-1))
            W3 = ABS(ZA(JX,JY-1)-ZA(JX-1,JY-1))
            SW = W2 + W3
            IF (SW.EQ.0.0) GO TO 430
            WX2 = W2/SW
            WX3 = W3/SW
            GO TO 440
  430       WX2 = 0.5
            WX3 = 0.5
  440       ZX(JX,JY) = WX2*ZA(JX,JY-1) + WX3*ZA(JX+1,JY-1)
            W2 = ABS(ZB(JX-1,JY+2)-ZB(JX-1,JY+1))
            W3 = ABS(ZB(JX-1,JY)-ZB(JX-1,JY-1))
            SW = W2 + W3
            IF (SW.EQ.0.0) GO TO 450
            WY2 = W2/SW
            WY3 = W3/SW
            GO TO 460
  450       WY2 = 0.5
            WY3 = 0.5
  460       ZY(JX,JY) = WY2*ZB(JX-1,JY) + WY3*ZB(JX-1,JY+1)
            ZXY(JX,JY) =
     *      WY2*(WX2*ZAB(JX-1,JY-1)+WX3*ZAB(JX,JY-1)) +
     *      WY3*(WX2*ZAB(JX-1,JY)+WX3*ZAB(JX,JY))
  470     CONTINUE
  480   CONTINUE
C WHEN (U(K).LT.X(1)).OR.(U(K).GT.X(LX))
        IF (IX.EQ.LXP1) GO TO 530
        IF (IX.NE.1) GO TO 590
        W2 = A4*(3.0*A3+A4)
        W1 = 2.0*A3*(A3-A4) + W2
        DO 500 JY=2,3
          ZX(1,JY) = (W1*ZA(1,JY-1)+W2*ZA(2,JY-1))/(W1+W2)
          ZY(1,JY) = ZY(2,JY) + ZY(2,JY) - ZY(3,JY)
          ZXY(1,JY) = ZXY(2,JY) + ZXY(2,JY) - ZXY(3,JY)
          DO 490 JX1=2,3
            JX = 5 - JX1
            ZX(JX,JY) = ZX(JX-1,JY)
            ZY(JX,JY) = ZY(JX-1,JY)
            ZXY(JX,JY) = ZXY(JX-1,JY)
  490     CONTINUE
  500   CONTINUE
        X3 = X3 - 1.0/A4
        Z33 = Z33 - Z3A2/A4
        DO 510 JY=1,5
          ZB(2,JY) = ZB(1,JY)
  510   CONTINUE
        DO 520 JY=2,4
          ZB(1,JY) = ZB(1,JY) - ZAB(1,JY-1)/A4
  520   CONTINUE
        A3 = A4
        JX = 1
        GO TO 570
  530   W4 = A2*(3.0*A3+A2)
        W5 = 2.0*A3*(A3-A2) + W4
        DO 550 JY=2,3
          ZX(4,JY) = (W4*ZA(4,JY-1)+W5*ZA(5,JY-1))/(W4+W5)
          ZY(4,JY) = ZY(3,JY) + ZY(3,JY) - ZY(2,JY)
          ZXY(4,JY) = ZXY(3,JY) + ZXY(3,JY) - ZXY(2,JY)
          DO 540 JX=2,3
            ZX(JX,JY) = ZX(JX+1,JY)
            ZY(JX,JY) = ZY(JX+1,JY)
            ZXY(JX,JY) = ZXY(JX+1,JY)
  540     CONTINUE
  550   CONTINUE
        X3 = X4
        Z33 = Z43
        DO 560 JY=1,5
          ZB(1,JY) = ZB(2,JY)
  560   CONTINUE
        A3 = A2
        JX = 3
  570   ZA(3,1) = ZA(JX+1,1)
        DO 580 JY=1,3
          ZAB(2,JY) = ZAB(JX,JY)
  580   CONTINUE
C WHEN (V(K).LT.Y(1)).OR.(V(K).GT.Y(LY))
  590   IF (IY.EQ.LYP1) GO TO 630
        IF (IY.NE.1) GO TO 680
        W2 = B4*(3.0*B3+B4)
        W1 = 2.0*B3*(B3-B4) + W2
        DO 620 JX=2,3
          IF (JX.EQ.3 .AND. IX.EQ.LXP1) GO TO 600
          IF (JX.EQ.2 .AND. IX.EQ.1) GO TO 600
          ZY(JX,1) = (W1*ZB(JX-1,1)+W2*ZB(JX-1,2))/(W1+W2)
          ZX(JX,1) = ZX(JX,2) + ZX(JX,2) - ZX(JX,3)
          ZXY(JX,1) = ZXY(JX,2) + ZXY(JX,2) - ZXY(JX,3)
  600     DO 610 JY1=2,3
            JY = 5 - JY1
            ZY(JX,JY) = ZY(JX,JY-1)
            ZX(JX,JY) = ZX(JX,JY-1)
            ZXY(JX,JY) = ZXY(JX,JY-1)
  610     CONTINUE
  620   CONTINUE
        Y3 = Y3 - 1.0/B4
        Z33 = Z33 - Z3B2/B4
        Z3A3 = Z3A3 - ZA3B2/B4
        Z3B3 = Z3B2
        ZA3B3 = ZA3B2
        B3 = B4
        GO TO 670
  630   W4 = B2*(3.0*B3+B2)
        W5 = 2.0*B3*(B3-B2) + W4
        DO 660 JX=2,3
          IF (JX.EQ.3 .AND. IX.EQ.LXP1) GO TO 640
          IF (JX.EQ.2 .AND. IX.EQ.1) GO TO 640
          ZY(JX,4) = (W4*ZB(JX-1,4)+W5*ZB(JX-1,5))/(W4+W5)
          ZX(JX,4) = ZX(JX,3) + ZX(JX,3) - ZX(JX,2)
          ZXY(JX,4) = ZXY(JX,3) + ZXY(JX,3) - ZXY(JX,2)
  640     DO 650 JY=2,3
            ZY(JX,JY) = ZY(JX,JY+1)
            ZX(JX,JY) = ZX(JX,JY+1)
            ZXY(JX,JY) = ZXY(JX,JY+1)
  650     CONTINUE
  660   CONTINUE
        Y3 = Y4
        Z33 = Z33 + Z3B3/B3
        Z3A3 = Z3A3 + ZA3B3/B3
        Z3B3 = Z3B4
        ZA3B3 = ZA3B4
        B3 = B2
  670   IF (IX.NE.1 .AND. IX.NE.LXP1) GO TO 680
        JX = IX/LXP1 + 2
        JX1 = 5 - JX
        JY = IY/LYP1 + 2
        JY1 = 5 - JY
        ZX(JX,JY) = ZX(JX1,JY) + ZX(JX,JY1) - ZX(JX1,JY1)
        ZY(JX,JY) = ZY(JX1,JY) + ZY(JX,JY1) - ZY(JX1,JY1)
        ZXY(JX,JY) = ZXY(JX1,JY) + ZXY(JX,JY1) - ZXY(JX1,JY1)
C DETERMINATION OF THE COEFFICIENTS OF THE POLYNOMIAL
  680   ZX3B3 = (ZX34-ZX33)*B3
        ZX4B3 = (ZX44-ZX43)*B3
        ZY3A3 = (ZY43-ZY33)*A3
        ZY4A3 = (ZY44-ZY34)*A3
        A = ZA3B3 - ZX3B3 - ZY3A3 + ZXY33
        B = ZX4B3 - ZX3B3 - ZXY43 + ZXY33
        C = ZY4A3 - ZY3A3 - ZXY34 + ZXY33
        D = ZXY44 - ZXY43 - ZXY34 + ZXY33
        E = A + A - B - C
        A3SQ = A3*A3
        B3SQ = B3*B3
        P02 = (2.0*(Z3B3-ZY33)+Z3B3-ZY34)*B3
        P03 = (-2.0*Z3B3+ZY34+ZY33)*B3SQ
        P12 = (2.0*(ZX3B3-ZXY33)+ZX3B3-ZXY34)*B3
        P13 = (-2.0*ZX3B3+ZXY34+ZXY33)*B3SQ
        P20 = (2.0*(Z3A3-ZX33)+Z3A3-ZX43)*A3
        P21 = (2.0*(ZY3A3-ZXY33)+ZY3A3-ZXY43)*A3
        P22 = (3.0*(A+E)+D)*A3*B3
        P23 = (-3.0*E-B-D)*A3*B3SQ
        P30 = (-2.0*Z3A3+ZX43+ZX33)*A3SQ
        P31 = (-2.0*ZY3A3+ZXY43+ZXY33)*A3SQ
        P32 = (-3.0*E-C-D)*B3*A3SQ
        P33 = (D+E+E)*A3SQ*B3SQ
C COMPUTATION OF THE POLYNOMIAL
  690   DY = VK - Y3
        Q0 = P00 + DY*(P01+DY*(P02+DY*P03))
        Q1 = P10 + DY*(P11+DY*(P12+DY*P13))
        Q2 = P20 + DY*(P21+DY*(P22+DY*P23))
        Q3 = P30 + DY*(P31+DY*(P32+DY*P33))
        DX = UK - X3
        W(K) = Q0 + DX*(Q1+DX*(Q2+DX*Q3))
  700 CONTINUE
C NORMAL EXIT
      RETURN
C ERROR EXIT
  710 WRITE (IU0,99999)
      GO TO 800
  720 WRITE (IU0,99998)
      GO TO 800
  730 WRITE (IU0,99997)
      GO TO 800
  740 WRITE (IU0,99996)
      GO TO 760
  750 WRITE (IU0,99995)
  760 WRITE (IU0,99994) IX, X(IX)
      GO TO 800
  770 WRITE (IU0,99993)
      GO TO 790
  780 WRITE (IU0,99992)
  790 WRITE (IU0,99991) IY, Y(IY)
  800 WRITE (IU0,99990) LX0, LY0, N0
      RETURN
C FORMAT STATEMENTS
99999 FORMAT(1X/23H  ***   LX = 1 OR LESS./)
99998 FORMAT(1X/23H  ***   LY = 1 OR LESS./)
99997 FORMAT(1X/22H  ***   N = 0 OR LESS./)
99996 FORMAT(1X/27H  ***   IDENTICAL X VALUES./)
99995 FORMAT(1X/33H  ***   X VALUES OUT OF SEQUENCE./)
99994 FORMAT(7H   IX =, I6, 10X, 7HX(IX) =, E12.3)
99993 FORMAT(1X/27H  ***   IDENTICAL Y VALUES./)
99992 FORMAT(1X/33H  ***   Y VALUES OUT OF SEQUENCE./)
99991 FORMAT(7H   IY =, I6, 10X, 7HY(IY) =, E12.3)
99990 FORMAT(7H   LX =, I6, 10X, 4HLY =, I6, 10X, 3HN =, I7/
     * 36H ERROR DETECTED IN ROUTINE    ITPLBV)
      END
      SUBROUTINE SFCFIT(IU, LX, LY, X, Y, Z, MX, MY, NU, NV, U,            B  10
     * V, W)
C SMOOTH SURFACE FITTING
C THIS SUBROUTINE FITS A SMOOTH SURFACE OF A SINGLE-VALUED
C BIVARIATE FUNCTION Z = Z(X,Y) TO A SET OF INPUT DATA
C POINTS GIVEN AT INPUT GRID POINTS IN AN X-Y PLANE.  IT
C GENERATES A SET OF OUTPUT GRID POINTS BY EQUALLY DIVIDING
C THE X AND Y COORDINATES IN EACH INTERVAL BETWEEN A PAIR
C OF INPUT GRID POINTS, INTERPOLATES THE Z VALUE FOR THE
C X AND Y VALUES OF EACH OUTPUT GRID POINT, AND GENERATES
C A SET OF OUTPUT POINTS CONSISTING OF INPUT DATA POINTS
C AND THE INTERPOLATED POINTS.
C THE METHOD IS BASED ON A PIECE-WISE FUNCTION COMPOSED OF
C A SET OF BICUBIC POLYNOMIALS IN X AND Y.  EACH POLYNOMIAL
C IS APPLICABLE TO A RECTANGLE OF THE INPUT GRID IN THE X-Y
C PLANE.  EACH POLYNOMIAL IS DETERMINED LOCALLY.
C THE INPUT PARAMETERS ARE
C IU  = LOGICAL UNIT NUMBER OF STANDARD OUTPUT UNIT
C LX  = NUMBER OF INPUT GRID POINTS IN THE X COORDINATE
C       (MUST BE 2 OR GREATER)
C LY  = NUMBER OF INPUT GRID POINTS IN THE Y COORDINATE
C       (MUST BE 2 OR GREATER)
C X   = ARRAY OF DIMENSION LX STORING THE X COORDINATES
C       OF INPUT GRID POINTS (IN ASCENDING OR DESCENDING
C       ORDER)
C Y   = ARRAY OF DIMENSION LY STORING THE Y COORDINATES
C       OF INPUT GRID POINTS (IN ASCENDING OR DESCENDING
C       ORDER)
C Z   = DOUBLY-DIMENSIONED ARRAY OF DIMENSION (LX,LY)
C       STORING THE VALUES OF THE FUNCTION AT INPUT
C       GRID POINTS
C MX  = NUMBER OF SUBINTERVALS BETWEEN EACH PAIR OF
C       INPUT GRID POINTS IN THE X COORDINATE
C       (MUST BE 2 OR GREATER)
C MY  = NUMBER OF SUBINTERVALS BETWEEN EACH PAIR OF
C       INPUT GRID POINTS IN THE Y COORDINATE
C       (MUST BE 2 OR GREATER)
C NU  = NUMBER OF OUTPUT GRID POINTS IN THE X COORDINATE
C     = (LX-1)*MX+1
C NV  = NUMBER OF OUTPUT GRID POINTS IN THE Y COORDINATE
C     = (LY-1)*MY+1
C THE OUTPUT PARAMETERS ARE
C U   = ARRAY OF DIMENSION NU WHERE THE X COORDINATES OF
C       OUTPUT POINTS ARE TO BE DISPLAYED
C V   = ARRAY OF DIMENSION NV WHERE THE Y COORDINATES OF
C       OUTPUT POINTS ARE TO BE DISPLAYED
C W   = DOUBLY-DIMENSIONED ARRAY OF DIMENSION (NU,NV)
C       WHERE THE Z COORDINATES OF OUTPUT POINTS ARE TO
C       BE DISPLAYED
C SOME VARIABLES INTERNALLY USED ARE
C ZA  = DIVIDED DIFFERENCE OF Z WITH RESPECT TO X
C ZB  = DIVIDED DIFFERENCE OF Z WITH RESPECT TO Y
C ZAB = SECOND ORDER DIVIDED DIFFERENCE OF Z WITH
C       RESPECT TO X AND Y
C ZX  = PARTIAL DERIVATIVE OF Z WITH RESPECT TO X
C ZY  = PARTIAL DERIVATIVE OF Z WITH RESPECT TO Y
C ZXY = SECOND ORDER PARTIAL DERIVATIVE OF Z WITH
C       RESPECT TO X AND Y
C DECLARATION STATEMENTS
      DIMENSION X(LX), Y(LY), Z(LX,LY), U(NU), V(NV), W(NU,NV)
      DIMENSION ZA(4,2), ZB(5), ZAB(2,3), ZX(2), ZY(2), ZXY(2)
      EQUIVALENCE (Z3A2,ZA(1,1)), (Z3A3,ZA(2,1)), (Z3A4,ZA(3,1)),
     * (Z3A5,ZA(4,1)), (Z4A2,ZA(1,2)), (Z4A3,ZA(2,2)), (Z4A4,ZA(3,2)),
     * (Z4A5,ZA(4,2)), (Z4B1,ZB(1)), (Z4B2,ZB(2)), (Z4B3,ZB(3)),
     * (Z4B4,ZB(4)), (Z4B5,ZB(5)), (ZA3B2,ZAB(1,1)),
     * (ZA4B2,ZAB(2,1)), (ZA3B3,ZAB(1,2)), (ZA4B3,ZAB(2,2)),
     * (ZA3B4,ZAB(1,3)), (ZA4B4,ZAB(2,3)), (ZX43,ZX(1)),
     * (ZX44,ZX(2)), (ZY43,ZY(1)), (ZY44,ZY(2)),
     * (ZXY43,ZXY(1)), (ZXY44,ZXY(2)), (P00,Z33), (P01,ZY33),
     * (P10,ZX33), (P11,ZXY33)
      EQUIVALENCE (IXM1,JX), (IXML,JY), (DU,DV,DX,DY),
     * (FMX,RMX,FMY,RMY,SW,E), (W2,WY2,A,Q0), (W3,WY3,B,Q1),
     * (WX2,C,Q2), (WX3,D,Q3), (Z3A2,P02), (Z4A2,P03),
     * (Z4B1,P12), (Z4B2,P13), (Z4B4,P20), (Z4B5,P21),
     * (ZA3B2,P22), (ZA3B4,P23)
C PRELIMINARY PROCESSING
C SETTING OF SOME INPUT PARAMETERS TO LOCAL VARIABLES
      IU0 = IU
      LX0 = LX
      LXM1 = LX0 - 1
      LXM2 = LXM1 - 1
      LY0 = LY
      LYM1 = LY0 - 1
      LYM2 = LYM1 - 1
      MX0 = MX
      MXP1 = MX0 + 1
      MXM1 = MX0 - 1
      MY0 = MY
      MYP1 = MY0 + 1
      MYM1 = MY0 - 1
      NU0 = NU
      NV0 = NV
C ERROR CHECK
      IF (LXM2.LT.0) GO TO 400
      IF (LYM2.LT.0) GO TO 410
      IF (MXM1.LE.0) GO TO 420
      IF (MYM1.LE.0) GO TO 430
      IF (NU0.NE.LXM1*MX0+1) GO TO 440
      IF (NV0.NE.LYM1*MY0+1) GO TO 450
      IX = 2
      IF (X(1)-X(2)) 10, 460, 30
   10 DO 20 IX=3,LX0
        IF (X(IX-1)-X(IX)) 20, 460, 470
   20 CONTINUE
      GO TO 50
   30 DO 40 IX=3,LX0
        IF (X(IX-1)-X(IX)) 470, 460, 40
   40 CONTINUE
   50 IY = 2
      IF (Y(1)-Y(2)) 60, 490, 80
   60 DO 70 IY=3,LY0
        IF (Y(IY-1)-Y(IY)) 70, 490, 500
   70 CONTINUE
      GO TO 100
   80 DO 90 IY=3,LY0
        IF (Y(IY-1)-Y(IY)) 500, 490, 90
   90 CONTINUE
C COMPUTATION OF THE U ARRAY
  100 FMX = MX0
      RMX = 1.0/FMX
      KU = 1
      X4 = X(1)
      U(1) = X4
      DO 120 IX=2,LX0
        X3 = X4
        X4 = X(IX)
        DU = (X4-X3)*RMX
        DO 110 JX=1,MXM1
          KU = KU + 1
          U(KU) = U(KU-1) + DU
  110   CONTINUE
        KU = KU + 1
        U(KU) = X4
  120 CONTINUE
C COMPUTATION OF THE V ARRAY
      FMY = MY0
      RMY = 1.0/FMY
      KV = 1
      Y4 = Y(1)
      V(1) = Y4
      DO 140 IY=2,LY0
        Y3 = Y4
        Y4 = Y(IY)
        DV = (Y4-Y3)*RMY
        DO 130 JY=1,MYM1
          KV = KV + 1
          V(KV) = V(KV-1) + DV
  130   CONTINUE
        KV = KV + 1
        V(KV) = Y4
  140 CONTINUE
C MAIN DO-LOOPS
      JYMX = MY0
      KV0 = 0
      DO 390 IY=2,LY0
        IYM2 = IY - 2
        IYM3 = IYM2 - 1
        IYML = IY - LY0
        IYML1 = IYML + 1
        IX6 = 0
        IF (IYML.EQ.0) JYMX = MYP1
        JXMX = MX0
        KU0 = 0
        DO 380 IX=1,LX0
          IXM1 = IX - 1
          IXML = IX - LX0
          IF (IXML.EQ.0) JXMX = MXP1
C ROUTINES TO PICK UP NECESSARY X, Y, AND Z VALUES, TO
C COMPUTE THE ZA, ZB, AND ZAB VALUES, AND TO ESTIMATE THEM
C WHEN NECESSARY
C PRELIMINARY WHEN IX.EQ.1
          IF (IXM1.NE.0) GO TO 150
          Y3 = Y(IY-1)
          Y4 = Y(IY)
          B3 = 1.0/(Y4-Y3)
          B3SQ = B3*B3
          IF (IYM2.GT.0) B2 = 1.0/(Y3-Y(IY-2))
          IF (IYM3.GT.0) B1 = 1.0/(Y(IY-2)-Y(IY-3))
          IF (IYML.LT.0) B4 = 1.0/(Y(IY+1)-Y4)
          IF (IYML1.LT.0) B5 = 1.0/(Y(IY+2)-Y(IY+1))
          GO TO 180
C TO SAVE THE OLD VALUES
  150     Z3A2 = Z3A3
          Z4A2 = Z4A3
          X3 = X4
          Z33 = Z43
          Z3B3 = Z4B3
          A3 = A4
          A3SQ = A3*A3
          Z3A3 = Z3A4
          Z4A3 = Z4A4
          ZA3B2 = ZA4B2
          ZA3B3 = ZA4B3
          ZA3B4 = ZA4B4
  160     X4 = X5
          Z43 = Z53
          Z4B1 = Z5B1
          Z4B2 = Z5B2
          Z4B3 = Z5B3
          Z4B4 = Z5B4
          Z4B5 = Z5B5
          A4 = A5
          Z3A4 = Z3A5
          Z4A4 = Z4A5
          ZA4B2 = ZA5B2
          ZA4B3 = ZA5B3
          ZA4B4 = ZA5B4
  170     X5 = X6
          Z53 = Z63
          Z54 = Z64
          Z5B1 = Z6B1
          Z5B2 = Z6B2
          Z5B3 = Z6B3
          Z5B4 = Z6B4
          Z5B5 = Z6B5
C TO COMPUTE THE ZA, ZB, AND ZAB VALUES AND
C TO ESTIMATE THE ZB VALUES
C WHEN (IY.LE.3).OR.(IY.GE.LY-1)
  180     IX6 = IX6 + 1
          IF (IX6.GT.LX0) GO TO 260
          X6 = X(IX6)
          Z63 = Z(IX6,IY-1)
          Z64 = Z(IX6,IY)
          Z6B3 = (Z64-Z63)*B3
          IF (LYM2.EQ.0) GO TO 200
          IF (IYM2.EQ.0) GO TO 190
          Z62 = Z(IX6,IY-2)
          Z6B2 = (Z63-Z62)*B2
          IF (IYML.NE.0) GO TO 190
          Z6B4 = Z6B3 + Z6B3 - Z6B2
          GO TO 210
  190     Z65 = Z(IX6,IY+1)
          Z6B4 = (Z65-Z64)*B4
          IF (IYM2.NE.0) GO TO 210
          Z6B2 = Z6B3 + Z6B3 - Z6B4
          GO TO 210
  200     Z6B2 = Z6B3
          Z6B4 = Z6B3
  210     IF (IYM3.LE.0) GO TO 220
          Z6B1 = (Z62-Z(IX6,IY-3))*B1
          GO TO 230
  220     Z6B1 = Z6B2 + Z6B2 - Z6B3
  230     IF (IYML1.GE.0) GO TO 240
          Z6B5 = (Z(IX6,IY+2)-Z65)*B5
          GO TO 250
  240     Z6B5 = Z6B4 + Z6B4 - Z6B3
  250     IF (IX6.EQ.1) GO TO 170
          A5 = 1.0/(X6-X5)
          Z3A5 = (Z63-Z53)*A5
          Z4A5 = (Z64-Z54)*A5
          ZA5B2 = (Z6B2-Z5B2)*A5
          ZA5B3 = (Z6B3-Z5B3)*A5
          ZA5B4 = (Z6B4-Z5B4)*A5
          IF (IX6.EQ.2) GO TO 160
          GO TO 280
C TO ESTIMATE THE ZA AND ZAB VALUES
C WHEN (IX.GE.LX-1).AND.(LX.GT.2)
  260     IF (LXM2.EQ.0) GO TO 270
          Z3A5 = Z3A4 + Z3A4 - Z3A3
          Z4A5 = Z4A4 + Z4A4 - Z4A3
          IF (IXML.EQ.0) GO TO 290
          ZA5B2 = ZA4B2 + ZA4B2 - ZA3B2
          ZA5B3 = ZA4B3 + ZA4B3 - ZA3B3
          ZA5B4 = ZA4B4 + ZA4B4 - ZA3B4
          GO TO 290
C TO ESTIMATE THE ZA AND ZAB VALUES
C WHEN (IX.GE.LX-1).AND.(LX.EQ.2)
  270     Z3A5 = Z3A4
          Z4A5 = Z4A4
          IF (IXML.EQ.0) GO TO 290
          ZA5B2 = ZA4B2
          ZA5B3 = ZA4B3
          ZA5B4 = ZA4B4
C TO ESTIMATE THE ZA AND ZAB VALUES
C WHEN IX.EQ.1
  280     IF (IXM1.NE.0) GO TO 290
          Z3A3 = Z3A4 + Z3A4 - Z3A5
          Z3A2 = Z3A3 + Z3A3 - Z3A4
          Z4A3 = Z4A4 + Z4A4 - Z4A5
          Z4A2 = Z4A3 + Z4A3 - Z4A4
          ZA3B2 = ZA4B2 + ZA4B2 - ZA5B2
          ZA3B3 = ZA4B3 + ZA4B3 - ZA5B3
          ZA3B4 = ZA4B4 + ZA4B4 - ZA5B4
          GO TO 300
C NUMERICAL DIFFERENTIATION   ---   TO DETERMINE PARTIAL
C DERIVATIVES ZX, ZY, AND ZXY AS WEIGHTED MEANS OF DIVIDED
C DIFFERENCES ZA, ZB, AND ZAB, RESPECTIVELY
C TO SAVE THE OLD VALUES WHEN IX.NE.1
  290     ZX33 = ZX43
          ZX34 = ZX44
          ZY33 = ZY43
          ZY34 = ZY44
          ZXY33 = ZXY43
          ZXY34 = ZXY44
C NEW COMPUTATION
  300     DO 350 JY=1,2
            W2 = ABS(ZA(4,JY)-ZA(3,JY))
            W3 = ABS(ZA(2,JY)-ZA(1,JY))
            SW = W2 + W3
            IF (SW.EQ.0.0) GO TO 310
            WX2 = W2/SW
            WX3 = W3/SW
            GO TO 320
  310       WX2 = 0.5
            WX3 = 0.5
  320       ZX(JY) = WX2*ZA(2,JY) + WX3*ZA(3,JY)
            W2 = ABS(ZB(JY+3)-ZB(JY+2))
            W3 = ABS(ZB(JY+1)-ZB(JY))
            SW = W2 + W3
            IF (SW.EQ.0.0) GO TO 330
            WY2 = W2/SW
            WY3 = W3/SW
            GO TO 340
  330       WY2 = 0.5
            WY3 = 0.5
  340       ZY(JY) = WY2*ZB(JY+1) + WY3*ZB(JY+2)
            ZXY(JY) = WY2*(WX2*ZAB(1,JY)+WX3*ZAB(2,JY)) +
     *      WY3*(WX2*ZAB(1,JY+1)+WX3*ZAB(2,JY+1))
  350     CONTINUE
          IF (IXM1.EQ.0) GO TO 380
C DETERMINATION OF THE COEFFICIENTS OF THE POLYNOMIAL
          ZX3B3 = (ZX34-ZX33)*B3
          ZX4B3 = (ZX44-ZX43)*B3
          ZY3A3 = (ZY43-ZY33)*A3
          ZY4A3 = (ZY44-ZY34)*A3
          A = ZA3B3 - ZX3B3 - ZY3A3 + ZXY33
          B = ZX4B3 - ZX3B3 - ZXY43 + ZXY33
          C = ZY4A3 - ZY3A3 - ZXY34 + ZXY33
          D = ZXY44 - ZXY43 - ZXY34 + ZXY33
          E = A + A - B - C
          P02 = (2.0*(Z3B3-ZY33)+Z3B3-ZY34)*B3
          P03 = (-2.0*Z3B3+ZY34+ZY33)*B3SQ
          P12 = (2.0*(ZX3B3-ZXY33)+ZX3B3-ZXY34)*B3
          P13 = (-2.0*ZX3B3+ZXY34+ZXY33)*B3SQ
          P20 = (2.0*(Z3A3-ZX33)+Z3A3-ZX43)*A3
          P21 = (2.0*(ZY3A3-ZXY33)+ZY3A3-ZXY43)*A3
          P22 = (3.0*(A+E)+D)*A3*B3
          P23 = (-3.0*E-B-D)*A3*B3SQ
          P30 = (-2.0*Z3A3+ZX43+ZX33)*A3SQ
          P31 = (-2.0*ZY3A3+ZXY43+ZXY33)*A3SQ
          P32 = (-3.0*E-C-D)*B3*A3SQ
          P33 = (D+E+E)*A3SQ*B3SQ
C COMPUTATION OF THE POLYNOMIAL
          DO 370 JY=1,JYMX
            KV = KV0 + JY
            DY = V(KV) - Y3
            Q0 = P00 + DY*(P01+DY*(P02+DY*P03))
            Q1 = P10 + DY*(P11+DY*(P12+DY*P13))
            Q2 = P20 + DY*(P21+DY*(P22+DY*P23))
            Q3 = P30 + DY*(P31+DY*(P32+DY*P33))
            DO 360 JX=1,JXMX
              KU = KU0 + JX
              DX = U(KU) - X3
              W(KU,KV) = Q0 + DX*(Q1+DX*(Q2+DX*Q3))
  360       CONTINUE
  370     CONTINUE
          KU0 = KU0 + MX0
  380   CONTINUE
        KV0 = KV0 + MY0
  390 CONTINUE
C NORMAL EXIT
      RETURN
C ERROR EXIT
  400 WRITE (IU0,99999)
      GO TO 520
  410 WRITE (IU0,99998)
      GO TO 520
  420 WRITE (IU0,99997)
      GO TO 520
  430 WRITE (IU0,99996)
      GO TO 520
  440 WRITE (IU0,99995)
      GO TO 520
  450 WRITE (IU0,99994)
      GO TO 520
  460 WRITE (IU0,99993)
      GO TO 480
  470 WRITE (IU0,99992)
  480 WRITE (IU0,99991) IX, X(IX)
      GO TO 520
  490 WRITE (IU0,99990)
      GO TO 510
  500 WRITE (IU0,99989)
  510 WRITE (IU0,99988) IY, Y(IY)
  520 WRITE (IU0,99987) LX0, MX0, NU0, LY0, MY0, NV0
      RETURN
C FORMAT STATEMENTS
99999 FORMAT(1X/23H  ***   LX = 1 OR LESS./)
99998 FORMAT(1X/23H  ***   LY = 1 OR LESS./)
99997 FORMAT(1X/23H  ***   MX = 1 OR LESS./)
99996 FORMAT(1X/23H  ***   MY = 1 OR LESS./)
99995 FORMAT(1X/26H  ***   IMPROPER NU VALUE./)
99994 FORMAT(1X/26H  ***   IMPROPER NV VALUE./)
99993 FORMAT(1X/27H  ***   IDENTICAL X VALUES./)
99992 FORMAT(1X/33H  ***   X VALUES OUT OF SEQUENCE./)
99991 FORMAT(7H   IX =, I6, 10X, 7HX(IX) =, E12.3)
99990 FORMAT(1X/27H  ***   IDENTICAL Y VALUES./)
99989 FORMAT(1X/33H  ***   Y VALUES OUT OF SEQUENCE./)
99988 FORMAT(7H   IY =, I6, 10X, 7HY(IY) =, E12.3)
99987 FORMAT(7H   LX =, I6, 10X, 4HMX =, I6, 10X, 4HNU =, I6/
     * 7H   LY =, I6, 10X, 4HMY =, I6, 10X, 4HNV =, I6/6H ERROR,
     * 30H DETECTED IN ROUTINE    SFCFIT)

      END

      subroutine intrpl(n,x,y,nn,u,v)
      real*8 x(n), y(n), b(n), c(n), d(n), u(nn), v(nn), seval, vold(nn)
      ii = 1

c$$$      call intrplold(n,x,y,nn,u,vold)
      call spline(n,x,y,b,c,d)
      do i = 1, nn
         v(i) = seval(ii, n, u(i), x, y, b, c, d)
c$$$         if (mod(i,100).eq.0) print'(i5,1p,5e12.4)',
c$$$     >      i,x(1),u(i),x(n),v(i),vold(i)
      enddo
      return
      
      end
      
      subroutine spline (n, x, y, b, c, d)
      integer n
      real*8 x(n), y(n), b(n), c(n), d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end

      double precision function seval(i, n, u, x, y, b, c, d)
      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n)
c
c  this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
      integer i, j, k
      double precision dx
c$$$      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end


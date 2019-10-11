c-----------------------------------------------------------------
c***************** Wigner coefficients ***************************
c-----------------------------------------------------------------
      function CGC0(WJ1,WJ2,WJ)
c     CGC0=CGC(WJ1,0.,WJ2,0.,WJ,0.)
      CGC0=((-1)**(NINT(WJ1-WJ2)))*(sqrt(2.*WJ+1.))*
     >   COF3J(WJ1,WJ2,WJ,0.,0.,0.)
      return
      end
c****************************************************************
      real function CGC(WJ1,WM1,WJ2,WM2,WJ,WM)
      CGC=((-1)**(NINT(WJ1-WJ2+WM)))*(sqrt(2.*WJ+1.))*
     >   COF3J(WJ1,WJ2,WJ,WM1,WM2,-WM)
      return
      end
c****************************************************************
C     FILE WIGNER.FOR
C
C     ================================================================
C
      FUNCTION COF12J ( X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4)
      REAL     COF12J, COF6J, X1, X2, X3, X4, Y1, Y2, Y3, Y4,
     >           Z1, Z2, Z3, Z4, A1, A2, A3, A4, B1, B2, B3, B4,
     >           G, G1, G2, D, PHASE, SUM, FUN
      INTEGER*4  NG1, NG2, NG3, I
C
C      *****************************************************************
C
C      THIS FUNCTION COMPUTES THE 12J SYMBOL IN TERMS OF 6J SYMBOL
C      SEE R. J. ORD-SMITH PHYS. REV. 94, 1227 (1954)
C
C     (X1  X2  X3  X4  )     S      -L       (X1 Z1 L ) (Z2 X2 L )
C     (  Y1  Y2  Y3  Y4) =(-) SUM(-)  (2L+1) (X4 Z4 Y4) (X1 Z1 Y1)
C     (Z1  Z2  Z3  Z4  )       L
C
C                                          * (X3 Z3 L ) (X4 Z4 L )
C                                            (Z2 X2 Y2) (Z3 X3 Y3)
C
C
C      *****************************************************************
C
      A1 = ABS( X1 - Z1 )
      A2 = ABS( X2 - Z2 )
      A3 = ABS( X3 - Z3 )
      A4 = ABS( X4 - Z4 )
      B1 = X1 + Z1
      B2 = X2 + Z2
      B3 = X3 + Z3
      B4 = X4 + Z4
      G1 = MAX(A1,A2,A3,A4) + 1.0D0
      G2 = MIN(B1,B2,B3,B4) + 1.0D0
      NG1 = G1 + 0.1
      NG2 = G2 + 0.1
      D = G1 - NG1
      PHASE = B1 + B2 + B3 + B4 + Y1 + Y2 + Y3 + Y4
      COF12J = 1.0D0
      SUM = 0.0D0
      IF ( NG2 .LT. NG1 )GO TO 2
  3   DO 1 I=NG1,NG2
      G = I + D - 1.0D0
      NG3 = PHASE + 3.0*G + 0.1
      FUN = ((-1.D0)**NG3)*(2.0*G+1.D0)*COF6J(X1,Z1,G,X4,Z4,Y4)*
     1       COF6J(Z2,X2,G,X1,Z1,Y1)*COF6J(X3,Z3,G,Z2,X2,Y2)*
     2       COF6J(X4,Z4,G,Z3,X3,Y3)
  1   SUM = SUM + FUN
  2   CONTINUE
      COF12J = COF12J * SUM
      RETURN
      END
C
C     ============================================================
C
      FUNCTION COF3J(A1,A2,A3,A4,A5,A6)
      REAL      A1, A2, A3, A4, A5, A6, A7, A8, A9,
     >            COF3J, CNJ
      real*8 sym3
      INTEGER*4   IC
C
      A7 = 0.0D0
      A8 = 0.0D0
      A9 = 0.0D0
      IC = 1
      call threej(nint(2*a1),nint(2*a2),nint(2*a3),
     >   nint(2*a4),nint(2*a5),nint(2*a6),sym3)
      cof3j = sym3
      if (a4.ne.0.0.or.a5.ne.0.0) then
         COF3Jold = CNJ(IC,A1,A2,A3,A4,A5,A6,A7,A8,A9)
         if (abs(cof3j-cof3jold).gt.1e-4) then
            print*,'Checking COF3J and THREEJ:',cof3j,cof3jold
         endif
      endif 
      RETURN
      END
C
C     =============================================================
C
      FUNCTION COF6J(A1,A2,A3,A4,A5,A6)
      REAL     A1, A2, A3, A4, A5, A6, A7, A8, A9,
     >           COF6J, CNJ
      INTEGER*4  IC
C      
      A7 = 0.0D0
      A8 = 0.0D0
      A9 = 0.0D0
      IC = 2
      COF6J = CNJ(IC,A1,A2,A3,A4,A5,A6,A7,A8,A9)
      RETURN
      END
C
C     =============================================================
C
      FUNCTION COF9J(A1,A2,A3,A4,A5,A6,A7,A8,A9)
      REAL     A1, A2, A3, A4, A5, A6, A7, A8, A9,
     >           COF9J, CNJ
      INTEGER*4  IC
C     
      IC = 3
      COF9J = CNJ(IC,A1,A2,A3,A4,A5,A6,A7,A8,A9)
      RETURN
      END
C
C     =============================================================
C
      FUNCTION CNJ(NJ,A1,A2,A3,A4,A5,A6,A7,A8,A9)
C
      REAL      CNJ, A1, A2, A3, A4, A5, A6, A7, A8, A9,
     >            AA, BB, CC, DD, EE, FF, GG, HH, PP,
     >            X, Y, Z, XX, YY, ZZ, XXX, YYY, ZZZ,
     >            FG01BD, FG01CD, FG01DD, FG02BD, FG02CD,
     >            FG03BD, CLEBSH, RACAH, ANINE, A, B, C,
     >            H, AY, RA, RB, UUU, VVV, WWW, Q, E, D,
     >            F, FACTOR
C
      INTEGER*4   NJ, K, K1, K2, K3, K4, K5, K6, KA, KB, KC,
     >            KK1, KK2, KK3, KK4, KK5, KK6, KK7, KK8, KK9,
     >            KUP, M1, M2, M3, M4, M5, M6, M7, M8, M9, M10,
     >            MM1, MM2, MM3, MM4, MM5, MA, MB, MC, MD, I,
     >            IJ, IX, IY, IZ, IYY, IAY, IJPAR, IERR, IERCT,
     >            N4, N5, N6, N5PAR, J, JJJ, JS, JSPAR, JQ, 
     >            KEY, KEYTRI, KEYRAC
C
C   -----------------------------------------------------------------
C
C  2. ARGUMENT LIST
C
C     ALL ARGUMENTS MUST BE REAL INTEGERS OR HALF-ODD-INTEGERS (REAL
C  IN THE SINGLE-LENGTH ROUTINES, DOUBLE PRECISION IN THE DOUBLE-LENGTH ONES).
C  THE NEAREST EXACT VALUE WILL BE ASSUMED, BUT ERRORS MAY OCCUR IF THE
C  VALUES GIVEN IN THE CALLING ROUTINE ARE +OR-0.05 OR MORE IN ERROR.
C
C
C  3. RESTRICTIONS
C
C     THE SUM OF THE THREE ANGULAR MOMENTA APPEARING IN ANY
C  "TRIANGULAR CONDITION" MUST NOT EXCEED 100.0.(MMAX) THIS LIMIT CAN BE RAISED
C  IF NECESSARY BY RECOMPILING WITH LARGER DIMENSIONS FOR THE ARRAYS H AND
C  J AND A CORRESPONDINGLY LARGER UPPER LIMIT ON THE INDEX OF THE FIRST DO
C  LOOP.
C
C     THE FOLLOWING "GEOMETRICAL" CONDITIONS ARE TESTED BY THE
C  PROGRAMS AND THE CORRECT VALUE OF ZERO RETURNED IN CASE OF VIOLATION:
C
C     (A) ALL TRIANGULAR CONDITIONS SATISFIED
C     (B) ALL ANGULAR MOMENTA NON-NEGATIVE
C     (C) IN COF3J AND FGO1B, (A+X), (B+Y) AND (C+Z) ARE INTEGRAL
C     (D) IN COF3J, X+Y+Z = 0.0; IN FGO1B, X+Y = Z
C     (E) IN FGO1C AND D, A,B AND C ALL INEGRAL WITH A+B+C EVEN.
C
C     SINCE A VIOLATION OF THESE CONDITIONS MAY BE THE RESULT OF AN
C  ERROR IN THE CALLING ROUTINE WHICH SETS UP THE ARGUMENTS, AN "ERROR
C  CHECK" IS PROVIDED.
C  A LABELLED COMMON AREA IS SET UP:
C
C     COMMON/FGERCM/IERR,IERCT
C
C  WHERE THE TWO MEMBERS ARE INTEGER*4. IERR WILL BE UNITY AT RETURN IF ANY
C  CONDITION WAS VIOLATED; ZERO IF NOT. THIS COUNT IS RESET AT EACH CALL
C  OF ANY OF THE ROUTINES IN THE PACKAGE. IERCT IS ZERO AT THE START OF THE
C  JOB AND IS INCREASED BY UNITY WHENEVER A CONDITION IS VIOLATED IN A CALL
C  OF ANY OF THE ROUTINES OF THE PACKAGE; THUS A SINGLE CHECK AT THE END OF
C  THE ENTIRE JOB CAN ENSURE THAT NO VIOLATIONS OCCURED DURING THE JOB.
C
C
C  4. GENERAL
C
C     THE 9 NAMES ARE DIFFERENT ENTRY POINTS TO A SINGLE ROUTINE
C  (CSECT NAME COF6J). THIS SAVES SPACE SINCE ALL USE THE SAME PRIVATE
C  TABLE OF FACTORIALS (1K BYTE). USERS WHO NEED ONLY ONE OF THESE
C  FUNCTIONS IN THEIR JOB AND WHO ARE PRESSED FOR SPACE MAY OF COURSE
C  RECOMPILE THE SOURCE ROUTINE AFTER REMOVING MOST OF THE CODE BELONGING
C  PRIVATELY TO THE OTHER ENTRIES (NOTE THAT C9J NEEDS MOST OF THE FGO2
C  CODE; AND THAT ALL ENTRIES NEED THE DO LOOP NEAR THE BEGINNING THAT SETS
C  UP THE FACTORIAL TABLES.
C
C     THE ROUTINES ARE SELF CONTAINED AND CAUSE NO OUTPUT; THEY SET
C  UP ONE NAMED COMMON AREA (DESCRIBED IN 3 ABOVE).
C
C     FGO3B IS COMPLETELY IDENTICAL TO FGO3A. IT IS INCLUDED FOR
C  COMPATIBILITY WITH EARLIER VERSIONS.
C
C
C
C##       FG01BD         08/01/74
C NAME FG01BD(R)                 CHECK



C
C     WIGNER AND RACAH COEFFICIENTS
C
C     FG01A - WIGNER 3-J SYMBOL
C     FG01B - CLEBSCH-GORDAN COEFFICIENT
C     FG01C & FG01D - SAME WITH ZERO MAGNETIC QUANTUM NUMBERS
C
C     FG02A - WIGNER 6-J SYMBOL
C     FG02B - RACAH COEFFICIENT
C     FG02C - U-FUNCTION (JAHN)
C
C     FG03A & FG03B - WIGNER 9-J SYMBOL
C
C
      parameter (mmax=501)
      COMMON / CNJSAVE / H(mmax), J(mmax)
      DIMENSION AY(4),IAY(4)
      COMMON/ FGERCM /IERR,IERCT
      DATA JJJ/0/
      INTPTF(Q) = nint(Q + Q)
      IPARF(I) = 4*(I/4) - I + 1
C      GO TO 77

C  COMPUTED GOTO 'SIMULATES' MULTIPLE ENTRY POINTS.

      GO TO (9993,9996,9999),NJ

C

9993  CONTINUE
      A = A1
      B = A2
      C = A3
      XX = A4
      YY = A5
      ZZ = A6

      ZZ=-ZZ
      KEY=1
      IERR=0
      GO TO 1
C
C
C
9996  CONTINUE
      UUU=A1
      VVV=A2
      WWW=A3
      XXX=A4
      YYY=A5
      ZZZ=A6


      KEY=11
      IERR=0
      GO TO 100
   77 KEY=2
      IERR=0
    1 K1=INTPTF(A)
      K2=INTPTF(B)
      K3=INTPTF(C)
      IF(KEY.GE.3) GO TO 100
      K4=INTPTF(XX)
      K5=INTPTF(YY)
      K6=INTPTF(ZZ)
C
  100 IF(JJJ.NE.0) GO TO 500
      JJJ=1
      IERCT=0
      H(1)=1.0D0
      J(1)=0
      X=0.D0
      DO 400 I=2,mmax
      X=X+1.0D0
      H(I)=H(I-1)*X
      J(I)=J(I-1)
  200 IF(H(I).LT.10.0D0) GO TO 400
      H(I)=0.01D0*H(I)
      J(I)=J(I)+2
      GO TO 200
  400 CONTINUE
C
  500 IF(KEY.LT.-5) GO TO 750
      IF(KEY.GE.3) GO TO 320
      IF((K4+K5-K6).NE.0) GO TO 710
      M1=K1+K2-K3
      M2=K2+K3-K1
      M3=K3+K1-K2
      M4=K1+K4
      M5=K1-K4
      M6=K2+K5
      M7=K2-K5
      M8=K3+K6
      M9=K3-K6
      M10=K1+K2+K3+2
C
      IF(M1.LT.0) GO TO 710
      IF(M2.LT.0) GO TO 710
      IF(M3.LT.0) GO TO 710
      IF(M4.LT.0) GO TO 710
      IF(M5.LT.0) GO TO 710
      IF(M6.LT.0) GO TO 710
      IF(M7.LT.0) GO TO 710
      IF(M8.LT.0) GO TO 710
      IF(M9.LT.0) GO TO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
      IF((M10-(M10/2)-(M10/2)).NE.0) GO TO 710
C
      Y=K3+1
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      M4=M4/2+1
      M5=M5/2+1
      M6=M6/2+1
      M7=M7/2+1
      M8=M8/2+1
      M9=M9/2+1
      M10=M10/2+1
      if (max(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10).gt.mmax) then
         print*,'Recompile wigner.f with mmax at least:',
     >      max(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
         stop 'Increse MMAX in CNJ (file wigner.f)'
      endif 
C
      Y  = SQRT(Y*H(M1)*H(M2)*H(M3)*H(M4)*H(M5)*
     >     H(M6)*H(M7)*H(M8)*H(M9)/H(M10))
      IY = (J(M1)+J(M2)+J(M3)+J(M4)+J(M5)+
     >     J(M6)+J(M7)+J(M8)+J(M9)-J(M10))/2
C
      N4=M1
      IF(N4.GT.M5)N4=M5
      IF(N4.GT.M6)N4=M6
      N4=N4-1
      M2=K2-K3-K4
      M3=K1+K5-K3
      N5=0
      IF(N5.LT.M2) N5=M2
      IF(N5.LT.M3) N5=M3
      N5PAR=IPARF(N5)
      N5=N5/2
      Z=0.0D0
      GO TO 610
C
  700 MM1=M1-N5
      MM2=M5-N5
      MM3=M6-N5
      MM4=N5-(M2/2)+1
      MM5=N5-(M3/2)+1
C
      X  = 1.D0/(H(MM1)*H(MM2)*H(MM3)*H(MM4)*H(MM5)*H(N5+1))
      IX = -J(MM1)-J(MM2)-J(MM3)-J(MM4)-J(MM5)-J(N5+1)
C
  800 IF(IX+IY)900,210,110
  900 X=0.1D0*X
      IX=IX+1
      GO TO 800
  110 X=10.0D0*X
      IX=IX-1
      GO TO 800
C
  210 IF(N5PAR.LT.0) X=-X
      Z=Z+X
  510 N5PAR=-N5PAR
      N5=N5+1
C
  610 IF(N5-N4)700,700,810
C
 710  CLEBSH=0.0D0
      IERR=1
      IERCT=IERCT+1
      GO TO 910
C
 810  CLEBSH=Z*Y
  910 GO TO(120,220),KEY
C
  220 FG01BD=CLEBSH
      RETURN
C
  120 JS=K1-K2+K6
      IF(JS.LT.0) JS=-JS
      JSPAR=IPARF(JS)
      CNJ=JSPAR*CLEBSH/SQRT(K3+1.0D0)
      ZZ=-ZZ
      RETURN
C
C
C
  320 IF(KEY.GE.10) GO TO 130
      KEY=KEY-2
      IF((K1-(K1/2)-(K1/2)).NE.0) GO TO 420
      IF((K2-(K2/2)-(K2/2)).NE.0) GO TO 420
      IF((K3-(K3/2)-(K3/2)).NE.0) GO TO 420
      IJ=K1+K2+K3
      IJPAR=IPARF(IJ)
      IF(IJPAR.LE.0) GO TO 420
      M1=IJ-K1-K1
      M2=IJ-K2-K2
      M3=IJ-K3-K3
      M4=IJ+2
      IF(M1.LT.0) GO TO 420
      IF(M2.LT.0) GO TO 420
      IF(M3.LT.0) GO TO 420
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      M4=IJ/2+2
      Y  = SQRT(H(M1)*H(M2)*H(M3)/H(M4))
      IY = (J(M1)+J(M2)+J(M3)-J(M4))/2
      IJ=IJ/2
      IJPAR=IPARF(IJ)
      IJ=IJ/2+1
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      Z=H(IJ)/(H(M1)*H(M2)*H(M3))
      IZ=J(IJ)-J(M1)-J(M2)-J(M3)
      IZ=IZ+IY
      CLEBSH=IJPAR*Y*Z*10.0D0**IZ
      GO TO(620,720),KEY
C
  620 FG01CD=CLEBSH
      RETURN
C
  720 JQ=K2-K1
      IF(JQ.LT.0) JQ=-JQ
      IJPAR=IPARF(JQ)
      FG01DD=CLEBSH*IJPAR*SQRT(K3+1.0D0)
      RETURN
C
  420 CLEBSH=0.0D0
      IERR=1
      IERCT=IERCT+1
      GO TO(620,720),KEY
C
  130 IF(KEY.EQ.11) GO TO 450
      IF(KEY.GT.19) GO TO 750
  550 stop'should not be here as D, E, and F are not defined (wigner.f)'
      K1=INTPTF(A)
      K2=INTPTF(B)
      K3=INTPTF(E)
      K4=INTPTF(D)
      K5=INTPTF(C)
      K6=INTPTF(F)
C
C
  750 KA=K1
      KB=K2
      KC=K3
      KEYTRI=1
      GO TO 630
C
  230 KA=K4
      KB=K5
      KEYTRI=2
      GO TO 630
C
  330 KB=K2
      KC=K6
      KEYTRI=3
      GO TO 630
C
  430 KA=K1
      KB=K5
      KEYTRI=4
      GO TO 630
C
  530 YY=AY(1)*AY(2)*AY(3)*AY(4)
      IYY=IAY(1)+IAY(2)+IAY(3)+IAY(4)
      M1=(K1+K2+K4+K5)/2+2
      M2=(K1+K2-K3)/2+1
      M3=(K4+K5-K3)/2+1
      M4=(K1+K5-K6)/2+1
      M5=(K2+K4-K6)/2+1
      M6=K1+K4-K3-K6
      M7=K2+K5-K3-K6
C
      N4=M1
      IF(N4.GT.M2) N4=M2
      IF(N4.GT.M3) N4=M3
      IF(N4.GT.M4) N4=M4
      IF(N4.GT.M5) N4=M5
      N4=N4-1
      N5=0
      IF(N5.LT.M6) N5=M6
      IF(N5.LT.M7) N5=M7
      N5PAR=IPARF(N5)
      N5=N5/2
      M6=M6/2-1
      M7=M7/2-1
      Z=0.0D0
      GO TO 730
C
  140 X  = H(M1-N5)/(H(N5+1)*H(M2-N5)*H(M3-N5)*H(M4-N5)
     >     *H(M5-N5)*H(N5-M6)*H(N5-M7))
      IX = J(M1-N5)-J(N5+1)-J(M2-N5)-J(M3-N5)-J(M4-N5)
     >    -J(M5-N5)-J(N5-M6)-J(N5-M7)
  240 IF(IX+IYY)340,440,540
  340 X=0.1D0*X
      IX=IX+1
      GO TO 240
  540 X=10.0D0*X
      IX=IX-1
      GO TO 240
  440 IF(N5PAR.LT.0) X=-X
      Z=Z+X
      N5PAR=-N5PAR
      N5=N5+1
C
  730 IF(N5.LE.N4) GO TO 140
C
      RACAH=Z*YY
  840 IF(KEY.LT.-5) GO TO 160
      KEY=KEY-10
      GO TO(150,250,350),KEY
C
  830 RACAH=0.0D0
      IERR=1
      IERCT=IERCT+1
      GO TO 840
C
  150 IJPAR=IPARF(K1+K2+K4+K5)
      IF(IJPAR.LT.0) RACAH=-RACAH
      CNJ=RACAH
      RETURN
C
  250 FG02BD=RACAH
      RETURN
C
  350 FACTOR = SQRT((K3+1.0D0)*(K6+1))
      FG02CD = FACTOR*RACAH
      RETURN
  450 K1 = INTPTF(UUU)
      K2 = INTPTF(VVV)
      K3=INTPTF(WWW)
      K4=INTPTF(XXX)
      K5=INTPTF(YYY)
      K6=INTPTF(ZZZ)
      GO TO 750
C
C
C     TRIANGLE FUNCTION
C
C
  630 MA=KA+KB-KC
      MB=KA-KB+KC
      MC=-KA+KB+KC
      MD=KA+KB+KC+2
      IF(MA.LT.0) GO TO 830
      IF(MB.LT.0) GO TO 830
      IF(MC.LT.0) GO TO 830
      IF((MD-(MD/2)-(MD/2)).NE.0) GO TO 830
      MA=MA/2+1
      MB=MB/2+1
      MC=MC/2+1
      MD=MD/2+1
      if (max(ma,mb,mc,md).gt.mmax) stop
     >   'Increse MMAX in CNJ (file wigner.f)'
      AY(KEYTRI) = SQRT(H(MA)*H(MB)*H(MC)/H(MD))
      IAY(KEYTRI) = (J(MA)+J(MB)+J(MC)-J(MD))/2
      GO TO(230,330,430,530),KEYTRI
C
C
C
C
C

9999  CONTINUE
      AA=A1
      BB=A2
      CC=A3
      DD=A4
      EE=A5
      FF=A6
      GG=A7
      HH=A8
      PP=A9


C      ENTRY FG03BD(AA,BB,CC,DD,EE,FF,GG,HH,PP)
C
      KEY=-10
      IERR=0
C
      KK1=INTPTF(AA)
      KK2=INTPTF(BB)
      KK3=INTPTF(CC)
      KK4=INTPTF(DD)
      KK5=INTPTF(EE)
      KK6=INTPTF(FF)
      KK7=INTPTF(GG)
      KK8=INTPTF(HH)
      KK9=INTPTF(PP)
C
      KUP=KK1+KK9
      M1=KK4+KK8
      M2=KK2+KK6
      IF(KUP.GT.M1) KUP=M1
      IF(KUP.GT.M2) KUP=M2
C
      K=KK1-KK9
      IF(K.LT.0) K=-K
      M1=KK4-KK8
      IF(M1.LT.0) M1=-M1
      M2=KK2-KK6
      IF(M2.LT.0) M2=-M2
      IF(K.LT.M1) K=M1
      IF(K.LT.M2) K=M2
C
      ANINE=0.0D0
C
  660 IF(K.GT.KUP) GO TO 260
      K1=KK1
      K2=KK4
      K3=KK7
      K4=KK8
      K5=KK9
      K6=K
      KEYRAC=1
      GO TO 100
C
  160 GO TO(360,460,560),KEYRAC
C
  360 RA=RACAH
      K1=KK2
      K2=KK8
      K3=KK5
      K4=KK4
      K5=KK6
      KEYRAC=2
      GO TO 750
C
  460 RB=RACAH
      K1=KK9
      K2=KK6
      K3=KK3
      K4=KK2
      K5=KK1
      KEYRAC=3
      GO TO 750
C
  560 ANINE=ANINE+RA*RB*RACAH*(K+1)
      K=K+2
      GO TO 660
C
  260 CNJ=ANINE
      FG03BD=ANINE
      RETURN
      END
C
C     ===============================================================
C
      FUNCTION BICO (A,B)
C
      REAL     BICO, A, B, FACTOR, X, AI, BJ
      INTEGER*4  IFIRST, I, J
C
C     THIS PROGRAM CAN ONLY HANDLE INTEGERS.
C
      COMMON / LOGFAC / FACTOR(230)
      DATA IFIRST / 0 /
      IF ( IFIRST .EQ. 0 ) CALL FACLOG
      IFIRST = 1
      X=A-B
      IF (X) 201,202,203
  201 BICO=0.D0
      GO TO 207
  202 BICO=1.D0
      GO TO 207
  203 IF (B) 201,202,204
 204  CONTINUE
      I=A+0.1
      J=B+0.1
      AI=I
      BJ=J
      IF(ABS(AI-A).GT.0.001) GOTO500
      IF(ABS(BJ-B).GT.0.001) GOTO500
      X=FACTOR(I)-FACTOR(J)-FACTOR(I-J)
      BICO= EXP(X)
      GO TO 207
  500 WRITE(9,501) A,B
  501 FORMAT(10X,'  A=',F10.3,'B=',F10.3,' NON INTEGER VALUES FOR BICO')
       BICO=0.0D0
 207   CONTINUE
       RETURN
      END
C
C     =============================================================
C
      SUBROUTINE FACLOG
C
      INTEGER*4         I 
      REAL            FACTOR, A
      COMMON / LOGFAC / FACTOR(230)
C
      FACTOR(1)=0.0D0
      DO 10 I=2,230
        A=I
        FACTOR(I)=FACTOR(I-1)+LOG(A)
10    CONTINUE
C
      RETURN
      END

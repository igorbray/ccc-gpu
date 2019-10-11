      SUBROUTINE COEF(GAM,DIR,S,T,KSU,N,LL,QQ,
     *P,LO1, S1,CH,KK,SIGM1,SIGM2,G,MU)
C
C  Calculation of direct and exchange Hartree-Fock coefficients 
C
c!!!
      include 'paratom.f'
      REAL*8 GAM(NSIZEqf),DIR(max10)
      INTEGER S,A,B,CH,G,T,N(max10),P(max10),LL(max10),
     *   QQ(max10),L1(max10),S1(max10),KK(NSIZEqf),SIGM1(NSIZEqf),
     *   SIGM2(NSIZEqf),PP(200),Q
      DO 10 I = 1,CH
         SIGM1(I) = 1
         SIGM2(I) = 1
         GAM(I) = 0.0
 10   continue
      mtmp = T
      I2 = 0
      NU = I2
      L1(1) = LO1
      NU1 = CH
      DO 20 I = 1,S
         PP(I) = 0
         IF((S1(1).LE.99.AND.(4*LL(I)+2).
     *      LT.QQ(I)).OR.(S1(1).GT.99.
     *      AND.(2*LL(I)+1).LT.QQ(I))) GO TO 1000
 20   CONTINUE
      IF(KSU.NE.1) GOTO 22
      DO 23 I = 1,S
         DO 230 J = I,S
            IF(I.EQ.J) GOTO 230
            II1 = IABS(LL(I)-LL(J))
            II2 = LL(I)+LL(J)
            DO 231 I1 = II1,II2,2
               NU = NU+1
               GAM(NU) = dble( G1(I1,LL(J),LL(I),QQ(I),QQ(J)) )
               IF(S1(1).LE.99) GOTO 25
               GAM(NU) = GAM(NU)*2.
               IF(DIR(I).NE.DIR(J)) GAM(NU) = 0.0
 25            KK(NU) = I1
               SIGM1(NU) = I
               SIGM2(NU) = J
 231        CONTINUE
 230     CONTINUE
 23   CONTINUE
 22   IF(KSU.NE.2.AND.KSU.NE.3) GOTO 26
      I3 = S
      IF(KSU.EQ.2) I3 = 1
      IF(KSU.EQ.3.AND.QQ(S).NE.1) write(66, 1001)
      DO 28 I = 1,S
         A = 2 *LL(I)+1
         IF(S1(1).LE.99) A = A*2
         DO 280 J = I3,S
            IF((J.LE.I.OR.KSU.NE.2).AND.
     *         (I.EQ.J.OR.KSU.NE.3)) GOTO 280
            B = 2 *LL(J)+1
            IF(S1(1).LE.99) B = B*2
            IF(A.NE.QQ(I).AND.B.NE.QQ(J)) GOTO 305
            IF(MU.NE.1) GOTO 280
            II3 = IABS(LL(I)- LL(J))
            II4 = LL(I)+ LL(J)
            DO 30 I1 = II3,II4,2
               NU = NU+1
               GAM(NU) = dble(G1(I1,LL(I),LL(J),QQ(I),QQ(J)) )
               IF(S1(1).LE.99) GOTO 31
               GAM(NU) = GAM(NU)*2.
               IF(DIR(I).NE.DIR(J)) GAM(NU) = 0.
 31            KK(NU) = I1
               SIGM1(NU) = I
               SIGM2(NU) = J
 30         continue
            GOTO 280
 305        IF(MU.NE.1) GOTO 32
            IF(S1(1).GT.99) GOTO 33
            IF((S1(1).EQ.1.OR.S1(1).EQ.3).AND.
     *         ((QQ(I).EQ.1.AND.QQ(J).EQ.B-1).
     *         OR.(QQ(J).EQ.1.AND.QQ(I).EQ.A-1))) GOTO 34
            GOTO 33
 34         II5 = IABS(LL(I)-LL(J))
            II6 = LL(I)+LL(J)
            DO 35 I1 = II5,II6,2
               NU = NU+1
               GAM(NU) = dble(G3(I1,LL(I),LL(J),L1(1),S1))
               KK(NU) = I1
               SIGM1(NU) = I
               SIGM2(NU) = J
 35         continue
            GOTO 32
 33         IF((QQ(I).EQ.1.AND.QQ(J).EQ.B-1).OR.(QQ(J)
     *         .EQ.1.AND.QQ(I).EQ.A-1)) GOTO 36
            GOTO 32
 36         II7 = IABS(LL(I)-LL(J))
            II8 = LL(I)+LL(J)
            DO 37 I1 = II7,II8,2
               NU = NU+1
               GAM(NU) = dble(G5(I1,LL(I),LL(J),L1(1)))
               IF(DIR(I).NE.DIR(J)) GAM(NU) = 0.
               KK(NU) = I1
               SIGM1(NU) = I
               SIGM2(NU) = J
 37         continue
 32         JS = 2 *LL(I)
            I4 = 2 *LL(J)
            I4 = MIN0(JS,I4)
            IF(I4.LT.2) GOTO 280
            DO 38 I1 = 2,I4,2
               I2 = I2+1
               GAM(NU1) = dble(F3(I1,LL(I),LL(J),L1(1)) )
               KK(NU1) = I1
               SIGM1(NU1) = I
               SIGM2(NU1) = J
               NU1 = NU1-1
 38         continue
 280     CONTINUE
 28   CONTINUE
 26   CONTINUE
      IF(KSU.NE.1.AND.KSU.NE.2) GOTO 391
      DO 39 I = 1,S
         IF(QQ(I).EQ.1) GOTO 39
         II9 = LL(I)+LL(I)
         IF(II9.LT.2) GOTO 39
         DO 392 I1 = 2,II9,2
            I2 = I2+1
            IF(S1(1)-99) 201,201,202
 201        GAM(NU1) = F1(I1,LL(I),QQ(I))
            GOTO 203
 202        GAM(NU1) = F4(I1,LL(I),QQ(I))
 203        KK(NU1) = I1
            SIGM1(NU1) = I
            SIGM2(NU1) = I
            NU1 = NU1-1
 392     continue
 39   CONTINUE
 391  I1 = 0
      N16 = MAX0(NU,CH)
      DO 293 N17 = 1,N16
         IF(QQ(SIGM1(N17)).EQ.0 .OR. QQ(SIGM2(N17)).EQ.0) GAM(N17) = 0.
 293  CONTINUE
      DO 40 I = 1,NU
         IF(GAM(I).NE.0.) GOTO 41
         J = I+1
 42      IF(J.GT.NU) GOTO 45
         IF(GAM(J).EQ.0.) GO TO 43
         GAM(I) = GAM(J)
         GAM(J) = 0.
         KK(I) = KK(J)
         SIGM1(I) = SIGM1(J)
         SIGM2(I) = SIGM2(J)
         GOTO 44
 43      J = J+1
      GOTO 42
 44   CONTINUE
 45   CONTINUE
41    IF(GAM(I).NE.0.) I1 = I1+1
40    CONTINUE
      G = 0
      DO 46 I = 1,CH
      IF(GAM(I).NE.0.) GOTO 47
      J = I+1
48    IF(J.GT.CH) GOTO 49
      IF(GAM(J).EQ.0.) GOTO 50
      GAM(I) = GAM(J)
      GAM(J) = 0.
      KK(I) = KK(J)
      SIGM1(I) = SIGM1(J)
      SIGM2(I) = SIGM2(J)
      GOTO 49
50    J = J+1
      GOTO 48
49    CONTINUE
47    IF(GAM(I).EQ.0.) GOTO 46
      G = G+1
46    CONTINUE
      CH = G
      G = NU
      G = I1
      IF(CH.NE.0) GOTO 53
      G = 1
      CH = 1
      GAM(1) = 0.
      KK(1) = 0
      SIGM1(1) = 1
      SIGM2(1) = 2
53    CONTINUE
      RETURN
1000  write(66, 1001)
1001  FORMAT(10X,'error in data')
      RETURN
      END

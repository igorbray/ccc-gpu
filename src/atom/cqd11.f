      SUBROUTINE CQD11(DELX,Q,OM,ALML,ALMM,
     * C,XO,DV1,DV2,   M,M1,IC1,JJ)
      REAL*8 K2,K3,KAP,OM,XO,DELX,ALML,XLAM,
     * DV1(20),DV2(20),C(20),Q(20),ALMM(20),
     * X1,X2,X3,A,U,V,B,T,X4,X5
      INTEGER P,W,F
      K3= DELX/3.0
      F=0
      DO 10 I=1,M
      Q(I)=0.0
   10 CONTINUE
      IF((OM+ALML).GT.0.0) GOTO 11
      XLAM=0.0
      GOTO 12
   11 XLAM=SQRT(ALMM(JJ))
   12 CONTINUE
C     IF(XLAM.EQ.0.0) GOTO 51
      IRA=M1+1
      IF(IRA.GT.M) GOTO 50
      DO 1 I=IRA,M,2
      X1=SQRT(ALMM(I))
      X2=X1+DELX
      X3=X2+DELX
      W=I
      IF(ABS(XLAM-X2).LT.1.D-04.OR.
     * ABS(XLAM-X3).LT.1.D-04)GOTO 50
1     CONTINUE
50    CONTINUE
      IF(ABS(XLAM-X2).GE.1.D-04) GOTO 200
      Q(W)=-2.*K3
      Q(W+1)=4.*K3
      Q(W+2)=-2.*K3
200   CONTINUE
      IF(ABS(XLAM-X3).GE.1.D-04) GOTO 201
      Q(W)=K3/3.
      Q(W+1)=-4.*K3/3.
      Q(W+2)=2.*K3
      Q(W+3)=-4.*K3/3.
      Q(W+4)=K3/3.
201   CONTINUE
51    CONTINUE
      IF(M1.LT.1) GOTO 1520
      DO 15 I=1,M1
 15   C(I)=1.0
 1520 CONTINUE
      C(M1+1)=K3
      IRA=M1+2
      IRA1=M-1
      IF(IRA.GT.IRA1) GOTO 40
      DO 4 I=IRA,IRA1,2
      C(I)=4.*K3
      C(I+1)=2.*K3
    4 CONTINUE
40    CONTINUE
      C(M)=K3
      P=1
      KAP=XO
      DO 5 W=1,M
      IF( W.EQ.JJ) GOTO 16
      IF(W.LE.M1.AND.W.NE.JJ) GOTO 18
      DV1(W)=KAP*2.0*DFLOAT(P)/
     * (OM-DFLOAT(P)*(KAP*KAP-ALML))
      GOTO 19
   18 DV1(W)=1.0*DFLOAT( P)/(OM-
     *  DFLOAT(    P)*(ALMM(W)-ALML))
      GOTO 19
   16 DV1(W)=0.0
   19 CONTINUE
      IF(W.GT.M1) KAP=KAP+DELX
    5 CONTINUE
      P=-1
      KAP=XO
      DO 6 W=1,M
      IF(W.LE.M1) GOTO 20
      DV2(W)=2.*DFLOAT(   P)*KAP/
     * (OM-DFLOAT(  P)*(KAP*KAP-ALML))
      GOTO 21
   20 DV2(W)=1.*DFLOAT(  P)/(OM-
     * DFLOAT(    P)*(ALMM(W)-ALML))
   21 CONTINUE
      IF(W.GT.M1) KAP=KAP+DELX
    6 CONTINUE
      RETURN
      END

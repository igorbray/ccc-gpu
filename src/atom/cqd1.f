      SUBROUTINE CQD1(DELX,Q,OM,ALML,ALMM,                              F4500190
     * C,XO,DV1,DV2,RE,IM,M,M1)                                         F4500200
      REAL*8 K2,K3,KAP,OM,XO,DELX,ALML,XLAM,                            F4500210
     * DV1( M),DV2( M),C(M ),Q(M ),ALMM( M),                            F4500220
     * X1,X2,X3,A,U,V,B,T,X4,X5,RE,IM,MIN                               F4500230
      INTEGER P,W,F                                                     F4500250
      K3= DELX/3.0                                                      F4500260
      F=0                                                               F4500270
      DO 10 I=1,M                                                       F4500280
      Q(I)=0.0                                                          F4500290
   10 CONTINUE                                                          F4500300
      IF((OM+ALML).GT.0.0) GOTO 11                                      F4500310
      XLAM=0.0                                                          F4500320
      GOTO 12                                                           F4500330
   11 XLAM= SQRT(OM+ALML)                                               F4500340
   12 CONTINUE                                                          F4500350
      IF(XLAM.EQ.0.0) GOTO 51                                           F4500360
      IRA=M1+1                                                          F4500370
      IF(IRA.GT.M) GOTO 50                                              F4500380
      DO 1 I=IRA,M,2                                                    F4500390
      X1= SQRT(ALMM(I))                                                 F4500400
      X2=X1+DELX                                                        F4500410
      X3=X2+DELX                                                        F4500420
      W=I                                                               F4500430
      IF(W.EQ.M) GOTO 51                                                F4500440
      IF(XLAM.GE.X1-DELX/2..AND.                                        F4500450
     * XLAM.LE.X3-DELX/2.) GOTO 50                                      F4500460
    1 CONTINUE                                                          F4500470
   50 IF(W.NE.(M1+1)) GOTO 13                                           F4500480
      A=(XLAM-X1)*(XLAM-X2)*(XLAM-X3)                                   F4500490
      U= ABS((X3-XLAM)/(XLAM-X1))                                       F4500500
      U= LOG(U)                                                         F4500510
      V=(1./(X1-XLAM)+4./(X2-XLAM)+                                     F4500520
     * 1./(X3-XLAM))*K3                                                 F4500530
      B=U-V                                                             F4500540
      T=B*A/DELX**2                                                     F4500550
      Q(W)=T/2.                                                         F4500560
      Q(W+1)=-T                                                         F4500570
      Q(W+2)=T/2.                                                       F4500580
      GOTO 14                                                           F4500590
   13 X4=X1-DELX                                                        F4500600
      X5=X4-DELX                                                        F4500610
      A=(XLAM-X1)*(XLAM-X2)*                                            F4500620
     *  (XLAM-X3)*(XLAM-X4)*(XLAM-X5)                                   F4500630
      U= ABS((X3-XLAM)/(XLAM-X5))                                       F4500640
      U= LOG(U)                                                         F4500650
      V=(1./(X5-XLAM)+4./(X4-XLAM)+2./                                  F4500660
     * (X1-XLAM)+4./(X2-XLAM)+1./(X3-XLAM))*K3                          F4500670
      B=U-V                                                             F4500680
      T=B*A/DELX**4                                                     F4500690
      Q(W-2)=T/24.                                                      F4500700
      Q(W-1)=-T/6.                                                      F4500710
      Q(W)=T/4.                                                         F4500720
      Q(W+1)=-T/6.                                                      F4500730
      Q(W+2)=T/24.                                                      F4500740
   14 CONTINUE                                                          F4500750
   51 CONTINUE                                                          F4500760
      IF(M1.EQ.0) GOTO 155                                              F4500770
      DO 15 I=1,M1                                                      F4500780
   15 C(I)=1.0                                                          F4500790
  155 CONTINUE                                                          F4500800
      C(M1+1)=K3                                                        F4500810
      IRA=M1+2                                                          F4500820
      IRA1=M-1                                                          F4500830
      IF(IRA.GT.IRA1) GOTO 40                                           F4500840
      DO 4 I=IRA,IRA1,2                                                 F4500850
      C(I)=4.*K3                                                        F4500860
      C(I+1)=2.*K3                                                      F4500870
    4 CONTINUE                                                          F4500880
   40 CONTINUE                                                          F4500890
      C(M)=K3                                                           F4500900
      RE=0.0                                                            F4500910
      KAP=XO                                                            F4500920
      DO 100 W=1,M                                                      F4500930
      DV1(W)=2.*KAP/(OM-KAP*KAP)                                        F4500940
      IF(W.LE.M1) DV1(W)=1./(OM-                                        F4500950
     * ALMM(W))                                                         F4500960
      IF(W.GT.M1) KAP=KAP+DELX                                          F4500970
100   CONTINUE                                                          F4500980
      DO 101 JJ=1,M                                                     F4500990
      DV1(JJ)=DV2(JJ)*(C(JJ)-Q(JJ))*DV1(JJ)                             F4501000
101   RE=RE+DV1(JJ)                                                     F4501010
      IM=0.0                                                            F4501020
      KAP=XO                                                            F4501030
      MIN=OM-KAP*KAP                                                    F4501040
      IF(MIN.LT.0.0) RETURN                                             F4501050
      IRA=M1+1                                                          F4501060
      IF(IRA.GT.M) GOTO 102                                             F4501070
      DO 103 I=IRA,M                                                    F4501080
      IF( ABS(OM-KAP*KAP).GT.MIN) GOTO 104                              F4501090
      IMIN=I                                                            F4501100
      MIN= ABS(OM-KAP*KAP)                                              F4501110
104   KAP=KAP+DELX                                                      F4501120
103   CONTINUE                                                          F4501130
102   CONTINUE                                                          F4501140
      IM=-DV2(IMIN)*3.1416D 00                                          F4501150
      RETURN                                                            F4501160
      END                                                               F4501170

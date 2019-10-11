      SUBROUTINE SELFEN(RR,RM,LAMK,                                     F4500010
     *  DELXM,DELXF,LAMI,LAML,LAMM,LAMF,                                F4500020
     *  MME4,MME5,MME1,MME2,R,RE,IM,RED,                                F4500030
     *  IMD,Q,C,DV1,AA1,AA2,AA3,AA4,                                    F4500040
     *  K,W,M,M1,F,F1,LV,W1,W2,LK,LI,KK1,FP,MP,KLL,IW1,IW2,YM,KM11)     F4500050
      INTEGER  K,W,M,M1,F,F1,LV,W2,W1,FP                                F4500060
      REAL*8 RR(K,K),RM(K,K) ,LAMK,AA1,AA2,AA3,AA4,                     F4500070
     *  DELXM,DELXF,LAMI,LAML,OM,RE1,IM1,Y1,Y2,Y3,Y4,YM,                F4500080
     *  LAMM(M),LAMF(F), A,B,SUM,SUM1,RE2,IM2,                          F4500090
     *  MME4(KK1,MP),MME5(KK1,MP),MME1(KK1,FP,MP),                      F4500100
     *  MME2(KK1,FP,MP),R(M),RE(F),IM(F),RED(F),                        F4500110
     *  IMD(F),Q(M),C(M),DV1(M),                                        F4500120
     *  RE3,RE4,IM3,IM4                                                 F4500130
      IF(W1*W2.EQ.1) WRITE(66,2) LAMK                                    F4500140
2     FORMAT(/' self-energy matrix',
     * ' is calculated at energy  ',E12.5)                                F4500160
3     FORMAT(1X,5I6/1X,4(E12.5,2X))                                     F4500190
      DO 6 I=1,F                                                        F4500200
      DO 7 I1=1,M                                                       F4500210
7     R(I1)=MME1(W2,I,I1)*MME1(W1,I,I1)*2                               F4500220
      OM=LAMK-LAMF(I)+LAML                                              F4500230
      CALL CQD1( SQRT(LAMM(M))- SQRT(LAMM(M-1)),                        F4500240
     * Q,OM,0.0,LAMM,C, SQRT(LAMM(M1+1)),DV1,R,                         F4500250
     * RE1,IM1,M,M1)                                                    F4500260
      RE(I)=RE1                                                         F4500270
      IM(I)=IM1                                                         F4500280
      IF(I.LE.F1) GOTO 8                                                F4500290
      A=2.* SQRT(LAMF(I))                                               F4500300
      RED(I-F1)=RE(I)*A                                                 F4500310
      IMD(I-F1)=IM(I)*A                                                 F4500320
8     CONTINUE                                                          F4500330
6     CONTINUE                                                          F4500340
      I1=F-F1                                                           F4500350
      IF(W.EQ.W1.AND.W.EQ.W2) write(66, 9)                                   F4500360
     * (RED(I),I=1,I1),(IMD(I),I=1,I1)                                  F4500370
9     FORMAT(10(2X,E10.4))                                              F4500380
      SUM=SIMPS1(RED,DELXF,F-F1)                                        F4500390
      SUM1=SIMPS1(IMD,DELXF,F-F1)                                       F4500400
      IF(F1.LE.0) GOTO 102                                              F4500410
      DO 10 I=1,F1                                                      F4500420
      SUM=SUM+RE(I)                                                     F4500430
10    SUM1=SUM1+IM(I)                                                   F4500440
102   B=1./DFLOAT((2*LV+1)*(2*LK+1))                                    F4500450
      RE1=SUM*B                                                         F4500460
      IM1=SUM1*B                                                        F4500470
      DO 16 I=1,F                                                       F4500480
      DO 17 I1=1,M                                                      F4500490
17    R(I1)=MME1(W2,I,I1)*MME2(W1,I,I1)                                 F4500500
      OM=LAMK-LAMF(I)+LAML                                              F4500510
      CALL CQD1( SQRT(LAMM(M))- SQRT(LAMM(M-1)),                        F4500520
     * Q,OM,0.0,LAMM,C, SQRT(LAMM(M1+1)),DV1,R,                         F4500530
     * RE2,IM2,M,M1)                                                    F4500540
      RE(I)=RE2                                                         F4500550
      IM(I)=IM2                                                         F4500560
      IF(I.LE.F1) GOTO 18                                               F4500570
      A=2.* SQRT(LAMF(I))                                               F4500580
      RED(I-F1)=RE(I)*A                                                 F4500590
      IMD(I-F1)=IM(I)*A                                                 F4500600
18    CONTINUE                                                          F4500610
16    CONTINUE                                                          F4500620
      I1=F-F1                                                           F4500630
      IF(W.EQ.W1.AND.W.EQ.W2) write(66, 9)                                   F4500640
     * (RED(I),I=1,I1),(IMD(I),I=1,I1)                                  F4500650
      SUM=SIMPS1(RED,DELXF,F-F1)                                        F4500660
      SUM1=SIMPS1(IMD,DELXF,F-F1)                                       F4500670
      IF(F1.LE.0) GOTO 202                                              F4500680
      DO 20 I=1,F1                                                      F4500690
      SUM=SUM+RE(I)                                                     F4500700
20    SUM1=SUM1+IM(I)                                                   F4500710
202   RE2=SUM*B                                                         F4500720
      IM2=SUM1*B                                                        F4500730
      RE3=0.0                                                           F4500740
      RE4=0.0                                                           F4500750
      IF(LV.GT.LK+LI.OR.LV.LT.                                          F4500760
     * IABS(LK-LI))  GOTO   23                                          F4500770
      DO 21 I=1,M                                                       F4500780
21    R(I )=MME4(W2,I)   *MME4(W1,I)*2.                                 F4500790
      OM=LAMI-LAMK   +LAML                                              F4500800
      CALL CQD1( SQRT(LAMM(M))- SQRT(LAMM(M-1)),                        F4500810
     * Q,OM,0.0,LAMM,C, SQRT(LAMM(M1+1)),DV1,R,                         F4500820
     * RE3,IM3,M,M1)                                                    F4500830
      RE3=-RE3*B                                                        F4500840
      IM3=-IM3*B                                                        F4500850
      DO 22 I=1,M                                                       F4500860
22    R(I )=MME4(W2,I)   *MME5(W1,I)                                    F4500870
      CALL CQD1( SQRT(LAMM(M))- SQRT(LAMM(M-1)),                        F4500880
     * Q,OM,0.0,LAMM,C, SQRT(LAMM(M1+1)),DV1,R,                         F4500890
     * RE4,IM4,M,M1)                                                    F4500900
      RE4=-RE4*B                                                        F4500910
      IM4=-IM4*B                                                        F4500920
23    SUM=RE1*AA1-RE2*AA2+RE3*AA3-RE4*AA4                               F4500930
      SUM1=IM1*AA1-IM2*AA2+IM3*AA3-IM4*AA4                              F4500940
      IF(W.NE.W1. OR.W.NE.W2)  GOTO 8748                                F4500950
      Y1=RE1*AA1                                                        F4500960
      Y2=RE2*AA2                                                        F4500970
      Y3=RE3*AA3                                                        F4500980
      Y4=RE4*AA4                                                        F4500990
      write(66, 8744) Y1,Y2,Y3,Y4                                            F4501000
8744  FORMAT('  BPEM. Bp. direct: ',E10.4,'     exch.',E10.4/             F4501010
     * '  BPEM. HAz. direct.: ',E10.4,'     exch.',E10.4)                  F4501020
8748  CONTINUE                                                          F4501030
      IF(W1.LE.KM11.OR.W2.LE.KM11) GOTO 8747                            F4501040
      IF( ABS(SUM).LE. ABS(YM)) GOTO 8745                               F4501050
      YM=SUM                                                            F4501060
      IW2=W2                                                            F4501070
      IW1=W1                                                            F4501080
8745  CONTINUE                                                          F4501090
8747  CONTINUE                                                          F4501100
      IF(W1.EQ.KLL. AND.W2.EQ.KLL ) write(66, 8746) IW2,IW1,YM               F4501110
8746  FORMAT('  max matrix element  RR(',I2,',',I2,')=',E10.4)          F4501120
      RR(W2,W1)=RR(W2,W1)+SUM                                           F4501130
      RM(W2,W1)=RM(W2,W1)+SUM1                                          F4501140
      RETURN                                                            F4501150
      END                                                               F4501170


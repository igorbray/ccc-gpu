      SUBROUTINE MAEL5(ALFA,BET,H,R1,
     * ME4,ME5,ME6,ME,R2,FUNC,R5,R6,PR,
     * EN,L,IT,NK,NM,M,LL,LI,LK,LM,
     * LV,NE5,KF)
      INTEGER   EN
      REAL*8 ALFA,BET,H,R1(NE5)     ,
     * ME4(  M),ME5(  M),ME6(  M),ME(1,M),
     * R2(NE5),FUNC(KF,NE5),R5(NE5),
     * R6(NE5),PR(1),A1,D1,B,  df
      CALL MATR(ALFA,BET,H,ME,R1,R2,FUNC,
     * R5,R6,
     * EN,1,M,NK,IT,L,NM,LV,NE5,KF)
      CALL IOT3(LK,LV,LI,0,0,0,A1)
      CALL IOT3(LM,LV,LL,0,0,0,D1)
      df=DFLOAT((2*LK+1)*
     * (2*LI+1)*(2*LL+1)*(2*LM+1))
      A1=A1*D1*2.*SQRT(df)
      DO 1 J=1,M
1     ME4(  J)=A1*ME(1,J)
      DO 2 J=1,M
2     ME5(  J)=0.0
      L1=IABS(LK-LL)
      IF(IABS(LK-LL).LE.IABS(LM-LI)) L1=
     * IABS(LM-LI)
      L2=LK+LL
      IF(LK+LL.GT.LM+LI) L2=LM+LI
      IF(L2.LT.L1) GOTO 4
      DO 3 J=L1,L2,2
      CALL MATR(ALFA,BET,H,ME,R1,R2,FUNC,
     * R5,R6,
     * EN,1,M,NK,L,IT,NM,J, NE5,KF)
      CALL IOT3(LK,J ,LL,0,0,0,A1)
      CALL IOT3(LM,J ,LI,0,0,0,D1)
      CALL IOT6(LK,LV,LI,LM,J,LL,B)
      PR(J)=A1*D1*B*(-1.D 00)**(LV+J)*
     * DFLOAT((2*LV+1)*2)
      df=DFLOAT((2*LK+1)*(2*LI+1)*
     * (2*LL+1)*(2*LM+1))
      A1=SQRT(df)
      DO 5 J1=1,M
5     ME5(  J1)=ME5(  J1)+PR(J)*ME(1,J1)*A1
3     CONTINUE
4     CONTINUE
      DO 6 J=1,M
6     ME6(  J)=2.*ME4(  J)-ME5(  J)
      RETURN
      END

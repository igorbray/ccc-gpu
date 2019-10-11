      SUBROUTINE MAEL4(ALFA,BET,H,R1,
     * ME1,ME2,ME3,ME,R2,FUNC,R5,R6,PR,
     * EN,NK,IT,F,NM,M,L,LL,LF,LM,LK,
     * LV,NE5,KF)
      INTEGER   EN,F
      REAL*8 ALFA,BET,H,R1(NE5)     ,
     * ME1(F,M),ME2(F,M),ME3(F,M),ME(F,M),
     * R2(NE5),FUNC(KF,NE5),R5(NE5),
     * R6(NE5),PR(1),A1,D1,B      , df
      CALL MATR(ALFA,BET,H,ME,R1,R2,FUNC,
     * R5,R6,
     * EN,M,F,L ,NM,NK ,IT,LV,NE5,KF)
      CALL IOT3(LK,LV,LF,0,0,0,A1)
      CALL IOT3(LL,LV,LM,0,0,0,D1)
      df=DFLOAT((2*LK+1)*
     * (2*LF+1)*(2*LL+1)*(2*LM+1))
      A1=A1*D1*2.*SQRT(df)
      DO 1 I=1,F
      DO 1 J=1,M
      JNN=(J+M*(I-1)-1)/F+1
      INN=J+M*(I-1)-(JNN-1)*F
1     ME1(I,J)=A1*ME(INN,JNN)
      DO 2 I=1,F
      DO 2 J=1,M
2     ME2(I,J)=0.0
      L1=IABS(LK-LM)
      IF(IABS(LK-LM).LE.IABS(LF-LL)) L1=
     * IABS(LF-LL)
      L2=LK+LM
      IF(LK+LM.GT.LF+LL) L2=LF+LL
      SUM=0.0
      IF(L2.LT.L1) GOTO 4
      DO 3 J=L1,L2,2
      CALL MATR(ALFA,BET,H,ME,R1,R2,FUNC,
     * R5,R6,
     * EN,F,M,L,IT,NK,NM,J, NE5,KF)
      CALL IOT3(LK,J ,LM,0,0,0,A1)
      CALL IOT3(LF,J ,LL,0,0,0,D1)
      CALL IOT6(LK,LV,LF,LL,J,LM,B)
      PR(J)=A1*D1*B*(-1.D 00)**(LV+J)*
     * DFLOAT((2*LV+1)*2)
      df=DFLOAT((2*LK+1)*(2*LF+1)*
     * (2*LL+1)*(2*LM+1))
      A1=SQRT(df)
      DO 5 I=1,F
      DO 5 J1=1,M
5     ME2(I,J1)=ME2(I,J1)+PR(J)*ME(I,J1)*A1
3     CONTINUE
4     CONTINUE
      DO 6 I=1,F
      DO 6 J=1,M
6     ME3(I,J)=2.*ME1(I,J)-ME2(I,J)
      RETURN
      END

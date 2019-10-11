c !!!
c     *CNST,
      SUBROUTINE MAEL3(EE,
     * BETA,
     * ALPHA,FUNC,H,R,MV,ML,YD,AU,QQ,
     * NCH,NES,HL,ND,LH,LP,  NN,
     *   NHL,      NWF,NPHS,LV,IH)
c !!!
c     *  NNL,
      real *8 simps1 ,  ddff
      INTEGER HLI,SI,NES(NCH),HL(NCH),ND(NCH),
     * LH(NCH),LP(NCH)
c !!!
c ,NNL(NCH)
      REAL*8FUNC(NWF,700),R(NN),YD(700),
     * EE(NPHS),F5,D,A,B,C,F,OM,CN,F4,QQ(IH),
     * AU(700),ML(IH,NPHS),MV(IH,NPHS),ER,BE,Y,A1,
     * CNST,BETA,ALPHA,H,  BET2,F1,F2,F3
      K1=0
      DO 118 I=1,NCH
      WRITE(6, 102) I,(QQ(J),J=1,IH)
102   FORMAT('    Channal       N   ',I2,/
     * 6X,'impulse:',5(6X,'Q=',F6.3,7X))
       WRITE(6, 2748)
 2748  FORMAT(14X,5('!',4X,'L',10X,'V',4X))
      NESI=NES(I)
      HLI=HL(I)
      NDI=ND(I)
      ILI=LH(I)
      ILM=LP(I)
      CALL IOT3(ILM,LV,ILI,0,0,0,CN)
      ddff=DFLOAT((2*ILM+1)*(2*ILI+1))
      CN=CN*SQRT(ddff)
      F1=FUNC(HLI,1)
      F2=FUNC(HLI,2)
      F3=FUNC(HLI,3)
      F4=FUNC(HLI,4)
      F5=FUNC(HLI,5)
      D=1./(12.*H)
      DO 103 J=1,NN
      IF(J.EQ.1) GOTO 104
      GOTO 105
104   C=-25.*F1+48.*F2-36.*F3+16.*F4-3.*F5
      F=F1
      GOTO 106
105   CONTINUE
      IF(J.EQ.NN) GOTO 107
      GOTO 108
107   C=25.*F5-48.*F4+36.*F3-16.*F2+3.*F1
      F=F5
      GOTO 106
108    IF(J.EQ.2) GOTO 109
      GOTO 110
109      C=-3.*F1-10.*F2+18.*F3-6.*F4+F5
       F=F2
       GOTO 106
110   IF(J.EQ.NN-1) GOTO 111
       GOTO  112
111   C=3*F5+10.*F4-18.*F3+6.*F2-F1
      F=F4
      GOTO 106
112   F1=F2
      F2=F3
      F=F4
      F3=F4
      F4=F5
      F5=FUNC(HLI,J+2)
      C=F1-8.*(F2-F4)-F5
106   C=C*D
103   YD(J)=C
      DO 118 JJ1=1,NESI
      K1=K1+1
      N1=NHL+K1
      OM=EE(K1)
       DO 1188 II11=1,IH
       Y=QQ(II11)
      DO 114 J=1,NN
      A=R(J)
      B=A/(ALPHA*A+BETA)
      ER=A
      CALL BESSEL(Y,ER,BE,LV)
      AU(J)=FUNC(HLI,J)*FUNC(N1,J)*BE*B*B
114    CONTINUE
      A=CN*SIMPS1(AU,H,NN)
      ML(II11,K1)=A
      IA=ILI*(ILI+1)-ILM*(ILM+1)
      A1=0.
      IF(IA.EQ.0) GOTO 500
      DO 200 IW=1,NN
      J1=LV-1
      ER=R(IW)
      CALL BESSEL(Y,ER,BE,J1)
      J1=LV+1
      CALL BESSEL(Y,ER,B ,J1)
      BE=BE/Y+Y*B
      AU(IW)=FUNC(HLI,IW)*FUNC(N1,IW)*BE*
     * R(IW)/(ALPHA*R(IW)+BETA)**2
200   CONTINUE
      A1=DFLOAT(IA)*SIMPS1(AU,H,NN)
500   CONTINUE
      F1=FUNC(N1 ,1)
      F2=FUNC(N1 ,2)
      F3=FUNC(N1 ,3)
      F4=FUNC(N1 ,4)
      F5=FUNC(N1 ,5)
      D=1./(12.*H)
      DO 10300  J=1,NN
      IF(J.EQ.1) GOTO 10400
      GOTO 10500
10400 C=-25.*F1+48.*F2-36.*F3+16.*F4-3.*F5
      F=F1
      GOTO 10600
10500 CONTINUE
      IF(J.EQ.NN) GOTO 10700
      GOTO 10800
10700 C=25.*F5-48.*F4+36.*F3-16.*F2+3.*F1
      F=F5
      GOTO 10600
10800  IF(J.EQ.2) GOTO 10900
      GOTO 11000
10900 C=-3.*F1-10.*F2+18.*F3-6.*F4+F5
       F=F2
       GOTO 10600
11000 IF(J.EQ.NN-1) GOTO 11100
       GOTO  11200
11100 C=3*F5+10.*F4-18.*F3+6.*F2-F1
      F=F4
      GOTO 10600
11200 F1=F2
      F2=F3
      F=F4
      F3=F4
      F4=F5
      F5=FUNC(N1 ,J+2)
      C=F1-8.*(F2-F4)-F5
10600 C=C*D
10300 AU(J)=C
      DO 300 IW=1,NN
      AU(IW)=FUNC(HLI,IW)*AU(IW)*R(IW)/
     * (ALPHA*R(IW)+BETA)
      J1=LV+1
      ER=R(IW)
      CALL BESSEL(Y,ER,BE,J1)
      J1=LV-1
      CALL BESSEL(Y,ER,B ,J1)
      BE=DFLOAT(LV+1)*Y*BE-DFLOAT(LV)/Y*B
      AU(IW)=AU(IW)*BE
300   CONTINUE
      A1=A1+SIMPS1(AU,H,NN)
      DO 400 IW=1,NN
      AU(IW)=FUNC(N1 ,IW)*YD(IW)*R(IW)/
     * (ALPHA*R(IW)+BETA)
      J1=LV+1
      ER=R(IW)
      CALL BESSEL(Y,ER,BE,J1)
      J1=LV-1
      CALL BESSEL(Y,ER,B ,J1)
      BE=DFLOAT(LV+1)*Y*BE-DFLOAT(LV)/Y*B
      AU(IW)=AU(IW)*BE
400   CONTINUE
      A1=(-A1+SIMPS1(AU,H,NN))*Y/DFLOAT(2*LV+1)
      MV(II11,K1)=A1*CN
 1188  CONTINUE
      WRITE(6, 117) K1,OM,(ML(J,K1),
     * MV(J,K1),J=1,IH)
117   FORMAT(1X,I3,2X,F7.4,1X,6('!',
     * 2(1X,F8.4,1X)))
 118  CONTINUE
      RETURN
      END

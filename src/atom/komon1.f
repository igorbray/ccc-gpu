      SUBROUTINE KOMON1(H,BETA,ALPHA,R,EE,FUNC,
     * PQN,EHL,KAPO,DKAP,PHS,F,Z,MT,NCH,HL,
     * EXPI,NES,ND,NNL,NN,LP,LH,NWF,NPHS,NHL,N0,NN5,ISTOP)
      INTEGER     BOO,PL,EXPI,HL(NCH),
     * PQN(NCH),LH(NCH),NNL(NCH),NES(NCH),
     * ND(NCH),LP(NCH),N0(NCH),ZNM
      REAL*8AR(15),R(NN),FUNC(NWF,700),
     * EHL(NCH),RCO,BETA,ALPHA,B,Z,X0,A,
     *KAPO(NCH),DKAP(NCH),EE(NPHS),
     * PHS(NPHS),F(705),XIS,XPHS,BETA1, C,
     * Z1,H1,H
      istop=0
      ZNM=NCH+1
      L7=0
      DO 413 J1=1,ZNM
      IF(L7.EQ.0)   READ(2 ) AR(1),AR(2),AR(3),AR(4),AR(5),
     * AR(6),AR(7),AR(8),AR(9),AR(10),
     * AR(11),AR(12),AR(13),AR(14),AR(15)
      if(l7.EQ.1)   READ(3 ) AR(1),AR(2),AR(3),AR(4),AR(5),
     * AR(6),AR(7),AR(8),AR(9),AR(10),
     * AR(11),AR(12),AR(13),AR(14),AR(15)
        IF(L7.EQ.2) READ(4 )(AR(II4),II4=1,15)
        IF(L7.EQ.3) READ(10)(AR(II4),II4=1,15)
        IF(L7.EQ.4) READ(8 )(AR(II4),II4=1,15)
        IF(L7.EQ.5) READ( 9)(AR(II4),II4=1,15)
       PL=int(((   AR(1))))
       N1=int(((   AR(3))))
       NF=int(((   AR(9))))
      IF(PL.EQ.1) GOTO 43
      IF(PL.NE.2) GOTO 42
      NF=int(((    AR(12)+AR(13))))
      GOTO 43
42    NF=1
43    continue
      IF(NF.NE.0) GOTO 44
      WRITE(66, 45) J1,PL,NF
45    FORMAT(//,' for J1=',I3,'  PL=',I3,
     * 1X,'  NF=',I3,'  --number of WF are in',
     * 'this file')
      ISTOP=1
      RETURN
44    BOO=1
      IF(J1.GT.1  ) GOTO 46
      IF(PL.EQ.1.OR.PL.EQ.3) GOTO 47
      WRITE(66, 48) PL,NHL,J1
48    FORMAT(2X,' PL=',I2,' NHL=',I2,' J1=',
     * I2,' entered WF of excited state ',
     * ' ')
      ISTOP=1
      RETURN
47    CONTINUE
      RCO=AR(2)
      NN=N1
      H=AR(4)
      BETA=AR(5)
      Z=AR(7)
      X0=-BETA*(10.+LOG(Z))
      A=X0+DFLOAT(  NN)*H
      ALPHA=(A-BETA*LOG(RCO))/RCO
      C=RCO
      DO 50 II=1,NN
      J=NN+1-II
      A=X0+DFLOAT(J)*H
1000  B=(ALPHA*C+BETA*LOG(C)-A)/(ALPHA*C+BETA)
      C=C*(1.-B)
      IF(ABS(B).GT.0.1D-06) GOTO 1000
50    R(J)=C
      WRITE(66, 51) R(1),RCO,NN,H,ALPHA,BETA,Z
51    FORMAT(' RMIN',
     * 24X,E15.8,/ , 1X,'RMAX',
     * 24X,E15.8,/ , 1X,'integration points',
     * 8X,I4, /, 1X,'step of integration',
     * 10X,F9.6,/ ,
     *  1X,'ALFA ',23X,F9.6,/ ,1X,'BETA',
     * 24X,F9.6,/ ,' nuclear charge ',19X,F4.1,6X/
     * 5X,'graund state configuration  :',
     * /9X,'N',9X,'L',9X,'Q')
      K1=0
      DO 60 K=1,NF
      CALL TAPE(F,XIS,XPHS,2,NN,IQ,2)
          II4=int((( XIS)))
          II3=int(((  F(2))))
          WRITE(66, 833) II4,II3,IQ
833      FORMAT(3(5X,I5))
      BOO=1
      DO 60 I=1,NCH
      IF(HL(I).NE.K) GOTO 6054
      IF(BOO.EQ.0) GOTO 58
      BOO=0
      K1=K1+1
      F(NN+2)=F(NN+1)
      DO 59 J=1,NN
59    FUNC(K1,J)=F(J+2)
58    HL(I)=-K1
      IF(EXPI.EQ.0) EHL(I)=F(1)
      PQN(I)=II4
      LH(I)=II3
      NNL(I)=IQ
6054  continue
60    CONTINUE
      WRITE(66, 6758)
6758  FORMAT(' wave functions of the graund state ',
     * /13X,'E',12X,'N',8X,'L')
      DO 6753 IO=1,NCH
      WRITE(66, 6754) IO,EHL(IO),PQN(IO),LH(IO)
6754  FORMAT(' channal',I2,2X,F9.4,6X,I2,8X,I1)
6753  CONTINUE
      REWIND 2
      DO 62 I=1,NCH
62    HL(I)=-HL(I)
      GOTO 4135
46    I=J1-1
      BOO=0
      IF(PL.EQ.2) GOTO 65
      WRITE(66, 66) J1,PL
66    FORMAT(2X,'J1=',I2,' PL=',I2,'reading',
     * ' WF of non-ground state')
      ISTOP=1
      RETURN
65    CONTINUE
      WRITE(66, 67)  I
67    FORMAT(/' Wave functions of excited states in the channel',I5)
      RCO=AR(2)
      H1=AR(4)
      BETA1=AR(5)
      Z1=AR(7)
       IF(ABS(RCO/R(NN)-1.).GT.0.1D-04.OR.
     * N1.NE.NN.OR.ABS(BETA1/BETA-1.).GT.
     *0.1D-04.OR.ABS(Z1/Z-1.).GT.0.1D-04) GOTO 200
      GOTO 69
200   CONTINUE
      WRITE(66, 70) RCO,N1,H1,BETA1,Z1
70    FORMAT(/,' Parameters of excited state',
     * ' : RMAX=',
     * F7.1,/,23X,'N    =',I7,2X,'H    =',
     * F7.3,/,22X,'BETA=',F7.4,/,' Z    =',I7,/)
      ISTOP=1
      RETURN
69    CONTINUE
      WRITE(66, 71)
71    FORMAT(5X,'E',7X,'N or KAP',4X,'L',
     *  6X,'face '   )
      N1=0
      II=I-1
      IF(II.LT.1) GOTO 572
      DO 72 I1=1,II
72    N1=N1+NES(I1)
572   K1=0
      IF(NES(J1-1  ).LT.NF) NF=NES(J1-1)
     *  -ND(J1-1)+int(((AR(12))))
      DO 75 K=1,NF
      IF(L7.EQ.1) CALL TAPE(F,XIS,XPHS,2,NN,IQ,3)
      IF(L7.EQ.2) CALL TAPE(F,XIS,XPHS,2,NN,IQ,4)
      IF(L7.EQ.3) CALL TAPE(F,XIS,XPHS,2,NN,IQ,10)
      IF(L7.EQ.4) CALL TAPE(F,XIS,XPHS,2,NN,IQ,8)
      IF(L7.EQ.5) CALL TAPE(F,XIS,XPHS,2,NN,IQ,9)
      IF(K1.LT.ND(I).AND.F(1).GT.0.) GOTO 74
      IF(K1.GE.ND(I).AND.F(1).LT.0.
     * .OR.K1.GE.NES(I)) GOTO 7554
      GOTO 76
74    ND(I)=K1
76    K1=K1+1
      N1=N1+1
      PHS(N1)=XPHS
      IF(K1.EQ.1) LP(I)=int(((F(2))))
      EE(N1)=F(1)-EHL(I)
      II4=NHL+N1
      DO 77 J=1,NN
77    FUNC( II4  ,J)=F(J+2)
      WRITE(66, 78) F(1),XIS ,LP(I ),XPHS
78    FORMAT(2X,2(F9.6,2X),I3,2X,F10.4)
      IF(K.NE.1) GOTO 79
      N0(I)=0
      IF(F(1).LT.0.) N0(I)=int(((XIS)))
79    continue
      IF(K1.EQ.ND(I)+1) KAPO(I)=XIS
      IF(K1.EQ.ND(I)+2) DKAP(I)=XIS -KAPO(I)
7554  continue
75    CONTINUE
          IF(L7.EQ.1) REWIND 3
          IF(L7.EQ.2) REWIND 4
          IF(L7.EQ.3) REWIND 10
          IF(L7.EQ.4) REWIND 8
4135             L7=L7+1
413       CONTINUE
      RETURN
      END

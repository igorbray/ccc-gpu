      SUBROUTINE MAIN1(ALFA,BET,H,Z,EPS,GAM,
     *TETA,R1,R2,R3,R4,R5,R6,R8,C3,ALA,ALA1,
     *ALA2,ALM3,MU,IS,NE,NE5,NE5IS,IQ,ISIGM,IN,
     *IL,IG,NU,ISIGM1,ISIGM2,KK,UD,KSUP,alpha,r0)
C  CALCULATION OF ONE DISCRETE SPECTRUM WAVE FUNCTION
C  FROM THE HARTREE-FOCK EQUATION (ORTHOG. BEFORE EXCHANGE)
      INTEGER IQ(IS),IN(IS),IL(IS),KK(NU),
     * ISIGM1(NU),ISIGM2(NU),UD(IS)
      REAL*8  ALFA,BET,H,Z,RO1,RO2,EPS,
     *GAM(NU ),TETA(IS),B,C,D,P,T,
     *  U,V,GAMA,TAU1,TAU2,ALM1,ALM2,
     *TETMAX,A1,B1,A2,OMEGA,ALM3(IS),
     * AKAP,DKAP,ALM4,ALM5,ALM6,ALM7,
     *C1,C5,ALA(IS),ALA1(IS),ALA2(IS),
     * C3(IS),R1(NE5 ),R2(NE5 ),R3(NE5 ),
     *R4(NE5 ),R5(NE5 ),
     *R6(NE5 ),R8(NE5IS),MU ,FN1,FN2,one,pone,temp
      one = 1d0
      pone=0.1d0
      ALM4=0.
      ALM5=0.
      ALM6=0.
      ALM7=0.
      DO 1 J=1,IS
      ALA(J)=0.
      ALA1(J)=0.
      ALA2(J)=0.
C     ALM3(J)=0.
   1  C3(J)=1.
c$$$      write(6, 8851)
8851  FORMAT('  VERSION 2. ORTOGONALIZATION, EXCHANGE ITERATIONS')
      KF=1
      KF1=1
      KSI1=0
      KSI =0
      DO 2 J=1,NE
   2  R3(J+1)=0.
      CALL DIRECT(ALFA,BET,H,GAM,R1,R5,
     * R3,R8,NE,NE5,NE5IS,IS  ,IQ,IG,NU,
     * ISIGM1,ISIGM2,KK,ISIGM,KSUP)
C  Patch for polarization potential put in by Igor 22/4/92
      DO 6 J=1,NE
         r = r1(j+2)
 6       R3(J)=R3(J  )-Z + r * polpot(r,alpha,r0)
      MMM=(ISIGM-1)*NE5
      DO 9221 J=1,NE5
9221  R2(J)=R8(MMM+J)
      IN1=IN(ISIGM)
      IQ1=IQ(ISIGM)
      IL1=IL(ISIGM)
      T=BET**2/4.
      P=ALFA*BET
      D=DFLOAT( IL1*(IL1+1))
      IAY=1
      DO 10 J=1,NE
      B=2.*R1(J+2)*R3(J)
      C=(ALFA*R1(J+2)+BET)**2
   10 R3(J)=(D+B+(P*R1(J+2)+T)/C)/C
      alm6_1 = 1.0
      alm6_2 = 2.0
  222 DO 11 J=1,NE5
   11 R4(J)=0.
      FN1=0.0
      FN2=0.0
      NPL=0
      NMI=0
      IY=IAY
      IF(ABS(MU-1.D 00).GT.1.D-04) GOTO 5700
      if (abs((alm6_1+alm6_2)/alm6_1).lt.1e-2) print*,'infinite loop?',
     >   alm6_1,alm6_2
      CALL EXCHA1(ALFA,BET,H,R1,R2,
     * R4,R5,R6,R8,ALM3,GAM,NE,NE5,
     * ISIGM,IG,IQ,KK,ISIGM1,ISIGM2,IS)
5700  IRA=NE+2
      DO 1254 J=1,NE5
 1254 R5(J)=R2(J)
9000  IY=IAY
       ISIS=IS-1
       DO 7826 IJJ=1,ISIS
 7826    C3(IJJ)=1.
  333 D=DFLOAT( IL1)+0.5
      M1=0
      M2=0
      IF(FN1*FN2.NE.0.0.AND.
     * ABS(FN1-FN2).LT.1.D-12) GOTO 4052
      GOTO 4053
4052  write(6, 4054) FN1,FN2
4054  FORMAT(' CAN NOT FIND TO ZEROS ',
     * ' FOR POTENTIAL P(X) WITHIN THE INTERVAL',2(3X,E22.15))
4053  CONTINUE
      IX=1
      DO 17 J=1,NE
      C= R3(J)+R2(1)*(R1(J+2) /(ALFA*R1(J+2)+BET))**2
      R6(J+1)=C*H**2/12.
      B=D*C
      IF(B.GE.0.0) GOTO 100
      GOTO 101
  100 D=C
      GOTO 17
  101 D=C
      GOTO (102,103),IX
  102 M1=J
      IX=2
      GOTO 17
  103 M2=J
   17 CONTINUE
      IF(M2.EQ.0.AND.M1.EQ.0) GOTO 205
      GOTO 206
  205 FN1=R2(1)
      R2(1)=R2(1)/2.
      IF(FN1*FN2.NE.0.0) R2(1)=(FN1+FN2)/2.
      GOTO 333
  206 IF(M1.NE.0.AND.M2.EQ.0) GOTO 207
      GOTO 208
  207 write(6,    18)
   18 FORMAT (20H ENERGY INCREASING   )
      FN2=R2(1)
      R2(1)=R2(1)*1.5
      IF(FN1*FN2.NE.0.0) R2(1)=(FN1+FN2)/2.
      GOTO 333
  208 DO 19 J=2,NE
      R6(J-1)=0.
      R2(J)=0.
      IA=J
      IF(R6(J+1).LT.1.) GOTO 20
   19 CONTINUE
   20 IRA=M1-1
      IF(IRA.LT.IA) GOTO 209
      DO 21 J=IA,IRA
      A1=1.-R6(J)
      B1=2.+10.*R6(J+1)
      A2=1.-R6(J+2)
      OMEGA=R4(J+2)+10.*R4(J+1)+R4(J)
      IF(J.EQ.IA) GOTO 22
      GOTO 23
   22 V=0.
   24 U=V
      V=A2/(B1-A1*V)
      IF(ABS(U-V).GT.1.D-07) GOTO 24
      U=0.
      GOTO 25
   23 B=B1-A1*V
      V=A2/B
      U=(A1*U+OMEGA)/B
   25 R6(J)=V
   21 R2(J+1)=U
  209 CONTINUE
      GAMA=1.
      RO1=U+V*GAMA
      KSI=0
      IRA=M2-1
      IF(IRA.LT.M1) GOTO 300
      gamam = 0.0
      DO 26 J=M1,IRA
      A1=1.-R6(J)
      B1=2.+10.*R6(J+1)
      A2=1.-R6(J+2)
      OMEGA=R4(J+2)+10.*R4(J+1)+R4(J)
      R2(J+2)=GAMA
      V=RO1
      RO1=GAMA
      if (abs(gama).gt.gamam) gamam = gama
      GAMA=((B1*RO1)-(A1*V)-OMEGA)/A2
      IF(RO1*GAMA.LT.0.0) KSI=KSI+1
   26 KSI1=KSI
      if (abs(gama).gt.1e1) print*,'KSI,gama:',
     >   ksi,gama
      if (abs(gama).gt.1e2) stop 'no convergence; reduce nnbtop?'
  300 CONTINUE
      IA=IN1-IL1-1
      IF(KSI.EQ.IA) GOTO 27
      IF(KSI.GT.IA) NMI=NMI+1
      IF(KSI.LT.IA) NPL=NPL+1
      temp = dble(KSI-IN1+IL1+1)
      A1=SIGN(pone , temp)
      IF((KSUP.LE.2.AND.IA.GT.3.OR.
     * KSUP.EQ.3).AND.IABS(KSI-IA).LE.1.
     * AND.NMI*NPL.NE.0) A1=A1*0.1
c      WRITE(3,9722) KSI,IA,M1,M2,NE,R2(1)
9722  FORMAT(' DIFF. NUMBER OF ZEROS  KSI=',I3,
     * ' INSTEAD OF =',I3, /, ' M1=',I4,' M2=',I4,
     * ' NE=',I4,' **** PREVIOUS ENERGY =',E15.8)
 444  R2(1)=R2(1)*(1.+A1/DFLOAT(IN1))
c      WRITE(3,9723)R2(1),ALM3(ISIGM)
9723  FORMAT('   NEW  = ',E15.8,
     * ' OFFDIAG. TERM =',E15.8)
      GOTO 333
   27 IRA=NE-1
      DO 28 IJ=1,IRA
      J=IRA-IJ+1
      R6(J+3)=0.
      R2(J+4)=0.
      IA=J
      IF(R6(J+1).LT.1.) GOTO 29
   28 CONTINUE
   29 IF(M2.GT.IA) GOTO 301
      DO 30 IJ=M2,IA
      J=IA-IJ+M2
      A1=1.-R6(J)
      B1=2.+10.*R6(J+1)
      A2=1.-R6(J+2)
      OMEGA=R4(J+2)+10.*R4(J+1)+R4(J)
      IF(J.EQ.IA) GOTO 31
      GOTO 32
   31 V=0.
   33 U=V
      V=A1/(B1-A2*V)
      IF(ABS(V-U).GT.1.D-07) GOTO 33
      U=0.
      GOTO 34
   32 B=B1-A2*V
      V=A1/B
      U=(A2*U+OMEGA)/B
   34 R6(J+2)=V
   30 R2(J+3)=U
  301 CONTINUE
      A1=(GAMA-(RO1*V+U))/H
      IF(KSUP.EQ.3) GOTO 70288
      IF(ABS(A1).LT.ABS(EPS*GAMA)) GOTO 107
      GOTO 70289
70288 CONTINUE
      IF(ABS(A1).LT.ABS(0.001*GAMA)) GOTO 107
70289 CONTINUE
c      WRITE(3,9724)A1,GAMA
9724  FORMAT(/'  DERIVATIVE JUMP=',E15.8,
     * ' GAMA(0.001)=',E15.8)
      IF(A1.GT.0.0) GOTO 35
      GOTO 36
   35 TAU1=A1
      ALM1=R2(1)
   36 IF(A1.LT.0.0) GOTO 37
      GOTO 38
   37 TAU2=A1
      ALM2=R2(1)
   38 GOTO(104,105,106,107),IY
  104 ISI=int(SIGN(one,A1))
      IY=2
      GOTO 39
  105 IF(ISI.NE.int(SIGN(one,
     *  A1))) GO TO 106
   39 A1=SIGN(pone,A1)
      IF(MOD(IN1-IL1,2).NE.0) A1=-A1
      GOTO 444
  106 IY=3
      R2(1)=(TAU2*ALM1-TAU1*ALM2)/(TAU2-TAU1)
      IF(ABS(TAU1+TAU2).GT.(TAU1-TAU2)/2.)
     * R2(1)=(ALM1+ALM2)/2.
      IF(ABS(R2(1)-ALM1).LT.1.D-20.OR.
     * ABS(R2(1)-ALM2).LT.1.D-20) write(6,
     *  4703) R2(1),ISIGM
4703  FORMAT(' ',88('*')/' DERIVATIVE JUMP CAN NOT ',
     * 'BE ELIMINATED AT E=',E18.11,
     * ' FOR ',I2,' SUBSHELL')
c      WRITE(3,9725)TAU1,TAU2, ALM1,ALM2,R2(1)
9725  FORMAT(' TAU1=',E15.8,
     *  '   TAU2=',E15.8,/' ALM1=',E15.8,
     * '    ALM2=',E15.8,/'    ****E=',E15.8)
      GOTO 333
  107 IRA=M1-1
      DO 40 IJ=1,IRA
      J=IRA-IJ+1
   40 R2(J+2)=R2(J+3)*R6(J)+R2(J+1)
      IRA=M2-1
      IRA1=NE-1
      IF(IRA.GT.IRA1) GOTO 302
      DO 41 J=IRA,IRA1
      IF(ABS(R2(J+2)).LT.1.D-30) R2(J+2)=0.0
   41 R2(J+3)=R2(J+2)*R6(J+3)+R2(J+4)
  302 CONTINUE
      IF(KSUP.EQ.3) GOTO 90071
      ALM6=0.0
      DO 99452 J=1,NE
      B=R2(J+2)-R5(J+2)
      IF(ABS(B).GT.ABS(ALM6)) ALM6=B
99452 CONTINUE
      ALM6=4.*ALM6
      GOTO 99453
99451 CONTINUE
      CALL ACCOBM(R2,R5,ALM4,ALM5,ALM6,
     * NE,NE5,-1,KF )
      KF=KF+1
99453 CONTINUE
      alm6_2 = alm6_1
      alm6_1 = alm6
      IF(ABS(ALM6).GT.EPS) GOTO 222
      IF(KSUP.NE.1) GOTO 90071
      MSK=NE5*(ISIGM-1)
      DO 2303 IW=1,NE5
2303  R5(IW)=R8(IW+MSK)
      T=(R2(1)-R5(1))/R2(1)
      DO 3377 J=1,NE
      B=R5(J+2)-R2(J+2)
      IF(ABS(R2(J+2)).GT.1.)  B=B/R2(J+2)
      IF(ABS(B).GT.ABS(T)) T=B
3377  R5(J+2)=R2(J+2)
      TETA(ISIGM)=T
      GOTO 90072
7903  ALM4=ABS(ALM6  )
      ALM5=ABS(1.-R5(1)/R2(1))
      IAY=4
      IF(ALM4.GT.ALM5) GOTO 5060
      TETA(ISIGM)=ALM5
      GOTO 5061
 5060 TETA(ISIGM)=ALM4
 5061 CONTINUE
90071 ALM3(IS)=0.0
      C=0.
      D=0.
      DO 5389 J=1,NE
      P=D
      D=C
      C=(R2(J+2) *R1(J+2) /(ALFA*R1(J+2)+BET))**2
      A1=5.*C
      B=8.*D
      A1=(A1+B-P)*H/12.
5389  ALM3(IS)=ALM3(IS)+A1
      IF(ABS(MU-1.D 00).GT.1.D-04) GOTO 92001
      CALL OFFDI1(ALFA,BET,H,R1,R2,R4,
     * R6,R8,ALA1,ALA2,ALM3,TETA,
     * C3,UD,NE,IS,NE5,ISIGM,IN,IL,KF1)
c      WRITE(3,887)R2(1),TETA(ISIGM)
      KF1=KF1+1
      IAY=1
      ISIS=IS-1
      DO 52 J=1,ISIS
      IF(ABS(TETA(J)).GT.EPS/10.)GOTO 9000
  52  CONTINUE
92001 CONTINUE
      CALL ACCOBM(R2,R5,ALM4,ALM5,ALM6,
     * NE,NE5,-1,KF)
      KF=KF+1
      C=ABS(1.-R5(1)/R2(1))
      IF(C.GT.ABS(ALM6)) ALM6=C
      KF1=1
      IF(ABS(MU-1.D 00).GT.1.D-04) GOTO 9010
      IF(ABS(ALM6).GT.EPS) GOTO 222
 9010  CONTINUE
90072 R5(2)=0.
887   FORMAT('  *****  ENERGY =',E15.8,/,
     * 6X,' MISMATCH =',E15.8,'  **************')
      C=0.
      D=0.
      DO 53 J=1,NE
      P=D
      D=C
      C=(R2(J+2) *R1(J+2) /(ALFA*R1(J+2)+BET))**2
      A1=5.*C
      B=8.*D
      A1=(A1+B-P)*H/12.
   53 R5(2)=R5(2)+A1
      R5(1)=R2(1)
      WRITE( 4,80020)ISIGM,R5(1),TETA(ISIGM),M1,M2
80020 FORMAT( '  ISIGM=',I3,' E=',
     * E15.8,' TETA=',E15.8,'M1=',I4,'  M2=',I4)
      IF(KSUP.NE.3) RETURN
      IRA=NE+2
      DO 54 IW=3,IRA
   54 R5(IW)=R2(IW)
      DO 55 J=1,NE
   55 R5(J+2)=R5(J+2)/SQRT(R5(2))
      R5(NE+3)=R5(2)
      R5(2)=DFLOAT( IL(IS))
      R5(1)=-R5(1)
      MMM=(ISIGM-1)*NE5
      DO 2005 I100=1,NE5
2005  R8(MMM+I100)=R5(I100)
      RETURN
      END
      SUBROUTINE MAIN2(ALFA,BET,H,Z,AKAP,EPS,TETA,
     *R1,R2,R3,R4,R5,R6,R8,ALA,ALA1,ALA2,ALM3,C3,
     *GAM,MU,IS,NE,NE5,NE5IS,IL,KK,ISIGM1,ISIGM2,
     *IN,IQ,NU,IG,ISIGM,UD,ksup,alpha,r0)
C   CALCULATION OF THE CONTINUOUS SPECTRUM
C   WAVE FUNCTION (ORTHOG. BEFORE EXCHANGE)
      REAL*8  ALFA,BET,H,Z,RO1,RO2,EPS,
     *GAM(140),TETA(IS),B,C,D,P,T,
     *  U,V,GAMA,TAU1,TAU2,ALM1,ALM2,
     *TETMAX,A1,B1,A2,OMEGA,ALM3(IS),
     * AKAP,DKAP,ALM4,ALM5,ALM6,ALM7,C1,
     *C5,ALA(IS),ALA1(IS),ALA2(IS),AKVINT,
     * C3(IS),R1(NE5 ),R2(NE5 ),R3(NE5 ),
     *R4(NE5),R5(NE5),R6(NE5 ),FA,R8(NE5IS),MU
      INTEGER IL(IS),KK(9),ISIGM1(9),IN(IS),IQ(IS),
     * ISIGM2( 9),UD(IS)
      IRA=IS-1
c$$$      write(6, 88513)
88513 FORMAT('  VERSION 2. ORTHOGONALIZATION, EXCHANGE ITERATIONS')
      DO 1 J=1,IRA
      ALA(J)=0.0
      ALA1(J)=0.0
      ALA2(J)=0.0
C     ALM3(J)=0.0
      C3(J)=1.0
    1 CONTINUE
      KF=1
      IN1=IN(ISIGM)
      IL1=IL(ISIGM)
      IQ1=IQ(ISIGM)
      KF1=1
      ALM4=0.0
      ALM5=0.0
      ALM6=0.0
      C1=0.0
      IRA=NE+2
      DO 301 J=1,IRA
      R2(J)=0.0
      R3(J)=0.0
  301 CONTINUE
      CALL DIRECT(ALFA,BET,H,GAM,R1,R5,R3,R8,
     *NE,NE5,NE5IS,IS-1,IQ,IG,NU,ISIGM1,ISIGM2,
     *KK,ISIGM,ksup)
C  Patch for polarization potential put in by Igor 22/4/92
      DO 302 J=1,NE
         r = r1(j+2)
 302     R3(J)=R3(J  )-Z + r * polpot(r,alpha,r0)
c$$$      DO 302 J=1,NE
c$$$      R3(J)=R3(J)-Z
c$$$  302 CONTINUE
      D=DFLOAT( IL1*(IL1+1))
      P=ALFA*BET
      T=BET**2/4.
      DO 2 J=1,NE
      B=2.0*R1(J+2)*R3(J)
      C=(ALFA*R1(J+2)+BET)**2
      B=B*MU
      R3(J)=(D+B-AKAP**2*R1(J+2)**2+(P*R1(J+2)+T)/C)/C
    2 CONTINUE
  444 IRA=NE+2
      DO 397 J=1,IRA
397   R4(J)=0.0
      IF(ABS(MU-1.D 00).GT.1.D-04) GOTO 5700
      CALL EXCHA1(ALFA,BET,H,R1,R2,R4,R5,R6,R8,
     *ALM3,GAM,NE,NE5,ISIGM,IG,IQ,KK,ISIGM1,
     *ISIGM2,IS)
 5700 CONTINUE
      DO 3973 J=1,IRA
 3973 R5(J)=R2(J)
9000  CONTINUE
      D=DFLOAT(   IL1)+ 0.5
       DO 2579 I33=1,IS
 2579  C3(I33)=1.
      P=H**2/12.
      M1=0
      DO 4 I3=1,NE
      J=NE-I3+1
      C=R3(J)
      R6(J+1)=C*P
      B=C*D
      IF(B.LT.0.0) M1=J
      D=C
    4 CONTINUE
      DO 5 J=2,NE
      R6(J-1)=0.0
      R2(J)=0.0
      IA=J
      IF(R6(J+1).LT.1.0) GOTO 120
    5 CONTINUE
  120 IRA=M1-1
      IF(IA.GT.IRA) GOTO 501
      DO 21 J=IA,IRA
      A1=1.-R6(J)
      B1=2.+10.*R6(J+1)
      A2=1.-R6(J+2)
      OMEGA=R4(J+2)+10.*R4(J+1)+R4(J)
      IF(J.EQ.IA) GOTO 6
      B=B1-A1*V
      V=A2/B
      U=(A1*U+OMEGA)/B
      GOTO 7
    6 V=0.0
  100 U=V
      V=A2/(B1-A1*V)
      IF(ABS(V-U).GT.1.D-7)GOTO 100
      U=0.0
    7 CONTINUE
      R6(J)=V
      R2(J+1)=U
   21 CONTINUE
  501 CONTINUE
      GAMA=1.0
      RO1=U+V*GAMA
      DO 8 J=M1,NE
      A1=1.-R6(J)
      B1=2.+10.*R6(J+1)
      A2=1.-R6(J+2)
      OMEGA= R4(J+2)+10.*R4(J+1)+R4(J)
      R2(J+2)=GAMA
      V=RO1
      RO1=GAMA
      GAMA=(B1*RO1-A1*V-OMEGA)/A2
    8 CONTINUE
  555 IRA=M1-1
      DO 306 I3=1,IRA
      J=IRA-I3+1
      R2(J+2)=R2(J+3)*R6(J)+R2(J+1)
  306 CONTINUE
      A1=0.0
      GOTO 9001
9002  CONTINUE
      CALL ACCOBM(R2,R5,ALM4,ALM5,ALM6,
     *NE,NE5,1,KF)
      KF=KF+1
      KF1=1
      IF(ABS(MU-1.D 00).GT.1.D-04)  GOTO 99991
      IF(ABS( ALM6)  .GT.EPS) GOTO 444
      GOTO 99991
9001  ALM3(IS)=1.
      IF(ABS(MU-1.D 00).GT.1.D-04)GOTO 3972
      CALL OFFDI1(ALFA,BET,H,
     * R1,R2,R4,R6,R8,ALA1,ALA2,ALM3,TETA,C3,
     * UD,NE,IS,NE5,ISIGM,IN,IL,KF1)
 3972 KF1=KF1+1
      IRA=IS-1
      DO 308 J=1,IRA
      IF(ABS(TETA(J)).GT.EPS/10.)GOTO 9000
  308 CONTINUE
      GOTO 9002
99991 IA=0
      TETMAX=0.0
      IRA=IS-1
      DO 309 J=1,IRA
      IA=IA+IQ(J)
  309 CONTINUE
      OMEGA=Z-DFLOAT(  IA)
      CALL PHASE(AKAP,C,OMEGA,R1,R2,IL1,NE,NE5,IM1)
      KSI=IM1
C******* Changes to the old version are made below (14/08/97) *******
      ALM1=R1(NE+1)
      ALM2=R1(NE+2)
      U=R2(NE+1)*SQRT(ALM1/(ALFA*ALM1+BET))
      V=R2(NE+2)*SQRT(ALM2/(ALFA*ALM2+BET))
      OMEGA=OMEGA*MU
C************* new subroutine ASYM is inserted here **************
      CALL ASYM(AKAP,OMEGA,IL1,ALM1,U,ALM2,V,1.D-5,KSI,FA,TETMAX)
C******** It produces the phase FA and amplitude TETMAX **********
C** corresponding to the asymptotic ampl/sqrt(k)*sin(...+phase) **
      C=SQRT(ABS(MU)/3.141592D0)/TETMAX
      DO 310 J=1,NE
      R5(J+2)=R2(J+2)*C
  310 CONTINUE
      R5(NE+4)=FA
      R5(NE+3)=C
      R5(1)=AKAP**2/ABS(MU)
      R5(2)=DFLOAT(  IL(ISIGM))
      MMM=(IS-1)*NE5
      DO 2009 I100=1,NE5
2009  R8(MMM+I100)=R5(I100)
c$$$      write(6,3333) R5(NE+4)
 3333 FORMAT(' PHASE =',25X,E15.8)
      RETURN
      END

      SUBROUTINE OFFDI1(ALFA,BET,H,R1,R2,R4,
     *R5,R8,ALA1,ALA2,ALM3,TETA,C3,UD,
     *NE,IS,NE5,ISIGM,IN,IL,KF)
      include 'paratom.f'
      REAL*8 ALFA,BET,H,R1(NE5),R2(NE5),R4(NE5),P1,
     *R5(NE5),ALA1(IS),ALA2(IS),R8(nsizeis*nsize),
     *ALM3(IS),TETA(IS),C3(IS),TETMAX,C,D,P,A1
      INTEGER UD(IS),IN(IS),IL(IS)
      IL1=IL(ISIGM)
      IN1=IN(ISIGM)
      P1=H*H/12.
      IRA=IS-1
      DO 30 J=1,IRA
      TETMAX=0.0
      IF(.NOT.(IL1.EQ.IL(J).AND.
     *IN1.NE.IN(J))) GOTO 301
      IF(UD(J).NE.1) GOTO 301
      IA=J
      MMM=(IA-1)*NE5
      DO 2002 I100=1,NE5
2002  R5(I100)=R8(MMM+I100)
      C=0.0
      D=0.0
      R5(NE+2)=C
      DO 10 IW=1,NE
      P=D
      D=C
      C=R2(IW+2)*R5(IW+2)*(R1(IW+2)/
     *(ALFA*R1(IW+2)+BET))**2
      C=C/SQRT(ALM3(IS))
      A1=(5.*C+8.*D-P)*H/12.
      TETMAX=A1+TETMAX
   10 CONTINUE
      WRITE( 4,987)J,TETMAX,ALM3(J)
987   FORMAT(' OVERLAP INTEGRAL FOR ',
     * I3,' SUBSHELL =',E15.8,'  OFFDIAG.PARAM. ',E15.8)
      C3(J)=ALM3(J)+TETMAX
      IF(KF.LT.2)  GOTO 6779
      C3(J)=(ALM3(J)*ALA1(J)-
     * ALA2(J)*TETMAX)/(ALA1(J)-TETMAX)
      IF(TETMAX*ALA1(J).GT.0.0.AND.
     * (ABS(ALA1(J)).GT.0.01.OR.
     * ABS(TETMAX).GT.0.01)) C3(J)=
     * ALM3(J)+MIN(ABS(2./TETMAX),ABS(ALA1(J)/(TETMAX-
     * ALA1(J))))*TETMAX
6779  ALA2(J)=ALM3(J)
      ALA1(J)=TETMAX
      ALM3(J)=C3(J)
      DO 6791 IW=1,NE
6791  R4(IW+1)=R4(IW+1)+(ALM3(J)-
     * ALA2(J))/R5(2)*R5(IW+2)*P1*
     * (R1(IW+2)/(ALFA*R1(IW+2)+BET))**2
301   TETA(J)=TETMAX
30    CONTINUE
      RETURN
      END

      SUBROUTINE EXCHA1(ALFA,BET,H,R1,
     *R2,R4,R5,R6,R8,ALM3,GAM,NE,NE5,ISIGM,
     *IG,IQ,KK,ISIGM1,ISIGM2,IS)
c
      include 'paratom.f'
c
      REAL*8 ALFA,BET,H,R1(NE5),R2(NE5),P,R8(nsize*nsizeis),
     *R4(NE5),R5(NE5),R6(NE5),ALM3(max10),GAM(NSIZEqf),GAMA
      INTEGER KK(NSIZEqf),ISIGM1(NSIZEqf),ISIGM2(NSIZEqf),IQ(max10)
c$$$      REAL*8 ALFA,BET,H,R1(NE5),R2(NE5),P,R8(4200),
c$$$     *R4(NE5),R5(NE5),R6(NE5),ALM3(1),GAM(1),GAMA
c$$$      INTEGER KK(1),ISIGM1(1),ISIGM2(1),IQ(1)
      IF(IG.LT.1) RETURN
      DO 3 IW=1,IG
      IF(ISIGM1(IW).EQ.ISIGM) GOTO 81
      IF(ISIGM2(IW).EQ.ISIGM) GOTO 82
      IA=0
      GOTO 83
81    IA=ISIGM2(IW)
      GOTO 83
82    IA=ISIGM1(IW)
83    CONTINUE
      IF(IA.EQ. 0 ) GOTO 110
      MMM=(IA    -1)*NE5
      DO 2001 I100=1,NE5
2001  R5(I100)=R8(MMM+I100)
      DO 50 JJ=1,NE
      R6(JJ+1)=R5(JJ+2)*R2(JJ+2)
   50 CONTINUE
      CALL POTATOM(ALFA,BET,H,R1,R6,NE,KK(IW),NE5)
      GAMA=2.*GAM(IW)/(R5(2)*DFLOAT(IQ(ISIGM)))
      DO 51 JJ=1,NE
      R4(JJ+1)=R4(JJ+1)+ GAMA*R6(JJ+1)*R5(JJ+2)
   51 CONTINUE
  110 CONTINUE
    3 CONTINUE
      ISIS=IS-1
      DO 210 J=1,ISIS
      MM=(J-1)*NE5
      DO 211 JJ=1,NE5
 211  R5(JJ)=R8(MM+JJ)
      DO 212 IW=1,NE
  212 R4(IW+1)=R4(IW+1)+ALM3(J)/R5(2)*
     * R1(IW+2)*R5(IW+2)
 210  CONTINUE
      P=H*H /12.
      DO 305 J=1,NE
      R4(J+1)=P*R1(J+2)*R4(J+1)/(ALFA*R1(J+2)+BET)**2
  305 CONTINUE
      RETURN
      END
      SUBROUTINE ASYM(AK,ZC,LOR,RL,PL,RR,PR,EPS,KSI,PHA,AMP)
C****************** This subroutine calculates the ********************
C*** asymptotic behaviour of the radial wave function and determines **
C* the phase and amplitude AMP/sqrt(k)*sin(kr+Z/k*ln(2kr)-l*Pi/2+PHA) *
C
C        Parameters transferred to ASYM are
C     AK - momentum
C     ZC - core charge (times mass of the projectile, 
C          <0 for positively charged particles)
C     LOR - orbital momentum
C     RL,PL,RR,PR - radii and wave functions to start from (RL<RR)
C     EPS - accuracy (asymptotic phase required to be stable to EPS/100) 
C     KSI - no. of nodes of the w.f. left of RR 
C
C       ASYM returnes the phase and the amplitude, PHA and AMP
C
      implicit REAL*8 (a-h,o-z)
c$$$      write(6,1)
 1    format('CALCULATION OF THE PHASE AND AMPLITUDE BY ASYM SUBROUTIE')
      PI=0.31415926D+01
      DR=RR-RL
      IF(DR.GT.0D0) GOTO 20
      write(6,10)
 10   format('RR must be greater than RL!!!')
      GOTO 1000
 20   CONTINUE
C******** Seeing whether the step should be made smaller ********
C**** (if AK*DR>1/32, i.e., less than 100 points per half wave) ***
      NIN=3.2D1*AK*DR+1
      DR=DR/NIN
c$$$      write(6,30) RL,RR,NIN,DR
 30   format('Interval RL=',D12.5,' RR=',D12.5,' divided into NIN=',
     * I4,' of DR=',D12.5)
C******** useful constants for Numerov integration ***********
      HA=DR*DR/12.D0
      HB=5.D0*HA
      AK2=AK*AK
      DLL1=DFLOAT(LOR*(LOR+1))
      R1=RL
      P1=PL
      R2=RR
      P2=PR
      IF (NIN.EQ.1) GOTO 50
C******* Calculation of P1,P2 by "progonka", initial values ********
      U=PL
      V=0.D0
      R1=RL
      R=RL+DR
      R2=R+DR
      RHS1=(DLL1/R1-2.D0*ZC)/R1-AK2
      RHS=(DLL1/R-2.D0*ZC)/R-AK2
      RHS2=(DLL1/R2-2.D0*ZC)/R2-AK2
C******* Calculation of P1,P2 by "progonka", loop of NIN-1 steps ******
      DO 60 NPR=2,NIN
      A1=1.D0-HA*RHS1
      B=1.D0+HB*RHS
      A2=1.D0-HA*RHS2
      DEN=2.D0*B-A1*V
      U=A1*U/DEN
      V=A2/DEN
      R1=R
      R=R2
      R2=R2+DR
      RHS1=RHS
      RHS=RHS2
      RHS2=(DLL1/R2-2.D0*ZC)/R2-AK2
 60   CONTINUE
      P1=U+V*PR
      R2=RR
      P2=PR
C******* end of "progonka" *******
 50   CONTINUE
c$$$      write(6,70) R1,P1,R2,P2
 70   format('Integration of the asymptotic radial w.f. starts with',/,
     * 'P(',D15.8,')=',D15.8,', P(',D15.8,')=',D15.8)
      KSC=0
      PHA=0.D0
      RA=RR
      KSIRR=KSI
      AMP=0.0D0
C********* Integration with DR stepsize beyond the R2 point ***********
C** terminated when the amplitude/phase are found or KSMAX nodes met **
      KSMAX=500
      R3=R2+DR
      RHS1=(DLL1/R1-2.D0*ZC)/R1-AK2
      RHS2=(DLL1/R2-2.D0*ZC)/R2-AK2
      RHS3=(DLL1/R3-2.D0*ZC)/R3-AK2
 80   CONTINUE
c**** The Numerov formula is written this way to keep errors small **** 
      P3=((2.D0*P2-P1)+2.D0*HB*RHS2*P2+HA*RHS3*P1)/(1.D0-HA*RHS1)
C**** Calculation of the amplitude near maximum using semiclassics ****
      IF((P3-P2)*(P2-P1).GT.0.D0) GOTO 180
      PM1=SQRT(ABS(RHS1))
      PM2=SQRT(ABS(RHS2))
      PM3=SQRT(ABS(RHS3))
      P22=P2*SQRT(PM2)
      P13=(P3*SQRT(PM3)-P1*SQRT(PM1))/(2.D0*DR)
      P13=P13/(PM2-0.5D0*ZC/(AK*R2)**3)
      AMP=SQRT(P22*P22+P13*P13)
c      write(6,190) R2,AMP,P2
 190  format('Amplitude at R=',D18.10,' is AMP=',D18.10,', P=',D18.10)
 180  CONTINUE
      R1=R2
      R2=R3
      R3=R3+DR
      RHS1=RHS2
      RHS2=RHS3
      RHS3=(DLL1/R3-2.D0*ZC)/R3-AK2
      P1=P2
      P2=P3
      IF(P1*P2.LT.0.D0) KSI=KSI+1
C*********** Printing the radial wave function *************
c      write(6,90) R2,P2
 90   format(3D18.10)
C**** Below we determine the asymptotic parameters (AMP, PHA) ****
C**** To determine them by comparison with the semiclassical  ****
C****     behaviour one must ensure that E-U>0, and mF<p^3    ****
      EU=-RHS2
      FRC=ABS(RHS3-RHS1)/(2.D0*DR)
      P3=(SQRT(ABS(EU)))**3
      IF(EU.LE.0.D0.OR.FRC.GE.P3) GOTO 110
      IF (KSC.EQ.1) GOTO 120
      KSC=1
c$$$      write(6,130) R2
 130  format('Semiclassical approximation is valid beyond R=',D12.5)
 120  CONTINUE
C*********** Checking for the node ***********
      IF(P1*P2.GT.0.D0) GOTO 110
      RND=(P2*R1-P1*R2)/(P2-P1)
      TKR=2.D0*AK*RND
C*** Matching with the semiclassical phase with up to 1/r^2 terms ***
      PHSC=AK*RND+ZC/AK*LOG(TKR)+(DLL1+ZC*ZC/AK2)/TKR
      PH2=KSI*PI-PHSC+0.5D0*PI*DFLOAT(LOR)
      PHSC=PHSC-ZC*(DLL1-1.D0+ZC**2/AK2)/(AK*TKR**2)
      PH1=KSI*PI-PHSC+0.5D0*PI*DFLOAT(LOR)
      DPH=ABS(PHA-PH1)
c      write(6,150) RND,PH1
 150  format('R=',D18.10,' PH1=',D18.10)
      RA=RND
      PHA=PH1
      IF(DPH.LT.0.1D-1*EPS.AND.AMP.NE.0.D0) GOTO 160
 110  CONTINUE
      KSNEW=KSI-KSIRR
      IF(KSNEW.LE.KSMAX) GOTO 80
      write(6,100) KSMAX,R2,EPS
 100  format('PHASE NOT STABLE TO EPS/100 AFTER ',I3,' NODES,
     * R=',D12.5,', EPS:',e10.3)
 160  continue 
c$$$      write(6,170) PHA,R2,KSI,AMP
 170  format('PHA=',D15.8,' at R=',D15.8,
     * ' (KSI=',I4,' nodes)',', AMP=',D15.8)
      RETURN
 1000 STOP
      END

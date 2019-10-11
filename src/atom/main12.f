      SUBROUTINE MAIN1(ALFA,BET,H,Z,EPS,GAM,
     *TETA,R1,R2,R3,R4,R5,R6,R8,C3,ALA,ALA1,
     *ALA2,ALM3,MU,IS,NE,NE5,NE5IS,IQ,ISIGM,IN,
     *IL,IG,NU,ISIGM1,ISIGM2,KK,UD,KSUP)
C
C  Calculation of one wave function
C  Solution of the Hartree-Fock equation    
C
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
     *R6(NE5 ),R8(NE5IS),MU ,FN1,FN2
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
      write(6, 8851)
8851  FORMAT('  BEPCàü 2. OPTOÉOHAãàáAñàü,CAMOCOÉãACOBAHàE')
      KF=1
      KF1=1
      KSI1=0
      KSI =0
      DO 2 J=1,NE
   2  R3(J+1)=0.
      CALL DIRECT(ALFA,BET,H,GAM,R1,R5,
     * R3,R8,NE,NE5,NE5IS,IS  ,IQ,IG,NU,
     * ISIGM1,ISIGM2,KK,ISIGM,KSUP)
      DO 6 J=1,NE
   6  R3(J)=R3(J  )-Z
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
  222 DO 11 J=1,NE5
   11 R4(J)=0.
      FN1=0.0
      FN2=0.0
      NPL=0
      NMI=0
      IY=IAY
      IF(DABS(MU-1.D 00).GT.1.D-04) GOTO 5700
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
     * DABS(FN1-FN2).LT.1.D-12) GOTO 4052
      GOTO 4053
4052  write(6, 4054) FN1,FN2
4054  FORMAT(' HE HAXOÑüTCü 2 KOPHü ',
     * 'îìHKñàà P(X) B èPOMEÜìTKE',2(3X,E22.15))
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
   18 FORMAT (20H ùHEPÉàü ìBEãàóàãACú)
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
      IF(DABS(U-V).GT.1.D-07) GOTO 24
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
      DO 26 J=M1,IRA
      A1=1.-R6(J)
      B1=2.+10.*R6(J+1)
      A2=1.-R6(J+2)
      OMEGA=R4(J+2)+10.*R4(J+1)+R4(J)
      R2(J+2)=GAMA
      V=RO1
      RO1=GAMA
      GAMA=((B1*RO1)-(A1*V)-OMEGA)/A2
      IF(RO1*GAMA.LT.0.0) KSI=KSI+1
   26 KSI1=KSI
  300 CONTINUE
      IA=IN1-IL1-1
      IF(KSI.EQ.IA) GOTO 27
      IF(KSI.GT.IA) NMI=NMI+1
      IF(KSI.LT.IA) NPL=NPL+1
      A1=DSIGN(1.D-01 ,DFLOAT(KSI-IN1+IL1+1))
      IF((KSUP.LE.2.AND.IA.GT.3.OR.
     * KSUP.EQ.3).AND.IABS(KSI-IA).LE.1.
     * AND.NMI*NPL.NE.0) A1=A1*0.1
c      WRITE(3,9722) KSI,IA,M1,M2,NE,R2(1)
9722  FORMAT(' PAáH. KOã.KOPHEâ  KSI=',I3,
     * ' Ñ. Å.=',I3, /, ' M1=',I4,' M2=',I4,
     * ' NE=',I4,' **** ùHEPÉàü ÅõãA =',E15.8)
 444  R2(1)=R2(1)*(1.+A1/DFLOAT(IN1))
c      WRITE(3,9723)R2(1),ALM3(ISIGM)
9723  FORMAT('   CTAãA = ',E15.8,
     * ' HEÑàAÉOHAã. óãEH=',E15.8)
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
      IF(DABS(V-U).GT.1.D-07) GOTO 33
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
      IF(DABS(A1).LT.DABS(EPS*GAMA)) GOTO 107
      GOTO 70289
70288 CONTINUE
      IF(DABS(A1).LT.DABS(0.001*GAMA)) GOTO 107
70289 CONTINUE
c      WRITE(3,9724)A1,GAMA
9724  FORMAT(/'  CKAó. èPOàá.=',E15.8,
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
  104 ISI=IFIX(SNGL(DSIGN(1.D 00,A1)))
      IY=2
      GOTO 39
  105 IF(ISI.NE.IFIX(SNGL(DSIGN(1.D 00,
     *  A1)))) GO TO 106
   39 A1=DSIGN(1.D-01,A1)
      IF(MOD(IN1-IL1,2).NE.0) A1=-A1
      GOTO 444
  106 IY=3
      R2(1)=(TAU2*ALM1-TAU1*ALM2)/(TAU2-TAU1)
      IF(DABS(TAU1+TAU2).GT.(TAU1-TAU2)/2.)
     * R2(1)=(ALM1+ALM2)/2.
      IF(DABS(R2(1)-ALM1).LT.1.D-20.OR.
     * DABS(R2(1)-ALM2).LT.1.D-20) write(6,
     *  4703) R2(1),ISIGM
4703  FORMAT(' ',88('*')/' HE èPOàCXOÑàT ',
     * 'CòàBAHàü èPOàáBOÑHOâ èPà E=',E18.11,
     * ' B ',I2,' OÅOãOóKE')
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
      IF(DABS(R2(J+2)).LT.1.D-30) R2(J+2)=0.0
   41 R2(J+3)=R2(J+2)*R6(J+3)+R2(J+4)
  302 CONTINUE
      IF(KSUP.EQ.3) GOTO 90071
      ALM6=0.0
      DO 99452 J=1,NE
      B=R2(J+2)-R5(J+2)
      IF(DABS(B).GT.DABS(ALM6)) ALM6=B
99452 CONTINUE
      ALM6=4.*ALM6
      GOTO 99453
99451 CONTINUE
      CALL ACCOBM(R2,R5,ALM4,ALM5,ALM6,
     * NE,NE5,NE5IS,-1,KF )
      KF=KF+1
99453 CONTINUE
      IF(DABS(ALM6).GT.EPS) GOTO 222
      IF(KSUP.NE.1) GOTO 90071
      MSK=NE5*(ISIGM-1)
      DO 2303 IW=1,NE5
2303  R5(IW)=R8(IW+MSK)
      T=(R2(1)-R5(1))/R2(1)
      DO 3377 J=1,NE
      B=R5(J+2)-R2(J+2)
      IF(DABS(R2(J+2)).GT.1.)  B=B/R2(J+2)
      IF(DABS(B).GT.DABS(T)) T=B
3377  R5(J+2)=R2(J+2)
      TETA(ISIGM)=T
      GOTO 90072
7903  ALM4=DABS(ALM6  )
      ALM5=DABS(1.-R5(1)/R2(1))
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
      IF(DABS(MU-1.D 00).GT.1.D-04) GOTO 92001
      CALL OFFDI1(ALFA,BET,H,R1,R2,R4,
     * R6,R8,ALA,ALA1,ALA2,ALM3,TETA,
     * C3,UD,NE,IS,NE5,ISIGM,IN,IL,KF1)
c      WRITE(3,887)R2(1),TETA(ISIGM)
      KF1=KF1+1
      IAY=1
      ISIS=IS-1
      DO 52 J=1,ISIS
      IF(DABS(TETA(J)).GT.EPS/10.)GOTO 9000
  52  CONTINUE
92001 CONTINUE
      CALL ACCOBM(R2,R5,ALM4,ALM5,ALM6,
     * NE,NE5,NE5IS,-1,KF)
      KF=KF+1
      C=DABS(1.-R5(1)/R2(1))
      IF(C.GT.DABS(ALM6)) ALM6=C
      KF1=1
      IF(DABS(MU-1.D 00).GT.1.D-04) GOTO 9010
      IF(DABS(ALM6).GT.EPS) GOTO 222
 9010  CONTINUE
90072 R5(2)=0.
887   FORMAT('  *****   ùHEPÉàü =',E15.8,/,
     * 6X,' HEBüáKA =',E15.8,'  ÑàAÉOHAãàáAñàü')
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
     * E15.8,' TETA=',E15.8,'M1=',I3,'  M2=',I3)
      IF(KSUP.NE.3) RETURN
      IRA=NE+2
      DO 54 IW=3,IRA
   54 R5(IW)=R2(IW)
      DO 55 J=1,NE
   55 R5(J+2)=R5(J+2)/DSQRT(R5(2))
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
     *IN,IQ,NU,IG,ISIGM,UD)
C   BõóàCãEHàE OÑHOâ BOãHOBOâ îìHKñàà
C   CèãOòHOÉO  CèEKTPA
      REAL*8  ALFA,BET,H,Z,RO1,RO2,EPS,
     *GAM(140),TETA(IS),B,C,D,P,T,
     *  U,V,GAMA,TAU1,TAU2,ALM1,ALM2,
     *TETMAX,A1,B1,A2,OMEGA,ALM3(IS),
     * AKAP,DKAP,ALM4,ALM5,ALM6,ALM7,C1,
     *C5,ALA(IS),ALA1(IS),ALA2(IS),
     * C3(IS),R1(NE5 ),R2(NE5 ),R3(NE5 ),
     *R4(NE5),R5(NE5),R6(NE5 ),FA,R8(NE5IS),MU
      INTEGER IL(IS),KK(9),ISIGM1(9),IN(IS),IQ(IS),
     * ISIGM2( 9),UD(IS)
      IRA=IS-1
      write(6, 88513)
88513 FORMAT('  BEPCàü 2. OPTOÉOHAãàáAñàü, CAMOCOÉãACOBAHàE')
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
     *KK,ISIGM)
      DO 302 J=1,NE
      R3(J)=R3(J)-Z
  302 CONTINUE
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
      IF(DABS(MU-1.D 00).GT.1.D-04) GOTO 5700
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
      IF(DABS(V-U).GT.1.D-7)GOTO 100
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
     *NE,NE5,NE5IS,1,KF)
      KF=KF+1
      KF1=1
      IF(DABS(MU-1.D 00).GT.1.D-04)  GOTO 99991
      IF(DABS( ALM6)  .GT.EPS) GOTO 444
      GOTO 99991
9001  ALM3(IS)=1.
      IF(DABS(MU-1.D 00).GT.1.D-04)GOTO 3972
      CALL OFFDI1(ALFA,BET,H,
     * R1,R2,R4,R6,R8,ALA,ALA1,ALA2,ALM3,TETA,C3,
     * UD,NE,IS,NE5,ISIGM,IN,IL,KF1)
 3972 KF1=KF1+1
      IRA=IS-1
      DO 308 J=1,IRA
      IF(DABS(TETA(J)).GT.EPS/10.)GOTO 9000
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
      ALM2=R1(NE+1)
      ALM1=ALM2
      V=R2(NE+1)*DSQRT(ALM1/(ALFA*ALM1+BET))
      OMEGA=OMEGA*MU
      A2=(DFLOAT( IL1*(IL1+1))/ALM1-2.*
     *OMEGA)/ALM1-AKAP**2
      TAU2=(DABS(A2))**0.25
      T=V *TAU2
      B1=A2
      ALM1=R1(NE+2)
      KF=0
      IF((ALM1-ALM2).LT.0.2/AKAP) GOTO 14
      ALM1=(ALM1+ALM2)/2.0
      R2(NE+2)=AKVINT(ALM1,R1,R2,NE+1,NE+2,NE5)
      write(6,   1111)
 1111 FORMAT( 20H òAÉ ÑEãàTCü èOèOãAM)
   14 CONTINUE
      U=R2(NE+2)*DSQRT(ALM1/(ALFA*ALM1+BET))
      A2=(DFLOAT(IL1*(IL1+1))/ALM1-2.*
     * OMEGA)/ALM1-AKAP**2
      TAU2=(DABS(A2))**0.25
      TAU2=TAU2*U
      J=0
      ALM2=ALM1-ALM2
      P=ALM2**2/12.
      IA=1
      R5(NE+3)=1.0
      IC4=0
      C5=0.0
  115 A1=B1
      B1=A2
      TAU1=T
      T=TAU2
      RO1=ALM1
      GAMA=V
      V=U
      ALM1=ALM1+ALM2
      A2=(DFLOAT( IL1*(IL1+1))/ALM1-2.*
     *OMEGA)/ALM1-AKAP**2
      TAU2=(DABS(A2))**0.25
      U=((B1*P*10.0+2.0)*V-(1.0-A1*P)*
     *  GAMA)/(1.0-A2*P)
      TAU2=TAU2*U
      IF( R2(NE+2)*U.LT.0.0.AND.KF.EQ.0.AND.ALM1.
     *LE.R1(NE+2)  )KSI=KSI-1
      KF=1
      IF(U*V.GE.0.0) GOTO 31
C     D=U/DSQRT(ALM1/(ALFA*ALM1+BET))
C     C=V/DSQRT(RO1/(ALFA*RO1+BET))
      C=(U*RO1-V*ALM1)/(U-V)
      KSI=KSI+1
      J=J+1
      IF(IC4.NE.1) GOTO 31
      FA=DFLOAT(KSI)*3.14159-AKAP*C-OMEGA/
     *AKAP*DLOG(2.0*C*AKAP)+3.14159/
     *2.*DFLOAT( IL1) -((OMEGA/AKAP)**2+
     *DFLOAT(IL1*(IL1+1)))/(2.*AKAP*C)
      IF(J.NE.6) GOTO 31
      J=0
      C5=R5(NE+3)
      R5(NE+3)=FA
   31 CONTINUE
      IF(IC4.EQ.1) GOTO 125
      IF( (TAU2-T)*(T-TAU1).GT.0.0) GOTO 115
      B=(GAMA-U)/(GAMA+U-V)
      C=B-1.0
      D=B+1.0
      B=DABS((GAMA*C*B+U*D*B)/2.0-D*C*V)
      C=TETMAX-B
      IF(IA.NE.4) GOTO 16
      TETMAX=B
      IA=0
   16 CONTINUE
      IA=IA+1
      J=0
      IF(DABS(C/10.0).GT.EPS) GOTO 115
      write(6,   2222)  KSI
 2222 FORMAT(' ', I5,'KOPHEâ ')
  125 IC4=1
      IF(J.NE.0) GOTO 115
      IF(R5(NE+3).LT.0.001.OR.KSI.GT.150.AND.
     *R5(NE+3).NE.1.) GOTO 130
      IF(DABS(C5-R5(NE+3)).GT.EPS) GOTO 115
  130 write(6,   3333)  ALM1,KSI,R5(NE+3)
 3333 FORMAT( 8H R KOHEó,
     *E15.8,5X,I3,7H KOPHEâ/' îAáA   =',25X,E15.8)
      C=DSQRT(DABS(MU)/(3.1416*AKAP))/TETMAX
      DO 310 J=1,NE
      R5(J+2)=R2(J+2)*C
  310 CONTINUE
      R5(NE+4)=R5(NE+3)
      R5(NE+3)=C
      R5(1)=AKAP**2/DABS(MU)
      R5(2)=DFLOAT(  IL(ISIGM))
      MMM=(IS-1)*NE5
      DO 2009 I100=1,NE5
2009  R8(MMM+I100)=R5(I100)
      RETURN
      END
      SUBROUTINE OFFDI1(ALFA,BET,H,R1,R2,R4,
     *R5,R8,ALA,ALA1,ALA2,ALM3,TETA,C3,UD,
     *NE,IS,NE5,ISIGM,IN,IL,KF)
      REAL*8 ALFA,BET,H,R1(NE5),R2(NE5),R4(1),P1,
     *R5(NE5),ALA(IS),ALA1(IS),ALA2(IS),R8(5000),
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
      C=C/DSQRT(ALM3(IS))
      A1=(5.*C+8.*D-P)*H/12.
      TETMAX=A1+TETMAX
   10 CONTINUE
      WRITE( 4,987)J,TETMAX,ALM3(J)
987   FORMAT('  àHTEÉPAã èEPEKP. Ñãü',
     * I3,' OÅOã.=',E15.8,'  HEÑ.óã. ',E15.8)
      C3(J)=ALM3(J)+TETMAX
      IF(KF.LT.2)  GOTO 6779
      C3(J)=(ALM3(J)*ALA1(J)-
     * ALA2(J)*TETMAX)/(ALA1(J)-TETMAX)
      IF(TETMAX*ALA1(J).GT.0.0.AND.
     * (DABS(ALA1(J)).GT.0.01.OR.
     * DABS(TETMAX).GT.0.01)) C3(J)=
     * ALM3(J)+DMIN1(DABS(2./TETMAX),DABS(ALA1(J)/(TETMAX-
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
      REAL*8 ALFA,BET,H,R1(NE5),R2(NE5),P,R8(4200),
     *R4(NE5),R5(NE5),R6(NE5),ALM3(1),GAM(1),GAMA
      INTEGER KK(1),ISIGM1(1),ISIGM2(1),IQ(1)
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


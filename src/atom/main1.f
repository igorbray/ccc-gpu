      SUBROUTINE MAIN1(ALFA,BET,H,Z,EPS,GAM,
     *TETA,R1,R2,R3,R4,R5,R6,R8,C3,ALA,ALA1,
     *ALA2,ALM3,MU,IS,NE,NE5,NE5IS,IQ,ISIGM,IN,
     *IL,IG,NU,ISIGM1,ISIGM2,KK,UD,KSUP,alpha,r0)
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
c$$$      write(*,*) ' MAIN1'
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
C
         
C         
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
   11 R4(J  )=0.
      NPL=0
      NMI=0
      FN1=0.0
      FN2=0.0
      IY=IAY
      IF(DABS(MU-1.D 00).GT.1.D-04) GOTO 5700
      CALL EXCHAN(ALFA,BET,H,R1,R2,
     * R4,R5,R6,R8,ALM3,GAM,NE,NE5,
     * ISIGM,IG,IQ,KK,ISIGM1,ISIGM2)
5700  IRA=NE+2
      DO 15 J=1,IRA
   15 R5(J)=R2(J)
  333 D=DFLOAT( IL1)+0.5
      M1=0
      M2=0
      IF(FN1*FN2.NE.0.0.AND.
     * DABS(FN1-FN2).LT.1.D-12) GOTO 4052
      GOTO 4053
4052  write(6, 4054) FN1,FN2
4054  FORMAT(' HE HAXOdIATCIA 2 KOPHIA ',
     * 'FYHKTSII P(X) B interval',2(3X,E22.15))
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
   18 FORMAT ('energy has been increased')
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
      A1=DFLOAT(KSI-IN1+IL1+1)
      A1=DSIGN(1.D 00 , A1)
      IF((KSUP.LE.2.AND.IA.GT.3.OR.
     * KSUP.EQ.3).AND.IABS(KSI-IA).LE.1.
     *  AND.NMI*NPL.NE.0) A1=A1*0.1
      WRITE( 3,9722) KSI,IA,M1,M2,NE,R2(1)
9722  FORMAT(' PAzH. KOL.KOPHEII  KSI=',I3,
     * ' d. b.=',I3, /, ' M1=',I4,' M2=',I4,
     * ' NE=',I4,' **** energy was =',E15.8)
      IF(IY.NE.1) GOTO 107
 444  R2(1)=R2(1)*(1.+A1/DFLOAT(IN1))
      WRITE( 3,9723)R2(1),ALM3(ISIGM)
9723  FORMAT('   CTALA = ',E15.8,
     * ' HEdIAgOHAL. chLEH=',E15.8)
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
      WRITE( 3,9724)A1,GAMA
9724  FORMAT(/'  CKAch. chPOIz.=',E15.8,
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
   39 A1=DSIGN(7.D-02,A1)
      IF(MOD(IN1-IL1,2).NE.0) A1=-A1
      GOTO 444
  106 IY=3
      R2(1)=(TAU2*ALM1-TAU1*ALM2)/(TAU2-TAU1)
      IF(DABS(TAU1+TAU2).GT.(TAU1-TAU2)/2.)
     * R2(1)=(ALM1+ALM2)/2.
      IF(DABS(R2(1)-ALM1).LT.1.D-20.OR.
     * DABS(R2(1)-ALM2).LT.1.D-20) write(6,
     *  4703) R2(1),ISIGM
4703  FORMAT(' ',88('*')/' HE pPOICXOdIT ',
     * 'CSHIBAHIIA pPOIzBOdHOII pPI E=',E18.11,
     * ' in ',I2,' shell')
      WRITE( 3,9725)TAU1,TAU2, ALM1,ALM2,R2(1)
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
C     write(6, 89099) NE,M1,M2, (R2(III),III=1,200,10)
89099 FORMAT(' R2',3I3  /200(6(2X,E12.5)/))
      IF(KSUP.EQ.3) GOTO 99451
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
      IF(DABS(ALM6).GT.EPS*0.1D 00) GOTO 222
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
90071 ALM4=DABS(ALM6  )
      ALM5=DABS(1.-R5(1)/R2(1))
      IAY=4
      IF(ALM4.GT.ALM5) GOTO 5060
      TETA(ISIGM)=ALM5
      GOTO 5061
 5060 TETA(ISIGM)=ALM4
 5061 CONTINUE
      ALM3(IS)=0.0
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
      CALL OFFDIA( ALFA,BET,H,R1,R2,
     * R5,R8,ALA,ALA1,ALA2,ALM3,TETA,
     * C3,UD,NE,IS,NE5,ISIGM,IN,IL,KF1)
      WRITE( 3,887)R2(1),TETA(ISIGM)
      KF1=KF1+1
      IAY=1
      KF=1
      DO 52 J=1,IS
      IF(DABS(TETA(J)).GT.EPS) GOTO 222
  52  CONTINUE
90072 R5(2)=0.
887   FORMAT('  *****   eHEPgIIA =',E15.8,/,
     * 6X,' HEBIAzKA =',E15.8,'  dIAgOHALIzATSIIA')
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





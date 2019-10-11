      SUBROUTINE MAIN2(ALFA,BET,H,Z,AKAP,EPS,TETA,
     *R1,R2,R3,R4,R5,R6,R8,ALA,ALA1,ALA2,ALM3,C3,
     *GAM,MU,IS,NE,NE5,NE5IS,IL,KK,ISIGM1,ISIGM2,
     *IN,IQ,NU,IG,ISIGM,UD,ksu,alpha,r0)
C
C  Calculation of one continuum wave function 
C
      REAL*8  ALFA,BET,H,Z,RO1,RO2,EPS,
     *GAM(140),TETA(IS),B,C,D,P,T,
     *  U,V,GAMA,TAU1,TAU2,ALM1,ALM2
      real*8 TETMAX,A1,B1,A2,OMEGA,ALM3(IS),
     * AKAP,DKAP,ALM4,ALM5,ALM6,ALM7,C1,
     *C5,ALA(IS),ALA1(IS),ALA2(IS),
     * C3(IS),R1(NE5 ),R2(NE5 ),R3(NE5 ),
     *R4(NE5),R5(NE5),R6(NE5 ),FA,R8(NE5IS),MU
      INTEGER IL(IS),KK(9),ISIGM1(9),IN(IS),IQ(IS),
     * ISIGM2( 9),UD(IS),ksu

      IRA=IS-1
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
     *kk,isigm,   ksu  )
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
      IF(DABS(MU-1.D 00).GT.1.D-04) GOTO 5700
      CALL EXCHAN(ALFA,BET,H,R1,R2,R4,R5,R6,R8,
     *ALM3,GAM,NE,NE5,ISIGM,IG,IQ,KK,ISIGM1,
     *ISIGM2)
5700  D=DFLOAT(   IL1)+ 0.5
      DO 3971 J=1,IRA
3971  R5(J)=R2(J)
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
      CALL ACCOBM(R2,R5,ALM4,ALM5,ALM6,
     *NE,NE5,NE5IS,1,KF)
      KF=KF+1
      IF(DABS(MU-1.D 00).GT.1.D-04)  GOTO 99991
      IF(DABS(ALM6).GT.EPS*1.D-01) GOTO 444
      ALM3(IS)=1.
      IF(DABS(MU-1.D 00).LT.1.D-04)
     * CALL OFFDIA(ALFA,BET,H,
     * R1,R2,R5,R8,ALA,ALA1,ALA2,ALM3,TETA,C3,
     * UD,NE,IS,NE5,ISIGM,IN,IL,KF1)
      KF1=KF1+1
      KF=1
      IRA=IS-1
      DO 308 J=1,IRA
      IF(DABS(TETA(J)).GT.EPS) GOTO 444
  308 CONTINUE
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
c$$$      write(6,   1111)
c$$$ 1111 FORMAT( 20H òAÉ ÑEãàTCü èOèOãAM)
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
c$$$      write(6,   2222) KSI
c$$$ 2222 FORMAT(' ', I5,'KOPHEâ ')
  125 IC4=1
      IF(J.NE.0) GOTO 115
      IF(R5(NE+3).LT.0.001.OR.KSI.GT.150.AND.
     *R5(NE+3).NE.1.) GOTO 130
      IF(DABS(C5-R5(NE+3)).GT.EPS) GOTO 115
 130  continue
c$$$      write(6,   3333) ALM1,KSI,R5(NE+3)
c$$$ 3333 FORMAT( 8H R KOHEó,
c$$$     *E15.8,5X,I3,7H KOPHEâ/' îAáA   =',25X,E15.8)
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

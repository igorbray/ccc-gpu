      SUBROUTINE FCHF(Z,GAM,DIR,R,BET,EPS,H,
     *TETA,ALA,ALA1,ALA2,ALM3,C3,R1,R2,R3,R4,R5,
     *R6,R7,R8,R88,AKAP,DKAP,MU,K1,K2,IS,IN,IL,
     *IQ,S1,SS1,L1,NU,IG,KK,ISIGM1,ISIGM2,
     *IP,PP,NE,NE5,NE5IS,KT,UD,KSU1,alpha,r0)
c
      include 'paratom.f'
c      
      REAL*8  ALFA,BET,H,Z,RO1,RO2,EPS,GAM( NU),
     *   TETA(IS),B,C,D,P,T,PRE,ROT,AA1,T2,AKSI,
     *   U,V,GAMA,TAU1,TAU2,ALM1,ALM2,TETMAX,
     *   A1,B1,A2,OMEGA,ALM3(IS),MU, UU,VV,
     *   AKAP,DKAP,ALM4,ALM5,ALM6,ALM7,C1,
     *   ALFA1,RN,RO1N,BET1,H1,R88(max10),
     *   C5,ALA(IS),ALA1(IS),ALA2(IS),
     *   C3(IS),R1(NE5 ),R2(NE5 ),R3(NE5 ),
     *   R4(NE5 ),R5(NE5 ), SOST,XN1,
     *   R6(NE5 ),R7(KT ),A(15),R,XNOBL,XPHS,
     *   XIS,QQ1,QQ2,QQ3,QQ4,QQ5,QQ6,ROH,
     *   DIR(max10),R8(NE5IS),RO2N
      INTEGER IN(IS),IL(IS),IQ(IS),
     *   KK(NU),ISIGM1(NU),
     *   IP(max10),PP(max10),SS1,ISIGM2(NU),CH,UD(max10),S1(max10)
      izon = 0
      izon1 = 0 ! defined by Igor, 10/4/2011
      IF(L1.EQ.111) GOTO 111
      DO 5 I=1,IS
5     IP(I)=1
      CH=1
      I1=IS
      S1(1)=SS1
      DO 1111 I=1,IS
      JJA=2*IL(I)+1
      IF(S1(1).LE.99) JJA=JJA*2
      DO 1111 J=I1,IS
      JJ=2*IL(J)+1
      IF(S1(1) .LE.99)  JJ=JJ*2
      IF(I.LT.J ) CH=CH+     ( IL(I)+IL(J)-
     *IABS(IL(I)-IL(J)))/2 +1
      IF(I.EQ.J.AND.IQ(I).NE.1) CH=CH+IL(I)
      IF(KSU1.NE.1.AND.(IQ(I).EQ.1.   !KSU changed to KSU1, Igor 10/4/2011
     * AND.IQ(J).EQ.JJ-1).OR.
     * (IQ(J).EQ.1.AND.IQ(I).EQ.JJA-1))
     * CH=CH+ MIN0(IL(I),IL(J))
 1111 CONTINUE
c
      ione = 1
c
      CALL COEF(GAM,DIR,IS,ione, KSU1,IN,IL,IQ,IP,
     *   L1,S1,CH,KK,ISIGM1,ISIGM2,IG,1 )
      NU=CH
111   write(66, 4) IS,NU,IG
4     FORMAT(' wave functions of excited state',
     * /' in Hartree-fock approximation',
     * /,
     * 'the number of sells',17X,
     *I5/,'the number of coeffitiants',12X, I5,
     * /' the number of exchange coeffitiants ',I5)
      READ(2 )(A(JJ1),JJ1=1,15)
      IF(A(1).EQ.1.) GOTO 1980
      write(66, 1981) A(1)
1981  FORMAT(5X,' Read  Bf from tape  ',
     *' not graund state ',E15.8)
      RETURN
1980  PRE=1.D-5
      IF(ABS(A(7)-Z).LT.PRE.AND.
     *ABS(DFLOAT(IS-1) -A(9)).LT.1.+PRE)  GOTO 1982
      write(66, 1983) (A(JJ1),JJ1=1,15),
     *R, Z,IS,PRE
1983     FORMAT(5X,'Some parameters',
     */,5X,'PL=',E15.8,3X,'R=',E15.8,
     * 3X,'N=',E15.8,3X,'H=',E15.8,3X,
     *'BET=',E15.8,/3X,'EPS=',E15.8,3X,
     * 'Z=',E15.8,3X,/,2(4(3X,E15.8),/),
     *3X, 2E12.5,I3,E15.8)
       RETURN
 1982   CONTINUE
      ALFA1=A(10)
      RN=A(2)
      RO1N=A(11)
      BET1=A(5)
      NE1=IFIX(sngl(A(3)))
      H1=A(4)
      IF(ABS(A(2)-R).LT.PRE.OR.K2.LT.1) GOTO 440
      ROH=ALFA1*R+BET1*LOG(R)
      NH=1+IFIX(sngl((ROH-A(11))/H1))
      ROH=A(11)+H1*DFLOAT(NH)
      IF(NE1.LT.NH) GOTO 7055
      write(66, 7054) RN,R
7054  FORMAT(' ',88('*')/' Error in ',
     * 'Calculation of the number of integration points '
     * /' R CTAPOE=',E15.8,' R HOBOE=',E15.8)
7055  CONTINUE
      CALL KOMON( Z,R,BET1,H1,GAM,EPS,ALFA,RO1N,
     * RO2,NH,IN,IL,IQ,NU,KK,ISIGM1,ISIGM2,IS,L1)
8921  FORMAT(' R first',22X,E15.8  )
      IF(NE5.GE.NH+5.AND.NE5IS.GE.(NH+5)*IS)
     * GOTO 7056
      write(66, 7057) NE5,NH,NE5IS
7057  FORMAT(' ',88('*'),/' the length of work arrayes ',
     * '=',I4,' is less then it is need.',
     * 'Set it as not less then',I4/' as well as',
     * 'Check the length of array R8, which is now  =',I6)
      RETURN
7056  CONTINUE
      IF(ABS(ALFA1-ALFA).LT.0.1.AND.ABS(RO2-ROH)
     * .LT.PRE) GOTO 7058
      write(66, 7059) RO1N,H1,BET1,RN,R,RO2,ROH,
     *ALFA1,ALFA
7059  FORMAT(' Error in expanding integral interval: ',
     * /' RO1=',E12.5,
     * ' H1=',E12.5,' BET=',E12.5,' R CTAPOE=',
     * E12.5,' R HOBOE=',E12.5/' RO final ',
     * 'Calculated in subroutine KOMON=',E15.8/' RO,',
     * 'Calculation in subr. FCHF=',E15.8/
     * ' ALFA old=',E15.8,' ALFA new=',E15.8)
      RETURN
7058  CONTINUE
      CALL VAR1(Z,R,BET1,H1,ALFA1,RO1N,
     * ROH,R1,NH,NH+5)
      write(66, 8921) R1(3)
      IRA=IS-1
      DO 55 IW=1,IRA
      IA=IW+IZON1-1
      IZCH=2
      CALL  TAPE(R3,XNOBL,XPHS,IZCH,NE1,KOLE,2)
      MMOM=IFIX(sngl(R3(2)))
      NOBL=IFIX(sngl(XNOBL))
      R3(2)=R3(NE1+3)
      R3(NE1+3)=0.
      DO 60006 J=1,NE1
60006 R3(J+2)=R3(J+2)*SQRT(R3(2))
      write(66, 7007) R3(1),R3(2),NOBL,MMOM,KOLE
 7007 FORMAT(/' Energy=',E15.8/,' HOPMA=',
     * E15.8,/,  ' Principal quantum number ',
     * I5,/,' Orbital moment',6X,I5 /
     *,' The number of electrons', 3X,I5)
      I101=NE1+1
      DO 402 LL=I101,NH
402   R3(LL+2)=0.0
      NH5=NH+5
      MMM=(IW-1)*NH5
      DO 1984 I100=1,NH5
1984  R8(MMM+I100)=R3(I100)
55    CONTINUE
      REWIND 2
      GOTO 44040
440   CONTINUE
      CALL KOMON( Z,RN,BET1,H1,GAM,EPS,ALFA1,
     * RO1N,RO2N,NE1,IN,IL,IQ,
     * NU,KK,ISIGM1,ISIGM2,IS,L1)
      CALL VAR1(Z,RN,BET1,H1,ALFA1,RO1N,
     * RO2N ,R1,NE1,NE1+5)
      write(66, 8921) R1(3)
      IRA=IS-1
      DO 596 IW=1,IRA
      IA=IW+IZON1-1
      IZCH=2
      CALL  TAPE(R5,XNOBL,XPHS,IZCH,NE1,KOLE,2)
      MMOM=IFIX(sngl(R5(2)))
      NOBL=IFIX(sngl(XNOBL))
      R5(2)=R5(NE1+3)
      R5(NE1+3)=0.
      DO 6 J=1,NE1
00006 R5(J+2)=R5(J+2)*SQRT(R5(2))
      write(66,    7) R5(1),R5(2),NOBL,MMOM,KOLE
00007 FORMAT(/' Energy=',E15.8/,' HOPMA=',
     * E15.8,/,  ' Principal quantum number ',
     * I5,/,' Orbital moment',6X,I5 /
     *,' The number of electrons', 3X,I5)
      NE5=NE1+5
      MMM=(IW-1)*NE5
      DO 1989 I100=1,NE5
1989  R8(MMM+I100)=R5(I100)
596   CONTINUE
      REWIND 2
      NH=NE1
44040 ISIGM=IS
      I=IS
      QQ3=DFLOAT(IS)
      QQ4=DFLOAT(IN(IS))
      QQ5=DFLOAT(K2)
      QQ6=DFLOAT(K1)
      SOST =2.D 00
      XN1=Z
      WRITE(1 )SOST,RN,A(3),H1,BET1,EPS,Z,XN1,
     *QQ3,QQ4,A(11),QQ5,QQ6,AKAP,DKAP
      IF(K2.LT.1) GOTO 20
      RO1=-BET*(10.D 00 + LOG(Z))
      RO2=RO1+DFLOAT(NE)*H
      ALFA=(RO2- BET*LOG(R))/R
      write(66, 6423)
6423  FORMAT(/)
      DO 8 J10=1,K2
      write(66, 9)   J10
   9  FORMAT('  Wave function of discret degenerate state',
     *'in Hartree-Fock approximation',I5)
      CALL ZERO( R5,Z,R8,TETA,NH ,
     * IS,IN,IQ,NH+5,NE5IS,KSU1)
      ISIGM=IS
       DO 9753 IT=1,IS
 9753  ALM3(IT)=0.0
      CALL MAIN1(ALFA1,BET1,H1,Z,EPS,GAM,TETA,
     *R1,R2,R3,R4,R5,R6,R8,C3,ALA,ALA1,ALA2,
     *ALM3,MU,IS,NH,NH+5,(NH+5)*IS, IQ,ISIGM,IN,IL,
     *IG,NU,ISIGM1,ISIGM2,KK,UD,KSU1,alpha,r0)
      IF(ABS(A(2)-R).LT.PRE) GOTO 44041
      C=0.
      R5(2)=DFLOAT(IL(IS))
      R5(NE1+3)=0.0
      D=0.
      DO 6658 J=1,NE1
      P=D
      D=C
      C=(R5(J+2) *R1(J+2) /(ALFA1*R1(J+2)+BET1))**2
      A1=5.*C
      B=8.*D
      A1=(A1+B-P)*H/12.
 6658 R5(NE1+3)=R5(NE1+3)+A1
      DO 429 J=1,NE1
429   R5(J+2)=R5(J+2)/SQRT(R5(NE1+3))
44041 CONTINUE
      IA=IZON
      IZCH=7
      XIS=DFLOAT(IN(IS))
      CALL  TAPE(R5,XIS,XPHS,IZCH,NE1,IQ(IS),1)
      CALL INTZZZ( R5,R8,RO1N,BET1,ALFA1,R7,
     *H1,R2,IS,NH,NH+5,NE5IS,KT,IN,KSU1,r1)
      IN(IS)=IN(IS)+1
    8 CONTINUE
      IN(IS)=IN(IS)-K2
      DO 6070 II=2,IS
      I101=(II-1)*(NE1+5)
      I102=(II-1)*(NH+5)
      I103=NE1+5
      DO 6071 J=1,I103
6071  R8(I101+J)=R8(I102+J)
6070  CONTINUE
   20 CONTINUE
      IF(K1.LT.1) GOTO 21
      NE5=NE1+5
      DO 11 II=1,K1
      write(66,   12) II
   12 FORMAT(    '  Continuum wave function', I5)
      CALL ZERO( R5,Z,R8 ,TETA,
     * NE1,IS,IN,IQ,NE5,NE5*IS,KSU1)
      ISIGM=IS
       DO 8642 IT=1,IS
 8642  ALM3(IT)=0.0
      CALL MAIN2(ALFA1,BET1,H1,Z,AKAP,EPS,
     *TETA,R1,R2,R3,R4,R5,R6,R8 ,ALA,ALA1,
     *ALA2,ALM3,C3,GAM,MU,IS,NE1,NE5,NE5*IS,
     *IL,KK,ISIGM1,ISIGM2,IN,IQ,NU,IG,ISIGM,UD,ksu1,alpha,r0)
      IA=IZON
      IZCH=7
      CALL TAPE(R5,AKAP,R5(NE1+4),
     * IZCH,NE1,IQ(IS),1)
      AKAP=AKAP+DKAP
      CALL INTZZZ( R5,R8 ,RO1N,BET1,ALFA1,
     * R7,H1,R2,IS,NE1,NE5,NE5*IS,KT,IN,KSU1,r1)
11       CONTINUE
21       CONTINUE
       RETURN
       END

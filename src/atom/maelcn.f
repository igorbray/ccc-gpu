      SUBROUTINE MAELCN(KAPO,DKAP,EE,U,V,
     * PHS,ML,MV,CNST, EHL, PF,PB,SC,S,
     * SUMO,SUML,SUMM,SUMV,R,Q,AA,ZA,SCH,
     * SPT,IP,           NCH,NPHS,
     * ND,NES,LP,HL,LH,NNL,NNNPHS,PQN)
c !!!
c     JJBEG)
      real *8 beta
      REAL*8 ZA,PF(  30),PB(  30),PHS(NPHS),
     * SC(30),S(NNNPHS,NNNPHS),SUMO(30),
     * SUML(30),SUMM(30),SUMV(30),
     * R(NNNPHS,NNNPHS),CNST,X,A,B,C,D,G,
     * GP,OM,OM1,CN,Y,Z,  RL,RV,RL1,
     * DKAP(NCH),EE(NPHS),EHL(NCH),
     * IV1,IL1,RV1,BETL,BETV,Q(NNNPHS,NNNPHS),AA(1)
      REAL*8 IL,IV,V(NPHS,NPHS),U(NPHS,NPHS),
     * AX,BX,RPL(200),RPV(200),BTL(200),BTV(200),
     * KAPO(NCH),ML(NPHS),MV(NPHS)
      INTEGER BOO,CONT,SCH(  3),SPT(  30),
     * ND(NCH),NES(NCH), IP(1),PQN(1),
     * LP(NCH),LH(NCH),HL(NCH),NNL(NCH)
      character  CI(5),  str1,  str2,  cie
       data      CI/'S','P','D','F','G'/,
     * STR1/'-'/,STR2/'>'/,CIE/'E'/
       data  btl/200*0./,  btv/200*0./
      J2=0
      DO 1 I=1,NCH
      J1=J2+ND(I)
      II=J2+1
      DO 2 J=II,J1
2     SC(J)=1.
      J2=J2+NES(I)
      IF(J1.EQ.J2) GOTO 55
      X=KAPO(I)
      A=X
      B=DKAP(I)
      C=A/B
      D=B/3.
      A=A*A/12.
      SC(J1+1)=2.*X*D
      X=X+B
      SC(J1+2)=8.*X*D
      II=J1+3
      LL=J2-2
      DO 3 J=II,LL,2
      X=X+B
      SC(J)=4.*X*D
      X=X+B
3     SC(J+1)=8.*X*D
      X=X+B
      SC(J2)=2.*X*D
55    CONTINUE
1     CONTINUE
      WRITE(1 ,10)
10    FORMAT(3X,'final coeff. ')
      WRITE(1 ,11)(SC(I),I=1,NPHS)
11    FORMAT(5X,E13.6)
      G=2./3.D 00
      N=2*NPHS
      JJ2=0
      GP=G*3.1415927D 00
      DO 12 I=1,NCH
12    SCH(I)=0
      DO 13 J=1,NPHS
13    SPT(J)=0
      DO 14 I=1,NCH
      SUMO(I)=0.0
      SUML(I)=0.0
      SUMV(I)=0.0
14    SUMM(I)=0.0
      DO 15 II=1,NCH
      JJ1=JJ2+1
      JJ2=JJ2+NES(II)
      LL=II+1
       IF(LL.GT.NCH)  GOTO 901
      DO 19 I=LL,NCH
      IF(HL(I).EQ.HL(II).AND.SCH(II).EQ.0) GOTO 18
      GOTO 19
18    SCH(I)=-II
      SCH(II)=I
19    CONTINUE
901     CONTINUE
      DO 40 JJ=JJ1,JJ2
      CONT=0
      IF(JJ.GE.JJ1+ND(II)) CONT=1
      OM=EE(JJ)
      OM1=OM
      CN=1.
      WRITE(6, 30) JJ,II
30    FORMAT(2X,20('*'),/,' point  N.',I2,
     * ' (channal ',I2,')')
      IF(CONT.EQ.0) GOTO 90
      I1=SCH(II)
      IF(I1.LT.0.AND.SPT(JJ).EQ.0) GOTO 90
      IF(I1.LT.0.AND.SPT(JJ).NE.0) GOTO 40
      IF(I1.EQ.0) GOTO 90
      J2=JJ2
      LL1=II+1
      LL2=I1-1
      IF(LL1.GT.LL2)   GOTO 902
      DO 31 I=LL1,LL2
31    J2=J2+NES(I)
902   J1=J2+ND(I1)+1
      J2=J2+NES(I1)
      DO 32 J=J1,J2
      IF(ABS(EE(J)/OM-1.).GE.1.D-05) GOTO 32
      SPT(JJ)=J
      SPT(J)=JJ
      GOTO 90
32    CONTINUE
90    CONTINUE
      IF(SPT(JJ).NE.0) WRITE(6, 33) SPT(JJ),I1
33    FORMAT('+',22X,' - point N.',I2,
     * ' (channal ',I2,')')
      WRITE(6, 3358) OM
 3358 FORMAT(' energy of the photon=',E16.6)
      CALL CQD2(OM,EE,G,SC,OM1,EHL,
     * KAPO,GP,U, V,DKAP,S,R,PF,PB,
     * Q,AA,IP,     NCH,
     * NPHS, NES, N,ND, CONT,
     *     BOO,ierror,NPHS*2)
      if ( ierror .eq. 1) goto 200
      IF(CONT.EQ.1) GOTO 35
      X=0.0
      Y=X
      Z=X
      DO 36 I=1,NPHS
      A=U(I,JJ)
      B=V(I,JJ)
      C=PF(I)
      D=PB(I)
      X=X+R(JJ,I)*A+R(JJ,I+NPHS)*B
      Y=Y+A*A*C+B*B*D
   36 Z=Z-(A*A*C*C-B*B*D*D)/ SC(I)
      CN=1.-Z
      A=G*X/CN
      OM1=OM+A
      C=Z/G
      WRITE(6, 37) OM1,X,Y,C
   37 FORMAT(' shifted photon energy=',
     *F8.5,5X,'correction to the energy=',E10.3,2X,
     *'A=',E10.3,
     * 2X,'B=',E10.3)
      IF(ABS(A/OM).LE.1.D-4) GOTO 38
      CALL CQD2(OM,EE,G,SC,OM1,EHL,
     * KAPO,GP,U,V,DKAP,S,R,PF,PB,
     *Q,AA,IP,         NCH,
     * NPHS,NES,N,ND, CONT,
     * BOO,ierror,NPHS*2)
      if ( ierror .eq. 1) goto 200
   38 OM=OM1
   35 CONTINUE
      IPR=1
      CALL FCROSS(R,ML,MV,S,OM,CN,SUML,
     * SUMV,SUMM,SUMO,SC,CNST,RL,RV,IL,
     * IV,AX,BX,II,JJ,NCH,NPHS,BOO,NNL,LH,
     * CONT,   NPHS*2,IPR)
      RPL(JJ)=AX
      RPV(JJ)=BX
      IF(SPT(JJ).EQ.0) GOTO 40
      RL1=RL
      IV1=IV
      IL1=IL
      RV1=RV
      K=SCH(II)
      K1=SPT(JJ)
      IPR=1
      CALL FCROSS(R,ML ,MV,S,OM,CN,SUML,
     * SUMV,SUMM,SUMO,SC ,CNST,RL,RV,IL,
     *  IV,AX,BX,K,K1,NCH,NPHS,BOO,NNL,LH,
     * CONT, NPHS*2,IPR)
      RPL(K1)=AX
      RPV(K1)=BX
      J1=LP(II)
      J2=LP(K)
      I1=LH(II)
      IF((J1-I1)*(I1-J2).EQ.1) GOTO 42
      GOTO 40
   42 CONTINUE
      IF(J1.GT.J2) GOTO 46
      J=I1+2
      J1=I1-1
      GOTO 47
   46 J=I1-1
      J1=I1+2
   47 C=PHS(JJ)-PHS(K1)
      A=COS(C)
      B=SIN(C)
      BETL=BETA(RL,IL,RL1,IL1,A,B, J,J1,I1)
      BETV=BETA(RV,IV,RV1,IV1,A,B, J,J1,I1)
      WRITE(6, 48) BETL,BETV
 48   FORMAT(' parameter of angular distribution'
     * ,2X,'BETA L=',E13.6,5X,'BETA V=',E13.6)
      BTL(JJ)=BETL
      BTV(JJ)=BETV
      BTL(K1)=BETL
      BTV(K1)=BETV
   40 CONTINUE
      WRITE(6, 49) SUMO(II),SUML(II),SUMV(II),SUMM(II)
49    FORMAT(2X,18('*'),/,5X,'SUMO=',E13.6,/,
     *5X,'SUML=',E13.6,/,5X,'SUMV=',E13.6,/,
     *5X,'SUMM=',E13.6)
      WRITE(6, 1122) ZA
 1122 FORMAT(//10X,'Z=',F4.0,5X,'transition')
      L17=LH(II)+1
      L18=LP(II)+1
      WRITE(6, 1123) PQN(II),CI(L17),STR1,
     * STR2,CIE,CI(L18),EHL(II)
1123  FORMAT('+',30X,I1,A1,2(1X,2A1),7X,
     * 'E=',F9.5)
      IF(NCH.EQ.1) GOTO 1134
      WRITE(6, 1125)
1125  FORMAT(/12X,'interaction with the channals')
      K9=NCH
      K8=II+1
      K6=1
7751  CONTINUE
      IF(K9.LT.K8) GOTO 1124
      DO 1127 K7=K8,K9
      L17=LH(K7)+1
      L18=LP(K7)+1
      IF(K6.NE.1) GOTO 2201
      WRITE(6, 1126) PQN(K7),CI(L17),STR1,
     * STR2,CIE,CI(L18)
1126  FORMAT('+',40X,I1,A1,2(1X,2A1))
      K6=K6+1
      GOTO 1127
2201  CONTINUE
      IF(K6.NE.2) GOTO 2202
      WRITE(6, 1136) PQN(K7),CI(L17),STR1,
     * STR2,CIE,CI(L18)
1136  FORMAT('+',52X,I1,A1,2(1X,2A1))
      K6=K6+1
      GOTO 1127
2202  CONTINUE
      IF(K6.NE.3) GOTO 2203
      WRITE(6, 1146) PQN(K7),CI(L17),STR1,
     * STR2,CIE,CI(L18)
1146  FORMAT('+',64X,I1,A1,2(1X,2A1))
      K6=K6+1
      GOTO 1127
2203  CONTINUE
      IF(K6.NE.4) GOTO 2204
      WRITE(6, 1156) PQN(K7),CI(L17),STR1,
     * STR2,CIE,CI(L18)
1156  FORMAT('+',76X,I1,A1,2(1X,2A1))
      K6=K6+1
2204  CONTINUE
1127  CONTINUE
1124  CONTINUE
      IF(II.EQ.1.OR.K6.EQ.NCH) GOTO 1134
      K8=1
      K9=II-1
      GOTO 7751
1134  CONTINUE
      WRITE(6, 3002)
3002  FORMAT('  Energy            cross section ',
     * '   photo-ionization',17X,'      face coefficiant',
     * '          '/
     *       '   photon       HF-dl        HF',
     * '-CK',11X,'pCfO-dl      pCfO-CK      shift',6X,
     * 'angular anizOTPOpii'/85X,'BETA L   BETA V')
      DO 3323 LY=JJ1,JJ2
      AX=DFLOAT(NNL(II))/DFLOAT(2*LH(II)+1)
      BX=AX*EE(LY)*ML(LY)*ML(LY)/3.
      AX=AX*4.*MV(LY)*MV(LY)/(3.*EE(LY))
      IF(LY.LE.JJ1+ND(II)-1) GOTO 3325
      BX=BX*CNST
      AX=AX*CNST
3325  CONTINUE
      WRITE(6, 3003) EE(LY),BX,AX,
     * RPL(LY),RPV(LY),PHS(LY),BTL(LY),BTV(LY)
3323  CONTINUE
3003  FORMAT(2X,F7.4,2(3X,2(E11.4,2X)),
     * 1X,F8.4,10X,2(F6.3,3X))
      WRITE(6, 3340)
3340  FORMAT(/)
   15 CONTINUE
  200 CONTINUE
      RETURN
      END


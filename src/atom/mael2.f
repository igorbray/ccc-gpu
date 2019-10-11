      SUBROUTINE MAEL2(H,ALPHA,BETA,R,FUNC,
     * U,V,COEF,Y,AU,L,CHA,CHB,LH,LP,HL,NES,
     * NNL,NN,NWF,NPHS,NHL, NCH,SPI)
      INTEGER CHA,CHB,EQOC,HLA,HLB,HL(NCH),
     * LH(NCH),LP(NCH),NES(NCH), NNL(NCH),SPI(NCH)
      REAL*8 R(NN),FUNC(NWF,700),CA,D,H,CB,ALPHA,
     * V(NPHS,NPHS),A,U(NPHS,NPHS),COEF(1),Y(1,1),
     * A1,A2,A3,B,B1,C,F,BETA,AU(1),  ddff
       real *8 rk
      ISPI=0
      IF(NNL(CHA).EQ.2*LH(CHA)+1.AND.
     * SPI(CHA)*SPI(CHB).GT.0) ISPI=1
      IF(NNL(CHA).EQ.2*LH(CHA)+1.AND.
     * SPI(CHA)*SPI(CHB).LT.0) ISPI=2
      IF(ISPI.NE.0) WRITE(6, 90845)
90845 FORMAT('  Calculation with WF calculated by spin-',
     * '-polalization method')
      INTRA=0
      IF(CHA.EQ.CHB) INTRA=1
      LHA=LH(CHA)
      LPA=LP(CHA)
      HLA=HL(CHA)
      NA=NES(CHA)
      MA=0
      J=CHA-1
      IF(J.LT.1) GOTO 293
      DO 2 I=1,J
2     MA=MA+NES(I)
293   CA=NNL(CHA)/DFLOAT(2 *LHA+1 )
      IF(INTRA.EQ.0) GOTO 3
      LHB=LHA
      LPB=LPA
      HLB=HLA
      NB=NA
      MB=MA
      I1=MA+1
      J1=I1
      I2=MA+NA
      J2=I2
      CALL MEZERO(U,I1,I2,J1,J2,NPHS)
c !!!
c      CALL MEZERO(U,I1,I2,J1,J2,NA,NB,NPHS)
      D=DFLOAT((2 *LHA+1 )*(2 *LPA+1))
      CB=CA
      EQOC=1
      GOTO 10
3     LHB=LH(CHB)
      LPB=LP(CHB)
      HLB=HL(CHB)
      NB=NES(CHB)
      MB=0
      J=CHB-1
      WRITE(1 ,1) NA,NB
1     FORMAT('  coulomb matrix elements ',
     * '      ',I3,'*',I3 )
      IF(J.LT.1) GOTO 275
      DO 4 I=1,J
4     MB=MB+NES(I)
275   CB=NNL(CHB)/DFLOAT (2 *LHB+1 )
      EQOC=0
      IF(ABS(CA-CB).LT.0.1D-5) EQOC=1
      ddff=DFLOAT((2 *LHA+1 )*(2 *LPA+1)*
     * (2 *LHB+1 )*(2 *LPB+1))
      D=SQRT(ddff)
      I=IABS(LHA-LHB)
      J=IABS(LPA-LPB)
      LMN=J
      IF(I.GT.J) LMN=I
      I=LHA+LHB
      J=LPA+LPB
      LMX=J
      IF(I.LT.J) LMX=I
      N=(LMX-LMN)/2+1
      A=DFLOAT(  IPARI(LMN+L)*(2 *L+1 )   )
      K1=0
      DO 101 K=LMN,LMX,2
      K1=K1+1
      CALL IOT3(  LHB,  K,  LHA,0,0,0,A1)
      CALL IOT3(  LPA,  K,  LPB,0,0,0,A2)
      CALL IOT6(LPA,L,LHA,LHB,K,LPB,A3)
      COEF(K1)=A*A1*A2*A3
      CALL YK(Y,H,ALPHA,BETA,R,FUNC,AU,
     * K,HLA,HLB,  K1+1, NN,       NWF)
101   CONTINUE
      CALL IOT3(  LPA,  L,  LHA,0,0,0,B)
      CALL IOT3(  LPB,  L,  LHB,0,0,0,B1)
      C=B*B1
      IF(INTRA.EQ.1) C=B*B
      CALL WRICO(COEF,N)
      I1=MA+1
      I2=MA+NA
      DO 102 I=I1,I2
      CALL YK(Y,H,ALPHA,BETA,R,FUNC,AU,
     * L,HLA,NHL+I,  1, NN,       NWF)
      J1=MB+1
      J2=MB+NB
      IF(INTRA.EQ.1) J2=I
      DO 102 J=J1,J2
      A=C*RK(Y,H,ALPHA,BETA, R,FUNC,AU,
     * HLB,NHL+J,  1, NN,       NWF)
      B=0
      DO 104 K1=1,N
104   B=B-COEF(K1)*RK(Y,H,ALPHA,BETA,R,
     * FUNC,AU, NHL+I,NHL+J,  K1+1,NN,NWF)
      IF(ISPI.NE.2) GOTO 6421
      F=D*CA*A
      U(I,J)=F
      U(J,I)=D*CB*A
      IF(EQOC.EQ.1) U(J,I)=F
      GOTO 6423
6421  CONTINUE
      F=D*(CA*A+B)
      U(I,J)=F
      U(J,I)=D*(CB*A+B)
      IF(EQOC.EQ.1 )U(J,I)=F
6423  CONTINUE
102   CONTINUE
      CALL WRITC(U,MA,NA,MB,NB,INTRA,NPHS)
10    CONTINUE
      WRITE(1 ,1)NA,NB
      IF(INTRA.EQ.1.AND.NNL(CHA).EQ.1) GOTO 12
      GOTO 13
c !!!
c 12    CALL MEZERO(V,MA+1,MA+NA,MB+1,J2,NA,NB,NPHS)
12    CALL MEZERO(V,MA+1,MA+NA,MB+1,J2,NPHS)
      GOTO 50
13    CONTINUE
      I=IABS(LPA-LHB)
      J=IABS(LHA-LPB)
      LMN=J
      IF(I.GT.J) LMN=I
      I=LPA+LHB
      J=LHA+LPB
      LMX=J
      IF(I.LT.J) LMX=I
      N=(LMX-LMN)/2+1
      A=DFLOAT( IPARI(LMN+L)*(2 *L+1))
      K1=0
      DO 201 K=LMN,LMX,2
      K1=K1+1
      CALL IOT3(  LPA,  K,  LHB,0,0,0,A1)
      CALL IOT3(  LHA,  K,  LPB,0,0,0,A2)
      CALL IOT6(LPA,L,LHA,LPB,K,LHB,A3)
      COEF(K1)=A*A1*A2*A3
201   CONTINUE
      CALL IOT3(  LPA,  L,  LHA,0,0,0,B)
      CALL IOT3(  LHB,  L,  LPB,0,0,0,B1)
      C=B*B1
      CALL WRICO(COEF,N)
      I1=MA+1
      I2=MA+NA
      DO 202 I=I1,I2
      CALL YK(Y,H,ALPHA,BETA,R,FUNC,AU,
     * L,HLA,NHL+I,  1, NN,       NWF)
      K1=0
      DO 203 K=LMN,LMX,2
      K1=K1+1
      CALL YK(Y,H,ALPHA,BETA,R,FUNC,AU,
     * K,HLB,NHL+I,  K1+1, NN       ,NWF)
203   CONTINUE
      J1=MB+1
      J2=MB+NB
      IF(INTRA.EQ.1) J2=I
      DO 202 J=J1,J2
      A=C*RK(Y,H,ALPHA,BETA,R,FUNC,AU,
     * HLB,NHL+J,  1, NN,       NWF)
      B=0.
      DO 204 K1=1,N
204   B=B-COEF(K1)*RK(Y,H,ALPHA,BETA,
     * R,FUNC,AU,HLA,NHL+J,  K1+1,NN,NWF)
      IF(ISPI.NE.2) GOTO 8643
      F=D*CA*A
      V(I,J)=F
      V(J,I)=D*CB*A
      IF(EQOC.EQ.1) V(J,I)=F
      GOTO 8645
8643  CONTINUE
      F=D*(CA*A+B)
      V(I,J)=F
      V(J,I)=D*(CB*A+B)
      IF(EQOC.EQ.1) V(J,I)=F
8645  CONTINUE
202   CONTINUE
      CALL WRITC(V,MA,NA,MB,NB,INTRA,NPHS)
50    CONTINUE
      WRITE(1 ,8)
8     FORMAT(2X,115('-'))
      RETURN
      END

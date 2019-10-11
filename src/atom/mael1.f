      SUBROUTINE MAEL1(EE,CNST,BETA,
     * ALPHA,FUNC,H,R,MV,ML,YD,AU,
     *   NNL,NCH,NES,HL,ND,LH,LP,  NN,
     *   NHL,      NWF,NPHS)
      INTEGER HLI,SI,NES(NCH),HL(NCH),ND(NCH),
     * LH(NCH),LP(NCH),NNL(NCH)
      REAL*8FUNC(NWF,700),R(NN),YD(700),
     * EE(NPHS),F5,D,A,B,C,F,OM,CN,F4 ,
     * AU(700),ML(NPHS),MV(NPHS),  ddff,
     * CNST,BETA,ALPHA,H,  BET2,F1,F2,F3
       real * 8 simps1
      K1=0
      BET2=BETA/2.
      DO 118 I=1,NCH
      WRITE(6, 102) I
102       FORMAT('    channal       N   ',I2,/
     * 9X,'energy ','   amplitude ',7X,
     * 'Cross section'/10X,'photon',4X,2('L',7X,'V',7X))
      NESI=NES(I)
      HLI=HL(I)
      NDI=ND(I)
      L1=LH(I)
      IF(LP(I).GT.LH(I)) L1=LP(I)
      ddff=DFLOAT(L1)
      CN=DFLOAT(   IPARI(L1))*SQRT(ddff)
      SI=-1
      IF(L1-LH(I).GT.0) SI=1
      F1=FUNC(HLI,1)
      F2=FUNC(HLI,2)
      F3=FUNC(HLI,3)
      F4=FUNC(HLI,4)
      F5=FUNC(HLI,5)
      D=1./(12.*H)
      DO 103 J=1,NN
      A=R(J)
      B=ALPHA*A+BETA
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
103   YD(J)=A/B*(-C+F*(DFLOAT(  SI*L1)-BET2/B)/B)
      DO 118 J1=1,NESI
      K1=K1+1
      N1=NHL+K1
      OM=EE(K1)
      DO 114 J=1,NN
      A=R(J)
      B=A/(ALPHA*A+BETA)
      AU(J)=FUNC(HLI,J)*FUNC(N1,J)*A*B*B
114    CONTINUE
      A=CN*SIMPS1(AU,H,NN)
      ML(K1)=A
      DO 115 J=1,NN
115   AU(J)=FUNC(N1,J)*YD(J)
      B=CN*SIMPS1(AU,H,NN)
      MV(K1)=B
      C=DFLOAT(    NNL(I))/DFLOAT(2 *LH(I)+1 )
      A=C*OM*A*A/3.
      B=4.*C*B*B/(3.*OM)
      IF(J1.LE.NDI) GOTO 116
      A=CNST*A
      B=CNST*B
116   CONTINUE
      WRITE(6, 117) K1,OM,ML(K1),
     * MV(K1),A,B
117   FORMAT(2X,I4,2X,5(F7.3,1X),
     * I4,1X, E12.5,1X,I2,1X,E15.8)
118   CONTINUE
      RETURN
      END
      SUBROUTINE WRICO(COEF,N)
      REAL*8    COEF(N)
      WRITE(1 ,1)
    1 FORMAT(/,'  coefficients:')
      WRITE(1 ,3)(COEF(I),I=1,N)
    3 FORMAT(5X,F8.5  )
      RETURN
      END

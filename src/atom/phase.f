      SUBROUTINE PHASE(X,P,T,R1,R2,L,NE,NE5,IM1)
C   OèPEÑEãEHàE  îAáOBOÉO  CÑBàÉA
      REAL*8X,P,T,R1(NE5 ),R2(NE5 ),C,B,D,F,H,G
      DO 3 IW=1,NE
      IF(R2(IW+2).LT.0.0005.AND.
     *R2(IW+2).GT.-0.05) GOTO 10
      GOTO 50
   10 R2(IW+2)=0.0
    3 CONTINUE
   50 C=R2(3)
      KSI=0
      IRA=NE-1
      DO 1 IW=1,IRA
      B=C
      C=R2(IW+3)
      IF(B*C.GE.0.0)  GOTO 1
      KSI=KSI+1
      D=R1(IW+2)
      F=R1(IW+3)
      H=R2(IW+2)
      G=R2(IW+3)
    1 CONTINUE
      IM1=KSI
      IF(KSI.EQ.0) RETURN
      B=(G*D-H*F)/(G-H)
      P=DFLOAT(KSI)*3.14159-X*B-T/(X*2.)*
     *LOG(2.*B*X)+3.14159*FLOAT(L)/2.
      IM1=KSI
      RETURN
      END

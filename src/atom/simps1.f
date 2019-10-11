      REAL  FUNCTION SIMPS1(F,H,NN)
      REAL*8    F(NN)    ,H,S,S1
      S=0.0
      S1=S
      DO 1 I=3,NN,2
      S=S+F(I-1)
1     S1=S1+F(I)
      S=F(1)+4.*S+2.*S1-F(NN)
      IF(IPARI(NN).GT.  0) S=S+(9.*F(NN)-
     *F(NN-2))/4.+F(NN-1)
      SIMPS1=S*H/3.
      RETURN
      END

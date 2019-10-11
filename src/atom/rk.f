      REAL    FUNCTION RK(Y,H,ALPHA,
     * BET,R,FUNC,AU,I,J,  K,NN, NWF)
      REAL*8    Y(20 ,700),FUNC(NWF,700),
     * H,ALPHA,BET, A,B, AU(700),R(NN)
      real *8 simps1
      IF(K.GT.21)  WRITE(6, 55) K
55    FORMAT(' ',88('*')/ '  K=',I5,
     * '  situation is not allowed K>20'/
     * '   increase the length of array Y ',
     * ' HA Y(40,700)  BO  BCEX subr. ',
     * ' where this array is used')
      DO 10 N=1,NN
      A=R(N)
      B=ALPHA*A+BET
10    AU(N)=FUNC(I,N)*FUNC(J,N)*Y(K,N)*A/(B*B)
      RK=SIMPS1(AU,H,NN)
      RETURN
      END

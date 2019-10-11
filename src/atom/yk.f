      SUBROUTINE YK(Y,H,ALPHA,BET,
     * R,FUNC,  AU, K,I,J,  M, NN, NWF)
      REAL*8    Y(20 ,700),R(NN),FUNC(NWF,700),
     * AU(700) ,H,ALPHA,BET  ,H12,Z,Y1,Y2,A,B
      H12=H/12.
      Z=0.0
      Y1=Z
      Y2=Z
      DO 10 N=1,NN
      A=R(N)
      B=ALPHA*A+BET
      A=A/B
      A=A*A*FUNC(I,N)*FUNC(J,N)
      B=DFLOAT(   K)/B
      Z=(Z+H12*(5.*A+8.*Y1-Y2))/(1.+5.*H12*B)
      Y2=Y1
      Y1=A-B*Z
10    AU(N)=Z
      Y(M,NN)=Z
      Y1=0.0
      Y2=Y1
      NN1=NN-1
      DO 20 L=1,NN1
      N=NN1+1-L
      A=R(N)
      B=ALPHA*A+BET
      A=-DFLOAT (2 *K+1 )*AU(N)/B
      B=DFLOAT(K+1 )/B
      Z=(Z-H12*(5.*A+8.*Y1-Y2))/(1.+5.*H12*B)
      Y2=Y1
      Y1=A+B*Z
20    Y(M,N)=Z
      RETURN
      END

      REAL FUNCTION F1(K,L,Q)
      INTEGER Q
      REAL*8A
      CALL IOT3(L,K,L,0,0,0     ,A)
      F1=DFLOAT(Q*(Q-1))/(2.*DFLOAT( 4 *L+1))*
     *A*A*DFLOAT(2 *L+1 )
      RETURN
      END

      REAL FUNCTION F4(K,L,Q)
      REAL*8 A
      INTEGER Q
      CALL IOT3(L,K,L,0,0,0,A)
      F4=DFLOAT(Q*(Q-1)*(2*L+1))/DFLOAT(4*L)*A*A
      RETURN
      END

      REAL FUNCTION F3(K,L,L2,L3)
      REAL*8 A,B,C
      CALL IOT3(L2,L2,K,0,0,0,B)
      CALL IOT6(L,L2,L3,L2,L,K,C)
      CALL IOT3(L,L,K,0,0,0,A)
      F3=(-1.)**L3*DFLOAT((2*L+1)*(2*L2+1))*A*B*C
      RETURN
      END

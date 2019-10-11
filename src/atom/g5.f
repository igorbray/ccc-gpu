      REAL FUNCTION G5(K,L,L2,L3)
      REAL*8 A
      CALL IOT3(L,K,L2,0,0,0,A)
      G5=DFLOAT(2*L+1)*A*A
      IF(L3.EQ.K) G5=G5*(1.-DFLOAT(2*L2+1)/
     *DFLOAT(2*L3+1))
      RETURN
      END

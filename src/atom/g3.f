      REAL FUNCTION G3(K,L,L2,L3,S1)
      REAL*8A,B
      INTEGER S1(1)
      CALL IOT3(L,K,L2,0,0,0,A)
      G3=DFLOAT((2*L+1)*(2*L2+1))*A*A
      B=0.0
      IF(L3.EQ.K.AND.S1(1).EQ.1) B=2./(2.*K+1.)
      G3=G3*(1./DFLOAT(2*L2+1)-B)
      RETURN
      END

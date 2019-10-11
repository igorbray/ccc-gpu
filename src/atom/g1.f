      FUNCTION G1(K,L,L2,Q,Q2)
      INTEGER Q,Q2
      REAL*8 A
      CALL IOT3(L,K,L2,0,0,0     ,A)
      G1 = FLOAT( Q*Q2)/2.*real(A*A)
      RETURN
      END

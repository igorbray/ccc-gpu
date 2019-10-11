      SUBROUTINE MXMLT(R,S,AA,N)
      REAL*8 R(N , N),S(N , N),AA(N  ) ,C
      DO 1 I=1,N
      DO 2 J=1,N
      C=0.
      DO 3 K=1,N
    3 C=C+R(I,K)*S(K,J)
    2 AA(J)=C
      DO 1 K=1,N
      R(I,K)=AA(K)
    1 CONTINUE
      RETURN
      END

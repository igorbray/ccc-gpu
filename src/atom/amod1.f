      REAL    FUNCTION AMOD1(I,J,K)
      real*8 ddff
      ddff=DFLOAT( IFACT(I+J-K)*IFACT(I+K-J)*
     *IFACT(K+J-I))/DFLOAT(   IFACT(I+J+K+1))
      AMOD1=SQRT(ddff)
      RETURN
      END

      SUBROUTINE  ZERO(R5,Z,R8,TETA,
     * NE,IS,IN,IQ,NE5,NE5IS,K)
      INTEGER IN(IS),IQ(IS)
      REAL*8Z,SUM, TETA(IS),R5(NE5 ),R8(NE5IS)
      DO 1 II=1,NE5
    1 R5(II)=0.
      R5(2)=1.
      SUM=0.
      DO 2 II=1,IS
      R5(1)=((Z-SUM)/DFLOAT(IN(II)))**2
      IF(K.EQ.3.AND.II.NE.IS) GOTO 455
      MSK= NE5 * (II-1)
      DO 3 IW=1,NE5
    3 R8(IW+MSK)=R5(IW)
455   SUM=SUM+DFLOAT(   IQ(II))
    2 TETA(II)=1.
      RETURN
      END

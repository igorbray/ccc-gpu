      SUBROUTINE INVERT(A,  EPS,DET,N,ierror)
      REAL*8    A(N , N),EPS,DET,Q
      integer ierror

      ierror=0
      DET=1.
      DO 1 I=1,N
      Q=A(I,I)
      DET=DET*Q
      IF(ABS(Q).GT.EPS)  goto 5
        ierror=1
        RETURN
5     A(I,I)=1.
      DO 2 K=1,N
2     A(I,K)=A(I,K)/Q
      DO 10 J=1,N
      IF(I.EQ.J) GOTO 10
      Q=A(J,I)
      A(J,I)=0.
      DO 100 L=1,N
      A(J,L)=A(J,L)-A(I,L)*Q
 100    CONTINUE
 10     CONTINUE
1     CONTINUE
      RETURN
       END

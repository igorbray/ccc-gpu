      SUBROUTINE INVER(A,  EPS,DET,N,*,NNN)                             F4503350
      REAL*8    A(NNN ,NNN),EPS,DET,Q                                   F4503360
      DET=1.                                                            F4503370
      DO 1 I=1,N                                                        F4503380
      Q=A(I,I)                                                          F4503390
      DET=DET*Q                                                         F4503400
      IF( ABS(Q).LT.EPS) RETURN1                                        F4503410
      A(I,I)=1.                                                         F4503420
      DO 2 K=1,N                                                        F4503430
C       IF(A(I,K).LT.1.D-20) PRINT 7,A(I,K),Q                           F4503440
C7       FORMAT(5X,'MA‹OE —ˆC‹O',5X,E15.8,'      „E‹ˆM HA ',E15.8)      F4503450
2     A(I,K)=A(I,K)/Q                                                   F4503460
      DO 10 J=1,N                                                       F4503470
      IF(I.EQ.J) GOTO 10                                                F4503480
      Q=A(J,I)                                                          F4503490
      A(J,I)=0.                                                         F4503500
      DO 100 L=1,N                                                      F4503510
      A(J,L)=A(J,L)-A(I,L)*Q                                            F4503520
 100    CONTINUE                                                        F4503530
 10     CONTINUE                                                        F4503540
1     CONTINUE                                                          F4503550
      RETURN                                                            F4503560
       END                                                              F4503580

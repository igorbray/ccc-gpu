      FUNCTION IC(NN,MM)
      IF(NN.GT.0.AND.MM.GT.0.AND.NN.LT.MM) write(6, 777) NN,MM
777   FORMAT(' ',88('*')/'  In subr.   IC ',
     * '-- Calculation of  chiClA COchETAHiii- ',
     * ' C from N  by  M ---is not allowed ',
     * ' Parameters        N=',I5,'   M=',I5)
      IF(NN.NE.0.OR.MM.NE.0) GOTO 788
      IC=1
      RETURN
788   CONTINUE
      IF(NN.LE.0.AND.MM.GT.0) write(6, 777) NN,MM
      IF(MM.LT.0) write(6, 777) NN,MM
      N=NN
      M=MM
      IF(M.LT.N-M) M=N-M
      N=N+1
      IA=1
      IF(M.LT.1) GOTO 2
      DO 1 I=1,M
1     IA=(IA*(N-I))/I
2     IC=IA
      RETURN
      END


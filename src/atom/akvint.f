      REAL   FUNCTION AKVINT (XX,X,Y,IN,IK,NE5)
C   KBAÑPATàóHAü  àHTEPèOãüñàü
      REAL*8XX,X(NE5 ),Y(NE5 ),T,Z
      K=IN
      L=IK-2
   10 I=(K+L-1)/2.0+0.6
      IF(I.EQ.K) GOTO 20
      IF(XX.GE.X(I)) GOTO 1
      L=I
      GOTO 10
    1 K=I
      GOTO 10
   20 T=XX-X(K)
      Z=((XX-X(K+1))/(X(K)-X(K+2)))*T
      AKVINT=Y(K)+(T+Z)/(X(K)-X(K+1))*(Y(K)-Y(K+1))-
     *(Z/(X(K+1)-  X(K+2)))*(Y(K+1)-Y(K+2))
      RETURN
      END

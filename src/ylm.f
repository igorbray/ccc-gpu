      function YYLM(l, m, theta )                                         

      REAL*8 plm,dcos
      real*8 fl

      common/fl0/fl(0:300)
      pi = acos(-1.0)
      if ((m .lt. 0) .or. (m .gt. l)) then
         YLM = 0.0
         return
      end if
      dcos = cos(theta)
      aux = PLM(l, m, dcos)
      tmp = sqrt(float(2*l+1)/4/pi) * aux * 
     >         exp((FL(l-m) - FL(l+m)) * 0.5)
      YYLM = tmp
      return
      end

      FUNCTION PLM(L,M,X)
C
C     CALCULATE ASSOCIATED LEGENDRE POLY FROM A RECURSION RELATION
C     IN MESSIAH. NOTE: M CAN'T BE NEGATIVE
C

      IMPLICIT REAL*8 (A-H,O-Z)
*      IMPLICIT integer (i-n)



      IF(ABS(X).GT.1.0D0) THEN
                             WRITE(6,256) X
256   FORMAT(' WARNING : ERROR IN PLM   X MUST BE BETWEEN [-1,1]',
     *        '   HERE X=',G22.16)
                            CALL EXIT
                        END IF
      IF(L.LT.0) WRITE(6,1) L
1     FORMAT(' ? L IN PLM = ',I6/)
      IF(M.LT.0) WRITE(6,2) M
2     FORMAT(' ? M IN PLM = ',I6/)
      IF(M.GT.L) WRITE(6,3) L,M
3     FORMAT(' ? L = ',I6,' & M = ',I6,' IN PLM'/)
C*********
      IF(L.NE.0) GO TO 4
      PLM=1.0D0
      RETURN
C*********
4     IF(L.NE.1) GO TO 6
      IF(M.NE.0) GO TO 5
      PLM=X
      RETURN
C*********
5     PLM=DSQRT(1.0D0-X*X)
      RETURN
C*********
C HERE  L.GE.2
C
6     IF(M.NE.0) GO TO 7
      A=1.0D0
      B=X
      DM=M
      GO TO 10
C*********
C    CALCULATES PLM(M,M,X)
C
7     DF=1.0D0
      ITMONE=2*M-1
      DO 8 LL=1,ITMONE,2
      DLL=LL
8     DF=DF*DLL
      OMX=DSQRT(1.0D0-X*X)**M
      A=OMX*DF
      IF(L.NE.M) GO TO 9
      PLM=A
      RETURN
C**********
9     DM=M
      B=A*(2.0D0*DM+1.D0)*X
      IF(L.NE.M+1) GO TO 10
      PLM=B
      RETURN
C**********
10      MPTWO=M+2
       DO 11 LL=MPTWO,L
      DLL = DBLE(LL)
      PLM=((2.0D0*DLL-1.0D0)*X*B-(DLL+DM-1.0D0)*A)/(DLL-DM)
      A=B
      B=PLM
11    CONTINUE


      RETURN
      END
      



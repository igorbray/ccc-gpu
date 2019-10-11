      SUBROUTINE PRINTM(AB,NPHS)
      REAL*8    AB(NPHS)
      WRITE(1 ,1)(AB(J),J=1,NPHS)
    1 FORMAT(5(2X,E13.6,3X))
      RETURN
      END
      REAL  FUNCTION BETA(RE,IM,RE1,
     * IM1,A,B,   J,J1,I1)
      REAL*8 IM,IM1,RE,RE1,A,B,Q,Q1,PR,PI,ddff
      Q=RE*RE+IM*IM
      Q1=RE1*RE1+IM1*IM1
      PR=RE*RE1+IM*IM1
      PI=RE*IM1-IM*RE1
      ddff=DFLOAT(     I1*(I1+1))
      BETA=(DFLOAT(J)*Q+DFLOAT(J1)*Q1+6.*
     * SQRT(ddff)*
     *(PR*A-PI*B))/(DFLOAT(   2 *I1+1 )*(Q+Q1))
      RETURN
      END

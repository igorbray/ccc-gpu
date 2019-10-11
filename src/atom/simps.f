      SUBROUTINE  SIMPS(ALFA,BET,H,P,R1,
     *R2,R5,R6,NE,IY)
      REAL*8 R2(ne+5),R1(ne+5),R5(ne+5),
     *R6(ne+5) ,ALFA,BET,H,P
      IF(IY.EQ.1) GOTO 10
      DO 1 II=1,NE
      R6(II+1)=R2(II+2)*R5(II+2)*R1(II+2)*
     *R6(II+1)/(ALFA*R1(II+2)+BET)**2
    1 CONTINUE
   10 P=R6(2)+4.*R6(3)
      IRA=NE-3
      DO 2 II=3,IRA,2
    2 P=P+R6(II+1)*2.+4.0*R6(II+2)
      P=(P+R6(NE))*H/3.
      RETURN
      END

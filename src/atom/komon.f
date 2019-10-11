      SUBROUTINE KOMON( Z,R,BET,H,GAM,EPS,ALFA,
     *RO1,RO2,NE,IN,IL,IQ,NU,KK,ISIGM1,ISIGM2,IS,L1)
C
C  Print of entered data,  Calculation of   ALFA
C
      REAL*8ALFA,BET,H,Z,RO1,RO2,EPS,GAM(NU),R
      INTEGER NE, KK(NU),ISIGM1(NU),ISIGM2(NU),
     *IN(IS),IL(IS),IQ(IS)
      write(66, 5) Z
    5 FORMAT('Nuclear charge' ,22X,F4.1,//
     */'Configuration:'/    ,
     *  '  Principal quantum number',
     */'',
     *12X,'electon moment'    )
      DO 6 IZCH=1,IS
c    6  7,(IN(IZCH),IL(IZCH),IQ(IZCH))
c!!!
    6 write(66, 7) IN(IZCH),IL(IZCH),IQ(IZCH)
    7 FORMAT(' ',3X,I5,12X,I5,7X,I5)
      IF(L1.EQ.111)  GOTO 111
      write(66, 8)
   8  FORMAT(//' ','Coefficiants: (Calculated)')
      GOTO 112
111   write(66, 88)
  88  FORMAT(//' ','Coefficiants: (Entered)')
112   write(66, 87)
   87 FORMAT(' value ',   16X,'Entetaction'
     */' coefficiant',10X,' between shells'
     */21X,'K',3X,
     *' ISIGM1','    ISIGM2')
      IF(NU.LT.1) GOTO 21
      DO 9 IZCH=1,NU
    9 write(66, 10) GAM(IZCH),KK(IZCH),
     *ISIGM1(IZCH),ISIGM2(IZCH)
   10 FORMAT(' ',E15.8,1X,I5,3X,I5,5X,I5)
   21 CONTINUE
   11 FORMAT(F15.8)
      write(66, 12) EPS
   12 FORMAT(/'precision',15X,F15.8)
      RO1=(-DFLOAT(10)-LOG(Z))*BET
      RO2=RO1+DFLOAT(  NE)*H
      ALFA=(RO2-BET*LOG(R))/R
      write(66, 13) R,NE,H,BET,ALFA
   13 FORMAT(' R MAX',22X,E15.8/
     *'the number of integration points' ,
     *I4/'integration step',9X,
     *E15.8/4H BET,24X,E15.8/5H ALFA,
     *23X,E15.8)
      RETURN
      END

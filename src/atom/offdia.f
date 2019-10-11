      SUBROUTINE OFFDIA(ALFA,BET,H,R1,R2,
     *R5,R8,ALA,ALA1,ALA2,ALM3,TETA,C3,UD,
     *NE,IS,NE5,ISIGM,IN,IL,KF)
      include 'paratom.f'
      REAL*8 ALFA,BET,H,R1(NE5),R2(NE5) ,
     *R5(NE5),ALA(IS),ALA1(IS),ALA2(IS),R8(nsizeis*nsize),
     *ALM3(IS),TETA(IS),C3(IS),TETMAX,C,D,P,A1,PRIN
      INTEGER UD(IS),IN(IS),IL(IS)
c$$$      write(*,*) ' OFFDIA'
      IL1=IL(ISIGM)
      IN1=IN(ISIGM)
      IRA=IS-1
      DO 30 J=1,IRA
      TETMAX=0.0
      IF(.NOT.(IL1.EQ.IL(J).AND.
     *IN1.NE.IN(J))) GOTO 301
      IF(UD(J).NE.1) GOTO 301
      IA=J
      MMM=(IA-1)*NE5
      DO 2002 I100=1,NE5
2002  R5(I100)=R8(MMM+I100)
      C=0.0
      D=0.0
      R5(NE+2)=C
      DO 10 IW=1,NE
      P=D
      D=C
      C=R2(IW+2)*R5(IW+2)*(R1(IW+2)/
     *(ALFA*R1(IW+2)+BET))**2
      C=C/SQRT(ALM3(IS))
      A1=(5.*C+8.*D-P)*H/12.
      TETMAX=A1+TETMAX
   10 CONTINUE
      WRITE( 4,987)J,TETMAX,ALM3(J)
987   FORMAT(' overlap integral for',
     * I3,' shell=',E15.8,' for non-diagonal terms',E12.5)
      IF(ABS(TETMAX).LT.1.D-07.AND.
     * KF.NE.1) GOTO 301
      IF(KF.NE.1)  GOTO 122
      C3(J)=ALM3(J)+TETMAX
      ALA(J)=0.005
      ALA1(J)=TETMAX
      GOTO 123
122   CONTINUE
      IF(TETMAX *TETA(J).GT.0.0) GOTO 124
      ALA(J)=ALA2(J)
      ALA1(J)=TETA(J)
      IF(ABS(TETMAX).LT.0.001) GOTO 124
      C3(J)=(ALM3(J)+ALA2(J))/2.
      GOTO 123
124   C3(J)=(ALM3(J)*TETA(J)-
     * ALA2(J)*TETMAX)/(TETA(J)-TETMAX)
      IF(ABS(C3(J)-ALM3(J)).GT.5.*
     * ABS(ALM3(J)-ALA2(J)).AND.ABS(TETMAX).GT.0.003)
     * C3(J)=3.*ALM3(J)-ALA2(J)
      IF(KF.LT.10) GOTO 370
      PRIN=(TETA(J)-TETMAX)/(ALA2(J)-ALM3(J))
      IF(PRIN.LT. -0.0005) GOTO 370
      C3(J)=ALM3(J)+SIGN(0.1d0,TETMAX)
      IF(ABS(ALA(J)-C3(J)).GT.0.1) C3(J)=
     * (C3(J)+ALA(J))/2.
      ALA(J)=ALA(J)-SIGN(0.1d0,TETA(J))
 370  CONTINUE
123   CONTINUE
      IF(ABS(C3(J)).GT.20.) C3(J)=-C3(J)/10.
      ALA2(J)=ALM3(J)
      ALM3(J)=C3(J)
301   TETA(J)=TETMAX
30    CONTINUE
      RETURN
      END

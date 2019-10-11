      SUBROUTINE ACCOBM(R2,R5,ALM4,
     *ALM5,ALM6,NE,NE5,NE5IS,IP,KF)
      REAL*8 R2(NE5),R5(NE5),ALM4,
     *ALM5,ALM6,U,V,A1,B,C1
c$$$      write(*,*) 'ACCOMB'
      U=0.0
      V=U
      A1=U
      DO 22 J=1,NE
      B=R2(J+2)-R5(J+2)
      IF(ABS(B).LE.ABS(A1)) GOTO 2277
      A1=B
      U=R2(J+2)
      V=R5(J+2)
      JJ=J+2
2277  continue
   22 CONTINUE
      ALM6=A1
      IF(ABS(U).GT.1.0) ALM6=ALM6/U
      IF(KF.NE.1) GOTO 200
      ALM4=ALM6
      ALM5=ALM6
200   CONTINUE
C     C1=0.7
C     GOTO 57
      IF(SIGN(1d0,ALM5-ALM4).
     * NE.SIGN(1d0,ALM6
     *-ALM5).AND.ABS(ALM4).GT.ABS(ALM5)) GOTO 56
      C1=0.01
      IF(ABS(ALM6-ALM4).LT.0.001.AND.ABS(ALM6).GT.0.05.
     * AND.KF.GT.5) C1=0.5
      GOTO 57
   56 C1=(ALM6-ALM5)/(ALM6-2.*ALM5+ALM4)
   57 CONTINUE
      DO 307 IW=1,NE
      R2(IW+2)=C1*R5(IW+2)+(1.-C1)*R2(IW+2)
  307 CONTINUE
577   ALM4=ALM5
      ALM5=ALM6
      IF(IP.EQ.-1) R2(1)=C1*R5(1)+(1.-C1)*R2(1)
c     WRITE( 4,98)KF,ALM6,U,V,JJ,C1
98    FORMAT(' *KF=',I3,'  MAX HEB. Bf =',E10.3,'  U=',E10.3,
     * '  V=',E10.3,' JJ=',I4,'  C1=',E10.3)
      RETURN
      END

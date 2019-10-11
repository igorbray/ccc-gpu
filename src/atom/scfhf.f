      SUBROUTINE SCFHF(Z,GAM,DIR,R,BET,EPS,H,TETA,
     *ALM3,R1,R2,R3,R4,R5,R6,R7,R8,MU,C3,ALA,ALA1,
     * ALA2,ALM6DU,
     *DZ,KZ,
     *IS,IN,IL,IQ,
     *SS1,L1,NU,IG,KK,ISIGM1,ISIGM2,IP,PP,NE,
     *NE5,NE5IS,KT,UD,KSU1,gamma,r0)
c !!!
      real *8 r
c
      include 'paratom.f'
c
      REAL*8ALFA,BET,H,Z,RO1,RO2,EPS,
     *GAM(NSIZEqf ),TETA(IS),B,C,D,P, DZ,ZF,DDZ,
     * T,U,V,GAMA,TAU1,TAU2,ALM1,ALM2,
     *TETMAX,ALM3(IS), MU,ALM6DU(IS),
     * R1( NE5),R2(NE5 ),R3( NE5),
     * R4( NE5),R5( NE5),R6( NE5),DIR(IS),
     * R7( KT),R8(NE5IS),A1,B1,A2,OMEGA ,
     * ALA(IS),C3(IS),ALA1(IS),ALA2(IS),EP
      INTEGER CH,KK(NSIZEqf),ISIGM1(NSIZEqf),PP(IS),
     *ISIGM2(NSIZEqf),IN(IS),IL(IS),IQ(IS), UD(IS),
     *SS1,S1(max10),TT,G,W,JJ,A,SS,IP(IS),L1

C Only one shell. No need for the coefficients.
      if (is .eq. 1) then
         ig = 0
         gam(1) = 0
         goto 111
      end if      
      IF(L1.EQ.111) GOTO 111
      DO 5 I=1,IS
5     IP(I)=1
      CH=1
      I1=1
      S1(1)=SS1
      DO 11 I=1,IS
      A=2*IL(I)+1
      IF(S1(1).LE.99) A=A*2
      DO 11 J=I1,IS
      JJ=2*IL(J)+1
      IF(S1(1) .LE.99)  JJ=JJ*2
      IF(I.LT.J ) CH=CH+     ( IL(I)+IL(J)-
     *IABS(IL(I)-IL(J)))/2 +1
      IF(I.EQ.J.AND.IQ(I).NE.1) CH=CH+IL(I)
      IF(KSU1.NE.1.AND.(IQ(I).EQ.1.AND. !KSU changed to KSU1, Igor 10/4/2011
     * IQ(J).EQ.JJ-1).OR.
     * (IQ(J).EQ.1.AND.IQ(I).EQ.A-1))
     * CH=CH+ MIN0(IL(I),IL(J))
11     CONTINUE
       CALL COEF(GAM,DIR,IS,1,KSU1,IN,IL,IQ,IP,
     * L1,S1,CH,KK,ISIGM1,ISIGM2,IG,1 )
c$$$      print*, 'gam', (gam(i), i=1,nu)
      NU=CH
111   write(66, 1)
    1 FORMAT('+','Wave functions of the ground state of the atom ',
     *''/
     *' in Hartree-Fock approximation')
      write(66, 4) IS,NU,IG
    4 FORMAT('   The number of shells           ',I5,
     * /     '   The number of coefficients     ',I5,
     * /     '   The number of exchange coeff.  ',I5)
      CALL KOMON( Z,R,BET,H,GAM,EPS,ALFA,RO1,
     *RO2,NE,IN,IL,IQ,NU,KK,ISIGM1,ISIGM2,IS,L1)
      CALL VAR1(Z,R,BET,H,ALFA,RO1,RO2,R1,NE,NE5)
      write(66, 777) R1(3)
777   FORMAT(' R first',22X,E15.8  )
      write(66, 778) DZ,KZ
778   FORMAT(' DZ' ,25X,F6.2/' KZ',25X,I4//)
      IF(KZ.NE.0) GOTO 750
      ZF=Z
      KZ=1
      DDZ=0.
      GOTO 751
750   ZF=Z+DZ
      DDZ=DZ/DFLOAT(KZ)
          KZ=KZ+1
751   CONTINUE
c!!!
      CALL ZERO( R5,ZF,R8,TETA,NE,IS,IN,IQ,
     * NE5,NE5IS,KSU1)
      ZF=ZF+DDZ
      DO 779 IKZ=1,KZ
      ZF=ZF-DDZ
      write(66, 760) ZF
760   FORMAT('  Z=',F6.2)
      EP=1.E-03
      IF(IKZ.EQ.KZ) EP=EPS
      CALL DH1(ALFA,BET,H,ZF,Z,RO1,RO2,EP ,R,
     *GAM,TETA,ALM3,R1,R2,R3,R4,R5,R6,R8,MU,
     * C3,ALA,ALA1,ALA2, ALM6DU,
     *NE,NE5,NE5IS,IS,NU,IZCH,IG,CH,KK,ISIGM1,
     *ISIGM2,IN,IL,IQ,UD,KSU1,gamma,r0)
779   CONTINUE
c!!!

      CALL INTZZZ( R5,R8,RO1,BET,ALFA,R7,H,R2,
     * IS,NE,NE5,NE5IS,KT,IN,KSU1,r1)
      RETURN
      END

      SUBROUTINE PHOABS(  EHL,  DKAP,KAPO,
     * EE,ML,MV,PHS,U,V,R,FUNC,F,AU,COEF,Y,
     * PF,PB,SC,S,SUMO,SUML,SUMM,SUMV,RD,Q,AA,
     * EXPI,NN,NCH,HL,ND,NES,
     * LH,LP,PQN,NNL,N0,SCH,SPT,IP,SPI,JJBEG)
      REAL*8 EHL(1),CNST, DKAP(1),KAPO(1),
     * EE(1),ML(1),MV(1),PHS(1),U(1,1),V(1,1),
     * R(1),FUNC(100,100),F(1),AU(1),COEF(1),Y(1,1),
     * PF(1),PB(1),SC(1),S(1,1),SUMO(1),SUML(1),
     * SUMM(1),SUMV(1),RD(1,1),Q(1,1),AA(1),
     * H,BETA,ALPHA,EH,Z
      INTEGER  EXPI,NN,NCH,HL(1),ND(1),NES(1),
     * LH(1),LP(1),PQN(1),NNL(1),N0(1),SCH(1),
     * SPT(1),IP(1),SPI(1)
      MT=3
      CNST=8.067
       NHL=0
      DO 6 I=1,NCH
      K=HL(I)
        IF(I.LT.2) GOTO 200
      DO 7 J=2,I
      IF(HL(J-1).EQ.K) GOTO 6
7     CONTINUE
200   NHL=NHL+1
6     CONTINUE
      NPHS=0
      DO 10 J=1,NCH
10    NPHS=NPHS+NES(J)
      NWF=NHL+NPHS
      WRITE(6, 11)
11    FORMAT(/10X,'calculation of the cross section ',
     * 'of photon absorption'//)
      ISTOP=0
      CALL KOMON1(H,BETA,ALPHA,R,EE,FUNC,PQN,EHL,
     * KAPO,DKAP,PHS,F,Z,MT,NCH,HL,EXPI,NES,
     *ND,    NNL,NN,LP,LH,NWF,NPHS,NHL,N0,NN+5,ISTOP)
      IF(ISTOP.EQ.1) RETURN
      WRITE(6, 90)
 90   FORMAT(/'Hartree- fock approximation')
      CALL MAEL1(EE,CNST,BETA,ALPHA,FUNC,H,R,MV,
     * ML,F,AU,
     * NNL,NCH, NES,HL,ND,LH,LP,NN,NHL,NWF,NPHS)
      WRITE(1,3)
3     FORMAT(5X,'matrix elements of ',
     * 'electron-hole coulomb interaction ')
      DO 9 IA=1,NCH
      DO 9 IB=1,IA
      CALL MAEL2(H,ALPHA,BETA,R,FUNC,
     *  U,V, COEF,Y,AU,1,IA,IB,LH,LP,HL,NES,NNL,
     * NN,NWF,NPHS,NHL,
     * NCH,SPI)
9     CONTINUE
      CALL MAELCN( KAPO,DKAP,EE,
     * U,V,PHS,ML,MV,CNST,EHL,
     * PF,PB,SC,S,SUMO,SUML,SUMM,SUMV,
     * RD,Q,AA,Z,SCH,SPT,IP,
     *  NCH,NPHS,ND,NES,LP,HL,LH,NNL,NPHS*2,PQN)
c !!!
c     ,JJBEG)
      RETURN
105   CONTINUE
       RETURN
      END

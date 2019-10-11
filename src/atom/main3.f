       PARAMETER (N30=20)
       PARAMETER (N2=5)
       PARAMETER (N60=N30 * 2)
       PARAMETER (N705=705)
       PARAMETER (N10=N30 + N2)
       PARAMETER (N200=200)
       PARAMETER (N20=20)
c
      character*15  con,  dum,  osso, voso
c
       integer  nch, nn, expi, JJBEG
       real*8 ehl(N30), dkap(N30), kapo(N30), sumo(N30),
     * suml(N30), summ(N30), sumv(N30)
       integer  hl(N2), nd(N2), nes(N2), lh(N2), pqn(N2),
     * nnl(N2), n0(N2), sch(N2), ip(N2), spi(N2), nphs, nphs22,
     *   nwf, ne5      , SPT(N60)
       real*8 EE(N30), ML(N30), MV(N30), PHS(N30), PF(N30), PB(N30),
     * SC(N30), U(N30,N30), V(N30,N30), R(N705), AU(N705), F(N705),
     * FUNC(N10,N705), COEF(N200), Y(N20,N705), S(N60,N60),RD(N60,N60),
     * Q(N60,N60), AA(N60)
c
      open(7,file='job3.dat')
      read(7,222) con
222   format(15a)
      read(7,222) dum
      read(7,222) osso
      READ(7,*) NCH
      READ(7,*) NN
      READ(7,*) EXPI
      READ(7,*) JJBEG
      ii = 100
      DO 1 I=1,NCH
      read(7,222) voso
      if(i . EQ . 1)  ii = 3
      if(i . EQ . 2)  ii = 4
      if(i . EQ . 3)  ii = 10
      if(i . EQ . 4)  ii = 8
      if(i . EQ . 5)  ii = 9
      open(ii,file=voso,form='unformatted')
      READ(7,*) HL(I)
      READ(7,*) nd(I)
      READ(7,*) NES(I)
      READ(7,*) SPI(I)
      IF(EXPI.NE.0) READ(7,*) EHL(I)
1     CONTINUE
c
      NPHS=0
      DO 2 I=1,NCH
2     NPHS=NPHS + NES(I)
      NPHS22= 2*NPHS
      NWF=NPHS+NCH
      NE5=NN+5
c
        close(7)
c
      open(2,file=osso,form='unformatted')
      open(1,file=dum)
      open(6,file=con)
c
      CALL PHOABS(EHL,DKAP,KAPO,EE,ML,MV,PHS,U,V,R,FUNC,F,AU,COEF,
     * Y,PF,PB,SC,S,SUMO,SUML,SUMM,SUMV,RD,Q,AA,
     * EXPI,NN,NCH,HL,ND,NES,LH,LP,PQN,NNL,N0,SCH,SPT,IP,SPI,JJBEG)
      STOP
      END



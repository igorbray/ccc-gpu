      SUBROUTINE CQD2(OM,EE,G,SC,OM1,EHL,
     * KAPO,GP,U,V,DKAP,S,R,PF,PB,Q,AA,IP,
     * NCH,NPHS,NES,N,ND,CONT,
     * BOO,ierror,NPHS02)
      integer ierror
      INTEGER BOO,CONT,IP(  3),
     * ND(NCH),NES(NCH)
      REAL*8 KAPO(NCH),EE(NPHS),
     * PF(NPHS),SC(NPHS)
      REAL*8 OM,G,OM1,GP,PB(NPHS),
     * EHL(NCH),  A,C,B,X,D,Y, AA(1),
     * S(NPHS02,NPHS02 ),U(NPHS ,NPHS) ,V(NPHS ,NPHS),
     * DKAP(NCH),Q(NPHS02,NPHS02) ,R(NPHS02,NPHS02)
      ierror=0
      DO 1 I=1,NCH
1     IP(I)=0
      DO 2 J=1,NPHS
      IF(ABS(OM-EE(J)).GE.0.1D-5) GOTO 3
      K=0
      I=0
4     K=K+1
      I=I+NES(K)
      IF(I.LT.J) GOTO 4
      IP(K)=J
      PF(J)=0.0
      GOTO 5
3     PF(J)=G*SC(J)/(OM1-EE(J))
5     PB(J)=-G*SC(J)/(OM1+EE(J))
2     CONTINUE
      WRITE(1 ,6)
6     FORMAT(3X,'BACK. PROPAGATOR')
      CALL PRINTM(PB,NPHS)
      WRITE(1 ,7)
7     FORMAT(3X,'FORW.PROPAGATOR')
      CALL PRINTM(PF,NPHS)
      DO 8 I=1,N
      DO 8 J=1,N
8     S(I,J)=0.0
      BOO=0
      J2=0
      DO 9 I=1,NCH
      J1=J2+ND(I)+1
      J2=J2+NES(I)
      IPI=IP(I)
      A=-EHL(I)
      C=KAPO(I)
      B=EE(J2)+C*C
      IF(OM1.LE.A.OR.OM1.GE.B) GOTO 50
      X=G*LOG((B-OM1)/(OM1-A))
        IF(J1.GT.J2) GOTO 571
      DO 10 J=J1,J2
10    X=X+PF(J)
571   IF(IPI.NE.0) GOTO 11
      C=ABS(OM-EE(J1))
      LL=J1+1
      DO 12 L=LL,J2
      D=C
      C=ABS(OM-EE(L))
      K1=L-1
      IF(C.GT.D) GOTO 13
12    CONTINUE
      GOTO 50
13    IP(I)=-K1
      BOO=1
      J3=K1-2
      IF(J3.LT.J1) J3=J1
      J4=K1+2
      IF(J4.GT.J2) J4=J2
      DO 14 J=J3,J4
      CALL INTER1(OM1,EE,Y,J,J3,J4)
      PF(J)=PF(J)-X*Y
      DO 14 K=1,NPHS
      K1=K+NPHS
      S(K,J)=S(K,J)+GP*Y*U(K,J)
      S(K1,J)=S(K1,J)+GP*Y*V(K,J)
14    CONTINUE
      GOTO 50
11    CONTINUE
      IF(IPI.LT.J1.OR.IPI.GT.J2)  GOTO 99
      BOO=1
      PF(IPI)=-X
          D=-G*SC(IPI)/(24.*DKAP(I)*(KAPO(I)+
     * DFLOAT (IPI-J1)*DKAP(I)))
      IF(IPI.EQ.J1.OR.IPI.EQ.J2) GOTO 15
      GOTO 16
15    I1=1
      IF(IPI.EQ.J2) I1=-1
      D=DFLOAT(   I1)*D
      J=IPI
      PF(J)=PF(J)-25.*D
      J=J+I1
      PF(J)=PF(J)+48.*D
      J=J+I1
      PF(J)=PF(J)-36.*D
      J=J+I1
      PF(J)=PF(J)+16.*D
      J=J+I1
      PF(J)=PF(J)-3.*D
      GOTO 17
16    CONTINUE
      IF(IPI.EQ.J1+1.OR.IPI.EQ.J2-1) GOTO 19
      GOTO 20
19    I1=1
      IF(IPI.EQ.J2-1) I1=-1
      D=DFLOAT(  I1)*D
      J=IPI-I1
      PF(J)=PF(J)-3.*D
      J=J+I1
      PF(J)=PF(J)-10.*D
      J=J+I1
      PF(J)=PF(J)+18.*D
      J=J+I1
      PF(J)=PF(J)-6.*D
      J=J+I1
      PF(J)=PF(J)+D
      GOTO 17
20    J=IPI-2
      PF(J)=PF(J)+D
      J=IPI-1
      PF(J)=PF(J) -8.*D
      J=IPI+1
      PF(J)=PF(J)+8.*D
      J=IPI+2
      PF(J)=PF(J)-D
17    CONTINUE
      DO 25 J=1,NPHS
      S(J,IPI)=S(J,IPI)+GP*U(J,IPI)
      J1=J+NPHS
25    S(J1,IPI)=S(J1,IPI)+GP*V(J,IPI)
50    CONTINUE
9     CONTINUE
      IF(CONT.EQ.1.AND.NCH.NE.1) WRITE(1 ,26)
26    FORMAT(' FORW.PROP. WITH PRINC. VAL. CORR.')
      CALL PRINTM(PF,NPHS)
      DO 27 J=1,NPHS
      A=PF(J)
      B=PB(J)
      J1=J+NPHS
      DO 27 I=1,NPHS
      C=0.
      IF(I.EQ.J) C=1.
      I1=I+NPHS
      X=C-A*U(I,J)
      Q( I,J)=X
      R( I,J)=X
      X=-A*V(I,J)
      Q(I1,J)=X
      R(I1,J)=X
      X=-B*V(I,J)
      Q(I,J1)=X
      R(I,J1)=X
      X=C-B*U(I,J)
      Q(I1,J1)=X
      R(I1,J1)=X
27    CONTINUE
      CALL INVERT(R  ,1.D-5,A,N,ierror)
      if( ierror .eq. 1 )goto 99
      WRITE(66, 30) A
30    FORMAT('+',43X,'DET(A)=',F8.5)
      IF(BOO.EQ.0) GOTO 31
      CALL MXMLT(R,S,AA,N)
      CALL MXMLT(S,R,AA,N)
      DO 32 I=1,N
      DO 32 J=1,N
      C=R(I,J)
      R(I,J)=Q(I,J)+S(I,J)
      S(I,J)=C
32    CONTINUE
      CALL INVERT(R,  1.D-5,A,N,ierror)
      if( ierror .eq. 1 )goto 99
      WRITE(66, 33) A
33    FORMAT('+',73X,'DET(1/R)=',F8.5)
      CALL MXMLT(S,R,AA,N)
31    CONTINUE
      RETURN
99    CONTINUE
        WRITE(66, 31311)
31311  FORMAT( 5X,15('*'),'Determinant  =0',15('*'))
      ierror=1
      RETURN
      END

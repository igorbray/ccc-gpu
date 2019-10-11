      program psdrive
      include 'par.f'
      real waveout(maxr,nnmax),erydout(nnmax)
      dimension jdouble(100)
*
      hmax = 0.03488889
      ra = 100.0
      njdouble = 11
      jdouble(1) = 1
      do 10 i=2,10
       jdouble(i) = (i-1)*32
 10   continue
      jdouble(11) = 8000
      zas = 1.0
      la =  1
      nbmax = 20
      call pseudo(jdouble,njdouble,hmax,zas,la,ra,nbmax,maxr,
     >                  erydout,waveout,lastpt)
      do 20 ne=1,nbmax
       write(6,*) ne,erydout(ne)
 20   continue
      print *,'lastpt = ',lastpt
      stop
      end
*======================================================================*
      subroutine pseudo(jdouble,njdouble,hmax,zas,la,ra,nbmax,maxr,
     >                  erydout,waveout,lastpt)
************************************************************************
*                                                                      *
*   THIS ROUTINE DETERMINES PSEUDO BOUND STATES OF A RADIAL POTENTIAL, *
*   I.E., STATES ORBITALS THAT VANISH AT R=0 AND R=A.                  *
*                                                                      *
*   IT IS THE BOUND-STATE PART OF THE CORE POTENTIAL ITERATION         *
*   PROGRAM PUBLISHED IN CHAPTER 2 OF                                  *
*                                                                      *
*            COMPUTATIONAL ATOMIC PHYSICS                              *
*            KLAUS BARTSCHAT  (ED.)                                    *
*            SPRINGER (1996)                                           *
*                                                                      *
*            WRITTEN BY:   KLAUS BARTSCHAT                             *
*                          PHYSICS DEPARTMENT                          *
*                          DRAKE UNIVERSITY                            *
*                          DES MOINES, IOWA 50311, U.S.A.              *
*                                                                      *
*            LAST UPDATE:  MARCH 17, 2003                              *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*
      PARAMETER (NDIM1=4,NDIM2=  9000,NDIM3=100,NNMAX=100,LLMAX=10)
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,SIX=6.0D0)
      PARAMETER (HUGE=1.0D30,SMALL=1.0D-20)
*
      real hmax,zas,ra
      real waveout(maxr,nnmax),erydout(nnmax)
      dimension jdouble(njdouble)
*
      COMMON / PATH / RFULL(NDIM2),RHALF(2*NDIM2),Y(NDIM1,NDIM2)
      COMMON / POT  / EGUESS,CORE(2*NDIM2),RDMASS
*
      DIMENSION IHX(NDIM3),IRX(NDIM3)
      DIMENSION NMINL(0:LLMAX),NMAXL(0:LLMAX)
      DIMENSION VSTART(NDIM1),EBOUND(NNMAX),
     +          OVRLAP(NNMAX,NNMAX)
      DIMENSION FVALUE(NDIM2)
      DIMENSION WAVES(NNMAX,2*NDIM2)
*
**********************  FORMAT STATEMENTS  *****************************
1000  FORMAT(' *******************************************************',
     +     /,' *                                                     *',
     +     /,' *            BOUND-STATE CALCULATION                  *',
     +     /,' *            -----------------------                  *',
     +     /,' *******************************************************',
     +       ///,' NUCLEAR CHARGE:                   ',I3,/,            
     +           ' ASYMPTOTIC CHARGE:                ',I3,/,            
     +           ' ERROR LIMIT FOR BOUND STATE ENERGY:',1P,D12.3,/,     
     +           ' LOWER LIMIT FOR BOUND STATE ENERGY:',1P,D12.3,/)
1001  FORMAT(/,' NEW BOUND STATE FOUND:   EBOUND(',I3,I3,') = ',          
     +       1P,D17.9,//,' NEXT GUESS = ',1P,D14.7,/)
1002  FORMAT(/,' WARNING: OVERLAP INTEGRAL BETWEEN FUNCTIONS WITH', 
     +         ' N1 = ',I2,' AND N2 = ',I2,' FOR L = ',I1,' :',/, 
     +           1P,D15.8,' +/- ',1P,D15.8)
1003  FORMAT(/,' DEBUG PRINT OUT OF FIRST 100 POINTS:',/,               
     +           8X,'R',10X,'FUNCTION',//,1X,1P,2D14.6)
1004  FORMAT(1X,1P,8D14.6)
1005  FORMAT(/,' INODE,NNODE,NBOUND,LBOUND = :',4I5,/,                  
     +         ' ELOW,EHIGH,EGUESS  = :',1P,3D16.8,/)
1006  FORMAT(/,' NITER = ',I5,'  GREATER THAN  NITMAX = ',I5,            
     +         '.  PROGRAM STOPPED.',/)
1007  FORMAT(//,' MESH:',/,' --------------',/,                
     +         ' HINT = ',1P,D14.6,'  NIX =',I2,'  NSTEP  =',I6,        
     +         '  LAST POINT = ',1P,D14.6)
1008  FORMAT(' IHX-ARRAY:',20I5)
1009  FORMAT(' IRX-ARRAY:',20I5)
1010  FORMAT(//,' STARTING BOUND STATE ITERATION FOR (NBOUND,LBOUND) =',
     +         ' (',I2,',',I2,')')
1011  FORMAT(/,' PROBLEMS WITH DEFINITION OF NINTEG. NEND = ',I6,/,     
     +         ' PROGRAM STOPPED.',/)
1012  FORMAT(/,' NORMALIZED ORBITAL FOR  (NBOUND,LBOUND) = (',I2,       
     +         ',',I2,') :',/,'        R          ORBITAL')
1013  FORMAT(1P101D14.6)
1014  FORMAT(/,' LMIN = ',I2,'  LMAX = ',I2,/)
1015  FORMAT(/,' NMIN = ',I2,'  NMAX = ',I2,'  FOR L = ',I2,/)
1016  FORMAT(/,' NITMAX = ',I3,'  IBUG = ',I2,'  ITOUT = ',I2,/)
1017  FORMAT(/,' ESTART = ',1PD14.6,'  ERROR = ',1PD14.6,/)
1018  FORMAT(/,' NUCLEAR CHARGE = ',F6.2,'  REDUCED MASS = ',F7.1,/)
1019  FORMAT(/,' OVERLAP INTEGRALS FOR L = ',I2,/)
1020  FORMAT(/,' NEW BOUND STATE FOUND:   EBOUND(',I3,I3,') = ',          
     +       1P,D17.9,/)
1021  FORMAT(//,' **********   BOUND STATE CALCULATION   **********',//,
     +       ' INPUT IPOT:  1 = COULOMB,  OTHERWISE = NUMERIC',//)
1022  FORMAT(' COULOMB PROBLEM',/)
1023  FORMAT(' NUMERICAL POTENTIAL',/)
1025  FORMAT(' IPOT = ',I2,' CURRENTLY INVALID.  PROGRAM STOPS.',//)
1026  FORMAT(' R0 = ',1PD14.6,'  D = ',1PD14.6,'  ALPHA = ',1PD14.6,0P,
     +       '  REDUCED MASS = ',F7.1,/)
1027  FORMAT(' V0 = ',1PD14.6,'  A = ',1PD14.6,
     +       '  REDUCED MASS = ',F7.1,/)
1028  FORMAT(//,'# BOUND-STATE ENERGIES FOUND FOR L = ',I2,/,
     +          '#      N       E (a.u)        E (eV)')
1029  FORMAT(I8,5X,F9.5,5X,F9.3)
2019  FORMAT(/,'WORST OVERLAP:',1PD12.4)
************************************************************************
*
      znuc = zas
      itout = 11
      ibug = 1
      ipot = 1
      rdmass = 1.0d0
      nix = njdouble-1
      hint = hmax*0.5d0**(nix-1)
      do 10 i=1,nix
       ihx(i) = 2**(i-1)
       irx(i) = jdouble(i+1)
 10   continue
      NSTEP = IRX(NIX)+1
*
*  CALCULATE ARRAYS FOR MESH
*
      RHALF(1)=ZERO
      RFULL(1)=ZERO
      DO 20 I=1,NIX
       HSTEP = HINT*IHX(I)
       HHALF=HALF*HSTEP
       IF (I.EQ.1) THEN
        JBEG = 1
       ELSE
        JBEG = IRX(I-1)+1
       ENDIF
       DO 20 J=JBEG,IRX(I)
        RFULL(J+1) = RFULL(J)+HSTEP
        I2=J+J
        RHALF(I2) = RFULL(J)+HHALF
        RHALF(I2+1) = RFULL(J+1)
20     CONTINUE
30    CONTINUE
*
*  redefine the last interval based on RA
*
      nstepsave = nstep
      do 35 j=1,nstep
*       print *,j,ra,rfull(j)
       if (ra.lt.rfull(j)) goto 36
 35   continue
      print *,' PROBLEM:  ra = ',ra,
     >        ' is greater than last point = ',rfull(nstep)
      stop
 36   irx(nix) = j/2*2
      nstep = irx(nix)+1
      print *,'mesh reset to  npts,ra = ',nstep,rfull(nstep)
      lastpt = nstep-1
      error = 1.0d-10
      nitmax = 100
      estart = -0.6d0
      elow = estart
*      if (la.lt.0) stop
*
*  PRINT MESH PARAMETERS
*
      WRITE(6,1007) HINT,NIX,NSTEP,RFULL(NSTEP)
      WRITE(6,1008) (IHX(I),I=1,NIX)
      WRITE(6,1009) (IRX(I),I=1,NIX)
*
*  SET POTENTIAL IN ATOMIC UNITS
*
      IF (IPOT.EQ.1) THEN
       CORE(1) = -HUGE
       DO 40 I=2,2*NSTEP-1
        CORE(I) = -ZNUC/RHALF(I) 
40     CONTINUE
      ENDIF
*
***      if (la.eq.-1) stop
*
*  HERE WE COULD PROVIDE ANY POTENTIAL ON THE RHALF-MESH
*
*
*  START ITERATION FOR BOUND STATES
*
      DO 170 LBOUND=la,la
       EGUESS = ESTART
       EHIGH = TWO
       NMIN = la+1
       NMAX = nbmax
       IBCNT = 0
       DO 140 NBOUND=NMIN,NMAX
        IF (IBUG.GT.0) WRITE(6,1010) NBOUND,LBOUND
        NNODE = NBOUND-LBOUND-1
        NITER = 0
        ELOW  = ESTART
        IF (NBOUND.GT.NMIN) ELOW = EBOUND(NBOUND-1)
*
*  SET STARTING VALUES FOR RUNGE-KUTTA INTEGRATION.
*  FUNCTION: 0.0  DERIVATIVE: 1.0
*
50      VSTART(1) = ZERO
        VSTART(2) = ONE
        NCLASS = 0
        CALL RKDUMB(VSTART,2,NSTEP,NEND,LBOUND,NCLASS)
        IF (IBUG.GT.3) WRITE(6,1003) RFULL(1),Y(1,1)
*
*  COUNT NODES
*
        INODE = 0
        DO 60 I=2,NEND
         IF (IBUG.GT.3.AND.I.LE.100) WRITE(6,1004) RFULL(I),Y(1,I)
         IF (Y(1,I)*Y(1,I-1).LT.ZERO) INODE = INODE+1
60      CONTINUE
*
*  CHECK WHETHER ENERGY IS TOO LOW OR TOO HIGH AND DEFINE
*  A NEW VALUE FOR EGUESS.
*
        IF (INODE.GT.NNODE) THEN
         EHIGH = EGUESS
        ELSE
         ELOW = EGUESS
        ENDIF
        EGUESS = HALF*(EHIGH+ELOW)
        IF (IBUG.GT.1) WRITE(6,1005) INODE,NNODE,NBOUND,LBOUND,       
     +                               ELOW,EHIGH,EGUESS
        IF (ABS(EHIGH-ELOW).GT.ERROR) THEN
         NITER = NITER+1
         IF (NITER.GT.NITMAX) THEN
          WRITE(6,1006) NITER,NITMAX
          STOP
         ENDIF
         GOTO 50
        ELSE
         IBCNT=IBCNT+1
         EBOUND(NBOUND) = EGUESS
*
*  NORMALIZE THE WAVEFUNCTION
*
         DO 70 I=1,NEND
          FVALUE(I) = Y(1,I)*Y(1,I)
70       CONTINUE
         DO 80 I=NEND-1,1,-1
          IF (FVALUE(I).GT.FVALUE(I+1)) GOTO 90
80       CONTINUE
         WRITE(6,1011) NEND
         STOP
90       NINTEG = I
         CALL INTEGR(RFULL,FVALUE,NINTEG,RESULT,ERRINT)
         FACTOR = ONE/DSQRT(RESULT)
         DO 100 I=1,NSTEP
          IF (I.LE.NINTEG) THEN
           WAVES(NBOUND,I) = Y(1,I)*FACTOR
          ELSE
           WAVES(NBOUND,I) = ZERO
          ENDIF
100      CONTINUE
         IF (IBUG.GT.2) THEN
          WRITE (6,1012) NBOUND, LBOUND
          DO 110 I=1,NINTEG+1
           WRITE(6,1004) RFULL(I),WAVES(NBOUND,I)
110       CONTINUE
         ENDIF
*
*  END OF NORMALISATION. SET THE NEXT GUESS IN SOME REASONABLE WAY.
*
         IF (EGUESS.LT.ZERO) THEN
          EGUESS = HALF*EGUESS
         ELSE 
          EGUESS = TWO*EGUESS
         ENDIF
         XKN = SQRT(TWO*ABS(EGUESS))
         XKNP1 = XKN+5.0/RFULL(NSTEP)
         EHIGH = XKNP1**2+ONE
         IF (NBOUND.LT.NMAX) THEN
          WRITE(6,1001) NBOUND,LBOUND,EBOUND(NBOUND),EGUESS
         ELSE
          WRITE(6,1020) NBOUND,LBOUND,EBOUND(NBOUND)
         ENDIF
        ENDIF
*
        IF (IBCNT.GE.NDIM3) GOTO 140
*
*  CALCULATE THE OVERLAP INTEGRALS
*
        OWORST = ZERO
        DO 130 NNN=NMIN,NBOUND
         DO 120 I=1,NDIM2
          FVALUE(I) = WAVES(NBOUND,I)*WAVES(NNN,I)
120      CONTINUE
         CALL INTEGR(RFULL,FVALUE,IRX(NIX)+1,RESULT,ERRINT)
         IF (ABS(RESULT).GT.1.0D-04.AND.ABS(RESULT-ONE).GT.1.0D-05) 
     +       WRITE(6,1002) NNN,NBOUND,LBOUND,RESULT,ERRINT
         OVRLAP(NNN,NBOUND) = RESULT
         IF (NNN.EQ.NBOUND) THEN
          OWORST = MAX(OWORST,ABS(RESULT-ONE))
         ELSE
          OWORST = MAX(OWORST,ABS(RESULT))
         ENDIF
130     CONTINUE
140    CONTINUE
*
*  PRINT THE OVERLAP INTEGRALS FOR THIS VALUE OF L
*
       WRITE(6,1019) LBOUND
       DO 150 NNN=NMIN,NMAX
        WRITE(6,1013) (OVRLAP(N,NNN),N=NMIN,NNN)
150    CONTINUE
       WRITE(6,2019) OWORST
*
*  WRITE THE ORBITALS TO A FILE (ONE UNIT PER L) IF ITOUT.NE.0
*
       IF (ITOUT.NE.0) THEN
        IT = ITOUT+LBOUND
        WRITE(IT,1013) RFULL(IRX(NIX)+1),
     +                 (EBOUND(NNN)*2.0D0,NNN=NMIN,NMAX)
        DO 160 I=1,IRX(NIX)+1
         WRITE(IT,1013) RFULL(I),(WAVES(NNN,I),
     +                            NNN=NMIN,NMAX)
160     CONTINUE
       ENDIF
170   CONTINUE  
*
*  PRINT ALL BOUND STATE ENERGIES FOUND
*
       WRITE(6,1028) LBOUND
       DO 180 N=nmin,nmax
        WRITE(6,1029) N,EBOUND(N),EBOUND(N)*27.21D0
        erydout(n-la) = EBOUND(N)*2.0d0
        do 185 i=2,IRX(NIX)+1
         waveout(i-1,n-la) = waves(n,i)
 185    continue
180    CONTINUE
*
      RETURN
      END
************************************************************************
      SUBROUTINE RKDUMB(VSTART,NVAR,NSTEP,NEND,LBOUND,NCLASS)
************************************************************************
*                                                                      *
*    THIS SUBROUTINE USES THE FOURTH ORDER RUNGE-KUTTA METHOD TO       *
*    PROPAGATE THE SOLUTION OF NVAR COUPLED DIFFERENTIAL EQUATIONS     *
*    KNOWN AT X1 (STORED IN VSTART) TO X2 BY TAKING NSTEP STEPS OF     *
*    EQUAL LENGTH (X2-X1)/NSTEP. THE RESULTS ARE STORED IN THE         *
*    COMMON BLOCK / PATH /.                                            *
*                                                                      *
************************************************************************
*                                                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER (NDIM1=4,NDIM2=  9000)
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,SIX=6.0D0)
*
      COMMON / PATH / RFULL(NDIM2),RHALF(2*NDIM2),Y(NDIM1,NDIM2)
      COMMON / POT  / EGUESS,CORE(2*NDIM2),RDMASS
*
      DIMENSION VSTART(NVAR),V(NDIM1),DV(NDIM1)
*
*  SET THE CLASSICAL TURNING POINT FOR ZERO ANGULAR MOMENTUM.
*  THIS ONLY MAKES SENSE FOR NEGATIVE ENERGIES
      IF (EGUESS.LT.ZERO) THEN
       DO 10 I=2*NSTEP-1,1,-1
        IF (EGUESS.GT.CORE(I)) THEN
         RCLASS = RHALF(I)
         NCLASS = I/2+1
         GOTO 20
        ENDIF
10     CONTINUE
*
*  THE CLASSICAL TURNING POINT IS AT THE ORIGIN. SOMETHING MUST 
*  BE WRONG. STOP. 
*
       WRITE(6,1000)
       STOP
      ELSE
       RCLASS = RFULL(NSTEP)
       NCLASS = NSTEP
      ENDIF
20    CONTINUE
*
*  LOAD STARTING VALUES
*
      DO 30 I=1,NVAR
       V(I) = VSTART(I)
       Y(I,1) = V(I)
30    CONTINUE
      NEND = NSTEP
*
      DO 50 K=1,NSTEP
       INDEX = K+K-1
       X = RHALF(INDEX)
       CALL DERIVS(INDEX,V,DV,LBOUND)
       H = RHALF(INDEX+2)-RHALF(INDEX)
       CALL RK4(V,DV,NVAR,INDEX,H,V,LBOUND)
       DO 40 I=1,NVAR
        Y(I,K+1) = V(I)
40     CONTINUE
*       IF (X.GT.RCLASS.AND.ABS(Y(1,K+1)).GT.ABS(Y(1,K))) THEN
*        NEND = K+1
*        RETURN
*       ENDIF
50    CONTINUE
*      IF (EGUESS.LT.ZERO) WRITE(6,1001)
*
1000  FORMAT('1',//,' PROGRAM STOPS IN SUBROUTINE RKDUMB, SINCE THE',
     +              ' CLASSICAL TURNING POINT COULD NOT BE DEFINED.',//)
1001  FORMAT(//,' WARNING : SOLUTION HAS NOT TURNED AROUND YET !',/)            
      RETURN
      END
************************************************************************
      SUBROUTINE RK4(Y,DYDX,N,INDEX,H,YOUT,LBOUND)
************************************************************************
*                                                                      *
*  THIS SUBROUTINE PROPAGATES THE SOLUTION FOR N VARIABLES Y BY        *
*  BY ONE STEP H USING THE FOURTH-ORDER RUNGE-KUTTA METHOD. DYDX       *
*  ARE THE DERIVATIVES AND Y THE VALUES OF THE FUNCTIONS AT X.         *
*                                                                      *
************************************************************************
*                                                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER (NDIM1=4,NDIM2=  9000)
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,SIX=6.0D0)
*
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NDIM1),DYT(NDIM1),DYM(NDIM1)
*
      HH = H*HALF
      H6 = H/SIX
*
*  FIRST STEP
*
      DO 10 I=1,N
       YT(I) = Y(I)+HH*DYDX(I)
10    CONTINUE
      CALL DERIVS(INDEX+1,YT,DYT,LBOUND)
*
*  SECOND STEP
*
      DO 20 I=1,N
       YT(I) = Y(I)+HH*DYT(I)
20    CONTINUE
      CALL DERIVS(INDEX+1,YT,DYM,LBOUND)
*
*  THIRD STEP
*
      DO 30 I=1,N
       YT(I) = Y(I)+H*DYM(I)
       DYM(I) = DYT(I)+DYM(I)
30    CONTINUE
      CALL DERIVS(INDEX+2,YT,DYT,LBOUND)
*
*  FOURTH STEP
*
      DO 40 I=1,N
       YOUT(I) = Y(I)+H6*(DYDX(I)+DYT(I)+TWO*DYM(I))
40    CONTINUE
      RETURN
      END
************************************************************************
      SUBROUTINE DERIVS(INDEX,YRK,F,LBOUND)
************************************************************************
*                                                                      *
*  THIS SUBROUTINE CALCULATES THE DERIVATIVES IN THE RUNGE-KUTTA       *
*  AT THE MESHPOINT NUMBER GIVEN BY INDEX.                             *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NDIM1=4,NDIM2=  9000)
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,SIX=6.0D0)
      PARAMETER (HUGE=1.0D30,SMALL=1.0D-20)
*
      COMMON / PATH / RFULL(NDIM2),RHALF(2*NDIM2),Y(NDIM1,NDIM2)
      COMMON / POT  / EGUESS,CORE(2*NDIM2),RDMASS
*
      DIMENSION YRK(NDIM1),F(NDIM1)
*
*  CALCULATE THE DERIVATIVES.
*
      X = RHALF(INDEX)
      F(1) = YRK(2)
      IF (X.NE.ZERO) THEN
       F(2) = (TWO*RDMASS*(CORE(INDEX)-EGUESS)
     +        +(LBOUND*(LBOUND+1)/X**2))*YRK(1)
      ELSE
       F(2) = -HUGE*YRK(1)
      ENDIF
      RETURN
      END
************************************************************************
      SUBROUTINE INTEGR(X,Y,N,RESULT,ERROR)
************************************************************************
*                                                                      *
*  THIS SUBROUTINE INTEGRATES ARBITRARILY SPACED DATA                  *
*                                                                      *
*                      INPUT                                           *
*                      -----                                           *
*                                                                      *
*  X .........     VECTOR CONTAINING THE MESHPOINTS                    *
*                  (EITHER IN ASCENDING OR DESCENDING ORDER)           *
*  Y .........     VECTOR CONTAINING THE FUNCTION VALUES               *
*  N .........     NUMBER OF MESHPOINTS   (AT LEAST 4)                 *
*                                                                      *
*                      OUTPUT                                          *
*                      ------                                          *
*                                                                      *
*  RESULT ....     RESULT OF THE INTEGRATION                           *
*  ERROR .....     ESTIMATED ERROR                                     *
*                                                                      *
************************************************************************
*                                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,                    
     +           FIVE=5.0D0,SIX=6.0D0,TEN=10.0D0,TWELVE=12.0D0,                 
     +           SIXTY=60.0D0,HUNTWE=120.0D0)
      PARAMETER (IREAD=5,IWRITE=6)
      DIMENSION X(N),Y(N)
*
      RESULT = ZERO
      ERROR = ZERO
*
*  CHECK THAT WE HAVE ENOUGH POINTS
*
      IF (N.LT.4) THEN
       WRITE(IWRITE,1000) N
       STOP
      ENDIF
*
*  CHECK THAT THE MESHPOINTS AR IN EITHER ASCENDING OR DESCENDING ORDER
*
      H2 = X(2)-X(1)
      DO 10 I=3,N
       H3 = X(I)-X(I-1)
       IF(H2*H3.LE.ZERO) THEN
        WRITE(IWRITE,1001)
        STOP
       ENDIF
10    CONTINUE
*
*  START THE INTEGRATION
*
      D3 = (Y(2)-Y(1))/H2
      H3 = X(3)-X(2)
      D1 = (Y(3)-Y(2))/H3
      H1 = H2+H3
      D2 = (D1-D3)/H1
      H4 = X(4)-X(3)
      R1 = (Y(4)-Y(3))/H4
      R2 = (R1-D1)/(H4+H3)
      H1 = H1+H4
      R3 = (R2-D2)/H1
      RESULT = H2*(Y(1)+H2*(D3/TWO-H2*(D2/SIX-(H2+TWO*H3)*R3/TWELVE)))
      S = -(H2**3)*(H2*(THREE*H2+FIVE*H4)+TEN*H3*H1)/SIXTY
      R4 = ZERO
      NN = N-1
*
*  LOOP OVER POINTS 2 TO N-1
*
      DO 20 I=3,NN
       RESULT = RESULT+H3*((Y(I)+Y(I-1))/TWO-H3*H3*(D2+R2+(H2-H4)*R3) 
     +         /TWELVE)
       C = H3**3*(TWO*H3*H3+FIVE*(H3*(H4+H2)+TWO*H4*H2))/HUNTWE
       ERROR = ERROR+(C+S)*R4
       IF (I.NE.3) THEN
        S = C
       ELSE
        S = S+TWO*C
       ENDIF
       IF (I.EQ.(N-1)) GOTO 30
       H1 = H2
       H2 = H3
       H3 = H4
       D1 = R1
       D2 = R2
       D3 = R3
       H4 = X(I+2)-X(I+1)
       R1 = (Y(I+2)-Y(I+1))/H4
       R4 = H4+H3
       R2 = (R1-D1)/R4
       R4 = R4+H2
       R3 = (R2-D2)/R4
       R4 = R4+H1
       R4 = (R3-D3)/R4
20    CONTINUE
30    CONTINUE
*
*  FINISH INTEGRATION
*
      RESULT = RESULT+H4*(Y(N)-H4*(R1/TWO+H4*(R2/SIX+(TWO*H3+H4)*R3             
     +        /TWELVE)))
      ERROR = ERROR-H4**3*R4*(H4*(THREE*H4+FIVE*H2)+TEN*H3*(H2+H3+H4)) 
     +        /SIXTY+S*R4
*
1000  FORMAT(//,' ERROR IN SUBROUTINE INTEGR . N = ',I2,3X,                     
     +      'PROGRAM STOPS.',//)                                                
1001  FORMAT(//,' ERROR IN SUBROUTINE INTEGR . MESHPOINTS ARE OUT ',            
     +      'OF ORDER. PROGRAM STOPS.',//)                                      
      RETURN
      END

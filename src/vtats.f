C  FILE 'VTAIL.FTN77'
C
      SUBROUTINE RDTAIL(L,LP,MRED,IGET)
      IMPLICIT REAL*8(A-H,K,O-Z)
      IMPLICIT INTEGER*4(I-J,L-N)
      include 'par.for'
c      PARAMETER ( NM=10 )
c     PARAMETER ( NQM=32 , NWM=NQM+1 , NTOTQ=NQM+NM )
c     PARAMETER ( LU=2 , LMAX=2*LU+1 )
      COMMON/KINT/ NQK,KI(NWM,NM),WT(NM,0:NWM,LMAX),ILOW(NM),IHIGH(NM)
     >   ,IHIGHB(NM),KALL(NTOTQ),IOPEN(NM),JOPEN(NM),NTOTK
      COMMON/TAILME/VTL(LMAX,NTOTQ,NTOTQ)
      COMMON/VTAIL/LOCN(120),IFLG(120),IUT,NRCD
      DATA I100/100/
C
C ****  THIS SUBROUTINE CHECKS WHETHER THE SPECIFIED BLOCK OF
C ****  RADIAL MATRIX ELEMENTS HAS ALREADY BEEN CALCULATED AND
C ****  IF THIS IS SO
C ****  DETERMINES THE POSITION OF THE BLOCK ON THE DIRECT
C ****  ACCESS FILE, AND IF IGET > 0
C ****  READS THE APPROPRIATE BLOCK OF MATRIX ELEMENTS
C ****  AND SETS A FLAG (MRED > 0) WHCIH IS PASSED BACK INTO
C ****  SETAIL (MRED IS REDUNDANT IN ATS) SO THAT THIS BLOCK
C ****  OF MATRIX ELEMENTS DOES NOT NEED TO BE RECOMPUTED.
C
      IF(NRCD.EQ.0) RETURN
      IF(L.LT.LP) THEN
         L1 = LP + 1
         L2 = L + 1
      ELSE
         L1 = L + 1
         L2 = LP + 1
      END IF
      LMN = IABS(L-LP)
      LMX = MIN(LMAX-1,L+LP)
      IF(LMN.EQ.0) LMN = 2
      IF(LMN.GT.LMX) RETURN
      DO 30 LM = LMN,LMX,2
         LAM = (LM-LMN)/2 + 1
         LNDX = (L1*I100 + L2)*I100 + LM
         DO 10 I = 1,NRCD
            IF(LOCN(I).NE.LNDX) GO TO 10
            IF(IGET.GT.0) THEN
               IF(L.GE.LP) THEN
                  READ(IUT,REC=I) ((VTL(LAM,I1,I2),I1=1,NTOTK),I2=1,
     >               NTOTK)
               ELSE
                  READ(IUT,REC=I) ((VTL(LAM,I2,I1),I1=1,NTOTK),I2=1,
     >               NTOTK)
               END IF
            END IF
            MRED = 1
            GO TO 20
10       CONTINUE
20       CONTINUE
30    CONTINUE
      RETURN
      END
      SUBROUTINE WRTAIL(L,LP,LMN,LMX)
      IMPLICIT REAL*8(A-H,K,O-Z)
      IMPLICIT INTEGER*4(I-J,L-N)
      include 'par.for'
c      PARAMETER ( NM=10 )
c     PARAMETER ( NQM=32 , NWM=NQM+1 )
c     PARAMETER ( NTOTQ=NQM+NM )
c     PARAMETER ( LU=2 , LMAX=2*LU+1 )
      COMMON/KINT/ NQK,KI(NWM,NM),WT(NM,0:NWM,LMAX),ILOW(NM),IHIGH(NM)
     >   ,IHIGHB(NM),KALL(NTOTQ),IOPEN(NM),JOPEN(NM),NTOTK
      COMMON/TAILME/VTL(LMAX,NTOTQ,NTOTQ)
      COMMON/VTAIL/LOCN(120),IFLG(120),IUT,NRCD
      DATA I100/100/
C
C ****  IF IT HAS BEEN NECESSARY TO COMPUTE A BLOCK OF RADIAL
C ****  MATRIX ELEMENTS THIS SUBROUTINE IS ENTERED AND THE
C ****  MATRIX ELEMENT BLOCK IS WRITTEN TO EITHER
C ****  A NEW RECORD,
C ****  OR AN OLD RECORD CONTAINING MATRIX ELEMENTS NO LONGER
C ****  NECESSARY IS OVERWRITTEN KEEPING THE SIZE OF THE
C ****  DIRECT ACCESS FILE FROM GROWING TOO LARGE.
C
C ****  OVERWRITE OLD RECORD
C
      IF(L.LT.LP) THEN
         L1 = LP + 1
         L2 = L + 1
      ELSE
         L1 = L + 1
         L2 = LP + 1
      END IF
      DO 40 LM = LMN,LMX,2
         LAM = (LM-LMN)/2 + 1
         LNDX = (L1*I100 + L2)*I100 + LM
         IF(NRCD.EQ.0) GO TO 20
         DO 10 I = 1,NRCD
            IF(IFLG(I).EQ.0) GO TO 10
            LOCN(I) = LNDX
            IFLG(I) = 0
            II = I
            GO TO 30
10       CONTINUE
20       CONTINUE
C
C ****  WRITE NEW RECORD
C
         NRCD = NRCD + 1
         IF(NRCD.GT.120) THEN
            WRITE(6,50)NRCD,L1,L2,LM
            STOP
         END IF
         LOCN(NRCD) = LNDX
         IFLG(NRCD) = 0
         II = NRCD
30       CONTINUE
         IF(L.GE.LP) THEN
            WRITE(IUT,REC=II) ((VTL(LAM,I1,I2),I1=1,NTOTK),I2=1,NTOTK)
         ELSE
            WRITE(IUT,REC=II) ((VTL(LAM,I2,I1),I1=1,NTOTK),I2=1,NTOTK)
         END IF
40       CONTINUE
         RETURN
50       FORMAT(1H0,' DIRECT ACCESS FILE TOO BIG',/,I4
     >      ,' TAIL MATRIX ELEMENTS',4X,' L =',I3,' LP =',I3,' LAM =',
     >      I3)
      END
      SUBROUTINE CHTAIL(J,LB,LT)
      IMPLICIT REAL*8(A-H,K,O-Z)
      IMPLICIT INTEGER*4(I-J,L-N)
      include 'par.for'
c      PARAMETER ( NM=10 , NQM=32 , NTOTQ=NQM+NM )
c     PARAMETER ( LU=2 , LMAX=2*LU+1 )
      COMMON/TAILME/VTL(LMAX,NTOTQ,NTOTQ)
      COMMON/VTAIL/LOCN(120),IFLG(120),IUT,NRCD
      DATA I100/100/
C
C ****  THIS SUBROUTINE SCANS THE CONTENTS OF THE LOOK-UP
C ****  TABLE TO DETERMINE WHETHER A PARTICULAR BLOCK WILL
C ****  BE NEEDED FOR THE NEXT J VALUE.
C
      LBMN = MAX0(J+1-LU,0)
      DO 10 I = 1,NRCD
         LNDX = LOCN(I)
         IF(LNDX.EQ.0) GO TO 10
         L = LNDX/(I100*I100) - 1
         LNDX = LNDX - (L+1)*I100*I100
         LP = LNDX/I100 - 1
         IF((L.LT.LBMN).OR.(LP.LT.LBMN)) THEN
            IFLG(I) = 1
         END IF
10    CONTINUE
      RETURN
      END

      SUBROUTINE SETAIL(J,LLB,LLT,INEW,NBND,NDPL,FPOS)
      IMPLICIT REAL*8(A-H,K,O-Z)
      IMPLICIT INTEGER*4(I-J,L-N)
      LOGICAL*4 LGDP
      COMPLEX*16 CLON,CLOFF,ALON,ALOFF,ALCOFF,TLON,TLOFF,ENLON,ENLOFF
      include 'par.for'
c      PARAMETER ( NR=1000 , NM=10 )
c     PARAMETER ( NQM=32 , NWM=NQM+1 )
c     PARAMETER ( NTOTQ=NQM+NM )
c     PARAMETER ( LU=2 , LU1=LU+1 , LMAX=2*LU+1 )
      COMMON/PHYCHN/NCHAN,NOPEN,LTT(NM),IST(NM),IPAR(NM),ET(NM)
      COMMON/CHANL/ ZZ,NORB,LC(NM),EC(NM)
      COMMON/KINT/ NQK,KI(NWM,NM),WT(NM,0:NWM,LMAX),ILOW(NM),IHIGH(NM)
     >   ,IHIGHB(NM),KALL(NTOTQ),IOPEN(NM),JOPEN(NM),NTOTK
      COMMON/DWPH/CLON(LMAX,NM),CLOFF(LMAX,NQM),ALON(LMAX,NM)
     >   ,ALOFF(LMAX,NQM),ALCOFF(LMAX,NQM),TLON(LMAX,NM)
     >   ,TLOFF(LMAX,NQM),ENLON(LMAX,NM),ENLOFF(LMAX,NQM)
      COMMON/WVCUT/ IBGCR(NM),IBGVL(NM),DELTWF,IBGT
      COMMON /RG/ NDW,RDW(NR),WDW(NR)
      COMMON/TAILME/VTL(LMAX,NTOTQ,NTOTQ)
      COMMON/VTAIL/LOCN(120),IFLG(120),IUT,NRCD
      COMMON/GMFN/FGS(34)
C
C ****  THIS SUBROUTINE COMPUTES THE VALUE OF THE INTEGRAL
C ****  OF TWO CONTINUUM FUNCTIONS AND AN INVERSE POWER OF
C ****  R FROM 0 TO INFINITY
C
C ****  INITIALISE THE LOOK-UP TABLE USED
C ****  TO DETERMINE THE POSITION OF THE RADIAL MATRIX
C ****  ELEMENT BLOCKS IN THE DIRECT ACCESS FILE.
C
      IUT = 27
      IF(INEW.GT.0) THEN
         NRCD = 0
c   double the record length
         IRCSZ = 2*(4*NTOTK*NTOTK)
c$$$         INQUIRE(FILE='/tmp/VTAIL.DMP',EXIST=LGDP)
c$$$         IF(LGDP) THEN
c$$$            CLOSE(UNIT=IUT,STATUS='DELETE')
c$$$         ELSE
c$$$            CLOSE(UNIT=IUT)
c$$$         END IF
c$$$         OPEN(UNIT=IUT,FILE='/tmp/VTAIL.DMP',FORM='UNFORMATTED'
c$$$     >      ,ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IRCSZ)
         CLOSE(UNIT=IUT)
         OPEN(UNIT=IUT,FORM='UNFORMATTED',ACCESS='DIRECT',
     >      STATUS='SCRATCH',RECL=IRCSZ)
C
C ****  CLEAR CONTENTS OF LOOK-UP TABLE AND SET FLAGS SO THAT
C ****  RECORDS CAN BE CREATED.
C ****  IF IFLG(I) > 0, THEN THE ITH RECORD MAY BE WRITTEN
C ****  OR OVERWRITTEN.
C
         CALL IPHASE
         CALL PSISET
         DO 10 I = 1,120
            LOCN(I) = 0
            IFLG(I) = 1
10       CONTINUE
      END IF
C
C ****  RFLEX1 AND RFLEX2 ARE THE INFLEXION POINTS OF THE
C ****  CONTINUUM FUNCTIONS.
C
      IF(INEW.EQ.0.AND.LLB.EQ.0) RETURN
      IGET = 0
      DO 60 L1 = LLB,LLT
         DO 50 L2 = LLB,L1
            LMN = IABS(L1-L2)
            LMX = MIN(LMAX-1,L1+L2)
            RLL1 = DBLE(L1*(L1+1))
            RLL2 = DBLE(L2*(L2+1))
            RFLEX1 = FPOS*ZZ + SQRT(ZZ*ZZ+RLL1)
            RFLEX2 = FPOS*ZZ + SQRT(ZZ*ZZ+RLL2)
            IF(LMN.EQ.0) LMN = 2
            IF(LMN.GT.LMX) GO TO 50
            DO 40 LM = LMN,LMX
               LM2 = (LM-LMN)/2 + 1
               IF(MOD(L1+L2+LM,2).NE.0) GO TO 40
               MRED = 0
               CALL RDTAIL(L1,L2,MRED,IGET)
               IF(MRED.GE.1) GO TO 50
C
C ****  JBLK IS SET EQUAL TO ONE IN DWTAIL AND IS USED TO
C ****  DETERMINE WHETHER QUANTITIES INDEPENDENT OF THE
C ****  WAVENUMBERS, K AND KP, NEED TO BE RECALCULATED.
C
               JBLK = 0
               DO 30 I = 1,NTOTK
                  K = KALL(I)
                  IPTOP = NTOTK
                  IF(L1.EQ.L2) IPTOP = I
C
                  IF(I.LE.NOPEN) THEN
                     x = real(ALON(L1-LLB+1,I))
                  ELSE
                     x = real(ALOFF(L1-LLB+1,I-NOPEN))
                  END IF
                  DO 20 IP = 1,IPTOP
                     KP = KALL(IP)
                     IF (KP.LT.0.OR.K.LT.0) THEN
C  At least one of the K's corresponds to a bound state.
                        VREND=0D0
                        GO TO 15
                     END IF
                     RKTOL = 3.0D0*ABS(K-KP)/DBLE(LM+2)
                     IF(IP.LE.NOPEN) THEN
                        xp = real(ALON(L2-LLB+1,IP))
                     ELSE
                        xp = real(ALOFF(L2-LLB+1,IP-NOPEN))
                     END IF
                     IF (NDPL.EQ.0.AND.RKTOL.LT.14.0D0) THEN
                        IF(1D0-min(x,xp).LT.1D-7) THEN
                           VREND = PWTAIL(L1,I,K,L2,IP,KP,LM,LLB,NBND,
     >                        JBLK)
                        ELSE
                           VREND = DWTAIL(L1,I,K,L2,IP,KP,LM,LLB,NBND,
     >                        JBLK)
                        END IF
                     ELSE
                        VREND = 0.0D0
                     END IF
15                   VTL(LM2,I,IP) = VREND
                     IF(L1.EQ.L2) VTL(LM2,IP,I) = VTL(LM2,I,IP)
20                CONTINUE
30             CONTINUE
40          CONTINUE
            CALL WRTAIL(L1,L2,LMN,LMX)
50       CONTINUE
60    CONTINUE
      RETURN
      END

      REAL*8 FUNCTION DWTAIL(L1,IK1,RK1,L2,IK2,RK2,LM,LLB,NBND,JBLK)
      IMPLICIT REAL*8(A-H,K,O-Z)
      COMPLEX*16 CLON,CLOFF,ALON,ALOFF,ALCOFF,TLON,TLOFF,ENLON,ENLOFF
      COMPLEX*16 AL1,AL2
      include 'par.for'
c     PARAMETER ( NR=1000 , NM=10 )
c     PARAMETER ( NQM=32 , NWM=NQM+1 )
c     PARAMETER ( NTOTQ=NQM+NM )
c     PARAMETER ( LU=2 , LU1=LU+1 , LMAX=2*LU+1 )
      COMMON/PHYCHN/NCHAN,NOPEN,LTT(NM),IST(NM),IPAR(NM),ET(NM)
      COMMON/DWPH/CLON(LMAX,NM),CLOFF(LMAX,NQM),ALON(LMAX,NM)
     >   ,ALOFF(LMAX,NQM),ALCOFF(LMAX,NQM),TLON(LMAX,NM)
     >   ,TLOFF(LMAX,NQM),ENLON(LMAX,NM),ENLOFF(LMAX,NQM)
      common/rdwcut/ rcton(nm),rctoff(nqm)
C
C ****  THIS SUBROUTINE RETURNS THE DEFINITE INTEGRAL FROM
C ****  RCT <= RGMAX TO INFINITY OF TWO CONTINUUM FUNCTIONS TIMES AN
C ****  INVERSE POWER OF R.
C
C ****  L1 IS L-VALUE OF THE IK1 QUADRATURE WITH K = RK1
C ****  L2 IS L-VALUE OF THE IK2 QUADRATURE WITH K = RK2
C ****  LM + 1  IS THE INVERSE POWER OF R.
C
      IF(IK1.LE.NOPEN) THEN
         AL1 = ALON(L1-LLB+1,IK1)
         rd1 = rcton(ik1)
      ELSE
         AL1 = ALOFF(L1-LLB+1,IK1-NOPEN)
         rd1 = rctoff(ik1-nopen)
      END IF
      IF(IK2.LE.NOPEN) THEN
         AL2 = ALON(L2-LLB+1,IK2)
         rd2 = rcton(ik2)
      ELSE
         AL2 = ALOFF(L2-LLB+1,IK2-NOPEN)
         rd2 = rctoff(ik2-nopen)
      END IF
      rgm=min(rd1,rd2)
      x=real(al1)
      y=dimag(al1)
      xp=real(al2)
      yp=dimag(al2)
      dwtail=ffgg(l1,rk1,l2,rk2,rgm,lm,x,y,xp,yp)+
     >   fggf(l1,rk1,l2,rk2,rgm,lm,x,y,xp,yp)
      RETURN
      END

      SUBROUTINE IPHASE
      include 'par.for'
      parameter(lll=lcoul*2+20)
      COMPLEX*16 CIPH
      COMMON/CPHAS/CIPH(-lll:lll)
C
C ****  CIPH(L) CONTAINS I**L
C
      CIPH(0) = DCMPLX(1.0D0,0.0D0)
      CIPH(1) = DCMPLX(0.0D0,1.0D0)
      DO 10 I = 2,lll
         CIPH(I) = CIPH(I-1)*CIPH(1)
10    CONTINUE
C
      DO 20 I = 1,lll
         CIPH(-I) = CIPH(-I+1)/CIPH(1)
20    CONTINUE
      RETURN
      END

      REAL*8 FUNCTION PWTAIL(L1,IK1,RK1,L2,IK2,RK2,LM,LLB,NBND,JBLK)
      IMPLICIT REAL*8(A-H,K,O-Z)
      IMPLICIT INTEGER*4(I-J,L-N)
      REAL UON,UOFF
      include 'par.for'
c      PARAMETER ( NR=1000 , NM=10 )
c     PARAMETER ( NQM=32 , NWM=NQM+1 , NTOTQ=NQM+NM )
c     PARAMETER ( LU=2 , LU1=LU+1 , LMAX=2*LU+1 )
      COMMON/GMFN/FGS(34)
      COMMON/KINT/ NQK,KI(NWM,NM),WT(NM,0:NWM,LMAX),ILOW(NM),IHIGH(NM)
     >   ,IHIGHB(NM),KALL(NTOTQ),IOPEN(NM),JOPEN(NM),NTOTK
      COMMON/PHYCHN/NCHAN,NOPEN,LTT(NM),IST(NM),IPAR(NM),ET(NM)
      COMMON/DW0/UON(LMAX,NM,NR)
      COMMON/DW1/UOFF(LMAX,NQM,NR)
      COMMON /WT/ FACT(0:20),GAM(LMAX),HART,PI,XMU,XMU2
      COMMON/MINI/IMON(LMAX,NM),IMOFF(LMAX,NQM)
      COMMON/DWCUT/ICTON(NM),ICTOFF(NQM)
      COMMON /RG/ NDW,RDW(NR),WDW(NR)
      COMMON/GRD2/ RGMAX,RGMX2,RGMIN,ZN,AE,ALP2,BE2,HB2,DE2
      COMMON/CM1/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,TWOL
      DIMENSION U1(NR),U2(NR),RLM(NR)
      SAVE RLM,U1,IKOLD,I1,I2
C
C ****  THIS SUBROUTINE COMPUTES THE VALUE OF THE INTEGRAL OF
C ****  TWO SPHERICAL BESSEL FUNCTIONS AND AN INVERSE POWER OF
C ****  R FROM A FINITE VALUE, MOSTLY RGMAX, TO INFINITY. THE ANALYTIC FORMULA
C ****  EQ. (6.574.2) OF GRADSTEYN AND RYZHIK IS USED TO WRITE THE
C ****  INTEGRAL AS A HYPERGEOMETRIC SERIES WHICH CONVERGES
C ****  MONOTONICALLY. WHEN THE RATIO RK1/RK2 IS CLOSE TO ONE THE
C ****  HYPERGEOMETRIC CONVERGES SLOWLY AND A KUMMER TRANSFORMATION,
C ****  EQ.(15.3.11) IS USED TO TRANSFORM THE HYPERGEOMETRIC INTO A
C ****  A FORM THAT CONVERGES MORE QUICKLY FOR SMALL L1 AND L2 VALUES.
C
C ****  EVALUATION OF INVARIANT CONSTANTS NEEDED IN PWTAIL
C
      IF(JBLK.EQ.0) THEN
         IZR = 0
         IONE = 1
         IKOLD = 0
         RT10 = SQRT(10.0D0)
         T1 = GAMX(IONE,2*LM)/(GAMX(IONE,L1-L2+LM+1)
     >      *GAMX(IONE,L2-L1+LM+1))
         T2 = GAMX(IZR,L1+L2-LM+2)/(GAMX(IZR,L1+L2+LM+2)*(1.0D1**LM))
         S3 = RT10**(L2-L1-LM-1)
         T3 = (GAMX(IZR,L1+L2+2-LM)*S3)/GAMX(IZR,2*L1+3)
         S4 = RT10**(L1-L2-LM-1)
         T4 = (GAMX(IZR,L1+L2+2-LM)*S4)/GAMX(IZR,2*L2+3)
         T5 = GAMX(IONE,L2-L1+LM+1)
         T6 = GAMX(IONE,L1-L2+LM+1)
         T7 = GAMX(IONE,2*LM)
         S8 = RT10**(L1-L2-LM+1)
         T8 = S8*GAMX(IZR,2*L1+3)/GAMX(IZR,L1+L2+LM+2)
         T8 = T8/GAMX(IONE,L1-L2+LM+1)
         S9 = RT10**(L2-L1-LM+1)
         T9 = S9*GAMX(IZR,2*L2+3)/GAMX(IZR,L1+L2+LM+2)
         T9 = T9/GAMX(IONE,L2-L1+LM+1)
         T10 = GAMX(IONE,L1-L2-LM+1)
         T11 = GAMX(IONE,L2-L1-LM+1)
         TWOL = DBLE(2**LM)
         DO 10 I = 1,NDW
            RLM(I) = WDW(I)/(RDW(I)**(LM+1))
10       CONTINUE
      END IF
      JBLK = 1
      PWTAIL = 0.0D0
      YINF = 0.0D0
      XRT = MAX(0.85D0,1.0D0-1.0D0/DBLE(L1+L2+1))
C
C ****  THIS SECTION EVALUATES THE DEFINITE INTEGRAL FROM
C ****  ZERO TO INFINITY. TWO DIFFERENT REPRESENTATIONS OF
C ****  THE HYPERGEOMETRIC FUNCTION ARE USED TO SPEED UP
C ****  THE CALCULATIONS. THE HYPERGEOMETRIC SERIES IS
C ****  NOTORIOUSLY SLOW TO CONVERGE WHEN THE ARGUMENT IS
C ****  CLOSE TO ONE.
C
C ****  FIRST SPECIAL CASE, RK1 = RK2
C
      TK = (RK1/RK2) - 1.0D0
      IF((IK1.EQ.IK2).OR.(ABS(TK).LT.1.0D-8)) THEN
         YINF = ((RK1*0.5D0)**LM)*T1*T2
C
C ****  NEXT CASES, RK1 < RK2  OR RK1 > RK2
C
      ELSE
C
C ****  RK1 < RK2
C
         WA = 0.5D0*DBLE(L1+L2+2-LM)
         IF(RK1.LT.RK2) THEN
            TK = RK1/RK2
            ARGF = TK*TK
            TK11 = (TK**(L1+1))
            VK = RK2**LM
            IF(ARGF.LE.XRT) THEN
               WB = 0.5D0*DBLE(L1-L2-LM+1)
               WC = 0.5D0*DBLE(2*L1+3)
               NTERM = -1
               F21 = FDHY(WA,WB,WC,NTERM,ARGF)
               YINF = F21*VK*TK11*T3/(T5*TWOL)
            ELSE
               ARGF = 1.0D0 - TK*TK
               WL = (-ARGF)**LM
               NTERM = LM - 1
               UB = 0.5D0*DBLE(L1-L2-LM+1)
               UC = DBLE(1-LM)
               F2X = T7*T8*FDHY(WA,UB,UC,NTERM,ARGF)
               VA = 0.5D0*DBLE(L1+L2+LM+2)
               VB = 0.5D0*DBLE(L1-L2+LM+1)
               VC = DBLE(LM)
               F2Y = WL*FDHY2(VA,VB,VC,ARGF)/(T3*T10)
               YINF = TK11*VK*T3*(F2X - F2Y)/(TWOL*T5)
            END IF
C
C ****  RK2 < RK1
C
         ELSE
            TK = RK2/RK1
            ARGF = TK*TK
            TK11 = (TK**(L2+1))
            VK = RK1**LM
            IF(ARGF.LE.XRT) THEN
               WB = 0.5D0*DBLE(L2-L1-LM+1)
               WC = 0.5D0*DBLE(2*L2+3)
               NTERM = -1
               F21 = FDHY(WA,WB,WC,NTERM,ARGF)
               YINF = F21*VK*TK11*T4/(T6*TWOL)
            ELSE
               ARGF = 1.0D0 - TK*TK
               WL = (-ARGF)**LM
               NTERM = LM - 1
               UB = 0.5D0*DBLE(L2-L1-LM+1)
               UC = DBLE(1-LM)
               F2X = T7*T9*FDHY(WA,UB,UC,NTERM,ARGF)
               VA = 0.5D0*DBLE(L1+L2+LM+2)
               VB = 0.5D0*DBLE(L2-L1+LM+1)
               VC = DBLE(LM)
               F2Y = WL*FDHY2(VA,VB,VC,ARGF)/(T4*T11)
               YINF = TK11*VK*T4*(F2X - F2Y)/(TWOL*T6)
            END IF
         END IF
      END IF
      YINF = PI*YINF/2.0D0
C
C ****  THIS NEXT SECTION EVALUATES THE INTEGRAL FROM
C ****  ZERO TO RGMAX. Actually, for large RK1 or RK2 the integral is
C ****  cut off before RGMAX. The cutoff point RDW(I) is determined
C ****  by ICTOFF(IK1) or ICTOFF(IK2). The cutoff is effected by setting
C ****  U1(I) or U2(I) to zero for I>ICTOFF.
C
      IF(IK1.EQ.IKOLD) GO TO 40
      IF(IK1.LE.NOPEN) THEN
         MCH1 = IOPEN(IK1)
         I1 = IMON(L1-LLB+1,MCH1)
         I2 = ICTON(MCH1)
         DO 20 I = 1,NDW
            U1(I) = UON(L1-LLB+1,MCH1,I)
20       CONTINUE
      ELSE
         I1 = IMOFF(L1-LLB+1,IK1-NOPEN)
         I2 = ICTOFF(IK1-NOPEN)
c         DO 30 I = 1,NDW
         DO 30 I = i1,i2
            U1(I) = UOFF(L1-LLB+1,IK1-NOPEN,I)
30       CONTINUE
      END IF
40    CONTINUE
      IF(IK2.LE.NOPEN) THEN
         MCH2 = IOPEN(IK2)
         I3 = IMON(L2-LLB+1,MCH2)
         I4 = ICTON(MCH2)
         DO 50 I = 1,NDW
            U2(I) = UON(L2-LLB+1,MCH2,I)
50       CONTINUE
      ELSE
         I3 = IMOFF(L2-LLB+1,IK2-NOPEN)
         I4 = ICTOFF(IK2-NOPEN)
c         DO 60 I = 1,NDW
         DO 60 I = i3,i4
            U2(I) = UOFF(L2-LLB+1,IK2-NOPEN,I)
60       CONTINUE
      END IF
      IKOLD = IK1
c      I5 = MIN(I1,I3)
      I5 = max(I1,I3)
      i6=min(i2,i4)
      Y0R = 0.0D0
c      DO 70 I = I5,NDW
      DO 70 I = I5,i6
         Y0R = Y0R + RLM(I)*U1(I)*U2(I)
70    CONTINUE
      PWTAIL = YINF - Y0R
      RETURN
      END

      DOUBLE PRECISION FUNCTION GAMX(ITEN,I)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      include 'par.for'
      parameter(lll=lcoul*2+20)
      COMMON/GAMMA/GS(lll)
C
C ****  THIS SUBROUTINE RETURNS THE VALUE OF THE GAMMA
C ****  FUNCTION FOR X = I*0.5D0
C ****  IF ITEN = 0,  RETURN GAMX(X) = GAMMA(X)/(10**X)
C ****  IF ITEN > 0,  RETURN GAMX(X) = GAMMA(X)
C ****  IF I < 0 AND EVEN, NO VALUE IS RETURNED SINCE THE GAMMA
C ****  FUNCTION HAS POLES AT NEGATIVE INTEGERS
C
      II = I
      IF(II.LE.0) THEN
         IF(MOD(II,2).EQ.0) THEN
            WRITE(6,20) II
            GAMX = 0.0D0
            RETURN
         END IF
         II = 1
      END IF
C
      IF(ITEN.EQ.0) THEN
         WGX = 1.0D0
      ELSE
         WGX = SQRT(10.0D0)
      END IF
C
      GAMX = GS(II)*(WGX**II)
C
C ****  RECURSION FORMULA ARE USED TO EXTEND THE LOOK-UP TABLE
C ****  FOR THE NEGATIVE 1/2 INTEGERS
C
      IF(I.NE.II) THEN
         IF(ITEN.EQ.0) THEN
            WGX = 10.0D0
         ELSE
            WGX = 1.0D0
         END IF
         II2 = II - 2
         DO 10 J = II2,I,-2
            GAMX = GAMX*WGX/(0.5D0*DBLE(J))
10       CONTINUE
      END IF
      RETURN
20    FORMAT(1X,' GAMMA FUNCTION UNDEFINED FOR X =',I3,'/2'
     >   ,'  GAMX = 0.0 RETURNED')
      END

      SUBROUTINE GAMSET
      IMPLICIT REAL*8(A-H,O-Z)
      include 'par.for'
      parameter(lll=lcoul*2+20)
      COMMON/GAMMA/GS(lll)
C
C ****  COMPUTES VALUES OF THE GAMMA FUNCTION TIMES A SCALE FACTOR.
C ****  GS(I)  CONTAINS GAMMA(I/2)/(10**(I/2))
C
      GS(2) = 1.0D0/10.0D0
      DO 10 I = 4,lll,2
         GS(I) = 0.5D0*DBLE(I-2)*GS(I-2)/10.0D0
10    CONTINUE
      GS(1) = 1.7724538509055D0/SQRT(10.0D0)
      DO 20 I = 3,lll,2
         GS(I) = 0.5D0*DBLE(I-2)*GS(I-2)/10.0D0
20    CONTINUE
      RETURN
      END

      SUBROUTINE PSISET
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PSVAL/PSI(600)
C
C ****  COMPUTES THE VALUE OF THE PSI(DIGAMMA) FUNCTION.
C ****  DIGAMMA(I/2) IS STORED IN PSI(I)
C
      PSI(1) = -1.96351002602143D0
      PSI(2) = -0.577215664901533D0
      DO 10 I = 3,600
         PSI(I) = PSI(I-2) + 2.0D0/DBLE(I-2)
10    CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FDHY(A,B,C,NTERM,X)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
C
C ****  CALCULATES THE HYPERGEOMETRIC POWER SERIES TO AN ACCURACY
C ****  OF 10**(-7.0). A SIMPLE POWER SERIES ALGORITHM IS USED TO
C ****  DETERMINE WHERE TO TRUNCATE THE SUMMATION WHEN THE SERIES
C ****  IS INFINITE. WHEN NTERM IS A POSITIVE INTEGER THE SERIES
C ****  TERMINATES AND A DIFFERENT SECTION OF CODE IS USED.
C
      ABCF = 1.0D0
      XI = 1.0D0
      S1 = 0.0D0
      s2=1d0
      i=0
      IF(NTERM.LT.0) THEN
c         while(abs(s2-s1).gt.1d-7) do
 9999     if(abs(s2-s1).gt.1d-7)then
            i=i+1
            s2=s1
            S1 = S1 + ABCF*XI
            SI = DBLE(I-1)
            ABCF = ABCF*(A+SI)*(B+SI)/((C+SI)*(SI+1.0D0))
            XI = XI*X
            go to 9999
           endif
c         end while
      ELSE
         NN = NTERM + 1
         DO 20 I = 1,NN
            S1 = S1 + ABCF*XI
            IF(I.EQ.NN) GO TO 30
            SI = DBLE(I-1)
            ABCF = ABCF*(A+SI)*(B+SI)/((C+SI)*(SI+1.0D0))
            XI = XI*X
20       CONTINUE
      END IF
30    CONTINUE
      FDHY = S1
      RETURN
      END

      DOUBLE PRECISION FUNCTION FDHY2(A,B,C,X)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/GMFN/FGS(34)
      COMMON/PSVAL/PSI(600)
C
C ****  THIS FUNCTION SUMS THE INFINITE PART OF THE ALTERNATE
C ****  FORM OF THE HYPERGEOMETRIC FUNCTION. THE ACCURACY TOLERANCE
C ****  IS 10**(-7.0D0)
C
      XLN = LOG(X)
      IA = INT(2.0D0*A+0.1D0)
      IB = INT(2.0D0*B+0.1D0)
      IC = INT(2.0D0*C+0.1D0)
      IM = INT(C+0.1D0)
      ABCF = 1.0D0/FGS(IM+1)
      XI = 1.0D0
      S1 = 0.0D0
      i=0
      s2=1d0
c      while(abs(s2-s1).gt.1d-7) do
 9997  if(abs(s2-s1).gt.1d-7)then
         i=i+1
         s2=s1
         I2 = I*2
         APSI = XLN - PSI(I2) - PSI(I2+IC) + PSI(I2-2+IA) + PSI(I2-2+IB)
         S1 = S1 + APSI*ABCF*XI
         SI = DBLE(I)
         SI1 = DBLE(I-1)
         ABCF = ABCF*(A+SI1)*(B+SI1)/((C+SI)*SI)
         XI = XI*X
c      end while
       go to 9997
       endif
      FDHY2 = S1
      IF(i+4.GT.(600/2)) THEN
         WRITE(6,20)A,B,C,X,i
      END IF
      RETURN
20    FORMAT(1H ,' PROBLEM IN FDHY2. WITH A=',G15.6,' B=',G15.6
     >   ,' C=',G15.6,' X=',G15.6,/
     >   ,' NN =',I5,' IS TOO LARGE')
      END

C  The following function returns x*xp*FF+y*yp*GG, where FF is
C  the integral from r=a to r=infinity of F(l1,r1*r)*F(l2,r2*r)/r**(lambda+1),
C  and where GG is the integral from r=a to r=infinity of
C  G(l1,r1*r)*G(l2,r2*r)/r**(lambda+1). F(l,r)=r*j(l,r), where j(l,r)
C  is the spherical Bessel function. G(l,r)=r*n(l,r), where n(l,r) is the
C  Neumann spherical function.
      function ffgg(l1,r1,l2,r2,a,lambda,x,y,xp,yp)
      implicit real*8 (a-h,o-z)
      logical*4 cosl
      include 'par.for'
      dimension c(0:2*lcoul),jstart(2),jstop(2),jstep(2)
      sum=0d0
      aer=1d-5
      jm=nint(-(r1*a+r2*a+1)+sqrt((2*r1*a-1)**2+4*l1*(l1+1))/2
     >   +sqrt((2*r2*a-1)**2+4*l2*(l2+1))/2)
      ltot=l1+l2
      if (jm.lt.0) jm=0
      if (jm.gt.ltot) jm=ltot
      jstart(1)=jm
      jstop(1)=0
      jstep(1)=-1
      jstart(2)=jm+1
      jstop(2)=ltot
      jstep(2)=1
      do 50 i=1,2
         do 10 j=jstart(i),jstop(i),jstep(i)
            n=j+ltot
            cosl=mod(n,2).eq.0
            mn=min(l1,j)
            mx=max(0,j-l2)
            tmp=0d0
            do 20 m=mx,mn
               c(m)=cnst(l1,m,r1,l2,j-m,r2)
               tmp=tmp+c(m)
20          continue
            app=tmp/(a**(lambda+j)*(lambda+j))
            rer=aer/app
            if (rer.gt.1d0) go to 50
            if (rer.lt.1d-12) then
               ffgg=0d0
c               print*,'Unable to handle this case of FF+GG integral '
c               print*,l1,'=l1',l2,'=l2',sngl(r1),'=k1',sngl(r2),'=k2'
c               print*,lambda,'=lambda. Tail set to zero.'
               return
            end if
            tmp1=csint(a,r1+r2,lambda+1+j,cosl,rer)
            tmp2=csint(a,r2-r1,lambda+1+j,cosl,rer)
            tmp=tmp*tmp1*(y*yp-x*xp)
            do 30 m=mx,mn
               tmp=tmp+c(m)*(-1)**(m+l1)*tmp2*(x*xp+y*yp)
30          continue
            if (cosl) then
               sum=sum+tmp*(-1)**(n/2)
            else
               sum=sum+tmp*(-1)**((n+1)/2)
            end if
10       continue
50    continue
      ffgg=(-1)**(l1+l2)*sum/2d0
      end

C  The following function returns x*yp*FG+y*xp*GF, where FG is
C  the integral from r=a to r=infinity of F(l1,r1*r)*G(l2,r2*r)/r**(lambda+1),
C  and where GF is the integral from r=a to r=infinity of
C  G(l1,r1*r)*F(l2,r2*r)/r**(lambda+1). F(l,r)=r*j(l,r), where j(l,r)
C  is the spherical Bessel function. G(l,r)=r*n(l,r), where n(l,r) is the
C  Neumann spherical function.
      function fggf(l1,r1,l2,r2,a,lambda,x,y,xp,yp)
      implicit real*8 (a-h,o-z)
      logical*4 cosl
      include 'par.for'
      dimension c(0:2*lcoul),jstart(2),jstop(2),jstep(2)
      sum=0d0
      aer=1d-5
      jm=nint(-(r1*a+r2*a+1)+sqrt((2*r1*a-1)**2+4*l1*(l1+1))/2
     >   +sqrt((2*r2*a-1)**2+4*l2*(l2+1))/2)
      ltot=l1+l2
      if (jm.lt.0) jm=0
      if (jm.gt.ltot) jm=ltot
      jstart(1)=jm
      jstop(1)=0
      jstep(1)=-1
      jstart(2)=jm+1
      jstop(2)=ltot
      jstep(2)=1
      yypmx=max(abs(y),abs(yp))
      do 50 i=1,2
         do 10 j=jstart(i),jstop(i),jstep(i)
            n=j+ltot
            cosl=mod(n,2).ne.0
            mn=min(l1,j)
            mx=max(0,j-l2)
            tmp=0d0
            do 20 m=mx,mn
               c(m)=cnst(l1,m,r1,l2,j-m,r2)
               tmp=tmp+c(m)
20          continue
            app=yypmx*tmp/(a**(lambda+j)*(lambda+j))
            rer=aer/app
            if (rer.gt.1d0) go to 50
            if (rer.lt.1d-12) then
               fggf=0d0
c               print*,'Unable to handle this case of FG+GF integral '
c               print*,l1,'=l1',l2,'=l2',sngl(r1),'=k1',sngl(r2),'=k2'
c               print*,lambda,'=lambda. Tail set to zero.'
               return
            end if
            tmp1=csint(a,r1+r2,lambda+1+j,cosl,rer)
            tmp2=csint(a,r2-r1,lambda+1+j,cosl,rer)
            tmp=tmp*tmp1*(x*yp+xp*y)
            do 30 m=mx,mn
               tmp=tmp-c(m)*(-1)**(m+l1)*tmp2*(x*yp-xp*y)
30          continue
            if (.not.cosl) then
               sum=sum-tmp*(-1)**(n/2)
            else
               sum=sum+tmp*(-1)**((n+1)/2)
            end if
10       continue
50    continue
      fggf=(-1)**(l1+l2+1)*sum/2d0
      end

      function cnst(l1,i1,r1,l2,i2,r2)
      implicit real*8 (a-h,o-z)
      include 'par.for'
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),res(63)
      cnst=exp(factl(l1+i1)+factl(l2+i2)-factl(i1)-factl(i2)
     >   -factl(l1-i1)-factl(l2-i2)-i1*dlog(2*r1)-i2*dlog(2*r2))
      end

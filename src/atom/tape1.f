      SUBROUTINE TAPE(R5,XI,PHS,IMA,IJ,NE,IQ,nn5)
C  áAèàCú (IJ=7)  àãà  CóàTõBAHàE  (IJ=2)
C   BOãHOBõX  îìHKñàâ
      REAL*8    R5(1004),PHS,XI2,XI
      IRA1=3
      NZ1=IMA
      KZ=1
      IRA2=NE+4
      DO 5 KI=1,KZ
      XI2=DFLOAT(IQ)
      IF(IJ.EQ.2) GOTO 3
      WRITE(1 )    R5(1),XI,R5(2),XI2,PHS,
     *  (R5(II),II=IRA1,IRA2)
      GOTO 5
      READ(2 ,   END=5) R5(1),XI,R5(2),XI2,PHS
3     continue
      do 55 ii=ira1,ira2
      read(2,  end=5)  R5(II)
55    continue
      IQ=XI2+0.5
    5 CONTINUE
      RETURN
      END
      SUBROUTINE TAPE1(R5,XI,PHS,IMA,IJ,NE,IQ,nn5)
C  áAèàCú (IJ=7)  àãà  CóàTõBAHàE  (IJ=2)
C   BOãHOBõX  îìHKñàâ
      REAL*8    R5(1004),PHS,XI2,XI
      IRA1=3
      NZ1=IMA
      KZ=1
      IRA2=NE+4
      DO 5 KI=1,KZ
      XI2=DFLOAT(IQ)
      IF(IJ.EQ.2) GOTO 3
      WRITE(1 )    R5(1),XI,R5(2),XI2,PHS,
     *  (R5(II),II=IRA1,IRA2)
      GOTO 5
      READ(3 ,   END=5) R5(1),XI,R5(2),XI2,PHS
3     continue
      do 55 ii=ira1,ira2
      read(3,  end=5)  R5(II)
55    continue
      IQ=XI2+0.5
    5 CONTINUE
      RETURN
C     DEBUG INIT
      END
      SUBROUTINE TAPE2(R5,XI,PHS,IMA,IJ,NE,IQ,nn5)
C  áAèàCú (IJ=7)  àãà  CóàTõBAHàE  (IJ=2)
C   BOãHOBõX  îìHKñàâ
      REAL*8    R5(1004),PHS,XI2,XI
      IRA1=3
      NZ1=IMA
      KZ=1
      IRA2=NE+4
      DO 5 KI=1,KZ
      XI2=DFLOAT(IQ)
      IF(IJ.EQ.2) GOTO 3
      WRITE(1 )    R5(1),XI,R5(2),XI2,PHS,
     *  (R5(II),II=IRA1,IRA2)
      GOTO 5
    3 READ( 4,   END=5) R5(1),XI,R5(2),XI2,PHS,
     *(R5(II),II=IRA1,IRA2)
      IQ=XI2+0.5
    5 CONTINUE
      RETURN
      END
      SUBROUTINE TAPE3(R5,XI,PHS,IMA,IJ,NE,IQ,nn5)
C  áAèàCú (IJ=7)  àãà  CóàTõBAHàE  (IJ=2)
C   BOãHOBõX  îìHKñàâ
      REAL*8    R5(1004),PHS,XI2,XI
      IRA1=3
      NZ1=IMA
      KZ=1
      IRA2=NE+4
      DO 5 KI=1,KZ
      XI2=DFLOAT(IQ)
      IF(IJ.EQ.2) GOTO 3
      WRITE(1 )    R5(1),XI,R5(2),XI2,PHS,
     *  (R5(II),II=IRA1,IRA2)
      GOTO 5
    3 READ(10,END=5) R5(1),XI,R5(2),XI2,PHS,
     *(R5(II),II=IRA1,IRA2)
      IQ=XI2+0.5
    5 CONTINUE
      RETURN
      END
      SUBROUTINE TAPE4(R5,XI,PHS,IMA,IJ,NE,IQ,nn5)
C  áAèàCú (IJ=7)  àãà  CóàTõBAHàE  (IJ=2)
C   BOãHOBõX  îìHKñàâ
      REAL*8    R5(1004),PHS,XI2,XI
      IRA1=3
      NZ1=IMA
      KZ=1
      IRA2=NE+4
      DO 5 KI=1,KZ
      XI2=DFLOAT(IQ)
      IF(IJ.EQ.2) GOTO 3
      WRITE(1 )    R5(1),XI,R5(2),XI2,PHS,
     *  (R5(II),II=IRA1,IRA2)
      GOTO 5
    3 READ( 8,   END=5) R5(1),XI,R5(2),XI2,PHS,
     *(R5(II),II=IRA1,IRA2)
      IQ=XI2+0.5
    5 CONTINUE
      RETURN
C     DEBUG INIT
      END
      SUBROUTINE TAPE5(R5,XI,PHS,IMA,IJ,NE,IQ,nn5)
C  áAèàCú (IJ=7)  àãà  CóàTõBAHàE  (IJ=2)
C   BOãHOBõX  îìHKñàâ
      REAL*8    R5(1004),PHS,XI2,XI
      IRA1=3
      NZ1=IMA
      KZ=1
      IRA2=NE+4
      DO 5 KI=1,KZ
      XI2=DFLOAT(IQ)
      IF(IJ.EQ.2) GOTO 3
      WRITE(1 )    R5(1),XI,R5(2),XI2,PHS,
     *  (R5(II),II=IRA1,IRA2)
      GOTO 5
    3 READ( 9,   END=5) R5(1),XI,R5(2),XI2,PHS,
     *(R5(II),II=IRA1,IRA2)
      IQ=XI2+0.5
    5 CONTINUE
      RETURN
      END

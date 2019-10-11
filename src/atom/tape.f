      SUBROUTINE TAPE(R5,XI,PHS,IJ,NE,IQ,nf)
      include 'paratom.f'
      REAL*8    R5(nsize),PHS,XI2,XI
      IRA1=3
      IRA2=NE+4
      XI2=DFLOAT(IQ)
      IF(IJ.EQ.2) GOTO 3
      WRITE(nf)    R5(1),XI,R5(2),XI2,PHS
      do 56 ii=ira1,ira2
      write(nf) r5(ii)
56    continue
      return
3     continue
      READ( nf ,   END=5) R5(1),XI,R5(2),XI2,PHS
      do 55 ii=ira1,ira2
      read(nf ,  end=5)  R5(II)
55    continue
      IQ=XI2+0.5
5     continue
      RETURN
      END

      SUBROUTINE JNT(
     >   R5,R8,RO1,BET,ALFA,R7,H,R2,RM,IS,NE,NE5,NE5IS,KT,IN,K)
!      REAL*8 ALFA,BET,H,Z,RO1,RO2,EPS,R2(KT ),
      REAL*8 ALFA,BET,H,Z,RO1,RO2,EPS,R2(KT+5 ),
     * R5(NE5 ),R7( KT),R8(NE5IS),ROT,
     * AA1,T2,AKSI,UU,VV ,R,RM
      INTEGER IN(IS)


      do iww = 1, kt + 5
         r2(iww) = 0.0
      enddo
      
      DO 1 JJ=1,IS
      IAA=JJ
      
* Write wave functions to file

c$$$      open(10,file='wf')
*      write(*,*) ' k=',k,' jj=',jj
*      write(*,*) ' kt=',kt
      IF(K.EQ.3.AND.JJ.NE.IS)  GOTO 1
      MSK= NE5 *  (IAA-1)
      DO 10 IW=1,NE5
   10 R5(IW)=R8(IW+MSK)
*      write(6,   2) R5(1),R5(NE+3)
    2 FORMAT(  ' ','Energy=',25X,E15.8,
     */' ','HOPMA  =' ,25X,E15.8)
*      IF(R5(1).LT.0.0)  write(6, 201) IN(JJ)
201   FORMAT(' Principal quantum number=',6X,I5)
*      write(6, 202) R5(2)
202   FORMAT(' Orbital moment=     ',
     *  9X,     E15.8// )
*      write(6, 5) JJ
    5 FORMAT( /  '   Wave function',I3,' shell'
     */3('  R',11X,'P(R)',18X))

      DO 3 IWW=1,KT
      ROT=ALFA*R7(IWW)+BET* LOG(R7(IWW))
      AA1=(ROT-RO1)/H
      IT1=AA1
      T2=AA1-IT1
      IF(IT1.EQ.0) GOTO 3
C
C Not too close to RM or NE
C
      IF(IT1.GT.NE .or. abs(R7(IWW)-RM).lt.1D-1) GOTO 1
      AKSI=1.-T2
      UU=T2*(1.-T2**2)/6.
      VV=AKSI*(1.-AKSI**2)/6.
      R2(IWW)=SQRT(R7(IWW)/(ALFA*
     *R7(IWW)+BET))*(R5(IT1+2)*AKSI
     *+R5(IT1+3)*T2-VV*(R5(IT1+3)-2.
     * *R5(IT1+2)+R5(IT1+1))-
     *UU*(R5(IT1
     *+4)-2.*R5(IT1+3)+R5(IT1+2)))
c$$$      write(6,'(I5,2x,2F9.4,2x,I7,2x,4F9.4)')
c$$$     >   IWW,R7(IWW),R2(IWW),it1,rot,aa1,uu,vv
    3 CONTINUE
    1 CONTINUE

*      write(6, 4) (R7(IW),R2(IW),IW=1,KT )
*      write(10, 40) (R7(IW),R2(IW),IW=1,KT )
    4 FORMAT(3(2X,F8.5,2X,E20.13,4X))
 40   FORMAT(2X,F9.4,2X,F11.6,4X)
*      write(6, 14)
14    FORMAT(/)
      RETURN
      END


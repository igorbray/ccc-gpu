      REAL*8 A(15),PRE        ,R5(1005),X(3),FUNC(4,1005),           
     * ALFA,BET,H,R,RO1,RO2,Z,R1(1005),R2(1005),R6(1005),ME,ME1,C1,C2
      INTEGER LM(4)
      character*15  fv
c
      pre=1.E-04
c
        open(5,file='culon.dat')
        open(6,file='con')
        open(7,file='culon.lst')
      WRITE(6,1)                                               
      WRITE(7,1)                                               
1     FORMAT(' Calculation of redused coulomb matrix elements'/ 
     *  //5X,'1-----2',/,2(8X,'!'/),5X,'3-----4'/)           
10    CONTINUE                                               
      N=11                                                   
      DO 11 I=1,4                                            
c      WRITE(6,2) I                                          
2     FORMAT('Enter the shell number  for ',I2,' Bf')        
      READ(5,*) J                                            
      read(5,222) fv
222     format(a15)
        open(n,file=fv,form='unformatted')
      IF(J.EQ.-1) STOP                                       
      READ(N) A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8),
     * A(9),A(10),A(11),A(12),A(13),A(14),A(15)
      IF(DABS(A(1)-1.D 00).LT.PRE) IS=IFIX(SNGL(A(9)))       
      IF(DABS(A(1)-2.D 00).LT.PRE) IS=IFIX(SNGL(A(12)+A(13)))
      NE=IFIX(SNGL(A(3)))                         
      NE4=NE+4                                    
      DO 12 L=1,IS                                
      READ(N) R5(1),X(1),R5(2),X(2),X(3)
      do 223 ii=3,ne4
      read(n) r5(ii)
223     continue
      IF(L.EQ.J) GOTO 13                          
 12   CONTINUE                                    
      WRITE(6,20) IS,J                            
 20   FORMAT(' *** ERROR ****   value of  Bf=',I4,
     * '  is less then entered sell number =',I4)
      GOTO 10                          
 13    CONTINUE                        
      NE5=NE+5                         
      LM(I)=IFIX(SNGL(R5(2)))          
      DO 14 II=3,NE5                   
 14   FUNC(I,II)=R5(II)                
      WRITE(6,3) I,R5(1),LM(I)         
 3    FORMAT(1X,I2,'   Function  E=',E12.5,3X,'L=',I1) 
      WRITE(7,3) I,R5(1),LM(I)                         
      REWIND N                                         
 11   N=N+1                                            
c      WRITE(6,4)                                      
 4    FORMAT(' Enter transfered orbital moment') 
      READ(5,*) LV                               
        Z=A( 7)                                  
       BET=A( 5)                                 
        H=A( 4)                                  
        R=A( 2)                                  
       RO1=A(11)                                 
       RO2=RO1+ A(3)*H                           
      ALFA=(RO2 - BET*DLOG(R))/R                 
       CALL VAR1(Z,R,BET,H,ALFA,RO1,RO2,R1,NE,NE5)
       CALL MATR(ALFA,BET,H,ME,R1,R2,FUNC,R5,R6,  
     *  NE,1,1,3,4,1,2,LV,NE5,4)                  
       CALL IOT3(LM(1),LV,LM(2),0,0,0,C1)         
       CALL IOT3(LM(3),LV,LM(4),0,0,0,C2)         
      C1=C1*C2*SQRT(FLOAT((2*LM(1)+1)*(2*LM(2)+1)*
     *  (2*LM(3)+1)*(2*LM(4)+1)))                 
      ME1=ME*C1                                   
      WRITE(6,5) LV,ME,ME1                        
 5    FORMAT(' Transfered orbital moment',I3/ 
     * '  Radial integral is equal ',E15.8/   
     * '  redused coulomb matrix element is  equal',E15.8/)
      WRITE(7,5) LV,ME,ME1                    
      GOTO 10                                 
      END                                     

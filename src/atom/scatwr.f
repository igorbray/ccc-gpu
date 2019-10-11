      SUBROUTINE SCATWR(LAML1,LAMM,LAMF,LAMK,LAMD,  
     * R1,R2,R6,FUNC,ME1,ME2,ME3,ME4,ME5,ME6,ME,    
     * R5,DV1,DV2,Q,C,PR,MME1,MME2,                 
     * MME4,MME5,RRR,RE,IM,RED,IMD,                 
     * F,NE,NE5,KS,NWF2,M,S,LD,KSU1,DL,IT,K,LV,     
     * F1,M1,   KK1,KKK1,K21,K22)                   
      INTEGER F,S,DL(S),KK1(S),W3,W1,W2,W,F1        
      REAL *8 LAMM(M),LAMF(F),R1(NE5),LAMK(KS),     
     * LAMD(S),PI,R,A(15),H,BET,Z,ALFA,RO1,EZZ,     
     * R2(NE5),R6(NE5),FUNC(NWF2,NE5),PRE,          
     *  EZP,ZZ1,ZZ2,ZZ3,                            
     * XNOBL,XPH,LAMI,LAML,LAML1,DELXM,DELXK,       
     * ME1(F,M),ME2(F,M),ME3(F,M),
     * DELXF,ME(F,M),ME4(M),R5(NE5),ME5(M),   
     * ME6(M),PR(1),DV1(K21),DV2(K21),C(K21),Q(K21),  
     * RE1,IM1,                                       
     * MME1(KKK1,F,M),MME2(KKK1,F,M),MME4(KKK1,M),    
     * MME5(KKK1,M),RRR(M),RE(F),AA,RE2,IM2,          
     * IM(F),RED(F),IMD(F)                            
      N1=1                                            
      PI=3.1416D 00                                   
      PRE=1.D-04                                      
      write(6, 2)                                     
2     FORMAT(/' calculation of coulomb matrixes', 
     * '')         
      READ( 3) A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8),
     * A(9),A(10),A(11),A(12),A(13),A(14),A(15)
      IF( ABS(A(1)-2.D 00).LT.PRE) GOTO 5      
      write(6, 6) (A(II),II=1,15)              
6     FORMAT(' CchiTAli Bf HE BOzbYgdEHHOgO',         
     *' COCTOiaHiia'/' MACCiB A'/,3(5(3X,E15.8)/))    
      RETURN                                          
5     CONTINUE                                        
      R=A(2)                                          
      H=A(4)                                          
      BET=A(5)                                        
      Z=A(7)                                          
      NE= IFIX(sngl(( A(3))))
      NE5=NE+5                                        
      DO 4 J=1,M                                      
      CALL TAPE(R5,XNOBL,XPH,2,NE,KOLE,3)             
      LM= IFIX(sngl((R5(2))))
      LAMM(J)=R5(1)                                   
      IA=J+S                                          
      DO 7 L=1,NE5                                    
7     FUNC(IA,L)=R5(L)                                
4     CONTINUE                                        
      REWIND 3                                        
      READ( 4) A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8),
     * A(9),A(10),A(11),A(12),A(13),A(14),A(15)
      IF( ABS(A(1)-2.D 00).LT.PRE.AND.                
     *  ABS(R-A(2)).LT.PRE.AND.                       
     *  ABS(H-A(4)).LT.PRE.AND.                       
     *  ABS(Z-A(7)).LT.PRE.AND.                       
     * NE.EQ. IFIX(sngl((A(3))))) GOTO 9                 
      write(6, 61) (A(II),II=1,15)                    
61    FORMAT(' CchiTAli Bf HE pPOMEgYTOchH',          
     *'OgO  COCTOiaHiia (22f)  ',                     
     * /' MACCiB A'/,3(5(3X,E15.8)/))                 
      RETURN                                          
9     CONTINUE                                        
      DO 8 J=1,F                                      
      CALL TAPE(R5,XNOBL,XPH,2,NE,KOLE,4)             
      LF= IFIX(sngl((R5(2))))
      LAMF(J)=R5(1)                                   
      IA=M+J+S                                        
      DO 10 L=1,NE5                                   
10    FUNC(IA,L)=R5(L)                                
8     CONTINUE                                        
      REWIND 4                                        
      READ( 2) A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8),
     * A(9),A(10),A(11),A(12),A(13),A(14),A(15)
      IF( ABS(A(1)-1.D 00).LT.PRE.AND.                 
     *  ABS(R-A(2)).LT.PRE.AND.                        
     *  ABS(H-A(4)).LT.PRE.AND.                        
     *  ABS(Z-A(7)).LT.PRE.AND.                        
     * NE.EQ. IFIX(sngl((A(3))))) GOTO 13                 
      write(6, 612) (A(II),II=1,15)                    
612   FORMAT(' CchiTAli Bf HE OCHOBHOgO  ',            
     *'  COCTOiaHiia (20f)  ',                         
     * /' MACCiB A'/,3(5(3X,E15.8)/))                  
      RETURN                                           
13    CONTINUE                                                          0800
      ALFA=A(10)                                                        0810
      RO1=A(11)                                                         0820
      write(6, 3) Z,R,NE,H                                                  0830
3     FORMAT(' zAPiad iadPA      ',12X,F5.1/                              0840
     * ' RMAX                        ',E9.2/                            0850
     * ' chiClO TOchEK iHTEgPiPOBAHiia'  ,I6  /                            0860
     * ' shAg iHTEgPiPOBAHiia          ',E9.2)                            0870
      CALL VAR1(Z,R,BET,H,ALFA,RO1,RO1+                                 0880
     * DFLOAT(NE)*H,R1,NE,NE5)                                          0890
      DO 14 J=1,S                                                       0900
      CALL TAPE( R5,XNOBL,XPH,2,NE,KOLE,2)                              0910
      IF(J.NE.LD) GOTO 15                                               0920
      LL= IFIX(sngl((R5(2))))
      LAML=LAML1                                                        0940
      IF(KSU1.EQ.0) LAML=R5(1)                                          0950
15    DL(J)= IFIX(sngl((R5(2))))
      IF(J.NE.IT) GOTO 16                                               0970
      LI= IFIX(sngl((R5(2))))
      LAMI=R5(1)                                                         0990
16    LAMD(J)=R5(1)                                                     1000
      DO 17 L=1,NE5                                                     1010
17    FUNC(J ,L)=R5(L)                                                 
14    CONTINUE                                                         
      REWIND 2                                                         
      READ( 8) A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8),
     * A(9),A(10),A(11),A(12),A(13),A(14),A(15)
      IF( ABS(A(1)-2.D 00).LT.PRE.AND.                                 
     *  ABS(R-A(2)).LT.PRE.AND.                                        
     *  ABS(H-A(4)).LT.PRE.AND.                                        
     *  ABS(Z-A(7)).LT.PRE.AND.                                        
     * NE.EQ. IFIX(sngl((A(3))))) GOTO 9787                               
      write(6, 611) (A(II),II=1,15)                                    
611   FORMAT(' CchiTAli Bf HE PACCEiaHHOgO',        
     *' elEKTPOHA     (23f)  ',                               
     * /' MACCiB A'/,3(5(3X,E15.8)/))                               
      RETURN                                                           
9787  CONTINUE                                                         
      DO 18 J=1,K                                                      
      CALL TAPE(R5,XNOBL,XPH,2,NE,KOLE,8)                              
      LK= IFIX(sngl((R5(2))))                                             
      LAMK(J)=R5(1)                                                    
      IA=M+J+S+F                                                       
      DO 19 L=1,NE5                                                    
19    FUNC(IA,L)=R5(L)                                                 
18    CONTINUE                                                         
      REWIND 8                                                         
      NM=S+1                                                           
      NF=S+1+M                                                         
      K1=0                                                             
      J=1                                                              
      DO 20 I=1,S                                                      
      IF(DL(I).NE.LK) GOTO 21                                          
      K1=K1+1                                                          
      KK1(J)=I                                                         
      J=J+1                                                            
21    CONTINUE                                                         
20    CONTINUE                                                         
      DO 22 I=1,K                                                      
      II=K-I+1                                                         
22    LAMK(II+K1)=LAMK(II)                                             
      IF(K1.LT.1) GOTO 2323                                            
      DO 23 I=1,K1                                                     
23    LAMK(I)=LAMD(KK1(I))                                             
2323  CONTINUE                                                         
      write(6, 24) LV                                                  
24    FORMAT(' pEPEdAHH@ii OPbiTAlbH@ii',  
     * ' MOMEHT  ',I2)                                                 
      IF1=F-F1                                                         
      write(6, 25) LF,F1,IF1                                           
25    FORMAT(' Bf pPOMEgYTOchHOgO',                   
     * ' COCTOiaHiia elEKTPOHA L ',I3,                    
     * /' ',I2,' diCKPET',                                       
     * 'H@X;   ',I2,' HEpPEP@BH@X fYHKtsiii') 
      write(6, 26) (LAMF(II),II=1,F)                                   
26    FORMAT(' eHEPgii:'/200(6(2X,E12.5)/))                
      write(6, 27) LI,LAMI                                             
27    FORMAT(' BPEMia HAzAd L ',I2,' E=',E12.5)                
      write(6, 28) LL,LAML                                             
28    FORMAT(' pOdObOlOchKA OCHOBHOgO ',              
     * 'COCTOiaHiia L ',I2,' E=',E12.5)                         
      IF1=M-M1                                                         
      IF(KSU1.EQ.0) GOTO 29                                            
      GOTO 30                                                          
29    write(6, 31)                                                     
 31   FORMAT('+',55X,'(TEOPETichECKOE)')                          
      GOTO 32                                                          
30    write(6, 33)                                                     
33    FORMAT('+',55X,'(eKCpEPiMEHTAlbHOE)')             
32    CONTINUE                                                         
      write(6, 333) LM,M1,IF1                                          
333   FORMAT(' Bf BOzbYgdEHHOgO   ',              
     * 'COCTOiaHiia L ',I2/' ',I2,                              
     * ' diCKPETH@X COCTOiaHiii   ',I2,                
     * ' HEpPEP@BH@X')                                        
      write(6, 34) (LAMM(II),II=1,M)                                   
34    FORMAT(' eHEPgiia='/200(6(2X,E12.5)/))                
      DELXM= SQRT(LAMM(M))- SQRT(LAMM(M-1))                            
      DELXF= SQRT(LAMF(F))- SQRT(LAMF(F-1))                            
      DELXK= SQRT(LAMK(K+K1))- SQRT(LAMK(K+K1-1))                      
      KM1=K+K1                                                         
      write(6, 35) LK,(LAMK(II),II=1,KM1)                              
35    FORMAT(' PACCEiaHH@ii elEKTPOH L ',I2,              
     * ' eHEPgiia='/  200(6(2X,E12.5)/))                    
      IF1=K+K1                                                         
      DO 37 W3=1,IF1                                                   
      NK=S+M+F+W3-K1                                                   
      IF(W3.LE.K1) NK=KK1(W3)                                          
      CALL MAEL4(ALFA,BET,H,R1,ME1,ME2,ME3,                            
     * ME,R2,FUNC,R5,R6,PR,                                            
     * NE,NK,NF,F,NM,M,LD,LL,LF,LM,LK,LV,NE5,NWF2)                     
      IF(LV.GT.LI+LK.OR.LV.LT.IABS(LK-LI))                             
     *  GOTO 40                                                        
      CALL MAEL5(ALFA,BET,H,R1,ME4,ME5,ME6,                            
     * ME,R2,FUNC,R5,R6,PR,                                            
     * NE,LD,IT,NK,NM,M,LL,LI,LK,LM,LV,NE5,NWF2)                       
40    CONTINUE                                                         
c     WRITE(9,9754) ((ME1(I9,J9),I9=1,F),J9=1,M)                       
9754  FORMAT(' ME1   DIRECT '/9(9(2X,E10.3)/))                         
c     WRITE(9,9755) ((ME2(I9,J9),I9=1,F),J9=1,M)                       
9755  FORMAT(' ME2 T. FORW. EXCHAN '/9(9(2X,E10.3)/))                  
c     WRITE(9,9756) ((ME3(I9,J9),I9=1,F),J9=1,M)                       
9756  FORMAT(' ME3    COMBINED'/9(9(2X,E10.3)/))                       
c     WRITE(9,9759) (ME4(J9),J9=1,M)                                   
9759  FORMAT(' ME4  DIRECT'/9(9(2X,E10.3)/))                           
c     WRITE(9,9757) (ME5(J9),J9=1,M)                                   
9757  FORMAT(' ME5 T. REV. EXCHAN'/9(9(2X,E10.3)/))                    
c     WRITE(9,9758) (ME6(J9),J9=1,M)                                   
9758  FORMAT(' ME6    COMBINED  '/9(9(2X,E10.3)/))                     
      DO 370 I=1,F                                                     
      DO 370 J=1,M                                                     
      MME1(W3,I,J)=ME1(I,J)                                            
370   MME2(W3,I,J)=ME2(I,J)                                            
      DO 371 J=1,M                                                     
      MME4(W3,J)=ME4(J)                                                
371   MME5(W3,J)=ME5(J)                                                
37    CONTINUE                                                         
      ZZ1=DFLOAT(IF1)                                                  
      ZZ2=DFLOAT(F)                                                    
      ZZ3=DFLOAT(M)                                                    
      WRITE(1) Z,ZZ1,ZZ2,ZZ3                                           
      do 2001 i1=1,if1
      do 2001 i2=1,f
      do 2001 i3=1,m
      WRITE(1)    MME1(I1,I2,I3)
2001  continue
      do 2011 i1=1,if1
      do 2011 i2=1,f
      do 2011 i3=1,m
      WRITE(1)    MME2(I1,I2,I3)
2011  continue
c     WRITE(1) (((MME2(I1,I2,I3),I3=1,M),I2=1,F),I1=1,IF1)             
      do 2021 i1=1,if1
      do 2021 i2=1,m
      WRITE(1)    MME4(I1,I2)
2021  continue
      do 2031 i1=1,if1
      do 2031 i2=1,m
      WRITE(1)    MME5(I1,I2)
2031  continue
c      WRITE(1)  ((MME5(I1,I2),I2=1,M),I1=1,IF1)                       
      ZZ1=DFLOAT(M1)                                                   
      ZZ2=DFLOAT(F1)                                                   
      ZZ3=DFLOAT(LV)                                                   
      WRITE(1) ZZ1,ZZ2,ZZ3                                             
      ZZ1=DFLOAT(LK)                                                   
      ZZ2=DFLOAT(LI)                                                   
      ZZ3=DFLOAT(KS)                                                   
      WRITE(1) ZZ1,ZZ2,ZZ3                                             
      WRITE(1) DELXM,DELXF,LAMI,LAML                                   
      do 2041 io=1,m
      WRITE(1)    lamm(io)
2041  continue
      do 2051 io=1,f
      WRITE(1)    lamf(io)
2051  continue
      do 2061 io=1,ks
      WRITE(1)    lamk(io)
2061  continue
c     WRITE(1) (LAMM(IO),IO=1,M),(LAMF(IO),IO=1,F),(LAMK(IO),IO=1,KS)  
      ZZ1=DFLOAT(K1)                                                   
      ZZ2=DFLOAT(K)                                                    
      ZZ3=DFLOAT(NE)                                                   
      WRITE(1) ZZ1,ZZ2,R,H,BET,ZZ3,ALFA,RO1                            
      do 2071 io=1,ne5
      WRITE(1)    r1(io)
2071  continue
c     WRITE(1) (R1(IO),IO=1,NE5)                                       
      write(6, 6787)                                                      
6787  FORMAT('  KYlOHOBCKAia MATPitsA  zApiCAHA')   
      RETURN                
      END                   


      SUBROUTINE  komon2(z,r,h,lamk,lamf,lami,laml,lamm,
     *  aa1,aa2,aa3,aa4,
     * ne,lk,if1,f1,f,ikis,lv,li,m1,m)
      integer ne,lk,if1,f1,f,ikis,lv,li,m1,m
      real*8  z,r,h,lamk(if1),lamf(f),lami,laml,lamm(m),
     *  aa1(1),aa2(1),aa3(1),aa4(1)
      write(66, 6051) Z,R,NE,H
6051  FORMAT(' nuclear charge',18X, 
     *F5.1/' RMAX',24X,F5.1/        
     *' number of integration points',I6/
     *' step of integration         ',E9.2) 
      write(66, 6053) LK,(LAMK(IO),IO=1,IF1)
6053  FORMAT(' scattering electron L=',
     *I1,' Energy:'/3(10(1X,E11.5)/))    
      IFN=F-F1                                        
      write(66, 6052) IKIS,LV,F1,IFN                   
6052  FORMAT(/' transition',I2/' transfered orbital moment',   
     * I3/' WF intermediate state:',              
     * I2,' discrete',I3,' continuum')                             
      write(66, 6054) (LAMF(IO),IO=1,F)
6054  FORMAT(' energy:',8(2X,E11.5)/2(10X,8(2X,E11.5)/))   
      write(66, 6055) LI,LAMI                               
6055  FORMAT(' time reverse L=',I1,4X,'E=',E12.5)  
      write(66, 6056) LAML   
6056  FORMAT(' sub-sell of ground state E=',   
     * E12.5)                                  
      MN=M-M1                                  
      write(66, 6057) M1,MN                     
6057  FORMAT(' WF excited state :',
     * I2,' discrete',I3,' continuum')
      write(66, 6054) (LAMM(IO),IO=1,M)
      write(66, 9444) AA1(IKIS),AA2(IKIS),AA3(IKIS),AA4(IKIS)    
9444  FORMAT(' Coefficients: ',4(2X,F8.5,3X)) 
      RETURN                                  
      END                                     

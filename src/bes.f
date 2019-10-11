      function BES(L,X)                        
                                               
C  Last update: 11 Jun 1998                    
                                               
C  This function calculates a modifyed         
C  spherical Bessel function                   
C                                              
C   j (X) = sqrt(pi/2*X) J  (X)                
C    L                    L+1/2                
C                                              
C  which at small |X|<1 and large |X|>1 is     
C  given by formulas (10.1.2) and (10.1.8)     
C  of the Abramovitz and Stegan's Handbook     
                                               
      IMPLICIT REAL*8(A-H,O-Z)                 
      IMPLICIT INTEGER*4(I-N)
      real*8 X,pi
      
      pi = dacos(-1.0d0)                       
                                               
      if (L .lt. 0) then                       
         bes = 0.0d0                           
         return                                
      end if                                   
                                               
      if (L .eq. 0) then                       
         if(x .eq. 0d0) then    ! j   = sin(x)/x  
            bes = 1.            !  L=0            
         else                                 
            bes = sin(x)/x
         end if                               
         return                               
      end if                                  
                                              
      if (L .eq. 1) then                      
         if(x .eq. 0d0) then    ! j   = sin(x)/x**2 - cos(x)/x 
            bes = 0.            !  L=1            
         else                                 
            bes = sin(x)/x**2 - cos(x)/x      
         end if                               
         return                               
      end if                                  
                                              
      if (L .eq. 2) then                      
         if(x .eq. 0d0) then    ! j   = sin(x)(3/x**3-1/x) 
            bes = 0.            !  L=2          - 3 cos(x)/x**2         
         else                                 
            bes = sin(x)*(3/x**3 - 1/x) - 3*cos(x)/x**2      
         end if                               
         return                               
      end if                                  
                                              
      m1 = int(L/2d0)                         
      m2 = int(abs((L-1)/2d0))                
      if ((X .gt. 1d0) .and. (X .gt. dble(L))) then  !Loop for large |X|
         c  = 1
         pp = 1
         d  = dble(L*(L+1)) / (2d0*X)
         qq = d
         if (m1 .ne. 0) then
            do i = 1, m1
               pp = - pp * dble(L-2*i+1) * dble(L-2*i+2)  /
     :                         (dble(  2*i)   * dble(  2*i-1)) *
     :                          dble(L+2*i)   * dble(L+2*i-1) / 
     :               (2d0*X)**2     
               c  = c + pp
            end do
         end if
         if (m2 .ne. 0) then
            do i = 1, m2
               qq = - qq * dble(L+2*i) * dble(L+2*i+1)  /
     :                         (dble(  2*i) * dble(  2*i+1)) *
     :                          dble(L-2*i) * dble(L-2*i+1) /
     :               (2d0*X)**2
               d = d + qq
            end do
         end if   
         bes = (c * sin(X - pi/2d0 * dble(L))  + 
     :          d * cos(X - pi/2d0 * dble(L))) / X
      else                                          !Loop for small |X|
         ii = 0
         bes = 0
         a1 = 1/dble(2*L+1)
         if (l .ne. 0) then
            do i = 1, L
               a1 = a1 * 2d0 / dble(i+L) * X
            end do
         end if   
         if (dabs(a1) .lt. 10d-11 ) goto 20
 10      bes = bes + a1
         ii = ii + 1
         a1 = a1 * (-1d0) * X**2 * dble(L+ii) / dble(ii) /
     :                       dble(2*(L+ii)+1) / dble(2*(L+ii))         
         if (dabs(a1) .gt. 10d-11) goto 10
 20      continue
      end if

      END            



















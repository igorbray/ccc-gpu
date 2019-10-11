      subroutine OVER(P, Q, result)                                   
                                                            
*     Radial overlap integral between two bound states      
*                   oo                                      
*                   f                                       
*              I =  | dr P(r) * Q(r)                        
*                   j                                       
*                   o

      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)

!Simpson's weights added
      result = 0.d0
      do i=1,meshr
         result = result + P(i)*Q(i) * rmesh(i,3)
      end do 

      RETURN
      END


      subroutine OVER1(maxNR,P, Q, result)
      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      COMMON           /TELNUM/KOUNT

!      kount = kount+1
!      print'(A,I3,I5)', '  OVER1', kount, maxNR
      result = 0.d0
      do i=1,maxNR
         rt = rmesh(i,1)                                
         result = result + P(i)*Q(i) * rmesh(i,3) 
!         write(300+kount,'(10F9.4)') rt,P(i),Q(i),result 
      end do 

      RETURN
      END




      subroutine rOVER(maxNR,P, Q, result)

C Radial overlap integral between two bound states      
                                                        
c                   oo                                  
c                   f                                   
c              I =  | r  dr P(r) * Q(r)                 
c                   j                                   
c                   o                                   
      include 'par.f'                                   
      DIMENSION  P(maxr), Q(maxr)                       
      common/meshrr/ meshr,rmesh(maxr,3)                
      COMMON           /TELNUM/KOUNT

!      kount = kount+1
!      print'(A,I3,I5)', '  rOVER', kount, maxNR
      result = 0.d0                                     
      do i=1,maxNR                                      
         rt = rmesh(i,1)                                
         result = result + P(i)*Q(i)*rt* rmesh(i,3) !Simded
!         write(300+kount,'(10F9.4)') rt,P(i),Q(i),result 
      end do                                            
                                                        
      RETURN                                            
      END                                               
                                                        
                                                        
                                                        
                                                        
      subroutine rOVER2(maxNR,P, Q, result)                     
C Radial overlap integral between two bound states      
                                                        
c                   oo                                  
c                   f  2                                
c              I =  | r  dr P(r) * Q(r)                 
c                   j                                   
c                   o  
                                        
      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      COMMON           /TELNUM/KOUNT
!      kount = kount+1
!      print'(A,I3,I5)', ' rOVER2', kount, maxNR
      result = 0.d0
      do i=1,maxNR
         rt = rmesh(i,1)
         result = result + P(i)*Q(i)*rt**2*rmesh(i,3)
!         write(300+kount,'(10F9.4)') rt,rt**2,P(i),Q(i),result 
      end do 

      RETURN
      END



      subroutine rOVER2r(maxNR,P, Q, result)                     
C Radial overlap integral between two bound states      
                                                        
c                   oo                                  
c                   f      2                                
c              I =  |j1(qr)  dr P(r) * Q(r)                 
c                   j                                   
c                   o  
      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      COMMON           /TELNUM/KOUNT
      COMMON /QGR/ fqr,fqr2
      DIMENSION    fqr(maxr),fqr2(maxr,maxr)
      result = 0.d0
      do i=1,maxNR
         rt = rmesh(i,1)
         result = result + P(i)*Q(i)*rt**2*fqr(i)*rmesh(i,3)
!         write(300+kount,'(10F9.4)') rt,rt**2,P(i),Q(i),result 
      end do 

      RETURN
      END


      subroutine rOVER22r(maxNR,P,Q,F,G, res)                     
C Radial overlap integral between two bound states      
                                                        
c                   oo             oo               
c                   f              f                   
c              I =  |rdr P(r)*Q(r) |r'dr' F(r')*G(r') Rho(r,r')
c                   j              j               
c                   o              o   
      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      DIMENSION  F(maxr), G(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      COMMON           /TELNUM/KOUNT
      COMMON /QGR/ fqr,fqr2
      DIMENSION    fqr(maxr),fqr2(maxr,maxr)

      SUM = 0.d0
      res = 0.d0
      do i=1,maxNR
         rt = rmesh(i,1)
         do j=1,maxNR
            rp = rmesh(j,1)
            res = res + P(i)*Q(i)*rt*rmesh(i,3)
     >                * F(j)*G(j)*rp*rmesh(j,3)
     >                * fqr2(i,j)
         end do 
      end do 

      RETURN
      END

      subroutine vOVER(maxNR,P,Q,lp,lq, result)

C Velocity operator taken between two bound states

c                   oo
c                   f
c              I =  |drP(r)(-v +- Lmax/r)Q(r) , v-differential operator
c                   j
c                   o

      include 'par.f'
      DIMENSION  P(maxr), Q(maxr), W(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)

      lm=lp                     !Lmax=max(lp,lq)
      if(lq.gt.lp) lm=lq

      il=-1                     !Lq=Lp-1
      if(lm-lp.gt.0) il=1       !Lq=Lp+1

      call  DERIVATIVE(P,W)


      result = 0.d0
      do i=1,maxNR
         result = result + Q(i)
     >      *(-W(i)+dble(il*Lm)*
     >      P(i)/rmesh(i,1))
     >      * rmesh(i,3)                         !Simpson's weights added
      end do

      RETURN
      END

      subroutine vOVER2(maxNR,P,Q,lp,lq, result)

C Two velocity operators between two bound states

c                   oo
c                   f       /             \2
c              I =  |drP(r) |-v +- Lmax/r)|  Q(r) , v-differential operator
c                   j       \             /
c                   o

      include 'par.f'
      DIMENSION  P(maxr), Q(maxr), W(maxr), U(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)

      lm=lp                     !Lmax=max(lp,lq)
      if(lq.gt.lp) lm=lq

      il=-1                     !Lq=Lp-1
      if(lm-lp.gt.0) il=1       !Lq=Lp+1

      call  DERIVATIVE(P,W)


      do i=1,maxNR
         U(i) = -W(i)+dble(il*Lm)* P(i)/rmesh(i,1)
      end do

      call  DERIVATIVE(U,W)

      result = 0.d0
      do i=1,maxNR
         result = result + Q(i)
     >      *(-W(i)+dble(il*Lm)*
     >      U(i)/rmesh(i,1))
     >      * rmesh(i,3)                         !Simpson's weights added
      end do

      RETURN
      END


      subroutine aOVER(maxNR,P, Q, result)

C Radial overlap integral between two bound states      
                                                        
c                   oo                                  
c                   f Z                                  
c              I =  |--2-  dr P(r) * Q(r)                 
c                   j r                                  
c                   o                                   
      include 'par.f'                                   
      DIMENSION  P(maxr), Q(maxr)                       
      common/meshrr/ meshr,rmesh(maxr,3)                
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym

      result = 0.d0                                     
      do i=1,maxNR                                      
         rt = rmesh(i,1)                                
         result = result + P(i)*Q(i)/rt**2 * nznuc*rmesh(i,3) !Simded
      end do                                            
                                                        
      RETURN                                            
      END                                               
                                                        
                                                        


      subroutine rOVERn(maxNR,P, Q, n, result)                     
C Radial overlap integral between two bound states      
                                                        
c                   oo                                  
c                   f  n                                
c              I =  | r  dr P(r) * Q(r)                 
c                   j                                   
c                   o  
                                        
      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      COMMON           /TELNUM/KOUNT
!      kount = kount+1

      result = 0.d0
      do i=1,maxNR
         rt = rmesh(i,1)
         result = result + P(i)*Q(i)*rt**n*rmesh(i,3)
!         write(300+kount,'(10F9.4)') rt,rt**n,P(i),Q(i),result 
      end do 

      RETURN
      END


      subroutine bOVER1(maxNR,P, Q, result)

C Radial integral with spherical Bessel function
c scaled to r-operator
      
c                   oo                                  
c                 3  f                                  
c              I =-  |j (q r)  dr P(r) * Q(r)                 
c                 q  j 1
c                    o             

      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      common/QS/ q0,q1
      real*8 x

      qav = (q0+q1)/2
      result = 0.d0
      do i=1,maxNR
         rt = rmesh(i,1)
         x = qav*RT
         tmp = bes(1,x)*3/qav
         result = result + P(i)*Q(i) * tmp * rmesh(i,3) 
!         write(33,'(10F9.4)') rt,tmp,rt,P(i),Q(i),result 
      end do 

      RETURN
      END



      subroutine bOVER2(maxNR,P, Q, result)                     

C Radial integral with two spherical Bessel functions 
c scaled to r^2-operator
      
c                      oo                                  
c                 3 3  f                                  
c              I =- -  |j (q r) j (q r)  dr P(r) * Q(r)                 
c                 q q  j 1  0    1  1                               
c                  0 1             

      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      common/QS/ q0,q1
      real*8 x0,x1
     
      result = 0.d0
      do i=1,maxNR
         rt = rmesh(i,1)
         x0 = q0*RT
         x1 = q1*RT
         tmp = bes(1,x0)* bes(1,x1)*3/q0*3/q1
         result = result + P(i)*Q(i)*tmp*rmesh(i,3)
!         write(34,'(10F9.4)') rt,tmp,rt**2,P(i),Q(i),result 
      end do 

      RETURN
      END



      subroutine bOVERl(maxNR,P, Q, l, q0, result)

C Radial integral with spherical Bessel function
      
c                  oo                                  
c                   f                                  
c              I =  |j (q0 r)  dr P(r) * Q(r)                 
c                   j l
c                   o             

      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      real*8 x
      COMMON           /TELNUM/KOUNT

      kount = kount+1
!      print'(A,10F9.4)', ' bOVERl', q0,(Q(i),i=1,7)
      result = 0.d0
      do i=1,maxNR
         rt = rmesh(i,1)
         x = q0*RT
         tmp = bes(l,x)
         if(l.eq.0) tmp = bes(l,x)-1.
!         if(l.eq.1) tmp = tmp/q0*3 !Back to old 2nd Born
!         tmp = q0*RT/3
         result = result + P(i)*Q(i) * tmp * rmesh(i,3) 
!         write(3000+kount,'(10F9.4)') rt,P(i),Q(i),result 
      end do 

      RETURN
      END


      subroutine bOVERll(maxNR,P, Q, l0,l1, q0,q1, result)

C Radial integral with two spherical Bessel function
      
c                  oo                                  
c                   f                                  
c              I =  |j (q0 r) j (q1 r)  dr P(r) * Q(r)                 
c                   j l        lp       
c                   o             

      include 'par.f'
      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      real*8 x0,x1
      COMMON           /TELNUM/KOUNT

      kount = kount+1
!      if(kount.eq.1)
!     >   print'(A,I3,I5,2F9.4)', 'bOVERll', kount, maxNR,q0,q1
      result = 0.d0
      do i=1,maxNR
         rt = rmesh(i,1)
         x0 = q0*RT
         x1 = q1*RT
         tmp0 = bes(l0,x0)
         tmp1 = bes(l1,x1)
         if(l0.eq.0) tmp0 = bes(l0,x0)-1.
         if(l1.eq.0) tmp1 = bes(l1,x1)-1.
!         if(l0 .eq. 1) tmp0 = tmp0/q0*3 !Back to old 2nd Born
!         if(l1 .eq. 1) tmp1 = tmp1/q1*3
!         tmp0 = q0*RT/3
!         tmp1 = q1*RT/3

         result = result + P(i)*Q(i) * tmp0 * tmp1 * rmesh(i,3) 
!      if(kount.eq.1)
!     >      write(500+kount,'(10E13.4)') rt,tmp0*tmp1,P(i),Q(i),result 
      end do 
      RETURN
      END









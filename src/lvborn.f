      subroutine lvBORN(maxNi,wi,wf,li,lf,lv,q,resl,resv)    

C  Last update: 11 Jun 1999          iqr     
C                                 \ e        
C  Calculates Born ME in both   i  \    f    
C  length and velocity forms     ---*---     
                                             
C  Velocity form defined according to        
C  Amusia, Chernysheva 1983 (9.11)           
C  Extra 1/2 is added to be consistent with  
C  photoionization subroutine VELOCITY       

      include 'par.f'                                  
      real*8 x,ang

      common/meshrr/ meshr,rmesh(maxr,3)
                                                       
      DIMENSION  wi(maxr), wf(maxr), yi(maxr),yf(maxr),yd(maxr)
                                                   
      call IOT3(LF,LV,LI,0,0,0,ang)
      ddff=hat(LI)*hat(LF)                         
      ang=ang*ddff                                 

C Born ME in length form  

      resl = 0.
      do J=1,maxNi
         RT=rmesh(J,1)                                  
         WT=rmesh(J,3)          !Simpson's weights added
         x = Q*RT
         BE= BES(LV,x)                          
         if(LV.eq.0)BE= BES(LV,x)-1.0
         resl = resl + WI(J)*WF(J) * BE * WT 
      end do                                       
      resl=resl*ang
                                                   
C Born ME in velocity form  
      
      A1=0.                                        
      LA=LI*(LI+1)-LF*(LF+1)                           
      IF(LA.EQ.0) GOTO 500                         
      DO I=1,maxNi                       
         RT=rmesh(i,1)                        
         WT=rmesh(i,3)                        ! wi*wf -term      
         x = Q*RT
         BE= BES(LV-1,x)                      !
         B = BES(LV+1,x)                      ! J (qr) + J (qr) 
         BE=BE+B                              !  l-1      l+1    
         if(RT.ne.0.)
     :      A1= A1 + WI(I)*WF(I)*BE /RT * WT
      end do                                             
      A1=A1*LA
 500  CONTINUE                                     
                    
C Initialize derivative. Otherwise crash

      do J=1,maxr
         YI(J)=0.      
         YF(J)=0.      
      end do                                       

      call DERIVATIVE(WI,YI)                           
      call DERIVATIVE(WF,YF)                           
      DO I=1,maxNi                        
         RT=rmesh(i,1)                                  
         WT=rmesh(i,3)                                  
         x = Q*RT
         YD(I)=WI(I)*YF(I)-WF(I)*YI(I)        ! wi*d/dr(wf) - wf*d/dr(wi)
         BE= BES(LV+1,x)                      !   
         B = BES(LV-1,x)                      !(l+1) J (qr) - l J (qr)
         BE=BE*(LV+1)-B*LV                    !       l+1        l-1
         YD(I)=YD(I)*BE                            
         A1=A1+YD(I) * WT
c$$$         PB=PB+YD(I) * WT
      end do                                       

      A1=A1*Q/hat(LV)**2                      
                                                  
      resv=ang*A1                             !Born ME in the VELOCITY form  
      resv=-resv/2                            !Extra -1/2    
                                                   
*      if(lv .eq. 1) resl = resl*3/q          !Back to photo
*      if(lv .eq. 1) resv = resv*3/q

      

      RETURN
      END


      function hat(l)

c$$$      real*8 hat
      integer l

      hat= sqrt(2.0 * l + 1.0)
      return
      end
      







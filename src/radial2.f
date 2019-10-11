      subroutine RADIAL2(w,w2,l,l1,J,r,v)
                                                       
      include 'par.f'                                            
      PARAMETER (Ncf = 20)
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr)
      
      COMMON /cat1/   CI(0:Ncf, Ncf, Ncf)
      COMMON /cat2/   Nc, nmax(0:Ncf), Limax
      COMMON /cat3/   wa(maxr, 0:Ncf, Ncf)
      COMMON /cat5/   maxNR(0:Ncf, Ncf)
      COMMON           /TELNUM/KOUNT
      kount=0

C Get rid of the Simpson's weights for continuum WF
C But don't corrupt them for later use. That's why w2, not w1

!      do i=1,maxr
!         if(rmesh(i,3) .ne. 0d0) w1(i) = w2(i)/rmesh(i,3)
!      end do
      do i=1,meshr
         w1(i) = w2(i)/rmesh(i,3)
      end do

C MCHF loop
      r=0.                                                  
      v=0.                                                  
      do Li=0,Limax                         
         do Ni=Li+1,nmax(Li)
            call sBORN2(maxNR(Li,Ni),
     :       wa(1,Li,Ni),w,w1,Li,l,l1,J,CI(Li,Ni,Ni),r,v)
         end do
      end do
                    
      RETURN
      END

      subroutine sBORN2(maxNi,wi,wj,wk,li,lj,lk,lv,C,r,v)     

      include 'par.f'                                              
      dimension  wi(maxr), wj(maxr), wk(maxr)                      
      real*8 six,threea,threeb
      character*1 hh(0:3)         !Spectroscopic labels
      data hh /'s','p','d','f'/
      COMMON           /TELNUM/KOUNT
!      goto 111 !D2b only
      
! D2a term                     !    i \  lk               
      t1 = triang(lk,1,li )    !   ---------
      t2 = triang(lj,1,li)     !     Lv=0,2    
      t3 =striang(lk,Lv,lj)    !   ---------
                               !    i /  lj
      
      if(t1*t2*t3 .ne. 0d0) then
!         call IOT6(1, Lv,1, lk,li,lj, six)
         six=COF6J(1.0, float(Lv),1.0, float(lk),float(li),float(lj))
         call IOT3(lk,1,li, 0,0,0,threea)
         call IOT3(lj,1,li, 0,0,0,threeb)
         call rover(maxNi,wi, wk, resa)                  
         call rover(maxNi,wi, wj, resb)                  
!         call aover(maxNi,wi, wk, vesa)                  
!         call aover(maxNi,wi, wj, vesb)                  
         call vover(maxNi,wi, wk,li,lk, vesa)                  
         call vover(maxNi,wi, wj,li,lj, vesb)                  
!        call rover22r(maxNi,wi,wk,wi,wj,resc)
         r = r + 2.0 * C  * (-1)**Lv * hat(lk)*hat(li)*hat(lj)
     >      * threea * threeb * six * resa * resb
         v = v + 2.0 * C  * (-1)**Lv * hat(lk)*hat(li)*hat(lj)
     >      * threea * threeb * six * vesa * vesb
!     >     * threea * threeb * six * resc
!         print'(/A))', 'D2a direct+exchange'
!         print'(2x,A,4x,4(A,8x))', 'li lj Lv lk','t1','t2','t3'
!         print'(4I3,4F9.4)', li,lj,Lv,lk,t1,t2,t3
!         print'(10F9.4)',C,threea, threeb, six, resa, resb,r
      end if
 111  continue
!      RETURN !D2a only
      
! D2b term - part 1                       !  i \/  k      
      t1 = triang(1, Lv,1)                ! --------- 
      t2 = triang(lk,Lv,li)               !  i Lv  j  
                                          ! ----x---- 
      call IOT3(1, 1,Lv, 0,0,0,threea) !Bug is fixed!
      if(t1*t2.ne.0 .and. li.eq.lj) then    
         call IOT3(lk,Lv,li,0,0,0,threeb)
         call rover2(maxNi,wi,wk,resa)
         call vover2(maxNi,wi,wk,li,lk,vesa)
         call over1 (maxNi,wi,wj,resb)
         r = r + C * (-1)**(li) / hat(li) * hat(lk) * hat(li)
     >      * threea * threeb  * resa * resb
         v = v + C * (-1)**(li) / hat(li) * hat(lk) * hat(li)
     >      * threea * threeb  * vesa * resb
!         print'(/A))', 'D2b direct'
!         print'(2x,A,4x,4(A,8x))', 'li lj Lv lk','t1','t2'
!         print'(4I3,4F9.4)', li,lj,Lv,lk,t1,t2
!         print'(10F9.4)',C,threea, threeb,  resa, resb,r
      end if                  
                                          ! ----x--- 
                                          !  i Lv k       
! D2b term - part 2                       ! --------      
      t2 = triang(lj,Lv,li)               !  i /\ j 
      if(t1*t2.ne.0 .and. li.eq.lk) then  
         call IOT3(lj,Lv,li,0,0,0,threeb)
         call over1 (maxNi,wi,wk,resa)
         call rover2(maxNi,wi,wj,resb)
         call vover2(maxNi,wi,wj,li,lj,vesb)
         r = r + C * (-1)**(li) / hat(li) * hat(lj) * hat(li)
     >      * threea * threeb * resa * resb
         v = v + C * (-1)**(li) / hat(li) * hat(lj) * hat(li)
     >      * threea * threeb * resa * vesb
!         print'(/A))', 'D2b exchange'
!         print'(2x,A,4x,4(A,8x))', 'li lj Lv lk','t1','t2'
!         print'(4I3,4F9.4)', li,lj,Lv,lk,t1,t2
!         print'(10F9.4)',C,threea, threeb,  resa, resb,r
      end if
      RETURN                                                  

      END
      

      function sTRIANG(l1, l2, l3)
    
C The L's must satisfy the triangular inequality
C "Soft" triangle without inforced parity
    
      striang = 1.
      if (l1+l2.lt.l3.or.l2+l3.lt.l1.or.l3+l1.lt.l2) striang = 0.

      return
      end

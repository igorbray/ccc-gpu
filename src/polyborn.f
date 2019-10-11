      subroutine polyBORN(w,w2,l,l1,J,Q,r1,v1,r,v)
                                                       
      include 'par.f'                                            
      PARAMETER (Ncf = 20)
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr)
      
      COMMON /cat1/   CI(0:Ncf, Ncf,Ncf)
      COMMON /cat2/   Nc, nmax(0:Ncf), Limax
      COMMON /cat3/   wa(maxr, 0:Ncf, Ncf)
      COMMON /cat5/   maxNR(0:Ncf, Ncf)

C Get rid of the Simpson's weights for continuum WF
C But don't corrupt them for later use. That's why w2, not w1

!      do i=1,maxr
!         if(rmesh(i,3) .ne. 0d0) w1(i) = w2(i)/rmesh(i,3)
!      end do
      do i=1,meshr
         w1(i) = w2(i)/rmesh(i,3)
      end do


*      if(J.eq.0) then
*         do i=1,maxr
*            write(213,'(2E12.4)')rmesh(i,1),w1(i)
*         end do
*      end if

      r=0.                                                  
      v=0.
      r1=0.
      v1=0.
C No correlation. S->J channel only

      call singleBORN(maxNR(0,1),wa(1,0,1),w,w1,0,l,l1,J,1.,Q,r1,v1)

C Angular momentum loop

      do Li=0,Limax                         
         do Ni=Li+1,nmax(Li)
            call singleBORN(
     :       maxNR(Li,Ni),wa(1,Li,Ni),w,w1,Li,l,l1,J,CI(Li,Ni,Ni),Q,r,v)
         end do
      end do
                              
                    
      RETURN
      END


      subroutine singleBORN(maxNi,wi,wj,wk,li,lj,lk,lv,C,Q,r,v)     

      include 'par.f'                                              
      dimension  wi(maxr), wj(maxr), wk(maxr)                      
      character*1 hh(0:3)         !Spectroscopic labels
      data hh /'s','p','d','f'/


C Phase factor: ang=1 for i=s

      t1 = triang(li, lv, lk)
      t2 = triang(li, lv, lj)

C Direct term: li=lj
      
      if(li.eq.lj .and. t1.ne.0.) then                                
         ang = (-1)**li / hat(li)                         ! \                  
         call OVER1(maxNi,wi, wj, ov)                     !   \ lv             
         call lvBORN(maxNi,wi,wk,li,lk,lv,q,resl,resv)    !i --*-- k           
         r = r + resl * ov * ang * C                      !i --x-- j           
         v = v + resv * ov * ang * C
!      write(6,'(A,I1,A1,7F9.4)')
!     :   'T1 ',1,hh(Li), ang,ov,C, resl,  resv ,r,v
      end if
      
C Exchange term: li = lk      

      if(li.eq.lk .and. t2.ne.0.) then              
         ang = (-1)**li / hat(li)                         ! \ 
         call OVER1(maxNi,wi, wk, ov)                     !   \ lv      
         call lvBORN(maxNi,wi,wj,li,lj,lv,q,resl,resv)    !i --*-- j    
         r = r + resl * ov * ang * C                      !i --x-- k           
         v = v + resv * ov * ang * C
!      write(6,'(A,I1,A1,7F9.4)')
!     :   'T2 ',1,hh(Li), ang,ov,C, resl,  resv ,r,v
      end if
                                                                   
      RETURN                                                                 
      END                                                                    



      function TRIANG(l1, l2, l3)
    
C The L's must satisfy the triangular inequality
C The sum of actual L's must be an integer
    
      triang = 1.
      if (l1+l2.lt.l3.or.l2+l3.lt.l1.or.l3+l1.lt.l2) triang = 0.
      if (mod(l1+l2+l3,2).ne.0) triang = 0.

      return
      end


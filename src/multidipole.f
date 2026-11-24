      subroutine multiDIPOLE(w,w2,l,l1,r1,v1,x1,r,v,x) 
                                                       
      include 'par.f'                                            
      PARAMETER (Ncf = 20)
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr), f(maxr), f1(maxr)
      
      COMMON /cat1/   CI(0:Ncf, Ncf, Ncf)
      COMMON /cat2/   Nc, nmax(0:Ncf), jmax
      COMMON /cat3/   wa(maxr, 0:Ncf, Ncf)
      COMMON /gstspin/   iSpin, meta

C Get rid of the Simpson's weights for continuum WF
C But don't corrupt them for later use. That's why w2, not w1

!      do i=1,maxr
!         if(rmesh(i,3) .ne. 0d0) w1(i) = w2(i)/rmesh(i,3)
!      end do
      do i=1,meshr
         w1(i) = w2(i)/rmesh(i,3)
      end do

      rm=0.                                                  
      vm=0.
      xm=0.
      rm1=0.
      vm1=0.
      xm1=0.

C Angular momentum loop
      
      do J=0,jmax
         
         I1=0
         I2=0
         
         if(l.eq. J .and. l1.eq.J+1) I1=1 !J,J+1 or J+1,J-combination
         if(l1.eq.J .and. l. eq.J+1) I2=1 !Always f=J, f1=J+1
         if(I1.eq.1)then
            call EQUAL(w,f)               !    \          
            call EQUAL(w1,f1)             ! J __\__ f1=J+1
         end if                 ! J __.__ f =J  
         if(I2.eq.1)then
            call EQUAL(w1,f)
            call EQUAL(w,f1)
         end if
         if(I1+I2.eq.0)goto 111

C J->J+1 transition

         do n=J+1,nmax(J)
            do n1=J+1,nmax(J)
               if(CI(J,n,n1).ne.0.)then
                  if(n.eq.n1)
     :               call singleDIPOLE
     :               (wa(1,J,n),f,f1,J,J,J+1,CI(J,n,n1),rm,vm,xm)
                  if(n.ne.n1)
     :               call doubleDIPOLE
     :    (wa(1,J,n),wa(1,J,n1),f,f1,J,J,J+1,CI(J,n,n1),rm,vm,xm)
               end if
           end do
        end do
                              
C No correlation. SP channel only

         if(J.eq.0 .and. l*l1 .eq. 0) then
            if(meta.eq.0)
     :         call singleDIPOLE
     :         (wa(1,0,1),f,f1,0,0,1,1.,rm1,vm1,xm1)
            if(meta.eq.1)
     :         call doubleDIPOLE
     :         (wa(1,0,1),wa(1,0,2),f,f1,0,0,1,1.,rm1,vm1,xm1)
         end if
                    
C J+1->J transition

         do n=J+2,nmax(J+1)
            do n1=J+2,nmax(J+1)
               if(CI(J,n,n1).ne.0.)then
                  if(n.eq.n1)
     :               call singleDIPOLE
     :              (wa(1,J+1,n),f1,f,J+1,J+1,J,CI(J+1,n,n1),rm,vm,xm)
                  if(n.ne.n1)
     :               call doubleDIPOLE
     : (wa(1,J+1,n),wa(1,J+1,n1),f1,f,J+1,J+1,J,CI(J+1,n,n1),rm,vm,xm)
               end if
            end do
         end do
         
 111     continue
      end do

      io = 1
      if(iSpin.eq.1) io=(-1)**l
                       
      
      r =rm *io
      v =vm *io
      x =xm *io
      r1=rm1*io
      v1=vm1*io
      x1=xm1*io
      
      RETURN
      END

      subroutine EQUAL(w,w1)                 !Substitute WF "w1" with "w"
      
      include 'par.f'
      
      dimension   w(maxr),w1(maxr)
      do i=1,maxr
         w1(i) = w(i)
      end do
      
      END
      
      subroutine singleDIPOLE(g,f,f1,lg,l,l1,c,r,v,x) 

C Input/output: r,v,x (length, velocity, accelration)
                                                       
      include 'par.f'                                            
           
      DIMENSION   g(maxr), f(maxr), f1(maxr)

C Angular momentum check

      if(lg .ne. l .or. iabs(lg-l1) .ne. 1) then
         write(6,'(3(A,I2))') 'Lo=', lg, ' L=', l, ' L1=', l1
         STOP 'WRONG ANGULAR MOMENTA'
      end if
      
      hat = dble(2*lg+1)                      !    l   ______   
      ang=(-1)**lg/sqrt(hat)                  !(-1) / V(2l+1)

      call OVER(g,f,ove)                       
      call LENGTH  (g,f1,lg,l1,ar)           
      call VELOCITY(g,f1,lg,l1,av)           
      call ACCEL   (g,f1,lg,l1,ax)           

      r = r + ar*ove*C*ang                    !     \      
      v = v + av*ove*C*ang                    ! g  __\__ f1
      x = x + ax*ove*C*ang                    ! g  __.__ f 

*      write(6,'(6(A,1p,E11.4))')'C=',C,' ANG=', ang, ' OVERLAP=', OVE,
*     :   ' L', r,
*     :   ' V', v,
*     :   ' A', x

      return
      end
      
      subroutine doubleDIPOLE(g1,g2,f,f1,lg,l,l1,c,r,v,x) 

C Input/output: r,v,x (length, velocity, accelration)
                                                       
      include 'par.f'                                            
      COMMON /gstspin/   iSpin, meta
           
      DIMENSION   g1(maxr), g2(maxr), f(maxr), f1(maxr)

C Angular momentum check

      if(lg .ne. l .or. iabs(lg-l1) .ne. 1) then
         write(6,'(3(A,I2))') 'Lo=', lg, ' L=', l, ' L1=', l1
         STOP 'WRONG ANGULAR MOMENTA'
      end if
      
      hat = dble(2*lg+1)                      !    l   ______   
      ang=(-1)**lg/sqrt(hat)                  !(-1) / V(2l+1)

      if(iSpin.eq.1) ang=1/sqrt(hat) !No reason for this. Just guess!
      
      call OVER(g2,f,ove2)                       
      call LENGTH  (g1,f1,lg,l1,ar1)          !     \      
      call VELOCITY(g1,f1,lg,l1,av1)          ! g1 __\__ f1
      call ACCEL   (g1,f1,lg,l1,ax1)          ! g2 __.__ f 

      call OVER(g1,f,ove1)                       
      call LENGTH  (g2,f1,lg,l1,ar2)          !     \      
      call VELOCITY(g2,f1,lg,l1,av2)          ! g2 __\__ f1
      call ACCEL   (g2,f1,lg,l1,ax2)          ! g1 __.__ f 
      
      is = (-1)**iSpin


      r = r + (ar1*ove2+ar2*ove1*is)*C*ang/sqrt(2.)
      v = v + (av1*ove2+av2*ove1*is)*C*ang/sqrt(2.)                    
      x = x + (ax1*ove2+ax2*ove1*is)*C*ang/sqrt(2.)                    

!      write(6,'(6(A,1p,E11.4))')'C=',C,' ANG=', ang, ' OVERLAP=', OVE,
!     :   ' L', r,
!     :   ' V', v,
!     :   ' A', x

      return
      end
      
      



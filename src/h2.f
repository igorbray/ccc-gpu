      subroutine H2 (rmax, nr, grid, gamma, r0)

C H2 molecule ground state

C INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
C          nr     - the number of points in r-grid
C          grid   - r-grid to calculate wavefunctions on it.

      
      include 'par.f'
      include 'paratom.f'
      PARAMETER (Ncf = 20)
      
      common/meshrr/ meshr,rmesh(maxr,3)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
     

      COMMON /ms/     wm(maxr, 0:Ncf, Ncf)
      COMMON /CIM1/   CM(0:Ncf, 0:Ncf, Ncf, Ncf)
      COMMON /CIM2/   Nc, jmax, nmin(0:Ncf), nmax(0:Ncf)
      COMMON /CIM3/   maxNR(0:Ncf, Ncf)
      dimension A(0:Ncf, 0:Ncf, Ncf, Ncf),
     >          B(0:Ncf, 0:Ncf, Ncf, Ncf),
     >         en(0:Ncf,Ncf),z(0:Ncf,Ncf)

      CHARACTER*4 a4
      CHARACTER*5 a5
      character *25  filnam
      
      character*1 hh(0:6)       !Spectroscopic labels
      data hh /'s','p','d','f','g','h','i'/

       
c  ---------------------------------------
c  Read information about the ground state
c  ---------------------------------------
       open(5,file='mol',STATUS='UNKNOWN')
       read (5,*)   IZ, iSpin, Nc, Etot
       write(6,101) IZ, iSpin, Nc, Etot
 101   FORMAT (////12x,' GROUND STATE OF H2 MOLECULE',
     :           /12x,' ----------------------------'/,
     :           /,   ' Nucleus charge           ',I3,
     :           /,   ' Spin                     ',I3,
     :           /,   ' Number of configurations ',I3,
     :           /,   ' Ground state energy, au  ',F10.6//)
       
       read(5,'(a25)') filnam

c     IZ - nuclear charge
c     Nc - number of two-electron configurations (nl)^2 in the GS
c     nmax(l) - maximum principle quantum number with orbital momentum l
c     CI(L, n) - contribution of the configuration (nl)^2 to the GS
       
c  Initialize CI variables
       
       do l = 0, Nc
          do l1 = 0, Nc
             do n = 1, Nc
                do n1 = 1, Nc
                   CM(l, l1, n, n1) = 0.
                end do
             end do
          end do
       end do
       
       jmax = 0
       S = 0
       do l = 0, Ncf
          nmax(l) = 0
          nmin(l) = Ncf
       end do
       
       meta = 0 !Ground state is assumed
       do i = 1, Nc
          read (5, 13) a5, x,x1,x2,x3,x4
          write(6, 13) a5, x,x1,x2,x3,x4
          call longLABEL(a5(1:3), n, l )
          call     LABEL(a5(4:5), n1,l1)

!          if(n.eq.n1 .and. l.eq.l1) STOP 'WRONG MOLECULAR BASIS'
!          if(n.eq.n1 .and. l.eq.l1) print*, 'WRONG MOLECULAR BASIS'
          if(l.gt.jmax.   or.l1.gt. jmax)    jmax    = max(l,l1)
          if(l.eq.l1) then
             if(n.gt.nmax(l).or.n1.gt.nmax(l1))     nmax(l) = max(n,n1)
            if(n.lt.nmin(l).or.n1.lt.nmin(l1))     nmin(l) = min(n,n1)
          else
             if(n. gt. nmax(l ))      nmax(l)  = n
             if(n1.gt. nmax(l1))      nmax(l1) = n1
             if(n. lt. nmin(l ))      nmin(l)  = n
             if(n1.lt. nmin(l1))      nmin(l1) = n1
          end if

          S = S + x**2
          CM(l,l1,n,n1) = x
          Z (l,n) = x1
          en(l,n) = x2
          Z (l1,n1) = x3
          en(l1,n1) = x4
       end do

       write(6, 11) 'Total', S
 11    FORMAT(A4, F20.8,4F13.7)
 13    FORMAT(A5, F20.8,4F13.7)
       
       print'(/A//A)', 'SLATER BASIS','nl      Zeta          En'
       do l=0,jmax
          do n=  nmin(l), nmax(l)
             do l1=0,jmax
                do n1=  nmin(l1), nmax(l1)
                   if(CM(l,l1,n,n1).ne.0) then
                      write(6,12) n,hh(l),z(l,n),en(l,n)
                      write(6,12) n1,hh(l1),z(l1,n1),en(l1,n1)
                   end if
                end do
             end do
          end do
       end do
! 12    FORMAT(I1,A1,2F12.6)
 12    FORMAT(I2,A1,2F12.6)
       
                          
       do i = 1, meshr
          rt = RMESH(i,1)
          do l=0,jmax
             do n=  nmin(l), nmax(l)
                do l1=0,jmax
                   do n1=  nmin(l1), nmax(l1)
                      if(CM(l,l1,n,n1).ne.0) then
                         tmp  = PR(en(l,n)  ,z(l,n)  ,rt)
                         tmp1 = PR(en(l1,n1),z(l1,n1),rt)
                         wm(i,l,n)   = tmp
                         wm(i,l1,n1) = tmp1
                         write(10,'(10E13.4)') rt,tmp,tmp1
                         phim = amax1(tmp,tmp1)
                         if(phim.gt.expcut)  maxphi=i
                         call OVER(wm(1,l,n),wm(1,l1,n1),A(l,l1,n,n1))
                         tn = 0.5*(en(l,n)+en(l1,n1))
                         tz = 0.5*( z(l,n)+z(l1,n1))
                         B(l,l1,n,n1) = PA(en(l,n)  ,z(l,n))*
     >                        PA(en(l1,n1),z(l1,n1))/PA(tn,tz)**2
                      end if
                   end do
                end do
             end do
          end do
       end do

       print'(/A/10x,A/)', 'MOLECULAR BASIS','CI           OVERLAPS'
       do l=0,jmax
          do n=  nmin(l), nmax(l)
             do l1=0,jmax
                do n1=  nmin(l1), nmax(l1)
                   if(CM(l,l1,n,n1).ne.0) then
!                      write(6,'(2(I1,A1),3F11.6)')
                      write(6,'(I2,A,I1,A1,5F11.6)')
     >                   n,hh(l),n1,hh(l1),CM(l,l1,n,n1),
     >                   A(l,l1,n,n1), B(l,l1,n,n1)
     >                   , PA(en(l,n),z(l,n)), PA(en(l1,n1),z(l1,n1))
                      if(l.eq.l1 .and. n.ne.n1 .and. NC.ne.22)
     >                   CM(l,l1,n,n1)=
     >                     CM(l,l1,n,n1)/(1+A(l,l1,n,n1)**2)**0.5
                   end if
                end do
             end do
          end do
       end do

c Only for MOL22
c Off-center overlap with hydrogen 1s

       if(NC. ne. 22) goto 222
       
       S1 = 0.0
       do n=  2, nmax(0)
          S1 = S1 + A(0,0,n,1)*CM(0,0,n,1)
       end do
       print'(A,F13.6)', 'Overlap with 1s', S1

c Barnett-Coulson function

       do i = 1, meshr
          rt = RMESH(i,1)
          do l=0,jmax
             S = 0.0
             do n=  nmin(l), nmax(l)
                S = S + wm(i,l,n)*CM(l,0,n,1)/rt/2*hat(l)
                if(i.eq.1 .and. CM(l,0,n,1) .ne.0)
     >             print'(I2,A1,F13.6)',n,hh(l),CM(l,0,n,1)
             end do
             write(400+l,'(2E13.4)') rt, S
          end do
       end do

c Normalization check
       
       S = 0.0
       do l=0,jmax
          do n=  max(2,nmin(l)), nmax(l)
             do n1=  max(2,nmin(l)), nmax(l)
                if( CM(l,0,n,1).ne.0 .and. CM(l,0,n1,1) .ne.0) then
                   tn = 0.5*(en(l,n)+en(l,n1))
                   tz = 0.5*( z(l,n)+z(l,n1))
                   B(l,l,n,n1) = PA(en(l,n)  ,z(l,n))*
     >                PA(en(l,n1),z(l,n1))/PA(tn,tz)**2
                   S = S + B(l,l,n,n1)*CM(l,0,n,1)*CM(l,0,n1,1)
                   print'(2(I2,A1),4F13.6)',n,hh(l),n1,hh(l),
     >                CM(l,0,n,1),CM(l,0,n1,1),B(l,l,n,n1)
                end if
             end do
          end do
       end do
       print'(A,F13.6)', 'Normalization', S

c CI modification

       do l=0,jmax
          do n=  nmin(l), nmax(l)
             if( CM(l,0,n,1).ne.0)
     >          CM(l,0,n,1) = CM(l,0,n,1)/(1.0+S1**2)
          end do
       end do
      
 222   return
       end

*********************************************************************
      FUNCTION PR(en,ze,r)
      external S14AAF
      REAL*8 S14AAF
      DOUBLE PRECISION X,A
      integer IFAIL
      real ze,en,r
      
! Normalization
      
!      A = (2.0*ze)**(en+0.5)
!      X = dble(2.0*en+1.0)
!      B = S14AAF(X,IFAIL)**(-0.5) !Gamma-function
!      A = A*B
      A = (2.0*ze)**(0.5*en+0.25)
      X = dble(2.0*en+1.0)
      b = fact(nint(x-1d0))**(-0.25)
c$$$  B = S14AAF(X,IFAIL)**(-0.25) !Gamma-function
      A = A*B
      A = A**2
!      PR = A * r**en * exp(-ze*r)
      tmp =   r**(0.5*en) * exp(-0.5*ze*r)
      PR = A * tmp**2
      RETURN
      END

*********************************************************************
      REAL*8 FUNCTION PA(en,ze)
      include 'par.f'
      external S14AAF
      REAL*8 S14AAF
      DOUBLE PRECISION X,B
      
! Normalization

      
!      A = (2.0*ze)**(en+0.5)
!      X = 2.0d0*en+1.0d0
!      B = S14AAF(X,IFAIL)**(-0.5)
!      PA = A*B
      A = (2.0*ze)**(0.5*en+0.25)
      X = 2.0d0*en+1.0d0
      b = fact(nint(x-1d0))**(-0.25)
c$$$      B = S14AAF(X,IFAIL)**(-0.25)
      PA = (A*B)**2
      RETURN
      END


*********************************************************************
      subroutine H2SAT(nchtop)
                                                                      
!    This subroutine calculates relative satellite intensities I(n)   
!    in hydrogen molecule with an arbitrary ground state                                         
      include 'par.f'
      include 'paratom.f'
      PARAMETER (Ncf = 20)

      common /meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym

      COMMON /ms/     wm(maxr, 0:Ncf, Ncf)
      COMMON /CIM1/   CM(0:Ncf, 0:Ncf, Ncf, Ncf)
      COMMON /CIM2/   Nc, jmax, nmin(0:Ncf), nmax(0:Ncf)
      COMMON /CIM3/   maxNR(0:Ncf, Ncf)

      dimension rint(Ncf)

      dimension Psi(maxr)
      character*1 hh(0:3)         !Spectroscopic labels
      data hh /'s','p','d','f'/

* Total cross-section

      n5=5

      do i=1,meshr
         Psi(i)= 0.
         do l1=0,0!jmax
            do n1=nmin(l1),nmax(l1)
               do l2=0,0!jmax
                  do n2=nmin(l2),nmax(l2)
                     if(CM(l1,l2,n1,n2).ne.0.)then
                        wr1 = wm(n5,l1,n1) 
                        wr2 = wm(n5,l2,n2)
                           if(l1.ne.l2 .or. n1.ne.n2) then
                              tmp=(wm(i,l1,n1)*wr2+wm(i,l2,n2)*wr1)
     >                           /sqrt(2.)
                           else
!                              tmp =  wm(i,l2,n2)*wr1
                              tmp =  wm(i,l1,n1)*wr1
                           end if
                     end if
                     Psi(i)= Psi(i) + tmp*CM(l1,l2,n1,n2)
                  end do
               end do
            end do
         end do
         write(141,'(2E13.4)') rmesh(i,1),Psi(i)
      end do

      call OVER(Psi, Psi, Sum)
!      print'(/A,3E11.4)', 'SUM', sum

* Partial ns cross-sections

      nn=0
      S = 0.
      do nch = 1, nchtop
         call  getchinfo (nch, nt, lg, psi, maxpsi, ea, la, na, l)
         if(la.eq.0 .and. na.gt.nn .and. ea .lt. 0.)nn=na 
         
!         write(6,1954) nch,Ea,na,hh(la),hh(l)
 1954    format(//11x,'CHANNEL ',I2/,
     >   11x,'Ion energy',F7.3/,
     >   11x,'Bound electron',I2,A1/,
     >   11x,'Free  electron k', A1)
      
         if(la.eq.0 .and. ea .lt. 0.)then !Bound s-states

*  Coulomb integrals

            ram = 0.
            do l1=0,0!jmax
               do n1=  nmin(l1), nmax(l1)
                  do l2=0,0!jmax
                     do n2=  nmin(l2), nmax(l2)
                        if(CM(l1,l2,n1,n2).ne.0) then

                           call OVER(psi, wm(1,l1,n1), ver1)
                           call OVER(psi, wm(1,l2,n2), ver2)
                           wr1 = wm(n5,l1,n1) 
                           wr2 = wm(n5,l2,n2) 
                           if(l1.ne.l2 .or. n1.ne.n2) then
                              tmp = (wr1*ver2 + wr2*ver1)/sqrt(2.)
                           else
!                              tmp =  wr2*ver1
                              tmp =  wr1*ver1
                           end if
                           ram = ram + tmp*CM(l1,l2,n1,n2)
                        end if
                        rint(na) = ram**2
                     end do
                  end do
               end do
            end do
            S = S + rint(na) 
         end if
      end do
      ratio = (sum - S)/S
      
      do n = 1, nn
         rint(n) = rint(n)/S
      end do

!      print'(/A,3E11.4)', 'SUM, INTEGRAL', S,sum,sum-S
      write(6,110)
      write(6,'(/4x,10(I2,A,6x),A)//') (i,hh(0), i =1,10),'oo'
      write(6,'(/20F9.4)//') (rint(n)*100,n=1,10),ratio*100
      
 20   FORMAT (2  F17.7)
 110  FORMAT(//10x,' SATELLITE INTENSITIES IN HYDROGEN MOLECULE'/,    
     :         10x,' ------------------------------------------') 
 120  FORMAT ( 20(F8.4,2x))


      
      END      

*********************************************************************
      subroutine hDIPOLE(w1,w,j1,j2,J,r,v,x) 
                                                       
      include 'par.f'                                            
      PARAMETER (Ncf = 20)
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr)
      DIMENSION   f1(maxr), f2(maxr), f3(maxr), f4(maxr)

      COMMON /ms/     wm(maxr, 0:Ncf, Ncf)
      COMMON /CIM1/   CM(0:Ncf, 0:Ncf, Ncf, Ncf)
      COMMON /CIM2/   Nc, jmax, nmin(0:Ncf), nmax(0:Ncf)
      COMMON /CIM3/   maxNR(0:Ncf, Ncf)

      character*1 hh(0:3)         !Spectroscopic labels
      data hh /'s','p','d','f'/


C Get rid of the Simpson's weights for continuum WF
C But don't corrupt them for later use. That's why w2, not w

!      do i=1,maxr
!         if(rmesh(i,3) .ne. 0d0) w2(i) = w(i)/rmesh(i,3)
!      end do
      do i=1,meshr
         w2(i) = w(i)/rmesh(i,3)
      end do

      r = 0.0
      v = 0.0
      x = 0.0

C Angular momentum loop. Always arrange:     
                                           
      do l1=0,jmax                         
         do n1=nmin(l1),nmax(l1)    
            do l2=0,jmax
               do n2=nmin(l2),nmax(l2)
                  anorm = 2.0
                  if(l1.ne.l2 .or. n1.ne.n2) anorm = sqrt(2.)
!                  anorm = 1.0 !No exchange
                  C = CM(l1,l2,n1,n2)/anorm
                  if(C.ne.0.)then
                     if(l2.eq.j2 .and. triang(l1,1,j1).eq.1) then       !     \         
                        call h2DIPOLE(wm(1,l1,n1),wm(1,l2,n2),w1,w2,    !l1  __\__ j1   
     :                     l1,l2,j1,j2,J,C,r,v,x)                       !l2  __.__ j2   
                     end if
                     if(l2.eq.j1 .and. triang(l1,1,j2).eq.1) then       !     \         
                        call h2DIPOLE(wm(1,l1,n1),wm(1,l2,n2),w2,w1,    !l1  __\__ j2   
     :                     l1,l2,j2,j1,J,C,r,v,x)                       !l2  __.__ j1   
                     end if               
                     if(l1.eq.j2 .and. triang(l2,1,j1).eq.1) then       !     \         
                        call h2DIPOLE(wm(1,l2,n2),wm(1,l1,n1),w1,w2,    !l2  __\__ j1   
     :                     l2,l1,j1,j2,J,C,r,v,x)                       !l1  __.__ j2   
                     end if
                     if(l1.eq.j1 .and. triang(l2,1,j2).eq.1) then       !     \          
                        call h2DIPOLE(wm(1,l2,n2),wm(1,l1,n1),w2,w1,    !l2  __\__ j2   
     :                     l2,l1,j2,j1,J,C,r,v,x)                       !l1  __.__ j1   
                     end if               
                  end if
               end do
            end do
         end do
      end do


      
      RETURN
      END
*********************************************************************
      subroutine h2DIPOLE(f1,f2,f3,f4,l1,l2,l3,l4,J,c,r,v,x) 

      include 'par.f'                                            
      DIMENSION   f1(maxr), f2(maxr), f3(maxr), f4(maxr)
      real*8 six,threej
      
!      J0 = l1+l2
!      Jmax = max(J,J0)
!      ang = sqrt(1.0*Jmax) * (-1)**Jmax  


C Polarization of light
      MP=1
      J0 = l1+l2

C Atomic ground state
!      J0 = 0
      
      call IOT3(J0,J,1,0,-MP,MP,threej)
      ang = hat(J)*hat(J0)*threej
      ang =-ang * sqrt(3.0) !To comply with atomic angular part

C Angular momentum check

      if(l2.ne.l4 .or. triang(l1,1,l3).ne.1) then
         print'(/A,I2)', 'l1=',l1    
         print'(A,I2) ', 'l2=',l2     
         print'(A,I2) ', 'l3=',l3   
         print'(A,I2/)', 'l4=',l4   
        STOP 'WRONG ANGULAR MOMENTA'
      end if
         
      call OVER(f2,f4,ove)                       
      call LENGTH  (f1,f3,l1,l3,ar)  !    \      
      call VELOCITY(f1,f3,l1,l3,av)  ! 1 __\__ 3
      call ACCEL   (f1,f3,l1,l3,ax)  ! 2 __.__ 4 


      call iot6(J0,J,1, l3,l1,l4,six)


      ang = ang * (-1)**l1 * six

      r = r + ar*ove*C*ang
      v = v + av*ove*C*ang               
      x = x + ax*ove*C*ang
      
!      print'(4I3,10F9.4)',l1,l2,l3,l4,six,ang,C,ar,ove,r
         

!      write(6,'(6(A,1p,E11.4))')'C=',C,' ANG=', ang, ' OVERLAP=', OVE,
!     :   ' L', r,
!     :   ' V', v,
!     :   ' A', x

      return
      end
      
**********************************************************************
      subroutine H2m 

      
      include 'par.f'
      include 'paratom.f'
      PARAMETER (Ncf = 50)
      
      common/meshrr/ meshr,rmesh(maxr,3)
      common /flogs/ faclog(1000)
     

      COMMON /ms/     wm(maxr, 0:Ncf, Ncf)
      COMMON /CIM1/   CM(0:Ncf, 0:Ncf, Ncf, Ncf)
      COMMON /CIM2/   Nc, Lm, nmin(0:Ncf), nmax(0:Ncf)
      COMMON /CIM3/   maxNR(0:Ncf, Ncf)
      COMMON /CIM4/   jmin(0:Ncf, Ncf), jmax(0:Ncf, Ncf)
      dimension dn(Ncf)

      real*8 faclog,dn,c,et2,x,pl
      CHARACTER*4 a4
      character *25  filnam
      
      character*1 hh(0:6)       !Spectroscopic labels
      data hh /'s','p','d','f','g','h','i'/

c  Read CI file 
      open(5,file='mol38',STATUS='UNKNOWN')
      read (5,*)   IZ, iSpin, Nc, Etot
      write(6,101) IZ, iSpin, Nc, Etot
 101  FORMAT (////12x,' GROUND STATE OF H2 MOLECULE',
     :   /12x,' ----------------------------'/,
     :   /,   ' Nucleus charge           ',I3,
     :   /,   ' Spin                     ',I3,
     :   /,   ' Number of configurations ',I3,
     :   /,   ' Ground state energy, au  ',F10.6//)
      
      read(5,'(a25)') filnam

c  Initialize CI variables
      do l = 0, Nc
         do l1 = 0, Nc
            do n = 1, Nc
               do n1 = 1, Nc
                  CM(l, l1, n, n1) = 0.
               end do
            end do
         end do
      end do
      
      Lm = 0
      S = 0
      do l = 0, Ncf
         nmax(l) = 0
         nmin(l) = Ncf
      end do
      
      do i = 1, Nc
         read (5, 11) a4, x
         write(6, 11) a4, x
         call label(a4(1:2), n, l )
         call label(a4(3:4), n1,l1)
         if(l.gt.Lm.   or.l1.gt. Lm)    Lm    = max(l,l1)
         if(l.eq.l1) then
            if(n.gt.nmax(l).or.n1.gt.nmax(l1))     nmax(l) = max(n,n1)
            if(n.lt.nmin(l).or.n1.lt.nmin(l1))     nmin(l) = min(n,n1)
         else
            if(n. gt. nmax(l ))      nmax(l)  = n
            if(n1.gt. nmax(l1))      nmax(l1) = n1
            if(n. lt. nmin(l ))      nmin(l)  = n
            if(n1.lt. nmin(l1))      nmin(l1) = n1
         end if
         S = S + x**2
         CM(l,l1,n,n1) = x
         if(l.eq.1 .and. l1.eq.1) CM(l,l1,n,n1) = -x/sqrt(3.0)
      end do
      write(6, 11) 'Total', S
 11   FORMAT(A4, F20.8,4F13.7)
       
      print'(/13x,A/,A/)','LAGUERRE BASIS',
     >   'nl     rmin  rmax    N O R M A     Orthogonality'


c Laguerre basis      
      eta = 2.9/1.4             !Hagstrom & Shull
      epscut = 1e-10
      
      do l=0,Lm
         do n=  nmin(l), nmax(l)
            if(l.eq.4)  eta =  8.7/1.4
            if(l.eq.6)  eta = 11.6/1.4
!            eta = 1./n         !Hydrogen test
!            eta = 2./n         !Helium test

            et2 = 2.0d0 * eta
            c = dsqrt(et2)**3
            dn(n) = c*dexp(0.5d0*faclog(n-l) - 1.5d0*faclog(n+l+2))

C Define Jmin(l,n)
            xfirst = exp(log(epscut/dn(nmin(l)))/real(l+1))
            rfirst = xfirst / et2
            i = 1
            do while ((rmesh(i,1).le.rfirst) .and. (i.lt.meshr))
               i = i+1
            end do 
            ifirst = i
            jmin(l,n) = ifirst
            jmax(l,n) = meshr

C Loop by  r-grid
            do i = 1, meshr
               r = rmesh(i,1)
               x = et2 * r
               c = exp(-0.5 * x  +  real(L) * log(x))
            
               call LAG(n+l+1,2*l+2,x,pl)
               wm(i,l,n) = 0.0
               if (i .le. jmax(l,n)) then
                  wm(i,l,n) = c  * real(dn(n)) * pl * r
                  if (i .gt. jmin(l,n)+3+100) then
                     if (     (abs(wm(i-3,l,n))  .lt.  epscut)
     >                  .and. (abs(wm(i-2,l,n))  .lt.  epscut)
     >                  .and. (abs(wm(i-1,l,n))  .lt.  epscut)
     >                  .and. (abs(wm(i,  l,n))  .lt.  epscut))
     >                  jmax(l,n) = i - 3
                  end if 
               end if 
               write(100+10*l+n,'(10E17.7)')  rmesh(i,1),wm(i,l,n),pl
            end do              !End r
            call OVER(wm(1,l,n),wm(1,l,n),ov)
            call OVER1(jmax(l,n),wm(1,l,n),wm(1,l,n),ov1)
            if(n.gt.nmin(l)) then
               call OVER(wm(1,l,n),wm(1,l,n-1),ov2)
               print'(I1,A1, 2I7,2F9.4,E13.4)',n,hh(l),
     >            jmin(l,n),jmax(l,n),ov,ov1,ov2
            else 
               print'(I1,A1, 2I7,2F9.4)',n,hh(l),
     >            jmin(l,n),jmax(l,n),ov,ov1
            end if
         end do                 !End n

      end do                    !End L

      return
      end

**********************************************************************
      subroutine LAG(n,m,x,cx)

!  Associated Laguerre polynomials
!  Recursion:
!    if N < M:    L(N,M)   = 0
!    if N = M:    L(N,M)   = (-1)**M * M! 
!    if N = M+1:  L(N,M) = (M+1)*(M+1-X)*L(M,M)
!    if N = M+2 or greater:
!                 L(N,M)  = -[(X+M-2*N+1)*L(N-1,M) + (N-1)*(N-1)*L(N-2,M)]*N /(N-M)
!  Input parameters:
!    Integer N, the highest order polynomial to compute.
!    Integer M, the parameter.  M must be nonnegative.
!    Real X, the point at which the polynomials are to be evaluated.

!  Output
!    Real CX(0:N), the values of the first N+1 Laguerre polynomials at the point X.
!    Note that polynomials 0 through N will be computed.
!  Source
!    http://www.psc.edu/~burkardt/src/polpak/polpak.f90
      
      common /flogs/ faclog(1000)

      integer n
      real*8 pl(0:n) !work array
      integer i
      integer ifact
      integer m
      real*8 x,cx,faclog

      if ( m .lt. 0 ) stop 'M must be nonnegative'
      if ( n .le. 0 ) stop 'N must be positive'

      do i=0,m-1
         pl(i) = 0.0E+00
      end do
      
      pl(m) = dexp(faclog(m+1))*(-1)**m
      pl(m+1) = dble ( m + 1 ) * ( dble ( m + 1 ) - x ) * pl(m)
      
      do i = m+2, n
         pl(i) = -  dble(i)/dble(i-m) * 
     >      ((x+dble(m-2*i+1))*pl(i-1)+ dble(i-1)**2*pl(i-2))
      end do
      cx = pl(n)
      
      return
      end
*********************************************************************
      subroutine H2SATm(nchtop)
                                                                      
!    This subroutine calculates relative satellite intensities I(n)   
!    in hydrogen molecule with an arbitrary ground state                                         
      include 'par.f'
      include 'paratom.f'
      PARAMETER (Ncf = 50)

      common /meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym

      COMMON /ms/     wm(maxr, 0:Ncf, Ncf)
      COMMON /CIM1/   CM(0:Ncf, 0:Ncf, Ncf, Ncf)
      COMMON /CIM2/   Nc, jmax, nmin(0:Ncf), nmax(0:Ncf)
      COMMON /CIM3/   maxNR(0:Ncf, Ncf)

      dimension rint(Ncf)

      dimension Psi(maxr)
      character*1 hh(0:3)         !Spectroscopic labels
      data hh /'s','p','d','f'/

* Total cross-section

      n5=1
      do i=1,meshr
         Psi(i)= 0.
         do l1=0,0
            do n1=nmin(l1),nmax(l1)
               do l2=0,0
                  do n2=nmin(l2),nmax(l2)
                     if(CM(l1,l2,n1,n2).ne.0.)then
                        wr1 = wm(n5,l1,n1)
                        wr2 = wm(n5,l2,n2)
                        anorm = sqrt(2.)
                        if(l1.eq.l2 .and. n1.eq.n2) anorm = 2.0
                        tmp=(wm(i,l1,n1)*wr2+wm(i,l2,n2)*wr1)/anorm
                     end if
                     Psi(i)= Psi(i) + tmp*CM(l1,l2,n1,n2)
                  end do
               end do
            end do
         end do
         write(41,'(2E13.4)') rmesh(i,1),Psi(i)
      end do

      call OVER(Psi, Psi, Sum)

* Partial ns cross-sections

      nn=0
      S = 0.
      do nch = 1, nchtop
         call  getchinfo (nch, nt, lg, psi, maxpsi, ea, la, na, l)
         if(la.eq.0 .and. na.gt.nn .and. ea .lt. 0.)nn=na 
         
!         write(6,1954) nch,Ea,na,hh(la),hh(l)
 1954    format(//11x,'CHANNEL ',I2/,
     >   11x,'Ion energy',F7.3/,
     >   11x,'Bound electron',I2,A1/,
     >   11x,'Free  electron k', A1)
      
         if(la.eq.0 .and. ea .lt. 0.)then !Bound s-states

*  Coulomb integrals

            ram = 0.
            do l1=0,0
               do n1=  nmin(l1), nmax(l1)
                  do l2=0,0
                     do n2=  nmin(l2), nmax(l2)
                        if(CM(l1,l2,n1,n2).ne.0) then
                           call OVER(psi, wm(1,l1,n1), ver1)
                           call OVER(psi, wm(1,l2,n2), ver2)
                           wr1 = wm(n5,l1,n1) 
                           wr2 = wm(n5,l2,n2) 
                           anorm = sqrt(2.)
                           if(l1.eq.l2 .and. n1.eq.n2) anorm = 2.0
                           tmp = (wr1*ver2 + wr2*ver1)/anorm
                           ram = ram + tmp*CM(l1,l2,n1,n2)
                        end if
                        rint(na) = ram**2
                     end do
                  end do
               end do
            end do
            S = S + rint(na) 
         end if
      end do
      ratio = (sum - S)/S
      
      do n = 1, nn
         rint(n) = rint(n)/S
      end do

      print'(/A,3E11.4)', 'SUM, INTEGRAL', S,sum,sum-S
      print'(A,E11.4)'  , 'Rfirst       ', rmesh(n5,1)
      write(6,110)
      write(6,'(/4x,12(I2,A1,6x))//')
     >   (i,'s', i =1,NN)
      write(6,120)
     >         (rint(n)*100,n=1,NN),ratio*100
      
 20   FORMAT (2  F17.7)
 110  FORMAT(//10x,' SATELLITE INTENSITIES IN HYDROGEN MOLECULE'/,    
     :         10x,' ------------------------------------------') 
 120  FORMAT ( 12(F7.4,2x))

      END      

*********************************************************************
      subroutine hDIPOLEm(w1,w,j1,j2,J,r,v,x) 
                                                       
      include 'par.f'                                            
      PARAMETER (Ncf = 50)
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr)
      DIMENSION   f1(maxr), f2(maxr), f3(maxr), f4(maxr)

      COMMON /ms/     wm(maxr, 0:Ncf, Ncf)
      COMMON /CIM1/   CM(0:Ncf, 0:Ncf, Ncf, Ncf)
      COMMON /CIM2/   Nc, jmax, nmin(0:Ncf), nmax(0:Ncf)
      COMMON /CIM3/   maxNR(0:Ncf, Ncf)

C Get rid of the Simpson's weights for continuum WF
C But don't corrupt them for later use. That's why w2, not w

!      do i=1,maxr
!         if(rmesh(i,3) .ne. 0d0) w2(i) = w(i)/rmesh(i,3)
!      end do
      do i=1,meshr
         w2(i) = w(i)/rmesh(i,3)
      end do

      r = 0.0
      v = 0.0
      x = 0.0

C Angular momentum loop. Always arrange:   !    \        
                                           !1  __\__ 3   
      do l1=0,jmax                         !2  __.__ 4   
         do n1=nmin(l1),nmax(l1)    
            do l2=0,jmax
               do n2=nmin(l2),nmax(l2)
                  anorm = sqrt(2.)
                  if(l1.eq.l2 .and. n1.eq.n2) anorm = 2.0
                  C = CM(l1,l2,n1,n2)/anorm
                  if(C.ne.0.)then
                     if(l2.eq.j2 .and. triang(l1,1,j1).eq.1) then
                        call h2DIPOLE(wm(1,l1,n1),wm(1,l2,n2),w1,w2,
     :                     l1,l2,j1,j2,J,C,r,v,x)
                     end if
                     if(l2.eq.j1 .and. triang(l1,1,j2).eq.1) then
                        call h2DIPOLE(wm(1,l1,n1),wm(1,l2,n2),w2,w1,
     :                     l1,l2,j2,j1,J,C,r,v,x)
                     end if               
                     if(l1.eq.j2 .and. triang(l2,1,j1).eq.1) then
                        call h2DIPOLE(wm(1,l2,n2),wm(1,l1,n1),w1,w2,
     :                     l2,l1,j1,j2,J,C,r,v,x)
                     end if
                     if(l1.eq.j1 .and. triang(l2,1,j2).eq.1) then
                        call h2DIPOLE(wm(1,l2,n2),wm(1,l1,n1),w2,w1,
     :                     l2,l1,j2,j1,J,C,r,v,x)
                     end if               
                  end if
               end do
            end do
         end do
      end do


      
      RETURN
      END


      subroutine longLABEL(aa, n, l)
      
      character*3 aa

      if(aa(2:2) .eq. '0') n = 0
      if(aa(2:2) .eq. '1') n = 1
      if(aa(2:2) .eq. '2') n = 2
      if(aa(2:2) .eq. '3') n = 3
      if(aa(2:2) .eq. '4') n = 4
      if(aa(2:2) .eq. '5') n = 5
      if(aa(2:2) .eq. '6') n = 6
      if(aa(2:2) .eq. '7') n = 7
      if(aa(2:2) .eq. '8') n = 8
      if(aa(2:2) .eq. '9') n = 9
      if(aa(1:2) .eq. '10') n = 10
      if(aa(1:2) .eq. '11') n = 11
      if(aa(1:2) .eq. '12') n = 12
      if(aa(1:2) .eq. '13') n = 13
      if(aa(1:2) .eq. '14') n = 14
      if(aa(1:2) .eq. '15') n = 15

      if(aa(3:3) .eq. 's') l = 0
      if(aa(3:3) .eq. 'p') l = 1
      if(aa(3:3) .eq. 'd') l = 2
      if(aa(3:3) .eq. 'f') l = 3
      if(aa(3:3) .eq. 'g') l = 4
      if(aa(3:3) .eq. 'h') l = 5
      if(aa(3:3) .eq. 'i') l = 6

      RETURN
      END

      

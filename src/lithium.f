      subroutine LiDIPOLEc(nst,w2,L,lp,lg,r,v,x) 

C Single-configuration one-electron dipole ME for Li <2s|r|kp>
C Ground state correlation is to be included
      
      use CI_MODULE
      include 'par.f'                                            

      PARAMETER (Ncf = 20)
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr), f(maxr), f1(maxr)
      
      COMMON /cat1/   CJ(0:Ncf, Ncf, Ncf)
      COMMON /cat2/   Nc, nmax(0:Ncf), jmax
      COMMON /cat3/   wa(maxr, 0:Ncf, Ncf)
      COMMON /gstspin/   iSpin, meta

      
      integer la, sa, lpar, np
      common /helium/ la(KNM), sa(KNM), lpar(KNM), np(KNM)
      double precision E
      common /CcoefsE/  E(KNM)

      double precision ortint
      common /ortog/  ortint(nspmax,nspmax)
      integer nspm, lo, ko, nset
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      real fl
      common /funLag/  fl(nmaxr,nspmax)
      integer maxf, minf
      common /minmaxf/ maxf(nspmax),minf(nspmax)

      integer na, nam
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer nicm, ncore
      common /corearray/ nicm, ncore(nspmCI)      

      double precision E_st
      integer L_st, S_st, Par_st

      character*3 chan(knm)
      common /charchan/ chan

      
      integer  jn1, jn2, n1, n2, l1, l2

      character lh(0:10)
c$$$      hat(l)=sqrt(2.0 * l + 1.0)
      data lh / 's','p','d','f','g','h','n','o',
     >          'x','x','x'/

C  Two-electron state number  nst  description
      E_st = E(nst)             ! two-electron energy
      L_st = la(nst)            ! orb. ang.  mom.
      S_st = sa(nst)            ! spin
      Par_st = lpar(nst)        ! parity
      AK = (1+(-1)**S_st)/2     ! deta_{S0}
      AL = (-1)**(S_st+1) * hat(S_st)

c$$$      print'(/A,I3,F9.4,2I3,3x,A)',
c$$$     >   'State', nst, E_st, L_st, S_st,chan(nst)  

C Get rid of the Simpson's weights for continuum WF
C But don't corrupt them for later use. That's why w2, not w1

!      do i=1,maxr
!         if(rmesh(i,3) .ne. 0d0) w1(i) = w2(i)/rmesh(i,3)
!      end do
      do i=1,meshr
         w1(i) = w2(i)/rmesh(i,3)
      end do

      
C Initialization
      r = 0.
      v = 0.
      x = 0.

      ang1=(-1)**L_st/hat(L_st)
      ang2=   hat(lp)/hat(L_st)/sqrt(3.0)

      
C MCHF loops

      do j=0,jmax
         do n=j+1, nmax(j)
            if(j.eq.0 .and. n.eq.2) CYCLE !Skip 2s in MCHF expansion

                                !    j   ______   
            ang =(-1)**j/hat(j) !(-1) / V(2j+1)

            CM = CJ(j,n,2)!*ang
            if(CM.eq.0) CYCLE
!            print'(I3,A,F9.4)', n,lh(j),  CM

            
C Target states loops            

      do jn1=1,nam(nst)
         n1 = na(jn1,nst)
         l1=lo(n1)
         do jn2=1,nam(nst)
            n2 = na(jn2,nst)
            l2=lo(n2)
            CC = C(nst,jn1,jn2)
            if(CC .eq. 0) cycle

c$$$            print'(A,2(I3,A),F9.4)',
c$$$     >         'Configuration', n1,lh(l1),n2,lh(l2),CC


C Direct matrices
            if(j.eq. l1) then
               call OVER(wa(1,j,n), fl(1,n1), res11) 
C DK-matrix
               if(j.eq.l2 .and. lp.eq.1) then
                  call OVER(wa(1,j,n), fl(1,n2), res12) 

c$$$                  print'(A,I1,A,F9.4)', '<1s|',n1,'s>', res11
c$$$                  print'(A,I1,A,F9.4)', '<1s|',n2,'s>', res12


                  call LENGTH  (wa(1,0,2),w1,0,lp,r2) !<2s|r|kp>       
                  call VELOCITY(wa(1,0,2),w1,0,lp,v2)           
                  call ACCEL   (wa(1,0,2),w1,0,lp,x2)           
                  
                  r = r + AK*res12*r2*res11*CC*CM*ang
                  v = v + AK*res12*v2*res11*CC*CM*ang
                  x = x + AK*res12*x2*res11*CC*CM*ang
               end if
C DL-matrix
               
               if(triang(j,lp,1).eq.1 .and. l2.eq.0) then
                  call OVER(wa(1,0,2), fl(1,n2), res22) 

c$$$                  print'(A,I1,A,F9.4)', '<2s|',n2,'s>', res22
                  
                  call LENGTH  (wa(1,j,n),w1,j,lp,r1) !<nl|r|kL>      
                  call VELOCITY(wa(1,j,n),w1,j,lp,v1)           
                  call ACCEL   (wa(1,j,n),w1,j,lp,x1)

                  r = r + AL*res22*r1*res11*CC*CM*ang
                  v = v + AL*res22*v1*res11*CC*CM*ang
                  x = x + AL*res22*x1*res11*CC*CM*ang
               end if
            end if

C Exchange matrices one
            if(j.eq. l1) then
               call OVER(wa(1,j,n), fl(1,n1), res11) 
C EK-matrix
               if(lp.eq.0 .and. triang(j,l2,1).eq.1) then

                  call OVER(wa(1,0,2), w1, res2k) !<2s|kL> 

                  call LENGTH  (wa(1,j,n),fl(1,n2),j,l2,r12) !<nl|r|n2l2> 
                  call VELOCITY(wa(1,j,n),fl(1,n2),j,l2,v12)           
                  call ACCEL   (wa(1,j,n),fl(1,n2),j,l2,x12)           
               
                  r = r + AK*r12*res2k*res11*CC*CM*ang*2
                  v = v + AK*v12*res2k*res11*CC*CM*ang*2
                  x = x + AK*x12*res2k*res11*CC*CM*ang*2
               end if
C EL2-matrix
               if(lp.eq.j .and. l2.eq.1) then

                  call OVER(wa(1,j,n), w1, res1k) !<nl|kL>

                  call LENGTH  (wa(1,0,2),fl(1,n2),0,1,r22) !<2s|r|n2p> 
                  call VELOCITY(wa(1,0,2),fl(1,n2),0,1,v22)           
                  call ACCEL   (wa(1,0,2),fl(1,n2),0,1,x22)           
                  
                  r = r + AL*r22*res1k*res11*CC*CM*ang
                  v = v + AL*v22*res1k*res11*CC*CM*ang
                  x = x + AL*x22*res1k*res11*CC*CM*ang
               end if
            end if
               
            if(triang(j,l1,1).eq.1) then
               call LENGTH  (wa(1,j,n),fl(1,n1),j,l1,r11) !<nl|r|n1l1> 
               call VELOCITY(wa(1,j,n),fl(1,n1),j,l1,v11)           
               call ACCEL   (wa(1,j,n),fl(1,n1),j,l1,x11)           
C EK2-matrix
               if(j.eq.l2 .and. lp.eq.0) then
                  call OVER(wa(1,j,n), fl(1,n2), res12) 
                  call OVER(wa(1,0,2), w1, res2k) !<2s|kL> 
                  
c$$$                  r = r + AK*res12*res2k*r11*CC*CM*ang
c$$$                  v = v + AK*res12*res2k*v11*CC*CM*ang
c$$$                  x = x + AK*res12*res2k*x11*CC*CM*ang
               end if
C EL-matrix
               if(l2.eq.0 .and. lp.eq.j) then
                  call OVER(wa(1,0,2), fl(1,n2), res22) 
                  call OVER(wa(1,j,n), w1, res1k) !<nl|kL>
                  
                  r = r + AL*res22*res1k*r11*CC*CM*ang
                  v = v + AL*res22*res1k*v11*CC*CM*ang
                  x = x + AL*res22*res1k*x11*CC*CM*ang
               end if
            endif
         enddo
      enddo

      enddo                     !MCHF loops
      enddo

      RETURN
      END

      subroutine LiSATEL(nchtop,lg)

C Asymptotic satellite intensity ratios
C Ground state correlation is to be included
      
      use CI_MODULE
      include 'par.f'                                            

      PARAMETER (Ncf = 20)
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr), f(maxr), f1(maxr)
      
      COMMON /cat1/   CJ(0:Ncf, Ncf, Ncf)
      COMMON /cat2/   Nc, nmax(0:Ncf), jmax
      COMMON /cat3/   wa(maxr, 0:Ncf, Ncf)
      COMMON /gstspin/   iSpin, meta

      character*3 chan(knm), chs(0:lamax), chat
      common /charchan/ chan
      
      integer la, sa, lpar, np
      common /helium/ la(KNM), sa(KNM), lpar(KNM), np(KNM)
      double precision E
      common /CcoefsE/  E(KNM)

      double precision ortint
      common /ortog/  ortint(nspmax,nspmax)
      integer nspm, lo, ko, nset
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      real fl
      common /funLag/  fl(nmaxr,nspmax)
      integer maxf, minf
      common /minmaxf/ maxf(nspmax),minf(nspmax)

      integer na, nam
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer nicm, ncore
      common /corearray/ nicm, ncore(nspmCI)      

      double precision E_st
      integer L_st, S_st, Par_st
      
      integer  jn1, jn2, n1, n2, l1, l2

      character lh(0:10)
      data lh / 's','p','d','f','g','h','n','o',
     >          'x','x','x'/
      character*1 hh(0:3)         !Spectroscopic labels
      data hh /'s','p','d','f'/
      dimension Psi(maxr)

      dimension SAT(0:1,Ncf)

      do n=1,10
         do i=0,1
            SAT(i,n)=0.0
         end do
      end do

      write(6,
     >   '(//,"ASYMPTOTIC SATELLITE INTENSITIES S(ns)/S(tot)"/,
     >   "---------------------------------------------" /,
     >   "  n     Singlet      Triplet"/)')

* Ground state overlap
      
      wrc = 0.0
      do n=1, nmax(0)
         if(n.eq.2) CYCLE       !Skip 2s in MCHF expansion
         CM = CJ(0,n,2)
         if(CM.eq.0) CYCLE
         wr1 = wa(5,0,n)/rmesh(5,1)
         wr2 = wa(5,0,2)/rmesh(5,1)
         wrc =wrc + CM**2*(wr1**2*2 + wr2**2)
      end do
      
c$$$      print'(A,4F9.4)', 'Radial 1s(0),2s(0)', wr1,wr2,wrc

* Partial ns cross-sections

      nn=0
      S = 0.
      do nch = 1, nchtop
         call  getchinfo (nch, nst, lg, psi, maxpsi, ea, laa, nna, l)

         E_st = E(nst)          ! two-electron energy
         L_st = la(nst)         ! orb. ang.  mom.
         S_st = sa(nst)         ! spin

         chat = chan(nst)
         call label(chat(2:3), nc, lc )


         
         if(L_st.eq.0 .and. nc.gt.nn
     >      .and. E_st .lt. 0. .and. S_st.eq.0) nn=nc !Singlet s-channels counter
         
c$$$         write(6,1954) nch,E_st,nna,hh(L_st),hh(l)
c$$$ 1954    format(//11x,'CHANNEL ',I2/,
c$$$     >      11x,'Ion energy',F7.3/,
c$$$     >      11x,'Bound electron',I2,A1/,
c$$$     >      11x,'Free  electron k',A1)
         
      
         if(L_st.eq.0 .and.  E_st .lt. 0.)then !Bound s-states
            
* MCHF loop
            TMP = 0.0
            do n=1, nmax(0)
               if(n.eq.2) CYCLE 
               CM = CJ(0,n,2)
               if(CM.eq.0) CYCLE
c$$$               print'(/A,I1,A1,E13.4)',  'MCHF  ', n, 's',  CM

               wr1 = wa(5,0,n)/rmesh(5,1)

* MCHF loop      
               
               over1 = 0.0
               over2 = 0.0

* Target states loops            

               do jn1=1,nam(nst)
                  n1 = na(jn1,nst)
                  l1=lo(n1)
                  do jn2=1,nam(nst)
                     n2 = na(jn2,nst)
                     l2=lo(n2)
                     CC = C(nst,jn1,jn2)
                     if(CC .eq. 0) cycle

c$$$                     print'(A,2(I3,A),F9.4)',
c$$$     >                  'Configuration', n1,lh(l1),n2,lh(l2),CC

                     if(l1.eq.0 .and. l2.eq.0)  then
                        call OVER(wa(1,0,n), fl(1,n1), res11) 
                        call OVER(wa(1,0,n), fl(1,n2), res12) 
                        call OVER(wa(1,0,2), fl(1,n2), res22)
                        
                        over1 = over1 + res11*res12*CC
                        over2 = over2 + res11*res22*CC*sqrt(2.0)

c$$$                 print'(A1,I1,A,I1,A,F9.4)', '<',n,'s|',n1,'s>', res11
c$$$                 print'(A1,I1,A,I1,A,F9.4)', '<',n,'s|',n2,'s>', res12
c$$$                 print'(A,I1,A,F9.4)', '<2s|',n2,'s>', res22

                     end if
                  enddo
               enddo
c$$$               print'(A,2E13.4)', 'Overlaps',  over1,over2

               print*,'S_st:',S_st
               if(S_st .eq. 0) then
!                  TMP = TMP + CM*(over2*wr1+over1*wr2*sqrt(2.)) !Its a bug?
                  TMP = TMP + CM*(-over2*wr1+over1*wr2*sqrt(2.))
               else
                  TMP = TMP + CM*(over2*wr1)
               end if
               
            end do

            Cns = TMP**2 * (2*S_st+1)/2./wrc
            Sat(S_st,nc) = TMP**2 * (2*S_st+1)/2./wrc
            
            print'(2A,2E13.4)', 'Satel  ', chan(nst),  Cns, E_st
         endif
      enddo

      if(nn.lt.10) print'(//4x,A/)',
     >   'WARNING: Less than 10 target s-states are used!'

      Sum0 = 0.0
      Sum1 = 0.0
      do n=1,9
         print'(I3,2E13.4)', n, Sat(0,n), Sat(1,n)
         Sum0 = Sum0 + Sat(0,n) 
         Sum1 = Sum1 + Sat(1,n)
      end do
      print'(A,2E13.4)', '1-9', Sum0, Sum1
      print'(A,2E13.4)', 'oo',  1.0-Sum0-Sum1      

!      STOP !!!!!!!!!
      RETURN
      END


      subroutine mainhe(Nmax,namax,pnewC,eproj,etot,
     >     lastop,nnbtop,ovlp,phasen,regcut,expcut,ry,enion,enlevel,
     >     enionry,nchanmax,ovlpnl,slowe,ne2e,vdcore)
      use vmat_module !for nodeid only
      include 'par.f'
      integer la, sa, lpar
      common /helium/ la(KNM), sa(KNM), lpar(KNM), np(KNM)
      double precision E, ortint, e1r, Z
      common/hame1/ e1r(nspmax,nspmax)
      common /CcoefsE/ E(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common /Zatom/ Z
      integer pnewC
      common /CIdata/ na_CI(nspmCI,KNM), nam_CI(KNM)

      integer:: na

c     should be uncommented also in osc.f file 
c     common /osc_strength/ oscstr(30,3,2), weight_cont(KNM)
c
C     Igor's stuff
      common /meshrr/ nr, gridr(nmaxr,3)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      character chan(knm)*3, opcl*6, enq*14
      common /charchan/ chan
      common /chanen/ enchan(knm)
      dimension psi(maxr),ovlp(ncmax,0:lamax),
     >   ovlpnl(ncmax,0:lamax,nnmax),slowe(ncmax)
      complex phasen(ncmax,0:lamax),sigc,phase,ovlpnl

      real:: vdcore(maxr,0:lamax)
      logical:: helium, exists, torf
      real:: temp(nr)
c      
      if (nodeid.eq.1) print*, 'MAINHE: nznuc, zasym:', nznuc, zasym
      torf = .false.

      temp(:) = 0.0
      itemp = -1
      if(nznuc - nint(zasym) .eq. 2) then
         helium = .true.
      else 
         helium = .false.
         pi = acos(-1.0)
         tmp = -(3d0/2d0)*(3d0/(4d0*pi*pi))**(1d0/3d0) / 2.0
         do l = labot, latop
            do k = 1,nabot(l)-1
               m1 = 1
               m2 = istoppsinb(k,l) 
               itemp = max(m2, itemp)
               temp(m1:m2) = temp(m1:m2) + 
     >              (4.0*l+2.0)*psinb(m1:m2,k,l)*psinb(m1:m2,k,l)
            enddo
         enddo
         temp(m1:m2) = tmp * (temp(m1:m2))**(1d0/3d0) 
      endif

      rmax=gridr(nr,1)
      gamma = 0.0
      r0 = 1.0
c$$$      Z = dble(nznuc)
      ry = 13.6058
C  Save the PSINB and ENPSINB arrays in Dmitry's arrays, and then redefine
C  the ENPSINB to contain the target state energies.
      do ni = 1, nspmax
         do nf = 1, nspmax
            ortint(nf,ni) = 0d0
            e1r(nf,ni) = 0d0
         enddo
      enddo       
c
      if (nodeid.eq.1) then
         print*,'These s.p.orbitals comes from diagonalization of ',
     >      'the ion Hamiltonian and go to the CI program as ',
     >      'the s.p.basis'
         print*, 'labot, latop:', labot, latop
      endif
      n = 0
      do l = labot, latop
         do k = nabot(l), natop(l)
            n = n + 1
            ko(n) = k - nabot(l) + 1
            lo(n) = l
            nset(n) = 1
            ortint(n,n) = 1d0
            e1r(n,n) = enpsinb(k,l) / 2.0  !a.u.
            minf(n) = 1
            maxf(n) = istoppsinb(k,l)
            if (nodeid.eq.1)
     >         write(*,'(3i5,1P,e12.4)') n,lo(n),ko(n),e1r(n,n)
            do i = minf(n), maxf(n)
               fl(i,n) = psinb(i,k,l)
            enddo 
         enddo
      enddo 
      nspm = n
c$$$c     Check starting and finishing values of one-electron fl(i,n) functions
c$$$      do n=1,nspm
c$$$         call minmaxi(fl(1,n),nr,i1,i2)
c$$$         if(abs(i2 - maxf(n)) .ge. 10)
c$$$     >      print*, '  ',n,i2,' changed to ',maxf(n)
c$$$         maxf(n) = i2
c$$$         minf(n) = i1
c$$$      end do
      
      
      if(nspm.gt.nspmCI) then
         print*, ' nspm > nspmCI, increase nspmCI in par.f'
         stop
      end if

      do l = labot, latop     
         if(natop(l)-nabot(l)+1.gt.komax) then
            print*, ' Warning: natop(l)-nabot(l)+1 > komax, ',
     >           'structure.f : array  Ns(0:lomax,0:1,komax,2)  ',
     >           'can be out of bounds'
            stop
         end if
      end do
      enionry = enpsinb(ninc,linc)


         
c      if (nznuc.eq.2) then
c         if (abs(enionry + 4.0).gt.1e-4) 
c     >      stop 'enpsinb(ninc,linc) should be - 4'
c
c     elseif (nznuc.eq.4) then
c         if (abs(enionry + 1.33).gt.1e-2)          
c     >      stop 'enpsinb(ninc,linc) should be - 1.33'
c      else
c         stop 'Have not done this yet'
c      endif 
      nspmW = max(nspm,10)
      call structure(Nmax,nspmW,namax,pnewC,E,enionry,eproj,vdcore,
     >   slowe(1)) 
c$$$      print*, 'Finish structure'
      lg = 0
      nchanmax = 0
c     This block is for calculation of the osc.str. for continuum states.
c     open(136,file='osc_str_cont')
c      write(136,'("# en(eV)/osc_str(a.u.): l        v")')

      if (ne2e.gt.0) then
         en = abs(slowe(1))/ry / 2.0
         do l = 0, 1
            do is = 0, 1
               do ip = 0, 1
                  write(enq,'(1p,e10.4,"_",3i1)')en*13.6058*2.0,
     >            l,is,ip
                  inquire(file=enq,exist=exists)
                  torf = torf.or.exists
               enddo
            enddo
         enddo
         print*,'Using the overlap ansatz is:',torf
      endif
      

      
      ovlpnl(:,:,:) = (0.0,0.0)
      minc=1
      maxc=nr
      iternum = 200
      accuracy = 1e-4
      do nst = 1, Nmax
c     find fc. orbital in fc. app. for state nst
         nfc_orb = na_CI(2,nst)
c
         call getchinfo (nst, ntmp, lg, psi, maxpsi, ea, latom, na, li)
         enpsinb(na,latom) = E(nst) * 2.0 - enionry
         if (nst.ne.ntmp)
     >      stop 'Expect NST = NTMP in mainhe'
         etot = eproj + enpsinb(ninc,linc)
         if (etot.ge.enpsinb(na,latom)) then
            onshellk = sqrt(etot-enpsinb(na,latom))
            opcl = 'open'
            if (enpsinb(na,latom).gt.0.0) then
c$$$               if (ne2e.ne.0) then
               q = sqrt(enpsinb(na,latom))
               en = q*q/2.0
               call contfunc(helium,temp,vdcore,
     >              latom,sa(nst),en,iternum,accuracy,
     >            fl(1,1),minf(1),maxf(1),psi,minc,maxc,phase,sigc,
     >              fl(1,nfc_orb),minf(nfc_orb),maxf(nfc_orb))
               call e2eovp(nst,nspm,lo,fl,maxf,minf,pnewC,Nmax,namax,
     >            E,en,latom,sa(nst),psi,minc,maxc,tmp)
               ovlp(na,latom) = tmp
               phasen(na,latom) = phase * sigc
c$$$               print*,'Overlap,phase,sigc',tmp,phase,sigc
c
c$$$               do nnn = 1, Nmax
c$$$                  if (ne2e.ne.0) then
c$$$                     if (nnn.eq.1) then
c$$$                        open(42,file=chan(nst))
c$$$                        write(42,'("# en =",f9.3)') en*ry*2.0
c$$$                     endif 
c$$$                  if (la(nnn).eq.latom.and.sa(nnn).eq.sa(nst)
c$$$     >                  .and.lpar(nst).eq.lpar(nnn)) then
c$$$                     call e2eovp(nnn,nspm,lo,fl,maxf,minf,pnewC,Nmax,
c$$$     >                  namax,E,en,latom,sa(nst),psi,minc,maxc,tmp)
c$$$c                     print*,'nst,nnn, overlap ',chan(nst),chan(nnn),tmp
c$$$                     write(*,'(2A3,1X,E9.2,2X)', ADVANCE='NO')
c$$$     >                    chan(nst),chan(nnn),tmp
c$$$                     write(42,*) (E(nnn)*2.0-enionry)*ry,tmp,chan(nnn)
c$$$                  endif
c$$$                  if (nnn.eq.Nmax) close(42)
c$$$                  endif 
c$$$               enddo 
               print*
c
c     This block is for calculation of the osc.str. for continuum states.
c               osc_cont_len = oscstr(nst,1,1) * ovlp(na,latom)**2/q
c               osc_cont_vel = oscstr(nst,1,2) * ovlp(na,latom)**2/q
c               osc_t = oscstr(nst,1,1) / weight_cont(nst)
c               osc_t = oscstr(nst,1,2) / weight_cont(nst)
               
c               write(136,*) q**2*27.2116/2.0,osc_cont_len,osc_cont_vel,
c     >         1.0/weight_cont(nst), ovlp(na,latom)**2/q  
c$$$               endif ! ne2e.ne.0
            else
               if (na.gt.nnbtop) then
                  print*,'Increase NNBTOP to at least NA:',nnbtop,na
                  stop 'Increase NNBTOP'
               endif 
               do nn = nabot(latom), nnbtop
                  if (na.eq.nn) then
                     ovlpnl(na,latom,nn) = 1.0
                  else
                     ovlpnl(na,latom,nn) = 0.0
                  endif
               enddo 
            endif 
            do ne = 1, ne2e * 2
               if (ne.le.ne2e) then
                  en = abs(slowe(ne))/ry / 2.0
               else
                  en = (etot - abs(slowe(ne-ne2e))/ry) / 2.0
               endif 
               if (en.le.0.0) then
                  print*,'Secondary energy <= 0;en,etot,ne',en,etot,ne
                  stop 'Secondary energy <= 0'
               else
                  print*,'en:',en
               endif
               if (lpar(nst).eq.(-1)**latom) then
                  ip = 0
               else
                  ip = 1
               endif 
               write(enq,'(1p,e10.4,"_",3i1)')en*13.6058*2.0,
     >            latom,sa(nst),ip
               inquire(file=enq,exist=exists)
               if (exists) then
                  print*,'reading:',enq
                  call get_e2eovp_gen(enq,nst,nspm,lo,fl,maxf,minf,
     >               namax,en,la(nst),sa(nst),lpar(nst),
     >               ovlpnl(na,latom,nnbtop+ne))
                  print*,'ovlp:',ovlpnl(na,latom,nnbtop+ne)
               else
                  if (torf) then
                     ovlpnl(na,latom,nnbtop+ne) = 0.0
                  else 
                     call contfunc(helium,temp,vdcore,
     >                  latom,sa(nst),en,iternum,accuracy,
     >                  fl(1,1),minf(1),maxf(1),psi,minc,maxc,phase,
     >                  sigc,fl(1,nfc_orb),minf(nfc_orb),maxf(nfc_orb))
                     call e2eovp(nst,nspm,lo,fl,maxf,minf,pnewC,Nmax,
     >                  namax,E,en,latom,sa(nst),psi,minc,maxc,tmp)
                     ovlpnl(na,latom,nnbtop+ne) = tmp * phase
                     print*,'FC ovlp:',ovlpnl(na,latom,nnbtop+ne)
                  endif 
               endif 
            enddo 
         else
            opcl = 'closed'
            onshellk = -sqrt(enpsinb(na,latom)-etot)
            ovlp(na,latom) = 0.0
         endif 
         elevel = (enpsinb(na,latom) * ry + enion) * 
     >      enlevel / enion
         if (nodeid.eq.1)
     >      print '(i3,a4,f12.4," Ry",f12.4," eV",f12.1," / cm",
     >      a8,f9.4,f8.2)', 
     >      nst, chan(nst),enpsinb(na,latom),enpsinb(na,latom)*ry,
     >      elevel, opcl, onshellk, ovlp(na,latom)
         nchanmax = nchanmax + latom + 1
      enddo 
C  EDELTA is used to shift the energy levels
      edelta = - enion - enpsinb(ninc,linc) * ry
      if (nodeid.eq.1) then
         write(*,'("finish structure calculation")')
         write(*,'("**********************************************")')
c
         print*,'Error in ground state energy and ',
     >      'effective incident energy (eV):', edelta, etot * ry + enion
      endif
c      close(4)
c      close(136)
      end




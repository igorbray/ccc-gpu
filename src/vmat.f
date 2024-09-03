
      subroutine initv(vmat,n,nsmax)
      real vmat(n,n+1)
      vmat(:,:) = 0.0
c$$$      do ki = 1, n
c$$$         do kf = ki, n
c$$$            vmat(kf,ki) = 0.0
c$$$            if (nsmax.eq.1) vmat(ki,kf+1) = 0.0
c$$$         end do
c$$$      end do
      return
      end
      
      subroutine initve2e(ve2e,n,nchtope2e)
      real ve2e(nchtope2e,n)
      do ki = 1, n
         do nchf = 1, nchtope2e
            ve2e(nchf,ki) = 0.0
         enddo
      enddo
      return
      end

      subroutine getche2e(nche2e,lg,lfast,lslow,lpsi,l)
      include 'par.f'
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym

      nch = 0
      do l = 0, lfast
         do lpsi = abs(lg-l), min(abs(lg+l), lslow)
            if ((-1)**(lpsi+l).eq.(-1)**(lg+ipar)) then
               nch = nch + 1
               if (nch.eq.nche2e) return
            endif
         enddo 
      enddo
      nche2e = 0
      end

      subroutine getchnl(ch,n,l,nc)
      character ch*3, c

      n = ichar(ch(2:2)) - ichar('0')
      c = ch(1:1)
      if (c .eq. 's') then
         nc = 1
      elseif (c .eq. ' ') then
         nc = 2
      elseif (c .eq. 't') then
         nc = 3
      elseif (c .eq. 'S') then
         nc = 4
      else
         nc = 5
      endif 
      if (ch(3:3).eq.'S') then
         l = 0
      elseif (ch(3:3).eq.'P') then
         l = 1
      elseif (ch(3:3).eq.'D') then
         l = 2
      elseif (ch(3:3).eq.'F') then
         l = 3
      elseif (ch(3:3).eq.'G') then         
         l = 4
      elseif (ch(3:3).eq.'H') then
         l = 5
      elseif (ch(3:3).eq.'I') then
         l = 6
      elseif (ch(3:3).eq.'p') then
         l = 13
      elseif (ch(3:3).eq.'s') then
         l = 16
      else
         l = ichar(ch(3:3)) - ichar('F') + 3         
      endif
      return
      end
      
      function getchchar(n,l)
      character ch*1, getchchar*2
      ch(n) = char(mod(n,75) + ichar('0'))
c$$$      ch(n) = char(n + ichar('0'))

      if (l.eq.0) then
         getchchar = ch(n)//'S'
      elseif (l.eq.1) then
         getchchar = ch(n)//'P'
      elseif (l.eq.2) then
         getchchar = ch(n)//'D'
      elseif (l.eq.3) then
         getchchar = ch(n)//'F'
      elseif (l.eq.4) then
         getchchar = ch(n)//'G'
      elseif (l.eq.13) then
         getchchar = ch(n)//'p'
      elseif (l.eq.16) then
         getchchar = ch(n)//'s'
      elseif (l.ge.5) then
         getchchar = ch(n)//char(l-4+ichar('G'))
      endif
      return
      end

      function positron(n,l,npos)
      include 'par.f'
      logical positron
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      positron = .false.
      if (l.le.lptop.and.n.ge.natop(l)-nptop(l)+npbot(l)) then
         positron = .true.
         npos = n - natop(l) + nptop(l)
      endif 
      return
      end
      
      subroutine getchinfo (nch,nchp,lg,psi,maxpsi,enpsi,la,na,lp)
      include 'par.f'
      dimension psi(maxr)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      character chan(knm)*3, getchchar*2, chspin(3)*1, chspinp(3)*1
      logical positron, first_time
      common /charchan/ chan
      common /chanen/ enchan(knm),enchandiff(knm)
      common /channels/ nfich(nchan)
      integer satom
      common /helium/ latom(KNM), satom(KNM), lpar(KNM), np(KNM)
C  Do not change the lower case t below as it is used for checking
C  the spin of the incident state in the SOLVET routine
      data chspin/'s',' ','t'/
C  Use the capital letters for the rare case of unnatural parity
      data chspinp/'S',' ','T'/
      data nchpmax/0/
      data first_time/.true./
c-----------------------------------------------------------
      common /statenumber_str_to_scat/ Nstate_str_to_scat, na_ng(KNM)
      common /noblgegas_switch/ i_sw_ng    ! set in  cc/main.f
      type state
         real radial(maxr), energy, energydiff
         integer na, la, maxpsi, lapar
         character chn*3
      end type state
      type(state) states(knm)
      save first_time, states, nchpmax

      if(i_sw_ng .eq. 1) then

         iparity = (-1)**(lg + ipar)
         ncht = 0
         nchp = 1
      
         do N=1,Nstate_str_to_scat
            la = latom(N)
            lapar = lpar(N)
            na = na_ng(N)
            
c     This is a horrible "hack", its purpose is to include "unnatural" parity states, 
c     like positive parity P-state, into channel looping for lg=0 symmetry.
c     It produces additional channels for lg=0 that in principel should not be there 
c     but they are coming out with zero V-matrix due to the parity and ang.mom. considerations anyway... 
            lapar = (-1)**la
            if (lpar(nchp).ne.0.and.lg.ne.0) lapar = lpar(N)
c
            do lp = abs(lg-la), lg+la
               if (((-1)**lp)*lapar.ne.iparity) cycle
               ncht = ncht + 1
               if (nch.eq.ncht) then
                  nchp = N
                  nfich(ncht) = nchp
                  enpsi = enpsinb(na,la)
                  enchan(nchp) = enpsi
                  if (na.gt.la+1) then
                     enchandiff(nchp)=enpsinb(na+1,la)-enpsinb(na-1,la)
                  else
                     enchandiff(nchp) = 1e0
                  endif
!            chan(nchp) = chspin(2*satom(nchp)+1)//getchchar(na,la)
!            if (lpar(nchp).eq.-(-1)**latom(nchp)) chan(nchp) 
!     >         = chspinp(2*satom(nchp)+1)//getchchar(na,la)
!
!                  print*
!                  print'("!chinfo:",10I5)',nch,ncht,N,la,lapar,ipar
                  return
               elseif(nch .lt. ncht) then
                  nch = 0
c                  print*, 'Could not find this channel: nch=',nch,
c     >                 ', lg,ipar:',lg,ipar
                  return
               endif
               
            enddo
         enddo
c         print*, 'Could not find this channel: nch=',nch, 
c     >        ', lg,ipar:',lg,ipar
         nch = 0
         return
      endif
c-----------------------------------------------------------

      nchp = 0
      iparity = (-1)**(lg + ipar)
      ncht = 0
c$$$      nchp = nchp + 1
c$$$      la = linc
c$$$      na = ninc
c$$$      lapar = (-1)**la
c$$$      if (lpar(nchp).ne.0.and.lg.ne.0) lapar = lpar(nchp)
c$$$      do i = 1, maxr
c$$$         psi(i) = 0.0
c$$$      enddo 
c$$$      if (ipar.eq.0) then
c$$$         enchan(nchp) = enpsinb(na,la)
c$$$         if (np(1).eq.0) then
c$$$            chan(nchp) = chspin(2)//getchchar(na,la)
c$$$            latom(nchp) = la
c$$$         else
c$$$            chan(nchp) = chspin(2*satom(nchp)+1)//getchchar(np(nchp),la)
c$$$            if (lpar(nchp).eq.-(-1)**latom(nchp)) chan(nchp) 
c$$$     >         = chspinp(2*satom(nchp)+1)//getchchar(np(nchp),la)
c$$$         endif
c$$$      endif 
c$$$      do 10 lp = abs(lg-la), lg+la
c$$$         if ((-1)**lp*lapar.ne.iparity) go to 10
c$$$         ncht = ncht + 1
c$$$         nfich(ncht) = nchp
c$$$         if (nch.eq.ncht) then
c$$$            enpsi = enpsinb(na,la)
c$$$            maxpsi = maxpsinb(na,la) 
c$$$            do i = 1, maxpsi
c$$$               psi(i) = psinb(i,na,la)
c$$$            end do
c$$$            return
c$$$         end if 
c$$$ 10   continue

      if (first_time.or.np(1).ne.0) then
cCCC      do la = labot, latop
cCCC         do na = nabot(la), natop(la)
C  The above looping is preferable, but the CROSS program expects it in
C  the order below. 
      nastart = nabot(labot)
      nastop  = natop(labot)
      do la = labot, latop
         if (nabot(la).le.la) then
            print*,'NABOT(LA) must be > LA, here NABOT(LA), LA:',
     >         nabot(la),la
         endif
         if (nabot(la).lt.nastart) nastart = nabot(la)
         if (natop(la).gt.nastop)  nastop  = natop(la)
      enddo
      
      do natmp = nastart, nastop
         do la = labot, min(natmp - 1, latop)
C  The following three lines are there for potassium to stop 3d being 
C  the first channel
            na = natmp
            if (la.le.1.and.nastart.lt.nabot(la).and.la+1.lt.nabot(la))
     >         na = natmp + nabot(la) - nastart
! The following statement ensures Ps states start after all He-like states
            if (np(1).ne.0.and.positron(natop(0),0,npos)
     >         .and.enpsinb(na,la).eq.0.0) cycle
            if (na.ge.nabot(la).and.na.le.natop(la)
c$$$     >         .and.(la.ne.linc.or.na.ne.ninc)
     >         ) then
               nchp = nchp + 1
               lapar = (-1)**la
               if (lpar(nchp).ne.0.and.lg.ne.0) lapar = lpar(nchp)
               if (ipar.eq.0) then
                  if (np(1).eq.0.or.positron(na,la,npos)) then
                     chan(nchp) = chspin(2)//getchchar(na,la)
                     latom(nchp) = la
                     if (positron(na,la,npos)) chan(nchp)='p'//
     >                  getchchar(npos,la)
                  else
                     chan(nchp) = chspin(2*satom(nchp)+1)//
     >                  getchchar(np(nchp),la)
                     if (lpar(nchp).eq.-(-1)**latom(nchp)) chan(nchp) =
     >                 chspinp(2*satom(nchp)+1)//getchchar(np(nchp),la)
                  endif
                  enchan(nchp) = enpsinb(na,la)
                  if (na.gt.la+1) then
                     enchandiff(nchp)=enpsinb(na+1,la)-enpsinb(na-1,la)
                  else
                     enchandiff(nchp) = 1e0
                  endif
                  states(nchp)%la = la
                  states(nchp)%na = na
                  states(nchp)%lapar = lapar
                  states(nchp)%chn = chan(nchp)
                  states(nchp)%energy = enchan(nchp)
                  states(nchp)%energydiff = enchandiff(nchp)
                  states(nchp)%radial(1:maxpsinb(na,la))=
     >                 psinb(1:maxpsinb(na,la),na,la)
                  states(nchp)%maxpsi=maxpsinb(na,la)
c$$$                  print*,'ENERGY:',enchan(nchp)
               endif 
               do 20 lp = abs(lg-la), lg+la
                  if ((-1)**lp*lapar.ne.iparity) go to 20
                  ncht = ncht + 1
                  nfich(ncht) = nchp
                  if (nch.eq.ncht) then
                     enpsi = enpsinb(na,la)
                     maxpsi = maxpsinb(na,la) 
                     do i = 1, maxpsi
                        psi(i) = psinb(i,na,la)
                     end do
                     return
                  end if 
 20            continue 
            end if 
         end do
      end do
      nchpmax = nchp
      if (np(1).eq.0) then
c$$$         print*,' reordering in getchinfo with NCHPMAX:',nchpmax
         call reorderstates(states,nchpmax,lptop.ge.0)
      endif 
      else ! not first_time
         if (nch.eq.0) stop 'nch = 0 in getchinfo'
         do nchp = 1, nchpmax
            la = states(nchp)%la
            lapar = states(nchp)%lapar
            latom(nchp) = la
            do 30 lp = abs(lg-la), lg+la
               if ((-1)**lp*lapar.ne.iparity) go to 30
               ncht = ncht + 1
               nfich(ncht) = nchp
               if (nch.eq.ncht) then
                  enpsi = states(nchp)%energy
                  enchan(nchp) = enpsi
                  enchandiff(nchp) = states(nchp)%energydiff
                  na = states(nchp)%na
                  chan(nchp) = states(nchp)%chn
                  maxpsi = states(nchp)%maxpsi
                  psi(1:maxpsi) = states(nchp)%radial(1:maxpsi)
                  psi(maxpsi+1:maxr) = 0.0
                  return
               end if 
 30         continue 
         enddo 
      endif 
      nch = 0
      first_time = .false.
      return
      end

      subroutine reorderstates(states,nchpmax,ps)
      include 'par.f'
      logical swapped,ps
      type state
         real radial(maxr), energy, energydiff
         integer na, la, maxpsi, lapar
         character chn*3
      end type state
      type(state) states(knm), staten

      if (ps) then
         nstart = 3
      else
         nstart = 2
      endif
      swapped = .true.
      do while (swapped)
         swapped = .false.
         do n = nstart, nchpmax !Leave the 1st alone to stop Ps(1s) becoming 1st
c$$$            if (states(n-1)%energy.gt.states(n)%energy) then !*0.99999) then
            e1 = states(n-1)%energy
            e2 = states(n)%energy
            if (e1.gt.e2) then
               swapped = .true.
c$$$           print*,'swapping:',abs((e1-e2)/(e1+e2)),n-1,n,states(n-1)%la,
c$$$     >              states(n-1)%energy,states(n)%la,states(n)%energy
               staten = states(n-1)
               states(n-1) = states(n)
               states(n) = staten
            endif
         enddo
      enddo 
      swapped = .true.
      do while (swapped)
         swapped = .false.
         do n = nstart, nchpmax !Leave the 1st alone to stop Ps(1s) becoming 1st
c$$$            if (states(n-1)%energy.gt.states(n)%energy) then !*0.99999) then
            e1 = states(n-1)%energy
            e2 = states(n)%energy
            if (abs((e1-e2)/(e1+e2)).lt.1e-5.and.
     >         states(n-1)%la.gt.states(n)%la) then
               swapped = .true.
c$$$           print*,'swapping:',abs((e1-e2)/(e1+e2)),n-1,n,states(n-1)%la,
c$$$     >              states(n-1)%energy,states(n)%la,states(n)%energy
               staten = states(n-1)
               states(n-1) = states(n)
               states(n) = staten
            endif
         enddo
      enddo 
      return
      end
 
      subroutine makev1e(nqmi,psii,maxpsii,ei,lia,li,
     >   chii,minchii,gki,ni,etot,thetain,ne2e,
     >   nqmf,psif,maxpsif,ef,lfa,lf,chif,minchif,gkf,nf,
!     >   lg,const,nqmfmax,uf,ui,nchf,nchi,nold,nznuci,npk,ve2e,vmatt,
     >   lg,const,uf,ui,nchf,nchi,nold,nznuci,npk,ve2e,nqmfmax,vmatt,
     >   nsmax,nchtop)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension chii(meshr,nqmi),minchii(nqmi),gki(kmax)
      dimension chif(meshr,nqmf),minchif(nqmf),gkf(kmax)
      dimension psii(maxr), psif(maxr), uf(maxr,nchan), ui(maxr), 
     >          npk(nchan+1)
!      real vmatt(kmax,kmax,0:1,nchtop),ve2e(nf,ni), temp(maxr),
      real vmatt(nqmfmax,nqmfmax,nchi:nchtop,0:nsmax),ve2e(nf,ni), 
     >   temp(maxr), ovlpf(kmax),ovlpkf(kmax), theta(0:lamax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
  
      do la = labot, latop
         theta(la) = thetain
      enddo
c$$$      if (int(thetain).eq.-1) then
c$$$         do la = labot, latop
c$$$            theta(la) = float(la+1) / 3.0
c$$$         enddo
c$$$      endif 
      if(li.eq.lf.and.lia.eq.lfa.and.li.le.latop.and.thetain.ne.0.0)then
         ovlp = 0.0
         do i = 1, min(maxpsif,maxpsii)
           ovlp = ovlp + psii(i) * psif(i) * rmesh(i,3)
         enddo
         tmp = - ovlp * const * etot * theta(li) / 2.0
         do n = nabot(li), natop(li)
            ep = enpsinb(n,li)
C  The following form isn't coded correctly, but can be made to work
cCCC  tmp = ovlp * const * (ei + ep - etot) / 2.0 * theta(lfa)
            do kf = 1, nqmf
               ovlpf(kf) = 0.0
               do i = 1, maxpsinb(n,li)
                  ovlpf(kf) = ovlpf(kf) + psinb(i,n,li) * chif(i,kf)
               enddo
            enddo 
            do ki = 1, nqmi
               kii = npk(nchi) + ki - 1
               ovlpi = 0.0
               do i = 1, maxpsinb(n,li)
                  ovlpi = ovlpi + psinb(i,n,li) * chii(i,ki)
               enddo
               do 5 kf = 1, nqmf

cCCC                  tmp = ovlp * const * (abs(gkf(kf)) * gkf(kf) + ep
cCCC     >               - etot) / 4.0
                  
c$$$c$$$                 if (ne2e.eq.0) then
c$$$                     kff = npk(nchf) + kf - 1
c$$$                     if (kff.lt.kii) go to 5
c$$$                     vmat(kii,kff+1) = vmat(kii,kff+1) +
c$$$     >                  tmp * ovlpf(kf) * ovlpi
c$$$c$$$                  else
c$$$c$$$                     kff = nchf
c$$$c$$$                  endif 
c$$$                  vmat(kff,kii) = vmat(kff,kii) +
c$$$     >               tmp * ovlpf(kf) * ovlpi
                  vmatt(kf,ki,nchf,0) = vmatt(kf,ki,nchf,0) 
     >               +tmp*ovlpf(kf)*ovlpi
                  if (nsmax.eq.1)
     >               vmatt(kf,ki,nchf,1) = vmatt(kf,ki,nchf,1) 
     >               +tmp*ovlpf(kf)*ovlpi
 5             continue 
            enddo
         enddo 
      endif
      if (li.eq.lfa.and.lf.eq.lia) then
         sign = (-1) ** (li + lia - lg)         
         do kf = 1, nqmf
            ekf = abs(gkf(kf)) * gkf(kf)
            ovlpf(kf) = 0.0
            ovlpkf(kf) = 0.0
            if (nznuc.gt.1) then
               do i = 1, meshr
                  temp(i) = 0.0
               enddo
               do i = 1, maxpsii
                  temp(i) = psii(i) * rmesh(i,3)
               enddo 
               call fcexch(temp,meshr,chif(1,kf),ovlpkf(kf),lf)
               ovlpkf(kf) = - ovlpkf(kf)
            endif
            do i = minchif(kf), maxpsii
               ovlpf(kf) = ovlpf(kf) + psii(i) * chif(i,kf)
               ovlpkf(kf) = ovlpkf(kf) + psii(i) * chif(i,kf) *
     >            (uf(i,nchf) / 2.0 + rpow2(i,0))
cCCC     >            (ei / 2.0 + uf(i) / 2.0 + 1.0 / rmesh(i,1))
            enddo
         enddo 
         do ki = 1, nqmi
            kii = npk(nchi) + ki - 1
            ovlpi = 0.0
            ovlpki = 0.0
            eki = abs(gki(ki)) * gki(ki)
            if (nznuc.gt.1) then
               do i = 1, meshr
                  temp(i) = 0.0
               enddo
               do i = 1, maxpsif
                  temp(i) = psif(i) * rmesh(i,3)
               enddo 
               call fcexch(temp,meshr,chii(1,ki),ovlpki,li)
               ovlpki = - ovlpki
            endif
C  For pure plane/Coulomb waves UI and UF both are  = - VDCORE * 2.0,
C  otherwise they are just the short ranged one-electron distorting potentials
C  The asymptotic part of V_{FC} (ZEFF/r) always cancels.
            do i = minchii(ki), maxpsif
               ovlpi = ovlpi + psif(i) * chii(i,ki)
               ovlpki = ovlpki + psif(i) * chii(i,ki) *
     >            (ui(i) / 2.0 + rpow2(i,0)) 
cCCC     >            (ef / 2.0 + ui(i) / 2.0 + 1.0 / rmesh(i,1)) 
            enddo
cCCC            print*,'ovlpki,ovlpi*e:',ki,ovlpki*2.0-ovlpi*eki,ovlpi*ef,
cCCC     >         (ovlpki*2.0-ovlpi*eki)/(ovlpi*ef+1e-30)
             do 10 kf = 1, nqmf
cCCC               if (ni.eq.nf) then
                kff = npk(nchf) + kf - 1
c$$$               if (ne2e.eq.0) then
c$$$                  kff = npk(nchf) + kf - 1
c$$$                  if (kff.lt.kii) go to 10
c$$$               else
c$$$                  kff = nchf
c$$$               endif 
               ekf = abs(gkf(kf)) * gkf(kf)

               etotnew = etot
C  The following line is Andris's first guess to implement step function SDCS
c$$$               if (nchi.eq.nchf.and.ei.gt.0.0.and.ei.lt.etot)
c$$$     >            etotnew = 0.0
C  The following line is Andris's second guess to implement step function SDCS
c$$$               if (nchi.eq.nchf.and.ei.gt.etot/2.0.and.ei.lt.etot)
c$$$     >            etotnew = 0.0
               eterm = 0.5 * (eki + ekf - etotnew * (1.0 - theta(lia)))

               
               tmp = eterm * ovlpi * ovlpf(kf) - ovlpki * ovlpf(kf) -
     >            ovlpkf(kf) * ovlpi

               if (nold.eq.0)
     >            tmp = (ei + ef - etot * (1.0 - theta(lfa))) / 2.0 *
     >            ovlpi * ovlpf(kf)    
cCCC     >            tmp = 0.5 * (ei + ef - etot) * ovlpi *
cCCC     >            ovlpf(kf) * (1.0 - theta(lfa))
               
c$$$               if (ne2e.ne.0) then
c$$$cCCC                  ve2e(kf,nchf,ki,nchi) = ve2e(kf,nchf,ki,nchi)
c$$$                  ve2e(kff,kii) = ve2e(kff,kii)
c$$$     >               + tmp * const * sign
c$$$               else

c$$$                  vmat(kii,kff+1) = vmat(kii,kff+1)
c$$$     >                - tmp * const * sign
c$$$                  vmat(kff,kii) = vmat(kff,kii) 
c$$$     >                + tmp * const * sign
               vmatt(kf,ki,nchf,0)=vmatt(kf,ki,nchf,0)+tmp*const*sign
               if (nsmax.eq.1)
     >            vmatt(kf,ki,nchf,1)=vmatt(kf,ki,nchf,1)-tmp*const*sign

c$$$  endif 
 10         continue 
         enddo 
      endif 
      end
               
!not used any more
c$$$      subroutine makev31d(dwpot,nqmi,psii,maxpsii,lia,li,nia,chii,
c$$$     >   minchii,gki,phasei,ni,nqmf,psif,maxpsif,lfa,lf,nfa,chif,
c$$$     >   minchif,gkf,phasef,nf,lg,rnorm,u,itail,nznuci,vdon,nchf,nchi,
c$$$     >   pos,second,npk,ne2e,vmatt)
c$$$      use ubb_module
c$$$      use apar
c$$$      use chil_module
c$$$      include 'par.f'
c$$$      common/meshrr/ meshr,rmesh(maxr,3)
c$$$      dimension chii(meshr,nqmi),minchii(nqmi),gki(kmax),
c$$$     >   phasei(kmax)
c$$$      dimension chif(meshr,nqmf),minchif(nqmf),gkf(kmax),
c$$$     >   phasef(kmax),npk(nchan+1)
c$$$      complex phasei, phasef
c$$$      dimension fun(maxr), temp(maxr), psii(maxr), psif(maxr),
c$$$     >   ovlpi(kmax), ovlpf(kmax), temp2(maxr), temp3(maxr),
c$$$     >   tail(kmax,kmax),u(maxr), dwpot(maxr,nchan)
c$$$      real vmatt(kmax,kmax,0:1),vdon(nchan,nchan,0:1)
c$$$      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
c$$$     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
c$$$      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
c$$$     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
c$$$      common /psinbc/ enpsinb(nnmax,0:lnabmax),
c$$$     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
c$$$      common/numericalpotential/ numericalv, lstoppos
c$$$      logical numericalv
c$$$      common/matchph/rphase(kmax,nchan)
c$$$*      real chitemp(maxr,kmax) !, rpowpos1(maxr), rpowpos2(maxr)
c$$$      allocatable :: chitemp(:,:)
c$$$      logical pos, second
c$$$      common/smallr/ formcut,regcut,expcut,fast,match,analyticd
c$$$      logical fast, match, analyticd
c$$$      real*8 xin(maxr),yin(maxr),xout(maxr),yout(maxr)
c$$$      character chan(knm)*3
c$$$      common /charchan/ chan
c$$$
c$$$      real prod(kmax,kmax)
c$$$c
c$$$      common /di_el_core_polarization/ gamma, r0, pol(maxr)   
c$$$      real*8, dimension (meshr) :: ubbb ! ANDREY
c$$$      character ch*1
c$$$      ch(n) = char(mod(n,75) + ichar('0'))
c$$$      
c$$$      hat(l)= sqrt(2.0 * l + 1.0)
c$$$
c$$$      minfun = 1
c$$$      maxfun = min(maxpsii,maxpsif)
c$$$      do i = minfun, maxfun
c$$$         fun(i) = psii(i) * psif(i) * rmesh(i,3)
c$$$      end do 
c$$$
c$$$      do i = 1, maxr
c$$$         temp3(i) = 0.0
c$$$      enddo
c$$$      mini = maxr
c$$$      maxi = 1
c$$$      ctemp = 0.0
c$$$
c$$$
c$$$      ltmin = 1000000
c$$$      ctemp = 0.0
c$$$      do 10 ilt = -lia, lia, 2
c$$$         lt = lfa + ilt
c$$$         if (lt.lt.0.or.lt.gt.ltmax) go to 10
c$$$         call cleb(2*lia,2*lt,2*lfa,0,0,0,c1)
c$$$         c1tmp = cgc0(float(lia),float(lt),float(lfa))
c$$$         if (abs((c1-c1tmp)/(c1tmp+1e-20)).gt.1e-3) then
c$$$            print*,'CGCs 1 do not agree:',c1, c1tmp,lia,lt,lfa
c$$$            stop 'CGCs 1 do not agree'
c$$$         endif 
c$$$         call cleb(2*li,2*lf,2*lt,0,0,0,c2)
c$$$         c2tmp = cgc0(float(li),float(lf),float(lt))         
c$$$         if (abs((c2-c2tmp)/(c2tmp+1e-20)).gt.1e-3) then
c$$$            c2tmp2 = cgcigor(2*li,2*lf,2*lt,0,0,0)
c$$$C$OMP critical(print)
c$$$            print*,'CGCs 2 do not agree:',c2, c2tmp,c2tmp2,li,lf,lt
c$$$C$OMP end critical(print)
c$$$            c2 = c2tmp2
c$$$         endif 
c$$$         call rac7(2*li,2*lf,2*lia,2*lfa,2*lt,2*lg,c3)
c$$$         c3tmp = cof6j(float(li),float(lf),float(lt),float(lfa),
c$$$     >      float(lia),float(lg))*(-1)**(li+lf+lia+lfa)
c$$$         if (abs((c3-c3tmp)/(c3tmp+1e-6)).gt.1e-2) then
c$$$C  below does happen for NPAR = 1, but not significant
c$$$C$OMP critical(print)
c$$$            print*,'WARNING: CJ6 and W do not agree in D:',c3, c3tmp,
c$$$     >         li,lf,lia,lfa,lt,lg
c$$$C$OMP end critical(print)
c$$$            c3 = c3tmp
c$$$c$$$            stop 'CJ6 and W do not agree in D'
c$$$         endif 
c$$$         
c$$$         c = c1 * c2 * c3
c$$$         if (abs(c).lt.1e-10) go to 10
c$$$         const = (-1)**(lg + lfa) * hat(li) *
c$$$     >      hat(lf) * hat(lia) / hat(lt) * c
c$$$         call form(fun,minfun,maxfun,rpow1(1,lt),rpow2(1,lt),
c$$$     >      minrp(lt),maxrp(lt),meshr,temp,i1,i2)
c$$$c$$$         if (lt.eq.0) print*,'lia,li,lg,lfa,lf,c:',lia,li,lg,lfa,lf,c
c$$$c
c$$$       if(lt.eq.1.and.gamma.ne.0.0) then
c$$$          sum1 = 0.0
c$$$          do i=minfun,maxfun
c$$$             sum1 = sum1 + fun(i)*pol(i)
c$$$          end do                       
c$$$          tmp = gamma * sum1 
c$$$          do i=i1,i2             
c$$$             temp(i) = temp(i) - tmp*pol(i)
c$$$          end do
c$$$       end if
c$$$       
c$$$c         
c$$$C  Subtract 1/r, but only for same atom-atom channels when lambda = 0
c$$$         if (lt.eq.0.and..not.pos.and.nchi.eq.nchf) then
c$$$cCCC  if (lt.eq.0) then
c$$$
c$$$c$$$            if (alkali) then
c$$$c$$$C ANDREY: subtract (ve(r)+1)/r for pos-alkali case             
c$$$c$$$               do i=1,i2
c$$$c$$$                  ! factor 2 is since u(i)/2 appears in nuclear(..)
c$$$c$$$                  u(i)=2.0d0*veint(rmesh(i,1))*rpow2(i ,0)
c$$$c$$$               end do
c$$$c$$$            end if                                 
c$$$            call nuclear(fun,.not.pos,minfun,maxfun,i2,u,nznuc,temp)
c$$$c$$$C     ANDREY: subtract (ve(r)+1)/r for pos-alkali atom scattering
c$$$c$$$            if (alkali) then           
c$$$c$$$               tmp=0.0
c$$$c$$$               do i=minfun,maxfun
c$$$c$$$                  tmp=tmp+fun(i)
c$$$c$$$               end do
c$$$c$$$               irstop = i2 
c$$$c$$$               do i=1,irstop
c$$$c$$$                  temp(i)=temp(i)-tmp*veint(rmesh(i,1))*rpow2(i ,0)
c$$$c$$$               end do
c$$$c$$$               do while (abs(temp(irstop)).lt.formcut.and.irstop.gt.1)
c$$$c$$$                  temp(irstop) = 0.0
c$$$c$$$                  irstop = irstop - 1
c$$$c$$$               end do
c$$$c$$$               i2=irstop
c$$$c$$$            end if             
c$$$         endif
c$$$
c$$$         if (pos) then !.and.lg.le.lstoppos) then  !.and.nqmi.gt.1
c$$$C     ANDREY: Hydrogen: ssalling + interpolation:          
c$$$            if (.not.alkali) then               
c$$$C  The factor of two below is the reduced mass. We are working with matrix
c$$$C  elements whose channels are multiplied by sqrt of the reduced mass. Here
c$$$C  we have positronium-positronium matrix element, hence a factor of 2 overall.
c$$$c$$$            const = const * (1-(-1)**lt)
c$$$               const = 2.0 * const * (1-(-1)**lt)
c$$$c$$$               const = - const !Rav's derivation, but not Alisher's
c$$$               do i = i1, i2
c$$$                  xin(i-i1+1) = rmesh(i,1)
c$$$                  yin(i-i1+1) = temp(i) * xin(i-i1+1) ** (lt+1)
c$$$                  xout(i-i1+1) = rmesh(i,1) * 2d0
c$$$               enddo 
c$$$               if (i2-i1+1.gt.maxr) then
c$$$                  print*,'I2,I1,LT:',i2,i1,lt
c$$$                  stop 'problem in call to intrpl'
c$$$               endif 
c$$$               call intrpl(i2-i1+1,xin,yin,i2-i1+1,xout,yout)
c$$$c$$$            nnn = nnn + 1
c$$$               do i = i1, i2
c$$$                  if (xout(i-i1+1).gt.xin(i2-i1+1))
c$$$     >                 yout(i-i1+1) = yin(i2-i1+1)
c$$$                  temp(i) = 2.0 * (yout(i-i1+1) / xout(i-i1+1)**(lt+1))
c$$$c$$$               write(50+nnn,'(1p,4e10.3,i5)') xout(i-i1+1),
c$$$c$$$     >            yout(i-i1+1) / xout(i-i1+1)**(lt+1), 
c$$$c$$$     >            xin(i-i1+1), yin(i-i1+1) / xin(i-i1+1)**(lt+1), i
c$$$               enddo
c$$$            else ! alkali              
c$$$C     ANDREY: pos-alkali case: radial integrals with Ubb               
c$$$               if (lt.gt.ubb_max3) stop
c$$$     $              'stopped: vmat: Ubb with larger lt are needed'
c$$$               const=2d0*const               
c$$$C              some constants from get_ubb: 
c$$$               rlgmin = log10(rmesh(1,1))
c$$$               rlgmax = log10(rmesh(meshr,1))
c$$$               drlg = (rlgmax-rlgmin)/real(ubb_max1-1)
c$$$
c$$$c---- find irmin and irmax needed for interpolation                  
c$$$               irmin = 1;
c$$$               rmin = rmesh(minfun, 1);  
c$$$               ir = 2;               
c$$$               do while (arho(ir).lt.rmin)                  
c$$$                  ir = ir+1
c$$$               end do
c$$$               irmin = ir-1               
c$$$C     ---------------------------------------------- 
c$$$               rmax = rmesh(maxfun, 1)  
c$$$               irmax = ubb_max1; ir = irmax;
c$$$               do while (arho(ir).gt.rmax)                  
c$$$                  ir = ir-1
c$$$               end do
c$$$               irmax = ir+1               
c$$$               maxir = irmax-irmin+1
c$$$C---------------------------------------------------
c$$$               
c$$$               do i = 1, i2
c$$$                  rho = rmesh(i,1);                                    
c$$$                  call get_ubb(minfun,maxfun,irmin,irmax,maxir,drlg,lt,
c$$$     $                 rho,ubbb)      
c$$$                  temp(i) = 0.0                 
c$$$                  do j = minfun, maxfun
c$$$                     r = rmesh(j,1)
c$$$c     explicit integration of Ubb (Eq. 23): 
c$$$C                    call Ubb0(lt, r, rho, res0, res)
c$$$C     Ubb integrals computed outside:
c$$$C                    res = ubb_res(j,i,lt)
c$$$C     interpolation: minus: alisher uses (1-(-1)**lt):
c$$$                     res = -ubbb(j)/max(rho,r/2.)
c$$$                     temp(i) = temp(i)+res*fun(j)
c$$$                  end do        ! j
c$$$               end do           ! i               
c$$$            end if              ! alkali or hydrogen                          
c$$$         else                   ! pos
c$$$C  Multiply by -1 for positron scattering as then NZE = 1, for non pos-pos
c$$$C  channels 
c$$$            const = - nze * const
c$$$         endif ! pos
c$$$         
c$$$         do i = i1, i2
c$$$            temp3(i) = const * temp(i) + temp3(i)
c$$$         enddo
c$$$         mini = min(mini,i1)
c$$$         maxi = max(maxi,i2)
c$$$C  The variable CTEMP is the asymptotic value of const * temp * r ** (lt+1)
c$$$c$$$         if ((lt.eq.1.or.lt.eq.2).and.i2.eq.meshr) then
c$$$            if (i2.eq.meshr.and.lt.ne.0.and.lt.lt.ltmin) then
c$$$               ctemp = rnorm * const * temp(i2)/rpow2(i2,lt)
c$$$               ltmin = lt
c$$$c$$$               print*,'ctemp,lt:',ctemp,lt
c$$$c$$$            ctemp = - float(nze) * rnorm * const * temp(i2)/rpow2(i2,lt)
c$$$         endif 
c$$$ 10   continue
c$$$c$$$      if (pos.and.maxi.eq.meshr) then
c$$$c$$$         print*,'Writing temp3 to disk'
c$$$c$$$         open(42,file='temp'//ch(lfa)//ch(lia)//ch(lt))
c$$$c$$$         write(42,*) '# lt,ltmin,ctemp:',lt,ltmin,ctemp
c$$$c$$$         do i = 1, meshr
c$$$c$$$            write(42,*) rmesh(i,1), temp3(i)
c$$$c$$$         enddo
c$$$c$$$         close(42)
c$$$c$$$      endif 
c$$$            
c$$$c$$$      if (ltmin.lt.10.and.maxi.eq.meshr)
c$$$c$$$     >   print*,'pos,lfa,lia,ltmin,maxi,ctemp,v(i):',pos,lfa,lia,ltmin,
c$$$c$$$     >   maxi,ctemp,(rnorm*temp3(i)/rpow2(i,ltmin),i=maxi-20,maxi,10)
c$$$      if (ltmin.lt.10.and.maxi.eq.meshr) then
c$$$         ctemp = rnorm * temp3(maxi) / rpow2(maxi,ltmin)
c$$$c$$$      else
c$$$c$$$         ctemp = 0.0
c$$$      endif 
c$$$      if (nchi.eq.nchf.and.nqmi.eq.1.and.u(1).eq.0.0) then
c$$$C  Define the channel dependent distorting potential when running the
c$$$C  Born case. The units are Rydbergs.
c$$$         do i = 1, meshr
c$$$            dwpot(i,nchi) = 0.0
c$$$         enddo
c$$$         do i = mini, maxi
c$$$            dwpot(i,nchi) = temp3(i) * 2.0
c$$$         enddo
c$$$      endif
c$$$      
c$$$C  As both CHII and CHIF contain the integration weights, we divide TEMP by   
c$$$C  them.
c$$$c$$$      open(42,file=chan(nchf)(2:3)//'.'//chan(nchi)(2:3))
c$$$c$$$      write(42,'("# ", i2, a4, "<-", a4, i3)')
c$$$c$$$     >   lf,chan(nchf),chan(nchi),li
c$$$      do i = mini, maxi
c$$$c$$$         write(42,*) i, rmesh(i,1), temp3(i)
c$$$         temp3(i) = rnorm * temp3(i) / rmesh(i,3)
c$$$      end do
c$$$c$$$      close(42)
c$$$      mini1 = mini
c$$$C put IF statement as chil(,,2) may not be allocated, but see below
c$$$      tail(:,:) = 0.0
c$$$      if (itail.eq.-1) then
c$$$      call maketail(itail,ctemp,chil(1,npk(nchi),2),minchil(npk(nchi),2)
c$$$     >   ,gki,phasei,li,nqmi,chil(1,npk(nchf),2),minchil(npk(nchf),2)
c$$$     >   ,gkf,phasef,lf,nqmf,nchf,nchi,ltmin,kmax,tail)
c$$$      elseif (itail.eq.1) then
c$$$      call maketail(itail,ctemp,chil(1,npk(nchi),1),minchil(npk(nchi),1)
c$$$     >   ,gki,phasei,li,nqmi,chil(1,npk(nchf),1),minchil(npk(nchf),1)
c$$$     >   ,gkf,phasef,lf,nqmf,nchf,nchi,ltmin,kmax,tail)
c$$$      endif
c$$$c$$$      call maketail(itail,ctemp,chii,minchii,gki,phasei,li,nqmi,
c$$$c$$$     >   chif,minchif,gkf,phasef,lf,nqmf,nchf,nchi,ltmin,tail)
c$$$
c$$$      allocate(chitemp(maxr,kmax))
c$$$      do ki = 1, nqmi
c$$$         minki = max(mini1, minchii(ki))
c$$$         do i = minki, maxi
c$$$            chitemp(i,ki) = temp3(i) * chii(i,ki)
c$$$         enddo
c$$$      enddo 
c$$$
c$$$C  The following gets correct prod(kmax,kmax) but the transfer takes all
c$$$C  the time! Need to transfer all of CHIL at the one time, but....
c$$$c$$$  if (nqmi.gt.1.and.nqmf.gt.1) 
c$$$c$$$     >     call cmmult(temp3,chif,chitemp,meshr,nqmf,nqmi,mini,maxi,
c$$$c$$$     >     prod)
c$$$
c$$$c$$$c$par doall      
c$$$      if (nchi.eq.nchf) then
c$$$      do ki = 1, nqmi
c$$$         qi = gki(ki)
c$$$         kii = npk(nchi) + ki - 1
c$$$         minki = max(mini1, minchii(ki))
c$$$         kfstop = nqmf
c$$$         if (second) then
c$$$            if (ki.eq.1) then
c$$$               kfstop = nqmf
c$$$            else
c$$$               kfstop = 1
c$$$            endif
c$$$         endif 
c$$$         do 20 kf = 1, kfstop
c$$$            qf = gkf(kf)
c$$$
c$$$            if (ne2e.eq.0) then
c$$$               kff = npk(nchf) + kf - 1
c$$$               if (kff.lt.kii) go to 20
c$$$            else
c$$$               kff = nchf
c$$$            endif 
c$$$            mini = max(minki, minchif(kf))
c$$$            n = maxi - mini + 1
c$$$            tmp = tail(kf,ki)
c$$$            tmp = tmp +
c$$$c     >         cuDDot(chif,chitemp,meshr,nqmf,maxr,kmax,mini,maxi,kf,ki)
c$$$     >         dot_product(chif(mini:maxi,kf),chitemp(mini:maxi,ki))
c$$$            vmatt(kf,ki,0) = vmatt(kf,ki,0) + tmp 
c$$$            vmatt(kf,ki,1) = vmatt(kf,ki,1) + tmp 
c$$$ 20      continue 
c$$$      end do 
c$$$      else
c$$$
c$$$      do ki = 1, nqmi
c$$$c         qi = gki(ki)
c$$$c         kii = npk(nchi) + ki - 1
c$$$         minki = max(mini1, minchii(ki))
c$$$c         kfstop = nqmf
c$$$         do kf = 1, nqmf
c$$$c            qf = gkf(kf)
c$$$c            kff = npk(nchf) + kf - 1
c$$$            mini = max(minki, minchif(kf))
c$$$            tmp = 
c$$$c     >         cuDDot(chif,chitemp,meshr,nqmf,maxr,kmax,mini,maxi,kf,ki)
c$$$     >         dot_product(chif(1:maxi,kf),chitemp(1:maxi,ki))
c$$$c            tmp = 0.0
c$$$c            do i = 1, maxi
c$$$c               tmp = tmp + chif(i,kf)*chii(i,ki)*temp3(i)
c$$$c            enddo
c$$$            vmatt(kf,ki,0) = vmatt(kf,ki,0) + tmp
c$$$            vmatt(kf,ki,1) = vmatt(kf,ki,1) + tmp
c$$$      end do
c$$$      end do
c$$$c      call cuDDOT(chif,chii,temp3,nqmf,nqmi,meshr,kmax,maxi,vmatt)
c$$$      end if
c$$$      deallocate(chitemp)
c$$$
c$$$      if (itail.gt.0.and..not.pos.and.abs(lf-li).eq.1) then
c$$$         do ki = 1, nqmi
c$$$            qi = gki(ki)
c$$$            kf = nqmf-1
c$$$            do while (kf.gt.1)
c$$$               qf = gkf(kf)
c$$$               if (qf*rmesh(meshr,1).lt.itail.and.qf.lt.qi) then
c$$$                  estimate = vmatt(kf+1,ki,0)*qf**2/gkf(kf+1)**2
c$$$c$$$                  print*,kf,qf,estimate,gkf(kf+1),vmatt(kf+1,ki,0)
c$$$                  vmatt(kf,ki,0) = estimate
c$$$                  vmatt(kf,ki,1) = estimate
c$$$               endif
c$$$               kf = kf - 1
c$$$            enddo
c$$$         enddo
c$$$      endif 
c$$$                  
c$$$C  Define the direct on shell V matrix used for analytic Born subtraction
c$$$C  in the cross program
c$$$      ki = 1
c$$$      kf = 1
c$$$      tmp = tail(kf,ki)
c$$$      mini = max(mini1, minchii(ki), minchif(kf))
c$$$      minki = max(mini1, minchii(ki))
c$$$      do i = minki, maxi
c$$$         temp2(i) = temp3(i) * chii(i,ki)
c$$$      enddo 
c$$$c$$$      n = maxi - mini + 1
c$$$c$$$      if(n.gt.0)tmp = tmp + sdot(n,chif(mini,kf),1,temp2(mini),1)
c$$$C The following If statement avoids a bug on the IBMs
c$$$      if (mini.lt.maxi) then
c$$$         tmp = tmp +
c$$$     >      dot_product(chif(mini:maxi,kf),temp2(mini:maxi))
c$$$      endif 
c$$$      
c$$$      
c$$$cCCC      do i = mini, maxi
c$$$cCCC         tmp = tmp + chii(i,ki) * chif(i,kf) * temp3(i)
c$$$cCCC      end do
c$$$      do ns = 0, 1
c$$$         vdon(nchf,nchi,ns) = vdon(nchf,nchi,ns) + tmp
c$$$         vdon(nchi,nchf,ns) = vdon(nchf,nchi,ns)
c$$$      enddo 
c$$$      end
      
      subroutine core(ifirst,nznuci,lg,etot,chilx,minchilx,nchtop,
     >   u,ldw,vdcore,minvdc,maxvdc,npk,vdon)!,vmat,vmatp)
      use vmat_module
      use chil_module
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      integer npk(nchtop+1)
      common/smallr/ formcut,regcut,expcut,fast,match,analyticd,packed
      logical fast, match, analyticd, packed
      real !vmat(npk(nchtop+1)-1,npk(nchtop+1)),
     >   vdon(nchan,nchan,0:1)
      dimension psic(maxr),ovlp(kmax),vdcore(maxr,0:lamax),temp(maxr),
     >   fun(maxr),u(maxr)!,vmatp((npk(nchtop+1)-1)*npk(nchtop+1)/2,0:1)
c$$$      dimension chil(meshr,npk(nchtop+1)-1),minchil(npk(nchtop+1)-1)
      data pi/3.1415927/
      logical posi, positron

      rnorm = 2.0/pi
      if (npk(2)-npk(1).eq.1) then
         nchii = 1
         nchif = nchtop
      else
         nchii = nchistart(nodeid)
         nchif = nchistop(nodeid)
      endif
      
      do nch = nchii, nchif
         call getchinfo (nch, nt, lg, temp, maxpsi, e, la, na, l)
         posi = positron(na,la,nposi)
         if(posi) cycle  ! core potential are included into
                         ! Ps->Ps channels via ve(r)  and va(r) 

         nqm = npk(nch+1) - npk(nch)
         if (l .gt. ldw) then
            do ki = 1, nqm
               kii = npk(nch) + ki - 1
               do i = max(minchil(kii,1),minvdc), maxvdc
                  temp(i) = chil(i,kii,1) * vdcore(i,la) / rmesh(i,3)
               enddo 
               do 5 kf = 1, nqm
                  kff = npk(nch) + kf - 1
                  if (kff.lt.kii) go to 5
                  tmp = 0.0
                  mini = max(minvdc, minchil(kii,1),minchil(kff,1))
                  maxi = maxvdc
                  do i = mini, maxi
                     tmp = tmp + chil(i,kff,1) * temp(i)
                  end do
c$$$                  if (packed) then
c$$$                     vmatp(kii+(kff-1)*kff/2,0) = - rnorm * tmp 
c$$$     >                  + vmatp(kii+(kff-1)*kff/2,0)
c$$$                     vmatp(kii+(kff-1)*kff/2,1) = - rnorm * tmp 
c$$$     >                  + vmatp(kii+(kff-1)*kff/2,1)
c$$$                  else
                  if (npk(2)-npk(1).eq.1.or.
     >               (.not.scalapack.and.nodeid.eq.1)) then
                     vmat(kff,kii) = vmat(kff,kii) + rnorm * tmp 
                     vmat(kii,kff+1) = vmat(kii,kff+1) + rnorm*tmp
                  else 
c$$$                     if (ifirst.eq.1) then !only for the second call
                     vmat01(kff,kii) = vmat01(kff,kii)+rnorm*tmp
                     vmat01(kii,kff+1)=vmat01(kii,kff+1)+rnorm*tmp
c$$$                  endif 
                     endif 
c$$$                  endif 
 5             continue 
            end do
         endif 
C  Define the direct on shell V matrix used for analytic Born subtraction
C  in the cross program
c$$$         if (packed) then
c$$$            vdon(nch,nch,0) = vmatp(npk(nch)+(npk(nch)-1)*npk(nch)/2,0)
c$$$         else
c$$$         if (npk(2)-npk(1).eq.1.and.nodeid.eq.1) then
c$$$            vdon(nch,nch,0) = vmat(npk(nch),npk(nch))
c$$$            vdon(nch,nch,1) = vdon(nch,nch,0)
c$$$        endif
c$$$        endif
         if (scalapack.and.ifirst.eq.1) then !had ifirst commented out, not sure why. Put back for Ga with DW.
            vdon(nch,nch,0) = vmat01(npk(nch),npk(nch))
            vdon(nch,nch,1) = vdon(nch,nch,0)
         else
            vdon(nch,nch,0) = vmat(npk(nch),npk(nch))
            vdon(nch,nch,1) = vdon(nch,nch,0)
         endif
C  Core exchange only for electron scattering and if not using local exchange
C  core potentials
            
         if (nze.eq.-1.and.ntype.ne.-3.and.ifirst.ne.0) then
            do lac = 0, lactop
               do nac = lac + 1, nabot(lac) - 1
                  if (lg.eq.0.and.nch.eq.1.and.nqm.eq.1)
     >               print*,'Core N, L:',nac,lac
                  if (lac.eq.linc.and.nac.eq.ninc) then
                     elecnum = float(2 * lac)
                  else
                     elecnum = float(2 * lac + 1)
                  endif 
                  maxpsic = istoppsinb(nac,lac)
                  ec = enpsinb(nac,lac)
                  do i = 1, maxpsic
                     psic(i) = psinb(i,nac,lac)
                  end do 
                  do 20 ilt = -lac, lac, 2
                     lt = l + ilt
                     if (lt.lt.0.or.lt.gt.ltmax) go to 20
                     call cleb(2*lac,2*lt,2*l,0,0,0,c)
                     const = - rnorm * elecnum * c * c /
     >                  float(2 * l + 1)
                     do ki = 1, nqm
                        kii = npk(nch) + ki - 1
                        minfun = minchil(kii,1)
                        maxfun = maxpsic
                        do i = minfun, maxfun
                           fun(i) = psic(i) * chil(i,kii,1)
                        end do
                        call form(fun,minfun,maxfun,rpow1(1,lt),
     >                     rpow2(1,lt),minrp(lt),maxrp(lt),maxfun,
     >                     temp,i1,i2)
                        mini1 = i1
                        maxi = min(i2,maxpsic)
                        do i = mini1, maxi
                           temp(i) = temp(i) * psic(i)
                        end do 
                        do 15 kf = 1, nqm
                           kff = npk(nch) + kf - 1
                           if (kff.lt.kii) go to 15
                           sum = 0.0
                           mini = max(mini1,minchil(kff,1))
                           do i = mini, maxi
                              sum = sum + chil(i,kff,1) * temp(i)
                           end do
c$$$                           if (packed) then
c$$$                              vmatp(kii+(kff-1)*kff/2,0) = const * sum
c$$$     >                           + vmatp(kii+(kff-1)*kff/2,0)
c$$$                              vmatp(kii+(kff-1)*kff/2,1) = const * sum
c$$$     >                           + vmatp(kii+(kff-1)*kff/2,1)
c$$$                           else
                           if (npk(2)-npk(1).eq.1.or.
     >                        (.not.scalapack.and.nodeid.eq.1)) then
                              vmat(kff,kii) = vmat(kff,kii)+
     >                           const * sum
                              vmat(kii,kff+1) = vmat(kii,kff+1) +
     >                           const * sum
                           else 
                              if (ifirst.eq.1) then
                                 vmat01(kff,kii) = vmat01(kff,kii)+
     >                              const * sum
                                 vmat01(kii,kff+1) = vmat01(kii,kff+1)+
     >                              const * sum
                              endif 
                              endif 
c$$$                           endif 
 15                     continue 
                     end do
 20               continue
                     
C  Core overlaps
cCCC                  if (l.eq.lac) then
cCCCC  Multiply by the energy term converted to a.u.
cCCCcCCC                     const = - rnorm * (ec + ec - etot) / 2.0
cCCCcCCC                     const = - rnorm * (e + ec - etot) / 2.0
cCCC                     const = 0.0
cCCC                     do k = 1, nqm
cCCC                        ovlp(k) = 0.0
cCCC                        do i = minchil(k,nch), maxpsic
cCCC                           ovlp(k) = ovlp(k)+chil(i,k,nch)*psic(i)
cCCC                        end do
cCCC                     end do
cCCC                     do ki = 1, nqm
cCCC                        do kf = 1, nqm
cCCC                           vmat(kf,ki,nch,nch)=vmat(kf,ki,nch,nch)
cCCC     >                        + const * ovlp(kf) * ovlp(ki)
cCCC                        end do
cCCC                     end do
cCCC                  end if

C  End of core loops
               enddo 
            end do
         end if 
C  End of channel loop
      end do 
      end
      
C  The following routine makes the core potential
      subroutine makevdcore(vdcore,minvdc,maxvdc,nznuci,u)
      use vmat_module, only: nodeid
      include 'par.f'
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/meshrr/ meshr,rmesh(maxr,3)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      dimension vdcore(maxr), vnl(maxr), fun(maxr), u(maxr)
c$$$      character ch
c$$$      ch(i)=char(i+ichar('0'))
      
      minvdc = meshr
      maxvdc = 0
      if (nodeid.eq.1)
     >   print '('' n l of the core states'')'
      do lac = 0, lactop
         do nac = lac + 1, nabot(lac) - 1
            if (nodeid.eq.1)
     >         print '(2i2)', nac, lac
            const = float(4 * lac + 2)
            if (lac.eq.linc.and.nac.eq.ninc) then
               print*,'number of core electrons for this l:',4 * lac
               const = float(4 * lac)
               print*, 'adding a*exp(-b*r)/r, give a and b:'
               read*,a,b
               i = 1
               expr = exp(-b*rmesh(i,1))
               do while (expr.gt.1e-10)
                  vdcore(i) =  vdcore(i) + a*expr/rmesh(i,1)
                  i = i + 1
                  expr = exp(-b*rmesh(i,1))
               enddo 
            endif 
            minfun = 1
            maxfun = istoppsinb(nac,lac)
            do i = minfun, maxfun
               fun(i) = psinb(i,nac,lac) * psinb(i,nac,lac) * rmesh(i,3)
            end do
            call form(fun,minfun,maxfun,rpow1(1,0),rpow2(1,0),
     >         minrp(0),maxrp(0),meshr,vnl,i1,i2)
            call nuclear(fun,.true.,minfun,maxfun,i2,u,nznuc,vnl)
c$$$            open(42,file='corepot'//ch(nac)//ch(lac))
c$$$            do i = 1, meshr
c$$$               write(42,*) rmesh(i,1), vnl(i)
c$$$            enddo 
c$$$            close(42)
            if (i1.lt.minvdc) minvdc = i1
            if (i2.gt.maxvdc) maxvdc = i2
            do i = minvdc, maxvdc
               vdcore(i) = vdcore(i) + const * vnl(i)
            end do
         end do
      end do
c$$$      do i = minvdc, maxvdc
c$$$         write(42,*) rmesh(i,1),vdcore(i)*rmesh(i,1)
c$$$      enddo
c$$$  stop
      return
      end 

c$$$      subroutine getsplines (psi, maxpsi, rmesh, der, diffder3)
c$$$      include 'par.f'
c$$$      real x(maxr), y(maxr), b(maxr), c(maxr), d(maxr), der(0:3,0:2),
c$$$     >   diffder3(0:2,maxr), psi(maxpsi), rmesh(maxpsi)
c$$$      x(1) = 0.0
c$$$      do i = 2, maxpsi
c$$$         x(i) = rmesh(i-1)
c$$$      end do 
c$$$      do l = 0, 2
c$$$         do i = 2, maxpsi
c$$$            y(i) = psi(i-1) / x(i)**l
c$$$         enddo 
c$$$         y(1) = 2.0 * y(2) - y(3)
c$$$         if (abs(y(1)).lt.1e-3) y(1) = 0.0
c$$$         call spline(maxpsi,x,y,b,c,d)
c$$$         der(0,l) = y(1)
c$$$         der(1,l) = b(1)
c$$$         der(2,l) = 2.0 * c(1)
c$$$         der(3,l) = 6.0 * d(1)
c$$$         do i = 1, maxpsi - 1
c$$$            diffder3(l,i) = (d(i+1) - d(i)) * 6.0
c$$$         end do
c$$$         diffder3(l,maxpsi) = 0.0
c$$$      end do
c$$$      return
c$$$      end
      
      subroutine getoverlap(psi,maxpsi,rmesh,l,rk,qcut,
     >   u,minchi,chi,diffder3,der,tmp)
      real psi(maxpsi), rmesh(maxpsi), u(maxpsi), chi(maxpsi),
     >   diffder3(0:2,maxpsi), der(0:3,0:2)
      tmp = 0.0
      if (l.gt.2.or.rk.lt.qcut/2.0.or.u(1).ne.0.0) then
         do i = minchi, maxpsi
            tmp = tmp + chi(i) * psi(i)
         end do
      else if (l.eq.0) then
         ts = 0.0
         do i = 1, maxpsi
            r = rmesh(i)
            ts = ts + sin(rk*r) * diffder3(l,i)
         enddo
         tmp = (ts / rk - der(2,l)) / rk**3 + der(0,l) / rk
cCCC                        print*,ki,a*(2.0*rk)/(a*a+rk*rk)**2/tmpi,tmpi
      else if (l.eq.1) then
         tc = 0.0
         ts = 0.0
         do i = 1, maxpsi
            r = rmesh(i)
            ts = ts + sin(rk*r) * diffder3(l,i)
            tc = tc + cos(rk*r) * diffder3(0,i)
         enddo
         tmp = (ts / rk - der(2,l)) / rk**4 + der(0,l) / rk**2
     >      + (der(1,0) - (der(3,0) + tc) / rk**2) / rk**2
      else if (l.eq.2) then
         tc = 0.0
         ts = 0.0
         do i = 1, maxpsi
            r = rmesh(i)
            ts = ts + sin(rk*r)*(3.0*diffder3(l,i)/rk**2-diffder3(0,i))
            tc = tc + cos(rk*r) * diffder3(1,i)
         enddo
         tmp = (ts / rk - 3.0*der(2,l)/rk**2 + der(2,0)) / rk**3
     >      + 3.0 * der(0,l) / rk**3 - der(0,0) / rk
     >      + (der(1,1) - (der(3,1) + tc) / rk**2) / rk**2 * 3.0 /rk
      endif 
      end

      function calctail(itail,r,rkf,phif,rki,phii)
      real*8 xp, xm, dci, dsi
      
      calctail = 0.0
      if (rkf.le.0.0.or.rki.le.0.0) return
      xp = (rkf + rki) * r
      xm = (rkf - rki) * r
      phip = phif + phii
      phim = phif - phii
      if (xm.lt.0.0) then
         xm = - xm
         phim = - phim
      endif 
      if (abs(rkf/rki-1.0).lt.1e-5) then
         partm = cos(phim)
      else
         partm = cos(xm+phim)+xm*(cos(phim)*dsi(xm)+sin(phim)*dci(xm))
      endif 
      calctail = (partm
     >   -(cos(xp+phip)+xp*(cos(phip)*dsi(xp)+sin(phip)*dci(xp))))/2./r
      return
      
c$$$      if (abs(rkf/rki-1.0).lt.1e-5) then
c$$$         argp = (rkf + rki) * r + phif + phii
c$$$         rkp = 2.0 * rkf 
c$$$         rp = rkp * r
c$$$         if (rp.le.itail) return
c$$$         calctail = 0.5 * cos(phif-phii) / r + 0.5 * rkp * (
c$$$     >      sin(argp) / rp ** 2
c$$$     >      - 2.0 * cos(argp) / rp ** 3
c$$$cCCC     >      - 6.0 * sin(argp) / rp ** 4
c$$$cCCC     >      + 24.0 * cos(argp) / rp ** 5
c$$$     >      )
c$$$      else
c$$$         argp = (rkf + rki) * r + phif + phii
c$$$         argm = (rkf - rki) * r + phif - phii
c$$$         rkp = (rkf + rki)
c$$$         rkm = (rkf - rki)
c$$$         rp = rkp * r
c$$$         rm = rkm * r
c$$$         if (rp.le.itail.or.abs(rm).le.itail) return
c$$$         calctail = 0.5 * rkm * (
c$$$     >      - sin(argm) / rm ** 2  
c$$$     >      + 2.0 * cos(argm) / rm ** 3
c$$$cCCC     >      + 6.0 * sin(argm) / rm ** 4 
c$$$cCCC     >      - 24.0 * cos(argm) / rm ** 5
c$$$     >      ) - 0.5 * rkp * (
c$$$     >      - sin(argp) / rp ** 2  
c$$$     >      + 2.0 * cos(argp) / rp ** 3
c$$$cCCC     >      + 6.0 * sin(argp) / rp ** 4 
c$$$cCCC     >      -  24.0 * cos(argp) / rp ** 5
c$$$     >      )
c$$$      endif
      return
      end
      

      function calctail2(a,rlam,r,rkf,phif,rki,phii)
C  Here we replace the 1/r**2 with a * exp(- r * rlam). Doesn't work well
C  enough.
      calctail2 = 0.0
      if (rkf.le.0.0.or.rki.le.0.0) return
      argp = (rkf + rki) * r + phif + phii
      argm = (rkf - rki) * r + phif - phii
      rkp = (rkf + rki)
      rkm = (rkf - rki)
      calctail2 = a * 0.5 * exp(- rlam * r) * (
     >   (rlam * cos(argm) - rkm * sin(argm)) / (rlam**2 + rkm**2) -
     >   (rlam * cos(argp) - rkp * sin(argp)) / (rlam**2 + rkp**2))
      return
      end

C The following integrates dr from 1 to oo the tail integrals that fall-off 
C as r**(-lt-1). The input projectiles are chi(Rkr), where R=rmesh(meshr,1).
      subroutine maketail(itail,ctemp,chii,minchii,gki,phasei,li,nqmi,
     >   chif,minchif,gkf,phasef,lf,nqmf,nchf,nchi,ltmin,nqmfmax,vmatt)
c$$$      use gf_module
      include 'par.f'
      implicit double precision (a-h,o-z)
      common/meshrr/ meshr,rmesh(maxr,3)
      common/matchph/rphase(kmax,nchan),trat
      real gki(kmax),gkf(kmax),vmatt(nqmfmax,nqmfmax),chii(meshr,nqmi),r
     >   ,chif(meshr,nqmf),temp(maxr),ctemp,rphase,rmesh,aimag,calctail
      integer minchii(nqmi),minchif(nqmf)
      complex phasei(kmax),phasef(kmax)
      data pi,small/3.141592653589793238,0.001/
c$$$      character ch
c$$$      ch(i)=char(i+ichar('0'))

c$$$      vmatt(:,:) = 0.0
!      if (itail.eq.0.or.ctemp.eq.0.0) return !.or.ltmin.ne.1) return

      if (itail.lt.0) then
         lt = ltmin
         const = ctemp * (trat / rmesh(meshr,1))**lt
         do ki = 1, nqmi
            rk1 = gki(ki)
            if (rk1.lt.0.0) cycle
            do i = minchii(ki), meshr
               temp(i) = chii(i,ki)/rmesh(i,3)/rmesh(i,1)**(lt+1)
            enddo 
            do kf = 1, nqmf
               mini = max(minchii(ki), minchif(kf))
c$$$               n = meshr - mini + 1
               rk2 = gkf(kf)
               if (rk2.lt.0.0) cycle
               vmatt(kf,ki) = const *
     >              dot_product(chif(mini:meshr,kf),temp(mini:meshr))
            enddo
         enddo
      elseif (itail.eq.1) then
c$$$      if (itail.ge.0.or.ctemp.eq.0.0.or.phasei(1).ne.(1.0,0.0).
c$$$     >   or.phasef(1).ne.(1.0,0.0)) return
c$$$
c$$$c$$$      analytic = nanalytic.le.-2.and.gki(1).ge.0.0.and.gkf(1).ge.0.0
c$$$c$$$      inquire(FILE="analytic",EXIST=analytic)
c$$$c$$$      if (analytic) print*,'tail integrals only for onshell points'
c$$$
         r = rmesh(meshr,1)
c$$$      if (itail.lt.0.and.ltmin.le.-itail) then 
         if (mod(li+lf,2).eq.0) then
            lt = 2
         else
            lt = 1
         endif
         lt = ltmin !need to modify He code for ltmin
         l1 = li
         l2 = lf
         lm = lt
         IZR = 0
         IONE = 1
         scale = 1d3
!         RT10 = SQRT(10.0D0)
         RT10 = SQRT(scale)
         T1 = GAMX(IONE,2*LM)/(GAMX(IONE,L1-L2+LM+1)
     >      *GAMX(IONE,L2-L1+LM+1))
!         T2 = GAMX(IZR,L1+L2-LM+2)/(GAMX(IZR,L1+L2+LM+2)*(1.0D1**LM))
         T2 = GAMX(IZR,L1+L2-LM+2)/(GAMX(IZR,L1+L2+LM+2)*(scale**LM))
         S3 = RT10**(L2-L1-LM-1)
         T3 = (GAMX(IZR,L1+L2+2-LM)*S3)/GAMX(IZR,2*L1+3)
         S4 = RT10**(L1-L2-LM-1)
         T4 = (GAMX(IZR,L1+L2+2-LM)*S4)/GAMX(IZR,2*L2+3)
         T5 = GAMX(IONE,L2-L1+LM+1)
         T6 = GAMX(IONE,L1-L2+LM+1)
         T7 = GAMX(IONE,2*LM)
         S8 = RT10**(L1-L2-LM+1)
         T8 = S8*GAMX(IZR,2*L1+3)/GAMX(IZR,L1+L2+LM+2)
         T8 = T8/GAMX(IONE,L1-L2+LM+1)
         S9 = RT10**(L2-L1-LM+1)
         T9 = S9*GAMX(IZR,2*L2+3)/GAMX(IZR,L1+L2+LM+2)
         T9 = T9/GAMX(IONE,L2-L1+LM+1)
         T10 = GAMX(IONE,L1-L2-LM+1)
         T11 = GAMX(IONE,L2-L1-LM+1)
         TWOL = DBLE(2**LM)

c$$$         open(42,file='chif'//ch(lf))
c$$$         write(42,'("#           ",200i12)') (minchif(k),k=1,nqmf)
c$$$         do i = 1, meshr
c$$$            write(42,'(200e12.4)') rmesh(i,1),(chif(i,k)/rmesh(i,3)
c$$$C     >         ,k=1,nqmf)
c$$$     >            /rmesh(i,1)**(lt+1),k=1,nqmf)
c$$$         enddo
c$$$         close(42)
         
c$$$         open(42,file='chif'//ch(lf))
c$$$         write(42,'("#           ",200i12)') lf,(minchif(k),k=1,1)
c$$$         do i = 1, meshr
c$$$            write(42,'(200e12.4)') rmesh(i,1),(chif(i,k)/rmesh(i,3),
c$$$     >         sin(gkf(k)*rmesh(i,1)+lf*pi/2d0*(-1)**lf),k=1,1)
c$$$         enddo
c$$$         close(42)
         do ki = 1, nqmi
            rk1 = gki(ki)
C  Have not got the following t work for rk1=0.0
            if (rk1.lt.0.0) cycle
c$$$            if (rk1.gt.0.0.and.abs(aimag(phasei(ki))).le.small) then
c$$$               n1 = rk1*rmesh(meshr,1)/2d0/pi
c$$$               d =  rk1*rmesh(meshr,1) - 2d0*n1*pi 
c$$$               phi1 = asin(chii(meshr,ki)/rmesh(meshr,3))
c$$$     >            - d
c$$$               print*,'li,rki,phii,sin(kr),chii:',li,rk1,phi1,
c$$$     >         sin(rk1*rmesh(meshr,1)+d),chii(meshr,ki)/rmesh(meshr,3)
            do i = minchii(ki), meshr
               temp(i) = chii(i,ki)/rmesh(i,3)/rmesh(i,1)**(lt+1)
            enddo 
            do kf = 1, nqmf
               mini = max(minchii(ki), minchif(kf))
               n = meshr - mini + 1
               rk2 = gkf(kf)
C  Have not got the following t work for rk2=0.0
               if (rk2.lt.0.0) cycle
c$$$                  if(rk2.gt.0.0.and.abs(aimag(phasef(kf))).le.small)then
c$$$               it = it + 1
               ztormax = 0.0
               if (n.gt.0) ztormax =
     >              dot_product(chif(mini:meshr,kf),temp(mini:meshr))
c$$$     >                  ddot(n,chif(mini,kf),1,temp(mini),1)
C
C ****  THIS SECTION EVALUATES THE DEFINITE INTEGRAL FROM
C ****  ZERO TO INFINITY. TWO DIFFERENT REPRESENTATIONS OF
C ****  THE HYPERGEOMETRIC FUNCTION ARE USED TO SPEED UP
C ****  THE CALCULATIONS. THE HYPERGEOMETRIC SERIES IS
C ****  NOTORIOUSLY SLOW TO CONVERGE WHEN THE ARGUMENT IS
C ****  CLOSE TO ONE.
               XRT = MAX(0.85D0,1.0D0-1.0D0/DBLE(L1+L2+1))
C
C ****  FIRST SPECIAL CASE, RK1 = RK2
C
               TK = (RK1/RK2) - 1.0D0
               IF((ABS(TK).LT.1.0D-8)) THEN
                  YINF = ((RK1*0.5D0)**LM)*T1*T2
c$$$                  print*,'yinf,rk1,lm,t1,t2',
c$$$     >                 yinf,rk1,lm,t1,t2
C
C ****  NEXT CASES, RK1 < RK2  OR RK1 > RK2
C
               ELSE
C
C ****  RK1 < RK2
C
                  WA = 0.5D0*DBLE(L1+L2+2-LM)
                  IF(RK1.LT.RK2) THEN
                     TK = RK1/RK2
                     ARGF = TK*TK
                     TK11 = (TK**(L1+1))
                     VK = RK2**LM
                     IF(ARGF.LE.XRT) THEN
                        WB = 0.5D0*DBLE(L1-L2-LM+1)
                        WC = 0.5D0*DBLE(2*L1+3)
                        NTERM = -1
                        F21 = FDHY(WA,WB,WC,NTERM,ARGF)
                        YINF = F21*VK*TK11*T3/(T5*TWOL)
c$$$                        print*,'YINF,F21,VK,TK11,T3,T5,TWOL',
c$$$     >                       YINF,F21,VK,TK11,T3,T5,TWOL
                        if (rk1.eq.0.0) YINF = F21*VK*T3/(T5*TWOL)
                     ELSE
                        ARGF = 1.0D0 - TK*TK
                        WL = (-ARGF)**LM
                        NTERM = LM - 1
                        UB = 0.5D0*DBLE(L1-L2-LM+1)
                        UC = DBLE(1-LM)
                        F2X = T7*T8*FDHY(WA,UB,UC,NTERM,ARGF)
                        VA = 0.5D0*DBLE(L1+L2+LM+2)
                        VB = 0.5D0*DBLE(L1-L2+LM+1)
                        VC = DBLE(LM)
                        F2Y = WL*FDHY2(VA,VB,VC,ARGF)/(T3*T10)
                        YINF = TK11*VK*T3*(F2X - F2Y)/(TWOL*T5)
c$$$                        print*,'YINF,TK11,VK,T3,F2X,F2Y,TWOL,T5',
c$$$     >                       YINF,TK11,VK,T3,F2X,F2Y,TWOL,T5
                     END IF
C     
C     ****  RK2 < RK1
C
                  ELSE
                     TK = RK2/RK1
                     ARGF = TK*TK
                     TK11 = (TK**(L2+1))
                     VK = RK1**LM
                     IF(ARGF.LE.XRT) THEN
                        WB = 0.5D0*DBLE(L2-L1-LM+1)
                        WC = 0.5D0*DBLE(2*L2+3)
                        NTERM = -1
                        F21 = FDHY(WA,WB,WC,NTERM,ARGF)
                        YINF = F21*VK*TK11*T4/(T6*TWOL)
c$$$                        print*,'YINF,F21,VK,TK11,T4,T6,TWOL',
c$$$     >                       YINF,F21,VK,TK11,T4,T6,TWOL
                        if (rk2.eq.0.0) YINF = F21*VK*T4/(T6*TWOL)
c$$$                        if (nchf.eq.6.and.nchi.eq.2.and.ki.eq.1)
c$$$     >                       print*,'kf,argf,xrt,yinf-ztormax:',
c$$$     >                       kf,argf,xrt,yinf*pi/2.0-ztormax
                     ELSE
                        ARGF = 1.0D0 - TK*TK
                        WL = (-ARGF)**LM
                        NTERM = LM - 1
                        UB = 0.5D0*DBLE(L2-L1-LM+1)
                        UC = DBLE(1-LM)
                        F2X = T7*T9*FDHY(WA,UB,UC,NTERM,ARGF)
                        VA = 0.5D0*DBLE(L1+L2+LM+2)
                        VB = 0.5D0*DBLE(L2-L1+LM+1)
                        VC = DBLE(LM)
                        F2Y = WL*FDHY2(VA,VB,VC,ARGF)/(T4*T11)
                        YINF = TK11*VK*T4*(F2X - F2Y)/(TWOL*T6)
c$$$                        print*,'YINF,TK11,VK,T4,F2X,F2Y,TWOL,T6',
c$$$     >                       YINF,TK11,VK,T4,F2X,F2Y,TWOL,T6
c$$$                        print*,'kf,argf,xrt,yinf-ztormax:',
c$$$     >                       kf,argf,xrt,yinf*pi/2.0-ztormax
                     END IF
                  END IF
               endif
c$$$                     write(50+ki,'(1p,4e12.3)') rk1,rk2,yinf*ctemp,
c$$$     >                  ztormax*ctemp
c$$$                     if (kf.eq.nqmf) write(50+ki,*)
                     
               yinf = yinf * pi / 2.0
c$$$                     if (lt.eq.3.and.(li.eq.0.or.lf.eq.0))
c$$$     >                  print*,'ztormax,yinf,rk1,rk2,li,lf,lt:',
c$$$     >                  ztormax,yinf,rk1,rk2,li,lf,lt
c$$$                     tail(kf,ki) = (YINF - ztormax
!                     if (ki.eq.1.and.kf.eq.1) print*,
!     >                    'yinf,ztormax,phasei,phasef,ctemp:',
!     >                    yinf,ztormax,phasei(ki),phasef(kf),ctemp
c$$$               if (nchf.eq.52.and.nchi.eq.43) print*,'tail:',
c$$$     >            yinf,ztormax,ctemp
               vmatt(kf,ki) = (YINF - ztormax
!     >                  *real(phasei(ki))*real(phasef(kf))
     >              ) * ctemp
c$$$                     if (abs((YINF-ztormax)/(YINF+ztormax+1e-20)).lt.
c$$$     >                  1e-4) tail(kf,ki) = 0d0
c$$$               print*,'lt,rk1,rk2,yinf,ztormax,ctemp:',
c$$$     >              lt,rk1,rk2,yinf,ztormax,ctemp
            enddo
         enddo
      endif
      
c$$$      else if (itail.eq.-1) then
c$$$         do ki = 1, nqmi
c$$$            do kf = 1, nqmf
c$$$               if (gkf(kf).gt.0.0.and.gki(ki).gt.0.0) then
c$$$                  xf = real(phasef(kf))
c$$$                  xi = real(phasei(ki))
c$$$                  yf = aimag(phasef(kf))
c$$$                  yi = aimag(phasei(ki))
c$$$                  dkf = gkf(kf)
c$$$                  dki = gki(ki)
c$$$                  drmax = r
c$$$                  lt = 1
c$$$                  ffggpart = ffgg(lf,dkf,li,dki,drmax,lt,xf,yf,xi,yi) 
c$$$                  fggfpart = fggf(lf,dkf,li,dki,drmax,lt,xf,yf,xi,yi)
c$$$cCCC                  appt = calctail(itail,r,gkf(kf),
c$$$cCCC     >               rphase(kf,nchf),gki(ki),rphase(ki,nchi))
c$$$cCCC                  if (abs(ffggpart).gt.100.0*abs(appt)) then
c$$$cCCC                     print*,ffggpart,appt,dkf,dki,lf,li,
c$$$cCCC     >                  'ffgg,appt,kf,ki,lf,li'
c$$$cCCC                  endif 
c$$$                  tail(kf,ki) = (ffggpart + fggfpart) * ctemp
c$$$                  if (abs(tail(kf,ki)).gt.1e0) then
c$$$                     print*,'nchf,nchi,kf,ki:',nchf,nchi,kf,ki
c$$$                     print*,'ffggpart, fggfpart:',ffggpart, fggfpart
c$$$c$$$                     print*,'maxfun,i2,meshr,lt,ctemp',
c$$$c$$$     >                  maxfun,i2,meshr,lt,ctemp
c$$$                  endif
c$$$               endif 
c$$$            end do
c$$$         end do
c$$$      else 
c$$$         do ki = 1, nqmi
c$$$            do kf = 1, nqmf
c$$$               tail(kf,ki) = calctail(itail,r,gkf(kf),
c$$$     >            rphase(kf,nchf),gki(ki),rphase(ki,nchi)) * ctemp
c$$$            enddo
c$$$         enddo
c$$$      endif
cCCC      print*,'Number of tail integral calculations:',it
      return
      end
      
      subroutine kgrid(ispeed,nkor,skor,etot,gk,wk,weightk,nbnd,nqm,lg,
     >   nchtop,nchopt,npk,nold,ndumm,luba,ifirst,npsbnd,nnbtop,hlike,
     >   uba,theta,nodes,nodeid)
      use gf_module !propagates analytic
      include 'par.f'
      complex wk(kmax*nchan), phase, sigc
      dimension nkor(5,0:lmax,2),skor(5,0:lmax,2),psin(maxr)
      dimension nk(10),sk(0:10),gridk(kmax),weightk(kmax),psi(maxr)
      real*8 xx(kmax), ww(kmax),wf(2*kmax),dstart,dstop,alpha,beta,p
      integer iwf(2*kmax),nbnd(0:lmax)
      dimension gk(kmax,nchan), gfixed(kmax), wfixed(kmax),npk(nchan+1)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax),itail
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      dimension reg(maxr),ucentr(maxr),tmpk(kmax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /double/njdouble,jdouble(22)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast,hlike,usetrapz,uba(nchan),positron,pos,lprint
      common/cont/ ce(ncmax),cint(ncmax),ncstates,energy
      character chan(knm)*3
      common /charchan/ chan
      common /chanen/ enchan(knm)

      data ucentr/maxr*0.0/

C     Added by alex
c$$$      logical analytic

c$$$      analytic = lg.le.(-nqm-2) !Use NQM in ccc.in as a switch.
      nanalytic = nqm
      pi = atan(1.0)*4.0

c$$$      nbox = 1

c$$$      inquire(FILE="analytic",EXIST=analytic)  
c$$$
c$$$      if (analytic) then
c$$$         open(42,file = "analytic")
c$$$         read(42,*) nbox
c$$$         close(42)
c$$$      endif
C     End added by Alex

C  NBAD stores the number of channels |phi(n)> such that
C  Int dk <phi(n)|k><k|phi(n)> .ne. 1
      nbad = 0
      lprint = nodeid.eq.1.and.lg.le.latop
      if (lprint) print '(''Kgrid quadrature set'')'
      nch = 1
      npk(nch) = 1
      call getchinfo (nch, ntmp, lg, psi, maxpsi, ea, la, na, li)
      if (nch.eq.0) stop 'NCH = 0 in KGRID'
C  The following is a loop of the form REPEAT ... UNTIL(condition)
 10   e = etot - ea ! Ry
      EeV = e * 13.6058
      posfac = 1.0
      if (positron(na,la,npos)) posfac = 2.0
      e = e * posfac
      if (e.ge.0.0) then
         rk = sqrt(e)
      else
         rk = - sqrt(-e)
      end if
      analytic = nanalytic.le.-2.and.e.ge.aenergyswitch !aenergyswitch is in modules.f
c$$$      print'(a3, " channel energy (eV) and k:",1p,2e15.6)', chan(ntmp),
c$$$     >   EeV,rk
c$$$      if (analytic.and.EeV.ge.0.0.and.EeV.lt.1e-3) then !+ve energies for extract_J
C analytic tail integrals have not been implemented for zero energy
c$$$      if (analytic.and.abs(EeV).lt.1e-4*(2.0*li+1.0).and.itail.ge.0)then
c$$$      if (analytic.and.abs(EeV).lt.1e-3*(2.0*li+1.0).and.ea.lt.0.0) then !.and.nch.gt.1)then         
      if (analytic.and.ea.lt.0.0.and.-e/ea.lt.1e-4) then !.and.nch.gt.1)then    
         print'("CAUTION: for NCH =",i3," L =",i2,
     >   " setting on-shell E (eV) to zero:",1p,e10.2)',nch,li,EeV
         rk = 0d0
         e = 0d0
      endif 
      if (analytic) print*,'Analytical GF used for channel:',nch
      nqk = 0
      lset = li
      nquad = 0
      do n = 1, 4
         nquad = nquad + max(nkor(n,lset,ispeed),0)
      enddo
C  Treat large target L using the UBA. Use negative LUBA so as not to 
C  make a mistake in the input file.        
      uba(nch) = la .ge. -luba .and. luba .lt. 0
      if (nqm.eq.1.or.uba(nch).or.nquad.eq.0) then
         print*,'nqm,uba(nch),nquad:',nqm,uba(nch),nquad
         go to 20
      endif 
c$$$      npoints = npointsl(lset)
c$$$      width = widthl(lset)
c$$$      midnp = abs(midnpl(lset))
c$$$      usetrapz = midnpl(lset).lt.0
c$$$      endk = endkl(lset)
c$$$      nendk = nendkl(lset)
c$$$      endp = endpl(lset)
C 1,3,   20,   0.8,    16,     2.5,    6,    6.0,   12,   0.09  
C     npoints, endk, npoints2, endk2, nendk, endp, midnp, width 
      npoints = abs(nkor(1,lset,ispeed))
      endk = skor(1,lset,ispeed)
      npoints2 = abs(nkor(2,lset,ispeed))
      endk2 = skor(2,lset,ispeed)
      nendk = abs(nkor(3,lset,ispeed))
      endp = skor(3,lset,ispeed)
      width = skor(4,lset,ispeed)
c$$$      midnp = abs(nkor(4,lset,ispeed))
      midnp = nkor(4,lset,ispeed)
c$$$      if (midnp.gt.0.and.rk.gt.endk) midnp = -midnp !sets -ve midnp if singularity is not in 1st interval
      if (midnp.lt.0.and.rk.gt.endk2) midnp = -midnp !allows +ve midnp at high energies
c$$$      if (rk.lt.0.0) then       !closed channel
c$$$         npoints2 = 2
c$$$         nendk = 2
c$$$      endif
c$$$  if (midnp.lt.0.and.rk.lt.width) midnp = -midnp !allows +ve midnp near thresholds
c$$$      usetrapz = nkor(4,lset,ispeed).lt.0
c$$$      if (width.lt.0.0) then
c$$$         dstart = 0d0
c$$$         nt = midnp
c$$$c$$$         nt = npoints + midnp
c$$$         nwf = 2 * nt
c$$$         niwf = 2 * nt
c$$$         dstop  = dble(-2.0 * width)
c$$$c$$$         dstop  = dble(endk)
c$$$         call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,
c$$$     >      0,nwf,wf,niwf,iwf,ier)
c$$$         do i = 1, nt
c$$$            gfixed(i) = real(xx(i))
c$$$            wfixed(i) = real(ww(i))
c$$$         enddo
c$$$      end if
      
      if (lprint)
     >   print '(''interval   points    start       stop          '//
     >   'test'')'
      enk = endk
 30   sumi=0.0
      sumip = 0.0
      nqk=0
      
      nk(1) = npoints
      sk(1) = endk
      nk(2) = npoints2
      sk(2) = endk2         
      mint = 3
      nk(mint) = nendk
      gridk(:) = 1.0

C Define the intervals and the number of points in each interval
c$$$      if (.not.analytic) call makeints(mint,sk,nk,rk,npoints,width,
      call makeints(mint,sk,nk,rk,npoints,width,
     >   midnp,enk,nendk,endp,npoints2,endk2)
c$$$      if (.not.analytic) then
      if (midnp.lt.0) then
         dstop = 0d0
         do j=1,mint-1
            nt=nk(j)
c$$$            if (mod(nt,2).ne.0) stop 'expect nt to be even'
            nwf=2*nt
            niwf=2*nt
            dstart = dstop
            dstop  = dble(sk(j))
            if (dstart.lt.rk.and.rk.lt.dstop) then
               if (rk.gt.dstop*0.9) then
                  if (lprint)
     >               print*,'rk, new dstop:',rk,dstop*1.1
                  dstop = dstop * 1.1
                  sk(j) = dstop
               endif 
               if (rk.lt.dstart*1.1) then
                  if (lprint)
     >               print*,'new dstart, rk:',dstart*0.9,rk
                  dstart = dstart*0.9
                  sk(j-1) = dstart
               endif 
c$$$               call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,
c$$$     >            0,nwf,wf,niwf,iwf,ier)
               npts = nt
               if (nkor(j,lset,ispeed).lt.0) npts = 2
               d = (dstop-dstart)/nt*npts
               do n = 1, nt/npts
                  call cgqf(npts,xx((n-1)*npts+1),ww((n-1)*npts+1),1,
     >               0d0,0d0,(n-1)*d+dstart,d*n+dstart,
     >               0,nwf,wf,niwf,iwf,ier)
               enddo

               test = 1.0
c$$$               if (nkor(j,lset,ispeed).lt.0) then !.and.rk.gt.dstop) then
c$$$c$$$               call trapez(dstart,dstop,nt,xx,ww)
c$$$                  print*,'Using Simpson rule'
c$$$                  call simpson(dstart,dstop,nt,xx,ww)
c$$$               endif
               startk = dstart
               stopk = dstop
              
               if (startk.eq.0.0.or.(stopk-rk)/(rk-startk).lt.40.0) then
                  dk = (stopk-startk)/nt/10
                  call getstopk(e,rk,startk,stopk,dk,test,xx,ww,nt)
                  dstop = stopk
                  sk(j) = stopk
c$$$                  print*,'stopk:',stopk,(stopk-rk)/(rk-startk)
               else 
                  dk = startk/nt/10
                  call getstartk(e,rk,startk,stopk,dk,test,xx,ww,nt)
                  dstart = startk
                  sk(j-1) = startk
c$$$                  print*,'startk:',startk,(stopk-rk)/(rk-startk)
               endif
            endif
         enddo 
      endif 
      if (width .le. 0.0.and.rk.lt.2.0*abs(width)) then
         print*,'WIDTH:', width
c$$$         if (.not.analytic)
         print*, 'Entering improvek'
         call improvek(e,sk(1),nk(1),gfixed,wfixed,endf,enk)
         sk(1) = endf
      endif 
      dstop = 0d0
C  Obtain Gaussian knots and weights within each interval
      do j=1,mint-1
         nt = nk(j)
c$$$         if (mod(nt,2).ne.0) stop 'expect nt to be even'
         nqk=nqk+nt
         nwf=2*nt
         niwf=2*nt
         if (nt.ge.0) then
C  Can use the following two lines here and in the i loop to make kgrid
C  symmetric about E rather than sqrt(E). A suitable change must be made
C  in the MAKEINTS routine.
c$$$            dstart = dstop**2
c$$$            dstop  = dble(sk(j))**2
            dstart = dstop
            dstop  = dble(sk(j))
c$$$            print*,'j,dstart,dstop:',j,dstart,dstop
c$$$            call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,
c$$$     >         0,nwf,wf,niwf,iwf,ier)
            if (nkor(j,lset,ispeed).gt.0) then
               npts = nt
c$$$               if (lprint) print*,'Using n-point quadratures',
c$$$     >            nt,nkor(j,lset,ispeed)
            else
               if (mod(nkor(j,lset,ispeed),2).eq.0) then
                  npts = 2
                  if (lprint) print*,'Using 2-point quadratures'
               else
                  npts = nt
               endif
            endif
            d = (dstop-dstart)/nt*npts
            do n = 1, nt/npts
               call cgqf(npts,xx((n-1)*npts+1),ww((n-1)*npts+1),1,
     >            0d0,0d0,(n-1)*d+dstart,d*n+dstart,
     >            0,nwf,wf,niwf,iwf,ier)
            enddo
c$$$            do i = 1, nt
c$$$               print*, i, xx(i), ww(i)
c$$$            enddo
c$$$            dstart = sqrt(dstart)
c$$$            dstop  = sqrt(dstop)

C  Can try using the trapezoidal rule, but it is not as good.
c$$$            if (j.ge.2.and.rk.gt.dstart.and.rk.lt.dstop.and.usetrapz)
c$$$     >         then
c$$$               call trapez(dstart,dstop,nt,xx,ww)
c$$$            print*,'Will use trapezoidal rule around the on-shell point'
c$$$            endif
c$$$            if (nkor(j,lset,ispeed).lt.0) then !.and.rk.gt.dstop) then
c$$$c$$$               call trapez(dstart,dstop,nt,xx,ww)
c$$$c$$$               print*,'Using trapezoidal rule'
c$$$               call simpson(dstart,dstop,nt,xx,ww)
c$$$               print*,'Using Simpson rule, dstart,dstop:',dstart,dstop
c$$$               do i = 1, nt
c$$$                  if (xx(i).lt.rk.and.rk.lt.xx(i+1)) 
c$$$     >               print*,i,xx(i),rk,xx(i+1)
c$$$               enddo 
c$$$            endif
            
c$$$            if (midnp.lt.0.and.dstart.lt.rk.and.rk.lt.dstop) then
c$$$               test = 1.0
c$$$               startk = dstart
c$$$               stopk = dstop
c$$$               if (startk.eq.0.0.or.stopk-rk.lt.rk-startk) then
c$$$                  dk = (stopk-startk)/nt/10
c$$$                  call getstopk(e,rk,startk,stopk,dk,test,xx,ww,nt)
c$$$                  dstop = stopk
c$$$                  sk(j) = stopk
c$$$               else 
c$$$                  dk = startk/nt/10
c$$$                  call getstartk(e,rk,startk,stopk,dk,test,xx,ww,nt)
c$$$                  toll = 1e-4
c$$$                  do while (abs(test).gt.toll)
c$$$                     startkold = startk
c$$$                     if (test.gt.0.0) then
c$$$                        startk = startk - dk
c$$$                     else
c$$$                        startk = startk + dk
c$$$                     endif
c$$$                     do n = 1, nt
c$$$                        xx(n) = (xx(n)-startkold)*(stopk-startk)/
c$$$     >                     (stopk-startkold) + startk
c$$$                        ww(n) = ww(n) * (stopk-startk)/(stopk-startkold)
c$$$                        if (xx(n).lt.rk) nmin = n
c$$$                     enddo 
c$$$                     testold = test
c$$$                     test = 2.0/rk*
c$$$     >                  (atanh(startk/rk)-acoth(stopk/rk))
c$$$                     do n = 1, nt
c$$$                        test = 2.0 * ww(n)/(e - xx(n)**2) + test
c$$$                     enddo 
c$$$                     print*,test,dk,startk,stopk,nmin
c$$$                     if (test*testold.lt.0.0) dk = dk / 2.0
c$$$                  enddo 
c$$$                  dstart = startk
c$$$                  sk(j-1) = startk
c$$$               endif 
c$$$  endif

c$$$           if (abs(dstart+dstop-2.0*rk).lt.1e-10) then
c$$$               h = (rk-dstart)/(nt/2-1)
c$$$               do i = 1, nt/2
c$$$                  xx(i) = dstart + (i-1)*h
c$$$                  ww(i) = h*(2.0*mod(i+1,2)+2.0)/3.0
c$$$                  print*,i,xx(i),ww(i)
c$$$               enddo
c$$$               ww(1) = ww(1)/2.0
c$$$               ww(nt/2) = ww(nt/2)/2.0
c$$$               do i = nt/2+1, nt
c$$$                  xx(i) = rk + (i-nt/2-1)*h
c$$$                  ww(i) = h*(2.0*mod(i,2)+2.0)/3.0
c$$$                  print*,i,xx(i),ww(i)
c$$$               enddo
c$$$               ww(nt/2+1) = ww(nt/2+1)/2.0
c$$$               ww(nt) = ww(nt)/2.0
c$$$               tmp = 0.0
c$$$               do i = 1, nt
c$$$                  tmp = tmp + xx(i)**3*ww(i)
c$$$               enddo
c$$$               print*,'Simpson rule test:',tmp,(dstop**4-dstart**4)/4.0
c$$$            endif
            do i=nqk-nk(j)+1,nqk
               gridk(i)=xx(i-nqk+nk(j))
               weightk(i)=ww(i-nqk+nk(j))
               if (gridk(i).eq.0.0) weightk(i) = 0.0
               ecmn = gridk(i) ** 2
               eta = 0.0
               if (la.eq.li.and.nze.eq.-1.and.nodes.eq.1) then
                  call regular(la,ecmn,eta,ucentr,cntfug(1,la),-1,
     >               rmesh(1,1),meshr,jdouble,njdouble,regcut,expcut,
     >               reg,jstart,jstop,phase,sigc)
                  call getprod(psi,maxpsi,jstart,reg,la,rmesh,
     >               meshr,ntmp,tmp)
                  sumi = sumi + tmp * weightk(i) * 2.0 / pi
               endif 
            end do
            tmpold = tmp
            if (lprint) then
               if (dstart.lt.rk.and.rk.lt.dstop) then
                  print '(i5,i10,f12.6,'' *'',f10.6,f13.5)', j,nk(j),
     >               dstart, dstop, sumi - sumip
               else 
                  print '(i5,i10,2f12.6,f13.5)', j,nk(j),dstart,dstop,
     >               sumi - sumip
               end if
            endif 
            sumip = sumi
         end if 
      end do 
C  Here we have the last interval
      nt=nk(mint)
      if (nt.gt.0) then
         if (dstop.lt.rk) then
            print*,'The last interval must start > than RK',nt,rk,dstop
            stop 'The last interval must start > than RK'
         end if  
         nqk=nqk+nt
         p=dble(sk(mint))
         if (p.gt.0.0) then
            dstart=0d0
            dstopp = dstop
            dstop=dstop**(1d0-p)
         else
            dstart = dstop
            dstopp = dstop
            dstop = -p
         endif
         nwf=2*nt
         niwf=2*nt
         call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,0,nwf,wf,niwf,iwf,
     >      ier)
         if (ier.ne.0) then
            print*,'KGRID IER, nt, dstart, dstop:',ier,
     >        nt, dstart, dstop
            error stop "KGRID IER not zero"
         endif
         jj = nqk - nt
         do j=nt,1,-1
            if (p.gt.0.0) then
               jj=nqk-j+1
               gridk(jj)=xx(j)**(1d0/(1d0-p))
               weightk(jj)=ww(j)/(p-1d0)*dble(gridk(jj))**p
            else
               jj = jj + 1
               gridk(jj)=xx(nt+1-j)
               weightk(jj)=ww(nt+1-j)
            endif
            ecmn = gridk(jj) ** 2
            if (la.eq.li.and.nze.eq.-1.and.nodes.eq.1) then
               eta = 0.0
               call regular(la,ecmn,eta,ucentr,cntfug(1,la),-1,
     >            rmesh,meshr,jdouble,njdouble,regcut,expcut,
     >            reg,jstart,jstop,phase,sigc)
               call getprod(psi,maxpsi,jstart,reg,la,rmesh,
     >            meshr,ntmp,tmp)
               if (abs(tmp).gt.abs(tmpold)*1.5) then
                  endkt = gridk(jj)
                  tmpold = tmp
               endif 
               sumi = sumi + tmp * weightk(jj) * 2.0 / pi
            endif 
         end do 
c$$$         if (endkt.gt.enk.and.enk.lt.5.9.and.abs(sumi-1.0).gt.0.1.and.
c$$$     >      nold.eq.-2) then
c$$$            enk = min(endkt + 0.2, 6.0)
c$$$            go to 30
c$$$         endif
         if (lprint) then
            if (p.gt.0.0) then
               print '(i5,i10,f12.6,''       oo'',f16.5)',mint,nk(mint),
     >            dstopp, sumi - sumip
               print'(''fall off power:'',f5.1,23x,''= '',f7.5)',
     >            sk(mint),sumi
            else
               print '(i5,i10,2f12.6,f13.5)', mint,nt,dstart,dstop,
     >            sumi - sumip
               print'(''fall off power:'',f5.1,23x,''= '',f7.5)',
     >            sk(mint),sumi
            endif
            print*,'last n, k:',jj,gridk(jj)
         endif   

c$$$         if (abs(sumi - 1.0).gt.1e-2.and.la.eq.li) nbad = nbad + 1
         if (abs(sumi - 1.0).gt.1e-2) nbad = nbad + 1
      end if
c$$$      close(60+nch)
      if (nqk+1.gt.kmax) then
         print*,'NQK + 1 > KMAX:',nqk+1,kmax
         stop 'Increase KMAX'
      end if 

 20   continue
      if (ispeed.eq.2) go to 40
C  Check that the integration rule will handle the principle value singularity
      j=1
      if (e.gt.0.0) then
         if (nqk.gt.0) then
            sum=0.0
            nt=0
            do while (sk(j).lt.rk-0.01.and.j.lt.mint)
               sum = 2.0 * atanh(sk(j)/rk)/rk
               nt=nt+nk(j)
               j=j+1
            end do
            j=mint-1
            tsum=0.0
            do while (j.ge.1.and.sk(j).gt.rk+0.01)
               tsum = - 2.0 * acoth(sk(j)/rk) / rk
               j=j-1
            end do
            sum=sum+tsum
            im=nt+1
            do while (gridk(im).le.(sk(j+1)).and.im.lt.kmax)
               gftmp = 1.0/(e-gridk(im)*gridk(im))
               if (abs(e-gridk(im)*gridk(im))/e.lt.1e-4) gftmp=0.25/e
               sum=sum + 2.0*weightk(im)*gftmp
               im=im+1
            end do
         else
            sum = 0.0
         endif 
         if (lprint)
     >      print '(''State, NCH, NST, NA, LA, L, K, EA:'',a4,2i4,
     >      3i3,1p,2e13.4)',chan(ntmp),nch,ntmp,na,la,li,rk,ea
      else
         if (nqm.gt.1) then
            sum = 0.0
            do i = 1, nk(1) + nk(2)
               sum = sum + weightk(i) * 2.0 / (e - gridk(i) ** 2)
            enddo
            print*,'Test of closed channel integral:',- 2.0 / sqrt(-e)*
     >         atan(sk(2)/sqrt(-e)) / sum, sum
         endif 
         if (lprint)
     >      print'(''State, NCH, NST, NA, LA, L, E:    '',a4,2i4,3i3,
     >      1p,e13.4,''    closed'')',chan(ntmp),nch,ntmp,na,la,li,e
      end if 
      if (nch.eq.1.and.analytic) then
         if (nqk.eq.-nanalytic) then
            nbox = 0
         else
            nbox = 1
         endif
         print*,"NBOX:", nbox!, nqk, nanalytic
      endif
      if (nqm.le.0.and.lg.le.1.and.ifirst.ne.0.and.hlike
     >   .and.theta.ne.0.0.and.nze.eq.-1.and.nodes.eq.1) then
         nchp = 1
         call getchinfo (nchp, ntm, lg, psin, mpsn, eat, lat,nt,lit)
         do while (nchp.ne.0)
            if (lat.eq.li) then
c$$$               call getchnl(chan(nchp),n,l,nc)
               print'(2x,a3,$)', chan(ntm)
            endif 
            nchp = nchp + 1
            call getchinfo (nchp, ntm, lg, psin, mpsn, eat, lat,nt,lit)
         enddo
         print*
         nchp = 1
         call getchinfo (nchp, ntm, lg, psin, mpsn, eat, lat,nt,lit)
         do while (nchp.ne.0)
            if (lat.eq.li) then
               sumi = 0.0
               eta = 0.0
C$doacross LOCAL(ecmn,jstart,jstop,tmp,reg,phase,sigc,i)
C$& SHARE(sumi,pi,mpsn,ntm,la,eta,weightk,rmesh,jdouble)
C$& SHARE(njdouble,regcut,expcut,cntfug,ucentr)
C$& SHARE(gridk,psin,lat,meshr,tmpk) 
C$PAR DOALL PRIVATE(ecmn,jstart,jstop,tmp,reg,phase,sigc,i)
C$PAR& SHARED(sumi,pi,mpsn,ntm,la,eta,weightk,rmesh,jdouble)
C$PAR& SHARED(njdouble,regcut,expcut,cntfug,ucentr)
C$PAR& SHARED(gridk,psin,lat,meshr,tmpk) SCHEDTYPE(SELF(1))
C  The following directives are for the IBM
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& PRIVATE(ecmn,jstart,jstop,tmp,reg,phase,sigc,i)
C$OMP& SHARED(sumi,pi,mpsn,ntm,la,eta,weightk,rmesh,jdouble)
C$OMP& SHARED(njdouble,regcut,expcut,cntfug,ucentr)
C$OMP& SHARED(gridk,psin,lat,meshr,tmpk) 
               do i = 1, nqk
                  ecmn = gridk(i) ** 2
                  call regular(la,ecmn,eta,ucentr,cntfug(1,la),-1,
     >               rmesh,mpsn,jdouble,njdouble,regcut,expcut,
     >               reg,jstart,jstop,phase,sigc)
                  call getprod(psin,mpsn,
     >               jstart,reg,lat,rmesh,meshr,ntm,tmpk(i))
c$$$                  sumi = sumi + tmp * weightk(i) * 2.0 / pi
               enddo
C$OMP END PARALLEL DO
               do i = 1, nqk
                  sumi = sumi + tmpk(i) * weightk(i) * 2.0 / pi
               enddo 
               print '(f5.2,$)',sumi
               if (abs(sumi - 1.0).gt.1e-2) nbad = nbad + 1
            endif
            nchp = nchp + 1
            call getchinfo (nchp, ntm, lg, psin, mpsn, eat, lat,nt,lit)
            call update(6)
         enddo 
         print*
      endif 

 40   continue
!      else !analytic
c$$$      if (analytic) then
c$$$        nqk = -nqm - nbnd(lset) !nbnd get's added below
c$$$         sum = 0.0
c$$$  endif

c$$$      sum = 0.0 ! temporary fix for rk in the k-grid
      
      do i = 1, nqk + 1
         kp = npk(nch) + i - 1
         if (i.eq.1) then
            gk(i,nch) = rk
            if (rk.eq.0.0) rk = 1.0
            if (e.ge.0.0) then
C     The T(kf,ki) matrix has been divided by KF and KI
               wk(kp) = - sum - cmplx(0.0, posfac * pi / rk)
C  Have multiplied the positronium-atom V-matrix elements by sqrt(2)
               wk(kp) = - sum - cmplx(0.0, pi / rk)
c$$$                        wk(kp) = - e * sum - cmplx(0.0, pi * rk)
            else
               wk(kp) = (0.0,0.0)
            end if 
            if (analytic) wk(kp)=cmplx(0.0,imag(wk(kp))) !needed for PHOTO (0.0,0.0)
            if (lprint) print*,'On-shell weight:',wk(kp)
            if (abs(real(wk(kp))).gt.1e-2) then
               error stop 'Real part of on shell weight must be < 1e-2'
            endif 
         else
            gk(i,nch) = gridk(i-1)
            ek = gridk(i-1) * gridk(i-1)
C     The T(kf,ki) matrix has been divided by KF and KI
            if (abs(rk-gridk(i-1)).gt.1e-5) then
               wk(kp) = 2.0 * weightk(i-1)/(e - ek) ! * posfac
            else
               wk(kp) = 2.0 * weightk(i-1)/(4.0*e)! * posfac
               print*,'Dropped GF for rk/k:',rk/gridk(i-1)
               print*,'test:',(1.0/(e-gridk(i-2)**2)+
     >            1.0/(e-gridk(i)**2))/(0.5/e)
            endif
C  Have multiplied the positronium-atom V-matrix elements by sqrt(2)

C  Added by Alex
            if (analytic) then
               if (nbox.EQ.0) then
                  wk(kp) = weightk(i-1)
               else
                  wk(kp) = 1.0
               endif
            else
c$$$               wk(kp) = 2.0 * weightk(i-1)/(e - ek) ! was overwriting the posfac above
c$$$               print*,kp,wk(kp)
            endif
C  End added by Alex

c$$$            if (ek.eq.0.0) wk(kp) = 1.0
c$$$                     wk(kp) = ek * 2.0 * weightk(i-1)/(e - ek)

c$$$            wk(kp) = wk(kp) * 2.0
            
         end if
c$$$         if (nold.eq.0.and.ntstop.eq.2.and.na.gt.nnbtop) then
c$$$            if (i.eq.1) print*,'Multiplying by CINT:', cint(na-nnbtop)
c$$$            wk(kp) = wk(kp) * cint(na-nnbtop)
c$$$         endif 
      end do !nqk loop

      if (analytic) then
        nqk = -nqm - nbnd(lset) !nbnd get's added below
         sum = 0.0
      endif 

      if (lprint) print*
      nchtop = nch
      if (nqm.gt.0) nqk = min(nqk,nqm-1)
      nqk = nqk + nbnd(lset) ! bound state wk are set in MAKECHIL
      if (nqk+1.gt.kmax) then
         print*,'Channel, NQK+1, KMAX:',nch,nqk+1,kmax
         stop 'Increase KMAX'
      endif
c$$$      if (lprint) print*,'NCH,NKQ,NPK(NCH):',nch,nqk,npk(nch)      
      npk(nchtop+1) = npk(nchtop) + nqk + 1

      ldumm = -1
      if (la.le.ldumm.and.na.le.ndumm) nchopt = nchtop
      nch = nch + 1
      call getchinfo (nch, ntmp, lg, psi, maxpsi, ea, la, na, li)
      if (nch.gt.nchan) then
         print*,'Recompile with NCHAN >=',nch
         stop 'Recompile with bigger NCHAN'
      end if 
      if (nch.ne.0) go to 10
      
      if (lprint) print*,'NBAD:', nbad
      return
      end
c-----------------------------------------------------------------
c***************** Wigner coefficients ***************************
c-----------------------------------------------------------------
      function CGC0(WJ1,WJ2,WJ)
c     CGC0=CGC(WJ1,0.,WJ2,0.,WJ,0.)
      CGC0=((-1)**(NINT(WJ1-WJ2)))*(sqrt(2.*WJ+1.))*
     >   COF3J(WJ1,WJ2,WJ,0.,0.,0.)
      return
      end
c****************************************************************
      real function CGC(WJ1,WM1,WJ2,WM2,WJ,WM)
      CGC=((-1)**(NINT(WJ1-WJ2+WM)))*(sqrt(2.*WJ+1.))*
     >   COF3J(WJ1,WJ2,WJ,WM1,WM2,-WM)
c$$$  print*,'WM1,WM2,-WM, sum:',WM1,WM2,-WM,wm1+wm2-wm
      return
      end
c****************************************************************
C     FILE WIGNER.FOR
C
C     ================================================================
C
      FUNCTION COF12J ( X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4)
      REAL     COF12J, COF6J, X1, X2, X3, X4, Y1, Y2, Y3, Y4,
     >           Z1, Z2, Z3, Z4, A1, A2, A3, A4, B1, B2, B3, B4,
     >           G, G1, G2, D, PHASE, SUM, FUN
      INTEGER*4  NG1, NG2, NG3, I
C
C      *****************************************************************
C
C      THIS FUNCTION COMPUTES THE 12J SYMBOL IN TERMS OF 6J SYMBOL
C      SEE R. J. ORD-SMITH PHYS. REV. 94, 1227 (1954)
C
C     (X1  X2  X3  X4  )     S      -L       (X1 Z1 L ) (Z2 X2 L )
C     (  Y1  Y2  Y3  Y4) =(-) SUM(-)  (2L+1) (X4 Z4 Y4) (X1 Z1 Y1)
C     (Z1  Z2  Z3  Z4  )       L
C
C                                          * (X3 Z3 L ) (X4 Z4 L )
C                                            (Z2 X2 Y2) (Z3 X3 Y3)
C
C
C      *****************************************************************
C
      A1 = ABS( X1 - Z1 )
      A2 = ABS( X2 - Z2 )
      A3 = ABS( X3 - Z3 )
      A4 = ABS( X4 - Z4 )
      B1 = X1 + Z1
      B2 = X2 + Z2
      B3 = X3 + Z3
      B4 = X4 + Z4
      G1 = MAX(A1,A2,A3,A4) + 1.0D0
      G2 = MIN(B1,B2,B3,B4) + 1.0D0
      NG1 = G1 + 0.1
      NG2 = G2 + 0.1
      D = G1 - NG1
      PHASE = B1 + B2 + B3 + B4 + Y1 + Y2 + Y3 + Y4
      COF12J = 1.0D0
      SUM = 0.0D0
      IF ( NG2 .LT. NG1 )GO TO 2
  3   DO 1 I=NG1,NG2
      G = I + D - 1.0D0
      NG3 = PHASE + 3.0*G + 0.1
      FUN = ((-1.D0)**NG3)*(2.0*G+1.D0)*COF6J(X1,Z1,G,X4,Z4,Y4)*
     1       COF6J(Z2,X2,G,X1,Z1,Y1)*COF6J(X3,Z3,G,Z2,X2,Y2)*
     2       COF6J(X4,Z4,G,Z3,X3,Y3)
  1   SUM = SUM + FUN
  2   CONTINUE
      COF12J = COF12J * SUM
      RETURN
      END
C
C     ============================================================
C
      FUNCTION COF3J(A1,A2,A3,A4,A5,A6)
      REAL      A1, A2, A3, A4, A5, A6, A7, A8, A9,
     >            COF3J, CNJ
      real*8 sym3
      INTEGER*4   IC
C
      A7 = 0.0D0
      A8 = 0.0D0
      A9 = 0.0D0
      IC = 1
      call threej(nint(2*a1),nint(2*a2),nint(2*a3),
     >   nint(2*a4),nint(2*a5),nint(2*a6),sym3)
      cof3j = sym3
      if (a4.ne.0.0.or.a5.ne.0.0) then
         COF3Jold = CNJ(IC,A1,A2,A3,A4,A5,A6,A7,A8,A9)
         if (abs(cof3j-cof3jold).gt.1e-4) then
            print*,'Checking COF3J and THREEJ:',cof3j,cof3jold
         endif
      endif 
      RETURN
      END
C
C     =============================================================
C
      FUNCTION COF6J(A1,A2,A3,A4,A5,A6)
      REAL     A1, A2, A3, A4, A5, A6, A7, A8, A9,
     >           COF6J, CNJ
      INTEGER*4  IC
C      
      A7 = 0.0D0
      A8 = 0.0D0
      A9 = 0.0D0
      IC = 2
      COF6J = CNJ(IC,A1,A2,A3,A4,A5,A6,A7,A8,A9)
      RETURN
      END
C
C     =============================================================
C
      FUNCTION COF9J(A1,A2,A3,A4,A5,A6,A7,A8,A9)
      REAL     A1, A2, A3, A4, A5, A6, A7, A8, A9,
     >           COF9J, CNJ
      INTEGER*4  IC
C     
      IC = 3
      COF9J = CNJ(IC,A1,A2,A3,A4,A5,A6,A7,A8,A9)
      RETURN
      END
C
C     =============================================================
C
      FUNCTION CNJ(NJ,A1,A2,A3,A4,A5,A6,A7,A8,A9)
C
      REAL      CNJ, A1, A2, A3, A4, A5, A6, A7, A8, A9,
     >            AA, BB, CC, DD, EE, FF, GG, HH, PP,
     >            X, Y, Z, XX, YY, ZZ, XXX, YYY, ZZZ,
     >            FG01BD, FG01CD, FG01DD, FG02BD, FG02CD,
     >            FG03BD, CLEBSH, RACAH, ANINE, A, B, C,
     >            H, AY, RA, RB, UUU, VVV, WWW, Q, E, D,
     >            F, FACTOR
C
      INTEGER*4   NJ, K, K1, K2, K3, K4, K5, K6, KA, KB, KC,
     >            KK1, KK2, KK3, KK4, KK5, KK6, KK7, KK8, KK9,
     >            KUP, M1, M2, M3, M4, M5, M6, M7, M8, M9, M10,
     >            MM1, MM2, MM3, MM4, MM5, MA, MB, MC, MD, I,
     >            IJ, IX, IY, IZ, IYY, IAY, IJPAR, IERR, IERCT,
     >            N4, N5, N6, N5PAR, J, JS, JSPAR, JQ, 
     >            KEY, KEYTRI, KEYRAC

      INTEGER*4,SAVE::   JJJ
C
C   -----------------------------------------------------------------
C
C  2. ARGUMENT LIST
C
C     ALL ARGUMENTS MUST BE REAL INTEGERS OR HALF-ODD-INTEGERS (REAL
C  IN THE SINGLE-LENGTH ROUTINES, DOUBLE PRECISION IN THE DOUBLE-LENGTH ONES).
C  THE NEAREST EXACT VALUE WILL BE ASSUMED, BUT ERRORS MAY OCCUR IF THE
C  VALUES GIVEN IN THE CALLING ROUTINE ARE +OR-0.05 OR MORE IN ERROR.
C
C
C  3. RESTRICTIONS
C
C     THE SUM OF THE THREE ANGULAR MOMENTA APPEARING IN ANY
C  "TRIANGULAR CONDITION" MUST NOT EXCEED 100.0.(MMAX) THIS LIMIT CAN BE RAISED
C  IF NECESSARY BY RECOMPILING WITH LARGER DIMENSIONS FOR THE ARRAYS H AND
C  J AND A CORRESPONDINGLY LARGER UPPER LIMIT ON THE INDEX OF THE FIRST DO
C  LOOP.
C
C     THE FOLLOWING "GEOMETRICAL" CONDITIONS ARE TESTED BY THE
C  PROGRAMS AND THE CORRECT VALUE OF ZERO RETURNED IN CASE OF VIOLATION:
C
C     (A) ALL TRIANGULAR CONDITIONS SATISFIED
C     (B) ALL ANGULAR MOMENTA NON-NEGATIVE
C     (C) IN COF3J AND FGO1B, (A+X), (B+Y) AND (C+Z) ARE INTEGRAL
C     (D) IN COF3J, X+Y+Z = 0.0; IN FGO1B, X+Y = Z
C     (E) IN FGO1C AND D, A,B AND C ALL INEGRAL WITH A+B+C EVEN.
C
C     SINCE A VIOLATION OF THESE CONDITIONS MAY BE THE RESULT OF AN
C  ERROR IN THE CALLING ROUTINE WHICH SETS UP THE ARGUMENTS, AN "ERROR
C  CHECK" IS PROVIDED.
C  A LABELLED COMMON AREA IS SET UP:
C
C     COMMON/FGERCM/IERR,IERCT
C
C  WHERE THE TWO MEMBERS ARE INTEGER*4. IERR WILL BE UNITY AT RETURN IF ANY
C  CONDITION WAS VIOLATED; ZERO IF NOT. THIS COUNT IS RESET AT EACH CALL
C  OF ANY OF THE ROUTINES IN THE PACKAGE. IERCT IS ZERO AT THE START OF THE
C  JOB AND IS INCREASED BY UNITY WHENEVER A CONDITION IS VIOLATED IN A CALL
C  OF ANY OF THE ROUTINES OF THE PACKAGE; THUS A SINGLE CHECK AT THE END OF
C  THE ENTIRE JOB CAN ENSURE THAT NO VIOLATIONS OCCURED DURING THE JOB.
C
C
C  4. GENERAL
C
C     THE 9 NAMES ARE DIFFERENT ENTRY POINTS TO A SINGLE ROUTINE
C  (CSECT NAME COF6J). THIS SAVES SPACE SINCE ALL USE THE SAME PRIVATE
C  TABLE OF FACTORIALS (1K BYTE). USERS WHO NEED ONLY ONE OF THESE
C  FUNCTIONS IN THEIR JOB AND WHO ARE PRESSED FOR SPACE MAY OF COURSE
C  RECOMPILE THE SOURCE ROUTINE AFTER REMOVING MOST OF THE CODE BELONGING
C  PRIVATELY TO THE OTHER ENTRIES (NOTE THAT C9J NEEDS MOST OF THE FGO2
C  CODE; AND THAT ALL ENTRIES NEED THE DO LOOP NEAR THE BEGINNING THAT SETS
C  UP THE FACTORIAL TABLES.
C
C     THE ROUTINES ARE SELF CONTAINED AND CAUSE NO OUTPUT; THEY SET
C  UP ONE NAMED COMMON AREA (DESCRIBED IN 3 ABOVE).
C
C     FGO3B IS COMPLETELY IDENTICAL TO FGO3A. IT IS INCLUDED FOR
C  COMPATIBILITY WITH EARLIER VERSIONS.
C
C
C
C##       FG01BD         08/01/74
C NAME FG01BD(R)                 CHECK



C
C     WIGNER AND RACAH COEFFICIENTS
C
C     FG01A - WIGNER 3-J SYMBOL
C     FG01B - CLEBSCH-GORDAN COEFFICIENT
C     FG01C & FG01D - SAME WITH ZERO MAGNETIC QUANTUM NUMBERS
C
C     FG02A - WIGNER 6-J SYMBOL
C     FG02B - RACAH COEFFICIENT
C     FG02C - U-FUNCTION (JAHN)
C
C     FG03A & FG03B - WIGNER 9-J SYMBOL
C
C
      parameter (mmax=1010)
      COMMON / CNJSAVE / H(mmax), J(mmax)
c$omp threadprivate(/CNJSAVE/)
      DIMENSION AY(4),IAY(4)
      COMMON/ FGERCM /IERR,IERCT
      DATA JJJ/0/
c$omp threadprivate(JJJ)
      INTPTF(Q) = nint(Q + Q)
      IPARF(I) = 4*(I/4) - I + 1
C      GO TO 77

C  COMPUTED GOTO 'SIMULATES' MULTIPLE ENTRY POINTS.

      GO TO (9993,9996,9999),NJ

C

9993  CONTINUE
      A = A1
      B = A2
      C = A3
      XX = A4
      YY = A5
      ZZ = A6

      ZZ=-ZZ
      KEY=1
      IERR=0
      GO TO 1
C
C
C
9996  CONTINUE
      UUU=A1
      VVV=A2
      WWW=A3
      XXX=A4
      YYY=A5
      ZZZ=A6


      KEY=11
      IERR=0
      GO TO 100
   77 KEY=2
      IERR=0
    1 K1=INTPTF(A)
      K2=INTPTF(B)
      K3=INTPTF(C)
      IF(KEY.GE.3) GO TO 100
      K4=INTPTF(XX)
      K5=INTPTF(YY)
      K6=INTPTF(ZZ)
C
  100 IF(JJJ.NE.0) GO TO 500
      JJJ=1
      IERCT=0
      H(1)=1.0D0
      J(1)=0
      X=0.D0
      DO 400 I=2,mmax
      X=X+1.0D0
      H(I)=H(I-1)*X
      J(I)=J(I-1)
  200 IF(H(I).LT.10.0D0) GO TO 400
      H(I)=0.01D0*H(I)
      J(I)=J(I)+2
      GO TO 200
  400 CONTINUE
C
  500 IF(KEY.LT.-5) GO TO 750
      IF(KEY.GE.3) GO TO 320
      IF((K4+K5-K6).NE.0) GO TO 710
      M1=K1+K2-K3
      M2=K2+K3-K1
      M3=K3+K1-K2
      M4=K1+K4
      M5=K1-K4
      M6=K2+K5
      M7=K2-K5
      M8=K3+K6
      M9=K3-K6
      M10=K1+K2+K3+2
C
      IF(M1.LT.0) GO TO 710
      IF(M2.LT.0) GO TO 710
      IF(M3.LT.0) GO TO 710
      IF(M4.LT.0) GO TO 710
      IF(M5.LT.0) GO TO 710
      IF(M6.LT.0) GO TO 710
      IF(M7.LT.0) GO TO 710
      IF(M8.LT.0) GO TO 710
      IF(M9.LT.0) GO TO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
      IF((M10-(M10/2)-(M10/2)).NE.0) GO TO 710
C
      Y=K3+1
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      M4=M4/2+1
      M5=M5/2+1
      M6=M6/2+1
      M7=M7/2+1
      M8=M8/2+1
      M9=M9/2+1
      M10=M10/2+1
      if (max(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10).gt.mmax) then
         print*,'Recompile wigner.f with new mmax at least:',
     >      max(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
         stop 'Increse MMAX in CNJ (file wigner.f)'
      endif 
C
      Y  = SQRT(Y*H(M1)*H(M2)*H(M3)*H(M4)*H(M5)*
     >     H(M6)*H(M7)*H(M8)*H(M9)/H(M10))
      IY = (J(M1)+J(M2)+J(M3)+J(M4)+J(M5)+
     >     J(M6)+J(M7)+J(M8)+J(M9)-J(M10))/2
C
      N4=M1
      IF(N4.GT.M5)N4=M5
      IF(N4.GT.M6)N4=M6
      N4=N4-1
      M2=K2-K3-K4
      M3=K1+K5-K3
      N5=0
      IF(N5.LT.M2) N5=M2
      IF(N5.LT.M3) N5=M3
      N5PAR=IPARF(N5)
      N5=N5/2
      Z=0.0D0
      GO TO 610
C
  700 MM1=M1-N5
      MM2=M5-N5
      MM3=M6-N5
      MM4=N5-(M2/2)+1
      MM5=N5-(M3/2)+1
C
      X  = 1.D0/(H(MM1)*H(MM2)*H(MM3)*H(MM4)*H(MM5)*H(N5+1))
      IX = -J(MM1)-J(MM2)-J(MM3)-J(MM4)-J(MM5)-J(N5+1)
C
  800 IF(IX+IY)900,210,110
  900 X=0.1D0*X
      IX=IX+1
      GO TO 800
  110 X=10.0D0*X
      IX=IX-1
      GO TO 800
C
  210 IF(N5PAR.LT.0) X=-X
      Z=Z+X
  510 N5PAR=-N5PAR
      N5=N5+1
C
  610 IF(N5-N4)700,700,810
C
 710  CLEBSH=0.0D0
      IERR=1
      IERCT=IERCT+1
      GO TO 910
C
 810  CLEBSH=Z*Y
  910 GO TO(120,220),KEY
C
  220 FG01BD=CLEBSH
      RETURN
C
  120 JS=K1-K2+K6
      IF(JS.LT.0) JS=-JS
      JSPAR=IPARF(JS)
      CNJ=JSPAR*CLEBSH/SQRT(K3+1.0D0)
      ZZ=-ZZ
      RETURN
C
C
C
  320 IF(KEY.GE.10) GO TO 130
      KEY=KEY-2
      IF((K1-(K1/2)-(K1/2)).NE.0) GO TO 420
      IF((K2-(K2/2)-(K2/2)).NE.0) GO TO 420
      IF((K3-(K3/2)-(K3/2)).NE.0) GO TO 420
      IJ=K1+K2+K3
      IJPAR=IPARF(IJ)
      IF(IJPAR.LE.0) GO TO 420
      M1=IJ-K1-K1
      M2=IJ-K2-K2
      M3=IJ-K3-K3
      M4=IJ+2
      IF(M1.LT.0) GO TO 420
      IF(M2.LT.0) GO TO 420
      IF(M3.LT.0) GO TO 420
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      M4=IJ/2+2
      Y  = SQRT(H(M1)*H(M2)*H(M3)/H(M4))
      IY = (J(M1)+J(M2)+J(M3)-J(M4))/2
      IJ=IJ/2
      IJPAR=IPARF(IJ)
      IJ=IJ/2+1
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      Z=H(IJ)/(H(M1)*H(M2)*H(M3))
      IZ=J(IJ)-J(M1)-J(M2)-J(M3)
      IZ=IZ+IY
      CLEBSH=IJPAR*Y*Z*10.0D0**IZ
      GO TO(620,720),KEY
C
  620 FG01CD=CLEBSH
      RETURN
C
  720 JQ=K2-K1
      IF(JQ.LT.0) JQ=-JQ
      IJPAR=IPARF(JQ)
      FG01DD=CLEBSH*IJPAR*SQRT(K3+1.0D0)
      RETURN
C
  420 CLEBSH=0.0D0
      IERR=1
      IERCT=IERCT+1
      GO TO(620,720),KEY
C
  130 IF(KEY.EQ.11) GO TO 450
      IF(KEY.GT.19) GO TO 750
  550 stop'should not be here as D, E, and F are not defined (wigner.f)'
      K1=INTPTF(A)
      K2=INTPTF(B)
      K3=INTPTF(E)
      K4=INTPTF(D)
      K5=INTPTF(C)
      K6=INTPTF(F)
C
C
  750 KA=K1
      KB=K2
      KC=K3
      KEYTRI=1
      GO TO 630
C
  230 KA=K4
      KB=K5
      KEYTRI=2
      GO TO 630
C
  330 KB=K2
      KC=K6
      KEYTRI=3
      GO TO 630
C
  430 KA=K1
      KB=K5
      KEYTRI=4
      GO TO 630
C
  530 YY=AY(1)*AY(2)*AY(3)*AY(4)
      IYY=IAY(1)+IAY(2)+IAY(3)+IAY(4)
      M1=(K1+K2+K4+K5)/2+2
      M2=(K1+K2-K3)/2+1
      M3=(K4+K5-K3)/2+1
      M4=(K1+K5-K6)/2+1
      M5=(K2+K4-K6)/2+1
      M6=K1+K4-K3-K6
      M7=K2+K5-K3-K6
C
      N4=M1
      IF(N4.GT.M2) N4=M2
      IF(N4.GT.M3) N4=M3
      IF(N4.GT.M4) N4=M4
      IF(N4.GT.M5) N4=M5
      N4=N4-1
      N5=0
      IF(N5.LT.M6) N5=M6
      IF(N5.LT.M7) N5=M7
      N5PAR=IPARF(N5)
      N5=N5/2
      M6=M6/2-1
      M7=M7/2-1
      Z=0.0D0
      GO TO 730
C
  140 X  = H(M1-N5)/(H(N5+1)*H(M2-N5)*H(M3-N5)*H(M4-N5)
     >     *H(M5-N5)*H(N5-M6)*H(N5-M7))
      IX = J(M1-N5)-J(N5+1)-J(M2-N5)-J(M3-N5)-J(M4-N5)
     >    -J(M5-N5)-J(N5-M6)-J(N5-M7)
  240 IF(IX+IYY)340,440,540
  340 X=0.1D0*X
      IX=IX+1
      GO TO 240
  540 X=10.0D0*X
      IX=IX-1
      GO TO 240
  440 IF(N5PAR.LT.0) X=-X
      Z=Z+X
      N5PAR=-N5PAR
      N5=N5+1
C
  730 IF(N5.LE.N4) GO TO 140
C
      RACAH=Z*YY
  840 IF(KEY.LT.-5) GO TO 160
      KEY=KEY-10
      GO TO(150,250,350),KEY
C
  830 RACAH=0.0D0
      IERR=1
      IERCT=IERCT+1
      GO TO 840
C
  150 IJPAR=IPARF(K1+K2+K4+K5)
      IF(IJPAR.LT.0) RACAH=-RACAH
      CNJ=RACAH
      RETURN
C
  250 FG02BD=RACAH
      RETURN
C
  350 FACTOR = SQRT((K3+1.0D0)*(K6+1))
      FG02CD = FACTOR*RACAH
      RETURN
  450 K1 = INTPTF(UUU)
      K2 = INTPTF(VVV)
      K3=INTPTF(WWW)
      K4=INTPTF(XXX)
      K5=INTPTF(YYY)
      K6=INTPTF(ZZZ)
      GO TO 750
C
C
C     TRIANGLE FUNCTION
C
C
  630 MA=KA+KB-KC
      MB=KA-KB+KC
      MC=-KA+KB+KC
      MD=KA+KB+KC+2
      IF(MA.LT.0) GO TO 830
      IF(MB.LT.0) GO TO 830
      IF(MC.LT.0) GO TO 830
      IF((MD-(MD/2)-(MD/2)).NE.0) GO TO 830
      MA=MA/2+1
      MB=MB/2+1
      MC=MC/2+1
      MD=MD/2+1
      if (max(ma,mb,mc,md).gt.mmax) stop
     >   'Increse MMAX in CNJ (file wigner.f)'
      AY(KEYTRI) = SQRT(H(MA)*H(MB)*H(MC)/H(MD))
      IAY(KEYTRI) = (J(MA)+J(MB)+J(MC)-J(MD))/2
      GO TO(230,330,430,530),KEYTRI
C
C
C
C
C

9999  CONTINUE
      AA=A1
      BB=A2
      CC=A3
      DD=A4
      EE=A5
      FF=A6
      GG=A7
      HH=A8
      PP=A9


C      ENTRY FG03BD(AA,BB,CC,DD,EE,FF,GG,HH,PP)
C
      KEY=-10
      IERR=0
C
      KK1=INTPTF(AA)
      KK2=INTPTF(BB)
      KK3=INTPTF(CC)
      KK4=INTPTF(DD)
      KK5=INTPTF(EE)
      KK6=INTPTF(FF)
      KK7=INTPTF(GG)
      KK8=INTPTF(HH)
      KK9=INTPTF(PP)
C
      KUP=KK1+KK9
      M1=KK4+KK8
      M2=KK2+KK6
      IF(KUP.GT.M1) KUP=M1
      IF(KUP.GT.M2) KUP=M2
C
      K=KK1-KK9
      IF(K.LT.0) K=-K
      M1=KK4-KK8
      IF(M1.LT.0) M1=-M1
      M2=KK2-KK6
      IF(M2.LT.0) M2=-M2
      IF(K.LT.M1) K=M1
      IF(K.LT.M2) K=M2
C
      ANINE=0.0D0
C
  660 IF(K.GT.KUP) GO TO 260
      K1=KK1
      K2=KK4
      K3=KK7
      K4=KK8
      K5=KK9
      K6=K
      KEYRAC=1
      GO TO 100
C
  160 GO TO(360,460,560),KEYRAC
C
  360 RA=RACAH
      K1=KK2
      K2=KK8
      K3=KK5
      K4=KK4
      K5=KK6
      KEYRAC=2
      GO TO 750
C
  460 RB=RACAH
      K1=KK9
      K2=KK6
      K3=KK3
      K4=KK2
      K5=KK1
      KEYRAC=3
      GO TO 750
C
  560 ANINE=ANINE+RA*RB*RACAH*(K+1)
      K=K+2
      GO TO 660
C
  260 CNJ=ANINE
      FG03BD=ANINE
      RETURN
      END
C
C     ===============================================================
C
      FUNCTION BICO (A,B)
C
      REAL     BICO, A, B, FACTOR, X, AI, BJ
      INTEGER*4  IFIRST, I, J
C
C     THIS PROGRAM CAN ONLY HANDLE INTEGERS.
C
      COMMON / LOGFAC / FACTOR(230)
      DATA IFIRST / 0 /
      IF ( IFIRST .EQ. 0 ) CALL FACLOG
      IFIRST = 1
      X=A-B
      IF (X) 201,202,203
  201 BICO=0.D0
      GO TO 207
  202 BICO=1.D0
      GO TO 207
  203 IF (B) 201,202,204
 204  CONTINUE
      I=A+0.1
      J=B+0.1
      AI=I
      BJ=J
      IF(ABS(AI-A).GT.0.001) GOTO500
      IF(ABS(BJ-B).GT.0.001) GOTO500
      X=FACTOR(I)-FACTOR(J)-FACTOR(I-J)
      BICO= EXP(X)
      GO TO 207
  500 WRITE(9,501) A,B
  501 FORMAT(10X,'  A=',F10.3,'B=',F10.3,' NON INTEGER VALUES FOR BICO')
       BICO=0.0D0
 207   CONTINUE
       RETURN
      END
C
C     =============================================================
C
      SUBROUTINE FACLOG
C
      INTEGER*4         I 
      REAL            FACTOR, A
      COMMON / LOGFAC / FACTOR(230)
C
      FACTOR(1)=0.0D0
      DO 10 I=2,230
        A=I
        FACTOR(I)=FACTOR(I-1)+LOG(A)
10    CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CLEB(IA,IB,IC,ID,IE,IF,RAC)
      COMMON /flogs/ FACLOG(1000)
      double precision faclog
      
      RAC=0.0

      IF( (ID.EQ.0) .AND. (IE.EQ.0) .AND. (IF.EQ.0) ) THEN
      LTL= (IA+IB+IC)/2
      IF( (LTL -2*(LTL/2) ).NE.0) RETURN
      ENDIF

      IF ((ID+IE-IF).NE.0) RETURN
      K1=IA+IB+IC
      IF ((K1-2*(K1/2)).NE.0) RETURN
      K1=IA+IB-IC
      K2=IC-IABS(IB-IA)
      K3=MIN0(K1,K2)
      IF (K3.LT.0) RETURN
      IF (((-1)**(IB+IE)).LT.0) RETURN
      IF (((-1)**(IC+IF)).LT.0) RETURN
      IF ((IA-IABS(ID)).LT.0) RETURN
      IF ((IB-IABS(IE)).LT.0) RETURN
      IF ((IC-IABS(IF)).LT.0) RETURN
      IF (IA.LT.0) RETURN
      IF (IB.LT.0) RETURN
      IF (IC.LT.0) RETURN
      RAC=1.0
      IF (IA.EQ.0) RETURN
      IF (IB.EQ.0) RETURN
      RAC=0.0
      IF (IC.GT.0) GO TO 130

      FB=IB+1
      RAC=((-1.0)**((IA-ID)/2))/SQRT (FB)
      RETURN

  130 FC2=IC+1
      IABCP=(IA+IB+IC)/2+1
      IABC=IABCP-IC
      ICAB=IABCP-IB
      IBCA=IABCP-IA
      IAPD=(IA+ID)/2+1
      IAMD=IAPD-ID
      IBPE=(IB+IE)/2+1
      IBME=IBPE-IE
      ICPF=(IC+IF)/2+1
      ICMF=ICPF-IF
      SQFCLG=0.5*(ALOG(FC2)-FACLOG(IABCP+1)
     *      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)
     *      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)
     *      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))
      NZMIC2 = (IB-IC-ID)/2
      NZMIC3 = (IA-IC+IE)/2
      NZMI= MAX0 (0,NZMIC2,NZMIC3)+1
      NZMX= MIN0 (IABC,IAMD,IBPE)
      S1=(-1.0)**(NZMI-1)
      DO 140 NZ=NZMI,NZMX
      NZM1=NZ-1
      NZT1=IABC-NZM1
      NZT2=IAMD-NZM1
      NZT3=IBPE-NZM1
      NZT4=NZ-NZMIC2
      NZT5=NZ-NZMIC3
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)
     *           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)
      SSTERM=S1*EXP (TERMLG)
      RAC=RAC+SSTERM
  140 S1=-S1
      RETURN
      END

      

      SUBROUTINE RAC7(IA,IB,IC,ID,IE,IF,RAC)
      COMMON /flogs/ FACLOG(1000)
      double precision faclog
      
      DIMENSION LT(6)

      K1=IA+IB-IE
      K3=IC+ID-IE
      K5=IA+IC-IF
      K7=IB+ID-IF
      K2=IE-IABS( IA-IB)
      K4=IE-IABS( IC-ID)
      K6=IF-IABS( IA-IC)
      K8=IF-IABS( IB-ID)
      K9= MIN0( K1,K2,K3,K4,K5,K6,K7,K8)
      RAC=0.0
      IF (K9.LT.0) RETURN
      K2=K1-2*(K1/2)
      K4=K3-2*(K3/2)
      K6=K5-2*(K5/2)
      K8=K7-2*(K7/2)
      IF (MAX0(K2,K4,K6,K8).NE.0) RETURN
      LTMIN= MIN0(IA,IB,IC,ID,IE,IF)
      IF (LTMIN.LT.0) RETURN
      IF (LTMIN.GT.0) GO TO 90

      LT(1)=IA
      LT(2)=IB
      LT(3)=IC
      LT(4)=ID
      LT(5)=IE
      LT(6)=IF
      LTMIN=LT(1)
      KMIN=1
      DO 50 N=2,6
      IF ((LT(N)-LTMIN).GE.0) GO TO 50
      LTMIN=LT(N)
      KMIN=N
   50 CONTINUE
      S1=1.0
      F1=IE
      F2=IF
      GO TO (80,80,80,80,60,70),KMIN
   60 F1=IA
      F2=IC
      S1=(-1.0)**(K5/2)
      GO TO 80
   70 F1=IA
      F2=IB
      S1=(-1.0)**(K1/2)
   80 RAC=S1/SQRT( (F1+1.)*(F2+1.))
      RETURN

   90 IABEP=(IA+IB+IE)/2+1
      ICDEP=(IC+ID+IE)/2+1
      IACFP=(IA+IC+IF)/2+1
      IBDFP=(IB+ID+IF)/2+1
      IABE=IABEP-IE
      IEAB=IABEP-IB
      IBEA=IABEP-IA
      ICDE=ICDEP-IE
      IECD=ICDEP-ID
      IDEC=ICDEP-IC
      IACF=IACFP-IF
      IFAC=IACFP-IC
      ICFA=IACFP-IA
      IBDF=IBDFP-IF
      IFBD=IBDFP-ID
      IDFB=IBDFP-IB
      NZMAX= MIN0( IABE,ICDE,IACF,IBDF)
      IABCD1=(IA+IB+IC+ID+4)/2
      IEFMAD=(IE+IF-IA-ID)/2
      IEFMBC=(IE+IF-IB-IC)/2
      NZMI1=-IEFMAD
      NZMI2=-IEFMBC
      NZMIN= MAX0( 0,NZMI1,NZMI2)+1
      SQLOG=0.5*(FACLOG(IABE)+FACLOG(IEAB)+FACLOG(IBEA)+FACLOG(ICDE)
     *          +FACLOG(IECD)+FACLOG(IDEC)+FACLOG(IACF)+FACLOG(IFAC)
     *          +FACLOG(ICFA)+FACLOG(IBDF)+FACLOG(IFBD)+FACLOG(IDFB)
     * -FACLOG(IABEP+1)-FACLOG(ICDEP+1)-FACLOG(IACFP+1)-FACLOG(IBDFP+1))
      DO 100 NZ=NZMIN,NZMAX
      NZM1=NZ-1
      K1=IABCD1-NZM1
      K2=IABE-NZM1
      K3=ICDE-NZM1
      K4=IACF-NZM1
      K5=IBDF-NZM1
      K6=NZ
      K7=IEFMAD+NZ
      K8=IEFMBC+NZ
      SSLOG=SQLOG+FACLOG(K1)-FACLOG(K2)-FACLOG(K3)-FACLOG(K4)
     *           -FACLOG(K5)-FACLOG(K6)-FACLOG(K7)-FACLOG(K8)
      SSTERM=((-1.0)**NZM1)*EXP( SSLOG)
      RAC=RAC+SSTERM
  100 CONTINUE
      RETURN
      END

      function polpot(r,corep,r0)
      include 'par.f'
      integer nin(0:1)
      real*8 rin(0:maxr,0:1),vin(0:maxr,0:1),rout,vout(0:1),dfac
      logical exists(0:1)
      character ch
      ch(i)=char(i+ichar('0'))
      data exists/.false.,.false./
      save exists,rin,vin,nin
      common /TMP_AB/ A,B

      if (corep.eq.0.0) then ! The following ensures the HF routines are called without polpot when using Bob McEachran's potentials
         polpot = 0.0
         return
      endif 

c$$$      polpot = - corep / (1.0 + (r/r0)**2)**2

      if (r.lt.0.01) then
         polpot = - corep * r*r / 2d0 / r0**6 
      else
         polpot = - corep / 2.0 / r**4 * (1d0 - exp(-(dble(r/r0))**6))
      end if
      
c$$$     >   - corepq / 2.0 / r**6 * (1d0 - exp(-(dble(r/r0))**10))

C The following reads Bob McEachran's ab initio pol pots.
      if (.not.exists(0)) then
         do iq = 0, 1
            inquire(file='polpotr'//ch(iq),exist=exists(iq))
            if (exists(iq)) then
               print*,'Reading polpotr'//ch(iq)
               open(42,file='polpotr'//ch(iq))
               read(42,*)
               read(42,*)
               read(42,*)
               rin(0,iq) = 0d0
               vin(0,iq) = 0d0
               i = 1
 10            read(42,*,end=20) rin(i,iq),vin(i,iq)
               if (rin(i,iq).lt.0.1) then
                  dfac = rin(i,iq)**(5+iq*2)
               else
                  dfac = 1d0-exp(-rin(i,iq)**(5+iq*2))
               endif 
               vin(i,iq) = vin(i,iq) / dfac
               i = i + 1
               go to 10
 20            close(42)
               nin(iq) = i
            endif 
         enddo
      endif 
      do iq = 0, 1
         if (exists(iq)) then
            rout = r
            nout = 1
            if (rout.gt.rin(nin(iq)-1,iq)) then
               vout(iq) = vin(nin(iq)-1,iq) / rout**(4+iq*2)
            else
               call intrpl(nin(iq),rin(0,iq),vin(0,iq),nout,rout,
     >            vout(iq))
               if (rout.lt.0.1) then
                  dfac = rout**(5+iq*2)
               else
                  dfac = 1d0-exp(-rout**(5+iq*2))
               endif 
               vout(iq) = vout(iq) * dfac / rout**(4+iq*2)
            endif
         endif 
      enddo
      if (exists(0)) polpot = vout(0)
      if (exists(1)) polpot = vout(0) - vout(1)
         
      end
      

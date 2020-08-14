      subroutine pwrite(nent,instate,nopen,energy,nznuci,
     >   zasymi,ry,noprint,ovlp,
     >   ovlpn,phasen,phaseq,jstart,jstop,projectile,target,nsmax,
     >   nchimax,nchanmax,npar,hlike,NTunit,vdcore,minvdc,maxvdc,ne2e,
     >   slowery,iborn,bornICS,tfile,nnbtop,ovlpnl,ovlpnn)
      implicit real (a-h,o-z)
      include 'par.f'
      parameter (ntype=2*(lamax+1),npoints=ncmax+10)
      character projectile*8,target*6,chunit(3)*8,chun*8,tfile*(*)
      dimension lea(nchan),onshellk(nchan),onshold(nchan),temp(maxr),
     >   psii(maxr),psif(maxr),ff(0:200,0:lamax),chi(maxr),ovlpn(knm),
     >   sweight(0:200,10),bornICS(knm,knm),ovlp(ncmax,0:lamax),
     >   ucentr(maxr),units(3),BornPCS(nchan,nchan),esum(0:1),
     >   partcs(nchan,nchan,0:1),sigtop(nchan,0:1),rcond(0:1),
     >   vdcore(maxr,0:lamax),test(6),slowery(ncmax),
     >   npt(npoints,-lamax:lamax),enf(npoints,ntype,-lamax:lamax),
     >   le2e(ntype),tin(npoints,ntype,-lamax:lamax,4),spinw(0:1),
     >   sdcs(ntype,npoints+100), finalcont(3),sigint(knm,ntype,0:1),
     >   sigintj(ntype,-lamax:lamax), tics(0:1),tcs(0:1),tnbcs(0:1),
     >   ovlpnn(knm,nnmax),ovlpnl(ncmax,0:lamax,nnmax),ntfun(0:lamax,3),
     >   nse2e(ntype),ovlpnold(knm),laold(knm),
     >   psdcs(ntype,npoints+100,0:1),res2(5)
      complex ovlpnl,ovlpnn
C  XIN,...,YOUT are used for (e,2e) and BORN interpolation
      real*8 rylm,xin(181),yin(181),xout(181*32),
     >   yout(181*32,4), wout(npoints + 100), wf(2*100), deta
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntyp,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      common /channels/ nfich(nchan)
      integer satom
      common /helium/ latom(KNM), satom(KNM), lpar(KNM), npp(KNM)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      logical noprint(nchan), exists, hlike, convergence, canstop,
     >   positron, posf, posi, Ps_start
      integer jlm(knm,0:1,0:lmax), non(knm,2*lamax+1,0:1,0:lmax),
     >   nonnew(knm,-lamax:lamax,0:1,0:lmax)
      integer nstep(10),nnset(3),nopen(0:1),nopenold(0:1),
     >   nchfi(knm,0:1), iwf(2*100), nttop(0:1), instate(100)
      parameter (lentr=8,lexit=lamax)
C ton(nchanmax,nchimax,0:nsmax,0:npar,0:jstop), 
      complex phase, sigc, !vdon(nchanmax,nchimax,0:nsmax,0:npar,0:jstop),
     >   phaseq(nchan),
     >   Bornamp(0:200,-lexit:lexit,-lentr:lentr),phasen(ncmax,0:lamax),
     >   phasediv(nchan),
     >   phasetim(ntype,npoints),te2ec(-latop:latop,ntype),fe2e(0:1),
     >   f0(2),f1(2),ce2e!te2e(0:nsmax,0:jstop,-latop:latop,ntype,ne2e*2)
      complex, allocatable :: te2e(:,:,:,:,:), ton(:,:,:,:,:),
     >   vdon(:,:,:,:,:)
      character chan(knm)*3,chanf(ntype)*3,csfile*30,sdcsfile*7,ch,
     >   cnode*3, ench*11, fmt*80
      complex*16 xx,eta1,zlmin,cfc(9),cgc(9),cfcp(9),cgcp(9),sig(9),ci,
     >   coulphase
      common /MPI_info/myid, ntasks, cnode, ench
      common /charchan/ chan
      common /chanen/ enchan(knm)
c$$$      ch(i)=char(i+ichar('0')) 
c$$$      pointer (ptrv, vdon)
c$$$      pointer (ptrt, ton)
c$$$      pointer (ptrte2e, te2e)
      data pi,chunit/3.1415927,' a0^2   ',' pi a0^2',' cm^2  '/
      data units/1.0,0.3183099,28.002e-18/
      if (ne2e*2.gt.npoints) stop 'INCREASE NPOINTS IN PWRITE'

      nunit = max(1,NTunit)
      if (nunit.gt.0) unit = units(nunit)

      le2e(:) = -1
      chanf(:) = 'XXX'
      do nn = 1, 10
         nstep(nn) = nn
         if (nn.eq.4) then
            nstep(nn) = 5
         elseif (nn.eq.5) then
            nstep(nn) = 6
         elseif (nn.eq.6) then
            nstep(nn) = 9
         elseif (nn.eq.7) then
            nstep(nn) = 10
         elseif (nn.eq.8) then
            nstep(nn) = 15
         elseif (nn.eq.9) then
            nstep(nn) = 18
         elseif (nn.eq.10) then
            nstep(nn) = 30
         endif             
         h = pi / 180.0 * nstep(nn)
         nprev = 1
         do n = 0, 200
            sweight(n,nn) = 0.0
         enddo 
         do n = 0, 180, nstep(nn)
            sweight(n,nn) = float(nprev) / 3.0 * h
            if (nprev.eq.4) then
               nprev = 2
            else
               nprev = 4
            endif 
         enddo
         sweight(0,nn) = 1.0/3.0 * h
         sweight(180,nn) = 1.0/3.0 * h
         sum = 0.0
         test(1) = 0.0
         do n = 0, 180, nstep(nn)
            th = n * pi / 180.0
            sum = sum + sin(th) * sweight(n,nn)
            test(1) = test(1) + th ** 2 * sweight(n,nn)
         enddo
         print*,'Nstep, Check of integration:', nstep(nn),
     >      sum / 2.0,test(1)*3.0/pi**3
      enddo
         
      npos = 0
      if (nze.eq.1) npos = 1
      zassymp = zasym
      xx = (1.0,0.0)
      ci = (0.0,1.0)
      nl = 1
      mode1 = 4
      kfn = 0

      do ip = 0, 1
         nttop(ip) = 0
         do nt = 1, ntype
            npt(nt,ip) = 0
         enddo
      enddo 

      Ps_start = .false.
      etot = energy / ry + enpsinb(ninc,linc)
      do ipar = 0, npar
         do n = 1, knm
            nchfi(n,ipar) = 0
         enddo 
         nch = 1
         nchpp = 0
         call getchinfo (nch, nchp, ipar, temp, maxpsi, ea, la, na, l)
         nopen(ipar) = 0
 10      continue 
         echan = (enpsinb(na,la)-enpsinb(ninc,linc)) * ry
         if (energy.ge.echan) then
c$$$         if ((energy.ge.echan.and.ea.lt.0.0).or.
c$$$     >      (energy.ge.echan.and.ne2e.ne.0)) then
c$$$            nchfi(nchp,ipar) = nch
C  Think about the role of NOPEN and NCH some more
c$$$  nopen = nch
C  The following line has the "if" statement as for IPAR = 1 and
C  unnatural parity states there is not 1:1 correspondence between
C  nch and nchp
            if (nchp.ne.nchpp) nopen(ipar) = nopen(ipar) + 1
            nchpp = nchp
            nchfi(nchp,ipar) = nopen(ipar)
            if (ipar.eq.0) then
c$$$               if (mod(nchfi(nchp,ipar),4).ne.0) then
               if (mod(nchp,4).ne.0) then
                  print '(i3,a4,"*",f10.3," eV,",$)',nch,chan(nch),ea*ry
               else
                  print '(i3,a4,"*",f10.3," eV")',nch, chan(nch),ea*ry
               endif
               posfac = 1.0
               if (positron(na,la,npos)) posfac = 2.0
               onshellk(nch) = sqrt(posfac*(energy - echan)/ry)
               lea(nchfi(nchp,ipar)) = 2 * la
C  The following change was done for the unnatural parity
               noprint(nchp) = .true.
            endif 
         else
            if (ipar.eq.0) then
               if (mod(nchp,4).ne.0) then
                  print '(i3,a4,f11.3," eV,",$)',nch,chan(nch),ea*ry
               else
                  print '(i3,a4,f11.3," eV")',nch, chan(nch),ea*ry
               endif
               onshellk(nch) = - sqrt(-(energy - echan)/ry)
C  The following change was done for unnatural parity
               noprint(nchp) = .false.
c$$$         noprint(nch) = .false.
            endif 
         end if
C  The following line puts the overlap of a pseudostate with the I^-
C  projection operator in to OVLPN which will be used in the CROSS program
         if (hlike) then
            if (nze.eq.1) then
               if (chan(nchp)(1:1).eq.'p') then
                  nt = 2 * la + 2
               else
                  nt = 2 * la + 1
               endif
            else 
               nt = la + 1
            endif 
         elseif(chan(nchp)(1:1).eq.'s'.or.chan(nchp)(1:1).eq.'S')then
            nt = 2 * la + 1
            if (Ps_start) stop "Must have Ps after all atomic states."
         else
            nt = 2 * la + 2
         endif
         if (ipar.eq.0) then
            ip = ipar
            if (chan(nchp)(1:1).eq.'T'.or.chan(nchp)(1:1).eq.'S') ip = 1
            if (nt.gt.nttop(ip)) nttop(ip) = nt
            npt(nt,ip) = npt(nt,ip) + 1
            if (npt(nt,ip).gt.npoints) stop
     >         'npt(nt,ip).gt.npoint: INCREASE NPOINTS IN PWRITE'
            enf(npt(nt,ip),nt,ip) = ea
            ovlpn(nch) = ovlp(na,la)
            if (chan(nchp)(1:1).eq.'p') then
               nestart = npbot(la)
               Ps_start = .true.
            else
               nestart = nabot(la)
            endif 
            do ne = nestart, max(nnbtop,nestart) + ne2e * 2
               ovlpnn(nch,ne) = ovlpnl(na,la,ne)
c$$$            print*,'Defining ovlpnn(nch,ne):',nch, ne,ovlpnn(nch,ne)
            enddo 
c$$$            print*,'Defining ovlpn(nch):',nch, ovlpn(nch)
            phaseq(nch) = phasen(na,la)
            eta1 = - (zasym + 1.0) / sqrt(abs(ea))
            zlmin = la
            call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp, sig,
     >         mode1,kfn,ifail)
            phasediv(nch) = exp(ci * sig(1))
            deta = - (zasym + 1.0) / sqrt(abs(ea))
            if (abs((phasediv(nch)-coulphase(deta,la))/
     >         coulphase(deta,la)).gt.1e-4) then
               print*,'PHASE checks:',phasediv(nch),coulphase(deta,la)
            endif 
c$$$            print*,'nt,n,phases:',nt,npt(nt,ip),phaseq(nch)/phasediv(nch)
c$$$     >         ,phasediv(nch),nch
            nchtop = nch
         endif
         nch = nch + 1
         call getchinfo (nch, nchp, ipar, temp, maxpsi, ea, la, na, l)
         if (nch.ne.0) go to 10
      enddo 
C  restore ipar to zero
      ipar = 0
      xerr = 0.1
      yerr = 0.0
      print*
      if (cnode.eq.'') then
         open(42,file='cfile')
         open(43,file='cfile1')
         open(44,file='dfile')
         open(45,file='dfile1')
      else
         open(42,file=cnode//'cfile'//ench)
         open(43,file=cnode//'cfile1'//ench)
         open(44,file=cnode//'dfile'//ench)
         open(45,file=cnode//'dfile1'//ench)
      endif 
      do ip = 0, 1
         do nt = 1, nttop(ip)
            do np = 1, npt(nt,ip)
               if (enf(np,nt,ip).gt.0.0) then
c$$$                  if (np.le.npt(nt,ip)-2) then
c$$$                     check = (enf(np,nt,ip)*4.0-enf(np+1,nt,ip)-
c$$$     >                  enf(np+2,nt,ip))/enf(np,nt,ip)
c$$$                  elseif (np.eq.npt(nt,ip)-1) then
c$$$                     check = (enf(np,nt,ip)*3.0-enf(np+1,nt,ip))/
c$$$     >                  enf(np,nt,ip)
c$$$                  else
c$$$                     check = 0.0
c$$$                  endif 
                  if (np.gt.1) then
                     write(42+ip,'(i2,1p,900e11.3)')nt,enf(np,nt,ip)*ry,
     >                  xerr,yerr,etot*ry,etot*ry/2.0,
     >                  (enf(np-1,nt,ip)+enf(np,nt,ip))/2.0*ry,
     >                  enf(np,nt,ip)*ry*2.0,
     >                  (slowery(ne)*ry,ne=1,ne2e)
                  else
                     write(42+ip,'(i2,1p,900e11.3)')nt,enf(np,nt,ip)*ry,
     >                  xerr,yerr,etot*ry,etot*ry/2.0,
     >                  enf(np,nt,ip)/2.0*ry,enf(np,nt,ip)*ry*2.0,
     >                  (slowery(ne)*ry,ne=1,ne2e)
                  endif 
               else
                  write(44+ip,'(i2,1p,100e13.5)')
     >               nt,enf(np,nt,ip)*ry,xerr,yerr,etot*ry
               endif 
            enddo
         enddo
      enddo 
      print*,'efiles written to disk'
      close(42)
      close(43)
      close(44)
      close(45)
      ip = 0
      
      nsp = nsmax + 1
      nparity = npar + 1

      ie2e = 1
      if (nchimax.gt.1.and.ne2e.gt.0) then
         print*,' More than one initial channel available, IE2E:'
         read*,ie2e
      endif 
      inquire(file=tfile,exist=exists)
      if (exists) then
         print'(''Will not calculate BornICS, potl file exists: '',
     >      a80)', tfile
         inquire(file='tcs',exist=exists)
         if (exists) then !tcs exists
            do ni = 1, knm
               do nf = 1, knm
                  BornICS(nf,ni) = 0.0
               enddo
            enddo 
C  Obtain the BornICS, after copying the 'totalcs' file to 'tcs',
C  from the 'tcs' file.
            csfile = 'tcs'
            esum(0:1) = 0.0
            rcond(0:1) = 0.0
            partcs(:,:,:) = 0.0
            sigtop(:,:) = 0.0
            BornPCS(:,:) = 0.0
            call wrtcs(partcs,sigtop,nchpmax,nent,instate,
     >         jstop,etot,nsmax,0,
     >         ovlpn,nunit,nznuc,zasym,projectile,target,canstop,
     >         BornICS,BornPCS,csfile,esum,rcond,nnbtop,ovlpnn)
            open(42,file='tcs')
            close(42,status='delete')
            open(42,file='tcs_J')
            close(42,status='delete')
            if (iborn.eq.0) then
               print*,'Setting BornICS=0 since iborn=0'
               BornICS(:,:) = 0.0
            endif 
         endif 

         irecl = max(100,28 * nchanmax * nchimax + 10)
         open(88,file=tfile,access='SEQUENTIAL',form='FORMATTED',
     >      recl=irecl)
         print*,'POTL and IRECL:',tfile,irecl
         call update(6)
c         open(88,file=tfile)
         read(88,*,err=25,end=25)
         read(88,*,err=25,end=25) enold, zasymold,nold
         if (abs(enold-energy)/energy.gt.1e-4) then
c$$$            print*,'Stopping: incident energy is not same in potl'
c$$$            call update(6)
c$$$            call mpi_finalize(ierr)
c$$$            print*,'MPI_FINALIZE returned:',ierr
c$$$            call update(6)
            stop 'incident energy is not same in potl'
         endif
         read(88,'(1000(i2,2(e12.4)))',err=25,end=25)
     >      (laold(n),ovlpnold(n),onshold(n),n=1,nold)
c$$$         read(88,*,err=25,end=25)
c$$$     >      (l,tmp,onshold(n),n=1,nold)
         do n = 1, min(nold,nopen(0))
            if (abs((onshold(n)-onshellk(n))/(onshellk(n)+1e-10))
     >         .gt.1e-4) then
               print*,n,onshold(n),onshellk(n)
               print*,'WARNING: onshellk are different'
c$$$               onshellk(n) = onshold(n)
c$$$               enchan(n) = etot - onshellk(n)**2
c$$$               enpsinb(n,laold(n)) = enchan(n)
c$$$               ovlpn(n) = ovlpnold(n)
            endif 
         enddo             
         read(88,*,err=25,end=25) (nopenold(np),np=0,npar), nparityold,
     >      nspold
         do np = 0, npar
            if (nopenold(np).ne.nopen(np)) then
               print*,' Caution: reset NOPEN(np) -> NOPENOLD(np)',
     >            nopen(np), nopenold(np), np
               nopen(np) = nopenold(np)
            endif
         enddo 
         if (nparityold.ne.nparity) stop 'Nparity must equal NPARITYOLD'
         if (nspold.ne.nsp) stop 'NSP must equal NSPOLD'

         mv = nchimax * nchanmax * nsp * nparity * (jstop + 1) * 8
         mt = mv
         me2e = (nsmax + 1) * (jstop + 1) * ( 2 * latop + 1) * ne2e * 2
     >      * ntype * 8
c$$$         call memalloc(ptrv,mv)
c$$$         call memalloc(ptrt,mt)
c$$$         call memalloc(ptrte2e,me2e)
         allocate(te2e(0:nsmax,0:jstop,-latop:latop,ntype,ne2e*2))
         allocate(vdon(nchanmax,nchimax,0:nsp-1,0:nparity-1,0:jstop))
         allocate( ton(nchanmax,nchimax,0:nsp-1,0:nparity-1,0:jstop))
c$$$         ton = 0.0
         do ne = 1, ne2e * 2
            do nt = 1, ntype
               do la = - latop, latop
                  do j = 0, jstop
                     do ns = 0, nsmax
                        te2e(ns,j,la,nt,ne) = (0.0,0.0)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         rnpow = 0
         rmpow = 0
         tics(0) = 0.0
         tics(1) = 0.0
         ticsi = 0.0
         ninterp = 0
         ninttype = 0
         if (ne2e.ne.0) then
            print*,'Enter interpolation type (-1, 0, 1, 2) and NINTERP'
            read*,ninttype,ninterp
            if (ninttype.eq.-1) then
               print*,'Will use the TCS section to define Te2e'
            else 
               print*,'Enter M for interpolation with e**M and lambda'
               read*, rmpow,alambda
               rnpow = 0.0
               if (ninttype.eq.0) then
                  print'("Interpolating real, imaginary and abs parts")'
               else if (ninttype.eq.1) then
                  print '(" Interpolating real and imaginary parts")'
               else if (ninttype.eq.2) then
                  print '(" Interpolating abs values and phases")'
               else if (ninttype.eq.3) then
                  print '("Will set T_1(2)=2/sqrt(3)(T_0(1)+T_0(2)/2)")'
               else if (ninttype.eq.4) then
                  print '("Will set T_0(2)=2/sqrt(3)(T_1(1)-T_1(2)/2)")'
               else if (ninttype.ge.5) then
                  print '(" Will sum singlet + triplet")'
               else 
                  stop 'Unknown type'
               endif
            endif 
            do ne = 1, ne2e
               xout(ne) = slowery(ne)
               xout(ne+ne2e) = etot - slowery(ne)
            enddo
c$$$         h = etot / ninterp
c$$$         do ne = 1, ninterp
c$$$            xout(ne + 2 * ne2e) = h * ne
c$$$            wout(ne + 2 * ne2e) = (mod(ne,2) * 2.0 + 2.0) / 3.0 * h
c$$$         enddo
c$$$         wout(ninterp + 2 * ne2e) = wout(ninterp + 2 * ne2e) / 2.0
            nwf = 2 * ninterp
            niwf = 2 * ninterp
            call cgqf(ninterp,xout(2*ne2e+1),wout(2*ne2e+1),1,0d0,0d0,
     >         0d0,dble(etot),0,nwf,wf,niwf,iwf,ier)
         
            do ne = 1, ninterp
               n = ne + 2 * ne2e
               if (abs((xout(n)+xout(ninterp+2*ne2e-ne+1))/etot-1d0)
     >            .gt.1d-5) print*,'XOUT not symmetric about E/2',
     >            xout(n)+xout(ninterp+2*ne2e-ne+1),etot,xout(n),wout(n)
            enddo 
         
            do ne = 1, ninterp + 2 * ne2e
               do nt = 1, ntype
                  sdcs(nt,ne) = 0.0
               enddo
            enddo 
            
            tests = 0.0
            do ne = 1, ninterp
               tests = tests + xout(ne + 2 * ne2e)**2*wout(ne + 2*ne2e)
            enddo
            print *, 'testing:',tests/etot**3*3.0
         endif 
         do ns = 0, 1 !nsmax
            call getspinw(chan(1)(1:1),ns,nze,spinw(ns),si)               
            tnbcs(ns) = 0.0
            tics(ns) = 0.0
            tcs(ns) = 0.0
            do nt = 1, ntype
               do nst = 1, knm
                  sigint(nst,nt,ns) = 0.0
               enddo
            enddo
         enddo 
         j = jstart

 15      if (j.gt.jstop) go to 20
         do np = 0, min(npar,j)
            ipar = np
            do n = 1, knm
               do la = 1, 2*lamax+1
                  non(n,la,np,j) = 0
               enddo
            enddo 
            do nchi = 1, nchan
               sigtop(nchi,0) = 0.0
               sigtop(nchi,1) = 0.0
               do nchf = 1, nchan
                  partcs(nchf,nchi,0) = 0.0
                  partcs(nchf,nchi,1) = 0.0
                  BornPCS(nchf,nchi) = 0.0
               enddo
            enddo 
            do ns = 0, nsmax
               if (ns.eq.0) read(88,'(1000i2)',err=20,end=20)
     >            (jlm(nch,np,j),nch=1,nopen(np))
               if (ns.eq.1) read(88,'(1000i2)',err=21,end=21)
     >            (jlm(nch,np,j),nch=1,nopen(np))
               if (nchanmax.lt.1000) then
               read(88,'(2i3,3x,10000i3)',err=22,end=22) lent,nymax,
     >            ((non(nch,l,np,j),l=1,jlm(nch,np,j)),nch=1,nopen(np))
               else
               read(88,'(2i4,3x,10000i4)',err=22,end=22) lent,nymax,
     >            ((non(nch,l,np,j),l=1,jlm(nch,np,j)),nch=1,nopen(np))
               endif
c$$$               if (ns.eq.0.and.np.eq.0) then
c$$$                  print*,'LENT, NCHIMAX:', lent, nchimax
c$$$                  print*,'NYMAX, NCHANMAX:',nymax,nchanmax
c$$$               endif
               write(fmt,'("(",i20,"(x,e12.5,x,e12.5,2x))")') lent*nymax
C  Can replace * with "fmt" below, see write(88 in solvet.
               read(88,*,err=23,end=23)
     >            ((vdon(nchf,nchi,ns,np,j),
     >            nchi=1,lent),nchf=1,nymax)
               read(88,*,err=24,end=24)
     >            ((ton(nchf,nchi,ns,np,j),
     >            nchi=1,lent),nchf=1,nymax) 
c$$$               if (j.eq.0.or.ns.eq.1) ton(:,:,ns,np,j)=(0.0,0.0) ! For Anatoli

c$$$  if (ntyp.lt.0) then
c$$$                  print*,'Will SET T = V for J => -1-NTYPE'
c$$$                  if (j.ge.-1-ntyp) then
c$$$                     print*,'setting T = V'
c$$$                     ton(1:nymax,1:lent,ns,np,j) =
c$$$     >                  vdon(1:nymax,1:lent,ns,np,j)
c$$$                  endif 
c$$$               endif 
C  Initialise the nfich and latom arrays
               do nch = 1, nchan
                  nfich(nch) = 0
               enddo 
               nch = 1          !nymax
               jltmp = -1 + ipar
               call getchinfo (nch, nchp, j, temp, maxpsi, ea,
     >            la, na, l)
               if (nchfi(nchp,ipar).gt.0) 
     >            nonnew(nchfi(nchp,ipar),j-l,ipar,j) = nch
               jltmp = jltmp + 2
               do while (nch.ne.0)
                  nchpold = nchp
                  nchtop = nch
                  nch = nch + 1
                  call getchinfo (nch, nchp, j, temp, maxpsi, ea,
     >               la, na, l)
                  if (nchp.eq.nchpold) then
                     jltmp = jltmp + 2
                  else
                     jltmp = -1 + ipar + 2
                  endif 
                  if (nchfi(nchp,ipar).gt.0.and.nch.gt.0) then
                     nonold = non(nchfi(nchp,ipar),jltmp,ipar,j)
                     nonnew(nchfi(nchp,ipar),j-l,ipar,j) = nonold
                     if (nonold.ne.nch.and.nch.gt.0) then
c$$$  print*,'Setting NON to NCH for NCHP:',
c$$$  >                     nonold,nch,nchp
                        non(nchfi(nchp,ipar),jltmp,ipar,j) = nch
                        nonnew(nchfi(nchp,ipar),j-l,ipar,j) = nch
c$$$  jlm(nchfi(nchp,ipar),ipar,j) = jltmp
                     endif
                  endif 
               enddo
               do nchi = 1, lent
                  do nchf = nymax + 1, nchtop
                     vdon(nchf,nchi,ns,np,j) = (0.0,0.0)
                     ton(nchf,nchi,ns,np,j)  = (0.0,0.0)
                  enddo 
               enddo 

               do nch = 1, nchtop
                  nchpf = nfich(nch)
                  if (onshellk(nchpf).lt.0.0) then
                     do nchi = 1, lent
                        do nchf = nchtop, nch+1, -1
                           vdon(nchf,nchi,ns,np,j) =
     >                        vdon(nchf-1,nchi,ns,np,j)
                           ton(nchf,nchi,ns,np,j) =
     >                        ton(nchf-1,nchi,ns,np,j)
                        enddo
                        vdon(nch,nchi,ns,np,j) = (0.0,0.0)
                        ton(nch,nchi,ns,np,j)  = (0.0,0.0)
                     enddo
                  endif
               enddo 
                     
               if (exists) then !tcs exists
               nchpiprev = 0
               do nchi = 1, lent
                  call getchinfo (nchi, nchp, j, temp, maxpsi, eia,
     >               la, na, l)
                  nchpi = nfich(nchi)
                  if (nchpi.ne.nchp) stop 'nchpi.ne.nchp'
                  lia = latom(nchpi)
                  if (lia.ne.la) stop 'lia.ne.la'
                  if (onshellk(nchpi).gt.0.0) then
                     deta = nze * zasym / onshellk(nchpi)
                     phase = coulphase(deta,l)
                     call getspinw(chan(nchpi)(1:1),ns,nze,spinw(ns),si)  
                     const = spinw(ns) * (2.0 * pi) ** 4 / (4.0 * pi) *
     >                  (2.0 * j + 1.0) / (2.0 * lia + 1.0)
                     sigtop(nchpi,ns) = sigtop(nchpi,ns) - const *
     >                  aimag(ton(nchi,nchi,ns,np,j)/phase**2) / 
     >                  (pi * onshellk(nchpi))
c$$$                     print*,'Caution T=T/phase**2 for Sonja',phase
c$$$                     ton(nchi,nchi,ns,np,j) = ton(nchi,nchi,ns,np,j) /
c$$$     >                  phase**2
                     if (nchpi.ne.nchpiprev) tcs(ns) = 0.0
                     nchpiprev = nchpi
                     do nchf = 1, nchtop
                        nchpf = nfich(nchf)
                        if (onshellk(nchpf).gt.0.0) then
                           if (nchpf.eq.nchpi) then
                              call getchinfo (nchf, nchp, j, temp, 
     >                           maxpsi, ef, la, na, lf)
c$$$                              if (j.eq.0.and.ns.eq.0) write(100+nchpi,
c$$$     >                           '(f10.4," eV on ",a3/
c$$$     >                           12x,"T",16x,"Phase(Li)",14x,
c$$$     >                           "Phase(Lf)",7x,"J Li Lf P S")')
c$$$     >                           onshellk(nchpi)**2*ry,chan(nchpi)
c$$$                              write(100+nchpi,'(1p,6e11.3,3i3,2i2)')
c$$$     >                        ton(nchf,nchi,ns,np,j),coulphase(deta,l),
c$$$     >                           coulphase(deta,lf),j,l,lf,np,ns
                           endif 
                        lfa = latom(nchpf)
                        const = spinw(ns) * (2.0 * pi) ** 4 / 
     >                     (4.0 * pi)*onshellk(nchpf)/onshellk(nchpi) *
     >                     (2.0 * j + 1.0) / (2.0 * lia + 1.0)
C  The following statement aims to get good extrapolation for positronium
C  channels after say J=20, even though these cross sections are set to zero
C  for J > 8. Note that this affects the accuracy of the TCS test using the
C     optical theorem, since the imaginary part of the elastic amplitude
C     is not altered. 
                        if (j.gt.8) then 
                           if (abs(ton(nchf,nchi,ns,np,j)).eq.0.0 .and.
     >                        abs(ton(nchf,nchi,ns,np,j-1)).ne.0.0)then
                              q = abs(ton(nchf,nchi,ns,np,j-1))**2/
     >                           abs(ton(nchf,nchi,ns,np,j-2)+1e-10)**2
                              if (q.gt.1.0) q = 0.0
                              ton(nchf,nchi,ns,np,j)=sqrt(q) *
     >                           ton(nchf,nchi,ns,np,j-1)
                           endif
                        endif 
                        partcs(nchpf,nchpi,ns) = partcs(nchpf,nchpi,ns) 
     >                     + abs(ton(nchf,nchi,ns,np,j))**2 * const
                        tcs(ns) = tcs(ns)
     >                     + abs(ton(nchf,nchi,ns,np,j))**2 * const
                        BornPCS(nchpf,nchpi) = BornPCS(nchpf,nchpi) +
     >                     abs(vdon(nchf,nchi,ns,np,j))**2 * const
c$$$                        if(nchpf  .eq. 18 .and. nchpi  .eq. 14) then
c$$$c !!!
c$$$                           write(*,*)
c$$$     >                        j,np,nchf,nchi,ns,
c$$$     >                        ton(nchf,nchi,ns,np,j),
c$$$     >                        vdon(nchf,nchi,ns,np,j)
c$$$                        end if
                     endif 
                     enddo
                     if (nfich(nchi+1).ne.nchpi.and.tcs(ns).gt.0.0) then
                       if(abs(sigtop(nchpi,ns)/tcs(ns)-1.0).gt.1e-4)then
                           print*,'ratio of tcs not 1 (Ps extrap?):',
     >                       sigtop(nchpi,ns)/tcs(ns),
     >                       tcs(ns),sigtop(nchpi,ns), nchpi
                           print*,'Setting SIGTOP=SIGT'
                           sigtop(nchpi,ns) = tcs(ns)
                        endif 
                     endif
                  endif 
               enddo

               do nch = 1, nchtop
                  call getchinfo (nch, nchp, j, temp, maxpsi, ea,
     >               la, na, l)
                  ch = chan(nchp)(1:1)
                  if (hlike) then
                     if (nze.eq.1) then
                        if (ch.eq.'p') then
                           nt = 2 * la + 2
                        else
                           nt = 2 * la + 1
                        endif
                     else 
                        nt = la + 1
                     endif 
                  else if (ch.eq.'s'.or.ch.eq.'S') then
                     nt = 2 * la + 1
                  else if (ch.eq.'t'.or.ch.eq.'T') then
                     nt = 2 * la + 2
                  endif
                           
c$$$                  if (ne2e.gt.0) print*,'Defining Te2e for ne2e:',ne2e
                  if (ninttype.eq.-1) then
                  do ne = 1, ne2e * 2
                     if (j.eq.0) then
                        eta1 = - (zasym + 1.0) / sqrt(xout(ne))
                        zlmin = la
                        call coulcc(xx,eta1,zlmin,nl,cfc,cgc,
     >                     cfcp,cgcp, sig,mode1,kfn,ifail)
                        phasetim(nt,ne) = exp(ci * sig(1))
                     endif 
                     te2e(ns,j,j-l,nt,ne) = te2e(ns,j,j-l,nt,ne) +
     >                  ovlpnn(nchp,nnbtop+ne) / sqrt(sqrt(xout(ne))) *
     >                  ton(nch,ie2e,ns,np,j) * phasetim(nt,ne)
                     if (ch.eq.'S'.or.ch.eq.'T') 
     >                  print*,'ns,nt,ne,|ovlp|,|t|:',ns,nt,ne,
     >                  abs(ovlpnn(nchp,nnbtop+ne)),
     >                  abs(ton(nch,ie2e,ns,np,j))
                     if (ne.eq.1.and.abs(slowery(ne)-ea).lt.1e-4) then
                        xx=ovlpnn(nchp,nnbtop+ne)/sqrt(sqrt(xout(ne)))*
     >                  ton(nch,ie2e,ns,np,j)*phasetim(nt,ne)
                        print'("J,S,L,l,|T|,arg(T)",
     >                     4i3,1x,e11.4,f7.3)',j,ns,l,la,abs(xx),
     >                     atan2(aimag(xx),real(xx))
                     endif 
                  enddo
                  endif ! ninttype.eq.-1
               enddo
               endif ! file tcs exists

               if (ne2e.ne.0.and.np.eq.0) then
                  call getspinw(chan(1)(1:1),ns,nze,spinw(ns),si)
                  do la = -latop, latop
                     do nt = 1, ntype
                        npt(nt,la) = 0
                     enddo
                  enddo 
                  ip = 0
                  nttop(ip) = 0
C  The following mess with JL, NCHF, NYMAX, NON etc is due to the fact
C  that there may be a mismatch between the read data and NCHF coming
C  from GETCHINFO. This happens when there is a closed channel in between
C  open ones.
                  nchprev = 1
                  jl = -1
                  nchf = 1
                  do la = - lamax, lamax
                     do nt = 1, ntype
                        sigintj(nt,la) = 0.0
                     enddo 
                  enddo
                  ntfun(:,:) = -1
                  lia = 0
                  call getchinfo (nchf, nchp, j, temp, maxpsi, ea,
     >               lfa, na, l)
                  do while (nchf.gt.0)  ! do nchf=1,nymax
                     if (ea.gt.0.0.and.onshellk(nchp).gt.0.0
     >                  .and.chan(nchp)(1:1).ne.'T'
     >                  .and.chan(nchp)(1:1).ne.'S') then
                        const = spinw(ns) * (2.0 * pi) ** 4/(4.0 * pi)*
     >                     (2.0 * j + 1.0) * onshellk(nchp)/onshellk(1) 
     >                     / (2.0 * lia + 1.0)
                        if (hlike) then
                           if (nze.eq.1) then
                              if (chan(nchp)(1:1).eq.'p') then
                                 nt = 2 * lfa + 2
                              else
                                 nt = 2 * lfa + 1
                              endif
                           else 
                              nt = lfa + 1
                           endif 
                           nse2e(nt) = 2
                        else if (chan(nchp)(1:1).eq.'s') then
                           nt = 2 * lfa + 1
                           nse2e(nt) = 1
                        else
                           nt = 2 * lfa + 2
                           nse2e(nt) = 3
                        endif
                        if (nchp.eq.nchprev) then
                           jl = jl + 2
                        else
                           jl = 1
                        endif
                        nchprev = nchp
                        nfold = non(nchfi(nchp,np),jl,np,j)
                        nf = nonnew(nchfi(nchp,np),j-l,np,j)
                        if (nf.ne.nfold) stop 'nonnew.ne.non, check!'
                        chanf(nt) = chan(nchp)
                        le2e(nt) = lfa
                        ntfun(lfa,nse2e(nt)) = nt
                        if (nt.gt.nttop(ip)) nttop(ip) = nt
                        npt(nt,j-l) = npt(nt,j-l) + 1
                        if (npt(nt,j-l).gt.npoints) stop
     >              'npt(nt,j-l).gt.npoints: INCREASE NPOINTS IN PWRITE'
                        csi = abs(ton(nf,ie2e,ns,np,j)) ** 2 * const
                        tics(ns) = tics(ns) + csi
                        sigint(nchp,nt,ns) = sigint(nchp,nt,ns) + csi
                        sigintj(nt,j-l) = sigintj(nt,j-l) + csi
                        enf(npt(nt,j-l),nt,j-l) = ea
                        treal = real(ton(nf,ie2e,ns,np,j)*ovlpn(nchp) *
     >                     phaseq(nchp)/phasediv(nchp)/sqrt(sqrt(ea)))
                        timag = aimag(ton(nf,ie2e,ns,np,j)*ovlpn(nchp)*
     >                     phaseq(nchp)/phasediv(nchp)/sqrt(sqrt(ea)))
                        tabs = abs(cmplx(treal,timag))
                        tin(npt(nt,j-l),nt,j-l,1) = treal
                        tin(npt(nt,j-l),nt,j-l,2) = timag
                        tin(npt(nt,j-l),nt,j-l,3) = tabs
                        if (tabs.eq.0.0) then
c$$$                           print*,'T input is zero!!!',ovlpn(nchp),
c$$$     >                        ton(nf,ie2e,ns,np,j),phaseq(nchp)
                        else 
                           tin(npt(nt,j-l),nt,j-l,4)=atan2(timag,treal)
                        endif 
                     endif
                     nchf = nchf + 1
                     call getchinfo (nchf, nchp, j, temp, maxpsi, ea,
     >                  lfa, na, l)
                  enddo ! end of nchf loop
                  do nt = 1, nttop(ip)
                     lfa = le2e(nt)
                     do l = abs(j - lfa), j + lfa, 2
                        xin(1) = 0d0
                        yin(1) = 0d0
                        nfinal = npt(nt,j-l) + 2
                        yin(nfinal) = 0d0
                        xin(nfinal) = etot
c$$$                        do n = 1, npt(nt,j-l)
c$$$                           print'("nt,n,abs(T),arg(T),T",2i3,e10.2,f6.2,
c$$$     >                        2e10.3)',nt,n,abs(
c$$$     >                        cmplx(tin(n,nt,j-l,1),tin(n,nt,j-l,2))),
c$$$     >                        atan2(tin(n,nt,j-l,2),tin(n,nt,j-l,1)),
c$$$     >                        tin(n,nt,j-l,1),tin(n,nt,j-l,2)
c$$$                        enddo
c$$$                        print*
                        do i = 1, 3
                           do n = 1, npt(nt,j-l)
                              xin(n+1) = enf(n,nt,j-l)
                              fac = enf(n,nt,j-l) ** rmpow *
     >                           (etot-enf(n,nt,j-l)) ** rnpow
                              yin(n+1) = tin(n,nt,j-l,i) * fac
                           enddo
                           if (alambda.gt.0.0) then
C  x = (e - lambda^2/8)/(e + lambda^2/8), where e is in a.u.
                              al = alambda
                              do n = 1, nfinal
                                 xin(n) = (xin(n)/2.0 - al*al/8.0)/
     >                              (xin(n)/2.0 + al*al/8.0)
                              enddo
                              do n = 1, 2*ne2e+ninterp
                                 xout(n) = (xout(n)/2.0 - al*al/8.0)/
     >                              (xout(n)/2.0 + al*al/8.0)
                              enddo
                           endif 
                           call intrpl(nfinal,xin,yin,2*ne2e+ninterp,
     >                        xout,yout(1,i))
                           if (alambda.gt.0.0) then
                              do n = 1, nfinal
                                 xin(n) = al*al/4.0*(xin(n) + 1.0)/
     >                              (1.0 - xin(n)/2.0)
                              enddo
                              do n = 1, 2*ne2e+ninterp
                                 xout(n) = al*al/4.0*(xout(n) + 1.0)/
     >                              (1.0 - xout(n))
                              enddo
                           endif 

c$$$                           if (i.eq.3.and.l.lt.10.and.j.lt.10) then
c$$$                              write(sdcsfile,
c$$$     >                           '(i1,".",i1,"-",i1,".",i1)') j,l,lfa,ns
c$$$                              open(17,file=sdcsfile)
c$$$                              do n = 1, nfinal
c$$$                                 write(17,*) xin(n),yin(n)
c$$$                              enddo
c$$$                              write(17,*)
c$$$                              do n = 2*ne2e+1, 2*ne2e+ninterp
c$$$                                 write(17,*) xout(n),yout(n,i)
c$$$                              enddo
c$$$                              close(17)
c$$$                           endif 
                        enddo 
                        if (ninttype.eq.2) then
                           i = 4
                           do n = 1, npt(nt,j-l)
                              xin(n) = enf(n,nt,j-l)
                              yin(n) = tin(n,nt,j-l,i)
                           enddo
                           call intrpl(npt(nt,j-l),xin,yin,
     >                        2*ne2e+ninterp,xout,yout(1,i))
                        endif 
                        do ne = 1, 2 * ne2e
                           if (j.eq.0) then
                              eta1 = - (zasym + 1.0) / sqrt(xout(ne))
                              zlmin = lfa
                              call coulcc(xx,eta1,zlmin,nl,cfc,cgc,
     >                           cfcp,cgcp, sig,mode1,kfn,ifail)
                              phasetim(nt,ne) = exp(ci * sig(1))
                           endif 
                           
                           fac = xout(ne)**rmpow*(etot-xout(ne))**rnpow
                           if (ninttype.eq.0.or.ninttype.gt.2) then
                              yabs = abs(cmplx(yout(ne,1),yout(ne,2)))
                              te2e(ns,j,j-l,nt,ne) = phasetim(nt,ne) *
     >                           cmplx(yout(ne,1)/fac,yout(ne,2)/fac) *
     >                           yout(ne,3) / (yabs + 1e-30)
                           else if (ninttype.eq.1) then
                              te2e(ns,j,j-l,nt,ne) = phasetim(nt,ne) *
     >                           cmplx(yout(ne,1)/fac,yout(ne,2)/fac)
                           else if (ninttype.eq.2) then
                              te2e(ns,j,j-l,nt,ne) = phasetim(nt,ne) *
     >                           yout(ne,3)/fac*exp(ci*yout(ne,4))
                           endif
                           if (ne.eq.1.and.
     >                        abs((etot-xout(ne))/etot-0.5).lt.1e-2)then
                              rmax = rmesh(meshr,1)
                              rk1 = sqrt(xout(ne))
                              rkappa = sqrt(etot)
c$$$                              te2e(ns,j,j-l,nt,ne)=te2e(ns,j,j-l,nt,ne)*
c$$$     >                          exp(ci/rk1*log(2.0*rk1*rk1*rmax/rkappa))
                              print'("J,ns,L,l,state,T:",4i3,a4,1p,
     >                           e13.3,e11.3)', j,ns,l,lfa,chanf(nt),
     >                           real(te2e(ns,j,j-l,nt,ne)),
     >                           aimag(te2e(ns,j,j-l,nt,ne))!,
c$$$     >                           atan2(aimag(te2e(ns,j,j-l,nt,ne)),
c$$$     >                           real(te2e(ns,j,j-l,nt,ne)))
c$$$                              if (j.eq.0) then
c$$$                                 print*,'Te2e = 0 for J, L, l, S:',
c$$$     >                              j,l,lfa,ns
c$$$                                 te2e(ns,j,j-l,nt,ne) = (0.0,0.0)
c$$$                              endif 
                           endif

                        enddo
c$$$                     do n = 1, npt(nt,j-l)
c$$$                        const  = spinw(ns) * (2.0 * pi) ** 4/(4.0 * pi)*
c$$$     >                     (2.0 * j + 1.0) * sqrt(etot-enf(n,nt,j-l))
c$$$     >                     / onshellk(1) / (2.0 * lia + 1.0)
c$$$                        print*,'NT,n,|T|',nt,n,units(nunit) * const *
c$$$     >                     abs(cmplx(tin(n,nt,j-l,1),
c$$$     >                     tin(n,nt,j-l,2)))**2*
c$$$     >                     sqrt(enf(n,nt,j-l))
c$$$                     enddo
                        testi = 0.0
                        do ne = 1, ninterp
                           n = ne + 2 * ne2e
                           fac = xout(n)**rmpow*(etot-xout(n))**rnpow
                           const  = spinw(ns) * (2.0 * pi)**4/(4.0*pi)*
     >                        (2.0 * j + 1.0) * sqrt(etot-xout(n))
     >                        / onshellk(1)/(2.0 * lia + 1.0)
                           if (ninttype.eq.1) then
                              abst = abs(cmplx(yout(n,1),yout(n,2)))
                           else 
                              abst = yout(n,3)
                           endif 
                           testi = testi+abst**2
     >                        * wout(n) * const / fac**2 / 2.0
                           sdcs(nt,n) = sdcs(nt,n) +
     >                        abst**2 * const / fac**2
                           psdcs(nt,n,ns) = abst**2 * const / fac**2 !atan2(yout(n,2),yout(n,1))
                        enddo
c$$$                     print*,'J, NS, NT, TICS:',j,ns,nt,
c$$$     >                  sigintj(nt,j-l),testi,
c$$$     >                  abs(sigintj(nt,j-l)-testi)/(testi+1e-30)
                     enddo ! end of the l loop
                  enddo ! end of the nt loop
                  if (abs((etot-xout(1))/etot-0.5).lt.1e-2
     >               .and.nze.eq.-1) then
                     if (ninttype.ge.5) te2e(ns,j,:,:,2) = (0.0,0.0)
                     te2ec(:,:) = (0.0,0.0)
                     do nt = 1, nttop(0)
                        if (hlike) then
                           lfa = nt - 1
                           do l = abs(j-lfa),j+lfa,2
                              if (l.le.latop) then     
                                 te2ec(j-l,nt) =
     >                              te2e(ns,j,j-l,nt,1)+(-1)**ns 
     >                              * te2e(ns,j,j-lfa,l+1,1)
                                 if (ninttype.eq.6) then
                                    if (lfa.lt.l) then
                                       te2ec(j-l,nt) =
     >                                    2.0*te2e(ns,j,j-l,nt,1)
                                    elseif (lfa.gt.l) then
                                       te2ec(j-l,nt) =
     >                                    2.0 * (-1)**ns 
     >                                    * te2e(ns,j,j-lfa,l+1,1)
                                    endif
                                 elseif (ninttype.eq.7) then
                                    if (lfa.gt.l) then
                                       te2ec(j-l,nt) =
     >                                    2.0*te2e(ns,j,j-l,nt,1)
                                    elseif (lfa.lt.l) then
                                       te2ec(j-l,nt) =
     >                                    2.0 * (-1)**ns 
     >                                    * te2e(ns,j,j-lfa,l+1,1)
                                    endif
                                 endif 
                              endif 
                           enddo
                        else
                           if (nse2e(nt).eq.1) then
                              nse2enp = 3
                              fac = sqrt(3.0)
                           else
                              nse2enp = 1
                              fac = 1.0/sqrt(3.0)
                           endif
                           lfa = le2e(nt)
                           do l = abs(j-lfa),j+lfa,2
                              if (l.gt.latop) cycle
                              ntp = ntfun(l,nse2enp)
                              nt2 = ntfun(l,nse2e(nt))
                              if (ntp.le.0.or.nt2.le.0) cycle
                              if (nse2enp.eq.3) then
                                 f0(1) = te2e(ns,j,j-l,nt,1)
                                 if (ntfun(l,1).gt.0)
     >                              f0(2)=te2e(ns,j,j-lfa,ntfun(l,1),1)
                                 if (ntfun(lfa,3).gt.0)
     >                              f1(1)=te2e(ns,j,j-l,ntfun(lfa,3),1)
                                 if (ntfun(l,3).gt.0)
     >                              f1(2)=te2e(ns,j,j-lfa,ntfun(l,3),1)
                                 ce2e = (te2e(ns,j,j-lfa,ntp,1) - 0.5 *
     >                              te2e(ns,j,j-l,ntp,1))*sqrt(4./3.)
                              else
                                 f1(1) = te2e(ns,j,j-l,nt,1)
                                 if (ntfun(l,3).gt.0)
     >                              f1(2)=te2e(ns,j,j-lfa,ntfun(l,3),1)
                                 if (ntfun(lfa,1).gt.0)
     >                              f0(1)=te2e(ns,j,j-l,ntfun(lfa,1),1)
                                 if (ntfun(l,1).gt.0)
     >                              f0(2)=te2e(ns,j,j-lfa,ntfun(l,1),1)
                                 ce2e = (te2e(ns,j,j-lfa,ntp,1) + 0.5 *
     >                              te2e(ns,j,j-l,ntp,1))*sqrt(4./3.)
                              endif 
                              fe2e(0)=f0(1)-f0(2)/2.0+sqrt(0.75)*f1(2)
                              fe2e(1)=f1(1)+sqrt(0.75)*f0(2)+f1(2)/2.0
                              testc = abs(fe2e(0))**2+abs(fe2e(1))**2
                              testi = 2.0*(abs(f0(1))**2+abs(f0(2))**2
     >                           +abs(f1(1))**2+abs(f1(2))**2)
c$$$                              print'("J,ns,L,l:",4i3,a4,
c$$$     >                           f7.2)', j,ns,l,lfa,
c$$$     >                           chanf(nt),abs((testc-testi)/testc)*100.
C  The following was a test of the symmetrisation inside the TDCS routine for
C  equal energy-sharing e-He ionization. It tested ok.
                              if (ninttype.eq.3.and.nse2enp.eq.3) then
                                 print*,'T_1(2), T_0(1), T_0(2)',
     >                              te2e(ns,j,j-lfa,ntp,1),
     >                              te2e(ns,j,j-l,nt,1),
     >                              te2e(ns,j,j-lfa,nt2,1)
                                 
c$$$                                 te2e(ns,j,j-lfa,ntp,1) =
c$$$     >                              2.0/sqrt(3.0) *
c$$$     >                              (te2e(ns,j,j-l,nt,1) + 0.5 *
c$$$     >                              te2e(ns,j,j-lfa,nt2,1))
c$$$                                 
c$$$                                 te2e(ns,j,j-lfa,ntp,2) =
c$$$     >                              te2e(ns,j,j-lfa,ntp,1)
c$$$                                 print*,'T_1(2) changed to',
c$$$     >                              te2e(ns,j,j-lfa,ntp,1)
                                 te2e(ns,j,j-l,nt,1) = 2.0/sqrt(3.0) *
     >                              (f1(2) - 0.5 * f1(1))
                                 te2e(ns,j,j-l,nt,2)=te2e(ns,j,j-l,nt,1)
                              endif
                              if (ninttype.eq.4.and.nse2enp.eq.1) then
                                 print*,'f0(2), f0(1)',f0(2),f0(1)
c$$$                                 
c$$$                                 te2e(ns,j,j-lfa,ntp,1) =
c$$$     >                              2.0/sqrt(3.0) *
c$$$     >                              (te2e(ns,j,j-l,nt,1) - 0.5 *
c$$$     >                              te2e(ns,j,j-lfa,nt2,1))
c$$$                                 
c$$$                                 te2e(ns,j,j-lfa,ntp,2) =
c$$$     >                              te2e(ns,j,j-lfa,ntp,1)
c$$$                                 print*,'T_0(2) changed to',
c$$$     >                              te2e(ns,j,j-lfa,ntp,1)
                                 te2e(ns,j,j-l,nt,1) = 2.0/sqrt(3.0) *
     >                              (f0(2) + 0.5 * f0(1))
                                 te2e(ns,j,j-l,nt,2)=te2e(ns,j,j-l,nt,1)
                                 
                              endif
                              if (ninttype.eq.5) then
c$$$                                 te2e(ns,j,j-l,nt,2) =
c$$$     >                              te2e(ns,j,j-l,nt,1)
c$$$     >                              + te2e(ns,j,j-lfa,ntp,1)
c$$$                                 sumc = abs(te2e(ns,j,j-l,nt,2))**2
c$$$                                 sumi = (abs(te2e(ns,j,j-l,nt,1))**2+
c$$$     >                              abs(te2e(ns,j,j-lfa,ntp,1))**2)
c$$$                                 print*,'1.9(|T(0)|^2+|T(1)|^2), ',
c$$$     >                              '|T(0)+T(1)|^2:',1.866*sumi,sumc,
c$$$     >                              100.0*abs((1.866*sumi-sumc)/sumc)
                                 if (nse2enp.eq.3) then
                                    te2e(ns,j,j-l,nt,2) =
     >                                 te2e(ns,j,j-l,nt,1) -
     >                                 0.5 * te2e(ns,j,j-lfa,nt2,1) +
     >                                 sqrt(0.75)*te2e(ns,j,j-lfa,ntp,1)
                                 else
                                    te2e(ns,j,j-l,nt,2) =
     >                                 te2e(ns,j,j-l,nt,1) +
     >                                 0.5 * te2e(ns,j,j-lfa,nt2,1) +
     >                                 sqrt(0.75)*te2e(ns,j,j-lfa,ntp,1)
                                 endif 
                              endif 
c$$$                              te2e(ns,j,j-lfa,ntp,1) =
c$$$     >                           te2e(ns,j,j-l,nt,1)
c$$$                              te2e(ns,j,j-lfa,ntp,2) =
c$$$     >                           te2e(ns,j,j-l,nt,2)
                           enddo  ! l loop
                        endif ! hlike
                     enddo ! nt loop
                     if (ninttype.ge.5.and.hlike) then
                        te2e(ns,j,:,:,:) = (0.0,0.0)
C  The following was a test of the symmetrisation inside the TDCS routine for
C  equal energy-sharing e-H ionization. It tested just fine!
                        do nt = 1, nttop(0)
                           lfa = nt - 1
                           do l = abs(j-lfa),j+lfa,2
                              if (l.le.latop) then
                                 print*,'J,ns,L,l:',j,ns,l,lfa,
     >                              te2e(ns,j,j-l,nt,1),
     >                              te2e(ns,j,j-lfa,l+1,1),
     >                              te2ec(j-l,nt)
                              endif 
                              te2e(ns,j,j-l,nt,1) = te2ec(j-l,nt)
                           enddo
                        enddo
                     endif 
                  endif         ! end of equal-energy-sharing
               endif            ! end of the ne2e .ne. 0 condition
            enddo ! end of ns loop
            if (exists) then
               csfile = 'tcs'
               call wrtcs(partcs,sigtop,nchpmax,nent,instate,
     >            j,etot,nsmax,np,
     >            ovlpn,nunit,nznuc,zasym,projectile,target,canstop,
     >            BornICS,BornPCS,csfile,esum,rcond,nnbtop,ovlpnn)
            endif 
         enddo ! end of np loop
         ipar = 0  ! reset ipar
         elrealt = real(ton(1,1,0,0,j))
         elimagt = aimag(ton(1,1,0,0,j))
c$$$         if (abs(elimagt/(elrealt+1e-30)).lt.0.5) then
            polfac = - elrealt * (2*j+3) * (2*j+1) * (2*j-1) /
     >         onshellk(1)
            if (target .eq. 'H  I') then
               exact = 4.5
            else if (target .eq. 'He I') then
               exact = 1.3832
            else if (target .eq. 'Li I') then
               exact = 161.77
            else if (target .eq. 'Na I') then
               exact = 164.7
            else if (target .eq. 'K I') then
               exact = 279.7
            else if (target .eq. 'Cs I') then
               exact = 401.0
            else
               exact = 0.0
            endif 
            print'(i3,2f10.3,f10.5,
     >         '' J, target polarization, exact, Imag(T)/Re(T)'')',
     >         j, polfac, exact, elimagt/(elrealt+1e-30)
c$$$         else
c$$$            print'(i3, '' read'')',j
c$$$         endif
         j = j + 1
         go to 15
 20      print*  ! end of J loop
         if (j.le.jstop) then
            close(88)
            open (88,file=tfile,position='append',status='old',
     >         recl=irecl)
            jstart = j
c$$$            call memfree(ptrv)
c$$$            call memfree(ptrt)
            deallocate(ton,stat=istat1)
            deallocate(vdon,stat=istat2)
            deallocate(te2e,stat=istat3)
            print*,'Deallocated Ton, Vdon, Te2e',istat1,istat2,istat3
         else
            print'("Interpolation powers M, N:",2f3.1)', rmpow, rnpow
            do nt = 1, nttop(ip)
               do nst = 1, knm
                  sum = 0.0
                  do ns = 0, nsmax
                     sum = sum + sigint(nst,nt,ns)
                  enddo 
                  if (sum.gt.0) print '(a4,1p,e10.3)',
     >               chan(nst), sum * units(nunit)
               enddo
            enddo 
            nn = max(1,iborn)
            open (47,file='sdcs.potl')
            open (48,file='psdcs.potl')
            ticsi = 0.0
            if (nunit.eq.3) then
               ufac = units(nunit) / ry / 2.0
            else
               ufac = 1.0
            endif 
            do ne = 1, ninterp
               n = ne + ne2e * 2
               sdcst = 0.0
               do nt = 1, nttop(ip)
                  sdcst = sdcst + sdcs(nt,n) 
     >               + sdcs(nt, ninterp + 2 * ne2e + 1 - ne)
               enddo
               ticsi = ticsi + sdcst * wout(n) / 2.0 * units(nunit)
               write(47,'(f8.2,1p,100e10.2)') xout(n)*ry, sdcst * ufac,
     >            (sdcs(nt,n)* ufac,  nt = 1, nttop(ip))
c$$$               write(48,'(f8.2,100f10.4)') xout(n)*ry,
c$$$     >            ((psdcs(nt,n,ns),  nt = 1, nttop(ip)), ns=0,nsmax)
c$$$     >            ((sdcs(nt,n)+sdcs(nt, ninterp + 2 * ne2e + 1 - ne))
c$$$     >            * ufac,  nt = 1, nttop(ip))
            enddo 
            write(47,*)
            print*,'Summed and integrated TICS:',
     >         (tics(0)+tics(1)) * units(nunit), ticsi/2.0
            do ne = 1, ne2e 
               do nt = 1, nttop(ip)
                  print*,'Writing amp for ip,nt,le2e(nt),latop:',
     >               ip,nt,le2e(nt),latop
                  
c     Alisher's changes               
                  if (chanf(nt)(1:1).eq.'p') then
                     rkf = sqrt(2.0*(etot-slowery(ne)))
                  else 
                     rkf = sqrt(etot-slowery(ne))
                  endif
                  
                  call makee2eamp(te2e(0,0,-latop,nt,ne),
     >               te2e(0,0,-latop,nt,ne+ne2e),nsmax,jstop,
     >               latop,le2e(nt),slowery(ne),rkf, ! sqrt(etot-slowery(ne))
     >               onshellk(ie2e),nunit,nt,ne,sweight(0,nn),nstep(nn),
     >               projectile,target,ninttype,rmpow,rnpow,chan(ie2e),
     >               chanf(nt),psdcs(nt,ne,0:1),psdcs(nt,ne+ne2e,0:1),
     >               chun,nze)
               enddo
C  Note that for e in Rydbergs we have dq = de / (2 * q), so to convert
C  atomic units to "/eV" we divide by (2*Ry), with the q division being
C  above. For e in atomic units dq = de / q, and so no division is necessary.
C  The division by 27.21 eV is done in the case we have cm^2/eV units in
C  the makee2eamp routine.
               print'(20(3x,a2,4x))',(chanf(nt)(1:1)//chanf(nt)(3:3),
     >            nt = 1, nttop(ip))
               print'(1p,20e9.2)',
     >            (psdcs(nt,ne,0)+psdcs(nt,ne,1), nt = 1, nttop(ip))
c$$$     >            (sdcs(nt,ne), nt = 1, nttop(ip))
c$$$     >            (sdcs(nt,ne)+sdcs(nt,ne+ne2e), nt = 1, nttop(ip))
               sdcst = 0.0
               do nt = 1, nttop(ip)
                  sdcst = sdcst + psdcs(nt,ne,0)+psdcs(nt,ne,1)!sdcs(nt,ne)
c$$$                  sdcst = sdcst + (sdcs(nt,ne)+sdcs(nt,ne+ne2e))
               enddo
               print'("At energy",f9.4,"eV SDCS:",
     >            1p,e10.3,a9)', xout(ne)*ry,
     >            sdcst,chun
c$$$               print'("At energy",f9.4,"eV and",f9.4,"eV  SDCS:",
c$$$     >            1p,e10.3,a9)', xout(ne)*ry, xout(ne+ne2e)*ry,
c$$$     >            sdcst,chun
c$$$               write(47,'(f6.2,1p,100e10.2)') xout(ne)*ry,sdcst,
c$$$     >            (sdcs(nt,ne), nt = 1, nttop(ip))
c$$$     >            (sdcs(nt,ne)+sdcs(nt,ne+ne2e), nt = 1, nttop(ip))
               write(48,'(f9.3,1p,100e10.3)') xout(ne)*ry,
     >            ((psdcs(nt,ne,ns),  nt = 1, nttop(ip)), ns=0,nsmax)
               write(48,'(f9.3,1p,100e10.3)') xout(ne2e+ne)*ry,
     >            ((psdcs(nt,ne2e+ne,ns), nt = 1, nttop(ip)),ns=0,nsmax)
            enddo
            close(47)
            close(48)
            if (iborn.gt.0) then
               bornsub = 1.0
            else
               bornsub = 0.0
            endif 
            call processt(ton,vdon,nchimax,nchanmax,jstart,jstop,
     >         npar,non,jlm,onshellk,knm,nopen,nze,etot,sweight(0,nn),
     >         nstep(nn),nsmax,nunit,projectile,target,hlike,
     >         vdcore,minvdc,maxvdc,bornsub,nchfi,ovlpn,phaseq,
     >         ne2e,slowery,nonnew,tfile)
            call MPI_FINALIZE(ierr)
            stop 'Job completed'
         endif 
      else
         irecl = max(100,28 * nchanmax * nchimax + 10)
         open(88,file=tfile,access='SEQUENTIAL',form='FORMATTED',
     >      recl=irecl)
         print*,'Opened POTL with IRECL:',irecl
c         open(88,file=tfile,access='SEQUENTIAL',form='FORMATTED',
c     >      recl=12000)
         if (jstart.eq.0) then
            write(88,"(a31,1p,e11.3,' eV')") projectile//' - '//target//
     >         ' scattering at',energy
            write(88,*) energy,zassymp,nopen(0),npos,jstop,jstop

            
c            do iorb=1,nopen(0)
c               write(88,'(i2,1p,2(e12.4))',ADVANCE='NO')
c     >            lea(iorb),ovlpn(iorb),onshellk(iorb)
c            end do
c            write(88,*)
            write(88,'(1p,1000(i2,2(e12.4)))')
     >         (lea(iorb),ovlpn(iorb),onshellk(iorb),
     >         iorb=1,nopen(0))
            write(88,*) (nopen(np),np=0,npar),nparity,nsp
c$$$  elseif (jstart.eq.1.and.ipar.eq.1) then
c$$$  write(88,*) nopen(np),nparity,nsp
         endif

         if (iborn.le.0) then
c$$$         if (zasym.ne.0.0.or.iborn.le.0) then
            print*,'Will set BornICS = 0.0'
            do incount = 1, nent
               ni = instate(incount)
               do nf = 1, knm
                  BornICS(nf,ni) = 0.0
               enddo
            enddo
            return
         endif 
            
C  Calculate Born form factors
         MSTEP = iborn
c$$$         print*,'Enter MSTEP (1-5) please'
c$$$         read*,MSTEP
         if (MSTEP.eq.1) then
            nnset(1) = 1
            nnset(2) = 2
            nnset(3) = 3
         elseif (MSTEP.eq.2) then
            nnset(1) = 2
            nnset(2) = 5
            nnset(3) = 7
         elseif (MSTEP.eq.3) then
            nnset(1) = 3
            nnset(2) = 5
            nnset(3) = 6
         elseif (MSTEP.eq.4) then
            nnset(1) = 4
            nnset(2) = 7
            nnset(3) = 8
         elseif (MSTEP.eq.5) then
            nnset(1) = 5
            nnset(2) = 9
            nnset(3) = 10
         endif 
         print'("Will use NSTEP:",3i3)',(nstep(nnset(n)),n=1,3)
         call clock(s1)
         jch = 0
         do incount = 1, nent
            ni = instate(incount)
            call getchinfo (ni,ntmpi,jch,psii,maxpsii,eai,lai,nai,li)
            if (ni.eq.0) cycle
            rli = lai
            posi = positron(nai,lai,nposi)
C The following has been commented out due to strange problems with two-electron targets. Waiting for nother version of compiler; 10.0.023 and lower were used.
C  The following directives are for OPENMP
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& PRIVATE(nf,ntmpf,psif,maxpsif,eaf,laf,naf,lf,rlf,posf,sf)
C$OMP& PRIVATE(const,Bornamp,mi,mf,res,res2,err1,err2,qmin,qmax,si)
            do nf = 1, nchtop
               call getchinfo (nf,ntmpf,jch,psif,maxpsif,eaf,
     >            laf,naf,lf)
               rlf = laf
               posf = positron(naf,laf,nposf)
               if (onshellk(ni).gt.0.0.and.onshellk(nf).gt.0.0) then
                  call getspinw(chan(ni)(1:1),0,nze,spinw(0),si)
                  call getspinw(chan(nf)(1:1),0,nze,spinw(0),sf)
                  const = (2.0 * pi) ** 4 * unit *
     >               onshellk(nf)/onshellk(ni) / (2.0 * lai + 1.0)
                  nn = nnset(1)
                  thfac = 0.0
                  thfac = 1.0
c$$$                  print*,'Enter THFAC'
c$$$                  read*,thfac
                  bornICS(nf,ni) = 0.0
                  if (si.eq.sf.and.posf.eqv.posi) then
                     call getBornamp(nf,laf,psif,maxpsif,onshellk(nf),
     >                  ni,lai,psii,maxpsii,onshellk(ni),hlike,nze,
     >                  nstep(nn),vdcore,minvdc,maxvdc,lexit,lentr,
     >                  thfac,Bornamp)
                     do mi = -lai, lai
                        do mf = -laf, laf
                           call intamp(Bornamp(0,mf,mi),thfac,
     >                        onshellk(nf),onshellk(ni),sweight,res,
     >                        res2,err1,err2)
                           bornICS(nf,ni) = bornICS(nf,ni) +
     >                        res * const
C$OMP critical(print)
                           print'(2X,a3,'' <-'',a3," MF, MI",2I3,
     >                        " Born ICS:",1p,e11.4,a8)', chan(nf),
     >                        chan(ni), mf, mi, res * const,
     >                        chunit(nunit)
C$OMP end critical(print)
                        enddo
                     enddo
c$$$                     pol = (sigm0 - sigm1) /
c$$$     >                  (2.375 * sigm0 + 3.749 * sigm1)
c$$$                     print '("polarization (%):",f7.3)', pol * 100.0
                     
                  endif
                  qmin = abs(onshellk(nf) - onshellk(ni))
                  qmax = onshellk(nf) + onshellk(ni)
c$$$                  print*,'testing Born:',bornICS(nf,ni)
c$$$                 
c$$$                  test(1) = 0.0
c$$$                  test(2) = 0.0
c$$$                  test(3) = 0.0
c$$$                  do nth = 0, 180
c$$$                     yin(nth+1) = 0d0
c$$$                     xin(nth+1) = nth * pi / 180.0
c$$$                     q = qmin + nth * (qmax - qmin) / 180.0 / thfac
c$$$                     xin(nth+1) = q
c$$$                     do mi = -lai, lai
c$$$                        do mf = -laf, laf
c$$$                           yin(nth+1) = yin(nth+1) +
c$$$     >                        abs(Bornamp(nth,mf,mi))**2
c$$$                        enddo
c$$$                     enddo
c$$$C  Note that sweight = 0.0 for some NTH, NN pairs
c$$$                     do n = 1, 3
c$$$                        nn = nnset(n)
c$$$                        test(n) = test(n) + sweight(nth,nn) * yin(nth+1)
c$$$     >                     * xin(nth+1) 
c$$$c$$$                        test(n) = test(n) + sweight(nth,nn) * yin(nth+1)
c$$$c$$$     >                     * sin(xin(nth+1)) * const * 2.0 * pi 
c$$$                     enddo 
c$$$                  enddo
c$$$                  do n = 1, 3
c$$$                     test(n) = test(n) * 2.0 * pi * const / thfac /
c$$$     >                  onshellk(nf) / onshellk(ni) * (qmax-qmin) / pi
c$$$                  enddo 
c$$$                  err1 = abs((test(1)-test(2))/(test(1)+1e-30))*100.0
c$$$                  err2 = abs((test(1)-test(3))/(test(1)+1e-30))*100.0
c$$$                  if (err1.gt.1.0.or.err2.gt.1.0) then
c$$$                     np = 2
c$$$                     if (ni.eq.nf) np = 0
c$$$                     do nth = 1, 181
c$$$                        yin(nth) = yin(nth) * xin(nth) ** np
c$$$                     enddo 
c$$$                     ntwo = 5
c$$$                     nthi = 180 * 2 ** ntwo + 1
c$$$                     do nth = 1, nthi
c$$$                        q = qmin + (nth-1)*(qmax - qmin)/(nthi-1)/thfac
c$$$c$$$                        xout(nth) = dble(nth - 1) / 2 ** ntwo * pi/180.0
c$$$                        xout(nth) = q
c$$$                     enddo
c$$$                     call intrpl(181,xin,yin,nthi,xout,yout(1,1))
c$$$c$$$                     open(42,file='borndcs.out')
c$$$c$$$                     do nth = 1, nthi
c$$$c$$$                        write(42,*) xout(nth),yout(nth,1)
c$$$c$$$                     enddo
c$$$c$$$                     close(42)
c$$$c$$$                     open(42,file='borndcs.in')
c$$$c$$$                     do nth = 1, 181
c$$$c$$$                        write(42,*) xin(nth),yin(nth)
c$$$c$$$                     enddo
c$$$c$$$                     close(42)
c$$$                     do n = 0, 5
c$$$                        tmp = 0.0
c$$$                        test1 = 0.0
c$$$                        h = xout(2**n+1)
c$$$                        h = (xout(2**n+1) - xout(1)) 
c$$$                        nprev = 1
c$$$                        do nth = 1, nthi, 2**n
c$$$                           weight = float(nprev) / 3.0 * h
c$$$c$$$                           contrib = weight*sin(xout(nth))*yout(nth,1)
c$$$                           contrib = weight*yout(nth,1) *
c$$$     >                        xout(nth)**(1-np)
c$$$                           tmp = tmp + contrib
c$$$                           contribt = xout(nth)**2 * weight
c$$$                           test1 = test1 + contribt
c$$$                           if (nprev.eq.4) then
c$$$                              nprev = 2
c$$$                           else
c$$$                              nprev = 4
c$$$                           endif
c$$$                        enddo
c$$$                        test(n+1) = (tmp - contrib / 2.0) * 2.0 * pi * const/thfac/
c$$$     >                     onshellk(nf) / onshellk(ni) ! * (qmax-qmin)/pi
c$$$                        test1 = test1 - contribt / 2.0
c$$$                        print*,'n,ICS, test:',n, test(n+1),
c$$$     >                     test1 * 3.0 / (qmax**3-qmin**3)
c$$$                     enddo 
c$$$                     err1 = abs((test(1)-test(2))/(test(1)+1e-30))*100.0
c$$$                     err2 = abs((test(1)-test(3))/(test(1)+1e-30))*100.0
c$$$                  endif
c$$$c$$$                  do n = 1, 3
c$$$c$$$                     nn = nnset(n)
c$$$c$$$                     test(n) = 0.0
c$$$c$$$                        finalcont(n) = 0.0
c$$$c$$$                        do mi = -lai, lai
c$$$c$$$                           do mf = -laf, laf
c$$$c$$$                              do nth = 0, 180, nstep(nn)
c$$$c$$$                                 thrad = nth * pi / 180.0 
c$$$c$$$                              fac=abs(Bornamp(nth,mf,mi))**2*sin(thrad)
c$$$c$$$c$$$                                 fac=abs(Bornamp(nth,mf,mi))**2 * q
c$$$c$$$                                 test(n) = test(n)+ sweight(nth,nn)*fac
c$$$c$$$                              enddo
c$$$c$$$                           enddo
c$$$c$$$                        enddo
c$$$c$$$                        finalcont(n) = finalcont(n) + const *2.0 * pi * 
c$$$c$$$     >                     fac * sweight(180,nn)
c$$$c$$$c$$$     >                     / thfac /
c$$$c$$$c$$$     >                     onshellk(nf) / onshellk(ni) * (qmax-qmin)/pi
c$$$c$$$                        test(n) = test(n) * 2.0 * pi * const
c$$$c$$$c$$$     >                     / thfac /
c$$$c$$$c$$$     >                     onshellk(nf) / onshellk(ni) * (qmax-qmin)/pi
c$$$c$$$                     enddo
c$$$c$$$                     err1 = abs((test(1)-test(2))/(test(1)+1e-30))*100.0
c$$$c$$$                     err2 = abs((test(1)-test(3))/(test(1)+1e-30))*100.0
c$$$c$$$                     if (thfac.ne.1.0) print'("CAUTION: thfac:",f5.1,
c$$$c$$$     >                  " for qratio",1p,1e9.1,". Res and last contr.:",
c$$$c$$$     >                  2e9.1)',thfac, qratio, test(1), finalcont(1)
c$$$c$$$                     print*,test(1), test(2), test(3)
c$$$                  bornICS(nf,ni) = test(1)
C$OMP critical(print)
                  print'(a3,'' <-'',a3," Born ICS:",1p,e11.4,a8,0p,
     >               2(f6.2,"%")," q(0),q(180):",2f7.3)', chan(nf),  
     >               chan(ni),BornICS(nf,ni),chunit(nunit),err1,err2,
     >               qmin,qmax
C$OMP end critical(print)
               endif            ! end of check of open channel
            enddo               ! end nf loop
C$OMP END PARALLEL DO
         enddo                  ! end ni loop
         call clock(s2)
         print '("Time (seconds) to calculate Born ICS:",f7.1)',s2-s1
      endif                     ! end of if exists condition
      return
 25   close(88)
      print*,'Delete potl and start again  '//tfile
      stop 'Delete potl and start again '
 21   close(88)
      print*,'Clean entry in potl and start again. J, NS, Par',j,NS,np,
     >   tfile
      stop 'Clean entry in potl and start again '
 22   close(88)
      print*,'Error reading NON line. J, NS, Par',j,NS,np,
     >   tfile
      stop 'Error reading NON line'
 23   close(88)
      print*,'Error reading VDON line. J, NS, Par',j,NS,np,
     >   tfile
      stop 'Error reading VDON line'
 24   close(88)
      print*,'Error reading TON line. J, NS, Par',j,NS,np,
     >   tfile
      stop 'Error reading TON line'

      end

      subroutine getspinw(ch,ns,nze,spinw,s)
      character ch
      if (ns.ne.0.and.ns.ne.1) stop 'GETSPINW assumes 0 or 1 for NS'
      if (ch.eq.'t'.or.ch.eq.'T') then
C  Incident on a triplet state of Helium-like target. Total spin S = ns + 0.5.
         s = 1.0
         Stot = float(ns) + 0.5
      elseif (ch.eq.'s'.or.ch.eq.'S') then
C  Incident on a singlet state of Helium-like target. Total spin S = 0.5.
C  For ns = 1 set Stot = -0.5 so that the spin weight comes out as zero
         s = 0.0
         Stot = 0.5
         if (ns.eq.1) Stot = - 0.5
      else
C  Incident on a hydrogenic target. Total spin S = ns.
         s = 0.5
         Stot = float(ns)
      endif
                  
      if (nze.eq.1) then
         spinw = 1.0
      else
         spinw = (2.0*Stot+1.0)/2.0/(2.0*s+1.0)
      endif
      return
      end
               
      
      subroutine processt(ton,vdon,nchimax,nchanmax,jstart,jstop,
     >   npar,non,jlm,onshellk,nomax,nopen,nze,etot,sweight,nstep,nsmax,
     >   nunit,projectile,target,hlike,vdcore,minvdc,maxvdc,bornsubin,
     >   nchfi,ovlpn,phaseq,ne2e,slowery,nonnew,tfile)
      include 'par.f'
      parameter (lentr=8,lexit=lamax,npoints=30,ntype=20)
      character projectile*(*), target*(*), tfile*(*)
      integer nchfi(nomax,0:1)
      dimension onshellk(nchan), temp(maxr), partcs(0:1), partbcs(0:1),
     >   ovlp(20,100), chi(maxr), sweight(0:200), spinw(0:1), units(3),
     >   psii(maxr), psif(maxr), vdcore(maxr,0:lamax),slowery(ne2e),
     >   ovlpn(knm),dcs(0:200,knm),sign(0:1),
     >   sigm(-lexit:lexit,-lentr:lentr),
     >   ticsm(-lentr:lentr),res2(5),csmt(5,knm),crossmt(5),sum(knm)
      logical hlike
      integer*4 jlm(nomax,0:1,0:lmax), non(nomax,2*lamax+1,0:1,0:lmax),
     >   nonnew(nomax,-lamax:lamax,0:1,0:lmax)
      complex ton(nchanmax,nchimax,0:nsmax,0:npar,0:jstop), cfac(0:1),
     >   vdon(nchanmax,nchimax,0:nsmax,0:npar,0:jstop), phase,
     >   f(0:200,0:1,-lexit:lexit,-lentr:lentr),ampt(0:200,0:1),cip,
     >   fB(0:200,-lexit:lexit,-lentr:lentr),phaseq(knm),ovlpq,ci,x,
     >   t, getT, getTextr, Tel(0:1), Vel(0:1),
     >   fbsum(0:200,-lexit:lexit,-lentr:lentr)
      real*8 rylm,xin(npoints),yin(npoints),xout(npoints),yout(npoints)
      character chan(knm)*3
      common /charchan/ chan
      common /chanen/ enchan(knm)
      common/meshrr/ meshr,rmesh(maxr,3)
      character ch*1, xunit(3)*8
      integer satom
      common /helium/ latom(KNM), satom(KNM), lpar(KNM), np(KNM)
      asym(sing,trip,fac) = (sing - trip / fac) / (sing + trip + 1e-30)
      ch(i)=char(i+ichar('0'))
      getT(j,js,t) = cmplx(real(t) * (2*js-1) * (2*js+1) * (2*js+3) /
     >   (2*j-1) / (2*j+1) / (2*j+3), aimag(t) * (2*js-1) * (2*js+1) *
     >   (2*js+3) * (2*js+5) * (2*js+7) /
     >   (2*j-1) / (2*j+1) / (2*j+3) / (2*j+5) / (2*j+7))
      DATA XUNIT/'A0**2   ','pi*A0**2','cm**2   '/
      data units/1.0,0.3183099,28.002e-18/
      data pi/3.1415927/
      data ci/(0.0,1.0)/
      

      if (nunit.gt.0) unit = units(nunit)
      fac = 2.0
      if (chan(1)(1:1).eq.' ') fac = 3.0
      if (hlike) then
         ich = 2
      else
         ich = 1
      endif 
      crossint = 0.0
      crossmt(:) = 0.0
      print*,
     >   'Enter max J for extrapolation' ! and bornsub (1 or 0)'
c$$$      read(*,err=100), jextrap
      read(*,*), jextrap
      bornsub = bornsubin
      print*, 'JEXTRAP, BORNSUBIN:',jextrap, bornsubin
 
 37   print*,'Please enter NISTART, NISTOP, NFSTART and NFSTOP:'
c$$$      read(*,err=100), NISTART, NISTOP, NFSTART, NFSTOP
      read(*,*), NISTART, NISTOP, NFSTART, NFSTOP
      print*,'NISTART, NISTOP, NFSTART, NFSTOP:',
     >   NISTART, NISTOP, NFSTART, NFSTOP
      if (nistart.le.0) then
c$$$         call MPI_FINALIZE(ierr)
         stop 'amplitude generation completed'
      endif 
      do ni = nistart, nistop
         if (onshellk(ni).le.0.0) cycle
         ndcs = 0
         ticsm(:) = 0.0
         do 10 nf = nfstart, nfstop
            ndcs = ndcs + 1
            if (ndcs.gt.knm) stop 'increase dimension of DCS'
            bornsub = bornsubin
C  Turn off Born subtraction for elastic channels is now done with extrapolation below..
C  At low energies forward amplitudes can be incorrect.
C  Also, we have a problem in that
C  analytical Born for P-states doesn't give zero amp for MF = 1, MI = -1.
C  Also turn off Born subtraction for atom-positronium pairs
            if ((chan(ni)(1:1).ne.chan(nf)(1:1).and.
     >         chan(nf)(1:1).eq.'p')) then
c$$$     >         chan(nf)(1:1).eq.'p').or.ni.eq.nf) then
               print*,'No Born subtraction for this channel pair'
               bornsub = 0.0
            endif 
            do nth = 0, 200
               dcs(nth,ndcs) = 0.0
            enddo
            rkf = onshellk(nf)
            if (onshellk(nf).le.0.0) then
               print*,'One of NI or NF channels is closed',ni,nf
               go to 10
            endif 
            do ns = 0, nsmax
               call getspinw(chan(ni)(1:1),ns,nze,spinw(ns),si)
               partcs(ns) = 0.0
               partbcs(ns) = 0.0
            enddo
            partcs(1) = 0.0
            partbcs(1) = 0.0
            
            do mi = - lentr, lentr
               do mf = - lexit, lexit
                  do ns = 0, nsmax
                     do nth = 0, 200
                        f(nth,ns,mf,mi) = (0.0,0.0)
                        fbsum(nth,mf,mi) = (0.0,0.0)
                     enddo 
                  enddo
               enddo
            enddo

            jch = 0
            c1 = onshellk(nf)/onshellk(ni) * (2.0 * pi) ** 4 * unit 
            call getchinfo (ni,nchip,jch,psii,maxpsii,ei,lia,nia,li)
            if (lia.gt.lentr) stop 'increase LENTR in PROCESST'
            rlia = lia
            liapar = (-1)**lia
            if (lpar(ni).ne.0) liapar = lpar(ni)
            
            call getchinfo (nf,nchp,jch,psif,maxpsif,ef,lfa,nfa,lf)
            if (lfa.gt.lexit) stop 'increase LEXIT in PROCESST'
            rlfa = lfa
            lfapar = (-1)**lfa
            if (lpar(nf).ne.0) lfapar = lpar(nf)
            
            if (ef.lt.0.0) then
               ovlpq = (1.0,0.0)
               nt = -1
            else
C An extra sqrt below is so that the TDCS program didn't have to worry about 
C Kb. It comes from the division by the Kb of the Coulomb wave combined with
C the multiplication by Kb in the definition of the TDCS resulting in
C the factor of 1/Kb multiplying the TDCS. Note OVLPQ will be squared for TDCS.
               ovlpq = phaseq(nf) * ovlpn(nf) / sqrt(sqrt(ef))
               print*, 'phase, ovlp:', phaseq(nf), ovlpn(nf)
               if (nunit.eq.3) ovlpq = ovlpq / sqrt(27.2116)
            endif
c$$$            jextrap = jstop
            rki = onshellk(ni)
            if (ni.eq.nf) then
c$$$               print*,'Elastic extrapolation is set to J = 80'
c$$$               jextrap = 80
               alpha = real(ton(1,1,0,0,jstop)) * (2*jstop-1) *
     >            (2*jstop+1) * (2*jstop+3)
               beta = aimag(ton(1,1,0,0,jstop)) * (2*jstop-1) *
     >            (2*jstop+1) * (2*jstop+3) * (2*jstop+5) * (2*jstop+7)
               print*,'ALPHA, BETA:',alpha/rki, beta/rki
               do j = jstop+1, jextrap
                  c2 = sqrt((2.0*j+1.0)) / sqrt(pi * 4.0)
                  treal = alpha / (2*j-1) / (2*j+1) / (2*j+3)
                  timag = beta / (2*j-1) / (2*j+1) / (2*j+3) / (2*j+5) /
     >               (2*j+7)
                  do ns = 0, nsmax
                     cfac(ns) = c2*cmplx(treal,timag)
                  enddo
                  do nth = 0, 200, nstep
                     theta = float(nth)
                     if (nth.gt.180) theta = (nth - 180) / 21.0
                     thrad = theta * pi / 180.0
                     sh=rylm(j,0,dble(thrad))
                     do ns = 0, nsmax
                        f(nth,ns,0,0) = sh * cfac(ns)+ f(nth,ns,0,0) 
                     enddo
                  enddo
               enddo
            endif 
c$$$            open (77,file='klaus.25ev.bsr')
c$$$            read(77,*)
            do j = 0, jextrap
               do ipar = 0, min(j,npar)
               rj = j
               iparity = (-1)**(j + ipar)
               jlmi = -1 + ipar
               do li = abs(j - lia), j + lia
                  if ((-1)**li*liapar.eq.iparity) then
                     jlmi = jlmi + 2 
                     if (j.le.jstop) then
                        nchiold = non(nchfi(ni,ipar),jlmi,ipar,j)
                        nchi = nonnew(nchfi(ni,ipar),j-li,ipar,j)
                        if (nchi.ne.nchiold) print*,
     >                     'WARNING: NCHI .ne. NCHIOLD'
                        if (nchi.eq.0) then
                           print*,'WARNING: NCHI = 0'
                           goto 10
                        endif
                     endif 
                     rli  = li
                     jlmf = -1 + ipar
                     if (lfapar.ne.(-1)**lfa.and.j.gt.0) jlmf = -ipar
                     do lf = abs(j - lfa), j + lfa
                        if ((-1)**lf*lfapar.eq.iparity) then
                           jlmf = jlmf + 2
                           if (j.le.jstop) then
                              nchf = non(nchfi(nf,ipar),jlmf,ipar,j)
                              if (nchf.ne.
     >                           nonnew(nchfi(nf,ipar),j-lf,ipar,j))
     >                           then
                                 print*,'Resetting NON to',
     >                              nonnew(nchfi(nf,ipar),j-lf,ipar,j)
                                 nchf =
     >                              nonnew(nchfi(nf,ipar),j-lf,ipar,j)
                              endif 
                              if (nchf.eq.0) then
                                 print*,'NCHF,NCHFI(nf,ip),jlmf,ipar,j',
     >                              NCHF,NCHFI(nf,ipar),jlmf,ipar,j
                                 print*,'WARNING: NCHF = 0'
                                 goto 10
                              endif
                              do ns = 0, nsmax
                                 Tel(ns) = ton(nchf,nchi,ns,ipar,j)
                                 Vel(ns) = vdon(nchf,nchi,ns,ipar,j)
c$$$                                 read(77,*) xx,kj,ks,kp,kli,klf,trk,tik
c$$$                                 write(*,'(5i5,1p,4e15.5)')
c$$$     >                              j,ipar,li,lf,ns,
c$$$     >                              -2.0*pi*ci*sqrt(rki*rkf)*Tel(ns),
c$$$     >                              trk,tik
c$$$                                 if (j.gt.11)
c$$$     >                              Tel(ns) =ci*cmplx(trk,tik)/
c$$$     >                              (2.0*pi*sqrt(rki*rkf))
                              enddo 
                           else
                              nchf = non(nchfi(nf,ipar),jlmf,ipar,jstop)
                              nchi = non(nchfi(ni,ipar),jlmi,ipar,jstop)
                              nchi =
     >                           nonnew(nchfi(ni,ipar),j-li,ipar,jstop)
                              nchf =
     >                           nonnew(nchfi(nf,ipar),j-lf,ipar,jstop)
                              Tel(:) = (0.0,0.0)
                              Vel(:) = 0.0
                              do ns = 0, nsmax
                                 if(spinw(ns) .eq. 0.0) cycle
!                                 Tel(ns) = getT(j,jstop,
!     >                              ton(nchf,nchi,ns,ipar,jstop))
!                                 Vel(ns) = getT(j,jstop,
!     >                              vdon(nchf,nchi,ns,ipar,jstop))
                                 Tel(ns) = getTextr(j,jstop,
     >                              ton(nchf,nchi,ns,ipar,jstop),
     >                                ton(nchf,nchi,ns,ipar,jstop-1))
                                 Vel(ns) = getTextr(j,jstop,
     >                              vdon(nchf,nchi,ns,ipar,jstop),
     >                              vdon(nchf,nchi,ns,ipar,jstop-1))
                              enddo 
                              if (j.eq.jstop+1) then
                                 t = getT(0,jstop,ton(nchf,nchi,0,
     >                              ipar,jstop))
                                 print*,'ALPHA, BETA:',-real(t)/rki
     >                              *3.0,-aimag(t)/rki*105.0
                                 do jpr = jstop - 3, jstop
                                    print*, jpr, getT(jpr,jstop,
     >                                 ton(nchf,nchi,0,ipar,jstop)),
     >                                 ton(nchf,nchi,0,ipar,jpr)
                                 enddo
                              endif 
                           endif 
                           rlf  = lf
                           cip = (0d0,1d0)**dfloat(li-lf)
                           do ns = 0, nsmax
                              if(spinw(ns) .eq. 0.0) cycle
                              const = spinw(ns) * c1 / (4.0 * pi) *
     >                           (2.0*j + 1.0) / (2.0 * rlia + 1.0)
                              partcs(ns) = partcs(ns) + const *
     >                           Tel(ns)*CONJG(Tel(ns))
!     >                           abs(Tel(ns))**2
                              partbcs(ns) = partbcs(ns) + const *
     >                           abs(Vel(ns))**2
                           enddo
                           c2 = sqrt((2.0*rli+1.0)) / sqrt(pi * 4.0)
                           do mi = - lia, lia
                              rmi = mi
                              c3=cgc(rli,0.0,rlia,rmi,rj,rmi)*c2
                              do mf = - lfa, lfa
                                 rmf = mf
                                 c4=cgc(rlf,rmi-rmf,rlfa,rmf,rj,rmi)*c3
                                 if (abs(c4).gt.1e-10) then
                                    x = c4 * cip * Vel(0)
                                    cfac(:) = (0.0,0.0)
                                    do ns = 0, nsmax
                                       if(spinw(ns) .eq. 0.0) cycle
                                       cfac(ns) = c4 * cip *
     >                                    (Tel(ns) - bornsub * Vel(ns))
                                    enddo
                                    do nth = 0, 200, nstep
                                       theta = float(nth)
                                       if (nth.gt.180) theta =
     >                                    (nth - 180) / 21.0
                                       thrad = theta * pi / 180.0
                                       sh=rylm(lf,mi-mf,dble(thrad))
                                       do ns = 0, nsmax
                                          f(nth,ns,mf,mi) = sh *
     >                                       cfac(ns) + f(nth,ns,mf,mi) 
                                       enddo
                                       fbsum(nth,mf,mi) = sh * x
     >                                    + fbsum(nth,mf,mi)
                                    enddo
                                 endif 
                              enddo ! end of mf loop
                           enddo !end of mi loop
                        endif   ! end of parity check for LF + LFA
                     enddo      ! end of lf loop
                  endif         ! end of parity check for LI + LIA
               enddo            ! end of li loop
            end do              ! end of parity loop
            enddo               ! end of J loop
            nth = 5
            x = (0.0,0.0)
            mff=0
            mii=0
            do mf = - lfa, lfa
               do mi = - lia, lia
c$$$                  print*,'mf,mi,born amp',mf,mi,fbsum(nth,mf,mi)
                  if (abs(fbsum(nth,mf,mi)).gt.abs(x)) then
                     x = fbsum(nth,mf,mi)
                     mff = mf
                     mii = mi
                  endif
               enddo
            enddo 
            thfac = 0.0
            fB(:,:,:) = (0.0,0.0)
            if (bornsub.ne.0.0) call getBornamp(nf,lfa,psif,maxpsif,
     >         onshellk(nf),ni,lia,psii,maxpsii,onshellk(ni),hlike,
     >         nze,nstep,vdcore,minvdc,maxvdc,lexit,lentr,thfac,fB)
c$$$            do mf = - lfa, lfa
c$$$               do mi = - lia, lia
c$$$                  print*,'mf,mi,born amp',mf,mi,fb(nth,mf,mi)
c$$$               enddo
c$$$            enddo 
c$$$            print*,'Born from partials and analytic',x,fb(nth,mff,mii)
            x = x / (fb(nth,mff,mii)+1e-30)
c$$$            print*,'enter sign'
c$$$            read*,x
            sign(0) = bornsub * real(x) / (abs(real(x)+1e-30))
            if (sign(0).lt.0.0) then
               print*,'CAUTION negative sign!, reset to plus one'
               sign(0) = 1.0
c$$$               print*,'Please enter new value (+/-1)'
c$$$               read*, sign(0)
            endif 
            
            sign(1) = sign(0)
C  The following line ensures that for 2-electron atoms there is no Born
C  subtraction for total spin = 3/2 if the transition involves singlet states.
            call getspinw(chan(nf)(1:1),0,nze,sw,sf)
            if (sf*si.lt.1e-5) sign(1) = 0.0
            if (bornsub.eq.0.0) then
               print*,'No Born subtraction will be used'
            else 
               print*,
     >            'X (should be complex +/-1), sign(0), analytic ',
     >            'and summed Born, mff,mii= , and matching NTH:',
     >            x, sign(0), fb(nth,mff,mii), fbsum(nth,mff,mii),
     >            mff, mii, nth 
            endif 
            y = sqrt(c1/(2.0 * rlia + 1.0)) 
            
            open(42,file='amp.'//chan(nf)(ich:3)//'-'//chan(ni)(ich:3))

            write(42,"(a8,' - ',a6,'scattering:',1p,e10.3,'eV on ',a3,
c$$$     >         ' -> ',a3,(2f10.5,'eV')/)") projectile,target,
     >         ' ->',a3,(2(e10.3,'eV'))/)") projectile,target,
     >         13.6058*(etot-enchan(ni)),chan(ni),chan(nf),
     >         13.6058 * enchan(nf), 13.6058 * (etot-enchan(nf))

            write(42,'(''When the amplitudes are squared the result '',
     >         ''contains the factor Kf/Ki/(2Li+1) ''/
     >         ''and is in units of: '',a8,
     >         ''The last amplitude is Born (analytic if calculated).''
     >         )') xunit(nunit)

            write(42,'(
     >         ''The CCC amplitudes may be read in the following way:''/
     >         ''      read(n,*) lf, li, nth, rkf, rki, tsp1, nsm,'',
     >         '' (spinw(ns), ns=0, nsm)''/
     >         ''      do mi = - li, li''/
     >         ''         do mf = - lf, lf''/
     >         ''            read(n,*) mfprev, miprev''/
     >         ''            if (mfprev.ne.mf.or.miprev.ne.mi) stop'',
     >         '' "incorrect format" '')')
            write(42,'(
     >         ''            do ith = 1, nth''/
     >         ''               read(n,*) theta(ith), (freal(ns),'',
     >         '' fimag(ns), ns = 0, nsm)''/
     >         ''               f(ith,mf,mi,0) = cmplx(freal(0),'',
     >         '' fimag(0))'')')
            write(42,'(
     >         ''               f(ith,mf,mi,1) = cmplx(freal(1),'',
     >         '' fimag(1))''/
     >         ''            end do''/
     >         ''         end do''/
     >         ''      end do'')')
            write(42,'(2i2,i4,2f8.5,2i2,2f8.5,
     >         " Lf,Li,Nth,Kf,Ki,2Sf+1,Nsm,SpinW(ns)")')
     >         lfa, lia, 181/1,  onshellk(nf), onshellk(ni),
     >         nint(2.0*sf+1.0),nsmax,(spinw(ns),ns=0,1)
            bornint = 0.0
            csmt(:,nf) = 0.0
            sum(nf) = 0.0
            do mli = -lia, lia
               do mlf = -lfa, lfa
                  sigm(mlf,mli) = 0.0
                  do ns = 0, nsmax
                     do nth = 0, 200, nstep
c$$$                        if (nth.lt.10) print*,
c$$$     >                     nth,f(nth,ns,mlf,mli),fb(nth,mlf,mli)
                        ampt(nth,ns)=y*(f(nth,ns,mlf,mli)+
     >                     sign(ns)*fB(nth,mlf,mli))
                        dcs(nth,ndcs) = dcs(nth,ndcs) + spinw(ns) *
     >                     abs(ampt(nth,ns)) ** 2
                     enddo
C  The ,0.0, below indicates that ampt is on a linear theta grid
                     call intamp(ampt(0,ns),0.0,onshellk(nf),
     >                  onshellk(ni),sweight,res,res2,err1,err2)
                     
c$$$                     call integrate(ampt(0,ns),sweight,nstep,res,res2)
                     sigm(mlf,mli) = sigm(mlf,mli) + spinw(ns) * res
                     csmt(:,nf) = csmt(:,nf) + spinw(ns) * res2(:)
                     sum(nf) = sum(nf) + spinw(ns) * res
                  enddo 
                  call intamp(fB(0,mlf,mli),0.0,onshellk(nf),
     >               onshellk(ni),sweight,res,res2,err1,err2)
c$$$                  call integrate(fB(0,mlf,mli),sweight,nstep,res,res2)
                  bornint = bornint + y ** 2 * res * bornsub
                  write(42,'(2i5,1p,e15.6,"  MF, MI, ICS")') 
     >               mlf,mli,sigm(mlf,mli)
                  print'(a3,'' <-'',a3,1p,e10.3,2i3,
     >               " magnetic sublevel integrated DCS, Mf, Mi")',
     >               chan(nf),chan(ni),sigm(mlf,mli),mlf,mli
                  if (enchan(nf).gt.0.0) ticsm(mli) =
     >                ticsm(mli) + sigm(mlf,mli)
                  if (enchan(nf).lt.0.0) then
                     do nth = 0, 180, nstep
                        write(42,'(f5.1,1p,4(1x,2e11.3))') float(nth),
     >                     (ampt(nth,ns) * ovlpq, ns = 0,nsmax),
     >                     sign(0) * y * fB(nth,mlf,mli) * ovlpq,
     >                     sign(0) * y * fBsum(nth,mlf,mli) * ovlpq
                     enddo
                  else if (abs((etot-2.0*enchan(nf))/etot).lt.0.01) then
                     print*,'AMP suitable for equal-energy (e,2e)'
                     do nth = 0, 180, nstep
                        write(42,'(f5.1,1p,4(1x,2e11.3))') float(nth),
     >                     (ampt(nth,ns) * ovlpq, ns = 0,nsmax),
     >                     (ampt(nth,ns) * ovlpq, ns = 0,nsmax)
                     enddo
                  else
                     print*,'AMP suitable for asymmetric-energy (e,2e)'
                     do nth = 0, 180, nstep
                        write(42,'(f5.1,1p,4(1x,2e11.3))') float(nth),
     >                     (ampt(nth,ns) * ovlpq, ns = 0,nsmax),
     >                     (ampt(nth,ns) * 0.0, ns = 0,nsmax)
                     enddo
                  endif 
               enddo
            enddo 
            close(42)
            if (target .eq. 'H  I') then
               if (lia.eq.0.and.lfa.eq.1) then
                  pol = (sigm(0,0) - sigm(1,0)) /
     >               (2.375 * sigm(0,0) + 3.749 * sigm(1,0))
                  print '(a3," polarization (%):",f7.3)',
     >               chan(nf), pol * 100.0
               elseif (lia.eq.0.and.lfa.eq.2) then
                  pol = 3.0 * (-2.0*sigm(2,0) + sigm(1,0) + sigm(0,0))/
     >               (6.0*sigm(2,0) + 9.0*sigm(1,0) + 5.0*sigm(0,0))
                  print '(a3," polarization (%):",f7.3)',
     >               chan(nf), pol * 100.0
               endif 
            endif
            print*,'written file: amp.'//chan(nf)(ich:3)//'-'//
     >         chan(ni)(ich:3)

            summedcs = partcs(0) + partcs(1)
            asymcs = asym(partcs(0),partcs(1),fac)
            print'(a3,'' <-'',a3,1p,2e10.3,
     >         " Summed and integrated Born CS")',
     >         chan(nf),chan(ni),partbcs(0) + partbcs(1), bornint
            bdiff = bornint - partbcs(0) - partbcs(1)
            print'(a3,'' <-'',a3,1p,2e10.3,0p,f8.4,
     >         " Summed CS, integrated CS, and spin-asym")',
     >         chan(nf),chan(ni), summedcs, sum(nf), asymcs
            tdiff = sum(nf) - summedcs
            if (abs(bdiff/(bornint+1e-30)).gt.0.01.and.abs(
     >         (tdiff-bdiff)/(tdiff+1e-30)).gt.0.1) print*,
     >         'CAUTION:  Borndiff and Tdiff are different',bdiff,
     >         tdiff
            print*
c$$$            if (ni.eq.nf) then
c$$$               crossint = sum
c$$$               crossmt(:) = csmt(:)
c$$$            endif 
 10      continue               ! end nf loop
         i = 1
         do while (tfile(i:i).ne.' ')
            i = i + 1
         enddo
         open(42,file=tfile(1:i-1)//'_dcs'//ch(ni))
         write(42,"('# ',f8.3,' eV ',a8,' - ',a6,a3,
     >      '   scattering DCS in units of ',a8)")
     >      13.6058*(etot-enchan(ni)),projectile,target,chan(ni),
     >      xunit(nunit)
         do n = 1, ndcs
            if (csmt(1,n+nfstart-1).gt.0.0)
     >         write(42,"('#',a4,' ICS:',1p,e10.3,'  MTCS(l):',5e11.3)") 
     >         chan(n+nfstart-1),sum(n+nfstart-1),
     >         (csmt(l,n+nfstart-1),l=1,5)
         enddo 
         write(42,"('#ICS',1p,1000e10.3)")(sum(n+nfstart-1),n=1,ndcs)
         write(42,"('#angle',1000(a4,6x))")(chan(n+nfstart-1),n=1,ndcs)
         
         do nth = 0, 180, nstep
            write(42,'(i4,1p,1000e10.3)') nth,(dcs(nth,n),n=1,ndcs)
         enddo
         close(42)
         do mli = -lia, lia
            print*,'chan, m, TICS(m):', chan(ni), mli, ticsm(mli)
         enddo 
      end do                    ! end ni loop
      go to 37
 100  stop 'incorrect input to amplitude generation'
      return
      end
     
C  The following routine integrates the complex amplitude over the theta
C  and phi angles
      subroutine integrate(amp,w,nstep,res1,res2)
      dimension amp(0:200),w(0:200,3)
      complex amp
      data pi/3.1415927/

      nn = 1
      res1 = 0.0
      res2 = 0.0
      do i = nstep, 180, nstep
         res1 = res1 + abs(amp(i)) ** 2 * sin(i * pi / 180.0) * w(i,nn)
         res2 = res2 + abs(amp(i)) ** 2 * sin(i * pi / 180.0) * w(i,nn)
     >      * (1.0 - cos(i * pi / 180.0))
      enddo
      res1 = res1 * 2.0 * pi
      res2 = res2 * 2.0 * pi
      end

      subroutine getBornamp(nf,lf,psif,maxpsif,rkf,ni,li,
     >   psii,maxpsii,rki,hlike,nze,nstep,vdcore,minvdc,maxvdc,
     >   lexit,lentr,thfac,Bornamp)
      use CI_MODULE 
      include 'par.f'
      complex Bornamp(0:200,-lexit:lexit,-lentr:lentr)
      logical hlike
      dimension psif(maxr), psii(maxr), vdcore(maxr,0:lamax), 
     >   temp(maxr), ucentr(maxr), q(0:200), qth(0:200), vc(0:200),
     >   ff(0:200,0:lamax)
      common/meshrr/ meshr,rmesh(maxr,3)
      real*8 rylm,dtheta,thrad,pi,dki,dkf,tmp
      common /double/njdouble,jdouble(22)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common/smallr/ formcut,regcut,expcut,fast,match
      logical fast,match,allocated
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      data ucentr/maxr*0.0/
      common /noblgegas_switch/ i_ng ! set in  cc/main.f
      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
c      common /dynamical_C/ Nmax, namax,pnewC
c      integer pnewC


      istat = 0
      if (.not.allocated(C)) allocate (C(1,1,1),stat=istat)
      if (istat.ne.0) then
         print*,'Allocation problems, Getbornamp ISTAT:', istat
c$$$         stop 'Allocation problems for C in Getbornamp'
      endif 
      
      Nmax = SIZE(C,1)
      namax = SIZE(C,2)

      pi = acos(-1d0)
      dki = rki
      dkf = rkf
      zeff = 0.0
      sqrt4pi = sqrt(4.0 * pi)
      do i = 1, maxr
         temp(i) = 0.0
      enddo
      if (lf.gt.lexit.or.li.gt.lentr.or.li+lf.gt.lamax) then
         print*,lf,lexit,li,lentr,li+lf,lamax
         stop 'lf.gt.lexit.or.li.gt.lentr.or.li+lf.gt.lamax'
      endif
      rlf = lf
      rli = li
      if (hlike) then
         ze = 1.0
         maxtemp = min(maxpsif,maxpsii)
         do i = 1, maxtemp
            temp(i) = psif(i) * psii(i) 
         enddo
      else
         ze = 2.0
         if(i_ng .eq. 1) then
            ze = 2*(2*l_ng+1)
         endif
         maxtemp = 0
         do n = 1, nspm
            maxtemp = max(maxtemp,maxf(n))
         enddo
      endif 
      do mi = - li, li
         do mf = - lf, lf
            do nth = 0, 200, nstep
               Bornamp(nth,mf,mi) = (0.0,0.0)
            enddo 
         enddo
      enddo

      qmin = abs(dkf - dki)
      qmax = dkf + dki

      do nth = 0, 200, nstep
         dtheta = dfloat(nth)
         if (nth.gt.180) dtheta = (nth - 180) / 21d0
         thrad = dtheta * pi / 180d0 
C  Q = Ki - Kf  (as vectors)               
         qsq = dkf**2 + dki**2 - 2d0 * dkf * dki * cos(thrad)
         qsq = dkf * dki * (dble(dkf/dki) + dble(dki/dkf) -
     >      dble(2d0 * cos(thrad)))
c$$$         qsq = (rkf - rki) ** 2 - 2.0 * rkf * rki * (cos(thrad) - 1d0)
         q(nth) = sqrt(abs(qsq))
c$$$         print*, dtheta, q(nth)
c$$$               tmp = -(rkf**2-q(nth)**2-rki**2)/
c$$$     >            2.0/q(nth)/rki

c$$$         if (dkf.eq.dki.and.nth.eq.0) print*,'Q at zero degrees:',q(nth)
         if (q(nth).eq.0.0) then
c$$$            q(nth) = 1e-10
            tmp = 0d0
C  The following ensures that for elastic P-P transitions we get 0 amplitude
C  for MI <> MF at theta = 0,  which leads to a step, unlike if tmp = 0 above.
c$$$            tmp = 1.0
         else 
            tmp = (dki-dkf*cos(thrad)) / q(nth)
         endif 

C The following is used to calculate the Born amp on q grid
         if (thfac.ne.0.0.and.qmin.gt.0.0.and.nth.le.180) then
c$$$            q(nth) = qmin + nth * (qmax - qmin) / 180.0 / thfac + 1e-10

C  Try logarithmic q-grid
            q(nth) = exp(
     >         log(qmin) + nth * (log(qmax) - log(qmin))/ 180.0)

            tmp = (dki ** 2 - dkf ** 2 + q(nth) ** 2)/(2d0*dki * q(nth))
            if (ni.eq.nf.and.nth.eq.0) tmp = 0d0
         endif 

         if (abs(tmp).gt.1.0) tmp = tmp/abs(tmp)         
C  Qth is the angle of the Q vector to the Ki (Z) axis                
         qth(nth) = acos(tmp)
         vc(nth) = 0.0
         do i = minvdc, maxvdc
            vc(nth) = vc(nth) + sin(q(nth)*rmesh(i,1)) *
     >         vdcore(i,li) * rmesh(i,1) * rmesh(i,3)
         enddo
         vc(nth) = - nze * vc(nth) * q(nth) / ze ! -nze here to cancel the cgc one below
      enddo

      nth = 0
      if (q(nth).le.0.0) then  !1e-5
C  In the following we have taken limit of VC/q**2 for q = 0
C  Note that VC(0) will not be continuous with others.
         vc(nth) = 0.0
         do i = minvdc, maxvdc
            vc(nth) = vc(nth) + rmesh(i,1) * rmesh(i,1) *
     >         vdcore(i,li) * rmesh(i,3)
         enddo
         vc(nth) = - nze * vc(nth) / ze ! -nze here to cancel the cgc one below
      endif



      do lam = abs(lf - li), lf + li
         call formfac(hlike,nf,ni,rlf,rli,temp,maxtemp,lam,nstep,
     >      q,vc,ff(0,lam),Nmax,namax,C,ze)
         rlam = lam
         do mai = -li, li
            rmi = mai
            do maf = -lf, lf
               rmf = maf
               mu = maf - mai
               rmu = mu
               cgcoeff = - nze * cgc(rli,rmi,rlam,rmu,rlf,rmf) 
               if (cgcoeff.ne.0.0) then
                  do nth = 0, 200, nstep
                     dtheta = dble(qth(nth))
                     srylm = sqrt4pi / sqrt(2.0 * rlf + 1.0) *
     >                  rylm(lam,mu,dtheta)*(-1)**mu
c$$$                     if (nth.eq.0.and.mu.ne.0) print*,
c$$$     >                  'dtheta,lam,mai,maf,ff(nth,lam),SRYLM:',
c$$$     >                  dtheta,lam,mai,maf,ff(nth,lam),srylm
                     Bornamp(nth,maf,mai) = Bornamp(nth,maf,mai) +
     >                  srylm * cgcoeff * ff(nth,lam) * (0.0,1.0) ** lam
                  enddo
               endif 
            enddo               ! end maf loop
         enddo                  ! end mai loop
      enddo                     ! end lam loop
      return
      end
      
      
      subroutine formfac(hlike,nf,ni,rlf,rli,temp,maxtemp,lam,nstep,
     >   q,vc,ff,Nmax,namax,C,ze)
      include 'par.f'
      logical hlike             !,zeroc
      real temp(maxr), chi(maxr), q(0:200), vc(0:200), ff(0:200),
     >   ucentr(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /double/njdouble,jdouble(22)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common/smallr/ formcut,regcut,expcut,fast,match
      logical fast,match
      complex phase, sigc
      common /CcoefsE/ E(KNM)
      double precision  C(Nmax,namax,namax), E, ortint
      common /helium/ la(KNM), sa(KNM), lpar(KNM), np(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      data sq2pi,ucentr/39.4784176,maxr*0.0/
      integer na, nam, sa
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      common /atompol/ atompol(maxr), ipolmin, ipolmax, ipol
      common /di_el_core_polarization/ gamma, r0, pol(nmaxr) 
      real chin(maxtemp,0:200)
      integer istartchi(0:200), istopchi(0:200)
      real tmp_hat, tmp_hat1, tmp_atom_pol(0:200)
      common /noblgegas_switch/ i_ng    ! set in  cc/main.f
      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      common /debye/ dbyexists, dblp, rmudeb, zdby
      logical  dbyexists
  
      hat(rl) = sqrt(2.0 * rl + 1.0)

      pi = acos(-1.0)
      eta = 0.0
      rlam = lam
      rm = 0.0
      rr = 0.0
      do nth = 0, 200
         ff(nth) = 0.0
      enddo 
      if (hlike) then
         const = cgc(rli,rm,rlam,rm,rlf,rm)
         if (const.ne.0.0) then
            do i = 1, maxtemp
               rr = rr + temp(i) * rmesh(i,1) * rmesh(i,1) * rmesh(i,3)
            enddo
            rr = rr * hat(rli) * hat(rlam) * const
            do nth = 0, 200, nstep
               if (q(nth).gt.0.0) then  !1e-5
                  qsq = q(nth) * q(nth)
                  call regular(lam,qsq,eta,ucentr,cntfug(1,lam),-1,
     >               rmesh,maxtemp,jdouble,njdouble,regcut,expcut,
     >               chi,istart,istop,phase,sigc)
                  do i = 1, maxtemp
                     chi(i) = chi(i) * rmesh(i,3) / rmesh(i,1)
                  enddo 
c     9.05.2000
c     modify array  chin(i,nth)  due to two-electron pol.potential
                  if(lam.eq.1.and.gamma.ne.0.0) then
                     sump = 1.0
                     do i=istart,istop
                        tmp = rmesh(i,1)*rmesh(i,1)*pol(i) - 1.0
                        if(abs(tmp).lt.expcut) go to 43
                        sump = sump + tmp*chi(i)
                     end do
c$$$  43               const2 = - gamma*q(nth)*sump/(4.0*pi)            
 43                  const2 = - gamma*qsq*sump/3.0
                     do i=istart,istop
                        chi(i) = chi(i) +
     >                     pol(i)*rmesh(i,3)*const2
                     end do
                  end if
c     ------------
                  sump = 0.0
                  do i = 1, maxtemp
                     sump = sump + temp(i) * chi(i)
                  enddo
                  ff(nth) = sump * hat(rli) * hat(rlam) * const / q(nth)
               endif 
            enddo 
         endif 
      else
         if (sa(ni).ne.sa(nf)) return
c     Make array of chin(i,nth) - it is the same inside s.p. loops. and 
c     can be made outside of them.
         do nth = 0, 200, nstep
            if (q(nth).gt.0.0) then  !1e-5
               qsq = q(nth) * q(nth)
               call regular(lam,qsq,eta,ucentr,cntfug(1,lam),-1,
     >              rmesh,maxtemp,jdouble,njdouble,regcut,expcut,
     >              chin(1,nth),istartchi(nth),istopchi(nth),
     >              phase,sigc)
               do i = 1, maxtemp
c$$$                  if (i.lt.10) print*,qsq,i,sin(q(nth)*rmesh(i,1)),
c$$$     >               chin(i,nth)
                  chin(i,nth) = chin(i,nth) * rmesh(i,3) / rmesh(i,1)
               enddo
c     modify array  chin(i,nth)  due to two-electron pol.potential
               if(lam.eq.1.and.gamma.ne.0.0) then
                  sump = 1.0
                  do i=istartchi(nth),istopchi(nth)
                     tmp = rmesh(i,1)*rmesh(i,1)*pol(i) - 1.0
                     if(abs(tmp).lt.expcut) go to 45
                     sump = sump + tmp*chin(i,nth)
                  end do
c$$$ 45               const2 = - gamma*q(nth)*sump/(4.0*pi)            
 45               const2 = - gamma*qsq*sump/3.0
                  do i=istartchi(nth),istopchi(nth)
                     chin(i,nth) = chin(i,nth) +
     >                    pol(i)*rmesh(i,3)*const2
                  end do
               end if
               if(ipol .eq.1 .and. ni .eq. nf .and. lam .eq. 0) then
                 minfun=max(istartchi(nth),ipolmin)
                 maxfun=min(istopchi(nth),ipolmax)
                 tmp_atom_pol(nth) = SUM(chin(minfun:maxfun,nth)*
     >              atompol(minfun:maxfun)*rmesh(i,1)*rmesh(i,1))
                 tmp_atom_pol(nth) = tmp_atom_pol(nth) * pi*q(nth)
              end if
            end if
         end do
c
         tmp_hat1 = hat(rli) * hat(rlf) * hat(rlam)

         if(i_ng .eq. 0) then
         li = nint(rli)
         do 67 jni1=1,nam(Ni)
            ni1 = na(jni1,Ni)
            li1=lo(ni1)
            rli1=li1
            tmp_hat = tmp_hat1 * hat(rli1) 
            do 68 jnf1=1,nam(Nf)
               nf1 = na(jnf1,Nf)
               lf1=lo(nf1)
               rlf1=lf1
               const = 0.0
               do 20 jni2=1,nam(Ni)
                  ni2 = na(jni2,Ni)
                  if(C(Ni,jni1,jni2).eq.0.0d0) goto 20
                  li2=lo(ni2)
                  rli2=li2
                  do 10 jnf2=1,nam(Nf)
                     nf2 = na(jnf2,Nf)
                     lf2=lo(nf2)
                     rlf2=lf2
                     if(C(Nf,jnf1,jnf2).eq.0.0d0.or.
     >                    ortint(ni2,nf2).eq.0d0) goto 10
                     const = const + C(Nf,jnf1,jnf2) * C(Ni,jni1,jni2)*
     >                  (-1)**(lf1+lf2+li+lam) * ortint(ni2,nf2) * 
     >                  COF6J(rlf2,rlf1,rlf,rlam,rli,rli1)
 10               continue 
 20            continue 
               const = const * cgc(rli1,rm,rlam,rm,rlf1,rm) * tmp_hat
               if (const.ne.0.0) then
                  istart = max(minf(ni1),minf(nf1))
                  istop = min(maxf(ni1),maxf(nf1))
                  sump = 0.0
                  do i = istart, istop
                     sump = sump + fl(i,ni1) * fl(i,nf1) *
     >                  rmesh(i,1) * rmesh(i,1) * rmesh(i,3)
                  enddo
                  rr = rr + sump * const
                  do nth = 0, 200, nstep
                     if (q(nth).gt.0.0) then  !1e-5
                        imax = min(istopchi(nth),maxf(ni1),maxf(nf1))
                        imin = max(istartchi(nth),minf(ni1),minf(nf1))
                                                        
                        sump = 0.0
                        do i = imin, imax
                           sump = sump + fl(i,ni1)*fl(i,nf1)*chin(i,nth) 
                        enddo
                        ff(nth) = ff(nth) + const * sump / q(nth)
                     endif
                  enddo    ! end nth loop
               end if
 68         continue 
 67      continue 
         elseif(i_ng .eq. 1) then
            nfcel = 4*l_ng + 1  ! number of 'frozen'-core electrons (5 : for Ne 2p^5)

         li = nint(rli)
         lf = nint(rlf)
         do  jni1=1,nam(Ni)
            ni1 = na(jni1,Ni)
            li1=lo(ni1)
            rli1=li1
            tmp_hat = tmp_hat1 * hat(rli1) 
            do  jnf1=1,nam(Nf)
               nf1 = na(jnf1,Nf)
               lf1=lo(nf1)
               rlf1=lf1
               const = 0.0
               do jni2=1,nam(Ni)
                  ni2 = na(jni2,Ni)
                  if(C(Ni,jni1,jni2).eq.0.0d0) cycle
                  li2=lo(ni2)
                  rli2=li2
                  do jnf2=1,nam(Nf)
                     nf2 = na(jnf2,Nf)
                     lf2=lo(nf2)
                     rlf2=lf2
                     if(C(Nf,jnf1,jnf2).eq.0.0d0) cycle
                     
                     const_ci = C(Nf,jnf1,jnf2) * C(Ni,jni1,jni2)

                     coefsym = 1.0
                     if(ni2 .eq. ni1 .and. nf2 .eq. nf1) then
                        i_12 = 1
                     elseif(ni2 .ne. ni1 .and. nf2 .ne. nf1) then
                        i_12 = 1
                     else
                        i_12 = 2                        
                        coefsym = sqrt(nfcel+1d0)   ! = sqrt(2) for lng=0, =sqrt(6) for lng=1, == sqrt(2*(2*lng+1))
                     endif

c     start <4l+1|f_k|4l+1> ME
                     if(ortint(ni2,nf2).ne.0d0 
     >                    .and. i_12.eq.1) then
                     
                     const1 = (-1)**(lf1+lf2+li+lam) * ortint(ni2,nf2) * 
     >                    COF6J(rlf2,rlf1,rlf,rlam,rli,rli1) *
     >                    cgc(rli1,rm,rlam,rm,rlf1,rm) * 
     >                    tmp_hat1 * hat(rli1) 
                     if(lam .eq. 0) then ! this is noble gas modification for ME 
                        const1 = const1 * nfcel
                     else
                        const1 = const1 * (-1)**(lam+1) 
                     endif

                     const = const1 * const_ci / ze
                     if (const.ne.0.0) then
                        istart = max(minf(ni1),minf(nf1))
                        istop = min(maxf(ni1),maxf(nf1))
                        sump = 0.0
                        do i = istart, istop
                           sump = sump + fl(i,ni1) * fl(i,nf1) *
     >                          rmesh(i,1) * rmesh(i,1) * rmesh(i,3)
                        enddo
                        rr = rr + sump * const
                        do nth = 0, 200, nstep
                           if (q(nth).gt.0.0) then !1e-5
                       imax = min(istopchi(nth),maxf(ni1),maxf(nf1))
                       imin = max(istartchi(nth),minf(ni1),minf(nf1))
                              
                              sump = 0.0
                              do i = imin, imax
                                 sump = sump + 
     >                           fl(i,ni1)*fl(i,nf1)*chin(i,nth) 
                              enddo
                              ff(nth) = ff(nth) + const*sump/q(nth)
                           endif
                        enddo   ! end nth loop
                     end if

                     endif
c     start <l2|f_k|l2'> ME
                     if(ortint(ni1,nf1).ne.0d0) then
                     const1 = (-1)**(li1+li2+lf+lam) * ortint(ni1,nf1) * 
     >                    COF6J(rlf1,rlf,rlf2,rlam,rli2,rli) *
     >                    cgc(rli2,rm,rlam,rm,rlf2,rm) 
     >                    * tmp_hat1 * hat(rli2) 
     >                    * coefsym 

                     const = const1 * const_ci / ze
                     if (const.ne.0.0) then
                        istart = max(minf(ni2),minf(nf2))
                        istop = min(maxf(ni2),maxf(nf2))
                        sump = 0.0
                        do i = istart, istop
                           sump = sump + fl(i,ni2) * fl(i,nf2) *
     >                          rmesh(i,1) * rmesh(i,1) * rmesh(i,3)
                        enddo
                        rr = rr + sump * const
                        do nth = 0, 200, nstep
                           if (q(nth).gt.0.0) then !1e-5
                       imax = min(istopchi(nth),maxf(ni1),maxf(nf1))
                       imin = max(istartchi(nth),minf(ni1),minf(nf1))
                              
                              sump = 0.0
                              do i = imin, imax
                                 sump = sump + 
     >                                fl(i,ni2)*fl(i,nf2)*chin(i,nth) 
                              enddo
                              ff(nth) = ff(nth) + const*sump/q(nth)
                           endif
                        enddo   ! end nth loop
                     end if

                     endif

                  enddo
               enddo

            enddo
         enddo

         else
            print*,'processt.f: wrong value for i_ng=',i_ng
            stop
         endif







      endif

      if (nf.eq.ni.and.lam.eq.0) then
         do nth = 0, 200, nstep
            ff(nth) = ff(nth) - hat(rli) * (1.0 - vc(nth))
c     account of additional atom polarazability:
            if(ipol .eq.1 ) then
               ff(nth) = ff(nth) + hat(rli) *  tmp_atom_pol(nth)
            endif
         enddo
      endif 
      rmusq = rmudeb * rmudeb
      do nth = 0, 200, nstep
         ff(nth) = ff(nth) * 2.0 / (q(nth) * q(nth) + rmusq + 1e-20)
      enddo 
      if (q(0).le.0.0) then  !1e-5
c$$$         if (ni.ne.nf) then
c$$$C  Here we have degeneracy as in the case of e-H 2s -> 2p. Set both
c$$$C  the zero and the next formfactor to large numbers. The latter is so
c$$$C  that the integrated CS came out to be large.
c$$$            ff(0) = 1e10
c$$$            ff(nstep) = 1e10
c$$$         else
         if (lam.eq.0) then
            ff(0) = - rr/3.0 + 2.0 * vc(0) * hat(rli)
         else if (lam.eq.2) then
            ff(0) = 2.0 / 15.0 * rr
         else
            ff(0) = 0.0
         endif
         ! Mark: Trying with 0 for elastic scattering.
         if (dbyexists) then
            ff(0) =0.0
         end if
      endif 
      do nth = 0, 200, nstep
         ff(nth) = ze * ff(nth) / sq2pi
      enddo 
      return
      end
      
      subroutine makee2eamp(te2es,te2ef,nsmax,jstop,latop,lfa,slowery,
     >   rkf,rki,nunit,nt,ne,sweight,nstep,projectile,target,
     >   ninttype,rmpow,rnpow,chani,chanf,sums,sumf,chun,nze)
      include 'par.f'
      complex te2es(0:nsmax,0:jstop,-latop:latop),cfacs(0:1),cfacf(0:1),
     >   fs(0:200,0:1,-lamax:lamax),ff(0:200,0:1,-lamax:lamax),
     >   te2ef(0:nsmax,0:jstop,-latop:latop)
      real*8 rylm
      dimension units(3),spinw(0:1),sweight(0:200),sums(0:1),sumf(0:1)
      character ch*1, xunit(3)*8,chanf*3,chani*3,fn*2,a3*3
      character projectile*(*), target*(*),chun*(*)
      ch(i)=char(i+ichar('0'))
      DATA XUNIT/' a. u.  ','PI*A0**2','cm**2/eV'/
c$$$      data units/1.0,0.3183099,28.002e-18/
C  Divide cm^2 by eV = 27.21 a.u.
      data units/1.0,0.3183099,1.0291e-18/
      data pi/3.1415927/

      if (nunit.eq.2) print*,'Warning: have pi a0**2 units in (e,2e)'
      chun = xunit(nunit)
      do mf = - lamax, lamax
         do ns = 0, nsmax
            do nth = 0, 200
               fs(nth,ns,mf) = (0.0,0.0)
               ff(nth,ns,mf) = (0.0,0.0)
            enddo 
         enddo
      enddo

      do ns = 0, nsmax
         call getspinw(chani(1:1),ns,nze,spinw(ns),si)
      enddo 
      if (nunit.gt.0) unit = units(nunit)

c Alisher's changes               
      if (chanf(1:1).eq.'p') then
         rks = sqrt(slowery/2.0)
         fastery = rkf**2/2.0
      else 
         rks = sqrt(slowery)
         fastery = rkf**2
      endif
      
      c1s = rkf/rki * (2.0 * pi) ** 4 * unit 
      c1f = rks/rki * (2.0 * pi) ** 4 * unit 

      lia = 0
      rlia = lia
      rlfa = lfa
      do j = 0, jstop
         do ipar = 0, min(j,0)
            rj = j
            iparity = (-1)**(j + ipar)
            do li = abs(j - lia), j + lia
               if ((-1)**(li+lia).eq.iparity) then
                  rli  = li
                  do lf = abs(j - lfa), j + lfa
c$$$                     if ((-1)**(lf+lfa).eq.iparity) then
                        rlf  = lf
                        c2 = sqrt((2.0*rli+1.0)) / sqrt(pi * 4.0)
                        do mi = - lia, lia
                           rmi = mi
                           c3=cgc(rli,0.0,rlia,rmi,rj,rmi)*c2
                           do mf = - lfa, lfa
                              rmf = mf
                              c4=cgc(rlf,rmi-rmf,rlfa,rmf,rj,rmi)*c3
                              if (abs(c4).gt.1e-10) then
                                 do ns = 0, nsmax
                                    cfacs(ns) = c4*
     >                                 (0d0,1d0)**dfloat(li-lf)*
     >                                 te2es(ns,j,j-lf)
                                    cfacf(ns) = c4*
     >                                 (0d0,1d0)**dfloat(li-lf)*
     >                                 te2ef(ns,j,j-lf)
c$$$                                       print*,'J,lfa,lf,ns,T',J,lfa,lf,
c$$$     >                                    ns,te2e(ns,j,j-lf)
                                 enddo
                                 do nth = 0, 180, nstep
                                    theta = float(nth)
                                    thrad = theta * pi / 180.0
                                    sh=rylm(lf,mi-mf,dble(thrad))
                                    do ns = 0, nsmax
                                       fs(nth,ns,mf)=sh*
     >                                    cfacs(ns)+fs(nth,ns,mf) 
                                       ff(nth,ns,mf)=sh*
     >                                    cfacf(ns)+ff(nth,ns,mf) 
                                    enddo
                                 enddo
                              endif 
                           enddo ! end of mf loop
                        enddo   ! end of mi loop
c$$$                     endif      ! end of parity check for LF + LFA
                  enddo         ! end of lf loop
               endif            ! end of parity check for LI + LIA
            enddo               ! end of li loop
         end do                 ! end of parity loop
      enddo                     ! end of J loop
      
      if (chanf(1:1).eq.' ') then
         fn = 'd'//chanf(3:3)
      else 
         fn = chanf(1:1)//chanf(3:3)
      endif
      a3 =  fn//ch(ne)
      print'("writing file: amp.",a3)', a3
      open(42,file='amp.'//fn//ch(ne))

      write(42,"(a8,' - ',a6,'ionization:',f8.3,'eV on ',a3,
     >   ' -> ',a2,2(f8.3,'eV'))") projectile,target,13.6058*
     >   rki**2,chani,fn,13.6058 * slowery,
     >   13.6058 * fastery !   rkf**2*13.6058

      write(42,'(''Amplitudes generated using NINTTYPE = '',i2,
     >   '', M = '',f3.1,'' and N = '',f3.1)') ninttype,rmpow,rnpow
      
      write(42,'(''When the amplitudes are squared the result '',
     >   ''contains the factor Kf/Ki/(2Li+1) ''/
     >   ''and is in units of: '',a8)') xunit(nunit)

      write(42,'(
     >   ''The CCC amplitudes may be read in the following way:''/
     >   ''      read(n,*) lf, li, nth, rkf, rki, tsp1, nsm,'',
     >   '' (spinw(ns), ns=0, nsm)''/
     >   ''      do mi = - li, li''/
     >   ''         do mf = - lf, lf''/
     >   ''            read(n,*) mfprev, miprev''/
     >   ''            if (mfprev.ne.mf.or.miprev.ne.mi) stop'',
     >   '' "incorrect format" '')')
      write(42,'(
     >   ''            do ith = 1, nth''/
     >   ''               read(n,*) theta(ith), (freal(ns),'',
     >   '' fimag(ns), ns = 0, nsm)''/
     >   ''               f(ith,mf,mi,0) = cmplx(freal(0),'',
     >   '' fimag(0))'')')
      write(42,'(
     >   ''               f(ith,mf,mi,1) = cmplx(freal(1),'',
     >   '' fimag(1))''/
     >   ''            end do''/
     >   ''         end do''/
     >   ''      end do'')')
      call getspinw(chanf(1:1),0,nze,sw,sf)
      write(42,'(2i2,i4,2f8.5,2i2,2f8.5,
     >   " Lf,Li,Nth,Kf,Ki,2Sf+1,Nsm,SpinW(ns)")')
     >   lfa, lia, nth,  rkf, rki,
     >   nint(2.0*sf+1.0),nsmax,(spinw(ns),ns=0,1)
      xs = sqrt(c1s/(2.0 * rlia + 1.0)) 
      xf = sqrt(c1f/(2.0 * rlia + 1.0))
      sums(:) = 0.0
      sumf(:) = 0.0
      do mli = -lia, lia
         do mlf = -lfa, lfa
            write(42,*) mlf, mli, '  MF,  MI'
            do nth = 0, 180, nstep
               write(42,'(f5.1,1p,2(2(1x,2e11.3),4x))') float(nth),
     >            (xs*fs(nth,ns,mlf),ns = 0,nsmax),
     >            (xf*ff(nth,ns,mlf),ns = 0,nsmax)
            enddo
            do ns = 0, nsmax
               call integrate(fs(0,ns,mlf),sweight,nstep,res,res2)
               sums(ns) = sums(ns) + spinw(ns) * res * xs ** 2
               call integrate(ff(0,ns,mlf),sweight,nstep,res,res2)
               sumf(ns) = sumf(ns) + spinw(ns) * res * xf ** 2
            enddo 
         enddo
      enddo
c$$$      print*,'The slow and fast TCS:',sums,sumf
      close(42)
      return
      end
      
      subroutine intamp(amp,thfac,rkf,rki,sweight,res,csmt,err1,err2)
      real sweight(0:200,10),test(6),csmt(5)
      complex amp(0:200)
      double precision yin(201),xin(201),yout(180*32+1),xout(180*32+1)
      data pi/3.141592654/

      do n = 1, 6
         test(n) = 0.0
      enddo
      csmt(:) = 0.0
      qmin = abs(rkf - rki)
      qmax = rkf + rki
      if (thfac.eq.0.0.or.qmin.eq.0.0) then
C  Here the amplitude is assumed to be on a linear theta grid
         thdiv = 1.0
         do n = 1, 3
            do i = n, 180, n
               test(n) = test(n) + sweight(i,n) *
     >            abs(amp(i)) ** 2 * sin(i * pi / 180.0)
            enddo
c$$$            print*,'testing theta:',test(n)
         enddo
         if (qmin.eq.0.0) then
            do l = 1, 5
               do i = 1, 180
                  csmt(l)  = csmt(l) + sweight(i,n) *
     >               abs(amp(i)) ** 2 * sin(i * pi / 180.0) *
     >               (1.0 - (cos(i * pi / 180.0))**l)
               enddo
            enddo 
         endif 

         do nth = 0, 200
            q = qmin + nth * (qmax - qmin) / 180.0 
            theta = float(nth)
            if (nth.gt.180) theta = (nth - 180) / 21.0
            thrad = theta * pi / 180.0 
C     Q = Ki - Kf  (as vectors)               
            qsq = rkf * rki * (dble(rkf/rki) + dble(rki/rkf) -
     >         dble(2.0 * cos(dble(thrad))))
            q = sqrt(qsq)
            if (nth.eq.0) then
               i = 1
            elseif (nth.le.180) then
               i = 21 + nth
            else
               i = nth - 179
            endif
            xin(i) = q
            yin(i) = abs(amp(nth))**2
         enddo
      else 
C  Here the amplitude is assumed to be on a linear q grid
         thdiv = thfac
         do nth = 0, 180
            yin(nth+1) = abs(amp(nth))**2
c$$$            q = qmin + nth * (qmax - qmin) / 180.0 / thfac
            q = exp(
     >         log(qmin) + nth * (log(qmax) - log(qmin))/ 180.0)
            xin(nth+1) = q
c$$$            xin(nth+1) = nth * pi / 180.0
C  Note that sweight = 0.0 for some NTH, NN pairs
            do n = 1, 3
               test(n) = test(n) +
     >            sweight(nth,n)*yin(nth+1)*xin(nth+1)**2
c$$$               test(n) = test(n) + sweight(nth,n)*yin(nth+1)*xin(nth+1)
            enddo 
         enddo
         do n = 1, 3
c$$$            test(n) = test(n) / rki / rkf * (qmax-qmin) / pi / thfac
            test(n) = test(n) / rki / rkf *
     >         (log(qmax)-log(qmin)) / pi 
c$$$            print*,'testing q:',test(n)
         enddo
      endif 
      err1 = abs((test(1)-test(2))/(test(1)+1e-30))*100.0
      err2 = abs((test(1)-test(3))/(test(1)+1e-30))*100.0
      res = test(1)

      if ((err1.gt.0.1.or.err2.gt.0.1).and.thfac.eq.0.0) then
C  Here the amplitude will be extrapolated on to a linear q grid
C  The value of NP is used to determine the interpolation. At high energies it
C  has some effect so some care needs to be taken!
         if (yin(2).gt.yin(1)) then
            np = 3
         else
            np = 4
         endif
         np = 2
         print*,'Interpolation with NP =',np
         if (rki.eq.rkf) then
C  This logic avoids the problem of 0**0 on the Cray.
            np = 0
         else 
            do nth = 1, 200
               yin(nth) = yin(nth) * xin(nth) ** np
c$$$               (1d0 - exp( - xin(nth) ** np ))
            enddo
         endif 
         ntwo = 5
         nthi = 180 * 2 ** ntwo + 1
         do nth = 1, nthi
            q = qmin + (nth-1) * (qmax - qmin) / (nthi-1) / thdiv
            xout(nth) = q
         enddo
c$$$         if (qmin.gt.0.0) then
c$$$            do nth = 1, 181
c$$$               write(67,*) xin(nth),yin(nth)
c$$$            enddo
c$$$            write(67,*)
c$$$         endif 
c$$$      
         call intrpl(200,xin,yin,nthi,xout,yout)
c$$$         if (qmin.gt.0.0) then
c$$$            do nth = 1, nthi
c$$$               write(68,*) xout(nth),yout(nth)
c$$$            enddo
c$$$            write(68,*)
c$$$         endif 
         do n = 0, 5
            tmp = 0.0
            test1 = 0.0
            h = xout(2**n+1)
            h = (xout(2**n+1) - xout(1)) 
            nprev = 1
            do nth = 1, nthi, 2**n
               weight = float(nprev) / 3.0 * h
               contrib = weight * yout(nth) * xout(nth) ** (1-np)
c$$$     >           / (1d0 - exp( - xout(nth) ** np ))
               tmp = tmp + contrib
               contribt = xout(nth)**2 * weight
               test1 = test1 + contribt
               if (nprev.eq.4) then
                  nprev = 2
               else
                  nprev = 4
               endif
            enddo
            test(n+1) = (tmp - contrib / 2.0) / rki / rkf
            test1 = test1 - contribt / 2.0
c$$$            print*,'n,ICS, test:',n, test(n+1),
c$$$     >         test1 * 3.0 / (qmax**3-qmin**3)
         enddo 
         err1 = abs((test(1)-test(2))/(test(1)+1e-30))*100.0
         err2 = abs((test(1)-test(3))/(test(1)+1e-30))*100.0
         if (abs((test(1)-res)/res).gt.0.2) print*,
     >      'CAUTION: extrapolation using NP had an effect > 20%',
     >      res, test(1), np
      endif 
      res = test(1) * 2.0 * pi
      csmt(:) = csmt(:) * 2.0 * pi
      return
      end
      
      complex function getTextr(j,js,t1,t2)
      integer j, js
      complex t1, t2
      complex res
      real*8  t1_r, t2_r, t1_i, t2_i, f1_r, f1_i, f2_r, f2_i
      real*8  rn, tmp
      

      t1_r = real(t1)
      t2_r = real(t2)
      t1_i = aimag(t1)
      t2_i = aimag(t2)

      
      f1_r = 0d0
      f2_r = 0d0
      f1_i = 0d0
      f2_i = 0d0


      f1_r = real(t1) * (2*js-1) * (2*js+1) * (2*js+3) /
     >     (2*j-1) / (2*j+1) / (2*j+3)
      f1_i = aimag(t1) * (2*js-1)*(2*js+1) *
     >     (2*js+3) * (2*js+5) * (2*js+7) /
     >     (2*j-1) / (2*j+1) / (2*j+3) / (2*j+5) / (2*j+7)
      

      if(t2_r .ne. 0.0)  then
         if(abs(t1_r) .lt. abs(t2_r)) then 
            f1_r = t1_r * (t1_r/t2_r)**(j-js)   
!            tmp = log(dble(js-1)/dble(js))
!            rn = log(t1_r/t2_r) / tmp
!            f2_r = 0.0 !t1_r * (dble(js)/dble(j))**(rn)            
!            f1_r = (f1_r + f2_r) /2.0
         endif

      endif



      if(t2_i .ne. 0.0)  then
         if(abs(t1_i) .lt. abs(t2_i)) then 
            f1_i = t1_i * (t1_i/t2_i)**(j-js)   
!            tmp = log(dble(js-1)/dble(js))
!            rn = log(t1_i/t2_i) / tmp
!           f2_i = 0.0 ! t1_i * (dble(js)/dble(j))**(rn)
!            f1_i = (f1_i + f2_i) /2.0
         endif

      endif
      
      res = cmplx( f1_r, f1_i )

      getTextr = res

      return
      end


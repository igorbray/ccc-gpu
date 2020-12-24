C  File rest.f
      subroutine solvet(iex,bb,vmat,gk,wk,weightk,nchtop,nqm,noprint,
     >   nopen,etot,lg,vdon,phasel,nent,instate,isecond,nunit,sigma,
     >   vmatop,nchopt,ovlpn,phaseq,ve2ed,ve2ee,dphasee2e,ephasee2e,
     >   nchtope2e,lfast,lslow,slowery,nznuci,zasymi,npk,ne2e,nchanmax,
     >   projectile,target,uba,nsmax,canstop,second,BornICS,tfile,
     >   csfile,e2efile,ichi,theta,nnbtop,ovlpnn,vmatp,scalapack,chil)
      use date_time_module
      include 'par.f'
      parameter (nchmaxe2e=1)
      integer npk(nchtop+1), nopen(0:1), instate(100)
      real vmat(ichi,ichi+1),vmatp(1,0:nsmax),
     >   bb(ichi,nchtop,2,0:nsmax), !vmatp(ichi*(ichi+1)/2,0:nsmax),
     >   ve2ed(nchmaxe2e,npk(nchtop+1)-1), ovlpn(knm), 
     >   vdon(nchan,nchan,0:1), ve2ee(nchmaxe2e,npk(nchtop+1)-1)
      complex ovlpnn
      complex vmatop(kmax,kmax,0:nchanop,nchanop),wk(kmax*nchan),
     >   cwork(66 * nchan),eigv(nchan),phaseq(knm),tcmplx,vs(1)
      complex*16 coulphase
      complex det(2),sigma(nchan),te2eH(-lamax:lamax,0:lamax),
     >   dphasee2e(nchane2e),ephasee2e(nchane2e),te2e(nchane2e),
     >   ve2e(nchane2e),phaseslow,sigc,ctmp,smat(nchan,nchan)
      complex ton(nchan,nchan,0:1),von(nchan,nchan),t2nd(nchan,nchan)
      complex wton(nchan,nchan), wvon(nchan,nchan), wt2nd(nchan,nchan),
     >   wvdon(nchan,nchan), phasel(kmax,nchan), tdist(nchan)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension jlm(nchan),non(nchan,2*lamax+1),partcs(nchan,nchan,0:1),
     >   weightk(kmax),gk(kmax,nchan), nychan(nchan), rwork(3*nchan),
     >   temp(maxr), err(nchan,nchan), psislow(maxr), slowery(ncmax),
     >   esum(0:1),rcond(0:1),ovlpnn(knm,nnmax)
      COMMON /gstspin/   iSpin, meta
      common/smallr/ formcut,regcut,expcut,fast,match,analyticd,packed
      common /double/id,jdouble(22)
      common /radpot/ ucentr(maxr)
      logical fast,match,analyticd,packed,scalapack
      real sigtop(nchan,0:1),imag
      real BornICS(knm,knm), BornPCS(nchan,nchan)
      logical same, sames, noprint(nchan), uba(nchan), canstop, second,
     >   instates,bwork(1)
      data pi,method,sames/3.1415927,1,.false./
      character ch,spin(0:1),chnopen*20,e2efile(10)*60,fmt*80,
     >   chan(knm)*3, target*6, tfile*80, csfile*80, projectile*8
      common /charchan/ chan
      common /chanen/ enchan(knm)
      save sames
      data spin/'+','-'/
      logical exists, select

      nch = 1
      ldvs=1
      n = 0
      do while (nch.le.nchtop)
         if (uba(nch)) then
            if (n.eq.0) then
               print '("NCHTOP, UBA channels:",i3,",",$)', nchtop
               n = 1
            endif
            nchs = nch
            do while (uba(nch).and.nch.le.nchtop)
               nch = nch + 1
            enddo
            print '(1x,i3," to ",i3,$)', nchs, nch - 1
         endif 
         nch = nch + 1
      enddo
      if (n.ne.0) print*
      
      do nch = 1, nchtop
         ai = aimag(phasel(1,nch))
c$$$         si = aimag(sigma(nch))
         if (lg.gt.20.and.abs(ai).gt.1e-8) then
            if (nch.eq.1)
     >         print*,'As J > 20 reset IM(phase) to 0 from:',ai
            ai = 0.0
         endif 
         tdist(nch) = - phasel(1,nch) * ai * sigma(nch) ** 2 /
     >      (pi * gk(1,nch)+1e-10)
         phasel(1,nch) = phasel(1,nch) * sigma(nch)
c$$$         tdist(nch) = - (phasel(1,nch) * ai - sigma(nch) * si) /
c$$$     >      (pi * gk(1,nch))
         jlm(nch) = 0
         do la = 1, 2*lamax+1
            non(nch,la) = 0
         enddo 
      enddo

      
      if (abs(tdist(1)).ne.0.0) print*,'TDIST:',tdist(1)
      
      same = .true.
      sames = .false.
      if (iex.eq.0.and.mod(nznuc - nint(zasym),2) .eq. 1) sames = .true.

      do nchi = 1, nchan
         sigtop(nchi,0) = 0.0
         sigtop(nchi,1) = 0.0
         do nchf = 1, nchan
            BornPCS(nchf,nchi) = 0.0
            partcs(nchf,nchi,0) = 0.0
            partcs(nchf,nchi,1) = 0.0
         enddo
      enddo
      rcond(:) = 0.0
      nsmin = 0
      if (iSpin.eq.1) nsmin = 1
      do ns = nsmin, nsmax
         call date_and_time(date,time,zone,valuesin)
         print '("GETTMAT entered at:",a10)',time
C  The following caused a bug for helium when two spins are used, but hopefully fixed with above
         if (sames.and.ns.eq.1) then
            do nchi = 1, nchtop
               do nchf = 1, nchtop
                  ton(nchf,nchi,1) = ton(nchf,nchi,0)
               enddo
            enddo
         else 
            call gettmat(lg,ns,nchtop,iex,npk,vmat,wk,gk,phasel,sigma,
     >      bb(1,1,1,ns),tdist,von,etot,t2nd,ton(1,1,ns),err,vmatop,
     >      nchopt,nsmax,ve2ed,ve2ee,dphasee2e,ephasee2e,nchtope2e,
     >      nchmaxe2e,ve2e,te2e,nds,second,uba,ovlpn,theta,ichi,csfile,
     >      projectile,slowery(1),det,rcond(ns),vmatp(1,ns),packed,
     >         phaseq,scalapack,chil)
         endif 
         call date_and_time(date,time,zone,valuesout)
         print '("GETTMAT exited at:",a10," diff (secs):",i5)',
     >      time,idiff(valuesin,valuesout)
         valuesin=valuesout
         call update(6)
         r = float((-1)**ns * min(1,iex))
         do ne = 1, abs(ne2e)
            do nchi = 1, 1
               call getchinfo (nchi,nt,lg,temp,maxpsii,ei,lia,nia,li)
               fastery = etot - slowery(ne)
               spinw = (2.0 * ns + 1.0)/4.0
               do nchf = 1, nchtope2e
                  call getche2e(nchf,lg,lfast,lslow,lfa,lf)
                  nchft = 2
                  nchftreq = -1
                  do while (nchft.gt.1)
                     call getchinfo (nchft,nft,lg,temp,maxpsif,ef,
     >                  lfat,nfa,lfp)
                     if (abs(ef-slowery(ne))/slowery(ne).lt.1e-2
     >                  .and.lfat.eq.lfa.and.lfp.eq.lf) then
c$$$                        print*,'nchft:',nchft
                        nchftreq = nchft
                        nftreq = nft
                        nchft = -10
                     endif
                     nchft = nchft + 1
                  enddo 

                  nche = (ne - 1) * nchtope2e + nchf
c$$$                  call getchinfo(nchf,ntp,lg,temp,maxpsif,ef,lfa,nfa,lf)
c$$$                  fastery = gk(1,nchf) ** 2
                  const = spinw * (2.0 * pi) ** 4 / (4.0 * pi) *
     >               (2.0 * lg + 1.0) * sqrt(fastery) / gk(1,nchi) /
     >               (2.0 * lia + 1.0)
                  if (nchftreq.gt.0) then
                     write(88+ne,'(4i3,1p,3(2e12.3,2x))') lg,lfa,lf,ns,
     >                  ve2e(nche), ton(nchftreq,nchi,ns) *
     >                  ovlpn(nftreq) * phaseq(nftreq) / sqrt(ef),
     >                  te2e(nche)
                  else
                     write(88+ne,'(4i3,1p,3(2e12.3,2x))') lg,lfa,lf,ns,
     >                  ve2e(nche), te2e(nche)
                  endif 
               enddo
            enddo
            close(88+ne)
            call openatend(88+ne,e2efile(ne))
c$$$            open(88+ne,file=e2efile(ne),status='old',fileopt='eof')
            call update(88+ne)
         enddo 

C  The following initialization uses NCHAN instead of NCHTOP because
C  The partial cross sections are labeled by the state number which
C  may be greater than than NCHTOP when IPAR = 1.

         do nchi = 1, nchan
            sigtop(nchi,ns) = 0.0
            do nchf = 1, nchan
               partcs(nchf,nchi,ns) = 0.0
               smat(nchf,nchi) = (0.0,0.0)
c$$$            print'(1p,10(2e10.3))',(von(nchf,nch),nchf=1,nchtop)
c$$$            print'(1p,10(2e10.3))',(ton(nchf,nch,ns) * phasel(1,nchf)
c$$$     >          * phasel(1,nch), nchf=1,nchtop)
c$$$            print*
            enddo
         enddo 

         do nchi = 1, nchtop
            smat(nchi,nchi) = (1.0,0.0)
            do nchf = 1, nchtop
               ctmp = conjg(sigma(nchf) * sigma(nchi))*(0.0,1.0)
               gkfac = sqrt(abs(gk(1,nchi) * gk(1,nchf)))*2.0*pi
               smat(nchf,nchi) = smat(nchf,nchi) - gkfac *
     >            (ton(nchf,nchi,ns) * ctmp)
            enddo
         enddo

         inquire(file='smatwrite',exist=exists)
         if(exists) then
            energy = (etot - enpsinb(ninc,linc)) * 13.6058
            call writesmat(lg,ns,ipar,nchtop,nchan,smat,maxr,energy)
         endif

c$$$         notun = 0
c$$$         do nchi = 1, nchtop
c$$$            do nchf = 1, nchtop
c$$$               if (nchi.eq.nchf) then
c$$$                  ctmp = (-1.0,0.0)
c$$$               else
c$$$                  ctmp = (0.0,0.0)
c$$$               endif
c$$$               do nchn = 1, nchtop
c$$$                  ctmp = ctmp + conjg(smat(nchf,nchn)) * smat(nchn,nchi)
c$$$               enddo
c$$$               err(nchf,nchi) = abs(ctmp)
c$$$               if (abs(ctmp).gt.1e-4.and.nchopt.eq.0) then
c$$$                  notun = notun + 1
c$$$c$$$                  print*,'nchf,nchi,err(nchf,nchi)',
c$$$c$$$     >               nchf,nchi,err(nchf,nchi)
c$$$               endif 
c$$$            enddo
c$$$         enddo
c$$$         if (notun.gt.0) print*,'S matrix not unitary', notun
c$$$         call date_and_time(date,time,zone,valuesout)
c$$$         print '("Unitarity checked at:  ",a10,", diff (secs):",i5)',
c$$$     >      time, idiff(valuesin,valuesout)
c$$$         call update(6)
         lwork = 66 * nchan
         ldvs = 1
         esum1 = 0.0
         esum2 = 0.0
         energyin = (etot - enpsinb(ninc,linc)) * 13.6058
         sum = 0d0
         n = 1
         do nch = 1, nchtop
            sum = sum + abs(smat(n,nch))**2*
     >         atan2(imag(smat(n,nch)),real(smat(n,nch)))
         enddo 
c$$$         print'("J,S,Ein, smat(1,1), abs and arg/2:",2i2,1p,4e12.4)',
c$$$     >      lg,ns,energyin,abs(smat(1,1)),
c$$$     >      atan2(imag(smat(1,1)),real(smat(1,1)))/2.0,sum
c$$$         print'("J,S,Ein,smat:",2i2,1p,32e12.4)',
c$$$     >      lg,ns,energyin,((abs(smat(n,np)),
c$$$     >      atan2(imag(smat(n,np)),real(smat(n,np))),np=n,4),n=1,4)

         eigv(:)=(1.0,0.0)

! On the SGI CGEES crashes due to a bug in the mkl libraries.
C     Diagonalize the S matrix to obtain the eigenphases. See p.153 of Bransden
c$$$!         include 'gees.f'
#if defined _single
      call cgees('n','n',select,nchtop,smat,nchan,ndim,eigv,vs,ldvs,
     >   cwork,lwork,rwork,bwork,info)
#elif defined _double
      call zgees('n','n',select,nchtop,smat,nchan,ndim,eigv,vs,ldvs,
     >   cwork,lwork,rwork,bwork,info)
#endif
         call date_and_time(date,time,zone,valuesout)
         print '("?GEES exited at:  ",a10,", diff (secs):",i5)',
     >      time, idiff(valuesin,valuesout)
         valuesin = valuesout
         call update(6)

c$$$  call zgees('n','n',select,nchtop,smat,nchan,ndim,eigv,vs,ldvs,
c$$$     >      cwork,lwork,rwork,bwork,info)
c$$$         if (info.ne.0) print*,'INFO from CGEES:',info
C  Note that eigenphases ere defined only modulo pi due to the 2i factor
         do nch = 1, nchtop
            esum1 = esum1 + real(log(eigv(nch))/2.0/(0.0,1.0))
            esum2 = esum2 +
     >         atan2(aimag(eigv(nch)),real(eigv(nch))) / 2.0
c$$$            esum1 = esum1 + atan(aimag(eigv(nch))/real(eigv(nch))) / 2.0
         enddo
         esum(ns) = esum2
         nchi = 0
         nchip = 0
         incount = 1
         do while (nchip.le.instate(nent).and.nchi.lt.nchtop)
            nchi = nchi + 1
c$$$         do nchi = 1, lent
            call getchinfo (nchi,nchip,lg,temp,maxpsi,ei,lia,nia,li)
            instates = .false.
            do n = 1, nent
               instates = instates.or.(instate(n).eq.nchip)
            enddo 
            
            if (instates) then
               call getspinw(chan(nchip)(1:1),ns,nze,spinw,si)               
               lent = nchi
               nymax = 0
               nfp = 0
               lfp = 0
               no = 0
               nchpp = 0
               do nchf = 1, nchtop
                  call getchinfo (nchf,nchp,lg,temp,maxpsi,ef,lfa,nfa,
     >               lf)
                  if (gk(1,nchi).eq.0.0) then
                     ratio = 1.0
                     if (li.gt.0) ratio = 0.0 ! plane waves have rho**(l+1)
                  else
                     ratio = gk(1,nchf) / gk(1,nchi)
                     if (gk(1,nchf).eq.0.0.and.zasym.ne.0.0)
     >                  ratio = 1.0/gk(1,nchi) ! needs checking at thresholds
                  endif 
                  const = spinw * (2.0 * pi) ** 4 / (4.0 * pi) *
     >               (2.0 * lg + 1.0) * ratio /
     >               (2.0 * lia + 1.0)
                  lt = lg + lfa
                  lb = abs(lg - lfa)
c$$$            print*,'ns,nchf,nchi,ton,vdon:',ns,nchf,nchi,
c$$$     >               ton(nchf,nchi,ns),vdon(nchf,nchi,ns)
                  if (noprint(nchp)) then
                     nymax = nymax + 1
                     nychan(nymax) = nchf
c$$$                     if (nfa.ne.nfp.or.lfa.ne.lfp) no = no + 1
c$$$                     nfp = nfa
c$$$                     lfp = lfa
                     if (nchp.ne.nchpp) no = no + 1
                     nchpp = nchp
                     jlm(no) = lt - lb + 1
                     non(no,lf-lb+1) = nymax
c$$$               non(nchp,lf-lb+1) = nymax
                     gkfac = gk(1,nchi) * gk(1,nchf)
                     if (gkfac.eq.0.0) gkfac = 1.0
                     if (isecond.eq.-3) ton(nchf,nchi,ns) =
     >                  cmplx(real(t2nd(nchf,nchi)),0.0)
c$$$                     if (nqm.eq.1) then
                     if (ntype.eq.-1) then
C  Born approximation
                        von(nchf,nchi) = vdon(nchf,nchi,ns)/gkfac
                        t2nd(nchf,nchi) = von(nchf,nchi)
                        ton(nchf,nchi,ns) =
     >                     vdon(nchf,nchi,ns)/gkfac
                     elseif (ntype.eq.-2) then
C  Born approximation with exchange
                        ton(nchf,nchi,ns) =  von(nchf,nchi)
                     endif
c$$$                     endif 
                     wvdon(nymax,nchi)= cmplx(vdon(nchf,nchi,ns)/gkfac)
                     wton(nymax,nchi) = ton(nchf,nchi,ns)
c$$$                     write(20+ns,*) ef,sqrt(const)*
c$$$     >                  real(ton(nchf,nchi,ns)*ovlpn(nchp)),
c$$$     >                  sqrt(const)*aimag(ton(nchf,nchi,ns)*ovlpn(nchp))
c$$$                     if (ef.gt.0.0) then
c$$$                        wton(nymax,nchi) = wton(nymax,nchi) *
c$$$     >                     ovlpn(nchp) / sqrt(ef)
c$$$                     endif 
                     wt2nd(nymax,nchi) = t2nd(nchf,nchi)
                     wvon(nymax,nchi) = von(nchf,nchi)
                  end if
                  if (gk(1,nchf).ge.0.0) then
                     at2 = abs(ton(nchf,nchi,ns))**2
                     ab2 = abs(vdon(nchf,nchi,ns)/gkfac)**2
                     partcs(nchp,nchip,ns) = partcs(nchp,nchip,ns) +
     >                  at2 * const

C  The Born partial cross section is infact independent of spin. We use
C  the spin weights to sum up to one, for convenience.
                     BornPCS(nchp,nchip) = BornPCS(nchp,nchip) +
     >                  ab2 * const
C  On occasion there will be a closed channel for nchp < nchpmax
C  This is no big deal, as its partcs = 0.0
                     nchpmax = max(nchp,nent)
                  endif
C  Define the total cross section from the imaginary part of the elastic
C  T matrix. Division by sigma**2 is necessary for charged targets.
                  if (nchp.eq.nchip.and.nchi.eq.nchf) then
                     sigtop(nchip,ns) = sigtop(nchip,ns) - const *
     >                  aimag(ton(nchi,nchi,ns)/sigma(nchi)**2) /
     >                  (pi * gk(1,nchf)+1e-10)
                  endif 
C  End of the "do nchf = 1, nchtop" loop
               end do
C  End of the "if (nchip.le.nent) then" statement
            endif 
            if (no.ne.nopen(ipar)) then
               print*,'NO and NOPEN are not equal',no,nopen(ipar),ipar
c$$$               stop 'NO and NOPEN are not equal'
            endif
         end do !while loop
         call date_and_time(date,time,zone,valuesout)
         print '("cross sections defined at:",a10," diff (secs):",i5)',
     >      time,idiff(valuesin,valuesout)
         valuesin=valuesout
         call update(6)

         write(chnopen,'(i5)') nopen(ipar)
c$$$         write(87) (jlm(nch),nch=1,nopen(ipar))
c$$$         write(87) lent,nymax,((non(nch,l),l=1,jlm(nch)),nch=1,nopen(ipar))
         fmt = '('//chnopen//'i2,"   J, S, Par:",3i4)'
c$$$         write(88,'('//chnopen//'i2,"   J, S, Par:",3i4)')
c$$$         print*,fmt
         write(88,fmt)
     >      (jlm(nch),nch=1,nopen(ipar)),lg,ns,ipar
         if (nchanmax.lt.1000) then
            write(88,'(2i3,3x,10000i3)')
     >         lent,nymax,((non(nch,l),l=1,jlm(nch)),nch=1,nopen(ipar))
         else 
            write(88,'(2i4,3x,10000i4)')
     >         lent,nymax,((non(nch,l),l=1,jlm(nch)),nch=1,nopen(ipar))
         endif       
C  Write the direct only 1st order V to POTL, used for analytic Born
C  subtraction.
C  Write the direct only 1st order V to POTL, used for analytic Born
C  subtraction.
c$$$         write(87) ((dcmplx(wvdon(nchf,nchi)),nchi=1,lent),
c$$$     >      nchf=1,nymax)
         inquire(88,recl=irecl)
c$$$         print*,'lg=',lg,', ns=',ns,', ipar=',ipar, ', recl=', irecl
         if(irecl .lt. 28*lent*nymax) then
            irecl_old = irecl
            irecl = 28*lent*nymax + 10
            print*,' Not sufficient record length for file potl'
            print*,' Set record length to be:', irecl 
c     call update(88)
            close(88)
            open(88,iostat=ios,file=tfile,status='OLD',
     >         access='SEQUENTIAL',form='FORMATTED',
     >         recl=irecl,position='APPEND')
            if(ios .ne. 0) then
               print*,'rest.f: ERROR on opening file ',tfile
               stop
            end if            
         end if
         
         write(chnopen,'(i20)') lent * nymax
         write(88,'(1p,'//chnopen//
     >      '(''('',e12.5,'','',e12.5,'') ''))')
c     >      '(''('',e12.5,'','',e12.5,'')'',1x))')
     >      ((wvdon(nchf,nchi),nchi=1,lent), nchf=1,nymax)

         write(88,'(1p,'//chnopen//
     >      '(''('',e12.5,'','',e12.5,'') ''))')
c    >      '(''('',e12.5,'','',e12.5,'')'',1x))')
     >      ((wton(nchf,nchi),nchi=1,lent),nchf=1,nymax)
c$$$  call openatend(88,tfile)
c$$$         open(88,file=tfile,status='old',fileopt='eof')
c$$$         call update(88)
         close(88)
         open(88,iostat=ios,file=tfile,status='OLD',
     >      access='SEQUENTIAL',form='FORMATTED',
     >      recl=irecl,position='APPEND')
         if(ios .ne. 0) then
            print*,'rest.f: ERROR on opening file ',tfile
            stop
         end if  


         sigt = 0.0
         sum = 0.0
         do nch = 1, nchpmax
            sigt = sigt + partcs(nch,1,ns)
         enddo
            
         
C  Print out the results
         j = mod(lg,100)
         do nchf = 1, nymax !min(nymax,50)
            call getchinfo(nychan(nchf),nchp,
     >         lg,temp,maxpsi,ef,lfa,nfa,lf)
            do nchi = 1, lent
               tcmplx = wton(nchf,nchi)
               if (nychan(nchf).eq.nchf) then
                  tcmplx = tcmplx * ovlpn(nchp)
               else
                  print*,'nychan(nchf),nchf:',nychan(nchf),nchf
               endif 
               tabs = abs(tcmplx)
               if (tabs.gt.0.0) then
                  targ = atan2(aimag(tcmplx),real(tcmplx))
               else
                  targ = 0.0
               endif 
               ctmp = (0.0,0.0)
               if (abs(eigv(nchf)).gt.0.0)
     >            ctmp = log(eigv(nchf))/2.0/(0.0,1.0)
               rdel = real(ctmp)
c$$$               rdel = atan2(aimag(eigv(nchf)),real(eigv(nchf))) / 2.0
               if (ef.le.etot)
c$$$     >            print '(i2,a,i3,i2,1p,e10.3,0p,f9.5,2x,
c$$$     >            1p,2(2(e10.3),1x),0p,f9.4)',
     >            print '(i2,a,2i2,1p,3(2(e10.3),1x),0p,f9.4)',
     >            j,spin(ns),nychan(nchf),nchi,
c$$$     >         tabs,targ,wt2nd(nchf,nchi), wvon(nchf,nchi),ef
     >         wton(nchf,nchi),wt2nd(nchf,nchi), wvon(nchf,nchi),ef
c$$$     >         rdel, nmpow(aimag(ctmp))
c$$$     >         abs(wton(nchf,nchi))
               same = same.and.abs((wvon(nchf,nchi)-wton(nchf,nchi))/
     >            (wton(nchf,nchi)+1e-30)).lt.1e-3
                  if (abs(ef-slowery(1))/abs(ef).lt.0.1) then
                     tcmplx = tcmplx / sqrt(sqrt(ef)) * phaseq(nchp)
c$$$     >                  coulphase(dble(-(zasym+1.0)/sqrt(ef)),lfa)
                     te2eH(lg-lf,lfa) = tcmplx
                  if (lfa.le.lf) then
                     print'("Ef,Lf,lfa,Te2e:",
     >                  f8.3,2i3,1p,2e10.3)',ef*13.6058,lf,lfa,tcmplx
                  else
                     print'("Ef,Lf,lfa,Te2e:",
     >                  f8.3,2i3,1p,2e10.3,i5)',ef*13.6058,lf,lfa,
     >                  tcmplx,nint(100.0*abs(
     >                  (tcmplx-(-1)**ns*te2eH(lg-lfa,lf))/
     >                  (tcmplx+(-1)**ns*te2eH(lg-lfa,lf))))
                  endif                      
               endif 
            end do 
         end do
         if (noprint(2).eqv..false..or.isecond.lt.0) then
            if (abs((sigt-sigtop(1,ns))/(sigtop(1,ns)+1e-30))
     >         .gt.1e-4) then
               print*,'Warning: optical theorem is only satisfied to',
     >            abs((sigt-sigtop(1,ns))/(sigtop(1,ns)+1e-30))
               print*,'Setting SIGTOP to SIGT'
               sigtop(1,ns) = sigt
            endif 
         endif
         eigmod = esum1
c$$$         n = 0
c$$$         do while (eigmod .lt. 0.0)
c$$$            n = n + 1
c$$$            eigmod = eigmod + pi / 2.0
c$$$         enddo 
c$$$         do while (eigmod .gt. pi / 2.0)
c$$$            n = n + 1
c$$$            eigmod = eigmod - pi / 2.0
c$$$         enddo 

         esum(ns) = eigmod
         do nchi = 1, lent
            print '(''N='',i6,'' Rcond:'',1p,e9.2, '', Eigphase sums:'',
     >         0p,2f9.4,i3,'' test:'',10i2)',
     >         nds,rcond(ns),eigmod,esum1,n,(nmpow(err(nchf,nchi)),
     >         nchf=1,min(10,nymax))
         enddo 
         call update(6)
         if (ns.eq.1) then
            sames = .true.
            do nchi = 1, lent
               do nchf = 1, nymax
                  sames = sames.and.
     >               abs((ton(nchf,nchi,0)-ton(nchf,nchi,1))/
     >               (ton(nchf,nchi,0)+1e-4)).lt.1e-3
               end do
            end do
            if (sames) print*,'Singlet and triplet are same'
         end if
C  End the spin loop
      end do 

C  Write out total cross sections
      if (nunit.gt.0) 
     >   call wrtcs(partcs,sigtop,nchpmax,nent,instate,
     >   lg,etot,nsmax,ipar,
     >   ovlpn,nunit,nznuc,zasym,projectile,target,canstop,
     >   BornICS,BornPCS,csfile,esum,rcond,nnbtop,ovlpnn)

C  The following code fails in the CCO part as there the imaginary part is
C  non-zero, but the real part is zero
c$$$      if (same.and.nqm.gt.1) then
c$$$         print*,'setting NQM to 1'
c$$$         nqm = 1
c$$$C  Make sure that the real parts of the on-shell super weights are 0 for UBA.
c$$$         do nch = 1, nchtop
c$$$            wk(1,nch) = 0.0
c$$$         end do 
c$$$      end if 
      end

      function nmpow(err)
      nmpow = 1
c$$$      do while (err * 10 ** nmpow .lt. 1.0.and.nmpow.lt.9)
c$$$         nmpow = nmpow + 1
c$$$      enddo
      return
      end

      subroutine getprod(psi,maxpsi,minchi,chi,l,rmesh,meshr,
     >   nch,ovlp)
      use CI_MODULE   
      include 'par.f'
      dimension psi(maxr), chi(maxr), rmesh(maxr,3)
C common blocks for He
      common /helium/ la(KNM), sa(KNM), lpar(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common /CcoefsE/  E(KNM)
      double precision E, ortint
      dimension ovlpn(nspmax)
      integer na, nam      
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      common /dynamical_C/ Nmax, namax,pnewC
      integer pnewC
ccc      double precision C 
ccc      pointer (pC, C(Nmax,namax,namax))
      
ccc      pC = pnewC
      


      ovlp = 0.0
      ovlp1 = 0.0
      if (nspm.eq.0) then
c$$$         do i = minchi, maxpsi
c$$$            ovlp = ovlp + chi(i) * psi(i) * rmesh(i,3)
c$$$         enddo
c$$$         print*,'ovlp1:',ovlp
         ovlp = sum(chi(minchi:maxpsi)*psi(minchi:maxpsi)*
     >      rmesh(minchi:maxpsi,3))
c$$$         print*,'ovlp2:',ovlp
         ovlp = ovlp * ovlp
      else
c$$$         do np = 1, nspm
c$$$            if (lo(np).eq.L) then
c$$$               if (inc.eq.0.or.np.eq.1.or.sa(np-1).eq.1) then
c$$$                  do i = max(minchi,minf(np)), maxf(np)
c$$$                     ovlpn(np) = ovlpn(np) + chi(i)*fl(i,np)*rmesh(i,3)
c$$$                  enddo
c$$$               endif
c$$$            endif
c$$$         enddo 
c$$$                  

         do 5 jni1 = 1, nam(nch)
            ovlpn(jni1) = 0.0
            ni1 = na(jni1,nch)
            if (lo(ni1).ne.l) go to 5
            do i = max(minchi,minf(ni1)), maxf(ni1)
               ovlpn(jni1) = ovlpn(jni1) + chi(i)*fl(i,ni1)*rmesh(i,3)
            enddo
 5       continue 
               
         do 20 jni1 = 1, nam(nch)
            ni1 = na(jni1,nch)
            if (lo(ni1).ne.l) go to 20
            do 10 jnf1 = 1, nam(nch)
               nf1 = na(jnf1,nch)
               if (lo(nf1).ne.l) go to 10
               do 30 jni2 = 1, nam(nch)
                  ni2 = na(jni2,nch)
                  if (c(nch,jni1,jni2).eq.0d0) go to 30
                  do jnf2 = 1, nam(nch)
                     nf2 = na(jnf2,nch)
                     ovlp = ovlp + c(nch,jni1,jni2) * c(nch,jnf1,jnf2) *
     >                  ortint(nf2,ni2) * ovlpn(jni1) * ovlpn(jnf1) 
c$$$     >                  * ortint(nf1,ni1)
                     ovlp1 = ovlp1 + c(nch,jni1,jni2)*c(nch,jnf1,jnf2) *
     >                  ortint(nf2,ni2) * ortint(nf1,ni1) 
                  enddo
 30            continue
 10         continue
 20      continue
C  A factor of two is necessary for He(1S) when l > 0 as part of the
C  integral disappears for L != l
         if (l.ne.0) ovlp = ovlp * 2.0
c$$$         print*,'He overlap:',ovlp1

c$$$         tmp = 0.0
c$$$         do ni1 = 1, nspm
c$$$            do 100 nf1 = 1, nspm
c$$$               if (lo(nf1).ne.lo(ni1)) go to 100
c$$$               do 300 ni2 = 1, nspm
c$$$                  if (c(nch,ni1,ni2).eq.0d0) go to 300
c$$$                  do 110 nf2 = 1, nspm
c$$$                     if (lo(nf2).ne.lo(ni2)) go to 110
c$$$                     tmp = tmp + c(nch,ni1,ni2) * c(nch,nf1,nf2) *
c$$$     >                  ortint(nf2,ni2) * ortint(nf1,ni1)
c$$$ 110              continue 
c$$$ 300           continue
c$$$ 100        continue
c$$$         enddo
c$$$         if (abs(tmp-1.0).gt.1e-3) print*,tmp
         
      endif
      return
      end
                     
      
C  This routine calculates the average potential of the intermediate state.
      subroutine getugb(psin,istartpsin,istoppsin,nznuc,ugb)
      include 'par.f'
      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
      dimension ff(maxr),temp(maxr),psin(maxr),ugb(maxr),u(maxr)
      data u/maxr*0.0/
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   istartrp(0:ltmax),istoprp(0:ltmax),cntfug(maxr,0:lmax)
      do j=1,maxr
         ugb(j)=0.0
      end do 
      do i=istartpsin,istoppsin
         ff(i)=psin(i)*psin(i)*rmesh(i,3)
      end do
C  Return the integral of FUNA * R<**LTI / R>**(LTI+1) for LTI=0 in TEMP
      call form(ff,istartpsin,istoppsin,rpow1(1,0),rpow2(1,0),
     >   istartrp(0),istoprp(0),meshr,temp,i1,i2)
C  Subtract the nuclear term.
      do i=i1,i2
         temp(i)=temp(i)-float(nznuc)/rmesh(i,1)
      end do 
C  We now define the potential on the XMESH.
      do j=1,meshr
         ugb(j)=temp(j)*2.0
      end do
      return
      end

C  This routine reads the input file 'ccc.in'
      subroutine readin(labot,latop,nabot,natop,lnabtop,nnbtop,energy,de
     >   ,ntstart,ntstop,lstart,lstop,lttop,nbnd,abnd,npsbnd,albnd,alpha
     >   ,ncstates,npot,lpot,lactop,formcut,regcut,expcut,theta,
     >   ifirst,isecond,nold,nq,qcut,rmax,fast,ldw,nk,sk,nznuc,nze,
     >   itail,corep,r0,ninc,linc,ipar,nent,zasym,nunit,ndbl,nps,iborn,
     >   ne2e,lslow,lfast,slowe,enion,enlevel,target,projectile,match,
     >   lpbot,lptop,npbot,nptop,npsp,alphap,luba,speed,
     >   erange,inoenergy,pint,igz,igp,analyticd,packed)
                               
      use ubb_module, only: ubb_max1
      use ql_module, only: maxql1, maxql2, dqlg, qlgmin, qlgmax
      use ftps_module, only: interpol, maxp
      
      include 'par.f'
      common/numericalpotential/ numericalv, lstoppos
      logical fast, match, speed, exists,analyticd,packed,numericalv
      dimension erange(kmax)
      dimension nk(5,0:lmax,2), sk(5,0:lmax,2),
     >   nbnd(0:lmax), abnd(0:lmax),
     >   nabot(0:lamax), natop(0:lamax), alpha(0:lamax), nps(0:lamax),
     >   slowe(ncmax), enions(20), enlevels(20), corep(0:lamax),
     >   r0(0:lamax),npbot(0:lamax),nptop(0:lamax),npsp(0:lamax),
     >   alphap(0:lamax)
      character targets(30)*6,roman(0:30)*10,target*(*),projectile*(*),
     >   infile*80
      data roman/' -',' I','II','III','IV',' V','VI','VII','VIII','IX',
     >   ' X','XI','XII','XIII','XIV','XV','XVI','XVII','XVIII','XIX',
     > 'XX','20+','21+','22+','23+','24+','25+','26+','27+','28+','29+'/
      data targets/'H I','He II','Li I','Be II','B III','C IV','N V',
     >   'O VI','F VII','Ne I','Na I','Mg II','Al III','Si IV',
     >   'P I','S I','Cl I','A VIII','K I','Ca II','Sc','Ti','V ','Cr',
     >   'Mn','Fe','Co','Ni','Cu','Zn'  /
      
C  The following energy ionization levels in Rydbergs are from
C  ATOMIC ENERGY LEVELS, Volume I, C. E. Moore, NSRDS (1971)
C  The data below corresponds to the Z of nearest hydrogen-like target.
C                H      HeII   Li    BeII   BIII   CIV    NV     OVI
      dataenions/13.595,54.400,5.390,18.206,37.920,64.476,97.863,138.08,
C        FVII    Ne     Na    MgII  AlIII SiIV  P    S      Cl    AVIII
     >   185.139,21.559,5.138,15.03,28.44,45.13,11.,10.357,13.01,143.46,
C        K     CaII
     >   4.339,11.87/
C   The data below corresponds to the Z of neutral target.
c$$$  data enions/13.595,24.580,5.390,9.320,8.296,11.264,14.54,13.614,
c$$$     >   17.42,21.559,5.138,7.644,5.984,8.149,11.0,10.357,13.01,15.755,
c$$$     >   4.339,6.111/
C  The energy levels are 0 for the ground state, and the given values below
C  corresponding to the ionization energy above.
C  The data below corresponds to the Z of nearest hydrogen-like target.
C                   H          HeII      Li       BeII     BIII
      data enlevels/109678.758,438889.04,43487.19,146881.7,305931.1,
C        CIV      NV       OVI       FVII    Ne       Na
     >   520177.8,789532.9,1113999.5,1493656,173931.7,41449.65,
C        MgII      AlIII     SiIV     P       S       Cl       AVIII
     >   121267.61,229453.99,364097.7,88560.0,83559.3,104991.0,1157400.,
C        K        CaII
     >   35009.78,95748.0/
C  The data below corresponds to the Z of neutral target.
c$$$      data enlevels/109678.758,198305.0,43487.19,75192.29,66930.0,
c$$$     >   90878.3,117345.0,109836.7,140553.5,173931.7,41449.65,
c$$$     >   61669.14,48279.16,65743.0,88560.0,83559.3,104991.0,127109.9,
c$$$     >   35009.78,49304.8/

      do la = 0, lamax
         corep(la) = 0.0
         r0(la) = 1.0
      enddo

      inquire(file='ccc.in',exist=exists)
      if (exists) then
         nin = 3
         open(nin,file='ccc.in')
         print*,'ccc.in found'
      else
         num_args = command_argument_count()
         if (num_args .eq. 1) then
            call get_command_argument(1,infile)
            nin=3
            open(nin,file=infile)
            print*,'ccc.in not found;will use: '//infile
         else
            nin = 5
            print*,'ccc.in not found;will use standard input'
         endif 
      endif 
      read(nin,*,end=13,err=13) energy,de,nznuc,zasym,nze,ninc,linc
      print '(''energy,de,nznuc,zasym,nze,ninc,linc:'',
     >   f8.3,f7.4,i7,f7.1,i7,i4,i3)',
     >   energy,de,nznuc,zasym,nze,ninc,linc
      lpbot = 0
      lptop = -1
      do l = 0, lamax
         nabot(l) = l + 1
         natop(l) = l
         npbot(l) = l + 1
         nptop(l) = l
      enddo 
      if (nze.eq.-1) then
         projectile = 'electron'
      else if (nze.eq.1) then
         projectile = 'positron'
C  Read positronium information
         read(nin,*) lpbot,lptop,(npbot(lp),nptop(lp),lp=lpbot,lptop)
         print '(''lpbot, lptop, npbot(l), nptop(l):    '',2i7,
     >   99(i4,i3))',lpbot, lptop,(npbot(l), nptop(l), l = lpbot, lptop)
         read(nin,*) (npsp(l),alphap(l), l = lpbot, lptop)
         print '(''npsp(l), alphap(l):                  '',
     >   99(i3,f5.2))',(npsp(l), alphap(l), l = lpbot, lptop)
         read(nin,*) igz, igp, analyticd, numericalv,lstoppos,
     1        interpol, maxp, ubb_max1, maxql1
         print'(''igz, igp, analyticd, numericalv:  '',2i7,2l7,i4)',
     >      igz, igp, analyticd, numericalv, lstoppos
!C     Andrey: interpolation parameters: 
         print'(''interpol, maxp, UBB_MAX, maxql1:    '', 1l7,3i7)',
     1        interpol, maxp, ubb_max1, maxql1
         dqlg =  (qlgmax-qlgmin)/dble(maxql1-1)
         maxql2 = maxql1
      else if (nze.eq.6) then
         projectile = ' pos6+'
      else if (nze.eq.2) then
         projectile = ' pos2+'
      else if (nze.eq.0) then
         projectile = 'photon'
         nze = -1
      else 
         print*,'NZE should be -1 for electron scattering'
         print*,'NZE should be +1 for positron scattering'
         print*,'NZE should be  0 for  photon  scattering'
         print*,'Here NZE:',nze
         stop 'Wrong value for NZE'
      end if

      nzasym = nint(zasym)
      n = nznuc - nzasym
      print*,"n, nznuc, nzasym:", n, nznuc, nzasym
c$$$      if (projectile.eq.'photon') n = n + 1
      target = ' ???? '
      if (nznuc.le.30) then
         target = targets(nznuc)
         enlevel = enlevels(nznuc)
         enion = enions(nznuc)
         target = targets(nznuc)(1:2)//roman(nint(zasym+1))
      endif 
      lactop = 0
      if (n.eq.1) then  ! H-like atoms and ions
         enion = 13.6058 * nznuc ** 2
         enlevel = enion * 8067.8
         if (nznuc.eq.30) then
            target = 'Zn XXX'
         elseif (nznuc.eq.26) then
            target = 'FeXXVI'
         endif 
      else if (n.eq.2) then ! He-like atoms and ions
         if (nzasym.eq.0) then ! He I
            enlevel = 198305.0
            enion  = 24.580
         elseif (nzasym.eq.8) then ! Ne IX
            enlevel = 9644957  ! cm-1
            enion = 1195.822   ! eV
         elseif (nzasym.eq.24) then ! FeXXV
            target = 'FeXXV'
            enlevel = 71204137.0 ! cm-1
            enion = 8828.1875    ! eV
         elseif (nzasym.eq.90) then
            target = 'U90+'
         endif 
      else if (n.eq.3) then ! Li-like atoms and ions
         if (nzasym.eq.7) then ! Ne VIII
            enlevel = 1928462.0
            enion = 239.0989
         elseif (nzasym.eq.23) then !FeXXIV
            target='FeXXIV'
            enlevel = 16500160.0
            enion = 2045.759
         endif 
      else if (n.eq.4) then ! Be-like atoms and ions
         if (nzasym.eq.0) then ! Be I
            enlevel = 75192.29
            enion  = 9.32
            do l = 0, lamax
               corep(l) = 0.05224
               r0(l) = 0.98
            enddo
            r0(1) = 0.9
         else if (nzasym.eq.2) then ! C III
            enlevel = 386241.0 ! cm-1
            enion = 47.8878 ! eV
         else if (nzasym.eq.4) then ! O V
            enlevel = 918702.0
            enion = 113.873
         else if (nzasym.eq.6) then ! Ne VII
            enlevel = 1671792.
            enion = 207.27
         else if (nzasym.eq.22) then ! Fe XXIII
            target = 'FeXXIII'
            enlevel = 15731000
            enion = 1950.4
         endif
      else if (n.eq.5) then ! B-like atoms and ions
         lactop = 0
         nabot(0) = 2
      else if (n.eq.6) then ! C-like atoms and ions
         lactop = 0
         nabot(0) = 3
      else if (n.eq.7) then ! N-like atoms and ions
         lactop = 0
         nabot(0) = 3
      else if (n.eq.8) then ! O-like atoms and ions
         lactop = 0
         nabot(0) = 3
      else if (n.eq.9) then ! F-like atoms and ions
         lactop = 1
         nabot(0) = 2
         nabot(1) = 2
         if (nzasym.eq.1) then ! Ne II
            enlevel = 330391.0
            enion = 40.96328
         endif 
      else if (n.eq.10) then    ! Ne-like atoms and ions
         lactop = 1
         nabot(0) = 2
         nabot(1) = 3
         enlevel =  173929.75   ! cm-1
         enion = 21.5645        ! eV         
      else if (n.eq.11) then ! Na-like atoms and ions
         lactop = 1
         nabot(0) = 3
         nabot(1) = 3
         if (nzasym.eq.0) then ! Na I
            do l = 0, lamax
               corep(l) = 0.99
               r0(l) = 1.439
            enddo
         else if (nzasym.eq.1) then ! MgII
            do l = 0, lamax
               corep(l) = 0.48
               r0(l) = 1.3
            enddo
            r0(0) = 1.1
            r0(1) = 1.1
         else if (nzasym.eq.5) then ! S VI
            enlevel = 710194.7
            enion = 88.05295
         endif 
      elseif (n.eq.12) then ! Mg-like atoms and ions
         lactop = 1
         nabot(0) = 3
         nabot(1) = 3         
         if (nzasym.eq.0) then ! Mg I
c$$$            enlevel = 61669.14
c$$$            enion  = 7.644
c$$$            do l = 0, lamax
c$$$               corep(l) = 0.48
c$$$               r0(l) = 1.115
c$$$            enddo
            print*,'Will be using Mg II parameters'
            do l = 0, lamax
               corep(l) = 0.4814
               r0(l) = 1.3
            enddo
            r0(0) = 1.1
            r0(1) = 1.1
            inquire(file='JMitroy_par',exist=exists)
            if (exists) then
               print*
               print*,"Set one-electron pol.pot parameters to ",
     >              "J.Mitroy choice: PRA 78(2008)012715"
               print*
               do l = 0, lamax
                  corep(l) = 0.4814
                  r0(l) = 1.361
               enddo
               r0(0) = 1.1795
               r0(1) = 1.302
               r0(2) = 1.442
               r0(3) = 1.52               
            endif

         elseif (nzasym.eq.4) then ! S V
            enlevel = 585514.1
            enion = 72.59452
         endif 
      else if (n.eq.17) then ! Cl-like atoms and ions
         lactop = 1
         nabot(0) = 3
         nabot(1) = 2
         if (nzasym.eq.1) then ! Ar II
            enlevel = 330391.0
            enion = 40.96328
         endif 
      elseif (n.eq.18) then ! Ar-like atoms and ions
         lactop = 1
         nabot(0) = 3
         nabot(1) = 4         
         if (nzasym.eq.0) then ! Ar I
            do l = 0, lamax
               corep(l) = 0.0
               r0(l) = 2.0
            enddo
            elevel = 127109.8
            enion = 15.7596
            target = 'Ar I'
         endif 
      elseif (n.eq.19) then ! K-like atoms and ions
         lactop = 1
         nabot(0) = 4
         nabot(1) = 4         
         if (nzasym.eq.0) then ! K I
            do l = 0, lamax
               corep(l) = 5.354
               r0(l) = 3.0
            enddo
            r0(0) = 2.0
            r0(1) = 2.0
            r0(2) = 2.5
c$$$            corep(3:lamax) = 0.0
         endif 
      elseif (n.eq.20) then ! Ca-like atoms and ions
         lactop = 1
         nabot(0) = 4
         nabot(1) = 4
         if (nzasym.eq.0) then
            elevel = 49304.8
            enion  = 6.111
            do l = 0, lamax
               corep(l) = 3.26
               r0(l) = 1.953
            enddo
            r0(0) = 1.71
            r0(1) = 1.71
         endif 
      elseif (n.eq.29) then ! Cu-like atoms and ions
         lactop = 2
         nabot(0) = 4
         nabot(1) = 4
         nabot(2) = 4
         if (nzasym.eq.0) then ! Cu I
            target = 'Cu I'
            enlevel = 62317.2
            enion  = 7.724
            do l = 0, lamax
               corep(l) = 5.0
               r0(l) = 2.0
            enddo
            corep(0) = 8.0
            corep(1) = 8.0            
            r0(0) = 2.3
            r0(1) = 2.5
         endif 
         if (nzasym.eq.2) then ! Ga III
            target = 'Ga III'
            print*, "Target: ",target
            enlevel = 247820.
            enion  = 30.7258 
            do l = 0, lamax
               corep(l) = 2.
               r0(l) = 2.0
            enddo
            r0(0) = 1.565
            r0(1) = 2.
         endif 

      elseif (n.eq.30) then ! Zn-like atoms and ions
         lactop = 2
         nabot(0) = 4
         nabot(1) = 4
         nabot(2) = 4
         if (nzasym.eq.0) then ! Zn I
            target = 'Zn I'
            enlevel = 75766.8
            enion  = 9.391
c            do l = 0, lamax
c               corep(l) = 2.6
c               r0(l) = 2.0
c            enddo
c            corep(0) = 8.0
c            corep(1) = 8.0            
c            r0(0) = 2.3
c            r0(1) = 2.5
         endif
         if (nzasym.eq.1) then ! Ga II
            target = 'Ga II'
            enlevel = 165465.8
            enion  = 20.51514
c            do l = 0, lamax
c               corep(l) = 2.6
c               r0(l) = 2.0
c            enddo
c            corep(0) = 8.0
c            corep(1) = 8.0            
c            r0(0) = 2.3
c            r0(1) = 2.5
         endif
 
      elseif (n.eq.31) then ! Ga-like atoms and ions
         lactop = 2
         nabot(0) = 5
         nabot(1) = 4
         nabot(2) = 4
         if (nzasym.eq.0) then ! Ga I
            target = 'Ga I'
            enlevel = 48387.634
            enion  = 5.9993
c            do l = 0, lamax
c               corep(l) = 2.6
c               r0(l) = 2.0
c            enddo
c            corep(0) = 8.0
c            corep(1) = 8.0            
c            r0(0) = 2.3
c            r0(1) = 2.5
         endif 
      elseif (n.eq.36) then ! Kr-like atoms and ions
         lactop = 2
         nabot(0) = 5
         nabot(1) = 4
         nabot(2) = 4

         enlevel = 112914.433
         enion = 13.9996
         target = 'Kr I'
      elseif (n.eq.37) then ! Rb-like atoms and ions
         lactop = 2
         nabot(0) = 5
         nabot(1) = 5
         nabot(2) = 4
         if (nzasym.eq.0) then ! Rb I
            target = 'Rb I'
            enlevel = 33691.0
            enion  = 4.176
            do l = 0, lamax
               corep(l) = 9.0
               r0(l) = 3.0
            enddo
            r0(0) = 2.05
            r0(1) = 2.25
         elseif (nzasym.eq.1) then ! Sr II
            target = 'Sr II'
            enlevel = 88964.0
            enion  = 11.027
            do l = 0, lamax
               corep(l) = 6.8
               r0(l) = 2.9
            enddo
            r0(0) = 1.95
            r0(1) = 2.25
            r0(2) = 2.75
         endif 
      elseif (n.eq.38) then ! Sr-like atoms and ions
         lactop = 2
         nabot(0) = 5
         nabot(1) = 5
         nabot(2) = 4
         if (nzasym.eq.0) then
            target = 'Sr I'
            enlevel = 45925.6
            enion  = 5.692
            do l = 0, lamax
               corep(l) = 6.8
               r0(l) = 2.9
            enddo
            r0(0) = 1.95
            r0(1) = 2.25
            r0(2) = 2.75
         endif 
      elseif (n.eq.47) then ! Ag-like atoms and ions
         lactop = 2
         nabot(0) = 5
         nabot(1) = 5
         nabot(2) = 4
         if (nzasym.eq.0) then
            target = 'Ag I'
            enion  = 7.576234
            enlevel = 8065.5447*enion
            do l = 0, lamax
               corep(l) = 7.66
               r0(l) = 2.5
            enddo
            r0(0) = 1.91
            r0(1) = 2.0
         endif 
      elseif (n.eq.48) then ! Cd-like atoms and ions
         lactop = 2
         nabot(0) = 5
         nabot(1) = 5
         nabot(2) = 4
         if (nzasym.eq.0) then
            target = 'Cd I'
            enlevel = 8065.5447*8.994
            enion  = 8.994
            do l = 0, lamax
               corep(l) = 4.971
               r0(l) = 1.9
            enddo
            r0(0) = 1.225
            r0(1) = 1.45
            r0(2) = 1.8
         endif 
!  Cd II  ionization energy 16.91 eV
      elseif (n.eq.54) then     ! Xe-like atoms and ions
         lactop = 2
         nabot(0) = 6
         nabot(1) = 5
         nabot(2) = 5
         target = 'Xe I'

         enlevel = 97833.787 
         enion = 12.1298

      elseif (n.eq.55) then ! Cs-like atoms and ions
         lactop = 2
         nabot(0) = 6
         nabot(1) = 6
         nabot(2) = 5
         if (nzasym.eq.0) then
            target = 'Cs I'
            enlevel = 31406.710
            enion  = 3.8937
            do l = 0, lamax
               corep(l) = 15.644
               r0(l)  = 4.0
            enddo
            r0(0) = 2.1
            r0(1) = 2.5
            r0(2) = 3.7
            r0(3) = 2.5
         elseif (nzasym.eq.1) then  
            target = 'Ba II'
            enlevel = 80686.87
            enion  = 10.001
            do l = 0, lamax
               corep(l) = 11.0
               r0(l)  = 3.4
            enddo
            r0(0) = 1.804
            r0(1) = 2.2
            r0(2)  = 3.4
            r0(3)  = 2.47
         endif 
      elseif (n.eq.56) then ! Ba-like atoms and ions
         lactop = 2
         nabot(0) = 6
         nabot(1) = 6
         nabot(2) = 5
         if (nzasym.eq.0) then
            target = 'Ba I'
            enlevel = 42032.4
            enion  = 5.210
            do l = 0, lamax
               corep(l) = 11.0
               r0(l)  = 3.4
            enddo
            r0(0) = 1.804
            r0(1) = 2.2
            r0(2)  = 3.4
            r0(3)  = 2.47
         endif 
      elseif (n.eq.69) then ! Tm-like atoms and ions
         lactop = 3
         nabot(0) = 6
         nabot(1) = 6
         nabot(2) = 5
         nabot(3) = 5
         if (nzasym.eq.0) then
            target = 'Tm I'
            enlevel = 0.0 
            enion  = 0.0 
         else if (nzasym.eq.1) then
            lactop = 3
            target = 'Yb II'
            enlevel = 98207
            enion  = 12.176
            do l = 0, lamax
               corep(l) = 28.4 
               r0(l)  = 3.0
            enddo
         endif 
      elseif (n.eq.70) then ! Yb-like atoms and ions
         lactop = 3
         nabot(0) = 6
         nabot(1) = 6
         nabot(2) = 5
         nabot(3) = 5
         if (nzasym.eq.0) then
            target = 'Yb I'
            enlevel = 50443.2
            enion  = 6.2542 
            print*,'Yb I, ionization energy 6.2542 eV'
            do l = 0, lamax
               corep(l) = 28.4 
               r0(l)  = 3.65   ! spd core:  2.86
            enddo
            r0(0)  = 3.0
            r0(1)  = 3.65
            r0(2)  = 3.65       ! 
         endif 
      elseif (n.eq.79) then ! Au-like atoms and ions
         lactop = 3
         nabot(0) = 6
         nabot(1) = 6
         nabot(2) = 6
         nabot(3) = 5
         if (nzasym.eq.0) then
            target = 'Au I'
            enlevel = 74410
            enion  = 9.2256
            do l = 0, lamax
               corep(l) = 12.7 
               r0(l)  = 2.1
            enddo
            r0(0) = 2.0
            print*, '********** Use pol. pot. for Au'
         else if (nzasym.eq.1) then
            lactop = 3
            nabot(0) = 6
            target = 'Hg II'
            enlevel = 151280.
            enion  = 18.751
         endif 
      elseif (n.eq.80) then ! Hg-like atoms and ions
         lactop = 3
         nabot(0) = 6
         nabot(1) = 6
         nabot(2) = 6
         nabot(3) = 5
         if (nzasym.eq.0) then
            target = 'Hg I'
            enlevel = 84184.1
            enion  = 10.43  
            do l = 0, lamax
               corep(l) = 8.4 
               r0(l)  = 2.3
            enddo
            r0(0) = 1.632
            r0(1) = 2.4
            r0(2)  = 2.8
            r0(3)  = 2.7
            r0(4)  = 2.65
         endif 
      else
         stop "Don't know how to calculate scattering on this target"
      endif
      
      read(nin,*) labot, latop, (nabot(l), natop(l), l = labot, latop)
      print '(''labot, latop, nabot(l), natop(l):    '',2i7,99(i4,i3))',
     >   labot, latop, (nabot(l), natop(l), l = labot, latop)

      if (latop.lt.lptop) stop 'expect LPTOP <= LATOP'
      if (latop.gt.lamax) then
         print*,'Must have LATOP <= LAMAX', latop, lamax
         stop 'Must have LATOP <= LAMAX'
      endif
      

      do l = labot, latop
         if (nabot(l).le.l) then
            print*,'Must have NABOT(L) >= L + 1, NABOT(L), L:',
     >         nabot(l),l
            stop 'Must have NABOT(L) >= L + 1'
         endif
         if (natop(l).gt.nnmax) then
            print*,'Must have NATOP(L) <= NNMAX',natop(l),nnmax
            stop 'Must have NATOP(L) <= NNMAX'
         endif
      end do

      read(nin,*) ntst,nunit,lnabtop,nnbtop,lttop,ncstates
      print '(''ntst,nunit,lnabtop,nnbtop,lttop,ncstates:'',
     >   i3,5i7)',ntst,nunit,lnabtop,nnbtop,lttop,ncstates
      ntstart = ntst
      ntstop = ntst

      read(nin,*) lstart,lstop,ipar,nent,iborn
      print '(''lstart, lstop, ipar, nent, iborn:    '',5i7)',
     >   lstart,lstop,ipar,nent,iborn
      if (lstop+latop.gt.lmax) stop 'Can not have LSTOP+LATOP > LMAX'

      read(nin,*) ndumm,luba,(nps(l),alpha(l),l=labot,latop)
      print '(''ndumm,luba,nps(l),alpha(l):         '',
     >   i4,i3,94(i3,f7.2))',
     >   ndumm,luba,(nps(l),alpha(l),l=labot,latop)
      if (alpha(latop).eq.0.0) stop 'ALPHA(LATOP) can''t be zero'
      
      read(nin,*) npot,lpot,ldw,npsbnd,albnd
      if (npot.eq.0) ldw = -1
      if (ldw.lt.0) then
         npot = 0
         lpot = 0
      endif 
      print '(''npot,lpot,ldw,npsbnd,albnd:          '',4i7,f7.2)',
     >   npot,lpot,ldw,npsbnd,albnd

C  FORMCUT cuts the form factors in FORM
C  REGCUT determines the minimum value at which regular solutions start
C  EXPCUT determines the smallest value of functions containing EXP fall off
      read(nin,*) formcut,regcut,expcut,gamm,rc
      print '(''formcut,regcut,expcut,gamma,rc:       '',1p,3e7.0,0p,
     >   2f7.3)', formcut,regcut,expcut,gamm,rc
      if (gamm.ge.0.0) then
         do la = 0, lamax
            corep(la) = gamm
            r0(la) = rc
            if (r0(la).eq.0.0) stop 'R0(la) can''t be 0.0'
         enddo
      endif 

      read(nin,*) ifirst,isecond,nold,itail,theta
      if (nze.ge.1) then
         if (ifirst.gt.0) ifirst = 0
         if (isecond.gt.0) isecond = 0
      end if 
      print '(''ifirst,isecond,nold,itail,theta:     '',4i7,f7.2)',
     >   ifirst,isecond,nold,itail,theta
      if (ndumm.eq.0.and.isecond.ge.0) then
         print*,'ISECOND >= 0 is inconsistent with NDUMM = 0'
         stop 'ISECOND >= 0 is inconsistent with NDUMM = 0'
      endif 

      slowe(:) = 0.0
      read(nin,*) ne2e, lslow, lfast, (slowe(n), n=1, abs(ne2e))
C  There is some bug on the SS10 machines that causes malloc to get a SEGV
C  if e2e is stopped before the CC calculation is completed
      lslow = max(latop,lslow)
      lfast = max(lstop,lfast)
      print '(''ne2e,lslow,lfast, slowe(n),n=1,|ne2e|'',3i7,100f7.3)',
     >   ne2e, lslow, lfast, (slowe(n), n=1, abs(ne2e))
      
      read(nin,*) nq,qcut,rmax,ndbl,fast,match,packed
      print '(''nq,qcut,rmax,ndbl,fast,match,packed: '',
     >   i7,2f7.1,i7,3l3)',nq,qcut,rmax,ndbl,fast,match,packed

      
      inquire(file='speed.in',exist=exists)
      if (exists) then
         ispeed = 2
         open(4,file='speed.in')
         read(4,*) speed
         read(4,*) nintsp, pint
         isp=0
         do n=1,nintsp
            read(4,*) estartsp, desp, nesp
            do ne=1,nesp
               isp=isp+1
               erange(isp)=estartsp+(ne-1)*desp
               print*,isp,erange(isp),'eV,',
     >          erange(isp)/13.6058,'Ryd'
            enddo
         enddo
         if (isp.gt.kmax) stop 'INCREASE DIMENSION OF ERANGE'
         inoenergy=isp
      else
         speed = .false.
         ispeed = 1
         inoenergy = 1
         erange(1) = energy
      endif
      
      do is = 1, ispeed
         lp = 0
         mint = 4
         do while (lp.le.lstop+latop+2) !+2 is for two-electron targets  
            if (is.eq.1) then
               read(nin-1+is,*)l,nbnd(l),
     >            (nk(j,l,is),sk(j,l,is),j=1,mint)
            else
               read(nin-1+is,*) l, ntmp,(nk(j,l,is),sk(j,l,is),j=1,mint)
            endif 
            if (l.gt.lmax) stop 'L > LMAX'
c$$$            if (albnd.lt.0.0) then
c$$$               abnd(l) = - albnd / float(l + 1)
c$$$            else
               abnd(l) = albnd
c$$$            endif 
            do while (lp.lt.l.and.lp.le.lstop+lamax)
               nbnd(lp) = max(0,nbnd(l)-lp)
c$$$            if (zasym.eq.0.0.and.lp.gt.ldw) nbnd(lp) = 0
c$$$               if (albnd.lt.0) then
c$$$                  abnd(lp) = - albnd / float(lp + 1)
c$$$               else
                  abnd(lp) =  albnd
c$$$               endif 
               do j = 1, mint
                  nk(j,lp,is) = nk(j,l,is)
                  sk(j,lp,is) = sk(j,l,is)
               enddo
               print '(''l,nbnd(l),(nk(j,l),sk(j,l),j=1,4)'',
     >            2i3,4(i3,f5.2))',
     >            lp,nbnd(lp),(nk(j,lp,is),sk(j,lp,is),j=1,mint)
               lp = lp + 1
            enddo
            nbnd(l) = max(0,nbnd(l)-l)
c$$$         if (zasym.eq.0.0.and.l.gt.ldw) nbnd(l) = 0
            lp = l + 1
            print '(''l,nbnd(l),(nk(j,l),sk(j,l),j=1,4)'',
     >         2i3,4(i3,f5.2))',
     >         l,nbnd(l),(nk(j,l,is),sk(j,l,is),j=1,mint)
            if (mod(nk(4,l,is),2).ne.0) stop 'NK(4,L,IS) must be even'
            nktot = nk(1,l,is)+nk(2,l,is)+nk(3,l,is)+max(0,nk(4,l,is))
     >         + nbnd(l)
            if (nktot.gt.kmax) stop 'KMAX too small in par.f'
         enddo 
         CLOSE(nin-1+is)
      enddo


      if (lttop.gt.ltmax) then
         print*,'LTTOP > LTMAX',lttop,ltmax
         print*,'Recompile with larger LTMAX'
         stop 'ABORTED'
      end if
      if (nq.gt.kmax) then
         print*,'Increase KMAX in par.f'
         stop 'Increase KMAX in par.f'
      endif 
      return
 13   stop 'no file ccc.in'
      end
      
      SUBROUTINE ARRANGE
C
C  this subroutine sets up the faclog array.
C
      INCLUDE 'par.f'
      double precision faclog
      COMMON /flogs/ FACLOG(1000)
      FACLOG(1)=0d0
      FACLOG(2)=0d0
      DO 10 N=3,1000
   10 FACLOG(N)=FACLOG(N-1)+log(dfloat(n-1))

      return
      end
      
C  This routine defines the continuum energy mesh CE. The integration range
C  goes from zero to infinity. The mesh is divided into three intervals
C  with the first two intervals having twice as many points as the last.
C  The first interval goes from zero to
C  CE0, the energy where the Green's function energy is zero, and the other two
C  after. CE is in Rydbergs. The integration is of the form 
C  dk F(k*k) = dE F(E)/2/sqrt(E). Note: the Jacobian k**2 has been cancelled
C  with the k appearing in the denominator of the two coulomb functions that
C  always appear in the integrand.
      subroutine egrid(npotgf,etot)
      include 'par.f'
      common/cont/ ce(ncmax),cint(ncmax),ncstates,energy
      real*8 t(ncmax),wts(ncmax),wf(2*ncmax),dstart,dstop
      integer iwf(2*ncmax)
      data ry/13.6058/
      if (ncstates.gt.ncmax) then
         print*,'NCSTATES cannot be greater than NCMAX'
         stop 'Job aborted in EGRID'
      end if
      if (mod(ncstates,5).ne.0) then
         ncstates=ncstates-mod(ncstates,5)
         print*,'NCSTATES not a multiple of five. Reset to:',ncstates
      end if
C  The Green's function energy is zero at CE0 = total energy.
      ce0 = etot
C  For incident energy less than 1 Ry. make sure CE0 is still +ve.
C  Requires further investigation.      
c$$$      if (ce0.lt.0.0) ce0 = 1.0
      if (ce0.lt.0.0) ce0 = 0.5
      
C  NT is the number of Gaussian points in each interval, NWF and NIWF are
C  dimensions of work arrays of routine CGQF.
      nt=2*ncstates/5
      nwf=2*ncmax
      niwf=2*ncmax

      if (npotgf.lt.0) then
C  First interval goes from zero to CE0.
         dstart=0d0
         dstop=ce0
C  Get plane Gaussian quadratures on the interval DSTART to DSTOP
         call cgqf(nt,t,wts,1,0d0,0d0,dstart,dstop,0,nwf,wf,niwf,iwf,i)
         print*,
     >      'Kn for which the Green''s function energy is positive'
         print*,'First interval'
         ns=0
         DO J=1,nt
            ns=ns+1
            CE(ns)= t(j)
            cint(ns)=wts(j)/2.0/sqrt(t(j))
            print*,ns,sqrt(ce(ns)),cint(ns)
         end do
      
C  Now for quadrature points that will lead to negative EPROJ, second interval.
         print*,
     >      'Kn for which the Green''s function energy is negative'
         print*,'Second interval'
         dstart=dstop
c$$$      print*,'enter first N please:'
c$$$      read*,n
C  Set the following N to 100000000 to get a zero energy for the G.F.s
         n=4
         n=1000000000
         t(1)=dstart+(1.0/n/n+1.0/(n+1)/(n+1))/2.0
c$$$      print*,'enter second N please:'
c$$$      read*,n
         n=7
c$$$         n=3
         t(2)=dstart+(1.0/n/n+1.0/(n+1)/(n+1))/2.0
c$$$      print*,'enter third N please:'
c$$$      read*,n
         n=1
         t(3)=dstart+(1.0/n/n+1.0/(n+1)/(n+1))/2.0
         t(3)=2.0*t(2)-t(1)
         do j=1,3
            ns=ns+1
            ce(ns)=t(j)
            cint(ns)=(t(2)-t(1))/3.0/2.0/sqrt(ce(ns))
            if (j.eq.2) cint(ns)=cint(ns)*4.0
            print*,ns,sqrt(ce(ns)),cint(ns)
         end do
         m=nt-3
         dstart=ce(ns)
         dstop=ce(ns)+2.2
c$$$         t(4)=dstart+1.3
c$$$      print*,'enter dstop please:'
c$$$      read*,dstop
c$$$         dstop=dstart+1.5
c$$$         call cliqf(4,t,wts,1,0d0,0d0,dstart,dstop,0,nwf,wf,niwf,iwf,i)
c$$$         DO J=1,4
c$$$            ns=ns+1
c$$$            CE(ns)= t(j)
c$$$            cint(ns)=wts(j)/2.0/sqrt(ce(ns))
c$$$            print*,ns,sqrt(ce(ns)),cint(ns)
c$$$         end do
c$$$         dstart=dstop
c$$$         dstop=dstart+1.0
c$$$         m=nt-4
         call cgqf(m,t,wts,1,0d0,0d0,sqrt(dstart),sqrt(dstop),0,nwf,wf,
     >      niwf,iwf,ier)
         DO J=1,m
            ns=ns+1
            CE(ns)= t(j)**2
            cint(ns)=wts(j)
            print*,ns,sqrt(ce(ns)),cint(ns)
         end do
      else
C  First interval goes from zero to CE0
         nt=ncstates*2/5
         ns=0
         dstart=0d0
         dstop=ce0
C  Get plane Gaussian quadratures on the interval DSTART to DSTOP
         call cgqf(nt,t,wts,1,0d0,0d0,dstart,dstop,0,nwf,wf,niwf,iwf,i)
         if (etot.gt.0.0) then
            print*,
     >         'Kn for which the Green''s function energy is positive'
         else 
            print*,
     >         'Kn for which the Green''s function energy is negative'
         endif 
         print*,'First interval'
         DO J=1,nt
            ns=ns+1
            CE(ns)= t(j)
            cint(ns)=wts(j)/2.0/sqrt(t(j))
            print*,ns,sqrt(ce(ns)),cint(ns)
         end do
         print*,
     >      'Kn for which the Green''s function energy is negative'
         print*,'Second interval'
         dstart=dstop
         dstop=ce0 + 1.0
C  Get plane Gaussian quadratures on the interval DSTART to DSTOP
         call cgqf(nt,t,wts,1,0d0,0d0,dstart,dstop,0,nwf,wf,niwf,iwf,i)
         DO J=1,nt
            ns=ns+1
            CE(ns)= t(j)
            cint(ns)=wts(j)/2.0/sqrt(t(j))
            print*,ns,sqrt(ce(ns)),cint(ns)
         end do
      end if       
C  For the last interval we assume that F(E) falls off as 1/E**P, where
C  by trial and error, I chose P to be two.
      print*,'Third interval'
      p=2.0
      dstart=0d0
      dstop=dstop**((1.0-2.0*p)/2.0)
      nt=ncstates/5
      call cgqf(nt,t,wts,1,0d0,0d0,dstart,dstop,0,nwf,wf,niwf,iwf,ier)
      nt=4*nt
      do j=nt+1,ncstates
         jj=ncstates-j+1+nt
         ce(jj)=t(j-nt)**(2.0/(1.0-2.0*p))
         cint(jj)=wts(j-nt)*ce(jj)**p/(2.0*p-1.0)
         print*,jj,sqrt(ce(jj)),cint(jj)
      end do
      do j = 1, ncstates
         cint(j) = cint(j) * ce(j)
      enddo 
      RETURN
      END

      subroutine setpow(lstop,lttop)
C  this subroutine sets up the arrays RPOW1=R**LNA and RPOW2=1/R**(LNA+1)
C  which are used in FORM, as well as CNTFUG=L(L+1)/XMESH(J)**2 which 
C  is used in Numerov integration.
      include 'par.f'
      common /debye/ dbyexists, dblp, rmudeb, zdby
      common/meshrr/ meshr,rmesh(maxr,3)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   istartrp(0:ltmax),istoprp(0:ltmax),cntfug(maxr,0:lmax)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast, dbyexists
      complex*16 :: XX, FC(0:ltmax), GC(0:ltmax), FCP(0:ltmax),
     >     GCP(0:ltmax), SIG(0:ltmax), CMPLXI=(0d0,1d0), ZLMIN, ETA1
      integer ::  NL, KFN, IFAIL, MODE1
      real*8 :: dblp_dp, FCR(0:ltmax),GCR(0:ltmax)  
      real*8 :: gfact
 

      do j=1,meshr
         cntfug(j,0)=0.0
         x=rmesh(j,1)
         x2=x*x
         cntfug(j,1)=2.0/x2
      end do
      do l=2,min(lstop+lttop,lmax)
         do j=1,meshr
            cntfug(j,l)=cntfug(j,l-1)*float(l+1)/float(l-1)
         end do
      end do

      istartrp(0)=1
      istoprp(0)=meshr

      rpow1(:,:)=0
      rpow2(:,:)=0

                 
      if (dbyexists) then
C     This is the partial wave expansion of the electron-electron
C     potential in a Debye plasma.

C     FC Spherical Bessel. FC(1) lna zero order
C     GC Hankel function first kind. GC(1) lna zero order
         
         istartrp(:)=1
         istoprp(:)=meshr
         dblp_dp = dblp

         ETA1=(0,0)    !Set to zero if extraction hankel and spherical bessel
         MODE1=12  !Extracting hankel function first order subroutine
         KFN=1     !Extracting spherical bessel function subroutine
         IFAIL=2   !Testing error subroutine. If okay IFAIL=0
         ZLMIN=(0,0)   !Starting order of l
         NL=lttop  !Max order of l

         do i=1, meshr
        
            XX=CMPLXI*rmesh(i,1)/dblp_dp

            call COULCC(XX,ETA1,ZLMIN,NL,FC,GC,FCP,GCP,
     >           SIG,MODE1,KFN,IFAIL)
            Debyelna: do lna=0,lttop
               
               FCR(lna)=  abs(real(FC(lna) / (CMPLXI ** (lna)))) 
               GCR(lna)=  abs(real((GC(lna) * CMPLXI ** (lna)))) 
               
               rpow1(i,lna)=  exp(gfact(2*lna+1))*
     >              FCR(lna)*dblp_dp**lna 
               rpow2(i,lna)= GCR(lna) / (exp(gfact(2 * lna - 1))*
     >              dblp_dp**(lna+1) ) 
               
               if (rpow1(i,lna) .lt. regcut*regcut .AND. lna>1) then
                  istartrp(lna)=i+1
!                  exit Debyelna
               end if
               
               if (rpow2(i,lna) .lt. expcut*expcut .AND. lna>1) then
                  istoprp(lna)=i-1
                  exit Debyelna
               end if
                           
            end do Debyelna             !end lna loop
         end do                 !end rmesh loop
         print*,'Debye screening ln,istart,istop'
         do lna = 0, lttop
            print*,lna,istartrp(lna),istoprp(lna)
         end do
c$$$C     Print to files
c$$$         open(89,name='rpow1')
c$$$         open(90,name='rpow2')
c$$$         open(91,name='rpow')
c$$$ 89      format(1P,E12.4,4x,100E15.4)
c$$$C     for l=1
c$$$         do i=1,meshr
c$$$            write(89,89)rmesh(i,1),rpow1(i,0),rpow1(i,1),rpow1(i,2),
c$$$     >           rpow1(i,3),rpow1(i,lttop-1),rpow1(i,lttop)
c$$$            write(90,89)rmesh(i,1),rpow2(i,0),rpow2(i,1),rpow2(i,2),
c$$$     >           rpow2(i,3),rpow2(i,lttop-1),rpow2(i,lttop)
c$$$            write(91,89) rmesh(i,1),(rpow2(i,lna)*rpow1(i,lna),
c$$$     >         lna=0,lttop)
c$$$         end do
c$$$         close(89)
c$$$         close(90)
c$$$         close(91)
      else
C     Calculation of electron-electron Coulomb Potential

      do i=1,meshr
         rpow1(i,0)=1.0
         rpow2(i,0)=1.0/rmesh(i,1)
      end do
      do lna=1,lttop
         istartrp(lna)=istartrp(lna-1)
         istoprp(lna)=istoprp(lna-1)
c$$$         print*,'lna,RPOW1(istoprp(lna-1),lna-1),',
c$$$     >      'RPOW2(istartrp(lna-1),lna-1)',
c$$$     >      lna,RPOW1(istoprp(lna-1),lna-1),
c$$$     >      RPOW2(istartrp(lna-1),lna-1),istartrp(lna-1),istoprp(lna-1)
         do i=istartrp(lna),istoprp(lna)
            rpow1(i,lna)=rpow1(i,lna-1)*rmesh(i,1)
            rpow2(i,lna)=rpow2(i,lna-1)/rmesh(i,1)
         end do
c$$$         do while (rpow1(istartrp(lna),lna).lt.regcut)
c$$$            istartrp(lna)=istartrp(lna)+1
c$$$         end do
c$$$         do while (rpow2(istoprp(lna),lna).lt.expcut*expcut)
c$$$            istoprp(lna)=istoprp(lna)-1
c$$$         end do
C The following stops overflows in single precision
         do while (rpow1(istoprp(lna),lna).gt.1.0/expcut)
            istoprp(lna)=istoprp(lna)-1
         end do
         do while (rpow2(istartrp(lna),lna).gt.1.0/expcut)
            istartrp(lna)=istartrp(lna)+1
         end do
      end do 
      end if
      end
         
C**************************************************************************
C  Debye screened electron-electron potential expansion functions         
      
      real*8 function gfact(n)
      implicit none
      integer::d,n
      real::p
      
      gfact=0
      if (n>0) then
         do d=1,n,2
            p=d
            gfact=gfact+log(p)
         end do
      end if
      end function gfact

C**************************************************************************

      SUBROUTINE POTENT(hlike,NB,LB,NZNUC,NZP,UCENTR)
      INCLUDE 'par.f'
      logical hlike
C  NOTE: the potential has been divided by X
CC
CC  IF NB= 0  THE POTENTIAL IS SET =0
CC
CC  IF NZNUC = 1 AND
CC  NB= 1   GROUND STATE POTENTIAL IS CALCULATED
CC    = 2   AND LB=0  2S EXCITED STATE POTENTIAL
CC    = 2   AND LB=1  2P EXCITED STATE POTENTIAL
CC
CC

      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      common/meshrr/ meshr,rmesh(maxr,3)
      DIMENSION UCENTR(MAXR),vdcore(maxr)
      DIMENSION VCENTR(440)

C     PREPARATION OF POTENTIALS

      do j = 1, meshr
         ucentr(j) = 0.0
      end do
      if (nb.eq.0) return

      small=1e-30
      i=0
      au=1.0
      if (nznuc.eq.-1) then
         if (nb.eq.1) then
            do while(i.lt.meshr.and.au.gt.small)
               i=i+1
               x=rmesh(i,1)
               au=2.0*exp(-2.0*x)*(1.0+1.0/x)
               ucentr(i)=-au
            end do
         else if (lb.eq.1) then
            do while(i.le.meshr.and.au.gt.small)
               i=i+1
               x=rmesh(i,1)
               x2=x*x
               x3=x2*x
               au=2.0*exp(-x)*(1.0/x+0.75+x/4.0+x2/24.0)
               ucentr(i)=-au
            end do
         else if (lb.eq.0) then
            do while(i.le.meshr.and.au.gt.small)
               i=i+1
               x=rmesh(i,1)
               x2=x*x
               x3=x2*x
               au=2.0*exp(-x)*(1.0/x+0.75+x/4.0+x2/8.0)
               ucentr(i)=-au
            end do
         else
            print*,'Invalid input to POTENT'
            stop 'aborted in POTENT'
         endif 
      else if (hlike) then
         jobh = 1
         la = lb
         l = lb
         call makehart(0, la, l, ltmax, psinb(1,nb,lb), 
     >      maxpsinb(nb,lb), enpsinb(nb,lb), jobh, -100000, ucentr)
c$$$         jobh = 2
c$$$         call makehart(0, 0, 0, ltmax, psinb(1,nb,lb), 
c$$$     >      maxpsinb(nb,lb), enpsinb(nb,lb), jobh, -100000, ucentr)
c$$$         i = meshr
c$$$         a = 0.0
c$$$         print*,'Enter the distorting potential fall off factor please'
c$$$         read*,a
c$$$         do while(i.lt.meshr.and.au.gt.small)
c$$$            i=i+1
c$$$            x=rmesh(i,1)
c$$$            au = - (2.0 * vdcore(i) + ucentr(i)) * exp(-a * x * x)
c$$$            ucentr(i) = - au
c$$$         end do
      else
C  Here for non-hydrogenic targets
         jobh = 1
         la = lb
         l = lb
         call makehart(0, la, l, ltmax, psinb(1,nb,lb), 
     >      maxpsinb(nb,lb), enpsinb(nb,lb), jobh, -100000, ucentr)
         ucentr(:) = ucentr(:)*2.0
c$$$         do while(i.lt.meshr.and.au.gt.small)
c$$$            i=i+1
c$$$            x=rmesh(i,1)
c$$$            au=2.0*exp(-2.0*x)*(1.0+1.0/x)
c$$$            ucentr(i)=-au * 2.0
c$$$         end do
      end if 
      i = meshr
      do while (abs(ucentr(i)).eq.0.0.and.i.gt.1)
         i = i - 1
      enddo
c$$$      do j = 1, i
c$$$         write(97,*) j,rmesh(j,1),ucentr(j)
c$$$      enddo 
      print*,'Calculated distorting potential up to R',rmesh(i,1),i
      return
      END

C======================================================================      
C  This subroutine sets up two grids: gridx and gridr. Gridx will be
C  used to solve the radial Schrodinger equation. Gridr, which is
C  typically every second point of gridx, is used to perform integrations
C  from zero to infinity using Simpson's rule. The structure is such
C  that a wave with momentum QMAX will have NPWAVE points per oscillation. This
C  determines HMAX, the largest step. At some point, XDBLE, the intervals,
C  dX, are progressively halved. This is done to accomodate irregular
C  solutions.
C  INPUT:
C   qmax  - the biggest momentum (a.u.) that can be calculated by this mesh.
C   rmax  - the largest "r=x" in the meshes.
C   nmaxx - array declaration parameter of first declaration of GRIDX
C   nmaxr - array declaration parameter of first declaration of GRIDR
C   nmaxj - array declaration parameter of first declaration of JDOUBLE
C  OUTPUT:
C   gridx - X grid
C   nx    - Number of X points
C   gridr - R grid
C   nr    - number of R points
C   jdouble - j points where dx doubles
C   njdouble - Number of doublings + 2
C======================================================================      
      subroutine grids(qmax,ndouble,rmax,gridr,nmaxr,nr,jdouble,
     >   nmaxj,njdouble)
C
C
      dimension gridr(nmaxr,3), jdouble(nmaxj)
C  jdouble stores the J index for GRIDR where DR doubles, with the first
C  point set to 1 and the last to NR. Used in the numerov integrations.

C  Set up GRIDR
C  The number of points per oscillation (half-wavelength): NPWAVE
      npwave = 6
C  The number of doubling of intervals is NDOUBLE
      njdouble = ndouble + 2
      if (njdouble.gt.nmaxj) then
         print*,'Increase NMAXJ in call to GRIDS to at least:',njdouble
         stop 'Increase NMAXJ in call to GRIDS'
      end if 
C  let NPDBL be the number of points with the same dx per interval
      npdbl = 40
      npdbl = 32

C  Make sure NPDBL is even
      npdbl=(npdbl/2) * 2
      if (npdbl.lt.4) then
         print*,'Need to have at least 4 equally spaced points in GRIDR'
         stop 'Need to have at least 4 equally spaced points in GRIDR'
      end if 
      hmax = 3.14/float(npwave)/qmax
      rdble = float(npdbl) * hmax * (2**ndouble - 1) / float(2**ndouble)
      rleft = rmax - rdble
      nleft = int(rleft / hmax) / 2 * 2
      ntot = nleft + npdbl * ndouble
      print*,'Estimated R max:',rdble + nleft * hmax,ntot,
     >   hmax * (float(npdbl) * (2**ndouble - 1) / float(2**ndouble)
     >   + nleft)
c$$$      if (ntot.le.nmaxr) then
c$$$         print*,'Old hmax:',hmax
c$$$         hmax = rmax / (float(npdbl) * (2**ndouble - 1) /
c$$$     >      float(2**ndouble) + nleft)
c$$$         print*,'New hmax:',hmax
c$$$      endif 
C  2**ndouble * hmin = hmax, therefore
      hmin = hmax/float(2**ndouble)
C  The value of the R point from which dR is constant, = HMAX, is RDBLE
C  RDBLE = NPDBL * hmin * (2**NDOUBLE-1)
      rdble = float(npdbl) * hmin * (2**ndouble-1)

      print*,'Grid R parameters:'
      print*,'NDOUBLE:',ndouble
      print*,'HMIN:',hmin
      print*,'HMAX:',hmax
      print*,'NPDBL:',npdbl
      print*,'RDBLE:',rdble
      
      dr = hmin
      jdouble(1) = 1
      do nd = 2, ndouble + 1
         jdouble(nd) = npdbl * (nd - 1)
      end do 

      r=0.0
      j = 0
      do nd = 1, ndouble
         do nj = 1, npdbl
            j = j + 1
            gridr(j,1) = r + float(nj) * dr
            gridr(j,2) = dr
C  Simpson's rule weights
            gridr(j,3) = float(mod(j,2) * 2 + 2) * dr / 3.0
         end do
         gridr(j,3) = dr
         r = gridr(j,1)
         dr = dr * 2.0
      end do

      if (abs((dr - hmax)/hmax).gt.1e-4) then
         print*,'DR and HMAX should be equal'
         stop 'DR and HMAX should be equal'
      end if 
      nj = nint((rmax-r)/hmax) / 2 * 2
c$$$      nj = nj / 32
c$$$      nj = nj * 32
      if (nj+j.gt.nmaxr) then
         nj = nmaxr - j
         print*,'Warning: NR had to be reduced to NMAXR'
      end if 
      do jj = 1, nj
         j = j + 1
         gridr(j,1) = r + float(jj) * dr
         gridr(j,2) = dr
C  Simpson's rule weights
         gridr(j,3) = float(mod(j,2) * 2 + 2) * dr / 3.0
      end do
C  Set the last weight so that the integrals ended exactly at rmax.
C  This is done so that the tail integrals could be done correctly.
      gridr(j,3) = dr / 3.0
      nr = j
      jdouble(ndouble+2) = j
      if (mod(j,2).ne.0)print*,'Warning expected J to be even in grids'

C  Redefine the Simpson's integration weights so that in the integration
C  all intervals were the same
c$$$      j = jdouble(ndouble+1)
c$$$      gridr(j,3) = gridr(j,3) / 3.0 * 4.0
c$$$      rlast = gridr(j,1)
c$$$      j = j - 1
c$$$      w = 1.0
c$$$      do while (j.gt.1)
c$$$         do while(abs((gridr(j,1)+hmax)/rlast-1.0).gt.1e-5.and.j.gt.1)
c$$$            gridr(j,3) = 1e-20
c$$$            j = j - 1
c$$$         enddo
c$$$         w = 3.0 - w
c$$$         rlast = gridr(j,1)
c$$$         gridr(j,3) = 2.0 * w * hmax / 3.0
c$$$         j = j - 1
c$$$      enddo
c$$$      gridr(1,3) = 1e-20
c$$$      
c$$$      j = 1
c$$$      do while (gridr(j,3).lt.1e-10)
c$$$         j = j + 1
c$$$      enddo
c$$$      dr1 = gridr(j,1)
c$$$      j1 = j
c$$$      j = j + 1
c$$$      do while (gridr(j,3).lt.1e-10)
c$$$         j = j + 1
c$$$      enddo
c$$$      j2 = j
c$$$      dr2 = gridr(j,1) - dr1

c$$$      if (gridr(j2,3).lt.gridr(j1,3)) then
c$$$         gridr(j1,3) = (dr1+dr2) / 2.0
c$$$         gridr(j2,3) = gridr(j2,3) / 2.0 + dr2 / 2.0
c$$$      else 
c$$$         gridr(j1,3) = gridr(j1,3) / 2.0 + dr1 / 2.0
c$$$      endif 
      
c$$$      r = 0.0
c$$$      do j = 1, npdbl * (ndouble+2)
c$$$         if (gridr(j,3).gt.1e-10) then
c$$$            print*,j,gridr(j,1) - r,gridr(j,3)/hmax
c$$$            r = gridr(j,1)
c$$$         endif
c$$$      enddo 
      do n = 0,0
         if (n.eq.0) then
            s = 2.0
         else
            s = 3.0
         endif
         w = 1.0
         test = 0.0
         do j = 2**n, nr, 2**n
            w = s - w
            test = test + gridr(j,1) ** 2 * gridr(j,3) * w
         enddo
         if (j-2**n.ne.nr)print*,'Warning last J should be NR',j-2**n,nr
         test = test * 2**n
         print*,2**n,test*3.0/gridr(nr,1)**3
      enddo 
      print*,'Last R and NR:', gridr(nr,1),nr

      return
      end

      subroutine getint(e,sk1,sk1p,gridk,weightk,nt,sum)
      real gridk(nt),weightk(nt)
      r = sk1p / sk1
      sum = 0.0
      do i = 1, nt
         sum = sum + 2.0 * weightk(i) * r / (e - (gridk(i) * r)**2)
      end do 
      rk = sqrt(e)
      tsum = - 2.0 * acoth(sk1p/rk) / rk
      sum = sum + tsum
      end
      
C  Define the k-grid integration intervals and the number of points in
C  each interval
C 1,3,   20,   0.8,    16,     2.5,    6,    6.0,   12,   0.09  
C     npoints, endk, npoints2, endk2, nendk, endp, midnp, width 
      subroutine makeints(mint,sk,nk,rk,np,width,midnp,endk,nendk,endp,
     >   np2,endk2)
      dimension sk(0:10), nk(10)
      nmin(npoints) = min(npoints/2,4)
      
      if (midnp.lt.0) then
         if (rk.gt.endk2) stop 'On shell point in third interval'
         nk(1) = np
         sk(1) = endk
         nk(2) = np2
         sk(2) = endk2         
         mint = 3
         nk(mint) = nendk
         sk(mint) = endp
         return
      endif 

      sk(0) = 0.0
c$$$      if (abs(width).lt.0.01) then
c$$$         print*,'Width must be greater than 0.01'
c$$$         stop 'Width must be greater than 0.01'
c$$$      endif
      if (endk+2.0*width.ge.endk2) then
         print*,'ENDK2 must be >= ENDK + 2 * WIDTH',endk2,endk+2.*width
         stop 'ENDK2 must be >= ENDK + 2 * WIDTH'
      endif 
      if (rk.lt.0.0.or.(width.lt.0.0.and.rk.lt.2.0*abs(width))) then
C  Here for closed channels
c$$$         nk(1) = midnp
c$$$         sk(1) = 2.0 * abs(width)
c$$$         nk(2) = np
         nk(1) = np
         sk(1) = endk
         nk(2) = np2
         sk(2) = endk2         
         mint = 3
      else if (rk.lt.width) then
C  Here if singularity is going to be in the first interval
         sk(1) = rk * 2.0
         nk(1) = midnp
         sk(2) = endk
         nk(2) = np
         sk(3) = endk2
         nk(3) = np2
         mint = 4
      else if (rk. lt. endk + 1.0 * width) then
C  Here if singularity is going to be in the second interval
         sk(1) = rk - width
         sk(2) = rk + width

C  Use the following to make the k grid symmetric about rk**2
c$$$            sk(1) = sqrt(rk**2 - 2.0 * rk * width)
c$$$            sk(2) = sqrt(rk**2 + 2.0 * rk * width)

         nk(2) = midnp
         if (rk .lt. endk - 1.0 * width ) then
            nk(1) = max(nmin(np),int(np * rk / (endk - 1.0 * width)))
            nk(1) = min(nk(1), np - nmin(np))
            npleft = np - nk(1)
            sk(3) = endk
            nk(3) = npleft
            sk(4) = endk2
            nk(4) = np2
            mint = 5
         else
            nk(1) = np
            sk(3) = endk2
            nk(3) = np2
            mint = 4
         endif 
      else if (rk. lt. endk2) then
C  Here the singularity is in the third interval
         sk(1) = endk
         nk(1) = np
         sk(2) = rk - width
         sk(3) = rk + width
         nk(3) = midnp
C  The two factors of 4.0 below put extra points in the first interval
         nk(2) = max(nmin(np2),int(np2 * (sk(2)-sk(1)) * 4.0 /
     >      (4.0 * (sk(2) - sk(1)) + (endk2-endk))))
c$$$         nk(2) = max(nmin(np2),int(np2 * (rk - endk) / (endk2-endk)))
         nk(2) = min(nk(2), np2 - nmin(np2))
         npleft = np2 - nk(2)
         if (rk .lt. endk2 - width) then
            sk(4) = endk2
            nk(4) = npleft
            mint = 5
         else
            nk(2) = np2
            mint = 4
         endif    
      else
C  This is for very large incident energies
         sk(1) = endk
         nk(1) = np
         sk(2) = rk - width
         nk(2) = np2
         sk(3) = rk + width
         nk(3) = midnp
         mint = 4
      endif
      nk(mint) = nendk
      sk(mint) = endp
      return
      end
      subroutine wrtcs(partcs,sigtop,nchpmax,nchipmax,instate,lg,etot,
     >   nsmax,ipar,ovlp,nunit,nznuci,zasymi,projectile,target,canstop,
     >   BornICS,BornPCS,csfile,esum,rcond,nnbtop,ovlpnn)
      include 'par.f'
      dimension partcs(nchan,nchan,0:1), oldpj(nchan,nchan,0:1),
     >   oldp(nchan,nchan,0:1), sigionold(nchan,0:1), sigb(nchan,0:1), 
     >   sigion(nchan,0:1), sigt(nchan,0:1), sigtop(nchan,0:1),
     >   sigtopold(nchan,0:1), sigtopoldj(nchan,0:1),
     >   sigione(nchan,0:1), sigtope(nchan,0:1), sigtJ(nchan,0:1), 
     >   sigionl(nchan,0:lamax,0:1),sigtl(nchan,0:lamax,0:1),
     >   sigbl(nchan,0:lamax,0:1),tbextra(knm),tnbBextra(knm),
     >   ticssum(nchan,0:1),ticsextra(knm),ovlpnn(knm,nnmax),
     >   PsTNBCS(nchan),PsTICS(nchan),PsTCS(nchan),
     >   PsTNBCSl(nchan,0:lamax),PsTICSl(nchan,0:lamax),
     >   PsTCSl(nchan,0:lamax),PsTICSJ(nchan),btics(2,nchan),
     >   PsTNBCSJ(nchan),SCSl(0:lamax,0:1,knm),SCS(0:1,knm),
     >   TNBCSJ(nchan,0:1),TICSJ(nchan,0:1)
      complex ovlpnn
     
      real units(3), BornICS(knm,knm),ovlp(knm),esum(0:1),rcond(0:1),
     >   BornPCS(nchan,nchan),oldBornPCS(nchan,nchan),ext(knm,knm)
      common /corepartarray/ corepart(KNM,nicmax),ecore(nicmax)
      common /corearray/ nicm, ncore(nspmCI)
      real*8 t(ncmax,0:lamax),wts(ncmax),wf(2*ncmax),dstart,dstop,
     >   y(ncmax,0:lamax),corepart,ecore,dt
      integer iwf(2*ncmax), nt(0:lamax), nlast(5,0:lmax), instate(100)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ip,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      common /MPI_info/myid, ntasks, cnode, ench
      logical canstop,exists,assigned
      character*8 chunit(3), chun
      character*3 chan(knm),chs(0:lamax)
      character projectile*(*), target*(*), csfile*(*), file_J*80,
     >   lockfile*80,cnode*3, ench*11
      common /charchan/ chan
      common /chanen/ enchan(knm)
      real sigblast(5,0:lamax),sigbprev(5,0:lamax),
     >   tiecs(nicmax,0:1,nchan),tiecse(nicmax,0:1,nchan)
      save ext
      asym(sing,trip,fac) = (sing - trip / fac) / (sing + trip + 1e-30)
      
      data pi,ry,chunit/3.1415927,13.6058,
     >   ' a0^2   ',' pi a0^2',' cm^2  '/
      data units/1.0,0.3183099,28.002e-18/

      tiecs(:,:,:) = 0.0
      tiecse(:,:,:) = 0.0
      btics(:,:) = 0.0
      ibl = 1
      do while (csfile(ibl:ibl).ne.' ')
         ibl = ibl + 1
      enddo
      lockfile = csfile(1:ibl-1)//'_LOCK'
      n = 1
 10   inquire(file = lockfile, exist=exists)
      if (exists) then
         print*,'Surprisingly lockfile exists:',lg,n,lockfile
         call sleep(lg)
         if (n.le.100) then
            n = n + 1
c$$$            call datetime(3)
         else
            stop 'LOCK file exists after N = 100'
         endif
         go to 10
      endif
      open(40,file = lockfile,err=10)
      write(40,*) lg,lockfile
      close(40,status='keep')
      inquire(file = lockfile, exist=exists)
      if (.not.exists) stop 'LOCK mechanism not working'
         
      call datetime(1)
      unit = units(nunit)
      nwf=2*ncmax
      niwf=2*ncmax

      fac = 2.0
      if (chan(1)(1:1).eq.' ') fac = 3.0
      small = 5e-2
      canstop = .true.
      do ns = 0, 1
         do nchip = 1, nchan
            sigtopold(nchip,ns) = 0.0
            do nchp = 1, nchan
               oldp(nchp,nchip,ns) = 0.0
               oldpj(nchp,nchip,ns) = 0.0
               oldBornPCS(nchp,nchip) = 0.0
            enddo
         enddo
      enddo

C  We redefine NCHPMAX and NCHIPMAX here because on occasion NCHPMAX for
C  the unnatural parity case comes out to be smaller than in the natural
C  parity case
      open(42,file=csfile)
      j=-1
c$$$      read(42,*,end=20) j, nchpmaxt, nchipmaxt
      read(42,"(3i4,7x,a8)",end=20) j, nchpmaxt, nchipmaxt, chun
      if (chun.eq.chunit(1)) then
         nunitold = 1
      elseif (chun.eq.chunit(2)) then
         nunitold = 2
      elseif (chun.eq.chunit(3)) then
         nunitold = 3
      else
         stop 'units do not match in WRTCS'
      endif
      conversion = units(nunit) / units(nunitold)
      if (nunit.ne.nunitold) print*,'WARNING: unit change in WRTCS'
         
      if (nchpmaxt.ne.nchpmax.or.nchipmaxt.ne.nchipmax) then
         print*,'WARNING have nchpmaxt<>nchpmax or nchipmaxt<>nchipmax',
     >      nchpmaxt,nchpmax,nchipmaxt,nchipmax
         print*,'Perhaps this is due to IPAR:',ipar
         nchipmax = nchipmaxt
         nchpmax  = nchpmaxt
         print*,'Set nchipmax and nchpmax to:', nchipmax, nchpmax
      endif
      do incount = 1, nchipmax
         nchip = instate(incount)
         read(42,*) 
         do nchp = 1, nchpmax
            if (BornICS(nchp,nchip).eq.0.0) then
               read(42,'(11x,e13.5,e15.5)') oldBornPCS(nchp,nchip)
     >            ,BornICS(nchp,nchip)
               oldBornPCS(nchp,nchip) = oldBornPCS(nchp,nchip) *
     >            conversion
               BornICS(nchp,nchip) = BornICS(nchp,nchip) * conversion
            else
               read(42,'(11x,e13.5,e15.5)') oldBornPCS(nchp,nchip)
               oldBornPCS(nchp,nchip) = oldBornPCS(nchp,nchip) *
     >            conversion
            endif
         enddo
      enddo 
      do ns = 0, nsmax
         do incount = 1, nchipmax
            nchip = instate(incount)
            sigtopold(nchip,1) = 0.0
            sigtopoldj(nchip,1) = 0.0
            read(42,*) sigtopold(nchip,ns), sigtopoldj(nchip,ns)
            sigtopold(nchip,ns) = sigtopold(nchip,ns) * conversion
            sigtopoldj(nchip,ns) = sigtopoldj(nchip,ns) * conversion
            do nchp = 1, nchpmax
               oldp(nchp,nchip,1) = 0.0
               oldpj(nchp,nchip,1) = 0.0
               read(42,'(3i3,4e15.5)') nsp,nf,ni,oldp(nchp,nchip,ns),
     >            en, overlap, dt ! dt is real*8 because it may be too small
               oldp(nchp,nchip,ns) = oldp(nchp,nchip,ns) * conversion
               oldpj(nchp,nchip,ns) = dt * conversion
            enddo
         enddo
      enddo
 20   close(42,status='delete')
      open(42,file=csfile)
      file_J = csfile(1:ibl-1)//'_J'
      open(43,file=file_J)
      close(43,status='keep')
      call openatend(43,file_J)
      write(42,"(3i4,' Units:',a8,2x,a8,' - ',a6,
     >   99(1p,e12.4,' eV on ',a3))")
     >   max(lg,j), nchpmax, nchipmax,chunit(nunit),
     >   projectile,target,(ry * (etot-enchan(instate(i))),
     >   chan(instate(i)),i=1,nchipmax) 
      write(43,"(3i4,' Units:',a8,2x,a8,' - ',a6,1p,
     >   99(/,e12.5,' eV on ',a3))")
     >   max(lg,j), nchpmax, nchipmax,chunit(nunit),
     >   projectile,target,(ry * (etot-enchan(instate(i))),
     >   chan(instate(i)),i=1,nchipmax) 
      do incount = 1, nchipmax
         nchip = instate(incount)
         tbextra(nchip) = 0.0
         tnbBextra(nchip) = 0.0
         ticsextra(nchip) = 0.0
         write(42,'("  transition   BornPCS        BornICS",
     >      "      last PCS(V)    last PCS(T) canstop")') 
         do nchp = 1, nchpmax
            partcsT = partcs(nchp,nchip,0) + partcs(nchp,nchip,1)
            if (BornICS(nchp,nchip).ne.0.0) then
               extra = max(0.0,BornICS(nchp,nchip) -
     >            BornPCS(nchp,nchip)*unit - oldBornPCS(nchp,nchip))
            else
               extra = 0.0
            endif 
            tbextra(nchip) = tbextra(nchip) + extra
            if (enchan(nchp).lt.0) then
               tnbBextra(nchip) = tnbBextra(nchip) +
     >            ovlp(nchp) * extra
            else 
               ticsextra(nchip) = ticsextra(nchip) + extra
               btics(1,nchip) = btics(1,nchip) + BornICS(nchp,nchip)
               if (enchan(nchp).lt.etot/2.0)
     >            btics(2,nchip) = btics(2,nchip) + BornICS(nchp,nchip)
            endif 
            diff = abs((partcsT - BornPCS(nchp,nchip))/(partcsT+1e-30))
            write(42,'(a5,'' <-'',a3,1p,e13.5,3e15.5,l5)') 
     >         chan(nchp),chan(nchip), BornPCS(nchp,nchip) * unit +
     >         oldBornPCS(nchp,nchip), BornICS(nchp,nchip),
     >         BornPCS(nchp,nchip) * unit, (partcs(nchp,nchip,0) +
     >         partcs(nchp,nchip,1)) * unit, diff.lt.small.or.
     >         abs(BornPCS(nchp,nchip)*unit+oldBornPCS(nchp,nchip)-
     >         BornICS(nchp,nchip))/(BornICS(nchp,nchip)+1e-30).lt.small
         enddo
      enddo 
      do ns = 0, nsmax
         do incount = 1, nchipmax
            nchip = instate(incount)
            ein = max(0.0,etot-enchan(nchip))
            SCS(ns,nchip) = 0.0
            sigt(nchip,ns) = 0.0
            sigtJ(nchip,ns) = 0.0
            sigb(nchip,ns) = 0.0
            sigt(nchip,1) = 0.0
            sigtJ(nchip,1) = 0.0
            sigb(nchip,1) = 0.0
            sigione(nchip,ns) = 0.0
            sigione(nchip,1) = 0.0
            sigion(nchip,ns) = 0.0
            sigion(nchip,1) = 0.0
            sigtope(nchip,ns) = 0.0
            sigtope(nchip,1) = 0.0
            ticssum(nchip,ns) = 0.0
            ticssum(nchip,1) = 0.0
            PsTNBCS(nchip) = 0.0
            PsTICS(nchip) = 0.0
            PsTICSJ(nchip) = 0.0
            TICSJ(nchip,0:1) = 0.0
            PsTNBCSJ(nchip) = 0.0
            TNBCSJ(nchip,ns) = 0.0
            PsTCS(nchip) = 0.0
            do l = 0, lamax
               SCSl(l,ns,nchip) = 0.0
               PsTNBCSl(nchip,l) = 0.0
               PsTICSl(nchip,l) = 0.0
               PsTCSl(nchip,l) = 0.0
               sigionl(nchip,l,ns) = 0.0
               sigtl(nchip,l,ns) = 0.0
               sigbl(nchip,l,ns) = 0.0
               sigionl(nchip,l,1) = 0.0
               sigtl(nchip,l,1) = 0.0
               sigbl(nchip,l,1) = 0.0
C  The following is used to estimate 1/n**3 scaling. NC is 1 for singlet,
C  2 for doublet, 3 for triplet, 4 for SINGLET and 5 for TRIPLET
               do nc = 1, 5
                  sigblast(nc,l) = 0.0
                  sigbprev(nc,l) = 0.0
                  nlast(nc,l) = 0
               enddo 
            enddo
            ltop = 0
            lPstop = 0
            sigbo = 0.0
            do nchp = 1, nchpmax
C  SIGT will be the total cross section for all J, calculated by summing 
C  individual cross sections. SIGTJ is same, but for the J only.
               sig = partcs(nchp,nchip,ns) * unit + oldp(nchp,nchip,ns)
               sigt(nchip,ns) = sigt(nchip,ns) + sig
               sigtJ(nchip,ns) = sigtJ(nchip,ns) +
     >            partcs(nchp,nchip,ns)
               call getchnl(chan(nchp),n,l,nc)
               if (chan(nchp)(1:1).eq.'p') then
                  PsTCS(nchip) = PsTCS(nchip) + sig
                  PsTCSl(nchip,l) = PsTCSl(nchip,l) + sig
                  if (enchan(nchp) .lt.0.0) then
                     PsTNBCS(nchip) = PsTNBCS(nchip) + sig * ovlp(nchp)
                     PsTNBCSJ(nchip) = PsTNBCSJ(nchip) +
     >                  partcs(nchp,nchip,ns) * unit
                     PsTNBCSl(nchip,l) = PsTNBCSl(nchip,l) +
     >                  sig * ovlp(nchp)
                  else
c$$$                     PsTICS(nchip) = PsTCS(nchip) - PsTNBCS(nchip)
c$$$                     PsTICSl(nchip,l) = 
c$$$     >                  PsTCSl(nchip,l) - PsTNBCSl(nchip,l)
                     PsTICS(nchip) = PsTICS(nchip) + sig
                     PsTICSl(nchip,l) = PsTICSl(nchip,l) + sig
                     PsTICSJ(nchip) = PsTICSJ(nchip) +
     >                  partcs(nchp,nchip,ns) * unit
                  endif                      
                  if (l.gt.lPstop) lPstop = l
               endif 
               if (l.gt.ltop) then
                  ltop = l
                  if (ltop.gt.lamax) then
                     print*,'chan,l,lamax:',chan(nchp),l,lamax
                     stop 'l > lamax'
                  endif
               endif 
               sigtl(nchip,l,ns) = sigtl(nchip,l,ns) + sig
               if (BornICS(nchp,nchip).ne.0.0) then
                  extra = max(0.0,BornICS(nchp,nchip) -
     >               BornPCS(nchp,nchip)*unit-oldBornPCS(nchp,nchip))
               else
                  extra = 0.0
               endif 
               SCS(ns,nchip) = SCS(ns,nchip) +
     >            0.5*(enchan(nchp)-enchan(nchip))*(sig+extra) !atomic units for energy
               SCSl(l,ns,nchip) = SCSl(l,ns,nchip) +
     >            0.5*(enchan(nchp)-enchan(nchip))*(sig+extra) !atomic units for energy

               if (enchan(nchp).lt.0.0.and.chan(nchp)(1:1).ne.'p') then
!               do nb = nabot(l), max(nnbtop,nabot(l))
                  tmp = sig !* ovlpnn(nchp,nb) ** 2
c$$$                  print*,nb,ovlpnn(nchp,nb)
                  sigbo = sigbo + tmp
                  sigb(nchip,ns) = sigb(nchip,ns) + tmp
                  sigbl(nchip,l,ns) = sigbl(nchip,l,ns) + tmp
                  TNBCSJ(nchip,ns) = TNBCSJ(nchip,ns) +
     >               partcs(nchp,nchip,ns) * unit
!               enddo
c$$$              print*,'nchp,sigt-sigb,sigt,sigb',nchp,sigtl(nchip,l,ns)-
c$$$     >            sigbl(nchip,l,ns),sigtl(nchip,l,ns),sigbl(nchip,l,ns)
               endif 
               if (enchan(nchp) .gt. 0.0) then
                  TICSJ(nchip,ns) = TICSJ(nchip,ns) +
     >                 partcs(nchp,nchip,ns) * unit
                  sigionl(nchip,l,ns) = sigionl(nchip,l,ns) + sig
                  sigion(nchip,ns) = sigion(nchip,ns) + sig
                  ticssum(nchip,ns) = ticssum(nchip,ns) + sig
                  ticsextra(nchip) = ticsextra(nchip) +
     >                 extrap((partcs(nchp,nchip,0)+
     >                 partcs(nchp,nchip,1))*unit,
     >                 oldpj(nchp,nchip,0) + oldpj(nchp,nchip,1), 0.0)

                  if (nicm.gt.nicmax) stop 'Increase NICMAX'
                  if (BornICS(nchp,nchip).ne.0.0) then
                     extra = max(0.0,BornICS(nchp,nchip) -
     >                  BornPCS(nchp,nchip)*unit-oldBornPCS(nchp,nchip))
                  else
                     extra = 0.0
                  endif 
                  do ic = 1, nicm ! Here for ionisation with excitation
                     contr = sig * corepart(nchp,ic)
                     contre = (sig+extra) * corepart(nchp,ic)
                     if (0.5*enchan(nchp).gt.ecore(ic)-ecore(1)) then ! a.u.
                        tiecs(ic,ns,nchip)=tiecs(ic,ns,nchip) + contr
                        tiecse(ic,ns,nchip)=tiecse(ic,ns,nchip) + contre
                     else  ! if below threshold put contr to nc=1 for now!
                        tiecs(1,ns,nchip)=tiecs(1,ns,nchip) + contr
                        tiecse(1,ns,nchip)=tiecse(1,ns,nchip) + contre
                     endif
                  enddo 
               else
c$$$                  sigbl(nchip,l,ns) = sigbl(nchip,l,ns) + sig *
c$$$     >               ovlp(nchp)
c$$$                  sigb(nchip,ns) = sigb(nchip,ns) + sig *
c$$$     >               ovlp(nchp)
                  sigbprev(nc,l) = sigblast(nc,l)
                  sigblast(nc,l) = sig
                  nlast(nc,l) = n
c$$$                  print*,'Compare:',chan(nchp),
c$$$     >               sigbprev(nc,l)*(n-1)**3/n**3, sig, n
               endif 
            enddo

c$$$            print*,'SIGB, SIGBO:',sigb(nchip,ns), sigbo
c$$$            print*,'TICS: sum, old proj, new proj',ticssum(nchip,ns),
c$$$     >         sigt(nchip,ns) - sigb(nchip,ns), sigt(nchip,ns)-sigbo
C  Redefine the ionization cross sections by SIGT - SIGB
            do l = 0, ltop
c$$$               sigionl(nchip,l,ns) = max(0.0,
c$$$     >            sigtl(nchip,l,ns)-sigbl(nchip,l,ns))
               do nc = 1, 5
                  if (nlast(nc,l).ne.0) then
c$$$                     print*,'Last bound excitation cross sections:',
c$$$     >                  sigbprev(nc,l),sigblast(nc,l),nlast(nc,l)
                  endif
               enddo 
            enddo
            
            write(42,'(1p,2e13.5,0p,f10.5,a4,
     >         '' TCS (op. th.), last cont., ratio to sum'')')
     >         sigtop(nchip,ns) * unit +
     >         sigtopold(nchip,ns), sigtop(nchip,ns) * unit,min(1e2,
     >         (sigtop(nchip,ns) * unit + sigtopold(nchip,ns))/
     >         (sigt(nchip,ns)+1e-30)), chan(nchip)
            write(43,'(1p,e10.4,"eV on ",a3," J =",i3,
     >         "  sum, opt, rat. TCS(",i1,",",i1,"): ",
     >         1p,2e10.3,0p,f7.4)') ry * ein,
     >         chan(nchip),lg, ns,ip,sigtJ(nchip,ns) * unit,
     >         sigtop(nchip,ns) * unit,
     >         sigtJ(nchip,ns)/(sigtop(nchip,ns)+1e-30)

C  The following statement tries to fix the problem of getting non-zero
C  ionization cross sections below threshold, due to the optical
C  theorem not being satisfied to full precision.
            if (abs((sigtop(nchip,ns) * unit + sigtopold(nchip,ns))/
     >         (sigt(nchip,ns)+1e-30) - 1.0).lt.0.01)
     >         sigtopold(nchip,ns)=sigt(nchip,ns)-sigtop(nchip,ns)*unit
            sum = 0.0
            sumold = 0.0
            sume = 0.0
            sumo = 0.0
            do nchp = 1, nchpmax
               sig = partcs(nchp,nchip,ns) * unit + oldp(nchp,nchip,ns)
               call getchnl(chan(nchp),n,l,nc)
               write(42,'(3i3,1p,4e15.5,a4,'' <-'',a3)') 
     >            ns,nchp,nchip, sig,enchan(nchp), ovlp(nchp),
     >            partcs(nchp,nchip,ns) * unit, chan(nchp),chan(nchip)
               if (chan(nchp)(1:1).eq.'p') then
                  nbstart = npbot(l)
               else
                  nbstart = nabot(l)
               endif 
!               do nb = nbstart, max(nbstart,nnbtop)
               if (enchan(nchp) .lt. 1e-10) then
C  SUM will be the total non-break-up cross section for all J
c$$$                  sum = sum + !ovlpnn(nchp,nb) ** 2 *
c$$$     >               sig
                  sum = sum + ovlp(nchp) * sig
C  SUMOLD will be the total non-break-up cross section for all previous J
c$$$                  sumold=sumold+!ovlpnn(nchp,nb)**2*
c$$$     >               oldp(nchp,nchip,ns)
                  sumold = sumold + ovlp(nchp)*oldp(nchp,nchip,ns)
C  SUMO will be the total non-break-up cross section for the current J only
c$$$                  sumo = sumo + !ovlpnn(nchp,nb) ** 2 *
c$$$     >               partcs(nchp,nchip,ns) * unit
                  sumo = sumo + ovlp(nchp) *
     >               partcs(nchp,nchip,ns) * unit
C  SUME will be the total non-break-up cross section for the previous J only
c$$$                  sume = sume + !ovlpnn(nchp,nb)**2*
c$$$     >               oldpj(nchp,nchip,ns)
                  sume = sume + ovlp(nchp) * oldpj(nchp,nchip,ns)
                  if (ovlp(nchp).gt.1.001) print*,
     >               'overlap too big, should be <= 1.0',n,l,ovlp(nchp)
               endif
!               enddo 
            enddo
C  Use the optical theorem to define the total and ionization cross sections
C  as then the code is suitable for both CCC and CCO. For CCC both 
C  forms should give much the same answer.
C  SIGION will be the total ionization cross section for all J
c$$$            sigion(nchip,ns) = max(sigtop(nchip,ns) * unit
c$$$     >         + sigtopold(nchip,ns) - sum,0.0)
c$$$C  SIGIONOLD will be the total ionization cross section for all previous J
c$$$            sigionold(nchip,ns) = max(sigtopold(nchip,ns) - sumold,0.0)
c$$$C  SIGIONE will be the extrapolated total ionization cross section
c$$$            sigione(nchip,ns) = sigionold(nchip,ns) +
c$$$     >         extrap(sigtop(nchip,ns) * unit - sumo,
c$$$     >         sigtopoldj(nchip,ns) - sume,0.0)
C  SIGTOPE will be the extrapolated total cross section
            sigtope(nchip,ns) = sigtopold(nchip,ns) +
     >         extrap(sigtop(nchip,ns) * unit,sigtopoldj(nchip,ns),0.0)
C  For pure Born approximation we get zero for the optical theorem. In this
C  case use the form below. Can't define sigionold in this case, need more work
C  Need to make SIGTOLD from the info above
c$$$            if (sigtop(nchip,ns).eq.0.0) then
c$$$               sigion(nchip,ns) = max(sigt(nchip,ns) - sum,0.0)
c$$$               sigionold(nchip,ns) = 0.0
c$$$            endif 
            write(43,'(1p,e10.4,"eV on ",a3," for partial wave J =",i3,
     >         " TICS(",i1,",",i1,"): ",1p,e11.3)')
     >         ry * max(0.0,ein),chan(nchip),lg,ns,ip,
     >         TICSJ(nchip,ns) !(sigtop(nchip,ns) * unit - sumo) 
            write(43,'(1p,e10.4,"eV on ",a3," for partial wave J =",i3,
     >         "   TNBCS(",i1,",",i1,"): ",1p,e11.3)')
     >         ry * ein,chan(nchip),lg,ns,ip,
     >         (TNBCSJ(nchip,ns)) 
            write(43,'(1p,e10.4,"eV on ",a3," for partial wave J =",i3,
     >         " PsTNBCS(",i1,",",i1,"): ",1p,e11.3)')
     >         ry * ein,chan(nchip),lg,ns,ip,
     >         (PsTNBCSJ(nchip)) 
            write(43,'(1p,e10.4,"eV on ",a3," for partial wave J =",i3,
     >         " PsTICS(",i1,",",i1,"): ",1p,e11.3)')
     >         ry * ein,chan(nchip),lg,ns,ip,
     >         (PsTICSJ(nchip)) 
         enddo ! incount
      enddo ! ns
      
      do ns = 0, nsmax
         write(43,'(1p,e10.4,"eV on ",a3," for J =",i3,
     >      " eigenphase sum and Rcond(",i1,",",i1,"):",1p,2e11.3)')
     >      ry * (etot-enchan(1)),chan(1),lg,ns,ip,
     >      esum(ns),rcond(ns)
         write(42,'(1p,e9.3,"eV on ",a3," for J =",i3,
     >      " eigenphase sum and Rcond(S=",i1,"):",1p,2e11.3)')
     >      ry * (etot-enchan(1)),chan(1),lg,ns,esum(ns),rcond(ns)
      enddo 
      do incount = 1, nchipmax
         nchip = instate(incount)
         ein = max(0.0,(etot-enchan(nchip)))
         sumion = 0.0
         do l = 0, ltop
            chs(l) = ' + '
            if (l.eq.0) chs(l) = ' = '
            do ns = 0, nsmax
               sumion = sumion + sigionl(nchip,l,ns)
            enddo
         enddo 
         write(42,'(20x,"CS        summed l=  ",i3,100i11)')
     >      (l, l = 0, ltop)
         if (projectile.eq.'positron') then
            write(42,'(1p,e9.3,"eV on ",a3," PsTNBCS:    ",
     >         1p,e8.2,19(a3,e8.2))') ry * ein,
     >         chan(nchip),PsTNBCS(nchip),
     >         (chs(l),PsTNBCSl(nchip,l),l=0,lPstop)
            write(42,'(1p,e9.3,"eV on ",a3,"  PsTICS:    ",
     >         1p,e8.2,19(a3,e8.2))') ry * ein,
     >         chan(nchip),PsTICS(nchip),
     >         (chs(l),PsTICSl(nchip,l),l=0,lPstop)
            write(42,'(1p,e9.3,"eV on ",a3,"   PsTCS:    ",
     >         1p,e8.2,19(a3,e8.2))') ry * ein,
     >         chan(nchip),PsTCS(nchip),
     >         (chs(l),PsTCSl(nchip,l),l=0,lPstop)
            write(42,'(79a)') ('-',i=1,79)
         else             
            do ns = 0, nsmax
               write(42,'(1p,e9.3,"eV on ",a3," TNBCS(S=",i1,"): ",
     >            1p,e8.2,19(a3,e8.2))') ry * ein,
     >            chan(nchip),ns,sigb(nchip,ns),
     >            (chs(l),sigbl(nchip,l,ns),l=0,ltop)
               write(42,'(1p,e9.3,"eV on ",a3,"  TICS(S=",i1,"): ",
     >            1p,e8.2,19(a3,e8.2))') ry * ein,
     >            chan(nchip),ns,sigion(nchip,ns),
     >            (chs(l),sigionl(nchip,l,ns),l=0,ltop)
               write(42,'(1p,e9.3,"eV on ",a3,"   TCS(S=",i1,"): ",
     >            1p,e8.2,19(a3,e8.2))') ry * ein,
     >            chan(nchip),ns,sigt(nchip,ns),
     >            (chs(l),sigtl(nchip,l,ns),l=0,ltop)
               write(42,'(1p,e9.3,"eV on ",a3,"   SCS(S=",i1,"):",
     >            1p,e9.2,19(a3,e8.2))') ry * ein,
     >            chan(nchip),ns,max(0.0,SCS(ns,nchip)),
     >            (chs(l),max(0.0,SCSl(l,ns,nchip)),l=0,ltop)              
               write(42,'(79a)') ('-',i=1,79)
            enddo
         endif 
         write(42,'(1p,e9.3,''eV on '',a3,
     >      '' BTICS: all, to E/2'',
     >      1p,4e11.3)') ry * ein,
     >      chan(nchip),btics(1,nchip),btics(2,nchip)
         write(42,'(1p,e9.3,"eV on ",a3," TNBCS:      ",
     >      1p,e8.2,19(a3,e8.2))') ry * ein,
     >      chan(nchip),sigb(nchip,0)+sigb(nchip,1),
     >      (chs(l),sigbl(nchip,l,0)+sigbl(nchip,l,1),l=0,ltop)
         write(42,'(1p,e9.3,"eV on ",a3,"  TICS:      ",
     >      1p,e8.2,19(a3,e8.2))') ry * ein,
     >      chan(nchip),sigion(nchip,0)+sigion(nchip,1),
     >      (chs(l),sigionl(nchip,l,0)+sigionl(nchip,l,1),l=0,ltop)
         write(42,'(1p,e9.3,"eV on ",a3,"   TCS:      ",
     >      1p,e8.2,19(a3,e8.2))') ry * ein,
     >      chan(nchip),sigt(nchip,0)+sigt(nchip,1),
     >      (chs(l),sigtl(nchip,l,0)+sigtl(nchip,l,1),l=0,ltop)
         write(42,'(79a)') ('-',i=1,79)
            
         write(42,'(1p,e9.3,''eV on '',a3,
     >      '' TNBCS,+extra, spin asymmetry:'',
     >      1p,4e11.3)') ry * ein,
     >      chan(nchip),sigb(nchip,0)+sigb(nchip,1),
     >      sigb(nchip,0)+sigb(nchip,1)+tnbBextra(nchip),
     >      asym(sigb(nchip,0),sigb(nchip,1),fac)

         write(42,'(1p,e9.3,''eV on '',a3,
     >      '' TICS: s, t-s, +, a'',
     >      1p,4e11.3)') ry * ein,
     >      chan(nchip),ticssum(nchip,0) + ticssum(nchip,1), 
     >      sigion(nchip,0) + sigion(nchip,1),
     >      sigion(nchip,0) + sigion(nchip,1) + ticsextra(nchip), 
     >      asym(sigion(nchip,0), sigion(nchip,1),fac)
         do ic = 1, nicm
            write(42,'(1p,e9.3,"eV on ",a3," TIECS(",i2,
     >         "), +extrap, and asym:",1p,3e11.3)')
     >         ry * ein,
     >         chan(nchip),ic,tiecs(ic,0,nchip)+tiecs(ic,1,nchip),
     >         tiecse(ic,0,nchip)+tiecse(ic,1,nchip),
     >         asym(tiecse(ic,0,nchip),tiecse(ic,1,nchip),fac)
         enddo 
         
         do ns = 0, nsmax
            write(43,'(i3,1p,e10.3,"eV on ",a3,
     >         "  TICS(",i1,",",i1,"): ",
     >         1p,e8.2,19(a3,e8.2))') lg, ry * ein,
     >         chan(nchip),ns,ip,sigion(nchip,ns),
     >         (chs(l),sigionl(nchip,l,ns),l=0,ltop)
         enddo 
         write(43,'("J=",i3,1p,e10.3,''eV on '',a3,
     >      '' TICS, +extra, spin asym:'',
     >      1p,4e11.3)')lg, ry * ein,
     >      chan(nchip),sigion(nchip,0) + sigion(nchip,1),
     >      sigion(nchip,0) + sigion(nchip,1) + ticsextra(nchip),
     >      asym(sigion(nchip,0), sigion(nchip,1),fac)
         
c$$$         write(42,'(1p,e9.3,''eV on '',a3,
c$$$     >      " TICS: ",1p,e9.3,19(a3,e8.2))') 
c$$$     >      ry * ein,chan(nchip),sumion,
c$$$     >      (chs(l),sigionl(nchip,l,0) + sigionl(nchip,l,1),l=0,ltop)
         
         sigtsum = sigt(nchip,0) + sigt(nchip,1)
         sigtopt = (sigtop(nchip,0) + sigtop(nchip,1)) * unit +
     >      sigtopold(nchip,0) + sigtopold(nchip,1)
         write(42,'(1p,e9.3,''eV on '',a3,
     >      '' TCS (op. th.), + extra, asym:'',
     >      1p,3e11.3)') ry * ein,chan(nchip),
     >      sigtopt, sigtopt + tbextra(nchip),
c$$$     >      sigtopt, sigtope(nchip,0) + sigtope(nchip,1),
     >      asym(sigtope(nchip,0),sigtope(nchip,1),fac)
         write(42,'(79a)') ('-',i=1,79)
         write(42,'(''transition  cross section  extrapolated'',
     >      ''     overlap       spin asym    energy'')')
         write(43,'('' J   trans  cross section  extrap    PCS(V) '',
     >      ''   PCS(T) S=0 PCS(T) S=1  energy     ovlp  ip'')')
         do l = 0, ltop
            nt(l) = 0
         enddo 
         do nchp = 1, nchpmax
            call getchnl(chan(nchp),n,l,nc)
            summedcs = partcs(nchp,nchip,0)*unit + oldp(nchp,nchip,0) +
     >         partcs(nchp,nchip,1) * unit + oldp(nchp,nchip,1)
            Borne = max(0.0,BornICS(nchp,nchip) - 
     >         oldBornPCS(nchp,nchip) - BornPCS(nchp,nchip) * unit)
C  The following is not right if Born extrapolation is used due to
C  the fact that the spin weights are not available here
c$$$            extrapcs0 = extrap(partcs(nchp,nchip,0) * unit,
c$$$     >         oldpj(nchp,nchip,0), Borne) + oldp(nchp,nchip,0)
c$$$            extrapcs1 = extrap(partcs(nchp,nchip,1) * unit,
c$$$     >         oldpj(nchp,nchip,1), Borne) + oldp(nchp,nchip,1)
c$$$            asymcs = asym(extrapcs0, extrapcs1,fac)
            asymcs = asym(partcs(nchp,nchip,0)*unit+oldp(nchp,nchip,0),
     >         partcs(nchp,nchip,1) * unit + oldp(nchp,nchip,1), fac)
            if (chan(nchip).eq.chan(nchp)) then ! elastic scattering
               if (target.ne.'H  I'.or.chan(nchp)(2:2).eq.'1') then
                  const = (partcs(nchp,nchip,0)+partcs(nchp,nchip,1)) *
     >               (2*lg-1)*(2*lg+1)*(2*lg+3)*(2*lg+5)*(2*lg+7)*unit
                  extrapcs = const/8.0/(2*lg+1)/(2*lg+3)/
     >               (2*lg+5)/(2*lg+7)
               else !excited H state elastic scattering
                  const = (partcs(nchp,nchip,0)+partcs(nchp,nchip,1)) *
     >               (2*lg-1)*(2*lg+1)*(2*lg+3)*unit
                  extrapcs = const/4.0/(4.0*(lg+1)**2-1.0)
               endif 
            else if (BornICS(nchp,nchip).ne.0.0) then
               extrapcs = Borne
            else 
c$$$               summedcs = oldp(nchp,nchip,0) + oldp(nchp,nchip,1) !need to sort this out, Igor
               extrapcs = extrap((partcs(nchp,nchip,0)+
     >            partcs(nchp,nchip,1))*unit,
     >            oldpj(nchp,nchip,0) + oldpj(nchp,nchip,1), 0.0)
            endif
! Below doubled up extrapolation when Borne was non zero for H-like targets
c$$$            if (ipar.eq.0) then 
c$$$               ext(nchp,nchip) = extrapcs
c$$$            else
c$$$               extrapcs = extrapcs + ext(nchp,nchip)
c$$$            endif
            extrapcs = summedcs + extrapcs
            partcsT = partcs(nchp,nchip,0) + partcs(nchp,nchip,1)
            diff = abs((partcsT - BornPCS(nchp,nchip))/(partcsT+1e-30))
            canstop = canstop.and.(diff.lt.small.or.
     >         abs(extrapcs-summedcs)/(summedcs+1e-30).lt.small)
            write(42,'(a3,'' <-'',a3,1p,4e15.5,e11.3)') chan(nchp),
     >         chan(nchip),summedcs, extrapcs, ovlp(nchp), asymcs,
     >         enchan(nchp)
c$$$     >         ,(diff.lt.small.or.
c$$$     >         abs(extrapcs-summedcs)/(summedcs+1e-30).lt.small)
            write(43,'(i3,a4,'' <-'',a3,1p,6e11.3,0p,f8.4,i2)') 
     >         lg,chan(nchp),chan(nchip),summedcs,extrapcs,
     >         BornPCS(nchp,nchip)*unit,
     >         partcs(nchp,nchip,0)*unit,partcs(nchp,nchip,1)*unit, 
     >         enchan(nchp),ovlp(nchp),ip
         enddo
      enddo 
      close(43)
      close(42)
      call datetime(2)
      open(40,file=lockfile)
      close(40,status='delete')
      inquire(file = lockfile, exist=exists)
      if (exists) print*,'LOCKFILE not deleted'
      return
      end
         
      
C  The EXTRAP function works in the following way. We assume that we have 
C  a geometric decreasing series, with last two elements being XP and X 
C  (XP > X). The function returns an estimate of the contribution to the 
C  total sum from elements X = X(J) > X(J+1) > X(J+2) ...
C  If, however, the Born approximation is known, then we use the integrated
C  BornICS - BornPCS as the extrapolation instead. Born-based
C  extrapolation works well for dipole transitions.
      function extrap(x,xp,Borne)
      extrap = x
      if (Borne.gt.0.0.and.x.gt.0.0) then
         extrap = max(0.0,Borne)
      else if (0.0.lt.x.and.x.lt.xp*0.9) then
         extrap = x / (1.0 - x / xp + 1e-20)
      endif
      end
      
      function extrapn3(c,na)
      extrapn3 = c/(na+1)**3
      n = na + 2
      do while (extrapn3 + c/n**3 .gt. extrapn3)
         extrapn3 = extrapn3 + c/n**3
         n = n + 1
      enddo
c$$$      print*,'exiting extrapn3 with n:',n
      return
      end
      
      subroutine simpson(dstart,dstop,nt,x,w)
      implicit real*8 (a-h,o-z)
      dimension x(nt),w(nt)

      k = 0
      if (dstart.lt.1e-10) k = 1

      if (mod(nt+k,2).eq.0) stop
     >   'Wrong number of points for the Simpson rule'
      h = (dstop - dstart) / (nt - 1 + k)
      do n = 1, nt
         x(n) = dstart + h * (n - 1 + k)
         w(n) = 2 * (1 + mod(n+1-k,2)) * h / 3.0
      enddo
      if (k.eq.0) w(1) = w(1)/2.0
      w(nt) = w(nt) / 2.0
      sum = 0d0
      do n = 1, nt
         sum = sum + x(n)**2 * w(n)
      enddo
      exact = (dstop**3 - dstart**3) / 3d0
      if (abs(exact/sum - 1d0).gt.1d-10) print*,'problems in SIMPS',
     >   exact/sum
      return
      end
      
      subroutine trapez(dstart,dstop,nt,x,w)
      implicit real*8 (a-h,o-z)
      dimension x(nt),w(nt)

      h = (dstop-dstart) / dfloat(nt - 1)
      do n = 1, nt
         x(n) = dstart + dfloat(n-1) * h
         w(n) = h
      enddo
      w(1) = h / 2d0
      w(nt) = h / 2d0
      if (abs(dstart).lt.1d-10) then
         x(1) = h / 2d0
         w(1) = h / 2d0
         w(2) = 0.75d0 * h
      endif
      sum = 0d0
      do n = 1, nt
         sum = sum + x(n) * w(n)
      enddo
      exact = (dstop**2 - dstart**2) / 2d0
      if (abs(exact/sum - 1d0).gt.1d-10) print*,'problems in TRAPEZ',
     >   exact/sum
      return
      end
      

      subroutine makehart(lg, lwell, l, lttop, fr, imaxfr, en, jobh,
     >   ntop, ucentx)
C
C         Common blocks
C
      include 'par.f'
      common /MESHRR/ meshr, rmesh(maxr, 3)
      Common /POWERS/ rpow1(maxr, 0:ltmax), rpow2(maxr, 0:ltmax),
     >   iminrp(0:ltmax), imaxrp(0:ltmax), cntfug(maxr, 0:lmax)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast, fastor
      COMMON/CONT/ CE(ncmax),CINT(ncmax),NCSTATES,ENERGY
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common /radpot/ ucentr(maxr)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)

      dimension fr(imaxfr), ucentx(maxr), pot0(maxr), temp(maxr),
     >   fun(maxr)

      if (imaxfr.eq.0) stop 'IMAXFR = 0 in MAKEHART'
      do i = 1, imaxfr
         temp(i) = fr(i)
      end do
c$$$      if (lwell.le.lnabmax) then
c$$$         do nq = max(lwell+1,ntop+1), nnmax
c$$$            do i = 1, min(istoppsinb(nq,lwell),imaxfr)
c$$$               temp(i) = temp(i) + psinb(i,nq,lwell)
c$$$            end do
c$$$         end do
c$$$      end if
c$$$ 10   print*,'Enter I please', imaxfr
c$$$      read*,i
c$$$      do nq = 1, ncstates
c$$$         print*, i,cint(nq) * psinc(i,nq,lwell) *
c$$$     >      en * en / ce(nq)
c$$$      end do
c$$$      go to 10
c$$$      
c$$$      do nq = 1, ncstates
c$$$         do i = istartpsinc(nq,lwell), imaxfr
c$$$            temp(i) = temp(i) + cint(nq) * psinc(i,nq,lwell) *
c$$$     >         en * en / ce(nq)
c$$$         end do
c$$$      end do
      if (lwell.ne.0.or.l.ne.0.or.lg.ne.0) print*,
     >   'WARNING: In MAKEHART the nuclear subtraction is outside the',
     >   ' LT loop'
C
C  Start the LAMBDA = LT loop
C
      do i = 1, meshr
         pot0(i) = 0.0
      end do
      do i = 1, imaxfr
         fun(i) = temp(i) * fr(i) * rmesh(i,3)
      end do 
      ltstop = 2 * min(l, lwell)
      if (jobh.eq.2) ltstop = 0
      fastor = fast
      fast = .false.
      do lt = 0, ltstop, 2
         call cleb(2*l, 2*l, 2*lt, 0, 0, 0, c1) 
         call cleb(2*lwell,2*lwell,2*lt, 0, 0, 0, c2) 
         call rac7(2*lwell, 2*l, 2*lwell, 2*l, 2*lg, 2*lt, c3) 
         c = c1 * c2 * c3 *
     >      float((2*lwell + 1) * (2*l + 1))/float(2*lt + 1) *
     >      float((-1)**lg)
         if (abs(c).gt.1e-10) then
            call form(fun,1,imaxfr,rpow1(1,lt),rpow2(1,lt),
     >         iminrp(lt),imaxrp(lt),meshr,temp,i1,i2)
            do i = i1, i2
               pot0(i) = pot0(i) + c * temp(i)
            end do
         end if
      end do
      fast = fastor
      
C  add the nuclear term
      do i = 1, meshr
c$$$         pot0(i) = pot0(i) + vnucl(i) / 2.0
c$$$         pot0(i) = pot0(i) + ucentr(i) / 2.0
         pot0(i) = pot0(i) - rpow2(i,0) 
c$$$         pot0(i) = pot0(i) - rpow2(i,0) - temp(i)
      end do
      j = meshr
      do while (abs(pot0(j)).lt.formcut.and.j.gt.1)
         pot0(j) = 0.0
         j = j - 1
      enddo
      if (j.eq.meshr) print*,'WARNING: the distorting potential is',
     >   ' non zero at last R point:', pot0(j)      
      do j=1,meshr
         ucentx(j)= pot0(j) * 2.0
      end do
      return
      end

      function oscil(ef,psif,maxf,ei,psii,maxi,rmesh,meshr)
      include 'par.f'
      dimension psif(maxr), psii(maxr), rmesh(maxr,3)

      sum = 0.0
      do i = 1, min(maxf,maxi)
         sum = sum + psif(i) * rmesh(i,1) * psii(i) * rmesh(i,3)
      enddo
      oscil = (ef - ei) * sum**2 / 3.0
      return
      end
      
      function oscil_modif(ef,psif,maxf,ei,psii,maxi,rmesh,meshr)
      include 'par.f'
      dimension psif(maxr), psii(maxr), rmesh(maxr,3)
      common /di_el_core_polarization/ gamma, r0_di_el, pol(maxr)

      sum = 0.0
      do i = 1, min(maxf,maxi)
         sum = sum + psif(i) * (rmesh(i,1) - gamma*pol(i))
     >      * psii(i) * rmesh(i,3)
      enddo
      oscil_modif = (ef - ei) * sum**2 / 3.0
      return
      end

      subroutine getstartk(e,rk,startk,stopk,dk,test,xx,ww,nt)
      real*8 xx(nt),ww(nt),err
#ifdef _single
      err = 1d-7
#elif defined _double
      err = 1d-14
#endif
      do while (dk/startk.gt.err)
         startkold = startk
         if (test.gt.0.0) then
            startk = startk - dk
         else
            startk = startk + dk
         endif
c$$$         print*,test,dk,startk,stopk,nmin
         do n = 1, nt
            xx(n) = (xx(n)-startkold)*(stopk-startk)/
     >         (stopk-startkold) + startk
            ww(n) = ww(n) * (stopk-startk)/(stopk-startkold)
            if (xx(n).lt.rk) nmin = n
         enddo 
         testold = test
         test = 2.0/rk*
     >      (atanh(startk/rk)-acoth(stopk/rk))
         do n = 1, nt
            test = 2.0 * ww(n)/(e - xx(n)**2) + test
         enddo 
         if (test*testold.lt.0.0) dk = dk / 2.0
      enddo 
      return
      end
      
      subroutine getstopk(e,rk,startk,stopk,dk,test,xx,ww,nt)
      real*8 xx(nt),ww(nt),err
#ifdef _single
      err = 1d-7
#elif defined _double
      err = 1d-14
#endif
      do while (dk/stopk.gt.err)
         stopkold = stopk
         if (test.gt.0.0) then
            stopk = stopk - dk
         else
            stopk = stopk + dk
         endif
c$$$         print*,test,dk,startk,stopk,nmin
         do n = 1, nt
            xx(n) = (xx(n)-startk)*(stopk-startk)/
     >         (stopkold-startk) + startk
            ww(n) = ww(n) * (stopk-startk)/(stopkold-startk)
            if (xx(n).lt.rk) nmin = n
         enddo 
         testold = test
         test = 2.0/rk*
     >      (atanh(startk/rk)-acoth(stopk/rk))
         do n = 1, nt
            test = 2.0 * ww(n)/(e - xx(n)**2) + test
         enddo 
         if (test*testold.lt.0.0) dk = dk / 2.0
      enddo 
      return
      end
      
      subroutine improvek(e,sk1,nt,gridk,weightk,endf,endk)
      real gridk(nt),weightk(nt)

      endf = sk1
      if (e.lt.0.0) return
      rk = sqrt(e)
      if (sk1.lt.rk) then
         endf = (2.0 * rk - endk)
         if (endf.lt.0.0) stop 'Problem in IMPOVEK; ENDK too big'
         return
      endif 
      k = 1
      do while (gridk(k) .lt. rk)
         k = k + 1
      enddo
      call getint(e,sk1,sk1,gridk,weightk,nt,sum)
      if (sum.gt.0.0) then
         sk1max = sk1
         sk1min = sk1 * (9.0 * rk + gridk(k)) / 10.0 / gridk(k)
         call getint(e,sk1,sk1min,gridk,weightk,nt,sum)
      else 
         sk1min = sk1
         sk1max = sk1 * (9.0 * rk + gridk(k-1)) / 10.0 / gridk(k-1)
         call getint(e,sk1,sk1max,gridk,weightk,nt,sum)
      endif 

      ncount = 0
      do while(abs((sk1max-sk1min)/sk1min).gt.1e-7.and.ncount.lt.90)
         ncount = ncount + 1
         endf = (sk1max + sk1min) / 2.0
         call getint(e,sk1,endf,gridk,weightk,nt,sum)
         if (sum.gt.0.0) then
            sk1max = endf
         else
            sk1min = endf
         endif 
         if (ncount.ge.90) print*,sk1min,sk1max, sum
      enddo 
      end

      subroutine makechil(lg,gk,wk,qcut,zasym,vdcore,npot,ui,ldw,dwpot,
     >   npk,minchilx,chilx,phasel,nchtop,etot,nbnd,abnd,npsbndin,albnd,
     >   sigma,nnbtop,pos,lnch)
      use gf_module
      use chil_module
      include 'par.f'
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   istartrp(0:ltmax),istoprp(0:ltmax),cntfug(maxr,0:lmax)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /double/id,jdouble(22)
      common/smallr/ formcut,regcut,expcut,fast, match
      logical fast, match
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasymfixed,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax),itail
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      dimension npk(nchtop+1),psi(maxr),natomps(nchtop),lnch(nchan,2),
     >   nbnd(0:lmax), abnd(0:lmax), vdcore(maxr,0:lamax), ui(maxr), 
     >   dwpot(maxr,nchan), utemp(maxr), chitemp(maxr), tempi(maxr)
      dimension gk(kmax,nchan), u(maxr), temp(maxr), phasel(kmax,nchan)
c$$$      dimension chil(:,:,:),minchil(:,:) ! gives SIGSEGV 
c$$$      dimension chil(meshr,npk(nchtop+1)-1,2),minchil(npk(nchtop+1)-1,2)
      real*8 dummy(ncmax,ncmax)
      logical torf, positron, pos(nchan),calculated(-lamax:lamax,2)
      complex phasel,sigc,sigma(nchan),wk(kmax*nchan),phase,sig
      common/matchph/rphase(kmax,nchan),trat
      character chan(knm)*3,ch,file*10
      character*20 stringtemp
      common /charchan/ chan
      common /chanen/ enchan(knm)
      common /debye/dbyexists, dblp, rmudeb, zdby
      logical::dbyexists,exists
C     Added by Alex
c$$$      logical analytic
      dimension vnucl(maxr)
      real waveout(maxr,nnmax),erydout(nnmax),Rpl(maxr),Rpg(maxr)
C     End add

      type wave
         real radial(maxr)
         integer min,max
         real en
      end type wave
#ifdef _CRAYFTN
      type(wave) ,allocatable :: psib(:,:), psibd(:,:)
#else
      type(wave) psib(nnmax,-lamax:lamax), psibd(nnmax,0:lamax)
#endif
      allocatable c(:,:)
      data pi/3.14159265358979/
      ch(i)=char(i+ichar('0'))
#ifdef _CRAYFTN
      allocate(psib(nnmax,-lamax:lamax), psibd(nnmax,0:lamax))
#endif
c$$$      inquire(FILE="analytic",EXIST=analytic)
      analytic = nanalytic.le.-2 !.and.e.ge.aenergyswitch
      if (analytic) then
c$$$         open(42,file = "analytic")
c$$$         read(42,*) nbox
c$$$         close(42)
         nkmax = 0
         do nch = 1, nchtop
            if (npk(nch+1)-npk(nch).gt.nkmax) nkmax=npk(nch+1)-npk(nch)
         enddo
         if (allocated(gf)) deallocate(gf)
         allocate(gf(nkmax,nkmax,nchtop))
         if (allocated(c)) deallocate(c)
         allocate(c(maxr,nkmax))    ! Form routine expects maxr
      endif 
      inquire(FILE="ndble",EXIST=torf)
      ndble = 0
      if (torf) then
         open(42,file="ndble")
         read(42,*) ndble
         close(42)
      endif 

      calculated(:,:) = .false. ! used to store psib in analytic call to pseudo
c$$$      allocate (psib(meshr,nnbtop,latop),epsib(nnbtop,latop))
      nz = nint(zasym)
      print*,'Entered MAKECHIL with ZASYM:',zasym
      alpha = 0.0
      utemp(:) = 0.0

      if (npot.lt.0.and.ldw.ge.0) ui(:)=dwpot(:,-npot) ! Redefine UI
C  First define the bound states, if any
c$$$      do l = max(0,lg-latop), ldw
      do l = max(0,lg-latop), min(lg+latop,lamax)         
         if (l.le.ldw) then
            utemp(1:meshr) = 2.0*vdcore(1:meshr,l)-nze * ui(1:meshr)
c$$$            print*,'Added VDCORE to DWPOT:',utemp(1)*rmesh(1,1)
         else
            utemp(:) = 0.0
         endif
         psib(:,l)%en = 0.0
         psib(:,l)%min = 1
         psib(:,l)%max = 0
         if (nbnd(l).gt.0) then
            do n = l+1, nbnd(l)+l !min(nnbtop,nnbm)
               psib(n,l)%radial(1:maxr) = 0.0
               call rnl(nz,n,l,psib(n,l)%radial,psib(n,l)%en,
     >            psib(n,l)%max)
               psib(n,l)%min = 1
               imax = psib(n,l)%max
               print'("n,l,imax,en,ovlp",3i6,1p,1e12.3,0p,f9.3)',
     >            n,l,psib(n,l)%max,psib(n,l)%en,dot_product(
     >               psib(n,l)%radial(1:imax)*rmesh(1:imax,3),
     >               psib(n,l)%radial(1:imax))
            enddo
            if (l.gt.ldw.and.npsbndin.lt.0) then
               psibd(:,l) = psib(:,l) ! use exact eigenstates
            else 
               nps = abs(npsbndin) !-l
               alpha = abnd(l)
               if (nbnd(l).gt.nps) then
                  print*,'CAUTION:l,nbnd(l),nps:',
     >               l,nbnd(l),nps
               endif 
c$$$               if (nz.eq.0) then
c$$$                  alpha = abnd(l)
c$$$               else 
c$$$                  alpha = max(abnd(l),alpha*1.1)
c$$$               endif
               torf = .false.
!               if (alpha.lt.10.0*(zasym+1)) then !.or.zasym.eq.-1.0) then
               niter = 20
               scale = 3.0
               if (alpha.gt.0.0) then
                  print*,"Projectile diagonalisation: N, L, al",nps,l,
     >                 alpha
                  call makeps(zasym, torf,alpha,l,expcut,nps,ps2,
     >               psen2,minps2, maxps2,rmesh,utemp,meshr,meshr,dummy)
                  if (nz.gt.0) then
                     do ntimes = 1, niter
                        npose = 1
                        do while (psen2(npose).lt.0.0)
                           npose = npose + 1
                        enddo
! Below ensures that the first +ve-energy E1 point integrates from zero to 2*E1 by making E2=3E1. Interestingly, the factor 3 also works for the two -ve-energy points, see below print.
                        test=(psen2(npose)*scale-psen2(npose+1))/
     >                       psen2(npose+1)
                        alpha = alpha * (1.0 - test /10.0)
                        print'("nt,npose,al,test,4psen:",2i3,6f10.6)',
     >                       ntimes,npose+l,alpha,test,
     >                       (psen2(n),n=npose-2,npose+1) 
                        if (abs(test).lt.1e-3) exit
                        call makeps(zasym,torf,alpha,l,expcut,nps,
     >                       ps2,psen2, minps2, maxps2, rmesh, utemp,
     >                       meshr,meshr,dummy)
                     enddo
                  endif
               else
                  hmax = rmesh(meshr,2)
                  ra = abs(abnd(l))
                  nbmax = abs(npsbndin) ! - l There is a subtraction of L inside
                  print*,'Projectile Box-based states; Z, N and R:',
     >               -nze*zasym-1.0,nbmax,ra !note that z+1 is used in pseudo
                  call pseudo(jdouble,id,hmax,-nze*zasym-1.0,l,ra,nbmax,
     >               maxr,utemp,psen2,ps2,jmax)
                  minps2(:) = 1
                  maxps2(:) = jmax
                  if (nz.gt.0) then
                     do ntimes = 1, niter
                        npose = 1
                        do while (psen2(npose).lt.0.0)
                           npose = npose + 1
                        enddo
c$$$                   test=(psen2(npose)*0.5+psen2(max(npose-1,1))*0.5)/
c$$$     >                       (psen2(npose) - psen2(max(npose-1,1)))
                        test=(psen2(npose)*scale-psen2(npose+1))/
     >                       psen2(npose+1)
                        ra = ra * (1.0 + test/10.0)
                        print'("nt,npose,Ra,test,4psen:",2i3,6f12.6)',
     >                       ntimes,npose+l,ra,test,
     >                       (psen2(n),n=npose-2,npose+1)
                        if (abs(test).lt.1e-3) exit
                        call pseudo(jdouble,id,hmax,-nze*zasym-1.0,l,ra,
     >                       nbmax,maxr,utemp,psen2,ps2,jmax)
                        minps2(:) = 1
                        maxps2(:) = jmax
                     enddo
                  endif
               endif
               n = l + 1
               psibd(n,l)%en = psen2(n-l)
               do while (n.le.nps+l) ! psibd(n,l)%en.lt.0.0.and.
                  proj = 0.0
                  psibd(n,l)%max = maxps2(n-l)
                  psibd(n,l)%radial(1:maxr) = 0.0
                  psibd(n,l)%radial(1:maxps2(n-l)) =
     >               ps2(1:maxps2(n-l),n-l)
                  imax = psibd(n,l)%max
                  do np = l+1, nbnd(l)+l !min(nnbtop,nnbm)
                     proj = proj + dot_product(
     >                  psibd(n,l)%radial(1:imax)*rmesh(1:imax,3),
     >                  psib(np,l)%radial(1:imax))**2
                  enddo
                  print'("n,l,imax,end,en,proj",3i6,1p,3e12.3)',
     >               n,l,psibd(n,l)%max,psibd(n,l)%en,psib(n,l)%en,proj
                  n = n + 1
                  psibd(n,l)%en = psen2(n-l)
               enddo
            endif 
c$$$            open(42,file='bndk'//ch(l))
c$$$            write(42,'("# radius        dpot  e(Ry):",1p,100e14.6)')
c$$$     >         (psibd(np,l)%en,np=l+1,n-1)
c$$$            do i = 1, psibd(n-1,l)%max
c$$$               write(42,'(1p,100e14.6)') rmesh(i,1),utemp(i),
c$$$     >            (psibd(np,l)%radial(i),np=l+1,n-1)
c$$$            enddo
c$$$            close(42)

c$$$            nend = min(nnmax,npsbndin)
c$$$            psibd(n+1:nend,l)%en = psen2(n+1-l:nend-l)
c$$$            print*,   'n,l,end        ',n,l,psibd(n,l)%en
         endif 
         call update(6)
      enddo
      trat = 1.0
      if (itail.lt.0) then
         trat = 5.0 !gk(1,1)/qcut*rmesh(meshr,1)
         print'("tail integral Rmax, Rtail, ratio:",3f6.1)',
     >      rmesh(meshr,1),rmesh(meshr,1)/trat,trat
      endif 

      do nch = 1, nchtop
C  NCH is the channel index
C  NT is the state index
C  LG is the total partial wave (J)
C  PSI is the state wave function
C  MAXPSI is the last point of the wave function
C  EA, LA and NA are the state's energy (Ry), angular momentum and
C    principle quantum number.
C  l is the projectile angular momentum (L)     
         call getchinfo (nch, nt, lg, psi, maxpsi, ea, la, na, l)
         pos(nch) = positron(na,la,npos)
         lnch(nch,1) = la
         lnch(nch,2) = l
         if (pos(nch)) then
            zeff = 0.0
            npcalc = 2
         else
            zeff = nze * zasym
            npcalc = 1
         endif
         if (l.le.ldw) then
c$$$            call makehart(lg,la,l,ltmax,psi,maxpsi,ea,1,-1,ui)
C  Distorted-waves
            dwpot(:,nch) = ui(:)
            if (ui(1).eq.0.0) stop 'DWPOT = 0.0'
c$$$            if (npot.ge.0) then
c$$$C  Define a single distorting potential (defined in POTENT) for all channels
c$$$               do i = 1, meshr
c$$$                  dwpot(i,nch) = ui(i)
c$$$               enddo
c$$$            else
c$$$C  Define a channel dependent distorting potential, defined in MAKEV31D or VDME
c$$$c$$$               alpha = log(expcut) / rmesh(meshr,1) ** (-npot)
c$$$c$$$               alpha = real(npot)/10.0 + 0.1
c$$$c$$$               m = - npot
c$$$c$$$               print*,'ALPHA:',alpha
c$$$               do i = 1, meshr
c$$$c$$$                  write(98,*) rmesh(i,1), dwpot(i,-npot), ui(i)
c$$$c$$$                  dwpot(i,nch) = dwpot(i,nch)*exp(alpha*rmesh(i,1)**m)
c$$$c$$$                  dwpot(i,nch) = dwpot(i,nch)*exp(alpha*rmesh(i,1))
c$$$                  dwpot(i,nch) = dwpot(i,-npot)
c$$$                  ui(i) = dwpot(i,-npot)
c$$$c$$$                  write(10*na+la,*) rmesh(i,1), - ui(i)
c$$$               enddo
c$$$c$$$               close(10*na+la)
c$$$               if (ui(1).eq.0.0) stop 'DWPOT = 0.0'
c$$$            endif
            ll = min(l,lamax)
            do i = 1, meshr
               u(i) = 2.0*vdcore(i,ll) - nze *
     >            (- zasym * 2.0 / rmesh(i,1) +ui(i))
               utemp(i) = - nze * (2.0*vdcore(i,ll)+ui(i))
            enddo
         else
C  Plane/Coulomb-waves
            do i = 1, meshr
c$$$               u(i) = - nze * (- zasym * 2.0 / rmesh(i,1))
               u(i) = zeff * 2.0 / rmesh(i,1)
               utemp(i) = 0.0
               dwpot(i,nch) = 0.0
            enddo
         endif 

C  LSET may either be L or LG, but must be the same as in routine KGRID
         lset = l
         nqm = npk(nch+1) - npk(nch)
         if (nqm.eq.1) then
C  Here when calculating the pure Born term for subtraction in cross, or 
C  a UBA run.
            nbndm = 0
         else
            nbndm = nbnd(lset)
         endif

         if (nbndm.gt.0.and.u(1).ne.0.0) then
C  We first find the bound states of the U potential.
c$$$            if (l.le.ldw) then
               npsbnd = nps
c$$$               call makeps(zasym,.false.,abnd(l), l, expcut,npsbnd, ps2,
c$$$     >            psen2, minps2, maxps2, rmesh, utemp,meshr,meshr,dummy)
               psen2(:)= 0.0
               n = 1
               do n = 1, nbndm  !while (psibd(n+l,l)%en.lt.0.0)
                  psen2(n) = psibd(n+l,l)%en
                  minps2(n) = 1
                  maxps2(n) = 0
                  ps2(1:maxr,n) = 0.0
                  if (psibd(n+l,l)%en.lt.0.0) then
                     maxps2(n) = psibd(n+l,l)%max
                     ps2(1:meshr,n) = psibd(n+l,l)%radial(1:meshr)
                  endif 
c$$$                  n = n + 1
               enddo
c$$$               psen2(n) = psibd(n+l,l)%en
               nneg = 1
               do while (psen2(nneg).lt.0.0.and.nneg.lt.npsbnd)
                  nneg = nneg + 1
               enddo
               nneg = nneg - 1
               if (nneg.gt.nbndm) then
                  print*,
     >               'Have generated more bound states than requested',
     >               ' NNEG, NBNDM, L:',NNEG, NBNDM, L
                  if (zasym.eq.0.0)
     >            stop'Have generated more bound states than requested'
               endif 
c$$$            else 
c$$$               nneg = nbndm
c$$$               do n = l+1, nbndm + l
c$$$                  call rnl(nz,n,l,temp,en,jstop)
c$$$c$$$                  sum = 0.0
c$$$c$$$                  do i = minps2(n-l),min(maxps2(n-l),jstop)
c$$$c$$$                     sum = sum + temp(i) * ps2(i,n-l) * rmesh(i,3)
c$$$c$$$                  enddo 
c$$$c$$$                  print*,'n, energies and overlap:',n,psen2(n-l),en,sum
c$$$                  minps2(n-l) = 1
c$$$                  maxps2(n-l) = jstop
c$$$                  psen2(n-l) = en
c$$$                  do i = minps2(n-l),maxps2(n-l)
c$$$                     ps2(i,n-l) = temp(i)
c$$$                  enddo
c$$$c$$$                  print*,n,l,en,psib(n,l)%en,
c$$$c$$$     >               dot_product(psib(n,l)%radial(1:psib(n,l)%max),
c$$$c$$$     >               temp(1:jstop)*rmesh(1:jstop,3))
c$$$               enddo
c$$$            endif 

            alpha = abnd(l)
            if (zasym.gt.0.0.and..false.) then
               do while (psen2(nneg+1).gt.0.0)
                  alpha = alpha - alpha/100.0
                  call makeps(zasym,.false.,alpha,l,1e-10,npsbnd,ps2,
     >               psen2, minps2, maxps2, rmesh, utemp, meshr, meshr,
     >               dummy)
               enddo
               alpha = alpha + alpha/100.0
               print*, 'ABND(l) changed to',alpha,l
               call makeps(zasym,.false.,alpha,l,1e-10,npsbnd,ps2,
     >            psen2, minps2, maxps2, rmesh, utemp, meshr, meshr,
     >            dummy)
            endif 
            print '(i2,2x,a3,i3,1p,100e12.5)', nneg, chan(nt), l,
     >         (psen2(i), i=1, nbndm)
         else
            do n = 1, nbndm
C  Any positive number will do
               psen2(n) = 1.0
            enddo 
         endif 

         testb = 0.0
         testc = 0.0
         do n = 1, nbndm
            k = npk(nch+1) - npk(nch) - nbndm + n
            kp = npk(nch+1) - nbndm - 1 + n
            wk(kp) = cmplx(2.0/(etot - ea - psen2(n)))
c$$$            if (etot - ea.lt.0.0) then
c$$$               print*,'Setting Green function to zero for bound',
c$$$     >            ' states in closed channels'
c$$$               wk(kp) = (1e,0.0)
c$$$            endif 
            e1 = etot-ea-psen2(n)
            e2 = etot-ea-psen2(n+1)
            if (n.lt.nbndm.and.e1*e2.lt.0.0) 
     >         print'("CAUTION: WK changes sign across bound states",
     >         i3,2f8.4,a5,i2)',n,e1,e2,chan(nt),l
C     The choice of phase makes no difference
            phasel(k,nch) = (1.0,0.0)
c$$$            phasel(k,nch) = (0.0,1.0)**(-l)
            if (psen2(n).lt.0.0) then
c$$$            if (psen2(n).lt.0.0.and.etot-ea.gt.0.0) then
               gk(k,nch) = -sqrt(-psen2(n))
               if (abs(etot - ea - psen2(n)).lt.0.1) print*,
     >            'Warning for bound state have etot - ea - psen2(n)',
     >            etot - ea - psen2(n), n
               minchil(kp,1) = 1
C The normalization sqrt(pi/2.0) is used as later all CHIL states are
C multiplied by sqrt(2.0/pi).
               do i = 1, maxps2(n)
                  chil(i,kp,1) = ps2(i,n) * rmesh(i,3) * sqrt(pi/2.0)
               enddo
               do i = maxps2(n) + 1, meshr
                  chil(i,kp,1) = 0.0
               enddo
c$$$               if (nbnd(lg).ne.0) then
                  call getprod(psi,maxpsi,1,ps2(1,n),l,rmesh,meshr,nt,
     >               tmp)
                  testb = testb + tmp
c$$$                  write(10*na+la,*) psen2(n),tmp
c$$$               endif 
            else
               gk(k,nch) =  sqrt(psen2(n))
               minchil(kp,1) = min(meshr + 1,maxr)
               do i = 1, meshr
                  chil(i,kp,1) = 0.0
               enddo
            endif 
         end do

C  Having defined the bound states we now define the distorted waves
         summax = 0.0
         ebmax = 0.0
C$OMP parallel do default(private)
C$OMP& shared(npk,nbndm,nch,zeff,gk,meshr,rmesh,u,cntfug,l,ldw,jdouble,
C$OMP& id,regcut,expcut,pi,etot,wk,ea,nze,zasym,vdcore,ui,nt,phasel,ll,
C$OMP& minchil,chil,rphase,nqm,testc,sigma,summax,itail,trat,qcut)
         do k = 1, npk(nch+1) - npk(nch) - nbndm
            kp = npk(nch) + k - 1
            if (gk(k,nch).eq.0.0) then
               eta = zeff
            else 
               eta = zeff / gk(k,nch)
            endif 
            en = gk(k,nch) * abs(gk(k,nch))
            jstart = min(meshr + 1,maxr)
            jstop = meshr
            if (en.ge.0.0) then
               call regular(l,en,eta,u,cntfug(1,l),ldw,
     >            rmesh,meshr,jdouble,id,regcut,expcut,
     >            temp,jstart,jstop,phasel(k,nch),sigc)
c$$$  if (dbyexists) sigc = (1.0,0.0)
c$$$               if (nbnd(lg).ne.0.and.k.gt.1) then
c$$$               if (k.lt.1) then     !Igor
c$$$               if (k.gt.1) then
c$$$                  call getprod(psi,maxpsi,jstart,temp,l,rmesh,
c$$$     >               meshr,nt,tmp)
c$$$                 tmp = tmp * wk(kp) * 2.0 / pi *
c$$$     >               (etot - ea - en) / 2.0
c$$$                  testc = testc + tmp
c$$$               endif 
c$$$               if (lg.gt.ldw) phasel(k,nch) = (1.0,0.0)

C  The following generates the contribution to the T matrix due
C  to the distorting potential.
               if (k.eq.1.and.l.le.ldw) then
                  do i = 1, meshr
c$$$                     utemp(i) = - nze * (- zasym * 2.0 / rmesh(i,1))
                     utemp(i) = zeff * 2.0 / rmesh(i,1)
                  enddo
                  call regular(l,en,eta,utemp,cntfug(1,l),ldw,
     >               rmesh,meshr,jdouble,id,regcut,expcut,
     >               chitemp,j1,j2,phase,sig)
                  do i = 1, meshr
                     utemp(i) = vdcore(i,ll) - nze * ui(i) / 2.0
                  enddo
c$$$                  do i = 1, meshr
c$$$                     utemp(i) = vdcore(i,ll) - nze *
c$$$     >                  (ui(i)-zasym*2.0/rmesh(i,1)) / 2.0
c$$$                  enddo

                  tmp = 0.0
                  do i = jstart, meshr
                     tmp = tmp + temp(i) * chitemp(i) * utemp(i) *
     >                  rmesh(i,3)
c$$$                     write(66,'(5e12.4)') rmesh(i,1),tmp,temp(i),
c$$$     >                  chitemp(i),utemp(i)
                  enddo
                  test = -aimag(phasel(k,nch)) * sqrt(en) / 2.0 / tmp
c$$$                  print*,'SIG, SIGC:',sig,sigc
c$$$                  print*,'phase,phasel:',phase,phasel(k,nch)
c$$$                  print'("NCH, l, En, test and T dist:",
c$$$     >               i4,i3,f7.3,f7.4,1p,2e13.4)', nch, l, en, test,
c$$$     >               tmp*2.0/pi/en*phasel(k,nch)*sig**2
                  if (abs(test-1.0).gt.0.1) 
     >               print*,'Caution distorted waves are inaccurate',
     >               tmp,-aimag(phasel(k,nch)) * sqrt(en) / 2.0
               endif 
            else
               do i = 1, meshr
                  temp(i) = 0.0
               enddo
               phasel(k,nch) = (1.0,0.0)
               sigc = (1.0,0.0)
            endif
            if (k.eq.1) sigma(nch) = sigc
c$$$            if (jstart.lt.-jstop) then   !Igor
c$$$            if (jstart.lt.jstop) then
c$$$               cc = cos(l * pi / 2.0)
c$$$               cs = sin(l * pi / 2.0)
c$$$               rp = real(phasel(k,nch))
c$$$               ap = aimag(phasel(k,nch))
c$$$               rphase(k,nch) = atan2(ap*cc - rp*cs, rp*cc - ap*cs)
c$$$C for large enough r*gk we should have temp(r) = sin(gk*r + rphase)
c$$$               i = meshr
c$$$               x = sin(gk(k,nch) * rmesh(i,1)+rphase(k,nch))
c$$$               y = sin(gk(k,nch) * rmesh(i,1)-l*pi/2.0)
c$$$               z = cos(gk(k,nch) * rmesh(i,1)+l*pi/2.0)
c$$$
c$$$c$$$               if (k.le.9) then
c$$$c$$$                  write(file,'(i1,"_",1p,e8.2)') l,en
c$$$c$$$                  open(50,file=file)
c$$$c$$$                  write(50,*) '# nch, gk, phase',
c$$$c$$$     >               nch,gk(k,nch),rphase(k,nch)
c$$$c$$$                  do i = 1, meshr
c$$$c$$$                     write(50,*) rmesh(i,1), temp(i),
c$$$c$$$     >                  sin(gk(k,nch) * rmesh(i,1)+rphase(k,nch))
c$$$c$$$                  enddo
c$$$c$$$               endif 
c$$$
c$$$               if (abs((x - rp * y - ap * z) / x).gt.1e-3.and.
c$$$     >            abs(x).gt.0.1) then
c$$$                  print*,'Matching process failed in MAKECHIL',
c$$$     >               x,rp * y + ap * z,l,gk(k,nch) 
c$$$                  i = meshr
c$$$                  do while(abs(sin(gk(k,nch)*rmesh(i,1)+rphase(k,nch))/
c$$$     >               (temp(i)+1e-30) - 1.0).lt.0.1.and.i.gt.meshr-10)
c$$$                     print*,i,sin(gk(k,nch) * rmesh(i,1)+rphase(k,nch)),
c$$$     >                  temp(i)
c$$$                     i = i - 1
c$$$                  enddo
c$$$                  print*,nch,k,i
c$$$C below we have temp(r) = rp * sin(k*r-l*pi/2) + ap * cos(k*l+l*pi/2)
c$$$                  i = meshr
c$$$                  do while (abs((rp*sin(gk(k,nch)*rmesh(i,1)-l*pi/2.0)+
c$$$     >               ap*cos(gk(k,nch)*rmesh(i,1)+l*pi/2.0)) /
c$$$     >               (temp(i)+1e-30) - 1.0).lt.0.1.and.i.gt.1)
c$$$                     i = i - 1
c$$$                  enddo
c$$$                  print*,nch,k,i
c$$$c$$$                  do i = 1, meshr
c$$$c$$$                     write(50,*) rmesh(i,1),
c$$$c$$$     >                  sin(gk(k,nch) * rmesh(i,1)+rphase(k,nch))
c$$$c$$$                     write(51,*) rmesh(i,1), temp(i)
c$$$c$$$                     write(52,*) rmesh(i,1),rp*sin(gk(k,nch)*rmesh(i,1)-
c$$$c$$$     >               l*pi/2.0) +ap*cos(gk(k,nch)*rmesh(i,1)+l*pi/2.0)
c$$$c$$$                  enddo
c$$$c$$$                  stop
c$$$               endif
c$$$            else
c$$$               rphase(k,nch) = 0.0
c$$$            endif
            rkmax=0.0
            do nb = 1, nbndm
               kpp = nqm - nbndm + nb
               kppp = npk(nch+1) - nbndm - 1 + nb
               if (gk(kpp,nch).lt.0.0) then
                  sum = 0.0
                  do i = jstart, meshr
                     sum = sum + chil(i,kppp,1) * temp(i)
                  enddo
                  sum = abs(sum)
                  if (sum.gt.summax) then
                     summax = sum
                     nbmax = nb
                     ebmax = - gk(kpp,nch) * gk(kpp,nch)
                     rkmax = gk(k,nch)
                  endif 
               endif 
               if (nb.eq.nbndm.and.k.eq.nqm-nbndm) print*,
     >            'L, Max overlap, bound state e, k:',
     >            l,summax,ebmax,rkmax
            enddo 
            minchil(kp,1) = jstart
C  The T(kf,ki) matrix has been divided by KF and KI so we no longer 
C  divide by sqrt(EN). Store the distorted waves multiplied by the integration
C  weights.
            alpha = log(sqrt(expcut)) / rmesh(meshr,1)
            do i = 1, meshr
               chil(i,kp,1) = 0.0
            end do
c$$$            print*,'l,en,kp,minchil(kp,1),temp(minchil(kp,1)):',
c$$$     >         l,en,kp,minchil(kp,1),temp(minchil(kp,1))
            do i = minchil(kp,1), meshr
               chil(i,kp,1) = temp(i) * rmesh(i,3)
c$$$     >            * exp(alpha * rmesh(i,1))
            end do
            if (itail.lt.0.and.en.ge.0.0) then
               entail = en * (rmesh(meshr,1)/trat)**2
               etatail = eta !zeff / sqrt(entail)
               do i = 1, meshr
                  chil(i,kp,2) = 0.0
               end do
               minchil(kp,2) = meshr
               if (entail.lt.1.1*qcut**2) then !qcut**2 too limited
                  do i = 1, meshr
!                     utemp(i) = - nze * (- zasym * 2.0 / rmesh(i,1))
                     utemp(i) = zeff*2.0*rmesh(meshr,1)/trat/rmesh(i,1)
                  enddo
c$$$                  call regular(l,entail,etatail,utemp,cntfug(1,l),ldw,
c$$$     >               rmesh,meshr,jdouble,id,regcut,expcut,
c$$$     >               temp,minchil(kp,2),jstop,phase,sigc)
c$$$                  if (abs(imag(phasel(k,nch))).gt.1e-3) then
                  call makegreen(l.le.ldw,l,entail,etatail,utemp,
     >               cntfug(1,l),rmesh,meshr,jdouble,id,regcut,expcut,
     >               temp,tempi,minchil(kp,2),jstop,phase)
c$$$                  endif
                  do while(rmesh(minchil(kp,2),1).lt.trat.and.
     >               minchil(kp,2).lt.meshr) !tail integrates from TRAT
                     minchil(kp,2) = minchil(kp,2) + 1
                  enddo 
                  do i = minchil(kp,2), meshr
                     chil(i,kp,2) = (real(phasel(k,nch))*temp(i) +
     >                  imag(phasel(k,nch))*tempi(i)) * rmesh(i,3)
                  enddo
                  if (en.eq.0.0) chil(minchil(kp,2):meshr,kp,2) =
     >               chil(minchil(kp,2):meshr,kp,2)*rmesh(meshr,1)/trat
c$$$                  if (k.eq.1)
c$$$     >               print'("K,ENTAIL,TRAT,MINCHIL:",1p,3e10.2,i7,a4)',
c$$$     >               gk(k,nch),entail,trat,minchil(kp,2),chan(nt)
               else
c$$$                  print*,'CHIL(:,kp,2)=0 for kp>qcut',sqrt(entail),qcut
               endif 
            end if !tails defined
         end do ! End of the loop over K
C$OMP end parallel do
         if (itail.lt.0) print*,'CHIL(:,kp,2)=0 for kp > qcut',qcut

         if (nqm.gt.1.and.nbnd(lg).ne.0) print 
     >      '("test integral for state:",a4,f10.6," =",f10.6," +",
     >      f10.6)',chan(nt),testb+testc,testb,testc
c$$$         close(10*na+la)

C Alex for box basis.
         if (gk(1,nch).GE.0.0) then
            eproj = gk(1,nch) * gk(1,nch)
         else 
            eproj = - gk(1,nch) * gk(1,nch)
         endif
         if (eproj.eq.0.0) then
            eta = zeff
         else 
            eta = zeff / sqrt(abs(eproj))
         endif
         nbmax=npk(nch+1)-npk(nch)-1
         if (analytic .AND. nbmax .GT. 1) then
            if (eproj.ge.aenergyswitch) then !aenergyswitch is in modules.f
               hmax=rmesh(meshr,2)
c$$$               zas=-nze*zasym-1.0 ! check this
               zas=-zeff-1.0 ! check this
               zas = -1.0 !potential is vnucl; hence no (z+1)/r in PSEUDO
               ra=rmesh(meshr,1)
C               vnucl(:)=0.0
               erydout(:)=0.0
               waveout(:,:)=0
               if (gk(1,nch).ne.0.0) then
                  do i = 1,maxr
c$$$                  vnucl(i) = (vdcore(i,l)-rpow2(i,0)*(zasym))*2.0
                     vnucl(i) = (vdcore(i,l)+zeff*rpow2(i,0))*2.0
                  enddo
!                  vnucl(1:meshr) = 0.0
               else !numerov routine wants below for calculation of J and Y
                  do i = 1,maxr
                     vnucl(i) = (vdcore(i,l)-abs(zeff)*rpow2(i,0))*2.0
                  enddo
               endif 
               if (calculated(lg-l,npcalc)) then
                  do k=2,npk(nch+1)-npk(nch)
                     waveout(:,k) = psib(k,lg-l)%radial(:)
                     erydout(k) = psib(k,lg-l)%en
                     maxps2(k) = psib(k,lg-l)%max
                     minps2(k) = psib(k,lg-l)%min
                  enddo 
               else
                  calculated(lg-l,npcalc) = .true. 
                  if (.true..or.nbmax.gt.npsbndin) then
                     zas=-zeff-1.0 ! check this
                     vnucl(:)=0.0
                     print*, "L, nbmax, npk(nch+1)-npk(nch), zas:",l,
     >                  nbmax,npk(nch+1)-npk(nch),zas
                     call pseudo(jdouble,id,hmax,zas,l,ra,nbmax+l,
     >                  maxr,vnucl,erydout(2),waveout(1,2),lastpt)
                     maxps2(:)=lastpt
                     minps2(:)=1
                     kstep = 0
c$$$                  print*,'l,erydout(1,...,3):',l,(erydout(i),i=1,3)
c$$$                  print*,'nbmax-ndble+2,nbmax+1:',nbmax-ndble+2,nbmax+1
                     do k = nbmax-ndble+2,nbmax+1
                     print*,'k,k+kstep,k+kstep+1:',k,k+kstep,k+kstep+1,
     >                     erydout(k+kstep),erydout(k+kstep+1)
                        waveout(1:lastpt,k) = (waveout(1:lastpt,k+kstep)
     >                     + waveout(1:lastpt,k+1+kstep))/sqrt(2d0)
                        erydout(k)=(erydout(k+kstep)+erydout(k+kstep+1))
     >                     / 2.0
                        kstep = kstep + 1
                     enddo 
                  else
                     temp(:) = 0d0
c$$$  nbmax = max(nbmax, npsbndin-la)
                     print*,'Diagonalising with L,N,a,zas,n:',
     >                  l,npsbndin,abnd(l),zas,nbmax
                     call makeps(zasym,.false.,abnd(l),l,expcut,
     >                  npsbndin,waveout(1,2),erydout(2),minps2(2),
     >                  maxps2(2),rmesh, temp, meshr, meshr, dummy)
                  endif 
                  do k = 2, npk(nch+1)-npk(nch)
                     psib(k,lg-l)%radial(:) = waveout(:,k)
                     psib(k,lg-l)%en = erydout(k)
                     psib(k,lg-l)%min = 1 !this needs to change for large L
                     psib(k,lg-l)%max= maxps2(k)
                  enddo 
               endif 
               wk(npk(nch)) = cmplx(0.0,imag(wk(npk(nch)))) !needed for PHOTO(0d0,0d0)
               do i = 1, meshr
c$$$                  utemp(i) = - nze * (- zasym * 2.0/rmesh(i,1))
                  utemp(i) = zeff * 2.0/rmesh(i,1)
               enddo
               do k=2,npk(nch+1)-npk(nch)
                  lastpt = maxps2(k)
                  if (erydout(k).LT.0) then
                     gk(k,nch) = -sqrt(abs(erydout(k)))
                  else
                     gk(k,nch) = sqrt(erydout(k))
                  endif
                  temp(:)=0d0
                  tmp = 1.0
                  if (itail.lt.0.or.nze.eq.1.and.lptop.ge.0) then ! for positrons or analytic tail integrals need projections
                     call regular(l,erydout(k),eta,vnucl,
     >                  cntfug(1,l),ldw,rmesh,meshr,jdouble,id,
     >                  regcut,expcut,temp,j1,j2,phase,sig)
                     imatch = lastpt
c$$$                  do while( abs(waveout(imatch-1,k)).GT.
c$$$     >               abs(waveout(imatch,k)).and.lastpt-imatch.lt.20)
c$$$                     print*,imatch, temp(imatch),waveout(imatch,k)
c$$$                     imatch = imatch - 1
c$$$                  enddo
c$$$                  psilast = waveout(imatch,k)
c$$$                  wk(k+npk(nch)-1) = (psilast / temp(imatch))**2
c$$$c$$$                  wk(k+npk(nch)-1)=1.0/(ra/2.0-(0.25)
c$$$c$$$     >               *sin(2.0*gk(k,nch)*ra)/gk(k,nch))
c$$$                  print*, 'analytic wk:',real(wk(k+npk(nch)-1))
                     temp(1:lastpt)=temp(1:lastpt)*rmesh(1:lastpt,3)
                     tmp=dot_product(waveout(1:lastpt,k),temp(1:lastpt))
                     if (tmp.eq.0.0) then
                        print*,'tmp is 0.0;j1,j2,lastpt,k,waveout,temp:'
     >                   ,j1,j2,lastpt,k,waveout(lastpt,k),temp(lastpt)
                        print*,'erydout(1,...,3):',(erydout(i),i=1,3)
                        stop 'tmp cannot be zero here in makechil'
                     endif
                     entail = erydout(k) * (rmesh(meshr,1)/trat)**2
                     kp = k + npk(nch)-1
                     if (itail.lt.0) then
                        do i = 1, meshr
                           chil(i,kp,2) = 0.0
                        end do
                        minchil(kp,2) = meshr
                     endif 
                     if (entail.lt.qcut**2*2.0.and.itail.lt.0) then
                        call regular(l,entail,eta,utemp,cntfug(1,l),
     >                     ldw,rmesh,meshr,jdouble,id,regcut,expcut,
     >                     temp,minchil(kp,2),j2,phase,sig)
                        do while(rmesh(minchil(kp,2),1).lt.trat) !tail integrates from TRAT
                           minchil(kp,2) = minchil(kp,2) + 1
                        enddo 
                        do i = minchil(kp,2), meshr
                           chil(i,kp,2) = temp(i) * rmesh(i,3)
                        enddo
                        print*,'Defined tail CHIL, MINCHIL for en:',
     >                     entail,minchil(kp,2)
                     endif 
                  else 
                     tmp = 1.0  ! electron scattering
                  endif
                  wk(k+npk(nch)-1) = (1.0,0.0)/tmp**2
c$$$                  print*, 'numeric  wk:',real(wk(k+npk(nch)-1))
                  do j=1,lastpt
c$$$                        chil(j,k+npk(nch)-1,1)=temp(j)*rmesh(j,3)
                     chil(j,k+npk(nch)-1,1)=waveout(j,k)*rmesh(j,3)
     >                  * tmp
!     >                  /sqrt(real(wk(k+npk(nch)-1)))
!     >                  *sqrt(pi/2d0)
                  enddo
                  chil(lastpt+1:meshr,k+npk(nch)-1,1) = 0d0
                  minchil(k+npk(nch)-1,1) = 1 !this needs to change for large L
!                  print*,npk(nch+1)-npk(nch),nps
                  if (npk(nch+1)-npk(nch).le.nps+1) then
                     print'("nch,k,gk,sqrt(abs(enb)):", 2i4,1p,2e12.3)', 
     >                  nch,k,gk(k,nch),sqrt(abs(psibd(k-1+l,l)%en))*
     >                  psibd(k-1+l,l)%en/abs(psibd(k-1+l,l)%en)
                     gk(k,nch) = sqrt(abs(psibd(k-1+l,l)%en))*
     >                  psibd(k-1+l,l)%en/abs(psibd(k-1+l,l)%en)
                     do j = 1, psibd(k-1+l,l)%max
                        chil(j,k+npk(nch)-1,1)=psibd(k-1+l,l)%radial(j)*
     >                     rmesh(j,3) * tmp
                     enddo 
                     chil(psibd(k-1+l,l)%max+1:meshr,k+npk(nch)-1,1)=0d0
                  else
                     if (k.eq.2.or.k.eq.npk(nch+1)-npk(nch)) then
                        print'("nch, k, gk, wk, 2/ra:", 2i4,1p,3e12.3)', 
     >                     nch,k,gk(k,nch),real(wk(k+npk(nch)-1)),2.0/ra
                     endif
                  endif
               enddo
               Rpl(:) = 0.0
               Rpg(:) = 0.0

               GF(:,:,nch)=0.0
               const = sqrt(abs(eproj)) !abs(gk(1,nch))
C  For zero energy GF goes as rho**(l+1)*rho**(-l) leaving a net k 
               if (const.eq.0.0) const = 1.0 !removing / k when it is zero
               vnucl(1:meshr)=-abs(zeff)*rpow2(1:meshr,0)*2.0
               call makegreen(.true.,l,eproj,eta,vnucl,cntfug(1,l),
     >            rmesh,meshr,jdouble,id,regcut,expcut,Rpl,Rpg,
     >            istart,istop,phase)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP& SCHEDULE(dynamic)
C$OMP& PRIVATE(i1,i2)
               do kpp=2,npk(nch+1)-npk(nch)
                  C(:,kpp) = 0.0
                  call form(chil(1,npk(nch)-1+kpp,1),
     >               minchil(kpp+npk(nch)-1,1),meshr,Rpl,Rpg,
     >               istart,istop,meshr,C(1,kpp),i1,i2)
               enddo
C$OMP END PARALLEL DO

c$$$            if (istop.lt.meshr) then
c$$$               do kp=2,npk(nch+1)-npk(nch)
c$$$                  temp(1:istop)=chil(1:istop,kp,1)/rmesh(1:istop,3)
c$$$                  call processform(istop,temp,ichiextreme)
c$$$               enddo
c$$$            endif 

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP& SCHEDULE(dynamic)
               do kpp=2,npk(nch+1)-npk(nch)
                  do kp=2,npk(nch+1)-npk(nch)
                     GF(kp,kpp,nch) =
     >                  - dot_product(chil(1:meshr,npk(nch)-1+kp,1),
     >                  C(1:meshr,kpp)) * pi / const * !*(2.0/pi)
     >                  real(wk(kp+npk(nch)-1))*real(wk(kpp+npk(nch)-1))
                  enddo
               enddo
C$OMP END PARALLEL DO
            else !negative energies
               do kp=2,npk(nch+1)-npk(nch)
                  GF(kp,kp,nch) = real(wk(kp+npk(nch)-1))
c$$$     >               *2.0/(eproj-gk(kp,nch)*abs(gk(kp,nch)))
c$$$                  GF(kp,kp,nch) = pi/(eproj-gk(kp,nch)*abs(gk(kp,nch)))
               enddo 
            endif 

            write(stringtemp,'("GFm",1P,SP,E10.3,"_",SS,I1)') eproj,l
            inquire(file='writegreen',exist=exists)
            if (exists) then
               open(42,file=stringtemp)
               do kp = 2, npk(nch+1)-npk(nch)
                  do kpp = kp,kp !2, npk(nch+1)-npk(nch)
                     write(42,'(1p,3e13.2)') 
     >                  gk(kpp,nch),GF(kp,kpp,nch),
     >                  pi/(eproj-gk(kp,nch)*abs(gk(kp,nch)))
c$$$     >                  gk(kp,nch),gk(kpp,nch),GF(kp,kpp,nch)
                  enddo
c$$$                  write(42,*)
               enddo
c$$$               write(42,'("#   r    ",1p,500e24.3)') (gk(kp,nch),
c$$$     >            kp=2,npk(nch+1)-npk(nch),10)
c$$$               do i = 1, istop
c$$$                  write(42,'(1p,500e12.3)') rmesh(i,1),(c(i,kp),
c$$$     >               chil(i,kp,1)/rmesh(i,3),kp=2,npk(nch+1)-npk(nch),10)
c$$$               enddo 
               close(42)
               print*,'Written:',stringtemp
            endif 
         endif                  !analytic
! As a check, plot 'chil.out' u 1:n, 'chiltail.out' u 1:n
! chiltail should start after chil ends
      inquire(file='waves',exist=exists)
      if (exists) then
         open(42,FILE='chil'//ch(lg)//'_'//ch(nch))
         write(42,'("#",11x,1000F6.3)')
     >      (phasel(k,nch),k=1,npk(nch+1)-npk(nch))
         write(42,'("#",1p,11x,1000E12.4)')(gk(k,nch)*abs(gk(k,nch)),
     >      k=1,npk(nch+1)-npk(nch))
         do i=1,meshr,1
            write(42,'(1p,1000E12.4)') rmesh(i,1),
     >         (chil(i,k,1)/rmesh(i,3),k=npk(nch),npk(nch+1)-1)
         enddo
         close(42)
         if (itail.ne.0) then
            open(42,FILE='chiltail'//ch(lg)//'_'//ch(nch))
            write(42,'("#",11x,1000F6.3)')
     >         (phasel(k,nch),k=1,npk(nch+1)-npk(nch))
            write(42,'("#",1p,11x,1000E12.4)')
     >      ((rmesh(meshr,1)/trat*gk(k,nch))**2,k=1,npk(nch+1)-npk(nch))
            do i=1,meshr,1
               write(42,'(1p,1000E12.4)')rmesh(i,1)*rmesh(meshr,1)/trat,
     >            (chil(i,k,2)/rmesh(i,3),k=npk(nch),npk(nch+1)-1)
            enddo
            close(42)
         endif 
         open(42,file='dwpot')
         do i = 1, meshr
            write(42,*) i,rmesh(i,1),ui(i)
         enddo
         close(42)
         print*,'Written chil for NCH, ZASYM:',nch,zasym
      endif 
      end do ! End of the loop over NCH
      if (allocated(c)) deallocate(c)
      
c$$$      do nchi = 1, nchtop
c$$$         natomps(nchi) = 0
c$$$         do nchf = nchi, nchtop
c$$$c$$$            if (pos(nchf).neqv.pos(nchi)) natomps(nchi)=natomps(nchi)!+1
c$$$c$$$     >           +(2*lnch(nchi,1))*(2*lnch(nchi,2))
c$$$c$$$     >           +2**lnch(nchi,1)
c$$$            if (pos(nchf).neqv.pos(nchi)) then
c$$$               if (pos(nchi)) then
c$$$                  natomps(nchi)=natomps(nchi) !+ 1
c$$$     >                 +1.8**lnch(nchf,1)*1.8**lnch(nchi,1)
c$$$               else
c$$$                  natomps(nchi)=natomps(nchi) !+ 1
c$$$     >                 +1.8**lnch(nchf,1)*1.8**lnch(nchi,1)
c$$$               endif
c$$$            endif
c$$$         enddo
c$$$         if (pos(nchi)) then
c$$$            natomps(nchi) = natomps(nchi)*1.8**lnch(nchi,1)
c$$$         else
c$$$            natomps(nchi) = natomps(nchi)*1.8**lnch(nchi,1)
c$$$         endif
c$$$      enddo 
#ifdef _CRAYFTN
      deallocate(psib,psibd)
#endif
      return
      end


      subroutine processform(istop,chi,ichiextreme)
      dimension chi(istop)

      i = istop
      dchi1 = chi(i)-chi(i-1)
      i = i - 1
      dchi2 = chi(i)-chi(i-1)
      do while (dchi1*dchi2.gt.0.0.and.i.gt.10)
         i = i - 1
         dchi1 = dchi2
         dchi2 = chi(i)-chi(i-1)
      enddo
      ichiextreme = i
      return
      end
      

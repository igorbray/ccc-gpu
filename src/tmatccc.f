C  File: solvet.f
C  The following routine solves the equation
C  T_fi(X_f,X_i) = V_fi(X_f,X_i) * phi(X_f) * phi(X_i) 
C   + Sum_l{ Sum_n w_nl*V_fn(X_f,X_n)*phi(X_f)*conjg(phi(X_n))*T_ni(X_n,X_i)
C   - i * pi * X_l * V_fl(X_f,X_l) * phi(X_f)*conjg(phi(X_l)) * T_fl(X_f,X_l)}.
C  Here sum over the l for the first term is the sum over all channels,
C  whereas the second term is summed over only the open channels.
C  We do this using real arithmetic by letting
C  Sum_i T_fi(X,X_i) * conjg(phi(X_f)) * conjg(phi(X_i)) *
C                 (I(i,i') + i * pi X_i * K_ii'(X_i,X_i')) = K_fi'(X,X_i'),
C  and then solving K(f,i) = V(f,i) + w(n) * V(f,n) * K(n,i),
C  where the sum over n is implied. We want to do this in a
C  symmetric manner so as w(f) = w(i) = 0 we solve for K(n,i) with n <> f
C  by writing (I(f',n)/w(n) - V(f',n)) (w(n) * K(n,i)) = V(f',i),
C  where f' ranges over all n <> f. This is solved for B(n,i) = w(n) * K(n,i),
C  and the required K(f,i) is then formed by K(f,i) = V(f,i) + V(f,n) * B(n,i).
C  This works fine, but gives enormous determinants, so instead we
C  solve (I(f',n) * sig(w(n)) - V'(f',n)) (sig(w(n)) * K'(n,i)) = V'(f',i),
C  where sig(w(n)) = w(n)/|w(n)|, V'(f',n) = V(f',n) * sqrt|w(f')*w(n)|,
C  K'(n,i) = K(n,i) * sqrt|w(n)| and V'(f',i) = V(f',i) * sqrt|w(f')|.
C  The required K(f,i) is recovered by the same equation as before.
      subroutine gettmat(lg,ns,nchtop,iex,npk,vmat,wk,gk,phasel,
     >   sigma,bb,tdist,von,etot,kon,ton,err,vmatop,nchopt,nsmax,ve2ed,
     >   ve2ee,dphasee2e,ephasee2e,nchtope2e,nchmaxe2e,ve2e,te2e,nds,
     >   second,uba,ovlpn,theta,ichi,csfile,projectile,slowery,det,
     >   rcond,vmatp,packed,phaseq,scalapack,chil)
      use gf_module
      use photo_module
      include 'par.f'
      parameter (nmax=kmax*nchan)
      integer npk(nchtop+1), iwork(nmax*5),ifail(nmax),valuesin(8),
     >   valuesout(8)
      real vmat(ichi,ichi+1),vmatp(*),bb(ichi,nchtop,2), !vmatp(ichi*(ichi+1)/2),
     >   ve2ed(nchmaxe2e,ichi),ve2ee(nchmaxe2e,ichi),
     >   tmp2(kmax,kmax),tmp3(kmax),tmp4(kmax,kmax),tmp5(kmax,kmax)
      complex vmatop(kmax,kmax,0:nchanop,nchanop),wk(kmax*nchan)
      complex cdet(2),dphasee2e(nchane2e),ephasee2e(nchane2e),
     >   te2e(nchane2e),ve2e(nchane2e), ke2e(nchane2e,nchan), cvec(nmax)
      real vec(nmax), offwk(nmax),
     >   ionshell(nchan+1),eigv(nmax), k2nd(nchan,nchan), ovlpn(knm)
c$$$      real v(ichi,nchtop,3),work((64+3)*ichi)
c$$$      pointer (ptrw,work)
c$$$      pointer (ptrv,v)
      real, allocatable :: work(:), v(:,:,:)

c$$$      pointer(ptrcv,cv)
c$$$      pointer(ptrcv2,cv2)
c$$$      complex cv(nchan,ichi), cv2(ichi,nchtop)
      complex, dimension(:,:), allocatable :: cv, cv2

      complex*16 coulphase
      real*8 deta,vgf
      complex ton(nchan,nchan),von(nchan,nchan), c,
     >   kon(nchan,nchan), ckern(nchan,nchan), cwork(nchan)
      complex phasel(kmax,nchan), tdist(nchan), sigma(nchan),
     >   phaseq(knm)
      dimension det(2), gk(kmax,nchan), temp(maxr),err(nchan,nchan)
      character ud(0:1), date*8,time*10,zone*5
!    fix for cray compiler
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      logical sprint,lprint,uba(nchan),second,packed,exists,scalapack,
     >   converged,photon,lastns
      character chan(knm)*3, ch, csfile*(*), projectile*(*), pfile*80
      common /uniquev/unv(kmax*nchan)
      common /charchan/ chan
      common /chanen/ enchan(knm),enchandiff(knm)
      COMMON /gstspin/   iSpin, meta
      common/meshrr/ meshr,rmesh(maxr,3)
c$$$      COMMON/dipole/  dr (kmax,nchan),dv (kmax,nchan),dx (kmax,nchan)

      ch(i)=char(i+ichar('0'))
      data pi/3.1415927/
      data ud/'L','U'/
c$$$      data sprint,lprint/.false.,.false./
c$$$      data sprint,lprint/.true.,.true./

C     Added by Alex
c$$$      real, allocatable :: kernel(:,:)
      logical :: box
C     End added by Alex

      call date_and_time(date,time,zone,valuesin)

      inquire (file='sprint',exist=sprint)
      inquire (file='lprint',exist=lprint)
      if (lprint) then
         open(42,file='lprint')
         read(42,*,err=10) nistart,nistop,nfstart,nfstop
         close(42)
         goto 20
 10      continue
         stop 'specfy start/stop values in lprint'
 20      continue
      endif
      r = float((-1)**ns * min(1,iex))
      nd = npk(nchtop+1) - 1
      if (nd .ne. ichi) STOP 'ND .NE. ICHI in GETTMAT'
      nds = npk(nchtop+1) - 1 ! - nchtop
      nk = 0
      do nch = 1, nchtop
c$$$         if (gk(1,nch).gt.0.0) then
c$$$            nk = 1
c$$$            do while
c$$$     >         (real(wk(npk(nch)+nk))*real(wk(npk(nch)+nk+1)).gt.0.0)
c$$$               nk = nk + 1
c$$$            enddo
c$$$            ionshell(nch) = npk(nch) + nk + 1
c$$$c$$$            print*,nch,gk(nk+1,nch),gk(1,nch),gk(nk+2,nch)
c$$$c$$$            print*,nch,real(wk(npk(nch)+nk)), real(wk(npk(nch)+nk+1))
c$$$         endif
!         do k = npk(nch)+1, npk(nch+1)-1
         do k = npk(nch), npk(nch+1)-1
            nk = nk + 1
            offwk(nk) = real(wk(k))
         enddo 
      enddo
      converged = .false.
      if (nk.ne.nds) stop 'NK .ne. NDS in tmatccc.f'
      inquire(file='iterate',exist=exists)
      if (exists) then
         open(42,file='iterate')
         read(42,*) maxit,rerr
         close(42)
         print '("ITERATE entered at:",a10," with MAXIT:",i3,
     >      " and RERR:",1p,e8.1)',time,maxit,rerr
         call date_and_time(date,time,zone,valuesin)
         call iterate(gk,kmax,ns,npk,nchtop,vmat,wk,nds,k2nd,nchan,
     >      maxit,rerr,converged)
         call date_and_time(date,time,zone,valuesout)
         print '("ITERATE exited at:  ",a10,", diff (secs):",i5)',
     >      time, idiff(valuesin,valuesout)
         call update(6)
      endif 
      
      if (ns.eq.0) then
c$$$         do ki = 1, npk(nchtop+1) - 1
c$$$            do kf = ki, npk(nchtop+1) - 1
c$$$               vns = vmat(kf,ki) + r * vmat(ki,kf+1) 
c$$$               vnt = vmat(kf,ki) - r * vmat(ki,kf+1)
c$$$               vmat(kf,ki) = vns 
c$$$               vmat(ki,kf+1) = vnt
c$$$            end do
c$$$         end do
      endif 
         
c$$$      do m = 1, min(nchan,nd)
c$$$         call makebasis(m,vmat,r,nqm,nchtop,wk,v)
c$$$         do nchf = 1, nchtop
c$$$            call makevec(vmat,1,nchf,nqm,nchtop,wk,r,vec)
c$$$            do nchi = 1, 1
c$$$               sum = sdot(nd,vec,1,v(1,nchi,2),1)
c$$$               sum = sum + onshellv(vmat,nqm,nchtop,nchf,nchi,r)
c$$$               err(nchf,nchi) = sum / gk(1,nchi) / gk(1,nchf)
c$$$            enddo
c$$$         enddo
c$$$         print'(i5,1p,5(e10.3,2x))',m,(err(nchf,1),nchf=1,min(nchtop,5))
c$$$      enddo 
 

c$$$
c$$$      nchi = 1
c$$$      kqi = 1
c$$$      do nchf=1,nchtop
c$$$         nqmf = npk(nchf+1)-npk(nchf) 
c$$$         do nqf=1,nqmf
c$$$            kqf = nqf + npk(nchf) - 1
c$$$            
c$$$            write(*,'("**",3I5,2E17.6)') nchf, kqf,kqi, 
c$$$     >           gk(nqf,nchf),vmat(kqf,kqi) 
c$$$            
c$$$         enddo
c$$$      enddo

C  Form the driving vector and the first Born
      mv = 4 * nd * nchtop * 3
c$$$      call memalloc(ptrv,mv)
c$$$      if (ptrv.eq.0) stop 'Not enough memory for V'
      allocate(v(nd,nchtop,4))
      v(:,:,:) = 0.0
c$$$      k2nd(:,:) = 0.0

C Set TON to VON before adding the 2nd term later
      do nchi = 1, nchtop
         kf = 1
         if (scalapack) then
            do nchf = 1, nchtop
               ton(nchf,nchi) = bb(npk(nchf),nchi,2)
            enddo
         else 
            call makev(vmat,kf,nchi,nchtop,ichi,npk,wk,ns,
     >         v(1,nchi,1),vmatp,packed)
            do n = 1, nd
               v(n,nchi,2:4) = v(n,nchi,1)
            enddo
            do nchf = 1, nchtop
               ton(nchf,nchi)=onshellk(vmat,nchtop,nchf,nchi,ns,npk,
     >            ichi,vmatp,packed)
c$$$               print*,nchf,nchi,ton(nchf,nchi)
c$$$               if (exists) print*,'nchf,nchi,ns,ton,kon:',
c$$$     >            ton(nchf,nchi),k2nd(nchf,nchi)
            enddo
         endif 
      enddo

C  The following is not coded correctly as ND should be NDS, with appropriate
C  mods.
c$$$      do nchi = 1, nchtop
c$$$         do ni = 1, nd
c$$$            s = (real(wk(ni))+1e-30) / abs(real(wk(ni))+1e-30)
c$$$            vec(ni) = v(ni,nchi,1) * s
c$$$         enddo         
c$$$         do nchn = 1, nchtop
c$$$            vec(npk(nchn)) = 0.0
c$$$         enddo 
c$$$         do nchf = nchi, nchtop
c$$$            tmp = 0.0
c$$$            if (nd.gt.0) tmp = sdot(nd,vec,1,v(1,nchf,1),1)
c$$$            k2nd(nchf,nchi) = real(ton(nchf,nchi)) + tmp
c$$$            k2nd(nchi,nchf) = k2nd(nchf,nchi)
c$$$         enddo 
c$$$      enddo

      epsil = 1e-4
C  Define the on-shell V, used for printing out only
      do nchi = 1, nchtop
         do nchf = 1, nchtop
            von(nchf,nchi) = ton(nchf,nchi) * phasel(1,nchf) *
     >         phasel(1,nchi) / (gk(1,nchi)+1e-10) / (gk(1,nchf)+1e-10)
         end do 
      end do 

      if (converged) then
         do nchi = 1, nchtop
            do nchf = nchi, nchtop
               ton(nchf,nchi) = k2nd(nchf,nchi)
               ton(nchi,nchf) = k2nd(nchf,nchi)
            enddo 
         enddo
      else
         rcond = 0.0
         det(:) = 0.0
         if (scalapack) then
            print*,'Using BB from main'
            do nchn = 1, nchtop
               do kn = 1, nd
                  v(kn,nchn,1)=bb(kn,nchn,1) ! solution
                  v(kn,nchn,2:4)=bb(kn,nchn,2) ! driving term
               enddo
            enddo
         else  
c$$$         do nchi = 1, nchtop
c$$$            do nchf = 1, nchtop
c$$$               do kf = npk(nchf) + 1, npk(nchf+1) - 1
c$$$                  v(kf-nchf,nchi,1:4) = v(kf,nchi,1)
c$$$               enddo 
c$$$            enddo
c$$$         enddo      

         
c$$$  if (packed) then
c$$$            do ki = 1, nd
c$$$               do kf = ki, nd
c$$$                  vmatp(ki+(kf-1)*kf/2) = - vmatp(ki+(kf-1)*kf/2) *
c$$$     >               sqrt(abs(real(wk(kf))*real(wk(ki))))
c$$$               enddo
c$$$            enddo
c$$$            do kf = 1, nd
c$$$               s = (real(wk(kf))+1e-30) / abs(real(wk(kf))+1e-30)
c$$$               vmatp(kf+(kf-1)*kf/2) = vmatp(kf+(kf-1)*kf/2) + s
c$$$            end do 
c$$$            do nchi = 1, nchtop
c$$$               do ki = npk(nchi) + 1, npk(nchi+1) - 1
c$$$                  do kf = ki, nd
c$$$                     vmatp(ki-nchi+(kf-1)*kf/2) = vmatp(ki+(kf-1)*kf/2)
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$            do nchf = 1, nchtop
c$$$               do kf = npk(nchf) + 1, npk(nchf+1) - 1
c$$$                  do ki = 1, kf - nchf
c$$$                     vmatp(ki+(kf-nchf-1)*(kf-nchf)/2) =
c$$$     >                  vmatp(ki+(kf-1)*kf/2)
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$            call date_and_time(date,time,zone,valuesin)
c$$$            print '("SSPSV entered at:",a10," with N:",i10)',time,nds
c$$$            call update(6)
c$$$            call sspsv('U',nds,nchtop,vmatp,iwork,v,nd,info)
c$$$            call date_and_time(date,time,zone,valuesout)
c$$$            print '("SSPSV exited at:",a10,", diff (secs):",i5)',
c$$$     >         time, idiff(valuesin,valuesout)
c$$$            call update(6)
c$$$            det = 0.0
c$$$            rcond = 0.0
c$$$         else
            delta = 0.0
            unv(:) = 0.0
            ndmax = 0
            inquire(file='delta',exist=exists)
            if (exists.and.lg.le.3) then
               open(42,file='delta')
               read(42,*) delta, ndmax
               print*,'Will be using DELTA, NDMAX:',delta, ndmax
               close(42)
               kn = 0
               do nch = 1, nchtop
                  do kg = 1, npk(nch+1)-npk(nch)
                     rk = gk(kg,nch)
                     kn = kn + 1
                     if (rk.gt.0.0) unv(kn) = 1.0 !rk * exp(-rk)
                  enddo
               enddo
               print*,'KN, ND:',kn,nd
            endif

            photon = projectile.eq.'photon'

c$$$            inquire(FILE="analytic",EXIST=analytic)  
c$$$
c$$$            print*,"analytic = ",analytic

            if (.false.) then !nbox.eq.1) then !analytic) then

c$$$               allocate( kernel(ichi,ichi))
c$$$
c$$$               open(42,file = "analytic")
c$$$               read(42,*) nbox
c$$$               close(42)
c$$$
c$$$               if (nbox .NE. 0) then
c$$$                  box = .TRUE.
c$$$                  print*, "box basis used"
c$$$               else
c$$$                  box = .FALSE.
c$$$                  print*, "true continuum functions used"
c$$$               endif
c$$$
c$$$               if (ns.eq.0) then
c$$$                  do ki = 1, nd
c$$$                     do kf = ki, nd
c$$$                        kernel(kf,ki) = vmat(kf,ki)
c$$$                        kernel(ki,kf) = vmat(kf,ki)
c$$$                     enddo
c$$$                  enddo
c$$$               else
c$$$                  do ki = 1, nd
c$$$                     do kf = ki, nd         
c$$$                        kernel(kf,ki) = vmat(ki,kf+1)
c$$$                        kernel(ki,kf) = vmat(ki,kf+1)
c$$$                     enddo
c$$$                  enddo
c$$$               endif

! Make V'=GVG
               do nchf = 1, nchtop
                  do nchn = 1, nchf !nchtop
! Extract V_{fn}
                     if (ns.eq.0) then
                        if (nchn.lt.nchf) then
                           do kf = npk(nchf)+1,npk(nchf+1)-1
                              do kn = npk(nchn)+1, npk(nchn+1)-1
                                 tmp4(kf-npk(nchf)+1,kn-npk(nchn)+1) =
     >                             vmat(kf,kn)
                              enddo
                           enddo
                        else ! nchn = nchf
                           do kf = npk(nchf)+1,npk(nchf+1)-1
                              do kn = npk(nchn)+1, kf
                                 tmp4(kf-npk(nchf)+1,kn-npk(nchn)+1) =
     >                             vmat(kf,kn)
                                 tmp4(kn-npk(nchn)+1,kf-npk(nchf)+1) =
     >                             vmat(kf,kn)
                              enddo
                           enddo
                        endif 
                     else
                        if (nchn.lt.nchf) then
                           do kf = npk(nchf)+1,npk(nchf+1)-1
                              do kn = npk(nchn)+1, npk(nchn+1)-1
                                 tmp4(kf-npk(nchf)+1,kn-npk(nchn)+1) =
     >                           vmat(kn,kf+1)
                              enddo
                           enddo
                        else ! nchn = nchf
                           do kf = npk(nchf)+1,npk(nchf+1)-1
                              do kn = npk(nchn)+1, kf
                                 tmp4(kf-npk(nchf)+1,kn-npk(nchn)+1) =
     >                           vmat(kn,kf+1)
                                 tmp4(kn-npk(nchn)+1,kf-npk(nchf)+1) =
     >                           vmat(kn,kf+1)
                              enddo
                           enddo
                        endif 
                     endif
! Define V'=GVG
c$$$                     do kf = npk(nchf)+1, npk(nchf+1)-1 !combine with below?
c$$$                        kgff = kf - npk(nchf) + 1
c$$$                        do kn = npk(nchn)+1, npk(nchn+1)-1
c$$$                           kgfn = kn - npk(nchn) + 1
c$$$                           do kp = 2, npk(nchn+1)-npk(nchn)
c$$$                              do kpp = 2, npk(nchf+1)-npk(nchf)
c$$$                                 tmp2(kgff,kgfn) = tmp2(kgff,kgfn) +
c$$$     >                              tmp4(kpp,kp)*gf(kgff,kpp,nchf)*
c$$$     >                              gf(kp,kgfn,nchn)
c$$$                              enddo
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP& SCHEDULE(dynamic)
C$OMP& PRIVATE(kgfn)
                     do kpp = 2, npk(nchf+1)-npk(nchf)
                        do kgfn = 2, npk(nchn+1)-npk(nchn)
                           tmp5(kpp,kgfn) = dot_product(
     >                        tmp4(kpp,2:npk(nchn+1)-npk(nchn)),
     >                        gf(2:npk(nchn+1)-npk(nchn),kgfn,nchn))
                        enddo
                     enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP& SCHEDULE(dynamic)
C$OMP& PRIVATE(kgfn)
                     do kgff = 2, npk(nchf+1)-npk(nchf)
                        do kgfn = 2, npk(nchn+1)-npk(nchn)
                           tmp2(kgff,kgfn) = dot_product(
     >                        gf(2:npk(nchf+1)-npk(nchf),kgff,nchf),
     >                        tmp5(2:npk(nchf+1)-npk(nchf),kgfn))
                        enddo
                     enddo
C$OMP END PARALLEL DO
! For each channel pair overwrite V with -V' for solving(G-V')K=V
                     if (ns.eq.0) then
                        do kf = npk(nchf)+1, npk(nchf+1)-1 
                           kgff = kf - npk(nchf) + 1
                           do kn = npk(nchn)+1, npk(nchn+1)-1
                              kgfn = kn - npk(nchn) + 1
                              if (kf.lt.kn) cycle
                              vmat(kf,kn) = - tmp2(kgff,kgfn)
                           enddo
                        enddo
                     else 
                        do kf = npk(nchf)+1, npk(nchf+1)-1 
                           kgff = kf - npk(nchf) + 1
                           do kn = npk(nchn)+1, npk(nchn+1)-1
                              kgfn = kn - npk(nchn) + 1
                              if (kf.lt.kn) cycle
                              vmat(kn,kf+1) = - tmp2(kgff,kgfn)
                           enddo
                        enddo
                     endif 
                  enddo
               enddo
C Redefine driving term V
               do nchi = 1, nchtop
                  do nchf = 1, nchtop
c$$$                     tmp3(:) = 0.0
                     do kf = npk(nchf)+1, npk(nchf+1)-1
                        kgff = kf - npk(nchf)+1
                        tmp3(kgff) = dot_product(  
     >                     gf(2:npk(nchf+1)-npk(nchf),kgff,nchf),
     >                     v(npk(nchf)+1:npk(nchf+1)-1,nchi,1))
c$$$                        do kp = npk(nchf)+1, npk(nchf+1)-1
c$$$                           kgfp = kp - npk(nchf) + 1
c$$$                           tmp3(kgff) = tmp3(kgff)+gf(kgff,kgfp,nchf)*
c$$$     >                        v(kp,nchi,1)
c$$$                        enddo
                     enddo
                     do kf = npk(nchf)+1, npk(nchf+1)-1
                        kgff = kf - npk(nchf)+1
                        v(kf,nchi,1:4) =  tmp3(kgff)
                     enddo
                     v(npk(nchf),nchi,1:4) = 0.0
                  enddo 
               enddo
C Add G to -V' on diagonal channels             
               do nch = 1, nchtop
                  if (ns.eq.0) then
                     vmat(npk(nch)+1:nd,npk(nch)) = 0.0
                     vmat(npk(nch),1:npk(nch)-1) = 0.0
                     vmat(npk(nch),npk(nch)) = 1.0
                     do kii = npk(nch)+1, npk(nch+1) - 1
                        do kff = npk(nch)+1, npk(nch+1) - 1
                           if (kff.lt.kii) cycle
                           vmat(kff,kii) = vmat(kff,kii) +
     >                        gf(kff-npk(nch)+1,kii-npk(nch)+1,nch)
                        enddo
                     enddo
                  else
                     vmat(npk(nch),npk(nch)+1+1:nd+1) = 0.0
                     vmat(1:npk(nch)-1,npk(nch)+1) = 0.0
                     vmat(npk(nch),npk(nch)+1) = 1.0
                     do kii = npk(nch)+1, npk(nch+1) - 1
                        do kff = npk(nch)+1, npk(nch+1) - 1
                           if (kff.lt.kii) cycle
                           vmat(kii,kff+1) = vmat(kii,kff+1) +
     >                        gf(kff-npk(nch)+1,kii-npk(nch)+1,nch)
                        enddo
                     enddo
                  endif 
               enddo
c$$$            call date_and_time(date,time,zone,valuesout)
c$$$            print '("added G at:",a10," diff (secs):",i5)',
c$$$     >         time,idiff(valuesin,valuesout)
c$$$            valuesin=valuesout
c$$$            call update(6)
               

c$$$C               call alexstatestest(gk,npk,nchtop,wk,weightk,nqm,ns
c$$$C     >                                                   ,kernel,box)
c$$$               zeff = zasym * nze
c$$$               temp(1:meshr) = zeff * 2.0 / rmesh(1:meshr,1)
c$$$c$$$                  temp(1:meshr) = - nze *
c$$$c$$$     >                 (- zasym * 2.0 / rmesh(1:meshr,1)) !+ui(i)) 
c$$$C                      + 2.0*vdcore(i,ll)
c$$$               print*,'Calling ALEXANALYTIC with Zeff:',zeff
c$$$               call alexanalytic(zeff,temp,lg,ichi,wk,gk,nchtop,npk,
c$$$     >            ns,v(1,1,2),kernel,box,chil,photon)
c$$$C               do nch = 1, nchtop
c$$$C                  do n = 1, npk(nchtop+1)-1
c$$$C                     v(n,nch,2) = -kernel(npk(nch),n)
c$$$C                  enddo
c$$$C               enddo 
c$$$C              Add the I matrix to -K
c$$$               do kf = 1, nd
c$$$!                  s = (real(wk(kf))+1e-30) / abs(real(wk(kf))+1e-30)
c$$$                  s = 1.0/(real(wk(kf))+1e-30)
c$$$                  kernel(kf,kf) = kernel(kf,kf) + s
c$$$               end do 
               if (sprint) then
                  do nchi = 1, min(9,nchtop)
                     do nchf = 1, min(9,nchtop)
                     open(42,file='splot.'//ch(nchf)//ch(nchi)//ch(ns))
            write(42,
     >         '("#   rkf       rki          vmat       kff kii")')
                     do ki = npk(nchi)+1, npk(nchi+1) - 1
                        kii = ki - npk(nchi) + 1
                        rki = gk(kii,nchi)
                        do kf = npk(nchf) + 1, npk(nchf+1) - 1
                           kff = kf - npk(nchf) + 1
                           rkf = gk(kff,nchf)
                           d = 1.0 ! rkf * rki
                           write(42,'(2f10.6,1p,e15.4,2i4)')
     >                        rkf,rki,vmat(kf,ki)/d,kff,kii
                        enddo
                        write(42,*)
c$$$                        write(42,'(1p,e10.2,4x,99e10.2)') v(ki,nchi,1),
c$$$     >                     (kernel(kf,ki),kf=npk(nchf),npk(nchf+1)-1)
                     enddo
                     close(42)
                     enddo 
                  enddo 
               endif 
            else ! original single node code
               if (ns.eq.0) then
                  do ki = 1, nd
                     do kf = ki, nd
                        vmat(kf,ki) = - vmat(kf,ki)
     >                 * (1.0+delta*unv(kf)*unv(ki))
!     >                 *  sqrt(abs(real(wk(kf))*real(wk(ki))))
C  Can comment out the above line and make a few more changes below
C  where sqrt(wk) occurs to get another way of solving the equations,
C  but this gives more ill-conditioned matrices
                     enddo
                  enddo
               else
                  do ki = 1, nd
                     do kf = ki, nd         
                        vmat(ki,kf+1) = - vmat(ki,kf+1)
     >                  * (1.0+delta*unv(kf+1)*unv(ki))
!     >                  * sqrt(abs(real(wk(kf))*real(wk(ki))))
                     enddo
                  enddo 
               endif
C  Add the I matrix to -K
               if (sprint) call splot(vmat,npk,nchtop,ns,gk)
               do kf = 1, nd
!                  s = (real(wk(kf))+1e-30) / abs(real(wk(kf))+1e-30)
                  s = 1.0/(real(wk(kf))+1e-30)
                  vmat(kf,kf+ns) = vmat(kf,kf+ns) + s
c$$$                  print*,'kf,v,vmat,wk:',kf,v(kf,1,1),vmat(kf,kf+ns),
c$$$     >                 real(wk(kf))
               end do 
c$$$               vmat(1,1+ns) = 1e30
c$$$               v(1,1,1) = 0.0
c$$$               print*,'vmat(1,1),v(1):',vmat(1,1+ns),v(1,1,1)
            endif               ! was for analytic, now set to .false. as immediate above works for all

            call date_and_time(date,time,zone,valuesout)
            print '("Kernel defined at:",a10," diff (secs):",i5)',
     >         time,idiff(valuesin,valuesout)
            valuesin=valuesout
            call update(6)


            call clock(s1)
c$$$            print*,'Mb allocated to work:', nw/1000000
c$$$            call update(6)
c$$$      call memalloc(ptrw,nw)
c$$$      if (ptrw.eq.0) stop 'Not enough memory for WORK'

c$$$            if (ns.eq.0) then
c$$$               do nchi = 1, nchtop
c$$$                  do ki = npk(nchi) + 1, npk(nchi+1) - 1
c$$$                     do kf = ki, nd
c$$$                        vmat(kf,ki-nchi) = vmat(kf,ki)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do nchf = 1, nchtop
c$$$                  do kf = npk(nchf) + 1, npk(nchf+1) - 1
c$$$                     do ki = 1, kf - nchf
c$$$                        vmat(kf-nchf,ki) = vmat(kf,ki)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            else
c$$$               do nchf = 1, nchtop
c$$$                  do kf = npk(nchf) + 1, npk(nchf+1) - 1
c$$$                     do ki = kf - nchf + 1, nd + 1
c$$$                        vmat(kf-nchf,ki) = vmat(kf,ki)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do nchi = 1, nchtop
c$$$                  do ki = npk(nchi) + 1, npk(nchi+1) - 1
c$$$                     do kf = 1, ki - 1
c$$$                        vmat(kf,ki-nchi+1) = vmat(kf,ki+1)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            endif
c$$$
c$$$            call date_and_time(date,time,zone,valuesin)
c$$$            print '("Kernel reduced at:",a10)',time
c$$$            call update(6)
            
c$$$      if (nsmax.eq.-10) then
c$$$         do ki = 1, nds
c$$$            do kf = ki, nds
c$$$              vmat(ki,kf) = vmat(kf,ki)
c$$$            enddo
c$$$         enddo
c$$$c$$$         call sgelss(nds,nds,nchtop,vmat,nd,v,nd,eigv,rcond,irank,
c$$$c$$$     >      work,lwork,info)
c$$$c$$$         do n = 1, nds
c$$$c$$$            print*,eigv(n)
c$$$c$$$         enddo 
c$$$c$$$
c$$$         do n = 1, nds
c$$$            iwork(n) = 0
c$$$         enddo
c$$$c$$$         print*,'enter RCOND please'
c$$$c$$$         read*,rcond
c$$$         rcond = 1e-14
c$$$         call sgelsx(nds,nds,nchtop,vmat,nd,v,nd,iwork,rcond,irank,
c$$$     >      work,info)
c$$$c$$$         print*,'INFO, RANK, NDS:',info,irank,nds
c$$$c$$$         call ssyev('v',ud(ns),nds,vmat,nd,eigv,work,lwork,info)
c$$$c$$$         do n = 1, nds
c$$$c$$$            if (abs(eigv(n)).lt.0.5) print*,n,eigv(n)
c$$$c$$$         enddo
c$$$c$$$         stop
c$$$c$$$c$$$         print*,'Enter Delta please'
c$$$c$$$c$$$         read*,delta
c$$$c$$$         delta = 1e-5
c$$$c$$$c$$$         call ssyevx('v','v',ud(ns),nds,vmat,nd,-delta,delta,1,1,0.0,
c$$$c$$$c$$$     >      neigv,eigv,v,nd,work,lwork,iwork,ifail,info)
c$$$c$$$         call clock(s2)
c$$$c$$$         print*,'time and INFO:', s2 - s1, info
c$$$c$$$         do i = 1, neigv
c$$$c$$$            print*,i,eigv(i),ifail(i)
c$$$c$$$         enddo
c$$$c$$$         nch = 1
c$$$c$$$         do while (gk(1,nch).gt.0.0.and.nch.le.nchtop)
c$$$c$$$            k = 2
c$$$c$$$            do while (gk(k,nch).lt.gk(1,nch))
c$$$c$$$               k = k + 1
c$$$c$$$            enddo
c$$$c$$$            iwork(nch) = npk(nch) + k - 1
c$$$c$$$            nch = nch + 1
c$$$c$$$         enddo
c$$$c$$$         nchtest = nch - 1
c$$$c$$$         do nchi = 1, nchtop
c$$$c$$$            do n = 1, nds
c$$$c$$$               v(n,nchi,1) = 0.0
c$$$c$$$            enddo 
c$$$c$$$            do nn = 1, nds
c$$$c$$$               ovlp = sdot(nds,v(1,nchi,2),1,vmat(1,nn),1)
c$$$c$$$               n = iwork(1)-1
c$$$c$$$               t1 = vmat(n,nn)/sqrt(abs(wk(n)))*real(wk(n))/abs(wk(n))
c$$$c$$$               n = iwork(1)
c$$$c$$$               t2 = vmat(n,nn)/sqrt(abs(wk(n)))*real(wk(n))/abs(wk(n))
c$$$c$$$               t = t1 * t2
c$$$c$$$               if ((abs(eigv(nn)).gt.delta.and.t.ge.0.0.and.
c$$$c$$$     >            abs(ovlp).gt.delta).or.abs(eigv(nn)).gt.0.1) then
c$$$c$$$                  coeff = ovlp / eigv(nn)
c$$$c$$$                  do n = 1, nds
c$$$c$$$                     v(n,nchi,1) = v(n,nchi,1) + vmat(n,nn) * coeff
c$$$c$$$                  enddo 
c$$$c$$$               else
c$$$c$$$                  if (ovlp.gt.delta) print*,'Skipping eigenvalue:',
c$$$c$$$     >               nchi,nn,eigv(nn),ovlp,t
c$$$c$$$c$$$                  if (t.gt.0.0) then
c$$$c$$$c$$$                     do nch = 1, nchtest                    
c$$$c$$$c$$$                        print*,(vmat(n,nn)/sqrt(abs(wk(n)))*real(wk(n))
c$$$c$$$c$$$     >                     /abs(wk(n)),n=iwork(nch)-2,iwork(nch)+1)
c$$$c$$$c$$$                     enddo
c$$$c$$$c$$$                  endif 
c$$$c$$$               endif
c$$$c$$$            enddo 
c$$$c$$$         enddo 
c$$$      else 

C     Solve the linear equations
         ifail = 0
         call date_and_time(date,time,zone,valuesin)
         print '("MATINV entered at:",a10," with N:",i10)',time,nds
         call update(6)
c$$$         call ssysv(ud(ns),nds,nchtop,vmat(1:nd,1+ns:nd+ns),nd,iwork,v,
c$$$     >      nd,work,lwork,info)

         lastns = ns.eq.nsmax
c$$$         if( analytic) then
c$$$            
c$$$            call alexmatinv(kernel,nd,nds,v,nchtop)
c$$$
c$$$            deallocate(kernel)
c$$$
c$$$            call date_and_time(date,time,zone,valuesout)
c$$$             print '("MATINV exited at:  ",a10,", diff (secs):",i5)',
c$$$     >            time, idiff(valuesin,valuesout)
c$$$            call update(6)
c$$$         else
            lwork = (256 + 3) * nd   !lwork = (64 + 3) * nd
            allocate(work(lwork))
            nw = 4 * lwork
            call matinv(vmat(1,1+ns),ud(ns),nd,nds,v,nchtop,work,
     >         lwork,delta,offwk,ndmax,det,rcond,ifail,lastns,packed)
            call date_and_time(date,time,zone,valuesout)
            print '("MATINV exited at:  ",a10,", diff (secs):",i5)',
     >         time, idiff(valuesin,valuesout)
            call update(6)
            deallocate(work)
c$$$         endif !analytic
c$$$      endif  ! packed

c$$$         do nchn = 1, nchtop
c$$$            do kn = 1, nd
c$$$               print*,kn,nchn,ns,v(kn,nchn,1),bb(kn,nchn),
c$$$     >            (v(kn,nchn,1)-bb(kn,nchn))/
c$$$     >            (v(kn,nchn,1)+bb(kn,nchn))
c$$$            enddo
c$$$         enddo
      endif ! scalapack: if true then using BB

c$$$  call memfree(ptrw)

c$$$  do nf = 2, nds
c$$$         print*,nf,v(nf,1,1)
c$$$      enddo 
      call clock(s2)
c$$$      print*,'Time for matinv:',s2-s1
c$$$      if (rc.gt.-1e-5) then
c$$$         call check(vmat,nchtop,nqm,wk,r,nds,work,v,err)
c$$$      endif
      do nchf = 1, nchtope2e
         nv = 0
         do kn = 1, nd
            cvec(kn) = !sqrt(abs(wk(kn))) *
     >         (ve2ed(nchf,kn) * dphasee2e(nchf) +
     >         r * ve2ee(nchf,kn) * ephasee2e(nchf) )
         enddo
         do nchn = 1, nchtop
            do kn = npk(nchn) + 1, npk(nchn+1) - 1
               cvec(kn-nchn) = cvec(kn)
            enddo 
         enddo 
         do nchi = 1, nchtop
            ke2e(nchf,nchi) = 
     >         ve2ed(nchf,npk(nchi)) * dphasee2e(nchf) +
     >         r * ve2ee(nchf,npk(nchi)) * ephasee2e(nchf)
            do kn = 1, nds
               ke2e(nchf,nchi) = ke2e(nchf,nchi) +
     >            cvec(kn) * v(kn,nchi,1)
            enddo 
            ke2e(nchf,nchi) = ke2e(nchf,nchi) / gk(1,nchi)
         enddo
         ve2e(nchf) = ve2ed(nchf,npk(1)) * dphasee2e(nchf)
     >      + r * ve2ee(nchf,npk(1)) * ephasee2e(nchf)
         ve2e(nchf) = ve2e(nchf) * phasel(1,1) / gk(1,1)
      enddo 

c$$$      call date_and_time(date,time,zone,valuesin)

      do nchf = 1, nchtop
c$$$         if (-lg.gt.latop) uba(nchf) = .true.
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
         do nchi = 1, nchtop
c$$$            if (nds.eq.0) then
c$$$               sum = 0.0
c$$$            else
c$$$c$$$               print*,'Calling sdot',nds,v(1,nchf,2),v(1,nchi,1)
c$$$               sum = sdot(nds,v(1,nchf,2),1,v(1,nchi,1),1)
c$$$c$$$               print*,'Exiting sdot', sum
c$$$            endif 
            ton(nchf,nchi) = ton(nchf,nchi) + 
     >         dot_product(v(1:nds,nchf,2),v(1:nds,nchi,1))
c$$$            print*,nchf,nchi,sum/gk(1,nchi)/gk(1,nchf),
c$$$     >         ton(nchf,nchi)/gk(1,nchi)/gk(1,nchf),
c$$$     >         abs(sum / (real(ton(nchf,nchi))+1e-5))
c$$$            if (nchi.lt.5) then
c$$$               uba(nchf) = uba(nchf) .and.
c$$$     >            abs(sum / (abs(ton(nchf,nchi))+1e-4)).lt.1e-2
c$$$            endif 

c$$$            if (gk(1,nchf).gt.0.0.and.gk(1,nchi).gt.0.0) then
c$$$            nk = ionshell(nchf)
c$$$            x1 = gk(nk-npk(nchf)-1,nchf)
c$$$            x2 = gk(nk-npk(nchf)+0,nchf)
c$$$            x3 = gk(nk-npk(nchf)+1,nchf)
c$$$            x4 = gk(nk-npk(nchf)+2,nchf)
c$$$            y1 = v(nk-nchf-2,nchi,1)/real(wk(nk-2))*sqrt(abs(wk(nk-2)))
c$$$            y2 = v(nk-nchf-1,nchi,1)/real(wk(nk-1))*sqrt(abs(wk(nk-1)))
c$$$            y3 = v(nk-nchf-0,nchi,1)/real(wk(nk-0))*sqrt(abs(wk(nk-0)))
c$$$            y4 = v(nk-nchf+1,nchi,1)/real(wk(nk+1))*sqrt(abs(wk(nk+1)))
c$$$            ton(nchf,nchi) = cubint(x1,y1,x2,y2,x3,y3,x4,y4,gk(1,nchf))
c$$$            endif 
c$$$            print*,nchf,nchi,x1,x2,gk(1,nchf),x3,x4
c$$$            print*,nchf,nchi,y1,y2,
c$$$     >         cubint(x1,y1,x2,y2,x3,y3,x4,y4,gk(1,nchf)),y3,y4
         enddo !nchi
C$OMP END PARALLEL DO

         if (lprint.and.nchf.le.74) then !74 corresponds to z
            if (nchf.lt.nfstart.or.nchf.gt.nfstop) cycle
            do nchi = 1, min(nchtop,74) !nchf
               if (nchi.lt.nistart.or.nchi.gt.nistop) cycle
C               write(lfile,'"lplot.",2i2) nchf,nchi
               open(52,file="lplot."//ch(nchf)//ch(nchi)//'_'//
     >            ch(mod(lg,10))//ch(ns)//ch(ipar))
            divk = gk(1,nchi) * gk(1,nchf)
            if (divk.eq.0.0) divk = 1.0
c$$$            write(nchi*100+nchf*10+ns,*) gk(1,nchf),
            write(52,'("# ",a3," <- ",a3,i4," <-",i3)') 
     >           chan(nchf), chan(nchi), nchf, nchi
            write(52,'("#    k               K(k,ki)         V(k,ki)")') 
            write(52,"(1p,5e16.5)") gk(1,nchf),             
     >         real(ton(nchf,nchi))/divk,
     >         real(von(nchf,nchi)/phasel(1,nchf)/phasel(1,nchi))
!     >        ,  real(wk(npk(nchf)))
c$$$     >         ,ovlpn(nchf), ovlpn(nchi)
C  The order is k, K(k), V(k), V(k) * K(k) / (E - En - k**2/2)
            do n = npk(nchf) + 1, npk(nchf+1) - 1
               kn = n - npk(nchf) + 1
               divk = gk(1,nchi) * gk(kn,nchf)
               vgf = dot_product(v(npk(nchf)+1:npk(nchf+1)-1,nchi,1),
     >            gf(2:npk(nchf+1)-npk(nchf),kn,nchf))
               if (divk.eq.0.0) divk = 1.0
               if (gk(kn,nchf).ne.0.0) then
c$$$                  write(nchi*100+nchf*10+ns,*) gk(kn,nchf),
                  if (scalapack) then
                     write(52,"(1p,5e16.5)") gk(kn,nchf),
c$$$     >                  v(n,nchi,1) / divk /
c$$$     >                  gf(kn-npk(nchf)+1,kn-npk(nchf)+1,nchf),!real(wk(n)),
!     >                  vgf / divk,
!     >                  v(n,nchi,1) / divk * gf(kn,kn,nchf),
     >                  vgf / divk,
     >                  v(n,nchi,2) / divk
     >                  ,gf(kn,kn,nchf)
                  else 
                     boxnorm = 1.0
                     if (nbox.eq.1) boxnorm = sqrt(2.0/rmesh(meshr,1))
                     write(52,"(1p,5e16.5)") gk(kn,nchf),
c$$$     >                  v(n,nchi,1) / divk / boxnorm
c$$$     >                  * gf(kn-npk(nchf)+1,kn-npk(nchf)+1,nchf) 
     >                  vgf / divk / boxnorm
!     >                  /  real(wk(n)) !* sqrt(abs(wk(n)))
     >                  ,v(n,nchi,2) /divk /boxnorm!/ sqrt(abs(wk(n)))
     >                  ,gf(kn,kn,nchf)
c$$$                  write(52,*) gk(kn,nchf),
c$$$     >               v(n-nchf,nchi,1) / divk /
c$$$     >               real(wk(n)) * sqrt(abs(wk(n))),
c$$$     >               v(n-nchf,nchi,2) / sqrt(abs(wk(n))) / divk
                  endif 
               else
                  divk = gk(1,nchi)
c$$$                  write(nchi*100+nchf*10+ns,*) gk(kn,nchf),
                  write(52,*) gk(kn,nchf),
     >               v(n,nchi,1) / divk,
     >               v(n,nchi,2) / divk
c$$$     >               v(n-nchf,nchi,1) / divk,
c$$$     >               v(n-nchf,nchi,2) / divk
               endif 
            enddo
            enddo 
c$$$            close(nchi*100+nchf*10+ns)
            close(52)
         endif 
      enddo
c$$$      call date_and_time(date,time,zone,valuesout)
c$$$      print '("On-shell K defined at:  ",a10,", diff (secs):",i5)',
c$$$     >   time, idiff(valuesin,valuesout)
c$$$      call update(6)
      endif ! end of not.converged

      nbad = 0
      do nchi = 1, nchtop
         if (gk(1,nchi).lt.0.0) cycle
         do nchf = nchi, nchtop
            if (gk(1,nchf).lt.0.0) cycle
            tmp =abs((ton(nchf,nchi) - ton(nchi,nchf))/
     >         (ton(nchf,nchi)+1e-30))
            if (tmp.gt.1e-2) then
               nbad = nbad + 1
c$$$               print*,'WARNING have a non symmetric',
c$$$     >            ' K matrix: K, F, I', real(ton(nchf,nchi)),
c$$$     >            real(ton(nchi,nchf)), nchf, nchi
            endif 
c$$$            ton(nchf,nchi) = (ton(nchf,nchi) + ton(nchi,nchf)) / 2.0
            ton(nchi,nchf) = ton(nchf,nchi)
c$$$            k1 = 0
c$$$            do n = npk(nchf) + 1, npk(nchf+1) - 1
c$$$               kn = n - npk(nchf) + 1
c$$$               divk = gk(1,nchi) * gk(kn,nchf)
c$$$               vgf = dot_product(v(npk(nchf)+1:npk(nchf+1)-1,nchi,1),
c$$$     >            gf(2:npk(nchf+1)-npk(nchf)-1,kn,nchf))
c$$$               if (gk(kn,nchf).lt.gk(1,nchf)) then
c$$$                  k1 = kn
c$$$                  offshellk1 = vgf / divk
c$$$               else
c$$$                  offshellk2 = vgf / divk
c$$$                  exit
c$$$               endif
c$$$            enddo
c$$$            f1 = dot_product(v(npk(nchf)+1:npk(nchf+1)-1,nchi,1),
c$$$     >         gf(2:npk(nchf+1)-npk(nchf)-1,kn-1,nchf))/
c$$$     >         gk(1,nchi)/ gk(kn-1,nchf)
c$$$            f2 = dot_product(v(npk(nchf)+1:npk(nchf+1)-1,nchi,1),
c$$$     >         gf(2:npk(nchf+1)-npk(nchf)-1,kn,nchf))/
c$$$     >         gk(1,nchi)/ gk(kn,nchf)
c$$$            f3 = dot_product(v(npk(nchf)+1:npk(nchf+1)-1,nchi,1),
c$$$     >         gf(2:npk(nchf+1)-npk(nchf)-1,kn+1,nchf))/
c$$$     >         gk(1,nchi)/ gk(kn+1,nchf)
c$$$            f4 = dot_product(v(npk(nchf)+1:npk(nchf+1)-1,nchi,1),
c$$$     >         gf(2:npk(nchf+1)-npk(nchf)-1,kn+2,nchf))/
c$$$     >         gk(1,nchi)/ gk(kn+2,nchf)
c$$$            call fourpointrule(
c$$$     >         gk(kn-1,nchf),f1,
c$$$     >         gk(kn  ,nchf),f2,
c$$$     >         gk(kn+1,nchf),f3,
c$$$     >         gk(kn+2,nchf),f4,       
c$$$     >         gk(1,nchf),f,df)
c$$$            slope = (offshellk2-offshellk1)/(gk(k1+1,nchf)-gk(k1,nchf))
c$$$            b = offshellk1-slope*gk(k1,nchf)
c$$$            est = slope*gk(1,nchf)+b
c$$$            print"('nchi,nchf,k1, offshellk1:',1p,2i4,2e11.3)",
c$$$     >         nchi,nchf,gk(k1,nchf),offshellk1
c$$$            print"('nchi,nchf,kon,  onshellk:',1p,2i4,6e11.3,e9.1,'%')",
c$$$     >         nchi,nchf,gk(1,nchf),
c$$$     >         real(ton(nchf,nchi))/(gk(1,nchi) * gk(1,nchf)),est,f,
c$$$     >         slope,df,abs((real(ton(nchf,nchi))/gk(1,nchi)/gk(1,nchf)-
c$$$     >         est)/est)*100.0
c$$$            print"('nchi,nchf,k2, offshellk2:',1p,2i4,2e11.3)",
c$$$     >         nchi,nchf,gk(k1+1,nchf),offshellk2
         enddo
      enddo 
      print'(
     >   "Nonsymmetric (to 1%) K elements:",i5," out of",i6,i3,"%")', 
     >   nbad, max(nds,1),nbad*100/max(nds,1)
      
      if (sprint) then
      open (42,file='pwkmat'//ch(ns))
      write(42,'(''# K-matrix elements for J ='',i3)') lg
      do nchi = 1, nchtop
         call getchinfo (nchi,nchip,lg,temp,maxpsi, ei, lai, ni, li)
         if (etot.gt.enchan(nchip)) then
         write(42,'(''#nchf Lf  lf  nf <-nchi  Li  li  ni Re(K-mat)'',
     >      ''        ef       ovlp'')')
            do nchf = nchi, nchtop
               call getchinfo (nchf,nchfp,lg,temp,maxpsf,ef,laf,nf,lf)
               if (etot.gt.enchan(nchfp)) write
     >            (42,'(4i4,'' <-'',4i4,1p,e12.4,0p,f12.5,f10.5)') 
     >            nchf,lf,laf,nf,nchi,li,lai,ni, real(ton(nchf,nchi))
     >            /gk(1,nchi)/gk(1,nchf),ef,enchandiff(nchfp)!ovlpn(nchfp)
c$$$               if (ei.gt.etot/2.0.or.ef.gt.etot/2.0) then
c$$$                  ton(nchf,nchi) = 0.0
c$$$                  ton(nchi,nchf) = 0.0
c$$$               endif 
            enddo
         endif 
      enddo
      close(42)
      endif

C  Check for Alisher
c$$$      print*,'Setting K=0 for e > E/2'
c$$$      do nchi = 1, nchtop
c$$$         call getchinfo (nchi,nchip,lg,temp,maxpsi, ei, lai, ni, li)
c$$$         do nchf = nchi, nchtop
c$$$            call getchinfo (nchf,nchfp,lg,temp,maxpsf,ef,laf,nf,lf)
c$$$            if (ei.gt.etot/2.0.or.ef.gt.etot/2.0) then
c$$$               ton(nchf,nchi) = 0.0
c$$$               ton(nchi,nchf) = 0.0
c$$$            endif 
c$$$         enddo
c$$$      enddo 
      
C  Get the T-matrix
c$$$      do nchi = 1, nchtop
c$$$         do nchf = 1, nchtop
c$$$            print"('ns,nchi,nchf,gki,gkf,Kon:',3i2,2f10.6,e12.4)", 
c$$$     >        ns,nchi,nchf,gk(1,nchi),gk(1,nchf),real(ton(nchf,nchi))
c$$$            if (gk(1,nchf).lt.0.0.or.gk(1,nchi).lt.0.0) then
c$$$c$$$               if (c.ne.0.0) print*,'C should be zero:',c
c$$$               c = 0.0
c$$$            else
c$$$               c = 1.0
c$$$            endif
c$$$            const = gk(1,nchi) * gk(1,nchf)
c$$$            if (const.eq.0.0) const = 1.0
c$$$            kon(nchf,nchi) = ton(nchf,nchi)/const
c$$$C  The variable KON is the K matrix K(f,i). To get the T matrix we use
c$$$C  Sum_i T_fi(X_f,X_i) * conjg(phi(X_f)) * conjg(phi(X_i)) *
c$$$C                 (I(i,i') + i * pi X_i * K_ii'(X_i,X_i')) = K_fi'(X_f,X_i').
c$$$C  As I, K and T are all symmetric we write this as (f <-> i')
c$$$C  Sum_i (I(f,i) + i * pi X_i * K_fi(X_f,X_i)) *
c$$$C   * T_ii'(X_i,X_i') * conjg(phi(X_i')) * conjg(phi(X_i)) = K_fi'(X_f,X_i'),
c$$$            ckern(nchf,nchi) = c * kon(nchf,nchi) *
c$$$     >         cmplx(0.0, pi*gk(1,nchi))
c$$$C  Maybe multiply sum by c in the next line also to get 0 for closed
c$$$C  channels when printing out NCHTOP T matrix elements
c$$$            ton(nchf,nchi) = c * kon(nchf,nchi)
c$$$         end do 
c$$$      end do

      ckern(:,:) = (0.0,0.0)
      if (gk(1,1).eq.0.0) then
         kon(1,1) = ton(1,1)
      else 
         do nchi = 1, nchtop
            gki = 1.0
            if (gk(1,nchi).gt.0.0) gki = sqrt(gk(1,nchi))
            do nchf = 1, nchtop
c$$$            print"('ns,nchi,nchf,gki,gkf,Kon:',3i2,2f10.6,e12.4)", 
c$$$     >        ns,nchi,nchf,gk(1,nchi),gk(1,nchf),real(ton(nchf,nchi))
               if (gk(1,nchf).lt.0.0.or.gk(1,nchi).lt.0.0) then
                  kon(nchf,nchi) = 0.0
               else
                  gkf = 1.0
                  if (gk(1,nchf).gt.0.0) gkf = sqrt(gk(1,nchf))
                  kon(nchf,nchi) = ton(nchf,nchi)/gki/gkf
               endif
C  The variable KON is the K matrix K(f,i). To get the T matrix we use
C  Sum_i T_fi(X_f,X_i) * conjg(phi(X_f)) * conjg(phi(X_i)) *
C                 (I(i,i') + i * pi X_i * K_ii'(X_i,X_i')) = K_fi'(X_f,X_i').
C  As I, K and T are all symmetric we write this as (f <-> i')
C  Sum_i (I(f,i) + i * pi X_i * K_fi(X_f,X_i)) *
C   * T_ii'(X_i,X_i') * conjg(phi(X_i')) * conjg(phi(X_i)) = K_fi'(X_f,X_i'),
               ckern(nchf,nchi) = kon(nchf,nchi) * cmplx(0.0, pi) 
c$$$     >            * gk(1,nchi)/(gk(1,nchi)+1e-10)
               ton(nchf,nchi) = kon(nchf,nchi)
            end do 
         end do
      endif 
      do nch = 1, nchtop
         ckern(nch,nch) = ckern(nch,nch) + (1.0,0.0)
      end do

C  Solve the linear equations
      if (projectile.eq.'photon') then
         call date_and_time(date,time,zone,valuesin)
         rmcv = 1e-6*8 * (ichi) * nchan
         rmcv2 = 1e-6*8 * (ichi) * nchtop
c$$$         call memalloc(ptrcv,mcv)
c$$$         call memalloc(ptrcv2,mcv2)

         print'(''Requested memory (Mb) for CV and CV2:'',f10.1)',
     >      (rmcv+rmcv2)
c$$$         if (ptrcv.eq.0) stop 'Not enough memory for CV'
c$$$         if (ptrcv2.eq.0) stop 'Not enough memory for CV2'
         allocate(cv(nchan,ichi))
         allocate(cv2(ichi,nchtop))
         if (analytic) then
            call analyticphoto(dr,npk,nchtop)
            call analyticphoto(dv,npk,nchtop)
            call analyticphoto(dx,npk,nchtop)
         endif 
         do nchf = 1, nchtop
            c = 1.0
            if (gk(1,nchf).lt.0.0) c = 0.0
            gkf = 1.0
            if (gk(1,nchf).gt.0.0) gkf = sqrt(gk(1,nchf))
            do nchi = 1, nchtop
               cv(nchf,npk(nchi)) = ton(nchf,nchi)
               nn = 1
               if (scalapack) then
                  do n = npk(nchi)+1, npk(nchi+1)-1
                     nn = nn + 1
c$$$                  cv(nchf,n) = v(n-nchi,nchf,1)/gk(nn,nchi)/gk(1,nchf)/
                     cv(nchf,n) = c * v(n,nchf,1)/sqrt(abs(gk(1,nchf)))
     >                  /real(wk(n))
                  enddo
               else 
                  do n = npk(nchi)+1, npk(nchi+1)-1
                     nn = nn + 1
c$$$  cv(nchf,n) = v(n-nchi,nchf,1)/gk(nn,nchi)/gk(1,nchf)/
c$$$                     cv(nchf,n) = c * v(n,nchf,1)/gk(nn,nchi)/gk(1,nchf)
                     cv(nchf,n) = c * v(n,nchf,1)/gkf !sqrt(abs(gk(1,nchf)))
     >                  /real(wk(n)) !* sqrt(abs(wk(n)))
                  enddo
               endif 
            enddo
         enddo 
         nv = ichi
         nd = nchtop
         call matinv2(ckern,nchan,nd,cv,nv,cwork,erfp,epsil)
         do nchi = 1, nchtop
            gki = 1.0
            if (gk(1,nchi).gt.0.0) gki = sqrt(gk(1,nchi))
            do nchf = 1, nchtop
               gkf = 1.0
               if (gk(1,nchf).gt.0.0) gkf = sqrt(gk(1,nchf))
               ton(nchf,nchi)=cv(nchf,npk(nchi))!*gk(1,nchi)*gk(1,nchf)
               nn = 0
               do n = npk(nchf), npk(nchf+1) - 1
                  nn = nn + 1
                  divk = sqrt(abs(gk(1,nchi)))
                  cv2(n,nchi) = cv(nchi,n) * gki !divk
               enddo
               cv2(npk(nchf),nchi)=cv2(npk(nchf),nchi)*gkf
c$$$     >            sqrt( abs(gk(1,nchf))) ! * gk(1,nchf)))
c$$$               print*,'nchi,nchf,ton:',nchi,nchf,
c$$$     >            ton(nchf,nchi) *sqrt( gk(1,nchi) * gk(1,nchf))
c$$$               print*,'nchi,nchf,cv2(1),ki,kf:',nchi,nchf,
c$$$     >            cv2(npk(nchf),nchi),gk(1,nchi),gk(1,nchf)
c$$$               print*,'nchi,nchf,cv2(2),ki,kf:',nchi,nchf,
c$$$     >            cv2(npk(nchf)+1,nchi)!, gk(1,nchi) , gk(1,nchf)
            enddo
         enddo
         
         if(lg.eq.1 .and. ns.eq.iSpin)
     >      call PHOTO(gk,npk,wk,cv2,nchtop,ovlpn,csfile,slowery,
     >      phasel,phaseq)

C (gamma,3e) by double shake-off
!         if(lg.eq.0 .and. ns.eq.ispin)
!     >      call shakeoff(gk,npk,wk,cv2,nchtop,ovlpn,csfile,
!     .      slowery,phasel)

C (e,3e) agenda
c$$$         call E2E(lg,gk,npk,wk,cv2,nchtop,ovlpn,csfile,slowery,phasel,
c$$$     >      phaseq)
c$$$         if(lg.eq.0 .or. lg.eq.2)
c$$$     >   call E2E2(lg,gk,npk,wk,cv2,nchtop,ovlpn,csfile,slowery,phasel,
c$$$     >      phaseq)
!         call E2Eq(lg,gk,npk,wk,cv2,nchtop,ovlpn,csfile,slowery,phasel)

         deallocate (cv, cv2)
c$$$         call memfree(ptrcv)
c$$$         call memfree(ptrcv2)
         call date_and_time(date,time,zone,valuesout)
         print '("Photo exited at:  ",a10,", diff (secs):",i5)',
     >      time, idiff(valuesin,valuesout)
         call update(6)
      else
         nv = nchtop
         nd = nchtop
         call matinv2(ckern,nchan,nd,ton,nv,cwork,erfp,epsil)
c$$$         call date_and_time(date,time,zone,valuesout)
c$$$         print '("MATINV2 returned at:  ",a10,", diff (secs):",i5)',
c$$$     >      time, idiff(valuesin,valuesout)
c$$$         call update(6)
      endif 
c$$$      call memfree(ptrv)

      deallocate (v, stat=istat)
      print*,'Deallocated V:', istat

      do nchi = 1, 1
         do nchf = 1, nchtope2e
            te2e(nchf) = ke2e(nchf,nchi)
            do nchn = 1, nchtop
               c = cmplx(0.0,pi * max(gk(1,nchn),0.0))
               te2e(nchf) = te2e(nchf) - c *
     >            ke2e(nchf,nchn) * ton(nchn,nchi)
            enddo
            te2e(nchf) = te2e(nchf) * phasel(1,nchi)
         enddo
      enddo
 
      do nchi = 1, nchtop
         gki = 1.0
         if (gk(1,nchi).gt.0.0) gki = sqrt(gk(1,nchi))
         do nchf = 1, nchtop
            gkf = 1.0
            if (gk(1,nchf).gt.0.0) gkf = sqrt(gk(1,nchf))
            ton(nchf,nchi) = ton(nchf,nchi) * phasel(1,nchf) *
     >         phasel(1,nchi)/gki/gkf
c$$$            print"('ns,nchi,nchf,gki,gkf,Ton:',3i2,2f10.6,2e12.4)", 
c$$$     >        ns,nchi,nchf,gk(1,nchi),gk(1,nchf),ton(nchf,nchi)
         enddo
C  In the case of Spin = 1.5 (ns=1) we musn't add TDIST to those TMATs which are 0 (singlet states in He-like)
C     
         if (abs(von(nchi,nchi)).gt.0.0 .or. ns.eq.0) then
c$$$            print*,'TON,TDIST:',ton(nchi,nchi), tdist(nchi)
            ton(nchi,nchi) = ton(nchi,nchi) + tdist(nchi)
            von(nchi,nchi) = von(nchi,nchi) + tdist(nchi)
         endif 
      enddo
c$$$      if (abs(slowery*2.0-etot).lt.1e-2)
c$$$     >   call gette2e(ton,kon,ovlpn,gk(1,1:nchtop),nchtop,lg,etot,
c$$$     >   slowery,ns)
C  The following checks unitarity relations for every T-matrix element!      D.Fursa
c$$$      do nchf = 1, nchtop
c$$$         do nchi = 1, nchtop
c$$$            c = (0.0,0.0)
c$$$            do nchn = 1, nchtop
c$$$               c = c + conjg(ton(nchf,nchn))*ton(nchn,nchi)*gk(1,nchn)
c$$$            enddo
c$$$            print*,nchf,nchi,c,-imag(ton(nchf,nchi))/pi
c$$$         enddo
c$$$      enddo 

c$$$C  The following inverts the T matrix to get the on-shell K matrix
c$$$C  It should always come out to be essentially real
c$$$      do nchi = 1, nchtop
c$$$         do nchf = 1, nchtop
c$$$            if (gk(1,nchf).lt.0.0.or.gk(1,nchi).lt.0.0) then
c$$$               c = 0.0
c$$$            else
c$$$               c = 1.0/sigma(nchf)/sigma(nchi)
c$$$            endif
c$$$            kon(nchf,nchi) = c * ton(nchf,nchi)
c$$$            ckern(nchf,nchi) = - c * ton(nchf,nchi) *
c$$$     >         cmplx(0.0, pi*gk(1,nchi))
c$$$         end do
c$$$         ckern(nchi,nchi) = ckern(nchi,nchi) + (1.0,0.0)
c$$$      enddo
c$$$      nd = nchtop
c$$$      nv = nchtop
c$$$      call matinv2(ckern,nchan,nd,kon,nv,cwork,erfp,epsil)
c$$$      if (sprint) then
c$$$         open (42,file='kmat'//ch(ns))
c$$$         write(42,'(''# K-matrix elements for J ='',i3)') lg
c$$$      endif 
c$$$      do nchi = 1, nchtop
c$$$         call getchinfo (nchi,nchip,lg,temp,maxpsi, ei, lai, ni, li)
c$$$         if (etot-enchan(nchip).ge.0.0) then
c$$$c$$$         write(42,'('' Lf lf nf <- Li li ni   Re(K-mat)   Im(K-mat)'',
c$$$c$$$     >      ''        ef       ovlp'')')
c$$$            do nchf = 1, nchtop
c$$$               call getchinfo (nchf,nchfp,lg,temp,maxpsf,ef,laf,nf,lf)
c$$$               deta = -1/sqrt(abs(ef))
c$$$               gkfac = gk(1,nchi)*gk(1,nchf)
c$$$               if (gkfac.eq.0.0) gkfac = 1.0
c$$$               if (etot-enchan(nchfp).ge.0.0.and.sprint) write
c$$$     >            (42,'(1p,10e10.2)') ef, ei, kon(nchf,nchi),
c$$$     >             ton(nchf,nchi),ovlpn(nchfp) * ovlpn(nchip)
c$$$c$$$     >            (42,'(3i3,'' <-'',3i3,1p,2e12.4,0p,f12.5,f10.5)') 
c$$$c$$$     >            lf,laf,nf,li,lai,ni, kon(nchf,nchi),ef,ovlpn(nchfp)
c$$$               if (nchf.ge.nchi) then
c$$$                  kon(nchf,nchi)=cmplx( real(kon(nchf,nchi)),
c$$$     >               k2nd(nchf,nchi)/gkfac)
c$$$                  kon(nchi,nchf) = kon(nchf,nchi)
c$$$               endif 
c$$$            enddo
c$$$            if (sprint) write(42,*)
c$$$         endif 
c$$$      enddo
c$$$      if (sprint) close(42)
c$$$      print'('' JS f i  abs(T)    arg(T)       K        K2nd   '
c$$$  >   //'   real(V)   imag(V)     Ef(Ry)'')'
      print'('' JS f i  real(T)    imag(T)     K        K2nd   '
     >   //'   real(V)   imag(V)     Ef(Ry)'')'
      return
      end
      
      subroutine fourpointrule(rk1,f1,rk2,f2,rk3,f3,rk4,f4,rk,f,df)
      c1 = f1/((rk1-rk2)*(rk1-rk3)*(rk1-rk4))
      c2 = f2/((rk2-rk1)*(rk2-rk3)*(rk2-rk4))
      c3 = f3/((rk3-rk1)*(rk3-rk2)*(rk3-rk4))
      c4 = f4/((rk4-rk1)*(rk4-rk3)*(rk4-rk2))
      f =  c1*(rk-rk2)*(rk-rk3)*(rk-rk4)
     >   + c2*(rk-rk1)*(rk-rk3)*(rk-rk4)
     >   + c3*(rk-rk1)*(rk-rk2)*(rk-rk4)
     >   + c4*(rk-rk1)*(rk-rk3)*(rk-rk2)      
      df = c1*(((rk-rk3)*(rk-rk4))+(rk-rk2)*(rk-rk4)+(rk-rk3)*(rk-rk2))+
     >   c2*(((rk-rk3)*(rk-rk4))+(rk-rk1)*(rk-rk4)+(rk-rk3)*(rk-rk1))+
     >   c3*(((rk-rk1)*(rk-rk4))+(rk-rk2)*(rk-rk4)+(rk-rk2)*(rk-rk1))+
     >   c4*(((rk-rk3)*(rk-rk2))+(rk-rk1)*(rk-rk2)+(rk-rk3)*(rk-rk1))      
      return
      end
      
      subroutine gette2e(ton,kon,ovlpn,gk,nchtop,lg,etot,slowery,ns)
      include 'par.f'
      complex ton(nchan,nchan),te2e(6*(lamax+1),nchan),
     >   kon(nchan,nchan),c,a(nchan,nchan),b(nchan,6*(lamax+1)),
     >   workc(nchan*64)
      complex*16 coulphase
      real gk(nchan),ovlpn(nchan),temp(maxr), imag, xsec(nchan),
     >   ate2e(6*(lamax+1),nchan),konr(nchan,nchan),w(nchan),
     >   work(nchan*3)
      real*8 xin(KNM,6*(lamax+1)),yin(KNM,6*(lamax+1),3),
     >   xout(1),yout(1,3)
      integer nin(6*(lamax+1)),ntc(0:lamax,-lamax:lamax,2)
      character chan(knm)*3,chanf*3
      common /charchan/ chan
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      logical canstop
      type channel
      real ea, k, ovlp
      integer la,l,isa,nt
      end type channel
      type(channel) chinfo(nchan)
      type lchannel
      integer la,l,is
      end type lchannel
      type(lchannel) lchan(6*(lamax+1))
      data pi/3.1415927/


      ntc(:,:,:) = 0
      nt = 0
      do nch = 1, nchtop
         call getchinfo (nch, nchp, lg, temp, maxpsi, ea, la, na, l)
         if (chan(nchp)(1:1).eq.' '.or.chan(nchp)(1:1).eq.'s') then
            isa = 1
         else if (chan(nchp)(1:1).eq.'t') then
            isa = 2
         else
            stop 'Have not coded for such states in gette2e'
         endif
         if (ntc(la,lg-l,isa).eq.0) then
            nt = nt + 1
            ntc(la,lg-l,isa) = nt
            lchan(ntc(la,lg-l,isa)) = lchannel(la,l,isa)
         endif 
         chinfo(nch) = channel(ea,gk(nch),ovlpn(nchp),la,l,isa,
     >      ntc(la,lg-l,isa))
      enddo
      ntmax = nt

c$$$      do nchi = 1, nchtop
c$$$         do nchf = 1, nchtop
c$$$            konr(nchf,nchi) = real(kon(nchf,nchi))
c$$$         enddo
c$$$      enddo
c$$$      lwork=nchan*3
c$$$      call ssyev('V','U',nchtop,konr,nchan,w,work,lwork,info)
c$$$      if (info.eq.0) print*,'used and optimal LWORK:',lwork,work(1)
c$$$
c$$$      if (ns.eq.0) then
c$$$         open(42,file='v0')
c$$$      else
c$$$         open(42,file='v1')
c$$$      endif 
c$$$      do nchi = 1, nchtop
c$$$         do nchf = 1, nchtop
c$$$            sum = 0.0
c$$$            do nchn = 1, nchtop
c$$$               sum = sum + konr(nchf,nchn)*w(nchn)*konr(nchi,nchn)
c$$$            enddo
c$$$            write(42,*) chinfo(nchf)%ea, nchi, real(ton(nchf,nchi)),
c$$$     >         imag(ton(nchf,nchi))
c$$$c$$$            print*,'KON,SUM:',kon(nchf,nchi),sum
c$$$         enddo
c$$$         write(42,*)
c$$$      enddo
c$$$      close(42)
      
      xout(1) = slowery
      nout = 1
      do nchi = 1, nchtop
         nin(:) = 1
         do nchf = 1, nchtop
            ea = chinfo(nchf)%ea
            nt = chinfo(nchf)%nt
            ovlp = chinfo(nchf)%ovlp
            if (ea.gt.0.0.and.ea.lt.etot) then
               nin(nt) = nin(nt) + 1
               xin(nin(nt),nt) = ea
               yin(nin(nt),nt,1) = real(ton(nchf,nchi)) * ovlp * ea
               yin(nin(nt),nt,2) = imag(ton(nchf,nchi)) * ovlp * ea
               yin(nin(nt),nt,3) =  abs(ton(nchf,nchi)) * ovlp * ea
            endif
         enddo
         do nt = 1, ntmax
            yin(1,nt,1) = 0.0
            yin(1,nt,2) = 0.0
            xin(1,nt) = 0.0
            nin(nt) = nin(nt) + 1
            xin(nin(nt),nt) = etot
            yin(nin(nt),nt,1) = 0.0
            yin(nin(nt),nt,2) = 0.0
            do i = 1, 3
               call intrpl(nin(nt),xin(1,nt),yin(1,nt,i),
     >            nout,xout,yout(1,i))
            enddo 
            te2e(nt,nchi) = cmplx(yout(1,1),yout(1,2)) / slowery
            ate2e(nt,nchi) = yout(1,3) / slowery 
         enddo
      enddo

      
c$$$      it = 0
c$$$      canstop = .false.
c$$$      do while (it.lt.20.and..not.canstop)
c$$$         it = it + 1
c$$$!         canstop = .true.
c$$$         do nt = 1, ntmax
c$$$            do nchi = 1, nchtop
c$$$               c(nt,nchi) = (0.0,0.0)
c$$$               do nchn = 1, nchtop
c$$$                  c(nt,nchi) = c(nt,nchi) + chinfo(nchn)%k *
c$$$     >               conjg(te2e(nt,nchn))*ton(nchn,nchi)
c$$$               enddo
c$$$            enddo
c$$$            sum1 = 0.0
c$$$            sum2 = 0.0
c$$$            r = 1.0
c$$$            do nchi = 1, nchtop
c$$$               sum1 = sum1 + abs(imag(c(nt,nchi)))
c$$$               sum2 = sum2 + imag(c(nt,nchi))
c$$$!               canstop = canstop.and.abs(imag(c)/real(c+1e-20)).lt.1e-4
c$$$               parti = -pi * real(c(nt,nchi))/r+
c$$$     >            (r-1.0)*imag(te2e(nt,nchi))/r
c$$$               diff = ate2e(nt,nchi)*ate2e(nt,nchi)-parti*parti
c$$$               if (diff.lt.0.0) then
c$$$                  partr = 0.0
c$$$               else 
c$$$                  partr = sqrt(diff) * real(te2e(nt,nchi)) /
c$$$     >               abs(real(te2e(nt,nchi))+1e-20)
c$$$               endif
c$$$               te2e(nt,nchi) = cmplx(partr,parti)
c$$$            enddo
c$$$            print*,'it,nt,sum1,sum2:',it,nt,sum1,sum2!,te2e(nt,1)
c$$$         enddo
c$$$      enddo 
c$$$      print*,'Iterations after exit:',it
      
C  The following checks unitarity relations for every Te2e-matrix element!
      do nt = 1, ntmax
         sum = 0.0
         do nchi = 1, nchtop
            c = (0.0,0.0)
            do nchn = 1, nchtop
               c = c + conjg(te2e(nt,nchn)) * ton(nchn,nchi) *
     >            chinfo(nchn)%k
            enddo
            sum = sum + abs(pi*c+imag(te2e(nt,nchi)))
!            print*,nt,nchi,-c*pi,imag(te2e(nt,nchi))
         enddo
         la = lchan(nt)%la
         l  = lchan(nt)%l
         print'("Unitarity check 1:",1p,e8.1,", J,ns,la,l,e2eT:",
     >      4i3,2e11.3)',sum/nchtop,lg,ns,la,l,te2e(nt,1)/
     >         sqrt(sqrt(slowery)) *
     >         coulphase(dble(-(zasym+1.0)/sqrt(slowery)),la)
      enddo 

      do nchi = 1, nchtop
         do nchf = 1, nchtop
            a(nchf,nchi) = - pi * ton(nchf,nchi) * chinfo(nchi)%k
         enddo
         do nt = 1, ntmax
            b(nchi,nt) = cmplx(imag(te2e(nt,nchi)),0.0)
         enddo 
      enddo
      
c$$$      do nchf = 1, 1 !nchtop
c$$$         do nchi = 1, nchtop
c$$$            c = (0.0,0.0)
c$$$            do nchn = 1, nchtop
c$$$               c= c +
c$$$     >            a(nchi,nchn) * conjg(te2e(nchf,nchn))
c$$$            enddo
c$$$            print*,nchf,nchi,c,
c$$$     >         b(nchi)!,c
c$$$         enddo
c$$$      enddo
      
      nrhs=ntmax
      lwork = nchan*64
c$$$      call cgels('N',nchtop,nchtop,nrhs,a,nchan,b,nchan,workc,
c$$$     >   lwork,info)
c$$$      if (info.eq.0) then
c$$$         print*,'used and optimal LWORK:',lwork,workc(1)
c$$$      else
c$$$         print*,'WARNING: INFO from CGELS is:',info
c$$$      endif 
c$$$      do nt = 1, ntmax
c$$$         sum = 0.0
c$$$         do nchi = 1, nchtop
c$$$            c = (0.0,0.0)
c$$$            do nchn = 1, nchtop
c$$$               c = c +
c$$$     >            b(nchn,nt)*ton(nchn,nchi)*
c$$$     >            chinfo(nchn)%k
c$$$            enddo
c$$$            sum = sum + abs(imag(b(nchi,nt))-c*pi)
c$$$c$$$            print*,nt,nchi,c*pi, imag(b(nchi,nt))
c$$$         enddo
c$$$         la = lchan(nt)%la
c$$$         l  = lchan(nt)%l
c$$$         print'("Unitarity check 2:",1p,e8.1,", J,ns,la,l,e2eT:",
c$$$     >      4i3,2e11.3)',sum/nchtop,lg,ns,la,l,conjg(b(1,nt))/
c$$$     >      sqrt(sqrt(slowery)) *
c$$$     >      coulphase(dble(-(zasym+1.0)/sqrt(slowery)),la)
c$$$      enddo 
c$$$      do nt = 1, ntmax
c$$$         la = lchan(nt)%la
c$$$         print*,'Interpolated T for channel 1:',la,
c$$$     >      te2e(nt,1) /
c$$$     >         sqrt(sqrt(slowery)) *
c$$$     >         coulphase(dble(-(zasym+1.0)/sqrt(slowery)),la)
c$$$      enddo 
               
      return
      end
      
      function cubint(x1,y1,x2,y2,x3,y3,x4,y4,x)
      cubint = y4*(x-x1)*(x-x2)*(x-x3)/(x4-x1)/(x4-x2)/(x4-x3) +
     >   y3*(x-x1)*(x-x2)*(x-x4)/(x3-x1)/(x3-x2)/(x3-x4) +
     >   y2*(x-x1)*(x-x4)*(x-x3)/(x2-x1)/(x2-x4)/(x2-x3) +
     >   y1*(x-x4)*(x-x2)*(x-x3)/(x1-x4)/(x1-x2)/(x1-x3)
      return
      end
      
      subroutine makev(vmat,kf,nchf,nchtop,ichi,npk,wk,ns,vec,
     >   vmatp,packed)
      include 'par.f'
      integer npk(nchtop+1)
      real vmat(ichi,ichi+1),vec(*),vmatp(*)!vmatp(ichi*(ichi+1)/2)
      complex wk(*)
      logical packed

      kff = npk(nchf) + kf - 1
c$$$      if (packed) then
c$$$         do nchn = nchf, nchtop
c$$$            do knn = npk(nchn) + 1, npk(nchn+1) - 1
c$$$               vec(knn) = vmatp(kff+(knn-1)*knn/2) * sqrt(abs(real(wk(knn)))) 
c$$$            end do
c$$$         end do
c$$$         do nchn = 1, nchf - 1
c$$$            do knn = npk(nchn) + 1, npk(nchn+1) - 1
c$$$               vec(knn) = vmatp(knn+(kff-1)*kff/2) * sqrt(abs(real(wk(knn)))) 
c$$$            end do
c$$$         end do
c$$$      else         

C        analytic
         if (ns.eq.0) then
            do nchn = nchf, nchtop
               do knn = npk(nchn), npk(nchn+1) - 1
                  vec(knn) = vmat(knn,kff) !* sqrt(abs(real(wk(knn)))) 
               end do
            end do
            do nchn = 1, nchf - 1
               do knn = npk(nchn), npk(nchn+1) - 1
                  vec(knn) = vmat(kff,knn) !* sqrt(abs(real(wk(knn)))) 
               end do
            end do
         else
            do nchn = nchf, nchtop 
               do knn = npk(nchn), npk(nchn+1) - 1
                  vec(knn) = vmat(kff,knn+1) !* sqrt(abs(real(wk(knn)))) 
               end do
            end do
            do nchn = 1, nchf - 1
               do knn = npk(nchn), npk(nchn+1) - 1
                  vec(knn) = vmat(knn,kff+1) !* sqrt(abs(real(wk(knn)))) 
               end do
            end do
         endif
 
C         if (ns.eq.0) then
C            do nchn = nchf, nchtop
C               do knn = npk(nchn) + 1, npk(nchn+1) - 1
C                  vec(knn) = vmat(knn,kff)! * sqrt(abs(real(wk(knn)))) 
C               end do
C            end do
C            do nchn = 1, nchf - 1
C               do knn = npk(nchn) + 1, npk(nchn+1) - 1
C                  vec(knn) = vmat(kff,knn)! * sqrt(abs(real(wk(knn)))) 
C               end do
C            end do
C         else
C            do nchn = nchf, nchtop 
C               do knn = npk(nchn) + 1, npk(nchn+1) - 1
C                  vec(knn) = vmat(kff,knn+1)! * sqrt(abs(real(wk(knn)))) 
C               end do
C            end do
C            do nchn = 1, nchf - 1
C               do knn = npk(nchn) + 1, npk(nchn+1) - 1
C                  vec(knn) = vmat(knn,kff+1)! * sqrt(abs(real(wk(knn)))) 
C               end do
C            end do
C         endif
c$$$      endif
      end
      
      function onshellk(vmat,nchtop,nchf,nchi,ns,npk,ichi,vmatp,packed)
      include 'par.f'
      integer npk(nchtop+1)
      real vmat(ichi,ichi+1),vmatp(*)!vmatp(ichi*(ichi+1)/2)
      logical packed

      kff = npk(nchf)
      kii = npk(nchi)
c$$$      if (packed) then
c$$$         if (nchf.ge.nchi) then
c$$$            onshellk = vmatp(kii+(kff-1)*kff/2)
c$$$         else 
c$$$            onshellk = vmatp(kff+(kii-1)*kii/2)
c$$$         endif 
c$$$      else       
         if (ns.eq.0) then
            if (nchf.ge.nchi) then
               onshellk = vmat(npk(nchf),npk(nchi))
            else 
               onshellk = vmat(npk(nchi),npk(nchf))
            endif 
         else
            if (nchf.ge.nchi) then
               onshellk = vmat(npk(nchi),npk(nchf)+1)
            else 
               onshellk = vmat(npk(nchf),npk(nchi)+1)
            endif
         endif 
c$$$      end if
      end
      
c$$$      function onshellv(vmat,nqm,nchtop,nchf,nchi,r)
c$$$      include 'par.f'
c$$$      real vmat(nqm,nqm,0:nchtop,nchtop)
c$$$c$$$      real vmat(kmax,kmax,0:nchan,nchan)
c$$$
c$$$      if (nchf.ge.nchi) then
c$$$         onshellv = vmat(1,1,nchf,nchi) + r * vmat(1,1,nchi-1,nchf) 
c$$$      else 
c$$$         onshellv = vmat(1,1,nchi,nchf) + r * vmat(1,1,nchf-1,nchi) 
c$$$      end if
c$$$      end
      
c$$$      subroutine check(vmat,nchtop,nqm,wk,r,nd,vec,v,err)
c$$$      include 'par.f'
c$$$      real vmat(nqm,nqm,0:nchtop,nchtop),v(nqm*nchtop,nchtop,2)
c$$$c$$$      real vmat(kmax,kmax,0:nchan,nchan),v(kmax*nchan,nchan,2)
c$$$     >   ,vec(kmax*nchan),work(kmax*nchan),
c$$$     >   err(nchan,nchan),anorm(nchan)
c$$$      complex wk(kmax,nchan)
c$$$      
c$$$      do nchi = 1, nchtop
c$$$         anorm(nchi) = 0.0
c$$$c$$$         call makevec(vmat,1,nchi,nqm,nchtop,wk,r,v(1,nchi,2))
c$$$      enddo 
c$$$      nn = 0
c$$$      do nchf = 1, nchtop
c$$$         do kf = 1, nqm
c$$$            call makevec(vmat,kf,nchf,nqm,nchtop,
c$$$     >         wk,r,vec)
c$$$            w = (real(wk(kf,nchf))+1e-30)
c$$$            nn = nn + 1
c$$$            do nchi = 1, nchtop
c$$$               sum = sdot(nd,vec,1,v(1,nchi,1),1)
c$$$c$$$               sum = 0.0
c$$$c$$$               do n = 1, nd
c$$$c$$$                  sum = sum + vec(n) * v(n,nchi)
c$$$c$$$               enddo
c$$$c$$$               anorm = anorm + abs(- v(nn,nchi)/w +
c$$$c$$$     >            sum + work(nn))
c$$$               if (kf.gt.1)
c$$$     >            anorm(nchi)=anorm(nchi)+(-v(nn,nchi,1)*sqrt(abs(w))/w
c$$$     >            + sum + v(nn,nchi,2) / sqrt(abs(w))) ** 2
c$$$            enddo 
c$$$         enddo
c$$$         do nchi = 1, nchtop
c$$$            err(nchf,nchi) = sqrt(anorm(nchi))
c$$$            anorm(nchi) = 0.0
c$$$         enddo 
c$$$      enddo 
c$$$      end
         
c$$$      subroutine makedriv(vmat,nchtop,nqm,wk,r,nd,c,vec,v)
c$$$      include 'par.f'
c$$$      real vmat(nqm,nqm,0:nchtop,nchtop),v(nqm*nchtop,nchtop,2)
c$$$c$$$      real vmat(kmax,kmax,0:nchan,nchan),v(kmax*nchan,nchan,2)
c$$$     >   ,vec(kmax*nchan)
c$$$      complex wk(kmax,nchan)
c$$$      
c$$$      nn = 0
c$$$      do nchf = 1, nchtop
c$$$         do kf = 2, nqm
c$$$            call makevec(vmat,kf,nchf,nqm,nchtop,wk,r,vec)
c$$$            nn = nn + 1
c$$$            do nchi = 1, nchtop
c$$$               sum = sdot(nd,vec,1,v(1,nchi,1),1)
c$$$c$$$               sum = 0.0
c$$$c$$$               do n = 1, nd
c$$$c$$$                  sum = sum + vec(n) * v(n,nchi,1) 
c$$$c$$$               enddo
c$$$c$$$               v(nn,nchi,2) = sum * c
c$$$               v(nn,nchi,2) = sum * sqrt(abs(wk(kf,nchf))) * c
c$$$            enddo
c$$$         enddo
c$$$      enddo 
c$$$      end
c$$$
C  The following routine solves the set of linear equations where the
      SUBROUTINE MATINV(KERNEL,ud,LDA,N,V,m,work,lwork,delta,offwk,
     >   ndmax,determ,rcond,i,lastns,packed)
      include 'par.f'
      parameter (nmax = kmax * nchan)
c$$$      real, dimension(n*n) :: afi
c$$$      pointer (ptrafi,afi)
c$$$      real, dimension(:), allocatable :: afi
c$$$      real kernel(lda,lda+1)
      real kernel(lda,lda)
      real*8 fact
      REAL V,work,DET,determ(2),offwk(lda)
      REAL ERFP,EPSIL
      DIMENSION V(lda,m,4),work(lwork),
     >   ipvt(nmax),inert(3),iwork(nmax),berr(nchan),ferr(nchan)
      common/smallr/ formcut,regcut,expcut,fast,match
      character ud, date*8,time*10,zone*5
      integer valuesin(8)
      logical fast, match, packed, lastns

      if (n.le.0) return
c$$$      if (lda.ne.nmax) stop 'LDA must equal NMAX in MATINV'
      IA=LDA
      IB=LDA
      IJOB=i
c$$$      if (fast) then
      print*,'packed,lastns,ndmax:',packed,lastns,ndmax
      if (packed.or..not.lastns.or.ndmax.gt.0) then
         call date_and_time(date,time,zone,valuesin)
         print '("XSYSV entered at: ",a10)',time
#ifdef _single
      
         call ssysv(ud,N,m,kernel,LDA,ipvt,v,lda,work,lwork,info)

#elif defined _double

         call dsysv(ud,N,m,kernel,LDA,ipvt,v,lda,work,lwork,info)

#endif
         call date_and_time(date,time,zone,valuesin)
         print '("XSYSV exited at:  ",a10)',time
         if (lwork.lt.work(1).or.nint(work(1)/n).ne.64)
     >      print*,'previous and optimal LWORK, NB:',lwork,work(1),
     >      work(1)/N
         if (delta.ne.0.0.and.ndmax.gt.0) then
            v(1:n,:,4) = v(1:n,:,1)
            do nd = 1, ndmax
               do i = 1, m
                  v(1:n,i,1) = (v(1:n,i,1) *
     >               abs(offwk(1:n))/(offwk(1:n)+1e-20) - v(1:n,i,2)) *
     >               nd / (1.0+delta)
                  v(1:n,i,2) = v(1:n,i,1)
               enddo 
#ifdef _single
      
               CALL SSYTRS(ud,N,m,kernel,LDA,ipvt,v(1,1,1),LDA,INFO)

#elif defined _double

               CALL DSYTRS(ud,N,m,kernel,LDA,ipvt,v(1,1,1),LDA,INFO)

#endif
               v(1:n,:,4) = v(1:n,:,4)+(-delta)**nd*v(1:n,:,1)/fact(nd)
            enddo
            v(:,:,1) = v(:,:,4)
            v(:,:,2) = v(:,:,3)
            call date_and_time(date,time,zone,valuesin)
            print '("DSYTRS exited at:  ",a10)',time
         endif 
         rcond = 0.0
         anorm =0.0
c$$$         print*,'Entered SLANSY'
c$$$         call update(6)
C         ANORM = SLANSY( 'I', ud, N, kernel, LDA, WORK )
c$$$         print*,'Entered SSYCON'
c$$$         call update(6)
C         CALL SSYCON( ud, N, kernel, LDA, ipvt, ANORM, RCOND, WORK, 
C     >      IWORK, INFO )
      else
         print*,'UD:',ud
c$$$         do i = 1, n
c$$$            print"(i3,1p,100e8.1)",i,(kernel(i,j),j=1,n+1)
c$$$         enddo
         if (ud.eq.'L'.or.ud.eq.'l') then
            do i = 2, n
               do j = 1, i - 1
                  kernel(j,i) = kernel(i,j)
               enddo
            enddo
         elseif (ud.eq.'U'.or.ud.eq.'u') then
            do i = 1, n
               do j = i+1, n
                  kernel(j,i) = kernel(i,j)
               enddo
            enddo
         else
            stop 'ud.ne. u/U/l/L'
         endif 
         call date_and_time(date,time,zone,valuesin)
         print '("XGESV entered at: ",a10)',time

#ifdef _single
      
         call sgesv(N,m,kernel,LDA,ipvt,v,lda,info)

#elif defined _double

         call dgesv(N,m,kernel,LDA,ipvt,v,lda,info)

#endif

         call date_and_time(date,time,zone,valuesin)
         print '("XGESV exited at:  ",a10)',time
c$$$         mafi = 4 * N * N
c$$$         call memalloc(ptrafi,mafi)
c$$$         if (ptrafi.eq.0) stop 'Not enough memory for AFI'

c$$$         call date_and_time(date,time,zone,valuesin)
c$$$         print '("SSYSVX entered at: ",a10)',time
c$$$         allocate (afi(n*n))
c$$$         call ssysvx('N',ud,N,m,kernel,LDA,AFI,N,ipvt,v(1,1,2),LDA,
c$$$     >      v,LDA,rcond,ferr,berr,work,lwork,iwork,info)
c$$$         print*,'Info,Rcond from SSYSVX:',info,rcond
c$$$         deallocate(afi)
c$$$         call date_and_time(date,time,zone,valuesin)
c$$$         print '("SSYSVX exited at:  ",a10)',time
c$$$  call memfree(ptrafi)
      endif       
      if (info.ne.0) print*,'Warning: INFO, N:', info, N
      END

      SUBROUTINE MATINV2(API,N1,NMAX,BHAR,nb,WK,ERFP,EPSIL)
      include 'par.f'
      COMPLEX API,BHAR,WK,DET(2)
      REAL ERFP,EPSIL
      DIMENSION API(N1,N1),BHAR(N1,nchan),WK(N1),kpvt(kmax*nchan)
      
      lda = n1
      n = nmax
! Due to occasional SEGV on Magnus, downloaded local copies. Igor 24/1/2018
#ifdef _single
c$$$      print*,'Calling own version of CGESV'      
      call cgesv(n,nb,api,lda,kpvt,bhar,lda,info)

#elif defined _double
c$$$      print*,'Calling own version of ZGESV'
      call zgesv(n,nb,api,lda,kpvt,bhar,lda,info)
#endif
      print*,'Exiting ?GESV'

c$$$      call cgeco(api,lda,n,kpvt,rcond,wk)
c$$$c$$$      call cgefa(api,lda,n,kpvt,info)
c$$$      do m = 1, nb
c$$$         call cgesl(api,lda,n,kpvt,bhar(1,m),0)
c$$$      enddo 
c$$$      call cgedi(api,lda,n,kpvt,det,wk,10)
C  The following Linpack routines assume that API is a symmetric kernel
c$$$      call csifa(api,lda,n,kpvt,info)
c$$$      call csisl(api,lda,n,kpvt,bhar)
c$$$      call csidi(api,lda,n,kpvt,det,wk,10)
      if (info.ne.0) then
         print*,'STOPPING since INFO in MATINV2:',info
         stop 'STOPPING since INFO in MATINV2 is nonzero'
      endif 
      return
C  The IMSL routine assumes a general kernel
c$$$      IA=N1
c$$$      IB=N1
c$$$      M=1
c$$$      IJOB=0
c$$$      CALL CLEQT(API,NMAX,IA,BHAR,M,IB,IJOB,WK,IER)
c$$$      DET(1) = (1.0,0.0)
c$$$      DET(2) = (0.0,0.0)
c$$$      DO 5 I = 1,NMAX
c$$$         IPVT = WK(I)
c$$$         IF (IPVT .NE. I) DET(1) = -DET(1)
c$$$         DET(1) = DET(1)*API(I,I)
c$$$         do while (abs(det(1)).gt.1e1)
c$$$            det(1) = det(1) / 1e1
c$$$            det(2) = det(2) + (1.0,0.0)
c$$$         enddo 
c$$$         do while (abs(det(1)).lt.1.0)
c$$$            det(1) = det(1) * 1e1
c$$$            det(2) = det(2) - (1.0,0.0)
c$$$         enddo 
c$$$ 5    CONTINUE
c$$$      if (ier.ne.0) print*,'IER not zero in CLEQT',ier
c$$$      IF (IER.EQ.129) THEN
c$$$         PRINT*,'MATRIX IS ALGORITHMICALLY SINGULAR.'
c$$$         RETURN
c$$$      END IF
c$$$      IF (IER.EQ.130) PRINT*,'MATRIX IS ILL CONDITIONED.'
      END

      real function sdotpr(n,sx,incx,sy,incy)
c
c     forms the dot product of two vectors with minimal precision loss
c
      real sx(n),sy(n),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
C Non-precision loosing dot product
      s = 0.0
      c = 0.0
      do j = 1, n
         y = c + sx(j) * sy(j)
         t = s + y
         f = 0.0
         if (y*s.ge.0.0) f = (0.46 * t - t) + t
         c = ((s - f) - (t - f)) + y
         s = t
      enddo
      sdotpr = s + c
      return
      end

c$$$      subroutine makebasis(m,vmat,r,nqm,nchtop,wk,z)
c$$$      include 'par.f'
c$$$      parameter (nmax = kmax * nchan)
c$$$      real vmat(nqm,nqm,0:nchtop,nchtop)
c$$$c$$$      real vmat(kmax,kmax,0:nchan,nchan),z(kmax*nchan,nchan,2)
c$$$     >   , z(nqm*nchtop,nchtop,2),
c$$$     >   vec(kmax*nchan), x(kmax*nchan,nchan),
c$$$     >   mat(nchan,nchan), b(nchan), kernel(nchan,nchan),
c$$$     >   bmat(nchan), vec2(kmax*nchan)
c$$$      complex wk(kmax,nchan), wdum(kmax,nchan)
c$$$      integer ipiv(nchan)
c$$$      data wdum /nmax * (1.0,0.0)/
c$$$      save mat, x, b
c$$$
c$$$      nd = (nqm - 1) * nchtop
c$$$      if (m.eq.1) then
c$$$         ki = 1
c$$$         nchi = 1
c$$$         call makevec(vmat,ki,nchi,nqm,nchtop,wdum,r,vec)
c$$$      else
c$$$         do n = 1, nd
c$$$            vec(n) = x(n,m-1)
c$$$         enddo 
c$$$         do j = 1, m - 1
c$$$            c = sdot(nd,z(1,j,1),1,x(1,m-1),1)
c$$$c$$$            c = 0.0
c$$$c$$$            do np = 1, nd
c$$$c$$$               c = c + z(np,j,1) * x(np,m-1)
c$$$c$$$            enddo
c$$$            do n = 1, nd
c$$$              vec(n) = vec(n) - z(n,j,1) * c
c$$$            enddo
c$$$         enddo
c$$$      endif 
c$$$
c$$$      veclen2 = sdot(nd,vec,1,vec,1)
c$$$      veclen = sqrt(veclen2)
c$$$      do n = 1, nd
c$$$         z(n,m,1) = vec(n) / veclen
c$$$      enddo
c$$$
c$$$      np = 0
c$$$      do nchn = 1, nchtop
c$$$         do kn = 2, nqm
c$$$            np = np + 1
c$$$            vec2(np) = real(wk(kn,nchn)) * z(np,m,1)
c$$$         enddo
c$$$      enddo
c$$$      
c$$$      n = 0
c$$$      do nchf = 1, nchtop
c$$$         do kf = 2, nqm
c$$$            call makevec(vmat,kf,nchf,nqm,nchtop,wdum,r,vec)
c$$$            n = n + 1
c$$$            x(n,m) = sdot(nd,vec2,1,vec,1)
c$$$         enddo
c$$$      enddo
c$$$      
c$$$      nchi = 1
c$$$      ki = 1
c$$$      call makevec(vmat,ki,nchi,nqm,nchtop,wdum,r,vec)
c$$$      b(m) = sdot(nd,z(1,m,1),1,vec,1)
c$$$c$$$      b(m) = 0.0
c$$$c$$$      do np = 1, nd
c$$$c$$$         b(m) = b(m) + z(np,m,1) * vec(np)
c$$$c$$$      enddo
c$$$      
c$$$      do mi = 1, m
c$$$         do mf = 1, m
c$$$            if (mf.eq.m.or.mi.eq.m) then
c$$$               mat(mf,mi) = sdot(nd,z(1,mf,1),1,x(1,mi),1)
c$$$c$$$               mat(mf,mi) = 0.0
c$$$c$$$               do n = 1, nd
c$$$c$$$                  mat(mf,mi) = mat(mf,mi) + z(n,mf,1) * x(n,mi)
c$$$c$$$               enddo
c$$$            endif
c$$$            kernel(mf,mi) = - mat(mf,mi)
c$$$         enddo
c$$$         kernel(mi,mi) = kernel(mi,mi) + 1.0
c$$$         bmat(mi) = b(mi)
c$$$      enddo
c$$$      nb = 1
c$$$      lda = nchan
c$$$      call sgesv(m,nb,kernel,lda,ipiv,bmat,lda,info)
c$$$      nchi = 1
c$$$      do n = 1, nd
c$$$         z(n,nchi,2) = 0.0
c$$$         do j = 1, m
c$$$            z(n,nchi,2) = z(n,nchi,2) + bmat(j) * z(n,j,1)
c$$$         enddo
c$$$      enddo
c$$$      n = 0
c$$$      do nchf = 1, nchtop
c$$$         do kf = 2, nqm
c$$$            n = n + 1
c$$$            z(n,nchi,2) = z(n,nchi,2) * sqrt(abs(wk(kf,nchf))) *
c$$$     >         real(wk(kf,nchf)) / abs(wk(kf,nchf))
c$$$         enddo
c$$$      enddo 
c$$$c$$$      mi = m
c$$$c$$$      do mf = max(mi - 2,1), mi
c$$$c$$$         c = 0.0
c$$$c$$$         do n = 1, nd
c$$$c$$$            c = c + z(n,mi) * z(n,mf)
c$$$c$$$         enddo
c$$$c$$$         print*,'mf,mi,ovlp',mf,mi,c
c$$$c$$$      enddo 
c$$$      return
c$$$      end
      
      subroutine splot(vmat,npk,nchtop,ns,gk)
      include 'par.f'
      integer npk(nchtop+1)
      dimension gk(kmax,nchan),vmat(npk(nchtop+1)-1,npk(nchtop+1))
      character ch*1
      ch(n) = char(n + ichar('0'))

      maxnch = 100
      do nchi = 1, min(nchtop,maxnch)
         do nchf = nchi, min(nchtop,maxnch)
            open(42,file='splot.'//ch(nchf)//ch(nchi)//ch(ns))
            write(42,
     >         '("#   rkf       rki          vmat       kff kii")')
            do ki = npk(nchi) + 0, npk(nchi+1) - 1
               kii = ki - npk(nchi) + 1
               rki = gk(kii,nchi)
               do kf = npk(nchf) + 0, npk(nchf+1) - 1
                  kff = kf - npk(nchf) + 1
                  rkf = gk(kff,nchf)
                  d = 1.0 ! rkf * rki
                  if (ns.eq.0) then
                     if (kf.ge.ki) then
                        write(42,'(2f10.6,e15.4,2i4)')
     >                     rkf,rki,vmat(kf,ki)/d,kff,kii
                     else
                        write(42,'(2f10.6,e15.4,2i4)')
     >                     rkf,rki,vmat(ki,kf)/d,kff,kii
                     endif
                  else 
                     if (kf.ge.ki) then
                        write(42,'(2f10.6,e15.4,2i4)')
     >                     rkf,rki,vmat(ki,kf+1)/d,kff,kii
                     else
                        write(42,'(2f10.6,e15.4,2i4)')
     >                     rkf,rki,vmat(kf,ki+1)/d,kff,kii
                     endif
                  endif
               enddo
               write(42,*)
            enddo
            close(42)
         enddo
      enddo
      return
      end
                     


      subroutine optpot
C  This routine is not used in CCC, and so is not loaded. For CCO the routine
C  is in the file de.f
      end
      
c$$$      SUBROUTINE CGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c$$$*
c$$$*  -- LAPACK driver routine (version 3.0) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c$$$*     Courant Institute, Argonne National Lab, and Rice University
c$$$*     March 31, 1993 
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            INFO, LDA, LDB, N, NRHS
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX            A( LDA, * ), B( LDB, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGESV computes the solution to a complex system of linear equations
c$$$*     A * X = B,
c$$$*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
c$$$*
c$$$*  The LU decomposition with partial pivoting and row interchanges is
c$$$*  used to factor A as
c$$$*     A = P * L * U,
c$$$*  where P is a permutation matrix, L is unit lower triangular, and U is
c$$$*  upper triangular.  The factored form of A is then used to solve the
c$$$*  system of equations A * X = B.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The number of linear equations, i.e., the order of the
c$$$*          matrix A.  N >= 0.
c$$$*
c$$$*  NRHS    (input) INTEGER
c$$$*          The number of right hand sides, i.e., the number of columns
c$$$*          of the matrix B.  NRHS >= 0.
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N)
c$$$*          On entry, the N-by-N coefficient matrix A.
c$$$*          On exit, the factors L and U from the factorization
c$$$*          A = P*L*U; the unit diagonal elements of L are not stored.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  IPIV    (output) INTEGER array, dimension (N)
c$$$*          The pivot indices that define the permutation matrix P;
c$$$*          row i of the matrix was interchanged with row IPIV(i).
c$$$*
c$$$*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
c$$$*          On entry, the N-by-NRHS matrix of right hand side matrix B.
c$$$*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
c$$$*
c$$$*  LDB     (input) INTEGER
c$$$*          The leading dimension of the array B.  LDB >= max(1,N).
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value
c$$$*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
c$$$*                has been completed, but the factor U is exactly
c$$$*                singular, so the solution could not be computed.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CGETRF, CGETRS, XERBLA
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      INFO = 0
c$$$      IF( N.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( NRHS.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -4
c$$$      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -7
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGESV ', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Compute the LU factorization of A.
c$$$*
c$$$      CALL CGETRF( N, N, A, LDA, IPIV, INFO )
c$$$      IF( INFO.EQ.0 ) THEN
c$$$*
c$$$*        Solve the system A*X = B, overwriting B with X.
c$$$*
c$$$         CALL CGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
c$$$     $                INFO )
c$$$      END IF
c$$$      RETURN
c$$$*
c$$$*     End of CGESV
c$$$*
c$$$      END
c$$$      SUBROUTINE CGETRF( M, N, A, LDA, IPIV, INFO )
c$$$*
c$$$*  -- LAPACK routine (version 3.0) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c$$$*     Courant Institute, Argonne National Lab, and Rice University
c$$$*     September 30, 1994
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            INFO, LDA, M, N
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX            A( LDA, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGETRF computes an LU factorization of a general M-by-N matrix A
c$$$*  using partial pivoting with row interchanges.
c$$$*
c$$$*  The factorization has the form
c$$$*     A = P * L * U
c$$$*  where P is a permutation matrix, L is lower triangular with unit
c$$$*  diagonal elements (lower trapezoidal if m > n), and U is upper
c$$$*  triangular (upper trapezoidal if m < n).
c$$$*
c$$$*  This is the right-looking Level 3 BLAS version of the algorithm.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  M       (input) INTEGER
c$$$*          The number of rows of the matrix A.  M >= 0.
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The number of columns of the matrix A.  N >= 0.
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N)
c$$$*          On entry, the M-by-N matrix to be factored.
c$$$*          On exit, the factors L and U from the factorization
c$$$*          A = P*L*U; the unit diagonal elements of L are not stored.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,M).
c$$$*
c$$$*  IPIV    (output) INTEGER array, dimension (min(M,N))
c$$$*          The pivot indices; for 1 <= i <= min(M,N), row i of the
c$$$*          matrix was interchanged with row IPIV(i).
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value
c$$$*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
c$$$*                has been completed, but the factor U is exactly
c$$$*                singular, and division by zero will occur if it is used
c$$$*                to solve a system of equations.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      COMPLEX            ONE
c$$$      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      INTEGER            I, IINFO, J, JB, NB
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CGEMM, CGETF2, CLASWP, CTRSM, XERBLA
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      INTEGER            ILAENV
c$$$      EXTERNAL           ILAENV
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX, MIN
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      INFO = 0
c$$$      IF( M.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
c$$$         INFO = -4
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGETRF', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      IF( M.EQ.0 .OR. N.EQ.0 )
c$$$     $   RETURN
c$$$*
c$$$*     Determine the block size for this environment.
c$$$*
c$$$      NB = ILAENV( 1, 'CGETRF', ' ', M, N, -1, -1 )
c$$$      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
c$$$*
c$$$*        Use unblocked code.
c$$$*
c$$$         CALL CGETF2( M, N, A, LDA, IPIV, INFO )
c$$$      ELSE
c$$$*
c$$$*        Use blocked code.
c$$$*
c$$$         DO 20 J = 1, MIN( M, N ), NB
c$$$            JB = MIN( MIN( M, N )-J+1, NB )
c$$$*
c$$$*           Factor diagonal and subdiagonal blocks and test for exact
c$$$*           singularity.
c$$$*
c$$$            CALL CGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
c$$$*
c$$$*           Adjust INFO and the pivot indices.
c$$$*
c$$$            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
c$$$     $         INFO = IINFO + J - 1
c$$$            DO 10 I = J, MIN( M, J+JB-1 )
c$$$               IPIV( I ) = J - 1 + IPIV( I )
c$$$   10       CONTINUE
c$$$*
c$$$*           Apply interchanges to columns 1:J-1.
c$$$*
c$$$            CALL CLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
c$$$*
c$$$            IF( J+JB.LE.N ) THEN
c$$$*
c$$$*              Apply interchanges to columns J+JB:N.
c$$$*
c$$$               CALL CLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
c$$$     $                      IPIV, 1 )
c$$$*
c$$$*              Compute block row of U.
c$$$*
c$$$               CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
c$$$     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
c$$$     $                     LDA )
c$$$               IF( J+JB.LE.M ) THEN
c$$$*
c$$$*                 Update trailing submatrix.
c$$$*
c$$$                  CALL CGEMM( 'No transpose', 'No transpose', M-J-JB+1,
c$$$     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
c$$$     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
c$$$     $                        LDA )
c$$$               END IF
c$$$            END IF
c$$$   20    CONTINUE
c$$$      END IF
c$$$      RETURN
c$$$*
c$$$*     End of CGETRF
c$$$*
c$$$      END
c$$$      SUBROUTINE CGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c$$$*
c$$$*  -- LAPACK routine (version 3.0) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c$$$*     Courant Institute, Argonne National Lab, and Rice University
c$$$*     September 30, 1994
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      CHARACTER          TRANS
c$$$      INTEGER            INFO, LDA, LDB, N, NRHS
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX            A( LDA, * ), B( LDB, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGETRS solves a system of linear equations
c$$$*     A * X = B,  A**T * X = B,  or  A**H * X = B
c$$$*  with a general N-by-N matrix A using the LU factorization computed
c$$$*  by CGETRF.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  TRANS   (input) CHARACTER*1
c$$$*          Specifies the form of the system of equations:
c$$$*          = 'N':  A * X = B     (No transpose)
c$$$*          = 'T':  A**T * X = B  (Transpose)
c$$$*          = 'C':  A**H * X = B  (Conjugate transpose)
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The order of the matrix A.  N >= 0.
c$$$*
c$$$*  NRHS    (input) INTEGER
c$$$*          The number of right hand sides, i.e., the number of columns
c$$$*          of the matrix B.  NRHS >= 0.
c$$$*
c$$$*  A       (input) COMPLEX array, dimension (LDA,N)
c$$$*          The factors L and U from the factorization A = P*L*U
c$$$*          as computed by CGETRF.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  IPIV    (input) INTEGER array, dimension (N)
c$$$*          The pivot indices from CGETRF; for 1<=i<=N, row i of the
c$$$*          matrix was interchanged with row IPIV(i).
c$$$*
c$$$*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
c$$$*          On entry, the right hand side matrix B.
c$$$*          On exit, the solution matrix X.
c$$$*
c$$$*  LDB     (input) INTEGER
c$$$*          The leading dimension of the array B.  LDB >= max(1,N).
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      COMPLEX            ONE
c$$$      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      LOGICAL            NOTRAN
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      LOGICAL            LSAME
c$$$      EXTERNAL           LSAME
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CLASWP, CTRSM, XERBLA
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      INFO = 0
c$$$      NOTRAN = LSAME( TRANS, 'N' )
c$$$      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
c$$$     $    LSAME( TRANS, 'C' ) ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( NRHS.LT.0 ) THEN
c$$$         INFO = -3
c$$$      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -5
c$$$      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -8
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGETRS', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      IF( N.EQ.0 .OR. NRHS.EQ.0 )
c$$$     $   RETURN
c$$$*
c$$$      IF( NOTRAN ) THEN
c$$$*
c$$$*        Solve A * X = B.
c$$$*
c$$$*        Apply row interchanges to the right hand sides.
c$$$*
c$$$         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
c$$$*
c$$$*        Solve L*X = B, overwriting B with X.
c$$$*
c$$$         CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
c$$$     $               ONE, A, LDA, B, LDB )
c$$$*
c$$$*        Solve U*X = B, overwriting B with X.
c$$$*
c$$$         CALL CTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
c$$$     $               NRHS, ONE, A, LDA, B, LDB )
c$$$      ELSE
c$$$*
c$$$*        Solve A**T * X = B  or A**H * X = B.
c$$$*
c$$$*        Solve U'*X = B, overwriting B with X.
c$$$*
c$$$         CALL CTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,
c$$$     $               A, LDA, B, LDB )
c$$$*
c$$$*        Solve L'*X = B, overwriting B with X.
c$$$*
c$$$         CALL CTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,
c$$$     $               LDA, B, LDB )
c$$$*
c$$$*        Apply row interchanges to the solution vectors.
c$$$*
c$$$         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
c$$$      END IF
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of CGETRS
c$$$*
c$$$      END
      SUBROUTINE CTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX            ALPHA
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  CTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX         .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - COMPLEX          array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B
*           or    B := alpha*inv( conjg( A' ) )*B.
*
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' )
*           or    B := alpha*B*inv( conjg( A' ) ).
*
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of CTRSM .
*
      END

c$$$      SUBROUTINE CGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,
c$$$     $                  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
c$$$*
c$$$*  -- LAPACK driver routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      CHARACTER          JOBVS, SORT
c$$$      INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      LOGICAL            BWORK( * )
c$$$      REAL               RWORK( * )
c$$$      COMPLEX            A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
c$$$*     ..
c$$$*     .. Function Arguments ..
c$$$      LOGICAL            SELECT
c$$$      EXTERNAL           SELECT
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGEES computes for an N-by-N complex nonsymmetric matrix A, the
c$$$*  eigenvalues, the Schur form T, and, optionally, the matrix of Schur
c$$$*  vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
c$$$*
c$$$*  Optionally, it also orders the eigenvalues on the diagonal of the
c$$$*  Schur form so that selected eigenvalues are at the top left.
c$$$*  The leading columns of Z then form an orthonormal basis for the
c$$$*  invariant subspace corresponding to the selected eigenvalues.
c$$$
c$$$*  A complex matrix is in Schur form if it is upper triangular.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  JOBVS   (input) CHARACTER*1
c$$$*          = 'N': Schur vectors are not computed;
c$$$*          = 'V': Schur vectors are computed.
c$$$*
c$$$*  SORT    (input) CHARACTER*1
c$$$*          Specifies whether or not to order the eigenvalues on the
c$$$*          diagonal of the Schur form.
c$$$*          = 'N': Eigenvalues are not ordered:
c$$$*          = 'S': Eigenvalues are ordered (see SELECT).
c$$$*
c$$$*  SELECT  (external procedure) LOGICAL FUNCTION of one COMPLEX argument
c$$$*          SELECT must be declared EXTERNAL in the calling subroutine.
c$$$*          If SORT = 'S', SELECT is used to select eigenvalues to order
c$$$*          to the top left of the Schur form.
c$$$*          IF SORT = 'N', SELECT is not referenced.
c$$$*          The eigenvalue W(j) is selected if SELECT(W(j)) is true.
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The order of the matrix A. N >= 0.
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N)
c$$$*          On entry, the N-by-N matrix A.
c$$$*          On exit, A has been overwritten by its Schur form T.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  SDIM    (output) INTEGER
c$$$*          If SORT = 'N', SDIM = 0.
c$$$*          If SORT = 'S', SDIM = number of eigenvalues for which
c$$$*                         SELECT is true.
c$$$*
c$$$*  W       (output) COMPLEX array, dimension (N)
c$$$*          W contains the computed eigenvalues, in the same order that
c$$$*          they appear on the diagonal of the output Schur form T.
c$$$*
c$$$*  VS      (output) COMPLEX array, dimension (LDVS,N)
c$$$*          If JOBVS = 'V', VS contains the unitary matrix Z of Schur
c$$$*          vectors.
c$$$*          If JOBVS = 'N', VS is not referenced.
c$$$*
c$$$*  LDVS    (input) INTEGER
c$$$*          The leading dimension of the array VS.  LDVS >= 1; if
c$$$*          JOBVS = 'V', LDVS >= N.
c$$$*
c$$$*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))
c$$$*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
c$$$*
c$$$*  LWORK   (input) INTEGER
c$$$*          The dimension of the array WORK.  LWORK >= max(1,2*N).
c$$$*          For good performance, LWORK must generally be larger.
c$$$*
c$$$*          If LWORK = -1, then a workspace query is assumed; the routine
c$$$*          only calculates the optimal size of the WORK array, returns
c$$$*          this value as the first entry of the WORK array, and no error
c$$$*          message related to LWORK is issued by XERBLA.
c$$$*
c$$$*  RWORK   (workspace) REAL array, dimension (N)
c$$$*
c$$$*  BWORK   (workspace) LOGICAL array, dimension (N)
c$$$*          Not referenced if SORT = 'N'.
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0: successful exit
c$$$*          < 0: if INFO = -i, the i-th argument had an illegal value.
c$$$*          > 0: if INFO = i, and i is
c$$$*               <= N:  the QR algorithm failed to compute all the
c$$$*                      eigenvalues; elements 1:ILO-1 and i+1:N of W
c$$$*                      contain those eigenvalues which have converged;
c$$$*                      if JOBVS = 'V', VS contains the matrix which
c$$$*                      reduces A to its partially converged Schur form.
c$$$*               = N+1: the eigenvalues could not be reordered because
c$$$*                      some eigenvalues were too close to separate (the
c$$$*                      problem is very ill-conditioned);
c$$$*               = N+2: after reordering, roundoff changed values of
c$$$*                      some complex eigenvalues so that leading
c$$$*                      eigenvalues in the Schur form no longer satisfy
c$$$*                      SELECT = .TRUE..  This could also be caused by
c$$$*                      underflow due to scaling.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      REAL               ZERO, ONE
c$$$      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      LOGICAL            LQUERY, SCALEA, WANTST, WANTVS
c$$$      INTEGER            HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO,
c$$$     $                   ITAU, IWRK, MAXWRK, MINWRK
c$$$      REAL               ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM
c$$$*     ..
c$$$*     .. Local Arrays ..
c$$$      REAL               DUM( 1 )
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CCOPY, CGEBAK, CGEBAL, CGEHRD, CHSEQR, CLACPY,
c$$$     $                   CLASCL, CTRSEN, CUNGHR, SLABAD, XERBLA
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      LOGICAL            LSAME
c$$$      INTEGER            ILAENV
c$$$      REAL               CLANGE, SLAMCH
c$$$      EXTERNAL           LSAME, ILAENV, CLANGE, SLAMCH
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX, SQRT
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input arguments
c$$$*
c$$$      INFO = 0
c$$$      LQUERY = ( LWORK.EQ.-1 )
c$$$      WANTVS = LSAME( JOBVS, 'V' )
c$$$      WANTST = LSAME( SORT, 'S' )
c$$$      IF( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -4
c$$$      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -6
c$$$      ELSE IF( LDVS.LT.1 .OR. ( WANTVS .AND. LDVS.LT.N ) ) THEN
c$$$         INFO = -10
c$$$      END IF
c$$$*
c$$$*     Compute workspace
c$$$*      (Note: Comments in the code beginning "Workspace:" describe the
c$$$*       minimal amount of workspace needed at that point in the code,
c$$$*       as well as the preferred amount for good performance.
c$$$*       CWorkspace refers to complex workspace, and RWorkspace to real
c$$$*       workspace. NB refers to the optimal block size for the
c$$$*       immediately following subroutine, as returned by ILAENV.
c$$$*       HSWORK refers to the workspace preferred by CHSEQR, as
c$$$*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
c$$$*       the worst case.)
c$$$*
c$$$      IF( INFO.EQ.0 ) THEN
c$$$         IF( N.EQ.0 ) THEN
c$$$            MINWRK = 1
c$$$            MAXWRK = 1
c$$$         ELSE
c$$$            MAXWRK = N + N*ILAENV( 1, 'CGEHRD', ' ', N, 1, N, 0 )
c$$$            MINWRK = 2*N
c$$$*
c$$$            CALL CHSEQR( 'S', JOBVS, N, 1, N, A, LDA, W, VS, LDVS,
c$$$     $             WORK, -1, IEVAL )
c$$$            HSWORK = WORK( 1 )
c$$$*
c$$$            IF( .NOT.WANTVS ) THEN
c$$$               MAXWRK = MAX( MAXWRK, HSWORK )
c$$$            ELSE
c$$$               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'CUNGHR',
c$$$     $                       ' ', N, 1, N, -1 ) )
c$$$               MAXWRK = MAX( MAXWRK, HSWORK )
c$$$            END IF
c$$$         END IF
c$$$         WORK( 1 ) = MAXWRK
c$$$*
c$$$         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
c$$$            INFO = -12
c$$$         END IF
c$$$      END IF
c$$$*
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGEES ', -INFO )
c$$$         RETURN
c$$$      ELSE IF( LQUERY ) THEN
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      IF( N.EQ.0 ) THEN
c$$$         SDIM = 0
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Get machine constants
c$$$*
c$$$      EPS = SLAMCH( 'P' )
c$$$      SMLNUM = SLAMCH( 'S' )
c$$$      BIGNUM = ONE / SMLNUM
c$$$      CALL SLABAD( SMLNUM, BIGNUM )
c$$$      SMLNUM = SQRT( SMLNUM ) / EPS
c$$$      BIGNUM = ONE / SMLNUM
c$$$*
c$$$*     Scale A if max element outside range [SMLNUM,BIGNUM]
c$$$*
c$$$      ANRM = CLANGE( 'M', N, N, A, LDA, DUM )
c$$$      SCALEA = .FALSE.
c$$$      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
c$$$         SCALEA = .TRUE.
c$$$         CSCALE = SMLNUM
c$$$      ELSE IF( ANRM.GT.BIGNUM ) THEN
c$$$         SCALEA = .TRUE.
c$$$         CSCALE = BIGNUM
c$$$      END IF
c$$$      IF( SCALEA )
c$$$     $   CALL CLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
c$$$*
c$$$*     Permute the matrix to make it more nearly triangular
c$$$*     (CWorkspace: none)
c$$$*     (RWorkspace: need N)
c$$$*
c$$$      IBAL = 1
c$$$      CALL CGEBAL( 'P', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )
c$$$*
c$$$*     Reduce to upper Hessenberg form
c$$$*     (CWorkspace: need 2*N, prefer N+N*NB)
c$$$*     (RWorkspace: none)
c$$$*
c$$$      ITAU = 1
c$$$      IWRK = N + ITAU
c$$$      CALL CGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ),
c$$$     $             LWORK-IWRK+1, IERR )
c$$$*
c$$$      IF( WANTVS ) THEN
c$$$*
c$$$*        Copy Householder vectors to VS
c$$$*
c$$$         CALL CLACPY( 'L', N, N, A, LDA, VS, LDVS )
c$$$*
c$$$*        Generate unitary matrix in VS
c$$$*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
c$$$*        (RWorkspace: none)
c$$$*
c$$$         CALL CUNGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ),
c$$$     $                LWORK-IWRK+1, IERR )
c$$$      END IF
c$$$*
c$$$      SDIM = 0
c$$$*
c$$$*     Perform QR iteration, accumulating Schur vectors in VS if desired
c$$$*     (CWorkspace: need 1, prefer HSWORK (see comments) )
c$$$*     (RWorkspace: none)
c$$$*
c$$$      IWRK = ITAU
c$$$      CALL CHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, W, VS, LDVS,
c$$$     $             WORK( IWRK ), LWORK-IWRK+1, IEVAL )
c$$$      IF( IEVAL.GT.0 )
c$$$     $   INFO = IEVAL
c$$$*
c$$$*     Sort eigenvalues if desired
c$$$*
c$$$      IF( WANTST .AND. INFO.EQ.0 ) THEN
c$$$         IF( SCALEA )
c$$$     $      CALL CLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, W, N, IERR )
c$$$         DO 10 I = 1, N
c$$$            BWORK( I ) = SELECT( W( I ) )
c$$$   10    CONTINUE
c$$$*
c$$$*        Reorder eigenvalues and transform Schur vectors
c$$$*        (CWorkspace: none)
c$$$*        (RWorkspace: none)
c$$$*
c$$$         CALL CTRSEN( 'N', JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM,
c$$$     $                S, SEP, WORK( IWRK ), LWORK-IWRK+1, ICOND )
c$$$      END IF
c$$$*
c$$$      IF( WANTVS ) THEN
c$$$*
c$$$*        Undo balancing
c$$$*        (CWorkspace: none)
c$$$*        (RWorkspace: need N)
c$$$*
c$$$         CALL CGEBAK( 'P', 'R', N, ILO, IHI, RWORK( IBAL ), N, VS, LDVS,
c$$$     $                IERR )
c$$$      END IF
c$$$*
c$$$      IF( SCALEA ) THEN
c$$$*
c$$$*        Undo scaling for the Schur form of A
c$$$*
c$$$         CALL CLASCL( 'U', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
c$$$         CALL CCOPY( N, A, LDA+1, W, 1 )
c$$$      END IF
c$$$*
c$$$      WORK( 1 ) = MAXWRK
c$$$      RETURN
c$$$*
c$$$*     End of CGEES
c$$$*
c$$$      END
c$$$
c$$$      SUBROUTINE CGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
c$$$*
c$$$*  -- LAPACK routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGEHRD reduces a complex general matrix A to upper Hessenberg form H by
c$$$*  an unitary similarity transformation:  Q' * A * Q = H .
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The order of the matrix A.  N >= 0.
c$$$*
c$$$*  ILO     (input) INTEGER
c$$$*  IHI     (input) INTEGER
c$$$*          It is assumed that A is already upper triangular in rows
c$$$*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
c$$$*          set by a previous call to CGEBAL; otherwise they should be
c$$$*          set to 1 and N respectively. See Further Details.
c$$$*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N)
c$$$*          On entry, the N-by-N general matrix to be reduced.
c$$$*          On exit, the upper triangle and the first subdiagonal of A
c$$$*          are overwritten with the upper Hessenberg matrix H, and the
c$$$*          elements below the first subdiagonal, with the array TAU,
c$$$*          represent the unitary matrix Q as a product of elementary
c$$$*          reflectors. See Further Details.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  TAU     (output) COMPLEX array, dimension (N-1)
c$$$*          The scalar factors of the elementary reflectors (see Further
c$$$*          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
c$$$*          zero.
c$$$*
c$$$*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
c$$$*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
c$$$*
c$$$*  LWORK   (input) INTEGER
c$$$*          The length of the array WORK.  LWORK >= max(1,N).
c$$$*          For optimum performance LWORK >= N*NB, where NB is the
c$$$*          optimal blocksize.
c$$$*
c$$$*          If LWORK = -1, then a workspace query is assumed; the routine
c$$$*          only calculates the optimal size of the WORK array, returns
c$$$*          this value as the first entry of the WORK array, and no error
c$$$*          message related to LWORK is issued by XERBLA.
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value.
c$$$*
c$$$*  Further Details
c$$$*  ===============
c$$$*
c$$$*  The matrix Q is represented as a product of (ihi-ilo) elementary
c$$$*  reflectors
c$$$*
c$$$*     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
c$$$*
c$$$*  Each H(i) has the form
c$$$*
c$$$*     H(i) = I - tau * v * v'
c$$$*
c$$$*  where tau is a complex scalar, and v is a complex vector with
c$$$*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
c$$$*  exit in A(i+2:ihi,i), and tau in TAU(i).
c$$$*
c$$$*  The contents of A are illustrated by the following example, with
c$$$*  n = 7, ilo = 2 and ihi = 6:
c$$$*
c$$$*  on entry,                        on exit,
c$$$*
c$$$*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
c$$$*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
c$$$*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
c$$$*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
c$$$*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
c$$$*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
c$$$*  (                         a )    (                          a )
c$$$*
c$$$*  where a denotes an element of the original matrix A, h denotes a
c$$$*  modified element of the upper Hessenberg matrix H, and vi denotes an
c$$$*  element of the vector defining H(i).
c$$$*
c$$$*  This file is a slight modification of LAPACK-3.0's CGEHRD
c$$$*  subroutine incorporating improvements proposed by Quintana-Orti and
c$$$*  Van de Geijn (2005). 
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      INTEGER            NBMAX, LDT
c$$$      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
c$$$      COMPLEX            ZERO, ONE
c$$$      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ), 
c$$$     $                     ONE = ( 1.0E+0, 0.0E+0 ) )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      LOGICAL            LQUERY
c$$$      INTEGER            I, IB, IINFO, IWS, J, LDWORK, LWKOPT, NB,
c$$$     $                   NBMIN, NH, NX
c$$$      COMPLEX            EI
c$$$*     ..
c$$$*     .. Local Arrays ..
c$$$      COMPLEX            T( LDT, NBMAX )
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CAXPY, CGEHD2, CGEMM, CLAHR2, CLARFB, CTRMM,
c$$$     $                   XERBLA
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX, MIN
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      INTEGER            ILAENV
c$$$      EXTERNAL           ILAENV
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters
c$$$*
c$$$      INFO = 0
c$$$      NB = MIN( NBMAX, ILAENV( 1, 'CGEHRD', ' ', N, ILO, IHI, -1 ) )
c$$$      LWKOPT = N*NB
c$$$      WORK( 1 ) = LWKOPT
c$$$      LQUERY = ( LWORK.EQ.-1 )
c$$$      IF( N.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
c$$$         INFO = -3
c$$$      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -5
c$$$      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
c$$$         INFO = -8
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGEHRD', -INFO )
c$$$         RETURN
c$$$      ELSE IF( LQUERY ) THEN
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
c$$$*
c$$$      DO 10 I = 1, ILO - 1
c$$$         TAU( I ) = ZERO
c$$$   10 CONTINUE
c$$$      DO 20 I = MAX( 1, IHI ), N - 1
c$$$         TAU( I ) = ZERO
c$$$   20 CONTINUE
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      NH = IHI - ILO + 1
c$$$      IF( NH.LE.1 ) THEN
c$$$         WORK( 1 ) = 1
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Determine the block size
c$$$*
c$$$      NB = MIN( NBMAX, ILAENV( 1, 'CGEHRD', ' ', N, ILO, IHI, -1 ) )
c$$$
c$$$      NBMIN = 2
c$$$      IWS = 1
c$$$      IF( NB.GT.1 .AND. NB.LT.NH ) THEN
c$$$*
c$$$*        Determine when to cross over from blocked to unblocked code
c$$$*        (last block is always handled by unblocked code)
c$$$*
c$$$         NX = MAX( NB, ILAENV( 3, 'CGEHRD', ' ', N, ILO, IHI, -1 ) )
c$$$         IF( NX.LT.NH ) THEN
c$$$*
c$$$*           Determine if workspace is large enough for blocked code
c$$$*
c$$$            IWS = N*NB
c$$$            IF( LWORK.LT.IWS ) THEN
c$$$*
c$$$*              Not enough workspace to use optimal NB:  determine the
c$$$*              minimum value of NB, and reduce NB or force use of
c$$$*              unblocked code
c$$$*
c$$$               NBMIN = MAX( 2, ILAENV( 2, 'CGEHRD', ' ', N, ILO, IHI,
c$$$     $                 -1 ) )
c$$$               IF( LWORK.GE.N*NBMIN ) THEN
c$$$                  NB = LWORK / N
c$$$               ELSE
c$$$                  NB = 1
c$$$               END IF
c$$$            END IF
c$$$         END IF
c$$$      END IF
c$$$      LDWORK = N
c$$$*
c$$$      IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
c$$$*
c$$$*        Use unblocked code below
c$$$*
c$$$         I = ILO
c$$$*
c$$$      ELSE
c$$$*
c$$$*        Use blocked code
c$$$*
c$$$         DO 40 I = ILO, IHI - 1 - NX, NB
c$$$            IB = MIN( NB, IHI-I )
c$$$*
c$$$*           Reduce columns i:i+ib-1 to Hessenberg form, returning the
c$$$*           matrices V and T of the block reflector H = I - V*T*V'
c$$$*           which performs the reduction, and also the matrix Y = A*V*T
c$$$*
c$$$            CALL CLAHR2( IHI, I, IB, A( 1, I ), LDA, TAU( I ), T, LDT,
c$$$     $                   WORK, LDWORK )
c$$$*
c$$$*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
c$$$*           right, computing  A := A - Y * V'. V(i+ib,ib-1) must be set
c$$$*           to 1
c$$$*
c$$$            EI = A( I+IB, I+IB-1 )
c$$$            A( I+IB, I+IB-1 ) = ONE
c$$$            CALL CGEMM( 'No transpose', 'Conjugate transpose', 
c$$$     $                  IHI, IHI-I-IB+1,
c$$$     $                  IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE,
c$$$     $                  A( 1, I+IB ), LDA )
c$$$            A( I+IB, I+IB-1 ) = EI
c$$$*
c$$$*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
c$$$*           right
c$$$*
c$$$            CALL CTRMM( 'Right', 'Lower', 'Conjugate transpose',
c$$$     $                  'Unit', I, IB-1,
c$$$     $                  ONE, A( I+1, I ), LDA, WORK, LDWORK )
c$$$            DO 30 J = 0, IB-2
c$$$               CALL CAXPY( I, -ONE, WORK( LDWORK*J+1 ), 1,
c$$$     $                     A( 1, I+J+1 ), 1 )
c$$$   30       CONTINUE
c$$$*
c$$$*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
c$$$*           left
c$$$*
c$$$            CALL CLARFB( 'Left', 'Conjugate transpose', 'Forward',
c$$$     $                   'Columnwise',
c$$$     $                   IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA, T, LDT,
c$$$     $                   A( I+1, I+IB ), LDA, WORK, LDWORK )
c$$$   40    CONTINUE
c$$$      END IF
c$$$*
c$$$*     Use unblocked code to reduce the rest of the matrix
c$$$*
c$$$      CALL CGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
c$$$      WORK( 1 ) = IWS
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of CGEHRD
c$$$*
c$$$      END
c$$$      
c$$$      SUBROUTINE CLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            K, LDA, LDT, LDY, N, NB
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            A( LDA, * ), T( LDT, NB ), TAU( NB ),
c$$$     $                   Y( LDY, NB )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CLAHR2 reduces the first NB columns of A complex general n-BY-(n-k+1)
c$$$*  matrix A so that elements below the k-th subdiagonal are zero. The
c$$$*  reduction is performed by an unitary similarity transformation
c$$$*  Q' * A * Q. The routine returns the matrices V and T which determine
c$$$*  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.
c$$$*
c$$$*  This is an auxiliary routine called by CGEHRD.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The order of the matrix A.
c$$$*
c$$$*  K       (input) INTEGER
c$$$*          The offset for the reduction. Elements below the k-th
c$$$*          subdiagonal in the first NB columns are reduced to zero.
c$$$*          K < N.
c$$$*
c$$$*  NB      (input) INTEGER
c$$$*          The number of columns to be reduced.
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N-K+1)
c$$$*          On entry, the n-by-(n-k+1) general matrix A.
c$$$*          On exit, the elements on and above the k-th subdiagonal in
c$$$*          the first NB columns are overwritten with the corresponding
c$$$*          elements of the reduced matrix; the elements below the k-th
c$$$*          subdiagonal, with the array TAU, represent the matrix Q as a
c$$$*          product of elementary reflectors. The other columns of A are
c$$$*          unchanged. See Further Details.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  TAU     (output) COMPLEX array, dimension (NB)
c$$$*          The scalar factors of the elementary reflectors. See Further
c$$$*          Details.
c$$$*
c$$$*  T       (output) COMPLEX array, dimension (LDT,NB)
c$$$*          The upper triangular matrix T.
c$$$*
c$$$*  LDT     (input) INTEGER
c$$$*          The leading dimension of the array T.  LDT >= NB.
c$$$*
c$$$*  Y       (output) COMPLEX array, dimension (LDY,NB)
c$$$*          The n-by-nb matrix Y.
c$$$*
c$$$*  LDY     (input) INTEGER
c$$$*          The leading dimension of the array Y. LDY >= N.
c$$$*
c$$$*  Further Details
c$$$*  ===============
c$$$*
c$$$*  The matrix Q is represented as a product of nb elementary reflectors
c$$$*
c$$$*     Q = H(1) H(2) . . . H(nb).
c$$$*
c$$$*  Each H(i) has the form
c$$$*
c$$$*     H(i) = I - tau * v * v'
c$$$*
c$$$*  where tau is a complex scalar, and v is a complex vector with
c$$$*  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
c$$$*  A(i+k+1:n,i), and tau in TAU(i).
c$$$*
c$$$*  The elements of the vectors v together form the (n-k+1)-by-nb matrix
c$$$*  V which is needed, with T and Y, to apply the transformation to the
c$$$*  unreduced part of the matrix, using an update of the form:
c$$$*  A := (I - V*T*V') * (A - Y*V').
c$$$*
c$$$*  The contents of A on exit are illustrated by the following example
c$$$*  with n = 7, k = 3 and nb = 2:
c$$$*
c$$$*     ( a   a   a   a   a )
c$$$*     ( a   a   a   a   a )
c$$$*     ( a   a   a   a   a )
c$$$*     ( h   h   a   a   a )
c$$$*     ( v1  h   a   a   a )
c$$$*     ( v1  v2  a   a   a )
c$$$*     ( v1  v2  a   a   a )
c$$$*
c$$$*  where a denotes an element of the original matrix A, h denotes a
c$$$*  modified element of the upper Hessenberg matrix H, and vi denotes an
c$$$*  element of the vector defining H(i).
c$$$*
c$$$*  This file is a slight modification of LAPACK-3.0's CLAHRD
c$$$*  incorporating improvements proposed by Quintana-Orti and Van de
c$$$*  Gejin. Note that the entries of A(1:K,2:NB) differ from those
c$$$*  returned by the original LAPACK routine. This function is
c$$$*  not backward compatible with LAPACK3.0.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      COMPLEX            ZERO, ONE
c$$$      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ), 
c$$$     $                     ONE = ( 1.0E+0, 0.0E+0 ) )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      INTEGER            I
c$$$      COMPLEX            EI
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CAXPY, CCOPY, CGEMM, CGEMV, CLACPY,
c$$$     $                   CLARFG, CSCAL, CTRMM, CTRMV, CLACGV
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MIN
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      IF( N.LE.1 )
c$$$     $   RETURN
c$$$*
c$$$      DO 10 I = 1, NB
c$$$         IF( I.GT.1 ) THEN
c$$$*
c$$$*           Update A(K+1:N,I)
c$$$*
c$$$*           Update I-th column of A - Y * V'
c$$$*
c$$$            CALL CLACGV( I-1, A( K+I-1, 1 ), LDA ) 
c$$$            CALL CGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, Y(K+1,1), LDY,
c$$$     $                  A( K+I-1, 1 ), LDA, ONE, A( K+1, I ), 1 )
c$$$            CALL CLACGV( I-1, A( K+I-1, 1 ), LDA ) 
c$$$*
c$$$*           Apply I - V * T' * V' to this column (call it b) from the
c$$$*           left, using the last column of T as workspace
c$$$*
c$$$*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
c$$$*                    ( V2 )             ( b2 )
c$$$*
c$$$*           where V1 is unit lower triangular
c$$$*
c$$$*           w := V1' * b1
c$$$*
c$$$            CALL CCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
c$$$            CALL CTRMV( 'Lower', 'Conjugate transpose', 'UNIT', 
c$$$     $                  I-1, A( K+1, 1 ),
c$$$     $                  LDA, T( 1, NB ), 1 )
c$$$*
c$$$*           w := w + V2'*b2
c$$$*
c$$$            CALL CGEMV( 'Conjugate transpose', N-K-I+1, I-1, 
c$$$     $                  ONE, A( K+I, 1 ),
c$$$     $                  LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
c$$$*
c$$$*           w := T'*w
c$$$*
c$$$            CALL CTRMV( 'Upper', 'Conjugate transpose', 'NON-UNIT', 
c$$$     $                  I-1, T, LDT,
c$$$     $                  T( 1, NB ), 1 )
c$$$*
c$$$*           b2 := b2 - V2*w
c$$$*
c$$$            CALL CGEMV( 'NO TRANSPOSE', N-K-I+1, I-1, -ONE, 
c$$$     $                  A( K+I, 1 ),
c$$$     $                  LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
c$$$*
c$$$*           b1 := b1 - V1*w
c$$$*
c$$$            CALL CTRMV( 'Lower', 'NO TRANSPOSE', 
c$$$     $                  'UNIT', I-1,
c$$$     $                  A( K+1, 1 ), LDA, T( 1, NB ), 1 )
c$$$            CALL CAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
c$$$*
c$$$            A( K+I-1, I-1 ) = EI
c$$$         END IF
c$$$*
c$$$*        Generate the elementary reflector H(I) to annihilate
c$$$*        A(K+I+1:N,I)
c$$$*
c$$$         CALL CLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1,
c$$$     $                TAU( I ) )
c$$$         EI = A( K+I, I )
c$$$         A( K+I, I ) = ONE
c$$$*
c$$$*        Compute  Y(K+1:N,I)
c$$$*
c$$$         CALL CGEMV( 'NO TRANSPOSE', N-K, N-K-I+1, 
c$$$     $               ONE, A( K+1, I+1 ),
c$$$     $               LDA, A( K+I, I ), 1, ZERO, Y( K+1, I ), 1 )
c$$$         CALL CGEMV( 'Conjugate transpose', N-K-I+1, I-1, 
c$$$     $               ONE, A( K+I, 1 ), LDA,
c$$$     $               A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
c$$$         CALL CGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, 
c$$$     $               Y( K+1, 1 ), LDY,
c$$$     $               T( 1, I ), 1, ONE, Y( K+1, I ), 1 )
c$$$         CALL CSCAL( N-K, TAU( I ), Y( K+1, I ), 1 )
c$$$*
c$$$*        Compute T(1:I,I)
c$$$*
c$$$         CALL CSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
c$$$         CALL CTRMV( 'Upper', 'No Transpose', 'NON-UNIT', 
c$$$     $               I-1, T, LDT,
c$$$     $               T( 1, I ), 1 )
c$$$         T( I, I ) = TAU( I )
c$$$*
c$$$   10 CONTINUE
c$$$      A( K+NB, NB ) = EI
c$$$*
c$$$*     Compute Y(1:K,1:NB)
c$$$*
c$$$      CALL CLACPY( 'ALL', K, NB, A( 1, 2 ), LDA, Y, LDY )
c$$$      CALL CTRMM( 'RIGHT', 'Lower', 'NO TRANSPOSE', 
c$$$     $            'UNIT', K, NB,
c$$$     $            ONE, A( K+1, 1 ), LDA, Y, LDY )
c$$$      IF( N.GT.K+NB )
c$$$     $   CALL CGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', K, 
c$$$     $               NB, N-K-NB, ONE,
c$$$     $               A( 1, 2+NB ), LDA, A( K+1+NB, 1 ), LDA, ONE, Y,
c$$$     $               LDY )
c$$$      CALL CTRMM( 'RIGHT', 'Upper', 'NO TRANSPOSE', 
c$$$     $            'NON-UNIT', K, NB,
c$$$     $            ONE, T, LDT, Y, LDY )
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of CLAHR2
c$$$*
c$$$      END
c$$$      SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            INCX, N
c$$$      COMPLEX            ALPHA, TAU
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            X( * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CLARFG generates a complex elementary reflector H of order n, such
c$$$*  that
c$$$*
c$$$*        H' * ( alpha ) = ( beta ),   H' * H = I.
c$$$*             (   x   )   (   0  )
c$$$*
c$$$*  where alpha and beta are scalars, with beta real, and x is an
c$$$*  (n-1)-element complex vector. H is represented in the form
c$$$*
c$$$*        H = I - tau * ( 1 ) * ( 1 v' ) ,
c$$$*                      ( v )
c$$$*
c$$$*  where tau is a complex scalar and v is a complex (n-1)-element
c$$$*  vector. Note that H is not hermitian.
c$$$*
c$$$*  If the elements of x are all zero and alpha is real, then tau = 0
c$$$*  and H is taken to be the unit matrix.
c$$$*
c$$$*  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The order of the elementary reflector.
c$$$*
c$$$*  ALPHA   (input/output) COMPLEX
c$$$*          On entry, the value alpha.
c$$$*          On exit, it is overwritten with the value beta.
c$$$*
c$$$*  X       (input/output) COMPLEX array, dimension
c$$$*                         (1+(N-2)*abs(INCX))
c$$$*          On entry, the vector x.
c$$$*          On exit, it is overwritten with the vector v.
c$$$*
c$$$*  INCX    (input) INTEGER
c$$$*          The increment between elements of X. INCX > 0.
c$$$*
c$$$*  TAU     (output) COMPLEX
c$$$*          The value tau.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      REAL               ONE, ZERO
c$$$      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      INTEGER            J, KNT
c$$$      REAL               ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      REAL               SCNRM2, SLAMCH, SLAPY3
c$$$      COMPLEX            CLADIV
c$$$      EXTERNAL           SCNRM2, SLAMCH, SLAPY3, CLADIV
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          ABS, AIMAG, CMPLX, REAL, SIGN
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CSCAL, CSSCAL
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$
c$$$      IF( N.LE.0 ) THEN
c$$$         TAU = ZERO
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$      XNORM = SCNRM2( N-1, X, INCX )
c$$$      ALPHR = REAL( ALPHA )
c$$$      ALPHI = AIMAG( ALPHA )
c$$$*
c$$$      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
c$$$*
c$$$*        H  =  I
c$$$*
c$$$         TAU = ZERO
c$$$      ELSE
c$$$*
c$$$*        general case
c$$$*
c$$$         BETA = -SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
c$$$         SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' )
c$$$         RSAFMN = ONE / SAFMIN
c$$$*
c$$$         IF( ABS( BETA ).LT.SAFMIN ) THEN
c$$$*
c$$$*           XNORM, BETA may be inaccurate; scale X and recompute them
c$$$*
c$$$            KNT = 0
c$$$   10       CONTINUE
c$$$            KNT = KNT + 1
c$$$            CALL CSSCAL( N-1, RSAFMN, X, INCX )
c$$$            BETA = BETA*RSAFMN
c$$$            ALPHI = ALPHI*RSAFMN
c$$$            ALPHR = ALPHR*RSAFMN
c$$$            IF( ABS( BETA ).LT.SAFMIN )
c$$$     $         GO TO 10
c$$$*
c$$$*           New BETA is at most 1, at least SAFMIN
c$$$*
c$$$            XNORM = SCNRM2( N-1, X, INCX )
c$$$            ALPHA = CMPLX( ALPHR, ALPHI )
c$$$            BETA = -SIGN( SLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
c$$$            TAU = CMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
c$$$            ALPHA = CLADIV( CMPLX( ONE ), ALPHA-BETA )
c$$$            CALL CSCAL( N-1, ALPHA, X, INCX )
c$$$*
c$$$*           If ALPHA is subnormal, it may lose relative accuracy
c$$$*
c$$$            ALPHA = BETA
c$$$            DO 20 J = 1, KNT
c$$$               ALPHA = ALPHA*SAFMIN
c$$$   20       CONTINUE
c$$$         ELSE
c$$$            TAU = CMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
c$$$            ALPHA = CLADIV( CMPLX( ONE ), ALPHA-BETA )
c$$$            CALL CSCAL( N-1, ALPHA, X, INCX )
c$$$            ALPHA = BETA
c$$$         END IF
c$$$      END IF
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of CLARFG
c$$$*
c$$$      END
c$$$      COMPLEX FUNCTION CLADIV( X, Y )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      COMPLEX            X, Y
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CLADIV := X / Y, where X and Y are complex.  The computation of X / Y
c$$$*  will not overflow on an intermediary step unless the results
c$$$*  overflows.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  X       (input) COMPLEX
c$$$*  Y       (input) COMPLEX
c$$$*          The complex scalars X and Y.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Local Scalars ..
c$$$      REAL               ZI, ZR
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           SLADIV
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          AIMAG, CMPLX, REAL
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$      CALL SLADIV( REAL( X ), AIMAG( X ), REAL( Y ), AIMAG( Y ), ZR,
c$$$     $             ZI )
c$$$      CLADIV = CMPLX( ZR, ZI )
c$$$
c$$$      RETURN
c$$$*
c$$$*     End of CLADIV
c$$$*
c$$$      END
c$$$      SUBROUTINE CHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,
c$$$     $                   WORK, LWORK, INFO )
c$$$*
c$$$*  -- LAPACK driver routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
c$$$      CHARACTER          COMPZ, JOB
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
c$$$*     ..
c$$$*     Purpose
c$$$*     =======
c$$$*
c$$$*     CHSEQR computes the eigenvalues of a Hessenberg matrix H
c$$$*     and, optionally, the matrices T and Z from the Schur decomposition
c$$$*     H = Z T Z**H, where T is an upper triangular matrix (the
c$$$*     Schur form), and Z is the unitary matrix of Schur vectors.
c$$$*
c$$$*     Optionally Z may be postmultiplied into an input unitary
c$$$*     matrix Q so that this routine can give the Schur factorization
c$$$*     of a matrix A which has been reduced to the Hessenberg form H
c$$$*     by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.
c$$$*
c$$$*     Arguments
c$$$*     =========
c$$$*
c$$$*     JOB   (input) CHARACTER*1
c$$$*           = 'E':  compute eigenvalues only;
c$$$*           = 'S':  compute eigenvalues and the Schur form T.
c$$$*
c$$$*     COMPZ (input) CHARACTER*1
c$$$*           = 'N':  no Schur vectors are computed;
c$$$*           = 'I':  Z is initialized to the unit matrix and the matrix Z
c$$$*                   of Schur vectors of H is returned;
c$$$*           = 'V':  Z must contain an unitary matrix Q on entry, and
c$$$*                   the product Q*Z is returned.
c$$$*
c$$$*     N     (input) INTEGER
c$$$*           The order of the matrix H.  N .GE. 0.
c$$$*
c$$$*     ILO   (input) INTEGER
c$$$*     IHI   (input) INTEGER
c$$$*           It is assumed that H is already upper triangular in rows
c$$$*           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
c$$$*           set by a previous call to CGEBAL, and then passed to CGEHRD
c$$$*           when the matrix output by CGEBAL is reduced to Hessenberg
c$$$*           form. Otherwise ILO and IHI should be set to 1 and N
c$$$*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
c$$$*           If N = 0, then ILO = 1 and IHI = 0.
c$$$*
c$$$*     H     (input/output) COMPLEX array, dimension (LDH,N)
c$$$*           On entry, the upper Hessenberg matrix H.
c$$$*           On exit, if INFO = 0 and JOB = 'S', H contains the upper
c$$$*           triangular matrix T from the Schur decomposition (the
c$$$*           Schur form). If INFO = 0 and JOB = 'E', the contents of
c$$$*           H are unspecified on exit.  (The output value of H when
c$$$*           INFO.GT.0 is given under the description of INFO below.)
c$$$*
c$$$*           Unlike earlier versions of CHSEQR, this subroutine may
c$$$*           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
c$$$*           or j = IHI+1, IHI+2, ... N.
c$$$*
c$$$*     LDH   (input) INTEGER
c$$$*           The leading dimension of the array H. LDH .GE. max(1,N).
c$$$*
c$$$*     W        (output) COMPLEX array, dimension (N)
c$$$*           The computed eigenvalues. If JOB = 'S', the eigenvalues are
c$$$*           stored in the same order as on the diagonal of the Schur
c$$$*           form returned in H, with W(i) = H(i,i).
c$$$*
c$$$*     Z     (input/output) COMPLEX array, dimension (LDZ,N)
c$$$*           If COMPZ = 'N', Z is not referenced.
c$$$*           If COMPZ = 'I', on entry Z need not be set and on exit,
c$$$*           if INFO = 0, Z contains the unitary matrix Z of the Schur
c$$$*           vectors of H.  If COMPZ = 'V', on entry Z must contain an
c$$$*           N-by-N matrix Q, which is assumed to be equal to the unit
c$$$*           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
c$$$*           if INFO = 0, Z contains Q*Z.
c$$$*           Normally Q is the unitary matrix generated by CUNGHR
c$$$*           after the call to CGEHRD which formed the Hessenberg matrix
c$$$*           H. (The output value of Z when INFO.GT.0 is given under
c$$$*           the description of INFO below.)
c$$$*
c$$$*     LDZ   (input) INTEGER
c$$$*           The leading dimension of the array Z.  if COMPZ = 'I' or
c$$$*           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.
c$$$*
c$$$*     WORK  (workspace/output) COMPLEX array, dimension (LWORK)
c$$$*           On exit, if INFO = 0, WORK(1) returns an estimate of
c$$$*           the optimal value for LWORK.
c$$$*
c$$$*     LWORK (input) INTEGER
c$$$*           The dimension of the array WORK.  LWORK .GE. max(1,N)
c$$$*           is sufficient, but LWORK typically as large as 6*N may
c$$$*           be required for optimal performance.  A workspace query
c$$$*           to determine the optimal workspace size is recommended.
c$$$*
c$$$*           If LWORK = -1, then CHSEQR does a workspace query.
c$$$*           In this case, CHSEQR checks the input parameters and
c$$$*           estimates the optimal workspace size for the given
c$$$*           values of N, ILO and IHI.  The estimate is returned
c$$$*           in WORK(1).  No error message related to LWORK is
c$$$*           issued by XERBLA.  Neither H nor Z are accessed.
c$$$*
c$$$*
c$$$*     INFO  (output) INTEGER
c$$$*             =  0:  successful exit
c$$$*           .LT. 0:  if INFO = -i, the i-th argument had an illegal
c$$$*                    value
c$$$*           .GT. 0:  if INFO = i, CHSEQR failed to compute all of
c$$$*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
c$$$*                and WI contain those eigenvalues which have been
c$$$*                successfully computed.  (Failures are rare.)
c$$$*
c$$$*                If INFO .GT. 0 and JOB = 'E', then on exit, the
c$$$*                remaining unconverged eigenvalues are the eigen-
c$$$*                values of the upper Hessenberg matrix rows and
c$$$*                columns ILO through INFO of the final, output
c$$$*                value of H.
c$$$*
c$$$*                If INFO .GT. 0 and JOB   = 'S', then on exit
c$$$*
c$$$*           (*)  (initial value of H)*U  = U*(final value of H)
c$$$*
c$$$*                where U is a unitary matrix.  The final
c$$$*                value of  H is upper Hessenberg and triangular in
c$$$*                rows and columns INFO+1 through IHI.
c$$$*
c$$$*                If INFO .GT. 0 and COMPZ = 'V', then on exit
c$$$*
c$$$*                  (final value of Z)  =  (initial value of Z)*U
c$$$*
c$$$*                where U is the unitary matrix in (*) (regard-
c$$$*                less of the value of JOB.)
c$$$*
c$$$*                If INFO .GT. 0 and COMPZ = 'I', then on exit
c$$$*                      (final value of Z)  = U
c$$$*                where U is the unitary matrix in (*) (regard-
c$$$*                less of the value of JOB.)
c$$$*
c$$$*                If INFO .GT. 0 and COMPZ = 'N', then Z is not
c$$$*                accessed.
c$$$*
c$$$*     ================================================================
c$$$*             Default values supplied by
c$$$*             ILAENV(ISPEC,'CHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
c$$$*             It is suggested that these defaults be adjusted in order
c$$$*             to attain best performance in each particular
c$$$*             computational environment.
c$$$*
c$$$*            ISPEC=1:  The CLAHQR vs CLAQR0 crossover point.
c$$$*                      Default: 75. (Must be at least 11.)
c$$$*
c$$$*            ISPEC=2:  Recommended deflation window size.
c$$$*                      This depends on ILO, IHI and NS.  NS is the
c$$$*                      number of simultaneous shifts returned
c$$$*                      by ILAENV(ISPEC=4).  (See ISPEC=4 below.)
c$$$*                      The default for (IHI-ILO+1).LE.500 is NS.
c$$$*                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.
c$$$*
c$$$*            ISPEC=3:  Nibble crossover point. (See ILAENV for
c$$$*                      details.)  Default: 14% of deflation window
c$$$*                      size.
c$$$*
c$$$*            ISPEC=4:  Number of simultaneous shifts, NS, in
c$$$*                      a multi-shift QR iteration.
c$$$*
c$$$*                      If IHI-ILO+1 is ...
c$$$*
c$$$*                      greater than      ...but less    ... the
c$$$*                      or equal to ...      than        default is
c$$$*
c$$$*                           1               30          NS -   2(+)
c$$$*                          30               60          NS -   4(+)
c$$$*                          60              150          NS =  10(+)
c$$$*                         150              590          NS =  **
c$$$*                         590             3000          NS =  64
c$$$*                        3000             6000          NS = 128
c$$$*                        6000             infinity      NS = 256
c$$$*
c$$$*                  (+)  By default some or all matrices of this order 
c$$$*                       are passed to the implicit double shift routine
c$$$*                       CLAHQR and NS is ignored.  See ISPEC=1 above 
c$$$*                       and comments in IPARM for details.
c$$$*
c$$$*                       The asterisks (**) indicate an ad-hoc
c$$$*                       function of N increasing from 10 to 64.
c$$$*
c$$$*            ISPEC=5:  Select structured matrix multiply.
c$$$*                      (See ILAENV for details.) Default: 3.
c$$$*
c$$$*     ================================================================
c$$$*     Based on contributions by
c$$$*        Karen Braman and Ralph Byers, Department of Mathematics,
c$$$*        University of Kansas, USA
c$$$*
c$$$*     ================================================================
c$$$*     References:
c$$$*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
c$$$*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
c$$$*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
c$$$*       929--947, 2002.
c$$$*
c$$$*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
c$$$*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
c$$$*       of Matrix Analysis, volume 23, pages 948--973, 2002.
c$$$*
c$$$*     ================================================================
c$$$*     .. Parameters ..
c$$$*
c$$$*     ==== Matrices of order NTINY or smaller must be processed by
c$$$*     .    CLAHQR because of insufficient subdiagonal scratch space.
c$$$*     .    (This is a hard limit.) ====
c$$$*
c$$$*     ==== NL allocates some local workspace to help small matrices
c$$$*     .    through a rare CLAHQR failure.  NL .GT. NTINY = 11 is
c$$$*     .    required and NL .LE. NMIN = ILAENV(ISPEC=1,...) is recom-
c$$$*     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
c$$$*     .    allows up to six simultaneous shifts and a 16-by-16
c$$$*     .    deflation window.  ====
c$$$*
c$$$      INTEGER            NTINY
c$$$      PARAMETER          ( NTINY = 11 )
c$$$      INTEGER            NL
c$$$      PARAMETER          ( NL = 49 )
c$$$      COMPLEX            ZERO, ONE
c$$$      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),
c$$$     $                   ONE = ( 1.0e0, 0.0e0 ) )
c$$$      REAL               RZERO
c$$$      PARAMETER          ( RZERO = 0.0e0 )
c$$$*     ..
c$$$*     .. Local Arrays ..
c$$$      COMPLEX            HL( NL, NL ), WORKL( NL )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      INTEGER            KBOT, NMIN
c$$$      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      INTEGER            ILAENV
c$$$      LOGICAL            LSAME
c$$$      EXTERNAL           ILAENV, LSAME
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CCOPY, CLACPY, CLAHQR, CLAQR0, CLASET, XERBLA
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          CMPLX, MAX, MIN, REAL
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     ==== Decode and check the input parameters. ====
c$$$*
c$$$      WANTT = LSAME( JOB, 'S' )
c$$$      INITZ = LSAME( COMPZ, 'I' )
c$$$      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
c$$$      WORK( 1 ) = CMPLX( REAL( MAX( 1, N ) ), RZERO )
c$$$      LQUERY = LWORK.EQ.-1
c$$$*
c$$$      INFO = 0
c$$$      IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -3
c$$$      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
c$$$         INFO = -4
c$$$      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
c$$$         INFO = -5
c$$$      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -7
c$$$      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
c$$$         INFO = -10
c$$$      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
c$$$         INFO = -12
c$$$      END IF
c$$$*
c$$$      IF( INFO.NE.0 ) THEN
c$$$*
c$$$*        ==== Quick return in case of invalid argument. ====
c$$$*
c$$$         CALL XERBLA( 'CHSEQR', -INFO )
c$$$         RETURN
c$$$*
c$$$      ELSE IF( N.EQ.0 ) THEN
c$$$*
c$$$*        ==== Quick return in case N = 0; nothing to do. ====
c$$$*
c$$$         RETURN
c$$$*
c$$$      ELSE IF( LQUERY ) THEN
c$$$*
c$$$*        ==== Quick return in case of a workspace query ====
c$$$*
c$$$         CALL CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z,
c$$$     $                LDZ, WORK, LWORK, INFO )
c$$$*        ==== Ensure reported workspace size is backward-compatible with
c$$$*        .    previous LAPACK versions. ====
c$$$         WORK( 1 ) = CMPLX( MAX( REAL( WORK( 1 ) ), REAL( MAX( 1,
c$$$     $               N ) ) ), RZERO )
c$$$         RETURN
c$$$*
c$$$      ELSE
c$$$*
c$$$*        ==== copy eigenvalues isolated by CGEBAL ====
c$$$*
c$$$         IF( ILO.GT.1 )
c$$$     $      CALL CCOPY( ILO-1, H, LDH+1, W, 1 )
c$$$         IF( IHI.LT.N )
c$$$     $      CALL CCOPY( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 )
c$$$*
c$$$*        ==== Initialize Z, if requested ====
c$$$*
c$$$         IF( INITZ )
c$$$     $      CALL CLASET( 'A', N, N, ZERO, ONE, Z, LDZ )
c$$$*
c$$$*        ==== Quick return if possible ====
c$$$*
c$$$         IF( ILO.EQ.IHI ) THEN
c$$$            W( ILO ) = H( ILO, ILO )
c$$$            RETURN
c$$$         END IF
c$$$*
c$$$*        ==== CLAHQR/CLAQR0 crossover point ====
c$$$*
c$$$         NMIN = ILAENV( 1, 'CHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO,
c$$$     $          IHI, LWORK )
c$$$         NMIN = MAX( NTINY, NMIN )
c$$$*
c$$$*        ==== CLAQR0 for big matrices; CLAHQR for small ones ====
c$$$*
c$$$         IF( N.GT.NMIN ) THEN
c$$$            CALL CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI,
c$$$     $                   Z, LDZ, WORK, LWORK, INFO )
c$$$         ELSE
c$$$*
c$$$*           ==== Small matrix ====
c$$$*
c$$$            CALL CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI,
c$$$     $                   Z, LDZ, INFO )
c$$$*
c$$$            IF( INFO.GT.0 ) THEN
c$$$*
c$$$*              ==== A rare CLAHQR failure!  CLAQR0 sometimes succeeds
c$$$*              .    when CLAHQR fails. ====
c$$$*
c$$$               KBOT = INFO
c$$$*
c$$$               IF( N.GE.NL ) THEN
c$$$*
c$$$*                 ==== Larger matrices have enough subdiagonal scratch
c$$$*                 .    space to call CLAQR0 directly. ====
c$$$*
c$$$                  CALL CLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, W,
c$$$     $                         ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
c$$$*
c$$$               ELSE
c$$$*
c$$$*                 ==== Tiny matrices don't have enough subdiagonal
c$$$*                 .    scratch space to benefit from CLAQR0.  Hence,
c$$$*                 .    tiny matrices must be copied into a larger
c$$$*                 .    array before calling CLAQR0. ====
c$$$*
c$$$                  CALL CLACPY( 'A', N, N, H, LDH, HL, NL )
c$$$                  HL( N+1, N ) = ZERO
c$$$                  CALL CLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ),
c$$$     $                         NL )
c$$$                  CALL CLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W,
c$$$     $                         ILO, IHI, Z, LDZ, WORKL, NL, INFO )
c$$$                  IF( WANTT .OR. INFO.NE.0 )
c$$$     $               CALL CLACPY( 'A', N, N, HL, NL, H, LDH )
c$$$               END IF
c$$$            END IF
c$$$         END IF
c$$$*
c$$$*        ==== Clear out the trash, if necessary. ====
c$$$*
c$$$         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 )
c$$$     $      CALL CLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )
c$$$*
c$$$*        ==== Ensure reported workspace size is backward-compatible with
c$$$*        .    previous LAPACK versions. ====
c$$$*
c$$$         WORK( 1 ) = CMPLX( MAX( REAL( MAX( 1, N ) ),
c$$$     $               REAL( WORK( 1 ) ) ), RZERO )
c$$$      END IF
c$$$*
c$$$*     ==== End of CHSEQR ====
c$$$*
c$$$      END
c$$$      SUBROUTINE CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
c$$$     $                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
c$$$      LOGICAL            WANTT, WANTZ
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
c$$$*     ..
c$$$*
c$$$*     Purpose
c$$$*     =======
c$$$*
c$$$*     CLAQR0 computes the eigenvalues of a Hessenberg matrix H
c$$$*     and, optionally, the matrices T and Z from the Schur decomposition
c$$$*     H = Z T Z**H, where T is an upper triangular matrix (the
c$$$*     Schur form), and Z is the unitary matrix of Schur vectors.
c$$$*
c$$$*     Optionally Z may be postmultiplied into an input unitary
c$$$*     matrix Q so that this routine can give the Schur factorization
c$$$*     of a matrix A which has been reduced to the Hessenberg form H
c$$$*     by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.
c$$$*
c$$$*     Arguments
c$$$*     =========
c$$$*
c$$$*     WANTT   (input) LOGICAL
c$$$*          = .TRUE. : the full Schur form T is required;
c$$$*          = .FALSE.: only eigenvalues are required.
c$$$*
c$$$*     WANTZ   (input) LOGICAL
c$$$*          = .TRUE. : the matrix of Schur vectors Z is required;
c$$$*          = .FALSE.: Schur vectors are not required.
c$$$*
c$$$*     N     (input) INTEGER
c$$$*           The order of the matrix H.  N .GE. 0.
c$$$*
c$$$*     ILO   (input) INTEGER
c$$$*     IHI   (input) INTEGER
c$$$*           It is assumed that H is already upper triangular in rows
c$$$*           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
c$$$*           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
c$$$*           previous call to CGEBAL, and then passed to CGEHRD when the
c$$$*           matrix output by CGEBAL is reduced to Hessenberg form.
c$$$*           Otherwise, ILO and IHI should be set to 1 and N,
c$$$*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
c$$$*           If N = 0, then ILO = 1 and IHI = 0.
c$$$*
c$$$*     H     (input/output) COMPLEX array, dimension (LDH,N)
c$$$*           On entry, the upper Hessenberg matrix H.
c$$$*           On exit, if INFO = 0 and WANTT is .TRUE., then H
c$$$*           contains the upper triangular matrix T from the Schur
c$$$*           decomposition (the Schur form). If INFO = 0 and WANT is
c$$$*           .FALSE., then the contents of H are unspecified on exit.
c$$$*           (The output value of H when INFO.GT.0 is given under the
c$$$*           description of INFO below.)
c$$$*
c$$$*           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
c$$$*           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
c$$$*
c$$$*     LDH   (input) INTEGER
c$$$*           The leading dimension of the array H. LDH .GE. max(1,N).
c$$$*
c$$$*     W        (output) COMPLEX array, dimension (N)
c$$$*           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored
c$$$*           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are
c$$$*           stored in the same order as on the diagonal of the Schur
c$$$*           form returned in H, with W(i) = H(i,i).
c$$$*
c$$$*     Z     (input/output) COMPLEX array, dimension (LDZ,IHI)
c$$$*           If WANTZ is .FALSE., then Z is not referenced.
c$$$*           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
c$$$*           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
c$$$*           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
c$$$*           (The output value of Z when INFO.GT.0 is given under
c$$$*           the description of INFO below.)
c$$$*
c$$$*     LDZ   (input) INTEGER
c$$$*           The leading dimension of the array Z.  if WANTZ is .TRUE.
c$$$*           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
c$$$*
c$$$*     WORK  (workspace/output) COMPLEX array, dimension LWORK
c$$$*           On exit, if LWORK = -1, WORK(1) returns an estimate of
c$$$*           the optimal value for LWORK.
c$$$*
c$$$*     LWORK (input) INTEGER
c$$$*           The dimension of the array WORK.  LWORK .GE. max(1,N)
c$$$*           is sufficient, but LWORK typically as large as 6*N may
c$$$*           be required for optimal performance.  A workspace query
c$$$*           to determine the optimal workspace size is recommended.
c$$$*
c$$$*           If LWORK = -1, then CLAQR0 does a workspace query.
c$$$*           In this case, CLAQR0 checks the input parameters and
c$$$*           estimates the optimal workspace size for the given
c$$$*           values of N, ILO and IHI.  The estimate is returned
c$$$*           in WORK(1).  No error message related to LWORK is
c$$$*           issued by XERBLA.  Neither H nor Z are accessed.
c$$$*
c$$$*
c$$$*     INFO  (output) INTEGER
c$$$*             =  0:  successful exit
c$$$*           .GT. 0:  if INFO = i, CLAQR0 failed to compute all of
c$$$*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
c$$$*                and WI contain those eigenvalues which have been
c$$$*                successfully computed.  (Failures are rare.)
c$$$*
c$$$*                If INFO .GT. 0 and WANT is .FALSE., then on exit,
c$$$*                the remaining unconverged eigenvalues are the eigen-
c$$$*                values of the upper Hessenberg matrix rows and
c$$$*                columns ILO through INFO of the final, output
c$$$*                value of H.
c$$$*
c$$$*                If INFO .GT. 0 and WANTT is .TRUE., then on exit
c$$$*
c$$$*           (*)  (initial value of H)*U  = U*(final value of H)
c$$$*
c$$$*                where U is a unitary matrix.  The final
c$$$*                value of  H is upper Hessenberg and triangular in
c$$$*                rows and columns INFO+1 through IHI.
c$$$*
c$$$*                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
c$$$*
c$$$*                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
c$$$*                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
c$$$*
c$$$*                where U is the unitary matrix in (*) (regard-
c$$$*                less of the value of WANTT.)
c$$$*
c$$$*                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
c$$$*                accessed.
c$$$*
c$$$*     ================================================================
c$$$*     Based on contributions by
c$$$*        Karen Braman and Ralph Byers, Department of Mathematics,
c$$$*        University of Kansas, USA
c$$$*
c$$$*     ================================================================
c$$$*     References:
c$$$*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
c$$$*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
c$$$*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
c$$$*       929--947, 2002.
c$$$*
c$$$*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
c$$$*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
c$$$*       of Matrix Analysis, volume 23, pages 948--973, 2002.
c$$$*
c$$$*     ================================================================
c$$$*     .. Parameters ..
c$$$*
c$$$*     ==== Matrices of order NTINY or smaller must be processed by
c$$$*     .    CLAHQR because of insufficient subdiagonal scratch space.
c$$$*     .    (This is a hard limit.) ====
c$$$*
c$$$*     ==== Exceptional deflation windows:  try to cure rare
c$$$*     .    slow convergence by increasing the size of the
c$$$*     .    deflation window after KEXNW iterations. =====
c$$$*
c$$$*     ==== Exceptional shifts: try to cure rare slow convergence
c$$$*     .    with ad-hoc exceptional shifts every KEXSH iterations.
c$$$*     .    The constants WILK1 and WILK2 are used to form the
c$$$*     .    exceptional shifts. ====
c$$$*
c$$$      INTEGER            NTINY
c$$$      PARAMETER          ( NTINY = 11 )
c$$$      INTEGER            KEXNW, KEXSH
c$$$      PARAMETER          ( KEXNW = 5, KEXSH = 6 )
c$$$      REAL               WILK1
c$$$      PARAMETER          ( WILK1 = 0.75e0 )
c$$$      COMPLEX            ZERO, ONE
c$$$      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),
c$$$     $                   ONE = ( 1.0e0, 0.0e0 ) )
c$$$      REAL               TWO
c$$$      PARAMETER          ( TWO = 2.0e0 )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      COMPLEX            AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
c$$$      REAL               S
c$$$      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
c$$$     $                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS,
c$$$     $                   LWKOPT, NDFL, NH, NHO, NIBBLE, NMIN, NS, NSMAX,
c$$$     $                   NSR, NVE, NW, NWMAX, NWR
c$$$      LOGICAL            NWINC, SORTED
c$$$      CHARACTER          JBCMPZ*2
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      INTEGER            ILAENV
c$$$      EXTERNAL           ILAENV
c$$$*     ..
c$$$*     .. Local Arrays ..
c$$$      COMPLEX            ZDUM( 1, 1 )
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CLACPY, CLAHQR, CLAQR3, CLAQR4, CLAQR5
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          ABS, AIMAG, CMPLX, INT, MAX, MIN, MOD, REAL,
c$$$     $                   SQRT
c$$$*     ..
c$$$*     .. Statement Functions ..
c$$$      REAL               CABS1
c$$$*     ..
c$$$*     .. Statement Function definitions ..
c$$$      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$      INFO = 0
c$$$*
c$$$*     ==== Quick return for N = 0: nothing to do. ====
c$$$*
c$$$      IF( N.EQ.0 ) THEN
c$$$         WORK( 1 ) = ONE
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     ==== Set up job flags for ILAENV. ====
c$$$*
c$$$      IF( WANTT ) THEN
c$$$         JBCMPZ( 1: 1 ) = 'S'
c$$$      ELSE
c$$$         JBCMPZ( 1: 1 ) = 'E'
c$$$      END IF
c$$$      IF( WANTZ ) THEN
c$$$         JBCMPZ( 2: 2 ) = 'V'
c$$$      ELSE
c$$$         JBCMPZ( 2: 2 ) = 'N'
c$$$      END IF
c$$$*
c$$$*     ==== Tiny matrices must use CLAHQR. ====
c$$$*
c$$$      IF( N.LE.NTINY ) THEN
c$$$*
c$$$*        ==== Estimate optimal workspace. ====
c$$$*
c$$$         LWKOPT = 1
c$$$         IF( LWORK.NE.-1 )
c$$$     $      CALL CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
c$$$     $                   IHIZ, Z, LDZ, INFO )
c$$$      ELSE
c$$$*
c$$$*        ==== Use small bulge multi-shift QR with aggressive early
c$$$*        .    deflation on larger-than-tiny matrices. ====
c$$$*
c$$$*        ==== Hope for the best. ====
c$$$*
c$$$         INFO = 0
c$$$*
c$$$*        ==== NWR = recommended deflation window size.  At this
c$$$*        .    point,  N .GT. NTINY = 11, so there is enough
c$$$*        .    subdiagonal workspace for NWR.GE.2 as required.
c$$$*        .    (In fact, there is enough subdiagonal space for
c$$$*        .    NWR.GE.3.) ====
c$$$*
c$$$         NWR = ILAENV( 13, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         NWR = MAX( 2, NWR )
c$$$         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
c$$$         NW = NWR
c$$$*
c$$$*        ==== NSR = recommended number of simultaneous shifts.
c$$$*        .    At this point N .GT. NTINY = 11, so there is at
c$$$*        .    enough subdiagonal workspace for NSR to be even
c$$$*        .    and greater than or equal to two as required. ====
c$$$*
c$$$         NSR = ILAENV( 15, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
c$$$         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
c$$$*
c$$$*        ==== Estimate optimal workspace ====
c$$$*
c$$$*        ==== Workspace query call to CLAQR3 ====
c$$$*
c$$$         CALL CLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ,
c$$$     $                IHIZ, Z, LDZ, LS, LD, W, H, LDH, N, H, LDH, N, H,
c$$$     $                LDH, WORK, -1 )
c$$$*
c$$$*        ==== Optimal workspace = MAX(CLAQR5, CLAQR3) ====
c$$$*
c$$$         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
c$$$*
c$$$*        ==== Quick return in case of workspace query. ====
c$$$*
c$$$         IF( LWORK.EQ.-1 ) THEN
c$$$            WORK( 1 ) = CMPLX( LWKOPT, 0 )
c$$$            RETURN
c$$$         END IF
c$$$*
c$$$*        ==== CLAHQR/CLAQR0 crossover point ====
c$$$*
c$$$         NMIN = ILAENV( 12, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         NMIN = MAX( NTINY, NMIN )
c$$$*
c$$$*        ==== Nibble crossover point ====
c$$$*
c$$$         NIBBLE = ILAENV( 14, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         NIBBLE = MAX( 0, NIBBLE )
c$$$*
c$$$*        ==== Accumulate reflections during ttswp?  Use block
c$$$*        .    2-by-2 structure during matrix-matrix multiply? ====
c$$$*
c$$$         KACC22 = ILAENV( 16, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         KACC22 = MAX( 0, KACC22 )
c$$$         KACC22 = MIN( 2, KACC22 )
c$$$*
c$$$*        ==== NWMAX = the largest possible deflation window for
c$$$*        .    which there is sufficient workspace. ====
c$$$*
c$$$         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
c$$$*
c$$$*        ==== NSMAX = the Largest number of simultaneous shifts
c$$$*        .    for which there is sufficient workspace. ====
c$$$*
c$$$         NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
c$$$         NSMAX = NSMAX - MOD( NSMAX, 2 )
c$$$*
c$$$*        ==== NDFL: an iteration count restarted at deflation. ====
c$$$*
c$$$         NDFL = 1
c$$$*
c$$$*        ==== ITMAX = iteration limit ====
c$$$*
c$$$         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
c$$$*
c$$$*        ==== Last row and column in the active block ====
c$$$*
c$$$         KBOT = IHI
c$$$*
c$$$*        ==== Main Loop ====
c$$$*
c$$$         DO 70 IT = 1, ITMAX
c$$$*
c$$$*           ==== Done when KBOT falls below ILO ====
c$$$*
c$$$            IF( KBOT.LT.ILO )
c$$$     $         GO TO 80
c$$$*
c$$$*           ==== Locate active block ====
c$$$*
c$$$            DO 10 K = KBOT, ILO + 1, -1
c$$$               IF( H( K, K-1 ).EQ.ZERO )
c$$$     $            GO TO 20
c$$$   10       CONTINUE
c$$$            K = ILO
c$$$   20       CONTINUE
c$$$            KTOP = K
c$$$*
c$$$*           ==== Select deflation window size ====
c$$$*
c$$$            NH = KBOT - KTOP + 1
c$$$            IF( NDFL.LT.KEXNW .OR. NH.LT.NW ) THEN
c$$$*
c$$$*              ==== Typical deflation window.  If possible and
c$$$*              .    advisable, nibble the entire active block.
c$$$*              .    If not, use size NWR or NWR+1 depending upon
c$$$*              .    which has the smaller corresponding subdiagonal
c$$$*              .    entry (a heuristic). ====
c$$$*
c$$$               NWINC = .TRUE.
c$$$               IF( NH.LE.MIN( NMIN, NWMAX ) ) THEN
c$$$                  NW = NH
c$$$               ELSE
c$$$                  NW = MIN( NWR, NH, NWMAX )
c$$$                  IF( NW.LT.NWMAX ) THEN
c$$$                     IF( NW.GE.NH-1 ) THEN
c$$$                        NW = NH
c$$$                     ELSE
c$$$                        KWTOP = KBOT - NW + 1
c$$$                        IF( CABS1( H( KWTOP, KWTOP-1 ) ).GT.
c$$$     $                      CABS1( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
c$$$                     END IF
c$$$                  END IF
c$$$               END IF
c$$$            ELSE
c$$$*
c$$$*              ==== Exceptional deflation window.  If there have
c$$$*              .    been no deflations in KEXNW or more iterations,
c$$$*              .    then vary the deflation window size.   At first,
c$$$*              .    because, larger windows are, in general, more
c$$$*              .    powerful than smaller ones, rapidly increase the
c$$$*              .    window up to the maximum reasonable and possible.
c$$$*              .    Then maybe try a slightly smaller window.  ====
c$$$*
c$$$               IF( NWINC .AND. NW.LT.MIN( NWMAX, NH ) ) THEN
c$$$                  NW = MIN( NWMAX, NH, 2*NW )
c$$$               ELSE
c$$$                  NWINC = .FALSE.
c$$$                  IF( NW.EQ.NH .AND. NH.GT.2 )
c$$$     $               NW = NH - 1
c$$$               END IF
c$$$            END IF
c$$$*
c$$$*           ==== Aggressive early deflation:
c$$$*           .    split workspace under the subdiagonal into
c$$$*           .      - an nw-by-nw work array V in the lower
c$$$*           .        left-hand-corner,
c$$$*           .      - an NW-by-at-least-NW-but-more-is-better
c$$$*           .        (NW-by-NHO) horizontal work array along
c$$$*           .        the bottom edge,
c$$$*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
c$$$*           .        vertical work array along the left-hand-edge.
c$$$*           .        ====
c$$$*
c$$$            KV = N - NW + 1
c$$$            KT = NW + 1
c$$$            NHO = ( N-NW-1 ) - KT + 1
c$$$            KWV = NW + 2
c$$$            NVE = ( N-NW ) - KWV + 1
c$$$*
c$$$*           ==== Aggressive early deflation ====
c$$$*
c$$$            CALL CLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
c$$$     $                   IHIZ, Z, LDZ, LS, LD, W, H( KV, 1 ), LDH, NHO,
c$$$     $                   H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, WORK,
c$$$     $                   LWORK )
c$$$*
c$$$*           ==== Adjust KBOT accounting for new deflations. ====
c$$$*
c$$$            KBOT = KBOT - LD
c$$$*
c$$$*           ==== KS points to the shifts. ====
c$$$*
c$$$            KS = KBOT - LS + 1
c$$$*
c$$$*           ==== Skip an expensive QR sweep if there is a (partly
c$$$*           .    heuristic) reason to expect that many eigenvalues
c$$$*           .    will deflate without it.  Here, the QR sweep is
c$$$*           .    skipped if many eigenvalues have just been deflated
c$$$*           .    or if the remaining active block is small.
c$$$*
c$$$            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT-
c$$$     $          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
c$$$*
c$$$*              ==== NS = nominal number of simultaneous shifts.
c$$$*              .    This may be lowered (slightly) if CLAQR3
c$$$*              .    did not provide that many shifts. ====
c$$$*
c$$$               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
c$$$               NS = NS - MOD( NS, 2 )
c$$$*
c$$$*              ==== If there have been no deflations
c$$$*              .    in a multiple of KEXSH iterations,
c$$$*              .    then try exceptional shifts.
c$$$*              .    Otherwise use shifts provided by
c$$$*              .    CLAQR3 above or from the eigenvalues
c$$$*              .    of a trailing principal submatrix. ====
c$$$*
c$$$               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
c$$$                  KS = KBOT - NS + 1
c$$$                  DO 30 I = KBOT, KS + 1, -2
c$$$                     W( I ) = H( I, I ) + WILK1*CABS1( H( I, I-1 ) )
c$$$                     W( I-1 ) = W( I )
c$$$   30             CONTINUE
c$$$               ELSE
c$$$*
c$$$*                 ==== Got NS/2 or fewer shifts? Use CLAQR4 or
c$$$*                 .    CLAHQR on a trailing principal submatrix to
c$$$*                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
c$$$*                 .    there is enough space below the subdiagonal
c$$$*                 .    to fit an NS-by-NS scratch array.) ====
c$$$*
c$$$                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
c$$$                     KS = KBOT - NS + 1
c$$$                     KT = N - NS + 1
c$$$                     CALL CLACPY( 'A', NS, NS, H( KS, KS ), LDH,
c$$$     $                            H( KT, 1 ), LDH )
c$$$                     IF( NS.GT.NMIN ) THEN
c$$$                        CALL CLAQR4( .false., .false., NS, 1, NS,
c$$$     $                               H( KT, 1 ), LDH, W( KS ), 1, 1,
c$$$     $                               ZDUM, 1, WORK, LWORK, INF )
c$$$                     ELSE
c$$$                        CALL CLAHQR( .false., .false., NS, 1, NS,
c$$$     $                               H( KT, 1 ), LDH, W( KS ), 1, 1,
c$$$     $                               ZDUM, 1, INF )
c$$$                     END IF
c$$$                     KS = KS + INF
c$$$*
c$$$*                    ==== In case of a rare QR failure use
c$$$*                    .    eigenvalues of the trailing 2-by-2
c$$$*                    .    principal submatrix.  Scale to avoid
c$$$*                    .    overflows, underflows and subnormals.
c$$$*                    .    (The scale factor S can not be zero,
c$$$*                    .    because H(KBOT,KBOT-1) is nonzero.) ====
c$$$*
c$$$                     IF( KS.GE.KBOT ) THEN
c$$$                        S = CABS1( H( KBOT-1, KBOT-1 ) ) +
c$$$     $                      CABS1( H( KBOT, KBOT-1 ) ) +
c$$$     $                      CABS1( H( KBOT-1, KBOT ) ) +
c$$$     $                      CABS1( H( KBOT, KBOT ) )
c$$$                        AA = H( KBOT-1, KBOT-1 ) / S
c$$$                        CC = H( KBOT, KBOT-1 ) / S
c$$$                        BB = H( KBOT-1, KBOT ) / S
c$$$                        DD = H( KBOT, KBOT ) / S
c$$$                        TR2 = ( AA+DD ) / TWO
c$$$                        DET = ( AA-TR2 )*( DD-TR2 ) - BB*CC
c$$$                        RTDISC = SQRT( -DET )
c$$$                        W( KBOT-1 ) = ( TR2+RTDISC )*S
c$$$                        W( KBOT ) = ( TR2-RTDISC )*S
c$$$*
c$$$                        KS = KBOT - 1
c$$$                     END IF
c$$$                  END IF
c$$$*
c$$$                  IF( KBOT-KS+1.GT.NS ) THEN
c$$$*
c$$$*                    ==== Sort the shifts (Helps a little) ====
c$$$*
c$$$                     SORTED = .false.
c$$$                     DO 50 K = KBOT, KS + 1, -1
c$$$                        IF( SORTED )
c$$$     $                     GO TO 60
c$$$                        SORTED = .true.
c$$$                        DO 40 I = KS, K - 1
c$$$                           IF( CABS1( W( I ) ).LT.CABS1( W( I+1 ) ) )
c$$$     $                          THEN
c$$$                              SORTED = .false.
c$$$                              SWAP = W( I )
c$$$                              W( I ) = W( I+1 )
c$$$                              W( I+1 ) = SWAP
c$$$                           END IF
c$$$   40                   CONTINUE
c$$$   50                CONTINUE
c$$$   60                CONTINUE
c$$$                  END IF
c$$$               END IF
c$$$*
c$$$*              ==== If there are only two shifts, then use
c$$$*              .    only one.  ====
c$$$*
c$$$               IF( KBOT-KS+1.EQ.2 ) THEN
c$$$                  IF( CABS1( W( KBOT )-H( KBOT, KBOT ) ).LT.
c$$$     $                CABS1( W( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
c$$$                     W( KBOT-1 ) = W( KBOT )
c$$$                  ELSE
c$$$                     W( KBOT ) = W( KBOT-1 )
c$$$                  END IF
c$$$               END IF
c$$$*
c$$$*              ==== Use up to NS of the the smallest magnatiude
c$$$*              .    shifts.  If there aren't NS shifts available,
c$$$*              .    then use them all, possibly dropping one to
c$$$*              .    make the number of shifts even. ====
c$$$*
c$$$               NS = MIN( NS, KBOT-KS+1 )
c$$$               NS = NS - MOD( NS, 2 )
c$$$               KS = KBOT - NS + 1
c$$$*
c$$$*              ==== Small-bulge multi-shift QR sweep:
c$$$*              .    split workspace under the subdiagonal into
c$$$*              .    - a KDU-by-KDU work array U in the lower
c$$$*              .      left-hand-corner,
c$$$*              .    - a KDU-by-at-least-KDU-but-more-is-better
c$$$*              .      (KDU-by-NHo) horizontal work array WH along
c$$$*              .      the bottom edge,
c$$$*              .    - and an at-least-KDU-but-more-is-better-by-KDU
c$$$*              .      (NVE-by-KDU) vertical work WV arrow along
c$$$*              .      the left-hand-edge. ====
c$$$*
c$$$               KDU = 3*NS - 3
c$$$               KU = N - KDU + 1
c$$$               KWH = KDU + 1
c$$$               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
c$$$               KWV = KDU + 4
c$$$               NVE = N - KDU - KWV + 1
c$$$*
c$$$*              ==== Small-bulge multi-shift QR sweep ====
c$$$*
c$$$               CALL CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS,
c$$$     $                      W( KS ), H, LDH, ILOZ, IHIZ, Z, LDZ, WORK,
c$$$     $                      3, H( KU, 1 ), LDH, NVE, H( KWV, 1 ), LDH,
c$$$     $                      NHO, H( KU, KWH ), LDH )
c$$$            END IF
c$$$*
c$$$*           ==== Note progress (or the lack of it). ====
c$$$*
c$$$            IF( LD.GT.0 ) THEN
c$$$               NDFL = 1
c$$$            ELSE
c$$$               NDFL = NDFL + 1
c$$$            END IF
c$$$*
c$$$*           ==== End of main loop ====
c$$$   70    CONTINUE
c$$$*
c$$$*        ==== Iteration limit exceeded.  Set INFO to show where
c$$$*        .    the problem occurred and exit. ====
c$$$*
c$$$         INFO = KBOT
c$$$   80    CONTINUE
c$$$      END IF
c$$$*
c$$$*     ==== Return the optimal value of LWORK. ====
c$$$*
c$$$      WORK( 1 ) = CMPLX( LWKOPT, 0 )
c$$$*
c$$$*     ==== End of CLAQR0 ====
c$$$*
c$$$      END
c$$$      SUBROUTINE SLADIV( A, B, C, D, P, Q )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      REAL               A, B, C, D, P, Q
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  SLADIV performs complex division in  real arithmetic
c$$$*
c$$$*                        a + i*b
c$$$*             p + i*q = ---------
c$$$*                        c + i*d
c$$$*
c$$$*  The algorithm is due to Robert L. Smith and can be found
c$$$*  in D. Knuth, The art of Computer Programming, Vol.2, p.195
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  A       (input) REAL
c$$$*  B       (input) REAL
c$$$*  C       (input) REAL
c$$$*  D       (input) REAL
c$$$*          The scalars a, b, c, and d in the above expression.
c$$$*
c$$$*  P       (output) REAL
c$$$*  Q       (output) REAL
c$$$*          The scalars p and q in the above expression.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Local Scalars ..
c$$$      REAL               E, F
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          ABS
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$      IF( ABS( D ).LT.ABS( C ) ) THEN
c$$$         E = D / C
c$$$         F = C + D*E
c$$$         P = ( A+B*E ) / F
c$$$         Q = ( B-A*E ) / F
c$$$      ELSE
c$$$         E = C / D
c$$$         F = D + C*E
c$$$         P = ( B+A*E ) / F
c$$$         Q = ( -A+B*E ) / F
c$$$      END IF
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of SLADIV
c$$$*
c$$$      END
c$$$      SUBROUTINE CLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
c$$$     $                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
c$$$     $                   NV, WV, LDWV, WORK, LWORK )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
c$$$     $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
c$$$      LOGICAL            WANTT, WANTZ
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
c$$$     $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
c$$$*     ..
c$$$*
c$$$*     ******************************************************************
c$$$*     Aggressive early deflation:
c$$$*
c$$$*     This subroutine accepts as input an upper Hessenberg matrix
c$$$*     H and performs an unitary similarity transformation
c$$$*     designed to detect and deflate fully converged eigenvalues from
c$$$*     a trailing principal submatrix.  On output H has been over-
c$$$*     written by a new Hessenberg matrix that is a perturbation of
c$$$*     an unitary similarity transformation of H.  It is to be
c$$$*     hoped that the final version of H has many zero subdiagonal
c$$$*     entries.
c$$$*
c$$$*     ******************************************************************
c$$$*     WANTT   (input) LOGICAL
c$$$*          If .TRUE., then the Hessenberg matrix H is fully updated
c$$$*          so that the triangular Schur factor may be
c$$$*          computed (in cooperation with the calling subroutine).
c$$$*          If .FALSE., then only enough of H is updated to preserve
c$$$*          the eigenvalues.
c$$$*
c$$$*     WANTZ   (input) LOGICAL
c$$$*          If .TRUE., then the unitary matrix Z is updated so
c$$$*          so that the unitary Schur factor may be computed
c$$$*          (in cooperation with the calling subroutine).
c$$$*          If .FALSE., then Z is not referenced.
c$$$*
c$$$*     N       (input) INTEGER
c$$$*          The order of the matrix H and (if WANTZ is .TRUE.) the
c$$$*          order of the unitary matrix Z.
c$$$*
c$$$*     KTOP    (input) INTEGER
c$$$*          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
c$$$*          KBOT and KTOP together determine an isolated block
c$$$*          along the diagonal of the Hessenberg matrix.
c$$$*
c$$$*     KBOT    (input) INTEGER
c$$$*          It is assumed without a check that either
c$$$*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
c$$$*          determine an isolated block along the diagonal of the
c$$$*          Hessenberg matrix.
c$$$*
c$$$*     NW      (input) INTEGER
c$$$*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
c$$$*
c$$$*     H       (input/output) COMPLEX array, dimension (LDH,N)
c$$$*          On input the initial N-by-N section of H stores the
c$$$*          Hessenberg matrix undergoing aggressive early deflation.
c$$$*          On output H has been transformed by a unitary
c$$$*          similarity transformation, perturbed, and the returned
c$$$*          to Hessenberg form that (it is to be hoped) has some
c$$$*          zero subdiagonal entries.
c$$$*
c$$$*     LDH     (input) integer
c$$$*          Leading dimension of H just as declared in the calling
c$$$*          subroutine.  N .LE. LDH
c$$$*
c$$$*     ILOZ    (input) INTEGER
c$$$*     IHIZ    (input) INTEGER
c$$$*          Specify the rows of Z to which transformations must be
c$$$*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
c$$$*
c$$$*     Z       (input/output) COMPLEX array, dimension (LDZ,IHI)
c$$$*          IF WANTZ is .TRUE., then on output, the unitary
c$$$*          similarity transformation mentioned above has been
c$$$*          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
c$$$*          If WANTZ is .FALSE., then Z is unreferenced.
c$$$*
c$$$*     LDZ     (input) integer
c$$$*          The leading dimension of Z just as declared in the
c$$$*          calling subroutine.  1 .LE. LDZ.
c$$$*
c$$$*     NS      (output) integer
c$$$*          The number of unconverged (ie approximate) eigenvalues
c$$$*          returned in SR and SI that may be used as shifts by the
c$$$*          calling subroutine.
c$$$*
c$$$*     ND      (output) integer
c$$$*          The number of converged eigenvalues uncovered by this
c$$$*          subroutine.
c$$$*
c$$$*     SH      (output) COMPLEX array, dimension KBOT
c$$$*          On output, approximate eigenvalues that may
c$$$*          be used for shifts are stored in SH(KBOT-ND-NS+1)
c$$$*          through SR(KBOT-ND).  Converged eigenvalues are
c$$$*          stored in SH(KBOT-ND+1) through SH(KBOT).
c$$$*
c$$$*     V       (workspace) COMPLEX array, dimension (LDV,NW)
c$$$*          An NW-by-NW work array.
c$$$*
c$$$*     LDV     (input) integer scalar
c$$$*          The leading dimension of V just as declared in the
c$$$*          calling subroutine.  NW .LE. LDV
c$$$*
c$$$*     NH      (input) integer scalar
c$$$*          The number of columns of T.  NH.GE.NW.
c$$$*
c$$$*     T       (workspace) COMPLEX array, dimension (LDT,NW)
c$$$*
c$$$*     LDT     (input) integer
c$$$*          The leading dimension of T just as declared in the
c$$$*          calling subroutine.  NW .LE. LDT
c$$$*
c$$$*     NV      (input) integer
c$$$*          The number of rows of work array WV available for
c$$$*          workspace.  NV.GE.NW.
c$$$*
c$$$*     WV      (workspace) COMPLEX array, dimension (LDWV,NW)
c$$$*
c$$$*     LDWV    (input) integer
c$$$*          The leading dimension of W just as declared in the
c$$$*          calling subroutine.  NW .LE. LDV
c$$$*
c$$$*     WORK    (workspace) COMPLEX array, dimension LWORK.
c$$$*          On exit, WORK(1) is set to an estimate of the optimal value
c$$$*          of LWORK for the given values of N, NW, KTOP and KBOT.
c$$$*
c$$$*     LWORK   (input) integer
c$$$*          The dimension of the work array WORK.  LWORK = 2*NW
c$$$*          suffices, but greater efficiency may result from larger
c$$$*          values of LWORK.
c$$$*
c$$$*          If LWORK = -1, then a workspace query is assumed; CLAQR3
c$$$*          only estimates the optimal workspace size for the given
c$$$*          values of N, NW, KTOP and KBOT.  The estimate is returned
c$$$*          in WORK(1).  No error message related to LWORK is issued
c$$$*          by XERBLA.  Neither H nor Z are accessed.
c$$$*
c$$$*     ================================================================
c$$$*     Based on contributions by
c$$$*        Karen Braman and Ralph Byers, Department of Mathematics,
c$$$*        University of Kansas, USA
c$$$*
c$$$*     ==================================================================
c$$$*     .. Parameters ..
c$$$      COMPLEX            ZERO, ONE
c$$$      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),
c$$$     $                   ONE = ( 1.0e0, 0.0e0 ) )
c$$$      REAL               RZERO, RONE
c$$$      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0 )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      COMPLEX            BETA, CDUM, S, TAU
c$$$      REAL               FOO, SAFMAX, SAFMIN, SMLNUM, ULP
c$$$      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN,
c$$$     $                   KNT, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3,
c$$$     $                   LWKOPT, NMIN
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      REAL               SLAMCH
c$$$      INTEGER            ILAENV
c$$$      EXTERNAL           SLAMCH, ILAENV
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CCOPY, CGEHRD, CGEMM, CLACPY, CLAHQR, CLAQR4,
c$$$     $                   CLARF, CLARFG, CLASET, CTREXC, CUNGHR, SLABAD
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          ABS, AIMAG, CMPLX, CONJG, INT, MAX, MIN, REAL
c$$$*     ..
c$$$*     .. Statement Functions ..
c$$$      REAL               CABS1
c$$$*     ..
c$$$*     .. Statement Function definitions ..
c$$$      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     ==== Estimate optimal workspace. ====
c$$$*
c$$$      JW = MIN( NW, KBOT-KTOP+1 )
c$$$      IF( JW.LE.2 ) THEN
c$$$         LWKOPT = 1
c$$$      ELSE
c$$$*
c$$$*        ==== Workspace query call to CGEHRD ====
c$$$*
c$$$         CALL CGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
c$$$         LWK1 = INT( WORK( 1 ) )
c$$$*
c$$$*        ==== Workspace query call to CUNGHR ====
c$$$*
c$$$         CALL CUNGHR( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
c$$$         LWK2 = INT( WORK( 1 ) )
c$$$*
c$$$*        ==== Workspace query call to CLAQR4 ====
c$$$*
c$$$         CALL CLAQR4( .true., .true., JW, 1, JW, T, LDT, SH, 1, JW, V,
c$$$     $                LDV, WORK, -1, INFQR )
c$$$         LWK3 = INT( WORK( 1 ) )
c$$$*
c$$$*        ==== Optimal workspace ====
c$$$*
c$$$         LWKOPT = MAX( JW+MAX( LWK1, LWK2 ), LWK3 )
c$$$      END IF
c$$$*
c$$$*     ==== Quick return in case of workspace query. ====
c$$$*
c$$$      IF( LWORK.EQ.-1 ) THEN
c$$$         WORK( 1 ) = CMPLX( LWKOPT, 0 )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     ==== Nothing to do ...
c$$$*     ... for an empty active block ... ====
c$$$      NS = 0
c$$$      ND = 0
c$$$      IF( KTOP.GT.KBOT )
c$$$     $   RETURN
c$$$*     ... nor for an empty deflation window. ====
c$$$      IF( NW.LT.1 )
c$$$     $   RETURN
c$$$*
c$$$*     ==== Machine constants ====
c$$$*
c$$$      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
c$$$      SAFMAX = RONE / SAFMIN
c$$$      CALL SLABAD( SAFMIN, SAFMAX )
c$$$      ULP = SLAMCH( 'PRECISION' )
c$$$      SMLNUM = SAFMIN*( REAL( N ) / ULP )
c$$$*
c$$$*     ==== Setup deflation window ====
c$$$*
c$$$      JW = MIN( NW, KBOT-KTOP+1 )
c$$$      KWTOP = KBOT - JW + 1
c$$$      IF( KWTOP.EQ.KTOP ) THEN
c$$$         S = ZERO
c$$$      ELSE
c$$$         S = H( KWTOP, KWTOP-1 )
c$$$      END IF
c$$$*
c$$$      IF( KBOT.EQ.KWTOP ) THEN
c$$$*
c$$$*        ==== 1-by-1 deflation window: not much to do ====
c$$$*
c$$$         SH( KWTOP ) = H( KWTOP, KWTOP )
c$$$         NS = 1
c$$$         ND = 0
c$$$         IF( CABS1( S ).LE.MAX( SMLNUM, ULP*CABS1( H( KWTOP,
c$$$     $       KWTOP ) ) ) ) THEN
c$$$
c$$$            NS = 0
c$$$            ND = 1
c$$$            IF( KWTOP.GT.KTOP )
c$$$     $         H( KWTOP, KWTOP-1 ) = ZERO
c$$$         END IF
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     ==== Convert to spike-triangular form.  (In case of a
c$$$*     .    rare QR failure, this routine continues to do
c$$$*     .    aggressive early deflation using that part of
c$$$*     .    the deflation window that converged using INFQR
c$$$*     .    here and there to keep track.) ====
c$$$*
c$$$      CALL CLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
c$$$      CALL CCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
c$$$*
c$$$      CALL CLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
c$$$      NMIN = ILAENV( 12, 'CLAQR3', 'SV', JW, 1, JW, LWORK )
c$$$      IF( JW.GT.NMIN ) THEN
c$$$         CALL CLAQR4( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1,
c$$$     $                JW, V, LDV, WORK, LWORK, INFQR )
c$$$      ELSE
c$$$         CALL CLAHQR( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1,
c$$$     $                JW, V, LDV, INFQR )
c$$$      END IF
c$$$*
c$$$*     ==== Deflation detection loop ====
c$$$*
c$$$      NS = JW
c$$$      ILST = INFQR + 1
c$$$      DO 10 KNT = INFQR + 1, JW
c$$$*
c$$$*        ==== Small spike tip deflation test ====
c$$$*
c$$$         FOO = CABS1( T( NS, NS ) )
c$$$         IF( FOO.EQ.RZERO )
c$$$     $      FOO = CABS1( S )
c$$$         IF( CABS1( S )*CABS1( V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) )
c$$$     $        THEN
c$$$*
c$$$*           ==== One more converged eigenvalue ====
c$$$*
c$$$            NS = NS - 1
c$$$         ELSE
c$$$*
c$$$*           ==== One undflatable eigenvalue.  Move it up out of the
c$$$*           .    way.   (CTREXC can not fail in this case.) ====
c$$$*
c$$$            IFST = NS
c$$$            CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
c$$$            ILST = ILST + 1
c$$$         END IF
c$$$   10 CONTINUE
c$$$*
c$$$*        ==== Return to Hessenberg form ====
c$$$*
c$$$      IF( NS.EQ.0 )
c$$$     $   S = ZERO
c$$$*
c$$$      IF( NS.LT.JW ) THEN
c$$$*
c$$$*        ==== sorting the diagonal of T improves accuracy for
c$$$*        .    graded matrices.  ====
c$$$*
c$$$         DO 30 I = INFQR + 1, NS
c$$$            IFST = I
c$$$            DO 20 J = I + 1, NS
c$$$               IF( CABS1( T( J, J ) ).GT.CABS1( T( IFST, IFST ) ) )
c$$$     $            IFST = J
c$$$   20       CONTINUE
c$$$            ILST = I
c$$$            IF( IFST.NE.ILST )
c$$$     $         CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
c$$$   30    CONTINUE
c$$$      END IF
c$$$*
c$$$*     ==== Restore shift/eigenvalue array from T ====
c$$$*
c$$$      DO 40 I = INFQR + 1, JW
c$$$         SH( KWTOP+I-1 ) = T( I, I )
c$$$   40 CONTINUE
c$$$*
c$$$*
c$$$      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
c$$$         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
c$$$*
c$$$*           ==== Reflect spike back into lower triangle ====
c$$$*
c$$$            CALL CCOPY( NS, V, LDV, WORK, 1 )
c$$$            DO 50 I = 1, NS
c$$$               WORK( I ) = CONJG( WORK( I ) )
c$$$   50       CONTINUE
c$$$            BETA = WORK( 1 )
c$$$            CALL CLARFG( NS, BETA, WORK( 2 ), 1, TAU )
c$$$            WORK( 1 ) = ONE
c$$$*
c$$$            CALL CLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
c$$$*
c$$$            CALL CLARF( 'L', NS, JW, WORK, 1, CONJG( TAU ), T, LDT,
c$$$     $                  WORK( JW+1 ) )
c$$$            CALL CLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT,
c$$$     $                  WORK( JW+1 ) )
c$$$            CALL CLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV,
c$$$     $                  WORK( JW+1 ) )
c$$$*
c$$$            CALL CGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),
c$$$     $                   LWORK-JW, INFO )
c$$$         END IF
c$$$*
c$$$*        ==== Copy updated reduced window into place ====
c$$$*
c$$$         IF( KWTOP.GT.1 )
c$$$     $      H( KWTOP, KWTOP-1 ) = S*CONJG( V( 1, 1 ) )
c$$$         CALL CLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
c$$$         CALL CCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ),
c$$$     $               LDH+1 )
c$$$*
c$$$*        ==== Accumulate orthogonal matrix in order update
c$$$*        .    H and Z, if requested.  (A modified version
c$$$*        .    of  CUNGHR that accumulates block Householder
c$$$*        .    transformations into V directly might be
c$$$*        .    marginally more efficient than the following.) ====
c$$$*
c$$$         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
c$$$            CALL CUNGHR( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),
c$$$     $                   LWORK-JW, INFO )
c$$$            CALL CGEMM( 'N', 'N', JW, NS, NS, ONE, V, LDV, T, LDT, ZERO,
c$$$     $                  WV, LDWV )
c$$$            CALL CLACPY( 'A', JW, NS, WV, LDWV, V, LDV )
c$$$         END IF
c$$$*
c$$$*        ==== Update vertical slab in H ====
c$$$*
c$$$         IF( WANTT ) THEN
c$$$            LTOP = 1
c$$$         ELSE
c$$$            LTOP = KTOP
c$$$         END IF
c$$$         DO 60 KROW = LTOP, KWTOP - 1, NV
c$$$            KLN = MIN( NV, KWTOP-KROW )
c$$$            CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ),
c$$$     $                  LDH, V, LDV, ZERO, WV, LDWV )
c$$$            CALL CLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
c$$$   60    CONTINUE
c$$$*
c$$$*        ==== Update horizontal slab in H ====
c$$$*
c$$$         IF( WANTT ) THEN
c$$$            DO 70 KCOL = KBOT + 1, N, NH
c$$$               KLN = MIN( NH, N-KCOL+1 )
c$$$               CALL CGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV,
c$$$     $                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
c$$$               CALL CLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ),
c$$$     $                      LDH )
c$$$   70       CONTINUE
c$$$         END IF
c$$$*
c$$$*        ==== Update vertical slab in Z ====
c$$$*
c$$$         IF( WANTZ ) THEN
c$$$            DO 80 KROW = ILOZ, IHIZ, NV
c$$$               KLN = MIN( NV, IHIZ-KROW+1 )
c$$$               CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ),
c$$$     $                     LDZ, V, LDV, ZERO, WV, LDWV )
c$$$               CALL CLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ),
c$$$     $                      LDZ )
c$$$   80       CONTINUE
c$$$         END IF
c$$$      END IF
c$$$*
c$$$*     ==== Return the number of deflations ... ====
c$$$*
c$$$      ND = JW - NS
c$$$*
c$$$*     ==== ... and the number of shifts. (Subtracting
c$$$*     .    INFQR from the spike length takes care
c$$$*     .    of the case of a rare QR failure while
c$$$*     .    calculating eigenvalues of the deflation
c$$$*     .    window.)  ====
c$$$*
c$$$      NS = NS - INFQR
c$$$*
c$$$*      ==== Return optimal workspace. ====
c$$$*
c$$$      WORK( 1 ) = CMPLX( LWKOPT, 0 )
c$$$*
c$$$*     ==== End of CLAQR3 ====
c$$$*
c$$$      END
c$$$      SUBROUTINE CLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
c$$$     $                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
c$$$      LOGICAL            WANTT, WANTZ
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
c$$$*     ..
c$$$*
c$$$*     This subroutine implements one level of recursion for CLAQR0.
c$$$*     It is a complete implementation of the small bulge multi-shift
c$$$*     QR algorithm.  It may be called by CLAQR0 and, for large enough
c$$$*     deflation window size, it may be called by CLAQR3.  This
c$$$*     subroutine is identical to CLAQR0 except that it calls CLAQR2
c$$$*     instead of CLAQR3.
c$$$*
c$$$*     Purpose
c$$$*     =======
c$$$*
c$$$*     CLAQR4 computes the eigenvalues of a Hessenberg matrix H
c$$$*     and, optionally, the matrices T and Z from the Schur decomposition
c$$$*     H = Z T Z**H, where T is an upper triangular matrix (the
c$$$*     Schur form), and Z is the unitary matrix of Schur vectors.
c$$$*
c$$$*     Optionally Z may be postmultiplied into an input unitary
c$$$*     matrix Q so that this routine can give the Schur factorization
c$$$*     of a matrix A which has been reduced to the Hessenberg form H
c$$$*     by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.
c$$$*
c$$$*     Arguments
c$$$*     =========
c$$$*
c$$$*     WANTT   (input) LOGICAL
c$$$*          = .TRUE. : the full Schur form T is required;
c$$$*          = .FALSE.: only eigenvalues are required.
c$$$*
c$$$*     WANTZ   (input) LOGICAL
c$$$*          = .TRUE. : the matrix of Schur vectors Z is required;
c$$$*          = .FALSE.: Schur vectors are not required.
c$$$*
c$$$*     N     (input) INTEGER
c$$$*           The order of the matrix H.  N .GE. 0.
c$$$*
c$$$*     ILO   (input) INTEGER
c$$$*     IHI   (input) INTEGER
c$$$*           It is assumed that H is already upper triangular in rows
c$$$*           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
c$$$*           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
c$$$*           previous call to CGEBAL, and then passed to CGEHRD when the
c$$$*           matrix output by CGEBAL is reduced to Hessenberg form.
c$$$*           Otherwise, ILO and IHI should be set to 1 and N,
c$$$*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
c$$$*           If N = 0, then ILO = 1 and IHI = 0.
c$$$*
c$$$*     H     (input/output) COMPLEX array, dimension (LDH,N)
c$$$*           On entry, the upper Hessenberg matrix H.
c$$$*           On exit, if INFO = 0 and WANTT is .TRUE., then H
c$$$*           contains the upper triangular matrix T from the Schur
c$$$*           decomposition (the Schur form). If INFO = 0 and WANT is
c$$$*           .FALSE., then the contents of H are unspecified on exit.
c$$$*           (The output value of H when INFO.GT.0 is given under the
c$$$*           description of INFO below.)
c$$$*
c$$$*           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
c$$$*           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
c$$$*
c$$$*     LDH   (input) INTEGER
c$$$*           The leading dimension of the array H. LDH .GE. max(1,N).
c$$$*
c$$$*     W        (output) COMPLEX array, dimension (N)
c$$$*           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored
c$$$*           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are
c$$$*           stored in the same order as on the diagonal of the Schur
c$$$*           form returned in H, with W(i) = H(i,i).
c$$$*
c$$$*     Z     (input/output) COMPLEX array, dimension (LDZ,IHI)
c$$$*           If WANTZ is .FALSE., then Z is not referenced.
c$$$*           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
c$$$*           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
c$$$*           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
c$$$*           (The output value of Z when INFO.GT.0 is given under
c$$$*           the description of INFO below.)
c$$$*
c$$$*     LDZ   (input) INTEGER
c$$$*           The leading dimension of the array Z.  if WANTZ is .TRUE.
c$$$*           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
c$$$*
c$$$*     WORK  (workspace/output) COMPLEX array, dimension LWORK
c$$$*           On exit, if LWORK = -1, WORK(1) returns an estimate of
c$$$*           the optimal value for LWORK.
c$$$*
c$$$*     LWORK (input) INTEGER
c$$$*           The dimension of the array WORK.  LWORK .GE. max(1,N)
c$$$*           is sufficient, but LWORK typically as large as 6*N may
c$$$*           be required for optimal performance.  A workspace query
c$$$*           to determine the optimal workspace size is recommended.
c$$$*
c$$$*           If LWORK = -1, then CLAQR4 does a workspace query.
c$$$*           In this case, CLAQR4 checks the input parameters and
c$$$*           estimates the optimal workspace size for the given
c$$$*           values of N, ILO and IHI.  The estimate is returned
c$$$*           in WORK(1).  No error message related to LWORK is
c$$$*           issued by XERBLA.  Neither H nor Z are accessed.
c$$$*
c$$$*
c$$$*     INFO  (output) INTEGER
c$$$*             =  0:  successful exit
c$$$*           .GT. 0:  if INFO = i, CLAQR4 failed to compute all of
c$$$*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
c$$$*                and WI contain those eigenvalues which have been
c$$$*                successfully computed.  (Failures are rare.)
c$$$*
c$$$*                If INFO .GT. 0 and WANT is .FALSE., then on exit,
c$$$*                the remaining unconverged eigenvalues are the eigen-
c$$$*                values of the upper Hessenberg matrix rows and
c$$$*                columns ILO through INFO of the final, output
c$$$*                value of H.
c$$$*
c$$$*                If INFO .GT. 0 and WANTT is .TRUE., then on exit
c$$$*
c$$$*           (*)  (initial value of H)*U  = U*(final value of H)
c$$$*
c$$$*                where U is a unitary matrix.  The final
c$$$*                value of  H is upper Hessenberg and triangular in
c$$$*                rows and columns INFO+1 through IHI.
c$$$*
c$$$*                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
c$$$*
c$$$*                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
c$$$*                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
c$$$*
c$$$*                where U is the unitary matrix in (*) (regard-
c$$$*                less of the value of WANTT.)
c$$$*
c$$$*                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
c$$$*                accessed.
c$$$*
c$$$*     ================================================================
c$$$*     Based on contributions by
c$$$*        Karen Braman and Ralph Byers, Department of Mathematics,
c$$$*        University of Kansas, USA
c$$$*
c$$$*     ================================================================
c$$$*     References:
c$$$*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
c$$$*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
c$$$*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
c$$$*       929--947, 2002.
c$$$*
c$$$*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
c$$$*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
c$$$*       of Matrix Analysis, volume 23, pages 948--973, 2002.
c$$$*
c$$$*     ================================================================
c$$$*     .. Parameters ..
c$$$*
c$$$*     ==== Matrices of order NTINY or smaller must be processed by
c$$$*     .    CLAHQR because of insufficient subdiagonal scratch space.
c$$$*     .    (This is a hard limit.) ====
c$$$*
c$$$*     ==== Exceptional deflation windows:  try to cure rare
c$$$*     .    slow convergence by increasing the size of the
c$$$*     .    deflation window after KEXNW iterations. =====
c$$$*
c$$$*     ==== Exceptional shifts: try to cure rare slow convergence
c$$$*     .    with ad-hoc exceptional shifts every KEXSH iterations.
c$$$*     .    The constants WILK1 and WILK2 are used to form the
c$$$*     .    exceptional shifts. ====
c$$$*
c$$$      INTEGER            NTINY
c$$$      PARAMETER          ( NTINY = 11 )
c$$$      INTEGER            KEXNW, KEXSH
c$$$      PARAMETER          ( KEXNW = 5, KEXSH = 6 )
c$$$      REAL               WILK1
c$$$      PARAMETER          ( WILK1 = 0.75e0 )
c$$$      COMPLEX            ZERO, ONE
c$$$      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),
c$$$     $                   ONE = ( 1.0e0, 0.0e0 ) )
c$$$      REAL               TWO
c$$$      PARAMETER          ( TWO = 2.0e0 )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      COMPLEX            AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
c$$$      REAL               S
c$$$      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
c$$$     $                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS,
c$$$     $                   LWKOPT, NDFL, NH, NHO, NIBBLE, NMIN, NS, NSMAX,
c$$$     $                   NSR, NVE, NW, NWMAX, NWR
c$$$      LOGICAL            NWINC, SORTED
c$$$      CHARACTER          JBCMPZ*2
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      INTEGER            ILAENV
c$$$      EXTERNAL           ILAENV
c$$$*     ..
c$$$*     .. Local Arrays ..
c$$$      COMPLEX            ZDUM( 1, 1 )
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CLACPY, CLAHQR, CLAQR2, CLAQR5
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          ABS, AIMAG, CMPLX, INT, MAX, MIN, MOD, REAL,
c$$$     $                   SQRT
c$$$*     ..
c$$$*     .. Statement Functions ..
c$$$      REAL               CABS1
c$$$*     ..
c$$$*     .. Statement Function definitions ..
c$$$      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$      INFO = 0
c$$$*
c$$$*     ==== Quick return for N = 0: nothing to do. ====
c$$$*
c$$$      IF( N.EQ.0 ) THEN
c$$$         WORK( 1 ) = ONE
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     ==== Set up job flags for ILAENV. ====
c$$$*
c$$$      IF( WANTT ) THEN
c$$$         JBCMPZ( 1: 1 ) = 'S'
c$$$      ELSE
c$$$         JBCMPZ( 1: 1 ) = 'E'
c$$$      END IF
c$$$      IF( WANTZ ) THEN
c$$$         JBCMPZ( 2: 2 ) = 'V'
c$$$      ELSE
c$$$         JBCMPZ( 2: 2 ) = 'N'
c$$$      END IF
c$$$*
c$$$*     ==== Tiny matrices must use CLAHQR. ====
c$$$*
c$$$      IF( N.LE.NTINY ) THEN
c$$$*
c$$$*        ==== Estimate optimal workspace. ====
c$$$*
c$$$         LWKOPT = 1
c$$$         IF( LWORK.NE.-1 )
c$$$     $      CALL CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
c$$$     $                   IHIZ, Z, LDZ, INFO )
c$$$      ELSE
c$$$*
c$$$*        ==== Use small bulge multi-shift QR with aggressive early
c$$$*        .    deflation on larger-than-tiny matrices. ====
c$$$*
c$$$*        ==== Hope for the best. ====
c$$$*
c$$$         INFO = 0
c$$$*
c$$$*        ==== NWR = recommended deflation window size.  At this
c$$$*        .    point,  N .GT. NTINY = 11, so there is enough
c$$$*        .    subdiagonal workspace for NWR.GE.2 as required.
c$$$*        .    (In fact, there is enough subdiagonal space for
c$$$*        .    NWR.GE.3.) ====
c$$$*
c$$$         NWR = ILAENV( 13, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         NWR = MAX( 2, NWR )
c$$$         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
c$$$         NW = NWR
c$$$*
c$$$*        ==== NSR = recommended number of simultaneous shifts.
c$$$*        .    At this point N .GT. NTINY = 11, so there is at
c$$$*        .    enough subdiagonal workspace for NSR to be even
c$$$*        .    and greater than or equal to two as required. ====
c$$$*
c$$$         NSR = ILAENV( 15, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
c$$$         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
c$$$*
c$$$*        ==== Estimate optimal workspace ====
c$$$*
c$$$*        ==== Workspace query call to CLAQR2 ====
c$$$*
c$$$         CALL CLAQR2( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ,
c$$$     $                IHIZ, Z, LDZ, LS, LD, W, H, LDH, N, H, LDH, N, H,
c$$$     $                LDH, WORK, -1 )
c$$$*
c$$$*        ==== Optimal workspace = MAX(CLAQR5, CLAQR2) ====
c$$$*
c$$$         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
c$$$*
c$$$*        ==== Quick return in case of workspace query. ====
c$$$*
c$$$         IF( LWORK.EQ.-1 ) THEN
c$$$            WORK( 1 ) = CMPLX( LWKOPT, 0 )
c$$$            RETURN
c$$$         END IF
c$$$*
c$$$*        ==== CLAHQR/CLAQR0 crossover point ====
c$$$*
c$$$         NMIN = ILAENV( 12, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         NMIN = MAX( NTINY, NMIN )
c$$$*
c$$$*        ==== Nibble crossover point ====
c$$$*
c$$$         NIBBLE = ILAENV( 14, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         NIBBLE = MAX( 0, NIBBLE )
c$$$*
c$$$*        ==== Accumulate reflections during ttswp?  Use block
c$$$*        .    2-by-2 structure during matrix-matrix multiply? ====
c$$$*
c$$$         KACC22 = ILAENV( 16, 'CLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
c$$$         KACC22 = MAX( 0, KACC22 )
c$$$         KACC22 = MIN( 2, KACC22 )
c$$$*
c$$$*        ==== NWMAX = the largest possible deflation window for
c$$$*        .    which there is sufficient workspace. ====
c$$$*
c$$$         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
c$$$*
c$$$*        ==== NSMAX = the Largest number of simultaneous shifts
c$$$*        .    for which there is sufficient workspace. ====
c$$$*
c$$$         NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
c$$$         NSMAX = NSMAX - MOD( NSMAX, 2 )
c$$$*
c$$$*        ==== NDFL: an iteration count restarted at deflation. ====
c$$$*
c$$$         NDFL = 1
c$$$*
c$$$*        ==== ITMAX = iteration limit ====
c$$$*
c$$$         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
c$$$*
c$$$*        ==== Last row and column in the active block ====
c$$$*
c$$$         KBOT = IHI
c$$$*
c$$$*        ==== Main Loop ====
c$$$*
c$$$         DO 70 IT = 1, ITMAX
c$$$*
c$$$*           ==== Done when KBOT falls below ILO ====
c$$$*
c$$$            IF( KBOT.LT.ILO )
c$$$     $         GO TO 80
c$$$*
c$$$*           ==== Locate active block ====
c$$$*
c$$$            DO 10 K = KBOT, ILO + 1, -1
c$$$               IF( H( K, K-1 ).EQ.ZERO )
c$$$     $            GO TO 20
c$$$   10       CONTINUE
c$$$            K = ILO
c$$$   20       CONTINUE
c$$$            KTOP = K
c$$$*
c$$$*           ==== Select deflation window size ====
c$$$*
c$$$            NH = KBOT - KTOP + 1
c$$$            IF( NDFL.LT.KEXNW .OR. NH.LT.NW ) THEN
c$$$*
c$$$*              ==== Typical deflation window.  If possible and
c$$$*              .    advisable, nibble the entire active block.
c$$$*              .    If not, use size NWR or NWR+1 depending upon
c$$$*              .    which has the smaller corresponding subdiagonal
c$$$*              .    entry (a heuristic). ====
c$$$*
c$$$               NWINC = .TRUE.
c$$$               IF( NH.LE.MIN( NMIN, NWMAX ) ) THEN
c$$$                  NW = NH
c$$$               ELSE
c$$$                  NW = MIN( NWR, NH, NWMAX )
c$$$                  IF( NW.LT.NWMAX ) THEN
c$$$                     IF( NW.GE.NH-1 ) THEN
c$$$                        NW = NH
c$$$                     ELSE
c$$$                        KWTOP = KBOT - NW + 1
c$$$                        IF( CABS1( H( KWTOP, KWTOP-1 ) ).GT.
c$$$     $                      CABS1( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
c$$$                     END IF
c$$$                  END IF
c$$$               END IF
c$$$            ELSE
c$$$*
c$$$*              ==== Exceptional deflation window.  If there have
c$$$*              .    been no deflations in KEXNW or more iterations,
c$$$*              .    then vary the deflation window size.   At first,
c$$$*              .    because, larger windows are, in general, more
c$$$*              .    powerful than smaller ones, rapidly increase the
c$$$*              .    window up to the maximum reasonable and possible.
c$$$*              .    Then maybe try a slightly smaller window.  ====
c$$$*
c$$$               IF( NWINC .AND. NW.LT.MIN( NWMAX, NH ) ) THEN
c$$$                  NW = MIN( NWMAX, NH, 2*NW )
c$$$               ELSE
c$$$                  NWINC = .FALSE.
c$$$                  IF( NW.EQ.NH .AND. NH.GT.2 )
c$$$     $               NW = NH - 1
c$$$               END IF
c$$$            END IF
c$$$*
c$$$*           ==== Aggressive early deflation:
c$$$*           .    split workspace under the subdiagonal into
c$$$*           .      - an nw-by-nw work array V in the lower
c$$$*           .        left-hand-corner,
c$$$*           .      - an NW-by-at-least-NW-but-more-is-better
c$$$*           .        (NW-by-NHO) horizontal work array along
c$$$*           .        the bottom edge,
c$$$*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
c$$$*           .        vertical work array along the left-hand-edge.
c$$$*           .        ====
c$$$*
c$$$            KV = N - NW + 1
c$$$            KT = NW + 1
c$$$            NHO = ( N-NW-1 ) - KT + 1
c$$$            KWV = NW + 2
c$$$            NVE = ( N-NW ) - KWV + 1
c$$$*
c$$$*           ==== Aggressive early deflation ====
c$$$*
c$$$            CALL CLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
c$$$     $                   IHIZ, Z, LDZ, LS, LD, W, H( KV, 1 ), LDH, NHO,
c$$$     $                   H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, WORK,
c$$$     $                   LWORK )
c$$$*
c$$$*           ==== Adjust KBOT accounting for new deflations. ====
c$$$*
c$$$            KBOT = KBOT - LD
c$$$*
c$$$*           ==== KS points to the shifts. ====
c$$$*
c$$$            KS = KBOT - LS + 1
c$$$*
c$$$*           ==== Skip an expensive QR sweep if there is a (partly
c$$$*           .    heuristic) reason to expect that many eigenvalues
c$$$*           .    will deflate without it.  Here, the QR sweep is
c$$$*           .    skipped if many eigenvalues have just been deflated
c$$$*           .    or if the remaining active block is small.
c$$$*
c$$$            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT-
c$$$     $          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
c$$$*
c$$$*              ==== NS = nominal number of simultaneous shifts.
c$$$*              .    This may be lowered (slightly) if CLAQR2
c$$$*              .    did not provide that many shifts. ====
c$$$*
c$$$               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
c$$$               NS = NS - MOD( NS, 2 )
c$$$*
c$$$*              ==== If there have been no deflations
c$$$*              .    in a multiple of KEXSH iterations,
c$$$*              .    then try exceptional shifts.
c$$$*              .    Otherwise use shifts provided by
c$$$*              .    CLAQR2 above or from the eigenvalues
c$$$*              .    of a trailing principal submatrix. ====
c$$$*
c$$$               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
c$$$                  KS = KBOT - NS + 1
c$$$                  DO 30 I = KBOT, KS + 1, -2
c$$$                     W( I ) = H( I, I ) + WILK1*CABS1( H( I, I-1 ) )
c$$$                     W( I-1 ) = W( I )
c$$$   30             CONTINUE
c$$$               ELSE
c$$$*
c$$$*                 ==== Got NS/2 or fewer shifts? Use CLAHQR
c$$$*                 .    on a trailing principal submatrix to
c$$$*                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
c$$$*                 .    there is enough space below the subdiagonal
c$$$*                 .    to fit an NS-by-NS scratch array.) ====
c$$$*
c$$$                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
c$$$                     KS = KBOT - NS + 1
c$$$                     KT = N - NS + 1
c$$$                     CALL CLACPY( 'A', NS, NS, H( KS, KS ), LDH,
c$$$     $                            H( KT, 1 ), LDH )
c$$$                     CALL CLAHQR( .false., .false., NS, 1, NS,
c$$$     $                            H( KT, 1 ), LDH, W( KS ), 1, 1, ZDUM,
c$$$     $                            1, INF )
c$$$                     KS = KS + INF
c$$$*
c$$$*                    ==== In case of a rare QR failure use
c$$$*                    .    eigenvalues of the trailing 2-by-2
c$$$*                    .    principal submatrix.  Scale to avoid
c$$$*                    .    overflows, underflows and subnormals.
c$$$*                    .    (The scale factor S can not be zero,
c$$$*                    .    because H(KBOT,KBOT-1) is nonzero.) ====
c$$$*
c$$$                     IF( KS.GE.KBOT ) THEN
c$$$                        S = CABS1( H( KBOT-1, KBOT-1 ) ) +
c$$$     $                      CABS1( H( KBOT, KBOT-1 ) ) +
c$$$     $                      CABS1( H( KBOT-1, KBOT ) ) +
c$$$     $                      CABS1( H( KBOT, KBOT ) )
c$$$                        AA = H( KBOT-1, KBOT-1 ) / S
c$$$                        CC = H( KBOT, KBOT-1 ) / S
c$$$                        BB = H( KBOT-1, KBOT ) / S
c$$$                        DD = H( KBOT, KBOT ) / S
c$$$                        TR2 = ( AA+DD ) / TWO
c$$$                        DET = ( AA-TR2 )*( DD-TR2 ) - BB*CC
c$$$                        RTDISC = SQRT( -DET )
c$$$                        W( KBOT-1 ) = ( TR2+RTDISC )*S
c$$$                        W( KBOT ) = ( TR2-RTDISC )*S
c$$$*
c$$$                        KS = KBOT - 1
c$$$                     END IF
c$$$                  END IF
c$$$*
c$$$                  IF( KBOT-KS+1.GT.NS ) THEN
c$$$*
c$$$*                    ==== Sort the shifts (Helps a little) ====
c$$$*
c$$$                     SORTED = .false.
c$$$                     DO 50 K = KBOT, KS + 1, -1
c$$$                        IF( SORTED )
c$$$     $                     GO TO 60
c$$$                        SORTED = .true.
c$$$                        DO 40 I = KS, K - 1
c$$$                           IF( CABS1( W( I ) ).LT.CABS1( W( I+1 ) ) )
c$$$     $                          THEN
c$$$                              SORTED = .false.
c$$$                              SWAP = W( I )
c$$$                              W( I ) = W( I+1 )
c$$$                              W( I+1 ) = SWAP
c$$$                           END IF
c$$$   40                   CONTINUE
c$$$   50                CONTINUE
c$$$   60                CONTINUE
c$$$                  END IF
c$$$               END IF
c$$$*
c$$$*              ==== If there are only two shifts, then use
c$$$*              .    only one.  ====
c$$$*
c$$$               IF( KBOT-KS+1.EQ.2 ) THEN
c$$$                  IF( CABS1( W( KBOT )-H( KBOT, KBOT ) ).LT.
c$$$     $                CABS1( W( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
c$$$                     W( KBOT-1 ) = W( KBOT )
c$$$                  ELSE
c$$$                     W( KBOT ) = W( KBOT-1 )
c$$$                  END IF
c$$$               END IF
c$$$*
c$$$*              ==== Use up to NS of the the smallest magnatiude
c$$$*              .    shifts.  If there aren't NS shifts available,
c$$$*              .    then use them all, possibly dropping one to
c$$$*              .    make the number of shifts even. ====
c$$$*
c$$$               NS = MIN( NS, KBOT-KS+1 )
c$$$               NS = NS - MOD( NS, 2 )
c$$$               KS = KBOT - NS + 1
c$$$*
c$$$*              ==== Small-bulge multi-shift QR sweep:
c$$$*              .    split workspace under the subdiagonal into
c$$$*              .    - a KDU-by-KDU work array U in the lower
c$$$*              .      left-hand-corner,
c$$$*              .    - a KDU-by-at-least-KDU-but-more-is-better
c$$$*              .      (KDU-by-NHo) horizontal work array WH along
c$$$*              .      the bottom edge,
c$$$*              .    - and an at-least-KDU-but-more-is-better-by-KDU
c$$$*              .      (NVE-by-KDU) vertical work WV arrow along
c$$$*              .      the left-hand-edge. ====
c$$$*
c$$$               KDU = 3*NS - 3
c$$$               KU = N - KDU + 1
c$$$               KWH = KDU + 1
c$$$               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
c$$$               KWV = KDU + 4
c$$$               NVE = N - KDU - KWV + 1
c$$$*
c$$$*              ==== Small-bulge multi-shift QR sweep ====
c$$$*
c$$$               CALL CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS,
c$$$     $                      W( KS ), H, LDH, ILOZ, IHIZ, Z, LDZ, WORK,
c$$$     $                      3, H( KU, 1 ), LDH, NVE, H( KWV, 1 ), LDH,
c$$$     $                      NHO, H( KU, KWH ), LDH )
c$$$            END IF
c$$$*
c$$$*           ==== Note progress (or the lack of it). ====
c$$$*
c$$$            IF( LD.GT.0 ) THEN
c$$$               NDFL = 1
c$$$            ELSE
c$$$               NDFL = NDFL + 1
c$$$            END IF
c$$$*
c$$$*           ==== End of main loop ====
c$$$   70    CONTINUE
c$$$*
c$$$*        ==== Iteration limit exceeded.  Set INFO to show where
c$$$*        .    the problem occurred and exit. ====
c$$$*
c$$$         INFO = KBOT
c$$$   80    CONTINUE
c$$$      END IF
c$$$*
c$$$*     ==== Return the optimal value of LWORK. ====
c$$$*
c$$$      WORK( 1 ) = CMPLX( LWKOPT, 0 )
c$$$*
c$$$*     ==== End of CLAQR4 ====
c$$$*
c$$$      END
c$$$      SUBROUTINE CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,
c$$$     $                   H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,
c$$$     $                   WV, LDWV, NH, WH, LDWH )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
c$$$     $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
c$$$      LOGICAL            WANTT, WANTZ
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ),
c$$$     $                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
c$$$*     ..
c$$$*
c$$$*     This auxiliary subroutine called by CLAQR0 performs a
c$$$*     single small-bulge multi-shift QR sweep.
c$$$*
c$$$*      WANTT  (input) logical scalar
c$$$*             WANTT = .true. if the triangular Schur factor
c$$$*             is being computed.  WANTT is set to .false. otherwise.
c$$$*
c$$$*      WANTZ  (input) logical scalar
c$$$*             WANTZ = .true. if the unitary Schur factor is being
c$$$*             computed.  WANTZ is set to .false. otherwise.
c$$$*
c$$$*      KACC22 (input) integer with value 0, 1, or 2.
c$$$*             Specifies the computation mode of far-from-diagonal
c$$$*             orthogonal updates.
c$$$*        = 0: CLAQR5 does not accumulate reflections and does not
c$$$*             use matrix-matrix multiply to update far-from-diagonal
c$$$*             matrix entries.
c$$$*        = 1: CLAQR5 accumulates reflections and uses matrix-matrix
c$$$*             multiply to update the far-from-diagonal matrix entries.
c$$$*        = 2: CLAQR5 accumulates reflections, uses matrix-matrix
c$$$*             multiply to update the far-from-diagonal matrix entries,
c$$$*             and takes advantage of 2-by-2 block structure during
c$$$*             matrix multiplies.
c$$$*
c$$$*      N      (input) integer scalar
c$$$*             N is the order of the Hessenberg matrix H upon which this
c$$$*             subroutine operates.
c$$$*
c$$$*      KTOP   (input) integer scalar
c$$$*      KBOT   (input) integer scalar
c$$$*             These are the first and last rows and columns of an
c$$$*             isolated diagonal block upon which the QR sweep is to be
c$$$*             applied. It is assumed without a check that
c$$$*                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0
c$$$*             and
c$$$*                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.
c$$$*
c$$$*      NSHFTS (input) integer scalar
c$$$*             NSHFTS gives the number of simultaneous shifts.  NSHFTS
c$$$*             must be positive and even.
c$$$*
c$$$*      S      (input) COMPLEX array of size (NSHFTS)
c$$$*             S contains the shifts of origin that define the multi-
c$$$*             shift QR sweep.
c$$$*
c$$$*      H      (input/output) COMPLEX array of size (LDH,N)
c$$$*             On input H contains a Hessenberg matrix.  On output a
c$$$*             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
c$$$*             to the isolated diagonal block in rows and columns KTOP
c$$$*             through KBOT.
c$$$*
c$$$*      LDH    (input) integer scalar
c$$$*             LDH is the leading dimension of H just as declared in the
c$$$*             calling procedure.  LDH.GE.MAX(1,N).
c$$$*
c$$$*      ILOZ   (input) INTEGER
c$$$*      IHIZ   (input) INTEGER
c$$$*             Specify the rows of Z to which transformations must be
c$$$*             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N
c$$$*
c$$$*      Z      (input/output) COMPLEX array of size (LDZ,IHI)
c$$$*             If WANTZ = .TRUE., then the QR Sweep unitary
c$$$*             similarity transformation is accumulated into
c$$$*             Z(ILOZ:IHIZ,ILO:IHI) from the right.
c$$$*             If WANTZ = .FALSE., then Z is unreferenced.
c$$$*
c$$$*      LDZ    (input) integer scalar
c$$$*             LDA is the leading dimension of Z just as declared in
c$$$*             the calling procedure. LDZ.GE.N.
c$$$*
c$$$*      V      (workspace) COMPLEX array of size (LDV,NSHFTS/2)
c$$$*
c$$$*      LDV    (input) integer scalar
c$$$*             LDV is the leading dimension of V as declared in the
c$$$*             calling procedure.  LDV.GE.3.
c$$$*
c$$$*      U      (workspace) COMPLEX array of size
c$$$*             (LDU,3*NSHFTS-3)
c$$$*
c$$$*      LDU    (input) integer scalar
c$$$*             LDU is the leading dimension of U just as declared in the
c$$$*             in the calling subroutine.  LDU.GE.3*NSHFTS-3.
c$$$*
c$$$*      NH     (input) integer scalar
c$$$*             NH is the number of columns in array WH available for
c$$$*             workspace. NH.GE.1.
c$$$*
c$$$*      WH     (workspace) COMPLEX array of size (LDWH,NH)
c$$$*
c$$$*      LDWH   (input) integer scalar
c$$$*             Leading dimension of WH just as declared in the
c$$$*             calling procedure.  LDWH.GE.3*NSHFTS-3.
c$$$*
c$$$*      NV     (input) integer scalar
c$$$*             NV is the number of rows in WV agailable for workspace.
c$$$*             NV.GE.1.
c$$$*
c$$$*      WV     (workspace) COMPLEX array of size
c$$$*             (LDWV,3*NSHFTS-3)
c$$$*
c$$$*      LDWV   (input) integer scalar
c$$$*             LDWV is the leading dimension of WV as declared in the
c$$$*             in the calling subroutine.  LDWV.GE.NV.
c$$$*
c$$$*
c$$$*     ================================================================
c$$$*     Based on contributions by
c$$$*        Karen Braman and Ralph Byers, Department of Mathematics,
c$$$*        University of Kansas, USA
c$$$*
c$$$*     ============================================================
c$$$*     Reference:
c$$$*
c$$$*     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
c$$$*     Algorithm Part I: Maintaining Well Focused Shifts, and
c$$$*     Level 3 Performance, SIAM Journal of Matrix Analysis,
c$$$*     volume 23, pages 929--947, 2002.
c$$$*
c$$$*     ============================================================
c$$$*     .. Parameters ..
c$$$      COMPLEX            ZERO, ONE
c$$$      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),
c$$$     $                   ONE = ( 1.0e0, 0.0e0 ) )
c$$$      REAL               RZERO, RONE
c$$$      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0 )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      COMPLEX            ALPHA, BETA, CDUM, REFSUM
c$$$      REAL               H11, H12, H21, H22, SAFMAX, SAFMIN, SCL,
c$$$     $                   SMLNUM, TST1, TST2, ULP
c$$$      INTEGER            I2, I4, INCOL, J, J2, J4, JBOT, JCOL, JLEN,
c$$$     $                   JROW, JTOP, K, K1, KDU, KMS, KNZ, KRCOL, KZS,
c$$$     $                   M, M22, MBOT, MEND, MSTART, MTOP, NBMPS, NDCOL,
c$$$     $                   NS, NU
c$$$      LOGICAL            ACCUM, BLK22, BMP22
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      REAL               SLAMCH
c$$$      EXTERNAL           SLAMCH
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$*
c$$$      INTRINSIC          ABS, AIMAG, CONJG, MAX, MIN, MOD, REAL
c$$$*     ..
c$$$*     .. Local Arrays ..
c$$$      COMPLEX            VT( 3 )
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CGEMM, CLACPY, CLAQR1, CLARFG, CLASET, CTRMM,
c$$$     $                   SLABAD
c$$$*     ..
c$$$*     .. Statement Functions ..
c$$$      REAL               CABS1
c$$$*     ..
c$$$*     .. Statement Function definitions ..
c$$$      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     ==== If there are no shifts, then there is nothing to do. ====
c$$$*
c$$$      IF( NSHFTS.LT.2 )
c$$$     $   RETURN
c$$$*
c$$$*     ==== If the active block is empty or 1-by-1, then there
c$$$*     .    is nothing to do. ====
c$$$*
c$$$      IF( KTOP.GE.KBOT )
c$$$     $   RETURN
c$$$*
c$$$*     ==== NSHFTS is supposed to be even, but if is odd,
c$$$*     .    then simply reduce it by one.  ====
c$$$*
c$$$      NS = NSHFTS - MOD( NSHFTS, 2 )
c$$$*
c$$$*     ==== Machine constants for deflation ====
c$$$*
c$$$      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
c$$$      SAFMAX = RONE / SAFMIN
c$$$      CALL SLABAD( SAFMIN, SAFMAX )
c$$$      ULP = SLAMCH( 'PRECISION' )
c$$$      SMLNUM = SAFMIN*( REAL( N ) / ULP )
c$$$*
c$$$*     ==== Use accumulated reflections to update far-from-diagonal
c$$$*     .    entries ? ====
c$$$*
c$$$      ACCUM = ( KACC22.EQ.1 ) .OR. ( KACC22.EQ.2 )
c$$$*
c$$$*     ==== If so, exploit the 2-by-2 block structure? ====
c$$$*
c$$$      BLK22 = ( NS.GT.2 ) .AND. ( KACC22.EQ.2 )
c$$$*
c$$$*     ==== clear trash ====
c$$$*
c$$$      IF( KTOP+2.LE.KBOT )
c$$$     $   H( KTOP+2, KTOP ) = ZERO
c$$$*
c$$$*     ==== NBMPS = number of 2-shift bulges in the chain ====
c$$$*
c$$$      NBMPS = NS / 2
c$$$*
c$$$*     ==== KDU = width of slab ====
c$$$*
c$$$      KDU = 6*NBMPS - 3
c$$$*
c$$$*     ==== Create and chase chains of NBMPS bulges ====
c$$$*
c$$$      DO 210 INCOL = 3*( 1-NBMPS ) + KTOP - 1, KBOT - 2, 3*NBMPS - 2
c$$$         NDCOL = INCOL + KDU
c$$$         IF( ACCUM )
c$$$     $      CALL CLASET( 'ALL', KDU, KDU, ZERO, ONE, U, LDU )
c$$$*
c$$$*        ==== Near-the-diagonal bulge chase.  The following loop
c$$$*        .    performs the near-the-diagonal part of a small bulge
c$$$*        .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
c$$$*        .    chunk extends from column INCOL to column NDCOL
c$$$*        .    (including both column INCOL and column NDCOL). The
c$$$*        .    following loop chases a 3*NBMPS column long chain of
c$$$*        .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
c$$$*        .    may be less than KTOP and and NDCOL may be greater than
c$$$*        .    KBOT indicating phantom columns from which to chase
c$$$*        .    bulges before they are actually introduced or to which
c$$$*        .    to chase bulges beyond column KBOT.)  ====
c$$$*
c$$$         DO 140 KRCOL = INCOL, MIN( INCOL+3*NBMPS-3, KBOT-2 )
c$$$*
c$$$*           ==== Bulges number MTOP to MBOT are active double implicit
c$$$*           .    shift bulges.  There may or may not also be small
c$$$*           .    2-by-2 bulge, if there is room.  The inactive bulges
c$$$*           .    (if any) must wait until the active bulges have moved
c$$$*           .    down the diagonal to make room.  The phantom matrix
c$$$*           .    paradigm described above helps keep track.  ====
c$$$*
c$$$            MTOP = MAX( 1, ( ( KTOP-1 )-KRCOL+2 ) / 3+1 )
c$$$            MBOT = MIN( NBMPS, ( KBOT-KRCOL ) / 3 )
c$$$            M22 = MBOT + 1
c$$$            BMP22 = ( MBOT.LT.NBMPS ) .AND. ( KRCOL+3*( M22-1 ) ).EQ.
c$$$     $              ( KBOT-2 )
c$$$*
c$$$*           ==== Generate reflections to chase the chain right
c$$$*           .    one column.  (The minimum value of K is KTOP-1.) ====
c$$$*
c$$$            DO 10 M = MTOP, MBOT
c$$$               K = KRCOL + 3*( M-1 )
c$$$               IF( K.EQ.KTOP-1 ) THEN
c$$$                  CALL CLAQR1( 3, H( KTOP, KTOP ), LDH, S( 2*M-1 ),
c$$$     $                         S( 2*M ), V( 1, M ) )
c$$$                  ALPHA = V( 1, M )
c$$$                  CALL CLARFG( 3, ALPHA, V( 2, M ), 1, V( 1, M ) )
c$$$               ELSE
c$$$                  BETA = H( K+1, K )
c$$$                  V( 2, M ) = H( K+2, K )
c$$$                  V( 3, M ) = H( K+3, K )
c$$$                  CALL CLARFG( 3, BETA, V( 2, M ), 1, V( 1, M ) )
c$$$*
c$$$*                 ==== A Bulge may collapse because of vigilant
c$$$*                 .    deflation or destructive underflow.  (The
c$$$*                 .    initial bulge is always collapsed.) Use
c$$$*                 .    the two-small-subdiagonals trick to try
c$$$*                 .    to get it started again. If V(2,M).NE.0 and
c$$$*                 .    V(3,M) = H(K+3,K+1) = H(K+3,K+2) = 0, then
c$$$*                 .    this bulge is collapsing into a zero
c$$$*                 .    subdiagonal.  It will be restarted next
c$$$*                 .    trip through the loop.)
c$$$*
c$$$                  IF( V( 1, M ).NE.ZERO .AND.
c$$$     $                ( V( 3, M ).NE.ZERO .OR. ( H( K+3,
c$$$     $                K+1 ).EQ.ZERO .AND. H( K+3, K+2 ).EQ.ZERO ) ) )
c$$$     $                 THEN
c$$$*
c$$$*                    ==== Typical case: not collapsed (yet). ====
c$$$*
c$$$                     H( K+1, K ) = BETA
c$$$                     H( K+2, K ) = ZERO
c$$$                     H( K+3, K ) = ZERO
c$$$                  ELSE
c$$$*
c$$$*                    ==== Atypical case: collapsed.  Attempt to
c$$$*                    .    reintroduce ignoring H(K+1,K).  If the
c$$$*                    .    fill resulting from the new reflector
c$$$*                    .    is too large, then abandon it.
c$$$*                    .    Otherwise, use the new one. ====
c$$$*
c$$$                     CALL CLAQR1( 3, H( K+1, K+1 ), LDH, S( 2*M-1 ),
c$$$     $                            S( 2*M ), VT )
c$$$                     SCL = CABS1( VT( 1 ) ) + CABS1( VT( 2 ) ) +
c$$$     $                     CABS1( VT( 3 ) )
c$$$                     IF( SCL.NE.RZERO ) THEN
c$$$                        VT( 1 ) = VT( 1 ) / SCL
c$$$                        VT( 2 ) = VT( 2 ) / SCL
c$$$                        VT( 3 ) = VT( 3 ) / SCL
c$$$                     END IF
c$$$*
c$$$*                    ==== The following is the traditional and
c$$$*                    .    conservative two-small-subdiagonals
c$$$*                    .    test.  ====
c$$$*                    .
c$$$                     IF( CABS1( H( K+1, K ) )*
c$$$     $                   ( CABS1( VT( 2 ) )+CABS1( VT( 3 ) ) ).GT.ULP*
c$$$     $                   CABS1( VT( 1 ) )*( CABS1( H( K,
c$$$     $                   K ) )+CABS1( H( K+1, K+1 ) )+CABS1( H( K+2,
c$$$     $                   K+2 ) ) ) ) THEN
c$$$*
c$$$*                       ==== Starting a new bulge here would
c$$$*                       .    create non-negligible fill.   If
c$$$*                       .    the old reflector is diagonal (only
c$$$*                       .    possible with underflows), then
c$$$*                       .    change it to I.  Otherwise, use
c$$$*                       .    it with trepidation. ====
c$$$*
c$$$                        IF( V( 2, M ).EQ.ZERO .AND. V( 3, M ).EQ.ZERO )
c$$$     $                       THEN
c$$$                           V( 1, M ) = ZERO
c$$$                        ELSE
c$$$                           H( K+1, K ) = BETA
c$$$                           H( K+2, K ) = ZERO
c$$$                           H( K+3, K ) = ZERO
c$$$                        END IF
c$$$                     ELSE
c$$$*
c$$$*                       ==== Stating a new bulge here would
c$$$*                       .    create only negligible fill.
c$$$*                       .    Replace the old reflector with
c$$$*                       .    the new one. ====
c$$$*
c$$$                        ALPHA = VT( 1 )
c$$$                        CALL CLARFG( 3, ALPHA, VT( 2 ), 1, VT( 1 ) )
c$$$                        REFSUM = H( K+1, K ) +
c$$$     $                           H( K+2, K )*CONJG( VT( 2 ) ) +
c$$$     $                           H( K+3, K )*CONJG( VT( 3 ) )
c$$$                        H( K+1, K ) = H( K+1, K ) -
c$$$     $                                CONJG( VT( 1 ) )*REFSUM
c$$$                        H( K+2, K ) = ZERO
c$$$                        H( K+3, K ) = ZERO
c$$$                        V( 1, M ) = VT( 1 )
c$$$                        V( 2, M ) = VT( 2 )
c$$$                        V( 3, M ) = VT( 3 )
c$$$                     END IF
c$$$                  END IF
c$$$               END IF
c$$$   10       CONTINUE
c$$$*
c$$$*           ==== Generate a 2-by-2 reflection, if needed. ====
c$$$*
c$$$            K = KRCOL + 3*( M22-1 )
c$$$            IF( BMP22 ) THEN
c$$$               IF( K.EQ.KTOP-1 ) THEN
c$$$                  CALL CLAQR1( 2, H( K+1, K+1 ), LDH, S( 2*M22-1 ),
c$$$     $                         S( 2*M22 ), V( 1, M22 ) )
c$$$                  BETA = V( 1, M22 )
c$$$                  CALL CLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
c$$$               ELSE
c$$$                  BETA = H( K+1, K )
c$$$                  V( 2, M22 ) = H( K+2, K )
c$$$                  CALL CLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
c$$$                  H( K+1, K ) = BETA
c$$$                  H( K+2, K ) = ZERO
c$$$               END IF
c$$$            ELSE
c$$$*
c$$$*              ==== Initialize V(1,M22) here to avoid possible undefined
c$$$*              .    variable problems later. ====
c$$$*
c$$$               V( 1, M22 ) = ZERO
c$$$            END IF
c$$$*
c$$$*           ==== Multiply H by reflections from the left ====
c$$$*
c$$$            IF( ACCUM ) THEN
c$$$               JBOT = MIN( NDCOL, KBOT )
c$$$            ELSE IF( WANTT ) THEN
c$$$               JBOT = N
c$$$            ELSE
c$$$               JBOT = KBOT
c$$$            END IF
c$$$            DO 30 J = MAX( KTOP, KRCOL ), JBOT
c$$$               MEND = MIN( MBOT, ( J-KRCOL+2 ) / 3 )
c$$$               DO 20 M = MTOP, MEND
c$$$                  K = KRCOL + 3*( M-1 )
c$$$                  REFSUM = CONJG( V( 1, M ) )*
c$$$     $                     ( H( K+1, J )+CONJG( V( 2, M ) )*H( K+2, J )+
c$$$     $                     CONJG( V( 3, M ) )*H( K+3, J ) )
c$$$                  H( K+1, J ) = H( K+1, J ) - REFSUM
c$$$                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M )
c$$$                  H( K+3, J ) = H( K+3, J ) - REFSUM*V( 3, M )
c$$$   20          CONTINUE
c$$$   30       CONTINUE
c$$$            IF( BMP22 ) THEN
c$$$               K = KRCOL + 3*( M22-1 )
c$$$               DO 40 J = MAX( K+1, KTOP ), JBOT
c$$$                  REFSUM = CONJG( V( 1, M22 ) )*
c$$$     $                     ( H( K+1, J )+CONJG( V( 2, M22 ) )*
c$$$     $                     H( K+2, J ) )
c$$$                  H( K+1, J ) = H( K+1, J ) - REFSUM
c$$$                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M22 )
c$$$   40          CONTINUE
c$$$            END IF
c$$$*
c$$$*           ==== Multiply H by reflections from the right.
c$$$*           .    Delay filling in the last row until the
c$$$*           .    vigilant deflation check is complete. ====
c$$$*
c$$$            IF( ACCUM ) THEN
c$$$               JTOP = MAX( KTOP, INCOL )
c$$$            ELSE IF( WANTT ) THEN
c$$$               JTOP = 1
c$$$            ELSE
c$$$               JTOP = KTOP
c$$$            END IF
c$$$            DO 80 M = MTOP, MBOT
c$$$               IF( V( 1, M ).NE.ZERO ) THEN
c$$$                  K = KRCOL + 3*( M-1 )
c$$$                  DO 50 J = JTOP, MIN( KBOT, K+3 )
c$$$                     REFSUM = V( 1, M )*( H( J, K+1 )+V( 2, M )*
c$$$     $                        H( J, K+2 )+V( 3, M )*H( J, K+3 ) )
c$$$                     H( J, K+1 ) = H( J, K+1 ) - REFSUM
c$$$                     H( J, K+2 ) = H( J, K+2 ) -
c$$$     $                             REFSUM*CONJG( V( 2, M ) )
c$$$                     H( J, K+3 ) = H( J, K+3 ) -
c$$$     $                             REFSUM*CONJG( V( 3, M ) )
c$$$   50             CONTINUE
c$$$*
c$$$                  IF( ACCUM ) THEN
c$$$*
c$$$*                    ==== Accumulate U. (If necessary, update Z later
c$$$*                    .    with with an efficient matrix-matrix
c$$$*                    .    multiply.) ====
c$$$*
c$$$                     KMS = K - INCOL
c$$$                     DO 60 J = MAX( 1, KTOP-INCOL ), KDU
c$$$                        REFSUM = V( 1, M )*( U( J, KMS+1 )+V( 2, M )*
c$$$     $                           U( J, KMS+2 )+V( 3, M )*U( J, KMS+3 ) )
c$$$                        U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
c$$$                        U( J, KMS+2 ) = U( J, KMS+2 ) -
c$$$     $                                  REFSUM*CONJG( V( 2, M ) )
c$$$                        U( J, KMS+3 ) = U( J, KMS+3 ) -
c$$$     $                                  REFSUM*CONJG( V( 3, M ) )
c$$$   60                CONTINUE
c$$$                  ELSE IF( WANTZ ) THEN
c$$$*
c$$$*                    ==== U is not accumulated, so update Z
c$$$*                    .    now by multiplying by reflections
c$$$*                    .    from the right. ====
c$$$*
c$$$                     DO 70 J = ILOZ, IHIZ
c$$$                        REFSUM = V( 1, M )*( Z( J, K+1 )+V( 2, M )*
c$$$     $                           Z( J, K+2 )+V( 3, M )*Z( J, K+3 ) )
c$$$                        Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
c$$$                        Z( J, K+2 ) = Z( J, K+2 ) -
c$$$     $                                REFSUM*CONJG( V( 2, M ) )
c$$$                        Z( J, K+3 ) = Z( J, K+3 ) -
c$$$     $                                REFSUM*CONJG( V( 3, M ) )
c$$$   70                CONTINUE
c$$$                  END IF
c$$$               END IF
c$$$   80       CONTINUE
c$$$*
c$$$*           ==== Special case: 2-by-2 reflection (if needed) ====
c$$$*
c$$$            K = KRCOL + 3*( M22-1 )
c$$$            IF( BMP22 .AND. ( V( 1, M22 ).NE.ZERO ) ) THEN
c$$$               DO 90 J = JTOP, MIN( KBOT, K+3 )
c$$$                  REFSUM = V( 1, M22 )*( H( J, K+1 )+V( 2, M22 )*
c$$$     $                     H( J, K+2 ) )
c$$$                  H( J, K+1 ) = H( J, K+1 ) - REFSUM
c$$$                  H( J, K+2 ) = H( J, K+2 ) -
c$$$     $                          REFSUM*CONJG( V( 2, M22 ) )
c$$$   90          CONTINUE
c$$$*
c$$$               IF( ACCUM ) THEN
c$$$                  KMS = K - INCOL
c$$$                  DO 100 J = MAX( 1, KTOP-INCOL ), KDU
c$$$                     REFSUM = V( 1, M22 )*( U( J, KMS+1 )+V( 2, M22 )*
c$$$     $                        U( J, KMS+2 ) )
c$$$                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
c$$$                     U( J, KMS+2 ) = U( J, KMS+2 ) -
c$$$     $                               REFSUM*CONJG( V( 2, M22 ) )
c$$$  100             CONTINUE
c$$$               ELSE IF( WANTZ ) THEN
c$$$                  DO 110 J = ILOZ, IHIZ
c$$$                     REFSUM = V( 1, M22 )*( Z( J, K+1 )+V( 2, M22 )*
c$$$     $                        Z( J, K+2 ) )
c$$$                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
c$$$                     Z( J, K+2 ) = Z( J, K+2 ) -
c$$$     $                             REFSUM*CONJG( V( 2, M22 ) )
c$$$  110             CONTINUE
c$$$               END IF
c$$$            END IF
c$$$*
c$$$*           ==== Vigilant deflation check ====
c$$$*
c$$$            MSTART = MTOP
c$$$            IF( KRCOL+3*( MSTART-1 ).LT.KTOP )
c$$$     $         MSTART = MSTART + 1
c$$$            MEND = MBOT
c$$$            IF( BMP22 )
c$$$     $         MEND = MEND + 1
c$$$            IF( KRCOL.EQ.KBOT-2 )
c$$$     $         MEND = MEND + 1
c$$$            DO 120 M = MSTART, MEND
c$$$               K = MIN( KBOT-1, KRCOL+3*( M-1 ) )
c$$$*
c$$$*              ==== The following convergence test requires that
c$$$*              .    the tradition small-compared-to-nearby-diagonals
c$$$*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
c$$$*              .    criteria both be satisfied.  The latter improves
c$$$*              .    accuracy in some examples. Falling back on an
c$$$*              .    alternate convergence criterion when TST1 or TST2
c$$$*              .    is zero (as done here) is traditional but probably
c$$$*              .    unnecessary. ====
c$$$*
c$$$               IF( H( K+1, K ).NE.ZERO ) THEN
c$$$                  TST1 = CABS1( H( K, K ) ) + CABS1( H( K+1, K+1 ) )
c$$$                  IF( TST1.EQ.RZERO ) THEN
c$$$                     IF( K.GE.KTOP+1 )
c$$$     $                  TST1 = TST1 + CABS1( H( K, K-1 ) )
c$$$                     IF( K.GE.KTOP+2 )
c$$$     $                  TST1 = TST1 + CABS1( H( K, K-2 ) )
c$$$                     IF( K.GE.KTOP+3 )
c$$$     $                  TST1 = TST1 + CABS1( H( K, K-3 ) )
c$$$                     IF( K.LE.KBOT-2 )
c$$$     $                  TST1 = TST1 + CABS1( H( K+2, K+1 ) )
c$$$                     IF( K.LE.KBOT-3 )
c$$$     $                  TST1 = TST1 + CABS1( H( K+3, K+1 ) )
c$$$                     IF( K.LE.KBOT-4 )
c$$$     $                  TST1 = TST1 + CABS1( H( K+4, K+1 ) )
c$$$                  END IF
c$$$                  IF( CABS1( H( K+1, K ) ).LE.MAX( SMLNUM, ULP*TST1 ) )
c$$$     $                 THEN
c$$$                     H12 = MAX( CABS1( H( K+1, K ) ),
c$$$     $                     CABS1( H( K, K+1 ) ) )
c$$$                     H21 = MIN( CABS1( H( K+1, K ) ),
c$$$     $                     CABS1( H( K, K+1 ) ) )
c$$$                     H11 = MAX( CABS1( H( K+1, K+1 ) ),
c$$$     $                     CABS1( H( K, K )-H( K+1, K+1 ) ) )
c$$$                     H22 = MIN( CABS1( H( K+1, K+1 ) ),
c$$$     $                     CABS1( H( K, K )-H( K+1, K+1 ) ) )
c$$$                     SCL = H11 + H12
c$$$                     TST2 = H22*( H11 / SCL )
c$$$*
c$$$                     IF( TST2.EQ.RZERO .OR. H21*( H12 / SCL ).LE.
c$$$     $                   MAX( SMLNUM, ULP*TST2 ) )H( K+1, K ) = ZERO
c$$$                  END IF
c$$$               END IF
c$$$  120       CONTINUE
c$$$*
c$$$*           ==== Fill in the last row of each bulge. ====
c$$$*
c$$$            MEND = MIN( NBMPS, ( KBOT-KRCOL-1 ) / 3 )
c$$$            DO 130 M = MTOP, MEND
c$$$               K = KRCOL + 3*( M-1 )
c$$$               REFSUM = V( 1, M )*V( 3, M )*H( K+4, K+3 )
c$$$               H( K+4, K+1 ) = -REFSUM
c$$$               H( K+4, K+2 ) = -REFSUM*CONJG( V( 2, M ) )
c$$$               H( K+4, K+3 ) = H( K+4, K+3 ) - REFSUM*CONJG( V( 3, M ) )
c$$$  130       CONTINUE
c$$$*
c$$$*           ==== End of near-the-diagonal bulge chase. ====
c$$$*
c$$$  140    CONTINUE
c$$$*
c$$$*        ==== Use U (if accumulated) to update far-from-diagonal
c$$$*        .    entries in H.  If required, use U to update Z as
c$$$*        .    well. ====
c$$$*
c$$$         IF( ACCUM ) THEN
c$$$            IF( WANTT ) THEN
c$$$               JTOP = 1
c$$$               JBOT = N
c$$$            ELSE
c$$$               JTOP = KTOP
c$$$               JBOT = KBOT
c$$$            END IF
c$$$            IF( ( .NOT.BLK22 ) .OR. ( INCOL.LT.KTOP ) .OR.
c$$$     $          ( NDCOL.GT.KBOT ) .OR. ( NS.LE.2 ) ) THEN
c$$$*
c$$$*              ==== Updates not exploiting the 2-by-2 block
c$$$*              .    structure of U.  K1 and NU keep track of
c$$$*              .    the location and size of U in the special
c$$$*              .    cases of introducing bulges and chasing
c$$$*              .    bulges off the bottom.  In these special
c$$$*              .    cases and in case the number of shifts
c$$$*              .    is NS = 2, there is no 2-by-2 block
c$$$*              .    structure to exploit.  ====
c$$$*
c$$$               K1 = MAX( 1, KTOP-INCOL )
c$$$               NU = ( KDU-MAX( 0, NDCOL-KBOT ) ) - K1 + 1
c$$$*
c$$$*              ==== Horizontal Multiply ====
c$$$*
c$$$               DO 150 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
c$$$                  JLEN = MIN( NH, JBOT-JCOL+1 )
c$$$                  CALL CGEMM( 'C', 'N', NU, JLEN, NU, ONE, U( K1, K1 ),
c$$$     $                        LDU, H( INCOL+K1, JCOL ), LDH, ZERO, WH,
c$$$     $                        LDWH )
c$$$                  CALL CLACPY( 'ALL', NU, JLEN, WH, LDWH,
c$$$     $                         H( INCOL+K1, JCOL ), LDH )
c$$$  150          CONTINUE
c$$$*
c$$$*              ==== Vertical multiply ====
c$$$*
c$$$               DO 160 JROW = JTOP, MAX( KTOP, INCOL ) - 1, NV
c$$$                  JLEN = MIN( NV, MAX( KTOP, INCOL )-JROW )
c$$$                  CALL CGEMM( 'N', 'N', JLEN, NU, NU, ONE,
c$$$     $                        H( JROW, INCOL+K1 ), LDH, U( K1, K1 ),
c$$$     $                        LDU, ZERO, WV, LDWV )
c$$$                  CALL CLACPY( 'ALL', JLEN, NU, WV, LDWV,
c$$$     $                         H( JROW, INCOL+K1 ), LDH )
c$$$  160          CONTINUE
c$$$*
c$$$*              ==== Z multiply (also vertical) ====
c$$$*
c$$$               IF( WANTZ ) THEN
c$$$                  DO 170 JROW = ILOZ, IHIZ, NV
c$$$                     JLEN = MIN( NV, IHIZ-JROW+1 )
c$$$                     CALL CGEMM( 'N', 'N', JLEN, NU, NU, ONE,
c$$$     $                           Z( JROW, INCOL+K1 ), LDZ, U( K1, K1 ),
c$$$     $                           LDU, ZERO, WV, LDWV )
c$$$                     CALL CLACPY( 'ALL', JLEN, NU, WV, LDWV,
c$$$     $                            Z( JROW, INCOL+K1 ), LDZ )
c$$$  170             CONTINUE
c$$$               END IF
c$$$            ELSE
c$$$*
c$$$*              ==== Updates exploiting U's 2-by-2 block structure.
c$$$*              .    (I2, I4, J2, J4 are the last rows and columns
c$$$*              .    of the blocks.) ====
c$$$*
c$$$               I2 = ( KDU+1 ) / 2
c$$$               I4 = KDU
c$$$               J2 = I4 - I2
c$$$               J4 = KDU
c$$$*
c$$$*              ==== KZS and KNZ deal with the band of zeros
c$$$*              .    along the diagonal of one of the triangular
c$$$*              .    blocks. ====
c$$$*
c$$$               KZS = ( J4-J2 ) - ( NS+1 )
c$$$               KNZ = NS + 1
c$$$*
c$$$*              ==== Horizontal multiply ====
c$$$*
c$$$               DO 180 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
c$$$                  JLEN = MIN( NH, JBOT-JCOL+1 )
c$$$*
c$$$*                 ==== Copy bottom of H to top+KZS of scratch ====
c$$$*                  (The first KZS rows get multiplied by zero.) ====
c$$$*
c$$$                  CALL CLACPY( 'ALL', KNZ, JLEN, H( INCOL+1+J2, JCOL ),
c$$$     $                         LDH, WH( KZS+1, 1 ), LDWH )
c$$$*
c$$$*                 ==== Multiply by U21' ====
c$$$*
c$$$                  CALL CLASET( 'ALL', KZS, JLEN, ZERO, ZERO, WH, LDWH )
c$$$                  CALL CTRMM( 'L', 'U', 'C', 'N', KNZ, JLEN, ONE,
c$$$     $                        U( J2+1, 1+KZS ), LDU, WH( KZS+1, 1 ),
c$$$     $                        LDWH )
c$$$*
c$$$*                 ==== Multiply top of H by U11' ====
c$$$*
c$$$                  CALL CGEMM( 'C', 'N', I2, JLEN, J2, ONE, U, LDU,
c$$$     $                        H( INCOL+1, JCOL ), LDH, ONE, WH, LDWH )
c$$$*
c$$$*                 ==== Copy top of H bottom of WH ====
c$$$*
c$$$                  CALL CLACPY( 'ALL', J2, JLEN, H( INCOL+1, JCOL ), LDH,
c$$$     $                         WH( I2+1, 1 ), LDWH )
c$$$*
c$$$*                 ==== Multiply by U21' ====
c$$$*
c$$$                  CALL CTRMM( 'L', 'L', 'C', 'N', J2, JLEN, ONE,
c$$$     $                        U( 1, I2+1 ), LDU, WH( I2+1, 1 ), LDWH )
c$$$*
c$$$*                 ==== Multiply by U22 ====
c$$$*
c$$$                  CALL CGEMM( 'C', 'N', I4-I2, JLEN, J4-J2, ONE,
c$$$     $                        U( J2+1, I2+1 ), LDU,
c$$$     $                        H( INCOL+1+J2, JCOL ), LDH, ONE,
c$$$     $                        WH( I2+1, 1 ), LDWH )
c$$$*
c$$$*                 ==== Copy it back ====
c$$$*
c$$$                  CALL CLACPY( 'ALL', KDU, JLEN, WH, LDWH,
c$$$     $                         H( INCOL+1, JCOL ), LDH )
c$$$  180          CONTINUE
c$$$*
c$$$*              ==== Vertical multiply ====
c$$$*
c$$$               DO 190 JROW = JTOP, MAX( INCOL, KTOP ) - 1, NV
c$$$                  JLEN = MIN( NV, MAX( INCOL, KTOP )-JROW )
c$$$*
c$$$*                 ==== Copy right of H to scratch (the first KZS
c$$$*                 .    columns get multiplied by zero) ====
c$$$*
c$$$                  CALL CLACPY( 'ALL', JLEN, KNZ, H( JROW, INCOL+1+J2 ),
c$$$     $                         LDH, WV( 1, 1+KZS ), LDWV )
c$$$*
c$$$*                 ==== Multiply by U21 ====
c$$$*
c$$$                  CALL CLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV, LDWV )
c$$$                  CALL CTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE,
c$$$     $                        U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ),
c$$$     $                        LDWV )
c$$$*
c$$$*                 ==== Multiply by U11 ====
c$$$*
c$$$                  CALL CGEMM( 'N', 'N', JLEN, I2, J2, ONE,
c$$$     $                        H( JROW, INCOL+1 ), LDH, U, LDU, ONE, WV,
c$$$     $                        LDWV )
c$$$*
c$$$*                 ==== Copy left of H to right of scratch ====
c$$$*
c$$$                  CALL CLACPY( 'ALL', JLEN, J2, H( JROW, INCOL+1 ), LDH,
c$$$     $                         WV( 1, 1+I2 ), LDWV )
c$$$*
c$$$*                 ==== Multiply by U21 ====
c$$$*
c$$$                  CALL CTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE,
c$$$     $                        U( 1, I2+1 ), LDU, WV( 1, 1+I2 ), LDWV )
c$$$*
c$$$*                 ==== Multiply by U22 ====
c$$$*
c$$$                  CALL CGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE,
c$$$     $                        H( JROW, INCOL+1+J2 ), LDH,
c$$$     $                        U( J2+1, I2+1 ), LDU, ONE, WV( 1, 1+I2 ),
c$$$     $                        LDWV )
c$$$*
c$$$*                 ==== Copy it back ====
c$$$*
c$$$                  CALL CLACPY( 'ALL', JLEN, KDU, WV, LDWV,
c$$$     $                         H( JROW, INCOL+1 ), LDH )
c$$$  190          CONTINUE
c$$$*
c$$$*              ==== Multiply Z (also vertical) ====
c$$$*
c$$$               IF( WANTZ ) THEN
c$$$                  DO 200 JROW = ILOZ, IHIZ, NV
c$$$                     JLEN = MIN( NV, IHIZ-JROW+1 )
c$$$*
c$$$*                    ==== Copy right of Z to left of scratch (first
c$$$*                    .     KZS columns get multiplied by zero) ====
c$$$*
c$$$                     CALL CLACPY( 'ALL', JLEN, KNZ,
c$$$     $                            Z( JROW, INCOL+1+J2 ), LDZ,
c$$$     $                            WV( 1, 1+KZS ), LDWV )
c$$$*
c$$$*                    ==== Multiply by U12 ====
c$$$*
c$$$                     CALL CLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV,
c$$$     $                            LDWV )
c$$$                     CALL CTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE,
c$$$     $                           U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ),
c$$$     $                           LDWV )
c$$$*
c$$$*                    ==== Multiply by U11 ====
c$$$*
c$$$                     CALL CGEMM( 'N', 'N', JLEN, I2, J2, ONE,
c$$$     $                           Z( JROW, INCOL+1 ), LDZ, U, LDU, ONE,
c$$$     $                           WV, LDWV )
c$$$*
c$$$*                    ==== Copy left of Z to right of scratch ====
c$$$*
c$$$                     CALL CLACPY( 'ALL', JLEN, J2, Z( JROW, INCOL+1 ),
c$$$     $                            LDZ, WV( 1, 1+I2 ), LDWV )
c$$$*
c$$$*                    ==== Multiply by U21 ====
c$$$*
c$$$                     CALL CTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE,
c$$$     $                           U( 1, I2+1 ), LDU, WV( 1, 1+I2 ),
c$$$     $                           LDWV )
c$$$*
c$$$*                    ==== Multiply by U22 ====
c$$$*
c$$$                     CALL CGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE,
c$$$     $                           Z( JROW, INCOL+1+J2 ), LDZ,
c$$$     $                           U( J2+1, I2+1 ), LDU, ONE,
c$$$     $                           WV( 1, 1+I2 ), LDWV )
c$$$*
c$$$*                    ==== Copy the result back to Z ====
c$$$*
c$$$                     CALL CLACPY( 'ALL', JLEN, KDU, WV, LDWV,
c$$$     $                            Z( JROW, INCOL+1 ), LDZ )
c$$$  200             CONTINUE
c$$$               END IF
c$$$            END IF
c$$$         END IF
c$$$  210 CONTINUE
c$$$*
c$$$*     ==== End of CLAQR5 ====
c$$$*
c$$$      END
c$$$      SUBROUTINE CLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
c$$$     $                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
c$$$     $                   NV, WV, LDWV, WORK, LWORK )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
c$$$     $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
c$$$      LOGICAL            WANTT, WANTZ
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
c$$$     $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
c$$$*     ..
c$$$*
c$$$*     This subroutine is identical to CLAQR3 except that it avoids
c$$$*     recursion by calling CLAHQR instead of CLAQR4.
c$$$*
c$$$*
c$$$*     ******************************************************************
c$$$*     Aggressive early deflation:
c$$$*
c$$$*     This subroutine accepts as input an upper Hessenberg matrix
c$$$*     H and performs an unitary similarity transformation
c$$$*     designed to detect and deflate fully converged eigenvalues from
c$$$*     a trailing principal submatrix.  On output H has been over-
c$$$*     written by a new Hessenberg matrix that is a perturbation of
c$$$*     an unitary similarity transformation of H.  It is to be
c$$$*     hoped that the final version of H has many zero subdiagonal
c$$$*     entries.
c$$$*
c$$$*     ******************************************************************
c$$$*     WANTT   (input) LOGICAL
c$$$*          If .TRUE., then the Hessenberg matrix H is fully updated
c$$$*          so that the triangular Schur factor may be
c$$$*          computed (in cooperation with the calling subroutine).
c$$$*          If .FALSE., then only enough of H is updated to preserve
c$$$*          the eigenvalues.
c$$$*
c$$$*     WANTZ   (input) LOGICAL
c$$$*          If .TRUE., then the unitary matrix Z is updated so
c$$$*          so that the unitary Schur factor may be computed
c$$$*          (in cooperation with the calling subroutine).
c$$$*          If .FALSE., then Z is not referenced.
c$$$*
c$$$*     N       (input) INTEGER
c$$$*          The order of the matrix H and (if WANTZ is .TRUE.) the
c$$$*          order of the unitary matrix Z.
c$$$*
c$$$*     KTOP    (input) INTEGER
c$$$*          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
c$$$*          KBOT and KTOP together determine an isolated block
c$$$*          along the diagonal of the Hessenberg matrix.
c$$$*
c$$$*     KBOT    (input) INTEGER
c$$$*          It is assumed without a check that either
c$$$*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
c$$$*          determine an isolated block along the diagonal of the
c$$$*          Hessenberg matrix.
c$$$*
c$$$*     NW      (input) INTEGER
c$$$*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
c$$$*
c$$$*     H       (input/output) COMPLEX array, dimension (LDH,N)
c$$$*          On input the initial N-by-N section of H stores the
c$$$*          Hessenberg matrix undergoing aggressive early deflation.
c$$$*          On output H has been transformed by a unitary
c$$$*          similarity transformation, perturbed, and the returned
c$$$*          to Hessenberg form that (it is to be hoped) has some
c$$$*          zero subdiagonal entries.
c$$$*
c$$$*     LDH     (input) integer
c$$$*          Leading dimension of H just as declared in the calling
c$$$*          subroutine.  N .LE. LDH
c$$$*
c$$$*     ILOZ    (input) INTEGER
c$$$*     IHIZ    (input) INTEGER
c$$$*          Specify the rows of Z to which transformations must be
c$$$*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
c$$$*
c$$$*     Z       (input/output) COMPLEX array, dimension (LDZ,IHI)
c$$$*          IF WANTZ is .TRUE., then on output, the unitary
c$$$*          similarity transformation mentioned above has been
c$$$*          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
c$$$*          If WANTZ is .FALSE., then Z is unreferenced.
c$$$*
c$$$*     LDZ     (input) integer
c$$$*          The leading dimension of Z just as declared in the
c$$$*          calling subroutine.  1 .LE. LDZ.
c$$$*
c$$$*     NS      (output) integer
c$$$*          The number of unconverged (ie approximate) eigenvalues
c$$$*          returned in SR and SI that may be used as shifts by the
c$$$*          calling subroutine.
c$$$*
c$$$*     ND      (output) integer
c$$$*          The number of converged eigenvalues uncovered by this
c$$$*          subroutine.
c$$$*
c$$$*     SH      (output) COMPLEX array, dimension KBOT
c$$$*          On output, approximate eigenvalues that may
c$$$*          be used for shifts are stored in SH(KBOT-ND-NS+1)
c$$$*          through SR(KBOT-ND).  Converged eigenvalues are
c$$$*          stored in SH(KBOT-ND+1) through SH(KBOT).
c$$$*
c$$$*     V       (workspace) COMPLEX array, dimension (LDV,NW)
c$$$*          An NW-by-NW work array.
c$$$*
c$$$*     LDV     (input) integer scalar
c$$$*          The leading dimension of V just as declared in the
c$$$*          calling subroutine.  NW .LE. LDV
c$$$*
c$$$*     NH      (input) integer scalar
c$$$*          The number of columns of T.  NH.GE.NW.
c$$$*
c$$$*     T       (workspace) COMPLEX array, dimension (LDT,NW)
c$$$*
c$$$*     LDT     (input) integer
c$$$*          The leading dimension of T just as declared in the
c$$$*          calling subroutine.  NW .LE. LDT
c$$$*
c$$$*     NV      (input) integer
c$$$*          The number of rows of work array WV available for
c$$$*          workspace.  NV.GE.NW.
c$$$*
c$$$*     WV      (workspace) COMPLEX array, dimension (LDWV,NW)
c$$$*
c$$$*     LDWV    (input) integer
c$$$*          The leading dimension of W just as declared in the
c$$$*          calling subroutine.  NW .LE. LDV
c$$$*
c$$$*     WORK    (workspace) COMPLEX array, dimension LWORK.
c$$$*          On exit, WORK(1) is set to an estimate of the optimal value
c$$$*          of LWORK for the given values of N, NW, KTOP and KBOT.
c$$$*
c$$$*     LWORK   (input) integer
c$$$*          The dimension of the work array WORK.  LWORK = 2*NW
c$$$*          suffices, but greater efficiency may result from larger
c$$$*          values of LWORK.
c$$$*
c$$$*          If LWORK = -1, then a workspace query is assumed; CLAQR2
c$$$*          only estimates the optimal workspace size for the given
c$$$*          values of N, NW, KTOP and KBOT.  The estimate is returned
c$$$*          in WORK(1).  No error message related to LWORK is issued
c$$$*          by XERBLA.  Neither H nor Z are accessed.
c$$$*
c$$$*     ================================================================
c$$$*     Based on contributions by
c$$$*        Karen Braman and Ralph Byers, Department of Mathematics,
c$$$*        University of Kansas, USA
c$$$*
c$$$*     ==================================================================
c$$$*     .. Parameters ..
c$$$      COMPLEX            ZERO, ONE
c$$$      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),
c$$$     $                   ONE = ( 1.0e0, 0.0e0 ) )
c$$$      REAL               RZERO, RONE
c$$$      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0 )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      COMPLEX            BETA, CDUM, S, TAU
c$$$      REAL               FOO, SAFMAX, SAFMIN, SMLNUM, ULP
c$$$      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN,
c$$$     $                   KNT, KROW, KWTOP, LTOP, LWK1, LWK2, LWKOPT
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      REAL               SLAMCH
c$$$      EXTERNAL           SLAMCH
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CCOPY, CGEHRD, CGEMM, CLACPY, CLAHQR, CLARF,
c$$$     $                   CLARFG, CLASET, CTREXC, CUNGHR, SLABAD
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          ABS, AIMAG, CMPLX, CONJG, INT, MAX, MIN, REAL
c$$$*     ..
c$$$*     .. Statement Functions ..
c$$$      REAL               CABS1
c$$$*     ..
c$$$*     .. Statement Function definitions ..
c$$$      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     ==== Estimate optimal workspace. ====
c$$$*
c$$$      JW = MIN( NW, KBOT-KTOP+1 )
c$$$      IF( JW.LE.2 ) THEN
c$$$         LWKOPT = 1
c$$$      ELSE
c$$$*
c$$$*        ==== Workspace query call to CGEHRD ====
c$$$*
c$$$         CALL CGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
c$$$         LWK1 = INT( WORK( 1 ) )
c$$$*
c$$$*        ==== Workspace query call to CUNGHR ====
c$$$*
c$$$         CALL CUNGHR( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
c$$$         LWK2 = INT( WORK( 1 ) )
c$$$*
c$$$*        ==== Optimal workspace ====
c$$$*
c$$$         LWKOPT = JW + MAX( LWK1, LWK2 )
c$$$      END IF
c$$$*
c$$$*     ==== Quick return in case of workspace query. ====
c$$$*
c$$$      IF( LWORK.EQ.-1 ) THEN
c$$$         WORK( 1 ) = CMPLX( LWKOPT, 0 )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     ==== Nothing to do ...
c$$$*     ... for an empty active block ... ====
c$$$      NS = 0
c$$$      ND = 0
c$$$      IF( KTOP.GT.KBOT )
c$$$     $   RETURN
c$$$*     ... nor for an empty deflation window. ====
c$$$      IF( NW.LT.1 )
c$$$     $   RETURN
c$$$*
c$$$*     ==== Machine constants ====
c$$$*
c$$$      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
c$$$      SAFMAX = RONE / SAFMIN
c$$$      CALL SLABAD( SAFMIN, SAFMAX )
c$$$      ULP = SLAMCH( 'PRECISION' )
c$$$      SMLNUM = SAFMIN*( REAL( N ) / ULP )
c$$$*
c$$$*     ==== Setup deflation window ====
c$$$*
c$$$      JW = MIN( NW, KBOT-KTOP+1 )
c$$$      KWTOP = KBOT - JW + 1
c$$$      IF( KWTOP.EQ.KTOP ) THEN
c$$$         S = ZERO
c$$$      ELSE
c$$$         S = H( KWTOP, KWTOP-1 )
c$$$      END IF
c$$$*
c$$$      IF( KBOT.EQ.KWTOP ) THEN
c$$$*
c$$$*        ==== 1-by-1 deflation window: not much to do ====
c$$$*
c$$$         SH( KWTOP ) = H( KWTOP, KWTOP )
c$$$         NS = 1
c$$$         ND = 0
c$$$         IF( CABS1( S ).LE.MAX( SMLNUM, ULP*CABS1( H( KWTOP,
c$$$     $       KWTOP ) ) ) ) THEN
c$$$
c$$$            NS = 0
c$$$            ND = 1
c$$$            IF( KWTOP.GT.KTOP )
c$$$     $         H( KWTOP, KWTOP-1 ) = ZERO
c$$$         END IF
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     ==== Convert to spike-triangular form.  (In case of a
c$$$*     .    rare QR failure, this routine continues to do
c$$$*     .    aggressive early deflation using that part of
c$$$*     .    the deflation window that converged using INFQR
c$$$*     .    here and there to keep track.) ====
c$$$*
c$$$      CALL CLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
c$$$      CALL CCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
c$$$*
c$$$      CALL CLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
c$$$      CALL CLAHQR( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1,
c$$$     $             JW, V, LDV, INFQR )
c$$$*
c$$$*     ==== Deflation detection loop ====
c$$$*
c$$$      NS = JW
c$$$      ILST = INFQR + 1
c$$$      DO 10 KNT = INFQR + 1, JW
c$$$*
c$$$*        ==== Small spike tip deflation test ====
c$$$*
c$$$         FOO = CABS1( T( NS, NS ) )
c$$$         IF( FOO.EQ.RZERO )
c$$$     $      FOO = CABS1( S )
c$$$         IF( CABS1( S )*CABS1( V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) )
c$$$     $        THEN
c$$$*
c$$$*           ==== One more converged eigenvalue ====
c$$$*
c$$$            NS = NS - 1
c$$$         ELSE
c$$$*
c$$$*           ==== One undflatable eigenvalue.  Move it up out of the
c$$$*           .    way.   (CTREXC can not fail in this case.) ====
c$$$*
c$$$            IFST = NS
c$$$            CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
c$$$            ILST = ILST + 1
c$$$         END IF
c$$$   10 CONTINUE
c$$$*
c$$$*        ==== Return to Hessenberg form ====
c$$$*
c$$$      IF( NS.EQ.0 )
c$$$     $   S = ZERO
c$$$*
c$$$      IF( NS.LT.JW ) THEN
c$$$*
c$$$*        ==== sorting the diagonal of T improves accuracy for
c$$$*        .    graded matrices.  ====
c$$$*
c$$$         DO 30 I = INFQR + 1, NS
c$$$            IFST = I
c$$$            DO 20 J = I + 1, NS
c$$$               IF( CABS1( T( J, J ) ).GT.CABS1( T( IFST, IFST ) ) )
c$$$     $            IFST = J
c$$$   20       CONTINUE
c$$$            ILST = I
c$$$            IF( IFST.NE.ILST )
c$$$     $         CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
c$$$   30    CONTINUE
c$$$      END IF
c$$$*
c$$$*     ==== Restore shift/eigenvalue array from T ====
c$$$*
c$$$      DO 40 I = INFQR + 1, JW
c$$$         SH( KWTOP+I-1 ) = T( I, I )
c$$$   40 CONTINUE
c$$$*
c$$$*
c$$$      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
c$$$         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
c$$$*
c$$$*           ==== Reflect spike back into lower triangle ====
c$$$*
c$$$            CALL CCOPY( NS, V, LDV, WORK, 1 )
c$$$            DO 50 I = 1, NS
c$$$               WORK( I ) = CONJG( WORK( I ) )
c$$$   50       CONTINUE
c$$$            BETA = WORK( 1 )
c$$$            CALL CLARFG( NS, BETA, WORK( 2 ), 1, TAU )
c$$$            WORK( 1 ) = ONE
c$$$*
c$$$            CALL CLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
c$$$*
c$$$            CALL CLARF( 'L', NS, JW, WORK, 1, CONJG( TAU ), T, LDT,
c$$$     $                  WORK( JW+1 ) )
c$$$            CALL CLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT,
c$$$     $                  WORK( JW+1 ) )
c$$$            CALL CLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV,
c$$$     $                  WORK( JW+1 ) )
c$$$*
c$$$            CALL CGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),
c$$$     $                   LWORK-JW, INFO )
c$$$         END IF
c$$$*
c$$$*        ==== Copy updated reduced window into place ====
c$$$*
c$$$         IF( KWTOP.GT.1 )
c$$$     $      H( KWTOP, KWTOP-1 ) = S*CONJG( V( 1, 1 ) )
c$$$         CALL CLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
c$$$         CALL CCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ),
c$$$     $               LDH+1 )
c$$$*
c$$$*        ==== Accumulate orthogonal matrix in order update
c$$$*        .    H and Z, if requested.  (A modified version
c$$$*        .    of  CUNGHR that accumulates block Householder
c$$$*        .    transformations into V directly might be
c$$$*        .    marginally more efficient than the following.) ====
c$$$*
c$$$         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
c$$$            CALL CUNGHR( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),
c$$$     $                   LWORK-JW, INFO )
c$$$            CALL CGEMM( 'N', 'N', JW, NS, NS, ONE, V, LDV, T, LDT, ZERO,
c$$$     $                  WV, LDWV )
c$$$            CALL CLACPY( 'A', JW, NS, WV, LDWV, V, LDV )
c$$$         END IF
c$$$*
c$$$*        ==== Update vertical slab in H ====
c$$$*
c$$$         IF( WANTT ) THEN
c$$$            LTOP = 1
c$$$         ELSE
c$$$            LTOP = KTOP
c$$$         END IF
c$$$         DO 60 KROW = LTOP, KWTOP - 1, NV
c$$$            KLN = MIN( NV, KWTOP-KROW )
c$$$            CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ),
c$$$     $                  LDH, V, LDV, ZERO, WV, LDWV )
c$$$            CALL CLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
c$$$   60    CONTINUE
c$$$*
c$$$*        ==== Update horizontal slab in H ====
c$$$*
c$$$         IF( WANTT ) THEN
c$$$            DO 70 KCOL = KBOT + 1, N, NH
c$$$               KLN = MIN( NH, N-KCOL+1 )
c$$$               CALL CGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV,
c$$$     $                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
c$$$               CALL CLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ),
c$$$     $                      LDH )
c$$$   70       CONTINUE
c$$$         END IF
c$$$*
c$$$*        ==== Update vertical slab in Z ====
c$$$*
c$$$         IF( WANTZ ) THEN
c$$$            DO 80 KROW = ILOZ, IHIZ, NV
c$$$               KLN = MIN( NV, IHIZ-KROW+1 )
c$$$               CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ),
c$$$     $                     LDZ, V, LDV, ZERO, WV, LDWV )
c$$$               CALL CLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ),
c$$$     $                      LDZ )
c$$$   80       CONTINUE
c$$$         END IF
c$$$      END IF
c$$$*
c$$$*     ==== Return the number of deflations ... ====
c$$$*
c$$$      ND = JW - NS
c$$$*
c$$$*     ==== ... and the number of shifts. (Subtracting
c$$$*     .    INFQR from the spike length takes care
c$$$*     .    of the case of a rare QR failure while
c$$$*     .    calculating eigenvalues of the deflation
c$$$*     .    window.)  ====
c$$$*
c$$$      NS = NS - INFQR
c$$$*
c$$$*      ==== Return optimal workspace. ====
c$$$*
c$$$      WORK( 1 ) = CMPLX( LWKOPT, 0 )
c$$$*
c$$$*     ==== End of CLAQR2 ====
c$$$*
c$$$      END
c$$$      SUBROUTINE CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
c$$$     $                   IHIZ, Z, LDZ, INFO )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
c$$$      LOGICAL            WANTT, WANTZ
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            H( LDH, * ), W( * ), Z( LDZ, * )
c$$$*     ..
c$$$*
c$$$*     Purpose
c$$$*     =======
c$$$*
c$$$*     CLAHQR is an auxiliary routine called by CHSEQR to update the
c$$$*     eigenvalues and Schur decomposition already computed by CHSEQR, by
c$$$*     dealing with the Hessenberg submatrix in rows and columns ILO to
c$$$*     IHI.
c$$$*
c$$$*     Arguments
c$$$*     =========
c$$$*
c$$$*     WANTT   (input) LOGICAL
c$$$*          = .TRUE. : the full Schur form T is required;
c$$$*          = .FALSE.: only eigenvalues are required.
c$$$*
c$$$*     WANTZ   (input) LOGICAL
c$$$*          = .TRUE. : the matrix of Schur vectors Z is required;
c$$$*          = .FALSE.: Schur vectors are not required.
c$$$*
c$$$*     N       (input) INTEGER
c$$$*          The order of the matrix H.  N >= 0.
c$$$*
c$$$*     ILO     (input) INTEGER
c$$$*     IHI     (input) INTEGER
c$$$*          It is assumed that H is already upper triangular in rows and
c$$$*          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
c$$$*          CLAHQR works primarily with the Hessenberg submatrix in rows
c$$$*          and columns ILO to IHI, but applies transformations to all of
c$$$*          H if WANTT is .TRUE..
c$$$*          1 <= ILO <= max(1,IHI); IHI <= N.
c$$$*
c$$$*     H       (input/output) COMPLEX array, dimension (LDH,N)
c$$$*          On entry, the upper Hessenberg matrix H.
c$$$*          On exit, if INFO is zero and if WANTT is .TRUE., then H
c$$$*          is upper triangular in rows and columns ILO:IHI.  If INFO
c$$$*          is zero and if WANTT is .FALSE., then the contents of H
c$$$*          are unspecified on exit.  The output state of H in case
c$$$*          INF is positive is below under the description of INFO.
c$$$*
c$$$*     LDH     (input) INTEGER
c$$$*          The leading dimension of the array H. LDH >= max(1,N).
c$$$*
c$$$*     W       (output) COMPLEX array, dimension (N)
c$$$*          The computed eigenvalues ILO to IHI are stored in the
c$$$*          corresponding elements of W. If WANTT is .TRUE., the
c$$$*          eigenvalues are stored in the same order as on the diagonal
c$$$*          of the Schur form returned in H, with W(i) = H(i,i).
c$$$*
c$$$*     ILOZ    (input) INTEGER
c$$$*     IHIZ    (input) INTEGER
c$$$*          Specify the rows of Z to which transformations must be
c$$$*          applied if WANTZ is .TRUE..
c$$$*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
c$$$*
c$$$*     Z       (input/output) COMPLEX array, dimension (LDZ,N)
c$$$*          If WANTZ is .TRUE., on entry Z must contain the current
c$$$*          matrix Z of transformations accumulated by CHSEQR, and on
c$$$*          exit Z has been updated; transformations are applied only to
c$$$*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
c$$$*          If WANTZ is .FALSE., Z is not referenced.
c$$$*
c$$$*     LDZ     (input) INTEGER
c$$$*          The leading dimension of the array Z. LDZ >= max(1,N).
c$$$*
c$$$*     INFO    (output) INTEGER
c$$$*           =   0: successful exit
c$$$*          .GT. 0: if INFO = i, CLAHQR failed to compute all the
c$$$*                  eigenvalues ILO to IHI in a total of 30 iterations
c$$$*                  per eigenvalue; elements i+1:ihi of W contain
c$$$*                  those eigenvalues which have been successfully
c$$$*                  computed.
c$$$*
c$$$*                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,
c$$$*                  the remaining unconverged eigenvalues are the
c$$$*                  eigenvalues of the upper Hessenberg matrix
c$$$*                  rows and columns ILO thorugh INFO of the final,
c$$$*                  output value of H.
c$$$*
c$$$*                  If INFO .GT. 0 and WANTT is .TRUE., then on exit
c$$$*          (*)       (initial value of H)*U  = U*(final value of H)
c$$$*                  where U is an orthognal matrix.    The final
c$$$*                  value of H is upper Hessenberg and triangular in
c$$$*                  rows and columns INFO+1 through IHI.
c$$$*
c$$$*                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit
c$$$*                      (final value of Z)  = (initial value of Z)*U
c$$$*                  where U is the orthogonal matrix in (*)
c$$$*                  (regardless of the value of WANTT.)
c$$$*
c$$$*     Further Details
c$$$*     ===============
c$$$*
c$$$*     02-96 Based on modifications by
c$$$*     David Day, Sandia National Laboratory, USA
c$$$*
c$$$*     12-04 Further modifications by
c$$$*     Ralph Byers, University of Kansas, USA
c$$$*
c$$$*       This is a modified version of CLAHQR from LAPACK version 3.0.
c$$$*       It is (1) more robust against overflow and underflow and
c$$$*       (2) adopts the more conservative Ahues & Tisseur stopping
c$$$*       criterion (LAWN 122, 1997).
c$$$*
c$$$*     =========================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      INTEGER            ITMAX
c$$$      PARAMETER          ( ITMAX = 30 )
c$$$      COMPLEX            ZERO, ONE
c$$$      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ),
c$$$     $                   ONE = ( 1.0e0, 0.0e0 ) )
c$$$      REAL               RZERO, RONE, HALF
c$$$      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0, HALF = 0.5e0 )
c$$$      REAL               DAT1
c$$$      PARAMETER          ( DAT1 = 3.0e0 / 4.0e0 )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      COMPLEX            CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U,
c$$$     $                   V2, X, Y
c$$$      REAL               AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX,
c$$$     $                   SAFMIN, SMLNUM, SX, T2, TST, ULP
c$$$      INTEGER            I, I1, I2, ITS, J, JHI, JLO, K, L, M, NH, NZ
c$$$*     ..
c$$$*     .. Local Arrays ..
c$$$      COMPLEX            V( 2 )
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      COMPLEX            CLADIV
c$$$      REAL               SLAMCH
c$$$      EXTERNAL           CLADIV, SLAMCH
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CCOPY, CLARFG, CSCAL, SLABAD
c$$$*     ..
c$$$*     .. Statement Functions ..
c$$$      REAL               CABS1
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          ABS, AIMAG, CONJG, MAX, MIN, REAL, SQRT
c$$$*     ..
c$$$*     .. Statement Function definitions ..
c$$$      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$      INFO = 0
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      IF( N.EQ.0 )
c$$$     $   RETURN
c$$$      IF( ILO.EQ.IHI ) THEN
c$$$         W( ILO ) = H( ILO, ILO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     ==== clear out the trash ====
c$$$      DO 10 J = ILO, IHI - 3
c$$$         H( J+2, J ) = ZERO
c$$$         H( J+3, J ) = ZERO
c$$$   10 CONTINUE
c$$$      IF( ILO.LE.IHI-2 )
c$$$     $   H( IHI, IHI-2 ) = ZERO
c$$$*     ==== ensure that subdiagonal entries are real ====
c$$$      DO 20 I = ILO + 1, IHI
c$$$         IF( AIMAG( H( I, I-1 ) ).NE.RZERO ) THEN
c$$$*           ==== The following redundant normalization
c$$$*           .    avoids problems with both gradual and
c$$$*           .    sudden underflow in ABS(H(I,I-1)) ====
c$$$            SC = H( I, I-1 ) / CABS1( H( I, I-1 ) )
c$$$            SC = CONJG( SC ) / ABS( SC )
c$$$            H( I, I-1 ) = ABS( H( I, I-1 ) )
c$$$            IF( WANTT ) THEN
c$$$               JLO = 1
c$$$               JHI = N
c$$$            ELSE
c$$$               JLO = ILO
c$$$               JHI = IHI
c$$$            END IF
c$$$            CALL CSCAL( JHI-I+1, SC, H( I, I ), LDH )
c$$$            CALL CSCAL( MIN( JHI, I+1 )-JLO+1, CONJG( SC ), H( JLO, I ),
c$$$     $                  1 )
c$$$            IF( WANTZ )
c$$$     $         CALL CSCAL( IHIZ-ILOZ+1, CONJG( SC ), Z( ILOZ, I ), 1 )
c$$$         END IF
c$$$   20 CONTINUE
c$$$*
c$$$      NH = IHI - ILO + 1
c$$$      NZ = IHIZ - ILOZ + 1
c$$$*
c$$$*     Set machine-dependent constants for the stopping criterion.
c$$$*
c$$$      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
c$$$      SAFMAX = RONE / SAFMIN
c$$$      CALL SLABAD( SAFMIN, SAFMAX )
c$$$      ULP = SLAMCH( 'PRECISION' )
c$$$      SMLNUM = SAFMIN*( REAL( NH ) / ULP )
c$$$*
c$$$*     I1 and I2 are the indices of the first row and last column of H
c$$$*     to which transformations must be applied. If eigenvalues only are
c$$$*     being computed, I1 and I2 are set inside the main loop.
c$$$*
c$$$      IF( WANTT ) THEN
c$$$         I1 = 1
c$$$         I2 = N
c$$$      END IF
c$$$*
c$$$*     The main loop begins here. I is the loop index and decreases from
c$$$*     IHI to ILO in steps of 1. Each iteration of the loop works
c$$$*     with the active submatrix in rows and columns L to I.
c$$$*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
c$$$*     H(L,L-1) is negligible so that the matrix splits.
c$$$*
c$$$      I = IHI
c$$$   30 CONTINUE
c$$$      IF( I.LT.ILO )
c$$$     $   GO TO 150
c$$$*
c$$$*     Perform QR iterations on rows and columns ILO to I until a
c$$$*     submatrix of order 1 splits off at the bottom because a
c$$$*     subdiagonal element has become negligible.
c$$$*
c$$$      L = ILO
c$$$      DO 130 ITS = 0, ITMAX
c$$$*
c$$$*        Look for a single small subdiagonal element.
c$$$*
c$$$         DO 40 K = I, L + 1, -1
c$$$            IF( CABS1( H( K, K-1 ) ).LE.SMLNUM )
c$$$     $         GO TO 50
c$$$            TST = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
c$$$            IF( TST.EQ.ZERO ) THEN
c$$$               IF( K-2.GE.ILO )
c$$$     $            TST = TST + ABS( REAL( H( K-1, K-2 ) ) )
c$$$               IF( K+1.LE.IHI )
c$$$     $            TST = TST + ABS( REAL( H( K+1, K ) ) )
c$$$            END IF
c$$$*           ==== The following is a conservative small subdiagonal
c$$$*           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
c$$$*           .    1997). It has better mathematical foundation and
c$$$*           .    improves accuracy in some examples.  ====
c$$$            IF( ABS( REAL( H( K, K-1 ) ) ).LE.ULP*TST ) THEN
c$$$               AB = MAX( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
c$$$               BA = MIN( CABS1( H( K, K-1 ) ), CABS1( H( K-1, K ) ) )
c$$$               AA = MAX( CABS1( H( K, K ) ),
c$$$     $              CABS1( H( K-1, K-1 )-H( K, K ) ) )
c$$$               BB = MIN( CABS1( H( K, K ) ),
c$$$     $              CABS1( H( K-1, K-1 )-H( K, K ) ) )
c$$$               S = AA + AB
c$$$               IF( BA*( AB / S ).LE.MAX( SMLNUM,
c$$$     $             ULP*( BB*( AA / S ) ) ) )GO TO 50
c$$$            END IF
c$$$   40    CONTINUE
c$$$   50    CONTINUE
c$$$         L = K
c$$$         IF( L.GT.ILO ) THEN
c$$$*
c$$$*           H(L,L-1) is negligible
c$$$*
c$$$            H( L, L-1 ) = ZERO
c$$$         END IF
c$$$*
c$$$*        Exit from loop if a submatrix of order 1 has split off.
c$$$*
c$$$         IF( L.GE.I )
c$$$     $      GO TO 140
c$$$*
c$$$*        Now the active submatrix is in rows and columns L to I. If
c$$$*        eigenvalues only are being computed, only the active submatrix
c$$$*        need be transformed.
c$$$*
c$$$         IF( .NOT.WANTT ) THEN
c$$$            I1 = L
c$$$            I2 = I
c$$$         END IF
c$$$*
c$$$         IF( ITS.EQ.10 .OR. ITS.EQ.20 ) THEN
c$$$*
c$$$*           Exceptional shift.
c$$$*
c$$$            S = DAT1*ABS( REAL( H( I, I-1 ) ) )
c$$$            T = S + H( I, I )
c$$$         ELSE
c$$$*
c$$$*           Wilkinson's shift.
c$$$*
c$$$            T = H( I, I )
c$$$            U = SQRT( H( I-1, I ) )*SQRT( H( I, I-1 ) )
c$$$            S = CABS1( U )
c$$$            IF( S.NE.RZERO ) THEN
c$$$               X = HALF*( H( I-1, I-1 )-T )
c$$$               SX = CABS1( X )
c$$$               S = MAX( S, CABS1( X ) )
c$$$               Y = S*SQRT( ( X / S )**2+( U / S )**2 )
c$$$               IF( SX.GT.RZERO ) THEN
c$$$                  IF( REAL( X / SX )*REAL( Y )+AIMAG( X / SX )*
c$$$     $                AIMAG( Y ).LT.RZERO )Y = -Y
c$$$               END IF
c$$$               T = T - U*CLADIV( U, ( X+Y ) )
c$$$            END IF
c$$$         END IF
c$$$*
c$$$*        Look for two consecutive small subdiagonal elements.
c$$$*
c$$$         DO 60 M = I - 1, L + 1, -1
c$$$*
c$$$*           Determine the effect of starting the single-shift QR
c$$$*           iteration at row M, and see if this would make H(M,M-1)
c$$$*           negligible.
c$$$*
c$$$            H11 = H( M, M )
c$$$            H22 = H( M+1, M+1 )
c$$$            H11S = H11 - T
c$$$            H21 = H( M+1, M )
c$$$            S = CABS1( H11S ) + ABS( H21 )
c$$$            H11S = H11S / S
c$$$            H21 = H21 / S
c$$$            V( 1 ) = H11S
c$$$            V( 2 ) = H21
c$$$            H10 = H( M, M-1 )
c$$$            IF( ABS( H10 )*ABS( H21 ).LE.ULP*
c$$$     $          ( CABS1( H11S )*( CABS1( H11 )+CABS1( H22 ) ) ) )
c$$$     $          GO TO 70
c$$$   60    CONTINUE
c$$$         H11 = H( L, L )
c$$$         H22 = H( L+1, L+1 )
c$$$         H11S = H11 - T
c$$$         H21 = H( L+1, L )
c$$$         S = CABS1( H11S ) + ABS( H21 )
c$$$         H11S = H11S / S
c$$$         H21 = H21 / S
c$$$         V( 1 ) = H11S
c$$$         V( 2 ) = H21
c$$$   70    CONTINUE
c$$$*
c$$$*        Single-shift QR step
c$$$*
c$$$         DO 120 K = M, I - 1
c$$$*
c$$$*           The first iteration of this loop determines a reflection G
c$$$*           from the vector V and applies it from left and right to H,
c$$$*           thus creating a nonzero bulge below the subdiagonal.
c$$$*
c$$$*           Each subsequent iteration determines a reflection G to
c$$$*           restore the Hessenberg form in the (K-1)th column, and thus
c$$$*           chases the bulge one step toward the bottom of the active
c$$$*           submatrix.
c$$$*
c$$$*           V(2) is always real before the call to CLARFG, and hence
c$$$*           after the call T2 ( = T1*V(2) ) is also real.
c$$$*
c$$$            IF( K.GT.M )
c$$$     $         CALL CCOPY( 2, H( K, K-1 ), 1, V, 1 )
c$$$            CALL CLARFG( 2, V( 1 ), V( 2 ), 1, T1 )
c$$$            IF( K.GT.M ) THEN
c$$$               H( K, K-1 ) = V( 1 )
c$$$               H( K+1, K-1 ) = ZERO
c$$$            END IF
c$$$            V2 = V( 2 )
c$$$            T2 = REAL( T1*V2 )
c$$$*
c$$$*           Apply G from the left to transform the rows of the matrix
c$$$*           in columns K to I2.
c$$$*
c$$$            DO 80 J = K, I2
c$$$               SUM = CONJG( T1 )*H( K, J ) + T2*H( K+1, J )
c$$$               H( K, J ) = H( K, J ) - SUM
c$$$               H( K+1, J ) = H( K+1, J ) - SUM*V2
c$$$   80       CONTINUE
c$$$*
c$$$*           Apply G from the right to transform the columns of the
c$$$*           matrix in rows I1 to min(K+2,I).
c$$$*
c$$$            DO 90 J = I1, MIN( K+2, I )
c$$$               SUM = T1*H( J, K ) + T2*H( J, K+1 )
c$$$               H( J, K ) = H( J, K ) - SUM
c$$$               H( J, K+1 ) = H( J, K+1 ) - SUM*CONJG( V2 )
c$$$   90       CONTINUE
c$$$*
c$$$            IF( WANTZ ) THEN
c$$$*
c$$$*              Accumulate transformations in the matrix Z
c$$$*
c$$$               DO 100 J = ILOZ, IHIZ
c$$$                  SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
c$$$                  Z( J, K ) = Z( J, K ) - SUM
c$$$                  Z( J, K+1 ) = Z( J, K+1 ) - SUM*CONJG( V2 )
c$$$  100          CONTINUE
c$$$            END IF
c$$$*
c$$$            IF( K.EQ.M .AND. M.GT.L ) THEN
c$$$*
c$$$*              If the QR step was started at row M > L because two
c$$$*              consecutive small subdiagonals were found, then extra
c$$$*              scaling must be performed to ensure that H(M,M-1) remains
c$$$*              real.
c$$$*
c$$$               TEMP = ONE - T1
c$$$               TEMP = TEMP / ABS( TEMP )
c$$$               H( M+1, M ) = H( M+1, M )*CONJG( TEMP )
c$$$               IF( M+2.LE.I )
c$$$     $            H( M+2, M+1 ) = H( M+2, M+1 )*TEMP
c$$$               DO 110 J = M, I
c$$$                  IF( J.NE.M+1 ) THEN
c$$$                     IF( I2.GT.J )
c$$$     $                  CALL CSCAL( I2-J, TEMP, H( J, J+1 ), LDH )
c$$$                     CALL CSCAL( J-I1, CONJG( TEMP ), H( I1, J ), 1 )
c$$$                     IF( WANTZ ) THEN
c$$$                        CALL CSCAL( NZ, CONJG( TEMP ), Z( ILOZ, J ), 1 )
c$$$                     END IF
c$$$                  END IF
c$$$  110          CONTINUE
c$$$            END IF
c$$$  120    CONTINUE
c$$$*
c$$$*        Ensure that H(I,I-1) is real.
c$$$*
c$$$         TEMP = H( I, I-1 )
c$$$         IF( AIMAG( TEMP ).NE.RZERO ) THEN
c$$$            RTEMP = ABS( TEMP )
c$$$            H( I, I-1 ) = RTEMP
c$$$            TEMP = TEMP / RTEMP
c$$$            IF( I2.GT.I )
c$$$     $         CALL CSCAL( I2-I, CONJG( TEMP ), H( I, I+1 ), LDH )
c$$$            CALL CSCAL( I-I1, TEMP, H( I1, I ), 1 )
c$$$            IF( WANTZ ) THEN
c$$$               CALL CSCAL( NZ, TEMP, Z( ILOZ, I ), 1 )
c$$$            END IF
c$$$         END IF
c$$$*
c$$$  130 CONTINUE
c$$$*
c$$$*     Failure to converge in remaining number of iterations
c$$$*
c$$$      INFO = I
c$$$      RETURN
c$$$*
c$$$  140 CONTINUE
c$$$*
c$$$*     H(I,I-1) is negligible: one eigenvalue has converged.
c$$$*
c$$$      W( I ) = H( I, I )
c$$$*
c$$$*     return to start of the main loop with new value of I.
c$$$*
c$$$      I = L - 1
c$$$      GO TO 30
c$$$*
c$$$  150 CONTINUE
c$$$      RETURN
c$$$*
c$$$*     End of CLAHQR
c$$$*
c$$$      END
c$$$      SUBROUTINE CLAQR1( N, H, LDH, S1, S2, V )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      COMPLEX            S1, S2
c$$$      INTEGER            LDH, N
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            H( LDH, * ), V( * )
c$$$*     ..
c$$$*
c$$$*       Given a 2-by-2 or 3-by-3 matrix H, CLAQR1 sets v to a
c$$$*       scalar multiple of the first column of the product
c$$$*
c$$$*       (*)  K = (H - s1*I)*(H - s2*I)
c$$$*
c$$$*       scaling to avoid overflows and most underflows.
c$$$*
c$$$*       This is useful for starting double implicit shift bulges
c$$$*       in the QR algorithm.
c$$$*
c$$$*
c$$$*       N      (input) integer
c$$$*              Order of the matrix H. N must be either 2 or 3.
c$$$*
c$$$*       H      (input) COMPLEX array of dimension (LDH,N)
c$$$*              The 2-by-2 or 3-by-3 matrix H in (*).
c$$$*
c$$$*       LDH    (input) integer
c$$$*              The leading dimension of H as declared in
c$$$*              the calling procedure.  LDH.GE.N
c$$$*
c$$$*       S1     (input) COMPLEX
c$$$*       S2     S1 and S2 are the shifts defining K in (*) above.
c$$$*
c$$$*       V      (output) COMPLEX array of dimension N
c$$$*              A scalar multiple of the first column of the
c$$$*              matrix K in (*).
c$$$*
c$$$*     ================================================================
c$$$*     Based on contributions by
c$$$*        Karen Braman and Ralph Byers, Department of Mathematics,
c$$$*        University of Kansas, USA
c$$$*
c$$$*     ================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      COMPLEX            ZERO
c$$$      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ) )
c$$$      REAL               RZERO
c$$$      PARAMETER          ( RZERO = 0.0e0 )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      COMPLEX            CDUM
c$$$      REAL               H21S, H31S, S
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          ABS, AIMAG, REAL
c$$$*     ..
c$$$*     .. Statement Functions ..
c$$$      REAL               CABS1
c$$$*     ..
c$$$*     .. Statement Function definitions ..
c$$$      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$      IF( N.EQ.2 ) THEN
c$$$         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) )
c$$$         IF( S.EQ.RZERO ) THEN
c$$$            V( 1 ) = ZERO
c$$$            V( 2 ) = ZERO
c$$$         ELSE
c$$$            H21S = H( 2, 1 ) / S
c$$$            V( 1 ) = H21S*H( 1, 2 ) + ( H( 1, 1 )-S1 )*
c$$$     $               ( ( H( 1, 1 )-S2 ) / S )
c$$$            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 )
c$$$         END IF
c$$$      ELSE
c$$$         S = CABS1( H( 1, 1 )-S2 ) + CABS1( H( 2, 1 ) ) +
c$$$     $       CABS1( H( 3, 1 ) )
c$$$         IF( S.EQ.ZERO ) THEN
c$$$            V( 1 ) = ZERO
c$$$            V( 2 ) = ZERO
c$$$            V( 3 ) = ZERO
c$$$         ELSE
c$$$            H21S = H( 2, 1 ) / S
c$$$            H31S = H( 3, 1 ) / S
c$$$            V( 1 ) = ( H( 1, 1 )-S1 )*( ( H( 1, 1 )-S2 ) / S ) +
c$$$     $               H( 1, 2 )*H21S + H( 1, 3 )*H31S
c$$$            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-S1-S2 ) + H( 2, 3 )*H31S
c$$$            V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-S1-S2 ) + H21S*H( 3, 2 )
c$$$         END IF
c$$$      END IF
c$$$      END
c$$$      SUBROUTINE CGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
c$$$*
c$$$*  -- LAPACK routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            IHI, ILO, INFO, LDA, N
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGEHD2 reduces a complex general matrix A to upper Hessenberg form H
c$$$*  by a unitary similarity transformation:  Q' * A * Q = H .
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The order of the matrix A.  N >= 0.
c$$$*
c$$$*  ILO     (input) INTEGER
c$$$*  IHI     (input) INTEGER
c$$$*          It is assumed that A is already upper triangular in rows
c$$$*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
c$$$*          set by a previous call to CGEBAL; otherwise they should be
c$$$*          set to 1 and N respectively. See Further Details.
c$$$*          1 <= ILO <= IHI <= max(1,N).
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N)
c$$$*          On entry, the n by n general matrix to be reduced.
c$$$*          On exit, the upper triangle and the first subdiagonal of A
c$$$*          are overwritten with the upper Hessenberg matrix H, and the
c$$$*          elements below the first subdiagonal, with the array TAU,
c$$$*          represent the unitary matrix Q as a product of elementary
c$$$*          reflectors. See Further Details.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  TAU     (output) COMPLEX array, dimension (N-1)
c$$$*          The scalar factors of the elementary reflectors (see Further
c$$$*          Details).
c$$$*
c$$$*  WORK    (workspace) COMPLEX array, dimension (N)
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value.
c$$$*
c$$$*  Further Details
c$$$*  ===============
c$$$*
c$$$*  The matrix Q is represented as a product of (ihi-ilo) elementary
c$$$*  reflectors
c$$$*
c$$$*     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
c$$$*
c$$$*  Each H(i) has the form
c$$$*
c$$$*     H(i) = I - tau * v * v'
c$$$*
c$$$*  where tau is a complex scalar, and v is a complex vector with
c$$$*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
c$$$*  exit in A(i+2:ihi,i), and tau in TAU(i).
c$$$*
c$$$*  The contents of A are illustrated by the following example, with
c$$$*  n = 7, ilo = 2 and ihi = 6:
c$$$*
c$$$*  on entry,                        on exit,
c$$$*
c$$$*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
c$$$*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
c$$$*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
c$$$*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
c$$$*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
c$$$*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
c$$$*  (                         a )    (                          a )
c$$$*
c$$$*  where a denotes an element of the original matrix A, h denotes a
c$$$*  modified element of the upper Hessenberg matrix H, and vi denotes an
c$$$*  element of the vector defining H(i).
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      COMPLEX            ONE
c$$$      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      INTEGER            I
c$$$      COMPLEX            ALPHA
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CLARF, CLARFG, XERBLA
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          CONJG, MAX, MIN
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters
c$$$*
c$$$      INFO = 0
c$$$      IF( N.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
c$$$         INFO = -3
c$$$      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -5
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGEHD2', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$      DO 10 I = ILO, IHI - 1
c$$$*
c$$$*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
c$$$*
c$$$         ALPHA = A( I+1, I )
c$$$         CALL CLARFG( IHI-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
c$$$         A( I+1, I ) = ONE
c$$$*
c$$$*        Apply H(i) to A(1:ihi,i+1:ihi) from the right
c$$$*
c$$$         CALL CLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ),
c$$$     $               A( 1, I+1 ), LDA, WORK )
c$$$*
c$$$*        Apply H(i)' to A(i+1:ihi,i+1:n) from the left
c$$$*
c$$$         CALL CLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1,
c$$$     $               CONJG( TAU( I ) ), A( I+1, I+1 ), LDA, WORK )
c$$$*
c$$$         A( I+1, I ) = ALPHA
c$$$   10 CONTINUE
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of CGEHD2
c$$$*
c$$$      END

      subroutine iterate(gk,kmax,ns,npk,nchtop,vmat,wk,nds,kon,nchan,
     >   maxit,rerr,converged)
      real vmat(nds,nds+1),kon(nchan,nchan),gk(kmax,nchan),vfn(nds)
      integer npk(nchtop+1)
      complex wk(nds)
      real, allocatable :: v(:,:,:)
      logical converged
      
      allocate(v(nds,nchtop,3))

      
c$$$      call makevni(ns,nds,npk,nchtop,vmat,v)
      do nchi = 1, nchtop
         call makevfn(ns,nds,npk,nchtop,vmat,npk(nchi),v(1,nchi,1))
         do nf = 1, nds
c$$$            print*,nf,vfn(nf),v(nf,nchi,1)
            v(nf,nchi,3) = v(nf,nchi,1)
            v(nf,nchi,1) = v(nf,nchi,1) * real(wk(nf))
         enddo
      enddo 

      it = 1
      converged = .false.
      do while (it.lt.maxit.and..not.converged)
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& PRIVATE(nchf,nchi,kf,vfn)
      do nchf = 1, nchtop
         do kf = npk(nchf), npk(nchf+1) - 1
            call makevfn(ns,nds,npk,nchtop,vmat,kf,vfn)
            do nchi = 1, nchtop
               v(kf,nchi,2) = dot_product(vfn(1:nds),v(1:nds,nchi,1))
            enddo
         enddo
      enddo 
C$OMP END PARALLEL DO
      converged = .true.
      do nchi = 1, nchtop
         do nf = 1, nds
            v(nf,nchi,3) = v(nf,nchi,3) + v(nf,nchi,2)
            v(nf,nchi,1) = v(nf,nchi,2) * real(wk(nf))
         enddo
         if (imag(wk(npk(nchi))).eq.0.0) cycle ! closed channel
         do nchf = nchi, nchtop
            if (imag(wk(npk(nchf))).eq.0.0) cycle ! closed channel
            if (.not.converged) cycle
            converged = converged .and. (rerr.gt.
     >         abs(v(npk(nchf),nchi,2)/(v(npk(nchf),nchi,3)+1e-20)))
            print*,it,converged,nchf,nchi,v(npk(nchf),nchi,3),
     >         abs(v(npk(nchf),nchi,2)/(v(npk(nchf),nchi,3)+1e-20))
         enddo 
      enddo
      if (converged) print*,'Converged afer iterations:',it
      it = it + 1
      enddo !it
      print*,'Setting CONVERGED to TRUE'
      converged = .true.

c$$$C  Test the solution of the L-S equation      
c$$$c$$$      call makevni(ns,nds,npk,nchtop,vmat,v)
c$$$      do nchi = 1, nchtop
c$$$         call makevfn(ns,nds,npk,nchtop,vmat,npk(nchi),v(1,nchi,1))
c$$$         do nf = 1, nds
c$$$            v(nf,nchi,2) = v(nf,nchi,3) * real(wk(nf))
c$$$         enddo
c$$$      enddo 
c$$$      
c$$$      do nchf = 1, nchtop
c$$$         do kf = npk(nchf), npk(nchf+1) - 1
c$$$            nchi = 1
c$$$            rkf = gk(kf-npk(nchf)+1,nchf)
c$$$            rki = gk(1,nchi)
c$$$            divk = rkf * rki
c$$$
c$$$            call makevfn(ns,nds,npk,nchtop,vmat,kf,vfn)
c$$$            do nchi = 1, nchtop
c$$$               v(kf,nchi,1) =  v(kf,nchi,1) +
c$$$     >            dot_product(vfn(1:nds),v(1:nds,nchi,2))
c$$$            enddo            
c$$$            if (kf.eq.npk(nchf).and.imag(wk(npk(nchf))).ne.0.0)
c$$$     >         print*,'kf,v(3)/v(1):',kf,v(kf,1,3)/v(kf,1,1)
c$$$         enddo
c$$$      enddo 
      
      do nchi = 1, nchtop
         if (imag(wk(npk(nchi))).eq.0.0) cycle ! closed channel
         do nchf = nchi, nchtop
            if (imag(wk(npk(nchf))).eq.0.0) cycle ! closed channel
            kon(nchf,nchi) = v(npk(nchf),nchi,3)
            kon(nchi,nchf) = kon(nchf,nchi)
         enddo
      enddo 
      

      
      deallocate(v)
      return
      end
      
      subroutine makevni(ns,nds,npk,nchtop,vmat,v)
      real vmat(nds,nds+1),v(nds,nchtop)
      integer npk(nchtop+1)

      do nchi = 1, nchtop
         if (ns.eq.0) then
            do nchf = nchi, nchtop
               do kf = npk(nchf), npk(nchf+1) - 1
                  v(kf,nchi) = vmat(kf,npk(nchi))
               enddo
            enddo
            do nchf = 1, nchi - 1
               do kf = npk(nchf), npk(nchf+1) - 1
                  v(kf,nchi) = vmat(npk(nchi),kf)
               enddo
            enddo 
         else 
            do nchf = nchi, nchtop
               do kf = npk(nchf), npk(nchf+1) - 1
                  v(kf,nchi) = vmat(npk(nchi),kf+1)
               enddo
            enddo 
            do nchf = 1, nchi - 1
               do kf = npk(nchf), npk(nchf+1) - 1
                  v(kf,nchi) = vmat(kf,npk(nchi)+1)
               enddo
            enddo 
         endif 
      enddo
      return
      end

      subroutine makevfn(ns,nds,npk,nchtop,vmat,kf,v)
      real vmat(nds,nds+1),v(nds)
      integer npk(nchtop+1)

      if (ns.eq.0) then
         do n = 1, kf
            v(n) = vmat(kf,n)
         enddo
         do n = kf+1, nds
            v(n) = vmat(n,kf)
         enddo 
      else 
         do n = 1, kf
            v(n) = vmat(n,kf+1)
         enddo
         do n = kf+1, nds
            v(n) = vmat(kf,n+1)
         enddo 
      endif 
      return
      end


c$$$      SUBROUTINE ZGESV_mine( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c$$$*
c$$$*  -- LAPACK driver routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            INFO, LDA, LDB, N, NRHS
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX*16         A( LDA, * ), B( LDB, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  ZGESV computes the solution to a complex system of linear equations
c$$$*     A * X = B,
c$$$*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
c$$$*
c$$$*  The LU decomposition with partial pivoting and row interchanges is
c$$$*  used to factor A as
c$$$*     A = P * L * U,
c$$$*  where P is a permutation matrix, L is unit lower triangular, and U is
c$$$*  upper triangular.  The factored form of A is then used to solve the
c$$$*  system of equations A * X = B.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The number of linear equations, i.e., the order of the
c$$$*          matrix A.  N >= 0.
c$$$*
c$$$*  NRHS    (input) INTEGER
c$$$*          The number of right hand sides, i.e., the number of columns
c$$$*          of the matrix B.  NRHS >= 0.
c$$$*
c$$$*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
c$$$*          On entry, the N-by-N coefficient matrix A.
c$$$*          On exit, the factors L and U from the factorization
c$$$*          A = P*L*U; the unit diagonal elements of L are not stored.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  IPIV    (output) INTEGER array, dimension (N)
c$$$*          The pivot indices that define the permutation matrix P;
c$$$*          row i of the matrix was interchanged with row IPIV(i).
c$$$*
c$$$*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
c$$$*          On entry, the N-by-NRHS matrix of right hand side matrix B.
c$$$*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
c$$$*
c$$$*  LDB     (input) INTEGER
c$$$*          The leading dimension of the array B.  LDB >= max(1,N).
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value
c$$$*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
c$$$*                has been completed, but the factor U is exactly
c$$$*                singular, so the solution could not be computed.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           XERBLA, ZGETRF, ZGETRS
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      INFO = 0
c$$$      IF( N.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( NRHS.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -4
c$$$      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -7
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'ZGESV ', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Compute the LU factorization of A.
c$$$*
c$$$      CALL ZGETRF( N, N, A, LDA, IPIV, INFO )
c$$$
c$$$      IF( INFO.EQ.0 ) THEN
c$$$*
c$$$*        Solve the system A*X = B, overwriting B with X.
c$$$*
c$$$         CALL ZGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
c$$$     $                INFO )
c$$$      END IF
c$$$      RETURN
c$$$*
c$$$*     End of ZGESV
c$$$*
c$$$      END
c$$$
c$$$
c$$$
c$$$      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c$$$*
c$$$*  -- LAPACK routine (version 3.1) --
c$$$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      CHARACTER          TRANS
c$$$      INTEGER            INFO, LDA, LDB, N, NRHS
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX*16         A( LDA, * ), B( LDB, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  ZGETRS solves a system of linear equations
c$$$*     A * X = B,  A**T * X = B,  or  A**H * X = B
c$$$*  with a general N-by-N matrix A using the LU factorization computed
c$$$*  by ZGETRF.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  TRANS   (input) CHARACTER*1
c$$$*          Specifies the form of the system of equations:
c$$$*          = 'N':  A * X = B     (No transpose)
c$$$*          = 'T':  A**T * X = B  (Transpose)
c$$$*          = 'C':  A**H * X = B  (Conjugate transpose)
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The order of the matrix A.  N >= 0.
c$$$*
c$$$*  NRHS    (input) INTEGER
c$$$*          The number of right hand sides, i.e., the number of columns
c$$$*          of the matrix B.  NRHS >= 0.
c$$$*
c$$$*  A       (input) COMPLEX*16 array, dimension (LDA,N)
c$$$*          The factors L and U from the factorization A = P*L*U
c$$$*          as computed by ZGETRF.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  IPIV    (input) INTEGER array, dimension (N)
c$$$*          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
c$$$*          matrix was interchanged with row IPIV(i).
c$$$*
c$$$*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
c$$$*          On entry, the right hand side matrix B.
c$$$*          On exit, the solution matrix X.
c$$$*
c$$$*  LDB     (input) INTEGER
c$$$*          The leading dimension of the array B.  LDB >= max(1,N).
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      COMPLEX*16         ONE
c$$$      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      LOGICAL            NOTRAN
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      LOGICAL            LSAME
c$$$      EXTERNAL           LSAME
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           XERBLA, ZLASWP, ZTRSM
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      INFO = 0
c$$$      NOTRAN = LSAME( TRANS, 'N' )
c$$$      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
c$$$     $    LSAME( TRANS, 'C' ) ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( NRHS.LT.0 ) THEN
c$$$         INFO = -3
c$$$      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -5
c$$$      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -8
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'ZGETRS', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      IF( N.EQ.0 .OR. NRHS.EQ.0 )
c$$$     $   RETURN
c$$$*
c$$$      IF( NOTRAN ) THEN
c$$$*
c$$$*        Solve A * X = B.
c$$$*
c$$$*        Apply row interchanges to the right hand sides.
c$$$*
c$$$         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
c$$$*
c$$$*        Solve L*X = B, overwriting B with X.
c$$$*
c$$$         CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
c$$$     $               ONE, A, LDA, B, LDB )
c$$$*
c$$$*        Solve U*X = B, overwriting B with X.
c$$$*
c$$$         CALL ZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
c$$$     $               NRHS, ONE, A, LDA, B, LDB )
c$$$      ELSE
c$$$*
c$$$*        Solve A**T * X = B  or A**H * X = B.
c$$$*
c$$$*        Solve U'*X = B, overwriting B with X.
c$$$*
c$$$         CALL ZTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,
c$$$     $               A, LDA, B, LDB )
c$$$*
c$$$*        Solve L'*X = B, overwriting B with X.
c$$$*
c$$$         CALL ZTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,
c$$$     $               LDA, B, LDB )
c$$$*
c$$$*        Apply row interchanges to the solution vectors.
c$$$*
c$$$         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
c$$$      END IF
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of ZGETRS
c$$$*
c$$$      END
c$$$
c$$$
c$$$
c$$$      SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
c$$$*     .. Scalar Arguments ..
c$$$      DOUBLE COMPLEX ALPHA
c$$$      INTEGER LDA,LDB,M,N
c$$$      CHARACTER DIAG,SIDE,TRANSA,UPLO
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      DOUBLE COMPLEX A(LDA,*),B(LDB,*)
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  ZTRSM  solves one of the matrix equations
c$$$*
c$$$*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
c$$$*
c$$$*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
c$$$*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c$$$*
c$$$*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
c$$$*
c$$$*  The matrix X is overwritten on B.
c$$$*
c$$$*  Arguments
c$$$*  ==========
c$$$*
c$$$*  SIDE   - CHARACTER*1.
c$$$*           On entry, SIDE specifies whether op( A ) appears on the left
c$$$*           or right of X as follows:
c$$$*
c$$$*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
c$$$*
c$$$*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
c$$$*
c$$$*           Unchanged on exit.
c$$$*
c$$$*  UPLO   - CHARACTER*1.
c$$$*           On entry, UPLO specifies whether the matrix A is an upper or
c$$$*           lower triangular matrix as follows:
c$$$*
c$$$*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c$$$*
c$$$*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c$$$*
c$$$*           Unchanged on exit.
c$$$*
c$$$*  TRANSA - CHARACTER*1.
c$$$*           On entry, TRANSA specifies the form of op( A ) to be used in
c$$$*           the matrix multiplication as follows:
c$$$*
c$$$*              TRANSA = 'N' or 'n'   op( A ) = A.
c$$$*
c$$$*              TRANSA = 'T' or 't'   op( A ) = A'.
c$$$*
c$$$*              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
c$$$*
c$$$*           Unchanged on exit.
c$$$*
c$$$*  DIAG   - CHARACTER*1.
c$$$*           On entry, DIAG specifies whether or not A is unit triangular
c$$$*           as follows:
c$$$*
c$$$*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c$$$*
c$$$*              DIAG = 'N' or 'n'   A is not assumed to be unit
c$$$*                                  triangular.
c$$$*
c$$$*           Unchanged on exit.
c$$$*
c$$$*  M      - INTEGER.
c$$$*           On entry, M specifies the number of rows of B. M must be at
c$$$*           least zero.
c$$$*           Unchanged on exit.
c$$$*
c$$$*  N      - INTEGER.
c$$$*           On entry, N specifies the number of columns of B.  N must be
c$$$*           at least zero.
c$$$*           Unchanged on exit.
c$$$*
c$$$*  ALPHA  - COMPLEX*16      .
c$$$*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c$$$*           zero then  A is not referenced and  B need not be set before
c$$$*           entry.
c$$$*           Unchanged on exit.
c$$$*
c$$$*  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
c$$$*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c$$$*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c$$$*           upper triangular part of the array  A must contain the upper
c$$$*           triangular matrix  and the strictly lower triangular part of
c$$$*           A is not referenced.
c$$$*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c$$$*           lower triangular part of the array  A must contain the lower
c$$$*           triangular matrix  and the strictly upper triangular part of
c$$$*           A is not referenced.
c$$$*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c$$$*           A  are not referenced either,  but are assumed to be  unity.
c$$$*           Unchanged on exit.
c$$$*
c$$$*  LDA    - INTEGER.
c$$$*           On entry, LDA specifies the first dimension of A as declared
c$$$*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c$$$*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c$$$*           then LDA must be at least max( 1, n ).
c$$$*           Unchanged on exit.
c$$$*
c$$$*  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
c$$$*           Before entry,  the leading  m by n part of the array  B must
c$$$*           contain  the  right-hand  side  matrix  B,  and  on exit  is
c$$$*           overwritten by the solution matrix  X.
c$$$*
c$$$*  LDB    - INTEGER.
c$$$*           On entry, LDB specifies the first dimension of B as declared
c$$$*           in  the  calling  (sub)  program.   LDB  must  be  at  least
c$$$*           max( 1, m ).
c$$$*           Unchanged on exit.
c$$$*
c$$$*
c$$$*  Level 3 Blas routine.
c$$$*
c$$$*  -- Written on 8-February-1989.
c$$$*     Jack Dongarra, Argonne National Laboratory.
c$$$*     Iain Duff, AERE Harwell.
c$$$*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c$$$*     Sven Hammarling, Numerical Algorithms Group Ltd.
c$$$*
c$$$*
c$$$*     .. External Functions ..
c$$$      LOGICAL LSAME
c$$$      EXTERNAL LSAME
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL XERBLA
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC DCONJG,MAX
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      DOUBLE COMPLEX TEMP
c$$$      INTEGER I,INFO,J,K,NROWA
c$$$      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
c$$$*     ..
c$$$*     .. Parameters ..
c$$$      DOUBLE COMPLEX ONE
c$$$      PARAMETER (ONE= (1.0D+0,0.0D+0))
c$$$      DOUBLE COMPLEX ZERO
c$$$      PARAMETER (ZERO= (0.0D+0,0.0D+0))
c$$$*     ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      LSIDE = LSAME(SIDE,'L')
c$$$      IF (LSIDE) THEN
c$$$          NROWA = M
c$$$      ELSE
c$$$          NROWA = N
c$$$      END IF
c$$$      NOCONJ = LSAME(TRANSA,'T')
c$$$      NOUNIT = LSAME(DIAG,'N')
c$$$      UPPER = LSAME(UPLO,'U')
c$$$*
c$$$      INFO = 0
c$$$      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
c$$$          INFO = 1
c$$$      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
c$$$          INFO = 2
c$$$      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
c$$$     +         (.NOT.LSAME(TRANSA,'T')) .AND.
c$$$     +         (.NOT.LSAME(TRANSA,'C'))) THEN
c$$$          INFO = 3
c$$$      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
c$$$          INFO = 4
c$$$      ELSE IF (M.LT.0) THEN
c$$$          INFO = 5
c$$$      ELSE IF (N.LT.0) THEN
c$$$          INFO = 6
c$$$      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
c$$$          INFO = 9
c$$$      ELSE IF (LDB.LT.MAX(1,M)) THEN
c$$$          INFO = 11
c$$$      END IF
c$$$      IF (INFO.NE.0) THEN
c$$$          CALL XERBLA('ZTRSM ',INFO)
c$$$          RETURN
c$$$      END IF
c$$$*
c$$$*     Quick return if possible.
c$$$*
c$$$      IF (N.EQ.0) RETURN
c$$$*
c$$$*     And when  alpha.eq.zero.
c$$$*
c$$$      IF (ALPHA.EQ.ZERO) THEN
c$$$          DO 20 J = 1,N
c$$$              DO 10 I = 1,M
c$$$                  B(I,J) = ZERO
c$$$   10         CONTINUE
c$$$   20     CONTINUE
c$$$          RETURN
c$$$      END IF
c$$$*
c$$$*     Start the operations.
c$$$*
c$$$      IF (LSIDE) THEN
c$$$          IF (LSAME(TRANSA,'N')) THEN
c$$$*
c$$$*           Form  B := alpha*inv( A )*B.
c$$$*
c$$$              IF (UPPER) THEN
c$$$                  DO 60 J = 1,N
c$$$                      IF (ALPHA.NE.ONE) THEN
c$$$                          DO 30 I = 1,M
c$$$                              B(I,J) = ALPHA*B(I,J)
c$$$   30                     CONTINUE
c$$$                      END IF
c$$$                      DO 50 K = M,1,-1
c$$$                          IF (B(K,J).NE.ZERO) THEN
c$$$                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
c$$$                              DO 40 I = 1,K - 1
c$$$                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
c$$$   40                         CONTINUE
c$$$                          END IF
c$$$   50                 CONTINUE
c$$$   60             CONTINUE
c$$$              ELSE
c$$$                  DO 100 J = 1,N
c$$$                      IF (ALPHA.NE.ONE) THEN
c$$$                          DO 70 I = 1,M
c$$$                              B(I,J) = ALPHA*B(I,J)
c$$$   70                     CONTINUE
c$$$                      END IF
c$$$                      DO 90 K = 1,M
c$$$                          IF (B(K,J).NE.ZERO) THEN
c$$$                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
c$$$                              DO 80 I = K + 1,M
c$$$                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
c$$$   80                         CONTINUE
c$$$                          END IF
c$$$   90                 CONTINUE
c$$$  100             CONTINUE
c$$$              END IF
c$$$          ELSE
c$$$*
c$$$*           Form  B := alpha*inv( A' )*B
c$$$*           or    B := alpha*inv( conjg( A' ) )*B.
c$$$*
c$$$              IF (UPPER) THEN
c$$$                  DO 140 J = 1,N
c$$$                      DO 130 I = 1,M
c$$$                          TEMP = ALPHA*B(I,J)
c$$$                          IF (NOCONJ) THEN
c$$$                              DO 110 K = 1,I - 1
c$$$                                  TEMP = TEMP - A(K,I)*B(K,J)
c$$$  110                         CONTINUE
c$$$                              IF (NOUNIT) TEMP = TEMP/A(I,I)
c$$$                          ELSE
c$$$                              DO 120 K = 1,I - 1
c$$$                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
c$$$  120                         CONTINUE
c$$$                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
c$$$                          END IF
c$$$                          B(I,J) = TEMP
c$$$  130                 CONTINUE
c$$$  140             CONTINUE
c$$$              ELSE
c$$$                  DO 180 J = 1,N
c$$$                      DO 170 I = M,1,-1
c$$$                          TEMP = ALPHA*B(I,J)
c$$$                          IF (NOCONJ) THEN
c$$$                              DO 150 K = I + 1,M
c$$$                                  TEMP = TEMP - A(K,I)*B(K,J)
c$$$  150                         CONTINUE
c$$$                              IF (NOUNIT) TEMP = TEMP/A(I,I)
c$$$                          ELSE
c$$$                              DO 160 K = I + 1,M
c$$$                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
c$$$  160                         CONTINUE
c$$$                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
c$$$                          END IF
c$$$                          B(I,J) = TEMP
c$$$  170                 CONTINUE
c$$$  180             CONTINUE
c$$$              END IF
c$$$          END IF
c$$$      ELSE
c$$$          IF (LSAME(TRANSA,'N')) THEN
c$$$*
c$$$*           Form  B := alpha*B*inv( A ).
c$$$*
c$$$              IF (UPPER) THEN
c$$$                  DO 230 J = 1,N
c$$$                      IF (ALPHA.NE.ONE) THEN
c$$$                          DO 190 I = 1,M
c$$$                              B(I,J) = ALPHA*B(I,J)
c$$$  190                     CONTINUE
c$$$                      END IF
c$$$                      DO 210 K = 1,J - 1
c$$$                          IF (A(K,J).NE.ZERO) THEN
c$$$                              DO 200 I = 1,M
c$$$                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
c$$$  200                         CONTINUE
c$$$                          END IF
c$$$  210                 CONTINUE
c$$$                      IF (NOUNIT) THEN
c$$$                          TEMP = ONE/A(J,J)
c$$$                          DO 220 I = 1,M
c$$$                              B(I,J) = TEMP*B(I,J)
c$$$  220                     CONTINUE
c$$$                      END IF
c$$$  230             CONTINUE
c$$$              ELSE
c$$$                  DO 280 J = N,1,-1
c$$$                      IF (ALPHA.NE.ONE) THEN
c$$$                          DO 240 I = 1,M
c$$$                              B(I,J) = ALPHA*B(I,J)
c$$$  240                     CONTINUE
c$$$                      END IF
c$$$                      DO 260 K = J + 1,N
c$$$                          IF (A(K,J).NE.ZERO) THEN
c$$$                              DO 250 I = 1,M
c$$$                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
c$$$  250                         CONTINUE
c$$$                          END IF
c$$$  260                 CONTINUE
c$$$                      IF (NOUNIT) THEN
c$$$                          TEMP = ONE/A(J,J)
c$$$                          DO 270 I = 1,M
c$$$                              B(I,J) = TEMP*B(I,J)
c$$$  270                     CONTINUE
c$$$                      END IF
c$$$  280             CONTINUE
c$$$              END IF
c$$$          ELSE
c$$$*
c$$$*           Form  B := alpha*B*inv( A' )
c$$$*           or    B := alpha*B*inv( conjg( A' ) ).
c$$$*
c$$$              IF (UPPER) THEN
c$$$                  DO 330 K = N,1,-1
c$$$                      IF (NOUNIT) THEN
c$$$                          IF (NOCONJ) THEN
c$$$                              TEMP = ONE/A(K,K)
c$$$                          ELSE
c$$$                              TEMP = ONE/DCONJG(A(K,K))
c$$$                          END IF
c$$$                          DO 290 I = 1,M
c$$$                              B(I,K) = TEMP*B(I,K)
c$$$  290                     CONTINUE
c$$$                      END IF
c$$$                      DO 310 J = 1,K - 1
c$$$                          IF (A(J,K).NE.ZERO) THEN
c$$$                              IF (NOCONJ) THEN
c$$$                                  TEMP = A(J,K)
c$$$                              ELSE
c$$$                                  TEMP = DCONJG(A(J,K))
c$$$                              END IF
c$$$                              DO 300 I = 1,M
c$$$                                  B(I,J) = B(I,J) - TEMP*B(I,K)
c$$$  300                         CONTINUE
c$$$                          END IF
c$$$  310                 CONTINUE
c$$$                      IF (ALPHA.NE.ONE) THEN
c$$$                          DO 320 I = 1,M
c$$$                              B(I,K) = ALPHA*B(I,K)
c$$$  320                     CONTINUE
c$$$                      END IF
c$$$  330             CONTINUE
c$$$              ELSE
c$$$                  DO 380 K = 1,N
c$$$                      IF (NOUNIT) THEN
c$$$                          IF (NOCONJ) THEN
c$$$                              TEMP = ONE/A(K,K)
c$$$                          ELSE
c$$$                              TEMP = ONE/DCONJG(A(K,K))
c$$$                          END IF
c$$$                          DO 340 I = 1,M
c$$$                              B(I,K) = TEMP*B(I,K)
c$$$  340                     CONTINUE
c$$$                      END IF
c$$$                      DO 360 J = K + 1,N
c$$$                          IF (A(J,K).NE.ZERO) THEN
c$$$                              IF (NOCONJ) THEN
c$$$                                  TEMP = A(J,K)
c$$$                              ELSE
c$$$                                  TEMP = DCONJG(A(J,K))
c$$$                              END IF
c$$$                              DO 350 I = 1,M
c$$$                                  B(I,J) = B(I,J) - TEMP*B(I,K)
c$$$  350                         CONTINUE
c$$$                          END IF
c$$$  360                 CONTINUE
c$$$                      IF (ALPHA.NE.ONE) THEN
c$$$                          DO 370 I = 1,M
c$$$                              B(I,K) = ALPHA*B(I,K)
c$$$  370                     CONTINUE
c$$$                      END IF
c$$$  380             CONTINUE
c$$$              END IF
c$$$          END IF
c$$$      END IF
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of ZTRSM .
c$$$*
c$$$      END
c$$$
c$$$      SUBROUTINE CGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c$$$*
c$$$*  -- LAPACK driver routine (version 3.2) --
c$$$*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c$$$*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            INFO, LDA, LDB, N, NRHS
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX            A( LDA, * ), B( LDB, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGESV computes the solution to a complex system of linear equations
c$$$*     A * X = B,
c$$$*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
c$$$*
c$$$*  The LU decomposition with partial pivoting and row interchanges is
c$$$*  used to factor A as
c$$$*     A = P * L * U,
c$$$*  where P is a permutation matrix, L is unit lower triangular, and U is
c$$$*  upper triangular.  The factored form of A is then used to solve the
c$$$*  system of equations A * X = B.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The number of linear equations, i.e., the order of the
c$$$*          matrix A.  N >= 0.
c$$$*
c$$$*  NRHS    (input) INTEGER
c$$$*          The number of right hand sides, i.e., the number of columns
c$$$*          of the matrix B.  NRHS >= 0.
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N)
c$$$*          On entry, the N-by-N coefficient matrix A.
c$$$*          On exit, the factors L and U from the factorization
c$$$*          A = P*L*U; the unit diagonal elements of L are not stored.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  IPIV    (output) INTEGER array, dimension (N)
c$$$*          The pivot indices that define the permutation matrix P;
c$$$*          row i of the matrix was interchanged with row IPIV(i).
c$$$*
c$$$*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
c$$$*          On entry, the N-by-NRHS matrix of right hand side matrix B.
c$$$*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
c$$$*
c$$$*  LDB     (input) INTEGER
c$$$*          The leading dimension of the array B.  LDB >= max(1,N).
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value
c$$$*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
c$$$*                has been completed, but the factor U is exactly
c$$$*                singular, so the solution could not be computed.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CGETRF, CGETRS, XERBLA
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      INFO = 0
c$$$      IF( N.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( NRHS.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -4
c$$$      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -7
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGESV ', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Compute the LU factorization of A.
c$$$*
c$$$      CALL CGETRF( N, N, A, LDA, IPIV, INFO )
c$$$      IF( INFO.EQ.0 ) THEN
c$$$*
c$$$*        Solve the system A*X = B, overwriting B with X.
c$$$*
c$$$         CALL CGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
c$$$     $                INFO )
c$$$      END IF
c$$$      RETURN
c$$$*
c$$$*     End of CGESV
c$$$*
c$$$      END
c$$$      SUBROUTINE CGETF2( M, N, A, LDA, IPIV, INFO )
c$$$*
c$$$*  -- LAPACK routine (version 3.2) --
c$$$*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c$$$*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            INFO, LDA, M, N
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX            A( LDA, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGETF2 computes an LU factorization of a general m-by-n matrix A
c$$$*  using partial pivoting with row interchanges.
c$$$*
c$$$*  The factorization has the form
c$$$*     A = P * L * U
c$$$*  where P is a permutation matrix, L is lower triangular with unit
c$$$*  diagonal elements (lower trapezoidal if m > n), and U is upper
c$$$*  triangular (upper trapezoidal if m < n).
c$$$*
c$$$*  This is the right-looking Level 2 BLAS version of the algorithm.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  M       (input) INTEGER
c$$$*          The number of rows of the matrix A.  M >= 0.
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The number of columns of the matrix A.  N >= 0.
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N)
c$$$*          On entry, the m by n matrix to be factored.
c$$$*          On exit, the factors L and U from the factorization
c$$$*          A = P*L*U; the unit diagonal elements of L are not stored.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,M).
c$$$*
c$$$*  IPIV    (output) INTEGER array, dimension (min(M,N))
c$$$*          The pivot indices; for 1 <= i <= min(M,N), row i of the
c$$$*          matrix was interchanged with row IPIV(i).
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0: successful exit
c$$$*          < 0: if INFO = -k, the k-th argument had an illegal value
c$$$*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
c$$$*               has been completed, but the factor U is exactly
c$$$*               singular, and division by zero will occur if it is used
c$$$*               to solve a system of equations.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      COMPLEX            ONE, ZERO
c$$$      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
c$$$     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      REAL               SFMIN
c$$$      INTEGER            I, J, JP
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      REAL               SLAMCH
c$$$      INTEGER            ICAMAX
c$$$      EXTERNAL           SLAMCH, ICAMAX
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CGERU, CSCAL, CSWAP, XERBLA
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX, MIN
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      INFO = 0
c$$$      IF( M.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
c$$$         INFO = -4
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGETF2', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      IF( M.EQ.0 .OR. N.EQ.0 )
c$$$     $   RETURN
c$$$*
c$$$*     Compute machine safe minimum
c$$$*
c$$$      SFMIN = SLAMCH('S') 
c$$$*
c$$$      DO 10 J = 1, MIN( M, N )
c$$$*
c$$$*        Find pivot and test for singularity.
c$$$*
c$$$         JP = J - 1 + ICAMAX( M-J+1, A( J, J ), 1 )
c$$$         IPIV( J ) = JP
c$$$         IF( A( JP, J ).NE.ZERO ) THEN
c$$$*
c$$$*           Apply the interchange to columns 1:N.
c$$$*
c$$$            IF( JP.NE.J )
c$$$     $         CALL CSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
c$$$*
c$$$*           Compute elements J+1:M of J-th column.
c$$$*
c$$$            IF( J.LT.M ) THEN
c$$$               IF( ABS(A( J, J )) .GE. SFMIN ) THEN
c$$$                  CALL CSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
c$$$               ELSE
c$$$                  DO 20 I = 1, M-J
c$$$                     A( J+I, J ) = A( J+I, J ) / A( J, J )
c$$$   20             CONTINUE
c$$$               END IF
c$$$            END IF
c$$$*
c$$$         ELSE IF( INFO.EQ.0 ) THEN
c$$$*
c$$$            INFO = J
c$$$         END IF
c$$$*
c$$$         IF( J.LT.MIN( M, N ) ) THEN
c$$$*
c$$$*           Update trailing submatrix.
c$$$*
c$$$            CALL CGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),
c$$$     $                  LDA, A( J+1, J+1 ), LDA )
c$$$         END IF
c$$$   10 CONTINUE
c$$$      RETURN
c$$$*
c$$$*     End of CGETF2
c$$$*
c$$$      END
c$$$      SUBROUTINE CGETRF( M, N, A, LDA, IPIV, INFO )
c$$$*
c$$$*  -- LAPACK routine (version 3.2) --
c$$$*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c$$$*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            INFO, LDA, M, N
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX            A( LDA, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGETRF computes an LU factorization of a general M-by-N matrix A
c$$$*  using partial pivoting with row interchanges.
c$$$*
c$$$*  The factorization has the form
c$$$*     A = P * L * U
c$$$*  where P is a permutation matrix, L is lower triangular with unit
c$$$*  diagonal elements (lower trapezoidal if m > n), and U is upper
c$$$*  triangular (upper trapezoidal if m < n).
c$$$*
c$$$*  This is the right-looking Level 3 BLAS version of the algorithm.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  M       (input) INTEGER
c$$$*          The number of rows of the matrix A.  M >= 0.
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The number of columns of the matrix A.  N >= 0.
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N)
c$$$*          On entry, the M-by-N matrix to be factored.
c$$$*          On exit, the factors L and U from the factorization
c$$$*          A = P*L*U; the unit diagonal elements of L are not stored.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,M).
c$$$*
c$$$*  IPIV    (output) INTEGER array, dimension (min(M,N))
c$$$*          The pivot indices; for 1 <= i <= min(M,N), row i of the
c$$$*          matrix was interchanged with row IPIV(i).
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value
c$$$*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
c$$$*                has been completed, but the factor U is exactly
c$$$*                singular, and division by zero will occur if it is used
c$$$*                to solve a system of equations.
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      COMPLEX            ONE
c$$$      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      INTEGER            I, IINFO, J, JB, NB
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CGEMM, CGETF2, CLASWP, CTRSM, XERBLA
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      INTEGER            ILAENV
c$$$      EXTERNAL           ILAENV
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX, MIN
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      INFO = 0
c$$$      IF( M.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
c$$$         INFO = -4
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGETRF', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      IF( M.EQ.0 .OR. N.EQ.0 )
c$$$     $   RETURN
c$$$*
c$$$*     Determine the block size for this environment.
c$$$*
c$$$      NB = ILAENV( 1, 'CGETRF', ' ', M, N, -1, -1 )
c$$$      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
c$$$*
c$$$*        Use unblocked code.
c$$$*
c$$$         CALL CGETF2( M, N, A, LDA, IPIV, INFO )
c$$$      ELSE
c$$$*
c$$$*        Use blocked code.
c$$$*
c$$$         DO 20 J = 1, MIN( M, N ), NB
c$$$            JB = MIN( MIN( M, N )-J+1, NB )
c$$$*
c$$$*           Factor diagonal and subdiagonal blocks and test for exact
c$$$*           singularity.
c$$$*
c$$$            CALL CGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
c$$$*
c$$$*           Adjust INFO and the pivot indices.
c$$$*
c$$$            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
c$$$     $         INFO = IINFO + J - 1
c$$$            DO 10 I = J, MIN( M, J+JB-1 )
c$$$               IPIV( I ) = J - 1 + IPIV( I )
c$$$   10       CONTINUE
c$$$*
c$$$*           Apply interchanges to columns 1:J-1.
c$$$*
c$$$            CALL CLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
c$$$*
c$$$            IF( J+JB.LE.N ) THEN
c$$$*
c$$$*              Apply interchanges to columns J+JB:N.
c$$$*
c$$$               CALL CLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
c$$$     $                      IPIV, 1 )
c$$$*
c$$$*              Compute block row of U.
c$$$*
c$$$               CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
c$$$     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
c$$$     $                     LDA )
c$$$               IF( J+JB.LE.M ) THEN
c$$$*
c$$$*                 Update trailing submatrix.
c$$$*
c$$$                  CALL CGEMM( 'No transpose', 'No transpose', M-J-JB+1,
c$$$     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
c$$$     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
c$$$     $                        LDA )
c$$$               END IF
c$$$            END IF
c$$$   20    CONTINUE
c$$$      END IF
c$$$      RETURN
c$$$*
c$$$*     End of CGETRF
c$$$*
c$$$      END
c$$$      SUBROUTINE CGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c$$$*
c$$$*  -- LAPACK routine (version 3.3.1) --
c$$$*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c$$$*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c$$$*  -- April 2011                                                      --
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      CHARACTER          TRANS
c$$$      INTEGER            INFO, LDA, LDB, N, NRHS
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX            A( LDA, * ), B( LDB, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CGETRS solves a system of linear equations
c$$$*     A * X = B,  A**T * X = B,  or  A**H * X = B
c$$$*  with a general N-by-N matrix A using the LU factorization computed
c$$$*  by CGETRF.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  TRANS   (input) CHARACTER*1
c$$$*          Specifies the form of the system of equations:
c$$$*          = 'N':  A * X = B     (No transpose)
c$$$*          = 'T':  A**T * X = B  (Transpose)
c$$$*          = 'C':  A**H * X = B  (Conjugate transpose)
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The order of the matrix A.  N >= 0.
c$$$*
c$$$*  NRHS    (input) INTEGER
c$$$*          The number of right hand sides, i.e., the number of columns
c$$$*          of the matrix B.  NRHS >= 0.
c$$$*
c$$$*  A       (input) COMPLEX array, dimension (LDA,N)
c$$$*          The factors L and U from the factorization A = P*L*U
c$$$*          as computed by CGETRF.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.  LDA >= max(1,N).
c$$$*
c$$$*  IPIV    (input) INTEGER array, dimension (N)
c$$$*          The pivot indices from CGETRF; for 1<=i<=N, row i of the
c$$$*          matrix was interchanged with row IPIV(i).
c$$$*
c$$$*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
c$$$*          On entry, the right hand side matrix B.
c$$$*          On exit, the solution matrix X.
c$$$*
c$$$*  LDB     (input) INTEGER
c$$$*          The leading dimension of the array B.  LDB >= max(1,N).
c$$$*
c$$$*  INFO    (output) INTEGER
c$$$*          = 0:  successful exit
c$$$*          < 0:  if INFO = -i, the i-th argument had an illegal value
c$$$*
c$$$*  =====================================================================
c$$$*
c$$$*     .. Parameters ..
c$$$      COMPLEX            ONE
c$$$      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
c$$$*     ..
c$$$*     .. Local Scalars ..
c$$$      LOGICAL            NOTRAN
c$$$*     ..
c$$$*     .. External Functions ..
c$$$      LOGICAL            LSAME
c$$$      EXTERNAL           LSAME
c$$$*     ..
c$$$*     .. External Subroutines ..
c$$$      EXTERNAL           CLASWP, CTRSM, XERBLA
c$$$*     ..
c$$$*     .. Intrinsic Functions ..
c$$$      INTRINSIC          MAX
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Test the input parameters.
c$$$*
c$$$      INFO = 0
c$$$      NOTRAN = LSAME( TRANS, 'N' )
c$$$      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
c$$$     $    LSAME( TRANS, 'C' ) ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( NRHS.LT.0 ) THEN
c$$$         INFO = -3
c$$$      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -5
c$$$      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
c$$$         INFO = -8
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'CGETRS', -INFO )
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$*     Quick return if possible
c$$$*
c$$$      IF( N.EQ.0 .OR. NRHS.EQ.0 )
c$$$     $   RETURN
c$$$*
c$$$      IF( NOTRAN ) THEN
c$$$*
c$$$*        Solve A * X = B.
c$$$*
c$$$*        Apply row interchanges to the right hand sides.
c$$$*
c$$$         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
c$$$*
c$$$*        Solve L*X = B, overwriting B with X.
c$$$*
c$$$         CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
c$$$     $               ONE, A, LDA, B, LDB )
c$$$*
c$$$*        Solve U*X = B, overwriting B with X.
c$$$*
c$$$         CALL CTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
c$$$     $               NRHS, ONE, A, LDA, B, LDB )
c$$$      ELSE
c$$$*
c$$$*        Solve A**T * X = B  or A**H * X = B.
c$$$*
c$$$*        Solve U**T *X = B or U**H *X = B, overwriting B with X.
c$$$*
c$$$         CALL CTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,
c$$$     $               A, LDA, B, LDB )
c$$$*
c$$$*        Solve L**T *X = B, or L**H *X = B overwriting B with X.
c$$$*
c$$$         CALL CTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,
c$$$     $               LDA, B, LDB )
c$$$*
c$$$*        Apply row interchanges to the solution vectors.
c$$$*
c$$$         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
c$$$      END IF
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of CGETRS
c$$$*
c$$$      END
c$$$      SUBROUTINE CLASWP( N, A, LDA, K1, K2, IPIV, INCX )
c$$$*
c$$$*  -- LAPACK auxiliary routine (version 3.2) --
c$$$*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c$$$*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c$$$*     November 2006
c$$$*
c$$$*     .. Scalar Arguments ..
c$$$      INTEGER            INCX, K1, K2, LDA, N
c$$$*     ..
c$$$*     .. Array Arguments ..
c$$$      INTEGER            IPIV( * )
c$$$      COMPLEX            A( LDA, * )
c$$$*     ..
c$$$*
c$$$*  Purpose
c$$$*  =======
c$$$*
c$$$*  CLASWP performs a series of row interchanges on the matrix A.
c$$$*  One row interchange is initiated for each of rows K1 through K2 of A.
c$$$*
c$$$*  Arguments
c$$$*  =========
c$$$*
c$$$*  N       (input) INTEGER
c$$$*          The number of columns of the matrix A.
c$$$*
c$$$*  A       (input/output) COMPLEX array, dimension (LDA,N)
c$$$*          On entry, the matrix of column dimension N to which the row
c$$$*          interchanges will be applied.
c$$$*          On exit, the permuted matrix.
c$$$*
c$$$*  LDA     (input) INTEGER
c$$$*          The leading dimension of the array A.
c$$$*
c$$$*  K1      (input) INTEGER
c$$$*          The first element of IPIV for which a row interchange will
c$$$*          be done.
c$$$*
c$$$*  K2      (input) INTEGER
c$$$*          The last element of IPIV for which a row interchange will
c$$$*          be done.
c$$$*
c$$$*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
c$$$*          The vector of pivot indices.  Only the elements in positions
c$$$*          K1 through K2 of IPIV are accessed.
c$$$*          IPIV(K) = L implies rows K and L are to be interchanged.
c$$$*
c$$$*  INCX    (input) INTEGER
c$$$*          The increment between successive values of IPIV.  If IPIV
c$$$*          is negative, the pivots are applied in reverse order.
c$$$*
c$$$*  Further Details
c$$$*  ===============
c$$$*
c$$$*  Modified by
c$$$*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
c$$$*
c$$$* =====================================================================
c$$$*
c$$$*     .. Local Scalars ..
c$$$      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
c$$$      COMPLEX            TEMP
c$$$*     ..
c$$$*     .. Executable Statements ..
c$$$*
c$$$*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
c$$$*
c$$$      IF( INCX.GT.0 ) THEN
c$$$         IX0 = K1
c$$$         I1 = K1
c$$$         I2 = K2
c$$$         INC = 1
c$$$      ELSE IF( INCX.LT.0 ) THEN
c$$$         IX0 = 1 + ( 1-K2 )*INCX
c$$$         I1 = K2
c$$$         I2 = K1
c$$$         INC = -1
c$$$      ELSE
c$$$         RETURN
c$$$      END IF
c$$$*
c$$$      N32 = ( N / 32 )*32
c$$$      IF( N32.NE.0 ) THEN
c$$$         DO 30 J = 1, N32, 32
c$$$            IX = IX0
c$$$            DO 20 I = I1, I2, INC
c$$$               IP = IPIV( IX )
c$$$               IF( IP.NE.I ) THEN
c$$$                  DO 10 K = J, J + 31
c$$$                     TEMP = A( I, K )
c$$$                     A( I, K ) = A( IP, K )
c$$$                     A( IP, K ) = TEMP
c$$$   10             CONTINUE
c$$$               END IF
c$$$               IX = IX + INCX
c$$$   20       CONTINUE
c$$$   30    CONTINUE
c$$$      END IF
c$$$      IF( N32.NE.N ) THEN
c$$$         N32 = N32 + 1
c$$$         IX = IX0
c$$$         DO 50 I = I1, I2, INC
c$$$            IP = IPIV( IX )
c$$$            IF( IP.NE.I ) THEN
c$$$               DO 40 K = N32, N
c$$$                  TEMP = A( I, K )
c$$$                  A( I, K ) = A( IP, K )
c$$$                  A( IP, K ) = TEMP
c$$$   40          CONTINUE
c$$$            END IF
c$$$            IX = IX + INCX
c$$$   50    CONTINUE
c$$$      END IF
c$$$*
c$$$      RETURN
c$$$*
c$$$*     End of CLASWP
c$$$*
c$$$      END

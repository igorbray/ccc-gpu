      subroutine first(ispeed,ifirst,second,nold,etot,lg,gk,npk,chil,
     >   minchil,u,ldw,dwpot,phasel,itail,nznuci,nchtop,nchtope2e,qcut,
     >   vdon,vmat,theta,vdcore,minvdc,maxvdc,lfast,lslow,slowery,td,
     >   te1,te2,ve2ed,ve2ee,dphasee2e,ephasee2e,ne2e,nchmaxe2e,
     >   vmatp,nsmax,
     >   nchistart,nchistop,nodeid,scalapack,
     >   vmat01,vmat0,vmat1,
     >   vni,vnf,vnd,nodes,myid)
#ifdef GPU
      use openacc
#endif
      use ubb_module
      use apar
!      use chil_module
C      use vmat_module
      include 'par.f'

      integer:: vni,vnf,vnd,nodes
      real, dimension(vnf+1:vnd,vni:vnf) :: vmat0
      real, dimension(vni:vnf,vni:vnf+1) :: vmat01
      real, dimension(vni:vnf,vnf+1+1:vnd+1) :: vmat1
      integer,dimension(nodes):: nchistart, nchistop
      integer::  nodeid
      logical:: scalapack

      integer npk(nchtop+1),mintemp3(nchan),maxtemp3(nchan),ltmin(nchan)
      dimension chil(meshr,npk(nchtop+1)-1,2),minchil(npk(nchtop+1)-1,2)
     >  ,u(maxr),gk(kmax,nchan),vdon(nchan,nchan,0:1),ui(maxr),
     >   uf(maxr,nchan)
     >   ,vdcore(maxr,0:lamax),u1e(maxr),dwpot(maxr,nchan),ctemp(nchan)
      complex phasel(kmax,nchan),phasefast,phaseslow,sigc,
     >   dphasee2e(nchane2e),ephasee2e(nchane2e)
      real vmat(npk(nchtop+1)-1,npk(nchtop+1)),
     >   vmatp((npk(nchtop+1)-1)*npk(nchtop+1)/2,0:nsmax)
      common/matchph/rphase(kmax,nchan),trat
      common/meshrr/ meshr,rmesh(maxr,3)
C      common/meshrr/ rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      common /double/id,jdouble(22)
      common /radpot/ ucentr(maxr)
      common /di_el_core_polarization/ gamma, r0, pol(maxr)
      common/smallr/ formcut,regcut,expcut,fast,match,analyticd,packed
      logical fast,second,posi,posf,positron,match,analyticd,packed
!      logical alkali
      dimension psii(maxr), psif(maxr), ovlp(kmax),temp(maxr),fun(maxr),
     >   psislow(maxr), psifast(maxr), slowery(ncmax)
     >   ,vmatt(kmax,kmax,0:1,nchtop)

      dimension psi_t(maxr,nchtop)
      dimension maxpsi_t(nchtop),la_t(nchtop),e_t(nchtop),
     >          na_t(nchtop),l_t(nchtop),npos_t(nchtop)

C  Note that ve2eon will not be defined correctly due to the declaration
C  in makev31d. However, ve2eon is not used.
      real ve2ed(nchmaxe2e,npk(nchtop+1)-1),
     >   ve2eon(nchane2e,nchan), uplane(maxr), ud(maxr),
     >   ve2ee(nchmaxe2e,npk(nchtop+1)-1)
      integer, allocatable :: nchansi(:), nchansf(:)
C      allocatable :: chitemp(:,:), temp3(:,:),vmati(:,:,:)
C      allocatable :: chitemp(:,:)
      real,allocatable :: vmati(:,:,:),temp3(:,:)!,vmatt(:,:,:,:)
!      data pi/3.14159265358979/
      data uplane/maxr*0.0/
 
      real*8, dimension (meshr) :: ubbb ! ANDREY
      real*8 xin(maxr),yin(maxr),xout(maxr),yout(maxr)
      logical, dimension (nchtop) :: pos

      integer ngpus,ntpg,nnt,gpunum,nthreads
      integer childim
      integer, external :: omp_get_max_threads
      integer, external :: omp_get_thread_num

C      integer nchistart(:),nchistop(:)
C      real vmat01(npk(nchistart(nodeid)):npk(nchistop(nodeid)+1)-1,
C     >     npk(nchistart(nodeid)):npk(nchistop(nodeid)+1))
C      real vmat0(npk(nchistop(nodeid)+1):npk(nchtop+1)-1,
C     >     npk(nchistart(nodeid)):npk(nchistop(nodeid)+1)-1)
C      real vmat1(npk(nchistart(nodeid)):npk(nchistop(nodeid)+1)-1,
C     >     npk(nchistop(nodeid)+1)+1:npk(nchtop+1))

      hat(l)= sqrt(2.0 * l + 1.0)

      ni=0
      nf=0

! set number of GPUs, OpenMP threads per GPU and nested threads
#ifdef GPU
      ngpus=max(1,acc_get_num_devices(acc_device_nvidia))
! 2 threads per GPU seems to be a good choice for P100 arch
      ntpg=1
      if(ngpus>0) then
        nnt=max(1,omp_get_max_threads()/(ngpus*ntpg))
        nthreads=ngpus*ntpg
      else
        nnt=1
        nthreads=max(1,omp_get_max_threads())
      endif
#else
      nnt=1
      nthreads=max(1,omp_get_max_threads())     
      write(*,*) "nnt,nthreads",nnt,nthreads
#endif

      td = 0.0
      te1 = 0.0
      te2 = 0.0
      rnorm = 2.0/pi
      if (npk(2)-npk(1).eq.1) then
         nchii = 1
         nchif = nchtop
      else
         nchii = nchistart(nodeid)
         nchif = nchistop(nodeid)
      endif 
C Unroll the nchi/nchf two loops into one over nch, for OpenMP efficiency.
      allocate(nchansi((nchif-nchii+1)*(nchtop-nchii+1)))
      allocate(nchansf((nchif-nchii+1)*(nchtop-nchii+1)))
      nch = 0
      do nchi = nchii, nchif
         do nchf = nchi, nchtop
            nch = nch + 1
            nchansi(nch) = nchi
            nchansf(nch) = nchf
         enddo
      enddo
      nchansmax = nch
      print'("nodeid, nchii, nchif, nchansmax, and allocation:",5i8)',
     >   nodeid,nchii,nchif,nchansmax,(nchif-nchii+1)*(nchtop-nchii+1)

      nch = 0

      if(itail.lt.0) then
        childim=2
      else
        childim=1
      endif

!      allocate(temp3(meshr,nchii:nchtop))
!      allocate(vmati(1:kmax,1:kmax,nchii:nchtop))
!      allocate(chitemp(meshr,kmax))


CC GPU ONLY
#ifdef GPU
      do gpunum=0,ngpus-1
         call acc_set_device_num(gpunum,acc_device_nvidia)

!$acc enter data copyin(chil(1:meshr,1:(npk(nchtop+1)-1),1:childim))
!$acc& copyin(nchtop,npk(1:nchtop+1))
! !$acc& create(chitemp(1:meshr,1:kmax))

       end do
#endif

!$omp parallel do private(nchf,nt) schedule(dynamic)
!$omp& shared(nchii,nchtop,lg,psi_t,maxpsi_t,e_t,la_t,na_t,l_t)
      do nchf = nchii, nchtop
         call getchinfo (nchf,nt,lg, psi_t(1,nchf), maxpsi_t(nchf), 
     >                 e_t(nchf), la_t(nchf), na_t(nchf),l_t(nchf))
         npos_t(nchf)=0
         pos(nchf)=positron(na_t(nchf),la_t(nchf),npos_t(nchf)) 
      end do
!$omp end parallel do


CC GPU ONLY
!$omp parallel do default(none) 
!$omp& num_threads(nthreads) schedule(dynamic)
!$omp& private(nchi,gpunum,nchf)
!$omp& private(nti,psii,maxpsii,ei,lia,nia,li,nposi,posi,ui,ud)
!$omp& private(nt,nposf,posf,nqmf,uf)
!$omp& private(td,s1,s2,s3,s4,te1,te2,nqmi)
!$omp& private(mini1,psif,kii,kff,maxpsif,ef,lfa,nfa,lf)
!$omp& private(nchtopf,nchtopi,minfun,maxfun,fun,c1tmp,c2tmp,c3tmp,mini)
!$omp& private(maxi,lt,c1,c2,c3,c2tmp2,c,const,i1,i2,i,ns,kf)
!$omp& shared(lg,nchii,nchif,npk,vdcore,minvdc,maxvdc,dwpot,meshr,ldw) 
!$omp& shared(rmesh,nsmax,rnorm,scalapack,nchtop,ltmin,second,ispeed)
!$omp& shared(ifirst,gk,nchistop,etot,nodeid,itail,nold,trat,nznuc)
!$omp& shared(vmat,vmat01,vmat0,vmat1,chil,minchil,theta,ctemp)
!$omp& shared(ve2ee,childim,ngpus,rpow1,rpow2,ni,pos,u,minrp,nf)
!$omp& shared(maxrp,vdon,nze)
!$omp& shared(formcut,gamma,pol)
!$omp$ shared(psi_t,maxpsi_t,e_t,la_t,na_t,l_t,nnt,npos_t)
!$omp$ shared(zasym,alkali,ubb_max3,ubb_max1,arho)
!$omp& private(temp3,mintemp3,maxtemp3)
!$omp& private(vmatt,temp,vmati,lm)

      do nchi = nchii, nchif
         gpunum=mod(omp_get_thread_num(),ngpus)
!         gpunum=1
#ifdef GPU
         call acc_set_device_num(gpunum,acc_device_nvidia)
#endif

!         call getchinfo (nchi,nt, lg, psii, maxpsii, ei, lia, nia, li)
         psii(:)=psi_t(:,nchi)
         maxpsii=maxpsi_t(nchi)
         ei=e_t(nchi)
         lia=la_t(nchi)
         nia=na_t(nchi)
         li=l_t(nchi) 

         nposi=0
         posi = positron(nia,lia,nposi)
         nqmi = npk(nchi+1) - npk(nchi)

         if (li.gt.ldw) then
C  UI has to be added to the kinetic energy operator K in order to
C  to extract the one-electron energy k**2/2. Associated with K is V
C  which contains VDCORE.
            ui(:) = 0.0
            do i = minvdc, maxvdc
               ui(i) = - vdcore(i,min(li,lamax)) * 2.0
            enddo
            ud(:) = 0.0
         else
            do i = 1, meshr
               ui(i) = dwpot(i,nchi)
               ud(i) = ui(i)
            enddo 
         endif             

       allocate(temp3(1:meshr,nchi:nchtop))
       allocate(vmati(1:kmax,1:nqmi,nchi:nchtop))
!$acc data create(vmati(1:kmax,1:nqmi,nchi:nchtop))

C$OMP PARALLEL DO DEFAULT(PRIVATE) num_threads(nnt)
C$OMP& SCHEDULE(dynamic)
C$OMP& SHARED(vdcore,npk,meshr,minvdc,maxvdc,dwpot,nchi,nchtop,lg)
C$OMP& SHARED(ldw,rmesh,temp3,mintemp3,maxtemp3,rpow1,rpow2,nznuc,nze)
C$OMP& SHARED(rnorm,nqmi,li,lia,psii,minrp,maxrp,u,maxpsii,ltmin,ctemp)
C$OMP& SHARED(pos,ni,nf)
C$OMP& SHARED(psi_t,maxpsi_t,e_t,la_t,na_t,l_t,uf)
C$OMP& SHARED(zasym,alkali,ubb_max3,ubb_max1,arho)
        do nchf = nchi, nchtop
!            call getchinfo (nchf,nt,lg, psif, maxpsif, ef, lfa, nfa,lf)
            ef=e_t(nchf)
            lfa=la_t(nchf)
            nfa=na_t(nchf)
            lf=l_t(nchf)

!            nposf=0
!            posf = positron(nfa,lfa,nposf)
            nqmf = npk(nchf+1) - npk(nchf)
            if(pos(nchf).neqv.pos(nchi)) cycle

            if (lf.gt.ldw) then
C  UF has to be added to the kinetic energy operator K in order
C  to extract the one-electron energy k**2/2. Associated with K is V
C  which contains VDCORE.
               uf(:,nchf)=0.0
               do i = minvdc, maxvdc
                  uf(i,nchf) = - vdcore(i,min(lf,lamax)) * 2.0
               enddo
            else 
               do i = 1, meshr
                  uf(i,nchf) = dwpot(i,nchf)
               enddo
            endif                

            nchtopf = nchtop
            nchtopi = nchtop
      minfun = 1
!      maxfun = min(maxpsii,maxpsif)
      maxfun = min(maxpsi_t(nchi),maxpsi_t(nchf))
      do i = minfun, maxfun
!         fun(i) = psii(i) * psif(i) * rmesh(i,3)
         fun(i) = psi_t(i,nchi)*psi_t(i,nchf)*rmesh(i,3)
      end do 

      c1tmp=0
      c2tmp=0
      c3tmp=0

      do i = 1, meshr
         temp3(i,nchf) = 0.0
      enddo
      mini = maxr
      maxi = 1
      ltmin(nchf) = 1000000
      ctemp(nchf) = 0.0
      do 10 ilt = -lia, lia, 2
         lt = lfa + ilt
         if (lt.lt.0.or.lt.gt.ltmax) go to 10
         call cleb(2*lia,2*lt,2*lfa,0,0,0,c1)
         c1tmp = cgc0(float(lia),float(lt),float(lfa))
         if (abs((c1-c1tmp)/(c1tmp+1e-20)).gt.1e-3) then
            print*,'CGCs 1 do not agree:',c1, c1tmp,lia,lt,lfa
            stop 'CGCs 1 do not agree'
         endif 
         call cleb(2*li,2*lf,2*lt,0,0,0,c2)
         c2tmp = cgc0(float(li),float(lf),float(lt))         
         if (abs((c2-c2tmp)/(c2tmp+1e-20)).gt.1e-3) then
            c2tmp2 = cgcigor(2*li,2*lf,2*lt,0,0,0)
C$OMP critical(print)
            print*,'CGCs 2 do not agree:',c2, c2tmp,c2tmp2,li,lf,lt
C$OMP end critical(print)
            c2 = c2tmp2
         endif 
         call rac7(2*li,2*lf,2*lia,2*lfa,2*lt,2*lg,c3)
         c3tmp = cof6j(float(li),float(lf),float(lt),float(lfa),
     >      float(lia),float(lg))*(-1)**(li+lf+lia+lfa)
         if (abs((c3-c3tmp)/(c3tmp+1e-6)).gt.1e-2) then
C  below does happen for NPAR = 1, but not significant
C$OMP critical(print)
            print*,'WARNING: CJ6 and W do not agree in D:',c3, c3tmp,
     >         li,lf,lia,lfa,lt,lg
C$OMP end critical(print)
            c3 = c3tmp
c$$$            stop 'CJ6 and W do not agree in D'
         endif 
         
         c = c1 * c2 * c3
         if (abs(c).lt.1e-10) go to 10
!          const = (-nze)*(-1)**(lg + lfa) * hat(li) *
         const = (-1)**(lg + lfa) * hat(li) *
     >      hat(lf) * hat(lia) / hat(lt) * c

         call form(fun,minfun,maxfun,rpow1(1,lt),rpow2(1,lt),
     >      minrp(lt),maxrp(lt),meshr,temp,i1,i2)

C  Subtract 1/r, but only for same atom-atom channels when lambda = 0
         if (lt.eq.0.and..not.pos(nchf).and.ni.eq.nf) then
          call nuclear(fun,.not.pos(nchf),minfun,maxfun,i2,u,nznuc,temp)
         endif

       if(pos(nchf)) then 
C     ANDREY: Hydrogen: ssalling + interpolation:
            if (.not.alkali) then
C  The factor of two below is the reduced mass. We are working with
C  matrix
C  elements whose channels are multiplied by sqrt of the reduced mass.
C  Here
C  we have positronium-positronium matrix element, hence a factor of 2
C  overall.
c$$$            const = const * (1-(-1)**lt)
              const = 2.0 * const * (1-(-1)**lt) * (zasym+1.0) !latter for He+
c$$$               const = - const !Rav's derivation, but not Alisher's
               do i = i1, i2
                  xin(i-i1+1) = rmesh(i,1)
                  yin(i-i1+1) = temp(i) * xin(i-i1+1) ** (lt+1)
                  xout(i-i1+1) = rmesh(i,1) * 2d0
               enddo
               if (i2-i1+1.gt.maxr) then
                  print*,'I2,I1,LT:',i2,i1,lt
                  stop 'problem in call to intrpl'
               endif
               call intrpl(i2-i1+1,xin,yin,i2-i1+1,xout,yout)
c$$$            nnn = nnn + 1
               do i = i1, i2
                  if (xout(i-i1+1).gt.xin(i2-i1+1))
     >                 yout(i-i1+1) = yin(i2-i1+1) 
                  temp(i) = 2.0 * (yout(i-i1+1) / xout(i-i1+1)**(lt+1))
                  
c$$$               write(50+nnn,'(1p,4e10.3,i5)') xout(i-i1+1),
c$$$     >            yout(i-i1+1) / xout(i-i1+1)**(lt+1),
c$$$     >            xin(i-i1+1), yin(i-i1+1) / xin(i-i1+1)**(lt+1), i
               enddo
            else ! alkali
C     ANDREY: pos-alkali case: radial integrals with Ubb
               if (lt.gt.ubb_max3) stop
     $              'stopped: vmat: Ubb with larger lt are needed'
               const=2d0*const
C              some constants from get_ubb:
               rlgmin = log10(rmesh(1,1))
               rlgmax = log10(rmesh(meshr,1))
               drlg = (rlgmax-rlgmin)/real(ubb_max1-1)

c---- find irmin and irmax needed for interpolation
               irmin = 1;
               rmin = rmesh(minfun, 1);
               ir = 2;
               do while (arho(ir).lt.rmin)
                  ir = ir+1
               end do
               irmin = ir-1
C     ----------------------------------------------
               rmax = rmesh(maxfun, 1)
               irmax = ubb_max1; ir = irmax;
               do while (arho(ir).gt.rmax)
                  ir = ir-1
               end do
               irmax = ir+1
               maxir = irmax-irmin+1
C---------------------------------------------------

               do i = 1, i2
                  rho = rmesh(i,1);
                  call get_ubb(minfun,maxfun,irmin,irmax,maxir,drlg,lt,
     $                 rho,ubbb)
                  temp(i) = 0.0
                  do j = minfun, maxfun
                     r = rmesh(j,1)
c     explicit integration of Ubb (Eq. 23):
C                    call Ubb0(lt, r, rho, res0, res)
C     Ubb integrals computed outside:
C                    res = ubb_res(j,i,lt)
C     interpolation: minus: alisher uses (1-(-1)**lt):
                     res = -ubbb(j)/max(rho,r/2.)
                     temp(i) = temp(i)+res*fun(j)
                  end do        ! j
               end do           ! i
            end if              ! alkali or hydrogen
         else                   ! pos
C  Multiply by -1 for positron scattering as then NZE = 1, for non
C  pos-pos
C  channels
            const = - nze * const
         endif ! pos

         do i = i1, i2
            temp3(i,nchf) = const * temp(i) + temp3(i,nchf)
         enddo
         mini = min(mini,i1)
         maxi = max(maxi,i2)
C  The variable CTEMP is the asymptotic value of const * temp * r ** (lt+1)
c$$$         if ((lt.eq.1.or.lt.eq.2).and.i2.eq.meshr) then
            if (i2.eq.meshr.and.lt.ne.0.and.lt.lt.ltmin(nchf)) then
               ctemp(nchf) = rnorm * const * temp(i2)/rpow2(i2,lt)
               ltmin(nchf) = lt
         endif 
 10   continue
      if (ltmin(nchf).lt.10.and.maxi.eq.meshr) then
         ctemp(nchf) = rnorm * temp3(maxi,nchf)/rpow2(maxi,ltmin(nchf))
      endif 
      if (nchi.eq.nchf.and.nqmi.eq.1.and.u(1).eq.0.0) then
C  Define the channel dependent distorting potential when running the
C  Born case. The units are Rydbergs.
         do i = 1, meshr
            dwpot(i,nchi) = 0.0
         enddo
         do i = mini, maxi
            dwpot(i,nchi) = temp3(i,nchf) * 2.0
         enddo
      endif
      
C  As both CHII and CHIF contain the integration weights, we divide TEMP by   
C  them.
      do i = mini, maxi
         temp3(i,nchf) = rnorm * temp3(i,nchf) / rmesh(i,3)
      end do
      mintemp3(nchf) = mini
      maxtemp3(nchf) = maxi
      end do !nchf
!$omp end parallel do

      mini1 = mintemp3(nchtop)

! create an array of pos to have all posf and 
      call gpuvdirect(maxr,meshr,rmesh,kmax,nqmi,nchi,nchtop,npk,
     >     mintemp3,maxtemp3,temp3,ltmin,minchil,chil,ctemp,itail,trat,
     >     nchan,vmati,childim,ngpus,nchii,second,pos)


! !$acc wait
! !$omp critical

!!$omp parallel do num_threads(2) schedule(dynamic)
      do nchf = nchi, nchtop
        if(pos(nchi).eqv.pos(nchf)) then
         vdon(nchf,nchi,0:1) = vdon(nchf,nchi,0:1) + vmati(1,1,nchf)
         vdon(nchi,nchf,0:1) = vdon(nchf,nchi,0:1)
!      end do
!!$omp end parallel do

C  Define the direct on shell V matrix used for analytic Born subtraction
C  tail integrals still need to be incorporated

!!$omp parallel do num_threads(2) schedule(dynamic)
!!$omp& shared(vdon,nchi,nchtop,vmatt,vmati,npk,nsmax,nqmi)
!!$omp& private(ns,ki,kf,nqmf)
!      do nchf=nchi,nchtop
         nqmf = npk(nchf+1) - npk(nchf)
         do ns = 0, nsmax
            do ki = 1, nqmi
               do kf = 1, nqmf
                  vmatt(kf,ki,ns,nchf) = vmati(kf,ki,nchf) 
               enddo
            enddo
         enddo
        else
               do ns = 0, nsmax
                do ki = 1, nqmi
                  do kf = 1, nqmf
                    vmatt(kf,ki,ns,nchf) = 0.0
                  enddo
                enddo
               enddo
               nqmf = npk(nchf+1) - npk(nchf)
               ef=e_t(nchf)
               lfa=la_t(nchf)
               nfa=na_t(nchf)
               lf=l_t(nchf)
               nposf=npos_t(nchf)
               lm=min0(Lf+lfa+lia,Li+lfa+lia) ! additional argument
               if (pos(nchf)) then
                call posvmat(nqmi,lia,li,nia,gk(1,nchi),
     >               npk(nchtop+1)-1,nchi,nqmf,lfa,lf,nposf,nfa,
     >               gk(1,nchf),nchf,lg,npk,etot,vmatt(1,1,0,nchf),lm)
               else
                call posvmat(nqmf,lfa,lf,nfa,gk(1,nchf),
     >               npk(nchtop+1)-1,nchf,nqmi,lia,li,nposi,nia,
     >               gk(1,nchi),nchi,lg,npk,etot,vmatt(1,1,0,nchf),lm)

!posvmat
               endif
        endif
      end do
!!$omp end parallel do
      call clock(s2)
      td = td + s2 - s1
 
!      do nchf = nchi, nchtop
!        nqmf = npk(nchf+1) - npk(nchf)
!          psif(:)=psi_t(:,nchf)
!          maxpsif=maxpsi_t(nchf)
!          lfa=la_t(nchf)
!          lf=l_t(nchf)
        
C Define exchange terms if IFIRST = 1
      if (ifirst.eq.1) then
         if (ispeed.ne.3) then
!            call makev3e(psii,maxpsii,lia,nchi,psif,maxpsif,lfa,
!     >         nchf,li,lf,chil(1,npk(nchi),1),minchil(npk(nchi),1),nqmi,
!     >         chil(1,npk(nchf),1),minchil(npk(nchf),1),nqmf,lg,rnorm,
!     >         second,npk,vmatt,nchtop)
!            call makev3e(chil(:,:,1),psii,maxpsii,lia,nchi,psi_t,


            call makev3e(chil,psii,maxpsii,lia,nchi,psi_t,
     >         maxpsi_t,la_t,li,l_t,minchil(npk(nchi),1),nqmi,
     >         lg,rnorm,second,npk,vmatt,nchtop,childim,nnt,gpunum)

            call clock(s3)
            te1 = te1 + s3 - s2
         endif 
      endif
!      end do
   
      do nchf = nchi, nchtop
        nqmf = npk(nchf+1) - npk(nchf)
!         call getchinfo (nchf,nt,lg, psif, maxpsif, ef, lfa, nfa, lf)
          psif(:)=psi_t(:,nchf)
          maxpsif=maxpsi_t(nchf)
          ef=e_t(nchf)
          lfa=la_t(nchf)
          lf=l_t(nchf)

      if (ifirst.eq.1) then
C  Define energy dependent exchange terms
         if (ispeed.ne.2) call makev1e(nqmi,psii,maxpsii,ei,lia,
     >      li,chil(1,npk(nchi),1),minchil(npk(nchi),1),gk(1,nchi),
     >      npk(nchtop+1)-1,etot,theta,0,nqmf,psif,maxpsif,
     >      ef,lfa,lf,chil(1,npk(nchf),1),minchil(npk(nchf),1),
     >      gk(1,nchf),npk(nchtop+1)-1,lg,rnorm,
     >      uf,ui,nchf,nchi,nold,nznuc,npk,ve2ee,vmatt,nchtop)
         call clock(s4)
         te2 = te2 + s4 - s3
C  End of exchange
      end if            
      end do
    
      do 200 nchf = nchi, nchtop 
         nqmf = npk(nchf+1) - npk(nchf)
            if (npk(2)-npk(1).eq.1.or.
     >         (.not.scalapack.and.nodeid.eq.1)) then
               do ki = 1, nqmi
                  kii = npk(nchi) + ki - 1
                  do kf = 1, nqmf
                     kff = npk(nchf) + kf - 1
                     if (kff.ge.kii) then
                        vmat(kff,kii)=vmat(kff,kii)+vmatt(kf,ki,0,nchf)
                        if (nsmax.eq.1) vmat(kii,kff+1) =
     >                     vmat(kii,kff+1)+vmatt(kf,ki,1,nchf)
                     endif
                  enddo
               enddo
            else 
               if (nchf.le.nchistop(nodeid)) then
                  do ki = 1, nqmi
                     kii = npk(nchi) + ki - 1
                     do kf = 1, nqmf
                        kff = npk(nchf) + kf - 1
                        if (kff.ge.kii) then
                           vmat01(kff,kii) = vmat01(kff,kii)
     >                        + vmatt(kf,ki,0,nchf)
                           if (nsmax.eq.1) vmat01(kii,kff+1) = 
     >                        vmat01(kii,kff+1) + vmatt(kf,ki,1,nchf)
                        endif
                     enddo
                  enddo
               else
                  do ki = 1, nqmi
                     kii = npk(nchi) + ki - 1
                     do kf = 1, nqmf
                        kff = npk(nchf) + kf - 1
                        vmat0(kff,kii) = vmat0(kff,kii)
     >                     + vmatt(kf,ki,0,nchf)
                        if (nsmax.eq.1) vmat1(kii,kff+1) =
     >                     vmat1(kii,kff+1) + vmatt(kf,ki,1,nchf)
                     enddo
                  enddo
               endif
            endif 

 200     continue

 
  
C  End of NCHI loop

! !$omp end critical

! !$omp end single

!$acc end data 
!      deallocate(vmati,temp3)
      deallocate(vmati,temp3)

      end do

!$omp end parallel do

! !$omp end parallel


CC GPU ONLY
#ifdef GPU
      do gpunum=0,ngpus-1
         call acc_set_device_num(gpunum,acc_device_nvidia)

!$acc exit data delete(chil(1:meshr,1:(npk(nchtop+1)-1),1:childim)) 
!$acc& delete(nchtop,npk(1:nchtop+1))
! !$acc& delete(vmati(1:kmax,1:kmax,nchii:nchtop))
! !$acc& delete(chitemp(1:meshr,1:kmax))
!$acc& finalize

      enddo
#endif

!      deallocate(temp3)
!      deallocate(vmati)
!      deallocate(chitemp)
      deallocate(nchansi,nchansf)
!      deallocate(temp3) 
      return
      end


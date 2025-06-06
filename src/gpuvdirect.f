      subroutine gpuvdirect(maxr,meshr,rmesh,kmax,nqmi,nchi,nchtop,npk,
     >   mintemp3,maxtemp3,temp3,!ltmin,minchilx,chilx,ctemp,itail,trat,
     >   nchan,nqmfmax,vmatt,ichildim_dum,ngpus,nnt,nchii_dummy,second,
     >   maxpsii,temp2,nqmi2,nchtop2,nsmax,ifirst)
      use chil_module
#ifdef GPU_ACC
      use openacc
#endif
      integer npk(nchtop+1)
      integer nqmf,kff,kii,maxi
      integer mintemp3(nchan),maxtemp3(nchan),ltmin(nchan)
      real rmesh(maxr,3),ctemp(nchan),temp(maxr)
c$$$      real tmp3(1:meshr,nchi:nchtop)
      real temp3(1:meshr,nchi:nchtop)
      real vmatt(nqmfmax,nqmfmax,nchi:nchtop,0:nsmax)
      real dotp1,dotp2
!      real,allocatable :: chitemp(:,:)!,tmp(:,:)
      real chitemp(meshr,nqmi)
      real tmp(nqmi,nqmfmax) 
      logical second
      integer maxpsii
      real temp2(1:maxpsii,1:nqmi2,nchi:nchtop2)

!      allocate(chitemp(meshr,nqmi))
      maxi2 = maxpsii
c$$$      if(nsmax.ne.1) then
c$$$        maxi2=0
c$$$      endif

      kii=0
      kff=0
      nqmf=0
      maxi=0
      dotp1=0.0
      dotp2=0.0
        
c$$$      temp3(1:meshr,nchi:nchtop)=tmp3(1:meshr,nchi:nchtop)

!#ifndef GPU_ACC
!#ifndef GPU_OMP
c$$$!$omp parallel num_threads(nnt) !default(private)
c$$$!$omp& private(gpunum,nchf,nqmf,maxi,mini,chitemp,ki,kf,i,kff,kii,tmp,
c$$$!$omp& tnum)
c$$$!$omp& shared(nchi,nchtop,npk,maxtemp3,temp3,chil,vmatt,ngpus,temp2)
c$$$!$omp&shared(nqmi,maxi2,nsmax,nnt,meshr,maxpsii,nchtop2,nqmfmax,minchil)
c$$$!$omp do schedule(dynamic)
!#endif
!#endif

#ifdef GPU_ACC
!$acc data if(nqmi>100) 
!$acc& copy(vmatt(1:nqmfmax,1:nqmi,nchi:nchtop,0:nsmax))
!$acc& present(npk(1:nchtop+1))
!$acc& present(chil(1:meshr,npkstart:npkstop,1))
!$acc& present(minchil(npkstart:npkstop,1))
!$acc& present(nchtop)
!$acc& copyin(nqmi,maxtemp3,temp3(1:meshr,nchi:nchtop))
!$acc& create(chitemp)
!$acc& create(tmp)
!$acc& copyin(temp2) !copyin(temp2(1:maxpsii,1:nqmi,nchi:nchtop)) 

#endif
#ifdef GPU_OMP
!$omp target data if(nqmi>100)
!$omp& map(to:temp3,temp2,maxtemp3)
!$omp& map(alloc:tmp,chitemp)
!$omp& map(tofrom:vmatt(1:nqmfmax,1:nqmi,nchi:nchtop,0:nsmax))
#endif
      do nchf = nchi, nchtop
         nqmf = npk(nchf+1) - npk(nchf)
         maxi = min(maxtemp3(nchf),meshr)
#ifdef GPU_ACC
!!!$acc parallel if(nqmi>100)
!$acc kernels if(nqmi>100)         
!$acc loop independent collapse(2)
#endif
#ifdef GPU_OMP
!$omp target teams distribute parallel do collapse(2) if(nqmi>100)
#endif
#ifndef GPU_ACC
#ifndef GPU_OMP
!$omp parallel do default(shared) private(ki,i,kii,minki) 
!$omp& collapse(1) schedule(static) !not dynamic!
#endif
#endif
         do ki = 1, nqmi
            kii = ki+npk(nchi)-1
            minki = minchil(kii,1)
            do i = 1, minki - 1 !not necessary since below: mini = max(minchil(kff,1),minchil(ki+npk(nchi)-1,1))
               chitemp(i,ki) = 0.0
            enddo
            do i = minki, maxi
               chitemp(i,ki) = temp3(i,nchf) * chil(i,kii,1)
            enddo
         enddo
#ifndef GPU_ACC
#ifndef GPU_OMP
!$omp end parallel do
#endif
#endif
#ifdef GPU_OMP
!$omp end target teams distribute parallel do
#endif
         if (ifirst.eq.1) then
#ifdef GPU_ACC                 
!$acc loop independent collapse(2)
#endif
#ifdef GPU_OMP
!$omp target teams distribute parallel do collapse(2) if(nqmi>100)
#endif
#ifndef GPU_ACC
#ifndef GPU_OMP
!$omp parallel do default(shared) private(kf,ki,mini,kff)
!$omp& collapse(2) schedule(static) !dynamic is ok too
#endif
#endif
            do ki = 1, nqmi
               do kf = 1, nqmf
                  kff = npk(nchf) + kf - 1
                  mini = minchil(kff,1)
                  tmp(ki,kf) = dot_product(chil(mini:maxi2,kff,1),
     >               temp2(mini:maxi2,ki,nchf))
               enddo
            enddo
#ifndef GPU_ACC
#ifndef GPU_OMP
!omp end parallel do
#endif
#endif
#ifdef GPU_OMP
!$omp end target teams distribute parallel do
#endif
         else
            tmp(1:nqmi,1:nqmf) = 0.0
         endif
#ifdef GPU_ACC
!$acc loop independent collapse(2)
#endif
#ifdef GPU_OMP
!$omp target teams distribute parallel do collapse(2) if(nqmi>100)
#endif
#ifndef GPU_ACC
#ifndef GPU_OMP
!$omp parallel do default(shared) private(kff,ki,kf,mini)
!$omp& collapse(2) schedule(static) !dynamic is ok too
#endif
#endif
         do ki = 1, nqmi
            do kf = 1, nqmf
               kff = npk(nchf) + kf - 1
               mini = max(minchil(kff,1),minchil(ki+npk(nchi)-1,1))
               vmatt(kf,ki,nchf,0) = tmp(ki,kf) + vmatt(kf,ki,nchf,0) +
     >            dot_product(chil(mini:maxi,kff,1),
     >            chitemp(mini:maxi,ki)) !the slow part
               if (nsmax.eq.1)  !does not slow things down
     >            vmatt(kf,ki,nchf,1)=vmatt(kf,ki,nchf,0)-2.0*tmp(ki,kf)
            end do
         end do
#ifndef GPU_ACC
#ifndef GPU_OMP
!$omp end parallel do
#endif
#endif
#ifdef GPU_ACC
!!!$acc end parallel
!$acc end kernels         
#endif
#ifdef GPU_OMP
!$omp end target teams distribute parallel do
#endif
       end do !nchf
#ifdef GPU_OMP
!$omp end target data
#endif

#ifndef GPU_ACC
#ifndef GPU_OMP
c$$$!$omp end do
c$$$!$omp end parallel
#endif
#endif

#ifdef GPU_ACC
!$acc wait
!$acc end data 
#endif

      return
      end


      subroutine makev3e(chilx,psii,maxpsii,lia,nchi,psif,maxpsif,lfa,
     >   li,lf,minchilx,nqmi,lg,rnorm,second,npk,nchtop,nnt,temp2)
      use chil_module
      include 'par.f'
      integer nnt
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension npk(nchtop+1),fun(maxr)
c$$$      dimension chil(meshr,npk(nchtop+1)-1)
c$$$      dimension minchil(npk(nchtop+1)-1)
      dimension psii(maxr),
     >   psif(maxr,nchtop),const(-lamax:lamax,nchtop)
      dimension maxpsif(nchtop), lfa(nchtop), lf(nchtop)
      real temp2(maxpsii,nqmi,nchi:nchtop) !temp2(meshr,nqmi,nchi:nchtop)
      real, allocatable :: temp(:)
!      real vmatt(1:kmax,1:kmax,0:1,1:nchtop)
!      real vmatt(nqmfmax,nqmfmax,nchi:nchtop,0:1)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      logical second
c
      common /di_el_core_polarization/ gamma, r0, pol(maxr)
      common/smallr/ formcut
      integer maxpsii
c      
      hat(l)= sqrt(2.0 * l + 1.0)

      allocate(temp(maxr))

!$omp parallel do default(private) num_threads(nnt) !collapse(2)
!$omp& schedule(dynamic)
!$omp& shared(rnorm,const,lf,lfa,li,lia,lg,nchi,nchtop)
      do nchf=nchi,nchtop
      do ilt = -lia, lia, 2
         const(ilt,nchf) = 0.0
         lt = lf(nchf) + ilt
         if (lt.ge.0.and.lt.le.ltmax) then
            call cleb(2*li,2*lt,2*lfa(nchf),0,0,0,c1)
            c1tmp = cgc0(float(li),float(lt),float(lfa(nchf)))
            if (abs((c1-c1tmp)/(c1tmp+1e-20)).gt.1e-3) then
               print*,'CGCs 1 do not agree:',c1, c1tmp,li,lt,lfa(nchf)
               stop 'CGCs 1 do not agree'
            endif 
            call cleb(2*lia,2*lf(nchf),2*lt,0,0,0,c2)
            c2tmp = cgc0(float(lia),float(lf(nchf)),float(lt))
            if (abs((c2-c2tmp)/(c2tmp+1e-20)).gt.1e-3) then
               c2tmp2 = cgcigor(2*lia,2*lf(nchf),2*lt,0,0,0)
! C$OMP critical(print)
               print*,'CGCs 2 do not agree:',c2, c2tmp,c2tmp2,lia,
     >                                      lf(nchf),lt
! C$OMP end critical(print)
               c2=c2tmp2
            endif 
            call rac7(2*lia,2*lf(nchf),2*li,2*lfa(nchf),2*lt,2*lg,c3)   
            c3tmp = cof6j(float(lia),float(lf(nchf)),float(lt),
     >              float(lfa(nchf)),
     >         float(li),float(lg))*(-1)**(lia+lf(nchf)+li+lfa(nchf))
            if (abs((c3-c3tmp)/(c3tmp+1e-6)).gt.1e-2) then
! C$OMP critical(print)
               print*,'WARNING:CJ6 and W do not agree in E:',c3, c3tmp
! C$OMP end critical(print)
               c3 = c3tmp
c$$$               stop 'CJ6 and W do not agree'
            endif 
            c = c1 * c2 * c3
            const(ilt,nchf) = rnorm * (-1)**(lf(nchf)) * hat(li) *
     >         hat(lf(nchf)) * hat(lia) / hat(lt) * c
         endif
      enddo ! ilt
      end do
!$omp end parallel do

!$omp parallel do default(private) num_threads(nnt) !collapse(2)
!$omp& schedule(dynamic)
!$omp& shared(chil,psif,psii,npk,lia,const,rpow1,meshr,nchi,lf,
!$omp& rpow2,temp2,nqmi,minchil,maxpsii,gamma,pol,minrp,maxrp,maxpsif,
!$omp& nchtop)
      do nchf=nchi,nchtop
      do ki = 1, nqmi
         kii = npk(nchi) + ki - 1
         minfun = minchil(kii,1)
         maxfun = maxpsif(nchf)

         do i = minfun, maxfun
            fun(i) = chil(i,kii,1) * psif(i,nchf)
         end do
         do i = 1, maxpsii !meshr
            temp2(i,ki,nchf) = 0.0
         enddo
         mini = meshr
         do 15 ilt = -lia, lia, 2
            lt = lf(nchf) + ilt
            if (lt.lt.0.or.lt.gt.ltmax.or.abs(const(ilt,nchf)).lt.1e-10)
     >         go to 15
            call form(fun,minfun,maxfun,rpow1(1,lt),
     >      rpow2(1,lt),minrp(lt),maxrp(lt),maxpsii,temp,i1,i2)
c
       if(lt.eq.1.and.gamma.ne.0.0) then
          sum1 = 0.0
          do i=minfun,maxfun
             sum1 = sum1 + fun(i)*pol(i)
          end do
          tmp = gamma * sum1
          do i=i1,i2
             temp(i) = temp(i) - tmp*pol(i)
          end do
       end if
       if (i2.gt.maxpsii) stop 'something is wrong'
            do i = i1, i2
               temp2(i,ki,nchf) = const(ilt,nchf) * temp(i) 
     >                          + temp2(i,ki,nchf)
            enddo
            mini = min(mini,i1)
 15      continue
         do i = mini, maxpsii
            temp2(i,ki,nchf) = temp2(i,ki,nchf) * psii(i)
         enddo
      enddo
      enddo
!$omp end parallel do

      deallocate(temp)
      end

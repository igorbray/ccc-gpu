      subroutine gpuvdirect(maxr,meshr,rmesh,kmax,nqmi,nchi,nchtop,npk,
     >     mintemp3,maxtemp3,temp3,ltmin,minchil,chil,ctemp,itail,trat,
     >     nchan,nqmfmax,vmatt,childim,ngpus,nnt,nchii,second)
#ifdef GPU
      use openacc
#endif
      integer npk(nchtop+1)
      integer childim
      integer mintemp3(nchan),maxtemp3(nchan),ltmin(nchan),
     >     minchil(npk(nchtop+1)-1)
      real chil(1:meshr,1:(npk(nchtop+1)-1))
      real rmesh(maxr,3),ctemp(nchan),temp(maxr)
      real temp3(1:meshr,nchi:nchtop)
      real vmatt(nqmfmax,nqmfmax,nchi:nchtop)
      allocatable :: chitemp(:,:)
      logical second
      integer gpunum,tnum
      integer, external :: omp_get_thread_num

      allocate(chitemp(meshr,nqmi))
#ifdef GPU
!$omp parallel num_threads(ngpus)
#else
!$omp parallel num_threads(nnt)
#endif
!$omp& private(gpunum,nchf,nqmf,maxi,mini,chitemp,ki,kf,i,kff,kii)
!$omp& shared(nchi,nchtop,npk,maxtemp3,temp3,chil,vmatt,ngpus)
!$omp& shared(nqmi)

#ifdef GPU
      tnum=omp_get_thread_num()
      gpunum=mod(tnum,ngpus)
      call acc_set_device_num(gpunum,acc_device_nvidia)
#endif

!$acc data 
!$acc& copyin(vmatt(1:nqmfmax,1:nqmi,nchi:nchtop))
!$acc& present(npk(1:nchtop+1))
!$acc& present(chil(1:meshr,1:(npk(nchtop+1)-1)))
!$acc& present(minchil(1:npk(nchtop+1)-1))
!$acc& present(nchtop)
!$acc& copyin(nqmi,maxtemp3,temp3(1:meshr,nchi:nchtop))
!$acc& create(chitemp)
!$omp do schedule(dynamic)
      do nchf = nchi, nchtop
         nqmf = npk(nchf+1) - npk(nchf)
         maxi = maxtemp3(nchf)
!$acc kernels async(nchf)
!$acc loop independent collapse(2)
         do ki = 1, nqmi
            do i = 1, maxi !minchil(ki+npk(nchi)-1), maxi !minki, maxi
               chitemp(i,ki) = temp3(i,nchf) * chil(i,ki+npk(nchi)-1)
            enddo
         enddo
!$acc loop independent collapse(2)
         do ki = 1, nqmi
            do kf=1,nqmf
               kii = npk(nchi) + ki - 1
               kff = npk(nchf) + kf - 1
c$$$               if (kff.ge.kii) then
                  mini = minchil(kff)
                  vmatt(kf,ki,nchf)=vmatt(kf,ki,nchf)+dot_product(
     >                 chil(mini:maxi,kff),
     >                 chitemp(mini:maxi,ki))
c$$$     >                 chil(minchil(kff):maxi,kff),
c$$$     >                 chitemp(minchil(kff):maxi,ki))
c$$$               endif
            end do
         end do
!$acc end kernels
!$acc update self(vmatt(1:nqmf,1:nqmi,nchf)) async(nchf)
!!$acc update self(vmati(1:nqmf,1:nqmi,nchf)) async(nchf)
       end do
!$omp end do

! START
! !$acc kernels
! !$acc loop independent collapse(2)
!         do ki = 1, nqmi
!            do kf = 1, nqmf
!               dotprod=0.0
!               do it = mini,maxi
!                chitemp=temp3(it,nchf) * chil(it,ki+npk(nchi)-1)
!                dotprod=dotprod+
!     >          (chil(it,kf+npk(nchf)-1)*chitemp)
!               end do
!               vmati(kf,ki,nchf)=dotprod
!            end do
!         end do
! !$acc end kernels
! END

c         if (itail.ne.0.and.ctemp(nchf).ne.0.0) then
c            lt = ltmin(nchf)
c            const = ctemp(nchf) * (trat / rmesh(meshr,1))**lt
c            do ki = 1, nqmi
c               do i = minchil(ki+npk(nchi)-1,2), meshr
c                  temp(i) = chil(i,ki+npk(nchi)-1,2)
c     >                 /rmesh(i,3)/rmesh(i,1)**(lt+1)
c               enddo 
c               do kf = 1, nqmf
c                  mini = max(
c     >              minchil(ki+npk(nchi)-1,2),minchil(kf+npk(nchf)-1,2))
c                  vmati(kf,ki,nchf)=vmati(kf,ki,nchf)+const*dot_product(
c     >               chil(mini:meshr,kf+npk(nchf)-1,2),temp(mini:meshr))
c               enddo
c            enddo
c         endif !itail
! check the following works with Na
c        vdon(nchf,nchi,0:1) = vdon(nchf,nchi,0:1) + vmati(1,1,nchf)
c         vdon(nchi,nchf,0:1) = vdon(nchf,nchi,0:1)
c$$$         vdon(nchf,nchi,0) = vdon(nchf,nchi,0) + vmati(1,1,nchf)
c$$$         vdon(nchi,nchf,0) = vdon(nchf,nchi,0)
!      enddo !nchf

!$acc wait
!$acc end data 
!$omp end parallel
      deallocate(chitemp)

      return
      end


      subroutine makev3e(chil,psii,maxpsii,lia,nchi,psif,maxpsif,lfa,
     >   li,lf,minchil,nqmi,lg,rnorm,second,npk,
     >   nqmfmax,vmatt,nchtop,nnt,ngpus)
#ifdef GPU
      use openacc
#endif
      include 'par.f'
      integer nnt,gpunum
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension npk(nchtop+1),fun(maxr)
      dimension chil(meshr,npk(nchtop+1)-1)
      dimension minchil(npk(nchtop+1)-1),
     >   const(-lamax:lamax,nchtop)
      dimension psii(maxr),
     >   psif(maxr,nchtop)
      dimension maxpsif(nchtop), lfa(nchtop), lf(nchtop)
      real, allocatable :: temp2(:,:,:)
      real, allocatable :: temp(:)
!      real vmatt(1:kmax,1:kmax,0:1,1:nchtop)
      real vmatt(nqmfmax,nqmfmax,nchi:nchtop,0:1)
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
      integer maxpsii,maxi
      integer, external :: omp_get_thread_num

c      
      hat(l)= sqrt(2.0 * l + 1.0)

#ifdef GPU
      do gpunum=0,ngpus-1
         call acc_set_device_num(gpunum,acc_device_nvidia)
!$acc enter data copyin(vmatt(1:nqmfmax,1:nqmi,nchi:nchtop,0:1)) async
      end do
#endif

      allocate(temp2(meshr,nqmi,nchtop))
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

      maxi=maxpsii

!$omp parallel do default(private) num_threads(nnt) !collapse(2)
!$omp& schedule(dynamic)
!$omp& shared(chil,psif,psii,npk,lia,const,rpow1,meshr,nchi,lf,
!$omp& rpow2,temp2,nqmi,minchil,maxi,gamma,pol,minrp,maxrp,maxpsif,
!$omp& nchtop,gpunum)
      do nchf=nchi,nchtop
      do ki = 1, nqmi
         kii = npk(nchi) + ki - 1
         minfun = minchil(kii)
         maxfun = maxpsif(nchf)

         do i = minfun, maxfun
            fun(i) = chil(i,kii) * psif(i,nchf)
         end do
         do i = 1, meshr
            temp2(i,ki,nchf) = 0.0
         enddo
         mini = meshr
         do 15 ilt = -lia, lia, 2
            lt = lf(nchf) + ilt
            if (lt.lt.0.or.lt.gt.ltmax.or.abs(const(ilt,nchf)).lt.1e-10)
     >         go to 15
            call form(fun,minfun,maxfun,rpow1(1,lt),
     >      rpow2(1,lt),minrp(lt),maxrp(lt),maxi,temp,i1,i2)
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

            do i = i1, i2
               temp2(i,ki,nchf) = const(ilt,nchf) * temp(i) 
     >                          + temp2(i,ki,nchf)
            enddo
            mini = min(mini,i1)
 15      continue
         do i = mini, maxi
            temp2(i,ki,nchf) = temp2(i,ki,nchf) * psii(i)
         enddo
      enddo
      enddo
!$omp end parallel do

! What is the purpose of the following loop?
#ifdef GPU
      do gpunum=0,ngpus-1
         call acc_set_device_num(gpunum,acc_device_nvidia)
!$acc wait
      end do
#endif

      call gpuvexchange(chil,temp2,nqmfmax,vmatt,ngpus,nnt,nchi,nchtop,
     >                        nqmi,minchil,npk,maxi,kmax)

#ifdef GPU
      do gpunum=0,ngpus-1
         call acc_set_device_num(gpunum,acc_device_nvidia)
!$acc exit data delete(vmatt(nqmfmax,nqmi,nchi:nchtop,0:1))
!$acc& delete(temp2(1:maxi,1:nqmi,nchi:nchtop))
      enddo
#endif

      deallocate(temp2)
      deallocate(temp)
      end


      subroutine gpuvexchange(chil,temp2, nqmfmax, vmatt, ngpus, nnt,
     >                        nchi,nchtop,nqmi,minchil,npk,maxi,kmax)
#ifdef GPU
      use openacc
#endif
      integer childim
      common/meshrr/ meshr
      dimension npk(nchtop+1)
      dimension chil(meshr,npk(nchtop+1)-1),minchil(npk(nchtop+1)-1)
      real temp2(meshr,nqmi,nchtop)
      real vmatt(nqmfmax,nqmfmax,nchi:nchtop,0:1)
      integer maxi,npktmp
      integer gpunum,tnum,ngpus
      integer nchf,nchtop,nchi,ki,kf,kff
      real tmp
      integer, external :: omp_get_thread_num

#ifdef GPU
!$omp parallel default(none) num_threads(ngpus)
#else
!$omp parallel default(none) num_threads(nnt)
#endif
!$omp& private(mini,gpunum,tnum,nchf,nqmf,ki,kf,kff,tmp)
!$omp& shared(nchi,nchtop,npk,temp2,chil,vmatt,ngpus,maxpsii)
!$omp& shared(nqmi,meshr,kmax,nqmfmax,minchil)
!$omp& firstprivate(maxi)

#ifdef GPU
      tnum=omp_get_thread_num()
      gpunum=mod(tnum,ngpus)
      call acc_set_device_num(gpunum,acc_device_nvidia)
#endif

!$acc data present(chil(1:meshr,1:npk(nchtop+1)-1))
!$acc& present(minchil(1:npk(nchtop+1)-1))
!$acc& copyin(temp2(1:maxi,1:nqmi,nchi:nchtop))
!$acc& present(vmatt(1:nqmfmax,1:nqmi,nchi:nchtop,0:1))
!$acc& present(npk(1:nchtop+1))

!$omp do schedule(dynamic)
      do nchf=nchi,nchtop
      nqmf=npk(nchf+1) - npk(nchf)
!$acc kernels async(nchf)
!$acc loop independent collapse(2)
       do ki = 1, nqmi
         do kf = 1, nqmf
            kff = npk(nchf) + kf - 1
            mini = minchil(kff)
            tmp = dot_product(chil(mini:maxi,kff)
     >            ,temp2(mini:maxi,ki,nchf))
            vmatt(kf,ki,nchf,0) = vmatt(kf,ki,nchf,0) + tmp
            vmatt(kf,ki,nchf,1) = vmatt(kf,ki,nchf,1) - tmp
        end do
       end do
!$acc end kernels
!$acc update self(vmatt(1:nqmf,1:nqmi,nchf,0:1)) async(nchf)
      end do
!$omp end do

!$acc wait
!$acc end data

!$omp end parallel

      end

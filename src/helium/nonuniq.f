      subroutine nonuniq(KJ,Tspin,nchm,iacc)
      use CI_MODULE
      use DM_MODULE
      include 'par.f'

      integer, intent(in):: KJ, nchm, iacc
      real, intent(in):: Tspin

      integer nspm,lo,ko,nset
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      integer la,sa,lpar,np
      common /helium/ la(KNM), sa(KNM), lpar(KNM), np(KNM)
      integer na, nam
      common /CIdata/ na(nspmCI,KNM), nam(KNM)

      double precision COF6J
      real temp(maxr)
      real*8 tmp, Ang1
!
!      
!     find the dimension of the DD matrix and set up arrays get_nch, get_na   and    get_nsp
      nchnsp = 0
      do nch = 1, nchm
         call getchinfo (nch, NN, KJ, temp, maxt, 
     >        eatmp, latmp, natmp, Lnch)
         do nsp=1,nspm
            if(lo(nsp) .eq. Lnch) then
               nchnsp = nchnsp + 1
               print*,'nchnsp,nch,nsp,Lnch', nchnsp,nch,nsp,Lnch
            endif
         enddo
      enddo
      nchnsp_max = nchnsp

      if(allocated(DD)) then
         deallocate(DD,get_nch,get_N,get_nsp)
      endif
      if(nchnsp_max .gt. 0) then
         allocate(DD(nchnsp_max,nchnsp_max),get_nch(nchnsp_max),
     >        get_N(nchnsp_max),get_nsp(nchnsp_max))
      else
         return
      endif
      nchnsp = 0
      do nch = 1, nchm
         call getchinfo(nch, NN, KJ, temp, maxt, 
     >        eatmp, latmp, natmp, Lnch)
         do nsp=1,nspm
            if(lo(nsp) .eq. Lnch) then
               nchnsp = nchnsp + 1
               get_nch(nchnsp) = nch
               get_N(nchnsp) = NN
               get_nsp(nchnsp) = nsp
            endif
         enddo
      enddo
!  finsh setting up arrays get_nch, get_na   and    get_nsp
!      print*,'nchnsp_max =', nchnsp_max
!  Start setting up array DD
      do nchnspi = 1, nchnsp_max
         nchi = get_nch(nchnspi)
         Ni = get_N(nchnspi)
         nspi = get_nsp(nchnspi)         
         linsp=lo(nspi)
         rlinsp=linsp
         Li = linsp
         do nchnspf = nchnspi, nchnsp_max
            nchf = get_nch(nchnspf)
            Nf = get_N(nchnspf)
            nspf = get_nsp(nchnspf)
            lfnsp=lo(nspf)
            rlfnsp=lfnsp
            Lf = lfnsp
           
!            print*,'Nf,Lf,nfnsp,Ni,Li,ninsp:',Nf,Lf,nfnsp,Ni,Li,ninsp
!            print*,'sa(Nf),sa(Ni):',sa(Nf),sa(Ni)
            rK=KJ
            lli=la(Ni)
            llf=la(Nf)
            lsi=sa(Ni)
            lsf=sa(Nf)
            rLi=Li
            rlli=lli
            rLf=Lf
            rllf=llf
            rsi=lsi
            rsf=lsf

!            do j=1,nam(Ni)
!            print*,j,na(j,Ni),(C(Ni,j,i),i=1,nam(Ni))
!            enddo

            Ang1 = 0.0
            do jni1=1,nam(Ni)
!               print*,'i:',jni1,na(jni1,Ni),nspf
               if(nspf .ne. na(jni1,Ni)) cycle
               do jnf1=1,nam(Nf)
!                  print*,'f:',jnf1,na(jnf1,Nf),nspi
                  if(nspi .ne. na(jnf1,Nf)) cycle
                  do jnf2=1,nam(Nf)
                     nf2 = na(jnf2,Nf)
!                     print*,'C(Nf,jnf1,jnf2)=',C(Nf,jnf1,jnf2)
                     if(C(Nf,jnf1,jnf2).eq.0.0D0) cycle
                     lf2=lo(nf2)
                     rlf2=lf2
                     do jni2=1,nam(Ni)
                        ni2 = na(jni2,Ni)
                        if(ni2 .ne. nf2) cycle ! assume orthogonal set
                        li2=lo(ni2)
                        rli2=li2                        
!                        print*,'C(Ni,jni1,jni2)=',C(Ni,jni1,jni2)
                        if(C(Ni,jni1,jni2).eq.0d0) cycle
                        tmp = C(Ni,jni1,jni2)*C(Nf,jnf1,jnf2)*
     >                       COF6J(rLi,rli2,rllf,rLf,rK,rlli)
!                        print*,'tmp=',tmp
                        Ang1 = Ang1 + tmp
                     enddo
                  enddo
               enddo
            enddo
            const = 2.0 * COF6J(0.5,0.5,rsf,0.5,Tspin,rsi)*
     >                       (-1)**(Lf+Li-lli+llf+1-lsi+lsf)*
     >                       sqrt(dble((2*lli+1)*(2*llf+1)*
     >                       (2*lsi+1)*(2*lsf+1)))
            
!            print*,'f,i,DD(f,i):',nchnspf,nchnspi,Ang1 * const
            DD(nchnspf,nchnspi) = Ang1 * const
            DD(nchnspi,nchnspf) = Ang1 * const
            
         enddo
      enddo


      call modifyDM(nchnsp_max,DD,iacc)

      return 
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine modifyDM(nchnsp_max,DD,iacc)
      integer:: nchnsp_max
      real*8, dimension(nchnsp_max,nchnsp_max) :: DD
      integer:: iacc
      real*8, dimension(nchnsp_max):: w,fv1,fv2
      real*8, dimension(nchnsp_max,nchnsp_max) :: z
      real*8:: tmp, acc
      integer:: matz,m,n,ierr
      integer, dimension(nchnsp_max):: ipoint
      integer:: ii, ipont_max, lam

      print*, 'Enter modify DD()'
      if(nchnsp_max .le. 0) return

      matz = 1

!      write(*,'(100I5)') (m, m=1,nchnsp_max)
!      write(*,'(100I5)') (get_nch(m), m=1,nchnsp_max)
!      write(*,'(100I5)') (get_nsp(m), m=1,nchnsp_max)

      do n=1,nchnsp_max
         write(*,'(100F12.5)') (DD(n,m), m=1,nchnsp_max)
      enddo

      call rs(nchnsp_max,nchnsp_max,DD,w,matz,z,fv1,fv2,ierr)
      write(*,'("ierr =",I3)') ierr
      print*,'eigenvalues: N=', nchnsp_max
      write(*,'(5F13.8)') (real(w(i)), i=1,nchnsp_max)
      print*
      do n=1,nchnsp_max
         write(*,'(100F12.5)') (z(n,m), m=1,nchnsp_max)
      enddo

      if(iacc .ge. 0) then
         acc = 10.0**(-iacc)
      else
         acc = -1.0
      endif
      print*,'nonuniq.f: acc=', acc

      ii = 0
      do m=1,nchnsp_max
         if(abs(w(m)-1.d0) .ge. acc) then
            ii = ii + 1
            ipoint(ii) = m
         endif
      enddo
      ipoint_max = ii

!      ipoint(2) = 3
!      ipoint(3) = 2
c
      do m=1,nchnsp_max
         do n=m,nchnsp_max
            tmp =0d0
            do ii=1,ipoint_max
               lam = ipoint(ii)
               tmp = tmp + z(m,lam)*w(lam)*z(n,lam)
            enddo
            DD(m,n) = tmp
            DD(n,m) = tmp
         enddo
      enddo
      
      print*
      print*,'D-matrix after modification'
      print*,'Number of eigenvalues is N=', nchnsp_max
      print*,'Number of excluded eigenvalues (=1) is:', 
     >     nchnsp_max-ipoint_max
      do n=1,nchnsp_max
         write(*,'(100F12.5)') (DD(n,m), m=1,nchnsp_max)
      enddo


!      stop
      return
      end 
!-------------------------------------------------------------------------------------
      function  iget_nchnsp(nch,nsp,nchnsp_max,get_nch,get_nsp)
      integer::  iget_nchnsp
      integer, intent(in):: nch,nsp,nchnsp_max
      integer, dimension(nchnsp_max), intent(in):: get_nch,get_nsp
      integer:: i
      
      iget_nchnsp = -1
      
      do i=1,nchnsp_max
         if(nch .eq. get_nch(i) .and. nsp .eq. get_nsp(i)) then
            iget_nchnsp = i
            exit
         endif         
      enddo
      if(iget_nchnsp .eq. -1) then
         print*,'nonuniq.f:iget_nchnsp = -1, nch,nsp:',nch,nsp
         stop
      endif

      return
      end function  iget_nchnsp
      
!------------------------------------------------------------------------------------------

      subroutine v_nu(nchm,nspm,npk,nchi,Ni,Li,si,nqmi,
     >     nchf,Nf,Lf,sf,nqmf,
     >     ortchil,Etot,KJ,na,nam,namax,lo,
     >     nchnsp_max,DD,get_nch,get_nsp,VEE)

      include 'par.f'

      integer,intent(in):: nchm,nspm
      integer, dimension(nchm+1),intent(in):: npk
      integer,intent(in):: nchi,Ni,Li,si,nqmi
      integer,intent(in):: nchf,Nf,Lf,sf,nqmf
      real, dimension(npk(nchm+1)-1,nspm),intent(in):: ortchil
      double precision,intent(in)::  Etot
      integer,intent(in):: KJ
      integer, dimension(nspmCI,KNM),intent(in):: na 
      integer, dimension(KNM),intent(in):: nam
      integer, dimension(nspmax),intent(in):: lo
      integer,intent(in):: nchnsp_max
      real*8, dimension(nchnsp_max,nchnsp_max),intent(in):: DD
      integer,dimension(nchnsp_max),intent(in):: get_nch,get_nsp
      real, dimension(npk(nchm+1)-1,npk(nchm+1)),intent(inout)::  VEE

      real*8:: pnorm, const
      real*8, dimension(nqmf,nqmi):: tmp
      integer:: loimax,Ke2,kf,kff,ki,kii
      integer:: jnf1,nf1,lf1,jni1,ni1,li1
      integer:: nchnspi,nchnspf
      real*8:: R0,R1

c     Ni,Nf  -   numbers of physical chanels
c      print*,'    '
c      print*,'Nf,Ni:',Nf,Ni

c
c     Lf=li1, therefore  Lf  could not be bigger then lomax. It makes
c     limitations on value of KJ: |KJ-llf|<=Lf<=KJ+llf (0<=llf<=2*lomax),
c     hence maximum value of KJ is: (the same as in the routine ve2me())
      loimax = 0
      do ni1 = 1, nspm
         loimax = max(loimax,lo(ni1))
      enddo
      Ke2 = 3*loimax
      if(KJ.gt.Ke2)  return
      if (Lf.gt.loimax) return
      if (Li.gt.loimax) return

      pnorm = 2.0d0/dacos(-1.0d0)
      const = Etot*pnorm 

      tmp(:,:) = 0d0
c     this is 1 s.p. loop on s.p. functions of the coordinate r1
      do jnf1=1,nam(Nf)
         nf1 = na(jnf1,Nf)
         lf1=lo(nf1)         
!         print '(A20,100I5)',
!     >        'jnf1,nf1,lf1,Li:',jnf1,nf1,lf1,Li,nam(Nf)
         if(lf1 .ne. Li) cycle
         nchnspi = iget_nchnsp(nchi,nf1,
     >        nchnsp_max,get_nch,get_nsp)
c     this is 1 s.p. loop on s.p. functions of the coordinate r0
         do jni1=1,nam(Ni)
            ni1 = na(jni1,Ni)
            li1=lo(ni1)       
!            print '(A20,100I5)',
!     >           'jni1,ni1,li1,Lf:',jni1,ni1,li1,Lf,nam(Ni)
            if(li1 .ne. Lf) cycle           
            nchnspf = iget_nchnsp(nchf,ni1,
     >           nchnsp_max,get_nch,get_nsp)
!            print*,'nchnspf,nchnspi:',nchnspf,nchnspi
            do kf = 1, nqmf
               kff = npk(nchf) + kf - 1
               R0 = ortchil(kff,ni1)
               do ki = 1, nqmi
                  kii = npk(nchi) + ki - 1
                  if (kff.lt.kii) cycle
                  R1 = ortchil(kii,nf1)   
!                  print*, kf,ki,R0,R1, DD(nchnspf,nchnspi)
                  tmp(kf,ki) = tmp(kf,ki) + DD(nchnspf,nchnspi)*R0*R1
               enddo            !     end ki loop
            enddo               !     end kf loop
            
         enddo                  !     end ni1 loop
      enddo                     !     end nf1 loop
      
      
      
      do kf = 1, nqmf
         kff = npk(nchf) + kf - 1
         do ki = 1, nqmi
            kii = npk(nchi) + ki - 1
            if (kff.lt.kii) cycle
            VEE(kff,kii) = VEE(kff,kii) + tmp(kf,ki) * const
c            VEE(kii,kff+1) = VEE(kii,kff+1) + veT * const2
         enddo
!         print '(100E15.8)',(VEE(kff,kii),kii=1,kff)
      enddo

c      print*,'********** Finish v_nu **********'
      return
      end subroutine v_nu        




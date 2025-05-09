c-----------------------------------------------------------------
c**************** hel(ll,ls) *************************************
c-----------------------------------------------------------------
c     here the set of s.p. functions "kl" consist of two sets:
c     a number of first  s.p. functions "kl" correspond to the alpha= alpha2,
c     while all other to the  alpha= alpha1.
c
c     it is regulated by the papameter "ndif" in input file F5 (if ndif=0, then
c     all  s.p. functions belong to the set with  alpha= alpha1. in other case
c     the functions that belong to the set with  alpha= alpha2 defined by
c     the parameters "ldif" and "nkd(l)"
c
c
c     subroutine inputdata 
c     subroutine checklomax(l1max,l2max)
c     subroutine config(ll,ls,lparity)
c     subroutine hammat   <--   call eigv(ncm,H,bb)
c     function sumlambda(n1,n2,n1p,n2p)
c     function RAD(lam,n1,n2,n1p,n2p)  <--  call form(fun,minfun,maxfun,
c     >   rpow1(1,lam),rpow2(1,lam),minrp(lam),maxrp(lam),meshr,temp,il,i2)
c     function fangint(lam,l1,l1p,l2,l2p)
c     subroutine setslatintfang
c     double precision function oneel(n1,n2,n1p,n2p)
c     - all one electron calculations in the file setmax.f
c     - Wigner coefficients - from civ3 in the file wigner.f
c     - subroutine eigv(n,H,bb)  - in the file eigv.f
c     - subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr) -from netlib
c                                                         in file rsg.f
c     
c     
c     
c     this subroutine is called to make CI calculation of helium states
c     with particular ll,ls,lparity.

      subroutine hel(enionry,ll,ls,lparity,
     >     CI_top,w_top,fl,maxf,minf,iter)

      include 'par.f'
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      dimension  nk1(0:lomax), nk2(0:lomax)  
      common/orbc/ncm,no1(ncmCI),no2(ncmCI)
      double precision CI_top(ncmCI,ncmCI), w_top(ncmCI)
      double precision, dimension(:,:), allocatable::  H, bb, CI
      double precision, dimension(:), allocatable::  w
      common /switch_manifold/  isw_manifold
      double precision, dimension(:,:), allocatable::  H_set, 
     >     bb_set, CI_set
      double precision, dimension(:), allocatable::  w_set, tmp_set

      if(iter.eq.0) then
         call  inputdata(ll,ls,lparity,l1max,l2max,nk1,nk2)
         call set_basis_optim(ll,ls,lparity,nk1)
         call config(ll,ls,lparity,l1max,l2max,nk1,nk2)
      end if
      if (ncm.eq.0) then
         print*,'WARNING, NCM=0'
         return
      endif 
      
      if(allocated(H)) then
         deallocate(H,bb,CI,w)
      else
         allocate(H(ncm,ncm),bb(ncm,ncm),CI(ncm,ncm),w(ncm))
      endif

      call hammat(H,bb,ll,ls,lparity,fl,maxf,minf)
c     start diagonalization
      call eigv(enionry,ncm,H,bb,CI,w,iter)
c     finish diagonalization


      CI_top(1:ncm,1:ncm) = CI(1:ncm,1:ncm)
      w_top(1:ncm) = w(1:ncm)


      if(isw_manifold .eq. 2 .and. 
     >     ll.eq.0 .and. ls.eq.0 .and. lparity.eq.1) then
         isw_old = isw_manifold
         isw_manifold = 0
         call hammat(H,bb,ll,ls,lparity,fl,maxf,minf)
c     start diagonalization
         call eigv(enionry,ncm,H,bb,CI,w,iter)
c     finish diagonalization
         isw_manifold =  isw_old
         
         CI_top(1:ncm,1) = CI(1:ncm,1)
         w_top(1) = w(1)

      endif


!  replace positive energy states with states that have no interaction between manifolds.
!  this is done by perfoming diagonalisation in  a set of target states (with energy larger then e_set) 
!  with intereaction betweenmanifolds removed.
                              
      if(isw_manifold .eq. 3) then   
         
         isw_old = isw_manifold
         isw_manifold = 1

         call hammat(H,bb,ll,ls,lparity,fl,maxf,minf)
! find how many states to be included here:
! set energy  e_set  to be defined later as ionization energy, 
! states with energyes above will be included in this rediagonalisation.
         e_set = enionry/2.0  ! (two-electron) 
         do n=1,ncm
            if(w_top(n) .gt. e_set) then
               n_set = n
               exit
            endif
         enddo
         nsc = ncm - n_set + 1

         print*, '>>> transforming positive energy states: nsc=', nsc

         if(allocated(H_set)) then
            deallocate(H_set,bb_set,CI_set,w_set,tmp_set)
         else
            allocate(H_set(nsc,nsc),bb_set(nsc,nsc),CI_set(nsc,nsc),
     >           w_set(nsc),tmp_set(ncm))
            bb_set(:,:) = 0d0
            H_set(:,:) = 0d0
            CI_set(:,:) = 0d0
            w_set(:) = 0d0
            tmp_set(:) = 0d0
         endif
         

         do n=1,nsc
            bb_set(n,n) = 1d0  ! states form orthonomal set
         enddo
         do n=1,nsc
            n_top = n + n_set - 1
            do np=1,n
               np_top = np + n_set - 1

               do ipr=1,ncm
                  tmp_set(ipr) = SUM( CI(1:ncm,n_top) * H(1:ncm,ipr) )
               enddo
               H_set(n,np) = SUM( CI(1:ncm,np_top) * tmp_set(1:ncm) )
               H_set(np,n) = H_set(n,np)               

            enddo
         enddo         

!         print*, '>>> H_set(1,1)=', H_set(1,1),enionry

c     start diagonalization         
         call eigv(enionry,nsc,H_set,bb_set,CI_set,w_set,iter)
c     finish diagonalization        
         
         do n=1,nsc
            n_top = n + n_set - 1
            do i=1,ncm
               CI_top(i,n_top) = SUM(CI_set(1:nsc,n)*CI(i,n_set:ncm))
            enddo
            w_top(n_top) = w_set(n)
         enddo

         deallocate(H_set,bb_set,CI_set,w_set,tmp_set)


         isw_manifold =  isw_old
         print*,'finish isw_manifold'
      endif

      deallocate(H,bb,CI,w)

      return
      end
c-----------------------------------------------------------------
c****************  Read input data  ******************************
c-----------------------------------------------------------------
c     This subroutine read input data form file F5, check them with par.f,
c     make factorials call factorials, and then return
c     to call subroutine config(ll,ls,lparity,l1max,l2max,nk1,nk2) 

      subroutine inputdata(ll,ls,lparity,l1max,l2max,nk1,nk2)
      use vmat_module, only: nodeid
      include 'par.f'
      dimension  nk1(0:lomax), nk2(0:lomax)  
      common /nstatearray/ nstate(0:lomax,0:1,2)
      character newsym*20
c      common /cut_off_energies/ en_max, global_en_max, sym_en_max



      nk1(:) = 0
      nk2(:) = 0
      
 1    read(3,'(a20)') newsym
c 1    read(3,FMT=*,ERR=20) newsym, sym_en_max
c      en_max = sym_en_max
c      go to 20
c 20   en_max = global_en_max
      read(3,*) ll,ls,lparity      
      ip = 1
      if(-(-1)**ll.eq.lparity) ip = 2
      read(3,*) l1max
      read(3,*) (nk1(l), l=0,l1max)
      read(3,*) l2max
      read(3,*) (nk2(l2), l2=0,l2max)

      if(nstate(ll,ls,ip).le.0) go to 1

      if (nodeid.eq.1) then
      write(4,'(a20)') newsym
      if (nodeid.eq.1) print*, 'Calculating target symmetry: ', newsym
      write(4,'("l =",I3,",  s =",I3,",  parity =",I3)')
     >   ll,ls,lparity
      write(4,'("l1max =",I5)') l1max
      write(4,'("nk1(l) =",10I5)') (nk1(l), l=0,l1max)
      write(4,'("l2max =",I5)') l2max
      write(4,'("nk2(l) =",10I5)') (nk2(l), l=0,l2max)
      endif
      call checklomax(l1max,l2max)
      return
      end
c-------------------------------------------------------------------
c     This is section to check consistency of the input data with par.f

      subroutine checklomax(l1max,l2max)

      include 'par.f'
      if(lomax.lt.l1max) then
         write(*,'("lomax<l1max, increase lomax in par.f")')
      end if
      if(lomax.lt.l2max) then
         write(*,'("lomax<l2max, increase lomax in par.f")')
      end if
      return
      end
c-----------------------------------------------------------------
c************** Set list of configurations  **********************
c-----------------------------------------------------------------
      subroutine config(ll,ls,lparity,l1max,l2max,nk1,nk2)
      use vmat_module, only: nodeid
      include 'par.f'
      dimension   nk1(0:lomax), nk2(0:lomax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common/orbc/ncm,no1(ncmCI),no2(ncmCI)
      common /corearray/ nicm, ncore(nspmCI)
      common /corearray_lsp/ nicm_lsp(0:lamax,0:1,-1:1), 
     >     ncore_lsp(nspmCI,0:lamax,0:1,-1:1)
      common /ionic_configurations/ l_ion_core, nk_ion_core(0:lomax)
      common /ngive_sym_opt/  ngivesym(nspmax), N_opt
      common /ngive_ionic_orb/  ngive_ion_orb(nspmax)
      common /include_opt/  include_opt(nspmCI,nspmCI)

      if (nodeid.eq.1)
     >   write(4,'("number of s.p. orbitals nspm=",I4)') nspm
!      print*, "l_ion_core:", l_ion_core,", nk_ion_core:",
!     >     (nk_ion_core(l2), l2=0,l2max)
!      print*, " nk2:", (nk2(l1), l1=0,l1max)
      
      nk2(l2max+1:) = 0
      nk1(l1max+1:) = 0

      nco=0

      do nsp2=1,nspm
         do nsp1=1,nspm
c     If nset(nsp)=-1 then this orbital is not used in the CI for this symmetry.
c            if(nset(nsp1).gt.0.and.nset(nsp2).gt.0) then

c     If N_opt.eq.0 - no opt.orb. have been built and this is standard CI run.
c     If N_opt.gt.0 - some opt.orb. have been built. In this case nsp2 will go
c     only over ionic core orbitals due to "if (l2.le.l2max.and.k2.le.nk2(l2)) then"
c     statement later on. Array include_opt(...) have information on which
c     s.p.orbital from the old set should be excluded for given ionic core ic=nsp2.
            if(N_opt.eq.0.or.
     >              (N_opt.gt.0.and.include_opt(nsp1,nsp2).ne.-1)) then
            l1=lo(nsp1)
            l2=lo(nsp2)
            k1=ko(nsp1)
            k2=ko(nsp2)
            if(l1.le.l1max.and.(k1.le.nk1(l1).or.nset(nsp1).eq.2)) then
c     for a given ionic core orbital (nsp2-->k2,l2) the configuration with outer 
c     electron  in (nsp1-->k1,l1) is included if: 
c         1. l_ion_core.lt.0  , all config. may be included
c         2. l_ion_core.ge.0.and.k1.le.nk2(l1), these are correlation config.
c         3. l_ion_core.ge.0.and.nset(nsp1).eq.2, these are optimized orbital config.
c         4. l_ion_core.ge.0.and.k2.le.nk_ion_core(l2) , continuum like config.
c                with large k1 are allowed for given inner electron orbital

               if(l_ion_core.lt.0.or.(l_ion_core.ge.0.and.
     >              (k1.le.nk2(l1)
c     >              (k1.le.nk2(l1).or.nset(nsp1).eq.2
     >              .or.k2.le.nk_ion_core(l2)) )) then

            if (l2.le.l2max.and.k2.le.nk2(l2)) then
               if(lparity.eq.(-1)**(l1+l2)) then
                  if(abs(l1-l2).le.ll.and.ll.le.l1+l2) then
                     lspar = 1
                     if(nsp1.eq.nsp2) then
                        lspar=(-1)**(ll+ls)
                     end if
                     if(lspar.eq.1) then
                        icheck = 0
                        do ncx=1,nco
                           if(no1(ncx).eq.nsp1.and.no2(ncx).eq.nsp2.or.
     >                        no1(ncx).eq.nsp2.and.no2(ncx).eq.nsp1)
     >                        then
                              icheck = 1
                           end if
                        end do
                        if(icheck.eq.0) then
                           nco=nco+1
      if(nco.gt.ncmCI) then         
         write(*,'("Stop in config(): check par.f: nco>ncmCI",2i5)')
     >        nco, ncmCI
c$$$         write(4,'("Stop in config(): check par.f: nco>ncmCI,2i5")')
c$$$     >        nco, ncmCI
         stop
      end if
                           no1(nco)=nsp1
                           no2(nco)=nsp2
c                   write(*,'("nco,no1(nco),no2(nco), l1,l2: ",8I5)')
c     >                        nco,no1(nco),no2(nco),lo(no1(nco)),
c     >                        lo(no2(nco)),nset(no1(nco)),
c     >                        nset(no2(nco)),ngivesym(no1(nco))
!                           write(*,'(3i5,"   ",2i5)') 
!     >                          k1,l1,nset(nsp1),k2,l2

                        end if
                     end if
                  end if
               end if
            end if
               end if
            end if
            end if
         end do                 ! end nsp1 loop
      end do                    ! end nsp2 loop
      ncm=nco
      if (nodeid.eq.1)
     >write(4,'("number of configurations ncm=",I5)') ncm
c     find core orbitals: they are in general given by array nk2(l2);
c     but not all of them will be used (due to selection rules). Therefore
c     it is better to find core orbitals after list of configurations is set.
c     This means that core orbitals is now correspond to array no2(nco) 
c     (because index nsp2 is placed in the outer loop and nsp1 in the inner loop)
      do nco=1,ncm
c         print*, 'nco,no1(nco),no2(nco), l1,l2: ',nco,no1(nco),no2(nco),
c     >        lo(no1(nco)),lo(no2(nco)),nset(no1(nco)),nset(no2(nco)), 
c     >        ngivesym(no1(nco))
         icheck = 0
         do ni=1,nicm
            if(ncore(ni).eq.no2(nco)) then
               icheck = 1
            end if
         end do
         if(icheck.eq.0) then
            nicm = nicm + 1
            ncore(nicm) = no2(nco)
c            print*, '      nicm,ncore(nicm)', nicm,ncore(nicm)
         end if

         icheck = 0
         do ni=1,nicm_lsp(ll,ls,lparity)
            if(ncore_lsp(ni,ll,ls,lparity).eq.no2(nco)) then
               icheck = 1
            end if
         end do
         if(icheck.eq.0) then
            nicm_lsp(ll,ls,lparity) = nicm_lsp(ll,ls,lparity) + 1
            n_t = nicm_lsp(ll,ls,lparity)
            ncore_lsp(n_t,ll,ls,lparity) = no2(nco)
c            print*, 'nicm_lsp,ncore_lsp', nicm_lsp(ll,ls,lparity),
c     >           ncore_lsp(n_t,ll,ls,lparity)
         end if
      end do
      return
      end
c-----------------------------------------------------------------
**** Calculation two-electron and one electron matrix elements ***
c-----------------------------------------------------------------
      subroutine hammat(H,bb,ll,ls,lparity,fl,maxf,minf)
      use vmat_module, only: nodeid
      include 'par.f'
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      common/orbc/ncm,no1(ncmCI),no2(ncmCI)
      double precision H(ncm,ncm), bb(ncm,ncm)
      double precision  H1,H2,H3,H4,b1,b2,b3,b4
      double precision  isym,isymp,sn,snp
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /switch_manifold/  isw_manifold

      if (nodeid.eq.1)
     >write(4,'("enter hammat")')
      do nc=1,ncm
         do ncp=1,nc
            n1=no1(nc)
            n2=no2(nc)
            n1p=no1(ncp)
            n2p=no2(ncp)
            l1=lo(n1)
            l2=lo(n2)
            l1p=lo(n1p)
            l2p=lo(n2p)
            isym= DBLE((-1)**(l1+l2-ll-ls))
            isymp= DBLE((-1)**(l1p+l2p-ll-ls))
            sn=dsqrt(2.D0)
            snp=sn
            isw = 0
            if(isw_manifold .eq. 1 .or. isw_manifold .eq. 2) then
               if(n2.ne.n2p .or. (n2.eq.n2p .and. l1.ne.l1p) ) then
                  isw = 1
               endif
            endif
            call hamconf(ll,l1,l2,l1p,l2p,n1,n2,n1p,n2p,
     >         H1,H2,H3,H4,b1,b2,b3,b4,sn,snp,fl,maxf,minf,isw)
            H(nc,ncp) = H1 + dble(isym)*H2 + dble(isymp)*H3 + 
     >           dble(isym*isymp)*H4
            bb(nc,ncp) = b1 + dble(isym)*b2 + dble(isymp)*b3 + 
     >           dble(isym*isymp)*b4
            H(nc,ncp) = H(nc,ncp)/(2.D0*sn*snp)
            bb(nc,ncp) = bb(nc,ncp)/(2.D0*sn*snp)
            H(ncp,nc) = H(nc,ncp)
            bb(ncp,nc) = bb(nc,ncp)
         end do
      end do
      return
      end
c----------------------------------------------------------------------
      subroutine hamconf(ll,l1,l2,l1p,l2p,n1,n2,n1p,n2p,
     >   H1,H2,H3,H4,b1,b2,b3,b4,sn,snp,fl,maxf,minf,isw)

      include 'par.f'
      double precision H11,H12,H1,H21,H22,H2,H31,H32,H3,H41,H42,H4
      double precision b1,b2,b3,b4
      double precision oneel,sumlambda,bmat,sn,snp
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)

      H12 = 0.0
      if(isw .eq. 0) then
         H12=sumlambda(ll,l1,l2,l1p,l2p,fl(1,n1),fl(1,n2),fl(1,n1p),
     >        fl(1,n2p),minf(n1),minf(n2),minf(n1p),minf(n2p),
     >        maxf(n1),maxf(n2),maxf(n1p),maxf(n2p))
      endif
      H11 = oneel(n1,n2,n1p,n2p)
      b1 = bmat(n1,n2,n1p,n2p)
      H1 = H11 + H12
      H22 = H12
      H21 = H11
      b2 = b1
c     if n1=n2, or n1p=n2p then all two-electron integrals are equal
c     hence only if all n1.ne.n2, and n1p.ne.n2p   we need calculate
c     one more two-electron integral (H12=H42, H22=H32  -  anyway)
      if(n1.ne.n2) then
         sn=1.D0
         H21 = oneel(n2,n1,n1p,n2p)
         if(n1p.ne.n2p) then
            if(isw .eq. 0) then
               H22 = sumlambda(ll,l2,l1,l1p,l2p,fl(1,n2),fl(1,n1),
     >              fl(1,n1p),fl(1,n2p),minf(n2),minf(n1),minf(n1p),
     >              minf(n2p),maxf(n2),maxf(n1),maxf(n1p),maxf(n2p))
            endif
            b2 = bmat(n2,n1,n1p,n2p)
         end if
      end if
      H2 = H21 + H22
      H31 = H11
      H41 = H21
      H32 = H22
      H42 = H12
      b3 = b2
      b4 = b1
      if(n1p.ne.n2p) then
         snp=1.D0
         H31 = oneel(n1,n2,n2p,n1p)
         H41 = H31
         if(n1.ne.n2) then
            H41 = oneel(n2,n1,n2p,n1p)
         end if 
      end if
      H3 = H31 + H32
      H4 = H41 + H42
      return
      end
c-----------------------------------------------------------------
c************** Calculation two-electron matrix elements *********
c-----------------------------------------------------------------
      double precision function sumlambda(ll,l1,l2,l1p,l2p,
     >   f1,f2,f1p,f2p,
     >   min1,min2,min1p,min2p,max1,max2,max1p,max2p)

      include 'par.f'
      dimension f1(nmaxr),f2(nmaxr),f1p(nmaxr),f2p(nmaxr)
      double precision sumlam, flam, Rlam, fangint, RAD
      lammin=max(abs(l1-l1p),abs(l2-l2p))
      lammax=min(l1+l1p,l2+l2p)
      sumlam=0.D0
      do lam=lammin,lammax
         Rlam=0.0D0
         flam = fangint(lam,l1,l1p,l2,l2p,ll)
         if(flam.ne.0.0D0) then
            Rlam = RAD(lam,f1,f1p,f2,f2p,min1,min1p,min2,min2p,
     >         max1,max1p,max2,max2p)
            sumlam = sumlam + flam * Rlam
         end if
      end do
      sumlambda = sumlam
      return
      end
c-----------------------------------------------------------
      double precision function RAD(lam,f1,f1p,f2,f2p,
     >   min1,min1p,min2,min2p,max1,max1p,max2,max2p)

      include 'par.f'
      double precision  tmp, sum1, sum2
      common /meshrr/ nr,gridr(nmaxr,3)
      dimension f1(nmaxr),f2(nmaxr),f1p(nmaxr),f2p(nmaxr)
      dimension temp(nmaxr), fun(nmaxr)
      common/powers/ rpow1(nmaxr,0:ltmax),rpow2(nmaxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(nmaxr,0:lmax)
      common /di_el_core_polarization/ gamma, r0, pol(nmaxr)
      minfun=max(min2,min2p)
      maxfun=min(max2,max2p)
      do i=minfun,maxfun
         fun(i) = f2(i)*f2p(i)*gridr(i,3)
      end do
      maxfm = min(max1,max1p)
      call form(fun,minfun,maxfun,rpow1(1,lam),rpow2(1,lam),
     >   minrp(lam),maxrp(lam),maxfm,temp,i1,i2)
      mini = max(i1,min1,min1p)
      maxi = min(i2,max1,max1p)
      tmp = 0d0
      do i=mini,maxi
         tmp = tmp + dble(temp(i)*f1(i)*f1p(i)*gridr(i,3))
      end do

      if(lam.eq.1.and.gamma.ne.0.0) then
         sum1 = 0.0
         sum2 = 0.0
         do i=max(min1,min1p),min(max1,max1p)
            dr = gridr(i,3)
            sum1 = sum1 + dble(f1(i)*f1p(i)*dr*pol(i))
         end do
         do i=minfun,maxfun
            sum2 = sum2 + dble(fun(i)*pol(i))
         end do

c     print*, sum1,sum2, gamma,r0
         tmp = tmp - sum1*sum2*dble(gamma)
      end if

      RAD = tmp
      return
      end
c-----------------------------------------------------------
      subroutine di_el_pol
      use vmat_module, only: nodeid
      include 'par.f'
      common /di_el_core_polarization/ gamma, r0, pol(nmaxr)
      common /meshrr/ nr,gridr(nmaxr,3)
      common/smallr/ formcut,regcut,expcut,fast,match
      double precision r,rat,rr0,ww,vv

      pol(:) = 0.0
      do i=1,nr
         r = gridr(i,1)
         rat = r/dble(r0)
         rr0 = rat*rat*rat
c     Here    ww = ww6
         vv =  dexp(-rr0*rr0)
         if(real(vv).lt.expcut) exit
         ww = 1d0 - vv
         pol(i) = dsqrt(ww)/(r*r)
c     Here    ww = ww3
c         vv = dexp(-rr0)
c         ww = 1d0 - vv
c         pol(i) = ww/(r*r)
      end do
      iCR = i-1
      do i=iCR,nr
         r = gridr(i,1)
         pol(i) = 1d0/(r*r)
      end do

      call minmaxi(pol,nr,i1,i2)
      if (nodeid.eq.1)
     >write(4,'("Di-electron pol.pot.: gamma =",F9.5,", r0 =",F9.5)') 
     >     gamma, r0
      if (nodeid.eq.1)
     >write(4,'("i1 =",I5,", i2 =",I5,", pol(i2)=",E12.6)') 
     >     i1, i2, pol(i2)
c      write(*,'("Di-electron pol.pot.: gamma =",F9.5,", r0 =",F9.5)') 
c     >     gamma, r0
c      write(*,'("i1 =",I5,", i2 =",I5)') i1, i2
      return 
      end

c-----------------------------------------------------------
      double precision function fangint(lam,l1,l1p,l2,l2p,ll)

c      include 'par.f'
      parameter (lomax = 4)
      double precision f1, f2, fang
      common /angint/ fang(0:2*lomax,0:lomax,0:lomax,0:lomax,
     >   0:lomax,0:lomax)
      if(lam.gt.2*lomax.or.l1.gt.lomax.or.l1p.gt.lomax.or.l2.gt.lomax
     >   .or.l2p.gt.lomax.or.ll.gt.lomax) then
         f2 = 1.11D0
         icheck = 1
      else
         icheck = 0
         f2=fang(lam,l1,l1p,l2,l2p,ll)
      end if
      fangint = f2
      if(f2.eq.1.11D0) then
         rlam=lam
         rl1=l1
         rl1p=l1p
         rl2=l2
         rl2p=l2p
         rll=ll
         rf1 =real((-1)**(l1p+l2+ll))*
     >      sqrt((2.*l1p+1.)*(2.*l2p+1.))*
     >      CGC0(rl1p,rlam,rl1)*CGC0(rl2p,rlam,rl2)*
     >      COF6J(rl1,rl2,rll,rl2p,rl1p,rlam)
         f1 = DBLE(rf1)
         fangint = f1
         if(icheck.eq.0) then
            fang(lam,l1,l1p,l2,l2p,ll) = f1
            fang(lam,l2,l2p,l1,l1p,ll) = f1
         end if
      end if
      return
      end
c-----------------------------------------------------------
c     this subroutine is called from structure.f

      subroutine setfang

c      include 'par.f'
      parameter (lomax = 4)
      double precision fang
      common /angint/ fang(0:2*lomax,0:lomax,0:lomax,0:lomax,
     >   0:lomax,0:lomax)
      do lam=0,2*lomax
         do l1=0,lomax
            do l1p=0,lomax
               do l2=0,lomax
                  do l2p=0,lomax
                     do ll=0,lomax
                        fang(lam,l1,l1p,l2,l2p,ll) = 1.11D0
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
      end
c-----------------------------------------------------------------
c*****************   b - matrix        ***************************
c-----------------------------------------------------------------
      double precision function bmat(n1,n2,n1p,n2p)

      include 'par.f'
      double precision  ortint, ort1, ort2      
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      bmat = 0.0D0
      if(lo(n2).eq.lo(n2p)) then
         ort2 = ortint(n2,n2p)
         if(lo(n1).eq.lo(n1p)) then
            ort1 = ortint(n1,n1p)
            bmat = ort1 * ort2
         end if
      end if
      return
      end
c-----------------------------------------------------------------
c*********** Calculation of one-electron matrix elements *********
c-----------------------------------------------------------------
      double precision function oneel(n1,n2,n1p,n2p)

      include 'par.f'
      common /hame1/ e1r(nspmax,nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      double precision  e1r, ortint, ort
      oneel = 0.0D0
      if(lo(n2).eq.lo(n2p)) then
         ort = ortint(n2,n2p)
         if(ort.ne.0.0D0) then
            if(lo(n1).eq.lo(n1p)) then
               oneel = e1r(n1,n1p) * ort
               oneel = 2.D0 * oneel
            end if
         end if
      end if
c      write(4,'("oneel = ",F10.5)') real(oneel)
      return
      end


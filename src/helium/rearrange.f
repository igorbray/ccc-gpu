c     This routine makes projection operator from s.p.orbitals to be used
c     in the scattering part to resolve nonuniqueness problem.
c      
      subroutine make_po(nspm,nr,fl,maxf,minf,lo,Nmax,nspmW,C)
      use po_module
      include 'par.f'
      integer nspm, nr
      dimension  fl(nmaxr,nspmax),lo(nspmax)
      dimension maxf(nspmax),minf(nspmax)
      integer Nmax, nspmW
      double precision, dimension(Nmax+1,nspmW,nspmW) :: C
      integer n, nst, i, j
      integer usedornot(nspm)
      common /corearray/ nicm, ncore(nspmCI)

      usedornot(:) = 0

c     find how many orbitals have been used in CI calculations:
c     it is often not all orbitals declared in ccc.in are used in CI calculations 
      j = 0
      do n=1,nspm

         do nst=1,Nmax
            do i=1,nspm
               if(C(nst,n,i) .ne. 0.0) then
                  j = j + 1                
                  usedornot(n) = 1
                  goto 10
               endif
            enddo
         enddo         
 10   enddo
      nspm_po = j
      if(allocated(fl_po)) then
         deallocate(fl_po,lo_po,minf_po,maxf_po)
      endif
      allocate(fl_po(nr,nspm_po),lo_po(nspm_po),
     >   minf_po(nspm_po),maxf_po(nspm_po))
      fl_po(:,:) = 0.0

      j = 0

c     first do core orbitals
      do  nic=1,nicm
         n=ncore(nic)
         if(usedornot(n) .eq. 1) then
            j = j + 1
            lo_po(j) = lo(n)
            minf_po(j) = minf(n)
            maxf_po(j) = maxf(n)
            fl_po(1:nr,j) = fl(1:nr,n)
         endif
      enddo
      
      do  n=1,nspm
c     core orbitals are already done, exlude them
         iscore = 0
          do  nic=1,nicm
             nc = ncore(nic)
             if(nc .eq. n) then
                iscore = 1
                exit
             endif
          enddo
          if(iscore .eq. 1) cycle
c
         if(usedornot(n) .eq. 1) then
            j = j + 1
            lo_po(j) = lo(n)
            minf_po(j) = minf(n)
            maxf_po(j) = maxf(n)
            fl_po(1:nr,j) = fl(1:nr,n)
         endif
      enddo
c$$$      print*,'Made s.p. projection operator with nspm_po=', nspm_po      
      if(j .ne. nspm_po) then
         print*,'rearramge.f: j .ne. nspm_po:', j, nspm_po
         stop
      endif

      return
      end
c      

c-------- This part of the code is to find the number of   -------------
c--------   frozen core func. N_fr_core = ionic core + frozen core orbitals ----------
c     Assume that core orbitals are orthogonal.
      subroutine find_nsp_tot(Nmax,nspmW,nspm,C,lo,nsp_tot)
      use vmat_module, only: nodeid
      include 'par.f'
      double precision  C(Nmax+1,nspmW,nspmW)
      common /corearray/ nicm, ncore(nspmCI)
      integer ltmp(lomax+1), lo(nspmax)
      integer itmp_core(nspmW)
     
c     The common block  /corearray/ nicm, ncore(nspmCI) must be set and
c     available at this stage.

c     set up array which is =1 if orbital is a core orbital, otherwise it is 0
      itmp_core(1:nspm) = 0
      itmp_core(ncore(1:nicm)) = 1
      


      nsp_tot = nicm
      do N=1,Nmax
         do nic=1,nicm
            ic=ncore(nic)
c     Here we find how many f.c. orbitals are built on the core orbital nic.
c     nlctmax is the number and ltmp(nlct) is the ang.mom. of these orbitals.
            nlct = 0
            do j=1,lomax+1
               ltmp(j) = -1
            end do
            do nsp=1,nspm
               if(C(N,ic,nsp).ne.0d0) then
c     This makes 'f.c.o. nspN based on the core orbital nic' orthogonal to all 
c     core orbitals.
                  if(itmp_core(nsp).eq.0) then                  
                     icheck = 1 ! check if new frozen-core orbital (new l) need to be added 
                     do il=1,max(nlct,1)
                        if(ltmp(il).eq.lo(nsp)) icheck = 0
                     end do
                     if(icheck.eq.1) then
                        nlct = nlct + 1
                        ltmp(nlct) = lo(nsp)
c                        print*, N,nic,nsp,lo(nsp),nlct,' * ',
c     >                       (ltmp(i), i=1,nlct)
                     end if                  
                  end if
               end if
            end do
            nlctmax = nlct
            nsp_tot = nsp_tot  + nlctmax
         end do
      end do
      if(nsp_tot.gt.nspmax) then
         print*,'stop in rearrange: nsp_tot > nspmax',nodeid       
         print*,'nsp_tot =', nsp_tot, ',  nspmax =',  nspmax
         stop
      end if
      if (nodeid.eq.1)
     >write(4,'("Number of frozen core functions is  nsp_tot =",I6)') 
     >      nsp_tot
      if (nodeid.eq.1)
     >   write(*,'("Number of frozen core functions is  nsp_tot =",I6)') 
     >      nsp_tot

      return 
      end
c********************************************************************************
      subroutine rearrange(Nmax,nspmW,nspm,nsp_tot,
     >     lo,ko,fl,maxf,minf,C,los)
      use vmat_module, only: nodeid
      include 'par.f'
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common/hame1/ e1r(nspmax,nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      common /meshrr/ nr, gridr(nmaxr,3)
      double precision  C(Nmax+1,nspmW,nspmW)
     >     , r1elk
      double precision  sym, sume, sumort, tmp, ortint, e1r
      dimension maxf(nspmax),minf(nspmax), koN(nsp_tot),
     >     los(nspmax),loN(nsp_tot),ngiveN(nsp_tot),ngiveic(nsp_tot),
     >     iorbmax(nsp_tot),ltmp(lomax+1)
      dimension  fl(nmaxr,nspmax),lo(nspmax),ko(nspmax)
      common /corearray/ nicm, ncore(nspmCI)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer na, nam
      integer npoint(nsp_tot),npt(nspm)
c
      double precision  e1rp(nsp_tot,nsp_tot)
      double precision  ortp(nsp_tot,nsp_tot)
      double precision  A(nspmW,nspmW)
      real p(nsp_tot)
      integer  iorbfc(nsp_tot,nspmW)
      double precision  ttt(nspmW,nspmW)
      integer itmp_core(nspmW)

c     In this routine  nspm  becomes the number of ionic core + frozen core orbitals
c     and is for a large calculation is a large number --> nspmax should be also 
c     set to larger number (700-1000).
      
c
c     The common block  /corearray/ nicm, ncore(nspmCI) must be set and
c     available at this stage.
c     In the start of this subroutine the nspm is number 
c     of s.p.functions (Laguerre or generated in Igor's programm) and
c     must be less then nspmCI.
c     After f.c.orbitals are made in this routine the nspm is number
c     of f.c.orbital plus number of core functions and 
c     must be less then nspmax. (Usually nspmCI<<nspmax)

      if (nodeid.eq.1)
     >write(20,'("******************************************")') 
      if (nodeid.eq.1)
     >write(20,'("rearranged s.p. functions ")') 
c     

c     restore the CI coef. to be as for antisymmetric configurations
c     C(N,{n1n2})|{n1n2}> rather than to be C(N,n1,n2)|n1,n2> where
c     |n1,n2> is not antisymmetrized
      do n1=1,nspm
         do n2=1,nspm
            sym = dsqrt(2.0d0)
            if(n1.eq.n2) sym=1.0d0
            do N=1,Nmax
               C(N,n1,n2) = C(N,n1,n2) * sym
            end do
         end do
      end do
      
c------------     
      imaxr = 0
      do n=1,nspm
         i1 = minf(n)
         i2 = maxf(n)
         imaxr = max(imaxr,i2)
         if(i1 .ne. 1)  fl(1:i1-1,n) = 0.0
         if(i2 .ne. nr) fl(i2+1:nr,n) = 0.0
      enddo

      A(:,:) = 0d0

      ngiveN(:) = 0
      ngiveic(:) = 0
      ortp(:,:) = 0d0
      e1rp(:,:) = 0d0

c     set up array which is =1 if orbital is a core orbital, otherwise it is 0
      itmp_core(1:nspm) = 0
      itmp_core(ncore(1:nicm)) = 1

c     Make frozen core orbitals
c     these are nicm core orbitals
      do nic=1,nicm
         ic=ncore(nic)
         iorbmax(nic) = 1
         iorbfc(nic,1) = ic
         loN(nic) = lo(ic)
         los(nic) = 1
         koN(nic) = ko(ic) 
      end do
      nspN = nicm
c     Here the frozen core orbitals are made.
c     Assume that core orbitals are orthogonal.
c    f.c.o. is orthogonal to all core orbital by construction.

c     We repeat the code from routine  find_nsp_tot()  just to avoid to declare arrays of N
      do N=1,Nmax
c         print*,'N=',N

c      print*, 'checking C'
c      do n1=1,nspm
c         write(*,'(I5,100F10.5)') n1, (C(N,n1,n2),n2=1,nspm)
c      enddo


         do nic=1,nicm
            ic=ncore(nic)
c     Here we find how many f.c. orbitals are built on the core orbital nic.
c     nlctmax is the number and ltmp(nlct) is the ang.mom. of these orbitals.
            nlct = 0
            do j=1,lomax+1
               ltmp(j) = -1
            end do
            do nsp=1,nspm
               if(C(N,ic,nsp).ne.0d0) then
                  if(itmp_core(nsp).eq.0) then  ! proceed if nsp is not a core orbital
                     icheck = 1 ! check if new frozen-core orbital (new l) need to be added 
                     do il=1,max(nlct,1)
                        if(ltmp(il).eq.lo(nsp)) icheck = 0
                     end do
                     if(icheck.eq.1) then
                        nlct = nlct + 1
                        ltmp(nlct) = lo(nsp)
c                        print*,nic,ic,nsp,lo(nsp),nlct,
c     >                       (ltmp(i), i=1,nlct)
                     end if                  
                  end if
               end if
            end do
            nlctmax = nlct
c            print*, N,nic,ic, nlctmax
c     if nlctmax stays  =0 then the core ic=ncore(ic) was not used
c     for the state N.
            do nlct=1,nlctmax
               nspN = nspN+1
               ngiveN(nspN) = N
               ngiveic(nspN) = ic
               los(nspN)= ls(N)
               loN(nspN) = ltmp(nlct)
               ltest = ltmp(nlct)
               ict = 0
               do nsp=1,nspm
                  if(ltest.eq.lo(nsp)) then
                     if(C(N,ic,nsp).ne.0d0) then                        
                        icheck1=0
                        do nicx=1,nicm ! nic
                           if(ncore(nicx).eq.nsp) icheck1 = 1
                        end do
                        if(icheck1.eq.0) then  
                           ict = ict + 1
                           iorbfc(nspN,ict) = nsp
                        end if
                     end if
                  end if
               end do
               iorbmax(nspN) = ict
            end do              ! end nlct loop
         end do                 ! end nic loop
      end do                    ! end N loop
      nspNm = nspN
      if(nspNm.ne.nsp_tot) then
         print*,'nspNm.ne.nsp_tot stop in rearrange.f '
         stop
      end if
      if(nspNm.gt.nspmax) then
         write(*,'("stop in rearrange(): nsp_tot=nspNm > nspmax, ",
     >      "increase  nspmax  in par.f")')
         write(*,'("nspNm =",I5,", nspmax =",I5)') nspNm, nspmax
         stop
      end if
c     Arrays   iorbfc(nspN,ict)  and   iorbmax(nspN)  give representation of 
c     f.c. functions in terms of original s.p.functions.  
      do i=1,imaxr
         do nspN=1,nicm
            ic=ncore(nspN)
            p(nspN) = fl(i,ic)
         end do
         do nspN=nicm+1,nspNm
            ic = ngiveic(nspN)
            N = ngiveN(nspN)
            sume = 0.0D0
            do ict=1,iorbmax(nspN)
               nsp = iorbfc(nspN,ict)
               tmp = 0d0
               do nicp=1,nicm  ! this makes f.c.o. nspN orthogonal to all core orbitals
                  icp = ncore(nicp)
                  if(lo(icp).eq.lo(nsp)) then
                     tmp = tmp + ortint(nsp,icp)*dble(fl(i,icp)) 
                  end if
               enddo
               if(i.ge.minf(nsp).and.i.le.maxf(nsp))
     >              sume = sume + C(N,ic,nsp)*(dble(fl(i,nsp)) - tmp)
            end do
            p(nspN) = sume
         end do
         do nspN=1,nspNm
            fl(i,nspN) = p(nspN)
         end do
      end do                    ! end i loop
      do nspN=1,nspNm
         call minmaxi(fl(1,nspN),nr,i1,i2)
         maxf(nspN) = i2
         minf(nspN) = i1
      end do
c     
c     form new array e1r(nsp1,nsp2) and ortint(nsp1,nsp2)
c     Three cases must be considered separately: (core-core),(core-fr.core),
c     (fr.core-fr.core) because core and fr.core orbitals are given by
c     different formulas.
c     Assume that core orbitals are orthogonal.
      do m1=1,nspm
         do m2=1,nspm
            ttt(m1,m2)=SUM(e1r(m1,ncore(1:nicm))*
     >           ortint(m2,ncore(1:nicm)))
         enddo
      enddo

c     core-core
      do nic1=1,nicm
         ic1=ncore(nic1)           
         do nic2=nic1,nicm
            ic2=ncore(nic2)           
            e1rp(nic1,nic2) = e1r(ic1,ic2)
            e1rp(nic2,nic1) = e1rp(nic1,nic2)
            ortp(nic1,nic2) = ortint(ic1,ic2)
            ortp(nic2,nic1) = ortp(nic1,nic2)
         end do
      end do
c     core-fr.core
      do nic=1,nicm
         ic=ncore(nic)
         do nspN=nicm+1,nspNm
            N = ngiveN(nspN)
            icN = ngiveic(nspN)
            sume = 0d0
            if(loN(nic).eq.loN(nspN)) then
               do ict=1,iorbmax(nspN)
                  nsp1 = iorbfc(nspN,ict) 
                  sume = sume + C(N,icN,nsp1)*(e1r(ic,nsp1)
     >                 - ttt(ic,nsp1))
               end do
            end if
            e1rp(nspN,nic) = sume
            e1rp(nic,nspN) = sume
            ortp(nspN,nic) = 0d0 ! orthogonal by construction of f.c.o.
            ortp(nic,nspN) = 0d0
         end do
      end do
c     fr.core-fr.core
      do nspN1=nicm+1,nspNm
         N1 = ngiveN(nspN1)
         ic1 = ngiveic(nspN1)
         do nspN2=nspN1,nspNm
            N2 = ngiveN(nspN2)
            ic2 = ngiveic(nspN2)
            sume = 0d0
            sumort = 0d0
            if(loN(nspN1).eq.loN(nspN2))
     >         then
               do ict1=1,iorbmax(nspN1)
                  nsp1 = iorbfc(nspN1,ict1)
                  do ict2=1,iorbmax(nspN2)
                     nsp2 = iorbfc(nspN2,ict2)

                     sume = sume + C(N1,ic1,nsp1)*C(N2,ic2,nsp2)*
     >                    (e1r(nsp1,nsp2)-
     >                    ttt(nsp2,nsp1)-ttt(nsp1,nsp2)+
     >                    SUM(ortint(nsp1,ncore(1:nicm))*
     >                    ttt(ncore(1:nicm),nsp2)))
                     
                     sumort = sumort + C(N1,ic1,nsp1)*C(N2,ic2,nsp2)*
     >                    (ortint(nsp1,nsp2) -
     >                    SUM(ortint(nsp1,ncore(1:nicm))*
     >                    ortint(nsp2,ncore(1:nicm))))
                  end do
               end do
            end if
            e1rp(nspN1,nspN2) = sume
            e1rp(nspN2,nspN1) = sume
            ortp(nspN1,nspN2) = sumort
            ortp(nspN2,nspN1) = sumort
         end do
      end do


c      print*, 'checking ortp'
c      do n=1,nspNm
c         write(*,'(I5,100F10.5)') n, (ortp(n,m),m=1,nspNm)
c      enddo

c-----------------------------------------------------------------     
c     form new coefficient A(:,:)=C(N,:,:) 
      
      C(:,:,:) =  C(:,:,:) / sqrt(2d0) 
      do n=1,nspm
         C(:,n,n) = C(:,n,n) * sqrt(2d0)
      enddo

      do N=1,Nmax

         A(:,:) = 0d0
         A(1:nicm,1:nicm) = C(N,ncore(1:nicm),ncore(1:nicm))

         do nic1=1,nicm
            ic1=ncore(nic1)           
            do nic2=1,nicm
               ic2=ncore(nic2)           
               
               do nsp=1,nspm
                  if(itmp_core(nsp).eq.0) then
                     A(nic1,nic2) =  A(nic1,nic2) + 
     >                    (C(N,ic1,nsp)*ortint(nsp,ic2) + 
     >                    C(N,nsp,ic2)*ortint(nsp,ic1))
                  end if
               end do
            end do
         enddo

         C(N,:,:)=0.0D0
         C(N,1:nicm,1:nicm) = A(1:nicm,1:nicm)  

      end do  ! end N loop
c-----------------------------------------------------------------     
      do N=1,Nmax
         A(1:nicm,1:nicm) =  C(N,1:nicm,1:nicm)
         C(N,1:nicm,1:nicm) = 0d0
         
c     arrays na(i,N) and nam(N) hold pointers to s.p.orbitals for state N
c     To be used to access CI coef. : C(N,na(i1,N),na(i2,N)).
         na(:,N) = 0
         nam(N) = 0
         
c     make new C() array
         sym = dsqrt(2.0D0)
         
c     print*, 'N and A'
c     print*, N,(real(A(i)), i=1,nicm)
         
         npoint(:) = 0
         i = 0
         do n1=1,nicm
            do n2=1,nicm  ! check if orbital n1 is used for state N
               if(A(n1,n2).ne.0d0) then
                  i = i + 1
                  npoint(n1)=i
                  npt(i) = n1
                  exit
               endif
            enddo
         enddo
         na(1:i,N) = npt(1:i)
         C(N,1:i,1:i) = A(npt(1:i),npt(1:i))

c         print*, 'N=',N, ',  i=',i
c         do k=1,i
c            write(*,'(100F10.5)') (C(N,k,j),j=1,i)
c         enddo
         
         do nic=1,nicm
            do nsp=nicm+1,nspNm
               if(ngiveN(nsp).eq.N.and.ngiveic(nsp).eq.ncore(nic)) then
                  
c                  print*, 'nic=',nic,', nsp=',nsp,', i=',i, 
c     >                 (na(ij,N), ij=1,i)

c     check if core orbital nic has been included in the list, if not include it.
                  itest = 0
                  do ij=1,i
                     if(na(ij,N).eq.nic) itest = 1
                  end do
                  if(itest.eq.0) then
                     i = i + 1
                     na(i,N) = nic
                     npoint(nic) = i
                  end if
c
c     check if orbital nsp has been included, if not include it.
                  itest = 0
                  do ij=1,i
                     if(na(ij,N).eq.nsp) itest = 1
                  end do
                  if(itest.eq.0) then
                     i = i + 1
                     na(i,N) = nsp
                     npoint(nsp) = i
                  end if
c                  print*, '!!! i=',i, (na(ij,N), ij=1,i)
                   
c     
                  if(npoint(nsp) .gt. nspmW) then
                     print*,'rearange.f: npoint(nsp) .gt. nspmW : ',
     >                    npoint(nsp),nspmW 
                  endif
c     
c                  print*,'nic,nsp:',nic,nsp,npoint(nic),npoint(nsp)
                 
                  C(N,npoint(nic),npoint(nsp))= 1.0D0/sym
                  C(N,npoint(nsp),npoint(nic))= 1.0D0/sym * 
     >                 dble((-1)**(ls(N)+ll(N)+loN(nic)+loN(nsp)))
               end if
            end do
         end do
         nam(N) = i

c         print*, 'N=',N, ',  nam(N)=',nam(N)
c         do i=1,nam(N)
c            write(*,'(100F10.5)') (C(N,i,j),j=1,nam(N))
c         enddo
         
      end do   ! end N loop
c
c     Important: Redefinitiion of number of s.p.o.
      nspm = nspNm
c
c-----------------------------------------------------------------           
c     Choose correct sign of the frozen-core orbitals - positive at start.
c     The sign of the second (i=minf(nsp)+1) after the starting minf(nsp)
c     r-grid point is used to determine the sign     
      do nsp=1,nspm
         p(nsp) = 1.0           ! this is array which record the sign of the f.c.o. nsp
         if(fl(minf(nsp)+1,nsp) .lt. 0.0) then
            fl(minf(nsp):maxf(nsp),nsp)=-fl(minf(nsp):maxf(nsp),nsp)
c     print*,'In rearrange:', nsp, lo(nsp)
            p(nsp) = -1.0
         end if
      end do
c     restore correct sign in e1r ortint and C arrays.
      do n1=1,nspm
         do n2=1,nspm
            e1r(n1,n2)=e1rp(n1,n2) * p(n1) * p(n2)
            ortint(n1,n2)=ortp(n1,n2) * p(n1) * p(n2)           
         end do
      end do
      do N=1,Nmax
         do jn1=1,nam(N)
            n1 = na(jn1,N)
            do jn2=1,nam(N)
               n2 = na(jn2,N)
               C(N,jn1,jn2) = C(N,jn1,jn2) * p(n1) * p(n2)
            end do
         end do
      end do
c$$$
c$$$ 137  print*, 'Enter state number, negative number to terminate'
c$$$      read(*,'(I5)') N_read
c$$$c      open(731,file='tmp')
c$$$c      N_read = 2
c$$$      if(N_read .gt. 0) then
c$$$         N = N_read
c$$$         do i=1,nam(N)
c$$$            write(*,'(2I5,100(E12.4))') i,na(i,N),
c$$$     >         (real(C(N,i,j)), j=1,nam(N))
c$$$         end do
c$$$         go to 137
c$$$      end if
         
c-----------------------------------------------------------------           
      do n1=1,-nspm
         if(n1 .eq. 1) then
            if (nodeid.eq.1)
     >write(20,'("rearranged s.p. state overlaps")') 
            if (nodeid.eq.1)
     >write(20,'(4X,"l1   l2   n1   n2   s1   s2  <n1|n2>")')
         endif
         do n2=n1,nspm
            if(loN(n1).eq.loN(n2)) then
               if(n1.ne.n2.and.dabs(ortint(n1,n2)).le.1.0D-08) then
                  ortint(n1,n2) = 0.0D0
                  ortint(n2,n1) = 0.0D0
c                 write(20,'(" <n1|n2> <= 1.0D-08 ")')
               end if
               if(dabs(ortint(n1,n2) - 1.0D0).le.1.0D-08) then
                  ortint(n1,n2) = 1.0D0
                  ortint(n2,n1) = 1.0D0
c                  write(20,'(" <n1|n2> - 1.0D0 <= 1.0D-08 ")')
               end if
               if (nodeid.eq.1)
     >write(20,'(6I5,1P,2E20.10E2)') loN(n1),loN(n2),n1,n2,
     >              los(n1),los(n2),real(ortint(n1,n2))
     >             , r1elk(0,fl(1,n1),fl(1,n2),
     >               minf(n1),minf(n2),maxf(n1),maxf(n2))

            end if
         end do
      end do
c     
      lo(1:nspm) = loN(1:nspm)         
      ko(1:nicm) = koN(1:nicm)         

c     After this routine the ion-core orbital are first nicm orbitals.
      do nic=1,nicm
         ncore(nic) = nic
      end do
      if (nodeid.eq.1)
     >write(20,'("rearranged s.p. functions:")')
      if (nodeid.eq.1)
     >write(20,'("if N=ngiveN(n) = 0, ic=ngiveic(n)=0  then ",
     >   "this is a core function")')
      if (nodeid.eq.1)
     >write(20,'(4X,"n    l   los   N   ic  minf maxf")') 
      do n=1,nspm
c     if ngiveN(n) = 0 and ngiveic(n)=0 then this is a core orbital
         if (nodeid.eq.1)
     >write(20,'(6I5,I6)') n,lo(n),los(n),ngiveN(n),ngiveic(n),
     >      minf(n),maxf(n)
      end do
c---------------------------------------------------------------------
c$$$      do n=1,nspm
c$$$         do m=n,nspm
c$$$            if(lo(n).eq.lo(m)) then
c$$$            tmp1 = r1elk(0,fl(1,n),fl(1,m),
c$$$     >         minf(n),minf(m),maxf(n),maxf(m))
c$$$            print*, 'n, m:', n,m,tmp1
c$$$            end if
c$$$         end do
c$$$      end do
c     
      return
      end
c------------------------------------------------------------------------
      subroutine subset(fl,minf,maxf,los,i_orb_subset)
      include 'par.f'
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common/hame1/ e1r(nspmax,nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      common /meshrr/ nr, gridr(nmaxr,3)
      double precision  ortint, e1r
      dimension maxf(nspmax),minf(nspmax),maxp(KNM),minp(KNM)
      dimension  p(KNM),fl(nmaxr,nspmax),los(nspmax)
      dimension ikko(0:lomax), nmms(KNM)

      double precision  e1rp(KNM,KNM)
      double precision  ortintp(KNM,KNM)


      do l=0,lomax
         ikko(l) = 0
      end do
      
      if(i_orb_subset .eq. 1) then
         i_orb_spin = 1
      else
         i_orb_spin = 0
      endif

      mm = 0
      do n=1,nspm
c         if(los(n) .eq. i_orb_spin .or. n.eq.1) then
         if(1 .eq. 1) then  ! use all orbitals generated
           l = lo(n)
            mm = mm + 1
            nmms(mm) = n
            ikko(l) = ikko(l) + 1
            ko(mm) = ikko(l)
            lo(mm) = lo(n)
            maxp(mm) = maxf(n)
            minp(mm) = minf(n)
         end if
      end do
      mmax = mm
      
      do i=1,nr
         do mm=1,mmax
            p(mm) = fl(i,nmms(mm))
         end do
         do mm=1,mmax
            fl(i,mm) = p(mm)
         end do         
      end do

c      
      do m2=1,mmax
         do m1=1,mmax
            e1rp(m1,m2) = e1r(nmms(m1),nmms(m2))
            ortintp(m1,m2) = ortint(nmms(m1),nmms(m2))
         end do
      end do
c
      do n2=1,nspm
         do n1=1,nspm
            ortint(n1,n2)=0.0D0
            e1r(n1,n2)=0.0D0
         end do
      end do
      nspm = mmax

      do n=1,nspm
         maxf(n) = maxp(n)
         minf(n) = minp(n)
         nset(n) = 1
      end do
c     
      do n2=1,nspm
         do n1=1,nspm
            e1r(n1,n2)=e1rp(n1,n2)
            ortint(n1,n2)=ortintp(n1,n2)
         end do
      end do
c     
      call gsort(nspm,fl,maxf,minf,ortint,e1r,lo)

      return
      end

c----------------------------------------------------------------------
      subroutine stNcore(Nmax,nspmW,N,C,la,is,lpar,En,enionry,
     >     inum,Nofgfcst_add,jpartN,partN,en_max,partstN,
     >     j_summed,Nadd,Num_ioncore)
      use vmat_module, only: nodeid
      include 'par.f'
      common /corepartarray/ corepart(KNM,nicmax), ecore(nicmax)
      common /corearray/ nicm, ncore(nspmCI)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      real partstN(KNM)
      double precision  C(Nmax+1,nspmW,nspmW),sym,symp,corepart,
     >     tmp1,En,ortint,ecore
      character inornot(2)*1
      data inornot /"+","-"/
      integer Num_ioncore(nspmCI)
      character chan(knm)*3
      common /charchan/ chan

      eV = 27.2116
c
c     core orbitals are determined in the config() subroutine (hel-3.f).
c     corepart(N,nic) is the probability to find state N in the manifold
c     based on the core orbital ic=ncore(nic).
c     Will be used for ionization with excitation where the excitation 
c     cross-section should be multiplied by corepart(N,nic) to get excitation 
c     corresponding to the manifold based on the core orbital ic=ncore(nic)
c     Note it is made before the frozen-core orbitals are made 
c     !!! It is assumed that s.p. functions that forms ionic core 
c     are orthogonal to each other and that ionic core functions are
c     orthogonal to the rest of the s.p. function, but the latter functions
c     may be nonorthoganal to each other. This is used when optimised orbitals
c     are constructed (they are not orthoganal to between itself and to 
c     the non-ionic core orbitals).
c
      if(nicm.gt.nicmax) then
         print*,'nicm.gt.nicmax, increase nicmax in file rearrange.f',
     >        'nicm =', nicm
         stop
      end if
c     
      tmp = 0.0
      do nnic=1,nicm
         ic = ncore(nnic)
         tmp1 = 0d0
         do nsp=1,nspm
            icheck=0
            do nicx=1,nnic-1
               if(ncore(nicx).eq.nsp) then
                  icheck = 1
               end if
            end do
            if(icheck.eq.0) then
               sym = dsqrt(2.0d0)
               if(ic.eq.nsp) sym=1.0d0
               do nspp=1,nspm
                  icheckp=0
                  do nicx=1,nnic-1
                     if(ncore(nicx).eq.nspp) then
                        icheckp = 1
                     end if
                  end do
                  if(icheckp.eq.0) then
                     symp = dsqrt(2d0)
                     if(ic.eq.nspp) symp=1.0d0
                     tmp1 = tmp1 + sym*symp*C(N,ic,nsp)*C(N,ic,nspp)*
     >                    ortint(nsp,nspp)
                  end if
               end do
            end if
         end do
         corepart(N,nnic) = tmp1
         tmp = tmp + corepart(N,nnic)
      end do
      if(abs(tmp).le.0.999.or.abs(tmp).ge.1.001) then
         if (nodeid.eq.1)
     >write(10,'("Warning: "," N =",i3,
     >        " in subroutine stNcore(..):",
     >         " tmp.ne.1.0; tmp =",F10.5)') N,tmp
      end if

      st_energy_eV = (En-enionry/2.0)*eV ! above frozen ionic core.

c     exclude or include states by ionic core orbitals
c     (1s, 1s+2s, 1s+2s+2p,...)
c     j_summed  - is how many ionic core orbital components should be summed,
c                 it should be less then nicm.
c     Nadd > 1 only if j_summed > 1.
      tmp = 0.0
c     go through all negative states, or 
c     when going first time (Nadd=0) go through all states
      if(st_energy_eV.lt.0.0 .or. Nadd.eq.0) then
         do ii=1,j_summed 
            nnic = Num_ioncore(ii)
            tmp = tmp + corepart(N,nnic)
         end do
         partstN(N) = tmp
      else
c     when going second (and all other) times first check if the state 
c     has been already included (negative energy states are always included, 
c     see code above), in this case the summation up to j_summed-1 
c     must satisfy the condition   partstN(N) > tmp
         do ii=1,j_summed - 1
            nnic = Num_ioncore(ii)
            tmp = tmp + corepart(N,nnic)
         end do
         partstN(N) = tmp
         if(partstN(N).le.partN) then
c     this state was not included before.
c     check if the summation up to j_summed will satisfy the condition 
c     partstN(N) > tmp. Only  Nadd  states can be included this way.
            if(Nofgfcst_add.lt.Nadd) then
               ii=j_summed
               nnic = Num_ioncore(ii)
               tmp = tmp + corepart(N,nnic)
               partstN(N) = tmp
               if(st_energy_eV.lt.en_max .and. partstN(N).gt.partN) 
     >              Nofgfcst_add = Nofgfcst_add + 1
c               print*, inum,N,la,partstN(N), Nadd, Nofgfcst_add
            end if
         end if
      end if
      
      if(jpartN.eq.0) then
         if (nodeid.eq.1)
     >write(15,'(1A,1X,I3,1X,"- -",1X,"- -",1X,2I2,I3,
     >      F8.3,2X,F5.3,1X,100(1X,F5.3))') 
     >      inornot(2),inum, la, is, lpar,st_energy_eV, 
     >      partstN(N), (real(corepart(N,i)), i=1,nicm)
      else        
         ii = 2
         jpartN = 0
         if (st_energy_eV.lt.0.0 .or. partstN(N).gt.partN) then
            ii = 1
            if(st_energy_eV.lt.en_max) jpartN = 1 ! state may be included     
         end if
         if(N.gt.Nmax.or.ii.eq.2) then         
            if (nodeid.eq.1)
     >write(15,'(1A,1X,I3,1X,"---",1X,"---",1X,2I2,I3,
     >         F8.3,2X,F5.3,1X,100(1X,F5.3))') 
     >         inornot(ii),inum, la, is, lpar,st_energy_eV, 
     >         partstN(N), (real(corepart(N,i)), i=1,nicm)
         else
            if (nodeid.eq.1)
     >write(15,'(1A,1X,I3,1X,I3,1X,A3,1X,2I2,I3,
     >         F8.3,2X,F5.3,1X,100(1X,F5.3))') 
     >         inornot(ii),inum, N, chan(N), la, is, lpar,st_energy_eV, 
     >         partstN(N), (real(corepart(N,i)), i=1,nicm)            
         end if
      end if
      return
      end
c************************************************************************   
      subroutine printcoreparts(partN,Nmax,nspmW,C,E,Nofgfcst,partstN)
      use vmat_module, only: nodeid
      include 'par.f'
      common /corepartarray/ corepart(KNM,nicmax),ecore(nicmax)
      common /corearray/ nicm, ncore(nspmCI)
      real partstN(KNM)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)      
      double precision  C(Nmax+1,nspmW,nspmW),E(KNM), corepart, e1r,
     >   ecore
      character chan(knm)*3
      common /charchan/ chan
      common/hame1/ e1r(nspmax,nspmax)
      if(nicm.gt.nicmax) then
         print*,'nicm.gt.nicmax, increase nicmax in file rearrange.f',
     >        'nicm =', nicm
         stop
      end if
c     
      if (nodeid.eq.1)
     >write(10,'("start printcoreparts(...)")')

      if (nodeid.eq.1) then
         write(*,'("core orbitals:   n   l(n)   ko(n)")')
         do nic=1,nicm
            ic = ncore(nic)
            write(*,'("             ",3I5)') ic, lo(ic), ko(ic)
         end do
      endif
      
      if (nodeid.eq.1)
     >write(10,'("corepart(N,nic) is the probability ",
     >   " to find state N")')
      if (nodeid.eq.1)
     >write(10,'("in the manifold based on core orbital i=1,nicm")')
      if (nodeid.eq.1)
     >write(10,'("partstN(N) is the probability to find state N")')
      if (nodeid.eq.1)
     >write(10,'("in the manifold based on the core orbitals ")') 
      if (nodeid.eq.1)
     >write(10,
     >   '("   nic  ic=ncore(nic) lo(ic)    ko(ic)   e1r(ic,ic)")')
      
      do nic=1,nicm
         ic = ncore(nic)
         ecore(nic) = e1r(ic,ic)
         if (nodeid.eq.1)
     >write(10,'(4(I5,5X),E12.5)') nic,ic,lo(ic),ko(ic),e1r(ic,ic)
      end do
      if (nodeid.eq.1)
     >write(10,'("  N      l s par   E(N)a.u.     partstN(N)",12X,
     >   "corepart(N,nic)")')
      if (nodeid.eq.1)
     >write(10,'(33X," nic= ",100I13)') (i, i=1,nicm)
      if (nodeid.eq.1)
     >write(10,'(33X," l(i)=",100I13)') (lo(ncore(i)), i=1,nicm)

      ii = 0
      do N=1,Nmax
         if (nodeid.eq.1)
     >write(10,'(I3,1x,A3,1X,2I2,I3,E15.5,2X,E12.5,1X,
     >        100(1X,E12.5))') N, chan(N), ll(N), ls(N), lparity(N),
     >        real(E(N)), partstN(N), (real(corepart(N,i)), i=1,nicm)
         if(partstN(N).gt.partN) ii = ii + 1
      end do
      if (nodeid.eq.1)
     >write(10,'("There were ",I3," states with  ",
     >     "partstN(N) > partN =",E10.5)') ii, partN
      if(ii.lt.Nofgfcst) then
         if (nodeid.eq.1)
     >write(10,'("But ",I3," of states with partstN(N) > partN",
     >        " were not included in the calculation.")') Nofgfcst - ii
         if (nodeid.eq.1)
     >write(10,'("To include those states increase  nstate(l,is,ip)",
     >        "  in file F5.")')
      end if
      
      
      return
      end

c------------------------------------------------------------------------
      subroutine natorb(fl,minf,maxf,Nmax,nspmW,C)
      use vmat_module, only: nodeid
      include 'par.f'
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common/hame1/ e1r(nspmax,nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      common /meshrr/ nr, gridr(nmaxr,3)

      dimension  fl(nmaxr,nspmax)
      dimension maxf(nspmax),minf(nspmax)
      integer Nmax, nspmW
      double precision, dimension(Nmax+1,nspmW,nspmW) :: C

      double precision, dimension(nspm,nspm) :: CI

      double precision  ortint, e1r
      integer nk(0:lomax), nkm(0:lomax)
      real tmpfl(nmaxr)
      real flp(nr,nspm)
      dimension maxfp(nspm),minfp(nspm)
      dimension  p(nspm)
      dimension lop(nspm),kp(nspm), lflp(nspm),kflp(nspm)
      dimension Num(nspm), Numinv(nspm), Num1(nspm)
      double precision, dimension(nspm,nspm) :: z, tm
      double precision, dimension(nspm):: w,fv1,fv2

      double precision  tmpar(nspm,nspm)
      double precision  tmp
!

      if (nodeid.eq.1) then
         write(*,'("making natural orbitals and associated work")')
         write(4,'("making natural orbitals and associated work")')
      endif
      tmpfl(:) = 0.0

      CI(1:nspm,1:nspm) = C(1,1:nspm,1:nspm)

      ltmpmax = 0
      do n=1,nspm
         ltmpmax = max(ltmpmax,lo(n))
      enddo


      imaxr = 0
      do n=1,nspm
         i1 = minf(n)
         i2 = maxf(n)
         imaxr = max(imaxr,i2)
         if(i1 .ne. 1)  fl(1:i1-1,n) = 0.0
         if(i2 .ne. nr) fl(i2+1:nr,n) = 0.0
!         print*,'n,i1,i2: ', n, lo(n), ko(n), i1, i2
      enddo
      print*,'imaxr= ', imaxr

      flp(:,:) = 0.0


!
!      do n=1,nspm
!         write(*,'(50F10.5)') (CI(n,m), m=1,nspm)         
!      enddo

! Assume that s.p.orbitals come from Igors routines and  orthonomal (ortint = I) 
! but not required to diagonalise one-electron Hamiltonian (e1r .not. diagonal matrix) 
! CI is an array if CI coef for the ground state : singlet S-state :
! It is a symmetric array and can be diagonalized.


! CI array is not necessary is in block-diagonal form: first block for 
! s-orbital then block for p-orbital, ... (ordering of s.p.o can be 
! arbitrary: 1s,2p,2s,3s,...))
! Sort according to orbital angular momentum.
c     This array will be sorted by insertion sort algorithm.
      do n=1,nspm
         Num(n) = n
      end do
      do n=2,nspm
         itmp = Num(n)
         j = n
         do while(j.ge.2 .and. lo(itmp).lt.lo(Num(j-1)))
            Num(j) = Num(j-1)
            j = j - 1
         end do
         Num(j) = itmp
      end do
      

      do n=1,nspm
         Numinv(Num(n)) = n
      end do      
      Num1(1:nspm) = Num(1:nspm)

!      print*,'New CI order'
      CI(1:nspm,1:nspm) = CI(Num(1:nspm),Num(1:nspm)) 
!      do n=1,nspm
!         write(*,'(50F10.5)') (CI(n,m), m=1,nspm)         
!      enddo

!      
      matz=2
      call  rs(nspm,nspm,CI,w,matz,z,fv1,fv2,ierr)
      print*, 'natorb: ierr =', ierr

c     This array will be sorted by insertion sort algorithm.
c     sorting by nat.orb. weights.
      do n=1,nspm
         Num(n) = n
      end do
      do n=2,nspm
         itmp = Num(n)
         j = n
         do while(j.ge.2)
            if(abs(w(itmp)).gt.abs(w(Num(j-1)))) then
               Num(j) = Num(j-1)
               j = j - 1
            else
               exit
            endif
         end do

!         print*,'!! j,itmp:',j,itmp, Num(j)

         Num(j) = itmp
      end do


      w(1:nspm) = w(Num(1:nspm))
      z(1:nspm,1:nspm) = z(Numinv(1:nspm),Num(1:nspm))

      do n=1,nspm
         tmp = SUM(z(1:nspm,n)*z(1:nspm,n) * dble(lo(1:nspm)))
         lop(n) = nint(tmp)
!         print*,n,lop(n), tmp
      enddo


      do l=0,ltmpmax
         kk = 0
         do n=1,nspm
            if(l .eq. lop(n)) then
               kk = kk + 1
               kp(n) = kk
            endif
         enddo
      enddo

      if (nodeid.eq.1)
     >write(4,'("Nat.orb. weights:")')
      if (nodeid.eq.1)
     >write(*,'("Nat.orb. weights:")')
      if (nodeid.eq.1)
     >write(4,'(" Sorting: n,lop(n),kp(n),abs(w(n))")')
      if (nodeid.eq.1)
     >write(*,'(" Sorting: n,lop(n),kp(n),abs(w(n))")')
      do n=1,nspm
         if (nodeid.eq.1)
     >write(4,'(3I5,5E15.6)') n,lop(n),kp(n),w(n)*w(n)
         if (nodeid.eq.1)
     >write(*,'(3I5,5E15.6)') n,lop(n),kp(n),w(n)*w(n)
      enddo


!      write(4,'("Nat.orb.    3s              4s            3p      ",
!     >     "      4p              3d")')
!      write(*,'("Nat.orb.    3s              4s            3p      ",
!     >     "      4p              3d")')
      do n=1,nspm
         if(lop(n) .eq. 0 .and. kp(n) .eq. 1 ) then 
            i3s=n
!            print*,'i3s=',i3s
         endif
         if(lop(n) .eq. 0 .and. kp(n) .eq. 2 ) then
            i4s=n
!            print*,'i4s=',i4s
         endif
         if(lop(n) .eq. 1 .and. kp(n) .eq. 1 ) then
            i3p=n
!            print*,'i3p=',i3p
         endif
         if(lop(n) .eq. 1 .and. kp(n) .eq. 2 ) then
            i4p=n
!            print*,'i4p=',i4p
         endif
         if(lop(n) .eq. 2 .and. kp(n) .eq. 1 ) then
            i3d=n
!            print*,'i3d=',i3d
         endif
      enddo
      do m=1,nspm
         n = Num1(m)
c$$$         write(*,'(I5,15E15.6)') n,z(n,i3s),z(n,i4s),
c$$$     >        z(n,i3p),z(n,i4p),z(n,i3d)
c$$$         write(4,'(I5,15E15.6)') n,z(n,i3s),z(n,i4s),
c$$$     >        z(n,i3p),z(n,i4p),z(n,i3d)
      enddo


c     read array nk(l) from F5 file, defining how many nat. orb. will be used 
c     to substitute for original orbitals
      read(3,*) lnk, (nk(l), l=0,lnk)
      write(*,'("Defining how many nat. orb. will be used", 
     >     " to substitute for original orbitals")')
      if (nodeid.eq.1)
     >write(4,'("Defining how many nat. orb. will be used", 
     >     " to substitute for original orbitals")')
      if (nodeid.eq.1)
     >write(*,*) lnk,  (nk(l), l=0,lnk), " */ lnk, nk(0:lnk)"      
      if (nodeid.eq.1)
     >write(4,*) lnk,  (nk(l), l=0,lnk), " */ lnk, nk(0:lnk)"      


c     read array nkm(l) from F5 file, defining how many orbitals wil be used for each l
      read(3,*) lnkm, (nkm(l), l=0,lnkm)
      if (nodeid.eq.1)
     >write(*,'("Defining how many orbitals will be used ",
     >     "for each l")')
      if (nodeid.eq.1)
     >write(4,'("Defining how many orbitals will be used ",
     >     "for each l")')
      if (nodeid.eq.1)
     >write(*,*) lnkm,  (nkm(l), l=0,lnkm), " */ lnkm, nkm(0:lnkm)"      
      if (nodeid.eq.1)
     >write(4,*) lnkm,  (nkm(l), l=0,lnkm), " */ lnkm, nkm(0:lnkm)"      

      ioffset = 0   ! for nk(l) < 0 :adding some original orbitals - need to reduce the number of nat.orbital to be included - hope to avoid linear depend. and would not go out of bounds due to larger number of orbitals than nspm.
      if(lnk .ge. 0) then
         do l=0,lnk
            if(nk(l) .lt. 0) then
c               ioffset = ioffset + abs(nk(l))
               nkm(l) = nkm(l) + nk(l)
            endif
         enddo
      endif
      

c     substitution of original orbitals is done if lnk>0 and 
c     if for orbital n the value of ko(n) is smaller then nk(l)
      nn = 0
      do l=0,ltmpmax
         kk = 0
         if(nk(l) .lt. 0) then
            do ij=1,abs(nk(l))
c$$$               m=0              ! finding nat.orb. with the same ko
c$$$               do j=1,nspm
c$$$                  if(lop(j) .eq. l .and. kp(j) .eq. ij) then
c$$$                     m=j
c$$$                     exit
c$$$                  endif
c$$$               enddo
c$$$               if(m.eq.0) then
c$$$                  print*, 'Error: rearange.f: natorb(): m=0,',
c$$$     >                 ' nk(l) <0'
c$$$                  stop
c$$$               endif
c$$$               nn = nn + 1
c$$$               kk = kk + 1
c$$$               kflp(nn) = kk
c$$$               lflp(nn) = l
c$$$               tm(1:nspm,nn) = z(1:nspm,m)
c$$$               p(1:nspm) = z(1:nspm,m)
c$$$               do i=1,imaxr
c$$$                  tmpfl(i) = SUM(fl(i,1:nspm) * p(1:nspm))
c$$$               enddo
c$$$               flp(1:imaxr,nn) = tmpfl(1:imaxr)
c$$$               call minmaxi(tmpfl,nr,i1,i2)
c$$$               maxfp(nn) = i2
c$$$               minfp(nn) = i1                     
               n=0              ! finding original.orb. fl() with the ko=ij
               do j=1,nspm
                  if(lo(j) .eq. l .and. ko(j) .eq. ij) then
                     n=j
                     exit
                  endif
               enddo
               if(n.eq.0) then
                  print*, 'Error: rearange.f: natorb(): n=0,',
     >                 ' nk(l) <0'
                  stop
               endif
               nn = nn + 1
               kk = kk + 1
               kflp(nn) = kk
               lflp(nn) = l
               tm(1:nspm,nn) = 0.0
               tm(n,nn) = 1.0
               flp(1:imaxr,nn) = fl(1:imaxr,n)                  
               maxfp(nn) = maxf(n)
               minfp(nn) = minf(n)
               
               nk(l) = abs(nk(l))
            enddo                 
         endif
         do n=1,nspm 
            if(l .eq. lo(n)) then
               if(ko(n) .le. nk(l) .or. lnk .lt. 0) then
                  m=0   ! finding nat.orb. with the same ko
                  do j=1,nspm
                     if(lop(j) .eq. l .and. kp(j) .eq. ko(n)) then
                        m=j
                        exit
                     endif
                  enddo
                  if(m.eq.0) then
                     print*, 'Error: rearange.f: natorb(): m=0,',
     >                    ' nk(l)>0'
                     stop
                  endif
                  nn = nn + 1
                  kk = kk + 1
                  kflp(nn) = kk
                  lflp(nn) = lop(m)
                  tm(1:nspm,nn) = z(1:nspm,m)
!                  print*, 'nn,m,nspm,imaxr:',nn,m,nspm,imaxr
                  p(1:nspm) = z(1:nspm,m)
                  do i=1,imaxr
                     tmpfl(i) = SUM(fl(i,1:nspm) * p(1:nspm))
                  enddo
                  flp(1:imaxr,nn) = tmpfl(1:imaxr)
                  call minmaxi(tmpfl,nr,i1,i2)
                  maxfp(nn) = i2
                  minfp(nn) = i1
               else
                  nn = nn + 1
                  kk = kk + 1
                  kflp(nn) = kk
                  lflp(nn) = l
                  tm(1:nspm,nn) = 0.0
                  tm(n,nn) = 1.0
                  flp(1:imaxr,nn) = fl(1:imaxr,n)                  
                  maxfp(nn) = maxf(n)
                  minfp(nn) = minf(n)
               endif
            endif
            if(kk .ge. nkm(l)) exit
         enddo

      enddo

      fl(1:imaxr,1:nn) = flp(1:imaxr,1:nn)
      maxf(1:nn) = maxfp(1:nn)
      minf(1:nn) = minfp(1:nn)
      lo(1:nn) = lflp(1:nn)
      ko(1:nn) = kflp(1:nn)

      tmpar(1:nspm,1:nspm) = e1r(1:nspm,1:nspm)
      do m=1,nn
         do n=1,nn
            tmp = 0d0
            do i=1,nspm
               do j=1,nspm
                  tmp =  tmp + tm(j,n) * tm(i,m) * tmpar(i,j)
               enddo
            enddo
            e1r(m,n) = tmp
         enddo
      enddo

      tmpar(1:nspm,1:nspm) = ortint(1:nspm,1:nspm)
      do m=1,nn
         do n=1,nn
            tmp = 0d0
            do i=1,nspm
               do j=1,nspm
                  tmp =  tmp + tm(j,n) * tm(i,m) * tmpar(i,j)
               enddo
            enddo
            ortint(m,n) = tmp
         enddo
      enddo

      nspm = nn

      do n=1,nspm
!         write(*,'(50F10.5)') (ortint(n,m), m=1,nspm)         
      enddo


      call gsort(nspm,fl,maxf,minf,ortint,e1r,lo)

      if (nodeid.eq.1)
     >write(*,*)
      if (nodeid.eq.1)
     >write(4,*)

      return
      end
c$$$      if(lnk .ge. 0) then
c$$$         m = 0
c$$$         do n=1,nspm
c$$$            if(ko(n) .le. nk(lo(n))) then
c$$$               m = m + 1
c$$$               Num(m) = n
c$$$            endif
c$$$         enddo
c$$$         
c$$$         nspm = m
c$$$      
c$$$         lo(1:nspm) = lo(Num(1:nspm))
c$$$         ko(1:nspm) = ko(Num(1:nspm))
c$$$         e1r(1:nspm,1:nspm) = e1r(Num(1:nspm),Num(1:nspm))
c$$$         ortint(1:nspm,1:nspm) = ortint(Num(1:nspm),Num(1:nspm))
c$$$         fl(1:imaxr,1:nspm) = fl(1:imaxr,Num(1:nspm))
c$$$         minf(1:nspm) = minf(Num(1:nspm))
c$$$         maxf(1:nspm) = maxf(Num(1:nspm))
c$$$
c$$$      endif

c----------------------------------------------------------------------


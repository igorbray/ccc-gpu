      subroutine noble1el(l,nps)
      include 'par.f'
      common/meshrr/ nr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >     ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >     npbot(0:lamax),nptop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >     psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common /worksp/
     >     ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax) 

      double precision ortint, e1r
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      common/hame1/ e1r(nspmax,nspmax)

      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)

      parameter (K_kb = 4)
      real pj(k_kb),cj(K_kb),tj(K_kb)
      real orb_kb(nr)
      integer min_kb, max_kb
      character ch_kb*2

      real fc(maxr)
      integer istopc
      integer im(1)
      real ovlp(nps), ovlp2(nps)
      integer ip(nps)
      real*8 tmpovlp(nps,nps), tmpe1r(nps,nps)
      real tmpps2(nr,nps)
      integer mintmp(nps), maxtmp(nps)


      Ry = 13.6058
      
      print*
      write(*,'("noble1el(): nps=",I5,", nspm=",I5,", l=",I5)') 
     >     nps,nspm,l

      itest = 0                 ! do not replace 2p orbital by hf orbital if itest = 1      
!!$   For helium only we do not call hf(), define here reqired values
      if(nznuc .eq. 2 .and. l .eq. 0) then
         itest = 1   ! no need to replace 1s orbital on hf orbital

      endif

c!!!!!!!!!!!      
      ioption_KB = 0
      if(ioption_KB .eq. 1) then
c     Use Klaus B. orbital from J.Phys.B 30(1997)4609
         K_kb_1 = K_kb
         ch_kb = '2p'
c     2p orbital
         if(ch_kb .eq. '2p') then
            pj(:) = 2
            cj(1) = 0.56762
            cj(2) = 0.28107
            cj(3) = 0.23310
            cj(4) = 0.01369
            tj(1) = 2.53704
            tj(2) = 4.75246
            tj(3) = 1.49352
            tj(4) = 9.59417
            K_kb_1 = 4
         elseif(ch_kb .eq. '3p') then
c     3p orbital
            pj(1) = 2
            pj(2) = 3
            cj(1) = 0.11864
            cj(2) = -0.99850
            tj(1) = 3.21475
            tj(2) = 0.54948
            K_kb_1 = 2
         elseif(ch_kb .eq. '3d') then
c     3d orbital
            pj(1) = 3
            cj(1) = 1.0
            tj(1) = 0.33636
            K_kb_1 = 1
         endif
         print*, '! Use Klaus.B  ', ch_kb,' orbital'
         
         
         orb_kb(:) = 0.0
         do i=1,nr
            r = rmesh(i,1)
            tmp = 0.0
            do k=1,K_kb_1
               tmp = tmp + cj(k)*(r**(pj(k)))*exp(-tj(k)*r)
            enddo
            orb_kb(i) = tmp
         enddo
         call minmaxi1(orb_kb,nr,i1,i2)
         print*,'!!!! i1,i2:', i1,i2
         max_kb = i2
         min_kb = i1      
         orbovlp = SUM(orb_kb(i1:i2)*orb_kb(i1:i2)*rmesh(i1:i2,3))
         print*,'test ovlp KB:', orbovlp
         
         orb_kb(i1:i2) =   orb_kb(i1:i2)/sqrt(orbovlp)
         call minmaxi1(orb_kb,nr,i1,i2)
         print*,'!!!! i1,i2:', i1,i2
         max_kb = i2
         min_kb = i1      
         orbovlp = SUM(orb_kb(i1:i2)*orb_kb(i1:i2)*rmesh(i1:i2,3))
         print*,'test ovlp KB:', orbovlp         
         
c     minf_ng = i1
c     maxf_ng = i2
c     f_ng(i1:i2) = orb_kb(i1:i2)
         
         open(168,file='tmp_ng-orb')
         
         i_max = min(max_kb, maxf_ng)
         do i=1,i_max,3
            write(168,'(F15.7,2E15.6)')  rmesh(i,1), orb_kb(i), f_ng(i)
         enddo
         
      endif
c!!!!!!!!!!
         

      if(l .ne. l_ng .or. itest .eq. 1) then
         do n=nspm+1,nspm+nps
            m = n - nspm
            lo(n) = l
            ko(n) = m
            i1 = minps2(m)
            i2 = maxps2(m)
            minf(n) = i1
            maxf(n) = i2
            fl(i1:i2,n) = ps2(i1:i2,m)
            e1r(n,n) = psen2(m)/2.0 ! in a.u.
            ortint(n,n) = 1.0
         enddo
         
         nspm = nspm + nps
         return

      endif
      if(itest .eq. 1) return

c     Note: get here only for l=l_ng, the task is to replace 'wrong'
c     2p orbital of Ne^{5+} obtained from one-electron diagonalization 
c     with 'correct' hf orbital of Ne^{+}. Need to orthogonlise 
c     the s.p.  basis.

c     Note l = l_ng below

      istopc = maxf_ng ! istoppsinb(n_ng,l)      
      fc(1:istopc) =  f_ng(1:istopc)  

      tmp = 0.0
      do n=1,nps
         i1 = minps2(n)
         i2 = min(istopc,maxps2(n))
         ovlp(n) = SUM(fc(i1:i2)*ps2(i1:i2,n)*rmesh(i1:i2,3))
         tmp = tmp + ovlp(n)*ovlp(n)
!         print*,'! n,psen2(n)',n,psen2(n)*Ry, ovlp(n), i1, i2
      enddo
      print*,'Laguere basis covers HF orbital:', tmp
      if(abs(tmp-1.0000000) .gt. 0.0001) then
         print*,'Really, you need to increase the size of Laguere ',
     >        'basis for l=',l
!         stop
      endif
!     make projection
      do n=1,nps
         i1 = minps2(n)
         i2 = min(istopc,maxps2(n))
         fc(i1:i2) = fc(i1:i2) + ps2(i1:i2,n)* ovlp(n)
      enddo
!     ensure that last few points are zero
      fc(istopc) = 0d0
      fc(istopc-1) = 0d0
      i1 = 1
      i2 = istopc      
!     normalize to unity      
      tmp = SUM(fc(i1:i2)*fc(i1:i2)*rmesh(i1:i2,3))
      tmp = sqrt(tmp)
      fc(i1:i2) = fc(i1:i2)/tmp
!     repeat overlap calculations
      tmp = 0.0
      do n=1,nps
         i1 = minps2(n)
         i2 = min(istopc,maxps2(n))
         ovlp(n) = SUM(fc(i1:i2)*ps2(i1:i2,n)*rmesh(i1:i2,3))
         tmp = tmp + ovlp(n)*ovlp(n)
!         print*,'! n,psen2(n)',n,psen2(n)*Ry, ovlp(n), i1, i2
      enddo
      print*,'Laguere basis covers HF orbital:', tmp
      if(abs(tmp-1.0000000) .gt. 0.0001) then
         print*,'Need to increase the size of Laguere basis',
     >        ' for l=',l
         stop
      endif





c     find max ovelap index
      im = MAXLOC(abs(ovlp))

c     place core orbital instead of this max-ovelap orbital and 
c     form new basis starting counting from this orbital. 
c     make ip(i) array that points to old basis.              
      ip(1) = im(1)
      i = 1
      do n=1,nps         
         if(n .eq. im(1)) cycle
         i = i + 1
         ip(i) = n
      enddo

c     deal with the hf core orbital: project it on s.p. basis
      i1 = 1
      i2 =  istopc
      do i=i1,i2
         fc(i) = SUM(ovlp(1:nps)*dble(ps2(i,1:nps)))
      enddo
      tmp = SUM(fc(i1:i2)*fc(i1:i2)*rmesh(i1:i2,3))
      fc(i1:i2) = fc(i1:i2)/sqrt(tmp)

c     correct the ovelaps (note: the change was only for the first orbital)
      do n=1,nps
         i1 = minps2(n)
         i2 = min(istopc,maxps2(n))
         ovlp(n) = SUM(fc(i1:i2)*ps2(i1:i2,n)*rmesh(i1:i2,3))
c         print*,'n,psen2(n)',n,psen2(n)/2.0,psen2(n)*Ry, ovlp(n)
      enddo
c      print*,'sum=',SUM(ovlp(1:nps)*ovlp(1:nps)) ! in Ry

      tmpovlp(:,:) = 0.0
      tmpe1r(:,:) = 0.0
      tmpps2(:,:) = 0.0

      tmpovlp(1,1) = 1.0
      tmpe1r(1,1) = SUM(ovlp(1:nps)*ovlp(1:nps)*dble(psen2(1:nps))) ! in Ry
      do n=2,nps
         tmpovlp(n,n) = 1.0
         tmpe1r(n,n) = psen2(ip(n))  ! in Ry
      enddo
      
      do n=2,nps
         m = ip(n)
         tmpovlp(1,n) = ovlp(m)
         tmpovlp(n,1) = tmpovlp(1,n)
         tmpe1r(1,n) = ovlp(m)*dble(psen2(m))
         tmpe1r(n,1) = tmpe1r(1,n)
      enddo

      mintmp(1) = 1
      maxtmp(1) = istopc
      tmpps2(1:istopc,1) = fc(1:istopc)
      do n=2,nps
         m = ip(n)
         i1 = minps2(m)
         i2 = maxps2(m)
         mintmp(n) = i1
         maxtmp(n) = i2
         tmpps2(i1:i2,n) = ps2(i1:i2,m)
      enddo


c      print*,'!  ', ' 11', tmpe1r(1,1),tmpovlp(1,1)
c      print*,'!  ', ' 22', tmpe1r(2,2),tmpovlp(2,2)
c      print*,'!  ', ' 12', tmpe1r(1,2),tmpovlp(1,2)


c     Note: tmpe1r in Ry
      call gsort1(nr,nps,tmpps2,maxtmp,mintmp,tmpovlp,tmpe1r)

c      print*,'!!!', ' 11', tmpe1r(1,1),tmpovlp(1,1)
c      print*,'!!!', ' 22', tmpe1r(2,2),tmpovlp(2,2)
c      print*,'!!!', ' 12', tmpe1r(1,2),tmpovlp(1,2)


      do n=1,nps
         i1 = mintmp(1)
         i2 = min(maxtmp(1),maxps2(n))
         ovlp(n) = SUM(tmpps2(i1:i2,1)*ps2(i1:i2,n)*rmesh(i1:i2,3))
      enddo
c      print*,'sum=',SUM(ovlp(1:nps)*ovlp(1:nps))
      do n=1,nps
         i1 = mintmp(2)
         i2 = min(maxtmp(2),maxps2(n))
         ovlp2(n) = SUM(tmpps2(i1:i2,2)*ps2(i1:i2,n)*rmesh(i1:i2,3))
      enddo
c      print*,'sum=',SUM(ovlp2(1:nps)*ovlp2(1:nps))

      e2r_22 = SUM(ovlp2(1:nps)*ovlp2(1:nps)*dble(psen2(1:nps)))
      e2r_12 = SUM(ovlp(1:nps)*ovlp2(1:nps)*dble(psen2(1:nps)))
c      print*, '!!22:',e2r_22,  '   12:', e2r_12 


      do n=nspm+1,nspm+nps
         m = n - nspm
         lo(n) = l
         ko(n) = m
         i1 = mintmp(m)
         i2 = maxtmp(m)
         minf(n) = i1
         maxf(n) = i2
         fl(i1:i2,n) = tmpps2(i1:i2,m)
         do np=nspm+1,n
            mp = np - nspm
            ortint(n,np) = tmpovlp(m,mp)
            ortint(np,n) = ortint(n,np) 
            e1r(n,np) = tmpe1r(m,mp)/2.0 ! in a.u.
            e1r(np,n) =  e1r(n,np)
            
            if(m .le. 4 ) then
               i1p = mintmp(mp)
               i2p = maxtmp(mp)
               i1 = max(i1,i1p)
               i2 = min(i2,i2p)
               tmp = SUM(tmpps2(i1:i2,m)*tmpps2(i1:i2,mp)*
     >              rmesh(i1:i2,3))
c               write(*,'("n,np,i1,i2,e1r(n,np),ortint(n,np)",4i5,
c     >         1P,4e12.4)') n,np, i1,i2,  e1r(n,np),ortint(n,np),tmp
            endif
         enddo
      enddo

      nspm = nspm + nps

      return
      end
c---------------------------------------------------------------------
      subroutine gsort1(nr,nspm,f,maxf,minf,ortint,e1r)
      
      real*8  ortint(nspm,nspm), e1r(nspm,nspm), sum
!     >   ,r1elk
      real f(nr,nspm)
      integer maxf(nspm),minf(nspm)

      print*,"!!! gsort1()"
c     form overlap array <u_j|v_k>, u_j - old nonorthogonal set,
c     v_k - new set but not yet normalized
c     Only elements with  j >= k  required.
      do j=1,nspm
         do k=1,j
            sum=0d0
            do n=1,k-1
               sum = sum + ortint(k,n)*ortint(j,n)/ortint(n,n)
            end do
            ortint(j,k) = ortint(j,k) - sum
c            write(20,'("j,k =",2I3,", <j|k> =",F10.5)')
c     >         j, k, real(ortint(j,k))
         end do
      end do
c     get new e1r(,):  <u_j|H_1|v_k>, j>=k
      do j=1,nspm
         do k=1,j
            sum=0d0
            do n=1,k-1
               sum = sum + ortint(k,n)*e1r(j,n)/ortint(n,n)
            end do
            e1r(j,k) = e1r(j,k) - sum
         end do
      end do
c     <vb_j|H_1|vb_k>, j>=k
      do j=1,nspm
         do k=1,j
            sum=0d0
            do n=1,j-1
               sum = sum + ortint(j,n)*e1r(n,k)/dsqrt(ortint(n,n))
            end do
            e1r(j,k) = (e1r(j,k)/dsqrt(ortint(k,k))-sum)/
     >           dsqrt(ortint(j,j))
            e1r(k,j) = e1r(j,k)
c            write(20,'("j,k =",2I3,", <j|H_1|k> =",F10.5)')
c     >         j, k, real(e1r(j,k))
         end do
      end do
c

      
c     form new orthonomal set vb_k, f_k = vb_k
      do k=1,nspm
         do i=1,nr
            sum=0d0
            do n=1,k-1
               sum = sum + dble(f(i,n))*ortint(k,n)/dsqrt(ortint(n,n))
            end do
            f(i,k) = (dble(f(i,k)) - sum)/dsqrt(ortint(k,k))
         end do
      end do
      do n=1,nspm
         call minmaxi1(f(1,n),nr,i1,i2)
c         print*,'!!!!n,i1,i2', n,i1,i2
         maxf(n) = i2
         minf(n) = i1
         do np=1,nspm
            ortint(n,np) = 0d0
         end do
         ortint(n,n) = 1d0
      end do
c$$$      do n1=1,nspm
c$$$         do n2=n1,nspm
c$$$            rnorm = 0.0D0
c$$$            if(lo(n1).eq.lo(n2)) rnorm =  r1elk(0,f(1,n1),f(1,n2),
c$$$     >         minf(n1),minf(n2),maxf(n1),maxf(n2))
c$$$            write(20,'("n1,n2 =",2I5,", <k1l|k2l> =",F10.5)') n1,n2,
c$$$     >         real(rnorm)
c$$$         end do 
c$$$      end do 

      return
      end

c-----------------------------------------------------------
      subroutine minmaxi1(f,nr,i1,i2)
      include 'par.f'
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      dimension  f(nr)
      i=1
      do while (i.lt.nr.and.abs(f(i)).lt.regcut)
         i=i+1
      end do
      i1=i
      i=nr
      do while (i.gt.1.and.abs(f(i)).lt.expcut)
         i=i-1
      end do
      i2=i
      return
      end

c-----------------------------------------------------------
      subroutine noble2el(Nmax,namax,nchanmax,eproj,etot,
     >     enion,enlevel)
      include 'par.f'
      parameter(nspmW = 1)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      double precision ortint, e1r
      common /ortog/  ortint(nspmax,nspmax)
      common/hame1/ e1r(nspmax,nspmax)
      common /nstatearray/ nstate(0:lomax,0:1,2)
      common /CcoefsE/ E(KNM)
      double precision E
      common /helium/ llst(KNM), lsst(KNM), lstparity(KNM), np(KNM)
      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      common /Nsarray/ Ns(0:lomax,0:1,komax,2)
      character chan(knm)*3,  opcl*6
      common /charchan/ chan
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer na, nam
      common /corearray/ nicm, ncore(nspmCI)
      common /di_el_core_polarization/ gamma, r0, pol(nmaxr)    
      common /noblgegas_switch/ i_sw_ng    ! set in  cc/main.f
      common /statenumber_str_to_scat/ Nstate_str_to_scat, na_ng(KNM)
      common /major_config/ l1orb(KNM),l2orb(KNM),n1orb(KNM),
     >   n2orb(KNM),config_maj(KNM)

      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >     ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >     npbot(0:lamax),nptop(0:lamax)
      common /MPI_info/myid, ntasks, cnode, ench
      character cnode*3, ench*11

      character label*3, getchchar*2, chnl*2
      character arr(0:5)*1
      data arr /"S","P","D","F","G","H"/

      character atom*20
      integer lmaxhe(0:1,2)
      integer nk1(0:lomax), nk2(0:lomax) 
      integer nkng(nspm)      
      real*8, dimension(:,:), allocatable:: CI, H, bb
      real*8, dimension(:), allocatable:: en
      common/meshrr/ nr,rmesh(maxr,3)
      real*8, dimension(:,:), allocatable:: Cst
      double precision, dimension(:,:,:),allocatable::  Ctmposc


      real  psi(maxr)   ! for getchinfo() 
      
      print*,'Start noble2el'


      na(:,:) = 0
      nam(:) = 0

c     find maximum value of l
      lspmax = -1
      do n=1,nspm
c         print*, n, lo(n)
         lspmax = max(lspmax,lo(n))
      enddo

      print*,'n_ng=',n_ng, ', l_ng=',l_ng, ', lspmax=',lspmax

      if(lspmax .eq. -1) then
         print*,'lspmax=-1, looks like no s.p. basis was calculated... '
         print*, 'nspm=', nspm
         stop
      endif

c     In a first go make s.p.orbitals
c     from diagonalization of Ne Hamiltonian (note: choice of 
c     singlet or triplet Ne states to create s.p. orbitals, or even both 
c     of them: like in Ca e2e code).
c      In a second go make target states.

      if (cnode.eq.'') then
         open(4,file='F6',action='WRITE')
      else 
         open(4,file=cnode//'F6'//ench,action='WRITE')
      endif 
      open(3,file='F5',iostat=iostatF5,action='READ')
      if(iostatF5.ne.0) then
         print*, '********  File  F5  was not found'
         stop
      end if
      open(10,file=cnode//'he-states'//ench,action='WRITE')

      print*, "!!! First make orbitals: FC orbitals",
     >     " for singlet states to form a basis."

      make_orbitals = 1  !  first go is to make orbitals
 110  read(3,'(a20)') atom
      write(4,'(a20)') atom


c     find maximum value of l
      lspmax = -1
      do n=1,nspm
c         print*, n, lo(n)
         lspmax = max(lspmax,lo(n))
      enddo

      do l=0,lspmax
         ntmp = 0
         do n=1,nspm
            if(lo(n) .eq. l) then
               ntmp = ntmp+1
            endif
         enddo
         print*, 'l=',l, ', number of orbitals=', ntmp
         write(4,'("l=",I5,", number of orbitals=",I5)') l, ntmp
         if(l .eq. l_ng ) then
            nsp_l_ng = ntmp   ! number of orbitals for l = l_ng
         endif
      enddo


      
      lmaxhe(:,:) = -1
      nstate(:,:,:) = 0
      Ns(:,:,:,:) = -1

      do ip=1,2
         do is=0,1
            read(3,*) lmaxhe(is,ip)
            write(4,'("lmaxhe(s=",I1,",ip=",I1,") =",I3)')
     >         is, ip, lmaxhe(is,ip)
            read(3,*) (nstate(i,is,ip), i=0,lmaxhe(is,ip))
            do i=0,lmaxhe(is,ip)
               if(nstate(i,is,ip) .lt. 0 ) nstate(i,is,ip) = 0
            enddo
            write(4,'("nstate(s=",I1,",ip=",I1,")=",10I3)')
     >         is, ip, (nstate(i,is,ip),i=0,lmaxhe(is,ip))
         end do
      end do
      
      read(3,*) inc, i_dipmom, gamma, r0
      write(4,'("inc =",I3)') inc
      write(4,'("i_dipmom =",I3)') i_dipmom
      write(4,'("Di-electron pol.pot.: gamma =",F10.5," r0 =",F10.5)') 
     >     gamma, r0

c     if i_dipmom = 0  no oscillator strengths are calculated
c     if i_dipmom = 1  absorption oscillator strengths are calculated
c     if i_dipmom = 2  lp >= l oscillator strengths are calculated
c     if i_dipmom = 3  dipole polarazability of the states (l=0)  is calculated
c                      in this case only S and P states should be included.
c     if i_dipmom = 4  singlet triplet mixing is calculated
c     if i_dipmom = 5  singlet triplet mixing and osc. strength are calculated
      i_stmix = 0
      if(i_dipmom .eq. 4) then
         i_dipmom = 1
         i_stmix = 1
      endif
      if(i_dipmom .eq. 5) then
         i_dipmom = 1
         i_stmix = 2
      endif

c     Set up the di-electron polarization array
      if(gamma.ne.0.0) 
     >     call di_el_pol
c

c     here  Nmax  is determined (number of states).
      isum = 0
      do ip=1,2
         do is=0,1
            do l=0,lmaxhe(is,ip)
               isum = isum + nstate(l,is,ip)
c               print*, ip,is,l,isum
            end do
         end do
      end do
      Nmax=isum
      Nstate_str_to_scat = Nmax
      write(4,'("Nmax =",I3)') Nmax
      print*, 'Nmax =', Nmax
      if (Nmax.gt.KNM) then
         print*,'Have too many target states, need to increse KNM to',
     >      ' at least',Nmax
         stop 'Have too many target states, need to increse KNM.'
      endif 

      if(allocated(Cst)) deallocate(Cst)
      allocate(Cst(Nmax,nspm))
      Cst(:,:) = 0.0

c     here we decide how many times helium cofiguration interaction
c     program should be called
      nslm = 0
      do ip=1,2
         do is=0,1
            do l=0,lmaxhe(is,ip)
               if(nstate(l,is,ip).gt.0) then
                  nslm = nslm + 1
               end if
            end do
         end do
      end do
      print*,'nslm=',nslm
c
      nicm = 1
      ncore(1) = 1
      ic = ncore(1)   ! core elctron first: 2p for Ne.
      lfcore = l_ng
      lfparity = (-1)**(l_ng)
      sfcore = 0.5
c     needed to have singlet S only for p6: this is basically due to fractional parantage coef.
      ll_gs = 0  
      s_gs = 0
      ipar_gs = 1

      ninc = 1   ! need for pwrite()
      linc = ll_gs
c!

      Nst = 0  ! current number of states
      do nsll=1,nslm
         
         print*
         write(4,*)
         print*,'nsll=',nsll

         call  inputdata(ll,ls,lparity,l1max,l2max,nk1,nk2)
         

         la = ll
         ls = ls
         lpar = lparity
         ip = 1
         if((-1)**ll.ne.lpar) ip = 2
         if(nstate(ll,ls,ip).le.0) cycle

         nstmax = nstate(ll,ls,ip)
        
         if(ls.eq.0) then
            sa = abs(sfcore - 0.5)
         elseif(ls.eq.1) then
            sa = sfcore + 0.5
         else
            print*,'noblegas.f: wrong value for spin index: ls=',ls
            stop
         endif
         
c     set up list of configurations
         nc = 0
         do nsp=1,nspm
            lnsp = lo(nsp)
            if(nsp.eq.ic .and. (la.ne.ll_gs .or. sa.ne.s_gs)) cycle
            if(lnsp .gt. l1max) cycle
            if(la.ge.abs(lfcore-lnsp) .and. la.le.(lfcore+lnsp)) then
               if(nint(2.0*sa).ge.abs(nint(2.0*(sfcore-0.5))) 
     >              .or. nint(2.0*sa).le.nint(2.0*(sfcore+0.5))) then
                  if(lpar.eq.lfparity*(-1)**(lnsp)) then
            
                     if(nk1(lnsp) .ge. ko(nsp)) then
                        nc = nc + 1
                        nkng(nc) = nsp
                        lnsp_1 = lnsp
                     endif
         
                  endif
               endif
            endif
         enddo
         ncm = nc

         print*,'la=',la,', sa=',sa,', lpar=',lpar,
     >        ', ncm=',ncm, ', nspm=',nspm
         if(ncm .eq. 0) then
            print*, 'Stop: ncm=0'
            stop
         endif

         if(make_orbitals .eq. 1) then
            lnsp = lnsp_1 

            i_subtruct_one = 0
            if(la .eq. 0 .and. ls .eq. 0) then
               i_subtruct_one = 1
C     "-1" is to account the core orbital (1s for He, 2p for Ne,...).
C     For example, with 2  s  orbitals there will be only 2 configurations -> states 
C     for singlet S symmetry but for two states there will be three orbitals (1s, 1s', 2s')
            endif

            itmp = nstmax-i_subtruct_one + (nabot(lnsp) - lnsp -1)
            if( ncm .lt. itmp ) then
               print*,'noblgas.f: too many states ordered,',
     >              ' increase number of Lag. func. in ccc.in file'
               print*, ' or reduce number of states in F5 file'
               print*, 'number of configurations minus ',
     >              'core orbitals:ncm=',ncm - (nabot(lnsp) - lnsp -1)
               print*, 'number of states ordered:', nstmax 
               stop
            endif
         else
            if(ncm.lt.nstmax) then
               print*,'Number of confugurations is less than ',
     >            'the number of ordered states'               
               print*,'ncm=', ncm, ',  nstmax=',nstmax
               print*,'decrease number of states in F5 file'
               stop
            endif

         endif                                
         
         allocate(H(ncm,ncm),bb(ncm,ncm))
         allocate(CI(ncm,ncm),en(ncm))
         CI(:,:) = 0.0

         call H12_ng(la,ls,ncm,nkng(1:ncm),H,bb)

c     start diagonalization
         enionry = 0.0          ! 2.*e1r(i_ng,i_ng)
         call eigv_ng(enionry,ncm,H,bb,CI,en)


         
         if(make_orbitals .eq. 1) then
            noffset = (nabot(lnsp) - lnsp) - 1 ! need to remove core orbitals
         else
            noffset = 0
         endif

c     make label
         if(ls .eq. 0) then
            if(ip .eq. 1) then
               label(1:1) = 's'
            elseif(ip .eq. 2) then
               label(1:1) = 'S'
            endif
         elseif(ls .eq. 1) then
            if(ip .eq. 1) then
               label(1:1) = 't'
            elseif(ip .eq. 2) then
               label(1:1) = 'T'
            endif
         else
            print*, 'nonlegas.f: Wrong value of spin ls=',ls
            stop
         endif
c$$$         label(3:3) = arr(la)
         chnl = getchchar(1,la)
         label(3:3) = chnl(2:2)
c     

         print*, 'noffset = ', noffset         
         write(4,'("noffset = ",I5)')  noffset
         
         nstmax_can_order = ncm - noffset
         print*,'Maximum number of states you can order is:',
     >        nstmax_can_order,
     >         " |  have actually ordered:",nstmax,
     >      nstmax_can_order-nstmax
         write(4,'("Maximum number of states you can order is:",
     >        I5, " |  have actually ordered:",2i5)') 
     >        nstmax_can_order, nstmax, nstmax_can_order-nstmax

         do n=1+noffset,nstmax+noffset
            m = n 
            Nst = Nst + 1
            label(2:2) = char(ichar('0') + n)
            chan(Nst) = label
            E(Nst) = en(m)
!            print*,"!!!", Nst,m,E(Nst)
            Ns(la,ls,n,ip) = Nst
            llst(Nst) = la
            lsst(Nst) = ls
            lstparity(Nst) = lpar
            nam(Nst) = ncm
            do nc=1,ncm
               nsp = nkng(nc)
               na(nc,Nst) = nsp
               Cst(Nst,nsp) = CI(nc,m)
c               print*,la,ls,Nst,nc,CI(nc,m)
            enddo
c     fix sign of Cst
            tmp = SUM(Cst(Nst,1:nspm))
            if(tmp. lt.0.0) then
               Cst(Nst,1:nspm) = -Cst(Nst,1:nspm)
            endif

         enddo
c         
         deallocate(CI,en,H,bb)
      enddo

      if(make_orbitals .eq. 1 .or. inc .eq. 1) then
         call rearange_ng(nr,nspm,Nmax,lo,ko,fl,minf,maxf,Cst)
      endif


      if(make_orbitals .eq. 1) then
         call gsort1(nr,nspm,fl(1:nr,1:nspm),
     >        maxf(1:nspm),minf(1:nspm),
     >        ortint(1:nspm,1:nspm),e1r(1:nspm,1:nspm))        
      endif
      
      
c      print*,'!! make_orbitals=',make_orbitals
      if(make_orbitals .eq. 1) then
         deallocate(Cst)
         make_orbitals = 0      ! orbitals are made go back for final run
         print*
         write(4,*)
         print*, "!!! FC orbitals are made go back for final run"
         print*
         write(4,*)
         goto 110
      endif

c     Compact Cst aray if inc .eq. 0, if inc = 1 then Cst array is already compacted.
      if(inc .eq. 0) then
         do Nst = 1,Nmax
            Cst(Nst,1:nam(Nst)) = Cst(Nst,na(1:nam(Nst),Nst))
         enddo
      endif

      call find_namax(Nmax,namax)
      write(4,*)
      write(4,'("namax =",I5)') namax
      call noblegas_make_new_C(Cst,Nmax,nspm,namax)
      deallocate(Cst)

      call energyorder(Nmax, KNM, nspmCI, llst, lsst, lstparity, np, 
     >   na, nam, nicm, E, chan, Ns,lomax,komax)

      if(atom.eq.'argon') then
         call  dataAr
      else  if(atom.eq.'neon' .or. atom.eq.'Neon') then
         call  dataNe
      endif

      call noblegas_find_major_config(Nmax,ll_gs,s_gs,ipar_gs)
      call print_energy(Nmax,E,enionry,atom)
      call printCcoef_ng(Nmax,E)


!      i_dipmom = 1
      allocate(Ctmposc(Nmax+1,nspmW,nspmW))
      if(i_dipmom .ne. 0) then 
         call dipmom(i_dipmom,atom,Nmax,nspmW,Ctmposc,E,fl,minf,
     >        maxf,lo,enionry)
      endif
      deallocate(Ctmposc)

      i_present_ng = i_sw_ng
      if(i_stmix .ne. 0) then
         call Spin_orbit(Nmax,E,enionry,i_stmix,i_present_ng)
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Set up for getchinfo() to work properly
      nchanmax = 0
      lg = 0

      lmax_tmp = -1
      lmin_tmp = 100
      do nst = 1, Nmax
         lmax_tmp = max(lmax_tmp,llst(nst))
         lmin_tmp = min(lmin_tmp,llst(nst))
      enddo
! set up arrays to get correct energies in getchinfo na_ng()  and enpsinb()
      do la=lmin_tmp,lmax_tmp
         
         icount = 0
         do nst = 1,Nmax
            if(llst(nst) .ne. la) cycle
            icount = icount + 1
            na_ng(nst) = icount
            enpsinb(na_ng(nst),la) = E(nst) * 2.0 - enionry         
            nchanmax = nchanmax + la + 1
!                  print*, '***la, nst,na_ng(nst):', la, nst,na_ng(nst)
         enddo
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      etot = eproj + E(1)*2.0



      print*, 'Finish structure'
      lg = 0
      ry = 13.6058

      ovlp_st = 0.0
      do nst = 1, Nmax
c     find fc. orbital in fc. app. for state nst

         na_st = na_ng(nst)
         latom = llst(nst)

         if (etot.ge.enpsinb(na_st,latom)) then
            onshellk = sqrt(etot-enpsinb(na_st,latom))
            opcl = 'open'
            if (enpsinb(na_st,latom).gt.0.0) then
c     put here e2e overlaps: see  helium/main.f
               ovlpnl_st = -1.0
            else
               ovlpnl_st = 1.0           
            endif 
c     put here e2e overlaps: see  heliuym/main.f
         else
            opcl = 'closed'
            onshellk = -sqrt(enpsinb(na_st,latom)-etot)
            ovlp_st = 0.0
         endif 
         elevel = (enpsinb(na_st,latom) * ry + enion) * 
     >      enlevel / enion
         print '(i3,a4,f12.4," Ry",f12.4," eV",f12.1," / cm",a8,f9.4,
     >      f8.2)', 
     >     nst, chan(nst),enpsinb(na_st,latom),
     >        enpsinb(na_st,latom)*ry,
     >      elevel, opcl, onshellk, ovlp_st 

      enddo 
      write(*,'("finish structure calculation")')
      write(*,'("**********************************************")')
c
C  EDELTA is used to shift the energy levels
      edelta = - enion - enpsinb(ninc,linc) * ry
      print*,'Error in ground state energy and ',
     >   'effective incident energy (eV):', edelta, etot * ry + enion


      print*,'Finish noble2el'
c      stop
      return
      end

c-----------------------------------------------------------
      subroutine H12_ng(la,ls,ncm,nkng,H,bb)
      include 'par.f'
      double precision e1r
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common/hame1/ e1r(nspmax,nspmax)
      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      
      integer ncm,la,ls      
      integer nkng(ncm)
      real*8  H(ncm,ncm), bb(ncm,ncm)
      real*8 coef, sumk, ang, radint, coef_a

      real c1, c2, c3, c123
      real rla, rl_ng, rlnsp, rlnspp, rk
      
      double precision RAD
      real CGC0, COF6J

      nfcel = 4*l_ng + 1   ! number of 'frozen'-core electrons (5 : for Ne 2p^5)
      coef_a = -1.0 
      rl_ng = l_ng


      do n=1,nspm
         if(lo(n) .eq. l_ng) then
            i_ng = n
            exit
         endif
      enddo
c      print*,'i_ng=',i_ng
      
      itest = 0 ! itest=1 for helium test case, itest=0 for Ne
      if(itest .eq. 1) then
         nfcel = 1              ! test case of Ne V
         coef_a = 1.0          ! fix sign to +1 for k!=0 for 2 electron config
      endif


      rla = la


      coef = 1.0     
      do nc=1,ncm
         nsp = nkng(nc)
         lnsp = lo(nsp)
         rlnsp = lnsp
c         print*,nsp,lo(nsp),e1r(nsp,nsp)
         do ncp = 1,nc
            nspp =nkng(ncp)
            lnspp = lo(nspp)
            rlnspp = lnspp

            if(nc.eq.1 .and. ncp.eq.1) then 
c               print*,'!!! nc,ncp', nc,ncp,e1r(nsp,nspp)
            endif


c     get one-electron term
            if( (nsp .eq. i_ng .and. nspp .ne. i_ng) .or. 
     >           (nsp .ne. i_ng .and. nspp .eq. i_ng)) then
               coef = sqrt(nfcel + 1.0)
            else
               coef = 1.0
            endif
            H(nc,ncp) = coef * e1r(nsp,nspp)
            if(itest.eq.1 .and. nsp.eq.nspp) then
               H(nc,ncp) = H(nc,ncp) + e1r(i_ng,i_ng)
            endif

c     get direct two-electron term : only even multipoles are allowed
            maxk = min(2*l_ng,lnsp+lnspp)
            sumk = 0.0
            do k=0,maxk,2
               rk = k   
               c1 = CGC0(rl_ng,rk,rl_ng)
               c2 = CGC0(rlnspp,rk,rlnsp)
               c3 = COF6J(rl_ng,rlnsp,rla,rlnspp,rl_ng,rk)
               ang = sqrt(2*l_ng+1.0)*CGC0(rl_ng,rk,rl_ng) * 
     >              sqrt(2*lnspp+1.0)*CGC0(rlnspp,rk,rlnsp) * 
     >              COF6J(rl_ng,rlnsp,rla,rlnspp,rl_ng,rk) * 
     >              dble((-1)**(lnsp+l_ng+la))
               coef_a = -1.0 
               if(itest .eq. 1) then
                  coef_a = 1.0 
               endif
               if(k .eq. 0) coef_a = nfcel
               ang = coef_a * ang
               if(ang .eq. 0.0) cycle
               
               n1 = i_ng
               n1p = i_ng
               n2 = nsp
               n2p = nspp
               radint = RAD(k,fl(1,n1),fl(1,n1p),
     >              fl(1,n2),fl(1,n2p),minf(n1),minf(n1p),minf(n2),
     >              minf(n2p),maxf(n1),maxf(n1p),maxf(n2),maxf(n2p))

               sumk = sumk + ang * radint
c               if(nc.eq.1 .and. ncp.eq.1) then
c                  print*,'k,ang,radint',k,ang,radint
c              endif
            enddo
            H(nc,ncp) =  H(nc,ncp) + coef*sumk
            H(ncp,nc) =  H(nc,ncp) 
c            if(nc .le. 2) then
c               print*,'direct: nc,ncp,H(nc,ncp)=',nc,ncp,H(nc,ncp) 
c            endif

!
c            H(ncp,nc) = H(nc,ncp)
c            cycle
!

c     get exchange two-electron term : only required when both nsp != 1 and nspp != 1
            if(nsp .eq. i_ng .or. nspp .eq. i_ng) then
c               print*,'noexch'
               cycle
            endif
            
            maxk = min(l_ng+lnsp,l_ng+lnspp)
            mink = max(abs(l_ng-lnsp),abs(l_ng-lnspp))
            sumk = 0.0
            do k=mink,maxk
               rk = k               
               
               c1 = 0.0
               if(lnsp .eq. lnspp) c1 = 1.0/(2*lnsp+1.0)

               c2 = 0.0
               if(ls .eq. 0) c2 = 1.0

               c3 = 0.0
               if(k .eq. la) c3 = 1.0

               c123 = c1 - 2.0*c2*c3*(-1)**(lnsp+lnspp)/(2*la+1.0)
               if(c123 .eq. 0.0) cycle

               ang = -sqrt(2*lnsp+1.0)*CGC0(rlnsp,rk,rl_ng) * 
     >              sqrt(2*lnspp+1.0)*CGC0(rlnspp,rk,rl_ng) 
               
               if(itest .eq. 1) then
                  tmpcof6j = COF6J(rl_ng,rlnsp,rla,rl_ng,rlnspp,rk)
                  ang = tmpcof6j * ang
               endif

               if(ang .eq. 0.0) cycle

               if(itest .eq. 0) then
                  ang = c123 * ang
               else
                  ang = ang  * (-1)**(ls)
               endif

               
               n1 = i_ng
               n1p = nspp
               n2 = nsp
               n2p = i_ng
               radint = RAD(k,fl(1,n1),fl(1,n1p),
     >              fl(1,n2),fl(1,n2p),minf(n1),minf(n1p),minf(n2),
     >              minf(n2p),maxf(n1),maxf(n1p),maxf(n2),maxf(n2p))

               sumk = sumk + ang * radint
c               print*,'exch: k,ang,radint',k,ang,radint
            enddo
            H(nc,ncp) =  H(nc,ncp) + coef * sumk

            H(ncp,nc) =  H(nc,ncp) 
c            print*,'exch: nc,ncp,H(nc,ncp)=',nc,ncp,H(nc,ncp) 

         enddo
         
      enddo

      do nc=1,-ncm
         write(*,'(1P,100E15.5)') (H(nc,ncp), ncp=1,ncm) 
      enddo

      bb(:,:) = 0.0
      do n=1,ncm     
         bb(n,n) = 1.0
      enddo

      return
      end

c-----------------------------------------------------------------
c*************** Find the eigenvalues and eigenvectors  **********
c********************* of the Hamiltonian matrix *****************
c-----------------------------------------------------------------
      subroutine eigv_ng(enionry,ncm,H,bb,CI,w)
      include 'par.f'
      double precision sum
      double precision H(ncm,ncm),bb(ncm,ncm),w(ncm),
     >   CI(ncm,ncm),fv1(ncm),fv2(ncm)
      write(*,'("start diagonalization")')

      matz=2
      call  rsg(ncm,ncm,H,bb,w,matz,CI,fv1,fv2,ierr)
      write(*,'("ierr =",I3)') ierr
      write(*,'("eigenvalues in a.u.")')
      write(*,'(5F12.6)') (real(w(i)), i=1,ncm)

c      write(*,'("eigenvalues in Ry")')
c      write(*,'(5F12.6)') (real(2d0*w(i)), i=1,ncm)

      write(*,'("eigenvalues - enionry/2.  -  in eV, enionry=",F10.5,
     >     " Ry")') enionry
      write(*,'(5F12.6)') ((real(w(i))-enionry/2.0)*27.2116,i=1,ncm)
      write(*,'("eigenvectors can be printed")')
      do i=1,-2 ! ncm
         write(*,*)
         write(*,'(5F12.6)') (real(CI(i,j)), j=1,ncm)
      end do
      write(4,'("number of configurations: ncm=",I10)') ncm
      write(4,'("eigenvalues - enionry/2.  -  in eV, enionry=",F10.5,
     >     " Ry")') enionry
      write(4,'(5F12.6)') ((real(w(i))-enionry/2.0)*27.2116,i=1,ncm)
      write(4,'("eigenvectors can be printed")')
      do i=1,-2 ! ncm
         write(4,*)
         write(4,'(5F12.6)') (real(CI(i,j)), j=1,ncm)
      end do
      return

      end
c-----------------------------------------------------------------
      subroutine rearange_ng(nr,nspm,Nmax,lo,ko,fl,minf,maxf,C)
      include 'par.f'
      integer  nr, nspm, Nmax
      integer lo(nspmax), ko(nspmax)
      real  fl(nmaxr,nspmax)
      integer maxf(nspmax),minf(nspmax)      
      real*8 C(Nmax,nspm)

      common /meshrr/ nrmeshr, gridr(nmaxr,3)

      double precision ortint, e1r
      common /ortog/  ortint(nspmax,nspmax)
      common/hame1/ e1r(nspmax,nspmax)

      common /corearray/ nicm, ncore(nspmCI)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer na, nam

      integer  nlstmax(Nmax), ltmp(Nmax,maxval(lo(1:nspm)))

      integer, dimension(:), allocatable:: ngiveN
      integer, dimension(:), allocatable:: loN
      integer, dimension(:,:), allocatable:: iorb
      integer, dimension(:), allocatable:: iorbmax
      real*8, dimension(:), allocatable:: p
      real*8, dimension(:,:), allocatable:: e1rp, ortp

      real*8 ovlpmin

      print*,'Start rearange_ng'

      ovlpmin=1.0d-7

      ic = ncore(1)! the first core orbital, 2p for Ne.

      imaxr = 0
      do n=1,nspm
         i1 = minf(n)
         i2 = maxf(n)
         imaxr = max(imaxr,i2)
         if(i1 .ne. 1)  fl(1:i1-1,n) = 0.0
         if(i2 .ne. nr) fl(i2+1:nr,n) = 0.0
      enddo

      nsp_tot = ic               ! include the first core orbital

      ltmp(:,:) = -1
      do N=1,Nmax
c     Here we find how many f.c. orbitals are built on the l^w core
c     nlstmax is the number and ltmp(nlct) is the ang.mom. of these orbitals.
         nlct = 0
         do nsp=ic+1,nspm          ! start from 2 as nsp=1 is the core orbital
            if(C(N,nsp).ne.0d0) then                              
               icheck = 1       ! check if new frozen-core orbital (new l) need to be added 
               do il=1,max(nlct,1)
                  if(ltmp(N,il).eq.lo(nsp)) icheck = 0
               end do
               if(icheck.eq.1) then
                  nlct = nlct + 1
                  ltmp(N,nlct) = lo(nsp)
c                  print*, N,nsp,lo(nsp),nlct,' ! ',
c     >                 (ltmp(N,i), i=1,nlct)
               end if                  
            end if

         end do
         nlstmax(N) = nlct
         nsp_tot = nsp_tot  + nlstmax(N)
      end do
      if(nsp_tot.gt.nspmax) then
         print*,'stop in rearrange: nsp_tot > nspmax'       
         print*,'nsp_tot =', nsp_tot, ',  nspmax =',  nspmax
         stop
      end if

      write(*,'("Number of frozen core functions is  nsp_tot =",I6)') 
     >      nsp_tot


      allocate(p(nsp_tot))
      allocate(ngiveN(nsp_tot),loN(nsp_tot))
      allocate(iorbmax(nsp_tot),iorb(nsp_tot,nspm))
      allocate(ortp(nsp_tot,nsp_tot),e1rp(nsp_tot,nsp_tot))
      p(:) = 0.0
      iorb(:,:) = 0
      iorbmax(:) = 0
      ortp(:,:) = 0.0
      e1rp(:,:) = 0.0

c     take care of core orbital
      iorbmax(ic) = ic
      iorb(ic,ic) = ic
      ngiveN(1) = 0
      loN(1) = lo(ic)

      nspN = 1      
      do N=1,Nmax
         do nlst=1,nlstmax(N)
            nspN = nspN + 1
            ngiveN(nspN) = N
            loN(nspN) = ltmp(N,nlst)
            ict = 0
            do nsp=ic+1,nspm
               if(lo(nsp) .eq. loN(nspN)) then
                  if(C(N,nsp).ne.0d0) then  
                     ict = ict + 1
                     iorb(nspN,ict) = nsp
                  endif
               endif
            enddo
            iorbmax(nspN) = ict
         enddo         

      enddo


      do i=1,imaxr
         do nspN=ic+1,nsp_tot
            Nst = ngiveN(nspN)
            sume = 0.0D0
            do ict=1,iorbmax(nspN)
               nsp = iorb(nspN,ict)
               if(i.ge.minf(nsp).and.i.le.maxf(nsp))
     >              sume = sume + C(Nst,nsp)*dble(fl(i,nsp))
            end do
            p(nspN) = sume
         enddo
         do nspN=ic+1,nsp_tot
            fl(i,nspN) = p(nspN)
         end do
         
      enddo
      do nspN=ic+1,nsp_tot
         call minmaxi1(fl(1,nspN),nr,i1,i2)
         maxf(nspN) = i2
         minf(nspN) = i1
      end do
 

c     form new array e1r(nsp1,nsp2) and ortint(nsp1,nsp2)
      e1rp(ic,ic) = e1r(ic,ic)
      ortp(ic,ic) = ortint(ic,ic)

      nsp2 = 1
      nspN2 = 1
      do nspN1=2,nsp_tot
         N1 = ngiveN(nspN1)
         sume = 0d0
         sumort = 0d0
         if(loN(nspN1).eq.loN(nspN2)) then
            do ict1=1,iorbmax(nspN1)
               nsp1 = iorb(nspN1,ict1)
               
               sume = sume + C(N1,nsp1)*e1r(nsp1,nsp2)
               
               sumort = sumort + C(N1,nsp1)*ortint(nsp1,nsp2)
            end do
         end if
         e1rp(nspN1,nspN2) = sume
         e1rp(nspN2,nspN1) = sume
         ortp(nspN1,nspN2) = sumort
         ortp(nspN2,nspN1) = sumort
      end do
      

      do nspN1=ic+1,nsp_tot
         N1 = ngiveN(nspN1)
         do nspN2=nspN1,nsp_tot
            N2 = ngiveN(nspN2)
            sume = 0d0
            sumort = 0d0
            if(loN(nspN1).eq.loN(nspN2)) then
               do ict1=1,iorbmax(nspN1)
                  nsp1 = iorb(nspN1,ict1)
                  do ict2=1,iorbmax(nspN2)
                     nsp2 = iorb(nspN2,ict2)

                     sume = sume + C(N1,nsp1)*C(N2,nsp2)*
     >                    e1r(nsp1,nsp2)
                     
                     sumort = sumort + C(N1,nsp1)*C(N2,nsp2)*
     >                    ortint(nsp1,nsp2)
                  end do
               end do
            end if
            e1rp(nspN1,nspN2) = sume
            e1rp(nspN2,nspN1) = sume
            ortp(nspN1,nspN2) = sumort
            ortp(nspN2,nspN1) = sumort
         end do
      end do

      nspm = nsp_tot
      do nsp1=1,nspm
         do nsp2=1,nsp1
            if(nspN1.ne.nspN2) then
               if(abs(ortp(nsp1,nsp2)) .lt. ovlpmin) then
                  ortp(nsp1,nsp2) = 0.0
                  ortp(nsp2,nsp1) = 0.0
               endif
            else
               if(abs(abs(ortp(nsp1,nsp2))-1.0d0) .lt. ovlpmin) then
                  ortp(nsp1,nsp2) = 1.0
               endif
            endif
         enddo
      enddo
      ortint(1:nspm,1:nspm) = ortp(1:nspm,1:nspm)
      e1r(1:nspm,1:nspm) = e1rp(1:nspm,1:nspm)
      lo(1:nspm) = loN(1:nspm)
      
      do l=0,lomax
         k=0
         do nsp=1,nspm
            if(lo(nsp) .eq. l) then
               k = k + 1
               ko(nsp) = k
            endif
         enddo
      enddo

      do N=1,Nmax
         A = C(N,ic)
         C(N,:) = 0.0
         C(N,ic) = A
         i = 0
         nam(N) =  nlstmax(N)         
         if(A .ne. 0.0) then
            i = ic
            nam(N) =  nlstmax(N) + 1
         endif
         do nsp=ic+1,nspm
            if(ngiveN(nsp).eq.N) then
               i = i + 1
               na(i,N) = nsp
               C(N,i) = 1.0
            endif
         enddo
         if(i .ne. nam(N)) then
            print*,'noblegas.f: i .ne. nam(N):',i, nam(N), N
            stop
         endif
         
      enddo


      do n=1,-nspm

         i1 = minf(n)
         i2 = maxf(n)

         tmp = sum(fl(i1:i2,n)*fl(i1:i2,n)*gridr(i1:i2,3))

         print*, 'l,N,n, ortint(n,n):', 
     >        lo(n), ngiveN(n),n, ortint(n,n),tmp


      enddo


      print*,'Finish  rearange_ng'
      return
      end
c-------------------------------------------------------------------
      subroutine printCcoef_ng(Nmax,E)
      use CI_MODULE             ! module with CI-array C
      include 'par.f'
      parameter (ic_20 = 20)
      double precision E(KNM), ortint, tmp, tmp1, sum
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      character chan(knm)*3
      common /charchan/ chan
      character lorb(0:10)*1
      data lorb /"s","p","d","f","g","h","i",
     >   "j","k","l","m"/      
      common /major_config/ l1orb(KNM),l2orb(KNM),n1orb(KNM),
     >   n2orb(KNM),config_maj(KNM)
      common /corearray/ nicm, ncore(nspmCI)
c     
      write(10,'("Nmax =",I3)') Nmax

      do N=1,Nmax
         write(10,'(A3,", N=",I3,"  l=",I3,"  s=",I3," parity=",I3,
     >        ", energy=",F12.5,", major config.:(",2(I2,A1),")",
     >        F10.4)')
     >      chan(N), N, ll(N), ls(N), lparity(N), E(N),
     >        n1orb(N), lorb(l1orb(N)), n2orb(N), lorb(l2orb(N)),
     >        config_maj(N)
         
         write(10,'("   n1   n2    l1  l2",8X,"-",7X,"CI weights")')

         sum = 0d0
         
         j = ncore(1)
         n1 = j
         do i=1,nam(N)
            tmp = C(N,j,i)
            n2 = na(i,N)
c            print*, '!', N, n1, n2, i, j
            if(n1.le.n2.and.tmp.ne.0.0D0) then                 
               tmp1 = tmp*dsqrt(ortint(n1,n1)*ortint(n2,n2))
               sum = sum + tmp1*tmp1     
               
               write(10,'(2I5,2X,2I4,2X,2E12.4)')
     >              n1,n2,lo(n1),lo(n2),
     >              real(tmp1),real(tmp1*tmp1)
            end if
         end do
         write(10,'("sum =",F10.6)') sum
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine noblegas_make_new_C(C_old,Nmax,nspm,namax)
      use CI_MODULE             ! module with CI-array C
      include 'par.f' 
      common /corearray/ nicm, ncore(nspmCI)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      real*8  C_old(Nmax,nspm)
      character chan(knm)*3
      common /charchan/ chan

     
      nic = 1
      ic = ncore(nic)
      
      itest = 1
      if(itest .eq. 1) then

         do N=1,Nmax
            
c     check if the core orbital  'ic'  is included in the na(i,Nst)  array:
            icheck = 0
            do nc=1,nam(N)
               if(ic .eq. na(nc,N)) then
                  icheck = 1
                  if(nc .ne. 1) then
                     print*, 'noblegas/f: possible problem: nc=',nc
                     stop
                  endif
                  exit
               endif
            enddo
            if(icheck .eq. 0) then
               nam(N) = nam(N) + 1
               do nc=nam(N),2,-1
                  na(nc,N) = na(nc-1,N)
                  C_old(N,nc) = C_old(N,nc-1)
               enddo
               na(1,N) = ic
               C_old(N,1) = 0.0
            endif
            namax = max(namax, nam(N))
         enddo
      endif

      allocate(C(Nmax,namax,namax))
      C(:,:,:) = 0.0

      do N=1,Nmax
         do j=1,nam(N)
            C(N,nic,j) =  C_old(N,j)
c            print'(i3,A5,4i3,e15.4)', N,  chan(N), nam(N), 
c     >           na(j,N),nic,j, real(C(N,nic,j))

         end do      
      end do      
      
      return
      end


c-------------------------------------------------------------------
      subroutine noblegas_find_major_config(Nmax,ll_gs,s_gs,ipar_gs)
      use CI_MODULE             ! module with CI-array C
      include 'par.f'
      double precision   ortint
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /major_config/ l1orb(KNM),l2orb(KNM),n1orb(KNM),
     >   n2orb(KNM),config_maj(KNM)
      common /corearray/ nicm, ncore(nspmCI)
      integer num(0:lomax,nspmCI)
      real  c_big(nspmCI), big(0:lomax)
      integer  l_big(nspmCI)

c     Find major configuration for each state.
      write(10,'("Find major configuration for each state.")')
      do il=0,lamax
         do is=0,1
            do ip=-1,1,2
               
               do nic=1,nicm
                  do l_orb=0,lamax
                     num(l_orb,nic) = 0
                  end do
               end do

               do N=1,Nmax
                  ila = ll(N)
                  isa = ls(N)
                  ipa = lparity(N)
                  if(ila.eq.il.and.isa.eq.is.and.ipa.eq.ip) then
                     do nic=1,nicm
                        c_big(nic) = 0d0
                        l_big(nic) = -1
                     end do
                     do nic=1,nicm ! go through all core orbitals
                        do l=0,lomax
                           big(l) = 0d0
                        end do
                        i = nic
                        n1 = ncore(nic)
                        if(ncore(nic).eq.n1) then ! n1 is a core orbital
                           do j=i,nam(N)
                              n2 = na(j,N)
                              ln2 = lo(n2)
                              tmp = C(N,i,j)                           
                              tmp=tmp*tmp*ortint(n1,n1)*ortint(n2,n2)
                              big(ln2) = big(ln2) + tmp
c                              print*,'N,n1,n2,big(l)',N,n1,n2,
c     >                             (big(l),l=0,lomax)
                           end do
                           tmp = 0.0 ! find the largest config. for given core orbital n1=ncore(nic)
                           do l=0,lomax
                              if(tmp.lt.big(l)) then
                                 tmp = big(l)
                                 l_tmp = l
c                                 print*,'N,n1,n2,l,tmp:',N,n1,n2,l,tmp
                              end if
                           end do
                           c_big(nic) = tmp
                           l_big(nic) = l_tmp
c                           print*, '&&', l_tmp, tmp
                        end if
                        
                     end do
c     We found that for ionic core orbital n1=ncore(nic)
c     the largest CI contribution c_big(nic)  comes 
c     from configurations with orbital angular momentum l_big(nic)
c     Find which core orbital has largest CI contribution
                     f_big = 0d0
                     m1 = 0
                     do nic=1,nicm
                        if(f_big.lt.c_big(nic)) then
                           f_big = c_big(nic)
                           m1 = nic
                        end if                        
                     end do

                     l1orb(N) = lo(ncore(m1))
                     l2orb(N) = l_big(m1)
                     lic = l_big(m1)
                     num(lic,m1) = num(lic,m1) + 1
                     
                     n1orb(N) =  ko(ncore(m1)) + nabot(l1orb(N)) - 1
                     n2orb(N) =  num(lic,m1) + nabot(l2orb(N)) - 1
c                     print*,'num:', N, num(lic,m1)

                     if(l1orb(N).eq.l2orb(N)) then   ! it is a core orbital (2p)
                        if(ila.ne.ll_gs .or. isa.ne.s_gs .or.    ! it is not singlet S-symmmetry
     >                       ipa.ne.ipar_gs) then
                           n2orb(N) = n2orb(N) + 1
                        endif
                     end if

c     We swap here the orbitals in the major configuration: here the first orbital
c     is ion-core orbital and its number n1 is smaller than n2 in the second orbital.
c     In all other code it is other way around.
                     n_tmp = n1orb(N)
                     n1orb(N) = n2orb(N) 
                     n2orb(N) = n_tmp
                     n_tmp = l1orb(N)
                     l1orb(N) = l2orb(N)
                     l2orb(N) = n_tmp   
c                     
                     
                     write(10,'("major config.:",5I5,F10.5)') N,
     >                    l1orb(N), n1orb(N), l2orb(N), n2orb(N), f_big
                     config_maj(N) = f_big
                     if(f_big .eq. 0.0) then
                        print*,'Attention: for state N=',N,', the ',
     >                     'major configuration is zero, wrong ',
     >                     'major config. labeling.'
                     end if
                  end if
               end do
               
            end do
         end do
      end do
      
      return
      end
c-------------------------------------------------------------------
c---------------------------------------------------------------------------
c     This subroutine sorts target states in the enrgy order.
      subroutine energyorder(Nmax, KNM, nspmCI, la, sa, lpar, np, 
     >   na, nam, nicm, E, chan,Ns,lomax,komax)
      use CI_MODULE
      IMPLICIT NONE
      integer Nmax, KNM, nspmCI
      integer la(KNM), sa(KNM), lpar(KNM), np(KNM)
      integer na(nspmCI,KNM), nam(KNM), nicm
      double precision E(KNM)
      character chan(knm)*3
      integer Ns(0:lomax,0:1,komax,2)
      integer lomax, komax

      integer Num(Nmax)
      integer i, j, itmp, nm
      integer natmp(nspmCI,KNM)
      integer il, is, ip, ik, Nst
      
      nm = MAXVAL(nam(1:Nmax))

c     This array will be sorted by insertion sort algorithm acording to the E(i) array
      do i=1,Nmax
         Num(i) = i
      end do
      do i=2,Nmax
         itmp = Num(i)
         j = i
         do while(j.ge.2) ! .and.E(itmp).lt.E(Num(j-1)))
            if( E(itmp).ge.E(Num(j-1)) ) exit 
            Num(j) = Num(j-1)
            j = j - 1
         end do
         Num(j) = itmp
      end do

      
      do i=1,Nmax
c         write(*,'(2I5,2X,A5,F15.5)') i, Num(i), chan(Num(i)), E(Num(i)) ! 27.2116 
      end do   


      la(1:Nmax) = la(Num(1:Nmax))
      sa(1:Nmax) = sa(Num(1:Nmax))
      lpar(1:Nmax) = lpar(Num(1:Nmax))
      np(1:Nmax) = np(Num(1:Nmax))
      chan(1:Nmax) = chan(Num(1:Nmax))
      nam(1:Nmax) = nam(Num(1:Nmax))
      E(1:Nmax) = E(Num(1:Nmax))
      

      natmp(:,:) = na(:,:)
      do i=1,Nmax
         na(:,i) = natmp(:,Num(i))                
      end do   
      C(1:Nmax,1:nicm,1:nm) = C(Num(1:Nmax),1:nicm,1:nm)

      do il=0,lomax
         do is=0,1
            do ip =1,2
               do ik=1,komax
                  Nst = Ns(il,is,ik,ip)
                  if(Nst .gt. 0) then
                     do i=1,Nmax
                        if(Nst .eq. Num(i)) exit
                     enddo
                     Ns(il,is,ik,ip) = i
c                     print*,'!',il, is, ip, ik, Nst, i
                  endif
               enddo
            enddo
         enddo
      enddo
         

      return
      end
c--------------------------------------------------------------------
c     < N || w*z_1 +z_2|| Np >, |N> = |l^w, n_2 l_2 : LS>, w = 4l+1
      subroutine   osc_ng(Nmax,N,Np,l,lp,fl,minf,maxf,lo,
     >     dip,deriv,dippol)
      use CI_MODULE             ! module with CI-array C
      include 'par.f'
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), npdrop(KNM)
      dimension fl(nmaxr,nspmax), lo(nspmax)
      dimension  minf(nspmax), maxf(nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      double precision  ortint,trm,trm1,r1elk,f1deriv,
     >   dip,deriv,sum,coef,tmp,ang,polsum,trmpol,dippol,dintegral
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer na, nam
      common /corearray/ nicm, ncore(nspmCI)
      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      common /di_el_core_polarization/ gamma, r0, pol(nmaxr)
      common /meshrr/ nr,gridr(nmaxr,3)
      double precision ttt
c
      dip = 0.0D0
      deriv = 0.0D0
      dippol = 0.0D0

      if(ls(N) .ne. ls(Np)) return

      nfcel = 4*l_ng + 1       ! number of 'frozen'-core electrons (5 : for Ne 2p^5)

      trm = 0.0D0
      trm1 = 0.0D0
      trmpol = 0.0D0
      rl = l
      rlp = lp
      ttt = dsqrt((2.0*l+1d0)*(2.0*lp+1d0))


      coef = 0.0      
      do nic=1,nicm             ! go through all core orbitals
         n1 = ncore(nic)
         jn1 = nic
         l1 = lo(n1)
         rl1 = l1

         do jn2=1,nam(N)
            n2 = na(jn2,N)
            l2 = lo(n2)
            rl2 = l2
            
            do nicp=1,nicm       ! go through all core orbitals
               n1p = ncore(nicp)
               jn1p = nicp
               l1p = lo(n1p)
               rl1p = l1p
               
               do jn2p=1,nam(Np)
                  n2p = na(jn2p,Np)
                  l2p = lo(n2p)
                  rl2p = l2p                  

                  coefsym = 1.0
                  if(n2 .eq. n1 .and. n2p .eq. n1p) then
                     coefsym = 1.0
                  elseif(n2 .ne. n1 .and. n2p .ne. n1p) then
                     coefsym = 1.0
                  else
                     coefsym = sqrt(nfcel+1d0)
                  endif
                  
                  ciN = C(N,jn1,jn2)
                  ciNp = C(Np,jn1p,jn2p)
                  ciprod = ciN * ciNp

c     deal with <l^w|z_1|l^w> - it should be zero by parity consideration, but keep it just to have a general case for a number of cores in future... ... needs to be checked
                  if( (abs(l1-1).le.l1p.or.l1p.le.l1+1 )
     >                 .and. ciprod .ne. 0) then
                     tmp = coefsym*ortint(n2,n2p)*ciprod  
                     sum1 = 0.0
                     if(tmp.ne.0d0) then
                        ang = -dble(nfcel)*ttt*dsqrt(2.0*l1p+1d0)
     >                       *dble(CGC0(rl1p,1.0,rl1))
     >                       *dble(COF6J(rl2,rl1p,rlp,1.,rl,rl1))
     >                       *dble((-1)**(l1+l2+lp+1))
                        sum1 = ang * tmp
                     end if
                     if(sum1.ne.0d0) then
                        coef = - l1p
                        if(l1p.eq.l1+1) coef = l1p+1
                        dintegral = sum1*r1elk(1,fl(1,n1),fl(1,n1p),
     >                       minf(n1),minf(n1p),maxf(n1),maxf(n1p))
                        trm = trm + dintegral
                        dderiv = sum1*( f1deriv(fl(1,n1),fl(1,n1p),
     >                       minf(n1),minf(n1p),maxf(n1),maxf(n1p)) 
     >                       + (coef-1d0)*r1elk(-1,fl(1,n1),fl(1,n1p),
     >                       minf(n1),minf(n1p),maxf(n1),maxf(n1p)) )                     
                        trm1 = trm1 + dderiv 
                        polsum = 0d0
                        if(gamma.ne.0.0) then
                           mini = max(minf(n1),minf(n1p))
                           maxi = min(maxf(n1),maxf(n1p))
                           do i=mini,maxi
                              polsum = polsum + gridr(i,3)*
     >                             pol(i)*fl(i,n1)*fl(i,n1p)
                           end do
                           polsum = polsum*gamma
                        end if
                        trmpol = trmpol + dintegral - sum1*polsum
                     end if               
                  end if
c                  print'(6i5,2e12.4)',N, Np, n1,n2,n1p,n2p,trm,trm1
c
c     deal with <l2|z_2|l2p>
                  if( (abs(l2-1).le.l2p.or.l2p.le.l2+1)
     >                 .and. ciprod .ne. 0) then
                     tmp = coefsym*ortint(n1,n1p)*ciprod 
                     sum2 = 0.0
                     if(tmp.ne.0d0) then
                        ang = ttt*dsqrt(2.0*l2p+1d0)
     >                       *dble(CGC0(rl2p,1.0,rl2))
     >                       *dble(COF6J(rl1,rl,rl2,1.,rl2p,rlp))
     >                       *dble((-1)**(l1+l+1+l2p))
                        sum2 = ang * tmp
                     end if
                     if(sum2.ne.0d0) then
                        coef = - l2p
                        if(l2p.eq.l2+1) coef = l2p+1
                        dintegral = sum2*r1elk(1,fl(1,n2),fl(1,n2p),
     >                       minf(n2),minf(n2p),maxf(n2),maxf(n2p))
                        trm = trm + dintegral
c     print*,n1,n2,n1p,n2p,dintegral,trm
c     print*,'norm=',n2p,r1elk(0,fl(1,n2p),fl(1,n2p),
c     >                       minf(n2p),minf(n2p),maxf(n2p),maxf(n2p))
                        dderiv = sum2*( f1deriv(fl(1,n2),fl(1,n2p),
     >                       minf(n2),minf(n2p),maxf(n2),maxf(n2p)) 
     >                       + (coef-1d0)*r1elk(-1,fl(1,n2),fl(1,n2p),
     >                       minf(n2),minf(n2p),maxf(n2),maxf(n2p)) )
                        trm1 = trm1 + dderiv
                        polsum = 0d0
                        if(gamma.ne.0.0) then
                           mini = max(minf(n2),minf(n2p))
                           maxi = min(maxf(n2),maxf(n2p))
                           do i=mini,maxi
                              polsum = polsum + gridr(i,3)*
     >                             pol(i)*fl(i,n2)*fl(i,n2p)
                           end do
                           polsum = polsum*gamma
                        end if
                        trmpol = trmpol + dintegral - sum2*polsum
                     end if
                  end if
c                  print'(6i5,2e12.4)',N, Np, n1,n2,n1p,n2p,trm,trm1
c

               
            end do
         end do
         end do
      end do
      dip = trm
      deriv = trm1
      dippol = trmpol
      return
      end



C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Matrix elements here are non-reduced (<|V|>) The difference is: <|V|> = sqrt(2J+1) <||V||>
      subroutine vdme_ng(nze,nchmi,Nmaxi,Ni,Li,lai,si,lpari,nspmi,loi,
     >   Ci,fli,maxfi,minfi,npki,chili,minci,
     >   nchmf,Nmaxf,Nf,Lf,laf,sf,lparf,nspmf,lof,
     >   Cf,flf,maxff,minff,npkf,chilf,mincf,ortint,KJ,dwpot,
     >   vdon,VVD,nchf,nchi,na,nam,namax)
      include 'par.f'
      double precision Z
      common /Zatom/ Z
      integer lai(Nmaxi), si(Nmaxi), lpari(Nmaxi)
      integer laf(Nmaxf), sf(Nmaxf), lparf(Nmaxf)
      real*8 Ci(Nmaxi,namax,namax),Cf(Nmaxf,namax,namax),
     >   ortint(nspmax,nspmax)
      dimension loi(nspmi), lof(nspmf),npki(nchmi+1),npkf(nchmf+1)
      dimension fli(maxr,nspmax),maxfi(nspmi),minfi(nspmi)
      dimension flf(maxr,nspmax),maxff(nspmf),minff(nspmf)
      dimension chili(nr,npki(nchmi+1)-1), minci(npki(nchmi+1)-1)
      dimension chilf(nr,npkf(nchmf+1)-1), mincf(npkf(nchmf+1)-1)
      dimension formout(maxr),dwpot(maxr,nchan)!,chiform(maxr,kmax)
      allocatable :: chiform(:,:)
      dimension VVD(kmax,kmax,0:1)
c$$$      dimension VVD(npkf(nchmf+1)-1,npki(nchmi+1))
      dimension  Vang(0:ltmax), upot(maxr), vdon(nchan,nchan,0:1)
      common /meshrr/ nr, gridr(nmaxr,3)
      dimension temp(maxr), zeroUpot(maxr)
c      logical zeroc
      logical diag_ME
      real tmp_const(0:ltmax), tmp_const2(0:ltmax)
      dimension na(nspmCI,KNM), nam(KNM)
      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      real*8:: tmp_Cf, tmp_Ci
      
c     Note : antisymmetrization is NOT included through CI coefficients.
c     explicit formulas per antysymmetric configuration
c     
      hat(l)=sqrt(2.0 * l + 1.0)
      ze = 2.0    !!!!!!!!!!!  to be corrected

      pnorm = 2.0/acos(-1.0)
      temp(:) = 0.0
      zeroUpot(:) = 0.0

      diag_ME = Ni.eq.Nf   ! set true if this is diagonal ME, to be used for subtraction of the nuclear term


      nfcel = 4*l_ng + 1        ! number of 'frozen'-core electrons (5 : for Ne 2p^5)
      
c     
c     These are 2 channels loops (nch = (L,la(N),s(N))
c     where L is projectile ang.mom., N is number of physical channel,
c     la(N),s(N)  are  target ang. and spin mom..
      nqmi = npki(nchi+1) - npki(nchi)
      do i = 1, nr
         upot(i) = dwpot(i,nchi) 
      enddo 
      nqmf = npkf(nchf+1) - npkf(nchf)

      formout(:) = 0.0

      minform = nr
      maxform = 0
      
      if(si(Ni).ne.sf(Nf)) return


      rK=KJ
      lli=lai(Ni)
      rLi=Li
      rlli=lli
      llf=laf(Nf)
      rLf=Lf
      rllf=llf
c     Ni,Nf  -   numbers of physical channels
c     s(Ni)=s(Nf) for direct VME

c     lambda loop: to avoid recalculations of 6J and CGC that do not depend
c     on ang. momenta of s.p.functions
      tmp_hat =  hat(lli)*hat(llf)*hat(Li)
      lamtop1=min(Li+Lf,lli+llf)
      lambot1=max(iabs(Li-Lf),iabs(lli-llf))
      do lam=lambot1,lamtop1
         rlam=lam
         tmp_const2(lam) =  tmp_hat *
     >        CGC(rLi,0.,rlam,0.,rLf,0.)*
     >        COF6J(rLf,rllf,rK,rlli,rLi,rlam)
         tmp_const(lam) = tmp_const2(lam)*(-1)**(llf+Li+KJ+lli+lam)
      end do



!!$ first deal with <l^w | V | l^w> term


c     these are 2 s.p. loops on s.p. functions of the coordinate r1
      do 68 jnf1=1,nam(Nf)
         nf1 = na(jnf1,Nf)
         lf1=lof(nf1)
         rlf1=lf1
         tmp2 = (-1)**(lf1)
         do 67 jni1=1,nam(Ni)
            ni1 = na(jni1,Ni)
            li1=loi(ni1)
            rli1=li1
            tmp1 = tmp2 * hat(li1) 
c     lambda loop:
            Vang(:) = 0.0
            lamtop=min(li1+lf1,lamtop1)
            lambot=max(iabs(li1-lf1),lambot1)
            
c     lambda loop for calculation of the angular coefficent:
            do 5 lam=lambot,lamtop
               rlam=lam
               
               const = tmp_const(lam) * tmp1 *
     >              CGC(rli1,0.,rlam,0.,rlf1,0.)
               if (const.eq.0.0) go to 5
               
c     These are loops on s.p. functions of the coordinate r2,
c     coordinate r2 is not affected by V(r0,r1), but  nf2 is not necessary 
c     equal to ni2 because nonorthogonal s.p. basis is used, but li2=lf2.
               
               lf2_old = -1
               do 10 jnf2=1,nam(Nf)
                  nf2 = na(jnf2,Nf)
                  lf2=lof(nf2)
                  rlf2=lf2
                  tmp_Cf = Cf(Nf,jnf1,jnf2)
                  do 20 jni2=1,nam(Ni)
                     ni2 = na(jni2,Ni)
                     tmp_Ci = Ci(Ni,jni1,jni2)
                     if(tmp_Ci.eq.0.0d0) goto 20
                     rli2 = loi(ni2)
                     const2 = ortint(ni2,nf2)*tmp_Cf*tmp_Ci
                     if (const2.eq.0.0) goto 20
                     
                     if(nf1 .eq. nf2 .and. ni1 .eq. ni2) then
                        const2 = const2 * (nfcel + 1) ! lam = 0 term only can be here, deal with all terms here - see if() statement for coord. 2 
                     elseif(nf1 .ne. nf2 .and. ni1 .ne. ni2) then
                        if(lam .eq. 0 ) then
                           const2 = const2 * nfcel
                        else
                           const2 = const2 * (-1)**(lam+1)                              
                        endif
                     else
                        const2 = 0d0 !  by orthogonality of the core orbital with any other orbital ortint(ni2,nf2)=0, should not really get to here...
                        print*,'something is wrong...'
                     endif
                     
                     if(lf2.ne.lf2_old) then
                        tmp_6J = (-1)**lf2 *
     >                       COF6J(rlf2,rlf1,rllf,rlam,rlli,rli1) 
                        lf2_old = lf2
                     end if
                     Vang(lam) = Vang(lam) + const2 * tmp_6J
                     
c     end ni2,nf2 loops
 20               continue
 10            continue
               
               Vang(lam) = Vang(lam) * const
c     end lambda loop                    
 5          continue 
            do lam=lambot,lamtop
               if (Vang(lam).ne.0.0) then
                  
                                !                        print*,"###",nf,ni,ni1,nf1,Vang(lam)
                  call getformout(lam,Vang(lam),zeroUpot,
     >                 fli(1,ni1),flf(1,nf1),
     >                 minfi(ni1),minff(nf1),
     >                 maxfi(ni1),maxff(nf1),nr,diag_ME,
     >                 formout,minform,maxform)
                  
               end if
            end do
                                            
c!!$ deal with < l2 | V | l2p > term. 
            if(ortint(ni1,nf1).eq.0d0) cycle
            do jnf2=1,nam(Nf)
               nf2 = na(jnf2,Nf)
               lf2=lof(nf2)
               rlf2=lf2
               tmp_Cf = Cf(Nf,jnf1,jnf2)
               if(Cf(Nf,jnf1,jnf2).eq.0.0d0) cycle
               
               do jni2=1,nam(Ni)
                  ni2 = na(jni2,Ni)
                  if(Ci(Ni,jni1,jni2).eq.0.0d0) cycle
                  li2 = loi(ni2)
                  rli2 = li2
                  const2 = ortint(ni1,nf1) * tmp_Cf * 
     >                 Ci(Ni,jni1,jni2)
                  if (const2.eq.0.0) cycle
                  
                  if(nf1 .eq. nf2 .and. ni1 .eq. ni2) cycle ! this case is accounted for above:   const2 = const2 * (nfcel + 1)
                  
                  if(nf1 .eq. nf2 .or. ni1 .eq. ni2) then
                     coefsym = sqrt(nfcel+1d0)
                  else
                     coefsym = 1d0
                  endif
                  
c     lambda loop:
                  lamtop=min(li2+lf2,lamtop1)
                  lambot=max(iabs(li2-lf2),lambot1)
                  
                  do lam=lambot,lamtop
                     rlam = lam
                     
                     tmp1 =  (-1)**(Li+KJ+li1+li2+lam)*hat(li2)
                     const = tmp_const2(lam) * const2 * tmp1 *
     >                    CGC(rli2,0.,rlam,0.,rlf2,0.)*coefsym
                     if (const.eq.0.0) cycle
                     tmp_6J=COF6J(rlf1,rllf,rlf2,rlam,rli2,rlli) 
                     
                     angCoef = const * tmp_6J
                                !                           print*,'###',nf,ni,ni1,nf1,nf2,ni2,
                                !     >                          "angCoef=",angCoef, tmp_6J, const
                     if (angCoef.ne.0.0) then
                        call getformout(lam,angCoef,zeroUpot,
     >                       fli(1,ni2),flf(1,nf2),
     >                       minfi(ni2),minff(nf2),
     >                       maxfi(ni2),maxff(nf2),nr,diag_ME,
     >                       formout,minform,maxform)
                     end if
                  end do        ! end lam loop
                  
               enddo            ! end ni2 loop
            enddo               ! end nf2 loop
            
c     end ni1,nf1 loops
 67      continue
 68   continue                                     
      
      
      if (nchi.eq.nchf.and.nqmi.eq.1.and.dwpot(1,nchi).eq.0.0)
     >     then
C     Define the channel dependent distorting potential when running the
C     Born case. The units are Rydbergs.
         do i = minform, maxform
            dwpot(i,nchi) =  formout(i) * 2.0
         enddo
                                !               print*,'Defined DPOT(i,nch):',nchi,dwpot(1,nchi)
      endif
      
c     Account for dist. potential for diagonal channel only
      if(nchf .eq. nchi .and. Ni. eq. Nf) then
         formout(minform:maxform) = formout(minform:maxform) 
     >        - upot(minform:maxform)/2.0
      endif

c$$$      if(nchf.eq.1 .and. nchi.eq.1 ) then
c$$$         read(*,*)
c$$$         kistart = npki(nchi)
c$$$         open(174,file='tmp_Ne_formout')
c$$$         write(174,'("#KJ=",I5)') KJ
c$$$         write(174,'("#minf,maxf=",2I5)') minfi(1),maxfi(1)
c$$$         do i=minform,maxform
c$$$            write(174,'(I5,F15.5,1P,10E15.5)') i, gridr(i,1),
c$$$     >           formout(i),chili(i,1), fli(i,1)
c$$$         enddo
c$$$         close(174)
c$$$      endif
                                !            print*,"!!",formout(100)
      
c     momentum loops
      kistep = npki(nchi) - 1
      const = - nze * pnorm
      do i = minform, maxform
         temp(i) = formout(i) / gridr(i,3)
      enddo
      
      allocate(chiform(maxr,kmax))
      do ki=npki(nchi), npki(nchi+1) - 1
         minki = max(minform, minci(ki))
         do i = minki, maxform
            chiform(i,ki-kistep) = temp(i) * chili(i,ki)
         enddo
      enddo
       kistart = npki(nchi)
      kistop = npki(nchi+1) - 1
      quart = 0.0
      if (si(Ni).eq.1.0) quart = 1.0
C     Note si(Ni) = sf(Nf) here
c$$$  c$par doall
      tmp = 0.0
      do ki=kistart, kistop
         minki = max(minform, minci(ki))
         kii = ki - kistep
         do 90 kf = npkf(nchf), npkf(nchf+1) - 1
            kff = kf - npkf(nchf) + 1
            mini = max(minki, mincf(kf))
            n = maxform - mini + 1
            if (kf.lt.ki.or.n.le.0) go to 90
                                !                  if(ki.eq.kistart .and. kf.eq.npkf(nchf)) then
                                !                     print*,"!!>",chilf(100,kf),chiform(100,kii)
                                !                  endif
            tmp = dot_product(chilf(mini:maxform,kf),
     >           chiform(mini:maxform,kii))
            
                                !                  if(nchi.eq.nchf) then
                                !                     print*, nchf,nchi,tmp
                                !                  endif
c$$$            if(nchf.eq.1 .and. nchi.eq.1 .and.nqmi.eq.1) then
c$$$               print*,'tmp=',tmp, const, VVD(kf,ki)
c$$$            endif
c$$$  tmp = sdot(n,chilf(mini,kf),1,chiform(mini,kii),1)
            vvd(kff,kii,0) = vvd(kff,kii,0) + tmp * const
            vvd(kff,kii,1) = (vvd(kff,kii,1) + tmp * const)*quart

c$$$            VVD(kf,ki) = VVD(kf,ki) + tmp * const
C     QUART is used to make sure that the quartet V matrix elements are
C     nonzero for triplet-triplet transitions only. The multiplication
C     by QUART zeros the contribution from the core.
c$$$            VVD(ki,kf+1) = (VVD(ki,kf+1) + tmp * const) * quart
c     end ki,kf loops
 90      continue 
      end do

      deallocate(chiform)

      vdon(nchf,nchi,0) = vdon(nchf,nchi,0) + tmp * const
      vdon(nchi,nchf,0) = vdon(nchf,nchi,0)
      vdon(nchf,nchi,1) = (vdon(nchf,nchi,1)+tmp*const) * quart
      vdon(nchi,nchf,1) = vdon(nchf,nchi,1)


      return
      end
c*********************************************************************
c     NOT USED
c     This routine is called only from the direct ME routine vdme(...)
      subroutine getformout_ng(lam,fli,flf,minfi,minff,
     >    maxfi,maxff,maxni1,temp,i1,i2)
       include 'par.f'
       common /meshrr/ nr, gridr(maxr,3)
       dimension  fli(maxr), flf(maxr)
       dimension temp(maxr), fun(maxr)
       common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >    minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
C the following is not used in the nuclear routine
       common /di_el_core_polarization/ gamma, r0, pol(nmaxr)
       
       minfun=max(minfi,minff)
       maxfun=min(maxfi,maxff)
       do i=minfun,maxfun
          fun(i) = fli(i) * flf(i) * gridr(i,3)
       end do

       call form(fun,minfun,maxfun,rpow1(1,lam),rpow2(1,lam),
     >    minrp(lam),maxrp(lam),maxni1,temp,i1,i2)
c     if this is diagonal ME then for lambda = 0  do subtraction of the nuclear term

!       if (lam.eq.0.and.torf) call nuclear(fun,.true.,minfun,maxfun,i2,
!     >      upot,nznuc,temp)

       if(lam.eq.1.and.gamma.ne.0.0) then
          sum1 = 0.0
          do i=minfun,maxfun
             sum1 = sum1 + fun(i)*pol(i)
          end do                       
          tmp = gamma * sum1 
          do i=i1,i2             
             temp(i) = temp(i) - tmp*pol(i)
          end do
       end if
       

       return
       end
c
c*********************************************************************
c
      subroutine makevexch_DF(lactop, lamax, nabot, lnabmax, nnmax, 
     $     ltmax, maxr, istoppsinb, psinb, rpow2, vdcore, nznuc, E, 
     $     maxvexch, vexch)
C=======================================================================
C     calculates the local-exchange core potential of the form suggested
C     by Furness & McCarthy (1973). This potential is used instead of
C     the exchage part of the Hartree-Fock potential.
C           
C INPUT:  E      - optimization parameter;      
C         vdcore - direct core (Hartree) potential;
C      
C     
C OUTPUT: maxvexch  - Vexch(i)=0 if i>maxvexch  
C         vexch     - LEA potential
C     
C     USE BEFORE THE POLARIZATION POTENTIAL IS ADDED TO VDCORE
C
C========================================================================
C
C     note that the electron-ion potential is Vst(r) = Vdcore(r) - 1/r
C
C======================================================================== 

      dimension nabot(0:lamax), istoppsinb(nnmax,0:lnabmax)
      dimension psinb(maxr,nnmax,0:lnabmax), rpow2(maxr,0:ltmax)
      dimension vdcore(1:maxr)
      dimension density(1:maxr)
      real vexch(1:maxr)
      integer maxvexch
            
C     LOCAL VARIABLES:
      real*8 density, const, Rnl, ve, ve2, rho

      vexch(:) = 0.0 

C     coef = 1.0
      if(nznuc .eq. 10) then
         coef = 0.6
      else
         stop 'Not defined for non Ne targets'
      endif
!      print*, 'Enter coef'
!      read(*,*) coef

C     calculate density
      maxvexch=0
      density(:)=0.0
      do lac = 0, lactop
         const = float(4 * lac + 2)
         do nac = lac + 1, nabot(lac) - 1
            maxvexch=max(maxvexch,istoppsinb(nac,lac))
            do i = 1, istoppsinb(nac,lac)
               Rnl=psinb(i,nac,lac)*rpow2(i,0) 
               density(i)=density(i) + const*Rnl*Rnl
            end do            
         end do
      end do
      
C     Slater exchange 
      
      pi = acos(-1d0)
      tmp = -(3d0/2d0)*(3d0/(4d0*pi*pi))**(1d0/3d0)  /2.0
c      tmp = -(3d0/2d0)*(24d0/pi)**(1d0/3d0) / 2d0
c      tmp = -(24d0/pi)**(1d0/3d0) / 2d0
      
      tmp = tmp * coef   ! adjustable coef

      do i = 1, maxvexch  
          vexch(i)  = tmp*(density(i))**(1d0/3d0)
       end do 

      return
      end subroutine makevexch_DF


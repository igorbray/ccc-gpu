C  File main.f
C  This program solves a set of coupled integral equations
C  written by Igor Bray for one-electron targets and extended to two-electron
C  targets by Dmitry Fursa.
c$$$      module psinc_module
c$$$      real, allocatable :: psinc(:,:,:)
c$$$      end module psinc_module
      program CCC      
      use vmat_module
c$$$      use chil_module      
      use ubb_module ! Andrey
      use ql_module  ! Andrey
      use apar                  ! Andrey: maxx, nmax, lmax1, maxq
      use ftps_module
      use gf_module
      use chil_module
c MPI Fortran header file
      include 'mpif.h'
      include 'par.f'
c MPI declarations
      integer ntasks, myid, ierr
      integer nlg, my_nlg, mylstart, mylstop
c$$$      pointer (ptrvmat,vmat)
      pointer (ptrvkeep,vkeep)
      pointer (ptre2ed,ve2ed)
      pointer (ptre2ee,ve2ee)
c$$$      real vmat(kmax,nchan,kmax,nchan+1),vdon(nchan,nchan,0:1)
      real, allocatable :: vmat(:,:), vmatp(:,:), AA(:,:), bb(:,:,:,:)
      integer, allocatable :: ipvt(:)
      real vkeep(1,1),vdon(nchan,nchan,0:1),
     >   vdondum(nchan,nchan,0:1), ve2ed(1,nchane2e,kmax,nchan),
     >   ve2ee(1,nchane2e,kmax,nchan)
      complex wk(kmax*nchan)
      complex*16 coulphase
      complex vmatop(kmax,kmax,0:nchanop,nchanop)
      pointer (ptrvopt,vmatop)
      include 'par.for'
      real*8 ffgg, factl, api, x2, w2, res, gamx, rylm, cknd, rlambda,
     >   Z, deta
      common /Zatom/ Z
      common /cnsts2/ factl(0:2*lcoul),api,x2(63,5),w2(63,5),res(63)
      complex phasen(ncmax,0:lamax),phaseq(knm)
      common/meshrr/ meshr,rmesh(maxr,3)
      COMMON/CONT/ CE(ncmax),CINT(ncmax),NCSTATES,ENERGY
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax),itail
      integer nptopin(0:lamax),lnch(nchan,2),nchtime(nchan)
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      dimension minps(ncmax),npk(nchan+1),npkb(nchan+1),npkeep(nchan+1),
     >   corep(0:lamax), r0(0:lamax), vdcore(maxr,0:lamax),vexch(maxr),
     >   vdcore_st(maxr,0:lamax), vdcore_pr(maxr,0:lamax)
      COMMON /gstspin/   iSpin, meta
      common /ininfo/isecond
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)
      common/smallr/ formcut,regcut,expcut,fast,match,analyticd,boxps
      real erange(kmax)
      logical fast, match, speed, analyticd, boxps, pos(nchan)
      integer satom, OMP_GET_MAX_THREADS, natomps(nchan)
      integer ,allocatable :: natompsnode(:)
      integer pnewC
      common /dynamical_C/ nmaxhe, namax,pnewC
      common /helium/ latom(KNM), satom(KNM), lpar(KNM), np(KNM)
      character opcl*10, e2efile(ncmax)*60,target*6,projectile*8,
     >   tfile*80,csfile*80,ench*11,ch*1,nodetfile*9
      character date*8,time*10,zone*5
      integer*8 npernode
      integer valuesin(8), valuesout(8), inc(1000,0:1), lgold(0:1),
     >     valuesinLG(8),incold(1000,0:1),nodesold(0:1)
      integer lmatch(0:lamax), nopen(0:1), instate(1000)
      common /radpot/ ucentr(maxr)
      common /double/id,jdouble(22)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      dimension chi(maxr), nchns(0:1,nchan,nchan),
     >   ovlp(ncmax,0:lamax), ovlpn(knm), ovlpnn(knm,nnmax),
     >   ovlpnl(ncmax,0:lamax,nnmax),
     >   slowe(ncmax), slowery(ncmax), BornICS(knm,knm)
      complex ovlpnl, ovlpnn
      real, allocatable :: kon(:,:) !,chil(:,:,:)
      complex, allocatable :: ct(:,:), ckern(:,:), cwork(:), ct2(:),
     $     psi12(:,:)
      integer, allocatable :: sendHandle(:), statusesR(:,:),
     >   statusesS(:,:), receiveHandle(:,:), receiveType(:,:)
      integer ,allocatable :: nrange(:) !,minchil(:,:)
      real ,allocatable :: tscale(:)
c$$$      pointer (ptrchi,chil)
c$$$      dimension chikeep(maxr,kmax,nchan),minchikeep(kmax,nchan)
      dimension chikeep(1),minchikeep(kmax,nchan),tave(0:1)
      pointer (ptrchikeep,chikeep)
      dimension nk(5,0:lmax,2), sk(5,0:lmax,2),
     >   npsp(0:lamax), alphap(0:lamax), nbnd0(0:lmax),
     >   alpha(0:lamax), nps(0:lamax), nbnd(0:lmax),abnd(0:lmax)
      complex phase,phasel(kmax,nchan),sigc,sigma(nchan),
     >   dphasee2e(nchane2e),ephasee2e(nchane2e),ctemp(maxr,30)
      dimension ui(maxr), uplane(maxr),qgrid(kmax),weightk(kmax)
      dimension gk(kmax,nchan),u(maxr),vnucl(maxr),gkeep(kmax,nchan),
     >   noprint(nchan),vasymp(maxr),temp(maxr), dwpot(maxr,nchan)
      character chan(knm)*3,cnode*3,chinstate(knm)*3,chann(knm)*5
      common /charchan/ chan
      common /chanen/ enchan(knm)
      common /formwarn/iw
      common /underflows/ nunder
      common /overflows/ nover
      common /MPI_info/myid, ntasks, cnode, ench
      logical small(0:1,nchan,nchan),noprint,hlike,uba(nchan),
     >   canstop(0:1), second, bothpar, allocated, exists,myflag
      ch(i)=char(i+ichar('0'))
      nchprs(n1,n2,n3) = (n2-n1+1)*(n2-n1+2)/2+(n2-n1+1)*(n3-n2)
      data uba/nchan*.false./
      common /relcorsetup/ iih
      real vrel(4000), rrel(0:lomax)
      common /atompol/ atompol(maxr), ipolmin, ipolmax, ipol
      common /di_el_core_polarization/ gamma, r0_di_el, pol(maxr)
      common /data_for_semirelBorn/ energy_for_semirelBorn, nunit
      common/nofcexch/ inofcexch, nfcmax
      common /TMP_AB/ A,B
      common /theta_scat/ theta
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /noblgegas_switch/ i_sw_ng
      common /statenumber_str_to_scat/ Nstate_str_to_scat, na_ng(KNM)
      common /debye/dbyexists, dblp, rmudeb, zdby
      logical::dbyexists,scaleexists, timeexists

c DPI of Lithium       
      common /LITHIUM/ lithium

      include 'par.pos'
      parameter (maxl=2*ltmax+1)

C     real *8 :: xgz1, wgz1, pl
C     common/gausz/xgz1(igzm),wgz1(igzm),pl(0:maxl,igzm),igz

      common/numericalpotential/ numericalv, lstoppos
      logical numericalv
      
C     ANDREY__________________________________________________
      real*8, dimension (1:2*(2**nmax-1)) ::  xgz, wgz
      real*8, dimension (0:lmax1, 1:2*(2**nmax-1)) :: pleg     
      common /gausc/ xgz, wgz, pleg
      
      real dx
      real, dimension (maxx) :: rmesh1      
      real, dimension (maxx) :: xmesh, ve, va, vaM
      common /potentials/ xmesh, ve, va, vaM, rmesh1, dx

C     fourier transforms of ps and ps*v (see makeftps)
      ! parameter maxp = 20000
      ! real, dimension (0:maxp) :: pmesh
      ! real*8, dimension (0:maxp,nnmax,0:lnabmax) :: ftps, ftpsv
      ! common /makeftpsc/ ftps, ftpsv, pmesh         
       
      real res0, res1
C      real veq(maxq), qmesh(maxq)
      integer, parameter :: dp=selected_real_kind(precision(1.0d0))
      real(kind=dp), dimension (nr) :: qmesh, veq
      common /veqc2/ qmesh, veq
      
      real, external :: electroncoreV

      real, dimension(3) :: c ! used in five
      common /veqc/ c      

      common /fivec/ veold(maxr),vaold(maxr)
      real ctail
      real*8 vefun, a1, b1
      common /funleg2c/ vefun(maxr,2), corepp, r0p, a1, b1, maxvdc
            
      integer icase(2)
      logical true_potential(2), exchange
                  
C     data for 2-D interpolation:    
      real, dimension(maxr) :: ubbb0, ubbb3            
      external positron
      logical positron, torf, reconstruct_psi

C     ql_module variables:
      real*8 pi8, lgqa, lgq, qq1, qa, dQlp

C     to test ftps1:
!      integer, parameter :: dp=selected_real_kind(15), dim=512
!     real(kind=dp), dimension(dim) :: kk, ft1

      character fname*5, fname1*7, enq*16

c     test for funleg2
      real*8 x, dx1, qq(5000), QL(5000)
      real*8 arg, Q0p, QlH
      real*8 qa2, pp, pp2, arg1
 
      dimension vdcore0(maxr,0:lamax)

      real, external :: vaint, veint, vaMint ! ANDREY: DEBUGING
#ifdef _CRAYFTN
      character*24 stacksize
#endif
      !integer, parameter :: dp=32
      integer, parameter :: qp = selected_real_kind(15, 307)
      !integer, parameter :: qp = selected_real_kind(33, 4931)
      !real (kind = qp), external :: tint
      real (kind = qp) :: QLvd, QLpp, part1, part2, part3 

c      real, allocatable :: vmat01(:,:), vmat0(:,:), vmat1(:,:)
c      integer, allocatable :: nchistart(:), nchistop(:)
c      integer  nodeid
c      logical scalapack

      
c$$$      integer :: provided, required
c$$$      required = MPI_THREAD_MULTIPLE
c$$$      call MPI_Init_thread(required,provided,ierr)
c$$$      if (provided .ne. required) then
c$$$         print*,'provided',provided
c$$$      endif
      va(:)=0.0
      ve(:)=0.0
      veold(:)=0.0
      vaold(:)=0.0
      xmesh(:)=0.0
      dx=0.0;
      c(:)=0.0
C  VDCORE will contain the direct plus l-dependent core polarization potential
      vdcore(:,:) = 0.0 !atomic units
      vdcore_st(:,:) = 0.0 ! for electrons, used for calculating the structure
      vdcore_pr(:,:) = 0.0 ! for projectile, electron or positron
      
C==== ANDREY ===================================================
      
      nspm = 0
      i_sw_ng = 0

      lstoppos = -1
      lptop = -1
      inc(:,:) = 0
      incold(:,:) = 0
      lgold(:) = -1
      nodesold(:)=-1
      myid = -1
      ntasks = 1
      nbnd0(:) = 0
#ifdef _CRAYFTN
      call pxfgetenv('OMP_STACKSIZE',0,stacksize,0,ierr)
      print*, 'OMP_STACKSIZE:',stacksize
#else
c      print*, 'KMP_STACKSIZE:',kmp_get_stacksize_s()
#endif
#ifdef _single
#define nbytes 4
#elif defined _double
#define nbytes 8
#endif

c MPI initialization and environment characteristics
c note: make these your first (or nearly so) executable statements
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, ntasks, ierr )
      nomporig = max(1,OMP_GET_MAX_THREADS())
      nomp = 1 ! nomporig ! revert for many tasks per node, see same comment below and in redistribute.f
c$$$      print*,'NUM_PROCESSES_PER_NODE:',NUM_PROCESSES_PER_NODE()
c$$$      if (ntasks.ge.8.and.nomp.ne.8) stop 'ntasks.ge.8.and.nomp.ne.8'
      nodes = max(1, ntasks / nomp)

      allocate(nchistart(nodes),nchistop(nodes))
      allocate(nchistartold(nodes,0:1),nchistopold(nodes,0:1))
      allocate(ntime(nodes,0:1),nntime(nodes))
      allocate(natompsnode(nodes),nrange(nodes),tscale(nodes))
      ntime(:,:)=0
      nntime(:)=0
C define the following here to stop "used without being defined" errors in do 770
      mylstart=0
      mylstop=0
      iparmin=0
      npar=0
!      inquire(file='scalapack',exist=scalapack)
      scalapack = nodes.gt.1 ! delete to use LAPACK, but with several nodes for VMAT
      if (myid.le.0) print*,'scalapack:',scalapack
      CALL CLOCK(START)

      if (mod(myid,nomp).eq.0) then
         nodeid = myid/nomp + 1
      allocate(receiveType(nodes,3))!These calls are for LAPACK with many nodes
      allocate(receiveHandle(nodes,3))
      allocate(statusesR(MPI_STATUS_SIZE,nodes))
      allocate(sendHandle(3))
      allocate(statusesS(MPI_STATUS_SIZE,3))

      
      do n = 1, knm
         latom(n) = 0
         satom(n) = 0
         lpar(n)  = 0
         np(n)    = 0
      enddo 
C  Define logs of factorials and reduced factorials
      call faclog0
      call fakred
      call dfset
      call gamset
      call psiset
C     initialize stuff for tail integrals
      call initialize

C

C  the following units correspond to files used by the program
C  3  ccc.in:  the input file
C  42 totalcs: integrated cross sections
C  88 potl:    T-matrix elements
C

!      call magmaf_init()


      nunder=0
      nover=0
C
C  SET UP SOME INITIAL CONSTANTS
C
      
      RY= 13.6058
      
      PI= ACOS(-1.)
      ANORM= 8./PI
      call readin(labot,latop,nabot,natop,lnabtop,nnbtop,energy,de,
     >   ntstart,ntstop,lstart,lstop,lttop,nbnd,abnd,npsbnd,albnd,alpha,
     >   ncstates,npot,lpot,lactop,formcut,regcut,expcut,theta,
     >   ifirst,isecond,nold,nqm,qcut,rmax,fast,ldw,nk,sk,nznuc,nze,
     >   itail,corep,r0,ninc,linc,npar,nent,zasym,nunit,ndbl,nps,iborn,
     >   ne2e,lslow,lfast,slowe,enion,enlevel,target,projectile,match,
     >   lpbot,lptop,npbot,nptop,npsp,alphap,luba,speed,
     >   erange,inoenergy,pint,igz,igp,analyticd,packed)
      cnode = ch(mod(lstart,10))//'_'//ch(nodeid)
      write(ench,'(1p,"_",e10.4)') energy
      if (nodes.eq.1) cnode = ''
      print*, 'nomporig, nodes, myid and nodeid are:',
     >   nomporig,nodes,myid,nodeid

c
      call date_and_time(date,time,zone,valuesin)
      print '(/,i4,": nodeid completed readin at: ",a10)',
     >      nodeid,time
      inofcexch = 0
      nfcmax = 100
!      print*, 'Enter 1 if NO fcexch() routine to be used'
!      read*,  inofcexch, nfcmax
!      print*,' inofcexch=', inofcexch
!      print*,'  nfcmax=', nfcmax
c
      energy_for_semirelBorn = energy
      Z = dble(nznuc)
      energy = energy + myid * de
      second = isecond.eq.-2
      if (second)print*,'Will be using the second order approximation'
      if (isecond.eq.-3) print*,'Will be using polarization only'
      ipar = 0
C  Some flags to use the Born or the Born with exchange approximations
      ntype = ntstart
      if (projectile.eq.'photon') then
         if (ntype.eq.0) then
            print*,'Will be using Hylleraas-20 helium ground state'
         else 
            print*,'Will be using an MCHF ground state'
         endif
         if(NZNUC.eq.3.and.ZASYM.eq.1.0) then
            lithium = 1
            print'(A)', 'DPI of lithium'
         else
            lithium = 0
         end if         
      else 
         if (ntype.eq.-1.and.nqm.eq.1) then
            print*,'Will use the Born approximation'
c$$$         elseif (ntype.eq.-2.and.nqm.eq.1) then
c$$$            print*,'Will use the Born approximation with exchange'
c$$$         elseif (ntype.eq.-3) then
c$$$            print*,'Will use a local core potential in file klaus.dat'
c$$$         elseif (ntype.eq.-4) then
c$$$            print*,'Will use only the Hartree core potential'
c$$$         elseif (ntype.eq.-5) then
c$$$  print*,'Will drop the ETOT term prportional to THETA'            
         endif
C ANDREY----------------------------------------------------------------
         if (nze.eq.1.and.(nznuc.eq.3.or.nznuc.eq.11.or.nznuc.eq.19.or
     $        .nznuc.eq.37.or.nznuc.eq.55.or.nznuc.eq.87)
     $        ) then
            alkali=.true.
c$$$            ntype = -3 
c$$$            print*," p-Alk case: ntype is set to -3 to block exchange"
         else            
            alkali=.false.
         end if         
      endif                     ! projectile

C  RAV--------
         IF (nznuc.eq.12) then
            helike=.true.
         else
            helike=.false.
         end if         
      if (myid.le.0) print*, 'He - like', helike
      
C  zeff is used in the calculation of the projectile plane/distorted/
C  coulomb waves.
      zeff = zasym * nze
      
C Define the projectile energy in Rydbergs
      ERY = ENERGY / RY
      
C  Setup the FACLOG array
      CALL ARRANGE
C  Setup the R integration grids
      call grids(qcut, ndbl, rmax, rmesh, maxr, meshr, jdouble, 40, id)

      rmudeb = 0.0
      minvdc = 1
      maxvdc = 0
      inquire(file='Debye',exist=dbyexists)
      if (dbyexists) then
         open(42,file='Debye')
         read(42,*) dblp, zdby
         rmudeb = 1d0/dblp
         print*,'Debye radius and Z:',dblp, zdby
         close(42)
         vdcore(1:meshr,0) = -(zdby-1.0)/rmesh(1:meshr,1)*            
     >      exp(-rmesh(1:meshr,1)/dblp)
         do la = 1, lamax
            vdcore(1:meshr,la) = vdcore(1:meshr,0)
         enddo
         vdcore_st(:,:) = vdcore(:,:)
         vdcore_pr(:,:) = - nze * vdcore(:,:)
         minvdc = 1
         maxvdc = meshr
      end if

      z=float(nznuc)
      z0 = 1.0
C DON'T Define the equivalent local potential felt by the valence electron
c$$$      call ipmhf(z,z0,vnucl,maxr,rmesh(1,1),meshr,expcut)
C  Above routine returns the potential * R, and in a.u. Convert it to Rydbergs
      if (isecond.lt.0) then
         lastop = latop
      else
         lastop = max(lnabtop, latop)
      endif 
      rmax = rmesh(meshr,1)
C  Define arrays containing the powers RMESH(I,1), as well as 
C  the CNTFUG array which is LN*(LN+1)/X/X.
      call setpow(lstop,ltmax)

      do i = 1, meshr
C  UCENTR is an array which is used only in the routine COULOMB which
C  makes continuum waves for the single electron of the atom/ion, e.g.,
C  for atomic hydrogen, He+, Li++,...
c$$$         ucentr(i) =  - 2.0 * (zasym + 1.0) * rpow2(i,0) ! Debye screening
         ucentr(i) =  - 2.0 * (zasym + 1.0) / rmesh(i,1)
c$$$         vnucl(i) = 2.0 * vnucl(i) / rmesh(i,1)
c$$$         vnucl(i) = 2.0 * (vnucl(i) - z0) / rmesh(i,1)
      end do
c$$$      i = meshr
c$$$      do while (i.gt.1.and.abs(vnucl(i)).lt.1e-10)
c$$$         i = i - 1
c$$$      end do
c$$$      maxpot = i
c$$$      print*,'maxpot:',maxpot

      if (projectile.eq.'photon') then
c$$$         call     ATOM(rmax, meshr, rmesh)
c$$$         if (ntype.eq.-3) then
c$$$            call corrATOM(rmax, meshr, rmesh)
c$$$         elseif (ntype.eq.-2) then
c$$$            call mcorrATOM(rmax, meshr, rmesh)
c$$$         elseif (ntype.eq.-1) then
c$$$            call ncorrATOM(rmax, meshr, rmesh)
c$$$         else 
c$$$            call xhylleraas
c$$$         endif

         if (ntype.eq.0) then
            iSpin = 0
            call xhylleraas
c$$$            call jones
c$$$            call lesech
         else if(ntype.eq.-1) then !H2 Slater basis
            call H2
         else if(ntype.eq.1) then
            call corrATOM(rmax, meshr, rmesh)
         endif 
      endif
      
      do i = 1, meshr
C  UPLANE will be used in calls to FIRST whenever there is no distorting
C  potential to subtract
         uplane(i) = 0.0
C  VASYMP will be used in calls to MAKECHIL. It will be used whenever
C  there is no distorting potential, otherwise U will be used. Both of
C  these may contain the Coulomb part for non-zero ZEFF.
c$$$         vasymp(i) = zeff * 2.0 * rpow2(i,0) ! generalized for Debye screening
         vasymp(i) = zeff * 2.0 / rmesh(i,1)
      end do
C  We assume that all targets are H like, except for He
      hlike = mod(nznuc - nint(zasym),2) .eq. 1
      print *, 'HLIKE:',hlike
C  For two-electron targets get the one-electron functions for the ion
      i_sw_ng = 0   ! for hlike==false targets make a separate case for noble gases
      inquire(file='noblegas_yes',exist=exists)
      if(exists) then
         i_sw_ng = 1
         i_noblegas_model_pot = 0  ! do not use model potential
         open(735,file='noblegas_yes')
         read(735,*) i_noblegas_model_pot
         close(735)
      endif
      
      orzasym = zasym
      if (.not.hlike) then
         if(i_sw_ng .eq. 1)  then 
            zasym = zasym + 2*(2*linc + 1) - 1
            print*,"In Dmitry's Neon part"
         else
            zasym = zasym + 1.0            
         endif
      endif
 
      if (alkali) call ubbgausleg(lmax1, nmax, xgz, wgz, pleg) ! Andrey         
            
      if (nznuc-nint(zasym).gt.2) then
C  Make the core states
         call makecorepsi(nznuc,zasym,ry,corep(0),r0(0))
C  Get the direct core potential.
         call makevdcore(temp,minvdc,maxvdc,nznuc,uplane)
         print*,'MAXVDC, MESHR:',maxvdc,meshr
         temp(maxvdc+1:meshr) = 0.0
           
C ANDREY: for pos-alkali and el-alkali without exchange
       ! if (alkali.or.ntype.eq.-3) then
         if (ntype.eq.-3) then
            if(i_sw_ng .eq. 0) then ! not noble gas   ! DFursa
               if (nznuc.eq.3) then
                                !     E_2s (eV)
C                               !-------------------------------
                  epar = -0.3064 ! -5.34156764 HF-calculation
                  epar = -0.3831 ! -5.39206788 experimental value 
               else if (nznuc.eq.11) then               
!     epar = -0.26405341010256356 ! 
                  epar = 0.0282819 ! E_3s = -5.139 eV
               else if (nznuc.eq.19) then               
                  epar=-0.14238185 ! K: E_4s = -4.341 eV   
               else if (nznuc.eq.37) then
                  epar=3.5      !in Ryd ==  Rb: E_5s = -4.177 eV
               else if (nznuc.eq.55) then
                  epar=-0.14    ! Cs: E_6s = -3.894 eV
               else 
                  if (nznuc.eq.12) then
                     epar = 0.681008
                  else 
                     print*, 'main: *** warning *** : specify epar'
                     print*, 'main: target structure may be inaccurate' 
                     epar = 0.0
                     stop
                  end if
               end if
            else                ! noble gas
               if(nznuc.eq.10) then ! Ne
                  epar = -0.1   ! not used (DF)
               endif
            endif
         end if

C RAV-----------         
         if (helike) then
            if (nznuc.eq.12) then
               epar = 0.681008
            else
               print*, 'main: *** warning *** : specify epar'
               print*, 'main: target structure may be inaccurate'     
               epar = 0.0
               stop
            end if                   
         end if


C ANDREY: for pos-alkali and el-alkali without exchange         
       ! if (alkali.or.ntype.eq.-3) then            
         if (ntype.eq.-3.and. i_sw_ng .eq. 0) then
               
C     extrapolate temp for i > maxvdc
c     andrey: makevdcore returns nonzero temp(i) 
c             when minvdc <i< maxvdc
            
            if (alkali) then
             
C     extrapolate temp for r > rmesh(maxvdc,1)
C     this function form is used for the reason that
C     the tail integral in functions Ql could be calculated
C     analytically      
               xx1 = rmesh(maxvdc-1,1)
               yy1 = temp(maxvdc-1) ! *rmesh(maxvdc-1,1)**2
               xx2 = rmesh(maxvdc,1)
               yy2 = temp(maxvdc)   ! *rmesh(maxvdc,1)**2              
               b1 = log(yy2/yy1)/(xx2-xx1)
               a1 = yy2/exp(b1*xx2)
               do i=maxvdc+1,meshr
                  r = rmesh(i,1)                  
                  temp(i) = a1*exp(b1*r)  ! /r/r
               end do
c$$$C----- TEST ----------------------------------------
c$$$!$omp critical          
c$$$      inquire(file='Ql',exist=dbyexists)      
c$$$      if (dbyexists) then
c$$$         print*,"a1, b1:", a1, b1          
c$$$         print*," maxvdc, rmesh(maxvdc,1): ",maxvdc,rmesh(maxvdc,1)
c$$$         open( unit = 55,file='tempr2.dat')
c$$$         do i = 1, meshr
c$$$            write (55,*) rmesh(i,1), temp(i) ! * rmesh(i,1) 
c$$$         end do
c$$$         i=maxvdc-1
c$$$         r = rmesh(i,1)
c$$$         print '("   r1 = ", f15.8, "; f1 = ", e16.10)',
c$$$     1        r, temp(i) ! * r**2
c$$$         i=maxvdc
c$$$         r = rmesh(i,1)
c$$$         print '("   r2 = ", f15.8, "; f2 = ", e16.10)',
c$$$     1        r, temp(i) ! * r**2
c$$$         close (55)
c$$$         stop ' CHECK tempr2.dat' 
c$$$      end if
c$$$!$omp end critical             
c$$$C----- TEST ----------------------------------------                 
            else
               if (nznuc.eq.12) then ! Mg
                  a1 = 0.0974031
                  b1 = 3.01811
c                 y(x) = -a*exp(-b*x)/x ;
                  do i=maxvdc+1,meshr
                     r = rmesh(i,1)
                     temp(i) = -a1*exp(-b1*r)/r
                  end do
               else
                  print*, 'main: specify temp(i) for i > maxvdc'
c$$$                  open(unit=56,file='temp')
c$$$                  print*, ' maxvdc = ', maxvdc
c$$$                  do i=1, maxvdc               
c$$$                    write(56,*) rmesh(i,1), temp(i), temp(i)* rmesh(i,1)
c$$$                  end do
c$$$                  close (56)                
                  stop
               end if           ! nznuc = 12
            end if              ! nznuc = 03 or else 11
         end if                 ! ntype

c$$$         !if (nznuc.eq.12) then
c$$$            open (unit=90,file='temp')
c$$$            do i= minvdc, meshr
c$$$               r = rmesh(i,1)
c$$$               write (90,*) r, temp(i) * r
c$$$            end do
c$$$            close (90)
c$$$            stop 'main.f: see temp'            
c$$$            !open (unit =55,file = 'epar')
c$$$            !read (55,*) epar
c$$$            !close (55)            
c$$$         !end if
         
         
C ANDREY: for pos-alkali and el-alkali without exchange         
         vexch(:)= 0.0
!        if ((ntype.eq.-3).or.alkali) then  
         if (ntype.eq.-3) then
            if(i_sw_ng .eq. 0) then ! not noble gas   ! DFursa            
            call makevexch(lactop, lamax, nabot, lnabmax, nnmax, ltmax,
     $         maxr, istoppsinb, psinb, rpow2, temp, nznuc, epar, 
     $           maxvexch, vexch)
            
c$$$            open(unit=56,file='vexch')
c$$$            print*, 'VEXCH: maxvexch,rmesh(maxvexch,1):'
c$$$            print*, maxvexch,rmesh(maxvexch,1)
c$$$            print*
c$$$            do i=1,meshr                  
c$$$               write(56,*) rmesh(i,1), vexch(i), temp(i)
c$$$            end do
c$$$            close (56)
c$$$            stop            
            else
               call makevexch_DF(lactop,lamax,nabot,lnabmax,nnmax,
     $              ltmax,maxr,istoppsinb,psinb,rpow2,temp,nznuc,epar, 
     $              maxvexch, vexch)
            endif
         endif ! ntype
         
         if (alkali.or.helike.and.nze.eq.1) then            
C           parameters for for common blocks: 
            maxr1=meshr;  znuc1=dble(nznuc-1) 
           if(helike)  znuc1=dble(nznuc-2)
            Pi4 = 4.0 * Pi;
C           select the  potential to use (see makepots)
            icase(1) = 1; true_potential(1) = .true.
            icase(2) = 1; true_potential(2) = .true.            
C           make tables of:
C           - va(r) and ve(r) for the e-core and p-core potentials 
C           - fourier transforms of ve(r)/r used in funlegnum
            veold(:)=0.0
            vaold(:)=0.0
            
            print'(a18,$)', '   main: makepots: '
c$$$            if (alkali) then
c$$$c     andrey: makevdcore returns zero temp(i) for i > maxvdc.
c$$$c             it results leads to the error when Ubb integrals
c$$$c             are calculated. i fix it bu extrapolatiing 
c$$$c             temp for this region
c$$$c               print*,'ANDREY: main: nznuc', nznuc
c$$$               if (nznuc.eq.3) then
c$$$                  print*,'extrapolating temp'
c$$$                  a1 = 1.65483437974677; ! this is for Li
c$$$                  b1 = 4.23186257292008
c$$$c              f(x)=-a*exp(-b*x)/x
c$$$                  do i=maxvdc+1, meshr
c$$$                     r = rmesh(i,1)
c$$$                     old = temp(i)
c$$$                     temp(i) = -a1*exp(-b1*r)/r
c$$$                     if (i.lt.maxvdc+10) print*, i, old, temp(i) 
c$$$                  end do
c$$$               else
c$$$                  print*, 'main: specify temp(i) for i > maxvdc'
c$$$                  stop
c$$$               end if               
c$$$            end if

                   
            call makepots(icase(1), true_potential(1), rmesh(1:meshr,:),
     $           meshr, temp(1:meshr),vexch(1:meshr),lamax,corep, r0,
     $           maxx,dx, xmesh, rmesh1, va, vaM, ve, veold(1:meshr)
     $           ,vaold(1:meshr)) !Rav:added vaold need for He-like targets           
c$$$c     andrey: now i restore temp(i) for i > maxvdc
c$$$            temp(minvdc+1:meshr) = 0.0            
            print*, 'done'
            
C           vefun: for funleg2 to calculate Q_l integrals
            ctail = corep(0)/2.0
            corepp = corep(0)
            r0p = r0(0)
            vefun = 0.0
            do i =1, meshr
               r = rmesh(i,1)
               pp =  polpot(r,corep(0),r0(0))
               !vefun(i) = 8d0*Pi*veold(i)*rmesh(i,1)*rmesh(i,3)
               r2dr = r*r*rmesh(i,3)
               vefun(i,1) = 8d0*Pi * (-temp(i)) * r2dr
               vefun(i,2) = 8d0*Pi * pp * r2dr
            end do            
            ! ntype=-3            ! to block exchange - do I need it?  
         end if                 ! alkali

         vdcore0(:,:) = vdcore(:,:) ! ANDREY            
         
         do la = labot, latop
            if (corep(la).gt.0.0) then
               minvdc = 1
               maxvdc = meshr
            endif 
            
            if (alkali.and.(.not.true_potential(2))) then
               stop 'Andrey: do not run for phen pots'
               do i = 1, meshr  ! phenomenological vdcore for Na
                  r = rmesh(i,1)                  
                  call phenpots(icase(2), r, va1, ve1)
                  vdcore(i,la) = va1
                  vdcore0(i,la) = va1
               end do               
            else               
               do i = 1, meshr
C  Phenomenological polarization potential at rmesh(i,1) in a.u. is polpotr
                  polpotr = polpot(rmesh(i,1),corep(la),r0(la))
c$$$               if (nznuc.eq.19) polpotr =
c$$$     >            -5.4605/2.0/rmesh(i,1)**4
c$$$  >            *(1.0-exp(-(rmesh(i,1)/2.22786)**2.871))

                  vdcore0(i,la) = - nze * temp(i) + polpotr ! ANDREY      
                  vdcore(i,la) = temp(i) + polpotr + vexch(i)
                  vdcore_st(i,la) = temp(i) + polpotr + vexch(i)
                  vdcore_pr(i,la) = - nze * temp(i) + polpotr 
               enddo
c$$$c$$$           if (ntype.eq.-3) then
c$$$c$$$               call makevexch(lactop, lamax, nabot, lnabmax, nnmax, 
c$$$c$$$     $            ltmax, maxr, istoppsinb, psinb, rpow2, vdcore(1,la), 
c$$$c$$$     $            nznuc, E, maxvexch, vexch)
c$$$c$$$               do i = 1, meshr
c$$$c$$$                  vdcore(i,la) = vdcore(i,la) + vexch(i)
c$$$c$$$               enddo
c$$$c$$$           endif
            end if            
         enddo
         do la = latop + 1, lamax
            vdcore0(:,la) = vdcore0(:,latop)
            vdcore(:,la) = vdcore(:,latop)
            vdcore_st(:,la) = vdcore_st(:,latop)
            vdcore_pr(:,la) = vdcore_pr(:,latop)
         enddo 
         if( i_sw_ng .eq. 1 .and. i_noblegas_model_pot .eq. 1) then

            if(nznuc.eq.10) then  ! Ne
               A1 = 2.91
               B1 = 5.5
c     print*,'Enter A, B  for neon model potential'
c     read*, A0, B0
               A0 = 3
               B0= 5
               A2 = A0
               B2 = B0
               do i = 1, meshr
                  rr = rmesh(i,1)
                  tmp1 = -A1*exp(-rr*B1)
                  tmp0 = -A0*exp(-rr*B0)
                  tmp2 = -A2*exp(-rr*B2)
                  vdcore(i,0) = vdcore(i,0) +  tmp0
                  vdcore(i,1) = vdcore(i,1) +  tmp1
                  vdcore(i,2) = vdcore(i,2) +  tmp2
               enddo
            endif
            if(nznuc.eq.18) then ! Ne
               A1 = 3.2
               B1 = 4.58
c     print*,'Enter A, B  for argon model potential'
c     read*, A0, B0
c     No point to fit all energy levels as singlet-triplet mixing is just too large for argon...
               A0 = A1
               B0= B1
               A2 = A0
               B2 = B0
               do i = 1, meshr
                  rr = rmesh(i,1)
                  tmp1 = -A1*exp(-rr*B1)
                  tmp0 = -A0*exp(-rr*B0)
                  tmp2 = -A2*exp(-rr*B2)
                  vdcore(i,0) = vdcore(i,0) +  tmp0
                  vdcore(i,1) = vdcore(i,1) +  tmp1
                  vdcore(i,2) = vdcore(i,2) +  tmp2
               enddo

            endif

         endif         
      else
C  Here for only H, He+, Li++,... 
         do la = labot, latop
            if (corep(la).ne.0.0) print*,'Setting corep(la) to zero'
            corep(la) = 0.0
         enddo 
      end if
c
c     Two-electron polarization potential for one-electron atoms.
c     Defined here and used in scattering part two modify Coulomb
c     potential. Requiered for heavy atoms to get correct high energy
c     behavior for S-P DCS (or GOS).
c     Write in the same common block  /di_el_core_polarization/ that is
c     used for two-electron atoms. However, this messes up the spin asymmetries
C     and we are not too sure how best to proceed.
      gamma = 0.0
      inquire(file='withdiel',exist=exists)
      if(hlike.and.exists) then
         gamma = corep(1)
         open(42,file='withdiel')
         read(42,*) r0_di_el         !r0_di_el = r0(1)-0.05
         close(42)
         print*,'Modifying POLPOT using di-electronic term, g,r0:',
     >      gamma, r0_di_el
         print*,'This messes up elastic spin asymmetries for e-Cs'
         print*,'See the 13 eV case on APAC' ! haven't played with r0_di_el
         do i = 1, meshr
c  2.0 in polpot(...) is to eliminate 1/2.0 in definition of polpot(...).
c  gamma is used in vmat.f via a common block
            tmp  = polpot(rmesh(i,1),2.0,r0_di_el) 
            pol(i) = sqrt(abs(tmp))
         enddo
      end if
c     
C  The following code makes sure that the L of the incident atomic state
C  gets calculated first, in order to get the correct total energy.
      do l = 0 , lastop
         lmatch(l) = l
      enddo
      lmatch(0) = linc
      lmatch(linc) = 0
      
      do l = 0 , lamax
         do n = 1, ncmax
            ovlp(n,l) = 1.0
            do nn = 1, nnmax
               ovlpnl(n,l,nn) = 0.0
            enddo
         enddo 
      enddo
C
      vrel(:) = 0.0
      rrel(:) = 0.0
      if(nznuc .eq. 80 .or. nznuc .eq. 79) then    ! Hg-target or Au         
         print*, 'Enter 0 if no rel.cor. to be included for Hg or Au'
         read(*,*) iih         
         if(iih.ne.0)  then
            open(648,file='TMP.A.B')
            read(648,*) A,B
            print*, 'A, B=', A, B
            read(648,*) rr1,rr2
            print*, 'rr1, rr2=', rr1, rr2
            rrel(0) = 1.0
            rrel(1) = rr1
            rrel(2) = rr2
            close(648)
            do i=1,meshr
               r = rmesh(i,1)
               vrel(i) = A*exp(dble(-B*r))!  *(1.0-exp(dble(-B*r)))
!               vrel(i) = polpot(r,corep(0),r0(0))
               if(abs(vrel(i)) .lt. expcut) then
                  print*, 'i=', i, r
                  exit
               endif
               if(i .ge. 4000-1) then
                  print*, 'Increase r-grid dimension in common ',
     >               'block /vrelcor/, files main.f, hevmat.f'
                  print*, 'i=', i, r
                  exit
                  stop
               end if
            end do
            iCR = i-1
         do l=0,2            
            vdcore(1:iCR,l) = vdcore(1:iCR,l) + rrel(l)*vrel(1:iCR)
         end do
         end if
      end if

C     Define the target states      
      nchanmax = 0
      nstmax = 0
      do la = 0, lastop
         if (hlike) then
            nnbin = nnbtop
         else
            nnbin = nabot(la) + nold + 1
            nnbin = min(8,nnbtop)  ! Igor added this to stop calls for large l
         endif
c$$$         nnbin = max(abs(nnbtop),nold)
         l = lmatch(la)
         if (i_sw_ng.eq.0) then
            print*,'Defining eigenstates with COREP and R0 to N:',
     >         corep(l),r0(l),nnbin
            call makeeigenpsi(nznuc,zasym,nnbin,l,corep(l),r0(l),
     >         ery,ry,enion,enlevel,ntstop,ovlpnl,nold)
         endif 
         npsstates(1,la) = 0

         etot = ery + enpsinb(ninc,linc)
         if (hlike.and.-slowe(1)/ry .gt. etot)
     >      slowe(1) = - etot * ry / 2.0
         
         if (nold.gt.0) then
            if (l.gt.0.and.alpha(l).lt.0.0) then
               al = -al
            else 
               al = alpha(l)
            endif 
            if (nps(l).gt.0) then
               npstat = nps(l) - l
            else
               npstat = natop(l) - l
            endif 
            print*,'Defining pseudostates with N and Alpha:',
     >         npstat, al
            if (ne2e.eq.0) slowe(1) = 0.0
            if (hlike) then
               ine2e = abs(ne2e)
            else
               ine2e = 0
            endif
            call makepspsi(nold,nznuc,zasym,npstat,nnbin,l,corep(l),
     >         r0(l),al,vdcore_st(1,l),ovlp,phasen,ery,ry,enion,enlevel,
     >         slowe,ine2e,ovlpnl,hlike,orzasym)
         endif
                  
         print*
c$$$         if (ipar.eq.0) then
            nchanmax = nchanmax + (natop(l) - nabot(l) + 1) * (l + 1)
            nstmax = nstmax + natop(l) - nabot(l) + 1
c$$$         else
c$$$            nchanmax = nchanmax + (natop(l) - nabot(l) + 1) * l
c$$$            if (l.gt.0) nstmax = nstmax + natop(l) - nabot(l) + 1
c$$$         endif 
         if (nchanmax.gt.nchan.and.l.eq.lstop.and.hlike) then
            print*,'NCHAN should be at least NCHANMAX',nchan,nchanmax
            stop 'NCHAN should be at least NCHANMAX'
         endif
            
c$$$         do n = nabot(l), min(natop(l),9)
c$$$            do i = 1, istoppsinb(n,l)
c$$$               write(70+l*10+n,*) rmesh(i,1), psinb(i,n,l)
c$$$            enddo
c$$$            close(70+l*10+n)
c$$$         enddo

c$$$         open(99,file=ch(l))
c$$$         write(99,'("#N =",i2,"  l =",i2," lam:",f4.1,
c$$$     >      "  E(n):",200e14.6)') natop(l)-nabot(l)+1,l,
c$$$     >      al*2.0,(enpsinb(n,l)*ry,n=nabot(l), natop(l))
c$$$         write(99,'(  "#     r(i)          dr(i)",
c$$$     >      "      psi(i,n), n = 1, N")')
c$$$         do i = istoppsinb(natop(l),l),1,-10
c$$$            write(99,'(200e14.6)') rmesh(i,1), rmesh(i,3),
c$$$     >         (psinb(i,n,l),n=nabot(l), natop(l))
c$$$         enddo
c$$$         close(99)

         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C______ANDREY: Get basis functions for TDSE calculations
         inquire(file='TDSE',exist=timeexists)
         if (timeexists.and.la.le.latop) then
            if (la.eq.0) print
     1           '(" Ps-states for TDSE calculations, latop = "
     2           ,I2.2," : ",$)'
            print'(I2.2,1x,$)', la            
            call savepsold(la)        
            if (la.eq.latop) stop ' ... done'
         end if
C______ANDREY:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         
      end do
      
C     ANDREY: vdcore has to be without the exchange part
C             after the target structure is done
C           For pos-alkali and el-alkali without exchange      
C     if (alkali.or.ntype.eq.-3) vdcore = vdcore0
c$$$      if (ntype.eq.-3) vdcore = vdcore0
      
C     The valence electron is assumed to be initially in the state
C     NINC, LINC.
      etot = ery + enpsinb(ninc,linc)
      if(nznuc .eq. 80) then
c         call  iterimp_rel(vdcore,corep,r0)
      end if
      
c     
      if (.not.hlike .or. nznuc .eq. 55) then
         call one_el_V_DERIV(vdcore) ! get fine structure splitting for positive ion 
      end if
c     
      if (.not.hlike) then
C  Here for helium-like targets. Will redefine ETOT.
c$$$         if (nze.eq.1) stop 'set OVLPNL properly'
         do l = 0 , lamax
            do n = 1, ncmax
               ovlp(n,l) = 1.0
               do nn = 1, nnmax
                  ovlpnl(n,l,nn) = 0.0
               enddo 
            enddo 
            corep(l)=0.0
            r0(l) =1.0
         enddo

         if(i_sw_ng .eq. 1) then
            print*,'Entering NOBLE2EL'
            zasym = orzasym
            call noble2el(nmaxhe,namax,nchanmax,ery,etot,
     >           enion,enlevel)
         else
            print*,'Entering MAINHE'
            zasym = orzasym
            rmax=rmesh(meshr,1)
            nptopin(0:lptop) = nptop(0:lptop)
            nptop(:) = 0 ! no Ps states when calling He-like structure
            call mainhe(nmaxhe,namax,pnewC,ery,etot,
     >           lastop,abs(nnbtop),ovlp,phasen,regcut,expcut,ry,
     >           enion,enlevel,enionry,nchanmax,
     >           ovlpnl,slowe,abs(ne2e),vdcore_st)
            nptop(0:lptop) = nptopin(0:lptop)
         endif
         nstmax = nmaxhe

      endif 
c$$$      if (projectile.eq.'photon'.and.nabot(1).eq.3) then
c$$$         nabot(1)=2             !anatoli
c$$$         print*,'putting 2p into the core'
c$$$      endif 

      
C  Define positronium states
      nposchmax = 0
      nposstmax = 0
      boxps = .false.
      natopi = natop(0)
      do lp=lpbot,lptop
         enpsinb(natop(lp)+1:nnmax,lp)=0.
      enddo     
      dppoltot = 0.0
      do lp = lpbot, lptop
         print*,'Defining positronium eigenstates for LP:', lp
         print*,' l Np N      e(Ry)     e(eV)     ovlp'
         nstart = natop(lp)
! Ps states must start after all He-states. When enpsinb=0.0 will cycle in getchinfo
         if (.not.hlike) nstart = max(natopi,natop(lp))

         n = nstart
         do npos = npbot(lp), abs(nnbtop)
            n = n + 1
            call rnl(0,npos,lp,psinb(1,n,lp),enpsinb(n,lp),
     >         istoppsinb(n,lp))
            tsum = 0.0
            do i = 1, istoppsinb(n,lp)
               tsum = tsum + psinb(i,n,lp) * psinb(i,n,lp) * rmesh(i,3)
            enddo
            ovlpnl(n,lp,npos) = tsum
            ovlp(n,lp) = tsum
            print'(3i3,3f10.5)',lp,npos,n,enpsinb(n,lp),
     >         enpsinb(n,lp)*ry,tsum
         enddo
         npsstates(2,lp) = 0
!         if (nptop(lp).le.0) then
            print*,'Defining positronium pseudostates for LP, NP:', lp,
     >         npsp(lp)-lp
            maxupl = 1
            npsstates(2,lp) = npsp(lp)-lp
            rlambda(2,lp) = abs(alphap(lp) * 2.0)
            ps2(:,:) = 0.0
            psen2(:) = 0.0
            if (rlambda(2,lp).gt.10.0) then
               print*,'Using Box-based states'
               boxps = .true.
               hmax = rmesh(meshr,2)
               zas = 1.0/sqrt(2.0)
               zas = 0.5
               ra = alphap(lp)
               nbmax = npsp(lp)
               call pseudo(jdouble,id,hmax,zas,lp,ra,nbmax,maxr,vnucl,
     >            psen2,ps2,jmax)
               minps2(:) = 1
               maxps2(:) = jmax
               psen2(:) = psen2(:) * 2.0
               print*,'R-check:',jmax,rmesh(jmax,1)
            else 
               call makeps(0.5, .true., abs(alphap(lp)), lp, expcut,
     >            npsp(lp)-lp, ps2, psen2, minps2, maxps2,rmesh,uplane,
     >            maxupl, meshr,cknd(1,nstart+1,lp))
               if (alphap(lp).lt.0.0.and.etot.gt.0.0) then
                  n = npsp(lp)-lp
                  do while (psen2(n).gt.etot)
                     n = n - 1
                  enddo
                  if (abs((psen2(n)+psen2(n+1))/2.0-etot).gt.
     >               abs((psen2(n-1)+psen2(n))/2.0-etot)) n = n - 1
                  test = (psen2(n)+psen2(n+1))/2.0 - etot
                  print*,test,psen2(n),etot,psen2(n+1)
                  alphanew = - alphap(lp)
                  iter = 1
                  do while (abs(test).gt.1e-4.and.iter.lt.100)
                     iter = iter + 1
                     alphanew = alphanew - test/etot/10.0
                     print*,'New alpha:',alphanew
                     ps2(:,:) = 0.0
                     psen2(:) = 0.0
                     call makeps(0.5, .true., alphanew, lp, expcut,
     >                  npsp(lp)-lp, ps2, psen2, minps2, maxps2, rmesh, 
     >                  uplane,maxupl, meshr,cknd(1,nstart+1,lp))
                     test = (psen2(n)+psen2(n+1))/2.0 - etot
                     print*,test,psen2(n),etot,psen2(n+1)
                  enddo
                  if (iter.ge.100) stop '100 iterations failed'
                  rlambda(2,lp) = alphanew * 2.0
               endif
            endif 
c$$$            if (nptop(lp).eq.0) then
c$$$               nptop(lp) = npsp(lp)
c$$$            elseif (nptop(lp).le.-100) then
c$$$               n = npsp(lp)-lp
c$$$               do while (psen2(n)*ry.gt.-nptop(lp).and.n.gt.1)
c$$$                  n = n - 1
c$$$               enddo
c$$$               nptop(lp) = n+lp
c$$$            else
c$$$               nptop(lp) = - nptop(lp)
c$$$            endif 
            if (nptop(lp).eq.0) then
               nptop(lp) = npsp(lp)
            elseif (nptop(lp).lt.0) then
               n = 1
               do while(n.le.npsp(lp).and.psen2(n).lt.etot)
                  n = n + 1
               enddo
               nptop(lp) = min(npsp(lp),n-nptop(lp)+lp-2) !-nptop(lp)
            endif
            n = nstart
            print*,' l Np N      e(Ry)     e(eV)     ovlp      proj'
            do npos = npbot(lp), npsp(lp)!nptop(lp)
               if (psen2(npos-lp).lt.etot) then
                  opcl = 'open'
               else
                  opcl = 'closed'
               endif
               if (npos.le.nptop(lp)) then
                  chin = '+'
               else
                  chin = '-'
               endif
               n = n + 1
               proj = 0.0
               if (minps2(npos-lp).ne.1) stop 'MINPS2 .NE. 1'
               if (psen2(npos-lp).lt.0.0) then
                  do npp = npbot(lp), abs(nnbtop)
c$$$  npt = npp + nnbtop - lp !This statement also was used, but why?
                     npt = npp + nstart - lp
                     tsum = 0.0
                     do i = minps2(npos-lp),
     >                  min(maxps2(npos-lp),istoppsinb(npt,lp))
                        tsum = tsum + ps2(i,npos-lp) * psinb(i,npt,lp) *
     >                     rmesh(i,3)
                     enddo
                     ovlpnl(n,lp,npp) = tsum
                     proj = proj + tsum * tsum
                  enddo
                  ovlp(n,lp) = proj
                  if (nnbtop.lt.0) then
                     ovlp(n,lp) = 1.0
                     ovlpnl(n,lp,npbot(lp)+1:abs(nnbtop)) = 0.0
                     ovlpnl(n,lp,npbot(lp)) = 1.0
                  endif 
c$$$                     do npp = npbot(lp), abs(nnbtop)
c$$$                        big = 0.0
c$$$                        nbig = npbot(lp)
c$$$                        if (ovlpnl(n,lp,npp).gt.big) nbi
c$$$                        print*,'ovlpnl(n,lp,npp)',n,lp,npp,
c$$$     >                     ovlpnl(n,lp,npp)
c$$$                        if (ovlpnl(n,lp,npp).gt.0.5) then
c$$$                           ovlpnl(n,lp,npp) = 1.0
c$$$                        else
c$$$                           ovlpnl(n,lp,npp) = 0.0
c$$$                        endif 
c$$$                        proj = 1.0
c$$$                     endif 

               else
                  en = psen2(npos-lp)
                  nznucpos=1
                  call hlikechi(nznucpos,0.5,en,lp,chi,phase,jstart,
     >               jstop,meshr,rmesh,expcut,regcut,corep,r0)
                  proj=sum(ps2(jstart:maxps2(npos-lp),npos-lp)*
     >               rmesh(jstart:maxps2(npos-lp),3)*
     >               chi(jstart:maxps2(npos-lp)))
                  ovlp(n,lp) = proj/sqrt(sqrt(2.0)) !compat with quadrature
                  ovlpnl(n,lp,:) = 0.0
                  phasen(n,lp) = phase
c$$$                  if (abs(en-1.0).lt.0.1)then
c$$$                     do i = 1, maxps2(npos-lp)
c$$$                        write(67,*) i, rmesh(i,1), ps2(i,npos-lp)*proj,
c$$$     >                     chi(i)
c$$$                     enddo
c$$$                  endif 
               endif
               tsum = 0.0
               do i = minps2(npos-lp), maxps2(npos-lp)
                  tsum = tsum+ps2(i,npos-lp)*ps2(i,npos-lp)*rmesh(i,3)
                  psinb(i,n,lp) = ps2(i,npos-lp)
               enddo
               do i = maxps2(npos-lp)+1, maxr
                  psinb(i,n,lp) = 0.0
               enddo
               if (npos.eq.1) nposn = n
               enpsinb(n,lp) = psen2(npos-lp)
               istoppsinb(n,lp) = maxps2(npos-lp)
               if (myid.le.0) 
     >            print'(3i3,4f10.5,2x,a8,a2)',lp,npos,n,psen2(npos-lp),
     >            psen2(npos-lp)*ry,tsum,proj,opcl,chin
         if (lp.eq.1) then
            l = lp
            osc = oscil(enpsinb(n,l),psinb(1,n,l),istoppsinb(n,l),
     >      enpsinb(nposn,0),psinb(1,nposn,0),
     >      istoppsinb(nposn,0),rmesh,meshr)
            dppol = 4.0 * osc / (enpsinb(nposn,0)-enpsinb(n,l))**2
            dppoltot = dppoltot + dppol
            if (myid.le.0) 
     >         print '("         Oscillator strength:",f6.3," dppol:",
     >         f7.2," dppoltot:",f7.2)', osc,dppol,dppoltot
         endif 
            enddo
!         endif
         npl = nptop(lp) - npbot(lp) + 1
         natop(lp) = nstart + npl
         nposchmax = nposchmax + (lp+1) * npl
         nposstmax = nposstmax + npl
      enddo 

* Alisher's addendum
      if (nze.eq.1.and.LPTOP.ge.LPBOT) then
* forms Qlarray with Leg.functions of the 2nd kind
         call Qltable(min(lstop+latop,lstoppos,ltmax),igz,igp)         
c$$$         call Qltable(lstop+latop,igz,igp)         
      endif                

c     Add atom polarization potential for Hg to account for
c     missing polarazability (should be about \alpha=34.4 au, but we get 
c     about \alpha = 21 au from the singlet P-states)
      ipol = 0
      if(nznuc .eq. 80 ) then    ! Hg-target
         print*,' Enter ipol (1 if you want to include polaris.pot.):'
         read(*,*) ipol
         if(ipol .eq. 1) then
         a =  -13.0 / 2.0
         r0Hg = 2.0
            print*,' Enter r0Hg'
            read(*,*) r0Hg
            print*,'Addition to the atom polarization potential: a =',
     >         a*2.0, ', r0Hg =', r0Hg
c     omega = 8.          ! for 25 eV
c            print*,' Enter omega'
c            read(*,*) omega
c            print*,'Addition to the atom polarization potential: a =',
c     >         a*2.0, ', omega =', omega
            r_max = r0Hg*(-log(expcut))**(1./6.)
c     r_max = sqrt(12.0*energy/(omega*omega*(-log(expcut))))
c            print*,'r_max = ', r_max
         atompol(:) = 0.0
c     atompol(i) = a*(1.0-exp(-dble(r/r0Hg)**6))/r**4
c     atompol(i) = a*exp(-12.0*energy/((omega*r)**2))/r**4   - as in  Klaus's paper, is used here.
c     atompol(i) = a/r**4        for  large r
         ipolmin = 1
         do i=1,meshr
            r = rmesh(i,1)
            atompol(i) = a*r*r/(r*r*r + 50.0)**2
c            if(r .lt. r_max) then
c               ipolmin = i+1
c            else
c               atompol(i) = a*exp(-12.0*energy/((omega*r)**2))/r**4
c               if(abs(atompol(i)) .lt. regcut) then
c                 ipolmax = i
c                 exit
c               end if
c            end if
         end do
            ipolmax = meshr
            print*, '  ipolmin =', ipolmin, 'ipolmax =', ipolmax
         else   ! set ipol = 0, for all other values of the entered number
            ipol = 0
         end if
      end if
c      


      
      print*,   'Number of      atomic states and max channels:',
     >   nstmax,nchanmax
      if (nze.eq.1) 
     >   print*,'Number of positronium states and max channels:',
     >   nposstmax,nposchmax
      nent = min(nent,nstmax+nposstmax)
      nchanmax = nchanmax + nposchmax
      if (nchanmax.gt.nchan) then
         print*,'NCHAN should be at least NCHANMAX',nchan,nchanmax
         stop 'NCHAN should be at least NCHANMAX'
      endif

C  For testing purposes it may be useful to run a calculation without
C  any distorting potential
      if (npot.eq.0) then
         do i = 1, meshr
            u(i) = vasymp(i)
            ui(i) = uplane(i)
         enddo
      else 
C  UI is the distorting potential in Rydbergs due to a single valence electron 
C  together with a nuclear term for one neucleon. It is added to the
C  distorting potential U of the projectile. 
C  UI contains that part of the distorting potential which will be subtracted
C  off in calculating the first order matrix elements. U contains the 
C  asymptotic, core, polarization, and ui potentials. It is used in
C  defining the distorted waves.
         if (npot.ge.0) then
            call potent(hlike,npot,lpot,nznuc,-1,ui)
         else
            call potent(hlike,ninc,linc,nznuc,-1,ui)
            alp = npot / 10.0 + 0.1
            print*,'ALPHA:',alp
            do i = 1, meshr
               ui(i) = ui(i) * exp(alp*rmesh(i,1))
            enddo
         endif 
C  for helium, VASYMP and VDCORE should be zero
         do i = 1, meshr
c$$$            u(i) = - nze *(vasymp(i) + ui(i) + 2.0 * vdcore(i,0))
c$$$            u(i) = vasymp(i) - nze *(ui(i) + 2.0 * vdcore(i,0))
            u(i) = vasymp(i) - nze*ui(i) + 2.0 * vdcore_pr(i,0)
         enddo
      endif

      inquire(file='vdcore',exist=exists)
      if (exists) then
         open(42,file='vdcore')
         do i=1, meshr
            write(42,'(30(e12.4))') rmesh(i,1), rpow2(i,0),
     >         (vdcore(i,la),la=0,lamax)
         enddo
         close(42)
      endif


      
      print'(''Total energy of the collision system:'',
     >   f13.4,'' eV ('',f11.4,'' Ryd )'')', ry * etot, etot
C  Define continuum grid
c$$$      if (natop(linc).gt.abs(nnbtop).and.ntstop.eq.2.and.nold.eq.0) then
c$$$         if (abs(nnbtop)+ncstates.gt.ncmax) stop 'Need to INCREASE NCMAX'
c$$$         call egrid(0,etot)
c$$$         print*,'Will be using exact continuum intermediate states'
c$$$         do l = labot, latop
c$$$            call makepsinc(zasym,nznuc,l,abs(nnbtop),
c$$$     >         psinb(1,1+abs(nnbtop),l),enpsinb(1+abs(nnbtop),l),minps,
c$$$     >         istoppsinb(1+abs(nnbtop),l),corep,r0)
c$$$         enddo 
c$$$      endif 
      nch = 0
      call getchinfo(nch,nchp,0,temp,maxpsi,enpsi,la,na,lp)
      do n = 1, nent
         instate(n) = n
      enddo
      mprev = 0
      if (nent.gt.1) then
         inquire(file='instates',exist=exists)
         if (exists) then
            open(42,file='instates')
            do n = 1, nent
               read(42,'(a3)') chinstate(n)
               m = 1
               do while(m.ne.0)
                  if (chan(m).eq.chinstate(n)) then
                     instate(n) = m
                     print*,'setting initial state to:',m,chan(m)
                     if (mprev.gt.m) then
                        stop 'need to have instates in order'
                     else 
                        mprev = m
                     endif 
                  endif 
                  m = m + 1
                  call getchinfo(m,nchp,0,temp,maxpsi,enpsi,la,na,lp)
               enddo 
            enddo
            close(42)
         endif
      endif 
      incount = 1
      n = instate(incount)
      nsmax = 0
      nchimax = 0
      nentin = nent
      bothpar = .false.
      do while(n.ne.0)
         call getchinfo (n,nchp,0,temp,maxpsi,enpsi,la,na,lp)
         if (n.ne.0) then
            print'(a8,'' scattering on '',a6,i3,a4,'' at'',f11.4,
     >         '' eV ('',f11.4,'' Ryd )'')', projectile,target,
     >         n,chan(n),ry * (etot - enpsi), etot - enpsi
            if (etot.gt.enpsi) then
               bothpar = bothpar .or. la .ge. 1
            endif 
            if (chan(n)(1:1).eq.' '.or.chan(n)(1:1).eq.'t') nsmax = 1
            if (nentin.le.0.and.enpsi.lt.etot) then
               nent = n
            else if (n.eq.nentin) then
               n = 0
            endif
            incount = incount + 1
            n = instate(incount)
            nchimax = nchimax + (la + 1)
         endif 
      enddo
      if (bothpar) then
         if (npar.eq.0) print*,'WARNING: Consider setting NPAR = 1'
      else
         if (npar.eq.1) then
            npar = 0
            print*,'WARNING: NPAR reset to 0'
         endif
      endif 
         
      print*,'Number of initial states:',nent
      if (nze.ge.1.or.projectile.eq.'photon') nsmax = 0
C The following is for photoionization from triplet states
      if (iSpin.eq.1) nsmax = 1
c determine values of lg task is to compute
      mylstart = lstart
      mylstop = lstop
      nlg = lstop - lstart + 1    
c$$$      if (ntasks.gt.1.and.de.eq.0.0) then
c$$$         if(ntasks.gt.nlg)
c$$$     >      stop'reduce the number of requested nodes to nlg'
c$$$         my_nlg = nlg/ntasks
c$$$         mylstart = lstart + my_nlg * myid
c$$$         if (myid .lt. ntasks-1) then
c$$$            mylstop = mylstart + my_nlg - 1
c$$$            print*,'FIRST part,mylstart,mylstop,myid,ntasks',
c$$$     >         mylstart,mylstop,myid,ntasks
c$$$         else
c$$$c$$$        mylstart = lstart + my_nlg*myid + mod(nlg,ntasks)
c$$$c$$$        mylstop = mylstart + my_nlg - 1
c$$$            mylstop = lstop
c$$$            print*,'SECOND part,mylstart,mylstop,myid,ntasks',
c$$$     >         mylstart,mylstop,myid,ntasks
c$$$         endif
c$$$      endif
      if (slowe(1)/ry.gt.etot) then
         print*,'Setting NE2E = 0'
         ne2e = 0
      endif 
      do n = 1, abs(ne2e)
         slowe(n) = abs(slowe(n))
         slowery(n) = slowe(n) / ry
         print'("The two outgoing electron energies (eV) are:",2f10.4)',
     >   slowe(n), (etot - slowery(n)) * ry
      end do 
      do icenergy=1,inoenergy
         if (de.eq.0.0) energy=erange(icenergy) ! energy defined above for MPI
         ery=energy/ry
         etot = ery + enpsinb(ninc,linc)
         write(ench,'(1p,"_",e10.4)') energy
c$$$         if (myid.eq.-1) then
            tfile = 'potl'//ch(lstart)//ench
            csfile = 'totalcs'//ench
c$$$         else 
c$$$            tfile = '/u/igor/potls/potl'//
c$$$     >         ch(mylstart)//ench
c$$$            csfile = '/u/igor/potls/totalcs'//ench
c$$$         endif 
         if (myid.eq.0) then
         ipstart = 0
         call pwrite(nent,instate,nopen,energy,nznuc,zasym,ry,noprint,
     >      ovlp,ovlpn,phasen,phaseq,mylstart,ipstart,lstop,projectile,
     >      target,nsmax,nchimax,nchanmax,abs(npar),hlike,
     >      nunit,vdcore_pr,minvdc,maxvdc,abs(ne2e),slowery,iborn,
     >      BornICS,tfile,abs(nnbtop),ovlpnl,ovlpnn)
         inquire(file='bornstop',exist=exists)
         if (exists) then
            open(138,file='bornstop')
            call Born_tiecs(energy,nent,nstmax,BornICS)
            close(138)
            stop 'bornstop found'
         endif
         endif ! myid.eq.0

      enddo
      if (dbyexists) then       !for shielded H-like ions
         ui(1:meshr) = ui(1:meshr) + 2.0 * vdcore(1:meshr,0)
         vdcore(:,:) = 0.0
      endif 

c$$$      call pwrite(nent,nopen,energy,nznuc,zasym,ry,noprint,ovlp,ovlpn,
c$$$     >   phasen,phaseq,mylstart,lstop,projectile,target,nsmax,nchimax,
c$$$     >   nchanmax,abs(npar),hlike,nunit,vdcore,minvdc,maxvdc,
c$$$     >   abs(ne2e),slowery,iborn,BornICS,tfile,abs(nnbtop),ovlpnl,ovlpnn)

c$$$C  Find the maximum of (e,2e) channels
c$$$      nchtope2e = 0
c$$$      if (ne2e.ne.0) then
c$$$         ne2e = 0
c$$$         print*,'Have set NE2E = 0'
c$$$      endif 
c$$$      if (ne2e.ne.0) then
c$$$         nche2e = 1
c$$$C  Set LG = LSLOW to get the max number of (e,2e) channels        
c$$$         call getche2e(nche2e,lslow,lfast,lslow,lfa,lf)
c$$$         do while (nche2e.ne.0)
c$$$            nchtope2e = nche2e
c$$$            nche2e = nche2e + 1
c$$$            call getche2e(nche2e,lslow,lfast,lslow,lfa,lf)
c$$$         enddo
c$$$         print*,'Max number of (e,2e) channels:', nchtope2e
c$$$         if (nchtope2e.gt.nchane2e) stop 'INCREASE NCHANE2E'
c$$$      endif 
      call date_and_time(date,time,zone,valuesout)
      print '(/,i4,": nodeid structure complete at: ",a10,
     >   ", diff (secs):",i5)',nodeid,time,
     >   idiff(valuesin,valuesout)
	valuesin = valuesout
      
      pi8 = 8d0 * pi * pi
      
C     LPTOP=-1 means that positronium channels are excluded
      if (alkali.and.lptop.ne.-1) then ! pos+alkali problem
         
C     I. COMPUTE Ubb(r1,r2,ll) - spher. harm. expansion for all ll         
            
         select case (interpolate_ubb)
         case (0)               ! no interpolation              
            ubb_max1 = meshr
            ubb_max2 = ubb_max1 
         case (1,2)             ! 1D or 2D interpolation
            ! ubb_max1 = 1000 defined in ccc.in
            ubb_max2 = meshr
            rlgmin = log10(rmesh(1,1))
            rlgmax = log10(rmesh(meshr,1))
            drlg = (rlgmax-rlgmin)/real(ubb_max1-1)                                             
            if (interpolate_ubb.eq.2) then
               ubb_max2 = ubb_max1                  
            end if                           
         end select
         ubb_min1 = 1; ubb_min2 = 1; 
         ubb_min3 = lpbot; ubb_max3 = 2*lptop;
         
         allocate(ubb_res(ubb_min1:ubb_max1,ubb_min2:ubb_max2
     $        ,ubb_min3:ubb_max3))
         
         if (.not.allocated(ubb_res)) stop 'ubb_res is not allocated'
         
         print*; print*,'calculate Ubb:'; print*
         
         if (interpolate_ubb.eq.2) then            
            do irho = 1, ubb_max2
               arholg(irho) = rlgmin+drlg * real(irho-1)
               arho(irho) = 10.0** arholg(irho)
            end do
         end if
                  
         do ll = ubb_min3, ubb_max3            
            print'(3x,a11,i2,$)',' ll = ',ll             
!$omp parallel do default(shared)
!$omp& private(irho,rho,rholg,ir,rlg,r,res0,res1)   
            do irho = 1, ubb_max2, 1                              
               if (interpolate_ubb.eq.2) then
C$$$                  rholg =  rlgmin+drlg*real(irho-1)
C$$$                  rho = 10.**rholg
C                  rholg =  arholg(irho) ! rlgmin+drlg*real(irho-1)
                  rho = arho(irho) ! 10.**rholg
               else
                  rho = rmesh(irho,1)
               end if
               do ir = 1, ubb_max1, 1
                  if (interpolate_ubb.eq.2) then
C$$$                     rlg = rlgmin+drlg*real(ir-1)
C$$$                     r = 10.**rlg                     
                     r = arho(ir)
                  else if (interpolate_ubb.eq.1) then
                     rlg = rlgmin+drlg*real(ir-1)
                     r = 10.**rlg
                  else
                     r = rmesh(ir,1)
                  end if                  
                  call Ubb0(.false., ll, r, rho, res0, res1)
                  ubb_res(ir, irho, ll) = res1 ! 1                     
               end do           ! ir1                  
            end do              ! ix
!$omp end parallel do            
         end do                 ! ll                        
         print'(a5)',' done'

         
C     select potential for atomic structure calculations

C     true_potential=.true.     ! for sodium only
C     icase(2) = 1              ! (1) SFG+POLPOT, (2) SFG, (3) HEWITT
         
C      end if                    ! alkali and.lptop.ne.-1      
C      if (alkali) then ! pos+alkali problem
      
C     II. COMPUTE QL(q,qa) _______________________________________         

         
         index_ql = 0
                         
C     find out what are correct limits minql3 and maxql3 for
c     given parameters
         minql3=min(lstart,abs(lstart - latop));
         minql3=0               ! for test
         maxql3=lstoppos+latop  ! minql3 = 0; maxql3 = 10;        
                  
         allocate(ql_res(maxql1, maxql2, minql3:maxql3))
         if (.not.allocated(ql_res)) stop 'Ql_res is not allocated'
         
         print*; print*,'calculate Ql: ', maxql3-minql3 ; print*;
         
c     qlgmax = log10(qcut);         
         do ll = minql3, maxql3 
            print'(3x,a11,i2,$)',' ll = ',ll
            inquire(file='qlres'//ch(ll),exist=exists)
            if (exists) then
               open(42,file='qlres'//ch(ll),form='unformatted')
               read(42)((ql_res(iqa,iq,ll),iqa=1,maxql1),iq=1,maxql2)
               close(42)         
               index_ql(ll) = 1
               cycle
            endif 
                           
            print '(i4,":",$)', maxql1
!$omp parallel do default(shared)
!$omp& SCHEDULE(dynamic)
!$omp& private(iqa,lgqa,qa,iq,lgq,qq1,dQlp)
            do iqa = 1, maxql1           
C$OMP critical(print)
               print '(i4,$)', iqa
C$OMP end critical(print)
               lgqa = qlgmin + dqlg * dble(iqa - 1)
               qa = 10d0 ** lgqa
!              qa=qqq(iqa)
               do iq = 1, maxql2                  
                  lgq = qlgmin + dqlg * dble(iq-1)
                  qq1 = 10d0 ** lgq
!                  q = qqq(iq)
C                  z = (qa*qa + q*q)/(2d0*qa*q)
C                  call funleg(z, ll, Q0, Qlp) ! i do not need to do it in posvmat now
                                ! dQlp = 0; if (alkali)              
	          if (ll.le.10) then
                     call funleg2(ll, qa, qq1, dQlp);
                  else
                     call funleg3(ll, qa, qq1, dQlp);
                  end if
C                  ql_res(iqa, iq, ll) =  Qlp+(qa * q/pi8) * dQlp                  
                  ql_res(iqa, iq, ll) =  (qa * qq1/pi8) * dQlp                  
               end do           ! iq                           
            end do              ! iqa
!$omp end parallel do                                     
            index_ql(ll) = 1
            if (nodes.eq.1) then
               open(42,file='qlres'//ch(ll),form='unformatted')
               write(42)((ql_res(iqa,iq,ll),iqa=1,maxql1),iq=1,maxql2)
               close(42)
            endif 
         end do                 ! ll                                 
         call date_and_time(date,time,zone,valuesout)
         print '(/,i4,": nodeid exited UBB calculation at: ",a10,
     >      ", diff (secs):",i5)',nodeid,time,
     >      idiff(valuesin,valuesout)
      end if                    ! alkali

            
C  Start the calculation
C  Define arrays that are used to speed up the calculation of the polorization
C  potential
      do nchf=1,nchan
         do nchi=1,nchan
            do i2 = 0, 1
               nchns(i2,nchi,nchf)=ncstates
               small(i2,nchi,nchf)=.false.
            end do
         end do
      end do
         
      if (isecond.ge.0) then
C  Define types and the ranges of intermediate states
         if (ntstart.eq.0) then
            print*,'Will be using continuum L2 states'
         else if (ntstart.eq.1) then
            print*,'Will be using exact bound intermediate states'
         end if
         if (ntstop.eq.2) then
C  Define continuum grid 
            print*,'Will be using exact continuum intermediate states'
            call egrid(0,etot)
c$$$  call egrid(npotgf,etot)
         else if (ntstop.eq.3) then
            print*,'Will be using Laguerre basis L2 intermediate states'
         end if
      end if 

      allocate(vmatp(1,0:1))
      vmatp(:,:) = 0.0
      nqmold = nqm              ! not used
      iparmin = 0
      if (npar.eq.-1) iparmin = 1
      call date_and_time(date,time,zone,valuesout)
      print '(/,i4,": nodeid starting partial waves at: ",a10,
     >   ", diff (secs):",i5)',nodeid,time,
     >   idiff(valuesin,valuesout)
      endif ! mod(myid,nomp).eq.0

      call mpi_bcast(mylstart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ipstart, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(mylstop, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(iparmin, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(npar,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

C     Start the partial wave LG (=J total orbital angular momentum) loop
      do 770 lg = mylstart, mylstop
         if (myid.le.0) then
            call date_and_time(date,time,zone,valuesinLG)
            print '("Partial wave LG:",i4," started at:",a10)',
     >           LG,time
         endif
         do ipar = iparmin, min(abs(npar),lg)
            if (lg.eq.mylstart.and.ipar.lt.ipstart) cycle !allow restart for ip=1
            if (mod(myid,nomp).eq.0) then
               nodeid = myid/nomp + 1
         nchtope2e = 0
         nchmaxe2e=1
         ne2e0 = 0
c$$$         if (ne2e.ne.0) then
c$$$            nche2e = 1
c$$$            call getche2e(nche2e,lg,lfast,lslow,lfa,lf)
c$$$            do while (nche2e.ne.0)
c$$$               nchtope2e = nche2e
c$$$               nche2e = nche2e + 1
c$$$               call getche2e(nche2e,lg,lfast,lslow,lfa,lf)
c$$$            enddo
c$$$C  The following line is very important so that the e2e code doesn't get
c$$$C  called when nchtope2e = 0
c$$$            if (nchtope2e.eq.0) then
c$$$               print*,'Resetting NE2E to zero'
c$$$               ne2e = 0
c$$$            endif 
c$$$            nchmaxe2e = nchtope2e*abs(ne2e)
c$$$            if (nchmaxe2e.gt.nchane2e) then
c$$$               print*,'nchmaxe2e>nchane2e; nchtop*abs(ne2e),nchane2e',
c$$$     >           nchmaxe2e , nchane2e
c$$$               stop 'nchmaxe2e.gt.nchane2e'
c$$$            endif 
c$$$         endif
         ndumm = 0

c$$$         if (speed) then
c$$$            ispeed = 2
c$$$            call kgrid(ispeed,nk,sk,etot,gkeep,wk,weightk,nbnd,nqm,lg,
c$$$     >         nchkeep,nchopt,npkeep,nold,ndumm,luba,ifirst,npsbnd,
c$$$     >         abs(nnbtop),hlike,uba,theta,nodes)
c$$$            mv = (npkeep(nchkeep+1)-1) * (npkeep(nchkeep+1)) ! * 4
c$$$            mchi =  meshr * (npkeep(nchkeep+1)-1) ! * 4
c$$$         
c$$$            print'(
c$$$     >         ''Memory (Mb) requested for'',
c$$$     >         '' VKEEP and CHIKEEP :'', 4i5)', 
c$$$     >         mv/250000, mchi/250000
c$$$            call memalloc(ptrvkeep,mv*4)
c$$$            if (ptrvkeep.eq.0) stop 'Not enough memory for VMAT'
c$$$         
c$$$            call memalloc(ptrchikeep,mchi*4)
c$$$            if (ptrchikeep.eq.0) stop 'Not enough memory for CHI'
c$$$C  Define the plane (or distorted by ui + vdcore) waves
c$$$            call clock(s1)
c$$$!            print*, 'calling makechil 1'
c$$$            call makechil(lg,gkeep,wk,qcut,zasym,vdcore_pr,npot,ui,ldw,
c$$$     >         dwpot,npkeep,minchikeep,chikeep,phasel,nchkeep,etot,
c$$$     >         nbnd0,abnd,npsbnd,albnd,sigma,abs(nnbtop),natomps)
c$$$            call clock(s2)
c$$$            print*,'Time to make the KEEP partial wave table:',s2-s1
c$$$            call update(6)
c$$$
c$$$C  Initialise the first order V matrix.
c$$$            call initv(vkeep,npkeep(nchkeep+1)-1,nsmax)
c$$$
c$$$ 
c$$$            
c$$$            call clock(v1)
c$$$            if (nabot(labot).gt.1) then
c$$$               call clock(s1)               
c$$$               call core(0,nznuc,lg,etot,chikeep,minchikeep,nchkeep,
c$$$     >            uplane,ldw,vdcore_pr,minvdc,maxvdc,npkeep,vdondum,
c$$$     >              vkeep,vmatp)            
c$$$               call clock(s2)
c$$$               tc = s2 - s1
c$$$            endif
c$$$
c$$$            if (hlike) then           
c$$$               call first(ispeed,ifirst,second,nold,etot,lg,gkeep,
c$$$     >            npkeep,chikeep,minchikeep,ui,ldw,dwpot,phasel,itail,
c$$$     >            nznuc,nchkeep,nchtope2e,qcut,vdondum,vkeep,theta,
c$$$     >            vdcore_pr,minvdc,maxvdc,lfast,lslow,slowery,td,te1,
c$$$     >            te2,ve2ed,ve2ee,dphasee2e,ephasee2e,ne2e0,nchmaxe2e,
c$$$     >            vmatp,nsmax,
c$$$     >            nchistart,nchistop,nodeid,scalapack,
c$$$     >            vmat01,vmat0,vmat1,ni,nf,nd,nodes,myid)
c$$$            else               
c$$$               call scattering(ispeed,ifirst,theta,nold,etot,lg,gkeep,
c$$$     >            enionry,npkeep,chikeep,minchikeep,vdcore_pr,dwpot,
c$$$     >            nchkeep,nmaxhe,namax,nze,td,te1,te2,t2nd,vdondum,
c$$$     >            vkeep,nsmax,itail,phasel)
c$$$            end if 
c$$$            call clock(v2)
c$$$            print*,'time for KEEP VMAT:',v2-v1            
c$$$            call update(6)
c$$$            call memfree(ptrchikeep)
c$$$         endif ! speed
c
c     This is the end of the speed routine.
c     
c     ispeed is set to 1, because it should be, this can be changed later
c
         ispeed = 1
c
c     this is the beginning of the energy loop
c
c
c$$$         do 775 icenergy=1,inoenergy
c$$$            if (de.eq.0.0) energy=erange(icenergy)
         energy = erange(1)
c
c     calling the clock, to get an initial time for the energy loop
c     and redefining the ery (energy in rydbergs) and the total
c     energy.
c
            call clock(timeenergy1)
            ery=energy/ry
            etot = ery + enpsinb(ninc,linc)

            
            
c$$$            if (lg.eq.lstart) then
               write(ench,'(1p,"_",e10.4)') energy
c$$$               if (myid.eq.-1) then
                  tfile = 'potl'//ch(lstart)//ench
                  csfile = 'totalcs'//ench
c$$$                  do n = 1, abs(ne2e)
c$$$                     if (n.lt.10) then
c$$$                        write(e2efile(n),'("te2e",i1)') n
c$$$                     else
c$$$                        write(e2efile(n),'("te2e",i2)') n
c$$$                     endif
c$$$                     open(88+n,file=e2efile(n))
c$$$                     if (mylstart.eq.0) write(88+n,*)lstop,energy,lslow, 
c$$$     >                  slowe(n), lfast, (etot - slowery(n)) * ry
c$$$                  enddo
c$$$               else 
c$$$                  tfile = '/u/igor/potls/potl'//
c$$$     >               ch(mylstart)//ench
c$$$                  csfile = '/u/igor/potls/totalcs'//ench
c$$$                  do n = 1, abs(ne2e)
c$$$                     if (n.lt.10) then
c$$$                        write(e2efile(n),
c$$$     >                     '("/u/igor/te2e",i1,".",a1,a11)') n,
c$$$     >                     ch(mylstart),ench
c$$$                     else
c$$$                        write(e2efile(n),
c$$$     >                     '("/u/igor/te2e",i2,".",a1,a11)') n,
c$$$     >                     ch(mylstart),ench
c$$$                     endif
c$$$                  enddo
c$$$                endif 
c$$$
c$$$               call pwrite(nent,nopen,energy,nznuc,zasym,ry,noprint,
c$$$     >         ovlp,ovlpn,phasen,phaseq,mylstart,lstop,projectile,
c$$$     >         target,nsmax,nchimax,nchanmax,abs(npar),hlike,
c$$$     >         nunit,vdcore,minvdc,maxvdc,abs(ne2e),slowery,iborn,
c$$$     >         BornICS,tfile,abs(nnbtop),ovlpnl,ovlpnn)
c$$$            endif 
         call date_and_time(date,time,zone,valuesin)
         print '(/,i2,": nodeid calling KGRID at: ",a10)',
     >      nodeid,time
         call kgrid(ispeed,nk,sk,etot,gk,wk,weightk,nbnd,nqm,lg,
     >      nchtop,nchopt,npk,nold,ndumm,luba,ifirst,npsbnd,abs(nnbtop),
     >      hlike,uba,theta,nodes)
         call date_and_time(date,time,zone,valuesout)
         print '(/,i2,": nodeid exited KGRID at: ",a10,
     >   ", diff (secs):",i5)',nodeid,time,
     >   idiff(valuesin,valuesout)

C---- ANDREY ---- make Fourier transforms of psinb and psinb*v -------
         if (nze.eq.1) then            
            !maxgk = npk(2)-npk(1) ! number of kgrid points of the 1st channel
            !pmax = 1.0d3          !      gk(maxgk,1)                        
            if (alkali.or.(interpol.and.nznuc.eq.1)) then              
               if (.not.allocated(pmesh)) allocate(pmesh(0:maxp))
               pmesh(0) = 0.0
               do ip = 1,maxp                            
                  pmesh(ip) = 10d0 ** (-3.0 + dble(ip)*6.0/dble(maxp))
               end do               
               if (.not.allocated(ftps)) allocate(ftps(0:maxp,nnmax
     $              ,0:lnabmax))
               if (.not.allocated(ftpsv)) allocate(ftpsv(0:maxp,nnmax
     $              ,0:lnabmax))
               !error = .false.
               call makeftps()  ! (qcut,pmax);  
               print*, ' makeftps: done ' 
            end if            

c$$$C test ftps ============================================================
c$$$!$omp critical             
c$$$            !if (interpol) then
c$$$            if (.true.) then
c$$$
c$$$c$$$               p = 0.1; p2 = p*p
c$$$c$$$               l = 11; n = 12
c$$$c$$$               call gnlp(1.0d0, l, p2, n, FT0, FT1)
c$$$c$$$               print'(" | num: ",2(e14.7,2x) \)', FT1, FT0
c$$$c$$$               call getftps(p, n, l, FT0, FT1)
c$$$c$$$               print'(" | int: ",2(e14.7,2x) )', FT0, FT1
c$$$c$$$               pause
c$$$
c$$$c$$$               ll = 0;
c$$$c$$$ 40            print*; print*,"ll =",ll
c$$$c$$$               print*,"=================================="
c$$$c$$$               continue               
c$$$c$$$               do iarg = 1,20
c$$$c$$$                  arg =0.4 *  (iarg-1) + 0.1
c$$$c$$$                  call sbessel(arg, ll, besj)
c$$$c$$$                  call SF32D(arg,ll,YYY,IERR); 
c$$$c$$$                  print*, arg, besj, YYY(ll)
c$$$c$$$               end do
c$$$c$$$               print*,"=================================="
c$$$c$$$               pause
c$$$c$$$               ll = ll+1
c$$$c$$$               go to 40
c$$$               
c$$$               icount = 0               
c$$$               !do ip = 1, maxp                                 
c$$$               !   pmesh(ip) = 10d0 ** (-3.0 +  dble(ip) * 6.0/maxp)
c$$$               !end do
c$$$               
c$$$            !if (.false.) then
c$$$
c$$$ 314           continue               
c$$$               print'("calculating for i = ",i2)', icount
c$$$               call system (' rm ftps?')
c$$$               
c$$$            
c$$$               do l = 0, latop
c$$$                                !l = 0
c$$$                  !do n = nabot(l), natop(l)
c$$$                  n = nabot(l)+icount
c$$$                  pos = positron(n,l,npos)
c$$$
c$$$                     if (.not.pos) then
c$$$                        print'("*** l  = ",i2,", n= ", i2 \)', l, n 
c$$$                        open (unit = 99, file='ftps1',position='append')
c$$$                        write(99,'("#  l = ",i2.2, ", n = ",i2.2)'),l,n                  
c$$$                        do i = 0, maxp
c$$$
c$$$                           if (i.eq.0) then
c$$$                              print'(" | ftps, ftpsV: ",2(e14.7,2x) \)',
c$$$     1                             ftps(i,n,l), ftpsv(i ,n,l)
c$$$                           end if
c$$$                           write (99,'(4(e14.7,3x))') pmesh(i),
c$$$     1                          ftps(i,n,l), ftpsv(i ,n,l)
c$$$                        end do
c$$$                        write (99,*)
c$$$                        write (99,*)                  
c$$$                        close (99)                           
c$$$
c$$$                  
c$$$c$$$                        open (unit=99,file='ftps2',position='append')
c$$$c$$$                        write(99,'("# l = ",i2.2,", n = ", i2.2)'),l,n
c$$$c$$$                        p = 1.0e-4
c$$$c$$$                        do i = 1, 70                     
c$$$c$$$                           call getftps(p, n, l, FT0, FT1)
c$$$c$$$                           write (99,'(4(e14.7,3x))') p, FT0, FT1
c$$$c$$$                           p = p * 10**(0.1)
c$$$c$$$                        end do
c$$$c$$$                        write (99,*)
c$$$c$$$                        write (99,*)    
c$$$c$$$                        
c$$$c$$$                        close (99)
c$$$                        
c$$$                        open (unit=99,file='ftps3',position='append')
c$$$                        write(99,'("# l = ",i2.2,", n = ",i2.2)'),l,n
c$$$                        p = 0.0001
c$$$                        do ip = 1, 7000                     
c$$$                           p2 = p*p
c$$$                           call gnlp(1.0d0, l, p2, n, FT0, FT1)
c$$$
c$$$                           if (ip.eq.1) then
c$$$                              print'(" | num: ",2(e14.7,2x) \)',
c$$$     1                             FT1, FT0
c$$$                           end if
c$$$                           
c$$$                           write (99,'(4(e14.7,3x))') p, FT1, FT0
c$$$                           p = p*10**(0.001)
c$$$                        end do
c$$$                        write (99,*)
c$$$                        write (99,*)
c$$$                        close (99)
c$$$                        
c$$$                        print'(" | enpsi, E(l,n): ",2(f10.6,2x))',
c$$$     1                             enpsinb(n,l)/2., -1./(2.0*n**2)
c$$$
c$$$                     end if     ! pos
c$$$                  
c$$$                  !end do
c$$$               end do
c$$$
c$$$               print*
c$$$               do l = 10,latop
c$$$                  n = nabot(l)+icount
c$$$                  pos = positron(n,l,npos)
c$$$                  print'("i = ",i2,", l = ",i2,", n = ",i2 \)',l-10,l,n
c$$$                  print*, pos
c$$$               end do 
c$$$               !print*               
c$$$               pause 'main @ makeftps: see '
c$$$               print*               
c$$$               icount = icount+1               
c$$$               go to 314
c$$$            end if
c$$$!$omp end critical
c$$$C test ftps ================================================================== 
            
         end if                 ! alkali

         
         
                  
C  Uncomment the following line to check that the (e,2e) code gives same
C  same answers as CCC. More changes are necessary in vmat.f and rest.f
c$$$         nchtope2e = nchtop

         if (isecond.lt.0) nchopt = 0
         if (nchopt.gt.nchanop) then
            print*,'NCHANOP, NCHOPT',nchanop,nchopt
            stop 'Increase NCHANOP'
         endif
         do nch = 1, nchtop + 1
            npkb(nch) = nch
         enddo 

         if (myid.le.0.and.lg.le.max(ldw,0).and.ipar.eq.0) then
            do k=1, max(1, npk(2) - npk(1) - nbnd(lg))
               qgrid(k) = gk(k,1)
            end do
            print*,' K     Q(K)    integration test'
            a=20.0/rmesh(meshr,1)
            l=lg
            jstart=1
            ione = 1
            izr = 0
            l1 = lg
            l2 = lg
            lm = 1
            T1 = GAMX(IONE,2*LM)/(GAMX(IONE,L1-L2+LM+1)
     >         *GAMX(IONE,L2-L1+LM+1))
            T2 = GAMX(IZR,L1+L2-LM+2)/(GAMX(IZR,L1+L2+LM+2)*(1.0D1**LM))
            do k=1, max(1, npk(2) - npk(1) - nbnd(lg))
               b=qgrid(k)
C  return SIN(X) in CHI(J)
               eta = zeff / (b + 1e-10)
c$$$               if (dbyexists) eta = 0.0

               call regular(l,b*b,eta,vasymp,cntfug(1,l),-1,rmesh,meshr,
c$$$               call regular(l,b*b,eta,ui,cntfug(1,l),ldw,rmesh,meshr,
     >            jdouble,id,regcut,expcut,chi,jstart,jstop,phase,sigc)
c$$$  if (k.le.9) then
c$$$                  do j = jstart, jstop
c$$$                     write(60+k,*) rmesh(j,1), chi(j)
c$$$                  enddo
c$$$                  close(60+k)
c$$$               endif 
               if (l.eq.0) then
                  tmp=b * b * rmesh(1,2) / 3.0
               else
                  tmp = 0.0
               endif 
               do i=jstart,meshr
                  r=rmesh(i,1)
                  tmp=tmp+chi(i)*chi(i)*rmesh(i,3)/r/r
               end do
               exact=b*pi/2.0/(2.0 * l + 1.0)
               cc = cos(l * pi / 2.0)
               cs = sin(l * pi / 2.0)
               rp = real(phase)
               ap = aimag(phase)
               x = atan2(ap * cc - rp * cs, rp * cc + ap * cs)
               tail = 0.0
               tailc = 0.0
               if (itail.ne.0) then
                  maxi = meshr
                  r1 = rmesh(maxi,1)
                  r2 = rmesh(maxi - 1,1)
                  rlam = log(rpow2(maxi,1)/rpow2(maxi-1,1))/(r2-r1)
                  a = rpow2(maxi,1) * exp(rlam*r1)
c$$$                  do i = 1, 100
c$$$                     r3 = rmesh(maxi-10,1) + i * rmesh(meshr,2)
c$$$                     print*,r3,1.0/r3**2,a*exp(-rlam*r3)
c$$$                  enddo 
                  tail2 = calctail2(a,rlam,r,b,x,b,x)
c$$$                  print*,tail2,a*exp(-rlam*r) * ( 0.5/rlam +
c$$$     >               (b*sin(2.0*b*r)-rlam*0.5*cos(2.0*r*b))/
c$$$     >               (rlam**2+4.0*b**2))
c$$$                  tailc = calctail(itail,r,b,x,b,x)
c$$$                  tail = ffgg(l,dble(b),l,dble(b),dble(r),1,
c$$$     >               1d0,0d0,1d0,0d0)
               endif
               rk1 = b
               sum2 = ((RK1*0.5D0)**LM)*T1*T2*pi/2.0
               tsum = tmp + tailc
               tmp = tmp + tail
c$$$               print*,(sin(b*rmesh(i,1)+x)/chi(i),i=meshr-3,meshr)
c$$$               print*,((real(phase) * sin(b*rmesh(i,1)-l*pi/2.0) + 
c$$$     >            imag(phase) * cos(b*rmesh(i,1)-l*pi/2.0))/
c$$$     >            chi(i),i=meshr-3,meshr)
               if (qgrid(k).gt.qcut)
     >            print*,'Warning QGRID point is > qcut'
               delta = real(log(phase)/(0.0,1.0))
               if (qgrid(k).gt.0.0)
     >            print'(i3,f10.4,7f12.5)',k,qgrid(k),tmp/exact,
     >            tsum/exact,sum2/exact,phase,delta,
     >            sin(delta)**2/qgrid(k)**2*4.0*pi
     >            
            end do
         end if 
         call update(6)
         
C  Define 1st Born matrix elements for subtraction in the cross program
         do ns = 0, 1
            do nchi = 1, nchtop
               do nchf = 1, nchtop
                  vdon(nchf,nchi,ns) = 0.0
               enddo
            enddo 
         enddo
         
         rmv = float((npk(nchtop+1)-1)) * (npk(nchtop+1)) ! * 4
         
         rmchi =  float(meshr) * (npk(nchtop+1)-1) ! * 4
         me2e = 0 * nchmaxe2e * (npk(nchtop+1)-1)
         mopt = 8 * kmax * kmax * (nchanop + 1) * nchanop
         
c$$$         mstep = 2 ** 20
c$$$         mtest = mstep
c$$$         ptrvmat = 1
c$$$         do while (ptrvmat.ne.0.and.mtest.lt.mv+mchi+me2e*2)
c$$$            mtest = mtest + mstep
c$$$            call memalloc(ptrvmat,mtest)
c$$$            call memfree(ptrvmat)
c$$$         enddo
         
c$$$         call memalloc(ptrvmat,mv*4)
c$$$         if (ptrvmat.eq.0) stop 'Not enough memory for VMAT'

c$$$         if (packed) then
c$$$            mv = mv * (nsmax+1) / 2
c$$$            print'(
c$$$     >         ''Memory (Mb) requested:'',i5,
c$$$     >         '' = VMATP + CHI:'', 2i5)', 
c$$$     >         (mv+mchi)/250000, mv/250000, mchi/250000
c$$$            deallocate(vmatp)
c$$$            allocate (vmatp((npk(nchtop+1)-1)*npk(nchtop+1)/2,0:nsmax))
c$$$            if (.not.allocated(vmatp))
c$$$     >         stop 'vmatp could not be allocated'
c$$$            vmatp(:,:) = 0.0
c$$$            allocate(vmat(1,1))
c$$$         else
         if (itail.lt.0) then
            allocate(chil(meshr,npkb(nchtop+1)-1,2))
            allocate(minchil(npkb(nchtop+1)-1,2))
            minchil(1:npkb(nchtop+1)-1,1:2) = min(meshr+1,maxr)
          print*,'Allocated CHIL(meshr,npkb,2):',meshr,npkb(nchtop+1)-1
         else
            allocate(chil(meshr,npkb(nchtop+1)-1,1))
            allocate(minchil(npkb(nchtop+1)-1,1))
          print*,'Allocated CHIL(meshr,npkb,1):',meshr,npkb(nchtop+1)-1
         endif 
         if (.not.allocated(chil)) stop 'chil could not be allocated'
C  Have set zasym to 0.0 and LDW to -1 for input to makechil so that
C  plane waves were always generated for Born subtraction
         zero=0.0
         m1=-1
         mnnbtop=abs(nnbtop)
         valuesin = valuesout
         call makechil(lg,gk,wk,qcut,zero,vdcore_pr,npot,ui,m1,dwpot,
     >      npkb,minchil,chil,phasel,nchtop,etot,nbnd0,abnd,npsbnd,
     >      albnd,sigma,mnnbtop,pos,lnch)
         call date_and_time(date,time,zone,valuesout)
         print '(/,i2,": nodeid exited first MAKECHIL at: ",a10,
     >   ", diff (secs):",i5)',nodeid,time,
     >   idiff(valuesin,valuesout)
         do nchi = 1, nchtop
            natomps(nchi) = 0
            do nchf = nchi, nchtop
c$$$            if (pos(nchf).neqv.pos(nchi)) natomps(nchi)=natomps(nchi)!+1
c$$$     >           +(2*lnch(nchi,1))*(2*lnch(nchi,2))
c$$$     >           +2**lnch(nchi,1)
               const = 1.0 
               if (pos(nchf).eqv.pos(nchi)) const = 0.0
               natomps(nchi)=natomps(nchi) + const * !+ 1
!     >              2.0**max(1,lnch(nchf,1))*2.0**max(1,lnch(nchi,1))*
     >              1.6**lnch(nchf,1)*1.6**lnch(nchi,1)*
     >              (npk(nchi+1)-npk(nchi))*
     >              (npk(nchf+1)-npk(nchf))/
     >              (npk(2)-npk(1))**2
c$$$     >           /1.15**abs(lnch(nchi,2)-lg)/1.15**abs(lnch(nchf,2)-lg)
c$$$               endif
            enddo
         enddo
         if (nodeid.eq.1.and.lptop.ge.0) then
            ntot = 0
            do nch = 1, nchtop
               ntot = ntot + natomps(nch)
               print*,'nch,natomps(nch),ntot:',nch,natomps(nch),ntot
            enddo
         endif 
         print'(
     >     ''Memory (Mb) requested:'',i7,
     >     '' = VMAT + CHI:'', 2i6)', 
     >     nint((rmv+rmchi)*nbytes/1e6), nint(rmv*nbytes/1e6),
     >     nint(rmchi*nbytes/1e6)

C Determine nchistart and nchistop for each node
            nchprst = (nchtop+1)*nchtop/2   !nch pairs total
            natompstot=0
            do nch= 1, nchtop
               natompstot = natompstot + natomps(nch)
            enddo
            natompspernode = natompstot/nodes
            sfactor=1.0            
            inquire(file='scale',exist=scaleexists)
            if (scaleexists) then
               open(42,file='scale')
               read(42,*) nn,(inc(n,0),n=1,nn)
               close(42)
               incsum = 0
               do n = 1, nodes - 1
                  incsum = incsum + inc(n,0)
               enddo
               inc(nodes,0) = - incsum ! just for output to time file purposes
            endif
            print"('nodeid, J, sfactor:',i4,i3,f5.2,2i5)",
     >         nodeid,lg,sfactor
            nchprst = nchprst / ((sfactor-1.0)/nodes+1.0)
            nchprspernode = nchprst / nodes
            nchistart(1) = 1
            ! made nodes.ge.1 below due to inefficiency with many Ps-states
!            if (natompstot.eq.0.or.nodes.ge.1.or.lg.gt.lstoppos) then ! no Ps states in the calculation or 1 node
c$$$            if (natompstot.eq.0.or.lg.gt.lstoppos) then ! no Ps states in the calculation or 1 node
            if (natompstot.eq.0.or.lg.gt.lstoppos) then ! no Ps states in the calculation or 1 node
               do nn = 1, nodes-1
                  print*,'nn,nchprspernode:',nn,nchprspernode
                  nchistop(nn) = nchistart(nn)
                  do while(nchprs(nchistart(nn),nchistop(nn),nchtop).le.
     >               nchprspernode) !was .le., but problems for many Ps states
                     nchistop(nn) = nchistop(nn) + 1
                  enddo
                  n1 = nchprs(nchistart(nn),nchistop(nn),nchtop)
                  n2 = nchprs(nchistart(nn),nchistop(nn)-1,nchtop)
                  if (n1-nchprspernode.gt.nchprspernode-n2)
     >               nchistop(nn) = nchistop(nn) - 1
                  nchprst = nchprst
     >               - nchprs(nchistart(nn),nchistop(nn),nchtop)
                  nchprspernode = nchprst / (nodes-nn)
                  nchistart(nn+1) = nchistop(nn) + 1
               enddo
               nchistop(nodes) = nchtop
            else
               natompsc = 0
               natompspernode = natompstot/nodes
               do nn = 1, nodes-1
                  natompsnode(nn)=natomps(nchistart(nn))
                  nchistop(nn) = nchistart(nn)
                  do while(natompsnode(nn).lt.natompspernode
     >                 -natomps(nchistop(nn)+1)/2
     >               .and.nchistop(nn).lt.nchtop)
                     nchistop(nn) = nchistop(nn) + 1
                     natompsnode(nn) =
     >                  natompsnode(nn) + natomps(nchistop(nn))
                  enddo
                  natompsc = natompsc + natompsnode(nn)
                  natompspernode = (natompstot-natompsc)/(nodes-nn)
                  nchistart(nn+1) = nchistop(nn) + 1
                  if (nodeid.eq.1)
     >               print"('nn,nchistart,natompsnode,natompspernode,',
     >               'natompsc,natompstot:',i3,i5,4i9)",
     >               nn,nchistart(nn), natompsnode(nn), natompspernode,
     >               natompsc, natompstot
               enddo
               nn = nodes
               natompsnode(nn) = natompstot - natompsc
               nchistop(nodes) = nchtop
               if (nodeid.eq.1)
     >            print"('nn,nchistart,natompsnode,natompspernode,',                                                                                                   
     >               'natompsc,natompstot:',i3,i5,4i9)",
     >               nn,nchistart(nn), natompsnode(nn),natompspernode,
     >               natompstot, natompstot
c$$$               if (natompsnode(nodes-1).lt.natompsnode(1)) then
c$$$                  nchistop(nodes-1) = nchistart(nodes)
c$$$                  nchistart(nodes) = nchistart(nodes) + 1
c$$$               endif
c$$$               if (nodes.gt.2) then ! put an extra nchi in the last two nodes
c$$$                  if (natompsnode(nodes-2).lt.natompsnode(1)) then
c$$$                     nchistop(nodes-2) = nchistop(nodes-2) + 1
c$$$                     nchistart(nodes-1) = nchistop(nodes-2) + 1
c$$$                     nchistop(nodes-1) = nchistop(nodes-1) + 1
c$$$                     nchistart(nodes) = nchistart(nodes) + 1
c$$$                  endif 
c$$$               endif
c$$$               do nn = 1, nodes
c$$$                  if (nodeid.eq.1) then
c$$$                     print"('nodeid,nchistart,nchistop,nchtop',4i5)",
c$$$     >               nn,nchistart(nn), nchistop(nn), nchtop
c$$$                     if (nchistart(nn).gt.nchtop) stop'nchistart>nchtop'
c$$$                  endif 
c$$$               enddo                  
            endif 
                  
            tscale(:) = 1.0
            tave(:)= 0.0
            tscale(nodes) = sfactor
            inquire(file='time'//ench,exist=exists)
            inquire(file='time_all',exist=timeexists)
            if (exists.or.timeexists) then !.and.lg.gt.mylstart) then
               if (exists) then
                  open(42,file='time'//ench)
               else
                  open(42,file='time_all')
               endif 
               ntime(:,:) = 0
               ip = 0
 10            read(42,*,end=20,err=20) lgp,n,ip,inc(n,ip),
     >            ntime(n,ip),nchistartold(n,ip),nchistopold(n,ip)
               lgold(ip) = lgp
               nodesold(ip) = n ! nodesold + 1
               if (ip.ne.ipar.or.lgp.ne.lg.or.n.ne.nodes) go to 10
c$$$               go to 10 !ensures read of last LG entry
 20            continue 
               close(42)
               print*,'Last LG read in time file:',lgold(ipar)
               if (nchistopold(nodesold(ip),ip).ne.nchtop) inc(:,:)=0  !reset to zero as sometimes non-zero for repeated low LG
               if (nodesold(ipar).eq.nodes.and.ntime(1,ipar).gt.0) then
                  ntimemin=10000000
                  ntimemax=0
                  ntimetot = 0
                  do n = 1, nodes
                     ntimetot = ntimetot + ntime(n,ipar)
                     if (ntime(n,ipar).lt.ntimemin) then
                        ntimemin = ntime(n,ipar)
                        nodemint = n
                     endif
                     if (ntime(n,ipar).gt.ntimemax) then
                        ntimemax = ntime(n,ipar)
                        nodemaxt = n
                     endif
                  enddo
                  tave(ipar) = float(ntimetot)/nodes
                  diffp = (ntimemax-ntimemin)/tave(ipar)
c$$$                  print*,'nodemint,nodemaxt,diffp:',nodemint,nodemaxt,
c$$$     >                 diffp
c$$$                  timeperi = float(ntimetot)/nchistopold(nodes,ipar)

            n = 1
            if (lgold(ipar).lt.10) then
               write(nodetfile,'(i3,"_",i1,"_",i1)') n,lgold(ipar),ipar
            elseif (lgold(ipar).lt.100) then
               write(nodetfile,'(i3,"_",i2,"_",i1)') n,lgold(ipar),ipar
            else
               write(nodetfile,'(i3,"_",i3,"_",i1)') n,lgold(ipar),ipar
            endif
            inquire(file=nodetfile,exist=exists)
            nchtimetot = 0
            do while (exists)
               open(42,file=nodetfile)
               nchistartold(n,ipar) = 10000000
               nchistopold(n,ipar) = 0
 13            read(42,'(i5,9x,i6)',end=14) nch,nchtime(nch)
               nchtimetot = nchtimetot + nchtime(nch)
               if (nch.lt.nchistartold(n,ipar)) nchistartold(n,ipar)=nch
               if (nch.gt.nchistopold(n,ipar)) nchistopold(n,ipar)=nch
               goto 13
 14            close(42)
               n = n + 1
               if (lgold(ipar).lt.10) then
                  write(nodetfile,'(i3,"_",i1,"_",i1)') 
     >                 n,lgold(ipar),ipar
               elseif (lgold(ipar).lt.100) then
                  write(nodetfile,'(i3,"_",i2,"_",i1)') 
     >                 n,lgold(ipar),ipar
               else
                  write(nodetfile,'(i3,"_",i3,"_",i1)') 
     >                 n,lgold(ipar),ipar
               endif
               inquire(file=nodetfile,exist=exists)
            enddo
            nodesprev = n - 1
            if (nodesprev.gt.1.and.
     >         nchistopold(nodesprev,ipar).eq.nchistop(nodes)) then
               tave(ipar) = float(nchtimetot)/float(nodesprev)
               if (myid.le.0) print'(
     >            "LGold,ipar,prev nodes,nchtimetot,ntmax,tave:",6i6)',
     >            lgold(ipar),ipar,nodesprev,nchtimetot,
     >              ntimemax,nint(tave(ipar))
               nt = 0
               n = 1
               nch = 1
               nistart = 1
               incsum = 0
               nodettot = 0 !tot node time before the current node
               nodetmax = 0 !max node time
               do n = 1, nodesprev - 1
                  nch = nistart
                  nodet = nchtime(nistart) !current node time
c$$$                  do while (nodet.lt.tave(ipar)) !.and.nch.lt.nchtopprev)
                  do while (nodet+nodettot.lt.n*tave(ipar).and.
     >                 nodet.lt.tave(ipar)*1.4)
                     nodetprev = nodet
                     nch = nch + 1
                     nodet = nodet + nchtime(nch)
                  enddo
c$$$                  if (nodet-tave(ipar).gt.tave(ipar)-nodetprev) then
                  if (nodet+nodettot-n*tave(ipar).gt.
     >                 n*tave(ipar)-nodetprev-nodettot
     >                 .or.nodet.gt.tave(ipar)*1.4) then
                     if (nch.gt.nistart) then
                        nch = nch - 1
                        nodet = nodetprev
                     endif
                  endif
                  if (nodet.gt.nodetmax) nodetmax = nodet
                  nodettot = nodettot + nodet
                  nistop = nch !max(nistart,nch)
                  incold(n,ipar) = nistop-nistart-
     >                 (nchistopold(n,ipar)-nchistartold(n,ipar))
                  if (nchistopold(nodesprev,ipar).eq.nchistop(nodes))
     >                 incold(n,ipar)=nistop-nistart-
     >                 (nchistop(n)-nchistart(n))-inc(n,ipar)
                  nistart = nistop + 1 !nch
                  if (myid.le.0) print'(
     >            "LG,ipar,nodeid,inc,nistop,nodet,nodettot,ntave:",
     >                 8i6)', lg,ipar,n,incold(n,ipar),
     >                 nistop,nodet,nodettot,nint(n*tave(ipar))
                  incsum = incsum + incold(n,ipar)
c$$$                  print"('node,nchistart,nchistop,time:',3i6,f6.1)", 
c$$$     >                 n,nchistart(n)inc(n,ipar),nodet
               enddo
               nodet = nchtimetot - nodettot
               if (nodet.gt.nodetmax) nodetmax = nodet
               neff = nint(100.0*nchtimetot/nodetmax/nodesprev)
               n = nodesprev
               nistop = nchistopold(n,ipar)
               if (myid.le.0)
     >            print'("LG,ipar,nodeid,inc,nistop,nodet:",
     >            15x,7i6,"%")',LG,ipar,n,-incsum,nistop,nodet,neff
            endif

                  incsum = 0
                  do n = 1, nodes - 1
c$$$                     timeperi = max(1.0 , float(ntime(n,ipar))/
c$$$     >                  (nchistopold(n,ipar)-nchistartold(n,ipar)+1))
! c$$$                     timeperi = 0.5 * ntime(n,ipar)/
! c$$$     >                  (nchistopold(n,ipar)-nchistartold(n,ipar)+1) +
! c$$$     >                  0.5 * ntime(n+1,ipar)/
! c$$$     >                  (nchistopold(n+1,ipar)-nchistartold(n+1,ipar)+1)
c$$$                     if (lptop.ge.0) then
c$$$                        if (diffp.gt.0.1) then
c$$$                           if (n.eq.nodemint) then
c$$$                              incstep = 1
c$$$                           elseif (n.eq.nodemaxt) then
c$$$                              incstep = -1
c$$$                           else
c$$$                              incstep = 0
c$$$                           endif
c$$$                        endif
c$$$                        incstep = 0
c$$$                        inc(n,ipar) = inc(n,ipar) + incold(n,ipar)
c$$$                     else
c$$$                        ni = nchistopold(n,ipar)-nchistartold(n,ipar)+1
c$$$                        timeperi = max(float(ntime(n,ipar))/ni, 1.0)
c$$$                       incstep=nint((tave(ipar)-ntime(n,ipar))/timeperi)
c$$$                        if (incstep.gt.ni/2) incstep = ni/2
c$$$                        if (incstep.lt.-ni/2) incstep = -ni/2
c$$$c$$$                     if ((tave-ntime(n,ipar))*(tave-ntime(n+1,ipar))
c$$$c$$$     >                  .lt.0.0) then
c$$$c$$$                        if (incstep.gt.1) incstep = 1
c$$$c$$$                        if (incstep.lt.-1) incstep = -1
c$$$c$$$                     else 
c$$$c$$$                        if (incstep.gt.ni) incstep = ni
c$$$c$$$                        if (incstep.lt.-ni/2) incstep = -ni/2
c$$$c$$$                     endif 
c$$$                     endif
c$$$                     inc(n,ipar) = inc(n,ipar) + incstep
                     inc(n,ipar) = inc(n,ipar) + incold(n,ipar)
                     incsum = incsum + inc(n,ipar)
c$$$                     print"('node,incstep,inc,timeperi:',3i6,f6.1)", 
c$$$     >                  n,incstep,inc(n,ipar),timeperi
                  enddo 
                  inc(nodes,ipar) = - incsum
c$$$                  tave = ntimetot / nodes
c$$$                  nchtot = 0
c$$$                  do n = 1, nodes !- 1
c$$$                     nrange(n) = nint(tave/tscale(n) * nchtop/nchtopold)
c$$$c$$$                     nrange(n) = nchistop(n)-nchistart(n)+1
c$$$c$$$                     if (ntime(n).lt.tave*0.9) then
c$$$c$$$                        nrange(n) = nrange(n) + 1
c$$$c$$$                     elseif (ntime(n).gt.tave*1.1) then
c$$$c$$$                        nrange(n) = nrange(n) - 1
c$$$c$$$                     endif
c$$$                     nchtot = nchtot + nrange(n)
c$$$                     if(nodeid.eq.1)print*,'node,nrange,nchtot,nchtop:',
c$$$     >                    n,nrange(n),nchtot,nchtop
c$$$                  enddo
c$$$c$$$                  if (nodes.gt.2) then
c$$$c$$$                     print*,'old nrange(nodes-1):',nrange(nodes-1)
c$$$c$$$                     nrange(nodes-1) = nrange(nodes-1)
c$$$c$$$     >                    - (nchtot-nchtop) / 2
c$$$c$$$                     print*,'new nrange(nodes-1):',nrange(nodes-1)
c$$$c$$$                  end if
c$$$c$$$                  if (nchtot.lt.nchtop) then
c$$$                     do n = 1, nodes - 1
c$$$                        nchistop(n) = nchistart(n) + nrange(n) - 1
c$$$                        nchistart(n+1) = nchistop(n) + 1
c$$$                     enddo
c$$$c$$$                  else
c$$$c$$$                     tscale(:)=1.0
c$$$c$$$                  endif 
               else
                  print*,
     >               'CAUTION: nodes<>nodesold',ipar,nodesold(ipar)
                  inc(:,:) = 0
               endif 
            endif
                  inct = 0
                  do n = 1, nodes - 1
                     inct = inct + inc(n,ipar)
c$$$                     print"('node,nchistopold:',2i6)",n,nchistop(n)
                     nchistop(n) = max(nchistop(n)+inct,nchistart(n))
                     nchistop(n) = min(nchistop(n),nchtop-nodes+n)
c$$$                     print"('node,nchistopnew:',2i6)",n,nchistop(n)
                     if (nchistart(n).gt.nchistop(n))
     >                  stop 'nchistart > nchistop'
                     nchistart(n+1) = nchistop(n) + 1
                  enddo


#ifdef _single
#define nbytes 4
#elif defined _double
#define nbytes 8
#endif

C Allocate the node-dependent VMAT arrays
            ni = npk(nchistart(nodeid))
            nf = npk(nchistop(nodeid)+1)-1
            nd = npk(nchtop+1)-1
            mb = 1000000
            mv01=0
            mv0=0
            mv1=0
            mchi = nint(1.0/mb*nd * meshr * nbytes)
            if (itail.ne.0) mchi = mchi*2
            npernode = (nf-ni+1)*(nf-ni+2)/2+(nf-ni+1)*(nd-nf)
c$$$            if (nodes.gt.nchtop) stop 'nodes > nchtop'	
c$$$            if (nodes.eq.1.and.scalapack)
c$$$     >         stop 'require at least 2 nodes with scalapack'
            if (nodeid.gt.1.or.scalapack) then 
               allocate(vmat01(ni:nf,ni:nf+1))
               mv01 = nint(1.0/mb*(nf-ni+1)*(nf+1-ni+1)*nbytes)
               mv0 = 0
               mv1 = 0
               if (.not.allocated(vmat01)) stop 'vmat01 not allocated'
c$$$               vmat01(:,:) = 0.0
               print*,'nodeid,ni,nf:',nodeid,ni,nf
               vmat01(ni:nf,ni:nf+1) = 0.0
               if (nf.lt.nd) then
                  allocate(vmat0(nf+1:nd,ni:nf))
                  mv0 = nint(1.0/mb*(nd-nf)*(nf-ni+1)*nbytes)
                  if (.not.allocated(vmat0)) stop 'vmat0 not allocated'
                  vmat0(:,:) = 0.0
                  if (nsmax.eq.1) then
                     allocate(vmat1(ni:nf,nf+1+1:nd+1))
                     mv1 = nint(1.0/mb*(nf-ni+1)*(nd-nf)*nbytes)
                     if (.not.allocated(vmat1))
     >                  stop 'vmat1 not allocated'
c$$$                     vmat1(:,:) = 0.0
                     vmat1(ni:nf,nf+1+1:nd+1) = 0.0
                  else
                     allocate(vmat1(1,1)) !avoid runtime errors
                  endif 
               else
                  allocate(vmat0(1,1),vmat1(1,1)) !avoid runtime errors
               endif
            else
               allocate(vmat01(1,2),vmat0(1,1),vmat1(1,1)) !as above
            endif 
            print"('nodeid,nchistart,nchistop,nchprs,mchi,mv01,mv0,mv1:'
     >         ,i4,':',2i5,i7,5(i7,'Mb'),'(tot)')",
     >         nodeid,nchistart(nodeid),
     >         nchistop(nodeid),nchprs(nchistart(nodeid),
     >         nchistop(nodeid),nchtop),mchi,mv01,mv0,mv1,
     >         mchi+mv01+mv0+mv1

            if (scalapack.or.nodeid.gt.1) then
               allocate(vmat(nchtop,nchtop+1))
            else
               allocate(vmat(npk(nchtop+1)-1,npk(nchtop+1)))
            endif 
            if (.not.allocated(vmat)) stop 'vmat could not be allocated'
            vmat(:,:) = 0.0
         
c$$$         call memalloc(ptrchi,mchi*4)
c$$$         if (ptrchi.eq.0) stop 'Not enough memory for CHI'

         if (isecond.ge.0) then
            call memalloc(ptrvopt,mvopt)
            if (ptrvopt.eq.0) stop 'Not enough memory for VMATOP'
         elseif (isecond.lt.-3) then
            second = lg .ge. -isecond
            if (second)print*,
     >         'Will be using the second order approximation'
         endif 

         call memalloc(ptre2ed,1)
         call memalloc(ptre2ee,1)
         call memalloc(ptrvopt,1)
         if (me2e.ne.0) then
            call memalloc(ptre2ed,me2e)
            if (ptre2ed.eq.0) stop 'Not enough memory for VE2ED'
            call memalloc(ptre2ee,me2e)
            if (ptre2ee.eq.0) stop 'Not enough memory for VE2EE'
            call initve2e(ve2ed,npkb(nchtop+1)-1,nchmaxe2e)
            call initve2e(ve2ee,npkb(nchtop+1)-1,nchmaxe2e)
         endif
         call update(6)
         td = 0.0
         te1 = 0.0
         te2 = 0.0
         t2nd = 0.0
                     
         call clock(s1)
         valuesin = valuesout
         if (nabot(labot).gt.1) call core(0,nznuc,lg,etot,chil,
     >      minchil,nchtop,uplane,-1,vdcore_pr,minvdc,maxvdc,npkb,
     >        vdon,vmat,vmatp)
         
         if (hlike) then
            ne2e1 = 0
            nchmaxe2e1 = 0
            nchtope2e1 = 0            

            call first(1,0,second,nold,etot,lg,gk,npkb,chil,minchil,
     >         uplane,-1,dwpot,phasel,itail,nznuc,nchtop,nchtope2e1,
     >         qcut,vdon,vmat,theta,vdcore_pr,minvdc,maxvdc,lfast,lslow,
     >         slowery,td,te1,te2,ve2ed,ve2ee,dphasee2e,ephasee2e,ne2e1,
     >         nchmaxe2e1,vmatp,nsmax,
     >         nchistart,nchistop,nodeid,scalapack,
     >         vmat01,vmat0,vmat1,ni,nf,nd,nodes,myid,natomps,lnch)
         else
            call scattering(myid,0,theta,nold,etot,lg,gk,enionry,npkb,
     >         chil,minchil,vdcore_pr,dwpot,nchtop,nmaxhe,namax,
     >         nze,td,te1,te2,t2nd,vdon,vmat,nsmax,itail,phasel)
         end if
         call clock(s2)
c$$$         print*,'Time to make the partial Born matrix elements:',s2-s1
         call date_and_time(date,time,zone,valuesout)
         print '(/,i2,": nodeid exited first VMAT routines at: ",a10,
     >      ", diff (secs):",i5)',nodeid,time,idiff(valuesin,valuesout)
         call update(6)
         
C  Define the plane (or distorted by ui + vdcore) waves
c$$$         if (dbyexists) then
c$$$            zas = 0.0
c$$$         else
c$$$            zas = zasym
c$$$         endif 
!            print*, 'calling makechil 3'
         deallocate(chil,stat=istat)
         deallocate(minchil,stat=istat)
         if (itail.lt.0) then
            allocate(chil(meshr,npk(nchtop+1)-1,2))
            allocate(minchil(npk(nchtop+1)-1,2))
            minchil(1:npk(nchtop+1)-1,1:2) = min(meshr+1,maxr)
            print*,'Allocated CHIL(meshr,npk,2):',meshr,npk(nchtop+1)-1
         else
            allocate(chil(meshr,npk(nchtop+1)-1,1))
            allocate(minchil(npk(nchtop+1)-1,1))
            print*,'Allocated CHIL(meshr,npk,1):',meshr,npk(nchtop+1)-1
         endif 
         if (.not.allocated(chil)) stop 'chil could not be allocated'

         call date_and_time(date,time,zone,valuesin)
         print '(/,i2,": nodeid calling second MAKECHIL at: ",a10)',
     >      nodeid,time
         mnnbtop=abs(nnbtop)
         call makechil(lg,gk,wk,qcut,zasym,vdcore_pr,npot,ui,ldw,dwpot,
     >      npk,minchil,chil,phasel,nchtop,etot,nbnd,abnd,npsbnd,albnd,
     >        sigma,mnnbtop,pos,lnch)
         call date_and_time(date,time,zone,valuesout)
         print '(/,i2,": nodeid exited  second MAKECHIL at: ",a10,
     >   ", diff (secs):",i5)',nodeid,time,
     >   idiff(valuesin,valuesout)

         call clock(s3)
c$$$         print*,'Time to make the partial wave table:',s3-s2
         call update(6)
         
C  Calculate the D matrix for photoionization
         if (projectile.eq.'photon') then
            if(ntype.eq.0)  call xSATEL(nchtop,lg)
            if(ntype.eq.1)  call mSATEL(nchtop,lg)
            call clock(s1)
            if (lg.eq.1) call dMATR(gk,npk,minchil,chil,nchtop)
            call clock(s2)
            print*,'time for dMATR:',s2-s1
!            if (mylstop.gt.1) then
C (gamma,3e) by double shake-off
!            if(lg.eq.0) call mMATR(gk,npk,minchil,chil,nchtop)

C (e,3e) agenda
               call clock(s1)
               call bMATR (lg,gk,npk,minchil,chil,nchtop)
               if(lg.eq.0 .or. lg.eq.2)
     >            call b2MATR(lg,gk,npk,minchil,chil,nchtop)
!               call qMATR(lg,gk,npk,minchil,chil,nchtop)
               call clock(s2)
               print*,'time for bMATR:',s2-s1
            endif 
!         endif 
         
C  Initialise the first order V matrix.
         call update(6)
c$$$         if (packed) then
c$$$            vmatp(:,:) = 0.0
c$$$         else 
            call update(6)
            vmat(:,:) = 0.0
C The following code can be commented out to yield previously working results
C Note that it does not touch on-shell wk (needed in ScaLAPACK).
            if (scalapack) then
               do nch = nchistart(nodeid), nchistop(nodeid)
                  do k = npk(nch)+1,npk(nch+1)-1
                     vmat01(k,k) = !-1.0/real(wk(k))
     >                    -1.0/gf(k-npk(nch)+1,k-npk(nch)+1,nch)
                     vmat01(k,k+1) = vmat01(k,k)
                  enddo
               enddo
            else
               do nch = 1, nchtop
                  do k = npk(nch)+1,npk(nch+1)-1
                     vmat(k,k) = !-1.0/real(wk(k))
     >                    -1.0/gf(k-npk(nch)+1,k-npk(nch)+1,nch)
c$$$                     print*,'nch,k,vmat(k,k):',nch,k,vmat(k,k),
c$$$     >                    -1.0/gf(k-npk(nch)+1,k-npk(nch)+1,nch)
                     vmat(k,k+1) = vmat(k,k) !-1.0/real(wk(k))
                  enddo
               enddo
            endif
            do nch = 1, nchtop 
               do k = npk(nch)+1,npk(nch+1)-1
                  wk(k)=cmplx(1e30,0.0) !offshell only, on every node
               enddo
            enddo
C-----------------------------------------------------------------------------
            
c$$$  call initv(vmat,npk(nchtop+1)-1,nsmax)
c$$$         endif 

c$$$  if (ne2e.ne.0) then
c$$$            call initve2e(ve2ed,npk(nchtop+1)-1,nchmaxe2e)
c$$$            call initve2e(ve2ee,npk(nchtop+1)-1,nchmaxe2e)
c$$$         endif 
c
c     if speed is true then we must obtain the parts of the v matrix
c     that have already been worked out for the input energy. This
c     is done via interpolation in the intrpv routine.
c         
         call clock(v1)
         if (speed) then
            stop 'should not be here'
            call intrpv(lg,nsmax,vkeep,gkeep,nbnd,npkeep,nchkeep,
     >         npkeep(nchkeep+1)-1,npk,nchtop,npk(nchtop+1)-1,gk,
     >         vmat,pint)
         endif
         tc=0.0
         if (nabot(labot).gt.1) then
            call clock(s1)
            call core(ifirst,nznuc,lg,etot,chil,minchil,nchtop,uplane,
     >         ldw,vdcore_pr,minvdc,maxvdc,npk,vdondum,vmat,vmatp)
            call clock(s2)
            tc = s2 - s1
         endif
c
c     if speed is true then ispeed needs to be set to 3 in order
c     for the first subroutine to only calculate the elements
c     that have not already been worked out via the interpolation
c     routine
c
         if (speed) ispeed=3

         call date_and_time(date,time,zone,valuesin)
         print '(/,i4,": nodeid entering VMAT routines at: ",a10)',
     >      nodeid,time
         if (hlike) then
            call first(ispeed,ifirst,second,nold,etot,lg,gk,npk,chil,
     >         minchil,
     >         ui,ldw,dwpot,phasel,itail,nznuc,nchtop,nchtope2e,qcut,
     >         vdondum,vmat,theta,vdcore_pr,minvdc,maxvdc,lfast,lslow,
     >         slowery,td,te1,te2,ve2ed,ve2ee,dphasee2e,ephasee2e,ne2e0,
     >         nchmaxe2e,vmatp,nsmax,
     >         nchistart,nchistop,nodeid,scalapack,
     >         vmat01,vmat0,vmat1,ni,nf,nd,nodes,myid,natomps,lnch)
            call clock(s1)
            if (isecond.ge.0) then
               stop 'Have not coded for LDW, NPK, or NQM'
               call optpot(ntstart,npsbnd,lnabtop,nnbtop,lstart,lg,nqm,
     >            nchopt,zasym,nznuc,ui,minchil,chil,lttop,nptop,lptop,
     >            etot,isecond,small,nchns,anorm,vdcore_pr,corep,r0,
     >            vmatop)
            end if 
            call clock(s2)
            t2nd = s2 - s1
         else
            call scattering(myid,ifirst,theta,nold,etot,lg,gk,enionry,
     >         npk,chil,minchil,vdcore_pr,dwpot,nchtop,nmaxhe,namax,
     >         nze,td,te1,te2,t2nd,vdondum,vmat,nsmax,itail,phasel)
         end if 
         call date_and_time(date,time,zone,valuesout)
cDIR$ SUPPRESS
         print '(/,i4,": nodeid exiting  VMAT routines at: ",a10,
     >      ", J=",i2,", nchprs:",i7,", diff (secs):",i5)',nodeid,time,
     >      lg,nchprs(nchistart(nodeid),nchistop(nodeid),nchtop),
     >      idiff(valuesin,valuesout)
         call update(6)
         ntime(:,ipar) = 0
         ntime(nodeid,ipar) = idiff(valuesin,valuesout)
         valuesin = valuesout

c     
c     resetting ispeed =1
c
         ispeed=1
         call clock(v2)
         print*,'time for VMAT:',v2-v1
         inquire(file='reconstruct_psi',exist=reconstruct_psi)
c$$$         reconstruct_psi = .false.
         if (.not.reconstruct_psi) then
            deallocate(chil,stat=istat)
            deallocate(minchil,stat=istat)
            allocate(chil(1,1,2)) !stop errors when passing as argument
            print*,'Deallocated CHIL',istat ! comment out this statement to use CHIL in tmatcco.f
         endif 
         endif ! mod(myid,nomp).eq.0 
#ifdef _single
#define MY_MPI_REAL MPI_REAL
#define MY_MPI_COMPLEX MPI_COMPLEX
#elif defined _double
#define MY_MPI_REAL MPI_DOUBLE_PRECISION
#define MY_MPI_COMPLEX MPI_DOUBLE_COMPLEX
#endif
#define TAG0 1
#define TAG1 2
#define TAG01 3
         if (.not. allocated(vmat0)) allocate(vmat0(1,1)) 
         if (.not. allocated(vmat1)) allocate(vmat1(1,1))
         if (.not. allocated(vmat01)) allocate(vmat01(1,2))
         if (scalapack) then
C The following allows scalapack to run efficiently with threaded libraries
#ifndef LIBSCI
c$$$            call mkl_set_num_threads(1) ! revert for many tasks per node
            call mkl_set_num_threads(nomporig) 
#endif
            call sleepy_barrier(MPI_COMM_WORLD)
            if (mod(myid,nomp).eq.0) then
               call date_and_time(date,time,zone,valuesout)
               print '(/,i4,": nodeid exiting Sleepy_Barrier at: ",a10,
     >            ", diff (secs):",i5)',nodeid,
     >            time,idiff(valuesin,valuesout)
               valuesin = valuesout
            endif
            call factor_two(ntasks,nrows,ncols)
            call blacs_get(0,0,nblacs_context)
            call blacs_gridinit(nblacs_context,'R',nrows,ncols)
            call mpi_bcast(nchistart,nodes,MPI_INTEGER,0,
     >         MPI_COMM_WORLD,ierr)
            call mpi_bcast(nchistop,nodes,MPI_INTEGER,0,
     >         MPI_COMM_WORLD,ierr)
            call mpi_bcast(nchtop,1,MPI_INTEGER,0,
     >         MPI_COMM_WORLD,ierr)
            call mpi_bcast(nsmax,1,MPI_INTEGER,0,
     >         MPI_COMM_WORLD,ierr)
            call mpi_bcast(npk,nchtop+1,MPI_INTEGER,0,
     >         MPI_COMM_WORLD,ierr)
            nd = npk(nchtop+1)-1
            if (mod(myid,nomp).eq.0) then
               ni = npk(nchistart(myid/nomp+1))
               nf = npk(nchistop(myid/nomp+1)+1)-1
            else
               ni=0
               nf=0
            end if
            call mpi_bcast(wk,nd,MY_MPI_COMPLEX,0,
     >         MPI_COMM_WORLD,ierr)

            if (myid.eq.0) then
               allocate(bb(nd,nchtop,2,0:nsmax))
            else
               if (.not.allocated(bb)) allocate(bb(1,1,2,0:1)) !not used
            endif
            do ns = 0, nsmax
               call redistributeAndSolve(ni,nf,nd
     >             ,nblacs_context,ns,nsmax
     >             ,nodes,nchtop+1,npk
     >             ,nchtop
     >             ,wk,nd,bb(1,1,1,ns))
            enddo
            call blacs_gridexit(nblacs_context)
         else ! .not.scalapack
C Reconstruct full VMAT on the first node to use with LAPACK
         if (myid.eq.0) then
            nd = npk(nchtop+1)-1
            allocate(bb(nd,nchtop,2,0:nsmax))
            bb(:,:,:,:) = 0.0
            receiveHandle(:,:)=MPI_REQUEST_NULL
            do nn = 2, nodes
               ni = npk(nchistart(nn))
               nf = npk(nchistop(nn)+1)-1
               call mpi_type_vector(nf-ni+2,nf-ni+1,nd,
     >         MY_MPI_REAL,
     >            receiveType(nn,3), ierr) ! Defines the strided memory access to receive vmat01
                  call mpi_type_commit(receiveType(nn,3),ierr)
                  call mpi_irecv(vmat(ni,ni), 1, receiveType(nn,3),
     >               (nn-1)*nomp, TAG01, MPI_COMM_WORLD,
     >               receiveHandle(nn,3),ierr)
               
               if (nf.lt.nd) then
                  call mpi_type_vector(nf-ni+1,nd-nf,nd,
     >            MY_MPI_REAL,
     >               receiveType(nn,1), ierr) ! Defines the strided memory access to receive vmat0
                  call mpi_type_commit(receiveType(nn,1),ierr)
                  call mpi_irecv(vmat(nf+1,ni), 1, receiveType(nn,1),
     >               (nn-1)*nomp, TAG0, MPI_COMM_WORLD,
     >               receiveHandle(nn,1),ierr)
                  if (nsmax.eq.1) then ! second spin
                     call mpi_type_vector(nd-nf,nf-ni+1,nd,
     >               MY_MPI_REAL,
     >                  receiveType(nn,2), ierr) ! Defines the strided memory access to receive vmat1
                     call mpi_type_commit(receiveType(nn,2),ierr)
                     call mpi_irecv(vmat(ni,nf+1+1), 1, 
     >                  receiveType(nn,2),(nn-1)*nomp, TAG1, 
     >                  MPI_COMM_WORLD,receiveHandle(nn,2),ierr)
                  endif 
               else
                  receiveHandle(nn,1) = MPI_REQUEST_NULL
               endif
            enddo
            do i = 1, 3
               call mpi_waitall(nodes,receiveHandle(:,i),
     >            statusesR(:,:),ierr)
            end do
            do nn=2,nodes
               ni = npk(nchistart(nn))
               nf = npk(nchistop(nn)+1)-1
               call mpi_type_free(receiveType(nn,3),ierr)
               if (nf.lt.nd) then
                  call mpi_type_free(receiveType(nn,1),ierr)
                  if (nsmax.eq.1)
     >               call mpi_type_free(receiveType(nn,2),ierr)
               end if
            end do
         else if (mod(myid,nomp).eq.0) then
            nd = npk(nchtop+1)-1
            ni = npk(nchistart(nodeid))
            nf = npk(nchistop(nodeid)+1)-1
            sendHandle(:)=MPI_REQUEST_NULL

c$$$            do i=ni,nf
c$$$               do j = nf+1+1,nd+1
c$$$                  vmat1(i,j)=float(i)+float(j)/1000.0
c$$$               enddo
c$$$            enddo 

            call mpi_isend(vmat01,(nf-ni+1)*(nf-ni+2),MY_MPI_REAL,0,
     >         TAG01,MPI_COMM_WORLD,sendHandle(3),ierr)
               
            if (nf.lt.nd) then
               call mpi_isend(vmat0,(nf-ni+1)*(nd-nf),MY_MPI_REAL,0,
     >            TAG0,MPI_COMM_WORLD,sendHandle(1),ierr)
               if (nsmax.eq.1)
     >            call mpi_isend(vmat1,(nf-ni+1)*(nd-nf),MY_MPI_REAL,
     >            0,TAG1,MPI_COMM_WORLD,sendHandle(2),ierr)
            endif 
            call mpi_waitall(3,sendHandle,statusesS(:,:))
         endif
         endif ! if scalapack
         if (allocated(vmat01)) deallocate(vmat01)
         if (allocated(vmat0)) deallocate(vmat0)
         if (allocated(vmat1)) deallocate(vmat1)

         nntime(:) = 0
! The following is causing problems on Frontera
         CALL MPI_REDUCE(ntime(1,ipar), nntime, nodes, MPI_INTEGER, 
     >      MPI_SUM, 0, MPI_COMM_WORLD, ierr)

C   Solve Ax=b inside the main routine; below will be replaced with scalapack
         if (myid.eq.0) then
            call date_and_time(date,time,zone,valuesout)
            if (scaleexists) then
               open(42,file='time'//ench,position='append')
               ntmax = 0
               timetot = 0.0
               do n=1,nodes
                  timetot = timetot + nntime(n)
                  if (nntime(n).gt.ntmax) ntmax = nntime(n)
               enddo
               neff = 0
               if (ntmax.gt.0) neff = 100.0*timetot/nodes/ntmax
               do n=1,nodes-1
                  naps = 0
                  do nch = nchistart(n), nchistop(n)
                     naps = naps + natomps(nch)
                  enddo
                  timeperi = float(nntime(n))/
     >               (nchistart(n+1)-nchistart(n))
c$$$                  incw = 0
c$$$                  if (n.eq.nodetimemax) inc(n,ipar) = inc(n,ipar)-inct
c$$$                  if (n.eq.nodetimemin) inc(n,ipar) = inc(n,ipar)+inct
c$$$                  write(42,'(4i4,3i7,f8.1,3i8," ave time of prev LG")')
                  write(42,'(4i4,3i7,f8.1,2i8)') 
     >               lg,n,ipar,inc(n,ipar),
     >               nntime(n),nchistart(n),nchistop(n),timeperi,
     >               nchprs(nchistart(n),nchistop(n),nchtop),
     >               naps
!                  write(42,'(4i4,7i7," ave time of prev LG")') 
!     >               lg,n,ipar,inc(n,ipar),
!     >               nntime(n),nchistart(n),nchistop(n),ntimeperi,
!     >               nchprs(nchistart(n),nchistop(n),nchtop),
!     >               naps, nint(tave)
               enddo
               n = nodes
               naps = 0
               do nch = nchistart(n), nchistop(n)
                  naps = naps + natomps(nch)
               enddo
               timeperi = float(nntime(n))/
     >            (nchistop(n)-nchistart(n)+1)
c$$$                  incw = 0
c$$$                  if (n.eq.nodetimemax) inc(n,ipar) = inc(n,ipar)-inct
c$$$                  if (n.eq.nodetimemin) inc(n,ipar) = inc(n,ipar)+inct
               write(42,'(4i4,3i7,f8.1,5i8,2i4,"% LG,node,ipar,inc,vt,",
     >            "i1,i2,tperi,nch,naps,NMAX,tave,mt,prev LG,eff",/)')
     >            lg,n,ipar,inc(n,ipar),nntime(n),nchistart(n),
     >            nchistop(n),timeperi,nchprs(nchistart(n),nchistop(n),
     >            nchtop),naps,npk(nchtop+1)-1,nint(timetot/nodes),
     >            idiff(valuesin,valuesout),lgold(ipar),neff
!               write(42,'(4i4,3i7,f8.1,3i8,2i4,"% LG,node,ipar,inc,vt,",
!     >            "i1,i2,tperi,nch,naps,mt,prev LG,eff",/)')
!     >            lg,n,ipar,inc(n,ipar),nntime(n),nchistart(n),
!     >            nchistop(n),timeperi,nchprs(nchistart(n),nchistop(n),
!     >            nchtop),naps,idiff(valuesin,valuesout),lgold(ipar),
!     >            neff
!               write(42,'(4i4,7i7,2i4,"% LG,node,ipar,inc,vt,",
!     >            "i1,i2,tperi,nch,naps,mt,LGold,eff")')
!     >            lg,n,ipar,inc(n,ipar),nntime(n),nchistart(n),
!     >            nchistop(n),ntimeperi,nchprs(nchistart(n),nchistop(n),
!     >            nchtop),naps,idiff(valuesin,valuesout),lgold(ipar),
!     >            neff
               close(42)
            endif 
 
c$$$         do ns = 0, nsmax
c$$$            allocate(AA(nd,nd))
c$$$            do nchi = 1, nchtop
c$$$               kf = 1
c$$$               call makev(vmat,kf,nchi,nchtop,nd,npk,wk,ns,
c$$$     >            bb(1,nchi,1,ns),vmatp,packed)
c$$$            enddo
c$$$            if (ns.eq.0) then
c$$$               do ki = 1, nd
c$$$                  do kf = ki, nd
c$$$                     AA(kf,ki) = - vmat(kf,ki) *
c$$$     >                  sqrt(abs(real(wk(kf))*real(wk(ki))))
c$$$                     AA(ki,kf) = AA(kf,ki)
c$$$                  enddo
c$$$               enddo
c$$$            else
c$$$               do ki = 1, nd
c$$$                  do kf = ki, nd         
c$$$                     AA(ki,kf) = - vmat(ki,kf+1) *
c$$$     >                  sqrt(abs(real(wk(kf))*real(wk(ki))))
c$$$                     AA(kf,ki) = AA(ki,kf)
c$$$                  enddo
c$$$               enddo 
c$$$            endif
c$$$            do nchi = 1, nchtop
c$$$               ki = npk(nchi)
c$$$               do nchf = 1, nchtop
c$$$                  kf = npk(nchf)
c$$$                  AA(kf,ki) = 0.0
c$$$               enddo
c$$$            enddo 
c$$$                  
c$$$C  Add the I matrix
c$$$            do kf = 1, nd
c$$$               s = (real(wk(kf))+1e-30) / abs(real(wk(kf))+1e-30)
c$$$               AA(kf,kf) = AA(kf,kf) + s
c$$$            end do 
c$$$            call date_and_time(date,time,zone,valuesin)
c$$$            print '("SGESV entered at: ",a10)',time
c$$$            allocate(ipvt(nd))
c$$$#ifdef _single
c$$$            call sgesv(nd,nchtop,AA,nd,ipvt,bb(1,1,1,ns),nd,info)
c$$$#elif defined _double
c$$$            call dgesv(nd,nchtop,AA,nd,ipvt,bb(1,1,1,ns),nd,info)
c$$$#endif
c$$$            deallocate(ipvt)
c$$$            deallocate(AA)
c$$$            call date_and_time(date,time,zone,valuesin)
c$$$            print '(".GESV exited at:  ",a10)',time
c$$$         enddo ! end of ns loop
            if (scalapack.and.reconstruct_psi) then
               allocate(kon(nchtop,nchtop))
               allocate(ckern(nchtop,nchtop))
               allocate(ct(nchtop,nd),cwork(nchtop),ct2(nd))
               allocate(psi12(meshr,meshr))
               
               do ns = 0, nsmax
                  kon(:,:) = 0.0
                  ckern(:,:) = (0.0,0.0)
                  do nchi = 1, nchtop
                     do nchf = 1, nchtop
                        if (gk(1,nchf).gt.0.0.and.gk(1,nchi).gt.0.0)then
                           kon(nchf,nchi) = (bb(npk(nchf),nchi,2,ns)+
     $                          dot_product(bb(1:nd,nchf,2,ns),
     $                          bb(1:nd,nchi,1,ns)))
c$$$     $                          /gk(1,nchf)/gk(1,nchi)
                           ckern(nchf,nchi) = cmplx(0.0,pi/gk(1,nchi))*
     $                          kon(nchf,nchi)
                        end if
                     end do
                     ckern(nchi,nchi) = ckern(nchi,nchi) + (1.0,0.0)
                  end do
                  do nchi = 1, nchtop
                     do nchf = 1, nchtop
                        ct(nchi,npk(nchf)) = kon(nchf,nchi) !* wk(npk(nchf)
                        nn = 1
c$$$                        print*,'gk, Kon:',gk(1,nchf),ct(nchi,npk(nchf))/
c$$$     $                       (gk(1,nchf)*gk(1,nchi))
                        do n = npk(nchf)+1, npk(nchf+1) - 1
                           nn = nn + 1
                           ct(nchi,n) = bb(n,nchi,1,ns) 
     $                          /(real(wk(n))+1e-20)!*gk(nn,nchf)*gk(1,nchi))
c$$$                           print*,'nchi,nn:',gk(nn,nchf),ct(nchi,n)/
c$$$     $                          (gk(nn,nchf)*gk(1,nchi))
                        end do
                     end do
                  end do
                  call matinv2(ckern,nchtop,nchtop,ct,nd,
     $                 cwork,0.0,0.0)
                  do nchi = 1, nchtop
                     ct2(1:nd) = ct(nchi,1:nd)*wk(1:nd)
                     do nchf = 1, nchtop
c$$$                        print*,'ns, nchi, nchf, On-shall T:',ns,nchi,
c$$$     $                    nchf,ct(nchf,npk(nchi))/gk(1,nchf)/gk(1,nchi),
c$$$     $                       (bb(npk(nchf),nchi,2,ns)+
c$$$     $                       dot_product(bb(1:nd,nchf,2,ns),
c$$$     $                       ct2(1:nd)))/gk(1,nchf)/gk(1,nchi)
                     end do
                  end do
                  nstep = 5
                  psi12(:,:) = (0.0,0.0)
                  if (nchtop.gt.30) stop 'increase ctemp array'
                  do nchi = 1, 1 !nchtop  ! 1  for psi12
                     ct2(1:nd) = ct(nchi,1:nd)*wk(1:nd)
                     do nch = 1, nchtop
                        call getchinfo (nch,nchp,lg,temp,maxpsi,enpsi,
     $                    la,na,lp)
                        chann(nch) = ch(na)//' '//ch(la)//' '//ch(lp)
                        ctemp(:,nch) = (0.0,0.0)
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
                        do n = npk(nch), npk(nch+1)-1
                           do i = 1, meshr
                              ctemp(i,nch) = ctemp(i,nch) +
     >                           chil(i,n,1)*ct2(n)/rmesh(i,3)
                           enddo 
                           do i = 1, meshr, nstep
                              do j = 1, maxpsi, nstep
                                 psi12(i,j) = psi12(i,j) + chil(i,n,1)/
     $                              rmesh(i,3)*temp(j)*ct2(n)
                              end do
                           end do
                        end do
C$OMP END PARALLEL DO
                     end do
                     print*,'Defn of psi12 complete'
                     call getchinfo(nchi,nchp,lg,temp,maxpsi,enpsi,la,
     $                 na,lp)
                     if (nznuc.eq.2.and.zasym.eq.1.0) then
                        write(enq,'(1p,e10.4,"_",3i1,"_",i1)')
     >                     energy-(4.0+enpsi)*ry,lg,ns,ipar,nchi
                        open(42,file=enq)
                        deta = nze * zasym / sqrt(ery)
                        write(42,'(" ",i6,i2,a11,a16,20a22)') meshr,
     >                     nchtop,chann(nchi),(chann(nch),nch=1,nchtop)
c$$$                     coulphase(deta,lg)
                        do i = 1, meshr
                           write(42,'(1p,42e11.3)') rmesh(i,1),
     >                        chil(i,npk(nchi),1)/rmesh(i,3),
     >                        (ctemp(i,nch),nch=1,nchtop)
                        enddo
                        close(42)
                     endif      !nznuc.eq.2...
                  enddo
                  nchi = 1
                  call getchinfo(nchi,nchp,lg,temp,maxpsi,enpsi,la,na,
     $                 lp)
                  open(42,file='psi_'//ch(ns))
                  write(42,'("# unsymmetrized psi(r1,r2) and phi(r1,r2)",
     >               " for S =",i2)') ns
                  write(42,'("#   r2        r1      Re(psi)   Im(psi)",
     >               "   Phi  ")')
c$$$                  do i = 1, meshr, nstep
c$$$                     do j = 1, maxpsi, nstep
c$$$                        psi12(i,j) = psi12(i,j) + chil(i,1,1)/
c$$$     $                       rmesh(i,3)*temp(j)
c$$$                     enddo
c$$$                  enddo
                  do i = 1, maxpsi, nstep
                     do j = 1, maxpsi, nstep
                        write(42,'(1p,5e10.2)') 
     $                     rmesh(j,1),rmesh(i,1),psi12(i,j),
     >                     chil(i,1,1)/rmesh(i,3)*temp(j)
c$$$     >                     (psi12(i,j)+(-1)**ns*psi12(j,i))
                     end do
                     write(42,*)
                  end do
                  close(42)
               end do           !ns
               deallocate(kon,ckern,ct,ct2,cwork,psi12)
            end if ! reconstruct_psi
            call clock(s2)
         call date_and_time(date,time,zone,valuesin)
         print '("SOLVET entered at:",a10)',time
 
         call solvet(ifirst,bb,vmat,gk,wk,weightk,nchtop,nqm,noprint,
     >      nopen,etot,lg,vdon,phasel,nent,instate,
     >      isecond,nunit,sigma,vmatop,
     >      nchopt,ovlpn,phaseq,ve2ed,ve2ee,dphasee2e,ephasee2e,
     >      nchtope2e,lfast,lslow,slowery,nznuc,zasym,npk,ne2e0,
     >      nchanmax,projectile,target,uba,nsmax,canstop(ipar),
     >      second,BornICS,tfile,csfile,e2efile,npk(nchtop+1)-1,theta,
     >      abs(nnbtop),ovlpnn,vmatp,scalapack,chil)
         call date_and_time(date,time,zone,valuesout)
         print '("SOLVET  exited at: ",a10,", diff (secs):",i5)',
     >      time, idiff(valuesin,valuesout)

         if (allocated(bb)) deallocate(bb)

         call clock(s3)
         tm = s3 - s2
c$$$         call gettime(utime,stime)
c$$$         call memfree(ptrvmat)

c$$$         if (packed) then
c$$$            deallocate(vmatp)
c$$$         else 
c$$$            deallocate(vmat)
c$$$            print*,'VMAT deallocated for myid, LG:',myid,lg
c$$$         endif
         
c$$$         call memfree(ptre2ed)
c$$$         call memfree(ptre2ee)
         print'(''J ='',i3,'' IP ='',i2,'' Tc:'',i5,'' Td:'',i6,
     >      '' Te1:'',i6,'' Te2:'',i5,'' T2nd:'',i5,'' Tm:'',i7)',lg,
     >      ipar,nint(tc),nint(td),nint(te1),nint(te2),nint(t2nd),
     >           nint(tm)
         if (nunder.gt.0) then
            print *,'Number of underflows was:',nunder
            nunder=0
         endif 
         if (nover.gt.0) then
            print *,'Number of overflows was:',nover
            nover=0
         endif 
         call update(6)

c
c     end of energy loop
c
         print*,'Energy completed - energy was ',energy
         call clock(timeenergy2)
         print*,'Time for energy loop was: ',timeenergy2-timeenergy1
      endif !  myid.eq.0
c$$$      print*,'LBOUNDs of VMAT:',lbound(vmat,1),lbound(vmat,2)
c$$$      print*,'UBOUNDs of VMAT:',ubound(vmat,1),ubound(vmat,2)
      
      if (mod(myid,nomp).eq.0) then
         deallocate(vmat,stat=istat)
         if (allocated(chil))deallocate(chil)
         if (allocated(minchil))deallocate(minchil)
c$$$ 775  continue
      else
c$$$         call sleepy_barrier(MPI_COMM_WORLD)! On Magnus this fails
      endif
      if (nodes.gt.1) call sleepy_barrier(MPI_COMM_WORLD)

c$$$      call sleepy_barrier(MPI_COMM_WORLD) ! This works perfectly on Magnus

      enddo                     !end of ipar loop
c$$$         if (canstop(0).and.canstop(abs(npar)).and.ntype.ge.0.and.
c$$$     >      ntasks.eq.-1) go to 780
      if (myid.eq.0) then
         call date_and_time(date,time,zone,valuesout)
         print '(/,"Partial wave LG: ",i4," completed at: ",a10,
     >   ", diff (secs):",i5)',lg,time,
     >        idiff(valuesinLG,valuesout)
      endif
 
 770  continue ! LG loop      
          
      if (alkali) then
         if (lptop.ne.-1) then
            deallocate(ubb_res)
            deallocate(ql_res)
         end if
         !if (.not.oldftps) then
            deallocate(ftps)
            deallocate(ftpsv)
         !end if
      end if
      deallocate(nchistart,nchistop)
      deallocate(ntime,nntime)
      deallocate(natompsnode,nrange,tscale)

!      call magmaf_finalize()       
     
 780  call clock(stop)
      print*,'Total CPU time:',stop-start
      call MPI_FINALIZE(ierr)
      STOP 'Job completed'
      END

      subroutine makecorepsi(nznuci,zasymi,ry,corepin,r0in)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/smallr/ formcut,regcut,expcut,fast,match
      common /noblgegas_switch/ i_sw_ng

      logical fast, match
      rmax = rmesh(meshr,1)
      r0 = 1.0
      corep = 0.0
      nnn = nznuc-nint(zasym)
      print*,'!!! nznuc, zasym:',nznuc,nint(zasym)
      if (nznuc-nint(zasym).eq.11) then
         call hfz11(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc-nint(zasym).eq.3) then
         call hfz3 (rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nnn .eq. 5  .and. i_sw_ng .eq. 0 ) then
         call hfz5(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
!         call hfz6(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nnn .eq. 7 .and. zasym .eq. 0) then
         call hfz7(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc .eq. 8 .and. nint(zasym) .eq. 1) then
         call hfz8(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc-nint(zasym).eq.9 .or.
     >        nznuc .eq. 10 .and. nint(zasym) .eq. 5) then
         if(i_sw_ng .eq. 0)  then 
            call hfz9 (rmax, meshr, rmesh, expcut, nznuc, corep, r0)
         elseif(i_sw_ng .eq. 1)  then 
            call hfz9_ng (rmax, meshr, rmesh, expcut, nznuc, corep, r0)
         else
            print*,'Wrong value for i_sw_ng',i_sw_ng
            stop
         endif
      else if (nznuc-nint(zasym).eq.17  .or.
     >        nznuc .eq. 18. and. nint(zasym) .eq. 5) then
         if(i_sw_ng .eq. 0)  then 
            call hfz17(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
         elseif(i_sw_ng .eq. 1)  then 
            call hfz17_ng(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
         else
            print*,'Wrong value for i_sw_ng',i_sw_ng
            stop
         endif 
      else if (nznuc-nint(zasym).eq.19 .and. i_sw_ng .eq. 0) then
         call hfz19(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc-nint(zasym).eq.29) then
         call hfz29(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc-nint(zasym).eq.31 .and. i_sw_ng .eq. 0) then
         call hfz31(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc.eq.36 .and. nint(zasym).eq.5) then
         call hfz35_ng(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc-nint(zasym).eq.37 .and. i_sw_ng .eq. 0) then
         call hfz37(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc-nint(zasym).eq.47 .and. i_sw_ng .eq. 0) then
         call hfz47(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc.eq.54 .and. nint(zasym).eq.5) then
         call hfz54_ng(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc-nint(zasym).eq.55 .and. i_sw_ng .eq. 0) then
         call hfz55(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc.eq.74.and.nint(zasym).eq.5) then
         call hfz74_5(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc-nint(zasym).eq.69) then
         call hfz69(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else if (nznuc-nint(zasym).eq.79) then
         call hfz79(rmax, meshr, rmesh, expcut, nznuc, corep, r0)
      else
         stop 'Unknown target'
      endif 
      print*,'Core wave functions'
      do lac = 0, lactop
         do nac = lac + 1, nabot(lac) - 1
            print'(2i3,i6,12x,f16.8,"  eV, ",f16.8,"  au")',
     >         lac,nac,istoppsinb(nac,lac),enpsinb(nac,lac)*ry,
     >         enpsinb(nac,lac)/2.0
         enddo 
      enddo 
      print*,'Initial state'
      nac = ninc !nabot(labot)
      lac = linc !0
      print'(2i3,i6,12x,f16.8,"  eV, ",f16.8,"  au")',
     >   lac,nac,istoppsinb(nac,lac),enpsinb(nac,lac)*ry,
     >         enpsinb(nac,lac)/2.0
      return
      end
      
      subroutine makeeigenpsi(nznuci,zasymi,nnbtop,l,corep,r0,
     >   ery,ry,enion,enlevel,ntstop,ovlpnl,nold)
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      COMMON/CONT/ CE(ncmax),CINT(ncmax),NCSTATES,ENERGY
      common/smallr/ formcut,regcut,expcut,fast, match
      logical fast, match
      dimension chi(maxr),ovlpnl(ncmax,0:lamax,nnmax)
      character opcl*10
      complex phase,ovlpnl
      common /there1Cowan/ there1
      logical there1

      rmax = rmesh(meshr,1)

      print*,'!!! nznuc, zasym:',nznuc,nint(zasym)
      print'(" LA NA  NRSTOP   E(1/cm)      E(eV)         ovlp")'
      if (nznuc-nint(zasym).eq.11) then
C  here for sodium like targets
c$$$      r0 = 1.439
c$$$      corep = 0.99
C  Si IV, i.e., Si 3+
c$$$            r0 = 1.3
c$$$            corep = 0.3
C  Al III, i.e., Al 2+
C  Mg II, i.e., Mg +
c$$$            r0 = 1.25
c$$$            corep = 0.5
C  Ar VIII, i.e., Ar 7+
c$$$            r0 = 1.3
c$$$            corep = 0.5
         nb = abs(nnbtop) - 2
         if (l.gt.2) nb = abs(nnbtop) - l
         if (nb.gt.0)
     >      call fcz11(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.3) then
C  Here for lithium like targets
C  Li I, i.e., Li
c$$$  r0 = 1.439
c$$$  corep = 0.1915
C  Be II, i.e., Be+
c$$$  r0 = 1.439  ?
c$$$  corep = 0.1915  ?
C  B III, i.e., B 2+
C  C IV, i.e., C 3+
C  N V, i.e., N 4+
c$$$  r0 = 1.439  ?
c$$$  corep = 0.1915  ?
C  O VI, i.e., O 5+
c$$$  r0 = 1.439  ?
c$$$  corep = 0.1915  ?
C  F VII, i.e., F 6+
c$$$  r0 = 1.439  ?
c$$$  corep = 0.1915  ?
         nb = abs(nnbtop) - 1
         if (l.gt.1) nb = abs(nnbtop) - l
         call fcz3(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.5) then
C     B, C+, N++
         nb = abs(nnbtop) - 1
         if (l.gt.1) nb = abs(nnbtop) - l
         call fcz5(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.9) then
         nb = abs(nnbtop) - 1
         if (l.gt.1) nb = abs(nnbtop) - l
         call fcz9(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.17) then
C  Here for argon II like targets
         nb = abs(nnbtop) - 3
         if (l.eq.2) nb = abs(nnbtop) - 2
         if (l.ge.3) nb = abs(nnbtop) - l 
         call fcz17(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.19) then
C  Here for potassium like targets
c$$$      r0 = 2.0
c$$$      corep = 5.354
C  Ca II, i.e., Ca+
c$$$      r0 = 2.0
c$$$      corep = 5.354
         nb = abs(nnbtop) - 3
         if (l.eq.2) nb = abs(nnbtop) - 2
         if (l.ge.3) nb = abs(nnbtop) - l 
         call fcz19(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.29) then
C  Here for Cu like targets
         nb = abs(nnbtop) - 3
         if (l.ge.3) nb = abs(nnbtop) - l 
         call fcz29(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.31) then
C  Here for Ga like targets
         nb = abs(nnbtop) - 3
         if (l.ge.3) nb = abs(nnbtop) - l 
         call fcz31(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.37) then
C  Here for rubidium like targets
c$$$      r0 = 2.0
c$$$      corep = 5.354
C  Ca II, i.e., Ca+
c$$$      r0 = 2.0
c$$$      corep = 5.354
         nb = abs(nnbtop) - 3 
         if (l.eq.2) nb = abs(nnbtop) - 2
         if (l.ge.3) nb = abs(nnbtop) - l 
         call fcz37(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)

      else if (nznuc-nint(zasym).eq.47) then
C  Here for cadmium like targets
         nb = abs(nnbtop) - 4
         if (l.eq.2) nb = abs(nnbtop) - 4
         if (l.ge.3) nb = abs(nnbtop) - l 

         call fcz47(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)

      else if (nznuc-nint(zasym).eq.55) then
C  Here for caesium like targets
c$$$      r0 = 2.0
c$$$      corep = 15.0
         nb = abs(nnbtop) - 5
         if (l.eq.2) nb = abs(nnbtop) - 4
         if (l.ge.3) nb = abs(nnbtop) - l 
         call fcz55(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc.eq.74.and.nint(zasym).eq.5) then  ! W VI target
c$$$      r0 = 2.0
c$$$      corep = 15.0
         nb = abs(nnbtop) - 5
         if (l.eq.2) nb = abs(nnbtop) - 4
         if (l.ge.3) nb = abs(nnbtop) - 4
         if (l.ge.4) nb = abs(nnbtop) - l
         call fcz74_5(rmax, meshr, rmesh, expcut, regcut, 0.0, 0,
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.69) then  ! Tm-like targets
c$$$      r0 = 2.0
c$$$      corep = 15.0
         nb = abs(nnbtop) - 5
         if (l.eq.2) nb = abs(nnbtop) - 4
         if (l.ge.3) nb = abs(nnbtop) - l 
         call fcz69(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >      chi, iwff, nb, l, nznuc, phase, corep, r0)
      else if (nznuc-nint(zasym).eq.79) then
C  Here for Hg-like targets
         nb = abs(nnbtop) - 5
         if (l.eq.3) nb = abs(nnbtop) - 4
         if (l.ge.4) nb = abs(nnbtop) - l 
         if(there1) then
            print*, 'We do not call ATOM frozen-core routines.'
            if(nold .ne. 1) then
               print*, 'nold=',nold,
     >            '   You should not call redefine routine, ',
     >            'nold is not 1, reset it to 1 in ccc.in'
               stop
            end if
            return
         else
            call fcz79(rmax, meshr, rmesh, expcut, regcut, 0.0, 0, 
     >         chi, iwff, nb, l, nznuc, phase, corep, r0)
         end if
      elseif (nznuc-nint(zasym).ne.1) then
         print*, 'Eigenstates for this target are not yet defined.'
      end if 

      nstart = max(nabot(l),l+1)
      if (l.eq.linc.and.ninc.lt.nstart) nstart = ninc
      dppoltot = 0.0
      dppoltotm = 0.0
      do n = nstart, min(abs(nnbtop), nnmax)
         if (nznuc-nint(zasym).eq.1) then
C  We are here for hydrogen like ions as well as hydrogen
            call rnl(nznuc,n,l,psinb(1,n,l),enpsinb(n,l),
     >         istoppsinb(n,l))
         end if
c$$$         if (n.eq.max(nabot(l),l+1)) then
c$$$            do i = 1, istoppsinb(n,l)               
c$$$               write(10*n+l,*) rmesh(i,1), psinb(i,n,l)
c$$$            enddo
c$$$            close(10*n+l)
c$$$         end if

         tsum = 0.0
         do i=1,istoppsinb(n,l)
            tsum = tsum + psinb(i,n,l)*psinb(i,n,l) * rmesh(i,3)
         end do
         ovlpnl(n,l,n) = tsum
         enion = 13.6058 * nznuc ** 2
         enlevel = enion * 8067.8
         elevel = (enpsinb(n,l) * ry + enion) * 
     >      enlevel / enion
         if (ery+enpsinb(ninc,linc)-enpsinb(n,l).gt.0.0)
     >      then
            opcl = 'open'
         else
            opcl = 'closed'
         endif 
         if (l.eq.1) then
            osc = oscil(enpsinb(n,l),psinb(1,n,l),istoppsinb(n,l),
     >      enpsinb(nabot(0),0),psinb(1,nabot(0),0),
     >      istoppsinb(nabot(0),0),rmesh,meshr)
            dppol = 4.0 * osc / (enpsinb(nabot(0),0)-enpsinb(n,l))**2
            dppoltot = dppoltot + dppol
            print '("         Oscillator strength:",f6.3," dppol:",
     >         f7.2," dppoltot:",f7.2)', osc,dppol,dppoltot
            osc = oscil_modif(enpsinb(n,l),psinb(1,n,l),istoppsinb(n,l),
     >      enpsinb(nabot(0),0),psinb(1,nabot(0),0),
     >      istoppsinb(nabot(0),0),rmesh,meshr)
            dppol = 4.0 * osc / (enpsinb(nabot(0),0)-enpsinb(n,l))**2
            dppoltotm = dppoltotm + dppol
            print '("Modified Oscillator strength:",f6.3," dppol:",
     >         f7.2," dppoltot:",f7.2)', osc,dppol,dppoltotm
         endif 
         print'(2i3,i6,f12.1,2f16.8,2x,a6)',
     >      l,n,istoppsinb(n,l),elevel,enpsinb(n,l)*ry,
     >      tsum,opcl
c$$$         if (abs(tsum-1.0).gt.1e-4.and.n.eq.nabot(l)) then
c$$$            if (nznuc-nint(zasym).ne.69)
c$$$     >         stop 'PROBLEM WITH WAVEFUNCTION'
c$$$         endif 
      end do
C General set of oscillator strengths, see Dmitry's notes
      if (l.gt.0) then
         print*,"nl->n'l' Oscillator strengths"
         do nm = nabot(l-1), 5
            do np = nabot(l), 5
               if (enpsinb(np,l)-enpsinb(nm,l-1).gt.0) then
                  print'(2i1,"->",2i1,1p,20e10.2)',nm,l-1,np,l,
     >               oscil(enpsinb(np,l),psinb(1,np,l),istoppsinb(np,l),
     >               enpsinb(nm,l-1),psinb(1,nm,l-1),
     >               istoppsinb(nm,l-1),rmesh,meshr)/(2.0*l-1)*l
               else
                  print'(2i1,"->",2i1,1p,20e10.2)',np,l,nm,l-1,
     >               oscil(
     >               enpsinb(nm,l-1),psinb(1,nm,l-1),istoppsinb(nm,l-1),
     >               enpsinb(np,l),psinb(1,np,l),istoppsinb(np,l),rmesh,
     >               meshr)/(2.0*l+1)*l
               endif
            enddo
         enddo
      endif
      return
      end
      
      subroutine makepspsi(nold,nznuci,zasymi,npsin,nnbtop,l,corep,r0,
     >   al,vdcore,ovlp,phasen,ery,ry,enion,enlevel,slowe,ne2e,ovlpnl,
     >   hlike,orzasym)
c$$$      use psinc_module
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      real*8 cknd, rlambda,rin(maxr),rout(maxr),vin(maxr),vout(maxr)
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      COMMON/CONT/ CE(ncmax),CINT(ncmax),NCSTATES,ENERGY
      common/smallr/ formcut,regcut,expcut,fast, match
      logical fast, match, exists, hlike
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      real, allocatable :: psinc(:,:,:)
      dimension chi(maxr),vdcore(maxr),als(ncmax),vasymp(maxr),
     >   ovlpnl(ncmax,0:lamax,nnmax), slowe(ncmax),
     >   vnucl(maxr),ovlp(ncmax,0:lamax),ovlpp(ncmax,ncmax),nsum(10)
      character opcl*10,chin*1,sumfile*4
      complex phasen(ncmax,0:lamax),phase,sigc,ovlpnl
      save etot
      common /relcorsetup/ iih_local_pot
      common /noblgegas_switch/ i_sw_ng
      common /double/id,jdouble(22)
      common /debye/dbyexists, dblp, rmudeb, zdby
      logical::dbyexists,jamalsv

      nps = npsin
      rmax = rmesh(meshr,1)
      do i = 1, meshr
c$$$         vnucl(i)=(vdcore(i)*rmesh(i,1)+float(nznuc-1-nint(zasym)))
c$$$     >      /rmesh(i,1)*2.0
c$$$         vnucl(i) = (vdcore(i)+float(nznuc-1)/rmesh(i,1)) * 2.0
         vnucl(i) = (vdcore(i)-rpow2(i,0)*(zasym+1.0))*2.0
      enddo
      maxpot = meshr
      if (zasym.ne.0.0) print*,'Check that the asym changes are working'
c$$$      z0 = float(nznuc) + zasym
      z0 = 0.0 !  zasym !float(nznuc)
c$$$      do i = 1, meshr
c$$$         write(42,*) rmesh(i,1), (vnucl(i)*rmesh(i,1)-2.0*z0)/2.0
c$$$      enddo
c$$$      stop
      alorig = al
      al = abs(al)
      slowery = abs(slowe(1)) / ry
      do n = 1, nps
         als(n) = al * 2.0
      enddo
c$$$      print*, 'Enter 1 for MAKEPS and 0 for HYD1EL'
c$$$      read*,nchoice
      nchoice = 1
      if (nchoice.eq.0) then
         do n = 1, nabot(l) + 1
            print*,'Enter Lambda for n =',n
            read*, als(n)
c$$$            als(n) = al * 2.0 * 1.2**(nabot(l) - n)
         enddo 
         call hyd1el(z0,l,nps,als,psen2,vnucl,ps2,minps2,maxps2)
      else 
         npsstates(1,l) = nps
         rlambda(1,l) = al * 2.0
         if (al.gt.4.0*(zasymi+1)) then
!            if (al*(zasymi+1)/rmesh(meshr,1).gt.5.0) then
            print*,'Using Box-based states'
            hmax = rmesh(meshr,2)
            zas = zasymi !+ 1.0
c$$$            if (dbyexists) then
c$$$               print*,'Debye radius and Z:',dblp, zdby
c$$$c$$$               vdcore(1:meshr) = -zas/rmesh(1:meshr,1)*
c$$$               vdcore(1:meshr) = -zdby/rmesh(1:meshr,1)*               
c$$$     >            exp(-rmesh(1:meshr,1)/dblp)
c$$$               zas = 0.0
c$$$            endif
c$$$            print*,'zas=',zas
            vnucl(:)=vdcore(:)*2.0 !Ry
            print*,'zas,vnucl(1)*rmesh(1,1):',zas,vnucl(1)*rmesh(1,1)
            ra = al
            nbmax = nps+l
            call pseudo(jdouble,id,hmax,zas,l,ra,nbmax,maxr,vnucl,
     >         psen2,ps2,jmax)
            minps2(:) = 1
            maxps2(:) = jmax
            print*,'R-check:',jmax,rmesh(jmax,1)
c$$$            if (l.eq.0) ps2(:,1)=0.0 ! For Anatoli
         else
c$$$            if (dbyexists) then
c$$$               print*,'Debye radius and Z:',dblp, zdby
c$$$c$$$               vnucl(1:meshr) = -2.0*z0/rmesh(1:meshr,1)*
c$$$               vnucl(1:meshr) = -2.0*zdby/rmesh(1:meshr,1)*
c$$$     >            exp(-rmesh(1:meshr,1)/dblp)
c$$$               z0 = 0.0 !Check that z0 is not needed again
c$$$            endif 
               
      inquire(file='vs_0.dat',exist=jamalsv)
      if (jamalsv) then
         open(42,file='vs_0.dat')
         i = 1
 10      read(42,*,end=20) rin(i),vin(i)
         i = i + 1
         go to 10
 20      continue
         close(42)
         rout(:) = rmesh(:,1)
         call intrpl(i-1,rin,vin,meshr,rout,vout)
         vnucl(:) = vout(:)*2.0
         z0 = 0.0
      endif 
c$$$            z0= 0.0
c$$$            vnucl(:)=-2.0/rmesh(:,1)
            z0= zasym + 1 !0.0
            vnucl(:)=vdcore(:)*2.0 !Ry
c$$$            vnucl(:)=-2.0/rmesh(:,1)
            print*,'z0,vnucl(1)*R(1):',z0,vnucl(1)*rmesh(1,1)
            call makeps(z0, .true., al, l, expcut, nps, ps2,
     >         psen2, minps2, maxps2, rmesh, vnucl,
     >         maxpot, meshr,cknd(1,nabot(l),l))

c$$$            if (l.eq.0) then !hollow Be+
c$$$               psen2(2) = psen2(1)
c$$$               ps2(:,2) = ps2(:,1)
c$$$            endif
c$$$            if (l.eq.0) then
c$$$               do n = 1, nps - 1
c$$$                  psen2(n) = psen2(n+1)
c$$$                  ps2(:,n) = ps2(:,n+1)
c$$$                  minps(n)=minps(n+1)
c$$$                  maxps(n)=maxps(n+1)
c$$$               enddo
c$$$            endif 
         endif 
         tmp = 1e10
         do n = 1, nps
            if (abs(slowery-psen2(n)).lt.tmp) then
               tmp = abs(slowery-psen2(n))
               nsmall = n
            endif
         enddo 
         
         if (nold.gt.1) call redefine(nold,ps2,psen2,minps2,maxps2,
     >      nps,psinb(1,1,l),enpsinb(1,l),istoppsinb(1,l),meshr,
     >      l,nabot,vdcore,rmesh)
         if(nznuc .eq. 80 .and. iih_local_pot .eq. 0) then ! Hg-target
            print*, 'Enter 0 if no rel.cor. to be included for Hg'
            read(*,*) iih         

            if(iih.ne.0)  call  redefine_rel(ps2,psen2,minps2,
     >         maxps2,nps,l,nznuc,nabot)      
         end if
  
         if (l.eq.linc) then
            etot = ery + psen2(ninc-l)
            if (nznuc.eq.6.and.zasym.eq.5.0) then
               etot = etot + 100.0/ry
               print*,'adjusting energy for CV to:',etot
            else 
               print*,'Setting total energy to (eV):',etot * ry
            endif 
c$$$            if (slowe.lt.0.0.and.etot.gt.0.0) slowe = etot / 2.0 * ry
         endif 
c$$$         if (slowe.ge.0.0.and.etot.gt.0.0.and.psen2(nps).gt.etot
c$$$     >      .and.alorig.lt.0.0) then
c$$$            n = 1
c$$$            do while (psen2(n).lt.etot)
c$$$               n = n + 1
c$$$            enddo
c$$$            distance = psen2(n) - psen2(n-1)
c$$$            smalle = min(0.4,distance/2.0)
c$$$            slowery = 0.0
c$$$            if (psen2(n-2).lt.0.0) then
c$$$               print*,'Setting one energy to:',(etot/2.0)*ry,n-1+l
c$$$               slowery = etot / 2.0
c$$$            else if (etot - psen2(n-1).lt.smalle) then
c$$$               print*,'Setting one energy to:',(etot - smalle)*ry,n-1+l
c$$$               slowery = etot - smalle
c$$$            else if (psen2(n) - etot .lt. smalle) then
c$$$               print*,'Setting one energy to:',(etot + smalle)*ry,n+l
c$$$               slowery = etot + smalle
c$$$            endif 
c$$$         endif 
         if (alorig.lt.0.0.and.etot.gt.0.0) then
            almin = al
            almax = al
            small = 1e10
            diffmin = -1e10
            diffmax = 1e10
            it = 0
            old = -1e1
            nf = 1
            alstep = al/50.0
            do while (small.gt.1e-5.and.it.lt.60)
c$$$            do while (abs((slowery-old)/old).gt.1e-5.and.it.lt.60)
c$$$            do while (abs(psen2(nf)+psen2(nf+1)-2.0*etot)/etot.gt.1e-4
c$$$     >         .and.it.lt.60)
               it = it + 1
               npsstates(1,l) = nps
               rlambda(1,l) = al * 2.0
               call makeps(z0, .true., al, l, expcut, nps, ps2,
     >            psen2, minps2, maxps2, rmesh, vnucl,
     >            maxpot, meshr,cknd(1,nabot(l),l))
               if (nold.gt.1)call redefine(nold,ps2,psen2,minps2,maxps2,
     >            nps,psinb(1,1,l),enpsinb(1,l),istoppsinb(1,l),meshr,l,
     >            nabot,vdcore,rmesh)
               if (slowe(1).ge.0.0) then
                  if (abs(slowe(1)/ry-etot*0.5)/etot.lt.1e-3) then
                     position = etot / 2.0
                  else 
                     print*,'Setting POSITION=ETOT'
                     position = etot
                     if (natop(l).eq.-99) position = etot / 2.0
                  endif 
                  n = 1
                  do while (psen2(n).lt.position.and.n.lt.nps)
                     n = n + 1
                  enddo
                  distance = psen2(n) - psen2(n-1)
                  nf = n - 1
                  if (psen2(n).lt.position) nf = n
                  old = slowery
                  slowery = 
     >               (position - distance / 2.0 + 4.0*psen2(n-1))/5.0
                  
c$$$                  if (natop(l).ge.-90.and.
c$$$c$$$     >               (psen2(n).lt.position))
c$$$     >               (slowery.lt.position/4.0.or.psen2(n).lt.position))
c$$$     >               slowery = (position / 4.0 + psen2(nf)) / 2.0
               endif                  
               if (it.lt.60) then
                  print*,it,
     >               ' slowery,alpha,almin,almax,small:',
     >               slowery*ry, al, almin,almax,small
               else
                  print*, 'CAUTION: it = 60'
               endif 
               small = 1e10               
               do n = 1, nps
                  diff = (slowery-psen2(n)) / slowery
                  if (abs(diff).lt.small) then
                     small = abs(diff)
                     nsmall = n
c$$$                     print*, n, diff, almin, al, almax, psen2(n) * ry
                  endif
               enddo
               
               if (psen2(nps).lt.etot.and.slowe(1).ge.0.0) then
                  al = al + alstep*it
                  print*,'New Lambda:', al * 2.0
               else 
                  diff = (slowery-psen2(nsmall)) / slowery
                  if (diff.lt.0.0) then
                     almax = al
                     diffmax = diff
                     if (diffmin.gt.0.0) then
c$$$  al = (almax + almin) / 2.0
                        al = (almax * diffmin - almin * diffmax) /
     >                     (diffmin - diffmax)
                     else
                        al = al - alstep
                        almin = al
                     endif 
                  else
                     almin = al
                     diffmin = diff
                     if (diffmax.lt.0.0) then
                        al = (almax * diffmin - almin * diffmax) /
     >                     (diffmin - diffmax)
c$$$  al = (almax + almin) / 2.0
                     else
                        al = al + alstep
                        almax = al
                     endif 
                  endif
               endif 
               if (abs(almax-almin)/almax.lt.1e-6) it = 100
            enddo                   
            if (slowe(1).ge.0.0.and.nsmall.ne.nf) 
     >         print*,'check energies for LA and ETOT (eV):',l,etot*ry
c$$$            print*, it, nsmall, almin, al, almax, psen2(nsmall)*ry,
c$$$     >         slowery * ry
            if (slowe(1).ge.0.0)
     >         print'(" Alpha redefined to be:",f10.6," after",i3,
     >         " iterations, Test:",f10.6)', al ,it,
     >         abs(psen2(nf)+psen2(nf+1)-2.0*position)/position
            if (nf.lt.nps.and.slowe(1).ge.0.0.and.
     >         abs(psen2(nf)+psen2(nf+1)-2.0*position)/position.gt.1e-4)
     >         print*, 'CAUTION: check test'
         endif 
      endif 
      if (etot.gt.0.0) then
         n = nps
         do while (psen2(n) .gt. etot .and. n.gt.1)
            n = n - 1
         enddo
         print*,'Energy test:',etot-psen2(n),psen2(min(n+1,nps))-etot
      endif

c$$$      do n = 1, nabot(l) - 1
c$$$         print*,'Energy and final I for n =',n,psen2(n) * ry,
c$$$     >      maxps2(n)
c$$$      enddo
      sumfile = 'sum'//char(l+ichar('0'))
      print*,'Checking for: ', sumfile
      inquire(file=sumfile,exist=exists)
      if (exists) then
         print*,'found'
         open(42,file=sumfile)
 100     read(42,*,end=200,err=200) n,(nsum(i),i=1,n)
         if (nsum(1).le.l) then
            nnn = l + 1 - nsum(1)
            n = n - nnn
            do i = 1, n
               nsum(i) = nsum(i + nnn)
            enddo
         endif 
         print*,'L and summed states:',l,(nsum(i),i=1,n)
         do i = 1, n
            nsum(i) = nsum(i) - l
         enddo 
         call sumup(nps,ps2,psen2,maxps2,n,nsum)
         go to 100
 200     continue
         close(42)
         n1 = 1
         n2 = 1
         do while (n1.lt.nps)
            n1 = n1 + 1
            if (maxps2(n1).gt.0) then
               n2 = n2 + 1
               psen2(n2) = psen2(n1)
               maxps2(n2) = maxps2(n1)
               do i = 1, maxps2(n1)
                  ps2(i,n2) = ps2(i,n1)
               enddo
            endif 
         enddo
         nps = n2
      else
         print*,'not found'
      endif 

      print'(" LA NA  NRSTOP   E(1/cm)      E(eV)         ovlp",
     >   "        proj/ovlp")'

      if (natop(l).eq.-99) then
         n = nps
         do while (psen2(n).gt.etot/2.0.and.n.gt.1)
            n = n - 1
         enddo
         natop(l) = n + l
      else if (natop(l).lt.-90) then
C  The following assumes that the number of summed states is -90-natop(l),
C  starting with the last
         n = nps
c$$$         do while (psen2(n) .gt. etot .and. nps-n .lt. -90-natop(l))
c$$$            n = n - 1
c$$$         enddo
C  The following assumes that the number of closed states is -90-natop(l)
         do while (psen2(n+90+natop(l)) .gt. etot)
            n = n - 1
         enddo
         
         nbeg = n
         en = psen2(n)
         do while (n .lt. nps)
            print*,'Summing state number:',n
            n = n + 1
            en = en + psen2(n)
         enddo
         nend = n
         en = en / (nend-nbeg+1)

         maxn = 0
         do j = 1, maxr
            chi(j) = 0.0
         enddo 
         do n = nbeg, nend
            if (maxps2(n).gt.maxn) maxn = maxps2(n)
            do j = 1, maxps2(n)
               chi(j) = chi(j) + ps2(j,n)
            enddo 
         enddo
         psen2(nbeg) = en
         do j = 1, maxn
            ps2(j,nbeg) = chi(j) / sqrt(float(nend-nbeg+1))
         enddo
         print*,'Redefining NPS (old,new)',nps, nps - (nend-nbeg)
         nps = nps - (nend-nbeg)
         natop(l) = 0
         do n = nbeg+1, nps
            psen2(n) = psen2(n+nend-nbeg)
            maxps2(n) = maxps2(n+nend-nbeg)
            do j = 1, maxps2(n+nend-nbeg)
               ps2(j,n) = ps2(j,n+nend-nbeg)
            enddo
         enddo
      endif

      ncontmin = 10000
      ncontmax = 0

      const = sqrt(2.0/acos(-1.0))
      ntop = min(nps + l,nnmax)
C  Define the projections of the pseudostate onto the discrete eigenstates
      do n = nabot(l), ntop
         ovlp(n,l) = 0.0
         ovlpnl(n,l,:) = 0.0
         ovlpnl(n,l,nabot(l)) = 1.0
         do ne = nabot(l), max(nabot(l),abs(nnbtop))
            sum = 0.0
            do i = 1, min(istoppsinb(ne,l),maxps2(n-l))
               sum = sum + psinb(i,ne,l) * ps2(i,n-l) * rmesh(i,3)
            enddo
            ovlpnl(n,l,ne) = sum
c$$$            if (etot.lt.0.0) ovlpnl(n,l,ne) = 1.0
            if (psen2(n-l).lt.0.0) ovlp(n,l) = ovlp(n,l) + sum * sum !1.0
            if (nnbtop.lt.0) then
               ovlp(n,l) = 1.0
               if (psen2(n-l).gt.0.0.or.ne.ne.nabot(l)) then
                  ovlpnl(n,l,ne) = 0.0
               else
                  ovlpnl(n,l,ne) = 1.0
               endif 
            endif 
         enddo
         if (abs(nnbtop)+2*ne2e.gt.nnmax) stop 'NNMAX is not big enough'
         do ne = 1, ne2e * 2
            if (ne.le.ne2e) then
               en = abs(slowe(ne)) / ry
            else
               en = etot - abs(slowe(ne-ne2e)) / ry
            endif
            if (en.lt.0.0) cycle
            if (ntype.eq.-3) then
               eta = -1.0 / sqrt(en)
               do i = 1, meshr
                  vasymp(i)=(vdcore(i)-1.0/rmesh(i,1))*2.0
               enddo 
               call regular(l,en,eta,vasymp,cntfug(1,l),l,rmesh,meshr,
     >            jdouble,id,regcut,expcut,chi,jstart,jstop,phase,sigc)
               do j = jstart, jstop
                  chi(j) = chi(j) * const
               enddo 
               phasen(n,l) = phase*sigc
            else            
               call hlikechi(nznuc,zasym,en,l,chi,phasen(n,l),jstart,
     >            jstop,meshr,rmesh,expcut,regcut,corep,r0)
            endif 
            sum = 0.0
            do i = jstart, maxps2(n-l)
               sum = sum + chi(i) * ps2(i,n-l) * rmesh(i,3)
            enddo
            ovlpnl(n,l,max(nabot(l),abs(nnbtop))+ne) = sum
         enddo 

c$$$         do ne = abs(nnbtop) + 1, abs(nnbtop) + ncstates
c$$$               sum = 0.0
c$$$               do i = 1, min(istoppsinb(ne,l),maxps2(n-l))
c$$$                  sum = sum + psinb(i,ne,l) * ps2(i,n-l) * rmesh(i,3)
c$$$               enddo
c$$$               ovlp(n,l) = ovlp(n,l) + sum * sum * cint(ne-abs(nnbtop))
c$$$          enddo
         if (etot.lt.0.0) ovlp(n,l) = 1.0
      enddo
C  If we are below the ionization threshold then set all overlaps to 1.0
C  so that we didn't get ionization cross sections below threshold!
      if (etot.lt.0.0)
     >  print*,'Note that the overlaps have been set to 1 for ETOT < 0'

C  Redefine the eigenstates with the pseudostates
      do n = l + 1, nabot(l) - 1
         elevel = (psen2(n-l) * ry + enion) *  enlevel / enion
         opcl = 'core'
         chin = '-'
         sum = 0.0
         print'(2i3,i6,f15.1,f15.5,2f14.8,2x,a6,1x,a)',
     >      l,n,maxps2(n-l),elevel,psen2(n-l)*ry,
     >      sum,sum,opcl,chin
      enddo

      dppoltot = 0.0
      dppoltotm = 0.0
      if (.not.allocated(psinc)) allocate(psinc(meshr,nnmax,0:lamax))
C Note that PSINC is available only for open states. Uncomment deallocate
C statement below to use in tmatcco.f
      psinc(:,:,:) = 0.0
      do n = nabot(l), ntop
         if (natop(l).eq.-99.and.psen2(n-l).gt.etot/2.0) then
            natop(l) = min(ntop, n)
            en = 0.0
            maxn = 0
            do j = 1, maxr
               psinb(j,n,l) = 0.0
            enddo 
            do i = n, ntop
               en = en + psen2(i-l)
               if (maxps2(i-l).gt.maxn) maxn = maxps2(i-l)
               do j = 1, maxps2(i-l)
                  psinb(j,n,l) = psinb(j,n,l) + ps2(j,i-l)
               enddo 
            enddo
            enpsinb(n,l) = en / (ntop - n + 1)
            istoppsinb(n,l) = maxn
            do j = 1, maxn
               psinb(j,n,l) = psinb(j,n,l) / sqrt(float(ntop - n + 1))
            enddo 
         else
            enpsinb(n,l) = psen2(n-l)
            istoppsinb(n,l) = maxps2(n-l)
            do i=1,istoppsinb(n,l)
               psinb(i,n,l)=ps2(i,n-l)
            enddo 
            do i = istoppsinb(n,l) + 1, maxr
               psinb(i,n,l) = 0.0
            enddo 
         endif 
         sum = 0.0
         do i = 1, istoppsinb(n,l)
            sum = sum + psinb(i,n,l)*psinb(i,n,l) * rmesh(i,3)
         end do
         en = enpsinb(n,l)
         if (ery+enpsinb(ninc,linc)-enpsinb(n,l).gt.0.0.and.en.gt.0e0
     >      .and.hlike) then
            if (ntype.eq.-3) then
               eta = -1.0 / sqrt(en)
               do i = 1, meshr
                  vasymp(i)=(vdcore(i)-1.0/rmesh(i,1))*2.0
               enddo 
               call regular(l,en,eta,vasymp,cntfug(1,l),l,rmesh,meshr,
     >            jdouble,id,regcut,expcut,chi,jstart,jstop,phase,sigc)
               do j = jstart, jstop
                  chi(j) = chi(j) * const
               enddo 
               phasen(n,l) = phase*sigc
            else
               call hlikechi(nznuc,zasym,en,l,chi,phasen(n,l),jstart,
     >            jstop,meshr,rmesh,expcut,regcut,corep,r0)
               print*,'phasen(n,l)',phasen(n,l)
            endif 
            j = jstart
            do while (j.lt.jstop.and.chi(j+1)/chi(j).gt.1.0)
               j = j + 1
            enddo
C  The following three methods for getting the correct normalization
C  factor all yield much the same OVLP(n,l). Given the T matrix <f|T|i>
C  for state <f| with energy Ef, then the corresponding "true" T matrix
C  T(Ef) = ovlp * <f|T|i>
            tmp1 = chi(j)/psinb(j,n,l)
            if (n-l.lt.nps.and.n-l-1.gt.0) then
               tmp2 = sqrt(sqrt(en)*4.0/
     >            (psen2(n-l+1)-psen2(n-l-1)))
            else
               tmp2 = 0.0
            endif 
            tmp = 0.0
            do i = jstart, istoppsinb(n,l)
               tmp = tmp + chi(i) * psinb(i,n,l) * rmesh(i,3)
            enddo
C  Note that the coulomb routine contains the CONST factor
            ovlp(n,l) = tmp
            print*,'ovlp,chi/psi,we,phase',tmp,tmp1,tmp2,phasen(n,l)
c$$$            print*,'OVLP,en,phase:',tmp,en,
c$$$     >         real(phasen(n,l)),imag(phasen(n,l))

            if (ncontmin .gt. n) ncontmin = n
            if (ncontmax .lt. n) ncontmax = n
            do np = nabot(l), ntop
               tmp = 0.0
               do i = jstart, maxps2(np-l)
                  tmp = tmp + chi(i) * ps2(i,np-l) * rmesh(i,3)
               enddo
               ovlpp(np,n) = tmp
               if (n.eq.np) psinc(1:meshr,n,l)=chi(1:meshr)/tmp
            enddo
c$$$            if (n.eq.nsmall+l) then
c$$$               open(42,file=char(l+ichar('0')))
c$$$               write(42,'("#",7x,100f12.5)') ry*enpsinb(n,l),ovlpp(n,n)
c$$$               do i = 1,istoppsinb(n,l)
c$$$                  write(42,'(1p,3e12.4)') rmesh(i,1),
c$$$     >               psinb(i,n,l)*ovlpp(n,n), chi(i)
c$$$               enddo
c$$$               close(42)
c$$$            endif 
         endif 
         elevel = (enpsinb(n,l) * ry + enion) * 
     >      enlevel / enion
         if (ery+enpsinb(ninc,linc)-enpsinb(n,l).gt.0.0) then
            opcl = 'open'
         else
            opcl = 'closed'
         endif
C  If natop(l) is -N, then open channels plus following N - 1 closed ones
C  will be included.
         if (natop(l).eq.-98.and.enpsinb(n,l).ge.etot/2.0) 
     >      natop(l) = min(ntop, n)
         if (natop(l).eq.-99.and.enpsinb(n,l).ge.etot/2.0) 
     >      natop(l) = min(ntop, n - 1)
         if (natop(l).lt.0.and.opcl.eq.'closed') 
     >      natop(l) = min(ntop, n - natop(l) - 2)
         if (n.le.natop(l).or.natop(l).le.0) then
            chin = '+'
         else
            chin = '-'
         endif      
         print'(2i3,i6,f12.1,f16.8,2f14.8,2x,a6,1x,a)',
     >      l,n,istoppsinb(n,l),elevel,enpsinb(n,l)*ry,
     >        sum,ovlp(n,l),opcl,chin
         if (l.eq.1) then
            osc = oscil(enpsinb(n,l),psinb(1,n,l),istoppsinb(n,l),
     >      enpsinb(nabot(0),0),psinb(1,nabot(0),0),
     >      istoppsinb(nabot(0),0),rmesh,meshr)
            dppol = 4.0 * osc / (enpsinb(nabot(0),0)-enpsinb(n,l))**2
            dppoltot = dppoltot + dppol
            print '("         Oscillator strength:",f6.3," dppol:",
     >         f7.2," dppoltot:",f7.2)', osc,dppol,dppoltot
         endif
C     ANDREY======================================================
c$$$         open (unit =55,file = 'epar')
c$$$         read (55,*) epar
c$$$         close (55)
c$$$         open (unit =55,file = 'e_3s',ACCESS = 'APPEND')
c$$$         write (55,*) epar, enpsinb(n,l)*ry
c$$$         close (55)
c$$$         if (l.eq.0.and.n.gt.8) stop '###ANDREY'
C     ANDREY======================================================         
      enddo 
      deallocate(psinc)
C  If natop(l) is 0 then all open and closed channels are included. The <=
C  is used in case all channels are open, but input latop < 0.
      if (natop(l).le.0) natop(l) = ntop
      if (ncontmin.le.ncontmax) then
         print'(''Overlaps between continuum functions and L2 states'')'
         print '(2x,900i10)', (n,n=ncontmin,ncontmax)      
         do np = nabot(l), ntop
            print'(i6,1p,900e10.2)',np,(ovlpp(np,n),n=ncontmin,ncontmax)
         enddo
         if (ne2e.gt.0) then
            print'(''Overlaps between e2e functions and L2 states'')'
            print '(2x,900i10)', (n,n=1,ne2e*2)      
            do np = nabot(l), ntop
               print '(i6,1p,900e10.2)', np,
     >         (ovlpnl(np,l,max(nabot(l),abs(nnbtop))+ne),ne=1,ne2e*2)
            enddo
         endif 
      endif
      if(i_sw_ng .eq. 1) then
         call noble1el(l,natop(l)-l)
      endif
      inquire(file='waves',exist=exists)
      if (exists) then
         open(42,file=char(l+ichar('0')))
         write(42,'("#",7x,100f12.5)')
     >      (ry*enpsinb(n,l),n=nabot(l),natop(l)) !nsmall+l,nsmall+l)
         do i = 1,istoppsinb(natop(l),l),10
            write(42,'(i5,1p,100e12.4)') i, rmesh(i,1),
     >         (psinb(i,n,l),n=nabot(l),natop(l)) !nsmall+l,nsmall+l)
         enddo
         close(42)
      endif 
      return
      end

      subroutine sumup(nps,ps2,psen2,maxps2,nt,nsum)
      include 'par.f'
      dimension ps2(maxr,ncmax),psen2(ncmax),maxps2(ncmax),nsum(10)

      do n = 2, nt
         do i = 1, maxps2(nsum(n))
            ps2(i,nsum(1)) = ps2(i,nsum(1)) + ps2(i,nsum(n))
            ps2(i,nsum(n)) = 0.0
         enddo
         psen2(nsum(1)) = psen2(nsum(1)) + psen2(nsum(n))
         maxps2(nsum(1)) = max(maxps2(nsum(1)),maxps2(nsum(n)))
         maxps2(nsum(n)) = 0
      enddo

      do i = 1, maxps2(nsum(1))
         ps2(i,nsum(1)) = 1.0/sqrt(float(nt))*ps2(i,nsum(1))
      enddo
      psen2(nsum(1)) = psen2(nsum(1)) / nt
      return
      end
      
      subroutine hlikechi(nznuc,zasym,en,l,chi,phasen,jstart,
     >   jstop,meshr,rmesh,expcut,regcut,corep,r0)
      include 'par.f'
      real chi(maxr), rmesh(maxr,3)
      complex phasen
            
      q = sqrt(en)
      nztest = nznuc-nint(zasym)
      phasen = (1.0,0.0)
      jstop = meshr
      rmax = rmesh(meshr,1)

      if (nztest.eq.1.or.zasym.eq.0.5) then
         call coulomb(zasym,en,l,chi,phasen,jstart,jstop)
      else if (nztest.eq.11) then
         call fcz11(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.3) then
         call fcz3(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.9) then
         call fcz9(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.17) then
         call fcz17(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.19) then
         call fcz19(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.29) then
         call fcz29(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.31) then
         call fcz31(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.37) then
         call fcz37(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.55) then
         call fcz55(rmax, meshr, rmesh, expcut, regcut, q, 1,
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nznuc.eq.74.and.nint(zasym).eq.5) then
         call fcz74_5(rmax, meshr, rmesh, expcut, regcut, q, 1,
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.69) then
         call fcz69(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      elseif (nztest.eq.79) then
         call fcz79(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >      chi, jstart, 0, l, nznuc, phasen, corep, r0)
      else
         stop 'no continuum wave defined for this target'
      end if
      return
      end


c
c     Below this point is interpolation routines only
c     (Justin's Stuff!)
c
      subroutine intrpv(lg,nsmax,vkeep,gkeep,nbnd,npkeep,nchkeep,
     >   nkeep,npk,nchtop,nnew,gk,vmat,pint)
      include 'par.f'
      dimension vkeep(nkeep,nkeep+1), vmat(nnew,nnew+1),
     >   gkeep(kmax,nchan), gk(kmax,nchan), npkeep(nchan+1),
     >   npk(nchan+1), nbnd(0:lmax), psii(maxr), psif(maxr)
      dimension vkeepch(200,200)
      dimension vmatch(200,200)
c$$$      character ch*1
c$$$      ch(n) = char(n + ichar('0'))

      call clock(timeint1)

c
c     this is the master interpolation subroutine.
c
      if (nchtop.ne.nchkeep) stop 'NCHTOP must be same as NCHKEEP'

      do nchi=1,nchtop
         do nchf=1,nchtop
c
c     checking to make sure the vkeepch and vmatch are big enough
c     (this error should not ever occure since at the current 
c     time (late 1999) there can not be more than 150 kgrid points
c     in the grid. 200 points have been allowed for here for 
c     future expansion of the kgrid)
c     In the event that this needs to be increased it also needs
c     to be increased in interpvd and interpvo.
c            
            if (npk(nchi+1)-npk(nchi).gt.200.or.
     >           npk(nchf+1)-npk(nchf).gt.200) then
               stop 'INCREASE size of vmatch in main.f'        
            endif
            
            if (npkeep(nchi+1)-npkeep(nchi).gt.200.or.
     >           npkeep(nchf+1)-npkeep(nchf).gt.200) then
               stop 'INCREASE size of vmatch in main.f'        
            endif

            if (nchi.eq.nchf) then
c
c     we have a diagonal vmatrix block element
c     we extract the appropriate part of the vkeep matrix
c     send it off to intrpvd to get interpolated, this then
c     returns a vmat block element which can be inserted 
c     appropriately in the vmat matrix.
c
               ii=0
               do i=npkeep(nchi),npkeep(nchi+1)-1
                  ii=ii+1
                  jj=0
                  do j=npkeep(nchf),npkeep(nchf+1)-1+1
                     jj=jj+1
                     vkeepch(ii,jj)=vkeep(i,j)
                  enddo
               enddo   
               
               call intrpvd(lg,nsmax,vkeepch,gkeep,nbnd,npkeep,
     >          nchkeep,nkeep,npk,nchtop,nnew,gk,vmatch,
     >          npkeep(nchi+1)-npkeep(nchi),
     >          npk(nchi+1)-npk(nchi),nchi,pint)

               ii=0
               do i=npk(nchi),npk(nchi+1)-1
                  ii=ii+1
                  jj=0
                  do j=npk(nchf),npk(nchf+1)-1+1
                     jj=jj+1
                     vmat(i,j)=vmatch(ii,jj)
                  enddo
               enddo   
               
            else
c
c     we have a non-diagonal vmatrix block element
c     in this case the block element is sent to interpvo to 
c     be interpolated.
c
              if (nchf.lt.nchi) then
                 ii=0
               do i=npkeep(nchi),npkeep(nchi+1)-1
                  ii=ii+1
                  jj=0
                  do j=npkeep(nchf),npkeep(nchf+1)-1+1
                     jj=jj+1
                     vkeepch(ii,jj)=vkeep(i,j)
                  enddo
               enddo
               else
                  ii=0
               do i=npkeep(nchi),npkeep(nchi+1)-1
                  ii=ii+1
                  jj=0
                  do j=npkeep(nchf)+1,npkeep(nchf+1)-1+1
                     jj=jj+1
                     vkeepch(ii,jj)=vkeep(i,j)
                  enddo
               enddo
              endif 
 
               call intrpvo(lg,nsmax,vkeepch,gkeep,nbnd,npkeep,
     >          nchkeep,nkeep,npk,nchtop,nnew,gk,vmatch,
     >          nchi,nchf,
     >          npkeep(nchi+1)-npkeep(nchi),
     >          npkeep(nchf+1)-npkeep(nchf),
     >          npk(nchi+1)-npk(nchi),
     >          npk(nchf+1)-npk(nchf),pint)
               
              if (nchf.lt.nchi) then
               ii=0
               do i=npk(nchi),npk(nchi+1)-1
                  ii=ii+1
                  jj=0
                  do j=npk(nchf),npk(nchf+1)-1+1
                     jj=jj+1
                     vmat(i,j)=vmatch(ii,jj)
                  enddo
               enddo
               else
                  ii=0
               do i=npk(nchi),npk(nchi+1)-1
                  ii=ii+1
                  jj=0
                  do j=npk(nchf)+1,npk(nchf+1)-1+1
                     jj=jj+1
                     vmat(i,j)=vmatch(ii,jj)
                  enddo
               enddo   
              endif
c
c     this next endif is the end of the diagonal element testing
c
            endif   
c
c     end do loops of nchi,nchf
c
         enddo
      enddo   

      call clock(timeint2)
      print*,'Time for interpolation was: ',timeint2-timeint1
c
c     this is a test section
c     only a valid test for dim(vkeep) = dim(vmat)
c     used for making sure interpolation is working ok.
c
c      open(57,file="check.dat")
c         do i=1,npk(nchtop+1)-1
c         do j=1,npk(nchtop+1)-1
c            write(57,*)'vkeep=',vkeep(i,j)
c            write(57,*)'vmat=',vmat(i,j)
c               if ((vkeep(i,j)-vmat(i,j)).ne.0.0.and.
c     >          vmat(i,j).ne.0.0) then
c            write(57,*)i,j,
c     >       (vmat(i,j)-vkeep(i,j))/vmat(i,j),
c     >       vmat(i,j),vkeep(i,j)
c               endif
c         enddo
c       enddo
c      close(57)

      return
      end

      subroutine intrpvo(lg,nsmax,vkeepb,gkeep,nbnd,npkeep,
     > nchkeep,nkeep,npk,nchtop,nnew,gk,vmat,nchi,nchf,
     > nkgridi,nkgridj,nvgridi,nvgridj,pint)
      include 'par.f'
      dimension vkeepb(200,200), vmat(200,200),
     >   gkeep(kmax,nchan), gk(kmax,nchan), npkeep(nchan+1),
     >   npk(nchan+1), nbnd(0:lmax), psii(maxr), psif(maxr)
      dimension vkgridi(nkgridi+1)
      dimension vkgridj(nkgridj+1)
      dimension tempi(nkgridj)
      dimension tempj(nkgridi)
      dimension xgrid(nvgridi*nvgridj) 
      dimension ygrid(nvgridi*nvgridj) 
      dimension vmatl(nvgridi*nvgridj)
      dimension vkeep(nkgridi+1,nkgridj+1)
      dimension vgridi(nvgridi)
      dimension vgridj(nvgridj)



c
c     this is a subroutine which performs the interpolation of 
c     the vkeep off-diagonal block elements to vmat.
c

c     making the vkeep matrix smaller from the bigger one that is
c     input

      do i=1,nkgridi
         do j=1,nkgridj
            vkeep(i,j)=vkeepb(i,j)
         enddo
      enddo

      do i=1,nkgridi
         vkeep(i,1)=0.
         vkeep(i,nkgridj+1)=0.
      enddo
      do j=1,nkgridj
         vkeep(1,j)=0.
         vkeep(nkgridi+1,j)=0.
      enddo
c
c     performing the mapping on the vkgrid/vkeep, if pint lt 0 then
c     this is a toggle to say that the mapping should no be
c     used. mapping is of the form
c     x -> (x/(1+x))**pint
c 
c     in the case of the mapping a point is added at infinity
c     in the case of no mapping a point is added at twice the
c     the maximum kgrid, so that we have interpolation rather
c     than extrapolation
c
      if (pint.lt.0.0) then
         do i=1,nkgridi
            vkgridi(i)=gkeep(i,nchi)
         enddo
         do j=1,nkgridj
            vkgridj(j)=gkeep(j,nchf)
         enddo
         vkgridi(nkgridi+1)=2*gkeep(nkgridi,nchi)
         vkgridj(nkgridj+1)=2*gkeep(nkgridj,nchf)
      else
         do i=1,nkgridi
            vkgridi(i)=(gkeep(i,nchi)/(1+gkeep(i,nchi)))**pint
         enddo
         do j=1,nkgridj
            vkgridj(j)=(gkeep(j,nchf)/(1+gkeep(j,nchf)))**pint
         enddo
         vkgridi(nkgridi+1)=1
         vkgridj(nkgridj+1)=1
      endif
      vkgridi(1)=0.
      vkgridj(1)=0.

c
c     the following commeneted out section is a reordering subroutin
c     that moved the onshell point to its appropriate location in 
c     the block element. It is commented out since its use in multi-
c     energy calculations is limited, however it is not deleted as
c     for a single energy speed improvement its use is valuable.
c      
c      oshelli=vkgridi(1)
c      itog=0
c      do i=2,nkgridi
c         if (oshelli.gt.vkgridi(i)) vkgridi(i-1) = vkgridi(i)
c         if (oshelli.lt.vkgridi(i).and.itog.eq.0) then
c            vkgridi(i-1)=oshelli
c            itog=1
c            ifit=i-1
c         endif
c      enddo
c
c      if (itog.eq.0) then
c         vkgridi(nkgridi)=oshelli
c         ifit=nkgridi
c      endif
c
c      oshellj=vkgridj(1)
c      itog=0
c      do i=2,nkgridj
c         if (oshellj.gt.vkgridj(i)) vkgridj(i-1) = vkgridj(i)
c         if (oshellj.lt.vkgridj(i).and.itog.eq.0) then
c            vkgridj(i-1)=oshellj
c            itog=1
c            jfit=i-1
c         endif
c      enddo
c
c      if (itog.eq.0) then
c         vkgridj(nkgridj)=oshellj
c         jfit=nkgridj
c      endif
c
c
c      do i=1,nkgridj
c         tempi(i)=vkeep(1,i)
c      enddo
c      do i=1,nkgridi
c         tempj(i)=vkeep(i,1)
c      enddo
c
c      do i=1,nkgridi
c         do j=1,nkgridj
c            if (i.lt.ifit.and.j.lt.jfit) then
c               vkeep(i,j)=vkeep(i+1,j+1)
c            elseif (i.lt.ifit.and.j.gt.jfit) then
c               vkeep(i,j)=vkeep(i+1,j)
c            elseif (i.gt.ifit.and.j.lt.jfit) then
c               vkeep(i,j)=vkeep(i,j+1)
c            endif
c         enddo
c      enddo
c
c      do i=1,nkgridj
c         if (i.lt.jfit) then
c            vkeep(ifit,i)=tempi(i+1)
c         elseif (i.eq.jfit) then
c            vkeep(ifit,jfit)=tempi(1)
c         else
c            vkeep(ifit,i)=tempi(i)
c         endif
c      enddo
c      
c      do j=1,nkgridi
c         if (j.lt.ifit) then
c            vkeep(j,jfit)=tempj(j+1)
c         elseif (j.eq.ifit) then
c            vkeep(ifit,jfit)=tempj(1)
c         else
c            vkeep(j,jfit)=tempj(j)
c         endif
c      enddo
c
c     end of reordering part.

c
c     firstly set up the arrays for the desired output points.
c     note that we must also set up the mapping or lack of it
c     on the vmat kgrid elements
c
      n=0
      do i=1,nvgridi
         vgridi(i)=gk(i,nchi)
      enddo
      if (vgridi(1).lt.0.0) vgridi(1)=0.0
      do j=1,nvgridj
         vgridj(j)=gk(j,nchf)
      enddo
      if (vgridj(1).lt.0.0) vgridj(1)=0.0

      if (pint.lt.0.0) then
         do i=1,nvgridi
            do j=1,nvgridj
               n=n+1
               xgrid(n)=vgridi(i)
               ygrid(n)=vgridj(j)
            enddo
         enddo
      else
         do i=1,nvgridi
            do j=1,nvgridj
               n=n+1
             xgrid(n)=(vgridi(i)/(1+vgridi(i)))**pint
             ygrid(n)=(vgridj(j)/(1+vgridj(j)))**pint
         enddo
      enddo
      endif

      call ITPLBV(6, nkgridi+1,nkgridj+1,vkgridi, vkgridj, 
     > vkeep,nvgridi*nvgridj,xgrid,ygrid,vmatl)

      do i=1,nvgridi
         do j=1,nvgridj
               vmat(i,j)=vmatl((i-1)*nvgridj+(j))
         enddo
      enddo

      return
      end


      subroutine intrpvd(lg,nsmax,vkeep,gkeep,nbnd,npkeep,nchkeep,
     >   nkeep,npk,nchtop,nnew,gk,vmat,nkgrid,nvgrid,ichan,pint)
      include 'par.f'
      dimension vkeep(200,200), vmat(200,200),
     >   gkeep(kmax,nchan), gk(kmax,nchan), npkeep(nchan+1),
     >   npk(nchan+1), nbnd(0:lmax), psii(maxr), psif(maxr)
      dimension vkeeps(nkgrid+1,nkgrid+1)
      dimension vkeept(nkgrid+1,nkgrid+1)
      dimension vkgrid(nkgrid+1)
      dimension tempa(nkgrid)
      dimension tempb(nkgrid)
      dimension vmats(nvgrid**2)
      dimension vmatt(nvgrid**2)
      dimension xgrid(nvgrid**2)
      dimension ygrid(nvgrid**2)
      dimension vgrid(nvgrid)
c$$$      character ch*1
c$$$      ch(n) = char(n + ichar('0'))

c
c     this is a subroutine which performs the interpolation of
c     the vkeep matrix diagonal elements to the vmat matrix.
c
c     changing the rectangular vkeep to two square matrices
c     vkeeps and vkeept (triplet(s) and singlet(t))
c
      do i=1,nkgrid
         do j=1,i
            vkeeps(i,j)=vkeep(i,j)
            vkeeps(j,i)=vkeep(i,j)
            vkeept(i,j)=vkeep(j,i+1)
            vkeept(j,i)=vkeep(j,i+1)
         enddo   
      enddo           

      do i=1,nkgrid
         vkeeps(1,i)=0.
         vkeeps(i,1)=0.
         vkeept(i,1)=0.
         vkeept(1,i)=0.
         vkeeps(nkgrid+1,i)=0.
         vkeeps(i,nkgrid+1)=0.
         vkeept(i,nkgrid+1)=0.
         vkeept(nkgrid+1,i)=0.        
      enddo

c
c     firstly we must get the kgrid for the vkeep.
c     again in this case we perform the mapping x->(1/(1+x))**pint
c     if pint is not negative
c
      if (pint.lt.0.0) then
         do i=1,nkgrid
            vkgrid(i)=gkeep(i,ichan)
         enddo
      vkgrid(nkgrid+1)=2*gkeep(nkgrid,ichan)
      else
         do i=1,nkgrid
            vkgrid(i)=(gkeep(i,ichan)/(1+gkeep(i,ichan)))**pint
         enddo
      vkgrid(nkgrid+1)=1.
      endif
      vkgrid(1)=0.
c
c     a small section of code to write the mapped kgrid parameters
c     to "kgrid.dat" this is important in seeing whether or not the
c     region is well spanned.
c
c      if (ichan.eq.1) then
c      open(57,file="kgrid.dat")
c      do i=1,nkgrid+1
c         write(57,*)i,vkgrid(i)
c      enddo
c      close(57)
c      endif

c     again the reordering of the vkeept matrix to include the
c     onshell point is commented out, see intrpvo for reasoning
c         now we perform the reodering process
c         firstly reorder the kgrid parameters
c         then reoder the matrices vkeeps and vkeept
c
c      oshell=vkgrid(1)
c      n=0
c      do i=2,nkgrid
c         if (oshell.gt.vkgrid(i)) vkgrid(i-1) = vkgrid(i)
c         if (oshell.lt.vkgrid(i).and.n.eq.0) then
c            vkgrid(i-1)=oshell
c            n=1
c            ifit=i-1
c         endif
c      enddo
c      
c      if (n.eq.0) then
c         vkgrid(nkgrid)=oshell
c         ifit=nkgrid
c      endif
c
c      do i=1,nkgrid
c         tempa(i)=vkeeps(1,i)
c         tempb(i)=vkeeps(i,1)
c      enddo
c
c      do i=1,nkgrid
c         do j=1,nkgrid
c            if (i.lt.ifit.and.j.lt.ifit) then
c               vkeeps(i,j)=vkeeps(i+1,j+1)
c            elseif (i.lt.ifit.and.j.gt.ifit) then
c               vkeeps(i,j)=vkeeps(i+1,j)
c            elseif (i.gt.ifit.and.j.lt.ifit) then
c               vkeeps(i,j)=vkeeps(i,j+1)
c            endif
c         enddo
c      enddo
c
c      do i=1,nkgrid
c         if (i.lt.ifit) then
c            vkeeps(ifit,i)=tempa(i+1)
c            vkeeps(i,ifit)=tempb(i+1)
c         elseif (i.eq.ifit) then
c            vkeeps(ifit,ifit)=tempa(1)
c         else
c            vkeeps(ifit,i)=tempa(i)
c            vkeeps(i,ifit)=tempb(i)
c         endif
c      enddo
c
c      do i=1,nkgrid
c         tempa(i)=vkeept(1,i)
c         tempb(i)=vkeept(i,1)
c      enddo
c
c      do i=1,nkgrid
c         do j=1,nkgrid
c            if (i.lt.ifit.and.j.lt.ifit) then
c               vkeept(i,j)=vkeept(i+1,j+1)
c            elseif (i.lt.ifit.and.j.gt.ifit) then
c               vkeept(i,j)=vkeept(i+1,j)
c            elseif (i.gt.ifit.and.j.lt.ifit) then
c               vkeept(i,j)=vkeept(i,j+1)
c            endif
c         enddo
c      enddo
c
c       do i=1,nkgrid
c         if (i.lt.ifit) then
c            vkeept(ifit,i)=tempa(i+1)
c            vkeept(i,ifit)=tempb(i+1)
c         elseif (i.eq.ifit) then
c            vkeept(ifit,ifit)=tempa(1)
c         else
c            vkeept(ifit,i)=tempa(i)
c            vkeept(i,ifit)=tempb(i)
c         endif
c      enddo

c
c     At this point the matrix vkeep has been broken into vkeeps and
c     vkeept 
c
c     firstly set up the arrays for the desired output points.
c     the vmat kgrid must also undergo the mapping
c
      n=0
      do i=1,nvgrid
         vgrid(i)=gk(i,ichan)
      enddo
      if (vgrid(1).lt.0.0) vgrid(1)=0.0
      if (pint.lt.0.0) then
         do i=1,nvgrid
            do j=1,nvgrid
               n=n+1
               xgrid(n)=vgrid(i)
               ygrid(n)=vgrid(j)
            enddo
         enddo
      else
         do i=1,nvgrid
            do j=1,nvgrid
             n=n+1
             xgrid(n)=(vgrid(i)/(1+vgrid(i)))**pint
             ygrid(n)=(vgrid(j)/(1+vgrid(j)))**pint
          enddo
       enddo
      endif
c
c     since we have two block elements (singlet and triplet)
c     we must call the interpolation subroutine twice.
c
      call ITPLBV(6, nkgrid+1,nkgrid+1,vkgrid, vkgrid, vkeeps,
     >    nvgrid**2,xgrid,ygrid,vmats)

      call ITPLBV(6, nkgrid+1,nkgrid+1,vkgrid, vkgrid, vkeept,
     >    nvgrid**2,xgrid,ygrid,vmatt)
c
c     remaking the vmat matrix
c

      do i=1,nvgrid
         do j=1,nvgrid+1
            if (j.gt.i) then
               vmat(i,j)=vmatt((i-1)*nvgrid+(j-1))
            else
               vmat(i,j)=vmats((i-1)*nvgrid+j)
            endif
         enddo
      enddo
c
c     this is a test subroutine
c     only a valid test for dim(vkeep) = dim(vmat)
c     WARNING this changes vkeep!!!
c
c      open(57,file="checkd.dat")
c         do i=1,nkgrid
c         do j=1,nkgrid
c            write(57,*)'vkeep=',vkeep(i,j)
c            write(57,*)'vmat=',vmat(i,j)
c            vkeep(i,j)=vmat(i,j)-vkeep(i,j)
c               if (vkeep(i,j).ne.0.0) then
c                  write(57,*)ichan,i,j,vkeep(i,j)
c               endif
c         enddo
c      enddo
c      close(57)

      return
      end

      subroutine  Born_tiecs(energy,nent,nstmax,BornICS)
      include 'par.f'

      integer nent,nstmax
      dimension BornICS(knm,knm)
      real*8  corepart, ecore, E
      common /corepartarray/ corepart(KNM,nicmax),ecore(nicmax)
      common /corearray/ nicm, ncore(nspmCI)
      common /helium/ latom(KNM), satom(KNM), lpar(KNM), np(KNM)
      integer latom, satom, lpar, np
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /CcoefsE/  E(KNM)


      write(138,'("energy=",F10.5,
     >     ", corepart(N,nic) is the probability ",
     >     " to find state N")') energy
      write(138,'("in the manifold based on core orbital i=1,nicm")')
      write(138,'("partstN(N) is the probability to find state N")')
      write(138,'("in the manifold based on the core orbitals ")') 
      write(138,
     >     '("   nic  ic=ncore(nic) lo(ic)    ko(ic)   e1r(ic,ic)")')
      
      do nic=1,nicm
         ic = ncore(nic)
         write(138,'(4(I5,5X),E12.5)') nic,ic,lo(ic),ko(ic),ecore(nic)
      end do
      write(138,'("  N      l s par   E(N)a.u.     ",12X,
     >     "corepart(N,nic)")')
      write(138,'(33X," nic= ",100(I10,3X))') (i, i=1,nicm)
      write(138,'(33X," l(i)=",100(I10,3X))') (lo(ncore(i)), i=1,nicm)


      write(138,'(" nf    ni  ICS  l(nf) s(nf) par(nf),  E   
     >     coreparts(nf)")')
      do ni=1,nent
         do nf=1,nstmax
            write(138,'(2I5,2X,1P,E12.4,3I5,0P,100(1X,E12.5))') 
     >           nf, ni, BornICS(nf,ni), latom(nf), satom(nf), 
     >           lpar(nf), real(E(nf)), (corepart(nf,k) , k=1,nicm)
         enddo
      enddo


      return 
      end

      
#define SLEEPY_BARRIER_TAG 678

      subroutine sleepy_barrier(comm)
      implicit none
      include 'mpif.h'
      integer :: comm
      integer, dimension(:), allocatable :: recvhandle,sendhandle,dataa
      integer, dimension(:,:), allocatable :: statuses
      integer :: n,me
      logical :: done
      integer :: i,ierr
#ifdef NEW_SLEEPY_BARRIER
      integer buf,req,n1,n2,nomp,OMP_GET_MAX_THREADS
#endif
c$$$      return
      call mpi_comm_size(comm,n,ierr)
      call mpi_comm_rank(comm,me,ierr)
#ifdef NEW_SLEEPY_BARRIER
      nomp=1 !max(1,OMP_GET_MAX_THREADS()) ! revert for many tasks per node
      n=(nomp*(1+me/nomp)-1)-(nomp*(me/nomp))+1
      n1=nomp*(me/nomp)
      n2=nomp*(1+me/nomp)-1
      allocate(sendhandle(n1:n2))
      if( mod(me,nomp).eq.0 ) then
         call mpi_irecv(buf,1,MPI_INTEGER,n1,SLEEPY_BARRIER_TAG,
     >         comm,req,ierr)
         do i=n1,n2
           call mpi_issend(1,1,MPI_INTEGER,i,SLEEPY_BARRIER_TAG,comm,
     >            sendhandle(i),ierr)
         end do
         call mpi_waitall(n,sendhandle,MPI_STATUSES_IGNORE,ierr)
         call mpi_wait(req,MPI_STATUSES_IGNORE,ierr)
      else
         call mpi_irecv(buf,1,MPI_INTEGER,n1,SLEEPY_BARRIER_TAG,
     >         comm,req,ierr)
         done=.false.
         call mpi_test(req,done,MPI_STATUSES_IGNORE,ierr)
         do while ( .not.done)
              call sleep(2)
              call mpi_test(req,done,MPI_STATUSES_IGNORE,ierr)
         enddo
         call mpi_test(req,done,MPI_STATUSES_IGNORE,ierr)
      endif
      deallocate(sendhandle)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
#else
      allocate(sendhandle(0:n-1))
      allocate(recvhandle(0:n-1))
      allocate(dataa(0:n-1))
      allocate(statuses(MPI_STATUS_SIZE,0:n-1))
      do i=0,n-1
        call mpi_irecv(dataa(i),1,MPI_INTEGER,i,SLEEPY_BARRIER_TAG,comm,
     >          recvhandle(i),ierr)
        call mpi_isend(1,1,MPI_INTEGER,i,SLEEPY_BARRIER_TAG,comm,
     >          sendhandle(i),ierr)
      end do
      done=.false.
      call mpi_testall(n,recvhandle,done,statuses,ierr)
      do while ( .not.done)
        call sleep(2)
        call mpi_testall(n,recvhandle,done,statuses,ierr)
      end do
      call mpi_waitall(n,sendhandle,statuses,ierr)
      deallocate(sendhandle,recvhandle,dataa,statuses)
#endif
      end subroutine sleepy_barrier

      subroutine date_and_time(date,time,zone,valuesin)
      character date*8,time*10,zone*5
      integer valuesin(8)
      valuesin(:)=0
      return
      end



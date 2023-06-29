C     File: structure.f
c-----------------------------------------------------------------
c**   main - configuration interaction program for helium **********
c***** numerical integration and numerical functions *************
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
c     subroutune  setmaxsp - set overall max s.p. state
c     subroutine inputdatamax
c     subroutine configsp
c     subroutine factorials
c     subroutine eigv(n,H,bb)
c     subroutine dipmom(lmaxhe,nstate,C,E)
c     subroutine r2s(Nmax,C)
      
      subroutine  structure(Nmax,nspmW,namax,pnewC,E,enionry,eproj,
     >   vdcore,slowe) 
c     Nmax  - out
c     nspmw - out
c     namax - out
c     pnewC - redundant
c     E     - out
c     enionry - in
c     eproj   - in
c     vdcore  - in
c     slowe   - in
      use vmat_module ! only for nodeid
      include 'par.f'
      character atom*20
      character string*40
      dimension lmaxhe(0:1,2)
      common /nstatearray/ nstate(0:lomax,0:1,2)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      double precision E(KNM)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /increarrange/ inc
      common /nhforb/ nhfo
      common /di_el_core_polarization/ gamma, r0, pol(nmaxr)      
      common /ionic_configurations/ l_ion_core, nk_ion_core(0:lomax)
      common /ngive_sym_opt/  ngivesym(nspmax), N_opt
      common /ngive_ionic_orb/  ngive_ion_orb(nspmax)
      common /include_opt/  include_opt(nspmCI,nspmCI)
c      common /cut_off_energies/  en_max, global_en_max, sym_en_max
      common /switch_manifold/  isw_manifold
      common /MPI_info/myid, ntasks, cnode, ench
      character cnode*3, ench*11
      real vdcore(maxr,0:lamax)
      logical exists
      integer pnewC
      
c     On entry to this routine 'nspm' is specified by Igor's routines.
c     If Igor's functions are not used and Laguerre functions are made then
c     'nspm' will be defined after call to setmaxsp(igorfunctions) routine.

c     Variable jsubset =1, =2) is used to differetiate between two possible modes
c     for the second go over F5 file: jsubset = 2 is for making optimized orbitals,
c     jsubset = 1 is for making f.c.o. corresponding to the states ordered in f.c.o. model
c     Variable nswitch is used for working on subset
      nswitch = 1
c     If nswitch = 0, then it indicates that program goes second time over F5 file.
      
c     Here we set number of optimized orbital to be zero. It is necessary for 
c     the routine set_basis(...)
      N_opt = 0
c     set 2-el. angular coef. array.
      call setfang

      open(15,file=adjustl(cnode//'states.core_parts'//ench),
     >   action='WRITE')
      open(4,file=adjustl(cnode//'F6'//ench),action='WRITE')
      open(10,file=adjustl(cnode//'he-states'//ench),action='WRITE')
      open(20,file=adjustl(cnode//'overlap.result'//ench),
     >   action='WRITE')
 10   open(3,file='F5',iostat=iostatF5,action='READ')

      if(iostatF5.ne.0) then
         print*, '********  File  F5  was not found'
         stop
      end if
 110  read(3,'(a20)') atom
      write(4,'(a20)') atom
      
      lmaxhe(:,:) = 0
      nstate(:,:,:) = 0
      
      l_ion_core = -1
      nk_ion_core(:) = 0

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
c     lmaxhe(0:1,2) -  maximum angular momentum of helium states
c     included in calculation (for singlet and triplet states separately).
c     nstate(l,s,ip=1) - number of required states with l=ll, s=ls,
c     lparity= (-1)**l.
c     nstate(l,s,ip=2) - number of required states with l=ll, s=ls,
c     lparity= - (-1)**l.
c
      gamma = 0.0
      rho = 1.0
c     if i_dipmom = 0  no oscillator strengths are calculated
c     if i_dipmom = 1  absorption oscillator strengths are calculated
c     if i_dipmom = 2  lp >= l oscillator strengths are calculated
c     if i_dipmom = 3  dipole polarazability of the states (l=0)  is calculated
c                      in this case only S and P states should be included.
c     if i_dipmom = 4  singlet triplet mixing is calculated
c     if i_dipmom = 5  singlet triplet mixing and osc. strength are calculated

      read(3,*)  nhfo, i_dipmom, i_r_sq, i_ort, i_gsort, gamma, r0
      write(4,'(" nhfo =",I3," i_dipmom =",I3," i_r_sq =",I3,
     >   " i_ort =",I3," , i_gsort =",I3)')
     >   nhfo, i_dipmom, i_r_sq, i_ort, i_gsort
      inquire(file='JMitroy_par',exist=exists)
      if (exists) then
         write(4,'("Set two-electron pol.pot parameters to ",
     >        "J.Mitroy choice: PRA 78(2008)012715")') 
         gamma = 0.4814
         r0 = 1.361
      endif
      write(4,'("Di-electron pol.pot.: gamma =",F10.5," r0 =",F10.5)') 
     >     gamma, r0

      if(nhfo.eq.1) then
         write(4,'(" numerical h.f. orbitals are used")')
      end if

      i_stmix = 0
      if(i_dipmom .eq. 4) then
         i_dipmom = 1
         i_stmix = 1
      endif
      if(i_dipmom .eq. 5) then
         i_dipmom = 1
         i_stmix = 2
      endif

c
      read(3,*)  inc, isubset, igorf, i_exclude, partN, en_max
      global_en_max = en_max
      
      i_orb_subset = 0
      if(isubset .eq. 3) then
         isubset = 1
         i_orb_subset = 1
      endif

      if(inc .eq. 3 .and. isubset .ne. 0) then    ! in F5 file: set i_exclude=0 and isubset=1
         nswitch = 1
      endif

c
c     If inc = 0, after CI calculation are performed the target states are
c                 represented  via ORIGINAL one-eletron functions.
c     If inc = 1, after CI calculation are performed the target states are
c                 represented  via NEW one-eletron functions which are built
c                 in the routine rearrange.f
c     If inc = 2, CI calculation are performed first for optimised states, 
c                 for which a set of optimised orbital is constructed. 
c                 Optimised orbitals and original one-electron functions are
c                 used then to perform second CI calculations.
c     If inc = 3, CI calculation are performed first for singlet-S state
c                 which is must be the first state and natural orbiatls are 
c                 formed for this state, then only a part of nat orbitals can 
c                 be used in following calculatiions (two more lines in F5 file)
c                !-> in F5 file: set i_exclude=0 and isubset=1
c                 [read array nk(l) from F5 file, defining how many nat. orb. will be used to substitute for original orbitals:      read(3,*) lnk, (nk(l), l=0,lnk) ]  
c                 [read array nkm(l) from F5 file, defining how many orbitals wil be used for each l:   read(3,*) lnkm, (nkm(l), l=0,lnkm) ]
c     If inc = 10, remove all interaction between frozen-core manifolds
c     If inc = 11, remove all interaction between frozen-core manifolds but make ground state with this interaction
c     If inc = 12, remove all interaction between frozen-core manifolds for states above set energy  e_set  but make all lower enery states with this interaction
      isw_manifold = 0
      if(inc. eq. 10) then
         inc = 1
         isw_manifold = 1
         write(4,'("isw_manifold = ",I5," setting inc=1")') 
     >        isw_manifold 
      elseif(inc. eq. 11) then
         inc = 1
         isw_manifold = 2
         write(4,'("isw_manifold = ",I5," setting inc=1")') 
     >        isw_manifold 
      elseif(inc. eq. 12) then
         inc = 1
         isw_manifold = 3
         write(4,'("isw_manifold = ",I5," setting inc=1")') 
     >        isw_manifold 
      endif


c     If nswitch=0 then programm goes over CI second time. This is the case when
c     F.C. orbitals for singlet states (isubset=1) or triplet states (isubset=3) 
c     are used to perform diagonalization or optimized orbitals have been made (inc=2, on 
c     the previous step of CI).  In this case we have set in routine 
c     start_structure(...)  nswitch=0 and igorfunctions=1. This way we avoid 
c     making new s.p.orbitals on this step.
      write(4,'("if inc=1 rearrange.f is on : inc =",I3)') inc
      if(nswitch.eq.1) igorfunctions = igorf
c
c     If isubset = 0  - standard calculation 
c     If isubset = 1(3)  - standard calculation , then second diagonalisation in the basis
c                       the singlet(triplet) state frozen-core  orbitals
c     If isubset = -1 - standard calculation +  the last state is summed up with all upper states
c     If isubset = 2  - standard calculation, then second diagonalisation in the basis
c                       the singlet state frozen-core  orbitals, then
c                       the last state is summed up with all upper states
c     If isubset = -1 or 2, Then the last state is summed up with all upper states
      if(isubset.eq.-1.or.isubset.eq.2) then
         write(4,'("isubset =",I3,", the last state for given ",
     >      "symmetry is summed up with all upper states")') isubset
c     
         if(isubset.eq.-1) then
            isubset = 0         ! no second diagonalization with usage of the singlet state
                                ! frozen-core  orbitals as new basis
         end if
         if(isubset.eq.2) then
            isubset = 1         ! second diagonalization with usage of the singlet state
                                ! frozen-core  orbitals as new basis, be aware of possible
                                !liniar dependencyof one-electron basis, make sure that
                                ! there were exluded states!
         end if
c            
         imake_summed_state = 1 ! this is a switch, setting it to be equal 1 will results in summation the last state  for given symmetry is summed up with all upper states:
         if(inc.eq.2) then
            print*,'Inconsistent use of  inc=2  and  isubset=-1'
            print*,'You probably need to set  inc=1'
            stop
         end if
      else
         imake_summed_state = 0 ! default
      end if     
      write(4,'("if isubset=1 (3) then subset of F.C. orbital is ",
     >   "used to perform diagonalization: isubset=",I3)') isubset
c
      write(4,'("if igorfunctions=1 then NO Laguerre functions are ",
     >   "made to perform diagonalization: igorfunctions=",I3)')
     >   igorfunctions
c      
      write(4,'("i_exclude =",I3,", partN =",F10.4,", en_max =",
     >   F10.5)') i_exclude, partN, en_max
c     if i_exclude=0, then no states are excluded due to partstN(N).le.partN
c                  = Ncore(.not.0), then coeffic. up to Ncore ionic orbital 
c                    will be summed to determine whether a state should be 
c                    to included or excluded. i_exclude will be mapped to 
c            j_summed = iexclude  - is how many ionic core orbital components 
c                                   should be summed.
c     partN - is the smallest allowed percentage/100 of the state in
c         the maniflod formed on s.p. functions that forms core (ionic core). 
c          The value of partN is compared with summed 
c         coeff. up to i_exclude of them to include or exclude  given state.
c     State N will be included if its partstN(N).gt.partN
c     if partN < 0, then all states are included.
c     en_max (eV) -  is the state energy (above frozen core)
c        above which states will be excluded,
c        if it is < 0 then en_max = (total energy in the incident channel = incident energy + 
c                                    incident state energy (above frozen core))
c     
      read(3,*) l_ion_core, (nk_ion_core(i), i=0,l_ion_core)
      write(4,'("l_ion_core=",I3,", nk_ion_core(i)=",20I3)') 
     >     l_ion_core,(nk_ion_core(i), i=0,l_ion_core)
c     l_ion_core < 0 - no limitations on outer electron configurations
c     l_ion_core >= 0 - outer electron allowed to be in continuum like config. 
c                       if inner electron (k2,l2) have k2 <= nk_ion_core(l2).
c                      for all other inner electron orbitals only correlation
c                      config. are allowed (k1,k2 <= nk2() )
 

c     here  Nmax  is determined (number of helium states).
      isum = 0
      do ip=1,2
         do is=0,1
            do l=0,lmaxhe(is,ip)
               isum = isum + nstate(l,is,ip)
            end do
         end do
      end do
      Nmax=isum
      write(4,'("Nmax =",I3)') Nmax
      if (nodeid.eq.1) print*, 'Nmax =', Nmax
      if (Nmax.gt.KNM) then
         print*,'Have too many target states, need to increse KNM to',
     >      ' at least',Nmax
         stop 'Have too many target states, need to increse KNM.'
      endif 


c     here we set all one-electron quantities - Igor's one-electron diag.
c     routine is not used if igorfunctions=0
c     if igorfunctions = 1 then in routine setmax Laguerre functions 
c     will not be made instead  functions from Igor programms 
c     or set of F.C. orbitals for triplet states will be used



c     If igorfunction=0 then Laguerre functions are made in this file,
      if(igorfunctions.le.0) then
         call setmaxsp(fl,maxf,minf)
      elseif(igorfunctions .eq. 2) then
         call set_dif_orb(fl,maxf,minf)
      else
c     If igorfunction=1 then Laguerre functions are NOT made and
c     s.p. functions from Igors programm or F.C. orbitals for triplet states
c     or optimized orbitals are used to perform diagonalization.
         write(4,'("igorfunctions=1, Laguere functions are not made")') 
         do n =1, 8
            read(3,'(a40)') string
            if (nodeid.eq.1) print '(''Skipping line: '',a40)', string
         enddo
      end if


c     if  igorfunctions.eq.0  then need to redefine  'nspmW'  to the new 
c     value of  'nspm'  that is available after call setmaxsp().
c     if  igorfunctions.eq.1  then  must be nspmW = nspm.
c     When optimized orbitals 
c     are made then after routine  optim_spf(...) is called the value of nspm 
c     is increased and on the second step of CI we need to set  nspmW = nspm and
c     reallocate memmory for array C(Nmax,nspmW,nspmW) in routine start_structure(...).

      nspmW = max(nspm,nspmW)
c
      call start_structure(atom,lmaxhe,i_dipmom,
     >   i_r_sq,i_ort,i_gsort,isubset,igorf,i_exclude, 
     >   partN,en_max,igorfunctions,nswitch,jsubset,i_orb_subset,
     >   enionry,eproj,Nmax,nspmW,namax,E,
     >   imake_summed_state,vdcore,slowe,i_stmix)
      
      if(jsubset.eq.1.and.nswitch.eq.0) then
         jsubset = 0
         go to 10
      end if
      if(jsubset.eq.2.and.nswitch.eq.0) then
         jsubset = 0
         go to 110
      end if
      
      return
      end
      
      subroutine start_structure(atom,lmaxhe,i_dipmom,
     >     i_r_sq,i_ort,i_gsort,isubset,igorf,i_exclude, 
     >     partN,en_max,igorfunctions,nswitch,jsubset,i_orb_subset,
     >     enionry,eproj,Nmax,nspmW,namax,E,
     >   imake_summed_state,vdcore,slowe,i_stmix)

      use vmat_module ! only for nodeid
      include 'par.f'
      dimension lmaxhe(0:1,2), nextsym(0:lamax)
      logical filled(0:3,0:lamax)
      common /nstatearray/ nstate(0:lomax,0:1,2)
      common /Nsarray/ Ns(0:lomax,0:1,komax,2)
      dimension los(nspmax),psen2(ncmax)
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      double precision E(KNM),CI(ncmCI,ncmCI),
     >   w(ncmCI),ortint
      common/orbc/ncm,no1(ncmCI),no2(ncmCI)
      common /ortog/  ortint(nspmax,nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common /increarrange/ inc
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      character atom*20
      double precision  Z, e1r
      common /Zatom/ Z
      common /nhforb/ nhfo
      common/hame1/ e1r(nspmax,nspmax)
      real partstN(KNM) 
      common /corearray/ nicm, ncore(nspmCI)      
      integer Num_ioncore(nspmCI)
      common /di_el_core_polarization/ gamma, r0, pol(nmaxr)      
      common /al_array/   al(0:lomax,2)
      double precision  al
      double precision, dimension(Nmax+1,nspmW,nspmW) :: C
      integer iorder_st(komax,0:lamax,0:1,2) ! used for setting state ordering based on ICS (Born)
      logical exists, done(0:lamax,0:1), fewer, equal
      double precision  coef,  En_summed_st
      double precision   sign_wf      
C     Igor's stuff
      real vdcore(maxr,0:lamax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      dimension nct(0:lamax,0:1,2)
      real psi(maxr)
      character chan(knm)*3
      common /charchan/ chan
      character  atom_label*3
      common /meshrr/ nr, gridr(nmaxr,3)
      common /theta_scat/ theta
       
      done(:,:) = .false.
      fewer = .false.
      equal = .true.
      Ns(:,:,:,:) = 0
      partstN(:) = 0.0

c     here helium states are enumerated
c
C  We first check that the input from optical.in is consistent with F5
c$$$      if (labot.ne.0.or.latop.ne.lmaxhe(0).or.latop.ne.lmaxhe(1)) then
c$$$         stop 'Ls do not match up between F5 and optical.in'
c$$$      endif
      
      labot = 0
      latop = max(lmaxhe(0,1),lmaxhe(1,1),lmaxhe(0,2),lmaxhe(1,2))
      do l = labot, latop
c$$$         if (natop(l) - nabot(l) + 1.ne.nstate(l,0)+nstate(l,1)) then
c$$$            stop 'Inconsistent use of nabot and natop'
c$$$         endif
C  nabot(l) must be set correctly in the ccc.in file. It is the principle
C  quantum number of the first active orbital. The core electrons range
C  from n = l + 1 to nabot(l) - 1.
         natop(l) = nabot(l) + nstate(l,0,1) + nstate(l,0,2) +
     >      nstate(l,1,1) + nstate(l,1,2) - 1
c         print*,'natop(l=',l,')=',natop(l)
         if (natop(l).gt.nnmax) then
            print*,'natop(l=',l,')=',natop(l),', nnmax=',nnmax
            stop 'INCREASE NNMAX'
         endif
      enddo
      
      do l=0,lamax
         nextsym(l) = 0
         nexts = 0
         do is=0,1
            do ip = 1, 2
               filled(nexts,l) = nstate(l,is,ip) .eq. 0
               nexts = nexts + 1
               nct(l,is,ip) = 0
            end do
         end do
      enddo
c$$$      n0 = 0
c$$$      n1 = 1
      nch = 0
      lg = 0
 100  nch = nch + 1
      call getchinfo (nch, ntmp, lg, psi, maxpsi, ea, la, na, li)
c
c$$$      print*,'la, nch : ', la, nch
      if (nch.ne.0) then
         nc = 1
         do while (filled(nextsym(la),la).and.nc.lt.5)
            nextsym(la) = mod(nextsym(la) + 1, 4)
            nc = nc + 1
         enddo 
         if(nc.gt.4) nch = 0
      endif 
      if (nch.ne.0) then
         if (nextsym(la).eq.0) then
            ip = 1
            is = 0
         elseif (nextsym(la).eq.1) then
            ip = 2
            is = 0
         elseif (nextsym(la).eq.2) then
            ip = 1
            is = 1
         elseif (nextsym(la).eq.3) then
            ip = 2
            is = 1
         endif
         nct(la,is,ip) = nct(la,is,ip) + 1
         nspin = is
         filled(nextsym(la),la) = nstate(la,is,ip) .eq. nct(la,is,ip)

C  Special case for helium (na.eq.1) so that s1S is followed by s2S.
         if (na.ne.1) nextsym(la) = mod(nextsym(la) + 1, 4)
         
         nc = 1         
         do while (filled(nextsym(la),la).and.nc.lt.5)
            nextsym(la) = mod(nextsym(la) + 1, 4)
            nc = nc + 1
         enddo 
         np(ntmp) = nabot(la) - 1 + nct(la,is,ip)

C  The following line makes sure that there the triplet S states begin with
C  a value one larger than the singlet S states.
         if (is.eq.1.and.la.eq.0) np(ntmp) = np(ntmp) + 1
         
C  The following line makes sure that there the singlet P states of even
C  parity begin with a value one larger than the triplet P states.
         if (is.eq.0.and.la.eq.1.and.ip.eq.2) np(ntmp) = np(ntmp) + 1
         
C  The following line makes sure that for Berillium the (2p)^2 ^1 D is
C  labaled as s2D, so that there was no confusion with the 2s3d ^{1,3}D states
         if (abs(z-4).lt.1e-4.and.is.eq.0.and.la.eq.2.and.ip.eq.1)
     >      np(ntmp) = np(ntmp) - 1
         ll(ntmp) = la
         ls(ntmp) = nspin               
c!!!  to be changed
         if (ip.eq.1) then
            lparity(ntmp) = (-1)**la
         else
            lparity(ntmp) = -(-1)**la
         endif
         if (nct(la,nspin,ip).gt.komax) stop "increase komax"
         Ns(la,nspin,nct(la,nspin,ip),ip) = ntmp         
c         print*,ntmp,na,la,nspin,' N, NA, LA, S'
         go to 100
      endif

      do la = labot, latop
         do nspin = 0, 1
            do ip = 1, 2
               if (nct(la,nspin,ip).ne.nstate(la,nspin,ip)) then
                  print*,"la,nspin,ip:", la,nspin,ip,
     >                 nct(la,nspin,ip),nstate(la,nspin,ip)
                  stop 'Enumeration of states is wrong (MAINHE)'
               endif
            enddo
         enddo
      enddo 
c
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
c
c     Gram-Schmidt orthogonalization of the s.p. basis
c     We should not call this routine when optimized orbitals have been made
c     This is case when nswitch=0 (do later: it also case when working on subset)
      if (i_gsort.eq.1.and.nswitch.ne.0) 
     >   call gsort(nspm,fl,maxf,minf,ortint,e1r,lo)
      
c
c     Set up the di-electron polarization array
      if(gamma.ne.0.0) 
     >     call di_el_pol
c
      if(Nmax+1.gt.KNM) then
         print*,'increase KNM to value of Nmax+1 =', Nmax+1
         stop
      end if
      write(15,'("states not included in the scattering calculation",
     >     " have  -  in first row  or N =",I3)') Nmax+1

      C(:,:,:) = 0d0

c     Here we set 'core' array: for frozen-core type calculations
      nicm = 0
      ncore(:) = 0
      Nofgfcst_total = 0
c-------------
C     Check if file Set.order_of_states is present. If it is, then inclusion and 
C     ordering of the target states for each symmetry is given in this file.

c     We should call this routine only when optimized orbitals have been made.
c     This is case when nswitch=0 (do later: it also case when working on subset)
      exists = 1 .eq. 0
      if(nswitch.eq.0) then
         inquire(FILE='Set.order_of_states',EXIST=exists)
         if(exists) then
            print*,'File Set.order_of_states exists, inclusion ',
     >         'and ordering of the target states for each ',
     >         'symmetry is given in this file.'
            print*,'Make sure that in file F5 i_exclude.ne.0, ',
     >         'and partN, en_max are chosen in such way that no ',
     >         'states are excluded.'
            write(4,'("File Set.order_of_states exists, inclusion ",
     >         "and ordering of the target states for each ",
     >         "symmetry is given in this file.")')
            write(4,'("Make sure that in file F5 i_exclude.ne.0, ",
     >         "and partN, en_max are chosen in such way that no ",
     >         "states are excluded.")')
            call read_Set_order_of_states(nstate,iorder_st)
         end if
      end if
c--------------------------------------------------------------------
c     here helium cofiguration interaction program is called
c     this programm is in the file  hel-3.f
      ry = 13.6058
      slowery = slowe/ry
      if (nodeid.eq.1)
     >   print*,'Calling configuration interaction routines'
      do nsll=1,nslm
         Nadd = 0
         write(4,*) nsll         
         iter = 0
         call hel(enionry,il,is,lpar,CI,w,fl,maxf,minf,iter)
         if (ncm.eq.0) cycle
         if (nsll.eq.1) etot = w(1) * 2.0 - enionry + eproj
         if (igorfunctions.le.-1.and..not.done(il,mod(is+1,2)).and.
     >      etot.gt.0.0) then
            if (-slowery.gt.etot) then
               slowery=etot/2.0
               slowe = slowery * ry
            endif 
            done(il,is) = .true.
            almin = al(il,1)
            almax = al(il,1)
            small = 1e10
            diffmin = -1e10
            diffmax = 1e10
            nf = 1
            alstep = al(il,1)/50.0
            position = etot
            if (slowery.gt.0.0) position = slowery
            do while (iter.lt.60.and.small.gt.1e-5)
               iter = iter + 1
               do n = 1, ncm
                  psen2(n) = w(n) * 2.0 - enionry
               enddo 
               call getnewal(small,diffmin,diffmax,alstep,position,
     >            slowery,psen2,ncm,diff,iter,almin,almax,al(il,1))
               call get_one_e_func(fl,maxf,minf)
               call hel(enionry,il,is,lpar,CI,w,fl,maxf,minf,iter)
               if (nsll.eq.1) etot = w(1) * 2.0 - enionry + eproj
            end do
         endif 
         ip = 1
         if((-1)**il.ne.lpar) ip = 2
         if (nodeid.eq.1)
     >      print*,' NCM, NSTATE(il,is,ip):', ncm,nstate(il,is,ip)
         if (ncm.lt.nstate(il,is,ip)) fewer = .true.
         equal = equal .and. ncm .eq. nstate(il,is,ip)
c     sort the ionic core orbitals by their energies
c     form Num_ioncore(i) array that keep the order number of the core orbitals
c     and will be passed to the routine stNcore(). 
c     This array will be sorted by insertion sort algorithm.
         do i=1,nicm
            Num_ioncore(i) = i
         end do
         if(i_exclude.ne.2) then
            do i=2,nicm
               itmp = Num_ioncore(i)
               j = i
               do while(j.ge.2.and.
     >            e1r(ncore(itmp),ncore(itmp)).lt.
     >            e1r(ncore(Num_ioncore(j-1)),ncore(Num_ioncore(j-1))))
                  Num_ioncore(j) = Num_ioncore(j-1)
                  j = j - 1
               end do
               Num_ioncore(j) = itmp
            end do
         end if
c$$$         write(4,'(" Energies of the ion-core orbitals:")')
c$$$         do j=1,nicm
c$$$            nc = ncore(Num_ioncore(j))
c$$$            write(4,'(" j,Num_ioncore(j),nc,lo(nc),e1r(nc,nc):",
c$$$     >         4i5,e15.5)') j,Num_ioncore(j),nc,lo(nc),e1r(nc,nc)
c$$$         end do

c     
         write(15,'("    Nc  N label l s par  E(eV)  partstN(N)",9X,
     >        "corepart(N,nic)")')
         write(15,'(29X," nic = ",10I6)') (i, i=1,nicm)
         write(15,'(29X," l(i) =",10I6)') (lo(ncore(i)), i=1,nicm)
         
         j_summed = i_exclude
 15      iplus = 0              ! if a state is not included then iplus = iplus + 1
         Nofgfcst = 0           ! how many good frozen core states have been included
         Nofgfcst_add = 0       ! how many states need to be added Nofgfcst_add 
c                                 goes from 1 to Nadd = nstate(il,is,ip) - Nofgfcst
         i_summed_states = 0    ! how many upper states have been added to the last state
c                               ! of the given symmetry
         En_summed_st    = 0d0  ! used to set correct energy of the summed state.
c         
         do num=1,ncm
            if(num.le.nstate(il,is,ip)) then
               N = Ns(il,is,num,ip)
               if(il.ne.ll(N).or.is.ne.ls(N).or.lpar.ne.lparity(N)) then
                  write(*,'("wrong enumeration of the states, stop ",
     >                 "in  structure()")')
                  print*,"il,is,lpar:", il,is,lpar
                  print*,"N=", N
                  print*,"ll,ls,lparity:", ll(N),ls(N),lparity(N)
                  stop
               end if
            else
               N = Nmax + 1
               do n2d=1,nspm
                  do n1d=1,nspm
                     C(N,n1d,n2d) = 0d0
                  end do
               end do
            end if
 20         inum = num + iplus
            
            if(num.le.nstate(il,is,ip).and.inum.gt.ncm) then
               if(j_summed.gt.nicm) then
                  print*,'Too many states have been excluded: inum=',
     >                 inum,'>ncm=',ncm,' decrease partN in file F5'
                  print*, 'There were', Nofgfcst, ' states with  ',
     >                 'partstN(N) > partN'
                  print*, 'for target symmetry',
     >                 '(l,s,par)=', il, is, lpar
                  stop
               else
c     At this run more states were odered then were included for given value
c     of partN compared with a value obtained when the number (j_summed) of 
c     coeffic. at ionic core functions were summed.
c     To include more states j_summed is increased by 1.
                  j_summed = j_summed + 1 
                  if(j_summed.gt.nicm) then
                     print*,'Too many states have been excluded: ',
     >                    'j_summed.gt.nicm'
                     stop
                  end if
                  Nadd = nstate(il,is,ip) - Nofgfcst
                  write(15,'("required number of states=",I3,
     >                 ", current number of states=",I3,", Nadd=",I3,
     >                 ", j_summed=",I3)') nstate(il,is,ip), Nofgfcst,
     >                 Nadd, j_summed
                  go to 15
               end if
            end if

            if(inum.gt.ncm) go to 25

            E(N) = w(inum)
c            sign = 1.0d0
c            if(CI(1,inum).lt.0.0d0) sign=-1.0d0
            call resolve_sign(Nmax,N,ncm,inum,CI(1,inum),sign_wf)

            do nc=1,ncm
               n1 = no1(nc)
               n2 = no2(nc)
c     Note : antisymmetrization is included through CI coefficients.
               sym = dsqrt(2.0D0)
               if(n1.eq.n2) sym = 1.0D0
               C(N,n1,n2) = CI(nc,inum)/sym*sign_wf
               C(N,n2,n1)=C(N,n1,n2)*
     >            (-1)**(lo(n1) + lo(n2) - il - is)
            end do
            if(i_exclude.ne.0) then
c     j_summed  - is how many ionic core orbital components should be summed,
c                 it should be less then nicm.
c     suppose the current state should be included
               jpartN = 1  
               if(exists) then  ! check if calling this routine is necessary
                  call check_state(iorder_st(1,il,is,ip),
     >               nstate(il,is,ip),inum,jpartN)
               end if
c     if the current state must be excluded due to routine  check_state() then jpartN = 0 is set in check_state().
c     pass this information to the routine stNcore().
               call stNcore(Nmax,nspmW,N,C,il,is,lpar,E(N),enionry,
     >              inum,Nofgfcst_add,jpartN,partN,en_max,partstN,
     >              j_summed,Nadd,Num_ioncore)
c--------------------------------------------------------------------------
c     This code is to sum all states of given symmetry in the last included
c     state of the given symmetry.
               if(inum.gt.nstate(il,is,ip).and.
     >            imake_summed_state.eq.1) then
c     find the state number for the last included state of the given symmetry
                  num_max = nstate(il,is,ip)
                  N_max = Ns(il,is,num_max,ip)
c     the current state here has the state number  N = Nmax + 1
                  i_summed_states = i_summed_states + 1
                  En_summed_st = En_summed_st + E(N)
                  do n2d=1,nspm
                     do n1d=1,nspm
c                        write(*,'(4i4,2e15.5)'),n1d,n2d,N_max,N,
c     >                     real(C(N_max,n1d,n2d)),real(C(N,n1d,n2d))
                        C(N_max,n1d,n2d) = C(N_max,n1d,n2d) +
     >                     C(N,n1d,n2d)
                     end do
                  end do
               endif
c--------------------------------------------------------------------------
               if(jpartN.eq.1) Nofgfcst = Nofgfcst + 1
               if(jpartN.eq.0) then
c     state excluded if jpartN=0
                  iplus = iplus + 1
                  C(N,:,:) = 0d0
                  go to 20
               end if
            end if
         end do                 !  end   num   loop
 25      continue   
c     print*, '***** Nofgfcst=', Nofgfcst
         if(Nofgfcst.gt.nstate(il,is,ip)) then
            write(10,'("Warning: ",I3,"  states were not included ",
     >           "for target symmetry (l,s,par)=", 3I3,
     >           " below energy (relative to ion core) en_max=",F10.5)')  
     >           Nofgfcst - nstate(il,is,ip), il, is, lpar, en_max
            write(10,'("To include those states increase  ",
     >           "nstate(il,is,ip)  in  file F5.")')
            write(10,*)
         end if
         Nofgfcst_total = Nofgfcst_total + Nofgfcst
c
c     If the last included state of the given symmetry was summed with
c     upper states then we need to restore correct normalizaton to this
c     new state : number of summed states is (i_summed_states+1)
         if(i_summed_states.gt.0) then
            print*,'i_summed_states=', i_summed_states
            coef = dsqrt(dble(i_summed_states+1))
            num_max = nstate(il,is,ip)
            N_max = Ns(il,is,num_max,ip)
            E(N_max) = (E(N_max) + En_summed_st)/dble(i_summed_states+1)
            do n2d=1,nspm
               do n1d=1,nspm
                  C(N_max,n1d,n2d) = C(N_max,n1d,n2d)/coef
c                  print*,n1d,n2d,C(N_max,n1d,n2d)
               end do
            end do
         end if
c         
      end do                    ! end  nsll  loop
      
      if (fewer) stop 'NCM < NSTATE for at least one symmetry'
      if (equal) then
         if (nodeid.eq.1) print*, 'NCM = NSTATE for all symmetries!'
      else 
         if (nodeid.eq.1)
     >      print*,'NCM>NSTATE for some symmetries, theta must be zero.'
      endif 
c     make optimized s.p.orbitals for each symmetry
      if(inc.eq.2.and.nswitch.eq.1) then
         call optim_spf(Nmax,nspmW,C,atom)
c     return and read further in F5 file new setup which is attached just below.
         igorfunctions=1 
         nswitch=0
         jsubset = 2
         return
      elseif(inc .eq. 3 .and. nswitch.eq.1) then
         call natorb(fl,minf,maxf,Nmax,nspmW,C)
c     return and read further in F5 file new setup which is attached just below.
         igorfunctions=1 
         nswitch=0
         jsubset = 2
         return 
      end if

      
      if(inc.eq.2) inc = 1
      

c     this subroutine must be called before all others subroutine if inc=0.
      if (inc.eq.0) 
     >     call  CItmp(Nmax,nspm,nspmW,C)

c     
      if (nodeid.eq.1) print*, 'Calling printcoreparts'
      call printcoreparts(partN,Nmax,nspmW,C,E,Nofgfcst_total,partstN)
      if (nodeid.eq.1) print*, 'Finish printcoreparts'
c-------------------------------------------------------------------
      if(i_ort.eq.1.and.inc.eq.0)
     >     call orthe(Nmax,nspmW,C)
c-------------------------------------------------------------------
c     set energy levels data
      call setdata(atom)
c     here dipole moments are calculated <|(r1+r2)/2|>
      if(i_dipmom.gt.0.and.inc.eq.0)
     >     call dipmom(i_dipmom,atom,Nmax,nspmW,C,E,fl,minf,
     >     maxf,lo,enionry)
c-------------------------------------------------------------------
c     here r**2 is calculated: 
      if(i_r_sq.eq.1.and.inc.eq.0)
     >     call r2s(Nmax,nspmW,C,fl,minf,maxf,lo)
c-------------------------------------------------------------------
c     Check orthogonality of the s.p. basis in case of theta!=0 scatetring calculations
c     This check is required to make sure that make_po() is called with orthogonal basis 
      if(theta .ne. 0 .and.nswitch.eq.0) then
         do n1=1,nspm
            do n2=1,n1-1
               if(abs(ortint(n1,n2)) .gt.1e-5 ) then
                  print*, 'nonortogonal s.p. basis: could be an ',
     >                 'error for theta != 0 scattering calculations'
                  print*, 'n1, n2:', n1, n2, ortint(n1,n2)
         print*, 'There is most likely an error in your calculations'
                  stop
               endif
               
            enddo
         enddo
      endif
c     call routine for making projection operator for nonuniqueness.
      call make_po(nspm,nr,fl,maxf,minf,lo,Nmax,nspmW,C)


      if (inc.eq.1) then
c     here s.p. states are rearranged in case of frozen core treatment of HE.
         write(4,'("here s.p. states are rearranged")')
         write(4,'("frozen core treatment of ",A20,"is assumed")') atom
         write(6,'("frozen core treatment of ",A20,"is assumed")') atom
         call find_nsp_tot(Nmax,nspmW,nspm,C,lo,nsp_tot)
         call rearrange(Nmax,nspmW,nspm,nsp_tot,lo,ko,
     >        fl,maxf,minf,C,los)
         write(10,'("******************************************")') 
         write(10,'("rearranged s.p. functions: HE state overlaps")') 
         if(i_ort.eq.1)
     >        call orthe(Nmax,nspmW,C)
         if(i_dipmom.gt.0)
     >        call dipmom(i_dipmom,atom,Nmax,nspmW,C,E,fl,minf,
     >        maxf,lo,enionry)
         if(i_r_sq.eq.1)
     >        call r2s(Nmax,nspmW,C,fl,minf,maxf,lo)
         if(isubset.eq.1.and.nswitch.eq.1) then
            write(4,*)
            write(4,'("*******************************************",
     >         "*********************")')
            write(4,*)
            if(i_orb_subset .eq. 1) then
               write(4,'("F.C. orbitals for triplet states are ",
     >              "used as single particle basis",/,
     >              "to perform diagonalization")')
            else
               write(4,'("F.C. orbitals for singlet states are ",
     >              "used as single particle basis",/,
     >              "to perform diagonalization")')

            endif
            call subset(fl,minf,maxf,los,i_orb_subset)

            if (nodeid.eq.1)
     >         print*,'F.C. orbitals for singlet states are used as',
     >         'single particle basis to perform diagonalization'
            igorfunctions=1  
            nswitch=0

c      need to recompile for the case of isub1 = 0

c     isub1 = 0
            isub1 = 1
c     if isub1=0 then the same set-up of F5 file is used for next CI calculation
            if(isub1.eq.0) then
               close(3)
               jsubset = 1
               return
c     if isub1=1 then new set-up of F5 file is used for next CI calculation.
c     This new set-up should be placed just after usual set up of the F5 file.
c     Note, in this case the F5 file must be read to the place where the new set-up
c     is placed (in first run better have only that 'tables for states' in 
c     F5 that are used)
            else  if(isub1.eq.1) then
               jsubset = 2
               return
            end if
            else
               jsubset = 0
         end if
      endif
c-------------------------------------------------------------------      
c     resolve sign of target states - might change the sign of a target state
c     only in the case when nonorthogonal basis was used for diagonalization.
      if (nodeid.eq.1) print*, 'Calling resolve_sign_final'
      call resolve_sign_final(Nmax,nspmW,C)
      if (nodeid.eq.1) print*, 'Finish resolve_sign_final'
c      
      if (nodeid.eq.1) print*, 'Calling find_major_config'
      call find_major_config(Nmax,nspmW,C)
      if (nodeid.eq.1) print*, 'Calling printCcoef'
      call  printCcoef(Nmax,nspmW,C,E)
c     
      if(i_dipmom.eq.0) then
         if(atom.eq.'helium') then
            call dataHe
         else  if(atom.eq.'beryllium') then
            call  dataBe
         else  if(atom.eq.'magnesium') then
            call  dataMg
         else  if(atom.eq.'calcium') then
            call  dataCa
         else  if(atom.eq.'zinc') then
            call  dataZn
         else  if(atom.eq.'barium') then
            call  dataBa
         else  if(atom.eq.'mercury') then
            call  dataHg
         end if
      end if

      if (nodeid.eq.1) print*, 'Calling print_energy'
      call print_energy(Nmax,E,enionry,atom)

c-------------------------------------------------------------------
      call find_namax(Nmax,namax)
      write(4,*)
      write(4,'("namax =",I5)') namax
      
c-------------------------------------------------------------------
      call F90_make_new_C(C,Nmax,nspmW,namax)

      if(i_stmix .ne. 0) then
         i_ng = 0
         call Spin_orbit(Nmax,E,enionry,i_stmix,i_ng)
      endif
c      if (atom.eq.'mercury'.or.atom.eq.'barium'
c     >     .or. atom .eq. 'ytterbium') then
c      endif 
c         call Spin_orbit(Nmax,E,vdcore,enionry)
c      call mix_s_t(maxr, nr, lamax, Nmax, KNM, atom_label, chan, 
c     >   ll, ls, vdcore)
      
c-------------------------------------------------------------------
      close(3)
      close(4)
      close(10)
      close(20)
      close(15)

      
C      call TwoElectWFValue(3,100,150)

c      STOP
      return
      end
      
      subroutine getnewal(small,diffmin,diffmax,alstep,position,
     >   slowe,psen2,nps,diff,it,almin,almax,al)
      real*8 al
      dimension psen2(nps)
      n = 1
      do while (psen2(n).lt.position.and.n.lt.nps)
         n = n + 1
      enddo
      if (n.eq.1) stop 'iteration failed in GETNEWAL'
      distance = psen2(n) - psen2(n-1)
      nf = n - 1
      if (psen2(n).lt.position) nf = n
      if (slowe.eq.0.0) then
         slowery = 
     >      (position + distance / 2.0 + 4.0*psen2(n))/5.0
      else
         slowery = slowe
      endif
      small = 1e10
      do n = 1, nps
         diff = (slowery-psen2(n)) / slowery
         if (abs(diff).lt.small) then
            small = abs(diff)
            nsmall = n
            if (small.lt.1e-5) then
               print*,'Alpha set at:',al
               return
            endif 
c$$$                     print*, n, diff, almin, al, almax, psen2(n) * ry
         endif
      enddo
               
      if (psen2(nps).lt.position) then
c$$$         al = al + alstep*it
         al = al - alstep
         small = 0.0
      else 
         diff = (slowery-psen2(nsmall)) / slowery
         if (diff.lt.0.0) then
            almax = al
            diffmax = diff
            if (diffmin.gt.0.0) then
c$$$  al = (almax + almin) / 2.0
               al = (almax * diffmin - almin * diffmax) /
     >            (diffmin - diffmax)
            else
               al = al - alstep
               almin = al
            endif 
         else
            almin = al
            diffmin = diff
            if (diffmax.lt.0.0) then
               al = (almax * diffmin - almin * diffmax) /
     >            (diffmin - diffmax)
c$$$  al = (almax + almin) / 2.0
            else
               al = al + alstep
               almax = al
            endif 
         endif
      endif
      print*,'reset AL:', al,psen2(nsmall)*13.6058,nsmall
      return
      end
      
c-------------------------------------------------------------------
      subroutine find_major_config(Nmax,nspmW,C)
      include 'par.f'
      double precision  C(Nmax+1,nspmW,nspmW),ortint
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
                        do i=1,nam(N)
                           n1 = na(i,N)
                           if(ncore(nic).eq.n1) then   ! n1 is a core orbital
                              do j=i,nam(N)
                                 n2 = na(j,N)
                                 ln2 = lo(n2)
                                 tmp = C(N,i,j)                           
                                 if(n1.ne.n2) tmp = tmp * dsqrt(2d0)
                                 tmp=tmp*tmp*ortint(n1,n1)*ortint(n2,n2)
                                 big(ln2) = big(ln2) + tmp
c                                 print*,'N,n1,n2,big(l)',N,n1,n2,
c     >                              (big(l),l=0,lomax)
                              end do
                              tmp = 0.0 ! find the largest config. for given core orbital n1=ncore(nic)
                              do l=0,lomax
                                 if(tmp.lt.big(l)) then
                                    tmp = big(l)
                                    l_tmp = l
c                                    print*,'N,n1,n2,l,tmp:',N,n1,n2,l,tmp
                                 end if
                              end do
                              c_big(nic) = tmp
                              l_big(nic) = l_tmp
                           end if
                        end do
                     end do
c     We found that for ionic core orbital n1=ncore(nic)
c     the largest CI contribution c_big(nic)  comes 
c     from configurations with orbital angular momentum l_big(nic)
c     Find which core orbital has largest CI contribution
                     f_big = 0d0
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
c     here the case (for Ba) of (5d ns) config. requires that s' orbital has n larger
c     than the largest n in previous ionic core orbitals
c     say, 2 orbitals with l=0 :  6s,7s  -->  n > 7
                     do nic1 = 1, m1 - 1
                        if(lo(ncore(nic1)).eq.l2orb(N).and.
     >                     num(lic,m1).le.ko(ncore(nic1))) then
                           num(lic,m1) = num(lic,m1) + 1
                        end if
                     end do
                     
                     n1orb(N) =  ko(ncore(m1)) + nabot(l1orb(N)) - 1
                     n2orb(N) =  num(lic,m1) + nabot(l2orb(N)) - 1
                     if((-1)**(il+is).eq.-1.and.
     >                    l1orb(N).eq.l2orb(N)) then
                        n2orb(N) = n2orb(N) + 1
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
      subroutine printCcoef(Nmax,nspmW,C,E)
      include 'par.f'
      parameter (ic_20 = 20)
      double precision  C(Nmax+1,nspmW,nspmW),  E(KNM),
     >     ortint, tmp, tmp1, sum
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      character chan(knm)*3
      common /charchan/ chan
      character lorb(0:23)*1 !(0:lamax)*1
      data lorb /"s","p","d","f","g","h","i",
     >     "j","k","l","m","n","o","P","q","r",
     >     "s","t","u","v","w","x","y","z"/      
      common /major_config/ l1orb(KNM),l2orb(KNM),n1orb(KNM),
     >   n2orb(KNM),config_maj(KNM)
      common /corearray/ nicm, ncore(nspmCI)
      real tmp_c(ic_20)
      integer icp_c(ic_20)
      common /increarrange/ inc
c     
      write(10,'("Nmax =",I3)') Nmax
c$$$      if (lamax.gt.23) then 
c$$$         print*, 'lorb(23) needs increasing'
c$$$         print*, 'at least to lamax=',lamax
c$$$         stop 
c$$$      endif
      do N=1,Nmax
         write(10,'(A3,", N=",I3,"  l=",I3,"  s=",I3," parity=",I3,
     >      ", energy=",F12.5,", major config.:(",2(I2,A1),")",
     >      F10.4)')
     >      chan(N), N, ll(N), ls(N), lparity(N), E(N),
     >      n1orb(N), lorb(l1orb(N)), n2orb(N), lorb(l2orb(N)),
     >      config_maj(N)

         write(10,'("   n1   n2    l1  l2",8X,"-",7X,"CI weights")')
         sum = 0d0
         do j=1,nam(N)
            do i=1,nam(N)
               tmp = C(N,i,j)
               n1 = na(i,N)
               n2 = na(j,N)
               if(n2.le.n1.and.tmp.ne.0.0D0) then                 
                  if(n1.ne.n2) tmp = tmp * dsqrt(2d0)
                  tmp1 = tmp*dsqrt(ortint(n1,n1)*ortint(n2,n2))
                  sum = sum + tmp1*tmp1     
                  if(n1.eq.n2) then
                     write(10,'(2I5,2X,2I4,2X,2E12.4)')
     >                  n1,n2,lo(n1),lo(n2),
     >                  real(tmp1),real(tmp1*tmp1)
                  else
c     here we find contributions to orbital n1 from ion core orbitals,
c     while (1.0 - sum_c)  will give contribution from all the rest orbitals.
                     sum_c = 0.0
                     kl = 0
                     if (inc.eq.1) then
                     do ic=n2+1,nicm
                        if(lo(ic).eq.lo(n1)) then
c     note: ion core orbitals (n2) form orthonomal set
                           kl = kl + 1
                           if(kl.gt.ic_20) then
                              print*,'increase size of arrays: ic_20,',
     >                           'file structure.f,",
     >                           " subroutine printCcoef'
                              stop
                           end if
                           tmp_c(kl)=ortint(n1,ic)/dsqrt(ortint(n1,n1))
                           icp_c(kl) = ic
                           sum_c = sum_c + tmp_c(kl) * tmp_c(kl)
                        end if
                     end do
                     end if
                     kl_m = kl
c
                     if(kl_m.eq.0) then
                        write(10,'(2I5,2X,2I4,2X,2E12.4)')
     >                     n1,n2,lo(n1),lo(n2),
     >                     real(tmp1),real(tmp1*tmp1)
                        else                           
                           write(10,'(2I5,2X,2I4,2X,2E12.4,
     >                        "; rest**2:",
     >                        E10.3,20(",",I3,":",E10.3:))')
     >                        n1,n2,lo(n1),lo(n2),
     >                        real(tmp1),real(tmp1*tmp1),1.0 - sum_c,
     >                        (icp_c(kl),tmp_c(kl), kl=1,kl_m)
                        end if
                     end if
               end if
            end do
        end do
        write(10,'("sum =",F10.6)') sum
      end do
      return
      end
c-------------------------------------------------------------------
      subroutine find_namax(Nmax,namax)
      include 'par.f' 
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      
      m = 0      
      do N=1,Nmax
           m = max(m,nam(N))
      end do
      namax = m

      return
      end
c-------------------------------------------------------------------
      subroutine CItmp(Nmax,nspm,nspmW,C)
      include 'par.f'
      double precision  C(Nmax+1,nspmW,nspmW)
      integer na, nam,  Nmax
      integer i,j,n1,n2      
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
c
      double precision  tmpC(nspmW,nspmW)
      
c      nspm = nspmW
      

c     initialize CItmp  arrays:
      do N=1,Nmax
         do i=1,nspm
            na(i,N) = 0
         end do
         nam(N) = 0
      end do

c     Reduced case:
c     Make all necessary CI arrays - case of only "a" electron,
c     only arrays nam(N), na(N,i) are defined.
      do N=1,Nmax
c     for each state N count how many s.p. orbitals are involved:
c     they are those which have at least one nonzero CI coefficient.
         i = 0
         do n1=1,nspm
            n1test = 0
            do n2=1,nspm
               if(C(N,n1,n2).ne.0d0) then
                  n1test = 1
                  go to 10
               end if
            end do
 10         i = i + n1test
            if(n1test.ne.0) then
               na(i,N) = n1
            end if
         end do
         nam(N) = i
c         print*,'******CItmp:  N,nam(N)', N,nam(N)
         if(nam(N) .eq. 0)  then
            print*, 'For N=', N, ' you have got nam(N)=0'
            do n1=1,nspm
               do n2=1,nspm
                  if(C(N,n1,n2).ne.0d0)
     >               print*,'n1,n2, C(N,n1,n2)', n1,n2, C(N,n1,n2)
               enddo
            enddo                  
            stop
         endif
c     
         do i=1,nam(N)
            n1 = na(i,N)
            do j=1,nam(N)
               n2 = na(j,N)
               tmpC(j,i) = C(N,n2,n1)
            end do      
         end do      
c     
         do n1=1,nspm
            do n2=1,nspm
               C(N,n1,n2) = 0d0
            end do
         end do

         do j=1,nam(N)
            do i=1,nam(N)
               C(N,j,i) =  tmpC(j,i)
            end do      
         end do      
c     
      end do

c
      
c$$$      do N=1,Nmax
c$$$         print*,N, ', nam(N) =', nam(N)  
c$$$         do i=1,nam(N)
c$$$            print*, '   ', i, ', na(i,N) = ', na(i,N) 
c$$$         end do
c$$$      end do

      return
      end
c-------------------------------------------------------------------
c     This routine makes optimized single particle orbitals for the description 
c     of the requested target states and then set a basis consisting of 
c     old basis plus all optimized orbitals (placed on top of the old basis).
c     N_opt  -  number of states to be optimized for given symmetry (il,is,ipar)
c     Nmax - the number of states (of all symmetries)
c     C(Nmax+1,nspmW,nspmW) - CI array
      subroutine optim_spf(Nmax,nspmW,C,atom)
      include 'par.f'
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common/hame1/ e1r(nspmax,nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      common /meshrr/ nr, gridr(nmaxr,3)
      double precision  C(Nmax+1,nspmW,nspmW)
      double precision  sym, sum, sumort, tmp, tmp2, ortint, e1r
      integer ltmp(lomax+1)
      common /corearray/ nicm, ncore(nspmCI)
      common /corearray_lsp/ nicm_lsp(0:lamax,0:1,-1:1), 
     >     ncore_lsp(nspmCI,0:lamax,0:1,-1:1)
      common /ngive_sym_opt/  ngivesym(nspmax), N_opt
      common /ngive_ionic_orb/  ngive_ion_orb(nspmax)
c     common /ionic_configurations/ l_ion_core, nk_ion_core(0:lomax)
      character atom*20
      integer n_state(0:lomax,0:1,-1:1), local_nstate(0:lomax,0:1,-1:1)
c     
      double precision, dimension(:,:), allocatable ::  A
      integer, dimension(:), allocatable :: nA
      real, dimension(:), allocatable :: sn_fl

c      
c     The common block  /corearray/ nicm, ncore(nspmCI) must be set and
c     available at this stage.
c     In the start of this subroutine the nspm is number 
c     of s.p.functions (Laguerre or generated in Igor's programm) and
c     must be less then nspmCI.
c     After f.c.orbitals are made in this routine the nspm is number
c     of f.c.orbital plus number of core functions and 
c     must be less then nspmax. (Usually nspmCI<<nspmax)

      write(20,'("******************************************")') 
      write(20,'("rearranged s.p. functions ")') 
      write(4,'("Make optimized s.p.orbitals")')
      print*, 'Making optimized s.p.orbitals and',
     >     ' starting new CI calculations'

c
      allocate(A(nspmax,nspmax))
      allocate(nA(nspmax))
      allocate(sn_fl(nspmax))
      A(:,:) = 0d0
      nA(:) = 0
      sn_fl(:) = 0.0

c     restore the CI coef. to be as for antisymmetric configurations
c     C(N,{n1n2})|{n1n2}> rather than to be C(N,n1,n2)|n1,n2> where
c     |n1,n2> is not antisymmetrized. Do it only for the current symmetry (il,is,ipar)
      tmp = dsqrt(2.0d0)
      do N=1,Nmax
         do n1=1,nspm
            do n2=1,nspm
               sym = tmp 
               if(n1.eq.n2) sym=1.0d0
               C(N,n1,n2) = C(N,n1,n2) * sym
            end do
         end do
      end do
      
c------------------------------------------------------------------------------------
      if(atom.eq."helium") then
         write(4,'("!!! Special choice for He optimised states",
     >      "have been made. See file structure.f, routine ",
     >      "optim_spf(Nmax,nspmW,C,atom), array local_nstate(,,)")')
         do ip=-1,1         
            do is=0,1
               do l=0,lomax
                  n_state(l,is,ip) = 0
                  local_nstate(l,is,ip) = -1
               end do
            end do
         end do
         local_nstate(0,0,1) = 2  ! sS
         local_nstate(1,0,-1) = 2  ! sP
         local_nstate(2,0,1) = 1  ! sD
         local_nstate(3,0,-1) = 1  ! sF
         local_nstate(4,0,1) = 1  ! sG         
         local_nstate(0,1,1) = 1 ! tS
         local_nstate(1,1,-1) = 2  ! tP
         local_nstate(2,1,1) = 1  ! tD
         local_nstate(3,1,-1) = 1  ! tF         
         local_nstate(4,1,1) = 1  ! tG         
      end if
c------------------------------------------------------------------------------------
c     Here frozen core orbitals are made. 
c     They are placed on the top of the current s.p.basis.
c     We assume that s.p. basis is  orthonomal. Then f.c.o. (coresponding to
c     a core orbital nc) is orthogonal to all core orbital involved for
c      given state.              ! ( ncp<=nc.)
      nspN = nspm
      do N=1,Nmax
         il_st = ll(N)
         is_st = ls(N)
         ip_st = lparity(N)
         nicm_st = nicm_lsp(il_st,is_st,ip_st)
c------------------------------------------------------------------------------------
         if(atom.eq."helium".and.
     >      local_nstate(il_st,is_st,ip_st).gt.0) then            
            n_state(il_st,is_st,ip_st)=n_state(il_st,is_st,ip_st) + 1
            if(n_state(il_st,is_st,ip_st).gt.
     >         local_nstate(il_st,is_st,ip_st))  nicm_st = 1
         end if
c------------------------------------------------------------------------------------
c     do nic=1,nicm_st
c     We assume here that nic=1 ionic core orbital always used for 
c     discretization of continuum in the consequent CI calculations.
c     Note: It would be better to first work out which ionic core orbitals will
c     be used for discretization of continuum in the consequent CI calculation
c     and only then exlcude these ionic core orbitals from building
c     associated with them optimized orbitals, but the array nk_ion_core(l)
c     is read later in the code.
c     Note: Exclusion of the optimised orbital builded on ionic core orbital nic=1
c     also assumes that all basis functions are used to build configurations.
c     If not all of them used then optimised orbitals with nic=1 should be
c     constructed and used.
C     Start from the first ionic core orbital:
c         do nic=2,nicm_st
         do nic=1,nicm_st
            ic = ncore_lsp(nic,il_st,is_st,ip_st)
c            if(ko(ic).gt.nk_ion_core(lo(ic))) then

c     Here we find how many f.c. orbitals are built on the core orbital nic.
c     nlctmax is the number and ltmp(nlct) is the ang.mom. of these orbitals.
               nlct = 0
               do j=1,lomax+1
                  ltmp(j) = -1
               end do
               do nsp=1,nspm
                  if(C(N,ic,nsp).ne.0d0) then
                     icheck1=0
                     do nicx=1,nicm_st
                        if(ncore_lsp(nicx,il_st,is_st,ip_st).eq.nsp) 
     >                     icheck1=1
                     end do
                     if(icheck1.eq.0) then                  
                        icheck = 1
                        do jl=1,max(nlct,1)
                           if(ltmp(jl).eq.lo(nsp)) icheck = 0
                        end do
                        if(icheck.eq.1) then
                           nlct = nlct + 1
                           ltmp(nlct) = lo(nsp)
c     print*, N,nic,nsp,lo(nsp),nlct,' * ',
c     >                       (ltmp(i), i=1,nlct)
                        end if                  
                     end if
                  end if
               end do
               nlctmax = nlct
c     print*, N,nic,ic, nlctmax
c     if nlctmax stays  =0 then the core ic=ncore(ic) was not used
c     for the state N.               
               do nlct=1,nlctmax
                  nspN = nspN + 1
                  if(nspN.gt.nspmax) then
                     print*,'routine optim_spf(): nspN.gt.nspmax'
                     stop
                  end if
                  lo(nspN) = ltmp(nlct)
                  ko(nspN) = 0  ! set to zero, this is essentually initialization
                  ngivesym(nspN) = 1000*il_st + 100*is_st + 10*(ip_st+1)
                  ngive_ion_orb(nspN) = ic                                    
                  ltest = ltmp(nlct)
                  do nsp=1,nspm
                     if(ltest.eq.lo(nsp).and.C(N,ic,nsp).ne.0d0) then
c     This makes 'f.c.o. nspN based on the core orbital nic' orthogonal 
c     to all core orbitals involved for given state N. 
c     We assume that s.p. basis is  orthogonomal.
                        icheck1=0
                        do nicx=1,nicm_st
                           if(ncore_lsp(nicx,il_st,is_st,ip_st).eq.nsp) 
     >                        icheck1=1
                        end do
                        if(icheck1.eq.0) A(nsp,nspN) = C(N,ic,nsp)
                     end if
                  end do
               end do           ! end nlct loop
c            end if
         end do                 ! end nic loop
      end do                    ! end N loop
      nspNm = nspN
c     Number of the opt.orbitals are: 
      N_opt = nspNm - nspm
c     
c     Array  A(nsp,nspN) gives representation of the new  basis
c     in terms of the original basis. Form new basis:
c     Just to have the same representation for all s.p. functions we set 
c     array A(nsp,nspN) also for the old s.p.basis.
      do n=1,nspm
         A(n,n) = 1d0
      end do

      do i=1,nr
         do nspN=nspm+1,nspm+N_opt
            sum = 0.0D0
            do  nsp=1,nspm
               if(i.ge.minf(nsp).and.i.le.maxf(nsp))
     >              sum = sum + A(nsp,nspN)*dble(fl(i,nsp))
            end do
            fl(i,nspN) = sum
         end do
      end do                    ! end i loop
      do nspN=nspm+1,nspm+N_opt
         call minmaxi(fl(1,nspN),nr,i1,i2)
         maxf(nspN) = i2
         minf(nspN) = i1
      end do
c
c     Choose the sign of the opt. functions - to be positive at start.
c     The sign of the second (i=minf(nspN)+1) after the starting minf(nspN)
c     r-grid point is used to determine the sign.
      do nspN=1,nspm
         sn_fl(nspN) = 1.0
      end do
      do nspN=nspm+1,nspm+N_opt
         sn_fl(nspN) = 1.0
            if(fl(minf(nspN)+1,nspN) .lt. 0.0) then
               do i=minf(nspN),maxf(nspN)
                  fl(i,nspN) = -fl(i,nspN)
                  sn_fl(nspN) = -1.0
               end do
c               print*,nspN, lo(nspN)
            end if
         end do
c
c     form new array e1r(nsp1,nsp2) and ortint(nsp1,nsp2)
      do nspN1=1,nspNm
         do nspN2=max(nspm+1,nspN1),nspNm
            sum = 0d0
            sumort = 0d0
            if(lo(nspN1).eq.lo(nspN2)) then
               do nsp1=1,nspm
                  do nsp2=1,nspm
                     sum = sum + A(nsp1,nspN1)*
     >                    A(nsp2,nspN2)*e1r(nsp1,nsp2)
                     sumort = sumort + A(nsp1,nspN1)*
     >                    A(nsp2,nspN2)*ortint(nsp1,nsp2) 
                  end do
               end do
            end if
            e1r(nspN1,nspN2) = sum * sn_fl(nspN1)* sn_fl(nspN2)
            e1r(nspN2,nspN1) = sum * sn_fl(nspN1)* sn_fl(nspN2)
            ortint(nspN1,nspN2) = sumort * sn_fl(nspN1)* sn_fl(nspN2)
            ortint(nspN2,nspN1) = sumort * sn_fl(nspN1)* sn_fl(nspN2)
         end do
      end do
c
c     Optimized orbitals are not normalized to unity. 
c     This can cause loss of accuracy in the consequent CI calculations. 
c     We normalize here opt.orb. to unity.
      do nsp=nspm+1,nspm+N_opt
         tmp = dsqrt(ortint(nsp,nsp))
         do i=minf(nsp),maxf(nsp)
            fl(i,nsp) = dble(fl(i,nsp))/tmp
         end do
         do n2=1,nspm
            if(lo(nsp).eq.lo(n2)) then
               ortint(n2,nsp) = ortint(n2,nsp)/tmp
               ortint(nsp,n2) = ortint(n2,nsp)
               e1r(n2,nsp) = e1r(n2,nsp)/tmp
               e1r(nsp,n2) = e1r(n2,nsp)
            end if
         end do
         do n2=nsp,nspm+N_opt
            if(lo(nsp).eq.lo(n2)) then
               tmp2 = dsqrt(ortint(n2,n2))
               ortint(n2,nsp) = ortint(n2,nsp)/tmp/tmp2
               ortint(nsp,n2) = ortint(n2,nsp)
               e1r(n2,nsp) = e1r(n2,nsp)/tmp/tmp2
               e1r(nsp,n2) = e1r(n2,nsp)
            end if
         end do
      end do
c--------------------------------------------------------------------------------
c     Here we exclude opt. orbitals that are too similar to avoid linear dependancy
c     find maximum angular momentum in the s.p. basis:
      lom = 0
      do nsp=1,nspm
         lom = max(lom,lo(nsp))
      end do
c     find maximum angular momentum in the target states:
      lmst = 0
      do n=1,Nmax
         lmst = max(lmst,ll(n))
      end do      
c     should be done separately for each symmetry
      do il=0,lmst
         do is=0,1
            do ip = -1,1,2
               Nsym = 1000*il + 100*is + 10*(ip+1)

               do l=0,lom
                  do n=nspm+1,nspm+N_opt
                     if(l.eq.lo(n).and.Nsym.eq.ngivesym(n)) then
                        ic = ngive_ion_orb(n)
c     Test overlaps between optimized orbitals of the same symmetry.
c     If the abs(overlap(n,nv)) is greater then 0.94 then these two optimised
c     orbitals are too similar and one of them should be excluded.
c     To identify such excluded orbitals in other parts of the code
c     set ko(n) = -1 for such orbitals.
c     The value of ko(n) for all 'good' opt.orbitals will be assigned
c     in the routine  set_basis_optim().
                        itest = 0
c     Set the maximum allowed value for overlaps: overlap_max, when orbitals is 
c     considered too similar to be included - if overlap=ortint(n,nv) > overlap_max 
c     then one of the orbitals must be excluded to avoid linear dependence problems.
                        overlap_max = 0.94
                        do nv=nspm+1,n-1
                           if(l.eq.lo(nv).and.Nsym.eq.ngivesym(nv).and.
     >                        ko(nv).ne.-1.and.ic.eq.ngive_ion_orb(nv)
     >                        .and.abs(ortint(n,nv)).gt.overlap_max)then
                              itest = 1
                              goto 20
                           end if
                        end do
 20                     if(itest.eq.1) then
                           ko(n) = -1
                           write(4,'("Excluded orbital: l=",i3,
     >                        ", Nsym=",i5,", ic=",i3,
     >                        ", ortint(n,nv)=",E10.4,
     >                        ", n,nv:",2i5)')
     >                        l, Nsym, ic, ortint(n,nv), n, nv
                        end if
                     end if
                  end do
               end do           ! end l loop
               
            end do
         end do
      end do                    !end il loop
c      

c     Exclude here all opt.orbitals that have ko(n) = -1
c     Initialize array
      do n=1,nspmax
         nA(n) = 0
      end do
c     All old basis functions will not be changed
      do n=1,nspm
         nA(n) = n
      end do
c     Here exclude opt.orbitals with ko(nspN) = -1
c     new basis --> nsp;  old basis --> nspN
      nsp = nspm
      do nspN=nspm+1,nspm+N_opt
         if(ko(nspN).ne.-1) then
            nsp = nsp + 1
            nA(nsp) = nspN
            if(nsp.ne.nspN) then
               lo(nsp) = lo(nspN)
               ko(nsp) = ko(nspN) ! should be zero
               ngivesym(nsp) = ngivesym(nspN)
               ngive_ion_orb(nsp) = ngive_ion_orb(nspN)
               do i=1,nr
                  fl(i,nsp) =  fl(i,nspN)
               end do
               maxf(nsp) = maxf(nspN)
               minf(nsp) = minf(nspN)
            end if
         end if
      end do
c     New total number of s.p.functions are
      nspNm = nsp
         
c     Number of the opt.orbitals are: 
      N_opt = nspNm - nspm
      if(nspNm.gt.nspmCI) then
         print*,'routine optim_spf(): nspNm.gt.nspmCI'
         print*,'nspNm=',nspNm, ', nspmCI=',nspmCI
         stop
      end if
c     form new array e1r(nsp1,nsp2) and ortint(nsp1,nsp2)
      do nspN1=1,nspNm
         do nspN2=max(nspm+1,nspN1),nspNm
            sum = 0d0
            sumort = 0d0
            if(lo(nspN1).eq.lo(nspN2)) then
               nsp1= nA(nspN1)
               nsp2= nA(nspN2)
               sum = e1r(nsp1,nsp2)
               sumort = ortint(nsp1,nsp2) 
            end if
            e1r(nspN1,nspN2) = sum
            e1r(nspN2,nspN1) = sum
            ortint(nspN1,nspN2) = sumort
            ortint(nspN2,nspN1) = sumort
         end do
      end do     
c-------------------------------------------------------------------------------- 
      write(4,'("Number of the optimized orbitals: N_opt =",I4)') N_opt  

c     reset to zero all necessary arrays:
      do nic=1,nicm
         ncore(nic) = 0
      end do
      nicm = 0
      do il=0,lmst
         do is=0,1
            do ip = -1,1,2
               do n=1,nicm_lsp(il,is,ip)
                  ncore_lsp(n,il,is,ip) = 0
               end do
               nicm_lsp(il,is,ip) = 0
            end do
         end do
      end do
c     array C(N,n1,n2) = 0d0 is set to zero in start_structure(...).

      nspm = nspNm

c$$$          do n1=1,nspm
c$$$             do n2=n1,nspm
c$$$                if(lo(n1).eq.lo(n2)) then
c$$$                   write(20,*) n1,n2, lo(n1),lo(n2),
c$$$     >                ortint(n1,n2), e1r(n1,n2)
c$$$                end if
c$$$             end do
c$$$          end do

c     call  G - S orthogonalization of the new s.p. basis
c      call gsort(nspm,fl,maxf,minf,ortint,e1r,lo)
c     

      deallocate(A)
      deallocate(nA)
      deallocate(sn_fl)

      return
      end
c------------------------------------------------------------------------
c     set optimized s.p. basis for the given symmetry. 
c     Called from file hel-3.f, routine  hel(...).
      subroutine set_basis_optim(il,is,ip,nk1)
      include 'par.f'
      dimension  nk1(0:lomax)   ! this array is passed here to exclude only those orbitals from old basis which will be used for given symmetry.
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /ngive_sym_opt/  ngivesym(nspmax), N_opt
      common /ortog/  ortint(nspmax,nspmax)
      double precision ortint, tmp
      common /ngive_ionic_orb/  ngive_ion_orb(nspmax)
      common /include_opt/  include_opt(nspmCI,nspmCI)
      common /ionic_configurations/ l_ion_core, nk_ion_core(0:lomax)


c     Return if N_opt=0, this is all cases except when the optimised orbitals 
c     have been built. For example programm returns from here when it goes 
c     first time to do  CI calculations to obtain states that 
c     have been requested to be optimised.     
      if(N_opt.eq.0) return
      Nsym_current = 1000*il + 100*is + 10*(ip+1)


c     Find maximum angular momentum in the s.p. basis:
      lom = 0
      do nsp=1,nspm
         lom = max(lom,lo(nsp))
      end do

      do nsp=1,nspm
         nset(nsp) = 1
      end do
c     
c     Array include_opt(n,ic) is used in routine config(...) to decide which
c     s.p.orbitals should be used in building configurations.
c     The second index of this array is the ionic core orbital index. For each ionic
c     core orbital we specify which s.p.orb. are allowed to be used for building
c     of configurations. If orbital nsp may not be used to build a configuration
c     with the ionic core orb. ic then include_opt(nsp,ic) = -1.
c     Any orbitals from the old basis can be used to form configuration with each other.
      do n2=1,nspm-N_opt
         do n1=1,nspm
            include_opt(n1,n2) = 0
         end do
      end do
c     Not all opt.orbitals are allowed to form configurations with old basis orbitals.
c     No configurations constructed only of opt.orbitals are allowed.
      do n2=nspm-N_opt+1,nspm
         do n1=1,nspm
            include_opt(n1,n2) = -1
         end do
      end do

c     This is a temporary coding. It assumes that CI expansion involving ionic core
c     orbitals (for which continuum is biult) comprise all s.p.orbitals from 
c     the old basis. In this case there no need to build opt.orbitals based on such
c     ionic core orbitals. This may be not the case if CI expansion involving ionic core
c     orbitals do not use all old s.p. orbitals or a new basis is built (to discretize
c     continuum with larger spacing between positive energy states).
C     Comment this codding
c$$$      do nsp=nspm-N_opt+1,nspm
c$$$         ic = ngive_ion_orb(nsp)
c$$$         if(l_ion_core.ge.lo(ic).and.nk_ion_core(lo(ic))
c$$$     >        .ge.ko(ic)) then
c$$$             ko(nsp) = -1
c$$$         end if
c$$$      end do



c     We go through all opt.orbital and find all corresponding ionic core orbitals.
c     We set  include_opt(ic1,ic2) for these orbitals to avoid excluding them when testing
c     on overlap with optimized orbitals (this can happen if in file F5 setting
c     for CI array nk1(l1) has been chosen to be very small).
c     We will restore nset(ic)=1  values in the end of the routine.
c     Go over optimized orbitals picking only those orbitals that correspond
c     to the current symmetry (il,is,ip).
c     Remember that opt.orbital n should be excluded if ko(n) = -1.
      do nsp=nspm - N_opt + 1,nspm
         if(Nsym_current.eq.ngivesym(nsp).and.ko(nsp).ne.-1) then
            ic2 = ngive_ion_orb(nsp)
            do nn=nspm - N_opt + 1,nspm
               if(Nsym_current.eq.ngivesym(nn).and.ko(nn).ne.-1) then
                  ic1 = ngive_ion_orb(nn)
                  include_opt(ic1,ic2) = -1
                  include_opt(ic2,ic1) = -1
               end if
            end do
         end if
      end do
c
c     Go over optimized orbitals picking only those orbitals that correspond
c     to the current symmetry (il,is,ip).
c     Remember that opt.orbital n should be excluded if ko(n) = -1.
      do nsp=nspm - N_opt + 1,nspm
         if(Nsym_current.eq.ngivesym(nsp).and.ko(nsp).ne.-1) then
            nset(nsp) = 2
c     find ionic-core orbital (ic) corresponding to the current optim. orbital (nsp)
c     and exclude all configurations except of (nsp,ic).
            ic = ngive_ion_orb(nsp)
            do jc=1,nspm
               include_opt(nsp,jc) = -1
            end do
            include_opt(nsp,ic) = 0
c     
c     find orbital (from the old set) that has largest overlap with 
c     current opt.orb.  and exclude former orbital. Ionic core orbital are 
c     orthogonal to the opt. orb. by construction of the latter ones therefore
c     the starting from nn=1 is safe, but we still exclude them by setting
c     include_opt(nn,ic) = -1 for ionic core orbitals in coding above.
            tmp = 0d0
            nn_find = -1
            do nn=1,nspm - N_opt
               if(lo(nn).eq.lo(nsp).and.ko(nn).le.nk1(lo(nsp))
     >            .and.include_opt(nn,ic).ne.-1) then
                  if(dabs(ortint(nn,nsp)).gt.tmp) then
                     tmp = dabs(ortint(nn,nsp))
                     nn_find = nn 
                  end if
               end if
            end do
            if(nn_find.ne.-1) then
               nset(nn_find) = -1
               include_opt(nn_find,ic) = -1
               ko(nsp) = ko(nn_find)
            else
               write(4,'("ATTENTION: Have not find old orbital",
     >            " to form overlap, symmetry: l,s,p:",3I4,
     >            ", orbital: n,l,ic:",3I4)')
     >            il, is, ip, nsp, lo(nsp), ic
               write(*,'("ATTENTION: Have not find old orbital",
     >            " to form overlap, symmetry: l,s,p:",3I4,
     >            ", orbital: n,l,ic:",3I4)')
     >            il, is, ip, nsp, lo(nsp), ic
            end if                        
         else
c     The optim. orbital  nsp does not correspond to the current symmetry and shell be excluded
            nset(nsp) = -1
            do jc=1,nspm
               include_opt(nsp,jc) = -1
            end do
         end if
      end do
      
c     Restore include_opt(ic1,ic2)=1  values for the ionic core orbitals.
      do nsp=nspm - N_opt + 1,nspm
         if(Nsym_current.eq.ngivesym(nsp).and.ko(nsp).ne.-1) then
            ic2 = ngive_ion_orb(nsp)
            do nn=nspm - N_opt + 1,nspm
               if(Nsym_current.eq.ngivesym(nn).and.ko(nn).ne.-1) then
                  ic1 = ngive_ion_orb(nn)
                  include_opt(ic1,ic2) = 1
                  include_opt(ic2,ic1) = 1
               end if
            end do
         end if
      end do
c      
      do n=nspm-N_opt+1,nspm
         Nsym = ngivesym(n)
         if(Nsym.eq.Nsym_current) then
            write(4,'(I5,", l s par:",3I3,
     >         ", n,lo(n),ko(n),ion_core: ",4I5)') 
     >         Nsym,il,is,ip,n,lo(n),ko(n), ngive_ion_orb(n)
         end if
      end do
         
      
      return 
      end
c------------------------------------------------------------------------
c     These routines deal with predetermined state ordering based on the value 
c     of ICS (Born). State ordering is recorded in the file Set.order_of_states

C     This routine is for given target symmetry (fixed il,is,ip)
      subroutine check_state(iorder_st,Max_st_sym,num,jpartN)
      include 'par.f'
      integer iorder_st(komax) ! used for setting state ordering based on ICS (Born)

      jpartN = 0
      do k=1,Max_st_sym
         if(num.eq.iorder_st(k)) then
            jpartN = 1
            goto 10
         end if
      end do
            
 10   return
      end
      
C     This routine reads array iorder_st(k,il,is,ip) from file Set.order_of_states.
      subroutine read_Set_order_of_states(nstate,iorder_st)
      include 'par.f'
      integer nstate(0:lomax,0:1,2)
      integer iorder_st(komax,0:lamax,0:1,2) ! used for setting state ordering based on ICS (Born)

      open(135,file='Set.order_of_states')

      lam = -1
      
 5    read(135, *, ERR = 20, END = 10) il,is,ip,Max

      ip = ip + 1

      lam = max0(lam,il)
      
      if(Max.lt.nstate(il,is,ip)) then
         print*,'Inconsistent input, (Max.lt.nstate(il,is,ip), in ',
     >      'file Set.order_of_states less states are present than ',
     >      'it was ordered in F5 file.'
         print*,'Program terminated'
         stop
      end if

      read(135, *, ERR = 25, END = 10)
     >   (iorder_st(k,il,is,ip), k=1,nstate(il,is,ip))

      goto 5


 20   print*,'20: Error when reading file Set.order_of_states ',
     >   'in routine read_Set.order_of_states(nstate,iorder_st), ',
     >   'file structure.f'
      print*,'Terminating program'
      stop
      
 25   print*,'25: Error when reading file Set.order_of_states ',
     >   'in routine read_Set.order_of_states(nstate,iorder_st) ,',
     >   'file structure.f'
      print*,'Terminating program'
      stop
      
 10   do il=0,lam
         do is=0,1
            do ip=1,2
               write(4,*) il,is,ip
               write(4,'(50I3)') (iorder_st(k,il,is,ip),
     >            k=1,nstate(il,is,ip))
            end do
         end do
      end do

         
      close(135)

      
      return
      end

c------------------------------------------------------------------------
c     This routine determine the sign of the wave-function by summing
c     all CI coefficients. The sign is positive if the sum is positive and
c     negative if  the sum is negative. In the last case all CI coefficients will
c     be multiplied by (-1.0).
      subroutine resolve_sign(Nmax,N,ncm,inum,CI,sign_wf)
      include 'par.f'
      
      double precision sign_wf, CI(ncmCI), tmp, small, CI_sum
      character chan(knm)*3
      common /charchan/ chan

      small = 0.001

      CI_sum = 0.0

      do i=1,ncm
         CI_sum = CI_sum + CI(i)
      end do

      if(CI_sum .lt. 0.0) then
         sign_wf = -1.
      else
         sign_wf = 1.
      end if
      
      return
      end
c------------------------------------------------------------------------
c     Resolve sign of target states - might change the sign of a target state
c     only in the case when nonorthogonal basis was used for diagonalization
c     This routine determine the sign of the wave-function by summing
c     all CI coefficients. The sign is positive if the sum is positive and
c     negative if  the sum is negative. In the last case all CI coefficients will
c     be multiplied by (-1.0).
      subroutine  resolve_sign_final(Nmax,nspmW,C)
      include 'par.f'
      
      double precision C(Nmax+1,nspmW,nspmW), ortint
      double precision   tmp, tmp1, small, CI_sum(KNM)
      real  sign_wf
      character chan(knm)*3, cnode*3, ench*11
      common /MPI_info/myid, ntasks, cnode, ench

      common /charchan/ chan
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common /meshrr/ nr, gridr(nmaxr,3)

      
c      print*,'***** Nmax=',Nmax
c      print*,'***** C(1,1,1)',C(1,1,1)

      small = 0.001

      do N=1,Nmax
c         print*, '*** N,nam(N):',N,nam(N)
         
         CI_sum(N) = 0.0
         
         do i=1,nam(N)
            do j=1,nam(N)
c               print*, '******** i,j,C(N,i,j):', i,j,C(N,i,j)
               tmp = C(N,i,j)
               n1 = na(i,N)
               n2 = na(j,N)
               if(n2.le.n1.and.tmp.ne.0.0D0) then                 
                  if(n1.ne.n2) tmp = tmp * dsqrt(2d0)
                  tmp1 = tmp*dsqrt(ortint(n1,n1)*ortint(n2,n2))
                  CI_sum(N) = CI_sum(N) + tmp1
               end if
            end do
         end do
         
         if(CI_sum(N) .lt. 0.0) then
c            print*, N, chan(N)
            sign_wf = -1.
         else
            sign_wf = 1.
         end if
c         print*,'***** N=',N,CI_sum(N)
         if(dabs(CI_sum(N)) .lt. small) then
            print*,'Attention: for state ', chan(N), ' N=',N,
     >         ' there can be a problem with the ',
     >         'sign of the wave function: CI coefficient sum ',
     >         'is too small in the absolute value:', CI_sum(N)
            print*,'nam(N)=', nam(N)
            do i=1,nam(N)
               do j=1,nam(N)
                  n1 = na(i,N)
                  n2 = na(j,N)
                  print*, ' n1, n2, C(N,i,j):', n1, n2, C(N,i,j)
               end do
            end do
         end if
c         print*,'*****111 N=',N,sign_wf
         if(sign_wf .lt. 0.0) then
c         print*,'***** N=',N,sign_wf
            do j=1,nam(N)
               do i=1,nam(N)
                  C(N,i,j) =  -C(N,i,j)
               end do
            end do            
         end if
         
      end do                    ! end  N  loop
      print*,'!!!!!!'
      open(534,file=adjustl(cnode//'save.sign.wf'//ench))
      write(534,'("lable         N   WaveFunc(i1,i2)        CI_sum")')
      do N=1,Nmax
         tmp1 = 0.0
         do i=1,nam(N)
            do j=1,nam(N)
               tmp = C(N,i,j)
               n1 = na(i,N)
               n2 = na(j,N)
               i1 = minf(n1) + 1
               i2 = minf(n2) + 1
               tmp1 = tmp1 + tmp*dble(fl(i1,n1)*fl(i2,n2))
            end do
         end do
c         write(534,'(A5,I5,2E12.4)') chan(N), N, tmp1, CI_sum(N)
c         write(534,*) chan(N), N, tmp1, CI_sum(N)
      end do                    ! end  N  loop
      close(534)
      print*,'-------'
      return
      end

c-----------------------------------------------------------------------
      subroutine writesmat(lg,ispin,ipar,nchtop,nchan,smat,maxr,
     >     energy)
      complex smat(nchan,nchan)
      dimension temp(maxr)
      character*20  stren
      
      write(stren,'(1p,E10.4)') energy

      nn = 175
      if(lg .eq. 0 .and. ispin .eq. 0) then
         open(nn,file='smat_out_'//stren)         
      endif

      tmp = 0.0
      write(nn,'("J S Par nchtop:",4i5)') lg, ispin, ipar, nchtop 
      do nchi = 1, nchtop
         call getchinfo (nchi,nchip,lg,temp,maxpsi,ei,lia,nia,li)
         tmp = 0.0
         do nchf = 1, nchtop
            call getchinfo (nchf,nchp,lg,temp,maxpsi,ef,lfa,nfa,lf)
            tmp = tmp + abs(smat(nchf,nchi))**2
            write(nn,'(6i4,2X,1p,3(1X,E12.5))') li,lia,nia,lfa,lf,
     >           nfa,abs(smat(nchf,nchi))**2,smat(nchf,nchi)

         enddo
      enddo

      return
      end



c-----------------------------------------------------------------------
C$$  Two-electron state is a product of two-electron angular function: |l1 l2 : L_st S_st Par_st > 
C$$  and corresponding radial function fst(r1,r2 | l1,l2) summed of all l1, l2
C$$
C$$    WF(r1,r2) = \sum_{l1,l2} fst(r1,r2 | l1,l2) *  |l1 l2 : L_st S_st Par_st >
C$$
C$$   NOTE: the terms like these: l1=1, l2=2   and  l1=2, l2=1
C$$   enter the sum  separately but they can be compbined to a single term using symmetry of
C$$   CI two-electron angular function:
C$$   |l1 l2 : L_st S_st Par_st >  = (-1)**(l1 + l2 - L_st + 1 - S_st)  * |l2 l1 : L_st S_st Par_st > 
      subroutine TwoElectWFValue(nst,i1,i2, res)
!      subroutine TwoElectWFValue(nst,i1,i2)

      use CI_MODULE      

      include 'par.f'
      
      integer nst    ! two-electron state number
      integer i1, i2   ! indexes to rgrid 
      real res(0:lamax,0:lamax)  ! array of values 
      
      integer la, sa, lpar, np
      common /helium/ la(KNM), sa(KNM), lpar(KNM), np(KNM)
      double precision E
      common /CcoefsE/  E(KNM)

      double precision ortint
      common /ortog/  ortint(nspmax,nspmax)
      integer nspm, lo, ko, nset
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      real fl
      common /funLag/  fl(nmaxr,nspmax)
      integer maxf, minf
      common /minmaxf/ maxf(nspmax),minf(nspmax)

      integer na, nam
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer nicm, ncore
      common /corearray/ nicm, ncore(nspmCI)      

      double precision E_st
      integer L_st, S_st, Par_st
      integer lomax_st
      
      integer  jn1, jn2, n1, n2, l1, l2

      lomax_st = -1
      res(:,:) = 0

C$$   two-electron state number  nst  description
      E_st = E(nst)  ! two-electron energy
      L_st = la(nst) ! orb. ang.  mom.
      S_st = sa(nst) ! spin
      Par_st = lpar(nst)  ! parity


      do jn1=1,nam(nst)
         n1 = na(jn1,nst)
         l1=lo(n1)
         
         if(i1 .lt. minf(n1) .or. i1 .gt. maxf(n1)) cycle
         
         do jn2=1,nam(nst)
            n2 = na(jn2,nst)
            l2=lo(n2)
            
            if( C(nst,jn1,jn2).eq. 0) cycle

            if(i2 .lt. minf(n2) .or. i2 .gt. maxf(n2)) cycle
            print*, l1,l2, n1, n2, fl(i1,n1), fl(i2,n2)!, C(nst,jn1,jn2)
            res(l1,l2) = res(l1,l2) + C(nst,jn1,jn2) * 
     >           fl(i1,n1) * fl(i2,n2)
            lomax_st = max( lomax_st,l1,l2)

         enddo
      enddo

      print*
      print*, '***>>  Test TwoElectWFValue(...)'
      print*, '*** >> lomax_st =',  lomax_st, lamax

      do l1=0,lomax_st
         do l2=0,lomax_st

            print*, '***>>', l1, l2, res(l1,l2)
            
         enddo
      enddo


      return
      end

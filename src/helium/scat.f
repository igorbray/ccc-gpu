      subroutine scattering(myid,ifirst,theta,nold,etotal,KJ,gridk,
     >   enionry,npk,vdcore,dwpot,nchm,Nmax,namax,
     >   nze,td,te1,te2,te3,vdon,
!     >   vmat,
     >   nsmax,itail,phasel)
      use CI_MODULE
      use chil_module
      use po_module
      use DM_MODULE
      use vmat_module
      use date_time_module
#ifdef GPU
      use openacc
#endif
c$$$      use mpi_f08 !avoids MPI argument missmatch warnings
      include 'mpif.h'
      
      include 'par.f'
      include 'par.pos'
c      implicit real*8 (a-h,o-z)       
      parameter (maxl=2*ltmax+1)
      integer my_status(MPI_STATUS_SIZE)
c$$$      type(MPI_Status) :: my_status
c$$$      type(MPI_Request) :: my_request
      integer sa, npk(nchm+1)
      common /helium/ la(KNM), sa(KNM), lpar(KNM), np(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      real*8, allocatable :: fli(:,:) ! fli(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common /CcoefsE/  E(KNM)
      double precision E, ortint, Etot, Etot_nu, Etot_ve
      common /meshrr/ nr, gridr(nmaxr,3)
c$$$      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
c$$$     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
c$$$     >   npbot(0:lamax),nptop(0:lamax)
      common /increarrange/ inc
      dimension minchiform(kmax),maxchiform(kmax)!,chiform(maxr,kmax)
      allocatable :: chiform(:,:),orts(:,:),fls(:,:),ortr(:,:),flr(:,:)
c$$$      dimension chil(nr,npk(nchm+1)-1),minc(npk(nchm+1)-1)
      dimension gridk(kmax,nchan), df(maxr), dwpot(maxr,nchan),
     >   vdcore(maxr,0:lamax),
     >   vdon(nchan,nchan,0:1), !vmat(npk(nchm+1)-1,npk(nchm+1)),
     >   temp(maxr),psii(maxr),psif(maxr)
c$$$      real  ortchil(npk(nchm+1)-1,nspm),flchil(npk(nchm+1)-1,nspm)
      real, allocatable :: ortchil(:,:),flchil(:,:)
      integer na, nam
      logical posi, posf, positron,ex
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      common /corearray/ nicm, ncore(nspmCI)      
      real, dimension(:,:), allocatable ::  ortchil_po
      real:: theta, theta_ve, vmatt(kmax,kmax,0:1)
      common /ionic_configurations/ l_ion_core, nk_ion_core(0:lomax)
      integer, dimension(:), allocatable :: is_core_orb  
      common /noblgegas_switch/ i_sw_ng  
crr -added by Rav
      real*8 xgz,wgz,pl,xgp,wgp,etotp 
      COMMON/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      common/gausp/xgp(igpm),wgp(igpm),igp
      real*8 Qlp(kmax,2*kmax*igp),Q0p(kmax,2*kmax*igp),
     >   Alp(kmax,2*kmax*igp)
      real*8 aq(0:kmax),qq(kmax),xp(2*kmax*igp+kmax)
     >                ,wp(2*kmax*igp)
      real*8, allocatable ::Qlpi(:,:,:),Q0pi(:,:,:),Alpi(:,:,:),
     >   aqi(:,:),xpi(:,:),wpi(:,:)
      integer*4 ionshi(nchm),nqgeni(nchm),nqimi(nchm),imaxi(nchm)
      real*8 ve2stat(maxr),va2stat(maxr)
      real*8, allocatable :: fd0(:,:,:,:)
      real  fd01(maxr,0:5) 
      real*8 gp(iglp),
     >   Bnlp(iglp,Nmax),Dnlp(iglp,Nmax),Unlp(iglp,Nmax),Amg(iglp)
      dimension nchps(nchm), nchat(nchm)
      complex phasel(kmax,nchan)
      integer, allocatable :: nchansi(:), nchansf(:) 
      real gridp(nmaxr,3)
      real*8 ya(10),xa(10),yr,ry,dy,formf0
      character*9 nodetfile
      integer ngpus, gpunum, omp_get_thread_num
      external omp_get_thread_num

        if (npk(2)-npk(1).eq.1) then
         nchii = 1
         nchif = nchm
      else
         nchii = nchistart(nodeid)
         nchif = nchistop(nodeid) 
      endif
c$$$      if (KJ.lt.10) then
c$$$         write(nodetfile,'(i3,"_",i1,"_",i1)') nodeid,KJ,ipar
c$$$      elseif (KJ.lt.100) then
c$$$         write(nodetfile,'(i3,"_",i2,"_",i1)') nodeid,KJ,ipar
c$$$      else
c$$$         write(nodetfile,'(i3,"_",i3,"_",i1)') nodeid,KJ,ipar
c$$$      endif
      
#ifdef GPU
      ngpus= max(1,acc_get_num_devices(acc_device_nvidia))
      gpunum = mod(myid,ngpus)
      call acc_set_device_num(gpunum,acc_device_nvidia)
      print*,'MYID associated with GPU:',myid,
     >     acc_get_device_num(acc_device_nvidia) 
!$acc enter data copyin(chil(1:nr,npkstart:npkstop,1))
!$acc& copyin(nchm,npk(1:nchm+1),minchil(npkstart:npkstop,1))
#endif 
 
C  One electron continuum of Helium assumes the target He+ is left in
C  the ground state. The ENIONRY is added in the mainhe routine, where the
C  energies are defined. Etot is in a.u.
      Etotp=etotal*0.5d0  ! 
      Etot = (etotal + enionry) / 2.0
      Etot_ve = Etot
      theta_ve = theta
!       PRINT*, nze,etotal, enionry
CRRRRRR - Rav's add
      npsch=0 
      nchatom=0
      nchps(:)=0
      nchat(:)=0

c$$$      if (npk(2)-npk(1).gt.1) open(42,file=nodetfile)

      IF(nze.eq.1) THEN
         nqm1= npk(1+1) - npk(1)
         nqmm=npk(nchm+1)-npk(nchm)

         IF(nqmm*nqm1.gt.1)  then
            lpsmax=0
            maxrps=0
            do i=nchii,nchm 
            call getchinfo (i, nt, KJ, psii, maxpsit, eit, lit, nit, Li)
            posf = positron(nit,lit,npost)
!        PRINT*, 'nch,Ps:: lp1,nion,lion', i,KJ,lit,posf, '::',lp1,lion,Li   
            if(posf) then
               npsch=npsch+1    ! max number of Ps-channels
               nchps(i)=npsch
!          PRINT*,'Ps: ei',eit
               if( maxrps.lt.maxpsit) maxrps=maxpsit
               if(lpsmax.lt.lit) lpsmax=lit ! max of Ps orb quant number
            elseif((.not.posf).and.i.le.nchif) then
               nchatom=nchatom+1
               nchat(i)=nchatom
!          PRINT*,'He: ei',eit
          ! nanchi(nchatom)=i
            endif         
            enddo

            IF(npsch.gt.0.or.nchatom.gt.0) THEN
               print*,'CAUTION: allocatable fli not yet tested'
               maxfliall = maxval(maxf)
               allocate(fli(maxfliall,nspm))
               DO ni1=1,nspm
                  maxfli=maxf(ni1)
                  do i=1,maxfli
                     fli(i,ni1)=dble(fl(i,ni1))
                  enddo
c$$$          fli(maxfli+1:nmaxr,ni1) = 0.d0
                  fli(maxfli+1:maxfliall,ni1) = 0.d0

                  maxflr=maxf(ni1)
                  r1=gridr(maxflr-50,1); r2=gridr(maxflr,1)
                  fl1=fl(maxflr-50,ni1); fl2=fl(maxflr,ni1)

                  alambda = 1./(r2-r1)*(Log(fl1/fl2)+ni1*Log(r2/r1))
c$$$          alambda=1./(r2-r1)*Log((fl1*r2**ni1)/(fl2*r1**ni1))
c$$$          print*,'alam1,alam:',alambda1,alambda
c$$$          print*,'nodeid,cfl1:',nodeid,fl1*exp(-ni1*log(r1)+alambda*r1)
c$$$  print*,'nodeid,cfl2:',nodeid,fl1*r1**(-ni1)*exp(alambda*r1)
c$$$          cfl=fl1*r1**(-ni1)*exp(alambda*r1)
                  cfl = fl1*exp(-ni1*log(r1)+alambda*r1)
!     WRITE(17,'(a1,3i8,2e20.10)')'#',ni1,maxflr,nr,alambda,cfl
!        WRITE(17,*) 

!        do i=1,nr
!           rr=gridr(i,1)
!           if(i.gt.maxflr) fli(i,ni1)=cfl*rr**(ni1)*exp(-alambda*rr)
!          WRITE(17,'(i5,3e20.10)') i,gridr(i,1),fl(i,ni1),fli(i,ni1)
!        enddo
               ENDDO
            ENDIF
      
            if(npsch.gt.0) then
               PRINT*,'number of Ps-channels:',npsch,'out of total',nchm
               PRINT*, 'nmaxr & lmax for Ps channels:', maxrps,lpsmax
               call makevstat(lo,fli,maxfliall,nspm,maxf,formf0,
     >            ve2stat,va2stat)

               if(nchatom.gt.0) then
                  Bnlp=0.; Dnlp=0.; Unlp=0.; Amg=0.;    
!!!!!!!!!!11
                  qmax =60.0 
                  ndoub=6
                  pmax = gridr(nr,1)        
                  call grids(qmax,ndoub,pmax,gridp,nmaxr,mpg,jdouble,
     >               40,njdouble) 
!!!!!!!!!!!!11
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& default(private)
C$OMP& shared(nze,Nmax,la,sa,lpar,nspm,lo,C,ortint,fli,maxf,minf)
C$OMP& shared(na,nam,namax,gp,Bnlp,Dnlp,Unlp,Amg,mpg,gridp)
                  DO Ni=1,Nmax
        
                     call makehepsff(nze,Nmax,Ni,la,sa,lpar,nspm,lo,
     >                  C,ortint,fli,maxfliall,nspm,maxf,minf,na,nam,
     >                  namax,gp,Bnlp,Dnlp,Unlp,Amg,mpg,gridp)

!        PRINT*,'He-pseudo FF:::', Bnlp(1,Ni),Unlp(1000,Ni),Amg(100)

                  ENDDO 
C$OMP END PARALLEL DO 
!         call clock(tff3) 
!           Print*, 'makehepsff time:', tff3-tff2
               endif            ! nchatom       
               if(.not.allocated(fd0)) then
                  allocate(fd0(maxr,0:lpsmax,npsch,npsch))
                  fd0 = 0.
                  PRINT*, 'fd0 is allocated'        
               endif 

               if(allocated(fd0)) then  
                  CALL hepsdir(gp,KJ,ve2stat,va2stat,nchm,npsch,nchps, 
     >               lpsmax,maxrps,nchii,nchif,formf0,fd0)
               endif

!          call clock(tff3) 
!          Print*, 'ps-ps fd0 cal time:',  tff3-tff2
            endif !npsch .gt.0
            if (allocated(fli)) deallocate(fli)
         ENDIF !nqmm*nqm1.gt.1
      ENDIF     !nze.eq.1
CRRR -End of Rav's add         
 
      if (ifirst.ne.0) then

c     make array which tells if a function is a core orbital (=1 or 2) or not (0), 
c     if it is a core orbital set value to 2 if this orbital is used for dicretization: l_ion_core, nk_ion_core(0:lomax)
         if (allocated(is_core_orb)) deallocate(is_core_orb) 
         allocate(is_core_orb(nspm))
         is_core_orb(:) = 0
         do nic=1,nicm
            nsp = ncore(nic)
            knsp = ko(nsp)
            lnsp = lo(nsp)
            if(l_ion_core .ge. 0) then
               is_core_orb(nsp) = 1
               if(lnsp .le. l_ion_core) then
                  if(knsp .le. nk_ion_core(lnsp)) then
                     is_core_orb(nsp) = 2
                  endif
               endif
               elseif(l_ion_core .lt. 0) then
                  is_core_orb(nsp) = 2
               endif
c     print*,'*** ion_core: ', nsp, is_core_orb(nsp), 
c     >           ko(nsp),nk_ion_core(lnsp)
         enddo
c

         if(abs(theta) .ge. 1000.0) then
            print*, 'scat.f: Setting up D matrix'
            Tspin = 0.5
            iacc = nint(abs(theta-1000.0))
            if(theta .le.-1000) then
               iacc = -1
            endif
            call nonuniq(KJ,Tspin,nchm,iacc)
         endif
!

         call date_and_time(date,time,zone,valuesin)
c$$$         print '("Making one-electron array at:",a10)',time
c$$$         call update(6)
         call clock(s1)
         
c     This block is to code the projection operator to fix nonuniqueness.
c     Note: only arrays lo_po(:) and ortchil_po(:,:) are required in vmat routine
c     in the part where nonuniqueness is dealt with.

         if(theta .ne. 0.0) then
            if (inc .eq. 0) then
               nspm_po = nspm
               if (allocated(lo_po)) deallocate(lo_po)                  
               allocate(lo_po(nspm_po))               
               if(nspm_po .ne. nspm) then
                  print*,'Please make sure that all s.p.orbitals ',
     >               'declared in ccc.in or F5 files are used in ',
     >               'CI calculations'
                  print*,'Can be changed later to account for more ',
     >               'compex case'
                  print*, 'scat.f: nspm_po .ne. nspm'
                  stop
               endif
               lo_po(1:nspm) = lo(1:nspm)
c     The aim is to copy ortchil(:,:)  array to  ortchil_po(:,:) array
c     after ortchil(:,:) array is calculated in routine ortchilnsp(...)
            else
c     ortchil_po(:,:) array has to be calculated directly
c     array lo_po(:) is already set up
            endif
            if (allocated(ortchil_po)) deallocate(ortchil_po)
            allocate(ortchil_po(npk(nchm+1)-1,nspm_po))
         else
            allocate(ortchil_po(1,1)) !to stop "unallocated" errors when used in arguments
         endif
         if (allocated(ortchil)) deallocate(ortchil)
         allocate(ortchil(npk(nchm+1)-1,nspm))
         if (allocated(flchil)) deallocate(flchil)
         allocate(flchil(npk(nchm+1)-1,nspm))
#ifdef _single
#define nbytes 4
#define MY_MPI_REAL MPI_REAL
#define MY_MPI_COMPLEX MPI_COMPLEX
#elif defined _double
#define nbytes 8
#define MY_MPI_REAL MPI_REAL8
!#define MY_MPI_COMPLEX MPI_COMPLEX16
!#define MY_MPI_REAL MPI_DOUBLE_PRECISION
#define MY_MPI_COMPLEX MPI_DOUBLE_COMPLEX
#endif
c     
c     
c     get overlap for projectile and s.p. functions
C  was from nchii to nchm, which is nchtop, now to nchif with MPI below         
         call ortchilnsp(KJ,npk,nchm,fl,maxf,
     >      minf,lo,nspm,vdcore,dwpot,inc,ortchil,flchil)
         call date_and_time(date,time,zone,valuesout)
         print '(i4,": nodeid, ortchilnsp call complete at: ",a10,
     >", diff (secs):",i5)',nodeid,time, idiff(valuesin,valuesout)
!MPI section begins         
! nodeid needs ortchil/flchil for nodeid to nodes, but calculates only for nodeid, and reads from nodeid+1 to nodes.
! nodeid sends ortchil/flchil: 1 to nodeid-1
         ntagort = 10
         ntagfl = 11
         npkstart = npk(nchii)
         npkstop = npk(nchif+1)-1
         nt = npkstop - npkstart + 1
         nsend = nt * nspm
         allocate(orts(nt,nspm),fls(nt,nspm))
         do nsp = 1, nspm
            orts(1:nt,nsp) =
     >         ortchil(npkstart:npkstop,nsp)
            fls(1:nt,nsp) =
     >         flchil(npkstart:npkstop,nsp)
         enddo
         do nn = 1, nodeid - 1
c$$$            print*,'Will send for nodeid to n:',nodeid,nn,npkstart,
c$$$     >         npkstop
            call MPI_ISEND(orts,nsend, 
     >         MY_MPI_REAL,nn-1,ntagort,MPI_COMM_WORLD,my_request,ierr)
            call MPI_ISEND(fls,nsend, 
     >         MY_MPI_REAL,nn-1,ntagfl,MPI_COMM_WORLD,my_request,ierr)
c$$$            print*,'Sent for nodeid to n:',nodeid,nn
         enddo
         do nn = nodeid+1, nodes
            nt = npk(nchistop(nn)+1) - npk(nchistart(nn))
            allocate(ortr(nt,nspm),flr(nt,nspm))
            nrcv = nt * nspm
c$$$            print*,'Will receive on nodeid from n:',nodeid,nn,
c$$$     >         npk(nchistart(nn)),npk(nchistop(nn)+1)-1
            call MPI_RECV(ortr,nrcv, 
     >         MY_MPI_REAL,nn-1,ntagort,MPI_COMM_WORLD,my_status,ierr)
            call MPI_RECV(flr,nrcv, 
     >         MY_MPI_REAL,nn-1,ntagfl,MPI_COMM_WORLD,my_status,ierr)
c$$$            print*,'Received on nodeid from n:',nodeid,nn
            do nsp = 1, nspm
               ortchil(npk(nchistart(nn)):npk(nchistop(nn)+1)-1,nsp) =
     >            ortr(1:nt,nsp)
               flchil(npk(nchistart(nn)):npk(nchistop(nn)+1)-1,nsp) =
     >            flr(1:nt,nsp)
            enddo
            deallocate(ortr,flr)
         enddo
         call mpi_barrier(MPI_COMM_WORLD,ierr)
         deallocate(orts,fls) 
!MPI section ends
         valuesin = valuesout
         if(theta .ne. 0.0) then
            if(inc .eq. 0) then
               ortchil_po = ortchil
            else
               call ortchilnsp_po(nr,KJ,npk,nchm,fl_po,
     >            maxf_po,minf_po,lo_po,nspm_po,ortchil_po) 
               call date_and_time(date,time,zone,valuesout)
               print '(i4,": nodeid, ortchilnsp_po call complete at: ",
     >  a10,", diff (secs):",i5)',nodeid,time,idiff(valuesin,valuesout)
            endif
         endif
         call clock(s2)
         call update(6)
      endif 

      td  = 0.0
      te1 = 0.0
      te2 = 0.0
      te3 = 0.0
      s1 = 0.0
      s2 = 0.0
      s3 = 0.0
      s4 = 0.0
      
c$$$C     The following is for the IBM
      if(.not.allocated(nchansi).and..not.allocated(nchansf)) then
         allocate(nchansi((nchif-nchii+1)*(nchm-nchii+1)))
         allocate(nchansf((nchif-nchii+1)*(nchm-nchii+1)))

         nch = 0
         do nchi = nchii, nchif
            do nchf = nchi, nchm
               nch = nch + 1
               nchansi(nch) = nchi
               nchansf(nch) = nchf
            enddo
         enddo
         nchansmax = nch
      endif
      IF(npsch.gt.0.and.nchatom.gt.0) THEN  !  if there is any Ps-state
         PRINT*,'Calculating Qlp(nchi,:,:) functions'
         if(.not.allocated(Qlpi)) then
            allocate(Qlpi(nchatom,kmax,2*kmax*igp))
            allocate(Q0pi(nchatom,kmax,2*kmax*igp))
            allocate(Alpi(nchatom,kmax,2*kmax*igp))
            allocate(aqi(nchatom,0:kmax))
            allocate(xpi(nchatom,2*kmax*igp+kmax))
            allocate(wpi(nchatom,2*kmax*igp))
            PRINT*, 'Qlpi, ... are allocated, :',nchii,nchatom
            Qlpi=0.; Q0pi=0.; aqi=0.; ionshi=0; nqgeni=0; nqimi=0
            xpi=0.; wpi=0.; imaxi=0; Alpi=0.
         endif
C$OMP PARALLEL DO ! why BornICS 0???
C$OMP& SCHEDULE(dynamic)
C$OMP& default(private)
C$OMP& shared(KJ,gridk,npk)
C$OMP& shared(Amg,gp)
C$OMP& shared(nchii,nchif,nchat,nchatom)
C$OMP& shared(Qlpi,Q0pi,aqi,ionshi,nqgeni,nqimi,xpi,wpi,imaxi,Alpi)
         DO nchi = nchii,nchif
            nqmi = npk(nchi+1) - npk(nchi)
            lioni=0
            call getchinfo(nchi, Ni, KJ, psii, maxpsii, ei, lia, nia,Li)
            posi = positron(nia,lia,nposi)
            if(nqmi.gt.1.and..not.posi)  then
               call qlandal(gp,Amg,nqmi,gridk,nchi,Li,Qlp,Q0p,aq,ionsh,
     *            nqgen,nqim ,xp,wp,imax,Alp)
               nchi1=nchat(nchi)
               Qlpi(nchi1,:,:)=Qlp(:,:)
               Q0pi(nchi1,:,:)=Q0p(:,:)
               aqi(nchi1,:)=aq(:)
               ionshi(nchi1)=ionsh
               nqgeni(nchi1)=nqgen
               nqimi(nchi1)=nqim
               xpi(nchi1,:)=xp(:)
               wpi(nchi1,:)=wp(:)
               imaxi(nchi1)=imax
               Alpi(nchi1,:,:)=Alp(:,:)        
            endif     
         ENDDO 
C$OMP END PARALLEL DO
         PRINT*,'Calculations of Qlp are done!'
      ENDIF
      print '("nodeid:",i4," nchii,nchif:",2i6)',
     >   nodeid,nchii,nchif

#ifdef GPU
!$omp parallel do num_threads(1)
#else
C$OMP PARALLEL DO
#endif
C$OMP& SCHEDULE(dynamic)
C$OMP& default(private)
C$OMP& shared(nchm,Nmax,la,sa,lpar,nspm,lo,itail,phasel)
C$OMP& shared(C,fl,maxf,minf,npk,chil,minchil,ortint,ortchil,flchil)
C$OMP& shared(KJ,gridk,vmat,vdon,ifirst,td,te1,te2,te3,na,nam,namax)
C$OMP& shared(is_core_orb)
C$OMP& shared(nicm,ncore,nspm_po,ortchil_po,lo_po)
C$OMP& shared(nchnsp_max,DD,get_nch,get_nsp,Etot_ve,theta_ve)
C$OMP& shared(nchps,maxrps,lpsmax,fd0)  
C$OMP& shared(nsmax,etot,etotp)
C$OMP& shared(Amg,gp)
C$OMP& shared(Bnlp,Dnlp,Unlp)      
C$OMP& shared(nchansi,nchansf,nchansmax,nchat)
C$OMP& shared(nodeid,scalapack,nze,nchistop)
C$OMP& shared(dwpot,i_sw_ng,nchii,nchif)
C$OMP& shared(vmat01,vmat0, vmat1,nr)
C$OMP& shared(ionshi,nqgeni,nqimi,imaxi,xpi,wpi,Qlpi,Q0pi,Alpi,aqi)
      do nch = 1, nchansmax
c$$$         print*,'OMP thread:',nch,omp_get_thread_num()
         nchi = nchansi(nch)
         nchf = nchansf(nch)
         nchi1=nchat(nchi)
         nchf1=nchat(nchf)

c$$$         if(MOD(nch*10,nchansmax).eq.0)
c$$$     >    print '(i4,"-",i3,"%",$)',nodeid,NINT(nch*100./nchansmax)
c$$$          call update(6)

         nqmi = npk(nchi+1) - npk(nchi)
         call getchinfo (nchi, Ni, KJ, psii, maxpsii, ei, lia, nia, Li)
         posi = positron(nia,lia,nposi)
         nqmf = npk(nchf+1) - npk(nchf)
c$$$           vmatt(:,:,:) = 0.0
         vmatt(1:nqmf,1:nqmi,0:nsmax) = 0.0
       call getchinfo (nchf,Nf,KJ,psif,maxpsif,ef,lfa,nfa,Lf)
            posf = positron(nfa,lfa,nposf)

            if (.not.posi.and..not.posf) then
               if(i_sw_ng.eq. 0) then  
                  call vdme(nze,nchm,Nmax,Ni,Li,la,sa,lpar,nspm,lo,
     >               C,fl,maxf,minf,npk,!chil,minchil,
     >               nchm,Nmax,Nf,Lf,la,sa,lpar,nspm,lo,
     >               C,fl,maxf,minf,npk,!chil,minchil,
     >               ortint,KJ,dwpot,
     >               vdon,vmatt,nchf,nchi,na,nam,namax,
     >               itail,gridk(1,nchi),phasel(1,nchi),gridk(1,nchf),
     >               phasel(1,nchf))
!       print*,'vdme is blocked for testing'
               elseif(i_sw_ng .eq. 1) then
                  call vdme_ng(nze,nchm,Nmax,Ni,Li,la,sa,lpar,nspm,lo,
     >                 C,fl,maxf,minf,npk,chil,minchil,
     >                 nchm,Nmax,Nf,Lf,la,sa,lpar,nspm,lo,
     >                 C,fl,maxf,minf,npk,chil,minchil,ortint,KJ,dwpot,
     >                 vdon,vmatt,nchf,nchi,na,nam,namax)
               else
                  print*,"Wrong value of i_sw_ng=",i_sw_ng
                  stop
               endif
            endif

c$$$            thetao = theta_ve 
c$$$            if (ei.gt.etot*2.0.or.ef.gt.etot*2.0) thetao = 0.0 !Igor
c!!        IF((posi.neqv.posf).and.ei.gt.0..and.ef.gt.0) CYCLE 
c!!        to test: no continuum-continuum rearrangment
 

            if (ifirst.ne.0) then
c     find maxumum value for maxf to be used in getformout2() and getformout3() routines.
               m_maxfi = 0
               do jni1=1,nam(Ni)   
                  ni1= na(jni1,Ni)
                  m_maxfi = max(m_maxfi,maxf(ni1))
               end do
               do jnf1=1,nam(Nf)   
                  nf1= na(jnf1,Nf)
                  m_maxfi = max(m_maxfi,maxf(nf1))
               end do

               nmr = m_maxfi !nr !maxr
               km = nqmi !kmax
               allocate(chiform(nmr,km))

c$$$  do k = 1, nqmi
c$$$                  do i = 1, maxr
c$$$                     chiform(i,k) = 0.0
c$$$                  enddo
c$$$                  minchiform(k) = maxr
c$$$                  maxchiform(k) = 0
c$$$               enddo 
               chiform(:,:) = 0.0
               minchiform(:) = maxr
               maxchiform(:) = 0
c$$$               call clock(s1)
               call ve2me(nchm,Nmax,Ni,Li,la,sa,lpar,nspm,lo,
     >            C,fl,maxf,minf,npk,chil,minchil,
     >            nchm,Nmax,Nf,Lf,la,sa,lpar,nspm,lo,
     >            C,fl,maxf,minf,npk,chil,minchil,
     >            ortint,ortchil,flchil,KJ,gridk,vmatt,nchf,nchi,
     >            chiform,minchiform,maxchiform,nmr,km,na,nam,namax)
c$$$               call clock(s2)
C  Store time for e2 matrix elements
c$$$               te2 = te2 + s2 - s1
               call ve1me(nchm,Nmax,Ni,Li,la,sa,lpar,nspm,lo,
     >            C,fl,maxf,minf,npk,chil,minchil,
     >            nchm,Nmax,Nf,Lf,la,sa,lpar,nspm,lo,
     >            C,fl,maxf,minf,npk,chil,minchil,
     >            ortint,KJ,vmatt,nchf,nchi,
     >            chiform,minchiform,maxchiform,nmr,km,na,nam,namax)
               deallocate(chiform)

c$$$  call clock(s3)
C  Store time for e1 matrix elements
c$$$               te1 = te1 + s3 - s2
c
c$$$               if(abs(theta) .ge. 1000.0) then
c$$$                  if(nchnsp_max .gt. 0) then
c$$$C vmatt still needs to be defined correctly in the routine below
c$$$                     call v_nu(nchm,nspm,npk,
c$$$     >                    nchi,Ni,Li,sa(Ni), nqmi,
c$$$     >                    nchf,Nf,Lf,sa(Nf),nqmf,
c$$$     >                    ortchil,Etot,KJ,na,nam,namax,lo,
c$$$     >                    nchnsp_max,DD,get_nch,get_nsp,vmatt)
c$$$                  endif
c$$$                  Etot_ve = 0d0
c$$$                  theta_ve = 0.0
c$$$               endif
c
               call ve2me12(nchm,Nmax,Ni,Li,la,sa,lpar,nspm,lo,
     >              C,fl,maxf,minf,npk,chil,minchil,
     >              nchm,Nmax,Nf,Lf,la,sa,lpar,nspm,lo,
     >              C,fl,maxf,minf,npk,chil,minchil,ortint,ortchil,
     >              flchil,Etot_ve,KJ,gridk,theta_ve,inc,vmatt,nchf,
c$$$     >              flchil,Etot_ve,KJ,gridk,thetao,inc,vmatt,nchf,
     >              nchi,na,nam,namax,nspm_po,ortchil_po,lo_po,nicm,
     >              ncore,is_core_orb)  
c$$$               call clock(s4)
C     Store time for both two electron exchange matrix elements
c$$$               te3 = te3 + s4 - s3
!               print*,'!!finish ve2me12'
            endif 
c$$$            call clock(s1)

CRRR adds by Rav: for Ps-formation and Ps-Ps transitions
         IF(posi.and.posf.and.nqmi*nqmf.gt.1) THEN
            npsi = nchps(nchi)  
            npsf = nchps(nchf) 
            fd01(:,:)=0.0
        maxps=min(maxpsii,maxpsif)
        if(allocated(fd0)) 
     >  fd01(1:maxps,0:lpsmax)=REAL(fd0(1:maxps,0:lpsmax,npsi,npsf)) 
        call pspsdme(nqmi,psii,maxpsii,lia,Li,nia,chil(1,npk(nchi),1),
     >     minchil(npk(nchi),1),gridk(1,nchi),nqmf,psif,
     >     maxpsif,lfa,Lf,nfa,chil(1,npk(nchf),1),
     >     minchil(npk(nchf),1),gridk(1,nchf),KJ,nchf,nchi,
     >     nchm,npk,posf,0,vdon,vmatt,
     >     maxrps,lpsmax,fd01,itail,phasel(1,nchi),
     >     phasel(1,nchf))
 
       Elseif (.not.posi.and.posf.and.nqmi.gt.1) then
       lm=min0(Lf+lfa+lia,Li+lfa+lia)
        Qlp(:,:)=Qlpi(nchi1,:,:)
        Q0p(:,:)=Q0pi(nchi1,:,:)
        aq(:)=aqi(nchi1,:)
        ionsh=ionshi(nchi1)
        nqgen=nqgeni(nchi1)
        nqim=nqimi(nchi1)
        xp(:)=xpi(nchi1,:)
        wp(:)=wpi(nchi1,:)
        imax=imaxi(nchi1)
        Alp(:,:)=Alpi(nchi1,:,:)   
      call posHeme(Nmax,nposf,gridk,npk,
     > nchi,nqmi,Li,Ni,la,sa,nchm,Bnlp(1,Ni),Dnlp(1,Ni),Unlp(1,Ni),
     > gp,nchf,nqmf,Lf,nfa,lfa,KJ,etotp,vmatt,lm,
     > Qlp,Q0p,Alp,aq,ionsh,nqgen,nqim ,xp,wp,imax)

!      PRINT*,'e+He ->He{+}+Ps ',nchi,nchf,vmatt(1,1,0),vmatt(10,10,0)
        

       Elseif (.not.posf.and.posi.and.nqmf.gt.1) then
       lm=min0(Li+lia+lfa,Li+lia+lfa)
        Qlp(:,:)=Qlpi(nchf1,:,:)
        Q0p(:,:)=Q0pi(nchf1,:,:)
        aq(:)=aqi(nchf1,:)
        ionsh=ionshi(nchf1)
        nqgen=nqgeni(nchf1)
        nqim=nqimi(nchf1)
        xp(:)=xpi(nchf1,:)
        wp(:)=wpi(nchf1,:)
        imax=imaxi(nchf1)
        Alp(:,:)=Alpi(nchf1,:,:)
      call posHeme(Nmax,nposi,gridk,npk,
     > nchf,nqmf,Lf,Nf,la,sa,nchm,Bnlp(1,Nf),Dnlp(1,Nf),Unlp(1,Nf),
     > gp,nchi,nqmi,Li,nia,lia,KJ,etotp,vmatt,lm,
     > Qlp,Q0p,Alp,aq,ionsh,nqgen,nqim,xp,wp,imax)
!      PRINT*,'He{+}+Ps ->e+He',nchf,nchi,vmatt(1,1,0),vmatt(10,10,0)

       ENDIF
CRRR end of adds by Rav

c$$$            call clock(s2)
C  Store time for direct matrix elements
c$$$            td = td + s2 - s1 
c$$$            vdon(nchf,nchi,0) = vmat(npk(nchf),npk(nchi)) +
c$$$     >         vdon(nchf,nchi,0)
c$$$            vdon(nchf,nchi,1) = vmat(npk(nchi),npk(nchf)+1) +
c$$$     >         vdon(nchf,nchi,1)
c$$$            vdon(nchi,nchf,0) = vdon(nchf,nchi,0)
c$$$            vdon(nchi,nchf,1) = vdon(nchf,nchi,1)
!            if (.not.posi.and.posf) then
!            endif
            if (npk(2)-npk(1).eq.1.or.
     >         (.not.scalapack.and.nodeid.eq.1)) then
!         write(65,*)
!         write(65,*)'#',nchi,'->',nchf

             do kf = 1, nqmf
                kff = npk(nchf) + kf - 1
!         write(65,*)
!         write(65,*) '# kf:',kff

               do ki = 1, nqmi
                  kii = npk(nchi) + ki - 1
!R                     if (kff.ge.kii) then
                        vmat(kff,kii) = vmat(kff,kii) + vmatt(kf,ki,0)
                        if (nsmax.eq.1) vmat(kii,kff+1) =
     >                     vmat(kii,kff+1)+vmatt(kf,ki,1)
!R                     endif

!                 write(65,'(4e20.10)') gridk(ki,nchi),
!     >          gridk(kf,nchf), vmat(kff,kii),vmatt(kf,ki,0)

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
     >                        + vmatt(kf,ki,0)
                           if (nsmax.eq.1) vmat01(kii,kff+1) = 
     >                        vmat01(kii,kff+1) + vmatt(kf,ki,1)
                        endif
                     enddo
                  enddo
               else
                  do ki = 1, nqmi
                     kii = npk(nchi) + ki - 1
                     do kf = 1, nqmf
                        kff = npk(nchf) + kf - 1
                        vmat0(kff,kii) = vmat0(kff,kii)
     >                     + vmatt(kf,ki,0)
                        if (nsmax.eq.1) vmat1(kii,kff+1) =
     >                     vmat1(kii,kff+1) + vmatt(kf,ki,1)
                     enddo
                  enddo
               endif
            endif 
c$$$         enddo
c$$$C$OMP END PARALLEL DO
         enddo 
C$OMP end parallel do
c     
      if(.not.scalapack.and.npk(2)-1.gt.1)then ! Do testing with LAPACK == 1 node
        if (nodeid.eq.1) then
        inquire(file='print_vmat',EXIST=ex)
        if(ex) then
           open(605,file='vmat_out')
           
           nchi = 1
           DO   nchf = 1,nchm
!           WRITE(605,*)
!           WRITE(605,*)'#', nchi,'->',nchf
!           nchf = 1
           kqi1 = npk(nchi) ! on shell
           kqi2 = npk(nchi) ! on shell
           kqf1 = npk(nchf)
           kqf2 = npk(nchf+1) - 1
           do kf=kqf1,kqf2
              kff=kf-kqf1+1
              do ki=kqi1,kqi2
                 kii=ki-kqi1+1
                 if(ki .gt. kf) cycle
!                 write(605,'(3e20.10)') gridk(kii,nchi),
!     >          gridk(kff,nchf), vmat(kf,ki)
              end do
           end do
           ENDDO 
          end if
        end if
      end if

#ifdef GPU
      call acc_set_device_num(gpunum,acc_device_nvidia)
!$acc exit data delete(chil(1:nr,npkstart:npkstop,1))
!$acc& delete(nchm,npk(1:nchm+1),minchil(npkstart:npkstop))
!$acc& finalize
#endif 
      if (allocated(fd0)) deallocate(fd0)
      if (allocated(flchil)) deallocate(flchil)
      if (allocated(ortchil)) deallocate(ortchil)
      if (allocated(ortchil_po)) deallocate(ortchil_po)
      if (allocated(nchansi)) deallocate(nchansi,nchansf)
        if(allocated(Qlpi)) then
         deallocate(Qlpi)
         deallocate(Q0pi) 
         deallocate(Alpi) 
         deallocate(aqi)  
         deallocate(xpi)  
         deallocate(wpi)  
        endif
      return
      end
c--------------------------------------------------------------------------
c 
      subroutine ortchilnsp(KJ,npk,nchm,fl,maxf,
     >   minf,lo,nspm,vdcore,dwpot,inc,ortchil,flchil)
      use chil_module
      include 'par.f'
      dimension  lo(nspmax), npk(nchm+1)
      dimension  fl(maxr,nspmax),maxf(nspmax),minf(nspmax),temp(maxr)
c$$$      dimension chil(nr,npk(nchm+1)-1), minchil(npk(nchm+1)-1)
      dimension ortchil(npk(nchm+1)-1,nspm),flchil(npk(nchm+1)-1,nspm)
      double precision Z
      common /Zatom/ Z
      common /meshrr/ nr, gridr(nmaxr,3)
      dimension dwpot(maxr,nchan),u(maxr),vdcore(maxr,0:lamax)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)      
      ze = 2.0

C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& private(nch,n,temp,maxt,ei,lia,nia,L,u,i,kqq,nsp,min1,max1)
C$OMP& private(tmp,tmp1)
C$OMP& shared(ze,ortchil,flchil,vdcore,dwpot,npk,nchii,nchif,nchm)
      do nch=nchii,nchif !nchm !return to nchm if removing MPI comms above
!         do nch=1,nchm
         call getchinfo (nch, N, KJ, temp, maxt, ei, lia, nia, L)
         if (dwpot(1,nch).eq.0.0) then
            do i = 1, nr
               u(i) = - 2.0 * vdcore(i,min(L,lamax)) ! * ze bug which stopped pw and dw results from being same 10/05/2016
            enddo
         else
            do i = 1, nr
               u(i) = dwpot(i,nch)
            enddo
         endif
         
         do kqq = npk(nch), npk(nch+1) - 1
            do nsp=1,nspm
               ortchil(kqq,nsp) = 0.0
               flchil(kqq,nsp) = 0.0
               if(lo(nsp).eq.L) then
                  temp(:) = 0.0
                  min1 = minf(nsp)
                  max1 = maxf(nsp)
                  do i = min1, max1
                     temp(i) = fl(i,nsp) * gridr(i,3)
                  enddo
                  call tmpfcexch(temp,min1,max1,nr,chil(1,kqq,1),tmp1,L)
c                  call fcexch(temp,nr,chil(1,kqq),tmp1,L)
                  tmp = 0.0
                  tmp1 = - tmp1 / ze
                  do i = max(minf(nsp),minchil(kqq,1)), maxf(nsp)
                     tmp = tmp + chil(i,kqq,1) * fl(i,nsp)
                     tmp1 = tmp1 + chil(i,kqq,1) * fl(i,nsp) *
     >                  (u(i)/2.0/ze + rpow2(i,0))   ! division by 2 is to go from Ry to au.
!     >                  (u(i)/2.0/ze + 1.0/gridr(i,1))   ! division by 2 is to go from Ry to au.
                  enddo 
                  ortchil(kqq,nsp)=tmp
                  flchil(kqq,nsp) = - ze * tmp1
                  if(Z.eq.-80.0) then
                     call relcorH12(lo(nsp),fl(1:nr,nsp),minf(nsp),
     >                  maxf(nsp),chil(1:nr,kqq,1),minchil(kqq,1),nr,
     >                  flchil(kqq,nsp))
                  end if

               end if
            end do

         end do
      end do
C$OMP END PARALLEL DO
      return
      end

c     Scattering from 2-electron atoms (2 el. above inert closed core).
c     The number of valence electrons: ze=2, 
c     Here we consider one electron in the field of a closed inert core.
c     The total core potential is :
c     V = V_{exchange} + Z/r + 2*V_{direct-core}.
c     The factor "2"  is to account two spin directions for closed inert core electrons.
c     The number of different electrons in the core 
c     is equal to 2*n_c.
c     The assymptotic charge of the target atom is 
c     Z_ass = Z - ze - 2*n_c.    (= 0 for neutral target).
c     The incident electron move in the potential of : (V(r) = V_{direct-core})
c     1. exchange potential
c     2. direct core potential = -Z/r + 2*V(r) = -(Z_ass +ze + 2*n_c)/r +2*V(r) =
c                              = -Z_ass/r - ze/r + 2(V(r) - n_c/r)
c     In Igor's programm: vdcore(r) = 2*(V(r) - n_c/r),
c                         exchange potential =  fcexch(temp,nr,fl(1,n1),tmp1,L) routine.
c     Potential -Z_ass/r is placed in Green function and therefore is removed from 
c     the potential.
c     This is valid for ions (Z_ass .not. 0) because Z_ass/r potential is accounted
c     by Green function.
    
c--------------------------------------------------------------------------
c
      subroutine ortchilnsp_po(nr,KJ,npk,nchm,fl,maxf,
     >   minf,lo,nspm,ortchil)
      use chil_module
      include 'par.f'
      dimension  lo(nspm),npk(nchm+1), LP(nchan)
      dimension  fl(nr,nspm),maxf(nspm),minf(nspm),
     >   temp(maxr)
c$$$      dimension chil(nr,npk(nchm+1)-1), minchil(npk(nchm+1)-1)
      dimension ortchil(npk(nchm+1)-1,nspm)
      
      do nch=nchii,nchm
!         do nch=1,nchm
         call getchinfo (nch, N, KJ, temp, maxt, ei, lia, nia, L)
         LP(nch) = L
      enddo

C$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC) COLLAPSE(1)
C$OMP&PRIVATE(nch,kqq,min1,max1,nsp)
      do nch=nchii,nchm
!         do nch=1,nchm
c$$$         call getchinfo (nch, N, KJ, temp, maxt, ei, lia, nia, L)   
         do kqq = npk(nch), npk(nch+1) - 1
            do nsp=1,nspm
               ortchil(kqq,nsp) = 0.0
c$$$               if(lo(nsp).eq.L) then
               if(lo(nsp).eq.LP(nch)) then
                  min1 = minf(nsp)
                  max1 = maxf(nsp)
                  ortchil(kqq,nsp) = dot_product(chil(min1:max1,kqq,1),
     >               fl(min1:max1,nsp))
c$$$     >                 SUM(chil(min1:max1,kqq,1)*fl(min1:max1,nsp))
!                  if(kqq .eq. 2) then
!                 print*,'po:', kqq,nsp,max1,
!     >                    fl(max1-2,nsp),fl(max1-1,nsp),
!     >                    fl(max1,nsp),ortchil(kqq,nsp)
!                  endif
               end if
            end do
         end do
      end do
C$OMP END PARALLEL DO
      return
      end

 
      subroutine makehepsff(nze,Nmaxi,Ni,lai,si,lpari,nspmi,loi,
     >   Ci,ortint,fli,maxfliall,nspm,maxfi,minfi,na,nam,namax,gp,Bnlp,
     >   Dnlp,Unlp,Amg,mpg,gridp)
   
      include 'par.f'
      include 'par.pos'

C!!! initial channel (target + e):
      implicit real*8 (a-h,o-z)
      double precision Z
      common /Zatom/ Z
      integer lai(Nmaxi), si(Nmaxi), lpari(Nmaxi)
      double precision Ci(Nmaxi,namax,namax),ortint(nspmax,nspmax)
      dimension loi(nspmi)
      real*8 fli(maxfliall,nspm) !fli(nmaxr,nspmax)
      dimension maxfi(nspmi),minfi(nspmi)
      real gridr
      common /meshrr/ nr, gridr(nmaxr,3)
      real tmp_const(0:ltmax)
      dimension na(nspmCI,KNM), nam(KNM)

C!!!  final channel (target-ion + Ps):
      parameter (maxl=2*ltmax+1)
      common/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      real veold, vaold
      common/fivec/ veold(maxr), vaold(maxr),vaoldne(maxr) 

      real*8 Fn2(nmaxr)
      real formf(nmaxr,nspmi),rpow1(nmaxr),rpow2(nmaxr),fun(nmaxr)
      real*8 Bnl(nmaxr),Dnl(nmaxr),dformf(nmaxr),Amg(iglp),gp(iglp),
     > Bnlp(iglp,Nmaxi),Dnlp(iglp,Nmaxi),Unlp(iglp,Nmaxi) 
      real cof3j,rli1,rli2,rlli,cl12,hat
      real gridp(nmaxr,3)
       real*8 dgridp(nmaxr,3),dgridr(nmaxr,3),formfn(nmaxr),Bnln(nmaxr)
     >  ,Dnln(nmaxr), veoldn(nmaxr)

!!*************  here pseudo-form-factors of He-like targets **********
! 
! Eqs. (2.134-136) in my thesis, Bnl, Dnl, uanl -on fine q^2 - mesh to be interpolated..

              dgridr(:,:)=1.d0*dble(gridr(:,:))
              dgridp(:,:)=1.d0*dble(gridp(:,:))


        do ip=1,iglp
         gp(ip)=1.e-8*exp(0.0065*ip) !0.00001*exp(0.00375*ip) !0.0001d0*exp(0.00625*ip)
        enddo 


********************************************************************
          Bnl(:)=0.d0
          Dnl(:)=0.d0
          formf(:,:)=0.
          maxBnl=1
          Fn2(:)=0.d0
            lli=lai(Ni)
            rlli=lli
            maxpsi=maxval(maxfi(:))
            minpsi = 1
c     these are 2 s.p. loops on s.p. functions of the coordinate r1
               do 10 jni1=1,nam(Ni)
                  ni1 = na(jni1,Ni)
                  li1=loi(ni1)
                  rli1=li1

!               if(maxrNi.lt.maxfi(ni1)) maxrNi=maxfi(ni1)

                     do 20 jni2=1,nam(Ni)
                           ni2 = na(jni2,Ni)
                           if(Ci(Ni,jni1,jni2).eq.0.0d0) goto 20
                         li2 = loi(ni2)
                         rli2 = li2
                         cl12 = hat(li1)*hat(li2)* (-1)**lli*
     ^                          cof3j(rli1,rli2,rlli,0.,0.,0.)
 
         const=Ci(Ni,jni1,jni2)*cl12*sqrt(2.0) ! 2./sqrt(2)

!        PRINT*, 'Ci::',Ni,ni1,ni2,const,ortint(ni2,1)

           if (const.eq.0.0) goto 20

!       write(17,'(a1,3i5,1e20.10)')'#', ni1,ni2,nam(Ni)
!     > ,Ci(Ni,jni1,jni2)
!      IF(ni1.ne.ni2) CYCLE 
!!! Integral for F_b(r_1) given in Eq.(2.104)
       fun(:)=0.
       maxpsi=min(maxfi(ni2),maxfi(ni1))
       maxBnl=max(maxBnl,maxpsi)
       do i=1,maxpsi
           r=dgridr(i,1)
!          rpow1(i)=r**li2        
!          rpow2(i)=1./r**(li2+1)
          hep=fli(i,1) !!   For He: 4.*sqrt(2.)*exp(-2.*r)
          fun(i)=hep*fli(i,ni2)*dgridr(i,3)/(2*li2+1)
       enddo
          rpow1(:)=gridr(:,1)**li2        
          rpow2(1:nr)=1./gridr(1:nr,1)**(li2+1)

      call form(fun,1,nr,rpow1,rpow2,1,nr,nr,
     >  formf(1,ni2),it1,it2)

!       PRINT*,'i1&i2', it1, it2,maxpsi,ni1,ni2
!       write(18,*) 
!       do i=1,nr,20
!          rr=gridr(i,1)
!        write(18,'(2i6,4e20.10)') ni1,ni2,rr, formf(i,ni2),
!     >      1./rr-(2.0+1./rr)*exp(-4.*rr),1./rr
!       enddo  



!      Fn2(:)=formf(:,ni2)*(1.d0-ortint(ni2,1)) !check ortint (must be 0 or 1)
!!!! below is the sum over configurations, as given in Eqs(2.102-103)
!      PRINT*, 'Sum over Ci=',Ni,ni1,ni2,ortint(ni2,1),cl12
      ort_wion = ortint(ni2,1)  !SUM(fli(1:maxfi(1),1)*fli(1:maxfi(1),ni2)
                                !     >              *gridr(1:maxfi(1),3))
      if(li2.ne.0) ort_wion=0.d0
      do ir=1,nr !min(minfi(ni1),minfi(ni2)),max(maxfi(ni1),maxfi(ni2))
       rr=dgridr(ir,1) 
      Fn2(ir)=(-2.0)/rr*ort_wion+formf(ir,ni2)
c$$$      Bnl(ir)=Bnl(ir)+const*ort_wion*fli(ir,ni1)
c$$$      Dnl(ir)=Dnl(ir)+const*Fn2(ir)*fli(ir,ni1) 
      enddo 
      do ir=1,maxfi(ni1)
         Bnl(ir)=Bnl(ir)+const*ort_wion*fli(ir,ni1)
         Dnl(ir)=Dnl(ir)+const*Fn2(ir)*fli(ir,ni1) 
      enddo 
    
 20                    continue
 10                  continue
!!   Make form-factors in momentum space: \int j_l(q*r)*Bnl(r) dr
!!   Results are in fine q-mesh wich to be interpolated later
!      do i=1,nam(Ni)
!      write(16,'(5i8)') Ni,lli,nam(Ni),maxfi(i),
!     >                  maxval(maxfi(1:nam(Ni)))
!      enddo
       IF(nam(Ni).gt.0) THEN 
         maxpsi = maxBnl !maxval(maxfi(1:nam(Ni)))
         minpsi=1  ! minval(minfi(1:nam(Ni)))
        ELSE
        STOP'stop!: nam(Ni) = 0, check F5 file lmax, l1max and etc.'
         maxpsi=nr
         minpsi=1
       ENDIF
!        print*, 'minpsi,maxpsi:',minpsi,maxpsi, Ni,nam(Ni),nspmi
! interpolation to a new rgrid 
        mpgmax = mpg
        mpgmin = 1
       do i=1,mpg
         if(dgridp(i,1).lt.dgridr(1,1)) mpgmin = mpgmin+1
         if(dgridp(i,1).gt.dgridr(maxpsi,1)) cycle
         mpgmax=i
       enddo
       mpgmin = max0(mpgmin,1)
       mpgmax = min0(mpgmax,mpg)

         do i=1,maxpsi
          dformf(i) = dble(formf(i,1))
         enddo
       call intrpl(maxpsi,dgridr(1,1),Bnl,mpgmax,dgridp(1,1),Bnln)
       call intrpl(maxpsi,dgridr(1,1),Dnl,mpgmax,dgridp(1,1),Dnln)
      call intrpl(maxpsi,dgridr(1,1),dformf,mpgmax,dgridp(1,1),formfn)
!       call intrpl(maxpsi,gridr(1,1),veold,mpgmax,gridp(1,1),veoldn)

   
      DO ip=1,iglp ! dense p2-grid for interpolation
  
         pp2=gp(ip)
         pp=sqrt(pp2)
         ppll=pp**lli

        res1=0.
        res2=0.
        res3=0.

        do i=mpgmin,mpgmax 
         rr=dgridp(i,1)
         arg = rr * pp
!         besskr=bes(lli,arg)
         call sbessel_qp(arg,lli,besskr) 
         rw=rr*dgridp(i,3)

!!!        f1s = 2.0/rr-formf(i,1)-vaold(i)/rr ! for He^{+1} -like ions  

        res1 = res1 + besskr*Bnln(i)*rw
        res2=res2 + besskr*Dnln(i)*rw
!!!        res3=res3 + besskr*f1s*Bnl(i)*rw
       enddo

        Bnlp(ip,Ni) = res1/ppll
        Dnlp(ip,Ni) = res2/ppll

       ENDDO  
!
      IF(Ni.eq.1) then
       Amg(:)=0.
       DO ip=1,iglp ! dense p2-grid for interpolation
 
         pp2=gp(ip)
         pp=sqrt(pp2)

        res4=0.0

        do i=mpgmin,mpgmax  
         rr=dgridp(i,1)
         arg = rr * pp

        f1sp=1.d0-formfn(i)*rr !+veoldn(i) ! f1sp=1/r-f1sp'; f1sp'=1/r-exp() for He{+1}
  
        res4=res4+ f1sp*sin(pp*rr)/pp*dgridp(i,3)
       enddo

        Amg(ip)=res4
        ENDDO        
       ENDIF

!!!!!!!
!        write(17,*) 
!        write(17,'(a1,5i5)'),'#',Ni,lli,maxBnl,maxpsi,minpsi  
!        do i=1, maxpsi 
!          rr=gridr(i,1) 
!        write(17,'(i5,5e20.10)') i,rr,Bnl(i),Dnl(i),formf(i,1),fli(i,1)    
!        enddo
! 
!        write(18,*)
!        write(18,*)'#', Ni,lli,mpgmin,mpgmax,gridp(mpgmax,1)
!        do i=mpgmin,mpgmax
!         write(18,'(i6,4e20.10)') i,gridp(i,1),Bnln(i),Dnln(i),formfn(i)
!        enddo
! 
!         write(19,*) 
!         write(19,'(a1,2i6)')'#',Ni,lli
!       DO ip=1,iglp ! dense p2-grid for interpolation
!          pp2=gp(ip)
!          pp=sqrt(pp2)
!          ppll=pp**lli
!        write(19,'(4e20.10)') pp,Bnlp(ip,Ni), Dnlp(ip,Ni),Amg(ip)   
!       ENDDO            
!!!!!!!!

  
       RETURN
       END


      subroutine makevstat(lo,fli,maxfliall,nspm,maxfi,formf0,
     >   ve2stat,va2stat)

      include 'par.f'
      include 'par.pos'
      implicit real*8 (a-h,o-z)       
      double precision Z
      common /Zatom/ Z
      dimension lo(nspm)
      real*8 fli(maxfliall,nspm) !fli(maxr,nspmax)
      integer  maxfi(nspm),minfi(nspm)
      real gridr,veold,vaold
      common /meshrr/ nr, gridr(nmaxr,3)
      common /fivec/ veold(maxr), vaold(maxr),vaoldne(maxr) 
      real formf(nmaxr,nspm),rpow1(nmaxr),rpow2(nmaxr),fun(nmaxr)
      real*8 ve2stat(maxr),va2stat(maxr)

       fun(:)=0.
       formf(:,:)=0.
        vaoldne(:) = 0.
!       maxpsi=max(maxfi(ni2),maxfi(1))
!        do ni2=1,nspm
           ni2=1
           li2=0
          maxpsi=maxfi(ni2) !nr !maxfi(1)
       do i=1,maxpsi
           r=gridr(i,1)
          hep=fli(i,1) !!   For He: 4.*sqrt(2.)*exp(-2.*r)
          fun(i)=hep*fli(i,ni2)*gridr(i,3)/(2*li2+1)
       enddo
          rpow1(:)=1.        
          rpow2(1:nr)=1./gridr(1:nr,1)

       call form(fun,1,maxpsi,rpow1,rpow2,1,nr,nr,
     >  formf(1,ni2),it1,it2)

               formf0=formf(1,1) !*dble(gridr(1,1))

        dzn = 0.d0 !Z-2.d0

       
        a_vee1 = 1.d0; b_vee1 = 0.4143; beta_vee1 = - 2.499;
 
        b_ver2 = 0.4; beta_ver2 = - 2.49;

        b_vep = 3.375; beta_vep = - 1.6875;

        frac = 1.0 ! 0.5
        do i=1,nr
           rr = gridr(i,1)

           vee1 = 1d0/rr - (1.d0/rr+b_vee1*rr)*exp(beta_vee1*rr) ! used by Igarashi etal.

!        vee2 =  - (1./rr + b_vee2*rr)*exp(beta_vee2*rr) ! used by Hewitt etal.

           vep = -1d0/rr + (1.d0/rr+b_vep*rr)*exp(beta_vep*rr) 

        ve2stat(i)= vep * rr*(1d0-frac) - rr*formf(i,1)*frac
!         ve2stat(i)= - rr * formf(i,1)

        va2stat(i)= vee1 * rr*(1d0-frac) + rr*formf(i,1)*frac
!         va2stat(i)= rr * formf(i,1)

!        vetest =  (1.d0/rr - (2.d0 + 1.d0/rr)*exp(-4.d0*rr))
!        vptest = - (1.d0/rr - (2.d0 + 1.d0/rr)*exp(-4.d0*rr))
!         write(17,'(6e20.10)') rr,formf(i,1), vee1,vep,vetest,vptest

        enddo

!       IF(Z.eq.2.) THEN
!        pcp1=2.499; pcp2=0.4143; pcp3=2.499;! pos-He^+1 potential parameters
!        pce1=2.499; pce2=0.4143; pce3=2.499;! electron-He^+1 pot. parameters
!        zmin1=1.d0
!        do i=1,nr
!        r=gridr(i,1)
!        va2stat(i)=1.d0-zmin1*exp(-pce1*r)-pce2*r*exp(-pce3*r) 
!        enddo
!        ENDIF

         RETURN
         END
        

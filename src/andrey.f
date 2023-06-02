      subroutine sd(max)
      real *8 PI_DP
      real, dimension (max) :: x
      real *8, dimension (max) :: y, z      
      dx = 1./real(max);
      PI = ACOS(-1.0)
      PI_DP =  ACOS(-1.0d0)
      do i = 1,max
         x(i) = sin(PI*dx*i)
         z(i) = sin(PI_DP*i/dble(max))
      end do
      y(:) =dble(x(:))
      print*
      print*, 'sd: x: ', x
      print*, 'sd: y: ', y
      print*, 'sd: z: ', z
      print*
      end subroutine sd
      
      
      real function electroncoreV(r)
C=====================================================================      
C     this is the electron-core potential Va = [va(r)-1]/r in r-space
      implicit none      
      real r      
      real , external :: vaint        
      electroncoreV  = (vaint(r) - 1.d0)/r      
      return
      end function electroncoreV

      function positroncoreVr(r)
C=====================================================================      
C     this is the electron-core potential Va = [ve(r)+1]/r in r-space
      implicit none      
      real, intent(in) :: r
      real :: positroncoreVr
      real, external :: veint         
      positroncoreVr = (veint(r) + 1.d0)/r      
      return
      end function positroncoreVr
      
*------------------------------------------------------

      function positroncoreV(q)
*     this is the positron-core potential in momentum space 
*     but w/o the pure Coulombic
*     positron-alkali_nucleus part (4*pi*NZNUC/q^2)      
      use apar                  ! maxx, nmax, lmax1, maxq, znuc1, Pi, Pi4
      use ftps_module, only: nr
C      use ql_module, only: dqlg
      implicit none
      real*8 q, positroncoreV
      integer, parameter :: io = 4 ! interpolation order            
      integer, parameter :: dp=selected_real_kind(precision(1.0d0))
      integer io1,j,k
      real(kind=dp), dimension (nr) :: qmesh, veq
      common /veqc2/ qmesh, veq
      real *8, dimension(io) :: xa, ya
      real *8 x, dx, xmin, q1, result
            
c$$$      real, dimension(3) :: c ! used in five
c$$$      common /veqc/ c                  
c$$$      if (q.lt.qmesh(1)) then
c$$$         positroncoreV = c(1)+c(2)*q**2+c(3)*q**4                  
c$$$      elseif (q.gt.qmesh(maxq)) then
c$$$         positroncoreV =Pi4*znuc1/q/q
c$$$      else                      ! ?????         
c$$$         call interpolate(io, q, maxq, qmesh, veq, positroncoreV)      
c$$$      end if      
c$$$      positroncoreV=positroncoreV + Pi4/q/q
      
      q1 = dble(q)      
      if (q.gt.qmesh(nr)) then
         positroncoreV =Pi4*znuc1/q/q                  
      else
         io1=io-1
         xmin = log(qmesh(1))
         ! dx = dqlg ;   
         dx = log(qmesh(2)/qmesh(1))
         x =  log(q1)                  
         j=max(int((x-xmin)/dx)+1,1)
!     k = min(max(j-io1/2,1),maxx-io1)
         k = min(max(j-io1/2,1),nr-io1)
         xa=dble(qmesh(k:k+io1))
         ya=dble(veq(k:k+io1))
         call intrpl(io, xa, ya, 1, q1, result)
         positroncoreV = result
      end if      
      return
      end function positroncoreV
*------------------------------------------------------

C 
C      
c$$$      subroutine  getf1z (pa, pb,  na, la, nb, lb, f1z)
c$$$!======================================================
c$$$!     computes F1z from Eq. 42 where
c$$$!        pa and pb  - momenta 
c$$$!        na, la, nb and lb - quantum numbers 
c$$$!------------------------------------------------------
c$$$      use ftps_module, only: oldftps
c$$$C     implicit real*8 (a-h,o-z)
c$$$      
c$$$* (1) interpolation
c$$$      if (oldftps) then      
c$$$         call getftps(pa, na, la, sk1, res0sk0)
c$$$         call getftps(pb,nb,lb,sk2,res0sk0)
c$$$      else
c$$$         call getftps3(pa, na, la, sk1, res0sk0)
c$$$         call getftps3(pb,nb,lb,sk2,res0sk0)
c$$$      end if
c$$$      
c$$$*     !      
c$$$* (2) direct numerical integration                                         
c$$$      ! call gnlp(bohr1,la, pa2, na ,res0sk0,sk1)
c$$$      ! call gnlp(bohr2,lb,pb2,nbi,res0sk0,sk2)      
c$$$*     V
c$$$* (3) analytical formula
c$$$      ! if(Nl1.eq.0) then
c$$$      !    call geigen(bohr1,na,la, pa2,sk1)
c$$$      ! else
c$$$      !    brap1=brp1*pa2+1.d0; x1=(brap1-2.d0)/brap1;
c$$$      !    sk1=pc1(Nl1-1,na,la)
c$$$      !    do k=Nl1-2,0,-1
c$$$      !       sk1=sk1*x1+pc1(k,na,la)
c$$$      !    end do
c$$$      !       sk1=sk1*bsp1/brap1**(la+2)
c$$$      ! end if                                                               
c$$$      !-------------------------------------------------
c$$$      ! if(Nl2.eq.0) then
c$$$      !    call geigen(bohr2,nb,lb, pb2,sk2)
c$$$      ! else
c$$$      !    brap2=brp2*pb2+1.d0
c$$$      !    x2=(brap2-2.d0)/brap2
c$$$      !    sk2=pc1(Nl2-1,nbi,lb)
c$$$      !    do k=Nl2-2,0,-1
c$$$      !        sk2=sk2*x2+pc1(k,nbi,lb)
c$$$      !    end do
c$$$      !    sk2=sk2*bsp2/brap2**(lb+2)
c$$$*     V end if
c$$$*--------------------------------------------------------
c$$$      f1z=sk2*sk1 
c$$$*--------------------------------------------------------
c$$$      return
c$$$      end 
      
      subroutine funleg2a(n, q, qa, rmax, Fmax, result)
c     --------------------------------------------------------     
c     estimates the integral of the form
c
c     dQlp = Integrate[f(r)*j(n,qa*r)j(n,q*r),{r,rmax,0}]
c      
c     which is omitted in funleg2
c      
c     this is done with the use of the asimptotic form of
c     the sherical bessels at large values of arguments
c     --------------------------------------------------------
C      print*,'1) funleg2a '
      z1 = q * rmax; z2 = qa * rmax
            
      a1 = q+qa ; a2 = q-qa;  
      
*     Cos(z1 + z2):
      C1 = (-1)**n*(Fc(n,z1)*Fc(n,z2) - Fs(n,z1)*Fs(n,z2))/2.
*     Sin(z1 + z2): 
      S1 = (-1)**n*(Fc(n,z2)*Fs(n,z1) + Fc(n,z1)*Fs(n,z2))/2.
*     Cos(z1 - z2):    
      C2 = (Fc(n,z1)*Fc(n,z2) + Fs(n,z1)*Fs(n,z2))/2.
*     Sin[z1 - z2]:
      S2 = (Fc(n,z2)*Fs(n,z1) - Fc(n,z1)*Fs(n,z2))/2.

      OIC1 = - sin(a1*rmax) * Fmax * C1 / a1;
      OIS1 =   cos(a1*rmax) * Fmax * S1 / a1;
      OIC2 = - sin(a2*rmax) * Fmax * C2 / a2;
      OIS2 =   cos(a2*rmax) * Fmax * S2 / a2;
      
      result =  OIC1+OIS1+OIC2+OIS2
C      print*,'2) funleg2a '
      return
      end

      function Fc(n,z)
      z2=z*z
      select case (n)
      case (0)
         Fc = 0.0
      case (1)        
         Fc = 1.0
      case (2)         
         Fc = 3.0
      case (3)         
         Fc = (3.*(-5.0 + 2.0*z2))/z2
      case (4)        
         Fc = (5.0*(-21.0 + 2.0*z2))/z2
      case default        
         a2 = float(n*(1 + n))/2.0
         a4 = float((-2+n)*(-1+n)*(2+n)*(3+n))/24.0
         a6 = float((-4 + n)*(-3 + n)*(4 + n)*(5 + n))/80.
         Fc = (1.0-(1.0-a6/z2)*a4/z2)*a2                  
      end select
      Fc = Fc/z2
      return
      end 

      function Fs(n,z)
      z2 = z*z
      select case (n)
      case (0:1)
         Fs = 1.0      
      case (2)         
         Fs =(-3.0 + z2)/z2
      case (3)
         Fs =(-15.0 + z2)/z2
      case default
         a2 = float((-1 + n)*n*(1 + n)*(2 + n))/8.0
         a4 = float((-3 + n)*(-2 + n)*(3 + n)*(4 + n))/48.0
         Fs = 1.0 - (1.0 - a4/z2) * a2/z2
      end select
      Fs = Fs/z            
      return
      end
            
      
      subroutine my_clock(time)
C     WALL-CLOCK TIME
      time = 0.0
      return
      call system_clock(COUNT=istart, COUNT_RATE=irate, COUNT_MAX=imax)      
      time = dble(istart)/dble(irate)
C     print'(a6,f14.5)',' time:',t1
      return
      end

      subroutine get_QL(ll, qa, maxq, q, QL)
      ! returns QL(q,qa) as an array for given ll and qa
      use ql_module
      implicit none
      integer, intent(in) :: ll, maxq
      real*8, intent(in) :: qa      
      integer, parameter :: io = 4, io1 = io-1
      real*8, dimension(maxq) :: q, QL
      real*8, dimension(maxql1) :: r1, ytmp
      !real*8, dimension(maxq) :: ytmp2
      real*8, dimension(io) :: x, ya2
      real*8 xlg, rlg
      
      integer ixlg,k,j
      
      xlg = log10(qa)
      ixlg = int((xlg-qlgmin)/dqlg)+1
      k = min(max(ixlg - io1/2,1), maxql2-io1) - 1
      do j = 1, io
         x(j) = qlgmin + dble(k+j-1)*dqlg
      end do                          
      do j = 1, maxql1
         rlg = qlgmin + dble(j-1)*dqlg
         r1(j) = 10.**rlg         
         ya2 = ql_res(k+1:k+io,j,ll)            
         call intrpl(io,x, ya2, 1, xlg, ytmp(j))  
      end do      
      call intrpl(maxql1,r1,ytmp,maxq,q,QL)
C      QL(:) = ytmp2(:)
c$$$      do j =1, maxq
c$$$         if (isnan(QL(j))) then
c$$$            print*, '(in) ll,qa,maxq: ', ll,qa,maxq
c$$$            print*,' qlgmin, dqlg:', qlgmin, dqlg
c$$$            print*,'j ', j
c$$$            print*,'QL', QL(j)
c$$$            STOP 'get_Ql: ERROR in QL'
c$$$         end if
c$$$      end do      
      return
      end subroutine get_QL
      
c$$$      subroutine get_QL(ll, qa, maxq, q, QL)
c$$$C     returns QL(q,qa) as an array for given ll and qa
c$$$      use ql_module      
c$$$      integer, parameter :: io = 4, io1 = io-1
c$$$      real, dimension(maxq) :: q, QL
c$$$      real, dimension(io) :: x, ya2
c$$$      real, dimension(io,io) :: ya
c$$$      real, dimension(maxql1) :: r1, ytmp      
c$$$            
c$$$      xlg = log10(qa)
c$$$      ixlg = int((xlg-qlgmin)/dqlg)+1
c$$$      k = min(max(ixlg - io1/2,1), maxql2-io1) - 1
c$$$      
c$$$      do j = 1, io
c$$$         x(j) = qlgmin + dble(k+j-1)*dqlg
c$$$      end do
c$$$                          
c$$$      do ir = 1, maxql1
c$$$c         rlg = qlgmin + dble(ir-1)*dqlg
c$$$c         r1(ir) = 10.**rlg         
c$$$         ya2 = ql_res(k+1:k+io,ir,ll)            
c$$$         call intrpl(io,x,ya2, 1, xlg, ytmp(ir))  
c$$$      end do
c$$$      
c$$$c      call intrpl(maxql1,r1,ytmp,maxq, q(1:maxq), QL(1:maxq))
c$$$      call intrpl(maxql1,qqq,ytmp,maxq, q(1:maxq), QL(1:maxq))
c$$$      
c$$$      return
c$$$      end subroutine get_QL

C=================================================
      subroutine get_ubb(minfun, maxfun, irmin, irmax, maxir, drlg, ll,
     $     rho,ubbb)
C     returns 1D array of ubb integrals for given ll and rho
      use ubb_module
      include 'par.f'

      real, intent(in)      :: drlg, rho
      integer, intent (in)  :: ll,minfun,maxfun,irmin,irmax
C      real*8,dimension(meshr),intent(out) :: ubbb
      real*8,dimension(minfun:maxfun), intent(out) :: ubbb
      
      
      common/meshrr/ meshr,rmesh(maxr,3)
      
      integer, parameter          :: io = 4, io1 = io-1      
      real*8, dimension(meshr)    :: tmp
      real*8, dimension(io)       :: xrho, ya2
C      real*8, dimension(ubb_max1) :: r1, ytmp
      real*8, dimension(irmin:irmax) :: r1, ytmp
      real*8 rholg, rlg
      
      rholg = dble(log10(rho))
C     irho = int((rholg-rlgmin)/drlg)+1
      irho = int((rholg-arholg(1))/drlg)+1
      krho = min(max(irho-io1/2,1),ubb_max2-io1)-1
      
c$$$      do j = 1, io         
c$$$         xrho(j) = dble(rlgmin) + dble(krho+j-1)*dble(drlg)
c$$$      end do
      
      xrho(1:io) = dble(arholg(krho+1:krho+io))

            
C     do ir = 1, ubb_max1
      do ir = irmin, irmax
c$$$         rlg = rlgmin + dble(ir-1)*drlg         
c$$$         r1(ir) = 10.0d0**rlg         
c$$$2         rlg = dble(arholg(ir))
c$$$2          r = dble(arho(ir))         
         ya2 = dble(ubb_res(ir,krho+1:krho+io,ll))            
         call intrpl(io,xrho, ya2, 1, rholg, ytmp(ir))         
      end do

      
      r1(irmin:irmax) = dble(arho(irmin:irmax))      
      
c$$$      tmp(1:meshr) = dble(rmesh(1:meshr,1))
c$$$      call intrpl(ubb_max1,r1,ytmp,meshr,tmp, ubbb)
c$$$C     ubbb(1:meshr) = tmp(1:meshr)

      maxtmp = maxfun - minfun + 1
      tmp(1:maxtmp) = dble(rmesh(minfun:maxfun,1))
C      call intrpl(ubb_max1,r1,ytmp,maxtmp, tmp, ubbb)      
      call intrpl(maxir,r1,ytmp,maxtmp,tmp,ubbb)

      
C     ubbb(1:meshr) = tmp(1:meshr)

C      stop ' get_ubb stopped'     
      return
      end subroutine get_ubb
      
c$$$      subroutine get_ubb_test(my_case, rlgmin, rlgmax, drlg, ll, rho,
c$$$     $     ubbb)
c$$$C     returns 1D array of ubb integrals for given ll and rho
c$$$C     different interpolation for comparison
c$$$      use ubb_module
c$$$      include 'par.f'
c$$$      common/meshrr/ meshr,rmesh(maxr,3)
c$$$      integer, parameter :: io = 4, io1 = io-1
c$$$      real, dimension(meshr) :: ubbb
c$$$      real, dimension(io) :: xrho, xr, ya2
c$$$      real, dimension(io,io) :: ya
c$$$      real, dimension(ubb_max1) :: r1, ytmp
c$$$      
c$$$C      open(unit=99, file = 'tmp.dat')
c$$$      
c$$$C     rlgmin = dlog10(rmesh(1,1))
c$$$C     rlgmax = dlog10(rmesh(meshr,1))
c$$$C     drlg = (rlgmax-rlgmin)/real(ubb_max1-1) 
c$$$C     m=io; n = io;
c$$$      
c$$$      rholg = log10(rho)
c$$$      irho = int((rholg-rlgmin)/drlg)+1
c$$$      krho = min(max(irho-io1/2,1),ubb_max2-io1)-1
c$$$      
c$$$      do j = 1, io
c$$$         xrho(j) = rlgmin + real(krho+j-1)*drlg
c$$$      end do
c$$$      
c$$$      select case (my_case)
c$$$      case(1)                   ! 2D polinomial interpolation
c$$$         do ir = 1, meshr               
c$$$            r = rmesh(ir,1)
c$$$            rlg = log10(r)               
c$$$            ir1 = int((rlg-rlgmin)/drlg)+1            
c$$$            kr1 = min(max(ir1-io1/2,1),ubb_max1-io1)-1
c$$$            do jrho = 1, io
c$$$               do j = 1, io     ! ky, ky+io1
c$$$                  xr(j) = rlgmin + real(kr1+j-1)*drlg
c$$$                  ya(jrho,j) = ubb_res(kr1+j, krho+jrho, ll)
c$$$               end do
c$$$            end do  
c$$$            call polin2(xrho,xr,ya,io,io,rholg,rlg,ubbb(ir),dy)
c$$$         end do                
c$$$      case (2)                  ! polint + Igor's routinue
c$$$         do ir = 1, ubb_max1
c$$$            rlg = rlgmin + real(ir-1)*drlg
c$$$            r1(ir) = 10.**rlg            
c$$$            ya2 = ubb_res(ir,krho+1:krho+io,ll)            
c$$$            call intrpl(io,xrho, ya2, 1, rholg, ytmp(ir))  
c$$$         end do
c$$$         call intrpl(ubb_max1,r1,ytmp,meshr,rmesh(1:meshr,1),
c$$$     $        ubbb(1:meshr))
c$$$      case (3)                  ! bilinear interplation
c$$$         i1 = irho
c$$$         i2 = irho+1
c$$$         x1 = rlgmin + real(i1-1)*drlg
c$$$         x2 = x1+drlg
c$$$         x = rholg
c$$$C         if ((x1.le.x).and.(x.le.x2)) then
c$$$         ir1 = 1
c$$$         y1 = rlgmin
c$$$         ir2 = ir1+1            
c$$$         f11 = ubb_res(ir1,i1,ll)
c$$$         f12 = ubb_res(ir2,i1,ll)
c$$$         f21 = ubb_res(ir1,i2,ll)
c$$$         f22 = ubb_res(ir2,i2,ll)         
c$$$         do i = 1, meshr               
c$$$            y = log10(rmesh(i,1))
c$$$ 10         y2 = y1+drlg               
c$$$            if (y.gt.y2) then 
c$$$               y1 = y2
c$$$               ir1 = ir2
c$$$               ir2 = ir1+1
c$$$               f11 = ubb_res(ir1,i1,ll)
c$$$               f12 = ubb_res(ir2,i1,ll)
c$$$               f21 = ubb_res(ir1,i2,ll)
c$$$               f22 = ubb_res(ir2,i2,ll)
c$$$               go to 10
c$$$            end if                                 
c$$$            ubbb(i) = (f11*(y2-y) + f12*(y-y1))*(x2-x)+ (f21*(y2-y) +
c$$$     $           f22*(y-y1))*(x-x1)
c$$$            ubbb(i) = ubbb(i)/(x2-x1)/(y2-y1)
c$$$         end do
c$$$C     else
c$$$C     stop 'get_ubb: irho is incorrect'
c$$$C     end if
c$$$      case default              ! exact calculation     
c$$$         do i = 1, meshr               
c$$$            r = rmesh(i,1)
c$$$            call Ubb0(ll, r, rho, res0, ubbb(i))
c$$$         end do
c$$$      end select
c$$$      return
c$$$      end subroutine get_ubb_test
            
      subroutine getubb(ir1, r1, meshr, rmesh, ll)
C     to check Ubb0
      use ubb_module
      real, dimension(meshr,3) :: rmesh
      
      logical exists
      character*100 command
      character*7 file1
      character*11 file2

      logical exchange
      
      file1 ='ubb_res'
      file2 ='ubb_res_old'

      exchange = .true.
      
      inquire(file=file1, exist=exists)      
      if (exists) then   
         write( command, '(3a,7a,1a,11a)') 'mv ', file1, ' ', file2
         print *, command
         call system(command)                
      end if
      open (88,file =file1)
      
      do ir2 = 1, ubb_max2, 10
         r2 = rmesh(ir2,1)/2.
         r = max(r1,r2)         
         call Ubb0(exchange,ll, r1, r2, res0, res1)
         write(88, '(4(e14.6,2x))') r2,res0/r,res1/r,ubb_res(ir2, ir1,
     $        ll)/r
C         write(88, '(2(e14.6,2x))') r2,ubb_res(ir2, ir1, ll) /r       
      end do
      close (88) 
      return
      end subroutine getubb
      
      
      subroutine makevexch(lactop, lamax, nabot, lnabmax, nnmax, ltmax,
     $     maxr, istoppsinb, psinb, rpow2, vdcore, nznuc, E, maxvexch,
     $     vexch)
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
c$$$      include 'par.f'
c$$$      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
c$$$     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
c$$$      common /psinbc/ enpsinb(nnmax,0:lnabmax),
c$$$     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
c$$$      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
c$$$     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
c$$$C     common/meshrr/ meshr,rmesh(maxr,3)

      dimension nabot(0:lamax), istoppsinb(nnmax,0:lnabmax)
      dimension psinb(maxr,nnmax,0:lnabmax), rpow2(maxr,0:ltmax)
      dimension vdcore(1:maxr)
      dimension density(1:maxr), vexch(1:maxr)
            
C     LOCAL VARIABLES:
      integer maxvexch
      real density, const, Rnl, ve, ve2, rho

      character fname*7  
      
C    for alkali targets only -------------------------------------------
      select case (nznuc)
C     case (3,11,19,37,55,87)
      case (3,11,19,37,55,87,12)  ! check if one active electron modeel works for magnesium 
         print '(3x,"calculate local exchange core potential")'
         print '(3x,"with energy parameter exchpar = ", f7.3," ... ",$)'
     $        , E 
      case default
         stop 'andrey.f: makevexch: not alkali target' 
      end select
      
C     calculate density ----------------------------------------------- 
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

c$$$C     Slater excahnge from Dima's Slater_exchange_nr
c$$$C     if(.true.) then
c$$$      pi = acos(-1d0)
c$$$      tmp = -(3d0/2d0)*(3d0/(4d0*pi*pi))**(1d0/3d0)  /2.0
c$$$c      tmp = -(3d0/2d0)*(24d0/pi)**(1d0/3d0) / 2d0
c$$$c      tmp = -(24d0/pi)**(1d0/3d0) / 2d0
c$$$      
c$$$      do i = 1, maxvexch  
c$$$          vexch(i)  = tmp*(density(i))**(1d0/3d0)
c$$$       end do 
c$$$       if (maxvexch.lt.maxr) vexch(maxvexch+1:maxr) = 0.0             
c$$$C=============================================================       
c$$$      open(unit=56,file='vexch')
c$$$      print*,'meshr', meshr
c$$$      do i=1,maxvexch+1
c$$$         r = 1./rpow2(i,0)         
c$$$         write(56,'(4(e12.5, 2x))') r, density(i), vexch(i) 
c$$$      end do
c$$$      close (56)
c$$$C      stop 'andrey: makevexch: see vexch' 
c$$$C==============================================================         
c$$$C       return
c$$$C       stop         
c$$$C      end if
         
C     make local-exchange potential ----------------------------------
C     - for details see Migdalek & Baylis PRA 24,649 (1981)
                  
      if (E.lt.0.0) then
         i0=1;
         do while (E + (rpow2(i0,0)-vdcore(i0)).gt.0.0)
            i0=i0+1
         end do
         i0=i0-2
      else 
         i0=maxvexch+1
      end if
                  
      do i = 1, maxvexch            
         ve = E + (rpow2(i,0)-vdcore(i)) ! rpow2(i,0) == 1/rmesh(i,1)
         ve2=ve*ve
         rho=density(i)                    
         if (i.lt.i0) then
            x= rho/ve2
            vexch(i) = -1.0*abs(0.5*rho/ve/(1.0+sqrt(1.0+x))) ! vexch(i) = 0.5*ve*x/(1.0+sqrt(1.0+rho/ve2))
         else
            !dv = vexch(i0-1)+0.5*sqrt(density(i0)) ! to fix  discontinuety at i0 ... it has to be sorted out
            vexch(i) = -0.5*sqrt(rho) ! + dv
         end if                     
      end do

c$$$      sgn = sign(1.0,E+(rpow2(1,0)-vdcore(1)))
c$$$      print *, ' ANDREY: sgn: ', sgn
c$$$      do i = 1, maxvexch
c$$$         ve = E + (rpow2(i,0)-vdcore(i)) ! rpow2(i,0) == 1/rmesh(i,1)
c$$$         rho=density(i)   
c$$$         If (sgn*ve.gt.0.0) then                    
c$$$            ve2=ve*ve                       
c$$$            x= rho/ve2
c$$$            vexch(i) = -1.0*abs(0.5*rho/ve/(1.0+sqrt(1.0+x))) ! vexch(i) = 0.5*ve*x/(1.0+sqrt(1.0+rho/ve2))
c$$$         else
c$$$            vexch(i) = -0.5*sqrt(rho)
c$$$         end if                     
c$$$      end do          
      
      if (maxvexch.lt.maxr) vexch(maxvexch+1:maxr) = 0.0      
      
c$$$C=============================================================       
c$$$      open(unit=56,file='vexch2')
c$$$      print*,'meshr', meshr
c$$$      do i=1,maxvexch
c$$$         r = 1./rpow2(i,0)
c$$$         ve = E + (rpow2(i,0)-vdcore(i))
c$$$         ve2 = ve*ve
c$$$         write(56,'(4(e12.5, 2x))') r, density(i), ve, vexch(i) 
c$$$      end do
c$$$      close (56)
c$$$      stop 'andrey: makevexch: see vexch' 
c$$$C==============================================================
      
      print '(" done")'
      
      return
      end subroutine makevexch    

      subroutine ubbgausleg(lmax, nmax, xgz, wgz, pleg)
C ================================================================
C     computes gauss. quadratures and legendre polinomials for Ubb
C     
C     calculates gausian quadratures (xgz and wgz) and legendre
c     polinomials (pleg) multiplied by wgz.  "lmax" is the maximum
c     orbital quantum number ll used for numerical calculation of Ubb
c     (see Eq.23) of pos-pos channels.
      
      implicit none       
      integer, intent(in) :: lmax, nmax
      real*8, dimension (1:2*(2**nmax-1)), intent(out) ::  xgz, wgz
      real*8, dimension (0:lmax, 1:2*(2**nmax-1)), intent(out) :: pleg
      
      integer i1, i2, i, ig, igz, n      
      real*8 z           
      real*8, dimension(0:lmax) :: pl8
      real*8, dimension (0:lmax,1:2**nmax) :: ptmp
                                ! DP due to cgqf in gausleg
      real *8, dimension (1:2**nmax) :: xtmp, wtmp      
      
C     gausian quadratures & legendre polinomials ---------------------
      i1 = 1
      do n=1, nmax                  ! for n = 64
         igz=2**n
         call gauleg(-1.d0,1.d0,xtmp(1:igz),wtmp(1:igz),igz)
         do ig=1,igz
            !z=real(xtmp(ig))
            z=xtmp(ig)
            pl8(0)=1.0; pl8(1)=z;
            !ptmp(0,ig) = real(wtmp(ig)); ptmp(1,ig) = z*real(wtmp(ig));
            ptmp(0,ig) = wtmp(ig); ptmp(1,ig) = z*wtmp(ig);
            do i = 2, lmax 
               pl8(i)=((2*i-1)*z*pl8(i-1)-(i-1)*pl8(i-2))/real(i)
               ptmp(i,ig)=pl8(i)*wtmp(ig)
            end do
         enddo         
C         call polleg2 (lmax, igz, xtmp(1:igz), wtmp(1:igz),
C     1        ptmp(0:lmax,1:igz))         
         i2=i1+igz-1
         xgz(i1:i2) = xtmp(1:igz)
         wgz(i1:i2) = wtmp(1:igz)
         pleg(:,i1:i2) = ptmp(:,1:igz)
         i1=i2+1
      enddo           
      return
      end subroutine ubbgausleg
      
      subroutine makepots(icase,true_potential,rmesh,maxr,vdcore,vexch,
     $     lamax,corep,r0,maxx,dx_sp,xmesh_sp,rmesh1_sp,va,vaM,ve,veold
     $     ,vaold)     
C
C     SP: va(r), ve(r) and vaM(r) are inaccurate at large value of r
C     
C=====================================================================
C     (1) calculates two arrays of va(r) and ve(r) at new rmesh. 
C         These functions  appear in the electron- and positron-core
C         potentials,
C     
C          Va(r) = [va(r) - 1]/r  and  Ve(r) = [1+ve(r)]/r
C
C         they are calculated (i) with the use of the HF method
C         (true_potential = true) or (ii) for some phenomenological
C         potentials specified with icase
C
C     (2) calculates the Fourier image of ve(r)/r needed for the re-
C         arrangement chanel matrix elements:
C
C                                 oo
C                          4 Pi   f
C                  veq  =  -----  | dr sin(q r) ve(r).
C                            q    J       
C                                 0
C     So,                    oo                           
C                            f     2                 4 Pi  
C              Ve(q) = 4 Pi  | dr r  j (q r) Ve(r) = ----  + veq.
C                            J        0              q*q             
C                            0        
C     INPUT
C     ----------------------------------------------------------------
C     true_potential - if true the HF-calculated  potentials are used
C     icase - select the pnenomenologialc potential
C     rmesh - old rmesh 
C     vdcore - a part of the static e-core interaction potenital 
C     corep and r0 - polarization-potential parameters
C     maxx - xmesh size
C     maxq - size of veq 
C      
C     OUTPUT
C     ----------------------------------------------------------------
C     xmesh, rmesh1 - new x and rmeshes of size maxx;
C     ve, va - see expression above;
C     veold  - table of ve but at the original rmesh
C     qmesh  - array of values q from 0 to 2Pi/rmash(maxr,1);
C     veq(q) - array of fourier image values at qmesh.
C
C=====================================================================
C      use ftps_module, only: nr, nexp 
      use apar, only: PI     ! maxx, nmax, lmax1, maxq, znuc1, Pi, Pi4 
      implicit none
      
C     input ---------------------------------------------------------
      integer maxr, lamax, maxx ! , maxq
      real, dimension  (maxr)   :: vdcore, vexch
      real, dimension  (maxr,3) :: rmesh, rm
      real, dimension (0:lamax) :: corep, r0
      real, external :: polpot
      integer icase
      logical true_potential
                  
C     output ---------------------------------------------------------
      real dx_sp
      real, dimension(maxx) :: va, ve, vaM, xmesh_sp, rmesh1_sp
      real, dimension(maxr) :: veold, vaold      
      
C     local varibales ------------------------------------------------
      real *8  dx
      real *8, dimension(maxx) :: xmesh, rmesh1, vtmp
      real *8, dimension (maxr) :: r, va1, ve1
!     real *8, dimension(maxq) :: veq, qmesh

C      real, dimension(nr) :: veq, qmesh
      
      real *8 pp, ri
      real *8 vaoldM(maxr) 
      integer la, i
      logical there
      real a,b,c

      real na, ca, ne, ce, naM, caM
      common /veva/ na, ca, ne, ce, naM, caM      
             
C     new meshes -----------------------------------------------------
      xmesh(1) = min(log10(rmesh(1,1)),1.e-10)      
      xmesh(maxx) = log10(rmesh(maxr,1))      
      dx = (xmesh(maxx) -  xmesh(1))/dble(maxx-1)      
      do i=1,maxx                
         xmesh(i) = (xmesh(1)+dble(i-1)*dx)  
         rmesh1(i) = 10.0d0**xmesh(i)  
      end do      
      
C     potentials on rmesh--------------------------------------------- 
      if (true_potential) then
         
!         print *, '*** full potentisls are used ***'           
         la=0                   ! ????
         do i=1, maxr           ! FULL POTENTIALS 
            r(i) = dble(rmesh(i,1))            
            pp = dble(polpot(rmesh(i,1),corep(la),r0(la)))
            va1(i) = (vdcore(i) + pp + vexch(i)) * r(i)
            ve1(i) = (-vdcore(i) + pp) * r(i)
c$$$  if (abs(ve1(i)).lt.1.e-16) print*, 'r, ve1: ', r(i), ve1(i)             
C     No exchange in e-ion potential (as Mitroy did): 
            vaoldM(i) = (vdcore(i) + pp ) * r(i)
         end do         
      else
!         print *, '*** phenomenological potentisls are used ***'
         do i=1, maxr          
            r(i) = dble(rmesh(i,1))             
            call phenpots(icase, r(i), va1(i), veold(i))             
         end do
      end if              

c$$$C====================================================================     
c$$$      open(unit=56,file='vave')
c$$$      print*, 'meshr, rmeshr(meshr,1): ', maxr, r(maxr)
c$$$      do i=1,maxr                                         
c$$$         write(56,'(4(e12.5, 2x))') r(i), va1(i), ve1(i)
c$$$      end do
c$$$      close (56)
c$$$      stop 'andrey: makepots: see vave: r(i), va1(i), ve1(i)'      
c$$$C====================================================================
      
c$$$C 10   continue
c$$$      la=0   
c$$$      do i=1, maxr              ! POTENTIALS WITHOUT HARTREE PART
c$$$         r(i)=dble(rmesh(i))
c$$$         pp=dble(polpot(rmesh(i),corep(la),r0(la)))
c$$$         vaold(i) = (vexch(i)+pp) * r(i)
c$$$         veold(i) = pp * r(i)    
c$$$      end do
            
C---- potentials on xmesh ---------------------------------------

      call intrpl(maxr, r, va1, maxx, rmesh1, vtmp)      
      va(:) = vtmp(:)
      
!     na and ca are needed to calculate va for r > rmax         
      na = log(va(maxx)/va(maxx-1))/log(rmesh1(maxx-1)/rmesh1(maxx))
      ca = va(maxx)*rmesh1(maxx)**na
            
      vaold(:) = va1(:)



      call intrpl(maxr, r, vaoldM, maxx, rmesh1, vtmp)
      vaM(:) = vtmp(:)      
!     naM and caM are needed to calculate vaM for r > rmax
      if (abs(vaM(maxx)).gt.1e-10) then
         naM=log(vaM(maxx)/vaM(maxx-1))/log(rmesh1(maxx-1)/rmesh1(maxx))
         caM=vaM(maxx)*rmesh1(maxx)**naM
      else
         caM = 0.0
         nam = 3
      endif
      
      call intrpl(maxr, r, ve1, maxx, rmesh1, vtmp)
      ve(:) = vtmp(:)
!     ne and ce are  needed to calculate ve for r > rmax
      if (abs(ve(maxx)).gt.1e-10) then
         ne = log(ve(maxx)/ve(maxx-1))/log(rmesh1(maxx-1)/rmesh1(maxx))
         ce = ve(maxx)*rmesh1(maxx)**na
      else
         ce = 0.0
         ne = 3
      endif 
      
      veold(:) = ve1(:)
      
      dx_sp = real(dx)
      rmesh1_sp(:) = rmesh1(:)
      xmesh_sp(:) = xmesh(:)

c$$$C=====================================================================C
c$$$      open(unit = 90, file = 'vexch')
c$$$      do i = 1, maxx 
c$$$C         write(90,'(4(E12.6,3x))') rmesh1_sp(i), va(i), ve(i), vaM(i)
c$$$         write(90,*) xmesh_sp(i), xmesh(i)
c$$$      end do
c$$$      close(90)
c$$$C=====================================================================C
      

C     va=0.                     ! FOR HYDROGEN 
C     ve=0.                     ! FOR HYDROGEN 
            
C     fourier image of ve(r)/r ------------------------------------
c$$$  call makeveq (maxr, veold, rmesh, maxq, qmesh, veq)
c$$$  call makeveq (maxr, rmesh, qmesh, veq)
      call makeveq (maxr, rmesh)
c$$$C=====================================================================C
c$$$      open(unit = 90, file = 'veq2')
c$$$      do i = 1, nr ! maxq 
c$$$         write(90,'(2(E12.6,3x))')  qmesh(i), veq(i)
c$$$      end do
c$$$      close(90)
c$$$C=====================================================================C      
c$$$      if (true_potential) then
c$$$         continue
c$$$      else
c$$$         if (icase.eq.0) then
c$$$            veq = 0.0
c$$$         end if
c$$$      end if      ! nexp
      return
      end subroutine makepots
      
      subroutine Ubb0(exchange,ll,r,rho,res0,res)
      !subroutine Ubb0(exchange,ll,corep,r0,r,rho,res0,res,resapp)
C------------------------------------------------------------      
C     Note that Ubb0(exchange, ll, r, rho, res0, res1) returns
C     
C     res0 = ((-1)**ll - 1) * X**ll
C
C     where X = min(r/2.,rho)/max(r/2.,rho). So, no scaling 
C     in makev31d is needed when Ubb0 is used.
C
C     no local exchange in Va(r) is set with the use of:
C           exchange = .false. 
C      
C------------------------------------------------------------ 
      implicit none       
      integer, intent(in) :: ll
      real, intent(in)  :: r,rho
      real, intent(out) :: res0, res
      real r1, r2
      logical, intent(in) :: exchange
      
      r1 = rho; r2 = 0.5 * r              
      if (r1.lt.r2) then       
         r1 = r2; r2 = rho
      end if

c$$$      !if (r2/r1.lt.0.1) then
c$$$         call UbbApp(ll,corep,r0,r1,r2,resapp)
c$$$         resapp = resapp * r1   ! it is devided by r1 in a different later
c$$$      !else
c$$$      !   resapp = 0.0
c$$$      !end if
         
      if (exchange) then       
         call Ubb3(ll, r1, r2, res0, res) ! with local exchange  
      else
         !if (r2/r1.lt.0.1.or.r1.lt.2) then
         !   res = resapp
         !else
            call Ubb3M(ll, r1, r2, res0, res) ! without local exchange
         !end if
      end if
            
C      res0=res0                 ! /r1 
C      res=res                   ! /r1
      return             
      end subroutine Ubb0


      subroutine UbbApp(ll,corep,r0,r1,r2,res)                  
C     calculates approximate values of Ubb's analytically
C     for small x = min(r1,r2)/max(r1,r2) and relatively small max(r1,r2)
C     ll must be even
C     corep and r0 are the polpot parameters for ll = 0
      !implicit none
      integer, intent(in) :: ll
      real, intent(in)  :: corep,r0,r1,r2
      real, intent(out) :: res
      
c$$$      r1 = rho; r2 = 0.5 * r              
c$$$      if (r1.lt.r2) then       
c$$$         r1 = r2; r2 = rho
c$$$      end if
      
      x = min(r1,r2)/max(r1,r2)
      rmax = max(r1,r2)

      rr0 = (rmax/r0)**6
      y = rr0
      expr = exp(-rr0)

      coef = 2d0/rmax**4

      !if (ll.gt.0) print'(1x," ll = ", i3.3)',ll
      select case (ll)
      case (0)
         if (rmax.lt.0.1) then
            res = (-1. - 1.*x**2)*y + (0.5 + 6.*x**2 + 12.6*x**4 + 6.*x
     $           **6 + 0.5*x**8)*y**2 +(-0.16666666666666666 - 5
     $           .833333333333333*x**2 - 45.5*x**4 - 119.16666666666667
     $           *x**6 - 119.16666666666667*x**8 - 45.5*x**10 -5
     $           .833333333333333*x**12 - 0.16666666666666666*x**14)*y
     $           **3
         else
            res =  -1. - 2.*x**2 - 3.*x**4 - 4.*x**6 - 5.*x**8 + (1. + x
     $           **2*(2. + y + 6.*y**2) + x**4*(3. + 3.*y + 14.1*y**2 -
     $           32.4*y**3 + 10.8*y**4) +x**6*(4. + 4.*y + 8.*y**2 - 112
     $           .5*y**3 +207.*y**4 - 87.94285714285714*y**5 + 9
     $           .257142857142858*y**6) +x**8*(5. + 5.*y + 3.*y**2 - 117
     $           .83333333333333*y**3+ 823.375*y**4 - 1173.*y**5 + 538.2
     $           *y**6 - 89.48571428571428*y**7 +4.628571428571429*y**8)
     $           )*expr
         end if
         res = coef*res
      case (2)
         if (rmax.lt.0.4) then
            res = (1.6*x**2 + 4.114285714285714*x**4 + 1.6*x**6)*y**2 +
     $           (-1.8666666666666667*x**2 - 20.*x**4 - 57
     $           .77777777777778*x**6 - 57.77777777777778*x**8 - 20.*x
     $           **10 - 1.8666666666666667*x**12)*y**3 + (x**2 + 24.*x
     $           **4 + 177.33333333333334*x**6 + 548.1212121212121*x**8
     $           + 790.5594405594405*x**10 + 548.1212121212121*x**12 +
     $           177.33333333333334*x**14 + 24.*x**16 + x**18)*y**4
         else
            res = -0.006926406926406926*x**2*(231. + 396.*x**2 + 550.*x
     $           **4 + 700.*x**6) +expr*(x**2*(1.6 + 1.6*y + 2.4*y**2) +
     $           x**4*(2.742857142857143 + 2.742857142857143*y + 5
     $           .485714285714286*y**2 -15.428571428571429*y**3 + 6
     $           .171428571428572*y**4) +x**6*(3.8095238095238093 + 3
     $           .8095238095238093*y + 3.5047619047619047*y**2 - 55
     $           .542857142857144*y**3 + 120.51428571428572*y**4 -55
     $           .542857142857144*y**5 + 6.171428571428572*y**6) +x**8
     $           *(4.848484848484849 + 4.848484848484849*y + 2
     $           .4242424242424243*y**2 - 56.96969696969697*y**3 + 490
     $           .54545454545456*y**4 -762.1558441558442*y**5 + 368
     $           .1350649350649*y**6 - 63.397402597402596*y**7 + 3
     $           .366233766233766*y**8))
         end if
         res = coef*res
      case (4)
         if (rmax.lt.0.8) then
            res =  0.20317460317460317*x**4*y**2 + (-2.3703703703703702
     $           *x**4 - 9.696969696969697*x**6 - 9.696969696969697*x**8
     $           -2.3703703703703702*x**10)*y**3 + (3.5555555555555554*x
     $           **4 + 40.72727272727273*x**6 + 148.8111888111888*x**8
     $           +224.87024087024088*x**10 + 148.8111888111888*x**12 +
     $           40.72727272727273*x**14 + 3.5555555555555554*x**16)*y
     $           **4
         else
            res = -0.014208014208014208*x**4*(143. + 234.*x**2 + 315.*x
     $           **4 + 392.*x**6) +expr*(x**4*(2.0317460317460316 + 2
     $           .0317460317460316*y + 1.2190476190476192*y**2 - 1
     $           .8285714285714285*y**3 + 1.3714285714285714*y**4) +x**6
     $           *(3.324675324675325 + 3.324675324675325*y + 1
     $           .6623376623376624*y**2 - 9.142857142857142*y**3 + 31
     $           .16883116883117*y**4 -17.57922077922078*y**5 + 2
     $           .244155844155844*y**6) +x**8*(4.475524475524476 + 4
     $           .475524475524476*y + 2.237762237762238*y**2 - 8
     $           .951048951048952*y**3 + 139.3006993006993*y**4 -267
     $           .42857142857144*y**5 + 146.34485514485513*y**6 - 27
     $           .44775224775225*y**7 + 1.5536463536463536*y**8) +x**10
     $           *(5.5695415695415695 + 5.5695415695415695*y + 2
     $           .7847707847707848*y**2 - 1.442113442113442*y**3 + 222
     $           .73193473193473*y**4 -1248.2685314685314*y**5 + 1765
     $           .4601398601399*y**6 - 921.0965034965035*y**7 + 203
     $           .7866133866134*y**8 - 19.16163836163836*y**9 +0
     $           .6214585414585415*y**10))
         end if
         res = coef*res
      case (6)
         if (rmax.lt.0.8) then
            res = (-0.3978243978243978*x**6 - 0.3978243978243978*x**8)*y
     $           **3 +(2.983682983682984*x**6 + 16.708624708624708*x**8
     $           + 28.011517893870835*x**10 + 16.708624708624708*x**12 +
     $           2.983682983682984*x**14)*y**4 + (-4.876190476190477*x
     $           **6 - 61.44*x**8 - 271.05882352941177*x**10 - 546
     $           .8730650154798*x**12 - 546.8730650154798*x**14 -271
     $           .05882352941177*x**16 - 61.44*x**18 - 4.876190476190477
     $           *x**20)*y**5
         else
            res = -0.0014779853789141715*x**6*(1615. + 2584.*x**2 + 3420
     $           .*x**4 + 4200.*x**6 + 4950.*x**8) +expr*(x**6*(2
     $           .386946386946387 + 2.386946386946387*y + 1
     $           .1934731934731935*y**2 + 2.6853146853146854*y**4 - 2
     $           .0715284715284716*y**5 +0.34525474525474525*y**6) + x
     $           **8*(3.819114219114219 + 3.819114219114219*y + 1
     $           .9095571095571096*y**2 + 0.2386946386946387*y**3 +16
     $           .46993006993007*y**4 - 44.89846153846154*y**5 + 30
     $           .29034965034965*y**6 - 6.55984015984016*y**7 + 0
     $           .4143056943056943*y**8) +x**10*(5.0547099958864665 + 5
     $           .0547099958864665*y + 2.5273549979432333*y**2 + 0
     $           .8424516659810778*y**3 + 28.222130810366103*y**4 -243
     $           .00518305224188*y**5 + 427.62846565199504*y**6 - 256
     $           .43085385203034*y**7 + 62.67592172533349*y**8 - 6
     $           .3608109537521305*y**9 +0.21933830875007346*y**10) + x
     $           **12*(6.20753859143952 + 6.20753859143952*y + 3
     $           .10376929571976*y**2 + 1.03458976523992*y**3 +16
     $           .96727214993469*y**4 - 530.1127108185932*y**5 + 2248
     $           .2352421572236*y**6 - 2956.616267942584*y**7 + 1616
     $           .4688215808958*y**8 -417.5460650804304*y**9 + 52
     $           .99713784667345*y**10 - 3.1515451730931607*y**11 + 0
     $           .06926472907897056*y**12) +x**14*(7.316027625625149 + 7
     $           .316027625625149*y + 3.6580138128125745*y**2 + 1
     $           .2193379376041915*y**3 + 3.2885174680840317*y**4 -543
     $           .8284151349167*y**5 + 5724.806027409123*y**6 - 15502
     $           .697218570902*y**7 + 16468.321503883115*y**8 - 8311
     $           .429984257229*y**9 +2195.471851987022*y**10 - 316
     $           .33394172391075*y**11 + 24.69864797740959*y**12 - 0
     $           .9697062071055879*y**13 +0.014842441945493692*y**14))
         end if
         res = coef*res   
      case (8)
          if (rmax.lt.1.0) then
             res = (0.5616344439873852*x**8 + 1.241507718287904*x**10 +
     $            0.5616344439873852*x**12)*y**4 +(-3.212549019607843*x
     $            **8 - 22.826006191950466*x**10 - 54.34763379035825*x
     $            **12 - 54.34763379035825*x**14 - 22.826006191950466*x
     $            **16 -3.212549019607843*x**18)*y**5 + (5
     $            .3542483660130715*x**8 + 74.39587203302374*x**10 + 384
     $            .3786721706226*x**12 +969.3027385172223*x**14 + 1308
     $            .55869699825*x**16 + 969.3027385172223*x**18 + 384
     $            .3786721706226*x**20 + 74.39587203302374*x**22 +5
     $            .3542483660130715*x**24)*y**6
         else
            res =  -0.0008812832073028601*x**8*(3059. + 4830.*x**2 +
     $           6325.*x**4 + 7700.*x**6) +expr*(x**8*(2
     $           .6958453311394486 + 2.6958453311394486*y + 1
     $           .3479226655697243*y**2 + 0.4493075551899081*y**3 + 0
     $           .6739613327848621*y**4 -2.6284491978609625*y**5 + 2
     $           .426260798025504*y**6 - 0.6498912851854028*y**7 + 0
     $           .04874184638890521*y**8) +x**10*(4.256597891272814 + 4
     $           .256597891272814*y + 2.128298945636407*y**2 + 0
     $           .7094329818788023*y**3 + 1.4188659637576047*y**4 -21
     $           .54902682456862*y**5 + 52.19653164173288*y**6 - 38
     $           .23413045159175*y**7 + 10.789905574301859*y**8 - 1
     $           .2236768803951465*y**9 +0.04617648605264704*y**10) + x
     $           **12*(5.5741162861905895 + 5.5741162861905895*y + 2
     $           .7870581430952948*y**2 +0.929019381031765*y**3 + 0
     $           .7938892892453264*y**4 - 53.739548377319274*y**5 + 330
     $           .31959743043336*y**6 - 536.2042878837218*y**7 +338
     $           .0320343080184*y**8 - 97.03108991848488*y**9 + 13
     $           .37743795346626*y**10 - 0.8509666715416384*y**11 + 0
     $           .01978992259399159*y**12)+ x**14*(6.785880696232022 + 6
     $           .785880696232022*y + 3.392940348116011*y**2 + 1
     $           .1309801160386703*y**3 + 0.2827450290096676*y**4 -54
     $           .29108478455632*y**5 + 914.9645295611643*y**6 - 3136
     $           .909069447177*y**7 + 3871.9662694016174*y**8 - 2174
     $           .7741866963297*y**9 +623.2955625341142*y**10 - 95
     $           .84749480160437*y**11 + 7.896035709766456*y**12 - 0
     $           .32438264425803603*y**13 +0.005162588502780415*y**14))
         end if
         res = coef*res  
      case (10)
          if (rmax.lt.0.4) then
             res =  0.01126084098220321*x**10*y**4 + (-0
     $            .6441201041820237*x**10 - 2.2684229755975616*x**12 - 2
     $            .2684229755975616*x**14 -0.6441201041820237*x**16)*y
     $            **5 + (3.0058938195161105*x**10 + 25.876825054964776*x
     $            **12 + 80.21815767039081*x**14 +114.88032456500413*x
     $            **16 + 80.21815767039081*x**18 + 25.876825054964776*x
     $            **20 + 3.0058938195161105*x**22)*y**6 +(-4
     $            .953590325018896*x**10 - 75.5960958296362*x**12 - 447
     $            .52888731144634*x**14 - 1353.6367579173377*x**16 -
     $            2310.517914376145*x**18 -2310.517914376145*x**20 -
     $            1353.6367579173377*x**22 - 447.52888731144634*x**24 -
     $            75.5960958296362*x**26 - 4.953590325018896*x**28)*y**7
         else
            res = -0.002872330453431544*x**10*(1035. + 1620.*x**2 + 2106
     $           .*x**4 + 2548.*x**6) +expr*(x**10*(2.972862019301648 +
     $           2.972862019301648*y + 1.486431009650824*y**2 + 0
     $           .4954770032169413*y**3 + 0.13513009178643853*y**4 -0
     $           .6080854130389735*y**5 + 2.3715331108519964*y**6 - 2
     $           .267289897188172*y**7 + 0.7769377018381883*y**8 -0
     $           .10261441345032676*y**9 + 0.004397760576442576*y**10)
     $           +x**12*(4.6531753345591005 + 4.6531753345591005*y + 2
     $           .3265876672795502*y**2 + 0.7755292224265168*y**3 + 0
     $           .1938823056066292*y**4 -2.2296465144762356*y**5 + 23
     $           .614864822887437*y**6 - 50.85255901339589*y**7 + 38
     $           .72833676332705*y**8 - 12.787922928368673*y**9 +1
     $           .9686670823935981*y**10 - 0.13709540579605767*y**11 + 0
     $           .0034417256685202765*y**12) +x**14*(6.049127934926831 +
     $           6.049127934926831*y + 3.0245639674634157*y**2 + 1
     $           .0081879891544718*y**3 + 0.25204699728861796*y**4 -2
     $           .218013576139838*y**5 + 77.95813626136953*y**6 - 368
     $           .44374090505767*y**7 + 554.8894760247348*y**8 - 358
     $           .73993843381396*y**9 +114.44013185883401*y**10 - 19
     $           .17018252528004*y**11 + 1.6949466399761806*y**12 - 0
     $           .07392826735981554*y**13 +0.0012390212406672995*y**14)
     $           + x**16*(7.318697995343573 + 7.318697995343573*y + 3
     $           .6593489976717866*y**2 +1.2197829992239289*y**3 + 0
     $           .3049457498059822*y**4 - 0.5831309542208273*y**5 + 114
     $           .24636931914897*y**6 -1239.0770412818065*y**7 + 3576
     $           .5573934837093*y**8 - 4164.874992596581*y**9 + 2371
     $           .0445406611157*y**10 -730.7314753099909*y**11 + 128
     $           .48496875230225*y**12 - 13.099290750746858*y**13 + 0
     $           .7573708540560453*y**14 -0.022807168763394364*y**15 + 0
     $           .0002753380534816221*y**16))
         end if
         res = coef*res  
      case (12) 
         if (rmax.lt.2) then
            res = (-0.02688501304411925*x**12 - 0.02688501304411925*x
     $           **14)*y**5 +(0.6273169710294492*x**12 + 3
     $           .066882969477307*x**14 + 4.91758820967913*x**16 + 3
     $           .066882969477307*x**18 + 0.6273169710294492*x**20)*y**6
     $           + (-2.481102632357291*x**12 - 25.08670439383483*x**14 -
     $           96.02152371433331*x**16 - 180.68566290331538*x**18 -180
     $           .68566290331538*x**20 - 96.02152371433331*x**22 - 25
     $           .08670439383483*x**24 - 2.481102632357291*x**26)*y**7
         else
            res = -0.0003987395334685836*x**12*(8091. + 12586.*x**2 +
     $           16275.*x**4 + 19600.*x**6) +expr*(x**12*(3
     $           .22620156529431 + 3.22620156529431*y + 1
     $           .613100782647155*y**2 + 0.537700260882385*y**3 + 0
     $           .13442506522059625*y**4 +0.6049127934926831*y**6 - 1
     $           .866588048491708*y**7 + 1.808257171976342*y**8 - 0
     $           .7134314896879365*y**9 +0.12653313213333212*y**10 - 0
     $           .009912169925338396*y**11 + 0.0002753380534816221*y**12
     $           ) +x**14*(5.018535768235593 + 5.018535768235593*y + 2
     $           .5092678841177967*y**2 + 0.8364226280392656*y**3 + 0
     $           .2091056570098164*y**4 +0.014936118357844027*y**5 + 3
     $           .0469681450001818*y**6 - 22.032268189655724*y**7 + 42
     $           .77741637982424*y**8 - 33.0282385247005*y**9 +12
     $           .058936538276638*y**10 - 2.250688032333799*y**11 + 0
     $           .21768022554143354*y**12 - 0.010248694212927045*y**13
     $           +0.00018355870232108143*y**14) + x**16*(6
     $           .489485907201198 + 6.489485907201198*y + 3
     $           .244742953600599*y**2 +1.0815809845335331*y**3 + 0
     $           .2703952461333833*y**4 + 0.054079049226676654*y**5 + 4
     $           .926601384550243*y**6 -91.10264790824402*y**7 + 349
     $           .0522382223307*y**8 - 489.9214788324546*y**9 + 319
     $           .7639425736048*y**10 - 109.61816254994211*y**11 +21
     $           .013785120981385*y**12 - 2.3028926664629505*y**13 + 0
     $           .14162978046417027*y**14 - 0.004500353012078927*y**15
     $           +0.00005696649382378389*y**16) + x**18*(7
     $           .815294855984239 + 7.815294855984239*y + 3
     $           .9076474279921194*y**2 +1.3025491426640399*y**3 + 0
     $           .32563728566600997*y**4 + 0.06512745713320199*y**5 + 3
     $           .077737545666174*y**6 -177.61722928009678*y**7 + 1381
     $           .8980439979323*y**8 - 3450.5389811048167*y**9 + 3804
     $           .0872947840758*y**10 -2174.451534067602*y**11 + 705
     $           .7017838810053*y**12 - 136.89614081204041*y**13 + 16
     $           .253315554750134*y**14 -1.1789728901547962*y**15 + 0
     $           .05049712151711837*y**16 - 0.0011632190513050065*y**17
     $           + 0.00001102577299815172*y**18))
         end if
         res = coef*res   
      case default
         !print'(1x," ll = ", i3.3)',ll
         stop 'UbbApp: no analytical expression for this ll'
      end select
      res = corep*res*(2.0*ll+1.0)/2.0
      
      return
      end subroutine UbbApp
      
      
      subroutine Ubb3(ll, r1, r2, res0, res)   
C======================================================================
C     is used to get 
C      
C     res = U (r1,r2) =                                      (Eq. 23)
C            ll                                 
C                (2ll+1)  f 1     _   _         _   _
C              = -------- |  (Ve(|r1+r2|) + Va(|r1-r2|)) P  (z) dz
C                   2     j-1                             ll
C      
C     for alkalies at r1 >= r2 . The integrand is divergent at r1 = r2.
C     P  (z) is Legendre function, res0 is what one has for hydrogen.
C      ll 
C======================================================================
      use apar, only: nmax, lmax1 ! maxx, nmax, lmax1, maxq 
      implicit none    
 
      integer, intent(in) :: ll
      real, intent(in)    :: r1, r2
      real, intent(out)   :: res0, res
      
      real*8, dimension (1:2*(2**nmax-1)) ::  xgaus, wgaus
      real*8, dimension (0:lmax1, 1:2*(2**nmax-1)) :: pleg 
      common /gausc/ xgaus, wgaus, pleg      

      real coef, coef1
      
!     integer i, i1, i2, ig, igz, j, jrange, k, kk,
      integer n, igz, i1, i2, ig, i
      real xgz(2**nmax), wgz(2**nmax), plz, plz0, plz1 
      real  X, res1, z, rr, y, zpl      
      real, external :: vaint, veint, vaMint
            
      X = r2/r1                 !  X must be less than  1
      
      n = nmax; igz=2**n; i1=igz-1; i2=2*igz-2;
C     n=4; igz=16; i1=15; i2=30;
C     n=5; igz=32; i1=31; i2=62;      
C     n=6; igz=64; i1=63; i2=126;
C      n=7; igz=128; i1=127; i2=254;

      
      xgz(1:igz) = xgaus(i1:i2)
      wgz(1:igz) = wgaus(i1:i2)    
      
      coef = (-1.0)**ll
      coef1 = coef-1.0d0
      
      res0=0.0d0                ! Hydrogen: Znuc = 1
      res1=0.0d0                ! Alkali Atom:  Znuc > 1
      select case (ll)
      case (0)
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaint(rr) + coef*veint(rr) + coef1
C           y = vaMint(rr) + coef*veint(rr) + coef1 ! no exchange: Mitroy == true
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                                                       
            res1=res1+y*wgz(ig) ! plz=1.0d0  ! case (0) 
            res0=res0+wgz(ig)              
         END DO
      case (1)
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaint(rr) + coef*veint(rr) + coef1
C           y = vaMint(rr) + coef*veint(rr) + coef1 ! no exchange: Mitroy == true
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz=zpl             ! case (1)                        
            res1=res1+y*plz*wgz(ig)
            res0=res0+plz*wgz(ig)
         END DO
      case (2)
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaint(rr) + coef*veint(rr) + coef1
C           y = vaMint(rr) + coef*veint(rr) + coef1 ! no exchange: Mitroy == true
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz=0.5d0*(3.0d0*zpl*zpl-1.0d0) ! case (2)                        
            res1=res1+y*plz*wgz(ig)
            res0=res0+plz*wgz(ig)
         END DO   
      case (3)         
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaint(rr) + coef*veint(rr) + coef1
C           y = vaMint(rr) + coef*veint(rr) + coef1 ! no exchange: Mitroy == true
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz=0.5d0*(5.0d0*zpl*zpl-3.0d0)*zpl ! case (3)                        
            res1=res1+y*plz*wgz(ig)
            res0=res0+plz*wgz(ig)
         END DO          
      case default         
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaint(rr) + coef*veint(rr) + coef1
C           y = vaMint(rr) + coef*veint(rr) + coef1 ! no exchange: Mitroy == true
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz0 = 1.0d0; plz1 = zpl; ! case default             
            do i = 2, ll
               plz=(dble(2*i-1)*zpl*plz1 - dble(i-1)*plz0)/dble(i)
               plz0=plz1
               plz1=plz
            end do                                     
            res1=res1+y*plz*wgz(ig)
            res0=res0+plz*wgz(ig)              
         END DO          
      end select            
      
      res0=((2.*ll+1.)/2.)*coef1*res0
      res=((2.*ll+1.)/2.)*res1
      
      return       
      end subroutine Ubb3


      subroutine Ubb3M(ll, r1, r2, res0, res)

C     Same as Ubb3 but without local exchange in Va(r)
      
C========================================================================
C     is used to get                                                    C 
C     
C     res = U (r1,r2) =                                      (Eq. 23)
C            ll                                 
C                (2ll+1)  f 1     _   _         _   _
C              = -------- |  (Ve(|r1+r2|) + Va(|r1-r2|)) P  (z) dz
C                   2     j-1                             ll
C      
C     for alkalies at r1 >= r2 . The integrand is divergent at r1 = r2.
C     P  (z) is Legendre function, res0 is what one has for hydrogen.
C      ll 
C======================================================================
      use apar, only: nmax, lmax1 ! maxx, nmax, lmax1, maxq 
      implicit none    
 
      integer, intent(in) :: ll
      real, intent(in)    :: r1, r2
      real, intent(out)   :: res0, res
      
      real*8, dimension (1:2*(2**nmax-1)) ::  xgaus, wgaus
      real*8, dimension (0:lmax1, 1:2*(2**nmax-1)) :: pleg 
      common /gausc/ xgaus, wgaus, pleg      

      real coef, coef1
      
!     integer i, i1, i2, ig, igz, j, jrange, k, kk,
      integer n, igz, i1, i2, ig, i
      real xgz(2**nmax), wgz(2**nmax), plz, plz0, plz1 
      real  X, res1, z, rr, y, zpl      
      real, external :: vaint, veint, vaMint

      real tc(0:1)
      real dr, r12, r10, y1, y2, y0, dy1

                  
      X = r2/r1                 !  X must be less than  1
      
      n = nmax; igz=2**n; i1=igz-1; i2=2*igz-2;
C     n=4; igz=16; i1=15; i2=30;
C     n=5; igz=32; i1=31; i2=62;      
C     n=6; igz=64; i1=63; i2=126;
C      n=7; igz=128; i1=127; i2=254;

      
c$$$      open(unit=88,file='vave')
c$$$      do i = 0, 100
c$$$         rr = 10.0 ** (-3.0+0.05*real(i))
c$$$         write (88,*) rr,vaMint(rr) + veint(rr)
c$$$      end do
c$$$      stop 'ubb3: see vave'

                  
      xgz(1:igz) = xgaus(i1:i2)      
      wgz(1:igz) = wgaus(i1:i2)
                 
      coef = (-1.0)**ll
      coef1 = coef-1.0d0
              
      res0=0.0d0                ! Hydrogen: Znuc = 1
      res1=0.0d0                ! Alkali Atom:  Znuc > 1
      select case (ll)
      case (0)
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaMint(rr) + coef*veint(rr) + coef1
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                                                       
            res1=res1+y*wgz(ig) ! plz=1.0d0  ! case (0) 
            res0=res0+wgz(ig)              
         END DO         
      case (1)
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaMint(rr) + coef*veint(rr) + coef1
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz=zpl             ! case (1)                        
            res1=res1+y*plz*wgz(ig)
            res0=res0+plz*wgz(ig)
         END DO
      case (2)
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaMint(rr) + coef*veint(rr) + coef1
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz=0.5d0*(3.0d0*zpl*zpl-1.0d0) ! case (2)                        
            res1=res1+y*plz*wgz(ig)
            res0=res0+plz*wgz(ig)
         END DO   
      case (3)         
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaMint(rr) + coef*veint(rr) + coef1
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz=0.5d0*(5.0d0*zpl*zpl-3.0d0)*zpl ! case (3)                        
            res1=res1+y*plz*wgz(ig)
            res0=res0+plz*wgz(ig)
         END DO          
      case default         
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             
            y = vaMint(rr) + coef*veint(rr) + coef1
            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz0 = 1.0d0; plz1 = zpl; ! case default             
            do i = 2, ll
               plz=(dble(2*i-1)*zpl*plz1 - dble(i-1)*plz0)/dble(i)
               plz0=plz1
               plz1=plz
            end do                                     
            res1=res1+y*plz*wgz(ig)
            res0=res0+plz*wgz(ig)              
         END DO          
      end select

C     factor 1/r1 is taken into account in vmat.f where pos-pos
C     matrix elements are calculated 
      res0=((2.*ll+1.)/2.)*coef1*res0; ! /r1
      res=((2.*ll+1.)/2.)*res1;        ! /r1

c$$$C     --------------------------------------------------------
c$$$!     One can factor out X**n (X = r2/r1) from res for odd ll.
c$$$!     It does not work for even ll.
c$$$      if (MOD(ll,2).eq.1) then
c$$$         n = ll                 ! odd ll         
c$$$      else
c$$$         n = ll + 3             ! even ll         
c$$$      end if     
c$$$      res0 = res0/X**n
c$$$      res  = res/X**n
c$$$C     ------------------------------------------------------
      return       
      end subroutine Ubb3M

      
C======================================================================C
C             ALL ROUTINES FOR REARRANGEMENT MATRIX ELEMENTS           C
C======================================================================C      

      subroutine makeveq (maxr, rmesh)
                                ! (maxr, rmesh, qmesh, veq)
                                !( maxr, ver, rmesh,maxq, qmesh, veq)
C======================================================================
C     calculates the Fourier image ve(q) of ve(r)/r = Ve(r)-1/r at some
C     qmash. Ve(r) is  the positron-core potential.
C
C                     oo
C             4 Pi   f
C     veq  =  -----  | dr sin(q r) ve(r).
C               q    J       
C                     0
C     So,        
C                   oo                           
C                   f     2                 4 Pi  
C     Ve(q) = 4 Pi  | dr r  j (q r) Ve(r) = ----  + veq.
C                   J        0              q*q             
C                   0                                                                        
C     veq is not divegent at q = 
C     
C     QUESTION: What mesh for q should be used?
C======================================================================
      use apar, only: PI, PI4
      use ftps_module, only: nr, nexp
      use ql_module
      implicit none
      
      integer, intent (in) :: maxr ! array sizes of ver    
      real, intent (in) :: rmesh(maxr,3)                                                 
C     real, intent(out) :: qmesh(nr), veq(nr)                             

      integer i, j, ll, nn, mr
      real :: pmax, r
      real, external :: veint

C     DP is required for NumSBT
      integer, parameter :: dp=selected_real_kind(precision(1.0d0))
      real(kind=dp) :: rmin,rmax,rhomin,rhomax,dr,kmin,cf
      real(kind=dp) :: fixedk,skk          
      real(kind=dp), dimension (nr) :: rr, kk, ver1, veq1     
      logical singlek

      common /veqc2/ kk, veq1
      
c     integer, intent (in) :: maxq ! array sizes of veq        
c     real, dimension (maxq):: qmesh, veq ! qmesh-points q(i) = i*dq      
c     real, external :: five
C     real ver(maxr)         ! array of ve(r) = Ve(r)*r - 1 at rmesh
c     real dpmin, pp
      
      pmax = min(2.0*PI/rmesh(1,1),100.0)
      
c$$$      if (.false.) then
c$$$  c     simpson integration
c$$$c       pmin = 1.0E-4      
c$$$         call ftcoef(rmesh,ver,maxr)    
c$$$         dpp = (pmax-pmin)/real(maxq-1)                  
c$$$         do i = 1, maxq
c$$$            qmesh(i) = pmin+(i-1) * dpp       
c$$$            veq(i)=five(qmesh(i))         
c$$$         end do
c$$$  else
      
c     sbt method
      ll = 0                    ! since we have integrand with bessel(0,r)
      nn = 2                    ! r**(-nn) - expected integrand behavior 
                                ! for large r
      mr = maxr                 ! maxr = mehsr            
!     setup new rmesh and kmesh parameters
         
      rmin = rmesh(1,1);
      rmax = rmesh(mr,1);      
      rhomin = log(rmin);
      rhomax = log(rmax)         
      dr = (rhomax-rhomin)/(nr-1)
      kmin = log(pmax)-rhomax+rhomin
      
!     set up the r and k space meshes
      rr(1)  = exp(rhomin)
      kk(1)  = exp(kmin)
      ver1(1) = veint(rr(1))/rr(1) 
      cf = exp(dr)

      do i = 2,nr
         rr(i)  = cf*rr(i-1)
         kk(i) = cf*kk(i-1)
         r = real(rr(i))                                    
         ver1(i) = veint(r)/rr(i)           
      end do

c     specify data for ql_module
c     maxql1 = nr;maxql2 = nr;
c      qlgmin = kmin; qlgmax = log(kk(nr));
c      dqlg = dr;                ! = (qlgmax-qlgmin)/dble(maxql1-1);

C     I am not sure that NumSBT should be used
C      maxql1 = 400; maxql2 = maxql1;
C      qlgmin = -8.0; qlgmax = 4.5;
C      dqlg = (qlgmax-qlgmin)/dble(maxql1-1);
         
      singlek = .false.
c                  i    o   i  i   i i   i   i     i    o    i
      call NumSBT(ver1,veq1,rr,kk,ll,nn,nexp,nr,fixedk,skk,singlek)         
      veq1(:) = PI4 * veq1(:)
      
C      veq(:) = veq1(:)
C      qmesh(:) = kk(:)

c$$$  end if
      
      return
      end subroutine makeveq
      
      
      real function five(q)
C     (q, ver, rmesh, maxr)
C=======================================================================
C     calculates the Fourier image of Ve(r)-1/r (simpsons integration)

C     INPUT -----------------------------------------------------------
      use apar, only : znuc1, PI4 
      include 'par.f'
      real q
      real, dimension (maxr,3) ::  rmesh
      real, dimension (maxr) :: veold
      common/meshrr/ meshr, rmesh
      real, dimension(3) :: c      
      common /veqc/ c      
      common /fivec/ veold      
                  
C     LOCAL VARIABLES -------------------------------------------------
      real, parameter :: delta = 1.E-10
      real *8 r, qr
      integer j, meshr
      
      five = 0.0d0
  
C     small q: Taylor expansion: 
      if (q.lt.delta) then 
         five = c(1)+c(2)*q**2+c(3)*q**4 
         return
      end if
      
C     0 < q < 100 : numerical integration 
      if (q.lt.100.0) then    
         do j = 1, meshr           
            five = five + sin(q*rmesh(j,1))*veold(j)*rmesh(j,3)
         end do
         five = PI4 * five/q
         return                  
      end if
                  
C     large q: analytical approximation:
         five = Pi4*znuc1/q/q       
      return
      end

      subroutine ftcoef(rmesh,ver,maxr)
      use apar, only: PI4
      implicit none
      real rmesh(maxr,3)
      real, dimension (maxr) :: ver
      real, dimension(3) :: c
      integer j, maxr
      common /veqc/ c
      real tmp, r2
      
      c(:) = 0.0
      do j = 1, maxr
C         r2=rmesh(j,1)*rmesh(j,1)
         tmp = ver(j)*rmesh(j,3)
         c(1) = c(1) + tmp*rmesh(j,1)
         c(2) = c(2) + tmp*rmesh(j,1)**3/6.
         c(3) = c(3) + tmp*rmesh(j,1)**5/120.
      end do
      c = Pi4 * c

      print*, 'TEST: andrey.f: tcoef: c(1-3):',  (c(j), j=1,3)
      
      return
      end subroutine ftcoef
      
c$$$      real*8 function five1(q)
c$$$C=======================================================================
c$$$C
c$$$C     Now i use the simpson rule to get five and this function is not
c$$$C     needed
c$$$C      
c$$$C     calculates the Fourier image of Ve(r)-1/r with the use of netlib
c$$$C     routine qdawfe.
c$$$C     
c$$$
c$$$C     INPUT -----------------------------------------------------------
c$$$      !implicit none
c$$$      implicit real*8 (a-h,o-z)
c$$$      real *8 q                             
c$$$      real *8, external :: veint ! array of ve(r) = Ve(r)*r - 1 at rmesh              
c$$$      real *8, external :: fun
c$$$      
c$$$C     LOCAL VARIABLES -------------------------------------------------
c$$$      real *8, parameter :: delta = 1.d-16
c$$$      real *8 Pi4, r, qr
c$$$      integer j 
c$$$      
c$$$      Pi4 = 4.d0*ACOS(-1.)
c$$$
c$$$C     check veint
c$$$c$$$      open(unit=80,file = "tmp.dat")
c$$$c$$$      do i =1, 2000
c$$$c$$$         r = 0.1*dble(i-1)+0.001
c$$$c$$$         y1=veint(r)
c$$$c$$$         write (80,*) r, y1
c$$$c$$$      end do
c$$$c$$$      close (80)
c$$$c$$$      stop
c$$$      
c$$$c$$$      call dqawfe(veint,a,omega,integr,epsabs,limlst,limit,maxp1, result
c$$$c$$$     $     ,abserr,neval,ier,result,erlst,ierlst,lst,alist,blist, rlist
c$$$c$$$     $     ,elist,iord,nnlog,chebmo)
c$$$
c$$$      integr = 2
c$$$      result = 0.0
c$$$c$$$  call ft(veint, integr, q, result, ier)
c$$$      
c$$$c$$$  five1 = Pi4*result/q      ! This is what has to be 
c$$$      five1 = Pi4*result/q      ! use this line to check if dqawfe works in five1
c$$$      
c$$$      return
c$$$      end      
c$$$
c$$$      double precision function fun(r)
c$$$      real *8 r
c$$$      fun = dexp(-r)
c$$$      return
c$$$      end 

C                     i    r4      r4       r4  r4 r4    
      function vint(maxx, xmesh, rmesh1max, dx, v, r)
C=====================================================================      
C     this is the electron-core potential for the target in r-space
C
C     this returns the polinomial interpolation of function v(r)
C     given  with the the arraies v and xmesh.  This procedure is
C     used in electroncoreV and positroncoreV routines.

      implicit none
                      
C     integers
      integer, intent(in) :: maxx      
      integer, parameter :: io = 4 ! polinomial-interpolation order
      integer, parameter :: io1 = io-1
      integer j, k

C     SP variables
      real, intent(in) :: r, rmesh1max, dx
      real vint
      
c     DP variables and arrays (due to intrpl)
      real  r1, x
      real*8 :: y
      real*8, dimension (io) :: xa, ya
      real, dimension (maxx), intent(in) :: xmesh, v        
      
      if (r.gt.rmesh1max) then
         vint = 0.0
      else
         x=log10(r)
         j=int((x-xmesh(1))/dx)+1 ! what if x <  xmesh(1) ?          
         k = min(max(j-io1/2,1),maxx-io1) 
         xa=dble(xmesh(k:k+io1))
         ya=dble(v(k:k+io1))

c$$$         print*,'vint:  x: ', x
c$$$         print*,'vint: xa:', (xa(k), k=1,io)
c$$$         print*,'vint: ya:', (ya(k), k=1,io)
         call intrpl(io,xa,ya,1,dble(x),y)
         
         vint = real(y)
      end if

c$$$      print*,'vint: ', r, vint
c$$$      stop
      
      return
      end function vint
    
      function veint(r)
C=====================================================================      
C     this function veint(r) appears in the potential
C      
C                   Ve(r) = [veint(r)+1]/r
C     
C     of positron - alkali ion interaction. One can see that
C     
C              Va(r -> 0) -> Z/r;  Va(r -> oo) -> 1/r
C---------------------------------------------------------------------
C     implicit none
      use apar, only:  maxx, alkali, znuc1,helike
      
      implicit none
C     SP variables      
      real r, veint, dx
      real, external :: vint
      
C     SP arrays
      real, dimension (maxx) :: rmesh1, xmesh, ve, va, vaM

      real na, ca, ne, ce, naM, caM
      common /veva/ na, ca, ne, ce, naM, caM
      
      common /potentials/ xmesh, ve, va, vaM, rmesh1, dx

      if (.not.alkali) then
         veint=0.0
      else
         if (r.lt.1.d-16) then
            veint=znuc1         ! z1 = znucl -1
         else                   ! vint(maxx, xmesh, rmesh1max, dx, v, r)
            if (r.lt.rmesh1(maxx)) then
               veint=vint(maxx, xmesh, rmesh1(maxx), dx, ve, r)
            else
               ! For potassium the tail is due to polpot.
               ! For large r it decreases like this:
               veint = ce/r**ne                             
            end if
         end if
      end if

      if (helike) then
       if (r.lt.1.d-16) then
         veint=znuc1               ! z1 = znucl -1
       else ! vint(maxx, xmesh, rmesh1max, dx, v, r)
         veint=vint(maxx, xmesh, rmesh1(maxx), dx, ve, r)
       end if
      endif  
      
      return
      end function veint
      
      function vaint(r)
C=====================================================================
C     this function veint(r) appears in the potential
C      
C                 Va(r) = [vaint(r)-1]/r ,
C     
C     of active electron - alkali ion interaction. One can see that
C     
C           Va(r -> 0) -> -Z/r;  Va(r -> oo) -> -1/r
C---------------------------------------------------------------------
C     implicit none
      use apar, only:  maxx, alkali, znuc1
      implicit none
      
      real r, vaint, dx
      real, external :: vint
      
      real, dimension (maxx) :: rmesh1          
      real , dimension (maxx) :: xmesh, ve, va, vaM
      real na, ca, ne, ce, naM, caM
      
      common /veva/ na, ca, ne, ce, naM, caM
      common /potentials/ xmesh, ve, va, vaM, rmesh1, dx

      if (.not.alkali) then
         vaint=0.0
      else
         if (r.lt.1.e-16) then
            vaint=-znuc1           
         else
            if (r.lt.rmesh1(maxx)) then
               vaint=vint(maxx, xmesh, rmesh1(maxx),  dx, va, r)
            else
               ! For potassium the tail is due to polpot.
               ! For large r it decreases like this:
               vaint = ca/r**na
            end if
         end if
      end if  
      return
      end function vaint

c$$$      double precision function vaMint(r)
      function vaMint(r)
C=====================================================================      
C     this is the electron-core potential Va = [va(r)-1]/r in r-space
C     implicit none
      use apar, only:  maxx, alkali, znuc1   
      implicit real*8 (a-h,o-z)   

C     SP variables            
      real r, dx, vaMint
      real, external :: vint
      
C     SP arrays
      real, dimension (maxx) :: rmesh1
      real, dimension (maxx) :: xmesh, ve, va, vaM
      real na, ca, ne, ce, naM, caM

      common /veva/ na, ca, ne, ce, naM, caM  
      common /potentials/ xmesh, ve, va, vaM, rmesh1, dx
      
      if (.not.alkali) then
         vaMint=0.0
      else
         if (r.lt.1.d-16) then
            vaMint=-znuc1       ! Has to be corrected
         else
            if (r.lt.rmesh1(maxx)) then
               vaMint=vint(maxx, xmesh, rmesh1(maxx),  dx, vaM, r)
            else
               ! For potassium the tail is due to polpot.
               ! For large r it decreases like this:
               vaMint = caM/r**naM
            end if
         end if
      end if  
      return
      end function vaMint
      
c$$$      subroutine interpolate(io, x, maxx, xtbl, ytbl, result)
c$$$C======================================================================
c$$$C     interpolates some table function given by two arrays xtbl(1:max)
c$$$C     and ytbl(1:max). The point of interpolation MUST satisfy
c$$$C     xtabl(1) <= x <= xtabl(maxx). It is not checked.
c$$$      
c$$$      integer io                ! interpolation order
c$$$      real *8 x                 ! point of interpolation
c$$$      real *8 xtbl(1:maxx)       ! 
c$$$      real *8 ytbl(1:maxx)       !
c$$$      
c$$$      real *8 result, accuracy
c$$$      
c$$$      integer io1, k
c$$$      real dx
c$$$      real *8 dy                ! interpolation accuracy
c$$$      real *8, dimension(io) :: xa, ya
c$$$      
c$$$      io1=io-1
c$$$      dx=xtbl(2)-xtbl(1)      
c$$$      j=int((x-xtbl(1))/dx)+1                
c$$$      k = min(max(j-io1/2,1),maxx-io1) 
c$$$      xa=xtbl(k:k+io1)
c$$$      ya=ytbl(k:k+io1)
c$$$      call intrpl(io, xa, ya, 1, x, result)
c$$$      return
c$$$      end       

      subroutine funleg2(L, q, qa, QL)
C==========================================================================
C     calculated integrals
C
C                      1
C                     f              - -    
C     Q (q,qa) = 2 Pi | dz P (z) V  (q-qa)
C      L              j     L     e
C                     -1
C     where
C        Ve(q) - Fourier image of the positron-core interaction
C                potential Ve(r) = calculated by makeveq, 
C        z = cos(angle_between_q_and_qa)
C        P (z) - Legendre polinomial with  L <= lmax - the maximum 
C          L     value of l defined in ubbgausleg
C         
C     for large q (max[q - qa] > qmax ~ 20-30) the integral is
C     highly oscillatory. The number of points must be increased
C     or another integration method must be used (see five1).
C
C     vefun from /funleg2c/ comes with integration weights
C      
C===========================================================================
                                
      use apar                  ! maxx, nmax, lmax1, maxq, Pi
      implicit real*8 (a-h,o-z)
C      implicit none      
      include 'par.f'
!      include 'par.for'
      
C     INPUT --------------------------------------------------------------
      integer, intent(in)  :: L
      real*8, intent (in)  :: q, qa
      real*8, intent(out)  :: QL
      
      real, external :: five
      real*8, external :: positroncoreV      
      
      real*8, dimension (1:2*(2**nmax-1)) ::  xgaus, wgaus
      real*8, dimension (0:lmax1, 1:2*(2**nmax-1)) :: pleg 
      common /gausc/ xgaus, wgaus, pleg
      
      real rmesh
      common/meshrr/ meshr,rmesh(maxr,3)
      
C     jdouble shows where the intervals  of rmesh are doubled
C     it happens njdouble times
C     common /double/njdouble,jdouble(22)
      real xgz, pl
      dimension xgz(2**nmax), pl(2**nmax)

      real ctail
      real*8 vefun, a1, b1
      !common /funleg2c/ vefun(maxr,2), ctail, maxvdc, a1, b1
      common /funleg2c/ vefun(maxr,2), corep, r0, a1, b1, maxvdc

      !real, dimension (maxr) :: coef
      

C     the second expression (direct integration) ------------------

C     NOTE: this works well for l <=2
C     for higher l the accuracy  can be lost on the border where
C     one method is replaced with the second one
      
      QL=0d0

      if (q.eq.qa) then
         q1 = q; q2 = qa
      else
         q1 = min(q,qa); q2 = max(q,qa)
      end if
      
!     if (2.*q1/q2.gt.5.e-2) then
!     if (((q.lt.50.0).and.(qa.lt.50.0)).or.(abs(q-qa).lt.50)) then   
      if (.false.) then    
C     small q & qa: direct integration from 0 to oo
         do i=1,meshr
            r = rmesh(i,1)            
            x = q*r
            call sbessel(x ,L, bes1)               
            x = qa*r
            call sbessel(x, L, bes2)            
C     QL = QL + veold(i) * r * bes1 * bes2 * rmesh(i,3)
            QL = QL + (vefun(i,1)+vefun(i,2)) * bes1 * bes2                       
         end do

c$$$            if ((qa.gt.800.).or.(q.gt.800.)) then
c$$$               open(unit=66,file='vefun')
c$$$               open(unit=67,file='besfun')
c$$$               print*,'funleg2: q, qa: ', q, qa               
c$$$               do i=1,meshr
c$$$                  r = rmesh(i,1)            
c$$$                  x = q*r
c$$$                  call sbessel(x ,L, bes1)               
c$$$                  x = qa*r
c$$$                  call sbessel(x, L, bes2)            
c$$$C     QL = QL + veold(i) * r * bes1 * bes2 * rmesh(i,3)
c$$$                  write(66,*)  r, vefun(i)/rmesh(i,3)
c$$$                  write(67,*)  r, bes1,  bes2
c$$$               end do
c$$$               close(66); close(67)
c$$$               stop
c$$$            end if
            
c$$$         else
                        
c     large q or qa or both: integration(0,r0) + approximation(r0,oo)
            
c$$$c           find r0 at the border of some simpson's interval and close to
c$$$            r0 = float(l+2)*2.0             
c$$$            j=2
c$$$            do while ((r0.gt.rmesh(jdouble(j),1)).and.(j.le.njdouble))
c$$$               j=j+1
c$$$            end do            
c$$$            jj = jdouble(j-1)
c$$$            do while (r0.gt.rmesh(jj,1))
c$$$               jj = jj+1
c$$$            end do
c$$$            if (mod(jj,2).eq.1) jj=jj-1
c$$$            if ((jj.eq.jdouble(j-1)).or.(jj.eq.jdouble(j))) jj=jj-2
c$$$            r0=rmesh(jj,1)
            
c$$$            print*, '(2) ro: ', r0
c$$$            print*, '(3) meshr, jj, rmesh(jj,1): ', meshr, jj, rmesh(jj
c$$$     $           ,1)            
c$$$            print*, '(4) rmesh: ', rmesh(jj-1,1), rmesh(jj,1), rmesh(jj
c$$$     $           +1,1)
c$$$            print*, '(5) dr(jj-1), dr(jj+1): ', rmesh(jj-1,3)/rmesh(jj,3
c$$$     $           ), rmesh(jj+1,3)/rmesh(jj,3)
c$$$            print*, '(6) jdouble(j-1), jdouble(j) ', jdouble(j-1),
c$$$     $           jdouble(j)
c$$$            print*, '(7) rmesh(jdouble(j-1 & j),1) ', rmesh(jdouble(j-1)
c$$$     $           ,1), rmesh(jdouble(j),1)            
c$$$            stop 'funleg2: second choice' 
            
c$$$c          (1) integral from 0 to r0            
c$$$            coef(:) = 1.0; coef(jj) = 0.5
c$$$            do i = jj, 1, -1
c$$$               r = rmesh(i,1)
c$$$               x = q*r;  call sbessel(x ,L, bes1)               
c$$$               x = qa*r; call sbessel(x, L, bes2)            
c$$$               QL = QL + coef(i) * vefun(i) * bes1 * bes2                        
c$$$            end do                                                                                
c$$$c           (2)integral from r0 to oo
c$$$c            * for slow oscillating integrands
c$$$            if (abs(q-qa).lt.10.0) then
c$$$               QlS = 0.0               
c$$$               do i = jj, meshr 
c$$$                  r = rmesh(i,1)
c$$$                  z1 = q*r; z2 = qa*r
c$$$                  FFC=(Fc(l,z1)*Fc(l,z2) + Fs(l,z1)*Fs(l,z2))/2.0
c$$$                  FFS =(Fc(l,z2)*Fs(l,z1) - Fc(l,z1)*Fs(l,z2))/2.0
c$$$                  bespart = (cos(z1-z2)*FFC+sin(z1-z2)*FFS)
c$$$                  QlS=QlS+coef(i)*vefun(i)* bespart
c$$$               end do
                                                                                
c$$$c     (2) integral for slow oscillating integrands
c$$$            QlS = 0.0  
c$$$            if (abs(q-qa).lt.20.0) then                           
c$$$               do i = 1, meshr 
c$$$                  r = rmesh(i,1)
c$$$                  z1 = q*r; z2 = qa*r
c$$$                  FFC=(Fc(l,z1)*Fc(l,z2) + Fs(l,z1)*Fs(l,z2))/2.0
c$$$                  FFS =(Fc(l,z2)*Fs(l,z1) - Fc(l,z1)*Fs(l,z2))/2.0
c$$$                  bespart = (cos(z1-z2)*FFC+sin(z1-z2)*FFS)
c$$$                  QlS=QlS+coef(i)*vefun(i)* bespart
c$$$               end do               
c$$$c$$$            else
c$$$c$$$               r = r0; z1 = q*r; z2 = qa*r
c$$$c$$$               FFC=(Fc(l,z1)*Fc(l,z2) + Fs(l,z1)*Fs(l,z2))/2.0
c$$$c$$$               FFS =(Fc(l,z2)*Fs(l,z1) - Fc(l,z1)*Fs(l,z2))/2.0
c$$$c$$$               bespart = (-sin(z1-z2)*FFC+cos(z1-z2)*FFS)/(q-qa)
c$$$c$$$               QlS = vefun(jj) * bespart/rmesh(jj,3)
c$$$            end if
c$$$c$$$c            * for fast oscillating integrands
c$$$c$$$            r = r0; z1 = q*r; z2 = qa*r
c$$$c$$$            FFC=(Fc(l,z1)*Fc(l,z2)-Fs(l,z1)*Fs(l,z2))/2.0
c$$$c$$$            FFS =(Fc(l,z2)*Fs(l,z1)+Fc(l,z1)*Fs(l,z2))/2.0
c$$$c$$$            bespart = (-sin(z1+z2)*FFC+cos(z1+z2)*FFS)/(q+qa)
c$$$c$$$            QlF = (-1.0)**l * vefun(jj) * bespart/rmesh(jj,3)
c$$$c$$$            
c$$$c$$$  QL = Qlf + QlS + QL
c$$$
c$$$            QL = QlS
c$$$            
c$$$c$$$c     correction taking into account integration from rmax to oo
c$$$c$$$         
c$$$c$$$         rmax = rmesh(meshr,1);         
c$$$c$$$         if (abs(qa-q)*rmax.gt.dble(20*L)) then                 
c$$$c$$$C           Fmax =  8d0*Pi*veold(maxr)*rmax
c$$$c$$$C           vefun(i) = 8d0*Pi*veold(i)*rmesh(i,1)*rmesh(i,3)
c$$$c$$$            Fmax =  vefun(meshr)/rmesh(meshr,3)
c$$$c$$$            call funleg2a(L, q, qa, rmax, Fmax, result)
c$$$c$$$            QL = QL + result
c$$$c$$$         end if
c$$$         end if
                  
      else
         
C     first expression: -------------------------------------------         
         q1 = q*q + qa*qa;
         q2 = 2.0d0 * q * qa;         
C     2**nmax-point gaussian quadrature
         igz=2**nmax; i1=igz-1; i2=2*igz-2;
         xgz(1:igz) = xgaus(i1:i2)
C        wgz(1:igz) = wgaus(i1:i2)
         pl(1:igz) = pleg(l,i1:i2)
                       
C     numerical integration (fourier image of ve is required)
         QL=0d0              
         do i=1,igz                
            dq = sqrt(q1-q2*xgz(i))        
C           QL = QL + five(dq) * pl(i)
            QL = QL + positroncoreV(dq) * pl(i)
         end do
         
         
      end if
      
      return
      end subroutine funleg2
      
C======================================================+==========

c$$$      subroutine makeftps3(pmin, pmax)
c$$$C----------------------------------------------------------------------------
c$$$C     returns the arraies of sperical bessel transforms of the basis w.f. and
c$$$C     their products with potentials which are C  Va(r) = [va(r) - 1]/r or 
c$$$C     or  Vb(r) = 1/r used with alkali or positronium wavefunctions
c$$$C     (see Eq. 48-49 and ftps1).  
c$$$C
c$$$C     (1) integration: NumSBT 
c$$$C     (2) new p-space mesh
c$$$C      
c$$$C     INPUT: pmin, pmax  specify the min and max values of p 
c$$$C----------------------------------------------------------------------------
c$$$      use apar, only:  alkali, maxp => maxpftps
c$$$      use ftps_module
c$$$C      implicit none
c$$$      real, intent(in) :: pmin, pmax
c$$$            
c$$$      include 'par.f' ! maxr, lamax, nnmax, lnabmax  are defined here
c$$$      common/meshrr/ meshr,rmesh(maxr,3)
c$$$      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
c$$$     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
c$$$     >   npbot(0:lamax),nptop(0:lamax)
c$$$      common /psinbc/ enpsinb(nnmax,0:lnabmax),
c$$$     >     psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)     
c$$$      
c$$$      real, external :: vaint
c$$$      logical, external :: positron             
c$$$             
c$$$!      real, dimension (nr, nnmax, 0:lnabmax) :: ftps, ftpsv
c$$$      
c$$$      integer LMIN, LMAX, l, n, mr
c$$$      integer MINL, MAXL, labot, lpbot, latop,lptop
c$$$      real r
c$$$!     real, dimension (nr) :: rr, Rnl, yc, yc1, kk1, kk2
c$$$      logical pos, singlek
c$$$!     real kappa0, kappa1
c$$$
c$$$      real*8, dimension (meshr) :: rold, fold
c$$$      real*8, dimension (nr) :: rr,Rnl,kk1,kk2,YC1,YC
c$$$      real*8 :: rmin, rmax, rhomin, rhomax, dr, cf
c$$$      real*8 :: kappa0,kappa1,kmin,fixedk,skk      
c$$$      
c$$$      print *,'makeftps: started'
c$$$                           
c$$$*     calculate ftps and ftps
c$$$      
c$$$      MINL = min(labot,lpbot);
c$$$      MAXL = max(latop,lptop);
c$$$
c$$$      mr = meshr
c$$$      
c$$$      kappa0 = log(dble(pmin));
c$$$      kappa1 = log(dble(pmax));
c$$$      dkk = (kappa1 - kappa0)/(nr-1)
c$$$      kk2(1) = pmin
c$$$      cf = exp((kappa1 - kappa0)/(nr-1))
c$$$      do i = 2,nr         
c$$$         kk2(i) = cf * kk2(i-1)
c$$$      end do
c$$$
c$$$      kk(:) = kk2(:)
c$$$      
c$$$      do l = MINL, MAXL
c$$$         do n = nabot(l), natop(l)
c$$$            mr = istoppsinb(n,l)
c$$$!     setup new rmesh and kmesh parameters
c$$$            rmin = dble(rmesh(1,1));
c$$$            rmax = dble(rmesh(mr,1));      
c$$$            rhomin = log(rmin);
c$$$            rhomax = log(rmax)
c$$$                      
c$$$            dr = (rhomax-rhomin)/(nr-1)
c$$$            kmin = log(dble(pmax))-rhomax+rhomin
c$$$                  
c$$$!     set up the r and k space meshes
c$$$            rr(1) = exp(rhomin)
c$$$            kk1(1) = exp(kmin)
c$$$            cf = exp(dr)
c$$$            do i = 2,nr
c$$$               rr(i) = cf*rr(i-1)
c$$$               kk1(i) = cf*kk1(i-1)
c$$$            end do
c$$$
c$$$            rold(1:mr) =  dble(rmesh(1:mr,1))
c$$$            fold(1:mr) = dble(psinb(1:mr,n,l)) ! wave function:
c$$$
c$$$C     I AM HERE ================================================
c$$$c            call intrpl(mr,rmesh(1:mr,1),psinb(1:mr,n,l),nr,rr,Rnl)
c$$$            call intrpl(mr,rold,fold,nr,rr,Rnl)            
c$$$            Rnl(1:nr) = Rnl(1:nr)/rr(1:nr)
c$$$
c$$$            
c$$$            singlek = .false.
c$$$            call NumSBT(Rnl,YC1,rr,kk1,l,n,nexp,nr,fixedk,skk ,singlek)            
c$$$            call intrpl(nr,kk1,YC1,nr,kk2,YC)
c$$$            ftps1(1:nr,n,l) = YC(1:nr)/kk2(1:nr)**l
c$$$
c$$$c$$$            print*,'makeftps3:  n, l: ', n, l
c$$$c$$$            print*,' pmin, pmax: ', pmin, pmax
c$$$c$$$            open (unit = 99, file = 'tmp1')
c$$$c$$$            do i = 1, nr
c$$$c$$$c$$$               write (99,'(4(e12.4,3x))') kk(i),ftps1(i,n,l), ftpsv1(i
c$$$c$$$c$$$     $              ,n,l)
c$$$c$$$               write (99,'(5(e12.4,3x))') kk1(i),yc1(i),kk(i),yc(i),
c$$$c$$$     $              ftps1(i,n,l)
c$$$c$$$            end do
c$$$c$$$            close (99)
c$$$c$$$!            stop
c$$$      
c$$$            pos = positron(n,l,npos)
c$$$            if (pos) then                                                        
c$$$               Rnl(1:nr) = Rnl(1:nr)/rr(1:nr)                          
c$$$            else                
c$$$               if (alkali) then                
c$$$                  do i = 1, nr
c$$$                     r = real(rr(i))
c$$$                     Rnl(i) = Rnl(i)*(1.0 - vaint(r))/r
c$$$                  end do
c$$$               else
c$$$                  Rnl(1:nr) = Rnl(1:nr)/rr(1:nr) ! for hydrogen
c$$$               end if
c$$$            end if              ! pos
c$$$            
c$$$            singlek = .false.
c$$$            call NumSBT(Rnl,YC1,rr,kk1,l,n,nexp,nr,fixedk ,skk ,singlek)
c$$$            call intrpl(nr,kk1,YC1,nr,kk2,YC)
c$$$            ftpsv1(1:nr,n,l) = YC(1:nr)/kk2(1:nr)**l
c$$$
c$$$c$$$            print*,'makeftps3:  n, l: ', n, l
c$$$c$$$            print*,' pmin, pmax: ', pmin, pmax
c$$$c$$$            open (unit = 99, file = 'tmp2')
c$$$c$$$            do i = 1, nr
c$$$c$$$               write (99,'(5(e12.4,3x))') kk1(i), yc1(i), kk(i), yc(i),
c$$$c$$$     $              ftpsv1(i,n,l)
c$$$c$$$            end do
c$$$c$$$            close (99)
c$$$c$$$            stop
c$$$
c$$$c$$$            print*,'makeftps3:  n, l: ', n, l
c$$$c$$$            print*,' pmin, pmax: ', pmin, pmax
c$$$c$$$            open (unit = 99, file = 'tmp2')
c$$$c$$$            do i = 1, nr
c$$$c$$$               write (99,'(5(e12.4,3x))') kk1(i), ftps1(i,n,l), ftpsv1(i
c$$$c$$$     $              ,n,l)
c$$$c$$$            end do
c$$$c$$$            close (99)
c$$$c$$$            stop              
c$$$
c$$$c$$$            do i = 1, 20
c$$$c$$$               print *, kk(i,l,n), ftps1(i,l,n), ftpsv1(i,l,n)
c$$$c$$$            end do   
c$$$c$$$            print *            
c$$$
c$$$         end do                 ! n                           
c$$$      end do                    ! l
c$$$         
c$$$      print*,'makeftps: done'
c$$$      return      
c$$$      end subroutine makeftps3
c$$$      

C=================================================================
      
      subroutine makeftps()
C----------------------------------------------------------------------------
C     returns the arrays of sperical bessel transforms of the basis w.f. and
C     their products with potentials which are Va(r) = [va(r) - 1]/r or 
C     or  Vb(r) = 1/r used with alkali or positronium wavefunctions
C     (see Eq. 48-49 and ftps1).  
C----------------------------------------------------------------------------
C     Functions P(p) and U(p) are tabulated. Interpolation is then used.
C
C     not accurate when  q > qcut 
C     
C----------------------------------------------------------------------------
      use apar, only:  pi, alkali ! , maxp => maxpftps
      use ftps_module
      
!      use apar, only:  maxx, alkali, znuc1  
      include 'par.f' ! maxr, lamax, nnmax, lnabmax  are defined here
      include 'par.pos'
      
C     implicit real*8 (a-h,o-z)      

      !real, intent(in) :: qcut, pmax
      logical numericalv
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /double/njdouble,jdouble(22)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >     psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common /laguercoeffs/ cknd(ncmax,ncmax,0:lnabmax),
     $    rlambda(2,0:lnabmax), npsstates(2,0:lnabmax)
      common/numericalpotential/ numericalv, lstoppos
      common/smallr/ formcut,regcut,expcut,fast,match,analyticd,packed
      logical fast, match, speed, analyticd, boxps
      
C     output: max value of the first index >= maxp
      
      !real, dimension (0:maxp) :: pmesh
      !real*8, dimension (0:maxp,nnmax,0:lnabmax) :: ftps, ftpsv
      !common /makeftpsc/ ftps, ftpsv, pmesh
      
C      common /makeftpsc/ pmesh(maxp), ftps(maxp,nnmax,0:lnabmax),
C     $     ftpsv(maxp,nnmax,0:lnabmax)
C     
      external vaint, positron
      integer npos ! pos. princ. quant. num:  npos = n - natop(l) + nptop(l) 
      logical pos, positron            
      real, dimension (maxr) :: vpot

      real*8 res0, res1, term, r, arg, besj, xmax, tail

      dimension pc0(0:ncmax-1,ncmax,0:6), pc1(0:ncmax-1,ncmax,0:6)
      dimension a(ncmax),f0(0:ncmax-1),f1(0:ncmax-1),reg(maxr),
     >   ucentr(maxr)
      complex phase,sigc
      ucentr(:) = 0.0
      
*     define pmesh: equidistant mesh for interpolation
      
c$$$      dp = pmax/real(maxp-1)           
c$$$      do ip = 1, maxp        
c$$$         pmesh(ip) = real(ip-1)*dp
c$$$      end do


c$$$!     linear scale: from pmesh(0)=0 to pmesh(1)=10**(-3)
c$$$!     logarithmic scale: from pmesh(1) to pmesh(maxp) = 10**3  
c$$$      pmesh(0) = 0.0
c$$$      do ip = 1,maxp                            
c$$$         pmesh(ip) = 10d0 ** (-3.0 +  dble(ip) * 6.0/dble(maxp))
c$$$      end do
      
      do i = 1,meshr
         vpot(i) = 1.0          ! for hydrogen
      end do

      vpot(1:meshr) = 1.0
      
      if (alkali) then
         do i = 1,meshr
            vpot(i) = vpot(i) - vaint(rmesh(i,1)) 
         end do
      end if
      eta = 0.0
      ldw = -1
      
*     ATOM and Ps: calculate ftps and ftpsv      
      do la =  labot, latop
         Nl1=npsstates(1,la); rlam1=rlambda(1,la)
c$$$            pos = positron(na,la,npos) ! npos = n - natop(l) + nptop(l)
c$$$            if (pos) cycle
         if (alkali.or.interpol) then !(.not.analyticd.and.nznuc.eq.1)) then
            jlarge = 1
            reglarge = 0.0
            plarge = 0.0
C$omp parallel do default(shared)
C$omp& private (ip,p,res0,res1,i,r,arg,besj,term,pl,reg,phase,sigc)
C$omp& private (jstart,jstop,na,j)
            do ip = 1, maxp+1
               p = pmesh(ip)               
               call regular(la,p*p,eta,ucentr,cntfug(1,la),ldw,rmesh,
     >            meshr,jdouble,njdouble,regcut,expcut,
     >            reg,jstart,jstop,phase,sigc)
c$$$               do j = 1, meshr ! can check against above
c$$$                  call sbessel(rmesh(j,1)*p,la,besj)
c$$$                  reg(j) = besj * rmesh(j,1)*p
c$$$               enddo
               do na = nabot(la), natop(la)            
                  res0=0.0d0; res1=0.d0;                                    
                  do i = 1, istoppsinb(na, la)                     
                     r = rmesh(i,1)
                     arg = r*p                     
                     term = psinb(i,na,la)*reg(i)*rmesh(i,3)/arg
                     res0 = res0 + term * vpot(i)
                     res1 = res1 + term * r             
                  end do        ! i                        
                  pl = p**la
C                 --------------------------
                  ftpsv(ip,na,la) = res0/pl
                  ftps(ip,na,la) = res1/pl
C                 --------------------------              
               enddo            !na
            end do              ! ip
C$omp end parallel do
            
            p = 0.0
            call regular(la,p*p,eta,ucentr,cntfug(1,la),ldw,rmesh,
     >         meshr,jdouble,njdouble,regcut,expcut,
     >         reg,jstart,jstop,phase,sigc)
            do na = nabot(la), natop(la)
               res0=0.0d0; res1=0.d0;
               do i = jstart, istoppsinb(na, la)
                  r = rmesh(i,1)
                  term = psinb(i,na,la)*reg(i)*rmesh(i,3)
                  res0 = res0 + term * vpot(i) / r
                  res1 = res1 + term 
               end do           ! i
               ftpsv(-1:0,na,la) = res0
               ftps(-1:0,na,la) = res1
            enddo               !na
         end if                 ! alkali.or.interpol
C     analytical calculations for Hydrogen only
         if ((nznuc.eq.1).and.analyticd) then
            do na = nabot(la), natop(la)
               do k=1,Nl1
                  a(k)=cknd(k,na,la)                        
               end do                     
               select case(la)
               case(0)
                  call coefsS(la,Nl1,a,f0,f1)
               case(1)
                  call coefsP(la,Nl1,a,f0,f1)
               case(2)
                  call coefsD(la,Nl1,a,f0,f1)
               case(3)
                  call coefsF(la,Nl1,a,f0,f1)
               case(4)
                  call coefsG(la,Nl1,a,f0,f1)
               case(5)
                  call coefsH(la,Nl1,a,f0,f1)
               case(6)
                  call coefsI(la,Nl1,a,f0,f1)
               case(7)
                  call newcoefs7(la,Nl1,a,f0,f1)
               case(8)
                  call newcoefs8(la,Nl1,a,f0,f1)
               case(9)
                  call newcoefs9(la,Nl1,a,f0,f1)
               case(10)
                  call newcoefs10(la,Nl1,a,f0,f1)
               case(11)
                  call newcoefs11(la,Nl1,a,f0,f1)
               case(12)
                  call newcoefs12(la,Nl1,a,f0,f1)
               case(13)
                  call newcoefs13(la,Nl1,a,f0,f1)
               case(14)
                  call newcoefs14(la,Nl1,a,f0,f1)
               case(15)
                  call newcoefs15(la,Nl1,a,f0,f1)
               case(16)
                  call newcoefs16(la,Nl1,a,f0,f1)
               case(17)
                  call newcoefs17(la,Nl1,a,f0,f1)
               case(18)
                  call newcoefs18(la,Nl1,a,f0,f1)
               case(19)
                  call newcoefs19(la,Nl1,a,f0,f1)
               case(20)
                  call newcoefs20(la,Nl1,a,f0,f1)   
               case default
                  stop 'coeffs are missing for la>20'
               end select       ! la               
               do k=0,Nl1-1
                  pc0(k,na,la)=f0(k)
                  pc1(k,na,la)=f1(k)
               end do           ! k                
               do ip = 0, maxp
                  p = pmesh(ip)                                   
                  p2 = p*p            
                  call f0zpart(Nl1,rlam1,1.0d0,na,la,na,p2,resH0
     $                    ,resH1,pc0,pc1)                  
                  ftpsv(ip,na,la) = resH0
                  ftps(ip,na,la) = resH1                  
               end do           ! ip               
            end do              ! na
         end if                 !analyticd
      end do                    ! la
c$$$      print*,'makeftps: target states: done'
      return

! below is no-longer used because incorporated above, Igor 29/04/2023
      
C     POSITRONIUM FTPS : FOR COMPARISON SINSE FOR THIS CASE
C              & FTPSV : THERE IS ANALYTIC EXPRESSION                   
C        
      do lp =  lpbot, lptop
         Nl2=npsstates(2,lp);
         if (Nl2.eq.0) cycle
         rlam2=rlambda(2,lp);             
         do np = nabot(lp), natop(lp)  
            pos = positron(np,lp,npos)
            if (.not.pos) cycle
c$$$            print*,'lp, npos, n: ', lp, npos, np 
            if (alkali.or.interpol) then!numericalv) then ! direct numerical integration
               if (lp.eq.0) then
                  ipstart = 0
               else
                  ipstart = 1
               endif
!$omp parallel do default(shared)
!$omp& private (ip,p,res0,res1,i,r,arg,besj,term,pl)             
               do ip = ipstart, maxp
                  p = pmesh(ip)   
                  res0=0.0d0; res1=0.d0;                              
                  do i = 1, istoppsinb(np, lp)
                     r = dble(rmesh(i,1))
                     arg = r * p
                     call sbessel(arg,lp,besj)         
                     term = dble(psinb(i,np,lp))*besj*dble(rmesh(i,3))
                     res0 = res0 + term 
                     res1 = res1 + term * r
                  end do        ! i                    
                  pl=p**lp
C                 --------------------------
                  ftpsv(ip,np,lp) = res0/pl
                  ftps(ip,np,lp) = res1/pl
C                 --------------------------
               end do           ! ip
!$omp end parallel do                  
               if (ipstart.eq.1) then
                  ftpsv(0,np,lp) = (pmesh(2)*ftpsv(1,np,lp)-
     >               pmesh(1)*ftpsv(2,np,lp))/(pmesh(2)-pmesh(1))
                  ftps(0,np,lp) = (pmesh(2)*ftps(1,np,lp)-
     >               pmesh(1)*ftps(2,np,lp))/(pmesh(2)-pmesh(1))
               endif

               ip = 0
               do while(pmesh(ip).lt.qcut.and.ip.lt.maxp)
                  ip = ip + 1
               enddo
               ipq = ip
               lp2 = 2*lp+2
               if (ip.lt.maxp) then
                  do while ((pmesh(ip-1)**(lp2)*ftpsv(ip-1,np,lp)-
     >               pmesh(ip)**(lp2)*ftpsv(ip,np,lp))*
     >               (pmesh(ip)**(lp2)*ftpsv(ip,np,lp)-
     >               pmesh(ip+1)**(lp2)*ftpsv(ip+1,np,lp)).gt.0.0
     >               .and.ip.lt.maxp)
                     ip = ip + 1
                  end do
               endif
c$$$               print*,'redefine ftpsv for p:',pmesh(ipq),pmesh(ip)
               do i = ip+1, maxp
                  ftpsv(i,np,lp) = ftpsv(ip,np,lp)
     >               *pmesh(ip)**(lp2)/pmesh(i)**(lp2)
                  ftps(i,np,lp) = ftps(ip,np,lp)
     >               *pmesh(ip)**(lp2)/pmesh(i)**(lp2)
               enddo
c$$$               do ip = 0, maxp
c$$$                  p = pmesh(ip)
c$$$                  write(np*10+lp,*) p,ftpsv(ip,np,lp),ftps(ip,np,lp)
c$$$               enddo 
            else                ! analytical calculations
               do k=1,Nl2
                  a(k)=cknd(k,np,lp)                        
               end do                                                                                             
               select case(lp)
               case(0)
                  call coefsS(lp,Nl2,a,f0,f1)
               case(1)
                  call coefsP(lp,Nl2,a,f0,f1)
               case(2)
                  call coefsD(lp,Nl2,a,f0,f1)
               case(3)
                  call coefsF(lp,Nl2,a,f0,f1)
               case(4)
                  call coefsG(lp,Nl2,a,f0,f1)
               case(5)
                  call coefsH(lp,Nl2,a,f0,f1)
               case(6)
                  call coefsI(lp,Nl2,a,f0,f1)
               case(7)
                  call newcoefs7(lp,Nl2,a,f0,f1)
               case(8)
                  call newcoefs8(lp,Nl2,a,f0,f1)
               case(9)
                  call newcoefs9(lp,Nl2,a,f0,f1)
               case(10)
                  call newcoefs10(lp,Nl2,a,f0,f1)
               case(11)
                  call newcoefs11(lp,Nl2,a,f0,f1)
               case(12)
                  call newcoefs12(lp,Nl2,a,f0,f1)
               case(13)
                  call newcoefs13(lp,Nl2,a,f0,f1)
               case(14)
                  call newcoefs14(lp,Nl2,a,f0,f1)
               case(15)
                  call newcoefs15(lp,Nl2,a,f0,f1)
               case(16)
                  call newcoefs16(lp,Nl2,a,f0,f1)
               case(17)
                  call newcoefs17(lp,Nl2,a,f0,f1)
               case(18)
                  call newcoefs18(lp,Nl2,a,f0,f1)
               case(19)
                  call newcoefs19(lp,Nl2,a,f0,f1)
               case(20)
                  call newcoefs20(lp,Nl2,a,f0,f1)   
               case default
                  stop 'coeffs are missing for lp>20'
               end select       ! lp
               do k=0,Nl2-1
                  pc0(k,np,lp)=f0(k)
                  pc1(k,np,lp)=f1(k)
               end do           ! k                  
c$$$!$omp critical  
c$$$               print*," Nl2, rlam2: ", Nl2, rlam2
c$$$               print*," npos,lp,np: ", npos, lp, np
c$$$!$omp end critical
               
               do ip = 0, maxp
                  p = pmesh(ip)                                   
                  p2 = p*p                                                            
                  call f0zpart(Nl2,rlam2,2.0d0,npos,lp,np,p2,resH0
     $                 ,resH1,pc0,pc1)                  
                  ftpsv(ip,np,lp) = resH0
                  ftps(ip,np,lp) = resH1                  
               end do           ! ip               
            end if            
         end do                 ! np
      end do                    ! lp
      
c$$$      print*,'makeftps: positronium states: done'
c$$$      if(Nl2.le.0) print*,'makeftps: positronium states: done'
      return      
      end subroutine makeftps
                
      subroutine getftps(p, n, l, FT0, FT1)
C---------------------------------------------------------------
C     interpolates ftps and ftpsv for p < pmax
C     
C     pmesh - equidistant mesh 
C     ftps and ftpsv - fourier transforms of ps and ps*v  
C     calculated with makeftps.
C---------------------------------------------------------------      
      use apar !,  maxp => maxpftps
      use ftps_module
      implicit  real*8 (a-h, o-z) !none
      integer, intent(in) :: n,l
      real*8, intent(in)  :: p
      real*8, intent(out) :: FT0, FT1
      include 'par.f'
      integer, parameter :: io = 4, io1=io-1 ! interpolation order 
      real dp, plg
      real*8, dimension(io) :: xa, y0, y1
      real*8 :: dy,res,pd
      integer j,j1,j2
!     parameter maxp = 2000
      
      !real, dimension (0:maxp) :: pmesh
      !real*8, dimension (0:maxp,nnmax,0:lnabmax) :: ftps, ftpsv
      !common /makeftpsc/ ftps, ftpsv, pmesh

      !dp = dble(pmesh(2))    ! pmesh(1) = 0
      !j = int(p/dp)+1  ! find where p could be in the table

!     pmesh(0)=0 < p < pmesh(1)

      !if (error) print*,'getftps: 1'
      if (p.lt.qcut) then
         i = int(p/pmesh(1))
         pm1 =p - pmesh(i-1)
         p0 = p - pmesh(i)
         p1 = p - pmesh(i+1)
         p2 = p - pmesh(i+2)
!     dp2i = 1.0/(pmesh(1)*pmesh(1)) !in module as is dp3i
! below is two-point interpolation
c$$$         FT0o = (ftps(i,n,l)*(pmesh(i+1)-p)+ftps(i+1,n,l)*(p-pmesh(i)))/
c$$$     >      pmesh(1) ! this is just dp
c$$$         FT1o=(ftpsv(i,n,l)*(pmesh(i+1)-p)+ftpsv(i+1,n,l)*(p-pmesh(i)))/
c$$$     >      pmesh(1)            ! this is just dp
! below is three-point interpolation
c$$$         c0 = p1*p2*0.5
c$$$         c1 = p0*p2
c$$$         c2 = p0*p1*0.5
c$$$         FT0o = (ftps(i,n,l)*c0
c$$$     >      -ftps(i+1,n,l)*c1
c$$$     >      +ftps(i+2,n,l)*c2)*dp2i
c$$$         FT1o = (ftpsv(i,n,l)*c0
c$$$     >      -ftpsv(i+1,n,l)*c1
c$$$     >      +ftpsv(i+2,n,l)*c2)*dp2i
! below is four-point interpolation
         cm1 = p0*p1*p2*0.16666666666
         c0 = pm1*p1*p2*0.5
         c1 = pm1*p0*p2*0.5
         c2 = pm1*p0*p1*0.16666666666
         FT0 = dp3i * (
     >      -ftps(i-1,n,l)*cm1
     >      +ftps(i,n,l)*c0
     >      -ftps(i+1,n,l)*c1
     >      +ftps(i+2,n,l)*c2)
         FT1 = dp3i * (
     >      -ftpsv(i-1,n,l)*cm1
     >      +ftpsv(i,n,l)*c0
     >      -ftpsv(i+1,n,l)*c1
     >      +ftpsv(i+2,n,l)*c2)
c$$$!$omp critical 
c$$$         if (abs((FT1o-FT1)/(FT1o+FT1)).gt.0.1.and.
c$$$     >      ftpsv(i,n,l)*ftpsv(i+1,n,l).gt.0.0) then
c$$$            print*,'l,n,p(i):',l,n,i,pmesh(i),p,pmesh(i+1)
c$$$            print*,'FTPSV:',ftpsv(i-1,n,l),ftpsv(i,n,l),ftpsv(i+1,n,l),
c$$$     >         ftpsv(i+2,n,l)
c$$$            print*,'FT1:',FT1o,FT1
c$$$         endif
c$$$!$omp end critical 
      else
         FT0 = ftps(maxp,n,l)*(pmesh(maxp)/p)**(2*l+4)
         FT1 = ftpsv(maxp,n,l)*(pmesh(maxp)/p)**(2*l+2)
      endif
c$$$      write(n*10+l,*) p,ft0,ft1
      return
         
      if (p.le.pmesh(1)) then         
         FT0 = ftps(0,n,l)+(ftps(1,n,l)-ftps(0,n,l))*(p/pmesh(1))
         FT1 = ftpsv(0,n,l)+(ftpsv(1,n,l)-ftpsv(0,n,l))*(p/pmesh(1))
         write(n*10+l,*) p,ft0,ft1
         return
      elseif (p.gt.qcut) then
         FT0 = ftps(maxp,n,l)*(pmesh(maxp)/p)**(2*l+2)
         FT1 = ftpsv(maxp,n,l)*(pmesh(maxp)/p)**(2*l+2)
         write(n*10+l,*) p,ft0,ft1
         return
      end if
            
!     p > pmesh(1)
      plg = log10(p)
      dp = 6.0/real(maxp) 
      irho = int((plg-log10(pmesh(1)))/dp)+1
      krho = min(max(irho-io1/2,1),maxp-io1)-1

       !if (error) print*,'getftps: 2'
      !xrho(1:io) = dble(arholg(krho+1:krho+io))
       xa(1:io) = dble(pmesh(krho+1:krho+io))

       !if (error) print*,'getftps: 3'
      
c$$$      select case(j)
c$$$      case(1)
c$$$         j1 = 1
c$$$         j2 = 4         
c$$$      case(2:maxp-2)
c$$$         j1 = j-1
c$$$         j2 = j+2         
c$$$      case(maxp-1)
c$$$         j1 = j-2
c$$$         j2 = j+1         
c$$$      case (maxp)
c$$$         j1 = maxp-3
c$$$         j2 = maxp 
c$$$      case(maxp+1:)
c$$$         FT0=0.0;  FT1= 0.0
c$$$         return
c$$$      end select

      j1 = krho+1; j2 = krho+io;
      
!      xa = dble(pmesh(j1:j2))
      
      y0 = ftps(j1:j2,n,l)
      y1 = ftpsv(j1:j2,n,l)
      
      call intrpl(io,xa,y0,1,p,res)
      FT0 = res
      call intrpl(io,xa,y1,1,p,res)
      FT1 = res
      write(n*10+l,*) p,ft0,ft1
      return
      end subroutine getftps
      
      subroutine getftps3(p, n, l, FT0, FT1)
C---------------------------------------------------------------
C     interpolates ftps and ftpsv for p < pmax
C     
C     pmesh - equidistant mesh 
C     ftps and ftpsv - fourier transforms of ps and ps*v  
C     calculated with makeftps.
C---------------------------------------------------------------      
      use ftps_module
      implicit none
      integer,intent(in) :: n,l
      real*8, intent (in)  :: p 
      real*8, intent (out) :: FT0, FT1
      
      real*8 res, pd
      real*8, dimension(4) :: xa, y0, y1

      integer j,j1,j2
      integer, parameter :: maxp1 = nr
!     find where p is positioned in the table
      
      j = int((log(p)-log(kk(1)))/dkk) +  1
      if (p.lt.kk(1)) j = 1 
      select case(j)
      case(1)
         j1 = 1; j2 = 4         
      case(2:maxp1-2)
         j1 = j-1; j2 = j+2         
      case(maxp1-1)
         j1 = j-2; j2 = j+1         
      case (maxp1)
         j1 = maxp1-3; j2 = maxp1 
      case(maxp1+1:)
         FT0=0.0d0;  FT1= 0.0d0
         return
      end select      
      xa = dble(kk(j1:j2))
      y0 = ftps(j1:j2,n,l)
      y1 = ftpsv(j1:j2,n,l)
      call intrpl(4,xa,y0,1,p,res)
      FT0 = res
      call intrpl(4,xa,y1,1,p,res)
      FT1 = res
      return
      end subroutine getftps3

      subroutine ftps2(p, n, l, res0, res1)
C-------------------------------------------------------------------
C    calculates spherical bessel transforms: 
C                           
C                      oo      
C                      f     2
C              F (Q) = | dX X  j (Q X) F(X) ,
C               L      j        L      
C                      0
C      
C     where F(r) is r * R(r) = P(r) => res1 or V(r) * P(r) => res0. 
C     R(r) is either target or positronium radial wavefunction and 
C     V(r)= (va(r)-1)/r [or 1/r] for target [or positronium] states
C     [see Eqs. (48-49)].   

      use apar, only : alkali 
      include 'par.f' ! maxr, lamax, nnmax, lnabmax  are defined here
C     implicit real*8 (a-h,o-z)            
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >     psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)     
      
      external vaint, positron           
      logical pos, positron
!      logical alkali
!      common /alk/ alkali
      real, dimension (maxr) :: vpot

C     Taylor-series expansion coefficients of sph. bes. functions:
      real cbesj (0:10)
      data cbesj / 1.000000000000000, 0.3333333333333333, 
     $     0.06666666666666667, 0.009523809523809524, 
     $     0.001058201058201058, 0.00009620009620009620, 
     $     7.400007400007400d-06, 4.933338266671600d-07, 
     $     2.901963686277412d-08, 1.527349308567059d-09,
     $     7.273091945557423d-11 /
                
      pos = positron(n,l,npos)
      if (pos) then
         vpot(:) = 1.0          ! positronium
      else 
         do i = 1, istoppsinb(n,l) ! atom
            rr = rmesh(i,1)
            vv = vaint(rr)           
            vpot(i) = 1.0 - vv           
         end do
      end if
          
      res0=0.0d0; res1=0.d0;      
C     p = 0 -------------------------------------------------
      if (p.lt.1.0d-16) then         
         do i = 1, istoppsinb(n,l)
            r = rmesh(i,1)
            besj = r ** l
            term = psinb(i,n,l) * besj * rmesh(i,3)
c$$$        res0 = res0 + term * (1.d0-vaint(r))
            res0 = res0 + term * vpot(i)
            res1 = res1 + term * r
         end do          
         res0 = res0 * cbesj(l) ! ftpsv(1,na,la)
         res1 = res1 * cbesj(l) ! ftps(1,na,la) 
         return
      end if                    ! p = 0            
c     p =!= 0  --------------------------------------------      
      do i = 1, istoppsinb(n,l)
         r = rmesh(i,1)
         arg = r*p
         call sbessel(arg, l, besj)         
         term = psinb(i,n,l) * besj * rmesh(i,3)
c     res0 = res0 + term * (1.d0-vaint(r))
         res0 = res0 + term * vpot(i)
         res1 = res1 + term * r         
      end do                    ! i      
      pl = p ** l
      res0 = res0/pl            ! ftpsv(ip,n,l) 
      res1 = res1/pl            ! ftps(ip,n,l)          
      return      
      end subroutine ftps2      

*------------------------------------------------------
         
      subroutine phenpots(icase, r, va, ve)
c=======================================================================
c     Here I calculate va(r) and ve(r) for different phenomenological
c     potentials used in the literature. va(r) and ve(r) appear in the
c     electron-core and positron-core potentails Va and Ve:
      
C     Va(r) = [va(r) - 1]/r  and  Ve(r) = [1+ve(r)]/r

C     the notation is slightly confusing but it was used by Alisher,
c     where a of va stands for alpha.
c      
c     now only for sodium
c      
c     A. Schweizer-Fabbinder-Gonzales-Ferez (SFG) potential
c        for electron-core [1,2] and positron-core [2] interactions.
c     -  for sodium and potasium
c     -  introduced by Klapisch in 1969
      
c     References:
c     -------------------------------------------------
C     [1] Han etal, PRA 77, 012721, 2008
c     [2] Le etal, PRA 71, 032713, 2005
c     [3] Kar & Ho, Eur. Phys. J. D 35, 453, 2005
c     [4] Hanssen et al, JPB 12, 899, 1979
c     [5] Schweizer etal, At. Data, Nucl. Data Tab. 72, 33, 1999
c     [6] Klapisch, PhD Thesis 1969 (see Hanssen [4])

      integer, intent(in) :: icase
      real, intent(in)  :: r
      real, intent(out) :: va, ve      
            
      select case (icase)
      case (0)                  ! NO CORE
         va = 0.0
         ve = 0.0
      case (1)                 ! SFG + POLPOT
         a1 = 3.32442452; a2 = 0.71372798; a3 = 1.83281815; 
         ac = 0.9457; rc = 0.52450638;
         W3 = 1.0 - exp(-(r/rc)**3)
         Vp = - ac*W3**2/(2.0d0*r**3)
      case (2)                 ! SFG
         a1 = 7.902; a2 = 23.51; a3 = 2.688;
         Vp = 0.0
      case (3)                  ! HEWITT
         a=3.23; b = 17.3;      ! sodium
C        a=2.22693; b = 3.43;      ! potasium
         va = - exp(-a*r)*(10.0+b*r)
         c = 3.733              ! from Cavaliere & Ferrante
         ve = 10.0*exp(-c*r)
         return         
      end select      
      Vd   = -(10.0*exp(-a1*r) + a2*r*exp(-a3*r))      
      va  = (Vp + Vd)
      ve = (Vp - Vd)      
      return
      end
      
      subroutine check_all(i_check)
C
C     all used tests are collected here 
C

C     case 5
      include 'par.f'           ! to get maxr
      common/meshrr/ meshr,rmesh(maxr,3)
      real, dimension(8) ::  sturm_e, sturm_a
      
      select case(i_check)
      case (0)
         return                 ! no test requested
      case (1)                  ! for vexch
         continue
c$$$            open(99,file='Vel.dat') 
c$$$            do i=1, meshr
c$$$               r = rmesh(i,1)
c$$$               tmp = (temp(i) + vexch(i)+polpot(rmesh(i,1),corep(0)
c$$$     $              ,r0(0)))*r            
c$$$               write (99,'(4(e14.6))') r, tmp
c$$$            end do
c$$$            close (99)
c$$$  stop
         continue
      case (2)                   ! potentials
c$$$            open(99,file='Vel.dat') 
c$$$            do i=1, maxx
c$$$               r = rmesh1(i)
c$$$C               tmp = electroncoreV(maxx, xmesh, dx, va, r)
c$$$               tmp1 = vint(maxx, xmesh, dx, ve, r)
c$$$               tmp2 = veint(r)
c$$$               tmp3 = vint(maxx, xmesh, dx, va, r)
c$$$               write (99,'(3(e20.12))') r, tmp1, tmp3
c$$$            end do
c$$$            close (99)
c$$$  stop
         continue
      case (3)                  ! funleg2               ! 
c$$$ 204        print*,'INPUT Q: '
c$$$            read(*,*) q
c$$$            do ll = 0, 7
c$$$               print *,'ll =', ll
c$$$               write(fname, 203) ll ! test
c$$$ 203           format('dQ_l.',I1) ! test
c$$$               open (unit=99, file=fname)
c$$$C               do nnn1 = -6,1
c$$$C                  do aaa1 = 1., 9., 0.5
c$$$C                     q = aaa1*10.00**nnn1
c$$$                     do nnn = -6,1
c$$$                        do aaa = 1., 9., 0.1
c$$$                           qa = aaa*10.00**nnn
c$$$                           call funleg2(ll, q, qa, QL)
c$$$                           write (99,*) qa, QL
c$$$C                           write (99,*) q, qa, QL
c$$$                        end do
c$$$                     end do
c$$$C                     write (99,*)
c$$$C                  end do
c$$$C               end do
c$$$               
c$$$               close(99)
c$$$            end do
c$$$            print*, 'CHECKING FUNLEG2: SEE dQ_l.n FILES'
c$$$            go to 204
c$$$            stop                  
         continue
      case (4)                  ! Ubb
C     test Ubb -----------------------------------------------
C     i_check_ubb = 0  - no test
C     i_check_ubb = 1  - X-dependence  for given r1
C     i_check_ubb = 2  - r1-dependence  for given X            
C            
C     ERROR: Ubb is calculated inacurately when X=r2/r1 < 1
C     is small and ll > 4.
C     
c$$$            X = 1.; r1 = 1.
c$$$            i_check_ubb = 0     !  select test here 
c$$$            select case (i_check_ubb)
c$$$            case (1)            !     X-dependence  for given r1
c$$$               do while (r1.gt.0.0)
c$$$                  print*,'INPUT r: '; read(*,*) r1;
c$$$                  do ll = 0, 6               
c$$$                     write(fname, '("ubb_",I1)') ll 
c$$$                     open (unit=99, file=fname)                                   
c$$$                     do xlg = -8.,0., 0.01
c$$$                        X = 10.0**xlg                  
c$$$                        r2 = X*r1                                   
c$$$                        call Ubb3(ll, r1, r2, res0, res1)                               
c$$$                        write (99,'(5(e14.7,3x))') X, res0, res1
c$$$                     end do     ! xr               1    2    3
c$$$                     close(99)
c$$$                  end do        ! ll                    
c$$$               end do
c$$$               stop
c$$$            case (2)            !    r1-dependence for given X    
c$$$               do while (X.gt.0.0)
c$$$                  print*,'INPUT X: '; read(*,*) X;
c$$$                  do ll = 0, 6               
c$$$                     write(fname, '("ubb_",I1)') ll 
c$$$                     open (unit=99, file=fname)                
c$$$                     do r1lg = -8.,2., 0.01
c$$$                        r1 = 10.0**r1lg                  
c$$$                        r2 = X*r1                                   
c$$$                        call Ubb3(ll, r1, r2, res0, res1)                               
c$$$                        write (99,'(5(e14.7,3x))') r1, res0, res1
c$$$                     end do     ! xr               1    2    3
c$$$                     close(99)
c$$$                  end do        ! ll                    
c$$$               end do
c$$$               stop
c$$$  end select
         continue
      case (5)

         print*, 'Sturmian basis expansion: '
            
C            open (77, file = 'min.dat')
C           do i = 1, 400
C               alamb = 3.2 +0.001*real(i-1)

         alamb = 3.388          ! best value 
               
C               print*,' lamb = ', alamb
! $omp parallel do default(shared)
! $omp& private(n,j,r,fa,fe)                
         do n=1, 8
            sturm_e(n)=0.0
            sturm_a(n)=0.0
                  do j = 1, meshr
                     r = rmesh(j,1)
                     fa = vaint(r)
                     fe = veint(r)
                     sturm_e(n)=sturm_e(n)+r*sturm(n,alamb,r)*fe*rmesh(j
     $                    ,3)
                     sturm_a(n)=sturm_a(n)+r*sturm(n,alamb,r)*fa*rmesh(j
     $                    ,3)
                  end do
C     print*, n, sturm_e(n), sturm_a(n)                
               end do
! $omp parallel do              
C     print*
               fa2 = 0.
               fe2 = 0.
               x = rmesh(1,1)
               do n = 1, 8               
                  fe2 = fe2 + sturm_e(n)*sturm(n,alamb,x)
                  fa2 = fa2 + sturm_a(n)*sturm(n,alamb,x)                   
               end do
               fa = vaint(x)
               fe = veint(x)
               pmin = (abs(fe2-fe)+abs(fa2-fa))/2.
C               write(77,*) alamb, pmin
C            end do
C            close(77)
            
            
            open(unit=88, file = 'fa.dat')            
            open(unit=99, file = 'fe.dat')

            j = 0; x=0.0            
            do j =1, meshr               
               x = rmesh(j,1)
               fa = vaint(x)
               fe = veint(x)
               fa2 = 0.
               fe2 = 0.
               do n = 1, 2 
                  fe2 = fe2 + sturm_e(n)*sturm(n,alamb,x)
                  fa2 = fa2 + sturm_a(n)*sturm(n,alamb,x)                   
               end do
               write(88,*) x, fa, fa2
               write(99,*) x, fe, fe2               
            end do
            close(88)
            close(99)
            stop 'fe,fa'
                     
      case default
         continue
      end select           
      return
      end subroutine check_all

      
      function sturm(n, lamb, r)
C      
C     eight sturmian functions
C     used as a basis for analytical representation
C     of the potentials Ve and Va
C     OUTCOME: one can do it but the basis itself is no good for
C     calculations of the spherical harmonic expanstion of the
C     potentials (see Ubb)
C     
      integer n
      real lamb, r
      x = lamb * r
      select case (n)
      case(1)
         sturm = 2.* exp(-x)
      case(2)
         sturm = (2.-2.*x)*exp(-x)*sqrt(2.)
      case(3)
         sturm = (3.-6.*x+2.*x*x)*exp(-x)*sqrt(1.3333333333333333)
      case(4)
         sturm = (-4.*(-3.+9.*x-6.*x*x+x**3)* exp(-x))/3.
      case(5)
         sturm = ((15.-60.*x+60.*x*x-20.*x**3+2.*x**4)*exp(-x)*sqrt(0.8
     $        ))/3.
      case(6)
         sturm = (-2.*(-45.+225.*x-300.*x*x+150.*x**3-30.*x**4+2.*x**5)
     $        *exp(-x)*sqrt(0.6666666666666666))/15.
      case(7)
         sturm = ((315.-1890.*x+3150.*x**2-2100.*x**3+630.*x**4-84.*x**5
     $        +4.*x**6)*exp(-x)*sqrt(0.5714285714285714))/45. 
      case(8)
         sturm = (-8.*(-315.+2205.*x-4410.*x*x+3675.*x**3-1470.*x**4 +
     $        294.*x**5-28.*x**6+x**7)*exp(-x)* sqrt(0.5))/315.
      case default
         stop 'more Sturmian functions are required'
      end select
      sturm = sturm*lamb
      return
      end function sturm      

      
      real*8 function FlegQ(l, z)
c------------------------------------------------------     
c     returns the Legendre function of the second kind 
c            
c            f 1
c         1  !      P(l,x) 
c     Q = -  | dx ---------- 
c         2  |     (z - x)
c            J -1     
c             
c     where P(n,x) is the Legendre polinomial 
c     and the argument z is large
c-------------------------------------------------------
      use apar                  ! maxx, nmax, lmax1, maxq 
!     implicit real*8 (a-h,o-z)     
      include 'par.f'
!      include 'par.for'
      real*8, intent(in) :: z      
      real *8 Q, zn
      parameter (jmaxC = 5, lmaxC = 20)
      real*8, dimension (0:lmaxC, 0:jmaxC) :: c
      data c/1.,0.3333333333333333,0.1333333333333333,0
     $     .05714285714285714,0.0253968253968254,0.01154401154401154,0
     $     .005328005328005328,0.002486402486402486,0.001170071758307052
     $     ,0.0005542445170928143,0.0002639259605203878,0
     $     .0001262254593793159,0.00006058822050207163,0
     $     .00002917210616766412,0.00001408308573611371,6
     $     .814396323925989d-6,3.303949732812601d-6,1.604775584508978d-6
     $     ,7.807016357070702d-7,3.803418225239573d-7,1.855325963531499d
     $     -7,0.3333333333333333,0.2,0.1142857142857143,0
     $     .06349206349206349,0.03463203463203463,0.01864801864801865,0
     $     .009945609945609946,0.005265322912381736,0.002771222585464072
     $     ,0.001451592782862133,0.0007573527562758953,0
     $     .0003938234332634656,0.0002042047431736488,0
     $     .0001056231430208528,0.00005451517059140791,0
     $     .00002808357272890711,0.0000144429802605808,7
     $     .416665539217167d-6,3.803418225239573d-6,1.948092261708074d-6
     $     ,9.966983664552936d-7,0.2,0.1428571428571429,0
     $     .09523809523809524,0.06060606060606061,0.0372960372960373,0
     $     .02237762237762238,0.01316330728095434,0.007620862110026197,0
     $     .004354778348586398,0.00246139645789666,0.001378382016422129
     $     ,0.0007657677869011831,0.0004224925720834113,0
     $     .0002316894750134836,0.000126376077280082,0
     $     .00006860415623775879,0.00003708332769608583,0
     $     .00001996794568250776,0.00001071450743939441,5
     $     .731015607117938d-6,3.056541657129567d-6,0.1428571428571429,0
     $     .1111111111111111,0.08080808080808081,0.05594405594405594,0
     $     .0372960372960373,0.02413273001508296,0.01524172422005239,0
     $     .009435353088603863,0.005743258401758873,0.003445955041055324
     $     ,0.002042047431736488,0.001197062287569665,0
     $     .0006950684250404509,0.0004001909113869263,0
     $     .0002286805207925293,0.0001297916469363004,0
     $     .00007321580083586177,0.00004107227851767856,0
     $     .00002292406242847175,0.00001273559023803986,7
     $     .045220131681626d-6,0.1111111111111111,0.09090909090909091,0
     $     .06993006993006993,0.05128205128205128,0.03619909502262443,0
     $     .02476780185758514,0.01651186790505676,0.01076860950329789,0
     $     .006891910082110647,0.004339350792440037,0.002693390147031747
     $     ,0.001650787509471071,0.001000477278467316,0
     $     .0006002863670803894,0.0003569270290748261,0
     $     .0002104954274031026,0.0001232168355530357,0
     $     .00007163769508897422,0.00004139066827362955,0
     $     .00002377761794442549,0.00001358721025395742,0
     $     .09090909090909091,0.07692307692307692,0.06153846153846154,0
     $     .04705882352941176,0.0346749226006192,0.02476780185758514,0
     $     .01722977520527662,0.0117162471395881,0.007810831426392067,0
     $     .00511744127936032,0.003301575018942142,0.002101002284781363
     $     ,0.001320630007576857,0.0008209321668721001,0
     $     .0005051890257674462,0.0003080420888825892,0
     $     .000186258007231333,0.0001117548043387998,0
     $     .00006657733024439136,0.00003940290973647652,0
     $     .00002317818219792737/

      Q = 0d0          
            
c     series expansion              
c         zn = 1/z**(l+1)
         do j = jmaxC, 0, -1                                
            n = l+2*j+1
            zn = exp(-n*log(z)) !1.0d0/z**n
!            print*,'zn:',zn,exp(-n*log(z))
            Q = Q + c(l,j) * zn
         end do         
c$$$    end if                    ! z
      FlegQ = Q      
      return
      end function FlegQ

c$$$        PROGRAM MLQNB
c$$$C
c$$$C       ===============================================================
c$$$C       Purpose: This program computes the Legendre functions Qn(x) 
c$$$C                and Qn'(x) using subroutine LQNB
c$$$C       Input :  x  --- Argument of Qn(x)
c$$$C                n  --- Degree of Qn(x)  ( n = 0,1,úúú)
c$$$C       Output:  QN(n) --- Qn(x)
c$$$C                QD(n) --- Qn'(x)
c$$$C       Examples:     x1 = 0.50,    x2 = 2.50
c$$$C
c$$$C       n      Qn(x1)        Qn'(x1)       Qn(x2)          Qn'(x2)
c$$$C     ----------------------------------------------------------------
c$$$C       0     .54930614    1.33333333   .42364893D+00  -.19047619D+00
c$$$C       1    -.72534693    1.21597281   .59122325D-01  -.52541546D-01
c$$$C       2    -.81866327    -.84270745   .98842555D-02  -.13109214D-01
c$$$C       3    -.19865477   -2.87734353   .17695141D-02  -.31202687D-02
c$$$C       4     .44017453   -2.23329085   .32843271D-03  -.72261513D-03
c$$$C       5     .55508089    1.08422720   .62335892D-04  -.16437427D-03
c$$$C       ===============================================================
c$$$C
c$$$        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$        DIMENSION QN(0:100),QD(0:100)
c$$$        WRITE(*,*)'Please enter Nmax and x '
c$$$        READ(*,*)N,X
c$$$        WRITE(*,40)X
c$$$        WRITE(*,*)
c$$$        WRITE(*,*)'  n          Qn(x)           Qn''(x)'
c$$$        WRITE(*,*)'--------------------------------------'
c$$$        CALL LQNB(N,X,QN,QD)
c$$$        DO 10 K=0,N
c$$$           IF (X.LE.1.0) THEN
c$$$              WRITE(*,20)K,QN(K),QD(K)
c$$$           ELSE
c$$$              WRITE(*,30)K,QN(K),QD(K)
c$$$           ENDIF
c$$$10      CONTINUE
c$$$20      FORMAT(1X,I3,2F17.8)
c$$$30      FORMAT(1X,I3,2D17.8)
c$$$40      FORMAT(3X,'x =',F5.2)
c$$$        END


        SUBROUTINE LQNB(N,X,QN) ! ,QD)
C
C       ====================================================
C       Purpose: Compute Legendre functions Qn(x) & Qn'(x)
C       Input :  x  --- Argument of Qn(x)
C                n  --- Degree of Qn(x)  ( n = 0,1,2,úúú)
C       Output:  QN(n) --- Qn(x)
C                QD(n) --- Qn'(x)
C       ====================================================
C
        IMPLICIT REAL *8 (A-H,O-Z)   
        DIMENSION QN(0:N) ! ,QD(0:N)
        EPS=1.0D-14
        IF (DABS(X).EQ.1.0D0) THEN
           DO 10 K=0,N
 10           QN(K)=1.0D+300
! 10           QD(K)=1.0D+300              
           RETURN
        ENDIF
        IF (X.LE.1.021D0) THEN
           X2=DABS((1.0D0+X)/(1.0D0-X))
           Q0=0.5D0*DLOG(X2)
           Q1=X*Q0-1.0D0
           QN(0)=Q0
           QN(1)=Q1
C           QD(0)=1.0D0/(1.0D0-X*X)
C           QD(1)=QN(0)+X*QD(0)
           DO 15 K=2,N
              QF=((2.0D0*K-1.0D0)*X*Q1-(K-1.0D0)*Q0)/K
              QN(K)=QF
C              QD(K)=(QN(K-1)-X*QF)*K/(1.0D0-X*X)
              Q0=Q1
15            Q1=QF
        ELSE
           QC2=1.0D0/X
           DO 20 J=1,N
              QC2=QC2*J/((2.0*J+1.0D0)*X)
              IF (J.EQ.N-1) QC1=QC2
20         CONTINUE
           DO 35 L=0,1
              NL=N+L
              QF=1.0D0
              QR=1.0D0
              DO 25 K=1,500
                 QR=QR*(0.5D0*NL+K-1.0D0)*(0.5D0*(NL-1)+K)
     &              /((NL+K-0.5D0)*K*X*X)
                 QF=QF+QR
                 IF (DABS(QR/QF).LT.EPS) GO TO 30
25            CONTINUE
30            IF (L.EQ.0) THEN
                 QN(N-1)=QF*QC1
              ELSE
                 QN(N)=QF*QC2
              ENDIF
35         CONTINUE
           QF2=QN(N)
           QF1=QN(N-1)
           DO 40 K=N,2,-1
              QF0=((2*K-1.0D0)*X*QF1-K*QF2)/(K-1.0D0)
              QN(K-2)=QF0
              QF2=QF1
40            QF1=QF0
C           QD(0)=1.0D0/(1.0D0-X*X)
C           DO 45 K=1,N
C45            QD(K)=K*(QN(K-1)-X*QN(K))/(1.0D0-X*X)
        ENDIF
        RETURN
        END

                
c$$$      subroutine funleg3(L, q, qa, QL)
c$$$C     ---------------------------------------------------------------------
c$$$C     same as funleg2 but with tail integral
c$$$C     
c$$$C==========================================================================
c$$$C     calculated integrals
c$$$C
c$$$C                      1
c$$$C                     f              - -    
c$$$C     Q (q,qa) = 2 Pi | dz P (z) V  (q-qa) = 
c$$$C      L              j     L     e
c$$$C                     -1
c$$$C
c$$$C                       oo    
c$$$C                     2 f
c$$$C              = 16 Pi  | dr u (r) r j (q r) j (q  r)   
c$$$C                       j     e       l       l  a
c$$$C                       0
c$$$C     where
c$$$C        Ve(q) - Fourier image of the positron-core interaction
c$$$C                potential Ve(r) = calculated by makeveq, 
c$$$C        z = cos(angle_between_q_and_qa)
c$$$C        P (z) - Legendre polinomial with  L <= lmax - the maximum 
c$$$C          L     value of l defined in ubbgausleg         
c$$$C
c$$$C     vefun from /funleg2c/ comes with integration weights
c$$$C      
c$$$C===========================================================================
c$$$                                
c$$$      use apar                  ! maxx, nmax, lmax1, maxq, Pi      
c$$$      implicit real*8 (a-h,o-z)
c$$$C     implicit none      
c$$$      include 'par.f'
c$$$      include 'par.for'
c$$$      
c$$$C     INPUT --------------------------------------------------------------
c$$$      integer, intent(in)  :: L
c$$$      real*8, intent (in)  :: q, qa
c$$$      real*8, intent(out)  :: QL
c$$$            
c$$$      real (kind = 16), external :: tint
c$$$      real ctail                  
c$$$      real rmesh
c$$$      common /meshrr/ meshr,rmesh(maxr,3)
c$$$             
c$$$      real*8 vefun
c$$$      common /funleg2c/ vefun(maxr), ctail      
c$$$      real (kind = 16) :: DQL, bes1, bes2, x
c$$$      
c$$$      QL=0d0;
c$$$      DQL = 0.0_16
c$$$      do i=1,meshr
c$$$         r = rmesh(i,1)
c$$$         if (r.gt.20.0) cycle
c$$$         x = q*r
c$$$         call sbessel_qp(x ,L, bes1)               
c$$$         x = qa*r
c$$$         call sbessel_qp(x, L, bes2)            
c$$$C        QL3 = QL3 + veint(r) * r * bes1 * bes2 * dr
c$$$         QL = QL + vefun(i) * bes1 * bes2 ! the weight is in vefun
c$$$         DQL = DQL + (bes1*bes2/rmesh(i,1)**2) * rmesh(i,3)                           
c$$$      end do                        
c$$$      QL = QL - 8.00_16 * pi * ctail * (tint(L,qa,q) - DQL)      
c$$$      return
c$$$      end subroutine funleg3

      
      subroutine funleg3(L, q, qa, QL) !, QL1, QL2, part1,part2,part3)
C     ---------------------------------------------------------------------
C     same as funleg2 but with tail integrals
C     
C==========================================================================
C     calculated integrals
C
C                      1
C                     f              - -    
C     Q (q,qa) = 2 Pi | dz P (z) V  (q-qa) = 
C      L              j     L     e
C                     -1
C
C                       oo    
C                     2 f
C              = 16 Pi  | dr u (r) r j (q r) j (q  r)   
C                       j     e       l       l  a
C                       0
C     where
C        Ve(q) - Fourier image of the positron-core interaction
C                potential Ve(r) = calculated by makeveq, 
C        z = cos(angle_between_q_and_qa)
C        P (z) - Legendre polinomial with  L <= lmax - the maximum 
C          L     value of l defined in ubbgausleg         
C
C     vefun from /funleg2c/ comes with integration weights
C      
C===========================================================================
C        QL3 = QL3 + veint(r) * r * bes1 * bes2 * dr
      
      use apar                  ! maxx, nmax, lmax1, maxq, Pi      
      implicit real*8 (a-h,o-z)
C     implicit none      
      include 'par.f'
!      include 'par.for'
      
C     INPUT --------------------------------------------------------------
      integer, intent(in)  :: L
      real*8, intent (in)  :: q, qa
      real*8, intent(out)  :: QL

      integer i

! PGI compiler didn't like that part
!      integer, parameter :: qp = selected_real_kind(15, 307)
      !integer, parameter :: qp = selected_real_kind(33, 4931)
      
!      real (kind = qp), external :: tint
      real*8 , external :: tint
      real corep, r0                  
      real rmesh
      common /meshrr/ meshr,rmesh(maxr,3)
             
      real*8 vefun, a1, b1
      common /funleg2c/ vefun(maxr,2), corep, r0, a1, b1,maxvdc
! PGI compiler didn't like that part
!      real (kind = qp) :: QL1,QL2,DQL1,DQL2, bes, bes1, bes2, x, r2, r4
!      real (kind = qp) :: part1, part2, part3, qq,qaa
      real*8 QL1,QL2,DQL1,DQL2, bes, bes1, bes2, x, r2, r4
      real*8 part1, part2, part3, qq,qaa
c     defining to allow gcc build
      real SinIntegral

      !real, dimension (0:10000) :: weight
      
      QL1=0d0
      DQL1 = 0d0
      QL2=0d0
      DQL2 = 0d0
      
      i = 1; r = 0.0;

      !print *, "funleg3: 1) starting calculating vdcore part " 
      
!     calculate vdcore contribution        
      !do while (i.le.meshr.and.(r.le.20.0.or.i.le.maxvdc))
      do i = 1,meshr
      !do i = 1,maxvdc
         r = rmesh(i,1)
         r2 = r * r
         x = q * r
         call sbessel_qp(x, L, bes1)               
         x = qa * r
         call sbessel_qp(x, L, bes2)
         bes = bes1 * bes2         
c$$$         if (i.le.maxvdc) then  ! calculate vdcore contribution                      
c$$$            QL1 = QL1 + vefun(i,1) * bes ! the weight is in vefun
c$$$            DQL1 = DQL1 + (a1*exp(-b1*r)*r)*bes * rmesh(i,3)
c$$$         end if                             
         QL1 = QL1 + vefun(i,1) * bes ! the weight is in vefun
!         DQL1 = DQL1 + (a1*exp(-b1*r)*r)*bes * rmesh(i,3)
      end do
C     tail integral:
      !if (.false.) then
      if ((q.gt.1.0).and.(qa.gt.1)) then
         x0=rmesh(meshr,1); DQL1 = -8.0*Pi*a1
     -        *((exp(b1*x0)*(-((b1**2 + (q + qa)**2)*
     -        (b1*Cos((q - qa)*x0) + (q - qa)*Sin((q - qa)*x0))) + 
     -        (-1.0)**l*(b1**2 + (q - qa)**2)*
     -        (b1*Cos((q + qa)*x0) + (q + qa)*Sin((q + qa)*x0))))/
     -        (2.*(b1**2 + (q - qa)**2)*(b1**2 + (q + qa)**2)))
         QL1 = QL1 + DQL1
      end if
      
      !print *, "funleg3: 2) starting calculating polpot part", r0 

!     calculate polpot contribution for the interval [0,0.01]

      qmax = max(qa,q); qmin = min(qa,q)

      if (.false.) then
      !if (qmin.gt.2.0.and.qmax.gt.20.0) then
         do i = 0, n1, 2        ! n must be even
            r = a + h * real(i);
            r2 = r*r; r4 = r2*r2
            x = q*r
            call sbessel_qp(x ,L, bes1)               
            x = qa*r
            call sbessel_qp(x, L, bes2)
            bes = bes1 * bes2
            weight = 2.0         
            if (i.eq.0.or.i.eq.n1) weight = 1.0
            part1 = part1 + (r4/r6) * bes * weight
         end do

! print *, "funleg3: 3) polpot part1 is computing"      
         weight = 4.0
         do i = 1, n1-1, 2      ! n 
            r = a + h * real(i);         
            x = q*r
            call sbessel_qp(x ,L, bes1)               
            x = qa*r
            call sbessel_qp(x, L, bes2)
            bes = bes1*bes2         
            part1 = part1 + (r**4/r6) * bes * weight
         end do
         part1 = part1 * h/3d0

      ! print *, "funleg3: 4) polpot part1 is done", part1 

!     calculate polpot contribution for the interval [0.01, 20]
         n1 = 80000
         a = 0.01;
         nnn = 1
         do while (b.lt.40.0)            
            b = real((nnn*int(qmax)+l+1))*pi/qmax/2;
            nnn = nnn + 1
         end do
         
         h = (b-a)/n1;
         r6 = r0**6
         part2 = 0.0
         do i = 0, n1, 2        ! n must be even
            r = a + h * real(i);
            r2= r*r
            x = q*r
            call sbessel_qp(x ,L, bes1)               
            x = qa*r
            call sbessel_qp(x, L, bes2)
            bes = bes1*bes2
            dr = 2.0         
            if (i.eq.0.or.i.eq.n1) dr = 1.0
            f = (1.-exp(-(r/r0)**6))/r2
            part2 = part2 + f * bes * dr
         end do
      !print *, "funleg3: 5) polpot part2 is computing"
      
         dr = 4.0
         do i = 1, n1-1, 2      ! n 
            r = a + h * real(i);
            r2 = r*r
            x = q*r
            call sbessel_qp(x ,L, bes1)               
            x = qa*r
            call sbessel_qp(x, L, bes2)
            bes = bes1*bes2
            f = (1.0-exp(-(r/r0)**6))/r2
            part2 = part2 + f * bes * dr         
         end do  
         part2 = part2 * h/3.0
      !print *, "funleg3: 6) polpot part2 is done", part2

         !qaa = qa; qq = q
         !x = qmin*b
         !call sbessel_qp(x, L, bes2)
         !f =  (1.0-exp(-(b/r0)**6))/b/b         
        ! part3 = f * bes2 * sin(qmax*b-real(l+1)/2.)/qmax               ! tail
         x0 = b
         part3 =   (4*Cos((q - qa)*x0) + 
     -        (-1)**l*(-(Pi*(q + qa)**3*x0**3) - 4*Cos((q + qa)*x0)) + 
     -        (q - qa)*x0*(-2*Sin((q - qa)*x0) + 
     -        (q - qa)*x0*(Pi*x0*Abs(q - qa) - 2*Cos((q - qa)*x0) + 
     -        2*(-q + qa)*x0*SinIntegral((q - qa)*x0))) + 
     -        2*(-1)**l*(q + qa)*x0*
     -        ((q + qa)*x0*Cos((q + qa)*x0) + Sin((q + qa)*x0) + 
     -        (q + qa)**2*x0**2*SinIntegral((q + qa)*x0)))/
     -        (24.*q*qa*x0**3)
         !part3 = 0.0
         
      !print *, "funleg3: 7) polpot part3 is done", part2      
         Ql2 = part1 + part2 + part3
      else
         n1 = 400 
         a = 0.00001; b = 0.01; h = (b-a)/real(n1);
         r6 = r0 ** 6
         part1 = 0.0

c$$$      weight(0) = 1.0
c$$$      do i = 1, n1-1, 2
c$$$         weight(i) = 4.0
c$$$         weight(i+1) = 2.0
c$$$      end do
c$$$      weight(n1) = 1.0
c$$$
c$$$      do i = 0, n1            ! n must be even
c$$$         r = a + h * real(i);
c$$$         r2 = r*r; r4 = r2*r2; 
c$$$         x = q*r; call sbessel(x, L, bes1)               
c$$$         x = qa*r; call sbessel(x, L, bes2)
c$$$         bes = bes1 * bes2 
c$$$         part1 = part1 + (r4/r6 - 1.0/r2) * bes * weight(i)
c$$$
c$$$         print'(" *** ", i3.3,5(2x,e15.7))',
c$$$     1        i, r, r4/r6 , 1.0/r2, bes1, bes2
c$$$      end do
c$$$      part1 = part1 * h/3.0
      
         do i = 0, n1, 2        ! n must be even
            r = a + h * real(i);
            r2 = r*r; r4 = r2*r2
            x = q*r
            call sbessel_qp(x ,L, bes1)               
            x = qa*r
            call sbessel_qp(x, L, bes2)
            bes = bes1 * bes2
            weight = 2.0         
            if (i.eq.0.or.i.eq.n1) weight = 1.0
            part1 = part1 + (r4/r6 - 1.0/r2) * bes * weight
         end do

! print *, "funleg3: 3) polpot part1 is computing"      
         weight = 4.0
         do i = 1, n1-1, 2      ! n 
            r = a + h * real(i);        
            x = q*r
            call sbessel_qp(x ,L, bes1)               
            x = qa*r
            call sbessel_qp(x, L, bes2)
            bes = bes1*bes2         
            part1 = part1 + (r**4/r6 - 1.0/r/r) * bes * weight
         end do
         part1 = part1 * h/3d0

      ! print *, "funleg3: 4) polpot part1 is done", part1 

!     calculate polpot contribution for the interval [0.01, 20]
         n1 = 100000 
         a = 0.01; b = 30.0; h = (b-a)/n1;
         r6 = r0**6
         part2 = 0.0
         do i = 0, n1, 2        ! n must be even
            r = a + h * real(i);
            r2= r*r
            x = q*r
            call sbessel_qp(x ,L, bes1)               
            x = qa*r
            call sbessel_qp(x, L, bes2)
            bes = bes1*bes2
            dr = 2.0         
            if (i.eq.0.or.i.eq.n1) dr = 1.0
            f = -exp(-(r/r0)**6)/r2
            part2 = part2 + f * bes * dr
         end do
      !print *, "funleg3: 5) polpot part2 is computing"
      
         dr = 4.0
         do i = 1, n1-1, 2      ! n 
            r = a + h * real(i);
            r2 = r*r
            x = q*r
            call sbessel_qp(x ,L, bes1)               
            x = qa*r
            call sbessel_qp(x, L, bes2)
            bes = bes1*bes2
            f = -exp(-(r/r0)**6)/r2
            part2 = part2 + f * bes * dr         
         end do  
         part2 = part2 * h/3.0
      !print *, "funleg3: 6) polpot part2 is done", part2

         qaa = qa; qq = q
         part3 = tint(L,qaa,qq)
      !print *, "funleg3: 7) polpot part3 is done", part2
      
         Ql2 = part1 + part2 + part3
      end if
      

      QL2 = 8.0*pi*(-corep/2.0) * QL2 ! tail is taken into account
      
c$$$         if (r.le.20.0) then    
c$$$            QL2 = QL2 + vefun(i,2) * bes ! the weight is in vefun
c$$$            s = min(qa,q)/max(qa,q)            
c$$$            if (s.gt.1.e-6) then
c$$$               DQL2 = DQL2 + (bes/r2) * rmesh(i,3) ! Direct integration
c$$$            else
c$$$               x = r            ! 
c$$$               call sbessel_qp(x ,L, bes1)               
c$$$               q2 = max(q,qa)
c$$$               x = s * r
c$$$               call sbessel_qp(x, L, bes2)
c$$$               bes = bes1*bes2
c$$$               DQL2 = DQL2 + q2*(bes/r2) * rmesh(i,3)
c$$$            end if
c$$$         end if                 
c$$$         i = i+1      
c$$$      end do
      
!      arg = (b1*b1+q*q+qa*qa)/(2.0*q*qa)
!      call funleg(arg,L,DQ0T,DQLT)

c$$$      qInt1 = QL1
c$$$      tail1 =  8.00_16 * pi * ((a1/(2.0*qa*q))*DQLT-DQL1)
c$$$      !tail1 = DQL1/a1
c$$$      qInt2 = QL2
c$$$      tail2 = - 8.00_16 * pi * ctail * (tint(L,qa,q) - DQL2)
c$$$      !tail2 = DQL2
            
      !QL1 = QL1  + 8.0 * pi * ((a1/(2.0*qa*q))*DQLT-DQL1)            
c$$$      QL2 = QL2 - 8.00_16 * pi * ctail * (tint(L,qa,q) - DQL2)      
      QL = QL1 + QL2
      
      return
      end subroutine funleg3
      

      function tint(l,qa,q)
C     integral from 0 to oo of j_l(q r)*j_l (q_a r)/r**2                                    
C     used for calculation of Q_l(q,qa) in funleg3.
      integer l
      !real :: qa, q
      integer, parameter :: qp = selected_real_kind(15, 307)
      !integer, parameter :: qp = selected_real_kind(33, 4931)
      real (kind = qp) :: qa, q, x, qmax, tint        
      real (kind=qp), parameter ::
     1     pi = 3.1415926535897932384626433832795	
      qmax = max(q,qa);
      x = min(q,qa)/qmax;
      tint = (pi*qmax/4d0) * x**l
     1     * (3.0+2.0*real(l)+(1.0-2.0*real(l))*x*x)
     2     /(2.0*real(l)+3.0)/(2.0*real(l)+1.0)/(2.0*real(l)-1.0)
	return
	end
	

      function SinIntegral(x)
!     Asymptotic expansion of Si(x) for rough evaluation        
      real, parameter :: pi = 3.1415926535897932384626433832795
      SinIntegral= Pi/2. + (87178291200.0/x**15 - 479001600/x**13 + 
     -     3628800/x**11 - 40320/x**9 + 720/x**7 - 24/x**5 + 
     -     2/x**3 - 1/x)*Cos(x) + 
     -     (-6227020800.0/x**14 + 39916800/x**12 - 362880/x**10 + 
     -     5040/x**8 - 120/x**6 + 6/x**4 - x**(-2))*Sin(x)
      return
      end




C ==============     

      subroutine savepsold(LL)
C      subroutine saveps(LL, alf, index)
      
C     save the ps-system information in a file 'info.par'
C     which is included in aspa.f
C
C     LL  - orbital quantum number      
C     alf - exponential fall-off parameter of the Laguerre basis
C           of given L
C     index = 1  - saves files in P
C     index = 2  - saves files in P+Q
C     
C     All I need is eigenvectors and eigenvalues which are defined with 
C     cknd(ncmax,ncmax,0:lnabmax) and enpsinb(nnmax, 0:lnabmax)
C
C     Note:
C
C     All information is saved in the following files:
C     - ps/info.par
C     - ps/rmesh.dat
C     - ps/n.f
C     - ps/E_lxx_nyy.dat, xx and yy are max values of l and n
C     - ps/ps_lxx_nyy.dat, xx and yy are values of l and n
C     - ps/cknd - coefficients needed to get ps-functions
C          as an expansion over the Laguerre basis
      
      include 'par.f'
      real*8  cknd, rlambda
      common/meshrr/ meshr, rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop, ntype
     $     ,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     $     npbot(0:lamax),nptop(0:lamax)
      character *64 fname
      character *6 dir
      character *6 path
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >     psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)      
      common /laguercoeffs/ cknd(ncmax,ncmax,0:lnabmax),
     $     rlambda(2 ,0:lnabmax), npsstates(2,0:lnabmax)

      real, parameter :: Ry = 13.6058
      logical, external ::  positron
      logical exists

      dimension grid(maxr), pl(nnmax), f(maxr, nnmax)
      
      
!     let unix take care of  folder ps where all ps-info is saved
C      call system( 'if [ -d ps ]; then echo " *** directory ps exist
C     $s"; else (mkdir ps; echo " *** directory ps has been created"); 
C     $  fi')


      print*, 'savepsold: la, latop, natop(l)', la, latop, natop(l)
      
c$$$      select case (index)
c$$$      case (1)
c$$$      call system( 'if [ -d P ]; then echo "directory P exist
c$$$     $     s"; else (mkdir P; echo "directory P has been created"); 
c$$$     $     fi')         
c$$$      case (2)
c$$$      call system( 'if [ -d P+Q ]; then echo "directory ps exist
c$$$     $     s"; else (mkdir P+Q; echo "directory ps has been create d"); 
c$$$     $     fi')         
c$$$      case default
c$$$         stop ' saveps: STOPPED: no case is specified'
c$$$      end select   
                     

!     save the ps-system information in info.par
!        which is used in aspa calculations

      print*, "savepsold: nmax:", natop(0)
      write(dir,'("ps",I3.3,"/")') natop(0)
      write(path,'("ps",I3.3)') natop(0)
      print*, "savepsold: dir:", dir

!     let unix take care of  folder ps where all ps-info is saved
      call system( 'if [ -d '// path //' ]; then echo " *** directory '
     $     // path // ' exists";
     $     else (mkdir '// path //' ; echo " *** directory '
     $     // path // ' has been created"); fi')      
      
C     inquire(file='ps/info.par', exist=exists)
      inquire(file=dir//'info.par', exist=exists)      
      if (.not.exists) then
         open(99,file=dir//'info.par') 
         nmax=0
         do l=labot, latop         
            if (nmax.lt.natop(l)) nmax=natop(l)
         end do
         write (99,5) labot
 5       format (6x, 'integer, parameter :: labot = ', I3)         
         write (99,10) latop
 10      format (6x, 'integer, parameter :: latop = ', I3)
         write (99,20) meshr
 20      format (6x, 'integer, parameter ::  maxr = ', I8)
         write (99,30) nmax
 30      format (6x, 'integer, parameter ::  nmax = ', I3 ) 
         close(99)
      end if                    ! exists

!     save rmesh
         inquire(file=dir//'rmesh.dat', exist=exists)
         if (.not.exists) then
            open (99,file=dir//'rmesh.dat',form='unformatted')           
            write(99) ((rmesh(i,l),i=1,meshr), l=1,3)
            !print*,' ANDREY: maxr,rmesh(1,l):',meshr,(rmesh(1,l),l=1,3)
            ! write(99) ((rmesh(i,l),l=1,3), i=1,maxr)                        
         !open (99,file='ps/rmesh.dat') 
         !do i = 1, meshr
         !   write(99,*) rmesh(i,1),rmesh(i,2),rmesh(i,3)
         !end do         
            close(99)
         end if                 ! exists
         
                     
!     save nabot(l) and natop(l) ---------------------------------- 
      open(99,file=dir//'n.f', ACCESS = 'APPEND')      
      write (99,12) LL, nabot(LL), LL, natop(LL)                 
 12   format (6x, 'nabot(',I2,') = ',I3,'; natop(',I2,') = ',I3,';') 
      close(99)
         
!     save pseudoenergies enpsinb(n,l)/2 given in atomic units!      
c$$$      if (LL.lt.10) then
c$$$         write(fname, 100) LL, natop(LL)
c$$$      else
c$$$         write(fname, 201) LL, natop(LL)
c$$$      end if
C 100  format('_l',I1.1,'_n',I2.2,'.dat')
 201  format('_l',I2.2,'_n',I2.2,'.dat')

      write(fname, 100) LL, natop(LL)
 100  format('_l',I2.2,'_n',I3.3,'.dat')
      
      print*, 'FILE FOR ENERGIES: ' , dir//'E'//fname
      open(99,file = dir//'E'//fname)   
      write(99,'("#   n   l    E(n,l), a.u.    E(n,l), eV ")')
      write(99,'("#---------------------------------------")')
      do n = nabot(LL), natop(LL)
         if (positron(n,LL,npos)) cycle              
         write(99,'(2(2x,I3.3),2(2x,e16.10))') n,LL,enpsinb(n,ll)/2.,
     $        enpsinb(n,LL) * Ry
      end do
      close(99)
      
!     save coefficients cknd -------------------------------------
      open(98,file=dir//'cknd', ACCESS = 'APPEND')
C     write(98,*) istoppsinb(n,l)               
      nmax1 = npsstates(1,LL)   ! for target 
      nmax2 = npsstates(2,LL)   ! for positronium
      print*, '   saveps: nmax1, nmax2', nmax1, nmax2
      write (98,*) '#   L = ', LL , ', nmax1 = ', nmax1
      write (98,*) '# =========================================' 
      do n = nabot(LL), natop(LL)
         if (positron(n, LL, npos)) cycle  
         do k=1, nmax1
            write(98,*) LL, n, k, cknd(k, n, LL)
         end do                 ! k         
      enddo                     ! n
      close(98)

         
c$$$C----------------------------------------------------------------------
c$$$C What is below is Igor's code for ps-states in Laguerre basis (psdlgr)
c$$$C----------------------------------------------------------------------      
c$$$
c$$$C            INPUT:
c$$$C  ----------------------------------      
c$$$C  L - orbital momentum,
c$$$C  alf -  basis parameter,
c$$$C  nmax - number of pseudostates,
c$$$C  nload, nld -  dimensions of arrays,
c$$$C  C    - matrix of eigenvectors,
c$$$C  grid - "R" - grid of "NR" points
c$$$
c$$$!     alf=1.5
c$$$      L  = LL              
c$$$      a2 = 2.0d0 * alf
c$$$      nr =meshr
c$$$      grid(:) = rmesh(:,1)      
c$$$
c$$$      print *, '*** saveps: nmax = ', nmax ! ANDREY
c$$$      
c$$$      L2 = 2 * L
c$$$      nmax = npsstates(1,l) 
c$$$C
c$$$C     Loop by  r-grid
c$$$C
c$$$      cc1 = a2**(L+1)
c$$$c$$$      do i = 1, nr
c$$$      do i = nr, 1, -1
c$$$         r = grid(i)
c$$$         x = a2 * r
c$$$         c1 = dexp(-0.5d0 * x  +  dble(L+1) * log(x))
c$$$C
c$$$C        Define Laguerre polynomials
c$$$C
c$$$         pl1   = 1.0d0
c$$$         pl(1) = pl1
c$$$         pl2   = dble(L2 + 3) - x
c$$$         if (nmax.gt.1) pl(2) = pl2
c$$$         do n = 3, nmax
c$$$            pl3 = ((dble(2*n-1+L2)-x)*pl2 - dble(n+L2)*pl1) /
c$$$     >         dble(n-1)
c$$$            pl(n) = pl3
c$$$            pl1 = pl2
c$$$            pl2 = pl3
c$$$         end do                 ! nmax
c$$$C
c$$$C        Loop by number of states.
c$$$C
c$$$         
c$$$         do n = nabot(l), natop(l)
c$$$            if (positron(n,l,npos)) cycle 
c$$$            
c$$$            f(i,n) = 0.0d0
c$$$c$$$            if (i .le. jmax(n)) then
c$$$            sum = 0.0d0
c$$$            
c$$$            do m = 1, nmax
c$$$!              sum = sum + c(m, n) * dnorm(m) * pl(m)
c$$$               sum = sum + cknd(m, n, L) * pl(m)                  
c$$$            end do              ! m
c$$$            f(i,  n) = c1  * sum
c$$$         end do                 ! n
c$$$      end do                    ! i

!     save pseudofunctions      
      l = LL

      open(70,file=dir//'istoppsinb', access ='append')     
      do n = nabot(l), natop(l)
         if (positron(n,l,npos)) cycle         
         write(70,*) l, n, istoppsinb(n,l)
      end do                    ! n
      close (70)
      
      do n = nabot(l), natop(l)            
         if (positron(n,l,npos)) cycle  
!     write(fname, 200) l,n
!     200        format('_l',I1.1,'_n',I2.2,'.dat')
         
         !if (l.lt.10) then
            write(fname, 100) l,n
         !else
         !   write(fname, 201) l,n
         !end if

         
         open(70,file=dir//'ps'//fname,form='unformatted')
         ! write(70+l*10+n,*) istoppsinb(n,l)         
c         write(70+l*10+n,*) '#  n, l,  En (au),  En (eV) : ', n,l,
c     $        enpsinb(n,l)/2., enpsinb(n,l) * Ry

!         do i = 1, istoppsinb(n,l)
!!            f has to be equal to psinb : works for H but not for alkalies
!!            write(70+l*10+n,*) rmesh(i,1), psinb(i,n,l), f(i,n+nabot(l)-l-1)
!            write(70+l*10+n) psinb(i,n,l)
!         end do

         write(70) (psinb(i,n,l),i=1,istoppsinb(n,l))          
         close(70)

c$$$C     this is for test:
c$$$         if (l.eq.1.and.n.eq.29) then
c$$$C            open(70,file='ps/ps'//fname,form='unformatted')
c$$$C            read (70) (psinb(i,n,l),i=1,istoppsinb(n,l))
c$$$C            close (70)
c$$$            open(70,file='ps/test_ps'//fname)
c$$$            do i = 1, istoppsinb(n,l)
c$$$               write(70,*) rmesh(i,1),  psinb(i,n,l)
c$$$            end do
c$$$            close (70)            
c$$$         end if            
         
      end do

c$$$      !     save pseudofunctions      
c$$$      l = LL
c$$$
c$$$      open(70,file='ps/ps',form='unformatted')
c$$$      do n = nabot(l), natop(l)            
c$$$         if (positron(n,l,npos)) cycle        
c$$$         write(70) istoppsinb(n,l)                  
c$$$         write(70) (psinb(i,n,l),i=1,istoppsinb(n,l))         
c$$$      end do
c$$$      close(70)
      
      end subroutine savepsold
      

C     SOME ROUTINE TO CALCULATE SPHERICAL BESSEL FUNCTIONS
      SUBROUTINE SF32D(X,N,Y,IERR)
      DIMENSION Y(N+2),A(100)
      INTEGER N,IERR,NMAX,J,NM1,NP1,NP2,NNP1,NPP,I,NU,NU2,
     1NU1,NP,LAM,MU,LAMP2,IRET,INT,MIN0
      DOUBLE PRECISION X,Y,A,SYS087,AA,B,ALOG2E,ALPHA,U,
     1DFLOAT,DSIN,DCOS,DSQRT
      DATA  NMAX/30/,SYS087/1.D-14/
      IERR=0
      IF(N.LT.0) GO TO 30
      IF(X.LT.0.D0) GO TO 29
      IF(N.GE.NMAX) GO TO 30
      NM1=N-1
      NP1=N+1
      NP2=N+2
      NNP1=N*NP1
      IF(X-SYS087) 1,3,3
    1 Y(1)=1.D0
      IF(NM1.EQ.0) GO TO 32
      DO 2 J=1,N
      Y(J+1)=0.D0
    2 CONTINUE
      GO TO 32
    3 CONTINUE
      IF(DFLOAT(NNP1)-X*X) 4,8,8
    4 Y(1)=DSIN(X)/X
      Y(2)=Y(1)/X-DCOS(X)/X
      IF(NM1) 7,7,5
    5 DO 6 J=1,NM1
      Y(J+2)=DFLOAT(2*J+1)*Y(J+1)/X-Y(J)
    6 CONTINUE
    7 GO TO 31
    8 AA=.1D0
      B=.35D0
      ALOG2E=48.D0
      NPP=N
      ASSIGN 9 TO IRET
      GO TO 28
    9 NU1=NU
      NP=INT(X-.5D0+DSQRT(B*X*ALOG2E))
      IF(NP-N) 12,12,10
   10 NPP=NP
      ASSIGN 11 TO IRET
      GO TO 28
   11 NU2=NU
      GO TO 13
   12 NU2=100000
   13 NU=MIN0(NU1,NU2,98)
      CONTINUE
      I=NU
      A(NU+2)=0.D0
   14 IF(I) 15,27,15
   15 A(I+1)=X/(DFLOAT(2*I+1)-X*A(I+2))
      IF(A(I+1)-1.D0) 16,16,17
   16 I=I-1
      GO TO 14
   17 LAM=I
      A(I)=1.D0
   18 I=I-1
      IF(I) 19,20,19
   19 A(I)=DFLOAT(2*I+1)*A(I+1)/X-A(I+2)
      GO TO 18
   20 ALPHA=1.D0/((A(1)-X*A(2))*DCOS(X)+X*A(1)*DSIN(X))
      MU=MIN0(LAM+1,NP1)
      DO 21 J=1,MU
      A(J)=ALPHA*A(J)
   21 CONTINUE
   22 CONTINUE
      IF(LAM+1-N) 23,23,25
   23 CONTINUE
      LAMP2=LAM+2
      DO 24 J=LAMP2,NP1
      A(J)=A(J)*A(J-1)
   24 CONTINUE
   25 DO 26 J=1,NP1
      Y(J)=A(J)
   26 CONTINUE
      GO TO 31
   27 A(1)=DSIN(X)/X
      LAM=0
      GO TO 22
   28 U=2.D0*X/DFLOAT(2*NPP+1)
      NU=NPP+INT(ALOG2E*(AA+B*U*(2.D0-U*U)/(2.D0*(1.D0-U*U))
     1))
      GO TO IRET,(9,11)
   29 IERR=65
      GO TO 31
   30 IERR=66
   31 IF(IERR.NE.0) CALL UTSF11(IERR,32)
   32 RETURN
      END

      SUBROUTINE UTSF11(IERR,N)
      IF(N.EQ.28.AND.IERR.EQ.1)GO TO 28
      IF(N.EQ.29.AND.IERR.EQ.1)GO TO 29
      I=IERR-64
      IF(N.EQ.10)GOTO 10
      IF(N.EQ.11)GO TO 11
      IF(N.EQ.13) GOTO 13
      IF(N.EQ.14) GO TO 14
      IF(N.EQ.15) GO TO 15
      IF (N. EQ. 16) GO TO 16
      IF (N .EQ. 17) GO TO 17
      IF (N .EQ. 18) GO TO 18
      IF (N.EQ.19) GO TO 19
      IF(N.EQ.20)GOTO 200
      IF(N.EQ.21) GO TO 21
      IF(N.EQ.22)GO TO 22
      IF(N.EQ.23)GO TO 23
      IF(N.EQ.28)GOTO 282
      IF(N.EQ.29)GO TO 292
      IF(N.EQ.31) GO TO 31
      IF (N .EQ. 32) GOTO 320
      IF(N.EQ.33)GO TO 330
      IF(N.EQ.34) GO TO 340
      IF(N.EQ.35) GO TO 350
      RETURN
   10 PRINT 101
      GO TO (105,106,107),I
  105 PRINT102
      RETURN
  106 PRINT 103
      RETURN
  107 PRINT 104
      RETURN
   11 PRINT 111
      GO TO (105,106,115),I
  115 PRINT 112
      RETURN
   13 IF (IERR .NE. 1) GOTO 135
      PRINT 131
  131 FORMAT(' OTEKA HB MY,YHK SF13D:HEATAHA OKA')
      PRINT 132
  132 FORMAT(' N 1 - APYMEHT AAH OTPATEHM HAEHEM '/
     *'       L=2,BCEH POOEH C COOBAHEM ACOTHOO'/
     *'       HAEH APYMEHTA')
      GOTO 9005
  135 I=IERR-65
      PRINT 136
  136 FORMAT(' OTEKA HB MY,YHK SF13D:ATAHA OKA')
      GOTO(137,138,139,1391),I
  137 PRINT 1392
 1392 FORMAT(' N 2 - HAEHE L MEHE 1  OE 3.'/
     *'       HAEHE SF13D OOEHO PABHM 1.7D308')
      GOTO 9005
  138 PRINT 1393
 1393 FORMAT(' N 3 - HAEHE APYMEHTA PABHO 0.HAEHE SF13D OOEH',
     *'0',
     *' PABHM'/7X,'1.7D308,EC',
     *' L=2, -1.7D308,EC L=1  3')
      GOTO 9005
  139 PRINT 1394
 1394 FORMAT(' N 4 - P L=1 HAEHE APYMEHTA PEBOCXOT MAKCMAH',
     *'0'/7X,
     *' OCTMOE 2.04D03.',
     *' HAEHE SF13D OOEHO PABHM 1.7D308')
      GO TO 9005
 1391 PRINT 1395
 1395 FORMAT(' N 5 - HAEHE APYMEHTA MEHE -2050 P L=1  2.'/
     *'       HAEHE SF13D OOEHO PABHM 1.7D308')
      GO TO 9005
   14 PRINT141
      GO TO(142,143,144,145,146),I
  142 PRINT 147
      RETURN
  143 PRINT148
      RETURN
  144 PRINT 149
      RETURN
  145 PRINT 1490
      RETURN
  146 PRINT 1491
      RETURN
   15 PRINT 151
      GOTO (142,143,152),I
  152 PRINT153
      RETURN
   16 IF (IERR .NE. 1) GO TO 161
  162 PRINT 163
      PRINT 164
      GOTO 9005
  161 PRINT 165
      PRINT 166
      GOTO 9005
   17 IF (IERR .NE. 1) GOTO 171
  172 PRINT 173
      PRINT 174
      GOTO 9005
  171 PRINT 175
      GOTO 9005
   18 IF (IERR .NE. 1) GO TO 181
  182 PRINT 183
      PRINT 184
      GOTO 9005
  181 PRINT 185
      PRINT 186
 9005 RETURN
   19 PRINT 191
      GOTO(195,196,197),I
  195 PRINT 192
      RETURN
  196 PRINT 193
      RETURN
  197 PRINT 194
      RETURN
  200 PRINT 201
      GOTO(202,203),I
  202 PRINT192
      RETURN
  203 PRINT 193
      RETURN
   21 PRINT 210
      GO TO(211,212,211),I
  211 PRINT 213
      RETURN
  212 PRINT 214
      RETURN     
   22 PRINT 220
      GOTO(221,224),I
  224 PRINT 222
      RETURN
  221 PRINT 223
      RETURN
   23 PRINT 230
      GO TO (231,335,232),I
  231 PRINT 233
      RETURN
  232 PRINT 234
      RETURN
   28 PRINT 281
      RETURN
   29 PRINT 291
      RETURN
  282 PRINT 283
      RETURN
  292 PRINT 293
      RETURN
   31 PRINT 311
      RETURN
  320 PRINT 321
      GO TO (322,323),I
  322 PRINT 324
      GOTO 326
  323 PRINT 325
  326 RETURN
  330 PRINT 331
      GO TO(332,333,334,335),I
  332 PRINT 338
      RETURN
  333 PRINT 336
      RETURN
  334 PRINT 337
         RETURN
  335 PRINT 339
      RETURN
  101 FORMAT (' OTEKA HB MY YHK SF10D: ATAHA OKA'/)
  111 FORMAT (' OTEKA HB MY YHK SF11D: ATAHA OKA'/)
  141 FORMAT (' OTEKA HB MY YHK SF14D: ATAHA OKA'/)
  151 FORMAT (' OTEKA HB MY YHK SF15D: ATAHA OKA'/)
  102 FORMAT(' N 1 - YCOBHOE CO L HE PABHO 1  HE PABHO 2.'/7X,
     *'HAEHE YHK OOEHO PABHM 1.7D308'/)
  103 FORMAT(' N 2 - HAEHE APYMEHTA MEHE  PABHO 0.'/7X,
     *'HAEHE YHK OOEHO PABHM 1.7D308'/)
  104 FORMAT (' N 3 - HAEHE APYMEHTA OE 2048.'/7X,
     *'HAEHE YHK OOEHO PABHM 0.'/)
  112 FORMAT (' N 3 - HAEHE APYMEHTA OE 2048.',
     *'HAEHE YHK'/7X,'OOEHO PABHM 0.'/)
  147 FORMAT(' N 1 - HAEHE L MEHE 1  OE 3.'/7X,
     *'HAEHE YHK OOEHO PABHM 1.7D308'/)
  148 FORMAT(' N 2 - HAEHE APYMEHTA MEHE 0 P L=1.'/7X,
     *'HAEHE YHK OOEHO PABHM 1.7D308'/)
  149 FORMAT(' N 3 - HAEHE APYMEHTA OE  PABHO 1',
     *' P L=1  L=2.'/7X,
     *'HAEHE YHK OOEHO PABHM 1.7D308'/)
 1490 FORMAT(' N 4 - HAEHE APYMEHTA OE 1 P L=3.'/7X,
     *'HAEHE YHK OOEHO PABHM 1.7D308'/)
 1491 FORMAT(' N 5 - HAEHE APYMEHTA MEHE  PABHO 0 P L=3.'/
     *7X,'HAEHE YHK OOEHO PABHM 1.7D308'/)
  153 FORMAT(' N 3 - HAEHE APYMEHTA OE 1 P L=1,2,3.'/7X,
     *'HAEHE YHK OOEHO PABHM 1.7D308'/)
  163 FORMAT(' OTEKA HB MY,OPOPAMMA SF16D:',
     *'HEATAHA OKA')
  164 FORMAT (' N 1 - HAEHE APYMEHTA MEHE 0.'/
     *'        HAEH AKER  AKEI OOEH PABHM 1.7D308')
  165 FORMAT(' OTEKA HB MY,OPOPAMMA SF16D:',
     *'ATAHA OKA')
  166 FORMAT(' N 2 - ACOTHOE HAEHE APYMEHTA OE 2.04D03.'/
     *'       HAEH BER  BEI OOEH PABHM 0.',
     *'B CYAE'/7X,'KOA APYMEHT HEOTPATEEH,',
     *' AKER  AKEI TAKE OOEH'/7X,'PABHM 0,',
     *'HAE AKER  AKEI OAATC PABHM 1.7D308')
  173 FORMAT(' OTEKA HB MY,OPOPAMMA SF17D:',
     *'HEATAHA OKA')
  174 FORMAT(' N 1 - HAEHE APYMEHTA MEHE 0.'/
     *'       HAEH AKERD  AKEID OOEH PABHM 1.7D308')
  175 FORMAT(' N 2 - ACOTHOE HAEHE APYMEHTA OE 2.04D03.'/
     *7X,'HAEH BERD  BEID OOEH PABHM 0,',
     *'B CYAE,KOA'/7X,'APYMEHT HEOTPATEEH,',
     *' AKERD  AKEID TAKE OOEH PABHM 0,'/7X,
     *'HAE AKERD  AKEID OAATC PAHM 1.7D308'/)
  183 FORMAT(' OTEKA HB MY OPOPAMMA SF18D:',
     *'HEATAHA OKA')
  184 FORMAT(' N 1 - HAEHE APYMEHTA MEHE 0.'/
     *'       HAEH AKER1  AKEI1 OOEH PABHM 1.7D308')
  185 FORMAT(' OTEKA HB MY,OPOPAMMA SF18D:',
     *'ATAHA OKA')
  186 FORMAT(' N 2 - ACOTHOE HAEHE APYMEHTA OE 2.04D03.'/
     *'       HAEH BER1  BEI1 OOEH PABHM 0.',
     *'B CYAE, KOA'/7X,'APYMEHT HEOTPATEEH,',
     *' AKER1  AKEI1 TAKE OOEH PABHM 0,'/7X,
     *'HAE AKER1,AKEI1 OAATC PABHM 1.7D308')
  191 FORMAT(' OTEKA HB MY, YHK SF19D: ATAHA OKA'/)
  201 FORMAT(' OTEKA HB MY, YHK SF20D: ATAHA OKA'/)
  210 FORMAT('   ,  SF21D:  '/)
  220 FORMAT(' OTEKA HB MY, YHK SF22D: ATAHA OKA'/)
  192 FORMAT(' N 1 - HEBEPHO AAH AKTECK APAMETP RK2: 0<RK2<1'/
     1'       HAEHE YHK OOEHO PABHM 1.7D308'/)
  193 FORMAT(' N 2- HEBEPHO AAH AKTECK APAMETP T: 0<T<1.5707',
     *'96326794896'/
     1'       HAEHE YHK OOEHO PABHM 1.7D308'/)
  194 FORMAT(' N 3 - P RK2=1  T=1.5707963267848 HAMEHATE ',
     *'OHTEPAHO'/7X,'YHK PABEH HY',
     1' HAEHE YHK OOEHO PABHM 1.7D308'/)
  213 FORMAT(' N 1 -      X   ',
     *' '/7X,'  .',
     1'  SF21D   1.7D308')
  214 FORMAT(' N 2 -    21.43, ',
     1'    1.7D308')
  223 FORMAT(' N 1 - TA-YHK HE OPEEEHA P X, Y,  X+Y,',
     1'PABHX'/7X,'HY ',
     *' OMY EOMY OTPATEHOMY CY.',
     1' HAEHE'/7X,'YHK OOEHO PABHM 1.7D308'/)
  230 FORMAT(1X,'OTEKA HB MY,',
     *'OPOPAMMA SF23P; ATAHA OKA'/)
  233 FORMAT(1X,' N 1 - BEECTBEHHA ACT ',
     *'APYMEHTA HE MOET T OE 21.4'/ )
  234 FORMAT(1X,' N 2 - AMMA-YHK ',
     *'HE OPEEEHA P X=0  OMY ',
     *'EOMY'/7X,'OTPATEHOMY CY. HAEHE ',
     1' YHK OOEHO PABHM 1.7D308'/)
  281 FORMAT (' OTEKA HB MY YHK SF28D: HEATAHA OKA'/
     *' N 1 - HAEHE APYMEHTA MEHE 0;'/
     *'       YHK BCETC  ACOTHOO HAEH APYMEHTA'/)
  283 FORMAT (' OTEKA HB MY YHK SF28D: ATAHA OKA'/
     *' N 2 - HAEHE APYMEHTA OE 2070'/
     *'       HAEHE YHK OOEHO PABHM 1.7D308'/)
  291 FORMAT (' OTEKA HB MY YHK SF29D: HEATAHA OKA'/
     *' N 1 - HAEHE APYMEHTA MEHE 0;'/
     *'       YHK BCETC  ACOTHOO HAEH APYMEHTA'/)
  293 FORMAT (' OTEKA HB MY YHK SF29D: ATAHA OKA'/
     *' N 2 - HAEHE APYMEHTA OE 2070'/
     *'       HAEHE YHK OOEHO PABHM 1.7D308'/)
  222 FORMAT (' N2 - HAEHE APMEHTA MEHE -961.199 ',
     1' HAEHE HK'/7X,'OOEHO PABHM 1.7D308'/)
  311 FORMAT(' OTEKA HB MY, YHK SF31D: ATAHA OKA'/
     *' N 1 - HTEPAH KOCHYC HE OPEEEH  HEOTPATEHX'/
     *7X,'APYMEHTOB',
     *' HAEHE YHK OOEHO PABHM 1.7D308'/)
  321 FORMAT(' OTEKA HB MY,OPOPAMMA SF32D:ATAHA OKA')
  324 FORMAT (' N 1 - HAEHE APYMEHTA OTPATEHO.')
  325 FORMAT (' N 2 - ADAHH OPOK CEPECKO YHK',
     *' ECCE HE PHAET'/7X,'HTEPBAY (0,29).')
  331 FORMAT(' OTEKA HB ',
     *'MY, OPOPAMMA SF33D: ',
     *'ATAHA OKA '/)
  336 FORMAT(' N 2 - NMAX (PAMEPHOCT',
     *' MACCBOB FUN1  FUN2) > MAX',
     *'(N,2*[X]+2)'/)
  337 FORMAT(' N 3 - HAEHE APYME',
     *'HTA X OHO T OE 0'/)
  338 FORMAT(' N 1 - N (CO HAEH ',
     *'OPKA,  KOTOPX HYHO OC',
     *'TAT YHK) > 2'/)
  339 FORMAT(' N 4 - HAEHE ',
     *'BCEMO YHK PEBOCX',
     *'OT HAOEE'/7X,'CO, PE',
     *'CTABMOE B MAHE (CM. COH',
     *'TEH APEC)'/)
  340 PRINT 341
  341 FORMAT('   , SFG4D:',
     *' ')
      GO TO 342
  342 GO TO (343,344,345),I
  343 PRINT 346
  346 FORMAT(' N 1 -   X  0,  ',
     *' '/7X,'    ',
     *' ')
      RETURN
  344 PRINT 347
  347 FORMAT(' N 2 -   A    0, ',
     *''/7X,'   ',
     *' '/7X,'  ')
      RETURN
  345 PRINT 348
  348 FORMAT(' N 3 -    ',
     *'  '/7X,'  A  EPS, ',
     *'  '/7X,'  ',
     *'   ')
      RETURN
  350 PRINT 351
  351 FORMAT('   , SFG5D:',
     *' ')
      GO TO 352
  352 GO TO (353,354,355),I
  353 PRINT 356
  356 FORMAT(' N 1 -   X  0,  ',
     *' '/7X,'    ',
     *' ')
      RETURN
  354 PRINT 357
  357 FORMAT(' N 2 -   A    0, ',
     *''/7X,'   ',
     *' '/7X,'  ')
      RETURN
  355 PRINT 358
  358 FORMAT(' N 3 -    ',
     *'  '/7X,'  A  EPS, ',
     *'  '/7X,'  ',
     *'    ')
      RETURN 
      END

c**************************************************************************
c     The aim of this file is to make one-electron functions, calculate
c     the overlap matrix ortint(n,m) and one-electron Hamiltonian matrix e1r(n,m)
c     for these functions 
c     subroutune  setmaxsp - set overall max s.p. functions
c     subroutine setgrids 
c     subroutine inputdatamax
c     subroutine configsp
c     subroutine factorials
c     subroutine spcoef  -  forms  Laguerre polinomials
c     call funlaguer -   forms Sturmian functians from Laguerre polinomials
c     subroutine minmaxi(fl,n,nr,i1,i2)
c     subroutine checkortsq
c     subroutine e1mat
c     subroutine setfang
c     double precision function oneelint(n1,n2)
c---  --------------------------------------------------------------------
c     this subroutine is called to make all necessary one-electron quantities:
c     enumerate s.p. states (call inputdatamax,call configsp),
c     set sq.int.basis (call  factorials,call spcoef),
c     overlap integrals between s.p. states (call checkortsq),
c     one electron integrals (call e1mat).
c     note - called only ones in calculation of HE states of all symmetries
      subroutine  setmaxsp(fl,maxf,minf)
      use vmat_module, only: nodeid
      include 'par.f'
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      dimension  nk1(0:lomax), nkd(0:lomax)
      common /al_array/ al(0:lomax,2)
      double precision  al
      common /nhforb/ nhfo
c     
      call setgrids
      call  factorials
      do l=0,lomax
         nk1(l) = 0
         nkd(l) = 0
      end do
c     read information on one electron orbitals
      if (nodeid.eq.1)
     >write(4,'(" description of s.p. orbitals")') 
      
      call inputdatamax(nk1,nkd,al,ndif)
      call configsp(nk1,nkd,ndif)

      call get_one_e_func(fl,maxf,minf)
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_one_e_func(fl,maxf,minf)
      include 'par.f'
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /al_array/   al(0:lomax,2)
      double precision  al, M(nspmCI)
      
      call spcoef(al,M)
c     if(nhfo.eq.1) then
c     call hforb
c     end if
      kom = 0
      do nsp=1,nspm
         kom = max(kom,ko(nsp))
      end do
      call make_spf(M,fl,maxf,minf,kom)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      subroutine make_spf(M,fl,maxf,minf,kom)
      
      include 'par.f'
      double precision  M(nspmCI)
      common /meshrr/ nr, gridr(nmaxr,3)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      integer  maxg(nspmCI),ming(nspmCI)
      real  gl(nr,nspm)
      double precision  f8(nr,kom)
      
c     from here all s.p. information depend on 'nsp'
      call funlaguer(M,fl,maxf,minf,gl,maxg,ming,f8,kom)
c     if(nhfo.eq.1) then
c     call hforbnum(gl)
c     end if

c     This subroutine calculate overlap integral between one-electron functions and defines the ortint arrays.
      call checkortsq(fl,maxf,minf)

      call e1mat(fl,maxf,minf,gl,maxg,ming,nr,nspm)
      
      return
      end
c-----------------------------------------------------------------
      subroutine setgrids
      include 'par.f'
      common /meshrr/ nr, gridr(nmaxr,3)
      double precision  nr1,grid
      common /gridr/ nr1,grid(nmaxr)
      nr1=DBLE(nr)
      do i=1,nr
         grid(i)=DBLE(gridr(i,1))
      end do
      return
      end
c-----------------------------------------------------------------
c****************  Read overall input data  **********************
c-----------------------------------------------------------------
      subroutine inputdatamax(nk1,nkd,al,ndif)
c     this is to read overall maximum on s.p. functions - for all HE states.
      use vmat_module, only: nodeid
      include 'par.f'
      dimension  nk1(0:lomax), nkd(0:lomax)
      dimension ral(0:lomax), raldif(0:lomax)
      double precision  al(0:lomax,2), Z
      common /Zatom/ Z
      read(3,*) rZ
      if (abs(rz-Z).gt.1e-4) then
         print*, 'Z(from F5) =', rz, ', Z(from ccc.in) = ', real(Z)
         stop 'Z in F5 not same as in ccc.in'
      end if
c--
      do l=0,lomax
         nk1(l) = 0
         nkd(l) = 0
      end do
      l2max = 0
c--
      if (nodeid.eq.1)
     >write(4,'(" description of s.p. states")') 
c     nset=1
      read(3,*) l1max
      if (nodeid.eq.1)
     >write(4,'("lmax =",I5)') l1max
      read(3,*) (ral(l), l=0,l1max)
      if (nodeid.eq.1)
     >write(4,'("al(l) =",10F10.5)') (ral(l), l=0,l1max)
      read(3,*) (nk1(l), l=0,l1max)
      if (nodeid.eq.1)
     >write(4,'("nk(l) =",10I5)') (nk1(l), l=0,l1max)
      call checklomax(l1max,l2max)
c     nset=2
      read(3,*) ndif
      if (nodeid.eq.1)
     >write(4,'("ndif =",I5)') ndif
      read(3,*) ldif
      if (nodeid.eq.1)
     >write(4,'("ldif =",I5)') ldif
      read(3,*) (raldif(l), l=0,ldif)
      if (nodeid.eq.1)
     >write(4,'("aldif(l) =",10F10.5)') (raldif(l), l=0,ldif)
      read(3,*) (nkd(l), l=0,ldif)
      if (nodeid.eq.1)
     >write(4,'("nkd(l) =",10I5)') (nkd(l), l=0,ldif)
      do l=0,lomax
         al(l,1) = DBLE(ral(l))
         al(l,2) = DBLE(raldif(l))
      end do
      return
      end
c-----------------------------------------------------------------
      subroutine configsp(nk1,nkd,ndif)
      include 'par.f'
      dimension  nk1(0:lomax), nkd(0:lomax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      nsp=0
      do l=0,lomax
         n = nk1(l)
         if(n.ne.0) then
            do k=1,n
               nsp=nsp+1
               lo(nsp)=l
               ko(nsp)=k
               nset(nsp)=1
               if(ndif.ne.0.and.k.le.nkd(l)) then
                  nset(nsp)=2
               end if
            end do
         end if
      end do
      nspm=nsp
      if (nspm.gt.nspmax) then
         print*,'Need to increase NSPMAX to at least:',nspm
         stop 'Need to increase NSPMAX'
      endif 
      return
      end
c-------------------------------------------------------------------
      subroutine factorials
c     n! = gamma(n+1)
      include 'par.f'
      double precision  gamma, Cgr
      common /factorial/ gamma(maxfac)
      common /factratio/ Cgr(0:lomax,komax)
c      write(*,'("set gamma")')
      gamma(1) = 1.D0
      do i=2,maxfac
         ii=i-1
         gamma(i) = DBLE(ii)*gamma(ii)
      end do

      do k=1,komax
         do l=0,lomax
            Cgr(l,k) = k
            do i=1,2*l+1
               Cgr(l,k) = Cgr(l,k) *dble(k+i)
            end do
         end do
      end do

c$$$      do k=1,komax
c$$$         do l=0,lomax
c$$$            print*, l,k,Cgr(l,k),gamma(2*l+k+2)/gamma(k)
c$$$         end do
c$$$      end do
     
      return
      end

c-----------------------------------------------------------------
c**************       Set s.p. basis              ****************
c-----------------------------------------------------------------
      subroutine spcoef(al,M)
c     This routine generates normalisation factors for Laguerre polinomials.
      use vmat_module, only: nodeid
      include 'par.f'
      double precision  als,M(nspmCI),al(0:lomax,2)
      double precision  Cgr ! ,gamma
      common /als/ als(nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
c      common /factorial/ gamma(maxfac)
      common /factratio/ Cgr(0:lomax,komax)
      if (nodeid.eq.1)
     >write(4,'("Set s.p. basis")')
      do nsp=1,nspm
         ns = nset(nsp)
         l = lo(nsp)
         k = ko(nsp)
         als(nsp) = al(l,ns)
c         M(nsp) = dsqrt(als(nsp)*gamma(k)/gamma(2*l+k+2))
         M(nsp) = dsqrt(als(nsp)/Cgr(l,k))
      end do
      return
      end
c-----------------------------------------------------------
      subroutine funlaguer(M,fl,maxf,minf,gl,maxg,ming,f8,kom)
c     this routine forms Sturmian functions (fl(i,nsp)) and second
c     derivatives (gl(i,nsp)) for each of them  from Laguerre polinomials.
      include 'par.f'
      common /meshrr/ nr,gridr(nmaxr,3)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      double precision  als, M(nspmCI), nr1, grid, x
      common /als/ als(nspmax)
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      real gl(nr,nspm)
      integer  maxg(nspmCI),ming(nspmCI)
      common /gridr/ nr1, grid(nmaxr)
      double precision pl(komax), f8(nr,kom)            

      do n=1,nspm
         l = lo(n)
         k = ko(n)
         il = 2*l+2
         call  lagpol8(il,als(n),f8,nr,k,grid,nr,pl)
         do i=1,nr
            x = (gridr(i,1))*als(n)
            fl(i,n) = M(n)*(x**(l+1))*exp(-x/2.0D0)*f8(i,k)
c            gl(i,n) = 0.0
            gl(i,n) = ((l*(l+1))/(x*x) - (k+l)/x +
     >         0.25d0)*(als(n)**2) * fl(i,n)
         end do
         call minmaxi(fl(1,n),nr,i1,i2)
         minf(n)=i1
         maxf(n)=i2
c
         if(k.gt.1) then
            call  lagpol8(il+1,als(n),f8,nr,k-1,grid,nr,pl)
            do i=1,nr
               x = DBLE(gridr(i,1))*als(n)
               gl(i,n) = dble(gl(i,n)) +
     >            (als(n)**2)*M(n)*(x**l)*dexp(-x/2.0D0)*f8(i,k-1)
            end do
         end if
         call minmaxi(gl(1,n),nr,i1,i2)
         ming(n)=i1
         maxg(n)=i2
      end do
      return
      end
c-----------------------------------------------------------
      subroutine minmaxi(f,nr,i1,i2)
      include 'par.f'
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      dimension  f(nmaxr)
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
      subroutine checkortsq(fl,maxf,minf)
      use vmat_module, only: nodeid
      include 'par.f'
      double precision  ortint
      common /ortog/  ortint(nspmax,nspmax)
      double precision  rnorm, rnorm1, r1elk
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
c      
      if (nodeid.eq.1)
     >write(4,'("check orthogonality of s.p. states")') 
      if (nodeid.eq.1)
     >write(20,'("check orthogonality of s.p. states")') 
      do n1=1,nspm
         do n2=n1,nspm
            rnorm = 0.0D0
            rnorm1 = 0.0D0
            if(lo(n1).eq.lo(n2)) then
               rnorm1 =  r1elk(0,fl(1,n1),fl(1,n2),
     >            minf(n1),minf(n2),maxf(n1),maxf(n2))
               if(n1.eq.n2) rnorm = 1.0D0
               if(nset(n1).ne.nset(n2)) then
                  rnorm =  r1elk(0,fl(1,n1),fl(1,n2),
     >               minf(n1),minf(n2),maxf(n1),maxf(n2))
               end if
            end if
c$$$            write(20,'("l1,l2,k1,k2,n1,n2 =",6I5,"   <k1l|k2l> =",
c$$$     >         F10.5,E15.5)')lo(n1),lo(n2),ko(n1),ko(n2),n1,n2,
c$$$     >         real(rnorm),rnorm1
            ortint(n1,n2) = rnorm
            ortint(n2,n1) = rnorm
         end do 
      end do 
      return
      end      
c-----------------------------------------------------------------
      subroutine e1mat(fl,maxf,minf,gl,maxg,ming,nr,nspm)
      use vmat_module, only: nodeid
      include 'par.f'
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      real gl(nr,nspm)
      integer  maxg(nspmCI),ming(nspmCI)
      double precision  re1,oneelint,e1r
      common/hame1/ e1r(nspmax,nspmax)
      if (nodeid.eq.1)
     >write(4,'("enter e1mat")')
      do n1=1,nspm
         do n1p=1,n1
            re1 = oneelint(n1,n1p,fl,maxf,minf,gl,maxg,ming) 
            e1r(n1,n1p) = re1
            e1r(n1p,n1) = re1
         end do
      end do
      return
      end
c-----------------------------------------------------------------
      double precision function oneelint(n1,n2,fl,maxf,minf,
     >     gl,maxg,ming)
      include 'par.f'
      common /meshrr/ nr, gridr(nmaxr,3)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      real  gl(nr,nspm)
      integer  maxg(nspmCI),ming(nspmCI)
      double precision  Z, sum, r
      common /Zatom/ Z
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      oneelint = 0.0D0
      if(lo(n1).ne.lo(n2)) return
      l = lo(n1)
      sum = 0d0

      maxrr=min(maxf(n1),maxg(n2))
      minrr=max(minf(n1),ming(n2))      
      do i=minrr,maxrr
         sum = sum + dble(gridr(i,3) * fl(i,n1) * gl(i,n2))
      end do

      maxrr=max(maxf(n1),maxf(n2))
      minrr=min(minf(n1),minf(n2))
      do i=minrr,maxrr
         r = gridr(i,1)
         sum = sum + dble(gridr(i,3) * fl(i,n1) * 
     >      fl(i,n2))*(2d0*Z*rpow2(i,0) - dble(l*(l+1.))/(r*r))
!     >      fl(i,n2))*(2d0*Z/r - dble(l*(l+1.))/(r*r))
      end do

      oneelint = - sum/2d0
      return
      end      
c************************************************************************

c$$$
c$$$      maxrr=max(maxf(n1),maxf(n2),maxg(n2))
c$$$      minrr=min(minf(n1),minf(n2),ming(n2))
c$$$      do i=minrr,maxrr
c$$$         r = gridr(i,1)
c$$$         sum = sum + gridr(i,3) * fl(i,n1) * ( gl(i,n2) +
c$$$     >      fl(i,n2)*(2.0*rZ/r - l*(l+1.)/(r*r)) )
c$$$      end do
c$$$

c-----------------------------------------------------------------
c** main - configuration interaction program for H **********
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
      subroutine maindiagH
      include 'par.f'
      character atom*8
      dimension fl(maxr,ncmax),maxf(ncmax),minf(ncmax)
      double precision  als(ncmax),al(2), E(ncmax)
      common /Zatom/ Z
      dimension Vpot(maxr)
c
      open(5,file='F5')
      open(6,file='F6')
      open(10,file='states')
      read(5,*) atom
      write(6,*) atom
      read(5,*) Z
      write(6,'("Z =",F10.5)') Z
      read(5,*) l
      write(6,'("l =",I3)') l
      read(5,*) nk
      write(6,'("number of H state nk=",10I3)') nk
c     nset=1
      read(5,*) ral
      write(6,'("al(set=1) =",F10.5)') ral
c     nset=2
      read(5,*) ndif
      write(6,'("ndif =",I5)') ndif
      read(5,*) raldif
      write(6,'("al(set=2) =",F10.5)') raldif
      read(5,*) nkd
      write(6,'("nkd =",I5)') nkd
      al(1) = ral
      al(2) = raldif
      do k=1,nk
         als(k) = al(1)
         if(ndif.ne.0.and.k.le.nkd) then
            als(k) = al(2)
         end if
      end do
c     here we call hydrogen CI routine
      call hyd1el(z,l,nk,als,E,Vpot,fl,minf,maxf)
      stop
      end
c**********************************************************************
      subroutine hyd1el(z,l,nk,als,E,Vpot,fl,minf,maxf)
      include 'par.f'
      double precision  CI(ncmax,ncmax)
      dimension fl(maxr,ncmax),maxf(ncmax),minf(ncmax)
      dimension gl(maxr,ncmax),maxg(ncmax),ming(ncmax)
      double precision  ortint(ncmax,ncmax)
      dimension Vpot(maxr),E(ncmax),als(ncmax)
c     here we set all one-electron quantities 
      call set1el(l,nk,als,fl,minf,maxf,gl,ming,maxg,ortint)
c     here H cofiguration interaction program is called
c     this programm is in the file  H-3.f

      call hyd(z,l,nk,CI,E,Vpot,fl,minf,maxf,gl,ming,maxg,ortint)
c-------------------------------------------------------------------
c$$$      do N=1,nk
c$$$         write(10,'("N=",I3)') N
c$$$         do n1=1,nk
c$$$            if(CI(n1,N).ne.0.0D0) then
c$$$               write(10,*) n1,real(CI(n1,N))
c$$$            end if
c$$$         end do
c$$$      end do
c$$$      call ort1el(nk,CI,ortint)
c-------------------------------------------------------------------
      call rear1el(nk,fl,maxf,minf,CI,ortint)
c-------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------
c****************   hyd(l,nk,CI,w)  ******************************
c-----------------------------------------------------------------
c     subroutine hyd(...)   <--   call eigv(ncm,H,bb)
c     - all one electron calculations in the file setmax.f
c     - subroutine eigv1el(nk,H,bb,CI,w) 
c     - subroutine rsg(nm,n,a,b,E,matz,z,fv1,fv2,ierr) -from netlib
c     double precision function el1int(...)
      subroutine hyd(z,l,nk,CI,E,Vpot,fl,minf,maxf,
     >   gl,ming,maxg,ortint)
      include 'par.f'
      double precision H(ncmax,ncmax), bb(ncmax,ncmax), 
     >  CI(ncmax,ncmax), w(ncmax)
      dimension fl(maxr,ncmax),maxf(ncmax),minf(ncmax)
      dimension gl(maxr,ncmax),maxg(ncmax),ming(ncmax)
      double precision  ortint(ncmax,ncmax)
      dimension Vpot(maxr),E(ncmax)

      do nc=1,nk
         do ncp=nc,nk
            H(nc,ncp) = dble(el1int(z,l,nc,ncp,Vpot,fl,minf,maxf,
     >         gl,ming,maxg))
            H(ncp,nc) = H(nc,ncp)
            bb(nc,ncp) =  ortint(nc,ncp)
            bb(ncp,nc) = bb(nc,ncp)
         end do
      end do
c     start diagonalization
      call eigv1el(nk,H,bb,CI,w)
      do n = 1, nk
         E(n) = sngl(w(n)) * 2.0
      enddo 
c     finish diagonalization
      return
      end
c-----------------------------------------------------------------
c*************** Find the eigenvalues and eigenvectors  **********
c********************* of the Hamiltonian matrix *****************
c-----------------------------------------------------------------
      subroutine eigv1el(n,H,bb,CI,w)
      include 'par.f'
      double precision sum, dnorm(ncmax,ncmax)
      double precision H(ncmax,ncmax),bb(ncmax,ncmax),w(ncmax),
     >   CI(ncmax,ncmax),fv1(ncmax),fv2(ncmax),b(ncmax,ncmax)
      matz=2
      do i=1,n
         do j=1,n
            b(j,i) = bb(j,i)
         end do
      end do
      call  rsg(ncmax,n,H,b,w,matz,CI,fv1,fv2,ierr)
      if (ierr.ne.0) write(6,'("WARNING ierr =",I3)') ierr
c$$$      write(6,'("eigenvalues in a.u.")')
c$$$      write(6,'(5F12.6)') (real(w(i)), i=1,n)
c$$$      write(6,'("eigenvectors  :  H-states")')
c$$$      do i=1,n
c$$$c         write(6,'(7F10.5)') (real(CI(i,j)), j=1,n)
c$$$      end do
c$$$      write(6,'("overlap matrix : overlap.result")')
c$$$      do i=1,n
c$$$c         write(6,'(7F10.5)') (real(bb(i,j)), j=1,n)
c$$$      end do
      return
      end
c************************************************************************
      subroutine rear1el(nk,fl,maxf,minf,CI,ortint)
      include 'par.f'
      double precision  ortint(ncmax,ncmax)
      double precision  CI(ncmax,ncmax)
      common /meshrr/ nr, gridr(maxr,3)
      dimension  p(maxr,ncmax),maxp(ncmax),minp(ncmax)
      dimension fl(maxr,ncmax),maxf(ncmax),minf(ncmax)
      double precision  rnorm,r1elk,sum
      
      do N=1,nk
         do i=1,nr
            sum = 0.0D0
            do nsp=1,nk
               sum = sum + CI(nsp,N) * fl(i,nsp)
            end do
            p(i,N) = sum
         end do
         call minmaxi(p(1,N),nr,i1,i2)
         maxp(N) = i2
         minp(N) = i1
      end do
      do n1=1,nk
         do n2=n1,nk
            rnorm = r1elk(0,p(1,n1),p(1,n2),
     >         minp(n1),minp(n2),maxp(n1),maxp(n2))
            ortint(n1,n2) = rnorm
            ortint(n2,n1) = ortint(n1,n2) 
         end do
      end do
      do n=1,nk
         do i=1,nr
            fl(i,n) = p(i,n)
         end do
         maxf(n) = maxp(n)
         minf(n) = minp(n)
      end do
      return
      end
c-----------------------------------------------------------------
      function el1int(Z,l,n1,n2,Vpot,
     >   fl,minf,maxf,gl,ming,maxg)
      include 'par.f'
      common /meshrr/ nr, gridr(maxr,3)
      dimension fl(maxr,ncmax),maxf(ncmax),minf(ncmax)
      dimension gl(maxr,ncmax),maxg(ncmax),ming(ncmax)
      dimension Vpot(maxr),f1(maxr),f2(maxr)
      sum = 0.0
      maxrr=max(maxf(n1),maxf(n2),maxg(n2))
      minrr=min(minf(n1),minf(n2),ming(n2))
C  Note that VPOT is in Rydbergs
      do i=minrr,maxrr
         r = gridr(i,1)
         sum = sum + gridr(i,3) * fl(i,n1) * ( gl(i,n2) +
     >      fl(i,n2)*(2.0*(Z/r) - l*(l+1.)/(r*r) - Vpot(i)) )
      end do
      do i = 1, nr
         f1(i) = fl(i,n1) * gridr(i,3)
         f2(i) = fl(i,n2) * gridr(i,3)
      enddo 
      call fcexch(f1,nr,f2,res,l)
      sum = sum - res * 2.0
      el1int = - dble(sum/2.0)
      return
      end      
c**************************************************************************
c     subroutune  set1el - set overall max s.p. state
c     call funlaguer1el -   forms Sturmian functians from Laguerre polinomials
c     subroutine checkort1el
c     subroutine  ort1el(nk,CI)
c---  --------------------------------------------------------------------
c     this subroutine is called to make all necessary one-electron quantities:
c     set sq.int.basis (call  factorials,call funlaguer1el),
c     overlap integrals between s.p. states (call checkortsq),
c     one electron integrals (call e1mat).
c     note - called only ones in calculation of HE states of all symmetries
      subroutine  set1el(l,nk,als,fl,minf,maxf,
     >   gl,ming,maxg,ortint)
      include 'par.f'
      double precision  M(ncmax)
      dimension fl(maxr,ncmax),maxf(ncmax),minf(ncmax),als(ncmax)
      dimension gl(maxr,ncmax),maxg(ncmax),ming(ncmax)
      double precision  ortint(ncmax,ncmax)
      double precision  gamman(ncmax + lmax * 2 + 2)
      call factor(nk,l,gamman)
      call funlaguer1el(l,nk,M,als,fl,minf,maxf,gl,ming,maxg,gamman)
      call checkort1el(nk,als,fl,minf,maxf,ortint)
      return
      end
c-----------------------------------------------------------------
c**************       Set s.p. basis              ****************
c-----------------------------------------------------------------
      subroutine funlaguer1el(l,nk,M,als,fl,minf,maxf,
     >   gl,ming,maxg,gamman)
c     this routine forms Sturmian functions (fl(i,nsp)) and second
c     derivatives (gl(i,nsp)) for each of them  from Laguerre polinomials.
      include 'par.f'
      common /meshrr/ nr,gridr(maxr,3)
      double precision  M(ncmax),grid(maxr),x
      dimension fl(maxr,ncmax),maxf(ncmax),minf(ncmax),als(ncmax)
      dimension gl(maxr,ncmax),maxg(ncmax),ming(ncmax)
      double precision pl(ncmax), f8(maxr,ncmax), dals
      double precision  gamman(nk+2+2*l)

      do i = 1, nr
         grid(i) = dble(gridr(i,1))
      enddo 
      il = 2*l+2
      do n=1,nk
         dals = dble(als(n))
         M(n) = dsqrt(dals*gamman(n)/gamman(2*l+n+2))
         call  lagpol8(il,dals,f8,maxr,n,grid,nr,pl)
         do i=1,nr
            x = DBLE(gridr(i,1)) * dals
            fl(i,n) = M(n)*(x**(l+1))*dexp(-x/2.0D0)*f8(i,n)
            gl(i,n) = (dble(l*(l+1))/(x*x) - dble(n+l)/x +
     >         0.25d0)*(dals**2) * fl(i,n)
         end do
         call minmaxi(fl(1,n),nr,i1,i2)
         minf(n)=i1
         maxf(n)=i2
c
         if(n.gt.1) then
            call  lagpol8(il+1,dals,f8,maxr,n-1,grid,nr,pl)
            do i=1,nr
               x = DBLE(gridr(i,1))*dals
               gl(i,n) = dble(gl(i,n)) +
     >            (dals**2)*M(n)*(x**l)*dexp(-x/2.0D0)*f8(i,n-1)
            end do
         end if
         call minmaxi(gl(1,n),nr,i1,i2)
         ming(n)=i1
         maxg(n)=i2
      end do
      return
      end
c-----------------------------------------------------------
      subroutine checkort1el(nk,als,fl,minf,maxf,ortint)
      include 'par.f'
      double precision  ortint(ncmax,ncmax)
      double precision  rnorm, r1elk
      dimension  fl(maxr,ncmax),maxf(ncmax),minf(ncmax), als(ncmax)
c      
      do n1=1,nk
         do n2=n1,nk
            rnorm = 0.0D0
            if(n1.eq.n2) rnorm = 1.0D0
            if(als(n1).ne.als(n2)) then
               rnorm =  r1elk(0,fl(1,n1),fl(1,n2),
     >            minf(n1),minf(n2),maxf(n1),maxf(n2))
            end if
            ortint(n1,n2) = rnorm
            ortint(n2,n1) = rnorm
         end do 
      end do 
      return
      end      
c*********************************************************************
      subroutine  ort1el(nk,CI,ortint)
      include 'par.f'
      double precision  ortint(ncmax,ncmax)
      double precision  CI(ncmax,ncmax)
      double precision sum, dr
      do N=1,nk
         do Np=1,N
            sum = 0.0D0
            do n1=1,nk
               do n1p=1,nk
                  dr = ortint(n1,n1p)
                  if(dr.ne.0.0D0) then    
                     sum = sum + dr*CI(n1,N)*CI(n1p,Np)
                  end if
               end do
            end do
         end do
      end do
      return
      end
c-------------------------------------------------------------------
      subroutine factor(nk,l,gamman)
c     n! = gamman(n+1)
      include 'par.f'
      double precision  gamman(nk+2+2*l)

      gamman(1) = 1.D0
      do i=2,nk+2+2*l
         ii=i-1
         gamman(i) = DBLE(ii)*gamman(ii)
      end do
      return
      end

      

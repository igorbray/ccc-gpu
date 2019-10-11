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
      double precision  als(nspmax), M(nspmax)
      dimension fl(nmaxr,nspmax),maxf(nspmax),minf(nspmax)
      dimension gl(nmaxr,nspmax),maxg(nspmax),ming(nspmax)
      double precision  ortint(nspmax,nspmax)
      double precision  gamma(maxfac)
      open(20,file='overlap.result')
      call setgrids
      call  factor(nk,l,gamma)
      call funlaguer1el(l,nk,M,als,fl,minf,maxf,gl,ming,maxg,gamma)
      call checkort1el(nk,als,fl,minf,maxf,ortint)
      return
      end
c-----------------------------------------------------------------
c**************       Set s.p. basis              ****************
c-----------------------------------------------------------------
      subroutine funlaguer1el(l,nk,M,als,fl,minf,maxf,
     >   gl,ming,maxg,gamma)
c     this routine forms Sturmian functions (fl(i,nsp)) and second
c     derivatives (gl(i,nsp)) for each of them  from Laguerre polinomials.
      include 'par.f'
      common /gridrr/ nr,gridr(nmaxr,3)
      double precision  M(nspmax),nr1,grid,x,als(nspmax)
      dimension fl(nmaxr,nspmax),maxf(nspmax),minf(nspmax)
      dimension gl(nmaxr,nspmax),maxg(nspmax),ming(nspmax)
      common /gridr/ nr1, grid(nmaxr)
      double precision pl(nspmax), f8(nmaxr,nspmax)
      double precision  gamma(maxfac)
      write(6,'("Set s.p. basis")')
      il = 2*l+2
      do n=1,nk
         M(n) = dsqrt(als(n)*gamma(n)/gamma(2*l+n+2))
         call  lagpol8(il,als(n),f8,nmaxr,n,grid,nr,pl)
         do i=1,nr
            x = DBLE(gridr(i,1))*als(n)
            fl(i,n) = M(n)*(x**(l+1))*dexp(-x/2.0D0)*f8(i,n)
            gl(i,n) = (dble(l*(l+1))/(x*x) - dble(n+l)/x +
     >         0.25d0)*(als(n)**2) * fl(i,n)
         end do
         call minmaxi(fl(1,n),nr,i1,i2)
         minf(n)=i1
         maxf(n)=i2
c
         if(n.gt.1) then
            call  lagpol8(il+1,als(n),f8,nmaxr,n-1,grid,nr,pl)
            do i=1,nr
               x = DBLE(gridr(i,1))*als(n)
               gl(i,n) = dble(gl(i,n)) +
     >            (als(n)**2)*M(n)*(x**l)*dexp(-x/2.0D0)*f8(i,n-1)
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
      double precision  ortint(nspmax,nspmax)
      double precision  rnorm, r1elk, als(nspmax)
      dimension  fl(nmaxr,nspmax),maxf(nspmax),minf(nspmax)
c      
      write(6,'("check ortoganality of s.p. states")') 
      write(20,'("check ortogonality of s.p. states")') 
      do n1=1,nk
         do n2=n1,nk
            rnorm = 0.0D0
            if(n1.eq.n2) rnorm = 1.0D0
            if(als(n1).ne.als(n2)) then
               rnorm =  r1elk(0,fl(1,n1),fl(1,n2),
     >            minf(n1),minf(n2),maxf(n1),maxf(n2))
            end if
            write(20,'("n1,n2 =",2I5)') n1,n2
            write(20,'("<n1|n2> =",F10.5)') real(rnorm)
            ortint(n1,n2) = rnorm
            ortint(n2,n1) = rnorm
         end do 
      end do 
      return
      end      
c*********************************************************************
      subroutine  ort1el(nk,CI,ortint)
      include 'par.f'
      double precision  ortint(nspmax,nspmax)
      double precision  CI(nspmax,nspmax)
      double precision sum, dr
      write(10,'("<N|Np>")') 
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
            write(10,'("N,Np=",2I3,",  <N|Np>=",F10.5)')
     >         N, Np, real(sum)
         end do
      end do
      return
      end
c-------------------------------------------------------------------
      subroutine factor(nk,l,gamma)
c     n! = gamma(n+1)
      include 'par.f'
      double precision  gamma(maxfac)
      write(6,'("set gamma")')
      gamma(1) = 1.D0
      do i=2,nk+2+2*l
         ii=i-1
         gamma(i) = DBLE(ii)*gamma(ii)
      end do
      return
      end

      

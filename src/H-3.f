c-----------------------------------------------------------------
c****************   hyd(l,nk,CI,w)  ******************************
c-----------------------------------------------------------------
c     subroutine hyd(...)   <--   call eigv(ncm,H,bb)
c     - all one electron calculations in the file setmax.f
c     - subroutine eigv1el(nk,H,bb,CI,w) 
c     - subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr) -from netlib
c     double precision function el1int(...)
      subroutine hyd(l,nk,CI,w,Vpot,fl,minf,maxf,
     >   gl,ming,maxg,ortint)
      include 'par.f'
      double precision H(nspmax,nspmax), bb(nspmax,nspmax), 
     >  CI(nspmax,nspmax), w(nspmax),el1int
      dimension fl(nmaxr,nspmax),maxf(nspmax),minf(nspmax)
      dimension gl(nmaxr,nspmax),maxg(nspmax),ming(nspmax)
      double precision  ortint(nspmax,nspmax)
      dimension Vpot(nmaxr)
      write(6,'("enter hyd")')
      do nc=1,nk
         do ncp=nc,nk
            H(nc,ncp) = el1int(l,nc,ncp,Vpot,fl,minf,maxf,
     >         gl,ming,maxg)
            H(ncp,nc) = H(nc,ncp)
            bb(nc,ncp) =  ortint(nc,ncp)
            bb(ncp,nc) = bb(nc,ncp)
         end do
      end do
c     start diagonalization
      call eigv1el(nk,H,bb,CI,w)
c     finish diagonalization
      return
      end
c-----------------------------------------------------------------
c*************** Find the eigenvalues and eigenvectors  **********
c********************* of the Hamiltonian matrix *****************
c-----------------------------------------------------------------
      subroutine eigv1el(n,H,bb,CI,w)
      include 'par.f'
      double precision sum, dnorm(nspmax,nspmax)
      double precision H(nspmax,nspmax),bb(nspmax,nspmax),w(nspmax),
     >   CI(nspmax,nspmax),fv1(nspmax),fv2(nspmax),b(nspmax,nspmax)
      double precision CI(nspmax,nspmax),w(nspmax)
      write(6,'("start diagonalization")')
      write(6,'("overlap matrix")')
      matz=2
      do i=1,n
         do j=1,n
            b(j,i) = bb(j,i)
         end do
      end do
      call  rsg(nspmax,n,H,b,w,matz,CI,fv1,fv2,ierr)
      write(6,'("ierr =",I3)') ierr
      write(6,'("eigenvalues in a.u.")')
      write(6,'(5F12.6)') (real(w(i)), i=1,n)
      write(6,'("eigenvectors  :  H-states")')
      do i=1,n
c         write(6,'(7F10.5)') (real(CI(i,j)), j=1,n)
      end do
      write(6,'("overlap matrix : overlap.result")')
      do i=1,n
c         write(6,'(7F10.5)') (real(bb(i,j)), j=1,n)
      end do
      do j=1,n
         do jp=1,n
            sum = 0.0D0
            do i=1,n
               do ip=1,n
                  sum = sum + CI(i,j)*CI(ip,jp)*bb(i,ip)
               end do
            end do
            dnorm(j,jp) = sum
         end do
      end do
      write(6,'("overlap matrix of H eigenstate :  H-states")')
      do i=1,n
c         write(6,'(7F10.5)') (real(dnorm(i,j)), j=1,n)
      end do
      return
      end
c************************************************************************
      subroutine rear1el(nk,fl,maxf,minf,CI,ortint)
      include 'par.f'
      double precision  ortint(nspmax,nspmax)
      double precision  CI(nspmax,nspmax)
      common /gridrr/ nr, gridr(nmaxr,3)
      dimension  p(nmaxr,nspmax),maxp(nspmax),minp(nspmax)
      dimension fl(nmaxr,nspmax),maxf(nspmax),minf(nspmax)
      double precision  rnorm,r1elk,sum
      
      write(10,'("******************************************")') 
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
      write(10,'("rearranged s.p. state overlaps")') 
      write(10,'("it is hydrogen eigenstates")') 
      do n1=1,nk
         do n2=n1,nk
            rnorm = r1elk(0,p(1,n1),p(1,n2),
     >         minp(n1),minp(n2),maxp(n1),maxp(n2))
            write(10,'("n1,n2 =",2I5)') n1,n2
            write(10,'("<n1|n2> =",F10.5)') real(rnorm)
            if(dabs(rnorm).le.1.0D-05) then
               rnorm = 0.0D0
               write(10,'(" <n1|n2> <= 1.0D-05 ")')
            end if
            if(dabs(rnorm - 1.0D0).le.1.0D-05) then
               rnorm = 1.0D0
               write(10,'(" <n1|n2> - 1.0D0 <= 1.0D-05 ")')
            end if
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
         write(10,'("n,minf(n),maxf(n)=",3I5)') n,minf(n),maxf(n)
      end do
      return
      end
c-----------------------------------------------------------------
      double precision function el1int(l,n1,n2,Vpot,
     >   fl,minf,maxf,gl,ming,maxg)
      include 'par.f'
      common /gridrr/ nr, gridr(nmaxr,3)
      dimension fl(nmaxr,nspmax),maxf(nspmax),minf(nspmax)
      dimension gl(nmaxr,nspmax),maxg(nspmax),ming(nspmax)
      common /Zatom/ Z
      dimension Vpot(nmaxr)
      sum = 0.0
      maxrr=max(maxf(n1),maxf(n2),maxg(n2))
      minrr=min(minf(n1),minf(n2),ming(n2))
      do i=minrr,maxrr
         r = gridr(i,1)
         sum = sum + gridr(i,3) * fl(i,n1) * ( gl(i,n2) +
     >      fl(i,n2)*(2.0*(Z/r + Vpot(i)) - l*(l+1.)/(r*r)) )
      end do
      el1int = - dble(sum/2.0)
      return
      end      

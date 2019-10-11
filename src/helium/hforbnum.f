c      subroutine hforb
c      subroutine hforbnum
c
c     this routine reads S.T.O representation of h.f. orbitals from
c     file F5hforb and form these orbitals.
c     Note: routines configsp and spcoef better be switched off if hforb works.
c     see file setmax.f where this routine is called.
      subroutine hforb
      include 'par.f'
      dimension  Nhf(0:lomax), NML(0:lomax), n(20), g(20), Chf(20)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      double precision   Cpol, gv
      common /polcoeffd/  Cpol(komax,nspmCI), gv(komax,nspmCI)
      common /polcoeffs/  im(nspmCI), nv(komax,nspmCI)
      double precision  gamma
      common /factorial/ gamma(maxfac)
c     
      open(55,file='F5hforb') 
      open(66,file='F6hforb') 
      nspc=0
      read(55,*) Lhfm
      write (66,'("Lhfm=",I5)') Lhfm
      do l=0,Lhfm
         read(55,*) Nhf(l)
         write(66,'("Nhf(l) =",I5)') Nhf(l)
         read(55,*) NML(l)
         write(66,'("NML(l) =",I5)') NML(l)
         do i=1,NML(l)
            read(55,*) rn, lv, rg
            n(i) = rn
            g(i) = rg
            write (66,*) n(i),lv,g(i)
         end do
         do nk=1,Nhf(l)
            nspc = nspc +1
            read(55,*) nsp
            write(66,*) nsp, l
            ko(nsp) = nk
            lo(nsp) = l
            nset(nsp) = nsp
            im(nsp) = NML(l)
            read(55,*) (Chf(j), j=1,NML(l))
            write(66,'(5F12.8)') (Chf(j), j=1,NML(l))
            do j=1,NML(l)
               nv(j,nsp) = n(j)
               gv(j,nsp) = dble(g(j))
               Cpol(j,nsp) =  Chf(j)*
     >            dsqrt( (2.0d0*gv(j,nsp))**(2*nv(j,nsp)+1)/
     >            gamma(2*nv(j,nsp)+1) ) 
            end do
         end do
      end do
      nspm = nspc
      close(55)
      close(66)
      return
      end
c-----------------------------------------------------------------------
c     this routine  use  S.T.O representation of h.f. orbitals to form
c     numerical h.f. orbitals (fl(i,nsp)) and second derivatives (gl(i,nsp))
c     for each of them.
      subroutine hforbnum(gl)
      include 'par.f'
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      double precision   Cpol, gv
      common /polcoeffd/  Cpol(komax,nspmCI), gv(komax,nspmCI)
      common /polcoeffs/  im(nspmCI), nv(komax,nspmCI)
      common /meshrr/ nr,gridr(nmaxr,3)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      real  gl(nr,nspm)
      common/minmaxg/  maxg(nspmCI),ming(nspmCI)
      double precision sum, sum2, r, g
c     
      do nsp=1,nspm
         do i=1,nr
            sum = 0.0d0
            sum2 = 0.0d0
            do mi=1,im(nsp)
               r = dble(gridr(i,1))
               n = nv(mi,nsp)
               g = gv(mi,nsp)
               sum = sum + Cpol(mi,nsp) * r**n * dexp(-g*r)
               sum2 = sum2 + Cpol(mi,nsp)*dexp(-g*r)*r**(n-2)*
     >            (dble(n*(n-1)) - dble(2.*n)*g*r + g**2*r**2)
            end do
            fl(i,nsp) = sum
            gl(i,nsp) = sum2
         end do
         call minmaxi(fl(1,nsp),nr,i1,i2)
         minf(nsp)=i1
         maxf(nsp)=i2
         call minmaxi(gl(1,nsp),nr,i1,i2)
         ming(nsp)=i1
         maxg(nsp)=i2
      end do
      return
      end


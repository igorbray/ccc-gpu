c======================================================================
c  The following routine diagonalizes the Hamiltonian
C             2
C            d         l (l + 1)      Zinf
C     H = - -----   +  ---------  -  ------  + pot0(r)
c               2          2            r
c            d r          r
C
C  in a Laguerre basis of size NPS and exponential fall off factor ALF.
C
C                          2l+2           l+1
C   Laguerre basis:       L   (2 alf r)  r    exp (- alf r)
C                          n-1
C
C  To get the first bound state for each l exactly set alf to be
C  Z of atom divided by (l+1) and take NPS sufficiently large.
C  This program has been tried with NPS up to 100.
C===================================================================
C                         i   i    i   i    i      i    o    o       o
      subroutine makeps(Zinf,exch,alf, l, expcut, nps, ps, energy, jmin, 
     >   jmax, gridr, pot0, maxpot, nr, cknd)
C         o      i     i       i    i
C  EXPCUT is used to avoid underflows at the start and end of the wave
C         functions PS.
C  ENERGY is the eigenvalues array in Rydbergs
C  GRIDR  is the grid containing NR R points and integration weights.
C  POT0(i) is defined for i = 1, maxpot is the input potential in Rydbergs
      include 'par.f'
c$$$      parameter (ndim=100,nld=3000)
      parameter (ndim=ncmax,nld=maxr)
      implicit real*8 (a-h, o-z)
      logical exch, exists
c
      real ps(maxr,ndim), gridr(maxr,3), pot0(maxr),energy(ndim),
     >   expcut, Zinf, alf
      real*8, allocatable :: rierm(:), wfierm(:,:)
c
      dimension c(ndim, ndim), h(ndim, ndim),ch(ndim,ndim)
      dimension enrgr(ndim), dnorm(ndim), work(ndim),
     >   dpot(nld), cknd(ndim,ndim)
      dimension grid(nld,3), fac(0:1000)
      dimension f(nld, ndim), f1(nld, ndim), jmin(ndim), jmax(ndim)
c
c$$$      real deltar
c$$$      common /delta/ deltar(maxr,maxr), ndelta
c$$$      common /pspace/ nabot,labot,natop,latop,ntype,ipar,nze,ninc,linc,
c$$$     >   lactop
c$$$      dimension nabot(0:lamax), natop(0:lamax)

      dimension potmtrx(ndim, ndim), pl(ndim), weight(nld)
      dimension dnorm2(ndim), work2(8*ndim), iwork(5*ndim)
      character ch1
      ch1(i)=char(i+ichar('0'))

c$$$      real*16 test
C
C     Parameter checking
C
      if (nr.gt.nld) then
         print*,'Increase NLD to at least:',nr
         stop 'Increase NLD'
      end if 
      if (nps.gt.ndim) then
         print*,'Increase NDIM to at least:',nps
         stop 'Increase NDIM'
      end if 
c
c     Define real*8 input arrays
c

      r0 = 0d0
      i = 1
      m = 0
      dr = gridr(i,1)
      do while (i.le.nr)
         grid(i,1) = r0 + dfloat(i-m) * dr
         grid(i,3) = dfloat(2 * mod(i,2) + 2) * dr / 3d0
         if (abs((grid(i,1)-gridr(i,1))/grid(i,1)).gt.1e-5) then
c$$$            print*,'Doubling at i:',i
            grid(i,1) = grid(i,1) + dr
            r0 = grid(i,1)
            dr = 2d0 * dr
            grid(i-1,3) = grid(i-1,3) * 1.5d0
            grid(i,3) = dfloat(2 * mod(i,2) + 2) * dr / 3d0
            m = i
         endif
         i = i + 1
      enddo
      grid(nr,3) = grid(nr,3) / 2d0

            
      inquire(file='iermwf'//ch1(l),exist=exists)
      if (exists) then
         print*,'Running with the IERM basis'
         open(42,file='iermwf'//ch1(l))
         read(42,*) nierm,lierm,npsierm
         if (lierm.ne.l) stop 'L and LIERM are not the same'
         if (npsierm-lierm.ne.nps)
     >      stop 'NPS and NPSIERM are not the same'
         read(42,*) r0ierm,(energy(n),n=1,npsierm-lierm)
         allocate(rierm(nierm))
         allocate(wfierm(nierm,npsierm-lierm))
         do i = 1, nierm
            read(42,*) rierm(i),(wfierm(i,n),n=1,npsierm-lierm)
         enddo
         close(42)
         i = 1
         do while (gridr(i,1).lt.r0ierm.and.i.lt.nr)
            i = i + 1
         enddo
         imax = i
         ps(:,:) = 0.0
         do n = 1, npsierm-lierm
            call intrpl(nierm,rierm,wfierm(1,n),imax,grid(1,1),f(1,n))
            do i = 1, imax-1
               ps(i,n) = real(f(i,n))
            enddo
            jmin(n)=1
            jmax(n)=imax
         enddo
         deallocate(wfierm)
         deallocate(rierm)
         return
      endif 
c$$$      do n=1,3,2
c$$$         do i=1,nr
c$$$            if (abs((grid(i,1)-gridr(i,1))/grid(i,1)).gt.1e-5.or.
c$$$     >         abs((grid(i,3)-gridr(i,3))/grid(i,3)).gt.1e-5) then
c$$$               print*,'Problem with grid for I, N:',i,n
c$$$            endif 
c$$$         end do
c$$$      end do
      test = 0d0
      do i=1,nr
c$$$         grid(i,1) = dfloat(i) * gridr(nr,1) / dfloat(nr)
c$$$         grid(i,3) = dfloat( 2 * mod(i,2) + 2 ) * grid(1,1) / 3d0
C  DPOT is in a.u.
         dpot(i) = dble(pot0(i)) / 2.0
         weight(i) = grid(i,3)
c$$$         psi1sr = 2d0 * grid(i,1) * exp ( - qreal(grid(i,1)))
c$$$         test = test +  psi1sr * psi1sr * weight(i)
      end do
c$$$      print*,'Analytic test:', abs(test - 1d0),test
      dalf = dble(alf)
      dz = dble(Zinf)
      dexpcut = dble(expcut)
c
c     Define Matrix of the potential
c
      call potl (exch, dpot, maxpot, L, dalf, potmtrx, nps, ndim, 
     >   f, pl, grid, weight, nr, nld)
c
      nfac = 2 * nps + 2 * L
C  Diagonalize the Hamiltonian
      call basisl (potmtrx, dz, L, dalf, nps, ndim, enrgr, c, h,
     >   dnorm, dnorm2, work2, iwork, fac, nfac)

      
C  Make the PS wave functions on GRID
      call psdlgr (L, dalf, nps, ndim, c, f, f1, dnorm, enrgr, work,
     >   grid, nr, nld, dexpcut, jmin, jmax)
C  Check that the the wave functions are orthonormal. Useful only
C  if the PS states are sufficiently small at last R.
      call ortogc (f, nld, nr, nps, ch, ndim, grid(1,3), jmin, jmax)

c$$$      if (zinf.gt.5d-1.and.l.eq.0) then
c$$$C  NK should be even
c$$$         nk = 100
c$$$         print*,'Enter the number of K points (even) and rkp'
c$$$         read*,nk,rkp
c$$$         dk = 5d0/dfloat(nk)
c$$$         testk = 0d0
c$$$         test2k = 0d0
c$$$c$$$         rkp = 2.5d0
c$$$         do k = 1, nk
c$$$            w = dfloat( 2 * mod(k,2) + 2) * dk / 3.0
c$$$            if (k.eq.nk) w = w / 2d0
c$$$            rk = dfloat(k) * dk
c$$$            testk = testk + w * sin(rk)
c$$$            sum = 0d0
c$$$            sum2 = 0d0
c$$$            do n = 1, nps, 1
c$$$               ovlpf = 0d0
c$$$               ovlpi = 0d0
c$$$               do i = jmin(n), jmax(n)
c$$$                  ovlpf = ovlpf + sin(rk*grid(i,1))*f(i,n) * grid(i,3)
c$$$                  ovlpi = ovlpi + sin(rkp*grid(i,1))*f(i,n) * grid(i,3)
c$$$               enddo
c$$$               sum = sum + ovlpf * ovlpi
c$$$               if (k.le.1) print*, n, enrgr(n), sum
c$$$               if (n.le.nps/2) sum2 = sum
c$$$            enddo
c$$$            test2k = test2k + sum * w
c$$$            write(51,*) real(rk),real(sum), real(sum2)
c$$$         enddo
c$$$         print*,'tests:',testk/(1d0-cos(rk)),test2k
c$$$         call update(51)
c$$$         stop
c$$$      endif 
c$$$         ndelta = jmax(nps)
c$$$         do j = 1, ndelta
c$$$            do i = 1, ndelta
c$$$               tmp = 0d0
c$$$               do n = 1, nps
c$$$c$$$               do n = nabot(0), natop(0)
c$$$                  tmp = tmp + f(i,n) * f(j,n)
c$$$               enddo
c$$$               deltar(i,j) = tmp
c$$$            enddo
c$$$         enddo
c$$$         do n = 1, nps
c$$$            do i = 1, jmax(n)
c$$$               write(40+n,*) real(grid(i,1)), real(f(i,n)/grid(i,1))
c$$$            enddo
c$$$            call update(50+n)
c$$$         enddo
c$$$         stop
c$$$      endif 
C  Define the output arrays
      do n = 1, nps
         do i = 1, maxr
            ps(i,n) = 0.0
         enddo 
         if (f(jmin(n),n).gt.0.0) then
c$$$         if ((-1)**L*f(jmin(n),n).gt.0) then
            if (f(jmin(n)+1,n).le.0.0)
     >         print*,'change of sign in the first two points (MAKEPS)'
            do i = jmin(n), jmax(n)
               ps(i,n) = real(f(i,n))
            end do
C  define the coefficients for Alisher's program
            do k = 1, nps
               cknd(k,n) = c(k,n) * dnorm(k)
            enddo
         else
c$$$            print*,'Changed sign of wave function for N, L:',n+l,l
            if (f(jmin(n)+1,n).gt.0.0)
     >         print*,'change of sign in the first two points (MAKEPS)'
            do i = jmin(n), jmax(n)
               ps(i,n) = - real(f(i,n))
            end do
C  define the coefficients for Alisher's program
            do k = 1, nps
               cknd(k,n) = - c(k,n) * dnorm(k)
            enddo
         end if 
C  Energy is in Rydbergs
         energy(n) = 2.0 * real(enrgr(n))
C  The following modification is for generating the positronium energies
         if (Zinf.eq.0.5) energy(n) = 2.0 * energy(n)
      end do 
      end

C=====================================================================
C  Define pseudostates in Laguerre basis using expansion coefficients C.
C=====================================================================
      subroutine psdlgr (L, alf, nmax, nload, c, f, f1, dnorm,enrgr,pl,
     >   grid, nr, nld, epscut, jmin, jmax)
      include 'par.f'
      implicit real*8 (a-h, o-z)
      dimension c(nload, nmax), f(nld, nmax), dnorm(nmax), pl(nmax),
     >   grid(nr), jmin(nmax), jmax(nmax), f1(nld, nmax), enrgr(nmax)
      common /flogs/ faclog(1000)
      dimension powr(ncmax,ncmax,0:lamax), en(ncmax,0:lamax),
     >   alpha(0:lamax), nps(0:lamax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      real ps2, psen2,zasym
      save nps,alpha,powr ! save is used for orbcc.dat
C
C INPUT:
C -----      
C  L - orbital momentum.
C  alf - parameter of basis.
C  nmax - number of pseudostates.
C  nload, nld -  dimensions of arrays
C  C    - matrix of eigenvectors.
C  grid - "R" - grid of "NR" points
C
C OUTPUT:
C ------
C  f(i,n) - "n"-radial wave functions.
C  f1(i,n) = f(i,n) * r**(-L-1) * exp(alf*r) for laguerre integration.
C     jmin(n) - the first point of "n"-state.
C     jmax(n) - the last  point of "n"-state.
C  Pl - work space for the program.
C
      if (l.gt.lamax) STOP 'L > LAMAX in PSDLGR (laguer.f)'
      do n = 1, nmax
         en(n,l) = enrgr(n)
      enddo 
      nps(l) = nmax
      alpha(l) = alf
      L2 = 2 * L
      alfl = alf 
      a2 = 2.0d0 * alfl
C
C     Do some checks
C
      if (nr .gt. nld) then
         print*, 'number of points is more than dimension of GRID'
         print*, 'NR=', nr, ' NLD=', nld
         stop    'number of points is more than dimension of GRID'
      end if 
C
C     Define JMIN(N)
C
      xfirst = dexp( log(epscut/dnorm(1)) / dble(l+1))
      rfirst = xfirst / a2
      i = 1
      do while ((grid(i) .le. rfirst) .and. (i .lt. nr))
         i = i + 1
      end do 
      ifirst = 1 ! i
      do n = 1, nmax
         jmin(n) = ifirst
         jmax(n) = nr
      end do 
C  Get the Slater type orbital coefficients      
      do n = 1, nmax
         do m = 1, nmax
            powr(m,n,l) = 0d0
         enddo
      enddo

      do n = 1, nmax
         do m = 1, nmax
            const = c(m,n) * dnorm(m)
            a2new = 1d0
            powr(:,n,l) = 0d0
            do k = 0, m - 1
               if (abs(a2new).gt.1d300) cycle ! Igor 22/03/2016
               powr(k+1,n,l) = powr(k+1,n,l) + const * a2new *
     >            exp(faclog(m+2+l2) - faclog(k+1) - faclog(m-k) -
     >            faclog(3+l2+k))
               a2new = - a2new * a2
            enddo
         enddo
      enddo 
C
C     Loop by  r-grid
C
      cc1 = a2**(L+1)
c$$$      do i = 1, nr
      do i = nr, 1, -1
         r = grid(i)
         x = a2 * r
         c1 = dexp(-0.5d0 * x  +  dble(L+1) * log(x))
C
C        Define Laguerre polynomials
C
         pl(:) = 0d0
         pl1   = 1.0d0
         pl(1) = pl1
         pl2   = dble(L2 + 3) - x
         if (nmax.gt.1) pl(2) = pl2
         pl3 = 0d0
         do n = 3, nmax
            if (abs(pl3).gt.1d300) cycle ! Igor 22/03/2016
            pl3 = ((dble(2*n-1+L2)-x)*pl2 - dble(n+L2)*pl1) /
     >         dble(n-1)
            pl(n) = pl3
            pl1 = pl2
            pl2 = pl3
         end do
C
C        Loop by number of states.
C
         do n = 1, nmax
            f(i,n) = 0.0d0
c$$$            if (i .le. jmax(n)) then
               sum = 0.0d0
               do m = 1, nmax
                  sum = sum + c(m, n) * dnorm(m) * pl(m)
               end do
               f(i,  n) = c1  * sum
               f1(i, n) = cc1 * sum
C  The commented code below evaluates the wave function using a Slater 
C  representation and compares with the Laguerre form
c$$$               sum = 0d0
c$$$               do m = 1, nmax
c$$$                  sum = sum + powr(m,n,l) * r ** (m+l)
c$$$               enddo
c$$$               test = sum * exp(-alf*r) * a2 ** (l+1)
c$$$               if (abs((test-f(i,n))/(test+f(i,n))).gt.1e-3) print*,
c$$$     >            l,n,i,test,f(i,n)

               if (jmax(n).eq.nr.and.abs(f(i,n)).gt.epscut)
     >            jmax(n) = i
               
c$$$               if (i .gt. jmin(n)+3+100) then
c$$$                  if (     (abs(f(i-3, n))  .lt.  epscut)
c$$$     >               .and. (abs(f(i-2, n))  .lt.  epscut)
c$$$     >               .and. (abs(f(i-1, n))  .lt.  epscut)
c$$$     >               .and. (abs(f(i,   n))  .lt.  epscut))
c$$$     >               jmax(n) = i - 3
c$$$               end if 
c$$$            end if 
         end do 
      end do 
      do n = 1, nmax
         if (f(jmin(n),n).lt.0d0) then
            do m = 1, nmax
               powr(m,n,l) = - powr(m,n,l)
            enddo 
         endif
      enddo
      lprint = latop
      lprint = -1
      if (nmax.gt.50) lprint = -1
      if (l.eq.lprint) then
         open(42, file='orbcc.dat')
         do lp = 0, lprint
            do n = 1, nps(lp)
               write(42,*) lp+n,alpha(lp)
            enddo
         enddo
         do nn = 1, nps(0)
            do lp = 0, min(nn-1,lamax)
               n = nn - lp
c$$$               if (n.le.nps(lp)) then
               if (en(n,lp).lt.0.0.and.n.le.nps(lp)) then
                  write(42,*) lp,'  1   ',-en(n,lp),'          ',nn,lp
                  do m = 1, nps(lp)
                     write(42,*) (2.0*alpha(lp))**(lp+1)*powr(m,n,lp) /
C  The crazy factor below multiplies the orbital in the cross program
     >                  sqrt((2.0*alpha(lp))**(2*(m+lp)+1)/
     >                  exp(faclog(2*(m+lp)+1))) 
                  enddo
               endif 
            enddo 
         enddo
         close(42)
      endif
      return
      end

C============================================================      
C  ch(i,j) = < f(i) I f(j) > - delta(i,j)
C============================================================      
      subroutine ortogc (f, nload, nr, nmax, ch, nch, weight,
     >   jmin, jmax)
      implicit real*8 (a-h, o-z)
      dimension f(nload,nmax), ch(nch,nmax), weight(nr),
     >   jmin(nmax), jmax(nmax)
      ncount = 0
      err = 0d0
      do m = 1, nmax 
         do n = 1, nmax
            sum = 0d0
            r = 0d0
            func = 0d0
            do i = max(jmin(m), jmin(n)), min(jmax(m), jmax(n))
               w = weight(i)
               func1 = f(i,m)
               func2 = f(i,n)
               sum = sum + func1 * func2 * w
            end do 
            diag = 0.0d0
            if (n .eq. m)  diag = 1.0d0
            ch(n,m) = abs(sum - diag)
c$$$            print*,n,m,ch(n,m)
            if (ch(n,m).gt.err) err = ch(n,m)
c$$$            err = err + ch(n,m)
         end do 
         do n = 1, nmax
            if (abs(ch(n,m)).gt.1e-3) ncount = ncount + 1
c$$$     >         print*,'Warning: two PS states are not orthonormal'//
c$$$     >         ' or not zero at last R',
c$$$     >         n,m,ch(n,m)
         end do 
      end do
c$$$      if (ncount.gt.0) then
         print*,'Number of non orthonormal (to 1e-3) states:', ncount
         print*,'Max orthogonality error:',err
c$$$      endif 
      return
      end

C=====================================================================
C        Define matrix of some potential in a Laguerre basis.
C        Normalization done later.
C=====================================================================
      subroutine potl (exch, pot0, jmax, L, alf, potmtrx, nmax, 
     >    nload, f, pl, grid, weight, nr, nld)
      implicit real*8 (a-h, o-z)
      include 'par.f'
      logical exch
      real f1(maxr), f2(maxr), sres
      dimension f(nld, nmax), grid(nr), pot0(nld)
      dimension potmtrx(nload, nmax), pl(nmax), weight(nr)
      integer maxf(nmax)
C
C INPUT:
C -----      
C  L     - orbital momentum.
C  alf   - parameter of basis.                     alfl = alf 
C  nmax  - number of basis functions.              ------------------
C  nld   - dimension of r-arrays.
C  nload - dimension of matrices.
C  f, pl - work arrays
C OUTPUT:
C ------
C  potmtrx - matrix of the potential. 
C
      L2 = 2 * L
      L22 = L2 + 2
      alfl = alf 
      a2 = 2.0d0 * alfl
C
C     Do some checks
C
      if (nr .gt. nld) then
         print*, 'number of points is more than dimension of GRID'
         print*, 'NR=', nr, ' NLD=', nld
         print*, 'Stop in potl'
         stop    'number of points is more than dimension of GRID'
      end if 
c
c     Define Laguerre polynomials on r-grid
c
      f(:,:) = 0.0
      call lagpol8 (L22, a2, f, nld, nmax, grid, nr, pl)
      do n = 1, nmax
         do i = 1, nr
            r = grid(i)
            f(i,n) = exp(- alfl * r + dble(l + 1) * log(r)) * f(i,n)
         end do
         call minmaxi(real(f(:,n)),nr,i1,i2) !don't set : to 1 Li fails!
         maxf(n)=i2
c         print*, n, i1, i2
      end do 
c
c     Define Matrix of the potential
c
      do n = 1, nmax
         in = maxf(n)
         f1(:)=0.0
         if (exch) f1(1:in) = f(1:in,n) * weight(1:in)
         do m = 1, n
            im = maxf(m)
c            call int3_8 (f(1,n), 1, in, f(1,m), 1, im, 
c     >         pot0, 1, jmax, res, weight, nr)
            j = min(im,in,jmax)
            res = SUM(dble(pot0(1:j)*weight(1:j))*f(1:j,n)*f(1:j,m))
            
            if (exch) then
               f2(:)=0.0
               f2(1:im) = f(1:im,m) * weight(1:im)
c               do i = 1, nr
c                  f1(i) = (f(i,n) * weight(i))
c                  f2(i) = (f(i,m) * weight(i))
c               enddo 
               call fcexch(f1,nr,f2,sres,l)
               res = res + dble(sres)
            endif
            
            potmtrx(n,m) = res * a2**L22
            potmtrx(m,n) = potmtrx(n,m)
         end do
      end do 
      return
      end

      subroutine fcexch(f1,nr,f2,res,l)
      include 'par.f'
      real f1(nr),f2(nr),res
      real ps2,psen2,enpsinb,psinb,rpow1,rpow2,cntfug,fun(maxr),
     >   vnl(maxr)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common/nofcexch/ inofcexch, nfcmax

      res = 0d0
      if(inofcexch .eq. 1) return
      if (ntype.eq.-3) return
      do lac = 0, lactop
         do nac = lac + 1, nabot(lac) - 1
            if( nfcmax .lt. nac) exit
            if (lac.eq.linc.and.nac.eq.ninc) then
               elecnum = float(2 * lac)
               stop 'Should not be here' ! Igor added 24/10/2008
            else
               elecnum = float(2 * lac + 1)
            endif 
C The above form fails for sodium if latop = 0.
C In the form below we assume that there are no d states in the core
            minfun = 1
            maxfun = istoppsinb(nac,lac)
            do i = minfun, maxfun
               fun(i) = psinb(i,nac,lac) * f2(i)
            end do 
            do 20 ilt = -lac, lac, 2
               lt = l + ilt
               if (lt.lt.0.or.lt.gt.ltmax) go to 20
               call cleb(2*lac,2*lt,2*l,0,0,0,c)
               const = - elecnum * c * c /
     >            float(2 * l + 1)
               call form(fun,minfun,maxfun,rpow1(1,lt),rpow2(1,lt),
     >            minrp(lt),maxrp(lt),maxfun,vnl,i1,i2)
               tmp = 0d0
               do i = i1, min(i2,maxfun)
                  tmp = tmp + psinb(i,nac,lac) * vnl(i) * f1(i)
               end do
               res = res + tmp * const
 20         continue 
         enddo 
      enddo 
      return
      end 

c     Note: radial integration weights (gridr(i,3)) are included in functions f1(i) and f2(i)
      subroutine tmpfcexch(f1,min1,max1,nr,f2,res,l)
      include 'par.f'
      real f1(nr),f2(nr),res
      real ps2,psen2,enpsinb,psinb,rpow1,rpow2,cntfug,fun(maxr),
     >   vnl(maxr)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common/nofcexch/ inofcexch, nfcmax

      res = 0d0
      if(inofcexch .eq. 1) return
      if (ntype.eq.-3) return
      do lac = 0, lactop
         do nac = lac + 1, nabot(lac) - 1
            if( nfcmax .lt. nac) exit
C The above form fails for sodium if latop = 0.
C In the form below we assume that there are no d states in the core
            minfun = 1
            maxfun = istoppsinb(nac,lac)
            do i = minfun, maxfun
               fun(i) = psinb(i,nac,lac) * f2(i)
            end do 
            do 20 ilt = -lac, lac, 2
               lt = l + ilt
               if (lt.lt.0.or.lt.gt.ltmax) go to 20
               call cleb(2*lac,2*lt,2*l,0,0,0,c)
               const = - float(2 * lac + 1) * c * c /
     >            float(2 * l + 1)
               call form(fun,minfun,maxfun,rpow1(1,lt),rpow2(1,lt),
     >            minrp(lt),maxrp(lt),maxfun,vnl,i1,i2)
               tmp = 0d0
               do i = max(min1,i1), min(i2,maxfun,max1)
                  tmp = tmp + psinb(i,nac,lac) * vnl(i) * f1(i)
               end do
               res = res + tmp * const
 20         continue 
         enddo 
      enddo 
      return
      end 

      
C=====================================================================
C          Define basis functions in local potential -Z0/r + potmtrx(m,n)
C=====================================================================
      subroutine basisl (potmtrx, Z0, L, alf, nmax, nload, enrg, c, h,
     >   dnorm, dnorm2, work, iwork, fac, nfac)
      implicit real*8 (a-h, o-z)
      dimension c(nload, nmax), h(nload, nmax), fac(0:nfac),
     >   potmtrx(nload, nmax)
      dimension enrg(nmax), dnorm(nmax), dnorm2(nmax), work(8*nmax)
      dimension iwork(5*nmax),ifail(1000)
C
C      f(n,l,r) = dnorm(n,l) * (2*alfl*r)**(l+1) * exp(-alfl*r) *
C                 * Laguerre(2*l+2;n-1;2*alfl*r)
C
C INPUT:
C -----      
C  potmtrx   - matrix of some potential. 
C  L     - orbital momentum.
C  alf   - parameter of basis.              alfl = alf 
C  nmax  - number of pseudostates.          ------------------
C  nload - dimension of arrays
c  Z0    - The charge of the coulomb potential.
C
C OUTPUT:
C ------
C  enrg - eigenvalues
C  H    - Hamiltonian matrix.
C  C    - matrix of eigenvectors.
C  Work, dnorm2 - work space for the program
C
      L2 = 2 * L
      alfl = alf 
      a2 = 2.0d0 * alfl
c
c     Do some checks
c
      if (Z0 .lt. 0d0) then
         print*, 'Z0 is out of range,  Z0=', Z0
         stop    'Stop in BASISL'
      end if 
C
C     Define factorials as n!=dexp(fac(n))
C
      fac(0) = 0.0d0
      fact = 0.0d0
      tmp = 1.0d0
      ntmp = 2 * nmax + L2
      if (nfac .lt. ntmp) then
         print*, 'nfac is not enough to store array of factorials,'
         print*, 'nfac has to be more than   2 * NMAX + 2 * L =', ntmp,
     >      ' but  nfac =', nfac
         stop 'nfac has to be more than   2 * NMAX + 2 * L'
      end if 
      do n = 1, ntmp
         fact   = fact + log(dble(n))
         fac(n) = fact
      end do
C
C     Define normalization coeff. Dnorm 
C
      c1 = sqrt(a2)
      c2 = 1.0d0 / a2
      do i = 1, nmax
         dnorm(i)  = c1 * exp(0.5d0 * (fac(i - 1)  -  fac(L2 + 1 + i)))
         dnorm2(i) = c2 * exp(fac(L2 + i)  -  fac(i - 1))
      end do 
C
C     Define Hamiltonian matrix 
C
      c2 = -a2 * a2 * 0.5d0
      do i = 1, nmax
         do j = 1, i
            sm2 = 0.0d0
            do jj = 1, j - 1
               do jjj = 1, min(i, jj)
                  sm2 = sm2 + dnorm2(jjj)
               end do
            end do
            sm1 = 0.0d0
            do jj = 1, min(i, j)
               sm1 = sm1 + dnorm2(jj)
            end do 
            diag = 0.d0
            if (i .eq. j)  diag = 0.25d0
            potlz = c2 * Z0 / alfl * sm1
            potl2  =  potmtrx(i,j)
            res = (dnorm(i) * dnorm(j) * (-dble(l+j) * sm1
     >         + sm2) + diag) * c2 
     >         + dnorm(i) * dnorm(j) * (potl2 + potlz)
            h(i, j) = res
            h(j, i) = res
         end do
      end do 
C         
C  if matc = 0, then  only eigenvalues.
      matc = 1
      call rs(nload, nmax, h, enrg, matc, c, dnorm2, work, ierr)
      lwork = 8 * nmax
c      call dsyevx('V','A','U',nmax,h,nload,0d0,0d0,0,0,0d0,nfound,
c     >   enrg,c,nload,work,lwork,iwork,ifail,info)
c      if (info.ne.0) then
c         print*,'INFO:',info
c         stop
c      endif 
      if (ierr .ne. 0)  then
         print*, 'Program "RS" finished abnormaly, ierr =', ierr
         stop    'Program "RS" finished abnormaly'
      end if 
      return
      end
C===================================================================
C                                     m
C   Laguerre's polinomials  f(i,n) = L   (dlambda * grid(i))
C                                     n-1
C===================================================================
      subroutine lagpol8 (m, dlambda, f, nload, nmax, grid, nr, pl)
      implicit real*8 (a-h, o-z)
      dimension f(nload, nmax), grid(nr),  pl(nmax)
C
C INPUT:
C -----
C  m - parameter of Laguerre polinomial.
C  nmax  - the max number of n.
C  nload - dimension of f.
C  grid  - "R" - grid of "NR" points.
C  pl(nmax) - work space for this program.
C
C OUTPUT:
C ------
C  f(i,n) - Laguerre polinomials.
C
      L2 = m - 2
C
C     Loop by  r-grid
C
c$$$      pl(:) = 1d0 !Igor 22/03/2016
      do i = 1, nr
         r = grid(i)
         x = dlambda * r
         f(i,:) = 0d0  !Igor 22/03/2016
C
C        Define Laguerre's polinomials
C        and store them in array 'f'
C
         pl1    = 1.0d0
         pl(1)  = pl1
         f(i,1) = pl1
c
         pl2    = dble(L2 + 3) - x
         if (nmax.gt.1) then
            pl(2)  = pl2
            f(i,2) = pl2
         endif 
c
         pl3 = 0d0
         do n = 3, nmax
            if (abs(pl3).gt.1d300) cycle !Igor 22/03/2016
            pl3 = ((dble(2*n-1+L2)-x)*pl2 - dble(n+L2)*pl1) /
     >         dble(n-1)
c$$$            pl(n) = pl3  !Igor 22/03/2016
            f(i,n) = pl3
            pl1 = pl2
            pl2 = pl3
         end do
      end do 
      return
      end
C
C==================================================================
C
C   int3_8    Program for integration of three real*8 functions.
C
C      call int3_8 (f1, nf1, nl1, f2, nf2, nl2, f3, nf3, nl3,
C     >   res, weight, nr)
C
C==================================================================
C   Real*8  Programe for integration of three functions
C==================================================================
      subroutine int3_8 (f1, nf1, nl1, f2, nf2, nl2, f3, nf3, nl3,
     >   res, weight, nr)
      implicit real*8 (a-h, o-z)
      dimension f1(nr), f2(nr), f3(nr), weight(nr)
      if ((nl1 .gt. nr)  .or.  (nl2 .gt. nr) .or.  (nl3 .gt. nr)) then
         print*, 'Number of points are not enough in r-mesh'
         print*, 'the last point of "F1"=', nl1
         print*, 'the last point of "F2"=', nl2
         print*, 'the last point of "F3"=', nl3
         print*, 'the last point of "WEIGHT"=', nr
         stop 'stop  in  "INT3_8"'
      end if 
      res = 0.0d0
      do i = max(nf1, nf2, nf3), min(nl1, nl2, nl3)
         res = res + f1(i) * f2(i) * f3(i) * weight(i)
      end do 
      return
      end

c$$$      program psdrive
c$$$      include 'par.f'
c$$$      real waveout(maxr,nnmax),erydout(nnmax)
c$$$      dimension jdouble(100)
c$$$*
c$$$      hmax = 0.03488889
c$$$      ra = 100.0
c$$$      njdouble = 11
c$$$      jdouble(1) = 1
c$$$      do 10 i=2,10
c$$$       jdouble(i) = (i-1)*32
c$$$ 10   continue
c$$$      jdouble(11) = 8000
c$$$      zas = 1.0
c$$$      la =  1
c$$$      nbmax = 20
c$$$      call pseudo(jdouble,njdouble,hmax,zas,la,ra,nbmax,maxr,
c$$$     >                  erydout,waveout,lastpt)
c$$$      do 20 ne=1,nbmax
c$$$       write(6,*) ne,erydout(ne)
c$$$ 20   continue
c$$$      print *,'lastpt = ',lastpt
c$$$      stop
c$$$      end
*======================================================================*
      subroutine pseudo(jdouble,njdouble,hmax,zas,la,ra,nbmax,maxrr,
     >                  vnucl,erydout,waveout,lastpt)
************************************************************************
*                                                                      *
*   THIS ROUTINE DETERMINES PSEUDO BOUND STATES OF A RADIAL POTENTIAL, *
*   I.E., STATES ORBITALS THAT VANISH AT R=0 AND R=A.                  *
*                                                                      *
*   IT IS THE BOUND-STATE PART OF THE CORE POTENTIAL ITERATION         *
*   PROGRAM PUBLISHED IN CHAPTER 2 OF                                  *
*                                                                      *
*            COMPUTATIONAL ATOMIC PHYSICS                              *
*            KLAUS BARTSCHAT  (ED.)                                    *
*            SPRINGER (1996)                                           *
*                                                                      *
*            WRITTEN BY:   KLAUS BARTSCHAT                             *
*                          PHYSICS DEPARTMENT                          *
*                          DRAKE UNIVERSITY                            *
*                          DES MOINES, IOWA 50311, U.S.A.              *
*                                                                      *
*            LAST UPDATE:  MARCH 17, 2003                              *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*
      include 'par.f'
      PARAMETER (NDIM1=4,NDIM2=maxr+1,NDIM3=100,LLMAX=lamax)
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,SIX=6.0D0)
      PARAMETER (HUGE=1.0D30,SMALL=1.0D-20)
*
      real hmax,zas,ra,vnucl(maxrr)
      real waveout(maxr,nnmax),erydout(nnmax)
      dimension jdouble(njdouble)
*
      real enpsinb,psinb,rpow1,rpow2,cntfug,fun(maxr),
     >   vnl(maxr),c
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      real zasym
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)

      COMMON / PATH / RFULL(NDIM2),RHALF(2*NDIM2),Y(NDIM1,NDIM2)
c$omp threadprivate(/PATH/)
      COMMON / POT  / EGUESS,CORE(2*NDIM2),RDMASS
c$omp threadprivate(/POT/)
*
      DIMENSION IHX(NDIM3),IRX(NDIM3)
      DIMENSION NMINL(0:LLMAX),NMAXL(0:LLMAX)
      DIMENSION VSTART(NDIM1),EBOUND(NNMAX),
     +          OVRLAP(NNMAX,NNMAX)
      DIMENSION FVALUE(NDIM2*2),vtemp(ndim2*2)
      DIMENSION WAVES(NNMAX,2*NDIM2), corei(2*ndim2)
*
**********************  FORMAT STATEMENTS  *****************************
1000  FORMAT(' *******************************************************',
     +     /,' *                                                     *',
     +     /,' *            BOUND-STATE CALCULATION                  *',
     +     /,' *            -----------------------                  *',
     +     /,' *******************************************************',
     +       ///,' NUCLEAR CHARGE:                   ',I3,/,            
     +           ' ASYMPTOTIC CHARGE:                ',I3,/,            
     +           ' ERROR LIMIT FOR BOUND STATE ENERGY:',1P,D12.3,/,     
     +           ' LOWER LIMIT FOR BOUND STATE ENERGY:',1P,D12.3,/)
1001  FORMAT(/,' NEW BOUND STATE FOUND:   EBOUND(',I3,I3,') = ',          
     +       1P,D17.9,//,' NEXT GUESS = ',1P,D14.7,/)
1002  FORMAT(/,' WARNING: OVERLAP INTEGRAL BETWEEN FUNCTIONS WITH', 
     +         ' N1 = ',I2,' AND N2 = ',I2,' FOR L = ',I1,' :',/, 
     +           1P,D15.8,' +/- ',1P,D15.8)
1003  FORMAT(/,' DEBUG PRINT OUT OF FIRST 100 POINTS:',/,               
     +           8X,'R',10X,'FUNCTION',//,1X,1P,2D14.6)
1004  FORMAT(1X,1P,8D14.6)
1005  FORMAT(/,' INODE,NNODE,NBOUND,LBOUND = :',4I5,/,                  
     +         ' ELOW,EHIGH,EGUESS  = :',1P,3D16.8,/)
1006  FORMAT(/,' NITER = ',I5,'  GREATER THAN  NITMAX = ',I5,            
     +         '.  PROGRAM STOPPED.',/)
1007  FORMAT(//,' MESH:',/,' --------------',/,                
     +         ' HINT = ',1P,D14.6,'  NIX =',I2,'  NSTEP  =',I6,        
     +         '  LAST POINT = ',1P,D14.6)
1008  FORMAT(' IHX-ARRAY:',20I5)
1009  FORMAT(' IRX-ARRAY:',20I5)
1010  FORMAT(//,' STARTING BOUND STATE ITERATION FOR (NBOUND,LBOUND) =',
     +         ' (',I2,',',I2,')')
1011  FORMAT(/,' PROBLEMS WITH DEFINITION OF NINTEG. NEND = ',I6,/,     
     +         ' PROGRAM STOPPED.',/)
1012  FORMAT(/,' NORMALIZED ORBITAL FOR  (NBOUND,LBOUND) = (',I2,       
     +         ',',I2,') :',/,'        R          ORBITAL')
1013  FORMAT(1P101D10.2)
1014  FORMAT(/,' LMIN = ',I2,'  LMAX = ',I2,/)
1015  FORMAT(/,' NMIN = ',I2,'  NMAX = ',I2,'  FOR L = ',I2,/)
1016  FORMAT(/,' NITMAX = ',I3,'  IBUG = ',I2,'  ITOUT = ',I2,/)
1017  FORMAT(/,' ESTART = ',1PD14.6,'  ERROR = ',1PD14.6,/)
1018  FORMAT(/,' NUCLEAR CHARGE = ',F6.2,'  REDUCED MASS = ',F7.1,/)
1019  FORMAT(/,' OVERLAP INTEGRALS FOR L = ',I2,/)
1020  FORMAT(/,' NEW BOUND STATE FOUND:   EBOUND(',I3,I3,') = ',          
     +       1P,D17.9,/)
1021  FORMAT(//,' **********   BOUND STATE CALCULATION   **********',//,
     +       ' INPUT IPOT:  1 = COULOMB,  OTHERWISE = NUMERIC',//)
1022  FORMAT(' COULOMB PROBLEM',/)
1023  FORMAT(' NUMERICAL POTENTIAL',/)
1025  FORMAT(' IPOT = ',I2,' CURRENTLY INVALID.  PROGRAM STOPS.',//)
1026  FORMAT(' R0 = ',1PD14.6,'  D = ',1PD14.6,'  ALPHA = ',1PD14.6,0P,
     +       '  REDUCED MASS = ',F7.1,/)
1027  FORMAT(' V0 = ',1PD14.6,'  A = ',1PD14.6,
     +       '  REDUCED MASS = ',F7.1,/)
1028  FORMAT(//,'# BOUND-STATE ENERGIES FOUND FOR L = ',I2,/,
     +          '#      N       E (a.u)        E (eV)')
1029  FORMAT(I8,5X,F9.5,5X,F9.3)
2019  FORMAT(/,'WORST OVERLAP:',1PD12.4)
************************************************************************
*
      oworst = zero
      znuc = zas + 1   ! restored "+1" for box-based target, but breaks DW
      print*,'ZNUC in PSEUDO:',znuc
      itout = 0
      ibug = 0
      ipot = 1
      rdmass = 1.0d0
      nix = njdouble-1
      hint = hmax*0.5d0**(nix-1)
      do 10 i=1,nix
       ihx(i) = 2**(i-1)
       irx(i) = jdouble(i+1)
!       print*, 'i,irx(i):',i,irx(i)
 10   continue
      y(:,:) = 0d0
      waves(:,:) = 0.0
      NSTEP = IRX(NIX)+1
*
*  CALCULATE ARRAYS FOR MESH
*
      RHALF(1)=ZERO
      RFULL(1)=ZERO
      DO 20 I=1,NIX
       HSTEP = HINT*IHX(I)
       HHALF=HALF*HSTEP
       IF (I.EQ.1) THEN
        JBEG = 1
       ELSE
        JBEG = IRX(I-1)+1
       ENDIF
       DO 20 J=JBEG,IRX(I)
        RFULL(J+1) = RFULL(J)+HSTEP
        I2=J+J
        RHALF(I2) = RFULL(J)+HHALF
        RHALF(I2+1) = RFULL(J+1)
20     CONTINUE
30    CONTINUE
*
*  redefine the last interval based on RA
*
      nstepsave = nstep
      do 35 j=1,nstep
*       print *,j,ra,rfull(j)
       if (ra.le.rfull(j)) goto 36
 35   continue
      if (j.eq.nstep+1) then
         ra = rfull(nstep)
         j = nstep
         print*,'CAUTION: reset R0 to:',ra
      endif 
 36   irx(nix) = j/2*2
      nstep = irx(nix)+1
      print *,'mesh reset to npts,ra = ',nstep-1,rfull(nstep)
      lastpt = nstep - 1
      error = 1.0d-12
      nitmax = 100
      estart = -0.5d0*znuc**2 - 0.1
      elow = estart
*      if (la.lt.0) stop
*
*  PRINT MESH PARAMETERS
*
c$$$      WRITE(6,1007) HINT,NIX,NSTEP,RFULL(NSTEP)
c$$$      WRITE(6,1008) (IHX(I),I=1,NIX)
c$$$      WRITE(6,1009) (IRX(I),I=1,NIX)
*
*  SET POTENTIAL IN ATOMIC UNITS
*
      IF (IPOT.EQ.1) THEN
         CORE(1) = -HUGE
         DO 40 I=2,2*NSTEP-1
            CORE(I) = -ZNUC/RHALF(I) 
 40      CONTINUE
c$$$       do i = 2, 2*NSTEP-1, 2*NSTEP-3
c$$$          print*,'i,rhalf(i),',i,rhalf(i)
c$$$       enddo
c$$$       do j = 1, 2
c$$$          print*,'J, RFULL(J), VNUCL(J)*RFULL(J+1):',
c$$$     >       j,rfull(j),vnucl(j)*rfull(j+1)
c$$$       enddo
! The following has been checked to work for bound states of distorting potentials
c$$$         if (vnucl(1).ne.0.0.and..false.) then !"Needed for DW, see above"
         if (vnucl(1).ne.0.0) then !"Needed for DW, see above"
            do j = 1, nstep-1
               vtemp(j+1) = vnucl(j)*rfull(j+1)
c$$$             write(75,*) rfull(j+1),vnucl(j)
            enddo
            vtemp(1) = dfloat(nint(vtemp(2)))
            call intrpl(nstep,rfull,vtemp,2*nstep-1,rhalf,fvalue)
            do j = 2, 2*nstep-1
               core(j) = core(j) + fvalue(j)/rhalf(j)/2.0 ! a.u.
            enddo
c$$$          do j = 2, 5
c$$$             print*,j,rfull(j)*vnucl(j-1),fvalue(2*j-1)
c$$$          enddo
            znapp = -vtemp(1)/2
            fvalue(:)=0d0
c$$$            estart = -0.5d0*(znuc+1)**2 - 0.1
            estart = -0.5d0*(znapp+1)**2 - 0.1
            print*,'check: ZNAPP,ZNUC,estart(ry)',znapp,znuc,estart*2.0
            elow = estart
c$$$          do j = 1, 2*nstep-1
c$$$             write(76,*) rhalf(j), core(j)
c$$$          enddo 
         endif 
      ENDIF
      do i = 2, 2*nstep-1
         corei(i) = core(i)  ! Intended for Hartree-Fock not yet used. 
      enddo 
 

*
***      if (la.eq.-1) stop
*
*  HERE WE COULD PROVIDE ANY POTENTIAL ON THE RHALF-MESH
*
*
*  START ITERATION FOR BOUND STATES
*
      DO 170 LBOUND=la,la
       EGUESS = ESTART
       EHIGH = TWO
       NMIN = la+1
       NMAX = nbmax
       IBCNT = 0
       DO 140 NBOUND=NMIN,NMAX
        IF (IBUG.GT.0) WRITE(6,1010) NBOUND,LBOUND
        NNODE = NBOUND-LBOUND-1
        NITER = 0
        ELOW  = ESTART
        IF (NBOUND.GT.NMIN) ELOW = EBOUND(NBOUND-1)
*
*  SET STARTING VALUES FOR RUNGE-KUTTA INTEGRATION.
*  FUNCTION: 0.0  DERIVATIVE: 1.0
*
50      VSTART(1) = ZERO
        VSTART(2) = ONE
        NCLASS = 0

! The following attempted to include the Fock part, but iteration is tricky
c$$$        NINTEG = 1
c$$$        factor = 1d0
c$$$        if (abs(y(1,ninteg)).lt.1e10*factor) then
c$$$           if (niter.gt.0) then
c$$$              i = nend-1
c$$$              do while (i.gt.1.and.abs(Y(1,I)).le.abs(Y(1,I+1)))
c$$$                 i = i - 1
c$$$              enddo 
c$$$              NINTEG = I
c$$$              DO I=1,ninteg
c$$$                 FVALUE(I) = Y(1,I)*Y(1,I)
c$$$              enddo
c$$$              y(1,ninteg:) = 0d0
c$$$              CALL INTEGR(RFULL,FVALUE,NINTEG,RESULT,ERRINT)
c$$$              FACTOR = ONE/DSQRT(RESULT)
c$$$           endif 
c$$$        vtemp(:) = 0d0
c$$$        do lac = 0, lactop
c$$$           do nac = lac + 1, nabot(lac) - 1
c$$$              minfun = 1
c$$$              maxfun = istoppsinb(nac,lac)
c$$$              do i = minfun, maxfun
c$$$                 fun(i) = psinb(i,nac,lac) * y(1,i+1)*factor
c$$$              end do 
c$$$              do ilt = -lac, lac, 2
c$$$                 lt = la + ilt
c$$$                 if (lt.lt.0.or.lt.gt.ltmax) cycle
c$$$                 call cleb(2*lac,2*lt,2*la,0,0,0,c)
c$$$                 const = - float(2 * lac + 1) * c * c /
c$$$     >              float(2 * la + 1)
c$$$                 call form(fun,minfun,maxfun,rpow1(1,lt),rpow2(1,lt),
c$$$     >              minrp(lt),maxrp(lt),maxfun,vnl,i1,i2)
c$$$                 do i = max(minfun,i1), min(i2,maxfun)
c$$$                    vtemp(i+1) = vtemp(i+1) +
c$$$     >                 const * psinb(i,nac,lac) * vnl(i)
c$$$                 end do
c$$$              enddo 
c$$$           enddo 
c$$$        enddo
c$$$        vtemp(1) = 0d0
c$$$        fvalue(:) = 0d0
c$$$        call intrpl(nstep,rfull,vtemp,2*nstep-1,rhalf,fvalue)
c$$$        do j = 2, 2*nstep-1
c$$$           core(j) = corei(j) + fvalue(j)
c$$$        enddo
c$$$        print*,'nbound,niter,elow,ehigh,factor,ninteg,core(100):',
c$$$     >     nbound,niter,elow,ehigh,factor,ninteg,core(100)
c$$$        endif 

        
        CALL RKDUMB(VSTART,2,NSTEP,NEND,LBOUND,NCLASS)
        IF (IBUG.GT.3) WRITE(6,1003) RFULL(1),Y(1,1)
*
*  COUNT NODES
*
        INODE = 0
        DO 60 I=2,NEND
         IF (IBUG.GT.3.AND.I.LE.100) WRITE(6,1004) RFULL(I),Y(1,I)
         IF (Y(1,I-1)/Y(1,I).LT.ZERO) INODE = INODE+1 ! Igor changed * to / to avoid overflows 19/4/04
60      CONTINUE
*
*  CHECK WHETHER ENERGY IS TOO LOW OR TOO HIGH AND DEFINE
*  A NEW VALUE FOR EGUESS.
*
        IF (INODE.GT.NNODE) THEN
         EHIGH = EGUESS
        ELSE
         ELOW = EGUESS
        ENDIF
        EGUESS = HALF*(EHIGH+ELOW)
        IF (IBUG.GT.1) WRITE(6,1005) INODE,NNODE,NBOUND,LBOUND,       
     +                               ELOW,EHIGH,EGUESS
c$$$        IF (ABS(EHIGH-ELOW).GT.ERROR) THEN
        IF (ABS(ELOW/EHIGH-1d0).GT.ERROR) THEN !Igor increased accuracy 07/2004
         NITER = NITER+1         
c$$$         print*,niter,elow,ehigh
         IF (NITER.GT.NITMAX) THEN
          WRITE(6,1006) NITER,NITMAX
          STOP
         ENDIF
         GOTO 50
        ELSE
         IBCNT=IBCNT+1
         EBOUND(NBOUND) = EGUESS
*
*  NORMALIZE THE WAVEFUNCTION
*
         DO 70 I=1,NEND
          FVALUE(I) = Y(1,I)*Y(1,I)
70       CONTINUE
         DO 80 I=NEND-1,1,-1
          IF (FVALUE(I).GT.FVALUE(I+1)) GOTO 90
80       CONTINUE
         print*,'i,fvalue(i),fvalue(i+1):',i,fvalue(i),fvalue(i+1)
         WRITE(6,1011) NEND
         STOP
90       NINTEG = I
         CALL INTEGR(RFULL,FVALUE,NINTEG,RESULT,ERRINT)
         FACTOR = ONE/DSQRT(RESULT)
         DO 100 I=1,NSTEP
          IF (I.LE.NINTEG) THEN
           WAVES(NBOUND,I) = Y(1,I)*FACTOR
          ELSE
           WAVES(NBOUND,I) = ZERO
          ENDIF
100      CONTINUE
         IF (IBUG.GT.2) THEN
          WRITE (6,1012) NBOUND, LBOUND
          DO 110 I=1,NINTEG+1
           WRITE(6,1004) RFULL(I),WAVES(NBOUND,I)
110       CONTINUE
         ENDIF
*
*  END OF NORMALISATION. SET THE NEXT GUESS IN SOME REASONABLE WAY.
*
         IF (EGUESS.LT.ZERO) THEN
          EGUESS = HALF*EGUESS
         ELSE 
          EGUESS = TWO*EGUESS
         ENDIF
         XKN = SQRT(TWO*ABS(EGUESS))
         XKNP1 = XKN+5.0/RFULL(NSTEP)
         EHIGH = XKNP1**2+ONE
c$$$         IF (NBOUND.LT.NMAX) THEN
c$$$          WRITE(6,1001) NBOUND,LBOUND,EBOUND(NBOUND),EGUESS
c$$$         ELSE
c$$$          WRITE(6,1020) NBOUND,LBOUND,EBOUND(NBOUND)
c$$$         ENDIF
        ENDIF
*
        IF (IBCNT.GE.NDIM3) GOTO 140
*
*  CALCULATE THE OVERLAP INTEGRALS
*
        OWORST = ZERO
        DO 130 NNN=NMIN,NBOUND
         DO 120 I=1,NDIM2
          FVALUE(I) = WAVES(NBOUND,I)*WAVES(NNN,I)
120      CONTINUE
         CALL INTEGR(RFULL,FVALUE,IRX(NIX)+1,RESULT,ERRINT)
         IF (ABS(RESULT).GT.1.0D-04.AND.ABS(RESULT-ONE).GT.1.0D-05) 
     +       WRITE(6,1002) NNN,NBOUND,LBOUND,RESULT,ERRINT
         OVRLAP(NNN,NBOUND) = RESULT
         IF (NNN.EQ.NBOUND) THEN
          OWORST = MAX(OWORST,ABS(RESULT-ONE))
         ELSE
          OWORST = MAX(OWORST,ABS(RESULT))
         ENDIF
130     CONTINUE
140    CONTINUE
*
*  PRINT THE OVERLAP INTEGRALS FOR THIS VALUE OF L
*
c$$$       WRITE(6,1019) la
c$$$       DO 150 NNN=NMIN,NMAX
c$$$        WRITE(6,1013) (OVRLAP(N,NNN),N=NMIN,NNN)
c$$$150    CONTINUE
       WRITE(6,2019) OWORST
       if (oworst.gt.1e-3) stop 'Box-basis orthonormality failed'
*
*  WRITE THE ORBITALS TO A FILE (ONE UNIT PER L) IF ITOUT.NE.0
*
       IF (ITOUT.NE.0) THEN
        IT = ITOUT+LBOUND
        WRITE(IT,1013) RFULL(IRX(NIX)+1),
     +                 (EBOUND(NNN)*2.0D0,NNN=NMIN,NMAX)
        DO 160 I=1,IRX(NIX)+1
         WRITE(IT,1013) RFULL(I),(WAVES(NNN,I),
     +                            NNN=NMIN,NMAX)
160     CONTINUE
       ENDIF
170   CONTINUE  
*
*  PRINT ALL BOUND STATE ENERGIES FOUND
*
c$$$       WRITE(6,1028) la
       DO 180 N=nmin,nmax
c$$$        WRITE(6,1029) N,EBOUND(N),EBOUND(N)*27.21D0
        erydout(n-la) = EBOUND(N)*2.0d0
        do 185 i=2,IRX(NIX)+1
         waveout(i-1,n-la) = waves(n,i)
 185    continue
180    CONTINUE
*
      RETURN
      END
************************************************************************
      SUBROUTINE RKDUMB(VSTART,NVAR,NSTEP,NEND,LBOUND,NCLASS)
************************************************************************
*                                                                      *
*    THIS SUBROUTINE USES THE FOURTH ORDER RUNGE-KUTTA METHOD TO       *
*    PROPAGATE THE SOLUTION OF NVAR COUPLED DIFFERENTIAL EQUATIONS     *
*    KNOWN AT X1 (STORED IN VSTART) TO X2 BY TAKING NSTEP STEPS OF     *
*    EQUAL LENGTH (X2-X1)/NSTEP. THE RESULTS ARE STORED IN THE         *
*    COMMON BLOCK / PATH /.                                            *
*                                                                      *
************************************************************************
*                                                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      include 'par.f'
      PARAMETER (NDIM1=4,NDIM2=maxr+1)
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,SIX=6.0D0)
*
      COMMON / PATH / RFULL(NDIM2),RHALF(2*NDIM2),Y(NDIM1,NDIM2)
      COMMON / POT  / EGUESS,CORE(2*NDIM2),RDMASS
*
      DIMENSION VSTART(NVAR),V(NDIM1),DV(NDIM1)
*
*  SET THE CLASSICAL TURNING POINT FOR ZERO ANGULAR MOMENTUM.
*  THIS ONLY MAKES SENSE FOR NEGATIVE ENERGIES
      IF (EGUESS.LT.ZERO) THEN
       DO 10 I=2*NSTEP-1,1,-1
        IF (EGUESS.GT.CORE(I)) THEN
         RCLASS = RHALF(I)
         NCLASS = I/2+1
         GOTO 20
        ENDIF
10     CONTINUE
*
*  THE CLASSICAL TURNING POINT IS AT THE ORIGIN. SOMETHING MUST 
*  BE WRONG. STOP. 
*
       WRITE(6,1000)
       STOP
      ELSE
       RCLASS = RFULL(NSTEP)
       NCLASS = NSTEP
      ENDIF
20    CONTINUE
*
*  LOAD STARTING VALUES
*
      DO 30 I=1,NVAR
       V(I) = VSTART(I)
       Y(I,1) = V(I)
30    CONTINUE
      NEND = NSTEP
*
      DO 50 K=1,NSTEP
       INDEX = K+K-1
       X = RHALF(INDEX)
       if (abs(v(1)).lt.1d100) then !Igor added this to avoid overflows 8/2004
          CALL DERIVS(INDEX,V,DV,LBOUND)
          H = RHALF(INDEX+2)-RHALF(INDEX)
          CALL RK4(V,DV,NVAR,INDEX,H,V,LBOUND)
       endif
       DO 40 I=1,NVAR
          Y(I,K+1) = V(I)
 40    CONTINUE
*     IF (X.GT.RCLASS.AND.ABS(Y(1,K+1)).GT.ABS(Y(1,K))) THEN
*        NEND = K+1
*        RETURN
*       ENDIF
50    CONTINUE
*      IF (EGUESS.LT.ZERO) WRITE(6,1001)
*
1000  FORMAT('1',//,' PROGRAM STOPS IN SUBROUTINE RKDUMB, SINCE THE',
     +              ' CLASSICAL TURNING POINT COULD NOT BE DEFINED.',//)
1001  FORMAT(//,' WARNING : SOLUTION HAS NOT TURNED AROUND YET !',/)            
      RETURN
      END
************************************************************************
      SUBROUTINE RK4(Y,DYDX,N,INDEX,H,YOUT,LBOUND)
************************************************************************
*                                                                      *
*  THIS SUBROUTINE PROPAGATES THE SOLUTION FOR N VARIABLES Y BY        *
*  BY ONE STEP H USING THE FOURTH-ORDER RUNGE-KUTTA METHOD. DYDX       *
*  ARE THE DERIVATIVES AND Y THE VALUES OF THE FUNCTIONS AT X.         *
*                                                                      *
************************************************************************
*                                                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      include 'par.f'
      PARAMETER (NDIM1=4,NDIM2=maxr+1)
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,SIX=6.0D0)
*
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NDIM1),DYT(NDIM1),DYM(NDIM1)
*
      HH = H*HALF
      H6 = H/SIX
*
*  FIRST STEP
*
      DO 10 I=1,N
       YT(I) = Y(I)+HH*DYDX(I)
10    CONTINUE
      CALL DERIVS(INDEX+1,YT,DYT,LBOUND)
*
*  SECOND STEP
*
      DO 20 I=1,N
       YT(I) = Y(I)+HH*DYT(I)
20    CONTINUE
      CALL DERIVS(INDEX+1,YT,DYM,LBOUND)
*
*  THIRD STEP
*
      DO 30 I=1,N
       YT(I) = Y(I)+H*DYM(I)
       DYM(I) = DYT(I)+DYM(I)
30    CONTINUE
      CALL DERIVS(INDEX+2,YT,DYT,LBOUND)
*
*  FOURTH STEP
*
      DO 40 I=1,N
       YOUT(I) = Y(I)+H6*(DYDX(I)+DYT(I)+TWO*DYM(I))
40    CONTINUE
      RETURN
      END
************************************************************************
      SUBROUTINE DERIVS(INDEX,YRK,F,LBOUND)
************************************************************************
*                                                                      *
*  THIS SUBROUTINE CALCULATES THE DERIVATIVES IN THE RUNGE-KUTTA       *
*  AT THE MESHPOINT NUMBER GIVEN BY INDEX.                             *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'par.f'
      PARAMETER (NDIM1=4,NDIM2=maxr+1)
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,SIX=6.0D0)
      PARAMETER (HUGE=1.0D30,SMALL=1.0D-20)
*
      COMMON / PATH / RFULL(NDIM2),RHALF(2*NDIM2),Y(NDIM1,NDIM2)
      COMMON / POT  / EGUESS,CORE(2*NDIM2),RDMASS
*
      DIMENSION YRK(NDIM1),F(NDIM1)
*
*  CALCULATE THE DERIVATIVES.
*
      X = RHALF(INDEX)
      F(1) = YRK(2)
      IF (X.NE.ZERO) THEN
       F(2) = (TWO*RDMASS*(CORE(INDEX)-EGUESS)
     +        +(LBOUND*(LBOUND+1)/X**2))*YRK(1)
      ELSE
       F(2) = -HUGE*YRK(1)
      ENDIF
      RETURN
      END
************************************************************************
      SUBROUTINE INTEGR(X,Y,N,RESULT,ERROR)
************************************************************************
*                                                                      *
*  THIS SUBROUTINE INTEGRATES ARBITRARILY SPACED DATA                  *
*                                                                      *
*                      INPUT                                           *
*                      -----                                           *
*                                                                      *
*  X .........     VECTOR CONTAINING THE MESHPOINTS                    *
*                  (EITHER IN ASCENDING OR DESCENDING ORDER)           *
*  Y .........     VECTOR CONTAINING THE FUNCTION VALUES               *
*  N .........     NUMBER OF MESHPOINTS   (AT LEAST 4)                 *
*                                                                      *
*                      OUTPUT                                          *
*                      ------                                          *
*                                                                      *
*  RESULT ....     RESULT OF THE INTEGRATION                           *
*  ERROR .....     ESTIMATED ERROR                                     *
*                                                                      *
************************************************************************
*                                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,                    
     +           FIVE=5.0D0,SIX=6.0D0,TEN=10.0D0,TWELVE=12.0D0,                 
     +           SIXTY=60.0D0,HUNTWE=120.0D0)
      PARAMETER (IREAD=5,IWRITE=6)
      DIMENSION X(N),Y(N)
*
      RESULT = ZERO
      ERROR = ZERO
*
*  CHECK THAT WE HAVE ENOUGH POINTS
*
      IF (N.LT.4) THEN
       WRITE(IWRITE,1000) N
       STOP
      ENDIF
*
*  CHECK THAT THE MESHPOINTS AR IN EITHER ASCENDING OR DESCENDING ORDER
*
      H2 = X(2)-X(1)
      DO 10 I=3,N
       H3 = X(I)-X(I-1)
       IF(H2*H3.LE.ZERO) THEN
        WRITE(IWRITE,1001)
        STOP
       ENDIF
10    CONTINUE
*
*  START THE INTEGRATION
*
      D3 = (Y(2)-Y(1))/H2
      H3 = X(3)-X(2)
      D1 = (Y(3)-Y(2))/H3
      H1 = H2+H3
      D2 = (D1-D3)/H1
      H4 = X(4)-X(3)
      R1 = (Y(4)-Y(3))/H4
      R2 = (R1-D1)/(H4+H3)
      H1 = H1+H4
      R3 = (R2-D2)/H1
      RESULT = H2*(Y(1)+H2*(D3/TWO-H2*(D2/SIX-(H2+TWO*H3)*R3/TWELVE)))
      S = -(H2**3)*(H2*(THREE*H2+FIVE*H4)+TEN*H3*H1)/SIXTY
      R4 = ZERO
      NN = N-1
*
*  LOOP OVER POINTS 2 TO N-1
*
      DO 20 I=3,NN
       RESULT = RESULT+H3*((Y(I)+Y(I-1))/TWO-H3*H3*(D2+R2+(H2-H4)*R3) 
     +         /TWELVE)
       C = H3**3*(TWO*H3*H3+FIVE*(H3*(H4+H2)+TWO*H4*H2))/HUNTWE
       ERROR = ERROR+(C+S)*R4
       IF (I.NE.3) THEN
        S = C
       ELSE
        S = S+TWO*C
       ENDIF
       IF (I.EQ.(N-1)) GOTO 30
       H1 = H2
       H2 = H3
       H3 = H4
       D1 = R1
       D2 = R2
       D3 = R3
       H4 = X(I+2)-X(I+1)
       R1 = (Y(I+2)-Y(I+1))/H4
       R4 = H4+H3
       R2 = (R1-D1)/R4
       R4 = R4+H2
       R3 = (R2-D2)/R4
       R4 = R4+H1
       R4 = (R3-D3)/R4
20    CONTINUE
30    CONTINUE
*
*  FINISH INTEGRATION
*
      RESULT = RESULT+H4*(Y(N)-H4*(R1/TWO+H4*(R2/SIX+(TWO*H3+H4)*R3             
     +        /TWELVE)))
      ERROR = ERROR-H4**3*R4*(H4*(THREE*H4+FIVE*H2)+TEN*H3*(H2+H3+H4)) 
     +        /SIXTY+S*R4
*
1000  FORMAT(//,' ERROR IN SUBROUTINE INTEGR . N = ',I2,3X,                     
     +      'PROGRAM STOPS.',//)                                                
1001  FORMAT(//,' ERROR IN SUBROUTINE INTEGR . MESHPOINTS ARE OUT ',            
     +      'OF ORDER. PROGRAM STOPS.',//)                                      
      RETURN
      END

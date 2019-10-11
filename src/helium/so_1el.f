c--------------------------------------------------------------------------
c     v = u/r;   v' = (u' - u/r) / r ; where v(:) = vdcore(:)

      subroutine  one_el_V_DERIV(vdcore)
      use  VDCORE_MODULE
      include 'par.f'
      real vdcore(maxr,0:lamax)
      double precision  funf(nr)

      common /meshrr/ nr, gridr(nmaxr,3)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      real Z_fit(0:lamax)      
      real v_exch(nr,0:lamax)   ! local equiv. exchange potential of the innert  core
      real rho(nr)      ! (electron density  of the innert  core) * 4 * pi
      real v_core(nr,0:lamax)
      real temp(nr)
      real fit_par(0:lamax)
      double precision sec_der(nr)
      real fun(nr)
      real  func_fit(nr,0:lamax), fit_exp(0:lamax), fit_coef(0:lamax)
      real A(0:lamax), d(0:lamax)
      double precision, dimension(:,:), allocatable :: Hso
c     
c     
      allocate(vd_deriv(nr,0:lamax))
      vd_deriv = 0.0
      allocate(Hso(nnmax,nnmax))
      Hso = 0.0
c
c      
c     ze = 1d0       ! this is charge of the core (for neutral 1-electron targets)
c      ze = 2d0      ! this is charge of the core (for A^{+} 1-electron targets)
      ze = zasym
      
      pi = acos(-1.)

      
c     define central direct core potential v_core:
      do l = labot, latop
         v_core(1:nr,l) = - ze/gridr(1:nr,1)
      end do
      v_core(1:nr,labot:latop) = v_core(1:nr,labot:latop) +
     >   vdcore(1:nr,labot:latop) 
      
c$$$c     calculate electron density in the core  (* 4 * pi)
c$$$      rho = 0.0      
c$$$      do l = 0, lactop
c$$$         const = real(2*(2*l+1))
c$$$         do n=l+1,nabot(l)-1
c$$$c            print*, n, l, gridr(1,1), psinb(1,n,l),psinb(1,n,l)/gridr(1,1)
c$$$            do i=1,istoppsinb(n,l)
c$$$               rho(i) = rho(i) + const*psinb(i,n,l)*psinb(i,n,l)/gridr(i,1)/gridr(i,1)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$
c$$$c     calculate equivalent local potential for the exchange.
c$$$      const = 6.0*(3.0/8.0/pi)**(1.0/3.0)
c$$$      fit_par = 0.0
c$$$      fit_par(0) = 1.0/17. !1.4196 
c$$$      fit_par(1) = 1.0/18. ! 1.36 
c$$$      fit_par(2) = 1.0/12. ! 1.48 
c$$$      fit_par(3) = 1.0/13. ! 1.09 
c$$$      v_exch = 0.0      
c$$$      do l = labot, latop
c$$$c         v_exch(:,l) = -fit_par(l)*(sqrt(v_core(:,l)*v_core(:,l) + rho(:)) + v_core(:,l))/2.0
c$$$         v_exch(:,l) = -fit_par(l)*const*(rho(:))**(1.0/3.0)
c$$$      end do
c$$$
c$$$      
c$$$c     Define total central local core potential:
c$$$c      v_core(1:nr,labot:latop) = v_core(1:nr,labot:latop) + v_exch(1:nr,labot:latop)
c$$$
c     

c     This is transformation to have more accurate numerical differentiation
      do l=0,latop
         v_core(1:nr,l) = v_core(1:nr,l) * gridr(1:nr,1)
      end do

      do l=1,latop
         call f_deriv(v_core(1:nr,l), 1, nr, gridr(1:nr,1),
     >      nr, 1, vd_deriv(1:nr,l))
      end do
      
      do l=1,latop
         vd_deriv(1:nr,l) = (vd_deriv(1:nr,l) -
     >      v_core(1:nr,l)/gridr(1:nr,1)) / gridr(1:nr,1)
      end do

c     restore back v_core() array
      do l=0,latop
         v_core(1:nr,l) = v_core(1:nr,l) / gridr(1:nr,1)
      end do
c     include  (1/2r)  factor
      do l=1,latop
         vd_deriv(1:nr,l) = vd_deriv(1:nr,l) / (2.0 * gridr(1:nr,1))
      end do

      Z_fit(:) = 1.0
      if(nznuc .eq. 80)  then  ! mercury
         print*, 'mercury fitting for spin-orbit term'
         Z_fit(0) = 0.0
         Z_fit(1) = 0.8  ! 2.4 ! 2.
         Z_fit(2) = 1.42  ! 1.9 ! 1.6
         Z_fit(3) = -10.  ! -3.0 !  -1.7
      else if(nznuc .eq. 56)  then ! barium
         print*, 'barium fitting for spin-orbit term'
         Z_fit(0) = 0.0
         Z_fit(1) = 1.465   ! 1.38
         Z_fit(2) =  0.93 ! 0.874
         Z_fit(3) = 0.95     ! 0.9
c         d(:) = 6.0
c         A(:) = 200.0
      else if(nznuc .eq. 70)  then ! Yb
         print*, 'barium fitting for spin-orbit term'
         Z_fit(0) = 0.0
         Z_fit(1) = 1.15
      end if
c      A(0) = 0.0  ! just to make sure that for l=0 states nothing is added.
      
      do l=1,latop      
         vd_deriv(1:nr,l) = Z_fit(l) * vd_deriv(1:nr,l)
c     >      + A(l)*(1.0 - exp(-(gridr(1:nr,1)/d(l))**6) )/gridr(1:nr,1)**2
      end do
            
      where(abs(vd_deriv) .lt. expcut)
         vd_deriv = 0.0
      end where

      
c$$$      go to 20
c$$$      Z_fit = 0.0
c$$$      Z_fit(1) = 72.5
c$$$      Z_fit(2) =  38.063 ! 41.8 ! 38.063
c$$$      Z_fit(3) = 22.74
c$$$      do l=1,latop
c$$$         vd_deriv(:,l) = Z_fit(l)/gridr(:,1)**3/2.0
c$$$      end do
c$$$ 20   continue
           
c     Get fine structure splitting for A^{+} (positive ion).      
      alpha = 1.0/137.04     ! fine structure constant
            
      do l = 0, latop ! max(1,labot), latop
         do k = nabot(l), nabot(l) + 3 !natop(l)
            m1 = 1
            m2 = istoppsinb(k,l)

            fun(1:nr) = 0.0
            fun(m1:m2) = psinb(m1:m2,k,l) * psinb(m1:m2,k,l) *
     >         gridr(m1:m2,3)
            funf(m1:m2) = fun(m1:m2) * vd_deriv(m1:m2,l)
            
            split = (l + 0.5) * SUM(funf(m1:m2))*alpha*alpha
            
            print*, 'k =',k,', l =', l,
     >         ', splitting =', split, '(au) ',split*27.2116, '(eV) ',
     >         split*27.2116 * 8068.7, '(cm^{-1})) ',' ; <r> =',
     >         SUM(fun(m1:m2)*gridr(m1:m2,1)), ' ; <r**3> =',
     >         SUM(fun(m1:m2)*gridr(m1:m2,1)**3)

         end do
      end do

c$$$      open(137,file='tmp.sp', recl=10000)
c$$$      write(137,'("#",15X,100(I8,8X))')
c$$$     >   ((k, k=nabot(l),nabot(l)+2), l=0,latop)
c$$$      do i=1,nr,3
c$$$         write(137,'(1p,100E16.8)')  gridr(i,1),((psinb(i,k,l),
c$$$     >      k=nabot(l), nabot(l) + 2), l = 0, latop)
c$$$      end do
c$$$      close(137)

c$$$      open(139,file='tmp.sp9', recl=10000)
c$$$      write(139,'("#",15X,100(I8,8X))')
c$$$     >   ((k, k=natop(l)-2,natop(l)), l=0,latop)
c$$$      do i=1,1100,2
c$$$         write(139,'(1p,100E16.8)')  gridr(i,1),((psinb(i,k,l),
c$$$     >       k=natop(l)-2,natop(l)), l = 0, latop)
c$$$      end do
c$$$      close(139)
      
c      open(137,file='tmp.pot')
c      do i=1,800
c         write(137,'(100E12.4)')  gridr(i,1), gridr(i,3),
c     >      (vd_deriv(i,l), l = 1, latop)
c      end do
c      close(137)


            
      do l = max(1,labot), latop
         do kp = nabot(l), natop(l)
            m1p = 1
            m2p = istoppsinb(kp,l)
            do k = kp, natop(l)
               m1 = 1
               m2 = istoppsinb(k,l)
               
               i1=max(m1p,m1)
               i2=min(m2p,m2)
               
c     fun(:) = 0.0
c     fun(i1:i2) = psinb(i1:i2,kp,l) * psinb(i1:i2,k,l) * gridr(i1:i2,3)
c     funf(i1:i2) = fun(i1:i2) * vd_deriv(i1:i2,l)
               
c     Hso(kp,k) = (l + 0.5) * SUM(funf(i1:i2))*alpha*alpha
            
            end do
         end do
      end do


      deallocate(Hso)
      
      return
      end  subroutine  one_el_V_DERIV

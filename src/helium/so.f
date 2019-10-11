c     array vd_deriv(nr,lamax) is set in subroutine one_el_V_DERIV, file so_1el.f
c      
c-------------------------------------------------------------------
      subroutine F90_make_new_C(C_old,Nmax,nspmW,namax)
      use CI_MODULE             ! module with CI-array C
      include 'par.f' 

      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      double precision  C_old(Nmax+1,nspmW,nspmW)

      allocate(C(Nmax,namax,namax))

      do N=1,Nmax
         do j=1,nam(N)
            do i=1,nam(N)
               C(N,j,i) =  C_old(N,j,i)
c               print'(4i3,2e15.4)', N,  nam(N), j, i, real(C(N,j,i)), real(C_old(N,j,i))
            end do      
         end do      
      end do      

      return
      end
c-------------------------------------------------------------------
c     spin-orbit term in the structure code.

      subroutine Spin_orbit(Nmax,Energy,enionry,i_stmix,i_ng)
c     Nmax is the number of states that comes from the norelativistic calculation
c
      use STATE_MODULE           ! module with definitions of the state type
      use CI_MODULE              ! module with CI-array C
      use  VDCORE_MODULE         ! module with array vd_deriv(nr,lamax)
      include 'par.f'
      double precision Energy(KNM)   ! state energy in a.au.
      character(LEN=3) chan(knm) ! state label 
      common /charchan/ chan     ! 
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer N, jn1,jn2, n1,n2, J2min,J2max,J2, ipa
      type(state), dimension(Nmax) :: States ! atomatic array of state structures
      type(configuration), pointer :: config
      
      integer:: i_sw, ncj
      integer:: nc_all, imix_max
      real, dimension(:,:), allocatable:: rmix_ci
      integer, dimension(:,:), allocatable:: imix_n
      type(state), dimension(:), allocatable :: StLSJ ! array of state structures: in intermediate coupling
      real, dimension(:), allocatable:: EnValues
      integer, dimension(:), allocatable:: Num
c
c
      if(i_ng .eq. 1) then
         print*, 'Start spin-orbit correction for noble gases'
      else
         print*, 'Start spin-orbit correction for 2-electron atoms'
      endif
c      
c     open file where mixing coef. are saved
      open(138,file='mixing_ordered.so',RECL=10000)
      write(138,'("# Ordered: Mixing coef. of nonrelativistic states",
     >   " due to the spin-orbit correction: ",
     >   " J, parity, rank_of_mixing_matrix")')

c      Find maximum value of orbital ang.mom. for s.p. orbitals
      lmax_orb = -1
      do nsp=1,nspmax
          lmax_orb = max(lmax_orb, lo(nsp))
      end do
      
      do N=1,Nmax
         States(N)%en = Energy(N)
         States(N)%label = chan(N)
         States(N)%pa =  lparity(N)
         States(N)%la =  ll(N)
         States(N)%sa =  ls(N)
         States(N)%ja =  -1.0 ! is not defined yet
         States(N)%nst =  N
c         print'(3i5,f8.2,A5,e12.5,A4,1P,e12.5,A4)',N,States(N)%pa,
c     >        States(N)%la,States(N)%sa,States(N)%label, 
c     >        real(States(N)%en),'(au)',
c     >        27.2116 * (real(States(N)%en)-enionry/2.0),'(eV)'


         inum = 0 
         
         allocate(States(N)%first_config)
         config => States(N)%first_config
         nullify (config%NEXT)
         do jn1=1,nam(N)
            n1 = na(jn1,N)
            do jn2=1,nam(N)
               n2 = na(jn2,N)
c               print'(i3,A5,3i3,e15.4)', N, States(N)%label, 
c     >              jn1, jn2, nam(N), real(C(N,jn1,jn2))
               if(C(N,jn1,jn2).eq.0d0) cycle             
               inum = inum + 1
               config%n1 = n1
               config%n2 = n2
               config%l1 = lo(n1)
               config%l2 = lo(n2)
               config%value = C(N,jn1,jn2)
               allocate(config%NEXT)
               config => config%NEXT
               nullify(config%NEXT)
            end do
         end do
         config => States(N)%first_config
         do
c            print*, associated(config),associated(config%NEXT), associated(config%NEXT%NEXT), config%n1,  config%n2
            if( .not. associated(config%NEXT%NEXT) ) then
c               print*,'removing list element'
               deallocate(config%NEXT)
               exit
            end if
            config => config%NEXT
         end do                  
         States(N)%num_config = inum
      end do                   !end N loop

c      do N_i=1,Nmax
c         do N_f=1,N_i
c            print'(2i5,e15.5)', N_f, N_i, real(st_overlap(States(N_f), States(N_i)))
c         end do
c      end do
      
c     Find maximum and minium possible values of J
      J2max = -1
      J2min = 10000
      do N=1,Nmax
         J2max = max(J2max, nint(2.0*(States(N)%la + States(N)%sa)))
         J2min=min(J2min,nint(2.0*(abs(States(N)%la - States(N)%sa))))
      end do

      print*,'J_min =', real(J2min)/2.0, ', J_max =', real(J2max)/2.0
c     Peform calculations for given J = nint( J2 / 2.0 ) and parity ipa
      nc_all = 0
      imix_max = 0
      allocate(imix_n(1,1),rmix_ci(1,1),STLSJ(1)) ! temporary allocation to avoid problems with passing of these variabls in subroutine Set_up_J()
      i_sw = 1   !  first find total number of states
      do J2=J2min, J2max,2
         do ipa=-1,1,2
            call Set_up_J(Nmax,States,J2,ipa,enionry,i_sw,ncj,
     >           nc_all,imix_max,imix_n,rmix_ci,StLSJ,Nst,i_ng)
            nc_all = nc_all + ncj
            imix_max = max(imix_max,ncj)
         end do
      end do
      i_sw = 0
      print*
      print*,'nc_all =', nc_all, ', imix_max=',  imix_max
      if(allocated(imix_n)) then
         deallocate(imix_n,rmix_ci,STLSJ)
      endif
      allocate(imix_n(nc_all,imix_max),rmix_ci(nc_all,imix_max))
      imix_n(:,:) = 0
      rmix_ci(:,:) = 0.0
      allocate(StLSJ(nc_all))
      allocate(EnValues(nc_all))
      allocate(Num(nc_all))

      Nst = 0
      do J2=J2min, J2max,2
         do ipa=-1,1,2
            call Set_up_J(Nmax,States,J2,ipa,enionry,i_sw,ncj,
     >           nc_all,imix_max,imix_n,rmix_ci,StLSJ,Nst,i_ng)
         end do
      end do

      EnValues(1:nc_all) = StLSJ(1:nc_all)%en
      do i=1,nc_all
         Num(i) = i
      enddo
      call Sort_subr(EnValues,Num,nc_all)
      do j=1,imix_max
         imix_n(1:nc_all,j) =  imix_n(Num(1:nc_all),j)
         rmix_ci(1:nc_all,j) = rmix_ci(Num(1:nc_all),j)
      enddo
      StLSJ(1:nc_all) =  StLSJ(Num(1:nc_all))


      imaxtmp = min(imix_max,5)
      print*,'   N label     E(N)       exit.en.  j    par  ',
     >     '( N-nonrel mix(N-nonrel) )' 
      do N=1,nc_all
         ncm = StLSJ(N)%nst
         ex_en =  (StLSJ(N)%en - StLSJ(1)%en) * 27.2116
         write(*,'(I5,1X,A3,1X,2F12.5,F5.1,I5,100(I5,F12.5))') 
     >        N, StLSJ(N)%label, StLSJ(N)%en, ex_en, StLSJ(N)%ja, 
     >        StLSJ(N)%pa, (imix_n(N,j),rmix_ci(N,j), j=1,imaxtmp)
         
      enddo


      if(i_stmix .eq. 2) then
         call osc_intcp(nc_all,imix_max,imix_n,rmix_ci,StLSJ,i_ng)
      endif

c     deallocate configuration structures.
      do N=1,Nmax
         config => States(N)%first_config
         do                     
            if( .not. associated(config%NEXT) ) then
               deallocate(config)
               exit
            end if
            config => config%NEXT
         end do
      end do

c     deallocate array vd_deriv(nr,lamax)
      deallocate(vd_deriv)
c
      close(138)

      deallocate(imix_n,rmix_ci)

c     stop

      return
      end  subroutine Spin_orbit

c--------------------------------------------------------------------------
      subroutine Set_up_J(Nmax,States,J2,ipa,enionry,i_sw,ncj,
     >     nc_all,imix_max,imix_n,rmix_ci,StLSJ,Nst,i_ng)
      use STATE_MODULE          ! module with stdefinitions of the state type
      include 'par.f'

      integer:: Nmax
      type(state), dimension(Nmax) :: States ! array of state structures
      integer:: J2, ipa
      real:: enionry
      integer:: i_sw, ncj
      integer:: nc_all, imix_max
      real, dimension(nc_all,imix_max):: rmix_ci
      integer, dimension(nc_all,imix_max):: imix_n
      type(state), dimension(nc_all):: StLSJ ! array of state structures: in intermediate coupling      
      integer:: Nst

      double precision H12_so, H12_so_ng, Mass_cor! function
      integer nco, ncm, nc1, nc2      
      integer, dimension(:), allocatable ::  N_nc
      double precision,  dimension(:,:), allocatable :: Ham
      double precision,  dimension(:), allocatable :: energy
      double precision alpha
      integer, dimension(:), allocatable :: Num
      double precision tmp, Hsotmp
      integer ik, i_max
c     integer i_max(1)
      
      alpha = 1d0/137.04d0
      
c     find number of configurations  for given J and parity
      nco = 0
      do N=1,Nmax
c     print*, N, States(N)%la, States(N)%sa, States(N)%pa
         if(J2 .le. nint(2.0*(States(N)%la + States(N)%sa)) .and.
     >        J2 .ge. nint(2.0*abs(States(N)%la - States(N)%sa)) .and.
     >        ipa .eq. States(N)%pa) then
            nco = nco + 1
         end if
      end do
      ncm = nco                 ! this is the number of configurations used to diagonalize spin-orbit Hamiltonian
      ncj = ncm
      if(ncm .eq. 0) then
         return                 ! this combination of J and parity is not present
      else
         print*
         write(*,'(" Calculating J = ", F4.1,", parity =",
     >        I2,", ncm =",I2)') J2/2.0, ipa, ncm
      end if
      if(i_sw .eq. 1) then
         return
      endif

      allocate(N_nc(ncm))
      allocate(energy(ncm))
      allocate(Ham(ncm,ncm))
      allocate(Num(ncm))
      N_nc = 0
      Ham = 0d0
      
      nco = 0
      do N=1,Nmax
         if(J2 .le. nint(2.0*(States(N)%la + States(N)%sa)) .and.
     >      J2 .ge. nint(2.0*abs(States(N)%la - States(N)%sa)) .and.
     >      ipa .eq. States(N)%pa) then
            nco = nco + 1
            N_nc(nco) = N
            States(N)%ja = real(J2)/2.0
         end if
      end do
      
c     Set diagonal matrix elements - energies of the nonrelativistic states.
      do nco=1,ncm
         Ham(nco,nco) =   States(N_nc(nco))%en
      end do
c     Calculate matrix elements for spin-orbit term. 
      do nc1=1,ncm
c     do nc2=nc1,nc1
         do nc2=1,nc1
            if(i_ng .eq. 1) then
               Hsotmp=H12_so_ng(States(N_nc(nc1)),
     >              States(N_nc(nc2)),J2)
            else
               Hsotmp=H12_so(States(N_nc(nc1)),States(N_nc(nc2)),J2)
            endif

            Ham(nc1,nc2) = Ham(nc1,nc2) +
     >         ( Hsotmp 
c     >       +  Mass_cor(States(N_nc(nc1)),States(N_nc(nc2))) 
     >         ) * alpha*alpha
            Ham(nc2,nc1) = Ham(nc1,nc2)
         end do
      end do
      

      call diagonalize_Ham(ncm, Ham, energy) 
      
c     write(*,'("eigenvectors")')
      write(138,'(F6.1,2I5," J, parity, number of config.")')
     >   real(J2)/2.0, ipa, ncm
      
      do i=1,ncm
 
c     write(138,'(F4.1,1X,I2,1X,A3,1X,I3,1X,2F12.6," : ",
c     >      1P,1000(e12.3,A4))')
c     >      real(J2)/2.0, ipa, States(i_max)%label, ncm, 
c     >      27.2116*(real(energy(i))-enionry/2.0), real(energy(i)),
c     >      (Ham(j,i),States(N_nc(j))%label, j=1,ncm)
      
c     set up Num() array - to be used in sorting CI coef. by their values by subroutine Sort_subr().
         do k=1,ncm
            Num(k) = k
         end do
c     Sort by the value of largest CI coef, note that the largest value CI coef. is the last
c     after sorting, hence inverse loop: j=ncm,1,-1
         call Sort_subr(abs(real(Ham(1:ncm,i))),Num,ncm)
         i_max = N_nc(Num(ncm))
c     write ordered CI coef.
         write(138,'(F4.1,1X,I2,1X,A3,1X,I3,1X,2F12.6," : ",
     >        1P,1000(e12.3,A4))')
     >        real(J2)/2.0, ipa, States(i_max)%label, ncm, 
     >        27.2116*(real(energy(i))-enionry/2.0), real(energy(i)),
     >        (Ham(Num(j),i),States(N_nc(Num(j)))%label, j=ncm,1,-1)
         
         
         Nst = Nst + 1
         StLSJ(Nst)%label = States(i_max)%label
         StLSJ(Nst)%en = real(energy(i))-enionry/2.0
         StLSJ(Nst)%pa = ipa
         StLSJ(Nst)%ja =  real(J2)/2.0
         StLSJ(Nst)%nst =  Nst
         StLSJ(Nst)%num_config =  ncm
         StLSJ(Nst)%la = States(i_max)%la
         StLSJ(Nst)%sa = States(i_max)%sa
c     StLSJ(Nst)%first_config = StLSJ(i_max)%first_config 
         
         do j=ncm,1,-1
            jj = ncm - j + 1
            imix_n(Nst,jj) = N_nc(Num(j))
            rmix_ci(Nst,jj) = Ham(Num(j),i)            
         enddo
      end do
      

      deallocate(N_nc)
      deallocate(energy)
      deallocate(Ham)
      deallocate(Num)
      
      return
      end  subroutine Set_up_J
c--------------------------------------------------------------------------
c     This is  M.E. of spin-orbit term. This is NOT reduced ME: <(f)J| H | J(i)>
c     Summation over all valence electrons is accounted here.
      double precision function H12_so(State_f, State_i, J2)
      use STATE_MODULE          ! module with definitions of the state type
      use  VDCORE_MODULE  ! array vd_deriv(nr,lamax) is set in subroutine one_el_V_DERIV, file so_1el.f
      include 'par.f'
      common /ortog/  ortint(nspmax,nspmax)

      integer J2
      type(state), intent (in) :: State_f, State_i !  state structures

      double precision ortint
      real rJ, sa_i, sa_f
      integer l1_i, l2_i, l1_f, l2_f
      double precision ang, tmp_hat, value_i, value_f
      double precision  Radint, Rel_Mass_cor ! function
      type(configuration), pointer :: config_i, config_f
      double precision  hat, hatr
      real rl
      integer l
      integer number_of_electrons
      data number_of_electrons /2/ ! set the number_of_electrons
c
      
      hat(l) = dsqrt(dble(2 * l + 1))
      hatr(rl) = dsqrt(2.0 * rl + 1d0)
      
      H12_so = 0d0
      
      sa_i = State_i%sa
      sa_f = State_f%sa
      if(sa_f .eq. 0.0 .and. sa_i .eq. 0.0) return ! H12_so is zero in this case
      la_i = State_i%la
      la_f = State_f%la
      if(la_f .eq. 0 .and. la_i .eq. 0) return ! H12_so is zero in this case
      rJ = State_i%ja
      if(nint(2.0*rJ) .ne. J2) then
         print*, 'Wrong value of J:', rJ, real(J2)/2.0
         stop
      end if

      tmp_hat =  -dsqrt(1.5d0) * (-1)**(nint(rJ + sa_i + sa_f)) *
c     >   hatr(rJ)*
     >     hat(la_i)*hat(la_f)*hatr(sa_i)*hatr(sa_f) *
     >   COF6J(0.5,0.5,real(sa_f),1.0,real(sa_i),0.5) *
     >   COF6J(real(la_f),sa_f,rJ,sa_i,real(la_i),1.0) *
     >   State_i%pa
c     Note: term  State_i%pa is for  (-1)**(l1_i + l2_i) which give parity of the initial state .
c     print*, 'tmp_hat =', tmp_hat
      if(tmp_hat .eq. 0d0) return      

      config_f => State_f%first_config
      do 
c         print*,  associated(config_f), associated(config_f%NEXT),
c     >      config_f%n1,config_f%n2, config_f%value
         if( .not. associated(config_f) ) exit         
         l1_f = config_f%l1
         if(l1_f .gt. 0) then  ! H12_so is zero in this case, otherwise we have:
            l2_f = config_f%l2
            n1_f = config_f%n1
            n2_f = config_f%n2
            value_f = config_f%value

            config_i => State_i%first_config
            do
               if( .not. associated(config_i) ) exit            
               n1_i = config_i%n1
               n2_i = config_i%n2
               l1_i = config_i%l1
               l2_i = config_i%l2
               value_i = config_i%value
               
               if(l1_f .eq. l1_i .and. l2_f .eq. l2_i) then
                  ang =  tmp_hat * ortint(n2_f,n2_i) * 
     >               dsqrt(dble(l1_i*(l1_i + 1)*(2*l1_i + 1))) * 
     >               COF6J(real(l2_i),real(l1_i),real(la_f),1.0,
     >               real(la_i),real(l1_i)) 
c     >               / dsqrt(J2 + 1d0)
                  if(ang .ne. 0d0) then
                     H12_so = H12_so + ang * value_f * value_i *
     >                  Radint(n1_f,n1_i,vd_deriv(1,l1_i)) 
                     
c                     print'(6I5,2e12.4))', State_f%nst, State_i%nst!, l1_f,l2_f,l1_i,l2_i, real(ang), real(H12_so)
c                  print'(8i4,1P,3e15.5)', n1_f,n1_i, n2_f,n2_i, l1_f,l1_i, l2_f,l2_i,
c     >                  real(ang), real(Radint(n1_f,n1_i,vd_deriv(:,l1_i)))
                  end if
               end if
               config_i => config_i%NEXT
            end do
         end if
         config_f => config_f%NEXT         
      end do
      H12_so = H12_so * dble(number_of_electrons) ! multiply by the number of valence electrons.
c      print*, real(H12_so)
      return
      end  function H12_so
c--------------------------------------------------------------------------
      double precision function  Radint(n1_f,n1_i,vd_deriv)
      include 'par.f'
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common /meshrr/ nr, gridr(nmaxr,3)
      double precision  vd_deriv(nr)
      double precision, dimension(nr) :: fun 
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
            
      Radint = 0d0

      fun(:) = 0d0
      
      minfun=max(minf(n1_f),minf(n1_i))
      maxfun=min(maxf(n1_f),maxf(n1_i))
      
      fun(minfun:maxfun) = dble(fl(minfun:maxfun,n1_f)) *
     >   dble(fl(minfun:maxfun,n1_i)) *
     >   dble(gridr(minfun:maxfun,3)) * vd_deriv(minfun:maxfun)
      

      Radint = SUM(fun(minfun:maxfun))

c       print*, 'Radint =', Radint

      return
      end  function   Radint
c--------------------------------------------------------------------------
c     Mass correction.  Summation over all valence electrons is accounted here.
      double precision function Mass_cor(State_f, State_i)
      use STATE_MODULE          ! module with definitions of the state type
      include 'par.f'
      type(state), intent (in) :: State_f, State_i !  state structures
      common /ortog/  ortint(nspmax,nspmax)
      double precision ortint
      real sa_i, sa_f
      integer l1_i, l2_i, l1_f, l2_f
      double precision value_i, value_f
      double precision  Rel_Mass_cor ! function
      type(configuration), pointer :: config_i, config_f
      integer number_of_electrons
      data number_of_electrons /2/ ! set the number_of_electrons
c     
      Mass_cor = 0d0
      
      sa_i = State_i%sa
      sa_f = State_f%sa
      la_i = State_i%la
      la_f = State_f%la
      if(sa_f .ne. sa_i .or. la_f .ne. la_i) return
      
      config_f => State_f%first_config
      do 
         if( .not. associated(config_f) ) exit         
         l1_f = config_f%l1
         l2_f = config_f%l2
         n1_f = config_f%n1
         n2_f = config_f%n2
         value_f = config_f%value
         
         config_i => State_i%first_config
         do
            if( .not. associated(config_i) ) exit            
            n1_i = config_i%n1
            n2_i = config_i%n2
            l1_i = config_i%l1
            l2_i = config_i%l2

            value_i = config_i%value
            if(l1_f .eq. l1_i .and. l2_f .eq. l2_i .and.
     >         ortint(n2_f,n2_i) .ne. 0d0) then
               Mass_cor = Mass_cor - ortint(n2_f,n2_i) *
     >            Rel_Mass_cor(n1_f,n1_i,l1_f) *
     >            value_f * value_i / 8d0
            end if
         
            config_i => config_i%NEXT
         end do
         config_f => config_f%NEXT         
      end do
      
      Mass_cor = Mass_cor * dble(number_of_electrons) ! multiply by the number of valence electrons.

      return
      end  function Mass_cor
c--------------------------------------------------------------------------
c     Mass term: (p^2)^2, Darwin term Z delta(r)
      double precision function  Rel_Mass_cor(n1_f,n1_i,lsp)
      include 'par.f'
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common /meshrr/ nr, gridr(nmaxr,3)
      double precision, dimension(nr) :: fun_f, fun_i, fun
      double precision Z
      common /Zatom/ Z

      fun = 0d0; fun_f = 0d0; fun_i = 0d0
      
      minfun=max(minf(n1_f),minf(n1_i))
      maxfun=min(maxf(n1_f),maxf(n1_i))

c     Mass term
      
      call f_deriv(fl(:,n1_f), minfun, maxfun, gridr(:,1),nr,2,fun_f)
      call f_deriv(fl(:,n1_i), minfun, maxfun, gridr(:,1),nr,2,fun_i)
      if(lsp .eq. 0) then
         fun_f(minfun:maxfun) = fun_f(minfun:maxfun) 
         fun_i(minfun:maxfun) = fun_i(minfun:maxfun) 
      else
         fun_f(minfun:maxfun) = fun_f(minfun:maxfun) -
     >      dble(lsp*(lsp + 1)) *
     >      dble(fl(minfun:maxfun,n1_f)/gridr(minfun:maxfun,1)/
     >      gridr(minfun:maxfun,1))
         fun_i(minfun:maxfun) = fun_i(minfun:maxfun) -
     >      dble(lsp*(lsp + 1)) *
     >      dble(fl(minfun:maxfun,n1_i)/gridr(minfun:maxfun,1)/
     >      gridr(minfun:maxfun,1))
      end if
      fun(minfun:maxfun) = fun_f(minfun:maxfun)*fun_i(minfun:maxfun)*
     >   dble(gridr(minfun:maxfun,3)) 

      Rel_Mass_cor = SUM(fun)

c     Darwin term
      tmp = 0.0
      if(lsp .eq. 0) then
         tmp = Z*(fl(1,n1_f)/gridr(1,1))*(fl(1,n1_i)/gridr(1,1))
      end if
      print*, n1_f, n1_i, lsp, real(Rel_Mass_cor), tmp
c
      Rel_Mass_cor =  Rel_Mass_cor  - tmp
      
      return
      end  function  Rel_Mass_cor
c--------------------------------------------------------------------------
c$$$c     array vd_deriv(nr,lamax) is set in subroutine one_el_V_DERIV
c$$$c     no need to call this routine now. 
c$$$      subroutine  set_VD_DERIV(vdcore, lmax_orb)
c$$$      use  VDCORE_MODULE
c$$$      include 'par.f'
c$$$      
c$$$      common /so_Z_eff/ Z_fit(1:lamax)
c$$$
c$$$      real vdcore(maxr,0:lamax)
c$$$      common/smallr/ formcut,regcut,expcut,fast
c$$$      logical fast
c$$$      common /meshrr/ nr, gridr(nmaxr,3)
c$$$      double precision  h, h1, h2, fun(nr,0:4), funf(nr), fun2(nr)
c$$$      common /funLag/  fl(nmaxr,nspmax)
c$$$      common /minmaxf/ maxf(nspmax),minf(nspmax)
c$$$      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
c$$$      common/hame1/ e1r(nspmax,nspmax)
c$$$      double precision  e1r
c$$$      double precision Z
c$$$      common /Zatom/ Z
c$$$      
c$$$
c$$$      alpha = 1d0/137.04d0
c$$$      do nsp=1,-nspm
c$$$         if(ko(nsp) .gt. 2) cycle
c$$$         m1 = minf(nsp)
c$$$         m2 = maxf(nsp)
c$$$         l = lo(nsp)
c$$$         funf(m1:m2) = fl(m1:m2,nsp)*fl(m1:m2,nsp)*gridr(m1:m2,3)
c$$$c     Mass correction
c$$$         fun(m1:m2,l) = dble(vdcore(m1:m2,l))  - ze/dble(gridr(m1:m2,1))
c$$$         result1 = SUM(fun(m1:m2,l)*funf(m1:m2))
c$$$         result2 = SUM(fun(m1:m2,l)*fun(m1:m2,l)*funf(m1:m2))
c$$$         result = e1r(nsp,nsp)**2 - 2d0*e1r(nsp,nsp)*result1 + result2
c$$$c
c$$$         tmp_D = 0.0; tmp_D2 = 0.0
c$$$         if(l .eq. 0) then
c$$$            do nsp2= nsp,nsp    ! 1,nspm
c$$$               if(lo(nsp2) .eq. 0) then
c$$$                  tmp_D = real(Z)*(fl(1,nsp)/gridr(1,1))*(fl(1,nsp2)/gridr(1,1)) ! Darwin term
c$$$                  m1 = max(2,m1,minf(nsp2))
c$$$                  m2 = min(m2,maxf(nsp2))
c$$$                  tmp_D2 = SUM(fun2(m1:m2)*fl(m1:m2,nsp)*fl(m1:m2,nsp2)*gridr(m1:m2,3))
c$$$               end if
c$$$            end do
c$$$         end if
c$$$      end do
c$$$      
c$$$      return
c$$$      end  subroutine  set_VD_DERIV
c--------------------------------------------------------------------------
      subroutine diagonalize_Ham(ncm,Ham,energy)
      double precision Ham(ncm,ncm), energy(ncm)
      integer ncm
      double precision  z(ncm,ncm), fv1(ncm),fv2(ncm)
      integer matz, ierr
      integer i, j
      double precision  CI_sum

c      write(*,'("Hamiltonian")')
c      do i=1,ncm
c         write(*,'(5E15.6)') (real(Ham(i,j)), j=1,ncm)
c      end do

      matz = 1
      call rs(ncm,ncm,Ham,energy,matz,z,fv1,fv2,ierr)
      write(*,'("ierr =",I3)') ierr
      if(ierr .ne. 0) stop

c      write(*,'("eigenvalues in a.u.")')
c      write(*,'(5E15.6)') (real(energy(i)), i=1,ncm)
      
c      write(*,'("eigenvectors")')
c     Fix sign of eigenvectors.
      do j=1,ncm
         CI_sum = SUM(z(:,j))
         if(CI_sum . lt. 0d0) z(:,j) = -z(:,j)
      end do

c      do i=1,ncm
c         write(*,'(5E15.6)') (real(z(i,j)), j=1,ncm)
c      end do

      Ham = z
      
      
      return
      end  subroutine  diagonalize_Ham
c--------------------------------------------------------------------------
      subroutine f_deriv(f, minf, maxf, grid, nr, m, result)
      real, intent(in) :: f(nr)
      integer minf, maxf
      real, intent(in) ::  grid(nr)
      integer, intent(in) :: nr
      integer, intent(in) :: m
      double precision, intent(out) :: result(nr)
      double precision  h1, h2, h
      
      result(:) = 0d0

      if(m .eq. 1) then         ! First derivative
c         print*,' get  First derivative' 
         if(minf .eq. 1) then
            h = dble(grid(2)) - dble(grid(1))
c     result(1) = dble(f(2))/2d0 /h  ! note - value of fl(r,n)  at r=0.0 is zero.
            result(1) = (2d0 * dble(f(2)) -1.5d0 * dble(f(1)) -
     >         dble(f(3))/2d0 )/h
         end if
         
         do i=max(2,minf),min(nr,maxf-1)
            h2 = dble(grid(i+1) - grid(i))
            h1 = dble(grid(i) - grid(i-1))
            if(h2.eq.h1) then
               result(i) = (dble(f(i+1)) - dble(f(i-1)))/(h1+h2)
            else
               result(i) = (dble(f(i+1)) * (h1/h2) -
     >            dble(f(i-1)) * (h2/h1))/(h1+h2)
     >            + dble(f(i)) * ((h2-h1)/h1)/h2
            end if
         end do         
      end if
      
      if(m .eq. 2) then         ! Second derivative
c         print*,' get   Second derivative' 
         if(minf.eq.1) then 
            h = dble(grid(2)) - dble(grid(1))
            result(1) = (dble(f(2)) - 2d0*dble(f(1)))/h/h    ! note - value of fl(r,n)  at r=0.0 is zero.
         end if
         
         do i=max(2,minf),min(nr,maxf-1)
            h2 = dble(grid(i+1) - grid(i))
            h1 = dble(grid(i) - grid(i-1))
            if(h2.eq.h1) then
               result(i) = (dble(f(i+1)) - 2d0*dble(f(i)) +
     >            dble(f(i-1)))/h1/h1
            else
               result(i) = (dble(f(i+1)) / h2 + dble(f(i-1)) / h1) *
     >            2d0/(h1+h2) - dble(f(i)) * 2d0 / h1 / h2
            end if
         end do      
         

      end if      
      
      return
      end  subroutine  f_deriv
c-----------------------------------------------------------------------------------
c     Function state overlap - should be zero or one.
      function st_overlap(State_f, State_i)
      use STATE_MODULE          ! module with definitions of the state type
      include 'par.f'
      type(state), intent (in) :: State_f, State_i !  state structures
      common /ortog/  ortint(nspmax,nspmax)
      double precision ortint
      real sa_i, sa_f
      integer l1_i, l2_i, l1_f, l2_f
      double precision value_i, value_f, tmp
      type(configuration), pointer :: config_i, config_f
c     
      st_overlap = 0d0

      tmp = 0d0
      
      sa_i = State_i%sa
      sa_f = State_f%sa
      la_i = State_i%la
      la_f = State_f%la
      if(sa_f .ne. sa_i .or. la_f .ne. la_i) return
      
      config_f => State_f%first_config
      do 
         if( .not. associated(config_f) ) exit         
         l1_f = config_f%l1
         l2_f = config_f%l2
         n1_f = config_f%n1
         n2_f = config_f%n2
         value_f = config_f%value
         
         config_i => State_i%first_config
         do
            if( .not. associated(config_i) ) exit            
            n1_i = config_i%n1
            n2_i = config_i%n2
            l1_i = config_i%l1
            l2_i = config_i%l2
            value_i = config_i%value
            
            if(l1_f .eq. l1_i .and. l2_f .eq. l2_i .and.
     >         ortint(n2_f,n2_i) .ne. 0d0 .and.
     >         ortint(n1_f,n1_i) .ne. 0d0) then
               tmp  = tmp + ortint(n2_f,n2_i) * ortint(n1_f,n1_i) *
     >            value_f * value_i
            end if
         
            config_i => config_i%NEXT
         end do
         config_f => config_f%NEXT         
      end do
      
      st_overlap = real(tmp)

      return
      end  function st_overlap

c---------------------------------------------------------------------------------------
c     The index Num(N) to the array  Values(N) will be sorted by insertion sort algorithm.
      subroutine Sort_subr(Values,Num,N)
      integer, intent(in)             :: N
      integer, dimension(N)           :: Num
      real, intent(in), dimension(N)  :: Values
      integer i, j, itmp
      
      do i=2,N
         itmp = Num(i)
         j = i
         do
            if(j. lt. 2) exit
            if(Values(itmp) .gt. Values(Num(j-1))) exit
            Num(j) = Num(j-1)
            j = j - 1
         end do
         Num(j) = itmp
      end do
           
      return
      end subroutine Sort_subr
      



c---------------------------------------------------------------------------------------
      subroutine osc_intcp(nc_all,imix_max,imix_n,rmix_ci,StLSJ,i_ng)
      use CI_MODULE             ! module with CI-array C
      use STATE_MODULE          ! module with definitions of the state type
      include 'par.f'
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)

      integer:: nc_all, imix_max
      real, dimension(nc_all,imix_max):: rmix_ci
      integer, dimension(nc_all,imix_max):: imix_n
      type(state), dimension(nc_all):: StLSJ ! array of state structures: in intermediate coupling      
      real COF6J
      double precision   dip, deriv, dippol

      Nmax = SIZE(C,1)
      namax = SIZE(C,2)
      

      do Ni=1,min(9,Nmax)
         rji = StLSJ(Ni)%ja
         Eni = StLSJ(Ni)%en
         ex_eni =  (Eni - StLSJ(1)%en) * 27.2116
         if(Eni .gt. 0.0) cycle

         do Nf=1,nc_all
            rjf = StLSJ(Nf)%ja
            Enf = StLSJ(Nf)%en
            ex_enf = (Enf - StLSJ(1)%en) * 27.2116
            if(Enf .gt. 0.0) cycle
            if(Enf-Eni .lt. 0.0) cycle

             if( StLSJ(Nf)%pa .eq. StLSJ(Ni)%pa) cycle
            if(nint(rjf - (rji+1.0)) .gt. 0) cycle
            if(nint(rjf - abs(rji-1.0)) .lt. 0) cycle

            tmp = 0.0
            tmpdip = 0.0
            tmpderiv = 0.0

            do kst_i=1,StLSJ(Ni)%num_config
               Nnri = imix_n(Ni,kst_i)
               li = ll(Nnri)
               rli = ll(Nnri)
               rsi = ls(Nnri)
               
               do kst_f=1,StLSJ(Nf)%num_config
                  Nnrf = imix_n(Nf,kst_f)
                  lf = ll(Nnrf)
                  rlf = ll(Nnrf)
                  rsf = ls(Nnrf)

                  if(lf .gt. li+1) cycle
                  if(lf .lt. abs(li-1)) cycle
                  if(ls(Nnrf) .ne. ls(Nnri)) cycle

                  coef = rmix_ci(Nf,kst_f) * rmix_ci(Ni,kst_i)
                  if(abs(coef) .lt. 1e-5) cycle

                  coef = coef * sqrt((2.0*rjf+1.0)*(2.0*rji+1.0))*
     >                 dble(-1**(nint(rlf+rsf+rji+1))) *
     >                 COF6J(rsf,rlf,rjf,1.0,rji,rli)
                  if(coef .ne.0.0) then
                     if(i_ng .eq. 0) then
                        call osc(Nmax-1,namax,Nnrf,Nnri,lf,li,C,
     >                       fl,minf,maxf,lo,dip,deriv,dippol)
                     elseif(i_ng .eq. 1) then
                        call osc_ng(Nmax,Nnrf,Nnri,lf,li,
     >                       fl,minf,maxf,lo,dip,deriv,dippol)
c                        print*, dip,deriv,dippol
                     else 
                        print*, 'osc.f: Wrong value of i_ng=',i_ng
                        stop
                     endif
c     print*,'!!!',Nnrf,Nnri,dip/sqrt(2*li+1.0),
c     >                    dippol/sqrt(2*li+1.0)
                     tmp = tmp + coef*dippol
                     tmpdip = tmpdip + coef*dip
                     tmpderiv = tmpderiv + coef*deriv
                  endif
                  
               enddo
            enddo
            tmp = tmp / sqrt(2.0*rji+1.0)
            tmpdip = tmpdip / sqrt(2.0*rji+1.0)
            tmpderiv = tmpderiv / sqrt(2.0*rji+1.0)
            flen_pol = 2.0*abs(Enf-Eni)*tmp*tmp/3.0
            flen=2.0*abs(Enf-Eni)*tmpdip*tmpdip/3.0
            fvel=2.0*tmpderiv*tmpderiv/abs(Enf-Eni)/3.0


            write(*,'("---",2(I5,1X,A3,1X,F5.1,I5,F10.3,2X),F12.5
     >           "(",F10.5,")",F10.5)') 
     >           Nf,  StLSJ(Nf)%label, rjf, StLSJ(Nf)%pa, ex_enf,Ni, 
     >           StLSJ(Ni)%label, rji, StLSJ(Ni)%pa, ex_eni, 
     >           flen,flen_pol,fvel
            

         enddo
      enddo


      end subroutine osc_intcp
c**********************************************************************************
c     Noble gases: |l^w l2>
c     This is  M.E. of spin-orbit term. This is NOT reduced ME: <J| H | J>
c     Summation over all valence electrons is accounted here.
      double precision function H12_so_ng(State_f, State_i, J2)
      use STATE_MODULE          ! module with definitions of the state type
      use  VDCORE_MODULE
      include 'par.f'
      double precision ortint
      common /ortog/  ortint(nspmax,nspmax)

      integer J2
      type(state), intent (in) :: State_f, State_i !  state structures

      real rJ, sa_i, sa_f
      integer l1_i, l2_i, l1_f, l2_f
      double precision ang, tmp_hat, value_i, value_f
      double precision  Radint, Rel_Mass_cor ! function
      type(configuration), pointer :: config_i, config_f
      double precision  hat, hatr
      real rl
      integer l
      integer number_of_electrons
      data number_of_electrons /2/ ! set the number_of_electrons
      real*8 tmprad
c
      
      hat(l) = dsqrt(dble(2 * l + 1))
      hatr(rl) = dsqrt(2.0 * rl + 1d0)
      
      H12_so_ng = 0d0
      
      sa_i = State_i%sa
      sa_f = State_f%sa
      if(sa_f .eq. 0.0 .and. sa_i .eq. 0.0) return ! H12_so is zero in this case
      la_i = State_i%la
      la_f = State_f%la
      if(la_f .eq. 0 .and. la_i .eq. 0) return ! H12_so is zero in this case
      rJ = State_i%ja
      if(nint(2.0*rJ) .ne. J2) then
         print*, 'Wrong value of J:', rJ, real(J2)/2.0
         stop
      end if

      phase1 = (-1)**(nint(rJ + sa_i + sa_f))
      phase2 = (-1)**(nint(rJ + la_i + la_f))
      c6j_1 = COF6J(0.5,0.5,real(sa_f),1.0,real(sa_i),0.5)
      c6j_2 = COF6J(real(la_f),sa_f,rJ,sa_i,real(la_i),1.0)
      tmp_hat =  -dsqrt(1.5d0) *
     >     hat(la_i)*hat(la_f)*hatr(sa_i)*hatr(sa_f) *
     >     c6j_1 * c6j_2 

c     Note: term  State_i%pa is for  (-1)**(l1_i + l2_i) which give parity of the initial state .
c     print*, 'tmp_hat =', tmp_hat
      if(tmp_hat .eq. 0d0) return      

      config_f => State_f%first_config
      do 
c         print*,  associated(config_f), associated(config_f%NEXT),
c     >      config_f%n1,config_f%n2, config_f%value
         if( .not. associated(config_f) ) exit         
         l1_f = config_f%l1         
         l2_f = config_f%l2
         n1_f = config_f%n1
         n2_f = config_f%n2
         value_f = config_f%value
         phase_12 = (-1)**(l1_f+l2_f) 
         
         config_i => State_i%first_config
         do
            if( .not. associated(config_i) ) exit            
            n1_i = config_i%n1
            n2_i = config_i%n2
            l1_i = config_i%l1
            l2_i = config_i%l2
            value_i = config_i%value
            prodCI = value_f*value_i
            if(l1_f.eq.l1_i .and. l2_f.eq.l2_i 
     >           .and. prodCI.ne.0) then
c     deal with l1^w  
               if(n1_i.ne.1 .or. n1_f.ne.1) then
                  print*,'so.f: noble gas: wrong core orbital:',
     >                 n1_f, n1_i
                  stop
               endif
               c6j_3 =  COF6J(real(l2_i),real(l1_i),real(la_i),1.0,
     >                 real(la_f),real(l1_i))
                  ang =  -tmp_hat * ortint(n2_f,n2_i) * 
     >                 dsqrt(dble(l1_i*(l1_i + 1)*(2*l1_i + 1))) * 
     >                c6j_3 * phase_12 * phase1
                  if(ang .ne. 0d0) then
                     tmprad = Radint(n1_f,n1_i,vd_deriv(1,l1_i))
                     H12_so_ng = H12_so_ng + ang * prodCI * tmprad
                     
c     print'(6I5,2e12.4))', State_f%nst, State_i%nst!, l1_f,l2_f,l1_i,l2_i, real(ang), real(H12_so)
c                     print'("*l1^w:",8i4,1P,3e15.5)',n1_f,n1_i,n2_f,
c     >                    n2_i,l1_f,l1_i, l2_f,l2_i,
c     >                    real(ang), real(tmprad)
                  end if
c     deal with l2
                  c6j_4 = COF6J(real(l1_i),real(l2_i),real(la_i),1.0,
     >                 real(la_f),real(l2_i))
                  ang =  tmp_hat * ortint(n1_f,n1_i) * 
     >                 dsqrt(dble(l2_i*(l2_i + 1)*(2*l2_i + 1))) *
     >                 c6j_4 * phase_12 * phase2
                  if(ang .ne. 0d0) then
                     tmprad = Radint(n2_f,n2_i,vd_deriv(1,l2_i))
                     H12_so_ng = H12_so_ng + ang * prodCI * tmprad
                     
c     print'(6I5,2e12.4))', State_f%nst, State_i%nst!, l1_f,l2_f,l1_i,l2_i, real(ang), real(H12_so)
c                     print'("*l2  :",8i4,1P,3e15.5)',n1_f,n1_i,n2_f,
c     >                    n2_i,l1_f,l1_i, l2_f,l2_i,
c     >                    real(ang), real(tmprad)
                  end if
                  
                  
                  
               end if
               config_i => config_i%NEXT
            end do

         config_f => config_f%NEXT         
      end do
      H12_so_ng = H12_so_ng
c     print*, real(H12_so)
      return
      end  function H12_so_ng


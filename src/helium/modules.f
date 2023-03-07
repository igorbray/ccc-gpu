      module CI_MODULE
      
      double precision, dimension(:,:,:), allocatable :: C
      
      end module CI_MODULE
      module DM_MODULE
      
      real*8, dimension(:,:), allocatable :: DD
      integer, dimension(:), allocatable :: get_nch,get_N,get_nsp
      integer:: nchnsp_max
      
    
      end module DM_MODULE
      
      module po_module
      integer:: nspm_po
      integer, allocatable:: lo_po(:), minf_po(:), maxf_po(:)
      real, allocatable :: fl_po(:,:)
      
      end module po_module
      module VDCORE_MODULE
      
      double precision, dimension(:,:), allocatable :: vd_deriv
      
      end module VDCORE_MODULE

      module photo_module
      real, allocatable  :: dr(:,:),dv(:,:),dx(:,:)
c$$$     >   ,dr1(:,:),dv1(:,:),dx1(:,:) ! not used
      end module photo_module

      module vmat_module
      real, allocatable :: vmat01(:,:), vmat0(:,:), vmat1(:,:)
      integer, allocatable :: nchistart(:), nchistop(:)
      integer, allocatable ::  ntime(:,:), nodet(:),
     >   nntime(:),nchistartold(:,:),nchistopold(:,:)
      integer  nodeid
      logical scalapack
      end module vmat_module
      
      module gf_module
      allocatable gf(:,:,:)
      logical analytic
      integer nanalytic,nbox,kmaxgf
      real aenergyswitch, gf
      data aenergyswitch/0.0/
      end module gf_module

      module date_time_module
      character date*8,time*10,zone*5
      integer valuesin(8),valuesout(8)
      end module date_time_module

      module chil_module
      real, allocatable :: chil(:,:,:)
      integer, allocatable :: minchil(:,:)
      integer meshrr,npkstart,npkstop,ichildim,nchii
      end module chil_module
      


!-----------------------------------------------------
! ANDREY: modules required for pos - alkali scattering
!----------------------------------------------------- 
      
      module apar
      save
      integer, parameter :: maxx = 50000
      !integer :: maxx
      integer, parameter :: nmax = 6
      integer, parameter :: lmax1 = 20
      integer, parameter :: maxq = 1000
      !integer, parameter :: maxpftps = 20000
      real :: PI, znuc1, Pi4
      logical alkali   
      logical helike 
      end module apar
      
      module ftps_module
      integer :: maxp
      logical :: interpol, error  
      real,   allocatable :: pmesh(:)
      real*8, allocatable :: ftps(:,:,:) , ftpsv(:,:,:)
      integer, parameter :: nexp = 8, nr = 2**nexp 
      real, dimension (nr) :: kk
      real dkk
c      logical, parameter :: oldftps = .false.
      logical, parameter :: oldftps = .true.
      end module ftps_module
      
      module ubb_module
                                ! table of spher.harm.expan.func. Ubb(r,rho,l) 
      real, allocatable :: ubb_res(:,:,:)     
      integer ubb_min1, ubb_min2, ubb_min3 ! min values of indeces
      integer ubb_max2, ubb_max3 ! max values of indeces
      integer :: ubb_max1 != 600
      integer, parameter:: interpolate_ubb = 2 ! No (0), 1D or 2D interpolation
      real, dimension (600) :: arholg, arho ! size = ubb_max1 
      end module ubb_module
      
      module ql_module
      implicit none
      real*8, allocatable :: ql_res(:,:,:)
      integer  minql3, maxql3
      integer, dimension (0:100) :: index_ql

!     parameters for the case of direct integration with f-n five:
!     -----------------------------------------------------------
!     for p-Li it was sufficient to use maxql1 = 300 
      !integer, parameter :: maxql1 = 300, maxql2 = maxql1
      integer :: maxql1, maxql2
      
!      integer, parameter :: maxql1 = 400, maxql2 = maxql1
      
      real*8, parameter :: qlgmin = -4.0, qlgmax = 4.0 ! qlgmax = log10(qcut)
C     it was
C     maxql1 = 400; maxql2 = maxql1;
C      qlgmin = -8.0; qlgmax = 4.5;      
      !real*8, parameter :: dqlg =  (qlgmax-qlgmin)/dble(maxql1-1)
       real*8 :: dqlg
      
      
!     real*8, parameter :: dqlg = (qlgmax-qlgmin)/dble(maxql1-1)
      
!     these parameters are calculated in andrey.f:makeveq:
!      integer :: maxql1, maxql2      
!      real*8  :: qlgmin, qlgmax, dqlg ! qlgmax = log10(qcut)
      
      end module ql_module

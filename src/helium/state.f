      module state_module
      public
      
c     This is  unsymmetrized configuration.
      type configuration
      double precision value
      integer n1                ! index to first one-electron function
      integer n2                ! index to second one-electron function
      integer l1                !  orbital angular momentum of the first one-electron function
      integer l2                !  orbital angular momentum of the second one-electron function
      type(configuration), pointer :: NEXT   ! pointer to next configuration      
      end type configuration
c      
c     This is a target state
      type state
      character(LEN=3) label    ! state label
      double precision en       ! state energy in a.u.
      integer pa                ! parity
      integer la                ! orbital angular momentum
      real sa                   ! spin angular momentum
      real ja                   ! total (la + sa - vectors) angular momentum
      integer nst               ! state number
      integer num_config        ! number of unsymmetrized configurations
      type(configuration), pointer ::  first_config ! pointer to first configuration structure, it is initialized to NULL.
      end type state
      
      interface assignment (=)
          module procedure copy_st
       end interface

      contains
      
      subroutine copy_st(st_l,st_r)
      type(state),intent(out):: st_l
      type(state), intent(in):: st_r
      
      
      st_l%label = st_r%label
      st_l%en = st_r%en
      st_l%pa = st_r%pa
      st_l%la = st_r%la
      st_l%sa = st_r%sa
      st_l%ja = st_r%ja
      st_l%nst = st_r%nst
      st_l%num_config = st_r%num_config
!      st_l%first_config = st_r%first_config

      end subroutine copy_st

      end module state_module
      

      

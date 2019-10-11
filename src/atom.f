      subroutine ATOM (rmax, nr, grid)
* Helium ground state Hartee-Fock wave function

* INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
*          nr     - the number of points in r-grid
*          grid   - r-grid to calculate wavefunctions on it.

* OUTPUT:  psi1s - contains helium ground state 1s orbital
      
*      IMPLICIT REAL*8(A-H, O-Z)
*      IMPLICIT INTEGER*4(I-N)
      include 'par.f'
      include 'paratom.f'

      common /schf/ psi1s(maxr)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      dimension wd(maxr)

      character *15  os, con,  du1,  du2

      logical there
      real  grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     * gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),r1(NSIZE),
     * r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     * r7(NSIZEk7),r8(NSIZE*NSIZEis),
     * mu, c3(NSIZEis),ala(NSIZEis),
     * ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)

       integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *  ss1,l1,nu,ig,kk(NSIZEqf),
     * is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     * pp(NSIZEis),ud(NSIZEis),kz

c      Quantuum numbers

       data  in(1) /1/
       data  il(1) /0/
       data  iq(1) /2/

C  Rmax = h * ne
C  BET is the factor of the logarithmic part of the grid transformation
C  EPS precision factor
C  H is a step in the rho grid
C  NE is the number of points in the rho grid

       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 1, 700/

C  KSU1 is 1 to calculate the complete set of SCFHF wave functions
C  KSU1 should be 3 for FCHF as only the single non-core state is needed
C  MU = 1 for electrons and -1 for positrons
C  DZ and KZ are steps to increment the charge of the nucleus

       data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/

**       ne = nr
**       h = grid(nr) - grid(nr-1)

c     Some protection

       if (nr .gt. nsizek7) then
          print*, 'Nr > Nsizek7 in hfz19.f to at least', nr
          stop    'Nr > Nsizek7 in hfz19.f'
       end if 
       if (ne .gt. nsize) then
          print*,'increase NSIZE in hfz19.f to at least',ne
          stop   'increase NSIZE in hfz19.f'
       end if 

C  Z is the nuclear charge

       Z=dble(nznuc)
       if(nznuc.eq.1) iq(1)=1   ! Hydrogen -
      
C  L1 is the total orbital angular momentum of the atom

       L1 = 0
       SS1 = 1

C  DIR is unknown here

       do i = 1, is
          dir(i) = 1
       end do 

c     The Grid
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do

C  Patch for polarization potential

      r0 = 1.0
      gamma = 0.0

c       
       r   = rmax + 1.0
c
       open(66,file=con)
       open(1,file=os,form='unformatted',access='sequential')
       open(3,file=du1)
       open(4,file=du2)
       inquire(file='hfwf',exist=there)
       if (there) go to 30
       open(10, file='hfwf',form='unformatted')
       open(20, file='hfen')
c
       ne5 = ne + 5
       ne5is = ne5 * is
c
       print*, 'We are about to call SCFHF'


       call scfhf(z,gam,dir,r,bet,eps,h,teta,
     * alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     * ala2,alm6du,
     * dz,kz,
     * is,in,il,iq,
     * ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
       close(10)
       close(20)

       print*, 'We finished SCFHF'

c     Find and Store wave functions

 30    open(10, file='hfwf', status='old',form='unformatted')
       open(20, file='hfen', status='old')

       read(20, 2)  R5(1)
       print*, 'Energy=', R5(1)
       
 2     FORMAT(  ' ','Energy=',25X,E15.8,/' ','HOPMA  =' ,25X,E15.8)
          enn(n) = R5(1)
          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40       FORMAT(2X, F8.5, 2X, F20.13, 4X)

       close(10)
       close(20)

       delta=1d-5
       print*, 'Rmax=', rmax
       print*, 'Delta=', delta
       do i = 1, nr
          psi1s(i) = r2(i)
          if(abs(psi1s(i)).gt.delta) write(17,*) grid(i), psi1s(i)
       end do
c$$$       close(17)
      
      call DERIVATIVE(psi1s,wd)

       do i = 1, nr
          if(abs(psi1s(i)).gt.delta)  write(18,*) grid(i), wd(i)
       end do

       close(66)
       close(1)
       close(3)
       close(4)


       return
       end






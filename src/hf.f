c===================================================================
c     Ground state Hartee-Fock wave functions of W 5+
c===================================================================
c     
      subroutine hfz74_5 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c
      common /there1Cowan/ there1
      logical there1
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10),in(11),in(12),in(13)
     >   /1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10),il(11),il(12),il(13)
     >   /0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10),iq(11),iq(12),iq(13)
     >   /2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 1.0d-5, 13, 900/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 5d0, 10/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz74_5.f to at least', nr
         stop    'Nr > Nsizek7 in hfz74_5.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz74_5.f to at least',ne
         stop   'increase NSIZE in hfz74_5.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
c     
      r   = rmax + 1.0
c     
c
      inquire(file='hfwf.formated',exist=there1)
      there1 = .false.
      if (there1) then
         call fit_rgrid(maxr,nnmax,lnabmax,psinb,istoppsinb,
     >      enpsinb,nr,nshells,grid,EXPCUT)
         return
      end if
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
      print*,'***Enter scfhf'
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      print*,'***Out of scfhf'
      close(10)
      close(20)
c     
c     Find and Store wave functions
c           
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      wff(:,:) = 0.0
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c
      istop_max = 0
      do j = 1, is
         nn = in(j)
         ll = il(j)
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
         istop_max = max(istop_max,istoppsinb(nn,ll))
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)

c
c
      inquire(file='hfwf.formated',exist=there1)
      there1 = .true.
      if (.not. there1) then
         open(129,iostat=iostat,file='hfwf.formated',
     >        status='new',form='formatted',recl=10000)
         if(iostat .ne. 0) then
            print*,'fit_rgrid: Can not open NEW file  hfwf.formated'
            stop
         end if
         write(129,'(2i10)') is, istop_max
         write(129,'(100i5)') (in(i), i=1,is)
         write(129,'(100i5)') (il(i), i=1,is)
         write(129,'(100i5)') (iq(i), i=1,is)
         write(129,'(100E16.8)') (enn(i), i=1,is)
         
         write(129,'("This line is ignored and is here to remind",
     >        " that wave function is zero at r=0.")')
         
         do i=1,istop_max
            write(129,*) grid(i),(wff(i,ndf), ndf=1, is)
         end do
         close(129)
      endif
c    
c


      return
      end

c===================================================================
c     Ground state Hartee-Fock wave functions of Yb
c===================================================================
c     
      subroutine hfz69 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c
      common /there1Cowan/ there1
      logical there1
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10),in(11),in(12),in(13)
     >   /1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10),il(11),il(12),il(13)
     >   /0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10),iq(11),iq(12),iq(13)
     >   /2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 1.0d-5, 13, 900/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 5d0, 10/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz69.f to at least', nr
         stop    'Nr > Nsizek7 in hfz69.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz69.f to at least',ne
         stop   'increase NSIZE in hfz69.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
c     
      r   = rmax + 1.0
c     
c
      inquire(file='hfwf.formated',exist=there1)
      there1 = .false.
      if (there1) then
         call fit_rgrid(maxr,nnmax,lnabmax,psinb,istoppsinb,
     >      enpsinb,nr,nshells,grid,EXPCUT)
         return
      end if
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
      print*,'***Enter scfhf'
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      print*,'***Out of scfhf'
      close(10)
      close(20)
c     
c     Find and Store wave functions
c           
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      wff(:,:) = 0.0
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c
      istop_max = 0
      do j = 1, is
         nn = in(j)
         ll = il(j)
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
         istop_max = max(istop_max,istoppsinb(nn,ll))
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)

c
c
      inquire(file='hfwf.formated',exist=there1)
      there1 = .true.
      if (.not. there1) then
         open(129,iostat=iostat,file='hfwf.formated',
     >        status='new',form='formatted',recl=10000)
         if(iostat .ne. 0) then
            print*,'fit_rgrid: Can not open NEW file  hfwf.formated'
            stop
         end if
         write(129,'(2i10)') is, istop_max
         write(129,'(100i5)') (in(i), i=1,is)
         write(129,'(100i5)') (il(i), i=1,is)
         write(129,'(100i5)') (iq(i), i=1,is)
         write(129,'(100E16.8)') (enn(i), i=1,is)
         
         write(129,'("This line is ignored and is here to remind",
     >        " that wave function is zero at r=0.")')
         
         do i=1,istop_max
            write(129,*) grid(i),(wff(i,ndf), ndf=1, is)
         end do
         close(129)
      endif
c    
c


      return
      end

c===================================================================
c     Ground state Hartee-Fock wave functions of Hg
c===================================================================
c     
      subroutine hfz79 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c
      common /there1Cowan/ there1
      logical there1
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10),in(11),in(12),in(13),in(14)
     >   /1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10),il(11),il(12),il(13),il(14)
     >   /0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10),iq(11),iq(12),iq(13),iq(14)
     >   /2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 14, 900/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz79.f to at least', nr
         stop    'Nr > Nsizek7 in hfz79.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz79.f to at least',ne
         stop   'increase NSIZE in hfz79.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
c     
      r   = rmax + 1.0
c     
c
      inquire(file='hfwf.formated',exist=there1)
      if (there1) then
         call fit_rgrid(maxr,nnmax,lnabmax,psinb,istoppsinb,
     >      enpsinb,nr,nshells,grid,EXPCUT)
         return
      end if
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = in(j)
         ll = il(j)
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end

c===================================================================
c     Ground state Hartee-Fock wave functions of Xe+
c===================================================================
c     
      subroutine hfz54_ng (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)

      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10),in(11)
     >   /1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10),il(11)
     >   /0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10),iq(11)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 5/

c old version the PGI compiler complained about
c      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
c     >   in(10),in(11),in(12)
c     >   /1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5/
c      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
c     >   il(10),il(11),il(12)
c     >   /0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1/
c      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
c     >   iq(10),iq(11),iq(12)
c     >   /2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 5/


c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 11, 900/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz54_ng.f to at least', nr
         stop    'Nr > Nsizek7 in hfz54_ng.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz54_ng.f to at least',ne
         stop   'increase NSIZE in hfz54_ng.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = in(j)
         ll = il(j)
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      l_ng = 1
      n_ng = 5
      minf_ng = 1
      maxf_ng = istoppsinb(5,1) 
      f_ng(1:maxf_ng) = psinb(1:maxf_ng,5,1) 
c     


      close(66)
      close(1)
      close(3)
      close(4)
      return
      end
      

c===================================================================
c     Ground state Hartee-Fock wave functions of caesium
c===================================================================
c     
      subroutine hfz55 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10),in(11),in(12)
     >   /1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10),il(11),il(12)
     >   /0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10),iq(11),iq(12)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 12, 900/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz55.f to at least', nr
         stop    'Nr > Nsizek7 in hfz55.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz55.f to at least',ne
         stop   'increase NSIZE in hfz55.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = in(j)
         ll = il(j)
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end
      

c===================================================================
c     Ground state Hartee-Fock wave functions of rubidium
c===================================================================
c     
      subroutine hfz37 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9)
     >   /1, 2, 2, 3, 3, 3, 4, 4, 5/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9)
     >   /0, 0, 1, 0, 1, 2, 0, 1, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 9, 700/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz37.f to at least', nr
         stop    'Nr > Nsizek7 in hfz37.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz37.f to at least',ne
         stop   'increase NSIZE in hfz37.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = in(j)
         ll = il(j)
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end


c===================================================================
c     Ground state Hartee-Fock wave functions of cadmium
c===================================================================
c     
      subroutine hfz47 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),
     >     in(9),in(10)
     >   /1, 2, 2, 3, 3, 3, 4, 4, 4,5/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),
     >     il(9),il(10)
     >   /0, 0, 1, 0, 1, 2, 0, 1, 2,0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),
     >     iq(9),iq(10)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 10,1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 10, 700/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz37.f to at least', nr
         stop    'Nr > Nsizek7 in hfz37.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz37.f to at least',ne
         stop   'increase NSIZE in hfz37.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = in(j)
         ll = il(j)
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end


c===================================================================
c     Ground state Hartee-Fock wave functions of Kr+
c===================================================================
c     
      subroutine hfz35_ng (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)

      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8)
     >   /1, 2, 2, 3, 3, 3, 4, 4/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8)
     >   /0, 0, 1, 0, 1, 2, 0, 1/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8)
     >   /2, 2, 6, 2, 6, 10, 2, 5/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 8, 700/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz37.f to at least', nr
         stop    'Nr > Nsizek7 in hfz37.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz37.f to at least',ne
         stop   'increase NSIZE in hfz37.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = in(j)
         ll = il(j)
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
      l_ng = 1
      n_ng = 4
      minf_ng = 1
      maxf_ng = istoppsinb(4,1) 
      f_ng(1:maxf_ng) = psinb(1:maxf_ng,4,1) 
c     
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end

      

c===================================================================
c     Ground state Hartee-Fock wave functions of potassium
c===================================================================
c     
      subroutine hfz19 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1), in(2), in(3), in(4), in(5),in(6) /1, 2, 2, 3, 3, 4/
      data  il(1), il(2), il(3), il(4), il(5),il(6) /0, 0, 1, 0, 1, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6) /2, 2, 6, 2, 6, 1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 6, 700/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz19.f to at least', nr
         stop    'Nr > Nsizek7 in hfz19.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz19.f to at least',ne
         stop   'increase NSIZE in hfz19.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      istop_max = 0
      do j = 1, is
         nn = in(j)
         ll = il(j)
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         istop_max = max(istop_max,istoppsinb(nn,ll))
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)


      return
      end


      
c===================================================================
c     Ground state Hartee-Fock wave functions of argon
c===================================================================
c     
      subroutine hfz17 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1), in(2), in(3), in(4), in(5),in(6) /1, 2, 2, 3, 3, 4/
      data  il(1), il(2), il(3), il(4), il(5),il(6) /0, 0, 1, 0, 1, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6) /2, 2, 6, 1, 6, 1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 5, 700/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz19.f to at least', nr
         stop    'Nr > Nsizek7 in hfz19.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz19.f to at least',ne
         stop   'increase NSIZE in hfz19.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      istop_max = 0
      do j = 1, is
         nn = in(j)
         ll = il(j)
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         istop_max = max(istop_max,istoppsinb(nn,ll))
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)


      return
      end

c===================================================================
c     Ground state Hartee-Fock wave functions of Ar+
c===================================================================
c     
      subroutine hfz17_ng (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)

      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1), in(2), in(3), in(4), in(5) /1, 2, 2, 3, 3/
      data  il(1), il(2), il(3), il(4), il(5) /0, 0, 1, 0, 1/
      data  iq(1), iq(2), iq(3), iq(4), iq(5) /2, 2, 6, 2, 5/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 5, 700/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz19.f to at least', nr
         stop    'Nr > Nsizek7 in hfz19.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz19.f to at least',ne
         stop   'increase NSIZE in hfz19.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      istop_max = 0
      do j = 1, is
         nn = in(j)
         ll = il(j)
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         istop_max = max(istop_max,istoppsinb(nn,ll))
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
      l_ng = 1
      n_ng = 3
      minf_ng = 1
      maxf_ng = istoppsinb(3,1) 
      f_ng(1:maxf_ng) = psinb(1:maxf_ng,3,1) 
c     
      close(66)
      close(1)
      close(3)
      close(4)


      return
      end

      
c===================================================================
c     Ground state Hartee-Fock wave functions of copper
c===================================================================
c     
      subroutine hfz29 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 3d, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7)/1,2,2,3,3, 3,4/
      data  il(1),il(2),il(3),il(4),il(5),il(6),il(7)/0,0,1,0,1, 2,0/
      data  iq(1),iq(2),iq(3),iq(4),iq(5),iq(6),iq(7)/2,2,6,2,6,10,1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 7, 700/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 1d0, 3/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz19.f to at least', nr
         stop    'Nr > Nsizek7 in hfz19.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz19.f to at least',ne
         stop   'increase NSIZE in hfz19.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if(j .eq. 1) then
            nn = 1
            ll = 0
         end if 
         if(j .eq. 2) then
            nn = 2
            ll = 0
         end if 
         if(j .eq. 3) then
            nn = 2
            ll = 1
         end if 
         if(j .eq. 4) then
            nn = 3
            ll = 0
         end if 
         if(j .eq. 5) then
            nn = 3
            ll = 1
         end if 
         if(j .eq. 6) then
            nn = 3
            ll = 2
         end if 
         if(j .eq. 7) then
            nn = 4
            ll = 0
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end
      
c===================================================================
c     Ground state Hartee-Fock wave functions of Galeum
c===================================================================
c     
      subroutine hfz31 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p, 3s, 3p, 3d, 4s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8)
     >   /1,2,2,3,3, 3,4,4/
      data il(1),il(2),il(3),il(4),il(5),il(6),il(7),il(8)
     >   /0,0,1,0,1, 2,0,1/
      data iq(1),iq(2),iq(3),iq(4),iq(5),iq(6),iq(7),iq(8)
     >   /2,2,6,2,6,10,2,1/
c     
C     Rmax = h * ne
C     BET is the factor of the logarithmic part of the grid transformation
C     EPS precision factor
C     H is a step in the rho grid
C     NE is the number of points in the rho grid
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 8, 700/
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
      data  ksu1, mu, dz, kz /1, 1, 1d0, 3/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz19.f to at least', nr
         stop    'Nr > Nsizek7 in hfz19.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz19.f to at least',ne
         stop   'increase NSIZE in hfz19.f'
      end if 
C     Z is the nuclear charge
c$$$  z   = 19d0
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
C     DIR is unknown here
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if (j .eq. 1) then
            nn = 1
            ll = 0
         elseif (j .eq. 2) then
            nn = 2
            ll = 0
         elseif (j .eq. 3) then
            nn = 2
            ll = 1
         elseif (j .eq. 4) then
            nn = 3
            ll = 0
         elseif (j .eq. 5) then
            nn = 3
            ll = 1
         elseif (j .eq. 6) then
            nn = 3
            ll = 2
         elseif (j .eq. 7) then
            nn = 4
            ll = 0
         elseif (j .eq. 8) then
            nn = 4
            ll = 1
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end
      
c===================================================================
c     Ground state Hartee-Fock wave functions of lithium
c===================================================================
c     
      subroutine hfz3 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
c     
c     Quantuum numbers
c     
      data  in(1), in(2), in(3), in(4) /1, 2, 2, 3/
      data  il(1), il(2), il(3), il(4) /0, 0, 1, 0/
      data  iq(1), iq(2), iq(3), iq(4) /2, 1, 0, 0/
c$$$      data  iq(1), iq(2), iq(3), iq(4) /1, 2, 0, 0/ ! hollow Be+
c     
C     Rmax = h * ne
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 2, 700/
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
c     
c     
c     Some protection
c     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz3.f to at least', nr
         stop    'Nr > Nsizek7 in hfz3.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz3.f to at least',ne
         stop   'increase NSIZE in hfz3.f'
      end if 
C     Z is the nuclear charge
      z = dfloat(nznuc)
c$$$  z   = 3d0
C     L1 is the total orbital angular momentum of the atom
      L1 = 0
      SS1 = 2
c     
      do i = 1, is
         dir(i) = 1
      end do 
c     
c     The Grid
c     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if(j .eq. 1) then  !change to 2 for hollow Be+
            nn = 1
            ll = 0
         end if 
         if(j .eq. 2) then  !change to 1 for hollow Be+
            nn = 2
            ll = 0
         end if 
         if(j .eq. 3) then
            nn = 2
            ll = 1
         end if 
         if(j .eq. 4) then
            nn = 3
            ll = 0
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end
      
c===================================================================
c        Ground state Hartee-Fock wave functions of sodium
c===================================================================
c
      subroutine hfz11 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c                        i    i    i     i   
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid to calculate wavefunctions on it.
c
c OUTPUT:  wff - contains 1s, 2s, 2p, 3s one electron     
c                 wave functions
c      
      include 'par.f'
      include 'paratom.f'
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)

      character *15  os, con,  du1,  du2
c
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
        real*8 z,r,bet,eps,h,dz,
     * gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),r1(NSIZE),
     * r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     * r7(NSIZEk7),r8(NSIZE*NSIZEis),
     * mu, c3(NSIZEis),ala(NSIZEis),
     * ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c
       integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *  ss1,l1,nu,ig,kk(NSIZEqf),
     * is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     * pp(NSIZEis),ud(NSIZEis),kz
c
c      Quantuum numbers
c
       data  in(1), in(2), in(3), in(4) /1, 2, 2, 3/
       data  il(1), il(2), il(3), il(4) /0, 0, 1, 0/
       data  iq(1), iq(2), iq(3), iq(4) /2, 2, 6, 1/
c
C  Rmax = h * ne
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 4, 700/
       data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
       ne = nr
       h = grid(nr) - grid(nr-1)
c
c
c     Some protection
c
       if (nr .gt. nsizek7) then
          print*, 'Nr > Nsizek7 in hfz11.f to at least', nr
          stop    'Nr > Nsizek7 in hfz11.f'
       end if 
       if (ne .gt. nsize) then
          print*,'increase NSIZE in hfz11.f to at least',ne
          stop   'increase NSIZE in hfz11.f'
       end if 
C  Z is the nuclear charge
       z = dfloat(nznuc)
c$$$       z   = 11d0
C  L1 is the total orbital angular momentum of the atom
       L1 = 0
       SS1 = 2
c
       do i = 1, is
          dir(i) = 1
       end do 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
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
       call scfhf(z,gam,dir,r,bet,eps,h,teta,
     * alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     * ala2,alm6du,
     * dz,kz,
     * is,in,il,iq,
     * ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
       close(10)
       close(20)
c
c     Find and Store wave functions
c
 30    open(10, file='hfwf', status='old',form='unformatted')
       open(20, file='hfen', status='old')
       do n = 1, is
          read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
          enn(n) = R5(1)
          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40       FORMAT(2X, F8.5, 2X, F20.13, 4X)
          do i = 1, kt
             wff(i,n) = r2(i)
          end do 
       end do   
       close(10)
       close(20)
c
       do j = 1, is
          nn = 0
          ll = 0
          if(j .eq. 1) then
             nn = 1
             ll = 0
          end if 
          if(j .eq. 2) then
             nn = 2
             ll = 0
          end if 
          if(j .eq. 3) then
             nn = 2
             ll = 1
          end if 
          if(j .eq. 4) then
             nn = 3
             ll = 0
          end if 
c
          istoppsinb(nn,ll) = nr
          do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
             istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
          end do 
c
          do i = 1, istoppsinb(nn, ll)
             psinb(i, nn, ll) = wff(i, j)
          end do 
c
          enpsinb(nn, ll) = enn(j)
c          
       end do 
c
       close(66)
       close(1)
       close(3)
       close(4)
       return
       end


c===================================================================
c     Ground state Hartee-Fock wave functions of helium+
c===================================================================
c     
      subroutine hfz2 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c                        i    i    i     i   
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid to calculate wavefunctions on it.
c
c OUTPUT:  wff - contains 1s
c                 wave functions
c      
      include 'par.f'
      include 'paratom.f'
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)

      character *15  os, con,  du1,  du2
c
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
        real*8 z,r,bet,eps,h,dz,
     * gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),r1(NSIZE),
     * r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     * r7(NSIZEk7),r8(NSIZE*NSIZEis),
     * mu, c3(NSIZEis),ala(NSIZEis),
     * ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c
       integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *  ss1,l1,nu,ig,kk(NSIZEqf),
     * is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     * pp(NSIZEis),ud(NSIZEis),kz
c
c      Quantuum numbers
c
       data  in(1), in(2), in(3), in(4) /1, 2, 2, 3/
       data  il(1), il(2), il(3), il(4) /0, 0, 1, 0/
       data  iq(1), iq(2), iq(3), iq(4) /1, 0, 0, 0/
c
C  Rmax = h * ne
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 1, 700/
       data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
       ne = nr
       h = grid(nr) - grid(nr-1)
c
c
c     Some protection
c
       if (nr .gt. nsizek7) then
          print*, 'Nr > Nsizek7 in hfz1.f to at least', nr
          stop    'Nr > Nsizek7 in hfz1.f'
       end if 
       if (ne .gt. nsize) then
          print*,'increase NSIZE in hfz1.f to at least',ne
          stop   'increase NSIZE in hfz1.f'
       end if 
C  Z is the nuclear charge
       z = dfloat(nznuc)
c$$$       z   = 2d0
C  L1 is the total orbital angular momentum of the atom
       L1 = 0
       SS1 = 2
c
       do i = 1, is
          dir(i) = 1
       end do 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
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
       call scfhf(z,gam,dir,r,bet,eps,h,teta,
     * alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     * ala2,alm6du,
     * dz,kz,
     * is,in,il,iq,
     * ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
       close(10)
       close(20)
c
c     Find and Store wave functions
c
 30    open(10, file='hfwf', status='old',form='unformatted')
       open(20, file='hfen', status='old')
       do n = 1, is
          read(20, *)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
          enn(n) = R5(1)
          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40       FORMAT(2X, F8.5, 2X, F20.13, 4X)
          do i = 1, kt
             wff(i,n) = r2(i)
          end do 
       end do   
       close(10)
       close(20)
c
       do j = 1, is
          nn = 0
          ll = 0
          if(j .eq. 1) then
             nn = 1
             ll = 0
          end if 
          if(j .eq. 2) then
             nn = 2
             ll = 0
          end if 
          if(j .eq. 3) then
             nn = 2
             ll = 1
          end if 
          if(j .eq. 4) then
             nn = 3
             ll = 0
          end if 
c
          istoppsinb(nn,ll) = nr
          do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
             istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
          end do 
c
          do i = 1, istoppsinb(nn, ll)
             psinb(i, nn, ll) = wff(i, j)
          end do 
c
          enpsinb(nn, ll) = enn(j)
c          
       end do 
c
       close(66)
       close(1)
       close(3)
       close(4)
       return
       end

c===================================================================
c     Ground state Hartee-Fock wave functions of neutral Ne
c===================================================================
c     
      subroutine hfz9 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)

      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
     
c     Quantuum numbers
     
      data  in(1), in(2), in(3), in(4) /1, 2, 2, 3/
      data  il(1), il(2), il(3), il(4) /0, 0, 1, 0/
      data  iq(1), iq(2), iq(3), iq(4) /2, 1, 6, 1/
c     
C     Rmax = h * ne
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 3, 700/
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
     
c     Some protection
     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz9.f to at least', nr
         stop    'Nr > Nsizek7 in hfz9.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz9.f to at least',ne
         stop   'increase NSIZE in hfz9.f'
      end if 
C     Z is the nuclear charge
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 1
      SS1 = 2
     
      do i = 1, is
         dir(i) = 1
      end do 
     
c     The Grid
     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
!         read(20, 2)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if(j .eq. 1) then
            nn = 1
            ll = 0
         end if 
         if(j .eq. 2) then
            nn = 2
            ll = 0
         end if 
         if(j .eq. 3) then
            nn = 2
            ll = 1
         end if 
         if(j .eq. 4) then
            nn = 3
            ll = 0
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      l_ng = 1
      n_ng = 2
      minf_ng = 1
      maxf_ng = istoppsinb(2,1) 
      f_ng(1:maxf_ng) = psinb(1:maxf_ng,2,1) 
c
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end
      
      
c===================================================================
c     Ground state Hartee-Fock wave functions of neutral Ne
c===================================================================
c     
      subroutine hfz9_ng (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)

      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
     
c     Quantuum numbers
     
      data  in(1), in(2), in(3), in(4) /1, 2, 2, 3/
      data  il(1), il(2), il(3), il(4) /0, 0, 1, 0/
      data  iq(1), iq(2), iq(3), iq(4) /2, 2, 5, 0/
c     
C     Rmax = h * ne
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 3, 700/
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
     

      print*,'Making Ne+ core orbitals: hfz9_ng'

c     Some protection
     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz9_ng.f to at least', nr
         stop    'Nr > Nsizek7 in hfz9_ng.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz9_ng.f to at least',ne
         stop   'increase NSIZE in hfz9_ng.f'
      end if 
C     Z is the nuclear charge
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 1
      SS1 = 2
     
      do i = 1, is
         dir(i) = 1
      end do 
     
c     The Grid
     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
!         read(20, 2)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if(j .eq. 1) then
            nn = 1
            ll = 0
         end if 
         if(j .eq. 2) then
            nn = 2
            ll = 0
         end if 
         if(j .eq. 3) then
            nn = 2
            ll = 1
         end if 
         if(j .eq. 4) then
            nn = 3
            ll = 0
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      l_ng = 1
      n_ng = 2
      minf_ng = 1
      maxf_ng = istoppsinb(2,1) 
      f_ng(1:maxf_ng) = psinb(1:maxf_ng,2,1) 
c
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end
      
      
c===================================================================
c     Ground state Hartee-Fock wave functions of doubly ionized Ne++
c===================================================================
c     
      subroutine hfz10_2 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
     
c     Quantuum numbers
     
      data  in(1), in(2), in(3), in(4) /1, 2, 2, 3/
      data  il(1), il(2), il(3), il(4) /0, 0, 1, 0/
      data  iq(1), iq(2), iq(3), iq(4) /2, 2, 4, 0/
c     
C     Rmax = h * ne
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 3, 700/
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
     
c     Some protection
     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz10_2.f to at least', nr
         stop    'Nr > Nsizek7 in hfz10_2.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz10_2.f to at least',ne
         stop   'increase NSIZE in hfz10_2.f'
      end if 
C     Z is the nuclear charge
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 1
      SS1 = 2
     
      do i = 1, is
         dir(i) = 1
      end do 
     
c     The Grid
     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
!         read(20, 2)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if(j .eq. 1) then
            nn = 1
            ll = 0
         end if 
         if(j .eq. 2) then
            nn = 2
            ll = 0
         end if 
         if(j .eq. 3) then
            nn = 2
            ll = 1
         end if 
         if(j .eq. 4) then
            nn = 3
            ll = 0
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end

      
      
c===================================================================
c     Ground state Hartee-Fock wave functions of neutral B
c===================================================================
      subroutine hfz5 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)

      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
     
c     Quantuum numbers
     
      data  in(1), in(2), in(3) /1, 2, 2/
      data  il(1), il(2), il(3) /0, 0, 1/
      data  iq(1), iq(2), iq(3) /2, 2, 1/

c old version that PGI compiler didn't like
c      data  in(1), in(2), in(3), in(4) /1, 2, 2/
c      data  il(1), il(2), il(3), il(4) /0, 0, 1/
c      data  iq(1), iq(2), iq(3), iq(4) /2, 2, 1/

c     
C     Rmax = h * ne
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 3, 700/
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
     
c     Some protection
     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz5.f to at least', nr
         stop    'Nr > Nsizek7 in hfz5.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz5.f to at least',ne
         stop   'increase NSIZE in hfz5.f'
      end if 
C     Z is the nuclear charge
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 1
      SS1 = 2
     
      do i = 1, is
         dir(i) = 1
      end do 
     
c     The Grid
     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
!         read(20, 2)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if(j .eq. 1) then
            nn = 1
            ll = 0
         end if 
         if(j .eq. 2) then
            nn = 2
            ll = 0
         end if 
         if(j .eq. 3) then
            nn = 2
            ll = 1
         end if 
         if(j .eq. 4) then
            nn = 3
            ll = 0
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      l_ng = 1
      n_ng = 2
      minf_ng = 1
      maxf_ng = istoppsinb(2,1) 
      f_ng(1:maxf_ng) = psinb(1:maxf_ng,2,1) 
c
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end
      
c===================================================================
c     Ground state Hartee-Fock wave functions of neutral C
c===================================================================
      subroutine hfz6 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)

      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
     
c     Quantuum numbers
     
      data  in(1), in(2), in(3) /1, 2, 2/
      data  il(1), il(2), il(3) /0, 0, 1/
      data  iq(1), iq(2), iq(3) /2, 2, 2/
c     
C     Rmax = h * ne
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 3, 700/
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
     
c     Some protection
     
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz6.f to at least', nr
         stop    'Nr > Nsizek7 in hfz6.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz6.f to at least',ne
         stop   'increase NSIZE in hfz6.f'
      end if 
C     Z is the nuclear charge
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 1
      SS1 = 2
     
      do i = 1, is
         dir(i) = 1
      end do 
     
c     The Grid
     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
!         read(20, 2)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if(j .eq. 1) then
            nn = 1
            ll = 0
         end if 
         if(j .eq. 2) then
            nn = 2
            ll = 0
         end if 
         if(j .eq. 3) then
            nn = 2
            ll = 1
         end if 
         if(j .eq. 4) then
            nn = 3
            ll = 0
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      l_ng = 1
      n_ng = 2
      minf_ng = 1
      maxf_ng = istoppsinb(2,1) 
      f_ng(1:maxf_ng) = psinb(1:maxf_ng,2,1) 
c
      close(66)
      close(1)
      close(3)
      close(4)
      return
      end
c===================================================================
c     Ground state Hartee-Fock wave functions of neutral N
c===================================================================
      subroutine hfz7 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)

      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
     
c     Quantuum numbers
     
      data  in(1), in(2), in(3) /1, 2, 2/
      data  il(1), il(2), il(3) /0, 0, 1/
      data  iq(1), iq(2), iq(3) /2, 2, 3/
c     
C     Rmax = h * ne
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 3, 700/
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
     
c     Some protection
      print*, 'start N: test'
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz7.f to at least', nr
         stop    'Nr > Nsizek7 in hfz7.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz7.f to at least',ne
         stop   'increase NSIZE in hfz7.f'
      end if 
C     Z is the nuclear charge
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 1
      SS1 = 2
     
      do i = 1, is
         dir(i) = 1
      end do 
     
c     The Grid
     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
!         read(20, 2)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if(j .eq. 1) then
            nn = 1
            ll = 0
         end if 
         if(j .eq. 2) then
            nn = 2
            ll = 0
         end if 
         if(j .eq. 3) then
            nn = 2
            ll = 1
         end if 
         if(j .eq. 4) then
            nn = 3
            ll = 0
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      l_ng = 1
      n_ng = 2
      minf_ng = 1
      maxf_ng = istoppsinb(2,1) 
      f_ng(1:maxf_ng) = psinb(1:maxf_ng,2,1) 
c
      close(66)
      close(1)
      close(3)
      close(4)

      print*, 'finish N: test'

      return
      end
      
c===================================================================
c     Ground state Hartee-Fock wave functions of neutral O
c===================================================================
      subroutine hfz8 (rmax, nr, grid, expcut, nznuc, gamma, r0)
c     i    i    i     i   
c     
c     INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c     nr     - the number of points in r-grid
c     grid   - r-grid to calculate wavefunctions on it.
c     
c     OUTPUT:  wff - contains 1s, 2s, 2p one-electron     
c     wave functions
c     
      include 'par.f'
      include 'paratom.f'
c     
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)

      common /noblegas/ l_ng, n_ng, minf_ng, maxf_ng, f_ng(maxr)
      
      character *15  os, con,  du1,  du2
c     
      logical there
      real wff(nsizek7, nshells) , grid(nr), enn(nshells)
      real*8 z,r,bet,eps,h,dz,
     *   gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),
     *   r1(NSIZE),r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     *   r7(NSIZEk7),r8(NSIZE*NSIZEis),
     *   mu, c3(NSIZEis),ala(NSIZEis),
     *   ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)
c     
      integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *   ss1,l1,nu,ig,kk(NSIZEqf),
     *   is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     *   pp(NSIZEis),ud(NSIZEis),kz
     
c     Quantuum numbers
     
      data  in(1), in(2), in(3) /1, 2, 2/
      data  il(1), il(2), il(3) /0, 0, 1/
      data  iq(1), iq(2), iq(3) /2, 2, 4/
c     
C     Rmax = h * ne
      data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 3, 700/
      data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
      data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/
      ne = nr
      h = grid(nr) - grid(nr-1)
     
c     Some protection
      print*, 'start N: test'
      if (nr .gt. nsizek7) then
         print*, 'Nr > Nsizek7 in hfz8.f to at least', nr
         stop    'Nr > Nsizek7 in hfz8.f'
      end if 
      if (ne .gt. nsize) then
         print*,'increase NSIZE in hfz8.f to at least',ne
         stop   'increase NSIZE in hfz8.f'
      end if 
C     Z is the nuclear charge
      z = dfloat(nznuc)
C     L1 is the total orbital angular momentum of the atom
      L1 = 1
      SS1 = 2
     
      do i = 1, is
         dir(i) = 1
      end do 
     
c     The Grid
     
      kt  = nr
      do i = 1, kt
         r7(i) = grid(i)
      end do
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
      call scfhf(z,gam,dir,r,bet,eps,h,teta,
     *   alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     *   ala2,alm6du,
     *   dz,kz,
     *   is,in,il,iq,
     *   ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     *   ne5,ne5is,kt,ud,ksu1,gamma,r0)
      close(10)
      close(20)
c     
c     Find and Store wave functions
c     
 30   open(10, file='hfwf', status='old',form='unformatted')
      open(20, file='hfen', status='old')
      do n = 1, is
         read(20, *)  R5(1)
!         read(20, 2)  R5(1)
 2       FORMAT(  ' ','Energy=',25X,E15.8,
     *      /' ','HOPMA  =' ,25X,E15.8)
         enn(n) = R5(1)
         read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$  read(10, 40) (R7(IW), R2(IW), IW = 1, KT )
 40      FORMAT(2X, F8.5, 2X, F20.13, 4X)
         do i = 1, kt
            wff(i,n) = r2(i)
         end do 
      end do   
      close(10)
      close(20)
c     
      do j = 1, is
         nn = 0
         ll = 0
         if(j .eq. 1) then
            nn = 1
            ll = 0
         end if 
         if(j .eq. 2) then
            nn = 2
            ll = 0
         end if 
         if(j .eq. 3) then
            nn = 2
            ll = 1
         end if 
         if(j .eq. 4) then
            nn = 3
            ll = 0
         end if 
c     
         istoppsinb(nn,ll) = nr
         do while (abs(wff(istoppsinb(nn,ll), j)) .lt. expcut)
            istoppsinb(nn,ll) = istoppsinb(nn,ll) - 1
         end do 
c     
         do i = 1, istoppsinb(nn, ll)
            psinb(i, nn, ll) = wff(i, j)
         end do 
c     
         enpsinb(nn, ll) = enn(j)
c     
      end do 
c     
      l_ng = 1
      n_ng = 2
      minf_ng = 1
      maxf_ng = istoppsinb(2,1) 
      f_ng(1:maxf_ng) = psinb(1:maxf_ng,2,1) 
c
      close(66)
      close(1)
      close(3)
      close(4)

      print*, 'finish N: test'

      return
      end
      

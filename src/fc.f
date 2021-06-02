c===================================================================
c       Frozen core Hartee-Fock wave functions of Yb
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz74_5(rmax, nr, grid, expcut, regcut, q, nc, wff, 
     >   iwff, nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10),in(11),in(12),in(13),in(14)
     >   /1, 2, 2, 3, 3,  3, 4, 4,  4, 4, 5, 5, 6, 7/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10),il(11),il(12),il(13),il(14)
     >   /0, 0, 1, 0, 1,  2, 0, 1,  2, 3, 0, 1, 0, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10),iq(11),iq(12),iq(13),iq(14)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 10,14, 2, 6, 1, 0/
c
C  See comments in the hfz11.f file
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 13, 900/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       iwff = 1
       wff(1) = 1.0
c$$$       print*,'fcz74_5 fails for unknown reason'
c$$$       return
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(6,lv+1)
       if (lv.eq.2) in(is) = 5
       if (lv.eq.3) in(is) = 5
       if (lv.eq.4) in(is) = 5
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz74_5'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='new')
c$$$       else
c$$$          write(ench,'(1p,"_",e10.4)') q*q
c$$$          open(10, file='fcwf'//ench, form='unformatted')
c$$$          open(20, file='fcen'//ench)
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       if (nc.eq.0) then
c$$$          close(10)
c$$$          close(20)
c$$$       else 
c$$$          close(10,status='delete')
c$$$          close(20,status='delete')
c$$$       endif 
c$$$c
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$          
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb
c$$$c
c$$$          nn = in(is) + n - 1
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end

c===================================================================
c       Frozen core Hartee-Fock wave functions of Yb
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz69(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10),in(11),in(12),in(13),in(14)
     >   /1, 2, 2, 3, 3,  3, 4, 4,  4, 4, 5, 5, 6, 7/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10),il(11),il(12),il(13),il(14)
     >   /0, 0, 1, 0, 1,  2, 0, 1,  2, 3, 0, 1, 0, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10),iq(11),iq(12),iq(13),iq(14)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 10,14, 2, 6, 0, 0/
c
C  See comments in the hfz11.f file
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 14, 900/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       iwff = 1
       wff(1) = 1.0
       print*,'fcz69 fails for unknown reason'
       return
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(6,lv+1)
       if (lv.eq.2) in(is) = 5
       if (lv.eq.3) in(is) = 5
       if (lv.eq.4) in(is) = 5
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz69'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='new')
c$$$       else
c$$$          write(ench,'(1p,"_",e10.4)') q*q
c$$$          open(10, file='fcwf'//ench, form='unformatted')
c$$$          open(20, file='fcen'//ench)
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       if (nc.eq.0) then
c$$$          close(10)
c$$$          close(20)
c$$$       else 
c$$$          close(10,status='delete')
c$$$          close(20,status='delete')
c$$$       endif 
c$$$c
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$          
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb
c$$$c
c$$$          nn = in(is) + n - 1
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end

c===================================================================
c       Frozen core Hartee-Fock wave functions of Hg
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz79(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10),in(11),in(12),in(13),in(14),in(15)
     >   /1, 2, 2, 3, 3,  3, 4, 4,  4, 4, 5, 5, 5, 6, 7/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10),il(11),il(12),il(13),il(14),il(15)
     >   /0, 0, 1, 0, 1,  2, 0, 1,  2, 3, 0, 1, 2, 0, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10),iq(11),iq(12),iq(13),iq(14),iq(15)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 10,14, 2, 6,10, 0, 1/
c
C  See comments in the hfz11.f file
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 15, 900/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(6,lv+1)
       if (lv.eq.2) in(is) = 6
       if (lv.eq.3) in(is) = 5
       if (lv.eq.4) in(is) = 5
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz79'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='new')
c$$$       else
c$$$          write(ench,'(1p,"_",e10.4)') q*q
c$$$          open(10, file='fcwf'//ench, form='unformatted')
c$$$          open(20, file='fcen'//ench)
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       if (nc.eq.0) then
c$$$          close(10)
c$$$          close(20)
c$$$       else 
c$$$          close(10,status='delete')
c$$$          close(20,status='delete')
c$$$       endif 
c$$$c
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$          
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb
c$$$c
c$$$          nn = in(is) + n - 1
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end

c===================================================================
c       Frozen core Hartee-Fock wave functions of caesium
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz55(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10),in(11),in(12),in(13)
     >   /1, 2, 2, 3, 3,  3, 4, 4,  4, 5, 5, 6, 7/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10),il(11),il(12),il(13)
     >   /0, 0, 1, 0, 1,  2, 0, 1,  2, 0, 1, 0, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10),iq(11),iq(12),iq(13)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 0, 1/
c
C  See comments in the hfz11.f file
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 13, 900/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(6,lv+1)
       if (lv.eq.2) in(is) = 5
       if (lv.eq.3) in(is) = 4
       if (lv.eq.4) in(is) = 5
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz55'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
       return
       end

c===================================================================
c       Frozen core Hartee-Fock wave functions of rubidium
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz37(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10)
     >   /1, 2, 2, 3, 3,  3, 4, 4,  5, 6/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10)
     >   /0, 0, 1, 0, 1,  2, 0, 1,  0, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 0, 1/
c
C  See comments in the hfz11.f file
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 10, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(5,lv+1)
       if (lv.eq.2) in(is) = 4
       if (lv.eq.3) in(is) = 4
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 1
C  Don't orthogonalise to the core states
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz11'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
       return
       end


c===================================================================
c       Frozen core Hartee-Fock wave functions of cadmium
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz47(rmax, nr, grid, expcut, regcut, q, nc,wff,iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),
     >   in(10), in(11)
     >   /1, 2, 2, 3, 3,  3, 4, 4,  4,5, 6/
      data  il(1), il(2), il(3), il(4), il(5),il(6),il(7),il(8),il(9),
     >   il(10), il(11)
     >   /0, 0, 1, 0, 1,  2, 0, 1,  2, 0, 0/
      data  iq(1), iq(2), iq(3), iq(4), iq(5),iq(6),iq(7),iq(8),iq(9),
     >   iq(10), iq(11)
     >   /2, 2, 6, 2, 6, 10, 2, 6, 10,0, 1/
c
C  See comments in the hfz11.f file
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 11, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(5,lv+1)
       if (lv.eq.2) in(is) = 4
       if (lv.eq.3) in(is) = 4
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 1
C  Don't orthogonalise to the core states
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz11'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
       return
       end


c===================================================================
c       Frozen core Hartee-Fock wave functions of potassium
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz19(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
       data  in(1),in(2),in(3),in(4),in(5),in(6),in(7)/1,2,2,3,3,4,5/
       data  il(1),il(2),il(3),il(4),il(5),il(6),il(7)/0,0,1,0,1,0,0/
       data  iq(1),iq(2),iq(3),iq(4),iq(5),iq(6),iq(7)/2,2,6,2,6,0,1/
c
C  See comments in the hfz11.f file
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 7, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
c$$$       z    = 19d0
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(4,lv+1)
       if (lv.eq.2) in(is) = 3
c$$$       if (nc.eq.0.and.lv .gt. 3) stop 'Have not coded for Latom > 3'
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
C  Don't orthogonalise to the core states
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states. May no longer be necessary.
       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz11'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='new')
c$$$       else
c$$$          write(ench,'(1p,"_",e10.4)') q*q
c$$$          open(10, file='fcwf'//ench, form='unformatted')
c$$$          open(20, file='fcen'//ench)
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       if (nc.eq.0) then
c$$$          close(10)
c$$$          close(20)
c$$$       else 
c$$$          close(10,status='delete')
c$$$          close(20,status='delete')
c$$$       endif 
c$$$c
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$          
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb
c$$$c
c$$$          nn = n + 3
c$$$          if (lv .eq. 2) nn = n + 2
c$$$          if (lv .ge. 3) nn = n + lv
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end

c===================================================================
c       Frozen core Hartee-Fock wave functions of argon
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz17(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
       data  in(1),in(2),in(3),in(4),in(5),in(6),in(7)/1,2,2,3,3,4,5/
       data  il(1),il(2),il(3),il(4),il(5),il(6),il(7)/0,0,1,0,1,0,0/
       data  iq(1),iq(2),iq(3),iq(4),iq(5),iq(6),iq(7)/2,2,6,0,6,1,0/
c
C  See comments in the hfz11.f file
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 6, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
c$$$       z    = 19d0
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(4,lv+1)
       if (lv.eq.2) in(is) = 3
       if (nc.eq.0.and.lv .eq. 2) then
          print*, 'bypassing Latom = 2'
          return
       endif 
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
C  Don't orthogonalise to the core states
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states. May no longer be necessary.
       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz11'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='new')
c$$$       else
c$$$          write(ench,'(1p,"_",e10.4)') q*q
c$$$          open(10, file='fcwf'//ench, form='unformatted')
c$$$          open(20, file='fcen'//ench)
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       if (nc.eq.0) then
c$$$          close(10)
c$$$          close(20)
c$$$       else 
c$$$          close(10,status='delete')
c$$$          close(20,status='delete')
c$$$       endif 
c$$$c
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$          
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb
c$$$c
c$$$          nn = n + 3
c$$$          if (lv .eq. 2) nn = n + 2
c$$$          if (lv .ge. 3) nn = n + lv
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end

c===================================================================
c       Frozen core Hartee-Fock wave functions of copper
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz29(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch, ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8)
     >   /1,2,2,3,3, 3,4,5/
      data  il(1),il(2),il(3),il(4),il(5),il(6),il(7),il(8)
     >   /0,0,1,0,1, 2,0,0/
      data  iq(1),iq(2),iq(3),iq(4),iq(5),iq(6),iq(7),iq(8)
     >   /2,2,6,2,6,10,0,1/
c
C  See comments in the hfz11.f file
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 8, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(4,lv+1)
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
C  Do orthogonalise to the core states
          ud(i)  = 1
C  Don't orthogonalise to the core states
c$$$          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
c$$$       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz29'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='new')
c$$$       else
c$$$          open(10, file='fcwf', form='unformatted')
c$$$          open(20, file='fcen')
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       close(10)
c$$$       close(20)
c$$$c
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$          
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb
c$$$c
c$$$          nn = n + 3
c$$$          if (lv .ge. 3) nn = n + lv
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end
c===================================================================
c       Frozen core Hartee-Fock wave functions of Ga
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz31(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch, ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
      data  in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8)
     >   /1,2,2,3,3, 3,4,5/
      data  il(1),il(2),il(3),il(4),il(5),il(6),il(7),il(8)
     >   /0,0,1,0,1, 2,0,1/
      data  iq(1),iq(2),iq(3),iq(4),iq(5),iq(6),iq(7),iq(8)
     >   /2,2,6,2,6,10,2,1/
c
C  See comments in the hfz11.f file
C     KSU1 is 1 to calculate the complete set of SCFHF wave functions
C     KSU1 should be 3 for FCHF as only the single non-core state is needed
C     MU = 1 for electrons and -1 for positrons
C     DZ and KZ are steps to increment the charge of the nucleus
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 8, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = max(4,lv+1)
       if (lv.eq.0) in(is)=5
       L1     = lv
       SS1    = 2
c
       print*,'FCZ31;in(is),lv,nb:',in(is),lv,nb
       do i = 1, is
          dir(i) = 1
C  Do orthogonalise to the core states
          ud(i)  = 1
C  Don't orthogonalise to the core states
c$$$          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(is-1) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz31'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
       return
       end
c===================================================================
c       Frozen core Hartee-Fock wave functions of lithium
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz3(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
       data  in(1), in(2), in(3), in(4) ,in(5) /1, 2, 3, 3, 4/
       data  il(1), il(2), il(3), il(4) ,il(5) /0, 0, 0, 0, 0/
       data  iq(1), iq(2), iq(3), iq(4) ,iq(5) /2, 0, 1, 0, 0/
c$$$       data  iq(1), iq(2), iq(3), iq(4) ,iq(5) /0, 2, 1, 0, 0/ !hollow Be+
c
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 3, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
c$$$       z    = 3d0
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = 2
c$$$       in(is) = 3 !hollow Be+
       if (lv .ge. 1) in(is) = lv + 1
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 1
C  Don't orthogonalise to the core states
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(2) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz11'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
       
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='new')
c$$$       else
c$$$          write(ench,'(1p,"_",e10.4)') q*q
c$$$          open(10, file='fcwf'//ench, form='unformatted')
c$$$          open(20, file='fcen'//ench)
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       if (nc.eq.0) then
c$$$          close(10)
c$$$          close(20)
c$$$       else 
c$$$          close(10,status='delete')
c$$$          close(20,status='delete')
c$$$       endif 
c$$$c
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$          
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb
c$$$c
c$$$          nn = n + 1
c$$$          if (lv .gt. 1) nn = lv + n
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end

c===================================================================
c       Frozen core Hartee-Fock wave functions of sodium
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz11(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
       data  in(1), in(2), in(3), in(4) ,in(5) /1, 2, 2, 3, 4/
       data  il(1), il(2), il(3), il(4) ,il(5) /0, 0, 1, 0, 0/
       data  iq(1), iq(2), iq(3), iq(4) ,iq(5) /2, 2, 6, 0, 1/
c$$$       data  in(1), in(2), in(3), in(4) /1, 2, 2, 3/
c$$$       data  il(1), il(2), il(3), il(4) /0, 0, 1, 0/
c$$$       data  iq(1), iq(2), iq(3), iq(4) /2, 2, 6, 1/
c
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 5, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
c$$$       z    = 11d0
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = 3
       if (lv .gt. 2) in(is) = lv + 1
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 1
C  Don't orthogonalise to the core states
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(4) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz11'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='new')
c$$$       else
c$$$          write(ench,'(1p,"_",e10.4)') q*q
c$$$          open(10, file='fcwf'//ench, form='unformatted')
c$$$          open(20, file='fcen'//ench)
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       if (nc.eq.0) then
c$$$          close(10)
c$$$          close(20)
c$$$       else 
c$$$          close(10,status='delete')
c$$$          close(20,status='delete')
c$$$       endif 
c$$$c     
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$          
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb
c$$$c
c$$$          nn = n + 2
c$$$          if (lv .gt. 2) nn = lv + n
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end

c===================================================================
c       Frozen core Hartee-Fock wave functions of helium
c===================================================================
c
c                        i    i    i      i      i     i  i    o    o 
      subroutine fcz2(rmax, nr, grid, expcut, regcut, q, nc, wff, iwff,
     >   nb, lv, nznuc, SS1, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch,ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
       data  in(1), in(2), in(3), in(4) ,in(5) /1, 2, 3, 3, 4/
       data  il(1), il(2), il(3), il(4) ,il(5) /0, 0, 0, 0, 0/
       data  iq(1), iq(2), iq(3), iq(4) ,iq(5) /1, 1, 0, 0, 0/
c
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 2, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c
       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = 2
       if (lv .gt. 1) in(is) = lv + 1
       L1     = lv
c$$$       SS1    = 1
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(2) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz11'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(SS1)//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(SS1)//ch(lv), status='new',
c$$$     >       form='unformatted')
c$$$          open(20, file='fcen'//ch(SS1)//ch(lv), status='new')
c$$$       else
c$$$          open(10, file='fcwf', form='unformatted')
c$$$          open(20, file='fcen')
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       close(10)
c$$$       close(20)
c$$$c
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(SS1)//ch(lv),status='old',
c$$$     >       form='unformatted')
c$$$          open(20, file='fcen'//ch(SS1)//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$          
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb * 0
c$$$c
c$$$          nn = n + 1
c$$$          if (lv .gt. 1) nn = lv + n
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end
c===================================================================
c     Excited state Hartee-Fock wave functions of singly ionized Ne+
c===================================================================
c
      subroutine fcz9 (rmax,nr,grid,expcut,regcut,q,nc,wff,iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch, ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
       data  in(1), in(2), in(3), in(4) ,in(5) /1, 2, 3, 3, 4/
       data  il(1), il(2), il(3), il(4) ,il(5) /0, 0, 1, 0, 0/
       data  iq(1), iq(2), iq(3), iq(4) ,iq(5) /2, 2, 4, 1, 0/
       
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 4, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c

       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = 3
       if (lv .gt. 1) in(is) = lv + 1
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 1
C  Don't orthogonalise to the core states
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(2) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz11'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
c$$$       open(2,file=os,form='unformatted',access='sequential')
c$$$       open(1,file=vo,
c$$$     >    form='unformatted',access='sequential')
c$$$       open(3,file=du1)
c$$$       open(4,file=du2)
c$$$       open(66,file=con)
c$$$       if (nc.eq.0) then
c$$$          inquire(file='fcwf'//ch(lv),exist=there)
c$$$          if (there) go to 10
c$$$          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='new')
c$$$       else
c$$$          open(10, file='fcwf', form='unformatted')
c$$$          open(20, file='fcen')
c$$$       end if 
c$$$          
c$$$c
c$$$       
c$$$       call fchf(z,gam,dir,r,bet,eps,h,teta,
c$$$     * ala,ala1,ala2,ala3,c3,
c$$$     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
c$$$     * is,in,il,iq,
c$$$     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
c$$$     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c$$$c
c$$$
c$$$       close(10)
c$$$       close(20)
c$$$c
c$$$c     Find and Store wave functions
c$$$c
c$$$ 10    continue
c$$$       if (nc.eq.0) then
c$$$          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
c$$$          open(20, file='fcen'//ch(lv), status='old')
c$$$       else
c$$$          open(10, file='fcwf', status='old', form='unformatted')
c$$$          open(20, file='fcen', status='old')
c$$$       end if 
c$$$
c$$$       if (nb.gt.0.and.nc.gt.0) then
c$$$          print *,'Either NB or NC must be zero'
c$$$          stop 'Either NB or NC must be zero'
c$$$       end if
c$$$
c$$$       if (nc.eq.1) then
c$$$          phase = exp((0.0,1.0) * r5(ne+4))
c$$$          read(10) (R7(IW), r2(IW), IW = 1, KT )
c$$$          iwff = 1
c$$$          do while (abs(r2(iwff)).lt.regcut)
c$$$             r2(iwff) = 0.0
c$$$             iwff = iwff + 1
c$$$          end do
c$$$C  Normalize the continuum wave function to range between +/- 1.0
c$$$c$$$          c = sqrt(acos(-1.0) * q)
c$$$C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
c$$$          c = sqrt(2.0 * q)
c$$$          do i = iwff, kt
c$$$             wff(i) = r2(i) * c
c$$$          end do 
c$$$       end if 
c$$$
c$$$       do n = 1, nb
c$$$c          nn = n + 1
c$$$          nn = n + 2
c$$$          if (lv .gt. 1) nn = lv + n
c$$$          if (nn.gt.nnmax) then
c$$$             print*,'NN > NNMAX',nn,nnmax
c$$$             stop 'increase NNMAX'
c$$$          end if
c$$$          
c$$$          read(20, *)  R5(1)
c$$$!          read(20, 22)  R5(1)
c$$$ 22       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$          enpsinb(nn, lv) = R5(1)
c$$$c
c$$$          read(10) (R7(IW), R2(IW), IW = 1, KT )
c$$$          do i = 1, kt
c$$$             psinb(i,nn,lv) = r2(i)
c$$$          end do 
c$$$          istoppsinb(nn,lv) = kt
c$$$          do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
c$$$             istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
c$$$          end do 
c$$$       end do   
c$$$ 2       FORMAT(  ' ','Energy=',25X,E15.8,
c$$$     *      /' ','HOPMA  =' ,25X,E15.8)
c$$$       close(10)
c$$$       close(20)
c$$$
c$$$       close(1)
c$$$       close(2)
c$$$       close(3)
c$$$       close(4)
c$$$       close(66)
       return
       end

c===================================================================
c     Excited state Hartee-Fock wave functions of neutral B, singly ionized C+
c===================================================================
c
      subroutine fcz5 (rmax,nr,grid,expcut,regcut,q,nc,wff,iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch, ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
       data  in(1), in(2), in(3) /1, 2, 2/
       data  il(1), il(2), il(3) /0, 0, 1/
       data  iq(1), iq(2), iq(3) /2, 2, 1/
       
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 3, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c

       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = 3
       if (lv .gt. 1) in(is) = lv + 1
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 1
C  Don't orthogonalise to the core states
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(2) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz5'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
       return
       end

c===================================================================
c     Excited state Hartee-Fock wave functions of neutral Al,..., Ar 5+
c===================================================================
c NOTE: need to change 
c     1. common bloks in(), il(), iq()
c     2. value of 'is'  in 'data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 5, 700/ '
c         to be equal to the maximum value in 1. common bloks (is =5 in this case).
      subroutine fcz13 (rmax,nr,grid,expcut,regcut,q,nc,wff,iwff,
     >   nb, lv, nznuc, phase, gamma, r0)
c        i   i
c
c INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
c          nr     - the number of points in r-grid
c          grid   - r-grid on which wavefunctions will be calculated
c          q      - the momentum of the continuum wavefunction
c          nc     - number of continuum states (0 or 1)
c          nb     - number of bound states
c          lv     - Orbital momentum of the valence electron
c
c OUTPUT:  wff - contains one electron continuum wavefunction if NC = 1.
c      
      include 'par.f'
      include 'paratom.f'
C  NS1 is the max number of core shells
      PARAMETER (ns1=15)
C  NS2 is the max number of shell coefficients
      PARAMETER (ns2=100)
      PARAMETER (ns3=maxr+5)
c
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      complex phase
c
c
       character*15   os, vo,  con,  du1,  du2
c
      real wff(nr) , grid(nr)
c
      logical there
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(ns3),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
      character ch, ench*11
      ch(i)=char(i+ichar('0'))
c
c      Quantuum numbers
c
       data  in(1), in(2), in(3), in(4) ,in(5) /1, 2, 2, 3, 3/
       data  il(1), il(2), il(3), il(4) ,il(5) /0, 0, 1, 0, 1/
       data  iq(1), iq(2), iq(3), iq(4) ,iq(5) /2, 2, 6, 2, 1/
       
       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 5, 700/
       data  ksu1, mu, dz, kz    /3, 1, 0d0, 0/
       data  os, con, du1, du2 /'a.pri', 'con', 'nul', 'cion'/
       data  vo /'h.pri'/
       ne = nr 
       h = grid(nr) - grid(nr-1)
c

       z = dfloat(nznuc)
       akap = q
       dkap = 0.0 !qstep
       k1   = nc
       k2   = nb 
       il(is) = lv
C  IN(IS) is the principal quantum number of the lowest excited state
       in(is) = 3
       if (lv .gt. 1) in(is) = lv + 1
       L1     = lv
       SS1    = 2
c
       do i = 1, is
          dir(i) = 1
          ud(i)  = 1
C  Don't orthogonalise to the core states
          ud(i)  = 0
       end do 
C  The following patch was put in to get the polarization potential working
C  in the main1 routine for discrete states.
       ud(2) = 0
c
c     Some protection
c
       if (nr .gt. ns3) then
          print*, 'Nr > ns3 in fcz13'
          stop
       end if 
c
c     The Grid
c
       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
c       
       r     = rmax + 1.0
       ne5   = ne + 5
       ne5is = ne5 * is
c
       nu = ns2
       do i = 1, nu
          is2(i) = is
       end do 
c
       include 'fcin.f'
       return
       end

C  File de.f
C  This routine makes and stores the intermediate continuum coulomb states PSIN
      subroutine makepsinc(zasym,nznuc,lna,ce,ncstates,
     >   ps,psen,minps,maxps,gamma,r0)
      INCLUDE 'par.f'
      parameter (ltm = ltmax / 2)
      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
      common/storeps/ psinc(maxrc,ncmax,0:ltm),minpsinc(ncmax,0:ltm)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),minpsinb(nnmax,0:lnabmax)
      common/smallr/ formcut,regcut,expcut,fast
      dimension ps(maxr,ncmax),psen(ncmax),minps(ncmax),maxps(ncmax)
      dimension chi(maxr),ce(ncmax)
      complex phasen
      data lcalc/-1/
      save lcalc
      
C  LCALC is used to keep track of which continuum wave functions have
C  already been calculated. The wave functions will only be saved if
C  LNA = LCALC + 1. This ensures that this routine will be called in 
C  the order LNA = 0, 1, 2, ...
      if (lna .gt. ltm) then
         print*,'Continuum wave functions are only stored up to L:',ltm
         stop 'Need larger LTM here'
      end if 
      jstop=meshr
      do ns = 1, ncstates
         en=ce(ns)
         psen(ns) = ce(ns)
         if (lna.gt.lcalc) then
            q = sqrt(psen(ns))
            if (nznuc.eq.11) then
               rmax = rmesh(meshr,1)
               call fcz11(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >            chi, jstart, 0, lna, nznuc, phasen, gamma, r0)
               jstop = meshr
            elseif (nznuc.eq.3) then
               rmax = rmesh(meshr,1)
               call fcz3(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >            chi, jstart, 0, lna, nznuc, phasen, gamma, r0)
               jstop = meshr
            elseif (nznuc.eq.19) then
               rmax = rmesh(meshr,1)
               call fcz19(rmax, meshr, rmesh, expcut, regcut, q, 1, 
     >            chi, jstart, 0, lna, nznuc, phasen, gamma, r0)
               jstop = meshr
            else
               call coulomb(zasym,en,lna,chi,phasen,jstart,jstop)
            end if 
            minps(ns) = jstart
            maxps(ns) = jstop
            do i = minps(ns), maxps(ns)
               ps(i,ns) = chi(i)
C  I added the following line o 10/12/93 and added some lines to the
C  end of the egrid routine
     >            / q
            end do
            if (lna.eq.-10) then
               do i = minps(ns), maxps(ns)
                  write(60+ns,*) rmesh(i,1), chi(i), lna
               end do
               close(60+ns)
            end if 
            if (lna.le.lnabmax) then
               do nn = lna + 1, min(nnmax,10)
                  sum = 0.0
                  do i = jstart, minpsinb(nn,lna)
                     sum = sum + chi(i) * psinb(i,nn,lna) * rmesh(i,3)
                  end do
                  if (nn.eq.lna+1) then
                     print '(i2,'':'',1p,e7.0,$)',ns,abs(sum)
                  else if (nn.lt.min(nnmax,10)) then
                     print '(1p,e7.0,$)',sum
                  else
                     print '(1p,e7.0)', sum
                  end if
               end do
            end if 
         else
C  read in the saved wave-functions
            minps(ns) = minpsinc(ns,lna)
            maxps(ns) = meshr
            do i = minps(ns), maxps(ns)
               ps(i,ns) = psinc(i,ns,lna)
            end do 
         end if 
      end do

C  Store the calculated wave-function
      if (lna.le.ltm .and. lna.eq.lcalc+1) then
         print*,'saving wave-functions for LNA:',lna
         lcalc = lna
         nwmax = 0
         do ns = 1, ncstates
            minpsinc(ns,lna) = minps(ns)
            do i = minps(ns), maxps(ns)
               psinc(i,ns,lna) = ps(i,ns)
            end do
            wmax = 0.0
            do i = max(minps(ns), maxps(ns) - 100), maxps(ns)
               if (abs(ps(i,ns)).gt.wmax) then
                  wmax=abs(ps(i,ns))
                  imax = i
               end if 
            end do
            if (abs(wmax - sqrt(1.0/acos(0.0))).gt.1e-2)
     >         nwmax = nwmax + 1
         end do
         if (nwmax.gt.0) print*,'WARNING: NWMAX in MAKEPSINC <> 0',
     >      nwmax
      end if 
      end


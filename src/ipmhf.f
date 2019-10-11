C=====Augest,1990===========================================================
C            Analytic Independent-Particle Model for Atoms.
C            Alex E. S. at al, Phys.Rev., v.184, 1969, p.1
C
C            -Z0/r is removed from this potential!
C             -----------------------------------
C
C=====By D.A.Konovalov======================MADE IN AUSTRALIA===============
      subroutine ipmhf (Z, Z0, v, imax, grid, nr, expcut)
      parameter (natom = 19)
      implicit real (a-h, o-z)
      dimension v(nr), grid(nr)
c
c INPUT:
c -----
c  Z     - charge of the atom
c  grid  - r-grid, nr - number of points in this grid
C  expcut- the smollest value of exp-functions in "IPMDHF.F"
C  Z0    - Coulomb potential removed from the potential V.
C          V(r) = Vatom(r)  -  Z0 / r     mult-ed by r if  Z.NE.Z0
c     
c OUTPUT:
c ------
c  V - Array of the values of the potential multiplyed by r  if  Z.NE.Z0
c  imax - the last non-zerro point in the potential
c
      dimension dn(natom)
      data (dn(n), n = 1, natom) /1.0d-30, 0.215e0, 0.563e0, 0.858e0
     >,0.979e0, 0.880e0, 0.776e0, 0.708e0, 0.575e0, 0.500e0, 0.561e0
     >,0.621e0, 0.729e0, 0.817e0, 0.868e0, 0.885e0, 0.881e0, 0.862e0
     >,1.006e0/
      NZ = nint(Z)
c
c     Do same check
c
      if ((NZ .lt. 1) .or. (NZ .gt. Natom)) then
         print*, 'This Z is not available yet,  Z=', Z
         print*, ' Max Z=', Natom
         stop    'Stop in "ipmdhf.f"'
      end if 
c
c     Some constants
c
      dd = dn(NZ)
      hh = dd * (Z - 1e0)**0.4e0
      imax = nr
c
c     Loop by r-grid
c     
      jump = 0
      do i = 1, nr
         r = grid(i)
         r1 = 1e0 / r
         v(i) = Z0 - 1e0
         if (jump .eq. 0) then
            dexpr = exp(- r / dd)
            if (dexpr .lt. expcut)  then
               jump = i
            else 
               teta = dexpr / (hh * (1e0 - dexpr) + dexpr)
               v(i) = - ((Z - 1e0)* teta + 1e0) + Z0 
            end if 
         end if 
      end do 
c
c     Put back 1/r if Z.EQ.Z0
c
      if (nint(Z) .eq. nint(Z0))  then
         do i = 1, imax
            v(i) = v(i) / grid(i)
         end do
      end if 
      return
      end

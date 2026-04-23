c-----------------------------------------------------------------
c*************** Find the eigenvalues and eigenvectors  **********
c********************* of the Hamiltonian matrix *****************
c-----------------------------------------------------------------
      subroutine eigv(enionry,n,H,bb,CI,w,iter)
      use vmat_module, only: nodeid
      include 'par.f'
      double precision sum
      double precision H(n,n),bb(n,n),w(n),
     >   CI(n,n),fv1(n),fv2(n)
      if (n.lt.1) stop 'n < 1 in eigv'
      if (nodeid.eq.1) write(4,'("start diagonalization")')
c      write(4,'("overlap matrix can be printed")')
c      do i=1,n
c         write(4,'(5F12.6)') (real(bb(i,j)), j=1,n)
c      end do

      matz=2
      call  rsg(n,n,H,bb,w,matz,CI,fv1,fv2,ierr)
      if (nodeid.eq.1) then
      write(4,'("ierr =",I3,", ITER =",i3)') ierr, iter
      write(4,'("eigenvalues in a.u.")')
      write(4,'(5F12.6)') (real(w(i)), i=1,n)

c      write(4,'("eigenvalues  - enionry/2. in Ry, enionry=",F10.5,
c     >   " Ry")') enionry
c      write(4,'(5F12.6)') ((real(2d0*w(i)) - enionry), i=1,n)

      write(4,'("eigenvalues in Ry")')
      write(4,'(5F12.6)') (real(2d0*w(i)), i=1,n)

      write(4,'("eigenvalues - enionry/2.  -  in eV, enionry=",F10.5,
     >     " Ry")') enionry
      write(4,'(5F12.6)') ((real(w(i)) - enionry/2.0)*27.2116, i=1,n)
      write(4,'("eigenvectors can be printed")')
c      do i=1,n
c         write(4,'(5F12.6)') (real(CI(i,j)), j=1,n)
c      end do
      endif
      return
c     
c     return
      end
c************************************************************************


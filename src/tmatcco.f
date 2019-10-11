C  File: tmatcco.f
C  The following routine solves the equations T(f,i) = V(f,i) * phi(f) * phi(i)
C  + w(n) * V(f,n) * phi(f) * conjg(phi(n)) * T(n,i) by writing them as
C  (I(f,n)/w(n) - V(f,n)) * (w(n) * conjg(phi(n)) * T(n,i)) = V(f,i) * phi(i).
C  The sum over n is implied. This form has a symmetric kernel.
C  Don't believe the above comment as problems arise when w(n) = 0.
      subroutine gettmat(lg,ns,nchtop,iex,npk,vmat,wk,gk,phasel,
     >   sigma,bb,tdist,von,etot,t2nd,ton,err,vmatop,nchopt,nsmax,ve2ed,
     >   ve2ee,dphasee2e,ephasee2e,nchtope2e,nchmaxe2e,ve2e,te2e,nd,
     >   second,uba,ovlpn,theta,ichi,csfile,projectile,slowery,det,
     >   rcond,vmatp,packed,phaseq)
      include 'par.f'
      parameter (nmax=kmax*nchan)
      integer npk(nchtop+1)
!      pointer (ptrk,kernel)
      complex, allocatable :: kernel(:,:), v(:,:), vdwb2(:,:)
      complex !kernel(ichi,ichi),
     >   dphasee2e(nchane2e),ephasee2e(nchane2e),te2e(nchane2e),
     >   ve2e(nchane2e)!, vdwb2(nchtop,ichi)
!      pointer (ptr2,vdwb2)
      real vmat(ichi,ichi+1),
     >   ve2ed(1,nchmaxe2e,ichi),
     >   ve2ee(1,nchmaxe2e,ichi)
      complex vmatop(kmax,kmax,0:nchanop,nchanop),wk(kmax*nchan)
!      pointer(ptrv,v)
      complex !v(ichi,nchtop),
     >   work(nmax),vns,von(nchan,nchan),
     >   t2nd(nchan,nchan), sigma(nchan),phaseq(knm)
      complex phasel(kmax,nchan),tdist(nchan),ton(nchan,nchan),phase
      dimension gk(kmax,nchan), err(nchan,nchan), ovlpn(knm)
      logical uba(nchan), second,packed
      character csfile*(*), projectile*(*), pfile*80
      
      
      r = float((-1)**ns * min(1,iex))
      r = float((-1)**ns)
c$$$      do m = 1, min(nchan,nchtop*nqm)
c$$$         call makebasis(m,vmat,vmatop,r,nqm,nchtop,nchopt,wk,v,ton)
c$$$         nchi = 1
c$$$         ntop = min(nchtop,3)
c$$$c$$$         if (m.gt.10) then
c$$$            print '(i2,5x,1p,3(2e10.3,1x))',m,
c$$$     >         (ton(nchf,nchi) * phasel(1,nchf) * phasel(1,nchi) /
c$$$     >         gk(1,nchi) / gk(1,nchf),nchf = 1, ntop)
c$$$c$$$         endif
c$$$      enddo
c$$$      print*

      if (ns.eq.0) then
         print'('' JS f i  real(t)   imag(t)    real(t2)  imag(t2)'
     >      //'   real(v)   imag(v) eigenphase'')'
         do nchi = 1, nchtop
            do nchf = nchi, nchtop
               do ki = 1, npk(nchi+1) - npk(nchi)
                  ni = npk(nchi) + ki - 1
                  do kf = 1, npk(nchf+1) - npk(nchf)
                     nf = npk(nchf) + kf - 1
                     if (nf.ge.ni) then
                        vd = (vmat(nf,ni) + vmat(ni,nf+1)) / 2.0
                        ve = (vmat(nf,ni) - vmat(ni,nf+1)) / 2.0
                        vmat(nf,ni) = vd
                        vmat(ni,nf+1) = ve
                     endif 
                  end do
               end do
            end do
         end do
      endif 

c$$$      mv = 8 * (npk(nchtop+1)-1) * nchtop
c$$$      call memalloc(ptrv,mv)
c$$$      if (ptrv.eq.0) stop 'Not enough memory for V'
      allocate(v(npk(nchtop+1)-1, nchtop))
C  Form the driving vector
      nchtopi = nchtop
      do nchi = 1, nchtopi
         call makevec(vmat,vmatop,1,nchi,nqm,npk,nchtop,nchopt,r,
     >      v(1,nchi),ichi)
C  Define the on-shell V, used mainly for printing and the check routine
         do nchf = 1, nchtop
            nf = npk(nchf)
            von(nchf,nchi) = v(nf, nchi)
         end do
      end do 

      if (second) then
c$$$         mdwb2 = 8 * nchtop * ichi
c$$$         call memalloc(ptr2,mdwb2)
c$$$         print'(''Memory (Mb) used for DWB2V:'',f6.1)', mdwb2 * 1e-6
c$$$         if (ptr2.eq.0) stop 'Not enough memory for DWB2V'
         allocate(vdwb2(nchtop, ichi))
         do nchi = 1, nchtop
            do nchf = nchi, nchtop
               do ki = 1, npk(nchi+1) - npk(nchi)
                  ni = npk(nchi) + ki - 1
                  do kf = 1, npk(nchf+1) - npk(nchf)
                     nf = npk(nchf) + kf - 1
                     if (nf.ge.ni) then
                        vns = vmat(nf,ni) + r * vmat(ni,nf+1)
                        if (kf.eq.1) vdwb2(nchf,ni) = - wk(ni) * vns 
                        if (ki.eq.1) vdwb2(nchi,nf) = - wk(nf) * vns
                     endif 
                  end do
               end do
            end do
         end do
         do nchi = 1, nchtop
            do nchf = 1, nchtop
               t2nd(nchf,nchi) = von(nchf,nchi)
               do nn = 1, npk(nchtop+1) - 1
                  t2nd(nchf,nchi) = t2nd(nchf,nchi) -
     >               vdwb2(nchf,nn) * v(nn,nchi) 
               end do
            enddo
         enddo 
         do nchi = 1, nchtop
            do nchf = 1, nchtop
               v(npk(nchf),nchi) = t2nd(nchf,nchi)
            enddo
         enddo 
c$$$         call memfree(ptr2)
      else 
c$$$         mk = 8 * (npk(nchtop+1)-1) ** 2 
c$$$         call memalloc(ptrk,mk)
c$$$         print'(''Memory (Mb) used for KERNEL:'',f6.1)', mk * 1e-6
c$$$         if (ptrk.eq.0) stop 'Not enough memory for KERNEL'
         allocate(kernel(npk(nchtop+1)-1,npk(nchtop+1)-1))
         do nchi = 1, nchopt
            do nchf = nchi, nchopt
               do ki = 1, npk(nchi+1) - npk(nchi)
                  ni = npk(nchi) + ki - 1
                  do kf = 1, npk(nchf+1) - npk(nchf)
                     nf = npk(nchf) + kf - 1
                     if (nf.ge.ni) then
                        vns = vmat(nf,ni)+r*vmat(ni,nf+1) +
     >                     vmatop(kf,ki,nchf,nchi) +
     >                     r * vmatop(kf,ki,nchi-1,nchf)
                        kernel(nf,ni) = - wk(ni) * vns 
                        kernel(ni,nf) = - wk(nf) * vns 
                     endif 
                  end do
               end do
            end do
         end do
C  Ravshan, the following is for you!
         do nchi = 1, nchtop
            do nchf = nchi, nchtop
               if (nchi.gt.nchopt.or.nchf.gt.nchopt) then
                  do ki = 1, npk(nchi+1) - npk(nchi)
                     ni = npk(nchi) + ki - 1
                     do kf = 1, npk(nchf+1) - npk(nchf)
                        nf = npk(nchf) + kf - 1
                        if (nf.ge.ni) then
                           vns = vmat(nf,ni) + r * vmat(ni,nf+1)
                           kernel(nf,ni) = - wk(ni) * vns 
                           kernel(ni,nf) = - wk(nf) * vns
                        endif 
                     end do
                  end do
               endif 
            end do
         end do
C  Get T2nd
         do nchi = 1, nchtopi
            do nchf = 1, nchtop
c$$$            nf = nqm * (nchf - 1) + 1
               nf = npk(nchf)
               t2nd(nchf,nchi) = von(nchf,nchi)
               do nn = 1, npk(nchtop+1) - 1
                  t2nd(nchf,nchi) = t2nd(nchf,nchi) -
     >               kernel(nf,nn) * v(nn,nchi) 
               end do
            end do 
         end do

         epsil = 1e-6
C  Add the I matrix to -K
         do nchf = 1, nchtop
            do kf = 1, npk(nchf+1) - npk(nchf)
               nf = npk(nchf) + kf - 1
               kernel(nf,nf) = kernel(nf,nf)  + (1e0,0e0) 
            end do
         end do 

C  Solve the linear equations
         nv = nchtopi
         nd = npk(nchtop+1) - 1
         call matinv2(kernel,nd,nd,v,nv,work,erfp,epsil)
         call check(vmat,vmatop,nchtop,nchopt,nqm,npk,wk,r,v,von,err,
     >      ichi)
c$$$         call memfree(ptrk)
      endif 


C Calculate photoionization for helium
      if (projectile.eq.'photon') then
         i = 1
         do while (csfile(i:i).ne.'_'.and.csfile(i:i).ne.' ')
            i = i + 1
         enddo 
         pfile = csfile(1:i-8)//'photocs'//csfile(i:80)
         call PHOTO(gk,npk,wk,v,nchtop,ovlpn,pfile,slowery,phasel,
     >      phaseq)
      end if
      
      pi = 3.141592654
C  Get the T-matrix         
      do nchi = 1, nchtopi
         do nchf = 1, nchtop
C  The following loop is a half-on-shell test of the T matrix
c$$$            do n = 1, nqm
c$$$               gkfac = gk(1,nchi) * gk(n,nchf)
c$$$               print*,imag(v(n + nqm * (nchf - 1), nchi)) / gkfac,
c$$$     >            - pi * gk(1,nchf) *
c$$$     >            v(n + nqm * (nchf - 1), nchi) / gkfac *
c$$$     >            conjg(v(1 + nqm * (nchf - 1), nchi)) /
c$$$     >            gk(1,nchi) / gk(1,nchf)
c$$$            enddo 
            gkfac = gk(1,nchi) * gk(1,nchf)
            phase = phasel(1,nchf) * phasel(1,nchi)
            ton(nchf,nchi)  = phase * v(npk(nchf),nchi) / gkfac 
            t2nd(nchf,nchi) = phase * t2nd(nchf,nchi) / gkfac 
            von(nchf,nchi)  = phase * von(nchf,nchi) / gkfac 
         end do
C  Add to the T matrix the part that comes from the distorting potential.
C  This applies only to the diagonal parts.
C  In the case of Spin = 1.5 we musn't add TDIST to those TMATs which are 0
         if (abs(ton(nchi,nchi)).gt.0.0) then
            ton(nchi,nchi)  = ton(nchi,nchi)  + tdist(nchi)
            t2nd(nchi,nchi) = t2nd(nchi,nchi) + tdist(nchi)
            von(nchi,nchi)  = von(nchi,nchi)  + tdist(nchi)
         endif 
c$$$            print'(1p,10(2e10.3))',(von(nchf,nchi),nchf=1,nchtop)
c$$$            print'(1p,10(2e10.3))',(ton(nchf,nchi),nchf=1,nchtop)
c$$$            print*
      end do
      if (nchtope2e.gt.0.and.second) print*,'Won''t get correct',
     >   ' e2e T matrix element as not all of T has been defined'
      do nchi = 1, 1
         ni = npk(nchi)
         do nchf = 1, nchtope2e
            ve2e(nchf) = ve2ed(1,nchf,ni) * dphasee2e(nchf) +
     >         r * ve2ee(1,nchf,ni) * ephasee2e(nchf)
            te2e(nchf) = ve2e(nchf)
            nn = 0
            do nchn = 1, nchtop
               do n = 1, npk(nchn+1) - npk(nchn)
                  nn = nn + 1
                  te2e(nchf) = te2e(nchf) +
     >               (ve2ed(1,nchf,nn) * dphasee2e(nchf)
     >               + r * ve2ee(1,nchf,nn) * ephasee2e(nchf))
     >               * wk(nn) * v(nn, nchi)
               enddo
            enddo
            te2e(nchf) = te2e(nchf) * phasel(1,nchi) / gk(1,nchi)
            ve2e(nchf) = ve2e(nchf) * phasel(1,nchi) / gk(1,nchi)
         enddo
      enddo      
      return
      end
      
      subroutine check(vmat,vmatop,nchtop,nchopt,nqm,npk,wk,r,v,von,err,
     >   ichi)
      include 'par.f'
      dimension npk(nchtop+1)
      real vmat(ichi,ichi+1)
      complex vmatop(kmax,kmax,0:nchanop,nchanop)
      complex v(npk(nchtop+1)-1,nchtop),
     >   wk(kmax*nchan),von(nchan,nchan),vec(kmax*nchan),sum
      real err(nchan,nchan)
      
      do nchf = 1, nchtop
         call makevec(vmat,vmatop,1,nchf,nqm,npk,nchtop,nchopt,r,vec,
     >      ichi)
         do nchi = 1, nchtop
            sum = (0.0,0.0)
            do n = 1, npk(nchtop+1) - 1
               sum = sum + vec(n) * v(n,nchi) * wk(n)
            enddo
            nt = npk(nchf)
            err(nchf,nchi) = abs(von(nchf,nchi) + sum - v(nt,nchi))
         enddo
      enddo
      return
      end
      
      subroutine makevec(vmat,vmatop,kf,nchf,nqm,npk,
     >   nchtop,nchopt,r,vec,ichi)
      include 'par.f'
      dimension npk(nchtop+1)
      real vmat(ichi,ichi+1)
      complex vmatop(kmax,kmax,0:nchanop,nchanop)
      complex vec(ichi)

      nf = npk(nchf) + kf - 1
      do nn = 1, nf
         vec(nn) = vmat(nf,nn) + r * vmat(nn,nf+1)
      enddo 

      do nn = nf + 1, npk(nchtop+1) - 1
         vec(nn) = vmat(nn,nf) +  r * vmat(nf,nn+1) 
      end do
               
      do nchn = 1, min(nchf,nchopt)
         do kn = 1,  npk(nchn+1) - npk(nchn)
            nn = npk(nchn) + kn - 1
            vec(nn) = vec(nn) + vmatop(kf,kn,nchf,nchn) + 
     >         r * vmatop(kf,kn,nchn-1,nchf)
         end do
      end do
      do nchn = nchf+1, nchopt
         do kn = 1, npk(nchn+1) - npk(nchn)
            nn = npk(nchn) + kn - 1
            vec(nn) = vec(nn) + vmatop(kn,kf,nchn,nchf) + 
     >         r * vmatop(kn,kf,nchf-1,nchn) 
         end do
      end do
      end

c$$$      subroutine makebasis(m,vmat,vmatop,r,nqm,nchtop,nchopt,wk,z,ton)
c$$$      include 'par.f'
c$$$      real vmat(nqm,nchtop,nqm,nchtop+1)
c$$$c$$$      real vmat(kmax,kmax,0:nchan,nchan)
c$$$      complex vmatop(kmax,kmax,0:nchanop,nchanop)
c$$$      complex z(kmax*nchan,nchan),
c$$$     >   wk(kmax,nchan), vec(kmax*nchan), x(kmax*nchan,nchan),
c$$$     >   mat(nchan,nchan), b(nchan), kernel(nchan,nchan),
c$$$     >   ton(nchan,nchan), bmat(nchan), c
c$$$      integer ipiv(nchan)
c$$$      save mat, x, b
c$$$
c$$$      nd = nqm * nchtop
c$$$      if (m.eq.1) then
c$$$         kf = 1
c$$$         nchf = 1
c$$$         call makevec(vmat,vmatop,kf,nchf,nqm,nchtop,nchopt,r,vec,ichi)
c$$$      else
c$$$         do n = 1, nd
c$$$            vec(n) = x(n,m-1)
c$$$         enddo 
c$$$         do j = 1, m-1
c$$$            c = (0.0,0.0)
c$$$            do np = 1, nd
c$$$               c = c + conjg(z(np,j)) * x(np,m-1)
c$$$            enddo
c$$$            do n = 1, nd
c$$$              vec(n) = vec(n) - z(n,j) * c
c$$$            enddo
c$$$         enddo
c$$$      endif 
c$$$
c$$$      veclen2 = 0.0
c$$$      do n = 1, nd
c$$$         veclen2 = veclen2 + vec(n) * conjg(vec(n))
c$$$      enddo
c$$$      veclen = sqrt(veclen2)
c$$$      do n = 1, nd
c$$$         z(n,m) = vec(n) / veclen
c$$$      enddo
c$$$      
c$$$      n = 0
c$$$      do nchf = 1, nchtop
c$$$         do kf = 1, nqm
c$$$            call makevec(vmat,vmatop,kf,nchf,nqm,nchtop,nchopt,r,vec,ichi)
c$$$            if (kf.eq.1.and.nchf.eq.1) then
c$$$               b(m) = (0.0,0.0)
c$$$               do np = 1, nd
c$$$                  b(m) = b(m) + conjg(z(np,m)) * vec(np)
c$$$               enddo
c$$$            endif 
c$$$            n = n + 1
c$$$            x(n,m) = (0.0,0.0)
c$$$            np = 0
c$$$            do nchn = 1, nchtop
c$$$               do kn = 1, nqm
c$$$                  np = np + 1
c$$$                  x(n,m) = x(n,m) + wk(kn,nchn) * vec(np) * z(np,m)
c$$$               enddo
c$$$            enddo 
c$$$         enddo
c$$$      enddo
c$$$      
c$$$      do mi = 1, m
c$$$         do mf = 1, m
c$$$            if (mf.eq.m.or.mi.eq.m) then
c$$$               mat(mf,mi) = (0.0,0.0)
c$$$               do n = 1, nd
c$$$                  mat(mf,mi) = mat(mf,mi) + conjg(z(n,mf)) * x(n,mi)
c$$$               enddo
c$$$            endif
c$$$            kernel(mf,mi) = - mat(mf,mi)
c$$$         enddo
c$$$         kernel(mi,mi) = kernel(mi,mi) + (1.0,0.0)
c$$$         bmat(mi) = b(mi)
c$$$      enddo
c$$$      nb = 1
c$$$      lda = nchan
c$$$      call cgesv(m,nb,kernel,lda,ipiv,bmat,lda,info)
c$$$      nchi = 1
c$$$      nn = 0
c$$$      do nchf = 1, nchtop
c$$$         ton(nchf,nchi) = (0.0,0.0)
c$$$         np = 1 + nqm * (nchf - 1)
c$$$         do j = 1, m
c$$$            ton(nchf,nchi) = ton(nchf,nchi) + bmat(j) * z(np,j)
c$$$         enddo
c$$$      enddo
c$$$
c$$$c$$$      do mi = 1, m
c$$$c$$$         do mf = 1, mi
c$$$c$$$            c = (0.0,0.0)
c$$$c$$$            do n = 1, nd
c$$$c$$$               c = c + z(n,mi) * conjg(z(n,mf))
c$$$c$$$            enddo
c$$$c$$$            print*,'mf,mi,ovlp',mf,mi,c
c$$$c$$$         enddo
c$$$c$$$      enddo 
c$$$      return
c$$$      end
c$$$      
      
         
         
      
      subroutine optpot
C  This routine is not used in CCC, and so is not loaded. For CCO the routine
C  is in the file de.f
      end
      
       
      SUBROUTINE MATINV2(API,N1,NMAX,BHAR,nb,WK,ERFP,EPSIL)
      include 'par.f'
      COMPLEX API,BHAR,WK,DET(2)
      REAL ERFP,EPSIL
      DIMENSION API(N1,N1),BHAR(N1,nchan),WK(N1),kpvt(kmax*nchan)
      
      lda = n1
      n = nmax
      call cgesv(n,nb,api,lda,kpvt,bhar,lda,info)
c$$$      call cgeco(api,lda,n,kpvt,rcond,wk)
c$$$c$$$      call cgefa(api,lda,n,kpvt,info)
c$$$      do m = 1, nb
c$$$         call cgesl(api,lda,n,kpvt,bhar(1,m),0)
c$$$      enddo 
c$$$      call cgedi(api,lda,n,kpvt,det,wk,10)
C  The following Linpack routines assume that API is a symmetric kernel
c$$$      call csifa(api,lda,n,kpvt,info)
c$$$      call csisl(api,lda,n,kpvt,bhar)
c$$$      call csidi(api,lda,n,kpvt,det,wk,10)
      if (info.ne.0) then
         print*,'STOPPING since INFO in MATINV2:',info
         stop 'STOPPING since INFO in MATINV2 is nonzero'
      endif 
      return
C  The IMSL routine assumes a general kernel
c$$$      IA=N1
c$$$      IB=N1
c$$$      M=1
c$$$      IJOB=0
c$$$      CALL CLEQT(API,NMAX,IA,BHAR,M,IB,IJOB,WK,IER)
c$$$      DET(1) = (1.0,0.0)
c$$$      DET(2) = (0.0,0.0)
c$$$      DO 5 I = 1,NMAX
c$$$         IPVT = WK(I)
c$$$         IF (IPVT .NE. I) DET(1) = -DET(1)
c$$$         DET(1) = DET(1)*API(I,I)
c$$$         do while (abs(det(1)).gt.1e1)
c$$$            det(1) = det(1) / 1e1
c$$$            det(2) = det(2) + (1.0,0.0)
c$$$         enddo 
c$$$         do while (abs(det(1)).lt.1.0)
c$$$            det(1) = det(1) * 1e1
c$$$            det(2) = det(2) - (1.0,0.0)
c$$$         enddo 
c$$$ 5    CONTINUE
c$$$      if (ier.ne.0) print*,'IER not zero in CLEQT',ier
c$$$      IF (IER.EQ.129) THEN
c$$$         PRINT*,'MATRIX IS ALGORITHMICALLY SINGULAR.'
c$$$         RETURN
c$$$      END IF
c$$$      IF (IER.EQ.130) PRINT*,'MATRIX IS ILL CONDITIONED.'
      END

c      
c        
      subroutine redefine(nold,ps,psen,minps,maxps,nps,
     >   psinb,enpsinb,istoppsinb,meshr,l,nabot,vdcore,gridr)
      include 'par.f'
      real vdcore(maxr), gridr(maxr,3)
      dimension ps(maxr,ncmax), psen(ncmax), minps(ncmax), maxps(ncmax),
     >   psinb(maxr,nnmax), enpsinb(nnmax), istoppsinb(nnmax),
     >   ps_tmp(nnmax)
      integer  nabot(0:lamax)
      double precision H(nps,nps),bb(nps,nps),w(nps),
     >     CI(nps,nps), r1elk, ovlp, tmp
      double precision fv1(nps),fv2(nps)
c      real G(nps,nps), gl(meshr)
      
c     nold is maximum princ. quant. number of the orbital that is replaced
c     by true one-electron atom bound state (all orbitals with n <= nold are replaced).

      if(l .eq. 0) then
         open(538,file='tmp.ovlp.ps')
      end if
      do n1=1,nps
         do n2=n1,nps
            ovlp =  r1elk(0,ps(1,n1),ps(1,n2),
     >         minps(n1),minps(n2),maxps(n1),maxps(n2))
            write(538,'(3I4,1p,E15.6)') l, n1,n2,ovlp
         end do
      end do
      if(l .eq. 3) then
         close(538)
      end if
     
c     Determine how many functions must be redefined:
      Nnew = min(nold-l,nps)
      if(Nnew .le. 0)  return
      do m = 1, Nnew
         n = m + l
         psen(m) = enpsinb(n) 
         minps(m) = 1
         maxps(m) = istoppsinb(n)
         ps(1:maxps(m),m) = psinb(1:maxps(m),n)
      enddo
      
      
      
      H(:,:) = 0d0
      bb(:,:) = 0d0
      CI(:,:) = 0d0

c     G(:,:) = 0.0

C     First Nnew orbitals are (approximate) eigenfunctions of the HF or FCHF
c     Hamiltonians:  H |n1> = psen(n1) |n1>
c     Note that HF and FCHF orbitals can be nonorthognal, and it would be
c     better to explisitly calculate <n1|H|n2>, but what is coded should be
c     accurate enough.
c     All other orbitals diagonalize FCHF Hamiltonian: <n1|H|n2> = psen(n1) delta(n1,n2)
      open(536,file='tmp.ovlp')
      do n1=1,Nnew      ! nps
         H(n1,n1) = psen(n1) 
         bb(n1,n1) = 1d0
         ovlp =  r1elk(0,ps(1,n1),ps(1,n1),
     >      minps(n1),minps(n1),maxps(n1),maxps(n1))
         write(536,'(3I4,1p,E15.6)') l, n1,n1,ovlp
         do n2=n1+1,nps         ! It is safer, because first Nnew orbitals can be 
                                ! from both hf and fch, and they are not orthogonal.
c     do n2=Nnew+1,nps
            ovlp =  r1elk(0,ps(1,n1),ps(1,n2),
     >         minps(n1),minps(n2),maxps(n1),maxps(n2))
            write(536,'(3I4,E15.6)') l, n1,n2,ovlp
            if(abs(ovlp) .lt. 5d-5) ovlp = 0d0
            H(n1,n2) = dble(psen(n1)) * ovlp
            H(n2,n1) = H(n1,n2)

c            il = max(minps(n1),minps(n2))
c            im = min(maxps(n1),maxps(n2))
c            call  Ham1el(meshr,l,vdcore(1:nr),il,im,ps(1:nr,n2),
c     >         ps(1:nr,n1),gl,gridr(1:nr,1),gridr(1:nr,3),res)
c            G(n1,n2) = 2.0*res
c            G(n2,n1) = G(n1,n2)
c            print*, n1,n2, H(n1,n2), G(n1,n2)
            bb(n1,n2) = ovlp
            bb(n2,n1) = bb(n1,n2)
         end do 
      end do

c      H(Nnew+1:nps,Nnew+1:nps) = 0.0

      do n1=Nnew+1, nps         
         H(n1,n1) = dble(psen(n1))
         bb(n1,n1) = 1d0
      end do
      
      matz=2
      call  rsg(nps,nps,H,bb,w,matz,CI,fv1,fv2,ierr)
      print*, 'ierr =', ierr, ', replaced n=',Nnew, '  lowest orbitals'
      
      psen(1:nps) = w(1:nps)

      i_max=0
      do n=1,nps
         i_max = max(i_max,maxps(n))
      end do
      do i=1,i_max              ! meshr
c     Core orbitals may not be overwritten if after this we go to 
c     scattering calculations. But the way the core orbitals are calculated 
c     ensures that they are the same after the diagonalization in this routine.
         do n=1,nps
            tmp = 0d0
            do m=1,nps
               if(minps(m).le.i.and.i.le.maxps(m)) then
                  tmp = tmp + CI(m,n)*dble(ps(i,m))
               end if
            end do
            ps_tmp(n) = tmp
         end do
         ps(i,1:nps) = ps_tmp(1:nps)
      end do
      do n=1,nps
         call minmaxi(ps(1,n),meshr,i1,i2)
         maxps(n) = i2
         minps(n) = i1
      end do
      
c     Fix sign of the one-electron functions
      do n=1,nps         
         if( sign(1.0,ps(minps(n)+1,n)) .lt. 0)  then
            do i=minps(n),maxps(n) 
               ps(i,n) = -ps(i,n)
            end do
         end if
      end do

     
      return 
      end
      
c----------------------------------------------------------------------------------
c     energy is in Ry.
      subroutine redefine_rel(ps,psen,minps,maxps,nps,lsp,nznuc,nabot)
      include 'par.f'
      common /meshrr/ meshr, gridr(maxr,3)
      real  ps(maxr,ncmax), psen(ncmax)
      integer  minps(ncmax), maxps(ncmax)
      integer  nps, lsp
      integer  nabot(0:lamax)
      real  ps_tmp(nps)
      double precision H(nps,nps),bb(nps,nps),w(nps),
     >   CI(nps,nps), tmp
      double precision fv1(nps),fv2(nps)
      double precision  Rel_cor   
      common /relcore/ pscore(3000,10,0:lamax), enpscore(10,0:lamax),
     >   istopcore(10,0:lamax)
      common/smallr/ formcut,regcut,expcut,fast,match
      double precision  CR(meshr), A, B
      integer iCR
      real tmppsen(nps)
      common /Hrel/ Hrel(10,10,0:lamax), nrel(0:lamax)
c     
      H = 0d0
      bb = 0d0
      CI = 0d0

      tmppsen(1:nps) = psen(1:nps)
      
c     <n2|H|n1>
      knmax = nabot(lsp) -lsp + 2
c      do n1=1,nps       
      do n1=1,knmax ! nps       
         do n2=1, n1
c            print*,n1,n2
            im = min(maxps(n1),maxps(n2),iCR)
c            H(n1,n2) = SUM(ps(1:im,n2)*ps(1:im,n1)*
c     >         gridr(1:im,3)*CR(1:im))
            H(n1,n2) = Rel_cor(meshr,ps(1:meshr,n2),ps(1:meshr,n1),
     >      lsp,minps(n2),maxps(n2),minps(n1),maxps(n1),
     >      gridr(1:meshr,1),gridr(1:meshr,3),nznuc) * 2d0 ! Rel_cor is in au but H is in Ry, so we multiply by 2.
            H(n2,n1) = H(n1,n2)
c            write(*,'("n1,n2,H(n1,n2)",2i5,e15.5)') n1,n2,H(n1,n2)
         end do 
      end do
      do n1=1,nps       
         H(n1,n1) =  H(n1,n1) +  dble(psen(n1))
         bb(n1,n1) = 1d0
      end do
      
      matz=2
      call  rsg(nps,nps,H,bb,w,matz,CI,fv1,fv2,ierr)
      print*, 'redefine_rel:  ierr =', ierr
      

      psen(1:nps) = w(1:nps)

      i_max=0
      do n=1,nps
         i_max = max(i_max,maxps(n))
      end do
      do i=1,i_max              ! meshr
         do n=1,nps
            tmp = 0d0
            do m=1,nps
               if(minps(m).le.i.and.i.le.maxps(m)) then
                  tmp = tmp + CI(m,n)*dble(ps(i,m))
               end if
            end do
            ps_tmp(n) = tmp
         end do
         ps(i,1:nps) = ps_tmp(1:nps)
      end do
      do n=1,nps
         call minmaxi(ps(1,n),meshr,i1,i2)
         maxps(n) = i2
         minps(n) = i1
      end do
      
c     Fix sign of the one-electron functions
      do n=1,nps         
         if( sign(1.0,ps(minps(n)+1,n)) .lt. 0)  then
            do i=minps(n),maxps(n) 
               ps(i,n) = -ps(i,n)
            end do
         end if
      end do
      print*,'Core wave-functions after account of Rel. Corrections'
      do k=1,6
         print'(2i3,i6,12x,f16.8,"  eV, ",f16.8,"  au")',
     >      lsp,k,maxps(k),w(k)*13.6058, w(k)/2.0
      end do

c$$$c     Save core orbitals for iterative improvement of orbitals
c$$$      kn = nabot(lsp) - lsp  - 1    ! last core orbital
c$$$      do m = 1, kn
c$$$         enpscore(m,lsp)  = psen(m) 
c$$$         istopcore(m,lsp) = maxps(m) 
c$$$         pscore(1:maxps(m),m,lsp) = ps(1:maxps(m),m) 
c$$$      enddo
c     Save orbitals for use in scattering calculations
      do m = 1, knmax
         enpscore(m,lsp)  = psen(m) 
         istopcore(m,lsp) = maxps(m)
         if(maxps(m) .gt. 3000) then
            print*,'Increase size of array ps(i,:,:) ',
     >         'in common /relcore/, to :', maxps(m)
            stop
         end if
            
         pscore(1:maxps(m),m,lsp) = ps(1:maxps(m),m) 
      enddo

c     get matrix elements of rel. corrections in for the newly calc. functions.
c
      do n1=1,nps       
         H(n1,n1) =  H(n1,n1) -  tmppsen(n1)
      end do
      Hrel(:,:,lsp) = 0.0
      do n1=1,knmax
         do n2=1,n1
            tmp = 0.0
            do i=1,knmax
               do j=1,knmax
                  Hrel(n1,n2,lsp) = Hrel(n1,n2,lsp) +
     >               H(i,j)*CI(i,n1)*CI(j,n2)
               end do
               tmp = tmp - CI(i,n1)*CI(i,n2)*tmppsen(i)
            end do
            if(n1.eq.n2) tmp = tmp + psen(n1)
            Hrel(n2,n1,lsp) = Hrel(n1,n2,lsp)
            write(*,'(3i5,1p,3e15.6)') lsp,n1,n2,
     >         Hrel(n1,n2,lsp),H(n1,n2), tmp
         end do
      end do
      nrel(lsp) = knmax

      return 
      end


c--------------------------------------------------------------------------
c     Mass term: (p^2)^2, Darwin term Z delta(r), atomic units
      double precision function  Rel_cor(nr,psf,psi,lsp,
     >   minf,maxf,mini,maxi,gridr,weightr,nznuc)
      real  psf(nr),  psi(nr)
      double precision, dimension(nr) :: fun_f, fun_i, fun !, firstderiv
      real, dimension(nr) :: gridr, weightr
      double precision  tmp  !, tmp1
      double precision  alpha
      data alpha / 7.29714d-3 /    ! 1d0/137.04d0
      
c     Mass term
      fun = 0d0; fun_f = 0d0; fun_i = 0d0
      
      minfun=max(mini,minf)
      maxfun=min(maxi,maxf)

      call secondderiv(psf(:), minfun, maxfun, gridr(:), nr, fun_f)
      call secondderiv(psi(:), minfun, maxfun, gridr(:), nr, fun_i)
c      open(128,file='file.tmp')
c      do i=1,300
c         write(128,'(10E15.5)') gridr(i),psf(i),fun_f(i),
c     >      psi(i),fun_i(i)
c      end do
c      close(128)
      if(lsp .ne. 0) then
         fun_f(minfun:maxfun) = fun_f(minfun:maxfun) -
     >      dble(lsp*(lsp + 1)) *dble(psf(minfun:maxfun)/
     >      gridr(minfun:maxfun)/gridr(minfun:maxfun))
         fun_i(minfun:maxfun) = fun_i(minfun:maxfun) -
     >      dble(lsp*(lsp + 1)) *dble(psi(minfun:maxfun)/
     >      gridr(minfun:maxfun)/gridr(minfun:maxfun))
      end if
      fun(minfun:maxfun) = fun_f(minfun:maxfun)*fun_i(minfun:maxfun)*
     >   dble(weightr(minfun:maxfun)) 

      tmp = ((psf(2)-2.0*psf(1))/gridr(1)/gridr(1)) *
     >   ((psi(2)-2.0*psi(1))/gridr(1)/gridr(1)) *
     >   weightr(1)/4.0
      Rel_cor = -(SUM(fun) + tmp)
c
      
c     Darwin term
      tmp = 0.0
      if(lsp .eq. 0) then
c     See  notes in the notebook for derivation.
         tmp = dble(nznuc)*((4d0*psf(1)-psf(2))/2d0/gridr(1))*
     >      ((4d0*psi(1)-psi(2))/2d0/gridr(1))

c         call f_deriv(psi,mini,maxi,gridr,nr,1,firstderiv)
         
c         fun(minfun:maxfun)=(psf(minfun:maxfun)/gridr(minfun:maxfun))
c     >      * (firstderiv(minfun:maxfun) -
c     >      psi(minfun:maxfun)/gridr(minfun:maxfun))
c     >      /gridr(minfun:maxfun)
         
c         tmp1 = -dble(nznuc)*2.0*SUM(fun(minfun:maxfun)
c     >      *dble(weightr(minfun:maxfun)) )         
c         tmp2 = dble(nznuc)*(psi(1)/gridr(1))*(psf(1)/gridr(1))         
c         print*, tmp, tmp1, tmp/tmp1,  tmp2
c         tmp = tmp1
      end if

c      print*,  lsp, real(Rel_cor*alpha*alpha/8d0),
c     >   real(tmp*alpha*alpha/8d0),
c     >   real((Rel_cor+tmp)*alpha*alpha/8d0)

      Rel_cor =  (Rel_cor  + tmp) * alpha*alpha/8d0

      return
      end  function Rel_cor
c--------------------------------------------------------------------------
      subroutine  secondderiv(f, minf, maxf, grid, nr, result)
      real, intent(in) :: f(nr)
      integer minf, maxf
      real, intent(in) ::  grid(nr)
      integer, intent(in) :: nr
      double precision, intent(out) :: result(nr)
      double precision  h1, h2, h, hh
      
      result(:) = 0d0
      
      
c     print*,' get   Second derivative' 
      if(minf.eq.1) then 
         h = dble(grid(2)) - dble(grid(1))
         result(1) = (dble(f(2)) - 2d0*dble(f(1)))/h/h ! note - value of fl(r,n)  at r=0.0 is zero.
      end if
      
      do i=max(2,minf),min(nr,maxf-1)
         h2 = dble(grid(i+1) - grid(i))
         h1 = dble(grid(i) - grid(i-1))
         if(abs(h2-h1) .le. 1e-4) then
            result(i) = (dble(f(i+1)) - 2d0*dble(f(i)) +
     >         dble(f(i-1)))/h1/h1
         else
            hh = h1 + h2
c            result(i) = (dble(f(i+1)) / h2 + dble(f(i-1)) / h1) *
c     >         2d0/(h1+h2)
c     >         - dble(f(i)) * 2d0 / h1 / h2
            result(i) = 2d0*( dble(f(i+1)) + dble(f(i-1))*(h2/h1)
     >         - dble(f(i))*(hh/h1))/ h2 / hh
         end if
      end do      
      
      
      return
      end  subroutine  secondderiv
      
c--------------------------------------------------------------------------

      subroutine relcorH12(l,fl,minf,maxf,chil,minc,maxc,result)
      include 'par.f'
      real fl(nr), chil(nr)
      real result, tmp
      real tc(10), tf(10)
      common /relcore/ ps(3000,10,0:lamax), enps(10,0:lamax),
     >   istop(10,0:lamax)
      common /Hrel/ Hrel(10,10,0:lamax), nrel(0:lamax)
      common /meshrr/ nr, gridr(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)

      kn = nrel(l)

      do n=1,kn
         in = min(istop(n,l),maxc)
         im = max(1,minc)
         tc(n) = SUM(ps(im:in,n,l)*chil(im:in))
         if(n .ge. nabot(l)-l) then
            tf(n) = SUM(ps(im:in,n,l)*fl(im:in)*gridr(im:in,3))
         else
            tf(n) = 0.0
         end if
      end do

      tmp = 0.0

      do n=nabot(l)-l,kn
         tmp = tmp + SUM(tc(1:kn)*Hrel(1:kn,n,l))*tf(n)
      end do
      
      result = result + tmp
      
      return
      end
      
c----------------------------------------------------------------------------------
c----------------------------------------------------------------------------------
c$$$c     Inclusion of relativistic mass and Darwin terms (step 2) will 
c$$$c     change pseudostates wave functions. They are not orthogonal to
c$$$c     nonrelativistic core w.f. anymore. This can be corrected (but not 100%)
c$$$c     by diagonalizing H(core w.f.) + H_{rel}, where in H() we use
c$$$c     core w.f. from diagonalization which inludes rel.corrections (step 2)
c$$$c     Note that the core w.f.  and pseudostate w.f. are still not orthogonal
c$$$c     but to a much less degree (at least rel.corrections are included in both).
c$$$c      
c$$$c     energy is in Ry.
c$$$      subroutine iterimp_rel(vdcore,corep,r0)
c$$$      include 'par.f'
c$$$      common/meshrr/ meshr,rmesh(maxr,3)
c$$$      real  vdcore(maxr,0:lamax), corep(0:lamax), r0(0:lamax)
c$$$      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
c$$$     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
c$$$     >   npbot(0:lamax),nptop(0:lamax)
c$$$      real  temp(maxr), uplane(maxr)
c$$$      real vdnewcore(meshr,0:latop)
c$$$      common /psinbc/ enpsinb(nnmax,0:lnabmax),
c$$$     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
c$$$      common /relcore/ pscore(2000,5,0:3), enpscore(5,0:3),
c$$$     >   istopcore(5,0:3)
c$$$      common /newrelcore/ pscoren(2000,5,0:3), enpscoren(5,0:3),
c$$$     >   istopcoren(5,0:3)
c$$$
c$$$      iter = 8
c$$$      print*,'***  Start iterimp_rel(): will do  ', iter,' iterations'
c$$$ 10   continue
c$$$      
c$$$      enpscoren =  0.0
c$$$      istopcoren =  0.0
c$$$      pscoren = 0.0
c$$$
c$$$
c$$$      
c$$$C     Get NEW direct core potential.
c$$$C     VDCORE will contain the direct plus l-dependent core polarization potential
c$$$      temp(1:meshr) = 0.0
c$$$      uplane(1:meshr) = 0.0
c$$$c
c$$$c     It uses  /relcore/ arrays
c$$$      call makenewvdcore(temp,minvdc,maxvdc,nznuc,uplane)  
c$$$      
c$$$      do la = labot, latop
c$$$         if (corep(la).gt.0.0) then
c$$$            minvdc = 1
c$$$            maxvdc = meshr
c$$$         endif 
c$$$         do i = 1, meshr
c$$$C     Phenomenological polarization potential at rmesh(i,1) in a.u. is polpotr
c$$$            polpotr = polpot(rmesh(i,1),corep(la),r0(la))
c$$$            vdnewcore(i,la) = temp(i) + polpotr
c$$$         enddo 
c$$$      enddo
c$$$      
c$$$      do la = labot, latop
c$$$         nps = natop(la) - nabot(la) + 1         
c$$$         call  iterimp_rel_la(la,nps,vdnewcore(:,la),vdcore(:,la),iter)
c$$$      end do
c$$$      
c$$$c     update vdcore array:
c$$$      vdcore(:,:) = vdnewcore(:,:)
c$$$      
c$$$c     update psinb - core arrays - to use FC exchange routine correctly:      
c$$$      do la = labot, latop
c$$$         do m = 1,nabot(la)-la-1
c$$$            n = m + la
c$$$            enpsinb(n,la) = enpscore(m,la)
c$$$            ist = istopcore(m,la)
c$$$            istoppsinb(n,la) = ist
c$$$            psinb(1:ist,n,la) = pscore(1:ist,m,la)
c$$$         enddo      
c$$$      enddo      
c$$$
c$$$      enpscore =  enpscoren
c$$$      istopcore =  istopcoren
c$$$      pscore = pscoren
c$$$
c$$$      iter = iter - 1
c$$$      if(iter .gt. 0)  goto 10
c$$$      
c$$$      return 
c$$$      end  subroutine iterimp_rel
c$$$c----------------------------------------------------------------------------------
c$$$c     energy is in Ry.
c$$$      subroutine iterimp_rel_la(la,nps,vdnewcore,vdcore,iter)
c$$$      include 'par.f'
c$$$      common /newrelcore/ pscoren(2000,5,0:3), enpscoren(5,0:3),
c$$$     >   istopcoren(5,0:3)
c$$$      common /relcore/ pscore(2000,5,0:3), enpscore(5,0:3),
c$$$     >   istopcore(5,0:3)
c$$$      common/meshrr/ meshr,rmesh(maxr,3)
c$$$      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
c$$$     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
c$$$     >   npbot(0:lamax),nptop(0:lamax)
c$$$      common /psinbc/ enpsinb(nnmax,0:lnabmax),
c$$$     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
c$$$      real vdnewcore(meshr)
c$$$      real  vdcore(maxr)
c$$$      real  ps(meshr,nps), psen(nps)
c$$$      integer  minps(nps), maxps(nps), nps, la
c$$$      real  ps_tmp(nps)
c$$$      real f1(meshr), f2(meshr)
c$$$      double precision H(nps,nps),bb(nps,nps),w(nps),
c$$$     >   CI(nps,nps), tmp
c$$$      double precision fv1(nps),fv2(nps)
c$$$c
c$$$      ps = 0.0
c$$$      psen = 0.0
c$$$c     Copy to local arrays:
c$$$c     First core orbitals
c$$$      do m = 1, nabot(la)-la-1
c$$$         psen(m) = enpscore(m,la) 
c$$$         minps(m) = 1
c$$$         maxps(m) = istopcore(m,la)
c$$$         ps(1:maxps(m),m) = pscore(1:maxps(m),m,la)
c$$$      enddo
c$$$c     And second all other orbitals      
c$$$      do m = nabot(la)-la, nps
c$$$         n = m + la
c$$$         psen(m) = enpsinb(n,la) 
c$$$         minps(m) = 1
c$$$         maxps(m) = istoppsinb(n,la)
c$$$         ps(1:maxps(m),m) = psinb(1:maxps(m),n,la)
c$$$      enddo
c$$$      
c$$$      
c$$$      H = 0d0
c$$$      bb = 0d0
c$$$      CI = 0d0
c$$$
c$$$c     <n2|H|n1>
c$$$      do n1=1,nps
c$$$         f1(1:meshr) = ps(1:meshr,n1)*rmesh(1:meshr,3)
c$$$         do n2=1,n1
c$$$            f2(1:meshr) = ps(1:meshr,n2)*rmesh(1:meshr,3)
c$$$            il = max(minps(n1),minps(n2))
c$$$            im = min(maxps(n1),maxps(n2))
c$$$c     It uses psinb() array for core functions
c$$$            call tmpfcexch(f1,minps(n1),maxps(n1),meshr,f2,res,la)
c$$$c     It uses pscore() array for core functions
c$$$            call newtmpfcexch(f1,minps(n1),maxps(n1),meshr,
c$$$     >         f2,resnew,la)
c$$$            f2(il:im) = ps(il:im,n1)*(vdnewcore(il:im) -
c$$$     >         vdcore(il:im))*ps(il:im,n2)*rmesh(1:meshr,3)
c$$$            tmp = SUM(f2(il:im))
c$$$            H(n1,n2) = 2d0*(tmp + resnew - res)
c$$$            H(n2,n1) = H(n1,n2)
c$$$         end do 
c$$$      end do
c$$$      do n1=1,nps       
c$$$         H(n1,n1) =  H(n1,n1) +  dble(psen(n1))
c$$$         bb(n1,n1) = 1d0
c$$$      end do
c$$$      
c$$$
c$$$      
c$$$      matz=2
c$$$      call  rsg(nps,nps,H,bb,w,matz,CI,fv1,fv2,ierr)
c$$$      if(ierr .ne. 0) then
c$$$         print*, ' iterimp_rel_la:  ierr =', ierr
c$$$      end if
c$$$
c$$$      psen(1:nps) = w(1:nps)
c$$$
c$$$      i_max=0
c$$$      do n=1,nps
c$$$         i_max = max(i_max,maxps(n))
c$$$      end do
c$$$      do i=1,i_max              ! meshr
c$$$         do n=1,nps
c$$$            tmp = 0d0
c$$$            do m=1,nps
c$$$               if(minps(m).le.i.and.i.le.maxps(m)) then
c$$$                  tmp = tmp + CI(m,n)*dble(ps(i,m))
c$$$               end if
c$$$            end do
c$$$            ps_tmp(n) = tmp
c$$$         end do
c$$$         ps(i,1:nps) = ps_tmp(1:nps)
c$$$      end do
c$$$      do n=1,nps
c$$$         call minmaxi(ps(1,n),meshr,i1,i2)
c$$$         maxps(n) = i2
c$$$         minps(n) = i1
c$$$      end do
c$$$      
c$$$c     Fix sign of the one-electron functions
c$$$      do n=1,nps         
c$$$         if( sign(1.0,ps(minps(n)+1,n)) .lt. 0)  then
c$$$            do i=minps(n),maxps(n) 
c$$$               ps(i,n) = -ps(i,n)
c$$$            end do
c$$$         end if
c$$$      end do
c$$$
c$$$c     update pseudostate arrays:      
c$$$      do m = nabot(la)-la, nps
c$$$         n = m + la
c$$$         enpsinb(n,la) = psen(m) 
c$$$         istoppsinb(n,la) = maxps(m)
c$$$         psinb(1:maxps(m),n,la) = ps(1:maxps(m),m)
c$$$      enddo      
c$$$
c$$$c     update new-pscore arrays:      
c$$$      do m = 1, nabot(la)-la-1
c$$$         enpscoren(m,la) = psen(m)
c$$$         istopcoren(m,la) = maxps(m) 
c$$$         pscoren(1:maxps(m),m,la) = ps(1:maxps(m),m) 
c$$$      enddo
c$$$c
c$$$      if(iter .eq. 1) then
c$$$         do k=1,min(10,nps)
c$$$            print'(2i3,i6,12x,f16.8,"  eV, ",f16.8,"  au")',
c$$$     >         la,k+la,maxps(k),w(k)*13.6058, w(k)/2.0
c$$$         end do
c$$$      end if
c$$$      
c$$$      return 
c$$$      end  subroutine iterimp_rel_la
c$$$
c$$$
c$$$c--------------------------------------------------------------------------
c$$$      subroutine makenewvdcore(vdcore,minvdc,maxvdc,nznuci,u)
c$$$      include 'par.f'
c$$$      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
c$$$     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
c$$$      common /relcore/ pscore(2000,5,0:3), enpscore(5,0:3),
c$$$     >   istopcore(5,0:3)
c$$$      common/meshrr/ meshr,rmesh(maxr,3)
c$$$      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
c$$$     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
c$$$      dimension vdcore(maxr), vnl(maxr), fun(maxr), u(maxr)
c$$$      
c$$$      minvdc = meshr
c$$$      maxvdc = 0
c$$$c      print '('' n l of the core states'')'
c$$$      do lac = 0, lactop
c$$$         do nac = 1, nabot(lac) - lac  - 1
c$$$c            print '(2i2)', nac, lac
c$$$            minfun = 1
c$$$            maxfun = istopcore(nac,lac)
c$$$            do i = minfun, maxfun
c$$$               fun(i) = pscore(i,nac,lac)*pscore(i,nac,lac)*rmesh(i,3)
c$$$            end do 
c$$$            call form(fun,minfun,maxfun,rpow1(1,0),rpow2(1,0),
c$$$     >         minrp(0),maxrp(0),meshr,vnl,i1,i2)
c$$$            call nuclear(fun,.true.,minfun,maxfun,i2,u,nznuc,vnl)
c$$$            if (i1.lt.minvdc) minvdc = i1
c$$$            if (i2.gt.maxvdc) maxvdc = i2
c$$$            do i = minvdc, maxvdc
c$$$               vdcore(i) = vdcore(i) + float(4 * lac + 2) * vnl(i)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      return
c$$$      end   subroutine makenewvdcore
c$$$c--------------------------------------------------------------------------
c$$$c     Note: radial integration weights (gridr(i,3)) are included in functions f1(i) and f2(i)
c$$$      subroutine newtmpfcexch(f1,min1,max1,nr,f2,res,l)
c$$$      include 'par.f'
c$$$      real f1(nr),f2(nr),res
c$$$      real ps2,psen2,enpsinb,psinb,rpow1,rpow2,cntfug,fun(maxr),
c$$$     >   vnl(maxr)
c$$$      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
c$$$     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
c$$$      common /worksp/
c$$$     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
c$$$      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
c$$$     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
c$$$      common /relcore/ pscore(2000,5,0:3), enpscore(5,0:3),
c$$$     >   istopcore(5,0:3)
c$$$
c$$$      res = 0d0
c$$$      if (ntype.eq.-3) return
c$$$      do lac = 0, lactop
c$$$         do nac = 1, nabot(lac) - lac - 1
c$$$C The above form fails for sodium if latop = 0.
c$$$C In the form below we assume that there are no d states in the core
c$$$            minfun = 1
c$$$            maxfun = istopcore(nac,lac)
c$$$            do i = minfun, maxfun
c$$$               fun(i) = pscore(i,nac,lac) * f2(i)
c$$$            end do 
c$$$            do 20 ilt = -lac, lac, 2
c$$$               lt = l + ilt
c$$$               if (lt.lt.0.or.lt.gt.ltmax) go to 20
c$$$               call cleb(2*lac,2*lt,2*l,0,0,0,c)
c$$$               const = - float(2 * lac + 1) * c * c /
c$$$     >            float(2 * l + 1)
c$$$               call form(fun,minfun,maxfun,rpow1(1,lt),rpow2(1,lt),
c$$$     >            minrp(lt),maxrp(lt),maxfun,vnl,i1,i2)
c$$$               tmp = 0d0
c$$$               do i = max(min1,i1), min(i2,maxfun,max1)
c$$$                  tmp = tmp + pscore(i,nac,lac) * vnl(i) * f1(i)
c$$$               end do
c$$$               res = res + tmp * const
c$$$ 20         continue 
c$$$         enddo 
c$$$      enddo 
c$$$      return
c$$$      end  subroutine newtmpfcexch
c$$$
c$$$
c$$$c--------------------------------------------------------------------------
c$$$c     res - is in atomic units
c$$$      subroutine Ham1el(nr,la,vdcore,il,im,f2,f1,gl1,grid,weight,res)
c$$$      real f1(nr),f2(nr),gl1(nr),vdcore(nr),grid(nr),weight(nr)
c$$$      real tmp, res, res1, tmpd
c$$$      real g(nr), ff1(nr), ff2(nr)
c$$$      double precision  fund(nr)
c$$$
c$$$      res1 = 0.0
c$$$      tmp = 0.0
c$$$      tmpd = 0.0
c$$$      tmpl = 0.0      
c$$$      
c$$$      ff1(:) = f1(:)*weight(:)
c$$$      ff2(:) = f2(:)*weight(:)
c$$$      
c$$$      g(il:im) = f2(il:im)*(vdcore(il:im)-2.0/grid(il:im))*f1(il:im)*
c$$$     >   weight(il:im)
c$$$      tmp = SUM(g(il:im))
c$$$c      call tmpfcexch(ff2,il,im,nr,ff1,res1,la)
c$$$      call secondderiv(f1, il, im, grid, nr, fund)
c$$$      
c$$$      if(l .ne. 0) then
c$$$         g(il:im) =  dble(la*(la + 1))*f2(il:im)*f1(il:im)*
c$$$     >      weight(il:im)/grid(il:im)/grid(il:im)
c$$$         tmpl = SUM(g(il:im))
c$$$      end if
c$$$c      tmpd = SUM(f2(il:im)*fund(il:im)*weight(il:im))
c$$$      tmpd = SUM(f2(il:im)*gl1(il:im)*weight(il:im))
c$$$      res = -(tmpd-tmpl)/2.0 + tmp + res1
c$$$      
c$$$      return
c$$$      end subroutine Ham1el
c$$$c--------------------------------------------------------------------------
c$$$      subroutine checkdiag(vdcore,l,al,nps,iih)
c$$$      include 'par.f'
c$$$      real vdcore(maxr)
c$$$      common/meshrr/ nr,gridr(maxr,3)
c$$$      double precision grid(nr), x(nr), pl(nps), alsp, DM
c$$$      real  fl(nr,nps), gl(nr,nps)
c$$$      integer minf(nps), maxf(nps), ming(nps), maxg(nps)
c$$$      double precision  f8(nr,nps), f8p(nr,nps)
c$$$      real en(nps)
c$$$      common /factratio/ Cgr(0:lomax,komax)
c$$$      double precision  Cgr
c$$$      double precision H(nps,nps),bb(nps,nps),w(nps),
c$$$     >   CI(nps,nps), tmp
c$$$      double precision  fv1(nps),fv2(nps)
c$$$      real  ps_tmp(nps)
c$$$      double precision  Rel_cor, Rel_cor_new
c$$$      double precision  alb
c$$$      
c$$$      nznuc = 80
c$$$
c$$$      vdcore(:) = -78.0/gridr(:,1)
c$$$      
c$$$      fl = 0.0
c$$$      gl = 0.0
c$$$      f8 = 0.0
c$$$      pl = 0.0
c$$$      H = 0.0
c$$$      bb = 0.0
c$$$      fv1 = 0.0
c$$$      fv2 = 0.0
c$$$      w = 0.0
c$$$      CI = 0.0
c$$$      
c$$$      al = al*2.0
c$$$      alsp = dble(al)
c$$$      
c$$$      grid(1:nr) = dble(gridr(1:nr,1))
c$$$
c$$$      call factorials
c$$$      
c$$$      il = 2*l+2
c$$$      call  lagpol8(il,alsp,f8,nr,nps,grid,nr,pl)
c$$$      call  lagpol8(il+1,alsp,f8p,nr,nps,grid,nr,pl)
c$$$      
c$$$      x(:) = grid(:)*alsp
c$$$
c$$$
c$$$      do k=1,nps
c$$$          DM = dsqrt(alsp/Cgr(l,k))
c$$$         fl(1:nr,k) = DM*(x(1:nr)**(l+1))*dexp(-x(1:nr)/2.0D0)*
c$$$     >       f8(1:nr,k)
c$$$         gl(1:nr,k) = (dble(l*(l+1))/(x(1:nr)*x(1:nr)) -
c$$$     >      dble(k+l)/x(1:nr) +
c$$$     >      0.25d0)*(alsp*alsp) * fl(1:nr,k)
c$$$
c$$$         call minmaxi(fl(1,k),nr,i1,i2)
c$$$         minf(k)=i1
c$$$         maxf(k)=i2
c$$$c
c$$$         if(k.gt.1) then
c$$$            gl(1:nr,k) = dble(gl(1:nr,k)) +
c$$$     >         (alsp*alsp)*DM*(x(1:nr)**l)*
c$$$     >         dexp(-x(1:nr)/2.0D0)*f8p(1:nr,k-1)
c$$$         end if
c$$$         call minmaxi(gl(1,k),nr,i1,i2)
c$$$         ming(k)=i1
c$$$         maxg(k)=i2
c$$$      end do
c$$$
c$$$c$$$      alb = 60.0
c$$$c$$$      fl(1:nr,1) = dsqrt(alb)*2d0*alb*grid(1:nr)*dexp(-alb*grid(1:nr))
c$$$c$$$      gl(1:nr,1) = dsqrt(alb)*2d0*alb*alb*(-2d0+alb*grid(1:nr))*
c$$$c$$$     >   dexp(-alb*grid(1:nr))
c$$$c$$$
c$$$c$$$      call minmaxi(fl(1,1),nr,i1,i2)
c$$$c$$$      minf(1)=i1
c$$$c$$$      maxf(1)=i2
c$$$c$$$      call minmaxi(gl(1,1),nr,i1,i2)
c$$$c$$$      ming(1)=i1
c$$$c$$$      maxg(1)=i2
c$$$c$$$     
c$$$c$$$      do n=2,nps
c$$$c$$$         il = max(minf(1),minf(n))
c$$$c$$$         im = min(maxf(1),maxf(n))
c$$$c$$$         bb(1,n) =  SUM(fl(il:im,1)*fl(il:im,n)*
c$$$c$$$     >      gridr(il:im,3))
c$$$c$$$         bb(n,1) = bb(1,n)
c$$$c$$$      end do
c$$$
c$$$      do n1=1,nps     
c$$$         do n2=1,n1            
c$$$            il = max(minf(n1),minf(n2))
c$$$            im = min(maxf(n1),maxf(n2))
c$$$            call  Ham1el(nr,l,vdcore(1:nr),il,im,fl(1:nr,n2),
c$$$     >         fl(1:nr,n1),gl(1:nr,n1),
c$$$     >         gridr(1:nr,1),gridr(1:nr,3),res)
c$$$            if(iih .ne. 0) then
c$$$c     resrel =  Rel_cor(nr,fl(1:nr,n2),fl(1:nr,n1),
c$$$c     >            l,minf(n2),maxf(n2),minf(n1),maxf(n1),
c$$$c     >            gridr(1:nr,1),gridr(1:nr,3),nznuc)
c$$$               resrelnew =  Rel_cor_new(nr,fl(1:nr,n2),fl(1:nr,n1),
c$$$     >            gl(1:nr,n2),gl(1:nr,n1),
c$$$     >            l,minf(n2),maxf(n2),minf(n1),maxf(n1),
c$$$     >            gridr(1:nr,1),gridr(1:nr,3),nznuc)
c$$$            else
c$$$               resrel = 0.0
c$$$               resrelnew = 0.0
c$$$            end if
c$$$            res = res +  resrelnew ! resrel
c$$$            H(n1,n2) = 2.0*res
c$$$            H(n2,n1) = H(n1,n2)
c$$$         end do
c$$$         
c$$$         bb(n1,n1) = 1d0
c$$$         
c$$$      end do
c$$$      
c$$$      matz=2
c$$$      call  rsg(nps,nps,H,bb,w,matz,CI,fv1,fv2,ierr)
c$$$      if(ierr .ne. 0) then
c$$$         print*, '********* ierr =', ierr
c$$$         stop
c$$$      end if
c$$$      
c$$$      en(1:nps) = w(1:nps)
c$$$
c$$$      i_max=0
c$$$      do n=1,nps
c$$$         i_max = max(i_max,maxf(n))
c$$$      end do
c$$$      do i=1,i_max
c$$$         do n=1,nps
c$$$            tmp = 0d0
c$$$            do m=1,nps
c$$$               if(minf(m).le.i.and.i.le.maxf(m)) then
c$$$                  tmp = tmp + CI(m,n)*dble(fl(i,m))
c$$$               end if
c$$$            end do
c$$$            ps_tmp(n) = tmp
c$$$         end do
c$$$         fl(i,1:nps) = ps_tmp(1:nps)
c$$$      end do
c$$$      do n=1,nps
c$$$         call minmaxi(fl(1,n),nr,i1,i2)
c$$$         maxf(n) = i2
c$$$         minf(n) = i1
c$$$      end do
c$$$      
c$$$c     Fix sign of the one-electron functions
c$$$      do n=1,nps         
c$$$         if( sign(1.0,fl(minf(n)+1,n)) .lt. 0)  then
c$$$            fl(minf(n):maxf(n),n) = -fl(minf(n):maxf(n),n)
c$$$         end if
c$$$      end do
c$$$
c$$$      
c$$$      print*,'New  wave-functions'
c$$$      do k=1,nps
c$$$         print'(2i3,i6,12x,f16.8,"  eV, ",f16.8,"  au")',
c$$$     >      l,k,maxf(k),en(k)*13.6058, en(k)/2.0
c$$$      end do
c$$$
c$$$      
c$$$      return
c$$$      end subroutine checkdiag
c$$$c--------------------------------------------------------------------------
c$$$c     Mass term: (p^2)^2, Darwin term Z delta(r), atomic units
c$$$      double precision function  Rel_cor_new(nr,psf,psi,gf,gi,lsp,
c$$$     >   minf,maxf,mini,maxi,gridr,weightr,nznuc)
c$$$      real  psf(nr),  psi(nr), gf(nr),  gi(nr)
c$$$      double precision, dimension(nr) :: fun_f, fun_i, fun
c$$$      real, dimension(nr) :: gridr, weightr
c$$$      double precision  tmp
c$$$      double precision  alpha
c$$$      data alpha / 7.29714d-3 /    ! 1d0/137.04d0
c$$$      
c$$$c     Mass term
c$$$      fun = 0d0; fun_f = 0d0; fun_i = 0d0
c$$$      
c$$$      minfun=max(mini,minf)
c$$$      maxfun=min(maxi,maxf)
c$$$      fun_f = gf
c$$$      fun_i = gi
c$$$      
c$$$c      call secondderiv(psf(:), minfun, maxfun, gridr(:), nr, fun_f)
c$$$c      call secondderiv(psi(:), minfun, maxfun, gridr(:), nr, fun_i)
c$$$      if(lsp .ne. 0) then
c$$$         fun_f(minfun:maxfun) = fun_f(minfun:maxfun) -
c$$$     >      dble(lsp*(lsp + 1)) *dble(psf(minfun:maxfun)/
c$$$     >      gridr(minfun:maxfun)/gridr(minfun:maxfun))
c$$$         fun_i(minfun:maxfun) = fun_i(minfun:maxfun) -
c$$$     >      dble(lsp*(lsp + 1)) *dble(psi(minfun:maxfun)/
c$$$     >      gridr(minfun:maxfun)/gridr(minfun:maxfun))
c$$$      end if
c$$$      fun(minfun:maxfun) = fun_f(minfun:maxfun)*fun_i(minfun:maxfun)*
c$$$     >   dble(weightr(minfun:maxfun)) 
c$$$
c$$$      tmp = ((psf(2)-2.0*psf(1))/gridr(1)/gridr(1)) *
c$$$     >   ((psi(2)-2.0*psi(1))/gridr(1)/gridr(1)) *
c$$$     >   weightr(1)/4.0
c$$$      Rel_cor_new = -(SUM(fun) + tmp)
c$$$c
c$$$      
c$$$
c$$$c     Darwin term
c$$$      tmp = 0.0
c$$$      if(lsp .eq. 0) then
c$$$         tmp = dble(nznuc)*((4d0*psf(1)-psf(2))/2d0/gridr(1))*
c$$$     >      ((4d0*psi(1)-psi(2))/2d0/gridr(1))
c$$$      end if
c$$$c      print*,  lsp, real(Rel_cor_new*alpha*alpha/8d0),
c$$$c     >   real(tmp*alpha*alpha/8d0),
c$$$c     >   real((Rel_cor_new+tmp)*alpha*alpha/8d0)
c$$$     
c$$$      Rel_cor_new =  (Rel_cor_new  + tmp) * alpha*alpha/8d0
c$$$      
c$$$      return
c$$$      end  function Rel_cor_new

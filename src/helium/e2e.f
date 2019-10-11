      subroutine e2eovp(n,nspm,lo,fl,maxf,minf,pnewC,Nmax,namax,
     >   E,en,lchi,lschi,chi, minc, maxc, ovlp)
      use CI_MODULE
      include 'par.f'
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      dimension maxf(nspmax),minf(nspmax),fl(nmaxr,nspmax),lo(nspmax)
      double precision  E(KNM),ortint,r1elk,orte2en2
      dimension  chi(nmaxr)
      data pi/3.141592654/
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer na, nam
      integer pnewC
      double precision tmp1, tmp
ccc      double precision  C
ccc      pointer (pC, C(Nmax,namax,namax))
      
ccc      pC = pnewC
!      n = n -1
c     calculate overlap : <(fl(1s) chi(k,l)):l s | n(n1n2),l_n s_n>
c     coulomb wave: (fl(1s) chi(k,l)) --> l = l_n, l = lo(n1) or l = lo(n2)
c     spin of the coulomb wave must be equal to the spin of the state n: s=s_n
c     n - is the number of the helium state we need overlap with.

      tmp = 0.0
      if(ll(n).eq.lchi.and.ls(n).eq.lschi) then
         do jn1=1,nam(n)
            n1 = na(jn1,n)
            if (ortint(1,n1) .eq. 0) cycle
            do jn2=1,nam(n)
               n2 = na(jn2,n)
               if(C(n,jn1,jn2).ne.0d0) then
                  orte2en2 = 0d0
                  if(lo(n2).eq.lchi) then
                     orte2en2 = r1elk(0,chi,fl(1,n2),
     >                  minc,minf(n2),maxc,maxf(n2))
                  end if
!                  print*, '--> old: orte2en2=', orte2en2, ortint(1,n1)
                  tmp1 = C(n,jn1,jn2)*ortint(1,n1)*orte2en2
                  tmp = tmp + tmp1
               end if
            end do
         end do
         tmp = tmp *sqrt(2.0)   !  sqrt(2.0) comes from an application of the antisym. operator A to the target state Phi_n
      end if
      ovlp = tmp * sqrt(2.0 / pi)
c$$$      write(*,'("l,s=",2I2,", en=",F7.4, " atom(n=",I3,", E_n=",
c$$$     >   F7.4,")",", la,sa=",2I3,", overlap=",F7.4)')
c$$$     >   lchi,lschi,en,n,E(n),ll(n),ls(n),overlap(n)
      return
      end

C--------------------------------------------------------------------------------
      subroutine get_e2eovp_gen(filename,n,nspm,lo,fl,maxf,minf,namax,
     >   en,lw,lsw,lparw,ovlp)
      
      include 'par.f'
      complex ovlp
      
      dimension maxf(nspmax),minf(nspmax),fl(nmaxr,nspmax),lo(nspmax)
      character  filename*16

      call limits_cw(en,lw,lsw,lparw,nchm,nrm,filename) ! get nchm, nrm, filename
!      print*, '-----> ', nchm, nrm, filename
      call e2eovp_gen(n,nspm,lo,fl,maxf,minf,namax,
     >     en,lw,lsw,lparw,nchm,nrm,filename,ovlp)

      return
      end
C--------------------------------------------------------------------------------
      subroutine limits_cw(ew,lw,lsw,lparw,nchm,nrm,filename)
            
      character  filename*16, tmpchar*1

c$$$      if( lparw .eq. (-1)**lw ) then
c$$$         write(filename,'(1p,E10.4,"_",I1,I1,"0")') ew,lw,lsw
c$$$      elseif ( lparw .eq. -(-1)**lw ) then
c$$$         write(filename,'(1p,E10.4,"_",I1,I1,"1")') ew,lw,lsw
c$$$      else
c$$$         print*,'wrong value of parity for the continuum wave'
c$$$         print*,'lparw = ', lparw
c$$$         stop
c$$$      endif
c$$$      
c$$$      filename = 'tmpfile'
      open(115,file=filename, status='OLD', err=10)
      
      read(115,*) nrm , nchm

      return

 10   print*, 'Errror openning file:', filename
      stop

 20   print*, 'Errror reading file:', filename
      stop

      end
C--------------------------------------------------------------------------------
      subroutine e2eovp_gen(n,nspm,lo,fl,maxf,minf,namax,
     >   en,lw,lsw,lparw,nchm,nrm,filename,ovlp)
      use CI_MODULE
      include 'par.f'
      common /meshrr/ nr, gridr(nmaxr,3)
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      dimension maxf(nspmax),minf(nspmax),fl(nmaxr,nspmax),lo(nspmax)
      double precision  E(KNM),ortint,r1elk, ort1
      integer lw, lsw, lparw
      complex  chi(nrm,0:nchm)
      dimension lchi(0:nchm),nchcore(0:nchm),minc(0:nchm),maxc(0:nchm)
      character  filename*16
      complex ovlp
      data pi/3.141592654/
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer na, nam
      complex orte2en2

c     calculate overlap : <(fl(i) chi(k,l,i)):lw sw | n(n1n2),l_n s_n>
c     CC-wave: 1/sqrt{2}[fl(i) chi(k,l,i) + (-1)^s*(-1)^{lw-l_i-l} chi(k,l,i) fl(i)]
c     spin of the coulomb wave must be equal to the spin of the state n: s=s_n
c     n - is the number of the helium state we need overlap with.
      
      ovlp = (0.0,0.0)

      if(ll(n).ne.lw) return
      if(ls(n).ne.lsw) return
      if(lparity(n).ne.lparw) return

      call read_cw(filename,nchm,nrm,lchi,nchcore,chi,minc,maxc)

      do i=0,nchm
         nc = nchcore(i)
         lnc = lo(nc)

         do jn1=1,nam(n)
            n1 = na(jn1,n)
            ort1 = ortint(nc,n1)
            if(ort1 .eq. 0) cycle
            do jn2=1,nam(n)
               n2 = na(jn2,n)
               Ctmp = C(n,jn1,jn2)
               if(Ctmp.eq.0d0) cycle

               if(lo(n2).ne.lchi(i)) cycle

               k1 = max(minc(i),minf(n2))
               k2 = min(maxc(i),maxf(n2))               
               if(k1 .gt. k2) return
!               orte2en2 = ovlp_comp(chi(1,i),minc(i),maxc(i),
!     >              fl(1,n2),minf(n2),maxf(n2))
               orte2en2 = SUM(chi(k1:k2,i)*fl(k1:k2,n2)*gridr(k1:k2,3))
!               print*, 'orte2en2 =' ,orte2en2
               ovlp = ovlp + Ctmp* ort1 * orte2en2 

            end do
         end do
      enddo                     ! end i loop
      ovlp = ovlp*dsqrt(2d0)
      ovlp = ovlp * sqrt(2.0 / pi)   

!      print*, '***> new ovlp =', ovlp
      return
      end
c
c--------------------------------------------------------------------------------------------------------
c     
      subroutine read_cw(filename,nchm,nrm,lchi,nchcore,chi,minc,maxc)

      include 'par.f'
      common /corearray/ nicm, ncore(nspmCI)  
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /meshrr/ nr, gridr(nmaxr,3)

      integer nchm, nrm
      complex  chi(nrm,0:nchm)
      dimension lchi(0:nchm),nchcore(0:nchm),minc(0:nchm),maxc(0:nchm)      
      character  filename*16      
      integer lst(0:nchm), nst(0:nchm)
      real gridw(nmaxr), ft, fr(nchm), fi(nchm)

      rewind(115)
      
      read(115,*) itmp, jtmp, nst(0),lst(0),lchi(0),
     >     ( nst(i),lst(i),lchi(i), i=1,nchm )

      chi(:,:) = (0.0,0.0)

      do j=1,nrm
         read(115,*) gridw(j), ft, ( fr(i),fi(i), i=1,nchm )
         chi(j,0) = cmplx(ft,0.0)
         do i=1,nchm
            chi(j,i) = cmplx(fr(i),fi(i))
         enddo
      enddo

      minc = 1
      maxc = nrm


c     Check that the r-grid for the continuum function read above from the disk is the same as the current r-grid
      if(nrm .ne. nr) then
         print*, 'e2e.f: read_cw(): the current r-grids is not the ',
     >        'same for the r-grid for the read cont.func.'
         print*, 'nr, nrm = ', nr, nrm
         stop 'e2e.f: 1. radial grids are not the same'
      endif
      if(nint(abs(gridw(1) - gridr(1,1))/gridw(1)) .gt. 1e-4 )then
         print*, 'e2e.f: read_cw(): the current r-grid is not the ',
     >        'same for the r-grid for the read cont.func.'
         print*, 'gridr(1,1), gridw(1) = ', gridr(1,1), gridw(1) 
         stop 'e2e.f: 2. radial grids are not the same'
      endif
      if(nint(abs(gridw(nr) - gridr(nr,1))/gridw(nr)) .gt. 1e-4 )then
         print*, 'e2e.f: read_cw(): the current r-grid is not the ',
     >        'same for the r-grid for the read cont.func.'
         print*, 'gridr(nr,1), gridw(nr) = ', gridr(nr,1), gridw(nr) 
         stop 'e2e.f: 3. radial grids are not the same'
      endif

c

      nchcore(:) = 0

      do i=0,nchm
         do nc=1,nicm
            nsp = ncore(nc)
            lnsp = lo(nsp)
            knsp = ko(nsp)
            if(lst(i) .eq. lnsp .and. nst(i) .eq. lnsp+knsp) then
               nchcore(i) = nsp
!               print*, 'i,nsp,lnsp,knsp:',  i,nsp,lnsp,knsp
               exit
            endif
         enddo
      enddo

      close(115)

      return 
      
      end
c
c---------------------------------------------------------------------------------------------------------------
c     
      subroutine write_he_st(N)

      use CI_MODULE

      include 'par.f'

      common /meshrr/ nr, gridr(nmaxr,3)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common /corearray/ nicm, ncore(nspmCI)      
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      
      integer  n1_ch(nspm),n2_ch(nspm),l1c(nspm),n1c(nspm),lfc(nspm)
      real Car(nspm)

      open(215,file='tmpfile')
      
      nrm = 0
      do j=1,nam(N)
         nsp = na(j,N)
         nrm = max(nrm,maxf(nsp))         
      enddo

      nch = 0

      do ic=1,nicm
         nic = ncore(ic)
         do j1=1,nam(N)
            n1 = na(j1,N)
            if(n1 .ne. nic) cycle 
            do j2=1,nam(N)
               n2 = na(j2,N)
               Ctmp = C(N,j1,j2)
               if(Ctmp .eq. 0.0) cycle
               nch = nch + 1
               Car(nch) = Ctmp
               n1_ch(nch) = nic
               n2_ch(nch) = n2
               l1c(nch) = lo(nic)
               n1c(nch) = ko(nic) + lo(nic)
               lfc(nch) = lo(n2)
            enddo
         enddo
      enddo

      nchm = nch
      
      write(215,'(I8,I5,4X,3i5,100(8X,3i5,7X))') nrm,nchm,
     >     n1c(1),l1c(1),lfc(1),
     >     (n1c(nch),l1c(nch),lfc(nch), nch=1,nchm)

      zero = 0.0
      do i=1,nrm
         write(215,'(E15.7,2X,E15.7,2X,100E15.7)') gridr(i,1), zero,
     >        (Car(nch)*fl(i,n2_ch(nch)),zero, nch=1,nchm)
      enddo

      close(215)

      return
      end
c---------------------------------------------------------------------------------------------------------------
c     
      subroutine write_fc_cw(nrm,lchi,chi,minc,maxc,phase,sigc)

      include 'par.f'

      common /meshrr/ nr, gridr(nmaxr,3)
      integer nchm, nrm, lchi, minc,maxc
      real  chi(nrm)
      data pi/3.141592654/
      complex phase, sigc

      open(315,file='tmpfile')

      lzero = 0
      ncst = 1

      zero = 0.0

      write(315,'(2i5,100(8X,3i5,7X))') nrm, 1, 
     >     ncst,lzero,lchi, ncst,lzero,lchi

      do j=1,nrm
         tmp1 = real(phase*chi(j))*0.3
         tmp2 = real(phase*chi(j))*0.7
         tmp3 = aimag(phase*chi(j))
         write(315,'(1P,E15.7,2X,100E15.7)') 
     >        gridr(j,1), tmp1, tmp2, tmp3
      enddo

      close(315)

      return
      end

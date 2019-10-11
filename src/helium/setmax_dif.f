      subroutine set_dif_orb(fl,maxf,minf)
      include 'par.f'
      common /meshrr/ nr, gridr(nmaxr,3)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)      

      integer  nkd(0:lomax)
      double precision  aldif(0:lomax)
      integer ndif
      integer kom
      integer, dimension(:), allocatable:: lodif,kodif
      
c     
      call  factorials

c     read information on one electron orbitals
      write(4,'(" description of  dif  s.p. orbitals")') 
      
      call inputdatamax_dif(nkd,aldif,ndif)
      if(ndif .eq. 0) then
         print*, "ndif = 0, No substittution of the orbitals",
     >        " will be done"
         write(4,'("ndif = 0, No substittution of the orbitals",
     >       " will be done")')
         return
      endif

c     find maximum value of  l  in the basis 
      lnspmax = 0
      do n=1,nspm
         lnspmax = max(lnspmax, lo(n))
      enddo

      nsp=0
      do l=0,min(lnspmax,lomax)
         n = nkd(l)
         if(n.ne.0) then
            do k=1,n
               nsp=nsp+1
            end do
         end if
      end do

      nspmdif = nsp
      allocate(lodif(nspmdif),kodif(nspmdif))

      nsp=0
      do l=0,lomax
         n = nkd(l)
         if(n.ne.0) then
            do k=1,n
               nsp=nsp+1
               lodif(nsp)=l
               kodif(nsp)=k
            end do
         end if
      end do
      nspmdif = nsp
      print*,'nspmdif=',nspmdif

      kom = 0
      do nsp=1,nspmdif
         kom = max(kom,kodif(nsp))
      end do

      call  set_sp_dif(fl,maxf,minf, 
     >     nkd,nspmdif,lodif,kodif,aldif,kom)

      deallocate(lodif,kodif)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  set_sp_dif(fl,maxf,minf, 
     >     nkd,nspmdif,lodif,kodif,aldif,kom)
      include 'par.f'

      common /meshrr/ nr,gridr(nmaxr,3)
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      integer  nkd(0:lomax)
      integer nspmdif
      integer lodif(nspmdif), kodif(nspmdif)
      double precision  aldif(0:lomax)
      integer nr, kom

      double precision  Cgr
      common /factratio/ Cgr(0:lomax,komax)
      double precision  ortint, e1r
      common /ortog/  ortint(nspmax,nspmax)
      common/hame1/ e1r(nspmax,nspmax)
      integer nspm, lo, ko, nset
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)

      double precision M(nspmdif),  als(nspmdif)
      real fldif(nr,nspmdif)
      integer maxfdif(nspmdif),minfdif(nspmdif)
      integer  maxg(nspmdif),ming(nspmdif)
      real  gl(nr,nspmdif)
      double precision  f8(nr,kom)
      double precision  rnorm, r1elk, re1, oneelint_A

      do nsp=1,nspmdif
         l = lodif(nsp)
         k = kodif(nsp)
         als(nsp) = aldif(l)
         M(nsp) = dsqrt(als(nsp)/Cgr(l,k))
c         print*, 'nsp,als(nsp)=',nsp,als(nsp)
      end do

  

      call funlaguer_dif(nspmdif,lodif,kodif,
     >     M,als,fldif,maxfdif,minfdif,gl,maxg,ming,f8,kom)



      do nsp=1,nspmdif
         print*, 'dif orb.: nsp=', nsp
         l = lodif(nsp)
         k = kodif(nsp)
         
C     find original orbital with lo = l, ko = k
         nsporig = 0
         do nn=1,nspm
            if(lo(nn) .eq. l .and. ko(nn).eq.k) then
               nsporig = nn
               exit
            endif
         enddo
         if(nsporig .eq. 0) then
            print*, "setmax_dif.f: set_sp_dif(): nsporig=0"
            stop
         endif
         print*, 'found nsporig=', nsporig, 
     >        '  will overwrite it with dif. orbital'

C     overwrite original orbital with dif orbital (only radial part)               
         i2 = maxfdif(nsp)
         i1 = minfdif(nsp)
         maxf(nsporig) = i2
         minf(nsporig) = i1
         fl(1:i1,nsporig) = 0.0
         fl(i2:nr,nsporig) = 0.0
         fl(i1:i2,nsporig) = fldif(i1:i2,nsp)

         print*,'nsporig, fl:', nsporig, fl(100,nsporig)

C     calculate overlap matrix for this overwritten orbital (nsporig)
         do nn=1,nspm
            if(lo(nn) .ne. lo(nsporig)) cycle
            if(nn .eq. nsporig) then
               ortint(nn,nn) = 1d0
            else
               n1=nn
               n2=nsporig
               i1=max(minf(n1),minf(n2))
               i2=min(maxf(n1),maxf(n2))
               rnorm =  SUM(fl(i1:i2,n1)*fl(i1:i2,n2)*gridr(i1:i2,3))
               ortint(n1,n2) = rnorm
               ortint(n2,n1) = rnorm
            endif                  
         enddo

C     calculate e1r matrix for this overwritten orbital (nsporig)
         do nn=1,nspm
            
            if(lo(nn) .ne. lo(nsporig)) cycle
            n1=nn
            n2=nsporig
            
            re1 = oneelint_A(n1,n2,fl,maxf,minf,
     >           gl(1,nsp),maxg(nsp),ming(nsp))
            e1r(n1,n2) = re1
            e1r(n2,n1) = re1
            
         enddo

      enddo ! end nsp loop

c      call gsort(nspm,fl,maxf,minf,ortint,e1r,lo)

c      print*,'ortint'
c      do n=1,nspm
c         write(*,'(100E12.4)') (ortint(n,np), np=1,nspm)
c      end do
c      print*,'e1r'
c      do n=1,nspm
c         write(*,'(100E12.4)') (e1r(n,np), np=1,nspm)
c      end do

      
      return
      end
c-----------------------------------------------------------------
c****************  Read overall input data for dif orbitals  **********************
c-----------------------------------------------------------------
      subroutine inputdatamax_dif(nkd,aldif,ndif)
c     this is to read overall maximum on s.p. functions - for all HE states: only for diff orbitals
      include 'par.f'
      integer  nkd(0:lomax)
      double precision aldif(0:lomax)
      integer ndif

      dimension raldif(0:lomax)
      double precision  Z
      common /Zatom/ Z
      real rZ
      character string*40


      read(3,*) rZ
      if (abs(rz-Z).gt.1e-4) then
         print*, 'Z(from F5) =', rz, ', Z(from ccc.in) = ', real(Z)
         stop 'Z in F5 not same as in ccc.in'
      end if
c--
      do l=0,lomax
         nkd(l) = 0
      end do
c--
      write(4,'(" description of s.p. states")') 
c     nset=1

      read(3,'(a40)') string
      print '(''Skipping line: '',a40)', string
      write(4,'("skip line l1max =")')

      read(3,'(a40)') string
      print '(''Skipping line: '',a40)', string
      write(4,'("skip line al(l) =")')

      read(3,'(a40)') string
      print '(''Skipping line: '',a40)', string
      write(4,'("skip line nk(l) =")')

      raldif(:) = 0.0
c     nset=2
      read(3,*) ndif
      write(4,'("ndif =",I5)') ndif
      read(3,*) ldif
      write(4,'("ldif =",I5)') ldif
      read(3,*) (raldif(l), l=0,ldif)
      write(4,'("aldif(l) =",10F10.5)') (raldif(l), l=0,ldif)
      read(3,*) (nkd(l), l=0,ldif)
      write(4,'("nkd(l) =",10I5)') (nkd(l), l=0,ldif)
      do l=0,lomax
         aldif(l) = DBLE(raldif(l))
      end do
      return
      end
c-----------------------------------------------------------------
c-----------------------------------------------------------
      subroutine funlaguer_dif(nspm,lo,ko,
     >     M,als,fl,maxf,minf,gl,maxg,ming,f8,kom)
c     this routine forms Sturmian functions (fl(i,nsp)) and second
c     derivatives (gl(i,nsp)) for each of them  from Laguerre polinomials.
      include 'par.f'
      common /meshrr/ nr,gridr(nmaxr,3)
      integer lo(nspm),ko(nspm)
      double precision  als(nspm), M(nspm)
      real  fl(nr,nspm)
      integer  maxf(nspm),minf(nspm)
      real gl(nr,nspm)
      integer  maxg(nspm),ming(nspm)
      double precision pl(komax), f8(nr,kom)            
      double precision grid(nr)
      double precision x

      grid(1:nr) = gridr(1:nr,1)
      
      do n=1,nspm
         l = lo(n)
         k = ko(n)
         il = 2*l+2
c         print*,'n,l,k=',n,l,k, 'M,als=',M(n), als(n)
         f8(:,:) = 0d0
         call  lagpol8(il,als(n),f8,nr,k,grid,nr,pl)
         do i=1,nr
            x = (gridr(i,1))*als(n)
            fl(i,n) = M(n)*(x**(l+1))*exp(-x/2.0D0)*f8(i,k)
c            gl(i,n) = 0.0
            gl(i,n) = ((l*(l+1))/(x*x) - (k+l)/x +
     >         0.25d0)*(als(n)**2) * fl(i,n)
         end do
         call minmaxi_nr(fl(1,n),nr,i1,i2)
         minf(n)=i1
         maxf(n)=i2
c
         if(k.gt.1) then
            call  lagpol8(il+1,als(n),f8,nr,k-1,grid,nr,pl)
            do i=1,nr
               x = DBLE(gridr(i,1))*als(n)
               gl(i,n) = dble(gl(i,n)) +
     >            (als(n)**2)*M(n)*(x**l)*dexp(-x/2.0D0)*f8(i,k-1)
            end do
         end if
         call minmaxi_nr(gl(1,n),nr,i1,i2)
         ming(n)=i1
         maxg(n)=i2
c         print*,'n=',n, 'i1,i2:',i1,i2, fl(100,n)
      end do

c      open(174,file='tmpfile')
c      do i=1,nr
c         write(174,'(i5,E15.4,1P,100E15.5)') i,gridr(i,1), 
c     >        (fl(i,n), n=1,nspm)
c      enddo
c      close(174)

      return
      end
c-----------------------------------------------------------
      subroutine minmaxi_nr(f,nr,i1,i2)
      include 'par.f'
      common/smallr/ formcut,regcut,expcut,fast
      logical fast
      dimension  f(nr)
      i=1
      do while (i.lt.nr.and.abs(f(i)).lt.regcut)
         i=i+1
      end do
      i1=i
      i=nr
      do while (i.gt.1.and.abs(f(i)).lt.expcut)
         i=i-1
      end do
      i2=i
      return
      end
c-----------------------------------------------------------
c-----------------------------------------------------------------
      double precision function oneelint_A(n1,n2,fl,maxf,minf,
     >     gl,maxg,ming)
      include 'par.f'
      common /meshrr/ nr, gridr(nmaxr,3)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      real  fl(nmaxr,nspmax)
      integer  maxf(nspmax),minf(nspmax)
      real  gl(nr)
      integer  maxg,ming
      double precision  Z, sum, r
      common /Zatom/ Z

      oneelint_A = 0.0D0
      if(lo(n1).ne.lo(n2)) return
      l = lo(n1)
      sum = 0d0

      maxrr=min(maxf(n1),maxg)
      minrr=max(minf(n1),ming)      
      do i=minrr,maxrr
         sum = sum + dble(gridr(i,3) * fl(i,n1) * gl(i))
      end do

      maxrr=max(maxf(n1),maxf(n2))
      minrr=min(minf(n1),minf(n2))
      do i=minrr,maxrr
         r = gridr(i,1)
         sum = sum + dble(gridr(i,3) * fl(i,n1) * 
     >      fl(i,n2))*(2d0*Z/r - dble(l*(l+1.))/(r*r))
      end do

      oneelint_A = - sum/2d0
      return
      end      
c************************************************************************


      subroutine corrATOM (rmax, nr, grid, gamma, r0)

C Helium ground state MCHartee-Fock wave function

C INPUT:   rmax   - maximal value of 'r'(the radial cut-off)
C          nr     - the number of points in r-grid
C          grid   - r-grid to calculate wavefunctions on it.

      
      include 'par.f'
      include 'paratom.f'
      PARAMETER (Ncf = 20)
      REAL*8 XN
      
      common/meshrr/ meshr,rmesh(maxr,3)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
     

      COMMON /cat1/   CI(0:Ncf, Ncf, Ncf)
      COMMON /cat2/   Nc, nmax(0:Ncf), jmax
      COMMON /cat3/   wa(maxr, 0:Ncf, Ncf)
      COMMON /cat4/   E(0:Ncf, Ncf)
      COMMON /cat5/   maxNR(0:Ncf, Ncf)
      COMMON /gstspin/   iSpin, meta
      dimension nmin(0:Ncf)

                        
      real  grid(nr)
*      real*8 over1,over2,over3
      real*8 z,r,alfa,bet,eps,h,ro1,!ro2,
     * r2(NSIZE),r5(NSIZE), r7(NSIZEk7),r8(NSIZE*NSIZEis),rm

      integer in(NSIZEis),kt,ne,ne5,ne5is,ksu

      CHARACTER*4 a4
      character *25  filnam
      
      character*1 hh(0:5)       !Spectroscopic labels
      data hh /'s','p','d','f','g','h'/

c DPI of Lithium       
      common /LITHIUM/ lithium,lithium2p
      
c Quantuum numbers

      data  in(1) /1/
      data  ksu   /1/

      real*8 a(15)

C  Rmax = h * ne
C  BET is the factor of the logarithmic part of the grid transformation
C  EPS precision factor
C  H is a step in the rho grid
C  NE is the number of points in the rho grid

       data  bet, eps, h, is, ne /1d0, 1.0d-5, 8.0d-2, 1, 1000/

C  KSU1 is 1 to calculate the complete set of SCFHF wave functions
C  KSU1 should be 3 for FCHF as only the single non-core state is needed

       data  ksu1, mu, dz, kz /1, 1, 0d0, 0/
c$$$       data  os, con, du1, du2 /'a.pri', 'list.lst', 'nul', 'con'/

c     Some protection

       if (nr .gt. nsizek7) then
          print*, 'Nr > Nsizek7 in corratom to at least', nr
          stop    'Nr > Nsizek7 in corratom'
       end if 
       if (ne .gt. nsize) then
          print*,'increase NSIZE in corratom to at least',ne
          stop   'increase NSIZE in corratom'
       end if 

* The grid

       kt  = nr
       do i = 1, kt
          r7(i) = grid(i)
       end do
       r   = rmax + 1.0
       
c  ---------------------------------------
c  Read information about the ground state
c  ---------------------------------------
       open(5,file='mchf',STATUS='UNKNOWN')
       read (5,*)   IZ, iSpin, Nc, Etot
       write(6,101) IZ, iSpin, Nc, Etot
 101   FORMAT (////12x,' MULTICONFIGURATION GROUND STATE   ',
     :           /12x,' -------------------------------   '/,
     :           /,   ' Nucleus charge           ',I3,
     :           /,   ' Spin                     ',I3,
     :           /,   ' Number of configurations ',I3,
     :           /,   ' Ground state energy, au  ',F10.6//)
       
       read(5,'(a25)') filnam

c     IZ - nuclear charge
c     Nc - number of two-electron configurations (nl)^2 in the GS
c     nmax(l) - maximum principle quantum number with orbital momentum l
c     CI(L, n) - contribution of the configuration (nl)^2 to the GS
       
c  Initialize CI variables
       
       do l = 0, Nc
          do n = 1, Nc
             do n1 = 1, Nc
                CI(l, n, n1) = 0
             end do
          end do
       end do
       
       jmax = 0
       S = 0
       do l = 0, Ncf
          nmax(l) = 0
          nmin(l) = 0
       end do
       
       meta = 0 !Ground state is assumed
       do i = 1, Nc
          read (5, 11) a4, x
          write(6, 11) a4, x
          call label(a4(1:2), n, l )
          call label(a4(3:4), n1,l1)
!          if(l.ne.l1) STOP 'NON-S MCHF is not implemented'
          if(n.ne.n1)  meta=1    !Metastable state
          if(l .gt. jmax)         jmax = l
          if(n .gt. nmax(l))      nmax(l) = n
          if(n1.gt. nmax(l))      nmax(l) = n1
          if(abs(x).gt.1e-10 .and. nmin(l).eq.0)    nmin(l) = n
          S = S + x**2
          CI(l,n,n1) = x
          if(i.eq.1) lin = l

C Lithium 2s/2p
          if(iz.eq.3 .and. lithium.eq.1) then
             if(i.eq.1 .and. n1.eq.2 .and. l1.eq.1) then
                lithium2p=1
                lithium  =0
                meta = 0
                print'(A)', 'Excited state Li 2p'
             end if
             if(i.eq.1 .and. n1.eq.2 .and. l1.eq.0) then
                lithium2p=0
                lithium  =1
                meta = 0
                print'(A)', 'Ground state Li 2s'
             end if
          end if
       end do


       
       write(6, 11) 'Total', S
 11    FORMAT(A4, E20.4)
                          
      OPEN(2,FILE=filnam,FORM='UNFORMATTED',STATUS='OLD')
      IZCH=2

      read(2) (a(jj1),jj1=1,15) 
      write(6,1954) (a(jj1),jj1=1,15)

      RM=a(2)
      NE=a(3)
      ne5 = ne + 5
      ne5is = ne5 * is
      H=a(4)
      bet=a(5)
      eps=a(6)
      Z=a(7)
      NW = A(9)
      alfa=a(10)
      RO1=a(11)

!      if (nint(Z).ne.IZ .or. IZ.ne.nznuc) stop 'Wrong Z in corratom.f'
      
      do j = 1, NW
c$$$         CALL TAPE(R5,XN,ZERO,IZCH,NE,IN,2)
         CALL TAPE(R5,XN,ZERO,2,NE,IN,2)
         n = nint(XN)
         l = nint(r5(2))
         E(l, n) = r5(1)
      call JNT(R8,R5,RO1,BET,ALFA,R7,H,R2,RM,IS,NE,NE5,NE5IS,KT,IN,KSU)
         do i = 1, kt
!            write(100+10*l+n,'(F9.4,E12.4)') r7(i),r2(i)
            wa(i, l, n) = r2(i)
         end do

C Significant number of radial points

         maxNR(l,n)=kt
         do i = kt,1,-1
            if(abs(wa(i, l, n)).gt.1e-10)  then
               maxNR(l,n)=i
               goto 111
            end if
         end do
 111     write(6,'(A1,I1,2A1,F9.4,A,I4)')'E',n,hh(l),'=',E(l, n),
     :      ' Max NR=', maxNR(l,n)

      end do
      close(2)
c$$$      open(2,file='test')
c$$$      print*,'Unit 2 opened'
c$$$      close(2)

!      write(6,'(/A,I1,A/)') ' Overlap with ',nmin(0),'s orbital'
      write(6,'(/A,I1,A1,A/)')'Overlap with ',nmin(lin),hh(lin),
     >   ' orbital'

      call OVER(wa(1,lin,nmin(lin)),wa(1,lin,nmin(lin)), overlap)
      print'(A,I1,A1,A,I1,A1,A,F9.4)',
     >   '<',nmin(lin),hh(lin),'|',nmin(lin),hh(lin),'> ', overlap
      if (abs(overlap-1.0).gt.1e-2) stop 'overlap not 1'

c$$$      write(6,'(/A/)') ' Overlap with 1s orbital'
c$$$
c$$$      call OVER(wa(1,0,1),wa(1,0,1), overlap)
c$$$      print*, '<1s~|1s~> ', overlap

c$$$  if (abs(overlap-1.0).gt.1e-2) stop 'overlap not 1'
      
      call OVER(wa(1,0,1),wa(1,0,2), overlap)
      print*, '<1s~|2s~> ', overlap
      
      call OVER(wa(1,0,1),wa(1,0,3), overlap)
      print*, '<1s~|3s~> ', overlap
      
!      call LENGTH(wa(1,0,1),wa(1,1,2),0,1,resl)
!      print*, '<1s|r|2p> ', resl
      call VELOCITY(wa(1,0,1),wa(1,1,2),0,1,resv)
      print*, '<1s|r|2p> ', resv

      
c Perturbation theory for CI
c Updated for calcium 17.05.2000
c Updated for neon     6.05.2003
 
      if(meta.eq.0 .and. lin.eq.0) then
      write(6,'(/A/)') ' Perturbation theory for CI'
      do L=0,jmax
!         do n=l+1, nmax(l)
         do n=nmin(l), nmax(l)
      call MATRcorr(wa(1,L,n),wa(1,L,n),wa(1,0,nmin(0)),
     :                              wa(1,0,nmin(0)),L,res)
!            if(l.eq.0 .and. n.eq.1) then
            if(l.eq.0 .and. n.eq.nmin(0)) then
               write(6,'(I1,A1,F10.4,F9.4)') n, hh(l), CI(L,n,n), E(L,n)
               goto 112
            end if
            write(6,'(I1,A1,F10.4,F9.4,2F10.4)') n, hh(l), CI(L, n,n),
     :   E(L,n),res,(-1)**l*res/(E(L,n)-E(0,nmin(0)))/dsqrt(dble(2*l+1))
 112        continue
         end do
      end do
      end if

 1954 format(/,'Some parameters',/,'SOST=',E10.4,/,'RM=',E10.4,/,
     :       'NE=',E10.4,/,'H0=',E10.4,/,'BET=',E10.4,/,
     :       'EPS=',E10.4,/,'Z=',E10.4,/,'Z=',E10.4,/,
     :       'NW=',E10.4,/,'ALFA=',E10.4,/,'RO1=',E10.4,/,4(1x,e10.4))

      
      return
      end

      subroutine LABEL(aa, n, l)
      
*  DECLARATION STATEMENTS

*      IMPLICIT REAL*8(A-H, O-Z)
*      IMPLICIT INTEGER*4(I-N)

      character*2 aa

      if(aa(1:1) .eq. '0') n = 0
      if(aa(1:1) .eq. '1') n = 1
      if(aa(1:1) .eq. '2') n = 2
      if(aa(1:1) .eq. '3') n = 3
      if(aa(1:1) .eq. '4') n = 4
      if(aa(1:1) .eq. '5') n = 5
      if(aa(1:1) .eq. '6') n = 6
      if(aa(1:1) .eq. '7') n = 7
      if(aa(1:1) .eq. '8') n = 8
      if(aa(1:1) .eq. '9') n = 9
      if(aa(1:1) .eq. ':') n = 10
      if(aa(1:1) .eq. ';') n = 11
      if(aa(1:1) .eq. '<') n = 12
      if(aa(1:1) .eq. '=') n = 13
      if(aa(1:1) .eq. '>') n = 14
      if(aa(1:1) .eq. '?') n = 15

      
      if(aa(2:2) .eq. 's') l = 0
      if(aa(2:2) .eq. 'p') l = 1
      if(aa(2:2) .eq. 'd') l = 2
      if(aa(2:2) .eq. 'f') l = 3
      if(aa(2:2) .eq. 'g') l = 4
      if(aa(2:2) .eq. 'h') l = 5

      if(aa(2:2) .eq. 'S') l = 0
      if(aa(2:2) .eq. 'P') l = 1
      if(aa(2:2) .eq. 'D') l = 2
      if(aa(2:2) .eq. 'F') l = 3
      if(aa(2:2) .eq. 'G') l = 4
      if(aa(2:2) .eq. 'H') l = 5


      RETURN
      END

      SUBROUTINE MATRcorr (fu1, fu2, fu3, fu4, LV,result)

C Calculates radial part of the matrix element with    
C Hylleraas ground state 1s^2

C
C   < f1, f2 || V  || f3, f4 > ~ R  (f1,f3; f2,f4)
C                Lv               Lv
C
C Functions f1-f4 should be free of any Simpsons wheigts

      include 'par.f'

      common /meshrr/ meshr,rmesh(maxr,3)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     :   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
           
      DIMENSION fu1(maxr), fu2(maxr)
      DIMENSION fu3(maxr), fu4(maxr)
      DIMENSION fun(maxr), temp(maxr)

      if(Lv.lt.0)then !No calculation for Lv<0
         result=0.
         return
      end if

      minfun = 1
      maxfun = meshr
      do i = 1, meshr
         fun(i) = fu1(i) * fu3(i) * rmesh(i,3)
      end do

c$$$      do i = 1, meshr
c$$$         write(201,'(F9.4,E12.4)') rmesh(i,1),rpow1(i,lv)
c$$$         write(202,'(F9.4,E12.4)') rmesh(i,1),rpow2(i,lv)
c$$$      end do
                                          


      call form(fun,minfun,maxfun,rpow1(1,lv),rpow2(1,lv),
     :      minrp(lv),maxrp(lv),meshr,temp,i1,i2)
                       
      result = 0.d0
      DO I = 1, meshr
         result = result + temp(i) * fu2(i) * fu4(i) * rmesh(i,3)
      end do
      
      RETURN
      END


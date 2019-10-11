      subroutine b2MATR(lg,gk,npk,minchil,chil,nchtop)

c Calculate 2nd Born matrix elements in dipole approximation
      
      include 'par.f'

      common/meshrr/ meshr,rmesh(maxr,3)

      COMMON/BM2/ br2(0:lamax,kmax,nchan),bv2(0:lamax,kmax,nchan)
      COMMON/MJ/  MJ

      COMMON/channel/ iopen(nchan)    !Indicates open channels (=1) 
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
           
      dimension npk(nchtop+1),psi(maxr)
      dimension gk(kmax,nchan)
      dimension chil(meshr,npk(nchtop+1)-1), minchil(npk(nchtop+1)-1)

      character*1 h(0:5)         !Spectroscopic labels
      data h /'s','p','d','f','g','h'/

      if (ntype.eq.0) then
         print*, '2nd Born is not implemented for Hylleraas'
         return
      end if

      if(nznuc.eq.1)GS  =1.05544 !H-   ion GS energy in Ry
      if(nznuc.eq.2)GS  =  5.807 !He  atom GS energy in Ry
      if(nznuc.eq.3)GS  = 14.559 !Li+  ion GS energy in Ry
      if(nznuc.eq.8)GS  =118.312 !O^6+ ion GS energy in Ry
      if(nznuc.eq.20)GS =1.3218  !CaIII-CaI   energy in Ry Radzig, Smirnov

      pi = acos(-1.)
      deg = 180/pi
      const=8.067
      Ry = 13.6058

!      MJ = 2

!      if(lg.eq.0) MJ=0
!      if(MJ.gt.lg) STOP "WRONG MJ"

      
      print '(i3,":",$)', nchtop
      write(6,1000)  
!      call QRGRIDn(lg,MJ)                                   


C  The following directives are for the IBM
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& SHARED(lg,nchtop,br2,bv2,GS,ntype,iopen,pi)
C$OMP& SHARED(npk,chil,gk)
C$OMP& PRIVATE(r0,v0,r2,kp,q,k,a)
C$OMP& PRIVATE(nch,nt,psi,maxpsi,ea,la,na,l,nqm,EI)
      do nch = 1, nchtop
C$OMP critical(print)
         print '(i3,$)', nch
C$OMP end critical(print)

!         IOPEN(nch)=1
         call  getchinfo (nch, nt, lg, psi, maxpsi, ea, la, na, l)
         nqm = npk(nch+1) - npk(nch)
         EI = GS + ea
         
!         write(6,1954) nch,EI,na,h(la),h(l)
!         write(6,170)


         do k = 1, nqm
            kp = npk(nch) + k - 1
            q=gk(k,nch)
             if (lg.ge.0 .and. lg.le.2 ) then
!              call RADIAL2x(psi,chil(1,kp),la,l,lg,r2)
               call RADIAL2 (psi,chil(1,kp),la,l,lg,r0,v0)
!               call INTEGRAL2x(res)
!               print'(A,10F9.4)', 'Ratio', r2,r0,res,r2/r0/res
C Normalization
               a=sqrt(pi/2d0)   !Bound states
               if (q.gt.0d0) a=sqrt(pi*q)
               br2 (lg,k,nch)=r0/a 
               bv2 (lg,k,nch)=v0/a 
!               br2 (lg,k,nch)=r2/res/a 
!               write(6,180) q, r0, v0
            else
               br2(lg,k,nch)=0. 
               bv2(lg,k,nch)=0. 
            end if
         end do
      end do
C$OMP END PARALLEL DO 
      write(6,'(A//)') '       '

      
      

 1000 format(//21x,'2nd BORN MATRIX ELEMENTS'/,21x,24('-'))
 1954 format(//11x,'CHANNEL ',I2/,
     >   11x,'Ion energy',F7.3/,
     >   11x,'Bound electron',I2,A1/,
     >   11x,'Free  electron k', A1)

 170  FORMAT(///'Momentum-',2x,'2nd Born amplitude'/,
     : 'PQ number '/,40('-'))

 180  FORMAT(F8.4, 1x, 2(2(E11.4),1x), 1x, 2(2(E11.4),1x))
        

      return
      end


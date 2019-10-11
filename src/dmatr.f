      subroutine dMATR(gk,npk,minchil,chil,nchtop)

* Calculate dipole matrix elements

      include 'par.f'

      common/meshrr/ meshr,rmesh(maxr,3)

      COMMON/dipole/  dr (kmax,nchan),dv (kmax,nchan),dx (kmax,nchan)
      COMMON/dipole1/ dr1(kmax,nchan),dv1(kmax,nchan),dx1(kmax,nchan)
      COMMON/channel/ iopen(nchan)    !Indicates open channels (=1) 
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
           
      dimension npk(nchtop+1),psi(maxr)
      dimension gk(kmax,nchan)
      dimension chil(meshr,npk(nchtop+1)-1), minchil(npk(nchtop+1)-1)

      character*1 h(0:3)         !Spectroscopic labels
      data h /'s','p','d','f'/

c DPI of Lithium       
      common /LITHIUM/ lithium

C  The following are experimental energies
      if(nznuc.eq.2)GS  =  5.807 !He  atom GS energy in Ry
      if(nznuc.eq.3)GS  = 14.559 !Li+  ion GS energy in Ry
      if(nznuc.eq.4) GS = 2.0236 !BeIII-BeI   energy in Ry Radzig, Smirnov
      if(nznuc.eq.8)GS  =118.312 !O^6+ ion GS energy in Ry
      if(nznuc.eq.1)GS  =1.05544 !H-   ion GS energy in Ry
      if(nznuc.eq.12)GS =1.6670  !MgIII-MgI   energy in Ry Radzig, Smirnov
      if(nznuc.eq.20)GS =1.3218  !CaIII-CaI   energy in Ry Radzig, Smirnov

      pi = dacos(-1d0)
      const=8.067
            
* Dipole channels only

      lg = 1
      write(6,1000)  

      print '(i3,":",$)', nchtop
C  The following directives are for the SGI
! c$doacross local(nch,nt,psi,maxpsi,ea,la,na,l,nqm,EI)
c$& local(r1,v1,x1,r,v,x,kp,q,k,a)
c$& share(lg,nchtop,dr,dv,dx,dr1,dv1,dx1,GS,ntype) 
C  The following directives are for the SUN
C$par DOALL SHARED(lg,nchtop,dr,dv,dx,dr1,dv1,dx1,GS,ntype,iopen,pi)
C$par& SHARED(const,npk,chil,gk)
C$par& PRIVATE(r1,v1,x1,r,v,x,kp,q,k,a,omega,cl,cv,cx,c)
C$par& PRIVATE(nch,nt,psi,maxpsi,ea,la,na,l,nqm,EI) SCHEDTYPE(SELF(1))
C  The following directives are for the IBM
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& SHARED(lg,nchtop,dr,dv,dx,dr1,dv1,dx1,GS,ntype,iopen,pi)
C$OMP& SHARED(const,npk,chil,gk)
C$OMP& PRIVATE(r1,v1,x1,r,v,x,kp,q,k,a,omega,cl,cv,cx,c)
C$OMP& PRIVATE(nch,nt,psi,maxpsi,ea,la,na,l,nqm,EI)
      do nch = 1, nchtop
C$OMP critical(print)
         print '(i3,$)', nch
C$OMP end critical(print)
         call update(6)
         IOPEN(nch)=1
         call  getchinfo (nch, nt, lg, psi, maxpsi, ea, la, na, l)
         nqm = npk(nch+1) - npk(nch)
         EI = GS + ea
         
c$$         write(6,1954) nch,EI,na,h(la),h(l)
c$$         write(6,170)

c$$         do ii=1,maxpsi
c$$            write(500+nch,'(2E13.4)') rmesh(ii,1), psi(ii)
c$$         end do

c$$         do ii=1,maxpsi
c$$         write(600+nch,'(2E13.4)')
c$$     :         rmesh(ii,1),chil(ii,1)/rmesh(ii,3)
c$$         end do

         do k = 1, nqm
            dr (k,nch)=0.0
            dv (k,nch)=0.0
            dx (k,nch)=0.0
            dr1(k,nch)=0.0
            dv1(k,nch)=0.0     
            dx1(k,nch)=0.0     
            kp = npk(nch) + k - 1
            q=gk(k,nch)
c$$$            if (q.eq.0.0) cycle
!            if (q.eq.0.0) PRINT*,'CAUTION: Q=0'
            if (ntype.eq.0) then
c$$$               call pdipole(psi,chil(1,kp),la,l,r1,v1,x1,r,v,x) ! for Lesech
               call xdipole(psi,chil(1,kp),la,l,r1,v1,x1,r,v,x)
            else if (ntype.eq.-1) then 
               call hDIPOLE(psi,chil(1,kp),la,l,lg,r,v,x)
               r1 = 0.0
               v1 = 0.0
               a1 = 0.0
            else if (ntype.eq.1 .and. lithium.eq.1) then 
               call LiDIPOLEc(nt,chil(1,kp),la,L,lg,r,v,x)
            else if (ntype.eq.1 .and. lithium.ne.1) then 
               call multiDIPOLE(psi,chil(1,kp),la,l,r1,v1,x1,r,v,x)
            endif 
* Restore normalization

            a=sqrt(pi/2d0) !Bound states
            if (q.gt.0d0) a=sqrt(pi*q)
            if(k.eq.1 .and. q. lt. 0.0) IOPEN(nch)=0 !Channel closed
            r  = r/a
            v  = v/a
            x  = x/a
            r1 = r1/a
            v1 = v1/a
            x1 = x1/a

            dr (k,nch)=r 
            dv (k,nch)=v 
            dx (k,nch)=x 
            dr1(k,nch)=r1
            dv1(k,nch)=v1      
            dx1(k,nch)=x1      

c$$$            energy = q**2
c$$$            if (q.lt.0d0) energy=-q**2 !Bound states
c$$$            omega = energy + EI       
c$$$            cl=   const*2/3*omega
c$$$            cv= 4*const*2/3/omega
c$$$            cx=16*const*2/3/omega**3
c$$$            csl1 = cl * r1**2
c$$$            csv1 = cv * v1**2
c$$$            csx1 = cx * x1**2
c$$$            csl  = cl * r**2
c$$$            csv  = cv * v**2
c$$$            csx  = cx * x**2
c$$$            write(6,180) q, r1,v1,r,v,
c$$            write(6,180) q, r1,v1,x1,v,
c$$     >                         csl1,csv1,csx1,csl,csv,csx
         end do
      end do
C$OMP END PARALLEL DO



 1000 format(//21x,'PHOTOIONIZATION MATRIX ELEMENTS'/,21x,31('-'))
 1954 format(//11x,'CHANNEL ',I2/,
     >   11x,'Ion energy',F7.3/,
     >   11x,'Bound electron',I2,A1/,
     >   11x,'Free  electron k', A1)

 170  FORMAT(///'Momentum-',7x,'Ionization amplitude',20x, 
     : 'Photoionization cross-section'/,
     : 'PQ number ',6x,'No GSC',12x,'GSC',
     :             16x,'No GSC',23x,'GSC'/,
     : 14x,2('r',7x,'v',8x),2x,2('r',7x,'v',8x,'a',9x)/,99('-'))  

 180  FORMAT(F8.4, 1x, 2(2(F8.4),1x), 1x, 2(3(F9.4),1x))
        

      return
      end


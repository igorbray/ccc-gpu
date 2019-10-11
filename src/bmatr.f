      subroutine bMATR(lg,gk,npk,minchil,chil,nchtop)

* Calculate Born matrix elements

      include 'par.f'

      common/meshrr/ meshr,rmesh(maxr,3)

      COMMON/BORN/ br(0:lamax,kmax,nchan),bv(0:lamax,kmax,nchan)
      COMMON/channel/ iopen(nchan)    !Indicates open channels (=1) 
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
           
      dimension npk(nchtop+1),psi(maxr)
      dimension gk(kmax,nchan)
      dimension chil(meshr,npk(nchtop+1)-1), minchil(npk(nchtop+1)-1)
      COMMON /moment/ qq,omega
      COMMON /KNMT/   E0,E1,theta
      COMMON /const/  pi, Ry, rad, SIP

      character*1 h(0:5)         !Spectroscopic labels
      data h /'s','p','d','f','g','h'/

      complex Hr(0:lamax,nchan),Hv(0:lamax,nchan)
      
      if(nznuc.eq.1)GS  =1.05544 !H-   ion GS energy in Ry
      if(nznuc.eq.2)GS  =  5.807 !He  atom GS energy in Ry
      if(nznuc.eq.3)GS  = 14.559 !Li+  ion GS energy in Ry
      if(nznuc.eq.8)GS  =118.312 !O^6+ ion GS energy in Ry
      if(nznuc.eq.20)GS =1.3218  !CaIII-CaI   energy in Ry Radzig, Smirnov

      pi = acos(-1.)
      deg = 180/pi
      rad = pi / 180.
      const=8.067
      Ry = 13.6058
      SIP=24.59       !Single ionization potential 

C Entered momentum transfer
      
      qq = 0.1290       !E1=10eV n=1
      qq = 0.3964       !E1=75eV n=1
      qq = 0.1343       !E1=5eV  n=1
      qq = 0.441        !Schlemmer Eo=400 E2=10eV  4deg
      qq = 0.955        !Schlemmer Eo=400 E2=10eV  10deg
      qq = 0.4356       !E1=75eV n=2   
      qq = 0.2241       !E1+E2=87eV
      qq = 0.1699       !E1=10eV OCW normalization
      qq = 0.447        !E1+E2=99eV E1=1000 eV
      qq = 0.298        !Azzedine 11eV
      qq = 0.4556       !Rouvellou n=1
      qq = 0.4412       !Franz & Altic
      qq = 0.4356       !E1=75eV n=2
      qq = 0.7420       !E1=40eV n=2   E0= 570eV 4deg
      qq = 0.2400       !E1+E2=99eV
      qq = 0.7994       !E1=20eV n=2   E0=1500eV 4deg
      qq = 0.1779       !E1=5eV  n=2
      qq = 0.6239       !E1=10eV n=2   E0= 570eV 4deg
      qq = 0.6784       !E1=9.25eV n=2   E0= 365.8eV 4.5deg
      QQ = 0.01         !Photo limit
      
c$$$      OPEN(UNIT=9,  FILE='kinematics',STATUS= 'unknown')
c$$$      READ(9, *) E0,E1,theta
c$$$
c$$$      E2 = E0-E1-SIP
c$$$      fk0 = sqrt(E0/Ry)
c$$$      fk1 = sqrt(E1/Ry)
c$$$      fk2 = sqrt(E2/Ry)
c$$$      QQ  = sqrt(fk0**2 + fk1**2 - 2*fk0*fk1*cos(theta/deg))
c$$$      thq = -asin(fk0/QQ*sin(theta/deg))
c$$$      write(6, 110) E0, E1, E2, fk0,fk1,fk2
c$$$      write(6,120) theta, QQ, thq*deg
c$$$      close(UNIT = 9)
c$$$
      
      write(6,1000)  
      print '(i3,":",$)', nchtop
C  The following directives are for the SGI
! c$doacross local(nch,nt,psi,maxpsi,ea,la,na,l,nqm,EI)
c$& local(r1,v1,x1,r,v,x,kp,q,k,a,energy)
C$& local(csl1,csv1,csl,csv)
c$& share(lg,nchtop,br,bv,dx,br1,bv1,dx1,GS,ntype,const,qq) 
C  The following directives are for the SUN
C$par DOALL SHARED(lg,nchtop,br,bv,dx,br1,bv1,dx1,GS,ntype,iopen,pi)
C$par& SHARED(const,npk,chil,gk,qq)
C$par& PRIVATE(r1,v1,x1,r,v,x,kp,q,k,a,omega,cl,cv,cx,c,energy)
C$par& PRIVATE(nch,nt,psi,maxpsi,ea,la,na,l,nqm,EI) SCHEDTYPE(SELF(1))
C$par& PRIVATE(csl1,csv1,csl,csv)
C  The following directives are for the IBM
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& SHARED(lg,nchtop,br,bv,dx,br1,bv1,dx1,GS,ntype,iopen,pi,qq)
C$OMP& SHARED(const,npk,chil,gk)
C$OMP& PRIVATE(r1,v1,x1,r,v,x,kp,q,k,a,omega,cl,cv,cx,c)
C$OMP& PRIVATE(csl1,csv1,csl,csv)
C$OMP& PRIVATE(nch,nt,psi,maxpsi,ea,la,na,l,nqm,EI,energy)
      do nch = 1, nchtop
C$OMP critical(print)
         print '(i3,$)', nch
C$OMP end critical(print)
         call update(6)
         IOPEN(nch)=1
         call  getchinfo (nch, nt, lg, psi, maxpsi, ea, la, na, l)
         nqm = npk(nch+1) - npk(nch)
         EI = GS + ea
         
c$$$         write(6,1954) nch,EI,na,h(la),h(l)
c$$$         write(6,170)

         do k = 1, nqm
            kp = npk(nch) + k - 1
            q=gk(k,nch)
            if (ntype.eq.0) then
               call xBORN(psi,chil(1,kp),la,l,lg,QQ,r1,v1,r,v)
!               call pBORN(psi,chil(1,kp),la,l,lg,QQ,r1,v1,r,v) !Lesech
            else
               call polyBORN(psi,chil(1,kp),la,l,lg,QQ,r1,v1,r,v)
            end if
            
C Restore normalization

            a=sqrt(pi/2d0) !Bound states
            if (q.gt.0d0) a=sqrt(pi*q)
            if(k.eq.1 .and. q. lt. 0.0) IOPEN(nch)=0 !Channel closed
            r  = r/a
            v  = v/a
            r1 = r1/a
            v1 = v1/a

            br (lg,k,nch)=r 
            bv (lg,k,nch)=v 

            energy = q**2
            if (q.lt.0d0) energy=-q**2 !Bound states
            om = energy + EI
            if(k.eq.1) omega=om
            cl=   const*2/3*om
            cv= 4*const*2/3/om
            csl1 = cl * r1**2
            csv1 = cv * v1**2
            csl  = cl * r**2
            csv  = cv * v**2
!            write(6,180) q, r1,v1,r,v,csl1,csv1,csl,csv
         end do
      end do
C$OMP END PARALLEL DO
      write(6,'(A//)') '       '

      
      

 110  format(//'            KINEMATICS OF THE REACTION',
     :       / '            --------------------------',
     :      // '                  PROJECTILE          ',
     :      // '                  Incoming          Outgoing',
     :       / '                  electron          electrons',
     :       / '                                Fast        Slow',
     :       / ' ===============================================',
     :      // ' Energy,   eV     ',    3(F7.2,5X),
     :       / ' Momentum, a.u.   ',    3(F7.4,5X)//)
      
 120  format(/ '             TRANSFERRED MOMENTUM      ',
     :      / '             --------------------      ',
     :      / 'Scattering    Momentum   Transfered momentum',
     :      / 'angle, deg   transfered      angle, deg    '/,
     :         F10.3,2x,F10.4,2x,F10.3,2x)
 1000 format(//21x,'BORN MATRIX ELEMENTS'/,21x,20('-'))
 1954 format(//11x,'CHANNEL ',I2/,
     >   11x,'Ion energy',F7.3/,
     >   11x,'Bound electron',I2,A1/,
     >   11x,'Free  electron k', A1)

 170  FORMAT(///'Momentum-',7x,'Ionization amplitude',13x, 
     : 'Photoionization cross-section'/,
     : 'PQ number ',6x,'No GSC',12x,'GSC',
     :             15x,'No GSC',13x,'GSC'/,
     : 14x,2('r',7x,'v',8x),2x,2('r',7x,'v',9x)/,81('-'))  

* 180  FORMAT(F8.4, 1x, 2(2(F8.4),1x), 1x, 2(3(F9.4),1x))
 180  FORMAT(F8.4, 1x, 2(2(E11.4),1x), 1x, 2(2(E11.4),1x))
        

      return
      end


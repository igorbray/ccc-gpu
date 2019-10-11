      subroutine CROSS(Lstart,Lstop)                                   
      real*8 rylm
C This subroutine calculates the TDCS from the shell b  
C in Born approximation in length and velocity forms as      
C a function of ejected electron angle.
C Coplanar geometry only - Phi=0
C Ejected electron is in the scattering plane
      
C      (3)      k' 1     q     2          
C  Sigma  =  8 --- -4  |T(b,k)|           
C               ko q                      
C                                     
C                                                       
C  q            iPhi(l) J-l          2    / l J lb \  r,v    
C T(b,k) = SUM e       i   Y(TH)  hat(J)  |        | D(El,i) 
C          Jlm               lm            \-m 0 m /  L      
                                
      include 'par.f'                                         
            
      COMMON/CONT/ CE(ncmax),CINT(ncmax),NCSTATES,ENERGY
      COMMON/channel/ iopen(nchan) !Indicates open channels (=1) 
      COMMON/BORN/ br(0:lamax,kmax,nchan),bv(0:lamax,kmax,nchan)
      COMMON/BORN2r/ Hr(0:lamax,nchan)
      COMMON/BORN2v/ Hv(0:lamax,nchan)
      COMMON /moment/ qq,omega
      COMMON /faza/ ph(0:lamax,nchan)
      dimension nopen(lamax+1,0:lamax)

      complex Tr,Tv,ci,F,br,bv
      complex r,v
      dimension CPr(lamax+1,0:lamax), CPv(lamax+1,0:lamax)
      dimension Cr(lamax+1), Cv(lamax+1)

      complex phase,ph
      dimension psi(maxr)
      complex Hr,Hv

      REAL*8 threej


      character*1 h(0:10) 
      data h / 's','p','d','f','g','h','n','o',
     >          'x','x','x'/
C Kinematic variables

*      data Eo,Ea,theta  /5640.41,5500.,-1.00/ !Eb=75eV, n=2
*      data Eo,Ea,theta  /5599.,5500.,-0.45/ !Ea+Eb=20 eV
      data Eo,Ea,theta  /5570.41,5500.,-0.35/ !Eb=5ev, n=2
*      data Eo,Ea,theta  /5529.59,5500.,-0.35/ !Eb=5eV, n=1
*      data Eo,Ea,theta  /5534.59,5500.,-0.32/ !Eb=10eV, n=1
*      data Eo,Ea,theta  /5599.59,5500.,-1.00/ !Eb=75eV, n=1
      pi = acos(-1.)
      deg= 180./pi
      SIP=24.59                 
      Ry = 13.6058
      
      Eb = Eo-Ea-SIP
      fko = sqrt(Eo/Ry)
      fka = sqrt(Ea/Ry)
      fkb = sqrt(Eb/Ry)
      FQ  = sqrt(fko**2 + fka**2 - 2*fko*fka*cos(theta/deg))  
      thb = -asin(fka/FQ*sin(theta/deg))      
      Eb = Eo-Ea-24.59
      write(6, 110) Eo, Ea, Eb, fko,fka,fkb
      write(6,120) theta, FQ, Thb*deg
      if(abs(Eb-ENERGY).gt.0.01) STOP "Energy is not conserved"
      if(abs(qq-fq).gt. 1e-2)
     :   write(6,'(A,F9.4)') 'Entered      ', QQ
      
      THo = 0.                  !- initial angle
      NTT=60                     !- number of points
      DT=2*pi/NTT                !- increment 
      
      const=8 * fka/fko /FQ**4  !TDCS constant
      const=const*2             !To compensate WF normalization on energy in Ry

*      print*, 'CONST', const, const*4/omega**2


      ci = (0.,1.)                                         

C Determine size nm,lm. No calculation
         
      nm = 1
      lm = 0
      do J = Lstart, Lstop     
         nch = 0
 10      nch = nch + 1
         call  getchinfo (nch, nt, J, psi, maxpsi, ea, la, na, l)
         if (nch.ne.0) then
            EI = GS + ea
            if (na.gt.nm) nm = na
            if (la.gt.lm) lm = la
            nopen(na,la)=0
            if(iopen(nch).eq.1) nopen(na,la)=1
            phase = ph(J,nch)
            phi = atan(aimag(phase)/real(phase))
*            write(6,1954) nch,EI,na,h(la),h(l),omega*Ry,phi
            go to 10
         end if
      end do

*      print*, 'nmax,lmax', nm,lm

C TDCS calculation
C Ejected electron angle loop   

      write(6,170)              !Title for TDCS
      
      do k=0,NTT
         TH = THo+k*DT
         TH = TH -THb           !Born angle transformation
*         write(6,'(A,F9.2)') 'ANGLE=', TH*deg


C Initialization

         do nb=1,nm
            Cr(nb) = 0.
            Cv(nb) = 0.
            do lb=0,lm
               CPr(nb,lb) = 0.
               CPv(nb,lb) = 0.
            end do
         end do
         
C Ion final state loop

         do nb=1,nm
            do lb=0,lm
               if(nopen(nb,lb) .eq. 0) goto 111
C Angular momentum projection loop
               
               do mb= -lb,lb
                  mm= iabs(mb)
                  mp= (mm+mb)/2

                  Tr = (0.,0.)
                  Tv = (0.,0.)

C Multipoles loop                  
                  do J = Lstart, Lstop     
                     nch = 0
 11                  nch = nch + 1
                     call  getchinfo (nch,nt,J,psi,maxpsi,ea,la,na,l)
                     if (nch.ne.0)then
                        if( na.eq.nb .and. la.eq.lb) then
*       write(6,1954) nch,EI,na,h(la),h(l),omega*Ry,phi
                           r = Hr(j,nch)
                           v = Hv(j,nch)
*      write(6,'(I1,A1,A2,A1,4(1p,E10.2)/)') nb,h(lb),'->',h(l),r,v
                           call IOT3(L,j,lb,-mb,0,mb,threej)             
                           ang= hat(j)**2 * threej                    
                           Y = RYLM(L,mm,dble(TH)) * (-1)**mp  
                           F = (ci)**(j-L)
                           Tr = Tr + F*Y*ang*r
                           Tv = Tv + F*Y*ang*v
*      write(6,'(2(A,I2),A,F9.4,A,F9.4,A,6E11.4)') 'mb=', mb, 
*     :  ' L=', L, ' ANG',ANG, ' Y=',Y,' Tv=', v, Tv
                        end if
                        go to 11
                     end if
                  end do               !J
                  CPr(nb,lb)= CPr(nb,lb) + abs(Tr)**2           
                  CPv(nb,lb)= CPv(nb,lb) + abs(Tv)**2           
               end do                  !mb 
               CPr(nb,lb)= CPr(nb,lb)* const           
               CPv(nb,lb)= CPv(nb,lb)* const*4/omega**2
*      write(6,'(A,I1,A1,A,2E11.4/)')'CP(',nb,h(lb),') ',CPr(nb,lb),
*     :            CPv(nb,lb)
               Cr(nb)= Cr(nb)+CPr(nb,lb)
               Cv(nb)= Cv(nb)+CPv(nb,lb)
               write(100+10*lb+nb,'(F6.1,2E12.4)')
     :                     (TH+THb)*deg, CPr(nb,lb), CPv(nb,lb)
 111           continue
            end do                     !lb 
*      write(6,'(A,I1,A,2E11.4/)') ' C(n=',nb,')',Cr(nb),Cv(nb)
            write(200+nb,'(F6.1,2E12.4)') (TH+THb)*deg, Cr(nb), Cv(nb)
         end do                        !nb
      end do                           !TH
      write(6,'(A//)')'   '
               
*      print*, 'CONST=', const*4/omega**2
         
 100  format(///10x,'TRIPLE DIFFERENTIAL CROSS SECTION CALCULATION'/,
     >       10x,45('-')// )
 110  format(//'            KINEMATICS OF THE REACTION',
     :       / '            --------------------------',
     :       / '         The Coplanar geometry of e2e setup',
     :      // '                  Incoming          Outgoing',
     :       / '                  electron          electrons',
     :       / '                                Fast        Slow',
     :       / ' ===============================================',
     :      // ' Energy,   eV     ',    3(F7.2,5X),
     :      // ' Momentum, a.u.   ',    3(F7.4,5X)//)
      
      
 120  format(/ '             TRANSFERRED MOMENTUM      ',
     :      / '             --------------------      ',
     :      / 'Scattering    Momentum   Transfered momentum',
     :      / 'angle, deg   transfered      angle, deg    '/,
     :         F10.3,2x,F10.4,2x,F10.3,2x)

 170    FORMAT(//' ANGLE',13X,'TDCS'/,
     >           ' DEGREE        r',12x,'v',12x/,35('-'))
 171    FORMAT(//'BOUND',9X,' GENERALIZED  OSCILLATOR  STRENGTH'/,
     >           'STATE       Direction  -Q',14x,'Direction  +Q',15x/,
     >    12x,2('r',12x,'v',12x)/,55('-'))
 200  format(1x,I1,A1, 2x,12(E11.4,2x))
 201  format(F9.2, 2x,12(F11.2,2x))
 202  format(F9.2, 2x,12(E11.3,2x))
 1954 format(//11x,'CHANNEL ',I2/,
     >   11x,'Ion energy',F7.3/,
     >   11x,'Bound electron',I2,A1/,
     >   11x,'Free  electron k', A1/,
     >   11x,'Photon energy', F7.3/,
     >   11x,'Coulomb phase', 2F7.3)

      return
      end


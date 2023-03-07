      subroutine PHOTO(gk,npk,wk,vt,nchtop,ovlpn,csfile,slowery,phasel,
     >   phaseq)
                                                             
C Integral form of the photoionization cross section                           
      use photo_module
      include 'par.f'
      parameter (knm6 = knm * 6, ncs = (lamax+1)*6)
                                                             
c$$$      COMMON/dipole/  dr (kmax,nchan),dv (kmax,nchan),dx (kmax,nchan)
      COMMON/channel/ iopen(nchan)    !Indicates open channels (=1) 
      COMMON/even/ ieven(nchan)    !Indicates even bound state channels (=0)
      COMMON /gstspin/   iSpin, meta
                                                             
      character*3 chan(knm), chs(0:lamax)
      common /charchan/ chan
      common /chanen/ enchan(knm)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      character lh(0:10), csfile*(*), file*50
      data lh / 's','p','d','f','g','h','n','o',
     >          'x','x','x'/
      integer npk(nchtop+1)
!      dimension gk(kmax,nchan), temp(maxr), ovlpn(knm), cs(knm,3),
!     >   tcs(0:lamax,3),tics(0:lamax,3),tnbcs(0:lamax,3)
      dimension gk(kmax,nchan), temp(maxr), ovlpn(knm), cs(knm,6),
     >   tcs(0:lamax,6),tics(0:lamax,6),tnbcs(0:lamax,6)
      complex wk(kmax*nchan)             !Greens function
      complex vt(npk(nchtop+1)-1,nchtop) !T-matrix
      complex GF,T,ci,TR,TV,TX,dTR,dTV,dTX,phase,phasel(kmax,nchan),
!     >   D(KNM,-1:1,3), To, SU
     >   D(KNM,-1:1,6), To, SU, T1,T2,T3, phaseq(knm)
      complex*16 coulphase
      real*8 deta
      data ltop,ry,cs,tnbcs,tics,tcs/0,13.6058,knm6*0.0,
     >   ncs*0.0, ncs*0.0, ncs*0.0/

c DPI of Lithium       
      common /LITHIUM/ lithium
      common /CcoefsE/  E(KNM)
      double precision E
            
      if(nznuc.eq.2)GS  =  5.807 !He  atom GS energy in Ry
!      if(nznuc.eq.3)GS  = 14.559 !Li+  ion GS energy in Ry
      if(nznuc.eq.3 .and. zasym.eq.0)GS  = 6.03/Ry !Li-  ion GS energy in Ry
      if(nznuc.eq.11 .and. zasym.eq.0)GS  = 5.492/Ry !Na-  ion GS energy in Ry
      if(nznuc.eq.8)GS  =118.312 !O^6+ ion GS energy in Ry
      if(nznuc.eq.1)GS  =1.05544 !H-   ion GS energy in Ry
      if(nznuc.eq.4) GS = 2.0236 !BeIII-BeI   energy in Ry Radzig, Smirnov
      if(nznuc.eq.12)GS = 1.6670 !MgIII-MgI   energy in Ry Radzig, Smirnov
      if(nznuc.eq.20)GS = 1.3218 !CaIII-CaI   energy in Ry Radzig, Smirnov
!      if(nznuc.eq.10)GS = 4.596  !NeIII-NeI   energy in Ry Radzig, Smirnov
      if(nznuc.eq.18)GS = 5.464  !ArIII ^1S   energy in Ry by Saloman
      if(nznuc.eq.4  .and. Zasym.eq.3) GS =27.308 !Be^2+
      if(nznuc.eq.5  .and. Zasym.eq.4) GS =44.07 !B^3+
      if(nznuc.eq.12 .and. Zasym.eq.11)GS =273.31 !Mg^10+ 
      if(nznuc.eq.13 .and. Zasym.eq.12)GS =322.49 !Al^11+ 
      if(nznuc.eq.14 .and. Zasym.eq.13)GS =376.28 !Si^12+ 
c$$$      if(nznuc.eq.10)GS = 5.103  !NeIII ^1S   energy in Ry Kilin et al.
      if(nznuc.eq.10)GS = 8.957 !NeIII ^1S   energy in Ry by Kramida and Nave
      if(meta.eq.1) then
         if(iSpin.eq.0) then
            GS   =4.283         !He atom 1s2s ^1S energy in Ry
            print*, 'Meta Singlet'
         else
            GS   =4.350         !He atom 1s2s ^3S energy in Ry
            print*, 'Meta Triplet'
         end if
      end if
      
      if(nznuc.eq.2 .and. ntype.lt.0)
     >               GS  =  3.77 !H2  molecule

      if(lithium.eq.1) GS =  14.955 !Li  atom GS energy in Ry
      if(lithium.eq.1 .and. NZNUC.eq.5) GS =  46.840 !B++
      if(lithium.eq.1 .and. NZNUC.eq.10) GS = 205.60 ! 205.5834 !Ne7+ Ry (all electrons) (NIST total binding energy)
      pi = acos(-1.0)                                       
      ci = cmplx(0.0, 1.0)
      do n = 1, knm
         do l = -1, 1
            do ng = 1, 3
!            do ng = 1, 6
               D(n,l,ng) = (0.0,0.0)
            enddo
         enddo
      enddo 

      open(11,file='tmp1')
      open(12,file='tmp2')
      open(13,file='tmp3')
      open(14,file='tmp4')
      open(20,file='tmp5')
      open(21,file='tmp6')
      open(23,file='tmp7')

      i = 1
      do while (csfile(i:i).ne.'_'.and.csfile(i:i).ne.' ')
         i = i + 1
      enddo 

      if (i.eq.8) then
         file = 'Lphotocs'//csfile(i:60)
         open(42,file=file)
         file = 'Vphotocs'//csfile(i:60)
         open(43,file=file)
         file = 'Aphotocs'//csfile(i:60)
         open(44,file=file)
         file = 'T1photocs'//csfile(i:60)
         open(45,file=file)
         file = 'T2photocs'//csfile(i:60)
         open(46,file=file)
         file = 'T3photocs'//csfile(i:60)
         open(47,file=file)
         file = 'photocs'//csfile(i:60)
         open(22,file=file)
      else 
         file = csfile(1:i-8)//'Lphotocs'//csfile(i:60)
         open(42,file=file)
         file = csfile(1:i-8)//'Vphotocs'//csfile(i:60)
         open(43,file=file)
         file = csfile(1:i-8)//'Aphotocs'//csfile(i:60)
         open(44,file=file)
         file = csfile(1:i-8)//'T1photocs'//csfile(i:60)
         open(45,file=file)
         file = csfile(1:i-8)//'T2photocs'//csfile(i:60)
         open(46,file=file)
         file = csfile(1:i-8)//'T3photocs'//csfile(i:60)
         open(47,file=file)
         file = csfile(1:i-8)//'photocs'//csfile(i:60)
         open(22,file=file)
      endif 

* Dipole channels only

      lg = 1

*      write(6,1000)

* Final channel loop
      
      sigma_l =0.0
      sigma_v =0.0
      sigma_x =0.0
      sigma1_l =0.0
      sigma1_v =0.0
      sigma1_x =0.0
    
      CS2l = 0.
      CS2v = 0.
      CS2x = 0.
      
      ntmax = 0
      do nchf = 1, nchtop
         call getchinfo (nchf,nt,lg,temp,maxpsif,ef,lf,nf,llf)
         EI=GS+ef
         if(lithium.eq.1)    EI=GS+E(nt)*2 !Li  patch     
         if (nchf.eq.1) then
            etot = ef + gk(1,nchf)**2
            om = EI + gk(1,nchf)**2
            do ng = 1, 3
!            do ng = 1, 6
               write(41+ng,"('Photon energy:',f10.5,' eV. Total '
     >            ,'electron energy:',f10.5,' eV.')") om*ry, etot*ry
               if (slowery.gt.0.0) write (41+ng,
     >            "('    L inner  L outer      Re(T)       Im(T) for ',
     >            'inner electron energy of',f6.2,' eV')") slowery*ry
            enddo
         endif 

         nqm = npk(nchf+1) - npk(nchf)
*         write(6,1954) nchf,EI,nf,lh(lf),lh(LLf),nqm            
*         write(6,'(11x,A,F9.4)') 'Bound state energy',ef
         ipseudo=0              !Eigenstates
         if(ef.gt.0.) ipseudo=1 !Pseudostates
         if(IOPEN(nchf).eq.0)then
*            write(6,'(11x,A,//)') 'Channel closed'
            goto 113
         end if
                                                             
         TR =(0.0,0.0)                                          
         TV =(0.0,0.0)                                          
         TX =(0.0,0.0)                                          
         
*         write(6,35)                                            
*         write(6,30) !Print T-matrix                                            

* Initial channel loop
      
         nv = 0
         const=8.067 !  4 pi^2 * alpha * a0^2 * 10^{18}, alpha = 1/137 
!         const = const * 3 ! Ne p^6 patch
         if(lithium.eq.1) const = const * 0.5 ! Li patch
         
         do nchi = 1, nchtop
            nqm = npk(nchi+1) - npk(nchi)
            dTR =(0.,0.)
            dTV =(0.,0.)
            dTX =(0.,0.)
            do k = 1, nqm
               nv = nv + 1
               q0 = gk(1,nchf)
               if (gk(1,nchf).eq.0.0) q0 = 1.0
               q  = gk(k,nchi)
               if (q.le.0.0) then
                  tmp = 1.0
               else
                  tmp = q
               end if 
               if(k.eq.1.and.q.lt.0.0) goto 112
               if(q.gt.0.0)then                              
                  rn=sqrt(2./q)                              
               else                                          
                  rn=1.0                                     
               end if                                        
               r =dr (k,NCHI)*rn                             
               v =dv (k,NCHI)*rn                             
               x =dx (k,NCHI)*rn                             
               GF= wk(nv) * tmp**2
               T =vt(nv,nchf)/q0/tmp
*               write(6,40) q,gf,T,r,v,x
               if(k.eq.1.and.NCHI.eq.NCHF)then   !Diagonal term f-f
                  p =q
                  ro=r
                  vo=v
                  xo=x
                  To= T
                  SU=1-2*pi*p*ci*To
               end if           
               T2 = T
               T3 = T
               if(k.eq.1.and.NCHI.eq.1)then   !electron-scattering ME
                  T1= T
               end if           
               if(k.eq.1.and.NCHI.eq.2)then   !electron-scattering ME
                  T2= T
               end if           
               if(k.eq.1.and.NCHI.eq.5)then   !electron-scattering ME
                  T3= T
               end if           

               dTR=dTR+r*GF*T                                
               dTV=dTV+v*GF*T                                
               dTX=dTX+x*GF*T                                
 112           continue                                         
            end do

!             io = iopen(NCHI)+1
!             io = 0
!             if(nchi.eq.24 .or. nchi.eq.27) io=1

!            io = 0
!            if(ieven(NCHI).ne.ieven(NCHF)) io=1
!            dTR = (-1)**io * dTR
!            dTV = (-1)**io * dTV
!            dTX = (-1)**io * dTX

*            write(6,'(I3,2x,6(E11.4,2x),I3)') NCHI,dTR,dTV,dTX,io  
*            write(6,41) NCHI,dTR,dTV,dTX 

            TR = (TR +  dTR)
            TV = (TV +  dTV)
            TX = (TX +  dTX)
         end do

*         write(6,42) 0, ro,0d0,vo,0d0,xo,0d0
*         write(6,43) TR+ro,TV+vo,TX+xo

         delta= 0.5*atan2(aimag(SU),real(SU))
*         write(6,1960) To                                                    
*         write(6,1962) SU,delta
 1960    FORMAT(/'Test on LS-integral equation'/, 2(2(E11.4,2x)/))
 1962    FORMAT( 'S-matrix'/,2(E11.4,2x),//
     :      'Phase shift respect to Coulomb'/,E11.4/)
                  

         
C Restore normalization

         if (p.eq.0.0) then
            const = const / 2.0
         else 
            const=const*p/2.0
         endif 
         om=EI+p**2

         TR = (TR+ro) 
         TV = (TV+vo) 
         TX = (TX+xo) 
         
         phase = phasel(1,nchf)
         if (ef.gt.0.0) then
            deta = - 2.0 / sqrt(ef)
c$$$            phase = phase * coulphase(deta,lf) * ovlpn(nt)
            phase = phase * phaseq(nt) * ovlpn(nt)
         endif 

         D(nt,llf-lf,1) = TR * sqrt(const*om*2.0/3.0) * phase
         D(nt,llf-lf,2) = TV * sqrt(const*8.0/3.0/om) * phase
         D(nt,llf-lf,3) = TX * sqrt(const*32.0/3.0/om**3) * phase
         D(nt,llf-lf,4) = T1 * phasel(1,1) *                phase
         D(nt,llf-lf,5) = T2 *                              phase
         D(nt,llf-lf,6) = T3 *                              phase

         
*         write(6,'(/A,F9.4/)') 'OMEGA=', om*ry

         CLo=  const*2/3*om*abs(ro)**2                          
         CVo=4*const*2/3/om*abs(vo)**2                          
         CXo=16*const*2/3/(om**3)*abs(xo)**2                          
         CL =  const*2/3*om*abs(TR)**2                          
         CV =4*const*2/3/om*abs(TV)**2                          
         CX=16*const*2/3/(om**3)*abs(TX)**2
         
         if (ef.gt.0.0.and.abs(ef-slowery)/ef.lt.1e-4) then
C  Here to write out the double ionization matrix elements.
C  The asymptotic charge seen by the inner target-space electron is +2. OVLP is
C  the overlap between the EF-energy state with the true continuum function
C  of same energy. This contains the continuum normalisation sqrt(2/pi), but
C  NOT division by sqrt(ef).
!            do ng = 1,3
            do ng = 1,6
               write(41+ng,70) lf,llf, D(nt,llf-lf,ng)
            enddo 
c$$$            deta = - 2.0 / sqrt(ef)
c$$$            phase = coulphase(deta,lf) * phasel(1,nchf)
c$$$            ovlp = ovlpn(nt)
c$$$            write (41+1,70) lf, llf, sqrt(const*2./3.*om) * 
c$$$     >         ovlp * phase * TR
c$$$            write (41+2,70) lf, llf, sqrt(const*8./3./om) *
c$$$     >         ovlp * phase * TV
c$$$            write (41+3,70) lf, llf, sqrt(const*32./3./om**3) *
c$$$     >         ovlp * phase * TX
 70         format(2i9,4x,1p,2e12.4)
         endif
            

*         write(6,1955)                                          
*         write(6,'(I2,A1,2x,10(3x,E11.4) )')nf,lh(min(10,lf)),
*     >                CLo,CVo,CXo,CL,CV,CX, ovlpn(nt)

         if(nchf.le.4)
     >      write(10+nchf,'(F7.1,3F11.4)') om*ry,CL,CV,CX  
*         write(6,'(//A//)') 

         if(nf.eq.2)then  !n=2 partial cross-section
            CS2l = CS2l+CL
            CS2v = CS2v+Cv
            CS2x = CS2x+Cx
         end if
         
         sigma_l = sigma_l+CL
         sigma_v = sigma_v+Cv
         sigma_x = sigma_x+Cx

         cs(nt,1) = cs(nt,1) + Cl
         cs(nt,2) = cs(nt,2) + Cv
         cs(nt,3) = cs(nt,3) + Cx
         tcs(lf,1) = tcs(lf,1) + Cl
         tcs(lf,2) = tcs(lf,2) + Cv
         tcs(lf,3) = tcs(lf,3) + Cx
         if (ef.lt.0.0) then
            tnbcs(lf,1) = tnbcs(lf,1) + Cl * ovlpn(nt)
            tnbcs(lf,2) = tnbcs(lf,2) + Cv * ovlpn(nt)
            tnbcs(lf,3) = tnbcs(lf,3) + Cx * ovlpn(nt)
         endif 
         if (nt.gt.ntmax) ntmax = nt
         if (lf.gt.ltop) ltop = lf
         
         if(ipseudo.eq.0)then       !Negative energy states only
            proj = ovlpn(nt)
            if (proj.lt.0.0.or.proj.gt.1.0001)
     >         print*, 'CAUTION: OVLPN for CHANNEL:',proj,nt
            sigma1_l = sigma1_l+CL * proj
            sigma1_v = sigma1_v+Cv * proj
            sigma1_x = sigma1_x+Cx * proj
         end if
 113     continue
      end do

      asymcs = 0.0
      do l = 0, ltop
         chs(l) = ' + '
         if (l.eq.0) chs(l) = ' = '
      enddo 
!      do ngauge = 1, 3
      do ngauge = 1, 6
         tcst = 0.0
         tnbcst = 0.0
         write(41+ngauge,'("state XSEC(Mb)       D(Lp-La=-1)       ",
     >      "    D(Lp-La=+1)     overlap beta   en(Ry)")')
         do nt = 1, ntmax
            call getchnl(chan(nt),n,la,nc)
            div = abs(D(nt,1,ngauge))**2+abs(D(nt,-1,ngauge))**2
            if (div.gt.0.0) then
               beta = ((la + 2) * abs(D(nt,1,ngauge))**2 +
     >            (la - 1) * abs(D(nt,-1,ngauge))**2 +
     >            6.0 * sqrt(float(la*(la+1))) *
     >            real(conjg(D(nt,1,ngauge))*D(nt,-1,ngauge)))/
     >            ((2*la+1)*div)
            else
               beta = 0.0
            endif 
            tcst = tcst + cs(nt,ngauge)
            if (enchan(nt) .lt.0.0)
     >         tnbcst = tnbcst + cs(nt,ngauge) * ovlpn(nt)
            if (ngauge.ge.4) cs(nt,ngauge) = div
            write(41+ngauge,'(a3,1p,5e11.3,0p,
     >         2f7.3,f10.6)') chan(nt),cs(nt,ngauge),D(nt,-1,ngauge),
     >         D(nt,1,ngauge), ovlpn(nt),beta,enchan(nt)
         enddo 
         write(41+ngauge,'(79a)') ('-',i=1,79)
         ticst = max(tcst - tnbcst,0.0)
         write(41+ngauge,'(f10.4,"eV on",a3," TNBCS: ",
     >      1p,e9.3,9(a3,e9.3))') om * ry,
     >      chan(1), tnbcst, (chs(l),tnbcs(l,ngauge),l=0,ltop)
         write(41+ngauge,'(f10.4,"eV on",a3,"  TICS: ",
     >      1p,e9.3,9(a3,e9.3))') om * ry,
     >      chan(1), max(0.0,tcst - tnbcst),
     >      (chs(l), max(0.0,tcs(l,ngauge)-tnbcs(l,ngauge)),l=0,ltop)
         write(41+ngauge,'(f10.4,"eV on",a3,"   TCS: ",
     >      1p,e9.3,9(a3,e9.3))') om * ry,
     >      chan(1), tcst, (chs(l), tcs(l,ngauge),l=0,ltop)
      enddo 
      
*      write(6,'(A,10(F9.4,A))')'n=2',CS2l,'&',CS2v,'&',CS2x,'&'
*      print*, ' '
*      print*, ' '

      ein = gk(1,1)**2*ry
      ddl = (sigma_l-sigma1_l)/sigma1_l*100.  
      ddv = (sigma_v-sigma1_v)/sigma1_v*100.  
      ddx = (sigma_x-sigma1_x)/sigma1_x*100.  
      
      write(20,'(F7.2,3F11.4)') om*ry,sigma_l, sigma_v, sigma_x   
      write(21,'(F7.2,3F11.4)') om*ry,sigma1_l,sigma1_v, sigma1_x
      write(22,'(''#ph energy  double L single     double V single'',
     >   ''     double A single   el energy'')') 
      write(22,'(F10.4,1x,1p,6E10.3,0p,f10.4)') om*ry, sigma_l-sigma1_l,
     >   sigma1_l,sigma_v-sigma1_v,sigma1_v,sigma_x-sigma1_x,sigma1_x,
     >   ein
      write(22,'(F10.4,3F11.4,2x,a20)') om*ry, ddl, ddv,ddx,file
      write(23,'(F8.2,3F11.4,2x,a20)') om*ry, ddl, ddv,ddx,file
 30   format(
     >   /' k final',9x,'ReG',11x,'ImG',10x,'ReT',11x,'ImT'/)
 35       format(/'Initial',8x,'Length',19x,'Velocity',
     >   17x,'Accelration'/,'channel',3x,
     >   3('Re            Im          ')/)     
 40   format(10(E11.4,2x))                                             
 41   format(I3,2x,6(E11.4,2x))                                        
 42   format(/I3,2x,6(E11.4,2x))                                       
 43   format(' Sum ',6(E11.4,2x))                                      
 100  format(10(3x,E11.4) )                                            
 1000 FORMAT(//15X,'INTEGRATION OF ',                                  
     >   'THE DIPOLE ME WITH HALF ON-SHELL T-MATRIX'/,15x,56('-')/)    
 1954 format(/11x,'FINAL CHANNEL ',I2/,                                
     >        11x,'Ion energy',F7.3/,                                  
     >        11x,'Bound electron',I2,A1/,                            
     >        11x,'Free  electron k', A1/,                             
     >        11x,'Total k     ', I2)                                  
 1955 format(/14x,'1st Cross Section',             
     >       26x,'TM Cross Section',/             
     >       2(6x,'Length',7x,'Velocty',6x,'Accelration')/,87('-'))
      close(42)
      close(43)
      close(44)
      close(22)
      close(11)
      close(12)
      close(13)
      close(14)
      close(20)
      close(21)
 
      return
      end


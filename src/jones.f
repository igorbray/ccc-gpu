      subroutine JONES

C  Pluvinage-type helium atom ground state               
C  Jones and Madison PRL 91, 073201 (2003)                               

      include 'par.f'
      include 'paratom.f'

      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
     
     
      common /meshrr/ meshr,rmesh(maxr,3)
      common /hyl/    phi1s(maxr),maxphi
      common /hyl3/   psi1s(maxr),en
      COMMON /hyl1/   z1, an, a1, a2, a3, a4, a5, a6, a7, a8, a9, 
     :        a10, a11, a12, a13,a14, a15, a16, a17, a18, a19
      COMMON /hyl2/   b6,b7,b8  ! For hydrogen-11 only
      common/smallr/ formcut,regcut,expcut,fast,match
      logical fast,match

      complex*16 xx,eta,zlmin,cfc(9),cgc(9),cfcp(9),cgcp(9),sig(9),ci
      complex*16 zero, one

      if(ntype.ne.0)Stop "This is not a Hylleraas ground state" 

C Initialization.

      a1 = 0.
      a2 = 0.
      a3 = 0.
      a4 = 0.
      a5 = 0.
      a6 = 0.
      a7 = 0.
      a8 = 0.
      a9 = 0.
      a10= 0.
      a11= 0.
      a12= 0.
      a13= 0.
      a14= 0.
      a15= 0.
      a16= 0.
      a17= 0.
      a18= 0.
      a19= 0.
      b6 = 0.
      b7 = 0.
      b8 = 0.

c Pluvinage ground state

      Z    =   2.               !              -Z(r1+r2)   -ikr12    
      ak   =   0.41             ! F(r r ) = N e           e
      aN   =   0.60337          !    1 2                 
      Etot =   2.903717         !           x (1+a1 r  + ... )   
      q    =   1.0              !                    12

      con  = Z**1.5*2
      con2 = con**2
      con4 = con**4

c Coulomb function Phi_0(1/2k,kr)expansion

!      q   =  0.7
!      aN  =  0.653 
      a1  =  q/2                                         !c1 u
      a5  =   (q**2/2 - ak**2)/6                         !c2 u^2
      a8  = q*(q**2   - ak**2*8)/144                     !c3 u^3
      a13 =   (q**4 - ak**2*q**2*20 + ak**4*24) /2880    !c4 u^4
      a14 = q*(q**4 - ak**2*q**2*40 + ak**4*184)/86400   !c5 u^5
      b6  = (q**6 - ak**2*q**4*70 + ak**4*q**2*784 - ak**6*720)/3628800 
      
c Exact Coulomb function

      eps = 1.E-8                                    
      pi =acos(-1.0)    
      zero = (0.0,0.0)
      one  = (1.0,0.0)
      ci   = (0.0,1.0)
      mode = 4
      kfn = 0
      ln = 0
      zlmin = cmplx(ln)
      nl = 1

      eta = one/ak/2
      
      Co  = pi*2*eta/(exp(pi*2*eta)-1.)
      Co  = sqrt(Co)

      
      do i = 1, meshr
         rt = RMESH(i,1)
         phi1s(i)  =  rt * con * exp(-Z*rt) * sqrt(an)
         if(phi1s(i).gt.expcut)  maxphi=i
         write(1001,'(3E13.4)') rt,  phi1s(i)
         rho = ak*rt
         XX = one*rho
         call coulcc(xx,eta,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode,
     >      kfn,ifai)
         psi1s(i)  =  phi1s(i) * real(cfc(1))/Co/rho
         write(1001,'(3E13.4)') rt,  phi1s(i), psi1s(i)


      end do
      
      en = 1
      maxphi = maxphi/en
      
      print*, 'Cut-off', expcut
      print*, 'Maxphi', maxphi
      print*, 'Rmax', RMESH(maxphi,1)


      write(6,101) nznuc,  Etot
 101  FORMAT (////12x,' JONES-MADISON GROUND STATE   ',
     :           /12x,' --------------------------   '/,
     :           /,   ' Nucleus charge           ',I3,
     :           /,   ' Ground state energy, au  ',F10.6,
     :           /,   ' Parameters Z, N, a1...a14 '//)

      write(6,'(F11.6)') z,an,a1,a5,a8,a13,a14,b6

C  Radial integrals

      call RnOVER(0,phi1s,phi1s,P0) 
      call RnOVER(1,phi1s,phi1s,P1) ! f n -ar        n!
      call RnOVER(2,phi1s,phi1s,P2) ! |r e    dr = ------
      call RnOVER(3,phi1s,phi1s,P3) ! j            a**(n+1) 
      call RnOVER(4,phi1s,phi1s,P4) 
      call RnOVER(6,phi1s,phi1s,P6) 

      write(6,111) '<r^0>=', P0, 'Err=', P0-AN*2  /(2*Z)**3*con2
      write(6,111) '<r^1>=', P1, 'Err=', P1-AN*6  /(2*Z)**4*con2
      write(6,111) '<r^2>=', P2, 'Err=', P2-AN*24 /(2*Z)**5*con2
      write(6,111) '<r^3>=', P3, 'Err=', P3-AN*120/(2*Z)**6*con2
      write(6,111) '<r^4>=', P4, 'Err=', P4-AN*720/(2*Z)**7*con2
      write(6,111)
      
 111  format(A,F11.6,3x,A,2x,E12.5)
      
C  Coulomb integrals

      call MATRnm (V020, phi1s, phi1s, 0, 2, 0)  
      call MATRnm (V111, phi1s, phi1s, 1, 1, 1)  
      call MATRnm (V040, phi1s, phi1s, 0, 4, 0)  
      call MATRnm (V022, phi1s, phi1s, 0, 2, 2)  
      call MATRnm (V222, phi1s, phi1s, 2, 2, 2)  
      call MATRnm (V060, phi1s, phi1s, 0, 6, 0)  
      call MATRnm (V042, phi1s, phi1s, 0, 4, 2)  
      call MATRnm (V133, phi1s, phi1s, 1, 3, 3)  
      call MATRnm (V333, phi1s, phi1s, 3, 3, 3)  

      write(6,111)'<F0 r^2 r^0>=',V020,'Err=',V020-AN**2*21  /2/(2*Z)**7
     >   * con4 
      write(6,111)'<F1 r^1 r^1>=',V111,'Err=',V111-AN**2*21  /4/(2*Z)**7
     >   * con4 
      write(6,111)'<F0 r^4 r^0>=',V040,'Err=',V040-AN**2*1845/8/(2*z)**9
     >   * con4 
      write(6,111)'<F0 r^2 r^2>=',V022,'Err=',V022-AN**2*837 /8/(2*z)**9
     >   * con4 
                                   
      F = P0**2 
     >  +    a1 *2 * (V020*2-V111*2/3)                !U^1
     >  +    a1**2 *  P0*P2*2                         !U^2 
     >  +    a5 *2 *  P0*P2*2                         !U^2 
     >  + a1*a5 *2 * (V040*2+V022*2-V222*4/5)         !U^3 
     >  +    a5**2 * (P0*P4*2+P2**2*10/3)             !U^4
     >  +    a8 *2 * (V040*2+V022*2-V222*4/5)         !U^3
     >  + a1*a8 *2 * (P0*P4*2+P2**2*10/3)             !U^4
     >  + a5*a8 *2 * (V060*2+V042*10-V133*2-V333*6/7) !U^5
     >  +    a8**2 * (P0*P6*2+P2*P4*14)               !U^6 

      write(6,'(A,F9.4)') 'Norma ', F

c Corrections for large distances
      
      call RmOVER(0,psi1s,psi1s,Q0) 
      F = F + P0*Q0*2 

      write(6,'(A,F9.4)') 'Norma ', F

      
C Asymptotics r2->0 r1->oo r12->r2

      
      do i = 1, meshr
         rt = RMESH(i,1)
         tmp =   exp(-z*rt) * (1 
     >      + a1*rt + a5*rt**2
     >      + a8*rt**3 + a13*rt**4
     >      + a14*rt**5)
         write(66,'(3E13.4)') rt, tmp, exp(-z*rt)
      end do

      

  20  FORMAT (2  F17.7)

      END      


      subroutine jDIPOLE(w,w2,l,l1,r1,v1,x1,r,v,x) 
                                                       
      include 'par.f'                                            
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
      common /schf/ psi1s(maxr)                    !w2-continuum
      common /hyl/  phi1s(maxr),maxphi
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr), f(maxr), f1(maxr)
      
C Get rid of the Simpson's weights for continuum WF
C But don't corrupt them for later use. That's why w2, not w1

!      do i=1,maxr
!         if(rmesh(i,3) .ne. 0d0) w1(i) = w2(i)/rmesh(i,3)
!      end do
      do i=1,meshr
         w1(i) = w2(i)/rmesh(i,3)
      end do

      rmatr=0d0                                                  
      dmatr=0d0
      xmatr=0d0
      rmatr1=0d0
      dmatr1=0d0
      xmatr1=0d0

C SP-channel                                  !     \        
                                              ! ns __\__ f1=p  
      Isp=0                                   ! ns __.__ f =s
      Ips=0                                   
      if(l.eq.0 .and. l1.eq.1) Isp=1          !SP or PS-combination
      if(l1.eq.0 .and. l.eq.1) Ips=1          !Always f=s, f1=p  
      if(Isp.eq.1)then                        
         call EQUAL(w,f)                      !f -discrete
         call EQUAL(w1,f1)                    !f1-continuous
      end if
      if(Ips.eq.1)then
         call EQUAL(w1,f)                     !f -continuous
         call EQUAL(w,f1)                     !f1-discrete
      end if
      if(Isp+Ips.eq.0)goto 111

      call jlENGTH  (f,f1, 0,1, rmatr)           
      call jvELOCITY(f,f1, 0,1, dmatr)           
      call jaCCEL   (f,f1, 0,1, xmatr)           

      goto 115
 111  continue


C PD-channel                                  !     \         
                                              ! np __\__ f1=d
      Ipd=0                                   ! np __.__ f =p
      Idp=0
      if(l.eq.1 .and. l1.eq.2) Ipd=1          !PD or DP-combination 
      if(l1.eq.1 .and. l.eq.2) Idp=1          !Always f=p, f1=d     
      if(Ipd.eq.1)then                        !This part has not ben tested
         call EQUAL(w,f)
         call EQUAL(w1,f1)
      end if
      if(Idp.eq.1)then
         call EQUAL(w1,f)
         call EQUAL(w,f1)
      end if
      if(Ipd+Idp.eq.0)goto 112

      call jlENGTH  (f,f1, 1,2, rmatr)           
      call jvELOCITY(f,f1, 1,2, dmatr)           
      call jaCCEL   (f,f1, 1,2, xmatr)           
      goto 115
 112  continue


C DF-channel                                  !     \        
                                              ! nd __\__ f1=f
      Idf=0                                   ! nd __.__ f =f
      Ifd=0
      if(l.eq.2 .and. l1.eq.3) Idf=1          !DF or FD-combination 
      if(l1.eq.2 .and. l.eq.3) Ifd=1          !Always f=d, f1=f     
      if(Idf.eq.1)then                                               
         call EQUAL(w,f)                      !f -discrete   
         call EQUAL(w1,f1)                    !f1-continuous 
      end if                                                 
      if(Ifd.eq.1)then                                       
         call EQUAL(w1,f)                     !f -continuous 
         call EQUAL(w,f1)                     !f1-discrete   
      end if
      if(Idf+Ifd.eq.0)goto 113

      call jlENGTH  (f,f1, 2,3, rmatr)           
      call jvELOCITY(f,f1, 2,3, dmatr)           
      call jaCCEL   (f,f1, 2,3, xmatr)           
      goto 115
 113  continue

C FG-channel                                  !     \        
                                              ! nf __\__ f1=g
      Ifg=0                                   ! nf __.__ f =f
      Igf=0
      if(l.eq.3 .and. l1.eq.4) Ifg=1          !FG or GF-combination 
      if(l1.eq.3 .and. l.eq.4) Igf=1          !Always f=f, f1=g     
      if(Ifg.eq.1)then                        
         call EQUAL(w,f)                      !f -discrete   
         call EQUAL(w1,f1)                    !f1-continuous 
      end if                                                 
      if(Igf.eq.1)then                                       
         call EQUAL(w1,f)                     !f -continuous 
         call EQUAL(w,f1)                     !f1-discrete   
      end if
      if(Ifg+Igf.eq.0)goto 114

      call jlENGTH  (f,f1, 3,4, rmatr)           
      call jvELOCITY(f,f1, 3,4, dmatr)           
      call xACCEL   (f,f1, 3,4, xmatr)           
      goto 115
 114  continue

C GH-channel                                  !     \        
                                              ! ng __\__ f1=h
      Igh=0                                   ! ng __.__ f =g
      Ihg=0
      if(l.eq.4 .and. l1.eq.5) Igh=1          !GH or GF-combination 
      if(l1.eq.4 .and. l.eq.5) Ihg=1          !Always f=g, f1=h     
      if(Igh.eq.1)then                        
         call EQUAL(w,f)                      !f -discrete   
         call EQUAL(w1,f1)                    !f1-continuous 
      end if                                                 
      if(Ihg.eq.1)then                                       
         call EQUAL(w1,f)                     !f -continuous 
         call EQUAL(w,f1)                     !f1-discrete   
      end if
      if(Igh+Ihg.eq.0)goto 115

      call jlENGTH  (f,f1, 4,5, rmatr)           
      call xVELOCITY(f,f1, 4,5, dmatr)           
      call xACCEL   (f,f1, 4,5, xmatr)           
 115  continue

      r =rmatr
      v =dmatr
      x =xmatr

      RETURN
      END

      subroutine jLENGTH(wg,wf,Lg,Lf,result)

C  Calculates dipole ME in length form     
C  from 10-term Hylleraas GS wavefunction for He     
C   r                               
C  D =<gf |r1 V| F > +  < gf |r2 V| F >
C                 L                  L 
C  Here F  is L-pole component of the Hylleraas      
C        L                                        

      include 'par.f'
      DIMENSION  wg(maxr), wf(maxr)
      common /hyl/    phi1s(maxr),maxphi
      COMMON /hyl1/   z1, an, a1, a2, a3, a4, a5, a6, a7, a8, a9, 
     :        a10, a11, a12, a13,a14, a15, a16, a17, a18, a19
      COMMON /hyl2/   b6,b7,b8  ! For hydrogen-11 only
      common /hyl3/   psi1s(maxr),en
      if(iabs(Lf-Lg).ne.1) STOP 'WRONG Lf,Lg' 
      Lm=Lg                                   
      if(Lf.gt.Lg) Lm=Lf                       
      const=float((-1)**(Lm))*sqrt(float(Lm)) 

c Large distance correction for Jones and Madison

      call RmOVER(0, psi1s, wg, Gm0) 
      call RmOVER(1, psi1s, wg, Gm1) 
      call RmOVER(0, psi1s, wf, Fm0) 
      call RmOVER(1, psi1s, wf, Fm1) 
      
      
C Initialization

      call RnOVER(0, phi1s, wg, G0) 
      call RnOVER(1, phi1s, wg, G1) 
      call RnOVER(2, phi1s, wg, G2) 
      call RnOVER(3, phi1s, wg, G3) 
      call RnOVER(4, phi1s, wg, G4) 
      call RnOVER(5, phi1s, wg, G5) 
      call RnOVER(6, phi1s, wg, G6) 
      call RnOVER(7, phi1s, wg, G7) 
      call RnOVER(0, phi1s, wf, F0) 
      call RnOVER(1, phi1s, wf, F1) 
      call RnOVER(2, phi1s, wf, F2) 
      call RnOVER(3, phi1s, wf, F3) 
      call RnOVER(4, phi1s, wf, F4) 
      call RnOVER(5, phi1s, wf, F5) 
      call RnOVER(6, phi1s, wf, F6) 
      call RnOVER(7, phi1s, wf, F7) 

      call MATRnm (p003, wg, wf, 0, 0, 3)
      call MATRnm (p004, wg, wf, 0, 0, 4)
      call MATRnm (p005, wg, wf, 0, 0, 5)
      call MATRnm (p006, wg, wf, 0, 0, 6)
      call MATRnm (p012, wg, wf, 0, 1, 2)
      call MATRnm (p013, wg, wf, 0, 1, 3)
      call MATRnm (p014, wg, wf, 0, 1, 4)
      call MATRnm (p015, wg, wf, 0, 1, 5)
      call MATRnm (p021, wg, wf, 0, 2, 1)
      call MATRnm (p022, wg, wf, 0, 2, 2)
      call MATRnm (p023, wg, wf, 0, 2, 3)
      call MATRnm (p024, wg, wf, 0, 2, 4)
      call MATRnm (p030, wg, wf, 0, 3, 0)
      call MATRnm (p031, wg, wf, 0, 3, 1)
      call MATRnm (p032, wg, wf, 0, 3, 2)
      call MATRnm (p040, wg, wf, 0, 4, 0)
      call MATRnm (p041, wg, wf, 0, 4, 1)
      call MATRnm (p042, wg, wf, 0, 4, 2)
      call MATRnm (p050, wg, wf, 0, 5, 0)
      call MATRnm (p051, wg, wf, 0, 5, 1)
      call MATRnm (p060, wg, wf, 0, 6, 0)
      call MATRnm (p112, wg, wf, 1, 1, 2)
      call MATRnm (p113, wg, wf, 1, 1, 3)
      call MATRnm (p114, wg, wf, 1, 1, 4)
      call MATRnm (p115, wg, wf, 1, 1, 5)
      call MATRnm (p121, wg, wf, 1, 2, 1)
      call MATRnm (p122, wg, wf, 1, 2, 2)
      call MATRnm (p123, wg, wf, 1, 2, 3)
      call MATRnm (p124, wg, wf, 1, 2, 4)
      call MATRnm (p131, wg, wf, 1, 3, 1)
      call MATRnm (p132, wg, wf, 1, 3, 2)
      call MATRnm (p133, wg, wf, 1, 3, 3)
      call MATRnm (p141, wg, wf, 1, 4, 1)
      call MATRnm (p142, wg, wf, 1, 4, 2)
      call MATRnm (p151, wg, wf, 1, 5, 1)
      call MATRnm (p223, wg, wf, 2, 2, 3)
      call MATRnm (p232, wg, wf, 2, 3, 2)
      call MATRnm (p332, wg, wf, 3, 3, 2)
      call MATRnm (p323, wg, wf, 3, 2, 3)

      call MATRnm (p070, wg, wf, 0, 7, 0)
      call MATRnm (p061, wg, wf, 0, 6, 1)
      call MATRnm (p052, wg, wf, 0, 5, 2)
      call MATRnm (p043, wg, wf, 0, 4, 3)
      call MATRnm (p034, wg, wf, 0, 3, 4)
      call MATRnm (p025, wg, wf, 0, 2, 5)
      call MATRnm (p016, wg, wf, 0, 1, 6)
      call MATRnm (p007, wg, wf, 0, 0, 7)

      call MATRnm (p152, wg, wf, 1, 5, 2)
      call MATRnm (p143, wg, wf, 1, 4, 3)
      call MATRnm (p134, wg, wf, 1, 3, 4)
      call MATRnm (p125, wg, wf, 1, 2, 5)

      call MATRnm (p252, wg, wf, 2, 5, 2)
      call MATRnm (p243, wg, wf, 2, 4, 3)
      call MATRnm (p234, wg, wf, 2, 3, 4)
      call MATRnm (p225, wg, wf, 2, 2, 5)

      call MATRnm (p352, wg, wf, 3, 5, 2)
      call MATRnm (p343, wg, wf, 3, 4, 3)
      call MATRnm (p334, wg, wf, 3, 3, 4)
      call MATRnm (p325, wg, wf, 3, 2, 5)

      call MATRnm (p452, wg, wf, 4, 5, 2)
      call MATRnm (p443, wg, wf, 4, 4, 3)
      call MATRnm (p434, wg, wf, 4, 3, 4)
      call MATRnm (p425, wg, wf, 4, 2, 5)

      call MATRnm (p543, wg, wf, 5, 4, 3)
      call MATRnm (p534, wg, wf, 5, 3, 4)

C First term D1: L=Lf->Lg transition

      L = Lf
      IF(L .eq. 0) then
         D1  = G1*F0 + a1 * (p030  + p012  - p121*2/3)
     :               + a2 * (G3*F0 + G1*F2 - G2*F1*2)
     :               + a4 * (G3*F0 + G1*F2 + G2*F1*2)
     :               + a5 * (G3*F0 + G1*F2)          
     :               + a3 * (G2*F0 + G1*F1)          
     : + a6 * (p040 + p022 + p031 + p013 - p131*2/3 - p122*2/3)
     : + a7 * (p050 + p032*2 + p014 - p041*2 - p023*2
     :                       - p141*2/3 - p123*2/3 + p132*4/3)
     : + a8 * (p050 + p014 + p032*2 - p232*4/5)
     : + a9 * (G5*F0 + G3*F2*2 + G1*F4 - G4*F1*2 - G2*F3*2)
     : + a10* (F3*G1 - F2*G2 - F1*G3 + F0*G4)
     : + a11* (F3*G1 + F2*G2*3 + F1*G3*3 + F0*G4)
     : + a12* (F6*G1 - F5*G2*2 + F4*G3*13/3 - F3*G4*20/3 + F2*G5*13/3
     :                - F1*G6*2 + F0*G7)
     : + a13* (F4*G1 + F2*G3*10/3 + F0*G5)
     : +  b6* (F6*G1 + F4*G3*7 + F2*G5*7 + F0*G7)
     : + a14 * (p070 + 5*p052 - 2*p143 - 6*p343/7 + 5*p034 + p016)
     : + a15 * (p070-2*p061+3*p052-4*p043+3*p034-2*p025+p016
     :       -4*p252/5 + 8*p243/5 - 4*p234/5)
     : + a16 * (G5*F0 - 2*G3*F2 + G1*F4)
     : + a17 * (G5*F0 + 4*G4*F1 + 6*G3*F2 + 4*G2*F3 + G1*F4)
     : + a18 * (p060 - p051 - p024 + p015
     :      -2*p151/3 + 2*p142/3 + 2*p133/3 - 2*p124/3)
     : + a19 * (G5*F0 - 4*G4*F1 + 6*G3*F2 - 4*G2*F3 + G1*F4)

         D1 = D1 + Gm1*F0 + G1*Fm0
      ELSE
         call MATRnm (s21, wg, wf, L+1, 2, 1)
         call MATRnm (s22, wg, wf, L+1, 2, 2)
         call MATRnm (s23, wg, wf, L+1, 2, 3)
         call MATRnm (s31, wg, wf, L+1, 3, 1)
         call MATRnm (s32, wg, wf, L+1, 3, 2)
         call MATRnm (s41, wg, wf, L+1, 4, 1)
         call MATRnm (t21, wg, wf, L-1, 2, 1)
         call MATRnm (t22, wg, wf, L-1, 2, 2)
         call MATRnm (t23, wg, wf, L-1, 2, 3)
         call MATRnm (t31, wg, wf, L-1, 3, 1)
         call MATRnm (t32, wg, wf, L-1, 3, 2)
         call MATRnm (t41, wg, wf, L-1, 4, 1)

         call MATRnm (P, wg, wf, L+2, 3, 2)
         call MATRnm (Q, wg, wf, L  , 3, 2)
         call MATRnm (R, wg, wf, L-2, 3, 2)
         
         D1 = a1 * (s21/(2*L+3) - t21/(2*L-1)) 
     :      + a6 *((s31+s22)/(2*L+3) - (t31+t22)/(2*L-1)) 
     :      + a7 *((s41+s23-s32*2)/(2*L+3) - (t41+t23-t32*2)/(2*L-1)) 
         B8 = a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))
         call MATRnm (P, wg, wf, L+3, 4, 3)
         call MATRnm (Q, wg, wf, L+1, 4, 3)
         call MATRnm (R, wg, wf, L-1, 4, 3)
         call MATRnm (S, wg, wf, L-3, 4, 3)

         B14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         call MATRnm (P52, wg, wf, L+2, 5, 2)
         call MATRnm (Q52, wg, wf, L  , 5, 2)
         call MATRnm (R52, wg, wf, L-2, 5, 2)
         call MATRnm (P43, wg, wf, L+2, 4, 3)
         call MATRnm (Q43, wg, wf, L  , 4, 3)
         call MATRnm (R43, wg, wf, L-2, 4, 3)
         call MATRnm (P34, wg, wf, L+2, 3, 4)
         call MATRnm (Q34, wg, wf, L  , 3, 4)
         call MATRnm (R34, wg, wf, L-2, 3, 4)

         B15= a15 * 3*( (P52-2*P43+P34)/( 3+2*L)/( 5+2*L)
     :              - 2*(Q52-2*Q43+Q34)/(-1+2*L)/( 3+2*L)
     :                 +(R52-2*R43+R34)/(-3+2*L)/(-1+2*L) )
     
         call MATRnm (S51, wg, wf, L+1, 5, 1)
         call MATRnm (S42, wg, wf, L+1, 4, 2)
         call MATRnm (S33, wg, wf, L+1, 3, 3)
         call MATRnm (S24, wg, wf, L+1, 2, 4)
         call MATRnm (T51, wg, wf, L-1, 5, 1)
         call MATRnm (T42, wg, wf, L-1, 4, 2)
         call MATRnm (T33, wg, wf, L-1, 3, 3)
         call MATRnm (T24, wg, wf, L-1, 2, 4)

         B18 = a18 *( (S51 - S42 - S33 + S24)/(2*L+3)
     :              - (T51 - T42 - T33 + T24)/(2*L-1) )

         IF(L .eq. 1) THEN
            D1  = D1 - a5 * G2*F1*2
     :               - a9 *(G4*F1 + G2*F3 - G3*F2*2)*2 
     :          + a12*4*(-F5*G2 + F4*G3*2 - F3*G4*2 + F2*G5*2 - F1*G6)
     :          + a13*4*(-F3*G2 - F1*G4)                              
     :          +  b6*6/5*(-F5*G2*5 - F3*G4*14 - F1*G6*5)                              
            B8 = -a8 *(p041 + p023 - p132*3/5 - p332/35)*3 
            B14= a14 *(-5*p061 - 10*p043 - 5*p025
     :               +  9*p152/5 + 3*p352/35 + 9*p134/5 + 3*p334/35
     :               -  4*p043/5 + 64*p243/35 - 4*p443/105)
            B15= a15 *(-3*p061 + 6*p052 - 6*p043 + 6*p034 - 3*p025
     :                + 9*p152/5 + 3*p352/35 - 18*p143/5
     :                - 6*p343/35 + 9*p134/5 + 3*p334/35)
         END IF

         IF(L .eq. 2) THEN
            D1  = D1
     :          + a12* 8/3 * (F4*G3 - 2*F3*G4 + F2*G5)
     :          + a13* 8/3 *  F2*G3
     :          +  b6* 8*(F4*G3 + F2*G5)
            B14 = a14 *(5*p052 - 2*p252/7 + p452/21
     :                + 5*p034 - 2*p234/7 + p434/21
     :              -  18*p143/7 - 2*p543/77)
         END IF

         IF(L .eq. 3) THEN
            D1  = D1
     :          -  b6* 16/5* F3*G4
         END IF
         D1  = D1 + B8 + B14 + B15 + B18
      END IF
         
C Second term D2: L=Lg->Lf transition

      L = Lg
      IF(L .eq. 0) then

         D2  = F1*G0 + a1 * (p021  + p003  - p112*2/3)
     :               + a2 * (G0*F3 + G2*F1 - G1*F2*2)   
     :               + a4 * (G0*F3 + G2*F1 + G1*F2*2)  
     :               + a5 * (G0*F3 + G2*F1)            
     :               + a3 * (G0*F2 + G1*F1)            
     : + a6 * (p031 + p013 + p022 + p004 - p122*2/3 - p113*2/3)
     : + a7 * (p041 + p023*2 + p005 - p032*2 - p014*2
     :                        - p132*2/3 - p114*2/3 + p123*4/3)
     : + a8 * (p041 + p005 + p023*2 - p223*4/5)
     : + a9 * (G4*F1 + G2*F3*2 + G0*F5 - G3*F2*2 - G1*F4*2)
     : + a10* (F4*G0 - F3*G1 - F2*G2 + F1*G3)                       
     : + a11* (F4*G0 + F3*G1*3 + F2*G2*3 + F1*G3)                   
     : + a12* (F7*G0 - F6*G1*2 + F5*G2*13/3 - F4*G3*20/3 + F3*G4*13/3
     :               - F2*G5*2 + F1*G6)
     : + a13* (F5*G0 + F3*G2*10/3 + F1*G4)
     : +  b6* (F7*G0 + F5*G2*7 + F3*G4*7 + F1*G6)
     : + a14 * (p061 + 5*p043 - 2*p134 - 6*p334/7 + 5*p025 + p007)
     : + a15 * (p061-2*p052+3*p043-4*p034+3*p025-2*p016+p007
     :       -4*p243/5 + 8*p234/5 - 4*p225/5)
     : + a16 * (G4*F1 - 2*G2*F3 + G0*F5)
     : + a17 * (G4*F1 + 4*G3*F2 + 6*G2*F3 + 4*G1*F4 + G0*F5)
     : + a18 * (p051 - p042 - p015 + p006
     :      -2*p142/3 + 2*p133/3 + 2*p124/3 - 2*p115/3)
     : + a19 * (G4*F1 - 4*G3*F2 + 6*G2*F3 - 4*G1*F4 + G0*F5)

         D2  = D2 + Fm1*G0 + F1*Gm0
      ELSE 

         call MATRnm (s12, wg, wf, L+1, 1, 2)
         call MATRnm (s13, wg, wf, L+1, 1, 3)
         call MATRnm (s14, wg, wf, L+1, 1, 4)
         call MATRnm (s22, wg, wf, L+1, 2, 2)
         call MATRnm (s23, wg, wf, L+1, 2, 3)
         call MATRnm (s32, wg, wf, L+1, 3, 2)
         call MATRnm (t12, wg, wf, L-1, 1, 2)
         call MATRnm (t13, wg, wf, L-1, 1, 3)
         call MATRnm (t14, wg, wf, L-1, 1, 4)
         call MATRnm (t22, wg, wf, L-1, 2, 2)
         call MATRnm (t23, wg, wf, L-1, 2, 3)
         call MATRnm (t32, wg, wf, L-1, 3, 2)

         call MATRnm (P, wg, wf, L+2, 2, 3)
         call MATRnm (Q, wg, wf, L  , 2, 3)
         call MATRnm (R, wg, wf, L-2, 2, 3)
         
         D2 = a1 * (s12/(2*L+3) - t12/(2*L-1)) 
     :      + a6 *((s22+s13)/(2*L+3) - (t22+t13)/(2*L-1)) 
     :      + a7 *((s32+s14-s23*2)/(2*L+3) - (t32+t14-t23*2)/(2*L-1)) 
         C8 = a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))
         call MATRnm (P, wg, wf, L+3, 3, 4)
         call MATRnm (Q, wg, wf, L+1, 3, 4)
         call MATRnm (R, wg, wf, L-1, 3, 4)
         call MATRnm (S, wg, wf, L-3, 3, 4)

         C14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         call MATRnm (P43, wg, wf, L+2, 4, 3)                
         call MATRnm (Q43, wg, wf, L  , 4, 3)                
         call MATRnm (R43, wg, wf, L-2, 4, 3)                
         call MATRnm (P34, wg, wf, L+2, 3, 4)                
         call MATRnm (Q34, wg, wf, L  , 3, 4)                
         call MATRnm (R34, wg, wf, L-2, 3, 4)                
         call MATRnm (P25, wg, wf, L+2, 2, 5)                
         call MATRnm (Q25, wg, wf, L  , 2, 5)                
         call MATRnm (R25, wg, wf, L-2, 2, 5)              
  
         C15= a15 * 3*( (P43-2*P34+P25)/( 3+2*L)/( 5+2*L) 
     :               -2*(Q43-2*Q34+Q25)/(-1+2*L)/( 3+2*L)
     :                 +(R43-2*R34+R25)/(-3+2*L)/(-1+2*L) )
                                                           
         call MATRnm (S42, wg, wf, L+1, 4, 2)
         call MATRnm (S33, wg, wf, L+1, 3, 3)
         call MATRnm (S24, wg, wf, L+1, 2, 4)
         call MATRnm (S15, wg, wf, L+1, 1, 5)
         call MATRnm (T42, wg, wf, L-1, 4, 2)
         call MATRnm (T33, wg, wf, L-1, 3, 3)
         call MATRnm (T24, wg, wf, L-1, 2, 4)
         call MATRnm (T15, wg, wf, L-1, 1, 5)

         C18 = a18 *( (S42 - S33 - S24 + S15)/(2*L+3)
     :              - (T42 - T33 - T24 + T15)/(2*L-1) )

         IF(L .eq. 1) THEN
            D2  = D2 - a5 * G1*F2*2
     :               - a9 *(G3*F2 + G1*F4 - G2*F3*2)*2 
     :        + a12*4 * (-F6*G1 + F5*G2*2 - F4*G3*2 + F3*G4*2 - F2*G5)
     :        + a13*4 * (-F4*G1 - F2*G3)                              
     :        +  b6*6/5 * (-F6*G1*5 - F4*G3*14 - F2*G5)                              
            C8 =     - a8 *(p032 + p014 - p123*3/5 - p323/35)*3 
           C14 =      a14 *(-5*p052 - 10*p034 - 5*p016
     :                   +  9*p143/5 + 3*p343/35 + 9*p125/5 + 3*p325/35
     :                   -  4*p034/5 + 64*p234/35 - 4*p434/105)
           C15=       a15 *(-3*p052 + 6*p043 - 6*p034 + 6*p025 - 3*p016
     :                     + 9*p143/5 + 3*p343/35 - 18*p134/5
     :                     - 6*p334/35 + 9*p125/5 + 3*p325/35)
         END IF

         IF(L .eq. 2) THEN
            D2  = D2
     :          + a12* 8/3 * (F5*G2 - F4*G3*2 + F3*G4)
     :          + a13* 8/3 *  F3*G2 
     :          +  b6* 8  * ( F5*G2 + F3*G4 )
            C14 = a14 *(5*p043 - 2*p243/7 + p443/21 
     :               +  5*p025 - 2*p225/7 + p425/21
     :               - 18*p134/7 - 2*p534/77)
         END IF

         IF(L .eq. 3) THEN
            D1  = D1
     :          -  b6* 16/5* F4*G3
         END IF

         D2  = D2 + C8 + C14 + C15 + C18
      END IF
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const


      RETURN
      END

      subroutine jVELOCITY(wg,wf,lg,lf,result)     
                                             
C  Calculates dipole ME in velocity form     
C  from Hylleraas GS wavefunction for He     
C                                    
C   V   dg                  df        
C  D =< -- f |V| F > +  < g -- |V| F >
C       dr        L         dr      L 
C  Here F  is L-pole component of the Hylleraas      
C        L             
      include 'par.f'
      DIMENSION  YD(maxr)
      DIMENSION  wg(MAXR), wf(MAXR), wd(MAXR)
      common /meshrr/ meshr,rmesh(maxr,3)
      common /hyl/    phi1s(maxr),maxphi
      COMMON /hyl1/   z1, an, a1, a2, a3, a4, a5, a6, a7, a8, a9, 
     :        a10, a11, a12, a13,a14, a15, a16, a17, a18, a19
      COMMON /hyl2/   b6,b7,b8  ! For hydrogen-11 only
      common /hyl3/   psi1s(maxr),en
      if(iabs(Lf-Lg).ne.1) STOP 'WRONG Lf,Lg'
      Lm=Lg                                   !Lmax=max(Lg,Lf) 
      if(Lf.gt.Lg) Lm=Lf                      !        Lmax  ____ 
      const=float((-1)**(Lm))*sqrt(float(Lm)) !const=(-1)   VLmax

C First term D1: L=Lf->Lg transition

      L = Lf
      if(Lg.gt.Lf) IL= 1                      !Lg=Lf+1
      if(Lg.lt.Lf) IL=-1                      !Lg=Lf-1
      call DERIVATIVE(wg,wd)
      DO i=1,meshr
         YD(i) = wd(i) + wg(i)/rmesh(i,1)*il*Lm
      end do

c Large distance correction for Jones and Madison

      call RmOVER(0, psi1s, YD, Gm0) 
      call RmOVER(1, psi1s, YD, Gm1) 
      call RmOVER(0, psi1s, wf, Fm0) 
      call RmOVER(1, psi1s, wf, Fm1) 
      
      call RnOVER(0, phi1s, YD, G0) 
      call RnOVER(1, phi1s, YD, G1) 
      call RnOVER(2, phi1s, YD, G2) 
      call RnOVER(3, phi1s, YD, G3) 
      call RnOVER(4, phi1s, YD, G4) 
      call RnOVER(5, phi1s, YD, G5) 
      call RnOVER(6, phi1s, YD, G6) 
      call RnOVER(0, phi1s, wf, F0) 
      call RnOVER(1, phi1s, wf, F1) 
      call RnOVER(2, phi1s, wf, F2) 
      call RnOVER(3, phi1s, wf, F3) 
      call RnOVER(4, phi1s, wf, F4)
      call RnOVER(5, phi1s, wf, F5)
      call RnOVER(6, phi1s, wf, F6)
      
      call MATRnm (p002, YD, wf, 0, 0, 2)
      call MATRnm (p003, YD, wf, 0, 0, 3)
      call MATRnm (p004, YD, wf, 0, 0, 4)
      call MATRnm (p005, YD, wf, 0, 0, 5)
      call MATRnm (p012, YD, wf, 0, 1, 2)
      call MATRnm (p013, YD, wf, 0, 1, 3)
      call MATRnm (p014, YD, wf, 0, 1, 4)
      call MATRnm (p020, YD, wf, 0, 2, 0)
      call MATRnm (p021, YD, wf, 0, 2, 1)
      call MATRnm (p022, YD, wf, 0, 2, 2)
      call MATRnm (p030, YD, wf, 0, 3, 0)
      call MATRnm (p031, YD, wf, 0, 3, 1)
      call MATRnm (p040, YD, wf, 0, 4, 0)
      call MATRnm (p041, YD, wf, 0, 4, 1)
      call MATRnm (p050, YD, wf, 0, 5, 0)
      call MATRnm (p111, YD, wf, 1, 1, 1)
      call MATRnm (p112, YD, wf, 1, 1, 2)
      call MATRnm (p113, YD, wf, 1, 1, 3)
      call MATRnm (p121, YD, wf, 1, 2, 1)
      call MATRnm (p122, YD, wf, 1, 2, 2)
      call MATRnm (p131, YD, wf, 1, 3, 1)
      call MATRnm (p222, YD, wf, 2, 2, 2)
      call MATRnm (p322, YD, wf, 3, 2, 2)
      
      call MATRnm (p006, YD, wf, 0, 0, 6)
      call MATRnm (p015, YD, wf, 0, 1, 5)
      call MATRnm (p024, YD, wf, 0, 2, 4)
      call MATRnm (p033, YD, wf, 0, 3, 3)
      call MATRnm (p042, YD, wf, 0, 4, 2)
      call MATRnm (p051, YD, wf, 0, 5, 1)
      call MATRnm (p060, YD, wf, 0, 6, 0)
      
      call MATRnm (p114, YD, wf, 1, 1, 4)
      call MATRnm (p132, YD, wf, 1, 3, 2)
      call MATRnm (p133, YD, wf, 1, 3, 3)
      call MATRnm (p123, YD, wf, 1, 2, 3)
      call MATRnm (p124, YD, wf, 1, 2, 4)
      call MATRnm (p141, YD, wf, 1, 4, 1)
      call MATRnm (p142, YD, wf, 1, 4, 2)
      
      call MATRnm (p233, YD, wf, 2, 3, 3)
      call MATRnm (p224, YD, wf, 2, 2, 4)
      call MATRnm (p242, YD, wf, 2, 4, 2)
      
      call MATRnm (p333, YD, wf, 3, 3, 3)
      call MATRnm (p324, YD, wf, 3, 2, 4)
      call MATRnm (p342, YD, wf, 3, 4, 2)
      
      call MATRnm (p433, YD, wf, 4, 3, 3)
      call MATRnm (p424, YD, wf, 4, 2, 4)
      call MATRnm (p442, YD, wf, 4, 4, 2)

      call MATRnm (p533, YD, wf, 5, 3, 3)
         
      IF(L .eq. 0) then

         D1  = G0*F0 + a1 * (p020 + p002 - p111*2/3)
     :               + a2 * (G2*F0 + G0*F2 - G1*F1*2)  
     :               + a4 * (G2*F0 + G0*F2 + G1*F1*2)  
     :               + a5 * (G2*F0 + G0*F2)            
     :               + a3 * (G1*F0 + G0*F1)            
     : + a6 * (p030 + p012 + p021 + p003 - p121*2/3 - p112*2/3)
     : + a7 * (p040 + p022*2 + p004 - p031*2 - p013*2
     :                       - p131*2/3 - p113*2/3 + p122*4/3)
     : + a8 * (p040 + p004 + p022*2 - p222*4/5)
     : + a9 * (G4*F0 + G2*F2*2 + G0*F4 - G3*F1*2 - G1*F3*2)
     : + a10* (F3*G0 - F2*G1 - F1*G2 + F0*G3)
     : + a11* (F3*G0 + F2*G1*3 + F1*G2*3 +F0*G3)
     : + a12* (F6*G0 - F5*G1*2 + F4*G2*13/3 - F3*G3*20/3 + F2*G4*13/3
     :               - F1*G5*2 + F0*G6)
     : + a13* (F4*G0 + F2*G2*10/3 + F0*G4)
     : + a14 * (p060 + 5*p042 + 5*p024 + p006
     :       -   2*p133 - 6*p333/7)
     : + a15 * (p060-2*p051+3*p042-4*p033+3*p024-2*p015+p006
     :       -4*p242/5 + 8*p233/5 - 4*p224/5)
     : + a16 * (G4*F0 - 2*G2*F2 + G0*F4)
     : + a17 * (G4*F0 + 4*G3*F1 + 6*G2*F2 + 4*G1*F3 + G0*F4)
     : + a18 * (p050 - p041 - p014 + p005
     :      - 2*p141/3 + 2*p132/3 + 2*p123/3 - 2*p114/3)
     : + a19 * (G4*F0 - 4*G3*F1 + 6*G2*F2 - 4*G1*F3 + G0*F4)

         D1  = D1 + Gm0*F0 + G0*Fm0

      ELSE
         call MATRnm (s11, YD, wf, L+1, 1, 1)
         call MATRnm (s12, YD, wf, L+1, 1, 2)
         call MATRnm (s21, YD, wf, L+1, 2, 1)
         call MATRnm (s13, YD, wf, L+1, 1, 3)
         call MATRnm (s31, YD, wf, L+1, 3, 1)
         call MATRnm (s22, YD, wf, L+1, 2, 2)
         call MATRnm (t11, YD, wf, L-1, 1, 1)
         call MATRnm (t12, YD, wf, L-1, 1, 2)
         call MATRnm (t21, YD, wf, L-1, 2, 1)
         call MATRnm (t13, YD, wf, L-1, 1, 3)
         call MATRnm (t31, YD, wf, L-1, 3, 1)
         call MATRnm (t22, YD, wf, L-1, 2, 2)

         call MATRnm (P, YD, wf, L+2, 2, 2)
         call MATRnm (Q, YD, wf, L  , 2, 2)
         call MATRnm (R, YD, wf, L-2, 2, 2)
         
         D1 = a1 * (s11/(2*L+3) - t11/(2*L-1)) 
     :      + a6 *((s21+s12)/(2*L+3) - (t21+t12)/(2*L-1)) 
     :      + a7 *((s31+s13-s22*2)/(2*L+3) - (t31+t13-t22*2)/(2*L-1)) 
         B8= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2
     :                                 + R/(2*L-1)/(2*L-3))

         call MATRnm (P, YD, wf, L+3, 3, 3)
         call MATRnm (Q, YD, wf, L+1, 3, 3)
         call MATRnm (R, YD, wf, L-1, 3, 3)
         call MATRnm (S, YD, wf, L-3, 3, 3)
                                    
         B14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         call MATRnm (P42, YD, wf, L+2, 4, 2)
         call MATRnm (Q42, YD, wf, L  , 4, 2)
         call MATRnm (R42, YD, wf, L-2, 4, 2)
         call MATRnm (P33, YD, wf, L+2, 3, 3)
         call MATRnm (Q33, YD, wf, L  , 3, 3)
         call MATRnm (R33, YD, wf, L-2, 3, 3)
         call MATRnm (P24, YD, wf, L+2, 2, 4)
         call MATRnm (Q24, YD, wf, L  , 2, 4)
         call MATRnm (R24, YD, wf, L-2, 2, 4)

         B15= a15 * 3*( (P42-2*P33+P24)/( 3+2*L)/( 5+2*L)
     :              - 2*(Q42-2*Q33+Q24)/(-1+2*L)/( 3+2*L)
     :                 +(R42-2*R33+R24)/(-3+2*L)/(-1+2*L) )

         call MATRnm (S41, YD, wf, L+1, 4, 1)
         call MATRnm (S32, YD, wf, L+1, 3, 2)
         call MATRnm (S23, YD, wf, L+1, 2, 3)
         call MATRnm (S14, YD, wf, L+1, 1, 4)
         call MATRnm (T41, YD, wf, L-1, 4, 1)
         call MATRnm (T32, YD, wf, L-1, 3, 2)
         call MATRnm (T23, YD, wf, L-1, 2, 3)
         call MATRnm (T14, YD, wf, L-1, 1, 4)

         B18 = a18 *( (S41 - S32 - S23 + S14)/(2*L+3)
     :              - (T41 - T32 - T23 + T14)/(2*L-1) )

         IF(L .eq. 1) THEN
            D1 =D1 - a5 * G1*F1*2
     :             - a9 *(G3*F1 + G1*F3 - G2*F2*2)*2 
     :             + a12 *4*(-F5*G1 + F4*G2*2 - F3*G3*2
     :                        + F2*G4*2 - F1*G5)
     :             + a13 *4*(-F3*G1 - F1*G3)
            B8=    - a8  *(p031 + p013 - p122*3/5 - p322/35)*3 
            B14=     a14 *(-5*p051 - 10*p033 - 5*p015
     :                 +  9*p142/5 + 3*p342/35 +  9*p124/5 + 3*p324/35
     :                   -  4*p033/5 + 64*p233/35 - 4*p433/105)
            B15=     a15 *(-3*p051 + 6*p042 - 6*p033 + 6*p024 - 3*p015
     :                    + 9*p142/5 + 3*p342/35 - 18*p133/5 
     :                    - 6*p333/35 + 9*p124/5 + 3*p324/35 )
         END IF

         IF(L .eq. 2) THEN
            D1  = D1
     :         +     a12 * 8/3 * (F4*G2 - F3*G3*2 + F2*G4)
     :         +     a13 * 8/3 *  F2*G2
            B14  =   a14 *(5*p042 - 2*p242/7 + p442/21
     :                    +5*p024 - 2*p224/7 + p424/21
     :                  - 18*p133/7 - 2*p533/77)
        END IF
        D1  = D1 + B8 + B14 + B15 + B18
      END IF
         
C Second term D2: L=Lg->Lf transition

      L = Lg
      if(Lf.gt.Lg) IL= 1                      !Lf=Lg+1
      if(Lf.lt.Lg) IL=-1                      !Lf=Lg-1
      call DERIVATIVE(wf,wd)
      DO i=1,meshr
         YD(i) = wd(i) + wf(i)/rmesh(i,1)*il*Lm
      end do
      
c Large distance correction for Jones and Madison

      call RmOVER(0, psi1s, wg, Gm0) 
      call RmOVER(1, psi1s, wg, Gm1) 
      call RmOVER(0, psi1s, YD, Fm0) 
      call RmOVER(1, psi1s, YD, Fm1) 

      call RnOVER(0, phi1s, wg, G0) 
      call RnOVER(1, phi1s, wg, G1) 
      call RnOVER(2, phi1s, wg, G2) 
      call RnOVER(3, phi1s, wg, G3) 
      call RnOVER(4, phi1s, wg, G4) 
      call RnOVER(5, phi1s, wg, G5) 
      call RnOVER(6, phi1s, wg, G6) 
      call RnOVER(0, phi1s, YD, F0) 
      call RnOVER(1, phi1s, YD, F1) 
      call RnOVER(2, phi1s, YD, F2) 
      call RnOVER(3, phi1s, YD, F3) 
      call RnOVER(4, phi1s, YD, F4)
      call RnOVER(5, phi1s, YD, F5)
      call RnOVER(6, phi1s, YD, F6)
 
      call MATRnm (p002, wg, YD, 0, 0, 2)
      call MATRnm (p003, wg, YD, 0, 0, 3)
      call MATRnm (p004, wg, YD, 0, 0, 4)
      call MATRnm (p005, wg, YD, 0, 0, 5)
      call MATRnm (p012, wg, YD, 0, 1, 2)
      call MATRnm (p013, wg, YD, 0, 1, 3)
      call MATRnm (p014, wg, YD, 0, 1, 4)
      call MATRnm (p020, wg, YD, 0, 2, 0)
      call MATRnm (p021, wg, YD, 0, 2, 1)
      call MATRnm (p022, wg, YD, 0, 2, 2)
      call MATRnm (p030, wg, YD, 0, 3, 0)
      call MATRnm (p031, wg, YD, 0, 3, 1)
      call MATRnm (p040, wg, YD, 0, 4, 0)
      call MATRnm (p041, wg, YD, 0, 4, 1)
      call MATRnm (p050, wg, YD, 0, 5, 0)
      call MATRnm (p111, wg, YD, 1, 1, 1)
      call MATRnm (p112, wg, YD, 1, 1, 2)
      call MATRnm (p113, wg, YD, 1, 1, 3)
      call MATRnm (p121, wg, YD, 1, 2, 1)
      call MATRnm (p122, wg, YD, 1, 2, 2)
      call MATRnm (p131, wg, YD, 1, 3, 1)
      call MATRnm (p222, wg, YD, 2, 2, 2)
      call MATRnm (p322, wg, YD, 3, 2, 2)
      
      call MATRnm (p006, wg, YD, 0, 0, 6)
      call MATRnm (p015, wg, YD, 0, 1, 5)
      call MATRnm (p024, wg, YD, 0, 2, 4)
      call MATRnm (p033, wg, YD, 0, 3, 3)
      call MATRnm (p042, wg, YD, 0, 4, 2)
      call MATRnm (p051, wg, YD, 0, 5, 1)
      call MATRnm (p060, wg, YD, 0, 6, 0)
      
      call MATRnm (p114, wg, YD, 1, 1, 4)
      call MATRnm (p132, wg, YD, 1, 3, 2)
      call MATRnm (p133, wg, YD, 1, 3, 3)
      call MATRnm (p123, wg, YD, 1, 2, 3)
      call MATRnm (p124, wg, YD, 1, 2, 4)
      call MATRnm (p141, wg, YD, 1, 4, 1)
      call MATRnm (p142, wg, YD, 1, 4, 2)
      
      call MATRnm (p233, wg, YD, 2, 3, 3)
      call MATRnm (p224, wg, YD, 2, 2, 4)
      call MATRnm (p242, wg, YD, 2, 4, 2)
      
      call MATRnm (p333, wg, YD, 3, 3, 3)
      call MATRnm (p324, wg, YD, 3, 2, 4)
      call MATRnm (p342, wg, YD, 3, 4, 2)
      
      call MATRnm (p433, wg, YD, 4, 3, 3)
      call MATRnm (p424, wg, YD, 4, 2, 4)
      call MATRnm (p442, wg, YD, 4, 4, 2)
      
      call MATRnm (p533, wg, YD, 5, 3, 3)

      IF(L .eq. 0) then

         D2  = G0*F0 + a1 * (p020 + p002 - p111*2/3)
     :               + a2 * (G2*F0 + G0*F2 - G1*F1*2)  
     :               + a4 * (G2*F0 + G0*F2 + G1*F1*2)  
     :               + a5 * (G2*F0 + G0*F2)            
     :               + a3 * (G1*F0 + G0*F1)            
     : + a6 * (p030 + p012 + p021 + p003 - p121*2/3 - p112*2/3)
     : + a7 * (p040 + p022*2 + p004 - p031*2 - p013*2
     :                        - p131*2/3 - p113*2/3 + p122*4/3)
     : + a8 * (p040 + p004 + p022*2 - p222*4/5)
     : + a9 * (G4*F0 + G2*F2*2 + G0*F4 - G3*F1*2 - G1*F3*2)
     : + a10* (F3*G0 - F2*G1 - F1*G2 + F0*G3)
     : + a11* (F3*G0 + F2*G1*3 + F1*G2*3 +F0*G3)
     : + a12* (F6*G0 - F5*G1*2 + F4*G2*13/3 - F3*G3*20/3 + F2*G4*13/3
     :               - F1*G5*2 + F0*G6)
     : + a13* (F4*G0 + F2*G2*10/3 + F0*G4)
     : + a14 * (p060 + 5*p042 + 5*p024 + p006
     :       -   2*p133 - 6*p333/7)
     : + a15 * (p060-2*p051+3*p042-4*p033+3*p024-2*p015+p006
     :       -4*p242/5 + 8*p233/5 - 4*p224/5)
     : + a16 * (G4*F0 - 2*G2*F2 + G0*F4)
     : + a17 * (G4*F0 + 4*G3*F1 + 6*G2*F2 + 4*G1*F3 + G0*F4)
     : + a18 * (p050 - p041 - p014 + p005
     :      - 2*p141/3 + 2*p132/3 + 2*p123/3 - 2*p114/3)
     : + a19 * (G4*F0 - 4*G3*F1 + 6*G2*F2 - 4*G1*F3 + G0*F4)

         D2  = D2 + Gm0*F0 + G0*Fm0
      ELSE
         call MATRnm (s11, wg, YD, L+1, 1, 1)
         call MATRnm (s12, wg, YD, L+1, 1, 2)
         call MATRnm (s21, wg, YD, L+1, 2, 1)
         call MATRnm (s13, wg, YD, L+1, 1, 3)
         call MATRnm (s31, wg, YD, L+1, 3, 1)
         call MATRnm (s22, wg, YD, L+1, 2, 2)
         call MATRnm (t11, wg, YD, L-1, 1, 1)
         call MATRnm (t12, wg, YD, L-1, 1, 2)
         call MATRnm (t21, wg, YD, L-1, 2, 1)
         call MATRnm (t13, wg, YD, L-1, 1, 3)
         call MATRnm (t31, wg, YD, L-1, 3, 1)
         call MATRnm (t22, wg, YD, L-1, 2, 2)

         call MATRnm (P, wg, YD, L+2, 2, 2)
         call MATRnm (Q, wg, YD, L  , 2, 2)
         call MATRnm (R, wg, YD, L-2, 2, 2)
         
         D2 = a1 * (s11/(2*L+3) - t11/(2*L-1)) 
     :      + a6 *((s21+s12)/(2*L+3) - (t21+t12)/(2*L-1)) 
     :      + a7 *((s31+s13-s22*2)/(2*L+3) - (t31+t13-t22*2)/(2*L-1)) 
         C8= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2
     :                                 + R/(2*L-1)/(2*L-3))
         call MATRnm (P, wg, YD, L+3, 3, 3)
         call MATRnm (Q, wg, YD, L+1, 3, 3)
         call MATRnm (R, wg, YD, L-1, 3, 3)
         call MATRnm (S, wg, YD, L-3, 3, 3)
                                    
         C14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         call MATRnm (P42, wg, YD, L+2, 4, 2)
         call MATRnm (Q42, wg, YD, L  , 4, 2)
         call MATRnm (R42, wg, YD, L-2, 4, 2)
         call MATRnm (P33, wg, YD, L+2, 3, 3)
         call MATRnm (Q33, wg, YD, L  , 3, 3)
         call MATRnm (R33, wg, YD, L-2, 3, 3)
         call MATRnm (P24, wg, YD, L+2, 2, 4)
         call MATRnm (Q24, wg, YD, L  , 2, 4)
         call MATRnm (R24, wg, YD, L-2, 2, 4)
         
         C15= a15 * 3*( (P42-2*P33+P24)/( 3+2*L)/( 5+2*L)
     :              - 2*(Q42-2*Q33+Q24)/(-1+2*L)/( 3+2*L)
     :                 +(R42-2*R33+R24)/(-3+2*L)/(-1+2*L) )

         call MATRnm (S41, wg, YD, L+1, 4, 1)
         call MATRnm (S32, wg, YD, L+1, 3, 2)
         call MATRnm (S23, wg, YD, L+1, 2, 3)
         call MATRnm (S14, wg, YD, L+1, 1, 4)
         call MATRnm (T41, wg, YD, L-1, 4, 1)
         call MATRnm (T32, wg, YD, L-1, 3, 2)
         call MATRnm (T23, wg, YD, L-1, 2, 3)
         call MATRnm (T14, wg, YD, L-1, 1, 4)

         C18 = a18 *( (S41 - S32 - S23 + S14)/(2*L+3)
     :              - (T41 - T32 - T23 + T14)/(2*L-1) )

         IF(L .eq. 1) THEN
            D2  = D2 - a5 * G1*F1*2
     :               - a9 *(G3*F1 + G1*F3 - G2*F2*2)*2 
     :  +  a12 * 4 * (-F5*G1 + F4*G2*2 - F3*G3*2 + F2*G4*2 - F1*G5)
     :  +  a13 * 4 * (-F3*G1 - F1*G3)
           C8=     - a8 *(p031 + p013 - p122*3/5 - p322/35)*3 
           C14=     a14 *(-5*p051 - 10*p033 - 5*p015
     :                  +  9*p142/5 + 3*p342/35 +  9*p124/5 + 3*p324/35
     :                  -  4*p033/5 + 64*p233/35 - 4*p433/105)
           C15=     a15 *(-3*p051 + 6*p042 - 6*p033 + 6*p024 - 3*p015
     :                   + 9*p142/5 + 3*p342/35 - 18*p133/5 
     :                   - 6*p333/35 + 9*p124/5 + 3*p324/35 )
         END IF
         IF(L .eq. 2) THEN
            D2  = D2
     :         +     a12 * 8/3 * (F4*G2 - F3*G3*2 + F2*G4)
     :         +     a13 * 8/3 *  F2*G2
            C14  =   a14 *(5*p042 - 2*p242/7 + p442/21
     :                    +5*p024 - 2*p224/7 + p424/21
     :                  - 18*p133/7 - 2*p533/77)
         END IF
         D2  = D2 + C8 + C14 + C15 + C18
      END IF
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const

      RETURN
      END

      subroutine jACCEL(wg,wf,Lg,Lf,result)

C  Calculates dipole ME in acceleration form     
C  from 10-term Hylleraas GS wavefunction for He     
C   .                                 
C   V           -2                   -2   
C  D = 2 <gf |r1  V| F > + 2 < gf |r2  V| F >
C                     L                    L 
C  Here F  is L-pole component of the Hylleraas      
C        L             
      include 'par.f'
      DIMENSION  wg(maxr), wf(maxr)
      common /hyl/    phi1s(maxr),maxphi
      COMMON /hyl1/   z1, an, a1, a2, a3, a4, a5, a6, a7, a8, a9, 
     :        a10, a11, a12, a13,a14, a15, a16, a17, a18, a19
      COMMON /hyl2/   b6,b7,b8  ! For hydrogen-11 only
      common /hyl3/   psi1s(maxr),en
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym

      if(iabs(Lf-Lg).ne.1) STOP 'WRONG Lf,Lg' 
      Lm=Lg                                   
      if(Lf.gt.Lg) Lm=Lf                       
      const=float((-1)**(Lm))*sqrt(float(Lm)) 

      const = const*nznuc

c Large distance correction for Jones and Madison

      call RmOVER(-2,psi1s, wg, Gm_2) 
      call RmOVER( 0,psi1s, wg, Gm0) 
      call RmOVER(-2,psi1s, wf, Fm_2) 
      call RmOVER( 0,psi1s, wf, Fm0) 

C Initialization

      call RnOVER(-2,phi1s, wg, G_2) 
      call RnOVER(-1,phi1s, wg, G_1) 
      call RnOVER( 0,phi1s, wg, G0) 
      call RnOVER( 1,phi1s, wg, G1) 
      call RnOVER( 2,phi1s, wg, G2) 
      call RnOVER( 3,phi1s, wg, G3) 
      call RnOVER( 4,phi1s, wg, G4) 
      call RnOVER( 5,phi1s, wg, G5) 
      call RnOVER( 6,phi1s, wg, G6) 
      call RnOVER(-2,phi1s, wf, F_2) 
      call RnOVER(-1,phi1s, wf, F_1) 
      call RnOVER( 0,phi1s, wf, F0) 
      call RnOVER( 1,phi1s, wf, F1) 
      call RnOVER( 2,phi1s, wf, F2) 
      call RnOVER( 3,phi1s, wf, F3) 
      call RnOVER( 4,phi1s, wf, F4) 
      call RnOVER( 5,phi1s, wf, F5) 
      call RnOVER( 6,phi1s, wf, F6) 

      call MATRnm (p0_22,wg, wf, 0,-2, 2)
      call MATRnm (p0_23,wg, wf, 0,-2, 3)
      call MATRnm (p0_24,wg, wf, 0,-2, 4)
      call MATRnm (p0_25,wg, wf, 0,-2, 5)
      call MATRnm (p0_12,wg, wf, 0,-1, 2)
      call MATRnm (p0_13,wg, wf, 0,-1, 3)
      call MATRnm (p0_14,wg, wf, 0,-1, 4)
      call MATRnm (p000, wg, wf, 0, 0, 0)
      call MATRnm (p001, wg, wf, 0, 0, 1)
      call MATRnm (p002, wg, wf, 0, 0, 2)
      call MATRnm (p003, wg, wf, 0, 0, 3)
      call MATRnm (p010, wg, wf, 0, 1, 0)
      call MATRnm (p011, wg, wf, 0, 1, 1)
      call MATRnm (p012, wg, wf, 0, 1, 2)
      call MATRnm (p02_2,wg, wf, 0, 2,-2)
      call MATRnm (p02_1,wg, wf, 0, 2,-1)
      call MATRnm (p020, wg, wf, 0, 2, 0)
      call MATRnm (p021, wg, wf, 0, 2, 1)
      call MATRnm (p030, wg, wf, 0, 3, 0)
      call MATRnm (p03_2,wg, wf, 0, 3,-2)
      call MATRnm (p03_1,wg, wf, 0, 3,-1)
      call MATRnm (p030, wg, wf, 0, 3, 0)
      call MATRnm (p04_1,wg, wf, 0, 4,-1)
      call MATRnm (p04_2,wg, wf, 0, 4,-2)
      call MATRnm (p05_2,wg, wf, 0, 5,-2)
      call MATRnm (p1_11,wg, wf, 1,-1, 1)
      call MATRnm (p1_12,wg, wf, 1,-1, 2)
      call MATRnm (p1_13,wg, wf, 1,-1, 3)
      call MATRnm (p1_14,wg, wf, 1,-1, 4)
      call MATRnm (p101, wg, wf, 1, 0, 1)
      call MATRnm (p102, wg, wf, 1, 0, 2)
      call MATRnm (p103, wg, wf, 1, 0, 3)
      call MATRnm (p11_1,wg, wf, 1, 1,-1)
      call MATRnm (p110, wg, wf, 1, 1, 0)
      call MATRnm (p111, wg, wf, 1, 1, 1)
      call MATRnm (p112, wg, wf, 1, 1, 2)
      call MATRnm (p12_1,wg, wf, 1, 2,-1)
      call MATRnm (p120, wg, wf, 1, 2, 0)
      call MATRnm (p121, wg, wf, 1, 2, 1)
      call MATRnm (p130, wg, wf, 1, 3, 0)
      call MATRnm (p13_1,wg, wf, 1, 3,-1)
      call MATRnm (p14_1,wg, wf, 1, 4,-1)
      call MATRnm (p202, wg, wf, 2, 0, 2)
      call MATRnm (p220, wg, wf, 2, 2, 0)
      call MATRnm (p302, wg, wf, 3, 0, 2)
      call MATRnm (p320, wg, wf, 3, 2, 0)

      call MATRnm (p0_26,wg, wf, 0,-2, 6)
      call MATRnm (p0_15,wg, wf, 0,-1, 5)
      call MATRnm (p004, wg, wf, 0, 0, 4)
      call MATRnm (p013, wg, wf, 0, 1, 3)
      call MATRnm (p022, wg, wf, 0, 2, 2)
      call MATRnm (p031, wg, wf, 0, 3, 1)
      call MATRnm (p040, wg, wf, 0, 4, 0)
      call MATRnm (p05_1,wg, wf, 0, 5,-1)
      call MATRnm (p06_2,wg, wf, 0, 6,-2)

      call MATRnm (p122, wg, wf, 1, 2, 2)
      call MATRnm (p113, wg, wf, 1, 1, 3)
      call MATRnm (p131, wg, wf, 1, 3, 1)
      call MATRnm (p104, wg, wf, 1, 0, 4)
      call MATRnm (p140, wg, wf, 1, 4, 0)

      call MATRnm (p222, wg, wf, 2, 2, 2)
      call MATRnm (p213, wg, wf, 2, 1, 3)
      call MATRnm (p231, wg, wf, 2, 3, 1)
      call MATRnm (p204, wg, wf, 2, 0, 4)
      call MATRnm (p240, wg, wf, 2, 4, 0)

      call MATRnm (p322, wg, wf, 3, 2, 2)
      call MATRnm (p313, wg, wf, 3, 1, 3)
      call MATRnm (p331, wg, wf, 3, 3, 1)
      call MATRnm (p304, wg, wf, 3, 0, 4)
      call MATRnm (p340, wg, wf, 3, 4, 0)

      call MATRnm (p422, wg, wf, 4, 2, 2)
      call MATRnm (p413, wg, wf, 4, 1, 3)
      call MATRnm (p431, wg, wf, 4, 3, 1)
      call MATRnm (p404, wg, wf, 4, 0, 4)
      call MATRnm (p440, wg, wf, 4, 4, 0)

      call MATRnm (p513, wg, wf, 5, 1, 3)
      call MATRnm (p531, wg, wf, 5, 3, 1)

C First term D1: L=Lf->Lg transition

      L = Lf
      IF(L .eq. 0) then

         D1  = G_2*F0 + a1 * (p000  + p0_22  - p1_11*2/3)
     :                + a2 * (G0*F0 + G_2*F2 - G_1*F1*2)
     :                + a4 * (G0*F0 + G_2*F2 + G_1*F1*2)
     :                + a5 * (G0*F0 + G_2*F2)          
     :                + a3 * (G_1*F0 + G_2*F1)          
     : + a6 * (p010 + p0_12 + p001 + p0_23 - p101*2/3 - p1_12*2/3)
     : + a7 * (p020 + p002*2 + p0_24 - p011*2 - p0_13*2
     :                       - p111*2/3 - p1_13*2/3 + p102*4/3)
     : + a8 * (p020 + p0_24 + p002*2 - p202*4/5)
     : + a9 * (G2*F0 + G0*F2*2 + G_2*F4 - G1*F1*2 - G_1*F3*2)
     : + a10* (F3*G_2 - F2*G_1 - F1*G0 + F0*G1)
     : + a11* (F3*G_2 + F2*G_1*3 + F1*G0*3 + F0*G1)
     : + a12* (F6*G_2 - F5*G_1*2 + F4*G0*13/3 - F3*G1*20/3 + F2*G2*13/3
     :                - F1*G3*2 + F0*G4)
     : + a13* (F4*G_2 + F2*G0*10/3 + F0*G2)
     : + a14 * (p040 + 5*p022 + 5*p004 + p0_26
     :       -  2*p113 - 6*p313/7)
     : + a15 * (p040-2*p031+3*p022-4*p013+3*p004-2*p0_15+p0_26
     :       -4*p222/5 + 8*p213/5 - 4*p204/5)
     : + a16 * (G2*F0 - 2*G0*F2 + G_2*F4)
     : + a17 * (G2*F0 + 4*G1*F1 + 6*G0*F2 + 4*G_1*F3+ G_2*F4)
     : + a18 * (p030 - p021 - p0_14 + p0_25
     :      -2*p121/3 + 2*p112/3 + 2*p103/3 - 2*p1_14/3) 
     : + a19 * (G2*F0 - 4*G1*F1 + 6*G0*F2 - 4*G_1*F3 + G_2*F4)

         D1  = D1 + Gm_2*F0 +  G_2*Fm0

      ELSE
         call MATRnm (s_11,wg, wf, L+1,-1, 1)
         call MATRnm (s_12,wg, wf, L+1,-1, 2)
         call MATRnm (s_13,wg, wf, L+1,-1, 3)
         call MATRnm (s01, wg, wf, L+1, 0, 1)
         call MATRnm (s02, wg, wf, L+1, 0, 2)
         call MATRnm (s11, wg, wf, L+1, 1, 1)
         call MATRnm (t_11,wg, wf, L-1,-1, 1)
         call MATRnm (t_12,wg, wf, L-1,-1, 2)
         call MATRnm (t_13,wg, wf, L-1,-1, 3)
         call MATRnm (t01, wg, wf, L-1, 0, 1)
         call MATRnm (t02, wg, wf, L-1, 0, 2)
         call MATRnm (t11, wg, wf, L-1, 1, 1)

         call MATRnm (P, wg, wf, L+2, 0, 2)
         call MATRnm (Q, wg, wf, L  , 0, 2)
         call MATRnm (R, wg, wf, L-2, 0, 2)
         
         D1 = a1 * (s_11/(2*L+3) - t_11/(2*L-1)) 
     : + a6 *((s01+s_12)/(2*L+3) - (t01+t_12)/(2*L-1)) 
     : + a7 *((s11+s_13-s02*2)/(2*L+3) - (t11+t_13-t02*2)/(2*L-1))
         B8 = a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))

         call MATRnm (P, wg, wf, L+3, 1, 3)
         call MATRnm (Q, wg, wf, L+1, 1, 3)
         call MATRnm (R, wg, wf, L-1, 1, 3)
         call MATRnm (S, wg, wf, L-3, 1, 3)

         B14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         call MATRnm (P22, wg, wf, L+2, 2, 2)              
         call MATRnm (Q22, wg, wf, L  , 2, 2)            
         call MATRnm (R22, wg, wf, L-2, 2, 2)            
         call MATRnm (P13, wg, wf, L+2, 1, 3)            
         call MATRnm (Q13, wg, wf, L  , 1, 3)            
         call MATRnm (R13, wg, wf, L-2, 1, 3)            
         call MATRnm (P04, wg, wf, L+2, 0, 4)            
         call MATRnm (Q04, wg, wf, L  , 0, 4)            
         call MATRnm (R04, wg, wf, L-2, 0, 4)            

         B15= a15 * 3*( (P22-2*P13+P04)/( 3+2*L)/( 5+2*L)
     :              - 2*(Q22-2*Q13+Q04)/(-1+2*L)/( 3+2*L)
     :                 +(R22-2*R13+R04)/(-3+2*L)/(-1+2*L) )

         call MATRnm (S 21, wg, wf, L+1, 2, 1)
         call MATRnm (S 12, wg, wf, L+1, 1, 2)
         call MATRnm (S 03, wg, wf, L+1, 0, 3)
         call MATRnm (S_14, wg, wf, L+1,-1, 4)
         call MATRnm (T 21, wg, wf, L-1, 2, 1)
         call MATRnm (T 12, wg, wf, L-1, 1, 2)
         call MATRnm (T 03, wg, wf, L-1, 0, 3)
         call MATRnm (T_14, wg, wf, L-1,-1, 4)

         B18 = a18 *( (S21 - S12 - S03 + S_14)/(2*L+3)
     :              - (T21 - T12 - T03 + T_14)/(2*L-1) )

         IF(L .eq. 1) THEN
            D1  = D1 - a5 * G_1*F1*2
     :               - a9 *(G1*F1 + G_1*F3 - G0*F2*2)*2 
     :          + a12*4*(-F5*G_1 + F4*G0*2 - F3*G1*2 + F2*G2*2 - F1*G3)
     :          + a13*4*(-F3*G_1 - F1*G1)                              
         B8=     - a8 *(p011 + p0_13 - p102*3/5 - p302/35)*3 
         B14=     a14 *(-5*p031 - 10*p013 - 5*p0_15
     :                 + 9*p122/5 + 3*p322/35 + 9*p104/5 + 3*p304/35
     :                 - 4*p013/5 + 64*p213/35 - 4*p413/105)
         B15=     a15 *(-3*p031 + 6*p022 - 6*p013 + 6*p004 - 3*p0_15
     :                 + 9*p122/5 + 3*p322/35 - 18*p113/5
     :                 - 6*p313/35 + 9*p104/5 + 3*p304/35)

         END IF

         IF(L .eq. 2) THEN
            D1  = D1
     :          + a12* 8/3 * (F4*G0 - F3*G1*2 + F2*G2)
     :          + a13* 8/3 *  F2*G0
            B14 = a14 *(5*p022 - 2*p222/7 + p422/21 
     :                + 5*p004 - 2*p204/7 + p404/21
     :               - 18*p113/7 - 2*p513/77)
           END IF
         D1  = D1 + B8 + B14 + B15 + B18
      END IF
         
C Second term D2: L=Lg->Lf transition

      L = Lg
      IF(L .eq. 0) then

         D2  = G0*F_2 + a1 * (p02_2  + p000  - p11_1*2/3)
     :                + a2 * (G0*F0 + G2*F_2 - G1*F_1*2)   
     :                + a4 * (G0*F0 + G2*F_2 + G1*F_1*2)  
     :                + a5 * (G0*F0 + G2*F_2)            
     :                + a3 * (G0*F_1+ G1*F_2)            
     : + a6 * (p03_2 + p010 + p02_1 + p001 - p12_1*2/3 - p110*2/3)
     : + a7 * (p04_2 + p020*2 + p002 - p03_1*2 - p011*2
     :                        - p13_1*2/3 - p111*2/3 + p120*4/3)
     : + a8 * (p04_2 + p002 + p020*2 - p220*4/5)
     : + a9 * (G4*F_2 + G2*F0*2 + G0*F2 - G3*F_1*2 - G1*F1*2)
     : + a10* (F1*G0 - F0*G1 - F_1*G2 + F_2*G3)                       
     : + a11* (F1*G0 + F0*G1*3 + F_1*G2*3 + F_2*G3)                   
     : + a12* (F4*G0 - F3*G1*2 + F2*G2*13/3 - F1*G3*20/3 + F0*G4*13/3
     :               - 2*F_1*G5 + F_2*G6)
     : + a13* (F2*G0 + F0*G2*10/3 + F_2*G4)
     : + a14 * (p06_2 + 5*p040 + 5*p022 + p004
     :       - 2*p131 - 6*p331/7)
     : + a15 * (p06_2-2*p05_1+3*p040-4*p031+3*p022-2*p013+p004
     :       -4*p240/5 + 8*p231/5 - 4*p222/5)
     : + a16 * (G4*F_2 -2*G2*F0 + G0*F2)
     : + a17 * (G4*F_2+ 4*G3*F_1+ 6*G2*F0 + 4*G1*F1 + G0*F2)
     : + a18 * (p05_2 - p04_1 - p012 + p003
     :      -2*p14_1/3 + 2*p121/3 - 2*p112/3 + 2*p130/3)
     : + a19 * (G4*F_2 - 4*G3*F_1+ 6*G2*F0 - 4*G1*F1 + G0*F2)

         D2  = D2 + Gm0*F_2 + G0*Fm_2
      ELSE 

         call MATRnm (s1_1, wg, wf, L+1, 1,-1)
         call MATRnm (s10,  wg, wf, L+1, 1, 0)
         call MATRnm (s11,  wg, wf, L+1, 1, 1)
         call MATRnm (s2_1, wg, wf, L+1, 2,-1)
         call MATRnm (s20,  wg, wf, L+1, 2, 0)
         call MATRnm (s3_1, wg, wf, L+1, 3,-1)
         call MATRnm (t1_1, wg, wf, L-1, 1,-1)
         call MATRnm (t10,  wg, wf, L-1, 1, 0)
         call MATRnm (t11,  wg, wf, L-1, 1, 1)
         call MATRnm (t2_1, wg, wf, L-1, 2,-1)
         call MATRnm (t20,  wg, wf, L-1, 2, 0)
         call MATRnm (t3_1, wg, wf, L-1, 3,-1)

         call MATRnm (P, wg, wf, L+2, 2, 0)
         call MATRnm (Q, wg, wf, L  , 2, 0)
         call MATRnm (R, wg, wf, L-2, 2, 0)
         
         D2 = a1 * (s1_1/(2*L+3) - t1_1/(2*L-1)) 
     :      + a6 *((s2_1+s10)/(2*L+3) - (t2_1+t10)/(2*L-1)) 
     :      + a7 *((s3_1+s11-s20*2)/(2*L+3) - (t3_1+t11-t20*2)/(2*L-1)) 
         C8 = a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))
         call MATRnm (P, wg, wf, L+3, 3, 1)
         call MATRnm (Q, wg, wf, L+1, 3, 1)
         call MATRnm (R, wg, wf, L-1, 3, 1)
         call MATRnm (S, wg, wf, L-3, 3, 1)

         C14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         call MATRnm (P40, wg, wf, L+2, 4, 0)       
         call MATRnm (Q40, wg, wf, L  , 4, 0)       
         call MATRnm (R40, wg, wf, L-2, 4, 0)       
         call MATRnm (P31, wg, wf, L+2, 3, 1)       
         call MATRnm (Q31, wg, wf, L  , 3, 1)       
         call MATRnm (R31, wg, wf, L-2, 3, 1)       
         call MATRnm (P22, wg, wf, L+2, 2, 2)       
         call MATRnm (Q22, wg, wf, L  , 2, 2)       
         call MATRnm (R22, wg, wf, L-2, 2, 2)       

         C15= a15 * 3*( (P40-2*P31+P22)/( 3+2*L)/( 5+2*L) 
     :               -2*(Q40-2*Q31+Q22)/(-1+2*L)/( 3+2*L)
     :                 +(R40-2*R31+R22)/(-3+2*L)/(-1+2*L) )
                                                           
         call MATRnm (S4_1, wg, wf, L+1, 4,-1)
         call MATRnm (S3 0, wg, wf, L+1, 3, 0)
         call MATRnm (S2 1, wg, wf, L+1, 2, 1)
         call MATRnm (S1 2, wg, wf, L+1, 1, 2)
         call MATRnm (T4_1, wg, wf, L-1, 4,-1)
         call MATRnm (T3 0, wg, wf, L-1, 3, 0)
         call MATRnm (T2 1, wg, wf, L-1, 2, 1)
         call MATRnm (T1 2, wg, wf, L-1, 1, 2)

         C18 = a18 *( (S4_1 - S30 - S21 + S12)/(2*L+3)
     :              - (T4_1 - T30 - T21 + T12)/(2*L-1) )

         IF(L .eq. 1) THEN
            D2  = D2 - a5 * G1*F_1*2
     :               - a9 *(G3*F_1 + G1*F1 - G2*F0*2)*2 
     :       + a12*4 * (-F3*G1 + F2*G2*2 - F1*G3*2 + F0*G4*2 - F_1*G5)
     :       + a13*4 * (-F1*G1 - F_1*G3)                              
             C8=     - a8 *(p03_1 + p011 - p120*3/5 - p320/35)*3 
             C14 =     a14*(-5*p05_1 - 10*p031 - 5*p013
     :                 +  9*p140/5 + 3*p340/35 + 9*p122/5 + 3*p322/35
     :                  - 4*p031/5 + 64*p231/35 - 4*p431/105)
             C15 =   a15 *(-3*p05_1+ 6*p040 - 6*p031 + 6*p022 - 3*p013
     :                    + 9*p140/5 + 3*p340/35 - 18*p131/5
     :                    - 6*p331/35 + 9*p122/5 + 3*p322/35)
         END IF

         IF(L .eq. 2) THEN
            D2  = D2
     :          + a12* 8/3 * (F2*G2 - F1*G3*2 + F0*G4)
     :          + a13* 8/3 *  F0*G2
            C14 = a14 *(5*p040 - 2*p240/7 + p440/21 
     :                + 5*p022 - 2*p222/7 + p422/21
     :               - 18*p131/7 - 2*p531/77)
         END IF
         D2  = D2 + C8 + C14 + C15 + C18
      END IF
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const
      RETURN
      END
      

      subroutine RmOVER(N, P, Q, result)

C Radial overlap integral between two bound states
                                                                        
C                   oo                                                   
C                   f  n                                                 
C              I =  | r  dr P(r) * Q(r)                                  
C                   j                                                    
C                   Rmax                                                    
      
      include 'par.f'

      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /hyl/    phi1s(maxr),maxphi
      common /hyl3/   psi1s(maxr),en

      result = 0.d0
      do i = maxphi, maxphi*en !It is assumed that P = phi1s
         rt = rmesh(i,1)
         sw = rmesh(i,3)  ! Simpson's weights added
         result = result + P(i)*Q(i)*rt**n * sw
      end do 

      RETURN
      END
      

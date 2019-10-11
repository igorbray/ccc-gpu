
      subroutine LeSech

C  Pluvinage-type helium atom ground state               
C  LeSech JPB 30, L47 (1997)                               

      include 'par.f'
      include 'paratom.f'

      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
     
C     psi1s   - 1s atomic ground states
C     phi1s   - 1s for 3-parameter Hylleraas
     
      common /meshrr/ meshr,rmesh(maxr,3)
      common /hyl/    psi1s(maxr),maxpsi !To pass to Hylleraas 
      common /plu/    phi1s(maxr),phip1s(maxr),maxphi
      COMMON /plu1/   Z, an, a1,a5,a8,a13,a14
      common/smallr/ formcut,regcut,expcut,fast,match
      logical fast,match
         

      if(ntype.ne.0)Stop "This is not a Hylleraas ground state" 

C Initialization.

      a1 = 0.
      a5 = 0.
      a8 = 0.
      a13= 0.
      a14= 0.

c Pluvinage ground state

      Z    =   2.               !            -Z(r1+r2)                     
      alam =   0.7              ! F(r r ) = Ne {cosh(lam r1)+cosh(lam r2)} 
      aN   =   0.70205          !    1 2                 
      a    =   0.17             !                     -ar12  
      Etot =   2.903717         !           x (1+0.5r e  )   
                                !                    12
c Taylor series

      a1  = 0.5                 !r12
      a5  =-a1 * a              !r12^2
      a8  = a1 * a**2/2         !r12^3
      a13 =-a1 * a**3/6         !r12^4
      a14 = a1 * a**4/24        !r12^5

c Jones and Madison      
c Coulomb function Phi_0(1/2k,kr)expansion

      ak  =  0.61
!      aN   = 0.57
!      alam = 1.12
      q   =  0.82
      a1  =  q/2                                         !c1 u
      a5  =   (q**2/2 - ak**2)/6                         !c2 u^2
      a8  = q*(q**2   - ak**2*8)/144                     !c3 u^3
      a13 =   (q**4 - ak**2*q**2*20 + ak**4*24) /2880    !c4 u^4
      a14 = q*(q**4 - ak**2*q**2*40 + ak**4*184)/86400   !c5 u^5
      b6  = (q**6 - ak**2*q**4*70 + ak**4*q**2*784 - ak**6*720)/3628800 
      
c  Test with 2-term Hylleraas by Green et al PR 112, 1187 (1952)

!      Z   =   1.8497            !              -Z(r1+r2)                    
!      aN  =   1.3900            ! F(r r ) = N e          {1 + a1  R12}      
!      a1  =   0.3658            !    1 2        
!      alam=   0.
!      Etot=   2.903717
      
c Jones and Madison

!      Z    =   2.               !              -Z(r1+r2)   -ikr12    
!      ak   =   0.41             ! F(r r ) = N e           e
!      aN   =   0.60337          !    1 2                 
!      Etot =   2.903717         !           x (1+a1 r  + ... )   
                                !                    12
!      alam = 0.0
!      aN   = aN * sqrt(1.62)
      
c Coulomb function Phi_0(1/2k,kr)expansion

!      a1  = 0.5                                  !c1 u
!      a5  = (0.5 - ak**2)/6                      !c2 u^2
!      a8  = (1.0 - ak**2*8)/144                  !c3 u^3
!      a13 = (1.0 - ak**2*20 + ak**4*24) /2880    !c4 u^4
!      a14 = (1.0 - ak**2*40 + ak**4*184)/86400   !c5 u^5
 


      eps = 1.E-8                                    
      pi =acos(-1.0)    

      do i = 1, meshr
         rt = RMESH(i,1)
         psi1s(i)  =  rt * exp(-Z*rt) * sqrt(an)
         phi1s(i)  =  rt * exp(-Z*rt) * sqrt(an)
!         phip1s(i) =  rt * exp(-Z*rt) * sqrt(an) * cosh(alam*rt)
         phip1s(i) =  rt * 0.5*(exp(-(Z+alam)*rt)+exp(-(Z-alam)*rt)) 
     >                                * sqrt(an) 
         if(phi1s(i).gt.expcut)  maxphi=i
         write(1001,'(3E13.4)') rt,  phi1s(i),  phip1s(i)
         write(2002,'(3E13.4)') rt, exp(-(Z-alam)*rt),  exp(-Z*rt)
      end do
      maxpsi = maxphi
      
      print*, 'Cut-off', expcut
      print*, 'Maxphi', maxphi


      write(6,101) nznuc,  Etot
 101  FORMAT (////12x,' PLUVINAGE GROUND STATE   ',
     :           /12x,' ----------------------   '/,
     :           /,   ' Nucleus charge           ',I3,
     :           /,   ' Ground state energy, au  ',F10.6,
     :           /,   ' Parameters Z, N, a1...a14 '//)

      write(6,'(F11.6)') z,an,a1,a5,a8,a13,a14

C  Radial integrals

      
 111  format(A,2E13.6,3x,A,2x,E12.5)
      

c 1st Order integrals
      
      call VMATRnm (V020, phip1s, phi1s, 0, 2, 0)  
      call VMATRnm (V002, phip1s, phi1s, 0, 0, 2)  
      call VMATRnm (V111, phip1s, phi1s, 1, 1, 1)  

      call MATRnm (Q020, phi1s, phi1s, 0, 2, 0)  
      call MATRnm (Q002, phi1s, phi1s, 0, 0, 2)  
      call MATRnm (Q111, phi1s, phi1s, 1, 1, 1)  


      write(6,111) "<1s1s'|r1^2 F0|Phi 0>  =", V020,Q020
      write(6,111) "<1s1s'|r2^2 F0|Phi 0>  =", V002,Q002
      write(6,111) "<1s1s'|r1r2 F1|Phi 0>  =", V111,Q111

      call UMATRnm (U00, phip1s, phi1s, 0, 0)  
      call UMATRnm (U02, phip1s, phi1s, 0, 2)  
      call UMATRnm (U20, phip1s, phi1s, 2, 0)  

      call RnOVER(0,phi1s,phi1s,P0)
      call RnOVER(2,phi1s,phi1s,P2)
      
      write(6,111) "<1s1s'|     |Phi 0>  =", U00,P0**2
      write(6,111) "<1s1s'|r1^2 |Phi 0>  =", U02,P0*P2
      write(6,111) "<1s1s'|r2^2 |Phi 0>  =", U20,P0*P2
            
      F = U00 +   a1*2*(V020+V002-V111*2/3) + a1**2 * (U02+U20)
!      F = P0**2 + a1*4*(Q020-Q111/3)      + a1**2 * (P0*P2+P2*P0)

c 2nd Order integrals

      call VMATRnm (V040, phi1s, phip1s, 0, 4, 0)  
      call VMATRnm (V004, phi1s, phip1s, 0, 0, 4)  
      call VMATRnm (V022, phi1s, phip1s, 0, 2, 2)  
      call VMATRnm (V222, phi1s, phip1s, 2, 2, 2)  

      call UMATRnm (U04, phi1s, phip1s, 0, 4)  
      call UMATRnm (U40, phi1s, phip1s, 4, 0)  
      call UMATRnm (U22, phi1s, phip1s, 2, 2)  

      F = F + a5*2*(U20 + U02) + a5*a1*2*(V040+V004+V022*2 - V222*4/5) 
     >   + a5**2*(U04 + U40 + U22*10/3)
      
c 3d Order integrals

      call VMATRnm (V060, phi1s, phip1s, 0, 6, 0)  
      call VMATRnm (V006, phi1s, phip1s, 0, 0, 6)  
      call VMATRnm (V042, phi1s, phip1s, 0, 4, 2)  
      call VMATRnm (V024, phi1s, phip1s, 0, 2, 4)  
      call VMATRnm (V133, phi1s, phip1s, 1, 3, 3)  
      call VMATRnm (V333, phi1s, phip1s, 3, 3, 3)  

      call UMATRnm (U06, phi1s, phip1s, 0, 6)  
      call UMATRnm (U60, phi1s, phip1s, 6, 0)  
      call UMATRnm (U24, phi1s, phip1s, 2, 4)  
      call UMATRnm (U42, phi1s, phip1s, 4, 2)  

      F = F + a8*2*(V040+V004+V022*2 - V222*4/5) 
     >   + a8*a1*2*(U04 + U40 + U22*10/3)
     >   + a8*a5*2*(V060+V042*5+V024*5 + V006- V133*2-V333*6/7) 
     >   + a8**2*(U60+U42*7+U24*7 + U06) 
      

      F = (pi*4)**2 * F

      write(6,'(A,4E13.4)') 'Norma ', F

C                             ___
C  Radial orbitals || 1s> =  V4pi  |1s>
                                              
      do i = 1, meshr
         phi1s(i)  = phi1s(i) *sqrt(pi*4)!/ sqrt(F)
         phip1s(i) = phip1s(i)*sqrt(pi*4)!/ sqrt(F)
      end do

C Asymptotics r2->0 r1->oo r12->r2

      
      do i = 1, meshr
         rt = RMESH(i,1)
         tmp =   exp(-(Z-alam)*rt) * (1
     >      + a1*rt + a5*rt**2
     >      + a8*rt**3 + a13*rt**4
     >      + a14*rt**5)
         write(77,'(3E13.4)') rt, tmp, exp(-z*rt)
      end do

      

  20  FORMAT (2  F17.7)

      END      

      SUBROUTINE VMATRnm (result, fu1, fu2, LV, N, M)

C Calculates radial part of the symmetrized matrix element     
C with Pluvenage ground state 1s1s'+1s'1s

C          N  M                       N         M
C <f1,f2||r  r  V ||1s'1s+1s1s'> ~ R(r f1,1s'; r f2,1s) + 1s<->1s'  
C          1  2  L                  L                                   
C
C Functions f1,f2 should be free of any Simpsons wheigts

      include 'par.f'

      common /meshrr/ meshr,rmesh(maxr,3)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     :   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /plu/    phi1s(maxr),phip1s(maxr),maxphi
           
      DIMENSION fu1(maxr), fu2(maxr)
      DIMENSION fun(maxr), temp(maxr)

      result=0.
      if(Lv.lt.0)then !No calculation for Lv<0
         return
      end if

! Direct term 1s'1s
      
      minfun = 1
      maxfun = maxphi
      do i = 1, maxphi
         fun(i) = fu1(i) * phip1s(i) * rmesh(i,3)
     :                              * rmesh(i,1)**N
      end do

      call form(fun,minfun,maxfun,rpow1(1,lv),rpow2(1,lv),
     :      minrp(lv),maxrp(lv),meshr,temp,i1,i2)
                       
      DO I = 1, maxphi
         result = result + temp(i) * fu2(i) * phi1s(i) * rmesh(i,3)
     :                                                  * rmesh(i,1)**M
      end do
      
! Exchange term 1s1s'

      do i = 1, maxphi
         fun(i) = fu1(i) * phi1s(i) * rmesh(i,3)
     :                              * rmesh(i,1)**N
      end do

      call form(fun,minfun,maxfun,rpow1(1,lv),rpow2(1,lv),
     :      minrp(lv),maxrp(lv),meshr,temp,i1,i2)
                       
      DO I = 1, maxphi
         result = result + temp(i) * fu2(i) * phip1s(i) * rmesh(i,3)
     :                                                  * rmesh(i,1)**M
      end do
      
!      result = 0.5*result
      
      RETURN
      END


      SUBROUTINE UMATRnm (result, fu1, fu2, N, M)

C Calculates radial part of the symmetrized matrix element     
C with Pluvenage ground state 1s1s'+1s'1s
C                                  oo           oo                   
C          N  M                    f  N         f  M               
C <f1,f2||r  r   ||1s'1s+1s1s'> ~  | r dr f1 1s | r dr f2 1s' + 1s<->1s'
C          1  2                    j            j                  
C                                  o            o                  
C
C Functions f1,f2 should be free of any Simpsons wheigts

      include 'par.f'

      DIMENSION fu1(maxr), fu2(maxr)
      DIMENSION fun(maxr), temp(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /plu/    phi1s(maxr),phip1s(maxr),maxphi
      
! Direct term 1s'1s
      tmp = 0.d0
      do i = 1, maxphi
         rt = rmesh(i,1)
         sw = rmesh(i,3) 
         tmp = tmp + fu1(i)*phi1s(i)*rt**N * sw
      end do 
      tmp1 = 0.d0
      do i = 1, maxphi 
         rt = rmesh(i,1)
         sw = rmesh(i,3)  
         tmp1 = tmp1 + fu2(i)*phip1s(i)*rt**M * sw
      end do 
      result = tmp*tmp1

! Exchange term 1s1s'
      tmp = 0.d0
      do i = 1, maxphi
         rt = rmesh(i,1)
         sw = rmesh(i,3) 
         tmp = tmp + fu1(i)*phip1s(i)*rt**N * sw
      end do 
      tmp1 = 0.d0
      do i = 1, maxphi 
         rt = rmesh(i,1)
         sw = rmesh(i,3)  
         tmp1 = tmp1 + fu2(i)*phi1s(i)*rt**M * sw
      end do 
      result = result + tmp*tmp1

!      result = 0.5*result
      
      RETURN
      END
      

      subroutine pDIPOLE(w,w2,l,l1,r1,v1,x1,r,v,x) 
                                                       
      include 'par.f'                                            
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
           
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

      call pLENGTH  (f,f1, 0,1, rmatr)           
      call pVELOCITY(f,f1, 0,1, dmatr)           
      call pACCEL   (f,f1, 0,1, xmatr)           

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

      call pLENGTH  (f,f1, 1,2, rmatr)           
      call pVELOCITY(f,f1, 1,2, dmatr)           
      call pACCEL   (f,f1, 1,2, xmatr)           
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

      call pLENGTH  (f,f1, 2,3, rmatr)           
      call pVELOCITY(f,f1, 2,3, dmatr)           
      call pACCEL   (f,f1, 2,3, xmatr)           
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

      call pLENGTH  (f,f1, 3,4, rmatr)           
      call pVELOCITY(f,f1, 3,4, dmatr)           
      call pACCEL   (f,f1, 3,4, xmatr)           
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

      call pLENGTH  (f,f1, 4,5, rmatr)           
      call pVELOCITY(f,f1, 4,5, dmatr)           
      call pACCEL   (f,f1, 4,5, xmatr)           
 115  continue

      r =rmatr
      v =dmatr
      x =xmatr

      RETURN
      END
      
      subroutine pLENGTH(wg,wf,Lg,Lf,result)

C  Calculates dipole ME in length form     
C  from 2-term Pluvinage GS wavefunction for He     
C   r                               
C  D =<gf |r1 V| F > +  < gf |r2 V| F >
C                 L                  L 
C  Here F  is L-pole component of the Hylleraas      
C        L                                        

      include 'par.f'
      DIMENSION  wg(maxr), wf(maxr)
      COMMON /plu1/   Z, an, a1,a5,a8,a13,a14
      if(iabs(Lf-Lg).ne.1) STOP 'WRONG Lf,Lg' 
      Lm=Lg                                   
      if(Lf.gt.Lg) Lm=Lf                       
      const=float((-1)**(Lm))*sqrt(float(Lm)) 

C Initialization
      
         call UMATRnm (U01,  wg, wf, 0, 1)  
         call UMATRnm (U03,  wg, wf, 0, 3)  
         call UMATRnm (U05,  wg, wf, 0, 5)  
         call UMATRnm (U10,  wg, wf, 1, 0)  
         call UMATRnm (U12,  wg, wf, 1, 2)  
         call UMATRnm (U14,  wg, wf, 1, 4)  
         call UMATRnm (U21,  wg, wf, 2, 1)  
         call UMATRnm (U23,  wg, wf, 2, 3)  
         call UMATRnm (U30,  wg, wf, 3, 0)  
         call UMATRnm (U32,  wg, wf, 3, 2)  
         call UMATRnm (U41,  wg, wf, 4, 1)  
         call UMATRnm (U50,  wg, wf, 5, 0)  

         call VMATRnm (p003, wg, wf, 0, 0, 3)
         call VMATRnm (p004, wg, wf, 0, 0, 4)
         call VMATRnm (p005, wg, wf, 0, 0, 5)
         call VMATRnm (p006, wg, wf, 0, 0, 6)
         call VMATRnm (p012, wg, wf, 0, 1, 2)
         call VMATRnm (p013, wg, wf, 0, 1, 3)
         call VMATRnm (p014, wg, wf, 0, 1, 4)
         call VMATRnm (p015, wg, wf, 0, 1, 5)
         call VMATRnm (p021, wg, wf, 0, 2, 1)
         call VMATRnm (p022, wg, wf, 0, 2, 2)
         call VMATRnm (p023, wg, wf, 0, 2, 3)
         call VMATRnm (p024, wg, wf, 0, 2, 4)
         call VMATRnm (p030, wg, wf, 0, 3, 0)
         call VMATRnm (p031, wg, wf, 0, 3, 1)
         call VMATRnm (p032, wg, wf, 0, 3, 2)
         call VMATRnm (p040, wg, wf, 0, 4, 0)
         call VMATRnm (p041, wg, wf, 0, 4, 1)
         call VMATRnm (p042, wg, wf, 0, 4, 2)
         call VMATRnm (p050, wg, wf, 0, 5, 0)
         call VMATRnm (p051, wg, wf, 0, 5, 1)
         call VMATRnm (p060, wg, wf, 0, 6, 0)
         call VMATRnm (p112, wg, wf, 1, 1, 2)
         call VMATRnm (p113, wg, wf, 1, 1, 3)
         call VMATRnm (p114, wg, wf, 1, 1, 4)
         call VMATRnm (p115, wg, wf, 1, 1, 5)
         call VMATRnm (p121, wg, wf, 1, 2, 1)
         call VMATRnm (p122, wg, wf, 1, 2, 2)
         call VMATRnm (p123, wg, wf, 1, 2, 3)
         call VMATRnm (p124, wg, wf, 1, 2, 4)
         call VMATRnm (p131, wg, wf, 1, 3, 1)
         call VMATRnm (p132, wg, wf, 1, 3, 2)
         call VMATRnm (p133, wg, wf, 1, 3, 3)
         call VMATRnm (p141, wg, wf, 1, 4, 1)
         call VMATRnm (p142, wg, wf, 1, 4, 2)
         call VMATRnm (p151, wg, wf, 1, 5, 1)
         call VMATRnm (p223, wg, wf, 2, 2, 3)
         call VMATRnm (p232, wg, wf, 2, 3, 2)
         call VMATRnm (p332, wg, wf, 3, 3, 2)
         call VMATRnm (p323, wg, wf, 3, 2, 3)
         
         call VMATRnm (p070, wg, wf, 0, 7, 0)
         call VMATRnm (p061, wg, wf, 0, 6, 1)
         call VMATRnm (p052, wg, wf, 0, 5, 2)
         call VMATRnm (p043, wg, wf, 0, 4, 3)
         call VMATRnm (p034, wg, wf, 0, 3, 4)
         call VMATRnm (p025, wg, wf, 0, 2, 5)
         call VMATRnm (p016, wg, wf, 0, 1, 6)
         call VMATRnm (p007, wg, wf, 0, 0, 7)
         
         call VMATRnm (p152, wg, wf, 1, 5, 2)
         call VMATRnm (p143, wg, wf, 1, 4, 3)
         call VMATRnm (p134, wg, wf, 1, 3, 4)
         call VMATRnm (p125, wg, wf, 1, 2, 5)
         
         call VMATRnm (p252, wg, wf, 2, 5, 2)
         call VMATRnm (p243, wg, wf, 2, 4, 3)
         call VMATRnm (p234, wg, wf, 2, 3, 4)
         call VMATRnm (p225, wg, wf, 2, 2, 5)
         
         call VMATRnm (p352, wg, wf, 3, 5, 2)
         call VMATRnm (p343, wg, wf, 3, 4, 3)
         call VMATRnm (p334, wg, wf, 3, 3, 4)
         call VMATRnm (p325, wg, wf, 3, 2, 5)
         
         call VMATRnm (p452, wg, wf, 4, 5, 2)
         call VMATRnm (p443, wg, wf, 4, 4, 3)
         call VMATRnm (p434, wg, wf, 4, 3, 4)
         call VMATRnm (p425, wg, wf, 4, 2, 5)

         call VMATRnm (p543, wg, wf, 5, 4, 3)
         call VMATRnm (p534, wg, wf, 5, 3, 4)
     
C First term D1: L=Lf->Lg transition

      L = Lf
      IF(L .eq. 0) then

         D1  = U10 + a1 * (p030  + p012  - p121*2/3)
     :             + a5 * (U30 + U12)          
     : + a8 * (p050 + p014 + p032*2 - p232*4/5)
     : + a13* (U14 + U32*10/3 + U50)
     : + a14 * (p070 + 5*p052 - 2*p143 - 6*p343/7 + 5*p034 + p016)

      ELSE

         call VMATRnm (s21, wg, wf, L+1, 2, 1)
         call VMATRnm (s22, wg, wf, L+1, 2, 2)
         call VMATRnm (s23, wg, wf, L+1, 2, 3)
         call VMATRnm (s31, wg, wf, L+1, 3, 1)
         call VMATRnm (s32, wg, wf, L+1, 3, 2)
         call VMATRnm (s41, wg, wf, L+1, 4, 1)
         call VMATRnm (t21, wg, wf, L-1, 2, 1)
         call VMATRnm (t22, wg, wf, L-1, 2, 2)
         call VMATRnm (t23, wg, wf, L-1, 2, 3)
         call VMATRnm (t31, wg, wf, L-1, 3, 1)
         call VMATRnm (t32, wg, wf, L-1, 3, 2)
         call VMATRnm (t41, wg, wf, L-1, 4, 1)

         call VMATRnm (P, wg, wf, L+2, 3, 2)
         call VMATRnm (Q, wg, wf, L  , 3, 2)
         call VMATRnm (R, wg, wf, L-2, 3, 2)
         
         D1 = a1 * (s21/(2*L+3) - t21/(2*L-1)) 
         B8 = a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))
         call VMATRnm (P, wg, wf, L+3, 4, 3)
         call VMATRnm (Q, wg, wf, L+1, 4, 3)
         call VMATRnm (R, wg, wf, L-1, 4, 3)
         call VMATRnm (S, wg, wf, L-3, 4, 3)

         B14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         IF(L .eq. 1) THEN
            D1  = D1 - a5 * U21*2
     :          + a13*4*(-U23 - U41)                              
            B8 = -a8 *(p041 + p023 - p132*3/5 - p332/35)*3 
            B14= a14 *(-5*p061 - 10*p043 - 5*p025
     :               +  9*p152/5 + 3*p352/35 + 9*p134/5 + 3*p334/35
     :               -  4*p043/5 + 64*p243/35 - 4*p443/105)
         END IF

         IF(L .eq. 2) THEN
            D1  = D1
     :          + a13* 8/3 *  U32
            B14 = a14 *(5*p052 - 2*p252/7 + p452/21
     :                + 5*p034 - 2*p234/7 + p434/21
     :              -  18*p143/7 - 2*p543/77)
         END IF
         D1  = D1 + B8 + B14 
      END IF
         
         
C Second term D2: L=Lg->Lf transition

      L = Lg
      IF(L .eq. 0) then
         D2  = U01 + a1 * (p021  + p003  - p112*2/3)
     :               + a5 * (U03 + U21)            
     : + a8 * (p041 + p005 + p023*2 - p223*4/5)
     : + a13* (U05 + U23*10/3 + U41)
     : + a14 * (p061 + 5*p043 - 2*p134 - 6*p334/7 + 5*p025 + p007)
      ELSE 
         call VMATRnm (s12, wg, wf, L+1, 1, 2)
         call VMATRnm (s13, wg, wf, L+1, 1, 3)
         call VMATRnm (s14, wg, wf, L+1, 1, 4)
         call VMATRnm (s22, wg, wf, L+1, 2, 2)
         call VMATRnm (s23, wg, wf, L+1, 2, 3)
         call VMATRnm (s32, wg, wf, L+1, 3, 2)
         call VMATRnm (t12, wg, wf, L-1, 1, 2)
         call VMATRnm (t13, wg, wf, L-1, 1, 3)
         call VMATRnm (t14, wg, wf, L-1, 1, 4)
         call VMATRnm (t22, wg, wf, L-1, 2, 2)
         call VMATRnm (t23, wg, wf, L-1, 2, 3)
         call VMATRnm (t32, wg, wf, L-1, 3, 2)

         call VMATRnm (P, wg, wf, L+2, 2, 3)
         call VMATRnm (Q, wg, wf, L  , 2, 3)
         call VMATRnm (R, wg, wf, L-2, 2, 3)
         
         D2 = a1 * (s12/(2*L+3) - t12/(2*L-1)) 
         C8 = a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))
         call VMATRnm (P, wg, wf, L+3, 3, 4)
         call VMATRnm (Q, wg, wf, L+1, 3, 4)
         call VMATRnm (R, wg, wf, L-1, 3, 4)
         call VMATRnm (S, wg, wf, L-3, 3, 4)

         C14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         IF(L .eq. 1) THEN
            D2  = D2 - a5 * U12*2
     :        + a13*4 * (-U14 - U32)                              
            C8 =     - a8 *(p032 + p014 - p123*3/5 - p323/35)*3 
           C14 =      a14 *(-5*p052 - 10*p034 - 5*p016
     :                   +  9*p143/5 + 3*p343/35 + 9*p125/5 + 3*p325/35
     :                   -  4*p034/5 + 64*p234/35 - 4*p434/105)
         END IF

         IF(L .eq. 2) THEN
            D2  = D2
     :          + a13* 8/3 *  U23 
            C14 = a14 *(5*p043 - 2*p243/7 + p443/21 
     :               +  5*p025 - 2*p225/7 + p425/21
     :               - 18*p134/7 - 2*p534/77)
         END IF
         D2  = D2 + C8 + C14 
      END IF

      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const


      RETURN
      END


      subroutine pVELOCITY(wg,wf,lg,lf,result)     
                                             
C  Calculates dipole ME in velocity form     
C  from 2-term Pluvinage GS wavefunction for He     
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
      COMMON /plu1/   Z, an, a1,a5,a8,a13,a14
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

      call UMATRnm (U00,  YD, wf, 0, 0)
      call UMATRnm (U20,  YD, wf, 2, 0)
      call UMATRnm (U02,  YD, wf, 0, 2)
      call UMATRnm (U40,  YD, wf, 4, 0)
      call UMATRnm (U22,  YD, wf, 2, 2)
      call UMATRnm (U04,  YD, wf, 0, 4)
      call UMATRnm (U11,  YD, wf, 1, 1)
      call UMATRnm (U31,  YD, wf, 3, 1)
      call UMATRnm (U13,  YD, wf, 1, 3)
      
      call VMATRnm (p002, YD, wf, 0, 0, 2)
      call VMATRnm (p003, YD, wf, 0, 0, 3)
      call VMATRnm (p004, YD, wf, 0, 0, 4)
      call VMATRnm (p005, YD, wf, 0, 0, 5)
      call VMATRnm (p012, YD, wf, 0, 1, 2)
      call VMATRnm (p013, YD, wf, 0, 1, 3)
      call VMATRnm (p014, YD, wf, 0, 1, 4)
      call VMATRnm (p020, YD, wf, 0, 2, 0)
      call VMATRnm (p021, YD, wf, 0, 2, 1)
      call VMATRnm (p022, YD, wf, 0, 2, 2)
      call VMATRnm (p030, YD, wf, 0, 3, 0)
      call VMATRnm (p031, YD, wf, 0, 3, 1)
      call VMATRnm (p040, YD, wf, 0, 4, 0)
      call VMATRnm (p041, YD, wf, 0, 4, 1)
      call VMATRnm (p050, YD, wf, 0, 5, 0)
      call VMATRnm (p111, YD, wf, 1, 1, 1)
      call VMATRnm (p112, YD, wf, 1, 1, 2)
      call VMATRnm (p113, YD, wf, 1, 1, 3)
      call VMATRnm (p121, YD, wf, 1, 2, 1)
      call VMATRnm (p122, YD, wf, 1, 2, 2)
      call VMATRnm (p131, YD, wf, 1, 3, 1)
      call VMATRnm (p222, YD, wf, 2, 2, 2)
      call VMATRnm (p322, YD, wf, 3, 2, 2)
      
      call VMATRnm (p006, YD, wf, 0, 0, 6)
      call VMATRnm (p015, YD, wf, 0, 1, 5)
      call VMATRnm (p024, YD, wf, 0, 2, 4)
      call VMATRnm (p033, YD, wf, 0, 3, 3)
      call VMATRnm (p042, YD, wf, 0, 4, 2)
      call VMATRnm (p051, YD, wf, 0, 5, 1)
      call VMATRnm (p060, YD, wf, 0, 6, 0)
      
      call VMATRnm (p114, YD, wf, 1, 1, 4)
      call VMATRnm (p132, YD, wf, 1, 3, 2)
      call VMATRnm (p133, YD, wf, 1, 3, 3)
      call VMATRnm (p123, YD, wf, 1, 2, 3)
      call VMATRnm (p124, YD, wf, 1, 2, 4)
      call VMATRnm (p141, YD, wf, 1, 4, 1)
      call VMATRnm (p142, YD, wf, 1, 4, 2)
      
      call VMATRnm (p233, YD, wf, 2, 3, 3)
      call VMATRnm (p224, YD, wf, 2, 2, 4)
      call VMATRnm (p242, YD, wf, 2, 4, 2)
      
      call VMATRnm (p333, YD, wf, 3, 3, 3)
      call VMATRnm (p324, YD, wf, 3, 2, 4)
      call VMATRnm (p342, YD, wf, 3, 4, 2)
      
      call VMATRnm (p433, YD, wf, 4, 3, 3)
      call VMATRnm (p424, YD, wf, 4, 2, 4)
      call VMATRnm (p442, YD, wf, 4, 4, 2)

      call VMATRnm (p533, YD, wf, 5, 3, 3)


      IF(L .eq. 0) then

         D1  = U00 + a1 * (p020 + p002 - p111*2/3)
     :               + a5 * (U20 + U02)            
     : + a8 * (p040 + p004 + p022*2 - p222*4/5)
     : + a13* (U40 + U22*10/3 + U04)
     : + a14 * (p060 + 5*p042 + 5*p024 + p006
     :       -   2*p133 - 6*p333/7)
      ELSE
         call VMATRnm (s11, YD, wf, L+1, 1, 1)
         call VMATRnm (s12, YD, wf, L+1, 1, 2)
         call VMATRnm (s21, YD, wf, L+1, 2, 1)
         call VMATRnm (s13, YD, wf, L+1, 1, 3)
         call VMATRnm (s31, YD, wf, L+1, 3, 1)
         call VMATRnm (s22, YD, wf, L+1, 2, 2)
         call VMATRnm (t11, YD, wf, L-1, 1, 1)
         call VMATRnm (t12, YD, wf, L-1, 1, 2)
         call VMATRnm (t21, YD, wf, L-1, 2, 1)
         call VMATRnm (t13, YD, wf, L-1, 1, 3)
         call VMATRnm (t31, YD, wf, L-1, 3, 1)
         call VMATRnm (t22, YD, wf, L-1, 2, 2)

         call VMATRnm (P, YD, wf, L+2, 2, 2)
         call VMATRnm (Q, YD, wf, L  , 2, 2)
         call VMATRnm (R, YD, wf, L-2, 2, 2)
         
         D1 = a1 * (s11/(2*L+3) - t11/(2*L-1)) 
         B8= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2
     :                                 + R/(2*L-1)/(2*L-3))

         call VMATRnm (P, YD, wf, L+3, 3, 3)
         call VMATRnm (Q, YD, wf, L+1, 3, 3)
         call VMATRnm (R, YD, wf, L-1, 3, 3)
         call VMATRnm (S, YD, wf, L-3, 3, 3)
                                    
         B14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 


         IF(L .eq. 1) THEN
            D1 =D1 - a5 * U11*2
     :             + a13 *4*(-U31 - U13)
            B8=    - a8  *(p031 + p013 - p122*3/5 - p322/35)*3 
            B14=     a14 *(-5*p051 - 10*p033 - 5*p015
     :                 +  9*p142/5 + 3*p342/35 +  9*p124/5 + 3*p324/35
     :                   -  4*p033/5 + 64*p233/35 - 4*p433/105)
         END IF

         IF(L .eq. 2) THEN
            D1  = D1
     :         +     a13 * 8/3 *  U22
            B14  =   a14 *(5*p042 - 2*p242/7 + p442/21
     :                    +5*p024 - 2*p224/7 + p424/21
     :                  - 18*p133/7 - 2*p533/77)
        END IF
        D1  = D1 + B8 + B14
      END IF
         
C Second term D2: L=Lg->Lf transition

      L = Lg
      if(Lf.gt.Lg) IL= 1                      !Lf=Lg+1
      if(Lf.lt.Lg) IL=-1                      !Lf=Lg-1
      call DERIVATIVE(wf,wd)
      DO i=1,meshr
         YD(i) = wd(i) + wf(i)/rmesh(i,1)*il*Lm
      end do

      call UMATRnm (U00,  wg, YD, 0, 0)
      call UMATRnm (U20,  wg, YD, 2, 0)
      call UMATRnm (U02,  wg, YD, 0, 2)
      call UMATRnm (U40,  wg, YD, 4, 0)
      call UMATRnm (U22,  wg, YD, 2, 2)
      call UMATRnm (U04,  wg, YD, 0, 4)
      call UMATRnm (U11,  wg, YD, 1, 1)
      call UMATRnm (U31,  wg, YD, 3, 1)
      call UMATRnm (U13,  wg, YD, 1, 3)
      
      call VMATRnm (p002, wg, YD, 0, 0, 2)
      call VMATRnm (p003, wg, YD, 0, 0, 3)
      call VMATRnm (p004, wg, YD, 0, 0, 4)
      call VMATRnm (p005, wg, YD, 0, 0, 5)
      call VMATRnm (p012, wg, YD, 0, 1, 2)
      call VMATRnm (p013, wg, YD, 0, 1, 3)
      call VMATRnm (p014, wg, YD, 0, 1, 4)
      call VMATRnm (p020, wg, YD, 0, 2, 0)
      call VMATRnm (p021, wg, YD, 0, 2, 1)
      call VMATRnm (p022, wg, YD, 0, 2, 2)
      call VMATRnm (p030, wg, YD, 0, 3, 0)
      call VMATRnm (p031, wg, YD, 0, 3, 1)
      call VMATRnm (p040, wg, YD, 0, 4, 0)
      call VMATRnm (p041, wg, YD, 0, 4, 1)
      call VMATRnm (p050, wg, YD, 0, 5, 0)
      call VMATRnm (p111, wg, YD, 1, 1, 1)
      call VMATRnm (p112, wg, YD, 1, 1, 2)
      call VMATRnm (p113, wg, YD, 1, 1, 3)
      call VMATRnm (p121, wg, YD, 1, 2, 1)
      call VMATRnm (p122, wg, YD, 1, 2, 2)
      call VMATRnm (p131, wg, YD, 1, 3, 1)
      call VMATRnm (p222, wg, YD, 2, 2, 2)
      call VMATRnm (p322, wg, YD, 3, 2, 2)
      
      call VMATRnm (p006, wg, YD, 0, 0, 6)
      call VMATRnm (p015, wg, YD, 0, 1, 5)
      call VMATRnm (p024, wg, YD, 0, 2, 4)
      call VMATRnm (p033, wg, YD, 0, 3, 3)
      call VMATRnm (p042, wg, YD, 0, 4, 2)
      call VMATRnm (p051, wg, YD, 0, 5, 1)
      call VMATRnm (p060, wg, YD, 0, 6, 0)
      
      call VMATRnm (p114, wg, YD, 1, 1, 4)
      call VMATRnm (p132, wg, YD, 1, 3, 2)
      call VMATRnm (p133, wg, YD, 1, 3, 3)
      call VMATRnm (p123, wg, YD, 1, 2, 3)
      call VMATRnm (p124, wg, YD, 1, 2, 4)
      call VMATRnm (p141, wg, YD, 1, 4, 1)
      call VMATRnm (p142, wg, YD, 1, 4, 2)
      
      call VMATRnm (p233, wg, YD, 2, 3, 3)
      call VMATRnm (p224, wg, YD, 2, 2, 4)
      call VMATRnm (p242, wg, YD, 2, 4, 2)
      
      call VMATRnm (p333, wg, YD, 3, 3, 3)
      call VMATRnm (p324, wg, YD, 3, 2, 4)
      call VMATRnm (p342, wg, YD, 3, 4, 2)
      
      call VMATRnm (p433, wg, YD, 4, 3, 3)
      call VMATRnm (p424, wg, YD, 4, 2, 4)
      call VMATRnm (p442, wg, YD, 4, 4, 2)
      
      call VMATRnm (p533, wg, YD, 5, 3, 3)
      
      IF(L .eq. 0) then

         D2  = U00 + a1 * (p020 + p002 - p111*2/3)
     :               + a5 * (U20 + U02)            
     : + a8 * (p040 + p004 + p022*2 - p222*4/5)
     : + a13* (U40 + U22*10/3 + U04)
     : + a14 * (p060 + 5*p042 + 5*p024 + p006
     :       -   2*p133 - 6*p333/7)
      ELSE
         call VMATRnm (s11, wg, YD, L+1, 1, 1)
         call VMATRnm (s12, wg, YD, L+1, 1, 2)
         call VMATRnm (s21, wg, YD, L+1, 2, 1)
         call VMATRnm (s13, wg, YD, L+1, 1, 3)
         call VMATRnm (s31, wg, YD, L+1, 3, 1)
         call VMATRnm (s22, wg, YD, L+1, 2, 2)
         call VMATRnm (t11, wg, YD, L-1, 1, 1)
         call VMATRnm (t12, wg, YD, L-1, 1, 2)
         call VMATRnm (t21, wg, YD, L-1, 2, 1)
         call VMATRnm (t13, wg, YD, L-1, 1, 3)
         call VMATRnm (t31, wg, YD, L-1, 3, 1)
         call VMATRnm (t22, wg, YD, L-1, 2, 2)

         call VMATRnm (P, wg, YD, L+2, 2, 2)
         call VMATRnm (Q, wg, YD, L  , 2, 2)
         call VMATRnm (R, wg, YD, L-2, 2, 2)
         
         D2 = a1 * (s11/(2*L+3) - t11/(2*L-1)) 
         C8= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2
     :                                 + R/(2*L-1)/(2*L-3))
         call VMATRnm (P, wg, YD, L+3, 3, 3)
         call VMATRnm (Q, wg, YD, L+1, 3, 3)
         call VMATRnm (R, wg, YD, L-1, 3, 3)
         call VMATRnm (S, wg, YD, L-3, 3, 3)
                                    
         C14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         IF(L .eq. 1) THEN
            D2  = D2 - a5 * U11*2
     :  +  a13 * 4 * (-U31 - U13)
           C8=     - a8 *(p031 + p013 - p122*3/5 - p322/35)*3 
           C14=     a14 *(-5*p051 - 10*p033 - 5*p015
     :                  +  9*p142/5 + 3*p342/35 +  9*p124/5 + 3*p324/35
     :                  -  4*p033/5 + 64*p233/35 - 4*p433/105)
         END IF
         IF(L .eq. 2) THEN
            D2  = D2
     :         +     a13 * 8/3 *  U22
            C14  =   a14 *(5*p042 - 2*p242/7 + p442/21
     :                    +5*p024 - 2*p224/7 + p424/21
     :                  - 18*p133/7 - 2*p533/77)
         END IF
         D2  = D2 + C8 + C14
      END IF
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const

      RETURN
      END

      subroutine pACCEL(wg,wf,Lg,Lf,result)

C  Calculates dipole ME in acceleration form     
C  from 2-term Pluvinage GS wavefunction for He     
C   .                                 
C   V           -2                   -2   
C  D = 2 <gf |r1  V| F > + 2 < gf |r2  V| F >
C                     L                    L 
C  Here F  is L-pole component of the Hylleraas      
C        L             
      include 'par.f'
      DIMENSION  wg(maxr), wf(maxr)
      COMMON /plu1/   Z, an, a1,a5,a8,a13,a14
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym

      if(iabs(Lf-Lg).ne.1) STOP 'WRONG Lf,Lg' 
      Lm=Lg                                   
      if(Lf.gt.Lg) Lm=Lf                       
      const=float((-1)**(Lm))*sqrt(float(Lm)) 

      const = const*nznuc

C Initialization

      call UMATRnm (U00,  wg, wf, 0, 0)
      call UMATRnm (U02,  wg, wf, 0, 2)
      call UMATRnm (U04,  wg, wf, 0, 4)
      call UMATRnm (U11,  wg, wf, 1, 1)
      call UMATRnm (U13,  wg, wf, 1, 3)
      call UMATRnm (U20,  wg, wf, 2, 0)
      call UMATRnm (U22,  wg, wf, 2, 2)
      call UMATRnm (U_11, wg, wf,-1, 1)
      call UMATRnm (U_13, wg, wf,-1, 3)
      call UMATRnm (U_20, wg, wf,-2, 0)
      call UMATRnm (U_22, wg, wf,-2, 2)
      call UMATRnm (U_24, wg, wf,-2, 4)
      call UMATRnm (U0_2, wg, wf, 0,-2)
      call UMATRnm (U1_1, wg, wf, 1,-1)
      call UMATRnm (U2_2, wg, wf, 2,-2)
      call UMATRnm (U3_1, wg, wf, 3,-1)
      call UMATRnm (U4_2, wg, wf, 4,-2)

      call VMATRnm (p0_22,wg, wf, 0,-2, 2)
      call VMATRnm (p0_23,wg, wf, 0,-2, 3)
      call VMATRnm (p0_24,wg, wf, 0,-2, 4)
      call VMATRnm (p0_25,wg, wf, 0,-2, 5)
      call VMATRnm (p0_12,wg, wf, 0,-1, 2)
      call VMATRnm (p0_13,wg, wf, 0,-1, 3)
      call VMATRnm (p0_14,wg, wf, 0,-1, 4)
      call VMATRnm (p000, wg, wf, 0, 0, 0)
      call VMATRnm (p001, wg, wf, 0, 0, 1)
      call VMATRnm (p002, wg, wf, 0, 0, 2)
      call VMATRnm (p003, wg, wf, 0, 0, 3)
      call VMATRnm (p010, wg, wf, 0, 1, 0)
      call VMATRnm (p011, wg, wf, 0, 1, 1)
      call VMATRnm (p012, wg, wf, 0, 1, 2)
      call VMATRnm (p02_2,wg, wf, 0, 2,-2)
      call VMATRnm (p02_1,wg, wf, 0, 2,-1)
      call VMATRnm (p020, wg, wf, 0, 2, 0)
      call VMATRnm (p021, wg, wf, 0, 2, 1)
      call VMATRnm (p030, wg, wf, 0, 3, 0)
      call VMATRnm (p03_2,wg, wf, 0, 3,-2)
      call VMATRnm (p03_1,wg, wf, 0, 3,-1)
      call VMATRnm (p030, wg, wf, 0, 3, 0)
      call VMATRnm (p04_1,wg, wf, 0, 4,-1)
      call VMATRnm (p04_2,wg, wf, 0, 4,-2)
      call VMATRnm (p05_2,wg, wf, 0, 5,-2)
      call VMATRnm (p1_11,wg, wf, 1,-1, 1)
      call VMATRnm (p1_12,wg, wf, 1,-1, 2)
      call VMATRnm (p1_13,wg, wf, 1,-1, 3)
      call VMATRnm (p1_14,wg, wf, 1,-1, 4)
      call VMATRnm (p101, wg, wf, 1, 0, 1)
      call VMATRnm (p102, wg, wf, 1, 0, 2)
      call VMATRnm (p103, wg, wf, 1, 0, 3)
      call VMATRnm (p11_1,wg, wf, 1, 1,-1)
      call VMATRnm (p110, wg, wf, 1, 1, 0)
      call VMATRnm (p111, wg, wf, 1, 1, 1)
      call VMATRnm (p112, wg, wf, 1, 1, 2)
      call VMATRnm (p12_1,wg, wf, 1, 2,-1)
      call VMATRnm (p120, wg, wf, 1, 2, 0)
      call VMATRnm (p121, wg, wf, 1, 2, 1)
      call VMATRnm (p130, wg, wf, 1, 3, 0)
      call VMATRnm (p13_1,wg, wf, 1, 3,-1)
      call VMATRnm (p14_1,wg, wf, 1, 4,-1)
      call VMATRnm (p202, wg, wf, 2, 0, 2)
      call VMATRnm (p220, wg, wf, 2, 2, 0)
      call VMATRnm (p302, wg, wf, 3, 0, 2)
      call VMATRnm (p320, wg, wf, 3, 2, 0)

      call VMATRnm (p0_26,wg, wf, 0,-2, 6)
      call VMATRnm (p0_15,wg, wf, 0,-1, 5)
      call VMATRnm (p004, wg, wf, 0, 0, 4)
      call VMATRnm (p013, wg, wf, 0, 1, 3)
      call VMATRnm (p022, wg, wf, 0, 2, 2)
      call VMATRnm (p031, wg, wf, 0, 3, 1)
      call VMATRnm (p040, wg, wf, 0, 4, 0)
      call VMATRnm (p05_1,wg, wf, 0, 5,-1)
      call VMATRnm (p06_2,wg, wf, 0, 6,-2)

      call VMATRnm (p122, wg, wf, 1, 2, 2)
      call VMATRnm (p113, wg, wf, 1, 1, 3)
      call VMATRnm (p131, wg, wf, 1, 3, 1)
      call VMATRnm (p104, wg, wf, 1, 0, 4)
      call VMATRnm (p140, wg, wf, 1, 4, 0)

      call VMATRnm (p222, wg, wf, 2, 2, 2)
      call VMATRnm (p213, wg, wf, 2, 1, 3)
      call VMATRnm (p231, wg, wf, 2, 3, 1)
      call VMATRnm (p204, wg, wf, 2, 0, 4)
      call VMATRnm (p240, wg, wf, 2, 4, 0)

      call VMATRnm (p322, wg, wf, 3, 2, 2)
      call VMATRnm (p313, wg, wf, 3, 1, 3)
      call VMATRnm (p331, wg, wf, 3, 3, 1)
      call VMATRnm (p304, wg, wf, 3, 0, 4)
      call VMATRnm (p340, wg, wf, 3, 4, 0)

      call VMATRnm (p422, wg, wf, 4, 2, 2)
      call VMATRnm (p413, wg, wf, 4, 1, 3)
      call VMATRnm (p431, wg, wf, 4, 3, 1)
      call VMATRnm (p404, wg, wf, 4, 0, 4)
      call VMATRnm (p440, wg, wf, 4, 4, 0)

      call VMATRnm (p513, wg, wf, 5, 1, 3)
      call VMATRnm (p531, wg, wf, 5, 3, 1)

C First term D1: L=Lf->Lg transition

      L = Lf
      IF(L .eq. 0) then
         D1  = U_20 + a1 * (p000  + p0_22  - p1_11*2/3)
     :                + a5 * (U00 + U_22)          
     : + a8 * (p020 + p0_24 + p002*2 - p202*4/5)
     : + a13* (U_24 + U02*10/3 + U20)
     : + a14 * (p040 + 5*p022 + 5*p004 + p0_26
     :       -  2*p113 - 6*p313/7)
      ELSE
         call VMATRnm (s_11,wg, wf, L+1,-1, 1)
         call VMATRnm (s_12,wg, wf, L+1,-1, 2)
         call VMATRnm (s_13,wg, wf, L+1,-1, 3)
         call VMATRnm (s01, wg, wf, L+1, 0, 1)
         call VMATRnm (s02, wg, wf, L+1, 0, 2)
         call VMATRnm (s11, wg, wf, L+1, 1, 1)
         call VMATRnm (t_11,wg, wf, L-1,-1, 1)
         call VMATRnm (t_12,wg, wf, L-1,-1, 2)
         call VMATRnm (t_13,wg, wf, L-1,-1, 3)
         call VMATRnm (t01, wg, wf, L-1, 0, 1)
         call VMATRnm (t02, wg, wf, L-1, 0, 2)
         call VMATRnm (t11, wg, wf, L-1, 1, 1)

         call VMATRnm (P, wg, wf, L+2, 0, 2)
         call VMATRnm (Q, wg, wf, L  , 0, 2)
         call VMATRnm (R, wg, wf, L-2, 0, 2)
         
         D1 = a1 * (s_11/(2*L+3) - t_11/(2*L-1)) 
         B8 = a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))

         call VMATRnm (P, wg, wf, L+3, 1, 3)
         call VMATRnm (Q, wg, wf, L+1, 1, 3)
         call VMATRnm (R, wg, wf, L-1, 1, 3)
         call VMATRnm (S, wg, wf, L-3, 1, 3)

         B14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         IF(L .eq. 1) THEN
            D1  = D1 - a5 * U_11*2
     :          + a13*4*(-U_13 - U11)                              
         B8=     - a8 *(p011 + p0_13 - p102*3/5 - p302/35)*3 
         B14=     a14 *(-5*p031 - 10*p013 - 5*p0_15
     :                 + 9*p122/5 + 3*p322/35 + 9*p104/5 + 3*p304/35
     :                 - 4*p013/5 + 64*p213/35 - 4*p413/105)

         END IF

         IF(L .eq. 2) THEN
            D1  = D1
     :          + a13* 8/3 *  U02
            B14 = a14 *(5*p022 - 2*p222/7 + p422/21 
     :                + 5*p004 - 2*p204/7 + p404/21
     :               - 18*p113/7 - 2*p513/77)
           END IF
         D1  = D1 + B8 + B14
      END IF
         
C Second term D2: L=Lg->Lf transition

      L = Lg
      IF(L .eq. 0) then
         D2  = U0_2 + a1 * (p02_2  + p000  - p11_1*2/3)
     :                + a5 * (U00 + U2_2)            
     : + a8 * (p04_2 + p002 + p020*2 - p220*4/5)
     : + a13* (U02 + U20*10/3 + U4_2)
     : + a14 * (p06_2 + 5*p040 + 5*p022 + p004
     :       - 2*p131 - 6*p331/7)
      ELSE 
         call VMATRnm (s1_1, wg, wf, L+1, 1,-1)
         call VMATRnm (s10,  wg, wf, L+1, 1, 0)
         call VMATRnm (s11,  wg, wf, L+1, 1, 1)
         call VMATRnm (s2_1, wg, wf, L+1, 2,-1)
         call VMATRnm (s20,  wg, wf, L+1, 2, 0)
         call VMATRnm (s3_1, wg, wf, L+1, 3,-1)
         call VMATRnm (t1_1, wg, wf, L-1, 1,-1)
         call VMATRnm (t10,  wg, wf, L-1, 1, 0)
         call VMATRnm (t11,  wg, wf, L-1, 1, 1)
         call VMATRnm (t2_1, wg, wf, L-1, 2,-1)
         call VMATRnm (t20,  wg, wf, L-1, 2, 0)
         call VMATRnm (t3_1, wg, wf, L-1, 3,-1)

         call VMATRnm (P, wg, wf, L+2, 2, 0)
         call VMATRnm (Q, wg, wf, L  , 2, 0)
         call VMATRnm (R, wg, wf, L-2, 2, 0)
         
         D2 = a1 * (s1_1/(2*L+3) - t1_1/(2*L-1)) 
         C8 = a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))
         call VMATRnm (P, wg, wf, L+3, 3, 1)
         call VMATRnm (Q, wg, wf, L+1, 3, 1)
         call VMATRnm (R, wg, wf, L-1, 3, 1)
         call VMATRnm (S, wg, wf, L-3, 3, 1)

         C14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         IF(L .eq. 1) THEN
            D2  = D2 - a5 * U1_1*2
     :       + a13*4 * (-U11 - U3_1)                              
             C8=     - a8 *(p03_1 + p011 - p120*3/5 - p320/35)*3 
             C14 =     a14*(-5*p05_1 - 10*p031 - 5*p013
     :                 +  9*p140/5 + 3*p340/35 + 9*p122/5 + 3*p322/35
     :                  - 4*p031/5 + 64*p231/35 - 4*p431/105)
         END IF

         IF(L .eq. 2) THEN
            D2  = D2
     :          + a13* 8/3 *  U20
            C14 = a14 *(5*p040 - 2*p240/7 + p440/21 
     :                + 5*p022 - 2*p222/7 + p422/21
     :               - 18*p131/7 - 2*p531/77)
         END IF
         D2  = D2 + C8 + C14 
      END IF

      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const
      RETURN
      END
      


      subroutine pbLENGTH(wg,wf,lg,lf,j,qq,result)     
                                             
C  Calculates Born ME in length form     
C  with Pluvinage GS wavefunction for He     
C                                    
C  D =< gf |M (qr_1) + M (qr_2)| F > 
C            j          j              
C
C  Here F  is L-pole component of the Hylleraas      
C        L             
C       M (qr) is the spherical Bessel function
C        j   

      include 'par.f'
      DIMENSION  YD(maxr)
      DIMENSION  wg(MAXR), wf(MAXR), wd(MAXR)
      common /meshrr/ meshr,rmesh(maxr,3)
      COMMON /plu1/   Z, an, a1,a5,a8,a13,a14
      real*8 x,ang
      
      T = triang(Lg, Lf, J)
      if(T.eq.0)  print *, 'WRONG Lf,Lg,J', Lf,Lg,J

      call IOT3(Lf,J,Lg,0,0,0,ang)
      ddff=hat(Lg)*hat(Lf)
      ang=ang*ddff
      
C First term D1: L=Lf->Lg transition

      L = Lf
      DO i=1,meshr
         RT=rmesh(i,1)
         x = QQ*RT
         BE= BES(J,x)
         if(J.eq.0)BE= BES(J,x)-1.0
         YD(i) = wg(i)*BE
      end do

      call UMATRnm (U00,  YD, wf, 0, 0)
      call UMATRnm (U20,  YD, wf, 2, 0)
      call UMATRnm (U02,  YD, wf, 0, 2)
      call UMATRnm (U40,  YD, wf, 4, 0)
      call UMATRnm (U22,  YD, wf, 2, 2)
      call UMATRnm (U04,  YD, wf, 0, 4)
      call UMATRnm (U11,  YD, wf, 1, 1)
      call UMATRnm (U31,  YD, wf, 3, 1)
      call UMATRnm (U13,  YD, wf, 1, 3)

      call VMATRnm (p002, YD, wf, 0, 0, 2)
      call VMATRnm (p003, YD, wf, 0, 0, 3)
      call VMATRnm (p004, YD, wf, 0, 0, 4)
      call VMATRnm (p005, YD, wf, 0, 0, 5)
      call VMATRnm (p012, YD, wf, 0, 1, 2)
      call VMATRnm (p013, YD, wf, 0, 1, 3)
      call VMATRnm (p014, YD, wf, 0, 1, 4)
      call VMATRnm (p020, YD, wf, 0, 2, 0)
      call VMATRnm (p021, YD, wf, 0, 2, 1)
      call VMATRnm (p022, YD, wf, 0, 2, 2)
      call VMATRnm (p030, YD, wf, 0, 3, 0)
      call VMATRnm (p031, YD, wf, 0, 3, 1)
      call VMATRnm (p040, YD, wf, 0, 4, 0)
      call VMATRnm (p041, YD, wf, 0, 4, 1)
      call VMATRnm (p050, YD, wf, 0, 5, 0)
      call VMATRnm (p111, YD, wf, 1, 1, 1)
      call VMATRnm (p112, YD, wf, 1, 1, 2)
      call VMATRnm (p113, YD, wf, 1, 1, 3)
      call VMATRnm (p121, YD, wf, 1, 2, 1)
      call VMATRnm (p122, YD, wf, 1, 2, 2)
      call VMATRnm (p131, YD, wf, 1, 3, 1)
      call VMATRnm (p222, YD, wf, 2, 2, 2)
      call VMATRnm (p322, YD, wf, 3, 2, 2)
      
      call VMATRnm (p006, YD, wf, 0, 0, 6)
      call VMATRnm (p015, YD, wf, 0, 1, 5)
      call VMATRnm (p024, YD, wf, 0, 2, 4)
      call VMATRnm (p033, YD, wf, 0, 3, 3)
      call VMATRnm (p042, YD, wf, 0, 4, 2)
      call VMATRnm (p051, YD, wf, 0, 5, 1)
      call VMATRnm (p060, YD, wf, 0, 6, 0)
      
      call VMATRnm (p114, YD, wf, 1, 1, 4)
      call VMATRnm (p132, YD, wf, 1, 3, 2)
      call VMATRnm (p133, YD, wf, 1, 3, 3)
      call VMATRnm (p123, YD, wf, 1, 2, 3)
      call VMATRnm (p124, YD, wf, 1, 2, 4)
      call VMATRnm (p141, YD, wf, 1, 4, 1)
      call VMATRnm (p142, YD, wf, 1, 4, 2)
      
      call VMATRnm (p233, YD, wf, 2, 3, 3)
      call VMATRnm (p224, YD, wf, 2, 2, 4)
      call VMATRnm (p242, YD, wf, 2, 4, 2)
      
      call VMATRnm (p333, YD, wf, 3, 3, 3)
      call VMATRnm (p324, YD, wf, 3, 2, 4)
      call VMATRnm (p342, YD, wf, 3, 4, 2)
      
      call VMATRnm (p433, YD, wf, 4, 3, 3)
      call VMATRnm (p424, YD, wf, 4, 2, 4)
      call VMATRnm (p442, YD, wf, 4, 4, 2)

      call VMATRnm (p533, YD, wf, 5, 3, 3)
         
      IF(L .eq. 0) then

         D1  = U00 + a1 * (p020 + p002 - p111*2/3)
     :               + a5 * (U20 + U02)            
     : + a8 * (p040 + p004 + p022*2 - p222*4/5)
     : + a13* (U40 + U22*10/3 + U04)
     : + a14 * (p060 + 5*p042 + 5*p024 + p006
     :       -   2*p133 - 6*p333/7)
      ELSE
         call VMATRnm (s11, YD, wf, L+1, 1, 1)
         call VMATRnm (s12, YD, wf, L+1, 1, 2)
         call VMATRnm (s21, YD, wf, L+1, 2, 1)
         call VMATRnm (s13, YD, wf, L+1, 1, 3)
         call VMATRnm (s31, YD, wf, L+1, 3, 1)
         call VMATRnm (s22, YD, wf, L+1, 2, 2)
         call VMATRnm (t11, YD, wf, L-1, 1, 1)
         call VMATRnm (t12, YD, wf, L-1, 1, 2)
         call VMATRnm (t21, YD, wf, L-1, 2, 1)
         call VMATRnm (t13, YD, wf, L-1, 1, 3)
         call VMATRnm (t31, YD, wf, L-1, 3, 1)
         call VMATRnm (t22, YD, wf, L-1, 2, 2)

         call VMATRnm (P, YD, wf, L+2, 2, 2)
         call VMATRnm (Q, YD, wf, L  , 2, 2)
         call VMATRnm (R, YD, wf, L-2, 2, 2)
         
         D1 = a1 * (s11/(2*L+3) - t11/(2*L-1)) 
         B8= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2
     :                                 + R/(2*L-1)/(2*L-3))

         call VMATRnm (P, YD, wf, L+3, 3, 3)
         call VMATRnm (Q, YD, wf, L+1, 3, 3)
         call VMATRnm (R, YD, wf, L-1, 3, 3)
         call VMATRnm (S, YD, wf, L-3, 3, 3)
                                    
         B14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 

         IF(L .eq. 1) THEN
            D1 =D1 - a5 * U11*2
     :             + a13 *4*(-U31 - U13)
            B8=    - a8  *(p031 + p013 - p122*3/5 - p322/35)*3 
            B14=     a14 *(-5*p051 - 10*p033 - 5*p015
     :                 +  9*p142/5 + 3*p342/35 +  9*p124/5 + 3*p324/35
     :                   -  4*p033/5 + 64*p233/35 - 4*p433/105)
         END IF

         IF(L .eq. 2) THEN
            D1  = D1
     :         +     a13 * 8/3 *  U22
            B14  =   a14 *(5*p042 - 2*p242/7 + p442/21
     :                    +5*p024 - 2*p224/7 + p424/21
     :                  - 18*p133/7 - 2*p533/77)
        END IF
        D1  = D1 + B8 + B14
      END IF
         
C Second term D2: L=Lg->Lf transition

      L = Lg
      DO i=1,meshr
         RT=rmesh(i,1)
         x = QQ*RT
         BE= BES(J,x)
         if(J.eq.0)BE= BES(J,x)-1.0
         YD(i) = wf(i)*BE
      end do
      
      call UMATRnm (U00,  wg, YD, 0, 0)
      call UMATRnm (U20,  wg, YD, 2, 0)
      call UMATRnm (U02,  wg, YD, 0, 2)
      call UMATRnm (U40,  wg, YD, 4, 0)
      call UMATRnm (U22,  wg, YD, 2, 2)
      call UMATRnm (U04,  wg, YD, 0, 4)
      call UMATRnm (U11,  wg, YD, 1, 1)
      call UMATRnm (U31,  wg, YD, 3, 1)
      call UMATRnm (U13,  wg, YD, 1, 3)
 
      call VMATRnm (p002, wg, YD, 0, 0, 2)
      call VMATRnm (p003, wg, YD, 0, 0, 3)
      call VMATRnm (p004, wg, YD, 0, 0, 4)
      call VMATRnm (p005, wg, YD, 0, 0, 5)
      call VMATRnm (p012, wg, YD, 0, 1, 2)
      call VMATRnm (p013, wg, YD, 0, 1, 3)
      call VMATRnm (p014, wg, YD, 0, 1, 4)
      call VMATRnm (p020, wg, YD, 0, 2, 0)
      call VMATRnm (p021, wg, YD, 0, 2, 1)
      call VMATRnm (p022, wg, YD, 0, 2, 2)
      call VMATRnm (p030, wg, YD, 0, 3, 0)
      call VMATRnm (p031, wg, YD, 0, 3, 1)
      call VMATRnm (p040, wg, YD, 0, 4, 0)
      call VMATRnm (p041, wg, YD, 0, 4, 1)
      call VMATRnm (p050, wg, YD, 0, 5, 0)
      call VMATRnm (p111, wg, YD, 1, 1, 1)
      call VMATRnm (p112, wg, YD, 1, 1, 2)
      call VMATRnm (p113, wg, YD, 1, 1, 3)
      call VMATRnm (p121, wg, YD, 1, 2, 1)
      call VMATRnm (p122, wg, YD, 1, 2, 2)
      call VMATRnm (p131, wg, YD, 1, 3, 1)
      call VMATRnm (p222, wg, YD, 2, 2, 2)
      call VMATRnm (p322, wg, YD, 3, 2, 2)
      
      call VMATRnm (p006, wg, YD, 0, 0, 6)
      call VMATRnm (p015, wg, YD, 0, 1, 5)
      call VMATRnm (p024, wg, YD, 0, 2, 4)
      call VMATRnm (p033, wg, YD, 0, 3, 3)
      call VMATRnm (p042, wg, YD, 0, 4, 2)
      call VMATRnm (p051, wg, YD, 0, 5, 1)
      call VMATRnm (p060, wg, YD, 0, 6, 0)
      
      call VMATRnm (p114, wg, YD, 1, 1, 4)
      call VMATRnm (p132, wg, YD, 1, 3, 2)
      call VMATRnm (p133, wg, YD, 1, 3, 3)
      call VMATRnm (p123, wg, YD, 1, 2, 3)
      call VMATRnm (p124, wg, YD, 1, 2, 4)
      call VMATRnm (p141, wg, YD, 1, 4, 1)
      call VMATRnm (p142, wg, YD, 1, 4, 2)
      
      call VMATRnm (p233, wg, YD, 2, 3, 3)
      call VMATRnm (p224, wg, YD, 2, 2, 4)
      call VMATRnm (p242, wg, YD, 2, 4, 2)
      
      call VMATRnm (p333, wg, YD, 3, 3, 3)
      call VMATRnm (p324, wg, YD, 3, 2, 4)
      call VMATRnm (p342, wg, YD, 3, 4, 2)
      
      call VMATRnm (p433, wg, YD, 4, 3, 3)
      call VMATRnm (p424, wg, YD, 4, 2, 4)
      call VMATRnm (p442, wg, YD, 4, 4, 2)
      
      call VMATRnm (p533, wg, YD, 5, 3, 3)

      IF(L .eq. 0) then

         D2  = U00 + a1 * (p020 + p002 - p111*2/3)
     :               + a5 * (U20 + U02)            
     : + a8 * (p040 + p004 + p022*2 - p222*4/5)
     : + a13* (U40 + U22*10/3 + U04)
     : + a14 * (p060 + 5*p042 + 5*p024 + p006
     :       -   2*p133 - 6*p333/7)
      ELSE
         call VMATRnm (s11, wg, YD, L+1, 1, 1)
         call VMATRnm (s12, wg, YD, L+1, 1, 2)
         call VMATRnm (s21, wg, YD, L+1, 2, 1)
         call VMATRnm (s13, wg, YD, L+1, 1, 3)
         call VMATRnm (s31, wg, YD, L+1, 3, 1)
         call VMATRnm (s22, wg, YD, L+1, 2, 2)
         call VMATRnm (t11, wg, YD, L-1, 1, 1)
         call VMATRnm (t12, wg, YD, L-1, 1, 2)
         call VMATRnm (t21, wg, YD, L-1, 2, 1)
         call VMATRnm (t13, wg, YD, L-1, 1, 3)
         call VMATRnm (t31, wg, YD, L-1, 3, 1)
         call VMATRnm (t22, wg, YD, L-1, 2, 2)

         call VMATRnm (P, wg, YD, L+2, 2, 2)
         call VMATRnm (Q, wg, YD, L  , 2, 2)
         call VMATRnm (R, wg, YD, L-2, 2, 2)
         
         D2 = a1 * (s11/(2*L+3) - t11/(2*L-1)) 
         C8= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2
     :                                 + R/(2*L-1)/(2*L-3))
         call VMATRnm (P, wg, YD, L+3, 3, 3)
         call VMATRnm (Q, wg, YD, L+1, 3, 3)
         call VMATRnm (R, wg, YD, L-1, 3, 3)
         call VMATRnm (S, wg, YD, L-3, 3, 3)
                                    
         C14 = a14 * 15*(  P/(( 3+2*L)*( 5+2*L)*( 7+2*L)) -
     :                   3*Q/((-1+2*L)*( 3+2*L)*( 5+2*L)) +
     :                   3*R/((-3+2*L)*(-1+2*L)*( 3+2*L)) -
     :                     S/((-5+2*L)*(-3+2*L)*(-1+2*L))) 
         IF(L .eq. 1) THEN
            D2  = D2 - a5 * U11*2
     :  +  a13 * 4 * (-U31 - U13)
           C8=     - a8 *(p031 + p013 - p122*3/5 - p322/35)*3 
           C14=     a14 *(-5*p051 - 10*p033 - 5*p015
     :                  +  9*p142/5 + 3*p342/35 +  9*p124/5 + 3*p324/35
     :                  -  4*p033/5 + 64*p233/35 - 4*p433/105)
         END IF
         IF(L .eq. 2) THEN
            D2  = D2
     :         +     a13 * 8/3 *  U22
            C14  =   a14 *(5*p042 - 2*p242/7 + p442/21
     :                    +5*p024 - 2*p224/7 + p424/21
     :                  - 18*p133/7 - 2*p533/77)
         END IF
         D2  = D2 + C8 + C14
      END IF
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * ang

      RETURN
      END

*********************************************************
      
      subroutine pBORN(w,w2,l,l1,J,qq,r1,v1,r,v) 
      include 'par.f'                                            
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
      DIMENSION   w(maxr),w1(maxr), w2(maxr)           
      
* Get rid of the Simpson's weights for continuum WF
* But don't corrupt them for later use. That's why w2, not w1

!      do i=1,maxr
!         if(rmesh(i,3) .ne. 0d0) w1(i) = w2(i)/rmesh(i,3)
!      end do
      do i=1,meshr
         w1(i) = w2(i)/rmesh(i,3)
      end do

      call pbLENGTH  (w,w1, l,l1, j,qq,rmatr)           
*      if(j.eq.1) call xLENGTH  (w,w1, l,l1, rmatr1)
      r =rmatr
*      r1 =rmatr1*qq/3
      r1 = 0.
      v = 0.
      v1 = 0.
      RETURN
      END





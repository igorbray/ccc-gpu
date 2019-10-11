      subroutine xACCEL(wg,wf,Lg,Lf,result)

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
     :                        a10, a11, a12, a13
      if(iabs(Lf-Lg).ne.1) STOP 'WRONG Lf,Lg' 
      Lm=Lg                                   
      if(Lf.gt.Lg) Lm=Lf                       
      const=float((-1)**(Lm))*sqrt(float(Lm)) 
      const = const*2

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

      IF(Lf .eq. 0 .or. Lg .eq. 0) then

         call MATRnm (p0_22,wg, wf, 0,-2, 2)
         call MATRnm (p0_23,wg, wf, 0,-2, 3)
         call MATRnm (p0_24,wg, wf, 0,-2, 4)
         call MATRnm (p0_12,wg, wf, 0,-1, 2)
         call MATRnm (p0_13,wg, wf, 0,-1, 3)
         call MATRnm (p000, wg, wf, 0, 0, 0)
         call MATRnm (p001, wg, wf, 0, 0, 1)
         call MATRnm (p002, wg, wf, 0, 0, 2)
         call MATRnm (p010, wg, wf, 0, 1, 0)
         call MATRnm (p011, wg, wf, 0, 1, 1)
         call MATRnm (p02_2,wg, wf, 0, 2,-2)
         call MATRnm (p02_1,wg, wf, 0, 2,-1)
         call MATRnm (p020, wg, wf, 0, 2, 0)
         call MATRnm (p03_2,wg, wf, 0, 3,-2)
         call MATRnm (p03_1,wg, wf, 0, 3,-1)
         call MATRnm (p030, wg, wf, 0, 3, 0)
         call MATRnm (p04_2,wg, wf, 0, 4,-2)
         call MATRnm (p1_11,wg, wf, 1,-1, 1)
         call MATRnm (p1_12,wg, wf, 1,-1, 2)
         call MATRnm (p1_13,wg, wf, 1,-1, 3)
         call MATRnm (p101, wg, wf, 1, 0, 1)
         call MATRnm (p102, wg, wf, 1, 0, 2)
         call MATRnm (p11_1,wg, wf, 1, 1,-1)
         call MATRnm (p110, wg, wf, 1, 1, 0)
         call MATRnm (p111, wg, wf, 1, 1, 1)
         call MATRnm (p12_1,wg, wf, 1, 2,-1)
         call MATRnm (p120, wg, wf, 1, 2, 0)
         call MATRnm (p13_1,wg, wf, 1, 3,-1)
         call MATRnm (p202, wg, wf, 2, 0, 2)
         call MATRnm (p220, wg, wf, 2, 2, 0)

      END IF

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
         D18= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))
         IF(L .eq. 1) THEN
            call MATRnm (p0_13,wg, wf, 0,-1, 3)
            call MATRnm (p011, wg, wf, 0, 1, 1)
            call MATRnm (p102, wg, wf, 1, 0, 2)
            call MATRnm (p302, wg, wf, 3, 0, 2)
            D1  = D1 - a5 * G_1*F1*2
     :               - a9 *(G1*F1 + G_1*F3 - G0*F2*2)*2 
     :          + a12*4*(-F5*G_1 + F4*G0*2 - F3*G1*2 + F2*G2*2 - F1*G3)
     :          + a13*4*(-F3*G_1 - F1*G1)                              
      D18=     - a8 *(p011 + p0_13 - p102*3/5 - p302/35)*3 
         END IF
         D1  = D1 + D18

         IF(L .eq. 2) THEN
            D1  = D1
     :          + a12* 8/3 * (F4*G0 - F3*G1*2 + F2*G2)
     :          + a13* 8/3 *  F2*G0
           END IF
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
         D28= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))
         IF(L .eq. 1) THEN
            call MATRnm (p011,  wg, wf, 0, 1, 1)
            call MATRnm (p03_1, wg, wf, 0, 3,-1)
            call MATRnm (p120,  wg, wf, 1, 2, 0)
            call MATRnm (p320,  wg, wf, 3, 2, 0)
            D2  = D2 - a5 * G1*F_1*2
     :               - a9 *(G3*F_1 + G1*F1 - G2*F0*2)*2 
     :       + a12*4 * (-F3*G1 + F2*G2*2 - F1*G3*2 + F0*G4*2 - F_1*G5)
     :       + a13*4 * (-F1*G1 - F_1*G3)                              
                 D28=     - a8 *(p03_1 + p011 - p120*3/5 - p320/35)*3 
         END IF
         D2  = D2 + D28
         IF(L .eq. 2) THEN
            D2  = D2
     :          + a12* 8/3 * (F2*G2 - F1*G3*2 + F0*G4)
     :          + a13* 8/3 *  F0*G2
         END IF
      END IF
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const


      RETURN
      END



      subroutine xDIPOLE(w,w2,l,l1,r1,v1,x1,r,v,x) 
                                                       
      include 'par.f'                                            
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
      common /schf/ psi1s(maxr)                    !w2-continuum
      common /hyl/  phi1s(maxr),maxphi
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr), f(maxr), f1(maxr)
      
* Get rid of the Simpson's weights for continuum WF
* But don't corrupt them for later use. That's why w2, not w1

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

* SP-channel                                  !     \        
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

      call xLENGTH  (f,f1, 0,1, rmatr)           
      call xVELOCITY(f,f1, 0,1, dmatr)           
      call xACCEL   (f,f1, 0,1, xmatr)           
      goto 115
 111  continue


* PD-channel                                  !     \         
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

      call xLENGTH  (f,f1, 1,2, rmatr)           
      call xVELOCITY(f,f1, 1,2, dmatr)           
      call xACCEL   (f,f1, 1,2, xmatr)           
      goto 115
 112  continue


* DF-channel                                  !     \        
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

      call xLENGTH  (f,f1, 2,3, rmatr)           
      call xVELOCITY(f,f1, 2,3, dmatr)           
      call xACCEL   (f,f1, 2,3, xmatr)           
      goto 115
 113  continue

* FG-channel                                  !     \        
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

      call xLENGTH  (f,f1, 3,4, rmatr)           
      call xVELOCITY(f,f1, 3,4, dmatr)           
      call xACCEL   (f,f1, 3,4, xmatr)           
      goto 115
 114  continue

* GH-channel                                  !     \        
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

      call xLENGTH  (f,f1, 4,5, rmatr)           
      call xVELOCITY(f,f1, 4,5, dmatr)           
      call xACCEL   (f,f1, 4,5, xmatr)           
 115  continue

      r =rmatr
      v =dmatr
      x =xmatr

      RETURN
      END

      subroutine xHYLLERAAS

C  Helium atom ground state               
C  Hylleraas 2nd order wave function.                               

      include 'par.f'
      include 'paratom.f'

     
C     psi1s   - 1s atomic ground states
C     phi1s   - 1s for 3-parameter Hylleraas
     
      common /schf/   psi1s(maxr)
      common /meshrr/ meshr,rmesh(maxr,3)
      common /hyl/    phi1s(maxr),maxphi
      COMMON /hyl1/   z1, an, a1, a2, a3, a4, a5, a6, a7, a8, a9,
     :                        a10, a11, a12, a13
      common/smallr/ formcut,regcut,expcut,fast,match
      logical fast,match
            
      DIMENSION  wd(maxr)
      write(6,100)

*      z1 = 1.816    !Green et. al PR 91, 35 (1953): 3-term Hylleraas
*      an = 1.32135  !
*      a1 = 0.30     !              -Z(r1+r2)                       2  
*      a2 = 0.13     ! F(r r ) = N e          {1 + a  R  + a (r1-r2) } 
*      a3 = 0.       !    1 2                       1  12   2          
*      a4 = 0.
*      a5 = 0.
*      a1 = 0.


*      z1 = 1.81    !Green et. al PR 91, 35 (1953): 6-term Hylleraas
*      an = 1.38189 !
*      a1 = 0.353   !              -Z(r1+r2)                          2   
*      a2 = 0.128   ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)    
*      a3 =-0.101   !    1 2                                 2         2  
*      a4 = 0.033   !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12  }
*      a5 =-0.032   !                                      
*      a6 = 0
*      a7 = 0
*      a8 = 0
*      a9 = 0
      

*      z1 = 1.755012  !Chandrasekhar PR 91, 1172 (1953): 10-term Hylleraas
*      an = 1.359625  !Set 1
*      a1 = 0.350563  !              -Z(r1+r2)                          2   
*      a2 = 0.157394  ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)    
*      a3 =-0.129341  !    1 2                                 2         2  
*      a4 = 0.013019  !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12 
*      a5 =-0.068133  !                                              2
*      a6 = 0.019238  !              + a6 R12(r1+r2) + a7 R12 (r1-r2) 
*      a7 =-0.033843  !                      3          2       2
*      a8 = 0.005575  !              + a8 R12   + a9 R12 (r1-r2)
*      a9 = 0.005342  !                                      
*      a10= 0
*      a11= 0
*      a12= 0
*      a13= 0

      z1 = 1.924965  !Chandrasekhar PR 98, 1050 (1955): 14-term Hylleraas
      an = 1.361717  !Set 1
      a1 = 0.398367  !              -Z(r1+r2)                          2   
      a2 = 0.177426  ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)    
      a3 = 0.011878  !    1 2                                 2         2  
      a4 = 0.020414  !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12 }
      a5 =-0.119940  !                                              2
      a6 = 0.077281  !              + a6 R12(r1+r2) + a7 R12 (r1-r2) 
      a7 =-0.084952  !                      3          2       2
      a8 = 0.022483  !              + a8 R12   + a9 R12 (r1-r2)
      a9 = 0.014528  !                                   2              3
      a10= 0.042902  !              + a10 (r1+r2) (r1-r2)  + a11 (r1+r2)
      a11= 0.001224  !                           2    4             4
      a12=-0.000100  !              + a12 (r1-r2)  R12     + a13 R12
      a13=-0.002061  !

      eps = 1.E-8                                    
      pi =acos(-1.0)

      do i = 1, meshr
         rt = RMESH(i,1)
         phi1s(i) =  rt * exp(-z1*rt) * sqrt(an)
         if(phi1s(i).gt.expcut)  maxphi=i
         wd(i) = phi1s(i) *rt**3 /2
      end do

*      write(6,'(4F9.4)'), a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13



*      print*, '    LENGTH'
*      print*, '    ------'
*      print*, ' '
*      call xLENGTH(phi1s,phi1s,1,0,result)
*      print*, 'RESULT= ',result
*      print*, ' '
*      print*, '    ACCELERATION'
*      print*, '    ------------'
*      print*, ' '
*      call  xACCEL(wd,phi1s,2,3,result) !Test D1
*      call  xACCEL(phi1s,wd,1,0,result) !Test D2
*      print*, 'RESULT= ',result
*
*      STOP!!!!!!

*  Radial integrals

      call RnOVER(0,phi1s,phi1s,P0)
      call RnOVER(1,phi1s,phi1s,P1)
      call RnOVER(2,phi1s,phi1s,P2)
      call RnOVER(3,phi1s,phi1s,P3)
      call RnOVER(4,phi1s,phi1s,P4)

      write(6,111) '<r^0>=', P0, 'Err=', P0-AN*2  /(2*z1)**3
      write(6,111) '<r^1>=', P1, 'Err=', P1-AN*6  /(2*z1)**4
      write(6,111) '<r^2>=', P2, 'Err=', P2-AN*24 /(2*z1)**5
      write(6,111) '<r^3>=', P3, 'Err=', P3-AN*120/(2*z1)**6
      write(6,111) '<r^4>=', P4, 'Err=', P4-AN*720/(2*z1)**7
      write(6,111)
      
 111  format(A,F11.6,3x,A,2x,E12.5)
      
*  Coulomb integrals

      call MATRnm (Q1, phi1s, phi1s, 0, 2, 0)  
      call MATRnm (Q2, phi1s, phi1s, 1, 1, 1)  
      call MATRnm (Q3, phi1s, phi1s, 0, 4, 0)  
      call MATRnm (Q4, phi1s, phi1s, 0, 3, 1)  
      call MATRnm (Q5, phi1s, phi1s, 0, 2, 2)  
      call MATRnm (Q6, phi1s, phi1s, 1, 3, 1)  
      call MATRnm (Q7, phi1s, phi1s, 1, 2, 2)  

      write(6,111) '<F0 r^2 r^0>=',Q1,'Err=', Q1-AN**2*21  /2/(2*z1)**7
      write(6,111) '<F1 r^1 r^1>=',Q2,'Err=', Q2-AN**2*21  /4/(2*z1)**7
      write(6,111) '<F0 r^4 r^0>=',Q3,'Err=', Q3-AN**2*1845/8/(2*z1)**9
      write(6,111) '<F0 r^3 r^1>=',Q4,'Err=', Q4-AN**2*1011/8/(2*z1)**9
      write(6,111) '<F0 r^2 r^2>=',Q5,'Err=', Q5-AN**2*837 /8/(2*z1)**9
      write(6,111) '<F1 r^3 r^1>=',Q6,'Err=', Q6-AN**2*621 /8/(2*z1)**9
      write(6,111) '<F1 r^2 r^2>=',Q7,'Err=', Q7-AN**2*555 /8/(2*z1)**9

      F1 = P0**2 + 4*a2*(P2*P0 - P1**2) 
     :       +  2*a2**2*(P4*P0 - 4* P3*P1 + 3*P2**2)

      F2 = 4*a1*(Q1-Q2*0.33333333) 
     :+ 4*a1*a2*(Q3-Q4*2+Q5 - 2*(Q6-Q7)*0.333333333)

      F3 = 2*a1**2*P0*P2
      F = F1 + F2 + F3
      F = (pi*4)**2 * F

      write(6,'(A,F9.4)') 'Norma ', F

*                          _  -Zr  ___
*  Radial orbitals |1s> = VN e    V4pi  
                                              
      do i = 1, meshr
         phi1s(i) = phi1s(i)*sqrt(pi*4)
      end do

      call RnOVER(0,phi1s,phi1s,over)
      print*, 'Norma 1s    ', over
      
* Overlap with HF ground state

      call RnOVER(0,phi1s,psi1s,over)
      print*, 'Overlap GS  ', over

  20  FORMAT (2  F17.7)
 100  FORMAT (//'++++++++++ Hylleraas 2nd order ++++++++++'//)
*
      END      














      
      subroutine xLENGTH(wg,wf,Lg,Lf,result)

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
     :                        a10, a11, a12, a13
      if(iabs(Lf-Lg).ne.1) STOP 'WRONG Lf,Lg' 
      Lm=Lg                                   
      if(Lf.gt.Lg) Lm=Lf                       
      const=float((-1)**(Lm))*sqrt(float(Lm)) 


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

      IF(Lf .eq. 0 .or. Lg .eq. 0) then

         call MATRnm (p003, wg, wf, 0, 0, 3)
         call MATRnm (p004, wg, wf, 0, 0, 4)
         call MATRnm (p005, wg, wf, 0, 0, 5)
         call MATRnm (p012, wg, wf, 0, 1, 2)
         call MATRnm (p013, wg, wf, 0, 1, 3)
         call MATRnm (p014, wg, wf, 0, 1, 4)
         call MATRnm (p021, wg, wf, 0, 2, 1)
         call MATRnm (p022, wg, wf, 0, 2, 2)
         call MATRnm (p023, wg, wf, 0, 2, 3)
         call MATRnm (p030, wg, wf, 0, 3, 0)
         call MATRnm (p031, wg, wf, 0, 3, 1)
         call MATRnm (p032, wg, wf, 0, 3, 2)
         call MATRnm (p040, wg, wf, 0, 4, 0)
         call MATRnm (p041, wg, wf, 0, 4, 1)
         call MATRnm (p050, wg, wf, 0, 5, 0)
         call MATRnm (p112, wg, wf, 1, 1, 2)
         call MATRnm (p113, wg, wf, 1, 1, 3)
         call MATRnm (p114, wg, wf, 1, 1, 4)
         call MATRnm (p121, wg, wf, 1, 2, 1)
         call MATRnm (p122, wg, wf, 1, 2, 2)
         call MATRnm (p123, wg, wf, 1, 2, 3)
         call MATRnm (p131, wg, wf, 1, 3, 1)
         call MATRnm (p132, wg, wf, 1, 3, 2)
         call MATRnm (p141, wg, wf, 1, 4, 1)
         call MATRnm (p223, wg, wf, 2, 2, 3)
         call MATRnm (p232, wg, wf, 2, 3, 2)
      END IF

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
         D18= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))

         IF(L .eq. 1) THEN
            call MATRnm (p023, wg, wf, 0, 2, 3)
            call MATRnm (p041, wg, wf, 0, 4, 1)
            call MATRnm (p132, wg, wf, 1, 3, 2)
            call MATRnm (p332, wg, wf, 3, 3, 2)
            D1  = D1 - a5 * G2*F1*2
     :               - a9 *(G4*F1 + G2*F3 - G3*F2*2)*2 
     :          + a12*4*(-F5*G2 + F4*G3*2 - F3*G4*2 + F2*G5*2 - F1*G6)
     :          + a13*4*(-F3*G2 - F1*G4)                              
            D18=     - a8 *(p041 + p023 - p132*3/5 - p332/35)*3 
         END IF
         D1  = D1 + D18

         IF(L .eq. 2) THEN
            D1  = D1
     :          + a12* 8/3 * (F4*G3 - 2*F3*G4 + F2*G5)
     :          + a13* 8/3 *  F2*G3
         END IF
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
         D28= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 
     :                                  + R/(2*L-1)/(2*L-3))
         IF(L .eq. 1) THEN
            call MATRnm (p014, wg, wf, 0, 1, 4)
            call MATRnm (p032, wg, wf, 0, 3, 2)
            call MATRnm (p123, wg, wf, 1, 2, 3)
            call MATRnm (p323, wg, wf, 3, 2, 3)
            D2  = D2 - a5 * G1*F2*2
     :               - a9 *(G3*F2 + G1*F4 - G2*F3*2)*2 
     :          + a12*4 * (-F6*G1 + F5*G2*2 - F4*G3*2 + F3*G4*2 - F2*G5)
     :          + a13*4 * (-F4*G1 - F2*G3)                              
            D28=     - a8 *(p032 + p014 - p123*3/5 - p323/35)*3 
         END IF
         D2  = D2 + D28
         IF(L .eq. 2) THEN
            D2  = D2
     :          + a12* 8/3 * (F5*G2 - F4*G3*2 + F3*G4)
     :          + a13* 8/3 *  F3*G2 
         END IF
      END IF

      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const
      RETURN
      END



      subroutine xVELOCITY(wg,wf,lg,lf,result)     
                                             
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
     :                        a10, a11, a12, a13
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
      
      IF(L .eq. 0) then

         call MATRnm (p020, YD, wf, 0, 2, 0)
         call MATRnm (p002, YD, wf, 0, 0, 2)
         call MATRnm (p111, YD, wf, 1, 1, 1)
         call MATRnm (p030, YD, wf, 0, 3, 0)
         call MATRnm (p003, YD, wf, 0, 0, 3)
         call MATRnm (p012, YD, wf, 0, 1, 2)
         call MATRnm (p021, YD, wf, 0, 2, 1)
         call MATRnm (p121, YD, wf, 1, 2, 1)
         call MATRnm (p112, YD, wf, 1, 1, 2)
         call MATRnm (p040, YD, wf, 0, 4, 0)
         call MATRnm (p004, YD, wf, 0, 0, 4)
         call MATRnm (p022, YD, wf, 0, 2, 2)
         call MATRnm (p031, YD, wf, 0, 3, 1)
         call MATRnm (p013, YD, wf, 0, 1, 3)
         call MATRnm (p131, YD, wf, 1, 3, 1)
         call MATRnm (p113, YD, wf, 1, 1, 3)
         call MATRnm (p122, YD, wf, 1, 2, 2)
         call MATRnm (p222, YD, wf, 0, 2, 2)
         
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
         D18= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 +
     :                                  + R/(2*L-1)/(2*L-3))

         IF(L .eq. 1) THEN
            call MATRnm (p013, YD, wf, 0, 1, 3)
            call MATRnm (p031, YD, wf, 0, 3, 1)
            call MATRnm (p122, YD, wf, 1, 2, 2)
            call MATRnm (p322, YD, wf, 3, 2, 2)
            D1  = D1 - a5 * G1*F1*2
     :               - a9 *(G3*F1 + G1*F3 - G2*F2*2)*2 
     :  +  a12 * 4 * (-F5*G1 + F4*G2*2 - F3*G3*2 + F2*G4*2 - F1*G5)
     :  +  a13 * 4 * (-F3*G1 - F1*G3)
            D18=     - a8 *(p031 + p013 - p122*3/5 - p322/35)*3 
         END IF
         D1  = D1 + D18

         IF(L .eq. 2) THEN
            D1  = D1
     :         +     a12 * 8/3 * (F4*G2 - F3*G3*2 + F2*G4)
     :         +     a13 * 8/3 *  F2*G2
         END IF
      END IF
         
C Second term D2: L=Lg->Lf transition

      L = Lg
      if(Lf.gt.Lg) IL= 1                      !Lf=Lg+1
      if(Lf.lt.Lg) IL=-1                      !Lf=Lg-1
      call DERIVATIVE(wf,wd)
      DO i=1,meshr
         YD(i) = wd(i) + wf(i)/rmesh(i,1)*il*Lm
      end do
      
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
 
      IF(L .eq. 0) then

         call MATRnm (p020, wg, YD, 0, 2, 0)
         call MATRnm (p002, wg, YD, 0, 0, 2)
         call MATRnm (p111, wg, YD, 1, 1, 1)
         call MATRnm (p030, wg, YD, 0, 3, 0)
         call MATRnm (p003, wg, YD, 0, 0, 3)
         call MATRnm (p012, wg, YD, 0, 1, 2)
         call MATRnm (p021, wg, YD, 0, 2, 1)
         call MATRnm (p121, wg, YD, 1, 2, 1)
         call MATRnm (p112, wg, YD, 1, 1, 2)
         call MATRnm (p040, wg, YD, 0, 4, 0)
         call MATRnm (p004, wg, YD, 0, 0, 4)
         call MATRnm (p022, wg, YD, 0, 2, 2)
         call MATRnm (p031, wg, YD, 0, 3, 1)
         call MATRnm (p013, wg, YD, 0, 1, 3)
         call MATRnm (p131, wg, YD, 1, 3, 1)
         call MATRnm (p113, wg, YD, 1, 1, 3)
         call MATRnm (p122, wg, YD, 1, 2, 2)
         call MATRnm (p222, wg, YD, 0, 2, 2)
      
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
         D28= a8 * 3*(P/(2*L+3)/(2*L+5) - Q/(2*L-1)/(2*L+3)*2 +
     :                                  + R/(2*L-1)/(2*L-3))

         IF(L .eq. 1) THEN
            call MATRnm (p013, wg, YD, 0, 1, 3)
            call MATRnm (p031, wg, YD, 0, 3, 1)
            call MATRnm (p122, wg, YD, 1, 2, 2)
            call MATRnm (p322, wg, YD, 3, 2, 2)
            D2  = D2 - a5 * G1*F1*2
     :               - a9 *(G3*F1 + G1*F3 - G2*F2*2)*2 
     :  +  a12 * 4 * (-F5*G1 + F4*G2*2 - F3*G3*2 + F2*G4*2 - F1*G5)
     :  +  a13 * 4 * (-F3*G1 - F1*G3)
            D28=     - a8 *(p031 + p013 - p122*3/5 - p322/35)*3 
         END IF
         D2  = D2 + D28
         IF(L .eq. 2) THEN
            D2  = D2
     :         +     a12 * 8/3 * (F4*G2 - F3*G3*2 + F2*G4)
     :         +     a13 * 8/3 *  F2*G2
         END IF
      END IF
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const

      RETURN
      END




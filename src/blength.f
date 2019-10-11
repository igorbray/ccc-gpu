      subroutine bLENGTH(wg,wf,lg,lf,j,qq,result)     
                                             
C  Calculates Born ME in length form     
C  with Hylleraas GS wavefunction for He     
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
      common /hyl/    phi1s(maxr),maxphi
      COMMON /hyl1/   z1, an, a1, a2, a3, a4, a5, a6, a7, a8, a9, 
     :           a10, a11, a12, a13,a14, a15, a16, a17, a18, a19
      real*8 x,ang,bes
      
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
      DO i=1,meshr
         RT=rmesh(i,1)
         x = QQ*RT
         BE= BES(J,x)
         if(J.eq.0)BE= BES(J,x)-1.0
         YD(i) = wf(i)*BE
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
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * ang

      RETURN
      END

*********************************************************
      
      subroutine xBORN(w,w2,l,l1,J,qq,r1,v1,r,v) 
      include 'par.f'                                            
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
      DIMENSION   w(maxr),w1(maxr), w2(maxr)           
      
* Get rid of the Simpson's weights for continuum WF
* But don't corrupt them for later use. That's why w2, not w1

!      do i=1,meshr!was maxr
!         if(rmesh(i,3) .ne. 0d0) w1(i) = w2(i)/rmesh(i,3)
!      end do
      do i=1,meshr
         w1(i) = w2(i)/rmesh(i,3)
      end do

      call bLENGTH  (w,w1, l,l1, j,qq,rmatr)           
*      if(j.eq.1) call xLENGTH  (w,w1, l,l1, rmatr1)
      r =rmatr
*      r1 =rmatr1*qq/3
      r1 = 0.
      v = 0.
      v1 = 0.
      RETURN
      END





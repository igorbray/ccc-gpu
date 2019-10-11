













      
      subroutine hLENGTH(wg,wf,Lg,Lf,result)

C  Calculates dipole ME in length form     
C  from 6-term Hylleraas GS wavefunction for He     
C   r                               
C  D =<gf |r1 V| F > +  < gf |r2 V| F >
C                 L                  L 
C  Here F  is L-pole component of the Hylleraas      
C        L                                        

      include 'par.f'
      DIMENSION  wg(maxr), wf(maxr)
      common /hyl/    phi1s(maxr),maxphi
      COMMON /hyl1/   z1, an, a1, a2, a3, a4, a5
      
      if(iabs(Lf-Lg).ne.1) STOP 'WRONG Lf,Lg'
      Lm=Lg                                   !Lmax=max(Lg,Lf) 
      if(Lf.gt.Lg) Lm=Lf                      !        Lmax  ____ 
      const=float((-1)**(Lm))*sqrt(float(Lm)) !const=(-1)   VLmax

C First term D1: L=Lf->Lg transition

      L = Lf
      IF(L .eq. 0) then
         call RnOVER(1, phi1s, wg, G1) 
         call RnOVER(2, phi1s, wg, G2) 
         call RnOVER(3, phi1s, wg, G3) 
         call RnOVER(0, phi1s, wf, F0) 
         call RnOVER(1, phi1s, wf, F1) 
         call RnOVER(2, phi1s, wf, F2) 
         call MATRnm (V1, wg, wf, 0, 3, 0)
         call MATRnm (V2, wg, wf, 0, 1, 2)
         call MATRnm (V3, wg, wf, 1, 2, 1)
         D1  = G1*F0 + a1 * (V1 + V2 - V3*2/3)
     :               + a2 * (G3*F0 + G1*F2 - G2*F1*2)
     :               + a4 * (G3*F0 + G1*F2 + G2*F1*2)
     :               + a5 * (G3*F0 + G1*F2)          
     :               + a3 * (G2*F0 + G1*F1)          

      ELSE
         call MATRnm (V1, wg, wf, L+1, 2, 1)
         call MATRnm (V2, wg, wf, L-1, 2, 1)
         D1 = a1 * (V1/(2*L+3) - V2/(2*L-1)) 
         IF(L .eq. 1) THEN
            call RnOVER(2, phi1s, wg, G2) 
            call RnOVER(1, phi1s, wf, F1) 
            D1  = D1 - a5 * G2*F1*2
         END IF
      END IF
         
C Second term D2: L=Lg->Lf transition

      L = Lg
      IF(L .eq. 0) then
         call RnOVER(0, phi1s, wg, G0) 
         call RnOVER(1, phi1s, wg, G1) 
         call RnOVER(2, phi1s, wg, G2) 
         call RnOVER(1, phi1s, wf, F1) 
         call RnOVER(2, phi1s, wf, F2) 
         call RnOVER(3, phi1s, wf, F3) 
         call MATRnm (V1, wg, wf, 0, 2, 1)
         call MATRnm (V2, wg, wf, 0, 0, 3)
         call MATRnm (V3, wg, wf, 1, 1, 2)
         D2 = F1*G0 + a1 * (V1 + V2 - V3*2/3)
     :              + a2 * (G0*F3 + G2*F1 - G1*F2*2)   
     :              + a4 * (G0*F3 + G2*F1 + G1*F2*2)  
     :              + a5 * (G0*F3 + G2*F1)            
     :              + a3 * (G0*F2 + G1*F1)            

      ELSE 
         call MATRnm (V1, wg, wf, L+1, 1, 2)
         call MATRnm (V2, wg, wf, L-1, 1, 2)
         D2 = a1 * (V1/(2*L+3) - V2/(2*L-1)) 
         IF(L .eq. 1) THEN
            call RnOVER(1, phi1s, wg, G1) 
            call RnOVER(2, phi1s, wf, F2) 
            D2  = D2 - a5 * G1*F2*2
         END IF
      END IF
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const

      RETURN
      END

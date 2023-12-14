      subroutine xHYLLERAAS

C  Helium atom ground state               
C  Hylleraas 2nd order wave function.                               

      include 'par.f'
      include 'paratom.f'

      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
     
C     psi1s   - 1s atomic ground states
C     phi1s   - 1s for 3-parameter Hylleraas
     
      common /schf/   psi1s(maxr)
      common /meshrr/ meshr,rmesh(maxr,3)
      common /hyl/    phi1s(maxr),maxphi
      COMMON /hyl1/   z1, an, a1, a2, a3, a4, a5, a6, a7, a8, a9, 
     :        a10, a11, a12, a13,a14, a15, a16, a17, a18, a19
      COMMON /hyl2/   b6,b7,b8  ! For hydrogen-11 only
      common/smallr/ formcut,regcut,expcut,fast,match
      logical fast,match
            
      DIMENSION  wd(maxr), wt(maxr)

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

      if(ntype.ne.0)Stop "This is not a Hylleraas ground state" 

C Generally a 20-term Hylleraas is used with ntype=0
C For inferior Hylleraas uncomment a reduced set of parameters

      if(nznuc.eq.1)then
         z1  =   0.675000       !Hart & Herzberg PR 106, 79 (1957): H- 20-term Hylleraas  
         aN  =   0.071640       !                                                      
         a1  =   0.337294       !              -Z(r1+r2)                          2    
         a2  =   0.080883       ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)     
         a3  =  -0.213129       !    1 2                                 2         2   
         a4  =   0.020038       !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12    
         a5  =  -0.028716       !                                              2       
         a6  =  -0.015438       !              + a6 R12(r1+r2) + a7 R12 (r1-r2)        
         a7  =  -9.218966d-3    !                      3          2       2            
         a8  =   4.329046d-3    !              + a8 R12   + a9 R12 (r1-r2)             
         a9  =   7.869764d-4    !                                   2              3   
         a10 =  -1.775586d-3    !              + a10 (r1+r2) (r1-r2)  + a11 (r1+r2)    
         a11 =  -7.408412d-4    !                           2    4             4       
         a12 =   1.630583d-6    !              + a12 (r1-r2)  R12     + a13 R12        
         a13 =  -2.731061d-4    !                                                      
         a14 =   6.274094d-6    !              + a_{14}u^5                             
         a15 =  -6.382931d-5    !              + a_{15}t^2u^3                          
         a16 =  -1.844231d-4    !              + a_{16}s^2t^2                          
         a17 =   1.558570d-5    !              + a_{17}s^4                             
         a18 =   6.483506d-4    !              + a_{18}st^2u                           
         a19 =   6.886024d-4    !              + a_{19}t^4                             
         Etot=   0.527644

*         z1 = 0.707735          !Henrich, Astrophys. J 99, 59 (1944): H- 11-term Hylleraas
*         an = 0.0728887         !
*         a1 = 0.276778          !              -Z(r1+r2)                          2   
*         a2 = 0.0645509         ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)    
*         a3 =-0.157708          !    1 2                                 2         2  
*         a4 = 0.0073342         !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12 
*         a5 =-0.0109206         !              4            6            4    2
*         b6 = 4.34469d-4        !  + b6 (r1-r2) + b7 (r1-r2)  + b8 (r1-r2) R12
*         b7 = 4.33717d-6        !          2       2              2    4          
*         b8 = 0.73220d-6        !  + a9 R12 (r1-r2)  + a12 (r1-r2)  R12           
*         a9 = 0.36386d-4        !                     
*         a12= 1.10006d-6        !
*         Etot=   0.527559
          
*         z1 = 0.7675            !Bethe Z. Phys. 57, 815 (1929): H- 3-term Hylleraas
*         an = 0.0617            !Bethe, Salpeter Quantum mechanics of 
*         a1 = 0.307             !  1 and 2 electron atoms, 1977, p.155
*         a2 = 0.117             !             
*         Etot=   0.5253
        
      end if

      if(nznuc.eq.2)then
         z1 = 1.935000          !Hart & Herzberg PR 106, 79 (1957): He 20-term Hylleraas 
         an = 1.366869          !                                                        
         a1 = 0.417070          !              -Z(r1+r2)                          2      
         a2 = 0.208366          ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)       
         a3 =-0.005234          !    1 2                                 2         2     
         a4 = 0.049090          !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12      
         a5 =-0.157224          !                                              2         
         a6 = 0.082327          !              + a6 R12(r1+r2) + a7 R12 (r1-r2)          
         a7 =-0.125878          !                      3          2       2              
         a8 = 0.047677          !              + a8 R12   + a9 R12 (r1-r2)               
         a9 = 0.040706          !                                   2              3     
         a10= 0.031908          !              + a10 (r1+r2) (r1-r2)  + a11 (r1+r2)      
         a11=-0.009497          !                           2    4             4         
         a12= 0.000512          !              + a12 (r1-r2)  R12     + a13 R12          
         a13=-0.009229          !                                                      
         a14= 0.000706          !              + a_{14}u^5                             
         a15=-0.007251          !              + a_{15}t^2u^3                          
         a16= 0.002303          !              + a_{16}s^2t^2                          
         a17= 0.001366          !              + a_{17}s^4                             
         a18= 0.001605          !              + a_{18}st^2u                           
         a19= 0.000758          !              + a_{19}t^4                             
         Etot=   2.903717

*         z1 = 1.924980          !Chandrasekhar PR 98, 1050 (1955): He 18-term Hylleraas 
*         an = 1.350463          !                                                        
*         a1 = 0.413896          !              -Z(r1+r2)                          2      
*         a2 = 0.211971          ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)       
*         a3 = 0.029010          !    1 2                                 2         2     
*         a4 = 0.005039          !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12      
*         a5 =-0.149093          !                                              2         
*         a6 = 0.079148          !              + a6 R12(r1+r2) + a7 R12 (r1-r2)          
*         a7 =-0.125874          !                      3          2       2              
*         a8 = 0.045441          !              + a8 R12   + a9 R12 (r1-r2)               
*         a9 = 0.043516          !                                   2              3     
*         a10= 0.028227          !              + a10 (r1+r2) (r1-r2)  + a11 (r1+r2)      
*         a11= 0.007138          !                           2    4             4         
*         a12= 0.000502          !              + a12 (r1-r2)  R12     + a13 R12          
*         a13=-0.009934          !                                                        
*         a14= 0.000930
*         a15=-0.007526
*         a16= 0.003074
*         a17=-0.000805
*         Etot=   2.903706
          
*         z1 = 1.924965          !Chandrasekhar PR 98, 1050 (1955): 14-term Hylleraas
*         an = 1.361717          !Set 1
*         a1 = 0.398367          !              -Z(r1+r2)                          2   
*         a2 = 0.177426          ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)    
*         a3 = 0.011878          !    1 2                                 2         2  
*         a4 = 0.020414          !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12 
*         a5 =-0.119940          !                                              2
*         a6 = 0.077281          !              + a6 R12(r1+r2) + a7 R12 (r1-r2)       
*         a7 =-0.084952          !                      3          2       2           
*         a8 = 0.022483          !              + a8 R12   + a9 R12 (r1-r2)            
*         a9 = 0.014528          !                                   2              3  
*         a10= 0.042902          !              + a10 (r1+r2) (r1-r2)  + a11 (r1+r2)   
*         a11= 0.001224          !                           2    4             4      
*         a12=-0.000100          !              + a12 (r1-r2)  R12     + a13 R12       
*         a13=-0.002061          !
*         Etot=   2.903700
          
*         z1 = 1.755012          !Chandrasekhar PR 91, 1172 (1953): He 10-term Hylleraas
*         an = 1.359625          !Set 1
*         a1 = 0.350563          !              -Z(r1+r2)                          2   
*         a2 = 0.157394          ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)    
*         a3 =-0.129341          !    1 2                                 2         2  
*         a4 = 0.013019          !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12 
*         a5 =-0.068133          !                                              2
*         a6 = 0.019238          !              + a6 R12(r1+r2) + a7 R12 (r1-r2) 
*         a7 =-0.033843          !                      3          2       2
*         a8 = 0.005575          !              + a8 R12   + a9 R12 (r1-r2)
*         a9 = 0.005342          !                                      
*         Etot=   2.903602
          
*         z1 = 1.81              !Green et. al PR 91, 35 (1953): He 6-term Hylleraas
*         an = 1.38189           !
*         a1 = 0.353             !              -Z(r1+r2)                          2   
*         a2 = 0.128             ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)    
*         a3 =-0.101             !    1 2                                 2         2  
*         a4 = 0.033             !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12  }
*         a5 =-0.032             !                                      
*         Etot=   2.903240

      endif

      if(nznuc.eq.3)then
         z1  =  3.100000        !Hart & Herzberg PR 106, 79 (1957): Li+ 20-term Hylleraas
         an  =  5.770574        !                                                      
         a1  =  0.427894        !              -Z(r1+r2)                          2    
         a2  =  0.320937        ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)     
         a3  =  0.153723        !    1 2                                 2         2   
         a4  =  0.102493        !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12    
         a5  = -0.291352        !                                              2       
         a6  =  0.178913        !              + a6 R12(r1+r2) + a7 R12 (r1-r2)        
         a7  = -0.338922        !                      3          2       2            
         a8  =  0.132405        !              + a8 R12   + a9 R12 (r1-r2)             
         a9  =  0.192820        !                                   2              3   
         a10 =  0.093489        !              + a10 (r1+r2) (r1-r2)  + a11 (r1+r2)    
         a11 = -0.022064        !                           2    4             4       
         a12 =  5.220202d-3     !              + a12 (r1-r2)  R12     + a13 R12        
         a13 = -0.041711        !                                                      
         a14 =  5.236457d-3     !              + a_{14}u^5                             
         a15 = -0.048497        !              + a_{15}t^2u^3                          
         a16 =  0.018400        !              + a_{16}s^2t^2                          
         a17 =  6.955161d-3     !              + a_{17}s^4                             
         a18 = -0.026701        !              + a_{18}st^2u                           
         a19 = -3.860425d-3     !              + a_{19}t^4                             
         Etot=  7.279905


*         z1 = 2.800028          !Chandrasekhar PR 98, 1050 (1955): Li+ 10-term Hylleraas
*         an = 5.750594          !Set 1
*         a1 = 0.359804          !              -Z(r1+r2)                          2   
*         a2 = 0.238563          ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)    
*         a3 =-0.088645          !    1 2                                 2         2  
*         a4 = 0.016782          !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12 
*         a5 =-0.130358          !                                              2
*         a6 = 0.055045          !              + a6 R12(r1+r2) + a7 R12 (r1-r2) 
*         a7 =-0.090609          !                      3          2       2
*         a8 = 0.016068          !              + a8 R12   + a9 R12 (r1-r2)
*         a9 = 0.019682          !                                      
*         Etot=  7.279762

      end if

      if(nznuc.eq.8)then
         z1  =  8.800000        !Hart & Herzberg PR 106, 79 (1957): O6+ 20-term Hylleraas
         an  =  141.2987        !                                                      
         a1  =  0.440607        !              -Z(r1+r2)                          2    
         a2  =  0.868595        ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)     
         a3  =  0.864501        !    1 2                                 2         2   
         a4  =  0.663489        !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12    
         a5  = -0.948511        !                                              2       
         a6  =  0.601062        !              + a6 R12(r1+r2) + a7 R12 (r1-r2)        
         a7  = -2.850145        !                      3          2       2            
         a8  =  1.177743        !              + a8 R12   + a9 R12 (r1-r2)             
         a9  =  4.998733        !                                   2              3   
         a10 =  0.807402        !              + a10 (r1+r2) (r1-r2)  + a11 (r1+r2)    
         a11 = -0.081824        !                           2    4             4       
         a12 =  0.923021        !              + a12 (r1-r2)  R12     + a13 R12        
         a13 = -1.076679        !                                                      
         a14 =  0.386667        !              + a_{14}u^5                             
         a15 = -3.163547        !              + a_{15}t^2u^3                          
         a16 =  0.553211        !              + a_{16}s^2t^2                          
         a17 =  0.227659        !              + a_{17}s^4                             
         a18 = -1.364837        !              + a_{18}st^2u                           
         a19 = -0.185285        !              + a_{19}t^4                             
         Etot=  59.156581
         
*         z1 = 8.200012          !Chandrasekhar PR 98, 1050 (1955): O^6+ 10-term Hylleraas 
*         an = 141.2470          !Set 1                                                    
*         a1 = 0.377172          !              -Z(r1+r2)                          2       
*         a2 = 0.675101          ! F(r r ) = N e          {1 + a1  R12 + a2 (r1-r2)        
*         a3 = 0.292231          !    1 2                                 2         2      
*         a4 = 0.134920          !              + a3 (r1+r2)  + a4 (r1+r2)  + a5 R12       
*         a5 =-0.506193          !                                              2          
*         a6 = 0.331757          !              + a6 R12(r1+r2) + a7 R12 (r1-r2)           
*         a7 =-0.784625          !                      3          2       2               
*         a8 = 0.168222          !              + a8 R12   + a9 R12 (r1-r2)                
*         a9 = 0.467885          !                                                         
*         Etot=  59.156413

      endif          
      eps = 1.E-8                                    
      pi =acos(-1.0)    

      do i = 1, meshr
         rt = RMESH(i,1)
         phi1s(i) =  rt * exp(-z1*rt) * sqrt(an)
         if(phi1s(i).gt.expcut)  maxphi=i
         wd(i) = phi1s(i) *rt**3 /2 !Test acceleration
         wt(i) =-phi1s(i)           !Test velocity  
     :          * (rt/z1+2./z1**2+2./rt/z1**3) 
      end do

      write(6,101) nznuc,  Etot
 101  FORMAT (////12x,' HYLLERAAS GROUND STATE   ',
     :           /12x,'  ----------------------   '/,
     :           /,   ' Nucleus charge           ',I3,
     :           /,   ' Ground state energy, au  ',F10.6,
     :           /,   ' Parameters Z1, N, a1 ... a20'//)

      write(6,'(F11.6)') z1,an,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,
     :         a11, a12, a13,a14, a15, a16, a17, a18, a19

C  Radial integrals

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
      
C  Coulomb integrals

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

C                          _  -Zr  ___
C  Radial orbitals |1s> = VN e    V4pi  
                                              
      do i = 1, meshr
         phi1s(i) = phi1s(i)*sqrt(pi*4)
      end do

      call RnOVER(0,phi1s,phi1s,over)
      print*, 'Norma 1s    ', over
      
C Overlap with HF ground state

      call RnOVER(0,phi1s,psi1s,over)
      print*, 'Overlap GS  ', over

  20  FORMAT (2  F17.7)

      END      

      
      subroutine xDIPOLE(w,w2,l,l1,r1,v1,x1,r,v,x) 
                                                       
      include 'par.f'                                            
      common/meshrr/ meshr,rmesh(maxr,3)           !w- discrete 
      common /schf/ psi1s(maxr)                    !w2-continuum
      common /hyl/  phi1s(maxr),maxphi
           
      DIMENSION   w(maxr),w1(maxr), w2(maxr), f(maxr), f1(maxr)
      
C Get rid of the Simpson's weights for continuum WF
C But don't corrupt them for later use. That's why w2, not w1

!      do i=1,meshr!was maxr
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

      call xLENGTH  (f,f1, 0,1, rmatr)           
      call xVELOCITY(f,f1, 0,1, dmatr)           
      call xACCEL   (f,f1, 0,1, xmatr)           

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

      call xLENGTH  (f,f1, 1,2, rmatr)           
      call xVELOCITY(f,f1, 1,2, dmatr)           
      call xACCEL   (f,f1, 1,2, xmatr)           
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

      call xLENGTH  (f,f1, 2,3, rmatr)           
      call xVELOCITY(f,f1, 2,3, dmatr)           
      call xACCEL   (f,f1, 2,3, xmatr)           
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

      call xLENGTH  (f,f1, 3,4, rmatr)           
      call xVELOCITY(f,f1, 3,4, dmatr)           
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

      call xLENGTH  (f,f1, 4,5, rmatr)           
      call xVELOCITY(f,f1, 4,5, dmatr)           
      call xACCEL   (f,f1, 4,5, xmatr)           
 115  continue

      r =rmatr
      v =dmatr
      x =xmatr

      RETURN
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
     :        a10, a11, a12, a13,a14, a15, a16, a17, a18, a19
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
     : + a14 * (p070 + 5*p052 - 2*p143 - 6*p343/7 + 5*p034 + p016)
     : + a15 * (p070-2*p061+3*p052-4*p043+3*p034-2*p025+p016
     :       -4*p252/5 + 8*p243/5 - 4*p234/5)
     : + a16 * (G5*F0 - 2*G3*F2 + G1*F4)
     : + a17 * (G5*F0 + 4*G4*F1 + 6*G3*F2 + 4*G2*F3 + G1*F4)
     : + a18 * (p060 - p051 - p024 + p015
     :      -2*p151/3 + 2*p142/3 + 2*p133/3 - 2*p124/3)
     : + a19 * (G5*F0 - 4*G4*F1 + 6*G3*F2 - 4*G2*F3 + G1*F4)
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
            B14 = a14 *(5*p052 - 2*p252/7 + p452/21
     :                + 5*p034 - 2*p234/7 + p434/21
     :              -  18*p143/7 - 2*p543/77)
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
     : + a14 * (p061 + 5*p043 - 2*p134 - 6*p334/7 + 5*p025 + p007)
     : + a15 * (p061-2*p052+3*p043-4*p034+3*p025-2*p016+p007
     :       -4*p243/5 + 8*p234/5 - 4*p225/5)
     : + a16 * (G4*F1 - 2*G2*F3 + G0*F5)
     : + a17 * (G4*F1 + 4*G3*F2 + 6*G2*F3 + 4*G1*F4 + G0*F5)
     : + a18 * (p051 - p042 - p015 + p006
     :      -2*p142/3 + 2*p133/3 + 2*p124/3 - 2*p115/3)
     : + a19 * (G4*F1 - 4*G3*F2 + 6*G2*F3 - 4*G1*F4 + G0*F5)
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
            C14 = a14 *(5*p043 - 2*p243/7 + p443/21 
     :               +  5*p025 - 2*p225/7 + p425/21
     :               - 18*p134/7 - 2*p534/77)
         END IF
         D2  = D2 + C8 + C14 + C15 + C18
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
     :        a10, a11, a12, a13,a14, a15, a16, a17, a18, a19
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
      result = (D1/(2*Lf+1) + D2/(2*Lg+1)) * const

      RETURN
      END

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
     :        a10, a11, a12, a13,a14, a15, a16, a17, a18, a19
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym

      if(iabs(Lf-Lg).ne.1) STOP 'WRONG Lf,Lg' 
      Lm=Lg                                   
      if(Lf.gt.Lg) Lm=Lf                       
      const=float((-1)**(Lm))*sqrt(float(Lm)) 

      const = const*nznuc

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
      


      SUBROUTINE MATRnm (result, fu1, fu2, LV, N, M)

C Calculates radial part of the matrix element with    
C Hylleraas ground state 1s^2

C               N  M         2           N        M
C   < f1, f2 ||r  r  V  || 1s   > ~ R  (r f1,1s; r f2,1s)
C               1  2  Lv             Lv
C
C Functions f1,f2 should be free of any Simpsons wheigts

      include 'par.f'

      common /meshrr/ meshr,rmesh(maxr,3)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     :   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      common /hyl/    phi1s(maxr),maxphi
           
      DIMENSION fu1(maxr), fu2(maxr)
      DIMENSION fun(maxr), temp(maxr)

      if(Lv.lt.0)then !No calculation for Lv<0
         result=0.
         return
      end if

      minfun = 1
      maxfun = maxphi
      do i = 1, maxphi
         fun(i) = fu1(i) * phi1s(i) * rmesh(i,3)
     :                              * rmesh(i,1)**N
      end do

      call form(fun,minfun,maxfun,rpow1(1,lv),rpow2(1,lv),
     :      minrp(lv),maxrp(lv),meshr,temp,i1,i2)
              
      result = 0.d0
      DO I = 1, maxphi
         result = result + temp(i) * fu2(i) * phi1s(i) * rmesh(i,3)
     :                                                 * rmesh(i,1)**M
      end do
      
      RETURN
      END

      subroutine RnOVER(N, P, Q, result)

C Radial overlap integral between two bound states
                                                                        
C                   oo                                                   
C                   f  n                                                 
C              I =  | r  dr P(r) * Q(r)                                  
C                   j                                                    
C                   o                                                    
      
      include 'par.f'

      DIMENSION  P(maxr), Q(maxr)
      common/meshrr/ meshr,rmesh(maxr,3)
      common /hyl/    phi1s(maxr),maxphi
      
      result = 0.d0
      do i = 1, maxphi !It is assumed that P = phi1s
         rt = rmesh(i,1)
         sw = rmesh(i,3)  ! Simpson's weights added
         result = result + P(i)*Q(i)*rt**n * sw
      end do 

      RETURN
      END
      
      subroutine xSATEL(nchtop,lg)
      
C  Satellite intensities in the high photon energy limit
C  according to Aberg (1970)

      include 'par.f'
      include 'paratom.f'

      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
     
      common /schf/   psi1s(maxr)
      common /meshrr/ meshr,rmesh(maxr,3)
      common /hyl/    phi1s(maxr),maxphi
      COMMON /hyl1/   z1, an, a1, a2, a3, a4, a5, a6, a7, a8, a9, 
     :        a10, a11, a12, a13,a14, a15, a16, a17, a18, a19

      COMMON /hyl2/   b6,b7,b8    ! For hydrogen-11 only
      character*1 hh(0:3)         !Spectroscopic labels
      data hh /'s','p','d','f'/
            
      dimension psi(maxr)
      DIMENSION  S(nnmax)
      pi = acos(-1.)
      
      write(6,
     : '(//11x,"ASYMPTOTIC SATELLITE INTENSITIES S(ns)/S(1s)x100%"/,
     :     11x,"-------------------------------------------------" )')

* Total cross-section

      call RnOVER( 0,phi1s, phi1s, G0) 
      call RnOVER( 1,phi1s, phi1s, G1) 
      call RnOVER( 2,phi1s, phi1s, G2) 
      call RnOVER( 3,phi1s, phi1s, G3) 
      call RnOVER( 4,phi1s, phi1s, G4) 
      call RnOVER( 5,phi1s, phi1s, G5) 
      call RnOVER( 6,phi1s, phi1s, G6) 
      call RnOVER( 7,phi1s, phi1s, G7) 
      call RnOVER( 8,phi1s, phi1s, G8) 
      call RnOVER( 9,phi1s, phi1s, G9) 
      call RnOVER(10,phi1s, phi1s, G10) 
      call RnOVER(11,phi1s, phi1s, G11) 
      call RnOVER(12,phi1s, phi1s, G12) 

      b = a1+a3                      !Linear
      c = a2+a4+a5+a6                !Quadratic 
      d = a7+a8+a10+a11              !Cubic
      f = a9+a13+a16+a17+a18+a19+b6  !Quadrupole
      g = a14+a15                    !5th
      h = a12+b7+b8                  !6th  

      ST =  G0                       !Mathematica script (rsphy4)
     :   + G1*2*b                                                            
     :   + G2*(b**2 + 2*c)           !x = 1+b*r+c*r^2+d*r^3+f*r^4+g*r^5+h*r^6  
     :   + G3*(2*b*c + 2*d)          !Collect[x^2,r]                     
     :   + G4*(c**2 + 2*b*d + 2*f)   !FortranForm[%]                     
     :   + G5*(2*c*d + 2*b*f + 2*g)        
     :   + G6*(d**2 + 2*c*f + 2*b*g + 2*h) 
     :   + G7*(2*d*f + 2*c*g + 2*b*h)      
     :   + G8*(f**2 + 2*d*g + 2*c*h)       
     :   + G9*(2*f*g + 2*d*h)              
     :   + G10*(g**2 + 2*f*h)              
     :   + G11*2*g*h                       
     :   + G12*h**2
      
* Partial ns cross-sections

      nn=0
      SP = 0.
      do nch = 1, nchtop
         call  getchinfo (nch, nt, lg, psi, maxpsi, ea, la, na, l)

         if(la.eq.0 .and. na.gt.nn .and. ea .lt. 0.)nn=na !Bound s-channels counter
         
c$$$         write(6,1954) nch,Ea,na,hh(la),hh(l)
c$$$ 1954 format(//11x,'CHANNEL ',I2/,
c$$$     >   11x,'Ion energy',F7.3/,
c$$$     >   11x,'Bound electron',I2,A1/,
c$$$     >   11x,'Free  electron k', A1)

         if(la.eq.0 .and. ea .lt. 0.)then !Bound s-states

* Store ns wave functions

c$$            do ii=1,maxpsi
c$$               write(800+na,'(2E13.4)') rmesh(ii,1), psi(ii)
c$$            end do

*  Coulomb integrals

            call RnOVER( 0,phi1s, psi, G0) 
            call RnOVER( 1,phi1s, psi, G1) 
            call RnOVER( 2,phi1s, psi, G2) 
            call RnOVER( 3,phi1s, psi, G3) 
            call RnOVER( 4,phi1s, psi, G4) 
            call RnOVER( 5,phi1s, psi, G5) 
            call RnOVER( 6,phi1s, psi, G6) 
            
            S(na) = G0 + G1*b + G2*c + G3*d + G4*f + G5*g + G6*h
            SP = SP + S(na)**2
         end if
      end do

      write(6,'(/4x,10(I2,A,6x))//') (i,hh(0), i =2,nn)
      write(6,'(/20F9.4)//') (S(i)**2/S(1)**2*100, i=2,nn)
                        
      write(6,'(/3x,A,F9.4,//)') 'Asymptotic ratio S(++)/S(+)x100%',
     : (ST-SP)/ST*100
                         

      END      

            subroutine mSATEL(nchtop,lg)
      
C  Satellite intensities in the high photon energy limit
C  according to Aberg (1970)
      use vmat_module, only: nodeid

      include 'par.f'
      include 'paratom.f'
      PARAMETER (Ncf = 20)

      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
     
      common /schf/   psi1s(maxr)
      common /meshrr/ meshr,rmesh(maxr,3)

      COMMON /cat1/   CI(0:Ncf, Ncf, Ncf)
      COMMON /cat2/   Nc, nmax(0:Ncf), jmax
      COMMON /cat3/   wa(maxr, 0:Ncf, Ncf)
      COMMON /cat5/   maxNR(0:Ncf, Ncf)
      COMMON /gstspin/   iSpin, meta

      dimension rint(Ncf)


      
      character*1 hh(0:3)         !Spectroscopic labels
      data hh /'s','p','d','f'/
            
      dimension Psi(maxr)
      pi = acos(-1.)
      is = (-1)**iSpin

      if (nodeid.eq.1) write(6,
     : '(//21x,"ASYMPTOTIC SATELLITE INTENSITIES S(ns)x100%"/,
     :     21x,"-------------------------------------------" )')


* Initialization

      do i=1,Ncf
         rint(i)=0
      end do
      
* Partial ns cross-sections

      nn=0
      S = 0.
      do nch = 1, nchtop
         call  getchinfo (nch, nt, lg, psi, maxpsi, ea, la, na, l)
         psi(maxpsi+1:meshr)= 0.0
         if(la.eq.0 .and. na.gt.nn .and. ea .lt. 0.)nn=na !Bound s-channels counter
         
!         write(6,1954) nch,Ea,na,hh(la),hh(l)
 1954 format(//11x,'CHANNEL ',I2/,
     >   11x,'Ion energy',F7.3/,
     >   11x,'Bound electron',I2,A1/,
     >   11x,'Free  electron k', A1)

         if(la.eq.0 .and. ea .lt. 0.)then !Bound s-states

*  Coulomb integrals


            ram = 0.
            do n1=1,nmax(0)
               do n2=1,nmax(0)
                  if(CI(0,n1,n2).ne.0.)then
                     call OVER(psi, wa(1,0,n1), ver1)
                     call OVER(psi, wa(1,0,n2), ver2)
                     wr1 = wa(5,0,n1) 
                     wr2 = wa(5,0,n2) 
                     if(n1.eq.n2) then
                        tmp = wr1*ver1
                     else
                        tmp = (wr1*ver2 + wr2*ver1*is)/sqrt(2.)
                     end if
                     ram = ram + tmp*CI(0,n1,n2)
                  end if
                  rint(na) = ram**2
               end do
            end do
            S = S + rint(na) 
         end if
      end do

!      do n = 1, nn
      if(nn.lt.ninc) STOP "Increase number of S-states"
      do n = ninc, nn
         rint(n) = rint(n)/S
      end do


* Total cross-section

      do i=1,meshr
        Psi(i)= 0.
         do n1=1,nmax(0)
            do n2=1,nmax(0)
               if(CI(0,n1,n2).ne.0.)then
                  wr1 = wa(5,0,n1) 
                  wr2 = wa(5,0,n2) 
                  if(n1.eq.n2) then
                     tmp= wa(i,0,n1)*wr1
                  else
                     tmp=(wa(i,0,n1)*wr2+wa(i,0,n2)*wr1*is)/sqrt(2.)
                  end if
                  Psi(i)= Psi(i) + tmp*CI(0,n1,n2)
               end if
            end do
         end do
      end do

      call OVER(Psi, Psi, Sum)
      
      ratio = (sum - S)/S

      
      if(meta.eq.1.and.nodeid.eq.1) then
         write(6,110) 2*iSpin+1
      write(6,'(/4x,10(I2,A,6x),A)//') (i,hh(0), i =1,10),'oo'
      write(6,'(/20F9.4)//') (rint(n),n=1,10),ratio*100
      elseif (nodeid.eq.1) then
!         write(6,111)
         write(6,'(/4x,10(I2,A,6x),A)//') (i,hh(0), i =1,10),'oo'
         write(6,120) (rint(n)*100,n=1,10),  ratio*100
      end if

      if(nn.lt.10.and.nodeid.eq.1) print'(//4x,A/)',
     >   'WARNING: Less than 10 target s-states are used!'

      
!      do i=1,7
!         print'(F9.4)', rint(i)
!      end do
!      print'(F9.4)', ratio*100

 20   FORMAT (2  F17.7)
 100  FORMAT (//'++++++++++ Satell ++++++++++'//)
 110  FORMAT(//'   SATELLITE INTENSITIES IN METASTABLE ATOM',I3,'S'/,    
     :         '   ----------------------------------------------')

 111  FORMAT(//'   SATELLITE INTENSITIES IN GROUND STATE ATOM'/,    
     :         '   ----------------------------------------------')
 120  FORMAT ( 12(F7.3,2x))
      
      
                        
                         

      END      

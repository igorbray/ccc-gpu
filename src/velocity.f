      subroutine VELOCITY(wi,wf,li,lf,result)
                                                 
*   Calculates dipole ME in velocity form:                        
*                        oo                                       
*              Lmax ___  f                                        
*   <f|r|i>=(-1)   VLmax |drP(r)(-v +- Lmax/r)P(r)                 
*                        j   f                i                   
*                        o
*    v-differential operator                                                  
                                             

      include 'par.f'

      common/meshrr/ meshr,rmesh(maxr,3)
      DIMENSION  wi(MAXR), wf(MAXR), wd(MAXR)
      

      lm=li                                       !Lmax=max(Li,Lf) 
      if(lf.gt.li) lm=lf                          !        Lmax  ____ 
      const=dfloat((-1)**(lm))*dsqrt(dfloat(lm))  !const=(-1)   VLmax

      il=-1                                       !Lf=Li-1
      if(lm-li.gt.0) il=1                         !Lf=Li+1

      call  DERIVATIVE(wi,wd)

      result = 0.d0
      do i=1,meshr
         result = result + wf(i)
     >      *(-wd(i)+dble(il*Lm)*
     >      wi(i)/rmesh(i,1))
     >      * rmesh(i,3)                         !Simpson's weights added
      end do 
      result=result*const

      RETURN
      END


      subroutine ACCEL(wi,wf,li,lf,result)             
                                                    
* Calculates dipole ME in acceleration form             
* Z=2 (helium)                                          
* Z=3 (lithium+)                                          
*                        oo                
*             Lmax  ___  f  Z                
* <f|r|i> = (-1)   VLmax | --- dr P(r) * P(r) 
*                        j r*r     f      i   
*                        o                 

      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      
      DIMENSION  wi(maxr), wf(maxr)
      lm=li                                              !Lmax=max(Li,Lf) 
      if(lf.gt.li) lm=lf                                 !        Lmax  ____ 
      const=dfloat((-1)**(lm))*dsqrt(dfloat(lm))         !const=(-1)   VLmax
      
      result = 0.d0

      do i=1,meshr
         result = result + wi(i)*wf(i) / rmesh(i,1) !1/R-added
     >                                 / rmesh(i,1) !2/R-added
     >                                 * nznuc
     >                                 * rmesh(i,3) !Simpson's weights added
      end do                                         

      result=result*const
      RETURN
      END















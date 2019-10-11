      subroutine LENGTH(wi,wf,li,lf,result)             
                                                    
* Calculates dipole ME in length form             
                                          
*                        oo                
*             Lmax  ___  f                 
* <f|r|i> = (-1)   VLmax | rdr P(r) * P(r) 
*                        j      f      i   
*                        o                 

      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      
      DIMENSION  wi(maxr), wf(maxr)
      lm=li                                              !Lmax=max(Li,Lf) 
      if(lf.gt.li) lm=lf                                 !        Lmax  ____ 
      const=dfloat((-1)**(lm))*dsqrt(dfloat(lm))         !const=(-1)   VLmax
      
      result = 0.d0

      do i=1,meshr
         result = result + wi(i)*wf(i) * rmesh(i,1) !R-added
     >                                 * rmesh(i,3) !Simpson's weights added
      end do                                         

      result=result*const
      RETURN
      END















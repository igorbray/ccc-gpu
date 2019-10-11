      subroutine STEP(h, ndouble, nstep, nmax)

* Returns h, ndouble, nstep and nmax      
                                                   
      include 'par.f'                              
      common/meshrr/ meshr,rmesh(maxr,3)
                                                   
      DIMENSION  n(10)                    
      
      h=rmesh(1,2)
      h0=h
      k=0
      do i=1,maxr
*         write(6,'(I5,F9.4)') i, rmesh(i,2)
         if(rmesh(i,2) .ne. h0) then
            k=k+1
            n(k)=i-1
            h0=rmesh(i,2)
            end if
      end do
      ndouble=k-1
      nstep=n(1)
      nmax=n(k)


      do k=1,ndouble
*         write(6,'(2I5)') k, n(k)
      end do
*         write(6,'(F9.4,3I5)') h,ndouble, nstep, nmax

                                             
      RETURN
      END










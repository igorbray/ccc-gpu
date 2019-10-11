      subroutine DERIVATIVE(w,w1)
                                                    
      include 'par.f'                               
                                                    
      dimension  w(MAXR), w1(MAXR)                     
      common/meshrr/ meshr,rmesh(maxr,3)

      common /double/id,jdouble(22)

      do k=1,id-1
         N1 = jdouble(k)+1
         if(k .eq. 1) N1=1
         N2 = jdouble(k+1)
         h=rmesh(N1,2)
         if(k.eq.ID+2) N2=meshr
*         write(6,'(3I5,2F9.4/)') k,N1,N2,h,rmesh(N2,1)
      
      F1=W(N1)                                       
      F2=W(N1+1)                                       
      F3=W(N1+2)                                       
      F4=W(N1+3)                                       
      F5=W(N1+4)                                       
      D=1./(12.*H)                                   
                                                    
      DO J=N1,N2                                                               
         IF(J.EQ.N1) GOTO 104                    !1st point in radial grid
         GOTO 105                                                             
 104     C=-25.*F1+48.*F2-36.*F3+16.*F4-3.*F5                                 
         GOTO 106                                                             
 105     CONTINUE                                                             
         IF(J.EQ.N2) GOTO 107                    !Nth point in radial grid    
         GOTO 108                                                             
 107     C=25.*F5-48.*F4+36.*F3-16.*F2+3.*F1                                  
         GOTO 106                                                             
 108     IF(J.EQ.N1+1) GOTO 109                  !2nd point in radial grid    
         GOTO 110                                                             
 109     C=-3.*F1-10.*F2+18.*F3-6.*F4+F5                                      
         GOTO 106                                                             
 110     IF(J.EQ.N2-1) GOTO 111                  !(N-1)th point in radial grid
         GOTO  112                                                            
 111     C=3*F5+10.*F4-18.*F3+6.*F2-F1                                        
         GOTO 106                                                             
 112     IF(J.EQ.N1+2) GOTO 113                  !3d point in radial grid
         F1=F2                                   !Take next 5 points          
         F2=F3                                                                
         F3=F4                                                                
         F4=F5                                                                
         F5=W(J+2)                                                            
 113     C=F1-8.*(F2-F4)-F5                       !Arbitrary point 
 106     C=C*D                                                   
         W1(J)=C
*         write(6,'(I4,6F9.4)')J, F1,F2,F3,F4,F5,C
      end do                                                     


      end do


      do i=1,meshr
*         write(119, '(2F9.4)') rmesh(i,1), w(i)
*         write(120, '(2F9.4)') rmesh(i,1), w1(i)
      end do

      
      RETURN
      END










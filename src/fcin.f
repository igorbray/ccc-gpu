       open(2,file=os,form='unformatted',access='sequential')
       open(1,file=vo,
     >    form='unformatted',access='sequential')
       open(3,file=du1)
       open(4,file=du2)
       open(66,file=con)
       if (nc.eq.0) then
          inquire(file='fcwf'//ch(lv),exist=there)
          if (there) go to 10
          open(10, file='fcwf'//ch(lv), status='new',form='unformatted')
          open(20, file='fcen'//ch(lv), status='new')
       else
          write(ench,'(1p,"_",e10.4)') q*q
          inquire(file='fcphase'//ch(lv)//ench,exist=there)
          if (there) go to 10
          open(10, file='fcwf'//ch(lv)//ench, form='unformatted')
          open(20, file='fcen'//ch(lv)//ench)
          open(30, file='fcphase'//ch(lv)//ench)
c$$$          open(10, file='fcwf', form='unformatted')
c$$$          open(20, file='fcen')
       end if 
          
c
       
       call fchf(z,gam,dir,r,bet,eps,h,teta,
     * ala,ala1,ala2,ala3,c3,
     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
     * is,in,il,iq,
     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     * ne5,ne5is,kt,ud,ksu1,gamma,r0)
c

       if (nc.eq.0) then
          close(10)
          close(20)
       else 
          phase = exp((0.0,1.0) * r5(ne+4))
          write(30,*) phase
          close(10)             !,status='delete')
          close(20)!,status='delete')
          close(30)!,status='delete')
       endif 
c
c     Find and Store wave functions
c
 10    continue
       if (nc.eq.0) then
          open(10,file='fcwf'//ch(lv),status='old',form='unformatted')
          open(20, file='fcen'//ch(lv), status='old')
       else
          open(10, file='fcwf'//ch(lv)//ench, status='old',
     >       form='unformatted')
          open(20, file='fcen'//ch(lv)//ench, status='old')
          open(30, file='fcphase'//ch(lv)//ench, status='old')
       end if 
          
       if (nb.gt.0.and.nc.gt.0) then
          print *,'Either NB or NC must be zero'
          stop 'Either NB or NC must be zero'
       end if

       if (nc.eq.1) then
          read(10) (R7(IW), r2(IW), IW = 1, KT )
          read(30,*) phase
          close(30)
          close(10)             !,status='delete')
          close(20)!,status='delete')
          iwff = 1
          do while (abs(r2(iwff)).lt.regcut)
             r2(iwff) = 0.0
             iwff = iwff + 1
          end do
C  Normalize the continuum wave function to range between +/- 1.0
c$$$          c = sqrt(acos(-1.0) * q)
C  Normalize the continuum wave function to range between +/- sqrt(2/pi)
          c = sqrt(2.0 * q)
          do i = iwff, kt
             wff(i) = r2(i) * c
          end do 
       else
          do n = 1, nb
c
             nn = in(is) + n - 1
             if (nn.gt.nnmax) then
                print*,'NN > NNMAX',nn,nnmax
                stop 'increase NNMAX'
             end if
          
             read(20, *)  R5(1)
             enpsinb(nn, lv) = R5(1)
c
             read(10) (R7(IW), R2(IW), IW = 1, KT )
             do i = 1, kt
                psinb(i,nn,lv) = r2(i)
             end do 
             istoppsinb(nn,lv) = kt
             do while (abs(psinb(istoppsinb(nn,lv),nn,lv)) .lt. expcut)
                istoppsinb(nn,lv) = istoppsinb(nn,lv) - 1
             end do 
          end do   
 2        FORMAT(  ' ','Energy=',25X,E15.8,
     *       /' ','HOPMA  =' ,25X,E15.8)
          close(10)
          close(20)
       endif 
       close(1)
       close(2)
       close(3)
       close(4)
       close(66)

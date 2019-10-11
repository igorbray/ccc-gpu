C  Compile as f90 alisher.f iqpackd.o
C     main program for routine kgrid
      implicit real*8 (a-h,o-z)
      parameter (kmax=100)
      dimension xk(kmax), wk(kmax)
      xk(:) = 0d0
      wk(:) = 0d0
      print*,'Enter x0, a, Na, b, Nb, p, Np'
      read*, x0, a, Na, b, Nb, p, Np
      nk = Na + Nb + Np
      call kgrid(x0,a,Na,b,Nb,p,Np,nk,xk,wk)
      do i = 1, nk
         print*,i, xk(i), wk(i)
      enddo 
      end
      
C The following routine sets up quadrature points with weights (gridk and
C weightk) for int dk f(k)/(k-x0). There are three intervals: [0,a], [a,b],
C and [c,oo). The number of points in the intervals are Na, Nb, and Np.
C The function f(k) is assumed to fall-off as k**(-p).    
      subroutine kgrid(x0,a,Na,b,Nb,p,Np,nk,gridk,weightk)
      implicit real*8 (a-h,o-z)
      dimension gridk(nk),weightk(nk),xx(nk),ww(nk),sk(0:10),nnk(10),
     >   wf(2*nk),iwf(2*nk)

      npoints = Na
      endk = a
      npoints2 = Nb
      endk2 = b
      nendk = Np
      endp = p
      midnp = -10  ! -ve means not used
      width = 0.1  ! not used when midnp < 0
      e = x0*x0  ! on-shell energy in Ry
      rk = x0
      
C  Define the intervals and the number of points in each interval
      call makeints(mint,sk,nnk,rk,npoints,width,midnp,enk,nendk,endp,
     >   npoints2,endk2)
      print '(''interval   points    start       stop          '//
     >   'test'')'
      enk = endk
 30   sumi=0.0
      sumip = 0.0
      nqk=0
C  Define the intervals and the number of points in each interval
      call makeints(mint,sk,nnk,rk,npoints,width,midnp,enk,nendk,endp,
     >   npoints2,endk2)
      if (midnp.lt.0) then
         dstop = 0d0
         do j=1,mint-1
            nt=nnk(j)
            nwf=2*nt
            niwf=2*nt
            dstart = dstop
            dstop  = dble(sk(j))
            if (dstart.lt.rk.and.rk.lt.dstop) then
               if (rk.gt.dstop*0.9) then
                  print*,'Increasing dstop:', dstop, dstop * 1.1
                  dstop = dstop * 1.1
               endif 
               if (rk.lt.dstart*1.1) then
                  print*,'Decreasing dstart:', dstart, dstart * 0.9
                  dstart = dstart * 0.9
               endif 
               call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,
     >            0,nwf,wf,niwf,iwf,ier)
               test = 1.0
               startk = dstart
               stopk = dstop
               if (startk.eq.0.0.or.(stopk-rk)/(rk-startk).lt.4.0) then
                  dk = (stopk-startk)/nt/10
                  call getstopk(e,rk,startk,stopk,dk,test,xx,ww,nt)
                  dstop = stopk
                  sk(j) = stopk
               else 
                  dk = startk/nt/10
                  call getstartk(e,rk,startk,stopk,dk,test,xx,ww,nt)
                  dstart = startk
                  sk(j-1) = startk
               endif
            endif
         enddo 
      endif 
c$$$      if (width .le. 0.0.and.rk.lt.2.0*abs(width)) then
c$$$         call improvek(e,sk(1),nnk(1),gfixed,wfixed,endf,enk)
c$$$         sk(1) = endf
c$$$      endif 
      dstop = 0d0
C  Obtain Gaussian knots and weights within each interval
      do j=1,mint-1
         nqk=nqk+nnk(j)
         nt=nnk(j)
         nwf=2*nt
         niwf=2*nt
         if (nt.gt.0) then
            dstart = dstop
            dstop  = dble(sk(j))
            call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,
     >         0,nwf,wf,niwf,iwf,ier)
            do i=nqk-nnk(j)+1,nqk
               gridk(i)=xx(i-nqk+nnk(j))
               weightk(i)=ww(i-nqk+nnk(j))
               if (gridk(i).eq.0.0) weightk(i) = 0.0
c$$$               ecmn = gridk(i) ** 2
c$$$               eta = 0.0
c$$$
c$$$               if (la.eq.li) then
c$$$                  call regular(la,ecmn,eta,ucentr,cntfug(1,la),-1,
c$$$     >               rmesh(1,1),meshr,jdouble,njdouble,regcut,expcut,
c$$$     >               reg,jstart,jstop,phase,sigc)
c$$$                  call getprod(psi,maxpsi,jstart,reg,la,rmesh,
c$$$     >               meshr,ntmp,tmp)
c$$$                  sumi = sumi + tmp * weightk(i) * 2.0 / pi
c$$$               endif 
            end do
c$$$            tmpold = tmp
            if (dstart.lt.rk.and.rk.lt.dstop) then
               print '(i5,i10,f12.6,'' *'',f10.6,f13.5)', j,nnk(j),
     >            dstart, dstop, sumi - sumip
            else 
               print '(i5,i10,2f12.6,f13.5)', j,nnk(j),dstart,dstop,
     >            sumi - sumip
            endif 
            sumip = sumi
         end if 
      end do 
C  Here we have the last interval
      nt=nnk(mint)
      if (nt.gt.0) then
         if (dstop.lt.rk) then
            print*,'The last interval must start > than RK'
            stop 'The last interval must start > than RK'
         end if  
         nqk=nqk+nt
         p=dble(sk(mint))
         dstart=0d0
         dstopp = dstop
         dstop=dstop**(1d0-p)
         nwf=2*nt
         niwf=2*nt
         call cgqf(nt,xx,ww,1,0d0,0d0,dstart,dstop,0,nwf,wf,niwf,iwf,
     >      ier)
         if (ier.ne.0) print*,'KGRID IER:',ier
         do j=nt,1,-1
            jj=nqk-j+1
            gridk(jj)=xx(j)**(1d0/(1d0-p))
            weightk(jj)=ww(j)/(p-1d0)*dble(gridk(jj))**p
            ecmn = gridk(jj) ** 2
c$$$            if (la.eq.li) then
c$$$               eta = 0.0
c$$$               call regular(la,ecmn,eta,ucentr,cntfug(1,la),-1,
c$$$     >            rmesh,meshr,jdouble,njdouble,regcut,expcut,
c$$$     >            reg,jstart,jstop,phase,sigc)
c$$$               call getprod(psi,maxpsi,jstart,reg,la,rmesh,
c$$$     >            meshr,ntmp,tmp)
c$$$               if (abs(tmp).gt.abs(tmpold)*1.5) then
c$$$                  endkt = gridk(jj)
c$$$                  tmpold = tmp
c$$$               endif 
c$$$               sumi = sumi + tmp * weightk(jj) * 2.0 / pi
c$$$            endif 
         end do 
         print '(i5,i10,f12.6,''       oo'',f16.5)', mint,nnk(mint),
     >      dstopp, sumi - sumip
         print'(''fall off power:'',f5.1,23x,''= '',f7.5)',sk(mint),
     >      sumi
c$$$         if (abs(sumi - 1.0).gt.1e-2.and.la.eq.li) nbad = nbad + 1
      end if
C  Check that the integration rule will handle the principle value singularity
      j=1
      if (e.gt.0.0) then
         if (nqk.gt.0) then
            sum=0.0
            nt=0
            do while (sk(j).lt.rk-0.01.and.j.lt.mint)
               sum = 2.0 * atanh(sk(j)/rk)/rk
               nt=nt+nnk(j)
               j=j+1
            end do
            j=mint-1
            tsum=0.0
            do while (sk(j).gt.rk+0.01.and.j.ge.1)
               tsum = - 2.0 * acoth(sk(j)/rk) / rk
               j=j-1
            end do
            sum=sum+tsum
            im=nt+1
            do while (gridk(im).le.(sk(j+1)).and.im.lt.nk)
               sum=sum + 2.0*weightk(im)/(e-gridk(im)*gridk(im))
               im=im+1
            end do
         else
            sum = 0.0
         endif 
         print*,'integral test of singularity:',sum
c$$$         print '(''State, NCH, NST, NA, LA, L, K, EA:'',a4,3i3,
c$$$     >      2i2,1p,2e13.4)',chan(ntmp),nch,ntmp,na,la,li,rk,ea
      else
c$$$         if (nqm.gt.1) then
            sum = 0.0
            do i = 1, nnk(1) + nnk(2)
               sum = sum + weightk(i) * 2.0 / (e - gridk(i) ** 2)
            enddo
            print*,'Test of closed channel integral:',- 2.0 / sqrt(-e)*
     >         atan(sk(2)/sqrt(-e)) / sum, sum
c$$$         endif 
c$$$         print'(''State, NCH, NST, NA, LA, L, E:    '',a4,3i3,2i2,
c$$$     >      1p,e13.4,''    closed'')',chan(ntmp),nch,ntmp,na,la,li,e
      end if
      if (nqk.ne.nk) stop 'inconsistency in kgrid'
      return
      end
      
      
      subroutine getstartk(e,rk,startk,stopk,dk,test,xx,ww,nt)
      implicit real*8 (a-h,o-z)
      dimension  xx(nt),ww(nt)

      do while (dk/startk.gt.1d-14)
         startkold = startk
         if (test.gt.0.0) then
            startk = startk - dk
         else
            startk = startk + dk
         endif
c$$$         print*,test,dk,startk,stopk,nmin
         do n = 1, nt
            xx(n) = (xx(n)-startkold)*(stopk-startk)/
     >         (stopk-startkold) + startk
            ww(n) = ww(n) * (stopk-startk)/(stopk-startkold)
            if (xx(n).lt.rk) nmin = n
         enddo 
         testold = test
         test = 2.0/rk*
     >      (atanh(startk/rk)-acoth(stopk/rk))
         do n = 1, nt
            test = 2.0 * ww(n)/(e - xx(n)**2) + test
         enddo 
         if (test*testold.lt.0.0) dk = dk / 2.0
      enddo 
      return
      end
      
      subroutine getstopk(e,rk,startk,stopk,dk,test,xx,ww,nt)
      implicit real*8 (a-h,o-z)
      dimension xx(nt),ww(nt)

      do while (dk/stopk.gt.1d-14)
         stopkold = stopk
         if (test.gt.0.0) then
            stopk = stopk - dk
         else
            stopk = stopk + dk
         endif
c$$$         print*,test,dk,startk,stopk,nmin
         do n = 1, nt
            xx(n) = (xx(n)-startk)*(stopk-startk)/
     >         (stopkold-startk) + startk
            ww(n) = ww(n) * (stopk-startk)/(stopkold-startk)
            if (xx(n).lt.rk) nmin = n
         enddo 
         testold = test
         test = 2.0/rk*
     >      (atanh(startk/rk)-acoth(stopk/rk))
         do n = 1, nt
            test = 2.0 * ww(n)/(e - xx(n)**2) + test
         enddo 
         if (test*testold.lt.0.0) dk = dk / 2.0
      enddo 
      return
      end

      
      subroutine improvek(e,sk1,nt,gridk,weightk,endf,endk)
      implicit real*8 (a-h,o-z)
      dimension  gridk(nt),weightk(nt)

      endf = sk1
      if (e.lt.0.0) return
      rk = sqrt(e)
      if (sk1.lt.rk) then
         endf = (2.0 * rk - endk)
         if (endf.lt.0.0) stop 'Problem in IMPOVEK; ENDK too big'
         return
      endif 
      k = 1
      do while (gridk(k) .lt. rk)
         k = k + 1
      enddo
      call getint(e,sk1,sk1,gridk,weightk,nt,sum)
      if (sum.gt.0.0) then
         sk1max = sk1
         sk1min = sk1 * (9.0 * rk + gridk(k)) / 10.0 / gridk(k)
         call getint(e,sk1,sk1min,gridk,weightk,nt,sum)
      else 
         sk1min = sk1
         sk1max = sk1 * (9.0 * rk + gridk(k-1)) / 10.0 / gridk(k-1)
         call getint(e,sk1,sk1max,gridk,weightk,nt,sum)
      endif 

      ncount = 0
      do while(abs((sk1max-sk1min)/sk1min).gt.1d-14.and.ncount.lt.50)
         ncount = ncount + 1
         endf = (sk1max + sk1min) / 2.0
         call getint(e,sk1,endf,gridk,weightk,nt,sum)
         if (sum.gt.0.0) then
            sk1max = endf
         else
            sk1min = endf
         endif 
         if (ncount.ge.90) print*,sk1min,sk1max, sum
      enddo 
      end

      subroutine getint(e,sk1,sk1p,gridk,weightk,nt,sum)
      implicit real*8 (a-h,o-z)
      dimension  gridk(nt),weightk(nt)
      r = sk1p / sk1
      sum = 0.0
      do i = 1, nt
         sum = sum + 2.0 * weightk(i) * r / (e - (gridk(i) * r)**2)
      end do 
      rk = sqrt(e)
      tsum = - 2.0 * acoth(sk1p/rk) / rk
      sum = sum + tsum
      end
      
C  Define the k-grid integration intervals and the number of points in
C  each interval
C 1,3,   20,   0.8,    16,     2.5,    6,    6.0,   12,   0.09  
C     npoints, endk, npoints2, endk2, nendk, endp, midnp, width 
      subroutine makeints(mint,sk,nk,rk,np,width,midnp,endk,nendk,endp,
     >   np2,endk2)
      implicit real*8 (a-h,o-z)
      dimension sk(0:10), nk(10)
      nmin(npoints) = min(npoints/2,4)
      
      if (midnp.lt.0) then
         if (rk.gt.endk2) stop 'On shell point in third interval'
         nk(1) = np
         sk(1) = endk
         nk(2) = np2
         sk(2) = endk2         
         mint = 3
         nk(mint) = nendk
         sk(mint) = endp
         return
      endif 

      sk(0) = 0.0
      if (abs(width).lt.0.01) then
         print*,'Width must be greater than 0.01'
         stop 'Width must be greater than 0.01'
      endif
      if (endk+2.0*width.ge.endk2) then
         print*,'ENDK2 must be >= ENDK + 2 * WIDTH',endk2,endk+2.*width
         stop 'ENDK2 must be >= ENDK + 2 * WIDTH'
      endif 
      if (rk.lt.0.0.or.(width.lt.0.0.and.rk.lt.2.0*abs(width))) then
C  Here for closed channels
c$$$         nk(1) = midnp
c$$$         sk(1) = 2.0 * abs(width)
c$$$         nk(2) = np
         nk(1) = np
         sk(1) = endk
         nk(2) = np2
         sk(2) = endk2         
         mint = 3
      else if (rk.lt.width) then
C  Here if singularity is going to be in the first interval
         sk(1) = rk * 2.0
         nk(1) = midnp
         sk(2) = endk
         nk(2) = np
         sk(3) = endk2
         nk(3) = np2
         mint = 4
      else if (rk. lt. endk + 1.0 * width) then
C  Here if singularity is going to be in the second interval
         sk(1) = rk - width
         sk(2) = rk + width

C  Use the following to make the k grid symmetric about rk**2
c$$$            sk(1) = sqrt(rk**2 - 2.0 * rk * width)
c$$$            sk(2) = sqrt(rk**2 + 2.0 * rk * width)

         nk(2) = midnp
         if (rk .lt. endk - 1.0 * width ) then
            nk(1) = max(nmin(np),int(np * rk / (endk - 1.0 * width)))
            nk(1) = min(nk(1), np - nmin(np))
            npleft = np - nk(1)
            sk(3) = endk
            nk(3) = npleft
            sk(4) = endk2
            nk(4) = np2
            mint = 5
         else
            nk(1) = np
            sk(3) = endk2
            nk(3) = np2
            mint = 4
         endif 
      else if (rk. lt. endk2) then
C  Here the singularity is in the third interval
         sk(1) = endk
         nk(1) = np
         sk(2) = rk - width
         sk(3) = rk + width
         nk(3) = midnp
C  The two factors of 4.0 below put extra points in the first interval
         nk(2) = max(nmin(np2),int(np2 * (sk(2)-sk(1)) * 4.0 /
     >      (4.0 * (sk(2) - sk(1)) + (endk2-endk))))
c$$$         nk(2) = max(nmin(np2),int(np2 * (rk - endk) / (endk2-endk)))
         nk(2) = min(nk(2), np2 - nmin(np2))
         npleft = np2 - nk(2)
         if (rk .lt. endk2 - width) then
            sk(4) = endk2
            nk(4) = npleft
            mint = 5
         else
            nk(2) = np2
            mint = 4
         endif    
      else
C  This is for very large incident energies
         sk(1) = endk
         nk(1) = np
         sk(2) = rk - width
         nk(2) = np2
         sk(3) = rk + width
         nk(3) = midnp
         mint = 4
      endif
      nk(mint) = nendk
      sk(mint) = endp
      return
      end

      function acoth(x)
      implicit real*8 (a-h,o-z)
      if (abs(x).le.1.0) stop 'ACOTH defined for |ARG| > 1 only'
      acoth=log((1.0+x)/(x-1.0))/2.0
      end

      function atanh(x)
      implicit real*8 (a-h,o-z)
      if (abs(x).ge.1.0) stop 'ATANH defined for |ARG| < 1 only'
      atanh=log((1.0+x)/(1.0-x))/2.0
      end

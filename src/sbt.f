      subroutine NumSBT(ff,gg,rr,kk,li,np,nexp,nr,fixedk,sgg,singlek)
      
      implicit none
      integer :: i,ll,nr2,nexpi,nexpp,it,ix,lk,kdiv
      integer, parameter :: dp=selected_real_kind(15),fdim=512,lmax=20
C      integer, parameter ::dp=selected_real_kind(precision(1.0d0)),
C     >   fdim=512,lmax=20
      integer, intent(in) :: li,np,nexp,nr
      real(kind=dp), intent(in) :: ff(1:nr),rr(1:nr),kk(1:nr),fixedk
      real(kind=dp), intent(out) :: gg(1:nr),sgg
      real(kind=dp), dimension(1:nr) :: ggs,ya
      real(kind=dp), dimension(1:2*fdim) :: premult
      real(kind=dp), dimension(1:fdim) :: smallr,postdiv
      real(kind=dp) :: xj(0:li),kmin,rmin,rpki,rhomin,kappamin,
     >   dr,factor,pi,dt,xx,C,aa,mindiff
      complex(kind=dp), dimension(1:2*fdim,0:lmax) :: mult_table2
      complex(kind=dp), dimension(1:fdim,0:lmax) :: mult_table1
      complex(kind=dp), dimension(1:2*nr)::ff2,temp1,temp2,temp3,temp4
      complex(kind=dp) :: ci,ca
      logical, intent(in) :: singlek
      logical :: first
      data first/.true./
      save premult,smallr,postdiv,rpki,nexpi,mult_table1,mult_table2,
     >   first
      
                                !   Verify that internal dimensions are adequate
      
                                !!$  open (unit = 99, file = 'tmp2')
                                !!$  do i = 1, nr
                                !!$     write (99,*) kk(i)
                                !!$  end do
                                !!$  close (99)
                                !!$  stop
      
      if ((li > lmax).or.(nr > fdim)) then
         write (6,*) 'stop - dimensions exceeded in sr NumSBT'
         write (6,101) li,lmax
         write (6,102) nr,fdim 
         stop
      endif
      
      rmin = rr(1)
      kmin = kk(1)  
      rhomin = log(rmin)
      kappamin = log(kmin)
      nexpp = nexp+1
      nr2 = 2*nr
      dr = log(rr(2)/rr(1)) 
      pi = 4.0*atan(1.0_dp) 
      dt = 2.0*pi/(nr2*dr)  
      ci = cmplx(0.0,1.0,kind=dp) 
      
                                !   If input meshes have changed, redo the initialization
      
      if (first.or.(rhomin+kappamin /= rpki).or.(nexp /= nexpi)) then
         
         rpki = rhomin+kappamin
         nexpi = nexp          
         first = .false.         
         call INITIALIZE    
      endif
      
                                !   Make the calculation for large k values
      
                                !   extend the input to the doubled mesh, extrapolating the input
                                !   as C r**(np+li)
      
      C = ff(1)/rr(1)**(np+li)   
      ff2(1:nr) = C*smallr(1:nr)**(np+li)
      ff2(nr+1:nr2) = ff   
      
                                !   multiply ff2(i) by r(i)**1.5 and store as input to the fft
      
      temp1 = premult(1:nr2)*ff2
      
      call NLOGN(nexpp,temp1,temp2)  
      
      temp1(1:nr) = temp2(1:nr)*mult_table1(1:nr,li) 
      temp1(nr+1:nr2) = 0.0  
      
                                !    obtain the large k results in the array gg
      
      call NLOGN(nexpp,temp1,temp2) 
      factor = (rmin/kmin)**1.5_dp 
      gg(1:nr) = factor*temp2(nr+1:nr2)*postdiv(1:nr)
      
                                !    obtain the small k results in the array ggs
      
      temp3(1:nr) = rr**3*ff 
      temp3(nr+1:nr2) = 0.0
      call NLOGN(nexpp,temp3,temp4) 
      temp3 = temp4*mult_table2(1:nr2,li)
      call NLOGN(nexpp,temp3,temp4)
      ggs(1:nr) = temp4(1:nr)*dr/nr2
      
                                !    determine the transition point and copy the small k results to gg
      
      kdiv = 0
      mindiff = huge(pi)
      do i = 1,nr
         aa = abs(gg(i)-ggs(i))
         if (aa < mindiff) then
         mindiff = aa
         kdiv = i
      endif
      enddo
      gg(1:kdiv) = ggs(1:kdiv)
      if (.not.singlek) return
      
                                !    obtain the result skk for a single k value fixedk
      
      if (fixedk == 0.0) then
         ix= 1
      else
         ix = (log(fixedk)-kappamin)/dr+1
      endif
      
      if (ix > kdiv) then
         temp4(1) = 1.0
         xx = (log(fixedk)-kappamin)*dt+pi
         ca = exp(ci*xx)
         do i = 2,nr
            temp4(i) = ca*temp4(i-1)
         enddo
         sgg = sum(temp1*temp4)*(rmin/fixedk)**1.5
      else
         do i = 1,nr
            aa = fixedk*rr(i)
            call XJL(aa,li,xj)
            ya(i) = xj(li)*ff(i)*rr(i)**3
         enddo
         sgg = sum(ya)*dr
      endif
      return
 101  format ('li=',i3,' lmax=',i3)
 102  format ('nr=',i4,' fdim=',i4)
      contains
      subroutine INITIALIZE
      implicit none
      integer :: i
      integer, parameter :: dp=selected_real_kind(15),lmax=20
      real(kind=dp), dimension(:,:), allocatable :: j_ltable
      real(kind=dp) :: factor,phi,phi1,phi2,phi3,rad,tt,xj(0:lmax)
      
                                !   Obtain the r values for the extended mesh, and the values r_i^1.5
                                !   in the arrays smallr and premult
                                !
      factor = exp(dr)
      smallr(nr) = rr(1)/factor 
      do i = nr-1,1,-1
         smallr(i) = smallr(i+1)/factor
      enddo
      factor = exp(1.5*dr) 
      premult(nr+1) = 1.0  
      do i = 2,nr
         premult(nr+i) = factor*premult(nr+i-1) 
      enddo  
      premult(nr) = 1.0/factor
      do i = 2,nr  
         premult(nr-i+1) = premult(nr-i+2)/factor
      enddo 
      
                                !   Obtain the values 1/k_i^1/5 in the array postdivide
      
      postdiv(1) = 1.0
      do i = 2,nr
         postdiv(i) = postdiv(i-1)/factor 
      enddo
      
                                !   construct the array of M_l(t) times the phase
      
      do it = 1,nr
         tt = (it-1)*dt         ! Define a t value
         phi3 = (kappamin+rhomin)*tt ! See Eq. (33)
         
         rad = sqrt(10.5_dp**2+tt**2)
         phi = atan((2*tt)/21)
         phi1 = -10*phi-log(rad)*tt+tt+sin(phi)/(12*rad)
     >      -sin(3*phi)/(360*rad**3)+sin(5*phi)/(1260*rad**5)
     >      -sin(7*phi)/(1680*rad**7)
         do ix = 1,10
            phi = 2*tt/(2*ix-1)
            phi1 = phi1+atan((2*tt)/(2*ix-1)) ! see Eqs. (27) and (28)
         enddo
         
         phi2 = -atan(sinh(pi*tt/2)/cosh(pi*tt/2)) ! see Eq. (20)
         phi = phi1+phi2+phi3
         mult_table1(it,0) = sqrt(pi/2)*exp(ci*phi)/nr ! Eq. (18)
         if (it.eq.1) mult_table1(it,0) = 0.5*mult_table1(it,0)
         phi = -phi2-atan(2*tt)
         mult_table1(it,1) = exp(2.0_dp*ci*phi)*mult_table1(it,0) ! See Eq. (21)
         
                                !    Apply Eq. (24)
         
         do lk = 1,lmax-1
            phi = -atan(2*tt/(2*lk+1))
            mult_table1(it,lk+1) = exp(2.0_dp*ci*phi)*
     >         mult_table1(it,lk-1)
         enddo
      enddo
      
                                !   make the initialization for the calculation at small k values
                                !   for 2N mesh values
      
      allocate (j_ltable(nr2,0:lmax))
      
                                !   construct a table of j_l values
      
      do i = 1,nr2
         xx = exp(rhomin+kappamin+(i-1)*dr)  
         call XJL(xx,lmax,xj)
         do ll = 0,lmax
            j_ltable(i,ll) = xj(ll)
         enddo
      enddo
      
      do ll = 0,lmax
         temp1(1:nr2) = j_ltable(1:nr2,ll)
         call NLOGN(nexpp,temp1,temp2) !  Make the FFT See Eq. (35)
         mult_table2(1:nr2,ll) = conjg(temp2)
      enddo  
      deallocate (j_ltable)
      end subroutine INITIALIZE
      
      subroutine XJL(xx,lc,xj)  ! Computes a table of j_l(x) for fixed xx, Eq. (39)
      
      implicit none
      integer :: k,l
      integer, parameter :: dp=selected_real_kind(15)
      integer, intent(in) :: lc
      real(kind=dp), intent(in) :: xx
      real(kind=dp), intent(out) :: xj(0:*)
      real(kind=dp) :: aam,aa,bbm,bb,sa,sb,qqm,aap,bbp,qq,cc
      do l = 0, lc
         xj(l) = 0.0
      enddo
      if (abs(xx) < 1.0d-10) then
         xj(0) = 1.0
         return
      endif
      if (xx.lt.0.75_dp*lc) then
         aam = 1.0
         aa = (2*lc+1)/xx
         bbm = 0.0
         bb = 1.0
         sa = -1.0
         qqm = 1.0d10
         do k = 1,50
            sb = (2*(lc+k)+1)/xx
            aap = sb*aa+sa*aam
            bbp = sb*bb+sa*bbm
            aam = aa
            bbm = bb
            aa = aap
            bb = bbp
            qq = aa/bb
            if (abs(qq-qqm).lt.1.0d-15) exit
            qqm = qq
         enddo
         xj(lc) = 1.0
         if (lc > 0) xj(lc-1) = qq
         if (lc > 1) then
            do l = lc-1,1,-1
               xj(l-1) = (2*l+1)*xj(l)/xx-xj(l+1)
            enddo
         endif
         cc = (sin(xx)/xx)/xj(0)
         do l = 0,lc
            xj(l) = cc*xj(l)
         enddo
      else
         xj(0) = sin(xx)/xx
         if (lc > 0) then
            xj(1) = xj(0)/xx-cos(xx)/xx
         endif
         if (lc > 1) then
            do l = 1,lc-1
               xj(l+1) = (2*l+1)*xj(l)/xx-xj(l-1)
            enddo
         endif
      endif
      return
      end subroutine XJL
      
      subroutine NLOGN(n,xxi,xxo) 
      
                                !   Fast Fourier transform routine.  See Eq. (30).
      
      implicit none  
      integer :: nr,nh,mm(10),lx,lblock,nblock,
     >   iblock,lbhalf,i,j,k,l,istart,jh,ii,nin
      integer, parameter :: dp=selected_real_kind(15)
      integer, intent(in) :: n
      real(kind=dp) :: pi
      complex(kind=dp), intent(in) :: xxi(*)
      complex(kind=dp), intent(out) :: xxo(*)
      complex(kind=dp) :: ww(1024),ci,wk,qq,hold,aa
      logical :: first
      save pi,ww,first,nin
      data first/.true./
      nr = 2**n
      xxo(1:nr) = xxi(1:nr)
      nh = nr/2
      ci = cmplx(0.0,1.0,kind=dp)
      do i = 1,n      
         mm(i) = 2**(n-i)  
      enddo
      lx = 2*mm(1)  
      if ((first).or.(n.ne.nin)) then
         pi = 4.0*atan(1.0_dp)
         do i = 1,nh 
            aa = ci*(i-1)*pi/nr
            ww(i+nh) = ci*exp(aa)
            aa = 2.0*aa
            ww(i) = exp(aa)
         enddo
         first = .false.
         nin = n
      end if
      do l = 1,n      
         nblock = 2**(l-1)   
         lblock = lx/nblock   
         lbhalf = lblock/2   
         k = 0       
         do iblock = 1,nblock 
            wk = ww(k+1)
            istart = lblock*(iblock-1) 
            do i = 1,lbhalf 
               j = istart+i  
               jh = j+lbhalf 
               qq = xxo(jh)*wk
               xxo(jh) = xxo(j)-qq
               xxo(j) = xxo(j)+qq
            enddo
            do i = 2,n 
               ii = i 
               if (k.lt.mm(i)) exit
               k = k-mm(i)  
            enddo
            k = k+mm(ii) 
         enddo
      enddo
      k = 0   
      do j = 1,lx 
         if (k >= j) then
            hold = xxo(j)
            xxo(j) = xxo(k+1)
            xxo(k+1) = hold
         endif
         do i = 1,n 
            ii = i 
            if (k.lt.mm(i)) exit
            k = k-mm(i)
         enddo
         k = k+mm(ii) 
      enddo
      return 
      end subroutine NLOGN
      end subroutine NumSBT
      
      

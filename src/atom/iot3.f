      function cgcigor(l1,l2,l3,m1,m2,m3)
C  The input numbers are integers, = 2 * actual L and M numbers.
      real*8 sym3
c$$$      CGC=((-1)**(NINT(WJ1-WJ2+WM)))*(sqrt(2.*WJ+1.))*
c$$$     >   COF3J(WJ1,WJ2,WJ,WM1,WM2,-WM)
      call threej(l1,l2,l3,m1,m2,-m3,sym3)
      cgcigor = (-1)**((l1-l2+m3)/2) * sqrt(float(l3+1)) * sym3
      return
      end
      
      
      SUBROUTINE IOT3(L1,L2,L3,M1,M2,M3,B)
      REAL*8 B
      call threej(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3,b)
      RETURN
      END

C The following subroutine calculates the Wigner 3J symbols.
C Written by Igor Bray, Flinders Uni, 16/07/90. The input numbers
C are integers, = 2 * actual L and M numbers.      
      subroutine threej(l1,l2,l3,m1,m2,m3,sym3)
      implicit double precision(a-h,o-z)
C FL(I) = log(I!)
      common/fl0/fl(0:2000)
C FAC is used for precision loss checks
      dimension fac(0:2000)

      sym3=0d0
C Check that the L, M pairs are sensible
      if (abs(m1).gt.l1.or.mod(l1-m1,2).ne.0) return
      if (abs(m2).gt.l2.or.mod(l2-m2,2).ne.0) return
      if (abs(m3).gt.l3.or.mod(l3-m3,2).ne.0) return

C Three selection rules
C Sum of the M's must be zero
      if (m1+m2+m3.ne.0) return
C The L's must satisfy the triangular inequality
      if (l1+l2.lt.l3.or.l2+l3.lt.l1.or.l3+l1.lt.l2) return
C The sum of actual L's must be an integer
      if (mod(l1+l2+l3,2).ne.0) return
      
C The often occuring special case where M1 = M2 = M3 = 0 is handled
C by a special formula which does not suffer from any precision loss.
      if (m1.eq.m2.and.m2.eq.m3.and.m3.eq.0) then
         j=(l1+l2+l3)/2
         if (mod(j,2).ne.0) return
         sym3=(-1)**(j/2)*exp((fl(j-l1)+fl(j-l2)+fl(j-l3)-fl(j+1))
     >      /2d0+fl(j/2)-fl(j/2-l1/2)-fl(j/2-l2/2)-fl(j/2-l3/2))
c$$$         print*,'exact result',sym3
         return
      end if 

C Now we use the standard result for the 3J symbols. This is an oscillating
C sum of large numbers and thus leads to precision loss.
      c1=fl((l1+l2-l3)/2)
      c2=fl((l1-l2+l3)/2)
      c3=fl((-l1+l2+l3)/2)
      c4=fl((l1-m1)/2)
      c5=fl((l1+m1)/2)
      c6=fl((l2-m2)/2)
      c7=fl((l2+m2)/2)
      c8=fl((l3-m3)/2)
      c9=fl((l3+m3)/2)
      c10=fl((l1+l2+l3)/2+1)

      const=(c1+c2+c3+c4+c5+c6+c7+c8+c9-c10)/2d0
      kstart=max(0,(l2-l3-m1)/2,(l1-l3+m2)/2)
      kstop=min((l1+l2-l3)/2,(l1-m1)/2,(l2+m2)/2)
      
C The positive terms are kept in SUMPOS and the negative in SUMNEG
      sumpos=0d0
      sumneg=0d0
      i1=(l1+l2-l3)/2
      i2=(l1-m1)/2
      i3=(l2+m2)/2
      i4=(l3-l2+m1)/2
      i5=(l3-l1-m2)/2
      do 10 k=kstart,kstop,1
         c1=fl(k)
         c2=fl(i1-k)
         c3=fl(i2-k)
         c4=fl(i3-k)
         c5=fl(i4+k)
         c6=fl(i5+k)
         fac(k)=(-1)**k/exp(c1+c2+c3+c4+c5+c6-const)
         if (fac(k).lt.0.0) then
            sumneg=sumneg+fac(k)
         else
            sumpos=sumpos+fac(k)
         end if
 10   continue
      
C Now we check for precision loss. If it is too great, we try another
C series that works a little better. If you have extended precision,
C then use that instead. In this case replace 1d7 below with a
C bigger number still. Check the results by comparing with the exact
C answer in the special case when M1 = M2 = M3 = 0.
      if (sumpos-sumneg.lt.1d7) then
         sym3=(-1)**((l1-l2-m3)/2)*(sumpos+sumneg)
c$$$         print*,'Method 1:',sym3
      else
         order=1d-6
         tmp=abs(fac(kstart))
         sum=fac(kstart)
         k=kstart+1
         do while (abs(fac(k)).lt.order.and.k.lt.kstop)
            sum=sum+fac(k)
            k=k+1
         end do
         k0=k
         do while (abs(fac(k)).ge.order.and.k.lt.kstop)
            k=k+1
         end do
         k1=k-1
         do 40 k=k1+1,kstop
            sum=sum+fac(k)
 40      continue
         
         temp=1d0
         do 50 k=k1,k0+1,-1
            temp=1d0-frac(i1,i2,i3,i4,i5,k-1)*temp
 50      continue
         sym3 = (-1)**((l1-l2-m3)/2)*(sum+temp*fac(k0))
c$$$         print*,'Method 2:',sym3
      end if 
      return
      end

      function frac(i1,i2,i3,i4,i5,k)
      real*8 frac
      frac=dfloat((i1-k)*(i2-k)*(i3-k))/dfloat((k+1)*(i4+k+1)*(i5+k+1))
      return
      end 

C  The following routine calculates and stores logs of factorials
      SUBROUTINE FACLOG0
      real*8 fl
      COMMON/FL0/FL(0:2000)
C
C ****  COMPUTE AND STORE, FL(I) = LOG((I)!)
C
      fl(0)=0d0
      do 10 i = 1, 2000
         fl(i) = fl(i-1) + log(dfloat(i))
 10   continue
      RETURN
      END


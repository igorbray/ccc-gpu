      subroutine gsort(nspm,f,maxf,minf,ortint,e1r,lo)
      include 'par.f'
      common /meshrr/ nr, gridr(nmaxr,3)
      
      double precision  ortint(nspmax,nspmax), e1r(nspmax,nspmax), 
     >     sum, r1elk
      dimension  f(nmaxr,nspmax), lo(nspmax)
      dimension maxf(nspmax),minf(nspmax)
      write(4,'(" G - S orthogonalization")') 
      write(20,'(" G - S orthogonalization")')

c     form overlap array <u_j|v_k>, u_j - old nonorthogonal set,
c     v_k - new set but not yet normalized
c     Only elements with  j >= k  required.
      do j=1,nspm
         do k=1,j
            sum=0d0
            do n=1,k-1
               sum = sum + ortint(k,n)*ortint(j,n)/ortint(n,n)
            end do
            ortint(j,k) = ortint(j,k) - sum
c            write(20,'("j,k =",2I3,", <j|k> =",F10.5)')
c     >         j, k, real(ortint(j,k))
         end do
      end do
c     get new e1r(,):  <u_j|H_1|v_k>, j>=k
      do j=1,nspm
         do k=1,j
            sum=0d0
            do n=1,k-1
               sum = sum + ortint(k,n)*e1r(j,n)/ortint(n,n)
            end do
            e1r(j,k) = e1r(j,k) - sum
         end do
      end do
c     <vb_j|H_1|vb_k>, j>=k
      do j=1,nspm
         do k=1,j
            sum=0d0
            do n=1,j-1
               sum = sum + ortint(j,n)*e1r(n,k)/dsqrt(ortint(n,n))
            end do
            e1r(j,k) = (e1r(j,k)/dsqrt(ortint(k,k))-sum)/
     >           dsqrt(ortint(j,j))
            e1r(k,j) = e1r(j,k)
c            write(20,'("j,k =",2I3,", <j|H_1|k> =",F10.5)')
c     >         j, k, real(e1r(j,k))
         end do
      end do
c

      
c     form new orthonomal set vb_k, f_k = vb_k
      do k=1,nspm
         do i=1,nr
            sum=0d0
            do n=1,k-1
               sum = sum + dble(f(i,n))*ortint(k,n)/dsqrt(ortint(n,n))
            end do
            f(i,k) = (dble(f(i,k)) - sum)/dsqrt(ortint(k,k))
         end do
      end do
      do n=1,nspm
         call minmaxi(f(1,n),nr,i1,i2)
         maxf(n) = i2
         minf(n) = i1
         do np=1,nspm
            ortint(n,np) = 0d0
         end do
         ortint(n,n) = 1d0
      end do
c$$$      do n1=1,nspm
c$$$         do n2=n1,nspm
c$$$            rnorm = 0.0D0
c$$$            if(lo(n1).eq.lo(n2)) rnorm =  r1elk(0,f(1,n1),f(1,n2),
c$$$     >         minf(n1),minf(n2),maxf(n1),maxf(n2))
c$$$            write(20,'("n1,n2 =",2I5,", <k1l|k2l> =",F10.5)') n1,n2,
c$$$     >         real(rnorm)
c$$$         end do 
c$$$      end do 

      return
      end













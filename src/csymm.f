      subroutine csidi(a,lda,n,kpvt,det,work,job)
      integer lda,n,job
      complex a(lda,1),det(2),work(1)
      integer kpvt(1)
c
c     csidi computes the determinant and inverse
c     of a complex symmetric matrix using the factors from csifa.
c
c     on entry
c
c        a       complex(lda,n)
c                the output from csifa.
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from csifa.
c
c        work    complex(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                job has the decimal expansion  ab  where
c                   if  b .ne. 0, the inverse is computed,
c                   if  a .ne. 0, the determinant is computed,
c
c                for example, job = 11  gives both.
c
c     on return
c
c        variables not requested by job are not used.
c
c        a      contains the upper triangle of the inverse of
c               the original matrix.  the strict lower triangle
c               is never referenced.
c
c        det    complex(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. abs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c     error condition
c
c        a division by zero may occur if the inverse is requested
c        and  csico  has set rcond .eq. 0.0
c        or  csifa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab
c
c     subroutines and functions
c
c     blas caxpy,ccopy,cdotu,cswap
c     fortran abs,cmplx,iabs,mod,real
c
c     internal variables.
c
      complex ak,akp1,akkp1,cdotu,d,t,temp
      real ten
      integer j,jb,k,km1,ks,kstep
      logical noinv,nodet
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
c
      if (nodet) go to 100
         det(1) = (1.0e0,0.0e0)
         det(2) = (0.0e0,0.0e0)
         ten = 10.0e0
         t = (0.0e0,0.0e0)
         do 90 k = 1, n
            d = a(k,k)
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 30
c
c              2 by 2 block
c              use det (d  t)  =  (d/t * c - t) * t
c                      (t  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (cabs1(t) .ne. 0.0e0) go to 10
                  t = a(k,k+1)
                  d = (d/t)*a(k+1,k+1) - t
               go to 20
   10          continue
                  d = t
                  t = (0.0e0,0.0e0)
   20          continue
   30       continue
c
            det(1) = d*det(1)
            if (cabs1(det(1)) .eq. 0.0e0) go to 80
   40          if (cabs1(det(1)) .ge. 1.0e0) go to 50
                  det(1) = cmplx(ten,0.0e0)*det(1)
                  det(2) = det(2) - (1.0e0,0.0e0)
               go to 40
   50          continue
   60          if (cabs1(det(1)) .lt. ten) go to 70
                  det(1) = det(1)/cmplx(ten,0.0e0)
                  det(2) = det(2) + (1.0e0,0.0e0)
               go to 60
   70          continue
   80       continue
   90    continue
  100 continue
c
c     compute inverse(a)
c
      if (noinv) go to 230
         k = 1
  110    if (k .gt. n) go to 220
            km1 = k - 1
            if (kpvt(k) .lt. 0) go to 140
c
c              1 by 1
c
               a(k,k) = (1.0e0,0.0e0)/a(k,k)
               if (km1 .lt. 1) go to 130
                  call ccopy(km1,a(1,k),1,work,1)
                  do 120 j = 1, km1
                     a(j,k) = cdotu(j,a(1,j),1,work,1)
                     call caxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  120             continue
                  a(k,k) = a(k,k) + cdotu(km1,work,1,a(1,k),1)
  130          continue
               kstep = 1
            go to 180
  140       continue
c
c              2 by 2
c
               t = a(k,k+1)
               ak = a(k,k)/t
               akp1 = a(k+1,k+1)/t
               akkp1 = a(k,k+1)/t
               d = t*(ak*akp1 - (1.0e0,0.0e0))
               a(k,k) = akp1/d
               a(k+1,k+1) = ak/d
               a(k,k+1) = -akkp1/d
               if (km1 .lt. 1) go to 170
                  call ccopy(km1,a(1,k+1),1,work,1)
                  do 150 j = 1, km1
                     a(j,k+1) = cdotu(j,a(1,j),1,work,1)
                     call caxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
  150             continue
                  a(k+1,k+1) = a(k+1,k+1)
     *                         + cdotu(km1,work,1,a(1,k+1),1)
                  a(k,k+1) = a(k,k+1) + cdotu(km1,a(1,k),1,a(1,k+1),1)
                  call ccopy(km1,a(1,k),1,work,1)
                  do 160 j = 1, km1
                     a(j,k) = cdotu(j,a(1,j),1,work,1)
                     call caxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  160             continue
                  a(k,k) = a(k,k) + cdotu(km1,work,1,a(1,k),1)
  170          continue
               kstep = 2
  180       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 210
               call cswap(ks,a(1,ks),1,a(1,k),1)
               do 190 jb = ks, k
                  j = k + ks - jb
                  temp = a(j,k)
                  a(j,k) = a(ks,j)
                  a(ks,j) = temp
  190          continue
               if (kstep .eq. 1) go to 200
                  temp = a(ks,k+1)
                  a(ks,k+1) = a(k,k+1)
                  a(k,k+1) = temp
  200          continue
  210       continue
            k = k + kstep
         go to 110
  220    continue
  230 continue
      return
      end
      subroutine  cswap (n,cx,incx,cy,incy)
c
c     interchanges two vectors.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1),cy(1),ctemp
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = cx(ix)
        cx(ix) = cy(iy)
        cy(iy) = ctemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
   20 do 30 i = 1,n
        ctemp = cx(i)
        cx(i) = cy(i)
        cy(i) = ctemp
   30 continue
      return
      end
      subroutine caxpy(n,ca,cx,incx,cy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1),cy(1),ca
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if (abs(real(ca)) + abs(aimag(ca)) .eq. 0.0 ) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cy(iy) + ca*cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        cy(i) = cy(i) + ca*cx(i)
   30 continue
      return
      end
      subroutine  ccopy(n,cx,incx,cy,incy)
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1),cy(1)
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        cy(i) = cx(i)
   30 continue
      return
      end
      complex function cdotu(n,cx,incx,cy,incy)
c
c     forms the dot product of two vectors.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1),cy(1),ctemp
      integer i,incx,incy,ix,iy,n
c
      ctemp = (0.0,0.0)
      cdotu = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = ctemp + cx(ix)*cy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      cdotu = ctemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ctemp = ctemp + cx(i)*cy(i)
   30 continue
      cdotu = ctemp
      return
      end

      subroutine csifa(a,lda,n,kpvt,info)
      integer lda,n,kpvt(1),info
      complex a(lda,1)
c
c     csifa factors a complex symmetric matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow csifa by csisl.
c     to compute  inverse(a)*c , follow csifa by csisl.
c     to compute  determinant(a) , follow csifa by csidi.
c     to compute  inverse(a) , follow csifa by csidi.
c
c     on entry
c
c        a       complex(lda,n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that csisl or csidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cswap,icamax
c     fortran abs,aimag,amax1,real,sqrt
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      real absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,icamax
      logical swap
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0e0 + sqrt(17.0e0))/8.0e0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (cabs1(a(1,1)) .eq. 0.0e0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = cabs1(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = icamax(k-1,a(1,k),1)
         colmax = cabs1(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0e0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = amax1(rowmax,cabs1(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = icamax(imax-1,a(1,imax),1)
               rowmax = amax1(rowmax,cabs1(a(jmax,imax)))
   50       continue
            if (cabs1(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (amax1(absakk,colmax) .ne. 0.0e0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call cswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call caxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call cswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0e0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call caxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call caxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end

      subroutine csisl(a,lda,n,kpvt,b)
      integer lda,n,kpvt(1)
      complex a(lda,1),b(1)
c
c     csisl solves the complex symmetric system
c     a * x = b
c     using the factors computed by csifa.
c
c     on entry
c
c        a       complex(lda,n)
c                the output from csifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from csifa.
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  csico  has set rcond .eq. 0.0
c        or  csifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call csifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call csisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotu
c     fortran iabs
c
c     internal variables.
c
      complex ak,akm1,bk,bkm1,cdotu,denom,temp
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call caxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call caxpy(k-2,b(k),a(1,k),1,b(1),1)
               call caxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0e0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + cdotu(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + cdotu(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + cdotu(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      integer function icamax(n,cx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      complex cx(1)
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      icamax = 0
      if( n .lt. 1 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = cabs1(cx(1))
      do 30 i = 2,n
         if(cabs1(cx(i)).le.smax) go to 30
         icamax = i
         smax = cabs1(cx(i))
   30 continue
      return
      end

      subroutine cgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      complex a(lda,1),det(2),work(1)
c
c     cgedi computes the determinant and inverse of a matrix
c     using the factors computed by cgeco or cgefa.
c
c     on entry
c
c        a       complex(lda, n)
c                the output from cgeco or cgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from cgeco or cgefa.
c
c        work    complex(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     complex(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if cgeco has set rcond .gt. 0.0 or cgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal,cswap
c     fortran abs,aimag,cmplx,mod,real
c
c     internal variables
c
      complex t
      real ten
      integer i,j,k,kb,kp1,l,nm1
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = (1.0e0,0.0e0)
         det(2) = (0.0e0,0.0e0)
         ten = 10.0e0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (cabs1(det(1)) .eq. 0.0e0) go to 60
   10       if (cabs1(det(1)) .ge. 1.0e0) go to 20
               det(1) = cmplx(ten,0.0e0)*det(1)
               det(2) = det(2) - (1.0e0,0.0e0)
            go to 10
   20       continue
   30       if (cabs1(det(1)) .lt. ten) go to 40
               det(1) = det(1)/cmplx(ten,0.0e0)
               det(2) = det(2) + (1.0e0,0.0e0)
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = (1.0e0,0.0e0)/a(k,k)
            t = -a(k,k)
            call cscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = (0.0e0,0.0e0)
               call caxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = (0.0e0,0.0e0)
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call caxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call cswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end
      subroutine cgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      complex a(lda,1)
c
c     cgefa factors a complex matrix by gaussian elimination.
c
c     cgefa is usually called by cgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for cgeco) = (1 + 9/n)*(time for cgefa) .
c
c     on entry
c
c        a       complex(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that cgesl or cgedi will divide by zero
c                     if called.  use  rcond  in cgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal,icamax
c     fortran abs,aimag,real
c
c     internal variables
c
      complex t
      integer icamax,j,k,kp1,l,nm1
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = icamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -(1.0e0,0.0e0)/a(k,k)
            call cscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call caxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0e0) info = n
      return
      end

      subroutine cgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      complex a(lda,1),b(1)
c
c     cgesl solves the complex system
c     a * x = b  or  ctrans(a) * x = b
c     using the factors computed by cgeco or cgefa.
c
c     on entry
c
c        a       complex(lda, n)
c                the output from cgeco or cgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from cgeco or cgefa.
c
c        b       complex(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  ctrans(a)*x = b  where
c                            ctrans(a)  is the conjugate transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if cgeco has set rcond .gt. 0.0
c        or cgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call cgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran conjg
c
c     internal variables
c
      complex cdotc,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call caxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call caxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  ctrans(a) * x = b
c        first solve  ctrans(u)*y = b
c
         do 60 k = 1, n
            t = cdotc(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/conjg(a(k,k))
   60    continue
c
c        now solve ctrans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + cdotc(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine  cscal(n,ca,cx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, linpack,  3/11/78.
c
      complex ca,cx(1)
      integer i,incx,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = ca*cx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
      complex function cdotc(n,cx,incx,cy,incy)
c
c     forms the dot product of two vectors, conjugating the first
c     vector.
c     jack dongarra, linpack,  3/11/78.
c
      complex cx(1),cy(1),ctemp
      integer i,incx,incy,ix,iy,n
c
      ctemp = (0.0,0.0)
      cdotc = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = ctemp + conjg(cx(ix))*cy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      cdotc = ctemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ctemp = ctemp + conjg(cx(i))*cy(i)
   30 continue
      cdotc = ctemp
      return
      end

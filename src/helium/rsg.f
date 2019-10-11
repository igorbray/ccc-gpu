      subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr)
c
      integer n,nm,ierr,matz
      double precision a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     for the real symmetric generalized eigenproblem  ax = (lambda)bx.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrices  a  and  b.
c
c        a  contains a real symmetric matrix.
c
c        b  contains a positive definite real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  reduc(nm,n,a,b,fv2,ierr)
      if (ierr .ne. 0) go to 50
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  rebak(nm,n,b,fv2,n,z)
   50 return
      end
      subroutine reduc(nm,n,a,b,dl,ierr)
c
      integer i,j,k,n,i1,j1,nm,nn,ierr
      double precision a(nm,n),b(nm,n),dl(n)
      double precision x,y
c
c     this subroutine is a translation of the algol procedure reduc1,
c     num. math. 11, 99-110(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c
c     this subroutine reduces the generalized symmetric eigenproblem
c     ax=(lambda)bx, where b is positive definite, to the standard
c     symmetric eigenproblem using the cholesky factorization of b.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices a and b.  if the cholesky
c          factor l of b is already available, n should be prefixed
c          with a minus sign.
c
c        a and b contain the real symmetric input matrices.  only the
c          full upper triangles of the matrices need be supplied.  if
c          n is negative, the strict lower triangle of b contains,
c          instead, the strict lower triangle of its cholesky factor l.
c
c        dl contains, if n is negative, the diagonal elements of l.
c
c     on output
c
c        a contains in its full lower triangle the full lower triangle
c          of the symmetric matrix derived from the reduction to the
c          standard form.  the strict upper triangle of a is unaltered.
c
c        b contains in its strict lower triangle the strict lower
c          triangle of its cholesky factor l.  the full upper
c          triangle of b is unaltered.
c
c        dl contains the diagonal elements of l.
c
c        ierr is set to
c          zero       for normal return,
c          7*n+1      if b is not positive definite.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      nn = iabs(n)
      if (n .lt. 0) go to 100
c     .......... form l in the arrays b and dl ..........
      do 80 i = 1, n
         i1 = i - 1
c
         do 80 j = i, n
            x = b(i,j)
            if (i .eq. 1) go to 40
c
            do 20 k = 1, i1
   20       x = x - b(i,k) * b(j,k)
c
   40       if (j .ne. i) go to 60
            if (x .le. 0.0d0) go to 1000
            y = dsqrt(x)
            dl(i) = y
            go to 80
   60       b(j,i) = x / y
   80 continue
c     .......... form the transpose of the upper triangle of inv(l)*a
c                in the lower triangle of the array a ..........
  100 do 200 i = 1, nn
         i1 = i - 1
         y = dl(i)
c
         do 200 j = i, nn
            x = a(i,j)
            if (i .eq. 1) go to 180
c
            do 160 k = 1, i1
  160       x = x - b(i,k) * a(j,k)
c
  180       a(j,i) = x / y
  200 continue
c     .......... pre-multiply by inv(l) and overwrite ..........
      do 300 j = 1, nn
         j1 = j - 1
c
         do 300 i = j, nn
            x = a(i,j)
            if (i .eq. j) go to 240
            i1 = i - 1
c
            do 220 k = j, i1
  220       x = x - a(k,j) * b(i,k)
c
  240       if (j .eq. 1) go to 280
c
            do 260 k = 1, j1
  260       x = x - a(j,k) * b(i,k)
c
  280       a(i,j) = x / dl(i)
  300 continue
c
      go to 1001
c     .......... set error -- b is not positive definite ..........
 1000 ierr = 7 * n + 1
 1001 return
      end
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      DOUBLE PRECISION D(N),E2(N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
C
C     This subroutine is a translation of the Algol procedure tqlrat,
C     Algorithm 464, Comm. ACM 16, 689(1973) by Reinsch.
C
C     This subroutine finds the eigenvalues of a symmetric
C     tridiagonal matrix by the rational QL method.
C
C     On input
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E2 contains the squares of the subdiagonal elements of the
C          input matrix in its last N-1 positions.  E2(1) is arbitrary.
C
C      On output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct and
C          ordered for indices 1,2,...IERR-1, but may not be
C          the smallest eigenvalues.
C
C        E2 has been destroyed.
C
C        IERR is set to
C          zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG for  DSQRT(A*A + B*B) .
C
C     Questions and comments should be directed to Burton S. Garbow,
C     Mathematics and Computer Science Div, Argonne National Laboratory
C
C     This version dated August 1987.
C     Modified by C. Moler to fix underflow/overflow difficulties,
C     especially on the VAX and other machines where epslon(1.0d0)**2
C     nearly underflows.  See the loop involving statement 102 and
C     the two statements just before statement 200.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0D0
      T = 0.0D0
      E2(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
         H = DABS(D(L)) + DSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
         if (c .ne. 0.0d0) go to 105
C        Spliting tolerance underflowed.  Look for larger value.
         do 102 i = l, n
            h = dabs(d(i)) + dsqrt(e2(i))
            if (h .gt. t) t = h
  102    continue
         b = epslon(t)
         c = b * b
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = PYTHAG(P,1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
C           Avoid division by zero on next pass
            if (g .eq. 0.0d0) g = epslon(d(i))
            h = g * (p / r)
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      subroutine rebak(nm,n,b,dl,m,z)
c
      integer i,j,k,m,n,i1,ii,nm
      double precision b(nm,n),dl(n),z(nm,m)
      double precision x
c
c     this subroutine is a translation of the algol procedure rebaka,
c     num. math. 11, 99-110(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c
c     this subroutine forms the eigenvectors of a generalized
c     symmetric eigensystem by back transforming those of the
c     derived symmetric matrix determined by  reduc.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix system.
c
c        b contains information about the similarity transformation
c          (cholesky decomposition) used in the reduction by  reduc
c          in its strict lower triangle.
c
c        dl contains further information about the transformation.
c
c        m is the number of eigenvectors to be back transformed.
c
c        z contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        z contains the transformed eigenvectors
c          in its first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
c
      do 100 j = 1, m
c     .......... for i=n step -1 until 1 do -- ..........
         do 100 ii = 1, n
            i = n + 1 - ii
            i1 = i + 1
            x = z(i,j)
            if (i .eq. n) go to 80
c
            do 60 k = i1, n
   60       x = x - b(k,i) * z(k,j)
c
   80       z(i,j) = x / dl(i)
  100 continue
c
  200 return
      end
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end

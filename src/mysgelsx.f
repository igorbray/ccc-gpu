      SUBROUTINE mySGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, 
     $                   RANK, WORK, INFO )
*
*  -- LAPACK driver routine (version 1.0b) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, M, N, NRHS, RANK
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SGELSX computes the minimum-norm solution to
*      min || A * X - B ||
*  using a complete orthogonal factorization. That is, first we compute
*  a QR factorization with column pivoting
*      A * P = Q * [ R11 R12 ]
*                  [  0  R22 ]
*  with the size RANK of R11 such that it is the largest leading
*  submatrix whose condition number is less than 1/RCOND.
*  Then, considering R22 to be negligible, we annihilate R12 through
*  orthogonal transformations from the right, arriving at
*     A * P = Q * [ T11 0 ] * Y
*                 [  0  0 ]
*  The minimum-norm solution is then
*     X = P * Y' [ inv(T11)*Q1'*B ]
*                [        0       ]
*  where Q1 are the first RANK columns of Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right-hand sides. NRHS >=0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the M by N matrix A.
*          On exit, A has been overwritten by its complete orthogonal
*          factorization.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          on entry, the M by NRHS right-hand side.
*          on exit, the first N rows of B contain the minimum-norm
*                   solution.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,M,N).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          on entry: If JPVT(i) <> 0, a(:,i) is permuted
*          to the front of AP, otherwise column i is a free column.
*          on exit: If JPVT(i) = k, then the i-th column of AP
*          was the k-th column of A.
*
*  RCOND   (input) REAL
*          The goal of the complete orthogonal factorization is to
*          identify a well-conditioned triangular matrix whose
*          condition number is less than 1/RCOND.
*
*  RANK    (output) INTEGER
*          The size of the largest leading triangular matrix in the
*          QR factorization (with pivoting) of A, whose condition number
*          is less than 1/RCOND.
*
*  WORK    (workspace) REAL array, dimension
*                      (max( MN+3*N, 2*MN+NRHS )),
*          where MN = min(M,N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            IMAX, IMIN
      PARAMETER          ( IMAX = 1, IMIN = 2 )
      REAL               ZERO, ONE, DONE, NTDONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, DONE = ZERO,
     $                   NTDONE = ONE )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IASCL, IBSCL, ISMAX, ISMIN, J, K, MN
      REAL               ANRM, BIGNUM, BNRM, C1, C2, S1, S2, SMAX,
     $                   SMAXPR, SMIN, SMINPR, SMLNUM, T1, T2
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEQPF, SLABAD, SLAIC1, SLASCL, SLASET, SLATZM,
     $                   SORM2R, STRSM, STZRQF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
      ISMIN = MN + 1
      ISMAX = 2*MN + 1
*
*     Test the input arguments.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGELSX', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF
*
*     Get machine parameters
*
      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
*
*     Scale A, B if max entries outside range [SMLNUM,BIGNUM]
*
      ANRM = SLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL SLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         RANK = 0
         GO TO 100
      END IF
*
      BNRM = SLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL SLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL SLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF
*
*     Compute QR factorization with column pivoting of A:
*        A * P = Q * R
*
      CALL SGEQPF( M, N, A, LDA, JPVT, WORK( 1 ), WORK( MN+1 ), INFO )
*
*     workspace 3*N. Details of Householder rotations stored
*     in WORK(1:MN).
*
*     Determine RANK using incremental condition estimation
*
      WORK( ISMIN ) = ONE
      WORK( ISMAX ) = ONE
      SMAX = ABS( A( 1, 1 ) )
      SMIN = SMAX
      IF( ABS( A( 1, 1 ) ).EQ.ZERO ) THEN
         RANK = 0
         CALL SLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         GO TO 100
      ELSE
         RANK = 1
      END IF
*
   10 CONTINUE
      IF( RANK.LT.MN ) THEN
         I = RANK + 1
         CALL SLAIC1( IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ),
     $                A( I, I ), SMINPR, S1, C1 )
         CALL SLAIC1( IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ),
     $                A( I, I ), SMAXPR, S2, C2 )
*
         if (i.ge.mn-30) print*,I,SMINPR/SMAXPR,'   RANK,SMINPR/SMAXPR'
         IF( SMAXPR*RCOND.LE.SMINPR ) THEN
            DO 20 I = 1, RANK
               WORK( ISMIN+I-1 ) = S1*WORK( ISMIN+I-1 )
               WORK( ISMAX+I-1 ) = S2*WORK( ISMAX+I-1 )
   20       CONTINUE
            WORK( ISMIN+RANK ) = C1
            WORK( ISMAX+RANK ) = C2
            SMIN = SMINPR
            SMAX = SMAXPR
            RANK = RANK + 1
            GO TO 10
         END IF
      END IF
*
*     Logically partition R = [ R11 R12 ]
*                             [  0  R22 ]
*     where R11 = R(1:RANK,1:RANK)
*
*     [R11,R12] = [ T11, 0 ] * Y
*
      IF( RANK.LT.N )
     $   CALL STZRQF( RANK, N, A, LDA, WORK( MN+1 ), INFO )
*
*     Details of Householder rotations stored in WORK(MN+1:2*MN)
*
*     B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS)
*
      CALL SORM2R( 'Left', 'Transpose', M, NRHS, MN, A, LDA, WORK( 1 ),
     $             B, LDB, WORK( 2*MN+1 ), INFO )
*
*     workspace NRHS
*
*     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
*
      CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', RANK,
     $            NRHS, ONE, A, LDA, B, LDB )
*
      DO 40 I = RANK + 1, N
         DO 30 J = 1, NRHS
            B( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
*
*     B(1:N,1:NRHS) := Y' * B(1:N,1:NRHS)
*
      IF( RANK.LT.N ) THEN
         DO 50 I = 1, RANK
            CALL SLATZM( 'Left', N-RANK+1, NRHS, A( I, RANK+1 ), LDA,
     $                   WORK( MN+I ), B( I, 1 ), B( RANK+1, 1 ), LDB,
     $                   WORK( 2*MN+1 ) )
   50    CONTINUE
      END IF
*
*     workspace NRHS
*
*     B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
*
      DO 90 J = 1, NRHS
         DO 60 I = 1, N
            WORK( 2*MN+I ) = NTDONE
   60    CONTINUE
         DO 80 I = 1, N
            IF( WORK( 2*MN+I ).EQ.NTDONE ) THEN
               IF( JPVT( I ).NE.I ) THEN
                  K = I
                  T1 = B( K, J )
                  T2 = B( JPVT( K ), J )
   70             CONTINUE
                  B( JPVT( K ), J ) = T1
                  WORK( 2*MN+K ) = DONE
                  T1 = T2
                  K = JPVT( K )
                  T2 = B( JPVT( K ), J )
                  IF( JPVT( K ).NE.I )
     $               GO TO 70
                  B( I, J ) = T1
                  WORK( 2*MN+K ) = DONE
               END IF
            END IF
   80    CONTINUE
   90 CONTINUE
*
*     Undo scaling
*
      IF( IASCL.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL SLASCL( 'U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA,
     $                INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL SLASCL( 'U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA,
     $                INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF
*
  100 CONTINUE
*
      RETURN
*
*     End of SGELSX
*
      END

      FUNCTION IFACT(N)
      IF(N.GE.13) write(6, 55) N
55    FORMAT(' ',88('*')/'  N=',I5,
     * '  In subr. IFACT  the argument is more then',
     * ' MAX value (13)'/ ' ',88('*'))
      J=1
      DO 1 I=1,N
    1 J=J*I
      IFACT=J
      IF(N .LE.-1) IFACT=0
      RETURN
      END

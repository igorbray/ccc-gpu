      SUBROUTINE WRICO(COEF,N)
      REAL*8    COEF(N)
      WRITE(1 ,1)
    1 FORMAT(/,'  coefficiants:')
      WRITE(1 ,3)(COEF(I),I=1,N)
    3 FORMAT(5X,F8.5  )
      RETURN
      END

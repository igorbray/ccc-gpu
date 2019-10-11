      SUBROUTINE POTATOM(ALFA,BET,H,R1,R6,NE,K,NE5)
      REAL*8 ALFA,BET,H,R1(NE5 ),R6(NE5 ),DD,PP,TT,CC

      DD=0.
      PP=0.
      TT=0.
      DO 1 II=1,NE
         CC=R6(II+1)*(R1(II+2)/(ALFA*R1(II+2)+BET))**2
         DD=(8.*PP-DD+5.*CC)*H/12.
         TT=TT+DD
         TT=TT/(1.+5.*DFLOAT(  K)*H/
     *      (12.*(ALFA*R1(II+2)+BET)))
         DD=PP
         PP=CC-TT*DFLOAT( K)/(ALFA*R1(II+2)+BET)
         R6(II+1)=TT
 1    continue
      
      PP=0.
      DD=0.
      DO 2 IJ=1,NE
         II=NE-IJ+1
         CC=DFLOAT(2 *K+1 )*R6(II+1)/
     *      (ALFA*R1(II+2)+BET)
         DD=(8.*PP-DD-5.*CC)*H/12.
         TT=TT-DD
         TT=TT/(5.*DFLOAT(K+1 )*H/(12.*
     *      (ALFA*R1(II+2)+BET))+1.)
         DD=PP
         PP=DFLOAT(K+1)*TT/(ALFA*R1(II+2)+BET)-CC
         R6(II+1)=TT
 2    continue
      RETURN
      END

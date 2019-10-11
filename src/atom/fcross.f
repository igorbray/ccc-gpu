      SUBROUTINE FCROSS(R,ML,MV,S,OM,CN,
     * SUML,SUMV,SUMM,SUMO,SC,CNST,RL,
     * RV,IL,IV,A,B,II,JJ,NCH,NPHS,BOO,NNL,LH,
     * CONT, NPHS02,IPR)
      REAL*8 R(NPHS02,NPHS02),S(NPHS02,NPHS02),
     * SUML(NPHS),SUMV(NPHS),
     * SUMM(NPHS),SUMO(NPHS),SC(NPHS),
     * OM,CN,CNST,RL,RV,C,A,CNSTB,CNSTA,B
      INTEGER BOO,CONT,NNL(NCH),LH(NCH)
      REAL*8 IL,IV,ML(NPHS),MV(NPHS)
      RL=0.0
      RV=0.0
      IL=0.0
      IV=0.0
      I1=NPHS
      DO 1 I=1,NPHS
      I1=I1+1
      RL=RL+(R(JJ,I)+R(JJ,I1))*ML(I)
1     RV=RV+(R(JJ,I)-R(JJ,I1))*MV(I)
      I1=NPHS
      IF(BOO.EQ.0) GOTO 2
      DO 3 I=1,NPHS
      I1=I1+1
      IL=IL-(S(JJ,I)+S(JJ,I1))*ML(I)
3     IV=IV-(S(JJ,I)-S(JJ,I1))*MV(I)
2     WRITE(6, 4) II,RL,IL,RV,IV
4     FORMAT(2X,'Dipol matrix element',
     * ' B  pCfO  (Channal ',I1,')',
     * /,2X,'RE length ',7X,'IM length ',7X,
     * 'RE rate',5X,'IM rate',/,1X,4(E13.6,3X))
      C=DFLOAT(  NNL(II))/DFLOAT (6 *LH(II)+3 )
      A=C*(RL*RL+IL*IL)*OM/CN
      B=4.*C*(RV*RV+IV*IV)/(OM*CN)
      IF(IPR.NE.1) GOTO 45
      SUML(II)=SUML(II)+SC(JJ)*A
      SUMV(II)=SUMV(II)+SC(JJ)*B
      SUMM(II)=SUMM(II)+2.*C*SC(JJ)*(RL*RV+
     * IL*IV)/CN
      SUMO(II)=SUMO(II)+2.*C*SC(JJ)*ML(JJ)*MV(JJ)
45    IF(CONT.EQ.0) GOTO 5
      CNSTA=CNST*A
      CNSTB=CNST*B
      WRITE(6, 6) CNSTA,CNSTB
6     FORMAT(' cross section in  L-form=',
     * E13.6, 5X,'B V-form=',E13.6)
      A=CNSTA
      B=CNSTB
      GOTO 7
5     CONTINUE
      WRITE(6, 8) A,B
8     FORMAT(2X,'Oscillator strength in  L-form=',
     * E13.6, 3X,'B V-form=',E13.6)
7     CONTINUE
      WRITE(1,9)SUMO(II),SUML(II),SUMV(II),SUMM(II)
9     FORMAT(/,' Total value:  HF MIX',
     * E13.6,/,27X,'L',2X,E13.6,/,27X,'V  ',
     * E13.6,/,25X,'MIX  ',E13.6)
      RETURN
      END

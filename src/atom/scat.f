      SUBROUTINE SCAT(LAML1,LAMM,LAMF,LAMK,LAMD,
     * R1,DV1,DV2,Q,C,SR,RR,RM,PR,MME1,MME2,
     * MME4,MME5,RRR,RE,IM,RED,IMD,AA1,AA2,AA3,AA4,
     * FP,NE,NE5,KS,NWF2,MP,S,LD,KSU1,DL,IT,K,LV,
     * F1,M1,   KK1,KKK1,K21,K22,KIS)
      INTEGER F,S,DL(S),KK1(S),W3,W1,W2,W,F1,FP,IW1(20),IW2(20)
      REAL *8 LAMM(MP),LAMF(FP),LAMK(KS),YMAX(20),
     * LAMD(S),PI,R,A(15),H,BET,Z,ALFA,RO1,EZZ,
     * PRE,R1(NE5),
     *  AA1(KIS),AA2(KIS),AA3(KIS),AA4(KIS),EZP,ZZ1,ZZ2,ZZ3,
     * XNOBL,XPH,LAMI,LAML,LAML1,DELXM,DELXK,
     * DELXF,
     * PR(1),DV1(K21),DV2(K21),C(K21),Q(K21),
     * SR(K22,K22),RR(K21,K21),RM(K21,K21),RE1,IM1,
     * MME1(KKK1,FP,MP),MME2(KKK1,FP,MP),MME4(KKK1,MP),
     * MME5(KKK1,MP),RRR(MP),RE(FP),AA,RE2,IM2,
     * IM(FP),RED(FP),IMD(FP)
      i1000=0
      N1=11
      PI=3.1416D 00
      PRE=1.D-04
      write(6, 2)
2     FORMAT(/' calculation of polirization corrections',
     * ' in scattering faces')
c
      DO 9403 IKIS=1,KIS
c
      call readm(mme1,mme2,
     * mme4,mme5,delxm,delxf,lami,laml,
     * lamm,lamf,lamk,r1,r,h,bet,alfa,ro1,ezp,z,
     * if1,lv,f1,m1,lk,li,ks,k1,k,m,f,jk1,jk,
     * kkk1,fp,mp,  n1, ikis)
c
      IF(IKIS.EQ.1)
     * call komon2(z,r,h,lamk,lamf,lami,laml,lamm,
     * aa1,aa2,aa3,aa4,
     * ne,lk,if1,f1,f,ikis,lv,li,m1,m)
c
      n1 = n1 + 1
9403  CONTINUE
      K=IF1
      K111=K1+1
      KM1=K-1
      DO 400 W=K111,KM1
       DO 5667 IO=1,K
       DO 5667 IOO=1,K
       RR(IO,IOO)=0.0
 5667  RM(IO,IOO)=0.0
      n1=11
      DO 9604 IKIS=1,KIS
c
      rewind n1
      call readm(mme1,mme2,
     * mme4,mme5,delxm,delxf,lami,laml,
     * lamm,lamf,lamk,r1,r,h,bet,alfa,ro1,ezp,z,
     * if1,lv,f1,m1,lk,li,ks,k1,k,m,f,jk1,jk,
     * kkk1,fp,mp,  n1, ikis)
c
        n1 = n1 + 1
      YMAX(IKIS)=0.
      IW1(IKIS)=99
      IW2(IKIS)=99
      DO 410 W2=1,K
      DO 410 W1=1,K
      CALL SELFEN(RR,RM,LAMK(W),DELXM,
     * DELXF,LAMI,LAML,LAMM,LAMF,MME4,
     * MME5,MME1,MME2,RRR,RE,IM,RED,IMD,Q,
     * C,DV1,AA1(IKIS),AA2(IKIS),AA3(IKIS),AA4(IKIS),
     * K21,W,M,M1,F,F1,LV,W1,W2,LK,LI,KKK1,FP,MP,K,IW1(IKIS),IW2(IKIS),
     *  YMAX(IKIS),K1)
410   CONTINUE
9604  CONTINUE
      AA=DSQRT(LAMK(K1+1))
      DELXK=DSQRT(LAMK(IF1 ))-DSQRT(LAMK(IF1 -1))
      CALL CQD11(DELXK,Q,LAMK(W),0.0,
     * LAMK,C,AA,DV1,DV2,K,K1,0,W)
      I1=2*K
      DO 201 I=1,I1
      DO 202 J=1,I1
202   SR(I,J)=0.0
201   SR(I,I)=1.0
      DO 203 I=1,K
      DO 203 J=1,K
      SR(I,J)=SR(I,J)-RR(I,J)*(C(J)-Q(J))*
     * DV1(J)
203   SR(I+K,J+K)=SR(I,J)
      DO 204 I=1,K
      DO 204 J=1,K
      SR(I+K,J)=SR(I+K,J)-RM(I,J)*(C(J)-
     * Q(J))*DV1(J)
204   SR(I,J+K)=-SR(I+K,J)
      WRITE(9,7324)I1,I1,((SR(I,J),I=1,I1),J=1,I1)
 7324 FORMAT(' SR ',I2,'*',I2/100(11(1X,E10.3)/))
      CALL INVER(SR,1.D-05,Y1,2*K,i1000,K22)
      if(i1000 . EQ . 1)  goto 1000
      WRITE(9,205)((RR(I,J),I=1,K),J=1,K)
205   FORMAT(' reak part:',4(2X,E12.5),
     * /200(6(2X,E12.5)/))
      WRITE(9,206) ((RM(I,J),I=1,K),J=1,K)
206   FORMAT('Im. part:',4(2X,E12.5),
     * /200(6(2X,E12.5)/))
      RE1=0.0
      IM1=0.0
      DO 207 I=1,K
      RE1=RE1+RR(I,W)*SR(W,I)+RM(I,W)*SR(W,I+K)
      IM1=IM1+RR(I,W)*SR(W+K,I)+RM(I,W)*SR(W+K,I+K)
      write(6, 208) RE1,IM1
208   FORMAT(' RE1=',E12.5,6X,'IM1=',E12.5)
207   CONTINUE
      RE1=-RE1*PI
      IM1=-IM1*PI
      RE2=PI*DSIGN(1.D 00,RE1)/2.
      IF(DABS(RE1).LT.1.D 00) RE2=
     * DATAN(2.*RE1/(1.-RE1*RE1-IM1*
     * IM1))/2.
      IF(DABS(RE1).GT.1.D 00) RE2=
     * DATAN(2.*RE1/(1.-RE1*RE1-IM1*
     * IM1))/2. + RE2
      IM2=0.25*DLOG(((1.+IM1)**2+RE1*RE1)/
     * ((1.-IM1)**2+RE1*RE1))
      write(6, 209) RE2,IM2
209   FORMAT(' polarization correction to scattering'/
     * ' RE2=',E12.5,6X,'IM2=',E12.5)
      DO 218 J=1,K
      RE1=0.
      DO 217 I=1,K
      RE1=RE1+RR(I,W)*SR(J,I)+RM(I,W)*SR(J,I+K)
217   CONTINUE
      C(J)=RE1
218   CONTINUE
      write(6, 219) (C(J),J=1,K)
219   FORMAT(' pPiB. CigMA:', 8(2X,E12.5))
400   CONTINUE
      RETURN
 1000 CONTINUE
      write(6, 7854)
 7854 FORMAT(' ',66('*')/
     * ' determinantd =0')
       RETURN
      END

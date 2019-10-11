      SUBROUTINE  readm(mme1,mme2,
     * mme4,mme5,delxm,delxf,lami,laml,
     * lamm,lamf,lamk,r1,r,h,bet,alfa,ro1,ezp,z,
     * if1,lv,f1,m1,lk,li,ks,k1,k,m,f,jk1,jk,
     * l1,l2,l3, n1, ikis)
      integer if1,lv,f1,m1,lk,li,ks,k1,k,m,f,jk1,jk
      real*8 ezp,zz1,zz2,zz3,mme1(l1,l2,l3),mme2(l1,l2,l3),
     * mme4(l1,l3),mme5(l1,l3),delxm,delxf,lami,laml,
     * lamm(m), lamf(f), lamk(ks), r1(ne5),
     * r,h,bet,alfa,ro1,  z
      READ(N1) EZP,ZZ1,ZZ2,ZZ3                                          F4500320
      IF1=IFIX(SNGL(ZZ1))                                               F4500340
      F=IFIX(SNGL(ZZ2))                                                 F4500350
      M=IFIX(SNGL(ZZ3))                                                 F4500360
      do 5020 i1=1,if1
      do 5020 i2=1,f
      do 5020 i3=1,m
       READ(n1) MME1(i1,i2,i3)
5020   continue
      do 5021 i1=1,if1
      do 5021 i2=1,f
      do 5021 i3=1,m
       READ(n1) MME2(i1,i2,i3)
5021   continue
      do 5022 i1=1,if1
      do 5022 i2=1,m
       READ(n1) MME4(i1,i2)
5022   continue
      do 5023 i1=1,if1
      do 5023 i2=1,m
       READ(n1) MME5(i1,i2)
5023   continue
       READ(N1) ZZ1,ZZ2,ZZ3                                             F4500410
      LV=IFIX(SNGL(ZZ3))                                                F4500420
      F1=IFIX(SNGL(ZZ2))                                                F4500430
      M1=IFIX(SNGL(ZZ1))                                                F4500440
       READ(N1) ZZ1,ZZ2,ZZ3                                             F4500500
      LK=IFIX(SNGL(ZZ1))                                                F4500510
      LI=IFIX(SNGL(ZZ2))                                                F4500520
      KS=IFIX(SNGL(ZZ3))                                                F4500540
       READ(N1) DELXM,DELXF,LAMI,LAML                                   F4500550
      do 5024 io=1,m
       READ(n1) lamm(io)
5024   continue
      do 5025 io=1,f
       READ(n1) lamf(io)
5025   continue
      do 5026 io=1,ks
       READ(n1) lamk(io)
5026   continue
      READ(N1) ZZ1,ZZ2,R,H,BET,ZZ3,ALFA,RO1            
      do 5027 io=1,ne5
       READ(n1) r1(io)
5027   continue
      K1=IFIX(SNGL(ZZ1))                               
      K=IFIX(SNGL(ZZ2))                                
      IF(IKIS.EQ.1) GOTO 9404                          
      IF(Z.EQ.EZP.AND.JK1.EQ.K1.AND.JK.EQ.K) GOTO 9404 
      write(6, 9405) IKIS,EZP,K1,K,Z,JK1,JK            
9405  FORMAT(' while reading  ',I3,' coulomb matrix',  
     * ' error is found'/
     * ' Z=',E15.8,'  K1=',I3,'  K=',I3,/              
     *  '  previous matrix had :  Z=',E15.8,
     *  '  K1=',I3,'  K=',I3)               
9404  CONTINUE                              
      Z=EZP                                 
      JK1=K1                                
      JK=K                                  
      REWIND  N1                            
      NE=IFIX(SNGL(ZZ3))                    
      RETURN                                
      END                                   

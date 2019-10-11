       PARAMETER (NS=6)
       PARAMETER (NM=14)
       PARAMETER (NK=17)
       PARAMETER (NF=13)
       PARAMETER (NNE=705)
       PARAMETER (NW=NM+NS+NF+NK)
       PARAMETER (NK1=NS+NK)
c
        real*8 laml1,ezz,
     * lamm(NM),me4(NM),me5(NM),me6(NM),r(NM),
     * dv1(NK1),dv2(NK1),q(NK1),d(NK1),lamk(NK1),
     * lamf(NF),re(NF),im(NF),red(NF),imd(NF),
     * lamd(NS),
     * r1(NNE),r2(NNE),r5(NNE),r6(NNE),
     * me1(NF,NM),me2(NF,NM),me3(NF,NM),me(NF,NM),
     * pr(300),
     * mme1(NK1,NF,NM),mme2(NK1,NF,NM),
     * mme4(NK1,NM),mme5(NK1,NM),
     *  func(NW,NNE)
c
       integer f1,f,ne,m1,m,s,k,ld,it,lv,ksu1,
     * nwf2,ne5,k21,k22,k25,k26,
     *  dl(NS),kk1(NS)
c
      character*15  osso,  nale,  prom,  voso,  matr, con
c
       open(2,file='job4.dat')
      read(2,222) con
222   format(15a)
      read(2,222) osso
      read(2,222) nale
      read(2,222) prom
      read(2,222) voso
      read(2,222) matr
c
       read(2,*) f1
       read(2,*) f
       read(2,*) ne
       read(2,*) m1
       read(2,*) m
       read(2,*) s
       read(2,*) k
       read(2,*) ld
       read(2,*) it
       read(2,*) lv
       read(2,*) ksu1
       if(ksu1 . NE . 0)  read(2,*) laml1
c
       ne5= ne+5
       nwf2=m+s+f+k
       k21=k+s
       k22=2*k21
       k25=k21
       k26=k21
c Откpытие файлов
c
        open(1,file=matr,form='unformatted',access='sequential')
c
c      open(6,file=con)
c
       open(2,file=osso,form='unformatted')
       open(3,file=voso,form='unformatted')
       open(4,file=prom,form='unformatted')
       open(8,file=nale,form='unformatted')
c       open(9,file='nul')
c
       call scatwr(laml1,lamm,lamf,lamk,lamd,
     * r1,r2,r6,func,me1,me2,me3,
     * me4,me5,me6,me,r5,dv1,dv2,
     * q,c,pr,mme1,mme2,
     * mme4,mme5,r,re,im,red,imd,
     * f,ne,ne5,k25,nwf2,m,s,ld,ksu1,dl,it,
     * k,lv,f1,m1,kk1,k26,k21,k22)
       end




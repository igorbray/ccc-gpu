       PARAMETER (NSIZEis=15)
       PARAMETER (NSIZEk7=100)
       PARAMETER (NSIZE=705)
       PARAMETER (NSIZEqf=600)
c
      character *15  os, con,  du1,  du2
c
        real*8 z,r,bet,eps,h,dz,
     * gam(NSIZEqf),dir(NSIZEis),teta(NSIZEis),alm3(NSIZEis),r1(NSIZE),
     * r2(NSIZE),r3(NSIZE),r4(NSIZE),r5(NSIZE),r6(NSIZE),
     * r7(NSIZEk7),r8(NSIZE,NSIZEis),
     * mu, c3(NSIZEis),ala(NSIZEis),
     * ala1(NSIZEis),ala2(NSIZEis),alm6du(NSIZEis)

       integer in(NSIZEis),il(NSIZEis),iq(NSIZEis),
     *  ss1,l1,nu,ig,kk(NSIZEqf),
     * is1(NSIZEqf),is2(NSIZEqf),ip(NSIZEis),
     * pp(NSIZEis),ud(NSIZEis),s1(NSIZEqf),kz

       open(2,file='job1.dat')
c
       read(2,222) os
222    format(a15)
       read(2,222) con
       read(2,222) du1
       read(2,222) du2
c
       read(2,*) z
       read(2,*) r
       read(2,*) bet
       read(2,*) eps
       read(2,*) h
       read(2,*) is
       read(2,*) ne
       read(2,*) l1
       read(2,*) ss1
c
       ne5= ne+5
       ne5is=ne5*is
       read(2,*) ksu1
       read(2,*) mu
       read(2,*) dz
       read(2,*) kz
c
       read(2,'(a)')
       read(2,*) (in(i),i=1,is)
       read(2,'(a)')
       read(2,*) (il(i),i=1,is)
       read(2,'(a)')
       read(2,*) (iq(i),i=1,is)
       if( ss1 .LT. 100 ) goto 444
       read(2,'(a)')
       read(2,*) (dir(i),i=1,is)
       goto 445
444    continue
       read(2,'(a)')
       read(2,'(a)')
445    continue
c
       read(2,*) kt
       read(2,'(a)')
       read(2,*) (r7(i),i=1,kt)
c
       if(l1.LT.100)  goto  555
       read(2,'(a)')
       read(2,*) nu
       read(2,*) ig
       read(2,'(a)')
       read(2,*) (is1(i),i=1,nu)
       read(2,'(a)')
       read(2,*) (is2(i),i=1,nu)
       read(2,'(a)')
       read(2,*) (kk(i),i=1,nu)
       read(2,'(a)')
       read(2,*) (gam(i),i=1,nu)
555    continue
c
c Откpытие файлов
       open(6,file=con)
       open(1,file=os,form='unformatted',access='sequential')
       open(3,file=du1)
       open(4,file=du2)
c
       call scfhf(z,gam,dir,r,bet,eps,h,teta,
     * alm3,r1,r2,r3,r4,r5,r6,r7,r8,mu,c3,ala,ala1,
     * ala2,alm6du,
     * dz,kz,
     * is,in,il,iq,
     * ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     * ne5,ne5is,kt,ud,ksu1)
       end




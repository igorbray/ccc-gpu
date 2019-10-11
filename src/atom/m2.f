       PARAMETER (ns1=15)
       PARAMETER (ns2=100)
       PARAMETER (ns3=705)
c
       character*15   os, vo,  con,  du1,  du2
c
       real*8 z,r,bet,eps,h,
     * r8(ns3,ns1)
       real*8  r88(ns3,ns1),
     * gam(ns2),dir(ns1),teta(ns1),r1(ns3),
     * r2(ns3),r3(ns3),r4(ns3),r5(ns3),r6(ns3),r7(100),
     * mu, c3(ns1),ala(ns1),ala3(ns1),
     * ala1(ns1),ala2(ns1),alm6du(ns1),akap,dkap
c
       integer in(ns1),il(ns1),iq(ns1),ss1,l1,nu,ig,kk(ns2),k1,k2,
     * is1(ns2),is2(ns2),ip(ns1),pp(ns1),ud(ns1),s1(ns2)
c
c
       open(5,file='job2.dat')
c
       read(5,222) os
222   format(a15)
       read(5,222) vo
       read(5,222) con
       read(5,222) du1
       read(5,222) du2
c
       read(5,*) z
       read(5,*) r
       read(5,*) bet
       read(5,*) eps
       read(5,*) h
       read(5,*) is
       read(5,*) ne
       read(5,*) l1
       read(5,*) ss1
       ne5= ne+5
       ne5is=ne5*is
       read(5,*) ksu1
       read(5,*) mu
       read(5,*) akap
       read(5,*) dkap
       read(5,*) k1
       read(5,*) k2
       read(5,'(a)')
       read(5,*) (in(i),i=1,is)
       read(5,'(a)')
       read(5,*) (il(i),i=1,is)
       read(5,'(a)')
       read(5,*) (iq(i),i=1,is)
       read(5,'(a)')
       read(5,*) (ud(i),i=1,is)
       if(l1.lt.100) goto 333
       read(5,'(a)')
       read(5,*) (dir(i),i=1,is)
       goto 334
333    continue
       read(5,'(a)')
       read(5,'(a)')
334    continue
       read(5,*) kt
       read(5,'(a)')
       read(5,*) (r7(i),i=1,kt)
c
       if(l1.LT.100) goto 555
       read(5,'(a)')
       read(5,*) nu
       read(5,*) ig
       read(5,'(a)')
       read(5,*) (is1(i),i=1,nu)
       read(5,'(a)')
       read(5,*) (kk(i),i=1,nu)
       read(5,'(a)')
       read(5,*) (gam(i),i=1,nu)
       do 556 i=1,nu
556    is2(i)=is
555    continue
c
c Откpытие файлов
       open(6,file=con)
c Входной файл(pезультат pаботы scfhf(job1))
       open(2,file=os,form='unformatted',access='sequential')
       open(1,file=vo,form='unformatted',access='sequential')
       open(3,file=du1)
       open(4,file=du2)
c
       call fchf(z,gam,dir,r,bet,eps,h,teta,
     * ala,ala1,ala2,ala3,c3,
     * r1,r2,r3,r4,r5,r6,r7,r8,r88,akap,dkap,mu,k1,k2,
     * is,in,il,iq,
     * s1,ss1,l1,nu,ig,kk,is1,is2,ip,pp,ne,
     * ne5,ne5is,kt,ud,ksu1)
       end




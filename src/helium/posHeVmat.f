c This subroutine is for Ps formation m.e. in pos+He.
c In initial channel -  Helim in lsp state + positron
c Final state - Ps + He^+(1S).
c* r in units of Bohr radius = 5.2917720859* 10^(-8)sm
       subroutine posHeme(Nmax,npos,gridk,npk,
     > nchi,nqmi,Li,Ni,lai,sa,nchm,Bnlp,Dnlp,Unlp,
     > gp, nchf,nqmf,Lf,nfa,lfa,KJ,etot,tmpv,lm,
     > Qlp,Q0p,Alp,aq,ionsh,nqgen,nqim ,xp,wp,imax)
  
      include 'par.f'
      include 'par.pos'
      implicit real*8 (a-h,o-z)       
      parameter (maxl=2*ltmax+1)
      double precision Z
      common /Zatom/ Z
      integer npk(nchm+1), lai(KNM),sa(KNM)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl)
c!     >  , Qlfactor(0:ltmax)
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)
      common/gauszR/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      common/gausp/xgp(igpm),wgp(igpm),igp
      common/numericalpotential/ numericalv, lstoppos
      logical numericalv
      dimension  gp(iglp),f0zex(igz)
c      dimension Vab(npk(nchm+1)-1,npk(nchm+1))
      real*8 Qlp(kmax,2*kmax*igp),Q0p(kmax,2*kmax*igp),
     >       fpqb(0:19,2*kmax*igp+kmax),Alp(kmax,2*kmax*igp)
*!!!     >  ,Alex(nqmi,1000,0:10,nspm),f1zex(igz,1000)
      dimension f0z(igz),f1z(igz,2*kmax*igp+kmax)
      real     tmpv(kmax,kmax,0:1),ressums
      real*8 aq(0:kmax),qq(kmax)
      real      gridk(kmax,nchan)
      real*8 xp(2*kmax*igp+kmax),wp(2*kmax*igp)
      dimension a(ncmax),f0(0:ncmax-1),f1(0:ncmax-1)
      dimension res0(0:19), res0ex(0:19),icalam(0:19)
      dimension pc0(0:ncmax-1,ncmax,0:3),pc1(0:ncmax-1,ncmax,0:3)
      dimension resPs0(igz),resPs1(igz)       !,resPs1f(igz,2*nqmi*igp+nqmi)
     >         ,resBnlp(igz),resDnlp(igz),resUnlp(igz)
     >         ,ptmp(igz) !,resPs1f2(igz,2*nqmi*igp+nqmi)
      data pi/3.1415926535897932384626434E0/
      dimension xx(20),yy(20)
      dimension Bnlp(iglp),Dnlp(iglp),Unlp(iglp) 
      real fla,flb,fla1,flb1,fla2,flb2,cof3j,
     > rlli,rK,rllf,rLi,rLf,fl2,flam,fl1,cof12j
       dimension dwp(2*kmax*igp)

          
       if(KJ.gt.lstoppos)  RETURN
         if(sa(2).eq.1.and.Ni.eq.1) return
         if(nqmi.eq.1.or.nqmf.eq.1) then
c          print*, 'no Born and on-shell m.e. for Ps-formation'
!          Beacuse cannot create xgp&wgp composite mesh when nqmi=1
          return
        endif

!        do i=1,iglp; write(15,'(3e20.10)') gp(i), Bnlp(i),Dnlp(i); enddo;

        ze = 2.0
        ielspin=sa(Ni)
!       PRINT*, 'SPIN of target electrons:', ielspin
            rK=KJ
            lli=lai(Ni)
            rLi=float(Li)
            rlli=float(lli)
            llf=lfa   
            rLf=float(Lf)
            rllf=float(llf)
c     Ni,Nf  -   numbers of physical channels

      iphase=llf-lli -(Li-Lf)
      ic0=(-1)**(iphase/2+KJ+Lf)*(2*llf+1)*(2*lli+1)
      c0=dble(ic0)*hat(Lf)*hat(Li)*sqrfct(2*llf)*sqrfct(2*lli)/pi 

      alfa=1.d0    !  -what is it and where it comes from?
      alfa2=alfa*alfa

c       print*,'i->f:', ic0,c0,lli,llf,nam(Ni),Ni
c  Ps:
       
      Nlf=npsstates(2,llf)
      rlam2=rlambda(2,llf)

      if(Nlf.ne.0) then
         do k=1,Nlf
            a(k)=cknd(k,nfa,llf)
c         print*, 'A(k)',llf, nlf,k,a(k)
          end do
c         call coef0(llf,Nlf,a,f0)
c         call coef1(llf,Nlf,a,f1)

         select case(llf)
         case(0)
            call coefsS(llf,Nlf,a,f0,f1)
         case(1)
            call coefsP(llf,Nlf,a,f0,f1)
         case(2)
            call coefsD(llf,Nlf,a,f0,f1)
         case(3)
            call coefsF(llf,Nlf,a,f0,f1)
         case default
            stop 'coeffs are missing for lb>3'
         end select ! llf

         do k=0,Nlf-1
            pc0(k,nfa,llf)=f0(k)
            pc1(k,nfa,llf)=f1(k)
         end do
       endif
******

       knqm=nqmi
       ntst=0  !     !-2 !  for full - 0,  for born -2
       nqm=-1 !  1  for on-shell!
c       call qmesh(gridk(1,nchi),knqm,imax,xp,wp,aq,ionsh,nqgen,nqim)
          imaxp=2*nqgen*igp !-igp
          mnqmi=nqmi
          mnqmf=nqmf    
        if(ntst.eq.-2 .or.nqm.eq.1) then
          mnqmi=1
          mnqmf=1
        endif
c        print*,'q-mesh::', imaxp,xp(imax),imax,xp(imaxp),nqgen,igp

c        do im=1,nqmi
c        print*, im,gridk(im,nchi)
c        enddo
!        do im=1,imax
!        if(xp(im).gt.50.) wp(im)=0.d0
!        enddo

       dp=1.7
c      print*, 'Ps:', nchf,mnqmf,mnqmi,nam(Ni)
      ddp=1.e+08
      pdp=0.0065d0
         nd1=2
         intp=2*nd1+1
       nnmin=INT(Log(1.E-06*ddp)/pdp) !800
       nnmax1=INT(Log(36.0e+02*ddp)/pdp)                  !iglp-2
       gpmin=gp(nnmin)
       gpmax=gp(nnmax1)      

 
! iqstep below is to speed up the q-integral when there are no singularities
             iqstep =1
             iqs1 = 0
             dwp(:) = wp(:)  !0.d0
!          IF(MOD(igp,3).eq.0) THEN
!             iqstep = 3
!              iqs1=1
!             do i=1+1,2*igp*nqim,iqstep
!                dwp(i) = wp(i-1) + wp(i) + wp(i+1)
!             enddo
!          ELSE
!             STOP'change igp to 6 or 9 or 12, it must be multiple of 3'
!          ENDIF

!       write(21,*)
!       write(21,'(a1,3i5)') '#', npos,nfa,llf
!      do ip = 1, iglp
!         pb2 = gp(ip)
!        call f0zpart(Nlf,rlam2,2.d0,npos,llf,nfa,pb2,
!     >        res0ps,res1ps,pc0,pc1)
!         write(21,'(3e20.10,i5)') pb2,res0ps,res1ps,ip
!      enddo

  
c$$$C$OMP PARALLEL DO
c$$$C$OMP& SCHEDULE(dynamic)
c$$$C$OMP& default(private)
c$$$C$OMP& shared(sqrfct,ortint)
c$$$C$OMP& shared(xgz,pl,igz,igp,ddp,pdp)
c$$$C$OMP& shared(tmpv,Nlf,rlam2,gridk,pi,Bnlp,Dnlp,Unlp,gp,etot,nfa)
c$$$C$OMP& shared(rK,KJ,Li,Lf,lli,rLi,rlli,llf,rLf,rllf,nchi,nchf,npos)
c$$$C$OMP& shared(c0,alfa,alfa2,a,pc0,pc1,Q0p,Qlp,Alp,hat)   !add Alex
c$$$C$OMP& shared(nqmi,nqim,mnqmi,mnqmf,imaxp,xp,wp,aq,ionsh,imax,maxf)
c$$$C$OMP& shared(nam,Ni,na,loi,Ci,dp,nspm,mpg,mpgmax,pg,fl,gridr)
       DO kf = 1, mnqmf      ! 1 - for born  ! kfstart,kfstop
c           print*, kf, mnqmf,gridk(kf,nchf)
            qb=dble(gridk(kf,nchf))
            if(qb.lt.0.) cycle !.or.qb.gt.20.) cycle
c          print '(i3,$)', kf
           iq=0
           sumf0z=0.
         do ki=1,mnqmi ! 1 - for born !kistart, kistop
            qa=dble(gridk(ki,nchi))
c         print*, 'gridk', nchf,npos,qb,qa
            qa2=qa*qa
            if(qa.lt.0.) cycle !.or.qa.gt.20) cycle
!!      IF((abs(qa-qb).lt.0.1.or.abs(qa-qb/2).lt.0.1).and.qa.gt.1) CYCLE
            iq = iq+1

        qb2=qb*qb
        qa2=qa*qa
        qbqa=qb*qa
        qbE=qa2 - etot 
                 
          f0z(:)=0.d0
!          f0zex(:)=0.d0
          if(iq.eq.1) f1z(:,:)=0.d0
cc          if(iq.eq.1) f1zex(:,:)=0.d0
          resPs0(:)=0.d0
          resPs1(:)=0.d0

       resDnlp(:)=0.d0
       do iz=1,igz
          xz=xgz(iz)
        pa2 =qb2 + qa2 - 2.d0*qbqa*xz
        nn=INT(Log(pa2*ddp)/pdp)
        if(nn.le.nd1+1.or.nn.gt.iglp-nd1) cycle

       IF(pa2.gt.gpmin.and.pa2.lt.gpmax) THEN
       yy(1:intp)=Bnlp(nn-nd1:nn+nd1)
       xx(1:intp)=gp(nn-nd1:nn+nd1)
       call intrpl(intp,xx,yy,1,pa2,ResBnlp(iz))
     
       yy(1:intp)=Unlp(nn-nd1:nn+nd1)
       call intrpl(intp,xx,yy,1,pa2,ResUnlp(iz))
        
       yy(1:intp)=Dnlp(nn-nd1:nn+nd1)
       call intrpl(intp,xx,yy,1,pa2,ResDnlp(iz))    
       ENDIF

       IF(pa2.le.gpmin) THEN
        ResBnlp(iz)=Bnlp(nnmin)
        ResUnlp(iz)=Unlp(nnmin)
        ResDnlp(iz)=Dnlp(nnmin)
       ELSEIF(pa2.ge.gpmax) THEN
        ResBnlp(iz)=0. !Bnlp(nnmax1)*gp(nnmax1)/pa2
        ResUnlp(iz)=0. !Unlp(nnmax1)*gp(nnmax1)/pa2 
        ResDnlp(iz)=0. !Dnlp(nnmax1)*gp(nnmax1)/pa2
       ENDIF

        enddo

       do iz=1,igz
          xz=xgz(iz)
        pb2 = 0.25d0*qb2 + qa2 - qbqa*xz
    
        call f0zpart(Nlf,rlam2,2.d0,npos,llf,nfa,pb2,
     >        resPs0(iz),resPs1(iz),pc0,pc1)
   
        efactor = 0.25d0*qb2 + pb2 - etot
c        efactor2 = 0.5d0*qa2 + 0.5d0*pa2 - etot  
        temp = (efactor*resPs1(iz)-resPs0(iz))*resBnlp(iz)
     >          + resPs1(iz)*resDnlp(iz) 
c!     >       - resPs1(iz)*resUnlp(iz)  
c!     >       - resPs0(iz)*resBnlp(iz)
        f0z(iz)=f0z(iz)+temp
c!        f0zex(iz)=(-1)**ielspin*resPs1(iz)*resDnlp(iz) 
        enddo
      
       if(iq.eq.1) then  ! change it back to iq.eq.1
         do ip=1+iqs1,2*nqim*igp, iqstep
            pp=xp(ip)
!         if(pp.gt.100.) cycle
            qbpp=qb*pp
            pp2=pp*pp

        ResBnlp(:) = 0.d0
       do iz=1,igz
          xz=xgz(iz)
        pa2 =qb2 + pp2 - 2.d0*qbpp*xz
          nn=INT(Log(pa2*ddp)/pdp)
          if(nn.le.nd1+1.or.nn.gt.iglp-nd1) cycle

       IF(pa2.gt.gpmin.and.pa2.lt.gpmax) THEN
        yy(1:intp)=Bnlp(nn-nd1:nn+nd1)
        xx(1:intp)=gp(nn-nd1:nn+nd1)
        call intrpl(intp,xx,yy,1,pa2,ResBnlp(iz))
       ENDIF

       IF(pa2.le.gpmin) THEN
        ResBnlp(iz)=Bnlp(nnmin)
       ELSEIF(pa2.ge.gpmax) THEN
        ResBnlp(iz)=0. !Bnlp(nnmax1)*gp(nnmax1)/pa2 !0.d0
       ENDIF
!
!      enddo
!           do iz=1,igz
!              xz=xgz(iz)

              pb2=0.25*qb2 + pp**2.d0 - qbpp*xz
        call f0zpart(Nlf,rlam2,2.d0,npos,llf,nfa,pb2,
     >         resPs01,resPs1f,pc0,pc1) 
       f1z(iz,ip)=resPs1f*resBnlp(iz)
!       f1z(iz,ip-1) = f1z(iz,ip)
!       f1z(iz,ip+1) = f1z(iz,ip)
       enddo
       enddo
!!!!
         do ip=2*igp*nqim+1, 2*igp*nqim+nqim
            pp=xp(ip)
         if(pp.gt.50.) cycle
            qbpp=qb*pp
            pp2=pp*pp

        ResBnlp(:) = 0.d0
       do iz=1,igz
          xz=xgz(iz)
        pa2 =qb2 + pp2 - 2.d0*qbpp*xz
          nn=INT(Log(pa2*ddp)/pdp)
          if(nn.le.nd1+1.or.nn.gt.iglp-nd1) cycle

       IF(pa2.gt.gpmin.and.pa2.lt.gpmax) THEN
        yy(1:intp)=Bnlp(nn-nd1:nn+nd1)
        xx(1:intp)=gp(nn-nd1:nn+nd1)
        call intrpl(intp,xx,yy,1,pa2,ResBnlp(iz))
       ENDIF

       IF(pa2.le.gpmin) THEN
        ResBnlp(iz)=Bnlp(nnmin)
       ELSEIF(pa2.ge.gpmax) THEN
        ResBnlp(iz)=0. !Bnlp(nnmax1)*gp(nnmax1)/pa2 !0.d0
       ENDIF
!
!      enddo
!           do iz=1,igz
!              xz=xgz(iz)

              pb2=0.25*qb2 + pp**2.d0 - qbpp*xz
        call f0zpart(Nlf,rlam2,2.d0,npos,llf,nfa,pb2,
     >         resPs01,resPs1f,pc0,pc1) 
       f1z(iz,ip)=resPs1f*resBnlp(iz)
       enddo
       enddo
        endif

c       if(iq.eq.1) then
c         write(15,*)
c          do i=1,imax
c          write(15,'(i6, 3e20.10)') i,xp(i),f1z(1,i),Qlp(ki,i)
c          enddo
c       endif
c*** Sum over l1,l2, l1',l1'',lambda ...
 
          resultex=0.
cc       if(ni1.gt.1) 
cc     >   call exchpart(ni1,lli,llf,Li,Lf,KJ,
cc     > nqmi,nchi,igp,imaxp,xp,wp,ki,qa,qb,gridk,nspm,
cc     > f1zex,Alex,resultex) 
  
           icalam(:) = 0
           res0(:) = 0.d0
c
           sumlb1=0.d0
           do lb1=0,llf 
               lb2=llf-lb1
               flb1=float(lb1)
               flb2=float(lb2)
               clb=sqrfct(2*lb1)*sqrfct(2*lb2)*2.d0**lb1
               sumla1=0.d0
            
               do la1=0,lli
                  la2=lli-la1
                  lab2=lb2+la2
                  qal2=qa**(lab2+1) ! here multiplyed by qa
                  fla1=float(la1)
                  fla2=float(la2)
                  cla=1.d0/sqrfct(2*la1)/sqrfct(2*la2)
                  cla=cla*qb**(lb1+la1)
                  suml1=0.d0

c******
                   do l1=iabs(lb1-la1),lb1+la1,2
                     fl1=float(l1)
                     w3j1=cof3j(fla1,flb1,fl1,0.,0.,0.)
                     if(w3j1.eq.0.) cycle
                     suml2=0.d0
                     do l2=iabs(lb2-la2),lb2+la2,2
                        fl2=float(l2)

                        w3j2=cof3j(fla2,flb2,fl2,0.,0.,0.) 
!     >                     * (-1)**(lb1-lb2)  ! this andi (-1)**llf above are different from pos-H code, 
!                                              ! they come from diff coordinate choice, but they cancel each other

                        if(w3j2.eq.0.) cycle
                        reslam=0.d0
                        do lam=max0(iabs(Lf-l1),iabs(Li-l2)),
     >                     min0(Lf+l1,Li+l2),2
                           flam=float(lam)
                           w3j3=cof3j(flam,fl1,rLf,0.,0.,0.)
                           w3j4=cof3j(fl2,flam,rLi,0.,0.,0.)
                           w3j34=dble(w3j3)*dble(w3j4)
                           if(w3j34.eq.0.d0) cycle
                           w12j=cof12j(fla1,rlli, rK,  rllf,
     >                                 fla2,rLi,rLf,flb1,
     >                                 flb2,fl2, flam,fl1)
                           wigner=w3j34*dble(w12j)
                           if(wigner.eq.0.d0) cycle
c      write(*,'(f4.1)') fla1,rlli,rK,rllf,fla2,rLi,rLf,flb1,flb2,fl2
c      PRINT*,'COEFS:', w12j,w3j3,w3j4,w3j34
c      PAUSE    
               if(icalam(lam).ne.1) then
                          resiz=0.d0
                       do iz=1,igz
                         resiz=resiz+f0z(iz)*pl(lam,iz)
                        enddo    
 
                         resizex=0.d0  
c              if(ni1.gt.1) then
c                    do iz=1,igz
c                       resizex=resizex+f0zex(iz)*pl(lam,iz)
c                    enddo
c               endif 
c               write(15,'(4e20.10)') qa,resiz,resizex
                        res0(lam)=resiz !+resizex
c              print*,'res0:',resiz,resizex
!      write(*,'(a13,2i4,4e20.10)')'kf,ki,vmatt:',kf,ki,qb,qa,resiz,etot

                        if(iq.eq.1) then
                    fpqb(lam,:)=0.d0
c                    sumfpqb=0.d0
                    do ipl=1+iqs1,2*nqim*igp,iqstep
                       pp=xp(ipl)
                       resi=0.d0
                       if(pp.gt.100.d0) cycle
                          do iz=1,igz
                             resi=resi+f1z(iz,ipl)*pl(lam,iz)
                          enddo
                          fpqb(lam,ipl)=resi
!                          fpqb(lam,ipl-1)=resi
!                          fpqb(lam,ipl+1)=resi                        
                    enddo
!!!
                      do ipl=2*nqim*igp+1,2*nqim*igp+nqim
                        pp=xp(ipl)
                         resi=0.d0
                      if(pp.gt.100.d0) cycle
                         do iz=1,igz
                             resi=resi+f1z(iz,ipl)*pl(lam,iz)
                            enddo
                         fpqb(lam,ipl)=resi
                        enddo
c        write(14,'(5i5,2e16.8)') nchi,nchf,ni1,kf,lam,sumfpqb,res0(lam)
                        endif
                  icalam(lam)=1
                  endif !icalam.ne.1

       res1 = 0.d0      
! Integration over composite mesh:
c             if(gridk(1,nchi).gt.0.d0) then! Rav: subtraction
!                           method must be used for both open and closed
!                           channels. Some changes have been made to
!                           make it work. 
             ising=ki    ! singular point
             if(ki.le.ionsh) ising=ki-1
             if(ki.eq.1) ising=ionsh
         if(gridk(1,nchi).lt.0.) ising=ki-1
c        print*,'Statrt Integr', nqmi,ising,ki,ionsh,nqgen 

            if(ising.ne.nqim) then   !i.e. singularity is not at last k-mesh
            res1=0.d0
*   integrals coming before singularity
            do i=1+iqs1,2*(ising-1)*igp,iqstep
             pp=xp(i)
             fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) !-Alp(ki,i))
             res1=res1+fp * dwp(i) 
             enddo

*   integral with singularities
             res2=0.d0
             do i=2*(ising-1)*igp+1,(2*ising-1)*igp
             pp=xp(i)
             fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) !-Alp(ki,i))
     >        - fpqb(lam,2*nqim*igp+ising)*qa**(lab2+1)*Q0p(ki,i)
             res2=res2+wp(i)*fp
             enddo

            res1b=2*qa*Log(2.)-qa*Log(qa-aq(ising-1))+2*qa*Log(qa)-
     >            qa*Log(aq(ising-1)+qa)-
     >          aq(ising-1)*Log((aq(ising-1)+qa)/(-aq(ising-1)+qa))

           res1b=res1b*fpqb(lam,2*nqim*igp+ising)*qa**(lab2+1)
           res2=res2+res1b
     
           res3=0.d0
           do i=(2*ising-1)*igp+1, 2*ising*igp
           pp=xp(i)
           fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) !-Alp(ki,i))
     >        - fpqb(lam,2*nqim*igp+ising)*qa**(lab2+1)*Q0p(ki,i)  ! check the Q0p and ...
           res3=res3+wp(i)*fp
           enddo

           res1c=-2*qa*Log(2.)+qa*Log(aq(ising)-qa)-
     >           2*qa*Log(qa)+qa*Log(aq(ising)+qa)
     >          +aq(ising)*Log((aq(ising)+qa)/(aq(ising)-qa))
           res1c=res1c*fpqb(lam,2*nqim*igp+ising)*qa**(lab2+1) 
           res3=res3+res1c
c           print*,'res3:', res3,res1c    
* Integrals coming after singularity
          res4=0.d0
          do i=2*ising*igp+1+iqs1,2*nqim*igp,iqstep
             pp=xp(i)
             fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) !-Alp(ki,i))
             res4=res4+fp * dwp(i)
          enddo

c          print*,'RES::', res1,res2,res3,res4
          res1=res1+res4+res2+res3
          else
*  singularity at last k-mesh point
          res1=0.d0
          do i=1+iqs1,2*(nqim-1)*igp,iqstep
             pp=xp(i)
             fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) !-Alp(ki,i))
             res1=res1+fp *dwp(i)
             enddo

*   integ with singularities
            res2=0.d0
          do i=2*(nqim-1)*igp+1,(2*nqim-1)*igp
            pp=xp(i)
            fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) !-Alp(ki,i))
     >        -fpqb(lam,2*nqim*igp+nqim)*qa**(lab2+1)*Q0p(ki,i)
           res2=res2+wp(i)*fp

          enddo

        res1b=2*qa*Log(2.)-qa*Log(qa-aq(nqim-1))
     >        +2*qa*Log(qa)-qa*Log(aq(nqim-1)+qa)
     >        -aq(nqim-1)*Log((aq(nqim-1)+qa)/(-aq(nqim-1)+qa))
        res1b=res1b*fpqb(lam,2*nqim*igp+nqim)*qa**(lab2+1)      ! check lab2+1??
        res2=res2+res1b

            res3=0.d0
            do i=(2*nqim-1)*igp+1,2*nqim*igp
               pp=xp(i)
           fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) !-Alp(ki,i))
     >      -fpqb(lam,2*nqim*igp+nqim)*qa**(lab2+1)*(qa2+alfa2)**2.d0/
     >         (pp*pp+alfa2)/(pp*pp+alfa2)*Q0p(ki,i)
           res3=res3+wp(i)*fp
                      enddo

          res1c=(-2*qa*ATan(qa/alfa)+Pi*qa)/(2.*alfa*(alfa2+qa2))
          res1c=res1c*fpqb(lam,2*nqim*igp+nqim)*
     >           qa**(lab2+1)*(qa2+alfa2)*(qa2+alfa2)
                   ! CHECK qa**lab2 or **lab2+1
          res3=res3+res1c

          res1=res1 +res2+res3
          endif
  

        resAlp=0.d0
***  this additional term with Alp doesn't contain singularity 1/q^2 
         do i=1+iqs1,2*nqim*igp,iqstep
           pp=xp(i)
           fp=fpqb(lam,i) *pp**(lab2+1)*Alp(ki,i)  
           resAlp=resAlp+fp * dwp(i) 
          enddo 

         res1=res1+resAlp 

!  Integration over q-mesh is done!
                resal = qal2 * res0(lam)
!                   resalex=qal2*res0ex(lam)
c               res1=0.   !!!!!!!!!!!!!!!!
c               resal=0.
                    result = res1/pi + resal 
c          PRINT*, 'I&II:',res1/pi, resal
                    reslam=reslam+(2*lam+1)*wigner*result

c       write(14,'(2i4,2f10.6,2e16.6)') nchi,nchf,qa,qb,result,wigner
c           print*,'coefs:', wigner,w3j2,w3j1,cla,clb,c0        
                     end do
                      suml2=suml2+(2*l2+1)*reslam*w3j2
                     end do
                     suml1=suml1+(2*l1+1)*suml2*w3j1
                  end do
                  sumla1=sumla1+suml1*cla
               end do
               sumlb1=sumlb1+sumla1/clb
            end do
        ressums = sumlb1*c0*qb * sqrt(2.d0) !!!+ resultex !! 
!       if(lli.eq.1) write(*,'(a13,3i4,a3,2i4,e20.10)')'kf,ki,vmatt:'
!     > ,lli,Li,KJ,'$$$',kf,ki,ressums
c          ressums=resultex
c        PRINT*, 'II:',qa,qb,ressums, resultex        
c            ressums0=ressums !-ressumex
c            resj1j2=ressums0 !const*ressums0         
* V-matrix is multiplied by qb*qa. Namely, qb*qa*V is used as driving term
* and in the kernal of integral equations. The solution is qb*qa*T.
* However, CCC prints  out  just V without qb*qa factor.
* (multiplification by qa is done earlier, see 'result')  -- check here carefully!
* V-matrix is also multiplied by sqrt(mu=2)

            if (nchf.ge.nchi) then
               tmpv(kf,ki,0) = tmpv(kf,ki,0) + ressums 
               tmpv(kf,ki,1) = tmpv(kf,ki,1) + ressums 
            else
               tmpv(ki,kf,0) = tmpv(ki,kf,0) + ressums
               tmpv(ki,kf,1) = tmpv(ki,kf,1) + ressums
            endif
c         sumj1j2=sumj1j2+resj1j2
        end do
c       write(14,*)
c       print*, 'kf:::',kf,qb
       END DO
!       write(*,*) 'lli,Li,KJ',lli,Li,KJ
c$$$C$OMP END PARALLEL DO

c        write(14,'(4i4, e16.6)') nchi,nchf, ki,kf,sumj1j2
c.
c.
c this is for rearranging the indexes as it was done in vmat.f.
c        print*, 'V-Ps::', tmpv(1,1,0), tmpv(1,1,1),
c     * tmpv(nqmi,nqmf,0),tmpv(nqmi,nqmf,1) 
 
c$$$        do ki = 1, nqmi
c$$$           kii = npk(nchi) + ki - 1
c$$$           do kf = 1, nqmf
c$$$          kff = npk(nchf) + kf - 1
c$$$      if (kff.ge.kii) then
c$$$                  Vab(kff,kii) = Vab(kff,kii) + tmpv(kf,ki,0)
c$$$                endif
c$$$                 enddo
c$$$c          print*, 'Ps:',ki,tmpv(ki,1,0) 
c$$$               enddo
       RETURN
       END  

c*****************************************************************
c  this subroutine is to create composite mesh
c  for q-integration in posvmat
c*****************************************************************
      subroutine qmesh(gki,nqmi,imax,xp,wp,aq,ionsh,nqgen,nqim)
      include 'par.f'
      include 'par.pos'
      implicit real*8 (a-h,o-z)       
      common/gausp/xgp(igpm),wgp(igpm),igp
      real gki(kmax)
      dimension aq(0:kmax),qq(kmax)
      dimension xp(2*kmax*igp+kmax),wp(2*kmax*igp)

      if(gki(1).gt.0.d0) then
*     putting on-shell point to its place

         ionsh=1
         do i=2,nqmi
            if(gki(i).lt.gki(1)) then
               qq(i-1)=dble(gki(i))
               ionsh=i
            else
               qq(i)=dble(gki(i))
            endif
         end do
         qq(ionsh)=dble(gki(1))

         do ji=1,nqmi
            j0=2*nqmi*igp+ji
            xp(j0)=qq(ji)
         end do
         nqgen=nqmi
         nqim=nqmi
         imax=2*nqmi*igp+nqmi

      else
*     on-shell point is negative (no on-shell point)

        ionsh=-1
        if (nqmi.eq.1) then
         nqim=1
         nqgen=nqim
         imax=2*nqim*igp+nqim
         else 
         do i=2,nqmi
           qq(i-1)=dble(gki(i))
           end do
              nqim=nqmi-1
              nqgen=nqim
              imax=2*nqim*igp+nqim
         do ji=1,nqim
            j0=2*nqim*igp+ji
            xp(j0)=qq(ji)
         end do
        endif

*### the way as below gives wrong nqim in some cases.. why?
!         do i=2,nqmi
!            qq(i-1)=gki(i)
!         end do
!         nqim=nqmi-1
!         nqgen=nqim
!         imax=2*nqim*igp
      endif

*     end points
      aq(0)=0.d0
      do i=1,nqgen-1
         aq(i)=(qq(i+1)+qq(i))*0.5d0 ! aq(i) are the middle points
*     forming composite mesh
         do ji=1,igp
            j1=2*(i-1)*igp+ji
            xp(j1)=aq(i-1)+
     >         (qq(i)-aq(i-1))*xgp(ji) ! linear scaling
            wp(j1)=(qq(i)-aq(i-1))*wgp(ji) ! linear scaling 
            j2=(2*i-1)*igp+ji
            xp(j2)=aq(i)-
     >         (aq(i)-qq(i))*xgp(igp+1-ji)
            wp(j2)=(aq(i)-qq(i))*
     >         wgp(igp+1-ji)
         end do
      end do
*     last 2 intervals
      do ji=1,igp
         j1=2*(nqgen-1)*igp+ji
         xp(j1)=aq(nqgen-1)+
     >      (qq(nqgen)-aq(nqgen-1))*xgp(ji)
         wp(j1)=(qq(i)-aq(i-1))*wgp(ji)
*     last interval
         j2=(2*nqgen-1)*igp+ji
         xp(j2)=qq(nqgen)/xgp(igp+1-ji)
         wp(j2)=qq(nqgen)*wgp(igp+1-ji)/
     >      xgp(igp+1-ji)/xgp(igp+1-ji)
      end do
*     composite mesh is ready 
      return
      end

c*****************************************
C this is for Ps-Ps direct V m.e.
!  uses fd0 for R < r_{max} of (wafe function)    
c*****************************************
      subroutine pspsdme(nqmi,psii,maxpsii,lia,li,nia,chii,
     >   minchii,gki,nqmf,psif,maxpsif,lfa,lf,nfa,chif,minchif,
     >   gkf,lg,nchf,nchi,nchm,npk,pos,ne2e,vdon,vmatt,
     >   maxrps,lpsmax,fd0,itail,phasei,phasef)
      use chil_module
      include 'par.f'
      include 'par.pos'
      parameter( maxl=2*ltmax+1)
c      implicit real*8 (a-h,o-z)
      real*8 xgz,wgz,pl
      real rmesh
      complex phasei(kmax),phasef(kmax)
      common/meshrr/ meshr,rmesh(maxr,3)
      real chii(meshr,nqmi),gki(kmax),
     >    chif(meshr,nqmf),gkf(kmax)
      dimension minchii(nqmi),minchif(nqmf),npk(nchm+1)
      dimension fun(maxr), temp(maxr), psii(maxr), psif(maxr)
     >   ,ovlpi(kmax), ovlpf(kmax), temp2(maxr), temp3(maxr)
     >   ,tail(kmax,kmax),u(maxr) !, dwpot(maxr,nchan)
      real vmatt(kmax,kmax,0:1),
     >  vdon(nchan,nchan,0:1) !,Vvd(npk(nchm+1)-1,npk(nchm+1))
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)
      real zasym
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),maxpsinb(nnmax,0:lnabmax)
      common/matchph/rphase(kmax,nchan)
      common/numericalpotential/ numericalv, lstoppos
      logical numericalv
*      real chitemp(maxr,kmax) !, rpowpos1(maxr), rpowpos2(maxr)
      allocatable :: chitemp(:,:)
      logical pos, second, packed
      common/smallr/formcut,regcut,expcut,fast,match,analyticd,packed
      logical fast, match, analyticd
      real*8 xin(maxr),yin(maxr),xout(maxr),yout(maxr)
      data pi/3.1415926535897932384626434E0/
      real*8 factrl, sqrfct,hat
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl)
c!     >  , Qlfactor(0:ltmax)
      COMMON /gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz 
      dimension  fd0(maxr,0:lpsmax) ! fd0(maxr,0:lpsmax,npsch,npsch)
      real c1,c2,c3,c3tmp,cof6j,c1tmp,cgc0,c2tmp2,cgcigor,c2tmp
      character ch*1
      logical posi,posf
      ch(n) = char(mod(n,75) + ichar('0'))
      if(lg.gt.lstoppos) return
!      RETURN noPsPs
    
      rnorm= 2./pi    
      minfun = 1
      maxfun = min(maxpsii,maxpsif)
      second=.false.                ! make it as input
       ifdstop=maxfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!          nr = meshr
!          maxfun=nr
!          posi = positron(nia,lia,nposi)
!
!          maxflr=maxpsii; ni1 = nposi
!          r1=rmesh(maxflr-50,1); r2=rmesh(maxflr,1)
!          fl1=psii(maxflr-50); fl2=psii(maxflr)
!          alambdai=1./(r2-r1)*Log((fl1*r2**ni1)/(fl2*r1**ni1))
!          cfli=fl1*r1**(-ni1)*exp(alambdai*r1)
!
!          do i=maxflr+1,nr
!             rr=rmesh(i,1)
!             psii(i)= cfli*rr**ni1*exp(-alambdai*rr)
!          enddo
!!        WRITE(17,'(a1,3i8,2e20.10)')'#',ni1,maxflr,nr,alambdai,cfli
!!        WRITE(17,*) 
!
!        posf = positron(nfa,lfa,nposf)
!
!          maxflr=maxpsif; ni1 = nposf
!          r1=rmesh(maxflr-50,1); r2=rmesh(maxflr,1)
!          fl1=psif(maxflr-50); fl2=psif(maxflr)
!         alambdaf=1./(r2-r1)*Log((fl1*r2**ni1)/(fl2*r1**ni1))
!         cflf=fl1*r1**(-ni1)*exp(alambdaf*r1)
!          do i=maxflr+1,nr
!             rr=rmesh(i,1)
!             psif(i)= cflf*rr**ni1*exp(-alambdaf*rr)
!          enddo
!!        WRITE(17,'(a1,3i8,2e20.10)')'#',ni1,maxflr,nr,alambdaf,cflf
!!        WRITE(17,*) 
!!        do i=1,nr
!!         rr=gridr(i,1)
!!         WRITE(17,'(i5,3e20.10)') i,rr,psii(i),psif(i)
!!        enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

      do i = minfun, maxfun
         fun(i) = psii(i) * psif(i) * rmesh(i,3)
      end do 
      do i = 1,meshr
         if(rmesh(i,1).le.(8.+4.*nznuc)) ifdstop = i
      end do 
      ifdstop=min(ifdstop,maxfun)
      IF(ifdstop.gt.maxfun) THEN
       PRINT*, '',ifdstop,maxfun
       STOP'check ifdstop'
      ENDIF
 
      do i = 1, maxr
         temp3(i) = 0.0
      enddo
      mini = maxr
      maxi = 1
      ltmin=1000000  
      ctemp = 0.0

!       IF(((-1)**(lia+lfa).eq.1).or.((-1)**(li+lf).eq.1)) RETURN

      do 10 ilt = -lia, lia, 2


         temp(:)=0.0
         lt = lfa + ilt
!         if((-1)**lt.eq.1) go to 10 !.and.nznuc.eq.2) go to 10                          
         if(lt.lt.abs(li-lf).or.lt.gt.abs(li+lf)) go to 10 
         if (lt.lt.0.or.lt.gt.ltmax) go to 10
         call cleb(2*lia,2*lt,2*lfa,0,0,0,c1)
         c1tmp = cgc0(float(lia),float(lt),float(lfa))
         if (abs((c1-c1tmp)/(c1tmp+1e-20)).gt.1e-3) then
            print*,'CGCs 1 do not agree:',c1, c1tmp,lia,lt,lfa
            stop 'CGCs 1 do not agree'
         endif 
         call cleb(2*li,2*lf,2*lt,0,0,0,c2)
         c2tmp = cgc0(float(li),float(lf),float(lt))         
         if (abs((c2-c2tmp)/(c2tmp+1e-20)).gt.1e-3) then
            c2tmp2 = cgcigor(2*li,2*lf,2*lt,0,0,0)
            print*,'CGCs 2 do not agree:',c2, c2tmp,c2tmp2,li,lf,lt
            c2 = c2tmp2
         endif 
         call rac7(2*li,2*lf,2*lia,2*lfa,2*lt,2*lg,c3)
         c3tmp = cof6j(float(li),float(lf),float(lt),float(lfa),
     >      float(lia),float(lg)) !*(-1)**(li+lf+lia+lfa)
         if (abs((c3-c3tmp)/(c3tmp+1e-6)).gt.1e-3) then
            print*,'WARNING: CJ6 and W do not agree in D:',c3, c3tmp
            c3 = c3tmp
         endif 
     
         c = c1 * c2 * c3
         if (abs(c).lt.1e-10) go to 10
         const = (-1)**(lg + lfa) * hat(li) *
     >      hat(lf) * hat(lia) / hat(lt) * c
 
         call form(fun,minfun,maxfun,rpow1(1,lt),rpow2(1,lt),
     >      minrp(lt),maxrp(lt),meshr,temp,i1,i2)

         const = 2.0 * const  !!  following factor must be added only for He
                               !!  and hepsdir must be changed.  * (1-(-1)**lt)   
              do i = i1, i2
               xin(i-i1+1) = dble(rmesh(i,1))
               yin(i-i1+1) = temp(i) * xin(i-i1+1) ** (lt+1)
               xout(i-i1+1) = dble(rmesh(i,1)) * 2d0
            enddo 
            if (i2-i1+1.gt.maxr) then
               print*,'I2,I1,LT:',i2,i1,lt
               stop 'problem in call to intrpl'
            endif 
            call intrpl(i2-i1+1,xin,yin,i2-i1+1,xout,yout)
            do i = i1, i2
               if (xout(i-i1+1).gt.xin(i2-i1+1))
     >            yout(i-i1+1) = yin(i2-i1+1)
               temp(i) = 2.0 * (yout(i-i1+1) / xout(i-i1+1)**(lt+1))
            enddo

c!       ELSE  ! below is for He+1 ion
             ilam=NINT((lt-abs(lia-lfa))*0.5)
!         write(17,*)
!         write(17,'(a1,a5,5i6)')'#','lam =',lt,nia,lia,nfa,lfa
!         do i =i1,i2
!            write(17,'(4e20.10)')rmesh(i,1)
!     >           ,temp(i)*((-1)**lt-1),fd0(i,ilam),rmesh(ifdstop,1)
!         enddo

         const = - const !Alisher's derivation, as the opposite is below  
           do i=i1,ifdstop-1 !maxfun-1 
           temp(i)=temp(i)*((-1)**lt - 1) + fd0(i,ilam)
           enddo
           do i=ifdstop,i2  !maxfun,i2
           temp(i)=temp(i)*((-1)**lt - 1) 
           enddo

!       write(17,'(i6,2e20.10)') i2,rmesh(i2,1), temp(i2)           
c*!  fd0(i,lt) is temp for He+1 ion, temp was for proton
c*!  in large distances (>maxfun) they must be the same (~ 1/r)!!!
c*!        ENDIF ! he1_ion

        do i = i1,i2
        temp3(i) = const * temp(i) + temp3(i)   
        enddo
         mini = min(mini,i1)
         maxi = max(maxi,i2)
!         if ((lt.eq.1.or.lt.eq.2).and.i2.eq.meshr) then
        if(i2.eq.meshr.and.lt.lt.ltmin.and.(-1)**lt.ne.1)then
            ctemp = rnorm * const * temp(i2)/rpow2(i2,lt)
            ltmin=lt
c$$$            print*,'lf,li,ctemp:',lf,li,(temp(i)/rpow2(i,lt),
c$$$     >         i=meshr-2,meshr)
         endif 
 10   continue

!      if (pos.and.maxi.eq.meshr) then
!         print*,'Writing temp3 to disk'
!         open(42,file='temp'//ch(lfa)//ch(lia)//ch(lt))
!         write(42,*) '# lt,ltmin,ctemp:',lt,ltmin,rmesh(i2,1),ctemp
!         do i = 1, meshr
!         write(42,*)rmesh(i,1),temp3(i),temp3(i)*rnorm/rpow2(i,ltmin)
!         enddo
!         close(42)
!      endif
!      if (ltmin.lt.10.and.maxi.eq.meshr)
!     >   print*,'pos,lfa,lia,ltmin,maxi,ctemp,v(i):',pos,lfa,lia,ltmin,
!     >   maxi,ctemp,(rnorm*temp3(i)/rpow2(i,ltmin),i=maxi-20,maxi,10)

       if (ltmin.lt.10.and.maxi.eq.meshr) then
         ctemp = rnorm * temp3(maxi) / rpow2(maxi,ltmin)
!      else
!         ctemp = 0.0
       endif  

C  As both CHII and CHIF contain the integration weights, we divide TEMP by them.
      
      do i =mini, maxi
         temp3(i) = rnorm * temp3(i) / rmesh(i,3)
      end do
      mini1 = mini

      tail(:,:) = 0.0
      if (itail.ne.0)
     >call maketail(itail,ctemp,chil(1,npk(nchi),2),minchil(npk(nchi),2)
     >   ,gki,phasei,li,nqmi,chil(1,npk(nchf),2),minchil(npk(nchf),2)
     >   ,gkf,phasef,lf,nqmf,nchf,nchi,ltmin,tail)
c$$$      call maketail(itail,ctemp,chii,minchii,gki,phasei,li,nqmi,
c$$$     >   chif,minchif,gkf,phasef,lf,nqmf,nchf,nchi,ltmin,tail)

      allocate(chitemp(maxr,kmax))
      do ki =1,nqmi         !in vdme: npk(nchi),npk(nchi+1)-1 
         minki = max(mini1, minchii(ki))
         do i = minki, maxi
            chitemp(i,ki) = temp3(i) * chii(i,ki)
c!         write(21,'(3e20.10)') gki(ki),rmesh(i,1), chitemp(i,ki)
         enddo
      enddo 



      do ki =1,nqmi     !in vdme: kistart, kistop
         qi = gki(ki)
!         if(qi.lt.0.) cycle
cr         kii = npk(nchi) + ki - 1  !vdme: ki-kistep
         minki = max(mini1, minchii(ki))
         kfstop = nqmf
         if (second) then
            if (ki.eq.1) then
               kfstop = nqmf
            else
               kfstop = 1
            endif
         endif 

         do 20 kf =1, kfstop  !in vdme: npkf(nchf), npkf(nchf+1)-1
            qf = gkf(kf)
!         if(qf.lt.0.) goto 20       
    
           mini = max(minki, minchif(kf))
            n = maxi - mini + 1 
            tmp = tail(kf,ki)

  
C The following If statement avoids a bug on the IBMs
          if (mini.lt.maxi) then
          tmp =tmp +
     >     dot_product(chif(mini:maxi,kf),chitemp(mini:maxi,ki))
           endif
  
        vmatt(kf,ki,0) = vmatt(kf,ki,0) + tmp 
        vmatt(kf,ki,1) = vmatt(kf,ki,1) + tmp 

 20      continue 
      end do 
      deallocate(chitemp)
C  Define the direct on shell V matrix used for analytic Born subtraction
C  in the cross section 
      ki = 1
      kf = 1
      tmp = tail(kf,ki)
      mini = max(mini1, minchii(ki), minchif(kf))
      minki = max(mini1, minchii(ki))
      do i = minki, maxi
         temp2(i) = temp3(i) * chii(i,ki)
      enddo 
C The following If statement avoids a bug on the IBMs
      if (mini.lt.maxi) then
         tmp =tmp +
     >      dot_product(chif(mini:maxi,kf),temp2(mini:maxi))
       endif 
      do ns = 0, 1
         vdon(nchf,nchi,ns) = vdon(nchf,nchi,ns) + tmp
         vdon(nchi,nchf,ns) = vdon(nchf,nchi,ns)
      enddo
      return 
      end
    

!     This subroutine calculates Q_L+A_L functions
!     used in poshevmat
       subroutine qlandal(gp,Amg,nqmi,gridk,nchi,Li,Qlp,Q0p,aq,ionsh,
     *           nqgen ,nqim ,xp,wp,imax,Alp)
      include 'par.f'
      include 'par.pos'
      implicit real*8 (a-h,o-z)       
      parameter (maxl=2*ltmax+1)
      common/gauszR/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      common /Zatom/ Z
      common/gausp/xgp(igpm),wgp(igpm),igp
      real zasym 
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common/numericalpotential/ numericalv, lstoppos
       real*8 Qlp(kmax,2*kmax*igp),Q0p(kmax,2*kmax*igp)
     >           ,gp(iglp),Amg(iglp),Alp(kmax,2*kmax*igp)
      real*8 aq(0:kmax),qq(kmax)
      real gridk(kmax,nchan)
      real*8 xp(2*kmax*igp+kmax),wp(2*kmax*igp)
      real*8 xx(20),yy(20)
   
       if(Li.gt.lstoppos+latop)  RETURN
       if(nqmi.eq.1) return
       knqm=nqmi
!!       nqm=-1 !  1  for on-shell!
       call qmesh(gridk(1,nchi),knqm,imax,xp,wp,aq,ionsh,nqgen,nqim)
          imaxp=2*nqgen*igp 
          mnqmi=nqmi

      ddp=1.e+08
      pdp=0.0065d0
      nd1=2
      intp=2*nd1+1
       nnmin=INT(Log(1.E-06*ddp)/pdp) !1000
       gpmin=gp(nnmin)
       gpmax=3600. !gp(iglp-5) 
       
       Q0p(:,:) = 0.d0;  Qlp(:,:) = 0.d0;
       Alp(:,:) = 0.d0 

       do i=1,imaxp
          pp=dble(xp(i))
!          if(pp.gt.40.) cycle
          pp2=pp*pp

!       write(17,*)

       do iq=1,nqmi  ! 1- for on-shell and born      
          qa=dble(gridk(iq,nchi))
          if(qa.lt.0.) cycle 
           qa2=qa*qa
           p2q2=pp2+qa2
           ppqa2=2.d0*pp*qa
           arg=p2q2/ppqa2
       call funleg(arg,Li,Q0p(iq,i),Qlp(iq,i))

            resA=0.d0

       IF(nznuc.eq.2) THEN
         do iz=1,igz
            qapp2=p2q2-ppqa2*xgz(iz)
!c             Akq1= (1-256./(16.+qapp2)**2.d0)/qapp2 ! for He only
             Akq=(32.+qapp2)/(256.+32.*qapp2+qapp2*qapp2)
             resA=resA+Akq*pl(Li,iz)
          enddo
        ENDIF
         
       IF(nznuc.gt.2) THEN
         do iz=1,igz
            qapp2=p2q2-ppqa2*xgz(iz)
!             Akq1= (1-256./(16.+qapp2)**2.d0)/qapp2 ! for He only
!!! Akq = \int e^{i|qa-pp|r} vcore_pos(r) dr 
!   Akq=positroncoreV(qapp) /8/pi**2
!! use veint(r) ? vcore_pos(r) = veint(r)
          nn=INT(Log(qapp2*ddp)/pdp)
          if(nn.lt.nd1+1.or.nn.gt.iglp-nd1) cycle

        IF(qapp2.gt.gpmin.and.qapp2.lt.gpmax) THEN
          yy(1:intp)=Amg(nn-nd1:nn+nd1)
          xx(1:intp)=gp(nn-nd1:nn+nd1)
          call intrpl(intp,xx,yy,1,qapp2,Akq) ! for both Mg and He
!         print*,'He-Akq&Qlp check:', qapp2, Akq1, Akq
!          intrpl(qapp2,) ! for Mg
        ELSEIF(qapp2.le.gpmin) THEN 
          Akq=Amg(nnmin)
        ELSEIF(qapp2.ge.gpmax) THEN 
          Akq=(Z-1.d0)/qapp2
        ENDIF
 
             resA=resA+Akq*pl(Li,iz)
             enddo
        ENDIF
 
!        write(17,'(4e20.10)') qa,pp,Qlp(iq,i),qa*pp*resA 

              Alp(iq,i)=qa*pp*resA   ! resA * or / (2*Li+1) ???

       enddo 
       enddo 
!       STOP'Ql are made'
        RETURN
        End


!!!                                   / 
!  This subroutine calculates I(R)=  | dr Phi_{\beta}(r)W_{lam}(r,R))Phi_{beta'}(r)
!!!                                  /
      subroutine hepsdir(gp,KJ,ve2stat,va2stat,nchm,npsch,nchps,
     >                   lpsmax,maxrps,nchii,nchif,formf0,fd0)
      include 'par.f'
      include 'par.pos'
      implicit real*8 (a-h,o-z)    
  
      parameter (maxl=2*ltmax+1)
      real gridr
      common /meshrr/ nr, gridr(nmaxr,3)
      real*8 dgridr(nmaxr,3)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      COMMON /gauszR/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      logical posi, posf,positron
      real psii, psif, ei,ef
      dimension temp(maxr),psii(maxr),psif(maxr),jdouble(40)
      real qmax,pmax, gridp(iglp,3) 
      dimension pg(iglp),wpg(iglp)
     >          ,fun(maxr),fung(iglp)
!       real, allocatable :: fd0(:,:,:,:)
      dimension ve2stat(maxr),va2stat(maxr),gp(iglp),ve2s_gp(iglp),
     >    va2s_gp(iglp),xx(5),yy(5),fd0(maxr,0:lpsmax,npsch,npsch) 
      real zasym 
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common/numericalpotential/ numericalv, lstoppos
       dimension nchps(nchm)
      integer, allocatable :: nchansi(:), nchansf(:)
! note that 0:3 will restrict using Ps states with l>3
! as max(lia+lfa)=6 and lam has 4 values: 0,2,4,6 saved as (0,1,2,3) 
! and with total Ps states = < 100... change if necessary  
! but it may reqire more kmp_stacksize
! we dont need so many Ps states anyway
      IF(KJ.gt.lstoppos) RETURN  
          do ip=1,iglp
          gp(ip)=1.e-8*exp(0.0065*ip)
          enddo

!!!        RETURN

!         If(.not. allocated(fd0))  allocate(fd0(maxr,0:3,npsch,npsch)) !STOP'fd0 is not allocated'
        dgridr(:,:)=1.d0*dble(gridr(:,:))
       qmax=6.0 
       ndoub=6
!       pmax=gridr(maxrps,1) 
       pmaxold=gridr(maxrps,1) 
       rmaxion=8.+6.*nznuc
       pmax=max(pmaxold,rmaxion+pmaxold/2)

       nmaxp=iglp 
        
       call grids(qmax,ndoub,pmax,gridp,nmaxp,mpg,jdouble,
     >   40,njdouble)

c!        print*,'gridp&gridr:', mpg,gridp(mpg,1),nr,dgridr(nr,1)
c        call gauleg(gridp(mpg,1),40.,pg,wpg,mpg)
c        if(mpg.gt.1000) mpg=1000
        do ip=1,mpg  
        pg(ip)=1.d0*dble(gridp(ip,1))
        wpg(ip)=1.d0*dble(gridp(ip,3))
        if(pg(ip).le.(rmaxion)) mpgmax=ip
        enddo
         mpgmax=mpg

         do i=1,nr
         if(dgridr(i,1).le.(rmaxion)) ifdstop = i
         enddo

         maxpsii = nr !maxval(maxf(:))
         maxgp=1
        do ip=1,iglp
        if(gp(ip).ge.dgridr(maxpsii-2,1)) cycle
        if(gp(ip).le.dgridr(1,1)) mingp=ip
         maxgp=ip
        enddo
         nd1=2
         intp=2*nd1+1
        gpmax=8.+5d0 * nznuc    
        gpmin=gp(mingp+2)
      print*,'intrpl:maxr->maxgp',maxpsii,dgridr(maxpsii,1),
     >         mingp,':',maxgp,gp(mingp),':',gp(maxgp)

        va2s_gp(:)=0.d0
      call intrpl(maxpsii,dgridr(1,1),va2stat,maxgp,gp,va2s_gp)
     
        ve2s_gp(:)=0.d0
      call intrpl(maxpsii,dgridr(1,1),ve2stat,maxgp,gp,ve2s_gp)

!       print*, 'intrpl: maxr -> maxgp is done, He-like Z=', nznuc
c!       STOP
!        do ir=1, maxpsii,20
!         rr=dgridr(ir,1)
!        write(68,'(4e20.10)')rr,ve2stat(ir),va2stat(ir),
!     >             -1./rr-(0.48+1/rr)*exp(-0.96*rr)
!        enddo
!        do ip=1, maxgp,20
!         pp=gp(ip)
!        write(69,'(4e20.10)')pp,ve2s_gp(ip),va2s_gp(ip),
!     >             -1./pp-(0.48+1./pp)*exp(-0.96*pp)
!        enddo
!       STOP
!!!
      ddp=1.e+08
      pdp=0.0065d0
       
       fd0(:,:,:,:)=0.d0
C Unroll the nchi/nchf two loops into one over nch, for OpenMP efficiency.
      if(.not.allocated(nchansi).and..not.allocated(nchansf)) then
         allocate(nchansi((nchif-nchii+1)*(nchm-nchii+1)))
         allocate(nchansf((nchif-nchii+1)*(nchm-nchii+1)))
      nch = 0
      do nchi = nchii, nchif
         do nchf = nchi, nchm
            nch = nch + 1
            nchansi(nch) = nchi
            nchansf(nch) = nchf
         enddo
      enddo
      nchansmax = nch
      endif
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& default(private)
C$OMP& shared(nr,nchm,KJ,mpg,pg,wpg,dgridr,igz,xgz,pl,wgz,nchii,nchif)
C$OMP& shared(nchps,fd0,gp,ve2s_gp,va2s_gp,ddp,pdp,iw19,nznuc)
C$OMP& shared(mpgmax,ifdstop,formf0,gpmax,gpmin,nd1,intp)
C$OMP& shared(nchansi,nchansf,nchansmax)
C$OMP& shared(gridr)
      do nch = 1, nchansmax
         nchi = nchansi(nch)
         nchf = nchansf(nch)
c$$$      DO nchi = nchii,nchif     !1,nchm 
        call getchinfo (nchi,Ni,KJ,psii,maxpsii,ei,lia,nia,Li)

        posi = positron(nia,lia,nposi)
        if (.not.posi) cycle
c!         print '(i4)', nchi
         npsi=nchps(nchi)
c$$$        do 10 nchf = nchi, nchm
c!         print '(i4,$)', nchf
        call getchinfo (nchf,Nf,KJ,psif,maxpsif,ef,lfa,nfa,Lf)
         posf = positron(nfa,lfa,nposf)
         if (.not.posf) cycle
         maxri=min(maxpsii,maxpsif)
         npsf=nchps(nchf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!          maxflr=maxpsii; ni1 = nposi
!          r1=gridr(maxflr-50,1); r2=gridr(maxflr,1)
!          fl1=psii(maxflr-50); fl2=psii(maxflr)
!          alambdai=1./(r2-r1)*Log((fl1*r2**ni1)/(fl2*r1**ni1))
!          cfli=fl1*r1**(-ni1)*exp(alambdai*r1)
!
!          do i=maxflr+1,nr
!             rr=gridr(i,1)
!             psii(i)= cfli*rr**ni1*exp(-alambdai*rr)
!          enddo
!!        WRITE(17,'(a1,3i8,2e20.10)')'#',ni1,maxflr,nr,alambdai,cfli
!!        WRITE(17,*) 
!!
!          maxflr=maxpsif; ni1 = nposf
!          r1=gridr(maxflr-50,1); r2=gridr(maxflr,1)
!          fl1=psif(maxflr-50); fl2=psif(maxflr)
!         alambdaf=1./(r2-r1)*Log((fl1*r2**ni1)/(fl2*r1**ni1))
!         cflf=fl1*r1**(-ni1)*exp(alambdaf*r1)
!         do i=maxflr+1,nr
!             rr=gridr(i,1)
!             psif(i)= cflf*rr**ni1*exp(-alambdaf*rr)
!          enddo
!        WRITE(17,'(a1,3i8,2e20.10)')'#',ni1,maxflr,nr,alambdaf,cflf
!        WRITE(17,*) 
!        do i=1,nr
!         rr=gridr(i,1)
!         WRITE(17,'(i5,3e20.10)') i,rr,psii(i),psif(i)
!        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         fun(:)=0.d0
         do i=1,maxri
          fun(i)=dble(psii(i)*psif(i)) 
         enddo



         mpg1=mpg
        do i=1,mpg
         if(pg(i).ge.dgridr(maxri,1)) cycle
          mpg1=i+1
        enddo    
        fung(:)=0.
        call intrpl(maxri,dgridr(1,1),fun,mpg1,pg,fung)
        do lam=abs(lfa-lia),lfa+lia,2 
           iclam= (-1)**lam - 1
!          if(iclam.eq.0) cycle !.and.nznuc.eq.2) cycle         
          ilam=NINT((lam-abs(lfa-lia))*0.5)
          cl22=(2.*lam+1)/2. 
          temp(:)=0.d0
        DO i=1,mpg 
           ResInt=0.0
           riz=pg(i)
       do j=1,mpg1  
       if(abs(fung(j)).le.1.e-12) cycle
           rjz=pg(j)     
           U2lam=0.d0
           Ue2lam=0.d0
           ULtest=0.d0
           UeLtest=0.d0

         rgt=max(riz,rjz/2.)
         rlt=min(riz,rjz/2.)

         if((rgt-rlt).lt.(8.+4.*nznuc)) then 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ll=lam
            r1=rgt
            r2=rlt
            X=r2/r1
            
            res0=0.0d0          ! Hydrogen: Znuc = 1
            res1=0.0d0          ! Alkali Atom:  Znuc > 1

            DO ig=1,igz            
               z=xgz(ig)
               rr2=r1*r1+r2*r2-2*r1*r2*z             
               rr=dsqrt(rr2)

               nn=INT(Log(rr*ddp)/pdp)
               if(rr.gt.gpmin.and.rr.lt.gpmax) then
                  yy(1:intp)=(-1)**lam * ve2s_gp(nn-nd1:nn+nd1)
     >                 + va2s_gp(nn-nd1:nn+nd1)
                  xx(1:intp)=gp(nn-nd1:nn+nd1)
                  call intrpl(intp,xx,yy,1,rr,y)
               else if(rr.le.gpmin) then
                  y=(nznuc-2-formf0*rr)*iclam 
               else if(rr.ge.gpmax) then
                  y=-1.*iclam
               endif

               res1=res1+y/rr*pl(lam,ig) !  
            END DO

    
            res=((2.*ll+1.)/2.)*res1
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subf=X**lam*iclam/rgt 

            U2lam= (res + subf)   

         else

            U2lam=0.d0          !X**lam*iclam

         endif 

                 ResInt=ResInt+(fung(j))*U2lam*wpg(j)
              enddo
              temp(i)= ResInt
           ENDDO
          do i=1,mpg
                if(abs(temp(i)).lt.1.e-15) temp(i) = 0.d0 
          enddo

           call intrpl(mpg,pg,temp,maxri,dgridr(1,1),
     >          fd0(1,ilam,npsi,npsf))  

          do i=1,maxri
       if(abs(fd0(i,ilam,npsi,npsf)).lt.1.e-15)
     >         fd0(i,ilam,npsi,npsf) = 0.d0 
          enddo


        enddo                   !lam
c$$$  10      continue
      ENDDO                     !nch
C$OMP END PARALLEL DO
      if(allocated(nchansi).and.allocated(nchansf)) then
       deallocate(nchansi,nchansf)
      endif
      RETURN
      END
 



!!!                                   / 
!  This subroutine calculates I(R)=  | dr Phi_{\beta}(r)W_{lam}(r,R))Phi_{beta'}(r)
!!!                                  /
      subroutine hepsdir2(gp,KJ,ve2stat,va2stat,nchm,npsch,nchps,
     >                   lpsmax,maxrps,nchii,nchif,formf0,fd0)
      include 'par.f'
      include 'par.pos'
      implicit real*8 (a-h,o-z)    
  
      parameter (maxl=2*ltmax+1)
      real gridr
      common /meshrr/ nr, gridr(nmaxr,3)
      real*8 dgridr(nmaxr,3)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      COMMON /gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      logical posi, posf,positron
      real psii, psif, ei,ef
      dimension temp(maxr),psii(maxr),psif(maxr),jdouble(40)
      real qmax,pmax, gridp(iglp,3) 
      dimension pg(iglp),wpg(iglp)
     >          ,fun(maxr),fung(iglp)
!       real, allocatable :: fd0(:,:,:,:)
      dimension ve2stat(maxr),va2stat(maxr),gp(iglp),ve2s_gp(iglp),
     >    va2s_gp(iglp),xx(5),yy(5),fd0(maxr,0:lpsmax,npsch,npsch) 
      real zasym 
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common/numericalpotential/ numericalv, lstoppos
       dimension nchps(nchm)
      integer, allocatable :: nchansi(:), nchansf(:)
! note that 0:3 will restrict using Ps states with l>3
! as max(lia+lfa)=6 and lam has 4 values: 0,2,4,6 saved as (0,1,2,3) 
! and with total Ps states = < 100... change if necessary  
! but it may reqire more kmp_stacksize
! we dont need so many Ps states anyway
      IF(KJ.gt.lstoppos) RETURN  
          do ip=1,iglp
          gp(ip)=1.e-8*exp(0.0065*ip)
          enddo

!!!        RETURN

!         If(.not. allocated(fd0))  allocate(fd0(maxr,0:3,npsch,npsch)) !STOP'fd0 is not allocated'
        dgridr(:,:)=1.d0*dble(gridr(:,:))
       qmax=4.0 
       ndoub=6
!       pmax=gridr(maxrps,1) 
       pmaxold=gridr(maxrps,1) 
       rmaxion=8.+5.*nznuc
       pmax=max(pmaxold,rmaxion+pmaxold/2)

       nmaxp=iglp 
        
       call grids(qmax,ndoub,pmax,gridp,nmaxp,mpg,jdouble,
     >   40,njdouble)

c!        print*,'gridp&gridr:', mpg,gridp(mpg,1),nr,dgridr(nr,1)
c        call gauleg(gridp(mpg,1),40.,pg,wpg,mpg)
c        if(mpg.gt.1000) mpg=1000
        do ip=1,mpg  
        pg(ip)=1.d0*dble(gridp(ip,1))
        wpg(ip)=1.d0*dble(gridp(ip,3))
        if(pg(ip).le.(8.+5.*nznuc)) mpgmax=ip
        enddo
         mpgmax=mpg

         do i=1,nr
         if(dgridr(i,1).le.(8.+5.*nznuc)) ifdstop = i
         enddo

         maxpsii = nr !maxval(maxf(:))
         maxgp=1
        do ip=1,iglp
        if(gp(ip).ge.dgridr(maxpsii-2,1)) cycle
        if(gp(ip).le.dgridr(1,1)) mingp=ip
         maxgp=ip
        enddo
         nd1=2
         intp=2*nd1+1
        gpmax=8+3*nznuc !   =gp(maxgp-2) 
        gpmin=gp(mingp+2)
      print*,'intrpl:maxr->maxgp',maxpsii,dgridr(maxpsii,1),
     >         mingp,':',maxgp,gp(mingp),':',gp(maxgp)

        va2s_gp(:)=0.d0
      call intrpl(maxpsii,dgridr(1,1),va2stat,maxgp,gp,va2s_gp)
     
        ve2s_gp(:)=0.d0
      call intrpl(maxpsii,dgridr(1,1),ve2stat,maxgp,gp,ve2s_gp)

!       print*, 'intrpl: maxr -> maxgp is done, He-like Z=', nznuc
c!       STOP
!        do ir=1, maxpsii,20
!         rr=dgridr(ir,1)
!        write(68,'(4e20.10)')rr,ve2stat(ir),va2stat(ir),
!     >             -1./rr-(0.48+1/rr)*exp(-0.96*rr)
!        enddo
!        do ip=1, maxgp,20
!         pp=gp(ip)
!        write(69,'(4e20.10)')pp,ve2s_gp(ip),va2s_gp(ip),
!     >             -1./pp-(0.48+1./pp)*exp(-0.96*pp)
!        enddo
!       STOP
!!!
      ddp=1.e+08
      pdp=0.0065d0
       
       fd0(:,:,:,:)=0.d0
C Unroll the nchi/nchf two loops into one over nch, for OpenMP efficiency.
      allocate(nchansi((nchif-nchii+1)*(nchm-nchii+1)))
      allocate(nchansf((nchif-nchii+1)*(nchm-nchii+1)))
      nch = 0
      do nchi = nchii, nchif
         do nchf = nchi, nchm
            nch = nch + 1
            nchansi(nch) = nchi
            nchansf(nch) = nchf
         enddo
      enddo
      nchansmax = nch

C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
C$OMP& default(private)
C$OMP& shared(nr,nchm,KJ,mpg,pg,wpg,dgridr,igz,xgz,pl,wgz,nchii,nchif)
C$OMP& shared(nchps,fd0,gp,ve2s_gp,va2s_gp,ddp,pdp,iw19,nznuc)
C$OMP& shared(mpgmax,ifdstop,formf0,gpmax,gpmin,nd1,intp)
C$OMP& shared(nchansi,nchansf,nchansmax)
      do nch = 1, nchansmax
         nchi = nchansi(nch)
         nchf = nchansf(nch)
c$$$      DO nchi = nchii,nchif     !1,nchm 
        call getchinfo (nchi,Ni,KJ,psii,maxpsii,ei,lia,nia,Li)

        posi = positron(nia,lia,nposi)
        if (.not.posi) cycle
c!         print '(i4)', nchi
         npsi=nchps(nchi)
c$$$        do 10 nchf = nchi, nchm
c!         print '(i4,$)', nchf
        call getchinfo (nchf,Nf,KJ,psif,maxpsif,ef,lfa,nfa,Lf)
         posf = positron(nfa,lfa,nposf)
         if (.not.posf) cycle
         maxri=min(maxpsii,maxpsif)
         npsf=nchps(nchf)

         fun(:)=0.d0
         do i=1,maxri
          fun(i)=dble(psii(i)*psif(i)) 
         enddo
         mpg1=mpg
        do i=1,mpg
         if(pg(i).ge.dgridr(maxri,1)) cycle
          mpg1=i+1
        enddo    
        fung(:)=0.
        call intrpl(maxri,dgridr(1,1),fun,mpg1,pg,fung)
        do lam=abs(lfa-lia),lfa+lia,2 
           iclam=1-(-1)**lam
          if(iclam.eq.0.and.nznuc.eq.2) cycle         
          ilam=NINT((lam-abs(lfa-lia))*0.5)
          cl22=(2.*lam+1)/2. 
          temp(:)=0.d0
        DO i=1,mpg 
           ResInt=0.0
           riz=pg(i)
       do j=1,mpg1  
       if(abs(fung(j)).le.1.e-12) cycle
           rjz=pg(j)     
           U2lam=0.d0
           Ue2lam=0.d0
           ULtest=0.d0
           UeLtest=0.d0

         rgt=max(riz,rjz/2.)
         rlt=min(riz,rjz/2.)

       if((rgt-rlt).lt.(8.+2.*nznuc)) then 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ll=lam
      r1=rgt
      r2=rlt
      X=r2/r1

      res0=0.0d0                ! Hydrogen: Znuc = 1
      res1=0.0d0                ! Alkali Atom:  Znuc > 1
      select case (ll)
      case (0)
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             

        nn=INT(Log(rr*ddp)/pdp)
        if(rr.gt.gpmin.and.rr.lt.gpmax) then
         yy(1:intp)=ve2s_gp(nn-nd1:nn+nd1)
     >              +(-1)**lam*va2s_gp(nn-nd1:nn+nd1)
         xx(1:intp)=gp(nn-nd1:nn+nd1)
         call intrpl(intp,xx,yy,1,rr,y)
        else if(rr.le.gpmin) then
         y=-formf0*rr*iclam 
        else if(rr.ge.gpmax) then
         y=-(nznuc-1.)*iclam
        endif

            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                                                       
            res1=res1+y*wgz(ig) ! plz=1.0d0  ! case (0) 
!            res0=res0+wgz(ig)              
         END DO

      case (1)
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             

        nn=INT(Log(rr*ddp)/pdp)
        if(rr.gt.gpmin.and.rr.lt.gpmax) then
         yy(1:intp)=ve2s_gp(nn-nd1:nn+nd1)
     >              +(-1)**lam*va2s_gp(nn-nd1:nn+nd1)
         xx(1:intp)=gp(nn-nd1:nn+nd1)
         call intrpl(intp,xx,yy,1,rr,y)
        else if(rr.le.gpmin) then
         y=-formf0*rr*iclam 
        else if(rr.ge.gpmax) then
         y=-(nznuc-1.)*iclam
        endif

            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz=zpl             ! case (1)                        
            res1=res1+y*plz*wgz(ig)
!            res0=res0+plz*wgz(ig)
         END DO
      case (2)
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             

        nn=INT(Log(rr*ddp)/pdp)
        if(rr.gt.gpmin.and.rr.lt.gpmax) then
         yy(1:intp)=ve2s_gp(nn-nd1:nn+nd1)
     >              +(-1)**lam*va2s_gp(nn-nd1:nn+nd1)
         xx(1:intp)=gp(nn-nd1:nn+nd1)
         call intrpl(intp,xx,yy,1,rr,y)
        else if(rr.le.gpmin) then
         y=-formf0*rr*iclam 
        else if(rr.ge.gpmax) then
         y=-(nznuc-1.)*iclam
        endif

            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz=0.5d0*(3.0d0*zpl*zpl-1.0d0) ! case (2)                        
            res1=res1+y*plz*wgz(ig)
!            res0=res0+plz*wgz(ig)
         END DO   
      case (3)         
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             

        nn=INT(Log(rr*ddp)/pdp)
        if(rr.gt.gpmin.and.rr.lt.gpmax) then
         yy(1:intp)=ve2s_gp(nn-nd1:nn+nd1)
     >              +(-1)**lam*va2s_gp(nn-nd1:nn+nd1)
         xx(1:intp)=gp(nn-nd1:nn+nd1)
         call intrpl(intp,xx,yy,1,rr,y)
        else if(rr.le.gpmin) then
         y=-formf0*rr*iclam 
        else if(rr.ge.gpmax) then
         y=-(nznuc-1.)*iclam
        endif

            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz=0.5d0*(5.0d0*zpl*zpl-3.0d0)*zpl ! case (3)                        
            res1=res1+y*plz*wgz(ig)
!            res0=res0+plz*wgz(ig)
         END DO          
      case default         
         DO ig=igz,1,-1            
            z=xgz(ig)
            rr=r1*(1.0D0+X*z)             

        nn=INT(Log(rr*ddp)/pdp)
        if(rr.gt.gpmin.and.rr.lt.gpmax) then
         yy(1:intp)=ve2s_gp(nn-nd1:nn+nd1)
     >              +(-1)**lam*va2s_gp(nn-nd1:nn+nd1)
         xx(1:intp)=gp(nn-nd1:nn+nd1)
         call intrpl(intp,xx,yy,1,rr,y)
        else if(rr.le.gpmin) then
         y=-formf0*rr*iclam 
        else if(rr.ge.gpmax) then
         y=-(nznuc-1.)*iclam
        endif

            zpl=0.5d0*X*(1.0d0-z*z)-z ! new argument of pollegs                        
            plz0 = 1.0d0; plz1 = zpl; ! case default             
            do ill = 2, ll
            plz=(dble(2*ill-1)*zpl*plz1 - dble(ill-1)*plz0)/dble(ill)
               plz0=plz1
               plz1=plz
            end do                                     
            res1=res1+y*plz*wgz(ig)
!            res0=res0+plz*wgz(ig)              
         END DO          
      end select            
      
!      res0=((2.*ll+1.)/2.)*coef1*res0
      res=((2.*ll+1.)/2.)*res1
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subf=X**lam*(nznuc-1)*iclam

         U2lam= (res + subf)   

        else

           U2lam=0.d0 !X**lam*iclam

        endif 

     	 ResInt=ResInt+(fung(j))*U2lam*wpg(j)/rgt
        enddo
        temp(i)=ResInt
        ENDDO
        call intrpl(mpg,pg,temp,maxri,dgridr(1,1),
     >       fd0(1,ilam,npsi,npsf))  

         enddo !lam
c$$$ 10      continue
        ENDDO !nch
C$OMP END PARALLEL DO
        deallocate(nchansi,nchansf)
        RETURN
        END
 



      DOUBLE PRECISION FUNCTION PLeg(L,X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C ****  CALCULATES LEGENDRE POLYNOMIAL, ORDER L, ARGUMENT X
C
      PLeg = 0.0D0
        IF((ABS(X).LE.1.0D0).AND.(L.GE.0)) THEN
        IF(L.GT.5) GO TO 200
          IF(L.LE.1) THEN
            IF(L.EQ.0) THEN
            PLeg = 1.0D0
            ELSE
            PLeg = X
            END IF
          ELSE
          X2 = X*X
            IF(L.LE.3) THEN
              IF(L.EQ.2) THEN
              PLeg = 1.5D0*(X2-0.33333333333333D0)
              ELSE
              PLeg = 2.5D0*X*(X2-0.6D0)
              END IF
            ELSE
       PLeg=4.375D0*((X2-0.857142857142857D0)*X2+0.857142857142857D-1)
           IF(L.EQ.4) RETURN
       PLeg=7.875D0*X*((X2-1.11111111111111D0)*X2+0.238095238095238D0)
            END IF
          END IF
          RETURN
C
C ****  EVALUATE THE LEGENDRE POLYNOMIAL USING A RECURSION FORMULA
C ****  IF ITS ORDER IS TOO LARGE FOR DIRECT EVALUATION.
C
  200   CONTINUE
        X2 = X*X
        P0=4.375D0*((X2-0.857142857142857D0)*X2+0.857142857142857D-1)
        P1=7.875D0*X*((X2-1.11111111111111D0)*X2+0.238095238095238D0)
        DO 100 I = 6,L
        YL = DBLE(I-1)/DBLE(I)
        P2 = X*P1*(1.0D0+YL) - YL*P0
        P0 = P1
        P1 = P2
  100   CONTINUE
        PLeg = P2
        ELSE
        WRITE(6,1000)L,X
        END IF
      RETURN
 1000 FORMAT(1H0,'ATTEMPT TO FIND P',I2,'(',F5.2,')')
      END



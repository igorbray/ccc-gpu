*     This program returns
*     positronium formation
*     matrix elements.
*     History:
*               10 Dec 99 composite mesh for p-integral of
*                         igp*2*nqmi points is introduced
*               20 Jan 00 memory rearrangements by Igor
*               10 Feb 00 optimization and serious rearragement
*               22 Mar 00 optimization as to f0z and f1z
*                  Jun 06 generilized to include numerical potentials
*
*     Important notes:
*      
*               - with the composite mesh calculation of Q_L
*                 became timeconsuming; when channel information
*                 and k-grid are ready the composite mesh can be
*                 formed for each channel and k-grid; then Q_L
*                 can be calculated before Vmat and saved for use
*                 at different E. (J must probably be fixed, alpha>0)
*     
*               - the most time consuming parts are E-independant
*                 (this can be used provided alphas are fixed).      
*               -              -//-                 J-independant.
*      
* in order parallelization to work posvmat should not have local
* common blocks; all common blocks come from outside; posvmat
* should not try renew elements of these common blocks either.      
      subroutine posVmat(nqmi,lia,li,nia,gki,ni,nchi,nqmf,lfa,lf,
     >   nposf,nfa,gkf,nchf,lg,npk,etot,nqmfmax,vmatt,lm)
      use ql_module, only: index_ql, qlgmax
      use apar, only: alkali
      use ftps_module, only: interpol ! , error
      include 'par.f'
      include 'par.pos'
      parameter (maxl=8*ltmax)
      implicit real*8 (a-h,o-z)

      common/gausp/xgp(igpm),wgp(igpm),igp
      common /laguercoeffs/ cknd(ncmax,ncmax,0:lnabmax),
     $    rlambda(2,0:lnabmax), npsstates(2,0:lnabmax)      
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      common/gauszR/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      common/gausz/xgzR(igzm),wgzR(igzm),plR(0:maxl,igzm),igzR 
       common/numericalpotential/ numericalv, lstoppos
      logical numericalv
      dimension fpqb(2*nqmi*igpm+nqmi,0:lmax)!ltmax+lamax)      
      dimension x(ipm),ylam(0:lm,ipm),y2lam(0:lm,ipm)
      real gki,gkf,etot,vmatt
      dimension gki(kmax),gkf(kmax),npk(nchan+1),vmatt(nqmfmax,nqmfmax)!vmatt(kmax,kmax,0:1)
      dimension pc0(0:ncmax-1,ncmax,0:lamax),
     >   pc1(0:ncmax-1,ncmax,0:lamax)
      dimension a(ncmax),f0(0:ncmax-1),f1(0:ncmax-1)
      dimension icalam(0:lm),res0(0:lm)
      dimension y(ipm),y2(ipm)
      real fJ,fLlb,fLla,flb,fla,flb1,flb2,fla1,fla2,fl1,fl2,flam,
     >     w3j1,w3j2,w3j3,w3j4,w12j,cof3j,cof12j

      dimension aq(0:nqmi),qq(nqmi)
      dimension xp(2*nqmi*igpm+nqmi),wp(2*nqmi*igpm+nqmi),
     >   xppl(2*nqmi*igpm+nqmi,0:2*lm+2)
      dimension Qlp(2*nqmi*igpm,nqmi),Q0p(nqmi,2*nqmi*igpm)
      dimension f0z(igzm),f1z(igzm,2*nqmi*igpm+nqmi)
     >   ,f1zC(igzm,2*nqmi*igpm+nqmi)
      
C     ANDREY: new variables -----------------------------------
C     dimension dQlp(nqmi,2*nqmi*igpm), dQlp2(nqmi, 2*nqmi*igpm)
      dimension dQL1(2*nqmi*igp)
      logical exists, Mitroy, found
      
      common/smallr/ formcut,regcut,expcut,fast,match,analyticd,boxps
      logical fast, match, speed, analyticd, boxps

      real psinb, enpsinb
      integer istoppsinb
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >     psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      
      real*8 lg2
!      real*8, dimension (2*nqmi*igpm+nqmi) :: xp2

      character(len=20) :: fname*20
      character fname1*5 
      character(len=1), dimension (0:9) :: ch      
      data ch /'0','1','2','3','4','5','6','7','8','9'/
C      real qa
C     ANDREY: end my variables ---------------------------------
      real tmp
      dimension dQ(nqmi,2*nqmi*igpm)
      common/cheb/xc(ich),wc(ich)
      common/funQl/Qlarray(0:ltmax,ich)
!Charlie Variables for Charged targets
      real zasym
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      dimension coulm(nqmi,2*nqmi*igpm),subc(nqmi,2*nqmi*igpm),!Cfactor(nqmi),
     >   fpqbC(2*nqmi*igpm+nqmi,0:lmax)!ltmax+lamax)
      dimension res2aC(1:2,1:nqmi),res2bC(1:2,1:nqmi)
      real, dimension (0:lfa,0:lia,0:(lfa+lia),0:(lfa+lia),
     >   max0(lf,li,iabs(lf-lfa-lia),iabs(li-lfa-li)):min0(
     >   lf+lfa+lia,li+lfa+lia)) :: w12jtest

      if (lg.gt.lstoppos) return
!      print*,'nchi,nqmi,nchf,nqmf:',nchi,nqmi,nchf,nqmf
c$$$      if(nqmi.eq.1) return      ! no need to calculate Born separately
     
      !PRINT*,'ENTERED POSVMAT with nqmi=',nqmi
      !STOP

      pi = acos(-1d0)

      if (alkali) then         
         if (index_ql(li).eq.0) then
          ! to avoid this block max(li) has to be used in main
            print*, ' posvmat: no QL is provided for li = ',li
C            go to 444
            stop ' STOPPED IN POSVMAT'
         end if
         qmax = 10. ** qlgmax   ! from ql_module
      end if

      
      J=lg                    ! J_total
      Llb=lf                  ! L_proton
      Lla=li                  ! L_positron
      lb=lfa                  ! positronium orb. mom. (OM)
      la=lia                  ! target OM
c$$$      print*,'1,J,Llb,Lla,lb,la:',J,Llb,Lla,lb,la
      fJ=float(J)
      fLlb=float(Llb)
      fLla=float(Lla)
      flb=float(lb)
      fla=float(la)
      nb=nposf                ! positronium princ. q. num. (PQN)
      na=nia                  ! atomic PQN
      nbi = nfa               ! used to get w. f. (see gnlp)
      En = dble(etot)/2.0d0   ! total three-particle energy (atomic units)
C     ANDREY: START: -----------------------------------------
      nqmi2 = 2 * nqmi
      lg2 = Log(2.)
      Mitroy=.false.    
C      Mitroy=.true.  
      If (Mitroy) then        
         EAEB = enpsinb(na,la) + enpsinb(nbi,lb) 
         EnM = dble(etot-EAEB)/2.0d0   ! atomic units are used
      end if
C     ANDREY: END: -------------------------------------------      

*     this is a max possible value for lam(bda).
*     look lam in transfer routine. Moved to vmat.f
*     lm=min0(Llb+lb+la,Lla+lb+la)

* iphase below is a result of division posVmat by i**(Lla-Llb).
* Namely this factor is omitted from direct channel matrix
* elements in dirVmat.
      
      iphase=lb-la -(Lla-Llb)
c$$$      print*,'2,J,Llb,Lla,lb,la:',J,Llb,Lla,lb,la
      ic0=(-1)**(iphase/2+J+Llb)*(2*lb+1)*(2*la+1)
      c0=dble(ic0)*hat(Llb)*hat(Lla)*sqrfct(2*lb)*sqrfct(2*la)/pi

      pi8=8.0d0*pi*pi
      
      alfa=1.d0
      alfa2=alfa*alfa
      if (zasym.eq.1.0) then !HeII
      bohr1=0.5d0
      else  !H
      bohr1=1.d0
      endif
      bohr2=2.d0
      Nl1=npsstates(1,la)
      Nl2=npsstates(2,lb); ! print*, 'posvmat, Nl2: ', Nl2
      rlam1=rlambda(1,la)
      rlam2=rlambda(2,lb)
      brp1=4d0/rlam1/rlam1
      brp2=4d0/rlam2/rlam2
      bsp1=(8.d0/rlam1)**(la+2)*factrl(la)*0.25d0
      bsp2=(8.d0/rlam2)**(lb+2)*factrl(lb)*0.25d0
      
      if (.not.alkali) then
         if(Nl1.ne.0) then      ! for pos+hydrogen problem
            do k=1,Nl1
               a(k)=cknd(k,na,la)
c            print*,'a(',k,')=',a(k)
            end do
c         call coef0(la,Nl1,a,f0)
c            call coef1(la,Nl1,a,f1)
            if (.not.numericalv.and..not.interpol) then
            select case(la)
            case(0)
               call coefsS(la,Nl1,a,f0,f1)
            case(1)
               call coefsP(la,Nl1,a,f0,f1)
            case(2)
               call coefsD(la,Nl1,a,f0,f1)
            case(3)
               call coefsF(la,Nl1,a,f0,f1)
            case(4)
               call coefsG(la,Nl1,a,f0,f1)
            case(5)
               call coefsH(la,Nl1,a,f0,f1)
            case(6)
               call coefsI(la,Nl1,a,f0,f1)
            case(7)
               call newcoefs7(la,Nl1,a,f0,f1)
            case(8)
               call newcoefs8(la,Nl1,a,f0,f1)
            case(9)
               call newcoefs9(la,Nl1,a,f0,f1)
            case(10)
               call newcoefs10(la,Nl1,a,f0,f1)
            case(11)
               call newcoefs11(la,Nl1,a,f0,f1)
            case(12)
               call newcoefs12(la,Nl1,a,f0,f1)
            case(13)
               call newcoefs13(la,Nl1,a,f0,f1)
            case(14)
               call newcoefs14(la,Nl1,a,f0,f1)
            case(15)
               call newcoefs15(la,Nl1,a,f0,f1)
            case(16)
               call newcoefs16(la,Nl1,a,f0,f1)
            case(17)
               call newcoefs17(la,Nl1,a,f0,f1)
            case(18)
               call newcoefs18(la,Nl1,a,f0,f1)
            case(19)
               call newcoefs19(la,Nl1,a,f0,f1)
            case(20)
               call newcoefs20(la,Nl1,a,f0,f1)
            case default
               stop 'coeffs are missing for la>20'
            end select        ! la               
            do k=0,Nl1-1
               pc0(k,na,la)=f0(k)
               pc1(k,na,la)=f1(k)
            end do
            endif
         endif                ! Nl1
      end if                    ! alkali
      
      if(Nl2.ne.0) then         ! this is for positronium
         do k=1,Nl2
            a(k)=cknd(k,nbi,lb)
         end do
c     call coef0(lb,Nl2,a,f0)
c     call coef1(lb,Nl2,a,f1)
         if (.not.numericalv) then
         select case(lb)
         case(0)
            call coefsS(lb,Nl2,a,f0,f1)
         case(1)
            call coefsP(lb,Nl2,a,f0,f1)
         case(2)
            call coefsD(lb,Nl2,a,f0,f1)
         case(3)
            call coefsF(lb,Nl2,a,f0,f1)
         case(4)
            call coefsG(lb,Nl2,a,f0,f1)
         case(5)
            call coefsH(lb,Nl2,a,f0,f1)
         case(6)
            call coefsI(lb,Nl2,a,f0,f1)
         case(7)
            call newcoefs7(lb,Nl2,a,f0,f1)
         case(8)
            call newcoefs8(lb,Nl2,a,f0,f1)
         case(9)
            call newcoefs9(lb,Nl2,a,f0,f1)
         case(10)
            call newcoefs10(lb,Nl2,a,f0,f1)
         case(11)
            call newcoefs11(lb,Nl2,a,f0,f1)
         case(12)
            call newcoefs12(lb,Nl2,a,f0,f1)
         case(13)
            call newcoefs13(lb,Nl2,a,f0,f1)
         case(14)
            call newcoefs14(lb,Nl2,a,f0,f1)
         case(15)
            call newcoefs15(lb,Nl2,a,f0,f1)
         case(16)
            call newcoefs16(lb,Nl2,a,f0,f1)
         case(17)
            call newcoefs17(lb,Nl2,a,f0,f1)
         case(18)
            call newcoefs18(lb,Nl2,a,f0,f1)
         case(19)
            call newcoefs19(lb,Nl2,a,f0,f1)
         case(20)
            call newcoefs20(lb,Nl2,a,f0,f1)
         case default
            stop 'coeffs are missing for lb>20'
         end select             ! lb
c     call coefs(lb,Nl2,a,f0,f1)
         do k=0,Nl2-1
            pc0(k,nbi,lb)=f0(k)
            pc1(k,nbi,lb)=f1(k)
         end do
         endif
      endif                     ! Nl2
        
      if(gki(1).gt.0.d0) then
*     putting on-shell point to its place

         ionsh=1
         do i=2,nqmi
            if(gki(i).lt.gki(1)) then
               qq(i-1)=gki(i)
               ionsh=i
            else
               qq(i)=gki(i)
            endif
         end do
         qq(ionsh)=gki(1)

         do ji=1,nqmi
            j0=nqmi2*igp+ji
            xp(j0)=qq(ji)
         end do
         nqgen=nqmi
         nqim=nqmi
         imax=nqmi2*igp+nqmi
      else
*     on-shell point is negative (no on-shell point)

C     ANDREY: ERROR when  gki(1) < 0 and nqmi=1
        ionsh=-1  
         if (nqmi.eq.1) then
            nqim=1
            nqgen=nqim
            imax=2*nqim*igp+nqim
         else 
            do i=2,nqmi
               qq(i-1)=gki(i)
            end do
            nqim=nqmi-1
            nqmi2=nqim*2
            nqgen=nqim
            imax=2*nqim*igp+nqim
         do ji=1,nqim
            j0=2*nqim*igp+ji
            xp(j0)=qq(ji)
         end do
          end if
            
      endif      

*     end points
      aq(0)=0.d0
      do i=1,nqgen-1
         aq(i)=(qq(i+1)+qq(i))*0.5d0
         i1 = i-1             
*     forming composite mesh
         do ji=1,igp            
            j1=2*i1*igp+ji
            xp(j1)=aq(i1)+ (qq(i)-aq(i1)) * xgp(ji)
            wp(j1)=(qq(i)-aq(i1)) * wgp(ji)
            j2=(2*i-1)*igp+ji
            xp(j2)=aq(i) - (aq(i)-qq(i)) * xgp(igp+1-ji)
            wp(j2)=(aq(i)-qq(i)) * wgp(igp+1-ji)
         end do
      end do
*     last 2 intervals
      nqgen1 =  nqgen -1
      do ji=1,igp
         j1=2*nqgen1*igp+ji
         xp(j1)=aq(nqgen1)+(qq(nqgen)-aq(nqgen1))*xgp(ji)
         wp(j1)=(qq(i)-aq(i-1))*wgp(ji)
*     last interval
         j2=(2*nqgen-1)*igp+ji
         xp(j2)=qq(nqgen)/xgp(igp+1-ji)
C     wp(j2)=qq(nqgen)*wgp(igp+1-ji)/xgp(igp+1-ji)/xgp(igp+1-ji)
         wp(j2)=xp(j2)*wgp(igp+1-ji)/xgp(igp+1-ji)
      end do
*     composite mesh is ready
      

C     Calculate Q0p(iqa,i) & Qlp(i,iqa) for integration in (44).
c     Very slow for direct integration so interpolation is used.
c     dQls are calculated in main.f     
!      xp2(:) = xp(:)*xp(:) ! resulted in an overflow for large LG and unnatural parity
      if (alkali) then
         iii = 2*nqgen*igp     
         do while (xp(iii).gt.qmax)
            iii = iii - 1        
         end do
      endif 
                 
      dQl1(:) = 0d0

      !print'("posvmat: test: 1",i6)', Lla
      
      do iqa=1, nqmi         
         qa = gki(iqa)
         if(qa.lt.0.0) cycle
c$$$         if (qa.eq.0.0) then
c$$$            print*,'iqa,nqmi,qa:',iqa,nqmi,qa
c$$$            stop 'cannot have qa=0'
c$$$         endif
                           ising=iqa
                           if(iqa.le.ionsh) ising=iqa-1
                           if(iqa.eq.1) ising=ionsh
                        if(gki(1).lt.0.) ising=iqa-1
                           ising1 = ising - 1

         qa2 = qa*qa         
         if (alkali) then       ! ANDREY
            ! if (iqa.eq.1) print'("posvmat: test: 2: iii = ", i6)', iii
            if (iii.ge.1) call get_QL(Lla,qa,iii,xp(1:iii),dQL1(1:iii))
            dQ(iqa,1:iii) = dQL1(1:iii) ! test
!           and/or when direct integration is needed
            do jj = iii+1, 2*nqgen*igp 
               print*
               print*, 'xp, qmax', xp(2*nqgen*igp), qmax
               print*, iii, 2*nqgen*igp
               stop 'posvmat: avoid additional calculation of dQl1'
               pp = xp(jj)      ! inserted
               call funleg2(Lla,qa,pp,ql2)
               dQl1(jj) = (qa*pp/pi8)*ql2              
            end do  
         endif                                  
         if (abs(zasym).gt.0.d0) then !Charged target (needs checking for alkalis)
            eta=zasym/qa
            if ((eta.gt.6.0d0).or.
     >           (eta.lt.-1000.d0)) then !where subtraction is inaccurate
               coulm(iqa,:)=0.d0
               subc(iqa,:)=0.d0
               res2aC(:,iqa)=0.d0
               res2bC(:,iqa)=0.d0
            else
               do i=1,2*nqgen*igp            
                  pp=xp(i)
!Wave function
                  call cwfnFull(pp,qa,Lla,eta,wavec)
                  coulm(iqa,i)=wavec
!Subtracting function
                  dp=pp-qa
                  call RegCharge(qa,dp,Lla,eta,subr,sub)
                  subc(iqa,i)=sub
        !print*,qa,pp,wavec
               enddo 
               call RegCharge(qa,aq(ising1)-qa,Lla,eta,res2a,subiii)
               call SubCharge(aq(ising1),qa,eta,Lla,res2aa)
               call RegCharge(qa,aq(ising)-qa,Lla,eta,res2b,subiiii)
               call SubCharge(aq(ising),qa,eta,Lla,res2bb)
               res2aC(1,iqa)=res2a
               res2aC(2,iqa)=res2aa
               res2bC(1,iqa)=res2b
               res2bC(2,iqa)=res2bb
            endif               !Coulomb accuracy
         else !here for all neutral targets, including alkalis
            do i=1,2*nqgen*igp         
c$$$               if (xp(i)*qa.ne.0.0) then
                  pp=xp(i)
                  pp2=pp * pp !xp2(i)
                  arg=(pp2+qa2)/(2.0*pp*qa)
                  call funleg(arg,Lla,Q0p(iqa,i),Qlp1)
                  Qlp(i,iqa) = Qlp1 + dQL1(i)
c$$$               else
c$$$                  print*,'xp(i)*qa=0.0:',i,xp(i),qa
c$$$               endif
            enddo
         endif
      enddo                     ! iqa


c$$$C     test ----------------------------------------------------------      
c$$$!$omp critical       
c$$$      write(fname1,'("Ql_",I2.2)') Lla     ! test
c$$$      inquire(file=fname1,exist=found)
c$$$      open (unit=99, file=fname1)
c$$$      if (found) go to 504
c$$$      print'(" posvmat: calculating  Q2_",I2.2,1x,$)',Lla
c$$$      print'(" nqmi = ",i6)', nqmi        
c$$$      do iqa=1, nqmi
c$$$         qa = gki(iqa)
c$$$         if(qa.lt.0.0) go to 505                        
c$$$!        if (iii.ge.1) call get_QL(Lla,qa,iii,xp(1:iii),dQL1(1:iii))
c$$$         do i = 1, iii, 10            
c$$$            call funleg2(lla, qa, xp(i), ql2)
c$$$            tmp = qa * xp(i)/pi8
c$$$            arg=(xp(i)*xp(i)+qa*qa)/(2.0*xp(i)*qa)
c$$$            call funleg(arg,Lla,Q0111,Qlp1)
c$$$            
c$$$            write (99,'(6(1x,e15.7))') qa, xp(i), Qlp1,
c$$$     1           dQ(iqa,i), tmp * ql2
c$$$         end do
c$$$         write (99,*)
c$$$         write (99,*)
c$$$ 505     continue
c$$$      enddo                     ! iqa
c$$$!      stop 'posvmat: test: see rename Q2_ll to Q2_ll_qcut'
c$$$ 504  continue
c$$$
c$$$      
c$$$      open (unit=99, file='ql.dat')      
c$$$      print'(" posvmat: calculating  ql.tmp_",I2.2,1x,$)',Lla
c$$$      print'(" nqmi = ",i6)', nqmi        
c$$$      do iqa=1, nqmi
c$$$         qa = gki(iqa)
c$$$         if(qa.lt.0.0) go to 508                        
c$$$!        if (iii.ge.1) call get_QL(Lla,qa,iii,xp(1:iii),dQL1(1:iii))
c$$$         dx = 0.01
c$$$         do xx = -5.0,3.0,0.025            
c$$$            pp = 10.0 ** xx           
c$$$            call funleg2(lla, qa, pp, ql2)
c$$$            tmp = qa * pp/pi8                        
c$$$            write (99,'(6(1x,e15.7))') qa, pp,
c$$$     1           tmp * ql2
c$$$         end do
c$$$         write (99,*)
c$$$         write (99,*)
c$$$ 508     continue
c$$$      enddo                     ! iqa
c$$$      stop 'posvmat: test: see rename Q2_ll to Q2_ll_qcut'
c$$$
c$$$      
c$$$!$omp end critical     
c$$$C     end test ------------------------------------------------------  

            
c$$$!$omp critical      
c$$$      if (alkali) then
c$$$         write(fname, 203) lLa  ! test
c$$$ 203     format('QL1_',I1)      ! test
c$$$      
c$$$         inquire(file=fname, exist=exists)      
c$$$         if (exists) go to 204  ! 
c$$$      
c$$$         open(unit=90, file = fname) !##############
c$$$         print*
c$$$         print*, 'see ', fname 
c$$$         print*, 'nqmi, 2*nqgen*igp: ', nqmi, 2*nqgen*igp
c$$$         print*
c$$$         indx = 0
c$$$      
c$$$         do iqa=1, nqmi
c$$$            write(90,*) '#                  INDEX: ', indx
c$$$            write(90,*) '#--------------------------------------------'         
c$$$            qa = gki(iqa)
c$$$            if(qa.lt.0.0) cycle
c$$$            qa2 = qa*qa                  
c$$$            call get_QL(Lla,qa,iii,xp(1:iii),dQL1(1:iii))                                             
c$$$         
c$$$            do i=1,2*nqgen*igp            
c$$$               pp=xp(i)
c$$$               pp2=xp2(i)
c$$$               arg=(pp2+qa2)/(2.0*pp*qa)            
c$$$               call funleg(arg,Lla,Q0p(iqa,i),Qlp1)         
c$$$               Qlp(i,iqa) = Qlp1 + dQL1(i)
c$$$               write(90,'(4(e13.4,2x))'), qa, pp, Qlp1, dQL1(i)
c$$$            enddo               ! i
c$$$            write(90,*)            
c$$$            write(90,*)            
c$$$            indx = indx + 1         
c$$$         enddo                  ! iqa
c$$$         print*, 'indx = ', indx
c$$$         close (90);
c$$$
c$$$ 204     continue
c$$$      end if
c$$$!$omp end critical
      
      
CCC   print'('' ready'')'
      

c$$$CCC   ANDREY: output for f0z:
c$$$!$omp critical      
c$$$      
c$$$      write(fname, 205) la,na,lb,nb                ! test
c$$$ 205  format('f0z_<',I1,',',I1, '|',I1,',',I1,'>') ! test
c$$$      
c$$$      inquire(file=fname, exist=exists)      
c$$$      if (exists) go to 206     ! 
c$$$      
c$$$      open(unit=90, file = fname) !##############
c$$$      print*
c$$$      print*, 'see ', fname 
c$$$          
c$$$      qbmin = 1.d-3; qbmax = 3.d+2;
c$$$      cf = (log(qbmax)-log(qbmin))/999.0
c$$$
c$$$      write(90,*) '# atom state: na,  la: ', na, la
c$$$      write(90,*) '#  pos state: nbi, lb: ', nbi,lb
c$$$      write(90,*) '# '
c$$$      do iqb=1, 1000
c$$$         qb = qbmin * exp(cf*(iqb-1))
c$$$         qb2 = qb*qb
c$$$         if (alkali) then
c$$$            if (oldftps) then
c$$$               call getftps(qb, na, la, resH1, resH0) ! interpolation                  
c$$$               call getftps(qb, nbi,lb, resP1, resP0) ! interpolation
c$$$            else
c$$$               call getftps3(qb, na, la, resH1, resH0) ! interpolation                  
c$$$               call getftps3(qb, nbi,lb, resP1, resP0) ! interpolation
c$$$            end if
c$$$         else
c$$$            call f0zpart(Nl1,rlam1,bohr1,na,la,na,qb2,resH0
c$$$     $           ,resH1,pc0,pc1)
c$$$            call f0zpart(Nl2,rlam2,bohr2,nb,lb,nbi,qb2,resP0,
c$$$     $           resP1, pc0,pc1)
c$$$         end if                                            
c$$$         write(90,'(5(e13.4,2x))') qb, resH1, resH0, resP1, resP0
c$$$      end do                    ! iqb
c$$$                 
c$$$ 206  continue      
c$$$!$omp end critical                  
                  
c$$$CCC   ANDREY: output for f0z:
c$$$!$omp critical      
c$$$
c$$$c      print*, 'lb, nb, nbi:', lb, nb, nbi
c$$$      fname='psi_l_'// ch(lb)//'_n_'//ch(nb)//'.dat'
c$$$      
c$$$      inquire(file=fname, exist=exists)      
c$$$      if (exists) go to 206     ! 
c$$$      
c$$$      open(unit=90, file = fname) !##############
c$$$      print*
c$$$      print*, ' posvmat: see ', fname 
c$$$          
c$$$      qbmin = 1.d-3; qbmax = 3.d+2;
c$$$      cf = (log(qbmax)-log(qbmin))/999.0
c$$$
c$$$      write(90,*) '#  pos state: nbi, lb: ', nbi,lb
c$$$      write(90,*) '# '
c$$$      
c$$$      do iqb=1, 1000
c$$$         qb = qbmin * exp(cf*(iqb-1))
c$$$         qb2 = qb*qb        
c$$$
c$$$         if (oldftps) then                               
c$$$            call getftps(qb, nbi,lb, resP1, resP0) ! interpolation
c$$$         else                               
c$$$            call getftps3(qb, nbi,lb, resP1, resP0) ! interpolation
c$$$         end if         
c$$$         
c$$$         call f0zpart(Nl2,rlam2,bohr2,nb,lb,nbi,qb2,rP0,
c$$$     $        rP1, pc0,pc1)                                          
c$$$         write(90,'(5(e13.4,2x))') qb, resP1, resP0, rP1, rP0
c$$$      end do                    ! iqb
c$$$                 
c$$$ 206  continue      
c$$$!$omp end critical        

c$$$CCC   ANDREY: output for f0z:
c$$$!$omp critical      
c$$$
c$$$      fname='psi_l_'// ch(la)//'_n_'//ch(na)//'.dat'      
c$$$      inquire(file=fname, exist=exists)      
c$$$      if (exists) go to 206; open(unit=90, file = fname) 
c$$$      print*, ' posvmat: see ', fname
c$$$C      print*, ' oldftps = ', oldftps
c$$$      
c$$$          
c$$$      qbmin = 1.d-3; qbmax = 3.d+2;
c$$$      cf = (log(qbmax)-log(qbmin))/999.0
c$$$
c$$$      write(90,*) '#  pos state: na, la: ', na,la
c$$$      write(90,*) '# '
c$$$      
c$$$      do iqb=1, 1000
c$$$         qb = qbmin * exp(cf*(iqb-1))
c$$$         qb2 = qb*qb        
c$$$
c$$$         if (oldftps) then                                           
c$$$            call getftps(qb, na, la, resH1, resH0) ! interpolation 
c$$$         else                                           
c$$$            call getftps3(qb, na, la, resH1, resH0) ! interpolation  
c$$$         end if         
c$$$c         call f0zpart(Nl1,rlam1,bohr1,na,la,na,qb2,rH0
c$$$c     $        ,rH1,pc0,pc1)
c$$$         
c$$$         write(90,'(5(e13.4,2x))') qb, resH1, resH0 ! , rH1, rH0
c$$$         
c$$$      end do                    ! iqb
c$$$                 
c$$$ 206  continue      
c$$$!$omp end critical        
!         WRITE(17,'(a2, 2i6)')'#H',na,la
!         WRITE(17,'(a2, 2i6)')'#P',nb,lb
!          WRITE(17,*) 
!        DO ipt=1,4000
!         pa2=0.0000001*ipt*ipt 
!         pb2=0.0000001*ipt*ipt 
!         call f0zpart(Nl1,rlam1,bohr1,na,la,na,pa2,resH0
!     $                ,resH1,pc0,pc1)
!          call f0zpart(Nl2,rlam2,bohr2,nb,lb,nbi,pb2,resP0,
!     $                 resP1, pc0,pc1)
!        WRITE(17,'(5e20.10)') pa2, resH0,resH1,resP0,resP1 
!        ENDDO

C-----------------Added by Ivan - w12j coefficient calculation--------
C All coefficient are: fla1, fla, fJ, flb, fla2, fLla, fLlb, flb1, flb2, fl2, flam, fl1

C The varying coefficients are:
C flb1(0,lb), flb2(lb, 0), fla1(0,la), fla2(la,0), fl1(iabs(lb1-la1),lb1+la1), fl2(iabs(lb2-la2),lb2+la2), flam(max0(iabs(Llb-l1),iabs(Lla-l2)),min0(Llb+l1,Lla+l2))

C The constant coefficienct are:
C fJ, fLlb, FLla, flb, fla

      do lb1=0,lb
         lb2=lb-lb1
         flb1=float(lb1)
         flb2=float(lb2)
         do la1=0,la
            la2=la-la1
            fla1=float(la1)
            fla2=float(la2)
            do l1=iabs(lb1-la1),lb1+la1
               fl1=float(l1)
               do l2=iabs(lb2-la2),lb2+la2
                  fl2=float(l2)
                  do lam=max0(iabs(Llb-l1),iabs(Lla-l2)),
     >               min0(Llb+l1,Lla+l2)
                     if (lam.gt.lmax) stop 'lam > lmax'
                     flam=float(lam)
                     w12jtest(lb1,la1,l1,l2,lam)=
     >                  cof12j(fla1,fla, fJ, flb,
     >                  fla2,fLla,fLlb,flb1,
     >                  flb2,fl2, flam,fl1)*
     >                  cof3j(fl2,flam,fLla,0.,0.,0.)*
     >                  cof3j(flam,fl1,fLlb,0.,0.,0.)
                  end do
               end do
            end do
         end do
      end do


      do iqb=1, nqmf 
         qb=dble(gkf(iqb));
         if(qb.lt.0.0) cycle
         qb2=qb*qb
         qb24En=qb2/4.-En
         !print*,iqb,qb,na,nb
         iq=0

         do iqa=1, nqmi 
            qa=dble(gki(iqa));  
            if(qa.lt.0.0) cycle
            qa2=qa*qa
            qbqa=qb*qa            
            qbqa2=2.d0*qbqa          

!        IF(abs(qa-qb).lt.0.001.or.abs(qa-qb/2.).lt.0.001) CYCLE 
!       IF(abs(qa-qb).le.0.001.and.qa.gt.2.and.lb.gt.2) THEN 
!         WRITE(17,'(a1,4i4,2e20.10)') '#', na,la, nb,lb, qa,qb
!          WRITE(17,*) 
!        ENDIF

            do iz=1,igz
               pb2=0.25d0*qb2 + qa2 - qbqa * xgz(iz)
               pa2=qb2 + qa2 - qbqa2 * xgz(iz)               
               
C     calculation of fourier transforms:
               !if (numericalv) then
c              if (alkali) then
               
               pa = sqrt(pa2)
               pb = sqrt(pb2)   ! for positronium
                  
c     - target states:              
               if (numericalv) then                 
                  call gnlp(bohr1,la,pa2,na,resH0,resH1)                                    
               else                                    
                  if ((.not.interpol).and.(.not.alkali)) then
                     !pause ' direct calculations'                   
                     call f0zpart(Nl1,rlam1,bohr1,na,la,na,pa2,resH0
     $                    ,resH1,pc0,pc1)
                     !if (alkali) STOP 'posVmat: 1'
                  else
                     !error = .true.
                     !print*,' *** pa, na, la, ERR :: ',pa,na,la,error
                     call getftps(pa, na, la, resH1, resH0)

                     !print*,pa, resH1, resH0
                     !STOP 'posVmat: 2'
c                    call getftps(pb, nbi,lb, resP1, resP0) 
                  end if 
               end if           ! numericalv               
               
c     - Ps states:
               if (numericalv) then                              
                  call gnlp(bohr2,lb,pb2,nbi,resP0,resP1)                     
               else             ! analytic formula
                  !if (interpol) then
                  call f0zpart(Nl2,rlam2,bohr2,nb,lb,nbi,pb2,resP0,
     $                    resP1, pc0,pc1)
                  !else         !interpolation
                  !   call getftps(pb, nbi,lb, resP1, resP0)
                  !end if       ! analytics 
               endif            !  numericalv
           

C-----ANDREY: mitroy simplification ---------------------------------                           
               if (Mitroy) then ! check if it works
                              ! for higher energy states
!                  EAEB = enpsinb(na,la) + enpsinb(nbi,lb)        
!                  EnM = dble(etot-EAEB)/2.0d0  
                  efactorM = qb2/4. - pa2/2. - EnM     
                  f0z(iz) = efactorM * resP1 * resH1      
               else           !  use Eq (34)
*     A general (off-shell, eigen/pseudostate) case           
                  efactor = qb24En+pb2 !qa2/2. - En   !  = qa2/2 + pa2/2 - En = qb2/4 + pb2 - En               

                  f0z(iz) = efactor*resP1*resH1-resP0*resH1-
     $                 resP1*resH0

!       IF(abs(qa-qb).le.0.001.and.qa.gt.2.and.lb.gt.2) 
!     >WRITE(17,'(i6,3e20.10)') iz,xgzR(iz), f0z(iz),plR(lb,iz)/wgzR(iz)

               end if
               
C-----ANDREY: end mitroy simplification ------------------------------
               
*     a pole dominance Faddeev approach; seems not to work at all. Why?
*     f0z(iz)=-efactor*resP1*resH1
            end do            !iz

            iq=iq+1
            
            if(iq.eq.1) then
               do i=1,imax
                 pp2 = xp(i)*xp(i)
!                  pp2 = xp2(i)
                  qbpp = qb*xp(i)
                  qbpp2 = 2.d0 * qbpp
                  bb = 0.25d0*qb2 + pp2
                  aa = qb2 + pp2        
                  do iz=1,igz
                     pb2=bb - qbpp * xgz(iz)
                     pa2 = aa - qbpp2 * xgz(iz)
!Charlie, Hamiltonian for Charged targets
                     efactorC=qb24En+pb2
                     pa = sqrt(pa2)
!     calculates w.f.                                         
                     pb = sqrt(pb2)                        
                     if (numericalv) then                           
                        call gnlp(bohr1,la,pa2,na,res0sk0,sk1) 
                        call gnlp(bohr2,lb,pb2,nbi,res0sk0,sk2)

                     else
                        if(interpol.or.alkali) then ! target
                           call getftps(pa, na, la, sk1, res0sk0)                           
c                          call getftps(pb,nbi,lb,sk2,res0sk0) 
                        elseif (abs(zasym).gt.0.d0) then ! Charged target
                        call f0zpart(Nl1,rlam1,bohr1,na,la,na,
     $                            pa2,sk11,sk1,pc0,pc1)
                        else    ! only for hydrogen
!                          call f0zpart(Nl1,rlam1,bohr1,na,la,na,
!     $                            pa2,res0sk0,sk1,pc0,pc1)                  
                           if(Nl1.eq.0) then
                              call geigen(bohr1,na,la,pa2,sk1)               
                           else
                              brap1=brp1*pa2+1.d0
                              x1=(brap1-2.d0)/brap1
                              sk1=pc1(Nl1-1,na,la)
                              do k=Nl1-2,0,-1
                                 sk1=sk1*x1+pc1(k,na,la)
                              end do
c$$$                              if (abs(brap1).gt.1d10) print*,
c$$$     >                           sk1,bsp1,brap1,la
                              sk1=sk1*bsp1/brap1**(la+2)
!                              sk1=sk1*exp(log(bsp1)-log(brap1)*(la+2))
                           end if ! Nl1                         
                        end if  ! interpol.or.alkali                          
c     this is for positronium
                        if (abs(zasym).gt.0.d0) then !Charged target
                       call f0zpart(Nl2,rlam2,bohr2,nb,lb,nbi,pb2,sk22,
     $                    sk2, pc0,pc1)
                        else
                        if(Nl2.eq.0) then
                           call geigen(bohr2,nb,lb, pb2,sk2)
                        else
                           brap2=brp2*pb2+1.d0
                           x2=(brap2-2.d0)/brap2
                           sk2=pc1(Nl2-1,nbi,lb)
                           do k=Nl2-2,0,-1
                              sk2=sk2*x2+pc1(k,nbi,lb)
                           end do
                           sk2=sk2*bsp2/brap2**(lb+2)
                        end if  ! Nl2
                        endif !Charge
                     end if     !numericalv                                                
                     if (abs(zasym).gt.0.d0) then !Charged target
                     f1zC(iz,i)=sk2*sk1*efactorC
     >          -dble(zasym+1d0)*sk2*sk11-sk22*sk1
                     else
                     f1zC(iz,i)=0.d0
                     endif
                     f1z(iz,i)=(zasym+1d0)*sk2*sk1 ! (Eq. 42)                     
                  end do        ! iz
               end do           ! i
            end if              !if(iq.eq.1) 
            
C-----ANDREY: calculate f1z for dQl -------------------------------            
            xppl(1:imax,0) = 1.0
            do ilam=1,lb+la+1 !2*lm+2
               xppl(1:imax,ilam) = xppl(1:imax,ilam-1)*xp(1:imax) 
            enddo

            do ilam=0,lm
               icalam(ilam)=0   ! or 1 if fpqb( : , lam) has been calculated
               res0(ilam)=0.
            enddo
            
            sumlb1=0.d0       ! Eq. (43)
            do lb1=0,lb
               lb2=lb-lb1
               flb1=float(lb1)
               flb2=float(lb2)
               clb=sqrfct(2*lb1)*sqrfct(2*lb2)*2.d0**lb1
               sumla1=0.d0
               do la1=0,la
                  la2=la-la1
                  lab2=lb2+la2
                  lab21 = lab2+1 ! ANDREY
                  qal2=qa**lab21
c$$$                  do i = 1, imax
c$$$                     xppl(i) = xp(i)**lab21
c$$$                  enddo
                  if (lab21.gt.lb+la+1) then !2*lm+2) then
                     print*,'lab21+1,la+lb+1:',lab21,la+lb+1 !2*lm+2
                     stop 'la+lb+1 too small'
                  endif
                  fla1=float(la1)
                  fla2=float(la2)
                  cla=1.d0/sqrfct(2*la1)/sqrfct(2*la2)
                  cla=cla*qb**(lb1+la1)
                  suml1=0.d0
                  do l1=iabs(lb1-la1),lb1+la1!,2 this and below for unnatural parity
                     fl1=float(l1)
                     w3j1=cof3j(fla1,flb1,fl1,0.,0.,0.)
                     if(w3j1.eq.0.) cycle                  
                     suml2=0.d0
                     do l2=iabs(lb2-la2),lb2+la2!,2
                        fl2=float(l2)
                        w3j2=cof3j(fla2,flb2,fl2,0.,0.,0.)
                        if(w3j2.eq.0.) cycle                     
                        reslam=0.d0
                        do lam=max0(iabs(Llb-l1),iabs(Lla-l2)),
     >                     min0(Llb+l1,Lla+l2)!,2
c$$$                           if (lam.gt.lmax) stop 'lam > lmax'
c$$$                           flam=float(lam)
c$$$                           w3j3=cof3j(flam,fl1,fLlb,0.,0.,0.)
c$$$                           w3j4=cof3j(fl2,flam,fLla,0.,0.,0.)
c$$$                           w3j34=dble(w3j3)*dble(w3j4)
c$$$                           if(w3j34.eq.0.d0) cycle                          
c$$$                           w12j=cof12j(fla1,fla, fJ,  flb,
c$$$     >                                 fla2,fLla,fLlb,flb1,
c$$$     >                                 flb2,fl2, flam,fl1)
c$$$                           wigner=w3j34*dble(w12j)
c$$$                           wigner=w3j34*dble(w12jtest(lb1,la1,
c$$$     >                        l1,l2,lam))
                           wigner = w12jtest(lb1,la1,l1,l2,lam)
                           if(wigner.eq.0.d0) cycle                         
                           if(icalam(lam).ne.1) then
                              
                              resz=0.d0 ! Eq.(37)
                              do iz=1,igz
                                 resz=resz + f0z(iz) *  pl(lam,iz)
                              end do
                              res0(lam)=resz                                  
                              if(iq.eq.1) then
*                                calculating f(p)----------------------------
                                 do i=1,imax !                        
                                    reszC=0.d0                       !   II
                                    resz=0.d0 ! same as Eq.(37) but for F
                                    do iz=1,igz !                        lam
                                       resz=resz+f1z(iz,i) * pl(lam,iz)
                                       reszC=reszC+f1zC(iz,i)*pl(lam,iz)
!Charlie: additional charge term
                                    end do
                                    fpqb(i,lam)=resz
                                    fpqbC(i,lam)=reszC
                                 end do
                              elseif(iq.eq.2) then
*     f(p) is ready
!      call zintegc(lam,x,y,y2,qb,lb,la,nb,na,nbi,pc1)
!      do ips=1,ipm
!      ylam(lam,ips)=y(ips)
!      y2lam(lam,ips)=y2(ips)
      !print*,x(ips),ylam(lam,ips),y2lam(lam,ips)  
!      enddo
                              endif
                              icalam(lam)=1
                           endif

*     sraightforward calculations (without using spline)
c     call transc(lam,lab2,res1,qb,qa,
c     >                  Lla,lb,la,nb,na,nbi,pc1)

c     ***** calculations of the p-integral using chebyshev mesh
c     open(33,file='fpold.dat')
!       call transc(lam,lab2,res1,qb,qa,Lla,lb,la,nb,na,nbi,pc1)
*     calculations using spline (10 to 20 times faster than trans)
!       y(:)=ylam(lam,:)
!       y2(:)=y2lam(lam,:)
!      call transpl(lam,lab2,res1,x,y,y2,qa,Lla,lm)
c      pause ' wrote to unit 10'
      !print*,' trans0 and trans1: ',res0(lam),res1,res1t
!        else
*     calculating p-integral using composite mesh
!                           if(gki(1).gt.0.d0) then ! Rav: subtraction
!                           method must be used for both open and closed
!                           channels. Some changes have been made to
!                           make it work. 
*     this is where singular point is
             if (abs(zasym).gt.0.d0) then


                           ising=iqa
                           if(iqa.le.ionsh) ising=iqa-1
                           if(iqa.eq.1) ising=ionsh
                        if(gki(1).lt.0.) ising=iqa-1
                           ising1 = ising - 1
                           
*     calculation of p integral                           
               eta=zasym/qa            
        
        i1=(2*ising-1)*igp
        i2=(2*ising-1)*igp+1
        q1=xp(i1)
        q2=xp(i2)

        Fqa =fpqbC(nqmi2*igp+ising,lam)*qa**(lab21) 
        Fqaq1=(fpqbC(i1,lam)+fpqb(i1,lam)*
     >          (qa**2.d0-q1**2.d0)/(2.d0*dble(zasym)))*q1**(lab21)
        Fqaq2=(fpqbC(i2,lam)+fpqb(i2,lam)*
     >          (qa**2.d0-q2**2.d0)/(2.d0*dble(zasym)))*q2**(lab21)
        FqaC=Fqa*2.0*pi!*exp(-pi*eta/2.0)/qa!2A*F*qa
        Fq1C=Fqaq1*2.0*pi
        Fq2C=Fqaq2*2.0*pi
        Fgrad=(Fq2C-Fq1C)/(q2-q1)

                           if(ising.ne.nqim) then
*     i.e. if singularity is not at last k-mesh point                              
                              res1 = 0.0d0
*     integrals coming before singularity
                              do i=1,2*ising1*igp
                                 pp=xp(i)
c                dp=pp-qa
c        call RegCharge(qa,dp,Lla,eta,subr,sub)
                          fpp=fpqbC(i,lam)+
     >            fpqb(i,lam)*(qa**2.d0-pp**2.d0)/(2.d0*dble(zasym))
                                fp=fpp*pp**(lab21+1)*coulm(iqa,i)!Coulomb wave
                                 res1 = res1 + wp(i) * fp
                              end do
                              
*     integrals with singularities
               res2=0.d0
               
                              do i=2*ising1*igp+1,(2*ising-1)*igp
                                 pp=xp(i)
        dp=pp-qa
c        call RegCharge(qa,dp,Lla,eta,subr,sub)
                fpp=fpqbC(i,lam)+
     >            fpqb(i,lam)*(qa**2.d0-pp**2.d0)/(2.d0*dble(zasym))
         
                                fp=fpp*pp**(lab21+1)*coulm(iqa,i)+
     >           2.d0*(Fgrad+FqaC/dp)*subc(iqa,i)/(pp+qa)!*exp(2.d0*pi*eta)/qa

                                 res2 = res2 + wp(i)*fp
                              end do                              

        res2a=2.d0*Fgrad*res2aC(2,iqa)+FqaC*res2aC(1,iqa)/dble(zasym)
        res2b=2.d0*Fgrad*res2bC(2,iqa)+FqaC*res2bC(1,iqa)/dble(zasym)

                              res3=0.d0
                              do i= (2*ising-1)*igp+1, 2*ising*igp
                                 pp=xp(i)
        dp=pp-qa
c        call RegCharge(qa,dp,Lla,eta,subr,sub)
                fpp=fpqbC(i,lam)+
     >            fpqb(i,lam)*(qa**2.d0-pp**2.d0)/(2.d0*dble(zasym))
         
                                fp=fpp*pp**(lab21+1)*coulm(iqa,i)+
     >                   2.d0*(Fgrad+FqaC/dp)*subc(iqa,i)/(pp+qa)
                                 res3 = res3 + wp(i)*fp
                              end do
                              
*     integrals coming after singularity
                              res4=0.d0
                              do i=2*ising*igp+1,nqmi2*igp
                                 pp=xp(i)
c        dp=pp-qa
c        call RegCharge(qa,dp,Lla,eta,subr,sub)
                        fpp=fpqbC(i,lam)+
     >            fpqb(i,lam)*(qa**2.d0-pp**2.d0)/(2.d0*dble(zasym))

                                fp=fpp*pp**(lab21+1)*coulm(iqa,i)!Coulomb wave
                                 res4 = res4 + wp(i) * fp
                              end do

!        print*,res1,res2,res3,res4,res2a,res2b
        res1=res1+res2+res3+res4
     >           +res2a-res2b
c        print*,res1,qa*qb*sqrt(2.d0)*res1/(2.d0*pi**3.d0)
c     >    ,'Current Idea'
        
                        else
*     i.e. if singularity is at last k-mesh point then

                              !nqmi1 = nqim-1
                              res1 = 0.0d0
                              
*     integrals coming before singularity
                              do i=1,2*(nqmi)*igp
                                 pp=xp(i)                              
                          fpp=fpqbC(i,lam)+
     >            fpqb(i,lam)*(qa**2.d0-pp**2.d0)/(2.d0*dble(zasym))
                                fp=fpp*pp**(lab21+1)*coulm(iqa,i)!Coulomb wave

                                 res1 = res1 + wp(i) * fp                                                                
                              end do
                           
*     integrals with singularities
                              Fqa = fpqbC(nqmi2*igp+nqim,lam)*qa**lab21
                              FqaC=Fqa*2.0*pi/qa!*exp(pi*eta/2.0)/qa
                                eta=zasym/qa

                              res2=0.d0
                              do i=2*(nqmi)*igp+1,(nqmi2)*igp
                                 pp=xp(i)
        dp=pp-qa
*        call RegCharge(qa,dp,Lla,eta,subr,sub)
                fpp=fpqbC(i,lam)+
     >            fpqb(i,lam)*(qa**2.d0-pp**2.d0)/(2.d0*dble(zasym))
         
                                fp=fpp*pp**(lab21+1)*coulm(iqa,i)+
     >                           (2.d0*FqaC/dp)*subc(iqa,i)/(qa+pp)!*exp(2.0*pi*eta)/qa

                                 res2 = res2 + wp(i)*fp                                 
                              end do
         dp=10E-6
        call RegCharge(qa,aq(nqim)-qa,Lla,eta,res2a,subiii)
        call RegCharge(qa,dp,Lla,eta,res2b,subiiii)

        res2a=res2a*FqaC/dble(zasym)
        res2b=res2b*FqaC/dble(zasym)

                                 
                              res1=res1+res2+res2a-res2b

                         endif
             else !neutral
                           ising=iqa
                           if(iqa.le.ionsh) ising=iqa-1
                           if(iqa.eq.1) ising=ionsh
                        if(gki(1).lt.0.) ising=iqa-1
                           ising1 = ising - 1
                           
*     calculation of p integral                           
                           
                           if(ising.ne.nqim) then
*     i.e. if singularity is not at last k-mesh point                              
                              
                              res1 = 0.0d0
*     integrals coming before singularity
                              if (imax.lt.2*ising1*igp) stop 'imax1'
                              do i=1,2*ising1*igp
                                 pp=xp(i)
!                                 fp=fpqb(i,lam)*pp**lab21*Qlp(i,iqa)
                                 fp=fpqb(i,lam)*xppl(i,lab21)*Qlp(i,iqa)
                                 res1 = res1 + wp(i) * fp
                              end do
                              
*     integrals with singularities
                              Fqa = fpqb(nqmi2*igp+ising,lam)*qa**lab21 
                              res2=0.d0
                              if (imax.lt.(2*ising-1)*igp) stop 'imax2'
                              do i=2*ising1*igp+1,(2*ising-1)*igp
                                 pp=xp(i)
!                                 fp=fpqb(i,lam)*pp**lab21*Qlp(i,iqa)
                                 fp=fpqb(i,lam)*xppl(i,lab21)*Qlp(i,iqa)
     $                                -Fqa * Q0p(iqa,i)
                                 res2 = res2 + wp(i)*fp
                                                                  
                              end do                              

                              res1b= - qa * Log(qa - aq(ising1)) - qa
     $                             * Log(aq(ising1)+qa) - aq(ising1)
     $                             * Log((aq(ising1) + qa)/ (-aq(ising1)
     $                             + qa))  

                              res3=0.d0
                              if (imax.lt.2*ising*igp) stop 'imax3'
                              do i= (2*ising-1)*igp+1, 2*ising*igp
                                 pp=xp(i)
!                                 fp=fpqb(i,lam)*pp**lab21*Qlp(i,iqa)
                                 fp=fpqb(i,lam)*xppl(i,lab21)*Qlp(i,iqa)
     $                                -Fqa * Q0p(iqa,i)
                                 res3 = res3 + wp(i)*fp
                              end do
!below gives zero if on-shell point is part of the k-grid.                            
!                              print*,'aq(ising)-qa,qa:',aq(ising)-qa,qa
                              res1c= qa*Log(aq(ising) - qa) + qa
     $                             *Log(aq(ising) + qa) + aq(ising)
     $                             *Log((aq(ising) + qa)/ (aq(ising) -
     $                             qa))

                              res3=res3+(res1b+res1c)*Fqa
                              
*     integrals coming after singularity
                              res4=0.d0
                              if (imax.lt.nqmi2*igp) stop 'imax4'
                              do i=2*ising*igp+1,nqmi2*igp
                                 pp=xp(i)
!                                 fp=fpqb(i,lam)*pp**lab21*Qlp(i,iqa)
                                 fp=fpqb(i,lam)*xppl(i,lab21)*Qlp(i,iqa)
                                 res4 = res4 + wp(i) * fp
                                 
c$$$C_TEST________________________________________________________________
c$$$!$omp critical 
c$$$            if(isnan(res4)) then
c$$$               open(unit=90, file = 'vmat', position='append')
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa
c$$$               write(90,*) 'i, res4', i,  res4
c$$$               write(90,*) 'fpqb, Qlp: ', fpqb(i,lam), Qlp(i,iqa) ! fpqb = NaN
c$$$               write(90,*) 'lam:', lam             
c$$$               close(90)
c$$$               stop 'stopped at posvmat -2-'
c$$$            end if
c$$$!$omp end critical 
c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^            
                              end do
                              res1=res1+res2+res3+res4
c$$$C_TEST________________________________________________________________
c$$$!$omp critical 
c$$$            if(isnan(res1).or.isnan(res2).or.isnan(res3).or.isnan(res4))
c$$$     1                             then
c$$$               open(unit=90, file = 'vmat', position='append')                                   
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              res1
c$$$               write(90,*) '(0) ----------------------------------'
c$$$               write(90,*) 'res2, res3, res4: ', res2, res3, res4
c$$$               write(90,*) 'lam:', lam 
c$$$               close(90)
c$$$               stop 'stopped at posvmat -1-'
c$$$            end if
c$$$!$omp end critical 
c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                                               
                           else
*     i.e. if singularity is at last k-mesh point then

                              nqmi1 = nqim-1
*     think of combining
                              res1 = 0.0d0
                              
*     integrals coming before singularity
                              if (imax.lt.2*nqmi1*igp) stop 'imax5'
                              do i=1,2*(nqmi1)*igp
                                 pp=xp(i)                           
!                                 fp=fpqb(i,lam)*pp**lab21*Qlp(i,iqa)   
                                 fp=fpqb(i,lam)*xppl(i,lab21)*Qlp(i,iqa)
                                 res1 = res1 + wp(i) * fp                                                                
                              end do
                           
*     integrals with singularities
                              Fqa = fpqb(nqmi2*igp+nqim,lam)*qa**lab21
                              res2=0.d0
                              if (imax.lt.(nqmi2-1)*igp) stop 'imax6'
                              do i=2*(nqmi1)*igp+1,(nqmi2-1)*igp
                                 pp=xp(i)
!                                 fp=fpqb(i,lam)*pp**lab21*Qlp(i,iqa)
                                 fp=fpqb(i,lam)*xppl(i,lab21)*Qlp(i,iqa)
     $                                -Fqa*Q0p(iqa,i)
                                 res2 = res2 + wp(i)*fp                                 
                              end do
                              
                              res1b=2*qa*lg2 - qa*Log(qa - aq(nqmi1)) +
     $                             2*qa*Log(qa) - qa*Log(aq(nqmi1) + qa)
     $                             - aq(nqmi1)*Log((aq(nqmi1) + qa)/
     $                             (qa-aq(nqmi1)))

                              res1b=res1b*Fqa
                              res2=res2+res1b

                              res3=0.d0
                              if (imax.lt.nqmi2*igp) stop 'imax7'
                              do i=(nqmi2-1)*igp+1,nqmi2*igp
                                 pp=xp(i)
                                 pp2 = pp * pp !xp2(i)
!                                 fp=fpqb(i,lam)*pp**lab21*Qlp(i,iqa)
                                 fp=fpqb(i,lam)*xppl(i,lab21)*Qlp(i,iqa)
     $                                -Fqa*(qa2+alfa2)*(qa2+alfa2) /(pp2
     $                                +alfa2)/(pp2+alfa2)*Q0p(iqa,i)
                                 res3 = res3 + wp(i)*fp
c$$$                              if (iqa.eq.1.and.iqb.eq.1) print*,
c$$$     >                              'i,fpqb,Qlp,Q0p:',i,fpqb(i,lam),
c$$$     >                              Qlp(i,iqa),Q0p(iqa,i)
                              end do
                              
                              res1c=(-2.*qa*ATan(qa/alfa) + 
     >                           Pi*qa)/(2.*alfa*(alfa2 + qa2))
                              
                              res1c=res1c*fpqb(2*nqim*igp+nqim,lam)*
     >                             qa**lab2*(qa2+alfa2)*(qa2+alfa2)
                              
                              res3=res3+res1c
c$$$                              if (iqa.eq.1.and.iqb.eq.1) print*,
c$$$     >                           'res1,res2,res3:',res1,res2,res3
                              res1=res1+res2+res3
c$$$C_TEST________________________________________________________________
c$$$!$omp critical 
c$$$            if(isnan(res1)) then
c$$$               open(unit=90, file = 'vmat', position='append')                                   
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              res1
c$$$               write(90,*) '(1) ---------------------------------------'
c$$$               write(90,*) 'res2,res3,Qlp(i,iqa): ',
c$$$     1              res2,res3,Qlp(iqa,(nqmi2-1)*igp+1) ! res2 = NaN, res3=NaN
c$$$               close(90)
c$$$               print*, 'ERROR: posvmat -3- '
c$$$            end if
c$$$!$omp end critical 
c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                                                                      
                                
                           endif                           
!                        else
                           
ccccccccccccccccccccccccccccccccccccccccccccccccccc
                           
c$$$*     calculation of p integral
c$$$
c$$$                           if(ising.ne.nqmi) then
c$$$*     i.e. if singularity is not at last k-mesh point
c$$$
c$$$                              res1 = 0.0d0
c$$$*     integrals coming before singularity
c$$$                              do i=1,2*(ising-1)*igp
c$$$                                 pp=xp(i)
c$$$                                 fp=fpqb(i,lam)*pp**(lab2+1)*Qlp(i,iqa)
c$$$                                 res1 = res1 + wp(i) * fp
c$$$                              end do
c$$$                              
c$$$*     integrals with singularities
c$$$                              res2=0.d0
c$$$                              do i=2*(ising-1)*igp+1,(2*ising-1)*igp
c$$$                                 pp=xp(i)
c$$$                                 fp=fpqb(i,lam)*pp**(lab2+1)*Qlp(i,iqa)
c$$$     >                              -fpqb(2*nqmi*igp+ising,lam)*
c$$$     >                              qa**(lab2+1)*Q0p(iqa,i)
c$$$                                 res2 = res2 + wp(i)*fp
c$$$                              end do
c$$$                              
c$$$                              res1b=2*qa*Log(2.) -
c$$$     >                           qa*Log(qa - aq(ising-1)) 
c$$$     >                           + 2*qa*Log(qa) -
c$$$     >                           qa*Log(aq(ising-1) + qa) 
c$$$     >                           - aq(ising-1)*Log((aq(ising-1) + qa)/
c$$$     >                           (-aq(ising-1) + qa))
c$$$
c$$$                              res1b=res1b*fpqb(2*nqmi*igp+ising,lam)*
c$$$     >                           qa**(lab2+1)
c$$$                              res2=res2+res1b
c$$$
c$$$                              res3=0.d0
c$$$                              do i=(2*ising-1)*igp+1,2*ising*igp
c$$$                                 pp=xp(i)
c$$$                                 fp=fpqb(i,lam)*pp**(lab2+1)*Qlp(i,iqa)
c$$$     >                              -fpqb(2*nqmi*igp+ising,lam)*
c$$$     >                              qa**(lab2+1)*Q0p(iqa,i)
c$$$                                 res3 = res3 + wp(i)*fp
c$$$                              end do
c$$$
c$$$                              res1c=-2*qa*Log(2.) +
c$$$     >                           qa*Log(aq(ising) - qa) -
c$$$     >                           2*qa*Log(qa) +  qa*Log(aq(ising) + qa)
c$$$     >                           + aq(ising)*Log((aq(ising) + qa)/
c$$$     >                           (aq(ising) - qa))
c$$$
c$$$                              res1c=res1c*fpqb(2*nqmi*igp+ising,lam)*
c$$$     >                           qa**(lab2+1)
c$$$                              res3=res3+res1c
c$$$                           
c$$$*     integrals coming after singularity
c$$$                              res4=0.d0
c$$$                              do i=2*ising*igp+1,2*nqmi*igp
c$$$                                 pp=xp(i)
c$$$                                 fp=fpqb(i,lam)*pp**(lab2+1)*Qlp(i,iqa)
c$$$                                 res4 = res4 + wp(i) * fp
c$$$                              end do
c$$$                              res1=res1+res2+res3+res4
c$$$                           else
c$$$*     i.e. if singularity is at last k-mesh point then
c$$$                              
c$$$*     think of combining
c$$$                              res1 = 0.0d0
c$$$
c$$$*     integrals coming before singularity
c$$$                              do i=1,2*(nqmi-1)*igp
c$$$                                 pp=xp(i)                              
c$$$                                 fp=fpqb(i,lam)*pp**(lab2+1)*Qlp(i,iqa)
c$$$                                 res1 = res1 + wp(i) * fp
c$$$                              end do
c$$$                           
c$$$*     integrals with singularities
c$$$                              res2=0.d0
c$$$                              do i=2*(nqmi-1)*igp+1,(2*nqmi-1)*igp
c$$$                                 pp=xp(i)
c$$$                                 fp=fpqb(i,lam)*pp**(lab2+1)*Qlp(i,iqa)
c$$$     >                              -fpqb(2*nqmi*igp+nqmi,lam)*
c$$$     >                              qa**(lab2+1)*Q0p(iqa,i)
c$$$                                 res2 = res2 + wp(i)*fp
c$$$                              end do
c$$$                              
c$$$                              res1b=2*qa*Log(2.) -
c$$$     >                           qa*Log(qa - aq(nqmi-1)) 
c$$$     >                           + 2*qa*Log(qa) -
c$$$     >                           qa*Log(aq(nqmi-1) + qa) 
c$$$     >                           - aq(nqmi-1)*Log((aq(nqmi-1) + qa)/
c$$$     >                           (-aq(nqmi-1) + qa))
c$$$
c$$$                              res1b=res1b*fpqb(2*nqmi*igp+nqmi,lam)*
c$$$     >                           qa**(lab2+1)
c$$$                              res2=res2+res1b
c$$$
c$$$                              res3=0.d0
c$$$                              do i=(2*nqmi-1)*igp+1,2*nqmi*igp
c$$$                                 pp=xp(i)
c$$$                                 fp=fpqb(i,lam)*pp**(lab2+1)*Qlp(i,iqa)
c$$$     >                              -fpqb(2*nqmi*igp+nqmi,lam)*
c$$$     >                              qa**(lab2+1)*(qa2+alfa2)*(qa2+alfa2)
c$$$     >                              /(pp*pp+alfa2)
c$$$     >                              /(pp*pp+alfa2)*Q0p(iqa,i)
c$$$                                 res3 = res3 + wp(i)*fp
c$$$                              end do
c$$$                              
c$$$                              res1c=(-2*qa*ATan(qa/alfa) + 
c$$$     >                           Pi*qa)/(2.*alfa*(alfa2 + qa2))
c$$$
c$$$                              res1c=res1c*fpqb(2*nqmi*igp+nqmi,lam)*
c$$$     >                           qa**lab2*(qa2+alfa2)*(qa2+alfa2)
c$$$                              res3=res3+res1c
c$$$
c$$$                              res1=res1+res2+res3
c$$$                           endif
c$$$                        else
                           
*     i.e. when on-shell point is negative

*     calculation of p integral

!                           res1 = 0.0d0
!                           do i=1,2*nqim*igp
!                              pp=xp(i)
!                              fp=fpqb(i,lam)*pp**lab21*Qlp(i,iqa)
!                              res1 = res1 + wp(i) * fp   
c$$$C_TEST________________________________________________________________
c$$$!$omp critical 
c$$$            if(isnan(res1)) then
c$$$               open(unit=90, file = 'vmat', position='append')                                   
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              res1
c$$$               write(90,*) '(2) ----------------------------------'
c$$$               write(90,*) ' i, lam, fp: ', i, lam, fp, fpqb(i,lam)
c$$$               write(90,*) ' fpqb(i,lam), Qlp(i,iqa): ',fpqb(i,lam),
c$$$     $              Qlp(i,iqa)
c$$$               close(90)
c$$$            end if
c$$$!$omp end critical 
c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                             
!                           end do                                    
!                        endif
                        
C-----Andrey --------------------------------------------------------
c$$$c     calculate contribution of dQl in the integral on the composite
c$$$c     mesh. 
c$$$                        if (alkali) then ! Gausian quadrature is used 
c$$$C     go to 10
c$$$                           open (unit=99, file='integrand1') ! TESTING
c$$$                           res5 = 0.0d0 
c$$$                           do i=1,2*nqmi*igp
c$$$                              pp=xp(i)
c$$$                              fp=fpqb(i,lam)*pp**(lab2+1)*dQlp(iqa,i)
c$$$                              
c$$$c     weights are in dQlp(iqa,i)
c$$$c                              res5 = res5 + wp(i) * fp 
c$$$                              res5 = res5 + fp                              
c$$$                              write(99,'(4(e12.6,1x))') pp, fp, res5
c$$$                           end do
c$$$                           close (99) ! TESTING
c$$$                           
c$$$C     integration on kgrid mesh
c$$$ 10                        continue                           
c$$$  end if

C     ANDREY: this block was used when contributin of dQlp was
C     calculated separately
                        
c$$$                        res5 = 0.0d0
c$$$                        go to 220 ! skip 
c$$$                        if (alkalitmp) then ! Fejer quadrature is used
c$$$c                           open (unit=99, file='integrand2') ! TESTING                           
c$$$                           do i=1,nqmi2*igp
c$$$                              pp=xp(i)
c$$$                              fp=fpqb(i,lam)*pp**lab21*dQlp2(iqa,i)
c$$$c                              res5 = res5 + wp(i) * fp
c$$$                              res5 = res5 + fp ! weights are in dQlp(iqa,i)
c$$$c$$$                          if ( dQlp2(iqa,i).gt.1.e-10) then
c$$$c$$$                             write(99,'(4(e12.6,1x))') pp, fp, res5
c$$$c$$$                          end if
c$$$                           end do
c$$$c                           close (99) ! TESTING                           
c$$$ 20                        continue                           
c$$$                        end if
c$$$ 220                    continue
c$$$C                        res5 = 0.0d0 ! TEST

                        
c$$$                        if(.false.) then ! rectangle rule                                                     
c$$$                           open (unit=99, file='integrand2') ! TESTING
c$$$                           res5 = 0.0d0
c$$$                           
c$$$*     - zero contribution for  [0,qq(1)] due to pp**(lab2+1)
c$$$                           
c$$$*     - [qq(1),qq(nqgen)] contribution:                           
c$$$                           do i = 1, nqgen-1
c$$$                              pp = aq(i)
c$$$                              dq = qq(i+1)-qq(i)
c$$$                              fp = gpqb(i,lam)*pp**(lab2+1)*dQlp2(iqa,i)
c$$$                              res5 = res5 + fp * dq
c$$$                              write(99,'(4(e12.6,1x))') pp,fp,res5
c$$$                           end do
c$$$                           
c$$$*     - [qq(nqgen),oo[ contribution:
c$$$                           do i = imax2,nqgen,-1
c$$$                              j2=(2*nqgen-1)*igp+i
c$$$                              pp = xp(j2)                              
c$$$                              fp = gpqb(i,lam)*pp**(lab2+1)*dQlp2(iqa,i)
c$$$                              res5 = res5 + wp(j2) * fp
c$$$                              write (99,'(4(e12.6,1x))') pp,fp,res5
c$$$                           end do
c$$$                           close (99) ! TESTING
c$$$                           continue
c$$$                        end if  ! alkali
                        
C     res1=res1+res5
        endif !Charge check
***** calculations of the p-integral using chebyshev mesh
c     open(33,file='fpold.dat')
c                           call transc(lam,lab2,res1test,qb,qa,Lla,lb,
c     >                        la,nb,na,nbi,pc1)
                           !print*,' chebyshev  =',res1test
c     close(33)
c     close(22)
c                           if(ising.eq.nqmi) pause'chto delat?'
                           
*     calculations using spline (10 to 20 times faster than trans)
c      call transpl(lam,lab2,res1t,x,ylam,y2lam,qa,Lla,lm)
      !print*,' trans0 and trans1: ',res0(lam),res1,res1t
        !if (iqb.eq.2) then
        !if (lam.gt.0) then
        !open(556, file='t-check')
        !write(556,*)qa,res1
        !endif
        !endif
        if (abs(zasym).gt.0.d0) then
           result=qa*res1/(2.d0*pi**2.d0)
           if ((eta.gt.6.0d0).or.(eta.lt.-1000.d0)) then !where subtraction is inaccurate
              result=0.d0
           endif
        else   
           result = qal2 * res0(lam) + res1 / pi
        endif
        !print*,qa,qb,result/qa
        reslam=reslam+(2*lam+1)*wigner*result
c$$$c$$$C_TEST________________________________________________________________
c$$$!$omp critical  
c$$$            if(isnan(reslam).or.isnan(result)) then
c$$$               open(unit=90, file = 'vmat', position='append')   
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              result, reslam
c$$$               write(90,*) 'reslam: lam, qal2',  lam, qal2 
c$$$                write(90,*) 'res0(lam), res1, pi: ', res0(lam), res1, pi ! res1 = NaN
c$$$               close (90)
c$$$            end if
c$$$!$omp end critical 
c$$$c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                           
      end do
                        
      suml2=suml2+(2*l2+1)*reslam*w3j2
                        
c$$$c$$$C_TEST________________________________________________________________
c$$$!$omp critical 
c$$$            if(isnan(suml2)) then
c$$$               open(unit=90, file = 'vmat', position='append')   
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              suml2
c$$$               write(90,*)'reslam: ', reslam
c$$$               close (90)
c$$$               stop
c$$$            end if
c$$$!$omp end critical 
c$$$c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                   
                      
                        
                     end do
                     suml1=suml1+(2*l1+1)*suml2*w3j1

c$$$c$$$C_TEST________________________________________________________________
c$$$!$omp critical 
c$$$            if(isnan(suml1)) then
c$$$               open(unit=90, file = 'vmat', position='append')                       
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              suml1
c$$$               write(90,*)'suml2 : ', suml2
c$$$               close (90)
c$$$               stop
c$$$            end if
c$$$!$omp end critical 
c$$$c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                   
                      
                  end do
                  sumla1=sumla1+suml1*cla

c$$$C_TEST________________________________________________________________
c$$$!$omp critical 
c$$$            if(isnan(sumla1)) then
c$$$               open(unit=90, file = 'vmat', position='append')               
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              sumla1
c$$$               write(90,*)'suml1, cla : ', suml1, cla
c$$$               close (90)               
c$$$            end if
c$$$!$omp end critical          
c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                   
 
                  
               end do
               sumlb1=sumlb1+sumla1/clb

c$$$C_TEST________________________________________________________________
c$$$!$omp critical                       
c$$$            if(isnan(sumlb1)) then               
c$$$               open(unit=90, file = 'vmat', position='append')
c$$$               write(90,'(6(2X, I4), 1(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              sumlb1
c$$$               write(90,*)'sumla1, clb : ', sumla1, clb
c$$$               close (90)
c$$$            end if
c$$$!$omp end critical          
c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                   
            end do

c            close (90)
c            stop
            result=sumlb1*c0*qb*sqrt(2.0d0)
c     print*,'posv',result

c$$$C_TEST________________________________________________________________
c$$$!$omp critical  
c$$$            if(isnan(result)) then
c$$$               open(unit=90, file = 'vmat', position='append')
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              result
c$$$               write(90,*)'sumlb1, c0, qb : ', sumlb1, c0, qb
c$$$               close (90)
c$$$            end if
c$$$!$omp end critical              
c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                   
                                      
* V-matrix is multiplied by qb*qa. Namely, qb*qa*V is used as driving term
* and in the kernal of integral equations. The solution is qb*qa*T.
* However, CCC prints  out  just V without qb*qa factor.
* (multiplification by qa is done earlier, see 'result')

* V-matrix is also multiplied by sqrt(mu=2)

c            write(99,'(''qb='',e23.6,'' qa='',e23.6,'' v='',e23.7)')
c     >      qb,qa,real(result)

c           write(6,*) 'at qb =',qb,' and qa =',qa
c           write(6,*) '  posVmat =',result/qb/qa/sqrt(2.0)
c      call dirvmat(dirres,qb,qa,J,Llb,Lla,lb,la,nb,na,nbi,En,bohr1)
c      print*,'  dirVmat   H->H:',dirres*2./pi
c      call dirvmat(dirres,qb,qa,J,Llb,Lla,lb,la,nb,na,nbi,En,bohr2)
c      print*,'  dirVmat Ps->Ps:',dirres*2./pi

*            vmat(kii,kff+1) = vmat(kii,kff+1) - real(result) ! old statement
            if (nchf.ge.nchi) then
c$$$               vmat(kff,kii) = vmat(kff,kii) + real(result)
               vmatt(iqb,iqa) = vmatt(iqb,iqa) + result
!               vmatt(iqb,iqa,1) = vmatt(iqb,iqa,1) + result
            else
c$$$               vmat(kii,kff) = vmat(kii,kff) + real(result)
               vmatt(iqa,iqb) = vmatt(iqa,iqb) + result
!               vmatt(iqa,iqb,1) = vmatt(iqa,iqb,1) + result
            endif

c$$$C_TEST________________________________________________________________
c$$$!$omp critical            
c$$$            if(isnan(result)) then
c$$$               open(unit=90, file = 'vmat', position='append') 
c$$$               write(90,'(6(2X, I4), 2(2X,E13.5))') lb,nb,la,na,iqb,iqa,
c$$$     $              vmatt(iqa,iqb,0), result
c$$$               !print*,vmatt(iqa,iqb,0), vmatt(iqa,iqb,1)
c$$$               close (90)
c$$$               !stop
c$$$            end if            
c$$$!$omp end critical   
c$$$C_TEST^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         end do
        !  print*, 'posvmat', iqa         
      end do
      return
      end subroutine posvmat

*------------------------------------------------------

      subroutine f0zpart(Nl,rlam,bohr,nn,ll,nni,pp2,res0,res1,pc0,pc1)

      include 'par.f'
      include 'par.pos'
      parameter (maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      dimension pc0(0:ncmax-1,ncmax,0:lamax),
     >     pc1(0:ncmax-1,ncmax,0:lamax)

      if(Nl.eq.0) then
         call geigen(bohr,nn,ll,pp2,res1)
         res0=res1*((bohr*nn)**2*pp2+1.d0)/2.d0/bohr/nn/nn
      else
c-- calculations without summation by mathematica
c   slow; for testing purposes only;
c         call gpseudo(ll,nni,Nl,rlam,pp2,res0,res1)
c--
         brap=4.d0*pp2/rlam/rlam+1.d0
         x=(brap-2.d0)/brap

         sk0=pc0(Nl-1,nni,ll)
         sk1=pc1(Nl-1,nni,ll)
         do k=Nl-2,0,-1
            sk0=sk0*x+pc0(k,nni,ll)
            sk1=sk1*x+pc1(k,nni,ll)
         end do

c-- or alternatively  (this compact form takes longer)
c         call sumk0an(x,ll,Nl,a,sk0)
c         call sumk1an(x,ll,Nl,a,sk1)
c--
         some=(8d0/(rlam*brap))**(ll+1)*factrl(ll)*0.5d0
         res0=some*sk0
         res1=some/rlam/brap*sk1*4.d0
      endif

c-- numerical calculations using Igor's rmesh;
c   extremely slow; for testing purposes only.
c      call gnlp(bohr,ll,pp2,nni,res0,res1) !,res2,res3)
c--
      return
      end

*---------------------------------------------

      subroutine f0zpartnum(Nl,rlam,bohr,nn,ll,nni,pp2,res0,res1) !,pc0,pc1)

      include 'par.f'
      include 'par.pos'
      parameter (maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
c$$$      dimension pc0(0:ncmax-1,ncmax,0:6),pc1(0:ncmax-1,ncmax,0:6)
      
      call gnlp(bohr,ll,pp2,nni,res0,res1) !,res2,res3)
      
      return
      end

*---------------------------------------------

      real*8 function f1z(pp,z,qb,lb,la,nb,na,nbi,pc1)
      include 'par.f'
      include 'par.pos'
      parameter (maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
c      common/pcfs/pc1(0:ncmax-1,ncmax,0:lnabmax)
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      dimension pc1(0:ncmax-1,ncmax,0:lamax)!(0:19,20,0:3)
        
      qb2=qb*qb
      pp2=pp*pp
      qbpp=qb*pp
      pb2=0.25d0*qb2 + pp2 - qbpp*z
      pa2=qb2 + pp2 - 2.d0*qbpp*z

      Nl=npsstates(1,la)
      if(Nl.eq.0) then
         bohr=1.d0
         call geigen(bohr,na,la,pa2,resH1)
      else
         rlam1=rlambda(1,la)
         brap1=4d0*pa2/rlam1/rlam1+1.d0
         x1=(brap1-2d0)/brap1

         sk1=pc1(Nl-1,na,la)
         do k=Nl-2,0,-1
            sk1=sk1*x1+pc1(k,na,la)
         end do
         resH1=(8d0/(rlam1*brap1))**(la+2)*factrl(la)*sk1*0.25d0
      endif

      Nl=npsstates(2,lb)
      if(Nl.eq.0) then
         bohr=2.d0
         call geigen(bohr,nb,lb,pb2,resP1)
      else
         rlam2=rlambda(2,lb)
         brap2=4d0*pb2/rlam2/rlam2+1d0
         x2=(brap2-2d0)/brap2

         sk2=pc1(Nl-1,nbi,lb)
         do k=Nl-2,0,-1
            sk2=sk2*x2+pc1(k,nbi,lb)
         end do
         resP1=(8d0/(rlam2*brap2))**(lb+2)*factrl(lb)*sk2*0.25d0
      endif

      f1z=resP1*resH1
      return
      end

*------------------------------------------------------


      subroutine zintegc(lam,x,y,y2,qb,lb,la,nb,na,nbi,pc1)

      include 'par.f'
      include 'par.pos'
      implicit real*8 (a-h,o-z)
      parameter (maxl=8*ltmax)
      common/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      dimension pc1(0:ncmax-1,ncmax,0:lamax)!(0:19,20,0:3)
      dimension x(ipm),y(ipm),y2(ipm),u(ipm)
      data yp1,ypn /0.d0,0.d0/

      Nl1=npsstates(1,la)
      Nl2=npsstates(2,lb)
      qb2=qb*qb
      x(1)=0.d0
      y(1)=0.d0
      do ips=2,ipm-1
         pp=dble(ips-1)/dble(ipm-ips)
         resz=0.d0
         pp2=pp*pp
         qbpp=qb*pp
         bb=0.25d0*qb2 + pp2
         aa=qb2 + pp2
         do iz=1,igz
            pb2=bb - qbpp*xgz(iz)
            pa2=aa - 2.d0*qbpp*xgz(iz)

            if(Nl1.eq.0) then
               bohr=1.d0
               call geigen(bohr,na,la,pa2,sk1)
            else
               rlam1=rlambda(1,la)
               brap1=4d0*pa2/rlam1/rlam1+1.d0
               x1=(brap1-2d0)/brap1
               sk1=pc1(Nl1-1,na,la)
               do k=Nl1-2,0,-1
                  sk1=sk1*x1+pc1(k,na,la)
               end do
               sk1=sk1*(8d0/rlam1/brap1)**(la+2)*factrl(la)*0.25d0
            endif

            if(Nl2.eq.0) then
               bohr=2.d0
               call geigen(bohr,nb,lb,pb2,sk2)
            else
               rlam2=rlambda(2,lb)
               brap2=4d0*pb2/rlam2/rlam2+1d0
               x2=(brap2-2d0)/brap2
               sk2=pc1(Nl2-1,nbi,lb)
               do k=Nl2-2,0,-1
                  sk2=sk2*x2+pc1(k,nbi,lb)
               end do
               sk2=sk2*(8d0/rlam2/brap2)**(lb+2)*factrl(lb)*0.25d0
            endif

* sk2*sk1 must be stored and re-used when lambda changes
* this to be done

            resz=resz + sk2*sk1*pl(lam,iz)
         end do

         x(ips)=pp
* function * argument**2 is gonna be splined,
* rather than function itself (for convinience)
         y(ips)=resz*pp2
      end do
      x(ipm)=1.0d30
      y(ipm)=0.0d0

* splining

      if (yp1.gt..99d30) then
         y2(1)=0.d0
         u(1)=0.d0
      else
         y2(1)=-0.5d0
         u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,ipm-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.d0)/p
         u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     >      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99d30) then
         qn=0.d0
         un=0.d0
      else
         qn=0.5d0
         un=(3.d0/(x(ipm)-x(ipm-1))) *
     >      (ypn-(y(ipm)-y(ipm-1))/(x(ipm)-x(ipm-1)))
      endif
      y2(ipm)=(un-qn*u(ipm-1))/(qn*y2(ipm-1)+1.d0)
      do k=ipm-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      end

*--------------------------------------------------------------------


      subroutine transpl(lam,lab2,result,xa,y,y2,qa,Lla,lm)

      include 'par.f'
      include 'par.pos'
      implicit real*8 (a-h,o-z)
      common/cheb/x(ich),w(ich)
      common/funQl/Qlarray(0:ltmax,ich)
      dimension xa(ipm),ylam(0:lm,ipm),y2lam(0:lm,ipm)
      dimension y(ipm),y2(ipm)

        !open(555,file="t-test")
        !write(555,*)"#qa",qa
      result = 0.0d0
      do i=1,ich
         t = x(i)
         ft=0.0d0
         QlQ0=Qlarray(Lla,i)
         p1=qa*t
         p2=qa/t
* spline interpolating

c         do ips=1,ipm
c            y(ips)=ylam(lam,ips)
c            y2(ips)=y2lam(lam,ips)
c            write(6,9) xa(ips),y(ips),ips,lam
c 9          format('po',2e16.7,2i4)
c         end do

         klo=1
         khi=ipm
 1       if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if(xa(k).gt.p1)then
               khi=k
            else
               klo=k
            endif
            goto 1
         endif
         h=xa(khi)-xa(klo)
         if (h.eq.0.d0) stop 'bad x input: check transpl'
         a=(xa(khi)-p1)/h
         b=(p1-xa(klo))/h
         fp1=a*y(klo)+b*y(khi)+
     >      ((a**3-a)*y2(klo)+(b**3-b)*y2(khi))*(h**2)/6.d0

         klo=1
         khi=ipm
 11      if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if(xa(k).gt.p2)then
               khi=k
            else
               klo=k
            endif
            goto 11
         endif
         h=xa(khi)-xa(klo)
         if (h.eq.0.d0) pause 'bad x input.'
         a=(xa(khi)-p2)/h
         b=(p2-xa(klo))/h
         fp2=a*y(klo)+b*y(khi)+
     >      ((a**3-a)*y2(klo)+(b**3-b)*y2(khi))*(h**2)/6.d0
         ft = (fp1*p1**lab2 + fp2*p2**lab2) * QlQ0 / t
        
        !write(555,*)sngl(t),sngl(p1),sngl(p2),sngl(fp1),sngl(fp2)
* function ft is regular (=0) when t is going to 0
c         write(10,*) t,ft
         result = result + w(i) * ft
      end do
c      pause ' wrote to unit 10'
      return
      end
                           
*---------------------------------------------

      subroutine transc(lam,lab2,result,qb,qa,
     >   Lla,lb,la,nb,na,nbi,pc1)

      include 'par.f'
      include 'par.pos'
      parameter (maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/cheb/x(ich),w(ich)
      common/funQl/Qlarray(0:ltmax,ich)
      common/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      common/gausp/xgp(igpm),wgp(igpm),igp
c$$$      common/const/pi
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      dimension pc1(0:ncmax-1,ncmax,0:lamax)!(0:19,20,0:3)
      external f1z
        
***** calculations of the p-integral using chebyshev mesh without spline
      result = 0.0d0
      pi = acos(-1d0)
      do i=1,ich
         t = x(i)
         QlQ0=Qlarray(Lla,i)
         p1=qa*t
         p2=qa/t
         resz1=0.d0
         resz2=0.d0
         do iz=1,igz
            z=xgz(iz)
            resz1=resz1+f1z(p1,z,qb,lb,la,nb,na,nbi,pc1)*pl(lam,iz)
            resz2=resz2+f1z(p2,z,qb,lb,la,nb,na,nbi,pc1)*pl(lam,iz)
        !print*,iz,f1z(p1,z,qb,lb,la,nb,na,nbi,pc1),pl(lam,iz)
         end do
* to be done: f1z doesn't depend on lambda - save!
         fp1=resz1*p1**(lab2+1)
         fp2=resz2*p2**(lab2+1)

c         write(33,*) p1,fp1!*log((p1+qa)/(qa-p1))
c         write(33,*) p2,fp2!*log((p2+qa)/(p2-qa))

         fts = fp1*p1 + fp2*p2
         result = result + w(i) * fts * QlQ0 / t
      end do

c      print*,' at  qa qb:',qa,qb
c      print*,'chebyshev result=',result,'is used'
      return
      end
*---------------------------------------------

      subroutine transcnew(lam,lab2,result,qb,qa,
     >   Lla,lb,la,nb,na,nbi,pc1)

      include 'par.f'
      include 'par.pos'
      parameter (maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/cheb/x(ich),w(ich)
      common/funQl/Qlarray(0:ltmax,ich)
      common/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
      common/gausp/xgp(igpm),wgp(igpm),igp
c$$$      common/const/pi
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)
      dimension pc1(0:ncmax-1,ncmax,0:lamax)!(0:19,20,0:3)
      external f1z
       !Charlie 
***** calculations of the p-integral using chebyshev mesh without spline
      result = 0.0d0
      pi = acos(-1d0)
      qb2=qb*qb

      do i=1,ich
         t = x(i)
         QlQ0=Qlarray(Lla,i)
         p1=qa*t
         p2=qa/t
         pp1=p1*p1
         pp2=p2*p2
         qbp1=qb*p1
         qbp2=qb*p2

         resz1=0.d0
         resz2=0.d0
         do iz=1,igz

            z=xgz(iz)

            
              pb12=0.25d0*qb2 + pp1 - qbp1*z
              pb22=0.25d0*qb2 + pp2 - qbp2*z
              pa12=qb2 + pp1 - 2.d0*qbp1*z
              pa22=qb2 + pp2 - 2.d0*qbp2*z

              Nl=npsstates(1,la)
              if(Nl.eq.0) then
                 bohr=1.d0
                 call geigen(bohr,na,la,pa12,resH11)
                call geigen(bohr,na,la,pa22,resH12)
              else
                 rlam1=rlambda(1,la)
                 brap11=4d0*pa12/rlam1/rlam1+1.d0
                 x11=(brap11-2d0)/brap11

                brap12=4d0*pa22/rlam1/rlam1+1.d0
                 x12=(brap12-2d0)/brap12

                 sk11=pc1(Nl-1,na,la)
                 sk12=sk11
                 do k=Nl-2,0,-1
                    sk11=sk11*x11+pc1(k,na,la)
                    sk12=sk12*x12+pc1(k,na,la)
                 end do
              resH11=(8d0/(rlam1*brap11))**(la+2)*factrl(la)*sk11*0.25d0
              resH12=(8d0/(rlam1*brap12))**(la+2)*factrl(la)*sk12*0.25d0
              endif

              Nl=npsstates(2,lb)
              if(Nl.eq.0) then
                 bohr=2.d0
                 call geigen(bohr,nb,lb,pb12,resP11)
                call geigen(bohr,nb,lb,pb22,resP12)
              else
                 rlam2=rlambda(2,lb)
                 brap21=4d0*pb12/rlam2/rlam2+1d0
                 x21=(brap21-2d0)/brap21

                brap22=4d0*pb22/rlam2/rlam2+1d0
                 x22=(brap22-2d0)/brap22

                 sk21=pc1(Nl-1,nbi,lb)
                 sk22=sk21

                 do k=Nl-2,0,-1
                    sk21=sk21*x21+pc1(k,nbi,lb)
                   sk22=sk22*x22+pc1(k,nbi,lb)
                 end do
              resP11=(8d0/(rlam2*brap21))**(lb+2)*factrl(lb)*sk21*0.25d0
              resP12=(8d0/(rlam2*brap22))**(lb+2)*factrl(lb)*sk22*0.25d0
              endif

              f1z1=resP11*resH11
              f1z2=resP12*resH12

            resz1=resz1+f1z1*pl(lam,iz)
            resz2=resz2+f1z2*pl(lam,iz)
        !print*,iz,f1z(p1,z,qb,lb,la,nb,na,nbi,pc1),pl(lam,iz)
         end do
* to be done: f1z doesn't depend on lambda - save!
         fp1=resz1*p1**(lab2+1)
         fp2=resz2*p2**(lab2+1)

c         write(33,*) p1,fp1!*log((p1+qa)/(qa-p1))
c         write(33,*) p2,fp2!*log((p2+qa)/(p2-qa))

         fts = fp1*p1 + fp2*p2
         result = result + w(i) * fts * QlQ0 / t
      end do

c      print*,' at  qa qb:',qa,qb
c      print*,'chebyshev result=',result,'is used'
      return
      end


*---------------------------------------------------

      subroutine geigen(bohr,nn,ll,pp2,result)

      include 'par.f'
      parameter(maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)

      bn=bohr*nn
      basis=bn*bn*pp2
      brap=basis+1.d0
      bram=basis-1.d0

      select case(nn)
      case(1)
         res0=2.d0*bohr*bohr/brap
         result=res0*2.d0*sqrt(dble(nn)/bohr)/brap
      case(2)
         if(ll.eq.0) then
            res0=8.d0*bohr*bohr*bram/brap/brap
         else
            res0=32.d0*bohr**3/sqrt(3.d0)/brap/brap
         endif
         result=res0*2.d0*sqrt(dble(nn)/bohr)/brap
      case(3)
         if(ll.eq.0) then
            bp2=pp2*bohr*bohr
            res0=18.d0*bohr*bohr*(81.d0*bp2*bp2-30.d0*bp2+1.d0)/brap**3
         else
            if(ll.eq.1) then
               res0=144.d0*bohr**3/sqrt(2.d0)*bram/brap**3
            else
               res0=864.d0*bohr**4/sqrt(10.d0)/brap**3
            endif
         endif
         result=res0*2.d0*sqrt(dble(nn)/bohr)/brap
      case default
         nr=nn-ll-1
         Anl=2.d0**(2*ll+2)*bn**(ll+1)*sqrt(bn)*factrl(ll)
         Anl=Anl*sqrt(dble(nn))*sqrfct(nr)/sqrfct(nn+ll)
         result=Anl/brap**(ll+2)
         if(nr.eq.0) return

*-- gegenbauer polynomial
         sum=1.d0
         bk=1.d0
         do k=1,nr
            bk=-bk*dble(nr+2*ll+1+k)*dble(nr+1-k)/k/(dble(ll+k)+0.5d0)
            sum=sum+bk/brap**k
         enddo
         gegen=sum*factrl(nr+2*ll+1)/factrl(nr)/factrl(2*ll+1)
*--
         result=result*gegen
      end select
      return
      end

*------------------------------------------------------

      subroutine gegenbauer(m,ll,brap,gegen)

      include 'par.f'
      parameter(maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)

c      gegen=1.
c      if(m.eq.0) return
      sum=1.d0
      bk=1.d0
      do k=1,m
         bk=-bk*dble(m+2*ll+1+k)*dble(m+1-k)/k/(dble(ll+k)+0.5d0)
         sum=sum+bk/brap**k
      enddo
      gegen=sum*factrl(m+2*ll+1)/factrl(m)/factrl(2*ll+1)
      return
      end

*------------------------------------------------------

      subroutine gpseudo(ll,nopt,Nl,rlam,pp2,res0,res1)

      include 'par.f'
      parameter(maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)
      dimension gegenm(15)

      brap=4.d0*pp2/rlam/rlam+1.d0
* Nl must be > 1
* k=Nl cycle
      summ0=1.d0
      summ1=dble(ll+1)
      do m=1,Nl-1

c         call gegenbauer(m,ll,brap,gegen)
c--   inlined
         sum=1.d0
         bk=1.d0
         do k=1,m
            bk=-bk*dble(m+2*ll+1+k)*dble(m+1-k)/k/(dble(ll+k)+0.5d0)
            sum=sum+bk/brap**k
         enddo
         gegen=sum*factrl(m+2*ll+1)/factrl(m)/factrl(2*ll+1)
c---

***** test for l=0 case
c
c         phi=acos((brap-2.d0)/brap)
c         gegen=sin((m+1)*phi)/sin(phi)
c
*****
         gegenm(m)=gegen
         summ0=summ0+gegen
         summ1=summ1+(m+ll+1)*gegen
      enddo
      ccc=cknd(Nl,nopt,ll)
      sumk0=summ0*ccc
      sumk1=summ1*ccc
      do k=1,Nl-1
         summ0=1.d0
         summ1=dble(ll+1)
         do m=1,k-1
            summ0=summ0+gegenm(m)
            summ1=summ1+(m+ll+1)*gegenm(m)
         enddo
         ccc=cknd(k,nopt,ll)
         sumk0=sumk0+summ0*ccc
         sumk1=sumk1+summ1*ccc
      enddo
      res0=2.d0**(3*ll+2)*factrl(ll)/rlam**(ll+1)/brap**(ll+1)*sumk0
      res1=2.d0**(3*ll+4)*factrl(ll)/rlam**(ll+2)/brap**(ll+2)*sumk1
      return
      end

*------------------------------------------------------

      subroutine gnlp(bohr,ll,pp2,nopt,res0,res1) 
C      
C     both for (hydrogen) atom and positronium (bohr)
C     both for R(r) and R(r)/r  (ipower) : Hydrogen target
C      
C     returns 
C     * res0 - Fourier transform of V(r) * R(r)
C     * res1 - Fourier transform of R(r) 
C
      use apar, only : alkali
      include 'par.f'
      implicit real*8 (a-h,o-z)
                         
      common/meshrr/ meshr,rmesh(maxr,3)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      real   enpsinb,psinb,rmesh                 
      pp=sqrt(pp2)
      res0=0.0d0
      res1=0.d0
      
      if (alkali.and.(bohr.eq.1)) then         
         do i = 1, istoppsinb(nopt, ll)
            rr = dble(rmesh(i,1))
            arg = rr * pp                                               
            call sbessel(arg, ll, besj)          
            term = dble(psinb(i,nopt,ll)) * besj * dble(rmesh(i,3))
            res0 = res0 + term * (1.d0-vaint(rr)) ! V * W.F.
            res1 = res1 + term * rr               ! W.F.
         enddo
      else                  
         do i = 1, istoppsinb(nopt, ll)
            rr = dble(rmesh(i,1))
            arg = rr * pp                        
            call sbessel(arg, ll, besj)           
            term = dble(psinb(i,nopt,ll)) * 
     1           besj * dble(rmesh(i,3))
            res0 = res0 + term  
            res1 = res1 + term * rr                        
c         call wpseudo(bohr,ll,rr,nopt,psir)
cc         print*,' wpseudo:  psir:',psir,rr
cc         print*,' Igor:   psinb:',psinb(i,nopt,ll)
cc         pause
c         res2 = res2 + psir * chi0 * dble(rmesh(i,3))
c         res3 = res3 + psir * rr * chi0 * dble(rmesh(i,3))            
         enddo                 
      end if
      
      ppll = pp ** ll
      res0 = res0 / ppll
      res1 = res1 / ppll
      return
      end subroutine gnlp
      

*------------------------------------------------------

      subroutine wpseudo(bohr,ll,rr,nopt,psir)

      include 'par.f'
      parameter(maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)

      if(bohr.eq.1.d0) then
         rlam=rlambda(1,ll)
         Nl=npsstates(1,ll)
      else
         rlam=rlambda(2,ll)
         Nl=npsstates(2,ll)
      endif
      arg=rlam*rr

c      rlam=rlam*2.
c      Nl=1

      psir=0.d0
      do k=1,Nl
         call laguerre(k-1,2*ll+2,arg,plag)
c         call laguer1(k-1,2*ll+2,arg,plag1)
c         if(cknd(k,nopt,ll).eq.0.) pause ' Cknd = 0'
         chlen=cknd(k,nopt,ll)*arg**(ll+1)*exp(-arg/2.d0)*plag
         psir=psir+chlen
      enddo
c      pause
      return
      end

*---------------------------------------------------

      subroutine laguerre(nn,ll,arg,result)

      include 'par.f'
      parameter(maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)

      sum=0.d0
      do mm=0,nn
         sum=sum+(-arg)**mm/factrl(nn-mm)/factrl(ll+mm)/factrl(mm)
      enddo
      result=sum*factrl(nn+ll)
c      if(result.eq.0.d0) pause 'Laguerre =0'
c      print*,' laguerre: res=',result,arg,sum
      return
      end

*---------------------------------------------------

      subroutine laguer1(nn,ll,arg,sum)

c      include 'par.f'
c      parameter(maxl=8*ltmax)
c      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
c     >   Qlfactor(0:maxl)
      implicit real*8 (a-h,o-z)

      sum=0.d0
      do mm=0,nn
         call laguerre(mm,ll-1,arg,result)
         sum=sum+result
      enddo
      print*,' laguer1: sum=',sum,arg
      print*,' '
      return
      end

*-----------------------------------------------------

      subroutine sbessel(X, L, BES)
      implicit real*8 (a-h,o-z)            
C     modifyed spherical Bessel function
      integer, intent(in) :: L
      real*8, intent(in) :: X
      real*8, intent(out) :: BES
      real*8, dimension(0:L) :: Y      
      call SF32D(X,L,Y,IERR)
      BES=Y(L)
      return
C     
C     inacurate for small X and L > 10
      if (L .lt. 0) then                       
         BES = 0d0                           
         return                                
      end if      
      if (X .eq. 0d0) then
         BES = 0d0 
         if (L .eq. 0)  BES = 1d0          
         return
      end if      
      select case(L)
      case(0)
         if (x.lt.1.e-2) then            
            x2 = x*x
            BES = 1.0 - (x2*(1.0 - ((1.0 - x2/42.)*x2)/20.))/6.
         else
            BES = sin(x)/x
         end if
      case(1)
         if (x.lt.1.e-2) then            
            x2 = x*x
            BES = (x*(1.-(x2*(1.-(x2*(1.-x2/54.))/28.))/10.))/3.
         else         
            BES = (sin(x)/x - cos(x))/x
         end if
      case(2)
         if (x.lt.1.e-2) then            
            x2 = x*x
            BES = (x2*(1.0 - ((1.0 - x2/36.)*x2)/14.))/15.
         else
            x2 = 3d0/x/x
            BES = sin(x)/x*(x2-1d0) - cos(x)*x2
         end if
      case(3)
         x2 = x*x
         if (x.lt.1.e-2) then
            BES =(x**3*(1.-((1.-x2/44.)*x2)/18.))/105. 
         else            
            BES = sin(x)/x2*(15d0/x2-6d0) - cos(x)/x*(15d0/x2-1d0)
         end if
      case default   ! this is from bes.f to save one call

         if (x.lt.2.0e-2.and.L.eq.4) then
            BES = x**4*(1.0 - x**2/22.0)/945.0
            return
         end if
         if (x.lt.0.1.and.L.eq.5) then
            BES = x**5*(1.0 - x**2/26.0)/10395.0
            return
         end if
         if (x.lt.0.2.and.L.eq.6) then
            BES = x**6*(1.0 - x**2/30.0)/135135.0
            return
         end if
         if (x.lt.0.4.and.L.eq.7) then
            BES = x**7*(1.0 - x**2/34.0)/2027025.0
            return
         end if
         if (x.lt.0.6.and.L.eq.8) then
            BES = x**8*(1.0 - x**2/38.0)/34459425.0
            return
         end if
         if (x.lt.0.8.and.L.eq.9) then
            BES = x**9*(1.0 - x**2/42.0)/654729075.0
            return
         end if
         if (x.lt.1.1.and.L.eq.10) then
            BES = x**10*(1.0 - x**2/46.0)/13749310575.0
            return
         end if
         
         
         arg =X - dacos(-1d0) * dble(L) /2d0 
         m1 = int(L/2d0)                         
         m2 = int(abs((L-1)/2d0))                
         if ((X .gt. 1d0) .and. (X .gt. dble(L))) then 
            c  = 1d0
            pp = 1d0
            X2 = 2d0 * X
            d  = dble(L*(L+1)) / X2
            qq = d
            if (m1 .ne. 0) then
               do i = 1, m1
                  i2 = 2*i
                  pp = - pp * dble(L-i2+1) * dble(L-i2+2)/
     :                 (dble(i2) * dble(i2-1)) * dble(L+i2) * 
     :                  dble(L+i2-1) / X2**2   
                  c  = c + pp
               end do
            end if
            if (m2 .ne. 0) then
               do i = 1, m2
                  i2 = 2*i
                  qq = - qq * dble(L+i2) * dble(L+i2+1) /
     :                 (dble(i2) * dble(i2+1)) * dble(L-i2) *
     :                 dble(L-i2+1) / X2**2
                  d = d + qq
               end do
            end if   
            BES = (c*sin(arg)  + d*cos(arg))/X
         else                 ! small |X|
            ii = 0
            BES = 0
            a1 = 1/dble(2*L+1)
            if (l .ne. 0) then
               do i = 1, L
                  a1 = a1 * 2d0 / dble(i+L) * X
               end do
            end if   
            if (abs(a1) .lt. 10e-11 ) goto 20
 10         BES = BES + a1
            ii = ii + 1
            Lii = L+ii
            a1 = a1 * (-1d0) * X**2 * dble(Lii) / dble(ii) /
     :           dble(2*Lii+1) / dble(2*Lii)         
            if (abs(a1) .gt. 10e-11) goto 10
 20         continue
         end if         
      end select  
      return
      end subroutine sbessel
      
      subroutine sbessel_qp(X, L, BES)
      integer, parameter :: qp = selected_real_kind(15, 307)
      !integer, parameter :: qp = selected_real_kind(33, 4931)
      real (kind = qp) ::  X, BES,arg,pp,c,X2,d,qq,a1
C     modifyed spherical Bessel function 
      if (L .lt. 0) then                       
         BES = 0d0                           
         return                                
      end if
      
      if (X .eq. 0d0) then
         BES = 0d0 
         if (L .eq. 0)  BES = 1d0          
         return
      end if
      
      select case(L)
      case(0)
         if (x.lt.1.e-2) then            
            x2 = x*x
            BES = 1.0 - (x2*(1.0 - ((1.0 - x2/42.)*x2)/20.))/6.
         else
            BES = sin(x)/x
         end if
      case(1)
         if (x.lt.1.e-2) then            
            x2 = x*x
            BES = (x*(1.-(x2*(1.-(x2*(1.-x2/54.))/28.))/10.))/3.
         else         
            BES = (sin(x)/x - cos(x))/x
         end if
      case(2)
         if (x.lt.1.e-2) then            
            x2 = x*x
            BES = (x2*(1.0 - ((1.0 - x2/36.)*x2)/14.))/15.
         else
            x2 = 3d0/x/x
            BES = sin(x)/x*(x2-1d0) - cos(x)*x2
         end if
      case(3)
         x2 = x*x
         if (x.lt.1.e-2) then
            BES =(x**3*(1.-((1.-x2/44.)*x2)/18.))/105. 
         else            
            BES = sin(x)/x2*(15d0/x2-6d0) - cos(x)/x*(15d0/x2-1d0)
         end if
      case default   ! this is from bes.f to save one call

         if (x.lt.2.0e-2.and.L.eq.4) then
            BES = x**4*(1.0 - x**2/22.0)/945.0
            return
         end if
         if (x.lt.0.1.and.L.eq.5) then
            BES = x**5*(1.0 - x**2/26.0)/10395.0
            return
         end if
         if (x.lt.0.2.and.L.eq.6) then
            BES = x**6*(1.0 - x**2/30.0)/135135.0
            return
         end if
         if (x.lt.0.4.and.L.eq.7) then
            BES = x**7*(1.0 - x**2/34.0)/2027025.0
            return
         end if
         if (x.lt.0.6.and.L.eq.8) then
            BES = x**8*(1.0 - x**2/38.0)/34459425.0
            return
         end if
         if (x.lt.0.8.and.L.eq.9) then
            BES = x**9*(1.0 - x**2/42.0)/654729075.0
            return
         end if
         if (x.lt.1.1.and.L.eq.10) then
            BES = x**10*(1.0 - x**2/46.0)/13749310575.0
            return
         end if
         
         
         arg =X - dacos(-1d0) * dble(L) /2d0 
         m1 = int(L/2d0)                         
         m2 = int(abs((L-1)/2d0))                
         if ((X .gt. 1d0) .and. (X .gt. dble(L))) then 
            c  = 1d0
            pp = 1d0
            X2 = 2d0 * X
            d  = dble(L*(L+1)) / X2
            qq = d
            if (m1 .ne. 0) then
               do i = 1, m1
                  i2 = 2*i
                  pp = - pp * dble(L-i2+1) * dble(L-i2+2)/
     :                 (dble(i2) * dble(i2-1)) * dble(L+i2) * 
     :                  dble(L+i2-1) / X2**2   
                  c  = c + pp
               end do
            end if
            if (m2 .ne. 0) then
               do i = 1, m2
                  i2 = 2*i
                  qq = - qq * dble(L+i2) * dble(L+i2+1) /
     :                 (dble(i2) * dble(i2+1)) * dble(L-i2) *
     :                 dble(L-i2+1) / X2**2
                  d = d + qq
               end do
            end if   
            BES = (c*sin(arg)  + d*cos(arg))/X
         else                 ! small |X|
            ii = 0
            BES = 0
            a1 = 1/dble(2*L+1)
            if (l .ne. 0) then
               do i = 1, L
                  a1 = a1 * 2d0 / dble(i+L) * X
               end do
            end if   
            if (abs(a1) .lt. 10e-11 ) goto 20
 10         BES = BES + a1
            ii = ii + 1
            Lii = L+ii
            a1 = a1 * (-1d0) * X**2 * dble(Lii) / dble(ii) /
     :           dble(2*Lii+1) / dble(2*Lii)         
            if (abs(a1) .gt. 10e-11) goto 10
 20         continue
         end if         
      end select  
      return
      end subroutine sbessel_qp
                  
*-------------------------------------------------

      subroutine Qltable(lstopm,igz,igp)

      include 'par.f'
      include 'par.pos'
      parameter( maxl=8*ltmax)
      implicit real*8 (a-h,o-z)    
c$$$      common/const/pi
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      common/fik/fik(0:19,0:19)
      common/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igza
      common/gauszR/xgzR(igzm),wgzR(igzm),plR(0:maxl,igzm),igzR
       common/gausp/xgp(igpm),wgp(igpm),igpa
      common/funQl/Qlarray(0:ltmax,ich)
      common/cheb/x(ich),w(ich)
      real*8 dfactrl(0:maxl),dfactrl2(0:maxl),arg(ich)
      
c      open(99,file='matritsa')
      pi = acos(-1d0)
      if(igz.gt.igzm .or. igp.gt.igpm) then
         print*, igz,igzm,igp,igpm
         stop'enlarge igzm and/or igpm in par.pos'
      endif
      igza=igz
      igpa=igp
      igzR=igza

      do i=1,ich
         t=x(i)
         arg(i)=(t*t+1.d0)/t/2.d0
      end do

      factrl(0)=1.d0
      hat(0)=1.d0
      sqrfct(0)=1.d0
      do i=1,max(40,lstopm) !maxl
         factrl(i)=factrl(i-1)*dble(i)
         hat(i)=sqrt(dble(2*i+1))
         sqrfct(i)=sqrt(factrl(i))
      end do

      do k=0,19
         do i=0,19
            fik(i,k)=factrl(k+i+1)/(factrl(i+2)*factrl(i))/
     >         (factrl(k+2)*factrl(k))
         enddo
      enddo
      
      print*
      print*,' posVmat: using ',igz,' points for z and ',
     >   igp,'*nqmi*2 for p integral'
      print*,' Forming tables of Ql and Pl'
      dfactrl(0)=1.d0
      dfactrl2(0)=1.d0
      Qlfactor(0)=1.d0
      do i=1,ich
         Qlarray(0,i)=1.d0
      end do
      do il=1,lstopm
         dfactrl(il)=dfactrl(il-1) * dble(il)
         dfactrl2(il)=dfactrl2(il-1) * dble(2*il+1)
         Qlfactor(il)=dfactrl(il)/dfactrl2(il)
*         do i=1,ich
*            Q0= 0.5d0*log((arg(i)+1.d0)/(arg(i)-1.d0))
*            call funleg(arg(i),il,resQl)
*            Qlarray(il,i)=resQl/Q0
*        end do
      enddo
      call polleg(2*lstopm+1,igz,igp)
      
      print*
      return
      end subroutine Qltable

*----------------------------------------------

      subroutine funleg(z, Lla, Q0, result)      
C     returns the Legendre polynomial of the second kind - Ql
C     for given Lla and z with the use of the recurrence
C     relation
            
      include 'par.f'
      include 'par.pos'
      parameter( maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
      real *8, dimension (0:Lla) :: QN

      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >     Qlfactor(0:maxl)      

*     arg of Legendre function is always .ge. 1
*     min[z]=1 at pp=qa
*     Mathematica cannot handle this region.
            
      Q0 = 0.5d0 * log((z+1d0)/(z-1d0+1d-10))
c$$$  Q0sp = real(Q0)      
      
      select case (Lla)
      case (0:4)
         zc = 4.0d0
      case (5)
         zc = 3.0d0
      case (6:7)
         zc = 2.5d0
      case (8:20)
         zc = 2.0d0
      case (21:)
         zc = 1d20
      end select
      
C     use analitical expression for small z
      if (z.lt.zc) then      
         select case(Lla)
         case(0)
            result = Q0
         case(1)
            result = z * Q0 - 1d0
         case(2)
            z3 = 3d0 * z
            result = ((z3 * z - 1d0) * Q0- z3)/2d0
         case(3)
            z2 = z*z
            result = ((5d0*z2 - 3d0)*z*Q0 - (5d0 * z2-4d0/3d0))/2d0
         case(4)
            z2 = z*z
            result = (((35d0*z2 - 30d0)*z2 + 3d0)*Q0 
     >           - (35d0*z2 - 55d0/3d0)*z)/8d0
         case(5)
            z2 = z * z
            z4 = z2 * z2
            result = -0.5333333333333333 + 6.125*z2 - 7.875*z4 + 0.9375
     $           *z*(1. - 4.666666666666666 * z2 + 4.2 * z4) * 2d0 * Q0
         case(6)
            z2 = z * z
            z4 = z2 * z2
            result= (-2.8875 + 14.875*z2 - 14.4375*z4)*z - 0.15625
     $           *(1d0-21d0*z2 + 62.99999999999999*z4-46.2*z2*z4)*2d0
     $           *Q0
         case(7)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            result=0.4571428571428571-10.6125*z2+34.375*z4 - 26
     $           .81249999999999*z6 - 1.09375*z*(1.-9.*z2+19.8*z4-12
     $           .25714285714285*z6)*2.d0*Q0
         case(8)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            result=(3.383705357142857 - 32.9140625*z2 + 77.0859375*z4 -
     $           50.2734375*z6)*z + 0.13671875* (1. - 36.*z2 + 198.*z4 -
     $           343.1999999999999*z6 + 183.8571428571428*z8)*2d0*Q0
         case(9)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            result=-0.4063492063492063 + 15.82477678571428*z2 - 92
     $           .7265625*z4 + 169.4401041666666*z6 - 94.9609375*z8 + 1
     $           .23046875*z* (1. - 14.66666666666666*z2 + 57.2*z4 - 81
     $           .7142857142857*z6 + 38.58730158730158*z8)*2d0*Q0
         case(10)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            result=(-3.817398313492063 + 59.68973214285714*z2 - 245
     $           .5578125*z4 + 367.1822916666666*z6 - 180.4257812499999
     $           *z8)*z - 0.123046875*(1. - 55.*z2 + 476.6666666666667
     $           *z4 - 1430.*z6 + 1736.428571428571*z8 - 733.15873015873
     $           *z8*z2)*2.d0*Q0

         case (11)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            z10 = z8 * z2
            r1 = 0.3694083694083694 - 344.44921875*z10 - 21
     $           .67392113095238*z2 +198.25*z4 - 622.8286458333333*z6 +
     $           787.3125*z8
            r2 =-2.70703125 + 344.44921875*z10 + 58.65234375*z2 - 351
     $           .9140625*z4 + 854.6484375*z6 - 902.12890625*z8
            result = r1 + r2 * z * Q0
         case (12)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            z10 = z8 * z2
            z12 = z6 * z6
            r1 = 4.207314495400433 - 660.1943359375*z10 - 96
     $           .25726996527777*z2 + 605.073828125*z4 - 1530.338671875
     $           *z6 + 1674.4059244791667*z8
            r2 = 0.2255859375 + 660.1943359375*z12 - 1894.470703125*z10
     $           -17.595703125*z2 + 219.9462890625*z4 - 997.08984375*z6
     $           +2029.7900390625*z8
            result = r1 * z + r2 * Q0
         case (13)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            z10 = z8 * z2
            z12 = z6 * z6
            r1 = -0.340992340992341 + 3537.9645182291665*z10 - 1269
     $           .6044921875*z12 + 28.09768584280303*z2 - 368
     $           .1101345486111*z4 + 1738.522265625*z6 - 3669.708984375
     $           *z8  
            r2 = 2.9326171875 - 3961.166015625*z10 + 1269.6044921875*z12
     $           -87.978515625*z2 + 747.8173828125*z4 - 2706.38671875*z6
     $           +4736.1767578125*z8
            result = r1 + r2 * z * Q0
         case (14)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            z10 = z8 * z2
            z12 = z6 * z6
            z14 = z10 * z4
            r1 = -4.564420117642774 + 7436.2548828125*z10 - 2448
     $           .52294921875*z12 + 143.5701448074495*z2 - 1271
     $           .78095703125*z4 + 4773.893136160715*z6 - 8632
     $           .101399739584*z8
            r2 = -0.20947265625 + 10893.20654296875*z10 - 8252
     $           .42919921875*z12 + 2448.52294921875*z14 + 21
     $           .99462890625*z2 - 373.90869140625*z4 + 2368.08837890625
     $           *z6 - 7104.26513671875*z8
            result = r1 * z + r2 * Q0
         case (15)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            z10 = z8 * z2
            z12 = z6 * z6
            z14 = z10 * z4
            r1 = 0.31825951825951826 - 19990.82958984375*z10 + 15561
     $           .7236328125*z12 - 4733.81103515625*z14 - 35
     $           .04905234739219*z2 + 621.1384055397727*z4 - 4081
     $           .3972981770835*z6 + 12654.588448660714*z8
            r2 =-3.14208984375 + 24757.28759765625*z10 - 17139
     $           .66064453125*z12 + 4733.81103515625*z14 + 124
     $           .63623046875*z2 - 1420.85302734375*z4 + 7104
     $           .26513671875*z6 - 18155.34423828125*z8
            result = r1 + r2 * z * Q0
         case (16)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            z10 = z8 * z2
            z12 = z6 * z6
            z14 = z10 * z4
            z16 = z8 * z8
            r1 = 4.895771676917917 - 45703.721282958984*z10 + 32446
     $           .329803466797*z12 - 9171.758880615234*z14 - 202
     $           .50454968005627*z2 + 2395.7503079501066*z4 - 12383
     $           .232080368769*z6 + 32610.860181535994*z8
            r2 = 0.196380615234375 - 45388.360595703125*z10 + 55703
     $           .89709472656*z12 - 35503.582763671875*z14 + 9171
     $           .758880615234*z16 - 26.707763671875*z2 + 592
     $           .0220947265625*z4 - 4972.985595703125*z6 + 20424
     $           .762268066406*z8
            result = r1 * z + r2 * Q0            
         case default
            
c$$$            if (z-1.0d0.ge.1.0D-6) then
            if (z-1.0d0.ge.0.021D0) then
               dLla=dble(Lla)
               a=dLla/2.d0+1.d0
               b=(dLla+1.d0)/2.d0
               c=dLla+1.5d0      
               ABCF = 1.0D0
               XI = 1.0D0
               S1 = 0.0D0
               s2=1.d0
               i=0
               zson = 1d0/(z*z)
               do while(abs((s1-s2)/(s1+s2)).gt.1d-6)
c$$$               do while(abs(s2-s1).gt.1d-7)
                  i=i+1
                  s2=s1
                  S1 = S1 + ABCF*XI
                  SI = DBLE(I-1)
                  ABCF = ABCF*(A+SI)*(B+SI)/((C+SI)*(SI+1.0D0))
! Below seems to be slower                  
c$$$                  lp2i = Lla+2*i
c$$$                  ABCF = ABCF*lp2i/(2*i)*(lp2i-1)/(lp2i+Lla+1)
                  XI = XI*zson
c$$$                  XI = XI/z/z
               end do
c$$$               print*,i,(A+SI)*(B+SI)/((C+SI)*(SI+1.0D0)),
c$$$     >            lp2i/(2.0*i)*(lp2i-1)/(lp2i+Lla+1)
c$$$               result=Qlfactor(Lla)/z**(Lla+1)*s1
c$$$               if (s1.lt.1e-10.or.Qlfactor(Lla).lt.1e-10) print*,
c$$$     >              'z,Lla,s1,Qlfactor(Lla):',z,Lla,s1,Qlfactor(Lla)
               result=exp(log(Qlfactor(Lla)) - log(z)*(Lla+1) + log(s1))
c$$$               print*,'i,s1,s2,Qlfactor:',i,s1,s2,Qlfactor(Lla)
            else
               call LQNB(Lla,z,QN)               
c$$$               print*,'z,Lla,QN:',z,Lla,QN(Lla)
               result = QN(Lla)
            end if
         end select        
      else           
!         result = FlegQ(min(Lla,20), z) !  series expansion for large z
         if (Lla.le.20) then
            result = FlegQ(Lla, z) !  series expansion for large z                        
c$$$            print*,'z,Lla,FleQ:',z,Lla,FlegQ(Lla, z)
         else
            stop 'must have Lla<=20 to call FlegQ'
         endif
      end if                    ! z
!     resultsp = real(result)      
      return
      end subroutine funleg     

      subroutine funlegnum(z,Lla,Q0,result,qa,qq) !result is Q_L
c$$$      use apar, only : alkali
c$$$      include 'par.f'
c$$$      include 'par.pos'
c$$$      parameter( maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
c$$$      real zasym
c$$$      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
c$$$     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
c$$$      common/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igz
c$$$      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
c$$$     >     Qlfactor(0:maxl)
      
*     arg of Legendre function is always .ge. 1
*     min[z]=1 at pp=qa
*     Mathematica cannot handle this region.
      
          
*     for hydrogen
      call funleg(z, Lla, Q0, result)

c$$$*     now this part is computed separately: 
c$$$      
c$$$*     for alkalies: V_{positron-core} = V_Coulomb + rest         
c$$$      if (alkali) then
c$$$         resz=0.d0
c$$$         do iz=1,igz
c$$$            deltaq=sqrt(qa*qa+qq*qq-2.d0*qa*qq*xgz(iz))
c$$$            resz=resz + positroncoreV(deltaq) * pl(Lla,iz)
c$$$         end do
c$$$         result=result+qa*qq*resz
c$$$      end if
      
      return
      end subroutine funlegnum            

*------------------------------------------------------

      subroutine polleg(limit,igz,igp)
      include 'par.f'
      include 'par.pos'
      parameter (maxl=8*ltmax)
      implicit real*8 (a-h,o-z)       
      common/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igza
      common/gauszR/xgzR(igzm),wgzR(igzm),plR(0:maxl,igzm),igzR
       common/gausp/xgp(igpm),wgp(igpm),igpa
      dimension endpts(2),b(igz)
      real*8 pl8(0:limit)

      real*8, dimension (igz) :: xga, wga
      real*8, dimension (igp) :: xgb, wgb
      real*8 :: x1, x2, z

      if (limit.gt.maxl) stop 'polleg: increase maxl'
***** for z integral
      x1=-1.d0; x2=1.d0

C     gauss-legendre quadrature is used
C     gauleg utilizes either cgqf or the code from Numerical recepies

C     print*, '   (1.1) polleg'
      call gauleg(x1,x2,xga,wga,igz)

      do ig = 1,igz
         z = xga(ig)
         xgz(ig) = z
         wgz(ig) = wga(ig)
         pl8(0)  = 1.0d0
         pl(0,ig) = 1.d0*wga(ig)
         pl(1,ig) = z*wga(ig)
         pl8(1)=z
         do i=2,limit
            pl8(i)=((2*i-1)*z*pl8(i-1)-(i-1)*pl8(i-2))/dble(i)
            pl(i,ig)=pl8(i)*wga(ig)
         end do
      enddo
cRRRR***** for z integral
      zeps=0.03d0
      x1=-1.d0; x2=1.d0-zeps
      igz1=igz/4*3
      igz2=igz/4
C     gauss-legendre quadrature is used
C     gauleg utilizes either cgqf or the code from Numerical recepies

C     print*, '   (1.1) polleg'
      call gauleg(x1,x2,xga,wga,igz1)
      x1=x2; x2=1.d0

      do ig = 1,igz1
         z = xga(ig)
         xgzR(ig) = z
         wgzR(ig) = wga(ig)
         pl8(0)  = 1.0d0
         plR(0,ig) = 1.d0*wga(ig)
         plR(1,ig) = z*wga(ig)
         pl8(1)=z
         do i=2,limit
            pl8(i)=((2*i-1)*z*pl8(i-1)-(i-1)*pl8(i-2))/dble(i)
            plR(i,ig)=pl8(i)*wga(ig)
         end do
      enddo

      call gauleg(x1,x2,xga,wga,igz2)
         ig=0
      do ig2 = igz1+1,igz
         ig=ig+1
         z = xga(ig)
         xgzR(ig2) = z
         wgzR(ig2) = wga(ig)
         pl8(0)  = 1.0d0
         plR(0,ig2) = 1.d0*wga(ig)
         plR(1,ig2) = z*wga(ig)
         pl8(1)=z
         do i=2,limit
            pl8(i)=((2*i-1)*z*pl8(i-1)-(i-1)*pl8(i-2))/dble(i)
            plR(i,ig2)=pl8(i)*wga(ig)
         end do
      enddo
 

c$$$      do ig = 1,igz
c$$$         print*,'   ' ,ig, xgz(ig),wgz(ig)
c$$$      end do
            
***** pl returns Leg. pols. * Gauss weights (for z integral)

***** for p integral
      x1=0.d0; x2=1.d0
      
      call gauleg(x1,x2,xgb,wgb,igp)
      
      do ig=1,igp
         xgp(ig)=xgb(ig)
c         xgp(igp2-ig+1)=1.d0/xgb(ig)
         wgp(ig)=wgb(ig)
c         wgp(igp2-ig+1)=wgb(ig)/xgb(ig)/xgb(ig)
      enddo
c$$$      print*, '   (1.2) polleg'      
c$$$      do ig = 1,igp
c$$$         print*,'   ' ,ig, xgp(ig),wgp(ig)
c$$$      end do
c$$$      stop      
      return
      end


*------------------------------------------------------

      subroutine sumk1an(x,ll,Nl,a,sk1) ! alternative entry !
      include 'par.f'
      implicit real*8 (a-h,o-z)
      dimension a(ncmax)
c
c this is the analytical form for sumk1 when l=0 only
c
      select case(ll)
      case(0)    ! for s states

      phi=acos(x)
      sk1 = 0.
      do k=1,Nl
         sk1=sk1+a(k)*((k+1.d0)*sin(k*phi)- k*sin((k+1.d0)*phi))
      enddo
      sk1=sk1/8.0d0/cos(phi/2.d0)/(sin(phi/2.d0))**3

c      print*,' x = ',x,' anal sk1 = ',sk1

      case default
         stop ' sumk1an: statement is missing for ll>0'
      end select
      return
      end

*------------------------------------------------------

      subroutine sumk0an(x,ll,Nl,a,sk0) ! alternative entry !
      include 'par.f'
      implicit real*8 (a-h,o-z)
      dimension a(ncmax)
c
c this is the analytical form for sumk0 when l=0 only
c
      select case(ll)
      case(0)    ! for s states

      phi=acos(x)
      cos2=cos(phi/2.d0)
      sk0 = 0.
      do k=1,Nl
         sk0=sk0+a(k)*(cos2- cos((k+0.5d0)*phi))
      enddo
      sk0=sk0/2.0d0/sin(phi/2.d0)/sin(phi)

c      print*,' x = ',x,' anal sk0 = ',sk0

      case default
         stop ' sumk0an: statement is missing for ll>0'
      end select
      return
      end


*------------------------------------------------------

      subroutine dirvswave(result,qfs,qis,nfa,nia,result2,ier)

* calculates s-wave (J=Llb=Lla=lb=la=lam=0).
* Clear that only atom --> atom term survive
*
*     Note about precision loss: pp. 1 2 3 4 5 9 and 8
*     are reliable (real*8)
*     pp. 6 and 7 are not (real*8)
      
      include 'par.f'
      parameter( maxl=8*ltmax)
      implicit real*8 (a-h,o-z)
c$$$      common/const/pi
      common/factors/factrl(0:maxl),sqrfct(0:maxl),hat(0:maxl),
     >   Qlfactor(0:maxl)
      common/fik/fik(0:19,0:19)
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)
      real result,qfs,qis,result2
      dimension sumi(2:40),trig(39)
      ier=0
      result=0.
      pi = acos(-1d0)
      Nl=npsstates(1,0)
      if(qfs.lt.0. .or. qis.lt.0. .or. Nl.eq.0) return
      qf=dble(qfs)
      qi=dble(qis)

* an alternative entry
* non-optimized, very slow, for testing purposes      
c      call dirvsnum(res2nu,qf,qi,nfa,nia)
c      result2=real(res2nu)
c      result=result2
c      return
* end of alternative entry     

c      pause' privet'
      rlam=rlambda(1,0)
      xp = (qf+qi)/rlam
      xm = (qf-qi)/rlam
      Rplus =sqrt(1.d0+xp*xp)
      Fplus =atan(xp)
      Rminus=sqrt(1.d0+xm*xm)
      Fminus=atan(xm)
      Rlog=log(Rplus/Rminus)
      summ=0.d0
      do m=1,2*Nl-1
         dm=dble(m)
* No 1         
c         if(abs((cos(dm*Fminus)/Rminus**m - cos(dm*Fplus)/Rplus**m)/
c     >      (cos(dm*Fminus)/Rminus**m + cos(dm*Fplus)/Rplus**m))
c     >      .lt. 1d-9) then
c            pause' problem No 1'
c         endif
* end No 1         
         trig(m)=(cos(dm*Fminus)/Rminus**m - cos(dm*Fplus)/Rplus**m)/dm
* No 2         
c         if(abs((summ+dble(2*Nl-m)*trig(m))/(summ-dble(2*Nl-m)*trig(m)))
c     >      .lt. 1d-9) then
c            pause' problem No 2'
c         endif
* end No 2         
         summ=summ+dble(2*Nl-m)*trig(m)
      end do
* No 3         
c      if(abs((dble(2*Nl) * Rlog + summ)/
c     >   (dble(2*Nl) * Rlog - summ)).lt. 1d-9) then
c         pause' problem No 3'
c      endif
* end No 3         
      sumi(2*Nl)= dble(2*Nl) * Rlog + summ
      do i=2,2*Nl-1
         summ=0.d0
         do m=1,i-1
* No 4         
c         if(abs((summ+dble(i-m)*trig(m))/(summ-dble(i-m)*trig(m)))
c     >      .lt. 1d-9) then
c            pause' problem No 4'
c         endif
* end No 4         
            summ=summ+dble(i-m)*trig(m)
         end do
* No 5         
c      if(abs((dble(i) * Rlog + summ)/(dble(i) * Rlog - summ))
c     >      .lt. 1d-9) then
c         pause' problem No 5'
c      endif
* end No 5         
         sumi(i)= dble(i) * Rlog + summ
      end do
      sumk1=0.d0
      do k1=1,Nl
         sumk2=0.d0
         do k2=1,Nl
            sumi1=0.d0
            do i1=0, k1-1

               sumi2a=0.d0
               do i2=0,k2-1,2
                  sumi2a=sumi2a+sumi(i1+i2+2)*fik(i1,i2)/
     >               factrl(k2-1-i2)/factrl(k1-1-i1)
c                  print*,'  sumi2a=',sumi2a,i2,
c     >               sumi(i1+i2+2),
c     >               fik(i1,i2),
c     >               factrl(k2-1-i2)*factrl(k1-1-i1)
               end do
               sumi2b=0.d0
               do i2=1,k2-1,2
                  sumi2b=sumi2b+sumi(i1+i2+2)*fik(i1,i2)/
     >               factrl(k2-1-i2)/factrl(k1-1-i1)
c                  print*,'  sumi2b=',sumi2b,i2,
c     >               sumi(i1+i2+2),
c     >               fik(i1,i2),
c     >               factrl(k2-1-i2)*factrl(k1-1-i1)

               end do
c               print*
* No 6         
c      if(abs((sumi2a-sumi2b)/(sumi2a+sumi2b)).lt.1d-9) then
c         ier=1
c         print'('' problem No 6: '',e23.16,'' -/+'')',sumi2a
c         print'(''               '',e23.16,'' ==>'',e24.16)',
c     >      sumi2b,sumi2a-sumi2b
c       print*,'k1 i1 k2:',k1,i1,k2
c         pause' problem No 6'
c      endif
* end No 6
    
                  sumi2= sumi2a-sumi2b
               if(real(i1/2).ne.real(i1)/2) sumi2=-sumi2
               
* No 7         
c      if(abs((sumi1+sumi2)/
c     >   (sumi1-sumi2)) .lt. 1d-9) then
c         ier=1
c         print'('' problem No 7: '',e23.16,'' -/+'')',sumi1
c         print'(''               '',e23.16,'' ==>'',e24.16)',
c     >    sumi2,sumi1+sumi2
c      pause' problem No 7'
c      endif
* end No 7
               sumi1=sumi1+sumi2
c               print*,'sumi1=',sumi1,i1
            end do
* No 8         
c      if(abs((sumk2+factrl(k2+1)*Cknd(k2,nfa,0)*sumi1)/
c     >   (sumk2-factrl(k2+1)*Cknd(k2,nfa,0)*sumi1)) .lt. 1d-9) then
c         pause' problem No 8'
c      endif
* end No 8         
            sumk2=sumk2+factrl(k2+1)*Cknd(k2,nfa,0)*sumi1
         end do
* No 9         
c      if(abs((sumk1+factrl(k1+1)*Cknd(k1,nia,0)*sumk2)/
c     >   (sumk1-factrl(k1+1)*Cknd(k1,nia,0)*sumk2)) .lt. 1d-9) then
c         pause' problem No 9'
c      endif
* end No 9         
         sumk1=sumk1+factrl(k1+1)*Cknd(k1,nia,0)*sumk2
      end do
* after multiplying by normalisation coef 2/pi and also by qf*qi one gets:
      result=real(sumk1/(rlam*pi))
c      if(ier.ne.0) then
c$$$         print'(''analytical res='',e23.10,'' used part='',e17.7)',
c$$$     >      sumk1/(rlam*pi),result
c$$$         print'(''numerical  res='',e23.10,'' used part='',e17.7)',
c$$$     >      res2nu,result2
c$$$         pause
c$$$c         print'(''qf r8='',e23.16,'' qf r4='',e17.7)',
c$$$c     >      qf,qfs
c$$$c         print'(''qi r8='',e23.16,'' qi r4='',e17.7)',
c$$$c     >      qi,qis
c         print*,' nfa nia:',nfa,nia
c         print*,' qfs qis:',qfs,qis
c      endif
      return
      end

*------------------------------------------------------
!Charlie Subtraction function for Coulomb wave
        subroutine RegCharge(q,dqq,l,eta,result,result2)
        implicit none

      double precision, intent(in) :: q,dqq
      integer, intent(in) :: l
      double precision, intent(in) :: eta
      double precision, intent(out) :: result,result2
      double complex :: value

      ! Other variables

      double precision :: qfactor,Phase
      double precision :: tmp,pi,tmp0,mag
      double complex :: ctmp,ctmp2
      double complex :: g

      interface
        double precision function cwfnSigmaL(lt, etat)
            integer, intent(in) :: lt
            double precision, intent(in) :: etat
        end function cwfnSigmaL
      end interface
        pi=4.d0*datan(1.d0)

        ctmp=dcmplx(1.d0,eta)
            tmp = cwfnSigmaL(l, eta)
            tmp0= cwfnSigmaL(0, eta)
            ctmp= dcmplx(0.d0,tmp0-tmp)
        call cwf_gamma_mag(0,eta,mag)!|gamma(1+in)|*exp(pi|n|/2)
        value=mag*exp(ctmp)

        ctmp=dcmplx(0.d0,eta)
         if (eta.gt.0.d0) then
                if (dqq.gt.0.d0) then
                Phase=exp(-pi*eta)
                else
                Phase=1.d0
                endif
         else
                if (dqq.gt.0.d0) then
                Phase=1.d0
                else
                Phase=exp(pi*eta)
                endif

        endif


        ctmp=dcmplx(0.d0,eta)
        !qfactor=abs(2.d0*q/dqq)
        qfactor=abs((2.d0*q+dqq)/dqq)
        value = value * dcmplx(qfactor, 0.d0)**ctmp
        value=value*Phase
        result=real(value)
        result2=imag(value)
        end
*------------------------------------------------------
!Charlie: Subtraction function
      subroutine SubCharge(q,qa,eta,l,ivalue)
      implicit none

      double precision, intent(in) :: q,qa,eta
      integer, intent(in) :: l
      double precision, intent(out) :: ivalue
      
      double complex :: ctmp, value,f21value,g,z2,ctmp2
      double precision :: z,dqq,qfactor,Phase,pi,tmp,mag,tmp0
        
      interface
        double precision function cwfnSigmaL(lt, etat)
            integer, intent(in) :: lt
            double precision, intent(in) :: etat
        end function cwfnSigmaL
      end interface
        pi=4.d0*datan(1.d0)
        dqq=q-qa

        ctmp = dcmplx(2.d0,-eta)
      z = (q-qa)/(qa+q)
      
      call cwfnf21(dcmplx(1.d0,0.d0),dcmplx(1.d0,-eta),ctmp,
     >    dcmplx(z,0.d0),f21value)
        
        ctmp=dcmplx(1.d0,eta)
            tmp = cwfnSigmaL(l, eta)
            tmp0= cwfnSigmaL(0, eta)
            ctmp= dcmplx(0.d0,tmp0-tmp)
        call cwf_gamma_mag(0,eta,mag)!|gamma(1+in)|*exp(pi|n|/2)
        value=mag*exp(ctmp)

        
        ctmp=dcmplx(-1.d0,eta)
         if (eta.gt.0.d0) then
                if (dqq.gt.0.d0) then
                Phase=exp(-pi*eta)
                else
                Phase=-1.d0
                endif
         else
                if (dqq.gt.0.d0) then
                Phase=1.d0
                else
                Phase=-exp(pi*eta)
                endif

        endif

        if ((mag*Phase).eq.0.d0) then
        f21value=dcmplx(1.d0,0.d0)
        endif


        qfactor=abs(1.d0/z)
        value = -value/ctmp * dcmplx(qfactor, 0.d0)**ctmp
        value=value*Phase*f21value
      ivalue=imag(value)
      end
*------------------------------------------------------

      subroutine dirvsnum(result,qf,qi,nfa,nia)

* calculates s-wave (J=Llb=Lla=lb=la=lam=0)
* numerically in double precision.

      include 'par.f'
      parameter( maxl=8*ltmax)
      parameter(limlst=100, limit=100, maxp1=100)
      implicit real*8 (a-h,o-z)
c$$$      common/const/pi
      common /laguercoeffs/
     >   cknd(ncmax,ncmax,0:lnabmax),rlambda(2,0:lnabmax),
     >   npsstates(2,0:lnabmax)
      common/ckif/cki(20),ckf(20),rlam,Nl
      common/omega/agemo
c      real result,qfs,qis

      dimension alist(limit),blist(limit),chebmo(maxp1,25),elist(limit),
     *  erlst(limlst),ierlst(limlst),iord(limit),nnlog(limit),psum(52),
     *  res3la(3),rlist(limit),rslst(limlst)
      external fr21,fr22

      pi = acos(-1d0)
      result=0.
      Nl=npsstates(1,0)
c      if(qfs.lt.0. .or. qis.lt.0. .or. Nl.eq.0) return
c      qf=dble(qfs)
c      qi=dble(qis)
      rlam=rlambda(1,0)
      do k=1,Nl
         cki(k)=Cknd(k,nia,0)
         ckf(k)=Cknd(k,nfa,0)
      enddo

      omegam = abs(qf-qi)
      omegap = qf+qi

      a=0.d0
      integr=1
      epsabs=1.d-10

      call dqawfe(fr21,a,omegam,integr,epsabs,limlst,limit,maxp1,
     *   res12m,abserr,neval,ier,rslst,erlst,ierlst,lst,alist,blist,
     *   rlist,elist,iord,nnlog,chebmo)

c      if(ier.ne.0) then
c         print*,'dqawfe m: ier=',ier
c      endif

      call dqawfe(fr21,a,omegap,integr,epsabs,limlst,limit,maxp1,
     *   res12p,abserr,neval,ier,rslst,erlst,ierlst,lst,alist,blist,
     *   rlist,elist,iord,nnlog,chebmo)

c      if(ier.ne.0) then
c         print*,'dqawfe p: ier=',ier
c      endif

      integr=2
      if(qf.ge.qi) then
         omega=qf
         agemo=qi
      else
         omega=qi
         agemo=qf
      endif

      call dqawfe(fr22,a,omega,integr,epsabs,limlst,limit,maxp1,
     *   res11,abserr,neval,ier,rslst,erlst,ierlst,lst,alist,blist,
     *   rlist,elist,iord,nnlog,chebmo)

c      if(ier.ne.0) then
c         print*,'dqawfe: ier=',ier
c      endif

      res=2d0*agemo*res11-res12m+res12p
* after multiplying by normalisation coef 2/pi and also by qf*qi one gets:
c      result=real(res*rlam*rlam/pi)
      result=res*rlam*rlam/pi
      return
      end

*------------------------------------------------------

      real*8 function fr21(r2)
      include 'par.f'
      parameter(maxl=8*ltmax)
      parameter(limit=100)
      implicit real*8 (a-h,o-z)
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  rlist(limit)
      external fr1

      fr21=0.d0
      if(r2.eq.0.d0) return

      inf=1
      epsabs=1.d-10
      epsrel=0.d0

      call dqagie(fr1,r2,inf,epsabs,epsrel,limit,fr21,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last)

c      if(ier.ne.0) then
c         print*,'dqagie: ier=',ier
c      endif

      return
      end

*------------------------------------------------------

      real*8 function fr22(r2)
      include 'par.f'
      parameter(maxl=8*ltmax)
      parameter(limit=100)
      implicit real*8 (a-h,o-z)
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  rlist(limit)
      common/omega/agemo
      external fr2

      fr22=0.d0
      if(r2.eq.0.d0) return

      inf=1
      epsabs=1.d-10
      epsrel=0.d0

      call dqagie(fr2,r2,inf,epsabs,epsrel,limit,fr22,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last)

c      if(ier.ne.0) then
c         print*,'dqagie: ier=',ier
c      endif

      fr22=fr22*(sin(agemo*r2)/(agemo*r2))

      return
      end

*-----------------------------------------------------------------------

      real*8 function fr1(r)
      implicit real*8 (a-h,o-z)
      common/ckif/cki(20),ckf(20),rlam,Nl

      arg=rlam*r
      nalfa=2
      pl0=1.d0
      pl1=dble(nalfa+1)-arg
      psii=cki(1)+cki(2)*pl1
      psif=ckf(1)+ckf(2)*pl1
      do k=3,Nl
         n=k-2
         pl2=((dble(2*n+nalfa+1)-arg)*pl1-dble(n+nalfa)*pl0)/dble(n+1)
         psii=psii+cki(k)*pl2
         psif=psif+ckf(k)*pl2
         pl0=pl1
         pl1=pl2
      enddo
      fr1=r*exp(-arg)*psif*psii
      return
      end

*------------------------------------------------------

      real*8 function fr2(r)
      implicit real*8 (a-h,o-z)
      external fr1
      fr2=r*fr1(r)
      return
      end
      
*------------------------------------------------------
* The end of the first part of posVmat.f
* This block has to be compiled
* together with positron.f (the second part)

      

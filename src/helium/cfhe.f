C     Note: There are a number of switches:
c     i_propagate = 1  ! use it to propagate solution from known square-integrable solotiion f_fc(r), note: for this case: i_cexch=1
c     i_cexch = 0  ! use it to smooth core-exchange at small r values
c
c
c     Note: for atoms with core electrons there is a problem of
c     calculating continuum functiions related to a unhom. term. 
c     Due to the oscillatory behaivior of core functions at small r
c     the unhom. term needs to be smoothed for small r values. This is
c     done below for case of Zn by setting unhom. term to zero if it is
c     too small comparing to ucentre(r). It will not necessary work
c     for other atoms. Probably the best way to deal with this problem is
c     to make use of accurate square-integrable function f_fc(r) for
c     small values of r and then propagate to larger r values 
c     using two values of f_fc(r) at some small values of r (0.3 - 2 au),
c     preferably where f_fc(r) has a node.
c
c     Calculation of a continuum wave for He in frozen-core approximation.
c     it is assumed here that frozen-core orbital is an eigenfunction of
c     the one-electron Hamiltonian
      subroutine  contfunc(helium,temp,vdcore,ll,ls,en,iternum,
     >   accuracy,ff1,min1,max1,ff2,min2,max2,phase,sigc,
     >     f_fc,minf_fc,maxf_fc)
      include 'par.f'
      common /meshrr/ nr,gridr(nmaxr,3)
      common/smallr/ formcut,regcut,expcut,fast,match      
      common /double/njdouble,jdouble(20) 
      common/powers/ rpow1(nmaxr,0:ltmax),rpow2(nmaxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(nmaxr,0:lmax)
      common /funLag/  fl(nmaxr,nspmax)
      common /minmaxf/ maxf(nspmax),minf(nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      logical:: helium
      real:: temp(nr)
      real:: vdcore(maxr,0:lamax)
      dimension ff1(nmaxr), ff2(nmaxr), f2tmp(nmaxr), f3tmp(nmaxr)
      real:: f_fc(nmaxr), ff2_ne(nmaxr)
      dimension temp0(nmaxr), templ(nmaxr), fun(nmaxr), ucentr(nmaxr),
     >   unhom(nmaxr)
      complex phase, sigc
      real:: res(nr)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)

c
c     en = q*q/2.0 (au) - is the energy of the continuum wave
c     ll - orbital angular momentum, ls - spin of helium atom.
c     ff1(r) is the frozen-core orbital of He+
c     ff2(r) is the calculated frozen-core orbital of He atom.
c     temp0(r) is the screening potential <ff1|(1/|r1-r2|)|ff1>
c     templ(r) is the calculated exchange 'potential' <ff1|(1/|r1-r2|)|ff2>
c     tmpss is the value of the 'nondiagonal lambda parameter'
c           =<ff1|<ff1|(1/|r1-r2|)|ff1>|ff2>
c
      eta = -(1.0+zasym)/sqrt(2.*en)
      znuc = 2.0 + zasym
      coef = real((-1)**ls)/(2.*ll+1.)
      incon = 0
      res(:) = 0.0
c     
c      do i=1,nr
c         ff1(i) = 4.0*sqrt(2.) * gridr(i,1) * exp(-2.0*gridr(i,1))
c      end do
c      call minmaxi(ff1,nr,i1,i2)
c      max1 = i2
c      min1 = i1
      
c     <ff1|(1/|r1-r2|)|ff1>  - calculated only once
      do i=min1,max1
         fun(i) = ff1(i)*ff1(i)*gridr(i,3)
      end do
      lam=0
      call form(fun,min1,max1,rpow1(1,lam),rpow2(1,lam),
     >   minrp(lam),maxrp(lam),nr,temp0,i10,i20)

      ucentr(:) = 0.0
      do i=i10,i20
         ucentr(i) = 2.0*temp0(i)
      end do
      do i=1,nr
         ucentr(i) = ucentr(i) - 2.0*znuc/gridr(i,1)
         unhom(i) = 0.0
      end do

      if(.not.helium) then
         ucentr(1:nr) = ucentr(1:nr) + 2.*vdcore(1:nr,ll)
      endif

            
c noexchange zero order solution for ff2(r)
      incon = 1
      ecmn=2.*en
      Min2=1
      Max2=nr

      call regularun(ll,ecmn,eta,ucentr,cntfug(1,ll),unhom,
     >   gridr(1,1),nr,jdouble,njdouble,regcut,expcut,
     >   ff2,min2,max2,phase,sigc,incon)
      incon = 0
      
      ff2_ne(1:nr) = ff2(1:nr)

c
c     For heavy atoms use f_fc(r) for small r and propagate it to large r.
      i_propagate = 0  ! for He: change here to 1 if you want to use propagation method
      if(.not. helium) then
         i_propagate = 0   ! change here to 1 if you want to use propagation method
      endif
      if(i_propagate .eq. 1) then
c     Find the "matching point", here it is set to the second maximum of f_fc(r).
c     Look for such maximum only on the first half of the r interval where
c     f_fc(r) is defined to avoid loss of accuracy due to exp. fall-off in f_fc(r)
c     Find extremum after second knot (zero) of f_fc(r)
         n_knots_max = 2
         n_knots = 0
         i_0 = 1
         do i=minf_fc+5,maxf_fc/2
            if(f_fc(i)*f_fc(i+1) .lt. 0.0) then
               i_0 = i          ! found  zero of f_fc(r)
               n_knots = n_knots + 1
               if(n_knots .eq. n_knots_max) exit
            endif         
         enddo
         i_m1 = i_0             ! found a point from which continuum function will be propagated
         ff2(:) = 0.0
         ff2(1:maxf_fc) = f_fc(1:maxf_fc)
         print*, 'i_m1=', i_m1, gridr(i_m1,1), 
     >        maxf_fc/2, gridr(maxf_fc/2,1)
      endif
c
c     


c           
c     Save ff2(r) to check on convergence after it is recalculated.
c     It is done here to avoid if(iter.eq.1) statements in iteration loop.
c     Other way the array f3tmp(i) will not be defined for the first iteration.
      f2tmp(:) = 0.0
      f3tmp(:) = 0.0
      min2tmp = min2
      max2tmp = max2
      do i=min2,max2
         f2tmp(i) = ff2(i)            
      end do
c
c     start iteration loop
      if(helium) isub_iter_max = 1
      ifirst_iter = 4  ! number of default iterations, always must be done.

      rksi_ret = 0.0 ! Check with Dmitry
      rksi_ret_tmp = 0.0
      c_ret = 0.0

      i_cexch = 0  ! use it to smooth core-exchange at small r values
      if(i_propagate .eq. 1)  i_cexch = 1  ! no need to do if small r values are excluded anyway
   


 10   lam = ll
      do iter=1,iternum
         
c     recalculated at each iteration-----------------------------------------
c     save f2tmp(i) to use in mixed ff2(i)
         do i=min2tmp,max2tmp
            f3tmp(i) = f2tmp(i)            
         end do
         
c     save ff2(r) to check on convergence after it is recalculated
         min2tmp = min2
         max2tmp = max2
         do i=min2,max2
            f2tmp(i) = ff2(i)            
         end do
         
c     define mixed ff2(i) to improve  stabilty of the convergence
         min2 = max(min2,min2tmp)
         max2 = min(max2,max2tmp)
c     print*, 'c_ret =', c_ret, rksi_ret, rksi_ret_tmp 
         if(helium) c_ret = 0.5
         do i=min2,max2
            ff2(i) =  (1.0-c_ret)*ff2(i) + c_ret*f3tmp(i)
         end do
         
c     <ff1|(1/|r1-r2|)|ff2>  
         minfun=max(min1,min2)
         maxfun=min(max1,max2)
         do i=minfun,maxfun
            fun(i) = ff1(i)*ff2(i)*gridr(i,3)
         end do
         call form(fun,minfun,maxfun,rpow1(1,lam),rpow2(1,lam),
     >        minrp(lam),maxrp(lam),nr,templ,i1l,i2l)
         
         if(.not.helium) then
            fun(:) = 0.0
            fun(1:nr) = ff2(1:nr)*gridr(1:nr,3)
            call cf_fcexch(nr,fun,res,imax_res,ll)
         endif
         
         coef1 = 2.0*coef
         do i=max(i1l,min1),min(i2l,max1)
            unhom(i) = coef1*templ(i)*ff1(i)
         end do
         if(.not.helium) then
            if(i_cexch .eq. 0) then              
               do i=1,imax_res
                  if(abs(ucentr(i))*1e-8 .lt. abs(res(i))) then
                     unhom(i) = unhom(i) + 2.*res(i)
                  endif
               enddo
            else
               unhom(1:imax_res)=unhom(1:imax_res) + 
     >              2.*res(1:imax_res)
            endif
         endif
         
c     <ff1|<ff1|(1/|r1-r2|)|ff1>|ff2>-nondiagonal lambda term: only for ll=0,ls=0
         if(ll.eq.0.and.ls.eq.0) then
            mini = max(i10,min1,min2)
            maxi = min(i20,max1,max2)
            tmpss = 0.0
            do i=mini,maxi
               tmpss = tmpss + temp0(i)*ff1(i)*ff2(i)*gridr(i,3)
            end do
            tmpss1 = 2.0*tmpss
            do i=min1,max1
               unhom(i) = unhom(i) - tmpss1*ff1(i)
            end do
         end if
         
c     exchange included solution for ff2(r) on the 'iter' step 
         min2 = 1
c     use min2 to propagation starting point i_m1
         if(i_propagate .eq. 1) then
            min2 = i_m1
            incon = -1
         endif
         max2=nr
         call regularun(ll,ecmn,eta,ucentr,cntfug(1,ll),unhom,
     >        gridr(1,1),nr,jdouble,njdouble,regcut,expcut,
     >        ff2,min2,max2,phase,sigc,incon)
         
         
         tmp_ret = 0.0
         tmp_ff2max = 0.0
         rksi_ret_tmp = rksi_ret
         i_ret = max(min2,min2tmp)
         do i=max(min2,min2tmp),min(max2,max2tmp)
            ttt = abs((f2tmp(i)-ff2(i)))
            if(ttt .gt. tmp_ret) then
               i_ret = i
               tmp_ret = ttt
            endif
            ttt = abs(ff2(i))
            if(ttt .gt. tmp_ff2max) then
               i_ff2max = i
               tmp_ff2max = ttt
            endif 
         end do
         tmp_ret = f2tmp(i_ret)-ff2(i_ret)
         rksi_ret = tmp_ret/tmp_ff2max
         if( rksi_ret * rksi_ret_tmp .ge. 0.0 
     >        .or. iter .eq. 1) then
            c_ret = 0.0
         else
            c_ret = rksi_ret_tmp/(rksi_ret_tmp - rksi_ret)
         endif
         
         if(iter .le. ifirst_iter) cycle
c     check on convergence of the calculated function ff2(r)
         naccurcy = 1
         if(abs(rksi_ret) .le. accuracy/100.) then
            naccurcy = 0
c$$$            print*,'abs(rksi_ret), accuracy:',
c$$$     >           abs(rksi_ret),accuracy
         endif         
         if(naccurcy.eq.0) then
c$$$            print*,'l=',ll, 
c$$$     >           ', iterations are convergent on step', iter
c     allow for matching in regularun() routine
            if(incon .eq. -1) then 
               incon = -2
            else
               incon = 1
               min2=1
            endif
            max2=nr
            call regularun(ll,ecmn,eta,ucentr,cntfug(1,ll),unhom,
     >           gridr(1,1),nr,jdouble,njdouble,regcut,expcut,
     >           ff2,min2,max2,phase,sigc,incon)   
            goto 1
         end if         
         
c     end iteration loop 
         
      end do

      tmp=0d0
      if (iter.eq.iternum .and. naccurcy .eq. 1) then ! not converged
         i = min(max2,max2tmp)
         tmp = abs((f2tmp(i)-ff2(i)))/(abs(ff2(i))+accuracy)
         print*,'!!!', tmp,f2tmp(i),ff2(i),i
      endif
      print*,'WARNING: iterations are not convergent; L, e(Ry), Acc',
     >     ll,ecmn,tmp
c     Still  allow for matching in regularun() routine to have correctly normalised output
      if(incon .eq. -1) then 
         incon = -2
      else
         incon = 1
         min2=1
      endif
      max2=nr
      call regularun(ll,ecmn,eta,ucentr,cntfug(1,ll),unhom,
     >   gridr(1,1),nr,jdouble,njdouble,regcut,expcut,
     >   ff2,min2,max2,phase,sigc,incon)
      
 1    continue
      
      idebug = 0
      if(idebug .eq. 1) then
         open(784,file='tmp_contfunc')
         if(.not.helium) then
            fun(:) = 0.0
            fun(1:nr) = ff2(1:nr)*gridr(1:nr,3)
            call cf_fcexch(nr,fun,res,imax_res,ll)
         endif
         write(784,'("# l=",I5)') ll
         do i=1,nr/3,1
            write(784,'(F15.8,6E15.6)') gridr(i,1), 
     >           ff2(i), f_fc(i), ff2_ne(i),res(i),
     >           ucentr(i)+cntfug(i,ll)
         enddo

         close(784)
         read(*,*)
      endif


      return
      end
c------------------------------------------------------------------------------
C  The following routine gets the regular solution of the following
C  differential equation. ECMN is the energy in Rydbergs.
C        2
C       d S   [ LN(LN+1)       V(X)          ]
C       --- - [ --------  + --------- - ECMN ] S = g(x)
C         2   [   X*X           X            ]
C       dX      
C
C  Note that UCENTR(X) = V(X)/X and CNTFUG(X) = LN(LN+1)/X/X, unhom(x)=g(x)
C incom - is the variable that idicate whether matching should be performed
C if incom=0 - no matching, incom=1 matching should be done
      subroutine regularun(ln,ecmn,eta,ucentr,cntfug,unhom,gridx,nx,
     >   jdouble,njdouble,regcut,expcut,reg,jstart,jstop,phase,sigc,
     >   incon)
      include 'par.f'
      
      dimension jdouble(njdouble),ucentr(nx),reg(nx),gridx(nx),
     >   cntfug(nx),f(maxr),g(maxr), unhom(nx)
      complex phase,sigc
      complex*16 xx,eta1,zlmin,cfc(9),cgc(9),cfcp(9),cgcp(9),sig(9)
      common/smallr/ formcut,regor,expor,fast,match
      logical match
      
      wnn=sqrt(abs(ecmn))
      if(incon .eq. 0 .or. incon .eq. 1) then
         reg(1:nx) = 0.0
      endif
      if (ecmn.le.0.0.or.ln.gt.3.or.ucentr(1).ne.0.0) then
         if (jstart.gt.jstop) return

         jstartorig = jstart
         if (ucentr(1) .eq. 0.0 .and. wnn .gt. 2.0) then
            rho = wnn * gridx(jstart)
            acc = 1.0
C  Need to take care that jstart stays less than jstop. For this reason
C  the wnn > 2 condition was added above.
            do while (rho.lt.ln.and.acc.gt.1e-6.and.jstart.lt.jstop)
               reg(jstart) = appf1(ln,wnn,gridx(jstart),acc)
               jstart = jstart + 1
               rho = wnn * gridx(jstart)
            enddo
C  make sure that we are not near a point where DX doubles
            do jd = 2, njdouble - 1
               if (abs(jstart-jdouble(jd)).lt.3) jstart=jdouble(jd)-3
            end do
            jstart = max(jstart,1)
         endif 
         
         if(incon .eq. 0 .or. incon .eq. 1) then
            if (jstart.eq.1) then
               s1 = 0.0
            else
               s1=appf1(ln,wnn,gridx(jstart-1),acc)
            end if 
            s2=appf1(ln,wnn,gridx(jstart),acc)
            do while (s2.eq.0.0.and.jstart.lt.jstop)
               jstart = jstart + 1
               s2=appf1(ln,wnn,gridx(jstart),acc)
            enddo 
         else
c     make sure that we are not near a point where DX doubles
            do jd = 2, njdouble - 1
               if (abs(jstart-jdouble(jd)).lt.3) jstart=jdouble(jd)-3
            end do
            jstart = max(jstart,1)
            s1 = reg(jstart-1)
            s2 = reg(jstart)
         endif
         
         call numerovfun(ln,ecmn,ucentr,cntfug,unhom,gridx,nx,
     >      jdouble,njdouble,s1,s2,reg,jstart,jstop)
         if(incon.eq.0 .or. incon .eq. -1) return
         if(incon .eq. -2) jstartorig = 1
         jstart = jstartorig
         jmatch = min(nx,jstop)
         tmp = (ecmn - ucentr(jmatch)) / ecmn - 1.0
         if (eta.eq.0.0) then
            do while (abs(tmp).lt.sqrt(regcut).and.jmatch.gt.jstart+40)
               jmatch = jmatch - 1
               tmp = (ecmn-ucentr(jmatch)) / ecmn - 1.0
            enddo 
            if (jmatch.eq.min(nx,jstop).and.jmatch.gt.jstart+40.and.
     >         ecmn.gt.1.0) then
               print*,'WARNING: asymptotic potential must be zero',
     >            ucentr(jmatch),ln,ecmn,cntfug(jmatch)
            endif 
         else
            asympot = 2.0*eta*ecmn/wnn/gridx(jmatch)
            tmp=abs(ucentr(jmatch)/asympot-1.0)
            do while (tmp.lt.sqrt(regcut).and.jmatch.gt.jstart+40)
               jmatch = jmatch - 1
               asympot = 2.0*eta*ecmn/wnn/gridx(jmatch)
               tmp = abs(ucentr(jmatch)/asympot-1.0)
            enddo 
c$$$            if (jmatch.eq.min(nx,jstop).and.jmatch.gt.jstart+40
c$$$     >         .and.ecmn.gt.1.0)
c$$$     >         print*,'WARNING: asymptotic potential must be Coulomb'
         endif             

         if (match.or.ucentr(1).ne.0.0.or.ecmn.lt.0.0) then
c$$$            print*,'JMATCH, RMATCH,U',jmatch,gridx(jmatch),
c$$$     >         ucentr(jmatch)
            if (abs(ecmn-34.527).lt.-1e-2) then
               do jm = jmatch-50, jmatch+400
                  do j = jstart, jstop
                     f(j) = reg(j)
                  enddo 
                  jmat = jm
                  call matchcc(ln,eta,ecmn,jstart,jstop,gridx,jmat,f,
     >               phase,sigc,j1,j2)
                  print*,j1,j2,jmatch,phase,ucentr(jmat)
               enddo
            else
               call matchcc(ln,eta,ecmn,jstart,jstop,gridx,jmatch,reg,
     >            phase,sigc,j1,j2)
            endif                
            if (eta.eq.0.0.and.ln.le.3.and.ucentr(1).ne.0.0) then
               call ricbessel(wnn,ln,gridx,jstop,jstart,f)
               call ricbessel(wnn,-ln-1,gridx,jstop,jstart,g)
               rph = real(phase)
               aph = aimag(phase)
               do j = jmatch, jstop
                  reg(j) = rph * f(j) + aph * g(j)
               enddo 
            endif 
         else
            phase = (1.0,0.0)
            sigc = (1.0,0.0)
         endif 
            

         do while (abs(reg(jstart)).lt.regcut)
            jstart = jstart + 1
            if (jstart.gt.jstop) return
         end do
         if (abs(ecmn-0.161599).lt.-1e-3) then
            mode1 = 12
            kfn = 0
            zlmin = cmplx(ln)
            nl = 1
            eta1 = (0.0,0.0)
            do j = jstart, jstop
 10            xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j))
               call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,
     >            kfn,ifai)
               tmp = appf1(ln,wnn,gridx(j),acc)
               write(*,'(10e14.6)') gridx(j), reg(j),
     >            real(cfc(1)),tmp,real(xx)
            enddo
            close(42)
            stop
         endif 

c$$$      if (abs(reg(jstart)/regcut).gt.1e+4.and.jstart.gt.1) then
c$$$         print*,'JSTART should be smaller. JSTART, REG(JSTART), LN,'
c$$$     >      //' ETA, ECMN, JSTOP',jstart,reg(jstart),ln,eta,ecmn,jstop
c$$$      end if 
      else
         phase = (1.0,0.0)
         sigc  = (1.0,0.0)
C  A lot of the sodium results were generated without setting jstop to nx,
C  but used from previous value. This goes wrong for NTYPE=1,2 for l > 0
C  after the negative energy Green's function has been invoked.
         jstart = 1
         jstop = nx
         call ricbessel(wnn,ln,gridx,jstop,jstart,reg)
      end if 
      return
      end
C----------------------------------------------------------------------------
      subroutine numerovfun(ln,en,ucentr,cntfug,unhom,gridx,nx,
     >   jdouble,njdouble,rs1,rs2,reg,jstart,jstop)
      implicit real*8 (a-h,o-z)
      real ucentr(nx),reg(nx),gridx(nx),cntfug(nx),en,rs1,rs2,
     >   appf1, acc, unhom(nx)
      dimension jdouble(njdouble)

      ecmn = en
      x = gridx(jstart)
      dx= gridx(jstart+1) - x
      f2 = ucentr(jstart)+cntfug(jstart)-ecmn
      g2 = unhom(jstart)
      h1=dx*dx
      h2=h1/12d0
      s1 = rs1
      s2 = rs2
      g1 = 0d0 !Check with Dmitry
      
      if (jstart.eq.1) then
C  We get here if the integration is going to start from the first X point.
C  This means that S1 is the solution of the differential equation at X=0
C  S2 is the solution at the first X point. 
         s1=0d0
         t1=0d0
         f1=0d0 !Check with Dmitry
C  LN = 1 is a special case
         if (ln.eq.1) t1=-h1/18d0
         if (ecmn.ne.0d0) t1=t1*ecmn
      else
         j=jstart-1
         f1=ucentr(j)+cntfug(j)-ecmn
         t1=s1*(1d0-h2*f1)
      end if

      reg(jstart) = s2
      t2=s2*(1d0-h2*f2) - h2*g2
      
      istart=2
      do while (jstart.gt.jdouble(istart).and.istart.lt.njdouble)
         istart=istart+1
      end do
      istart=istart-1
      istop=njdouble-1
C  JDOUBLE(ISTART) points to the first doubling of DX that happens after JSTART
C  JDOUBLE(ISTOP) points to the last doubling of DX that happens before JSTOP
      
C    integration loop
      s3 = 0.0
      f3 = 0.0
      g3 = 0.0
      g0 = 0.0
      f0 = 0.0
      s0 = 0.0
      do i=istart,istop
         j1=max(jstart,jdouble(i))+1
         j2=min(jstop,jdouble(i+1))
         do j=j1,j2
            g3 = unhom(j)
            f3 = ucentr(j)+cntfug(j)-ecmn
            t3 = 2.*t2-t1+h1*(f2*s2 + g2)
            s3 = (t3 + h2*g3)/(1d0-h2*f3)
c$$$            if (j.lt.jstart+10) then
c$$$               test=appf1(ln,sqrt(abs(en)),gridx(j),acc)
c$$$               print*,test/s3
c$$$               s3 = test
c$$$            endif 
            reg(j)=s3      

            t1=t2
            t2=t3
            f0=f1
            f1=f2
            f2=f3
            g0 = g1
            g1 = g2
            g2 = g3
            s0=s1
            s1=s2
            s2=s3
         end do
         dx=2d0*dx
         h1=dx*dx
         h2=h1/12d0
         t2=s3*(1d0-h2*f3) - h2*g3
         t1=s0*(1d0-h2*f0) - h2*g0
      end do
      return
      end      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     Note: radial integration weights (gridr(i,3)) are included in functions f1(i) and f2(i)
      subroutine cf_fcexch(nr,f2,res,imax_res,l)
      include 'par.f'
      real f2(nr),res(nr)
      real ps2,psen2,enpsinb,psinb,rpow1,rpow2,cntfug,fun(maxr),
     >   vnl(maxr)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      common /worksp/
     >   ps2(maxr,ncmax),psen2(ncmax),minps2(ncmax),maxps2(ncmax)
      common /psinbc/ enpsinb(nnmax,0:lnabmax),
     >   psinb(maxr,nnmax,0:lnabmax),istoppsinb(nnmax,0:lnabmax)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   minrp(0:ltmax),maxrp(0:ltmax),cntfug(maxr,0:lmax)


      res(:) = 0d0
      imax_res = 0
      if (ntype.eq.-3) return
      do lac = 0, lactop
         do nac = lac + 1, nabot(lac) - 1
C The above form fails for sodium if latop = 0.
C In the form below we assume that there are no d states in the core
            minfun = 1
            maxfun = istoppsinb(nac,lac)
            do i = minfun, maxfun
               fun(i) = psinb(i,nac,lac) * f2(i)
            end do 
            do 20 ilt = -lac, lac, 2
               lt = l + ilt
               if (lt.lt.0.or.lt.gt.ltmax) go to 20
               call cleb(2*lac,2*lt,2*l,0,0,0,c)
               const = - float(2 * lac + 1) * c * c /
     >            float(2 * l + 1)
               call form(fun,minfun,maxfun,rpow1(1,lt),rpow2(1,lt),
     >            minrp(lt),maxrp(lt),maxfun,vnl,i1,i2)
               imax_res =  max(imax_res,min(i2,maxfun))
               do i = i1, min(i2,maxfun)
                  res(i) = res(i) + psinb(i,nac,lac) * vnl(i) * const
               end do
 20         continue 
         enddo 
      enddo 
      return
      end 

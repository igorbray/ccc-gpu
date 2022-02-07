C  The following routine gets the regular solution of the following
C  differential equation. ECMN is the energy in Rydbergs.
C        2
C       d S   [ LN(LN+1)       V(X)          ]
C       --- - [ --------  + --------- - ECMN ] S =0
C         2   [   X*X           X            ]
C       dX      
C
C  Note that UCENTR(X) = V(X)/X and CNTFUG(X) = LN(LN+1)/X/X
      
      subroutine regular(ln,ecmn,eta,ucentr,cntfug,ldw,gridx,nx,
     >   jdouble,njdouble,regcut,expcut,reg,jstart,jstop,phase,sigc)
      include 'par.f'
      
      dimension jdouble(njdouble),ucentr(nx),reg(nx),gridx(nx),
     >   cntfug(nx),f(maxr),g(maxr)
      complex phase,sigc
      complex*16 xx,eta1,zlmin,cfc(9),cgc(9),cfcp(9),cgcp(9),sig(9),ci,
     >   coulphase
      common/smallr/ formcut,regor,expor,fast,match
      logical fast,match
      character file*10
      common /debye/dbyexists, dblp, rmudeb, zdby
      logical::dbyexists

c$$$      if (eta.gt.0.0) stop  'routine REGULAR not set up for +ve ETA'
      s1 = 0.0
      wnn=sqrt(abs(ecmn))
      phase = (1.0,0.0)
      sigc = coulphase(dble(eta),ln)

      do j = 1, nx
         reg(j) = 0.0
      enddo

!      s3=appf1(ln,ecmn,gridx(1),acc)
!      s4 = coulrho(ln,0.0,ecmn,gridx(1),acc)
!      print*,'ecmn,s3,s4',ecmn,s3,s4
      if (ecmn.le.0.0.or.ln.gt.3.or.ucentr(1).ne.0.0) then
         jstart = jfirst1(ln,ecmn,gridx,nx,jdouble,njdouble,regcut)
         if (jstart.ge.nx) return
 20      if (jstart.gt.100) then
            if (eta.ge.0.0) then
               s1=appf1(ln,ecmn,gridx(jstart-1),acc)
            else
               s1 = coulrho(ln,eta,ecmn,gridx(jstart-1),acc)
            endif 
            if (abs(s1).gt.regcut.or.acc.lt.1d-6) then
c$$$     >           .and.wnn*gridx(jstart-1).gt.1.0) then
c$$$               mode1 = 12
c$$$               kfn = 0
c$$$               zlmin = cmplx(ln)
c$$$               nl = 1
c$$$               eta1 = eta
c$$$               xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(jstart-1))
c$$$               call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,
c$$$     >            sig,mode1,kfn,ifai)
c$$$               print'("reducing jstart, ln, ecmn, s1, acc:",
c$$$     >            2i5,1p,3d14.3)',jstart,ln,ecmn,s1,acc
               jstart = jstart - 100
C  make sure that we are not near a point where DX doubles
               do jd = 2, njdouble - 1
                  if (abs(jstart-jdouble(jd)).lt.3) jstart=jdouble(jd)+3
               end do       
               go to 20
            endif 
         endif 
         jstop  = jlast1(ecmn,gridx,nx,jdouble,njdouble,expcut)
c$$$         if (ecmn.gt.200.0) then
c$$$            jstop = nx
c$$$            jstart = nx + 1
c$$$         endif 
         jstartorig = jstart

         if (ucentr(1) .eq. 0.0.and.jstart.ge.1) then
c$$$         if (ucentr(1) .eq. 0.0 .and. wnn .gt. 2.0) then
            rho = wnn * gridx(jstart)
            acc = 1.0
            reg(jstart) = appf1(ln,ecmn,gridx(jstart),acc)
c$$$C  Need to take care that jstart stays less than jstop. For this reason
c$$$C  the wnn > 2 condition was added above.
            do while (rho.lt.ln.and.acc.gt.1e-6.and.jstart.lt.jstop)
c$$$            do while (abs(reg(jstart)).lt.regcut.and.jstart.lt.jstop)
c$$$               if (ln.eq.2.and.jstart.gt.2) return
               jstart = jstart + 1
               rho = wnn * gridx(jstart)
               reg(jstart) = appf1(ln,ecmn,gridx(jstart),acc)
c$$$               print*,'increased jstart:',
c$$$     >            ln,rho,reg(jstart),acc,jstart 
            enddo
C  make sure that we are not near a point where DX doubles
            do jd = 2, njdouble - 1
               if (abs(jstart-jdouble(jd)).lt.3) jstart=jdouble(jd)-3
            end do
            jstart = max(jstart,1)
         endif 

         if (jstart.gt.jstop) return
         if (eta.ne.0.0.and.jstart.eq.1.and.ln.eq.1) jstart = 2 !avoids inaccuracy
         if (jstart.eq.1) then
            s1 = 0.0
         else
            s1=appf1(ln,ecmn,gridx(jstart-1),acc)
c$$$            if (eta.lt.0.0) s1=coulrho(ln,eta,ecmn,gridx(jstart-1),acc)
            if (eta.ne.0.0) s1=coulrho(ln,eta,ecmn,gridx(jstart-1),acc)
            if (s1.eq.0.0) then
               print*,
     >       "STOP: S1=0.0 in numerov.f; ln,jstartorig,jstart,rho,acc:",
     >       ln,jstartorig,jstart,rho,acc
               stop "STOP: cannot have S1=0.0 in numerov.f"
            endif
         end if 
         s2=appf1(ln,ecmn,gridx(jstart),acc)
         if (eta.ne.0.0) s2=coulrho(ln,eta,ecmn,gridx(jstart),acc)
         if (s2.eq.0.0) then
            print*,
     >       "STOP: S2=0.0 in numerov.f; ln,jstartorig,jstart,rho,acc:",
     >       ln,jstartorig,jstart,rho,acc
            stop "STOP: cannot have S2=0.0 in numerov.f"
         endif
c$$$         print*,'jstart,S1,S2:',jstart,s1,s2 !was used to check the importance of jstart=2 above
c$$$  if (eta.lt.0.0) then
c$$$            s2 = coulrho(ln,eta,ecmn,gridx(jstart),acc)
c$$$         elseif (eta.gt.0.0) then
c$$$            mode1 = 2
c$$$            kfn = 0
c$$$            zlmin = cmplx(ln)
c$$$            nl = 1
c$$$            eta1 = eta
c$$$            xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(jstart+1))
c$$$            call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,
c$$$     >         sig,mode1,kfn,ifai)
c$$$            s2 = real(cfc(1))
c$$$            print*,'ECMN,ETA,S2:',
c$$$     >         ecmn,eta,s2!,coulrho(ln,eta,ecmn,gridx(jstart),acc)
c$$$c$$$            s2 = 1e-30 ! this requires more investigation, Igor 30/4/2011
c$$$         endif 
c$$$         if (eta.eq.0.0.and.abs((s2-s2old)/(s2old+1e-10)).gt.1e-4)then
c$$$            print*,'S2,S2OLD,rho,l,jstart',s2,s2old,wnn*gridx(jstart),
c$$$     >         jstart
c$$$         endif 

         if (s1.eq.s2.and.s2.eq.0.0) then
            print*,'S1=S2=0 for ln,ecmn:',ln,ecmn
            jstart = jstop+1
         endif 
         call numerovf(ln,ecmn,ucentr,cntfug,gridx,nx,
     >      jdouble,njdouble,s1,s2,reg,jstart,jstop)
         if (eta.ne.0.0.and.jstart.eq.2.and.ln.eq.1) then
            jstart = 1
            reg(1) = s1
         endif 

c$$$         if (abs(reg(jstop)).gt.1.5) print*,'WARNING reg(j)>>1:',
c$$$     >      reg(jstop),ln,ecmn,jstartorig,jstart
         jstart = jstartorig
         jmatch = min(nx,jstop)
         rho = wnn * gridx(jmatch)
         tmp = (ecmn - ucentr(jmatch)) / (ecmn+1e-10) - 1.0
         if (eta.eq.0.0) then
            do while (abs(tmp).lt.sqrt(regcut).and.rho.gt.1.0)
c$$$            do while (abs(tmp).lt.sqrt(regcut).and.jmatch.gt.jstart+40)
               jmatch = jmatch - 1
               rho = wnn * gridx(jmatch)
               tmp = (ecmn-ucentr(jmatch)) /
     >            (ecmn) - 1.0
            enddo

c$$$            if (jmatch.eq.min(nx,jstop).and.jmatch.gt.jstart+40.and.
c$$$     >         ecmn.gt.1.0.and.abs(ucentr(jmatch)/
c$$$     >         (cntfug(jmatch)+1e-30)).gt.2e-2) then
c$$$               print*,'WARNING: abs(ucentr(jmatch)/cntfug(jmatch)) > ',
c$$$     >            '2e-2',ucentr(jmatch)/(cntfug(jmatch)+1e-30)
c$$$            endif 
         else
!            asympot = 2.0*eta*ecmn/wnn/gridx(jmatch)
            asympot = 2.0*eta*abs(ecmn)/(wnn+1e-10)/gridx(jmatch)
            tmp=abs(ucentr(jmatch)/(asympot+1e-10)-1.0)
c$$$            print*,'jmatch,asympot,tmp,ucentr:',jmatch,asympot,tmp,
c$$$     >         ucentr(jmatch)
            do while (tmp.lt.1e-5.and.rho.gt.0.1.and.jmatch.gt.jstart)
c$$$            do while (tmp.lt.1e-5.and.jmatch.gt.jstart+40)
               jmatch = jmatch - 1
               rho = wnn * gridx(jmatch)
!               asympot = 2.0*eta*ecmn/wnn/gridx(jmatch)
               asympot = 2.0*eta*abs(ecmn)/wnn/gridx(jmatch)
               tmp = abs(ucentr(jmatch)/asympot-1.0)
            enddo 
c$$$            if (jmatch.eq.min(nx,jstop))
c$$$     >         print*,'WARNING: asymptotic potential must be Coulomb',
c$$$     >         'ECMN,jmatch,tmp:',ecmn,jmatch,tmp
         endif

         if (match.or.ln.le.ldw) then !.or.ucentr(1).ne.0.0) then !.or.ecmn.lt.0.0) then
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
               call matchcc(ln,eta,ecmn,jstart,jstop,gridx,jmatch,
     >            reg,phase,sigc,j1,j2)
c$$$         print*,sqrt(ecmn),jmatch,phase,reg(jstart),reg(jstart+1)
               diff = abs(sigc-coulphase(dble(eta),ln))
               if (diff.gt.1e-5) then
                  print*,'Coulphase difference:',sigc,
     >               coulphase(dble(eta),ln)
c$$$                  stop 'Coulphase difference > 1e-5'
               endif 
c$$$               if (ln.gt.ldw.and..not.dbyexists) then
c$$$                  diff = abs(real(phase)-1.0)
c$$$                  if (diff.gt.1e-5) then
c$$$                     print*,'Expected phase to be 1 for L, ECMN:',
c$$$     >                  phase,ln,ecmn
c$$$                     write(file,'(i1,"_",1p,e8.2)') ln,ecmn
c$$$                     mode1 = 12
c$$$                     kfn = 0
c$$$                     zlmin = cmplx(ln)
c$$$                     nl = 1
c$$$                     eta1 = eta
c$$$                     open(57,file=file)
c$$$                     do j = jstart, jstop
c$$$                        xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j))
c$$$                        call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,
c$$$     >                     sig,mode1,kfn,ifai)
c$$$                        write(57,'(5e14.6,i8)') gridx(j), reg(j),
c$$$     >                     real(cfc(1)),
c$$$     >                     coulrho(ln,eta,ecmn,gridx(j),acc),acc,j
c$$$c$$$                        reg(j) = real(cfc(1))
c$$$                     enddo
c$$$                     close(57)
c$$$                  endif 
c$$$               endif 
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
            if (ecmn.eq.0.0) then
               sigc = (1d0,0d0)
            else
               sigc = coulphase(dble(eta),ln)
            endif 
         endif 
            
         do while (abs(reg(jstart)).lt.regcut)
            jstart = jstart + 1
            if (jstart.gt.jstop) return
         end do

         if (match.and.ln.gt.ldw) then
c$$$         if (match.and.ln.gt.ldw.and.eta.ne.0.0) then
            if (abs(aimag(phase)).gt.1e-5)
     >         print*,'Old phase (NOT reset to 1) for L, ecmn:',
     >         phase,ln,ecmn
c$$$            phase = (1.0,0.0)
            mode1 = 12
            kfn = 0
            zlmin = cmplx(ln)
            nl = 1
            eta1 = eta

c$$$            write(file,'(i1,"_",1p,e8.2)') ln,ecmn
c$$$            open(57,file=file)
            acc = 1.0
            nerrs = 0
            do j = jstart, jstop
               if (acc.gt.1e-6) then
                  cfc(1) = cmplx(coulrho(ln,eta,ecmn,gridx(j),acc))
               else 
                  xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j))
                  call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,
     >               mode1,kfn,ifai)
               endif 
c$$$               if (wnn*gridx(j).lt.1.0) then
c$$$                  tmp = appf1(ln,ecmn,gridx(j),acc)
c$$$               else
c$$$                  tmp = 0.0
c$$$               endif 
c$$$               write(57,'(10e14.6)') gridx(j), reg(j),
c$$$     >            real(cfc(1)), tmp
               if (abs((reg(j)-real(cfc(1)))/
     >            (reg(j)+real(cfc(1))+regcut)).gt.1e-2) then
                  print'("ln,R0,k,reg(j),cfc:",i3,f12.3,1p,3e12.3)',
     >                 ln,gridx(j),sqrt(ecmn), reg(j), real(cfc(1))
c$$$                  cfc(1) = abs(real(cfc(1))) * reg(j)/abs(reg(j))
c$$$                  print*,'cfc of wrong sign?:',reg(j-2),reg(j-1),
c$$$     >               reg(j),reg(j+1),reg(j+2)
               endif
               if (abs((reg(j)-real(cfc(1)))/(reg(j)+1e-20)).gt.0.1.
     >            and.abs(reg(j)).gt.0.1)
     >            nerrs = nerrs + 1
c$$$     >            ,print*,'>10% error in reg(j):',ln,ecmn,j,reg(j),
c$$$     >            real(cfc(1))
               reg(j) = real(cfc(1))
            enddo
            if (nerrs.gt.0) print*,'WARNING; L, ECMN, NERRS:',
     >         ln,ecmn,nerrs
c$$$  close(57)
         endif 

c$$$      if (abs(reg(jstart)/regcut).gt.1e+4.and.jstart.gt.1) then
c$$$         print*,'JSTART should be smaller. JSTART, REG(JSTART), LN,'
c$$$     >      //' ETA, ECMN, JSTOP',jstart,reg(jstart),ln,eta,ecmn,jstop
c$$$      end if 
      else
         phase = (1.0,0.0)
         sigc  = (1.0,0.0)
         jstart = 1
         jstop = nx
         call ricbessel(wnn,ln,gridx,jstop,jstart,reg)
      end if 
      return
      end
      
      subroutine ricbesselg(ecmn,ln,gridx,jstop,jstart,reg)
      real gridx(jstop), reg(jstop)

      rhosmall = 10**(-3.0/(ln+0.5))
      wnn = sqrt(abs(ecmn))
      if (ln.eq.0.or.ln.eq.-1) jstart = 1
      j = jstart

      if (ln.gt.0) then
         rk = wnn * gridx(j)
         do while (rk.lt.1d0.and.j.lt.jstop)
            reg(j) = appf1(ln,ecmn,gridx(j),acc)
            j = j + 1
            rk = wnn * gridx(j)
         enddo
      endif 

c$$$      select case (ln)
c$$$         case (1)
c$$$            rk = wnn * gridx(j)
c$$$            do while (rk.lt.rhosmall)
c$$$               reg(j) = appf1(ln,ecmn,gridx(j),acc) !rk * rk / 3d0 - rk * rk * rk * rk / 30d0
c$$$               j = j + 1
c$$$               rk = wnn * gridx(j)
c$$$            enddo
c$$$         case (2)
c$$$            rk = wnn * gridx(j)
c$$$            do while (rk.lt.rhosmall)
c$$$               reg(j) =  rk * rk * rk / 15d0
c$$$               j = j + 1
c$$$               rk = wnn * gridx(j)
c$$$            enddo
c$$$         case (3)
c$$$            rk = wnn * gridx(j)
c$$$            do while (rk.lt.rhosmall)
c$$$               reg(j) =  rk * rk * rk * rk / 105d0
c$$$               j = j + 1
c$$$               rk = wnn * gridx(j)
c$$$            enddo
c$$$      end select
      js = j
      if (ecmn.gt.0.0) then
         select case (ln)
            case (0)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = sin(rk)
               enddo
            case (1)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = sin(rk) / rk - cos(rk)
               enddo
            case (2)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = sin(rk) * (3d0/rk/rk - 1d0) - 3d0*cos(rk)/rk
               enddo
            case (3)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = ((15d0/rk/rk/rk-6d0/rk) * sin(rk) -
     >               (15d0/rk/rk - 1d0) * cos(rk))
               enddo
            case (-1)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = cos(rk)
               enddo
            case (-2)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = sin(rk) + cos(rk) / rk
               enddo
            case (-3)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = 3d0 * sin(rk) / rk + (3d0/rk/rk-1d0)*cos(rk)
               enddo
            case (-4)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = - sin(rk) * (1d0 - 15d0/rk/rk)
     >               - cos(rk) * (6d0/rk - 15d0/rk/rk/rk)
               enddo
            case default
               stop 'ricbesselg not coded for L>3'
         end select
      elseif (ecmn.eq.0.0) then
         if (ln.ge.0) then
            do j = jstart, jstop
               reg(j) = gridx(j)**(ln+1)/(2*ln+1) !check for ln > 0
            enddo
         else 
            do j = jstart, jstop
               reg(j) = gridx(j)**(-(-ln-1))/(2*(-ln-1)+1) !check for ln > 0
            enddo
         endif 
      else                      !negative energies
         select case (ln)
            case (0)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = -(exp(-rk)-exp(rk))/2d0
               enddo
            case (1)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = (exp(-rk)*(rk+1d0)+exp(rk)*(rk-1d0))/(2d0*rk)
               enddo
            case (2)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = -(exp(-rk)*(3d0*(1d0+rk)+rk*rk)-exp(rk)*
     >               (3d0*(1d0-rk)+rk*rk))/(2d0*rk*rk)
               enddo
            case (3)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = (exp(-rk)*(15d0+15d0*rk+6d0*rk*rk+rk**3)+
     >               exp(rk)*(15d0*rk-15d0-6d0*rk*rk+rk**3))/(2d0*rk**3)
               enddo
            case (-1)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = exp(-rk)
               enddo
            case (-2)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = exp(-rk)*(rk+1d0)/rk
               enddo
            case (-3)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = exp(-rk)*(3d0+3d0*rk+rk*rk)/(rk*rk)
               enddo
            case (-4)
               do j = js, jstop
                  rk = wnn * gridx(j)
                  reg(j) = exp(-rk)*(15d0+15d0*rk+6d0*rk*rk+rk**3)/
     >               (rk**3)
               enddo
            case default
               stop 'ricbesselg not coded for L>3'
         end select
      endif
c$$$      if (ln.gt.0) then
c$$$         print"('ln,jstart,js,jstop,ecmn,reg(jstart),error:',
c$$$     >3i4,i6,f9.3,1p,3e13.4,e10.1)",
c$$$     >      ln,jstart,js,jstop,ecmn,reg(jstart),(reg(js)-
c$$$     >      appf1(ln,ecmn,gridx(js),acc))/appf1(ln,ecmn,gridx(js),acc)
c$$$      endif 

c$$$      real rho
c$$$      complex*16 xx,eta1,zlmin,cfc(1),cgc(1),cfcp(1),cgcp(1),sig(1)
c$$$
c$$$      xx = rho
c$$$      eta1 = (0.0,0.0)
c$$$      zlmin = dcmplx(l)
c$$$      nl = 1
c$$$      ifai = 0
c$$$      mode1 = 4
c$$$      kfn = 0
c$$$      call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,kfn,ifai)
c$$$      ricbessel = cfc(1)
      return
      end

      subroutine ricbessel(wnn,ln,gridx,jstop,jstart,reg)
      real gridx(jstop), reg(jstop)
      real*8 rk,rhosmall
      if (ln.gt.0) rhosmall = 10**(-3.0/ln)
      if (ln.eq.0) then
         jstart = 1
         do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = sin(rk)
         enddo
      else if (ln.eq.1) then
         jstart = 1
         do i = jstart, jstop
            rk = wnn * gridx(i)
            if (rk.lt.rhosmall) then
               reg(i) = rk * rk / 3d0 - rk * rk * rk * rk / 30d0
            else 
               reg(i) = sin(rk) / rk - cos(rk)
            endif 
         enddo
      else if (ln.eq.2) then
         jstart = 1
         do i = jstart, jstop
            rk = wnn * gridx(i)
            if (rk.lt.rhosmall) then
               reg(i) = rk * rk * rk / 15d0
            else 
               reg(i) = sin(rk) * (3d0/rk/rk - 1d0) - 3d0 * cos(rk)/rk
            endif 
         enddo
      else if (ln.eq.3) then
         jstart = 1
         do i = jstart, jstop
            rk = wnn * gridx(i)
            if (rk.lt.rhosmall) then
               reg(i) = rk * rk * rk * rk / 105d0
            else 
               reg(i) = ((15d0/rk/rk/rk-6d0/rk) * sin(rk) -
     >            (15d0/rk/rk - 1d0) * cos(rk))
            endif
         enddo 
      else if (ln.eq.-1) then
         do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = cos(rk)
         enddo
      elseif (ln.eq.-2) then
          do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = sin(rk) + cos(rk) / rk
         enddo
      elseif (ln.eq.-3) then
          do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = 3d0 * sin(rk) / rk + (3d0/rk/rk - 1d0) * cos(rk)
         enddo
      elseif (ln.eq.-4) then
          do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = - sin(rk) * (1d0 - 15d0/rk/rk)
     >         - cos(rk) * (6d0/rk - 15d0/rk/rk/rk)
         enddo
      else 
         print*,'Riccati-Bessel functions are not coded for |L| > 3.'
      end if
c$$$      real rho
c$$$      complex*16 xx,eta1,zlmin,cfc(1),cgc(1),cfcp(1),cgcp(1),sig(1)
c$$$
c$$$      xx = rho
c$$$      eta1 = (0.0,0.0)
c$$$      zlmin = dcmplx(l)
c$$$      nl = 1
c$$$      ifai = 0
c$$$      mode1 = 4
c$$$      kfn = 0
c$$$      call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,kfn,ifai)
c$$$      ricbessel = cfc(1)
      return
      end

      subroutine numerovf(ln,en,ucentr,cntfug,gridx,nx,
     >   jdouble,njdouble,rs1,rs2,reg,jstart,jstop)
      implicit real*8 (a-h,o-z)
      real ucentr(nx),reg(nx),gridx(nx),cntfug(nx),en,rs1,rs2,
     >   appf1, acc
      dimension jdouble(njdouble)

      if (jstart.ge.jstop) return
      ecmn = en
      x = gridx(jstart)
      dx= gridx(min(jstart+1,nx)) - x
      f1 = 0.0
      f2 = ucentr(jstart)+cntfug(jstart)-ecmn
      h1=dx*dx
      h2=h1/12d0
      s1 = rs1
      s2 = rs2
      s3 = 0.0
      f3 = 0.0
      s0 = 0.0
      f0 = 0.0
      if (jstart.eq.1) then
C  We get here if the integration is going to start from the first X point.
C  This means that S1 is the solution of the differential equation at X=0
C  S2 is the solution at the first X point. 
         s1=0d0
         t1=0d0
C  LN = 1 is a special case
         if (ln.eq.1) t1 = h1/18d0     !there was a - sign here until 11/04/2016.
c$$$  if (ln.eq.1) t1=-h1/18d0 !seems to affect only -ve energies
         if (ecmn.ne.0d0) t1=t1*ecmn
      else
         j=jstart-1
         f1=ucentr(j)+cntfug(j)-ecmn
         t1=s1*(1d0-h2*f1)
      end if

      reg(jstart) = s2
      t2=s2*(1d0-h2*f2)
      
      istart=2
      do while (jstart.gt.jdouble(istart).and.istart.lt.njdouble)
         istart=istart+1
      end do
      istart=istart-1
      istop=njdouble-1
C  JDOUBLE(ISTART) points to the first doubling of DX that happens after JSTART
C  JDOUBLE(ISTOP) points to the last doubling of DX that happens before JSTOP

C    integration loop
      do i=istart,istop
         j1=max(jstart,jdouble(i))+1
         j2=min(jstop,jdouble(i+1))
         do j=j1,j2
            f3 = ucentr(j)+cntfug(j)-ecmn
            t3 = 2.*t2-t1+h1*f2*s2
            s3 = t3/(1d0-h2*f3)
c$$$            if (j.lt.jstart+10) then
c$$$               test=appf1(ln,en,gridx(j),acc)
c$$$               print*,j,test/s3,f3,t3
c$$$c$$$               s3 = test
c$$$            endif 
            reg(j)=s3      

            t1=t2
            t2=t3
            f0=f1
            f1=f2
            f2=f3
            s0=s1
            s1=s2
            s2=s3
         end do
         dx=2d0*dx
         h1=dx*dx
         h2=h1/12d0
         t2=s3*(1d0-h2*f3)
         t1=s0*(1d0-h2*f0)
      end do
      return
      end

      function coulrho(l,seta,ecmn,x,acc)
      include 'par.f'
!      include 'par.for'
      parameter (nmax=1000)
      implicit real*8 (a-h,o-z)
      real seta, ecmn, x, coulrho, acc
      dimension C(0:lmax),A(nmax,0:lmax)
      logical converged
      real*8 ffgg, factl, api, x2, w2, res
      common /cnsts2/ factl(0:2*lcoul),api,x2(63,5),w2(63,5),res(63)
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
      data delta /1d-10/

      coulrho = 0.0
c$$$      if (seta.gt.30.0) return ! failed for small negative ecmn
      if (ecmn.eq.0.0) then
         fact = 2**l/exp(factl(2*l+1))
         if (seta.eq.0.0) then
            coulrho = sqrt(2d0*pi) * x ** (l+1)*fact
         else
            coulrho = sqrt(2d0*pi/abs(seta)) * (-seta * x) ** (l+1)*fact
         endif
         return
      endif 
      eta = seta
      if (ecmn.lt.0.0) eta = - abs(seta) !test for repulsive cases
      if (eta.eq.0d0) then
         C(0) = 1d0
      else
c$$$         print*,'2 pi eta:',2d0*pi*eta
c$$$         tmp = exp(2d0*pi*eta)
c$$$         print*,'exp(2 pi eta):',tmp
         if (eta.gt.100.0) then
            coulrho = 0.0
            return
         else 
            C(0) = sqrt(2d0*pi*eta/(exp(2d0*pi*eta)-1d0))
         endif 
      endif
      do lp = 1, l
         C(lp) = C(lp-1)*sqrt(lp*lp+eta*eta)/lp/(2*lp+1)
      enddo 
      rho = sqrt(abs(ecmn))*x
      s = ecmn/abs(ecmn)
      converged = .false.
      A(l+1,l) = 1d0
      A(l+2,l) = eta/(l+1)
      n = l+1
      sum = 0d0
      termp = 0.0
      termn = 0.0
      rhop = rho**l
      do while (.not.converged.and.n.lt.nmax-2)
         rhop = rhop * rho
         term = A(n,l) * rhop * C(l)
         if (term.ge.0d0) then
            termp = termp + term
         else
            termn = termn + term
         endif 
         if (termp.gt.1d30.or.
     >      abs((termp+termn)/(termp-termn)).lt.delta) then
            coulrho = sum
c$$$            print*,'No convergence;l,n,termp,termn',l,n,termp,termn
            n = nmax-2
         endif 
         sum = sum + term
         converged = abs(term/sum).lt.delta.and.abs(A(n,l)).gt.0d0
c$$$         print*,'l,n,rho,termp,termn,sum:',l,n,rho,termp,termn,sum
         n = n + 1
         A(n+1,l) = (2d0*eta*A(n,l)-s*A(n-1,l))/(n-l)/(n+1+l)
      enddo
      if (n.eq.nmax.or.termp.eq.-termn) then
         acc = 1e-30
      else
         acc = abs(termn/termp+1d0)
      endif 
      coulrho = sum
c$$$      print*,'ecmn,x,seta,C(0),coulrho:',ecmn,x,seta,C(0),coulrho
c$$$      if (ecmn.eq.1d-10) then
c$$$         coulrho = sum / 1d-5**(l+1)
c$$$         ecmn = 0.0
c$$$      endif 
c$$$      if (abs(sum) .le. 1e-30) then
c$$$         stop 'COULRHO <= 1e-30'
c$$$      endif 
      return
      end

C  Riccati-Bessel function for small rho. Same as spherical Bessel * rho that works for positive and negative energies
      function appf1(ln,ecmn,x,acc)
      implicit real*8 (a-h,o-z)
      real appf1, ecmn, x, acc
      w = sqrt(abs(ecmn))
      s = ecmn / abs(ecmn+1d-30)
      acc = 0.0
      appf1 = 0d0
      if (w.ne.0.0) then
C  We use the expansion for the Riccati-Bessel function j(l,rho)
C  j(l,rho) = rho^(l+1) sum(k=0,oo) (-1)^k 2^k (rho/2)^(2k)/k!/(2(k+l)+1)!!
         kmax = 1000
         rho=w*x
         k = 1
         if (rho.gt.1e3) return
         ff = 1d0
         do i=1,ln
            ff = ff * rho / float(2*i+1)
         end do
         sum = ff
         summ = 0d0
         sump = ff
         zo2k = 1d0
         sumold = 0d0
         if (abs(ff).lt.1d-300) return
         do while (k.lt.kmax.and.abs(sumold/sum-1d0).gt.1d-6)
            zo2k = zo2k  * rho * rho / 2d0 / float(k) /
     >         float(2*(ln+k)+1)
            sumold = sum
            if (mod(k,2).eq.0) then
               sump = sump + zo2k * ff
            else
               summ = summ + zo2k * ff
            endif 
            sum = sump - s*summ
c$$$            print '(2f9.3,i4,4e24.14)', w,x,k,sump,-s*summ, sum
            k = k + 1
            if (abs(s*summ/sump-1d0).lt.1d-12.and.(k.ge.kmax
     >           .or.abs(sum).lt.1d-300)) then
c$$$               k = kmax
               acc = 0.0
               appf1 = 0.0
c$$$               print'("Precision loss in APPF1 for ln, k, ecmn:",
c$$$     >            2i5,1p,e11.3,3e18.10)',ln,k,ecmn,-s*summ,sump,sum
               return
            endif 
         enddo
         acc = abs(s*summ/sump-1d0)
         appf1 = rho * sum
c$$$         print*,rho,k,s*summ,sump
c$$$         if (k.eq.kmax) acc = 0.0
c$$$     >      print'("Possible precision loss in APPF1;result and error",
c$$$     >      1p,2e14.4)',sum,acc
      else
C  The following is the limiting case of the Riccati-Bessel/w**(ln+1) for w=0.0
         iprod = 2 * ln + 1
         iterm = 2 * ln + 1
         do while (iterm.gt.3)
            iterm = iterm - 2
            iprod = iprod * iterm
         enddo 
         appf1=x**(ln+1)/float(iprod)
      end if 
      return
      end

C  This function returns the index of the first X for which the regular
C  solution is < REGSTART.
      function jfirst1(l,ecmn,gridx,nx,jdouble,njdouble,regcut)
      dimension gridx(nx),jdouble(njdouble)
      w=sqrt(abs(ecmn))
      if (w.ne.0.0) then
         tmp=0.0
         do i=1,l
            tmp=tmp+log(float(2*i+1))
         end do
      else
         tmp=log(float(2*l+1))
         w=1.0
      end if
 10   xstart = 0
      if (regcut.ne.0)
     >     xstart=min(exp((tmp+log(regcut))/(l+1))/w,gridx(nx))
      j=max(int(xstart*nx/gridx(nx)),1)
      j = min(j,nx)
      do while (gridx(j).lt.xstart.and.j.lt.nx)
         j=j+1
      end do
c$$$      s1 = appf1(l,ecmn,gridx(j),acc)
c$$$      if (acc.lt.1e-10) then
c$$$         regcut = regcut / 10.0
c$$$         print '("regcut redefined. S1, REGCUT, acc:",1p,3e13.4)',
c$$$     $      s1,regcut,acc
c$$$         go to 10
c$$$      endif 

C  make sure that we are not near a point where DX doubles
      do jd = 2, njdouble - 1
         if (abs(j-jdouble(jd)).lt.3) j = jdouble(jd) + 3
      end do       
c$$$      if (j.eq.nx) j = nx + 1
      jfirst1=j
      return
      end

C  This function returns the index of the last X for
C  which EXP(-WNN*X) < EXPCUT
      function jlast1(ecmn,gridx,nx,jdouble,njdouble,expcut)
      dimension gridx(nx),jdouble(njdouble)
      if (ecmn.ge.-0.005) then
         jlast1 = nx
      else
         wnn = sqrt(-ecmn)
         xmax=-log(expcut)/wnn
         j=min(int(xmax*nx/gridx(nx)),nx)
         do while(gridx(j).gt.xmax.and.j.gt.1)
            j=j-1
         end do
         do while (gridx(j)**(1.0/wnn)*exp(-gridx(j)*wnn).gt.expcut
     >      .and.j.lt.nx)
            j=j+1
         end do 

C  make sure that we are not near a point where DX doubles
         do jd = 2, njdouble - 1
            if (abs(j-jdouble(jd)).lt.3) j = jdouble(jd) + 3
         end do       

         jlast1=j
      end if 
      return
      end
      
C  The following routine solves the differential equation below backwards
C  from S(JSTOP)=S3 and S(JSTOP-1)=S2
C        2
C       D S   [ LN(LN+1)       V(X)          ]
C       --- - [ --------  + --------- - ECMN ] S =0
C         2   [   X*X           X            ]
C       DX      
C
C  Note that V(X)/X = UCENTR(X)
      subroutine numerovb(ln,ecmn,ucentr,cntfug,gridx,nx,
     >   jdouble,njdouble,s2,s3,chi,jstart,jstop)

      dimension ucentr(nx),cntfug(nx),gridx(nx),jdouble(njdouble),
     >   chi(nx)

C  We have to be careful not to lose precission in defining DX for large X.
      dx= gridx(jstop)-gridx(jstop-1)
      n = nint(dx/gridx(1))
      dx = gridx(1) * float(n)
      
      h2= dx*dx
      h2d= h2/12.
      xl= gridx(jstop)
      xlm1= gridx(jstop-1)
      wnn=sqrt(abs(ecmn))

      f3= ucentr(jstop)+cntfug(jstop)-ecmn
      f2= ucentr(jstop-1)+cntfug(jstop-1)-ecmn
      t3= (1. -h2d*f3)*s3
      t2= (1. -h2d*f2)*s2

      chi(jstop)= s3
      chi(jstop-1)= s2

      istart=njdouble-1
      do while (istart.gt.1.and.jdouble(istart).gt.jstop)
         istart=istart-1
      end do
      istop=istart
      istart=istart+1
      do while (istop.gt.1.and.jdouble(istop).gt.jstart)
         istop=istop-1
      end do
      istop=istop+1
C  JDOUBLE(ISTOP) points to the first doubling of DX that happens after JSTOP
C  JDOUBLE(ISTART) points to the last doubling of DX that happens before JSTART
      do i=istart,istop,-1
         j1=min(jstop,jdouble(i))-2
         j2=max(jstart,jdouble(i-1))
         do j=j1,j2,-1
            f1= ucentr(j)+cntfug(j)-ecmn
            t1= 2.0*t2 +h2*f2*s2 -t3
            s1= t1/(1.0d0- h2d*f1)
            t3=t2
            t2=t1
            f3=f2
            f2=f1
            s3=s2
            s2=s1
            chi(j)= s1
!            if (j.le.5) print*,'j,f1,t1,s1:',j,f1,t1,s1
         end do
         if (j2.eq.jstart) return
         j=j2-1
         dx= dx/2.0
         h2= dx*dx
         h2d= h2/12.
         f1= ucentr(j)+cntfug(j)-ecmn
         s1= s2*(36. +33.*h2*f2) +s3*(-12.+5.*h2*f3)
         s1= s1/(24. +2.*h2*f1)
         t2= s2*(1. -h2d*f2)
         t1= s1*(1. -h2d*f1)
         t3=t2
         t2=t1
         f3=f2
         f2=f1
         s3=s2
         s2=s1
         chi(j)= s1
      end do
      
      return
      end

C  This routine matches the input REG solution to the apropriate
C  boundary conditions. It matches at two points of REG to the coulomb
C  functions F and H = G+iF. The two points are chosen where REG is not
C  too small.
      subroutine matchcc(l,eta,ecmn,jstart,jstop,gridx,jm,reg,
     >   phase,sigc,j1,j2)
      dimension reg(jstop),gridx(jstop)
      real*8 dx, deta, xlmin, fc(0:1), gc(0:1), fcp(0:1), gcp(0:1)
      complex phase,sigc
      complex*16 xx,eta1,zlmin,cfc(9),cgc(9),cfcp(9),cgcp(9),sig(9),ci,
     >   f1,f2,f3,h1,h2,h3
      data pi,ci/3.1415927,(0d0,1d0)/
      if (jstart.eq.jstop) then
         phase = (1.0,0.0)
         sigc = (1.0,0.0)
         return
      endif
c$$$      if (reg(jstart).eq.0.0.and.reg(jstop).eq.0.0.and.reg(jm).eq.0.0)
c$$$     >   then
c$$$         phase = (1.0,0.0)
c$$$         sigc = (1.0,0.0)
c$$$         return
c$$$      endif
      mode1=12
      kfn=0
      deta = dble(eta)
      xlmin = dble(0)
      lrange = l
      if (ecmn.gt.0.0) then
         eta1=cmplx(dble(eta),0d0)
      else
         eta1=cmplx(0d0,-dble(eta))
      end if 
      zlmin=cmplx(dble(l),0d0)
      nl=3
      wnn=sqrt(abs(ecmn))
      if (wnn.lt.1e-15) wnn=1.0
      dx=5.0*pi/4.0/wnn

c$$$      j=min(jm+40,jstop)
c$$$C  Find the smallest value of REG closest to REG(JM)
c$$$      do while (abs(reg(j)).gt.abs(reg(j-1)).and.j-jm.gt.20)
c$$$         j=j-1
c$$$      end do
c$$$
c$$$C  Find the previous largest value of REG
c$$$      do while (abs(reg(j)).lt.abs(reg(j-1)).and.jm-j.lt.40)
c$$$         j=j-1
c$$$      end do
      j = jm
      j1 = jm
      j2 = jm
      regmax = reg(jm)
      regmin = abs(reg(jm))
      do while (j.lt.min(jm+40,jstop))
         j = j + 1
         if (reg(j).gt.regmax) then
            regmax = reg(j)
            j1 = j
         endif
         if (abs(reg(j)).lt.regmin) then
            regmin = abs(reg(j))
            j2 = j
         endif
      enddo
      if (j1.eq.j2) j1 = j2 - 2
c$$$      print*,'JM, J1, J2:', jm,j1,j2
C  First matching point
      if (j1.lt.jstart) j1=min(jstart+1,jstop)
      ifai=0
 10   xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j1))
      dx = real(xx)
      if (abs(ecmn).lt.1e-30) then
C  Here to evaluate the Coulomb G.F. at zero energy
         xx=dcmplx(sqrt(8.0*gridx(j1)))
         kfn=2
         zlmin=2*l+1
      end if 
      call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,kfn,ifai)
c$$$      call coul90(dx, deta, xlmin, lrange, fc, gc, fcp, gcp, kfn, ifai)
c$$$      print*,'CFC,fc,CGC,GC:',cfc(1),fc(l),cgc(1),gc(l)
c$$$      print*,'|CGC-(GC+iFC)|:',abs(cgc(1)-cmplx(gc(l),fc(l)))
      if (ifai.ne.0.and.j1.ne.jstop) then
         ifai=0
         j1=jstop
         j = j1
         goto 10
      else if (ifai.ne.0) then
         print*,'IFAIL, ECMN, J1 in routine MATCH:',ifai,ecmn,j1
         jstart=jstop
         jstop =jstop -1
         phase = (1.0,0.0)
         do j = 1, jstop
            reg(j) = 0.0
         enddo 
         return
      end if 
      if (abs(ecmn).gt.1e-30) then
         f1=cfc(1)
         h1=cgc(1)
c$$$         print*,j1,cgc(1) * cfcp(1) - cgcp(1) * cfc(1)
      else if (eta.ne.0.0) then
         f1=cfc(1)*sqrt(gridx(j1))*pi
         h1=cgc(1)*sqrt(gridx(j1))*pi
      else
         f1=appf1(l,0.0,gridx(j1),acc)  
         h1=cmplx(gridx(j1)**(-l),real(f1))
      end if 

c$$$C  Choose the other matching point
c$$$      do while (gridx(j)+dx.gt.gridx(j1).and.j1-j.lt.20)
c$$$         j=j-1
c$$$      end do
c$$$      j2=j
      if (j2.lt.jstart.or.j1.eq.j2) j2=j1-1
 20   xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j2))
      if (abs(ecmn).lt.1e-30) then
C  Here to evaluate the Coulomb G.F. at zero energy
         xx=dcmplx(sqrt(8.0*gridx(j2)))
         kfn=2
         zlmin=2*l+1
      end if 
      call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,kfn,ifai)
      if (ifai.ne.0) then
         if (j2.ne.j1-1) then
            ifai=0
            j2 = j1 - 1
            j = j2
            goto 20
         else if (j1.ne.jstop) then
            ifai=0
            j1 = jstop
            j = j1
            goto 10
         else
            print*,'IFAIL, ECMN, J1, J2 in routine MATCH:',
     >         ifai,ecmn,j1,j2
            jstart=jstop
            jstop =jstop -1
            phase = (1.0,0.0)
            do j = 1, jstop
               reg(j) = 0.0
            enddo 
            return
         endif 
      end if 
      if (abs(ecmn).gt.1e-30) then
         f2=cfc(1)
         h2=cgc(1)
c$$$         print*,j2,cgc(1) * cfcp(1) - cgcp(1) * cfc(1)
      else if (eta.ne.0.0) then
         f2=cfc(1)*sqrt(gridx(j2))*pi
         h2=cgc(1)*sqrt(gridx(j2))*pi
      else
         f2=appf1(l,0.0,gridx(j2),acc)
         h2=cmplx(gridx(j2)**(-l),real(f2))
      end if 
      
C  Find real and imaginary parts of the phase
      if (abs(reg(j1)*h2-reg(j2)*h1).eq.0.0) then
         print*,'Denominator is zero; j1,j2,ecmn,reg1,reg2,h1,h2:',
     >   j1,j2,ecmn,reg(j1),reg(j2),h1,h2
         phase = (1.0,0.0)
         const = 0.0
      else 
         phase=(f1*h2-f2*h1)/(reg(j1)*h2-reg(j2)*h1)
         const=abs(phase)
         phase=phase/const
      endif 
c$$$      if (abs(phase+(1.0,0.0)).lt.1e-3) then
c$$$         phase = - phase
c$$$         const = - const
c$$$      endif
      if (reg(j1).eq.0.0.or.reg(j2).eq.0.0) then
         print*,'REG(J1) or REG(J2) = 0 in MATCH',ecmn,J1,J2
         reg(:) = 0.0
         return
      end if
      do j=jstart,jstop
         reg(j)=reg(j)*const
      end do
      
      y =  (reg(j1)*f2-reg(j2)*f1) / (f2*h1-f1*h2)
c$$$      tmp1 = (reg(j1)*real(f2)-reg(j2)*real(f1))
c$$$      tmp2 = (real(f2)*real(h1)-real(f1)*real(h2))
c$$$      y =  (reg(j1)*real(f2)-reg(j2)*real(f1)) /
c$$$     >   (real(f2)*real(h1)-real(f1)*real(h2))
      phase = cmplx(real(phase),y)
      if (abs((y-aimag(phase))/(abs(y)+1e-2)).gt.1e-2) then
         print*,'Warning different imaginary parts of phase',
     >      aimag(phase),y,(reg(j1)*f2-reg(j2)*f1) / (f2*h1-f1*h2),
     >      f2,h1,f1,h2,reg(j1),reg(j2),(f2*h1-f1*h2)
      endif 
c$$$      print*,'l, ecmn, const, phase:',l,ecmn, const, phase
c$$$      if (abs(phase-(1.0,0.0)).lt.1e-5) then
c$$$         if (abs(const-1.0).gt.0.2) then
c$$$            print*,'CAUTION: CONST <> 1.0; const,phase,ecmn,l',
c$$$     >         const,phase,ecmn,l
c$$$c$$$            stop 'CONST <> 1.0'
c$$$         endif 
c$$$      endif

            
      sigc = exp(ci * sig(1))
c$$$      phase = phase * sigc
c$$$      jm = j1
C  We now test the matching procedure by calling COULCC again
      j3=(j1 + j2) / 2
      xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j3))
      if (abs(ecmn).lt.1e-30) then
C  Here to evaluate the Coulomb G.F. at zero energy
         xx=dcmplx(sqrt(8.0*gridx(j3)))
         kfn=2
         zlmin=2*l+1
      end if 
      call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,kfn,ifai)
      if (j3.ne.j1.and.j3.ne.j2.and.ifai.eq.0) then
         if (abs(ecmn).gt.1e-30) then
            f3=cfc(1)
            h3=cgc(1)
c$$$  print*,j3,cgc(1) * cfcp(1) - cgcp(1) * cfc(1)
         else if (eta.ne.0.0) then
            f3=cfc(1)*sqrt(gridx(j3))*pi
            h3=cgc(1)*sqrt(gridx(j3))*pi
         else
            f3=appf1(l,0.0,gridx(j3),acc)
            h3=cmplx(gridx(j3)**(-l),real(f3))
         end if
         test = abs(phase-(f1*h3-f3*h1)/(reg(j1)*h3-reg(j3)*h1))
         if (test.gt.1e-2) then
            print*,'Matching process has problems:',
     >         test, l, ecmn,j1,j2,j3,phase,
     >         (f1*h3-f3*h1)/(reg(j1)*h3-reg(j3)*h1),
     >        f1,h3,f3,h1,reg(j1),reg(j3)
         endif 
      endif
      return
      end

      SUBROUTINE COULOMB(ZASYM,ECMIN,L,WF,phase,jstart,jstop)
      INCLUDE 'par.f'
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasymor,lpbot,lptop,
     >   npbot(0:lamax),nptop(0:lamax)
      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   istartrp(0:ltmax),istoprp(0:ltmax),cntfug(maxr,0:lmax)
      common /double/id,jdouble(22)
      common/smallr/ formcut,regcut,expcut,fast, match
      logical fast, match
      DIMENSION WF(MAXR), utemp(maxr)
      COMPLEX phase,FAC,TTI,sigc
      common /radpot/ ucentr(maxr)
CC
CC THIS PROGRAM CALCULATES REGULAR COULOMB WAVE FUNCTIONS
CC      FOR AN ANGULAR MOMENTUM L OVER A MESH STORED IN
CC      XMESH(..,1) WHICH HAS MESH NUMBER OF POINTS
CC
C
CC     ECM IS THE ENERGY IN RYDBERGS OF THE WAVE FUNCTION
CC
c$$$      WN= SQRT(ECM)
      RMAS=1.0
c$$$      zeff = nze * (zasym + 1.0)
      zeff = - (zasym + 1.0) ! The form above fails for positron scattering
      if (zasym.eq.0.5) then ! positronium
         RMAS=0.5
         zeff = -1.0
      endif
      utemp(:) = ucentr(:)*rmas   
C  As ZASYM is the asymptotic Z of the atom, and here we need to find
C  the coulomb wave of the single electron moving in the potential of
C  the bare nucleus, the effective Z is - (ZASYM + 1.0)
      ETA = ZEFF * SQRT(RMAS/ECMIN)
      ECM = ECMIN * RMAS

c$$$      if (eta-10.0.lt.0.0) then
c$$$         ETA2=ETA*ETA
c$$$         ETA2A=2.0*ETA
c$$$         ETA6=ETA2+16.0
c$$$         SIGMAZ=-(ETA/(12.*ETA6))*(1.+(ETA2-48.)/(30.*(ETA6**2))+((ETA2
c$$$     *      -160.)*ETA2+1280.)/(105.*(ETA6**4)))-ETA+(ETA/2.)*ALOG(ETA6)
c$$$     *      +3.5*ATAN(0.25*ETA)-(ATAN(ETA)+ATAN(0.5*ETA)+ATAN(ETA/3.))
c$$$      else
c$$$         EINV1=1.0/ETA
c$$$         EINV2=EINV1*EINV1
c$$$         EINV3=EINV1*EINV2
c$$$         EINV5=EINV3*EINV2
c$$$         EINV7=EINV5*EINV2
c$$$         EINV9=EINV7*EINV2
c$$$         SIGMAZ=0.7853981634+ETA*ALOG(ETA)-ETA
c$$$     *      -(0.08333333333*EINV1+0.00277777777*EINV3
c$$$     *      +0.00079365079*EINV5+0.00059523810*EINV7
c$$$     *      +0.00084175084*EINV9)
c$$$      end if
c$$$      ETASQ=ETA*ETA
c$$$      EXSGMR=COS(SIGMAZ)
c$$$      EXSGMI=SIN(SIGMAZ)
c$$$      FL=1.0
c$$$      DO Lp=1,L
c$$$         DENOM=SQRT(1.0/(ETASQ+FL*FL))
c$$$         ETR=EXSGMR
c$$$         ETI=EXSGMI
c$$$         EXSGMR=(FL*ETR-ETA*ETI)*DENOM
c$$$         EXSGMI=(FL*ETI+ETA*ETR)*DENOM
c$$$         FL=FL+1.0
c$$$      end do
c$$$
      ldw = -1
      call regular(l,ecm,eta,utemp,cntfug(1,l),ldw,rmesh,meshr,
     >   jdouble,id,regcut,expcut,wf,jstart,jstop,phase,sigc)
      if (real(phase).lt.0.9) 
     >   print*,'Phase from Coulomb for L and E:',phase,l,ecm

c         print*,' RMAS = ',rmas,' Coulomb phase  SIGC = ',sigc
c 
c      if (rmas.eq.0.5) then
c         print*,' Coulomb phase  sigc = ',sigc
c         sigc = (1.0,0.0)       !Alisher
c      endif 
c 
      phase = phase * sigc
C  SIGC should be same as cmplx(EXSGMR,EXSGMI)
      
      if (jstart.ge.jstop) return
c$$$      phase=phase*cmplx(EXSGMR,EXSGMI)
      TTI= (0.0,1.0)
      PI= 3.1415926536
C  The following line has been changed as the /WN has been cancelled
C  with the Jacobian of the continuum integration
c$$$      FAC= SQRT(2./PI)*phase*TTI**L/WN
      FAC= SQRT(2./PI)*phase
      absfac=abs(fac)
      phase=fac/absfac
c$$$      print'(2f10.3)', phase
      DO j=jstart,jstop
         WF(j)= WF(j)*absfac
      end do
      RETURN
      END

C  Make the Green's function on GRIDX
      subroutine makegreen(both,ln,eproj,eta,ucentr,cntfug,gridx,nx,
     >   jdouble,njdouble,regcut,expcut,reg,beta,jstart,jstop,phase)

      dimension ucentr(nx),cntfug(nx),gridx(nx),jdouble(njdouble),
     >   reg(nx),beta(nx)
      complex phase,sigc
      complex*16 xx,eta1,zlmin,cfc(9),cgc(9),cfcp(9),
     >   cgcp(9),sig(9)
      character*20 stringtemp
      logical exists,both

c$$$      print*,'Entered MAKEGREEN with eta, eproj, ucentr(1):',
c$$$     >   eta,eproj,ucentr(1)
c$$$      ucentr(:)=0.0 !-ucentr(:)
      ldw = -1
      pi = acos(-1d0)
      reg(:)=0.0
      beta(:) = 0.0
      phase = (1.0,0.0)
      jstart = jfirst1(ln,eproj,gridx,nx,jdouble,njdouble,regcut)
      jstop = jlast1(eproj,gridx,nx,jdouble,njdouble,expcut)

      if (jstart+3.gt.jstop) return

      if (ln.le.3 .AND. abs(eta).LT.1e-10) then
         call ricbesselg(eproj,ln,gridx,jstop,jstart,reg) !regular
         if (.not.both) return ! if both then calculate the irregular part as well
         call ricbesselg(eproj,-ln-1,gridx,jstop,jstart,beta) !irregular
c$$$         tmp = reg(jstop)
c$$$         call regular(ln,eproj,eta,ucentr,cntfug,ldw,gridx,nx,
c$$$     >      jdouble,njdouble,small,small,reg,jstart,jstop,phase,sigc)
c$$$         print"('eproj,ln,reg(jstop),error:',i2,f10.3,1p,2e11.3)",
c$$$     >      ln,eproj,reg(jstop),(tmp-reg(jstop))/(tmp+reg(jstop))
c$$$         print"('eproj,ln,jstart,jstop,reg,beta',
c$$$     >      f10.3,3i6,1p,4e11.3)",eproj,ln,jstart,jstop,
c$$$     >      reg(jstart),reg(jstop),beta(jstart),beta(jstop)
      else ! non-zero eta or ln > 3
         call regular(ln,eproj,eta,ucentr,cntfug,ldw,gridx,nx,
     >      jdouble,njdouble,regcut,expcut,reg,jstart,jstop,phase,sigc)
         if (.not.both) return ! if both then calculate the irregular part as well
c$$$         print*,'jstart,jstop,nx,ln,eproj,reg(jstop):',
c$$$     >      jstart,jstop,nx,ln,eproj,reg(jstop)
C  Now to find the irregular solution. We take care so as not to start
C  at values that are too small. Use routine COULCC to start the G function
         xlmin=ln
         zlmin=cmplx(dble(xlmin),0d0)
         ifail=1
         nl = 1
      
         j=jstop
         if (eproj.lt.0.0) then
!         phase=phase/(0.0,1.0)
            mode1=12
            kfn=0
            eta1=eta/(0.0,1.0)
            jcount = 0
            cgc(1)= (0.0,0.0)
            do while (j.ge.jstop-1.or.jcount.lt.2.or.
     >         abs(cgc(1)).lt.expcut)
               xx=sqrt(cmplx(dble(eproj)))*dble(gridx(j))
               call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,
     >            kfn,ifail)
               if (ifail.eq.0) then !new cc.f gives errors when f < 1e-60
                  jcount = jcount + 1
                  beta(j)=imag(cfc(1)*cgc(1))/reg(j) !abs(cgc(1))
               endif
               j=j-1
            end do
            jstop = min(j+2,nx)
c$$$            do jj = j+2, jstop
c$$$               beta(jj) = beta(j+1)/exp(-sqrt(-eproj)*gridx(j+1))*
c$$$     >            exp(-sqrt(-eproj)*gridx(jj))
c$$$            enddo 

c$$$            if (eta.ne.0.0) then
c$$$               phase=cfc(1)/abs(cfc(1))*cgc(1)/abs(cgc(1))
c$$$               const = 0d0
c$$$               if (reg(j+1).ne.0d0) const = abs(cfc(1))/abs(reg(j+1))
c$$$               print'("Green phase (should be i), const and j+1:",
c$$$     >            1p,3e12.3,i5)',phase,const,j+1
c$$$               do jj = jstart, jstop
c$$$                  reg(jj) = reg(jj) * const
c$$$               enddo 
c$$$            endif 

c$$$            rph=real(phase)
c$$$            aph=aimag(phase)
c$$$            if (abs(aph).gt.1e-5) then
c$$$               print*,'Problems in making Green''s functions'
c$$$               print*,'Imaginary part of phase is too big, expect 0',aph
c$$$            end if 
         else if (eproj.eq.0.0) then
            do while (j.ge.jstop-1)
               if (eta.eq.0.0) then
                  beta(j)=gridx(j)**(-ln)
               else
C  Here to evaluate the Coulomb G.F. at zero energy
                  xx=dcmplx(sqrt(abs(eta)*8.0*gridx(j)))
                  mode1=2
                  kfn=2
                  xlmin=2*ln+1
                  zlmin=cmplx(dble(xlmin),0d0)
                  call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,
     >               mode1,kfn,ifail)
c$$$                  beta(j)=-real(cgc(1))*sqrt(gridx(j))
                  beta(j)=-real(cgc(1))*sqrt(gridx(j)*pi)
               end if
               j=j-1
            end do 
         else ! eproj positive
            jj=jstop
            do while (j.ge.jstop-1)!.or.abs(beta(jj)).lt.0.1.and.j.gt.1)
               xx=sqrt(cmplx(dble(eproj)))*dble(gridx(j))
               eta1=cmplx(dble(eta),0d0)
               mode1=2
               kfn=0
               call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,
     >            kfn,ifail)
c$$$            print*,'G.F. test',cfc(1)*cgc(2)-cfc(2)*cgc(1)
               beta(j)=real(cgc(1))
               if (ifail.ne.0) then
                  print*,'ln,eproj,cgc:',ln,eproj,beta(j)
                  stop 'COULCC failed'
               endif 
               jj=j
               j=j-1
            end do
         end if 
         j=j+2
         s2=beta(j-1)
         s3=beta(j)
         call numerovb(ln,eproj,ucentr,cntfug,gridx,nx,
     >      jdouble,njdouble,s2,s3,beta,jstart,j)
         
c$$$         rp = real(phase)
c$$$         t2=1.0/(1d-14+rp)
c$$$         t=t2*aimag(phase)
c$$$         if (abs(t).gt.1e-4) then
c$$$            DO J=JSTART,JSTOP
c$$$               BETA(J)=beta(j)*t2-t*reg(j)
c$$$            end do
c$$$         else
c$$$            print*,'PHASE,rp:',phase, rp
c$$$C  T2 may be -1.0
c$$$            if (rp.lt.1d0) then
c$$$               print*,'Expected phase to be 1:',phase
c$$$               phase = phase/rp
c$$$            endif
c$$$            DO J=JSTART,JSTOP
c$$$               BETA(J)=beta(j)  !*t2*rp
c$$$            end do
c$$$         end if 
      endif ! ln <= 3 loop
      write(stringtemp,'("GFr",1P,SP,E10.3,"_",SS,I1)') eproj,ln
      if (eproj.ne.0.0) then
C First ensure that beta(j) has rho**(-ln) behavior at small rho.
         print"('file,j,beta(j)*rho(j)**l:',a17,1p,5(i4,e10.2))",
     >      stringtemp,(j,beta(j)*(
     >      gridx(j)*sqrt(abs(eproj)))**ln,j=jstart,min(jstart+4,jstop))
         j = min(jstart + 4,jstop)
         tmp = beta(j)* (gridx(j)*sqrt(abs(eproj)))**ln
         do j = jstart, min(jstart+3,jstop)
            beta(j) = tmp / (gridx(j)*sqrt(abs(eproj)))**ln
         enddo
C These two print statements should be the same for ln, irrespective of eproj.
         print"('file,j,beta(j)*rho(j)**l:',a17,1p,5(i4,e10.2))",
     >      stringtemp,(j,beta(j)*(
     >      gridx(j)*sqrt(abs(eproj)))**ln,j=jstart,min(jstart+4,jstop))
c$$$      jstart = 1 ! This seems best for the FORM routine
      else
         print"('file:',a17)", stringtemp
      endif                     !eproj.ne.0.0
C Now ensure that the exp(rho) and exp(-rho) behaviour is not too big.
c$$$      if (eproj.lt.0.0.and.jstop.lt.nx) 
c$$$     >   call scalerb(eproj,jstart,jstop,gridx,nx,reg,beta)
c$$$      jstop = nx ! This again seems best for the FORM routine
      inquire(file='writegreen',exist=exists)
      if (exists) then
         open(72,file=stringtemp)
         if (eproj.lt.0.0) then
            eta1 = eta / (0.0,1.0)
            mode1 = 12
            kfn = 0
         else 
            eta1 = eta
            mode1 = 2
            kfn = 0
         endif 
         xlmin = ln
         zlmin=cmplx(dble(xlmin),0d0)
         write(72,"('#    R           reg         ireg       coulrho',
     >      '     Re(cfc)   	Im(cfc)	    Re(cgc)   	Im(cgc)',
     >      '      APPF      phase(cfc)*phase(cgc)',
     >      '   Im(fg)/ireg Im(fg)/reg   j')")
         do j = jstart,jstop
            xx=sqrt(cmplx(dble(eproj)))*dble(gridx(j))
            if (eproj.eq.0.0.and.eta.ne.0.0) then
               xlmin=2*ln+1
               xx=dcmplx(sqrt(-eta*8.0*gridx(j)))
               mode1=2
               kfn=2
               zlmin=cmplx(dble(xlmin),0d0)
C reg=cfc(1)*sqrt(r*pi), beta= -cgc(1)*sqrt(r*pi)
            endif
            if (eta.eq.0.0) then !coulcc crashes for zero eta
               cfc(1) = (1.0,0.0)
               cgc(1) = (1.0,0.0)
            else 
               call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,
     >         kfn,ifail)
            endif 
            if (eproj.eq.0.0.and.eta.ne.0.0) then
               cfc(1) = cfc(1)*sqrt(gridx(j)*pi)
               cgc(1) = -cgc(1)*sqrt(gridx(j)*pi)
c$$$               reg(j) = cfc(1)
c$$$               beta(j)=cgc(1)
            endif 
!            beta(j)=real(cgc(1))
            write(72,'(1p,13e12.3,i5)') gridx(j),reg(j),beta(j)
     >         ,coulrho(ln,eta,eproj,gridx(j),acc),cfc(1),cgc(1)
     >         ,appf1(ln,eproj,gridx(j),acc)
     >         ,cfc(1)/abs(cfc(1))*cgc(1)/abs(cgc(1)),
     >         imag(cfc(1)*cgc(1))/beta(j),imag(cfc(1)*cgc(1))/reg(j),j
c$$$     >         abs(cfc(1)),abs(cgc(1)),j            
         enddo
         close(72)
         print*,'Written:',stringtemp
      endif 
      RETURN
      END

      subroutine scalerb(eproj,jstart,jstop,gridx,nx,reg,beta)
      real reg(nx),beta(nx),gridx(nx)
      w = sqrt(-eproj) !eproj is negative
      rhodiff = w*gridx(nx) - 240.0
      if (rhodiff.lt.0.0) return
      if (rhodiff.lt.240.0) then
         factor = exp(rhodiff)
         do j = jstart, jstop
            reg(j) = reg(j)/factor
            beta(j) = beta(j)*factor
         enddo 
         constr = reg(jstop)*exp(-w*gridx(jstop))
         constb = beta(jstop)*exp(w*gridx(jstop))
         do j = jstop+1, nx
            reg(j) = constr*exp(w*gridx(j))
            beta(j) = constb*exp(-w*gridx(j))
         enddo
         print"('eproj,factor,reg(nx),beta(nx):',f10.3,1p,3e12.2)",
     >      eproj,factor,reg(nx),beta(nx)
      else
         print"('eproj,jstop,reg(jstop),beta(jstop):',
     >      f10.3,i7,1p,2e12.2)", eproj,jstop,reg(jstop),beta(jstop)
      endif 
      return
      end

      subroutine getnewal2(small,diffmin,diffmax,alstep,position,
     >   etot,psen2,nps,diff,it,al)
      dimension psen2(nps)
      n = 1
      small = 1e10
      do while (psen2(n).lt.position.and.n.lt.nps)
         n = n + 1
      enddo
      distance = psen2(n) - psen2(n-1)
      nf = n - 1
      if (psen2(n).lt.position) nf = n
      slowery = 
     >   (position + distance / 2.0 + 4.0*psen2(n))/5.0
      do n = 1, nps
         diff = (slowery-psen2(n)) / slowery
         if (abs(diff).lt.small) then
            small = abs(diff)
            nsmall = n
c$$$                     print*, n, diff, almin, al, almax, psen2(n) * ry
         endif
      enddo
               
      if (psen2(nps).lt.etot) then
c$$$         al = al + alstep*it
         al = al - alstep
         small = 0.0
      else 
         diff = (slowery-psen2(nsmall)) / slowery
         if (diff.lt.0.0) then
            almax = al
            diffmax = diff
            if (diffmin.gt.0.0) then
c$$$  al = (almax + almin) / 2.0
               al = (almax * diffmin - almin * diffmax) /
     >            (diffmin - diffmax)
            else
               al = al - alstep
               almin = al
            endif 
         else
            almin = al
            diffmin = diff
            if (diffmax.lt.0.0) then
               al = (almax * diffmin - almin * diffmax) /
     >            (diffmin - diffmax)
c$$$  al = (almax + almin) / 2.0
            else
               al = al + alstep
               almax = al
            endif 
         endif
      endif
      print*,'reset AL, slowery to:', al,slowery,it
      return
      end
      

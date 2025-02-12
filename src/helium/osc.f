c********** Transition moments between states of helium-like atoms ***********
c*****************************************************************************
      subroutine setdata(atom)
      character atom*20

      if(atom.eq.'helium') then
         call dataHe
      else  if(atom.eq.'LiII') then  
         call dataLiII
      else  if(atom.eq.'beryllium') then  
         call  dataBe
      else  if(atom.eq.'oxygenV') then  
         call  dataOV
      else  if(atom.eq.'magnesium') then
         call  dataMg
      else  if(atom.eq.'calcium') then
         call  dataCa
      else  if(atom.eq.'zinc') then
         call  dataZn
      else  if(atom.eq.'cadmium') then
         call  dataCd
      else  if(atom.eq.'galliumII') then
         call  dataGaII
      else  if(atom.eq.'strontium') then
         call  dataSr
      else  if(atom.eq.'barium') then
         call  dataBa
      else  if(atom.eq.'ytterbium') then
         call  dataYb
      else  if(atom.eq.'mercury') then
         call  dataHg
      end if

      return
      end
C-----------------------------------------------------------------------------
      subroutine dipmom(i_dipmom,atom,Nmax,nspmW,C,E,fl,minf,
     >     maxf,lo,enionry)
      use vmat_module, only: nodeid
      include 'par.f'
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      dimension fl(nmaxr,nspmax), numcore(0:lomax), ncount(KNM)
      dimension  minf(nspmax), maxf(nspmax), lo(nspmax)
      double precision  C(Nmax+1,nspmW,nspmW),E(KNM),dip,deriv,dippol
      character arr(0:5)*1, atom*20, parity(-1:1)*1
      data arr /"S","P","D","F","G","H"/,
     >   parity /"o"," ","e"/
      common /Nsarray/ Ns(0:lomax,0:1,komax,2)
      common /nstatearray/ nstate(0:lomax,0:1,2)
c      common /osc_strength/ oscstr(30,3,2), weight_cont(KNM)
c
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      character chan(knm)*3
      common /charchan/ chan
      real br(KNM,KNM)
      common /noblgegas_switch/ i_sw_ng  
      common /MPI_info/myid, ntasks, cnode, ench
      character cnode*3, ench*11
c
      eV = 27.2116
      
      i_ng = i_sw_ng 
      
      do NNp=1,Nmax
         do NN=1,Nmax
            br(NNp,NN) = 0.0
         end do
      end do

      do ip=1,2
         do is=0,1
            do l=0,lomax
               nmm = 0
               do N=1,Nmax
                  lN = ll(N)
                  isN = ls(N)
                  lparN = lparity(N)
                  ipN = 1
                  if(lparN.ne.(-1)**lN) ipN = 2
                  if(ip.eq.ipN.and.is.eq.isN.and.l.eq.lN) then
                     nmm = nmm + 1
                     ncount(N) = nmm
                  end if
               end do
            end do
         end do
      end do
      go to 1
c     calculation of the weights for positive en. states to place them
c     on continuum scale
c$$$      E_ion = -2.0
c$$$      do ip=1,2
c$$$         do is=0,1
c$$$            do l=0,lomax
c$$$               do num1=1,nstate(l,is,ip)
c$$$                  N = Ns(l,is,num1,ip)
c$$$                  if (E(N)-E_ion.gt.0.0) then
c$$$                     Np1 = Ns(l,is,num1+1,ip)
c$$$                     weight_cont(N) = (E(Np1) + E(N) - 2.0*E_ion)/2.0
c$$$                     nfirst = num1
c$$$                     do num=nfirst+1,nstate(l,is,ip)-1
c$$$                        N = Ns(l,is,num,ip)
c$$$                        Np1 = Ns(l,is,num+1,ip)
c$$$                        Nm1 = Ns(l,is,num-1,ip)
c$$$                        weight_cont(N) = (E(Np1) - E(Nm1))/2.0
c$$$                     end do
c$$$                     num = nstate(l,is,ip)
c$$$                     N = Ns(l,is,num,ip)
c$$$                     Nm1 = Ns(l,is,num-1,ip)
c$$$                     weight_cont(N) = (E(Nm1) + E(N) - 2.0*E_ion)/2.0
c$$$                     go to 15
c$$$                  end if
c$$$               end do
c$$$ 15         end do
c$$$         end do
c$$$      end do
c
 1    do l=0,lomax
         numcore(l) = nabot(l) - l - 1
c         print*, l, numcore(l)
      end do
c
      if (nodeid.eq.1) then
         write(4,'("dipole moments and oscillator strength ",
     >"are calculated and written in files ",
     >"osc.strength and dip.moment")')
         close(197)
         open(197,file=adjustl(cnode//'osc.strength'//ench),
     >      action='WRITE')
c$$$         open(197,file='osc.strength')

         write(197,'("oscillator strength for ",A20)') atom
         if(i_dipmom.eq.1) then
            write(197,'("E(np) > E(n): This is absorption oscillator ",
     >         "strengths")')
         end if
         if(i_dipmom.eq.2) then
            write(197,'("lp > l: This is sometimes absorption and ",
     >"sometimes emission oscillator strengths")')
         end if
         if(atom.eq.'helium') then
            write(197,'(14X,"   & FCM(l)& FCM(v)&",
     >        "Perimetric & Hylleraas & wavelen(nm) // /hline")')
         else if(atom.eq.'beryllium') then
            write(197,'(14X,"& CCC-l & CCC-v &",
     >       " RM-l  & RM-v  & wavelen(nm)  // /hline")')
         else if(atom.eq.'magnesium') then
            write(197,'(14X,"& CCC-l & CCC-v &",
     >      "MCHF-l & MCHF-v & wavelen(nm) // /hline")')
         else
         write(197,'(14X,"& length     & velocity &",
     >         "    wavelen(nm) // /hline")')
         end if
c      
         if(i_dipmom.eq.3) then
c$$$            open(232,file='dip.mom')
c$$$            open(233,file='moments_S(k)')
            open(232,file=adjustl(cnode//'dip.mom'//ench),
     >         action='WRITE')
            open(233,file=adjustl(cnode//'moments_S(k)'//ench),
     >         action='WRITE')
               write(233,'("Moments S(k), L(-1), see Eqs. 2 and 3 in ",
     >         "J.S.Bigs and Y-k.Kim Phys.Rev.A 3(971)1342")')
         end if
      endif

c     
      do ispin=0,1
         do NN=1,Nmax
            if((2.0*E(NN) - enionry).ge.0.0.and.i_dipmom.ne.3) go to 20
            l = ll(NN)
            rl = l
            is = ls(NN)
            if(ispin.ne.is) go to 20
            lpar = lparity(NN)
c     
            if(i_dipmom.eq.3) then
               dipsum = 0.0
               dipsum_v = 0.0
               dipsum_pol = 0.0
               pos_sum = 0.0
               pos_sum_pol = 0.0
               fl_sum = 0.0
               fv_sum = 0.0
               fl_sum_pol = 0.0
               sm1 = 0.0
               sp1 = 0.0
               rlm1 = 0.0
               sm1_m = 0.0
               sp1_m = 0.0
               rlm1_m = 0.0
               sm1_v = 0.0
               sp1_v = 0.0
               rlm1_v = 0.0

               if (nodeid.eq.1)
     >write(232,'("dipole polarazability of ",A3,
     >            " state N =",I3)') chan(NN),NN
               if (nodeid.eq.1)
     >write(232,'("E(Np),E(N) - eV; lp, l, modif.osc.str.:",
     >            "current state Np,summed; ",
     >            "standard osc.str.: current state Np,summed")')
               if (nodeid.eq.1)
     >write(233,'(6X,"Energy, Moments S(k), k=-2,-1,0,1,",
     >            " L(-1)  for  ",A3,
     >            " state N =",I3)') chan(NN),NN
            end if
c     
            i_br = 0
            sum_branch = 0.0
            
            do NNp=1,Nmax
               if((2.0*E(NNp)-enionry).ge.0.0.and.i_dipmom.ne.3) goto 10
               lp = ll(NNp)
               rlp = lp
               isp = ls(NNp)
               lparp = lparity(NNp)
               if(i_dipmom.eq.1.and.E(NNp).lt.E(NN)) go to 10
               if(i_dipmom.eq.2.and.
     >            (lp.lt.l.or.(l.eq.lp.and.NN.gt.NNp))) go to 10
c     if(i_dipmom.eq.3.and.lp.lt.l) go to 10
c     request lp >= l to avoid calculation of the same thing twice
               if(is.eq.isp.and.lpar.ne.lparp.and.
     >            (abs(l-1).le.lp.and.lp.le.l+1)) then               
c     >         lp.ge.l.and.
c     >         E(NNp).ge.E(NN).and.

                  if(i_ng .eq. 0) then
                     call osc(Nmax,nspmW,NNp,NN,lp,l,C,fl,minf,maxf,
     >                    lo,dip,deriv,dippol)
                  elseif(i_ng .eq. 1) then
                     call osc_ng(Nmax,NNp,NN,lp,l,fl,minf,maxf,lo,
     >                    dip,deriv,dippol)
                  else
                     print*, 'osc.f: Wrong value of i_ng=',i_ng
                     stop
                  endif

                  dip = dip/sqrt(2*l+1.0)
                  deriv = deriv/sqrt(2*l+1.0)
                  dippol = dippol/sqrt(2*l+1.0)
c
                  ae = abs(E(NNp)-E(NN)) + 1e-10 !avoid divs by zero
                  flen=2.0*ae*dip*dip/3.0
                  fvel=2.0*deriv*deriv/ae/3.0
                  flen_pol=2.0*ae*dippol*dippol/3.0
c
                  if(i_dipmom.eq.3) then
                     en_dif = ae ! Note, energy in au, not in R as in many other papers
                     tmp = flen/(en_dif*en_dif)
                     tmp_v = fvel/(en_dif*en_dif)
                     tmp_pol = flen_pol/(en_dif*en_dif)

                     dipsum = dipsum + tmp
                     dipsum_v = dipsum_v + tmp_v
                     dipsum_pol = dipsum_pol + tmp_pol
                     if((2.0*E(NNp) - enionry).lt.0.0) then
                        pos_sum = pos_sum + tmp
                        pos_sum_pol = pos_sum_pol + tmp_pol
                     end if
                     fl_sum = fl_sum + flen
                     fv_sum = fv_sum + fvel
                     fl_sum_pol = fl_sum_pol + flen_pol
                     sm1 = sm1 + flen / (en_dif * 2.0)
                     sp1 = sp1 + flen * en_dif * 2.0
                     rlm1 = rlm1 + flen *
     >                  log(2.0*abs(en_dif))/(2.0*en_dif)
                     sm1_m = sm1_m + flen_pol / (en_dif * 2.0)
                     sp1_m = sp1_m + flen_pol * en_dif * 2.0
                     rlm1_m = rlm1_m + flen_pol *
     >                  log(2.0*en_dif)/(2.0*en_dif)
                     sm1_v = sm1_v + fvel / (en_dif * 2.0)
                     sp1_v = sp1_v + fvel * en_dif * 2.0
                     rlm1_v = rlm1_v + fvel *
     >                  log(2.0*abs(en_dif))/(2.0*en_dif)

                     if (nodeid.eq.1)
     >write(232,'(A3,2E13.4,2I3,1P,4E13.4)')
     >                  chan(NNp),
     >                  (E(NNp) - enionry/2.)*eV,
     >                  (E(NN) - enionry/2.)*eV,
     >                  lp,l, tmp_pol, dipsum_pol, tmp, dipsum

                     if (nodeid.eq.1)
     >write(233,'(A3,1P,6E12.3)')
     >                  chan(NNp),
     >                  (E(NNp) - enionry/2.)*eV,
     >                  dipsum/4., sm1, fl_sum, sp1, rlm1
                  end if
c     

c     Get emmision osc.str. to calculate branching ratious
                  if(E(NNp).lt.E(NN)) then
                     br(NNp,NN) = flen
                     br(NN,NNp) = flen_pol
                  else
                     br(NNp,NN) = flen * (2*l+1.0)/(2*lp+1.0)
                     br(NN,NNp) = flen_pol * (2*l+1.0)/(2*lp+1.0)
                  end if                  
         
                  if (nodeid.eq.1)
     >write(4,'(A3,"-",A3,"   f-len =",
     >               F8.4,"(",F8.4,")","   f-vel =",F8.4,
     >               "   dipmom =",F8.4,"  deriv =",F8.4)')
     >               chan(NNp),chan(NN),
     >               flen,flen_pol,fvel,dip, deriv
               
                  if(atom.eq.'helium') then               
                     call exp_osc(atom,ncount(NN),ncount(NNp),l,lp,is,
     >                  lpar,exp_osc_l,exp_osc_v)
                     if (nodeid.eq.1)
     >               write(197,'(" $ ",A3," - ",A3," $ & ",F5.3," & ",
     >                 F5.3," & " ,F5.3," & ", F12.3," & // ")')
     >                 chan(NNp),chan(NN),flen,fvel,exp_osc_l,45.58/ae
                  else if(atom.eq.'beryllium'.or.
     >                  atom.eq.'magnesium') then
                     call exp_osc(atom,ncount(NN),ncount(NNp),l,lp,is,
     >                  lpar,exp_osc_l,exp_osc_v)
                     if (nodeid.eq.1)
     >               write(197,'("$ ",A3," - ",A3," $ & ",F5.3,
     >                  " (",F5.3,")"," & ",F5.3," & " ,F5.3,
     >                  " & " ,F5.3," & ", F12.3," &  // ")')
     >                  chan(NNp),chan(NN),flen,flen_pol,
     >                  fvel,exp_osc_l,exp_osc_v,45.58/ae
                  else
                     if (nodeid.eq.1)
     >               write(197,'("$ ",A3," - ",A3," $ & ",F5.3,
     >             " (",F5.3,")"," & " ,F5.3," & ", F12.3," &   // ")')
     >                  chan(NNp),chan(NN),flen,flen_pol,fvel,45.58/ae
                  end if
c     dipci=sqrt( cefs(num,nump,l,is) /
c     >               (2.0*abs(enci(nump,lp,is)-enci(num,l,is))) )
               end if
 10         end do
            if(i_dipmom.eq.3) then
               if (nodeid.eq.1)
     >write(232,'(5X,A3,": full polarizability: modified:",
     >            1P,E12.4,", stand.:",E12.4,
     >            ", from discreet spectrum: modified:",E12.4,
     >            ", stand.:",E12.4)')
     >            chan(NN),dipsum_pol,dipsum,pos_sum_pol, pos_sum
               if (nodeid.eq.1)
     >write(232,'("Summed osc.strength: length:",F10.4,
     >            ", length modif.:",F10.4,", velocity :",F10.4)')
     >            fl_sum, fl_sum_pol, fv_sum

               if (nodeid.eq.1)
     >write(233,'("Standard length, Summed values for ",
     >            A3,2X,1P,5E12.3)')
     >            chan(NN), dipsum/4., sm1, fl_sum, sp1, rlm1               
               if (nodeid.eq.1)
     >write(233,'("Modified length, Summed values for ",
     >            A3,2X,1P,5E12.3)')
     >            chan(NN), dipsum_pol/4., sm1_m, fl_sum_pol,
     >            sp1_m, rlm1_m
               if (nodeid.eq.1)
     >write(233,'("Stand. velocity, Summed values for ",
     >            A3,2X,1P,5E12.3)')
     >            chan(NN), dipsum_v/4., sm1_v, fv_sum, sp1_v, rlm1_v

               if (nodeid.eq.1)
     >write(232,*)
               if (nodeid.eq.1)
     >write(233,*)
               
            end if
 20      end do
      end do
      close(197)

      if (nodeid.eq.1)
     >write(4,'("Branching ratios and emmision osc.str. f_l: ",
     >   "ordinary and modified transition operator")')
      do NN=1,Nmax
         if((2.0*E(NN) - enionry).lt.0d0) then
            sum = 0.0
            sumpol = 0.0
            do NNp=1,Nmax
               if(E(NNp).lt.E(NN)) then
                  sum = sum + br(NNp,NN) * (E(NN)-E(NNp))**2
                  sumpol = sumpol + br(NN,NNp) * (E(NN)-E(NNp))**2
               end if
            end do
            if(sum.ne.0.0.or.sumpol.ne.0.0) then
               if (nodeid.eq.1)
     >write(4,*) chan(NN)
               do NNp=1,Nmax
                  if(E(NNp).lt.E(NN).and.nodeid.eq.1) then
                     if(br(NNp,NN).ne.0.0) 
     >write(4,'(3X,A3,3X,2F10.5,5X,2F10.5)')chan(NNp),
     >                  br(NNp,NN) * (E(NN)-E(NNp))**2/sum, br(NNp,NN),
     >                  br(NN,NNp) * (E(NN)-E(NNp))**2/sumpol,br(NN,NNp)
                  end if
               end do
            end if
         end if
      end do

      close(232)
      close(233)

      return
      end
c--------------------------------------------------------------------
c     < N || 2*z_1 || Np >
      subroutine   osc(Nmax,nspmW,N,Np,l,lp,C,fl,minf,maxf,lo,
     >     dip,deriv,dippol)
      include 'par.f'
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), npdrop(KNM)
      dimension fl(nmaxr,nspmax), lo(nspmax)
      dimension  minf(nspmax), maxf(nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      double precision  C(Nmax+1,nspmW,nspmW)
      double precision  ortint,trm,trm1,r1elk,f1deriv,
     >   dip,deriv,sum,coef,tmp,ang,polsum,trmpol,dippol,dintegral
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer na, nam
      common /di_el_core_polarization/ gamma, r0, pol(nmaxr)
      common /meshrr/ nr,gridr(nmaxr,3)
      double precision ttt
c
      dip = 0.0D0
      deriv = 0.0D0
      dippol = 0.0D0

      if(ls(N) .ne. ls(Np)) return

      trm = 0.0D0
      trm1 = 0.0D0
      trmpol = 0.0D0
      rl = l
      rlp = lp
      ttt = dsqrt((2.0*l+1d0)*(2.0*lp+1d0))
      do jn1=1,nam(N)
         n1 = na(jn1,N)
         l1 = lo(n1)
         rl1 = l1
         do jn1p=1,nam(Np)
            n1p = na(jn1p,Np)
            l1p = lo(n1p)
            rl1p = l1p
            if(abs(l1-1).le.l1p.or.l1p.le.l1+1) then
               sum=0d0
               do jn2=1,nam(N)
                  n2 = na(jn2,N)
                  l2 = lo(n2)
                  rl2 = l2
                  do jn2p=1,nam(Np)
                     n2p = na(jn2p,Np)
                     l2p = lo(n2p)
                     rl2p = l2p
                     if(l2.eq.l2p) then
                     tmp = ortint(n2,n2p)*C(N,jn1,jn2)*C(Np,jn1p,jn2p)
                     if(tmp.ne.0d0) then
                        ang = ttt*dsqrt(2.0*l1p+1d0)
     >                     *dble(CGC0(rl1p,1.0,rl1))
     >                     *dble(COF6J(rl2,rl1p,rlp,1.,rl,rl1))
     >                     *dble((-1)**(l1+l2+lp+1))
                        sum = sum + ang * tmp
                     end if
                     end if
                  end do
               end do
               if(sum.ne.0d0) then
                  coef = - l1p
                  if(l1p.eq.l1+1) coef = l1p+1
                  dintegral = sum*r1elk(1,fl(1,n1),fl(1,n1p),
     >               minf(n1),minf(n1p),maxf(n1),maxf(n1p))
                  trm = trm + dintegral
                  trm1 = trm1 + sum*( f1deriv(fl(1,n1),fl(1,n1p),
     >               minf(n1),minf(n1p),maxf(n1),maxf(n1p)) 
     >               + (coef-1d0)*r1elk(-1,fl(1,n1),fl(1,n1p),
     >               minf(n1),minf(n1p),maxf(n1),maxf(n1p)) )

                  polsum = 0d0
                  if(gamma.ne.0.0) then
                     mini = max(minf(n1),minf(n1p))
                     maxi = min(maxf(n1),maxf(n1p))
                     do i=mini,maxi
                        polsum = polsum + gridr(i,3)*
     >                     pol(i)*fl(i,n1)*fl(i,n1p)
                     end do
                     polsum = polsum*gamma
                  end if
                  trmpol = trmpol + dintegral - sum*polsum
               end if
            end if
         end do
      end do
      dip = trm*2.0D0
      deriv = trm1*2.0D0
      dippol = trmpol*2.0D0
      return
      end
c*******************************************************************
c     pay attention to the value of   k (k>-2, better to have k>=0)
      double precision function r1elk(k,f1,f2,l1,l2,m1,m2)
      include 'par.f'
      real*8  sum,r
      common /meshrr/ nr, gridr(nmaxr,3)
      dimension  f1(nmaxr),f2(nmaxr)

      sum = 0.0

      minfun=max(l1,l2)
      maxfun=min(m1,m2)

      if(k.eq.0) then
         do i=minfun,maxfun
            sum = sum + dble( f1(i) * f2(i) * gridr(i,3) )
         end do
      else if(k.eq.-1) then
         do i=minfun,maxfun
            r = dble(gridr(i,1))         
            sum = sum + dble( f1(i) * f2(i) * gridr(i,3) ) / r
         end do
      else if(k.eq.1) then
         do i=minfun,maxfun
            r = dble(gridr(i,1))         
            sum = sum + r * dble( f1(i) * f2(i) * gridr(i,3) )
         end do
      else if(k.eq.2) then
         do i=minfun,maxfun
            r = dble(gridr(i,1))         
            sum = sum + r*r * dble( f1(i) * f2(i) * gridr(i,3) )
         end do
      else
         do i=minfun,maxfun
            r = dble(gridr(i,1))         
            sum = sum + r**(k) * dble( f1(i) * f2(i) * gridr(i,3) )
         end do
      end if
         
      r1elk = sum
      return
      end
c------------------------------------------------------------------------
      double precision function f1deriv(f1,f2,l1,l2,m1,m2)
      include 'par.f'
      real*8  sum, h1,h2,fdf2
      common /meshrr/ nr, gridr(nmaxr,3)
      dimension  f1(nmaxr),f2(nmaxr)
      sum = 0.0d0
      minfun=max(l1,l2)
      maxfun=min(m1,m2)
      if(minfun.eq.1) then
         i = minfun
         sum = sum +  dble( f1(i) *( f2(i+1) - f2(i)) * gridr(i,3) ) /
     >      (gridr(i+1,1) - gridr(i,1))
         minfun = 2
      end if
      do i=minfun,maxfun
         h2 = gridr(i+1,1) - gridr(i,1)
         h1 = gridr(i,1) - gridr(i-1,1)
         if(h2.eq.h1) then
            fdf2 =  dble(f2(i+1) - f2(i-1))/(h1+h2)
         else
            fdf2 = (dble(f2(i+1))*(h1/h2) -
     >         dble(f2(i-1))*(h2/h1))/(h1+h2)
     >         + dble(f2(i))*((h2-h1)/h1)/h2
         end if
         sum = sum +  dble( gridr(i,3) * f1(i)) * fdf2
      end do
      f1deriv = sum
      return
      end
c------------------------------------------------------------------------
      subroutine dataHe
c     This is data for Helium
c     index 'ip' is not used here: it is assumed that lparity = (-1)**l
      parameter (nmax=8,llmax=5)
      common /oscHe/ wiesef(nmax,nmax,0:llmax,0:1),
     >     cefs(nmax,nmax,0:llmax,0:1)
      common /enHe_th/ enci(nmax,0:llmax,0:1)
      common /enHe/ enHe(nmax,0:llmax,0:1,-1:1)
      do n=1,nmax
         do l=0,llmax
            do is=0,1
               do np=1,nmax
                  wiesef(n,np,l,is)=1e30
                  cefs(n,np,l,is)=1e30
               end do
               enci(n,l,is)=1e30
            end do
         end do
      end do
c     Oscillator strenghts data from Weise (1966) wiesef(n,np,l,is)
c     It is assumed that n --> l, and np --> lp=l+1
      wiesef(1,1,0,0)=0.2762
      wiesef(1,2,0,0)=0.0734
      wiesef(1,3,0,0)=0.0302
      wiesef(2,1,0,0)=0.3764
      wiesef(2,2,0,0)=0.1514
      wiesef(2,3,0,0)=0.0507
      wiesef(3,1,0,0)=0.0480
      wiesef(3,2,0,0)=0.629
      wiesef(3,3,0,0)=0.140
      wiesef(4,1,0,0)=0.00834
      wiesef(4,2,0,0)=0.103
      wiesef(4,3,0,0)=0.853
      wiesef(1,1,1,0)=0.711
      wiesef(1,2,1,0)=0.122
      wiesef(2,1,1,0)=0.0139
      wiesef(2,2,1,0)=0.647
      wiesef(3,1,1,0)=0.00858
      wiesef(3,2,1,0)=0.0240
      wiesef(1,1,0,1)=0.5391
      wiesef(1,2,0,1)=0.06446
      wiesef(1,3,0,1)=0.0231
      wiesef(2,1,0,1)=0.0693
      wiesef(2,2,0,1)=0.896
      wiesef(2,3,0,1)=0.0429
      wiesef(3,1,0,1)=0.0118
      wiesef(3,2,0,1)=0.145
      wiesef(3,3,0,1)=1.21
      wiesef(1,1,1,1)=0.609
      wiesef(1,2,1,1)=0.125
      wiesef(2,1,1,1)=0.111
      wiesef(2,2,1,1)=0.482
      wiesef(3,1,1,1)=0.0205
      wiesef(3,2,1,1)=0.200
      
c     Oscillator strenghts data from Schieff \etal  (1971) and 
c     Weiss (1967) ($3^{1,3}D-2^{1,3}P,3^{1,3}D-3^{1,3}P$)
c     cefs(n,np,l,is)
c     It is assumed that n --> l, and np --> lp=l+1
      cefs(1,1,0,0)=0.2762
      cefs(1,2,0,0)=0.073
      cefs(1,3,0,0)=0.030
      cefs(2,1,0,0)=0.3764
      cefs(2,2,0,0)=0.1514
      cefs(2,3,0,0)=0.049
      cefs(3,1,0,0)=0.1454
      cefs(3,2,0,0)=0.626
      cefs(3,3,0,0)=0.144
      cefs(4,1,0,0)=0.0258
      cefs(4,2,0,0)=0.306
      cefs(4,3,0,0)=0.85
      cefs(1,1,1,0)=0.711
      cefs(2,1,1,0)=0.023
      cefs(1,1,0,1)=0.539086
      cefs(1,2,0,1)=0.06446
      cefs(1,3,0,1)=0.02577
      cefs(2,1,0,1)=0.20854
      cefs(2,2,0,1)=0.8909
      cefs(2,3,0,1)=0.0501
      cefs(3,1,0,1)=0.0317
      cefs(3,2,0,1)=0.4357
      cefs(3,3,0,1)=1.2153
      cefs(1,1,1,1)=0.610
      cefs(2,1,1,1)=0.111
c     energies from Accad et al (1971) and
c     Davis and Chung (1982) ($3^1D,3^3D$)    enci(n,l,is)
c     FCM ($4^1D,4^3D$)    
      do is=0,1
         do l=0,llmax               
            do n=1,nmax
               enci(n,l,is)= 1e10
            end do
         end do
      end do

      enci(1,0,0)=2.903724377
      enci(2,0,0)=2.145974044
      enci(3,0,0)=2.06127198
      enci(4,0,0)=2.0335866
      enci(1,1,0)=2.1238430858
      enci(2,1,0)=2.0551463554
      enci(3,1,0)=2.0310698
      enci(1,0,1)=2.175229378
      enci(2,0,1)=2.06868906
      enci(3,0,1)=2.03651208
      enci(1,1,1)=2.1331641905
      enci(2,1,1)=2.05808108
      enci(3,1,1)=2.0323243
      enci(1,2,0)=2.055620
      enci(2,2,0)=2.0313
      enci(1,2,1)=2.055636
      enci(2,2,1)=2.0313



c     He energy data from C.Moore  enHe(n,l,is,lpar)
c     HJ69 : J-E Holmstr\"{o}m and L. Johansson Arkiv Fysik 40,133-138
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enHe(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
            
c     sS
      enHe(1,0,0,1)= 0.0
      enHe(2,0,0,1)= 166271.70
      enHe(3,0,0,1)= 184859.06
      enHe(4,0,0,1)= 190934.5
      enHe(5,0,0,1)= 193657.78
      enHe(6,0,0,1)= 195109.17
      enHe(7,0,0,1)= 195973.19
      enHe(8,0,0,1)= 196529.03
c     sP
      enHe(1,1,0,-1)= 171129.148
      enHe(2,1,0,-1)= 186203.62
      enHe(3,1,0,-1)= 191486.95
      enHe(4,1,0,-1)= 193936.75
      enHe(5,1,0,-1)= 195269.17
      enHe(6,1,0,-1)= 196073.41
      enHe(7,1,0,-1)= 196595.56
c     sD
      enHe(1,2,0,1)= 186099.22
      enHe(2,2,0,1)= 191440.71
      enHe(3,2,0,1)= 193912.54
      enHe(4,2,0,1)= 195255.02
      enHe(5,2,0,1)= 196064.31
      enHe(6,2,0,1)= 196589.73
c     sF
      enHe(1,3,0,-1)= 191447.24        
      enHe(2,3,0,-1)= 193914.31        
      enHe(3,3,0,-1)= 195256.7        
      enHe(4,3,0,-1)= 196065.4        
      enHe(5,3,0,-1)= 196590.3  
c     tS
      enHe(1,0,1,1)= 159850.318
      enHe(2,0,1,1)= 183231.08
      enHe(3,0,1,1)= 190292.46
      enHe(4,0,1,1)= 193341.33
      enHe(5,0,1,1)= 194930.46
      enHe(6,0,1,1)= 195862.63
      enHe(7,0,1,1)= 196455.79
c     tP
      enHe(1,1,1,-1)= 169081.256 
      enHe(2,1,1,-1)= 185559.01466
      enHe(3,1,1,-1)= 191211.42
      enHe(4,1,1,-1)= 193795.07
      enHe(5,1,1,-1)= 195187.21
      enHe(6,1,1,-1)= 196021.72
      enHe(7,1,1,-1)= 196561.08
c     tD
      enHe(1,2,1,1)= 186095.90
      enHe(2,2,1,1)= 191438.83
      enHe(3,2,1,1)= 193911.48
      enHe(4,2,1,1)= 195254.37
      enHe(5,2,1,1)= 196064.00
      enHe(6,2,1,1)= 196589.42
c     tF 
      enHe(1,3,1,-1)= 191446.61
      enHe(2,3,1,-1)= 193915.79
      enHe(3,3,1,-1)= 195256.82
      enHe(4,3,1,-1)= 196065.51
      enHe(5,3,1,-1)= 196590.42
c     Convert to ionization energies (He^+ 1s) in eV.
c     Ground state ionization energy is : 24.580 eV = 198305. cm^{-1}
      gs_wave_number =  198305
c      gs_wave_number =  0.0     ! to print exciation energies
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enHe(n,l,is,lpar)=(enHe(n,l,is,lpar)-gs_wave_number)/
     >                 8065.7
               end do
            end do
         end do
      end do
            


      return
      end
      subroutine dataLiII
c     This is data for LiII
c     index 'ip' is not used here: it is assumed that lparity = (-1)**l
      parameter (nmax=8,llmax=5)
      common /enLiII/ enLiII(nmax,0:llmax,0:1,-1:1)
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enLiII(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
            
c     sS
      enLiII(1,0,0,1)= 0.0
      enLiII(2,0,0,1)= 491373.71
      enLiII(3,0,0,1)= 558777.01
      enLiII(4,0,0,1)= 581595.90
      enLiII(5,0,0,1)= 591988.68
      enLiII(6,0,0,1)= 597579.66
      enLiII(7,0,0,1)= 600929.7
      enLiII(8,0,0,1)= 603091.84
c     sP
      enLiII(1,1,0,-1)= 501807.72
      enLiII(2,1,0,-1)= 561751.95
      enLiII(3,1,0,-1)= 582836.5
      enLiII(4,1,0,-1)= 592618.4
      enLiII(5,1,0,-1)= 597936.4
      enLiII(6,1,0,-1)= 601169
      enLiII(7,1,0,-1)= 603244.7
c     sD
      enLiII(1,2,0,1)= 561272.75
      enLiII(2,2,0,1)= 582630.08
      enLiII(3,2,0,1)= 592513.56
      enLiII(4,2,0,1)= 597881.65
      enLiII(5,2,0,1)= 601118.15
      enLiII(6,2,0,1)= 603218.63
c     sF
      enLiII(1,3,0,-1)= 582643.17
      enLiII(2,3,0,-1)= 592520.24
      enLiII(3,3,0,-1)= 597885.61
      enLiII(4,3,0,-1)= 601120.68
      enLiII(5,3,0,-1)= 603220.34
c     tS
      enLiII(1,0,1,1)= 476034.11
      enLiII(2,0,1,1)= 554753.58
      enLiII(3,0,1,1)= 579980.46
      enLiII(4,0,1,1)= 591183.39
      enLiII(5,0,1,1)= 597121.08
      enLiII(6,0,1,1)= 600643.03
      enLiII(7,0,1,1)= 602902.11
c     tP
      enLiII(1,1,1,-1)= 494262.57
      enLiII(2,1,1,-1)= 559500.55
      enLiII(3,1,1,-1)= 581885.11
      enLiII(4,1,1,-1)= 592133.16
      enLiII(5,1,1,-1)= 597661.86
      enLiII(6,1,1,-1)= 600979.87
      enLiII(7,1,1,-1)= 603125.99
c     tD
      enLiII(1,2,1,1)= 561242.90
      enLiII(2,2,1,1)= 582612.54
      enLiII(3,2,1,1)= 592503.22
      enLiII(4,2,1,1)= 597875.07
      enLiII(5,2,1,1)= 601113.58
      enLiII(6,2,1,1)= 603215.34
c     tF 
      enLiII(1,3,1,-1)= 582642.10
      enLiII(2,3,1,-1)= 592519.24
      enLiII(3,3,1,-1)= 597884.56
      enLiII(4,3,1,-1)= 601119.54
      enLiII(5,3,1,-1)= 603219.16
c     Convert to ionization energies (LiIII) in eV.
c     Ground state ionization energy is 610078.5260 cm^{-1}
      gs_wave_number = 610078.5260
c      gs_wave_number =  0.0     ! to print exciation energies
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enLiII(n,l,is,lpar)=
     >               (enLiII(n,l,is,lpar)-gs_wave_number)/8065.7
               end do
            end do
         end do
      end do
            


      return
      end
c**************************************************************************
      subroutine dataBe
      parameter (nmax=8,llmax=5)
      common /oscBe/ oscBe_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
     >   oscBe_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enBe/ enBe(nmax,0:llmax,0:1,-1:1)
      do n=1,nmax
         do np=1,nmax
            do l=0,llmax
               do lp=0,llmax
                  do is=0,1
                     do lpar=-1,1
                        oscBe_l(n,np,l,lp,is,lpar)=1e30
                        oscBe_v(n,np,l,lp,is,lpar)=1e30
                     end do
                  end do
               end do
            end do
         end do
      end do
               
c  Be Oscillator strenghts data from Fon etal (1992) oscBe(n,np,l,lp,is,lpar)
c     It is assumed that n --> l, and np --> lp, and n --> lpar.
c     E(n) < E(np)   --> This is absorption oscillator strengths
      oscBe_l(1,1,0,1,0,1)= 1.408 
      oscBe_v(1,1,0,1,0,1)= 1.453 
      oscBe_l(1,2,0,1,0,1)=0.019
      oscBe_v(1,2,0,1,0,1)=0.023
      oscBe_l(1,3,0,1,0,1)=0.000
      oscBe_v(1,3,0,1,0,1)=0.002
      oscBe_l(1,2,1,0,0,-1)=0.127
      oscBe_v(1,2,1,0,0,-1)=0.121
      oscBe_l(1,3,1,0,0,-1)=0.012
      oscBe_v(1,3,1,0,0,-1)=0.010
      oscBe_l(2,2,0,1,0,1)=0.986
      oscBe_v(2,2,0,1,0,1)=0.948
      oscBe_l(2,3,0,1,0,1)=0.018
      oscBe_v(2,3,0,1,0,1)=0.001
      oscBe_l(2,3,1,0,0,-1)=0.254
      oscBe_v(2,3,1,0,0,-1)=0.233
      oscBe_l(3,3,0,1,0,1)=1.302
      oscBe_v(3,3,0,1,0,1)=1.700
      oscBe_l(1,1,1,2,0,-1)=0.002
      oscBe_v(1,1,1,2,0,-1)=0.001
      oscBe_l(1,2,1,2,0,-1)=0.395
      oscBe_v(1,2,1,2,0,-1)=0.412
      oscBe_l(1,3,1,2,0,-1)=0.428
      oscBe_v(1,3,1,2,0,-1)=0.459
      oscBe_l(1,2,2,1,0,1)=0.086
      oscBe_v(1,2,2,1,0,1)=0.040
      oscBe_l(1,3,2,1,0,1)=0.021
      oscBe_v(1,3,2,1,0,1)=0.001
      oscBe_l(2,2,1,2,0,-1)=0.697
      oscBe_v(2,2,1,2,0,-1)=0.551
      oscBe_l(2,3,1,2,0,-1)=0.005
      oscBe_v(2,3,1,2,0,-1)=0.039
      oscBe_l(2,3,2,1,0,1)=0.343
      oscBe_v(2,3,2,1,0,1)=0.099
      oscBe_l(3,3,1,2,0,-1)=1.420
      oscBe_v(3,3,1,2,0,-1)=0.599

      oscBe_l(1,1,1,0,1,-1)=0.082
      oscBe_v(1,1,1,0,1,-1)=0.086
      oscBe_l(1,2,1,0,1,-1)=0.011
      oscBe_v(1,2,1,0,1,-1)=0.011
      oscBe_l(1,2,0,1,1,1)=1.152
      oscBe_v(1,2,0,1,1,1)=1.103
      oscBe_l(1,3,0,1,1,1)=0.011
      oscBe_v(1,3,0,1,1,1)=0.005
      oscBe_l(2,2,1,0,1,-1)=0.233
      oscBe_v(2,2,1,0,1,-1)=0.203
      oscBe_l(2,3,0,1,1,1)=1.753
      oscBe_v(2,3,0,1,1,1)=1.398
      oscBe_l(1,1,1,1,1,-1)=0.457
      oscBe_v(1,1,1,1,1,-1)=0.496
      oscBe_l(2,1,1,1,1,-1)=0.000
      oscBe_v(2,1,1,1,1,-1)=0.001
      oscBe_l(1,3,1,1,1,1)=0.00
      oscBe_v(1,3,1,1,1,1)=0.00
      oscBe_l(1,1,1,2,1,-1)=0.270
      oscBe_v(1,1,1,2,1,-1)=0.255
      oscBe_l(1,2,1,2,1,-1)=0.090
      oscBe_v(1,2,1,2,1,-1)=0.087
      oscBe_l(2,1,1,2,1,-1)=0.556
      oscBe_v(2,1,1,2,1,-1)=0.505
      oscBe_l(2,2,1,2,1,-1)=0.085
      oscBe_v(2,2,1,2,1,-1)=0.123
      oscBe_l(1,3,1,2,1,1)=0.082
      oscBe_v(1,3,1,2,1,1)=0.073
      oscBe_l(3,2,1,2,1,-1)=0.863
      oscBe_v(3,2,1,2,1,-1)=0.748


      do n=1,nmax
         do np=1,nmax
            do l=0,llmax
               do lp=0,llmax
                  tmp = (2.0*l + 1.0)/(2.0*lp + 1.0)
                  do is=0,1
                     do lpar=-1,1
                        if(oscBe_l(n,np,l,lp,is,lpar).lt.1e30) then
                           oscBe_l(np,n,lp,l,is,-lpar)=
     >                          oscBe_l(n,np,l,lp,is,lpar)*tmp
                           oscBe_v(np,n,lp,l,is,-lpar)=
     >                          oscBe_v(n,np,l,lp,is,lpar)*tmp
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
c     
c     Be energy data from C.Moore  enBe(n,l,is,lpar)
c     HJ69 : J-E Holmstr\"{o}m and L. Johansson Arkiv Fysik 40,133-138
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enBe(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
            
c     sS
      enBe(1,0,0,1)= 0.0
      enBe(2,0,0,1)= 54677.2
      enBe(3,0,0,1)= 65245.4
      enBe(4,0,0,1)= 69322.3
      enBe(5,0,0,1)= 71320.7
      enBe(6,0,0,1)= 72448.3
      enBe(7,0,0,1)= 73146.7
      enBe(8,0,0,1)= 73608.5
c     sP
      enBe(1,1,0,-1)= 42565.3
      enBe(2,1,0,-1)= 60187     !
      enBe(3,1,0,-1)= 67228     !
      enBe(5,1,0,-1)= 71746.09  ! HJ69
c     sD
      enBe(1,2,0,1)= 56882.43   ! HJ69
      enBe(2,2,0,1)= 64428.15
      enBe(3,2,0,1)= 68781.2
      enBe(4,2,0,1)= 71002.3
      enBe(5,2,0,1)= 72251.1
      enBe(6,2,0,1)= 73017.2
      enBe(7,2,0,1)= 73519.7
c     tS
      enBe(1,0,1,1)= 52082.07
      enBe(2,0,1,1)= 64507.7
      enBe(3,0,1,1)= 69009.3
      enBe(4,0,1,1)= 71161.9
      enBe(5,0,1,1)= 72355.4
      enBe(6,0,1,1)= 73089.1
c     tP
      enBe(1,1,1,-1)= 21980.67
      enBe(2,1,1,-1)= 58791.6
      enBe(3,1,1,-1)= 65949     !
      enBe(4,1,1,-1)= 69634.5   !
      enBe(5,1,1,-1)= 71482.9   !
c     tD
      enBe(1,2,1,1)= 62054.8
      enBe(2,2,1,1)= 67943.6
      enBe(3,2,1,1)= 70606.7
      enBe(4,2,1,1)= 72030.6
      enBe(5,2,1,1)= 72881.9
      enBe(6,2,1,1)= 73429.6
c     TP
      enBe(1,1,1,1)= 59696.22
c     no data in C.Moore: F-states     
c     sF
      enBe(1,3,0,-1)= 68241.02  ! Fon etal 1992
      
c     tF : HJ69
      enBe(1,3,1,-1)= 68241.02
      enBe(2,3,1,-1)= 70749.72
      enBe(3,3,1,-1)= 72111.55
      enBe(4,3,1,-1)= 72931.60
c     Convert to ionization energies (Be^+ 2s) in eV.
c     Ground state ionization energy is : 9.320 eV = 75192.29 cm^{-1}
      gs_wave_number =  75192.29
c      gs_wave_number =  0.0     ! to print exciation energies
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enBe(n,l,is,lpar)=(enBe(n,l,is,lpar)-gs_wave_number)/
     >                 8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end

c**************************************************************************
      subroutine dataOV
      parameter (nmax=8,llmax=5)
      common /enOV/ enOV(nmax,0:llmax,0:1,-1:1)
c     
c     OV energy data from C.Moore  enOV(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enOV(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
            
c     sS
      enOV(1,0,0,1)= 0.0
      enOV(2,0,0,1)= 287909.
      enOV(3,0,0,1)= 561278.
      enOV(4,0,0,1)= 707630.
      enOV(5,0,0,1)= 731667.
c     sP
      enOV(1,1,0,-1)= 158798.
      enOV(2,1,0,-1)= 580826.
      enOV(3,1,0,-1)= 664486.
      enOV(4,1,0,-1)= 719277.
      enOV(5,1,0,-1)= 737883.
c     sD
      enOV(1,2,0,1)= 231722.
      enOV(2,2,0,1)= 612617.
      enOV(3,2,0,1)= 697170.
      enOV(4,2,0,1)= 746280.
c     sF
      enOV(1,3,0,-1)= 712967.
      enOV(2,3,0,-1)= 749857.
c     tS
      enOV(1,0,1,1)= 547150.
      enOV(2,0,1,1)= 684124.
      enOV(3,0,1,1)= 722666.
      enOV(4,0,1,1)= 796263.
c     tP
      enOV(1,1,1,-1)= 82413.
      enOV(2,1,1,-1)= 583058.
      enOV(3,1,1,-1)= 689793.
      enOV(4,1,1,-1)= 708326.
      enOV(5,1,1,-1)= 736125.
c     tD
      enOV(1,2,1,1)= 600943.4
      enOV(2,2,1,1)= 677639.2
      enOV(3,2,1,1)= 742412.3
c     SP
      enOV(1,1,0,1)= 672695.
c     TP
      enOV(1,1,1,1)= 213929.
      enOV(2,1,1,1)= 698793.
c     SD
      enOV(1,2,0,-1)= 694646.
c     TD
      enOV(1,2,1,-1)= 704459.
c     Convert to ionization energies (OV^+ 2s) in eV.
c     Ground state ionization energy is : 113.873 eV = 918702. cm^{-1}
      gs_wave_numOVr =  918702.
c      gs_wave_numOVr =  0.0     ! to print exciation energies
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enOV(n,l,is,lpar)=(enOV(n,l,is,lpar)-gs_wave_numOVr)/
     >                 8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end
c**************************************************************************
      subroutine dataMg
      parameter (nmax=8,llmax=5)
      common /oscMg/ oscMg_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
     >   oscMg_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enMg/ enMg(nmax,0:llmax,0:1,-1:1)

      do n=1,nmax
         do np=1,nmax
            do l=0,llmax
               do lp=0,llmax
                  do is=0,1
                     do lpar=-1,1
                        oscMg_l(n,np,l,lp,is,lpar)=1e30
                        oscMg_v(n,np,l,lp,is,lpar)=1e30
                     end do
                  end do
               end do
            end do
         end do
      end do
               
c     Mg Oscillator strenghts data from C.Froese-Fischer (1975) 
c     oscMg(n,np,l,lp,is,lpar)
c     It is assumed that n --> l, and np --> lp, and n --> lpar.
c     <np lp |O|n l>, these values do not have some factors yet, 
c     and it is symmetric  <np lp |O|n l> = <n l |O|np lp>
c     S-P singlet
      oscMg_l(1,1,0,1,0,1)=1.757
      oscMg_v(1,1,0,1,0,1)=1.736
      oscMg_l(1,2,0,1,0,1)=0.117
      oscMg_v(1,2,0,1,0,1)=0.114
      oscMg_l(1,3,0,1,0,1)=0.0260
      oscMg_v(1,3,0,1,0,1)=0.0249
      oscMg_l(1,4,0,1,0,1)=0.0095
      oscMg_v(1,4,0,1,0,1)=0.0091

      oscMg_l(2,1,0,1,0,1)=0.4658
      oscMg_v(2,1,0,1,0,1)=0.4659
      oscMg_l(2,2,0,1,0,1)=1.236
      oscMg_v(2,2,0,1,0,1)=1.240
      oscMg_l(2,3,0,1,0,1)=0.0222
      oscMg_v(2,3,0,1,0,1)=0.0223
      oscMg_l(2,4,0,1,0,1)=0.00143
      oscMg_v(2,4,0,1,0,1)=0.00140

      oscMg_l(3,1,0,1,0,1)=0.0182
      oscMg_v(3,1,0,1,0,1)=0.0183
      oscMg_l(3,2,0,1,0,1)=0.900
      oscMg_v(3,2,0,1,0,1)=0.882
      oscMg_l(3,3,0,1,0,1)=1.758
      oscMg_v(3,3,0,1,0,1)=1.751
      oscMg_l(3,4,0,1,0,1)=0.0503
      oscMg_v(3,4,0,1,0,1)=0.0491

      oscMg_l(4,1,0,1,0,1)=0.00418
      oscMg_v(4,1,0,1,0,1)=0.00421
      oscMg_l(4,2,0,1,0,1)=0.0591
      oscMg_v(4,2,0,1,0,1)=0.0569
      oscMg_l(4,3,0,1,0,1)=1.300
      oscMg_v(4,3,0,1,0,1)=1.275
      oscMg_l(4,4,0,1,0,1)=2.210
      oscMg_v(4,4,0,1,0,1)=2.195

c     S-P triplet
      oscMg_l(1,1,0,1,1,1)=1.221
      oscMg_v(1,1,0,1,1,1)=1.201
      oscMg_l(1,2,0,1,1,1)=3.942
      oscMg_v(1,2,0,1,1,1)=3.972
      oscMg_l(1,3,0,1,1,1)=0.094
      oscMg_v(1,3,0,1,1,1)=0.097
      oscMg_l(1,4,0,1,1,1)=0.019
      oscMg_v(1,4,0,1,1,1)=0.020

      oscMg_l(2,1,0,1,1,1)=0.141
      oscMg_v(2,1,0,1,1,1)=0.139
      oscMg_l(2,2,0,1,1,1)=2.522
      oscMg_v(2,2,0,1,1,1)=2.535
      oscMg_l(2,3,0,1,1,1)=5.342
      oscMg_v(2,3,0,1,1,1)=5.380
      oscMg_l(2,4,0,1,1,1)=0.178
      oscMg_v(2,4,0,1,1,1)=0.181

      oscMg_l(3,1,0,1,1,1)=0.0481
      oscMg_v(3,1,0,1,1,1)=0.0473
      oscMg_l(3,2,0,1,1,1)=0.209
      oscMg_v(3,2,0,1,1,1)=0.211
      oscMg_l(3,3,0,1,1,1)=3.787
      oscMg_v(3,3,0,1,1,1)=3.797
      oscMg_l(3,4,0,1,1,1)=6.06
      oscMg_v(3,4,0,1,1,1)=6.662

c     P-D singlet
      oscMg_l(1,1,1,2,0,-1)=0.619
      oscMg_v(1,1,1,2,0,-1)=0.636
      oscMg_l(1,2,1,2,0,-1)=0.456
      oscMg_v(1,2,1,2,0,-1)=0.445
      oscMg_l(1,3,1,2,0,-1)=0.424
      oscMg_v(1,3,1,2,0,-1)=0.416
      oscMg_l(1,4,1,2,0,-1)=0.195
      oscMg_v(1,4,1,2,0,-1)=0.192

      oscMg_l(2,1,1,2,0,-1)=0.685
      oscMg_v(2,1,1,2,0,-1)=0.699
      oscMg_l(2,2,1,2,0,-1)=2.72
      oscMg_v(2,2,1,2,0,-1)=2.77
      oscMg_l(2,3,1,2,0,-1)=0.0007
      oscMg_v(2,3,1,2,0,-1)=0.0004
      oscMg_l(2,4,1,2,0,-1)=0.0244
      oscMg_v(2,4,1,2,0,-1)=0.0241

      oscMg_l(3,1,1,2,0,-1)=0.0424
      oscMg_v(3,1,1,2,0,-1)=0.0431
      oscMg_l(3,2,1,2,0,-1)=1.38
      oscMg_v(3,2,1,2,0,-1)=1.40
      oscMg_l(3,3,1,2,0,-1)=3.99
      oscMg_v(3,3,1,2,0,-1)=4.01
      oscMg_l(3,4,1,2,0,-1)=0.0251
      oscMg_v(3,4,1,2,0,-1)=0.0247

      oscMg_l(4,1,1,2,0,-1)=0.0160
      oscMg_v(4,1,1,2,0,-1)=0.0164
      oscMg_l(4,2,1,2,0,-1)=0.0817
      oscMg_v(4,2,1,2,0,-1)=0.0832
      oscMg_l(4,3,1,2,0,-1)=2.11
      oscMg_v(4,3,1,2,0,-1)=2.12
      oscMg_l(4,4,1,2,0,-1)=5.05
      oscMg_v(4,4,1,2,0,-1)=5.02

c     P-D triplet
      oscMg_l(1,1,1,2,1,-1)=5.68
      oscMg_v(1,1,1,2,1,-1)=5.71
      oscMg_l(1,2,1,2,1,-1)=1.129
      oscMg_v(1,2,1,2,1,-1)=1.140
      oscMg_l(1,3,1,2,1,-1)=0.427
      oscMg_v(1,3,1,2,1,-1)=0.432
      oscMg_l(1,4,1,2,1,-1)=0.209
      oscMg_v(1,4,1,2,1,-1)=0.212

      oscMg_l(2,1,1,2,1,-1)=0.0502
      oscMg_v(2,1,1,2,1,-1)=0.0408
      oscMg_l(2,2,1,2,1,-1)=5.69
      oscMg_v(2,2,1,2,1,-1)=5.71
      oscMg_l(2,3,1,2,1,-1)=1.241
      oscMg_v(2,3,1,2,1,-1)=1.249
      oscMg_l(2,4,1,2,1,-1)=0.488
      oscMg_v(2,4,1,2,1,-1)=0.491

      oscMg_l(3,1,1,2,1,-1)=0.121
      oscMg_v(3,1,1,2,1,-1)=0.118
      oscMg_l(3,2,1,2,1,-1)=0.532
      oscMg_v(3,2,1,2,1,-1)=0.566
      oscMg_l(3,3,1,2,1,-1)=6.09
      oscMg_v(3,3,1,2,1,-1)=6.11
      oscMg_l(3,4,1,2,1,-1)=1.395
      oscMg_v(3,4,1,2,1,-1)=1.398

      oscMg_l(4,1,1,2,1,-1)=0.0241
      oscMg_v(4,1,1,2,1,-1)=0.0236
      oscMg_l(4,2,1,2,1,-1)=0.255
      oscMg_v(4,2,1,2,1,-1)=0.253
      oscMg_l(4,3,1,2,1,-1)=1.13
      oscMg_v(4,3,1,2,1,-1)=1.19
      oscMg_l(4,4,1,2,1,-1)=6.61
      oscMg_v(4,4,1,2,1,-1)=6.60

c     here we use the symmetry of the oscMg_l().
      do n=1,nmax
         do np=1,nmax
            do is=0,1
               do lpar=-1,1
                  oscMg_l(np,n,1,0,is,-1)= oscMg_l(n,np,0,1,is,1)
                  oscMg_l(np,n,2,1,is,1)= oscMg_l(n,np,1,2,is,-1)
                  oscMg_v(np,n,1,0,is,-1)= oscMg_v(n,np,0,1,is,1)
                  oscMg_v(np,n,2,1,is,1)= oscMg_v(n,np,1,2,is,-1)
               end do
            end do
         end do
      end do

c  Mg energy data from C.Moore  enMg(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enMg(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enMg(1,0,0,1)= 0.0
      enMg(2,0,0,1)= 43503.0
      enMg(3,0,0,1)= 52556.37
      enMg(4,0,0,1)= 56187.03
      enMg(5,0,0,1)= 58009.46
c     sP
      enMg(1,1,0,-1)= 35051.36
      enMg(2,1,0,-1)= 49346.6
      enMg(3,1,0,-1)= 54699.4
c     sD
      enMg(1,2,0,1)= 46403.14
      enMg(2,2,0,1)= 53134.70
      enMg(3,2,0,1)= 56308.43
      enMg(4,2,0,1)= 58023.27
      enMg(5,2,0,1)= 59041.09
c     sF
      enMg(1,3,0,-1)= 54676.38
      enMg(2,3,0,-1)= 57204.22
      enMg(3,3,0,-1)= 58575.54
      enMg(4,3,0,-1)= 59400.77
      enMg(5,3,0,-1)= 59935.38
c     tS
      enMg(1,0,1,1)= 41197.37
      enMg(2,0,1,1)= 51872.36
      enMg(3,0,1,1)= 55891.83
      enMg(4,0,1,1)= 57853.5
      enMg(5,0,1,1)= 58962.49
      enMg(6,0,1,1)= 59648.2
c     tP
      enMg(1,1,1,-1)= 21877.31
      enMg(2,1,1,-1)= 47849.75
      enMg(3,1,1,-1)= 54252.6
      enMg(4,1,1,-1)= 57019.45
      enMg(5,1,1,-1)= 58478.4
c     tD
      enMg(1,2,1,1)= 47957.03
      enMg(2,2,1,1)= 54192.16
      enMg(3,2,1,1)= 56968.31
      enMg(4,2,1,1)= 58442.62
      enMg(5,2,1,1)= 59317.4
c     TP
      enMg(1,1,1,1)= 57839.96      
c     tF
      enMg(1,3,1,-1)= 54676.38
      enMg(2,3,1,-1)= 57204.22
      enMg(3,3,1,-1)= 58575.54
      enMg(4,3,1,-1)= 59400.77
      enMg(5,3,1,-1)= 59935.38
c     Convert to ionization energies (Mg^+ 3s) in eV.
c     Ground state ionization energy is : 7.644 eV = 61669.14 cm^{-1}
      gs_wave_number =  61669.14
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enMg(n,l,is,lpar)=(enMg(n,l,is,lpar)-gs_wave_number)/
     >                 8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end
c**************************************************************************
      subroutine dataCa
      parameter (nmax=8,llmax=5)
c      common /oscCa/ oscCa_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscCa_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enCa/ enCa(nmax,0:llmax,0:1,-1:1)

c  Ca energy data from C.Moore  enCa(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enCa(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enCa(1,0,0,1)= 0.0
      enCa(2,0,0,1)= 33317.25
      enCa(3,0,0,1)= 41786.312
      enCa(4,0,0,1)= 44276.638
      enCa(5,0,0,1)= 45887.31
      enCa(6,0,0,1)= 46835.2
c     sP
      enCa(1,1,0,-1)= 23652.324
      enCa(2,1,0,-1)= 36731.622
      enCa(3,1,0,-1)= 41678.997
      enCa(4,1,0,-1)= 43933.341
      enCa(5,1,0,-1)= 45425.283
      enCa(6,1,0,-1)= 46479.95
c     sD
      enCa(1,2,0,1)= 21849.61
      enCa(2,2,0,1)= 37298.312
      enCa(3,2,0,1)= 40719.867
      enCa(4,2,0,1)= 42919.074
      enCa(5,2,0,1)= 44989.882
      enCa(6,2,0,1)= 46309.9
c     sF
      enCa(1,3,0,-1)= 40537.860
      enCa(2,3,0,-1)= 42343.554
      enCa(3,3,0,-1)= 44804.786
      enCa(4,3,0,-1)= 46182.23
      enCa(5,3,0,-1)= 47015.137
c     tS
      enCa(1,0,1,1)= 31539.510
      enCa(2,0,1,1)= 40474.275
      enCa(3,0,1,1)= 43980.798
      enCa(4,0,1,1)= 45738.732
      enCa(5,0,1,1)= 46748.21
c     tP
      enCa(1,1,1,-1)= 15227.975
      enCa(2,1,1,-1)= 36560.
      enCa(3,1,1,-1)= 39337.
      enCa(4,1,1,-1)= 42520.012
      enCa(5,1,1,-1)= 44960.
c     tD
      enCa(1,2,1,1)= 20351.86
      enCa(2,2,1,1)= 37752.513
      enCa(3,2,1,1)= 42745.1
      enCa(4,2,1,1)= 45050.6
      enCa(5,2,1,1)= 46305.
c     tF
      enCa(1,3,1,-1)= 35816.
      enCa(2,3,1,-1)= 42170.5
      enCa(3,3,1,-1)= 44762.9
      enCa(4,3,1,-1)= 46164.82
      enCa(5,3,1,-1)= 47006.11
c     SD
      enCa(1,2,0,-1)= 35835.4
c     TD
      enCa(1,2,1,-1)= 38224.
c     TP
      enCa(1,1,1,1)= 38478.      

c     Convert to ionization energies (Ca^+ 4s) in eV.
c     Ground state ionization energy is : 6.11 eV = 49304.80 cm^{-1}
      gs_wave_number =  49304.80
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enCa(n,l,is,lpar)=(enCa(n,l,is,lpar)-gs_wave_number)/
     >                 8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end
      
c**************************************************************************
      subroutine dataZn
      parameter (nmax=8,llmax=5)
c      common /oscZn/ oscZn_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscZn_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enZn/ enZn(nmax,0:llmax,0:1,-1:1)

c  Znnergy data from C.Moore  enZnn,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enZn(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enZn(1,0,0,1)= 0.0
      enZn(2,0,0,1)= 55789.22
      enZn(3,0,0,1)= 66037.60
      enZn(4,0,0,1)= 70003.72
      enZn(5,0,0,1)= 71822.5
c     sP
      enZn(1,1,0,-1)= 46745.37
      enZn(2,1,0,-1)= 62910.0
      enZn(3,1,0,-1)= 68607.26
      enZn(4,1,0,-1)= 71219.08
      enZn(5,1,0,-1)= 72626.2
c     sD
      enZn(1,2,0,1)= 62458.51
      enZn(2,2,0,1)= 68338.48
      enZn(3,2,0,1)= 71050.45
      enZn(4,2,0,1)= 72489.13
c     sF
      enZn(1,3,0,-1)= 68834.25  
      enZn(2,3,0,-1)= 71336.15 
c     tS
      enZn(1,0,1,1)= 53672.241
      enZn(2,0,1,1)= 65432.333  
      enZn(3,0,1,1)= 69745.96 
      enZn(4,0,1,1)= 71822.5
c     tP
      enZn(1,1,1,-1)= 32696.34
      enZn(2,1,1,-1)= 61302.2
      enZn(3,1,1,-1)= 68091.4
      enZn(4,1,1,-1)= 70985.
      enZn(5,1,1,-1)= 72501
c     tD
      enZn(1,2,1,1)= 62774.
      enZn(2,2,1,1)= 68582.
      enZn(3,2,1,1)= 71213.5
      enZn(4,2,1,1)= 72627.9
c     tF
      enZn(1,3,1,-1)= 68834.4   
      enZn(2,3,1,-1)= 71335.6
      enZn(3,3,1,-1)= 72690.8
      enZn(4,3,1,-1)= 73499.5 

c     Convert to ionization energies (Zn^+ 4s) in eV.
c     Ground state ionization energy is : 9.391 eV = 75766.8 cm^{-1}
      gs_wave_number =  75766.8
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enZn(n,l,is,lpar)=(enZn(n,l,is,lpar) - 
     >                 gs_wave_number)/8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end
c**************************************************************************
      subroutine dataCd
      parameter (nmax=8,llmax=5)
c      common /oscCd/ oscCd_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscCd_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enCd/ enCd(nmax,0:llmax,0:1,-1:1)

c  Cd energy data from NIST   enCdn,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enCd(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enCd(1,0,0,1)= 0.0
      enCd(2,0,0,1)= 53310.101
c     sP
      enCd(1,1,0,-1)= 43692.384
c     sD
      enCd(1,2,0,1)= 59219.734
c     sF
!      enCd(1,3,0,-1)= 68834.25  
c     tS
      enCd(1,0,1,1)= 51483.98
c     tP
      enCd(1,1,1,-1)= 30865.0
c     tD
      enCd(1,2,1,1)= 59499.9
c     tF
!      enCd(1,3,1,-1)= 68834.4   

c     Convert to ionization energies (Cd^+ 4s) in eV.
c     Ground state ionization energy is : 8.99382 eV = 72540.07 cm^{-1}
      gs_wave_number =   72540.07
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enCd(n,l,is,lpar)=(enCd(n,l,is,lpar) - 
     >                 gs_wave_number)/8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end
c**************************************************************************
      subroutine dataGaII
      parameter (nmax=8,llmax=5)
c      common /oscGaII/ oscGaII_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscGaII_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enGaII/ enGaII(nmax,0:llmax,0:1,-1:1)


c  GaII energy data from NIST web  enGaII(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enGaII(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enGaII(1,0,0,1)= 0.0
      enGaII(2,0,0,1)= 106662.379
      enGaII(3,0,0,1)= 133741.20
      enGaII(4,0,0,1)= 135639.40
      enGaII(5,0,0,1)= 146 023.600
      enGaII(6,0,0,1)= 152 200.866
c     sP
      enGaII(1,1,0,-1)= 70701.427
      enGaII(2,1,0,-1)= 120550.431
!      enGaII(3,1,0,-1)= 68607.26

c     sD
      enGaII(1,2,0,1)= 107720.716
      enGaII(2,2,0,1)= 126187.61
      enGaII(3,2,0,1)= 139703.442
      enGaII(4,2,0,1)= 148446.729

c     sF
      enGaII(1,3,0,-1)= 137342.570 
      enGaII(2,3,0,-1)= 147493.327
      enGaII(3,3,0,-1)= 153 009.57 

c     tS
      enGaII(1,0,1,1)= 102944.595
      enGaII(2,0,1,1)= 133010.30
      enGaII(3,0,1,1)= 145494.205
      enGaII(4,0,1,1)= 151923.93

c     tP
      enGaII(1,1,1,-1)= 47814.114
      enGaII(2,1,1,-1)= 118518.461

c     tD
      enGaII(1,2,1,1)= 113842.301
      enGaII(2,2,1,1)= 137168.592
      enGaII(3,2,1,1)= 147526.09 
      enGaII(4,2,1,1)= 153068.30

c     tF
      enGaII(1,3,1,-1)= 137333.474
      enGaII(2,3,1,-1)= 147 485.466
      enGaII(3,3,1,-1)= 153001.63

c     TP
      enGaII(1,1,1,1)= 115224.47

c     Convert to ionization energies (GaII^+ 4s) in eV.
c     Ground state ionization energy is : 20.51514 eV = 165465.8 cm^{-1}
      gs_wave_number =  165465.8
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enGaII(n,l,is,lpar)=(enGaII(n,l,is,lpar) - 
     >                 gs_wave_number)/8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end
c**************************************************************************
      subroutine dataSr
      parameter (nmax=8,llmax=5)
c      common /oscSr/ oscSr_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscSr_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enSr/ enSr(nmax,0:llmax,0:1,-1:1)

c  Sr energy data from C.Moore  enSr(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enSr(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enSr(1,0,0,1)= 0.0
      enSr(2,0,0,1)= 30591.8     
      enSr(3,0,0,1)= 37160.278
      enSr(4,0,0,1)= 38444.054
      enSr(5,0,0,1)= 41052.5
      enSr(6,0,0,1)= 42596.0
c     sP
      enSr(1,1,0,-1)= 21698.482
      enSr(2,1,0,-1)= 34098.44
      enSr(3,1,0,-1)= 38906.90
      enSr(4,1,0,-1)= 41172.15
      enSr(5,1,0,-1)= 42462.36
c     sD
      enSr(1,2,0,1)= 20149.7
      enSr(2,2,0,1)= 34727.483
      
      enSr(3,2,0,1)= 36960.881
      enSr(4,2,0,1)= 39733.114
      enSr(5,2,0,1)= 41831.7
      enSr(6,2,0,1)= 43020.9
c     sF
      enSr(1,3,0,-1)= 38008.0      
      enSr(2,3,0,-1)=  39539.04
      enSr(3,3,0,-1)= 41518.91
c     tS
      enSr(1,0,1,1)= 29038.795
      enSr(2,0,1,1)= 37424.713      
      enSr(3,0,1,1)= 40761.44
      enSr(4,0,1,1)= 42451.2
c     tP
      enSr(1,1,1,-1)= 14702.987
      enSr(2,1,1,-1)= 33924.88      
      enSr(3,1,1,-1)= 39442.02
c     tD
      enSr(1,2,1,1)= 18267.98
      enSr(2,2,1,1)= 35033.14      
      enSr(3,2,1,1)= 39695.
      enSr(4,2,1,1)= 41870.
c     tF
      enSr(1,3,1,-1)= 33654.
      enSr(2,3,1,-1)= 38753.
      enSr(3,3,1,-1)= 41366.
c     SD
      enSr(1,2,0,-1)= 33826.927
c     TD
      enSr(1,2,1,-1)= 36441.2
c     TP
      enSr(1,1,1,1)= 35529.69

c     Convert to ionization energies (Sr^+ 4s) in eV.
c     Ground state ionization energy is : 5.692 eV = 45925.6 cm^{-1}
      gs_wave_number =  45925.6
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enSr(n,l,is,lpar)=(enSr(n,l,is,lpar)-gs_wave_number)/
     >                 8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end
      
c**************************************************************************

      subroutine dataBa
      parameter (nmax=8,llmax=5)
c      common /oscBa/ oscBa_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscBa_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enBa/ enBa(nmax,0:llmax,0:1,-1:1)

c  Ba energy data from C.Moore  enBa(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enBa(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enBa(1,0,0,1)= 0.0
      enBa(2,0,0,1)= 26757.4 ! this is 5d5d level which is absent in Moore tables but given in H.P.Palenius Phys.Let. 56A(1976)451-452
      enBa(3,0,0,1)= 28230.08
      enBa(4,0,0,1)= 34370.78
      enBa(5,0,0,1)= 37041.00
      enBa(6,0,0,1)= 38267.59
      enBa(7,0,0,1)= 38662.91
c     sP
      enBa(1,1,0,-1)= 18060.264
      enBa(2,1,0,-1)= 28554.257
      enBa(3,1,0,-1)= 32547.076
      enBa(4,1,0,-1)= 35892.518
      enBa(5,1,0,-1)= 37774.50
      enBa(6,1,0,-1)= 38499.852
c     sD
      enBa(1,2,0,1)= 11395.382
      enBa(2,2,0,1)= 23062.06
      enBa(3,2,0,1)= 30236.815
      enBa(4,2,0,1)= 33795.84
      enBa(5,2,0,1)= 35344.423
      enBa(6,2,0,1)= 37434.956
      enBa(7,2,0,1)= 37837.40
      enBa(8,2,0,1)= 38556.18      
c     sF
      enBa(1,3,0,-1)= 26816.293
      enBa(2,3,0,-1)= 34736.423
      enBa(3,3,0,-1)= 37282.175
      enBa(4,3,0,-1)= 37739.95
      enBa(5,3,0,-1)= 38884.00
c     sG
      enBa(1,4,0,1)= 24300. ! this is 5d5d level which is absent in Moore tables but given in H.P.Palenius Phys.Let. 56A(1976)451-452
      enBa(2,4,0,1)= 38177.10
c     tS
      enBa(1,0,1,1)= 26160.284
      enBa(2,0,1,1)= 33905.349
      enBa(3,0,1,1)= 36446.583
      enBa(4,0,1,1)= 37095.492
      enBa(5,0,1,1)= 38662.91
c     tP
      enBa(1,1,1,-1)= 13082.669
      enBa(2,1,1,-1)= 25837.48
      enBa(3,1,1,-1)= 30902.96
      enBa(4,1,1,-1)= 37143.16
c     tD
      enBa(1,2,1,1)= 9357.03
      enBa(2,2,1,1)= 30771.13
      enBa(3,2,1,1)= 33187.66
      enBa(4,2,1,1)= 35762.41
      enBa(5,2,1,1)= 36347.05
      enBa(6,2,1,1)= 37978.41
c     tF
      enBa(1,3,1,-1)= 23084.24
      enBa(2,3,1,-1)= 34619.43
      enBa(3,3,1,-1)= 36711.52
      enBa(4,3,1,-1)= 37458.34
      enBa(5,3,1,-1)= 38820.71
c     tG
      enBa(1,4,1,1)= 36430.12
      
c     SD
      enBa(1,2,0,-1)= 23074.416
      enBa(2,2,0,-1)= 37077.50
c     TD
      enBa(1,2,1,-1)= 24672.85
      enBa(2,2,1,-1)= 37172.37
c     TP
      enBa(1,1,1,1)= 23693.76
      enBa(2,1,1,1)= 35227.655
      enBa(3,1,1,1)= 38120.43
c     SP
      enBa(1,1,0,1)= 36902.55
c     SF
      enBa(1,3,0,1)= 36165.311
c     TF
      enBa(1,3,1,1)= 21335.012 ! this is 5d5d level which is absent in Moore tables but given in H.P.Palenius Phys.Let. 56A(1976)451-452
      enBa(2,3,1,1)= 37502.97
      

c     Convert to ionization energies (Ba^+ 6s) in eV.
c     Ground state ionization energy is : 5.21 eV = 42032.40 cm^{-1}
      gs_wave_number =  42032.40
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enBa(n,l,is,lpar)=(enBa(n,l,is,lpar)-gs_wave_number)/
     >                 8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end

c**************************************************************************
      subroutine dataYb
      parameter (nmax=8,llmax=5)
c      common /oscHg/ oscHg_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscHg_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enYb/ enYb(nmax,0:llmax,0:1,-1:1)

c  Yb energy data from C.Moore  enYb(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enYb(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enYb(1,0,0,1)= 0.0
      enYb(2,0,0,1)= 34350.65  
      enYb(3,0,0,1)= 41939.90  
c     sP
      enYb(1,1,0,-1)= 25068.222  
      enYb(2,1,0,-1)= 40563.97  
c     sD
      enYb(1,2,0,1)= 27677.665  
      enYb(2,2,0,1)= 40061.51  
c     sF
c      enYb(1,3,0,-1)= 77241.68
c     sG
c      enYb(1,4,0,1)= 
c     tS
      enYb(1,0,1,1)= 32694.692  
      enYb(2,0,1,1)= 41615.04  
c     tP
      enYb(1,1,1,-1)= 19500.0  ! 17288.439, 17992.007,  19710.388  
      enYb(2,1,1,-1)= 38400.0  ! 38090.71,  38174.17,  38551.93  
c     tD
      enYb(1,2,1,1)= 25000.0   ! 24489.102, 24751.948, 25270.902 
      enYb(2,2,1,1)= 39966.09  ! 39808.72, 39838.04, 39966.09   
c     tF
c      enYb(1,3,1,-1)= 77259.17 ! 77237.04, 77239.2, 77286.99
c     tG
c      enYb(1,4,1,1)= 36430.12
      
c     SD
c      enYb(1,2,0,-1)= 23074.416
c      enYb(2,2,0,-1)= 37077.50
c     TD
c      enYb(1,2,1,-1)= 24672.85
c      enYb(2,2,1,-1)= 37172.37
c     TP
      enYb(1,1,1,1)= 44400.0  ! 42436.91, 43805.42, 44760.37  

c     SP
c      enYb(1,1,0,1)= 36902.55
c     SF
c      enYb(1,3,0,1)= 36165.311
c     TF
c      enYb(1,3,1,1)= 21335.012 
c      enYb(2,3,1,1)= 37502.97
      

c     Convert to ionization energies (Yb^+ 6s) in eV.
c     Ground state ionization energy is : 6.2542 eV =  50441.0  cm^{-1}
      gs_wave_number = 50441.0 ! 52778.5
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enYb(n,l,is,lpar)=(enYb(n,l,is,lpar)-gs_wave_number)/
     >                 8065.7
               end do
            end do
         end do
      end do
            
      return
      end
c**************************************************************************
      subroutine dataHg
      parameter (nmax=8,llmax=5)
c      common /oscHg/ oscHg_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscHg_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enHg/ enHg(nmax,0:llmax,0:1,-1:1)

c  Hg energy data from C.Moore  enHg(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enHg(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enHg(1,0,0,1)= 0.0
      enHg(2,0,0,1)= 63928.243
      enHg(3,0,0,1)= 74404.590
      enHg(4,0,0,1)= 78404.387
      enHg(5,0,0,1)= 80365.653
      enHg(6,0,0,1)= 81473.4
      enHg(7,0,0,1)= 82160.8
      enHg(8,0,0,1)= 82616.2
c     sP
      enHg(1,1,0,-1)= 54068.781
      enHg(2,1,0,-1)= 71295.15
      enHg(3,1,0,-1)= 76863.264
      enHg(4,1,0,-1)= 79964.1
      enHg(5,1,0,-1)= 81153.614
      enHg(6,1,0,-1)= 81942.444
      enHg(7,1,0,-1)= 82464.05
c     sD
      enHg(1,2,0,1)= 71333.182
      enHg(2,2,0,1)= 77064.097
      enHg(3,2,0,1)= 79660.785
      enHg(4,2,0,1)= 81057.749
      enHg(5,2,0,1)= 81895.0
      enHg(6,2,0,1)= 82436.2
      enHg(7,2,0,1)= 82807.421
      enHg(8,2,0,1)= 83069.2     
c     sF
      enHg(1,3,0,-1)= 77241.68
      enHg(2,3,0,-1)= 79745.3
      enHg(3,3,0,-1)= 81106.5
      enHg(4,3,0,-1)= 81925.3
      enHg(5,3,0,-1)= 82455.3
      enHg(6,3,0,-1)= 82820.4
c     sG
c      enHg(1,4,0,1)= 
c      enHg(2,4,0,1)= 
c     tS
      enHg(1,0,1,1)= 62350.456
      enHg(2,0,1,1)= 73961.298
      enHg(3,0,1,1)= 78216.261
      enHg(4,0,1,1)= 80268.056
      enHg(5,0,1,1)= 81416.352
      enHg(6,0,1,1)= 82124.081
      enHg(7,0,1,1)= 82591.3
c     tP
      enHg(1,1,1,-1)= 41788.5  ! 37645.08, 39412.3, 44042.977
      enHg(2,1,1,-1)= 70504.43 ! 69516.66, 69661.89, 71207.51
      enHg(3,1,1,-1)= 76662.88 ! 76447.243, 76467.067, 76823.5
      enHg(4,1,1,-1)= 79520.06 ! 79375.783, 79412.745, 79613.3
      enHg(5,1,1,-1)= 80974.09 ! 80902.27, 80916.686, 81022.9
c     tD
      enHg(1,2,1,1)= 71400.58  ! 71336.164, 71396.22, 71431.311
      enHg(2,2,1,1)= 77113.35  ! 77084.632, 77107.917, 77129.535
      enHg(3,2,1,1)= 79693.74  ! 79678.708, 79690.3, 79702.634
      enHg(4,2,1,1)= 81079.86  ! 81071.027, 81077.8, 81085.126
      enHg(5,2,1,1)= 81910.16  ! 81904.5, 81908.7, 81913.632
      enHg(6,2,1,1)= 82446.86  ! 82443.0, 82445.9, 82449.2
c     tF
      enHg(1,3,1,-1)= 77259.17 ! 77237.04, 77239.2, 77286.99
      enHg(2,3,1,-1)= 79746.1  ! 79743.7, 79745.0, 79748.3
      enHg(3,3,1,-1)= 81105.5  ! 81103.9, 81104.6, 81107.1
      enHg(4,3,1,-1)= 81925.6  ! 81923.5, 81924.3, 81927.8
      enHg(5,3,1,-1)= 82455.8  ! 82454.9, 82455.2, 82456.7
c     tG
c      enHg(1,4,1,1)= 36430.12
      
c     SD
c      enHg(1,2,0,-1)= 23074.416
c      enHg(2,2,0,-1)= 37077.50
c     TD
c      enHg(1,2,1,-1)= 24672.85
c      enHg(2,2,1,-1)= 37172.37
c     TP
c      enHg(1,1,1,1)= 23693.76
c      enHg(2,1,1,1)= 35227.655
c      enHg(3,1,1,1)= 38120.43
c     SP
c      enHg(1,1,0,1)= 36902.55
c     SF
c      enHg(1,3,0,1)= 36165.311
c     TF
c      enHg(1,3,1,1)= 21335.012 
c      enHg(2,3,1,1)= 37502.97
      

c     Convert to ionization energies (Hg^+ 6s) in eV.
c     Ground state ionization energy is : 10.43 eV = 84184.1 cm^{-1}
      gs_wave_number = 84184.1 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enHg(n,l,is,lpar)=(enHg(n,l,is,lpar)-gs_wave_number)/
     >                 8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end
c**************************************************************************
      subroutine dataNe
      parameter (nmax=8,llmax=5)
c      common /oscNe/ oscNe_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscNe_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enNe/ enNe(nmax,0:llmax,0:1,-1:1)

c  Ne energy data from C.Moore  enNe(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enNe(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enNe(1,0,0,1)= 0.0           ! 2p2p
      enNe(2,0,0,1)= 152970.7328   ! 2p3p
      enNe(3,0,0,1)= 164285.8872   ! 2p4p 
      enNe(4,0,0,1)= 168586.8304   ! 2p5p
c     sP
      enNe(1,1,0,-1)= 135888.7173  ! 2p3s 
      enNe(2,1,0,-1)= 159534.6196  ! 2p4s
      enNe(3,1,0,-1)= 161636.6175  ! 2p3d
      enNe(4,1,0,-1)= 166975.3424  ! 2p5s

c$$$c     sD
      enNe(1,2,0,1)= 150315.8612   ! 2p3p
      enNe(2,2,0,1)= 163708.6029   ! 2p4p

c$$$c     sF
      enNe(1,3,0,-1)= 161 701.4486 ! 2p3d

c$$$c     sG
c$$$c      enNe(1,4,0,1)= 
c$$$c      enNe(2,4,0,1)= 

c$$$c     tS
      enNe(1,0,1,1)= 150917.4307    ! 2p3p
      enNe(2,0,1,1)= 163401.3061    ! 2p4p
      enNe(3,0,1,1)= 167867.1941    ! 2p5p

c$$$c     tP
      enNe(1,1,1,-1)= 134459.2871  ! 2p3s
      enNe(2,1,1,-1)= 158795.9924  ! 2p4s
      enNe(3,1,1,-1)= 161524.1739  ! 2p3d
      enNe(4,1,1,-1)= 165912.7861  ! 2p5s

c$$$c     tD
      enNe(1,2,1,1)= 149 824.2215 ! 2p3p
      enNe(2,2,1,1)= 162899.1169  ! 2p4p

c$$$c     tF
      enNe(1,3,1,-1)= 161 592.1200 ! 2p3d

c     tG
c      enNe(1,4,1,1)= 36430.12
      
c     SD
      enNe(1,2,0,-1)= 161699.6613   ! 2p3d


c     TD
      enNe(1,2,1,-1)= 161607.2609   ! 2p3d

c     TP
      enNe(1,1,1,1)= 148257.7898    !  2p3p
      enNe(2,1,1,1)= 162517.8755    !  2p3p


c     SP
      enNe(1,1,0,1)= 150121.5922    !  2p3p
      enNe(2,1,0,1)= 163012.6247    !  2p4p

c     SF
c      enNe(1,3,0,1)= 36165.311
c     TF
c      enNe(1,3,1,1)= 21335.012 
c      enNe(2,3,1,1)= 37502.97
      

c     Convert to ionization energies (Ne^+ 62p^5) in eV.
c     Ground state ionization energy is : 21.5645 eV = 173929.75  cm^{-1}
      gs_wave_number = 173929.75 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enNe(n,l,is,lpar)=(enNe(n,l,is,lpar)-gs_wave_number)/
     >                   8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end
c**************************************************************************
      subroutine dataAr
      parameter (nmax=8,llmax=5)
c      common /oscAr/ oscAr_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
c     >   oscAr_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /enAr/ enAr(nmax,0:llmax,0:1,-1:1)

c  Ar energy data from C.Moore  enAr(n,l,is,lpar)
c     excitation energies (cm^{-1}) 
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax               
               do n=1,nmax
                  enAr(n,l,is,lpar)= 1e10
               end do
            end do
         end do
      end do
c	sS
      enAr(1,0,0,1)= 0.0
c$$$      enAr(2,0,0,1)= 63928.243
c$$$      enAr(3,0,0,1)= 74404.590
c$$$      enAr(4,0,0,1)= 78404.387
c$$$      enAr(5,0,0,1)= 80365.653
c$$$      enAr(6,0,0,1)= 81473.4
c$$$      enAr(7,0,0,1)= 82160.8
c$$$      enAr(8,0,0,1)= 82616.2
c     sP
      enAr(1,1,0,-1)= 95399.8276
c$$$      enAr(2,1,0,-1)= 225809.07
c$$$      enAr(3,1,0,-1)= 229950.46
c$$$      enAr(4,1,0,-1)= 79964.1
c$$$      enAr(5,1,0,-1)= 81153.614
c$$$      enAr(6,1,0,-1)= 81942.444
c$$$      enAr(7,1,0,-1)= 82464.05
c$$$c     sD
c$$$      enAr(1,2,0,1)= 71333.182
c$$$      enAr(2,2,0,1)= 77064.097
c$$$      enAr(3,2,0,1)= 79660.785
c$$$      enAr(4,2,0,1)= 81057.749
c$$$      enAr(5,2,0,1)= 81895.0
c$$$      enAr(6,2,0,1)= 82436.2
c$$$      enAr(7,2,0,1)= 82807.421
c$$$      enAr(8,2,0,1)= 83069.2     
c$$$c     sF
c$$$      enAr(1,3,0,-1)= 77241.68
c$$$      enAr(2,3,0,-1)= 79745.3
c$$$      enAr(3,3,0,-1)= 81106.5
c$$$      enAr(4,3,0,-1)= 81925.3
c$$$      enAr(5,3,0,-1)= 82455.3
c$$$      enAr(6,3,0,-1)= 82820.4
c$$$c     sG
c$$$c      enAr(1,4,0,1)= 
c$$$c      enAr(2,4,0,1)= 
c$$$c     tS
c$$$      enAr(1,0,1,1)= 62350.456
c$$$      enAr(2,0,1,1)= 73961.298
c$$$      enAr(3,0,1,1)= 78216.261
c$$$      enAr(4,0,1,1)= 80268.056
c$$$      enAr(5,0,1,1)= 81416.352
c$$$      enAr(6,0,1,1)= 82124.081
c$$$      enAr(7,0,1,1)= 82591.3
c$$$c     tP
c$$$      enAr(1,1,1,-1)= 41788.5  ! 37645.08, 39412.3, 44042.977
c$$$      enAr(2,1,1,-1)= 70504.43 ! 69516.66, 69661.89, 71207.51
c$$$      enAr(3,1,1,-1)= 76662.88 ! 76447.243, 76467.067, 76823.5
c$$$      enAr(4,1,1,-1)= 79520.06 ! 79375.783, 79412.745, 79613.3
c$$$      enAr(5,1,1,-1)= 80974.09 ! 80902.27, 80916.686, 81022.9
c$$$c     tD
c$$$      enAr(1,2,1,1)= 71400.58  ! 71336.164, 71396.22, 71431.311
c$$$      enAr(2,2,1,1)= 77113.35  ! 77084.632, 77107.917, 77129.535
c$$$      enAr(3,2,1,1)= 79693.74  ! 79678.708, 79690.3, 79702.634
c$$$      enAr(4,2,1,1)= 81079.86  ! 81071.027, 81077.8, 81085.126
c$$$      enAr(5,2,1,1)= 81910.16  ! 81904.5, 81908.7, 81913.632
c$$$      enAr(6,2,1,1)= 82446.86  ! 82443.0, 82445.9, 82449.2
c$$$c     tF
c$$$      enAr(1,3,1,-1)= 77259.17 ! 77237.04, 77239.2, 77286.99
c$$$      enAr(2,3,1,-1)= 79746.1  ! 79743.7, 79745.0, 79748.3
c$$$      enAr(3,3,1,-1)= 81105.5  ! 81103.9, 81104.6, 81107.1
c$$$      enAr(4,3,1,-1)= 81925.6  ! 81923.5, 81924.3, 81927.8
c$$$      enAr(5,3,1,-1)= 82455.8  ! 82454.9, 82455.2, 82456.7
c     tG
c      enAr(1,4,1,1)= 36430.12
      
c     SD
c      enAr(1,2,0,-1)= 23074.416
c      enAr(2,2,0,-1)= 37077.50
c     TD
c      enAr(1,2,1,-1)= 24672.85
c      enAr(2,2,1,-1)= 37172.37
c     TP
c      enAr(1,1,1,1)= 23693.76
c      enAr(2,1,1,1)= 35227.655
c      enAr(3,1,1,1)= 38120.43
c     SP
c      enAr(1,1,0,1)= 36902.55
c     SF
c      enAr(1,3,0,1)= 36165.311
c     TF
c      enAr(1,3,1,1)= 21335.012 
c      enAr(2,3,1,1)= 37502.97
      

c     Convert to ionization energies (Ar^+ 6s) in eV.
c     Ground state ionization energy is : 15.76 eV = 127109.8 cm^{-1}
      gs_wave_number = 127109.8
      do lpar=-1,2,2
         do is=0,1
            do l=0,llmax         
               do n=1,nmax
                  enAr(n,l,is,lpar)=(enAr(n,l,is,lpar)-gs_wave_number)/
     >                 8065.7
               end do
            end do
         end do
      end do
            
      
      return
      end

c**************************************************************************
      subroutine exp_osc(atom,n,np,l,lp,is,lpar,exp_osc_l,exp_osc_v)
      parameter (nmax=8,llmax=5)
      common /oscMg/ oscMg_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
     >   oscMg_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /oscBe/  oscBe_l(nmax,nmax,0:llmax,0:llmax,0:1,-1:1),
     >   oscBe_v(nmax,nmax,0:llmax,0:llmax,0:1,-1:1)
      common /oscHe/ wiesef(nmax,nmax,0:llmax,0:1),
     >   cefs(nmax,nmax,0:llmax,0:1)
      character  atom*20
      
      exp_osc_l = 1e30
      exp_osc_v = 1e30
      if(n.gt.nmax.or.np.gt.nmax.or.l.gt.llmax.or.lp.gt.llmax) 
     >   return
      
      if(atom.eq.'helium') then
         if(lp.ge.l) then
            exp_osc_l = cefs(n,np,l,is)
            if(exp_osc_l.eq.0.0) exp_osc_l = wiesef(n,np,l,is)
         end if
         if(lp.lt.l) then
            tmp = (2.0*lp+1.0)/(2.0*l+1.0)
            exp_osc_l = cefs(np,n,lp,is)*tmp
            if(exp_osc_l.eq.0.0) exp_osc_l = wiesef(np,n,lp,is)*tmp
         end if
      else if(atom.eq.'beryllium') then
         exp_osc_l = oscBe_l(n,np,l,lp,is,lpar)
         exp_osc_v = oscBe_v(n,np,l,lp,is,lpar)
      else if(atom.eq.'magnesium') then
         tmp = (2.0*is + 1.0)*(2.0*l+1.0)
         exp_osc_l = oscMg_l(n,np,l,lp,is,lpar)/tmp
         exp_osc_v = oscMg_v(n,np,l,lp,is,lpar)/tmp
      end if
      
      return
      end
c**************************************************************************
      subroutine print_energy(Nstmax,E,enionry,atom)
      use vmat_module, only: nodeid
      include 'par.f'     
      parameter (nmax=8,llmax=5)
      double precision  E(KNM),tmp
      common /enNe/ enNe(nmax,0:llmax,0:1,-1:1)
      common /enAr/ enAr(nmax,0:llmax,0:1,-1:1)
      common /enHg/ enHg(nmax,0:llmax,0:1,-1:1)
      common /enYb/ enYb(nmax,0:llmax,0:1,-1:1)
      common /enBa/ enBa(nmax,0:llmax,0:1,-1:1)
      common /enZn/ enZn(nmax,0:llmax,0:1,-1:1)
      common /enCd/ enCd(nmax,0:llmax,0:1,-1:1)
      common /enGaII/ enGaII(nmax,0:llmax,0:1,-1:1)
      common /enCa/ enCa(nmax,0:llmax,0:1,-1:1)
      common /enSr/ enSr(nmax,0:llmax,0:1,-1:1)
      common /enMg/ enMg(nmax,0:llmax,0:1,-1:1)
      common /enOV/ enOV(nmax,0:llmax,0:1,-1:1)
      common /enBe/ enBe(nmax,0:llmax,0:1,-1:1)
      common /enHe/ enHe(nmax,0:llmax,0:1,-1:1)
      common /enLiII/ enLiII(nmax,0:llmax,0:1,-1:1)
      common /enHe_th/ enci(nmax,0:llmax,0:1)
      common /Nsarray/ Ns(0:lomax,0:1,komax,2)
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)    
      character atom*20
      character chan(knm)*3
      common /charchan/ chan
      common /major_config/ l1orb(KNM),l2orb(KNM),n1orb(KNM),
     >   n2orb(KNM),config_maj(KNM)
      character lorb(0:23)*1
      data lorb /"s","p","d","f","g","h","i",
     >     "j","k","l","m","n","o","P","q","r",
     >     "s","t","u","v","w","x","y","z"/      

      
c     find  minimum energy of the states:
      tmp =0d0
      lmax_a = -1
      do Nst=1,Nstmax
         tmp = min(E(Nst),tmp)
         lmax_a = max(lmax_a,ll(Nst))
      end do
      if(lmax_a .eq. -1) then
         print*, 'Stop in print_energy(), file osc.f, ',
     >      'value of lmax_a = -1'
         stop
      end if
      if (nodeid.eq.1)
     >write(4,'(4X,"excitation energy in au and in eV,",
     >     " calculated and experimental ionization energies (eV)")')
      do is=0,1
         do il=0,lmax_a ! llmax
            do ip=-1,1,2
               jp = 1
               if((-1)**il.ne.ip) jp = 2
               do n=1,nmax
                  Nst = Ns(il,is,n,jp)
                  if(Nst.gt.0) then
                     en_eV = 27.2116 * (E(Nst) - enionry/2.0)
                     ex_en = E(Nst) - tmp
                     ex_en_eV = ex_en * 27.2116
                     en_exp = 1e10  ! set to large value to print **** if it is not available in the data base (enBe() or enMg() ...)
                     l1 = l1orb(Nst)
                     l2 = l2orb(Nst)
                     n1 = n1orb(Nst)
                     n2 = n2orb(Nst)
                     imode = 1
                     if(en_eV.lt.0.0 .and. il .le. llmax) then
                        if(atom.eq.'helium') then
                           en_exp = enHe(n,il,is,ip) ! -(enci(n,il,is)-2.0)*27.2116
                        else  if(atom.eq.'LiII') then
                           en_exp = enLiII(n,il,is,ip)
                        else  if(atom.eq.'beryllium') then
                           en_exp = enBe(n,il,is,ip)
                        else  if(atom.eq.'oxygenV') then
                           en_exp = enOV(n,il,is,ip)
                        else  if(atom.eq.'neon'.or.atom.eq.'Neon')then
                           en_exp = enNe(n,il,is,ip)
                        else  if(atom.eq.'magnesium') then
                           en_exp = enMg(n,il,is,ip)
                        else  if(atom.eq.'calcium') then
                           en_exp = enCa(n,il,is,ip)
                        else  if(atom.eq.'zinc') then
                           en_exp = enZn(n,il,is,ip)
                        else  if(atom.eq.'cadmium') then
                           en_exp = enCd(n,il,is,ip)
                        else  if(atom.eq.'galliumII') then
                           en_exp = enGaII(n,il,is,ip)
                        else  if(atom.eq.'strontium') then
                        else  if(atom.eq.'strontium') then
                           en_exp = enSr(n,il,is,ip)
                        else  if(atom.eq.'barium') then
                           en_exp = enBa(n,il,is,ip)
                        else  if(atom.eq.'ytterbium') then
                           en_exp = enYb(n,il,is,ip)
                        else  if(atom.eq.'mercury') then
                           en_exp = enHg(n,il,is,ip)
                        else  if(atom.eq.'argon') then
                           en_exp = enAr(n,il,is,ip)
                        else
                           imode = 0
                        end if
                        if(imode .eq. 1) then
                           if (nodeid.eq.1)
     >write(4,'(I3,I5,A7," (",2(I2,A1),")",
     >                          4(2X,F10.5))')
     >                          n,Nst,chan(Nst),n1,
     >                          lorb(l1),n2,lorb(l2),
     >                          ex_en,ex_en_eV,en_eV,en_exp
                        else
                           if (nodeid.eq.1)
     >write(4,'(I3,I5,A7," (",2(I2,A1),")",
     >                          3(2X,F10.5),
     >                          (2X,"**********"))')
     >                          n,Nst,chan(Nst),n1,
     >                          lorb(l1),n2,lorb(l2),
     >                          ex_en,ex_en_eV,en_eV
                           
                        endif
                     end if
                  end if
               end do
            end do
         end do
      end do

      return
      end




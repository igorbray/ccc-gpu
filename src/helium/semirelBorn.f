c     This program will use amplitudes from non-relativistic CCC calculations and
c     mixing coef. between (non-relativistic)  states associated with those amplitudes
c     in order to calculate amplitude in intermediate coupling. This amplitude
c     can be used than to calculates all sort of scattering parameters (DCS,ICS, EICPs, ...).
     
c     THE PROGRAM CHECK if file "input.semirelBorn"  exist and if yes
c     than proceeds to calculate FBA DCS and ICS in intermediate coupling.
c     Program reads from file  "input.semirelBorn" input data for each channel pair
c     1. Initial state angular momentum (J), parity (+1 or -1), and label (s5D)
c     and then on the same line
c     2. Final state angular momentum (J), parity (+1 or -1), and label (t6P)
c     3. Program returns when END OF FILE is detected while reading from file "input.semirelBorn"
      
c     Mixing coef. of nonrelativistic states due to the spin-orbit correction are in file "mixing_ordered.so"
c     The structure of the  "mixing_ordered.so"  file is:
c     # Ordered: Mixing coef. of nonrelativistic states due to the spin-orbit correction:  J, parity, rank_of_mixing_matrix
c     0.0   -1    9 J, parity, number of config.
c     0.0 -1 t6P   9    -3.747211   -0.505400 :    9.996E-01 t6P   1.877E-02 t8P   1.333E-02 t;P   1.164E-02 t9P   7.898E-03 t>P   5.606E-03 t=P   5.415E-03 t:P  -1.529E-03 t7P   1.118E-03 t<P

c     PROGRAM OUTPUT:
c     1. file cr_s-t.T6P5-s6P3 - the results of the mixing for cross sections.
c     2. file: amp.T6P5-s6P3, where T6P and s6P are labels of the final and initial states, and the last letter is multiplicity of the final and initial states (2J+1).

c     NOTE: 
c     From file with mixing coef. the program will determine which nonrelativistic
c     state amplitudes  are required for calculation. It is you rresponsibility to include all requirerd states in the calculations.
      
c     SOBROUTINES CALLED FROM THIS  PROGRAM:
c     
c     get_LSJM(lmax,f,g,rJf,rJi,amp_file,title,units,iflag,nth,th) -
c            form (LSJM) coupling amp. from non-relativistic amp. in (L M_L S M_S) coupling.
c     amp_write(mmax,fst,nth,theta,rJf,rJi,rkf,rki,state_f,state_i,title,units,cs)
c            opens amp.xxx-yyyy file.
c
c     THIS PROGRAM opens files cr_s-t.xxxx-yyyy  and  amp.xxxx-yyyy
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mix_s_t(maxr, nr, lmax, Nmax, KNM, atom_label, chan, 
     >   la, sa, vdcore)
      IMPLICIT NONE
      integer maxr, nr, lmax, Nmax, KNM
      character chan(KNM)*3, atom_label*3
      integer  la(KNM), sa(KNM)
      real vdcore(maxr,0:lmax)

      real energy
      integer NUNIT
      common /data_for_semirelBorn/ energy, NUNIT
      
      complex  g(0:200,2*lmax,2*lmax,0:1,0:1)
      complex ctmp
      real tmp
      real af(Nmax), ai(Nmax)
      real cr_dcs(0:200), cs(2*lmax,2*lmax), cr_ics
      character  title*80, units*80, file*20
      complex fst(0:200,2*lmax,2*lmax,0:1,0:1)
      integer  KJMi, KJMf, ltf, lti   
      integer nth
      real th(0:200), tmpa(0:200), tmpcs(0:200)
      character label_f(Nmax)*3, label_i(Nmax)*3
      character  state_i*3, state_f*3, hi*1,hf*1
      character  label*3, col*1,  label_transition*9
      integer  Ngive_m(Nmax), Ngive_n(Nmax)
      logical there1
      double precision weight(0:200)
      integer iostatIN, iostatMIX, iostatCR
      integer iparity_i, iparity_f, ipar
      real rki, rkf, rJi, rJf,  rj, sf, si
      integer itest
      real e_fi, en_eV_f, en_eV_i, en_au_f, en_au_i, en_incident_eV
      integer N_tmp, N_amp_f, N_amp_i, N_amp
      integer kMf, kMi, isigf, isigi
      real const, pi
      integer nfst, nist
      integer it, k, nstep, m, n, i
      real CROSSINTr, CROSSINTr1, cr_ics_1,  rintitp,  cr_ics_2,
     >   cr_ics_3
      real ttt(201), aaa(201)
c
      atom_label = 'Ba'
c     Check if file input.semirelBorn exist and if yes open it.
      inquire(file='input.semirelBorn',exist=there1)
      if(.not.there1) return
      print*
      print*, 'Start simirel. Born calculation for channel pairs ',
     >   'defined in the file input.semirelBorn'
      open(300,file='input.semirelBorn',iostat=iostatIN,action='READ')
      if(iostatIN.ne.0) then
         print*, '********  File  input.semirelBorn  was not found'
         stop
      end if
c     Open file with mixing coef.
         open(310,file='mixing_ordered.so',
     >      iostat=iostatMIX,action='READ')
         if(iostatMIX.ne.0) then
            print*, '********  File  mixing_ordered.so  was not found'
            stop
         end if    

c     set rki from other routines
      rki =  sqrt(2.0*energy/27.2116)  ! 1.27909  ! 22.26 eV

      do         
c     Please enter total ang.mom. j, parity and label (3 letters) for the initial state.
c     Please enter  total ang.mom. j, parity and label (3 letters),for the final state.
         read(300,*,iostat=iostatIN,ERR=160,END=170)
     >      rJi, iparity_i, state_i, rJf, iparity_f, state_f
 160     if(iostatIN .ne. 0) then
            print*,  'ERROR reading from input.semirelBorn file.'
            stop
         endif

         write(*,'("Initial state: ",F5.1,I5,A6)')
     >      rJi, iparity_i, state_i
         write(*,'("Final   state: ",F5.1,I5,A6)')
     >      rJf, iparity_f, state_f
         
         if(state_f .eq. "") then
            print*,'empty final state label', state_f
            stop
         end if
         if(state_i .eq. "") then
            print*,'empty initial state label', state_i
            stop
         end if
         
         
c     find state with given label state_f, and rJf, iparity_f in the file with mixing coef. and read these coef. and corresponding labels. 
         itest = 0              ! set value itest=0 to check if it is changed later
         read(310,*)
 12      read(310,*,end=13) rj, ipar, N_tmp
         do k=1,N_tmp
            if((rj .eq. rJf) .and. (ipar .eq. iparity_f)) then
               read(310,*,end=13) rj,ipar,label,N_amp_f,en_eV_f,en_au_f,
     >            col,(af(i), label_f(i), i=1,N_amp_f)
               if(label .eq. state_f)   then
                  write(*,*) rj,ipar,' ',label,N_amp_f,
     >               (af(i),' ',label_f(i),i=1,N_amp_f)
                  itest = 1
                  goto 13       ! found the final state
               end if
            else
               read(310,*)
            end if
         end do
         goto 12      
 13      if(itest .eq. 0) then
            print*,'Could not find final state ', state_f,
     >         ' J=', rJf, '  parity=',iparity_f,
     >         ' Wrong reading from file mixing_ordered.so'
            stop
         end if
c     
         rewind(310)
c     find state with given label state_i, and rJi, iparity_i in the file with mixing coef. and read these coef. and corresponding labels. 
         itest = 0              ! set value itest=0 to check if it is changed later
         read(310,*) 
 112     read(310,*,end=113) rj, ipar, N_tmp
c         print*, rj, ipar, N_tmp
         do k=1,N_tmp
            if((rj .eq. rJi) .and. (ipar .eq. iparity_i)) then
               read(310,*,end=113) rj,ipar,label,N_amp_i,en_eV_i,
     >            en_au_i,col,(ai(i), label_i(i), i=1,N_amp_i)
c               read(310,'(F4.1,I3,A4,I4)',end=113,ADVANCE='NO') 
c     >            rj,ipar,label,N_amp_i
c               print*, label,N_amp_i
c               read(310,*,end=113) en_eV_i,en_au_i,col,
c     >            (ai(i), label_i(i), i=1,N_amp_i)
               if(label .eq. state_i) then
                  write(*,*) rj,ipar,' ',label,N_amp_i,
     >               (ai(i),' ',label_i(i),i=1,N_amp_i)
                  itest = 1
                  goto 113      ! found the initial state
               end if
            else
               read(310,*)
            end if
         end do
         goto 112
 113     if(itest .eq. 0) then
            print*,'Could not find initial state ', state_i,
     >         ' J=', rJi, '  parity=',iparity_i,
     >         ' Wrong reading from file  mixing_ordered.so'
            stop
         end if
         
         rewind(310)            ! set file to the begining for he next channel pair.
         
c     Energy difference (in eV) between initial and final states. 
         e_fi = en_eV_f - en_eV_i
         
c     find value of rkf
         tmp = rki*rki - 2.0*abs(en_au_f - en_au_i)
         if( tmp .le.0 ) then
            print*, 'This state is closed at given incident energy'
                                ! goto next channel pair
            go to 150
         endif
         rkf = sqrt(tmp)
         
c     
         do k=1,Nmax
            Ngive_m(k) = 0
            Ngive_n(k) = 0
         end do
c     
         N_amp = 0
         tmp = 0.0
         
         KJMf = nint(2.0*rJf + 1)
         KJMi = nint(2.0*rJi + 1)
         
         const = (2.0*rJi + 1.0)*2.0
         fst = (0.0,0.0)        ! set amp. to zero         
         cr_dcs = 0.0 ! set crsection to zero
         
         
         do n=1,N_amp_i
c     find state number for label label_i
            itest = 0
            do nist=1,Nmax
c               print*, label_i(n),' ',chan(nist)
               if(label_i(n) .eq. chan(nist)) then
                  itest = 1
                  exit
               endif
            end do
            if(itest .eq. 0) then
               print*, 'Could not find label label_i(',n,')=',
     >            label_i(n)
               stop
            endif
c            print*, 'Found:', chan(nist), ' ', nist
            do m=1,N_amp_f               
               g = 0            ! set Born amp. for given nonrel. component to zero
c     find state number for label label_f
               itest = 0
               do nfst=1,Nmax
                  if(label_f(m) .eq. chan(nfst)) then
                     itest = 1
                     exit
                  endif
               end do
            if(itest .eq. 0) then
               print*, 'Could not find label label_f(',m,')=',
     >            label_f(m)
               stop
            endif
c               print*, 'Found:', chan(nfst), ' ', nfst
c     Here  nfst  and  nist  give state numbers for the channel pair nonrel. components
c     Calculate Born amplitude for this channel pair: only if they have the same spin
c               print*, sa(nfst), ' ',sa(nist)
               if(sa(nfst) .ne. sa(nist)) cycle
               ltf = la(nfst)
               lti = la(nist)
               sf = sa(nfst)
               si = sa(nist)
c               print*,'ready for get_LSJM'
               call  get_LSJM(KNM,lmax, g, rJf, rJi, nfst, nist,
     >            ltf, lti, sf, si, rkf, rki, nth, th,maxr,nr,vdcore)
               
               N_amp = N_amp + 1
               Ngive_n(N_amp) = n
               Ngive_m(N_amp) = m               
c     
               do isigi=0,1
                  do isigf=0,1
                     do kMi=1,KJMi
                        do kMf=1,KJMf
                           do it=0,200
                              ctmp = g(it,kMf,kMi,isigf,isigi)
                              fst(it,kMf,kMi,isigf,isigi) =
     >                           fst(it,kMf,kMi,isigf,isigi) +
     >                           af(m)*ai(n)*ctmp
                           end do
                        end do
                     end do
                  end do
               end do
               
            end do
         end do
         
c     Fix units here. Currently you have atomic units.
         units =  "and is in units of: a.u. (a_0^2)"
         if(NUNIT .eq. 2) then
c     To change to \pi a_0^2
            units =  "and is in units of:  \pi a_0^2"
            pi = acos(-1.0)
            fst = fst / pi
         elseif(NUNIT .eq. 3) then
c     To change to cm**2
            units =  "and is in units of:  cm**2"
            fst = fst * 0.28e-16
         endif
c     
c     Prepare weights for integration over scattering angles
         nth = 201
         nstep = 1      
         do it = 0, 200, nstep
            th(it) = dfloat(it)
            if(it.gt.180) th(it) = (it - 180) / 21d0
         enddo
         
         weight(0) = 1.0/21.0
         k = 2
         do it = 181, 200, nstep
            k = 6 - k
            weight(it) = real(k)/21.0
         enddo
         weight(1) = 1.0 + 1.0/21.0
         do it = 2, 178, nstep
            k = 6 - k
            weight(it) = k
         enddo
         weight(179) = 2.5
         weight(180) = 1.5
c
         aaa(1) = th(0)
         aaa(2:21) = th(181:200)
         aaa(22:201) = th(1:180)
c         
         do kMi=1,KJMi
            do kMf=1,KJMf
c     print*, 'kMf=', kMf
               tmpcs = 0.0
               do isigi=0,1
                  do isigf=0,1
                     do it=0,200
                        tmp = cabs(fst(it,kMf,kMi,isigf,isigi))*
     >                     cabs(fst(it,kMf,kMi,isigf,isigi))/const
                        cr_dcs(it)  = cr_dcs(it) + tmp
                        tmpcs(it) = tmpcs(it) + tmp
                     end do
                  end do
               end do
c               cs(kMf,kMi) =  CROSSINTr(tmpcs,th,weight)
               ttt(1) = tmpcs(0)
               ttt(2:21) = tmpcs(181:200)
               ttt(22:201) = tmpcs(1:180)
               cs(kMf,kMi) = rintitp(ttt,rki,rkf,nth,aaa)
            end do
         end do
c
         ttt(1) = cr_dcs(0)
         ttt(2:21) = cr_dcs(181:200)
         ttt(22:201) = cr_dcs(1:180)
         cr_ics = rintitp(ttt,rki,rkf,nth,aaa)
c

         if(KJMi .lt.10) then
            hi = CHAR(ICHAR('0')  + KJMi)
         else
            hi = CHAR(ICHAR('a') - 10 + KJMi)
         end if
         if(KJMf .lt.10) then
            hf = CHAR(ICHAR('0')  + KJMf)
         else
            hf = CHAR(ICHAR('a') - 10 + KJMf)
         end if
         label_transition = state_f//hf//'-'//state_i//hi
         
         file = 'cr_s-t.'//label_transition
         open(20,iostat=iostatCR,file=file)
         if(iostatCR.ne.0) then
            print*,' can not open file: ', file
            stop
         else
            print*,' open file: ', file                  
         end if
         
         
         en_incident_eV = 27.2116*rki*rki/2.0
         write(title,'("# electron - ",A4,"  scattering: ",
     >      F8.3,"eV on ",A3," ->",A3,F8.3,"eV ",F8.3,"eV")')
     >      atom_label, en_incident_eV, state_i, state_f,
     >      en_eV_i, en_incident_eV +  en_eV_i
         write(20,'("# ",A80)') title
         write(20,'("#",A80)')  units
         
         write(20,'("# Diff. cross section between fine-structure ",
     >      "levels J_f=",F5.1," <-- J_i=",F5.1)')  rJf, rJi
         write(20,'("# ICS:",1p,E12.4)')  cr_ics
         do it=0,180
            write(20,'(F5.1,1p,E12.4)') real(it), 
     >         cr_dcs(it)
         end do
         close(20)
         
         nth = 181
         call amp_write(lmax,fst(0:181,:,:,:,:),nth,th(0:180),rJf,rJi,
     >      rkf,rki,state_f,state_i,title,units,cs)

 150     print*
         
      enddo

      
 170  print*, 'Finished simirel. Born calculation'   
      print*
      
      return
      end

c     This subroutine form (LSJM) coupling amp. from
c     non-relativistic amp. in (L M_L S M_S) coupling.
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  get_LSJM(KNM, lmax, g, rJf, rJi, nf, ni,
     >   ltf, lti, sf, si, rkf, rki, nth, th, maxr, nr, vdcore)
c     lmax - in
c     g(0:200,2*lmax,2*lmax,0:1,0:1) - out,
c                this is non-relativistic amp. in (LSJM) coupling
c     rJf - in, ang.mom of the final state
c     rJi - in, ang.mom of the intial state
c     nfst, nist - in, are state numbers of the nonrelativistic states
c     ltf, lti - in, orbital angular moemntum of the nonrelativistic states
c     sf, si - in, orbital angular moemntum of the nonrelativistic states
c     rki, rkf - in, are initial and final electron energies for the semirel. states
c     nth - out, must be 201 for consistency with ogetBornamp routines.
c     th(0:200) -out, array of angles at which amplitudes are calculated, should be 1 degree 
c                   interval from 0 to 180, for consistency with other routines.
c     
c     This program is written for conversion of amplitudes in $L M_L S M_S$
c     coupling scheme to $(LS)J M_J$ representation.
c      
c     Note: factor  sqrt(ki/kf)  is included in the amplitudes.

c     SOUBROTINES CALLED FROM THIS PROGRAM:
c     getBornamp(nf,lf,psif,maxpsif,rkf,ni,li,
c     >   psii,maxpsii,rki,hlike,nze,nstep,vdcore,minvdc,maxvdc,
c     >   lexit,lentr,thfac,Bornamp)
c     this is where non-relativistic amp. is calculated
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      integer KNM, lmax
      complex  g(0:200,2*lmax,2*lmax,0:1,0:1)
      real rJf, rJi
      integer nf, ni
      integer ltf, lti
      real rki,rkf
      integer nth
      real  th(0:200)
      integer maxr, nr
      real vdcore(maxr,0:lmax)

      complex  Bornamp(0:200,-lmax:lmax,-lmax:lmax)       
      real  si, rli, sf, rlf, rjmin, rjmax
      real rms, rmf, rmi, sig, rJMf, rJMi
      integer KJMf, KJMi, kMf, kMi, isig, is, maxis, mf, mi, it
      real tmp, const
      real CGC
c     This is data required to perform call to routine  getBornamp
      logical hlike
      real psif(maxr), psii(maxr)
      integer minvdc,maxvdc
      integer maxpsif, maxpsii   ! is not used here for hlike .eq. .false.
      integer nze, nstep
      integer lexit, lentr
      real thfac
      real pi, y
c
      integer ns, MaxS, iMs, nsm, isigi, isigf
      real Stot(0:1), sigf, sigi, rmsi, rmsf
      character  states*7
c     
      
      pi = acos(-1.0)

      Bornamp(:,:,:) = (0.0,0.0)
      
      hlike = .false.
      lexit = lmax
      lentr = lmax
      nze = -1   ! incident electron
      nstep = 1
      thfac = 0.0
      minvdc = 1
      maxvdc = nr
c----------------------------------------------------- ---------            
      call getBornamp(nf,ltf,psif,maxpsif,rkf,ni,lti,
     >   psii,maxpsii,rki,hlike,nze,nstep,vdcore,minvdc,maxvdc,
     >   lexit,lentr,thfac,Bornamp)
c------------------------------------------- ---------- ---------

      y = sqrt(rkf/rki)  * (2.0 * pi) ** 2 
      Bornamp = y * Bornamp

      
      rli=lti
      rlf=ltf

      rjmin = abs(sf - rlf)
      rjmax = sf + rlf
      if(rJf .gt. rjmax .or. rJf .lt. rjmin) then
         print*, "Can not have this value of j_f=", rJf,
     >      " for amplitude with l_f=",ltf," and s_f=",sf
         stop
      end if
      rjmin = abs(si - rli)
      rjmax = si + rli
      if(rJi .gt. rjmax .or. rJi .lt. rjmin) then
         print*, "Can not have this value of j_i=", rJi,
     >      " for amplitude with l_i=",lti," and s_i=",si
         stop
      end if

c
     
c     Total spin:
      nsm = 1
      Stot(0) = abs(si-0.5)
      Stot(1) = abs(si+0.5)
      if(si.eq.0.0) nsm = 0

      
      states = 't7S-t6P'
      call ampNR_write(lmax,Bornamp(0:180,:,:),181,
     >   Nsm,sf,si,rlf,rli,rkf,rki,states)
      
c      
      
c     fixed:
c     rJf, rJMf, sigf
c     rJi, rJMi, sigi

c      isigi = nint(sigi + 0.5)
c      isigf = nint(sigf + 0.5)
c     NOTE for Born: sigf = sigi = sig,  and   rmsf = rmsi = rms
c     rJMi and rJMf are magnetic sublevels for target total angular momentum Ji and Jf.
      KJMf = nint(2.0*rJf + 1)
      KJMi = nint(2.0*rJi + 1)

      const = real(KJMi)*2.0
      
      do kMi=1,KJMi
         rJMi = rJi - (kMi - 1)
         do kMf=1,KJMf
            rJMf = rJf - (kMf - 1)

            
            do isigi=0,1
               sigi = isigi - 0.5
               do isigf=0,1                        
                  sigf = isigf - 0.5
                 write(*,'(4F5.1)')
     >               sigf, sigi, rJMf, rJMi
c     do isig=0,1
c     sig = isig - 0.5
                  
                  do ns=0,nsm   ! total spin loop
                     maxS = nint(2.0*Stot(ns)+ 1.0)
                     do iMs=1,maxS
                        rMs = Stot(ns) - iMs + 1
                        rmi = rJMi + sigi - rMs 
                        rmf = rJMf + sigf - rMs
                        mi = nint(rmi)
                        mf = nint(rmf)
                        rmsi = rMs - sigi
                        rmsf = rMs - sigf
                        tmp = CGC(rlf,rmf,sf,rmsf,rJf,rJMf) *
     >                     CGC(rli,rmi,si,rmsi,rJi,rJMi) *
     >                     CGC(0.5,sigf,sf,rmsf,Stot(ns),rMs) *
     >                     CGC(0.5,sigi,si,rmsi,Stot(ns),rMs)
c                        print*,CGC(rlf,rmf,sf,rmsf,rJf,rJMf),
c     >                     CGC(rli,rmi,si,rmsi,rJi,rJMi),
c     >                     CGC(0.5,sigf,sf,rmsf,Stot(ns),rMs),
c     >                     CGC(0.5,sigi,si,rmsi,Stot(ns),rMs)
c                        write(*,'(10F5.1)')  rlf,rmf,sf,rmsf,rJf,rJMf
                        write(*,'("  ",F10.5,10F5.1)')
     >                     tmp,rlf,rmf,sf,rJf,
     >                           rli,rmi,si,rJi
c     print*, '   ',rMs,sigi,rmsi,sigf,rmsf
c     read(*,*)
                        if(tmp.ne.0.0) then
                           do it=0,200
                              g(it,kMf,kMi,isigf,isigi) =
     >                           g(it,kMf,kMi,isigf,isigi) + 
     >                           tmp*Bornamp(it,mf,mi)
                           end do
                        end if
                     end do
                     
                  end do
                  
c$$$  
c$$$               maxis = nint(2*sf) + 1
c$$$               do is=1,maxis
c$$$                  rms = sf - is + 1
c$$$                  rmi = rJMi - rms
c$$$                  rmf = rJMf - rms
c$$$                  mi = nint(rmi)
c$$$                  mf = nint(rmf)
c$$$                  tmp = CGC(rlf,rmf,sf,rms,rJf,rJMf) *
c$$$     >               CGC(rli,rmi,si,rms,rJi,rJMi)
c$$$                  print*, tmp,rlf,rmf,sf,rJf,rJMf,
c$$$     >               rli,rmi,si,rJi,rJMi
c$$$c     print*, '   ',rMs,sigi,rmsi,sigf,rmsf
c$$$c     read(*,*)
c$$$                  if(tmp.ne.0.0) then
c$$$                     do it=0,200
c$$$                        g(it,kMf,kMi,isig,isig) =
c$$$     >                     g(it,kMf,kMi,isig,isig) + 
c$$$     >                     tmp*Bornamp(it,mf,mi)
c$$$                     end do
c$$$                     print*, kMf,kMi,isig,isig,
c$$$     >                  g(10,kMf,kMi,isig,isig) 
c$$$                  end if
c$$$               end do                     

               end do
            end do
            
         end do
      end do
      
      return
      end

c****************************************************************
      real function CROSSINTr(xx,th,weight)
      real th(0:200), xx(0:200)
      double precision weight(0:200)

      pi = acos(-1.0)
      
      h=pi/180.
      S = 0.0
      do i=0,200
         S = S + xx(i)*sin(pi*th(i)/180.) * weight(i)
      enddo
      
      CROSSINTr=2.*pi*h*S/3.
            
      return
      end
c****************************************************************
      real function CROSSINTr1(xx,th)
      dimension  xx(181)
      real th(181)
      
      pi = acos(-1.0)
      
      h=pi/180.
      k=2
      S = 0.
      do i=1,179
         k=6-k
         S = S + k*xx(1+i)*sin(pi*th(1+i)/180.)
      end do
      CROSSINTr1 = 2.*pi*h*S/3.
      return
      end

c****************************************************************
      subroutine amp_write(mmax,fst,nth,theta,rJf,rJi,
     >   rkf,rki,state_f,state_i,title,units,cs)
      complex fst(181,2*mmax,2*mmax,0:1,0:1)
      character  file*60, title*80, units*80
      character  state_f*3, state_i*3, label*9, hi*1,hf*1
      real theta(nth), tmp(nth)
      real cs(2*mmax,2*mmax)
c      common /scat_stoks/  KJMi, KJMf
c      
      KJMf = nint(2.0*rJf+1.0)
      KJMi = nint(2.0*rJi+1.0)

      i_Born = 0

      if(KJMi .lt.10) then
         hi = CHAR(ICHAR('0')  + KJMi)
      else
         hi = CHAR(ICHAR('a') - 10 + KJMi)
      end if
      if(KJMf .lt.10) then
         hf = CHAR(ICHAR('0')  + KJMf)
      else
         hf = CHAR(ICHAR('a') - 10 + KJMf)
      end if
      label = state_f//hf//'-'//state_i//hi
      file = 'amp.'//label

      n = 42
      open(n, FILE=file,IOSTAT=iost)
      if(iost.ne.0) then
         print*, 'File amp.',states,' was not created.'
         stop
      else
         print*,'Writing file ', file
      end if

c     Get integrated cross sections
      summed_cs = 0.0
      do kmf=1,KJMf
         do kmi=1,KJMi          
            summed_cs = summed_cs +  cs(kmf,kmi) 
         end do
      end do

      write(n,*) title
      write(n,'(10X," Jf=",F5.1," <-- Ji=",F5.1)')  rJf, rJi
      write(n,'("When the amplitudes are squared the result contains",
     >   " the factor Kf/Ki")')
      write(n,'(A80)') units
      write(n,'("The CCC amplitudes may be read in the following",
     >   " way:")')
      write(n,'(6X,"read(n,*) KJMf, KJMi, nth, rkf, rki")')
      write(n,'(6X,"do kMi=1,KJMi")')
      write(n,'(6X,"   do kMf=1,KJMf")')
      write(n,'(6X,"      read(n,*) kMfprev, kMiprev")')
      write(n,'(6X,"      if (kMfprev.ne.kMf.or.kMiprev.ne.kMi) ",
     >   "stop ""incorrect format"" ")')
      write(n,'(6X,"      do ith = 1, nth")')
      write(n,'(6X,"         read(n,*) theta(ith), (freal(i,j),",
     >   " fimag(i,j), i=0,1), j=0,1)")')
      write(n,'(6X,"         do isigi=0,1")')
      write(n,'(6X,"            do isigf=0,1")')
      write(n,'(6X,"               f(ith,kMf,kMi,isigf,isigi) = ",
     >   "cmplx(freal(isigf,isigi),fimag(isigf,isigi))")')
      write(n,'(6X,"            end do")')
      write(n,'(6X,"         end do")')
      write(n,'(6X,"      end do")')
      write(n,'(6X,"   end do")')
      write(n,'(6X,"end do")')
      
      write(n,'(2i4,i4,2f8.5,
     >   " KJMf, KJMi, nth, rkf, rki;  ICS = ", 1p,e11.3,
     >   "; KJMf=(2.0*Jf+1.0), Mf=(Jf-kMf+1)")')
     >   KJMf, KJMi, nth, rkf, rki,  summed_cs 
      if (KJMi.gt.mmax.or.KJMf.gt.mmax) stop
      do kMi=1,KJMi
         do kMf=1,KJMf
            write(n,'(2i4,2X,1p,e12.5," kMf, kMi, ICS,  Mf,Mi:",
     >         0p,2f5.1)') 
     >         kMf, kMi, cs(kMf,kMi), (rJf-kMf+1), (rJi-kMi+1)
            do ith = 1, nth
               write(n,'(F5.1,1p,4(1X,2E11.3))')
     >            theta(ith), ((real(fst(ith,kMf,kMi,isigi,isigf)),
     >            aimag(fst(ith,kMf,kMi,isigi,isigf)),
     >            isigi=0,1), isigf=0,1)
            end do
         end do
      end do
      close(42)

      return
      end



c     This function is used to integrate DCS. It uses interpolation of
c     DCS * q**2 which allow for better accuracy at small angles at high energies.
c     It uses library interpolation subroutine INTRPL()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      real function rintitp(ttt,rki,rkf,nth,theta)
c     ttt(nth) - in, this is DCS
c     rki, rkf - in, initial and final projectile momentum.
c     nth - in, number of angles
c     theta(nth) -in, array of angles at which DCS are calculated.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      integer nth
      real theta(nth)
      real  ttt(nth)
      real rkf, rki
      real  a(0:10)
      real*8  x(nth), y(nth), u(600), v(600), FINT,q2,th,pi
      integer n(10)
      COMMON/INTEG/ INT,FINT
      
      rintitp = 0.0
      if(nth .eq. 0) return

      pi = acos(-1.)

c     For elastic transition the value of q can be equal to zero.
c     No division by q**2 should be done in this case, and we do not use
c     any extrapolation technique for elastic transitions. Simple Simpson rule
c     is used in routine   CROSSINTr1(ttt,nth,theta).
      if(rki .eq. rkf) then
         nmax = nth
         u = theta
         v = ttt
      else
c
         INT = 0

         do i = 1,nth
            th = pi * real(theta(i))/180.0
c     th = pi * real(i-1)/180.0
            x(i) = th
            q2 = dble(rki*rki) + dble(rkf*rkf) -
     $         dble(2.0*rki*rkf*dcos(th))
            y(i) = dble(ttt(i))  * q2
         end do
         
         mmax = 5.
         a(0) = 0.
         a(1) = 2.
         a(2) = 5.
         a(3) = 10.
         a(4) = 30.
         a(5) = 180.
         n(1) = 40
         n(2) = 30
         n(3) = 25
         n(4) = 40
         n(5) = 150
         i = 1
         u(1) = 0.0
         do m=1,mmax
            do j=1,n(m)
               i = i+1
               u(i) = pi*(a(m-1) + dble(j)*dble(a(m)-
     $            a(m-1))/dble(n(m)))/180d0
            end do
         end do
         nmax = 1
         do m=1,mmax
            nmax = nmax + n(m)
         end do
         call INTRPL(nth,x,y,nmax,u,v)
c     open(1, file='tmp1')
         do i=1,nmax
            q2 =  dble(rki*rki) + dble(rkf*rkf) - 
     $         dble(2.0*rki*rkf*dcos(u(i)))
c     write(1,'(I5,1p,2E12.4)') i, real(u(i))*180.0/pi, real(v(i)/q2)
            v(i) = dsin(u(i)) * v(i) /  q2
         end do
c     close(1)
         
      end if

      
      INT = 1
      call INTRPL(nmax,u,v,nmax,u,v)
      rintitp = 2.0 * pi * FINT

      return
      end


c--------------------------------------------------------------------------------
      subroutine  ampNR_write(mmax,fn,nth,Nsm,sf,si,rlf,rli,
     >   rkf,rki,states)
      complex fn(181,-mmax:mmax,-mmax:mmax)
      character  states*7, file*60, title*65, units*80
      real spinw(0:1)
      real theta(181)

      do i=1,181
         theta(i) = i-1
      enddo
      
      spinw(0) = (2.0*abs(si-0.5)+1)/(2.0*(2.0*si+1))
      spinw(1) = (2.0*abs(si+0.5)+1)/(2.0*(2.0*si+1))
      if(si .eq. 0.0)  spinw(1) = 0.0
      
      lf = nint(rlf)
      li = nint(rli)

      i_Born = 0
      
      n = 42
      file = 'ampNR.'//states(5:7)//'-'//states(1:3)
      open(n, FILE=file,IOSTAT=iost)
      if(iost.ne.0) then
         print*, 'File amp.',states,' was not created.'
         stop
      else
         print*,'Writing file ', file
      end if


      exen = 27.2116*(rki**2 - rkf**2)/2.0
      ii = INDEX(title,states(5:7))
      if(ii .eq. 0) then
         print*,'Can not find label ',states(5:7),' in the title string'
         stop
      end if
      ii = ii + 9
      iif = INDEX(title(ii:),'eV')
      if(iif .eq. 0 ) then
         print*,'Can not find string  eV  in the title string'         
         print*, 'iif=', iif, '  ii=', ii
         print*, title(ii:)
         stop
      end if
      iif= ii -1 + iif - 1
c      print*,ii,iif
c      print*,title
c      print*, title(ii:iif)
      read(UNIT=title(ii:iif),*) en_i

      en_f = en_i + exen

      
      write(n,'(A28,F9.3,"eV on ",A3," ->",A3,2(F9.3,"eV"))')
     >   title(1:28),27.2116*rki**2/2.0, 
     >   states(1:3),states(5:7), en_f, 27.2116*rkf**2/2.0
      
      write(n,*)
      write(n,'("When the amplitudes are squared the result contains",
     >   " the factor Kf/Ki/(2Li+1)")')
      write(n,'(A80)') units
      write(n,'("The CCC amplitudes may be read in the following",
     >   " way:")')
      write(n,'(6X,"read(n,*) lf, li, nth, rkf, rki, tsp1, nsm, ",
     >   "(spinw(ns), ns=0, nsm)")')
      write(n,'(6X,"do mi = - li, li")')
      write(n,'(6X,"   do mf = - lf, lf")')
      write(n,'(6X,"      read(n,*) mfprev, miprev")')
      write(n,'(6X,"      if (mfprev.ne.mf.or.miprev.ne.mi) ",
     >   "stop ""incorrect format"" ")')
      write(n,'(6X,"      do ith = 1, nth")')
      write(n,'(6X,"         read(n,*) theta(ith), (freal(ns),",
     >   " fimag(ns), ns = 0, nsm)")')
      write(n,'(6X,"         f(ith,mf,mi,0) = cmplx(freal(0),",
     >   " fimag(0))")')
      write(n,'(6X,"         f(ith,mf,mi,1) = cmplx(freal(1),",
     >   " fimag(1))")')
      write(n,'(6X,"      end do")')
      write(n,'(6X,"   end do")')
      write(n,'(6X,"end do")')
      
      write(n,'(2i2,i4,2f8.5,2i2,2f8.5,
     >         " Lf,Li,Nth,Kf,Ki,2Sf+1,Nsm,SpinW(ns)")')
     >   lf, li, nth, rkf,rki,KSM,Nsm,(spinw(ns), ns=0,Nsm)
      if (li.gt.mmax.or.lf.gt.mmax) stop
      do mi = - li, li
         do mf = - lf, lf
            write(n,*) mf, mi, 0.0, ' MF, MI, ICS'
            do ith = 1, nth
               if(i_Born.eq.1) then
                  if(Nsm.eq.0) then
c                     write(n,*) theta(ith), a,b,freal(0),fimag(0)
                  else if(Nsm.eq.1) then
c                     write(n,*) theta(ith), a,b,c,d,freal(0),fimag(0)
                  else
                     print*,'Wrong value for Nsm in amp. file'
                     stop
                  end if
               else
                  write(n,'(F5.1,1p,4(1X,2E11.3))')
     >               theta(ith), (real(fn(ith,mf,mi)),
     >               aimag(fn(ith,mf,mi)), ns=0,Nsm)
               end if
            end do
         end do
      end do
      close(42)

      return
      end

c****************************************************************

C  Only the top analyticphoto subroutine is being used
      subroutine analyticphoto(d,npk,nchtop)
      use gf_module
      include 'par.f'
      real d(kmax,nchan),photo(kmax,nchan)
      integer npk(nchtop+1)
      do nch = 1, nchtop
         do kpp=2,npk(nch+1)-npk(nch)
            photo(kpp,nch)=0.0
            do kp=2,npk(nch+1)-npk(nch)
               photo(kpp,nch)=photo(kpp,nch)+d(kp,nch)
     >            *gf(kp,kpp,nch)
            enddo
         enddo
         do kpp=2,npk(nch+1)-npk(nch)
            d(kpp,nch) = photo(kpp,nch)
         enddo
      enddo
      return
      end
      
      SUBROUTINE  alexanalytic(zeff,ucentr,lg,ichi,wk,gk,nchtop,npk,ns,
     >              v,kernel,box,chil,photon)
      use gf_module
      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      common/smallr/ formcut,regcut,expcut,fast
      common /double/id,jdouble(22)
      common/powers/ rpow1(maxr,0:ltmax),rpow2(maxr,0:ltmax),
     >   istartrp(0:ltmax),istoprp(0:ltmax),cntfug(maxr,0:lmax)
      COMMON/dipole/  dr (kmax,nchan),dv (kmax,nchan),dx (kmax,nchan)
      logical fast,positron,photon
      dimension gk(kmax,nchan)
      complex wk(kmax*nchan)
      dimension npk(nchan+1)      
      real kernel(ichi,ichi)
c$$$      real, allocatable ::C(:,:)!A(:,:),B(:,:)
c$$$      real, allocatable ::Rpl(:),Rpg(:)
      real R(maxr),C(maxr,kmax)
      integer :: rpp,rp,rstop,i,j,n,m
      logical closed, box
      real v(ichi,nchtop)
      character :: ch, ch2
      real psi(maxr),renorm,chil(meshr,ichi),ucentr(maxr),Rpl(maxr),
     >   Rpg(maxr),B(kmax,kmax),A(kmax,kmax),vphoto(kmax,nchan)
      complex phase

      ch(n) = char(n + ichar('0'))

      pi= acos(-1d0)
      small = regcut
      istart = 1
C      ucentr(:) = 0d0
      if (npk(nchtop+1)-1.ne.ichi) stop 'ichi ne npk(nchtop+1)-1' 

      print*, "Entering alexanalytic"

      do nch=1,nchtop
!     fix for cray compiler
         ch2 = ch(nch)

         rstop=meshr

         call getchinfo( nch, nt, lg, psi, maxpsi, e, la, na, l)
!     fix for cray compiler
!         print*, "nch: ", ch(nch)," l = ", l," gk(1,"//ch(nch)//"): "
!     >      ,gk(1,nch) 
         print*, "nch: ", ch2," l = ", l," gk(1,"//ch2//"): "
     >      ,gk(1,nch)

c$$$         allocate(!Rpl(maxr),Rpg(maxr),
c$$$     >      C(maxr,npk(nch+1)-npk(nch)))
         Rpl(:) = 0.0
         Rpg(:) = 0.0

         if (gk(1,nch).GE.0.0) then
            eproj = gk(1,nch) * gk(1,nch)
         else 
            eproj = - gk(1,nch) * gk(1,nch)
         endif
         if (eproj.eq.0.0) then
            eta = zeff
         else 
            eta = zeff / sqrt(abs(eproj))
         endif 
         call makegreen(1,l,eproj,eta,ucentr,cntfug(1,l),rmesh,meshr,
     >      jdouble,id,regcut,expcut,Rpl,Rpg,istart,rstop,phase)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP& SCHEDULE(dynamic)
C$OMP& PRIVATE(R,n,m)
         do kpp=2,npk(nch+1)-npk(nch)
            R(1:meshr) = chil(1:meshr,npk(nch)-1+kpp)
            C(:,kpp)=0.0
            call form(R,istart,meshr,Rpl,Rpg,istart,rstop,meshr,
     >         C(istart,kpp),n,m)
         enddo
C$OMP END PARALLEL DO

c$$$         deallocate(Rpl,Rpg)
c$$$         print*,'deallocated RPL and RPG'
         
c$$$         print*,'allocating B with:',npk(nch+1)-npk(nch),
c$$$     >      npk(nch+1)-npk(nch)
c$$$         allocate(B(npk(nch+1)-npk(nch),npk(nch+1)-npk(nch)),stat=istat)
c$$$         print*,'istat:',istat

         B(:,:)=0.0
         const = abs(gk(1,nch))
C  For zero energy GF goes as rho**(l+1)*rho**(-l) leaving a net k 
         if (const.eq.0.0) const = 1.0 !removing / k when it is zero

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP& SCHEDULE(dynamic)
C$OMP& PRIVATE(R,kp)
         do kpp=2,npk(nch+1)-npk(nch)
            do kp=2,npk(nch+1)-npk(nch)
c$$$               call ricbessel(gk(kp,nch),l,rmesh,meshr,jstart,R)
               B(kp,kpp) = - dot_product(chil(1:meshr,npk(nch)-1+kp),
     >            C(1:meshr,kpp)) * pi / const !*(2.0/pi)
!     >            *real(wk(kp+npk(nch)-1))!*real(wk(kpp+npk(nch)-1))
c$$$               R(1:meshr) = chil(1:meshr,npk(nch)-1+kp)
c$$$               B(kp,kpp)=0.0               
c$$$               do rp=1,rstop
c$$$                  B(kp,kpp)=B(kp,kpp)+R(rp) !sin(rmesh(rp,1)*gk(kp,nch))
c$$$     >               *C(rp,kpp) !*rmesh(rp,3)
c$$$               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
C write out the G(k,kp) surface and check symmetry
c$$$         open(42,file='test')
c$$$         do kpp=2,npk(nch+1)-npk(nch)
c$$$            do kp=2,npk(nch+1)-npk(nch)
c$$$               write(42,'(1p,4e12.4,2i4)') gk(kp,nch),gk(kpp,nch),
c$$$     >            B(kp,kpp),(gf(kp,kpp,nch)-B(kpp,kp))/
c$$$     >            (gf(kp,kpp,nch)+B(kpp,kp)),kp,kpp
c$$$            enddo
c$$$            write(42,*)
c$$$         enddo
c$$$         close(42)

c$$$         write(60+nch,'("#",10x,1000f11.6)') (gk(kpp,nch),kpp=
c$$$     >      2,npk(nch+1)-npk(nch))
c$$$         do kp = 2,npk(nch+1)-npk(nch)
c$$$            write(60+nch,'(f11.6,1p,1000e11.3)') gk(kp,nch),(B(kp,kpp)
c$$$     >         ,kpp=2,npk(nch+1)-npk(nch))
c$$$         enddo
c$$$         B(2:npk(nch+1)-npk(nch),2:npk(nch+1)-npk(nch)) =
c$$$     >      GF(2:npk(nch+1)-npk(nch),2:npk(nch+1)-npk(nch),nch)
         do nchf=1,nchtop
c$$$            allocate(A(npk(nchf+1)-npk(nchf),npk(nch+1)-npk(nch)))
            a(:,:)=0.0
            if (box) then
               const = abs(gk(1,nch))
C  For zero energy GF goes as rho**(l+1)*rho**(-l) leaving a net k 
               if (const.eq.0.0) const = 1.0 !removing / k when it is zero
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
               do kf=1,npk(nchf+1)-npk(nchf)
                  do kpp=2,npk(nch+1)-npk(nch)
                     A(kf,kpp)=0.0
                     do kp=2,npk(nch+1)-npk(nch) !off shell only
                        A(kf,kpp)=A(kf,kpp)+kernel(kf+npk(nchf)-1,
     >                     kp+npk(nch)-1)
     >                     *B(kp,kpp)
!     >                     *(-pi)/const ! Box basis
     >                     *real(wk(kp+npk(nch)-1)) ! wk = 1.0
C     >                                       *(-4.0)/(gk(1,nch)*pi)
                     enddo
                  enddo
               enddo
C$OMP END PARALLEL DO
            else
C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
               do kf=1,npk(nchf+1)-npk(nchf)
                  do kpp=1,npk(nch+1)-npk(nch)
                     A(kf,kpp)=0.0
                     do kp=2,npk(nch+1)-npk(nch) !off shell only
                        A(kf,kpp)=A(kf,kpp)+kernel(kf+npk(nchf)-1,
     >                     kp+npk(nch)-1)
     >                     *B(kp,kpp)
!     >                     *real(wk(kp+npk(nch)-1)) !wk = 1.0
C     >                                       *(-pi)/(gk(1,nch))! Box basis
     >                     *(-4.0)/abs(gk(1,nch)*pi)
                     enddo
                  enddo
               enddo
C$OMP END PARALLEL DO

            endif               !box
! Below is used to form the dot product with the solution to get on-shell K    
            v(npk(nch),nchf) = 0.0
            do n=npk(nch)+1,npk(nch+1)-1
               v(n,nchf) = A(1,n-npk(nch)+1)!*sqrt(abs(real(wk(n))))
            enddo

C$OMP PARALLEL DO
C$OMP& SCHEDULE(dynamic)
            do kpp=2,npk(nch+1)-npk(nch)
               do kf=2,npk(nchf+1)-npk(nchf)
                  kernel(kf+npk(nchf)-1,kpp+npk(nch)-1)=-A(kf,kpp)
!     >               *sqrt(abs(real(wk(kf+npk(nchf)-1))*
!     >               real(wk(kpp+npk(nch)-1))))
C     >      *real(wk(kpp+npk(nch)-1))
               enddo
            enddo
C$OMP END PARALLEL DO
c$$$            deallocate(A)
c$$$            kpp=2
c$$$            kf =2
c$$$            print"('nchf,nch,A,B:',2i3,1p,3e11.3)",nchf,nch,
c$$$     >         A(kf,kpp),B(kf,kpp)
         enddo                  !nchf loop

c$$$         print*,'Deallocating C for nch:',nch
c$$$         deallocate(C)          !,B)
c$$$         print*,'C deallocated for nch:', nch
         if (photon) then
            const = abs(gk(1,nch))
C  For zero energy GF goes as rho**(l+1)*rho**(-l) leaving a net k 
            if (const.eq.0.0) const = 1.0 !removing / k when it is zero
C Velocity gauge
            do kpp=2,npk(nch+1)-npk(nch)
               vphoto(kpp,nch)=0.0
               do kp=2,npk(nch+1)-npk(nch)
                  vphoto(kpp,nch)=vphoto(kpp,nch)+dv(kp,nch)
     >               *B(kp,kpp)
!     >               *(-pi)/const ! Box basis
!     >               *real(wk(kp+npk(nch)-1)) ! wk = 1.0
               enddo
            enddo
            do kpp=2,npk(nch+1)-npk(nch)
               dv(kpp,nch) = vphoto(kpp,nch)
            enddo
C Acceleration gauge            
            do kpp=2,npk(nch+1)-npk(nch)
               vphoto(kpp,nch)=0.0
               do kp=2,npk(nch+1)-npk(nch)
                  vphoto(kpp,nch)=vphoto(kpp,nch)+dx(kp,nch)
     >               *B(kp,kpp)
!     >               *(-pi)/const ! Box basis
!     >               *real(wk(kp+npk(nch)-1)) ! wk = 1.0
               enddo
            enddo
            do kpp=2,npk(nch+1)-npk(nch)
               dx(kpp,nch) = vphoto(kpp,nch)
            enddo
C Length gauge            
            do kpp=2,npk(nch+1)-npk(nch)
               vphoto(kpp,nch)=0.0
               do kp=2,npk(nch+1)-npk(nch)
                  vphoto(kpp,nch)=vphoto(kpp,nch)+dr(kp,nch)
     >               *B(kp,kpp)
!     >               *(-pi)/const ! Box basis
!     >               *real(wk(kp+npk(nch)-1)) ! wk = 1.0
               enddo
            enddo
            do kpp=2,npk(nch+1)-npk(nch)
               dr(kpp,nch) = vphoto(kpp,nch)
            enddo
         endif !photon
      enddo                     !nch loop

      print*, "Exited alexanalytic"

      return
      END

C     For testing if the complete set of states is working as
C     intended

      SUBROUTINE alexstatestest(gk,npk,nchtop,wk,weightk,nqm,ns,kernel
     >                                                            ,box)
      include 'par.f'
      complex wk(kmax*nchan)
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension npk(nchan+1)      
      dimension gk(kmax,nchan), weightk(kmax)
      real kernel(npk(nchtop+1)-1,npk(nchtop+1)-1),Vsum(npk(nchtop+1)-1
     >                                                 ,npk(nchtop+1)-1)
      logical box
      character ch*1,ch2

      ch(n) = char(n + ichar('0'))
      pi= acos(-1.)

      print*, "Entering alexstatestest"

      print*, "nchtop = ",nchtop

      do nch=1,nchtop

!     fix for cray compiler
          ch2=ch(nch)
!      print*, "gk(1,"//ch(nch)//")= ",gk(1,nch)
          print*, "gk(1,"//ch2//")= ",gk(1,nch) 

              if (ns.eq.0) then
!     fix for cray compiler                      
!               open(42,file="statestestsingletnch"//ch(nch))
               open(42,file="statestestsingletnch"//ch2)        
              else
!     fix for cray compiler                      
!               open(42,file="statestesttripletnch"//ch(nch))
               open(42,file="statestesttripletnch"//ch2)
              endif

C           Calculates the complete set of states integrals and writes
C           the results to see if it is working well

               if (box) then

                  do kf=1,1
                        do k=1,npk(nch+1)-npk(nch)
                              Vsum(kf,k)=0.0
                              do i=1,meshr
                                    do kp=1,npk(nch+1)-npk(nch)
                                          Vsum(kf,k)=Vsum(kf,k)
     >                       +kernel(kf+npk(nch)-1,kp+npk(nch)-1)
     >                       *sin(rmesh(i,1)*gk(kp,nch))*rmesh(i,3)!   *chil(i,kp+npk(nch)-1)*chil(i,k+npk(nch)-1)
     >                       *sin(rmesh(i,1)*gk(k,nch))
     >                       *real(wk(kp+npk(nch)-1))
C     >                       *2/pi !only for non box basis
                                 enddo
                           enddo                           
                           write(42,*) gk(kf,nch),gk(k,nch)
     >                              ,kernel(kf+npk(nch)-1,k+npk(nch)-1)!/real(wk(j))
     >                              ,Vsum(kf,k)
                        enddo
                  enddo

               else

                  do kf=1,1
                        do k=1,npk(nch+1)-npk(nch)
                              Vsum(kf,k)=0.0
                              do i=1,meshr
                                    do kp=1,npk(nch+1)-npk(nch)
                                          Vsum(kf,k)=Vsum(kf,k)
     >                       +kernel(kf+npk(nch)-1,kp+npk(nch)-1)
     >                       *sin(rmesh(i,1)*gk(kp,nch))*rmesh(i,3)!   *chil(i,kp+npk(nch)-1)*chil(i,k+npk(nch)-1)
     >                       *sin(rmesh(i,1)*gk(k,nch))
     >                       *real(wk(kp+npk(nch)-1))
     >                       *2/pi !only for non box basis
                                 enddo
                           enddo
                           write(42,*) gk(kf,nch),gk(k,nch)
     >                              ,kernel(kf+npk(nch)-1,k+npk(nch)-1)!/real(wk(j))
     >                              ,Vsum(kf,k)
                        enddo
                  enddo

               endif

               close(42)

         enddo

      print*, "Exiting alexstatestest"

      return
      end

      subroutine alexmatinv(kernel,lda,n,v,m)
      include 'par.f'
      parameter (nmax = kmax * nchan)
      real kernel(lda,lda)
      real v
      dimension v(lda,m,2), ipvt(nmax)

      print*,'In ALEXMATINV with LDA, M:',lda,m
      if (n.le.0) return

#ifdef _single
      
      call sgesv(n,m,kernel,lda,ipvt,v,lda,info)

#elif defined _double

      call dgesv(n,m,kernel,lda,ipvt,v,lda,info)

#endif

      if (info.ne.0) print*,'Warning: INFO, N:', info, N

      return
      end

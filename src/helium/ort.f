      subroutine  orthe(Nmax,nspmW,C)
      include 'par.f'
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), nprinc(KNM)
      common /ortog/  ortint(nspmax,nspmax)
      common/orbsp/nspm,lo(nspmax),ko(nspmax),nset(nspmax)
      double precision  C(Nmax+1,nspmW,nspmW), ortint, sum, dr1, dr2
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer na, nam
c
      write(10,'("<N|Np>")') 
      do N=1,Nmax
         do Np=1,N
            sum = 0.0D0
            if(ll(N).eq.ll(Np)) then
               do jn1=1,nam(N)
                  n1 = na(jn1,N)
                  do jn1p=1,nam(Np)
                     n1p = na(jn1p,Np)
                     dr1 = ortint(n1,n1p)
                     if(lo(n1).eq.lo(n1p).and.dr1.ne.0.0D0) then
                        do jn2=1,nam(N)
                           n2 = na(jn2,N)     
                           do jn2p=1,nam(Np)
                              n2p = na(jn2p,Np)          
                              dr2 = ortint(n2,n2p)
                              if(lo(n2).eq.lo(n2p).and.dr2.ne.0d0) then
                                 sum = sum + dr2*dr1*
     >                              C(N,jn1,jn2)*C(Np,jn1p,jn2p)
                              end if
                           end do
                        end do
                     end if
                  end do
               end do
            end if
            write(10,'("N,Np=",2I3,",  <N|Np>=",F10.5,E15.5)')
     >         N, Np, real(sum), real(sum)
         end do
      end do
      
      return
      end


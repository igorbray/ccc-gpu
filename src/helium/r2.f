c************** r**2 for states of HE ****************************************
c*****************************************************************************
      subroutine r2s(Nmax,nspmW,C,fl,minf,maxf,lo)
      include 'par.f'
      dimension  fl(nmaxr,nspmax)
      dimension  minf(nspmax), maxf(nspmax), lo(nspmax)
      common /helium/ ll(KNM), ls(KNM), lparity(KNM), np(KNM)
      double precision  C(Nmax+1,nspmW,nspmW)
      double precision r2
c
      write(4,'("here <N|r**2|N> are calculated")')
      do N=1,Nmax
         rsq = real(r2(Nmax,nspmW,N,C,fl,minf,maxf,lo))/2.
         write(4,'("N=",I3,", l=",I3,", s=",I3,", r**2 =",F10.5)')
     >      N,ll(N),ls(N),rsq
      end do
      return
      end
c
      double precision function r2(Nmax,nspmW,N,C,fl,minf,maxf,lo)
      include 'par.f'
      dimension  fl(nmaxr,nspmax), lo(nspmax)
      dimension  minf(nspmax), maxf(nspmax)
      common /ortog/  ortint(nspmax,nspmax)
      double precision  trm, dr, r1elk, ortint
      double precision  C(Nmax+1,nspmW,nspmW)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      integer na, nam
c     
      r2 = 0.0D0
      trm = 0.0D0
      do jn1=1,nam(N)
         n1 = na(jn1,N)
         do jn1p=1,nam(N)
            n1p = na(jn1p,N)
            if(lo(n1).eq.lo(n1p)) then
               dr = ortint(n1,n1p)
               if(dr.ne.0.0D0) then
c     write(4,'(2I3,",  dr=",F10.5)') n1,n1p,real(dr)
                  do jn2=1,nam(N)
                     n2 = na(jn2,N)
                     do jn2p=1,nam(N)
                        n2p = na(jn2p,N)
                        if(C(N,jn1,jn2).ne.0.0D0.
     >                     and.C(N,jn1p,jn2p).ne.0.0D0.
     >                     and.lo(n2).eq.lo(n2p)) then
                           trm = trm + r1elk(2,fl(1,n2),fl(1,n2p),
     >                        minf(n2),minf(n2p),maxf(n2),maxf(n2p))
     >                        * dr * C(N,jn1,jn2)*C(N,jn1p,jn2p)
                        end if
                     end do
                  end do
               end if
            end if
         end do
      end do
      r2 = trm*2.0D0
      return
      end
c*******************************************************************

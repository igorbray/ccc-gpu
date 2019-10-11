c     This is Fortran 90 file.
      subroutine fit_rgrid(maxr,nnmax,lnabmax,psinb,istoppsinb,
     >   enpsinb,nr,nshells,grid,EXPCUT)
      real psinb(maxr,nnmax,0:lnabmax)
      real enpsinb(nnmax,0:lnabmax)
      integer istoppsinb(nnmax,0:lnabmax)
      integer  nr, nshells, ncspvs, imesh, i, ndf
      integer nnn(nshells), l(nshells)
      real wwnl(nshells)
      double precision  ee(nshells)
      double precision  r(0:nr), pnl(0:nr,nshells)
      real   grid(nr)
      double precision  dgrid(nr)
      double precision  v(0:nr)
      integer ir
      real EXPCUT
      double precision r1elk
c
      r(:) = 0.0
      pnl(:,:) = 0.0
      v(:) = 0.0
      dgrid(:) = 0.0
      print*,'Read HF core wave functions from formated file'
c
      open(129,iostat=iostat,file='hfwf.formated',
     >   status='old',form='formatted',recl=10000)
      if(iostat .ne. 0) then
         print*,'fit_rgrid: Can not open file  hfwf.formated'
         stop
      end if
      read(129,*) ncspvs, imesh
      read(129,*) (nnn(i), i=1,ncspvs)
      read(129,*) (l(i), i=1,ncspvs)
      read(129,*) (wwnl(i), i=1,ncspvs)
      read(129,*) (ee(i), i=1,ncspvs)

      read(129,*)
      imesh = imesh   !   have been changed from imesh = imesh -1
      do i=1,imesh
         read(129,*) r(i),(pnl(i,ndf), ndf=1, ncspvs)
      end do
      close(129)

!     find index to the largest r point in new r-grid.
      rmax = r(imesh)
      do i=1,nr
         if(rmax .lt. grid(i)) then
            imax = i
            exit
         endif
      enddo

      dgrid(1:nr) = dble(grid(1:nr))           
      do ndf=1, ncspvs
         do i=imesh,2,-1
            if(abs(pnl(i,ndf)) .ge. EXPCUT) exit
         end do
         ir = i
         call INTRPL(ir+1,r(0:ir),pnl(0:ir,ndf),
     >        imax,dgrid(1:imax),v(1:imax))
         nn = nnn(ndf)
         ll = l(ndf)
         enpsinb(nn, ll) = ee(ndf)   ! in Ry
         psinb(1:imax, nn, ll) = v(1:imax)
         do i=imax,2,-1
            if(abs(v(i)) .ge. EXPCUT) exit
         end do
         istoppsinb(nn, ll) = i
      end do
      
       do ndf=1, -ncspvs
          do ndfp=1,ndf
             if(l(ndf) .eq. l(ndfp)) then
                nn = nnn(ndf)
                ll = l(ndf)
                i2 = istoppsinb(nn, ll)
                nnp = nnn(ndfp)
                llp = l(ndfp)
                i2p = istoppsinb(nn, ll)
                i2 = min(i2,i2p)
                tmp = r1elk(0,psinb(1,nn,ll),psinb(1,nnp,llp),
     >               1,1,i2,i2)
                print*,ndf,ndfp,ll,tmp
             endif
          enddo
       enddo

      return
      end

      

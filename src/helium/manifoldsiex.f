      subroutine writesmat(lg,ispin,ipar,nchan,smat,maxr,etot)
      complex smat(nchan,nchan)
      dimension temp(maxr)
      character*20  stren
      
      write(stren,'(E10.3)') etot

      nn = 175
      if(lg .eq. 0 .and ispin .eq. 0) then
         open(nn,file='smat_out_'//etot)         
      endif

      write(nn,'("J S Par:",3i5)') lg, ispin, ipar 
      do nchi = 1, nchan
         call getchinfo (nchi,nchip,lg,temp,maxpsi,ei,lia,nia,li)
         do nchf = 1, nchan
            call getchinfo (nchf,nchp,lg,temp,maxpsi,ef,lfa,nfa,lf)
            write(nn,'(6i5,E12.5)') li,lia,nia,lfa,nfa,lf, smat(nchf,nchi)

         enddo
      enddo

      return
      end

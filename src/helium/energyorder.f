      subroutine energyorder(Nmax)
      IMPLICIT NONE
      use CI_MODULE
      integer Nmax
      integer la,sa,lpar,np, na, nam
      common /helium/ la(KNM), sa(KNM), lpar(KNM), np(KNM)
      common /CIdata/ na(nspmCI,KNM), nam(KNM)
      double precision E
      common /CcoefsE/  E(KNM)
      character chan(knm)*3
      common /charchan/ chan

      integer Num(Nmax)
      integer i, j, itmp


c     This array will be sorted by insertion sort algorithm acording to the E(i) array
      do i=1,Nmax
         Num(i) = i
      end do
      do i=2,Nmax
         itmp = Num(i)
         j = i
         do while(j.ge.2.and.E(itmp).lt.E(Num(j-1)))
            Num(j) = Num(j-1)
            j = j - 1
         end do
         Num(j) = itmp
      end do

      
      do i=1,Nmax
         write(*,'(2I5,F15.5)') i, Num(i), E(Num(i))
      end do   



      return
      end
c---------------------------------------------------------------------------
c     This subroutine sorts target states in the enrgy order.
      subroutine energyorder(Nmax, KNM, nspmCI, la, sa, lpar, np, 
     >   na, nam, E, chan)
      use CI_MODULE
      IMPLICIT NONE
      integer Nmax, KNM, nspmCI
      integer la(KNM), sa(KNM), lpar(KNM), np(KNM)
      integer na(nspmCI,KNM), nam(KNM)
      double precision E(KNM)
      character chan(knm)*3

      integer Num(Nmax)
      integer i, j, itmp
      integer natmp(nspmCI,KNM)


c     This array will be sorted by insertion sort algorithm acording to the E(i) array
      do i=1,Nmax
         Num(i) = i
      end do
      do i=2,Nmax
         itmp = Num(i)
         j = i
         do while(j.ge.2.and.E(itmp).lt.E(Num(j-1)))
            Num(j) = Num(j-1)
            j = j - 1
         end do
         Num(j) = itmp
      end do

      
      do i=1,Nmax
         write(*,'(2I5,2X,A5,F15.5)') i, Num(i), chan(Num(i)), E(Num(i))
      end do   


      la(1:Nmax) = la(Num(1:Nmax))
      sa(1:Nmax) = sa(Num(1:Nmax))
      lpar(1:Nmax) = lpar(Num(1:Nmax))
      np(1:Nmax) = np(Num(1:Nmax))
      chan(1:Nmax) = chan(Num(1:Nmax))
      nam(1:Nmax) = nam(Num(1:Nmax))

      natmp(:,:) = na(:,:)
      do i=1,Nmax
         na(:,i) = natmp(:,Num(i))                
      end do   


      return
      end

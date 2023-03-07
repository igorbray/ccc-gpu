#define REQUESTTAG 5
#define DATATAG 6
#if defined _single
#define MY_MPI_REAL MPI_REAL
#define MY_MPI_COMPLEX MPI_COMPLEX
#define MY_GEMR2D psgemr2d
#elif defined _double
#define MY_MPI_REAL MPI_DOUBLE_PRECISION
#define MY_MPI_COMPLEX MPI_DOUBLE_COMPLEX
#define MY_GEMR2D pdgemr2d
#endif
#define MATRIX 0
#define VECTOR 1

      subroutine factor_two(n, a, b)
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: a, b
      integer :: i, j
      loop_a : do i=0, n-1
         do j=0, n-1
            if((n-i) * (n-j) .eq. n) then
               a = n-i
               b = n-j
               if(a .le. b) exit loop_a
            end if
         end do
      end do loop_a
      end subroutine factor_two

      
      subroutine redistributeAndSolve(ni,nf,nd
     >             ,blacs_ctx,ns,nsmax
     >             ,nodes,npklen,npk
     >             ,nchtop
     >             ,wk,wklen,soln)
      use vmat_module
      implicit none
      include 'mpif.h'
c$$$      real, dimension(nf+1:nd,ni:nf) :: vmat0
c$$$      real, dimension(ni:nf,ni:nf+1) :: vmat01
c$$$      real, dimension(ni:nf,nf+1+1:nd+1) :: vmat1
      integer :: ni,nf,nd
      integer :: blacs_ctx
      real, dimension(:,:), allocatable :: kl,vl
      integer, dimension(2) :: kblock,vblock
      integer :: ns 
      integer :: rows,cols,myrow,mycol
      integer, dimension(2) :: kldim,vldim
      integer, dimension(9) :: desc_kl, desc_vl, desc_soln
      integer :: ierr
      integer :: numroc
      integer :: wklen
      complex, dimension(wklen) :: wk
      integer :: i, nsmax
      integer :: nodes,matrix,vector
      integer :: npklen
c$$$      integer, dimension(nodes) :: nchistart
c$$$      integer, dimension(nodes) :: nchistop
      integer, dimension(npklen) :: npk
      integer :: nchtop,omp_get_num_threads,nomp
      real, dimension(nd,nchtop,2) :: soln
      character date*8,time*10,zone*5
      integer valuesin(8), valuesout(8), idiff, myid

!      wk(1) = cmplx(0.0,imag(wk(1)))

      call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
      call blacs_gridinfo(blacs_ctx,rows,cols,myrow,mycol)
      kblock(1)=min(64,nd/rows)
      kblock(2)=min(64,nd/cols)
      kblock(1)=min(kblock(1),kblock(2))
      kblock(2)=kblock(1)
!      vblock(1)=min(64,nd/rows)
!      vblock(2)=1
      vblock=kblock
      vblock(2)=1
      kldim(1)=numroc(nd,kblock(1),myrow,0,rows)+1
      kldim(2)=numroc(nd,kblock(2),mycol,0,cols)+1 
      vldim(1)=numroc(nd,vblock(1),myrow,0,rows)+1
      vldim(2)=numroc(nchtop,vblock(2),mycol,0,cols)+1
      call descinit(desc_kl,nd,nd,kblock(1),kblock(2),
     > 0,0,blacs_ctx,max(kldim(1),1),ierr)
      call descinit(desc_vl,nd,nchtop,vblock(1),vblock(2),
     > 0,0,blacs_ctx,max(vldim(1),1),ierr)
      call descinit(desc_soln,nd,nchtop,nd,nchtop,
     > 0,0,blacs_ctx,nd,ierr)
      allocate(kl(kldim(1),kldim(2)))
      allocate(vl(vldim(1),vldim(2)))

      call date_and_time(date,time,zone,valuesin)
      call distribute(vmat0,vmat1,vmat01,ni,nf,nd,vl,vldim,desc_vl,
     >             nodes,nchistart,nchistop,npklen,npk,
     > blacs_ctx,ns,VECTOR,1,wk,wklen)
      call date_and_time(date,time,zone,valuesout)
      if (myid.eq.0) print '(" exited DISTRIBUTE VECTOR at: ",a10,
     >      ", diff (secs):",i5)', time, idiff(valuesin,valuesout)
      valuesin = valuesout
      call distribute(vmat0,vmat1,vmat01,ni,nf,nd,kl,kldim,desc_kl,
     >             nodes,nchistart,nchistop,npklen,npk,
     > blacs_ctx,ns,MATRIX,1,wk,wklen)
      call date_and_time(date,time,zone,valuesout)
      if (myid.eq.0) print '(" exited DISTRIBUTE MATRIX at: ",a10,
     >      ", diff (secs):",i5)', time, idiff(valuesin,valuesout)
      valuesin = valuesout

C  Save the driving term on to the zeroth process      
      call MY_GEMR2D(nd,nchtop,vl,1,1,desc_vl,soln(1,1,2),1,1,
     >   desc_soln,blacs_ctx)
      call date_and_time(date,time,zone,valuesout)
      if (myid.eq.0) print '(" exited first P*GEMR2D at: ",a10,
     >      ", diff (secs):",i5)', time, idiff(valuesin,valuesout)
      valuesin = valuesout

      if (ns.eq.nsmax) then
         if (allocated(vmat01)) deallocate(vmat01)
         if (allocated(vmat0)) deallocate(vmat0)
         if (allocated(vmat1)) deallocate(vmat1)
c$$$      else
c$$$         if (allocated(vmat0)) deallocate(vmat0)
      endif 
C  Solve the linear equations
C  The following was an unsuccessful attempt to generate single executables
C  for both scalapack and lapack linear equation solvers.
c$$$  nomp = OMP_GET_NUM_THREADS()
c$$$      call OMP_SET_NUM_THREADS(1)
c$$$      call sleepy_barrier(MPI_COMM_WORLD)
      call scalapackSolve(kl,kldim,desc_kl,vl,vldim,desc_vl)
c$$$      call OMP_SET_NUM_THREADS(nomp)
c$$$      call sleepy_barrier(MPI_COMM_WORLD)
      
      call date_and_time(date,time,zone,valuesout)
      if (myid.eq.0) print '(" exited SCALAPACK SOLVE at: ",a10,
     >      ", diff (secs):",i5)', time, idiff(valuesin,valuesout)
      valuesin = valuesout

c$$$      nomp=omp_get_num_threads()
c$$$      call omp_set_num_threads(1)
c$$$      call omp_set_num_threads(nomp)
      call MY_GEMR2D(nd,nchtop,vl,1,1,desc_vl,soln(1,1,1),1,1,
     >   desc_soln,blacs_ctx)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call date_and_time(date,time,zone,valuesout)
      if (myid.eq.0) print '(" exited second P*GEMR2D at: ",a10,
     >      ", diff (secs):",i5)', time, idiff(valuesin,valuesout)
      valuesin = valuesout
c$$$      call sleepy_barrier(MPI_COMM_WORLD)
      deallocate(kl)
      deallocate(vl)
c$$$      call date_and_time(date,time,zone,valuesout)
c$$$      if (myid.eq.0) print '(" exited deallocate(kl) at: ",a10,
c$$$     >      ", diff (secs):",i5)', time, idiff(valuesin,valuesout)

      return
      end subroutine redistributeAndSolve


      
      subroutine genreqcord(lidx,lidxend,gidx,gidxend,nproc,bctx,nodes,
     > nchistart,nchistop,npklen,npk,kblock,reqcord,ns,
     > part,kf)
      ! map a local index to a global index
      ! then determine which processor that global index is on
      ! then determine the extent of a block that is limited by
      ! a) the block size in the distributed matrix
      ! b) the blocking in the undistributed matrix
      ! to this end, I assume that groups rows are distributed over
      ! the mod(nproc,OMP_NUM_THREADS)==0 processors, that rowmap is of
      ! length mpisize and that rowmap(nproc) contains the global index
      ! of the first row stored on nproc. If no rows are stored (i.e.
      ! mod(nproc,OMP_NUM_THREDS).ne.0 then rowmap(nproc) should be very
      ! large
      implicit none
      integer :: nd
      integer, dimension(2) :: lidx
      integer, dimension(2) :: lidxend
      integer, dimension(2) :: gidx
      integer, dimension(2) :: gidxend
      integer :: nproc
      integer :: nodes
      integer :: npklen
      integer, dimension(nodes) :: nchistart
      integer, dimension(nodes) :: nchistop
      integer, dimension(npklen) :: npk
      integer, dimension(2) :: kblock
      integer, dimension(5) :: reqcord
      integer :: bctx
      integer :: ns
      integer, dimension(2) :: pgriddim, pcord
      integer :: endp
      integer :: indxl2g, indxg2p
      integer :: i
      integer :: nf
      integer :: nomp
      integer :: mgidx
      integer, dimension(2) :: limit
      integer :: part
      integer :: kf
      integer :: OMP_GET_MAX_THREADS

      nomp=1 !max(1,OMP_GET_MAX_THREADS()) ! revert for many tasks per node
!      nomp = 8  ! temporary fix for epic scalapack problem'

      call blacs_gridinfo(bctx,pgriddim(1),pgriddim(2)
     >                    ,pcord(1),pcord(2))
      gidx(2)=indxl2g(lidx(2),kblock(2),pcord(2),0,pgriddim(2))
      if (part.eq.VECTOR) then
        gidx(2)=npk(gidx(2))+kf-1
      end if
      gidx(1)=indxl2g(lidx(1),kblock(1),pcord(1),0,pgriddim(1))

      mgidx=min(gidx(1),gidx(2))
      do i=1,nodes
        if (npk(nchistart(i)).le.mgidx .and. 
     >          npk(nchistop(i)+1)-1.ge.mgidx) then
          nproc=i
        end if
      end do
      nf=npk(nchistop(nproc)+1)-1
      nd=npk(npklen)-1
      if (gidx(1).gt.nf) then
        limit(1)=nd
      else
        limit(1)=nf
      end if
      if (gidx(2).gt.nf) then
        limit(2)=nd
      else 
        limit(2)=nf
      end if

      gidxend(2)=gidx(2)
      if (part.ne.VECTOR) then ! if we are distributing the vecotor, we
                               ! we will only want data from the same
                               ! column, if we are distributing the
                               ! matrix, we can take a block of data
        endp=indxg2p(gidxend(2),kblock(2),0,0,pgriddim(2))
        do while (endp.eq.pcord(2).and.
     >     (gidxend(2).le.limit(2)) .and.
     >     (gidxend(2)-gidx(2)).lt.kblock(2))
          gidxend(2)=gidxend(2)+1
          endp=indxg2p(gidxend(2),kblock(2),0,0,pgriddim(2))
        end do
        gidxend(2)=gidxend(2)-1
      end if
      lidxend(2)=lidx(2)+(gidxend(2)-gidx(2))
      gidxend(1)=gidx(1)
      endp=indxg2p(gidxend(1),kblock(1),0,0,pgriddim(1))
      do while (endp.eq.pcord(1).and.
     >   (gidxend(1).le.limit(1)) .and.
     >   (gidxend(1)-gidx(1)).lt.kblock(1))
        gidxend(1)=gidxend(1)+1
        endp=indxg2p(gidxend(1),kblock(1),0,0,pgriddim(1))
      end do
      gidxend(1)=gidxend(1)-1

      lidxend(1)=lidx(1)+(gidxend(1)-gidx(1))

      if (ns.ne.0) then
          reqcord(1)=gidx(1)
          reqcord(2)=gidx(2)
          reqcord(3)=gidxend(1)
          reqcord(4)=gidxend(2)
      else
          reqcord(1)=gidx(1)
          reqcord(2)=gidx(2)
          reqcord(3)=gidxend(1)
          reqcord(4)=gidxend(2)
      end if
      reqcord(5)=ns
      nproc=(nproc-1)*nomp

      end subroutine genreqcord


      subroutine distribute(vmat0,vmat1,vmat01,ni,nf,nd,
     >                      kdist,kldim,desc_kl,
     >             nodes,nchistart,nchistop,npklen,npk,
     >                      bctx,ns,part,kf,wk,wklen)
      implicit none
      include 'mpif.h'

      real, dimension(nf+1:nd,ni:nf) :: vmat0
      real, dimension(ni:nf,ni:nf+1) :: vmat01
      real, dimension(ni:nf,nf+1+1:nd+1) :: vmat1
      integer :: ni,nf,nd
      integer,dimension(2) :: kldim
      real, dimension(kldim(1),kldim(2)) :: kdist
      integer :: ns  
      integer :: bctx
      integer, dimension(9) :: desc_kl, desc_vl
      integer, dimension(9) :: desc_vmaster, desc_kmaster
      integer :: ierr
      integer :: myid
      integer :: i,j
      integer :: OMP_GET_MAX_THREADS
      integer :: nopenmp
      integer :: mpisize
      integer :: myrow, mycol, rows, cols
      integer :: rnctx
      integer :: blacs_pnum
      integer :: blacs_k_block
      integer :: finished
      integer, dimension(2) :: lidx,gidx,gidxend,lidxend
      integer :: lidx1,lidx2,lidxend1,lidxend2
      integer, dimension(2) :: pgriddim,pcord
      integer, dimension(2) :: kblock
      integer :: indxg2l, indxl2g
      real, dimension(:,:), allocatable :: kbuff
      integer, dimension(4) :: buff
      integer :: idx
      integer :: wklen
      complex, dimension(wklen) :: wk
      integer, dimension(:), allocatable :: sends, requests
      integer, dimension(:,:), allocatable :: reqbuffer
      integer, dimension(5) :: reqcord
      logical :: waiting, reqdone
      integer :: gproc
      integer, dimension(2) :: fillrequests
      integer, dimension(MPI_STATUS_SIZE,2) :: fillstats
      integer, dimension(20) :: rowmap
      real, dimension(:,:,:), allocatable :: sendbuffer
      integer :: nodes
      integer :: npklen
      integer, dimension(nodes) :: nchistart
      integer, dimension(nodes) :: nchistop
      integer, dimension(npklen) :: npk
      integer :: part
      integer :: kf 



      kblock(1)=desc_kl(5)
      kblock(2)=desc_kl(6)
      allocate(kbuff(kblock(1),kblock(2)))

      call mpi_comm_size(MPI_COMM_WORLD,mpisize,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)

      allocate(sends(mpisize))
      allocate(requests(mpisize))
      allocate(reqbuffer(5,mpisize))
      allocate(sendbuffer(kblock(1),kblock(2),mpisize))
      do idx=1,mpisize
        call mpi_irecv(reqbuffer(:,idx),size(reqbuffer,1),MPI_INTEGER,
     >     idx-1,REQUESTTAG,MPI_COMM_WORLD
     >     ,requests(idx),ierr)
        sends(idx)=MPI_REQUEST_NULL
      end do
      lidx(1)=1
      lidx(2)=1
      kdist=0.0

      do while (lidx(2).lt.kldim(2))
        lidx(1)=1
        do while (lidx(1).lt.kldim(1))

          ! turn a local coordinate into a global coordinate
          ! work out the dimensions of the block, and request it
          call genreqcord(lidx,lidxend,gidx,gidxend,gproc,bctx,
     >             nodes,nchistart,nchistop,npklen,npk,
     >                kblock,reqcord,ns,part,kf)
          call mpi_isend(reqcord,size(reqcord),MPI_INTEGER,gproc,
     >            REQUESTTAG,MPI_COMM_WORLD,fillrequests(1),ierr)
        kbuff=0.0
          call mpi_irecv(kbuff,kblock(1)*kblock(2)
     >              ,MY_MPI_REAL,gproc,
     >              DATATAG,MPI_COMM_WORLD,fillrequests(2),ierr)
          waiting=.true.
          do while (waiting)
            call service_requests(vmat0,vmat1,vmat01,ni,nf,nd,sends,
     >              kblock,mpisize,requests,reqbuffer,sendbuffer)
            call mpi_test(fillrequests(2),reqdone,fillstats(:,2),ierr)
            if (reqdone) then
              waiting=.false.
            end if
          end do
          if (part.eq.MATRIX) then ! if distributing the matrix
            do j=lidx(2),lidxend(2)
              do i=lidx(1),lidxend(1)
                
                kdist(i,j) = -kbuff(i-lidx(1)+1,j-lidx(2)+1)
c$$$     >                       *sqrt(abs(real(wk(gidx(1)+i-lidx(1)))
c$$$     >                                *real(wk(gidx(2)+j-lidx(2)))))

                if ((gidx(1)+i-lidx(1)).eq.(gidx(2)+j-lidx(2))) then
                  kdist(i,j) = kdist(i,j)
c$$$     >                              +real(wk(gidx(1)+i-lidx(1))+1e-30)
c$$$     >                         /abs(real(wk(gidx(1)+i-lidx(1))+1e-30))
     >                + 1.0/real(wk(gidx(1)+i-lidx(1))+1e-30)
c$$$                  print*,'gidx(1),i,lidx(1),wk:',gidx(1),i,lidx(1),
c$$$     >               real(wk(gidx(1)+i-lidx(1)))
               end if
              end do
            end do
          else ! if distributing the vector
            do j=lidx(2),lidxend(2)
              do i=lidx(1),lidxend(1)
                kdist(i,j) = kbuff(i-lidx(1)+1,j-lidx(2)+1)
c$$$     >                       *sqrt(abs(real(wk(gidx(1)+i-lidx(1)))))
              end do
            end do
          end if
          lidx(1)=lidxend(1)+1
        end do
        lidx(2)=lidxend(2)+1
      end do
      ! signal that this processor has recieved all the data it needs
      reqcord(:)=-1
      do idx=0,mpisize-1
        call mpi_isend(reqcord,size(reqcord),MPI_INTEGER,idx,
     >            REQUESTTAG,MPI_COMM_WORLD,fillrequests(1),ierr)
      end do
      ! spin waiting for all processors to signal completion
      idx=1
      do while (idx.lt.mpisize+1)
        if (requests(idx).ne.MPI_REQUEST_NULL) then
          call service_requests(vmat0,vmat1,vmat01,ni,nf,nd,sends,
     >              kblock,mpisize,requests,reqbuffer,sendbuffer)
          idx=1
        else
          idx=idx+1
        end if
      end do
      call mpi_wait(fillrequests(1),fillstats(:,1),ierr)
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      deallocate(sends)
      deallocate(requests)
      deallocate(reqbuffer)
      deallocate(kbuff)
      deallocate(sendbuffer)
      
      end subroutine distribute

      subroutine service_requests(vmat0,vmat1,vmat01,
     >                          ni,nf,nd,sends,kblock,
     >                          nprocs,reqhandle,reqbuffer,sendbuffer)
      ! check if any processors have requested a block of data
      ! if they have check if the last block of data they requested has
      ! finished sending and send a new block
      ! check if any more processors have requested a block of data and
      ! repeat
      ! if a processor requests a block at -1, its done
      implicit none
      include 'mpif.h'
      real, dimension(nf+1:nd,ni:nf) :: vmat0
      real, dimension(ni:nf,ni:nf+1) :: vmat01
      real, dimension(ni:nf,nf+1+1:nd+1) :: vmat1
      integer :: ni,nf,nd
      integer, dimension(2) :: kblock
      integer :: finished
      integer, dimension(MPI_STATUS_SIZE) :: mpistat
      integer :: ierr
      integer :: nprocs
      integer, dimension(nprocs) :: sends
      integer, dimension(5,nprocs) :: reqbuffer
      integer :: idxsize
      integer, dimension(nprocs) :: reqhandle
      integer :: i,j
      integer :: idx
      logical :: flag
      real, dimension(kblock(1),kblock(2),nprocs) :: sendbuffer
      integer :: ns
      integer, dimension(2) :: gidx, gidxend


      call mpi_testany(nprocs,reqhandle,idx,flag,mpistat,ierr)
      do while (idx.ne.MPI_UNDEFINED)
         if (reqbuffer(1,idx).eq.-1) then
            reqhandle(idx)=MPI_REQUEST_NULL
            call mpi_wait(sends(idx),mpistat,ierr)
            call mpi_wait(reqhandle(idx),mpistat,ierr)
         else
            call mpi_wait(sends(idx),mpistat,ierr)
            call mpi_wait(reqhandle(idx),mpistat,ierr)
            ns=reqbuffer(5,idx) ! spin
            gidx=reqbuffer(1:2,idx) ! top left corner of the region
            gidxend=reqbuffer(3:4,idx) ! bottom right corner
            if (ns .eq. 0) then
               if (gidx(1).gt.nf) then
               sendbuffer(1:gidxend(1)-gidx(1)+1,1:gidxend(2)-gidx(2)+1,
     >               idx) =
     >               vmat0(gidx(1):gidxend(1),gidx(2):gidxend(2))
               else
!              if (gidx(2).gt.nf+1) then
                  if (gidx(2).gt.nf) then
                     sendbuffer(1:gidxend(1)-gidx(1)+1,
     >                  1:gidxend(2)-gidx(2)+1,idx) =
     >                  transpose(vmat0(gidx(2):gidxend(2),
     >                  gidx(1):gidxend(1)))
                  else
                     do j=gidx(2),gidxend(2)
                        do i=gidx(1),gidxend(1)
                           if (i.ge.j) then
                              sendbuffer(i-gidx(1)+1,j-gidx(2)+1,idx)= 
     >                           vmat01(i,j)
                           else
c$$$                       print*,i-gidx(1)+1,j-gidx(2)+1,idx,j,i
                              sendbuffer(i-gidx(1)+1,j-gidx(2)+1,idx)= 
     >                           vmat01(j,i)
                           end if
                        end do
                     end do
                  end if
               end if
            else                !ns eq 1
               if (gidx(1).gt.nf) then
               sendbuffer(1:gidxend(1)-gidx(1)+1,1:gidxend(2)-gidx(2)+1,
     >               idx) =
     >      transpose(vmat1(gidx(2):gidxend(2),gidx(1)+1:gidxend(1)+1))
               else
!              if (gidx(2).gt.nf+1) then
                  if (gidx(2).gt.nf) then
                     sendbuffer(1:gidxend(1)-gidx(1)+1,
     >                  1:gidxend(2)-gidx(2)+1,idx) =
     >                  vmat1(gidx(1):gidxend(1),gidx(2)+1:gidxend(2)+1)
                  else 
                     do j=gidx(2),gidxend(2)
                        do i=gidx(1),gidxend(1)
                           if (i.ge.j) then
                              sendbuffer(i-gidx(1)+1,j-gidx(2)+1,idx) = 
     >                           vmat01(j,i+1)
                           else
                              sendbuffer(i-gidx(1)+1,j-gidx(2)+1,idx) = 
     >                           vmat01(i,j+1)
                           end if
                        end do
                     end do
                  end if
               end if
            end if

            call mpi_isend(sendbuffer(:,:,idx),kblock(1)*kblock(2)
     >         ,MY_MPI_REAL,idx-1,DATATAG
     >         ,MPI_COMM_WORLD,sends(idx),ierr)
            
            call mpi_irecv(reqbuffer(:,idx),size(reqbuffer,1),
     >         MPI_INTEGER,idx-1,REQUESTTAG,MPI_COMM_WORLD
     >         ,reqhandle(idx),ierr)
         end if
         call mpi_testany(nprocs,reqhandle,idx,flag,mpistat,ierr)
      end do

      end subroutine service_requests

      subroutine scalapackSolve(a,dima,desca,b,dimb,descb)
      implicit none
      include 'mpif.h'
      integer, dimension(2) :: dima, dimb
      real, dimension(dima(1),dima(2)) :: a
      real, dimension(dimb(1),dimb(2)) :: b
      integer, dimension(9) :: desca, descb
      integer :: bctx
      integer :: ierr
      integer, dimension(:), allocatable :: ipiv
      allocate(ipiv(dima(1)+desca(4)))
#ifdef _double
!      print*,'Entering PDGESV'
      call pdgesv(desca(4),descb(4),a,1,1,desca,ipiv,b,1,1,descb,ierr)
!      print*,'Exiting PDGESV'
#elif defined _single
!      print*,'Entering PSGESV'
      call psgesv(desca(4),descb(4),a,1,1,desca,ipiv,b,1,1,descb,ierr)
!      print*,'Exiting PSGESV'
#endif
      if (ierr.ne.0) then
        print*,'p*gesv failed with error value',ierr
        stop
      end if
      deallocate(ipiv)

      end subroutine scalapackSolve

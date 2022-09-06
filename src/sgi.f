      subroutine datetime(n)
      return
      end
            
      subroutine mysleep(n)
      return
      end
      
c$$$      subroutine MPI_INIT( ierr )
c$$$      return
c$$$      end
c$$$
c$$$      subroutine  MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
c$$$      return
c$$$      end
c$$$      
c$$$      subroutine  MPI_COMM_SIZE( MPI_COMM_WORLD, ntasks, ierr )
c$$$      return
c$$$      end
c$$$
c$$$      subroutine MPI_FINALIZE(ierr)
c$$$      return
c$$$      end
      
      subroutine openatend(n,file)
      character file*(*)
      open (n,file=file,position='append',status='old')
      return
      end
      
      subroutine memalloc(ptr,m)
      pointer (ptr,v)
      ptr = malloc(m)
      return
      end

      subroutine memfree(ptr)
      pointer (ptr,v)
      call free(ptr)
      return
      end

c$$$      subroutine gettime(utime,stime)
c$$$c     changing to allow gcc build
c$$$c      real tarray(2)
c$$$      real*4 tarray(2)
c$$$      tt = dtime(tarray)
c$$$      utime = tarray(1)
c$$$      stime = tarray(2)
c$$$      end
      
C  This routine returns the CPU time in seconds since last call to CLOCK
      subroutine clock( ctime )
      real*4 tarray(2), etime
      ctime=0.0 !etime(tarray)
      end

      subroutine update(n)
      return
c$$$      if (n.eq.6) then
c$$$         call flush(101,istat)
c$$$      else
c$$$         call flush(n,istat)
c$$$      endif
      end
      function idiff(valuesin,valuesout)
      integer valuesin(8), valuesout(8)
      idiff = (valuesout(3)-valuesin(3))*86400+
     >   (valuesout(5)-valuesin(5))*3600 +
     >   (valuesout(6)-valuesin(6))*60 +
     >   (valuesout(7)-valuesin(7))
      return
      end

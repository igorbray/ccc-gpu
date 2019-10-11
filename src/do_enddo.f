      character line*80, file*20
      integer nend(100)
      np = 0
      n = 0
      print*,'Enter file please'
      read '(a20)', file
      open (20,file=file)
 10   read(20,'(a80)',end=20) line
      i = 1
      do while (i.lt.80.and.line(i:i).eq.' ')
         i = i + 1
      enddo
      if (line(i:i+2).eq.'do '.and.ichar(line(i+3:i+3)).gt.64) then
         np = np + 1
         n = n + 10
         nend(np) = n
         if (n.lt.100) then
            do j = 77, i + 3, -1
               line(j+3:j+3) = line(j:j)
            enddo
            line(i+3:i+3) = char(n/10+ichar('0'))
            line(i+4:i+4) = char(0+ichar('0'))
            line(i+5:i+5) = ' '
         else if (n.lt.1000) then
            do j = 76, i + 3, -1
               line(j+4:j+4) = line(j:j)
            enddo
            line(i+3:i+3) = char(n/100+ichar('0'))
            line(i+4:i+4) = char((n-n/100*100)/10+ichar('0'))
            line(i+5:i+5) = char(0+ichar('0'))
            line(i+6:i+6) = ' '
         else
            stop 'Labels got bigger than 990'
         endif 
      else if (line(i:i+5).eq.'enddo'.or.line(i:i+6).eq.'end do') then
         line(i:i+8) = 'continue'
         if (nend(np).lt.100) then
            line(2:2) = char(nend(np)/10+ichar('0'))
            line(3:3) = '0'
         else if (nend(np).lt.1000) then
            line(2:2) = char(nend(np)/100+ichar('0'))
            line(3:3) = char((nend(np)-nend(np)/100*100)/10+ichar('0'))
            line(4:4) = '0'
         endif 
         np = np - 1
      end if
      i = 80
      do while (line(i:i).eq.' '.and.i.gt.1)
         i = i - 1
      enddo 
      print '(a)', line(1:i)
      go to 10
 20   continue 
      end
      

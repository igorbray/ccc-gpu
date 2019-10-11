	subroutine cuDDOT(chif,chii,temp3,nqmf,nqmi,meshr,kmax,maxi,vmatt)
	integer		:: nqmf,nqmi,meshr,maxi
	double precision :: chif(meshr,nqmf),chii(meshr,nqmi),temp3(meshr)
	double precision :: vmatt(kmax,kmax,0:1),tmp
	integer :: i
!$acc data copyin(chif,chii,vmatt,nqmf,nqmi,meshr,kmax,maxi) copy(vmatt)
!$acc parallel loop reduction(+:tmp)
	do ki=1,nqmi
		do kf=1,nqmf	
			tmp=0.0
			do i = 1, maxi
				tmp=tmp+chif(i,kf)*chii(i,ki)*temp3(i)
			enddo
            		vmatt(kf,ki,0) = vmatt(kf,ki,0) + tmp
            		vmatt(kf,ki,1) = vmatt(kf,ki,1) + tmp
		enddo
	enddo
!$acc end data
	return
	end subroutine


	function cuDOT(a,b,mymin,mymax)
	integer :: n
	double precision :: a(1:(mymax-mymin+1)),b(1:(mymax-mymin+1))
	double precision :: icnt,dummy
	integer 	::i
	n=mymax-mymin+1
	icnt=0
!$acc data present_or_copyin(a(:),b(:),n,i) copyout(icnt)
!$acc parallel loop
	do i=1,n
		icnt=icnt+a(i)*b(i)
	enddo
!$acc end data
!dummy=dot_product(a(:),b(:))
	cuDOT=icnt
	return
	end function


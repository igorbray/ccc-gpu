C  File form.f
C  This routine returns  FORM(R) = int dR' FUN(R') * F(R<) * G(R>), where
C  the range of integration is zero to infinity.
      subroutine form(fun,ifuns1,ifuns2,f,g,ifstart,igstop,irstop,
     >   formf,i1,i2)
c!$acc routine
      include 'par.f'
      dimension fun(maxr),formf(maxr),f(maxr),g(maxr),temp1(maxr),
     >   temp4(0:maxr),temp5(maxr+1),temp2(maxr),temp3(maxr)
!c$omp threadprivate(temp1,temp2,temp3,temp4,temp5)
!c$omp threadprivate(ifunstart,ifunstop,istop,i)
      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
      common/smallr/ formcut,regcut,expcut,fast
      logical fast

      ifunstart=max(ifuns1,ifstart)
      ifunstop=min(ifuns2,igstop)
      formf(:) = 0.0
      if (ifunstop.le.ifunstart+2) then
         i1=2
         i2=1
         return
      end if 
      istop=min(ifunstop,irstop)
      if (fast) then
C  Find the integral for FUN(R') * F(R') from 0 to IRSTOP. The function FUN 
C  already contains the Simson's integration weights. 
         temp4(ifunstart-1)=0.0
         do i=ifunstart,istop
            temp4(i)=temp4(i-1)+fun(i)*f(i)
         end do 
C  Find the integral of FUN(R') * G(R') from infinity to IFUNSTART
         temp5(ifunstop+1)=0.0
         do i=ifunstop,ifunstart,-1
            temp5(i)=temp5(i+1)+fun(i)*g(i)
         end do 

C  Make the form factor
         do i=ifunstart,istop
            formf(i)=temp4(i)*g(i)+temp5(i+1)*f(i)
         end do
      else
C  This method is more accurate but takes longer.
         do i=ifunstart,istop
            temp1(i)=fun(i)*f(i)/rmesh(i,3)
         end do
         call maketemp(ifunstart,istop,temp1,temp2)
         do i=ifunstart,ifunstop
            temp1(i)=fun(i)*g(i)/rmesh(i,3)
         end do
         call maketempb(ifunstart,ifunstop,temp1,temp3)
         do i=ifunstart,istop
            formf(i)=temp2(i)*g(i)+temp3(i)*f(i)
         end do
         temp4(istop)=temp2(istop)
         temp5(ifunstart)=temp3(ifunstart)
      end if 
            
C  IRSTOP may be greater than IFUNSTOP
      i=istop
      do while (abs(formf(i)).gt.formcut.and.i.lt.min(irstop,igstop))
         i=i+1
         formf(i)=temp4(istop)*g(i)
      end do 
      i2=i
      do while (abs(formf(i2)).lt.formcut.and.i2.gt.ifunstart)
         i2=i2-1
      end do 

      i=ifunstart
      do while (i.gt.ifstart.and.abs(formf(i)).gt.formcut)
         i=i-1
         formf(i)=f(i)*temp5(ifunstart)
      end do
      i1=i
      do while (abs(formf(i1)).lt.formcut.and.i1.le.i2)
         i1=i1+1
      end do 
      end

      subroutine nuclear(fun,direct,ifuns1,ifuns2,irstop,u,nz,form)
      include 'par.f'
      dimension fun(maxr),form(maxr),u(maxr)
      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
      Common /POWERS/ rpow1(maxr, 0:ltmax), rpow2(maxr, 0:ltmax),
     >   iminrp(0:ltmax), imaxrp(0:ltmax), cntfug(maxr, 0:lmax)
      common/smallr/ formcut,regcut,expcut,fast
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym

      logical direct
      z=float(nz)
      tmp=0.0
      if (direct) then
         do i=ifuns1,ifuns2
            tmp=tmp+fun(i)
         end do
c$$$         print*, '********* tmp=', tmp
C  The following statement was sensible for one-electron targets, but
C  makes no sense for two-electron targets. Leave commented out or delete!
c$$$  if (tmp.gt.0.5) then
            do i=1,irstop
               form(i) = form(i) - tmp * (u(i) / 2.0 + rpow2(i,0))
            end do
            do while (abs(form(irstop)).lt.formcut.and.irstop.gt.1)
               form(irstop) = 0.0
               irstop = irstop - 1
            end do 
c$$$            if (abs(tmp-1.0).gt.1e-3) print*,
c$$$     >         'Wave functions are not orthogonal, expected 1.0:',tmp
c$$$         else if (abs(tmp).gt.1e-3) then
c$$$            print*,'Wave functions are not orthogonal, expected 0.0:',
c$$$     >         tmp
c$$$         end if             
      else
C  Exchange
         alpha = log(sqrt(expcut)) / rmesh(meshr,1) * 0.0
         do i=ifuns1,ifuns2
            tmp=tmp+fun(i)*(rpow2(i,0)+u(i)/2.0)
         end do
         do i=1,irstop
            form(i)=(form(i) - tmp) * exp(alpha * rmesh(i,1))
         end do
      end if
      end
      
C  The following routine returns the integral of FUN(I) from 
C  RMESH(I) to RMESH(ISTOP) for I=ISTART,ISTOP.
      subroutine maketempb(iistart,iistop,fun,temp)
      include 'par.f'
      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
      common /double/id,jdouble(22)
      dimension fun(maxr),temp(maxr)

      istop=iistop/2*2
      istart=(iistart+1)/2*2
      do i=1,meshr
         temp(i)=0.0
      end do
      if (istop-istart.lt.4) return
      
C  Define the integral at even points of RMESH using Simpson's rule.
C  Define the integral at odd points of RMESH using a combination of
C  Simpson's rule and Newton's 3/8 rule.
      temp(istop)=0.0
      do i=istop,istart+2,-2
         temp(i-2)=temp(i)+(fun(i-2)+4.*fun(i-1)+fun(i))*rmesh(i-1,2)/3.
         temp(i-3)=temp(i)+(fun(i-3)+3.0*fun(i-2)+3.0*fun(i-1)+fun(i))*
     >      0.375*rmesh(i,2)
      end do

C  The above Newton's 3/8 rule fails where dR doubles
C  I is the next point before the doubling. It should be odd.
      do n=2,id-1
         i=jdouble(n)-1
         if (i.lt.istart) cycle
         temp(i)=temp(i-2)-(fun(i-2)+4.0*fun(i-1)+fun(i))*rmesh(i,2)/3.0
         if (i.eq.istart+1) temp(i)=temp(istart)
      end do

C  The following point had been left out
      temp(istop-1)=temp(istop-3)-(fun(istop-3)+4.0*fun(istop-2)+
     >   fun(istop-1))*rmesh(istop-1,2)/3.0

C  FUN(R) from R=RMESH(istop,1) to 'infinity' is assumed to be zero.
      do i=istop,meshr
         temp(i)=0.0
      end do

      if (iistart.gt.1) then
C  FUN(R) from R=0 to RMESH(istart-1,1) is assumed to be zero.
         do i=1,istart-1
            temp(i)=temp(istart)
         end do
      end if
      end

C  The following routine returns the integral of FUN(RP) from 0 to
C  RMESH(IRP) for IRP=1,IRPSTOP.
      subroutine maketemp(istart,istop,fun,temp)
      include 'par.f'
      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
      common /double/id,jdouble(22)
      dimension fun(maxr),temp(maxr)
      irpstart=(istart+1)/2*2
      irpstop=istop/2*2
      if (irpstart.eq.istart) fun(istart-1)=0.0
      
      do i=1,meshr
         temp(i)=0.0
      end do
      if (irpstop-irpstart.lt.4) return
         
      if (istart.eq.1) then
C  We assume that fun(0)=0. 
         temp(2)=(4.0*fun(1)+fun(2))*rmesh(1,2)/3.0
         temp(3)=(3.0*fun(1)+3.0*fun(2)+fun(3))*rmesh(1,2)*0.375
         temp(1)=temp(3)-(fun(1)+4.0*fun(2)+fun(3))*rmesh(1,2)/3.0
      else
         temp(irpstart)=(4.0*fun(irpstart-1)+fun(irpstart))*
     >      rmesh(irpstart-1,2)/3.0
         temp(irpstart+1)=(2.0*fun(irpstart-1)+
     >      4.0*fun(irpstart)+fun(irpstart+1))*rmesh(irpstart,2)/3.0
c$$$         temp(irpstart+1)=(4.0*fun(irpstart)+fun(irpstart+1))*
c$$$     >      rmesh(irpstart,2)/3.0
      end if
      
C  Define the integral at even points of RMESH using Simpson's rule.
C  Define the integral at odd points of RMESH using a combination of
C  Simpson's rule and Newton's 3/8 rule.
      do irp=irpstart+2,irpstop-2,2
         temp(irp)=temp(irp-2)+(fun(irp-2)+4.0*fun(irp-1)+fun(irp))
     >      *rmesh(irp-1,2)/3.0
         temp(irp+1)=temp(irp-2)+(fun(irp-2)+3.0*fun(irp-1)+
     >      3.0*fun(irp)+fun(irp+1))*0.375*rmesh(irp,2)
      end do

C  The above Newton's 3/8 rule fails where dR doubles.
C  I is the next point after the doubling. It should be odd. 
      do n=2,id-1
         i=jdouble(n)+1
         if (i.lt.istart) cycle
         temp(i)=temp(i+2)-(fun(i)+4.0*fun(i+1)+fun(i+2))*rmesh(i,2)/3.0
         if (i.eq.irpstop-1) temp(i)=temp(i-1)
      end do

      if (istart.gt.1) then
C  FUN(R) from R=0 to RMESH(irpstart,1) is assumed to be zero.
         do irp=1,irpstart
            temp(irp)=0.0
         end do
      end if
C  FUN(R) from R=RMESH(irpstop,1) to 'infinity' is assumed to be zero.
      do i=irpstop,meshr
         temp(i)=temp(irpstop-1)
      end do
      end

      SUBROUTINE FUNC(N,L,UU,jstop)
      INCLUDE 'par.f'

      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
      DIMENSION U(MAXR),UU(MAXR)
      call rnl(1,n,l,uu,en,jstop)
      return
      I=10*N + L
      IF(I.EQ.10) ASSIGN 10 TO I
      IF(I.EQ.20) ASSIGN 20 TO I
      IF(I.EQ.21) ASSIGN 21 TO I
      IF(I.EQ.30) ASSIGN 30 TO I
      IF(I.EQ.31) ASSIGN 31 TO I
      IF(I.EQ.32) ASSIGN 32 TO I
      IF(I.EQ.40) ASSIGN 40 TO I
      IF(I.EQ.41) ASSIGN 41 TO I
      IF(I.EQ.42) ASSIGN 42 TO I
      IF(I.EQ.43) ASSIGN 43 TO I
      IF(I.EQ.50) ASSIGN 50 TO I
      IF(I.EQ.51) ASSIGN 51 TO I
      IF(I.EQ.52) ASSIGN 52 TO I
      IF(I.EQ.53) ASSIGN 53 TO I
      IF(I.EQ.54) ASSIGN 54 TO I
      IF(I.EQ.60) ASSIGN 60 TO I
      IF(I.EQ.61) ASSIGN 61 TO I
      IF(I.EQ.62) ASSIGN 62 TO I
      IF(I.EQ.63) ASSIGN 63 TO I
      IF(I.EQ.64) ASSIGN 64 TO I
      IF(I.EQ.65) ASSIGN 65 TO I

      rn=float(n)
      do j=1,jstop
         R=rmesh(J,1)
         t=r/rn
         X=2.0*t
         ae=exp(-t)
         GO TO I
10    CONTINUE
      U(J)= 2*AE
      GO TO 100
20    U(J)= SQRT(3.)*1./(4.89898)*(2.-X)*AE
      GO TO 100
21    U(J)= 1./(4.89898)*X*AE
      GO TO 100
30    U(J)= 1./(9.*SQRT(3.))*(6.-6.*X+X**2)*AE
      GO TO 100
31    U(J)= 1./(9.*SQRT(6.))*(4.-X)*X*AE
      GO TO 100
32    U(J)= 1./(9.*SQRT(30.))*X**2*AE
      GO TO 100
40    U(J)= 1./(96.)*(24.-36*X+12.*X*X-X**3)*AE
      GO TO 100
41    U(J)=1./(32.*SQRT(15.))*(20.-10.*X+X**2)*X*AE
      GO TO 100
42    U(J)=1./(96.*SQRT(5.))*(6.-X)*X**2*AE
      GO TO 100
43    U(J)= 1./(96.*SQRT(35.))*X**3*AE
      GO TO 100
50    U(J)= 1./(300.*SQRT(5.))*(120.-240*X+120*X**2-20*X**3+X**4)
     &    *AE
      GO TO 100
51    U(J)= 1./(150.*SQRT(30.))*(120.-90.*X+18*X**2-X**3)*X*
     &   AE
      GO TO 100
52    U(J)= 1./(150.*SQRT(70.))*(42.-14.*X+X**2)*X**2*AE
      GO TO 100
53    U(J)= 1./(300.*SQRT(70.))*(8.-X)*X**3*AE
      GO TO 100
54    U(J)= 1./(900.*SQRT(70.))*X**4*AE
      GO TO 100
60    U(J)= 1./(2160.*SQRT(6.))*(720.-1800.*X+1200.*X**2-300.*X**3
     &  +30.*X**4-X**5)*AE
      GO TO 100
61    U(J)= 1./(432.*SQRT(210.))*(840.-840.*X+252.*X**2-28.*X**3+
     &    X**4)*X*AE
      GO TO 100
62    U(J)= 1./(864.*SQRT(105.))*(336.-168.*X+24.*X**2-X**3)*X**2*
     &   AE
      GO TO 100
63    U(J)= 1./(2592.*SQRT(35.))*(72.-18.*X+X**2)*X**3*AE
      GO TO 100
64    U(J)= 1./(12960.*SQRT(7.))*(10.-X)*X**4*AE
      GO TO 100
65    U(J)= 1./(12960.*SQRT(77.))*X**5*AE
100   CONTINUE
      UU(J)=U(J)*R
      end do
      RETURN
      END

      function fact(k)
      real*8 fact
      fact=1d0
      do i=2, k
         fact=fact*float(i)
      end do
      end

C  The following routine returns the hydrogenic radial functions.
      subroutine rnl(nznuc,n,l,chi,en,jstop)
      include 'par.f'
!      include 'par.for'
      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
      common/smallr/ formcut,regcut,expcut,fast
      common /cnsts2/ factl(0:2*lcoul),api,x2(63,5),w2(63,5),res(63)
      logical fast
      dimension chi(maxr)
      real*8 fact, rlrho, rho, const2(0:2*lcoul), sumpos, sumneg,
     >   exprho,factl,api,x2,w2,res
      if (nznuc.gt.0) then
         Z = float(nznuc)
         en = - (Z / float(n)) ** 2  ! Ry
      else
C  Here we define positronium          
         Z = 0.5
         en = - 0.5 / float(n)**2  ! Ry
      endif 
      do j=1,meshr
         chi(j)=0.0
      end do

c$$$      const1=2.0/float(n*n)*sqrt(fact(n-l-1)*fact(n+l))
c$$$      do k=0,n-l-1
c$$$         const2=(-1)**k/fact(n-l-1-k)/fact(2*l+1+k)/fact(k)
c$$$         do j=1,jstop
c$$$            rho = 2.0 * Z * rmesh(j,1) / float(n)
c$$$            chi(j)=chi(j)+const1*const2*exp(-rho/2.0)*rho**(l+k)
c$$$         end do
c$$$      end do
c$$$      do j=1,jstop
c$$$         chi(j)=chi(j)*rmesh(j,1)
c$$$      end do
c$$$      const1 = sqrt(Z * fact(n-l-1) * fact(n+l)) / float(n)
      const1 = sqrt(Z) / float(n)
      rho = 0.0
      exprho = exp(- rho / 2.0)
      do k = 0, n - l - 1
         const2(k)=(factl(n-l-1)+factl(n+l))/2.0
     >      -(factl(n-l-1-k)+factl(2*l+1+k)+factl(k))
      enddo 
      j = 0
      do while (exprho.gt.expcut.and.j.lt.meshr)
         j = j + 1
c$$$      do j = 1, jstop
         rlrho = 0.0
         sumpos = 0.0
         sumneg = 0.0
         rho = 2d0 * Z * rmesh(j,1) / float(n)
         exprho = exp(- rho / 2.0)
c$$$         do k = 0, n - l - 1
c$$$c$$$            const2=(-1)**k/fact(n-l-1-k)/fact(2*l+1+k)/fact(k)
c$$$            const2=(-1)**k*exp(
c$$$     >         (factl(n-l-1)+factl(n+l))/2.0
c$$$     >         -(factl(n-l-1-k)+factl(2*l+1+k)+factl(k)))
c$$$            rlrho = rlrho + const2 * rho ** k
c$$$         enddo
         do k = 0, n - l - 1, 2
            sumpos = sumpos + exp(const2(k) + k * log(rho))
         enddo
         do k = 1, n - l - 1, 2
            sumneg = sumneg + exp(const2(k) + k * log(rho))
         enddo
         rlrho = sumpos-sumneg
         if (abs(sumpos-sumneg)/(sumpos+sumneg).lt.1e-13) rlrho=0.0
         tmp = chi(j)
c$$$         chi(j) = const1 * rlrho * rho ** (l + 1) * exprho
         chi(j) = const1 * rlrho *exp((l+1)*log(rho)-rho/2.0)
c$$$         if (abs((chi(j)-tmp)/tmp).gt.1e-4) print*,n,l,j,chi(j),tmp
      enddo
      jstop = j
      return
      end
      
      function ei(x)
      real*8 fact
      data gamma/0.5772156649/
      tmp=0.0
      n=1
      xn=x
      tmp1=1.0
      do while (abs(tmp-tmp1).gt.0.0)
         tmp1=tmp
         tmp=tmp+xn/float(n)/fact(n)
         xn=xn*x
         n=n+1
      end do
      ei=gamma+log(abs(x))+tmp
      end

      function acoth(x)
      if (abs(x).le.1.0) stop 'ACOTH defined for |ARG| > 1 only'
      acoth=log((1.0+x)/(x-1.0))/2.0
      end

      function atanh(x)
      if (abs(x).ge.1.0) stop 'ATANH defined for |ARG| < 1 only'
      atanh=log((1.0+x)/(1.0-x))/2.0
      end
      

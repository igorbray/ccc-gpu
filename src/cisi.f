C  FILE 'CISI.FTN77'
C
C  The following function returns the value of the integral from
C  A to infinity of COS(R*X)/X**N dX, if COSL is true, or
C  SIN(R*X)/X**N dX, if otherwise. RER is the relative error required.
      function csint(a,rpm,n,cosl,rer)
      implicit real*8 (a-h,o-z)
      include 'par.for'
      logical*4 cosl
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),res(63)
      if (n.le.1) then
         print*,'N is',n,' in CSINT; should be greater than 1.'
         call exit
      end if
      r=rpm
      if (r.lt.0d0) r=-r
      if (r.lt.1d-10) then
         if (cosl) then
            csint=a**(1-n)/(n-1)
         else
            csint=0d0
         end if
         return
      end if
      app=appcs(a,r,n,cosl)
      if (rer.gt.1d-1) then
         csint=app
         if (rpm.lt.0d0.and..not.cosl) csint=-csint
         return
      end if
      if (cosl) then
         fterm=ciorsi(n-1,a,r)
      else
         fterm=ciorsi(n-2,a,r)
      end if
      if (abs(fterm*1d-13/app).gt.rer) then
         csint=csass(a,r,n,cosl,rer)
         if (rpm.lt.0d0.and..not.cosl) csint=-csint
         return
      end if
      tmp=0d0
      if (cosl) then
         do 30 i=0,n-2
            temp=(-1)**i*der(i,a,r)*rint(i+1,a,n)
            tmp=tmp+temp
30       continue
         temp=(-1)**(n-1)*fterm
      else
         do 58 i=0,n-3
            temp=(-1)**i*der(i,a,r)*rint(i+1,a,n-1)
            tmp=tmp+temp
58       continue
         temp=(-1)**(n-2)*fterm
         tmp=(tmp+temp)*r/(n-1)
         temp=sin(r*a)*a**(-n+1)/(n-1)
      end if
      csint=tmp+temp
      if (rpm.lt.0d0.and..not.cosl) csint=-csint
      end

C  The following function is used when the integration by parts of eg.
C  SIN(RX)/X**N to reduce the power N, fails to give the required accuracy.
C  This function integrates by parts the above function to get a series
C  with increasing powers of N. This is the assymptotic expansion. If
C  this series diverges before RER is satisfied then the last term is
C  evaluated numerically with CSNUM.
      real*8 function csass(a,r,n,cosl,rer)
      implicit real*8 (a-h,o-z)
      include 'par.for'
      logical*4 cosl
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),res(63)
      sum=0d0
      en=1d10
      ep=1d11
      app=appcs(a,r,n,cosl)
      if (cosl) then
         i=0
 11      if (abs(en).lt.abs(ep)) then
c         while(abs(en).lt.abs(ep)) do
            ep=en
            en=(-1)**i*exp(factl(n-1+2*i)-factl(n-1)-
     >         (2*i+1)*log(r)-(n+2*i)*log(a))*
     >         (sin(r*a)-cos(r*a)*(n+2*i)/(r*a))
            sum=sum+en
            i=i+1
            if (abs(en/sum).lt.rer) then
               csass=-sum
               return
            end if
            go to 11
         end if
c         end while
         j=2*(i-2)+1
         tmp3=exp(factl(n+j)-factl(n-1)-(j+1)*log(r))
         tmp2=csnum(a,r,n+j+1,cosl,rer*abs(sum/en))*(-1)**i
      else
         i=1
 12      if (abs(en).lt.abs(ep)) then
c         while(abs(en).lt.abs(ep)) do
            ep=en
            en=(-1)**i*exp(factl(n-3+2*i)-factl(n-1)-
     >         (2*i-1)*log(r)-(n+2*i-2)*log(a))*
     >         (cos(r*a)+sin(r*a)*(n+2*i-2)/(r*a))
            sum=sum+en
            i=i+1
            if (abs(en/sum).lt.rer) then
               csass=-sum
               return
            end if
            go to 12
         end if
c        end while
         j=2*(i-2)-1
         tmp3=exp(factl(n+j)-factl(n-1)-(j+1)*log(r))
         tmp2=csnum(a,r,n+j+1,cosl,rer*abs(sum/en))*(-1)**(i-1)
      end if
      csass=-sum-tmp2*tmp3+en
      end

C  The following function approximates x**(-n) by aa*exp(al*x) and hence
C  gives a good approximation to the integral of eg. COS(r*x)*x**(-n)
C  from A to infinity.
      real*8 function appcs(a,r,n,cosl)
      implicit real*8 (a-h,o-z)
      logical*4 cosl
      aa=(exp(1d0)/a)**n
      al=-n/a
      if (cosl) then
         appcs=-aa*exp(-dble(n))*(al*cos(r*a)+r*sin(r*a))/(al*al+r*r)
      else
         appcs=-aa*exp(-dble(n))*(al*sin(r*a)-r*cos(r*a))/(al*al+r*r)
      end if
      end

C  The function DER(I,A) returns the negative of the Ith derivative
C  of COS(R*X) evaluated at X=A.
      function der(i,a,r)
      implicit real*8 (a-h,o-z)
      m=mod(i,4)
      der=-((m-3)*m+1)*r**i*(mod(i+1,2)*cos(a*r)+mod(i,2)*sin(a*r))
      end

C  The function RINT(I,A) returns the Ith derivative of X**-N (I<N),
C  evaluated at X=A.
      function rint(i,a,n)
      implicit real*8 (a-h,o-z)
      include 'par.for'
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),res(63)
      rint=(-1)**i*a**(-n+i)*exp(factl(n-1-i)-factl(n-1))
      end

C  The function CIORSI(I,A) returns the value of the integral from A to
C  infinity of the Ith derivative of COS((R1+R2)*X) times the Ith derivative
C  of X**-N. Called once only for I=N-1.
      function ciorsi(i,a,r)
      implicit real*8 (a-h,o-z)
      include 'par.for'
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),res(63)
      if (mod(i,2).eq.0) then
         ciorsi=dci(a*r)
      else
         ciorsi=dsi(a*r)
      end if
      m=mod(i,4)
      ciorsi=-((m-3)*m+1)*(-r)**i*ciorsi/exp(factl(i))
      end

C  The following function evaluates the integral of eg. SIN(RX)/X**N
C  from A to infinity numerically.
      real*8 function csnum(a,r,n,cosl,rer)
      implicit real*8 (a-h,o-z)
      logical*4 cosl
      include 'par.for'
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),result(63)
C  In the following line we multiply RER by 100 as APP is accurate to at
C  least 2 significant figures.
      rerp=rer*1d2
      m=int(a*r/pi)+1
      if (cosl) then
         sum=cosint(a,pi*m/r,r,n)
         app=appcs(a,r,n,cosl)
         i=0
 13      if (abs(app/sum).gt.rerp) then
c         while(abs(app/sum).gt.rerp) do
            as=pi*(m+i)/r
            af=as+pi/r
            sum=sum+cosint(as,af,r,n)
            app=appcs(af,r,n,cosl)
            i=i+1
            go to 13
         end if
      else
         sum=sinint(a,pi*m/r,r,n)
         app=appcs(a,r,n,cosl)
         i=0
 14      if (abs(app/sum).gt.rerp) then
c         while(abs(app/sum).gt.rerp) do
            as=pi*(m+i)/r
            af=as+pi/r
            sum=sum+sinint(as,af,r,n)
            app=appcs(af,r,n,cosl)
            i=i+1
            goto 14
         end if
      end if
      csnum=sum+app
      end

C  The following function evaluates the integral COS(RX)/X**N from A to B to
C  near machine precision.
      real*8 function cosint(a,b,r,n)
      implicit real*8 (a-h,o-z)
      dimension i2(0:5)
      include 'par.for'
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),result(63)
      data i2/0,3,7,15,31,63/
      tp=1d9
      do 100 j=1,5
         tn=0d0
         do 80 i=i2(j-1)+1,i2(j)
            x=((b-a)*x2(i,j)+b+a)/2
            result(i)=cos(r*x)/x**n*(b-a)/2
            tn=tn+w2(i,j)*result(i)
80       continue
         do 90 i=1,i2(j-1)
            tn=tn+w2(i,j)*result(i)
90       continue
         if (abs(1d0-tn/tp).lt.1d-7) go to 110
         tp=tn
100   continue
110   cosint=tn
      end

C  The following function evaluates the integral SIN(RX)/X**N from A to B to
C  near machine precision.
      real*8 function sinint(a,b,r,n)
      implicit real*8 (a-h,o-z)
      dimension i2(0:5)
      include 'par.for'
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),result(63)
      data i2/0,3,7,15,31,63/
      tp=1d9
      do 100 j=1,5
         tn=0d0
         do 80 i=i2(j-1)+1,i2(j)
            x=((b-a)*x2(i,j)+b+a)/2
            result(i)=sin(r*x)/x**n*(b-a)/2
            tn=tn+w2(i,j)*result(i)
80       continue
         do 90 i=1,i2(j-1)
            tn=tn+w2(i,j)*result(i)
90       continue
         if (abs(1d0-tn/tp).lt.1d-7) go to 110
         tp=tn
100   continue
110   sinint=tn
      end

      
Compliments of netlib   Tue Nov  4 20:43:48 CST 1986
      double precision function dsi (x)
c december 1980 edition, w. fullerton, bell labs.
      double precision x, sics(18), pi2, xsml, absx, f, g, cosx,
     1  dcsevl, d1mach, cos, sin, dsqrt
      external d1mach, dcsevl, initds
c
c series for si   on the interval  0.00000e+00 to  1.60000e+01
c                                        with weighted error   8.58e-32
c                                         log weighted error  31.07
c                               significant figures required  30.53
c                                    decimal places required  31.69
c
      data si  cs(  1) / -0.1315646598 1848419289 0427517300 0457d0/
      data si  cs(  2) / -0.2776578526 9736018920 4828766015 7299d0/
      data si  cs(  3) /  0.0354414054 8666591797 4913546471 0086d0/
      data si  cs(  4) / -0.0025631631 4479339776 5875278836 1530d0/
      data si  cs(  5) /  0.0001162365 3904970092 8126492148 2985d0/
      data si  cs(  6) / -0.0000035904 3272416060 4267000434 7148d0/
      data si  cs(  7) /  0.0000000802 3421237057 1016230865 2976d0/
      data si  cs(  8) / -0.0000000013 5629976925 4025064993 1846d0/
      data si  cs(  9) /  0.0000000000 1794407215 9973677556 7759d0/
      data si  cs( 10) / -0.0000000000 0019083873 4308714549 0737d0/
      data si  cs( 11) /  0.0000000000 0000166699 8958682433 0853d0/
      data si  cs( 12) / -0.0000000000 0000001217 3098836850 3042d0/
      data si  cs( 13) /  0.0000000000 0000000007 5418186699 3865d0/
      data si  cs( 14) / -0.0000000000 0000000000 0401417884 2446d0/
      data si  cs( 15) /  0.0000000000 0000000000 0001855369 0716d0/
      data si  cs( 16) / -0.0000000000 0000000000 0000007516 6966d0/
      data si  cs( 17) /  0.0000000000 0000000000 0000000026 9113d0/
      data si  cs( 18) / -0.0000000000 0000000000 0000000000 0858d0/
c
      data pi2 / 1.5707963267 9489661923 1321691639 75 d0 /
      data nsi, xsml /0, 0.0d0/
c
      if (nsi.ne.0) go to 10
c$$$      print*,'d1mach(3):',d1mach(3)
      nsi = initds (sics, 18, 1e-12) !0.1*sngl(d1mach(3)))
      xsml = dsqrt(d1mach(3))
c
 10   absx = dabs(x)
      if (absx.gt.4.0d0) go to 20
      dsi = x-pi2
      if (absx.lt.xsml) return
c
      dsi = x*(0.75d0 + dcsevl ((x*x-8.0d0)*.125d0, sics, nsi))-pi2
      return
c
 20   call d9sifg (absx, f, g)
      cosx = cos (absx)
      call erroff
c      dsi = pi2 - f*cosx - g*sin(x)
      dsi = - f*cosx - g*sin(x)
      if (x.lt.0.0d0) dsi = -dsi
c
      return
      end

      subroutine erroff
c
c  turns off the error state off by setting lerror=0.
c
      external i8save
c
      i=i8save(1,0,.true.)
      return
c
      end
      integer function i8save(isw,ivalue,set)
c
c  if (isw = 1) i8save returns the current error number and
c               sets it to ivalue if set = .true. .
c
c  if (isw = 2) i8save returns the current recovery switch and
c               sets it to ivalue if set = .true. .
c
      logical set
c
      integer iparam(2)
c  iparam(1) is the error number and iparam(2) is the recovery switch.
c
c  start execution error free and with recovery turned off.
c
      data iparam(1) /0/,  iparam(2) /2/
c
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
c
      return
c
      end
      function initds (dos, nos, eta)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
c
c initialize the double precision orthogonal series dos so that initds
c is the number of terms needed to insure the error is no larger than
c eta.  ordinarily eta will be chosen to be one-tenth machine precision.
c
c             input arguments --
c dos    dble prec array of nos coefficients in an orthogonal series.
c nos    number of coefficients in dos.
c eta    requested accuracy of series.
c
      double precision dos(nos)
c
      if (nos.lt.1) call seteru (
     1  'initds  number of coefficients lt 1', 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call seteru ('initds  eta may be too small', 28,
     1  1, 2)
      initds = i
c
      return
      end
      subroutine seteru (messg, nmessg, nerr, iopt)
      common /cseter/ iunflo
      integer messg(1)
      iunflo= 0
c
      if (iopt.ne.0) call seterr (messg, nmessg, nerr, iopt)
      if (iopt.ne.0) return
c
      if (iunflo.le.0) return
      call seterr (messg, nmessg, nerr, 1)
c
      return
      end

      function csevl (x, cs, n)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate the n-term chebyshev series cs at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
c and parker, chebyshev polys in numerical analysis, oxford press, p.56.
c
c             input arguments --
c x      value at which the series is to be evaluated.
c cs     array of n terms of a chebyshev series.  in eval-
c        uating cs, only half the first coef is summed.
c n      number of terms in array cs.
c
      dimension cs(1)
c
      if (n.lt.1) call seteru ('csevl   number of terms le 0', 28, 2,2)
      if (n.gt.1000) call seteru ('csevl   number of terms gt 1000',
     1  31, 3, 2)
      if (x.lt.(-1.1) .or. x.gt.1.1) call seteru (
     1  'csevl   x outside (-1,+1)', 25, 1, 1)
c
      b1 = 0.
      b0 = 0.
      twox = 2.*x
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
 10   continue
c
      csevl = 0.5 * (b0-b2)
c
      return
      end
      double precision function d9pak (y, n)
c december 1979 edition. w. fullerton, c3, los alamos scientific lab.
c
c pack a base 2 exponent into floating point number x.  this routine is
c almost the inverse of d9upak.  it is not exactly the inverse, because
c dabs(x) need not be between 0.5 and 1.0.  if both d9pak and 2.d0**n
c were known to be in range we could compute
c                d9pak = x * 2.0d0**n
c
      double precision y, aln2b, aln210, d1mach
      external d1mach, i1mach
      data nmin, nmax / 2*0 /
      data aln210 / 3.321928094 8873623478 7031942948 9 d0 /
c
      if (nmin.ne.0) go to 10
      aln2b = 1.0d0
      if (i1mach(10).ne.2) aln2b = d1mach(5)*aln210
      nmin = aln2b*dble(float(i1mach(15)))
      nmax = aln2b*dble(float(i1mach(16)))
c
 10   call d9upak (y, d9pak, ny)
c
      nsum = n + ny
      if (nsum.lt.nmin) go to 40
      if (nsum.gt.nmax) call seteru (
     1  'd9pak   packed number overflows', 31, 1, 2)
c
      if (nsum.eq.0) return
      if (nsum.gt.0) go to 30
c
 20   d9pak = 0.5d0*d9pak
      nsum = nsum + 1
      if (nsum.ne.0) go to 20
      return
c
 30   d9pak = 2.0d0*d9pak
      nsum = nsum - 1
      if (nsum.ne.0) go to 30
      return
c
 40   call seteru ('d9pak   packed number underflows', 32, 1, 0)
      d9pak = 0.0d0
      return
c
      end
      subroutine d9sifg (x, f, g)
c december 1980 edition.  w. fullerton, bell labs.
      double precision x, f, g, f1cs(43), f2cs(99), g1cs(44),
     1  g2cs(44), g3cs(56), xbnd, xbndg, xbig, xmaxf, xmaxg,
     2  dcsevl, d1mach, exp, log, dsqrt
      external d1mach, dcsevl, initds
c
c series for f1   on the interval  2.00000e-02 to  6.25000e-02
c                                        with weighted error   2.45e-32
c                                         log weighted error  31.61
c                               significant figures required  30.42
c                                    decimal places required  32.43
c
      data f1  cs(  1) / -0.1191081969 0513636103 4820196582 8918d0/
      data f1  cs(  2) / -0.0247823144 9962362475 9007415082 3133d0/
      data f1  cs(  3) /  0.0011910281 4533578212 6812036305 4457d0/
      data f1  cs(  4) / -0.0000927027 7143885617 4830860036 0706d0/
      data f1  cs(  5) /  0.0000093373 1415682709 9686820458 2766d0/
      data f1  cs(  6) / -0.0000011058 2878205571 4393897942 6306d0/
      data f1  cs(  7) /  0.0000001464 7720714601 6216933655 0799d0/
      data f1  cs(  8) / -0.0000000210 6944962876 8953260122 7548d0/
      data f1  cs(  9) /  0.0000000032 2934923668 4823638285 7374d0/
      data f1  cs( 10) / -0.0000000005 2065296175 2937582801 4986d0/
      data f1  cs( 11) /  0.0000000000 8748788845 7027875026 8316d0/
      data f1  cs( 12) / -0.0000000000 1521761870 5612366829 4574d0/
      data f1  cs( 13) /  0.0000000000 0272571924 0541957390 0583d0/
      data f1  cs( 14) / -0.0000000000 0050070530 7596855629 0255d0/
      data f1  cs( 15) /  0.0000000000 0009402409 0272606851 1779d0/
      data f1  cs( 16) / -0.0000000000 0001800144 4479180367 8336d0/
      data f1  cs( 17) /  0.0000000000 0000350626 2143274178 5826d0/
      data f1  cs( 18) / -0.0000000000 0000069352 8292676914 9709d0/
      data f1  cs( 19) /  0.0000000000 0000013909 2513645421 6568d0/
      data f1  cs( 20) / -0.0000000000 0000002824 8688507417 0585d0/
      data f1  cs( 21) /  0.0000000000 0000000580 3130569357 9081d0/
      data f1  cs( 22) / -0.0000000000 0000000120 4690157337 5820d0/
      data f1  cs( 23) /  0.0000000000 0000000025 2505244365 5940d0/
      data f1  cs( 24) / -0.0000000000 0000000005 3398026880 5594d0/
      data f1  cs( 25) /  0.0000000000 0000000001 1385578627 4122d0/
      data f1  cs( 26) / -0.0000000000 0000000000 2446286150 5259d0/
      data f1  cs( 27) /  0.0000000000 0000000000 0529365932 0439d0/
      data f1  cs( 28) / -0.0000000000 0000000000 0115318494 0277d0/
      data f1  cs( 29) /  0.0000000000 0000000000 0025278656 8318d0/
      data f1  cs( 30) / -0.0000000000 0000000000 0005573864 5378d0/
      data f1  cs( 31) /  0.0000000000 0000000000 0001235824 5621d0/
      data f1  cs( 32) / -0.0000000000 0000000000 0000275435 0842d0/
      data f1  cs( 33) /  0.0000000000 0000000000 0000061690 6808d0/
      data f1  cs( 34) / -0.0000000000 0000000000 0000013881 7443d0/
      data f1  cs( 35) /  0.0000000000 0000000000 0000003137 5329d0/
      data f1  cs( 36) / -0.0000000000 0000000000 0000000712 1249d0/
      data f1  cs( 37) /  0.0000000000 0000000000 0000000162 2778d0/
      data f1  cs( 38) / -0.0000000000 0000000000 0000000037 1206d0/
      data f1  cs( 39) /  0.0000000000 0000000000 0000000008 5221d0/
      data f1  cs( 40) / -0.0000000000 0000000000 0000000001 9633d0/
      data f1  cs( 41) /  0.0000000000 0000000000 0000000000 4538d0/
      data f1  cs( 42) / -0.0000000000 0000000000 0000000000 1052d0/
      data f1  cs( 43) /  0.0000000000 0000000000 0000000000 0245d0/
c
c series for f2   on the interval  0.00000e+00 to  2.00000e-02
c                                        with weighted error   2.38e-32
c                                         log weighted error  31.62
c                               significant figures required  30.01
c                                    decimal places required  32.62
c
      data f2  cs(  1) / -0.0348409253 8970132330 8360497337 45577d0/
      data f2  cs(  2) / -0.0166842205 6779596873 2467863122 78676d0/
      data f2  cs(  3) /  0.0006752901 2412377385 0452078592 39727d0/
      data f2  cs(  4) / -0.0000535066 6225447013 6287855775 57429d0/
      data f2  cs(  5) /  0.0000062693 4217790075 2670507594 31626d0/
      data f2  cs(  6) / -0.0000009526 6388019916 6806777904 14293d0/
      data f2  cs(  7) /  0.0000001745 6292242509 8804255044 27666d0/
      data f2  cs(  8) / -0.0000000368 7954030653 0933070976 46628d0/
      data f2  cs(  9) /  0.0000000087 2026777051 3952640758 16938d0/
      data f2  cs( 10) / -0.0000000022 6019703919 7387485304 23167d0/
      data f2  cs( 11) /  0.0000000006 3246249765 2506125204 44877d0/
      data f2  cs( 12) / -0.0000000001 8889118884 7178692409 11480d0/
      data f2  cs( 13) /  0.0000000000 5967746729 9978133726 20472d0/
      data f2  cs( 14) / -0.0000000000 1980443117 3722390111 96007d0/
      data f2  cs( 15) /  0.0000000000 0686413954 7721033837 13264d0/
      data f2  cs( 16) / -0.0000000000 0247310193 0701991060 74890d0/
      data f2  cs( 17) /  0.0000000000 0092263594 5499414041 96042d0/
      data f2  cs( 18) / -0.0000000000 0035523634 9992617844 97297d0/
      data f2  cs( 19) /  0.0000000000 0014076049 6253515914 61820d0/
      data f2  cs( 20) / -0.0000000000 0005726228 4997476527 94311d0/
      data f2  cs( 21) /  0.0000000000 0002386537 5454131718 10106d0/
      data f2  cs( 22) / -0.0000000000 0001017141 8907645971 42232d0/
      data f2  cs( 23) /  0.0000000000 0000442594 5310783644 24968d0/
      data f2  cs( 24) / -0.0000000000 0000196344 9330491897 61979d0/
      data f2  cs( 25) /  0.0000000000 0000088688 7483148104 61024d0/
      data f2  cs( 26) / -0.0000000000 0000040743 3450273115 46948d0/
      data f2  cs( 27) /  0.0000000000 0000019016 8372156753 39859d0/
      data f2  cs( 28) / -0.0000000000 0000009009 7072974780 42442d0/
      data f2  cs( 29) /  0.0000000000 0000004329 2112740956 68667d0/
      data f2  cs( 30) / -0.0000000000 0000002108 1444653224 79526d0/
      data f2  cs( 31) /  0.0000000000 0000001039 6379070264 52274d0/
      data f2  cs( 32) / -0.0000000000 0000000518 8910079489 31936d0/
      data f2  cs( 33) /  0.0000000000 0000000261 9553248698 99371d0/
      data f2  cs( 34) / -0.0000000000 0000000133 6903999513 01570d0/
      data f2  cs( 35) /  0.0000000000 0000000068 9410577029 31664d0/
      data f2  cs( 36) / -0.0000000000 0000000035 9053626104 37250d0/
      data f2  cs( 37) /  0.0000000000 0000000018 8780772557 91706d0/
      data f2  cs( 38) / -0.0000000000 0000000010 0161252655 94380d0/
      data f2  cs( 39) /  0.0000000000 0000000005 3607256915 78228d0/
      data f2  cs( 40) / -0.0000000000 0000000002 8931989749 44827d0/
      data f2  cs( 41) /  0.0000000000 0000000001 5740651002 02625d0/
      data f2  cs( 42) / -0.0000000000 0000000000 8630271064 31206d0/
      data f2  cs( 43) /  0.0000000000 0000000000 4767156028 62288d0/
      data f2  cs( 44) / -0.0000000000 0000000000 2652227399 98504d0/
      data f2  cs( 45) /  0.0000000000 0000000000 1485828650 63866d0/
      data f2  cs( 46) / -0.0000000000 0000000000 0837972359 23135d0/
      data f2  cs( 47) /  0.0000000000 0000000000 0475659164 22711d0/
      data f2  cs( 48) / -0.0000000000 0000000000 0271690733 53112d0/
      data f2  cs( 49) /  0.0000000000 0000000000 0156127388 81686d0/
      data f2  cs( 50) / -0.0000000000 0000000000 0090245550 78347d0/
      data f2  cs( 51) /  0.0000000000 0000000000 0052460970 49119d0/
      data f2  cs( 52) / -0.0000000000 0000000000 0030664508 18697d0/
      data f2  cs( 53) /  0.0000000000 0000000000 0018019962 50957d0/
      data f2  cs( 54) / -0.0000000000 0000000000 0010644430 50752d0/
      data f2  cs( 55) /  0.0000000000 0000000000 0006319421 58881d0/
      data f2  cs( 56) / -0.0000000000 0000000000 0003770138 12246d0/
      data f2  cs( 57) /  0.0000000000 0000000000 0002259975 42918d0/
      data f2  cs( 58) / -0.0000000000 0000000000 0001361008 44814d0/
      data f2  cs( 59) /  0.0000000000 0000000000 0000823332 32003d0/
      data f2  cs( 60) / -0.0000000000 0000000000 0000500259 86091d0/
      data f2  cs( 61) /  0.0000000000 0000000000 0000305262 45684d0/
      data f2  cs( 62) / -0.0000000000 0000000000 0000187051 64021d0/
      data f2  cs( 63) /  0.0000000000 0000000000 0000115084 04393d0/
      data f2  cs( 64) / -0.0000000000 0000000000 0000071087 14611d0/
      data f2  cs( 65) /  0.0000000000 0000000000 0000044080 65533d0/
      data f2  cs( 66) / -0.0000000000 0000000000 0000027437 60867d0/
      data f2  cs( 67) /  0.0000000000 0000000000 0000017141 44851d0/
      data f2  cs( 68) / -0.0000000000 0000000000 0000010747 68860d0/
      data f2  cs( 69) /  0.0000000000 0000000000 0000006762 59777d0/
      data f2  cs( 70) / -0.0000000000 0000000000 0000004269 81348d0/
      data f2  cs( 71) /  0.0000000000 0000000000 0000002705 00637d0/
      data f2  cs( 72) / -0.0000000000 0000000000 0000001719 33331d0/
      data f2  cs( 73) /  0.0000000000 0000000000 0000001096 36138d0/
      data f2  cs( 74) / -0.0000000000 0000000000 0000000701 32573d0/
      data f2  cs( 75) /  0.0000000000 0000000000 0000000450 01784d0/
      data f2  cs( 76) / -0.0000000000 0000000000 0000000289 63835d0/
      data f2  cs( 77) /  0.0000000000 0000000000 0000000186 97009d0/
      data f2  cs( 78) / -0.0000000000 0000000000 0000000121 04646d0/
      data f2  cs( 79) /  0.0000000000 0000000000 0000000078 59065d0/
      data f2  cs( 80) / -0.0000000000 0000000000 0000000051 16867d0/
      data f2  cs( 81) /  0.0000000000 0000000000 0000000033 40627d0/
      data f2  cs( 82) / -0.0000000000 0000000000 0000000021 86851d0/
      data f2  cs( 83) /  0.0000000000 0000000000 0000000014 35340d0/
      data f2  cs( 84) / -0.0000000000 0000000000 0000000009 44523d0/
      data f2  cs( 85) /  0.0000000000 0000000000 0000000006 23117d0/
      data f2  cs( 86) / -0.0000000000 0000000000 0000000004 12101d0/
      data f2  cs( 87) /  0.0000000000 0000000000 0000000002 73208d0/
      data f2  cs( 88) / -0.0000000000 0000000000 0000000001 81558d0/
      data f2  cs( 89) /  0.0000000000 0000000000 0000000001 20934d0/
      data f2  cs( 90) / -0.0000000000 0000000000 0000000000 80737d0/
      data f2  cs( 91) /  0.0000000000 0000000000 0000000000 54022d0/
      data f2  cs( 92) / -0.0000000000 0000000000 0000000000 36227d0/
      data f2  cs( 93) /  0.0000000000 0000000000 0000000000 24348d0/
      data f2  cs( 94) / -0.0000000000 0000000000 0000000000 16401d0/
      data f2  cs( 95) /  0.0000000000 0000000000 0000000000 11074d0/
      data f2  cs( 96) / -0.0000000000 0000000000 0000000000 07497d0/
      data f2  cs( 97) /  0.0000000000 0000000000 0000000000 05091d0/
      data f2  cs( 98) / -0.0000000000 0000000000 0000000000 03470d0/
      data f2  cs( 99) /  0.0000000000 0000000000 0000000000 02377d0/
c
c series for g1   on the interval  2.00000e-02 to  6.25000e-02
c                                        with weighted error   7.23e-32
c                                         log weighted error  31.14
c                               significant figures required  30.35
c                                    decimal places required  31.96
c
      data g1  cs(  1) / -0.3040578798 2534959544 9972668209 1083d0/
      data g1  cs(  2) / -0.0566890984 5971205877 3133915611 8269d0/
      data g1  cs(  3) /  0.0039046158 1732756439 1998407155 4082d0/
      data g1  cs(  4) / -0.0003746075 9592022606 1861933986 7489d0/
      data g1  cs(  5) /  0.0000435431 5565598436 7955222084 0065d0/
      data g1  cs(  6) / -0.0000057417 2944530250 4656197072 3475d0/
      data g1  cs(  7) /  0.0000008282 5521045026 2974193761 6492d0/
      data g1  cs(  8) / -0.0000001278 2458925946 4272788391 3223d0/
      data g1  cs(  9) /  0.0000000207 9783529486 8788443925 7529d0/
      data g1  cs( 10) / -0.0000000035 3132059219 9079804203 2682d0/
      data g1  cs( 11) /  0.0000000006 2108242363 0895106863 1449d0/
      data g1  cs( 12) / -0.0000000001 1252154744 4629264933 6987d0/
      data g1  cs( 13) /  0.0000000000 2090889176 8442160526 7019d0/
      data g1  cs( 14) / -0.0000000000 0397158317 3768172768 9158d0/
      data g1  cs( 15) /  0.0000000000 0076904313 1427208993 9005d0/
      data g1  cs( 16) / -0.0000000000 0015146967 4273161351 9826d0/
      data g1  cs( 17) /  0.0000000000 0003028921 4655235968 4119d0/
      data g1  cs( 18) / -0.0000000000 0000613997 0383470882 5400d0/
      data g1  cs( 19) /  0.0000000000 0000126006 0582951093 3553d0/
      data g1  cs( 20) / -0.0000000000 0000026150 2925093948 3683d0/
      data g1  cs( 21) /  0.0000000000 0000005482 7884489179 6821d0/
      data g1  cs( 22) / -0.0000000000 0000001160 3818212952 6571d0/
      data g1  cs( 23) /  0.0000000000 0000000247 7165410712 9795d0/
      data g1  cs( 24) / -0.0000000000 0000000053 3067275322 3389d0/
      data g1  cs( 25) /  0.0000000000 0000000011 5566607559 8465d0/
      data g1  cs( 26) / -0.0000000000 0000000002 5228054774 4957d0/
      data g1  cs( 27) /  0.0000000000 0000000000 5542903855 0786d0/
      data g1  cs( 28) / -0.0000000000 0000000000 1225220842 1297d0/
      data g1  cs( 29) /  0.0000000000 0000000000 0272366431 8684d0/
      data g1  cs( 30) / -0.0000000000 0000000000 0060870783 1422d0/
      data g1  cs( 31) /  0.0000000000 0000000000 0013672487 4476d0/
      data g1  cs( 32) / -0.0000000000 0000000000 0003085662 6806d0/
      data g1  cs( 33) /  0.0000000000 0000000000 0000699521 2319d0/
      data g1  cs( 34) / -0.0000000000 0000000000 0000159258 7569d0/
      data g1  cs( 35) /  0.0000000000 0000000000 0000036405 1056d0/
      data g1  cs( 36) / -0.0000000000 0000000000 0000008353 9465d0/
      data g1  cs( 37) /  0.0000000000 0000000000 0000001924 0303d0/
      data g1  cs( 38) / -0.0000000000 0000000000 0000000444 6816d0/
      data g1  cs( 39) /  0.0000000000 0000000000 0000000103 1182d0/
      data g1  cs( 40) / -0.0000000000 0000000000 0000000023 9887d0/
      data g1  cs( 41) /  0.0000000000 0000000000 0000000005 5976d0/
      data g1  cs( 42) / -0.0000000000 0000000000 0000000001 3100d0/
      data g1  cs( 43) /  0.0000000000 0000000000 0000000000 3074d0/
      data g1  cs( 44) / -0.0000000000 0000000000 0000000000 0723d0/
c
c series for g2   on the interval  5.00000e-03 to  2.00000e-02
c                                        with weighted error   3.25e-32
c                                         log weighted error  31.49
c                               significant figures required  30.32
c                                    decimal places required  32.31
c
      data g2  cs(  1) / -0.1211802894 7316462635 4183404685 8267d0/
      data g2  cs(  2) / -0.0316761386 3949502867 0140792350 5610d0/
      data g2  cs(  3) /  0.0013383199 7788626801 6381942949 2182d0/
      data g2  cs(  4) / -0.0000895511 0113922524 2553190506 9518d0/
      data g2  cs(  5) /  0.0000079155 5629617182 1311524946 7924d0/
      data g2  cs(  6) / -0.0000008438 7933222415 2018141898 2080d0/
      data g2  cs(  7) /  0.0000001029 9804256775 3014664722 7274d0/
      data g2  cs(  8) / -0.0000000139 2957506051 8383579583 4444d0/
      data g2  cs(  9) /  0.0000000020 4227039598 7598040067 7594d0/
      data g2  cs( 10) / -0.0000000003 1965346942 0642703543 4752d0/
      data g2  cs( 11) /  0.0000000000 5281478326 5726769861 5312d0/
      data g2  cs( 12) / -0.0000000000 0913395546 7267103373 5289d0/
      data g2  cs( 13) /  0.0000000000 0164262512 3896776044 4819d0/
      data g2  cs( 14) / -0.0000000000 0030558970 3932266000 2410d0/
      data g2  cs( 15) /  0.0000000000 0005856558 2578577971 7892d0/
      data g2  cs( 16) / -0.0000000000 0001152291 9773094012 0563d0/
      data g2  cs( 17) /  0.0000000000 0000232094 6911998853 7310d0/
      data g2  cs( 18) / -0.0000000000 0000047743 5583417753 5025d0/
      data g2  cs( 19) /  0.0000000000 0000010009 9676580018 0573d0/
      data g2  cs( 20) / -0.0000000000 0000002135 3377808225 6704d0/
      data g2  cs( 21) /  0.0000000000 0000000462 7719077736 7671d0/
      data g2  cs( 22) / -0.0000000000 0000000101 7580741022 7657d0/
      data g2  cs( 23) /  0.0000000000 0000000022 6765739988 4672d0/
      data g2  cs( 24) / -0.0000000000 0000000005 1163077607 6426d0/
      data g2  cs( 25) /  0.0000000000 0000000001 1676701491 3108d0/
      data g2  cs( 26) / -0.0000000000 0000000000 2693542767 2470d0/
      data g2  cs( 27) /  0.0000000000 0000000000 0627566584 1146d0/
      data g2  cs( 28) / -0.0000000000 0000000000 0147588055 7531d0/
      data g2  cs( 29) /  0.0000000000 0000000000 0035014531 4739d0/
      data g2  cs( 30) / -0.0000000000 0000000000 0008375773 2152d0/
      data g2  cs( 31) /  0.0000000000 0000000000 0002019181 5152d0/
      data g2  cs( 32) / -0.0000000000 0000000000 0000490356 7705d0/
      data g2  cs( 33) /  0.0000000000 0000000000 0000119912 3348d0/
      data g2  cs( 34) / -0.0000000000 0000000000 0000029517 0610d0/
      data g2  cs( 35) /  0.0000000000 0000000000 0000007311 3112d0/
      data g2  cs( 36) / -0.0000000000 0000000000 0000001821 7843d0/
      data g2  cs( 37) /  0.0000000000 0000000000 0000000456 5148d0/
      data g2  cs( 38) / -0.0000000000 0000000000 0000000115 0151d0/
      data g2  cs( 39) /  0.0000000000 0000000000 0000000029 1267d0/
      data g2  cs( 40) / -0.0000000000 0000000000 0000000007 4125d0/
      data g2  cs( 41) /  0.0000000000 0000000000 0000000001 8953d0/
      data g2  cs( 42) / -0.0000000000 0000000000 0000000000 4868d0/
      data g2  cs( 43) /  0.0000000000 0000000000 0000000000 1256d0/
      data g2  cs( 44) / -0.0000000000 0000000000 0000000000 0325d0/
c
c series for g3   on the interval  0.00000e+00 to  5.00000e-03
c                                        with weighted error   3.83e-32
c                                         log weighted error  31.42
c                               significant figures required  29.71
c                                    decimal places required  32.29
c
      data g3  cs(  1) / -0.0280574367 8094729284 0281526433 5299d0/
      data g3  cs(  2) / -0.0137271597 1622369754 0910050808 9556d0/
      data g3  cs(  3) /  0.0002894032 6387602960 2744894127 3751d0/
      data g3  cs(  4) / -0.0000114129 2393911971 4590874362 2517d0/
      data g3  cs(  5) /  0.0000006813 9655907262 4299772020 7302d0/
      data g3  cs(  6) / -0.0000000547 9522896046 5236366905 8052d0/
      data g3  cs(  7) /  0.0000000055 2074299182 1252910940 6521d0/
      data g3  cs(  8) / -0.0000000006 6414641993 2292002249 1428d0/
      data g3  cs(  9) /  0.0000000000 9223736634 8704110856 4960d0/
      data g3  cs( 10) / -0.0000000000 1442990888 8668286261 1718d0/
      data g3  cs( 11) /  0.0000000000 0249639048 9203071024 8705d0/
      data g3  cs( 12) / -0.0000000000 0047082406 7587524472 2971d0/
      data g3  cs( 13) /  0.0000000000 0009572176 5921675998 8140d0/
      data g3  cs( 14) / -0.0000000000 0002078899 6609580903 0537d0/
      data g3  cs( 15) /  0.0000000000 0000478750 9997087743 1627d0/
      data g3  cs( 16) / -0.0000000000 0000116190 7058337717 3759d0/
      data g3  cs( 17) /  0.0000000000 0000029565 0896926783 6974d0/
      data g3  cs( 18) / -0.0000000000 0000007852 9498825649 2025d0/
      data g3  cs( 19) /  0.0000000000 0000002169 2226436825 6612d0/
      data g3  cs( 20) / -0.0000000000 0000000621 1351583167 6342d0/
      data g3  cs( 21) /  0.0000000000 0000000183 8456883845 0977d0/
      data g3  cs( 22) / -0.0000000000 0000000056 1088748213 7276d0/
      data g3  cs( 23) /  0.0000000000 0000000017 6186280528 0062d0/
      data g3  cs( 24) / -0.0000000000 0000000005 6811105054 1451d0/
      data g3  cs( 25) /  0.0000000000 0000000001 8778627958 2313d0/
      data g3  cs( 26) / -0.0000000000 0000000000 6353169415 1124d0/
      data g3  cs( 27) /  0.0000000000 0000000000 2196880236 8238d0/
      data g3  cs( 28) / -0.0000000000 0000000000 0775466655 0395d0/
      data g3  cs( 29) /  0.0000000000 0000000000 0279101835 6581d0/
      data g3  cs( 30) / -0.0000000000 0000000000 0102317852 5247d0/
      data g3  cs( 31) /  0.0000000000 0000000000 0038169340 3919d0/
      data g3  cs( 32) / -0.0000000000 0000000000 0014476789 5606d0/
      data g3  cs( 33) /  0.0000000000 0000000000 0005577951 2634d0/
      data g3  cs( 34) / -0.0000000000 0000000000 0002181723 9071d0/
      data g3  cs( 35) /  0.0000000000 0000000000 0000865664 6309d0/
      data g3  cs( 36) / -0.0000000000 0000000000 0000348215 7895d0/
      data g3  cs( 37) /  0.0000000000 0000000000 0000141918 8130d0/
      data g3  cs( 38) / -0.0000000000 0000000000 0000058571 4314d0/
      data g3  cs( 39) /  0.0000000000 0000000000 0000024466 0482d0/
      data g3  cs( 40) / -0.0000000000 0000000000 0000010338 7099d0/
      data g3  cs( 41) /  0.0000000000 0000000000 0000004417 7299d0/
      data g3  cs( 42) / -0.0000000000 0000000000 0000001908 0079d0/
      data g3  cs( 43) /  0.0000000000 0000000000 0000000832 6038d0/
      data g3  cs( 44) / -0.0000000000 0000000000 0000000366 9553d0/
      data g3  cs( 45) /  0.0000000000 0000000000 0000000163 2875d0/
      data g3  cs( 46) / -0.0000000000 0000000000 0000000073 3357d0/
      data g3  cs( 47) /  0.0000000000 0000000000 0000000033 2327d0/
      data g3  cs( 48) / -0.0000000000 0000000000 0000000015 1906d0/
      data g3  cs( 49) /  0.0000000000 0000000000 0000000007 0020d0/
      data g3  cs( 50) / -0.0000000000 0000000000 0000000003 2539d0/
      data g3  cs( 51) /  0.0000000000 0000000000 0000000001 5240d0/
      data g3  cs( 52) / -0.0000000000 0000000000 0000000000 7193d0/
      data g3  cs( 53) /  0.0000000000 0000000000 0000000000 3420d0/
      data g3  cs( 54) / -0.0000000000 0000000000 0000000000 1638d0/
      data g3  cs( 55) /  0.0000000000 0000000000 0000000000 0790d0/
      data g3  cs( 56) / -0.0000000000 0000000000 0000000000 0383d0/
c
      data nf1, nf2, ng1, ng2, ng3 / 5*0 /
      data xbnd, xbndg, xbig, xmaxf, xmaxg / 5*0.0d0 /
c
      if (nf1.ne.0) go to 10
      tol = 1e-12 !0.1d0*d1mach(3)
      nf1 = initds (f1cs, 43, tol)
      nf2 = initds (f2cs, 99, tol)
      ng1 = initds (g1cs, 44, tol)
      ng2 = initds (g2cs, 44, tol)
      ng3 = initds (g3cs, 56, tol)
c
      xbig = dsqrt(1.0d0/d1mach(3))
      xmaxf = exp (dmin1 (-log(d1mach(1)), log(d1mach(2))) - .01d0)
      xmaxg = 1.0d0/dsqrt(d1mach(1))
      xbnd = dsqrt (50.0d0)
      xbndg = dsqrt (200.d0)
c
 10   if (x.lt.4.0d0) call seteru (
     1  'd9sifg  approxs invalid for x lt 4', 34, 1, 2)
c
      if (x.gt.xbnd) go to 20
      f = (1.0d0 + dcsevl ((1.d0/x**2-0.04125d0)/.02125d0, f1cs, nf1))/x
      g = (1.0d0 + dcsevl((1.d0/x**2-.04125d0)/.02125d0, g1cs,ng1))/x**2
      return
c
 20   if (x.gt.xbig) go to 30
      f = (1.0d0 + dcsevl (100.d0/x**2-1.d0, f2cs, nf2))/x
      if (x.le.xbndg) g = (1.0d0 + dcsevl ((10000.d0/x**2-125.d0)/75.d0,
     1  g2cs, ng2))/x**2
      if (x.gt.xbndg) g = (1.0d0 + dcsevl (400.d0/x**2-1.d0, g3cs,
     1  ng3))/x**2
      return
c
 30   f = 0.d0
      if (x.lt.xmaxf) f = 1.0d0/x
      g = 0.d0
      if (x.lt.xmaxg) g = 1.0d0/x**2
      return
c
      end
      subroutine d9upak (x, y, n)
c august 1980 portable edition.  w. fullerton, los alamos scientific lab
c
c unpack floating point number x so that x = y * 2.0**n, where
c 0.5 .le. abs(y) .lt. 1.0 .
c
      double precision x, y, absx
c
      absx = dabs(x)
      n = 0
      y = 0.0d0
      if (absx.lt.1d-100) return
c
 10   if (absx.ge.0.5d0) go to 20
      n = n - 1
      absx = absx*2.0d0
      go to 10
c
 20   if (absx.lt.1.0d0) go to 30
      n = n + 1
      absx = absx*0.5d0
      go to 20
c
 30   y = dsign (absx, x)
      return
c
      end

      double precision function dcsevl (x, a, n)
c
c evaluate the n-term chebyshev series a at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
c
c             input arguments --
c x      dble prec value at which the series is to be evaluated.
c a      dble prec array of n terms of a chebyshev series.  in eval-
c        uating a, only half the first coef is summed.
c n      number of terms in array a.
c
      double precision a(n), x, twox, b0, b1, b2
c
      if (n.lt.1) call seteru ('dcsevl  number of terms le 0', 28, 2,2)
      if (n.gt.1000) call seteru ('dcsevl  number of terms gt 1000',
     1  31, 3, 2)
      if (x.lt.(-1.1d0) .or. x.gt.1.1d0) call seteru (
     1  'dcsevl  x outside (-1,+1)', 25, 1, 1)
c
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
c
      dcsevl = 0.5d0 * (b0-b2)
c
      return
      end

      function inits (os, nos, eta)
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c initialize the orthogonal series so that inits is the number of terms
c needed to insure the error is no larger than eta.  ordinarily, eta
c will be chosen to be one-tenth machine precision.
c
c             input arguments --
c os     array of nos coefficients in an orthogonal series.
c nos    number of coefficients in os.
c eta    requested accuracy of series.
c
      dimension os(nos)
c
      if (nos.lt.1) call seteru (
     1  'inits   number of coefficients lt 1', 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call seteru ('inits   eta may be too small', 28,
     1  1, 2)
      inits = i
c
      return
      end
      subroutine r9upak (x, y, n)
c august 1980 portable edition.  w. fullerton, los alamos scientific lab
c
c unpack floating point number x so that x = y * 2.0**n, where
c 0.5 .le. abs(y) .lt. 1.0 .
c
      absx = abs(x)
      n = 0
      y = 0.0
      if (absx.lt.1d-100) return
c
 10   if (absx.ge.0.5) go to 20
      n = n - 1
      absx = absx*2.0
      go to 10
c
 20   if (absx.lt.1.0) go to 30
      n = n + 1
      absx = absx*0.5
      go to 20
c
 30   y = sign (absx, x)
      return
c
      end
      subroutine seterr (messg, nmessg, nerr, iopt)
c
c  this version modified by w. fullerton to dump if iopt = 1 and
c  not recovering.
c  seterr sets lerror = nerr, optionally prints the message and dumps
c  according to the following rules...
c
c    if iopt = 1 and recovering      - just remember the error.
c    if iopt = 1 and not recovering  - print, dump and stop.
c    if iopt = 2                     - print, dump and stop.
c
c  input
c
c    messg  - the error message.
c    nmessg - the length of the message, in characters.
c    nerr   - the error number. must have nerr non-zero.
c    iopt   - the option. must have iopt=1 or 2.
c
c  error states -
c
c    1 - message length not positive.
c    2 - cannot have nerr=0.
c    3 - an unrecovered error followed by another error.
c    4 - bad value for iopt.
c
c  only the first 72 characters of the message are printed.
c
c  the error handler calls a subroutine named fdump to produce a
c  symbolic dump. to complete the package, a dummy version of fdump
c  is supplied, but it should be replaced by a locally written version
c  which at least gives a trace-back.
c
      integer messg(1)
      external i1mach, i8save
c
c  the unit for error messages.
c
      iwunit=i1mach(4)
c
      if (nmessg.ge.1) go to 10
c
c  a message of non-positive length is fatal.
c
        write(iwunit,9000)
 9000   format(52h1error    1 in seterr - message length not positive.)
        go to 60
c
c  nw is the number of words the message occupies.
c
 10   nw=(min0(nmessg,72)-1)/i1mach(6)+1
c
      if (nerr.ne.0) go to 20
c
c  cannot turn the error state off using seterr.
c
        write(iwunit,9001)
 9001   format(42h1error    2 in seterr - cannot have nerr=0//
     1         34h the current error message follows///)
        call e9rint(messg,nw,nerr,.true.)
        itemp=i8save(1,1,.true.)
        go to 50
c
c  set lerror and test for a previous unrecovered error.
c
 20   if (i8save(1,nerr,.true.).eq.0) go to 30
c
        write(iwunit,9002)
 9002   format(23h1error    3 in seterr -,
     1         48h an unrecovered error followed by another error.//
     2         48h the previous and current error messages follow.///)
        call eprint
        call e9rint(messg,nw,nerr,.true.)
        go to 50
c
c  save this message in case it is not recovered from properly.
c
 30   call e9rint(messg,nw,nerr,.true.)
c
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
c
c  must have iopt = 1 or 2.
c
        write(iwunit,9003)
 9003   format(42h1error    4 in seterr - bad value for iopt//
     1         34h the current error message follows///)
        go to 50
c
c  test for recovery.
c
 40   if (iopt.eq.2) go to 50
c
      if (i8save(2,0,.false.).eq.1) return
c
c     call eprint
c     stop
c
 50   call eprint
 60   call fdump
      stop
c
      end
      subroutine fdump
      return
      end

      subroutine e9rint(messg,nw,nerr,save)
c
c  this routine stores the current error message or prints the old one,
c  if any, depending on whether or not save = .true. .
c
      integer messg(nw)
      logical save
      external i1mach, i8save
c
c  messgp stores at least the first 72 characters of the previous
c  message. its length is machine dependent and must be at least
c
c       1 + 71/(the number of characters stored per integer word).
c
      integer messgp(36),fmt(14),ccplus
c
c  start with no previous message.
c
      data messgp(1)/1h1/, nwp/0/, nerrp/0/
c
c  set up the format for printing the error message.
c  the format is simply (a1,14x,72axx) where xx=i1mach(6) is the
c  number of characters stored per integer word.
c
      data ccplus  / 1h+ /
c
      data fmt( 1) / 1h( /
      data fmt( 2) / 1ha /
      data fmt( 3) / 1h1 /
      data fmt( 4) / 1h, /
      data fmt( 5) / 1h1 /
      data fmt( 6) / 1h4 /
      data fmt( 7) / 1hx /
      data fmt( 8) / 1h, /
      data fmt( 9) / 1h7 /
      data fmt(10) / 1h2 /
      data fmt(11) / 1ha /
      data fmt(12) / 1hx /
      data fmt(13) / 1hx /
      data fmt(14) / 1h) /
c
      if (.not.save) go to 20
c
c  save the message.
c
        nwp=nw
        nerrp=nerr
        do 10 i=1,nw
 10     messgp(i)=messg(i)
c
        go to 30
c
 20   if (i8save(1,0,.false.).eq.0) go to 30
c
c  print the message.
c
        iwunit=i1mach(4)
        write(iwunit,9000) nerrp
 9000   format(' error ' ,i4,' in ')
c
        call s88fmt(2,i1mach(6),fmt(12))
        write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
c
 30   return
c
      end
      subroutine eprint
c
c  this subroutine prints the last error message, if any.
c
      integer messg(1)
c
      call e9rint(messg,1,1,.false.)
      return
c
      end
      subroutine s88fmt( n, w, ifmt )
c
c  s88fmt  replaces ifmt(1), ... , ifmt(n) with
c  the characters corresponding to the n least significant
c  digits of w.
c
      integer n,w,ifmt(n)
c
      integer nt,wt,digits(10)
c
      data digits( 1) / 1h0 /
      data digits( 2) / 1h1 /
      data digits( 3) / 1h2 /
      data digits( 4) / 1h3 /
      data digits( 5) / 1h4 /
      data digits( 6) / 1h5 /
      data digits( 7) / 1h6 /
      data digits( 8) / 1h7 /
      data digits( 9) / 1h8 /
      data digits(10) / 1h9 /
c
      nt = n
      wt = w
c
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
c
      end
      double precision function dci (x)
c december 1980 edition, w. fullerton, bell labs.
      double precision x, cics(19), xsml, y, f, g, sinx, dcsevl,
     1  d1mach, cos, log, sin, dsqrt
      external d1mach, dcsevl, initds
c
c series for ci   on the interval  0.00000e+00 to  1.60000e+01
c                                        with weighted error   9.26e-33
c                                         log weighted error  32.03
c                               significant figures required  32.06
c                                    decimal places required  32.67
c
      data ci  cs(  1) / -0.3400428185 6055363156 2810766331 29873d0/
      data ci  cs(  2) / -1.0330216640 1177456807 1592710401 63751d0/
      data ci  cs(  3) /  0.1938822265 9917082876 7158746060 81709d0/
      data ci  cs(  4) / -0.0191826043 6019865893 9463462701 75301d0/
      data ci  cs(  5) /  0.0011078925 2584784967 1840980992 66118d0/
      data ci  cs(  6) / -0.0000415723 4558247208 8038402318 14601d0/
      data ci  cs(  7) /  0.0000010927 8524300228 7152955789 66285d0/
      data ci  cs(  8) / -0.0000000212 3285954183 4652196012 80329d0/
      data ci  cs(  9) /  0.0000000003 1733482164 3485448651 29873d0/
      data ci  cs( 10) / -0.0000000000 0376141547 9876836993 81798d0/
      data ci  cs( 11) /  0.0000000000 0003622653 4884839643 36956d0/
      data ci  cs( 12) / -0.0000000000 0000028911 5284936518 52433d0/
      data ci  cs( 13) /  0.0000000000 0000000194 3278606764 94420d0/
      data ci  cs( 14) / -0.0000000000 0000000001 1151831826 50184d0/
      data ci  cs( 15) /  0.0000000000 0000000000 0055278588 87706d0/
      data ci  cs( 16) / -0.0000000000 0000000000 0000239070 13943d0/
      data ci  cs( 17) /  0.0000000000 0000000000 0000000910 01612d0/
      data ci  cs( 18) / -0.0000000000 0000000000 0000000003 07233d0/
      data ci  cs( 19) /  0.0000000000 0000000000 0000000000 00926d0/
c
      data nci, xsml /0, 0.0d0/
c
      if (nci.ne.0) go to 10
      nci = initds (cics, 19, 1e-12) !0.1*sngl(d1mach(3)))
      xsml = dsqrt (d1mach(3))
c
 10   if (x.le.0.0d0) call seteru ('dci     x is le 0', 17, 1, 2)
c
      if (x.gt.4.0d0) go to 20
      y = -1.0d0
      if (x.gt.xsml) y = (x*x-8.d0)*0.125d0
c
      dci = log(x) - 0.5d0 + dcsevl (y, cics, nci)
      return
c
 20   call d9sifg (x, f, g)
      sinx = sin (x)
      call erroff
      dci = f*sinx - g*cos(x)
c
      return
      end

      INTEGER FUNCTION I1MACH(I)
C
C  I/O UNIT NUMBERS.
C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS.
C
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER CHARACTER STORAGE UNIT.
C                 FOR FORTRAN 77, THIS IS ALWAYS 1.  FOR FORTRAN 66,
C                 CHARACTER STORAGE UNIT = INTEGER STORAGE UNIT.
C
C  INTEGERS.
C
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C
C    I1MACH( 7) = A, THE BASE.
C
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C
C    I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.  FOR FORTRAN 77, YOU MAY WISH
C  TO ADJUST THE DATA STATEMENT SO IMACH(6) IS SET TO 1, AND
C  THEN TO COMMENT OUT THE EXECUTABLE TEST ON I .EQ. 6 BELOW.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE, EXCEPT PERHAPS
C  FOR IMACH(1) - IMACH(4).
C
      INTEGER IMACH(16),OUTPUT,SANITY
C
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
      DATA IMACH( 1) /    5 /
      DATA IMACH( 2) /    6 /
      DATA IMACH( 3) /    7 /
      DATA IMACH( 4) /    6 /
      DATA IMACH( 5) /   32 /
      DATA IMACH( 6) /    4 /
      DATA IMACH( 7) /    2 /
      DATA IMACH( 8) /   31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /    2 /
      DATA IMACH(11) /   24 /
      DATA IMACH(12) / -125 /
      DATA IMACH(13) /  128 /
      DATA IMACH(14) /   53 /
      DATA IMACH(15) / -1021 /
      DATA IMACH(16) /  1024 /, SANITY/987/
C  ***  ISSUE STOP 777 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SANITY .NE. 987) STOP 777
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 999
      I1MACH=IMACH(I)
C/6S
C/7S
      IF(I.EQ.6) I1MACH=1
C/
      RETURN
  999 WRITE(OUTPUT,1999) I
 1999 FORMAT(' I1MACH - I OUT OF BOUNDS',I10)
      STOP
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES AND MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.
C
      DATA SMALL(1),SMALL(2) /    1048576,          0 /
      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
      DATA DIVER(1),DIVER(2) / 1018167296,          0 /
      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /
      IF (I .LT. 1  .OR.  I .GT. 5) GOTO 999
      D1MACH = DMACH(I)
      RETURN
 999  WRITE(I1MACH(2),1999) I
 1999 FORMAT(' D1MACH - I OUT OF BOUNDS',I10)
      STOP
      END

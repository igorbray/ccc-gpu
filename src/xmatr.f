      subroutine QRGRIDn(J,MJ)                                   

c Full spherical integration
      
!      IMPLICIT REAL*8(A-H,O-Z)

      include 'par.f'
      common/meshrr/ meshr,rmesh(maxr,3)
      dimension ii(maxr)
      PARAMETER (Ncf = 20)
      COMMON /cat2/   Nc, nmax(0:Ncf), jmax
      COMMON /cat5/   maxNR(0:Ncf, Ncf)

      COMMON /JZ/  Jt
      COMMON /RR/  rt,rp
      COMMON /QGR/ fqr,fqr2
      DIMENSION    fqr(maxr),fqr2(maxr,1)!maxr)
      COMMON /const/  pi, Ry, rad, SIP
      COMMON /KNMT/   E0,E1,theta1
      common /KK/     k0,k1,p,teta1
      real          k0,k1

      logical there
      character nh(0:9)
      data nh / '0','1','2','3','4','5','6','7','8','9'/


      PARAMETER(NN=1000)
      real*8 AGS(NN),WGS(NN),pid
      stop 'redefine fqr2(maxr,maxr) in QGRIDn of xmatr.f'
      pid = dacos(-1d0)
      rad = pid/180
      teta1 = theta1*rad
      
      k0= sqrt(E0/Ry)
      k1= sqrt(E1/Ry)
      p = sqrt(k0*k1)
      
      NG  = 50
      x1=cos(teta1)
      y1=sin(teta1)
      q = sqrt(k0**2 + k1**2-2*k0*k1*x1)
      print'(A,3F9.4)', 'Ks', k0,k1,p
      Maxn = 1
      do l = 0, jmax
         do n = l+1, nmax(l)
            if(maxNR(l,n).gt.Maxn) maxn = maxNR(l,n)
         end do
      end do
      print'(3(A,I3))', 'QGRID J=', J,' MJ=',MJ, ' Maxn=', maxn
c$$$      if(maxn.gt.Nax) STOP "INCREASE ARRAY SIZE"
      
      inquire(file='matrix.J='//nh(J)//nh(MJ),exist=there)
      if (there) then
         open(301,file='matrix.J='//nh(J)//nh(MJ),
     >      STATUS='OLD',FORM='UNFORMATTED')
         print'(A)', 'MATRIX EXIST'
         read(301) p0,p1,pp
         eps=(p0-k0)**2 + (p1-k1)**2 + (p-pp)**2
         read(301)   j1,  maxn1
         if(eps.gt.1.e-4 .or. j.ne.j1 .or. maxn.ne.maxn1)
     >      STOP "WRONG MATRIX"
         do i=1,Maxn
            do ip=1,Maxn
               read(301) fqr2(i,ip)
            end do
            rt =  rmesh(i,1)
            fqr(i) =  fqr2(i,i)
c$$$            write(101,'(F9.4,5E13.4)') rt, fqr(i)
            write(501,'(5000E13.4)') (fqr2(i,ip),ip=1,Maxn)
         end do
         goto 111
      else 
         open(301,file='matrix.J='//nh(J)//nh(MJ),
     >      STATUS='NEW',FORM='UNFORMATTED')
         open(201,file='vector',FORM='FORMATTED')
         print'(A)', 'MATRIX NOT FOUND'
         call clock(s1)
         write(301) k0,k1,p
         write(301)   J, maxn
      end if

      do it=1,Maxn
         rt =  rmesh(it,1)
         do ip=1,Maxn
            rp =  rmesh(ip,1)
            sum = 0d0


            call WAGAUS(0d0,pid,AGS,WGS,NG)
C$OMP PARALLEL DO
C$OMP& default(SHARED)
C$OMP& reduction(+: sum)
C$OMP& PRIVATE(i,theta,x,y,q0,w)
            do i = 1,NG
               theta = ags(i)
               x=cos(theta)
               y=sin(theta)
               q0 = (p**2+k0**2-2.*p*k0*x)**0.5
               w = wgs(i)
               
               if(J.eq.0) then
                  sum = sum + y*FUQ0(x) * w
     >              *(FUQ1(theta)+(q0**2-q**2)*FUQ2(theta))
               else if(J.eq.2 .and. MJ.eq.0) then
                  sum = sum + y*FUQ0(x) * w
     >               * (FUQ1(theta) + (q0**2-q**2 +  6*(k0-p*x)
     >               *     (p*x-k1*x1))*FUQ2(theta))
               else if(J.eq.2 .and. MJ.eq.1) then
                  sum = sum + y*FUQ0(x) * w
     >               * (p*y*(2*p*x-k1*x1-k0)*FUQ3(theta)
     >               + k1*y1*(k0-p*x)*FUQ2(theta))
               else if(J.eq.2 .and. MJ.eq.2) then
                  sum = sum + y*FUQ0(x) * w
     >               * (k1*y*y1*FUQ3(theta) - p*y**2*FUQ4(theta))
               end if
            end do
C$OMP END PARALLEL DO 

            if(J.eq.0) then
               sum = sum*9*sqrt(3.0)*p/64/pi**2
            else if(J.eq.2 .and. MJ.eq.0) then
               sum = sum*9*sqrt(6.0)*p/128/pi**2
            else if(J.eq.2 .and. MJ.eq.1) then
               sum = sum*27*p   /64/pi**2
            else if(J.eq.2 .and. MJ.eq.2) then
               sum = sum*27*p**2/64/pi**2
            end if
            if(it.eq.1 .and. ip.eq.1) sum1=sum
            fqr2(it,ip) = sum!/sum1
            write(301) fqr2(it,ip)
         end do
         fqr(it)= fqr2(it,it)
         if(it.eq.1) print'(2(A,I3),E13.4)','IKF J=', J,' MJ=',MJ,sum1 
c$$$         write(101,'(F9.4,5E13.4)') rt, fqr(it)
         write(201,'(100E13.4)') rt
      end do
      call clock(s2)
      print*,'Time for QRgrid:',s2-s1


      
 111  continue
      close(301)

      return
      end

******************************************************************

      function FUQ0(x)


      common /KK/     k0,k1,p,theta1
      common /RR/     r,rp
      
      real k0,k1
      real*8 z
      
      q0 = (p**2+k0**2-2.*p*k0*x)**0.5
      z = q0*r
      FUQ0 = bes(1,z)/r/q0**3

      return
      end
      
******************************************************************

      function FUQ1(theta)

      PARAMETER(NN=1000)
      
      real*8 AGS(NN),WGS(NN),twopi

      common /KK/     k0,k1,p,theta1
      common /RR/     r,rp
      real k0,k1
      real*8 z
      pi = dacos(-1d0)
      twopi = pi*2

      N  = 10
      sum = 0d0
      call WAGAUS(0d0,twopi,AGS,WGS,N)

      do i = 1,N
         phi = ags(i)
         w = wgs(i)
         x = cos(theta)*cos(theta1)+sin(theta)*sin(theta1)*cos(phi)
         q1= (p**2+k1**2-2.*p*k1*x)**0.5
         z = q1*rp
         sum = sum + bes(1,z)/rp/q1*w
!         print'(2F9.4,3E13.4)', theta,sin(theta),x,w,sum
      end do

      FUQ1 = sum

      return
      end
      
******************************************************************

      function FUQ2(theta)

      PARAMETER(NN=1000)
      
      real*8 AGS(NN),WGS(NN)

      common /KK/     k0,k1,p,theta1
      common /RR/     r,rp
      real k0,k1
      real*8 z,twopi
      pi = dacos(-1d0)
      twopi = pi*2

      N  = 10
      sum = 0d0
      call WAGAUS(0d0,twopi,AGS,WGS,N)

      do i = 1,N
         phi = ags(i)
         w = wgs(i)
         x = cos(theta)*cos(theta1)+sin(theta)*sin(theta1)*cos(phi)
         q1= (p**2+k1**2-2.*p*k1*x)**0.5
         z = q1*rp
          
         sum = sum + bes(1,z)/rp/q1**3*w
!         print'(2F9.4,3E13.4)', theta,sin(theta),x,w,sum
      end do

      FUQ2 = sum

      return
      end
      
******************************************************************

      function FUQ3(theta)

      PARAMETER(NN=1000)
      
      real*8 AGS(NN),WGS(NN)

      common /KK/     k0,k1,p,theta1
      common /RR/     r,rp
      real k0,k1
      real*8 z,twopi
      pi = dacos(-1d0)
      twopi = pi*2

      q = sqrt(k0**2 + k1**2-2*k0*k1*cos(theta1))
      q0 = (p**2+k0**2-2.*p*k0*cos(theta))**0.5
      
      N  = 10
      sum = 0d0
      call WAGAUS(0d0,twopi,AGS,WGS,N)

      do i = 1,N
         phi = ags(i)
         w = wgs(i)
         x = cos(theta)*cos(theta1)+sin(theta)*sin(theta1)*cos(phi)
         q1= (p**2+k1**2-2.*p*k1*x)**0.5
         z = q1*rp
          
         sum = sum + bes(1,z)/rp*cos(phi)/q1**3*w
!         print'(2F9.4,3E13.4)', theta,sin(theta),x,w,sum
      end do

      FUQ3 = sum

      return
      end
      
******************************************************************

      function FUQ4(theta)

      PARAMETER(NN=1000)
      
      real*8 AGS(NN),WGS(NN)

      common /KK/     k0,k1,p,theta1
      common /RR/     r,rp
      real k0,k1
      real*8 z,twopi
      pi = dacos(-1d0)
      twopi = pi*2

      q = sqrt(k0**2 + k1**2-2*k0*k1*cos(theta1))
      q0 = (p**2+k0**2-2.*p*k0*cos(theta))**0.5
      
      N  = 10
      sum = 0d0
      call WAGAUS(0d0,twopi,AGS,WGS,N)

      do i = 1,N
         phi = ags(i)
         w = wgs(i)
         x = cos(theta)*cos(theta1)+sin(theta)*sin(theta1)*cos(phi)
         q1= (p**2+k1**2-2.*p*k1*x)**0.5
         z = q1*rp
          
         sum = sum + bes(1,z)/rp*cos(phi*2)/q1**3*w
!         print'(2F9.4,3E13.4)', theta,sin(theta),x,w,sum
      end do

      FUQ4 = sum

      return
      end
      




       
      SUBROUTINE WAGAUS(A,B,AGS,WGS,NGAUSS)
C   ******************************************************************
C   *                                                                *
C   *    Obtain abscissas and weight factors for Gauss-Legendre      *
C   *    n-point integration in the range [A,B]. To see abscissas    *
C   *    and weigths in the reduced range [-1,1] set IDBG = 1.       *
C   *                                                                *
C   *    A     : Lower integration limit                             *
C   *    B     : Upper integration limit                             *
C   *    NGAUSS: Number of points in the range [A,B]                 *
C   *    AGS   : Abscissas in the range [A,B]                        *
C   *    WGS   : Corresponding weight factors                        *
C   *    EPS   : Floating precission in the root-finding             *
C   *                                                                *
C   *    HLS: 30-JAN-87                                              *
C   ******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ONE=1,EPS=1.D-15,IDBG=0)
      DIMENSION AGS(NGAUSS),WGS(NGAUSS)
  101 FORMAT('0',10X,'Mesh and weights found in WAGAUSS',//,3X,'I',9X,
     1       'X',17X,'WEIGHT',/)
  102 FORMAT(' ',I3,2F20.15)
  103 FORMAT('0WAGAUS: EPS too small for machine precission',/,8X,
     1       'EPS = ',E15.8)
C
      PI=ACOS(-ONE)
      NHALF=(NGAUSS+1)/2
C
C     Guess at the I'th root of the Legendre polynomial and refine
C     by Newton to obtain abscissas and weights in the normalized
C     range [-1,1]
C
      DO 11 I=1,NHALF
      Z=COS(PI*(I-.25D0)/(NGAUSS+0.5D0))
      NEWTON=0
   10    NEWTON=NEWTON+1
         IF(NEWTON.GT.20) THEN
            WRITE(6,103) EPS
            STOP
         ENDIF
         P1=1.D0
         P2=0.D0
         DO 12 J=1,NGAUSS
         P3=P2
         P2=P1
         P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
   12    CONTINUE
C
C        Newton root-finding using the derivative PP
C
         PP=NGAUSS*(Z*P1-P2)/(Z*Z-1.D0)
         Z1=Z
         Z=Z1-P1/PP
      IF(ABS(Z-Z1).GT.EPS) GO TO 10
      AGS(I)=-Z
      AGS(NGAUSS+1-I)=Z
      WGS(I)=2.D0/((1.D0-Z*Z)*PP*PP)
      WGS(NGAUSS+1-I)=WGS(I)
   11 CONTINUE
      IF(IDBG.EQ.1) THEN
         WRITE(6,101)
         DO 16 I=1,NGAUSS
   16    WRITE(6,102) I,AGS(I),WGS(I)
      ENDIF
C
C     Transform abscissas and weights to the range [A,B]
C
      C=(B-A)/2.0D0
      D=(B+A)/2.0D0
      DO 20 I=1,NGAUSS
      AGS(I)=C*AGS(I)+D
   20 WGS(I)=C*WGS(I)
      RETURN
      END

      block data
      implicit real*8 (a-h,o-z)
      include 'par.f'
      COMMON/GMFN/FS(34)
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),res(63)
      DATA X2( 1,1),W2( 1,1)/ 0.7745966692415D+00, 0.5555555555555D+00/
      DATA X2( 2,1),W2( 2,1)/-0.5421010862427D-18, 0.8888888888889D+00/
      DATA X2( 3,1),W2( 3,1)/-0.7745966692415D+00, 0.5555555555555D+00/
      DATA X2( 1,2),W2( 1,2)/ 0.7745966692415D+00, 0.2684880898683D+00/
      DATA X2( 2,2),W2( 2,2)/-0.5421010862427D-18, 0.4509165386585D+00/
      DATA X2( 3,2),W2( 3,2)/-0.7745966692415D+00, 0.2684880898683D+00/
      DATA X2( 4,2),W2( 4,2)/ 0.9604912687080D+00, 0.1046562260265D+00/
      DATA X2( 5,2),W2( 5,2)/ 0.4342437493468D+00, 0.4013974147759D+00/
      DATA X2( 6,2),W2( 6,2)/-0.4342437493468D+00, 0.4013974147759D+00/
      DATA X2( 7,2),W2( 7,2)/-0.9604912687080D+00, 0.1046562260265D+00/
      DATA X2( 1,3),W2( 1,3)/ 0.7745966692415D+00, 0.1344152552438D+00/
      DATA X2( 2,3),W2( 2,3)/-0.5421010862427D-18, 0.2255104997982D+00/
      DATA X2( 3,3),W2( 3,3)/-0.7745966692415D+00, 0.1344152552438D+00/
      DATA X2( 4,3),W2( 4,3)/ 0.9604912687080D+00, 0.5160328299708D-01/
      DATA X2( 5,3),W2( 5,3)/ 0.4342437493468D+00, 0.2006285293770D+00/
      DATA X2( 6,3),W2( 6,3)/-0.4342437493468D+00, 0.2006285293770D+00/
      DATA X2( 7,3),W2( 7,3)/-0.9604912687080D+00, 0.5160328299708D-01/
      DATA X2( 8,3),W2( 8,3)/ 0.9938319632127D+00, 0.1700171962994D-01/
      DATA X2( 9,3),W2( 9,3)/ 0.8884592328722D+00, 0.9292719531512D-01/
      DATA X2(10,3),W2(10,3)/ 0.6211029467372D+00, 0.1715119091364D+00/
      DATA X2(11,3),W2(11,3)/ 0.2233866864290D+00, 0.2191568584016D+00/
      DATA X2(12,3),W2(12,3)/-0.2233866864290D+00, 0.2191568584016D+00/
      DATA X2(13,3),W2(13,3)/-0.6211029467372D+00, 0.1715119091364D+00/
      DATA X2(14,3),W2(14,3)/-0.8884592328722D+00, 0.9292719531512D-01/
      DATA X2(15,3),W2(15,3)/-0.9938319632127D+00, 0.1700171962994D-01/
      DATA X2( 1,4),W2( 1,4)/ 0.7745966692415D+00, 0.6720775429600D-01/
      DATA X2( 2,4),W2( 2,4)/-0.5421010862427D-18, 0.1127552567208D+00/
      DATA X2( 3,4),W2( 3,4)/-0.7745966692415D+00, 0.6720775429599D-01/
      DATA X2( 4,4),W2( 4,4)/ 0.9604912687080D+00, 0.2580759809619D-01/
      DATA X2( 5,4),W2( 5,4)/ 0.4342437493468D+00, 0.1003142786118D+00/
      DATA X2( 6,4),W2( 6,4)/-0.4342437493468D+00, 0.1003142786118D+00/
      DATA X2( 7,4),W2( 7,4)/-0.9604912687080D+00, 0.2580759809620D-01/
      DATA X2( 8,4),W2( 8,4)/ 0.9938319632127D+00, 0.8434565739273D-02/
      DATA X2( 9,4),W2( 9,4)/ 0.8884592328722D+00, 0.4646289326174D-01/
      DATA X2(10,4),W2(10,4)/ 0.6211029467372D+00, 0.8575592004998D-01/
      DATA X2(11,4),W2(11,4)/ 0.2233866864290D+00, 0.1095784210559D+00/
      DATA X2(12,4),W2(12,4)/-0.2233866864290D+00, 0.1095784210559D+00/
      DATA X2(13,4),W2(13,4)/-0.6211029467372D+00, 0.8575592004998D-01/
      DATA X2(14,4),W2(14,4)/-0.8884592328722D+00, 0.4646289326175D-01/
      DATA X2(15,4),W2(15,4)/-0.9938319632127D+00, 0.8434565739234D-02/
      DATA X2(16,4),W2(16,4)/ 0.9990981249676D+00, 0.2544780791599D-02/
      DATA X2(17,4),W2(17,4)/ 0.9815311495537D+00, 0.1644604985439D-01/
      DATA X2(18,4),W2(18,4)/ 0.9296548574297D+00, 0.3595710330713D-01/
      DATA X2(19,4),W2(19,4)/ 0.8367259381688D+00, 0.5697950949412D-01/
      DATA X2(20,4),W2(20,4)/ 0.7024962064915D+00, 0.7687962049900D-01/
      DATA X2(21,4),W2(21,4)/ 0.5313197436444D+00, 0.9362710998126D-01/
      DATA X2(22,4),W2(22,4)/ 0.3311353932580D+00, 0.1056698935802D+00/
      DATA X2(23,4),W2(23,4)/ 0.1124889431332D+00, 0.1119568730210D+00/
      DATA X2(24,4),W2(24,4)/-0.1124889431332D+00, 0.1119568730210D+00/
      DATA X2(25,4),W2(25,4)/-0.3311353932580D+00, 0.1056698935802D+00/
      DATA X2(26,4),W2(26,4)/-0.5313197436444D+00, 0.9362710998126D-01/
      DATA X2(27,4),W2(27,4)/-0.7024962064915D+00, 0.7687962049900D-01/
      DATA X2(28,4),W2(28,4)/-0.8367259381688D+00, 0.5697950949412D-01/
      DATA X2(29,4),W2(29,4)/-0.9296548574297D+00, 0.3595710330712D-01/
      DATA X2(30,4),W2(30,4)/-0.9815311495538D+00, 0.1644604985441D-01/
      DATA X2(31,4),W2(31,4)/-0.9990981249676D+00, 0.2544780791612D-02/
      DATA X2( 1,5),W2( 1,5)/ 0.7745966692415D+00, 0.3360387714799D-01/
      DATA X2( 2,5),W2( 2,5)/-0.5421010862427D-18, 0.5637762836039D-01/
      DATA X2( 3,5),W2( 3,5)/-0.7745966692415D+00, 0.3360387714713D-01/
      DATA X2( 4,5),W2( 4,5)/ 0.9604912687080D+00, 0.1290380022626D-01/
      DATA X2( 5,5),W2( 5,5)/ 0.4342437493468D+00, 0.5015713930586D-01/
      DATA X2( 6,5),W2( 6,5)/-0.4342437493468D+00, 0.5015713930582D-01/
      DATA X2( 7,5),W2( 7,5)/-0.9604912687080D+00, 0.1290379999680D-01/
      DATA X2( 8,5),W2( 8,5)/ 0.9938319632127D+00, 0.4217731327956D-02/
      DATA X2( 9,5),W2( 9,5)/ 0.8884592328722D+00, 0.2323144663721D-01/
      DATA X2(10,5),W2(10,5)/ 0.6211029467372D+00, 0.4287796002535D-01/
      DATA X2(11,5),W2(11,5)/ 0.2233866864290D+00, 0.5478921052795D-01/
      DATA X2(12,5),W2(12,5)/-0.2233866864290D+00, 0.5478921052795D-01/
      DATA X2(13,5),W2(13,5)/-0.6211029467372D+00, 0.4287796002531D-01/
      DATA X2(14,5),W2(14,5)/-0.8884592328722D+00, 0.2323144664337D-01/
      DATA X2(15,5),W2(15,5)/-0.9938319632127D+00, 0.4217614778805D-02/
      DATA X2(16,5),W2(16,5)/ 0.9990981249676D+00, 0.1264660920956D-02/
      DATA X2(17,5),W2(17,5)/ 0.9815311495537D+00, 0.8223012874950D-02/
      DATA X2(18,5),W2(18,5)/ 0.9296548574297D+00, 0.1797855141828D-01/
      DATA X2(19,5),W2(19,5)/ 0.8367259381688D+00, 0.2848975473281D-01/
      DATA X2(20,5),W2(20,5)/ 0.7024962064915D+00, 0.3843981024901D-01/
      DATA X2(21,5),W2(21,5)/ 0.5313197436444D+00, 0.4681355499089D-01/
      DATA X2(22,5),W2(22,5)/ 0.3311353932580D+00, 0.5283494679008D-01/
      DATA X2(23,5),W2(23,5)/ 0.1124889431332D+00, 0.5597843651048D-01/
      DATA X2(24,5),W2(24,5)/-0.1124889431332D+00, 0.5597843651051D-01/
      DATA X2(25,5),W2(25,5)/-0.3311353932580D+00, 0.5283494679007D-01/
      DATA X2(26,5),W2(26,5)/-0.5313197436444D+00, 0.4681355499065D-01/
      DATA X2(27,5),W2(27,5)/-0.7024962064915D+00, 0.3843981024984D-01/
      DATA X2(28,5),W2(28,5)/-0.8367259381688D+00, 0.2848975474162D-01/
      DATA X2(29,5),W2(29,5)/-0.9296548574297D+00, 0.1797855160991D-01/
      DATA X2(30,5),W2(30,5)/-0.9815311495538D+00, 0.8223009147808D-02/
      DATA X2(31,5),W2(31,5)/-0.9990981249676D+00, 0.1265258789067D-02/
      DATA X2(32,5),W2(32,5)/ 0.9998726561183D+00, 0.3634916215097D-03/
      DATA X2(33,5),W2(33,5)/ 0.9972064490897D+00, 0.2579198894668D-02/
      DATA X2(34,5),W2(34,5)/ 0.9886847496087D+00, 0.6115479896454D-02/
      DATA X2(35,5),W2(35,5)/ 0.9721828732593D+00, 0.1049824455438D-01/
      DATA X2(36,5),W2(36,5)/ 0.9463428580969D+00, 0.1540675033430D-01/
      DATA X2(37,5),W2(37,5)/ 0.9103711569593D+00, 0.2059423394988D-01/
      DATA X2(38,5),W2(38,5)/ 0.8639079381850D+00, 0.2586967933249D-01/
      DATA X2(39,5),W2(39,5)/ 0.8069405319532D+00, 0.3107355111651D-01/
      DATA X2(40,5),W2(40,5)/ 0.7397560443530D+00, 0.3606443278134D-01/
      DATA X2(41,5),W2(41,5)/ 0.6629096600252D+00, 0.4071551011704D-01/
      DATA X2(42,5),W2(42,5)/ 0.5771957100521D+00, 0.4491453165343D-01/
      DATA X2(43,5),W2(43,5)/ 0.4836180269457D+00, 0.4856433040658D-01/
      DATA X2(44,5),W2(44,5)/ 0.3833593241987D+00, 0.5158325395207D-01/
      DATA X2(45,5),W2(45,5)/ 0.2777498220218D+00, 0.5390549933528D-01/
      DATA X2(46,5),W2(46,5)/ 0.1682352515522D+00, 0.5548140435656D-01/
      DATA X2(47,5),W2(47,5)/ 0.5634431304660D-01, 0.5627769983125D-01/
      DATA X2(48,5),W2(48,5)/-0.5634431304658D-01, 0.5627769983124D-01/
      DATA X2(49,5),W2(49,5)/-0.1682352515522D+00, 0.5548140435655D-01/
      DATA X2(50,5),W2(50,5)/-0.2777498220218D+00, 0.5390549933529D-01/
      DATA X2(51,5),W2(51,5)/-0.3833593241987D+00, 0.5158325395209D-01/
      DATA X2(52,5),W2(52,5)/-0.4836180269458D+00, 0.4856433040668D-01/
      DATA X2(53,5),W2(53,5)/-0.5771957100519D+00, 0.4491453165350D-01/
      DATA X2(54,5),W2(54,5)/-0.6629096600249D+00, 0.4071551011672D-01/
      DATA X2(55,5),W2(55,5)/-0.7397560443534D+00, 0.3606443278112D-01/
      DATA X2(56,5),W2(56,5)/-0.8069405319506D+00, 0.3107355111347D-01/
      DATA X2(57,5),W2(57,5)/-0.8639079381891D+00, 0.2586967932669D-01/
      DATA X2(58,5),W2(58,5)/-0.9103711569487D+00, 0.2059423390077D-01/
      DATA X2(59,5),W2(59,5)/-0.9463428584183D+00, 0.1540675048263D-01/
      DATA X2(60,5),W2(60,5)/-0.9721828745704D+00, 0.1049824666711D-01/
      DATA X2(61,5),W2(61,5)/-0.9886847605064D+00, 0.6115510620291D-02/
      DATA X2(62,5),W2(62,5)/-0.9972062302276D+00, 0.2579021304445D-02/
      DATA X2(63,5),W2(63,5)/-0.9998729514522D+00, 0.3631587165155D-03/
      DATA FS/1.0D0,1.0D0,2.0D0,6.0D0,2.4D1,1.2D2,7.2D2,5.04D3,
     1 4032.0D1,36288.0D1,36288.0D2,399168.0D2,4790016.0D2,
     2 62270208.0D2,871782912.0D2,1307674368.0D3,
     3 20922789888.0D03,355687428096.0D03,640237370572.8D04,
     4 12164510040883.2D04,24329020081766.4D05,0.510909421717094D20,
     5 0.112400072777761D22,0.258520167388850D23,0.620448401733239D24,
     6 0.155112100433310D26,0.403291461126606D27,0.108888694504184D29,
     7 0.304888344611714D30,0.884176199373970D31,0.265252859812191D33,
     8 0.822283865417792D34,0.263130836933694D36,0.868331761881189D37/
      end
      
      subroutine initialize
      implicit real*8 (a-h,o-z)
      include 'par.f'
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),res(63)
      pi = acos(-1d0)
      factl(0)=0d0
      do 15 i=1,2*lcoul
         factl(i)=factl(i-1)+log(dfloat(i))
15    continue
      return
      end
      
C  The following function returns x*xp*FF+y*yp*GG, where FF is
C  the integral from r=a to r=infinity of F(l1,r1*r)*F(l2,r2*r)/r**(lambda+1),
C  and where GG is the integral from r=a to r=infinity of
C  G(l1,r1*r)*G(l2,r2*r)/r**(lambda+1). F(l,r)=r*j(l,r), where j(l,r)
C  is the spherical Bessel function. G(l,r)=r*n(l,r), where n(l,r) is the
C  Neumann spherical function.
      function ffgg(l1,r1,l2,r2,a,lambda,x,y,xp,yp)
      implicit real*8 (a-h,o-z)
      logical*4 cosl
      include 'par.f'
      dimension c(0:2*lcoul),jstart(2),jstop(2),jstep(2)
      sum=0d0
c$$$      aer=1d-5
      rer = 1d-4
      jm=nint(-(r1*a+r2*a+1)+sqrt((2*r1*a-1)**2+4*l1*(l1+1))/2
     >   +sqrt((2*r2*a-1)**2+4*l2*(l2+1))/2)
      ltot=l1+l2
      if (jm.lt.0) jm=0
      if (jm.gt.ltot) jm=ltot
      jstart(1)=jm
      jstop(1)=0
      jstep(1)=-1
      jstart(2)=jm+1
      jstop(2)=ltot
      jstep(2)=1
      do 50 i=1,2
         do 10 j=jstart(i),jstop(i),jstep(i)
            n=j+ltot
            cosl=mod(n,2).eq.0
            mn=min(l1,j)
            mx=max(0,j-l2)
            tmp=0d0
            do 20 m=mx,mn
               c(m)=cnst(l1,m,r1,l2,j-m,r2)
               tmp=tmp+c(m)
20          continue
            app=tmp/(a**(lambda+j)*(lambda+j))
c$$$            rer=aer/app
            if (rer.gt.1d0) go to 50
            if (rer.lt.1d-12) then
               ffgg=0d0
c               print*,'Unable to handle this case of FF+GG integral '
c               print*,l1,'=l1',l2,'=l2',sngl(r1),'=k1',sngl(r2),'=k2'
c               print*,lambda,'=lambda. Tail set to zero.'
               return
            end if
            tmp1=csint(a,r1+r2,lambda+1+j,cosl,rer)
            if (abs(r1/r2-1d0).gt.1d-4) then
               tmp2=csint(a,r2-r1,lambda+1+j,cosl,rer)
            else 
               rrr = 0d0
               tmp2=csint(a,rrr,lambda+1+j,cosl,rer)
            endif 
ci            print*,tmp1,a,r1+r2,lambda+1+j,cosl,rer
            tmp=tmp*tmp1*(y*yp-x*xp)
            do 30 m=mx,mn
               tmp=tmp+c(m)*(-1)**(m+l1)*tmp2*(x*xp+y*yp)
30          continue
            if (cosl) then
               sum=sum+tmp*(-1)**(n/2)
            else
               sum=sum+tmp*(-1)**((n+1)/2)
            end if
10       continue
50    continue
      ffgg=(-1)**(l1+l2)*sum/2d0
      end

C  The following function returns x*yp*FG+y*xp*GF, where FG is
C  the integral from r=a to r=infinity of F(l1,r1*r)*G(l2,r2*r)/r**(lambda+1),
C  and where GF is the integral from r=a to r=infinity of
C  G(l1,r1*r)*F(l2,r2*r)/r**(lambda+1). F(l,r)=r*j(l,r), where j(l,r)
C  is the spherical Bessel function. G(l,r)=r*n(l,r), where n(l,r) is the
C  Neumann spherical function.
      function fggf(l1,r1,l2,r2,a,lambda,x,y,xp,yp)
      implicit real*8 (a-h,o-z)
      logical*4 cosl
      include 'par.f'
      dimension c(0:2*lcoul),jstart(2),jstop(2),jstep(2)
      if (abs(y*yp).lt.1d-30) then
         fggf=0d0
         return
      end if
      sum=0d0
c$$$      aer=1d-5
      rer = 1d-4
      jm=nint(-(r1*a+r2*a+1)+sqrt((2*r1*a-1)**2+4*l1*(l1+1))/2
     >   +sqrt((2*r2*a-1)**2+4*l2*(l2+1))/2)
      ltot=l1+l2
      if (jm.lt.0) jm=0
      if (jm.gt.ltot) jm=ltot
      jstart(1)=jm
      jstop(1)=0
      jstep(1)=-1
      jstart(2)=jm+1
      jstop(2)=ltot
      jstep(2)=1
      yypmx=max(abs(y),abs(yp))
      do 50 i=1,2
         do 10 j=jstart(i),jstop(i),jstep(i)
            n=j+ltot
            cosl=mod(n,2).ne.0
            mn=min(l1,j)
            mx=max(0,j-l2)
            tmp=0d0
            do 20 m=mx,mn
               c(m)=cnst(l1,m,r1,l2,j-m,r2)
               tmp=tmp+c(m)
20          continue
            app=yypmx*tmp/(a**(lambda+j)*(lambda+j))
c$$$            rer=aer/app
            if (rer.gt.1d0) go to 50
            if (rer.lt.1d-12) then
               fggf=0d0
c               print*,'Unable to handle this case of FG+GF integral '
c               print*,l1,'=l1',l2,'=l2',sngl(r1),'=k1',sngl(r2),'=k2'
c               print*,lambda,'=lambda. Tail set to zero.'
               return
            end if
            tmp1=csint(a,r1+r2,lambda+1+j,cosl,rer)
            if (abs(r1/r2-1d0).gt.1d-4) then
               tmp2=csint(a,r2-r1,lambda+1+j,cosl,rer)
            else
               rrr = 0d0
               tmp2=csint(a,rrr,lambda+1+j,cosl,rer)
            endif                
            tmp=tmp*tmp1*(x*yp+xp*y)
            do 30 m=mx,mn
               tmp=tmp-c(m)*(-1)**(m+l1)*tmp2*(x*yp-xp*y)
30          continue
            if (.not.cosl) then
               sum=sum-tmp*(-1)**(n/2)
            else
               sum=sum+tmp*(-1)**((n+1)/2)
            end if
10       continue
50    continue
      fggf=(-1)**(l1+l2+1)*sum/2d0
      end

      function cnst(l1,i1,r1,l2,i2,r2)
      implicit real*8 (a-h,o-z)
      include 'par.f'
      common /cnsts2/ factl(0:2*lcoul),pi,x2(63,5),w2(63,5),res(63)
      cnst=exp(factl(l1+i1)+factl(l2+i2)-factl(i1)-factl(i2)
     >   -factl(l1-i1)-factl(l2-i2)-i1*log(2*r1)-i2*log(2*r2))
      end

      DOUBLE PRECISION FUNCTION GAMX(ITEN,I)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      include 'par.f'
      parameter(lll=lcoul*2+20)
      COMMON/GAMMA/GS(lll)
C
C ****  THIS SUBROUTINE RETURNS THE VALUE OF THE GAMMA
C ****  FUNCTION FOR X = I*0.5D0
C ****  IF ITEN = 0,  RETURN GAMX(X) = GAMMA(X)/(10**X)
C ****  IF ITEN > 0,  RETURN GAMX(X) = GAMMA(X)
C ****  IF I < 0 AND EVEN, NO VALUE IS RETURNED SINCE THE GAMMA
C ****  FUNCTION HAS POLES AT NEGATIVE INTEGERS
C
      scale = 1d3
      II = I
      IF(II.LE.0) THEN
         IF(MOD(II,2).EQ.0) THEN
            WRITE(6,20) II
            GAMX = 0.0D0
            RETURN
         END IF
         II = 1
      END IF
C
      IF(ITEN.EQ.0) THEN
         WGX = 1.0D0
      ELSE
!         WGX = SQRT(10.0D0)
         WGX = SQRT(scale)
      END IF
C
      GAMX = GS(II)*(WGX**II)
C
C ****  RECURSION FORMULA ARE USED TO EXTEND THE LOOK-UP TABLE
C ****  FOR THE NEGATIVE 1/2 INTEGERS
C
      IF(I.NE.II) THEN
         IF(ITEN.EQ.0) THEN
!            WGX = 10.0D0
            WGX = scale
         ELSE
            WGX = 1.0D0
         END IF
         II2 = II - 2
         DO 10 J = II2,I,-2
            GAMX = GAMX*WGX/(0.5D0*DBLE(J))
10       CONTINUE
      END IF
      RETURN
20    FORMAT(1X,' GAMMA FUNCTION UNDEFINED FOR X =',I3,'/2'
     >   ,'  GAMX = 0.0 RETURNED')
      END

      SUBROUTINE GAMSET
      IMPLICIT REAL*8(A-H,O-Z)
      include 'par.f'
      parameter(lll=lcoul*2+20)
      COMMON/GAMMA/GS(lll)
C
C ****  COMPUTES VALUES OF THE GAMMA FUNCTION TIMES A SCALE FACTOR.
C ****  GS(I)  CONTAINS GAMMA(I/2)/(10**(I/2))
C
      scale = 1d3
!      GS(2) = 1.0D0/10.0D0
      GS(2) = 1.0D0/scale
      DO 10 I = 4,lll,2
!         GS(I) = 0.5D0*DBLE(I-2)*GS(I-2)/10.0D0
         GS(I) = 0.5D0*DBLE(I-2)*GS(I-2)/scale
10    CONTINUE
!      GS(1) = 1.7724538509055D0/SQRT(10.0D0)
      GS(1) = 1.7724538509055D0/SQRT(scale)
      DO 20 I = 3,lll,2
!         GS(I) = 0.5D0*DBLE(I-2)*GS(I-2)/10.0D0
         GS(I) = 0.5D0*DBLE(I-2)*GS(I-2)/scale
20    CONTINUE
      RETURN
      END

      SUBROUTINE PSISET
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PSVAL/PSI(600)
C
C ****  COMPUTES THE VALUE OF THE PSI(DIGAMMA) FUNCTION.
C ****  DIGAMMA(I/2) IS STORED IN PSI(I)
C
      PSI(1) = -1.96351002602143D0
      PSI(2) = -0.577215664901533D0
      DO 10 I = 3,600
         PSI(I) = PSI(I-2) + 2.0D0/DBLE(I-2)
10    CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FDHY(A,B,C,NTERM,X)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
C
C ****  CALCULATES THE HYPERGEOMETRIC POWER SERIES TO AN ACCURACY
C ****  OF 10**(-7.0). A SIMPLE POWER SERIES ALGORITHM IS USED TO
C ****  DETERMINE WHERE TO TRUNCATE THE SUMMATION WHEN THE SERIES
C ****  IS INFINITE. WHEN NTERM IS A POSITIVE INTEGER THE SERIES
C ****  TERMINATES AND A DIFFERENT SECTION OF CODE IS USED.
C
      ABCF = 1.0D0
      XI = 1.0D0
      S1 = 0.0D0
      s2=1d0
      i=0
      IF(NTERM.LT.0) THEN
c         while(abs(s2-s1).gt.1d-7) do
 9999     if(abs(s2-s1).gt.1d-7)then
            i=i+1
            s2=s1
            S1 = S1 + ABCF*XI
            SI = DBLE(I-1)
            ABCF = ABCF*(A+SI)*(B+SI)/((C+SI)*(SI+1.0D0))
            XI = XI*X
            go to 9999
           endif
c         end while
      ELSE
         NN = NTERM + 1
         DO 20 I = 1,NN
            S1 = S1 + ABCF*XI
            IF(I.EQ.NN) GO TO 30
            SI = DBLE(I-1)
            ABCF = ABCF*(A+SI)*(B+SI)/((C+SI)*(SI+1.0D0))
            XI = XI*X
20       CONTINUE
      END IF
30    CONTINUE
      FDHY = S1
      RETURN
      END

      DOUBLE PRECISION FUNCTION FDHY2(A,B,C,X)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      COMMON/GMFN/FGS(34)
      COMMON/PSVAL/PSI(600)
C
C ****  THIS FUNCTION SUMS THE INFINITE PART OF THE ALTERNATE
C ****  FORM OF THE HYPERGEOMETRIC FUNCTION. THE ACCURACY TOLERANCE
C ****  IS 10**(-7.0D0)
C
      XLN = LOG(X)
      IA = INT(2.0D0*A+0.1D0)
      IB = INT(2.0D0*B+0.1D0)
      IC = INT(2.0D0*C+0.1D0)
      IM = INT(C+0.1D0)
      ABCF = 1.0D0/FGS(IM+1)
      XI = 1.0D0
      S1 = 0.0D0
      i=0
      s2=1d0
c      while(abs(s2-s1).gt.1d-7) do
 9997  if(abs(s2-s1).gt.1d-7)then
         i=i+1
         s2=s1
         I2 = I*2
         APSI = XLN - PSI(I2) - PSI(I2+IC) + PSI(I2-2+IA) + PSI(I2-2+IB)
         S1 = S1 + APSI*ABCF*XI
         SI = DBLE(I)
         SI1 = DBLE(I-1)
         ABCF = ABCF*(A+SI1)*(B+SI1)/((C+SI)*SI)
         XI = XI*X
c      end while
       go to 9997
       endif
      FDHY2 = S1
      IF(i+4.GT.(600/2)) THEN
         WRITE(6,20)A,B,C,X,i
      END IF
      RETURN
20    FORMAT(1H ,' PROBLEM IN FDHY2. WITH A=',G15.6,' B=',G15.6
     >   ,' C=',G15.6,' X=',G15.6,/
     >   ,' NN =',I5,' IS TOO LARGE')
      END


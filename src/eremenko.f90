!Charlie: gets real term for regB
subroutine cwf_asym(p,q,l,eta,result)
implicit none


    ! Variables
    double precision :: p
    double precision :: q
    integer :: l
    double precision :: eta,sigma,result,aa
    double complex :: bb,denom,ctmp


!interface
!        double precision function cwfnSigmaL(lt, etat)
!            integer, intent(in) :: lt
!            double precision, intent(in) :: etat
!        end function cwfnSigmaL
!    end interface

call cwfnrega(q,l,eta,aa)
call cwfnregb(p,q,l,eta,bb)
    ! The denominator
    denom = dcmplx((p - q), 0.d0)
    ctmp = dcmplx(0.d0, eta)
    denom = denom**ctmp

bb=bb/denom

result=aa*dreal(bb)*(p+q)

end subroutine cwf_asym
!Charlie
!Calculates magnitude of complex gamma function*exp(pi*|eta|/2.d0)
subroutine cwf_gamma_mag(l,eta,mag)
implicit none
integer, intent(in) :: l
    double precision, intent(in) :: eta
    double precision, intent(out) :: mag

double precision::denom,pi
integer :: i
! Working...
pi = 4.d0 * datan(1.d0)
!Zero value
!denom=sinh(pi*eta)
if (eta.gt.0.d0) then
denom=1.d0-exp(-eta*pi*2.d0)
else
denom=exp(eta*pi*2.d0)-1.d0
endif

mag=sqrt(2.d0*eta*pi/denom)
!print*,denom,mag

!recurence
if (l.gt.0) then
do i=1,l
mag=mag*sqrt(eta**2.d0+i**2.d0)
enddo
endif

end subroutine cwf_gamma_mag
! Computes natural logarithm of Gamma function of
! double complex argument ln[\Gamma(z)].
!
! 2013-07-15 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration
!
!
!
!       prof. Charlotte Elster: this routine seems to originate in
!           CERNLib.
!
!       Initial version of the code is published at
!       S. Zhang & JM. Jin "Computation of Special Functions"
!       (Wiley, 1996).
!       @book{zhang1996computation,
!       title={Computation of special functions},
!       author={Zhang, S. and Jin, J.M.},
!       lccn={95049544},
!       series={A Wiley-Interscience Publication},
!       url={http://books.google.com/books?id=ASfvAAAAMAAJ},
!       year={1996},
!       publisher={John Wiley}
!       }
!
!       But it was *completely refactored* by Vasily Eremenko.
!
!
!
!
!        See http://dlmf.nist.gov/5 Digital Library of
!        Mathematical Functions (DLMF NIST), Chapter 5,
!        or Abramowitz-Stegun, 1964, pages 256--293
!
!
!        See testing results at the end of file.
!
!
!     This code in its core uses Asymptotic Expansion from
!     A-S, eq.(6.1.40), or
!     DLMF, eq.(5.11.1), http://dlmf.nist.gov/5.11.E1 .
!
!     S. Zhang & JM. Jin call this expansion as Stirling's formula,
!     at A-S there is not any names around it,
!     at DLMF the expansion is at the section Poincare-Type Expansions.
!
!     It is unclear, why the authors stated that asymptotic starts
!     at Re(z) >= 7.

subroutine cwfncgammaln(z, value)
    implicit none

    double complex, intent(in) :: z
    double complex, intent(out) :: value

    ! Other variables and even constants
    double precision :: pi

    double precision :: tmp
    double complex :: ctmp
    ! Z with POSitive REal part
    double complex :: zrepos
    integer :: n
    ! zrepos, but the real part is extended to be > 7
    double complex :: z0
    integer :: k


    ! Working...
pi = 4.d0 * datan(1.d0)
    ! Checking the argument's prohibited values
    tmp = dreal(z)
    if (dimag(z).eq.0.0d0 &
            .and. tmp.eq.int(tmp) &
            .and. tmp.le.0.0d0) then
        ! Stop the program!
        write(0, *) "Argument is 0 or negative integer!"
        stop 1
    end if

    ! Checking the sign(Re(z))
    if (tmp.lt.0.0d0) then
        ! Re(z) < 0
        zrepos = dcmplx(1.d0, 0.d0) - z
    else
        zrepos = z
    end if
    ! Re(zrepos) > 0

    ! Constructing z0.
    ! Making x0 > 7
    ! S. Zhang & JM. Jin stated, that Re(z) >= 7 is already
    ! asymptotic region.
    tmp = dreal(zrepos)
    if (tmp.le.7.0d0) then
        n = ceiling(7.d0 - tmp)
        tmp = tmp + n
    endif
    z0 = dcmplx(tmp, dimag(zrepos))

    ! The core computations...
    ! value = ln[Gamma(z0)]
    call cwfncgammalncore(z0, value)
    ! Done!

    ! Way back...
    !
    ! If Re(zrepos) < 7,
    ! ln[Gamma(zrepos)] = ln[Gamma(z0)] - lz1;
    ! lz1 = \Sum_{k = 0}^{n - 1} ln(z_k);
    ! z_k = zrepos + k;
    ! z_{n} = z0.
    ! Sorry for confusion notations.
    ! This formulas are obtained from
    ! A-S, eqs.(6.1.15) and (6.1.16), or
    ! http://dlmf.nist.gov/5.5.E1 .
    tmp = dreal(zrepos)
    if (tmp.le.7.0d0) then
        n = ceiling(7.d0 - tmp)

        ctmp = dcmplx(0.d0, 0.d0)
        do k = 0, (n - 1)
            ctmp = ctmp + log(zrepos + k)
        end do

        value = value - ctmp
    endif
    ! Now value = ln[Gamma(zrepos)].

    ! If Re(z) < 0,
    ! apply reflection formula A-S, eq.(6.1.17) or
    ! http://dlmf.nist.gov/5.5.E3 .
    tmp = dreal(z)
    if (tmp.lt.0.0d0) then
        ! \ln \pi - \ln \sin (\pi z)
        tmp = log(pi)
        ctmp = z * pi
        ctmp = sin(ctmp)
        ctmp = log(ctmp)
        ctmp = tmp - ctmp

        ! - \ln \Gamma (1 - z)
        value = ctmp - value
    endif
    ! Done!
end subroutine cwfncgammaln







! The computational core of cwfncgammaln().
! Computes the ln[Gamma(z)].
! *NO ARGUMENT CHECKING!*
subroutine cwfncgammalncore(z, value)
    implicit none

    double complex, intent(in) :: z
    double complex, intent(out) :: value

    ! Coefficients at asymptotic expansion of ln[Gamma(z)] in case
    ! z -> \infty and |\arg(z)| < \pi.
    ! See DLMF, eq.(5.11.1), http://dlmf.nist.gov/5.11.E1 .
    !
    ! These values are explicitly published by S. Zhang & JM. Jin
    ! "Computation of Special Functions" (Wiley, 1996).
    ! Every A(i) could be computed using A-S, eq.(6.1.40) or
    ! http://dlmf.nist.gov/5.11.E1
    ! and tables A-S, page 810.
    ! This coefficients are:
    ! a_k = B_{2k} / [2 k (2 k - 1)];
    ! PLEASE, MODIFY THE EXPANSION CORRESPONDINGLY!
    double precision, parameter, dimension(10) :: a = (/ &
          8.333333333333333d-02  &
        , -2.777777777777778d-03 &
        , 7.936507936507937d-04  &
        , -5.952380952380952d-04 &
        , 8.417508417508418d-04  &
        , -1.917526917526918d-03 &
        , 6.410256410256410d-03  &
        , -2.955065359477124d-02 &
        , 1.796443723688307d-01  &
        , -1.39243221690590d+00  &
    /)

    ! Other variables and even constants
    double precision :: pi

    double complex :: ctmp
    double precision :: tmp
    integer :: k
    integer :: n

    ! Working...
pi = 4.d0 * datan(1.d0)
    ! Calculate ln[Gamma(z)] using asymptotic expansion
    ! A-S, eq.(6.1.40) or http://dlmf.nist.gov/5.11.E1,

    ! The first term
    ctmp = dcmplx(0.5d0, 0.d0)
    ctmp = (z - ctmp)
    value = ctmp * log(z)

    ! The second term
    value = value - z

    ! And the third one
    tmp = 0.5d0 * log(2.0d0 * pi)
    ctmp = dcmplx(tmp, 0.d0)
    value = value + ctmp

    ! The sum from the formula
    do k = 1, 10
        ctmp = z**(2 * k - 1)
        ctmp = a(k) / ctmp
        value = value + ctmp
    end do
    ! Done!
end subroutine cwfncgammalncore







! ***
! TESTING RESULTS
!
!
!
! Vasily Eremenko <to@vsl.name>
! TORUS collaboration
!
! 2013-07-15
!
!
!
! The comparison is performed by rejudge utility with default tolerance
! 1.e-10 (related discrepancy).
!
!
!
! TESTS
!
! 1. Agree with tables from Abramowitz-Stegan, 1964, pp.267--287,
!    the relative tolerance mentioned above. All the set of the real
!    argument Gamma(x) testing points, plus 24 complex points were picked
!    up from tables and value of Gamma(z) tested in this points.
!
! 2. Agree with independently generated values with the same related
!    tolerance for the set of points from square (-40 - i 90) to
!    (150 + i 100), all the grid points have reasonably `big' imaginary
!    part.
!  *WARNING:* In case of `small' argument's imaginary part
!    (|Im(z)| ~ 10^{-26} or so), the results, generated by this code,
!    differ from independently generated, but
!    (a)there are not any referees in this area,
!    (b)these values of the argument should be treated as real values
!       in our case.
!
!
!
! Platform/compiler independence
!
! All the tests
! are passed with the same tolerance on 32-bit machine (f95 only) and
! 64-bit machines (compilers: f95, pgf90, ifort, gfortran) with options:
! <no options>, -O2, -O4, -C -g.
!
!
!
!
!
! OVERALL OUTCOME
!
! This code works for us (TORUS Collaboration, Few-Body Group, AGS Project).
! Subroutine cwfnchi(p, q, value)
! Computes the linear transformed \chi(p, q) = (1 - \zeta(p, q))
! 4-th argument
! of the hypergeometric 2F1()
! for the `regular' representation of Coulomb WF.
! See Eq.(3.12) at libcwfn-doc.pdf.
!
! This subroutine is computing \chi(p, q) = (1 - \zeta(p, q))
! with full precision, even in case \zeta(p, q) is limiting to 1.
!
! 2013-08-01 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnchi(p, q, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    double precision, intent(out) :: value

    ! Other variables
    double precision :: p2
    double precision :: q2
    double precision :: denom

    ! Working...
    p2 = p * p
    q2 = q * q

    ! \chi(p, q) = [(p^2 - q^2) / (p^2 + q^2)]^2
    value = p2 - q2
    denom = p2 + q2
    value = value / denom
    value = value * value
end subroutine cwfnchi
! Computes hypergeometric function 2F1(a, b; c; z) of complex arguments.
!
! This is a simple wrapper around the "provider" (the actual code to compute
! hypergeometric function).
!
! The wrapper is used to prevent possible accuracy degradation and to simplify
! the possible switch from one "provider" to another.
!
! 2014-08-21 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnf21(aa,bb,cc,zz,fun)
    implicit none

    double complex, intent(in) :: aa
    double complex, intent(in) :: bb
    double complex, intent(in) :: cc
    double complex, intent(in) :: zz
    double complex, intent(out) :: fun

    interface
        function cwfnMiSt2F1(a, b, c, z)
            double complex, intent(in) :: a
            double complex, intent(in) :: b
            double complex, intent(in) :: c
            double complex, intent(in) :: z
            double complex :: cwfnMiSt2F1
        end function cwfnMiSt2F1
    end interface

    ! Working...
    fun = cwfnMiSt2F1(aa,bb,cc,zz)
end subroutine cwfnf21
! Subroutine cwfnFull(p, q, l, eta, value)
! Computes `full' Coulomb WF \phi(p, q, l, \eta) from the Chapt.3
! at libcwfn-doc.pdf.
!
! WARNING: Be careful! It is your responsibility to to understand, what
! `phone numbers' (prof.Charlotte Elster) the code will generate.
!
! 2013-08-10 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

!Charlie: value now defined as real
subroutine cwfnFull(p, q, l, eta, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double complex, intent(out) :: value
!    double precision, intent(out) :: value

    ! Parameters (computed later) to decide the best representation
    double precision:: pMin
    double precision:: pMax

    ! The threshould of |Im[Xi()]/Re[Xi()]|.
    ! We want to have at least 7 significant digits of \Im[\Xi()]
    double precision, parameter :: xiThreshould = 1.d-6

    ! Variables
    double complex :: xi
    double precision :: imxi
    double precision :: ratio


    ! Working...
    ! Parameters
    pMin = 0.3d0 * q
    pMax = 3.d0 * q

    if(p.lt.pMin .or. p.gt.pMax) then
        ! Non-pole representation
        call cwfnfullnonpole(p, q, l, eta, value)
        !print*,p,q,'nonpole',value
    else
        ! ratio = |Im[Xi()]/Re[Xi()]|
        call cwfnxi(p, q, l, eta, xi)
        !print*,xi
        ratio = dimag(xi)
        ratio = ratio / dreal(xi)
        ratio = abs(ratio)
        if(ratio.le.xiThreshould) then
            ! Non-pole representation
            call cwfnfullnonpole(p, q, l, eta, value)
        !print*,p,q,'nonpole',value
        else
            ! Pole representation
            imxi = dimag(xi)
            call cwfnfullpolecached(p, q, l, eta, imxi, value)
        end if
    end if
end subroutine cwfnFull
! Subroutine cwfnfullnonpole(p, q, l, eta, value)
! Computes \Psi(p, q, l, \eta) from Eq.(3.6) at libcwfn-doc.pdf.
!
! 2013-07-11 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

!Charlie: Value declared as real
subroutine cwfnfullnonpole(p, q, l, eta, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double complex, intent(out) :: value
!    double precision, intent(out) :: value
    ! Other variables
!    double complex :: upsilon
!    double complex :: omega
    double precision :: upsilon
    double precision :: omega

    ! Working...
    call cwfnupsilon(q, l, eta, upsilon)
    call cwfnomega(p, q, l, eta, omega)
     !print*,q,p,upsilon,omega
    value = upsilon * omega
end subroutine cwfnfullnonpole
! Subroutine cwfnfullpole(p, q, l, eta, value)
! Computes \Phi(p, q, l, \eta) from the eq.(3.18) at libcwfn-doc.pdf.
!
! 2013-07-11 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnfullpole(p, q, l, eta, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
!    double complex, intent(out) :: value
    double precision, intent(out) :: value
    ! Other variables
!    double complex :: theta
!    double complex :: lambda
     double precision :: theta, lambda

    ! Working...
    call cwfntheta(q, l, eta, theta)
    call cwfnlambda(p, q, l, eta, lambda)
    value = theta * lambda
end subroutine cwfnfullpole
! Subroutine cwfnfullpolecached(p, q, l, eta, imxi, value)
! Computes \Phi(p, q, l, \eta) from the eq.(3.18) at libcwfn-doc.pdf
! *USING* cwfnlambdacached().
!
! 2013-07-11 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

!Charlie
subroutine cwfnfullpolecached(p, q, l, eta, imxi, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double precision, intent(in) :: imxi
    double complex, intent(out) :: value
!    double precision, intent(out) :: value
    ! Other variables
!    double complex :: theta
!    double complex :: lambda
    double precision :: theta
    double precision :: lambda

    ! Working...
    call cwfntheta(q, l, eta, theta)
    call cwfnlambdacached(p, q, l, eta, imxi, lambda)
    value = theta * lambda
end subroutine cwfnfullpolecached
! Subroutine cwfnFullV(n, pv, q, l, eta, values)
! Computes `full' Coulomb WF \phi(pv(i), q, l, eta)
! for every i = 1, 2,... n;
! results are stored at values.
!
! 2014-02-20 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnFullV(n, pv, q, l, eta, values)
    implicit none

    ! I/O units
    integer, parameter :: stderr = 0

    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: pv
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double complex, dimension(n), intent(out) :: values

    ! Other variables
    integer :: i
    double precision :: p
    double complex :: val


    ! Working...
    do i = 1, n
        p = pv(i)
        call cwfnfull(p, q, l, eta, val)
        values(i) = val
    end do
end subroutine cwfnFullV
! Subroutine cwfnlambda(p, q, l, eta, value)
! Computes \Lambda(p, q, l, \eta) from the eq.(3.20) at libcwfn-doc.pdf.
!
! 2014-05-15 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnlambda(p, q, l, eta, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
!    double complex, intent(out) :: value
    double precision, intent(out) :: value

    ! Other variables
    double precision :: tmp
    double complex :: xi

    ! Working...
    ! Computing Xi(p, q, l, eta)
    call cwfnxi(p, q, l, eta, xi)

    ! Other stuff
    tmp = (p + q)**(l * 2 - 1)
    tmp = tmp / p**(l + 1)
    tmp = tmp * dimag(xi)
    value = dcmplx(tmp, 0.d0)
end subroutine cwfnlambda
! Subroutine cwfnlambdacached(p, q, l, eta, imxi, value)
! Computes \Lambda(p, q, l, \eta) from Eq.(3.20) at libcwfn-doc.pdf
! *USING PROVIDED* Im[Xi()] value.
!
! 2014-05-15 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

!Charlie
subroutine cwfnlambdacached(p, q, l, eta, imxi, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double precision, intent(in) :: imxi
!    double complex, intent(out) :: value
    double precision, intent(out) :: value
    ! Other variables
    double precision :: tmp

    ! Working...
    tmp = (p + q)**(l * 2 - 1)
    tmp = tmp / p**(l + 1)
    tmp = tmp * imxi
!    value = dcmplx(tmp, 0.d0)
    value=tmp
end subroutine cwfnlambdacached
! Everything required to compute the hypergeometric function 2F1(a, b, c; z)
! with complex arguments.
!
! The initial code was taken from http://cpc.cs.qub.ac.uk/summaries/AEAE
! [N.Michel, M.V.Stoitsov. // Comput. Phys. Commun. 178 (2008) 535].
!
! The actual code was slightly modified by Vasily Eremenko <to@vsl.name> in
! the name of TORUS Collaboration ( http://reactiontheory.org ).  The
! modifications are about to rename the module and functions, plus few other
! minor changes to increase the code portability and prevent naming confusion.
! All the comments below "***" are preserved as it is.
!
! 2014-08-21 Vasily Eremenko <to@vsl.name> in the name of TORUS Collaboration
!
! ***
!
! Gamma_inv denotes the entire inverse of the Gamma function.
! F(z) means 2F1(a,b,c,z) with the a, b, c and z given as inputs 
! in the routine.
!
! Elementary functions and standard constants 
! are defined in the module.
! See N.J.~Higham, ``Accuracy and Stability of Numerical Algorithms'',
! SIAM, Philadelphia, 1996 for expm1 implementation.
! log1p follows instantly.
!
!----------------------------------------------------------------------
MODULE CWFN_MIST
  !--------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, PARAMETER :: PR=KIND(1.0D0),IPR=KIND(1)
  REAL(PR)     :: EPS15=1.0D-15
  REAL(PR)     :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.50D0
  REAL(PR)     :: M_PI=3.14159265358979323846D0
  REAL(PR)     :: M_PI_2=1.57079632679489661923D0
  REAL(PR)     :: M_1_PI=0.31830988618379067154D0 
CONTAINS
  !
  FUNCTION INF_NORM(Z)
    COMPLEX(PR),INTENT(IN) :: Z
    REAL(PR)  :: INF_NORM
    INF_NORM=MAX(ABS(REAL(Z,PR)),ABS(AIMAG(Z)))
    RETURN
  END FUNCTION INF_NORM
  !
  FUNCTION TANZ(Z)
    COMPLEX(PR),INTENT(IN) :: Z
    COMPLEX(PR) :: TANZ
    TANZ=SIN(Z)/COS(Z)
    RETURN
  END FUNCTION TANZ
  !
  FUNCTION LOG1P(Z)
    COMPLEX(PR),INTENT(IN) :: Z
    REAL(PR) :: X,XP1,LOG1P_X
    REAL(PR) :: Y,YX,YX2,YX2P1,LOG1P_YX2
    REAL(PR) :: RE_LOG1P,IM_LOG1P
    COMPLEX(PR) :: LOG1P
    IF(INF_NORM(Z).LT.ONE) THEN
       X = REAL(Z,PR); XP1 = X+ONE
       IF(XP1.EQ.ONE) THEN
          LOG1P_X = X
       ELSE
          LOG1P_X = LOG(XP1)*X/(XP1-ONE)
       ENDIF
       Y = AIMAG(Z)
       YX = Y/XP1; YX2 = YX*YX; YX2P1 = YX2+ONE
       IF(YX2P1.EQ.ONE) THEN
          LOG1P_YX2 = YX2
       ELSE
          LOG1P_YX2 = LOG(YX2P1)*YX2/(YX2P1-ONE)
       ENDIF
       RE_LOG1P = LOG1P_X + HALF*LOG1P_YX2
       IM_LOG1P = ATAN2(Y,XP1)
       LOG1P = CMPLX(RE_LOG1P,IM_LOG1P,PR)
       RETURN
    ELSE
       LOG1P=LOG(ONE+Z)
       RETURN
    ENDIF
  END FUNCTION LOG1P
  !
  FUNCTION EXPM1(Z)
    COMPLEX(PR),INTENT(IN) :: Z
    REAL(PR) :: X,EXPM1_X,EXP_X,Y,SIN_HALF_Y
    REAL(PR) :: RE_EXPM1,IM_EXPM1
    COMPLEX(PR) :: EXPM1
    IF(INF_NORM(Z).LT.ONE) THEN
       X = REAL(Z,PR); EXP_X = EXP(X)
       Y = AIMAG(Z); SIN_HALF_Y=SIN(HALF*Y)
       IF(EXP_X.EQ.ONE) THEN
          EXPM1_X = X
       ELSE 
          EXPM1_X = (EXP_X-ONE)*X/LOG(EXP_X)
       ENDIF
       RE_EXPM1 = EXPM1_X-TWO*EXP_X*SIN_HALF_Y*SIN_HALF_Y 
       IM_EXPM1 = EXP_X*SIN(Y)
       EXPM1 = CMPLX(RE_EXPM1,IM_EXPM1,PR)
       RETURN
    ELSE
       EXPM1=EXP(Z)-ONE
       RETURN
    ENDIF
  END FUNCTION EXPM1
  !
END MODULE CWFN_MIST
!
!----------------------------------------------------------------------
RECURSIVE FUNCTION CWFN_MIST_LOG_GAMMA(Z) RESULT(RES)
!----------------------------------------------------------------------
! Logarithm of Gamma[z] and Gamma inverse function
! ------------------------------------------------
!
! For log[Gamma[z]],if z is not finite 
! or is a negative integer, the program 
! returns an error message and stops.
! The Lanczos method is used. Precision : ~ 1E-15
! The method works for Re[z]>0.5 .
! If Re[z]<=0.5, one uses the formula Gamma[z].Gamma[1-z]=Pi/sin(Pi.z)
! log[sin(Pi.z)] is calculated with the Kolbig method 
! (K.S. Kolbig, Comp. Phys. Comm., Vol. 4, p.221(1972)): 
! If z=x+iy and y>=0, log[sin(Pi.z)]=log[sin(Pi.eps)]-i.Pi.n, 
! with z=n+eps so 0<=Re[eps]< 1 and n integer.
! If y>110, log[sin(Pi.z)]=-i.Pi.z+log[0.5]+i.Pi/2 
! numerically so that no overflow can occur.
! If z=x+iy and y< 0, log[Gamma(z)]=[log[Gamma(z*)]]*, 
! so that one can use the previous formula with z*.
!
! For Gamma inverse, Lanczos method is also used 
! with Euler reflection formula.
! sin (Pi.z) is calculated as sin (Pi.(z-n)) 
! to avoid inaccuracy with z = n + eps 
! with n integer and |eps| as small as possible.
!
!
! Variables:
! ----------
! x,y: Re[z], Im[z]
! log_sqrt_2Pi,log_Pi : log[sqrt(2.Pi)], log(Pi).
! sum : Rational function in the Lanczos method
! log_Gamma_z : log[Gamma(z)] value.
! c : table containing the fifteen coefficients in the expansion 
! used in the Lanczos method.
! eps,n : z=n+eps so 0<=Re[eps]< 1 and n integer for Log[Gamma].
!         z=n+eps and n integer 
!         so |eps| is as small as possible for Gamma_inv.
! log_const : log[0.5]+i.Pi/2
! g : coefficient used in the Lanczos formula. It is here 607/128.
! z,z_m_0p5,z_p_g_m0p5,zm1 : argument of the Gamma function, 
! z-0.5, z-0.5+g, z-1 
! res: returned value
!----------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: Z
  INTEGER(IPR) :: N,I
  REAL(PR)     :: X,Y,LOG_SQRT_2PI,G,LOG_PI,M_LN2,C(0:14)
  COMPLEX(PR)  :: GAMMA_SUM,Z_M_0P5,Z_P_G_M0P5,ZM1
  COMPLEX(PR)  :: LOG_CONST,I_PI,EPS,LOG_SIN_PI_Z,RES
  !
  M_LN2=0.69314718055994530942D0; X=REAL(Z,PR); Y=AIMAG(Z)
  IF((Z.EQ.NINT(X)).AND.(X.LE.ZERO)) &
       STOP 'Z IS NEGATIVE INTEGER IN CWFN_MIST_LOG_GAMMA'
  IF(X.GE.HALF) THEN
     LOG_SQRT_2PI=0.91893853320467274177D0; G=4.7421875D0
     Z_M_0P5=Z-HALF; Z_P_G_M0P5=Z_M_0P5+G; ZM1=Z-ONE
     C=(/ 0.99999999999999709182D0,57.156235665862923517D0,       &
          -59.597960355475491248D0,  14.136097974741747174D0,     &
          -0.49191381609762019978D0, 0.33994649984811888699D-4,   &
          0.46523628927048575665D-4, -0.98374475304879564677D-4,  &
          0.15808870322491248884D-3, -0.21026444172410488319D-3,  &
          0.21743961811521264320D-3, -0.16431810653676389022D-3,  &
          0.84418223983852743293D-4, -0.26190838401581408670D-4,  &
          0.36899182659531622704D-5 /)

     GAMMA_SUM=C(0)
     DO I=1,14
        GAMMA_SUM=GAMMA_SUM+C(I)/(ZM1+I)
     ENDDO
     RES=LOG_SQRT_2PI+LOG(GAMMA_SUM)+Z_M_0P5*LOG(Z_P_G_M0P5) &
          -Z_P_G_M0P5
     RETURN
  ELSE IF(Y.GE.ZERO) THEN
     IF(X.LT.NINT(X)) THEN
        N=NINT(X)-1
     ELSE
        N=NINT(X)
     ENDIF
     LOG_PI=1.1447298858494002D0
     LOG_CONST=CMPLX(-M_LN2,M_PI_2,PR); I_PI=CMPLX(ZERO,M_PI,PR)
     EPS=Z-N
     IF(Y.GT.110.0D0) THEN
        LOG_SIN_PI_Z=-I_PI*Z+LOG_CONST
     ELSE
        LOG_SIN_PI_Z=LOG(SIN(M_PI*EPS))-I_PI*N
     ENDIF
     RES=LOG_PI-LOG_SIN_PI_Z-CWFN_MIST_LOG_GAMMA(ONE-Z);
     RETURN
  ELSE
     RES=CONJG(CWFN_MIST_LOG_GAMMA(CONJG(Z)))
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_LOG_GAMMA
!
!----------------------------------------------------------------------
! Inverse of the Gamma function [1/Gamma](z)
! ------------------------------------------
! It is calculated with the Lanczos method for Re[z] >= 0.5 
! and is precise up to 10^{-15}.
! If Re[z] <= 0.5, one uses the formula 
! Gamma[z].Gamma[1-z] = Pi/sin (Pi.z).
! sin (Pi.z) is calculated as sin (Pi.(z-n)) to avoid inaccuracy,
! with z = n + eps with n integer and |eps| as small as possible.
! 
! Variables 
! ---------
! z : argument of the function
! x: Re[z]
! eps,n : z = n + eps with n integer and |eps| as small as possible.
! res: returned value
!----------------------------------------------------------------------
RECURSIVE FUNCTION CWFN_MIST_GAMMA_INV(Z) RESULT(RES)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: Z
  INTEGER(IPR) :: N,I
  REAL(PR)     :: X,LOG_SQRT_2PI,G,C(0:14)
  COMPLEX(PR)  :: RES,GAMMA_SUM,Z_M_0P5,Z_P_G_M0P5,ZM1,EPS
  !
  X=REAL(Z,PR)
  IF(X.GE.HALF) THEN
     LOG_SQRT_2PI=0.91893853320467274177D0; G=4.7421875D0
     Z_M_0P5=Z-HALF; Z_P_G_M0P5=Z_M_0P5+G; ZM1=Z-ONE
     C=(/ 0.99999999999999709182D0,57.156235665862923517D0,       &
          -59.597960355475491248D0,  14.136097974741747174D0,     &
          -0.49191381609762019978D0, 0.33994649984811888699D-4,   &
          0.46523628927048575665D-4, -0.98374475304879564677D-4,  &
          0.15808870322491248884D-3, -0.21026444172410488319D-3,  &
          0.21743961811521264320D-3, -0.16431810653676389022D-3,  &
          0.84418223983852743293D-4, -0.26190838401581408670D-4,  &
          0.36899182659531622704D-5 /)

     GAMMA_SUM=C(0)
     DO I=1,14
        GAMMA_SUM=GAMMA_SUM+C(I)/(ZM1+I);
     ENDDO
     RES=EXP(Z_P_G_M0P5-Z_M_0P5*LOG(Z_P_G_M0P5)-LOG_SQRT_2PI) &
          /GAMMA_SUM
     RETURN
  ELSE
     X=REAL(Z,PR); N=NINT(X)
     EPS=Z-N
     IF(MOD(N,2).EQ.0) THEN
        RES=SIN(M_PI*EPS)*M_1_PI/CWFN_MIST_GAMMA_INV (ONE-Z)
        RETURN
     ELSE
        RES=-SIN(M_PI*EPS)*M_1_PI/CWFN_MIST_GAMMA_INV (ONE-Z)
        RETURN
     ENDIF
  ENDIF
END FUNCTION CWFN_MIST_GAMMA_INV
!----------------------------------------------------------------------
!
! Calculation of H(z,eps) = [Gamma(z+eps)/Gamma(z) - 1]/eps, with e and
! ---------------------------------------------------------------------
! z complex so z,z+eps are not negative integers and 0 <= |eps|oo < 0.1
! ---------------------------------------------------------------------
! The function H(z,eps) = [Gamma(z+eps)/Gamma(z) - 1]/e is calculated 
! here with the Lanczos method.
! For the Lanczos method, the gamma parameter, denoted as g, 
! is 4.7421875 and one uses a sum of 15 numbers with the table c[15], 
! so that it is precise up to machine accuracy.
! The H(z,eps) function is used in formulas occuring in1-z and 1/z 
! transformations (see Comp. Phys. Comm. paper).
!
! One must have z and z+eps not negative integers as otherwise 
! it is clearly not defined.
! As this function is meant to be precise for small |eps|oo, 
! one has to have 0 <= |eps|oo < 0.1 .
! Indeed, a direct implementation of H(z,eps) with Gamma_inv or 
! log_Gamma for |eps|oo >= 0.1 is numerically stable.
! The returned function has full numerical accuracy 
! even if |eps|oo is very small.
!
! eps not equal to zero
! ---------------------
! If Re(z) >= 0.5 or Re(z+eps) >= 0.5, one clearly has Re(z) > 0.4 
! and Re(z+eps) > 0.4, 
! so that the Lanczos summation can be used for both Gamma(z) 
! and Gamma(z+eps).
! One then has:
! log[Gamma(z+eps)/Gamma(z)] = 
! (z-0.5) log1p[eps/(z+g-0.5)] + eps log(z+g-0.5+eps) - eps 
! + log1p[-eps \sum_{i=1}^{14} c[i]/((z-1+i)(z-1+i+eps)) 
! / (c[0] + \sum_{i=1}^{14} c[i]/(z-1+i))]
! H(z,eps) = expm1[log[Gamma(z+eps)/Gamma(z)]]/eps .
!
! If Re(z) < 0.5 and Re(z+eps) < 0.5, 
! Euler reflection formula is used for both Gamma(z) and Gamma(z+eps).
! One then has: 
! H(z+eps,-eps) = [cos(pi.eps) + sin(pi.eps)/tan(pi(z-n))].H(1-z,-eps) 
! + (2/eps).sin^2(eps.pi/2) - sin(pi.eps)/(eps.tan(pi.(z-n)))
! H(1-z,-eps) is calculated with the Lanczos summation 
! as Re(1-z) >= 0.5 and Re(1-z-eps) >= 0.5 .
! z-n is used in tan(pi.z) instead of z to avoid inaccuracies 
! due the finite number of digits of pi.
! H(z,eps) = H(z+eps,-eps)/(1 - eps.H(z+eps,-eps)) 
! provides the final result.
!
! eps equal to zero
! -----------------
! It is obtained with the previous case and eps -> 0 :
! If Re(z) >= 0.5, one has:
! H(z,eps) = (z-0.5)/(z+g-0.5) + log(z+g-0.5) - 1 -
! \sum_{i=1}^{14} c[i]/((z-1+i)^2)/(c[0]+\sum_{i=1}^{14} c[i]/(z-1+i))
!
! If Re(z) < 0.5, one has:
! H(z,0) = H(1-z,0) - pi/tan(pi.(z-n))
!
! Variables
! ---------
! z,eps: input variables of the function H(z,eps)
! g,c[15]: double and table of 15 doubles defining the Lanczos sum 
! so that it provides the Gamma function 
! precise up to machine accuracy.
! eps_pz,z_m_0p5,z_pg_m0p5,eps_pz_pg_m0p5,zm1,zm1_p_eps: 
! z+eps,z-0.5,z+g-0.5,z+eps+g-0.5,z-1,z-1+eps
! x,eps_px: real parts of z and z+eps.
! n,m: closest integer ot the real part of z, same for z+eps.
! sum_num,sum_den: \sum_{i=1}^{14} c[i]/((z-1+i)(z-1+i+eps)) 
! and (c[0] + \sum_{i=1}^{14} c[i]/(z-1+i)). 
! They appear respectively as numerator and denominator in formulas.
! Pi_eps,term,T1_eps_z: pi.eps, sin (pi.eps)/tan(pi.(z-n)), 
! [cos(pi.eps) + sin(pi.eps)/tan(pi(z-n))].H(1-z,-eps)
! sin_Pi_2_eps,T2_eps_z,T_eps_z: sin^2(eps.pi/2), 
! (2/eps).sin^2(eps.pi/2) - sin(pi.eps)/(eps.tan(pi.(z-n))), 
! H(z+eps,-eps)
! res: returned value
!----------------------------------------------------------------------
RECURSIVE FUNCTION CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS(Z,EPS) RESULT(RES)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: Z,EPS
  INTEGER(IPR) :: N,M,I
  REAL(PR)     :: G,X,EPS_PX,C(0:14)
  COMPLEX(PR)  :: RES,SUM_NUM,SUM_DEN
  COMPLEX(PR)  :: EPS_PZ,Z_M_0P5,Z_PG_M0P5,EPS_PZ_PG_M0P5,ZM1
  COMPLEX(PR)  :: CI_ZM1_PI_INV,PI_EPS,TT,T1_EPS_Z,SIN_PI_2_EPS
  COMPLEX(PR)  :: ZM1_P_EPS,T2_EPS_Z,T_EPS_Z
  !
  G=4.74218750D0
  IF(INF_NORM(EPS).GT.0.1D0) &
       STOP 'ONE MUST HAVE |EPS|< 0.1 IN CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS'
  EPS_PZ=Z+EPS; Z_M_0P5=Z-HALF; Z_PG_M0P5=Z_M_0P5+G
  EPS_PZ_PG_M0P5=Z_PG_M0P5+EPS; ZM1=Z-ONE; ZM1_P_EPS=ZM1+EPS
  X=REAL(Z,PR); EPS_PX=REAL(EPS_PZ,PR); N=NINT(X); M=NINT(EPS_PX)
  IF((Z.EQ.N).AND.(N.LE.0)) THEN
     STOP 'Z IS NEGATIVE INTEGER IN CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS'
  ENDIF
  IF((EPS_PZ.EQ.M).AND.(M.LE.0)) THEN
     STOP 'Z+EPS IS NEGATIVE INTEGER IN CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS'
  ENDIF
  C=(/ 0.99999999999999709182D0,57.156235665862923517D0,     &
       -59.597960355475491248D0,14.136097974741747174D0,     &
       -0.49191381609762019978D0,0.33994649984811888699D-4,  &
       0.46523628927048575665D-4,-0.98374475304879564677D-4, &
       0.15808870322491248884D-3,-0.21026444172410488319D-3, &
       0.21743961811521264320D-3,-0.16431810653676389022D-3, &
       0.84418223983852743293D-4,-0.26190838401581408670D-4, &
       0.36899182659531622704D-5 /)
  IF((X.GE.HALF).OR.(EPS_PX.GE.HALF)) THEN
     SUM_NUM=ZERO;SUM_DEN=C(0)
     DO I=1,14
        CI_ZM1_PI_INV=C(I)/(ZM1+I)
        SUM_NUM=SUM_NUM+CI_ZM1_PI_INV/(ZM1_P_EPS+I)
        SUM_DEN=SUM_DEN+CI_ZM1_PI_INV
     ENDDO
     IF(EPS.NE.ZERO) THEN
        RES=EXPM1(Z_M_0P5*LOG1P(EPS/Z_PG_M0P5) &
             +EPS*LOG(EPS_PZ_PG_M0P5)-EPS+LOG1P(-EPS*SUM_NUM/SUM_DEN))&
             /EPS
        RETURN
     ELSE
        RES=Z_M_0P5/Z_PG_M0P5 &
             +LOG(EPS_PZ_PG_M0P5)-ONE-SUM_NUM/SUM_DEN
        RETURN
     ENDIF
  ELSE
     IF(EPS.NE.ZERO) THEN
        PI_EPS=M_PI*EPS
        TT=SIN(PI_EPS)/TANZ(M_PI*(Z-N))
        T1_EPS_Z=(COS(PI_EPS)+TT)*& 
             CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS(ONE-Z,-EPS)
        SIN_PI_2_EPS=SIN(M_PI_2*EPS)
        T2_EPS_Z=(TWO*SIN_PI_2_EPS*SIN_PI_2_EPS-TT)/EPS
        T_EPS_Z=T1_EPS_Z+T2_EPS_Z
        RES=(T_EPS_Z/(ONE-EPS*T_EPS_Z))
        RETURN
     ELSE
        RES=CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS(ONE-Z,-EPS) &
             -M_PI/TANZ(M_PI*(Z-N))
        RETURN
     ENDIF
  ENDIF
END FUNCTION CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS
!
!----------------------------------------------------------------------
! Calculation of G(z,eps) = [Gamma_inv(z) - Gamma_inv(z+eps)]/eps 
! ---------------------------------------------------------------
! with e and z complex
!---------------------
! The G(z,eps) function is used in formulas occuring in 1-z 
! and 1/z transformations (see Comp. Phys. Comm. paper).
! Several case have to be considered for its evaluation. 
! eps is considered equal to zero 
! if z+eps and z are equal numerically.
!
! |eps|oo > 0.1
! -------------
! A direct evaluation with the values Gamma_inv(z) 
! and Gamma_inv(z+eps) is stable and returned.
!
! |eps|oo <= 0.1 with z+eps and z numerically different
! -----------------------------------------------------
! If z is a negative integer, z+eps is not, 
! so that G(z,eps) = -Gamma_inv(z+eps)/eps, 
! for which a direct evaluation is precise and returned.
! If z+eps is a negative integer, z is not, 
! so that G(z,eps) = Gamma_inv(z)/eps, 
! for which a direct evaluation is precise and returned.
! If both of them are not negative integers, 
! one looks for the one of z and z+eps 
! which is the closest to a negative integer.
! If it is z, one returns H(z,eps).Gamma_inv(z+eps). 
! If it is z+eps, one returns H(z+eps,-eps).Gamma_inv(z).
! Both values are equal, so that one chooses the one 
! which makes the Gamma ratio Gamma(z+eps)/Gamma(z) 
! in H(z,eps) the smallest in modulus.
!
! z+eps and z numerically equal
! -----------------------------
! If z is negative integer, G(z,0) = (-1)^(n+1) n!, 
! where z = -n, n integer, which is returned.
! If z is not negative integer, one returns H(z,eps).Gamma_inv(z+eps)
!
! Variables
! ---------
! z,eps: input variables of the function G(z,eps)
! eps_pz,x,eps_px: z+eps,real parts of z and z+eps.
! n,m: closest integer ot the real part of z, same for z+eps.
! fact,k: (-1)^(n+1) n!, returned when z = -n, n integer 
! and z and z+eps identical numerically (eps ~ 0). 
! It is calculated with integer index k.
! is_z_negative_integer,is_eps_pz_negative_integer: 
! true if z is a negative integer, false if not, same for z+eps.
! z_neg_int_distance, eps_pz_neg_int_distance: 
! |z + |n||oo, |z + eps + |m||oo. 
! If |z + |n||oo < |z + eps + |m||oo, 
! z is closer to the set of negative integers than z+eps.
! Gamma_inv(z+eps) is then of moderate modulus 
! if Gamma_inv(z) is very small. 
! If z ~ n, H(z,eps) ~ -1/eps, 
! that so returning 
! G(z,eps) = H(z,eps).Gamma_inv(z+eps) here is preferred.
! Same for |z + |n||oo > |z + eps + |m||oo with z <-> z+eps.
!
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_GAMMA_INV_DIFF_EPS(Z,EPS)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: Z,EPS
  INTEGER(IPR) :: M,N,K
  REAL(PR)     :: X,EPS_PX,FACT
  REAL(PR)     :: Z_NEG_INT_DISTANCE
  REAL(PR)     :: EPS_PZ_NEG_INT_DISTANCE
  COMPLEX(PR)  :: CWFN_MIST_GAMMA_INV_DIFF_EPS,EPS_PZ,CWFN_MIST_GAMMA_INV
  COMPLEX(PR)  :: CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS
  LOGICAL      :: IS_Z_NEG_INT,IS_EPS_PZ_NEG_INT

  EPS_PZ=Z+EPS; X=REAL(Z,PR); EPS_PX=REAL(EPS_PZ,PR)
  N=NINT(X); M=NINT(EPS_PX)
  IS_Z_NEG_INT=(Z.EQ.N).AND.(N.LE.0)
  IS_EPS_PZ_NEG_INT=(EPS_PZ.EQ.M).AND.(M.LE.0)
  IF(INF_NORM(EPS).GT.0.10D0) THEN
     CWFN_MIST_GAMMA_INV_DIFF_EPS = (CWFN_MIST_GAMMA_INV (Z) - CWFN_MIST_GAMMA_INV (EPS_PZ))/EPS
     RETURN
  ELSE IF(EPS_PZ.NE.Z) THEN 
     IF(IS_Z_NEG_INT) THEN
        CWFN_MIST_GAMMA_INV_DIFF_EPS = (-CWFN_MIST_GAMMA_INV (EPS_PZ)/EPS)
        RETURN
     ELSE IF(IS_EPS_PZ_NEG_INT) THEN
        CWFN_MIST_GAMMA_INV_DIFF_EPS = (CWFN_MIST_GAMMA_INV (Z)/EPS)
        RETURN
     ELSE
        Z_NEG_INT_DISTANCE = INF_NORM (Z + ABS (N))
        EPS_PZ_NEG_INT_DISTANCE = INF_NORM (EPS_PZ + ABS (M))
        IF(Z_NEG_INT_DISTANCE.LT.EPS_PZ_NEG_INT_DISTANCE) THEN
           CWFN_MIST_GAMMA_INV_DIFF_EPS= &
                CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS (Z,EPS)*CWFN_MIST_GAMMA_INV (EPS_PZ)
           RETURN
        ELSE
           CWFN_MIST_GAMMA_INV_DIFF_EPS= &
                CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS (EPS_PZ,-EPS)*CWFN_MIST_GAMMA_INV (Z)
           RETURN
        ENDIF
     ENDIF
  ELSE IF(IS_Z_NEG_INT.AND.IS_EPS_PZ_NEG_INT) THEN
     FACT = -ONE;K=-1
     DO WHILE (K.GE.N) 
        FACT=FACT*K
        K=K-1 
     ENDDO
     CWFN_MIST_GAMMA_INV_DIFF_EPS = FACT
     RETURN
  ELSE
     CWFN_MIST_GAMMA_INV_DIFF_EPS = &
          CWFN_MIST_GAMMA_RATIO_DIFF_SMALL_EPS (Z,EPS)*CWFN_MIST_GAMMA_INV (EPS_PZ)
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_GAMMA_INV_DIFF_EPS
!----------------------------------------------------------------------
!
! Calculation of Gamma_inv(1-m-eps)/eps of the A(z) polynomial in 1-z
! -------------------------------------------------------------------
! and 1/z transformations
! -----------------------
! This value occurs in A(z) in 1-z and 1/z transformations 
! (see Comp. Phys. Comm. paper) for m > 0.
! Both cases of 1-m-eps numerically negative integer 
! or not have to be considered
! 
! 1-eps-m and 1-m numerically different
! -------------------------------------
! One returns Gamma_inv(1-m-eps)/eps directly 
! as its value is accurate.
! To calculate Gamma_inv(1-m-eps), 
! one uses the value Gamma_inv(1-eps), 
! needed in considered transformations,
! and one uses the equality 
! Gamma_inv(1-m-eps) = Gamma_inv(1-eps) \prod_{i=1}^{m} (1-eps-i) 
! for m > 0.
! It is trivially demonstrated 
! from the equality Gamma(x+1) = x.Gamma(x). 
! One Gamma function evaluation is removed this way 
! from the calculation.
! 
! 1-eps-m and 1-m numerically equal
! ---------------------------------
! This implies that 1-m-eps is negative integer numerically.
! Here, eps~0, so that one returns the limit of Gamma_inv(1-m-eps)/eps
! for eps -> 0, which is (-1)^m (m-1)!
!
! Variables
! ---------
! m,eps: variable inputs of the function 
! (m,eps) -> Gamma_inv(1-m-eps)/eps
! Gamma_inv_one_meps: Gamma_inv(1-eps), 
! previously calculated and here recycled 
! to quickly calculate Gamma_inv(1-m-eps).
! one_meps: 1-eps
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_A_SUM_INIT(M,EPS,GAMMA_INV_ONE_MEPS)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  INTEGER(IPR),INTENT(IN) :: M
  COMPLEX(PR),INTENT(IN) :: EPS,GAMMA_INV_ONE_MEPS
  INTEGER(IPR) :: N,I
  REAL(PR)     :: FACT
  COMPLEX(PR)  :: CWFN_MIST_A_SUM_INIT,ONE_MEPS
  COMPLEX(PR)  :: GAMMA_INV_ONE_MEPS_MM
  !
  ONE_MEPS = ONE - EPS
  IF(ONE_MEPS-M.NE.1-M) THEN
     GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS
     DO I=1,M
        GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS_MM*(ONE_MEPS-I)
     ENDDO
     CWFN_MIST_A_SUM_INIT=GAMMA_INV_ONE_MEPS_MM/EPS
     RETURN
  ELSE
     FACT=ONE
     DO N=2,M-1
        FACT=FACT*N
     ENDDO
     IF(MOD(M,2).EQ.0) THEN
        CWFN_MIST_A_SUM_INIT=FACT
     ELSE
        CWFN_MIST_A_SUM_INIT=-FACT
     ENDIF
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_A_SUM_INIT
!
!----------------------------------------------------------------------
! Calculation of the log of Gamma_inv(1-m-eps)/eps
! ------------------------------------------------
! See previous function. 
! It is used in case Gamma_inv(1-m-eps)/eps might overflow.
!
! Variables
! ---------
! m,eps: variable inputs of the function 
! (m,eps) -> log[Gamma_inv(1-m-eps)/eps]
! one_meps_mm: 1-eps-m
! i_Pi: i.Pi
! log_fact: logarithm of (-1)^m (m-1)!, 
! here defined as log((m-1)!) + i.Pi if m is odd.
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_LOG_A_SUM_INIT(M,EPS)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  INTEGER(IPR),INTENT(IN) :: M
  COMPLEX(PR),INTENT(IN) :: EPS
  INTEGER(IPR) :: N
  REAL(PR)     :: LOG_FACT
  COMPLEX(PR)  :: ONE_MEPS_MM,CWFN_MIST_LOG_A_SUM_INIT,CWFN_MIST_LOG_GAMMA
  !
  ONE_MEPS_MM=ONE-EPS-M
  IF(ONE_MEPS_MM.NE.1-M) THEN
     CWFN_MIST_LOG_A_SUM_INIT=(-CWFN_MIST_LOG_GAMMA(ONE_MEPS_MM) - LOG(EPS))
     RETURN
  ELSE
     LOG_FACT=ZERO
     DO N=2,M-1
        LOG_FACT=LOG_FACT + LOG(DBLE(N))
     ENDDO
     IF(MOD(M,2).EQ.0) THEN
        CWFN_MIST_LOG_A_SUM_INIT=LOG_FACT
     ELSE
        CWFN_MIST_LOG_A_SUM_INIT=CMPLX(LOG_FACT,M_PI,PR)
     ENDIF
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_LOG_A_SUM_INIT
!----------------------------------------------------------------------
! Calculation of the first term of the B(z) power series
! ------------------------------------------------------
! in the 1-z transformation, divided by (1-z)^m
! ----------------------------------------------
! In the 1-z transformation, 
! the power series B(z) = \sum_{n=0}^{+oo} \beta_n (1-z)^n occurs 
! (see Comp. Phys. Comm. paper).
! The first term \beta_0, divided by (1-z)^m, is calculated here. 
! m is the closest integer to Re(c-a-b) >= 0 and eps = c-a-b-m.
!
! One has to consider |eps|oo > 0.1 and |eps|oo <= 0.1, 
! where 1-m-eps and 1-m can be different or equal numerically, 
! leading to some changes in this last case.
!
! |eps|oo > 0.1
! -------------
! One has \beta_0/(1-z)^m = [(a)_m (b)_m Gamma_inv(1-eps) 
! Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) Gamma_inv(m+1)
! - (1-z)^eps Gamma_inv(a) Gamma_inv(b) Gamma_inv(1+m+eps)]
! [Gamma(c)/eps], stable in this regime for a direct evaluation.
!
! The values of Gamma(c), Gamma_inv(a+m+eps) 
! and Gamma_inv(b+m+eps) were already calculated and recycled here.
! Gamma_inv(m+1) is calculated as 1/(m!).
!
! Gamma_inv(1+m+eps) is calculated from Gamma_inv(1-eps), 
! using the equalities:
! Gamma_inv(1-m-eps) = Gamma_inv(1-eps) \prod_{i=1}^{m} (1-eps-i), 
! where the product is 1 by definition if m = 0,
! Gamma_inv(1+m+eps) = (-1)^m sin (pi.eps)
! /[pi.(eps+m).Gamma_inv(1-m-eps)] 
! from Euler reflection formula, Gamma(x+1) = x.Gamma(x) equality, 
! and m+eps no zero.
! This scheme is much faster than 
! to recalculate Gamma_inv(1+m+eps) directly.
! 
! |eps|oo <= 0.1
! --------------
! The \beta_0/(1-z)^m expression is rewritten 
! so that it contains no instabilities:
! \beta_0/(1-z)^m = Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) 
! [(G(1,-eps) Gamma_inv(m+1) + G(m+1,eps))
! - Gamma_inv(1+m+eps) (G(a+m,eps) Gamma_inv(b+m+eps) 
! + G(b+m,eps) Gamma_inv(a+m)) 
! - E(log(1-z),eps) Gamma_inv(a+m) Gamma_inv(b+m) Gamma_inv(1+m+eps)] 
! (a)_m (b)_m Gamma(c)
!
! E(log(1-z),eps) is [(1-z)^eps - 1]/eps 
! if 1-m-eps and 1-m are different numerically, 
! and log(1-z) otherwise (eps ~ 0).
! If 1-m-eps and 1-m are equal numerically, 
! Gamma_inv(1+m+eps) is numerically equal to Gamma_inv(1+m), 
! already calculated as 1/(m!).
! See |eps|oo > 0.1 case for data recycling of other values 
! or for 1-m-eps and 1-m different numerically.
!
!----------------------------------------------------------------------
! Variables
! ---------
! a,b,c,one_minus_z: a,b,c and 1-z parameters and arguments 
! of the 2F1(a,b,c,z) function.
! m,eps: closest integer to c-a-b, with Re(c-a-b) >= 0 
! and eps = c-a-b-m
! Gamma_c,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm, Gamma_inv_eps_pb_pm: 
! recycled values of Gamma(c), Gamma_inv(1-eps), 
! Gamma_inv(a+m+eps) and Gamma_inv(b+m+eps).
! inf_norm_eps,phase,a_pm,b_pm,one_meps,Pi_eps,Pi_eps_pm: 
! |eps|oo,(-1)^m,a+m,b+m,1-eps,pi.eps,pi.(eps+m)
! Gamma_inv_one_meps_mm,Gamma_inv_eps_pm_p1: 
! Gamma_inv(1-m-eps) and Gamma_inv(1+m+eps) 
! calculated with the recycling scheme.
! prod1: (a)_m (b)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) 
! x Gamma_inv(b+m+eps) Gamma_inv(m+1) in |eps|oo > 0.1 case.
! prod2: (1-z)^eps Gamma_inv(a) Gamma_inv(b) Gamma_inv(1+m+eps) 
! in |eps|oo > 0.1 case.
! Gamma_inv_mp1,prod_ab: Gamma_inv(m+1) calculated as 1/(m!) 
! and (a)_m (b)_m in |eps|oo <= 0.1 case.
! is_eps_non_zero: true if 1-m-eps and 1-m are different numerically,
! false if not.
! Gamma_inv_a_pm,Gamma_inv_b_pm,z_term: Gamma_inv(a+m),Gamma_inv(b+m),
! E(eps,log(1-z))
! prod1: Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) 
! x [(G(1,-eps) Gamma_inv(m+1) + G(m+1,eps)) in |eps|oo <= 0.1 case.
! prod2: Gamma_inv(1+m+eps) (G(a+m,eps) Gamma_inv(b+m+eps) 
! + G(b+m,eps) Gamma_inv(a+m))
! prod3: E(eps,log(1-z)) Gamma_inv(a+m) Gamma_inv(b+m) 
! Gamma_inv(1+m+eps) 
! res: returned \beta_0/(1-z)^m value in all cases.
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_B_SUM_INIT_PS_ONE(A,B,GAMMA_C,GAMMA_INV_ONE_MEPS, &
     GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM,MZP1,M,EPS)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  INTEGER(IPR),INTENT(IN) :: M
  COMPLEX(PR),INTENT(IN) :: A,B,GAMMA_C,GAMMA_INV_ONE_MEPS, &
       GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM,MZP1,EPS
  INTEGER(IPR) :: M_M1,N,I,PHASE
  REAL(PR)     :: INF_NORM_EPS,GAMMA_INV_MP1
  COMPLEX(PR)  :: A_PM,CWFN_MIST_B_SUM_INIT_PS_ONE,PI_EPS,GAMMA_INV_ONE_MEPS_MM
  COMPLEX(PR)  :: B_PM,TMP1,TMP2
  COMPLEX(PR)  :: Z_TERM,PROD1,PROD2,PROD3,ONE_MEPS,PI_EPS_PM
  COMPLEX(PR)  :: GAMMA_INV_A_PM,PROD_AB,CWFN_MIST_GAMMA_INV,GAMMA_INV_B_PM
  COMPLEX(PR)  :: CWFN_MIST_GAMMA_INV_DIFF_EPS,GAMMA_INV_EPS_PM_P1
  !       
  INF_NORM_EPS=INF_NORM(EPS); M_M1=M-1; A_PM=A+M; B_PM=B+M
  ONE_MEPS=ONE-EPS; PI_EPS=M_PI*EPS; PI_EPS_PM = M_PI*(EPS+M)
  IF(MOD(M,2).EQ.0) THEN
     PHASE = 1
  ELSE
     PHASE = -1
  ENDIF
  GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS
  DO I=1,M
     GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS_MM*(ONE_MEPS - I)
  ENDDO
  IF(INF_NORM_EPS.GT.0.10D0) THEN
     GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
          /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
     PROD1=GAMMA_INV_ONE_MEPS*GAMMA_INV_EPS_PA_PM*GAMMA_INV_EPS_PB_PM
     DO N=0,M_M1
        PROD1=PROD1*(A+N)*(B+N)/(N+ONE)
     ENDDO
     PROD2=CWFN_MIST_GAMMA_INV(A)*CWFN_MIST_GAMMA_INV(B)*GAMMA_INV_EPS_PM_P1*(MZP1**EPS)
     CWFN_MIST_B_SUM_INIT_PS_ONE=GAMMA_C*(PROD1-PROD2)/EPS
     RETURN
  ELSE
     GAMMA_INV_MP1=ONE;PROD_AB=ONE
     DO N=0,M_M1
        GAMMA_INV_MP1 = GAMMA_INV_MP1/(N+ONE)
        PROD_AB = PROD_AB*(A+N)*(B+N)
     ENDDO
     IF(ONE_MEPS-M.NE.1-M) THEN
        Z_TERM=EXPM1(EPS*LOG(MZP1))/EPS
        GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
             /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
     ELSE
        Z_TERM=LOG(MZP1)
        GAMMA_INV_EPS_PM_P1 = GAMMA_INV_MP1
     ENDIF
     GAMMA_INV_A_PM=CWFN_MIST_GAMMA_INV(A_PM);GAMMA_INV_B_PM=CWFN_MIST_GAMMA_INV(B_PM)
     TMP1=ONE; TMP2=M+1;
     PROD1 = GAMMA_INV_EPS_PA_PM*GAMMA_INV_EPS_PB_PM    &
          *(GAMMA_INV_MP1*CWFN_MIST_GAMMA_INV_DIFF_EPS(TMP1,-EPS) &
          +CWFN_MIST_GAMMA_INV_DIFF_EPS(TMP2,EPS))
     PROD2 = GAMMA_INV_EPS_PM_P1 &
          *(GAMMA_INV_EPS_PB_PM*CWFN_MIST_GAMMA_INV_DIFF_EPS(A_PM,EPS) &
          +GAMMA_INV_A_PM*CWFN_MIST_GAMMA_INV_DIFF_EPS (B_PM,EPS))
     PROD3 = GAMMA_INV_A_PM*GAMMA_INV_B_PM*GAMMA_INV_EPS_PM_P1*Z_TERM
     CWFN_MIST_B_SUM_INIT_PS_ONE=GAMMA_C*PROD_AB*(PROD1-PROD2-PROD3)
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_B_SUM_INIT_PS_ONE
!
!----------------------------------------------------------------------
! Calculation of the first term of the B(z) power series 
! ------------------------------------------------------
! in the 1/z transformation, divided by z^{-m}
!---------------------------------------------
! In the 1/z transformation, the power series 
! B(z) = \sum_{n=0}^{+oo} \beta_n z^{-n} occurs 
! (see Comp. Phys. Comm. paper).
! The first term \beta_0, divided by z^{-m}, is calculated here. 
! m is the closest integer to Re(b-a) >= 0 and eps = b-a-m.
!
! One has to consider |eps|oo > 0.1 and |eps|oo <= 0.1, 
! where 1-m-eps and 1-m can be different or equal numerically, 
! leading to some changes in this last case.
!
! |eps|oo > 0.1
! -------------
! One has \beta_0/z^{-m} = [(a)_m (1-c+a)_m Gamma_inv(1-eps) 
! Gamma_inv(a+m+eps) Gamma_inv(c-a) Gamma_inv(m+1)
! - (-z)^{-eps} (1-c+a+eps)_m Gamma_inv(a) Gamma_inv(c-a-eps) 
! Gamma_inv(1+m+eps)].[Gamma(c)/eps], 
! stable in this regime for a direct evaluation.
!
! The values of Gamma(c), Gamma_inv(c-a) and Gamma_inv(a+m+eps) 
! were already calculated and recycled here.
! Gamma_inv(m+1) is calculated as 1/(m!). 
! Gamma_inv(1+m+eps) is calculated from Gamma_inv(1-eps) 
! as in the 1-z transformation routine.
! 
! |eps|oo <= 0.1
! --------------
! The \beta_0/z^{-m} expression is rewritten 
! so that it contains no instabilities:
! \beta_0/z^{-m} = [((1-c+a+eps)_m G(1,-eps) - P(m,eps,1-c+a) 
! Gamma_inv(1-eps)) Gamma_inv(c-a) Gamma_inv(a+m+eps) Gamma_inv(m+1)
! + (1-c+a+eps)_m [G(m+1,eps) Gamma_inv(c-a) Gamma_inv(a+m+eps) 
! - G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)]
! - (G(c-a,-eps) - E(log(-z),-eps)) Gamma_inv(m+1+eps) 
! Gamma_inv(a+m)]] (a)_m Gamma(c)
!
! Definitions and method are the same 
! as in the 1-z transformation routine, except for P(m,eps,1-c+a).
! P(m,eps,s) = [(s+eps)_m - (s)_m]/eps 
! for eps non zero and has a limit for eps -> 0.
! Let n0 be the closest integer to -Re(s) for s complex. 
! A stable formula available for eps -> 0 for P(m,eps,s) is:
! P(m,eps,s) = (s)_m E(\sum_{n=0}^{m-1} L(1/(s+n),eps),eps) 
! if n0 is not in [0:m-1],
! P(m,eps,s) = \prod_{n=0, n not equal to n0}^{m-1} (s+eps+n) 
! + (s)_m E(\sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps),eps) 
! if n0 is in [0:m-1].
! L(s,eps) is log1p(s eps)/eps if eps is not zero, 
! and L(s,0) = s.
! This expression is used in the code.
!
! Variables
! ---------
! a,b,c,z: a,b,c and z parameters 
! and arguments of the 2F1(a,b,c,z) function.
! m,eps: closest integer to b-a, with Re(b-a) >= 0 and eps = b-a-m.
! Gamma_c,Gamma_inv_cma,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm: 
! recycled values of Gamma(c), Gamma_inv(c-a), Gamma_inv(1-eps) 
! and Gamma_inv(a+m+eps).
! inf_norm_eps,phase,cma,a_mc_p1,a_mc_p1_pm,cma_eps,eps_pa_mc_p1,a_pm: 
! |eps|oo,(-1)^m,c-a,1-c+a+m,c-a-eps,1-c+a+eps,a+m
! Gamma_inv_cma_meps,one_meps,Pi_eps,Pi_eps_pm: 
! Gamma_inv(c-a-eps),1-eps,pi.eps,pi.(eps+m)
! Gamma_inv_one_meps_mm,Gamma_inv_eps_pm_p1: Gamma_inv(1-m-eps) 
! and Gamma_inv(1+m+eps) calculated with the recycling scheme.
! prod1: (a)_m (1-c+a)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) 
! x Gamma_inv(c-a) Gamma_inv(m+1) in |eps|oo > 0.1 case.
! prod2: (-z)^{-eps} (1-c+a+eps)_m Gamma_inv(a) 
! x Gamma_inv(c-a-eps) Gamma_inv(1+m+eps) in |eps|oo > 0.1 case.
! n0: closest integer to -Re(1-c+a)
! is_n0_here: true is n0 belongs to [0:m-1], false if not.
! is_eps_non_zero: true if 1-m-eps and 1-m are different numerically, 
! false if not.
! Gamma_inv_mp1,prod_a,prod_a_mc_p1: 
! Gamma_inv(m+1) calculated as 1/(m!), 
! (a)_m and (1-c+a)_m in |eps|oo <= 0.1 case.
! prod_eps_pa_mc_p1_n0: 
! \prod_{n=0, n not equal to n0}^{m-1} (1-c+a+eps+n) 
! if n0 belongs to [0:m-1], 0.0 if not, in |eps|oo <= 0.1 case.
! prod_eps_pa_mc_p1: (1-c+a+eps)_m in |eps|oo <= 0.1 case.
! sum: \sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps) if 1-m-eps 
! and 1-m are different numerically, 
! \sum_{n=0, n not equal to n0}^{m-1} 1/(s+n) if not.
! a_pn,a_mc_p1_pn,eps_pa_mc_p1_pn: a+n,1-c+a+n,1-c+a+eps+n values 
! used in (a)_m, (1-c+a)_m and (1-c+a+eps)_m evaluations.
! sum_term,prod_diff_eps,z_term: 
! E(\sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps),eps), 
! P(m,eps,1-c+a), -E(-eps,log(-z))
! Gamma_inv_a_pm,Gamma_prod1: Gamma_inv(a+m), 
! Gamma_inv(c-a).Gamma_inv(a+m+eps)
! prod1: ((1-c+a+eps)_m G(1,-eps) 
! - P(m,eps,1-c+a) Gamma_inv(1-eps)) Gamma_inv(c-a) 
! x Gamma_inv(a+m+eps) Gamma_inv(m+1)
! prod_2a: Gamma_inv(c-a).Gamma_inv(a+m+eps).G(m+1,eps)
! prod_2b: G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)
! prod_2c: (G(c-a,-eps) 
! - E(log(-z),-eps)) Gamma_inv(m+1+eps) Gamma_inv(a+m)
! prod2: (1-c+a+eps)_m [G(m+1,eps) Gamma_inv(c-a) Gamma_inv(a+m+eps) 
! - G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)] 
! - (G(c-a,-eps) - E(log(-z),-eps)) 
! x Gamma_inv(m+1+eps) Gamma_inv(a+m)]]
! res: returned \beta_0/z^{-m} value in all cases.
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_B_SUM_INIT_PS_INFINITY(A,C,GAMMA_C,GAMMA_INV_CMA, &
     GAMMA_INV_ONE_MEPS,GAMMA_INV_EPS_PA_PM,Z,M,EPS)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  INTEGER(IPR),INTENT(IN) :: M
  COMPLEX(PR),INTENT(IN) :: A,C,GAMMA_C,GAMMA_INV_CMA,Z,EPS
  COMPLEX(PR),INTENT(IN) :: GAMMA_INV_ONE_MEPS,GAMMA_INV_EPS_PA_PM
  INTEGER(IPR) :: M_M1,I,N,N0,PHASE
  LOGICAL      :: IS_N0_HERE,IS_EPS_NON_ZERO
  REAL(PR)     :: INF_NORM_EPS,NP1,GAMMA_INV_MP1
  COMPLEX(PR)  :: CWFN_MIST_B_SUM_INIT_PS_INFINITY,CWFN_MIST_GAMMA_INV,TMP1
  COMPLEX(PR)  :: CMA,A_MC_P1,A_MC_P1_PM,CMA_MEPS,EPS_PA_MC_P1,A_PM
  COMPLEX(PR)  :: GAMMA_INV_EPS_PM_P1,GAMMA_INV_CMA_MEPS,PI_EPS
  COMPLEX(PR)  :: PROD1,PROD2,A_PN,A_MC_P1_PN,ONE_MEPS
  COMPLEX(PR)  :: PROD_A,PROD_A_MC_P1,PROD_EPS_PA_MC_P1_N0,PI_EPS_PM
  COMPLEX(PR)  :: PROD_EPS_PA_MC_P1,SUM_N0,Z_TERM,SUM_TERM
  COMPLEX(PR)  :: PROD_DIFF_EPS,GAMMA_INV_A_PM,GAMMA_PROD1
  COMPLEX(PR)  :: PROD_2A,PROD_2B,PROD_2C,CWFN_MIST_GAMMA_INV_DIFF_EPS
  COMPLEX(PR)  :: EPS_PA_MC_P1_PN,GAMMA_INV_ONE_MEPS_MM
  !
  INF_NORM_EPS=INF_NORM(EPS); CMA=C-A; A_MC_P1=A-C+ONE
  A_MC_P1_PM=A_MC_P1+M; CMA_MEPS=CMA-EPS; EPS_PA_MC_P1=EPS+A_MC_P1
  A_PM=A+M; M_M1=M-1; ONE_MEPS=ONE-EPS; PI_EPS=M_PI*EPS
  PI_EPS_PM=M_PI*(EPS+M); GAMMA_INV_CMA_MEPS=CWFN_MIST_GAMMA_INV(CMA_MEPS)
  IF(MOD(M,2).EQ.0) THEN
     PHASE = 1
  ELSE
     PHASE = -1
  ENDIF
  GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS
  DO I=1,M
     GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS_MM*(ONE_MEPS - I)
  ENDDO
  IF(INF_NORM_EPS.GT.0.1D0) THEN
     GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
          /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
     PROD1 = GAMMA_INV_CMA*GAMMA_INV_EPS_PA_PM*GAMMA_INV_ONE_MEPS
     PROD2 = CWFN_MIST_GAMMA_INV(A)*GAMMA_INV_CMA_MEPS*GAMMA_INV_EPS_PM_P1 &
          *((-Z)**(-EPS))
     DO N=0,M_M1
        A_PN=A+N; A_MC_P1_PN=A_MC_P1+N
        EPS_PA_MC_P1_PN=EPS+A_MC_P1_PN;NP1=N+ONE
        PROD1 = PROD1*A_PN*A_MC_P1_PN/NP1
        PROD2 = PROD2*EPS_PA_MC_P1_PN
     ENDDO
     CWFN_MIST_B_SUM_INIT_PS_INFINITY = GAMMA_C*(PROD1-PROD2)/EPS
     RETURN
  ELSE
     N0=-NINT(REAL(A_MC_P1,PR))
     IS_EPS_NON_ZERO=ONE_MEPS-M.NE.1-M
     IS_N0_HERE=(N0.GE.0).AND.(N0.LT.M)     
     GAMMA_INV_MP1=ONE; PROD_A=ONE; PROD_A_MC_P1=ONE
     PROD_EPS_PA_MC_P1=ONE; SUM_N0=ZERO
     IF(IS_N0_HERE) THEN
        PROD_EPS_PA_MC_P1_N0 = ONE
     ELSE
        PROD_EPS_PA_MC_P1_N0 = ZERO
     ENDIF
     DO N=0,M_M1
        A_PN=A+N; A_MC_P1_PN=A_MC_P1+N
        EPS_PA_MC_P1_PN=EPS+A_MC_P1_PN; NP1=N+ONE
        PROD_A = PROD_A*A_PN
        PROD_A_MC_P1 = PROD_A_MC_P1*A_MC_P1_PN
        PROD_EPS_PA_MC_P1 = PROD_EPS_PA_MC_P1*EPS_PA_MC_P1_PN
        GAMMA_INV_MP1 = GAMMA_INV_MP1/NP1
        IF(N.NE.N0) THEN
           IF(IS_N0_HERE) THEN
              PROD_EPS_PA_MC_P1_N0=PROD_EPS_PA_MC_P1_N0 &
                   *EPS_PA_MC_P1_PN
           ENDIF
           IF(IS_EPS_NON_ZERO) THEN
              SUM_N0 = SUM_N0 + LOG1P(EPS/A_MC_P1_PN)
           ELSE
              SUM_N0 = SUM_N0 + ONE/A_MC_P1_PN
           ENDIF
        ENDIF
     ENDDO
     IF(IS_EPS_NON_ZERO) THEN
        GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
             /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
        SUM_TERM = EXPM1(SUM_N0)/EPS
        Z_TERM = EXPM1(-EPS*LOG(-Z))/EPS
     ELSE
        GAMMA_INV_EPS_PM_P1 = GAMMA_INV_MP1
        SUM_TERM = SUM_N0
        Z_TERM = -LOG(-Z)
     ENDIF
     PROD_DIFF_EPS = PROD_EPS_PA_MC_P1_N0 + PROD_A_MC_P1*SUM_TERM
     GAMMA_INV_A_PM = CWFN_MIST_GAMMA_INV(A_PM)
     GAMMA_PROD1=GAMMA_INV_CMA*GAMMA_INV_EPS_PA_PM
     TMP1=ONE
     PROD1 = GAMMA_PROD1*GAMMA_INV_MP1*(CWFN_MIST_GAMMA_INV_DIFF_EPS(TMP1,-EPS) &
          *PROD_EPS_PA_MC_P1 - GAMMA_INV_ONE_MEPS*PROD_DIFF_EPS)
     TMP1=M+1
     PROD_2A = GAMMA_PROD1*CWFN_MIST_GAMMA_INV_DIFF_EPS(TMP1,EPS) 
     PROD_2B = GAMMA_INV_CMA*GAMMA_INV_EPS_PM_P1  &
          *CWFN_MIST_GAMMA_INV_DIFF_EPS(A_PM,EPS)
     PROD_2C = GAMMA_INV_EPS_PM_P1*GAMMA_INV_A_PM &
          *(CWFN_MIST_GAMMA_INV_DIFF_EPS(CMA,-EPS) + GAMMA_INV_CMA_MEPS*Z_TERM)
     PROD2 = PROD_EPS_PA_MC_P1*(PROD_2A - PROD_2B - PROD_2C)
     CWFN_MIST_B_SUM_INIT_PS_INFINITY = GAMMA_C*PROD_A*(PROD1+PROD2)
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_B_SUM_INIT_PS_INFINITY
!
!----------------------------------------------------------------------
! Calculation of the derivative of the polynomial P(X) 
! ----------------------------------------------------
! testing power series convergence
! --------------------------------
! P(X) = |z(a+X)(b+X)|^2 - |(c+X)(X+1)|^2 
!      = \sum_{i=0}^{4} c[i] X^{i}, for |z| < 1.
! It is positive when the power series term modulus increases 
! and negative when it decreases, 
! so that its derivative provides information on its convergence 
! (see Comp. Phys. Comm. paper).
! Its derivative components cv_poly_der_tab[i] = (i+1) c[i+1] 
! for i in [0:3] 
! so that P'(X) = \sum_{i=0}^{3} cv_poly_der_tab[i] X^{i} 
! are calculated.
!
! Variables:
! ----------
! a,b,c,z: a,b,c and z parameters and arguments 
! of the 2F1(a,b,c,z) function.
! cv_poly_der_tab[3]: table of four doubles 
! containing the P'(X) components.
! mod_a2,mod_b2,mod_c2,mod_z2,R_a,Re_b,Re_c: |a|^2, |b|^2, |c|^2, 
! |z|^2, Re(a), Re(b), Re(c), with which P(X) can be expressed.
!----------------------------------------------------------------------
SUBROUTINE CWFN_MIST_CV_POLY_DER_TAB_CALC(A,B,C,Z,CV_POLY_DER_TAB)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  REAL(PR),INTENT(OUT) :: CV_POLY_DER_TAB(0:3)
  REAL(PR)     :: MOD_A2,MOD_B2,MOD_C2,MOD_Z2
  REAL(PR)     :: RE_A,RE_B,RE_C,IM_A,IM_B,IM_C,RE_Z,IM_Z
  !
  RE_A=REAL(A,PR); IM_A=AIMAG(A); MOD_A2=RE_A*RE_A+IM_A*IM_A
  RE_B=REAL(B,PR); IM_B=AIMAG(B); MOD_B2=RE_B*RE_B+IM_B*IM_B
  RE_C=REAL(C,PR); IM_C=AIMAG(C); MOD_C2=RE_C*RE_C+IM_C*IM_C
  RE_Z=REAL(Z,PR); IM_Z=AIMAG(Z); MOD_Z2=RE_Z*RE_Z+IM_Z*IM_Z
  CV_POLY_DER_TAB(0)=TWO*((RE_A*MOD_B2+RE_B*MOD_A2)*MOD_Z2-RE_C-MOD_C2)
  CV_POLY_DER_TAB(1)=TWO*((MOD_A2+MOD_B2+4.0D0*RE_A*RE_B)*MOD_Z2 &
       -ONE-4.0D0*RE_C-MOD_C2)
  CV_POLY_DER_TAB(2)=6.0D0*((RE_A+RE_B)*MOD_Z2-RE_C-ONE)
  CV_POLY_DER_TAB(3)=4.0D0*(MOD_Z2-ONE)
END SUBROUTINE CWFN_MIST_CV_POLY_DER_TAB_CALC
!
!----------------------------------------------------------------------
! Calculation of the derivative of the polynomial P(X) 
! ----------------------------------------------------
! testing power series convergence at one x value
! -----------------------------------------------
! P'(x) is calculated for a real x. 
! See P'(X) components calculation routine for definitions.
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_CV_POLY_DER_CALC(CV_POLY_DER_TAB,X)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  REAL(PR),INTENT(IN) :: X
  REAL(PR),INTENT(IN) :: CV_POLY_DER_TAB(0:3)
  REAL(PR) :: CWFN_MIST_CV_POLY_DER_CALC
  !
  CWFN_MIST_CV_POLY_DER_CALC=CV_POLY_DER_TAB(0)+X*(CV_POLY_DER_TAB(1) &
       +X*(CV_POLY_DER_TAB(2)+X*CV_POLY_DER_TAB(3)))
  RETURN
END FUNCTION CWFN_MIST_CV_POLY_DER_CALC
!
!----------------------------------------------------------------------
! Calculation of an integer after which false convergence cannot occur
! --------------------------------------------------------------------
! See cv_poly_der_tab_calc routine for definitions.
! If P'(x) < 0 and P''(x) < 0 for x > xc, it will be so for all x > xc 
! as P(x) -> -oo for x -> +oo 
! and P(x) can have at most one maximum for x > xc. 
! It means that the 2F1 power series term modulus will increase 
! or decrease to 0 for n > nc, 
! with nc the smallest positive integer larger than xc.
!
! If P'(X) = C0 + C1.X + C2.X^2 + C3.X^3, 
! the discriminant of P''(X) is Delta = C2^2 - 3 C1 C3.
!
! If Delta > 0, P''(X) has two different real roots 
! and its largest root is -(C2 + sqrt(Delta))/(3 C3), 
! because C3 = 4(|z|^2 - 1) < 0.
! One can take xc = -(C2 + sqrt(Delta))/(3 C3) 
! and one returns its associated nc integer.
!
! If Delta <= 0, P''(X) has at most one real root, 
! so that P'(X) has only one root and then P(X) only one maximum.
! In this case, one can choose xc = nc = 0, which is returned.
!
! Variables
! ---------
! cv_poly_der_tab: table of four doubles 
! containing the P'(X) coefficients
! C1,C2,three_C3: cv_poly_der_tab[1], cv_poly_der_tab[2] 
! and 3.0*cv_poly_der_tab[3], so that P''(X) = C1 + 2.C2.x + three_C3.x^2
! Delta: discriminant of P''(X), equal to C2^2 - 3 C1 C3.
! largest_root: if Delta > 0, 
! P''(X) largest real root equal to -(C2 + sqrt(Delta))/(3 C3).
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_MIN_N_CALC(CV_POLY_DER_TAB)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  REAL(PR),INTENT(IN) :: CV_POLY_DER_TAB(0:3)
  INTEGER(IPR) :: CWFN_MIST_MIN_N_CALC
  REAL(PR)     :: C1,C2,THREE_C3,DELTA,LARGEST_ROOT
  !
  C1=CV_POLY_DER_TAB(1); C2=CV_POLY_DER_TAB(2)
  THREE_C3=3.0D0*CV_POLY_DER_TAB(3); DELTA = C2*C2 - THREE_C3*C1
  IF(DELTA.LE.ZERO) THEN
     CWFN_MIST_MIN_N_CALC = 0
     RETURN
  ELSE
     LARGEST_ROOT = -(C2 + SQRT (DELTA))/THREE_C3
     CWFN_MIST_MIN_N_CALC = MAX(CEILING(LARGEST_ROOT),0)
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_MIN_N_CALC
!
!----------------------------------------------------------------------
! Calculation of the 2F1 power series converging for |z| < 1
! ----------------------------------------------------------
! One has 2F1(a,b,c,z) 
! = \sum_{n = 0}^{+oo} (a)_n (b)_n / ((c)_n n!) z^n,
! so that 2F1(a,b,c,z) = \sum_{n = 0}^{+oo} t[n] z^n, 
! with t[0] = 1 and t[n+1] = (a+n)(b+n)/((c+n)(n+1)) t[n] for n >= 0.
! If a or b are negative integers, 
! F(z) is a polynomial of degree -a or -b, evaluated directly.
! If not, one uses the test of convergence |t[n] z^n|oo < 1E-15 
! to truncate the series after it was checked 
! that false convergence cannot occur.
! Variables:
! ----------
! a,b,c,z: a,b,c and z parameters and arguments 
! of the 2F1(a,b,c,z) function. One must have here |z| < 1.
! term,sum: term of the 2F1 power series equal to t[n] z^n, 
! truncated sum at given n of the 2F1 power series.
! na,nb: absolute values of the closest integers to Re(a) and Re(b). 
! a = -na or b = -nb means one is in the polynomial case.
! cv_poly_der_tab: coefficients of the derivative 
! of the polynomial P(X) = |z(a+X)(b+X)|^2 - |(c+X)(X+1)|^2
! min_n: smallest integer after which false convergence cannot occur. 
! It is calculated in min_n_calc.
! possible_false_cv: always true if n < min_n. 
! If n >= min_n, it is true if P'(n) > 0. 
! If n >= min_n and P'(n) < 0, 
! it becomes false and remains as such for the rest of the calculation. 
! One can then check if |t[n] z^n|oo < 1E-15 to truncate the series.
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_HYP_PS_ZERO(A,B,C,Z)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  INTEGER(IPR) :: N,NA,NB,MIN_N,CWFN_MIST_MIN_N_CALC
  COMPLEX(PR)  :: CWFN_MIST_HYP_PS_ZERO,TERM
  LOGICAL :: POSSIBLE_FALSE_CV
  REAL(PR) :: CV_POLY_DER_TAB(0:3)
  REAL(PR) :: CWFN_MIST_CV_POLY_DER_CALC
  !
  NA = ABS(NINT(REAL(A,PR)))
  NB = ABS(NINT(REAL(B,PR)))
  TERM=ONE; CWFN_MIST_HYP_PS_ZERO=ONE  
  IF(A.EQ.(-NA)) THEN
     DO N=0,NA-1
        TERM = TERM*Z*(A+N)*(B+N)/((N+ONE)*(C+N))
        CWFN_MIST_HYP_PS_ZERO = CWFN_MIST_HYP_PS_ZERO + TERM
     ENDDO
     RETURN
  ELSE IF(B.EQ.(-NB)) THEN
     DO N=0,NB-1
        TERM = TERM*Z*(A+N)*(B+N)/((N+ONE)*(C+N))
        CWFN_MIST_HYP_PS_ZERO = CWFN_MIST_HYP_PS_ZERO + TERM
     ENDDO
     RETURN
  ELSE
     CALL CWFN_MIST_CV_POLY_DER_TAB_CALC(A,B,C,Z,CV_POLY_DER_TAB)
     POSSIBLE_FALSE_CV=.TRUE.
     MIN_N=CWFN_MIST_MIN_N_CALC(CV_POLY_DER_TAB);N=0
     DO WHILE(POSSIBLE_FALSE_CV.OR.(INF_NORM(TERM).GT.EPS15))
        TERM = TERM*Z*(A+N)*(B+N)/((N+ONE)*(C+N))
        CWFN_MIST_HYP_PS_ZERO = CWFN_MIST_HYP_PS_ZERO + TERM
        IF(POSSIBLE_FALSE_CV.AND.(N.GT.MIN_N)) THEN
           POSSIBLE_FALSE_CV = &
                (CWFN_MIST_CV_POLY_DER_CALC (CV_POLY_DER_TAB,DBLE(N)).GT.ZERO)
        ENDIF
        N=N+1 
     ENDDO
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_HYP_PS_ZERO
!
!----------------------------------------------------------------------
! Calculation of the 2F1 power series 
! -----------------------------------
! converging with the 1-z transformation
! --------------------------------------
! The formula for F(z) in the 1-z transformation holds:
! F(z) = (-1)^m (pi.eps)/sin (pi.eps) [A(z) + B(z)] 
! for eps not equal to zero, F(z) = (-1)^m [A(z) + B(z)] for eps = 0
! where m = |Re(c-a-b)], eps = c-a-b-m, 
! A(z) = \sum_{n=0}^{m-1} alpha[n] (1-z)^n, 
! B(z) = \sum_{n=0}^{+oo} beta[n] (1-z)^n, and:
!
! alpha[0] = [Gamma_inv(1-m-eps)/eps] Gamma_inv(a+m+eps) 
!          x Gamma_inv(b+m+eps) Gamma(c)
! [Gamma_inv(1-m-eps)/eps] is calculated in A_sum_init. 
! alpha[0] is calculated with log[Gamma] 
! if the previous expression might overflow, 
! and its imaginary part removed if a, b and c are real.
! alpha[n+1] = (a+n)(b+n)/[(n+1)(1-m-eps+n)] alpha[n], n in [0:m-2].
!
! beta[0] is defined in B_sum_init_PS_one function comments.
! gamma[0] = Gamma(c) (a)_m (b)_m (1-z)^m Gamma_inv(a+m+eps) 
!          x Gamma_inv(b+m+eps) Gamma_inv(m+1) Gamma_inv(1-eps)
!
! beta[n+1] = (a+m+n+eps)(b+m+n+eps)/[(m+n+1+eps)(n+1)] beta[n]
! + [(a+m+n)(b+m+n)/(m+n+1) - (a+m+n) - (b+m+n) - eps 
! + (a+m+n+eps)(b+m+n+eps)/(n+1)]
!             x gamma[n]/[(n+m+1+eps)(n+1+eps)], n >= 0.
! gamma[n+1] = (a+m+n)(b+m+n)/[(m+n+1)(n+1-eps)] gamma[n], n >= 0.
!
! B(z) converges <=> |1-z| < 1
! The test of convergence is |beta[n] (1-z)^n|oo < 1E-15 |beta[0]|oo
! for n large enough so that false convergence cannot occur.
!
! Variables
! ---------
! a,b,c,one_minus_z: a,b,c parameters 
! and 1-z from z argument of 2F1(a,b,c,z)
! m,phase,m_p1,eps,eps_pm,eps_pm_p1,
! a_pm,b_pm,one_meps,one_meps_pm: 
! |Re(c-a-b)], (-1)^m, m+1, c-a-b-m, 
! eps+m, eps+m+1, a+m, b+m, 1-eps, 1-eps-m
! eps_pa,eps_pb,eps_pa_pm,eps_pb_pm,Pi_eps,Gamma_c: 
! eps+a, eps+b, eps+a+m, eps+b+m, pi.eps, Gamma(c)
! Gamma_inv_eps_pa_pm,Gamma_inv_eps_pb_pm,Gamma_prod: 
! Gamma_inv(eps+a+m), Gamma_inv(eps+b+m), 
! Gamma(c).Gamma_inv(eps+a+m).Gamma_inv(eps+b+m)
! Gamma_inv_one_meps,A_first_term,A_sum,A_term: 
! Gamma_inv(1-eps), alpha[0], A(z), alpha[n] (1-z)^n
! pow_mzp1_m,B_first_term,prod_B,ratio: (1-z)^m, beta[0], 
! (a)_m (b)_m (1-z)^m, (a+n)(b+n)/(n+1) for n in [0:m-2].
! B_extra_term,B_term,B_sum,B_prec: 
! gamma[n], beta[n] (1-z)^n, B(z), 1E-15 |beta[0|oo
! cv_poly1_der_tab,cv_poly2_der_tab: P1'(X) and P2'(X) coefficients 
! of the potentials derivatives of P1(X) and P2(X) 
! defined in cv_poly_der_tab_calc with parameters 
! a1 = a, b1 = b, c1 = 1-m-eps, z1 = 1-z 
! and a2 = eps+b+m, b2 = eps+a+m,c2 = eps+m+1, z2 = 1-z.
! min_n: smallest integer after which false convergence cannot occur. 
! It is calculated in min_n_calc with both P1'(X) and P2'(X), 
! so one takes the largest integer coming from both calculations.
! possible_false_cv: always true if n < min_n. 
! If n >= min_n, it is true if P1'(n) > 0 or P2'(n) > 0. 
! If n >= min_n and P1'(n) < 0 and P2'(n) < 0, 
! it becomes false and remains as such for the rest of the calculation.
! One can then check if |beta[n] z^n|oo < 1E-15 to truncate the series.
! n,n_pm_p1,n_p1,a_pm_pn,b_pm_pn,eps_pm_p1_pn,n_p1_meps,eps_pa_pm_pn,
! eps_pb_pm_pn,eps_pm_pn: index of power series, n+m+1, n+1, 
! a+m+n, b+m+n, eps+m+n+1, n+1-eps, eps+a+m+n, eps+b+m+n, eps+m+n,
! prod1,prod2,prod3: (eps+a+m+n)(eps+b+m+n), 
! (eps+m+1+n)(n+1), (a+m+n)(b+m+n)
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_HYP_PS_ONE(A,B,C,MZP1)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,MZP1
  INTEGER(IPR) :: N,M,PHASE,M_M2,MIN_N,CWFN_MIST_MIN_N_CALC,M_P1
  REAL(PR)     :: B_PREC,N_P1,N_PM_P1,CWFN_MIST_CV_POLY_DER_CALC 
  COMPLEX(PR)  :: CWFN_MIST_HYP_PS_ONE,EPS,EPS_PM,EPS_PM_P1,A_PM
  COMPLEX(PR)  :: B_PM,ONE_MEPS_MM,EPS_PA,EPS_PB,PI_EPS,GAMMA_PROD
  COMPLEX(PR)  :: EPS_PA_PM,EPS_PB_PM,CWFN_MIST_GAMMA_INV,CWFN_MIST_B_SUM_INIT_PS_ONE
  COMPLEX(PR)  :: CWFN_MIST_A_SUM_INIT,CWFN_MIST_LOG_A_SUM_INIT,A_SUM,A_TERM,ONE_MEPS
  COMPLEX(PR)  :: B_EXTRA_TERM,B_TERM,B_SUM,GAMMA_C,CWFN_MIST_LOG_GAMMA,RATIO
  COMPLEX(PR)  :: A_PM_PN,B_PM_PN,EPS_PM_P1_PN,N_P1_MEPS
  COMPLEX(PR)  :: PROD1,PROD2,PROD3
  COMPLEX(PR)  :: EPS_PA_PM_PN,EPS_PB_PM_PN,EPS_PM_PN,PROD_B,POW_MZP1_M
  COMPLEX(PR)  :: GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM
  COMPLEX(PR)  :: GAMMA_INV_ONE_MEPS
  LOGICAL :: POSSIBLE_FALSE_CV
  REAL(PR) :: CV_POLY1_DER_TAB(0:3),CV_POLY2_DER_TAB(0:3)
  !
  M=NINT(REAL(C-A-B,PR)); M_M2=M-2; M_P1=M+1
  IF(MOD(M,2).EQ.0) THEN
     PHASE=1
  ELSE
     PHASE=-1
  ENDIF
  EPS=C-A-B-M; EPS_PM=EPS+M; EPS_PM_P1=EPS_PM+ONE; A_PM=A+M;B_PM=B+M
  ONE_MEPS=ONE-EPS; ONE_MEPS_MM=ONE_MEPS-M; EPS_PA=EPS+A; EPS_PB=EPS+B 
  PI_EPS=M_PI*EPS; EPS_PA_PM=EPS_PA+M; EPS_PB_PM=EPS_PB+M
  GAMMA_C=ONE/CWFN_MIST_GAMMA_INV(C)
  GAMMA_INV_EPS_PA_PM=CWFN_MIST_GAMMA_INV(EPS_PA_PM)
  GAMMA_INV_EPS_PB_PM=CWFN_MIST_GAMMA_INV(EPS_PB_PM)
  GAMMA_PROD=GAMMA_C*GAMMA_INV_EPS_PA_PM*GAMMA_INV_EPS_PB_PM
  GAMMA_INV_ONE_MEPS=CWFN_MIST_GAMMA_INV(ONE_MEPS)
  IF(M.EQ.0) THEN
     A_TERM=ZERO
  ELSE IF(INF_NORM(ONE_MEPS_MM &
       *(LOG(ONE + ABS(ONE_MEPS_MM))-ONE)).LT.300.0d0) THEN
     A_TERM=GAMMA_PROD*CWFN_MIST_A_SUM_INIT(M,EPS,GAMMA_INV_ONE_MEPS)
  ELSE
     A_TERM=EXP(CWFN_MIST_LOG_GAMMA(C)-CWFN_MIST_LOG_GAMMA(EPS_PA_PM)&
          -CWFN_MIST_LOG_GAMMA(EPS_PB_PM)+CWFN_MIST_LOG_A_SUM_INIT(M,EPS))
     IF((AIMAG(A).EQ.ZERO).AND.(AIMAG(B).EQ.ZERO)&
          .AND.(AIMAG(C).EQ.ZERO)) THEN
        A_TERM=REAL(A_TERM,PR)
     ENDIF
  ENDIF
  A_SUM=A_TERM
  POW_MZP1_M = MZP1**M
  B_TERM=CWFN_MIST_B_SUM_INIT_PS_ONE(A,B,GAMMA_C,GAMMA_INV_ONE_MEPS, &
       GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM,MZP1,M,EPS)*POW_MZP1_M
  PROD_B=POW_MZP1_M
  DO N=0,M_M2
     RATIO=(A+N)*(B+N)/(N+ONE)
     A_TERM=A_TERM*MZP1*RATIO/(N+ONE_MEPS_MM)
     A_SUM=A_SUM+A_TERM
     PROD_B = PROD_B*RATIO
  ENDDO
  IF(M.GT.0) THEN
     PROD_B = PROD_B*(A+M-ONE)*(B+M-ONE)/DBLE(M)
  ENDIF
  B_EXTRA_TERM = PROD_B*GAMMA_PROD*GAMMA_INV_ONE_MEPS; B_SUM=B_TERM
  B_PREC=EPS15*INF_NORM(B_TERM)
  CALL CWFN_MIST_CV_POLY_DER_TAB_CALC(A,B,ONE_MEPS_MM,MZP1,CV_POLY1_DER_TAB)
  CALL CWFN_MIST_CV_POLY_DER_TAB_CALC(EPS_PB_PM,EPS_PA_PM,EPS_PM_P1,MZP1, &
       CV_POLY2_DER_TAB)
  MIN_N=MAX(CWFN_MIST_MIN_N_CALC(CV_POLY1_DER_TAB),CWFN_MIST_MIN_N_CALC(CV_POLY2_DER_TAB))
  POSSIBLE_FALSE_CV=.TRUE.; N=0
  DO WHILE(POSSIBLE_FALSE_CV.OR.(INF_NORM(B_TERM).GT.B_PREC))
     N_PM_P1=N+M_P1; N_P1=N+ONE; A_PM_PN=A_PM+N; B_PM_PN=B_PM+N
     EPS_PM_P1_PN=EPS_PM_P1+N; N_P1_MEPS=ONE_MEPS+N
     EPS_PM_PN=EPS_PM+N; EPS_PA_PM_PN=EPS_PA_PM+N 
     EPS_PB_PM_PN=EPS_PB_PM+N
     PROD1=EPS_PA_PM_PN*EPS_PB_PM_PN
     PROD2=EPS_PM_P1_PN*N_P1
     PROD3=A_PM_PN*B_PM_PN
     B_TERM = MZP1*(B_TERM*PROD1/PROD2+B_EXTRA_TERM*(PROD3/N_PM_P1 &
          -A_PM_PN-B_PM_PN-EPS+PROD1/N_P1)/(EPS_PM_P1_PN*N_P1_MEPS))
     B_SUM=B_SUM+B_TERM
     B_EXTRA_TERM=B_EXTRA_TERM*MZP1*PROD3/(N_PM_P1*N_P1_MEPS)
     IF(POSSIBLE_FALSE_CV.AND.(N.GT.MIN_N)) THEN
        POSSIBLE_FALSE_CV = &
             (CWFN_MIST_CV_POLY_DER_CALC(CV_POLY1_DER_TAB,DBLE(N)).GT.ZERO).OR. &
             (CWFN_MIST_CV_POLY_DER_CALC(CV_POLY2_DER_TAB,DBLE(N)).GT.ZERO)
     ENDIF
     N=N+1
  ENDDO
  IF(EPS.EQ.ZERO) THEN
     CWFN_MIST_HYP_PS_ONE=PHASE*(A_SUM+B_SUM)
     RETURN
  ELSE
     CWFN_MIST_HYP_PS_ONE=PHASE*(A_SUM+B_SUM)*PI_EPS/SIN(PI_EPS)
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_HYP_PS_ONE
!
!----------------------------------------------------------------------
! Calculation of the 2F1 power series 
! -----------------------------------
! converging with the 1/z transformation
! --------------------------------------
! The formula for F(z) in the 1/z transformation holds:
! F(z) = (-1)^m (pi.eps)/sin (pi.eps) [A(z) + B(z)] 
! for eps not equal to zero, 
! F(z) = (-1)^m [A(z) + B(z)] for eps = 0
! where m = |Re(b-a)], eps = b-a-m, 
! A(z) = \sum_{n=0}^{m-1} alpha[n] z^{-n}, 
! B(z) = \sum_{n=0}^{+oo} beta[n] z^{-n}, and:
!
! alpha[0] = [Gamma_inv(1-m-eps)/eps] Gamma_inv(c-a) 
!          x Gamma_inv(a+m+eps) Gamma(c)
! [Gamma_inv(1-m-eps)/eps] is calculated in A_sum_init. 
! alpha[0] is calculated with log[Gamma] 
! if the previous expression might overflow, 
! and its imaginary part removed if a, b and c are real.
! alpha[n+1] = (a+n)(1-c+a+n)/[(n+1)(1-m-eps+n)] alpha[n], 
! n in [0:m-2].
!
! beta[0] is defined in B_sum_init_PS_infinity function comments.
! gamma[0] = Gamma(c) (a)_m (1-c+a)_m z^{-m} Gamma_inv(a+m+eps) 
!          x Gamma_inv(c-a) Gamma_inv(m+1) Gamma_inv(1-eps)
!
! beta[n+1] = (a+m+n+eps)(1-c+a+m+n+eps)/[(m+n+1+eps)(n+1)] beta[n] 
! + [(a+m+n)(1-c+a+m+n)/(m+n+1) - (a+m+n) - (1-c+a+m+n) 
! - eps + (a+m+n+eps)(1-c+a+m+n+eps)/(n+1)]
! x gamma[n]/[(n+m+1+eps)(n+1+eps)], n >= 0.
! gamma[n+1] = (a+m+n)(b+m+n)/[(m+n+1)(n+1-eps)] gamma[n], n >= 0.
!
! B(z) converges <=> |z| > 1
! The test of convergence is |beta[n] z^{-n}|oo < 1E-15 |beta[0]|oo
! for n large enough so that false convergence cannot occur.
!
! Variables
! ---------
! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
! m,phase,m_p1,eps,a_mc_p1,one_meps,
! one_meps_pm,a_pm,a_mc_p1_pm,cma: |Re(b-a)], (-1)^m, m+1, b-a-m, 
! 1-c+a, 1-eps, 1-eps-m, a+m, 1-c+a+m, c-a
! eps_pa,eps_pm_p1,eps_pa_mc_p1_pm,Pi_eps,eps_pa_pm,eps_pm,Gamma_c: 
! eps+a, eps+m+1, eps+1-c+a+m, pi.eps, eps+a+m, eps+m, Gamma(c)
! Gamma_inv_eps_pa_pm,Gamma_inv_cma,z_inv,pow_mz_ma,
! Gamma_inv_one_meps,Gamma_prod: Gamma_inv(eps+a+m), Gamma_inv(c-a), 
! 1/z, (-z)^(-a), Gamma_inv(1-eps), 
! Gamma(c) Gamma_inv(c-a) Gamma_inv(eps+a+m)
! A_first_term,A_sum,A_term: alpha[0], A(z), alpha[n] z^{-n}
! pow_z_inv_m,B_first_term,prod_B,ratio: z^{-m}, beta[0], 
! (a)_m (1-c+a)_m z^{-m}, (a+n)(1-c+a+n)/(n+1) for n in [0:m-2].
! B_extra_term,B_term,B_sum,B_prec: 
! gamma[n], beta[n] z^{-n}, B(z), 1E-15 |beta[0|oo
! cv_poly1_der_tab,cv_poly2_der_tab: P1'(X) and P2'(X) coefficients 
! of the potentials derivatives of P1(X) and P2(X) 
! defined in cv_poly_der_tab_calc 
! with parameters a1 = a, b1 = 1-c+a, c1 = 1-m-eps, z1 = 1/z 
! and a2 = b, b2 = eps+1-c+a+m,c2 = eps+m+1, z2 = 1/z.
! min_n: smallest integer after which false convergence cannot occur. 
!        It is calculated in min_n_calc with both P1'(X) and P2'(X), 
! so one takes the largest integer coming from both calculations.
! possible_false_cv: always true if n < min_n. If n >= min_n, 
! it is true if P1'(n) > 0 or P2'(n) > 0. 
! If n >= min_n and P1'(n) < 0 and P2'(n) < 0, 
! it becomes false and remains as such for the rest of the calculation. 
! One can then check if |beta[n] z^n|oo < 1E-15 to truncate the series.
! n,n_pm_p1,n_p1,a_pm_pn,a_mc_p1_pm_pn,eps_pm_p1_pn,n_p1_meps,
! eps_pa_pm_pn,eps_pa_mc_p1_pm_pn,eps_pm_pn: 
! index of power series, n+m+1, n+1, a+m+n, 1-c+a+m+n, eps+m+n+1,
! n+1-eps, eps+a+m+n, eps+1-c+a+m+n, eps+m+n,
! prod1,prod2,prod3: (eps+a+m+n)(eps+1-c+a+m+n),
! (eps+m+1+n)(n+1), (a+m+n)(1-c+a+m+n)
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_HYP_PS_INFINITY(A,B,C,Z)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  INTEGER(IPR) :: N,M,PHASE,M_M2,MIN_N,CWFN_MIST_MIN_N_CALC,M_P1
  REAL(PR)     :: B_PREC,N_P1,N_PM_P1,CWFN_MIST_CV_POLY_DER_CALC
  COMPLEX(PR)  :: CWFN_MIST_B_SUM_INIT_PS_INFINITY,CWFN_MIST_LOG_GAMMA,POW_Z_INV_M
  COMPLEX(PR)  :: CWFN_MIST_HYP_PS_INFINITY,Z_INV,CWFN_MIST_GAMMA_INV,RATIO
  COMPLEX(PR)  :: EPS,A_MC_P1,ONE_MEPS,ONE_MEPS_MM,A_PM,A_MC_P1_PM
  COMPLEX(PR)  :: CMA,EPS_PA,EPS_PM_P1,EPS_PA_MC_P1_PM,PI_EPS
  COMPLEX(PR)  :: EPS_PA_PM,EPS_PM,GAMMA_C,GAMMA_INV_CMA,POW_MZ_MA
  COMPLEX(PR)  :: CWFN_MIST_A_SUM_INIT,CWFN_MIST_LOG_A_SUM_INIT,A_SUM,A_TERM
  COMPLEX(PR)  :: GAMMA_INV_EPS_PA_PM,GAMMA_INV_ONE_MEPS
  COMPLEX(PR)  :: PROD_B,B_EXTRA_TERM,B_TERM,B_SUM,PROD1
  COMPLEX(PR)  :: A_PM_PN,A_MC_P1_PM_PN,EPS_PM_P1_PN,N_P1_MEPS
  COMPLEX(PR)  ::   PROD2,PROD3,GAMMA_PROD
  COMPLEX(PR)  :: EPS_PA_PM_PN,EPS_PA_MC_P1_PM_PN,EPS_PM_PN
  LOGICAL :: POSSIBLE_FALSE_CV
  REAL(PR) :: CV_POLY1_DER_TAB(0:3),CV_POLY2_DER_TAB(0:3)
  !
  M=NINT(REAL(B-A,PR)); M_M2=M-2;M_P1=M+1
  IF(MOD(M,2).EQ.0) THEN
     PHASE=1
  ELSE
     PHASE=-1
  ENDIF
  EPS=B-A-M; A_MC_P1=ONE-C+A; ONE_MEPS=ONE-EPS; ONE_MEPS_MM=ONE_MEPS-M
  A_PM=A+M; A_MC_P1_PM=A_MC_P1+M; CMA=C-A; EPS_PA=EPS+A
  EPS_PM=EPS+M; EPS_PM_P1=EPS_PM+ONE; EPS_PA_MC_P1_PM=EPS+A_MC_P1_PM
  PI_EPS=M_PI*EPS; EPS_PA_PM=EPS_PA+M
  GAMMA_C=ONE/CWFN_MIST_GAMMA_INV(C); GAMMA_INV_EPS_PA_PM = CWFN_MIST_GAMMA_INV(EPS_PA_PM)
  GAMMA_INV_ONE_MEPS = CWFN_MIST_GAMMA_INV(ONE_MEPS)
  GAMMA_INV_CMA=CWFN_MIST_GAMMA_INV(CMA); Z_INV=ONE/Z;POW_MZ_MA=(-Z)**(-A)
  GAMMA_PROD=GAMMA_C*GAMMA_INV_CMA*GAMMA_INV_EPS_PA_PM
  IF(M.EQ.0) THEN
     A_TERM=ZERO
  ELSE IF(INF_NORM(ONE_MEPS_MM &
       *(LOG(ONE + ABS(ONE_MEPS_MM))-ONE)).LT.300.0d0) THEN
     A_TERM=GAMMA_PROD*CWFN_MIST_A_SUM_INIT(M,EPS,GAMMA_INV_ONE_MEPS)
  ELSE
     A_TERM=EXP(CWFN_MIST_LOG_GAMMA(C)-CWFN_MIST_LOG_GAMMA(CMA) &
          -CWFN_MIST_LOG_GAMMA(B) &
          + CWFN_MIST_LOG_A_SUM_INIT(M,EPS))
     IF((AIMAG(A).EQ.ZERO).AND.(AIMAG(B).EQ.ZERO).AND.     &
          (AIMAG(C).EQ.ZERO)) THEN
        A_TERM=REAL(A_TERM,PR)
     ENDIF
  ENDIF
  A_SUM=A_TERM
  POW_Z_INV_M=Z_INV**M
  B_TERM=CWFN_MIST_B_SUM_INIT_PS_INFINITY(A,C,GAMMA_C,GAMMA_INV_CMA, &
       GAMMA_INV_ONE_MEPS,GAMMA_INV_EPS_PA_PM,Z,M,EPS)*POW_Z_INV_M
  PROD_B=POW_Z_INV_M
  DO N=0,M_M2
     RATIO=(A+N)*(A_MC_P1+N)/(N+ONE)
     A_TERM = A_TERM*Z_INV*RATIO/(N+ONE_MEPS_MM)
     A_SUM = A_SUM+A_TERM
     PROD_B = PROD_B*RATIO
  ENDDO
  IF (M.GT.0) THEN
     PROD_B=PROD_B*(A+M-ONE)*(A_MC_P1+M-ONE)/DBLE(M)
  ENDIF
  B_EXTRA_TERM = PROD_B*GAMMA_PROD*GAMMA_INV_ONE_MEPS
  B_SUM=B_TERM
  B_PREC=EPS15*INF_NORM(B_TERM)
  CALL CWFN_MIST_CV_POLY_DER_TAB_CALC(A,A_MC_P1,ONE_MEPS_MM,Z_INV, &
       CV_POLY1_DER_TAB)
  CALL CWFN_MIST_CV_POLY_DER_TAB_CALC(B,EPS_PA_MC_P1_PM,EPS_PM_P1, &
       Z_INV,CV_POLY2_DER_TAB)
  MIN_N=MAX(CWFN_MIST_MIN_N_CALC(CV_POLY1_DER_TAB),CWFN_MIST_MIN_N_CALC(CV_POLY2_DER_TAB))
  POSSIBLE_FALSE_CV=.TRUE.; N=0
  DO WHILE(POSSIBLE_FALSE_CV.OR.(INF_NORM(B_TERM).GT.B_PREC))
     N_PM_P1=N+M_P1; N_P1=N+ONE; A_PM_PN=A_PM+N
     A_MC_P1_PM_PN=A_MC_P1_PM+N; EPS_PM_P1_PN=EPS_PM_P1+N
     N_P1_MEPS=N_P1-EPS; EPS_PA_PM_PN=EPS_PA_PM+N
     EPS_PA_MC_P1_PM_PN=EPS_PA_MC_P1_PM+N; EPS_PM_PN=EPS_PM+N
     PROD1=EPS_PA_PM_PN*EPS_PA_MC_P1_PM_PN; PROD2=EPS_PM_P1_PN*N_P1
     PROD3=A_PM_PN*A_MC_P1_PM_PN
     B_TERM = Z_INV*(B_TERM*PROD1/PROD2+B_EXTRA_TERM*(PROD3/N_PM_P1 &
          -A_PM_PN-A_MC_P1_PM_PN-EPS+PROD1/N_P1)                    &
          /(EPS_PM_P1_PN*N_P1_MEPS))
     B_SUM=B_SUM+B_TERM
     B_EXTRA_TERM=B_EXTRA_TERM*Z_INV*PROD3/(N_PM_P1*N_P1_MEPS)
     IF(POSSIBLE_FALSE_CV.AND.(N.GT.MIN_N)) THEN
        POSSIBLE_FALSE_CV = (CWFN_MIST_CV_POLY_DER_CALC( &
             CV_POLY1_DER_TAB,DBLE(N)).GT.ZERO).OR.(&
             CWFN_MIST_CV_POLY_DER_CALC(CV_POLY2_DER_TAB,DBLE(N)).GT.ZERO)
     ENDIF
     N=N+1
  ENDDO
  IF(EPS.EQ.ZERO) THEN
     CWFN_MIST_HYP_PS_INFINITY=PHASE*POW_MZ_MA*(A_SUM+B_SUM)
     RETURN
  ELSE
     CWFN_MIST_HYP_PS_INFINITY=PHASE*POW_MZ_MA*(A_SUM+B_SUM)*PI_EPS &
          /SIN(PI_EPS)
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_HYP_PS_INFINITY
!
!----------------------------------------------------------------------
! Calculation of F(z) in transformation theory missing zones 
! ----------------------------------------------------------
! of the complex plane with a Taylor series
! -----------------------------------------
! If z is close to exp(+/- i.pi/3), no transformation in 1-z, z, 
! z/(z-1) or combination of them can transform z in a complex number 
! of modulus smaller than a given Rmax < 1 .
! Rmax is a radius for which one considers power series summation 
! for |z| > Rmax is too slow to be processed. One takes Rmax = 0.9 .
! Nevertheless, for Rmax = 0.9, 
! these zones are small enough to be handled 
! with a Taylor series expansion around a point z0 close to z 
! where transformation theory can be used to calculate F(z).
! One then chooses z0 to be 0.9 z/|z| if |z| < 1, and 1.1 z/|z| 
! if |z| > 1, 
! so that hyp_PS_zero or hyp_PS_infinity can be used 
! (see comments of these functions above).
! For this z0, F(z) = \sum_{n=0}^{+oo} q[n] (z-z0)^n, with:
! q[0] = F(z0), q[1] = F'(z0) = (a b/c) 2F1(a+1,b+1,c+1,z0)
! q[n+2] = [q[n+1] (n (2 z0 - 1) - c + (a+b+c+1) z0) 
! + q[n] (a+n)(b+n)/(n+1)]/(z0(1-z0)(n+2))
! As |z-z0| < 0.1, it converges with around 15 terms, 
! so that no instability can occur for moderate a, b and c.
! Convergence is tested 
! with |q[n] (z-z0)^n|oo + |q[n+1] (z-z0)^{n+1}|oo. 
! Series is truncated when this test is smaller 
! than 1E-15 (|q[0]|oo + |q[1] (z-z0)|oo).
! No false convergence can happen here 
! as q[n] behaves smoothly for n -> +oo.
!
! Variables
! ---------
! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
! abs_z,is_abs_z_small: |z|, true if |z| < 1 and false if not.
! z0,zc_z0_ratio,z0_term1,z0_term2: 0.9 z/|z| if |z| < 1, 
! and 1.1 z/|z| if |z| > 1, (z-z0)/(z0 (1-z0)), 
! 2 z0 - 1, c - (a+b+c+1) z0
! hyp_PS_z0,dhyp_PS_z0,prec: F(z0), F'(z0) calculated with 2F1 
! as F'(z0) = (a b/c) 2F1(a+1,b+1,c+1,z0), 
! precision demanded for series truncation 
! equal to 1E-15 (|q[0]|oo + |q[1] (z-z0)|oo).
! n,an,anp1,anp2,sum: index of the series, q[n] (z-z0)^n, 
! q[n+1] (z-z0)^{n+1}, q[n+2] (z-z0)^{n+2}, 
! truncated sum of the power series.
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_HYP_PS_COMPLEX_PLANE_REST(A,B,C,Z)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  INTEGER(IPR) :: N
  REAL(PR)     :: ABS_Z,PREC
  COMPLEX(PR)  :: CWFN_MIST_HYP_PS_COMPLEX_PLANE_REST
  COMPLEX(PR)  :: Z0,ZC,ZC_Z0_RATIO,Z0_TERM1,Z0_TERM2
  COMPLEX(PR)  :: HYP_PS_Z0,DHYP_PS_Z0,AN,ANP1,ANP2
  COMPLEX(PR)  :: CWFN_MIST_HYP_PS_ZERO,CWFN_MIST_HYP_PS_INFINITY
  !
  ABS_Z=ABS(Z)
  IF(ABS_Z.LT.ONE) THEN
     Z0=0.9D0*Z/ABS_Z; ZC=Z-Z0; ZC_Z0_RATIO=ZC/(Z0*(ONE-Z0))
     Z0_TERM1=TWO*Z0 - ONE; Z0_TERM2=C-(A+B+ONE)*Z0
     HYP_PS_Z0=CWFN_MIST_HYP_PS_ZERO(A,B,C,Z0)
     DHYP_PS_Z0=CWFN_MIST_HYP_PS_ZERO(A+ONE,B+ONE,C+ONE,Z0)*A*B/C 
  ELSE
     Z0=1.1D0*Z/ABS_Z; ZC=Z-Z0; ZC_Z0_RATIO=ZC/(Z0*(ONE-Z0))
     Z0_TERM1=TWO*Z0 - ONE; Z0_TERM2=C-(A+B+ONE)*Z0
     HYP_PS_Z0=CWFN_MIST_HYP_PS_INFINITY(A,B,C,Z0)
     DHYP_PS_Z0=CWFN_MIST_HYP_PS_INFINITY(A+ONE,B+ONE,C+ONE,Z0)*A*B/C 
  ENDIF
  AN=HYP_PS_Z0;ANP1=ZC*DHYP_PS_Z0;CWFN_MIST_HYP_PS_COMPLEX_PLANE_REST=AN+ANP1
  PREC=EPS15*(INF_NORM(AN)+INF_NORM(ANP1)); N=0
  DO WHILE(INF_NORM(AN).GT.PREC)
     ANP2=ZC_Z0_RATIO*(ANP1*(N*Z0_TERM1-Z0_TERM2)+AN*ZC*(A+N)*(B+N) &
          /(N+ONE))/(N+TWO)
     CWFN_MIST_HYP_PS_COMPLEX_PLANE_REST = CWFN_MIST_HYP_PS_COMPLEX_PLANE_REST + ANP2
     N=N+1
     AN=ANP1
     ANP1=ANP2
  ENDDO
  RETURN
END FUNCTION CWFN_MIST_HYP_PS_COMPLEX_PLANE_REST
!
!----------------------------------------------------------------------
! Calculation of F(z) for arbitrary z using previous routines
! -----------------------------------------------------------
! Firstly, it is checked if a,b and c are negative integers.
! If neither a nor b is negative integer but c is, 
! F(z) is undefined so that the program stops with an error message.
! If a and c are negative integers with c < a, 
! or b and c are negative integers with b < a, 
! or c is not negative integer integer but a or b is, 
! one is in the polynomial case.
! In this case, if |z| < |z/(z-1)| or z = 1, 
! hyp_PS_zero is used directly, as then |z| <= 2 
! and no instability arises with hyp_PS_zero 
! as long the degree of the polynomial is small (<= 10 typically).
! If not, one uses the transformation 
! F(z) = (1-z)^{-a} 2F1(a,c-b,c,z/(z-1)) if a is negative integer 
! or F(z) = (1-z)^{-b} 2F1(b,c-a,c,z/(z-1)) if b is negative integer 
! along with hyp_PS_zero.
! Indeed, 2F1(a,c-b,c,X) is a polynomial if a is negative integer, 
! and so is 2F1(b,c-a,c,X) if b is negative integer, 
! so that one has here |z/(z-1)| <= 2 
! and the stability of the method is the same 
! as for the |z| < |z/(z-1)| case.
! If one is in the non-polynomial case, one checks if z >= 1. 
! If it is, one is the cut of F(z) 
! so that z is replaced by z - 10^{-307}i.
! Then, using F(z) = 2F1(b,a,c,z) 
! and F(z) = (1-z)^{c-a-b} 2F1(c-a,c-b,c,z), 
! one replaces a,b,c parameters by combinations of them 
! so that Re(b-a) >= 0 and Re(c-a-b) >= 0.
! Exchanging a and b does not change convergence properties, 
! while having Re(c-a-b) >= 0 accelerates it 
! (In hyp_PS_zero, t[n] z^n ~ z^n/(n^{c-a-b}) for n -> +oo).
! If |1-z| < 1E-5, one uses hyp_PS_one 
! as the vicinity of the singular point z = 1 is treated properly.
! After that, one compares |z| and |z/(z-1)| 
! to R in {0.5,0.6,0.7,0.8,0.9}. 
! If one of them is smaller than R, 
! one uses hyp_PS_zero without transformation
! or with the transformation F(z) = (1-z)^{-a} 2F1(a,c-b,c,z/(z-1)).
! Then, if both of them are larger than 0.9, 
! one compares |1/z|, |(z-1)/z|, |1-z| and |1/(1-z)| 
! to R in {0.5,0.6,0.7,0.8,0.9}. 
! If one of them is found smaller than R, 
! with the condition that |c-b|oo < 5 for (z-1)/z transformation, 
! |a,b,c|oo < 5 for |1-z| transformation 
! and |a,c-b,c|oo < 5 for |1/(1-z)| transformation,
! the corresponding transformation is used. 
! If none of them was smaller than 0.9, 
! one is in the missing zones of transformation theory 
! so that the Taylor series of hyp_PS_complex_plane_rest is used.
!
! Variables
! ---------
! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
! Re_a,Re_b,Re_c,na,nb,nc,is_a_neg_int,is_b_neg_int,is_c_neg_int: 
! real parts of a,b,c, closest integers to a,b,c, 
! true if a,b,c is negative integers and false if not.
! zm1,z_over_zm1,z_shift: z-1, z/(z-1), z - 10^{-307}i in case z >= 1.
! ab_condition, cab_condition: true if Re(b-a) >= 0 and false if not, 
! true if Re(c-a-b) >= 0 and false if not.
! abs_zm1,abz_z,abs_z_inv,abs_z_over_zm1,abs_zm1_inv,abs_zm1_over_z: 
! |z-1|, |z|, |1/z|, |z/(z-1)|, |1/(z-1)|, |(z-1)/z|
! are_ac_small: true if |a|oo < 5 and |c|oo < 5, false if not.
! is_cmb_small: true if |c-b|oo < 5, false if not.
! are_abc_small: true if |a|oo < 5, |b|oo < 5 and |c|oo < 5, 
! false if not.
! are_a_cmb_c_small: true if |a|oo < 5, |c-b|oo < 5 and |c|oo < 5, 
! false if not.
! R_tab,R: table of radii {0.5,0.6,0.7,0.8,0.9}, one of these radii.
! res: returned result
!----------------------------------------------------------------------
RECURSIVE FUNCTION cwfnMiSt2F1(A,B,C,Z) RESULT(RES)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  INTEGER(IPR) :: NA,NB,NC,I
  REAL(PR)     :: RE_A,RE_B,RE_C,ABS_Z,ABS_ZM1,ABS_Z_OVER_ZM1
  REAL(PR)     :: ABS_ZM1_OVER_Z,ABS_ZM1_INV,R_TABLE(1:5),R,ABS_Z_INV
  COMPLEX(PR)  :: RES,CWFN_MIST_HYP_PS_INFINITY,CWFN_MIST_HYP_PS_ZERO,Z_SHIFT
  COMPLEX(PR)  :: CWFN_MIST_HYP_PS_COMPLEX_PLANE_REST,CWFN_MIST_HYP_PS_ONE,Z_OVER_ZM1,ZM1
  LOGICAL      :: IS_A_NEG_INT,IS_B_NEG_INT,IS_C_NEG_INT
  LOGICAL      :: AB_CONDITION,CAB_CONDITION,ARE_A_CMB_C_SMALL
  LOGICAL      :: IS_CMB_SMALL,ARE_AC_SMALL,ARE_ABC_SMALL
  !
  RE_A=REAL(A,PR); RE_B=REAL(B,PR); RE_C=REAL(C,PR);
  NA=NINT(RE_A); NB=NINT(RE_B); NC=NINT(RE_C);
  IS_A_NEG_INT=A.EQ.NA.AND.NA.LE.0
  IS_B_NEG_INT=B.EQ.NB.AND.NB.LE.0
  IS_C_NEG_INT=C.EQ.NC.AND.NC.LE.0
  ZM1=Z-ONE
  IF(IS_C_NEG_INT) THEN
     ABS_Z=ABS(Z); Z_OVER_ZM1 = Z/ZM1
     ABS_Z_OVER_ZM1=ABS(Z_OVER_ZM1)
     IF(IS_A_NEG_INT.AND.(NC.LT.NA)) THEN
        IF((Z.EQ.ONE).OR.(ABS_Z.LT.ABS_Z_OVER_ZM1)) THEN
           RES=CWFN_MIST_HYP_PS_ZERO(A,B,C,Z)
           RETURN
        ELSE
           RES=((-ZM1)**(-A))*CWFN_MIST_HYP_PS_ZERO(A,C-B,C,Z_OVER_ZM1)
           RETURN
        ENDIF
     ELSE IF(IS_B_NEG_INT.AND.(NC.LT.NB)) THEN
        IF((Z.EQ.ONE).OR.(ABS_Z.LT.ABS_Z_OVER_ZM1)) THEN
           RES=CWFN_MIST_HYP_PS_ZERO(A,B,C,Z)
           RETURN
        ELSE
           RES=((-ZM1)**(-B))*CWFN_MIST_HYP_PS_ZERO(B,C-A,C,Z_OVER_ZM1)
           RETURN
        ENDIF
     ELSE
        STOP '2F1 UNDEFINED'
     ENDIF
  ENDIF
  IF(IS_A_NEG_INT) THEN
     ABS_Z=ABS(Z); Z_OVER_ZM1 = Z/ZM1
     ABS_Z_OVER_ZM1=ABS(Z_OVER_ZM1)
     IF((Z.EQ.ONE).OR.(ABS_Z.LT.ABS_Z_OVER_ZM1)) THEN
        RES=CWFN_MIST_HYP_PS_ZERO(A,B,C,Z)
        RETURN
     ELSE
        RES=((-ZM1)**(-A))*CWFN_MIST_HYP_PS_ZERO(A,C-B,C,Z_OVER_ZM1)
        RETURN
     ENDIF
  ELSE IF(IS_B_NEG_INT) THEN
     ABS_Z=ABS(Z); Z_OVER_ZM1 = Z/ZM1
     ABS_Z_OVER_ZM1=ABS(Z_OVER_ZM1)
     IF((Z.EQ.ONE).OR.(ABS_Z.LT.ABS_Z_OVER_ZM1)) THEN
        RES=CWFN_MIST_HYP_PS_ZERO(A,B,C,Z)
        RETURN
     ELSE
        RES=((-ZM1)**(-B))*CWFN_MIST_HYP_PS_ZERO(B,C-A,C,Z_OVER_ZM1)
        RETURN
     ENDIF
  ENDIF
  IF((REAL(Z,PR).GE.ONE).AND.(AIMAG(Z).EQ.ZERO)) THEN
     Z_SHIFT=CMPLX(REAL(Z,PR),-1.0D-307,PR)
     RES=CWFNMIST2F1(A,B,C,Z_SHIFT)
     RETURN
  ENDIF
  AB_CONDITION = (RE_B.GE.RE_A)
  CAB_CONDITION = (RE_C.GE.RE_A + RE_B)
  IF ((.NOT.AB_CONDITION).OR.(.NOT.CAB_CONDITION)) THEN
     IF ((.NOT.AB_CONDITION).AND.(CAB_CONDITION)) THEN
        RES=CWFNMIST2F1(B,A,C,Z)
        RETURN
     ELSE IF((.NOT.CAB_CONDITION).AND.(AB_CONDITION)) THEN
        RES=((-ZM1)**(C-A-B))*CWFNMIST2F1(C-B,C-A,C,Z)
        RETURN
     ELSE
        RES=((-ZM1)**(C-A-B))*CWFNMIST2F1(C-A,C-B,C,Z)
        RETURN 
     ENDIF
  ENDIF
  ABS_ZM1=ABS(ZM1)
  IF(ABS_ZM1.LT.1D-5) THEN 
     RES=CWFN_MIST_HYP_PS_ONE (A,B,C,-ZM1)
     RETURN
  ENDIF
  ABS_Z=ABS(Z); ABS_Z_OVER_ZM1=ABS_Z/ABS_ZM1; ABS_Z_INV=ONE/ABS_Z
  ABS_ZM1_OVER_Z=ONE/ABS_Z_OVER_ZM1; ABS_ZM1_INV=ONE/ABS_ZM1
  IS_CMB_SMALL = INF_NORM(C-B).LT.5.0D0; 
  ARE_AC_SMALL = (INF_NORM(A).LT.5.0D0).AND.(INF_NORM(C).LT.5.0D0)
  ARE_ABC_SMALL = ARE_AC_SMALL.AND.(INF_NORM(B).LT.5.0D0)
  ARE_A_CMB_C_SMALL = ARE_AC_SMALL.AND.IS_CMB_SMALL
  R_TABLE=(/0.5D0,0.6D0,0.7D0,0.8D0,0.9D0/)
  DO I=1,5
     R=R_TABLE(I)
     IF(ABS_Z.LE.R) THEN 
        RES=CWFN_MIST_HYP_PS_ZERO (A,B,C,Z)
        RETURN
     ENDIF
     IF(IS_CMB_SMALL.AND.(ABS_Z_OVER_ZM1.LE.R)) THEN
        RES=((-ZM1)**(-A))*CWFN_MIST_HYP_PS_ZERO (A,C-B,C,Z/ZM1)
        RETURN
     ENDIF
  ENDDO
  DO I=1,5
     R=R_TABLE(I)
     IF(ABS_Z_INV.LE.R) THEN 
        RES=CWFN_MIST_HYP_PS_INFINITY (A,B,C,Z)
        RETURN 
     ENDIF
     IF(IS_CMB_SMALL.AND.(ABS_ZM1_OVER_Z.LE.R)) THEN 
        RES=((-ZM1)**(-A))*CWFN_MIST_HYP_PS_INFINITY (A,C-B,C,Z/ZM1)
        RETURN
     ENDIF
     IF(ARE_ABC_SMALL.AND.(ABS_ZM1.LE.R)) THEN 
        RES=CWFN_MIST_HYP_PS_ONE (A,B,C,-ZM1)
        RETURN
     ENDIF
     IF(ARE_A_CMB_C_SMALL.AND.(ABS_ZM1_INV.LE.R)) THEN 
        RES=((-ZM1)**(-A))*CWFN_MIST_HYP_PS_ONE (A,C-B,C,-ONE/ZM1)
        RETURN
     ENDIF
  ENDDO
  RES=CWFN_MIST_HYP_PS_COMPLEX_PLANE_REST (A,B,C,Z)
  RETURN
END FUNCTION cwfnMiSt2F1
!
!----------------------------------------------------------------------
! Test of 2F1 numerical accuracy 
! ------------------------------
! using hypergeometric differential equation
! ------------------------------------------
! If z = 0, F(z) = 1 so that this value is trivially tested.
! To test otherwise if the value of F(z) is accurate, 
! one uses the fact that 
! z(z-1) F''(z) + (c - (a+b+1) z) F'(z) - a b F(z) = 0.
! If z is not equal to one, a relative precision test is provided 
! by |F''(z) + [(c - (a+b+1) z) F'(z) - a b F(z)]/[z(z-1)]|oo
! /(|F(z)|oo + F'(z)|oo + |F''(z)|oo).
! If z is equal to one, one uses |(c - (a+b+1)) F'(z) - a b F(z)|oo
! /(|F(z)|oo + F'(z)|oo + 1E-307).
! F'(z) and F''(z) are calculated using equalities 
! F'(z) = (a b/c) 2F1(a+1,b+1,c+1,z) 
! and F'(z) = ((a+1)(b+1)/(c+1)) (a b/c) 2F1(a+2,b+2,c+2,z).
!
! Variables
! ---------
! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
! F,dF,d2F: F(z), F'(z) and F''(z) calculated with hyp_2F1 
! using F'(z) = (a b/c) 2F1(a+1,b+1,c+1,z) 
! and F'(z) = ((a+1)(b+1)/(c+1)) (a b/c) 2F1(a+2,b+2,c+2,z).
!----------------------------------------------------------------------
FUNCTION CWFN_MIST_TEST_2F1(A,B,C,Z,F)
  !--------------------------------------------------------------------
  USE CWFN_MIST
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  REAL(PR)    :: CWFN_MIST_TEST_2F1
  COMPLEX(PR) :: F,DF,D2F,CWFNMIST2F1  
  !
  IF(Z.EQ.ZERO) THEN
     CWFN_MIST_TEST_2F1=INF_NORM(F-ONE)
     RETURN
  ELSE IF(Z.EQ.ONE) THEN
     DF = CWFNMIST2F1(A+ONE,B+ONE,C+ONE,Z)*A*B/C
     CWFN_MIST_TEST_2F1=INF_NORM((C-(A+B+ONE))*DF-A*B*F) &
          /(INF_NORM (F)+INF_NORM(DF)+1D-307)
     RETURN
  ELSE
     DF = CWFNMIST2F1(A+ONE,B+ONE,C+ONE,Z)*A*B/C
     D2F = CWFNMIST2F1(A+TWO,B+TWO,C+TWO,Z)*A*(A+ONE)*B*(B+ONE) &
          /(C*(C+ONE))
     CWFN_MIST_TEST_2F1=INF_NORM(D2F+((C-(A+B+ONE)*Z)*DF-A*B*F)/(Z*(ONE-Z))) &
          /(INF_NORM(F)+INF_NORM(DF)+INF_NORM(D2F))
     RETURN
  ENDIF
END FUNCTION CWFN_MIST_TEST_2F1
! double precision function cwfnnormphase(phase)
! Computes normalized phase, i.e. phase between -\pi and \pi.
!
! 2013-08-05 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

double precision function cwfnnormphase(phase)
    implicit none

    double precision, intent(in) :: phase

    ! Other variables
    double precision :: pi
    double precision :: dblpi

    double precision :: ph
    double precision :: tmp

    ! Working...

    pi = 4.d0*datan(1.d0)
    dblpi = 2.d0 * pi

    ph = phase

    ! Eliminating (2 \pi n)
    tmp = abs(ph)
    if(tmp.ge.dblpi) then
        tmp = floor( tmp / dblpi )
        tmp = tmp * dblpi
        if(ph.lt.0.d0) then
            tmp = -tmp
        end if
        ph = ph - tmp
    end if

    ! Eliminating \pi
    tmp = abs(ph)
    if(tmp.gt.pi) then
        tmp = tmp - pi
        if(ph.gt.0.d0) then
            ph = -pi + tmp
        else
            ! ph < 0, but |ph| > \pi
            ph = pi - tmp
        end if
    end if

    cwfnnormphase = ph
end function cwfnnormphase
! Subroutine cwfnomega(p, q, l, eta, value)
! Computes \Omega(p, q, l, \eta) from Eq.(3.8) at libcwfn-doc.pdf.
!
! 2013-08-05 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

!Charlie: Changed value to real
subroutine cwfnomega(p, q, l, eta, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double precision, intent(out) :: value
    ! Other variables
    double precision :: p2
    double precision :: q2
    double complex :: nom
    double complex :: denom
    double complex :: varphi

    double precision :: tmp
    double complex :: ctmp,valuec

    ! Working...
    ! The factor in front of \varphi()
    ! The nominator
    call cwfntau(p, q, eta, nom)
    tmp = dble(l)
    nom = nom * p**tmp

    ! The denominator
    ctmp = dcmplx(1.d0 + dble(l), eta)
    p2 = p * p
    q2 = q * q
    tmp = p2 + q2
    denom = dcmplx(tmp, 0.d0)
    denom = denom**ctmp
        
    if (nom.eq.0.d0) then
    varphi=1.d0
    else
    ! \varphi()
    call cwfnvarphi(p, q, l, eta, varphi)
    endif

    ! Combining all together
    valuec = nom * varphi / denom
    value=real(valuec)
end subroutine cwfnomega
! Subroutine cwfnphi(a, b, c, z, value)
! Computes \phi(a, b, c, z) from Eq.(3.13) at libcwfn-doc.pdf.
!
! 2013-08-01 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnphi(a, b, c, z, value)
    implicit none

    double complex, intent(in) :: a
    double complex, intent(in) :: b
    double complex, intent(in) :: c
    double complex, intent(in) :: z
    double complex, intent(out) :: value

    ! Working..
    ! Computing hypergeometric 2F1()
    call cwfnf21(a, b, c, z, value)
end subroutine cwfnphi
! Subroutine cwfnphitilde(a, b, c, h, value)
! Computes \widetilde{\phi}(a, b, c, h) from Eq.(3.14) at libcwfn-doc.pdf.
!
! In other words, computes hypergeometric function, applying linear
! transformation.
! See libcwfn-doc.pdf for details.
!
! 2013-08-05 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnphitilde(a, b, c, h, value)
    implicit none

    ! Interface variables
    double complex, intent(in) :: a
    double complex, intent(in) :: b
    double complex, intent(in) :: c
    double complex, intent(in) :: h
    double complex, intent(out) :: value

    ! Other variables
    double precision :: tmp
    double complex :: ctmp
    double complex :: cab
    double complex :: abc

    double precision :: phase
    double complex :: f21
    double complex :: lng1
    double complex :: lng2
    double complex :: lng3
    double complex :: lng4
    double complex :: fct
    double complex :: brackets
    double complex :: hpow

    interface
        double precision function cwfnnormphase(ph)
            double precision, intent(in) :: ph
        end function cwfnnormphase
    end interface

    ! Working...
    ! Some preparations...
    abc = a + b - c

    ! The front factor
    ! [(\Gamma(c) * \Gamma(a + b - c)) / (\Gamma(a) \Gamma(b))]
    ! * F(a, b, a + b - c + 1, h)
    ctmp = abc + dcmplx(1.d0, 0.0)
    call cwfnf21(a, b, ctmp, h, f21)

    ctmp = c
    call cwfncgammaln(ctmp, lng1)

    ctmp = abc
    call cwfncgammaln(ctmp, lng2)

    ctmp = a
    call cwfncgammaln(ctmp, lng3)

    ctmp = b
    call cwfncgammaln(ctmp, lng4)

    ctmp = lng1 + lng2 - lng3 - lng4

    ctmp = exp(ctmp)
    fct = f21 * ctmp

    ! The brackets.
    ! The first term
    ! \exp\{ 2 i [ -\arg(\Gamma(a + b - c)) + \arg(\Gamma(a)) + \arg(\Gamma(b)) ] \}
    phase = -dimag(lng2) + dimag(lng3) + dimag(lng4)
    phase = 2.d0 * phase
    phase = cwfnnormphase(phase)

    ctmp = dcmplx(0.d0, phase)
    brackets = exp(ctmp)

    ! The second term
    ! h^{c - a - b} * \exp[ -2 i \arg(F()) ]
    phase = -2.d0 * dimag(log(f21))
    phase = cwfnnormphase(phase)

    ctmp = dcmplx(0.d0, phase)
    ctmp = exp(ctmp)
    ! h^{c - a - b}
    cab = c - a - b
    hpow = h**cab
    brackets = brackets + hpow * ctmp

    value = fct * brackets
end subroutine cwfnphitilde
! Subroutine cwfnrega(q, l, eta, value)
! Computes \mathcal{A}(q, l, \eta) from Eq.(3.24) in libcwfn-doc.pdf.
!
! 2014-05-15 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration


!Charlie: Changed to remove complex phase factor
subroutine cwfnrega(q, l, eta, value)
    implicit none

    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
!    double complex, intent(out) :: value
    double precision, intent(out) :: value

    ! Other variables (and even constants)
    double complex :: icompl = dcmplx(0.d0, 1.d0) ! complex i
    double precision :: pi

    double precision :: tmp
    double precision ::ctmp
!    double complex :: ctmp
    double complex :: z

!    interface
!        double precision function cwfnSigmaL(lt, etat)
!            integer, intent(in) :: lt
!            double precision, intent(in) :: etat
!        end function cwfnSigmaL
!    end interface

    ! Working...
    pi = 4.d0 * datan(1.d0)

    ! Real factor
    tmp = 2.d0 * pi / (4.d0 * q)**l
!    value = dcmplx(tmp, 0.d0)
    value=tmp

    ! Complex factor
!    tmp = cwfnSigmaL(l, eta)
!    ctmp = dcmplx(-pi * eta / 2.d0, tmp)
     ctmp=exp(-pi*eta/2.d0)
!    ctmp = icompl * exp(ctmp)

    value = value * ctmp
end subroutine cwfnrega
! Subroutine cwfnregb(p, q, l, eta, value)
! Computes \mathcal{B}(p, q, l, \eta) from Eqs.(3.25)
! and (3.26) in libcwfn-doc.pdf.
!
! 2014-05-15 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnregb(p, q, l, eta, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double complex, intent(out) :: value

    ! Other variables
    double precision, parameter :: zeroThreshould = 1.e-13

    double complex :: nom
    double complex :: denom

    double precision :: rho
    double complex :: a
    double complex :: b
    double complex :: c
    double complex :: z
    double complex :: f21

    double precision :: pq
    double precision :: tmp
    double complex :: ctmp
    double complex :: g

    interface
        double precision function cwfnSigmaL(lt, etat)
            integer, intent(in) :: lt
            double precision, intent(in) :: etat
        end function cwfnSigmaL
    end interface

    ! Working...
    ! Factors in front of all
    ! \Gamma(1 + i \eta) e^{-i \sigma_l}
    ctmp = dcmplx(1.d0, eta)
    call cwfncgammaln(ctmp, g)
    tmp = cwfnSigmaL(l, eta)
    ctmp = dcmplx(0.d0, -tmp)
    value = exp(g + ctmp)

    ! Is p = q?
    call cwfnrho(p, q, rho)
    if(rho < zeroThreshould) then
        ! p = q
        ! The simplified expression from Eq.(3.26) is used
        ! 2^{2l - 1 + i \eta}
        tmp = dble(2 * l - 1)
        ctmp = dcmplx(tmp, eta)
        value = value * dcmplx(2.d0, 0.d0)**ctmp

        ! q^{l - 2 + i \eta}
        tmp = dble(l - 2)
        ctmp = dcmplx(tmp, eta)
        value = value * dcmplx(q, 0.d0)**ctmp

        return
    end if

    ! Full-scale case, Eq.(3.24)

    ! Computing hypergeometric 2F1()
    ! Constructing arguments
    a = dcmplx(-dble(l), 0.d0)
    b = dcmplx(-dble(l), -eta)
    c = dcmplx(1.d0, -eta)
    z = dcmplx(rho, 0.d0)

    ! Computing 2F2() itself
    call cwfnf21(a, b, c, z, f21)

    ! Computing the factor in front of 2F1()
    ! The nominator
    g = dcmplx(p + q, 0.d0)
    ctmp = dcmplx(dble(2 * l - 1), eta)
    nom = g**ctmp

    ! The denominator
    tmp = dble(l) + 1.d0
    denom = p**tmp

    ! Combining all together
    value = value * nom * f21 / denom
end subroutine cwfnregb
! Subroutine cwfnrho(p, q, value)
! Computes the 4-th argument of the hypergeometric 2F1()
! for the `pole proximity' representation of Coulomb WF.
! See Eq.(3.22) at libcwfn-doc.pdf.
!
! 2013-07-05 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnrho(p, q, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    double precision, intent(out) :: value

    ! Other variables
    double precision :: nom
    double precision :: denom

    ! Working...
    nom = p - q
    nom = nom * nom
    denom = p + q
    denom = denom * denom

    value = nom / denom
end subroutine cwfnrho
! double precision function cwfnSigmaL(l, eta)
! Computes Coulomb phase shift \sigma_l (normalized) using
! Eq.(33.2.10) in http://dlmf.nist.gov:
! \sigma_l = \arg\Gamma(1 + l + i \eta) = \Im\ln\Gamma(1 + l + i \eta).
!
! 2014-05-15 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

double precision function cwfnSigmaL(l, eta)
    implicit none

    integer, intent(in) :: l
    double precision, intent(in) :: eta

    ! Other variables
    double complex :: z
    double complex :: lng

    interface
        double precision function cwfnnormphase(ph)
            double precision, intent(in) :: ph
        end function cwfnnormphase
    end interface

    ! Working...
    ! cwfnSigmaL = Im{ln[Gamma(l + 1 + i \eta)]}
    z = dcmplx(dble(l + 1), eta)
    call cwfncgammaln(z, lng)
    cwfnSigmaL = dimag(lng)
    cwfnSigmaL = cwfnnormphase(cwfnSigmaL)
end function cwfnSigmaL
! Subroutine cwfntau(p, q, eta, value)
! Computes \tau(p, q, \eta) using Eq.(3.9) at libcwfn-doc.pdf.
!
! 2013-07-10 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

!Charlie: Moved exp(-pi*eta/2.0) from upsilon to here to calculate low eta
!values
!Update: Further moved exp(pi*|eta|/2.0) to gamma_mag
subroutine cwfntau(p, q, eta, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    double precision, intent(in) :: eta
    double complex, intent(out) :: value

    ! Other variables and constants
    double precision :: pi
    double precision :: a
    double precision :: phaseFactor
    double complex :: tmp

    ! Working...
    pi = 4.d0 * datan(1.d0)

    if(abs(p - q) < 1.d-7) then
        write(0, *) "cwfntau(): Coulomb wave function does not exist at p = q!"
        stop 1
    end if

    a = abs(p*p - q*q)

    ! Computing phase shift factor
    if(p.gt.q) then
        if (eta.gt.0.d0) then
        phaseFactor = 1.d0 * exp(-pi * eta)
        else
        phaseFactor = 1.d0! * exp(-pi * eta/2.d0)
        endif
    else
        if (eta.gt.0.d0) then
        phaseFactor = -1.d0! * exp(pi * eta/2.d0)
        else
        phaseFactor = -1.d0 * exp(pi*eta)
        endif
    end if

    ! Computing tau
    value = dcmplx(a, 0.d0)
    tmp = dcmplx(-1.d0, eta)
    value = value**tmp
    value = value * phaseFactor
end subroutine cwfntau
! Subroutine cwfntheta(q, l, eta, value)
! Computes \Theta(q, l, \eta) from the eq.(3.19) at libcwfn-doc.pdf.
!
! 2014-05-15 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

!Charlie: Changes made to remove sigmaL
!       removed exp(-pieta/2) term and moved to xi
subroutine cwfntheta(q, l, eta, value)
    implicit none

    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
!    double complex, intent(out) :: value
    double precision, intent(out) :: value
    ! Other variables (and even constants)
    double precision :: pi
    double precision :: tmp
    double precision :: sigma
!    double complex :: ctmp
    double precision :: ctmp
!    interface
!        double precision function cwfnSigmaL(lt, etat)
!            integer, intent(in) :: lt
!            double precision, intent(in) :: etat
!        end function cwfnSigmaL
!    end interface


    ! Working...
    pi = 4.d0 * datan(1.d0)
!    sigma = cwfnSigmaL(l, eta)

    ! Real factor
    tmp = 4.d0
    tmp = -(tmp**(1-l)) * pi / q**l

    ! Complex exponent
!    ctmp = dcmplx(-pi * eta / 2.d0, sigma)
!    ctmp = exp(ctmp)
!     ctmp=exp(-pi*eta/2.d0)

    ! Combining
!    value = dcmplx(tmp, 0.d0) * ctmp
    value=tmp!*ctmp
end subroutine cwfntheta
! Subroutine cwfnupsilon(q, l, eta, value)
! Computes \Upsilon(q, l, \eta) from Eq.(3.7) at libcwfn-doc.pdf.
!
! 2014-04-01 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

!Charlie: Gamma subroutine replaced to just calculate magnitude
!Moved exp(-pi*eta/2) to tau
subroutine cwfnupsilon(q, l, eta, value)
    implicit none

    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
!    double complex, intent(out) :: value
    double precision,intent(out)::value

    ! Other variables (and even constants)
    integer :: i
    double precision :: pi
    double precision :: tmp
!    double complex :: ctmp
    double precision :: ctmp
    double complex :: z

    ! Working...
    pi = 4.d0 * datan(1.d0)

    ! Everything outside of the square brackets
    tmp = -4.d0 * pi * eta
    tmp = tmp * q**(l+1)!*exp(-pi*eta/2.d0)
    !value = dcmplx(tmp, 0.d0)
    value=tmp

    ! Square brackets
    ! Gamma(z)
!    z = dcmplx(dble(l + 1), eta)
!    call cwfncgammaln(z, ctmp)
!    ctmp = exp(ctmp)

     call cwf_gamma_mag(l,eta,ctmp)
     !print*,"Gamma",ctmp
    ! Inverted Pochhammer's symbol 1 / (1/2)_{l + 1}
    tmp = 0.5d0
    do i = 1, l
       tmp = tmp * (0.5d0 + i)
    end do
    tmp = 1.d0 / tmp
        
        !print*,tmp,ctmp,value
    ! Combining together
    !ctmp = ctmp * dcmplx(tmp, 0.d0)
    ctmp=ctmp*tmp
    ! The last operation...
    value = value * ctmp
end subroutine cwfnupsilon
! Subroutine cwfnvarphi(p, q, l, eta, value)
! Computes \varphi(p, q, l, \eta) from Eq.(3.10) at libcwfn-doc.pdf.
!
! 2013-08-05 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnvarphi(p, q, l, eta, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double complex, intent(out) :: value

    ! Other variables
    double precision :: zeta
    double precision :: chi
    double complex :: a
    double complex :: b
    double complex :: c

    double complex :: ctmp

    ! Working...
    ! Constructing arguments
    a = dcmplx(1.d0 + dble(l) / 2.d0, eta / 2.d0)
    b = dcmplx(0.5d0 + dble(l) / 2.d0, eta / 2.d0)
    c = dcmplx(dble(l) + 1.5d0, 0.d0)

    call cwfnzeta(p, q, zeta)
    call cwfnchi(p, q, chi)

    ! Compute \phi() or \widetilde{\phi}() correspondingly
    if(zeta.le.0.5d0) then
        ! Compute \phi()
        ctmp = dcmplx(zeta, 0.d0)
        call cwfnphi(a, b, c, ctmp, value)
    else
        ! Compute \widetilde{\phi}()
        ctmp = dcmplx(chi, 0.d0)
        call cwfnphitilde(a, b, c, ctmp, value)
    end if
end subroutine cwfnvarphi
! Subroutine cwfnWaveFunction(p, q, l, eta, value)
! Computes `full' Coulomb WF \phi(p, q, l, \eta) from the Chapt.3
! at libcwfn-doc.pdf, but checking the applicability of the library
! using the eq.(3.4) from the same paper.
! The execution will be STOPPED, if l or \eta are outside of the
! region of applicability.
! From the other point of view, it is a wrapper around cwfnFull().
!
! 2014-03-28 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnWaveFunction(p, q, l, eta, value)
    implicit none

    ! I/O units
    integer, parameter :: stderr = 0

    ! Interface variables
    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double complex, intent(out) :: value

    ! Parameters to check the applicability using the eq.(3.4--3.5)
    double precision, parameter:: a = 0.001
    double precision, parameter:: b = 100
    double precision, parameter:: f = 0.1

    ! Variables
    double precision :: ratio
    double precision :: dbl ! l converted to the double precision
    double precision :: rhs ! right hand side of the eq.(3.4)


    ! Working...
    dbl = dble(l)

    ! Checks and exceptions
    if(l > 0) then
       ratio = eta / dbl
    else
       ratio = 1.
    end if
    rhs = a * (b + dbl * dbl) - f

    if(l > 5 .and. ratio < rhs) then
        write(stderr, *) &
             "cwfnWaveFunction(): ERROR: L or eta are outside of the"
        write(stderr, *) "region of applicability!"
        write(stderr, *) "cwfnWaveFunction(): l = ", l, ", eta = ", eta
        write(stderr, *) "cwfnWaveFunction(): eta / l = ", ratio
        write(stderr, *) "cwfnWaveFunction(): right hand side (rhs) = ", rhs
        write(stderr, *) "cwfnWaveFunction(): eta / l >= rhs must be satisfied!"
        write(stderr, *) "cwfnWaveFunction(): Execution halted!"
        stop 1
    end if

    ! Compute the value
    call cwfnFull(p, q, l, eta, value)
end subroutine cwfnWaveFunction
! Subroutine cwfnWaveFunctionV(n, pv, q, l, eta, values) Computes `full'
! Coulomb WF \phi(pv(i), q, l, \eta) for every i = 1, 2,... n; results are
! stored at values.  Computations and checks are performed using
! cwfnWaveFunction() subroutine.
! 
! The execution will be STOPPED, if l or \eta are outside of the
! region of applicability.
! From the other point of view, it is a wrapper around cwfnWaveFunction().
!
! 2014-03-31 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnWaveFunctionV(n, pv, q, l, eta, values)
    implicit none

    ! Interface variables
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: pv
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double complex, dimension(n), intent(out) :: values

    ! Other variables
    integer :: i
    double precision :: p
    double complex :: val


    ! Working...
    do i = 1, n
        p = pv(i)
        call cwfnWaveFunction(p, q, l, eta, val)
        values(i) = val
    end do
end subroutine cwfnWaveFunctionV
! Subroutine cwfnxi(p, q, l, eta, value)
! Computes \Xi(p, q, l, \eta) from Eq.(3.21) at libcwfn-doc.pdf.
!
! 2014-05-15 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

!Charlie: Added exp(-pieta/2) from theta
subroutine cwfnxi(p, q, l, eta, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    integer, intent(in) :: l
    double precision, intent(in) :: eta
    double complex, intent(out) :: value

    ! Other variables
    double precision :: rho,phase,pi
    double complex :: nom
    double complex :: denom
    double complex :: a
    double complex :: b
    double complex :: c
    double complex :: z
    double complex :: f21
!    double complex :: g
    double precision :: sigma,sigma0,g

    double complex :: ctmp

    interface
        double precision function cwfnSigmaL(lt, etat)
            integer, intent(in) :: lt
            double precision, intent(in) :: etat
        end function cwfnSigmaL
    end interface
    
    pi = 4.d0 * datan(1.d0)
    ! Working..
    ! Computing hypergeometric 2F1()
    ! Constructing arguments
    a = dcmplx(-l, 0.d0)
    b = dcmplx(-l, -eta)
    c = dcmplx(1.d0, -eta)

    call cwfnrho(p, q, rho)
    z = dcmplx(rho, 0.d0)

    ! Computing 2F2() itself
    call cwfnf21(a, b, c, z, f21)

    ! Computing the factor in front of 2F1()
    ! The nominator
    nom = dcmplx((p + q), 0.d0)
    ctmp = dcmplx(0.d0, eta)
    nom = nom**ctmp
        

    sigma=cwfnSigmaL(l,eta)
    sigma0=cwfnSigmaL(0,eta)
    ctmp=dcmplx(0.d0,sigma0-sigma)
    call cwf_gamma_mag(0,eta,g)
    if (eta.gt.0.d0) then
            if ((p-q).gt.0.d0) then
            phase=1.d0*exp(-pi*eta)
            else
            phase=-1.d0
            endif
    else
            if ((p-q).gt.0.d0) then
            phase=1.d0
            else
            phase=-1.d0*exp(pi*eta)
            endif
           
    endif
    nom = nom * g * exp(ctmp)

    ! The denominator
    denom = dcmplx(dabs(p - q), 0.d0)
    ctmp = dcmplx(1.d0, eta)
    denom = denom**ctmp
    
    ! Combining all together
    value = nom * f21 / denom*phase
end subroutine cwfnxi
! Subroutine cwfnzeta(p, q, value)
! Computes the 4-th argument of the hypergeometric 2F1()
! for the `regular' representation of Coulomb WF.
! See Eq.(3.11) at libcwfn-doc.pdf.
!
! 2013-07-04 Vasily Eremenko <to@vsl.name> in the name of
! TORUS Collaboration

subroutine cwfnzeta(p, q, value)
    implicit none

    double precision, intent(in) :: p
    double precision, intent(in) :: q
    double precision, intent(out) :: value

    ! Other variables
    double precision :: p2
    double precision :: q2
    double precision :: denom

    ! Working...
    p2 = p * p
    q2 = q * q
    denom = p2 + q2
    denom = denom * denom

    value = 4.d0 * q2 * p2 / denom
end subroutine cwfnzeta

! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
recursive subroutine quicksort(a, first, last)
  implicit none
  real*8  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort

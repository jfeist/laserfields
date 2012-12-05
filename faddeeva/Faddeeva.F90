!!$ Copyright (c) 2012, Johannes Feist
!!$ licensed under the MIT open source license, see LICENSE file
!!$ The open-source C++ code by Steven G. Johnson this is based on is
!!$ is provided in Faddeeva.cc and is available
!!$ at: http://ab-initio.mit.edu/Faddeeva
!!$ I have translated it to Fortran 95 in Oct 2012.
!!$ the file faddeeva.f90 used for the main code can be produced
!!$ with "cpp -P -DUSENRTYPE Faadeeva_w.F90 > faddeeva.f90"
!!$   Johannes Feist, Universidad Autonoma de Madrid

!!$   Compute the Faddeeva function w(z) = exp(-z^2) * erfc(-i*z),
!!$   to a desired relative accuracy relerr, for arbitrary complex
!!$   arguments z.
!!$
!!$   For sufficiently large |z|, we use a continued-fraction expansion
!!$   for w(z) similar to those described in:
!!$
!!$      Walter Gautschi, "Efficient computation of the complex error
!!$      function," SIAM J. Numer. Anal. 7(1), pp. 187-198 (1970)
!!$
!!$      G. P. M. Poppe and C. M. J. Wijers, "More efficient computation
!!$      of the complex error function," ACM Trans. Math. Soft. 16(1),
!!$      pp. 38-46 (1990).
!!$
!!$   Unlike those papers, however, we switch to a completely different
!!$   algorithm for smaller |z|:
!!$
!!$      Mofreh R. Zaghloul and Ahmed N. Ali, "Algorithm 916: Computing the
!!$      Faddeyeva and Voigt Functions," ACM Trans. Math. Soft. 38(2), 15
!!$      (2011).
!!$
!!$   (I initially used this algorithm for all z, but it turned out to be
!!$    significantly slower than the continued-fraction expansion for
!!$    larger |z|.  On the other hand, it is competitive for smaller |z|,
!!$    and is significantly more accurate than the Poppe & Wijers code
!!$    in some regions, e.g. in the vicinity of z=1+1i.)
!!$
!!$   Note that this is an INDEPENDENT RE-IMPLEMENTATION of these algorithms,
!!$   based on the description in the papers ONLY.  In particular, I did
!!$   not refer to the authors' Fortran or Matlab implementations, respectively,
!!$   (which are under restrictive ACM copyright terms and therefore unusable
!!$    in free/open-source software).
!!$
!!$   Steven G. Johnson, Massachusetts Institute of Technology
!!$   http://math.mit.edu/~stevenj
!!$   October 2012.
!!$
!!$    -- Note that Algorithm 916 assumes that the erfc(x) function,
!!$       or rather the scaled function erfcx(x) = exp(x**2)*erfc(x),
!!$       is supplied for REAL arguments x. I originally used an
!!$       erfcx routine derived from DERFC in SLATEC, but I have
!!$       since replaced it with a much faster routine written by
!!$       me which uses a combination of continued-fraction expansions
!!$       and a lookup table of Chebyshev polynomials.
!!$
!!$   A small test program is included the end, which checks
!!$   the w(z) results against several known values.  To compile
!!$   the test function, compile with -DFADDEEVA_W_TEST (that is,
!!$   #define FADDEEVA_W_TEST).
!!$
!!$   REVISION HISTORY:
!!$       4 October 2012: Initial public release (SGJ)
!!$       5 October 2012: Revised (SGJ) to fix spelling error,
!!$                       start summation for large x at round(x/a) (> 1)
!!$		       rather than ceil(x/a) as in the original
!!$		       paper, which should slightly improve performance
!!$     		       (and, apparently, slightly improves accuracy)
!!$      19 October 2012: Revised (SGJ) to fix bugs for large x, large -y,
!!$                       and 15<x<26. Performance improvements. Prototype
!!$		       now supplies default value for relerr.
!!$      24 October 2012: Switch to continued-fraction expansion for
!!$                       sufficiently large z, for performance reasons.
!!$		       Also, avoid spurious overflow for |z| > 1e154.
!!$		       Set relerr argument to min(relerr,0.1).
!!$      27 October 2012: Enhance accuracy in Re[w(z)] taken by itself,
!!$                       by switching to Alg. 916 in a region near
!!$		       the real-z axis where continued fractions
!!$		       have poor relative accuracy in Re[w(z)].  Thanks
!!$		       to M. Zaghloul for the tip.
!!$      29 October 2012: Replace SLATEC-derived erfcx routine with
!!$                       completely rewritten code by me, using a very
!!$                      different algorithm which is much faster.
!!$      30 October 2012: Implemented special-case code for real z
!!$                       (where real part is exp(-x^2) and imag part is
!!$                       Dawson integral), using algorithm similar to erfx.
!!$                      Export ImFaddeeva_w function to make Dawson's
!!$                      integral directly accessible.

module faddeeva
#ifdef USENRTYPE
  use nrtype
#endif
  implicit none
  private
  public :: faddeeva_w, erf
#ifndef USENRTYPE
  ! include these explicitly so we can compile as standalone if we want
  integer, parameter :: dp  = kind(1.d0)
  integer, parameter :: dpc = kind((1.d0,1.d0))
  real(dp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582d0
  complex(dpc), parameter :: IU = (0.d0,1.d0)
#endif
  real(dp), parameter :: ispi = 0.56418958354775628694807945156d0 ! 1 / sqrt(pi)

  interface erf
     module procedure cderf
  end interface erf
  private :: cderf
contains
  !-----------------------------------------------------------------------
  complex(dpc) pure elemental function cderf(z)
    complex(dpc), intent(in) :: z
    ! using recommendations from http://ab-initio.mit.edu/Faddeeva_w
    if (real(z)>=0) then
       cderf = 1.d0 - faddeeva_w( IU*z) * exp(-z**2)
    else
       cderf = faddeeva_w(-IU*z) * exp(-z**2) - 1.d0
    end if
  end function cderf
  !-----------------------------------------------------------------------
  ! return sinc(x) = sin(x)/x, given both x and sin(x)
  ! [since we only use this in cases where sin(x) has already been computed]
  real(dp) pure elemental function sinc(x,sinx)
    real(dp), intent(in) :: x, sinx
    if (abs(x)<1.d-4) then
       sinc = 1 - (0.1666666666666666666667d0)*x**2
    else
       sinc = sinx / x
    end if
  end function sinc

  ! sinh(x) via Taylor series, accurate to machine precision for |x| < 1e-2
  real(dp) pure elemental function sinh_taylor(x)
    real(dp), intent(in) :: x
    sinh_taylor = x * (1 + x**2 * (0.1666666666666666666667d0 &
         &                       + 0.00833333333333333333333d0 * x**2))
  end function sinh_taylor

  complex(dpc) pure elemental function polar(r,phi)
    real(dp), intent(in) :: r, phi
    polar = r * cmplx(cos(phi),sin(phi),dpc)
  end function polar

  complex(dpc) pure elemental function Faddeeva_w(z, relerr_in) result(ret)
    real(dp), parameter :: c0=3.9, c1=11.398, c2=0.08254, c3=0.1421, c4=0.2023 ! fit
    real(dp), parameter :: DBL_EPSILON = epsilon(1.d0)
    complex(dpc), intent(in) :: z
    real(dp), intent(in), optional :: relerr_in
    real(dp) :: relerr, a, a2, c
    real(dp) :: x, y, ya, x2, ax2, xs, xya, yax
    real(dp) :: sum1, sum2, sum3, sum4, sum5
    real(dp) :: coef, coef1, coef2, denom
    real(dp) :: prod2ax, prodm2ax
    real(dp) :: dr, di, dx, exp1, exp1dn, exp2ax, expm2ax
    real(dp) :: cos2xy, sin2xy, sinxy, expx2, expx2erfcxy
    real(dp) :: n0, nm, np, nu, tm, tp, wr, wi
    integer :: dn, n
    ! precomputed table of expa2n2[n-1] = exp(-a2*n*n) for double-precision a2 = 0.26865... below
    real(dp), parameter :: expa2n2(52) = (/ &
         7.64405281671221563d-01, 3.41424527166548425d-01, 8.91072646929412548d-02, &
         1.35887299055460086d-02, 1.21085455253437481d-03, 6.30452613933449404d-05, &
         1.91805156577114683d-06, 3.40969447714832381d-08, 3.54175089099469393d-10, &
         2.14965079583260682d-12, 7.62368911833724354d-15, 1.57982797110681093d-17, &
         1.91294189103582677d-20, 1.35344656764205340d-23, 5.59535712428588720d-27, &
         1.35164257972401769d-30, 1.90784582843501167d-34, 1.57351920291442930d-38, &
         7.58312432328032845d-43, 2.13536275438697082d-47, 3.51352063787195769d-52, &
         3.37800830266396920d-57, 1.89769439468301000d-62, 6.22929926072668851d-68, &
         1.19481172006938722d-73, 1.33908181133005953d-79, 8.76924303483223939d-86, &
         3.35555576166254986d-92, 7.50264110688173024d-99, 9.80192200745410268d-106, &
         7.48265412822268959d-113, 3.33770122566809425d-120, 8.69934598159861140d-128, &
         1.32486951484088852d-135, 1.17898144201315253d-143, 6.13039120236180012d-152, &
         1.86258785950822098d-160, 3.30668408201432783d-169, 3.43017280887946235d-178, &
         2.07915397775808219d-187, 7.36384545323984966d-197, 1.52394760394085741d-206, &
         1.84281935046532100d-216, 1.30209553802992923d-226, 5.37588903521080531d-237, &
         1.29689584599763145d-247, 1.82813078022866562d-258, 1.50576355348684241d-269, &
         7.24692320799294194d-281, 2.03797051314726829d-292, 3.34880215927873807d-304, &
         0.d0 /) ! underflow (also prevents reads past array end, below)

    if (real(z)==0.d0) then
       ret = cmplx(erfcx(aimag(z)),real(z),dpc) ! give correct sign of 0 in imag(w)
       return
    else if (aimag(z) == 0) then
       ret = cmplx(exp(-real(z)**2),ImFaddeeva_w(real(z)),dpc)
       return
    end if

    relerr = DBL_EPSILON
    if (present(relerr_in)) relerr = relerr_in

    if (relerr <= DBL_EPSILON) then
       relerr = DBL_EPSILON
       a = 0.518321480430085929872d0 ! pi / sqrt(-log(eps*0.5))
       c = 0.329973702884629072537d0 ! (2/pi) * a
       a2 = 0.268657157075235951582d0 ! a^2
    else
       if (relerr > 0.1d0) relerr = 0.1d0 ! not sensible to compute < 1 digit
       a = pi / sqrt(-log(relerr*0.5d0))
       c = (2/pi)*a
       a2 = a*a
    end if
    x = abs(real(z))
    y = aimag(z); ya = abs(y)

    ret = 0 ! return value

    sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0; sum5 = 0

    ! 1 to use continued fraction for large |z|
#define USE_CONTINUED_FRACTION 1

#if USE_CONTINUED_FRACTION
    ! continued fraction is faster
    ! As pointed out by M. Zaghloul, the continued fraction seems to give a large relative error in
    ! Re w(z) for |x| ~ 6 and small |y|, so use algorithm 816 in this region:
    if (ya > 7 .or. (x > 6 .and. (ya > 0.1d0 .or. (x > 8 .and. ya > 1.d-10) .or. x > 28))) then
       ! Poppe & Wijers suggest using a number of terms
       !    nu = 3 + 1442 / (26*rho + 77)
       ! where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
       ! (They only use this expansion for rho >= 1, but rho a little less
       ! than 1 seems okay too.)
       ! Instead, I did my own fit to a slightly different function
       ! that avoids the hypotenuse calculation, using NLopt to minimize
       ! the sum of the squares of the errors in nu with the constraint
       ! that the estimated nu be >= minimum nu to attain machine precision.
       ! I also separate the regions where nu == 2 and nu == 1.
       xs = real(z)
       if (y<0) xs = -xs ! compute for -z if y < 0
       if (x + ya > 4000) then ! nu <= 2
          if (x + ya > 1.d7) then ! nu == 1, w(z) = i/sqrt(pi) / z
             ! scale to avoid overflow
             if (x > ya) then
                yax = ya / xs
                denom = ispi / (xs + yax*ya)
                ret = cmplx(denom*yax, denom, dpc)
             else
                xya = xs / ya
                denom = ispi / (xya*xs + ya)
                ret = cmplx(denom, denom*xya, dpc)
             end if
          else ! nu == 2, w(z) = i/sqrt(pi) * z / (z*z - 0.5)
             dr = xs*xs - ya*ya - 0.5d0
             di = 2*xs*ya
             denom = ispi / (dr*dr + di*di)
             ret = cmplx(denom * (xs*di-ya*dr), denom * (xs*dr+ya*di), dpc)
          end if
       else ! compute nu(z) estimate and do general continued fraction
          nu = floor(c0 + c1 / (c2*x + c3*ya + c4))
          wr = xs
          wi = ya
          nu = 0.5d0 * (nu-1)
          do while (nu>0.4d0)
             ! w <- z - nu/w:
             denom = nu / (wr*wr + wi*wi)
             wr = xs - wr * denom
             wi = ya + wi * denom
             nu = nu - 0.5d0
          end do
          ! w(z) = i/sqrt(pi) / w:
          denom = ispi / (wr*wr + wi*wi)
          ret = cmplx(denom*wi, denom*wr, dpc)
       end if
       if (y < 0) then
          ! use w(z) = 2.0*exp(-z*z) - w(-z),
          ! but be careful of overflow in exp(-z*z)
          !                                = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
          ret = 2.d0*exp(cmplx((ya-xs)*(xs+ya), 2*xs*y, dpc)) - ret
       end if
       return
#else
       ! !USE_CONTINUED_FRACTION
       if (x + ya > 1.d7) then ! w(z) = i/sqrt(pi) / z, to machine precision
          xs = real(z)
          if (y<0) xs = -xs ! compute for -z if y < 0
          ! scale to avoid overflow
          if (x > ya) then
             yax = ya / xs
             denom = ispi / (xs + yax*ya)
             ret = cmplx(denom*yax, denom, dpc)
          else
             xya = xs / ya
             denom = ispi / (xya*xs + ya)
             ret = cmplx(denom, denom*xya, dpc)
          end if
          if (y < 0) then
             ! use w(z) = 2.0*exp(-z*z) - w(-z),
             ! but be careful of overflow in exp(-z*z)
             !                                = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
             ret = 2.d0*exp(cmplx((ya-xs)*(xs+ya), 2*xs*y, dpc)) - ret
          end if
          return
#endif
          ! USE_CONTINUED_FRACTION
       else if (x < 10) then
          ! Note: The test that seems to be suggested in the paper is x <
          ! sqrt(-log(DBL_MIN)), about 26.6, since otherwise exp(-x^2)
          ! underflows to zero and sum1,sum2,sum4 are zero.  However, long
          ! before this occurs, the sum1,sum2,sum4 contributions are
          ! negligible in double precision; I find that this happens for x >
          ! about 6, for all y.  On the other hand, I find that the case
          ! where we compute all of the sums is faster (at least with the
          ! precomputed expa2n2 table) until about x=10.  Furthermore, if we
          ! try to compute all of the sums for x > 20, I find that we
          ! sometimes run into numerical problems because underflow/overflow
          ! problems start to appear in the various coefficients of the sums,
          ! below.  Therefore, we use x < 10 here.
          prod2ax = 1
          prodm2ax = 1

          ! Somewhat ugly copy-and-paste duplication here, but I see significant
          ! speedups from using the special-case code with the precomputed
          ! exponential, and the x < 5e-4 special case is needed for accuracy.
          if (relerr == DBL_EPSILON) then ! use precomputed exp(-a2*(n*n)) table
             if (x < 5.d-4) then ! compute sum4 and sum5 together as sum5-sum4
                x2 = x**2
                expx2 = 1 - x2 * (1 - 0.5d0*x2) ! exp(-x**2) via Taylor
                ! compute exp(2*a*x) and exp(-2*a*x) via Taylor, to double precision
                ax2 = 1.036642960860171859744d0*x ! 2*a*x
                exp2ax  = 1 + ax2 * (1 + ax2 * (0.5 + 0.166666666666666666667d0*ax2))
                expm2ax = 1 - ax2 * (1 - ax2 * (0.5 - 0.166666666666666666667d0*ax2))
                do n = 1, size(expa2n2)
                   coef = expa2n2(n) * expx2 / (a2*n**2 + y**2)
                   prod2ax  = prod2ax  * exp2ax
                   prodm2ax = prodm2ax * expm2ax
                   sum1 = sum1 + coef
                   sum2 = sum2 + coef * prodm2ax
                   sum3 = sum3 + coef * prod2ax
                   ! really = sum5 - sum4
                   sum5 = sum5 + coef * (2*a) * n * sinh_taylor((2*a)*n*x)
                   ! test convergence via sum3
                   if (coef * prod2ax < relerr * sum3) exit
                end do
             else ! x > 5e-4, compute sum4 and sum5 separately
                expx2  = exp(-x**2)
                exp2ax = exp((2*a)*x)
                expm2ax = 1 / exp2ax
                do n = 1, size(expa2n2)
                   coef = expa2n2(n) * expx2 / (a2*n**2 + y**2)
                   prod2ax  = prod2ax  * exp2ax
                   prodm2ax = prodm2ax * expm2ax
                   sum1 = sum1 + coef
                   sum2 = sum2 + coef * prodm2ax
                   sum4 = sum4 + (coef * prodm2ax) * (a*n)
                   sum3 = sum3 + coef * prod2ax
                   sum5 = sum5 + (coef * prod2ax) * (a*n)
                   ! test convergence via sum5, since this sum has the slowest decay
                   if ((coef * prod2ax) * (a*n) < relerr * sum5) exit
                end do
             end if
          else ! relerr != DBL_EPSILON, compute exp(-a2*(n*n)) on the fly
             exp2ax = exp((2*a)*x)
             expm2ax = 1 / exp2ax
             if (x < 5e-4) then ! compute sum4 and sum5 together as sum5-sum4
                x2 = x**2
                expx2 = 1 - x2 * (1 - 0.5*x2) ! exp(-x**2) via Taylor
                n = 1
                do
                   coef = exp(-a2*n**2) * expx2 / (a2*n**2 + y**2)
                   prod2ax  = prod2ax  * exp2ax
                   prodm2ax = prodm2ax * expm2ax
                   sum1 = sum1 + coef
                   sum2 = sum2 + coef * prodm2ax
                   sum3 = sum3 + coef * prod2ax

                   ! really = sum5 - sum4
                   sum5 = sum5 + coef * (2*a) * n * sinh_taylor((2*a)*n*x)

                   ! test convergence via sum3
                   if (coef * prod2ax < relerr * sum3) exit
                   n = n+1
                end do
             else ! x > 5e-4, compute sum4 and sum5 separately
                expx2 = exp(-x**2)
                n = 1
                do
                   coef = exp(-a2*n**2) * expx2 / (a2*n**2 + y**2)
                   prod2ax  = prod2ax *  exp2ax
                   prodm2ax = prodm2ax * expm2ax
                   sum1 = sum1 + coef
                   sum2 = sum2 + coef * prodm2ax
                   sum4 = sum4 + (coef * prodm2ax) * (a*n)
                   sum3 = sum3 + coef * prod2ax
                   sum5 = sum5 + (coef * prod2ax) * (a*n)
                   ! test convergence via sum5, since this sum has the slowest decay
                   if ((coef * prod2ax) * (a*n) < relerr * sum5) exit
                   n = n+1
                end do
             end if
          end if
          if (y>-6) then ! avoid spurious overflow for large negative y
             expx2erfcxy = expx2*erfcx(y)
          else ! for y < -6, erfcx(y) = 2*exp(y**2) to double precision
             expx2erfcxy = 2*exp(y**2-x**2)
          end if
          if (y > 5) then ! imaginary terms cancel
             sinxy = sin(x*y)
             ret = (expx2erfcxy - c*y*sum1) * cos(2*x*y) &
                  + (c*x*expx2) * sinxy * sinc(x*y, sinxy)
          else
             xs = real(z)
             sinxy = sin(xs*y)
             sin2xy = sin(2*xs*y)
             cos2xy = cos(2*xs*y)
             coef1 = expx2erfcxy - c*y*sum1
             coef2 = c*xs*expx2
             ret = cmplx(coef1 * cos2xy + coef2 * sinxy * sinc(xs*y, sinxy), &
                  &        coef2 * sinc(2*xs*y, sin2xy) - coef1 * sin2xy, dpc)
          end if
       else ! x large: only sum3 & sum5 contribute (see above note)
#if USE_CONTINUED_FRACTION
          ret = exp(-x**2) ! |y| < 1e-10, so we only need exp(-x**2) term
#else
          if (y < 0) then
             ! erfcx(y) ~ 2*exp(y**2) + (< 1) if y < 0, so
             ! erfcx(y)*exp(-x**2) ~ 2*exp(y**2-x**2) term may not be negligible
             ! if y**2 - x**2 > -36 or so.  So, compute this term just in case.
             ! We also need the -exp(-x**2) term to compute Re[w] accurately
             ! in the case where y is very small.
             ret = polar(2*exp(y**2-x**2) - exp(-x**2), -2*real(z)*y)
          else
             ret = exp(-x**2) ! not negligible in real part if y very small
          end if
#endif
          ! (round instead of ceil as in original paper; note that x/a > 1 here)
          n0 = floor(x/a + 0.5) ! sum in both directions, starting at n0
          dx = a*n0 - x
          sum3 = exp(-dx**2) / (a2*n0**2 + y**2)
          sum5 = a*n0 * sum3
          exp1 = exp(4*a*dx)
          exp1dn = 1
          dn = 1
          do while (n0-dn>0)
             np = n0 + dn
             nm = n0 - dn
             tp = exp(-(a*dn+dx)**2)
             exp1dn = exp1dn * exp1
             tm = tp * exp1dn ! trick to get tm from tp
             tp = tp / (a2*np**2 + y**2)
             tm = tm / (a2*nm**2 + y**2)
             sum3 = sum3 + tp + tm
             sum5 = sum5 + a * (np * tp + nm * tm)
             if (a * (np * tp + nm * tm) < relerr * sum5) goto 100
             dn = dn + 1
          end do
          do ! loop over n0+dn terms only (since n0-dn <= 0)
             np = n0 + dn
             dn = dn + 1
             tp = exp(-(a*dn+dx)**2) / (a2*np**2 + y**2)
             sum3 = sum3 + tp
             sum5 = sum5 + a * np * tp
             if (a * np * tp < relerr * sum5) goto 100
          end do
       end if
100    ret = ret + cmplx((0.5*c)*y*(sum2+sum3), &
            &            (0.5*c)*sign(sum5-sum4, real(z)), dpc)
     end function Faddeeva_w
     !-----------------------------------------------------------------------
     !-----------------------------------------------------------------------
     ! erfcx(x) = exp(x^2) erfc(x) function, for real x, written by
     ! Steven G. Johnson, October 2012.
     !
     ! This function combines a few different ideas.
     !
     ! First, for x > 50, it uses a continued-fraction expansion (same as
     ! for the Faddeeva function, but with algebraic simplifications for z=i*x).
     !
     ! Second, for 0 <= x <= 50, it uses Chebyshev polynomial approximations,
     ! but with two twists:
     !
     !    a) It maps x to y = 4 / (4+x) in [0,1].  This simple transformation,
     !       inspired by a similar transformation in the octave-forge/specfun
     !    erfcx by Soren Hauberg, results in much faster Chebyshev convergence
     !    than other simple transformations I have examined.
     !
     !    b) Instead of using a single Chebyshev polynomial for the entire
     !       [0,1] y interval, we break the interval up into 100 equal
     !    subintervals, with a switch/lookup table, and use much lower
     !    degree Chebyshev polynomials in each subinterval. This greatly
     !    improves performance in my tests.
     !
     ! For x < 0, we use the relationship erfcx(-x) = 2 exp(x^2) - erfc(x),
     ! with the usual checks for overflow etcetera.
     !
     ! Performance-wise, it seems to be substantially faster than either
     ! the SLATEC DERFC function [or an erfcx function derived therefrom]
     ! or Cody's CALERF function (from netlib.org/specfun), while
     ! retaining near machine precision in accuracy.

     ! Given y100=100*y, where y = 4/(4+x) for x >= 0, compute erfc(x).
     ! Uses a look-up table of 100 different Chebyshev polynomials
     ! for y intervals [0,0.01], [0.01,0.02], ...., [0.99,1], generated
     ! with the help of Maple and a little shell script.   This allows
     ! the Chebyshev polynomials to be of significantly lower degree (about 1/4)
     ! compared to fitting the whole [0,1] interval with a single polynomial.
     real(dp) pure elemental function erfcx_y100(y100)
       real(dp), intent(in) :: y100
       real(dp) :: t
       select case(int(y100))
       case(0)
          t = 2*y100 - 1
          erfcx_y100 = 0.70878032454106438663d-3 + (0.71234091047026302958d-3 + (0.35779077297597742384d-5 +&
               & (0.17403143962587937815d-7 + (0.81710660047307788845d-10 + (0.36885022360434957634d-12 +&
               & 0.15917038551111111111d-14 * t) * t) * t) * t) * t) * t
       case(1)
          t = 2*y100 - 3
          erfcx_y100 = 0.21479143208285144230d-2 + (0.72686402367379996033d-3 + (0.36843175430938995552d-5 +&
               & (0.18071841272149201685d-7 + (0.85496449296040325555d-10 + (0.38852037518534291510d-12 +&
               & 0.16868473576888888889d-14 * t) * t) * t) * t) * t) * t
       case(2)
          t = 2*y100 - 5
          erfcx_y100 = 0.36165255935630175090d-2 + (0.74182092323555510862d-3 + (0.37948319957528242260d-5 +&
               & (0.18771627021793087350d-7 + (0.89484715122415089123d-10 + (0.40935858517772440862d-12 +&
               & 0.17872061464888888889d-14 * t) * t) * t) * t) * t) * t
       case(3)
          t = 2*y100 - 7
          erfcx_y100 = 0.51154983860031979264d-2 + (0.75722840734791660540d-3 + (0.39096425726735703941d-5 +&
               & (0.19504168704300468210d-7 + (0.93687503063178993915d-10 + (0.43143925959079664747d-12 +&
               & 0.18939926435555555556d-14 * t) * t) * t) * t) * t) * t
       case(4)
          t = 2*y100 - 9
          erfcx_y100 = 0.66457513172673049824d-2 + (0.77310406054447454920d-3 + (0.40289510589399439385d-5 +&
               & (0.20271233238288381092d-7 + (0.98117631321709100264d-10 + (0.45484207406017752971d-12 +&
               & 0.20076352213333333333d-14 * t) * t) * t) * t) * t) * t
       case(5)
          t = 2*y100 - 11
          erfcx_y100 = 0.82082389970241207883d-2 + (0.78946629611881710721d-3 + (0.41529701552622656574d-5 +&
               & (0.21074693344544655714d-7 + (0.10278874108587317989d-9 + (0.47965201390613339638d-12 +&
               & 0.21285907413333333333d-14 * t) * t) * t) * t) * t) * t
       case(6)
          t = 2*y100 - 13
          erfcx_y100 = 0.98039537275352193165d-2 + (0.80633440108342840956d-3 + (0.42819241329736982942d-5 +&
               & (0.21916534346907168612d-7 + (0.10771535136565470914d-9 + (0.50595972623692822410d-12 +&
               & 0.22573462684444444444d-14 * t) * t) * t) * t) * t) * t
       case(7)
          t = 2*y100 - 15
          erfcx_y100 = 0.11433927298290302370d-1 + (0.82372858383196561209d-3 + (0.44160495311765438816d-5 +&
               & (0.22798861426211986056d-7 + (0.11291291745879239736d-9 + (0.53386189365816880454d-12 +&
               & 0.23944209546666666667d-14 * t) * t) * t) * t) * t) * t
       case(8)
          t = 2*y100 - 17
          erfcx_y100 = 0.13099232878814653979d-1 + (0.84167002467906968214d-3 + (0.45555958988457506002d-5 +&
               & (0.23723907357214175198d-7 + (0.11839789326602695603d-9 + (0.56346163067550237877d-12 +&
               & 0.25403679644444444444d-14 * t) * t) * t) * t) * t) * t
       case(9)
          t = 2*y100 - 19
          erfcx_y100 = 0.14800987015587535621d-1 + (0.86018092946345943214d-3 + (0.47008265848816866105d-5 +&
               & (0.24694040760197315333d-7 + (0.12418779768752299093d-9 + (0.59486890370320261949d-12 +&
               & 0.26957764568888888889d-14 * t) * t) * t) * t) * t) * t
       case(10)
          t = 2*y100 - 21
          erfcx_y100 = 0.16540351739394069380d-1 + (0.87928458641241463952d-3 + (0.48520195793001753903d-5 +&
               & (0.25711774900881709176d-7 + (0.13030128534230822419d-9 + (0.62820097586874779402d-12 +&
               & 0.28612737351111111111d-14 * t) * t) * t) * t) * t) * t
       case(11)
          t = 2*y100 - 23
          erfcx_y100 = 0.18318536789842392647d-1 + (0.89900542647891721692d-3 + (0.50094684089553365810d-5 +&
               & (0.26779777074218070482d-7 + (0.13675822186304615566d-9 + (0.66358287745352705725d-12 +&
               & 0.30375273884444444444d-14 * t) * t) * t) * t) * t) * t
       case(12)
          t = 2*y100 - 25
          erfcx_y100 = 0.20136801964214276775d-1 + (0.91936908737673676012d-3 + (0.51734830914104276820d-5 +&
               & (0.27900878609710432673d-7 + (0.14357976402809042257d-9 + (0.70114790311043728387d-12 +&
               & 0.32252476000000000000d-14 * t) * t) * t) * t) * t) * t
       case(13)
          t = 2*y100 - 27
          erfcx_y100 = 0.21996459598282740954d-1 + (0.94040248155366777784d-3 + (0.53443911508041164739d-5 +&
               & (0.29078085538049374673d-7 + (0.15078844500329731137d-9 + (0.74103813647499204269d-12 +&
               & 0.34251892320000000000d-14 * t) * t) * t) * t) * t) * t
       case(14)
          t = 2*y100 - 29
          erfcx_y100 = 0.23898877187226319502d-1 + (0.96213386835900177540d-3 + (0.55225386998049012752d-5 +&
               & (0.30314589961047687059d-7 + (0.15840826497296335264d-9 + (0.78340500472414454395d-12 +&
               & 0.36381553564444444445d-14 * t) * t) * t) * t) * t) * t
       case(15)
          t = 2*y100 - 31
          erfcx_y100 = 0.25845480155298518485d-1 + (0.98459293067820123389d-3 + (0.57082915920051843672d-5 +&
               & (0.31613782169164830118d-7 + (0.16646478745529630813d-9 + (0.82840985928785407942d-12 +&
               & 0.38649975768888888890d-14 * t) * t) * t) * t) * t) * t
       case(16)
          t = 2*y100 - 33
          erfcx_y100 = 0.27837754783474696598d-1 + (0.10078108563256892757d-2 + (0.59020366493792212221d-5 +&
               & (0.32979263553246520417d-7 + (0.17498524159268458073d-9 + (0.87622459124842525110d-12 +&
               & 0.41066206488888888890d-14 * t) * t) * t) * t) * t) * t
       case(17)
          t = 2*y100 - 35
          erfcx_y100 = 0.29877251304899307550d-1 + (0.10318204245057349310d-2 + (0.61041829697162055093d-5 +&
               & (0.34414860359542720579d-7 + (0.18399863072934089607d-9 + (0.92703227366365046533d-12 +&
               & 0.43639844053333333334d-14 * t) * t) * t) * t) * t) * t
       case(18)
          t = 2*y100 - 37
          erfcx_y100 = 0.31965587178596443475d-1 + (0.10566560976716574401d-2 + (0.63151633192414586770d-5 +&
               & (0.35924638339521924242d-7 + (0.19353584758781174038d-9 + (0.98102783859889264382d-12 +&
               & 0.46381060817777777779d-14 * t) * t) * t) * t) * t) * t
       case(19)
          t = 2*y100 - 39
          erfcx_y100 = 0.34104450552588334840d-1 + (0.10823541191350532574d-2 + (0.65354356159553934436d-5 +&
               & (0.37512918348533521149d-7 + (0.20362979635817883229d-9 + (0.10384187833037282363d-11 +&
               & 0.49300625262222222221d-14 * t) * t) * t) * t) * t) * t
       case(20)
          t = 2*y100 - 41
          erfcx_y100 = 0.36295603928292425716d-1 + (0.11089526167995268200d-2 + (0.67654845095518363577d-5 +&
               & (0.39184292949913591646d-7 + (0.21431552202133775150d-9 + (0.10994259106646731797d-11 +&
               & 0.52409949102222222221d-14 * t) * t) * t) * t) * t) * t
       case(21)
          t = 2*y100 - 43
          erfcx_y100 = 0.38540888038840509795d-1 + (0.11364917134175420009d-2 + (0.70058230641246312003d-5 +&
               & (0.40943644083718586939d-7 + (0.22563034723692881631d-9 + (0.11642841011361992885d-11 +&
               & 0.55721092871111111110d-14 * t) * t) * t) * t) * t) * t
       case(22)
          t = 2*y100 - 45
          erfcx_y100 = 0.40842225954785960651d-1 + (0.11650136437945673891d-2 + (0.72569945502343006619d-5 +&
               & (0.42796161861855042273d-7 + (0.23761401711005024162d-9 + (0.12332431172381557035d-11 +&
               & 0.59246802364444444445d-14 * t) * t) * t) * t) * t) * t
       case(23)
          t = 2*y100 - 47
          erfcx_y100 = 0.43201627431540222422d-1 + (0.11945628793917272199d-2 + (0.75195743532849206263d-5 +&
               & (0.44747364553960993492d-7 + (0.25030885216472953674d-9 + (0.13065684400300476484d-11 +&
               & 0.63000532853333333334d-14 * t) * t) * t) * t) * t) * t
       case(24)
          t = 2*y100 - 49
          erfcx_y100 = 0.45621193513810471438d-1 + (0.12251862608067529503d-2 + (0.77941720055551920319d-5 +&
               & (0.46803119830954460212d-7 + (0.26375990983978426273d-9 + (0.13845421370977119765d-11 +&
               & 0.66996477404444444445d-14 * t) * t) * t) * t) * t) * t
       case(25)
          t = 2*y100 - 51
          erfcx_y100 = 0.48103121413299865517d-1 + (0.12569331386432195113d-2 + (0.80814333496367673980d-5 +&
               & (0.48969667335682018324d-7 + (0.27801515481905748484d-9 + (0.14674637611609884208d-11 +&
               & 0.71249589351111111110d-14 * t) * t) * t) * t) * t) * t
       case(26)
          t = 2*y100 - 53
          erfcx_y100 = 0.50649709676983338501d-1 + (0.12898555233099055810d-2 + (0.83820428414568799654d-5 +&
               & (0.51253642652551838659d-7 + (0.29312563849675507232d-9 + (0.15556512782814827846d-11 +&
               & 0.75775607822222222221d-14 * t) * t) * t) * t) * t) * t
       case(27)
          t = 2*y100 - 55
          erfcx_y100 = 0.53263363664388864181d-1 + (0.13240082443256975769d-2 + (0.86967260015007658418d-5 +&
               & (0.53662102750396795566d-7 + (0.30914568786634796807d-9 + (0.16494420240828493176d-11 +&
               & 0.80591079644444444445d-14 * t) * t) * t) * t) * t) * t
       case(28)
          t = 2*y100 - 57
          erfcx_y100 = 0.55946601353500013794d-1 + (0.13594491197408190706d-2 + (0.90262520233016380987d-5 +&
               & (0.56202552975056695376d-7 + (0.32613310410503135996d-9 + (0.17491936862246367398d-11 +&
               & 0.85713381688888888890d-14 * t) * t) * t) * t) * t) * t
       case(29)
          t = 2*y100 - 59
          erfcx_y100 = 0.58702059496154081813d-1 + (0.13962391363223647892d-2 + (0.93714365487312784270d-5 +&
               & (0.58882975670265286526d-7 + (0.34414937110591753387d-9 + (0.18552853109751857859d-11 +&
               & 0.91160736711111111110d-14 * t) * t) * t) * t) * t) * t
       case(30)
          t = 2*y100 - 61
          erfcx_y100 = 0.61532500145144778048d-1 + (0.14344426411912015247d-2 + (0.97331446201016809696d-5 +&
               & (0.61711860507347175097d-7 + (0.36325987418295300221d-9 + (0.19681183310134518232d-11 +&
               & 0.96952238400000000000d-14 * t) * t) * t) * t) * t) * t
       case(31)
          t = 2*y100 - 63
          erfcx_y100 = 0.64440817576653297993d-1 + (0.14741275456383131151d-2 + (0.10112293819576437838d-4 +&
               & (0.64698236605933246196d-7 + (0.38353412915303665586d-9 + (0.20881176114385120186d-11 +&
               & 0.10310784480000000000d-13 * t) * t) * t) * t) * t) * t
       case(32)
          t = 2*y100 - 65
          erfcx_y100 = 0.67430045633130393282d-1 + (0.15153655418916540370d-2 + (0.10509857606888328667d-4 +&
               & (0.67851706529363332855d-7 + (0.40504602194811140006d-9 + (0.22157325110542534469d-11 +&
               & 0.10964842115555555556d-13 * t) * t) * t) * t) * t) * t
       case(33)
          t = 2*y100 - 67
          erfcx_y100 = 0.70503365513338850709d-1 + (0.15582323336495709827d-2 + (0.10926868866865231089d-4 +&
               & (0.71182482239613507542d-7 + (0.42787405890153386710d-9 + (0.23514379522274416437d-11 +&
               & 0.11659571751111111111d-13 * t) * t) * t) * t) * t) * t
       case(34)
          t = 2*y100 - 69
          erfcx_y100 = 0.73664114037944596353d-1 + (0.16028078812438820413d-2 + (0.11364423678778207991d-4 +&
               & (0.74701423097423182009d-7 + (0.45210162777476488324d-9 + (0.24957355004088569134d-11 +&
               & 0.12397238257777777778d-13 * t) * t) * t) * t) * t) * t
       case(35)
          t = 2*y100 - 71
          erfcx_y100 = 0.76915792420819562379d-1 + (0.16491766623447889354d-2 + (0.11823685320041302169d-4 +&
               & (0.78420075993781544386d-7 + (0.47781726956916478925d-9 + (0.26491544403815724749d-11 +&
               & 0.13180196462222222222d-13 * t) * t) * t) * t) * t) * t
       case(36)
          t = 2*y100 - 73
          erfcx_y100 = 0.80262075578094612819d-1 + (0.16974279491709504117d-2 + (0.12305888517309891674d-4 +&
               & (0.82350717698979042290d-7 + (0.50511496109857113929d-9 + (0.28122528497626897696d-11 +&
               & 0.14010889635555555556d-13 * t) * t) * t) * t) * t) * t
       case(37)
          t = 2*y100 - 75
          erfcx_y100 = 0.83706822008980357446d-1 + (0.17476561032212656962d-2 + (0.12812343958540763368d-4 +&
               & (0.86506399515036435592d-7 + (0.53409440823869467453d-9 + (0.29856186620887555043d-11 +&
               & 0.14891851591111111111d-13 * t) * t) * t) * t) * t) * t
       case(38)
          t = 2*y100 - 77
          erfcx_y100 = 0.87254084284461718231d-1 + (0.17999608886001962327d-2 + (0.13344443080089492218d-4 +&
               & (0.90900994316429008631d-7 + (0.56486134972616465316d-9 + (0.31698707080033956934d-11 +&
               & 0.15825697795555555556d-13 * t) * t) * t) * t) * t) * t
       case(39)
          t = 2*y100 - 79
          erfcx_y100 = 0.90908120182172748487d-1 + (0.18544478050657699758d-2 + (0.13903663143426120077d-4 +&
               & (0.95549246062549906177d-7 + (0.59752787125242054315d-9 + (0.33656597366099099413d-11 +&
               & 0.16815130613333333333d-13 * t) * t) * t) * t) * t) * t
       case(40)
          t = 2*y100 - 81
          erfcx_y100 = 0.94673404508075481121d-1 + (0.19112284419887303347d-2 + (0.14491572616545004930d-4 +&
               & (0.10046682186333613697d-6 + (0.63221272959791000515d-9 + (0.35736693975589130818d-11 +&
               & 0.17862931591111111111d-13 * t) * t) * t) * t) * t) * t
       case(41)
          t = 2*y100 - 83
          erfcx_y100 = 0.98554641648004456555d-1 + (0.19704208544725622126d-2 + (0.15109836875625443935d-4 +&
               & (0.10567036667675984067d-6 + (0.66904168640019354565d-9 + (0.37946171850824333014d-11 +&
               & 0.18971959040000000000d-13 * t) * t) * t) * t) * t) * t
       case(42)
          t = 2*y100 - 85
          erfcx_y100 = 0.10255677889470089531d0 + (0.20321499629472857418d-2 + (0.15760224242962179564d-4 +&
               & (0.11117756071353507391d-6 + (0.70814785110097658502d-9 + (0.40292553276632563925d-11 +&
               & 0.20145143075555555556d-13 * t) * t) * t) * t) * t) * t
       case(43)
          t = 2*y100 - 87
          erfcx_y100 = 0.10668502059865093318d0 + (0.20965479776148731610d-2 + (0.16444612377624983565d-4 +&
               & (0.11700717962026152749d-6 + (0.74967203250938418991d-9 + (0.42783716186085922176d-11 +&
               & 0.21385479360000000000d-13 * t) * t) * t) * t) * t) * t
       case(44)
          t = 2*y100 - 89
          erfcx_y100 = 0.11094484319386444474d0 + (0.21637548491908170841d-2 + (0.17164995035719657111d-4 +&
               & (0.12317915750735938089d-6 + (0.79376309831499633734d-9 + (0.45427901763106353914d-11 +&
               & 0.22696025653333333333d-13 * t) * t) * t) * t) * t) * t
       case(45)
          t = 2*y100 - 91
          erfcx_y100 = 0.11534201115268804714d0 + (0.22339187474546420375d-2 + (0.17923489217504226813d-4 +&
               & (0.12971465288245997681d-6 + (0.84057834180389073587d-9 + (0.48233721206418027227d-11 +&
               & 0.24079890062222222222d-13 * t) * t) * t) * t) * t) * t
       case(46)
          t = 2*y100 - 93
          erfcx_y100 = 0.11988259392684094740d0 + (0.23071965691918689601d-2 + (0.18722342718958935446d-4 +&
               & (0.13663611754337957520d-6 + (0.89028385488493287005d-9 + (0.51210161569225846701d-11 +&
               & 0.25540227111111111111d-13 * t) * t) * t) * t) * t) * t
       case(47)
          t = 2*y100 - 95
          erfcx_y100 = 0.12457298393509812907d0 + (0.23837544771809575380d-2 + (0.19563942105711612475d-4 +&
               & (0.14396736847739470782d-6 + (0.94305490646459247016d-9 + (0.54366590583134218096d-11 +&
               & 0.27080225920000000000d-13 * t) * t) * t) * t) * t) * t
       case(48)
          t = 2*y100 - 97
          erfcx_y100 = 0.12941991566142438816d0 + (0.24637684719508859484d-2 + (0.20450821127475879816d-4 +&
               & (0.15173366280523906622d-6 + (0.99907632506389027739d-9 + (0.57712760311351625221d-11 +&
               & 0.28703099555555555556d-13 * t) * t) * t) * t) * t) * t
       case(49)
          t = 2*y100 - 99
          erfcx_y100 = 0.13443048593088696613d0 + (0.25474249981080823877d-2 + (0.21385669591362915223d-4 +&
               & (0.15996177579900443030d-6 + (0.10585428844575134013d-8 + (0.61258809536787882989d-11 +&
               & 0.30412080142222222222d-13 * t) * t) * t) * t) * t) * t
       case(50)
          t = 2*y100 - 101
          erfcx_y100 = 0.13961217543434561353d0 + (0.26349215871051761416d-2 + (0.22371342712572567744d-4 +&
               & (0.16868008199296822247d-6 + (0.11216596910444996246d-8 + (0.65015264753090890662d-11 +&
               & 0.32210394506666666666d-13 * t) * t) * t) * t) * t) * t
       case(51)
          t = 2*y100 - 103
          erfcx_y100 = 0.14497287157673800690d0 + (0.27264675383982439814d-2 + (0.23410870961050950197d-4 +&
               & (0.17791863939526376477d-6 + (0.11886425714330958106d-8 + (0.68993039665054288034d-11 +&
               & 0.34101266222222222221d-13 * t) * t) * t) * t) * t) * t
       case(52)
          t = 2*y100 - 105
          erfcx_y100 = 0.15052089272774618151d0 + (0.28222846410136238008d-2 + (0.24507470422713397006d-4 +&
               & (0.18770927679626136909d-6 + (0.12597184587583370712d-8 + (0.73203433049229821618d-11 +&
               & 0.36087889048888888890d-13 * t) * t) * t) * t) * t) * t
       case(53)
          t = 2*y100 - 107
          erfcx_y100 = 0.15626501395774612325d0 + (0.29226079376196624949d-2 + (0.25664553693768450545d-4 +&
               & (0.19808568415654461964d-6 + (0.13351257759815557897d-8 + (0.77658124891046760667d-11 +&
               & 0.38173420035555555555d-13 * t) * t) * t) * t) * t) * t
       case(54)
          t = 2*y100 - 109
          erfcx_y100 = 0.16221449434620737567d0 + (0.30276865332726475672d-2 + (0.26885741326534564336d-4 +&
               & (0.20908350604346384143d-6 + (0.14151148144240728728d-8 + (0.82369170665974313027d-11 +&
               & 0.40360957457777777779d-13 * t) * t) * t) * t) * t) * t
       case(55)
          t = 2*y100 - 111
          erfcx_y100 = 0.16837910595412130659d0 + (0.31377844510793082301d-2 + (0.28174873844911175026d-4 +&
               & (0.22074043807045782387d-6 + (0.14999481055996090039d-8 + (0.87348993661930809254d-11 +&
               & 0.42653528977777777779d-13 * t) * t) * t) * t) * t) * t
       case(56)
          t = 2*y100 - 113
          erfcx_y100 = 0.17476916455659369953d0 + (0.32531815370903068316d-2 + (0.29536024347344364074d-4 +&
               & (0.23309632627767074202d-6 + (0.15899007843582444846d-8 + (0.92610375235427359475d-11 +&
               & 0.45054073102222222221d-13 * t) * t) * t) * t) * t) * t
       case(57)
          t = 2*y100 - 115
          erfcx_y100 = 0.18139556223643701364d0 + (0.33741744168096996041d-2 + (0.30973511714709500836d-4 +&
               & (0.24619326937592290996d-6 + (0.16852609412267750744d-8 + (0.98166442942854895573d-11 +&
               & 0.47565418097777777779d-13 * t) * t) * t) * t) * t) * t
       case(58)
          t = 2*y100 - 117
          erfcx_y100 = 0.18826980194443664549d0 + (0.35010775057740317997d-2 + (0.32491914440014267480d-4 +&
               & (0.26007572375886319028d-6 + (0.17863299617388376116d-8 + (0.10403065638343878679d-10 +&
               & 0.50190265831111111110d-13 * t) * t) * t) * t) * t) * t
       case(59)
          t = 2*y100 - 119
          erfcx_y100 = 0.19540403413693967350d0 + (0.36342240767211326315d-2 + (0.34096085096200907289d-4 +&
               & (0.27479061117017637474d-6 + (0.18934228504790032826d-8 + (0.11021679075323598664d-10 +&
               & 0.52931171733333333334d-13 * t) * t) * t) * t) * t) * t
       case(60)
          t = 2*y100 - 121
          erfcx_y100 = 0.20281109560651886959d0 + (0.37739673859323597060d-2 + (0.35791165457592409054d-4 +&
               & (0.29038742889416172404d-6 + (0.20068685374849001770d-8 + (0.11673891799578381999d-10 +&
               & 0.55790523093333333334d-13 * t) * t) * t) * t) * t) * t
       case(61)
          t = 2*y100 - 123
          erfcx_y100 = 0.21050455062669334978d0 + (0.39206818613925652425d-2 + (0.37582602289680101704d-4 +&
               & (0.30691836231886877385d-6 + (0.21270101645763677824d-8 + (0.12361138551062899455d-10 +&
               & 0.58770520160000000000d-13 * t) * t) * t) * t) * t) * t
       case(62)
          t = 2*y100 - 125
          erfcx_y100 = 0.21849873453703332479d0 + (0.40747643554689586041d-2 + (0.39476163820986711501d-4 +&
               & (0.32443839970139918836d-6 + (0.22542053491518680200d-8 + (0.13084879235290858490d-10 +&
               & 0.61873153262222222221d-13 * t) * t) * t) * t) * t) * t
       case(63)
          t = 2*y100 - 127
          erfcx_y100 = 0.22680879990043229327d0 + (0.42366354648628516935d-2 + (0.41477956909656896779d-4 +&
               & (0.34300544894502810002d-6 + (0.23888264229264067658d-8 + (0.13846596292818514601d-10 +&
               & 0.65100183751111111110d-13 * t) * t) * t) * t) * t) * t
       case(64)
          t = 2*y100 - 129
          erfcx_y100 = 0.23545076536988703937d0 + (0.44067409206365170888d-2 + (0.43594444916224700881d-4 +&
               & (0.36268045617760415178d-6 + (0.25312606430853202748d-8 + (0.14647791812837903061d-10 +&
               & 0.68453122631111111110d-13 * t) * t) * t) * t) * t) * t
       case(65)
          t = 2*y100 - 131
          erfcx_y100 = 0.24444156740777432838d0 + (0.45855530511605787178d-2 + (0.45832466292683085475d-4 +&
               & (0.38352752590033030472d-6 + (0.26819103733055603460d-8 + (0.15489984390884756993d-10 +&
               & 0.71933206364444444445d-13 * t) * t) * t) * t) * t) * t
       case(66)
          t = 2*y100 - 133
          erfcx_y100 = 0.25379911500634264643d0 + (0.47735723208650032167d-2 + (0.48199253896534185372d-4 +&
               & (0.40561404245564732314d-6 + (0.28411932320871165585d-8 + (0.16374705736458320149d-10 +&
               & 0.75541379822222222221d-13 * t) * t) * t) * t) * t) * t
       case(67)
          t = 2*y100 - 135
          erfcx_y100 = 0.26354234756393613032d0 + (0.49713289477083781266d-2 + (0.50702455036930367504d-4 +&
               & (0.42901079254268185722d-6 + (0.30095422058900481753d-8 + (0.17303497025347342498d-10 +&
               & 0.79278273368888888890d-13 * t) * t) * t) * t) * t) * t
       case(68)
          t = 2*y100 - 137
          erfcx_y100 = 0.27369129607732343398d0 + (0.51793846023052643767d-2 + (0.53350152258326602629d-4 +&
               & (0.45379208848865015485d-6 + (0.31874057245814381257d-8 + (0.18277905010245111046d-10 +&
               & 0.83144182364444444445d-13 * t) * t) * t) * t) * t) * t
       case(69)
          t = 2*y100 - 139
          erfcx_y100 = 0.28426714781640316172d0 + (0.53983341916695141966d-2 + (0.56150884865255810638d-4 +&
               & (0.48003589196494734238d-6 + (0.33752476967570796349d-8 + (0.19299477888083469086d-10 +&
               & 0.87139049137777777779d-13 * t) * t) * t) * t) * t) * t
       case(70)
          t = 2*y100 - 141
          erfcx_y100 = 0.29529231465348519920d0 + (0.56288077305420795663d-2 + (0.59113671189913307427d-4 +&
               & (0.50782393781744840482d-6 + (0.35735475025851713168d-8 + (0.20369760937017070382d-10 +&
               & 0.91262442613333333334d-13 * t) * t) * t) * t) * t) * t
       case(71)
          t = 2*y100 - 143
          erfcx_y100 = 0.30679050522528838613d0 + (0.58714723032745403331d-2 + (0.62248031602197686791d-4 +&
               & (0.53724185766200945789d-6 + (0.37827999418960232678d-8 + (0.21490291930444538307d-10 +&
               & 0.95513539182222222221d-13 * t) * t) * t) * t) * t) * t
       case(72)
          t = 2*y100 - 145
          erfcx_y100 = 0.31878680111173319425d0 + (0.61270341192339103514d-2 + (0.65564012259707640976d-4 +&
               & (0.56837930287837738996d-6 + (0.40035151353392378882d-8 + (0.22662596341239294792d-10 +&
               & 0.99891109760000000000d-13 * t) * t) * t) * t) * t) * t
       case(73)
          t = 2*y100 - 147
          erfcx_y100 = 0.33130773722152622027d0 + (0.63962406646798080903d-2 + (0.69072209592942396666d-4 +&
               & (0.60133006661885941812d-6 + (0.42362183765883466691d-8 + (0.23888182347073698382d-10 +&
               & 0.10439349811555555556d-12 * t) * t) * t) * t) * t) * t
       case(74)
          t = 2*y100 - 149
          erfcx_y100 = 0.34438138658041336523d0 + (0.66798829540414007258d-2 + (0.72783795518603561144d-4 +&
               & (0.63619220443228800680d-6 + (0.44814499336514453364d-8 + (0.25168535651285475274d-10 +&
               & 0.10901861383111111111d-12 * t) * t) * t) * t) * t) * t
       case(75)
          t = 2*y100 - 151
          erfcx_y100 = 0.35803744972380175583d0 + (0.69787978834882685031d-2 + (0.76710543371454822497d-4 +&
               & (0.67306815308917386747d-6 + (0.47397647975845228205d-8 + (0.26505114141143050509d-10 +&
               & 0.11376390933333333333d-12 * t) * t) * t) * t) * t) * t
       case(76)
          t = 2*y100 - 153
          erfcx_y100 = 0.37230734890119724188d0 + (0.72938706896461381003d-2 + (0.80864854542670714092d-4 +&
               & (0.71206484718062688779d-6 + (0.50117323769745883805d-8 + (0.27899342394100074165d-10 +&
               & 0.11862637614222222222d-12 * t) * t) * t) * t) * t) * t
       case(77)
          t = 2*y100 - 155
          erfcx_y100 = 0.38722432730555448223d0 + (0.76260375162549802745d-2 + (0.85259785810004603848d-4 +&
               & (0.75329383305171327677d-6 + (0.52979361368388119355d-8 + (0.29352606054164086709d-10 +&
               & 0.12360253370666666667d-12 * t) * t) * t) * t) * t) * t
       case(78)
          t = 2*y100 - 157
          erfcx_y100 = 0.40282355354616940667d0 + (0.79762880915029728079d-2 + (0.89909077342438246452d-4 +&
               & (0.79687137961956194579d-6 + (0.55989731807360403195d-8 + (0.30866246101464869050d-10 +&
               & 0.12868841946666666667d-12 * t) * t) * t) * t) * t) * t
       case(79)
          t = 2*y100 - 159
          erfcx_y100 = 0.41914223158913787649d0 + (0.83456685186950463538d-2 + (0.94827181359250161335d-4 +&
               & (0.84291858561783141014d-6 + (0.59154537751083485684d-8 + (0.32441553034347469291d-10 +&
               & 0.13387957943111111111d-12 * t) * t) * t) * t) * t) * t
       case(80)
          t = 2*y100 - 161
          erfcx_y100 = 0.43621971639463786896d0 + (0.87352841828289495773d-2 + (0.10002929142066799966d-3 +&
               & (0.89156148280219880024d-6 + (0.62480008150788597147d-8 + (0.34079760983458878910d-10 +&
               & 0.13917107176888888889d-12 * t) * t) * t) * t) * t) * t
       case(81)
          t = 2*y100 - 163
          erfcx_y100 = 0.45409763548534330981d0 + (0.91463027755548240654d-2 + (0.10553137232446167258d-3 +&
               & (0.94293113464638623798d-6 + (0.65972492312219959885d-8 + (0.35782041795476563662d-10 +&
               & 0.14455745872000000000d-12 * t) * t) * t) * t) * t) * t
       case(82)
          t = 2*y100 - 165
          erfcx_y100 = 0.47282001668512331468d0 + (0.95799574408860463394d-2 + (0.11135019058000067469d-3 +&
               & (0.99716373005509038080d-6 + (0.69638453369956970347d-8 + (0.37549499088161345850d-10 +&
               & 0.15003280712888888889d-12 * t) * t) * t) * t) * t) * t
       case(83)
          t = 2*y100 - 167
          erfcx_y100 = 0.49243342227179841649d0 + (0.10037550043909497071d-1 + (0.11750334542845234952d-3 +&
               & (0.10544006716188967172d-5 + (0.73484461168242224872d-8 + (0.39383162326435752965d-10 +&
               & 0.15559069118222222222d-12 * t) * t) * t) * t) * t) * t
       case(84)
          t = 2*y100 - 169
          erfcx_y100 = 0.51298708979209258326d0 + (0.10520454564612427224d-1 + (0.12400930037494996655d-3 +&
               & (0.11147886579371265246d-5 + (0.77517184550568711454d-8 + (0.41283980931872622611d-10 +&
               & 0.16122419680000000000d-12 * t) * t) * t) * t) * t) * t
       case(85)
          t = 2*y100 - 171
          erfcx_y100 = 0.53453307979101369843d0 + (0.11030120618800726938d-1 + (0.13088741519572269581d-3 +&
               & (0.11784797595374515432d-5 + (0.81743383063044825400d-8 + (0.43252818449517081051d-10 +&
               & 0.16692592640000000000d-12 * t) * t) * t) * t) * t) * t
       case(86)
          t = 2*y100 - 173
          erfcx_y100 = 0.55712643071169299478d0 + (0.11568077107929735233d-1 + (0.13815797838036651289d-3 +&
               & (0.12456314879260904558d-5 + (0.86169898078969313597d-8 + (0.45290446811539652525d-10 +&
               & 0.17268801084444444444d-12 * t) * t) * t) * t) * t) * t
       case(87)
          t = 2*y100 - 175
          erfcx_y100 = 0.58082532122519320968d0 + (0.12135935999503877077d-1 + (0.14584223996665838559d-3 +&
               & (0.13164068573095710742d-5 + (0.90803643355106020163d-8 + (0.47397540713124619155d-10 +&
               & 0.17850211608888888889d-12 * t) * t) * t) * t) * t) * t
       case(88)
          t = 2*y100 - 177
          erfcx_y100 = 0.60569124025293375554d0 + (0.12735396239525550361d-1 + (0.15396244472258863344d-3 +&
               & (0.13909744385382818253d-5 + (0.95651595032306228245d-8 + (0.49574672127669041550d-10 +&
               & 0.18435945564444444444d-12 * t) * t) * t) * t) * t) * t
       case(89)
          t = 2*y100 - 179
          erfcx_y100 = 0.63178916494715716894d0 + (0.13368247798287030927d-1 + (0.16254186562762076141d-3 +&
               & (0.14695084048334056083d-5 + (0.10072078109604152350d-7 + (0.51822304995680707483d-10 +&
               & 0.19025081422222222222d-12 * t) * t) * t) * t) * t) * t
       case(90)
          t = 2*y100 - 181
          erfcx_y100 = 0.65918774689725319200d0 + (0.14036375850601992063d-1 + (0.17160483760259706354d-3 +&
               & (0.15521885688723188371d-5 + (0.10601827031535280590d-7 + (0.54140790105837520499d-10 +&
               & 0.19616655146666666667d-12 * t) * t) * t) * t) * t) * t
       case(91)
          t = 2*y100 - 183
          erfcx_y100 = 0.68795950683174433822d0 + (0.14741765091365869084d-1 + (0.18117679143520433835d-3 +&
               & (0.16392004108230585213d-5 + (0.11155116068018043001d-7 + (0.56530360194925690374d-10 +&
               & 0.20209663662222222222d-12 * t) * t) * t) * t) * t) * t
       case(92)
          t = 2*y100 - 185
          erfcx_y100 = 0.71818103808729967036d0 + (0.15486504187117112279d-1 + (0.19128428784550923217d-3 +&
               & (0.17307350969359975848d-5 + (0.11732656736113607751d-7 + (0.58991125287563833603d-10 +&
               & 0.20803065333333333333d-12 * t) * t) * t) * t) * t) * t
       case(93)
          t = 2*y100 - 187
          erfcx_y100 = 0.74993321911726254661d0 + (0.16272790364044783382d-1 + (0.20195505163377912645d-3 +&
               & (0.18269894883203346953d-5 + (0.12335161021630225535d-7 + (0.61523068312169087227d-10 +&
               & 0.21395783431111111111d-12 * t) * t) * t) * t) * t) * t
       case(94)
          t = 2*y100 - 189
          erfcx_y100 = 0.78330143531283492729d0 + (0.17102934132652429240d-1 + (0.21321800585063327041d-3 +&
               & (0.19281661395543913713d-5 + (0.12963340087354341574d-7 + (0.64126040998066348872d-10 +&
               & 0.21986708942222222222d-12 * t) * t) * t) * t) * t) * t
       case(95)
          t = 2*y100 - 191
          erfcx_y100 = 0.81837581041023811832d0 + (0.17979364149044223802d-1 + (0.22510330592753129006d-3 +&
               & (0.20344732868018175389d-5 + (0.13617902941839949718d-7 + (0.66799760083972474642d-10 +&
               & 0.22574701262222222222d-12 * t) * t) * t) * t) * t) * t
       case(96)
          t = 2*y100 - 193
          erfcx_y100 = 0.85525144775685126237d0 + (0.18904632212547561026d-1 + (0.23764237370371255638d-3 +&
               & (0.21461248251306387979d-5 + (0.14299555071870523786d-7 + (0.69543803864694171934d-10 +&
               & 0.23158593688888888889d-12 * t) * t) * t) * t) * t) * t
       case(97)
          t = 2*y100 - 195
          erfcx_y100 = 0.89402868170849933734d0 + (0.19881418399127202569d-1 + (0.25086793128395995798d-3 +&
               & (0.22633402747585233180d-5 + (0.15008997042116532283d-7 + (0.72357609075043941261d-10 +&
               & 0.23737194737777777778d-12 * t) * t) * t) * t) * t) * t
       case(98)
          t = 2*y100 - 197
          erfcx_y100 = 0.93481333942870796363d0 + (0.20912536329780368893d-1 + (0.26481403465998477969d-3 +&
               & (0.23863447359754921676d-5 + (0.15746923065472184451d-7 + (0.75240468141720143653d-10 +&
               & 0.24309291271111111111d-12 * t) * t) * t) * t) * t) * t
       case(99)
          t = 2*y100 - 199
          erfcx_y100 = 0.97771701335885035464d0 + (0.22000938572830479551d-1 + (0.27951610702682383001d-3 +&
               & (0.25153688325245314530d-5 + (0.16514019547822821453d-7 + (0.78191526829368231251d-10 +&
               & 0.24873652355555555556d-12 * t) * t) * t) * t) * t) * t
       case default
          ! we only get here if y = 1, i.e. |x| < 4*eps, in which case
          ! erfcx is within 1d-15 of 1..
          erfcx_y100 = 1.d0
       end select
     end function erfcx_y100
     !-----------------------------------------------------------------------
     real(dp) pure elemental function erfcx(x)
       real(dp), intent(in) :: x
       if (x >= 0) then
          if (x > 50) then ! continued-fraction expansion is faster
             if (x > 5e7) then ! 1-term expansion, important to avoid overflow
                erfcx = ispi / x
                return
             end if
             ! 5-term expansion (rely on compiler for CSE), simplified from:
             !     ispi / (x+0.5/(x+1/(x+1.5/(x+2/x))))
             erfcx = ispi*(x**2 * (x**2+4.5) + 2) / (x * (x**2 * (x**2+5) + 3.75))
             return
          end if
          erfcx = erfcx_y100(400/(4+x))
       else if (x < -26.7) then
          erfcx = HUGE(1.d0)
       else if (x < -6.1) then
          erfcx = 2*exp(x**2)
       else
          erfcx = 2*exp(x**2) - erfcx_y100(400/(4-x))
       end if
     end function erfcx
     !-----------------------------------------------------------------------
     ! Compute a scaled Dawson integral
     !           ImFaddeeva_w(x) = 2*Dawson(x)/sqrt(pi)
     !  equivalent to the imaginary part w(x) for real x.
     !
     !  Uses methods similar to the erfcx calculation above: continued fractions
     !  for large |x|, a lookup table of Chebyshev polynomials for smaller |x|,
     !  and finally a Taylor expansion for |x|<0.01.
     !
     !  Steven G. Johnson, October 2012.

     !  Given y100=100*y, where y = 1/(1+x) for x >= 0, compute ImFaddeeva_w(x).
     !  Uses a look-up table of 100 different Chebyshev polynomials
     !  for y intervals [0,0.01], [0.01,0.02], ...., [0.99,1], generated
     !  with the help of Maple and a little shell script.   This allows
     !  the Chebyshev polynomials to be of significantly lower degree (about 1/30)
     !  compared to fitting the whole [0,1] interval with a single polynomial.
     !-----------------------------------------------------------------------
     real(dp) pure elemental function ImFaddeeva_w_y100(y100,x) result(ImF)
       real(dp), intent(in) :: y100, x
       real(dp) :: t, x2
       select case(int(y100))
       case(0)
          t = 2*y100 - 1
          ImF = 0.28351593328822191546d-2 + (0.28494783221378400759d-2 + (0.14427470563276734183d-4 + (0.10939723080231588129d-6 +&
               & (0.92474307943275042045d-9 + (0.89128907666450075245d-11 + 0.92974121935111111110d-13 * t) * t) * t) * t) * t) * t
       case(1)
          t = 2*y100 - 3
          ImF = 0.85927161243940350562d-2 + (0.29085312941641339862d-2 + (0.15106783707725582090d-4 + (0.11716709978531327367d-6 +&
               & (0.10197387816021040024d-8 + (0.10122678863073360769d-10 + 0.10917479678400000000d-12 * t) * t) * t) * t) * t) * t
       case(2)
          t = 2*y100 - 5
          ImF = 0.14471159831187703054d-1 + (0.29703978970263836210d-2 + (0.15835096760173030976d-4 + (0.12574803383199211596d-6 +&
               & (0.11278672159518415848d-8 + (0.11547462300333495797d-10 + 0.12894535335111111111d-12 * t) * t) * t) * t) * t) * t
       case(3)
          t = 2*y100 - 7
          ImF = 0.20476320420324610618d-1 + (0.30352843012898665856d-2 + (0.16617609387003727409d-4 + (0.13525429711163116103d-6 +&
               & (0.12515095552507169013d-8 + (0.13235687543603382345d-10 + 0.15326595042666666667d-12 * t) * t) * t) * t) * t) * t
       case(4)
          t = 2*y100 - 9
          ImF = 0.26614461952489004566d-1 + (0.31034189276234947088d-2 + (0.17460268109986214274d-4 + (0.14582130824485709573d-6 +&
               & (0.13935959083809746345d-8 + (0.15249438072998932900d-10 + 0.18344741882133333333d-12 * t) * t) * t) * t) * t) * t
       case(5)
          t = 2*y100 - 11
          ImF = 0.32892330248093586215d-1 + (0.31750557067975068584d-2 + (0.18369907582308672632d-4 + (0.15761063702089457882d-6 +&
               & (0.15577638230480894382d-8 + (0.17663868462699097951d-10 + (0.22126732680711111111d-12 +&
               & 0.30273474177737853668d-14 * t) * t) * t) * t) * t) * t) * t
       case(6)
          t = 2*y100 - 13
          ImF = 0.39317207681134336024d-1 + (0.32504779701937539333d-2 + (0.19354426046513400534d-4 + (0.17081646971321290539d-6 +&
               & (0.17485733959327106250d-8 + (0.20593687304921961410d-10 + (0.26917401949155555556d-12 +&
               & 0.38562123837725712270d-14 * t) * t) * t) * t) * t) * t) * t
       case(7)
          t = 2*y100 - 15
          ImF = 0.45896976511367738235d-1 + (0.33300031273110976165d-2 + (0.20423005398039037313d-4 + (0.18567412470376467303d-6 +&
               & (0.19718038363586588213d-8 + (0.24175006536781219807d-10 + (0.33059982791466666666d-12 +&
               & 0.49756574284439426165d-14 * t) * t) * t) * t) * t) * t) * t
       case(8)
          t = 2*y100 - 17
          ImF = 0.52640192524848962855d-1 + (0.34139883358846720806d-2 + (0.21586390240603337337d-4 + (0.20247136501568904646d-6 +&
               & (0.22348696948197102935d-8 + (0.28597516301950162548d-10 + (0.41045502119111111110d-12 +&
               & 0.65151614515238361946d-14 * t) * t) * t) * t) * t) * t) * t
       case(9)
          t = 2*y100 - 19
          ImF = 0.59556171228656770456d-1 + (0.35028374386648914444d-2 + (0.22857246150998562824d-4 + (0.22156372146525190679d-6 +&
               & (0.25474171590893813583d-8 + (0.34122390890697400584d-10 + (0.51593189879111111110d-12 +&
               & 0.86775076853908006938d-14 * t) * t) * t) * t) * t) * t) * t
       case(10)
          t = 2*y100 - 21
          ImF = 0.66655089485108212551d-1 + (0.35970095381271285568d-2 + (0.24250626164318672928d-4 + (0.24339561521785040536d-6 +&
               & (0.29221990406518411415d-8 + (0.41117013527967776467d-10 + (0.65786450716444444445d-12 +&
               & 0.11791885745450623331d-13 * t) * t) * t) * t) * t) * t) * t
       case(11)
          t = 2*y100 - 23
          ImF = 0.73948106345519174661d-1 + (0.36970297216569341748d-2 + (0.25784588137312868792d-4 + (0.26853012002366752770d-6 +&
               & (0.33763958861206729592d-8 + (0.50111549981376976397d-10 + (0.85313857496888888890d-12 +&
               & 0.16417079927706899860d-13 * t) * t) * t) * t) * t) * t) * t
       case(12)
          t = 2*y100 - 25
          ImF = 0.81447508065002963203d-1 + (0.38035026606492705117d-2 + (0.27481027572231851896d-4 + (0.29769200731832331364d-6 +&
               & (0.39336816287457655076d-8 + (0.61895471132038157624d-10 + (0.11292303213511111111d-11 +&
               & 0.23558532213703884304d-13 * t) * t) * t) * t) * t) * t) * t
       case(13)
          t = 2*y100 - 27
          ImF = 0.89166884027582716628d-1 + (0.39171301322438946014d-2 + (0.29366827260422311668d-4 + (0.33183204390350724895d-6 +&
               & (0.46276006281647330524d-8 + (0.77692631378169813324d-10 + (0.15335153258844444444d-11 +&
               & 0.35183103415916026911d-13 * t) * t) * t) * t) * t) * t) * t
       case(14)
          t = 2*y100 - 29
          ImF = 0.97121342888032322019d-1 + (0.40387340353207909514d-2 + (0.31475490395950776930d-4 + (0.37222714227125135042d-6 +&
               & (0.55074373178613809996d-8 + (0.99509175283990337944d-10 + (0.21552645758222222222d-11 +&
               & 0.55728651431872687605d-13 * t) * t) * t) * t) * t) * t) * t
       case(15)
          t = 2*y100 - 31
          ImF = 0.10532778218603311137d0 + (0.41692873614065380607d-2 + (0.33849549774889456984d-4 + (0.42064596193692630143d-6 +&
               & (0.66494579697622432987d-8 + (0.13094103581931802337d-9 + (0.31896187409777777778d-11 +&
               & 0.97271974184476560742d-13 * t) * t) * t) * t) * t) * t) * t
       case(16)
          t = 2*y100 - 33
          ImF = 0.11380523107427108222d0 + (0.43099572287871821013d-2 + (0.36544324341565929930d-4 + (0.47965044028581857764d-6 +&
               & (0.81819034238463698796d-8 + (0.17934133239549647357d-9 + (0.50956666166186293627d-11 +&
               & (0.18850487318190638010d-12 + 0.79697813173519853340d-14 * t) * t) * t) * t) * t) * t) * t) * t
       case(17)
          t = 2*y100 - 35
          ImF = 0.12257529703447467345d0 + (0.44621675710026986366d-2 + (0.39634304721292440285d-4 + (0.55321553769873381819d-6 +&
               & (0.10343619428848520870d-7 + (0.26033830170470368088d-9 + (0.87743837749108025357d-11 +&
               & (0.34427092430230063401d-12 + 0.10205506615709843189d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(18)
          t = 2*y100 - 37
          ImF = 0.13166276955656699478d0 + (0.46276970481783001803d-2 + (0.43225026380496399310d-4 + (0.64799164020016902656d-6 +&
               & (0.13580082794704641782d-7 + (0.39839800853954313927d-9 + (0.14431142411840000000d-10 +&
               & 0.42193457308830027541d-12 * t) * t) * t) * t) * t) * t) * t
       case(19)
          t = 2*y100 - 39
          ImF = 0.14109647869803356475d0 + (0.48088424418545347758d-2 + (0.47474504753352150205d-4 + (0.77509866468724360352d-6 +&
               & (0.18536851570794291724d-7 + (0.60146623257887570439d-9 + (0.18533978397305276318d-10 +&
               & (0.41033845938901048380d-13 - 0.46160680279304825485d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(20)
          t = 2*y100 - 41
          ImF = 0.15091057940548936603d0 + (0.50086864672004685703d-2 + (0.52622482832192230762d-4 + (0.95034664722040355212d-6 +&
               & (0.25614261331144718769d-7 + (0.80183196716888606252d-9 + (0.12282524750534352272d-10 +&
               & (-0.10531774117332273617d-11 - 0.86157181395039646412d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(21)
          t = 2*y100 - 43
          ImF = 0.16114648116017010770d0 + (0.52314661581655369795d-2 + (0.59005534545908331315d-4 + (0.11885518333915387760d-5 +&
               & (0.33975801443239949256d-7 + (0.82111547144080388610d-9 + (-0.12357674017312854138d-10 +&
               & (-0.24355112256914479176d-11 - 0.75155506863572930844d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(22)
          t = 2*y100 - 45
          ImF = 0.17185551279680451144d0 + (0.54829002967599420860d-2 + (0.67013226658738082118d-4 + (0.14897400671425088807d-5 +&
               & (0.40690283917126153701d-7 + (0.44060872913473778318d-9 + (-0.52641873433280000000d-10 -&
               & 0.30940587864543343124d-11 * t) * t) * t) * t) * t) * t) * t
       case(23)
          t = 2*y100 - 47
          ImF = 0.18310194559815257381d0 + (0.57701559375966953174d-2 + (0.76948789401735193483d-4 + (0.18227569842290822512d-5 +&
               & (0.41092208344387212276d-7 + (-0.44009499965694442143d-9 + (-0.92195414685628803451d-10 +&
               & (-0.22657389705721753299d-11 + 0.10004784908106839254d-12 * t) * t) * t) * t) * t) * t) * t) * t
       case(24)
          t = 2*y100 - 49
          ImF = 0.19496527191546630345d0 + (0.61010853144364724856d-2 + (0.88812881056342004864d-4 + (0.21180686746360261031d-5 +&
               & (0.30652145555130049203d-7 + (-0.16841328574105890409d-8 + (-0.11008129460612823934d-9 +&
               & (-0.12180794204544515779d-12 + 0.15703325634590334097d-12 * t) * t) * t) * t) * t) * t) * t) * t
       case(25)
          t = 2*y100 - 51
          ImF = 0.20754006813966575720d0 + (0.64825787724922073908d-2 + (0.10209599627522311893d-3 + (0.22785233392557600468d-5 +&
               & (0.73495224449907568402d-8 + (-0.29442705974150112783d-8 + (-0.94082603434315016546d-10 +&
               & (0.23609990400179321267d-11 + 0.14141908654269023788d-12 * t) * t) * t) * t) * t) * t) * t) * t
       case(26)
          t = 2*y100 - 53
          ImF = 0.22093185554845172146d0 + (0.69182878150187964499d-2 + (0.11568723331156335712d-3 + (0.22060577946323627739d-5 +&
               & (-0.26929730679360840096d-7 + (-0.38176506152362058013d-8 + (-0.47399503861054459243d-10 +&
               & (0.40953700187172127264d-11 + 0.69157730376118511127d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(27)
          t = 2*y100 - 55
          ImF = 0.23524827304057813918d0 + (0.74063350762008734520d-2 + (0.12796333874615790348d-3 + (0.18327267316171054273d-5 +&
               & (-0.66742910737957100098d-7 + (-0.40204740975496797870d-8 + (0.14515984139495745330d-10 +&
               & (0.44921608954536047975d-11 - 0.18583341338983776219d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(28)
          t = 2*y100 - 57
          ImF = 0.25058626331812744775d0 + (0.79377285151602061328d-2 + (0.13704268650417478346d-3 + (0.11427511739544695861d-5 +&
               & (-0.10485442447768377485d-6 + (-0.34850364756499369763d-8 + (0.72656453829502179208d-10 +&
               & (0.36195460197779299406d-11 - 0.84882136022200714710d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(29)
          t = 2*y100 - 59
          ImF = 0.26701724900280689785d0 + (0.84959936119625864274d-2 + (0.14112359443938883232d-3 + (0.17800427288596909634d-6 +&
               & (-0.13443492107643109071d-6 + (-0.23512456315677680293d-8 + (0.11245846264695936769d-9 +&
               & (0.19850501334649565404d-11 - 0.11284666134635050832d-12 * t) * t) * t) * t) * t) * t) * t) * t
       case(30)
          t = 2*y100 - 61
          ImF = 0.28457293586253654144d0 + (0.90581563892650431899d-2 + (0.13880520331140646738d-3 + (-0.97262302362522896157d-6 +&
               & (-0.15077100040254187366d-6 + (-0.88574317464577116689d-9 + (0.12760311125637474581d-9 +&
               & (0.20155151018282695055d-12 - 0.10514169375181734921d-12 * t) * t) * t) * t) * t) * t) * t) * t
       case(31)
          t = 2*y100 - 63
          ImF = 0.30323425595617385705d0 + (0.95968346790597422934d-2 + (0.12931067776725883939d-3 + (-0.21938741702795543986d-5 +&
               & (-0.15202888584907373963d-6 + (0.61788350541116331411d-9 + (0.11957835742791248256d-9 +&
               & (-0.12598179834007710908d-11 - 0.75151817129574614194d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(32)
          t = 2*y100 - 65
          ImF = 0.32292521181517384379d0 + (0.10082957727001199408d-1 + (0.11257589426154962226d-3 + (-0.33670890319327881129d-5 +&
               & (-0.13910529040004008158d-6 + (0.19170714373047512945d-8 + (0.94840222377720494290d-10 +&
               & (-0.21650018351795353201d-11 - 0.37875211678024922689d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(33)
          t = 2*y100 - 67
          ImF = 0.34351233557911753862d0 + (0.10488575435572745309d-1 + (0.89209444197248726614d-4 + (-0.43893459576483345364d-5 +&
               & (-0.11488595830450424419d-6 + (0.28599494117122464806d-8 + (0.61537542799857777779d-10 -&
               & 0.24935749227658002212d-11 * t) * t) * t) * t) * t) * t) * t
       case(34)
          t = 2*y100 - 69
          ImF = 0.36480946642143669093d0 + (0.10789304203431861366d-1 + (0.60357993745283076834d-4 + (-0.51855862174130669389d-5 +&
               & (-0.83291664087289801313d-7 + (0.33898011178582671546d-8 + (0.27082948188277716482d-10 +&
               & (-0.23603379397408694974d-11 + 0.19328087692252869842d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(35)
          t = 2*y100 - 71
          ImF = 0.38658679935694939199d0 + (0.10966119158288804999d-1 + (0.27521612041849561426d-4 + (-0.57132774537670953638d-5 +&
               & (-0.48404772799207914899d-7 + (0.35268354132474570493d-8 + (-0.32383477652514618094d-11 +&
               & (-0.19334202915190442501d-11 + 0.32333189861286460270d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(36)
          t = 2*y100 - 73
          ImF = 0.40858275583808707870d0 + (0.11006378016848466550d-1 + (-0.76396376685213286033d-5 + (-0.59609835484245791439d-5&
               & + (-0.13834610033859313213d-7 + (0.33406952974861448790d-8 + (-0.26474915974296612559d-10 +&
               & (-0.13750229270354351983d-11 + 0.36169366979417390637d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(37)
          t = 2*y100 - 75
          ImF = 0.43051714914006682977d0 + (0.10904106549500816155d-1 + (-0.43477527256787216909d-4 + (-0.59429739547798343948d-5&
               & + (0.17639200194091885949d-7 + (0.29235991689639918688d-8 + (-0.41718791216277812879d-10 +&
               & (-0.81023337739508049606d-12 + 0.33618915934461994428d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(38)
          t = 2*y100 - 77
          ImF = 0.45210428135559607406d0 + (0.10659670756384400554d-1 + (-0.78488639913256978087d-4 + (-0.56919860886214735936d-5&
               & + (0.44181850467477733407d-7 + (0.23694306174312688151d-8 + (-0.49492621596685443247d-10 +&
               & (-0.31827275712126287222d-12 + 0.27494438742721623654d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(39)
          t = 2*y100 - 79
          ImF = 0.47306491195005224077d0 + (0.10279006119745977570d-1 + (-0.11140268171830478306d-3 + (-0.52518035247451432069d-5&
               & + (0.64846898158889479518d-7 + (0.17603624837787337662d-8 + (-0.51129481592926104316d-10 +&
               & (0.62674584974141049511d-13 + 0.20055478560829935356d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(40)
          t = 2*y100 - 81
          ImF = 0.49313638965719857647d0 + (0.97725799114772017662d-2 + (-0.14122854267291533334d-3 + (-0.46707252568834951907d-5&
               & + (0.79421347979319449524d-7 + (0.11603027184324708643d-8 + (-0.48269605844397175946d-10 +&
               & (0.32477251431748571219d-12 + 0.12831052634143527985d-13 * t) * t) * t) * t) * t) * t) * t) * t
       case(41)
          t = 2*y100 - 83
          ImF = 0.51208057433416004042d0 + (0.91542422354009224951d-2 + (-0.16726530230228647275d-3 + (-0.39964621752527649409d-5&
               & + (0.88232252903213171454d-7 + (0.61343113364949928501d-9 + (-0.42516755603130443051d-10 +&
               & (0.47910437172240209262d-12 + 0.66784341874437478953d-14 * t) * t) * t) * t) * t) * t) * t) * t
       case(42)
          t = 2*y100 - 85
          ImF = 0.52968945458607484524d0 + (0.84400880445116786088d-2 + (-0.18908729783854258774d-3 + (-0.32725905467782951931d-5&
               & + (0.91956190588652090659d-7 + (0.14593989152420122909d-9 + (-0.35239490687644444445d-10 +&
               & 0.54613829888448694898d-12 * t) * t) * t) * t) * t) * t) * t
       case(43)
          t = 2*y100 - 87
          ImF = 0.54578857454330070965d0 + (0.76474155195880295311d-2 + (-0.20651230590808213884d-3 + (-0.25364339140543131706d-5&
               & + (0.91455367999510681979d-7 + (-0.23061359005297528898d-9 + (-0.27512928625244444444d-10 +&
               & 0.54895806008493285579d-12 * t) * t) * t) * t) * t) * t) * t
       case(44)
          t = 2*y100 - 89
          ImF = 0.56023851910298493910d0 + (0.67938321739997196804d-2 + (-0.21956066613331411760d-3 + (-0.18181127670443266395d-5&
               & + (0.87650335075416845987d-7 + (-0.51548062050366615977d-9 + (-0.20068462174044444444d-10 +&
               & 0.50912654909758187264d-12 * t) * t) * t) * t) * t) * t) * t
       case(45)
          t = 2*y100 - 91
          ImF = 0.57293478057455721150d0 + (0.58965321010394044087d-2 + (-0.22841145229276575597d-3 + (-0.11404605562013443659d-5&
               & + (0.81430290992322326296d-7 + (-0.71512447242755357629d-9 + (-0.13372664928000000000d-10 +&
               & 0.44461498336689298148d-12 * t) * t) * t) * t) * t) * t) * t
       case(46)
          t = 2*y100 - 93
          ImF = 0.58380635448407827360d0 + (0.49717469530842831182d-2 + (-0.23336001540009645365d-3 + (-0.51952064448608850822d-6&
               & + (0.73596577815411080511d-7 + (-0.84020916763091566035d-9 + (-0.76700972702222222221d-11 +&
               & 0.36914462807972467044d-12 * t) * t) * t) * t) * t) * t) * t
       case(47)
          t = 2*y100 - 95
          ImF = 0.59281340237769489597d0 + (0.40343592069379730568d-2 + (-0.23477963738658326185d-3 + (0.34615944987790224234d-7 +&
               & (0.64832803248395814574d-7 + (-0.90329163587627007971d-9 + (-0.30421940400000000000d-11 +&
               & 0.29237386653743536669d-12 * t) * t) * t) * t) * t) * t) * t
       case(48)
          t = 2*y100 - 97
          ImF = 0.59994428743114271918d0 + (0.30976579788271744329d-2 + (-0.23308875765700082835d-3 + (0.51681681023846925160d-6 +&
               & (0.55694594264948268169d-7 + (-0.91719117313243464652d-9 + (0.53982743680000000000d-12 +&
               & 0.22050829296187771142d-12 * t) * t) * t) * t) * t) * t) * t
       case(49)
          t = 2*y100 - 99
          ImF = 0.60521224471819875444d0 + (0.21732138012345456060d-2 + (-0.22872428969625997456d-3 + (0.92588959922653404233d-6 +&
               & (0.46612665806531930684d-7 + (-0.89393722514414153351d-9 + (0.31718550353777777778d-11 +&
               & 0.15705458816080549117d-12 * t) * t) * t) * t) * t) * t) * t
       case(50)
          t = 2*y100 - 101
          ImF = 0.60865189969791123620d0 + (0.12708480848877451719d-2 + (-0.22212090111534847166d-3 + (0.12636236031532793467d-5 +&
               & (0.37904037100232937574d-7 + (-0.84417089968101223519d-9 + (0.49843180828444444445d-11 +&
               & 0.10355439441049048273d-12 * t) * t) * t) * t) * t) * t) * t
       case(51)
          t = 2*y100 - 103
          ImF = 0.61031580103499200191d0 + (0.39867436055861038223d-3 + (-0.21369573439579869291d-3 + (0.15339402129026183670d-5 +&
               & (0.29787479206646594442d-7 + (-0.77687792914228632974d-9 + (0.61192452741333333334d-11 +&
               & 0.60216691829459295780d-13 * t) * t) * t) * t) * t) * t) * t
       case(52)
          t = 2*y100 - 105
          ImF = 0.61027109047879835868d0 + (-0.43680904508059878254d-3 + (-0.20383783788303894442d-3 + (0.17421743090883439959d-5&
               & + (0.22400425572175715576d-7 + (-0.69934719320045128997d-9 + (0.67152759655111111110d-11 +&
               & 0.26419960042578359995d-13 * t) * t) * t) * t) * t) * t) * t
       case(53)
          t = 2*y100 - 107
          ImF = 0.60859639489217430521d0 + (-0.12305921390962936873d-2 + (-0.19290150253894682629d-3 + (0.18944904654478310128d-5&
               & + (0.15815530398618149110d-7 + (-0.61726850580964876070d-9 + 0.68987888999111111110d-11 * t) * t) * t) * t) * t) &
               &* t
       case(54)
          t = 2*y100 - 109
          ImF = 0.60537899426486075181d0 + (-0.19790062241395705751d-2 + (-0.18120271393047062253d-3 + (0.19974264162313241405d-5&
               & + (0.10055795094298172492d-7 + (-0.53491997919318263593d-9 + (0.67794550295111111110d-11 -&
               & 0.17059208095741511603d-13 * t) * t) * t) * t) * t) * t) * t
       case(55)
          t = 2*y100 - 111
          ImF = 0.60071229457904110537d0 + (-0.26795676776166354354d-2 + (-0.16901799553627508781d-3 + (0.20575498324332621581d-5&
               & + (0.51077165074461745053d-8 + (-0.45536079828057221858d-9 + (0.64488005516444444445d-11 -&
               & 0.29311677573152766338d-13 * t) * t) * t) * t) * t) * t) * t
       case(56)
          t = 2*y100 - 113
          ImF = 0.59469361520112714738d0 + (-0.33308208190600993470d-2 + (-0.15658501295912405679d-3 + (0.20812116912895417272d-5&
               & + (0.93227468760614182021d-9 + (-0.38066673740116080415d-9 + (0.59806790359111111110d-11 -&
               & 0.36887077278950440597d-13 * t) * t) * t) * t) * t) * t) * t
       case(57)
          t = 2*y100 - 115
          ImF = 0.58742228631775388268d0 + (-0.39321858196059227251d-2 + (-0.14410441141450122535d-3 + (0.20743790018404020716d-5&
               & + (-0.25261903811221913762d-8 + (-0.31212416519526924318d-9 + (0.54328422462222222221d-11 -&
               & 0.40864152484979815972d-13 * t) * t) * t) * t) * t) * t) * t
       case(58)
          t = 2*y100 - 117
          ImF = 0.57899804200033018447d0 + (-0.44838157005618913447d-2 + (-0.13174245966501437965d-3 + (0.20425306888294362674d-5&
               & + (-0.53330296023875447782d-8 + (-0.25041289435539821014d-9 + (0.48490437205333333334d-11 -&
               & 0.42162206939169045177d-13 * t) * t) * t) * t) * t) * t) * t
       case(59)
          t = 2*y100 - 119
          ImF = 0.56951968796931245974d0 + (-0.49864649488074868952d-2 + (-0.11963416583477567125d-3 + (0.19906021780991036425d-5&
               & + (-0.75580140299436494248d-8 + (-0.19576060961919820491d-9 + (0.42613011928888888890d-11 -&
               & 0.41539443304115604377d-13 * t) * t) * t) * t) * t) * t) * t
       case(60)
          t = 2*y100 - 121
          ImF = 0.55908401930063918964d0 + (-0.54413711036826877753d-2 + (-0.10788661102511914628d-3 + (0.19229663322982839331d-5&
               & + (-0.92714731195118129616d-8 + (-0.14807038677197394186d-9 + (0.36920870298666666666d-11 -&
               & 0.39603726688419162617d-13 * t) * t) * t) * t) * t) * t) * t
       case(61)
          t = 2*y100 - 123
          ImF = 0.54778496152925675315d0 + (-0.58501497933213396670d-2 + (-0.96582314317855227421d-4 + (0.18434405235069270228d-5&
               & + (-0.10541580254317078711d-7 + (-0.10702303407788943498d-9 + (0.31563175582222222222d-11 -&
               & 0.36829748079110481422d-13 * t) * t) * t) * t) * t) * t) * t
       case(62)
          t = 2*y100 - 125
          ImF = 0.53571290831682823999d0 + (-0.62147030670760791791d-2 + (-0.85782497917111760790d-4 + (0.17553116363443470478d-5&
               & + (-0.11432547349815541084d-7 + (-0.72157091369041330520d-10 + (0.26630811607111111111d-11 -&
               & 0.33578660425893164084d-13 * t) * t) * t) * t) * t) * t) * t
       case(63)
          t = 2*y100 - 127
          ImF = 0.52295422962048434978d0 + (-0.65371404367776320720d-2 + (-0.75530164941473343780d-4 + (0.16613725797181276790d-5&
               & + (-0.12003521296598910761d-7 + (-0.42929753689181106171d-10 + (0.22170894940444444444d-11 -&
               & 0.30117697501065110505d-13 * t) * t) * t) * t) * t) * t) * t
       case(64)
          t = 2*y100 - 129
          ImF = 0.50959092577577886140d0 + (-0.68197117603118591766d-2 + (-0.65852936198953623307d-4 + (0.15639654113906716939d-5&
               & + (-0.12308007991056524902d-7 + (-0.18761997536910939570d-10 + (0.18198628922666666667d-11 -&
               & 0.26638355362285200932d-13 * t) * t) * t) * t) * t) * t) * t
       case(65)
          t = 2*y100 - 131
          ImF = 0.49570040481823167970d0 + (-0.70647509397614398066d-2 + (-0.56765617728962588218d-4 + (0.14650274449141448497d-5&
               & + (-0.12393681471984051132d-7 + (0.92904351801168955424d-12 + (0.14706755960177777778d-11 -&
               & 0.23272455351266325318d-13 * t) * t) * t) * t) * t) * t) * t
       case(66)
          t = 2*y100 - 133
          ImF = 0.48135536250935238066d0 + (-0.72746293327402359783d-2 + (-0.48272489495730030780d-4 + (0.13661377309113939689d-5&
               & + (-0.12302464447599382189d-7 + (0.16707760028737074907d-10 + (0.11672928324444444444d-11 -&
               & 0.20105801424709924499d-13 * t) * t) * t) * t) * t) * t) * t
       case(67)
          t = 2*y100 - 135
          ImF = 0.46662374675511439448d0 + (-0.74517177649528487002d-2 + (-0.40369318744279128718d-4 + (0.12685621118898535407d-5&
               & + (-0.12070791463315156250d-7 + (0.29105507892605823871d-10 + (0.90653314645333333334d-12 -&
               & 0.17189503312102982646d-13 * t) * t) * t) * t) * t) * t) * t
       case(68)
          t = 2*y100 - 137
          ImF = 0.45156879030168268778d0 + (-0.75983560650033817497d-2 + (-0.33045110380705139759d-4 + (0.11732956732035040896d-5&
               & + (-0.11729986947158201869d-7 + (0.38611905704166441308d-10 + (0.68468768305777777779d-12 -&
               & 0.14549134330396754575d-13 * t) * t) * t) * t) * t) * t) * t
       case(69)
          t = 2*y100 - 139
          ImF = 0.43624909769330896904d0 + (-0.77168291040309554679d-2 + (-0.26283612321339907756d-4 + (0.10811018836893550820d-5&
               & + (-0.11306707563739851552d-7 + (0.45670446788529607380d-10 + (0.49782492549333333334d-12 -&
               & 0.12191983967561779442d-13 * t) * t) * t) * t) * t) * t) * t
       case(70)
          t = 2*y100 - 141
          ImF = 0.42071877443548481181d0 + (-0.78093484015052730097d-2 + (-0.20064596897224934705d-4 + (0.99254806680671890766d-6&
               & + (-0.10823412088884741451d-7 + (0.50677203326904716247d-10 + (0.34200547594666666666d-12 -&
               & 0.10112698698356194618d-13 * t) * t) * t) * t) * t) * t) * t
       case(71)
          t = 2*y100 - 143
          ImF = 0.40502758809710844280d0 + (-0.78780384460872937555d-2 + (-0.14364940764532853112d-4 + (0.90803709228265217384d-6&
               & + (-0.10298832847014466907d-7 + (0.53981671221969478551d-10 + (0.21342751381333333333d-12 -&
               & 0.82975901848387729274d-14 * t) * t) * t) * t) * t) * t) * t
       case(72)
          t = 2*y100 - 145
          ImF = 0.38922115269731446690d0 + (-0.79249269708242064120d-2 + (-0.91595258799106970453d-5 + (0.82783535102217576495d-6&
               & + (-0.97484311059617744437d-8 + (0.55889029041660225629d-10 + (0.10851981336888888889d-12 -&
               & 0.67278553237853459757d-14 * t) * t) * t) * t) * t) * t) * t
       case(73)
          t = 2*y100 - 147
          ImF = 0.37334112915460307335d0 + (-0.79519385109223148791d-2 + (-0.44219833548840469752d-5 + (0.75209719038240314732d-6&
               & + (-0.91848251458553190451d-8 + (0.56663266668051433844d-10 + (0.23995894257777777778d-13 -&
               & 0.53819475285389344313d-14 * t) * t) * t) * t) * t) * t) * t
       case(74)
          t = 2*y100 - 149
          ImF = 0.35742543583374223085d0 + (-0.79608906571527956177d-2 + (-0.12530071050975781198d-6 + (0.68088605744900552505d-6&
               & + (-0.86181844090844164075d-8 + (0.56530784203816176153d-10 + (-0.43120012248888888890d-13 -&
               & 0.42372603392496813810d-14 * t) * t) * t) * t) * t) * t) * t
       case(75)
          t = 2*y100 - 151
          ImF = 0.34150846431979618536d0 + (-0.79534924968773806029d-2 + (0.37576885610891515813d-5 + (0.61419263633090524326d-6 +&
               & (-0.80565865409945960125d-8 + (0.55684175248749269411d-10 + (-0.95486860764444444445d-13 -&
               & 0.32712946432984510595d-14 * t) * t) * t) * t) * t) * t) * t
       case(76)
          t = 2*y100 - 153
          ImF = 0.32562129649136346824d0 + (-0.79313448067948884309d-2 + (0.72539159933545300034d-5 + (0.55195028297415503083d-6 +&
               & (-0.75063365335570475258d-8 + (0.54281686749699595941d-10 - 0.13545424295111111111d-12 * t) * t) * t) * t) * t) *&
               & t
       case(77)
          t = 2*y100 - 155
          ImF = 0.30979191977078391864d0 + (-0.78959416264207333695d-2 + (0.10389774377677210794d-4 + (0.49404804463196316464d-6 +&
               & (-0.69722488229411164685d-8 + (0.52469254655951393842d-10 - 0.16507860650666666667d-12 * t) * t) * t) * t) * t) *&
               & t
       case(78)
          t = 2*y100 - 157
          ImF = 0.29404543811214459904d0 + (-0.78486728990364155356d-2 + (0.13190885683106990459d-4 + (0.44034158861387909694d-6 +&
               & (-0.64578942561562616481d-8 + (0.50354306498006928984d-10 - 0.18614473550222222222d-12 * t) * t) * t) * t) * t) *&
               & t
       case(79)
          t = 2*y100 - 159
          ImF = 0.27840427686253660515d0 + (-0.77908279176252742013d-2 + (0.15681928798708548349d-4 + (0.39066226205099807573d-6 +&
               & (-0.59658144820660420814d-8 + (0.48030086420373141763d-10 - 0.20018995173333333333d-12 * t) * t) * t) * t) * t) *&
               & t
       case(80)
          t = 2*y100 - 161
          ImF = 0.26288838011163800908d0 + (-0.77235993576119469018d-2 + (0.17886516796198660969d-4 + (0.34482457073472497720d-6 +&
               & (-0.54977066551955420066d-8 + (0.45572749379147269213d-10 - 0.20852924954666666667d-12 * t) * t) * t) * t) * t) *&
               & t
       case(81)
          t = 2*y100 - 163
          ImF = 0.24751539954181029717d0 + (-0.76480877165290370975d-2 + (0.19827114835033977049d-4 + (0.30263228619976332110d-6 +&
               & (-0.50545814570120129947d-8 + (0.43043879374212005966d-10 - 0.21228012028444444444d-12 * t) * t) * t) * t) * t) *&
               & t
       case(82)
          t = 2*y100 - 165
          ImF = 0.23230087411688914593d0 + (-0.75653060136384041587d-2 + (0.21524991113020016415d-4 + (0.26388338542539382413d-6 +&
               & (-0.46368974069671446622d-8 + (0.40492715758206515307d-10 - 0.21238627815111111111d-12 * t) * t) * t) * t) * t) *&
               & t
       case(83)
          t = 2*y100 - 167
          ImF = 0.21725840021297341931d0 + (-0.74761846305979730439d-2 + (0.23000194404129495243d-4 + (0.22837400135642906796d-6 +&
               & (-0.42446743058417541277d-8 + (0.37958104071765923728d-10 - 0.20963978568888888889d-12 * t) * t) * t) * t) * t) *&
               & t
       case(84)
          t = 2*y100 - 169
          ImF = 0.20239979200788191491d0 + (-0.73815761980493466516d-2 + (0.24271552727631854013d-4 + (0.19590154043390012843d-6 +&
               & (-0.38775884642456551753d-8 + (0.35470192372162901168d-10 - 0.20470131678222222222d-12 * t) * t) * t) * t) * t) *&
               & t
       case(85)
          t = 2*y100 - 171
          ImF = 0.18773523211558098962d0 + (-0.72822604530339834448d-2 + (0.25356688567841293697d-4 + (0.16626710297744290016d-6 +&
               & (-0.35350521468015310830d-8 + (0.33051896213898864306d-10 - 0.19811844544000000000d-12 * t) * t) * t) * t) * t) *&
               & t
       case(86)
          t = 2*y100 - 173
          ImF = 0.17327341258479649442d0 + (-0.71789490089142761950d-2 + (0.26272046822383820476d-4 + (0.13927732375657362345d-6 +&
               & (-0.32162794266956859603d-8 + (0.30720156036105652035d-10 - 0.19034196304000000000d-12 * t) * t) * t) * t) * t) *&
               & t
       case(87)
          t = 2*y100 - 175
          ImF = 0.15902166648328672043d0 + (-0.70722899934245504034d-2 + (0.27032932310132226025d-4 + (0.11474573347816568279d-6 +&
               & (-0.29203404091754665063d-8 + (0.28487010262547971859d-10 - 0.18174029063111111111d-12 * t) * t) * t) * t) * t) *&
               & t
       case(88)
          t = 2*y100 - 177
          ImF = 0.14498609036610283865d0 + (-0.69628725220045029273d-2 + (0.27653554229160596221d-4 + (0.92493727167393036470d-7 +&
               & (-0.26462055548683583849d-8 + (0.26360506250989943739d-10 - 0.17261211260444444444d-12 * t) * t) * t) * t) * t) *&
               & t
       case(89)
          t = 2*y100 - 179
          ImF = 0.13117165798208050667d0 + (-0.68512309830281084723d-2 + (0.28147075431133863774d-4 + (0.72351212437979583441d-7 +&
               & (-0.23927816200314358570d-8 + (0.24345469651209833155d-10 - 0.16319736960000000000d-12 * t) * t) * t) * t) * t) *&
               & t
       case(90)
          t = 2*y100 - 181
          ImF = 0.11758232561160626306d0 + (-0.67378491192463392927d-2 + (0.28525664781722907847d-4 + (0.54156999310046790024d-7 +&
               & (-0.21589405340123827823d-8 + (0.22444150951727334619d-10 - 0.15368675584000000000d-12 * t) * t) * t) * t) * t) *&
               & t
       case(91)
          t = 2*y100 - 183
          ImF = 0.10422112945361673560d0 + (-0.66231638959845581564d-2 + (0.28800551216363918088d-4 + (0.37758983397952149613d-7 +&
               & (-0.19435423557038933431d-8 + (0.20656766125421362458d-10 - 0.14422990012444444444d-12 * t) * t) * t) * t) * t) *&
               & t
       case(92)
          t = 2*y100 - 185
          ImF = 0.91090275493541084785d-1 + (-0.65075691516115160062d-2 + (0.28982078385527224867d-4 + (0.23014165807643012781d-7&
               & + (-0.17454532910249875958d-8 + (0.18981946442680092373d-10 - 0.13494234691555555556d-12 * t) * t) * t) * t) * t)&
               & * t
       case(93)
          t = 2*y100 - 187
          ImF = 0.78191222288771379358d-1 + (-0.63914190297303976434d-2 + (0.29079759021299682675d-4 + (0.97885458059415717014d-8&
               & + (-0.15635596116134296819d-8 + (0.17417110744051331974d-10 - 0.12591151763555555556d-12 * t) * t) * t) * t) * t)&
               & * t
       case(94)
          t = 2*y100 - 189
          ImF = 0.65524757106147402224d-1 + (-0.62750311956082444159d-2 + (0.29102328354323449795d-4 + (-0.20430838882727954582d-8&
               & + (-0.13967781903855367270d-8 + (0.15958771833747057569d-10 - 0.11720175765333333333d-12 * t) * t) * t) * t) * t)&
               & * t
       case(95)
          t = 2*y100 - 191
          ImF = 0.53091065838453612773d-1 + (-0.61586898417077043662d-2 + (0.29057796072960100710d-4 + (-0.12597414620517987536d-7&
               & + (-0.12440642607426861943d-8 + (0.14602787128447932137d-10 - 0.10885859114666666667d-12 * t) * t) * t) * t) * t)&
               & * t
       case(96)
          t = 2*y100 - 193
          ImF = 0.40889797115352738582d-1 + (-0.60426484889413678200d-2 + (0.28953496450191694606d-4 + (-0.21982952021823718400d-7&
               & + (-0.11044169117553026211d-8 + (0.13344562332430552171d-10 - 0.10091231402844444444d-12 * t) * t) * t) * t) * t)&
               & * t
       case(97)
          t = 2*y100 - 195
          ImF = 0.28920121009594899986d-1 + (-0.59271325915413781788d-2 + (0.28796136372768177423d-4 + (-0.30300382596279568642d-7&
               & + (-0.97688275022802329749d-9 + (0.12179215701512592356d-10 - 0.93380988481777777779d-13 * t) * t) * t) * t) * t)&
               & * t
       case(98)
          t = 2*y100 - 197
          ImF = 0.17180782722617876655d-1 + (-0.58123419543161127769d-2 + (0.28591841095380959666d-4 + (-0.37642963496443667043d-7&
               & + (-0.86055809047367300024d-9 + (0.11101709356762665578d-10 - 0.86272947493333333334d-13 * t) * t) * t) * t) * t)&
               & * t
       case(99,100)
          ! use Taylor expansion for small x (|x| <= 0.010101...)
          !  (2/sqrt(pi)) * (x - 2/3 x^3  + 4/15 x^5  - 8/105 x^7)
          x2 = x*x
          ImF = x * (1.1283791670955125739d0    &
               - x2 * (0.75225277806367504925d0 &
               - x2 * (0.30090111122547001970d0 &
               - x2 * 0.085971746064420005629d0)))
       case default
          ImF = 0.d0 ! never reached for 0 <= y < 101
       end select
     end function ImFaddeeva_w_y100
     !-----------------------------------------------------------------------
     real(dp) pure elemental function ImFaddeeva_w(x)
       real(dp), intent(in) :: x
       if (x >= 0) then
          if (x > 45) then ! continued-fraction expansion is faster
             if (x > 5d7) then ! 1-term expansion, important to avoid overflow
                ImFaddeeva_w = ispi / x
             else
                ! 5-term expansion (rely on compiler for CSE), simplified from:
                !         ispi / (x+0.5/(x+1/(x+1.5/(x+2/x))))
                ImFaddeeva_w = ispi * (x**2 * (x**2-4.5d0) + 2) / (x * (x**2 * (x**2-5) + 3.75d0))
             end if
          else
             ImFaddeeva_w = ImFaddeeva_w_y100(100/(1+x), x)
          end if
       else ! = -ImFaddeeva_w(-x)
          if (x < -45) then ! continued-fraction expansion is faster
             if (x < -5e7) then ! 1-term expansion, important to avoid overflow
                ImFaddeeva_w = ispi / x
             else
                ! 5-term expansion (rely on compiler for CSE), simplified from:
                !         ispi / (x+0.5/(x+1/(x+1.5/(x+2/x))))
                ImFaddeeva_w = ispi * (x**2 * (x**2-4.5d0) + 2) / (x * (x**2 * (x**2-5) + 3.75d0))
             end if
          else
             ImFaddeeva_w = -ImFaddeeva_w_y100(100/(1-x), -x)
          end if
       end if
     end function ImFaddeeva_w
     !-----------------------------------------------------------------------
   end module faddeeva

#ifdef FADDEEVA_W_TEST
! Compile with -DFADDEEVA_W_TEST to compile a little test program

program main
  use faddeeva
  implicit none
  integer, parameter :: dp  = kind(1.d0)
  integer, parameter :: dpc = kind((1.d0,1.d0))
  integer, parameter :: NTST = 46
  complex(dpc), parameter :: z(NTST) = (/ &
       (624.2d0,-0.26123d0), &
       (-0.4d0,3.d0), &
       (0.6d0,2.d0), &
       (-1.d0,1.d0), &
       (-1.d0,-9.d0), &
       (-1.d0,9.d0), &
       (-0.0000000234545d0,1.1234d0), &
       (-3.d0,5.1d0), &
       (-53d0,30.1d0), &
       (0.0d0,0.12345d0), &
       (11d0,1d0), &
       (-22d0,-2d0), &
       (9d0,-28d0), &
       (21d0,-33d0), &
       (1d5,1d5), &
       (1d14,1d14), &
       (-3001d0,-1000d0), &
       (1d160,-1d159), &
       (-6.01d0,0.01d0), &
       (-0.7d0,-0.7d0), &
       (2.611780000000000d+01, 4.540909610972489d+03), &
       (0.8d7,0.3d7), &
       (-20d0,-19.8081d0), &
       (1d-16,-1.1d-16), &
       (2.3d-8,1.3d-8), &
       (6.3d0,-1d-13), &
       (6.3d0,1d-20), &
       (1d-20,6.3d0), &
       (1d-20,16.3d0), &
       (9d0,1d-300), &
       (6.01d0,0.11d0), &
       (8.01d0,1.01d-10), &
       (28.01d0,1d-300), &
       (10.01d0,1d-200), &
       (10.01d0,-1d-200), &
       (10.01d0,0.99d-10), &
       (10.01d0,-0.99d-10), &
       (1d-20,7.01d0), &
       (-1d0,7.01d0), &
       (5.99d0,7.01d0), &
       (1d0,0d0), &
       (55d0,0d0), &
       (-0.1d0,0d0), &
       (1d-20,0d0), &
       (0d0,5d-14), &
       (0d0,51d0) /)
  complex(dpc), parameter :: w(NTST) = (/ &
       ! w(z), computed with WolframAlpha
       ! ... note that WolframAlpha is problematic
       ! some of the above inputs, so I had to
       ! use the continued-fraction expansion
       ! in WolframAlpha in some cases, or switch to Maple
       (-3.78270245518980507452677445620103199303131110d-7,0.000903861276433172057331093754199933411710053155d0), &
       (0.1764906227004816847297495349730234591778719532788d0,-0.02146550539468457616788719893991501311573031095617d0), &
       (0.2410250715772692146133539023007113781272362309451d0,0.06087579663428089745895459735240964093522265589350d0), &
       (0.30474420525691259245713884106959496013413834051768d0,-0.20821893820283162728743734725471561394145872072738d0), &
       (7.317131068972378096865595229600561710140617977d34,8.321873499714402777186848353320412813066170427d34), &
       (0.0615698507236323685519612934241429530190806818395d0,-0.00676005783716575013073036218018565206070072304635d0), &
       (0.3960793007699874918961319170187598400134746631d0,-5.593152259116644920546186222529802777409274656d-9), &
       (0.08217199226739447943295069917990417630675021771804d0,-0.04701291087643609891018366143118110965272615832184d0), &
       (0.00457246000350281640952328010227885008541748668738d0,-0.00804900791411691821818731763401840373998654987934d0), &
       (0.8746342859608052666092782112565360755791467973338452d0,0.d0), &
       (0.00468190164965444174367477874864366058339647648741d0,0.0510735563901306197993676329845149741675029197050d0), &
       (-0.0023193175200187620902125853834909543869428763219d0,-0.025460054739731556004902057663500272721780776336d0), &
       (9.11463368405637174660562096516414499772662584d304,3.97101807145263333769664875189354358563218932d305), &
       (-4.4927207857715598976165541011143706155432296d281,-2.8019591213423077494444700357168707775769028d281), &
       (2.820947917809305132678577516325951485807107151d-6,2.820947917668257736791638444590253942253354058d-6), &
       (2.82094791773878143474039725787438662716372268d-15,2.82094791773878143474039725773333923127678361d-15), &
       (-0.0000563851289696244350147899376081488003110150498d0,-0.000169211755126812174631861529808288295454992688d0), &
       (-5.586035480670854326218608431294778077663867d-162,5.586035480670854326218608431294778077663867d-161), &
       (0.00016318325137140451888255634399123461580248456d0,-0.095232456573009287370728788146686162555021209999d0), &
       (0.69504753678406939989115375989939096800793577783885d0,-1.8916411171103639136680830887017670616339912024317d0), &
       (0.0001242418269653279656612334210746733213167234822d0,7.145975826320186888508563111992099992116786763d-7), &
       (2.318587329648353318615800865959225429377529825d-8,6.182899545728857485721417893323317843200933380d-8), &
       (-0.0133426877243506022053521927604277115767311800303d0,-0.0148087097143220769493341484176979826888871576145d0), &
       (1.00000000000000012412170838050638522857747934d0,1.12837916709551279389615890312156495593616433d-16), &
       (0.9999999853310704677583504063775310832036830015d0,2.595272024519678881897196435157270184030360773d-8), &
       (-1.4731421795638279504242963027196663601154624d-15,0.090727659684127365236479098488823462473074709d0), &
       (5.79246077884410284575834156425396800754409308d-18,0.0907276596841273652364790985059772809093822374d0), &
       (0.0884658993528521953466533278764830881245144368d0,1.37088352495749125283269718778582613192166760d-22), &
       (0.0345480845419190424370085249304184266813447878d0,2.11161102895179044968099038990446187626075258d-23), &
       (6.63967719958073440070225527042829242391918213d-36,0.0630820900592582863713653132559743161572639353d0), &
       (0.00179435233208702644891092397579091030658500743634d0,0.0951983814805270647939647438459699953990788064762d0), &
       (9.09760377102097999924241322094863528771095448d-13,0.0709979210725138550986782242355007611074966717d0), &
       (7.2049510279742166460047102593255688682910274423d-304,0.0201552956479526953866611812593266285000876784321d0), &
       (3.04543604652250734193622967873276113872279682d-44,0.0566481651760675042930042117726713294607499165d0), &
       (3.04543604652250734193622967873276113872279682d-44,0.0566481651760675042930042117726713294607499165d0), &
       (0.5659928732065273429286988428080855057102069081d-12,0.056648165176067504292998527162143030538756683302d0), &
       (-0.56599287320652734292869884280802459698927645d-12,0.0566481651760675042929985271621430305387566833029d0), &
       (0.0796884251721652215687859778119964009569455462d0,1.11474461817561675017794941973556302717225126d-22), &
       (0.07817195821247357458545539935996687005781943386550d0,-0.01093913670103576690766705513142246633056714279654d0), &
       (0.04670032980990449912809326141164730850466208439937d0,0.03944038961933534137558064191650437353429669886545d0), &
       (0.36787944117144232159552377016146086744581113103176d0,0.60715770584139372911503823580074492116122092866515d0), &
       (0.d0,0.010259688805536830986089913987516716056946786526145d0), &
       (0.99004983374916805357390597718003655777207908125383d0,-0.11208866436449538036721343053869621153527769495574d0), &
       (0.99999999999999999999999999999999999999990000d0,1.12837916709551257389615890312154517168802603d-20), &
       (0.999999999999943581041645226871305192054749891144158d0,0d0), &
       (0.0110604154853277201542582159216317923453996211744250d0,0d0) /)
  real(dp) :: errmax, err, re_err, im_err
  complex(dpc) :: fw
  integer :: i

  errmax = 0.d0
  do i = 1, NTST
     fw = Faddeeva_w(z(i),0.d0)
     err = abs((fw - w(i)) / w(i))
     re_err = abs(real(fw) - real(w(i)))
     if (real(w(i))/=0) re_err = re_err / abs(real(w(i)))
     im_err = abs(aimag(fw) - aimag(w(i)))
     if (aimag(w(i))/=0) im_err = im_err / abs(aimag(w(i)))
     write(6,'(a,99g24.15e3)') "z, w, wref, err, re_err, im_err = ", z(i), fw, w(i), err, re_err, im_err
     errmax = max(errmax,err)
     errmax = max(errmax,re_err)
     errmax = max(errmax,im_err)
  end do
  if (errmax > 1d-13) then
     write(6,'(a,g15.5,a)') "FAILURE -- relative error ", errmax, " too large!"
     STOP 1
  end if
  write(6,'(a,g15.5,a)') "SUCCESS (max relative error = ", errmax, ")"
end program main

#endif

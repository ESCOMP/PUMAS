!Note:  This code was taken directly from the shared code
!       library used by the NCAR Community Earth System Model (CESM).
!       A version can be found on Github here:
!       https://github.com/ESCOMP/CESM_share/blob/main/src/shr_spfn_mod.F90

! Define flags for compilers supporting Fortran 2008 intrinsics
! HAVE_GAMMA_INTRINSICS: gamma and log_gamma

! These compilers have the intrinsics.
! Intel also has them (and Cray), but as of mid-2015, our implementation is
! actually faster, in part because they do not properly vectorize, so we
! pretend that the compiler version doesn't exist.
#if defined CPRIBM || defined __GFORTRAN__
#define HAVE_GAMMA_INTRINSICS
#endif

! As of 5.3.1, NAG does not have any of these.

module pumas_gamma_function
  ! Module containg PUMAS gamma function

  use pumas_kinds, only: r8 => kind_r8

  implicit none
  private
  save

  real(r8), parameter :: pi = 3.1415926535897932384626434E0_r8

  ! Gamma functions
  ! Note that we lack an implementation of log_gamma, but we do have an
  ! implementation of the upper incomplete gamma function, which is not in
  ! Fortran 2008.

  ! Note also that this gamma function is only for double precision. We
  ! haven't needed an r4 version yet.

  public :: pumas_gamma

  interface pumas_gamma
     module procedure pumas_gamma_r8
  end interface pumas_gamma

  ! Mathematical constants
  ! sqrt(pi)
  real(r8), parameter :: sqrtpi = 1.77245385090551602729_r8

  ! Define machine-specific constants needed in this module.
  ! These were used by the original gamma and calerf functions to guarantee
  ! safety against overflow, and precision, on many different machines.

  ! By defining the constants in this way, we assume that 1/xmin is
  ! representable (i.e. does not overflow the real type). This assumption was
  ! not in the original code, but is valid for IEEE single and double
  ! precision.

  ! Double precision
  !---------------------------------------------------------------------
  ! Machine epsilon
  real(r8), parameter :: epsr8 = epsilon(1._r8)
  ! "Huge" value is returned when actual value would be infinite.
  real(r8), parameter :: xinfr8 = huge(1._r8)
  ! Smallest normal value.
  real(r8), parameter :: xminr8 = tiny(1._r8)
  ! Largest number that, when added to 1., yields 1.
  real(r8), parameter :: xsmallr8 = epsr8/2._r8
  ! Largest argument for which erfcx > 0.
  real(r8), parameter :: xmaxr8 = 1._r8/(sqrtpi*xminr8)

  ! For gamma
  ! Approximate value of largest acceptable argument to gamma,
  ! for IEEE double-precision.
  real(r8), parameter :: xbig_gamma = 171.624_r8

  !$acc declare copyin(xinfr8,epsr8,xminr8,xbig_gamma)

contains

  elemental function pumas_gamma_r8(x) result(res)
    !$acc routine seq
    real(r8), intent(in) :: x
    real(r8) :: res

#if defined HAVE_GAMMA_INTRINSICS
    ! Call intrinsic gamma.
    intrinsic gamma
    res = gamma(x)
#else
    ! No intrinsic
    res = pumas_gamma_nonintrinsic_r8(x)
#endif

  end function pumas_gamma_r8

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  pure function pumas_gamma_nonintrinsic_r8(X) result(gamma)
    !$acc routine seq

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    ! 7 Feb 2013 -- S. Santos
    ! The following comments are from the original version. Changes have
    ! been made to update syntax and allow inclusion into this module.
    !
    !----------------------------------------------------------------------
    !
    ! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
    !   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
    !   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
    !   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
    !   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
    !   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
    !   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
    !   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
    !   MACHINE-DEPENDENT CONSTANTS.
    !
    !
    !*******************************************************************
    !*******************************************************************
    !
    ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
    !
    ! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
    ! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
    ! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
    !          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
    !                  GAMMA(XBIG) = BETA**MAXEXP
    ! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
    !          APPROXIMATELY BETA**MAXEXP
    ! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
    !          1.0+EPS .GT. 1.0
    ! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
    !          1/XMININ IS MACHINE REPRESENTABLE
    !
    !     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
    !
    !                            BETA       MAXEXP        XBIG
    !
    ! CRAY-1         (S.P.)        2         8191        966.961
    ! CYBER 180/855
    !   UNDER NOS    (S.P.)        2         1070        177.803
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (S.P.)        2          128        35.040
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (D.P.)        2         1024        171.624
    ! IBM 3033       (D.P.)       16           63        57.574
    ! VAX D-FORMAT   (D.P.)        2          127        34.844
    ! VAX G-FORMAT   (D.P.)        2         1023        171.489
    !
    !                            XINF         EPS        XMININ
    !
    ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
    ! CYBER 180/855
    !   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
    ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
    ! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
    ! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
    !
    !*******************************************************************
    !*******************************************************************
    !
    ! ERROR RETURNS
    !
    !  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
    !     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
    !     TO BE FREE OF UNDERFLOW AND OVERFLOW.
    !
    !
    !  INTRINSIC FUNCTIONS REQUIRED ARE:
    !
    !     INT, DBLE, EXP, LOG, REAL, SIN
    !
    !
    ! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
    !              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
    !              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
    !              (ED.), SPRINGER VERLAG, BERLIN, 1976.
    !
    !              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
    !              SONS, NEW YORK, 1968.
    !
    !  LATEST MODIFICATION: OCTOBER 12, 1989
    !
    !  AUTHORS: W. J. CODY AND L. STOLTZ
    !           APPLIED MATHEMATICS DIVISION
    !           ARGONNE NATIONAL LABORATORY
    !           ARGONNE, IL 60439
    !
    !----------------------------------------------------------------------

    real(r8), intent(in) :: x
    real(r8) :: gamma
    real(r8) :: fact, res, sum, xden, xnum, y, y1, ysq, z

    integer :: i, n
    logical :: negative_odd

    ! log(2*pi)/2
    real(r8), parameter :: logsqrt2pi = 0.9189385332046727417803297E0_r8

    !----------------------------------------------------------------------
    !  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
    !     APPROXIMATION OVER (1,2).
    !----------------------------------------------------------------------
    real(r8), parameter :: P(8) = &
         (/-1.71618513886549492533811E+0_r8, 2.47656508055759199108314E+1_r8, &
         -3.79804256470945635097577E+2_r8, 6.29331155312818442661052E+2_r8, &
         8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8, &
         -3.61444134186911729807069E+4_r8, 6.64561438202405440627855E+4_r8 /)
    real(r8), parameter :: Q(8) = &
         (/-3.08402300119738975254353E+1_r8, 3.15350626979604161529144E+2_r8, &
         -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8, &
         2.25381184209801510330112E+4_r8, 4.75584627752788110767815E+3_r8, &
         -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8 /)
    !----------------------------------------------------------------------
    !  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
    !----------------------------------------------------------------------
    real(r8), parameter :: C(7) = &
         (/-1.910444077728E-03_r8,          8.4171387781295E-04_r8, &
         -5.952379913043012E-04_r8,       7.93650793500350248E-04_r8, &
         -2.777777777777681622553E-03_r8, 8.333333333333333331554247E-02_r8, &
         5.7083835261E-03_r8 /)

    negative_odd = .false.
    fact = 1._r8
    n = 0
    y = x
    if (y <= 0._r8) then
       !----------------------------------------------------------------------
       !  ARGUMENT IS NEGATIVE
       !----------------------------------------------------------------------
       y = -x
       y1 = aint(y)
       res = y - y1
       if (res /= 0._r8) then
          negative_odd = (y1 /= aint(y1*0.5_r8)*2._r8)
          fact = -pi/sin(pi*res)
          y = y + 1._r8
       else
          gamma = xinfr8
          return
       end if
    end if
    !----------------------------------------------------------------------
    !  ARGUMENT IS POSITIVE
    !----------------------------------------------------------------------
    if (y < epsr8) then
       !----------------------------------------------------------------------
       !  ARGUMENT .LT. EPS
       !----------------------------------------------------------------------
       if (y >= xminr8) then
          res = 1._r8/y
       else
          gamma = xinfr8
          return
       end if
    elseif (y < 12._r8) then
       y1 = y
       if (y < 1._r8) then
          !----------------------------------------------------------------------
          !  0.0 .LT. ARGUMENT .LT. 1.0
          !----------------------------------------------------------------------
          z = y
          y = y + 1._r8
       else
          !----------------------------------------------------------------------
          !  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
          !----------------------------------------------------------------------
          n = int(y) - 1
          y = y - real(n, r8)
          z = y - 1._r8
       end if
       !----------------------------------------------------------------------
       !  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
       !----------------------------------------------------------------------
       xnum = 0._r8
       xden = 1._r8
       do i=1,8
          xnum = (xnum+P(i))*z
          xden = xden*z + Q(i)
       end do
       res = xnum/xden + 1._r8
       if (y1 < y) then
          !----------------------------------------------------------------------
          !  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
          !----------------------------------------------------------------------
          res = res/y1
       elseif (y1 > y) then
          !----------------------------------------------------------------------
          !  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
          !----------------------------------------------------------------------
          do i = 1,n
             res = res*y
             y = y + 1._r8
          end do
       end if
    else
       !----------------------------------------------------------------------
       !  EVALUATE FOR ARGUMENT .GE. 12.0,
       !----------------------------------------------------------------------
       if (y <= xbig_gamma) then
          ysq = y*y
          sum = C(7)
          do i=1,6
             sum = sum/ysq + C(i)
          end do
          sum = sum/y - y + logsqrt2pi
          sum = sum + (y-0.5_r8)*log(y)
          res = exp(sum)
       else
          gamma = xinfr8
          return
       end if
    end if
    !----------------------------------------------------------------------
    !  FINAL ADJUSTMENTS AND RETURN
    !----------------------------------------------------------------------
    if (negative_odd)  res = -res
    if (fact /= 1._r8) res = fact/res
    gamma = res
    ! ---------- LAST LINE OF GAMMA ----------
  end function pumas_gamma_nonintrinsic_r8

end module pumas_gamma_function

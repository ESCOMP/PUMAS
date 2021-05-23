module micro_mg_utils

!--------------------------------------------------------------------------
!
! This module contains process rates and utility functions used by the MG
! microphysics.
!
! Original MG authors: Andrew Gettelman, Hugh Morrison
! Contributions from: Peter Caldwell, Xiaohong Liu and Steve Ghan
!
! Separated from MG 1.5 by B. Eaton.
! Separated module switched to MG 2.0 and further changes by S. Santos.
!
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!
!--------------------------------------------------------------------------
!
! List of required external functions that must be supplied:
!   gamma --> standard mathematical gamma function (if gamma is an
!       intrinsic, define HAVE_GAMMA_INTRINSICS)
!
!--------------------------------------------------------------------------
!
! Constants that must be specified in the "init" method (module variables):
!
! kind            kind of reals (to verify correct linkage only) -
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                   J kg-1 K-1
! rh2o            gas constant for water vapor                   J kg-1 K-1
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! tmelt           temperature of melting point for water         K
! latvap          latent heat of vaporization                    J kg-1
! latice          latent heat of fusion                          J kg-1
!
!--------------------------------------------------------------------------

#ifndef HAVE_GAMMA_INTRINSICS
use shr_spfn_mod, only: gamma => shr_spfn_gamma
#endif

implicit none
private
save

public :: &
     micro_mg_utils_init, &
     size_dist_param_liq, &
     size_dist_param_basic, &
     avg_diameter, &
     rising_factorial, &
     ice_deposition_sublimation, &
     ice_deposition_sublimation_mg4, &  !! ktc 
     sb2001v2_liq_autoconversion,&
     sb2001v2_accre_cld_water_rain,&       
     kk2000_liq_autoconversion, &
     ice_autoconversion, &
     immersion_freezing, &
     contact_freezing, &
     snow_self_aggregation, &
     ice_self_aggregation, &  !! mg4 added for using instead of snow_self_agg
     accrete_cloud_water_snow, &
     accrete_cloud_water_ice, &  !! mg4 added for using instead of cloud_water_snow
     secondary_ice_production, &
     accrete_rain_snow, &
     accrete_rain_ice, &  !! mg4 added for using instead of rain_snow
     heterogeneous_rain_freezing, &
     accrete_cloud_water_rain, &
     self_collection_rain, &
     accrete_cloud_ice_snow, &
     evaporate_sublimate_precip, &
     bergeron_process_snow, &
     graupel_collecting_snow, &
     graupel_collecting_rain, &
     graupel_collecting_cld_water, &
     graupel_riming_liquid_snow, &
     graupel_rain_riming_snow, &
     graupel_rime_splintering, &
     evaporate_sublimate_precip_graupel, &
     evaporate_sublimate_precip_mg4, &
     evaporate_sublimate_precip_graupel_mg4, &
     access_lookup_table, & !! mg4
     access_lookup_table_coll, & !! mg4
     init_lookup_table, &      !! mg4
     avg_diameter_vec

! 8 byte real and integer
integer, parameter, public :: r8 = selected_real_kind(12)
integer, parameter, public :: i8 = selected_int_kind(18)
integer, parameter         :: VLENS = 128  ! vector length of a GPU compute kernel

public :: MGHydrometeorProps

type :: MGHydrometeorProps
   ! Density (kg/m^3)
   real(r8) :: rho
   ! Information for size calculations.
   ! Basic calculation of mean size is:
   !     lambda = (shape_coef*nic/qic)^(1/eff_dim)
   ! Then lambda is constrained by bounds.
   real(r8) :: eff_dim
   real(r8) :: shape_coef
   real(r8) :: lambda_bounds(2)
   ! Minimum average particle mass (kg).
   ! Limit is applied at the beginning of the size distribution calculations.
   real(r8) :: min_mean_mass
end type MGHydrometeorProps

interface MGHydrometeorProps
   module procedure NewMGHydrometeorProps
end interface

type(MGHydrometeorProps), public :: mg_liq_props
type(MGHydrometeorProps), public :: mg_ice_props
type(MGHydrometeorProps), public :: mg_rain_props
type(MGHydrometeorProps), public :: mg_snow_props
type(MGHydrometeorProps), public :: mg_graupel_props
type(MGHydrometeorProps), public :: mg_hail_props

interface size_dist_param_liq
  module procedure size_dist_param_liq_2D
  module procedure size_dist_param_liq_vect
  module procedure size_dist_param_liq_line
end interface
interface size_dist_param_basic
  module procedure size_dist_param_basic_2D
  module procedure size_dist_param_basic_vect
  module procedure size_dist_param_basic_line
end interface
interface calc_ab
  module procedure calc_ab_line
  module procedure calc_ab_vect
end interface

!=================================================
! Public module parameters (mostly for MG itself)
!=================================================

! Pi to 20 digits; more than enough to reach the limit of double precision.
real(r8), parameter, public :: pi = 3.14159265358979323846_r8

! "One minus small number": number near unity for round-off issues.
real(r8), parameter, public :: omsm   = 1._r8 - 1.e-5_r8

! Smallest mixing ratio considered in microphysics.
real(r8), parameter, public :: qsmall = 1.e-18_r8

! minimum allowed cloud fraction
real(r8), parameter, public :: mincld = 0.0001_r8

real(r8), parameter, public :: rhosn = 250._r8  ! bulk density snow
real(r8), parameter, public :: rhoi = 500._r8   ! bulk density ice
real(r8), parameter, public :: rhow = 1000._r8  ! bulk density liquid
real(r8), parameter, public :: rhows = 917._r8  ! bulk density water solid

!Hail and Graupel (set in MG3)
real(r8), parameter, public :: rhog = 500._r8 
real(r8), parameter, public :: rhoh = 900._r8 

! fall speed parameters, V = aD^b (V is in m/s)
! droplets
real(r8), parameter, public :: ac = 3.e7_r8
real(r8), parameter, public :: bc = 2._r8
! snow
real(r8), parameter, public :: as = 11.72_r8
real(r8), parameter, public :: bs = 0.41_r8
! cloud ice
real(r8), parameter, public :: ai = 700._r8
real(r8), parameter, public :: bi = 1._r8
! small cloud ice (r< 10 um) - sphere, bulk density
real(r8), parameter, public :: aj = ac*((rhoi/rhows)**(bc/3._r8))*rhows/rhow
real(r8), parameter, public :: bj = bc
! rain
real(r8), parameter, public :: ar = 841.99667_r8
real(r8), parameter, public :: br = 0.8_r8
! graupel
real(r8), parameter, public :: ag = 19.3_r8
real(r8), parameter, public :: bg = 0.37_r8
! hail
real(r8), parameter, public :: ah = 114.5_r8 
real(r8), parameter, public :: bh = 0.5_r8

! mass of new crystal due to aerosol freezing and growth (kg)
! Make this consistent with the lower bound, to support UTLS and
! stratospheric ice, and the smaller ice size limit.
real(r8), parameter, public :: mi0 = 4._r8/3._r8*pi*rhoi*(1.e-6_r8)**3

! mass of new graupel particle  (assume same as mi0 for now, may want to make bigger?)
!real(r8), parameter, public :: mg0 = 4._r8/3._r8*pi*rhoi*(1.e-6_r8)**3
!or set based on M2005:
real(r8), parameter, public :: mg0 = 1.6e-10_r8
! radius of contact nuclei
real(r8), parameter, public :: mmult = 4._r8/3._r8*pi*rhoi*(5.e-6_r8)**3

!=================================================
! Private module parameters
!=================================================

! Signaling NaN bit pattern that represents a limiter that's turned off.
integer(i8), parameter :: limiter_off = int(Z'7FF1111111111111', i8)

! alternate threshold used for some in-cloud mmr
real(r8), parameter, public :: icsmall = 1.e-8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d
! exponent
real(r8), parameter :: dsph = 3._r8

! Bounds for mean diameter for different constituents.
real(r8), parameter :: lam_bnd_rain(2) = 1._r8/[500.e-6_r8, 20.e-6_r8]
real(r8), parameter :: lam_bnd_snow(2) = 1._r8/[2000.e-6_r8, 10.e-6_r8]

! Minimum average mass of particles.
real(r8), parameter :: min_mean_mass_liq = 1.e-20_r8
real(r8), parameter :: min_mean_mass_ice = 1.e-20_r8

! ventilation parameters
! for snow
real(r8), parameter :: f1s = 0.86_r8
real(r8), parameter :: f2s = 0.28_r8
! for rain
real(r8), parameter :: f1r = 0.78_r8
real(r8), parameter :: f2r = 0.308_r8

! collection efficiencies
! aggregation of cloud ice and snow
real(r8), parameter :: eii = 0.5_r8
! collection efficiency, ice-droplet collisions
real(r8), parameter, public :: ecid = 0.7_r8
! collection efficiency between droplets/rain and snow/rain
real(r8), parameter, public :: ecr = 1.0_r8

! immersion freezing parameters, bigg 1953
real(r8), parameter :: bimm = 100._r8
real(r8), parameter :: aimm = 0.66_r8

! Mass of each raindrop created from autoconversion.
real(r8), parameter :: droplet_mass_25um = 4._r8/3._r8*pi*rhow*(25.e-6_r8)**3
real(r8), parameter :: droplet_mass_40um = 4._r8/3._r8*pi*rhow*(40.e-6_r8)**3

!!!!!!
! MG4 look up table parameters
!!!!!!
! ice microphysics lookup table array dimensions
integer, parameter, public :: tsize     = 3
integer, parameter, public :: isize     = 20
integer, parameter, public :: jsize     = 20
integer, parameter, public :: rcollsize = 30

! number of ice microphysical quantities used from lookup table
integer, parameter, public :: tabsize   = 14

! number of ice-rain collection microphysical quantities used from lookup table
integer, parameter, public :: colltabsize = 2


!ice lookup table values
real(r8) :: itab(tsize,isize,jsize,tabsize)

!ice lookup table values for ice-rain collision/collection
real(r8) :: itabcoll(tsize,isize,jsize,rcollsize,colltabsize)

private :: itab,itabcoll   

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_mg_init
real(r8) :: rv          ! water vapor gas constant
real(r8) :: cpp         ! specific heat of dry air
real(r8) :: tmelt       ! freezing point of water (K)

real(r8) :: ra        ! dry air gas constant

! latent heats of:
real(r8) :: xxlv        ! vaporization
real(r8) :: xlf         ! freezing
real(r8) :: xxls        ! sublimation

! additional constants to help speed up code
real(r8) :: gamma_bs_plus3
real(r8) :: gamma_half_br_plus5
real(r8) :: gamma_half_bs_plus5
real(r8) :: gamma_2bs_plus2

!$acc declare create (rv,cpp,tmelt,xxlv,xxls,gamma_bs_plus3,   &
!$acc                 gamma_half_br_plus5,gamma_half_bs_plus5, &
!$acc                 gamma_2bs_plus2)

!=========================================================
! Utilities that are cheaper if the compiler knows that
! some argument is an integer.
!=========================================================

interface rising_factorial
   module procedure rising_factorial_r8
   module procedure rising_factorial_r8_vec
   module procedure rising_factorial_integer
   module procedure rising_factorial_integer_vec
end interface rising_factorial

interface var_coef
   module procedure var_coef_r8
   module procedure var_coef_r8_vect
   module procedure var_coef_integer
   module procedure var_coef_integer_vect
end interface var_coef

!==========================================================================
contains
!==========================================================================

! Initialize module variables.
!
! "kind" serves no purpose here except to check for unlikely linking
! issues; always pass in the kind for a double precision real.
!
! "errstring" is the only output; it is blank if there is no error, or set
! to a message if there is an error.
!
! Check the list at the top of this module for descriptions of all other
! arguments.
subroutine micro_mg_utils_init( kind, rair, rh2o, cpair, tmelt_in, latvap, &
     latice, dcs, errstring)

  integer,  intent(in)  :: kind
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: tmelt_in
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice
  real(r8), intent(in)  :: dcs

  character(128), intent(out) :: errstring

  ! Name this array to workaround an XLF bug (otherwise could just use the
  ! expression that sets it).
  real(r8) :: ice_lambda_bounds(2)

  !-----------------------------------------------------------------------

  errstring = ' '

  if( kind .ne. r8 ) then
     errstring = 'micro_mg_init: KIND of reals does not match'
     return
  endif

  ! declarations for MG code (transforms variable names)

  rv    = rh2o          ! water vapor gas constant
  cpp   = cpair         ! specific heat of dry air
  tmelt = tmelt_in
  ra    = rair          ! dry air gas constant

  ! latent heats

  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

  ! Define constants to help speed up code (this limits calls to gamma function)
  gamma_bs_plus3      = gamma(3._r8+bs)
  gamma_half_br_plus5 = gamma(5._r8/2._r8+br/2._r8)
  gamma_half_bs_plus5 = gamma(5._r8/2._r8+bs/2._r8)
  gamma_2bs_plus2     = gamma(2._r8*bs+2._r8)  

  ! Don't specify lambda bounds for cloud liquid, as they are determined by
  ! pgam dynamically.
  mg_liq_props = MGHydrometeorProps(rhow, dsph, &
       min_mean_mass=min_mean_mass_liq)

  ! Mean ice diameter can not grow bigger than twice the autoconversion
  ! threshold for snow.
  ice_lambda_bounds = 1._r8/[2._r8*dcs, 1.e-6_r8]

  mg_ice_props = MGHydrometeorProps(rhoi, dsph, &
       ice_lambda_bounds, min_mean_mass_ice)

  mg_rain_props = MGHydrometeorProps(rhow, dsph, lam_bnd_rain)
  mg_snow_props = MGHydrometeorProps(rhosn, dsph, lam_bnd_snow)
  mg_graupel_props = MGHydrometeorProps(rhog, dsph, lam_bnd_snow)
  mg_hail_props = MGHydrometeorProps(rhoh, dsph, lam_bnd_snow)

  !$acc update device (rv,cpp,tmelt,xxlv,xxls,gamma_bs_plus3,   &
  !$acc                gamma_half_br_plus5,gamma_half_bs_plus5, &
  !$acc                gamma_2bs_plus2)

end subroutine micro_mg_utils_init

! Constructor for a constituent property object.
function NewMGHydrometeorProps(rho, eff_dim, lambda_bounds, min_mean_mass) &
     result(res)
  real(r8), intent(in) :: rho, eff_dim
  real(r8), intent(in), optional :: lambda_bounds(2), min_mean_mass
  type(MGHydrometeorProps) :: res

  res%rho = rho
  res%eff_dim = eff_dim
  if (present(lambda_bounds)) then
     res%lambda_bounds = lambda_bounds
  else
     res%lambda_bounds = no_limiter()
  end if
  if (present(min_mean_mass)) then
     res%min_mean_mass = min_mean_mass
  else
     res%min_mean_mass = no_limiter()
  end if

  res%shape_coef = rho*pi*gamma(eff_dim+1._r8)/6._r8

end function NewMGHydrometeorProps

!========================================================================
!FORMULAS
!========================================================================

! Use gamma function to implement rising factorial extended to the reals.
subroutine rising_factorial_r8(x, n, res)
  !$acc routine seq
  real(r8), intent(in) :: x, n
  real(r8), intent(out) :: res

  !$acc data present (x,res)

  res = gamma(x+n)/gamma(x)
  
  !$acc end data
end subroutine rising_factorial_r8

subroutine rising_factorial_r8_vec(x, n, res,vlen)
  integer, intent(in)   :: vlen
  real(r8), intent(in)  :: x(vlen), n
  real(r8), intent(out) :: res(vlen)
  integer :: i
  real(r8) :: tmp

  !$acc data present (x,res)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(tmp)
  do i=1,vlen
     tmp = x(i)+n
     res(i) = gamma(tmp)
     tmp = gamma(x(i))
     res(i) = res(i)/tmp
  end do
  !$acc end parallel

  !$acc end data

end subroutine rising_factorial_r8_vec

! Rising factorial can be performed much cheaper if n is a small integer.
subroutine rising_factorial_integer(x, n, res)
  !$acc routine seq
  real(r8), intent(in)  :: x
  integer,  intent(in)  :: n
  real(r8), intent(out) :: res

  integer :: i
  real(r8) :: factor

  !$acc data present (x,res)

  res = 1._r8
  factor = x

  !$acc loop seq
  do i = 1, n
     res = res * factor
     factor = factor + 1._r8
  end do

  !$acc end data
end subroutine rising_factorial_integer

subroutine rising_factorial_integer_vec(x, n, res,vlen)
  integer,  intent(in)  :: vlen
  real(r8), intent(in)  :: x(vlen)
  integer,  intent(in)  :: n
  real(r8), intent(out) :: res(vlen)

  integer  :: i,j
  real(r8) :: factor(vlen)

  !$acc data present (x,res) &
  !$acc      create  (factor)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     res(i)    = 1._r8
     factor(i) = x(i)
  end do

  if (n == 3) then
    !$acc loop gang vector
    do i=1,vlen
       res(i)    = res(i) * factor(i)
       factor(i) = factor(i) + 1._r8
       res(i)    = res(i) * factor(i)
       factor(i) = factor(i) + 1._r8
       res(i)    = res(i) * factor(i)
    end do
  elseif (n == 2) then
    !$acc loop gang vector
    do i=1,vlen
       res(i)    = res(i) * factor(i)
       factor(i) = factor(i) + 1._r8
       res(i)    = res(i) * factor(i)
    end do
  else
    !$acc loop seq
    do j = 1, n
       !$acc loop gang vector
       do i = 1, vlen
          res(i)    = res(i) * factor(i)
          factor(i) = factor(i) + 1._r8
       end do
    end do
  end if
  !$acc end parallel

  !$acc end data
end subroutine rising_factorial_integer_vec

! Calculate correction due to latent heat for evaporation/sublimation
subroutine calc_ab_line(t, qv, xxl, ab)
  !$acc routine seq
  real(r8), intent(in)  :: t     ! Temperature
  real(r8), intent(in)  :: qv    ! Saturation vapor pressure
  real(r8), intent(in)  :: xxl   ! Latent heat
  real(r8), intent(out) :: ab

  real(r8) :: dqsdt

  !$acc data present (t,qv,xxl,ab)

  dqsdt = xxl*qv / (rv * t**2)
  ab = 1._r8 + dqsdt*xxl/cpp

  !$acc end data
end subroutine calc_ab_line

! Calculate correction due to latent heat for evaporation/sublimation
subroutine calc_ab_vect(t, qv, xxl, ab, vlen)
  integer,  intent(in) :: vlen
  real(r8), intent(in) :: t(vlen)     ! Temperature
  real(r8), intent(in) :: qv(vlen)    ! Saturation vapor pressure
  real(r8), intent(in) :: xxl         ! Latent heat

  real(r8), intent(out) :: ab(vlen)
  real(r8) :: dqsdt
  integer :: i

  !$acc data present (t,qv,xxl,ab)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(dqsdt)
  do i=1,vlen
     dqsdt = xxl*qv(i) / (rv * t(i)**2)
     ab(i) = 1._r8 + dqsdt*xxl/cpp
  end do
  !$acc end parallel

  !$acc end data
end subroutine calc_ab_vect

! get cloud droplet size distribution parameters
subroutine size_dist_param_liq_line(props, qcic, ncic, rho, pgam, lamc)
  !$acc routine seq
  type(MGHydrometeorProps), intent(in)    :: props
  real(r8),                 intent(in)    :: qcic
  real(r8),                 intent(inout) :: ncic
  real(r8),                 intent(in)    :: rho
  real(r8),                 intent(out)   :: pgam
  real(r8),                 intent(out)   :: lamc

  ! local variables
  type(MGHydrometeorProps)                :: props_loc
  real(r8)                                :: tmp

  !$acc data present (props,qcic,ncic,rho,pgam,lamc)

  if (qcic > qsmall) then
     ! Local copy of properties that can be modified.
     ! (Elemental routines that operate on arrays can't modify scalar
     ! arguments.)
     props_loc = props

     ! Get pgam from fit Rotstayn and Liu 2003 (changed from Martin 1994 for CAM6)
     pgam = 1.0_r8 - 0.7_r8 * exp(-0.008_r8*1.e-6_r8*ncic*rho)
     pgam = 1._r8/(pgam**2) - 1._r8
     pgam = max(pgam, 2._r8)

     ! Set coefficient for use in size_dist_param_basic.
     ! The 3D case is so common and optimizable that we specialize it:
     if (props_loc%eff_dim == 3._r8) then
        call rising_factorial(pgam+1._r8, 3, tmp)
        props_loc%shape_coef = pi / 6._r8 * props_loc%rho * tmp
     else
        call rising_factorial(pgam+1._r8, props_loc%eff_dim, tmp)
        props_loc%shape_coef = pi / 6._r8 * props_loc%rho * tmp
     end if

     ! Limit to between 2 and 50 microns mean size.
     props_loc%lambda_bounds = (pgam+1._r8)*1._r8/[50.e-6_r8, 2.e-6_r8]

     call size_dist_param_basic(props_loc, qcic, ncic, lamc)
  else
     ! pgam not calculated in this case, so set it to a value likely to
     ! cause an error if it is accidentally used
     ! (gamma function undefined for negative integers)
     pgam = -100._r8
     lamc = 0._r8
  end if

  !$acc end data
end subroutine size_dist_param_liq_line

! get cloud droplet size distribution parameters

subroutine size_dist_param_liq_2D(props, qcic, ncic, rho, pgam, lamc, dim1, dim2)

  type(mghydrometeorprops),       intent(in)    :: props
  integer,                        intent(in)    :: dim1, dim2
  real(r8), dimension(dim1,dim2), intent(in)    :: qcic
  real(r8), dimension(dim1,dim2), intent(inout) :: ncic
  real(r8), dimension(dim1,dim2), intent(in)    :: rho
  real(r8), dimension(dim1,dim2), intent(out)   :: pgam
  real(r8), dimension(dim1,dim2), intent(out)   :: lamc

  ! local variables
  integer  :: i, k
  real(r8) :: tmp(dim1,dim2),pgamp1(dim1,dim2)
  real(r8) :: shapeC(dim1,dim2),lbnd(dim1,dim2),ubnd(dim1,dim2)

  !$acc data present (props,qcic,ncic,rho,pgam,lamc) &
  !$acc      create  (tmp,pgamp1,shapeC,lbnd,ubnd)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k = 1, dim2
     do i = 1, dim1
        if (qcic(i,k) > qsmall) then
           ! Get pgam from fit Rotstayn and Liu 2003 (changed from Martin 1994 for CAM6)
           pgam(i,k) = 1.0_r8 - 0.7_r8 * exp(-0.008_r8*1.e-6_r8*ncic(i,k)*rho(i,k))
           pgam(i,k) = 1._r8/(pgam(i,k)**2) - 1._r8
           pgam(i,k) = max(pgam(i,k), 2._r8)
           pgamp1(i,k) = pgam(i,k)+1._r8
        else
           pgamp1(i,k) = 0._r8
        end if
     end do
  end do
  !$acc end parallel

  ! Set coefficient for use in size_dist_param_basic.
  ! The 3D case is so common and optimizable that we specialize it:
  if (props%eff_dim == 3._r8) then
     call rising_factorial_integer_vec(pgamp1,3,tmp,dim1*dim2)
  else
     call rising_factorial_r8_vec(pgamp1, props%eff_dim,tmp,dim1*dim2)
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k = 1, dim2
     do i = 1, dim1
        if (qcic(i,k) > qsmall) then
           shapeC(i,k) = pi / 6._r8 * props%rho * tmp(i,k)
           ! Limit to between 2 and 50 microns mean size.
           lbnd(i,k)   = pgamp1(i,k)*1._r8/50.e-6_r8
           ubnd(i,k)   = pgamp1(i,k)*1._r8/2.e-6_r8
        else
           shapeC(i,k) = 0._r8
           lbnd(i,k) = 0._r8
           ubnd(i,k) = 0._r8
        end if
     end do
  end do
  !$acc end parallel 

  call size_dist_param_basic_vect2(props, qcic, ncic, shapeC, lbnd, ubnd, lamc, dim1*dim2)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k = 1, dim2
     do i = 1, dim1
        if (qcic(i,k) <= qsmall) then
           ! pgam not calculated in this case, so set it to a value likely to
           ! cause an error if it is accidentally used
           ! (gamma function undefined for negative integers)
           pgam(i,k) = -100._r8
           lamc(i,k) = 0._r8
        end if
     end do
  end do
  !$acc end parallel

  !$acc end data
end subroutine size_dist_param_liq_2D

! get cloud droplet size distribution parameters

subroutine size_dist_param_liq_vect(props, qcic, ncic, rho, pgam, lamc, vlen)

  type(mghydrometeorprops),  intent(in)    :: props
  integer,                   intent(in)    :: vlen 
  real(r8), dimension(vlen), intent(in)    :: qcic
  real(r8), dimension(vlen), intent(inout) :: ncic
  real(r8), dimension(vlen), intent(in)    :: rho
  real(r8), dimension(vlen), intent(out)   :: pgam
  real(r8), dimension(vlen), intent(out)   :: lamc

  ! local variables
  integer  :: i
  real(r8) :: tmp(vlen),pgamp1(vlen)
  real(r8) :: shapeC(vlen),lbnd(vlen),ubnd(vlen)

  !$acc data present (props,qcic,ncic,rho,pgam,lamc) &
  !$acc      create  (tmp,pgamp1,shapeC,lbnd,ubnd)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1, vlen
     if (qcic(i) > qsmall) then
        ! Get pgam from fit Rotstayn and Liu 2003 (changed from Martin 1994 for CAM6)
        pgam(i) = 1.0_r8 - 0.7_r8 * exp(-0.008_r8*1.e-6_r8*ncic(i)*rho(i))
        pgam(i) = 1._r8/(pgam(i)**2) - 1._r8
        pgam(i) = max(pgam(i), 2._r8)
        pgamp1(i) = pgam(i)+1._r8
     else
        pgamp1(i) = 0._r8   
     end if
  end do
  !$acc end parallel

  ! Set coefficient for use in size_dist_param_basic.
  ! The 3D case is so common and optimizable that we specialize it:
  if (props%eff_dim == 3._r8) then
     call rising_factorial_integer_vec(pgamp1,3,tmp,vlen)
  else
     call rising_factorial_r8_vec(pgamp1, props%eff_dim,tmp,vlen)
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1, vlen 
     if (qcic(i) > qsmall) then
        shapeC(i) = pi / 6._r8 * props%rho * tmp(i)
        ! Limit to between 2 and 50 microns mean size.
        lbnd(i)   = pgamp1(i)*1._r8/50.e-6_r8
        ubnd(i)   = pgamp1(i)*1._r8/2.e-6_r8
     else
        shapeC(i) = 0._r8
        lbnd(i) = 0._r8
        ubnd(i) = 0._r8
     end if
  end do
  !$acc end parallel 

  call size_dist_param_basic_vect2(props, qcic, ncic, shapeC, lbnd, ubnd, lamc, vlen)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1, vlen
     if (qcic(i) <= qsmall) then
        ! pgam not calculated in this case, so set it to a value likely to
        ! cause an error if it is accidentally used
        ! (gamma function undefined for negative integers)
        pgam(i) = -100._r8
        lamc(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine size_dist_param_liq_vect

! Basic routine for getting size distribution parameters.
subroutine size_dist_param_basic_line(props, qic, nic, lam, n0)
  !$acc routine seq
  type(MGHydrometeorProps), intent(in)    :: props
  real(r8),                 intent(in)    :: qic
  real(r8),                 intent(inout) :: nic
  real(r8),                 intent(out)           :: lam
  real(r8),                 intent(out), optional :: n0
  
  logical :: present_n0 
  present_n0 = present(n0)

  !$acc data present (props,qic,nic,lam,n0)

  if (qic > qsmall) then
     ! add upper limit to in-cloud number concentration to prevent
     ! numerical error
     if (limiter_is_on(props%min_mean_mass)) then
        nic = min(nic, qic / props%min_mean_mass)
     end if

     ! lambda = (c n/q)^(1/d)
     lam = (props%shape_coef * nic/qic)**(1._r8/props%eff_dim)

     ! check for slope
     ! adjust vars
     if (lam < props%lambda_bounds(1)) then
        lam = props%lambda_bounds(1)
        nic = lam**(props%eff_dim) * qic/props%shape_coef
     else if (lam > props%lambda_bounds(2)) then
        lam = props%lambda_bounds(2)
        nic = lam**(props%eff_dim) * qic/props%shape_coef
     end if
  else
     lam = 0._r8
  end if

  if (present_n0) n0 = nic * lam

  !$acc end data
end subroutine size_dist_param_basic_line

subroutine size_dist_param_basic_vect(props, qic, nic, lam, vlen, n0)

  type (mghydrometeorprops), intent(in)    :: props
  integer,                   intent(in)    :: vlen
  real(r8), dimension(vlen), intent(in)    :: qic
  real(r8), dimension(vlen), intent(inout) :: nic
  real(r8), dimension(vlen), intent(out)   :: lam
  real(r8), dimension(vlen), intent(out), optional :: n0
  integer  :: i
  logical  :: limiterActive, present_n0
  real(r8) :: effDim,shapeCoef,ubnd,lbnd, minMass

  !$acc data present (props,qic,nic,lam,n0)

  limiterActive = limiter_is_on(props%min_mean_mass)
  effDim    = props%eff_dim
  shapeCoef = props%shape_coef
  lbnd      = props%lambda_bounds(1)
  ubnd      = props%lambda_bounds(2)
  minMass   = props%min_mean_mass
  present_n0 = present(n0)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1, vlen
     if (qic(i) > qsmall) then
        ! add upper limit to in-cloud number concentration to prevent
        ! numerical error
        if (limiterActive) then
           nic(i) = min(nic(i), qic(i) / minMass)
        end if

        ! lambda = (c n/q)^(1/d)
        lam(i) = (shapeCoef * nic(i)/qic(i))**(1._r8/effDim)

        ! check for slope
        ! adjust vars
        if (lam(i) < lbnd) then
           lam(i) = lbnd
           nic(i) = lam(i)**(effDim) * qic(i)/shapeCoef
        else if (lam(i) > ubnd) then
           lam(i) = ubnd
           nic(i) = lam(i)**(effDim) * qic(i)/shapeCoef
        end if

     else
        lam(i) = 0._r8
     end if
  end do

  if (present_n0) then
     !$acc loop gang vector
     do i = 1, vlen 
        n0(i) = nic(i) * lam(i)
     end do
  end if
  !$acc end parallel

  !$acc end data
end subroutine size_dist_param_basic_vect

subroutine size_dist_param_basic_2D(props, qic, nic, lam, dim1, dim2, n0)

  type (mghydrometeorprops),      intent(in)    :: props
  integer,                        intent(in)    :: dim1, dim2
  real(r8), dimension(dim1,dim2), intent(in)    :: qic
  real(r8), dimension(dim1,dim2), intent(inout) :: nic
  real(r8), dimension(dim1,dim2), intent(out)   :: lam
  real(r8), dimension(dim1,dim2), intent(out), optional :: n0
  integer  :: i, k 
  logical  :: limiterActive, present_n0
  real(r8) :: effDim,shapeCoef,ubnd,lbnd, minMass

  limiterActive = limiter_is_on(props%min_mean_mass)
  effDim    = props%eff_dim
  shapeCoef = props%shape_coef
  lbnd      = props%lambda_bounds(1)
  ubnd      = props%lambda_bounds(2)
  minMass   = props%min_mean_mass
  present_n0 = present(n0)

  !$acc data present (props,qic,nic,lam,n0)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector collapse(2)
  do k = 1, dim2
     do i = 1, dim1
        if (qic(i,k) > qsmall) then
           ! add upper limit to in-cloud number concentration to prevent
           ! numerical error
           if (limiterActive) then
              nic(i,k) = min(nic(i,k), qic(i,k) / minMass)
           end if
   
           ! lambda = (c n/q)^(1/d)
           lam(i,k) = (shapeCoef * nic(i,k)/qic(i,k))**(1._r8/effDim)
   
           ! check for slope
           ! adjust vars
           if (lam(i,k) < lbnd) then
              lam(i,k) = lbnd
              nic(i,k) = lam(i,k)**(effDim) * qic(i,k)/shapeCoef
           else if (lam(i,k) > ubnd) then
              lam(i,k) = ubnd
              nic(i,k) = lam(i,k)**(effDim) * qic(i,k)/shapeCoef
           end if
   
        else
           lam(i,k) = 0._r8
        end if
     end do
  end do

  if (present_n0) then
     !$acc loop gang vector collapse(2)
     do k = 1, dim2
        do i = 1, dim1
           n0(i,k) = nic(i,k) * lam(i,k)
        end do
     end do
  end if
  !$acc end parallel

  !$acc end data
end subroutine size_dist_param_basic_2D

subroutine size_dist_param_basic_vect2(props, qic, nic, shapeC, lbnd, ubnd, lam, vlen, n0)

  type (mghydrometeorprops), intent(in)    :: props
  integer,                   intent(in)    :: vlen
  real(r8), dimension(vlen), intent(in)    :: qic
  real(r8), dimension(vlen), intent(inout) :: nic
  real(r8), dimension(vlen), intent(in)    :: shapeC,lbnd,ubnd
  real(r8), dimension(vlen), intent(out)   :: lam
  real(r8), dimension(vlen), intent(out), optional :: n0
  integer  :: i
  integer  :: cnt
  logical  :: limiterActive, present_n0
  real(r8) :: effDim,shapeCoef, minMass

  limiterActive = limiter_is_on(props%min_mean_mass)
  effDim        = props%eff_dim
  minMass       = props%min_mean_mass
  present_n0    = present(n0)

  !$acc data present (props,qic,nic,shapeC,lbnd,ubnd,lam,n0)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen

     if (qic(i) > qsmall) then
        ! add upper limit to in-cloud number concentration to prevent
        ! numerical error

        if (limiterActive) then
           nic(i) = min(nic(i), qic(i) / minMass)
        end if
        ! lambda = (c n/q)^(1/d)

        lam(i) = (shapeC(i) * nic(i)/qic(i))**(1._r8/effDim)
        ! check for slope
        ! adjust vars

        if (lam(i) < lbnd(i)) then
           lam(i) = lbnd(i)
           nic(i) = lam(i)**(effDim) * qic(i)/shapeC(i)
        else if (lam(i) > ubnd(i)) then
           lam(i) = ubnd(i)
           nic(i) = lam(i)**(effDim) * qic(i)/shapeC(i)
        end if
     else
        lam(i) = 0._r8
     end if

  end do

  if (present_n0) then
     !$acc loop gang vector
     do i = 1,vlen
        n0(i) = nic(i) * lam(i)
     end do
  end if
  !$acc end parallel

  !$acc end data
end subroutine size_dist_param_basic_vect2

real(r8) elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  ! Assumes that diameter follows an exponential distribution.
  real(r8), intent(in) :: q         ! mass mixing ratio
  real(r8), intent(in) :: n         ! number concentration (per volume)
  real(r8), intent(in) :: rho_air   ! local density of the air
  real(r8), intent(in) :: rho_sub   ! density of the particle substance

  avg_diameter = (q*rho_air/(pi * rho_sub * n))**(1._r8/3._r8)

end function avg_diameter

subroutine avg_diameter_vec (q, n, rho_air, rho_sub, avg_diameter, vlen)
   ! Finds the average diameter of particles given their density, and
   ! mass/number concentrations in the air.
   ! Assumes that diameter follows an exponential distribution.
   integer,  intent(in)  :: vlen
   real(r8), intent(in)  :: q(vlen)         ! mass mixing ratio
   real(r8), intent(in)  :: n(vlen)         ! number concentration (per volume)
   real(r8), intent(in)  :: rho_air(vlen)   ! local density of the air
   real(r8), intent(in)  :: rho_sub   ! density of the particle substance
   real(r8), intent(out) :: avg_diameter(vlen)
   integer :: i

   !$acc data present (q,n,rho_air,avg_diameter)

   !$acc parallel vector_length(VLENS) default(present)
   !$acc loop gang vector
   do i=1,vlen
      avg_diameter(i) = (q(i)*rho_air(i)/(pi * rho_sub * n(i)))**(1._r8/3._r8)
   end do
   !$acc end parallel

   !$acc end data
end subroutine avg_diameter_vec

subroutine var_coef_r8(relvar, a, res)
  !$acc routine seq

  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(r8), intent(in)  :: relvar
  real(r8), intent(in)  :: a
  real(r8), intent(out) :: res

  !$acc data present (relvar,res)

  call rising_factorial(relvar, a, res) 
  res = res / relvar**a

  !$acc end data
end subroutine var_coef_r8

subroutine var_coef_r8_vect(relvar, a, res, vlen)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  integer,  intent(in)  :: vlen
  real(r8), intent(in)  :: relvar(vlen)
  real(r8), intent(in)  :: a
  real(r8), intent(out) :: res(vlen)
  integer  :: i
  real(r8) :: tmpA(vlen)

  !$acc data present (relvar,res) &
  !$acc      create  (tmpA)

  call rising_factorial(relvar,a,tmpA,vlen)
  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     res(i) = tmpA(i)/relvar(i)**a
  end do
  !$acc end parallel

  !$acc end data
end subroutine var_coef_r8_vect

subroutine var_coef_integer(relvar, a, res)
  !$acc routine seq

  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(r8), intent(in) :: relvar
  integer, intent(in) :: a
  real(r8), intent(out) :: res

  !$acc data present (relvar,res)

  call rising_factorial(relvar, a, res)
  res = res / relvar**a

  !$acc end data
end subroutine var_coef_integer

subroutine var_coef_integer_vect(relvar, a, res, vlen)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  integer, intent(in)   :: vlen
  real(r8), intent(in)  :: relvar(vlen)
  integer, intent(in)   :: a
  real(r8), intent(out) :: res(vlen)
  integer  :: i
  real(r8) :: tmp(vlen)

  !$acc data present (relvar,res) &
  !$acc      create  (tmp)

  call rising_factorial(relvar, a,tmp,vlen)
  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     res(i) = tmp(i) / relvar(i)**a
  end do
  !$acc end parallel

  !$acc end data
end subroutine var_coef_integer_vect

!========================================================================
!MICROPHYSICAL PROCESS CALCULATIONS
!========================================================================
!========================================================================
! Initial ice deposition and sublimation loop.
! Run before the main loop
! This subroutine written by Peter Caldwell

subroutine ice_deposition_sublimation(t, qv, qi, ni, &
                                      icldm, rho, dv,qvl, qvi, &
                                      berg, vap_dep, ice_sublim, vlen)

  !INPUT VARS:
  !===============================================
  integer,  intent(in)                  :: vlen
  real(r8), dimension(vlen), intent(in) :: t
  real(r8), dimension(vlen), intent(in) :: qv
  real(r8), dimension(vlen), intent(in) :: qi
  real(r8), dimension(vlen), intent(in) :: ni
  real(r8), dimension(vlen), intent(in) :: icldm
  real(r8), dimension(vlen), intent(in) :: rho
  real(r8), dimension(vlen), intent(in) :: dv
  real(r8), dimension(vlen), intent(in) :: qvl
  real(r8), dimension(vlen), intent(in) :: qvi

  !OUTPUT VARS:
  !===============================================
  real(r8), dimension(vlen), intent(out) :: vap_dep !ice deposition (cell-ave value)
  real(r8), dimension(vlen), intent(out) :: ice_sublim !ice sublimation (cell-ave value)
  real(r8), dimension(vlen), intent(out) :: berg !bergeron enhancement (cell-ave value)

  !INTERNAL VARS:
  !===============================================
  real(r8) :: ab(vlen)
  real(r8) :: epsi
  real(r8) :: qiic(vlen)
  real(r8) :: niic(vlen)
  real(r8) :: lami(vlen)
  real(r8) :: n0i(vlen)
  integer :: i

  !$acc data present (t,qv,qi,ni,icldm,rho,dv,qvl) &
  !$acc      present (qvi,vap_dep,ice_sublim,berg) &
  !$acc      create  (ab,qiic,niic,lami,n0i)

  !GET IN-CLOUD qi, ni
  !===============================================

  !Compute linearized condensational heating correction
  call calc_ab(t, qvi, xxls, ab, vlen)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1,vlen
     qiic(i) = qi(i)/icldm(i)
     niic(i) = ni(i)/icldm(i)
     if (qi(i)<qsmall) then
        qiic(i) = 0._r8
        niic(i) = 0._r8
        ab(i)   = 0._r8
     end if
  end do
  !$acc end parallel

  !Get slope and intercept of gamma distn for ice.
  call size_dist_param_basic_vect(mg_ice_props, qiic, niic, lami, vlen, n0i)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(epsi)
  do i=1,vlen
     if (qi(i)>=qsmall) then
        !Get depletion timescale=1/eps
        epsi = 2._r8*pi*n0i(i)*rho(i)*Dv(i)/(lami(i)*lami(i))

        !Compute deposition/sublimation
        vap_dep(i) = epsi/ab(i)*(qv(i) - qvi(i))

        !Make this a grid-averaged quantity
        vap_dep(i) = vap_dep(i)*icldm(i)

        !Split into deposition or sublimation.
        if (t(i) < tmelt .and. vap_dep(i) > 0._r8) then
           ice_sublim(i) = 0._r8
        else
        ! make ice_sublim negative for consistency with other evap/sub processes
           ice_sublim(i) = min(vap_dep(i),0._r8)
           vap_dep(i)    = 0._r8
        end if

        !sublimation occurs @ any T. Not so for berg.
        if (t(i) < tmelt) then

           !Compute bergeron rate assuming cloud for whole step.
           berg(i) = max(epsi/ab(i)*(qvl(i) - qvi(i)), 0._r8)
        else   !T>frz
           berg(i)=0._r8
        end if !T<frz
     else      !where qi<qsmall
        berg(i)       = 0._r8
        vap_dep(i)    = 0._r8
        ice_sublim(i) = 0._r8
     end if    !qi>qsmall
  end do
  !$acc end parallel

  !$acc end data
end subroutine ice_deposition_sublimation

subroutine ice_deposition_sublimation_mg4(t, qv, qi, niic, &
                                      icldm, rho, dv,qvl, qvi, &
                                      berg, vap_dep, ice_sublim, &
                                      af1pr5, af1pr14, rhof, mu, sc, &
                                      vlen)

  !INPUT VARS:
  !===============================================
  integer,  intent(in)                  :: vlen
  real(r8), dimension(vlen), intent(in) :: t
  real(r8), dimension(vlen), intent(in) :: qv
  real(r8), dimension(vlen), intent(in) :: qi
  real(r8), dimension(vlen), intent(in) :: niic ! trude: input nicc as other routines
  real(r8), dimension(vlen), intent(in) :: icldm
  real(r8), dimension(vlen), intent(in) :: rho
  real(r8), dimension(vlen), intent(in) :: dv
  real(r8), dimension(vlen), intent(in) :: qvl
  real(r8), dimension(vlen), intent(in) :: qvi
  real(r8), dimension(vlen), intent(in) :: af1pr5, af1pr14
  real(r8), dimension(vlen), intent(in) :: rhof
  real(r8), dimension(vlen), intent(in) :: mu
  real(r8), dimension(vlen), intent(in) :: sc

  !OUTPUT VARS:
  !===============================================
  real(r8), dimension(vlen), intent(out) :: vap_dep !ice deposition (cell-ave value)
  real(r8), dimension(vlen), intent(out) :: ice_sublim !ice sublimation (cell-ave value)
  real(r8), dimension(vlen), intent(out) :: berg !bergeron enhancement (cell-ave value)

  !INTERNAL VARS:
  !===============================================
  real(r8) :: ab(vlen)
  real(r8) :: epsi
  real(r8) :: qiic(vlen)
  real(r8) :: lami(vlen)
  real(r8) :: n0i(vlen)
  real(r8), parameter :: thrd = 1._r8/3._r8
  integer :: i

  !$acc data present (t,qv,qi,icldm,rho,dv,qvl) &
  !$acc      present (qvi,vap_dep,ice_sublim,berg) &
  !$acc      present (niic,af1pr5,af1pr14,rhof,mu,sc) &
  !$acc      create  (ab,qiic,lami,n0i)

  !GET IN-CLOUD qi, ni
  !===============================================
  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i = 1,vlen
     if (qi(i)>=qsmall) then
        qiic(i) = qi(i)/icldm(i)
        !Compute linearized condensational heating correction
        call calc_ab(t(i), qvi(i), xxls, ab(i))
     end if
  end do
  !$acc end parallel

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(epsi)
  do i=1,vlen
     if (qi(i)>=qsmall) then
        if ( t(i) .lt. tmelt ) then
           epsi = (af1pr5(i)+af1pr14(i)*sc(i)**thrd*(rhof(i)*rho(i)/mu(i))**0.5_r8)*2._r8*pi* &
             rho(i)*dv(i)
        else
           epsi = 0._r8
        end if

        !Compute deposition/sublimation
        vap_dep(i) = epsi/ab(i)*(qv(i) - qvi(i))

        !Make this a grid-averaged quantity
        vap_dep(i) = vap_dep(i)*icldm(i)

        !Split into deposition or sublimation.
        if (t(i) < tmelt .and. vap_dep(i) > 0._r8) then
           ice_sublim(i)=0._r8
        else
           ! make ice_sublim negative for consistency with other evap/sub processes
           ice_sublim(i)=min(vap_dep(i),0._r8)
           vap_dep(i)=0._r8
        end if

        !sublimation occurs @ any T. Not so for berg.
        if (t(i) < tmelt) then
           !Compute bergeron rate assuming cloud for whole step.
           berg(i) = max(epsi/ab(i)*(qvl(i) - qvi(i)), 0._r8)
        else   !T>frz
           berg(i) = 0._r8
        end if !T<frz
     else      !where qi<qsmall
        berg(i)       = 0._r8
        vap_dep(i)    = 0._r8
        ice_sublim(i) = 0._r8
     end if    !qi>qsmall
  end do
  !$acc end parallel

  !$acc end data
end subroutine ice_deposition_sublimation_mg4

!========================================================================
! autoconversion of cloud liquid water to rain
! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
! minimum qc of 1 x 10^-8 prevents floating point error

subroutine kk2000_liq_autoconversion(microp_uniform, qcic, &
     ncic, rho, relvar, prc, nprc, nprc1, vlen)

  integer, intent(in) :: vlen
  logical, intent(in) :: microp_uniform

  real(r8), dimension(vlen), intent(in) :: qcic
  real(r8), dimension(vlen), intent(in) :: ncic
  real(r8), dimension(vlen), intent(in) :: rho

  real(r8), dimension(vlen), intent(in) :: relvar

  real(r8), dimension(vlen), intent(out) :: prc
  real(r8), dimension(vlen), intent(out) :: nprc
  real(r8), dimension(vlen), intent(out) :: nprc1

  real(r8), dimension(vlen) :: prc_coef
  integer :: i

  !$acc data present (qcic,ncic,rho,relvar,prc,nprc,nprc1) &
  !$acc      create  (prc_coef)

  ! Take variance into account, or use uniform value.
  if (.not. microp_uniform) then
     call var_coef(relvar, 2.47_r8, prc_coef,vlen)
  else
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector
     do i = 1,vlen
        prc_coef(i) = 1._r8
     end do
     !$acc end parallel
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if (qcic(i) >= icsmall) then
        ! nprc is increase in rain number conc due to autoconversion
        ! nprc1 is decrease in cloud droplet conc due to autoconversion

        ! assume exponential sub-grid distribution of qc, resulting in additional
        ! factor related to qcvar below
        ! switch for sub-columns, don't include sub-grid qc
        prc(i)   = prc_coef(i) * &
             0.01_r8 * 1350._r8 * qcic(i)**2.47_r8 * (ncic(i)*1.e-6_r8*rho(i))**(-1.1_r8)
        nprc(i)  = prc(i) * (1._r8/droplet_mass_25um)
        nprc1(i) = prc(i)*ncic(i)/qcic(i)

     else
        prc(i)   = 0._r8
        nprc(i)  = 0._r8
        nprc1(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine kk2000_liq_autoconversion
  
!========================================================================

subroutine sb2001v2_liq_autoconversion(pgam,qc,nc,qr,rho,relvar,au,nprc,nprc1,vlen)
  !
  ! ---------------------------------------------------------------------
  ! AUTO_SB:  calculates the evolution of mass- and number mxg-ratio for
  ! drizzle drops due to autoconversion. The autoconversion rate assumes
  ! f(x)=A*x**(nu_c)*exp(-Bx) in drop MASS x. 

  ! Code from Hugh Morrison, Sept 2014

  ! autoconversion
  ! use simple lookup table of dnu values to get mass spectral shape parameter
  ! equivalent to the size spectral shape parameter pgam
    
  integer,                   intent(in)     :: vlen  
  real(r8), dimension(vlen), intent (in)    :: pgam
  real(r8), dimension(vlen), intent (in)    :: qc  ! = qc (cld water mixing ratio)
  real(r8), dimension(vlen), intent (in)    :: nc  ! = nc (cld water number conc /kg)    
  real(r8), dimension(vlen), intent (in)    :: qr  ! = qr (rain water mixing ratio)
  real(r8), dimension(vlen), intent (in)    :: rho ! = rho : density profile
  real(r8), dimension(vlen), intent (in)    :: relvar 
  
  real(r8), dimension(vlen), intent (out)   :: au ! = prc autoconversion rate
  real(r8), dimension(vlen), intent (out)   :: nprc1 ! = number tendency
  real(r8), dimension(vlen), intent (out)   :: nprc ! = number tendency fixed size for rain
 
  ! parameters for droplet mass spectral shape, 
  !used by Seifert and Beheng (2001)                             
  ! warm rain scheme only (iparam = 1)                                                                        
  real(r8), parameter :: dnu(16) = [0._r8,-0.557_r8,-0.430_r8,-0.307_r8, & 
     -0.186_r8,-0.067_r8,0.050_r8,0.167_r8,0.282_r8,0.397_r8,0.512_r8, &
     0.626_r8,0.739_r8,0.853_r8,0.966_r8,0.966_r8]

  ! parameters for Seifert and Beheng (2001) autoconversion/accretion                                         
  real(r8), parameter :: kc = 9.44e9_r8
  real(r8), parameter :: kr = 5.78e3_r8
  real(r8) :: dum, dum1, nu
  integer  :: dumi, i

  !$acc data present (pgam,qc,nc,qr,rho,relvar,au,nprc1,nprc)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(dumi,nu,dum,dum1)
  do i=1,vlen
     if (qc(i) > qsmall) then
       dumi = int(pgam(i))
       nu   = dnu(dumi)+(dnu(dumi+1)-dnu(dumi))* &
               (pgam(i)-dumi)

       dum  = 1._r8-qc(i)/(qc(i)+qr(i))
       dum1 = 600._r8*dum**0.68_r8*(1._r8-dum**0.68_r8)**3

       au(i) = kc/(20._r8*2.6e-7_r8)* &
         (nu+2._r8)*(nu+4._r8)/(nu+1._r8)**2._r8* &
         (rho(i)*qc(i)/1000._r8)**4._r8/(rho(i)*nc(i)/1.e6_r8)**2._r8* &
         (1._r8+dum1/(1._r8-dum)**2)*1000._r8 / rho(i)

       nprc1(i) = au(i)*2._r8/2.6e-7_r8*1000._r8
       nprc(i)  = au(i)/droplet_mass_40um
     else
       au(i)    = 0._r8
       nprc1(i) = 0._r8
       nprc(i)  = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine sb2001v2_liq_autoconversion 
  
!========================================================================
!SB2001 Accretion V2

subroutine sb2001v2_accre_cld_water_rain(qc,nc,qr,rho,relvar,pra,npra,vlen)
  !
  ! ---------------------------------------------------------------------
  ! ACCR_SB calculates the evolution of mass mxng-ratio due to accretion
  ! and self collection following Seifert & Beheng (2001).  
  !
  
  integer, intent(in) :: vlen
  
  real(r8), dimension(vlen), intent (in)    :: qc  ! = qc (cld water mixing ratio)
  real(r8), dimension(vlen), intent (in)    :: nc  ! = nc (cld water number conc /kg)    
  real(r8), dimension(vlen), intent (in)    :: qr  ! = qr (rain water mixing ratio)
  real(r8), dimension(vlen), intent (in)    :: rho ! = rho : density profile
  real(r8), dimension(vlen), intent (in)    :: relvar

  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: pra  ! MMR
  real(r8), dimension(vlen), intent(out) :: npra ! Number

  ! parameters for Seifert and Beheng (2001) autoconversion/accretion                                         
  real(r8), parameter :: kc = 9.44e9_r8
  real(r8), parameter :: kr = 5.78e3_r8

  real(r8) :: dum, dum1
  integer :: i

  ! accretion

  !$acc data present (qc,nc,qr,rho,relvar,pra,npra)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(dum,dum1)
  do i = 1,vlen
    if (qc(i) > qsmall) then
      dum     = 1._r8-qc(i)/(qc(i)+qr(i))
      dum1    = (dum/(dum+5.e-4_r8))**4._r8
      pra(i)  = kr*rho(i)*0.001_r8*qc(i)*qr(i)*dum1
      npra(i) = pra(i)*rho(i)*0.001_r8*(nc(i)*rho(i)*1.e-6_r8)/ &
           (qc(i)*rho(i)*0.001_r8)*1.e6_r8 / rho(i)
    else
      pra(i)  = 0._r8
      npra(i) = 0._r8
    end if 
  end do
  !$acc end parallel

  !$acc end data
end subroutine sb2001v2_accre_cld_water_rain   

!========================================================================
! Autoconversion of cloud ice to snow
! similar to Ferrier (1994)

subroutine ice_autoconversion(t, qiic, lami, n0i, dcs, prci, nprci, vlen)

  integer, intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t
  real(r8), dimension(vlen), intent(in) :: qiic
  real(r8), dimension(vlen), intent(in) :: lami
  real(r8), dimension(vlen), intent(in) :: n0i
  real(r8),                    intent(in) :: dcs

  real(r8), dimension(vlen), intent(out) :: prci
  real(r8), dimension(vlen), intent(out) :: nprci

  ! Assume autoconversion timescale of 180 seconds.
  real(r8), parameter :: ac_time = 180._r8

  ! Average mass of an ice particle.
  real(r8) :: m_ip
  ! Ratio of autoconversion diameter to average diameter.
  real(r8) :: d_rat
  integer :: i

  !$acc data present (t,qiic,lami,n0i,prci,nprci)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(d_rat,m_ip)
  do i=1,vlen
     if (t(i) <= tmelt .and. qiic(i) >= qsmall) then
        d_rat = lami(i)*dcs

        ! Rate of ice particle conversion (number).
        nprci(i) = n0i(i)/(lami(i)*ac_time)*exp(-d_rat)
        m_ip = (rhoi*pi/6._r8) / lami(i)**3

        ! Rate of mass conversion.
        ! Note that this is:
        ! m n (d^3 + 3 d^2 + 6 d + 6)
        prci(i) = m_ip * nprci(i) * &
             (((d_rat + 3._r8)*d_rat + 6._r8)*d_rat + 6._r8)
     else
        prci(i)  = 0._r8
        nprci(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine ice_autoconversion

! immersion freezing (Bigg, 1953)
!===================================

subroutine immersion_freezing(microp_uniform, t, pgam, lamc, &
     qcic, ncic, relvar, mnuccc, nnuccc, vlen)

  integer, intent(in) :: vlen
  logical, intent(in) :: microp_uniform

  ! Temperature
  real(r8), dimension(vlen), intent(in) :: t

  ! Cloud droplet size distribution parameters
  real(r8), dimension(vlen), intent(in) :: pgam
  real(r8), dimension(vlen), intent(in) :: lamc

  ! MMR and number concentration of in-cloud liquid water
  real(r8), dimension(vlen), intent(in) :: qcic
  real(r8), dimension(vlen), intent(in) :: ncic

  ! Relative variance of cloud water
  real(r8), dimension(vlen), intent(in) :: relvar

  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: mnuccc ! MMR
  real(r8), dimension(vlen), intent(out) :: nnuccc ! Number

  ! Coefficients that will be omitted for sub-columns
  real(r8), dimension(vlen) :: dum
  integer  :: i
  real(r8) :: tmp

  !$acc data present (t,pgam,lamc,qcic,ncic) &
  !$acc      present (relvar,mnuccc,nnuccc)  &
  !$acc      create  (dum)

  if (.not. microp_uniform) then
     call var_coef(relvar, 2, dum, vlen)
  else
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector
     do i =1,vlen
        dum(i) = 1._r8
     end do
     !$acc end parallel
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(tmp)
  do i=1,vlen
     if (qcic(i) >= qsmall .and. t(i) < 269.15_r8) then
        call rising_factorial(pgam(i)+1._r8, 3, tmp)
        nnuccc(i) = &
             pi/6._r8*ncic(i)*tmp* &
             bimm*(exp(aimm*(tmelt - t(i)))-1._r8)/lamc(i)**3

        call rising_factorial(pgam(i)+4._r8, 3, tmp)
        mnuccc(i) = dum(i) * nnuccc(i) * &
             pi/6._r8*rhow* &
             tmp/lamc(i)**3
     else
        mnuccc(i) = 0._r8
        nnuccc(i) = 0._r8
     end if ! qcic > qsmall and t < 4 deg C
  end do
  !$acc end parallel

  !$acc end data
end subroutine immersion_freezing

! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
!===================================================================
! dust size and number in multiple bins are read in from companion routine

subroutine contact_freezing (microp_uniform, t, p, rndst, nacon, &
     pgam, lamc, qcic, ncic, relvar, mnucct, nnucct, vlen, mdust)

  logical, intent(in) :: microp_uniform

  integer, intent(in) :: vlen
  integer, intent(in) :: mdust

  real(r8), dimension(vlen), intent(in) :: t            ! Temperature
  real(r8), dimension(vlen), intent(in) :: p            ! Pressure
  real(r8), dimension(vlen, mdust), intent(in) :: rndst ! Radius (for multiple dust bins)
  real(r8), dimension(vlen, mdust), intent(in) :: nacon ! Number (for multiple dust bins)

  ! Size distribution parameters for cloud droplets
  real(r8), dimension(vlen), intent(in) :: pgam
  real(r8), dimension(vlen), intent(in) :: lamc

  ! MMR and number concentration of in-cloud liquid water
  real(r8), dimension(vlen), intent(in) :: qcic
  real(r8), dimension(vlen), intent(in) :: ncic

  ! Relative cloud water variance
  real(r8), dimension(vlen), intent(in) :: relvar

  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: mnucct ! MMR
  real(r8), dimension(vlen), intent(out) :: nnucct ! Number

  real(r8) :: tcnt                  ! scaled relative temperature
  real(r8) :: viscosity             ! temperature-specific viscosity (kg/m/s)
  real(r8) :: mfp                   ! temperature-specific mean free path (m)

  ! Dimension these according to number of dust bins, inferred from rndst size
  real(r8) :: nslip(size(rndst,2))  ! slip correction factors
  real(r8) :: ndfaer(size(rndst,2)) ! aerosol diffusivities (m^2/sec)

  ! Coefficients not used for subcolumns
  real(r8) :: dum, dum1, tmp

  ! Common factor between mass and number.
  real(r8) :: contact_factor

  integer  :: i, j

  !$acc data present (t,p,rndst,nacon,pgam,lamc) &
  !$acc      present (qcic,ncic,relvar,mnucct,nnucct) &
  !$acc      create  (nslip,ndfaer)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(dum,dum1,tcnt,viscosity,mfp,contact_factor,tmp)
  do i = 1,vlen
     if (qcic(i) >= qsmall .and. t(i) < 269.15_r8) then
        if (microp_uniform) then
           dum  = 1._r8
           dum1 = 1._r8
        else
           call var_coef(relvar(i), 4._r8/3._r8, dum)
           call var_coef(relvar(i), 1._r8/3._r8, dum1)
        end if

        tcnt=(270.16_r8-t(i))**1.3_r8
        viscosity = 1.8e-5_r8*(t(i)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
        mfp = 2.0_r8*viscosity/ &                         ! Mean free path (m)
                     (p(i)*sqrt( 8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i)) ))

        do j = 1, mdust
           ! Note that these two are vectors.
           nslip(j) = 1.0_r8+(mfp/rndst(i,j))*(1.257_r8+(0.4_r8*exp(-(1.1_r8*rndst(i,j)/mfp))))! Slip correction factor
   
           ndfaer(j) = 1.381e-23_r8*t(i)*nslip(j)/(6._r8*pi*viscosity*rndst(i,j))  ! aerosol diffusivity (m2/s)
        end do

        contact_factor = dot_product(ndfaer,nacon(i,:)*tcnt) * pi * &
                                     ncic(i) * (pgam(i) + 1._r8) / lamc(i)
        call rising_factorial(pgam(i)+2._r8, 3, tmp)
        mnucct(i) = dum * contact_factor * pi/3._r8*rhow*tmp/lamc(i)**3
        nnucct(i) = dum1 * 2._r8 * contact_factor
     else
        mnucct(i) = 0._r8
        nnucct(i) = 0._r8
     end if ! qcic > qsmall and t < 4 deg C
  end do
  !$acc end parallel

  !$acc end data
end subroutine contact_freezing

! snow self-aggregation from passarelli, 1978, used by reisner, 1998
!===================================================================
! this is hard-wired for bs = 0.4 for now
! ignore self-collection of cloud ice

subroutine snow_self_aggregation(t, rho, asn, rhosn, qsic, nsic, nsagg, vlen)

  integer,                   intent(in) :: vlen

  real(r8), dimension(vlen), intent(in) :: t     ! Temperature
  real(r8), dimension(vlen), intent(in) :: rho   ! Density
  real(r8), dimension(vlen), intent(in) :: asn   ! fall speed parameter for snow
  real(r8),                  intent(in) :: rhosn ! density of snow

  ! In-cloud snow
  real(r8), dimension(vlen), intent(in) :: qsic ! MMR
  real(r8), dimension(vlen), intent(in) :: nsic ! Number

  ! Output number tendency
  real(r8), dimension(vlen), intent(out) :: nsagg

  integer :: i

  !$acc data present (t,rho,asn,qsic,nsic,nsagg)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if (qsic(i) >= qsmall .and. t(i) <= tmelt) then
        nsagg(i) = -1108._r8*eii/(4._r8*720._r8*rhosn)*asn(i)*qsic(i)*nsic(i)*rho(i)*&
             ((qsic(i)/nsic(i))*(1._r8/(rhosn*pi)))**((bs-1._r8)/3._r8)
     else
        nsagg(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine snow_self_aggregation

! ice self-aggregation used in mg4 
!===================================================================

subroutine ice_self_aggregation(t, rho, rhof, af1pr3, qiic, nsagg, vlen)

  integer,                   intent(in) :: vlen

  real(r8), dimension(vlen), intent(in) :: t        ! Temperature
  real(r8), dimension(vlen), intent(in) :: rho      ! Density
  real(r8), dimension(vlen), intent(in) :: rhof     ! Density correction
  
  real(r8), dimension(vlen), intent(in) :: af1pr3   ! process rate 

  ! In-cloud ice
  real(r8), dimension(vlen), intent(in) :: qiic ! MMR

  ! Collection efficiency
  real(r8), parameter :: eii = 0.1_r8

  ! Output number tendency
  real(r8), dimension(vlen), intent(out) :: nsagg

  integer :: i

  !$acc data present (t,rho,rhof,af1pr3,qiic,nsagg)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen  
     if (qiic(i) >= qsmall .and. t(i) <= tmelt) then
        nsagg(i) = af1pr3(i)*rho(i)*eii*rhof(i)
     else
        nsagg(i)=0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine ice_self_aggregation

! accretion of cloud droplets onto snow/graupel
!===================================================================
! here use continuous collection equation with
! simple gravitational collection kernel
! ignore collisions between droplets/cloud ice
! since minimum size ice particle for accretion is 50 - 150 micron

subroutine accrete_cloud_water_snow(t, rho, asn, uns, mu, qcic, ncic, qsic, &
     pgam, lamc, lams, n0s, psacws, npsacws, vlen)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t   ! Temperature
  real(r8), dimension(vlen), intent(in) :: rho ! Density
  real(r8), dimension(vlen), intent(in) :: asn ! Fallspeed parameter (snow)
  real(r8), dimension(vlen), intent(in) :: uns ! Current fallspeed   (snow)
  real(r8), dimension(vlen), intent(in) :: mu  ! Viscosity

  ! In-cloud liquid water
  real(r8), dimension(vlen), intent(in) :: qcic ! MMR
  real(r8), dimension(vlen), intent(in) :: ncic ! Number

  ! In-cloud snow
  real(r8), dimension(vlen), intent(in) :: qsic ! MMR

  ! Cloud droplet size parameters
  real(r8), dimension(vlen), intent(in) :: pgam
  real(r8), dimension(vlen), intent(in) :: lamc

  ! Snow size parameters
  real(r8), dimension(vlen), intent(in) :: lams
  real(r8), dimension(vlen), intent(in) :: n0s

  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: psacws  ! Mass mixing ratio
  real(r8), dimension(vlen), intent(out) :: npsacws ! Number concentration

  real(r8) :: dc0 ! Provisional mean droplet size
  real(r8) :: dum
  real(r8) :: eci ! collection efficiency for riming of snow by droplets

  ! Fraction of cloud droplets accreted per second
  real(r8) :: accrete_rate
  integer :: i

  ! ignore collision of snow with droplets above freezing

  !$acc data present (t,rho,asn,uns,mu,qcic,ncic,qsic) &
  !$acc      present (pgam,lamc,lams,n0s,psacws,npsacws)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(dc0,dum,eci,accrete_rate)
  do i=1,vlen
     if (qsic(i) >= qsmall .and. t(i) <= tmelt .and. qcic(i) >= qsmall) then
        ! put in size dependent collection efficiency
        ! mean diameter of snow is area-weighted, since
        ! accretion is function of crystal geometric area
        ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)
        dc0 = (pgam(i)+1._r8)/lamc(i)
        dum = dc0*dc0*uns(i)*rhow*lams(i)/(9._r8*mu(i))
        eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))
        eci = max(eci,0._r8)
        eci = min(eci,1._r8)

        ! no impact of sub-grid distribution of qc since psacws
        ! is linear in qc
        accrete_rate = pi/4._r8*asn(i)*rho(i)*n0s(i)*eci*gamma_bs_plus3 / lams(i)**(bs+3._r8)
        psacws(i)  = accrete_rate*qcic(i)
        npsacws(i) = accrete_rate*ncic(i)
     else
        psacws(i)  = 0._r8
        npsacws(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine accrete_cloud_water_snow

! Version used in MG4
subroutine accrete_cloud_water_ice(t, rho, rhof, af1pr4, qcic, ncic, qiic, psacws, npsacws, vlen)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t        ! Temperature
  real(r8), dimension(vlen), intent(in) :: rho      ! Density
  real(r8), dimension(vlen), intent(in) :: rhof     ! Density correction factor 
  real(r8), dimension(vlen), intent(in) :: af1pr4   ! Process rate from look-up table
  ! In-cloud liquid water
  real(r8), dimension(vlen), intent(in) :: qcic     ! MMR
  real(r8), dimension(vlen), intent(in) :: ncic     ! Number
  ! In-cloud ice
  real(r8), dimension(vlen), intent(in) :: qiic     ! MMR
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: psacws  ! Mass mixing ratio
  real(r8), dimension(vlen), intent(out) :: npsacws ! Number concentration

  real(r8) :: eci     ! collection efficiency for riming of snow by droplets
  integer  :: i

  ! ignore collision of snow with droplets above freezing
  eci = 0.5_r8 ! from p3 program
 
  !$acc data present (t,rho,rhof,af1pr4,qcic) &
  !$acc      present (ncic,qiic,psacws,npsacws)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if (qiic(i) >= qsmall .and. t(i) <= tmelt .and. qcic(i) >= qsmall) then
         psacws(i)  = rhof(i)*af1pr4(i)*qcic(i)*eci*rho(i)
         npsacws(i) = rhof(i)*af1pr4(i)*ncic(i)*eci*rho(i)
     else
         psacws(i)  = 0._r8
         npsacws(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine accrete_cloud_water_ice

! add secondary ice production due to accretion of droplets by snow
!===================================================================
! (Hallet-Mossop process) (from Cotton et al., 1986)

subroutine secondary_ice_production(t, psacws, msacwi, nsacwi, vlen)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t ! Temperature

  ! Accretion of cloud water to snow tendencies
  real(r8), dimension(vlen), intent(inout) :: psacws ! MMR

  ! Output (ice) tendencies
  real(r8), dimension(vlen), intent(out) :: msacwi ! MMR
  real(r8), dimension(vlen), intent(out) :: nsacwi ! Number
  integer :: i

  !$acc data present (t,psacws,msacwi,nsacwi)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if((t(i) < 270.16_r8) .and. (t(i) >= 268.16_r8)) then
        nsacwi(i) = 3.5e8_r8*(270.16_r8-t(i))/2.0_r8*psacws(i)
     else if((t(i) < 268.16_r8) .and. (t(i) >= 265.16_r8)) then
        nsacwi(i) = 3.5e8_r8*(t(i)-265.16_r8)/3.0_r8*psacws(i)
     else
        nsacwi(i) = 0.0_r8
     endif
  end do

  !$acc loop gang vector
  do i=1,vlen
     msacwi(i) = min(nsacwi(i)*mi0, psacws(i))
     psacws(i) = psacws(i) - msacwi(i)
  end do
  !$acc end parallel

  !$acc end data
end subroutine secondary_ice_production

! accretion of rain water by snow
!===================================================================
! formula from ikawa and saito, 1991, used by reisner et al., 1998

subroutine accrete_rain_snow(t, rho, umr, ums, unr, uns, qric, qsic, &
     lamr, n0r, lams, n0s, pracs, npracs, vlen)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t   ! Temperature
  real(r8), dimension(vlen), intent(in) :: rho ! Density
  ! Fallspeeds
  ! mass-weighted
  real(r8), dimension(vlen), intent(in) :: umr ! rain
  real(r8), dimension(vlen), intent(in) :: ums ! snow
  ! number-weighted
  real(r8), dimension(vlen), intent(in) :: unr ! rain
  real(r8), dimension(vlen), intent(in) :: uns ! snow
  ! In cloud MMRs
  real(r8), dimension(vlen), intent(in) :: qric ! rain
  real(r8), dimension(vlen), intent(in) :: qsic ! snow
  ! Size distribution parameters
  ! rain
  real(r8), dimension(vlen), intent(in) :: lamr
  real(r8), dimension(vlen), intent(in) :: n0r
  ! snow
  real(r8), dimension(vlen), intent(in) :: lams
  real(r8), dimension(vlen), intent(in) :: n0s
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: pracs  ! MMR
  real(r8), dimension(vlen), intent(out) :: npracs ! Number
  ! Collection efficiency for accretion of rain by snow
  real(r8), parameter :: ecr = 1.0_r8
  ! Ratio of average snow diameter to average rain diameter.
  real(r8) :: d_rat
  ! Common factor between mass and number expressions
  real(r8) :: common_factor
  integer  :: i

  !$acc data present (t,rho,umr,ums,unr,uns,qric,qsic) &
  !$acc      present (lamr,n0r,lams,n0s,pracs,npracs)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(common_factor,d_rat)
  do i=1,vlen
     if (qric(i) >= icsmall .and. qsic(i) >= icsmall .and. t(i) <= tmelt) then
        common_factor = pi*ecr*rho(i)*n0r(i)*n0s(i)/(lamr(i)**3 * lams(i))
        d_rat = lamr(i)/lams(i)
        pracs(i) = common_factor*pi*rhow* &
             sqrt((1.2_r8*umr(i)-0.95_r8*ums(i))**2 + 0.08_r8*ums(i)*umr(i)) / lamr(i)**3 * &
             ((0.5_r8*d_rat + 2._r8)*d_rat + 5._r8)
        npracs(i) = common_factor*0.5_r8* &
             sqrt(1.7_r8*(unr(i)-uns(i))**2 + 0.3_r8*unr(i)*uns(i)) * &
             ((d_rat + 1._r8)*d_rat + 1._r8)
     else
        pracs(i) = 0._r8
        npracs(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine accrete_rain_snow

! Version used in MG4
subroutine accrete_rain_ice(t, rho, rhof, af1pr8, af1pr7, qric, qiic, &
      n0r, pracs, npracs, vlen )

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: rhof ! Density correction factor 
  real(r8), dimension(vlen), intent(in) :: rho  ! Density
  real(r8), dimension(vlen), intent(in) :: t    ! Temperature
  ! In cloud MMRs
  real(r8), dimension(vlen), intent(in) :: qric ! rain
  real(r8), dimension(vlen), intent(in) :: qiic ! combined snow and ice
  ! Size distribution parameters
  ! rain
  real(r8), dimension(vlen), intent(in) :: n0r
  ! Processrates 
  real(r8), dimension(vlen), intent(in) :: af1pr8  
  real(r8), dimension(vlen), intent(in) :: af1pr7  
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: pracs  ! MMR
  real(r8), dimension(vlen), intent(out) :: npracs ! Number
  ! Collection efficiency for accretion of rain by snow
  real(r8), parameter :: ecr = 1.0_r8
  integer :: i

  !$acc data present (t,rho,rhof,qric,qiic,n0r) &
  !$acc      present (af1pr8,af1pr7,pracs,npracs)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if (qric(i) >= icsmall .and. qiic(i) >= icsmall .and. t(i) <= tmelt) then
        pracs(i)  = 10._r8**(af1pr8(i)+dlog10(n0r(i)))*rho(i)*rhof(i)*ecr
        npracs(i) = 10._r8**(af1pr7(i)+dlog10(n0r(i)))*rho(i)*rhof(i)*ecr
     else
        pracs(i)  = 0._r8
        npracs(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine accrete_rain_ice

! heterogeneous freezing of rain drops
!===================================================================
! follows from Bigg (1953)

subroutine heterogeneous_rain_freezing(t, qric, nric, lamr, mnuccr, nnuccr, vlen)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t    ! Temperature
  ! In-cloud rain
  real(r8), dimension(vlen), intent(in) :: qric ! MMR
  real(r8), dimension(vlen), intent(in) :: nric ! Number
  real(r8), dimension(vlen), intent(in) :: lamr ! size parameter
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: mnuccr ! MMR
  real(r8), dimension(vlen), intent(out) :: nnuccr ! Number
  integer :: i

  !$acc data present (t,qric,nric,lamr,mnuccr,nnuccr)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if (t(i) < 269.15_r8 .and. qric(i) >= qsmall) then
        nnuccr(i) = pi*nric(i)*bimm* &
             (exp(aimm*(tmelt - t(i)))-1._r8)/lamr(i)**3
        mnuccr(i) = nnuccr(i) * 20._r8*pi*rhow/lamr(i)**3
     else
        mnuccr(i) = 0._r8
        nnuccr(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine heterogeneous_rain_freezing

! accretion of cloud liquid water by rain
!===================================================================
! formula from Khrouditnov and Kogan (2000)
! gravitational collection kernel, droplet fall speed neglected

subroutine accrete_cloud_water_rain(microp_uniform, qric, qcic, &
     ncic, relvar, accre_enhan, pra, npra, vlen)

  logical, intent(in) :: microp_uniform
  integer, intent(in) :: vlen
  ! In-cloud rain
  real(r8), dimension(vlen), intent(in) :: qric ! MMR
  ! Cloud droplets
  real(r8), dimension(vlen), intent(in) :: qcic ! MMR
  real(r8), dimension(vlen), intent(in) :: ncic ! Number
  ! SGS variability
  real(r8), dimension(vlen), intent(in) :: relvar
  real(r8), dimension(vlen), intent(in) :: accre_enhan
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: pra  ! MMR
  real(r8), dimension(vlen), intent(out) :: npra ! Number
  ! Coefficient that varies for subcolumns
  real(r8), dimension(vlen) :: pra_coef
  integer :: i

  !$acc data present (qric,qcic,ncic,relvar,accre_enhan,pra,npra) &
  !$acc      create  (pra_coef)

  if (.not. microp_uniform) then
     call var_coef(relvar, 1.15_r8, pra_coef, vlen)
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector
     do i = 1, vlen
        pra_coef(i) = accre_enhan(i) * pra_coef(i)
     end do
     !$acc end parallel
  else
     !$acc parallel vector_length(VLENS) default(present)
     !$acc loop gang vector
     do i = 1, vlen
        pra_coef(i) = 1._r8
     end do
     !$acc end parallel
  end if

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if (qric(i) >= qsmall .and. qcic(i) >= qsmall) then
        ! include sub-grid distribution of cloud water
        pra(i)  = pra_coef(i) * 67._r8*(qcic(i)*qric(i))**1.15_r8
        npra(i) = pra(i)*ncic(i)/qcic(i)
     else
        pra(i)  = 0._r8
        npra(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine accrete_cloud_water_rain

! Self-collection of rain drops
!===================================================================
! from Beheng(1994)

subroutine self_collection_rain(rho, qric, nric, nragg, vlen)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: rho  ! Air density
  ! Rain
  real(r8), dimension(vlen), intent(in) :: qric ! MMR
  real(r8), dimension(vlen), intent(in) :: nric ! Number
  ! Output number tendency
  real(r8), dimension(vlen), intent(out) :: nragg
  integer :: i

  !$acc data present (rho,qric,nric,nragg)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if (qric(i) >= qsmall) then
        nragg(i) = -8._r8*nric(i)*qric(i)*rho(i)
     else
        nragg(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine self_collection_rain

! Accretion of cloud ice by snow
!===================================================================
! For this calculation, it is assumed that the Vs >> Vi
! and Ds >> Di for continuous collection

subroutine accrete_cloud_ice_snow(t, rho, asn, qiic, niic, qsic, &
     lams, n0s, prai, nprai, vlen)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t    ! Temperature
  real(r8), dimension(vlen), intent(in) :: rho   ! Density
  real(r8), dimension(vlen), intent(in) :: asn  ! Snow fallspeed parameter
  ! Cloud ice
  real(r8), dimension(vlen), intent(in) :: qiic ! MMR
  real(r8), dimension(vlen), intent(in) :: niic ! Number
  real(r8), dimension(vlen), intent(in) :: qsic ! Snow MMR
  ! Snow size parameters
  real(r8), dimension(vlen), intent(in) :: lams
  real(r8), dimension(vlen), intent(in) :: n0s
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: prai ! MMR
  real(r8), dimension(vlen), intent(out) :: nprai ! Number
  ! Fraction of cloud ice particles accreted per second
  real(r8) :: accrete_rate
  integer  :: i

  !$acc data present (t,rho,asn,qiic,niic,qsic,lams,n0s,prai,nprai)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(accrete_rate)
  do i=1,vlen
     if (qsic(i) >= qsmall .and. qiic(i) >= qsmall .and. t(i) <= tmelt) then
        accrete_rate = pi/4._r8 * eii * asn(i) * rho(i) * n0s(i) * gamma_bs_plus3/ &
             lams(i)**(bs+3._r8)
        prai(i)  = accrete_rate * qiic(i)
        nprai(i) = accrete_rate * niic(i)
     else
        prai(i)  = 0._r8
        nprai(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine accrete_cloud_ice_snow

! calculate evaporation/sublimation of rain and snow
!===================================================================
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

subroutine evaporate_sublimate_precip(t, rho, dv, mu, sc, q, qvl, qvi, &
     lcldm, precip_frac, arn, asn, qcic, qiic, qric, qsic, lamr, n0r, lams, n0s, &
     pre, prds, am_evp_st, vlen, evap_rhthrsh)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t           ! temperature
  real(r8), dimension(vlen), intent(in) :: rho         ! air density
  real(r8), dimension(vlen), intent(in) :: dv          ! water vapor diffusivity
  real(r8), dimension(vlen), intent(in) :: mu          ! viscosity
  real(r8), dimension(vlen), intent(in) :: sc          ! schmidt number
  real(r8), dimension(vlen), intent(in) :: q           ! humidity
  real(r8), dimension(vlen), intent(in) :: qvl         ! saturation humidity (water)
  real(r8), dimension(vlen), intent(in) :: qvi         ! saturation humidity (ice)
  real(r8), dimension(vlen), intent(in) :: lcldm       ! liquid cloud fraction
  real(r8), dimension(vlen), intent(in) :: precip_frac ! precipitation fraction (maximum overlap)
  ! fallspeed parameters
  real(r8), dimension(vlen), intent(in) :: arn         ! rain
  real(r8), dimension(vlen), intent(in) :: asn         ! snow
  ! In-cloud MMRs
  real(r8), dimension(vlen), intent(in) :: qcic        ! cloud liquid
  real(r8), dimension(vlen), intent(in) :: qiic        ! cloud ice
  real(r8), dimension(vlen), intent(in) :: qric        ! rain
  real(r8), dimension(vlen), intent(in) :: qsic        ! snow
  ! Size parameters
  ! rain
  real(r8), dimension(vlen), intent(in) :: lamr
  real(r8), dimension(vlen), intent(in) :: n0r
  ! snow
  real(r8), dimension(vlen), intent(in) :: lams
  real(r8), dimension(vlen), intent(in) :: n0s
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: pre
  real(r8), dimension(vlen), intent(out) :: prds
  real(r8), dimension(vlen), intent(out) :: am_evp_st  ! Fractional area where rain evaporates.
  ! Switch
  logical, intent(in)                    :: evap_rhthrsh 

  real(r8) :: qclr   ! water vapor mixing ratio in clear air
  real(r8) :: ab,abr ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale
  real(r8), dimension(vlen) :: dum
  integer :: i

  !$acc data present (t,rho,dv,mu,sc,q,qvl,qvi,lcldm,precip_frac,arn,asn) &
  !$acc      present (qcic,qiic,qric,qsic,lamr,n0r,lams,n0s,pre,prds,am_evp_st) &
  !$acc      create  (dum)

  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     am_evp_st(i) = 0._r8
     if (qcic(i)+qiic(i) < 1.e-6_r8) then
        dum(i) = 0._r8
     else
        dum(i) = lcldm(i)
     end if
  end do

  !$acc loop gang vector private(qclr,eps,abr,ab)
  do i=1,vlen
     ! only calculate if there is some precip fraction > cloud fraction
     if (precip_frac(i) > dum(i)) then
        if (qric(i) >= qsmall .or. qsic(i) >= qsmall) then
           am_evp_st(i) = precip_frac(i) - dum(i)
           ! calculate q for out-of-cloud region
           qclr=(q(i)-dum(i)*qvl(i))/(1._r8-dum(i))
        end if

        ! evaporation of rain
        if (qric(i) >= qsmall) then
           if (.not.evap_rhthrsh.or.(evap_rhthrsh.and.qclr/qvl(i).le.0.9_r8)) then
              call calc_ab(t(i), qvl(i), xxlv, abr)
              eps = 2._r8*pi*n0r(i)*rho(i)*Dv(i)* &
                    (f1r/(lamr(i)*lamr(i))+ &
                    f2r*(arn(i)*rho(i)/mu(i))**0.5_r8* &
                    sc(i)**(1._r8/3._r8)*gamma_half_br_plus5/ &
                    (lamr(i)**(5._r8/2._r8+br/2._r8)))
              pre(i) = eps*(qclr-qvl(i))/abr
              ! only evaporate in out-of-cloud region
              ! and distribute across precip_frac
              pre(i)=min(pre(i)*am_evp_st(i),0._r8)
              pre(i)=pre(i)/precip_frac(i)
           else
              pre(i) = 0._r8
           end if
        else
           pre(i) = 0._r8
        end if

        ! sublimation of snow
        if (qsic(i) >= qsmall) then
           call calc_ab(t(i), qvi(i), xxls, ab)
           eps = 2._r8*pi*n0s(i)*rho(i)*Dv(i)* &
                (f1s/(lams(i)*lams(i))+ &
                f2s*(asn(i)*rho(i)/mu(i))**0.5_r8* &
                sc(i)**(1._r8/3._r8)*gamma_half_bs_plus5/ &
                (lams(i)**(5._r8/2._r8+bs/2._r8)))
           prds(i) = eps*(qclr-qvi(i))/ab
           ! only sublimate in out-of-cloud region and distribute over precip_frac
           prds(i) = min(prds(i)*am_evp_st(i),0._r8)
           prds(i) = prds(i)/precip_frac(i)
        else
           prds(i) = 0._r8
        end if
     else
        prds(i) = 0._r8
        pre(i)  = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine evaporate_sublimate_precip

! calculate evaporation/sublimation of rain (no snow in mg4)
!===================================================================
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

subroutine evaporate_sublimate_precip_mg4(t, rho, dv, mu, sc, q, qvl, qvi, &
     lcldm, precip_frac, arn, qcic, qiic, qric, lamr, n0r, &
     pre, am_evp_st, vlen, evap_rhthrsh)

  integer,  intent(in)                   :: vlen
  real(r8), dimension(vlen), intent(in)  :: t           ! temperature
  real(r8), dimension(vlen), intent(in)  :: rho         ! air density
  real(r8), dimension(vlen), intent(in)  :: dv          ! water vapor diffusivity
  real(r8), dimension(vlen), intent(in)  :: mu          ! viscosity
  real(r8), dimension(vlen), intent(in)  :: sc          ! schmidt number
  real(r8), dimension(vlen), intent(in)  :: q           ! humidity
  real(r8), dimension(vlen), intent(in)  :: qvl         ! saturation humidity (water)
  real(r8), dimension(vlen), intent(in)  :: qvi         ! saturation humidity (ice)
  real(r8), dimension(vlen), intent(in)  :: lcldm       ! liquid cloud fraction
  real(r8), dimension(vlen), intent(in)  :: precip_frac ! precipitation fraction (maximum overlap)
  ! fallspeed parameters
  real(r8), dimension(vlen), intent(in)  :: arn         ! rain
  ! In-cloud MMRs
  real(r8), dimension(vlen), intent(in)  :: qcic        ! cloud liquid
  real(r8), dimension(vlen), intent(in)  :: qiic        ! cloud ice
  real(r8), dimension(vlen), intent(in)  :: qric        ! rain
  ! Size parameters
  ! rain
  real(r8), dimension(vlen), intent(in)  :: lamr
  real(r8), dimension(vlen), intent(in)  :: n0r
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: pre
  real(r8), dimension(vlen), intent(out) :: am_evp_st   ! Fractional area where rain evaporates.
  ! Switch
  logical, intent(in)                    :: evap_rhthrsh 
  real(r8) :: qclr   ! water vapor mixing ratio in clear air
  real(r8) :: abr    ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale
  real(r8), dimension(vlen) :: dum
  integer :: i

  !$acc data present (t,rho,dv,mu,sc,q,qvl,qvi,lcldm,precip_frac,arn) &
  !$acc      present (qcic,qiic,qric,lamr,n0r,pre,am_evp_st) &
  !$acc      create  (dum)

  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     am_evp_st(i) = 0._r8
     if (qcic(i)+qiic(i) < 1.e-6_r8) then
        dum(i) = 0._r8
     else
        dum(i) = lcldm(i)
     end if
  end do

  !$acc loop gang vector private(qclr,eps,abr)
  do i=1,vlen
     ! only calculate if there is some precip fraction > cloud fraction
     if (precip_frac(i) > dum(i)) then
        if (qric(i) >= qsmall) then
           am_evp_st(i) = precip_frac(i) - dum(i)
           ! calculate q for out-of-cloud region
           qclr=(q(i)-dum(i)*qvl(i))/(1._r8-dum(i))
        end if

        ! evaporation of rain
        if (qric(i) >= qsmall) then
           if (.not.evap_rhthrsh.or.(evap_rhthrsh.and.qclr/qvl(i).le.0.9_r8)) then
              call calc_ab(t(i), qvl(i), xxlv, abr)
              eps = 2._r8*pi*n0r(i)*rho(i)*Dv(i)* &
                 (f1r/(lamr(i)*lamr(i))+ &
                 f2r*(arn(i)*rho(i)/mu(i))**0.5_r8* &
                 sc(i)**(1._r8/3._r8)*gamma_half_br_plus5/ &
                 (lamr(i)**(5._r8/2._r8+br/2._r8)))
              pre(i) = eps*(qclr-qvl(i))/abr
              ! only evaporate in out-of-cloud region
              ! and distribute across precip_frac
              pre(i) = min(pre(i)*am_evp_st(i),0._r8)
              pre(i) = pre(i)/precip_frac(i)
           else
              pre(i) = 0._r8
           end if
        else
           pre(i) = 0._r8
        end if
     else
        pre(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine evaporate_sublimate_precip_mg4

! evaporation/sublimation of rain, snow and graupel
!===================================================================
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

subroutine evaporate_sublimate_precip_graupel(t, rho, dv, mu, sc, q, qvl, qvi, &
     lcldm, precip_frac, arn, asn, agn, bg, qcic, qiic, qric, qsic, qgic, lamr, n0r, lams, n0s, lamg, n0g, &
     pre, prds, prdg, am_evp_st, vlen, evap_rhthrsh)

  integer,                   intent(in)  :: vlen
  real(r8), dimension(vlen), intent(in)  :: t           ! temperature
  real(r8), dimension(vlen), intent(in)  :: rho         ! air density
  real(r8), dimension(vlen), intent(in)  :: dv          ! water vapor diffusivity
  real(r8), dimension(vlen), intent(in)  :: mu          ! viscosity
  real(r8), dimension(vlen), intent(in)  :: sc          ! schmidt number
  real(r8), dimension(vlen), intent(in)  :: q           ! humidity
  real(r8), dimension(vlen), intent(in)  :: qvl         ! saturation humidity (water)
  real(r8), dimension(vlen), intent(in)  :: qvi         ! saturation humidity (ice)
  real(r8), dimension(vlen), intent(in)  :: lcldm       ! liquid cloud fraction
  real(r8), dimension(vlen), intent(in)  :: precip_frac ! precipitation fraction (maximum overlap)
  ! fallspeed parameters
  real(r8), dimension(vlen), intent(in)  :: arn         ! rain
  real(r8), dimension(vlen), intent(in)  :: asn         ! snow
  real(r8), dimension(vlen), intent(in)  :: agn         ! graupel
  real(r8),                  intent(in)  :: bg 
  ! In-cloud MMRs
  real(r8), dimension(vlen), intent(in)  :: qcic        ! cloud liquid
  real(r8), dimension(vlen), intent(in)  :: qiic        ! cloud ice
  real(r8), dimension(vlen), intent(in)  :: qric        ! rain
  real(r8), dimension(vlen), intent(in)  :: qsic        ! snow
  real(r8), dimension(vlen), intent(in)  :: qgic        ! graupel
  ! Size parameters
  ! rain
  real(r8), dimension(vlen), intent(in)  :: lamr
  real(r8), dimension(vlen), intent(in)  :: n0r
  ! snow
  real(r8), dimension(vlen), intent(in)  :: lams
  real(r8), dimension(vlen), intent(in)  :: n0s
  ! graupel
  real(r8), dimension(vlen), intent(in)  :: lamg
  real(r8), dimension(vlen), intent(in)  :: n0g
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: pre
  real(r8), dimension(vlen), intent(out) :: prds
  real(r8), dimension(vlen), intent(out) :: prdg
  real(r8), dimension(vlen), intent(out) :: am_evp_st   ! Fractional area where rain evaporates.
  ! Switch
  logical, intent(in) :: evap_rhthrsh 

  real(r8) :: qclr   ! water vapor mixing ratio in clear air
  real(r8) :: ab, abr, abg    ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale
  real(r8), dimension(vlen) :: dum
  integer :: i

  !$acc data present (t,rho,dv,mu,sc,q,qvl,qvi,lcldm,precip_frac) &
  !$acc      present (arn,asn,agn,qcic,qiic,qric,qsic,qgic,lamr)  &
  !$acc      present (n0r,lams,n0s,lamg,n0g,pre,prds,prdg,am_evp_st) &
  !$acc      create  (dum)

  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     am_evp_st(i) = 0._r8
     if (qcic(i)+qiic(i) < 1.e-6_r8) then
        dum(i) = 0._r8
     else
        dum(i) = lcldm(i)
     end if
  end do

  !$acc loop gang vector private(qclr,eps,abr,ab,abg)
  do i=1,vlen
     ! only calculate if there is some precip fraction > cloud fraction
     if (precip_frac(i) > dum(i)) then
        if (qric(i) >= qsmall .or. qsic(i) >= qsmall .or. qgic(i) >= qsmall) then
           am_evp_st(i) = precip_frac(i) - dum(i)
           ! calculate q for out-of-cloud region
           qclr=(q(i)-dum(i)*qvl(i))/(1._r8-dum(i))
        end if

        ! evaporation of rain
        if (qric(i) >= qsmall) then
           if (.not.evap_rhthrsh.or.(evap_rhthrsh.and.qclr/qvl(i).le.0.9_r8)) then
              call calc_ab(t(i), qvl(i), xxlv, abr)
              eps = 2._r8*pi*n0r(i)*rho(i)*Dv(i)* &
                 (f1r/(lamr(i)*lamr(i))+ &
                 f2r*(arn(i)*rho(i)/mu(i))**0.5_r8* &
                 sc(i)**(1._r8/3._r8)*gamma_half_br_plus5/ &
                 (lamr(i)**(5._r8/2._r8+br/2._r8)))
              pre(i) = eps*(qclr-qvl(i))/abr
              ! only evaporate in out-of-cloud region
              ! and distribute across precip_frac
              pre(i) = min(pre(i)*am_evp_st(i),0._r8)
              pre(i) = pre(i)/precip_frac(i)
           else
              pre(i) = 0._r8
           end if
        else
           pre(i) = 0._r8
        end if

        ! sublimation of snow
        if (qsic(i) >= qsmall) then
           call calc_ab(t(i), qvi(i), xxls, ab)
           eps = 2._r8*pi*n0s(i)*rho(i)*Dv(i)* &
                (f1s/(lams(i)*lams(i))+ &
                f2s*(asn(i)*rho(i)/mu(i))**0.5_r8* &
                sc(i)**(1._r8/3._r8)*gamma_half_bs_plus5/ &
                (lams(i)**(5._r8/2._r8+bs/2._r8)))
           prds(i) = eps*(qclr-qvi(i))/ab
           ! only sublimate in out-of-cloud region and distribute over precip_frac
           prds(i) = min(prds(i)*am_evp_st(i),0._r8)
           prds(i) = prds(i)/precip_frac(i)
        else
           prds(i) = 0._r8
        end if

        ! ADD GRAUPEL, do Same with prdg.
        if (qgic(i).ge.qsmall) then
           call calc_ab(t(i), qvi(i), xxls, abg)
           eps = 2._r8*pi*n0g(i)*rho(i)*Dv(i)*                    &
                (f1s/(lamg(i)*lamg(i))+                           &
                f2s*(agn(i)*rho(i)/mu(i))**0.5_r8*                &
                sc(i)**(1._r8/3._r8)*gamma(5._r8/2._r8+bg/2._r8)/ &
                (lamg(i)**(5._r8/2._r8+bs/2._r8)))
           prdg(i) = eps*(qclr-qvi(i))/abg
           ! only sublimate in out-of-cloud region and distribute over precip_frac
           prdg(i) = min(prdg(i)*am_evp_st(i),0._r8)
           prdg(i) = prdg(i)/precip_frac(i)
        else
           prdg(i) = 0._r8
        end if

     else
        prds(i) = 0._r8
        pre(i)  = 0._r8
        prdg(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine evaporate_sublimate_precip_graupel

! evaporation/sublimation of rain and graupel (no snow in mg4)
!===================================================================
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and graupel is neglected
! except for transfer of cloud water to snow through bergeron process

subroutine evaporate_sublimate_precip_graupel_mg4(t, rho, dv, mu, sc, q, qvl, qvi, &
     lcldm, precip_frac, arn, agn, bg, qcic, qiic, qric, qgic, lamr, n0r, lamg, n0g, &
     pre, prdg, am_evp_st, vlen)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t            ! temperature
  real(r8), dimension(vlen), intent(in) :: rho          ! air density
  real(r8), dimension(vlen), intent(in) :: dv           ! water vapor diffusivity
  real(r8), dimension(vlen), intent(in) :: mu           ! viscosity
  real(r8), dimension(vlen), intent(in) :: sc           ! schmidt number
  real(r8), dimension(vlen), intent(in) :: q            ! humidity
  real(r8), dimension(vlen), intent(in) :: qvl          ! saturation humidity (water)
  real(r8), dimension(vlen), intent(in) :: qvi          ! saturation humidity (ice)
  real(r8), dimension(vlen), intent(in) :: lcldm        ! liquid cloud fraction
  real(r8), dimension(vlen), intent(in) :: precip_frac  ! precipitation fraction (maximum overlap)
  ! fallspeed parameters
  real(r8), dimension(vlen), intent(in) :: arn          ! rain
  real(r8), dimension(vlen), intent(in) :: agn          ! graupel
  real(r8),                  intent(in) :: bg 
  ! In-cloud MMRs
  real(r8), dimension(vlen), intent(in) :: qcic         ! cloud liquid
  real(r8), dimension(vlen), intent(in) :: qiic         ! cloud ice
  real(r8), dimension(vlen), intent(in) :: qric         ! rain
  real(r8), dimension(vlen), intent(in) :: qgic         ! graupel
  ! Size parameters
  ! rain
  real(r8), dimension(vlen), intent(in) :: lamr
  real(r8), dimension(vlen), intent(in) :: n0r
  ! graupel
  real(r8), dimension(vlen), intent(in) :: lamg
  real(r8), dimension(vlen), intent(in) :: n0g
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: pre
  real(r8), dimension(vlen), intent(out) :: prdg
  real(r8), dimension(vlen), intent(out) :: am_evp_st   ! Fractional area where rain evaporates.

  real(r8) :: qclr                                      ! water vapor mixing ratio in clear air
  real(r8) :: abr, abg                                  ! correction to account for latent heat
  real(r8) :: eps                                       ! 1/ sat relaxation timescale
  real(r8), dimension(vlen) :: dum
  integer :: i

  !$acc data present (t,rho,dv,mu,sc,q,qvl,qvi,lcldm,precip_frac) &
  !$acc      present (arn,agn,qcic,qiic,qric,qgic,lamr,n0r,lamg)  &
  !$acc      present (n0g,pre,prdg,am_evp_st) &
  !$acc      create  (dum)

  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise
  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     am_evp_st(i) = 0._r8
     if (qcic(i)+qiic(i) < 1.e-6_r8) then
        dum(i) = 0._r8
     else
        dum(i) = lcldm(i)
     end if
  end do

  !$acc loop gang vector private(qclr,eps,abr,abg)
  do i=1,vlen
     ! only calculate if there is some precip fraction > cloud fraction
     if (precip_frac(i) > dum(i)) then
        if (qric(i) >= qsmall .or. qgic(i) >= qsmall) then
           am_evp_st(i) = precip_frac(i) - dum(i)
           ! calculate q for out-of-cloud region
           qclr = (q(i)-dum(i)*qvl(i))/(1._r8-dum(i))
        end if

        ! evaporation of rain
        if (qric(i) >= qsmall) then
           call calc_ab(t(i), qvl(i), xxlv, abr)
           eps = 2._r8*pi*n0r(i)*rho(i)*Dv(i)* &
                (f1r/(lamr(i)*lamr(i))+ &
                f2r*(arn(i)*rho(i)/mu(i))**0.5_r8* &
                sc(i)**(1._r8/3._r8)*gamma_half_br_plus5/ &
                (lamr(i)**(5._r8/2._r8+br/2._r8)))
           pre(i) = eps*(qclr-qvl(i))/abr
           ! only evaporate in out-of-cloud region
           ! and distribute across precip_frac
           pre(i) = min(pre(i)*am_evp_st(i),0._r8)
           pre(i) = pre(i)/precip_frac(i)
        else
           pre(i) = 0._r8
        end if

        !++AG ADD GRAUPEL, do Same with prdg.
        if (qgic(i).ge.qsmall) then
           call calc_ab(t(i), qvi(i), xxls, abg)
           eps = 2._r8*pi*n0g(i)*rho(i)*Dv(i)*                    &
                (f1s/(lamg(i)*lamg(i))+                           &
                f2s*(agn(i)*rho(i)/mu(i))**0.5_r8*                &
                sc(i)**(1._r8/3._r8)*gamma(5._r8/2._r8+bg/2._r8)/ &
                (lamg(i)**(5._r8/2._r8+bs/2._r8)))
           prdg(i) = eps*(qclr-qvi(i))/abg
           ! only sublimate in out-of-cloud region and distribute over precip_frac
           prdg(i) = min(prdg(i)*am_evp_st(i),0._r8)
           prdg(i) = prdg(i)/precip_frac(i)
        else
           prdg(i) = 0._r8
        end if
     else
        pre(i) = 0._r8
!++ag
        prdg(i) = 0._r8
!--ag
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine evaporate_sublimate_precip_graupel_mg4

! bergeron process - evaporation of droplets and deposition onto snow
!===================================================================

subroutine bergeron_process_snow(t, rho, dv, mu, sc, qvl, qvi, asn, &
     qcic, qsic, lams, n0s, bergs, vlen)

  integer,                   intent(in) :: vlen
  real(r8), dimension(vlen), intent(in) :: t    ! temperature
  real(r8), dimension(vlen), intent(in) :: rho  ! air density
  real(r8), dimension(vlen), intent(in) :: dv   ! water vapor diffusivity
  real(r8), dimension(vlen), intent(in) :: mu   ! viscosity
  real(r8), dimension(vlen), intent(in) :: sc   ! schmidt number
  real(r8), dimension(vlen), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), dimension(vlen), intent(in) :: qvi  ! saturation humidity (ice)
  ! fallspeed parameter for snow
  real(r8), dimension(vlen), intent(in) :: asn
  ! In-cloud MMRs
  real(r8), dimension(vlen), intent(in) :: qcic ! cloud liquid mixing ratio
  real(r8), dimension(vlen), intent(in) :: qsic ! snow mixing ratio
  ! Size parameters for snow
  real(r8), dimension(vlen), intent(in) :: lams
  real(r8), dimension(vlen), intent(in) :: n0s
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: bergs

  real(r8) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale
  integer :: i

  !$acc data present (t,rho,dv,mu,sc,qvl,qvi) &
  !$acc      present (asn,qcic,qsic,lams,n0s,bergs)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(eps,ab)
  do i=1,vlen
     if (qsic(i) >= qsmall.and. qcic(i) >= qsmall .and. t(i) < tmelt) then
        call calc_ab(t(i), qvi(i), xxls, ab)
        eps = 2._r8*pi*n0s(i)*rho(i)*Dv(i)* &
             (f1s/(lams(i)*lams(i))+ &
             f2s*(asn(i)*rho(i)/mu(i))**0.5_r8* &
             sc(i)**(1._r8/3._r8)*gamma_half_bs_plus5/ &
             (lams(i)**(5._r8/2._r8+bs/2._r8)))
        bergs(i) = eps*(qvl(i)-qvi(i))/ab
     else
        bergs(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine bergeron_process_snow

!========================================================================
! Collection of snow by rain to form graupel
!========================================================================

subroutine graupel_collecting_snow(qsic,qric,umr,ums,rho,lamr,n0r,lams,n0s, &
     psacr, vlen)
  
  integer,                   intent(in)  :: vlen
  ! In-cloud MMRs
  real(r8), dimension(vlen), intent(in)  :: qsic  ! snow
  real(r8), dimension(vlen), intent(in)  :: qric  ! rain
  ! mass-weighted fall speeds
  real(r8), dimension(vlen), intent(in)  :: umr   ! rain
  real(r8), dimension(vlen), intent(in)  :: ums   ! snow
  real(r8), dimension(vlen), intent(in)  :: rho   ! air density
  ! Size parameters for rain
  real(r8), dimension(vlen), intent(in)  :: lamr
  real(r8), dimension(vlen), intent(in)  :: n0r
  ! Size parameters for snow
  real(r8), dimension(vlen), intent(in)  :: lams
  real(r8), dimension(vlen), intent(in)  :: n0s
  real(r8), dimension(vlen), intent(out) :: psacr ! conversion due to coll of snow by rain

  real(r8), parameter :: cons31 = pi*pi*ecr*rhosn
  integer :: i

  !$acc data present (qsic,qric,umr,ums,rho,lamr) &
  !$acc      present (n0r,lams,n0s,psacr)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if (qsic(i).ge.0.1e-3_r8 .and. qric(i).ge.0.1e-3_r8) then
        psacr(i) = cons31*(((1.2_r8*umr(i)-0.95_r8*ums(i))**2+ &
             0.08_r8*ums(i)*umr(i))**0.5_r8*rho(i)*            &
             n0r(i)*n0s(i)/lams(i)**3*                         &
             (5._r8/(lams(i)**3*lamr(i))+                      &
             2._r8/(lams(i)**2*lamr(i)**2)+                    &
             0.5_r8/(lams(i)*lamr(i)**3)))            
     else
        psacr(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine graupel_collecting_snow

!========================================================================
! Collection of cloud water by graupel
!========================================================================

subroutine graupel_collecting_cld_water(qgic,qcic,ncic,rho,n0g,lamg,bg,agn, &
     psacwg, npsacwg, vlen)

  integer,                   intent(in) :: vlen
  ! In-cloud MMRs
  real(r8), dimension(vlen), intent(in) :: qgic ! graupel
  real(r8), dimension(vlen), intent(in) :: qcic ! cloud water
  real(r8), dimension(vlen), intent(in) :: ncic ! cloud water number conc
  real(r8), dimension(vlen), intent(in) :: rho  ! air density
  ! Size parameters for graupel
  real(r8), dimension(vlen), intent(in) :: lamg
  real(r8), dimension(vlen), intent(in) :: n0g
  ! fallspeed parameters for graupel
  ! Set AGN and BG  as input (in micro_mg3_0.F90)  (select hail or graupel as appropriate)
  real(r8),                  intent(in) :: bg
  real(r8), dimension(vlen), intent(in) :: agn
  ! Output tendencies
  real(r8), dimension(vlen), intent(out) :: psacwg
  real(r8), dimension(vlen), intent(out) :: npsacwg

  real(r8) :: cons
  integer  :: i 

  cons = gamma(bg + 3._r8)*pi/4._r8 * ecid

  !$acc data present (qgic,qcic,ncic,rho,lamg,n0g,agn,psacwg,npsacwg)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector
  do i=1,vlen
     if (qgic(i).ge.1.e-8_r8 .and. qcic(i).ge.qsmall) then
        psacwg(i) = cons*agn(i)*qcic(i)*rho(i)*  &
               n0g(i)/                           &
               lamg(i)**(bg+3._r8)
        npsacwg(i) = cons*agn(i)*ncic(i)*rho(i)* &
               n0g(i)/                           &
               lamg(i)**(bg+3._r8)
     else
        psacwg(i)  = 0._r8
        npsacwg(i) = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine graupel_collecting_cld_water

!========================================================================
! Conversion of rimed cloud water onto snow to graupel/hail
!========================================================================

subroutine graupel_riming_liquid_snow(psacws,qsic,qcic,nsic,rho,rhosn,rhog,asn,lams,n0s,dtime, &
     pgsacw,nscng,vlen)

  integer,                   intent(in)    :: vlen
  ! Accretion of cloud water to snow tendency (modified)
  real(r8), dimension(vlen), intent(inout) :: psacws
  real(r8), dimension(vlen), intent(in)    :: qsic    ! snow mixing ratio
  real(r8), dimension(vlen), intent(in)    :: qcic    ! cloud liquid mixing ratio
  real(r8), dimension(vlen), intent(in)    :: nsic    ! snow number concentration
  real(r8), dimension(vlen), intent(in)    :: rho     ! air density
  real(r8),                  intent(in)    :: rhosn   ! snow density
  real(r8),                  intent(in)    :: rhog    ! graupel density
  real(r8), dimension(vlen), intent(in)    :: asn     ! fall speed parameter for snow
  ! Size parameters for snow
  real(r8), dimension(vlen), intent(in)    :: lams
  real(r8), dimension(vlen), intent(in)    :: n0s
  real(r8),                  intent(in)    :: dtime
  !Output process rates
  real(r8), dimension(vlen), intent(out)   :: pgsacw  ! dQ graupel due to collection droplets by snow
  real(r8), dimension(vlen), intent(out)   :: nscng   ! dN graupel due to collection droplets by snow

  real(r8) :: cons
  real(r8) :: rhosu
  real(r8) :: dum
  integer  :: i 

!........................................................................
!Input: PSACWS,qs,qc,n0s,rho,lams,rhos,rhog
!Output:PSACWS,PGSACW,NSCNG

  rhosu = 85000._r8/(ra * tmelt)    ! typical air density at 850 mb

  cons  = 4._r8 *2._r8 *3._r8*rhosu*pi*ecid*ecid*gamma_2bs_plus2/(8._r8*(rhog-rhosn))

  !$acc data present (psacws,qsic,qcic,nsic,rho,asn,lams,n0s,pgsacw,nscng)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(dum)
  do i=1,vlen
     if (psacws(i).gt.0._r8 .and. qsic(i).GE.0.1e-3_r8 .AND. qcic(i).GE.0.5E-3_r8) then
        ! Only allow conversion if qni > 0.1 and qc > 0.5 g/kg following Rutledge and Hobbs (1984)
        ! if (qsic(i).GE.0.1e-3_r8 .AND. qcic(i).GE.0.5E-3_r8) then

        ! portion of riming converted to graupel (Reisner et al. 1998, originally IS1991)
        ! dtime here is correct.
        pgsacw(i) = min(psacws(i),cons*dtime*n0s(i)*qcic(i)*qcic(i)* &
             asn(i)*asn(i)/ &
             (rho(i)*lams(i)**(2._r8*bs+2._r8))) 

        ! Mix rat converted into graupel as embryo (Reisner et al. 1998, orig M1990)
        dum       = max(rhosn/(rhog-rhosn)*pgsacw(i),0._r8) 

        ! Number concentraiton of embryo graupel from riming of snow
        nscng(i)  = dum/mg0*rho(i)
        ! Limit max number converted to snow number  (dtime here correct? We think yes)
        nscng(i)  = min(nscng(i),nsic(i)/dtime)

        ! Portion of riming left for snow
        psacws(i) = psacws(i) - pgsacw(i)
     else
        pgsacw(i) = 0._r8
        nscng(i)  = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine graupel_riming_liquid_snow

!========================================================================
!CHANGE IN Q,N COLLECTION RAIN BY GRAUPEL
!========================================================================

subroutine graupel_collecting_rain(qric,qgic,umg,umr,ung,unr,rho,n0r,lamr,n0g,lamg,&
     pracg,npracg,vlen)

  integer,                   intent(in) :: vlen
  !MMR
  real(r8), dimension(vlen), intent(in) :: qric  !rain mixing ratio
  real(r8), dimension(vlen), intent(in) :: qgic  !graupel mixing ratio
  !Mass weighted Fall speeds
  real(r8), dimension(vlen), intent(in) :: umg ! graupel fall speed
  real(r8), dimension(vlen), intent(in) :: umr ! rain fall speed
  !Number weighted fall speeds
  real(r8), dimension(vlen), intent(in) :: ung ! graupel fall speed
  real(r8), dimension(vlen), intent(in) :: unr ! rain fall speed
  real(r8), dimension(vlen), intent(in) :: rho   ! air density
  ! Size parameters for rain
  real(r8), dimension(vlen), intent(in) :: n0r
  real(r8), dimension(vlen), intent(in) :: lamr
  ! Size parameters for graupel
  real(r8), dimension(vlen), intent(in) :: n0g
  real(r8), dimension(vlen), intent(in) :: lamg
  !Output process rates
  real(r8), dimension(vlen), intent(out) :: pracg   ! Q collection rain by graupel
  real(r8), dimension(vlen), intent(out) :: npracg  ! N collection rain by graupel

! Add collection of graupel by rain above freezing
! assume all rain collection by graupel above freezing is shed
! assume shed drops are 1 mm in size
  integer  :: i 
  real(r8), parameter :: cons41 = pi*pi*ecr*rhow
  real(r8), parameter :: cons32 = pi/2._r8*ecr
  real(r8) :: dum

  !$acc data present (qric,qgic,umg,umr,ung,unr,rho) &
  !$acc      present (n0r,lamr,n0g,lamg,pracg,npracg)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(dum)
  do i=1,vlen
     if (qric(i).ge.1.e-8_r8.and.qgic(i).ge.1.e-8_r8) then
        ! pracg is mixing ratio of rain per sec collected by graupel/hail
        pracg(i) = cons41*(((1.2_r8*umr(i)-0.95_r8*umg(i))**2._r8+  &
             0.08_r8*umg(i)*umr(i))**0.5_r8*rho(i)*                 &
             n0r(i)*n0g(i)/lamr(i)**3._r8*                          &
             (5._r8/(lamr(i)**3._r8*lamg(i))+                       &
             2._r8/(lamr(i)**2._r8*lamg(i)**2._r8)+                 &
             0.5_r8/(lamr(i)*lamg(i)**3._r8)))

! assume 1 mm drops are shed, get number shed per sec
        dum = pracg(i)/5.2e-7_r8
        npracg(i) = cons32*rho(i)*(1.7_r8*(unr(i)-ung(i))**2._r8+   &
             0.3_r8*unr(i)*ung(i))**0.5_r8*n0r(i)*n0g(i)*           &
             (1._r8/(lamr(i)**3._r8*lamg(i))+                       &
             1._r8/(lamr(i)**2._r8*lamg(i)**2._r8)+                 &
             1._r8/(lamr(i)*lamg(i)**3._r8))
        
! hm 7/15/13, remove limit so that the number of collected drops can smaller than 
! number of shed drops
!            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)
        npracg(i) = npracg(i) - dum
     else
        npracg(i) = 0._r8
        pracg(i)  = 0._r8
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine graupel_collecting_rain

!========================================================================
! Rain riming snow to graupel
!========================================================================
! Conversion of rimed rainwater onto snow converted to graupel

subroutine graupel_rain_riming_snow(pracs,npracs,psacr,qsic,qric,nric,nsic,n0s,lams,n0r,lamr,dtime,&
     pgracs,ngracs,vlen)

  integer,                   intent(in)    :: vlen
  ! Accretion of rain by snow
  real(r8), dimension(vlen), intent(inout) :: pracs
  real(r8), dimension(vlen), intent(inout) :: npracs
  real(r8), dimension(vlen), intent(inout) :: psacr  ! conversion due to coll of snow by rain
  !MMR
  real(r8), dimension(vlen), intent(in)    :: qsic   ! snow mixing ratio
  real(r8), dimension(vlen), intent(in)    :: qric   ! rain mixing ratio
  real(r8), dimension(vlen), intent(in)    :: nric   ! rain number concentration
  real(r8), dimension(vlen), intent(in)    :: nsic   ! snow number concentration
  ! Size parameters for snow
  real(r8), dimension(vlen), intent(in)    :: n0s
  real(r8), dimension(vlen), intent(in)    :: lams
  ! Size parameters for rain
  real(r8), dimension(vlen), intent(in)    :: n0r
  real(r8), dimension(vlen), intent(in)    :: lamr
  real(r8),                  intent(in)    :: dtime
  !Output process rates
  real(r8), dimension(vlen), intent(out)   :: pgracs  ! Q graupel due to collection rain by snow
  real(r8), dimension(vlen), intent(out)   :: ngracs  ! N graupel due to collection rain by snow

!Input: PRACS,NPRACS,PSACR,qni,qr,lams,lamr,nr,ns  Note: No PSACR in MG2
!Output:PGRACS,NGRACS,PRACS,PSACR
  integer  :: i 
  real(r8), parameter :: cons18 = rhosn*rhosn
  real(r8), parameter :: cons19 = rhow*rhow
  real(r8) :: dum

  !$acc data present (pracs,npracs,psacr,qsic,qric,nric)  &
  !$acc      present (nsic,n0s,lams,n0r,lamr,pgracs,ngracs)

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(dum)
  do i=1,vlen
     if (pracs(i).gt.0._r8.and.qsic(i).ge.0.1e-3_r8.and.qric(i).ge.0.1e-3_r8) then
        ! only allow conversion if qs > 0.1 and qr > 0.1 g/kg following rutledge and hobbs (1984)
        !if (qsic(i).ge.0.1e-3_r8.and.qric(i).ge.0.1e-3_r8) then
           ! portion of collected rainwater converted to graupel (reisner et al. 1998)
        dum = cons18*(4._r8/lams(i))**3*(4._r8/lams(i))**3 &    
             /(cons18*(4._r8/lams(i))**3*(4._r8/lams(i))**3+ &  
             cons19*(4._r8/lamr(i))**3*(4._r8/lamr(i))**3)
        dum = min(dum,1._r8)
        dum = max(dum,0._r8)
        pgracs(i) = (1._r8-dum)*pracs(i)
        ngracs(i) = (1._r8-dum)*npracs(i)
        ! limit max number converted to min of either rain or snow number concentration
        ngracs(i) = min(ngracs(i),nric(i)/dtime)
        ngracs(i) = min(ngracs(i),nsic(i)/dtime)
        
        ! amount left for snow production
        pracs(i)  = pracs(i) - pgracs(i)
        npracs(i) = npracs(i) - ngracs(i)
        
        ! conversion to graupel due to collection of snow by rain
        psacr(i)=psacr(i)*(1._r8-dum)
     else
        pgracs(i) = 0._r8
        ngracs(i) = 0._r8
     end if
  end do 
  !$acc end parallel

  !$acc end data
end subroutine graupel_rain_riming_snow

!========================================================================
! Rime Splintering
!========================================================================
subroutine graupel_rime_splintering(t,qcic,qric,qgic,psacwg,pracg,&
     qmultg,nmultg,qmultrg,nmultrg,vlen)

  integer,                   intent(in)    :: vlen
  real(r8), dimension(vlen), intent(in)    :: t        ! temperature
  !MMR
  real(r8), dimension(vlen), intent(in)    :: qcic     ! liquid mixing ratio
  real(r8), dimension(vlen), intent(in)    :: qric     ! rain mixing ratio
  real(r8), dimension(vlen), intent(in)    :: qgic     ! graupel mixing ratio
  ! Already calculated terms for collection 
  real(r8), dimension(vlen), intent(inout) :: psacwg   ! collection droplets by graupel
  real(r8), dimension(vlen), intent(inout) :: pracg    ! collection rain by graupel
  !Output process rates for splintering
  real(r8), dimension(vlen), intent(out)   :: qmultg   ! Q ice mult of droplets/graupel
  real(r8), dimension(vlen), intent(out)   :: nmultg   ! N ice mult of droplets/graupel
  real(r8), dimension(vlen), intent(out)   :: qmultrg  ! Q ice mult of rain/graupel
  real(r8), dimension(vlen), intent(out)   :: nmultrg  ! N ice mult of rain/graupel

!Input: qg,qc,qr, PSACWG,PRACG,T
!Output: NMULTG,QMULTG,NMULTRG,QMULTRG,PSACWG,PRACG
  integer  :: i 
  real(r8) :: fmult
  real(r8) :: tm_3, tm_5, tm_8
 
  tm_3 = tmelt - 3._r8
  tm_5 = tmelt - 5._r8
  tm_8 = tmelt - 8._r8

  !$acc data present (t,qcic,qric,qgic,psacwg,pracg) &
  !$acc      present (qmultg,nmultg,qmultrg,nmultrg)

!nmultg,qmultg                                                                             .
!========================================================================

  !$acc parallel vector_length(VLENS) default(present)
  !$acc loop gang vector private(fmult)
  do i=1,vlen
     nmultrg(i) = 0._r8
     qmultrg(i) = 0._r8
     nmultg(i)  = 0._r8
     qmultg(i)  = 0._r8
     if (qgic(i).ge.0.1e-3_r8) then
        if (qcic(i).ge.0.5e-3_r8.or.qric(i).ge.0.1e-3_r8) then
           if (psacwg(i).gt.0._r8.or.pracg(i).gt.0._r8) then
              if (t(i).lt.tm_3 .and. t(i).gt.tm_8) then
                 if (t(i).gt.tm_3) then
                    fmult = 0._r8
                 else if (t(i).le.tm_3.and.t(i).gt.tm_5)  then
                    fmult = (tm_3-t(i))/2._r8
                 else if (t(i).ge.tm_8.and.t(i).le.tm_5)   then
                    fmult = (t(i)-tm_8)/3._r8
                 else if (t(i).lt.tm_8) then
                    fmult = 0._r8
                 end if
! 1000 is to convert from kg to g  (Check if needed for MG)
! splintering from droplets accreted onto graupel
                 if (psacwg(i).gt.0._r8) then
                    nmultg(i) = 35.e4_r8*psacwg(i)*fmult*1000._r8
                    qmultg(i) = nmultg(i)*mmult
! constrain so that transfer of mass from graupel to ice cannot be more mass
! than was rimed onto graupel
                    qmultg(i) = min(qmultg(i),psacwg(i))
                    psacwg(i) = psacwg(i)-qmultg(i)
                 end if

!nmultrg,qmultrg                                                                             .
!========================================================================
! riming and splintering from accreted raindrops
! Factor of 1000. again (Check)
                 if (pracg(i).gt.0._r8) then
                    nmultrg(i) = 35.e4_r8*pracg(i)*fmult*1000._r8
                    qmultrg(i) = nmultrg(i)*mmult
! constrain so that transfer of mass from graupel to ice cannot be more mass
! than was rimed onto graupel
                    qmultrg(i) = min(qmultrg(i),pracg(i))
                    pracg(i)   = pracg(i)-qmultrg(i)
                 end if
              end if
           end if
        end if
     end if
  end do
  !$acc end parallel

  !$acc end data
end subroutine graupel_rime_splintering

!========================================================================
! Evaporation and sublimation of graupel  
!========================================================================

!MERGE WITH RAIN AND SNOW EVAP
!
!subroutine graupel_sublimate_evap(t,q,qgic,rho,n0g,lamg,qvi,dv,mu,sc,bg,agn,&
!     prdg,eprdg,vlen)
!
!  integer, intent(in) :: vlen
!  
!  real(r8), dimension(vlen), intent(in) :: t  !temperature
!  real(r8), dimension(vlen), intent(in) :: q  !specific humidity (mixing ratio)
!
!  !MMR
!  real(r8), dimension(vlen), intent(in) :: qgic  !graupel mixing ratio
!
!  real(r8), dimension(vlen), intent(in) :: rho   ! air density
!
! ! Size parameters for graupel
!  real(r8), dimension(vlen), intent(in) :: n0g
!  real(r8), dimension(vlen), intent(in) :: lamg
!
!  real(r8), dimension(vlen), intent(in) :: qvi  !saturation humidity (ice)
!
!  real(r8), dimension(vlen), intent(in) :: dv   ! water vapor diffusivity
!  real(r8), dimension(vlen), intent(in) :: mu   ! viscosity
!  real(r8), dimension(vlen), intent(in) :: sc   ! schmidt number
!
!  ! fallspeed parameters for graupel
!  ! Set AGN and BG  as input (in micro_mg3_0.F90)  (select hail or graupel as appropriate)
!  real(r8),                    intent(in) :: bg
!  real(r8), dimension(vlen), intent(in) :: agn
!
!  ! Output tendencies (sublimation or evaporation of graupel)
!  real(r8), dimension(vlen), intent(out) :: prdg
!  real(r8), dimension(vlen), intent(out) :: eprdg
!
!  real(r8) :: cons11,cons36
!  real(r8) :: abi
!  real(r8) :: epsg
!  integer :: i 
!
!  cons11=gamma(2.5_r8+bg/2._r8)  !bg will be different for graupel (bg) or hail (bh)
!  cons36=(2.5_r8+bg/2._r8)  
!
!  
!  do i=1,vlen
! 
!     abi = calc_ab(t(i), qvi(i), xxls)
!
!      if (qgic(i).ge.qsmall) then
!         epsg = 2._r8*pi*n0g(i)*rho(i)*dv(i)*                                &
!              (f1s/(lamg(i)*lamg(i))+                               &
!              f2s*(agn(i)*rho(i)/mu(i))**0.5_r8*                      &
!              sc(i)**(1._r8/3._r8)*cons11/                   &
!              (lamg(i)**cons36))
!      else
!         epsg = 0.
!      end if
!
!! vapor dpeosition on graupel
!      prdg(i) = epsg*(q(i)-qvi(i))/abi
!
!! make sure not pushed into ice supersat/subsat
!! put this in main mg3 code…..check for it…
!! formula from reisner 2 scheme  
!!
!!          dum = (qv3d(k)-qvi(k))/dt
!!
!!           fudgef = 0.9999
!!           sum_dep = prd(k)+prds(k)+mnuccd(k)+prdg(k)
!!
!!           if( (dum.gt.0. .and. sum_dep.gt.dum*fudgef) .or.                      &
!!               (dum.lt.0. .and. sum_dep.lt.dum*fudgef) ) then
!!	       prdg(k) = fudgef*prdg(k)*dum/sum_dep
!!           endif
!
!! if cloud ice/snow/graupel vap deposition is neg, then assign to sublimation processes
!
!      eprdg(i)=0._r8
!
!      if (prdg(i).lt.0._r8) then
!         eprdg(i)=prdg(i)
!         prdg(i)=0.
!      end if
!
!   enddo
!
!end subroutine graupel_sublimate_evap

!========================================================================
!MG4 Lookup Table
!========================================================================

SUBROUTINE init_lookup_table(lkuptable_filename)

  implicit none
    !character(100) :: lkuptable_filename = "/home/katec/mg4work/lookup_table.dat"
    !character(100) :: lkuptable_filename = "/glade/p/work/katec/mg3work/lookup_table.dat"
    !character(100) :: lkuptable_filename = "./lookup_table.dat"
    character(len=256), intent(in) :: lkuptable_filename 
    integer :: i,j,k,ii,jj,kk, tt
    real(r8)    :: dum

    ! ++ trude
    open(unit=10,file=lkuptable_filename,status='old')

    !------------------------------------------------------------------------------------------!
    ! read in ice microphysics table
     do tt = 1, tsize
       do i = 1,isize
          do k = 1,jsize
             read(10,*) dum,dum,itab(tt,i,k,1),itab(tt,i,k,2),             &
                  itab(tt,i,k,3),itab(tt,i,k,4),itab(tt,i,k,5),            &
                  itab(tt,i,k,6),itab(tt,i,k,7),itab(tt,i,k,8),            &
                  itab(tt,i,k,13),itab(tt,i,k,9),itab(tt,i,k,10),          &
                  itab(tt,i,k,11),itab(tt,i,k,12),itab(tt,i,k,14) 
            enddo
       enddo
       ! read in table for ice-rain collection
       do i = 1,isize
          do k = 1,jsize
             do j = 1,rcollsize
                read(10,*) dum,dum,dum,itabcoll(tt,i,k,j,1),                  &
                     itabcoll(tt,i,k,j,2),dum
                itabcoll(tt,i,k,j,1) = dlog10(itabcoll(tt,i,k,j,1))
                itabcoll(tt,i,k,j,2) = dlog10(itabcoll(tt,i,k,j,2))
             end do
          end do
       end do
    end do

    close(unit=10)

END SUBROUTINE init_lookup_table

pure SUBROUTINE access_lookup_table(dumii,dumi,dumk,index,dum1,dum2,dum4,proc)

  implicit none

  real(r8), intent(in) :: dum1,dum2,dum4
  real(r8), intent(out) :: proc
  real(r8) :: dproc1,dproc2,iproc1,gproc1,tmp1,tmp2
  integer, intent(in) :: dumi,dumk,dumii,index


  !(trude,  dumii could be the same as dumt)
  ! get value at current density index
  !***** without temperature dependence *******
  ! first interpolate for current rimed fraction index
  !   dproc1 = itab(dumi,dumk,index)+(dum1-real(dumi))*(itab(       &
  !            dumi+1,dumk,index)-itab(dumi,dumk,index))

  !   dproc2 = itab(dumi,dumk+1,index)+(dum1-real(dumi))*(itab(     &
  !          dumi+1,dumk+1,index)-itab(dumi,dumk+1,index))

  !   iproc1 = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

  !   proc    = iproc1
  !! ***** end without temperature dependence ***********
  ! ******************* from orig P3 ********************

  ! trude comment, dumjj is for density index, dumii is for rime, which is subsituted by temperature. All dumjj should be removed. 
  ! first interpolate for current temperature 
  dproc1 = itab(dumii,dumi,dumk,index)+(dum1-real(dumi))*(itab(dumii,       &
       dumi+1,dumk,index)-itab(dumii,dumi,dumk,index))

  dproc2 = itab(dumii,dumi,dumk+1,index)+(dum1-real(dumi))*(itab(dumii,     &
       dumi+1,dumk+1,index)-itab(dumii,dumi,dumk+1,index))

  iproc1 = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

  ! linearly interpolate to get process rates for temperature index + 1
  dproc1 = itab(dumii+1,dumi,dumk,index)+(dum1-real(dumi))*(itab(dumii+1,   &
       dumi+1,dumk,index)-itab(dumii+1,dumi,dumk,index))

  dproc2 = itab(dumii+1,dumi,dumk+1,index)+(dum1-real(dumi))*(itab(dumii+1, &
       dumi+1,dumk+1,index)-itab(dumii+1,dumi,dumk+1,index))

  gproc1 = dproc1+(dum2-real(dumk))*(dproc2-dproc1)
  tmp1   = iproc1+(dum4-real(dumii))*(gproc1-iproc1)
  proc   = tmp1

  ! ********** end from orig P3 **************

END SUBROUTINE access_lookup_table

!------------------------------------------------------------------------------------------!
pure SUBROUTINE access_lookup_table_coll(dumii,dumj,dumi,dumk,index,dum1,dum2,dum3,     &
     dum4,proc)

  implicit none

  real(r8),intent(in)    :: dum1,dum2,dum3,dum4
  real(r8), intent(out) :: proc
  real(r8) :: dproc1,dproc2,iproc1,gproc1,tmp1,tmp2,dproc11, &
       dproc12,dproc21,dproc22
  integer, intent(in) :: dumii,dumj,dumi,dumk,index


  ! This subroutine interpolates lookup table values for rain/ice collection processes
  !!******for without temperature dependeny *****
  !   dproc11 = itabcoll(dumi,dumk,dumj,index)+(dum1-real(dumi))*               &
  !             (itabcoll(dumi+1,dumk,dumj,index)-itabcoll(dumi,    &
  !             dumk,dumj,index))

  !   dproc21 = itabcoll(dumi,dumk+1,dumj,index)+(dum1-real(dumi))*             &
  !             (itabcoll(dumi+1,dumk+1,dumj,index)-itabcoll(dumi,  &
  !             dumk+1,dumj,index))

  !   dproc1  = dproc11+(dum2-real(dumk))*(dproc21-dproc11)

  !   dproc12 = itabcoll(dumi,dumk,dumj+1,index)+(dum1-real(dumi))*             &
  !             (itabcoll(dumi+1,dumk,dumj+1,index)-itabcoll(dumi,  &
  !             dumk,dumj+1,index))

  !   dproc22 = itabcoll(dumi,dumk+1,dumj+1,index)+(dum1-real(dumi))*           &
  !             (itabcoll(dumi+1,dumk+1,dumj+1,index)-itabcoll(     &
  !             dumi,dumk+1,dumj+1,index))

  !   dproc2  = dproc12+(dum2-real(dumk))*(dproc22-dproc12)
  !   iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

  !   proc    = iproc1 
  !!******end for without temperature dependeny *****

  ! ************ from original P3 ************
  ! current density index

  ! current temp index
  dproc11 = itabcoll(dumii,dumi,dumk,dumj,index)+(dum1-real(dumi))*               &
       (itabcoll(dumii,dumi+1,dumk,dumj,index)-itabcoll(dumii,dumi,    &
       dumk,dumj,index))

  dproc21 = itabcoll(dumii,dumi,dumk+1,dumj,index)+(dum1-real(dumi))*             &
       (itabcoll(dumii,dumi+1,dumk+1,dumj,index)-itabcoll(dumii,dumi,  &
       dumk+1,dumj,index))

  dproc1  = dproc11+(dum2-real(dumk))*(dproc21-dproc11)
  dproc12 = itabcoll(dumii,dumi,dumk,dumj+1,index)+(dum1-real(dumi))*             &
       (itabcoll(dumii,dumi+1,dumk,dumj+1,index)-itabcoll(dumii,dumi,  &
       dumk,dumj+1,index))

  dproc22 = itabcoll(dumii,dumi,dumk+1,dumj+1,index)+(dum1-real(dumi))*           &
       (itabcoll(dumii,dumi+1,dumk+1,dumj+1,index)-itabcoll(dumii,     &
       dumi,dumk+1,dumj+1,index))

  dproc2  = dproc12+(dum2-real(dumk))*(dproc22-dproc12)
  iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

  ! temp index + 1

  dproc11 = itabcoll(dumii+1,dumi,dumk,dumj,index)+(dum1-real(dumi))*             &
       (itabcoll(dumii+1,dumi+1,dumk,dumj,index)-itabcoll(dumii+1,     &
       dumi,dumk,dumj,index))

  dproc21 = itabcoll(dumii+1,dumi,dumk+1,dumj,index)+(dum1-real(dumi))*           &
       (itabcoll(dumii+1,dumi+1,dumk+1,dumj,index)-itabcoll(dumii+1,   &
       dumi,dumk+1,dumj,index))

  dproc1  = dproc11+(dum2-real(dumk))*(dproc21-dproc11)

  dproc12 = itabcoll(dumii+1,dumi,dumk,dumj+1,index)+(dum1-real(dumi))*           &
       (itabcoll(dumii+1,dumi+1,dumk,dumj+1,index)-itabcoll(dumii+1,   &
       dumi,dumk,dumj+1,index))

  dproc22 = itabcoll(dumii+1,dumi,dumk+1,dumj+1,index)+(dum1-real(dumi))*         &
       (itabcoll(dumii+1,dumi+1,dumk+1,dumj+1,index)-itabcoll(dumii+1, &
       dumi,dumk+1,dumj+1,index))

  dproc2  = dproc12+(dum2-real(dumk))*(dproc22-dproc12)

  gproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)
  tmp1    = iproc1+(dum4-real(dumii))*(gproc1-iproc1)
  proc    = tmp1

  ! ********** end orig P3 *******
END SUBROUTINE access_lookup_table_coll

!========================================================================
!UTILITIES
!========================================================================

pure function no_limiter()
  real(r8) :: no_limiter

  no_limiter = transfer(limiter_off, no_limiter)

end function no_limiter

pure function limiter_is_on(lim)
  !$acc routine seq

  real(r8), intent(in) :: lim
  logical :: limiter_is_on

  limiter_is_on = transfer(lim, limiter_off) /= limiter_off

end function limiter_is_on

end module micro_mg_utils

module micro_mg4_0
!---------------------------------------------------------------------------------
! Purpose:
!   MG microphysics version 4.0 - Update of MG microphysics with
!                                 the Eidhammer et al. ice microphysics scheme
!
! Author: Andrew Gettelman, Hugh Morrison
!
! Version 3 history: Sep 2016: development begun for hail, graupel 
!
! Version 2 history: Sep 2011: Development begun.
!                    Feb 2013: Added of prognostic precipitation.
!                    Aug 2015: Published and released version
! Contributions from:  Sean Santos, Peter Caldwell, Xiaohong Liu and Steve Ghan
!
! invoked in CAM by specifying -microphys=mg4
!
! References: 
!
!           Gettelman, A. and H. Morrison, Advanced Two-Moment Microphysics for Global Models. 
!
!           Part I: Off line tests and comparisons with other schemes. 
!
!           J. Climate, 28, 1268-1287. doi: 10.1175/JCLI-D-14-00102.1, 2015. 
!
!
!
!           Gettelman, A., H. Morrison, S. Santos, P. Bogenschutz and P. H. Caldwell 
!
!           Advanced Two-Moment Microphysics for Global Models. 
!
!           Part II: Global model solutions and Aerosol-Cloud Interactions. 
!
!           J. Climate, 28, 1288-1307. doi:10.1175/JCLI-D-14-00103.1 , 2015. 
!
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------
!
! NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
! microphysics in cooperation with the MG liquid microphysics. This is
! controlled by the do_cldice variable.
!
! If do_cldice is false, then MG microphysics should not update CLDICE or
! NUMICE; it is assumed that the other microphysics scheme will have updated
! CLDICE and NUMICE. The other microphysics should handle the following
! processes that would have been done by MG:
!   - Detrainment (liquid and ice)
!   - Homogeneous ice nucleation
!   - Heterogeneous ice nucleation
!   - Bergeron process
!   - Melting of ice
!   - Freezing of cloud drops
!   - Autoconversion (ice -> snow)
!   - Growth/Sublimation of ice
!   - Sedimentation of ice
!
! This option has not been updated since the introduction of prognostic
! precipitation, and probably should be adjusted to cover snow as well.
!
!---------------------------------------------------------------------------------
! Version 4.O based on micro_mg3_0.F90 
!---------------------------------------------------------------------------------
! Based on micro_mg (restructuring of former cldwat2m_micro)
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------
! Code comments added by HM, 093011
! General code structure:
!
! Code is divided into two main subroutines:
!   subroutine micro_mg_init --> initializes microphysics routine, should be called
!                                  once at start of simulation
!   subroutine micro_mg_tend --> main microphysics routine to be called each time step
!                                this also calls several smaller subroutines to calculate
!                                microphysical processes and other utilities
!
! List of external functions:
!   qsat_water --> for calculating saturation vapor pressure with respect to liquid water
!   qsat_ice --> for calculating saturation vapor pressure with respect to ice
!   gamma   --> standard mathematical gamma function
! .........................................................................
! List of inputs through use statement in fortran90:
! Variable Name                      Description                Units
! .........................................................................
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                  J kg-1 K-1
! tmelt           temperature of melting point for water          K
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! rh2o            gas constant for water vapor                  J kg-1 K-1
! latvap          latent heat of vaporization                   J kg-1
! latice          latent heat of fusion                         J kg-1
! qsat_water      external function for calculating liquid water
!                 saturation vapor pressure/humidity              -
! qsat_ice        external function for calculating ice
!                 saturation vapor pressure/humidity              pa
! rhmini          relative humidity threshold parameter for
!                 nucleating ice                                  -
! .........................................................................
! NOTE: List of all inputs/outputs passed through the call/subroutine statement
!       for micro_mg_tend is given below at the start of subroutine micro_mg_tend.
!---------------------------------------------------------------------------------

! Procedures required:
! 1) An implementation of the gamma function (if not intrinsic).
! 2) saturation vapor pressure and specific humidity over water
! 3) svp over ice

#ifndef HAVE_GAMMA_INTRINSICS
use shr_spfn_mod, only: gamma => shr_spfn_gamma
#endif

!++ag DEBUG: remove later
use cam_logfile,    only: iulog
!--ag

use wv_sat_methods, only: &
     qsat_water => wv_sat_qsat_water, &
     qsat_ice => wv_sat_qsat_ice

! Parameters from the utilities module.
use micro_mg_utils, only: &
     r8, &
     pi, &
     omsm, &
     qsmall, &
     icsmall, &
     mincld, &
     rhosn, &
     rhoi, &
     rhow, &
     rhows, &
     ac, bc, &
     ai, bi, &
     aj, bj, &
     ar, br, &
     as, bs, &
!++ag
     ag, bg, &
     ah, bh, &
     rhog,rhoh, &
!--ag
     mi0, &
     rising_factorial

implicit none
private
save

public :: &
     micro_mg_init, &
     micro_mg_get_cols, &
     micro_mg_tend_pre_adjust, &
     micro_mg_tend

! Switches for specification rather than prediction of droplet and crystal number
! note: number will be adjusted as needed to keep mean size within bounds,
! even when specified droplet or ice number is used
!
! If constant cloud ice number is set (nicons = .true.),
! then all microphysical processes except mass transfer due to ice nucleation
! (mnuccd) are based on the fixed cloud ice number. Calculation of
! mnuccd follows from the prognosed ice crystal number ni.

logical :: nccons ! nccons = .true. to specify constant cloud droplet number
logical :: nicons ! nicons = .true. to specify constant cloud ice number
!++ag kt
logical :: ngcons = .false. ! ngcons = .true. to specify constant graupel number
!--ag kt
! specified ice and droplet number concentrations
! note: these are local in-clouse values, not grid-mean
real(r8) :: ncnst    ! droplet num concentration when nccons=.true. (m-3)
real(r8) :: ninst    ! ice num concentration when nicons=.true. (m-3)
!++ag
real(r8) :: ngnst = 0.1e6_r8     ! graupel num concentration when ngcons=.true. (m-3)
!--ag

!=========================================================
! Private module parameters
!=========================================================

!Range of cloudsat reflectivities (dBz) for analytic simulator
real(r8), parameter :: csmin = -30._r8
real(r8), parameter :: csmax = 26._r8
real(r8), parameter :: mindbz = -99._r8
real(r8), parameter :: minrefl = 1.26e-10_r8    ! minrefl = 10._r8**(mindbz/10._r8)

! autoconversion size threshold for cloud ice to snow (m)
real(r8) :: dcs

! minimum mass of new crystal due to freezing of cloud droplets done
! externally (kg)
real(r8), parameter :: mi0l_min = 4._r8/3._r8*pi*rhow*(4.e-6_r8)**3

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_mg_init
real(r8) :: g           ! gravity
real(r8) :: r           ! dry air gas constant
real(r8) :: rv          ! water vapor gas constant
real(r8) :: cpp         ! specific heat of dry air
real(r8) :: tmelt       ! freezing point of water (K)

! latent heats of:
real(r8) :: xxlv        ! vaporization
real(r8) :: xlf         ! freezing
real(r8) :: xxls        ! sublimation

real(r8) :: rhmini      ! Minimum rh for ice cloud fraction > 0.

! flags
logical :: microp_uniform
logical :: do_cldice
logical :: use_hetfrz_classnuc
!++ag
logical :: do_hail
logical :: do_graupel
!--ag

real(r8) :: rhosu       ! typical 850mn air density

real(r8) :: icenuct     ! ice nucleation temperature: currently -5 degrees C

real(r8) :: snowmelt    ! what temp to melt all snow: currently 2 degrees C
real(r8) :: rainfrze    ! what temp to freeze all rain: currently -5 degrees C

! additional constants to help speed up code
real(r8) :: gamma_br_plus1
real(r8) :: gamma_br_plus4
real(r8) :: gamma_bs_plus1
real(r8) :: gamma_bs_plus4
real(r8) :: gamma_bi_plus1
real(r8) :: gamma_bi_plus4
real(r8) :: gamma_bj_plus1
real(r8) :: gamma_bj_plus4
real(r8) :: xxlv_squared
real(r8) :: xxls_squared

character(len=16)  :: micro_mg_precip_frac_method  ! type of precipitation fraction method
real(r8)           :: micro_mg_berg_eff_factor     ! berg efficiency factor
! ++ Trude
real(r8)           :: micro_mg_accre_enhan_fact     ! accretion enhancment factor
real(r8)           :: micro_mg_autocon_fact     ! autoconversion prefactor
real(r8)           :: micro_mg_autocon_exp     ! autoconversion exponent factor
real(r8)           :: micro_mg_homog_size  ! size of freezing homogeneous ice
real(r8)           :: micro_mg_vtrmi_factor
real(r8)           :: micro_mg_effi_factor
real(r8)           :: micro_mg_iaccr_factor
real(r8)           :: micro_mg_max_nicons
! -- Trude
logical  :: allow_sed_supersat ! Allow supersaturated conditions after sedimentation loop
logical  :: do_sb_physics ! do SB 2001 autoconversion or accretion physics

!===============================================================================
contains
!===============================================================================

subroutine micro_mg_init( &
     kind, gravit, rair, rh2o, cpair,    &
     tmelt_in, latvap, latice,           &
     rhmini_in, micro_mg_dcs,            &
!++ag
     micro_mg_do_hail_in,micro_mg_do_graupel_in, &
!--ag
     microp_uniform_in, do_cldice_in, use_hetfrz_classnuc_in, &
     micro_mg_precip_frac_method_in, micro_mg_berg_eff_factor_in, &
! ++ trude
     micro_mg_accre_enhan_fact_in, micro_mg_autocon_fact_in, & !++ trude
     micro_mg_autocon_exp_in, micro_mg_homog_size_in, & !++ trude
     micro_mg_vtrmi_factor_in, micro_mg_effi_factor_in,  micro_mg_iaccr_factor_in,&
     micro_mg_max_nicons_in, &
!-- trude
     allow_sed_supersat_in, do_sb_physics_in, &
     nccons_in, nicons_in, ncnst_in, ninst_in, errstring)

  use micro_mg_utils, only: micro_mg_utils_init

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! initialize constants for MG microphysics
  !
  ! Author: Andrew Gettelman Dec 2005
  !
  !-----------------------------------------------------------------------

  integer,  intent(in)  :: kind         ! Kind used for reals
  real(r8), intent(in)  :: gravit
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: tmelt_in     ! Freezing point of water (K)
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice
  real(r8), intent(in)  :: rhmini_in    ! Minimum rh for ice cloud fraction > 0.
  real(r8), intent(in)  :: micro_mg_dcs

!++ag
!MG3 dense precipitating ice. Note, only 1 can be true, or both false.
  logical,  intent(in)  :: micro_mg_do_graupel_in    ! .true. = configure with graupel
                                                   ! .false. = no graupel (hail possible)
  logical,  intent(in)  :: micro_mg_do_hail_in    ! .true. = configure with hail
                                                   ! .false. = no hail (graupel possible)
!--ag

  logical,  intent(in)  :: microp_uniform_in    ! .true. = configure uniform for sub-columns
                                            ! .false. = use w/o sub-columns (standard)
  logical,  intent(in)  :: do_cldice_in     ! .true. = do all processes (standard)
                                            ! .false. = skip all processes affecting
                                            !           cloud ice
  logical,  intent(in)  :: use_hetfrz_classnuc_in ! use heterogeneous freezing

  character(len=16),intent(in)  :: micro_mg_precip_frac_method_in  ! type of precipitation fraction method
  real(r8),         intent(in)  :: micro_mg_berg_eff_factor_in     ! berg efficiency factor
!++ trude
  real(r8),         intent(in)  :: micro_mg_accre_enhan_fact_in     !accretion enhancment factor
  real(r8),         intent(in) ::  micro_mg_autocon_fact_in    !autconversion prefactor
  real(r8),         intent(in) ::  micro_mg_autocon_exp_in    !autconversion exponent factor
  real(r8),         intent(in) ::  micro_mg_homog_size_in  ! size of homoegenous freezing ice
  real(r8),         intent(in)  :: micro_mg_vtrmi_factor_in    !factor for ice fall velocity
  real(r8),         intent(in)  :: micro_mg_effi_factor_in    !factor for ice effective radius
  real(r8),         intent(in)  :: micro_mg_iaccr_factor_in  ! ice accretion factor
  real(r8),         intent(in)  :: micro_mg_max_nicons_in ! maximum number ice crystal allowed 
! -- trude

  logical,  intent(in)  ::  allow_sed_supersat_in ! allow supersaturated conditions after sedimentation loop
  logical,  intent(in)  ::  do_sb_physics_in ! do SB autoconversion and accretion physics

  logical,  intent(in)  :: nccons_in
  logical,  intent(in)  :: nicons_in
  real(r8), intent(in)  :: ncnst_in
  real(r8), intent(in)  :: ninst_in

  character(128), intent(out) :: errstring    ! Output status (non-blank for error return)

  !-----------------------------------------------------------------------
 
  dcs = micro_mg_dcs

  ! Initialize subordinate utilities module.
  call micro_mg_utils_init(kind, rair, rh2o, cpair, tmelt_in, latvap, latice, &
       dcs, errstring)

  if (trim(errstring) /= "") return

  ! declarations for MG code (transforms variable names)

  g= gravit                 ! gravity
  r= rair                   ! dry air gas constant: note units(phys_constants are in J/K/kmol)
  rv= rh2o                  ! water vapor gas constant
  cpp = cpair               ! specific heat of dry air
  tmelt = tmelt_in
  rhmini = rhmini_in
  micro_mg_precip_frac_method = micro_mg_precip_frac_method_in
  micro_mg_berg_eff_factor    = micro_mg_berg_eff_factor_in
  
  !++ trude
  micro_mg_accre_enhan_fact   =  micro_mg_accre_enhan_fact_in
  micro_mg_autocon_fact  = micro_mg_autocon_fact_in
  micro_mg_autocon_exp = micro_mg_autocon_exp_in
  micro_mg_homog_size = micro_mg_homog_size_in
  micro_mg_vtrmi_factor = micro_mg_vtrmi_factor_in
  micro_mg_effi_factor = micro_mg_effi_factor_in
  micro_mg_iaccr_factor=micro_mg_iaccr_factor_in
  micro_mg_max_nicons = micro_mg_max_nicons_in
   ! -- trude
  allow_sed_supersat          = allow_sed_supersat_in
  do_sb_physics               = do_sb_physics_in

 
  nccons = nccons_in
  nicons = nicons_in
  ncnst = ncnst_in
  ninst = ninst_in

  ! latent heats

  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

  ! flags
  microp_uniform = microp_uniform_in
  do_cldice  = do_cldice_in
  use_hetfrz_classnuc = use_hetfrz_classnuc_in
!++ag
  do_hail          = micro_mg_do_hail_in
  do_graupel       = micro_mg_do_graupel_in
!--ag

  ! typical air density at 850 mb

  rhosu = 85000._r8/(rair * tmelt)

  ! Maximum temperature at which snow is allowed to exist
  snowmelt = tmelt + 2._r8
  ! Minimum temperature at which rain is allowed to exist
  rainfrze = tmelt - 40._r8

  ! Ice nucleation temperature
  icenuct  = tmelt - 5._r8

  ! Define constants to help speed up code (this limits calls to gamma function)
  gamma_br_plus1=gamma(1._r8+br)
  gamma_br_plus4=gamma(4._r8+br)
  gamma_bs_plus1=gamma(1._r8+bs)
  gamma_bs_plus4=gamma(4._r8+bs)
  gamma_bi_plus1=gamma(1._r8+bi)
  gamma_bi_plus4=gamma(4._r8+bi)
  gamma_bj_plus1=gamma(1._r8+bj)
  gamma_bj_plus4=gamma(4._r8+bj)

  xxlv_squared=xxlv**2
  xxls_squared=xxls**2

end subroutine micro_mg_init

!===============================================================================
!microphysics routine for each timestep goes here...

! First adjust input mixing ratios to instantaneously adjust to initial temperatures
subroutine micro_mg_tend_pre_adjust( &
     mgncol,             nlev,             deltatin,           &
     t,                                  &
     qcn,                          qin,                          &
     ncn,                          nin,                          &
     qrn,                                &
     nrn,                                &
     qgr,                          ngr,                          &
     tlat,                               &
     qctend,                       qitend,                       &
     nctend,                       nitend,                       &
     qrtend,                       nrtend,                       &
     qgtend,                       ngtend                        &
     )

  ! input arguments
  integer,  intent(in) :: mgncol         ! number of microphysics columns
  integer,  intent(in) :: nlev           ! number of layers
  real(r8), intent(in) :: deltatin       ! time step (s)
  real(r8), intent(inout) :: t(mgncol,nlev) ! input temperature (K)

  ! note: all input cloud variables are grid-averaged
  real(r8), intent(inout) :: qcn(mgncol,nlev)       ! cloud water mixing ratio (kg/kg)
  real(r8), intent(inout) :: qin(mgncol,nlev)       ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(inout) :: ncn(mgncol,nlev)       ! cloud water number conc (1/kg)
  real(r8), intent(inout) :: nin(mgncol,nlev)       ! cloud ice number conc (1/kg)

  real(r8), intent(inout) :: qrn(mgncol,nlev)       ! rain mixing ratio (kg/kg)
  real(r8), intent(inout) :: nrn(mgncol,nlev)       ! rain number conc (1/kg)
  real(r8), intent(inout) :: qgr(mgncol,nlev)       ! graupel/hail mixing ratio (kg/kg)
  real(r8), intent(inout) :: ngr(mgncol,nlev)       ! graupel/hail number conc (1/kg)

  real(r8), intent(out) :: tlat(mgncol,nlev)         ! latent heating rate       (W/kg)
  real(r8), intent(out) :: qctend(mgncol,nlev)       ! microphysical tendency qc (1/s)
  real(r8), intent(out) :: qitend(mgncol,nlev)       ! microphysical tendency qi (1/s)
  real(r8), intent(out) :: nctend(mgncol,nlev)       ! microphysical tendency nc (1/(kg*s))
  real(r8), intent(out) :: nitend(mgncol,nlev)       ! microphysical tendency ni (1/(kg*s))
  real(r8), intent(out) :: qrtend(mgncol,nlev)       ! microphysical tendency qr (1/s)
  real(r8), intent(out) :: nrtend(mgncol,nlev)       ! microphysical tendency nr (1/(kg*s))
  real(r8), intent(out) :: qgtend(mgncol,nlev)       ! microphysical tendency qg (1/s)
  real(r8), intent(out) :: ngtend(mgncol,nlev)       ! microphysical tendency ng (1/(kg*s))

  ! Locally defined variables

  ! loop array variables
  ! "i" and "k" are column/level iterators for internal (MG) variables
  integer i, k
  real(r8) :: deltat            ! sub-time step (s)
  ! dummy variables
  real(r8) :: dum
  real(r8) :: dum1

  ! local copies of input variables
  real(r8) :: qc(mgncol,nlev)      ! cloud liquid mixing ratio (kg/kg)
  real(r8) :: qi(mgncol,nlev)      ! cloud ice mixing ratio (kg/kg)
  real(r8) :: nc(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: ni(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: qr(mgncol,nlev)      ! rain mixing ratio (kg/kg)
  real(r8) :: nr(mgncol,nlev)      ! rain number concentration (1/kg)
  real(r8) :: qg(mgncol,nlev)      ! graupel mixing ratio (kg/kg)
  real(r8) :: ng(mgncol,nlev)      ! graupel number concentration (1/kg)

  ! Rates/tendencies due to:
  ! Instantaneous snow melting
  real(r8) :: minstsm(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstsm(mgncol,nlev)    ! number concentration
  ! Instantaneous graupel melting
  real(r8) :: minstgm(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstgm(mgncol,nlev)    ! number concentration
  ! Instantaneous rain freezing
  real(r8) :: minstrf(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstrf(mgncol,nlev)    ! number concentration

  ! Copies of input concentrations that may be changed internally.
  deltat = deltatin
  qc = qcn
  nc = ncn
  qi = qin
  ni = nin
  qr = qrn
  nr = nrn
  qg = qgr
  ng = ngr
  qctend = 0._r8
  nctend = 0._r8
  qitend = 0._r8
  nitend = 0._r8
  qrtend = 0._r8
  nrtend = 0._r8
  qgtend = 0._r8
  ngtend = 0._r8
  minstsm = 0._r8
  ninstsm = 0._r8
  minstgm = 0._r8
  ninstgm = 0._r8
  minstrf = 0._r8
  ninstrf = 0._r8
  tlat=0._r8

  !=============================================================================
  do k=1,nlev

     do i=1,mgncol

        ! calculate instantaneous precip processes (melting and homogeneous freezing)
        !++kt
        ! melting of ice 

        if (t(i,k) > snowmelt) then
           if (qi(i,k) > 0._r8) then

              ! make sure melting ice/snow doesn't reduce temperature below threshold
              dum = -xlf/cpp*qi(i,k)
              if (t(i,k)+dum < snowmelt) then
                 dum = (t(i,k)-snowmelt)*cpp/xlf
                 dum = dum/qi(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              minstsm(i,k) = dum*qi(i,k)
              ninstsm(i,k) = dum*ni(i,k)

              dum1=-xlf*minstsm(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1

              qi(i,k) = max(qi(i,k) - minstsm(i,k), 0._r8)
              ni(i,k) = max(ni(i,k) - ninstsm(i,k), 0._r8)
       !--kt
              qr(i,k) = max(qr(i,k) + minstsm(i,k), 0._r8)
              nr(i,k) = max(nr(i,k) + ninstsm(i,k), 0._r8)
           end if
        end if

     end do
  end do 

!++ag

        ! melting of graupel at +2 C

  do k=1,nlev
     do i=1,mgncol

        if (t(i,k) > snowmelt) then
           if (qg(i,k) > 0._r8) then

              ! make sure melting graupel doesn't reduce temperature below threshold
              dum = -xlf/cpp*qg(i,k)
              if (t(i,k)+dum < snowmelt) then
                 dum = (t(i,k)-snowmelt)*cpp/xlf
                 dum = dum/qg(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              minstgm(i,k) = dum*qg(i,k)
              ninstgm(i,k) = dum*ng(i,k)

              dum1=-xlf*minstgm(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1

              qg(i,k) = max(qg(i,k) - minstgm(i,k), 0._r8)
              ng(i,k) = max(ng(i,k) - ninstgm(i,k), 0._r8)
              qr(i,k) = max(qr(i,k) + minstgm(i,k), 0._r8)
              nr(i,k) = max(nr(i,k) + ninstgm(i,k), 0._r8)
           end if
        end if

     end do
  end do 
!--ag

  do k=1,nlev
    do i=1,mgncol
        ! freezing of rain at -5 C

        if (t(i,k) < rainfrze) then

           if (qr(i,k) > 0._r8) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*qr(i,k)
              if (t(i,k)+dum > rainfrze) then
                 dum = -(t(i,k)-rainfrze)*cpp/xlf
                 dum = dum/qr(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              minstrf(i,k) = dum*qr(i,k)
              ninstrf(i,k) = dum*nr(i,k)

              ! heating tendency
              dum1 = xlf*minstrf(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1

              qr(i,k) = max(qr(i,k) - minstrf(i,k), 0._r8)
              nr(i,k) = max(nr(i,k) - ninstrf(i,k), 0._r8)
!++ag
! freeze rain to graupel not snow.
              if(do_hail.or.do_graupel) then
                 qg(i,k) = max(qg(i,k) + minstrf(i,k), 0._r8)
                 ng(i,k) = max(ng(i,k) + ninstrf(i,k), 0._r8)
              else
!++kt for mg4 freeze rain to ice, not snow, if no graupel
                 qi(i,k) = max(qi(i,k) + minstrf(i,k), 0._r8)
                 ni(i,k) = max(ni(i,k) + ninstrf(i,k), 0._r8)
!--kt
              end if
!--ag
           end if
        end if
     end do
  end do 

  ! Calculate tendencies from these instantaneous changes
  qctend = (qc-qcn)/deltat
  nctend = (nc-ncn)/deltat
  qitend = (qi-qin)/deltat
  nitend = (ni-nin)/deltat
  qrtend = (qr-qrn)/deltat
  nrtend = (nr-nrn)/deltat
  qgtend = (qg-qgr)/deltat
  ngtend = (ng-ngr)/deltat

  ! Output updated variables for the next microphysics step
  ! ktc-note: these should probably be seperate output variables

  ! Update temperature profile
  ! This could be dangerous because of a lack of checks to make
  ! sure that t does not go negative (for now)
  do k=1,nlev
     do i=1,mgncol
        t(i,k) = t(i,k)+tlat(i,k)/cpp*deltat
     end do
  end do

  ! Update hydrometeors
  qcn = qc
  ncn = nc
  qin = qi
  nin = ni
  qrn = qr
  nrn = nr
  qgr = qg
  ngr = ng

  if(minval(qc).lt.0._r8) & 
     write(iulog,*) "KTC inst qc < 0 : min=",minval(qc)
  if(minval(qi).lt.0._r8) & 
     write(iulog,*) "KTC inst qi < 0 : min=",minval(qi)
  if(minval(qr).lt.0._r8) & 
     write(iulog,*) "KTC inst qr < 0 : min=",minval(qr)
  if(minval(qg).lt.0._r8) & 
     write(iulog,*) "KTC inst qg < 0 : min=",minval(qg)
  if(minval(nc).lt.0._r8) & 
     write(iulog,*) "KTC inst nc < 0 : min=",minval(nc)
  if(minval(ni).lt.0._r8) & 
     write(iulog,*) "KTC inst ni < 0 : min=",minval(ni)
  if(minval(nr).lt.0._r8) & 
     write(iulog,*) "KTC inst nr < 0 : min=",minval(nr)
  if(minval(ng).lt.0._r8) & 
     write(iulog,*) "KTC inst ng < 0 : min=",minval(ng)

end subroutine micro_mg_tend_pre_adjust

! Then call full mg_tend

subroutine micro_mg_tend ( &
     mgncol,             nlev,               deltatin,           &
     t,                            q,                            &
     qcn,                          qin,                          &
     ncn,                          nin,                          &
     qrn,                          qsn,                          &
     nrn,                          nsn,                          &
!++ag
     qgr,                          ngr,                          &
!--ag
     relvar,                       accre_enhan,                  &
     p,                            pdel,                         &
     cldn,    liqcldf,        icecldf,       qsatfac,            &
     qcsinksum_rate1ord,                                         &
     naai,                         npccn,                        &
     rndst,                        nacon,                        &
     tlat,                         qvlat,                        &
     qctend,                       qitend,   qitendresid,        &
     nctend,                       nitend,                       &
     qrtend,                       qstend,                       &
     nrtend,                       nstend,                       &
!++ag
     qgtend,                       ngtend,                       &
!--ag
     effc,               effc_fn,            effi,               &
     sadice,                       sadsnow,                      &
     prect,                        preci,                        &
     nevapr,                       evapsnow,                     &
     am_evp_st,                                                  &
     prain,                        prodsnow,                     &
     cmeout,                       deffi,                        &
     pgamrad,                      lamcrad,                      &
     qsout,                        dsout,                        &
!++ag
     qgout,     ngout,             dgout,                        &
!--ag
     lflx,               iflx,                                   &
!++ag
     gflx,                                                       &
!--ag
     rflx,               sflx,               qrout,              &
!++ag
     reff_rain,          reff_snow,          reff_grau,          &
!--ag
     qcsevap,            qisevap,            qvres,              &
     cmeitot,            vtrmc,              vtrmi,              &
     umr,                          ums,                          &
!++ag
     umg,                qgsedten,                               &
!--ag
     qcsedten,                     qisedten,                     &
     qrsedten,                     qssedten,                     &
     pratot,                       prctot,                       &
     mnuccctot,          mnuccttot,          msacwitot,          &
     psacwstot,          bergstot,           bergtot,            &
!++ktc
     npsacwstot,                                                 &
!--ktc
     melttot,                      homotot,                      &
     qcrestot,           prcitot,            praitot,            &
!++ag
     qirestot,           mnuccrtot,          mnuccritot, pracstot,           &
!--ag
     meltsdttot,         frzrdttot,          mnuccdtot,          &
!++ag
     pracgtot,           psacwgtot,          pgsacwtot,          &
     pgracstot,          prdgtot,           &
     qmultgtot,          qmultrgtot,         psacrtot,           &
     npracgtot,          nscngtot,           ngracstot,          &
     nmultgtot,          nmultrgtot,         npsacwgtot,         & 
!--ag
     nrout,                        nsout,                        &
     refl,               arefl,              areflz,             &
     frefl,              csrfl,              acsrfl,             &
     fcsrfl,                       rercld,                       &
     ncai,                         ncal,                         &
     qrout2,                       qsout2,                       &
     nrout2,                       nsout2,                       &
     drout2,                       dsout2,                       &
!++ag     
     qgout2,        ngout2,        dgout2,    freqg,                   &
!--ag
     freqs,                        freqr,                        &
     nfice,                        qcrat,                        &
     errstring, & ! Below arguments are "optional" (pass null pointers to omit).
     tnd_qsnow,          tnd_nsnow,          re_ice,             &
     prer_evap,                                                      &
     frzimm,             frzcnt,             frzdep,             &
     lamsout,            n0sout ,        lamrout,          n0rout)

  ! Constituent properties.
  use micro_mg_utils, only: &
       mg_liq_props, &
       mg_ice_props, &
       mg_rain_props, &
!++ag
       mg_graupel_props,&
!--ag
       mg_snow_props

  ! Size calculation functions.
  use micro_mg_utils, only: &
       size_dist_param_liq, &
       size_dist_param_basic, &
       avg_diameter

  ! Microphysical processes.
  use micro_mg_utils, only: &
       ice_deposition_sublimation, &
       sb2001v2_liq_autoconversion,&
       sb2001v2_accre_cld_water_rain,&       
       kk2000_liq_autoconversion, &
       ice_autoconversion, &
       immersion_freezing, &
       contact_freezing, &
       snow_self_aggregation, &
       accrete_cloud_water_snow, &
       secondary_ice_production, &
       accrete_rain_snow, &
       heterogeneous_rain_freezing, &
       accrete_cloud_water_rain, &
       self_collection_rain, &
       accrete_cloud_ice_snow, &
       bergeron_process_snow, &
!++ag
     graupel_collecting_snow, &
     graupel_collecting_rain, &
     graupel_collecting_cld_water, &
     graupel_riming_liquid_snow, &
     graupel_rain_riming_snow, &
     graupel_rime_splintering, &
!     graupel_sublimate_evap
!--ag
!++kt
     ice_deposition_sublimation_mg4, &
     evaporate_sublimate_precip_mg4, &
     evaporate_sublimate_precip_graupel_mg4, &
     ice_self_aggregation, &
     accrete_cloud_water_ice, &
     accrete_rain_ice, &
     access_lookup_table, &
     access_lookup_table_coll, &
     tsize, isize, jsize, rcollsize
!--kt

  !Authors: Hugh Morrison, Andrew Gettelman, NCAR, Peter Caldwell, LLNL
  ! e-mail: morrison@ucar.edu, andrew@ucar.edu

  ! input arguments
  integer,  intent(in) :: mgncol         ! number of microphysics columns
  integer,  intent(in) :: nlev           ! number of layers
  real(r8), intent(in) :: deltatin       ! time step (s)
  real(r8), intent(in) :: t(mgncol,nlev) ! input temperature (K)
  real(r8), intent(in) :: q(mgncol,nlev) ! input h20 vapor mixing ratio (kg/kg)

  ! note: all input cloud variables are grid-averaged
  real(r8), intent(in) :: qcn(mgncol,nlev)       ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(mgncol,nlev)       ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: ncn(mgncol,nlev)       ! cloud water number conc (1/kg)
  real(r8), intent(in) :: nin(mgncol,nlev)       ! cloud ice number conc (1/kg)

  real(r8), intent(in) :: qrn(mgncol,nlev)       ! rain mixing ratio (kg/kg)
  real(r8), intent(in) :: qsn(mgncol,nlev)       ! snow mixing ratio (kg/kg)
  real(r8), intent(in) :: nrn(mgncol,nlev)       ! rain number conc (1/kg)
  real(r8), intent(in) :: nsn(mgncol,nlev)       ! snow number conc (1/kg)

!++ag
  real(r8), intent(in) :: qgr(mgncol,nlev)       ! graupel/hail mixing ratio (kg/kg)
  real(r8), intent(in) :: ngr(mgncol,nlev)       ! graupel/hail number conc (1/kg)
!--ag

  real(r8), intent(in) :: relvar(mgncol,nlev)      ! cloud water relative variance (-)
  real(r8), intent(in) :: accre_enhan(mgncol,nlev) ! optional accretion
                                             ! enhancement factor (-)

  real(r8), intent(in) :: p(mgncol,nlev)        ! air pressure (pa)
  real(r8), intent(in) :: pdel(mgncol,nlev)     ! pressure difference across level (pa)

  real(r8), intent(in) :: cldn(mgncol,nlev)      ! cloud fraction (no units)
  real(r8), intent(in) :: liqcldf(mgncol,nlev)   ! liquid cloud fraction (no units)
  real(r8), intent(in) :: icecldf(mgncol,nlev)   ! ice cloud fraction (no units)
  real(r8), intent(in) :: qsatfac(mgncol,nlev)   ! subgrid cloud water saturation scaling factor (no units)

  ! used for scavenging
  ! Inputs for aerosol activation
  real(r8), intent(in) :: naai(mgncol,nlev)     ! ice nucleation number (from microp_aero_ts) (1/kg)
  real(r8), intent(in) :: npccn(mgncol,nlev)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  ! Note that for these variables, the dust bin is assumed to be the last index.
  ! (For example, in CAM, the last dimension is always size 4.)
  real(r8), intent(in) :: rndst(:,:,:)  ! radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
  real(r8), intent(in) :: nacon(:,:,:) ! number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)
  
  ! output arguments

  real(r8), intent(out) :: qcsinksum_rate1ord(mgncol,nlev) ! 1st order rate for
  ! direct cw to precip conversion
  real(r8), intent(out) :: tlat(mgncol,nlev)         ! latent heating rate       (W/kg)
  real(r8), intent(out) :: qvlat(mgncol,nlev)        ! microphysical tendency qv (1/s)
  real(r8), intent(out) :: qctend(mgncol,nlev)       ! microphysical tendency qc (1/s)
  real(r8), intent(out) :: qitend(mgncol,nlev)       ! microphysical tendency qi (1/s)
  real(r8), intent(out) :: qitendresid(mgncol,nlev)
  real(r8), intent(out) :: nctend(mgncol,nlev)       ! microphysical tendency nc (1/(kg*s))
  real(r8), intent(out) :: nitend(mgncol,nlev)       ! microphysical tendency ni (1/(kg*s))

  real(r8), intent(out) :: qrtend(mgncol,nlev)       ! microphysical tendency qr (1/s)
  real(r8), intent(out) :: qstend(mgncol,nlev)       ! microphysical tendency qs (1/s)
  real(r8), intent(out) :: nrtend(mgncol,nlev)       ! microphysical tendency nr (1/(kg*s))
  real(r8), intent(out) :: nstend(mgncol,nlev)       ! microphysical tendency ns (1/(kg*s))
!++ag
  real(r8), intent(out) :: qgtend(mgncol,nlev)       ! microphysical tendency qg (1/s)
  real(r8), intent(out) :: ngtend(mgncol,nlev)       ! microphysical tendency ng (1/(kg*s))
!--ag
  real(r8), intent(out) :: effc(mgncol,nlev)         ! droplet effective radius (micron)
  real(r8), intent(out) :: effc_fn(mgncol,nlev)      ! droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8), intent(out) :: effi(mgncol,nlev)         ! cloud ice effective radius (micron)
  real(r8), intent(out) :: sadice(mgncol,nlev)       ! cloud ice surface area density (cm2/cm3)
  real(r8), intent(out) :: sadsnow(mgncol,nlev)      ! cloud snow surface area density (cm2/cm3)
  real(r8), intent(out) :: prect(mgncol)             ! surface precip rate (m/s)
  real(r8), intent(out) :: preci(mgncol)             ! cloud ice/snow precip rate (m/s)
  real(r8), intent(out) :: nevapr(mgncol,nlev)       ! evaporation rate of rain + snow (1/s)
  real(r8), intent(out) :: evapsnow(mgncol,nlev)     ! sublimation rate of snow (1/s)
  real(r8), intent(out) :: am_evp_st(mgncol,nlev)    ! stratiform evaporation area (frac)
  real(r8), intent(out) :: prain(mgncol,nlev)        ! production of rain + snow (1/s)
  real(r8), intent(out) :: prodsnow(mgncol,nlev)     ! production of snow (1/s)
  real(r8), intent(out) :: cmeout(mgncol,nlev)       ! evap/sub of cloud (1/s)
  real(r8), intent(out) :: deffi(mgncol,nlev)        ! ice effective diameter for optics (radiation) (micron)
  real(r8), intent(out) :: pgamrad(mgncol,nlev)      ! ice gamma parameter for optics (radiation) (no units)
  real(r8), intent(out) :: lamcrad(mgncol,nlev)      ! slope of droplet distribution for optics (radiation) (1/m)
  real(r8), intent(out) :: qsout(mgncol,nlev)        ! snow mixing ratio (kg/kg)
  real(r8), intent(out) :: dsout(mgncol,nlev)        ! snow diameter (m)
  real(r8), intent(out) :: lflx(mgncol,nlev+1)       ! grid-box average liquid condensate flux (kg m^-2 s^-1)
  real(r8), intent(out) :: iflx(mgncol,nlev+1)       ! grid-box average ice condensate flux (kg m^-2 s^-1)
  real(r8), intent(out) :: rflx(mgncol,nlev+1)       ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8), intent(out) :: sflx(mgncol,nlev+1)       ! grid-box average snow flux (kg m^-2 s^-1)
!++ag
  real(r8), intent(out) :: gflx(mgncol,nlev+1)       ! grid-box average graupel/hail flux (kg m^-2 s^-1)
!--ag
  real(r8), intent(out) :: qrout(mgncol,nlev)        ! grid-box average rain mixing ratio (kg/kg)
  real(r8), intent(out) :: reff_rain(mgncol,nlev)    ! rain effective radius (micron)
  real(r8), intent(out) :: reff_snow(mgncol,nlev)    ! snow effective radius (micron)
!++ag
  real(r8), intent(out) :: reff_grau(mgncol,nlev)    ! graupel effective radius (micron)
!--ag
  real(r8), intent(out) :: qcsevap(mgncol,nlev)      ! cloud water evaporation due to sedimentation (1/s)
  real(r8), intent(out) :: qisevap(mgncol,nlev)      ! cloud ice sublimation due to sublimation (1/s)
  real(r8), intent(out) :: qvres(mgncol,nlev)        ! residual condensation term to ensure RH < 100% (1/s)
  real(r8), intent(out) :: cmeitot(mgncol,nlev)      ! grid-mean cloud ice sub/dep (1/s)
  real(r8), intent(out) :: vtrmc(mgncol,nlev)        ! mass-weighted cloud water fallspeed (m/s)
  real(r8), intent(out) :: vtrmi(mgncol,nlev)        ! mass-weighted cloud ice fallspeed (m/s)
  real(r8), intent(out) :: umr(mgncol,nlev)          ! mass weighted rain fallspeed (m/s)
  real(r8), intent(out) :: ums(mgncol,nlev)          ! mass weighted snow fallspeed (m/s)
!++ag
  real(r8), intent(out) :: umg(mgncol,nlev)          ! mass weighted graupel/hail fallspeed (m/s)
  real(r8), intent(out) :: qgsedten(mgncol,nlev)     ! qg sedimentation tendency (1/s)
!--ag
  real(r8), intent(out) :: qcsedten(mgncol,nlev)     ! qc sedimentation tendency (1/s)
  real(r8), intent(out) :: qisedten(mgncol,nlev)     ! qi sedimentation tendency (1/s)
  real(r8), intent(out) :: qrsedten(mgncol,nlev)     ! qr sedimentation tendency (1/s)
  real(r8), intent(out) :: qssedten(mgncol,nlev)     ! qs sedimentation tendency (1/s)

  ! microphysical process rates for output (mixing ratio tendencies) (all have units of 1/s)
  real(r8), intent(out) :: pratot(mgncol,nlev)          ! accretion of cloud by rain
  real(r8), intent(out) :: prctot(mgncol,nlev)          ! autoconversion of cloud to rain
  real(r8), intent(out) :: mnuccctot(mgncol,nlev)       ! mixing ratio tend due to immersion freezing
  real(r8), intent(out) :: mnuccttot(mgncol,nlev)       ! mixing ratio tend due to contact freezing
  real(r8), intent(out) :: msacwitot(mgncol,nlev)       ! mixing ratio tend due to H-M splintering
  real(r8), intent(out) :: psacwstot(mgncol,nlev)       ! collection of cloud water by snow
!++ktc
  real(r8), intent(out) :: npsacwstot(mgncol,nlev)      ! number collection of coud water by snow
!--ktc
  real(r8), intent(out) :: bergstot(mgncol,nlev)        ! bergeron process on snow
  real(r8), intent(out) :: bergtot(mgncol,nlev)         ! bergeron process on cloud ice
  real(r8), intent(out) :: melttot(mgncol,nlev)         ! melting of cloud ice
  real(r8), intent(out) :: homotot(mgncol,nlev)         ! homogeneous freezing cloud water
  real(r8), intent(out) :: qcrestot(mgncol,nlev)        ! residual cloud condensation due to removal of excess supersat
  real(r8), intent(out) :: prcitot(mgncol,nlev)         ! autoconversion of cloud ice to snow
  real(r8), intent(out) :: praitot(mgncol,nlev)         ! accretion of cloud ice by snow
  real(r8), intent(out) :: qirestot(mgncol,nlev)        ! residual ice deposition due to removal of excess supersat
  real(r8), intent(out) :: mnuccrtot(mgncol,nlev)       ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
  real(r8), intent(out) :: mnuccritot(mgncol,nlev)      ! mixing ratio tendency due to heterogeneous freezing of rain to ice (1/s)
  real(r8), intent(out) :: pracstot(mgncol,nlev)        ! mixing ratio tendency due to accretion of rain by snow (1/s)
  real(r8), intent(out) :: meltsdttot(mgncol,nlev)      ! latent heating rate due to melting of snow  (W/kg)
  real(r8), intent(out) :: frzrdttot(mgncol,nlev)       ! latent heating rate due to homogeneous freezing of rain (W/kg)
  real(r8), intent(out) :: mnuccdtot(mgncol,nlev)       ! mass tendency from ice nucleation
!++ag Hail/Graupel Tendencies
  real(r8), intent(out) :: pracgtot(mgncol,nlev)        ! change in q collection rain by graupel  (precipf)
  real(r8), intent(out) :: psacwgtot(mgncol,nlev)       ! change in q collection droplets by graupel (lcldm)
  real(r8), intent(out) :: pgsacwtot(mgncol,nlev)       ! conversion q to graupel due to collection droplets by snow  (lcldm)
  real(r8), intent(out) :: pgracstot(mgncol,nlev)       ! conversion q to graupel due to collection rain by snow (precipf)
  real(r8), intent(out) :: prdgtot(mgncol,nlev)         ! dep of graupel (precipf)
!  real(r8), intent(out) :: eprdgtot(mgncol,nlev)        ! sub of graupel (precipf)
  real(r8), intent(out) :: qmultgtot(mgncol,nlev)       ! change q due to ice mult droplets/graupel  (lcldm)
  real(r8), intent(out) :: qmultrgtot(mgncol,nlev)      ! change q due to ice mult rain/graupel (precipf)
  real(r8), intent(out) :: psacrtot(mgncol,nlev)        ! conversion due to coll of snow by rain (precipf)
  real(r8), intent(out) :: npracgtot(mgncol,nlev)       ! change n collection rain by graupel  (precipf)
  real(r8), intent(out) :: nscngtot(mgncol,nlev)        ! change n conversion to graupel due to collection droplets by snow (lcldm)
  real(r8), intent(out) :: ngracstot(mgncol,nlev)       ! change n conversion to graupel due to collection rain by snow (precipf)
  real(r8), intent(out) :: nmultgtot(mgncol,nlev)       ! ice mult due to acc droplets by graupel  (lcldm)
  real(r8), intent(out) :: nmultrgtot(mgncol,nlev)      ! ice mult due to acc rain by graupel  (precipf)
  real(r8), intent(out) :: npsacwgtot(mgncol,nlev)      ! change n collection droplets by graupel (lcldm?)
!--ag
  real(r8), intent(out) :: nrout(mgncol,nlev)        ! rain number concentration (1/m3)
  real(r8), intent(out) :: nsout(mgncol,nlev)        ! snow number concentration (1/m3)
  real(r8), intent(out) :: refl(mgncol,nlev)         ! analytic radar reflectivity
  real(r8), intent(out) :: arefl(mgncol,nlev)        ! average reflectivity will zero points outside valid range
  real(r8), intent(out) :: areflz(mgncol,nlev)       ! average reflectivity in z.
  real(r8), intent(out) :: frefl(mgncol,nlev)        ! fractional occurrence of radar reflectivity
  real(r8), intent(out) :: csrfl(mgncol,nlev)        ! cloudsat reflectivity
  real(r8), intent(out) :: acsrfl(mgncol,nlev)       ! cloudsat average
  real(r8), intent(out) :: fcsrfl(mgncol,nlev)       ! cloudsat fractional occurrence of radar reflectivity
  real(r8), intent(out) :: rercld(mgncol,nlev)       ! effective radius calculation for rain + cloud
  real(r8), intent(out) :: ncai(mgncol,nlev)         ! output number conc of ice nuclei available (1/m3)
  real(r8), intent(out) :: ncal(mgncol,nlev)         ! output number conc of CCN (1/m3)
  real(r8), intent(out) :: qrout2(mgncol,nlev)       ! copy of qrout as used to compute drout2
  real(r8), intent(out) :: qsout2(mgncol,nlev)       ! copy of qsout as used to compute dsout2
  real(r8), intent(out) :: nrout2(mgncol,nlev)       ! copy of nrout as used to compute drout2
  real(r8), intent(out) :: nsout2(mgncol,nlev)       ! copy of nsout as used to compute dsout2
  real(r8), intent(out) :: drout2(mgncol,nlev)       ! mean rain particle diameter (m)
  real(r8), intent(out) :: dsout2(mgncol,nlev)       ! mean snow particle diameter (m)
  real(r8), intent(out) :: freqs(mgncol,nlev)        ! fractional occurrence of snow
  real(r8), intent(out) :: freqr(mgncol,nlev)        ! fractional occurrence of rain
  real(r8), intent(out) :: nfice(mgncol,nlev)        ! fractional occurrence of ice
  real(r8), intent(out) :: qcrat(mgncol,nlev)        ! limiter for qc process rates (1=no limit --> 0. no qc)
!++ag
  real(r8), intent(out) :: qgout(mgncol,nlev)        ! graupel/hail mixing ratio (kg/kg)
  real(r8), intent(out) :: dgout(mgncol,nlev)        ! graupel/hail diameter (m)
  real(r8), intent(out) :: ngout(mgncol,nlev)        ! graupel/hail number concentration (1/m3)
!Not sure if these are needed since graupel/hail is prognostic?
  real(r8), intent(out) :: qgout2(mgncol,nlev)       ! copy of qgout as used to compute dgout2
  real(r8), intent(out) :: ngout2(mgncol,nlev)       ! copy of ngout as used to compute dgout2
  real(r8), intent(out) :: dgout2(mgncol,nlev)       ! mean graupel/hail particle diameter (m)
  real(r8), intent(out) :: freqg(mgncol,nlev)        ! fractional occurrence of graupel  

  real(r8), intent(out) :: lamsout(mgncol,nlev)
  real(r8), intent(out) :: n0sout(mgncol,nlev)
  real(r8), intent(out) :: lamrout(mgncol,nlev)
  real(r8), intent(out) :: n0rout(mgncol,nlev)

!--ag

  real(r8), intent(out) :: prer_evap(mgncol,nlev)

  character(128),   intent(out) :: errstring  ! output status (non-blank for error return)

  ! Tendencies calculated by external schemes that can replace MG's native
  ! process tendencies.

  ! Used with CARMA cirrus microphysics
  ! (or similar external microphysics model)
  real(r8), intent(in) :: tnd_qsnow(:,:) ! snow mass tendency (kg/kg/s)
  real(r8), intent(in) :: tnd_nsnow(:,:) ! snow number tendency (#/kg/s)
  real(r8), intent(in) :: re_ice(:,:)    ! ice effective radius (m)

  ! From external ice nucleation.
  real(r8), intent(in) :: frzimm(:,:) ! Number tendency due to immersion freezing (1/cm3)
  real(r8), intent(in) :: frzcnt(:,:) ! Number tendency due to contact freezing (1/cm3)
  real(r8), intent(in) :: frzdep(:,:) ! Number tendency due to deposition nucleation (1/cm3)

  ! local workspace
  ! all units mks unless otherwise stated

  ! local copies of input variables
  real(r8) :: qc(mgncol,nlev)      ! cloud liquid mixing ratio (kg/kg)
  real(r8) :: qi(mgncol,nlev)      ! cloud ice mixing ratio (kg/kg)
  real(r8) :: nc(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: ni(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: qr(mgncol,nlev)      ! rain mixing ratio (kg/kg)
  real(r8) :: qs(mgncol,nlev)      ! snow mixing ratio (kg/kg)
  real(r8) :: nr(mgncol,nlev)      ! rain number concentration (1/kg)
  real(r8) :: ns(mgncol,nlev)      ! snow number concentration (1/kg)
!++ag
  real(r8) :: qg(mgncol,nlev)      ! graupel mixing ratio (kg/kg)
  real(r8) :: ng(mgncol,nlev)      ! graupel number concentration (1/kg)
  real(r8) :: rhogtmp              ! hail or graupel density (kg m-3)

!--ag

  ! general purpose variables
  real(r8) :: deltat            ! sub-time step (s)
  real(r8) :: mtime             ! the assumed ice nucleation timescale

  ! physical properties of the air at a given point
  real(r8) :: rho(mgncol,nlev)    ! density (kg m-3)
  real(r8) :: dv(mgncol,nlev)     ! diffusivity of water vapor
  real(r8) :: mu(mgncol,nlev)     ! viscosity
  real(r8) :: sc(mgncol,nlev)     ! schmidt number
  real(r8) :: rhof(mgncol,nlev)   ! density correction factor for fallspeed

  ! cloud fractions
  real(r8) :: precip_frac(mgncol,nlev) ! precip fraction assuming maximum overlap
  real(r8) :: cldm(mgncol,nlev)   ! cloud fraction
  real(r8) :: icldm(mgncol,nlev)  ! ice cloud fraction
  real(r8) :: lcldm(mgncol,nlev)  ! liq cloud fraction
  real(r8) :: qsfm(mgncol,nlev)   ! subgrid cloud water saturation scaling factor

  ! mass mixing ratios
  real(r8) :: qcic(mgncol,nlev)   ! in-cloud cloud liquid
  real(r8) :: qiic(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: qsic(mgncol,nlev)   ! in-precip snow
  real(r8) :: qric(mgncol,nlev)   ! in-precip rain
!++ag
  real(r8) :: qgic(mgncol,nlev)   ! in-precip graupel/hail
!++ag

  ! number concentrations
  real(r8) :: ncic(mgncol,nlev)   ! in-cloud droplet
  real(r8) :: niic(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: nsic(mgncol,nlev)   ! in-precip snow
  real(r8) :: nric(mgncol,nlev)   ! in-precip rain
!++ag
  real(r8) :: ngic(mgncol,nlev)   ! in-precip graupel/hail
!++ag

  ! maximum allowed ni value
  real(r8) :: nimax(mgncol,nlev)

  ! Size distribution parameters for:
  ! cloud ice
  real(r8) :: lami(mgncol,nlev)   ! slope
  real(r8) :: n0i(mgncol,nlev)    ! intercept
  ! cloud liquid
  real(r8) :: lamc(mgncol,nlev)   ! slope
  real(r8) :: pgam(mgncol,nlev)   ! spectral width parameter
  ! snow
  real(r8) :: lams(mgncol,nlev)   ! slope
  real(r8) :: n0s(mgncol,nlev)    ! intercept
  ! rain
  real(r8) :: lamr(mgncol,nlev)   ! slope
  real(r8) :: n0r(mgncol,nlev)    ! intercept
!++ag 
  ! graupel/hail
  real(r8) :: lamg(mgncol,nlev)   ! slope
  real(r8) :: n0g(mgncol,nlev)    ! intercept
  real(r8) :: bgtmp               ! tmp fall speed parameter
!--ag

  ! Rates/tendencies due to:

  ! deposition of cloud ice
  real(r8) :: vap_dep(mgncol,nlev)    ! deposition from vapor to ice PMC 12/3/12
  ! sublimation of cloud ice
  real(r8) :: ice_sublim(mgncol,nlev) ! sublimation from ice to vapor PMC 12/3/12
  ! ice nucleation
  real(r8) :: nnuccd(mgncol,nlev) ! number rate from deposition/cond.-freezing
  real(r8) :: mnuccd(mgncol,nlev) ! mass mixing ratio
  ! freezing of cloud water
  real(r8) :: mnuccc(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnuccc(mgncol,nlev) ! number concentration
  ! contact freezing of cloud water
  real(r8) :: mnucct(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnucct(mgncol,nlev) ! number concentration
  ! deposition nucleation in mixed-phase clouds (from external scheme)
  real(r8) :: mnudep(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnudep(mgncol,nlev) ! number concentration
  ! ice multiplication
  real(r8) :: msacwi(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nsacwi(mgncol,nlev) ! number concentration
  ! autoconversion of cloud droplets
  real(r8) :: prc(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: nprc(mgncol,nlev)   ! number concentration (rain)
  real(r8) :: nprc1(mgncol,nlev)  ! number concentration (cloud droplets)
  ! self-aggregation of snow & ice in mg4 (kt)
  real(r8) :: niagg(mgncol,nlev)  ! number concentration
  ! self-collection of rain
  real(r8) :: nragg(mgncol,nlev)  ! number concentration
  ! collection of droplets by snow
  real(r8) :: psacws(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: npsacws(mgncol,nlev)    ! number concentration
  ! collection of rain by snow
  real(r8) :: pracs(mgncol,nlev)  ! mass mixing ratio
  real(r8) :: npracs(mgncol,nlev) ! number concentration
  ! freezing of rain
  real(r8) :: mnuccr(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnuccr(mgncol,nlev) ! number concentration
  ! freezing of rain to form ice (mg add 4/26/13)
  real(r8) :: mnuccri(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: nnuccri(mgncol,nlev)    ! number concentration
  ! accretion of droplets by rain
  real(r8) :: pra(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: npra(mgncol,nlev)   ! number concentration
  ! evaporation of rain
  real(r8) :: pre(mgncol,nlev)    ! mass mixing ratio
  ! number evaporation
  real(r8) :: nsubi(mgncol,nlev)  ! cloud ice
  real(r8) :: nsubc(mgncol,nlev)  ! droplet
  real(r8) :: nsubr(mgncol,nlev)  ! rain
  ! bergeron process
  real(r8) :: berg(mgncol,nlev)   ! mass mixing ratio (cloud ice)

!++ag
  !graupel/hail processes
  real(r8) :: npracg(mgncol,nlev)  ! change n collection rain by graupel  (precipf)
  real(r8) :: nscng(mgncol,nlev)   ! change n conversion to graupel due to collection droplets by snow (lcldm)
  real(r8) :: ngracs(mgncol,nlev)  ! change n conversion to graupel due to collection rain by snow (precipf)
  real(r8) :: nmultg(mgncol,nlev)  ! ice mult due to acc droplets by graupel  (lcldm)
  real(r8) :: nmultrg(mgncol,nlev) ! ice mult due to acc rain by graupel  (precipf)
  real(r8) :: npsacwg(mgncol,nlev) ! change n collection droplets by graupel (lcldm)

  real(r8) :: psacr(mgncol,nlev)   ! conversion due to coll of snow by rain (precipf)
  real(r8) :: pracg(mgncol,nlev)   ! change in q collection rain by graupel  (precipf)
  real(r8) :: psacwg(mgncol,nlev)  ! change in q collection droplets by graupel (lcldm)
  real(r8) :: pgsacw(mgncol,nlev)  ! conversion q to graupel due to collection droplets by snow  (lcldm)
  real(r8) :: pgracs(mgncol,nlev)  ! conversion q to graupel due to collection rain by snow (precipf)
  real(r8) :: prdg(mgncol,nlev)    ! dep of graupel (precipf)
!  real(r8) :: eprdg(mgncol,nlev)   ! evap/sub of graupel (precipf)
  real(r8) :: qmultg(mgncol,nlev)  ! change q due to ice mult droplets/graupel  (lcldm)
  real(r8) :: qmultrg(mgncol,nlev) ! change q due to ice mult rain/graupel (precipf)
!--ag

  ! fallspeeds
  ! number-weighted
  real(r8) :: unr(mgncol,nlev)    ! rain
!++ag
  real(r8) :: ung(mgncol,nlev)    ! graupel/hail
!--ag
  ! air density corrected fallspeed parameters
  real(r8) :: arn(mgncol,nlev)    ! rain
  real(r8) :: asn(mgncol,nlev)    ! snow
!++a
  real(r8) :: agn(mgncol,nlev)    ! graupel
!--ag
  real(r8) :: acn(mgncol,nlev)    ! cloud droplet
  real(r8) :: ain(mgncol,nlev)    ! cloud ice
  real(r8) :: ajn(mgncol,nlev)    ! cloud small ice

  ! Mass of liquid droplets used with external heterogeneous freezing.
  real(r8) :: mi0l(mgncol)

  ! saturation vapor pressures
  real(r8) :: esl(mgncol,nlev)    ! liquid
  real(r8) :: esi(mgncol,nlev)    ! ice
  real(r8) :: esn                 ! checking for RH after rain evap

  ! saturation vapor mixing ratios
  real(r8) :: qvl(mgncol,nlev)    ! liquid
  real(r8) :: qvi(mgncol,nlev)    ! ice
  real(r8) :: qvn                 ! checking for RH after rain evap

  ! relative humidity
  real(r8) :: relhum(mgncol,nlev)

  ! parameters for cloud water and cloud ice sedimentation calculations
  real(r8) :: fc(mgncol,nlev)
  real(r8) :: fnc(mgncol,nlev)
  real(r8) :: fi(mgncol,nlev)
  real(r8) :: fni(mgncol,nlev)

!++ag
  real(r8) :: fg(mgncol,nlev)
  real(r8) :: fng(mgncol,nlev)
!--ag

  real(r8) :: fr(mgncol,nlev)
  real(r8) :: fnr(mgncol,nlev)
  real(r8) :: fs(mgncol,nlev)
  real(r8) :: fns(mgncol,nlev)

  real(r8) :: faloutc(nlev)
  real(r8) :: faloutnc(nlev)
  real(r8) :: falouti(nlev)
  real(r8) :: faloutni(nlev)

  real(r8) :: faloutr(nlev)
  real(r8) :: faloutnr(nlev)
  real(r8) :: falouts(nlev)
  real(r8) :: faloutns(nlev)

  real(r8) :: faltndc
  real(r8) :: faltndnc
  real(r8) :: faltndi
  real(r8) :: faltndni
  real(r8) :: faltndqie
  real(r8) :: faltndqce

  real(r8) :: faltndr
  real(r8) :: faltndnr
  real(r8) :: faltnds
  real(r8) :: faltndns

!++ag
  real(r8) :: faloutg(nlev)
  real(r8) :: faloutng(nlev)
  real(r8) :: faltndg
  real(r8) :: faltndng
!--ag

  real(r8) :: rainrt(mgncol,nlev)     ! rain rate for reflectivity calculation

  ! dummy variables
  real(r8) :: dum
  real(r8) :: dum1
  real(r8) :: dum2
!++ag
  real(r8) :: dum3
!--ag
  real(r8) :: dumni0
  real(r8) :: dumns0
  ! dummies for checking RH
  real(r8) :: qtmp
  real(r8) :: ttmp
  ! dummies for conservation check
  real(r8) :: ratio
  real(r8) :: tmpfrz
  ! dummies for in-cloud variables
  real(r8) :: dumc(mgncol,nlev)   ! qc
  real(r8) :: dumnc(mgncol,nlev)  ! nc
  real(r8) :: dumi(mgncol,nlev)   ! qi
  real(r8) :: dumqitend(mgncol,nlev)
  real(r8) :: rainfrztend(mgncol,nlev)
  real(r8) :: icemelttend(mgncol,nlev)
  real(r8) :: dumni(mgncol,nlev)  ! ni
  real(r8) :: dumr(mgncol,nlev)   ! rain mixing ratio
  real(r8) :: dumnr(mgncol,nlev)  ! rain number concentration
!++ag
  real(r8) :: dumg(mgncol,nlev)   ! graupel mixing ratio
  real(r8) :: dumng(mgncol,nlev)  ! graupel number concentration
!--ag
  ! Array dummy variable
  real(r8) :: dum_2D(mgncol,nlev)
  real(r8) :: pdel_inv(mgncol,nlev)

  ! loop array variables
  ! "i" and "k" are column/level iterators for internal (MG) variables
  ! "n" is used for other looping (currently just sedimentation)
  integer i, k, n

  ! number of sub-steps for loops over "n" (for sedimentation)
  integer nstep
  integer mdust

  ! Varaibles to scale fall velocity between small and regular ice regimes.
  real(r8) :: irad
  real(r8) :: ifrac

!++kt
  ! parameters for lookup table
  real(r8)    :: dum11,dum22,dum5,dum4, dumlr, kttmp
  integer :: dumit, dumk,dumj, dumtt
  real(r8), parameter :: thrd = 1._r8/3._r8
  ! interpolated quantities from ice lookup table
  real(r8) :: f1pr1,f1pr2,f1pr3,f1pr4,f1pr5,f1pr6,f1pr7,f1pr8,f1pr9,f1pr10,f1pr13,     &
       f1pr14,f1pr15,f1pr16, f1pr12
  ! values for lams and n0s
  real(r8) :: f1pr11lams, f1pr17n0s
  real(r8) :: af1pr11lams(mgncol,nlev), af1pr17n0s(mgncol,nlev)

  real(r8) :: af1pr1(mgncol,nlev),af1pr2(mgncol,nlev),af1pr3(mgncol,nlev),af1pr4(mgncol,nlev),&
       af1pr5(mgncol,nlev),af1pr6(mgncol,nlev),af1pr7(mgncol,nlev),af1pr8(mgncol,nlev),af1pr9(mgncol,nlev), &
       af1pr10(mgncol,nlev),af1pr13(mgncol,nlev),     &
       af1pr14(mgncol,nlev),af1pr15(mgncol,nlev),af1pr16(mgncol,nlev), &
       af1pr12(mgncol,nlev)

!--kt

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! Return error message
  errstring = ' '

  ! Process inputs

  ! assign variable deltat to deltatin
  deltat = deltatin

  ! Copies of input concentrations that may be changed internally.
  qc = qcn
  nc = ncn
!++kt MG4 combines snow and ice... might be better to remove qsin and nsn from the
!++kt    parameter list
!  qi = qin+qsn
!  ni = nin+nsn
  qi = qin
  ni = nin
!--kt
  qr = qrn
  nr = nrn
!++ag
  qg = qgr
  ng = ngr

  ! cldn: used to set cldm, unused for subcolumns
  ! liqcldf: used to set lcldm, unused for subcolumns
  ! icecldf: used to set icldm, unused for subcolumns

  if (microp_uniform) then
     ! subcolumns, set cloud fraction variables to one
     ! if cloud water or ice is present, if not present
     ! set to mincld (mincld used instead of zero, to prevent
     ! possible division by zero errors).

     where (qc >= qsmall)
        lcldm = 1._r8
     elsewhere
        lcldm = mincld
     end where

     where (qi >= qsmall)
        icldm = 1._r8
     elsewhere
        icldm = mincld
     end where

     cldm = max(icldm, lcldm)
     qsfm = 1._r8

  else
     ! get cloud fraction, check for minimum
     cldm = max(cldn,mincld)
     lcldm = max(liqcldf,mincld)
     qsfm = qsatfac
!++kt
     icldm = cldm
     do k=1,nlev
        ! calculate precip fraction based on maximum overlap assumption
        ! if rain or snow mix ratios are smaller than threshold,
        ! then leave cldmax as cloud fraction at current level
        if (k /= 1) then
           where (qr(:,k-1) >= qsmall .or. qi(:,k-1) >= qsmall)
              icldm(:,k)=max(icldm(:,k-1),icldm(:,k))
           end where
        end if
      end do
!--kt
  end if

  ! Initialize local variables

  ! local physical properties
  rho = p/(r*t)
  dv = 8.794E-5_r8 * t**1.81_r8 / p
  mu = 1.496E-6_r8 * t**1.5_r8 / (t + 120._r8)
  sc = mu/(rho*dv)

  ! air density adjustment for fallspeed parameters
  ! includes air density correction factor to the
  ! power of 0.54 following Heymsfield and Bansemer 2007

  rhof=(rhosu/rho)**0.54_r8

  arn=ar*rhof
  asn=as*rhof
!++ag if do hail then agn = ah *rhof else ag*rhof
  if (do_hail) then
     agn = ah * rhof
  end if
  if (do_graupel) then
     agn=ag*rhof
  end if
!--ag
  acn=g*rhow/(18._r8*mu)
  ain=ai*(rhosu/rho)**0.35_r8
  ajn=aj*(rhosu/rho)**0.35_r8

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Get humidity and saturation vapor pressures

  do k=1,nlev
     do i=1,mgncol

        call qsat_water(t(i,k), p(i,k), esl(i,k), qvl(i,k))

        ! make sure when above freezing that esi=esl, not active yet
        if (t(i,k) >= tmelt) then
           esi(i,k)=esl(i,k)
           qvi(i,k)=qvl(i,k)
        else
           call qsat_ice(t(i,k), p(i,k), esi(i,k), qvi(i,k))

           ! Scale the water saturation values to reflect subgrid scale
           ! ice cloud fraction, where ice clouds begin forming at a
           ! gridbox average relative humidity of rhmini (not 1).
           !
           ! NOTE: For subcolumns and other non-subgrid clouds, qsfm willi
           ! be 1.
           qvi(i,k) = qsfm(i,k) * qvi(i,k)
           esi(i,k) = qsfm(i,k) * esi(i,k)
           qvl(i,k) = qsfm(i,k) * qvl(i,k)
           esl(i,k) = qsfm(i,k) * esl(i,k)
        end if

     end do
  end do

  relhum = q / max(qvl, qsmall)

  !===============================================

  ! set mtime here to avoid answer-changing
  mtime=deltat

  ! initialize microphysics output
  qcsevap=0._r8
  qisevap=0._r8
  qvres  =0._r8
  cmeitot =0._r8
  vtrmc =0._r8
  vtrmi =0._r8
  qcsedten =0._r8
  qisedten =0._r8
  qrsedten =0._r8 
  qssedten =0._r8 !++kt This is not used, should be removed from parameters list
!++ag
  qgsedten =0._r8
!--ag

  pratot=0._r8
  prctot=0._r8
  mnuccctot=0._r8
  mnuccttot=0._r8
  msacwitot=0._r8
  psacwstot=0._r8
  npsacwstot=0._r8
  bergstot=0._r8
  bergtot=0._r8
  melttot=0._r8
  homotot=0._r8
  qcrestot=0._r8
  prcitot=0._r8
  praitot=0._r8
  qirestot=0._r8
  mnuccrtot=0._r8
!++ag
  mnuccritot=0._r8
!--ag
  pracstot=0._r8
  meltsdttot=0._r8
  frzrdttot=0._r8
  mnuccdtot=0._r8
!++ag
  psacrtot=0._r8
  pracgtot=0._r8
  psacwgtot=0._r8
  pgsacwtot=0._r8
  pgracstot=0._r8
  prdgtot=0._r8
!  eprdgtot=0._r8
  qmultgtot=0._r8
  qmultrgtot=0._r8
  npracgtot =0._r8
  nscngtot =0._r8
  ngracstot=0._r8
  nmultgtot=0._r8
  nmultrgtot=0._r8
  npsacwgtot=0._r8
!need to zero these out to be totally switchable (for conservation)
  psacr=0._r8
  pracg=0._r8
  psacwg=0._r8
  pgsacw=0._r8
  pgracs=0._r8
  prdg=0._r8
!  eprdg=0._r8
  qmultg=0._r8
  qmultrg=0._r8
  npracg=0._r8
  nscng=0._r8
  ngracs=0._r8
  nmultg=0._r8
  nmultrg=0._r8
  npsacwg=0._r8
!--ag

  rflx=0._r8
  sflx=0._r8
  lflx=0._r8
  iflx=0._r8
!++ag
  gflx=0._r8
  bgtmp=0._r8
  rhogtmp=0._r8
!--ag

  ! initialize precip output

  qrout=0._r8
  qsout=0._r8
  nrout=0._r8
  nsout=0._r8
!++ag
  qgout=0._r8
  ngout=0._r8
!--ag

  ! for refl calc
  rainrt = 0._r8

  ! initialize rain size
  rercld=0._r8

  qcsinksum_rate1ord = 0._r8

  ! initialize variables for trop_mozart
  nevapr = 0._r8
  prer_evap = 0._r8
  evapsnow = 0._r8
  am_evp_st = 0._r8
  prain = 0._r8
  prodsnow = 0._r8
  cmeout = 0._r8

  precip_frac = mincld

  lamc=0._r8

  ! initialize microphysical tendencies

  tlat=0._r8
  qvlat=0._r8
  qctend=0._r8
  qitend=0._r8
  qitendresid=0._r8
  rainfrztend=0._r8
  icemelttend=0._r8
  qstend = 0._r8
  qrtend = 0._r8
  nctend=0._r8
  nitend=0._r8
  nrtend = 0._r8
  nstend = 0._r8
!++ag
  qgtend = 0._r8
  ngtend = 0._r8
!--ag

  ! initialize in-cloud and in-precip quantities to zero
  qcic  = 0._r8
  qiic  = 0._r8
  qsic  = 0._r8
  qric  = 0._r8
!++ag
  qgic  = 0._r8
!--ag

  ncic  = 0._r8
  niic  = 0._r8
  nsic  = 0._r8
  nric  = 0._r8
!++ag
  ngic  = 0._r8
!--ag

  ! initialize precip at surface

  prect = 0._r8
  preci = 0._r8

  ! initialize precip fallspeeds to zero
  ums = 0._r8  
  umr = 0._r8
  unr = 0._r8

!++trude
  n0s = 0._r8
  lams = 0._r8
!--trude
!++ag
  umg = 0._r8
  ung = 0._r8
!--ag

  ! initialize limiter for output
  qcrat = 1._r8

  ! Many outputs have to be initialized here at the top to work around
  ! ifort problems, even if they are always overwritten later.
  effc = 10._r8
  lamcrad = 0._r8
  pgamrad = 0._r8
  effc_fn = 10._r8
  effi = 25._r8
  !++ trude
  effi = effi*micro_mg_effi_factor
  !-- trude
  sadice = 0._r8
  sadsnow = 0._r8
  deffi = 50._r8

  qrout2 = 0._r8
  nrout2 = 0._r8
  drout2 = 0._r8
  qsout2 = 0._r8
  nsout2 = 0._r8
  dsout = 0._r8
  dsout2 = 0._r8
!++ag 
  qgout2 =0._r8
  ngout2= 0._r8
  freqg = 0._r8
!--ag

  freqr = 0._r8
  freqs = 0._r8

  reff_rain = 0._r8
  reff_snow = 0._r8
!++ag
  reff_grau = 0._r8
!--ag

  refl = -9999._r8
  arefl = 0._r8
  areflz = 0._r8
  frefl = 0._r8
  csrfl = 0._r8
  acsrfl = 0._r8
  fcsrfl = 0._r8

  ncal = 0._r8
  ncai = 0._r8

  nfice = 0._r8

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! droplet activation
  ! get provisional droplet number after activation. This is used for
  ! all microphysical process calculations, for consistency with update of
  ! droplet mass before microphysics

  ! calculate potential for droplet activation if cloud water is present
  ! tendency from activation (npccn) is read in from companion routine

  ! output activated liquid and ice (convert from #/kg -> #/m3)
  !--------------------------------------------------
  where (qc >= qsmall)
     nc = max(nc + npccn*deltat, 0._r8)
     ncal = nc*rho/lcldm ! sghan minimum in #/cm3
  elsewhere
     ncal = 0._r8
  end where

  where (t < icenuct)
     ncai = naai*rho
  elsewhere
     ncai = 0._r8
  end where

  !===============================================

  ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
  !
  ! NOTE: If using gridbox average values, condensation will not occur until rh=1,
  ! so the threshold seems like it should be 1.05 and not rhmini + 0.05. For subgrid
  ! clouds (using rhmini and qsfacm), the relhum has already been adjusted, and thus
  ! the nucleation threshold should also be 1.05 and not rhmini + 0.05.

  !-------------------------------------------------------

  if (do_cldice) then
     where (naai > 0._r8 .and. t < icenuct .and. &
          relhum*esl/esi > 1.05_r8)

        !if NAAI > 0. then set numice = naai (as before)
        !note: this is gridbox averaged
        nnuccd = (naai-ni/icldm)/mtime*icldm
        nnuccd = max(nnuccd,0._r8)
        nimax = naai*icldm

        !Calc mass of new particles using new crystal mass...
        !also this will be multiplied by mtime as nnuccd is...

        mnuccd = nnuccd * mi0

     elsewhere
        nnuccd = 0._r8
        nimax = 0._r8
        mnuccd = 0._r8
     end where

  end if

  !KTC Removed instantaneous adjustment code

  do k=1,nlev
    do i=1,mgncol
        ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
        !-------------------------------------------------------
        ! for microphysical process calculations
        ! units are kg/kg for mixing ratio, 1/kg for number conc

        if (qc(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           qcic(i,k)=min(qc(i,k)/lcldm(i,k),5.e-3_r8)
           ncic(i,k)=max(nc(i,k)/lcldm(i,k),0._r8)

           ! specify droplet concentration
           if (nccons) then
              ncic(i,k)=ncnst/rho(i,k)
           end if
        else
           qcic(i,k)=0._r8
           ncic(i,k)=0._r8
        end if

        if (qi(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           qiic(i,k)=min(qi(i,k)/icldm(i,k),5.e-3_r8)
           niic(i,k)=max(ni(i,k)/icldm(i,k),0._r8)
!++kt
           ! mg4 max niic = 500 L-1
           if (niic(i,k).ge.1.e-20_r8) then	
!              dum=500.e+3_r8/(rho(i,k)*niic(i,k))
              dum=micro_mg_max_nicons/(rho(i,k)*niic(i,k))
              niic(i,k)=niic(i,k)*min(dum,1._r8)
!++trude The assignment below was added to MG2 single ice category. 
!              However, by changing ni(i,k) inside this 
!              routine I belive we create problems regarding ice number concentration
!              ni(i,k) = niic(i,k)*icldm(i,k)
!--trude
           endif
!--kt
           ! switch for specification of cloud ice number
           if (nicons) then
              niic(i,k)=ninst/rho(i,k)
           end if
        else
           qiic(i,k)=0._r8
           niic(i,k)=0._r8
        end if

     end do
  end do

  !========================================================================

  ! for sub-columns cldm has already been set to 1 if cloud
  ! water or ice is present, so precip_frac will be correctly set below
  ! and nothing extra needs to be done here

  precip_frac = cldm
  dumi = qiic
  dumni = niic

  micro_vert_loop: do k=1,nlev
      if(mgncol.gt.0) then
         if (trim(micro_mg_precip_frac_method) == 'in_cloud') then

            if (k /= 1) then
               where (qc(:,k) < qsmall .and. qi(:,k) < qsmall)
                  precip_frac(:,k) = precip_frac(:,k-1)
               end where
            endif

         else if (trim(micro_mg_precip_frac_method) == 'max_overlap') then

!++ag  add graupel to precip frac?

        ! calculate precip fraction based on maximum overlap assumption

        ! if rain or snow mix ratios are smaller than threshold,
        ! then leave precip_frac as cloud fraction at current level
            if (k /= 1) then
!++ag
!           where (qr(:,k-1) >= qsmall .or. qs(:,k-1) >= qsmall .or. qg(:,k-1) >= qsmall)
!--ag
!++kt   ! mg4 use qi instead of qs
               where (qr(:,k-1) >= qsmall .or. qi(:,k-1) >= qsmall)
                  precip_frac(:,k)=max(precip_frac(:,k-1),precip_frac(:,k))
               end where
!--kt
            end if

         endif ! max_overlap


     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ! get size distribution parameters based on in-cloud cloud water
     ! these calculations also ensure consistency between number and mixing ratio
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     ! cloud liquid
     !-------------------------------------------

          call size_dist_param_liq(mg_liq_props, qcic(:,k), ncic(:,k), rho(:,k), &
                pgam(:,k), lamc(:,k), mgncol)

     !========================================================================
     ! autoconversion of cloud liquid water to rain
     ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
     ! minimum qc of 1 x 10^-8 prevents floating point error

           if (.not. do_sb_physics) then
              call kk2000_liq_autoconversion(microp_uniform, qcic(:,k), &
              ncic(:,k), rho(:,k), relvar(:,k), prc(:,k), nprc(:,k), nprc1(:,k), micro_mg_autocon_fact, micro_mg_autocon_exp, mgncol)
           endif

     ! assign qric based on prognostic qr, using assumed precip fraction
     ! note: this could be moved above for consistency with qcic and qiic calculations
           qric(:,k) = qr(:,k)/precip_frac(:,k)
           nric(:,k) = nr(:,k)/precip_frac(:,k)

     ! limit in-precip mixing ratios to 10 g/kg
           qric(:,k)=min(qric(:,k),0.01_r8)

     ! add autoconversion to precip from above to get provisional rain mixing ratio
     ! and number concentration (qric and nric)

           where (qric(:,k).lt.qsmall)
              qric(:,k)=0._r8
              nric(:,k)=0._r8
           end where

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

           nric(:,k)=max(nric(:,k),0._r8)
	  
     ! Alternative autoconversion 
           if (do_sb_physics) then
              call sb2001v2_liq_autoconversion(pgam(:,k),qcic(:,k),ncic(:,k), &
              qric(:,k),rho(:,k),relvar(:,k),prc(:,k),nprc(:,k),nprc1(:,k), mgncol)     
           endif

!++kt get ice processrates from 3D lookuptable
           do i = 1, mgncol

        ! find indices in 4D ice lookup table
        !------------------------------------------------------------------------------------------!
        ! ice lookup table for ice micro processes
         if (qiic(i,k).ge.qsmall) then  
            f1pr1 = 0._r8
            f1pr2 = 0._r8
            f1pr3 = 0._r8
            f1pr4 = 0._r8
            f1pr5 = 0._r8
            f1pr6 = 0._r8
            f1pr7 = 0._r8
            f1pr8 = 0._r8
            f1pr9 = 0._r8
            f1pr10 = 0._r8
            f1pr13 = 0._r8
            f1pr14 = 0._r8
            f1pr15 = 0._r8
            f1pr16 = 0._r8

           ! find index for qi (total ice mass mixing ratio)
            dum11 = (dlog10(qiic(i,k))+16._r8)*1.41328_r8
            dumit = int(dum11)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum11 = min(dum11,real(isize))
            dum11 = min(dum11,dble(isize))
            dum11 = max(dum11,1._r8)
            dumit = max(1,dumit)
            dumit = min(isize-1,dumit)

           ! find index for Ni (ice number mixing ratio)
            dum22 = (dlog10(niic(i,k))+10._r8)*1.10731_r8
            dumk = int(dum22)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum22 = min(dum22,real(jsize))
            dum22 = min(dum22,dble(jsize))
            dum22 = max(dum22,1._r8)
            dumk = max(1,dumk)
            dumk = min(jsize-1,dumk)

           ! find index for scaled mean rain size
           ! if no rain, then just choose dumj = 1 and don't calculate rain-ice
           ! collection processes

            if (qric(i,k).ge.qsmall) then
              ! calculate scaled mean size for consistency with ice lookup table
               dumlr = (qric(i,k)/(pi*rhow*nric(i,k)))**thrd
               dum5  = (dlog10(1._r8*dumlr)+5._r8)*10.70415_r8
               dumj  = int(dum5)
              ! set limits
               dum5  = min(dum5,dble(rcollsize))
               dum5  = max(dum5,1._r8)
               dumj  = max(1,dumj)
               dumj  = min(rcollsize-1,dumj)
            else
               dumj  = 1
               dum5  = 1._r8
            endif

           ! find index for temperature
            ! Temp ranges: -65 - -55,  -55 - -40, -40 - -20
           ! Midpoint :           -60              -47.5        -30

            dum4  = t(i,k)      ! trude, and something else here for finding index qitot(i,k,iice)*3. + 1.
            if (t(i,k) .lt. 213.15_r8) then ! (-60)
               dum4 = 3         !trude, need to calculate the "interger" value
            else if (t(i,k).ge. 243.15_r8) then  
               dum4 = 1         ! (-30) 
            else
               dum4 = 17.21_r8-(t(i,k)/15._r8) 
            endif
            dumtt = int(dum4)
           ! set limits
!             dum4  = min(dum4,real(tsize))
            dum4  = min(dum4,dble(tsize))
            dum4  = max(dum4,1._r8)
            dumtt = max(1,dumtt)
            dumtt = min(tsize-1,dumtt)

           ! call to lookup table interpolation subroutines to get process rates
            call access_lookup_table(dumtt,dumit,dumk,7,dum11,dum22,dum4,f1pr9)
            call access_lookup_table(dumtt,dumit,dumk,8,dum11,dum22,dum4,f1pr10)

           ! adjust Ni if needed to make sure mean size is in bounds (i.e., apply lambda limiters) 
            niic(i,k) = min(niic(i,k),f1pr9)
            niic(i,k) = max(niic(i,k),f1pr10)

!            niic(i,k) = min(niic(i,k),500.e3_r8/rho(i,k))
            niic(i,k) = min(niic(i,k),micro_mg_max_nicons/rho(i,k))

!++trude The assignment below was added to MG2 single ice category. However, by changing ni(i,k) inside this 
!                  routine I belive we create problems regarding ice number concentration
!            ni(i,k) = niic(i,k)*icldm(i,k)
!--trude

           ! find index for Ni (ice number mixing ratio)
            dum22 = (dlog10(niic(i,k))+10._r8)*1.10731_r8
            dumk = int(dum22)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum22 = min(dum22,real(jsize))
            dum22 = min(dum22,dble(jsize))
            dum22 = max(dum22,1._r8)
            dumk = max(1,dumk)
            dumk = min(jsize-1,dumk)

           ! call to lookup table interpolation subroutines to get process rates
            call access_lookup_table(dumtt,dumit,dumk,7,dum11,dum22,dum4,f1pr9)
            call access_lookup_table(dumtt,dumit,dumk,8,dum11,dum22,dum4,f1pr10)
            call access_lookup_table(dumtt,dumit,dumk,2,dum11,dum22,dum4,f1pr2)
            call access_lookup_table(dumtt,dumit,dumk,3,dum11,dum22,dum4,f1pr3)
            call access_lookup_table(dumtt,dumit,dumk,4,dum11,dum22,dum4,f1pr4)
            call access_lookup_table(dumtt,dumit,dumk,5,dum11,dum22,dum4,f1pr5)
            call access_lookup_table(dumtt,dumit,dumk,10,dum11,dum22,dum4,f1pr14)
            call access_lookup_table(dumtt,dumit,dumk,12,dum11,dum22,dum4,f1pr12)

            af1pr1(i,k) = f1pr1  ! mass weighed ice number concentration fall speed
            af1pr2(i,k) = f1pr2  ! mass weighted ice mixing ratio fall speed
            af1pr3(i,k) = f1pr3  ! ice self accregation
            af1pr4(i,k) = f1pr4  ! accretion of cloud water by ice
            af1pr5(i,k) = f1pr5  ! ice deposition/sublimation
            af1pr6(i,k) = f1pr6  ! effective radius
            af1pr7(i,k) = f1pr7  ! accrete rain by ice
            af1pr8(i,k) = f1pr8  ! accrete rain by ice
            af1pr9(i,k) = f1pr9  ! max ice number concentration 
            af1pr10(i,k) = f1pr10 ! min ice number concentration 
            af1pr13(i,k) = f1pr13 ! not used
            af1pr14(i,k) = f1pr14 ! ice deposition/sublimation
            af1pr15(i,k) = f1pr15 ! not used
            af1pr16(i,k) = f1pr16 ! not used
            af1pr12(i,k) = f1pr12 ! not used

         endif                  ! qitot > qsmall
      enddo                     ! end get process rates from lookup table

!--kt
!++ag kt If precip mix ratio is zero, so should number concentration for 
!++ graupel, which is assumed to be 'precip_frac'
      qgic(:,k) = qg(:,k)/precip_frac(:,k)
      ngic(:,k) = ng(:,k)/precip_frac(:,k)

     ! limit in-precip mixing ratios to 10 g/kg
      qgic(:,k)=min(qgic(:,k),0.01_r8)

     ! if precip mix ratio is zero so should number concentration
      where (qgic(:,k) < qsmall)
         qgic(:,k)=0._r8
         ngic(:,k)=0._r8
      end where

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

      ngic(:,k)=max(ngic(:,k),0._r8)

!--ag    

     !.......................................................................
     ! get size distribution parameters for precip
     !......................................................................
     ! rain

      call size_dist_param_basic(mg_rain_props, qric(:,k), nric(:,k), &
      lamr(:,k), mgncol, n0=n0r(:,k))

      where (lamr(:,k) >= qsmall)

        ! provisional rain number and mass weighted mean fallspeed (m/s)

         unr(:,k) = min(arn(:,k)*gamma_br_plus1/lamr(:,k)**br,9.1_r8*rhof(:,k))
         umr(:,k) = min(arn(:,k)*gamma_br_plus4/(6._r8*lamr(:,k)**br),9.1_r8*rhof(:,k))

      elsewhere
         umr(:,k) = 0._r8
         unr(:,k) = 0._r8
      end where

!++kt
     !......................................................................
     ! "snow" - updated for mg4 as falling ice
      ums(:,k) = 0._r8

!++trude  updated dumni
      dumi(:,k) = qiic(:,k)
      dumni(:,k) = niic(:,k)
!--trude     
      do i = 1, mgncol

        if ((dumi(i,k).ge.qsmall).and.(dumni(i,k).ge.icsmall)) then
           !! ++ trude. Since ice has been updated, call access_lookup_table again for fall speed.            
           ! find index for qi (total ice mass mixing ratio)
           dum11 = (dlog10(dumi(i,k))+16._r8)*1.41328_r8
           dumit = int(dum11)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum11 = min(dum11,real(isize))
           dum11 = min(dum11,dble(isize))
           dum11 = max(dum11,1._r8)
           dumit = max(1,dumit)
           dumit = min(isize-1,dumit)

           ! find index for Ni (ice number mixing ratio)
           dum22 = (dlog10(dumni(i,k))+10._r8)*1.10731_r8
           dumk = int(dum22)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum22 = min(dum22,real(jsize))
           dum22 = min(dum22,dble(jsize))
           dum22 = max(dum22,1._r8)
           dumk = max(1,dumk)
           dumk = min(jsize-1,dumk)
           ! find index for temperature
           dum4  = t(i,k)   ! trude, and something else here for finding index qitot(i,k,iice)*3. + 1.
           if (t(i,k) .lt. 213.15_r8) then    ! (-60)
              dum4 = 3     !trude, need to calculate the "interger" value
           else if (t(i,k).ge. 243.15_r8) then  
              dum4 = 1                      ! (-30) 
           else
              dum4 = 17.21_r8-(t(i,k)/15._r8) 
           endif
           dumtt = int(dum4)
!     set limits
           dum4  = min(dum4,dble(tsize))
           dum4  = max(dum4,1._r8)
           dumtt = max(1,dumtt)
           dumtt = min(tsize-1,dumtt)

           call access_lookup_table(dumtt,dumit,dumk,2,dum11,dum22,dum4,f1pr2)
           ! lams
           call access_lookup_table(dumtt,dumit,dumk,13,dum11,dum22,dum4,f1pr11lams)
           ! n0s
           call access_lookup_table(dumtt,dumit,dumk,14,dum11,dum22,dum4,f1pr17n0s)
           
            lams(i,k) = f1pr11lams
! ++ trude, to test changes to fall speed
            f1pr2=f1pr2*micro_mg_vtrmi_factor
! -- trude
            n0s(i,k)  = f1pr17n0s
            n0sout(i,k) = n0s(i,k)
            vtrmi(i,k) = min(f1pr2*rhof(i,k),1.2_r8*rhof(i,k))
            ums(i,k) = vtrmi(i,k)
           
        endif
      enddo                     !mgncol
!--kt
!++ag
!use correct bg or bh  (bgtmp=bg or bh)
     !......................................................................
     ! graupel/hail

!++AG SET rhog here and for mg_graupel_props?
! For now: rhog is constant. Set to same in micro_mg_utils.F90
! Ideally: find a method to set once. (Hail = 400, Graupel = 500 from M2005)

     if (do_hail) then 
        bgtmp = bh 
        rhogtmp = rhoh
     end if
     if (do_graupel) then 
        bgtmp = bg
        rhogtmp = rhog
     end if

!mg,snow_props or mg_graupel props?

     call size_dist_param_basic(mg_graupel_props, qgic(:,k), ngic(:,k), &
          lamg(:,k), mgncol, n0=n0g(:,k))

     where (lamg(:,k) > 0._r8)

        ! provisional graupel/hail number and mass weighted mean fallspeed (m/s)
        !++KTC note: if do_hail and do_graupel are both false, why is lamg > 0??
        umg(:,k) = min(agn(:,k)*gamma(4._r8+bgtmp)/(6._r8*lamg(:,k)**bgtmp),20._r8*rhof(:,k))
        ung(:,k) = min(agn(:,k)*gamma(1._r8+bgtmp)/lamg(:,k)**bgtmp,20._r8*rhof(:,k))

     elsewhere
        umg(:,k) = 0._r8
        ung(:,k) = 0._r8
     end where
!--ag


     if (do_cldice) then
        if (.not. use_hetfrz_classnuc) then

           ! heterogeneous freezing of cloud water
           !----------------------------------------------

           call immersion_freezing(microp_uniform, t(:,k), pgam(:,k), lamc(:,k), &
                qcic(:,k), ncic(:,k), relvar(:,k), mnuccc(:,k), nnuccc(:,k), mgncol)

           ! make sure number of droplets frozen does not exceed available ice nuclei concentration
           ! this prevents 'runaway' droplet freezing

           where (qcic(:,k).ge.qsmall .and. t(:,k).lt.269.15_r8)
              where (nnuccc(:,k)*lcldm(:,k).gt.nnuccd(:,k))
                 ! scale mixing ratio of droplet freezing with limit
                 mnuccc(:,k)=mnuccc(:,k)*(nnuccd(:,k)/(nnuccc(:,k)*lcldm(:,k)))
                 nnuccc(:,k)=nnuccd(:,k)/lcldm(:,k)
              end where
           end where

           mdust = size(rndst,3)
           call contact_freezing(microp_uniform, t(:,k), p(:,k), rndst(:,k,:), &
                nacon(:,k,:), pgam(:,k), lamc(:,k), qcic(:,k), ncic(:,k), &
                relvar(:,k), mnucct(:,k), nnucct(:,k), mgncol, mdust)

           mnudep(:,k)=0._r8
           nnudep(:,k)=0._r8

        else

           ! Mass of droplets frozen is the average droplet mass, except
           ! with two limiters: concentration must be at least 1/cm^-3, and
           ! mass must be at least the minimum defined above.
           mi0l = qcic(:,k)/max(ncic(:,k), 1.0e6_r8/rho(:,k))
           mi0l = max(mi0l_min, mi0l)

           where (qcic(:,k) >= qsmall)
              nnuccc(:,k) = frzimm(:,k)*1.0e6_r8/rho(:,k)
              mnuccc(:,k) = nnuccc(:,k)*mi0l

              nnucct(:,k) = frzcnt(:,k)*1.0e6_r8/rho(:,k)
              mnucct(:,k) = nnucct(:,k)*mi0l

              nnudep(:,k) = frzdep(:,k)*1.0e6_r8/rho(:,k)
              mnudep(:,k) = nnudep(:,k)*mi0
           elsewhere
              nnuccc(:,k) = 0._r8
              mnuccc(:,k) = 0._r8

              nnucct(:,k) = 0._r8
              mnucct(:,k) = 0._r8

              nnudep(:,k) = 0._r8
              mnudep(:,k) = 0._r8
           end where

        end if

     else
        mnuccc(:,k)=0._r8
        nnuccc(:,k)=0._r8
        mnucct(:,k)=0._r8
        nnucct(:,k)=0._r8
        mnudep(:,k)=0._r8
        nnudep(:,k)=0._r8
     end if

!++kt
     call ice_self_aggregation(t(:,k), rho(:,k), rhof(:,k), af1pr3(:,k), qiic(:,k), niagg(:,k))
     niagg(:,k) = -niagg(:,k)
 
!++ trude, test for changes in ice collection cloud water
    af1pr4(:,k)=af1pr4(:,k)*micro_mg_iaccr_factor
!--trude

    call accrete_cloud_water_ice(t(:,k), rho(:,k), rhof(:,k), af1pr4(:,k), &
          qcic(:,k), ncic(:,k), qiic(:,k), &
          psacws(:,k), npsacws(:,k))

 
!--kt

     if (do_cldice) then
        call secondary_ice_production(t(:,k), psacws(:,k), msacwi(:,k), nsacwi(:,k), mgncol)
     else
        nsacwi(:,k) = 0.0_r8
        msacwi(:,k) = 0.0_r8
     end if

!++kt
     !  We do the access_lookup_table_coll here, because nric is modified in  size_dist_param_basic

     af1pr7(:,k) = 0._r8
     af1pr8(:,k) = 0._r8

     do i = 1, mgncol
        if(qiic(i,k).ge.qsmall) then
           !     find index for qi (total ice mass mixing ratio)
           dum11 = (dlog10(qiic(i,k))+16._r8)*1.41328_r8
           dumit = int(dum11)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum11 = min(dum11,real(isize))
           dum11 = min(dum11,dble(isize))
           dum11 = max(dum11,1._r8)
           dumit = max(1,dumit)
           dumit = min(isize-1,dumit)

           ! find index for Ni (ice number mixing ratio)
           dum22 = (dlog10(niic(i,k))+10._r8)*1.10731_r8
           dumk = int(dum22)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum22 = min(dum22,real(jsize))
           dum22 = min(dum22,dble(jsize))
           dum22 = max(dum22,1._r8)
           dumk = max(1,dumk)
           dumk = min(jsize-1,dumk)

           if (qric(i,k).ge.qsmall) then
              ! calculate scaled mean size for consistency with ice lookup table
              dumlr = (qric(i,k)/(pi*rhow*nric(i,k)))**thrd
              dum5  = (dlog10(1._r8*dumlr)+5._r8)*10.70415_r8
              dumj  = int(dum5)
              ! set limits
              
              dum5  = min(dum5,dble(rcollsize))
              dum5  = max(dum5,1._r8)
              dumj  = max(1,dumj)
              dumj  = min(rcollsize-1,dumj)
           else
              dumj  = 1
              dum5  = 1._r8
           endif
           ! find index for temperature
           dum4  = t(i,k)       ! trude, and something else here for finding index qitot(i,k,iice)*3. + 1.
           if (t(i,k) .lt. 213.15_r8) then ! (-60)
              dum4 = 3          !trude, need to calculate the "interger" value
           else if (t(i,k).ge. 243.15_r8) then  
              dum4 = 1          ! (-30) 
           else
              dum4 = 17.21_r8-(t(i,k)/15._r8) 
           endif
           dumtt = int(dum4)
           ! set limits
           dum4  = min(dum4,dble(tsize))
           dum4  = max(dum4,1._r8)
           dumtt = max(1,dumtt)
           dumtt = min(tsize-1,dumtt)
           ! ice-rain collection processes
           if (qric(i,k).ge.qsmall) then
              call access_lookup_table_coll(dumtt,dumj,dumit,dumk,1,dum11,dum22,    &
              dum5,dum4,f1pr7)
              call access_lookup_table_coll(dumtt,dumj,dumit,dumk,2,dum11,dum22,    &
              dum5,dum4,f1pr8)
           else
              f1pr7 = 0._r8
              f1pr8 = 0._r8
           endif
           af1pr7(i,k) = f1pr7
           af1pr8(i,k) = f1pr8
        end if ! qiic > qsmall
     end do

!++ trude, test for changes in ice collection rain
    af1pr7(:,k)=af1pr7(:,k)*micro_mg_iaccr_factor
    af1pr8(:,k)=af1pr8(:,k)*micro_mg_iaccr_factor
!--trude

     call accrete_rain_ice(t(:,k), rho(:,k), rhof(:,k), af1pr8(:,k), af1pr7(:,k), &
          qric(:,k), qiic(:,k), n0r(:,k), pracs(:,k), npracs(:,k))
!--kt

     call heterogeneous_rain_freezing(t(:,k), qric(:,k), nric(:,k), lamr(:,k), &
          mnuccr(:,k), nnuccr(:,k), mgncol)

!     mnuccr(:,k) =  0._r8

     if (do_sb_physics) then
       call sb2001v2_accre_cld_water_rain(qcic(:,k), ncic(:,k), qric(:,k), &
            rho(:,k), relvar(:,k), pra(:,k), npra(:,k), mgncol)     
     else
       call accrete_cloud_water_rain(microp_uniform, qric(:,k), qcic(:,k), &
            ncic(:,k), relvar(:,k), accre_enhan(:,k), pra(:,k), npra(:,k), mgncol)
     endif
!++ trude
     pra(:,k)=pra(:,k)*micro_mg_accre_enhan_fact
     npra(:,k)=npra(:,k)*micro_mg_accre_enhan_fact
! -- trude

     call self_collection_rain(rho(:,k), qric(:,k), nric(:,k), nragg(:,k), mgncol)

!++ag Moved below graupel conditional, now two different versions
!     call evaporate_sublimate_precip(t(:,k), rho(:,k), &
!          dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k), &
!          lcldm(:,k), precip_frac(:,k), arn(:,k), asn(:,k), qcic(:,k), qiic(:,k), &
!          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
!          pre(:,k), prds(:,k), am_evp_st(:,k), mgncol)
!--ag

     !+++PMC 12/3/12 - NEW VAPOR DEP/SUBLIMATION GOES HERE!!!
     if (do_cldice) then
!++kt
        call ice_deposition_sublimation_mg4(t(:,k), q(:,k), qi(:,k), niic(:,k), &
             icldm(:,k), rho(:,k), dv(:,k), qvl(:,k), qvi(:,k), &
             berg(:,k), vap_dep(:,k), ice_sublim(:,k), &
             af1pr5(:,k),af1pr14(:,k), rhof(:,k), mu(:,k), sc(:,k), &
             mgncol)
!--kt
        berg(:,k)=berg(:,k)*micro_mg_berg_eff_factor

        where (vap_dep(:,k) < 0._r8 .and. qi(:,k) > qsmall .and. icldm(:,k) > mincld)
!++kt use niic instead of ni
           nsubi(:,k) = vap_dep(:,k) / qi(:,k) * niic(:,k) / icldm(:,k)
!--kt
        elsewhere
           nsubi(:,k) = 0._r8
        end where

        ! bergeron process should not reduce nc unless
        ! all ql is removed (which is handled elsewhere)
        !in fact, nothing in this entire file makes nsubc nonzero.
        nsubc(:,k) = 0._r8

     end if !do_cldice
     !---PMC 12/3/12

!++ag Process rate calls for graupel here.  
!   (Should this be in do_cldice loop?)
!===================================================================
!++kt Instead of snow, use ice
        !write(*,*) "KTC lams(:,k)=",lams(:,k)," k=",k
        !write(*,*) "KTC n0s(:,k)=",n0s(:,k)," k=",k

        lamsout(:,k) = lams(:,k)
        n0sout(:,k)  = n0s(:,k)
        lamrout(:,k) = lamr(:,k)
        n0rout(:,k)  = n0r(:,k)

     if(do_hail.or.do_graupel) then

        call graupel_collecting_snow(qiic(:,k),qric(:,k),umr(:,k),ums(:,k), &
             rho(:,k),lamr(:,k),n0r(:,k),lams(:,k),n0s(:,k), &
             psacr(:,k), mgncol)

       call graupel_collecting_cld_water(qgic(:,k),qcic(:,k),ncic(:,k),rho(:,k), &
             n0g(:,k),lamg(:,k),bgtmp,agn(:,k), psacwg(:,k), npsacwg(:,k), mgncol)
        
        call graupel_riming_liquid_snow(psacws(:,k),qiic(:,k),qcic(:,k),niic(:,k), &
             rho(:,k),rhosn,rhogtmp,asn(:,k),lams(:,k),n0s(:,k),deltat, &
             pgsacw(:,k),nscng(:,k),mgncol)

        call graupel_collecting_rain(qric(:,k),qgic(:,k),umg(:,k), &
             umr(:,k),ung(:,k),unr(:,k),rho(:,k),n0r(:,k),lamr(:,k),n0g(:,k), &
             lamg(:,k), pracg(:,k),npracg(:,k),mgncol)
        
!AG note: Graupel rain riming snow changes  
!    pracs, npracs, (accretion of rain by snow)  psacr (collection of snow by rain)

       call graupel_rain_riming_snow(pracs(:,k),npracs(:,k),psacr(:,k),qiic(:,k), &
             qric(:,k),nric(:,k),niic(:,k),n0s(:,k),lams(:,k),n0r(:,k),lamr(:,k), &
             deltat,pgracs(:,k),ngracs(:,k),mgncol)
       
        call graupel_rime_splintering(t(:,k),qcic(:,k),qric(:,k),qgic(:,k), &
             psacwg(:,k),pracg(:,k),qmultg(:,k),nmultg(:,k),qmultrg(:,k), &
             nmultrg(:,k),mgncol)

        call evaporate_sublimate_precip_graupel_mg4(t(:,k), rho(:,k), &
             dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k), &
             lcldm(:,k), precip_frac(:,k), arn(:,k), agn(:,k), bgtmp, &
             qcic(:,k), qiic(:,k), qric(:,k), qgic(:,k), &
             lamr(:,k), n0r(:,k), lamg(:,k), n0g(:,k), &
             pre(:,k), prdg(:,k), am_evp_st(:,k), mgncol)   

!!Not used: part of above
!!        call graupel_sublimate_evap(t(:,k),q(:,k),qgic(:,k),rho(:,k),n0g(:,k), &
!!             lamg(:,k),qvi(:,k),dv(:,k),mu(:,k),sc(:,k),bgtmp,agn(:,k), &
!!             prdg(:,k),eprdg(:,k),mgncol)


! Add a tentative check for process rates being positive....(cannot find output)?

!Check to make sure this works...
!        if (maxval(qmultg(:,k)).gt.0._r8) & 
!             write(iulog,*) "TEST, qmultg > 0 : max=",maxval(qmultg(:,k))
!        if (maxval(pgracs(:,k)).gt.0._r8) & 
!             write(iulog,*) "TEST, pgracs > 0 : max=",maxval(pgracs(:,k))

        if (minval(qmultg(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, qmultg < 0 : min=",minval(qmultg(:,k))

        if (minval(qmultrg(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, qmultrg < 0 : min=",minval(qmultrg(:,k))

        if (minval(pgracs(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, pgracs < 0 : min=",minval(pgracs(:,k))

        if (minval(psacwg(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, psacwg < 0 : min=",minval(psacwg(:,k))

        if (minval(npsacwg(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, npsacwg < 0 : min=",minval(npsacwg(:,k))

        if (minval(pracg(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, pracg < 0 : min=",minval(pracg(:,k))

        if (maxval(prdg(:,k)).gt.0._r8) & 
             write(iulog,*) "OOPS, prdg > 0 : max=",maxval(prdg(:,k))

        if (minval(nmultg(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, nmultg < 0 : min=",minval(nmultg(:,k))

        if (minval(nmultrg(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, nmultrg < 0 : min=",minval(nmultrg(:,k))

        if (minval(ngracs(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, ngracs < 0 : min=",minval(ngracs(:,k))

        if (minval(psacr(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, psacr < 0 : min=",minval(psacr(:,k))

        if (minval(nscng(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, nscng < 0 : min=",minval(nscng(:,k))

        if (minval(npracg(:,k)).lt.0._r8) & 
             write(iulog,*) "OOPS, npracg < 0 : min=",minval(npracg(:,k))


     else
! Routine without Graupel (original) - mg4: still use ice instead of snow        
        call evaporate_sublimate_precip_mg4(t(:,k), rho(:,k), &
          dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k), &
          lcldm(:,k), precip_frac(:,k), arn(:,k), qcic(:,k), qiic(:,k), &
          qric(:,k), lamr(:,k), n0r(:,k), &
          pre(:,k), am_evp_st(:,k), mgncol)

     end if ! end do_graupel/hail loop
!--kt
!--ag

     do i=1,mgncol

        ! conservation to ensure no negative values of cloud water/precipitation
        ! in case microphysical process rates are large
        !===================================================================

        ! note: for check on conservation, processes are multiplied by omsm
        ! to prevent problems due to round off error

        ! conservation of qc
        !-------------------------------------------------------------------

!++ag Add graupel tendencies for qc to equation ON
!        dum = ((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+ &
!             psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k))*deltat
!++kt remove bergs in mg4
        dum = ((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+ &
             psacws(i,k)+qmultg(i,k)+psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)+ &
             berg(i,k))*deltat 
!--ag

        if (dum.gt.qc(i,k)) then
!++ag
!           ratio = qc(i,k)/deltat/((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+ &
!                msacwi(i,k)+psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k))*omsm
           ratio = qc(i,k)/deltat/((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+ &
                msacwi(i,k)+psacws(i,k)+qmultg(i,k)+psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)+&
                berg(i,k))*omsm

           qmultg(i,k)=qmultg(i,k)*ratio
           psacwg(i,k)=psacwg(i,k)*ratio
           pgsacw(i,k)=pgsacw(i,k)*ratio
!--ag
           prc(i,k) = prc(i,k)*ratio
           pra(i,k) = pra(i,k)*ratio
           mnuccc(i,k) = mnuccc(i,k)*ratio
           mnucct(i,k) = mnucct(i,k)*ratio
           msacwi(i,k) = msacwi(i,k)*ratio
           psacws(i,k) = psacws(i,k)*ratio
           berg(i,k) = berg(i,k)*ratio
           qcrat(i,k) = ratio
        else
           qcrat(i,k) = 1._r8
        end if

        !PMC 12/3/12: ratio is also frac of step w/ liquid.
        !thus we apply berg for "ratio" of timestep and vapor
        !deposition for the remaining frac of the timestep.
        if (qc(i,k) >= qsmall) then
           vap_dep(i,k) = vap_dep(i,k)*(1._r8-qcrat(i,k))
        end if

     end do

     do i=1,mgncol

        !=================================================================
        ! apply limiter to ensure that ice/snow sublimation and rain evap
        ! don't push conditions into supersaturation, and ice deposition/nucleation don't
        ! push conditions into sub-saturation
        ! note this is done after qc conservation since we don't know how large
        ! vap_dep is before then
        ! estimates are only approximate since other process terms haven't been limited
        ! for conservation yet

        ! first limit ice deposition/nucleation vap_dep + mnuccd
        dum1 = vap_dep(i,k) + mnuccd(i,k)
        if (dum1 > 1.e-20_r8) then
           dum = (q(i,k)-qvi(i,k))/(1._r8 + xxls_squared*qvi(i,k)/(cpp*rv*t(i,k)**2))/deltat
           dum = max(dum,0._r8)
           if (dum1 > dum) then
              ! Allocate the limited "dum" tendency to mnuccd and vap_dep
              ! processes. Don't divide by cloud fraction; these are grid-
              ! mean rates.
              dum1 = mnuccd(i,k) / (vap_dep(i,k)+mnuccd(i,k))
              mnuccd(i,k) = dum*dum1
              vap_dep(i,k) = dum - mnuccd(i,k)
           end if
        end if

     end do

     do i=1,mgncol

        !===================================================================
        ! conservation of nc
        !-------------------------------------------------------------------
!++ag NEW ONE ON
!        dum = (nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+ &
!             npsacws(i,k)-nsubc(i,k))*lcldm(i,k)*deltat
        dum = (nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+ &
             npsacws(i,k)-nsubc(i,k)+npsacwg(i,k))*lcldm(i,k)*deltat
!--ag
        if (dum.gt.nc(i,k)) then
!++ag
!           ratio = nc(i,k)/deltat/((nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+&
!                npsacws(i,k)-nsubc(i,k))*lcldm(i,k))*omsm
           ratio = nc(i,k)/deltat/((nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+&
                npsacws(i,k)-nsubc(i,k)+npsacwg(i,k))*lcldm(i,k))*omsm
           npsacwg(i,k)=npsacwg(i,k)*ratio
!--ag
           nprc1(i,k) = nprc1(i,k)*ratio
           npra(i,k) = npra(i,k)*ratio
           nnuccc(i,k) = nnuccc(i,k)*ratio
           nnucct(i,k) = nnucct(i,k)*ratio
           npsacws(i,k) = npsacws(i,k)*ratio
           nsubc(i,k)=nsubc(i,k)*ratio
        end if

        mnuccri(i,k)=0._r8
        nnuccri(i,k)=0._r8

        if (do_cldice) then

           ! freezing of rain to produce ice if mean rain size is smaller than Dcs
           ! mg4 move all nuccr to nucci -> move all snow to ice
!++kt      !if (lamr(i,k) > qsmall .and. 1._r8/lamr(i,k) < Dcs) then
!++trude This test must be added again when we run with graupel
           if(do_hail.or.do_graupel) then
              if (lamr(i,k) > qsmall .and. 1._r8/lamr(i,k) < Dcs) then
                 mnuccri(i,k)=mnuccr(i,k)
                 nnuccri(i,k)=nnuccr(i,k)
                 mnuccr(i,k)=0._r8
                 nnuccr(i,k)=0._r8
              end if
           else 
              mnuccri(i,k)=mnuccr(i,k)
              nnuccri(i,k)=nnuccr(i,k)
              mnuccr(i,k)=0._r8
              nnuccr(i,k)=0._r8
           endif
        end if
!--trude

     end do

     do i=1,mgncol

        ! conservation of rain mixing ratio
        !-------------------------------------------------------------------
!++ag Implemented change for graupel
!        dum = ((-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k))*precip_frac(i,k)- &
!             (pra(i,k)+prc(i,k))*lcldm(i,k))*deltat

        dum = ((-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k) &
             +qmultrg(i,k)+pracg(i,k)+pgracs(i,k))*precip_frac(i,k)- &
             (pra(i,k)+prc(i,k))*lcldm(i,k))*deltat
!--ag

        ! note that qrtend is included below because of instantaneous freezing/melt

!++ag
!        if (dum.gt.qr(i,k).and. &
!             (-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k)).ge.qsmall) then
!           ratio = (qr(i,k)/deltat+(pra(i,k)+prc(i,k))*lcldm(i,k))/   &
!                precip_frac(i,k)/(-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k))*omsm
        if (dum.gt.qr(i,k).and. &
             (-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k)+qmultrg(i,k)+pracg(i,k)+pgracs(i,k)).ge.qsmall) then

           ratio = (qr(i,k)/deltat+(pra(i,k)+prc(i,k))*lcldm(i,k))/   &
                precip_frac(i,k)/(-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k) &
                +qmultrg(i,k)+pracg(i,k)+pgracs(i,k))*omsm

           qmultrg(i,k)= qmultrg(i,k)*ratio
           pracg(i,k)=pracg(i,k)*ratio
           pgracs(i,k)=pgracs(i,k)*ratio
!--ag
           pre(i,k)=pre(i,k)*ratio
           pracs(i,k)=pracs(i,k)*ratio
           mnuccr(i,k)=mnuccr(i,k)*ratio
           mnuccri(i,k)=mnuccri(i,k)*ratio
        end if

     end do

     do i=1,mgncol

        ! conservation of rain number
        !-------------------------------------------------------------------

        ! Add evaporation of rain number.
        if (pre(i,k) < 0._r8) then
           dum = pre(i,k)*deltat/qr(i,k)
           dum = max(-1._r8,dum)
           nsubr(i,k) = dum*nr(i,k)/deltat
        else
           nsubr(i,k) = 0._r8
        end if

     end do

     do i=1,mgncol

!++ag IMplemented change for graupel
!       dum = ((-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k))*precip_frac(i,k)- &
!             nprc(i,k)*lcldm(i,k))*deltat
        dum = ((-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k)+npracg(i,k)+ngracs(i,k)) &
             *precip_frac(i,k)- nprc(i,k)*lcldm(i,k))*deltat

!--ag

        if (dum.gt.nr(i,k)) then
!++ag
!            ratio = (nr(i,k)/deltat+nprc(i,k)*lcldm(i,k)/precip_frac(i,k))/ &
!                (-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k))*omsm
           ratio = (nr(i,k)/deltat+nprc(i,k)*lcldm(i,k))/precip_frac(i,k)/ &
                (-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k)+npracg(i,k)+ngracs(i,k))*omsm

           npracg(i,k)=npracg(i,k)*ratio
           ngracs(i,k)=ngracs(i,k)*ratio
!--ag
           nragg(i,k)=nragg(i,k)*ratio
           npracs(i,k)=npracs(i,k)*ratio
           nnuccr(i,k)=nnuccr(i,k)*ratio
           nsubr(i,k)=nsubr(i,k)*ratio
           nnuccri(i,k)=nnuccri(i,k)*ratio
        end if

     end do

     if (do_cldice) then

        do i=1,mgncol

           ! conservation of qi
           !-------------------------------------------------------------------

!++ag
!           dum = ((-mnuccc(i,k)-mnucct(i,k)-mnudep(i,k)-msacwi(i,k))*lcldm(i,k)+(prci(i,k)+ &
!                prai(i,k))*icldm(i,k)-mnuccri(i,k)*precip_frac(i,k) &
!                -ice_sublim(i,k)-vap_dep(i,k)-berg(i,k)-mnuccd(i,k))*deltat
!++kt
           dum = ((-mnuccc(i,k)-mnucct(i,k)-msacwi(i,k)-qmultg(i,k)-psacws(i,k))*lcldm(i,k)+ &
                (-mnudep(i,k))*icldm(i,k)+&            
                (-qmultrg(i,k)-mnuccri(i,k)-pracs(i,k))*precip_frac(i,k) &
                -ice_sublim(i,k)-vap_dep(i,k)-berg(i,k)-mnuccd(i,k))*deltat
!--kt
!-ag

           if (dum.gt.qi(i,k)) then
!++ag
!              ratio = (qi(i,k)/deltat+vap_dep(i,k)+berg(i,k)+mnuccd(i,k)+ &
!                   (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k)+ &
!                   mnuccri(i,k)*precip_frac(i,k))/ &
!                   ((prci(i,k)+prai(i,k))*icldm(i,k)-ice_sublim(i,k))*omsm
!++kt
              ratio = (qi(i,k)/deltat+vap_dep(i,k)+berg(i,k)+mnuccd(i,k)+ &
                   (mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+qmultg(i,k)+psacws(i,k))*lcldm(i,k)+ &
                   (qmultrg(i,k)+mnuccri(i,k)+pracs(i,k))*precip_frac(i,k)+mnudep(i,k)*icldm(i,k))/ &
                   (-ice_sublim(i,k))*omsm

! Only sink terms are limited. 
!              qmultg(i,k)=qmultg(i,k)*ratio
!              qmultrg(i,k)=qmultrg(i,k)*ratio
!--ag
              ice_sublim(i,k) = ice_sublim(i,k)*ratio
!--kt
           end if
!++kt
           ! KT: Adding in some "snow" conservation checks to ensure graupel processes are scaled correctly
!++ag
           if (do_hail .or. do_graupel) then
              !NOTE: mnuccr is moved to graupel when active
              !psacr is a positive value, but a loss for snow
              !HM: psacr is positive in dum (two negatives)

              dum = (-(pracs(i,k)-psacr(i,k))*precip_frac(i,k) &
                   -(psacws(i,k))*lcldm(i,k))*deltat
              
              if (dum.gt.qi(i,k).and.psacr(i,k).ge.qsmall) then

                 ratio = (qi(i,k)/deltat+ &
                      (psacws(i,k))*lcldm(i,k)+pracs(i,k)*precip_frac(i,k))/ &
                      precip_frac(i,k)/(psacr(i,k))*omsm
            
                 psacr(i,k)=psacr(i,k)*ratio
!--ag
              end if
           end if
!--kt
        end do

     end if

     if (do_cldice) then

        do i=1,mgncol

           ! conservation of ni
           !-------------------------------------------------------------------
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = 0._r8
           end if
!++ag
!           dum = ((-nnucct(i,k)-tmpfrz-nnudep(i,k)-nsacwi(i,k))*lcldm(i,k)+(nprci(i,k)+ &
!                nprai(i,k)-nsubi(i,k))*icldm(i,k)-nnuccri(i,k)*precip_frac(i,k)- &
!                nnuccd(i,k))*deltat
!++kt
!--ag
!++trude add graupel
           if (do_hail .or. do_graupel) then        
              dum = ((-nnucct(i,k)-tmpfrz-nsacwi(i,k)-nmultg(i,k))*lcldm(i,k)+&
                   (-nsubi(i,k)-niagg(i,k)-nnudep(i,k))*icldm(i,k)+(-nmultrg(i,k)+ngracs(i,k)-nnuccri(i,k))*precip_frac(i,k)- &
                   nnuccd(i,k)+nscng(i,k)*lcldm(i,k))*deltat
           else
              dum = ((-nnucct(i,k)-tmpfrz-nsacwi(i,k)-nmultg(i,k))*lcldm(i,k)+&
                   (-nsubi(i,k)-niagg(i,k)-nnudep(i,k))*icldm(i,k)+(-nmultrg(i,k)-nnuccri(i,k))*precip_frac(i,k)- &
                   nnuccd(i,k))*deltat
           end if
!--trude end add graupel

           if (dum.gt.ni(i,k)) then
!++ag
!              ratio = (ni(i,k)/deltat+nnuccd(i,k)+ &
!                   (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k)+ &
!                   nnuccri(i,k)*precip_frac(i,k))/ &
!                   ((nprci(i,k)+nprai(i,k)-nsubi(i,k))*icldm(i,k))*omsm
!++trude add graupel
              if (do_hail .or. do_graupel) then        
                 ratio = (ni(i,k)/deltat+nnuccd(i,k)+ &
                      (nnucct(i,k)+tmpfrz+nsacwi(i,k)+nmultg(i,k))*lcldm(i,k)+ &
                      (nnuccri(i,k)+nmultrg(i,k))*precip_frac(i,k)+nnudep(i,k)*icldm(i,k))/ &
                      (ngracs(i,k)*precip_frac(i,k)+nscng(i,k)*lcldm(i,k)+&
                      (-nsubi(i,k)-niagg(i,k))*icldm(i,k))*omsm
                 nscng(i,k)=nscng(i,k)*ratio
                 ngracs(i,k)=ngracs(i,k)*ratio
              else
                 ratio = (ni(i,k)/deltat+nnuccd(i,k)+ &
                      (nnucct(i,k)+tmpfrz+nsacwi(i,k)+nmultg(i,k))*lcldm(i,k)+ &
                      (nnuccri(i,k)+nmultrg(i,k))*precip_frac(i,k)+nnudep(i,k)*icldm(i,k))/ &
                      ((-nsubi(i,k)-niagg(i,k))*icldm(i,k))*omsm
              endif
!--trude               
!--ag
              nsubi(i,k) = nsubi(i,k)*ratio
              niagg(i,k) = niagg(i,k)*ratio
!--kt
           end if

        end do

     end if

!++ag Graupel Conservation Checks
!-------------------------------------------------------------------
     if(do_hail.or.do_graupel) then
! conservation of graupel mass
!-------------------------------------------------------------------
        do i=1,mgncol

           dum= ((-pracg(i,k)-pgracs(i,k)-prdg(i,k)-psacr(i,k)-mnuccr(i,k))*precip_frac(i,k) &
                + (-psacwg(i,k)-pgsacw(i,k))*lcldm(i,k))*deltat
           
           if (dum.gt.qg(i,k)) then

! hm added
! note: prdg is always negative (like prds), so it needs to be subtracted in ratio
              ratio = (qg(i,k)/deltat + (pracg(i,k)+pgracs(i,k)+psacr(i,k)+mnuccr(i,k))*precip_frac(i,k) &
                   + (psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)) / ((-prdg(i,k))*precip_frac(i,k))  *omsm

              prdg(i,k)= prdg(i,k)*ratio

           end if

        end do

! conservation of graupel number: not needed, no sinks
!-------------------------------------------------------------------
     end if
!--ag

     do i=1,mgncol

        ! next limit ice and snow sublimation and rain evaporation
        ! get estimate of q and t at end of time step
        ! don't include other microphysical processes since they haven't
        ! been limited via conservation checks yet

!++ag need to add graupel sublimation/evap here too (prdg)? May not need eprdg?
!++ag
!        if ((pre(i,k)+prds(i,k))*precip_frac(i,k)+ice_sublim(i,k) < -1.e-20_r8) then
!
!           qtmp=q(i,k)-(ice_sublim(i,k)+vap_dep(i,k)+mnuccd(i,k)+ &
!                (pre(i,k)+prds(i,k))*precip_frac(i,k))*deltat
!           ttmp=t(i,k)+((pre(i,k)*precip_frac(i,k))*xxlv+ &
!                (prds(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k))*xxls)*deltat/cpp
!++ktc removed prds
        if ((pre(i,k)+prdg(i,k))*precip_frac(i,k)+ice_sublim(i,k) < -1.e-20_r8) then
           qtmp=q(i,k)-(ice_sublim(i,k)+vap_dep(i,k)+mnuccd(i,k)+ &
                (pre(i,k)+prdg(i,k))*precip_frac(i,k))*deltat
           ttmp=t(i,k)+((pre(i,k)*precip_frac(i,k))*xxlv+ &
                (prdg(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k))*xxls)*deltat/cpp
!--ktc
!--ag
           ! use rhw to allow ice supersaturation
           call qsat_water(ttmp, p(i,k), esn, qvn)

           ! modify ice/precip evaporation rate if q > qsat
           if (qtmp > qvn) then

!++ag 
!              dum1=pre(i,k)*precip_frac(i,k)/((pre(i,k)+prds(i,k))*precip_frac(i,k)+ice_sublim(i,k))
!              dum2=prds(i,k)*precip_frac(i,k)/((pre(i,k)+prds(i,k))*precip_frac(i,k)+ice_sublim(i,k))
!++kt removed prds so dum2 = 0
              dum1=pre(i,k)*precip_frac(i,k)/((pre(i,k)+prdg(i,k))*precip_frac(i,k)+ice_sublim(i,k))
              dum2=0._r8
              dum3=prdg(i,k)*precip_frac(i,k)/((pre(i,k)+prdg(i,k))*precip_frac(i,k)+ice_sublim(i,k))
!--ag kt
              ! recalculate q and t after vap_dep and mnuccd but without evap or sublim
              qtmp=q(i,k)-(vap_dep(i,k)+mnuccd(i,k))*deltat
              ttmp=t(i,k)+((vap_dep(i,k)+mnuccd(i,k))*xxls)*deltat/cpp

              ! use rhw to allow ice supersaturation
              call qsat_water(ttmp, p(i,k), esn, qvn)

              dum=(qtmp-qvn)/(1._r8 + xxlv_squared*qvn/(cpp*rv*ttmp**2))
              dum=min(dum,0._r8)

              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              pre(i,k)=dum*dum1/deltat/precip_frac(i,k)

              ! do separately using RHI for prds and ice_sublim
              call qsat_ice(ttmp, p(i,k), esn, qvn)

              dum=(qtmp-qvn)/(1._r8 + xxls_squared*qvn/(cpp*rv*ttmp**2))
              dum=min(dum,0._r8)

!++ag
              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              prdg(i,k) = dum*dum3/deltat/precip_frac(i,k)
!--ag
              ! don't divide ice_sublim by cloud fraction since it is grid-averaged
!++ag
!              dum1 = (1._r8-dum1-dum2)
              dum1 = (1._r8-dum1-dum2-dum3)
!--ag
              ice_sublim(i,k) = dum*dum1/deltat
           end if
        end if

     end do

     ! Big "administration" loop enforces conservation, updates variables
     ! that accumulate over substeps, and sets output variables.

     do i=1,mgncol

        ! get tendencies due to microphysical conversion processes
        !==========================================================
        ! note: tendencies are multiplied by appropriate cloud/precip
        ! fraction to get grid-scale values
        ! note: vap_dep is already grid-average values

        ! The net tendencies need to be added to rather than overwritten,
        ! because they may have a value already set for instantaneous
        ! melting/freezing.

!++ag
!        qvlat(i,k) = qvlat(i,k)-(pre(i,k)+prds(i,k))*precip_frac(i,k)-&
!             vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*lcldm(i,k)
!++kt remove prds
! ++ trude. mnudep*icldm instead of mnudep*lcldm
        qvlat(i,k) = qvlat(i,k)-(pre(i,k)*precip_frac(i,k))-&
             vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*icldm(i,k) &
             -prdg(i,k)*precip_frac(i,k) 

!        tlat(i,k) = tlat(i,k)+((pre(i,k)*precip_frac(i,k)) &
!             *xxlv+(prds(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k)+mnudep(i,k)*lcldm(i,k))*xxls+ &
!             ((bergs(i,k)+psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k))*lcldm(i,k)+(mnuccr(i,k)+ &
!             pracs(i,k)+mnuccri(i,k))*precip_frac(i,k)+berg(i,k))*xlf)
!++kt remove bergs
! ++ trude. mnudep*icldm instead of mnudep*lcldm
         tlat(i,k) = tlat(i,k)+((pre(i,k)*precip_frac(i,k))*xxlv+ &
             ((prdg(i,k))*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+ &
                mnuccd(i,k)+mnudep(i,k)*icldm(i,k))*xxls+ &
             ((psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+psacwg(i,k)+ &
                qmultg(i,k)+pgsacw(i,k))*lcldm(i,k)+ &
             (mnuccr(i,k)+pracs(i,k)+mnuccri(i,k)+pracg(i,k)+pgracs(i,k)+qmultrg(i,k))*precip_frac(i,k)+ &
             berg(i,k))*xlf ) 

!       qctend(i,k) = qctend(i,k)+ &
!             (-pra(i,k)-prc(i,k)-mnuccc(i,k)-mnucct(i,k)-msacwi(i,k)- &
!             psacws(i,k)-bergs(i,k))*lcldm(i,k)-berg(i,k)
!++kt remove bergs
       qctend(i,k) = qctend(i,k)+ &
             (-pra(i,k)-prc(i,k)-mnuccc(i,k)-mnucct(i,k)-msacwi(i,k)- &
             psacws(i,k)-qmultg(i,k)-psacwg(i,k)-pgsacw(i,k))*lcldm(i,k)-berg(i,k)


        if (do_cldice) then
!           qitend(i,k) = qitend(i,k)+ &
!                (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k)+(-prci(i,k)- &
!                prai(i,k))*icldm(i,k)+vap_dep(i,k)+berg(i,k)+ice_sublim(i,k)+ &
!                mnuccd(i,k)+mnuccri(i,k)*precip_frac(i,k)
!++kt move psacws and pracs from snow to ice
           qitend(i,k) = qitend(i,k)+ &
                (mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+psacws(i,k)+ &
                qmultg(i,k))*lcldm(i,k)+ mnudep(i,k)*icldm(i,k)+&
                vap_dep(i,k)+berg(i,k)+ice_sublim(i,k)+ &
                mnuccd(i,k)+(mnuccri(i,k)+pracs(i,k)+qmultrg(i,k))*precip_frac(i,k)
        end if

!        qrtend(i,k) = qrtend(i,k)+ &
!             (pra(i,k)+prc(i,k))*lcldm(i,k)+(pre(i,k)-pracs(i,k)- &
!             mnuccr(i,k)-mnuccri(i,k))*precip_frac(i,k)

        qrtend(i,k) = qrtend(i,k)+ &
             (pra(i,k)+prc(i,k))*lcldm(i,k)+(pre(i,k)-pracs(i,k)- &
             mnuccr(i,k)-mnuccri(i,k)-qmultrg(i,k)-pracg(i,k)-pgracs(i,k))*precip_frac(i,k)

        if (do_hail.or.do_graupel) then
           qgtend(i,k) = qgtend(i,k) + (pracg(i,k)+pgracs(i,k)+prdg(i,k)+psacr(i,k)+mnuccr(i,k))*precip_frac(i,k) &
                + (psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)

           qitend(i,k) = qitend(i,k)-psacr(i,k)*precip_frac(i,k)
        else
           !necessary since mnuccr moved to graupel
           qitend(i,k) = qitend(i,k)+mnuccr(i,k)*precip_frac(i,k)  
        end if
!--ag        
!++kt
        cmeout(i,k) = cmeout(i,k) + vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

        ! add output for cmei (accumulate)
        cmeitot(i,k) = cmeitot(i,k) + vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)
!--kt
        ! assign variables for trop_mozart, these are grid-average
        !-------------------------------------------------------------------
        ! evaporation/sublimation is stored here as positive term

!++add evaporation/sublimation of graupel too? YES: After conservation checks. 

!++ag
!ADD GRAUPEL to evapsnow: prdg. (sign? same as prds: negative, so this is a positive number)
!        evapsnow(i,k) = -prds(i,k)*precip_frac(i,k)
!++kt  remove prds
        evapsnow(i,k) = -prdg(i,k)*precip_frac(i,k)
!--ag
        nevapr(i,k) = -pre(i,k)*precip_frac(i,k)
        prer_evap(i,k) = -pre(i,k)*precip_frac(i,k)

        ! change to make sure prain is positive: do not remove snow from
        ! prain used for wet deposition

!++AG NEED TO MAKE CONSISTENT WITH BUDGETS
        prain(i,k) = (pra(i,k)+prc(i,k))*lcldm(i,k)+(-pracs(i,k)- &
             mnuccr(i,k)-mnuccri(i,k))*precip_frac(i,k)
!++kt remove prai and bergs and prci from prodsnow
        if (do_hail .or. do_graupel) then
! Subtract PSACR here or not? Ask Hugh
           prodsnow(i,k) = (psacws(i,k))*lcldm(i,k)+(pracs(i,k))*precip_frac(i,k)
        else
           prodsnow(i,k) = (psacws(i,k))*lcldm(i,k)+(pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)
        end if

        ! following are used to calculate 1st order conversion rate of cloud water
        !    to rain and snow (1/s), for later use in aerosol wet removal routine
        ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
        !    used to calculate pra, prc, ... in this routine
        ! qcsinksum_rate1ord = { rate of direct transfer of cloud water to rain & snow }
        !                      (no cloud ice or bergeron terms)

!++AG NEED TO MAKE CONSITANT: PGSACW, PSACWG (check budgets)? More sink terms? Check. No. Just loss to precip. 
!Ask Hugh
!        qcsinksum_rate1ord(i,k) = (pra(i,k)+prc(i,k)+psacws(i,k))*lcldm(i,k) 
        qcsinksum_rate1ord(i,k) = (pra(i,k)+prc(i,k)+psacws(i,k)+psacwg(i,k)+pgsacw(i,k))*lcldm(i,k) 
!--ag
        ! Avoid zero/near-zero division.
        qcsinksum_rate1ord(i,k) = qcsinksum_rate1ord(i,k) / &
             max(qc(i,k),1.0e-30_r8)


        ! microphysics output, note this is grid-averaged
        pratot(i,k) = pra(i,k)*lcldm(i,k)
        prctot(i,k) = prc(i,k)*lcldm(i,k)
        mnuccctot(i,k) = mnuccc(i,k)*lcldm(i,k)
        mnuccttot(i,k) = mnucct(i,k)*lcldm(i,k)
        msacwitot(i,k) = msacwi(i,k)*lcldm(i,k)
        psacwstot(i,k) = psacws(i,k)*lcldm(i,k)
!++kt
        npsacwstot(i,k) = npsacws(i,k)*lcldm(i,k)
        bergtot(i,k) = berg(i,k)
        mnuccdtot(i,k) = mnuccd(i,k)*icldm(i,k)
!--kt
        pracstot(i,k) = pracs(i,k)*precip_frac(i,k)
        mnuccrtot(i,k) = mnuccr(i,k)*precip_frac(i,k)
!++ag
        mnuccritot(i,k) = mnuccri(i,k)*precip_frac(i,k)
!--ag

!++ag Hail/Graupel tendencies for output
        psacrtot(i,k) = psacr(i,k)*precip_frac(i,k)
        pracgtot(i,k) = pracg(i,k)*precip_frac(i,k)
        psacwgtot(i,k) = psacwg(i,k)*lcldm(i,k)
        pgsacwtot(i,k) = pgsacw(i,k)*lcldm(i,k)
        pgracstot(i,k) = pgracs(i,k)*precip_frac(i,k)
        prdgtot(i,k) = prdg(i,k)*precip_frac(i,k)
        qmultgtot(i,k) = qmultg(i,k)*lcldm(i,k)
        qmultrgtot(i,k) = qmultrg(i,k)*precip_frac(i,k)
        npracgtot(i,k) = npracg(i,k)*precip_frac(i,k) 
        nscngtot(i,k) = nscng(i,k)*lcldm(i,k) 
        ngracstot(i,k) = ngracs(i,k)*precip_frac(i,k)
        nmultgtot(i,k) = nmultg(i,k)*lcldm(i,k)
        nmultrgtot(i,k) = nmultrg(i,k)*precip_frac(i,k)
        npsacwgtot(i,k) = npsacwg(i,k)*lcldm(i,k)
!--ag

!++ag
!        nctend(i,k) = nctend(i,k)+&
!             (-nnuccc(i,k)-nnucct(i,k)-npsacws(i,k)+nsubc(i,k) &
!             -npra(i,k)-nprc1(i,k))*lcldm(i,k)

        nctend(i,k) = nctend(i,k)+&
             (-nnuccc(i,k)-nnucct(i,k)-npsacws(i,k)+nsubc(i,k) &
             -npra(i,k)-nprc1(i,k)-npsacwg(i,k))*lcldm(i,k)


        if (do_cldice) then
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = 0._r8
           end if
!           nitend(i,k) = nitend(i,k)+ nnuccd(i,k)+ &
!                (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k)+(nsubi(i,k)-nprci(i,k)- &
!                nprai(i,k))*icldm(i,k)+nnuccri(i,k)*precip_frac(i,k)

          nitend(i,k) = nitend(i,k)+ nnuccd(i,k)+ &
                (nnucct(i,k)+tmpfrz+nsacwi(i,k)+nmultg(i,k))*lcldm(i,k)&
                +(nsubi(i,k)+niagg(i,k)+nnudep(i,k))*icldm(i,k) &
                +(nnuccri(i,k)+nmultrg(i,k))*precip_frac(i,k)

        end if

!++kt removed nstend
!        if(do_graupel.or.do_hail) then
!           ngtend(i,k) = ngtend(i,k)+nscng(i,k)*lcldm(i,k)+(ngracs(i,k)+nnuccr(i,k))*precip_frac(i,k)
!     end if

        if(do_graupel.or.do_hail) then
!           nstend(i,k) = nstend(i,k)+(nsubs(i,k)+ &
!                nsagg(i,k)+nnuccr(i,k))*precip_frac(i,k)+nprci(i,k)*icldm(i,k)
!++trude nsubi(i,k) is alredy in nitend(i,k) equtaiton above. Remove the nsubi process below
           nitend(i,k) = nitend(i,k)+(-ngracs(i,k))*precip_frac(i,k)-nscng(i,k)*lcldm(i,k)
!--trude

           ngtend(i,k) = ngtend(i,k)+nscng(i,k)*lcldm(i,k)+(ngracs(i,k)+nnuccr(i,k))*precip_frac(i,k)

        else
           !necessary since mnuccr moved to graupel
!++trude, remove nsubi, alreadi included above.
           nitend(i,k) = nitend(i,k)+(nnuccr(i,k))*precip_frac(i,k)

        end if

!        nrtend(i,k) = nrtend(i,k)+ &
!             nprc(i,k)*lcldm(i,k)+(nsubr(i,k)-npracs(i,k)-nnuccr(i,k) &
!             -nnuccri(i,k)+nragg(i,k))*precip_frac(i,k)

        nrtend(i,k) = nrtend(i,k)+ &
             nprc(i,k)*lcldm(i,k)+(nsubr(i,k)-npracs(i,k)-nnuccr(i,k) &
             -nnuccri(i,k)+nragg(i,k)-npracg(i,k)-ngracs(i,k))*precip_frac(i,k)
!--ag

        ! make sure that ni at advanced time step does not exceed
        ! maximum (existing N + source terms*dt), which is possible if mtime < deltat
        ! note that currently mtime = deltat
        !================================================================

        if (do_cldice .and. nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax(i,k)) then
           nitend(i,k)=max(0._r8,(nimax(i,k)-ni(i,k))/deltat)
        end if

     end do

     ! End of "administration" loop
      end if ! mgncol > 0
  end do micro_vert_loop ! end k loop

  !-----------------------------------------------------
  ! convert rain/snow q and N for output to history, note,
  ! output is for gridbox average

  qrout = qr
  nrout = nr * rho
!++kt remove qsout
!++ag
  qgout = qg
  ngout = ng * rho
!--ag


  ! calculate n0r and lamr from rain mass and number
  ! divide by precip fraction to get in-precip (local) values of
  ! rain mass and number, divide by rhow to get rain number in kg^-1

  do k=1,nlev
      if( mgncol > 0.0_r8) then
         call size_dist_param_basic(mg_rain_props, qric(:,k), nric(:,k), lamr(:,k), mgncol, n0=n0r(:,k))

     !        Calculate rercld

     !        calculate mean size of combined rain and cloud water

         call calc_rercld(lamr(:,k), n0r(:,k), lamc(:,k), pgam(:,k), qric(:,k), qcic(:,k), ncic(:,k), &
         rercld(:,k), mgncol)
      end if
  enddo

  ! Assign variables back to start-of-timestep values
  ! Some state variables are changed before the main microphysics loop
  ! to make "instantaneous" adjustments. Afterward, we must move those changes
  ! back into the tendencies.
  ! These processes:
  !  - Droplet activation (npccn, impacts nc)
  !  - Instantaneous snow melting  (minstsm/ninstsm, impacts qr/qs/nr/ns)
  !  - Instantaneous rain freezing (minstfr/ninstrf, impacts qr/qs/nr/ns)
  !================================================================================

  ! Re-apply droplet activation tendency
  nc = ncn
  nctend = nctend + npccn


  !.............................................................................

  !================================================================================

  ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
  nevapr = nevapr + evapsnow
  prain = prain + prodsnow

  do k=1,nlev

     do i=1,mgncol

        ! calculate sedimentation for cloud water and ice
!++ag   ! and Graupel (mg3)
        !================================================================================

        ! update in-cloud cloud mixing ratio and number concentration
        ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
        ! note: these are in-cloud values***, hence we divide by cloud fraction

        dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)
        dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)/icldm(i,k)
        dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k),0._r8)
        dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)/icldm(i,k),0._r8)

        dumr(i,k) = (qr(i,k)+qrtend(i,k)*deltat)/precip_frac(i,k)
        dumnr(i,k) = max((nr(i,k)+nrtend(i,k)*deltat)/precip_frac(i,k),0._r8)
!++kt remove dums and dumns

!++ag Add graupel
        dumg(i,k) = (qg(i,k)+qgtend(i,k)*deltat)/precip_frac(i,k)
        dumng(i,k) = max((ng(i,k)+ngtend(i,k)*deltat)/precip_frac(i,k),0._r8)
        ! switch for specification of droplet and crystal number
        if (ngcons) then
           dumng(i,k)=ngnst/rho(i,k)
        end if
!--ag
        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)
        end if
     enddo
  enddo

  do k=1,nlev

     ! obtain new slope parameter to avoid possible singularity
!++kt
      do i=1,mgncol

         !if ((qiic(i,k).ge.qsmall).and.(dumni(i,k).lt.icsmall)) then
         !   write(*,*) "MG4 WARNING 2 ice present in cloud with no number conc at i=",i," k=",k
         !end if

         if ((qiic(i,k).ge.qsmall).and.(dumi(i,k).ge.qsmall).and.(dumni(i,k).ge.icsmall)) then
!     call size_dist_param_basic(mg_ice_props,dumi(:,k), dumni(:,k), &
!     lami(:,k), mgncol)
! find index for qi (total ice mass mixing ratio)
            dum11 = (dlog10(dumi(i,k))+16._r8)*1.41328_r8
            dumit = int(dum11)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum11 = min(dum11,real(isize))
            dum11 = min(dum11,dble(isize))
            dum11 = max(dum11,1._r8)
            dumit = max(1,dumit)
            dumit = min(isize-1,dumit)

           ! find index for Ni (ice number mixing ratio)
            dum22 = (dlog10(dumni(i,k))+10._r8)*1.10731_r8
            dumk = int(dum22)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum22 = min(dum22,real(jsize))
            dum22 = min(dum22,dble(jsize))
            dum22 = max(dum22,1._r8)
            dumk = max(1,dumk)
            dumk = min(jsize-1,dumk)

           ! find index for scaled mean rain size
           ! if no rain, then just choose dumj = 1 and don't calculate rain-ice
           ! collection processes

           ! find index for temperature
            dum4  = t(i,k)      ! trude, and something else here for finding index qitot(i,k,iice)*3. + 1.
            if (t(i,k) .lt. 213.15_r8) then ! (-60)
               dum4 = 3         !trude, need to calculate the "interger" value
            else if (t(i,k).ge. 243.15_r8) then  
               dum4 = 1         ! (-30) 
            else
               dum4 = 17.21_r8-(t(i,k)/15._r8) 
            endif
            dumtt = int(dum4)
           ! set limits
!             dum4  = min(dum4,real(tsize))
            dum4  = min(dum4,dble(tsize))
            dum4  = max(dum4,1._r8)
            dumtt = max(1,dumtt)
            dumtt = min(tsize-1,dumtt)

           ! call to lookup table interpolation subroutines to get process rates

            call access_lookup_table(dumtt,dumit,dumk,7,dum11,dum22,dum4,f1pr9)
            call access_lookup_table(dumtt,dumit,dumk,8,dum11,dum22,dum4,f1pr10)

            dumni(i,k) = min(dumni(i,k),f1pr9)
            dumni(i,k) = max(dumni(i,k),f1pr10)
!            dumni(i,k) = min(dumni(i,k), 500.e3_r8/rho(i,k))
            dumni(i,k) = min(dumni(i,k), micro_mg_max_nicons/rho(i,k))
         end if                 !qiic > qsmall
        end do !mgncol
!--kt
        if (mgncol > 0) then
           call size_dist_param_liq(mg_liq_props, dumc(:,k), dumnc(:,k), rho(:,k), &
           pgam(:,k), lamc(:,k), mgncol)
        end if

  enddo !k=nlev

  do k=1,nlev
     do i=1,mgncol

        ! calculate number and mass weighted fall velocity for droplets and cloud ice
        !-------------------------------------------------------------------


        if (dumc(i,k).ge.qsmall) then

           vtrmc(i,k)=acn(i,k)*gamma(4._r8+bc+pgam(i,k))/ &
                (lamc(i,k)**bc*gamma(pgam(i,k)+4._r8))

           fc(i,k) = g*rho(i,k)*vtrmc(i,k)

           fnc(i,k) = g*rho(i,k)* &
                acn(i,k)*gamma(1._r8+bc+pgam(i,k))/ &
                (lamc(i,k)**bc*gamma(pgam(i,k)+1._r8))
        else
           fc(i,k) = 0._r8
           fnc(i,k)= 0._r8
        end if

        ! calculate number and mass weighted fall velocity for cloud ice

        !if ((dumi(i,k).ge.qsmall).and.(dumni(i,k).lt.icsmall)) then
        !   write(*,*) "MG4 WARNING 3 ice present in cloud with no number conc at i=",i," k=",k
        !end if

        if ((dumi(i,k).ge.qsmall).and.(dumni(i,k).ge.icsmall)) then
!++kt Since ice has been updated, call access_lookup_table again for fall speed. 
           
           ! find index for qi (total ice mass mixing ratio)
           dum11 = (dlog10(dumi(i,k))+16._r8)*1.41328_r8
           dumit = int(dum11)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum11 = min(dum11,real(isize))
           dum11 = min(dum11,dble(isize))
           dum11 = max(dum11,1._r8)
           dumit = max(1,dumit)
           dumit = min(isize-1,dumit)

           ! find index for Ni (ice number mixing ratio)
           dum22 = (dlog10(dumni(i,k))+10._r8)*1.10731_r8
           dumk = int(dum22)

           ! set limits to make sure the calculated index doesn't exceed range of lookup table
           !             dum22 = min(dum22,real(jsize))
           dum22 = min(dum22,dble(jsize))
           dum22 = max(dum22,1._r8)
           dumk = max(1,dumk)
           dumk = min(jsize-1,dumk)
           ! find index for temperature
           dum4  = t(i,k)   ! trude, and something else here for finding index qitot(i,k,iice)*3. + 1.
           if (t(i,k) .lt. 213.15_r8) then    ! (-60)
              dum4 = 3     !trude, need to calculate the "interger" value
           else if (t(i,k).ge. 243.15_r8) then  
              dum4 = 1                      ! (-30) 
           else
              dum4 = 17.21_r8-(t(i,k)/15._r8) 
           endif
           dumtt = int(dum4)
           ! set limits
           dum4  = min(dum4,dble(tsize))
           dum4  = max(dum4,1._r8)
           dumtt = max(1,dumtt)
           dumtt = min(tsize-1,dumtt)

           call access_lookup_table(dumtt,dumit,dumk,1,dum11,dum22,dum4,f1pr1)
           call access_lookup_table(dumtt,dumit,dumk,2,dum11,dum22,dum4,f1pr2)

           ! ++ trude, to test changes to fall speed
           f1pr2=f1pr2*micro_mg_vtrmi_factor
           ! -- trude

           vtrmi(i,k) = min(f1pr2*rhof(i,k),1.2_r8*rhof(i,k))
           fi(i,k) = g*rho(i,k)*vtrmi(i,k)

           ! fni(i,k) = g*rho(i,k)* &
           !     min(ain(i,k)*gamma_bi_plus1/lami(i,k)**bi,1.2_r8*rhof(i,k))
           fni(i,k) = g*rho(i,k)* &
                min(f1pr1*rhof(i,k),1.2_r8*rhof(i,k))

           ! adjust the ice fall velocity for smaller (r < 20 um) ice
           ! particles (blend over 8-20 um)
           ! kt: remove this for now...
           ! irad = 1.5_r8 / lami(i,k) * 1e6_r8
           ! ifrac = min(1._r8, max(0._r8, (irad - 18._r8) / 2._r8))
 
           ! if (ifrac .lt. 1._r8) then
           !    vtrmi(i,k) = ifrac * vtrmi(i,k) + & 
           !       (1._r8 - ifrac) * &
           !       min(ajn(i,k)*gamma_bj_plus4/(6._r8*lami(i,k)**bj), &
           !       1.2_r8*rhof(i,k))

           !    fi(i,k) = g*rho(i,k)*vtrmi(i,k)
           !    fni(i,k) = ifrac * fni(i,k) + & 
           !       (1._r8 - ifrac) * &
           !       g*rho(i,k)* &
           !       min(ajn(i,k)*gamma_bj_plus1/lami(i,k)**bj,1.2_r8*rhof(i,k))
           ! end if
           !--
!--kt
        else
           fi(i,k) = 0._r8
           fni(i,k)= 0._r8
        end if

     enddo

  enddo

  do k=1,nlev

        ! fallspeed for rain
      if(mgncol > 0) then
         call size_dist_param_basic(mg_rain_props, dumr(:,k), dumnr(:,k), &
         lamr(:,k), mgncol)
      end if
  enddo

  do k=1,nlev

     do i=1,mgncol
        if (lamr(i,k).ge.qsmall) then

           ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

           unr(i,k) = min(arn(i,k)*gamma_br_plus1/lamr(i,k)**br,9.1_r8*rhof(i,k))
           umr(i,k) = min(arn(i,k)*gamma_br_plus4/(6._r8*lamr(i,k)**br),9.1_r8*rhof(i,k))

           fr(i,k) = g*rho(i,k)*umr(i,k)
           fnr(i,k) = g*rho(i,k)*unr(i,k)

        else
           fr(i,k)=0._r8
           fnr(i,k)=0._r8
        end if

!++kt remove fallspeed for snow

!++ag
        ! fallspeed for graupel

        call size_dist_param_basic(mg_graupel_props, dumg(i,k), dumng(i,k), &
             lamg(i,k))

        if (lamg(i,k).ge.qsmall) then

           ! 'final' values of number and mass weighted mean fallspeed for graupel (m/s)
           umg(i,k) = min(agn(i,k)*gamma(4._r8+bgtmp)/(6._r8*lamg(i,k)**bgtmp),20._r8*rhof(i,k))
           ung(i,k) = min(agn(i,k)*gamma(1._r8+bgtmp)/lamg(i,k)**bgtmp,20._r8*rhof(i,k))

           fg(i,k) = g*rho(i,k)*umg(i,k)
           fng(i,k) = g*rho(i,k)*ung(i,k)

        else
           fg(i,k)=0._r8
           fng(i,k)=0._r8
        end if

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation

        dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
        dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat),0._r8)
        dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
        dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat),0._r8)
        dumr(i,k) = (qr(i,k)+qrtend(i,k)*deltat)
        dumnr(i,k) = max((nr(i,k)+nrtend(i,k)*deltat),0._r8)
!++kt remove dums
!++ag
        dumg(i,k) = (qg(i,k)+qgtend(i,k)*deltat)
        dumng(i,k) = max((ng(i,k)+ngtend(i,k)*deltat),0._r8)
!--ag

        if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
        if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
        if (dumr(i,k).lt.qsmall) dumnr(i,k)=0._r8
!++ag
        if (dumg(i,k).lt.qsmall) dumng(i,k)=0._r8
!--ag

     enddo
  end do       !!! vertical loop

  do k=1,nlev
     do i=1,mgncol
       pdel_inv(i,k) = 1._r8/pdel(i,k)
     enddo
  enddo

  ! initialize nstep for sedimentation sub-steps

  ! calculate number of split time steps to ensure courant stability criteria
  ! for sedimentation calculations
  !-------------------------------------------------------------------
  do i=1,mgncol
     nstep = 1 + int(max( &
          maxval( fi(i,:)*pdel_inv(i,:)), &
          maxval(fni(i,:)*pdel_inv(i,:))) &
          * deltat)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        if (do_cldice) then
           falouti  = fi(i,:)  * dumi(i,:)
           faloutni = fni(i,:) * dumni(i,:)
        else
           falouti  = 0._r8
           faloutni = 0._r8
        end if

        ! top of model

        k = 1

        ! add fallout terms to microphysical tendencies
        faltndi = falouti(k)/pdel(i,k)
        faltndni = faloutni(k)/pdel(i,k)
        qitend(i,k) = qitend(i,k)-faltndi/nstep
        nitend(i,k) = nitend(i,k)-faltndni/nstep

        ! sedimentation tendency for output
        qisedten(i,k)=qisedten(i,k)-faltndi/nstep

        dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
        dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep

        do k = 2,nlev

           ! for cloud liquid and ice, if cloud fraction increases with height
           ! then add flux from above to both vapor and cloud water of current level
           ! this means that flux entering clear portion of cell from above evaporates
           ! instantly

           ! note: this is not an issue with precip, since we assume max overlap
           dum1=icldm(i,k)/icldm(i,k-1)
           dum1=min(dum1,1._r8)

           faltndqie=(falouti(k)-falouti(k-1))/pdel(i,k)
!++ trude. Why is faltndi with dum1 commented out?? Why is dum1 not included? In MG3 they are included
           !faltndi=(falouti(k)-dum1*falouti(k-1))/pdel(i,k)
           !faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel(i,k)
           faltndi=(falouti(k)-falouti(k-1))/pdel(i,k)
           faltndni=(faloutni(k)-faloutni(k-1))/pdel(i,k)
!-- trude

           ! add fallout terms to eulerian tendencies

           qitend(i,k) = qitend(i,k)-faltndi/nstep
           nitend(i,k) = nitend(i,k)-faltndni/nstep

           ! sedimentation tendency for output
           qisedten(i,k)=qisedten(i,k)-faltndi/nstep

           ! add terms to to evap/sub of cloud water

!++   ktc
          qvlat(i,k)=qvlat(i,k)-(faltndqie-faltndi)/nstep
!     for output
!++trude test. Uncomment qisevap 01/10/20
!           qisevap(i,k)=0._r8
           qisevap(i,k)=qisevap(i,k)-(faltndqie-faltndi)/nstep
!++   ktc
           !ask trude why??
           tlat(i,k)=tlat(i,k)+(faltndqie-faltndi)*xxls/nstep
!-- ktc
           dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
           dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep

        end do

        ! Ice flux
        do k = 1,nlev
          iflx(i,k+1) = iflx(i,k+1) + falouti(k) / g / real(nstep)
        end do

        ! units below are m/s
        ! sedimentation flux at surface is added to precip flux at surface
        ! to get total precip (cloud + precip water) rate

        prect(i) = prect(i)+falouti(nlev)/g/real(nstep)/1000._r8
        preci(i) = preci(i)+falouti(nlev)/g/real(nstep)/1000._r8

     end do

     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + int(max( &
          maxval( fc(i,:)*pdel_inv(i,:)), &
          maxval(fnc(i,:)*pdel_inv(i,:))) &
          * deltat)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        faloutc  = fc(i,:)  * dumc(i,:)
        faloutnc = fnc(i,:) * dumnc(i,:)

        ! top of model
        k = 1

        ! add fallout terms to microphysical tendencies
        faltndc = faloutc(k)/pdel(i,k)
        faltndnc = faloutnc(k)/pdel(i,k)
        qctend(i,k) = qctend(i,k)-faltndc/nstep
        nctend(i,k) = nctend(i,k)-faltndnc/nstep

        ! sedimentation tendency for output
        qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep

        dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
        dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

        do k = 2,nlev

           dum=lcldm(i,k)/lcldm(i,k-1)
           dum=min(dum,1._r8)
           faltndqce=(faloutc(k)-faloutc(k-1))/pdel(i,k)
           faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(i,k)
           faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(i,k)

           ! add fallout terms to eulerian tendencies
           qctend(i,k) = qctend(i,k)-faltndc/nstep
           nctend(i,k) = nctend(i,k)-faltndnc/nstep

           ! sedimentation tendency for output
           qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep

           ! add terms to to evap/sub of cloud water
           qvlat(i,k)=qvlat(i,k)-(faltndqce-faltndc)/nstep
           ! for output
           qcsevap(i,k)=qcsevap(i,k)-(faltndqce-faltndc)/nstep

           tlat(i,k)=tlat(i,k)+(faltndqce-faltndc)*xxlv/nstep

           dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
           dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

        end do

        !Liquid condensate flux here
        do k = 1,nlev
           lflx(i,k+1) = lflx(i,k+1) + faloutc(k) / g / real(nstep)
        end do

        prect(i) = prect(i)+faloutc(nlev)/g/real(nstep)/1000._r8

     end do

     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + int(max( &
          maxval( fr(i,:)*pdel_inv(i,:)), &
          maxval(fnr(i,:)*pdel_inv(i,:))) &
          * deltat)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        faloutr  = fr(i,:)  * dumr(i,:)
        faloutnr = fnr(i,:) * dumnr(i,:)

        ! top of model
        k = 1

        ! add fallout terms to microphysical tendencies
        faltndr = faloutr(k)/pdel(i,k)
        faltndnr = faloutnr(k)/pdel(i,k)
        qrtend(i,k) = qrtend(i,k)-faltndr/nstep
        nrtend(i,k) = nrtend(i,k)-faltndnr/nstep

        ! sedimentation tendency for output
        qrsedten(i,k)=qrsedten(i,k)-faltndr/nstep

        dumr(i,k) = dumr(i,k)-faltndr*deltat/real(nstep)
        dumnr(i,k) = dumnr(i,k)-faltndnr*deltat/real(nstep)

        do k = 2,nlev

           faltndr=(faloutr(k)-faloutr(k-1))/pdel(i,k)
           faltndnr=(faloutnr(k)-faloutnr(k-1))/pdel(i,k)

           ! add fallout terms to eulerian tendencies
           qrtend(i,k) = qrtend(i,k)-faltndr/nstep
           nrtend(i,k) = nrtend(i,k)-faltndnr/nstep

           ! sedimentation tendency for output
           qrsedten(i,k)=qrsedten(i,k)-faltndr/nstep

           dumr(i,k) = dumr(i,k)-faltndr*deltat/real(nstep)
           dumnr(i,k) = dumnr(i,k)-faltndnr*deltat/real(nstep)

        end do

        ! Rain Flux
        do k = 1,nlev
           rflx(i,k+1) = rflx(i,k+1) + faloutr(k) / g / real(nstep)
        end do

        prect(i) = prect(i)+faloutr(nlev)/g/real(nstep)/1000._r8

     end do

!++kt remove falouts and faloutns
     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     !nstep = 1 + int(max( &
     !     maxval( fs(i,:)*pdel_inv(i,:)), &
     !     maxval(fns(i,:)*pdel_inv(i,:))) &
     !     * deltat)



!++ag Graupel Sedimentation
    ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + int(max( &
          maxval( fg(i,:)*pdel_inv(i,:)), &
          maxval(fng(i,:)*pdel_inv(i,:))) &
          * deltat)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        faloutg  = fg(i,:)  * dumg(i,:)
        faloutng = fng(i,:) * dumng(i,:)

        ! top of model
        k = 1

        ! add fallout terms to microphysical tendencies
        faltndg = faloutg(k)/pdel(i,k)
        faltndng = faloutng(k)/pdel(i,k)
        qgtend(i,k) = qgtend(i,k)-faltndg/nstep
        ngtend(i,k) = ngtend(i,k)-faltndng/nstep

        ! sedimentation tendency for output
        qgsedten(i,k)=qgsedten(i,k)-faltndg/nstep

           dumg(i,k) = dumg(i,k)-faltndg*deltat/real(nstep)
           dumng(i,k) = dumng(i,k)-faltndng*deltat/real(nstep)

        do k = 2,nlev

           faltndg=(faloutg(k)-faloutg(k-1))/pdel(i,k)
           faltndng=(faloutng(k)-faloutng(k-1))/pdel(i,k)

           ! add fallout terms to eulerian tendencies
           qgtend(i,k) = qgtend(i,k)-faltndg/nstep
           ngtend(i,k) = ngtend(i,k)-faltndng/nstep

           ! sedimentation tendency for output
           qgsedten(i,k)=qgsedten(i,k)-faltndg/nstep

           dumg(i,k) = dumg(i,k)-faltndg*deltat/real(nstep)
           dumng(i,k) = dumng(i,k)-faltndng*deltat/real(nstep)

        end do   !! k loop

        ! Graupel Flux
        do k = 1,nlev
           gflx(i,k+1) = gflx(i,k+1) + faloutg(k) / g / real(nstep)
        end do

        ! Add to snow flux at surface
        prect(i) = prect(i)+faloutg(nlev)/g/real(nstep)/1000._r8
        preci(i) = preci(i)+faloutg(nlev)/g/real(nstep)/1000._r8

     end do   !! nstep loop
!--ag
  enddo
  ! end sedimentation

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! get new update for variables that includes sedimentation tendency
  ! note : here dum variables are grid-average, NOT in-cloud

  do k=1,nlev
     do i=1,mgncol
        dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)
        dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)

        dumr(i,k) = max(qr(i,k)+qrtend(i,k)*deltat,0._r8)
        dumnr(i,k) = max(nr(i,k)+nrtend(i,k)*deltat,0._r8)
!++kt remove dums
!++ag
        dumg(i,k) = max(qg(i,k)+qgtend(i,k)*deltat,0._r8)
        dumng(i,k) = max(ng(i,k)+ngtend(i,k)*deltat,0._r8)
!--ag

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)*lcldm(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)*icldm(i,k)
        end if

!++ag
        ! switch for specification of graupel number
        if (ngcons) then
           dumng(i,k)=ngnst/rho(i,k)*precip_frac(i,k)
        end if
!--ag

        if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
        if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
        if (dumr(i,k).lt.qsmall) dumnr(i,k)=0._r8
!--kt
!++ag
        if (dumg(i,k).lt.qsmall) dumng(i,k)=0._r8
!--ag
     enddo

  enddo

  ! calculate instantaneous processes (melting, homogeneous freezing)
  !====================================================================

  ! melting of ice at +2 C
  do k=1,nlev

     do i=1,mgncol
!++ kt melting of ice instead of snow
        if (t(i,k)+tlat(i,k)/cpp*deltat > snowmelt) then
           if (dumi(i,k) > 0._r8) then

              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -xlf/cpp*dumi(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt. snowmelt) then
                 dum = (t(i,k)+tlat(i,k)/cpp*deltat-snowmelt)*cpp/xlf
                 dum = dum/dumi(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qitend(i,k)=qitend(i,k)-dum*dumi(i,k)/deltat
              nitend(i,k)=nitend(i,k)-dum*dumni(i,k)/deltat
              qrtend(i,k)=qrtend(i,k)+dum*dumi(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)+dum*dumni(i,k)/deltat

              icemelttend(i,k) = -1.0_r8*dum*dumi(i,k)/deltat

              dum1=-xlf*dum*dumi(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              meltsdttot(i,k)=meltsdttot(i,k) + dum1
           end if
        end if
!--kt
     enddo
  enddo
!++ag

  ! melting of graupel at +2 C

  do k=1,nlev

     do i=1,mgncol

        if (t(i,k)+tlat(i,k)/cpp*deltat > snowmelt) then
           if (dumg(i,k) > 0._r8) then

              ! make sure melting graupel doesn't reduce temperature below threshold
              dum = -xlf/cpp*dumg(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt. snowmelt) then
                 dum = (t(i,k)+tlat(i,k)/cpp*deltat-snowmelt)*cpp/xlf
                 dum = dum/dumg(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if


              qgtend(i,k)=qgtend(i,k)-dum*dumg(i,k)/deltat
              ngtend(i,k)=ngtend(i,k)-dum*dumng(i,k)/deltat
              qrtend(i,k)=qrtend(i,k)+dum*dumg(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)+dum*dumng(i,k)/deltat

              dum1=-xlf*dum*dumg(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              meltsdttot(i,k)=meltsdttot(i,k) + dum1
           end if
        end if
     enddo
  enddo
!--ag



   do k=1,nlev
      do i=1,mgncol

        ! freezing of rain at -5 C

        if (t(i,k)+tlat(i,k)/cpp*deltat < rainfrze) then

           if (dumr(i,k) > 0._r8) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*dumr(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.rainfrze) then
                 dum = -(t(i,k)+tlat(i,k)/cpp*deltat-rainfrze)*cpp/xlf
                 dum = dum/dumr(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qrtend(i,k)=qrtend(i,k)-dum*dumr(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)-dum*dumnr(i,k)/deltat

              ! get mean size of rain = 1/lamr, add frozen rain to either snow or cloud ice
              ! depending on mean rain size
              ! add to graupel if using that option....

              call size_dist_param_basic(mg_rain_props, dumr(i,k), dumnr(i,k), &
                   lamr(i,k))

              if (lamr(i,k) < 1._r8/Dcs) then
!++ag freeze rain to graupel
!++kt always send to ice
                 if(do_hail.or.do_graupel) then
                    qgtend(i,k)=qgtend(i,k)+dum*dumr(i,k)/deltat
                    ngtend(i,k)=ngtend(i,k)+dum*dumnr(i,k)/deltat
                 else
                    qitend(i,k)=qitend(i,k)+dum*dumr(i,k)/deltat
                    nitend(i,k)=nitend(i,k)+dum*dumnr(i,k)/deltat
                 end if
!--ag kt
              else
                 qitend(i,k)=qitend(i,k)+dum*dumr(i,k)/deltat
                 nitend(i,k)=nitend(i,k)+dum*dumnr(i,k)/deltat
              end if

              rainfrztend(i,k)=dum*dumr(i,k)/deltat
              ! heating tendency
              dum1 = xlf*dum*dumr(i,k)/deltat
              frzrdttot(i,k)=frzrdttot(i,k) + dum1
              tlat(i,k)=tlat(i,k)+dum1

           end if
        end if

      enddo
   enddo

   if (do_cldice) then
      !-- kt removed

     ! homogeneously freeze droplets at -40 C
     !-----------------------------------------------------------------

     do k=1,nlev
        do i=1,mgncol
           if (t(i,k)+tlat(i,k)/cpp*deltat < 233.15_r8) then
              if (dumc(i,k) > 0._r8) then

                 ! limit so that freezing does not push temperature above threshold
                 dum = dumc(i,k)*xlf/cpp
                 if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.233.15_r8) then
                    dum = -(t(i,k)+tlat(i,k)/cpp*deltat-233.15_r8)*cpp/xlf
                    dum = dum/dumc(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qitend(i,k)=qitend(i,k)+dum*dumc(i,k)/deltat
                 ! for output
                 homotot(i,k)=dum*dumc(i,k)/deltat

!++ kt reduce these values by half (deltat*8)
                 ! assume 25 micron mean volume radius of homogeneously frozen droplets
                 ! consistent with size of detrained ice in stratiform.F90
!                 nitend(i,k)=nitend(i,k)+dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*1.563e-14_r8* &
 ! ++ trude
                 nitend(i,k)=nitend(i,k)+dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*micro_mg_homog_size**3._r8* &
                 !-- trude
                 500._r8)/(deltat*8)
!--kt
                 qctend(i,k)=((1._r8-dum)*dumc(i,k)-qc(i,k))/deltat
                 nctend(i,k)=((1._r8-dum)*dumnc(i,k)-nc(i,k))/deltat
                 tlat(i,k)=tlat(i,k)+xlf*dum*dumc(i,k)/deltat
              end if
           end if
        enddo 
     enddo 

     ! remove any excess over-saturation, which is possible due to non-linearity when adding
     ! together all microphysical processes
     !-----------------------------------------------------------------
     ! follow code similar to old CAM scheme
     do k=1,nlev
        do i=1,mgncol

           qtmp=q(i,k)+qvlat(i,k)*deltat
           ttmp=t(i,k)+tlat(i,k)/cpp*deltat

           ! use rhw to allow ice supersaturation
           call qsat_water(ttmp, p(i,k), esn, qvn)

           if (qtmp > qvn .and. qvn > 0 .and. allow_sed_supersat) then
              ! expression below is approximate since there may be ice deposition
              dum = (qtmp-qvn)/(1._r8+xxlv_squared*qvn/(cpp*rv*ttmp**2))/deltat
              ! add to output cme
              cmeout(i,k) = cmeout(i,k)+dum
              ! now add to tendencies, partition between liquid and ice based on temperature
              if (ttmp > 268.15_r8) then
                 dum1=0.0_r8
                 ! now add to tendencies, partition between liquid and ice based on te
                 !-------------------------------------------------------
              else if (ttmp < 238.15_r8) then
                 dum1=1.0_r8
              else
                 dum1=(268.15_r8-ttmp)/30._r8
              end if

              dum = (qtmp-qvn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                   *qvn/(cpp*rv*ttmp**2))/deltat
              qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
              ! for output
              qcrestot(i,k)=dum*(1._r8-dum1)
              qitend(i,k)=qitend(i,k)+dum*dum1
              qirestot(i,k)=dum*dum1
              qvlat(i,k)=qvlat(i,k)-dum
              ! for output
              qvres(i,k)=-dum
              tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
           end if
        enddo 
     enddo 

  end if

  ! calculate effective radius for pass to radiation code
  !=========================================================
  ! if no cloud water, default value is 10 micron for droplets,
  ! 25 micron for cloud ice

  ! update cloud variables after instantaneous processes to get effective radius
  ! variables are in-cloud to calculate size dist parameters
  do k=1,nlev
     do i=1,mgncol
        dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)/lcldm(i,k)
        dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)/icldm(i,k)
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)/lcldm(i,k)
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)/icldm(i,k)

        !if ((dumi(i,k).ge.qsmall).and.(dumni(i,k).lt.icsmall)) then
        !   write(iulog,*) "*** qi=",qi(i,k)," qitend=",qitend(i,k)
        !   write(iulog,*) "*** ni=",ni(i,k)," nitend=",nitend(i,k)
        !end if

        dumr(i,k) = max(qr(i,k)+qrtend(i,k)*deltat,0._r8)/precip_frac(i,k)
        dumnr(i,k) = max(nr(i,k)+nrtend(i,k)*deltat,0._r8)/precip_frac(i,k)
!--kt remove dums
!++ag
        dumg(i,k) = max(qg(i,k)+qgtend(i,k)*deltat,0._r8)
        dumng(i,k) = max(ng(i,k)+ngtend(i,k)*deltat,0._r8)
!--ag

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)
        end if

!++ag
        ! switch for specification of graupel number
        if (ngcons) then
           dumng(i,k)=ngnst/rho(i,k)*precip_frac(i,k)
        end if
!--ag


        ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
        dumc(i,k)=min(dumc(i,k),5.e-3_r8)
        dumi(i,k)=min(dumi(i,k),5.e-3_r8)
        ! limit in-precip mixing ratios
        dumr(i,k)=min(dumr(i,k),10.e-3_r8)
!--kt remove dums
!++ag
        dumg(i,k)=min(dumg(i,k),10.e-3_r8)
!--ag

     enddo
  enddo
  ! cloud ice effective radius
  !-----------------------------------------------------------------

  if (do_cldice) then
     do k=1,nlev
        do i=1,mgncol

           !if ((dumi(i,k).ge.qsmall).and.(dumni(i,k).lt.icsmall)) then
           !   write(*,*) "MG4 WARNING 4 ice present in cloud with no number conc at i=",i," k=",k
           !end if

           if ((dumi(i,k).ge.qsmall).and.(dumni(i,k).ge.icsmall)) then

              ! find index for qi (total ice mass mixing ratio)
              dum11 = (dlog10(dumi(i,k))+16._r8)*1.41328_r8
              dumit = int(dum11)

              ! set limits to make sure the calculated index doesn't exceed range of lookup table
              !             dum11 = min(dum11,real(isize))
              dum11 = min(dum11,dble(isize))
              dum11 = max(dum11,1._r8)
              dumit = max(1,dumit)
              dumit = min(isize-1,dumit)

              ! find index for Ni (ice number mixing ratio)
              dum22 = (dlog10(dumni(i,k))+10._r8)*1.10731_r8
              dumk = int(dum22)

              ! set limits to make sure the calculated index doesn't exceed range of lookup table
              !             dum22 = min(dum22,real(jsize))
              dum22 = min(dum22,dble(jsize))
              dum22 = max(dum22,1._r8)
              dumk = max(1,dumk)
              dumk = min(jsize-1,dumk)

              ! find index for scaled mean rain size
              ! if no rain, then just choose dumj = 1 and don't calculate rain-ice
              ! collection processes

              ! find index for temperature
              dum4  = t(i,k)   ! trude, and something else here for finding index qitot(i,k,iice)*3. + 1.
              if (t(i,k) .lt. 213.15_r8) then    ! (-60)
                 dum4 = 3     !trude, need to calculate the "interger" value
              else if (t(i,k).ge. 243.15_r8) then  
                 dum4 = 1                      ! (-30) 
              else
                 dum4 = 17.21_r8-(t(i,k)/15._r8) 
              endif
              dumtt = int(dum4)
              ! set limits
!             dum4  = min(dum4,real(tsize))
              dum4  = min(dum4,dble(tsize))
              dum4  = max(dum4,1._r8)
              dumtt = max(1,dumtt)
              dumtt = min(tsize-1,dumtt)

              ! call to lookup table interpolation subroutines to get process rates

              call access_lookup_table(dumtt,dumit,dumk,7,dum11,dum22,dum4,f1pr9)
              call access_lookup_table(dumtt,dumit,dumk,8,dum11,dum22,dum4,f1pr10)

              dum_2D(i,k) = dumni(i,k)
              dumni(i,k) = min(dumni(i,k),f1pr9)
              dumni(i,k) = max(dumni(i,k),f1pr10)
!              dumni(i,k) = min(dumni(i,k), 500.e3_r8/rho(i,k))
              dumni(i,k) = min(dumni(i,k), micro_mg_max_nicons/rho(i,k))
              ! dum_2D(i,k) = dumni(i,k)
              ! call size_dist_param_basic(mg_ice_props, dumi(i,k), dumni(i,k), &
              !     lami(i,k), dumni0)

              if (dumni(i,k) /=dum_2D(i,k)) then
                 ! adjust number conc if needed to keep mean size in reasonable range
                 nitend(i,k)=(dumni(i,k)*icldm(i,k)-ni(i,k))/deltat
              end if

!++ kt Trude suggests "do something about effi" Not sure about "sadice"
              !effi(i,k) = 1.5_r8/lami(i,k)*1.e6_r8
              !++kt will probably end up removing sadice.. but for now
              !if(lami(i,k).gt.0.0) then
              !   sadice(i,k) = 2._r8*pi*(lami(i,k)**(-3))*dumni0*rho(i,k)*1.e-2_r8 ! m2/m3 -> cm2/cm3
              !else
                 sadice(i,k) = 0._r8
              !end if !++kt How will this impact radiation??
              dum11 = (dlog10(dumi(i,k))+16._r8)*1.41328_r8
              dumit = int(dum11)

              ! set limits to make sure the calculated index doesn't exceed range of lookup table
              !             dum11 = min(dum11,real(isize))
              dum11 = min(dum11,dble(isize))
              dum11 = max(dum11,1._r8)
              dumit = max(1,dumit)
              dumit = min(isize-1,dumit)

              ! find index for Ni (ice number mixing ratio)
              dum22 = (dlog10(dumni(i,k))+10._r8)*1.10731_r8
              dumk = int(dum22)

              ! set limits to make sure the calculated index doesn't exceed range of lookup table
              !             dum22 = min(dum22,real(jsize))
              dum22 = min(dum22,dble(jsize))
              dum22 = max(dum22,1._r8)
              dumk = max(1,dumk)
              dumk = min(jsize-1,dumk)
              ! find index for temperature
              dum4  = t(i,k)   ! trude, and something else here for finding index qitot(i,k,iice)*3. + 1.
           if (t(i,k) .lt. 213.15_r8) then    ! (-60)
              dum4 = 3     !trude, need to calculate the "interger" value
           else if (t(i,k).ge. 243.15_r8) then  
              dum4 = 1                      ! (-30) 
           else
              dum4 = 17.21_r8-(t(i,k)/15._r8) 
           endif
              dumtt = int(dum4)
              ! set limits
              dum4  = min(dum4,dble(tsize))
              dum4  = max(dum4,1._r8)
              dumtt = max(1,dumtt)
              dumtt = min(tsize-1,dumtt)

              call access_lookup_table(dumtt,dumit,dumk,6,dum11,dum22,dum4,f1pr6)
              call access_lookup_table(dumtt,dumit,dumk,12,dum11,dum22,dum4,f1pr12)

              effi(i,k) = f1pr6*1.e6_r8   ! f1pr6 is in meter, effi is in micrometer 
!++ trude
              effi = effi*micro_mg_effi_factor
!-- trude
              !++trude, not sure if sadice is needed
             sadice(i,k)=4._r8*pi*(effi(i,k)**2)*ni(i,k)*rho(i,k)*1e-2_r8
          else
!--trude
              effi(i,k) = 25._r8
!++ trude
              effi = effi*micro_mg_effi_factor
!-- trude
              sadice(i,k) = 0._r8
           end if
!--kt
           ! ice effective diameter for david mitchell's optics
           deffi(i,k)=effi(i,k)*rhoi/rhows*2._r8
        enddo
     enddo
  else
     do k=1,nlev
        do i=1,mgncol
           ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
           ! radius has already been determined from the size distribution.
           effi(i,k) = re_ice(i,k) * 1.e6_r8      ! m -> um
!++ trude
           effi = effi*micro_mg_effi_factor
!-- trude
           deffi(i,k)=effi(i,k) * 2._r8
           sadice(i,k) = 4._r8*pi*(effi(i,k)**2)*ni(i,k)*rho(i,k)*1e-2_r8
        enddo
     enddo
  end if

  ! cloud droplet effective radius
  !-----------------------------------------------------------------
  do k=1,nlev
     do i=1,mgncol
        if (dumc(i,k).ge.qsmall) then


           ! switch for specification of droplet and crystal number
           if (nccons) then
              ! make sure nc is consistence with the constant N by adjusting tendency, need
              ! to multiply by cloud fraction
              ! note that nctend may be further adjusted below if mean droplet size is
              ! out of bounds

              nctend(i,k)=(ncnst/rho(i,k)*lcldm(i,k)-nc(i,k))/deltat

           end if

           dum = dumnc(i,k)

           call size_dist_param_liq(mg_liq_props, dumc(i,k), dumnc(i,k), rho(i,k), &
                pgam(i,k), lamc(i,k))

           if (dum /= dumnc(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nctend(i,k)=(dumnc(i,k)*lcldm(i,k)-nc(i,k))/deltat
           end if

           effc(i,k) = (pgam(i,k)+3._r8)/lamc(i,k)/2._r8*1.e6_r8
           !assign output fields for shape here
           lamcrad(i,k)=lamc(i,k)
           pgamrad(i,k)=pgam(i,k)


           ! recalculate effective radius for constant number, in order to separate
           ! first and second indirect effects
           !======================================
           ! assume constant number of 10^8 kg-1

           dumnc(i,k)=1.e8_r8

           ! Pass in "false" adjust flag to prevent number from being changed within
           ! size distribution subroutine.
           call size_dist_param_liq(mg_liq_props, dumc(i,k), dumnc(i,k), rho(i,k), &
                pgam(i,k), lamc(i,k))

           effc_fn(i,k) = (pgam(i,k)+3._r8)/lamc(i,k)/2._r8*1.e6_r8

        else
           effc(i,k) = 10._r8
           lamcrad(i,k)=0._r8
           pgamrad(i,k)=0._r8
           effc_fn(i,k) = 10._r8
        end if
     enddo
  enddo
  ! recalculate 'final' rain size distribution parameters
  ! to ensure that rain size is in bounds, adjust rain number if needed
  do k=1,nlev
     do i=1,mgncol

        if (dumr(i,k).ge.qsmall) then

           dum = dumnr(i,k)

           call size_dist_param_basic(mg_rain_props, dumr(i,k), dumnr(i,k), &
                lamr(i,k))

           if (dum /= dumnr(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nrtend(i,k)=(dumnr(i,k)*precip_frac(i,k)-nr(i,k))/deltat
           end if

        end if
     enddo
  enddo
!--kt remove snow size distribution parameters
!++ag
  ! recalculate 'final' graupel size distribution parameters
  ! to ensure that  size is in bounds, addjust number if needed
  do k=1,nlev
     do i=1,mgncol

        if (dumg(i,k).ge.qsmall) then

           dum = dumng(i,k)

           call size_dist_param_basic(mg_graupel_props, dumg(i,k), dumng(i,k), &
                lamg(i,k))

           if (dum /= dumng(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              ngtend(i,k)=(dumng(i,k)*precip_frac(i,k)-ng(i,k))/deltat
           end if

        end if
     enddo
  enddo
!--ag

  do k=1,nlev
     do i=1,mgncol
        ! if updated q (after microphysics) is zero, then ensure updated n is also zero
        !=================================================================================
        if (qc(i,k)+qctend(i,k)*deltat.lt.qsmall) nctend(i,k)=-nc(i,k)/deltat
        if (do_cldice .and. qi(i,k)+qitend(i,k)*deltat.lt.qsmall) nitend(i,k)=-ni(i,k)/deltat
        if (qr(i,k)+qrtend(i,k)*deltat.lt.qsmall) nrtend(i,k)=-nr(i,k)/deltat
!--kt remove qs
!++ag
        if (qg(i,k)+qgtend(i,k)*deltat.lt.qsmall) ngtend(i,k)=-ng(i,k)/deltat
!--ag
     end do

  end do

!  write(iulog,*) "Max qgtend: last = ",maxval(qgtend)
!  write(iulog,*) "Min qgtend: last = ",minval(qgtend)
!  write(iulog,*) "Max qvtend: last = ",maxval(qvlat)
!  write(iulog,*) "Min qvtend: last = ",minval(qvlat)


  ! DO STUFF FOR OUTPUT:
  !==================================================

  ! qc and qi are only used for output calculations past here,
  ! so add qctend and qitend back in one more time
  qc = qc + qctend*deltat
  qi = qi + qitend*deltat

  ! averaging for snow and rain number and diameter
  !--------------------------------------------------

  ! drout2/dsout2:
  ! diameter of rain and snow
  ! dsout:
  ! scaled diameter of snow (passed to radiation in CAM)
  ! reff_rain/reff_snow:
  ! calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual

  where (qrout .gt. 1.e-7_r8 &
       .and. nrout.gt.0._r8)
     qrout2 = qrout * precip_frac
     nrout2 = nrout * precip_frac
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just drout2 times constants.
     drout2 = avg_diameter(qrout, nrout, rho, rhow)
     freqr = precip_frac

     reff_rain=1.5_r8*drout2*1.e6_r8
  elsewhere
     qrout2 = 0._r8
     nrout2 = 0._r8
     drout2 = 0._r8
     freqr = 0._r8
     reff_rain = 0._r8
  end where
!--kt remove qsout

!++ag
  where (qgout .gt. 1.e-7_r8 &
       .and. ngout.gt.0._r8)
     qgout2 = qgout * precip_frac
     ngout2 = ngout * precip_frac
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just dsout2 times constants.
     dgout2 = avg_diameter(qgout, ngout, rho, rhogtmp)
     freqg = precip_frac

     dgout=3._r8*rhogtmp/rhows*dsout2

     reff_grau=1.5_r8*dsout2*1.e6_r8
  elsewhere
     dgout  = 0._r8
     qgout2 = 0._r8
     ngout2 = 0._r8
     dgout2 = 0._r8
     freqg  = 0._r8
     reff_grau=0._r8
  end where

!--ag

  ! analytic radar reflectivity
  !--------------------------------------------------
  ! formulas from Matthew Shupe, NOAA/CERES
  ! *****note: radar reflectivity is local (in-precip average)
  ! units of mm^6/m^3

  do i = 1,mgncol
     do k=1,nlev
        if (qc(i,k).ge.qsmall .and. (nc(i,k)+nctend(i,k)*deltat).gt.10._r8) then
           dum=(qc(i,k)/lcldm(i,k)*rho(i,k)*1000._r8)**2 &
                /(0.109_r8*(nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)/1.e6_r8)*lcldm(i,k)/precip_frac(i,k)
        else
           dum=0._r8
        end if
        if (qi(i,k).ge.qsmall) then
           dum1=(qi(i,k)*rho(i,k)/icldm(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm(i,k)/precip_frac(i,k)
        else
           dum1=0._r8
        end if
!--kt remove qsout

        refl(i,k)=dum+dum1

        ! add rain rate, but for 37 GHz formulation instead of 94 GHz
        ! formula approximated from data of Matrasov (2007)
        ! rainrt is the rain rate in mm/hr
        ! reflectivity (dum) is in DBz

        if (rainrt(i,k).ge.0.001_r8) then
           dum=log10(rainrt(i,k)**6._r8)+16._r8

           ! convert from DBz to mm^6/m^3

           dum = 10._r8**(dum/10._r8)
        else
           ! don't include rain rate in R calculation for values less than 0.001 mm/hr
           dum=0._r8
        end if

        ! add to refl

        refl(i,k)=refl(i,k)+dum

        !output reflectivity in Z.
        areflz(i,k)=refl(i,k) * precip_frac(i,k)

        ! convert back to DBz

        if (refl(i,k).gt.minrefl) then
           refl(i,k)=10._r8*log10(refl(i,k))
        else
           refl(i,k)=-9999._r8
        end if

        !set averaging flag
        if (refl(i,k).gt.mindbz) then
           arefl(i,k)=refl(i,k) * precip_frac(i,k)
           frefl(i,k)=precip_frac(i,k)
        else
           arefl(i,k)=0._r8
           areflz(i,k)=0._r8
           frefl(i,k)=0._r8
        end if

        ! bound cloudsat reflectivity

        csrfl(i,k)=min(csmax,refl(i,k))

        !set averaging flag
        if (csrfl(i,k).gt.csmin) then
           acsrfl(i,k)=refl(i,k) * precip_frac(i,k)
           fcsrfl(i,k)=precip_frac(i,k)
        else
           acsrfl(i,k)=0._r8
           fcsrfl(i,k)=0._r8
        end if

     end do
  end do

!--kt calculate qitendresid
  do i = 1,mgncol
     do k=1,nlev
        dumqitend(i,k) = (mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+psacws(i,k)+ &
             qmultg(i,k))*lcldm(i,k)+mnudep(i,k)*icldm(i,k)+ &
             vap_dep(i,k)+berg(i,k)+ice_sublim(i,k)+ &
             mnuccd(i,k)+(mnuccri(i,k)+pracs(i,k)+qmultrg(i,k))*precip_frac(i,k)
        dumqitend(i,k) = dumqitend(i,k)+(mnuccr(i,k)*precip_frac(i,k))+qisedten(i,k)+ &
             homotot(i,k)+qirestot(i,k)
        !--kt below are test tendencies to see where the residual is coming from
        dumqitend(i,k) = dumqitend(i,k)+rainfrztend(i,k)
        dumqitend(i,k) = dumqitend(i,k)+icemelttend(i,k)
        qitendresid(i,k) = rainfrztend(i,k)+icemelttend(i,k)
        !qitendresid(i,k) = qitend(i,k)-dumqitend(i,k)
     end do
  end do

  !redefine fice here....
!--kt remove qsout
  dum_2D = qrout + qc + qi
  dumi = qi
  where (dumi .gt. qsmall .and. dum_2D .gt. qsmall)
     nfice=min(dumi/dum_2D,1._r8)
  elsewhere
     nfice=0._r8
  end where

  !do i = 1,mgncol
  !   do k=1,nlev
  !      kttmp = qin(i,k)+qitend(i,k)*deltat
  !      if(kttmp < -1.0e-20_r8) then
  !         write(iulog,*) "**Mg4KT qi<0 qin=",qin(i,k)," tend=",qitend(i,k)*deltat 
  !      end if
  !   end do
  !end do

end subroutine micro_mg_tend

!========================================================================
!OUTPUT CALCULATIONS
!========================================================================

subroutine calc_rercld(lamr, n0r, lamc, pgam, qric, qcic, ncic, rercld, mgncol)
  integer, intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: lamr          ! rain size parameter (slope)
  real(r8), dimension(mgncol), intent(in) :: n0r           ! rain size parameter (intercept)
  real(r8), dimension(mgncol), intent(in) :: lamc          ! size distribution parameter (slope)
  real(r8), dimension(mgncol), intent(in) :: pgam          ! droplet size parameter
  real(r8), dimension(mgncol), intent(in) :: qric          ! in-cloud rain mass mixing ratio
  real(r8), dimension(mgncol), intent(in) :: qcic          ! in-cloud cloud liquid
  real(r8), dimension(mgncol), intent(in) :: ncic          ! in-cloud droplet number concentration

  real(r8), dimension(mgncol), intent(inout) :: rercld     ! effective radius calculation for rain + cloud

  ! combined size of precip & cloud drops
  real(r8) :: Atmp

  integer :: i

  do i=1,mgncol
     ! Rain drops
     if (lamr(i) > 0._r8) then
        Atmp = n0r(i) * pi / (2._r8 * lamr(i)**3._r8)
     else
        Atmp = 0._r8
     end if

     ! Add cloud drops
     if (lamc(i) > 0._r8) then
        Atmp = Atmp + &
             ncic(i) * pi * rising_factorial(pgam(i)+1._r8, 2)/(4._r8 * lamc(i)**2._r8)
     end if

     if (Atmp > 0._r8) then
        rercld(i) = rercld(i) + 3._r8 *(qric(i) + qcic(i)) / (4._r8 * rhow * Atmp)
     end if
  enddo
end subroutine calc_rercld

!========================================================================
!UTILITIES
!========================================================================


pure subroutine micro_mg_get_cols(ncol, nlev, top_lev, qcn, qin, &
!++ag
     qrn, qsn, qgr,mgncol, mgcols)
!--ag

  ! Determines which columns microphysics should operate over by
  ! checking for non-zero cloud water/ice.

  integer, intent(in) :: ncol      ! Number of columns with meaningful data
  integer, intent(in) :: nlev      ! Number of levels to use
  integer, intent(in) :: top_lev   ! Top level for microphysics

  real(r8), intent(in) :: qcn(:,:) ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:) ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: qrn(:,:) ! rain mixing ratio (kg/kg)
  real(r8), intent(in) :: qsn(:,:) ! snow mixing ratio (kg/kg)
!++ag
  real(r8), intent(in) :: qgr(:,:) ! graupel mixing ratio (kg/kg)
!--ag

  integer, intent(out) :: mgncol   ! Number of columns MG will use
  integer, allocatable, intent(out) :: mgcols(:) ! column indices

  integer :: lev_offset  ! top_lev - 1 (defined here for consistency)
  logical :: ltrue(ncol) ! store tests for each column

  integer :: i, ii ! column indices

  if (allocated(mgcols)) deallocate(mgcols)

  lev_offset = top_lev - 1

  ! Using "any" along dimension 2 collapses across levels, but
  ! not columns, so we know if water is present at any level
  ! in each column.

  ltrue = any(qcn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qin(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qrn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qsn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
!++ag
  ltrue = ltrue .or. any(qgr(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
!--ag

  ! Scan for true values to get a usable list of indices.

  mgncol = count(ltrue)
  allocate(mgcols(mgncol))
  i = 0
  do ii = 1,ncol
     if (ltrue(ii)) then
        i = i + 1
        mgcols(i) = ii
     end if
  end do

end subroutine micro_mg_get_cols

end module micro_mg4_0

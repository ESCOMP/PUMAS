! *****************************COPYRIGHT*******************************
! based on mphys_morr_two_moment
! *****************************COPYRIGHT*******************************
!
! Module containing interface to Morrison & Gettelman version 2 (mg2)
!  Andrew Gettelman, February 2013
!
module mphys_mg1_5

  Use parameters, only : num_h_moments, num_h_bins, nspecies, nz, dt &
       , h_names, mom_units, max_char_len, mom_names, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi

  Use module_mp_mg1_5, only: micro_mg_init, micro_mg_tend
  Use diagnostics, only: save_dg, i_dgtime
  Use common_physics, only : qsaturation, qisaturation

  Use wv_sat_methods, only: wv_sat_qsat_water, wv_sat_qsat_ice, wv_sat_methods_init

  Implicit None

  !Logical switches 
  logical :: micro_unset=.True.
  integer:: ih, imom
  character(max_char_len) :: name, units

contains
  
  Subroutine mphys_mg1_5_interface
    
    real :: t1d(nz), p1d(nz), dz1d(nz),qv1d(nz),qc1d(nz) &
         , qi1d(nz), ni1d(nz) & 
         , qr1d(nz), qs1d(nz)         &
         , ns1d(nz), nr1d(nz)        &
         , w1d(nz)

    real :: wvar1d(nz)

    real :: qc_tend1d(nz), ni_tend1d(nz)  &
         ,  qi_tend1d(nz), qni_tend1d(nz) & 
         ,  qv_tend1d(nz), t_tend1d(nz)   &
         ,  qr_tend1d(nz), nr_tend1d(nz) &
         ,  qs_tend1d(nz), ns_tend1d(nz) 
!         qg_tend1d(nz), ng_tend1d(nz)

    real :: precprt1d, snowrt1d

    real :: effc1d(nz), effi1d(nz), effs1d(nz),        &
         effr1d(nz), effg1d(nz)

    real :: qrcu1d(nz), qscu1d(nz), qicu1d(nz)

    real :: qgsten(nz), qrsten(nz), qisten(nz),  &
         qnisten(nz), qcsten(nz)

    ! KiD_2D diag arrays
    real :: qrsten_2d(nz,nx),precprt2d(nx), snowrt2d(nx) 

    integer :: kts, kte, i, j, k

!++ag
    
!Init Variables
    character(128) :: errstring
    logical :: microp_uniform_in,do_cldice_in     
    integer, parameter :: kind = selected_real_kind(12) ! 8 byte real
    integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
    real(R8) :: gravit,rair,rh2o,cpair,tmelt_in,latvap,latice,rhmini_in

    real(R8) :: h2otrip_in, tboil_in, ttrice_in, epsilo_in 

    ! From shr_const_mod.F90
    real(R8),parameter :: SHR_CONST_G       = 9.80616_R8      ! acceleration of gravity ~ m/s^2
    real(R8),parameter :: SHR_CONST_BOLTZ   = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
    real(R8),parameter :: SHR_CONST_AVOGAD  = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole
    real(R8),parameter :: SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
    real(R8),parameter :: SHR_CONST_MWDAIR  = 28.966_R8       ! molecular weight dry air ~ kg/kmole
    real(R8),parameter :: SHR_CONST_MWWV    = 18.016_R8       ! molecular weight water vapor
    real(R8),parameter :: SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant     ~ J/K/kg
    real(R8),parameter :: SHR_CONST_RWV     = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
    real(R8),parameter :: SHR_CONST_TKFRZ   = 273.15_R8       ! freezing T of fresh water          ~ K 
    real(R8),parameter :: SHR_CONST_LATICE  = 3.337e5_R8      ! latent heat of fusion      ~ J/kg
    real(R8),parameter :: SHR_CONST_LATVAP  = 2.501e6_R8      ! latent heat of evaporation ~ J/kg
    real(R8),parameter :: SHR_CONST_CPDAIR  = 1.00464e3_R8    ! specific heat of dry air   ~ J/kg/K

    real(R8),parameter :: qsmall = 1.e-12_R8

! DUMMY VARIABLES

     REAL(R8) DUMT,DUMQV,DUMQSS,DUMQC,PCC,DUMS, DNC,DUMNC,ESL,QVL,CPM

!Inputs and outputs.
  ! input arguments
  integer :: mgncol                ! number of microphysics columns
  integer :: mgcols(nx)             ! list of microphysics columns
  integer :: nlev                  ! number of layers
  integer :: top_lev               ! top level to do microphysics
  real(r8) :: deltatin              ! time step (s)
  real(r8) :: tn(nx,nz)        ! input temperature (K)
  real(r8) :: qn(nx,nz)        ! input h20 vapor mixing ratio (kg/kg)

  ! note: all input cloud variables are grid-averaged
  real(r8) :: qcn(nx,nz)       ! cloud water mixing ratio (kg/kg)
  real(r8) :: qin(nx,nz)       ! cloud ice mixing ratio (kg/kg)
  real(r8) :: ncn(nx,nz)       ! cloud water number conc (1/kg)
  real(r8) :: nin(nx,nz)       ! cloud ice number conc (1/kg)
  real(r8) :: pn(nx,nz)         ! air pressure (pa)
  real(r8) :: pdeln(nx,nz)      ! pressure difference across level (pa)
  ! hm add 11-16-11, interface pressure
  real(r8) :: pint(nx,0:nz)    ! level interface pressure (pa)
  real(r8) :: cldn(nx,nz)      ! cloud fraction (no units)
  real(r8) :: liqcldf(nx,nz)   ! liquid cloud fraction (no units)
  real(r8) :: icecldf(nx,nz)   ! ice cloud fraction (no units)
  ! used for scavenging
  ! Inputs for aerosol activation
  real(r8) :: naain(nx,nz)     ! ice nucleation number (from microp_aero_ts) (1/kg)
  real(r8) :: npccnin(nx,nz)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  ! Note that for these variables, the dust bin is assumed to be the last index.
  ! (For example, in CAM, the last dimension is always size 4.)
  real(r8) :: rndstn(nx,nz,4)  ! radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
  real(r8) :: naconin(nx,nz,4) ! number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)

  ! Used with CARMA cirrus microphysics
  ! (or similar external microphysics model)
  real(r8) :: tnd_qsnown(nx,nz) ! snow mass tendency (kg/kg/s)
  real(r8) :: tnd_nsnown(nx,nz) ! snow number tendency (#/kg/s)
  real(r8) :: re_icen(nx,nz)    ! ice effective radius (m)

  ! output arguments

  real(r8) :: rate1ord_cw2pr_st(nx,nz)    ! 1st order rate for
  ! direct cw to precip conversion
  real(r8) :: tlato(nx,nz)         ! latent heating rate       (W/kg)
  real(r8) :: qvlato(nx,nz)        ! microphysical tendency qv (1/s)
  real(r8) :: qctendo(nx,nz)       ! microphysical tendency qc (1/s)
  real(r8) :: qitendo(nx,nz)       ! microphysical tendency qi (1/s)
  real(r8) :: nctendo(nx,nz)       ! microphysical tendency nc (1/(kg*s))
  real(r8) :: nitendo(nx,nz)       ! microphysical tendency ni (1/(kg*s))
  real(r8) :: effco(nx,nz)         ! droplet effective radius (micron)
  real(r8) :: effco_fn(nx,nz)      ! droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8) :: effio(nx,nz)         ! cloud ice effective radius (micron)
  real(r8) :: precto(nx)          ! surface precip rate (m/s)
  real(r8) :: precio(nx)          ! cloud ice/snow precip rate (m/s)
  real(r8) :: nevapro(nx,nz)       ! evaporation rate of rain + snow (1/s)
  real(r8) :: evapsnowo(nx,nz)     ! sublimation rate of snow (1/s)
  real(r8) :: praino(nx,nz)        ! production of rain + snow (1/s)
  real(r8) :: prodsnowo(nx,nz)     ! production of snow (1/s)
  real(r8) :: cmeouto(nx,nz)       ! evap/sub of cloud (1/s)
  real(r8) :: deffio(nx,nz)        ! ice effective diameter for optics (radiation) (micron)
  real(r8) :: pgamrado(nx,nz)      ! ice gamma parameter for optics (radiation) (no units)
  real(r8) :: lamcrado(nx,nz)      ! slope of droplet distribution for optics (radiation) (1/m)
  real(r8) :: qsouto(nx,nz)        ! snow mixing ratio (kg/kg)
  real(r8) :: dsouto(nx,nz)        ! snow diameter (m)
  real(r8) :: rflxo(nx,nz)         ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8) :: sflxo(nx,nz)         ! grid-box average snow flux (kg m^-2 s^-1)
  real(r8) :: qrouto(nx,nz)        ! grid-box average rain mixing ratio (kg/kg)
  real(r8) :: reff_raino(nx,nz)    ! rain effective radius (micron)
  real(r8) :: reff_snowo(nx,nz)    ! snow effective radius (micron)
  real(r8) :: qcsevapo(nx,nz)      ! cloud water evaporation due to sedimentation (1/s)
  real(r8) :: qisevapo(nx,nz)      ! cloud ice sublimation due to sublimation (1/s)
  real(r8) :: qvreso(nx,nz)        ! residual condensation term to ensure RH < 100% (1/s)
  real(r8) :: cmeiout(nx,nz)       ! grid-mean cloud ice sub/dep (1/s)
  real(r8) :: vtrmco(nx,nz)        ! mass-weighted cloud water fallspeed (m/s)
  real(r8) :: vtrmio(nx,nz)        ! mass-weighted cloud ice fallspeed (m/s)
  real(r8) :: qcsedteno(nx,nz)     ! qc sedimentation tendency (1/s)
  real(r8) :: qisedteno(nx,nz)     ! qi sedimentation tendency (1/s)

  ! microphysical process rates for output (mixing ratio tendencies) (all have units of 1/s)
  real(r8) :: prao(nx,nz)         ! accretion of cloud by rain 
  real(r8) :: prco(nx,nz)         ! autoconversion of cloud to rain
  real(r8) :: mnuccco(nx,nz)      ! mixing ratio tend due to immersion freezing
  real(r8) :: mnuccto(nx,nz)      ! mixing ratio tend due to contact freezing
  real(r8) :: msacwio(nx,nz)      ! mixing ratio tend due to H-M splintering
  real(r8) :: psacwso(nx,nz)      ! collection of cloud water by snow
  real(r8) :: bergso(nx,nz)       ! bergeron process on snow
  real(r8) :: bergo(nx,nz)        ! bergeron process on cloud ice
  real(r8) :: melto(nx,nz)        ! melting of cloud ice
  real(r8) :: homoo(nx,nz)        ! homogeneous freezing cloud water
  real(r8) :: qcreso(nx,nz)       ! residual cloud condensation due to removal of excess supersat
  real(r8) :: prcio(nx,nz)        ! autoconversion of cloud ice to snow
  real(r8) :: praio(nx,nz)        ! accretion of cloud ice by snow
  real(r8) :: qireso(nx,nz)       ! residual ice deposition due to removal of excess supersat
  real(r8) :: mnuccro(nx,nz)      ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
  real(r8) :: pracso(nx,nz)       ! mixing ratio tendency due to accretion of rain by snow (1/s)
  real(r8) :: meltsdto(nx,nz)     ! latent heating rate due to melting of snow  (W/kg)
  real(r8) :: frzrdto(nx,nz)      ! latent heating rate due to homogeneous freezing of rain (W/kg)
  real(r8) :: mnuccdo(nx,nz)      ! mass tendency from ice nucleation
  real(r8) :: nrouto(nx,nz)        ! rain number concentration (1/m3)
  real(r8) :: nsouto(nx,nz)        ! snow number concentration (1/m3)
  real(r8) :: reflo(nx,nz)         ! analytic radar reflectivity        
  real(r8) :: areflo(nx,nz)        ! average reflectivity will zero points outside valid range
  real(r8) :: areflzo(nx,nz)       ! average reflectivity in z.
  real(r8) :: freflo(nx,nz)        ! fractional occurrence of radar reflectivity
  real(r8) :: csrflo(nx,nz)        ! cloudsat reflectivity 
  real(r8) :: acsrflo(nx,nz)       ! cloudsat average
  real(r8) :: fcsrflo(nx,nz)       ! cloudsat fractional occurrence of radar reflectivity
  real(r8) :: rercldo(nx,nz)       ! effective radius calculation for rain + cloud
  real(r8) :: ncaio(nx,nz)        ! output number conc of ice nuclei available (1/m3)
  real(r8) :: ncalo(nx,nz)        ! output number conc of CCN (1/m3)
  real(r8) :: qrouto2(nx,nz)       ! copy of qrout as used to compute drout2
  real(r8) :: qsouto2(nx,nz)       ! copy of qsout as used to compute dsout2
  real(r8) :: nrouto2(nx,nz)       ! copy of nrout as used to compute drout2
  real(r8) :: nsouto2(nx,nz)       ! copy of nsout as used to compute dsout2
  real(r8) :: drouto2(nx,nz)       ! mean rain particle diameter (m)
  real(r8) :: dsouto2(nx,nz)       ! mean snow particle diameter (m)
  real(r8) :: freqso(nx,nz)        ! fractional occurrence of snow
  real(r8) :: freqro(nx,nz)        ! fractional occurrence of rain
  real(r8) :: nficeo(nx,nz)        ! fractional occurrence of ice
  real(r8) :: qcrato(nx,nz)        ! limiter for qc process rates (1=no limit --> 0. no qc)

! Flipped arrays Input:  tn,qn,qcn,ncn,qin,nin,cldn,liqcldf,icecldf
  real(r8) ::flip_pn(nx,nz)
  real(r8) ::flip_pdeln(nx,nz)
  real(r8) ::flip_pint(nx,0:nz)
  real(r8) ::flip_tn(nx,nz)
  real(r8) ::flip_qn(nx,nz)
  real(r8) ::flip_qcn(nx,nz)
  real(r8) ::flip_ncn(nx,nz)
  real(r8) ::flip_qin(nx,nz)
  real(r8) ::flip_nin(nx,nz) 
  real(r8) ::flip_cldn(nx,nz)
  real(r8) ::flip_liqcldf(nx,nz)
  real(r8) ::flip_icecldf(nx,nz)  
  real(r8) ::flip_naain(nx,nz) 
  real(r8) ::flip_npccnin(nx,nz)
  real(r8) ::flip_relvarn(nx,nz)
  real(r8) ::flip_accre_enhann(nx,nz)

!Flipped arrays output
!tlato,qvlato,qctendo,qitendo,nctendo,nitendo

!tendencies
  real(r8) ::flip_tlato(nx,nz)
  real(r8) ::flip_qvlato(nx,nz)
  real(r8) ::flip_qctendo(nx,nz)
  real(r8) ::flip_qitendo(nx,nz)
  real(r8) ::flip_nitendo(nx,nz)
  real(r8) ::flip_nctendo(nx,nz)
!other (for output)
  real(r8) ::flip_qrouto(nx,nz)
  real(r8) ::flip_nrouto(nx,nz) 
  real(r8) ::flip_qrouto2(nx,nz)
  real(r8) ::flip_nrouto2(nx,nz)  
  real(r8) ::flip_prao(nx,nz)
  real(r8) ::flip_prco(nx,nz)
  real(r8) ::flip_cmeiout(nx,nz)
!--ag

!note different order of i,k for MG...
!also: MG is top down, kid arrays are bottom up.    

    kts=1
    kte=nz
    j=1

    qrouto(:,:)= 0.
    nrouto(:,:)= 0.
    flip_qrouto2(:,:)=0.
    flip_nrouto2(:,:)=0.

    do i=1,nx
       !zero out precipitation
       precto(i) = 0.
       precio(i) = 0.

       !interface pressures
       pint(i,:) = 100. * pmb_half(:,i) 

       do k=1,nz
       ! zero some of these for safety
          tlato(i,k)   = 0.
          qvlato(i,k)  = 0.
          qctendo(i,k)  = 0.
          qitendo(i,k)  = 0.
          nctendo(i,k) = 0.
          nitendo(i,k)  = 0.

          effco(i,k)     = 0.
          effio(i,k)     = 0.
          effco_fn(i,k)     = 0.

          !set input arguments

          !temperature

          tn(i,k) = (theta(k,i) + (dtheta_adv(k,i)+dtheta_div(k,i))*dt )*exner(k,i)
          !layer pressure
          pn(i,k) = p0*exner(k,i)**(1./r_on_cp)

          !interface pressure (assume we are going from the bottom up)
          pdeln(i,k)= pint(i,k-1) - pint(i,k)
          
          !corrected: use dp = rho * g * dz = P/RdT * g *dz
          pdeln(i,k)= pn(i,k)/(287.15*tn(i,k)) * 9.81 * dz(k)

          qn(i,k) = qv(k,i)+ (dqv_adv(k,i)+dqv_div(k,i))*dt
          if (qn(i,k).lt.qsmall) qn(i,k)=0.
          qcn(i,k) = hydrometeors(k,i,1)%moments(1,1) & 
               + (dhydrometeors_adv(k,i,1)%moments(1,1) &
               + dhydrometeors_div(k,i,1)%moments(1,1))*dt
          if (qcn(i,k).lt.qsmall) qcn(i,k)=0.
          if (num_h_moments(1) >= 2) &
               ncn(i,k) = hydrometeors(k,i,1)%moments(1,2)& 
               + (dhydrometeors_adv(k,i,1)%moments(1,2) &
               + dhydrometeors_div(k,i,1)%moments(1,2))*dt
          if (qcn(i,k).lt.qsmall) ncn(i,k) = 0.
          if (num_h_moments(3) >= 1) &
               qin(i,k) = hydrometeors(k,i,3)%moments(1,1)& 
               + (dhydrometeors_adv(k,i,3)%moments(1,1) &
               + dhydrometeors_div(k,i,3)%moments(1,1))*dt
          if (qin(i,k).lt.qsmall) qin(i,k)=0.
          if (num_h_moments(3) >= 2) &
               nin(i,k) = hydrometeors(k,i,3)%moments(1,2)& 
               + (dhydrometeors_adv(k,i,3)%moments(1,2) &
               + dhydrometeors_div(k,i,3)%moments(1,2))*dt
          if (qin(i,k).lt.qsmall) nin(i,k) = 0.
!set cloud information
          liqcldf(i,k)= 0.
          icecldf(i,k)= 0.
          cldn(i,k)=0.

          if (qcn(i,k) > qsmall) liqcldf(i,k) = 1.0_r8
          if (qin(i,k) > qsmall) icecldf(i,k) = 1.0_r8
          if (qcn(i,k) > qsmall .or. qin(i,k) > qsmall) cldn(i,k) = 1.0_r8

!set activated aerosol number (uniform for now)
!          if (tn(i,k) > tmelt_in) &
               npccnin = 10.e6   ! #/m3   
!          if (tn(i,k) < tmelt_in - 10.) &
               naain = 0.e6   ! #/m3             

!          if (num_h_moments(3) >= 1) &
!               qr1d(k) = hydrometeors(k,i,3)%moments(1,1)& 
!               + (dhydrometeors_adv(k,i,3)%moments(1,1) &
!               + dhydrometeors_div(k,i,3)%moments(1,1))*dt
!          if (num_h_moments(3) >= 2) &
!               nr1d(k) = hydrometeors(k,i,3)%moments(1,2)& 
!               + (dhydrometeors_adv(k,i,3)%moments(1,2) &
!               + dhydrometeors_div(k,i,3)%moments(1,2))*dt
!          if (num_h_moments(4) >= 1) &
!               qs1d(k) = hydrometeors(k,i,4)%moments(1,1)& 
!               + (dhydrometeors_adv(k,i,4)%moments(1,1) &
!               + dhydrometeors_div(k,i,4)%moments(1,1))*dt
!          if (num_h_moments(4) >= 2) &
!               ns1d(k) = hydrometeors(k,i,4)%moments(1,2)& 
!               + (dhydrometeors_adv(k,i,4)%moments(1,2) &
!               + dhydrometeors_div(k,i,4)%moments(1,2))*dt
!          if (num_h_moments(5) >= 1) &
!               qg1d(k) = hydrometeors(k,i,5)%moments(1,1)& 
!               + (dhydrometeors_adv(k,i,5)%moments(1,1) &
!               + dhydrometeors_div(k,i,5)%moments(1,1))*dt
!          if (num_h_moments(5) >= 2) &
!               ng1d(k) = hydrometeors(k,i,5)%moments(1,2)& 
!               + (dhydrometeors_adv(k,i,5)%moments(1,2) &
!               + dhydrometeors_div(k,i,5)%moments(1,2))*dt         
!          wvar1d(k) = 0.5 ! hard-wired not coupled to forcing!
!          w1d(k) = w_half(k,i)
       end do
    end do


!print*,'pn(1),pn(nz),pint(0),pint(nz),pdel(1),pdel(nz)=',pn(1,1),pn(1,nz),pint(1,0),pint(1,nz),pdeln(1,1),pdeln(1,nz)

!set up constants for microphysics here...
          gravit =  SHR_CONST_G !Gravity
          rair =    SHR_CONST_RDAIR ! dry air gas constant: note units (phys_constants are in J/K/kmol)
          rh2o =    SHR_CONST_RWV ! water vapor gas constant
          cpair =    SHR_CONST_CPDAIR ! specific heat of dry air J/kg/K
          tmelt_in = SHR_CONST_TKFRZ !Freezing point of Water (K)
          latvap =   SHR_CONST_LATVAP ! latent heat vaporization J/kg
          latice=    SHR_CONST_LATICE ! latent heat freezing J/kg

       ! Initialise microphysics 
       if (micro_unset)then
!++ag
          rhmini_in =  0.7 ! Minimum rh for ice cloud fraction > 0.
          microp_uniform_in = .true.
          do_cldice_in = .true.

          h2otrip_in = 273.16
          tboil_in = 373.16
          ttrice_in = 20.
          epsilo_in =  SHR_CONST_MWWV / SHR_CONST_MWDAIR  

          call wv_sat_methods_init(kind, tmelt_in, h2otrip_in, tboil_in, &
     ttrice_in, epsilo_in, errstring)
          
          call micro_mg_init( &
     kind, gravit, rair, rh2o, cpair,    &
     tmelt_in, latvap, latice,           &
     rhmini_in, microp_uniform_in, do_cldice_in, &
     errstring)
!--ag
          micro_unset=.False.
       end if


!       call morr_two_moment_micro(qc_tend1d, qi_tend1d, qni_tend1d, &
!            qr_tend1d, ni_tend1d, ns_tend1d, nr_tend1d,             &
!            qc1d, qi1d, qs1d, qr1d, ni1d, ns1d, nr1d,               &
!            t_tend1d, qv_tend1d, t1d, qv1d, p1d, dz1d, w1d, wvar1d, &
!            precprt1d, snowrt1d,                                    &
!            effc1d, effi1d, effs1d, effr1d, dt,                     &
!            i,i,j,j,kts,kte,                                        &
!            i,i,j,j,kts,kte,                                        &
!            qg_tend1d, ng_tend1d, qg1d, ng1d, effg1d,               &
!            qrcu1d, qscu1d, qicu1d,                                 &
!            qgsten, qrsten, qisten, qnisten, qcsten)


!microphysics routine for each timestep goes here...

       mgncol=nx
       mgcols(:) = 1
       nlev=nz
       top_lev=1  
       deltatin=dt

!set dust information (for contact freezing):
!size
       rndstn(:,:,1)  = 0.5e-6      ! radius (m)
       rndstn(:,:,2)  = 1.e-6      ! radius (m)
       rndstn(:,:,3)  = 2.e-6      ! radius (m)
       rndstn(:,:,4)  = 10.e-6      ! radius (m)
!number
       naconin(:,:,:) = 1000.  ! 100 m-3

!now: need to flip all inputs in the vertical: MG assumes 1=top , nz=surf
! Flipped arrays Input:  tn,qn,qcn,ncn,qin,nin,cldn,liqcldf,icecldf

       flip_pn(:,:)=pn(:,nz:1:-1)
       flip_pdeln(:,:)=pdeln(:,nz:1:-1)
       flip_pint(:,:)=pint(:,nz:0:-1)
       flip_tn(:,:)=tn(:,nz:1:-1)
       flip_qn(:,:)=qn(:,nz:1:-1)
       flip_qcn(:,:)=qcn(:,nz:1:-1)
       flip_qin(:,:)=qin(:,nz:1:-1)
       flip_ncn(:,:)=ncn(:,nz:1:-1)
       flip_nin(:,:)=nin(:,nz:1:-1)
       flip_cldn(:,:)=cldn(:,nz:1:-1)
       flip_liqcldf(:,:)=liqcldf(:,nz:1:-1)
       flip_icecldf(:,:)=icecldf(:,nz:1:-1)

!Flipped activated aerosol number concentration
       flip_naain(:,:)=naain(:,nz:1:-1)
       flip_npccnin(:,:)=npccnin(:,nz:1:-1)

!Accretion and Relative variance for sub-grid processes...now passed in
!for KiD: set to 1 for accretion, 0. for relative variance.
       flip_relvarn(:,:)= 0.
       flip_accre_enhann(:,:)= 1.

!      i=1
!       k=25
!      do k=1,nz 
!         print*,'z,pmb,tn,qn,qcn,ncn,qin,nin=',&
!              z(k),pmb(k,i),tn(i,k),qn(i,k),qcn(i,k),ncn(i,k),qin(i,k),&
!              nin(i,k)
!      end do

!,qn(i,k),qcn(i,k),ncn(i,k),qin(i,k),&
!                  nin(i,k)

call micro_mg_tend ( &
!Input
     mgncol,   mgcols,   nlev,     top_lev,  deltatin,           &
     flip_tn,                 flip_qn,               &
     flip_qcn,                          flip_qin,    &
     flip_ncn,                          flip_nin,    &
     flip_relvarn,            flip_accre_enhann,     &     
     flip_pn,                 flip_pdeln,              flip_pint,      &
     flip_cldn,               flip_liqcldf,            flip_icecldf,   &
!!Output
     rate1ord_cw2pr_st,  flip_naain,    flip_npccnin,  rndstn,   naconin,  &
     flip_tlato,    flip_qvlato,   flip_qctendo,  flip_qitendo,  flip_nctendo,  flip_nitendo,  &
     effco,    effco_fn, effio,              precto,   precio,   &
     nevapro, evapsnowo, praino,  prodsnowo, cmeouto,  deffio,   &
     pgamrado, lamcrado, qsouto,   dsouto,   rflxo,    sflxo,    &
     flip_qrouto,             reff_raino,         reff_snowo,         &
     qcsevapo, qisevapo, qvreso,   flip_cmeiout,  vtrmco,   vtrmio,   &
     qcsedteno,qisedteno,flip_prao,     flip_prco,     mnuccco,  mnuccto,  &
     msacwio,  psacwso,  bergso,   bergo,    melto,    homoo,    &
     qcreso,             prcio,    praio,    qireso,             &
     mnuccro,  pracso,   meltsdto, frzrdto,  mnuccdo,            &
     flip_nrouto,   nsouto,   reflo,    areflo,   areflzo,  freflo,   &
     csrflo,   acsrflo,  fcsrflo,            rercldo,            &
     ncaio,    ncalo,    flip_qrouto2,  qsouto2,  flip_nrouto2,  nsouto2,  &
     drouto2,  dsouto2,  freqso,   freqro,   nficeo,  qcrato,    &
     tnd_qsnown,         tnd_nsnown,         re_icen,            &
     errstring)

  !Authors: Hugh Morrison, Andrew Gettelman, NCAR, Peter Caldwell, LLNL
  ! e-mail: morrison@ucar.edu, andrew@ucar.edu

!Rembember: need to flip all output in the vertical

! Key tendencies to hook up prognostically:
!tlato,qvlato,qctendo,qitendo,nctendo,nitendo
       
!note: flip_tlato is dry static energy: W/kg divide by cpair

       tlato(:,:)=flip_tlato(:,nz:1:-1) / cpair
       qvlato(:,:)=flip_qvlato(:,nz:1:-1)
       qctendo(:,:)=flip_qctendo(:,nz:1:-1)
       qitendo(:,:)=flip_qitendo(:,nz:1:-1)
       nctendo(:,:)=flip_nctendo(:,nz:1:-1)
       nitendo(:,:)=flip_nitendo(:,nz:1:-1)

!       i=1
!       k=25
!       print*,'z,pmb,t,tlato,qvlato,qctendo,nctendo,qitendo,nitendo=',&
!                  z(k),pmb(k,i),tn(i,k),tlato(i,k),qvlato(i,k),qctendo(i,k),nctendo(i,k),&
!                  qitendo(i,k),nitendo(i,k)


! Flip for output

       nrouto(:,:)=flip_nrouto2(:,nz:1:-1)
       qrouto(:,:)=flip_qrouto2(:,nz:1:-1)
       prco(:,:)=flip_prco(:,nz:1:-1)
       prao(:,:)=flip_prao(:,nz:1:-1)
       cmeiout(:,:)=flip_cmeiout(:,nz:1:-1)

!Add saturation adjustment...based on m2005
       do i=1,nx
          do k=1,nz
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION  

             DUMT = tn(i,k) + dt*tlato(i,k)
             DUMQV = qn(i,k) + dt*qvlato(i,k)

             call wv_sat_qsat_water(DUMT, pn(i,k), ESL, QVL, 1)

!               DUMQSS = qsaturation(tn(k,i),pmb(k,i))
             DUMQSS = QVL

             DUMQC = qcn(i,k) + dt*qctendo(i,k)
             DUMQC = MAX(DUMQC,0.)
             DUMNC = ncn(i,k) + dt*nctendo(i,k)

! SATURATION ADJUSTMENT FOR LIQUID
             PCC= 0.

             DUMS = DUMQV-DUMQSS

             CPM = cpair*(1.+0.887*DUMQC) 

!Do condensation only (doesnt seem to work without a positive limit warm2: why doesn't evap work?)
!             if (DUMS.GT.qsmall) &
                  PCC = DUMS / (1.+ latvap**2 * DUMQSS/(CPM * rh2o * DUMT**2)) / dt

!limit evap to qc
             IF (PCC*dt + DUMQC.LT.qsmall) THEN
                PCC = -DUMQC/dt
             END IF

!Figure out additional number concentration (assume size of 8 microns...)
             DNC =  3._r8 * PCC / (4._r8*3.14_r8* 8.e-6_r8**3*997._r8)

!make sure it does not get less than zero
             IF (DNC*dt + DUMNC.LT.qsmall) THEN
                DNC = -DUMNC/dt
             END IF         

            
!             if (i.eq.1.and.k.eq. 25) &
!                  print*, 'z,pmb,DUMT,DUMqv,DUMQS,DUMQC,pcc,dnc,tlato=',&
!                  z(k),pmb(k,i),DUMT,DUMQV,QVL,DUMQC,PCC,DNC,tlato(k,i)

!apply tendencies
             qvlato(i,k) = qvlato(i,k) - PCC
             tlato(i,k) = tlato(i,k) + PCC * latvap/CPM
             qctendo(i,k) = qctendo(i,k) + PCC
             nctendo(i,k) = nctendo(i,k) + DNC 

!limters to make sure if (a) non negative mass and number and (b) no mass, then no number
             if (qcn(i,k)+qctendo(i,k).lt.qsmall) then
                qctendo(i,k)=-qcn(i,k)/dt
                nctendo(i,k)=-ncn(i,k)/dt
             end if

             if (qin(i,k)+qitendo(i,k).lt.qsmall) then
                qitendo(i,k)=-qin(i,k)/dt
                nitendo(i,k)=-nin(i,k)/dt
             end if


             if (i.eq.1.and.k.eq. 25) &
                  print*,'z,pmb,t,tlato,qvlato,qctendo,nctendo,qitendo,nitendo=',&
                  z(k),pmb(k,i),tn(i,k),tlato(i,k),qvlato(i,k),qctendo(i,k),nctendo(i,k),&
                  qitendo(i,k),nitendo(i,k)

             
          end do
       end do

!everything else is diagnostic....
      
!       qrsten_2d(:,i)=qrsten(:)
       precprt2d = precto
       snowrt2d = precio
   
!       print*,'precto,precio',precto(1),precio(1)
   
       ! save tendencies
  do i=1,nx
       do k=1,nz
          dtheta_mphys(k,i)=tlato(i,k)/exner(k,i)
          dqv_mphys(k,i)=qvlato(i,k)
          dhydrometeors_mphys(k,i,1)%moments(1,1)= qctendo(i,k)
          dhydrometeors_mphys(k,i,1)%moments(1,2)= nctendo(i,k)
          dhydrometeors_mphys(k,i,3)%moments(1,1)= qitendo(i,k)
          dhydrometeors_mphys(k,i,3)%moments(1,2)= nitendo(i,k)
       end do   
    end do

    ! Save some diagnostics

!Rain/Snow (what are units?)

    name='total_surface_ppt'
    units='m s-1'

    if (nx == 1) then
       call save_dg(precto(1), name, i_dgtime,  units, dim='time')
    else
       call save_dg(precto(1:nx), name, i_dgtime,  units, dim='time')
    endif
 
    name='surface_ppt_for_snow'
    units='m s-1'
    if (nx == 1) then
       call save_dg(precio(1), name, i_dgtime,  units, dim='time')
    else
       call save_dg(precio(1:nx), name, i_dgtime,  units, dim='time')
    endif

    name='rain_mass'
    units='kg kg-1'
    if (nx == 1) then
       call save_dg(qrouto(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(qrouto(1:nx,:), name, i_dgtime,  units, dim='z')
    endif   

    name='rain_number'
    units='kg-1'
    if (nx == 1) then
       call save_dg(nrouto(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(nrouto(1:nx,:), name, i_dgtime,  units, dim='z')
    endif   

    name='snow_mass'
    units='kg kg-1'
    if (nx == 1) then
       call save_dg(qsouto(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(qsouto(1:nx,:), name, i_dgtime,  units, dim='z')
    endif   

    name='snow_number'
    units='kg-1'
    if (nx == 1) then
       call save_dg(nsouto(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(nsouto(1:nx,:), name, i_dgtime,  units, dim='z')
    endif   


    name='autoconversion_rate'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(prco(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(prco(1:nx,:), name, i_dgtime,  units, dim='z')
    endif   

    name='accretion_rate'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(prao(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(prao(1:nx,:), name, i_dgtime,  units, dim='z')
    endif   

    name='ice_condensation_rate'
    units='kg kg-1 s-1'
    if (nx == 1) then
       call save_dg(cmeiout(1,:), name, i_dgtime,  units, dim='z')
    else
       call save_dg(cmeiout(1:nx,:), name, i_dgtime,  units, dim='z')
    endif   


  end Subroutine mphys_mg1_5_interface

end module mphys_mg1_5

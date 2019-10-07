! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Switches
!
Module switches

  Use typeKind
  Use parameters, only : nspecies,naerosol  & 
     , num_aero_moments, num_aero_bins  
  Implicit None

  logical ::                     &
       l_mphys=.True.            &   ! Switch to use microphysics
      ,l_advect=.True.           &   ! Switch to use advection
      ,l_diverge=.False.         &   ! Switch to use separate
                                     ! divergence code
      ,l_diverge_advection=.False.    ! Switch to add divergent 
                                     ! u velocity to advection calculation
  
  ! Switches to allow/prevent feedback onto background field  
  logical ::                         &
        l_fix_qv        =.False.     & ! Fix qv field
       ,l_nomphys_qv    =.False.     & ! Don't feedback from mphys on qv 
       ,l_noadv_qv      =.False.     & ! Don't feedback from advect on qv
       ,l_posadv_qv     =.False.       ! Only allow positive increments
                                       ! from advect to qv
  logical ::                         &
        l_fix_theta     =.False.     & ! Fix theta field
       ,l_nomphys_theta =.False.     & ! Don't feedback from mphys on theta 
       ,l_noadv_theta   =.False.       ! Don't feedback from advect on
                                       ! theta
  logical ::                         &
        l_noadv_hydrometeors=.False.   ! Don't advect hydrometeors
                                       ! NB This doesn't include
                                       ! sedimentation
  logical ::                         & 
        l_nodiv_hydrometeors=.False.   ! Don't apply divergence terms to
                                       ! hydrometeors
  logical ::                         &
        l_noadv_aerosols=.False.   ! Don't advect aerosols
                                       ! NB This doesn't include
                                       ! sedimentation
  logical ::                         & 
        l_nodiv_aerosols=.False.   ! Don't apply divergence terms to
                                       ! aerosols
  logical ::                         &
       l_fix_aerosols=.True.      ! Fix aerosol field.

  logical ::                      &
       l_pupdate=.False.   ! Switch to update pressure field

  ! Switches to go in mphys routines

  logical ::                      &
       l_sediment=.True. ! Allow sedimentation of
                         ! hydrometeors
  integer,parameter ::   &
       isurface_fixed=1    ! Fix surface values in advection
                           ! scheme (modification still allowed
                           ! through mphys)
  integer,parameter ::   &
       isurface_flux=2   ! Flux values for surface in advection
                         ! scheme 
  logical ::    &
       l_periodic_bound=.True. ! True for periodic boundaries (default)
                              ! False for fixed boundaries

  logical :: l_force_positive=.false. ! Switch to put on checking of negative number

                             
  ! Switches for hydrometeor initiation
  logical ::                      &
       l_hinit(nspecies)=.False.  & ! Switch to initialize hydrometeor field
       ,l_ainit(naerosol)=.False.   ! Switch to initialize aerosol field

  ! Integer switches for surface condition
  integer :: isurface=isurface_fixed 
  !integer :: isurface=isurface_flux

  ! Integer switches for microphysics schemes
  integer, parameter ::             &
        imphys_lem2_4           =1  & ! UKMO LEM version 2.4 bulk scheme
       ,imphys_tau_bin          =2  & ! UKMO LEM version 2.4 bin scheme
                                      ! (warm only)
       ,imphys_lem2_4_icebin    =3  & ! UKMO LEM version 3.0 3-phase 
                                      ! bin scheme
       ,imphys_um7_1            =4  & ! UKMO UM version 7.1
       ,imphys_thompson06       =5  & ! Greg Thompson's (2006) scheme
       ,imphys_morr_two_moment  =6  & ! Hugh Morrisons's 2M scheme
       ,imphys_thompson07       =7  & ! Greg Thompson's (2007) scheme
       ,imphys_thompson09       =8  & ! Greg Thompson's (2009) scheme
       ,imphys_um7_3            =9  & ! UM test: Jonathan Wilkinson
       ,imphys_wsm6             =10 & ! WSM6 scheme
       ,imphys_wdm6             =11 & ! WDM6 scheme
       ,imphys_4A               =12 & ! Shipway scheme
       ,imphys_mg1_5            =13 & ! Morrison Gettelman, v1_5  
       ,imphys_mg2              =14 & ! Morrison Gettelman, v2    
       ,imphys_mg3              =15   ! Morrison Gettelman, v3    
 
  ! Integer switches to choose input data type and/or test case
  integer, parameter ::    & 
        iukmo_lem     =1   &  ! From UKMO LEM
       ,iukmo_um      =2   &  ! From UKMO UM
       ,iecmwf        =3   &  ! From ECMWF
       ,iarm          =4      ! From ARM

       ! GCSS Warm rain cases
  integer, parameter ::    & 
        igcss_warm1   =101 &  
       ,igcss_warm2   =102 &  
       ,igcss_warm3   =103 &  
!       ,igcss_warm4   =104 &  !removed - now use 101 with different namelist
       ,igcss_warm5   =105 &  
       ,igcss_warm6   =106 &  
       ,igcss_warm7   =107 &  
       ,igcss_warm8   =108 &   
!++ag
       ,igcss_warmX   =109 &   
       ,igcss_warm1A  =110 &
       ,igcss_warm00  =198 &
       ,igcss_warm0   =199
!--ag
       ! GCSS mixed-phase cases
  integer, parameter ::    & 
        igcss_mixed1  =201 &  
!       ,igcss_mixed2  =202 &  !removed - now use 201 with different namelist
       ,igcss_mixed3  =203 &
       ,igcss_mixed4  =204
   
       ! GCSS Deep cases
  integer, parameter ::    & 
        igcss_deep1   =301 &  
       ,igcss_deep2   =302 &  
       ,igcss_deep3   =303
       ! Other cases 
  integer, parameter ::    & 
        iicel         =401 &  ! ICE-L-type case
!++ag
       ,ice1          =402 &  ! Cirrus ice case (higher)
!--ag
       ,itest_case    =998 &  ! Test case
       ,itest_case1   =999    ! Test case 2 
  integer, parameter :: &
        ipassive_sed_bulk1 = 501 & ! Passive sedimentation
       ,ipassive_sed_bulk2 = 502   ! Passive sedimentation
  integer, parameter :: &
       itest_2D = 600       &
       ,igcss_2d_Cu = 601   &  ! 2-D Cu case based on Morrison & Grabowski 06
       ,igcss_2d_Sc = 602   &  ! 2-D Sc case based on modified Morrison & Grabowski 06
       ,igcss_2d_ISDAC = 604 &  ! 2-D ISDAC case based on modified Morrison & Grabowski 06
       ,igcss_2d_squall = 605   ! 2-D Ideallised squall case

  integer, parameter :: iwmo_case1 = 650 ! 2-D Sc case for WMO workshop

  integer :: icase=igcss_warm1

  ! Various switches and their default values
  character(20) :: mphys_scheme='lem2.4'
  ! variant on default microphysics scheme
  ! will need documenting for each scheme
  integer :: mphys_var=0  
  
  character(100) :: input_file=''
  logical :: l_input_file=.True.
  integer :: ifiletype=iecmwf
  
  integer :: imphys=imphys_lem2_4

  ! Switches to be used for interfacing to python driver
  ! Otherwise, no need to change these.
  logical :: l_namelists=.true.
  logical :: l_write_dgs=.true.

  ! Switches to control test cases
  integer, parameter :: nctrl = 8
  real(wp), dimension(nctrl) :: &
         wctrl=0.           & ! control w values
       , zctrl=0.           & ! control z values
       , xctrl=0.           & ! control x values
       , tctrl=0.           & ! control time values
       , pctrl_z=0.         & ! control profile values for height
       , pctrl_T=0.         & ! control profile values for temperature
       , pctrl_v=0.         & ! control profile values for vapour
       , lhf_ctrl=0.        & ! control latent heat forcing
       , shf_ctrl=0.          ! control sensible heat forcing
  integer :: ipctrl=0 ! control standard profile choices
  
  integer :: iforce_method=0 ! Method of applying forcing
                             !  0 = tendencies applied
                             !  1 = relaxation applied

  
  logical :: l_cu_cold=.false. ! make the cumulus case 20C colder

  !---------------------------------------
  ! microphysics scheme specific switches
  !---------------------------------------

  ! Thompson09...
  logical :: l_reuse_thompson_lookup=.false. ! If true, reuse previously 
                                             ! calculated tables
  ! Shipway 4A...
  integer :: set_4A_option=-999 ! sets 4A scheme option
  logical :: l_evap=.true.   ! evaporate rain
  logical :: l_rsource=.true. ! rain sources
  logical :: l_raut=.true. ! rain sources
  logical :: l_cond=.true.   ! condensation
  logical :: l_bouss=.false. ! set rho=1 everywhere
  logical :: l_3msed=.false. ! set 3m sedimentation
  real(wp):: dt_sub=-999. ! 4A substep length (s)
  real(wp):: dt_sed=-999. ! 4A sedimentation substep length (s)
  integer :: diag_mu=-999 ! select diagnostic mu
  real(wp):: set_p1=-999. ! set p1
  real(wp):: set_sp1=-999. ! set sp1
  real(wp):: set_p2=-999. ! set p2
  real(wp):: set_sp2=-999. ! set sp2
  real(wp):: set_p3=-999. ! set p2
  real(wp):: set_sp3=-999. ! set sp2
  logical :: set_l_abelshipway=.false. ! use abelshipway fallspeed
  logical :: set_l_cons=.false. ! conserve number


end Module switches

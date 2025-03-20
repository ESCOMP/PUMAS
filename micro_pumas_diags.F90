module micro_pumas_diags

!----------------------------------------
! PUMAS diagnostics support package
!----------------------------------------

use pumas_kinds,   only: r8=>kind_r8

!> \section arg_table_proc_rates_type  Argument Table
!! \htmlinclude proc_rates_type.html
  type, public :: proc_rates_type

    real(r8), allocatable :: prodsnow(:,:)     ! production of snow (kg/kg/s)
    real(r8), allocatable :: evapsnow(:,:)     ! sublimation rate of snow (1/s)
    real(r8), allocatable :: qcsevap(:,:)      ! cloud water evaporation due to sedimentation (1/s)
    real(r8), allocatable :: qisevap(:,:)      ! cloud ice sublimation due to sublimation (1/s)
    real(r8), allocatable :: qvres(:,:)        ! residual condensation term to ensure RH < 100% (1/s)
    real(r8), allocatable :: cmeitot(:,:)      ! grid-mean cloud ice sub/dep (1/s)
    real(r8), allocatable :: vtrmc(:,:)        ! mass-weighted cloud water fallspeed (m/s)
    real(r8), allocatable :: vtrmi(:,:)        ! mass-weighted cloud ice fallspeed (m/s)
    real(r8), allocatable :: umr(:,:)          ! mass weighted rain fallspeed (m/s)
    real(r8), allocatable :: ums(:,:)          ! mass weighted snow fallspeed (m/s)
    real(r8), allocatable :: umg(:,:)          ! mass weighted graupel/hail fallspeed (m/s)
    real(r8), allocatable :: qgsedten(:,:)     ! qg sedimentation tendency (1/s)
    real(r8), allocatable :: qcsedten(:,:)     ! qc sedimentation tendency (1/s)
    real(r8), allocatable :: qisedten(:,:)     ! qi sedimentation tendency (1/s)
    real(r8), allocatable :: qrsedten(:,:)     ! qr sedimentation tendency (1/s)
    real(r8), allocatable :: qssedten(:,:)     ! qs sedimentation tendency (1/s)

    real(r8), allocatable :: pratot(:,:)
    real(r8), allocatable :: prctot(:,:)
    real(r8), allocatable :: mnuccctot(:,:)
    real(r8), allocatable :: mnuccttot(:,:)
    real(r8), allocatable :: msacwitot(:,:)
    real(r8), allocatable :: psacwstot(:,:)
    real(r8), allocatable :: bergstot(:,:)
    real(r8), allocatable :: vapdepstot(:,:)
    real(r8), allocatable :: bergtot(:,:)
    real(r8), allocatable :: melttot(:,:)
    real(r8), allocatable :: meltstot(:,:)
    real(r8), allocatable :: meltgtot(:,:)
    real(r8), allocatable :: homotot(:,:)
    real(r8), allocatable :: qcrestot(:,:)
    real(r8), allocatable :: prcitot(:,:)
    real(r8), allocatable :: praitot(:,:)
    real(r8), allocatable :: qirestot(:,:)
    real(r8), allocatable :: mnuccrtot(:,:)
    real(r8), allocatable :: mnudeptot(:,:)
    real(r8), allocatable :: mnuccritot(:,:)
    real(r8), allocatable :: pracstot(:,:)
    real(r8), allocatable :: meltsdttot(:,:)
    real(r8), allocatable :: frzrdttot(:,:)
    real(r8), allocatable :: mnuccdtot(:,:)
    real(r8), allocatable :: pracgtot(:,:)
    real(r8), allocatable :: psacwgtot(:,:)
    real(r8), allocatable :: pgsacwtot(:,:)
    real(r8), allocatable :: pgracstot(:,:)
    real(r8), allocatable :: prdgtot(:,:)
    real(r8), allocatable :: qmultgtot(:,:)
    real(r8), allocatable :: qmultrgtot(:,:)
    real(r8), allocatable :: psacrtot(:,:)
    real(r8), allocatable :: npracgtot(:,:)
    real(r8), allocatable :: nscngtot(:,:)
    real(r8), allocatable :: ngracstot(:,:)
    real(r8), allocatable :: nmultgtot(:,:)
    real(r8), allocatable :: nmultrgtot(:,:)
    real(r8), allocatable :: npsacwgtot(:,:)

  real(r8), allocatable :: nnuccctot(:,:)        ! change n  due to Immersion freezing of cloud water
  real(r8), allocatable :: nnuccttot(:,:)        ! change n  due to Contact freezing of cloud water
  real(r8), allocatable :: nnuccdtot(:,:)        ! change n  due to Ice nucleation
  real(r8), allocatable :: nnudeptot(:,:)        ! change n  due to Deposition Nucleation
  real(r8), allocatable :: nhomotot(:,:)         ! change n  due to Homogeneous freezing of cloud water
  real(r8), allocatable :: nnuccrtot(:,:)        ! change n  due to heterogeneous freezing of rain to snow (1/s)
  real(r8), allocatable :: nnuccritot(:,:)       ! change n  due to Heterogeneous freezing of rain to ice
  real(r8), allocatable :: nsacwitot(:,:)        ! change n  due to Conversion of cloud water [to cloud ice]
                                                         !                  from rime-splintering
  real(r8), allocatable :: npratot(:,:)          ! change n  due to Accretion of cloud water by rain
  real(r8), allocatable :: npsacwstot(:,:)       ! change n  due to Accretion of cloud water by snow
  real(r8), allocatable :: npraitot(:,:)         ! change n  due to Accretion of cloud ice to snow
  real(r8), allocatable :: npracstot(:,:)        ! change n  due to Accretion of rain by snow
  real(r8), allocatable :: nprctot(:,:)          ! change nr  due to Autoconversion of cloud water [to rain]
  real(r8), allocatable :: nprcitot(:,:)         ! change n  due to Autoconversion of cloud ice to snow
  real(r8), allocatable :: ncsedten(:,:)         ! change n  due to cloud liquid sedimentation
  real(r8), allocatable :: nisedten(:,:)         ! change n  due to cloud ice sedimentation
  real(r8), allocatable :: nrsedten(:,:)         ! change n  due to rain sedimentation
  real(r8), allocatable :: nssedten(:,:)         ! change n  due to snow sedimentation
  real(r8), allocatable :: ngsedten(:,:)         ! change n  due to graupel sedimentation
  real(r8), allocatable :: nmelttot(:,:)         ! change n  due to Melting of cloud ice
  real(r8), allocatable :: nmeltstot(:,:)        ! change n  due to Melting of snow
  real(r8), allocatable :: nmeltgtot(:,:)        ! change n  due to Melting of graupel

  ! TAU diagnostic variables
  real(r8), allocatable :: nraggtot(:,:)          ! change nr  due to self collection of rain


  real(r8), allocatable :: pgam_out(:,:)      ! Liquid Size distribution parameter Mu for output
  real(r8), allocatable :: lamc_out(:,:)      ! Liquid Size distribution parameter Lambda for output
  real(r8), allocatable :: lamr_out(:,:)      ! Rain Size distribution parameter Lambda for output
  real(r8), allocatable :: n0r_out(:,:)       ! Size distribution parameter n0 for output
  real(r8), allocatable :: scale_qc(:,:)        !TAU scaling factor for liquid mass to ensure conservation
  real(r8), allocatable :: scale_nc(:,:)        !TAU scaling factor for liquid number to ensure conservation
  real(r8), allocatable :: scale_qr(:,:)        !TAU scaling factor for rain mass to ensure conservation
  real(r8), allocatable :: scale_nr(:,:)        !TAU scaling factor for rain_number to ensure conservation
  real(r8), allocatable :: amk_c(:,:,:)         !TAU cloud liquid mass from bins
  real(r8), allocatable :: ank_c(:,:,:)         !TAU cloud liquid number from bins
  real(r8), allocatable :: amk_r(:,:,:)         !TAU cloud rain mass from bins
  real(r8), allocatable :: ank_r(:,:,:)         !TAU cloud rain number from bins
  real(r8), allocatable :: amk(:,:,:)           !TAU all liquid mass from bins
  real(r8), allocatable :: ank(:,:,:)           !TAU all liquid number from bins
  real(r8), allocatable :: amk_out(:,:,:)       !TAU all liquid number from bins output
  real(r8), allocatable :: ank_out(:,:,:)       !TAU all liquid mass from bins output
  real(r8), allocatable :: qc_out_TAU(:,:)      !TAU: output total cloud liquid mass
  real(r8), allocatable :: nc_out_TAU(:,:)      !TAU: output total cloud liquid number
  real(r8), allocatable :: qr_out_TAU(:,:)      !TAU: output total rain mass
  real(r8), allocatable :: nr_out_TAU(:,:)      !TAU: output total cloud rain number
  real(r8), allocatable :: qc_in_TAU(:,:)       !TAU: input total cloud liquid mass
  real(r8), allocatable :: nc_in_TAU(:,:)       !TAU: input total cloud liquid number
  real(r8), allocatable :: qr_in_TAU(:,:)       !TAU: input total rain mass
  real(r8), allocatable :: nr_in_TAU(:,:)       !TAU: input total cloud rain number
  real(r8), allocatable :: qctend_KK2000(:,:)   !cloud liquid mass tendency due to autoconversion  & accretion from KK2000
  real(r8), allocatable :: nctend_KK2000(:,:)   !cloud liquid number tendency due to autoconversion  & accretion from KK2000
  real(r8), allocatable :: qrtend_KK2000(:,:)   !rain mass tendency due to autoconversion  & accretion from KK2000
  real(r8), allocatable :: nrtend_KK2000(:,:)   !rain number tendency due to autoconversion   & accretion from KK2000
  real(r8), allocatable :: qctend_SB2001(:,:)   !cloud liquid mass tendency due to autoconversion  & accretion from SB2001
  real(r8), allocatable :: nctend_SB2001(:,:)   !cloud liquid number tendency due to autoconversion  & accretion from SB2001
  real(r8), allocatable :: qrtend_SB2001(:,:)   !rain mass tendency due to autoconversion  & accretion from SB2001
  real(r8), allocatable :: nrtend_SB2001(:,:)   !rain number tendency due to autoconversion  & accretion from SB2001
  real(r8), allocatable :: qctend_TAU(:,:)      !cloud liquid mass tendency due to autoconversion & accretion from TAU or Emulator code
  real(r8), allocatable :: nctend_TAU(:,:)      !cloud liquid number tendency due to autoconversion & accretion from TAU or Emulator code
  real(r8), allocatable :: qrtend_TAU(:,:)      !rain mass tendency due to autoconversion & accretion from TAU or Emulator code
  real(r8), allocatable :: nrtend_TAU(:,:)      !rain number tendency due to autoconversion & accretion from TAU or Emulatorcode
  real(r8), allocatable :: gmnnn_lmnnn_TAU(:,:) ! TAU sum of mass gain and loss from bin code
  real(r8), allocatable :: ML_fixer(:,:)     !Emulated: frequency of ML fixer is activated
  real(r8), allocatable :: QC_fixer(:,:)     !Emulated: change in cloud liquid mass due to ML fixer
  real(r8), allocatable :: NC_fixer(:,:)     !Emulated: change in cloud number number due to ML fixer
  real(r8), allocatable :: QR_fixer(:,:)     !Emulated: change in rain mass due to ML fixer
  real(r8), allocatable :: NR_fixer(:,:)     !Emulated: change in rain number due to ML fixer

    contains
      procedure :: allocate => proc_rates_allocate
      procedure :: deallocate => proc_rates_deallocate
  end type proc_rates_type

contains

   subroutine proc_rates_allocate(this, psetcols, nlev, ncd, warm_rain, errstring)
   !--------------------------------------------------------------
   ! Routine to allocate the elements of the proc_rates DDT
   !--------------------------------------------------------------

      implicit none

      class(proc_rates_type) :: this

      integer,           intent(in) :: psetcols, nlev
      integer,           intent(in) :: ncd
      character(len=16), intent(in) :: warm_rain            ! 'tau','emulated','sb2001' or 'kk2000'
      character(128),   intent(out) :: errstring

      integer :: ierr

      errstring=' '

      allocate(this%prodsnow(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%prodsnow'
      end if
      allocate(this%evapsnow(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%evapsnow'
      end if
      allocate(this%qcsevap(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qcsevap'
      end if
      allocate(this%qisevap(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qisevap'
      end if
      allocate(this%qvres(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qvres'
      end if
      allocate(this%cmeitot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%cmeitot'
      end if
      allocate(this%vtrmc(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%vtrmc'
      end if
      allocate(this%vtrmi(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%vtrmi'
      end if
      allocate(this%umr(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%umr'
      end if
      allocate(this%ums(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%ums'
      end if
      allocate(this%umg(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%umg'
      end if
      allocate(this%qgsedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qgsedten'
      end if
      allocate(this%qcsedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qcsedten'
      end if
      allocate(this%qisedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qisedten'
      end if
      allocate(this%qrsedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qrsedten'
      end if
      allocate(this%qssedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qssedten'
      end if
      allocate(this%pratot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pratot'
      end if
      allocate(this%prctot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%prctot'
      end if
      allocate(this%mnuccctot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccctot'
      end if
      allocate(this%mnuccttot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccttot'
      end if
      allocate(this%msacwitot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%msacwitot'
      end if
      allocate(this%psacwstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%psacwstot'
      end if
      allocate(this%bergstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%bergstot'
      end if
      allocate(this%vapdepstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%vapdepstot'
      end if
      allocate(this%bergtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%bergtot'
      end if
      allocate(this%melttot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%melttot'
      end if
      allocate(this%meltstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%meltstot'
      end if
      allocate(this%meltgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%meltgtot'
      end if
      allocate(this%homotot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%homotot'
      end if
      allocate(this%qcrestot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qcrestot'
      end if
      allocate(this%prcitot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%prcitot'
      end if
      allocate(this%praitot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%praitot'
      end if
      allocate(this%qirestot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qirestot'
      end if
      allocate(this%mnuccrtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccrtot'
      end if
      allocate(this%mnudeptot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnudeptot'
      end if
      allocate(this%mnuccritot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccritot'
      end if
      allocate(this%pracstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pracstot'
      end if
      allocate(this%meltsdttot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%meltsdttot'
      end if
      allocate(this%frzrdttot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%frzrdttot'
      end if
      allocate(this%mnuccdtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%mnuccdtot'
      end if
      allocate(this%pracgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pracgtot'
      end if
      allocate(this%psacwgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%psacwgtot'
      end if
      allocate(this%pgsacwtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pgsacwtot'
      end if
      allocate(this%pgracstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%pgracstot'
      end if
      allocate(this%prdgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%prdgtot'
      end if
      allocate(this%qmultgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qmultgtot'
      end if
      allocate(this%qmultrgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qmultrgtot'
      end if
      allocate(this%psacrtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%psacrtot'
      end if
      allocate(this%npracgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npracgtot'
      end if
      allocate(this%nscngtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nscngtot'
      end if
      allocate(this%ngracstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%ngracstot'
      end if
      allocate(this%nmultgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmultgtot'
      end if
      allocate(this%nmultrgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmultrgtot'
      end if
      allocate(this%npsacwgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npsacwgtot'
      end if
      allocate(this%nnuccctot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccctot'
      end if
      allocate(this%nnuccttot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccttot'
      end if
      allocate(this%nnuccdtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccdtot'
      end if
      allocate(this%nnudeptot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnudeptot'
      end if
      allocate(this%nhomotot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nhomotot'
      end if
      allocate(this%nnuccrtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccrtot'
      end if
      allocate(this%nnuccritot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nnuccritot'
      end if
      allocate(this%nsacwitot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nsacwitot'
      end if
      allocate(this%npratot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npratot'
      end if
      allocate(this%npsacwstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npsacwstot'
      end if
      allocate(this%npraitot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npraitot'
      end if
      allocate(this%npracstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%npracstot'
      end if
      allocate(this%nprctot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nprctot'
      end if
      allocate(this%nraggtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nraggtot'
      end if
      allocate(this%nprcitot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nprcitot'
      end if
      allocate(this%ncsedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%ncsedten'
      end if
      allocate(this%nisedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nisedten'
      end if
      allocate(this%nrsedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nrsedten'
      end if
      allocate(this%nssedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nssedten'
      end if
      allocate(this%ngsedten(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%ngsedten'
      end if
      allocate(this%nmelttot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmelttot'
      end if
      allocate(this%nmeltstot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmeltstot'
      end if
      allocate(this%nmeltgtot(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nmeltgtot'
      end if
      allocate(this%lamc_out(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
         errstring='Error allocating this%lamc_out'
      end if
      allocate(this%lamr_out(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
         errstring='Error allocating this%lamr_out'
      end if
      allocate(this%pgam_out(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
         errstring='Error allocating this%pgam_out'
      end if
      allocate(this%n0r_out(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
         errstring='Error allocating this%n0r_out'
      end if

      ! Only allocate these variables if machine learning turned on

      if (trim(warm_rain) == 'tau' .or. trim(warm_rain) == 'emulated') then
         allocate(this%scale_qc(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%scale_qc'
         end if
         allocate(this%scale_nc(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%scale_nc'
         end if
         allocate(this%scale_qr(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%scale_qr'
         end if
         allocate(this%scale_nr(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%scale_nr'
         end if
         allocate(this%amk_c(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%amk_c'
         end if
         allocate(this%ank_c(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ank_c'
         end if
         allocate(this%amk_r(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%amk_r'
         end if
         allocate(this%ank_r(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ank_r'
         end if
         allocate(this%amk(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%amk'
         end if
         allocate(this%ank(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ank'
         end if
         allocate(this%amk_out(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%amk_out'
         end if
         allocate(this%ank_out(psetcols,nlev,ncd), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ank_out'
         end if
         allocate(this%qc_out_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qc_out_TAU'
         end if
         allocate(this%nc_out_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nc_out_TAU'
         end if
         allocate(this%qr_out_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qr_out_TAU'
         end if
         allocate(this%nr_out_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nr_out_TAU'
         end if
         allocate(this%qc_in_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qc_in_TAU'
         end if
         allocate(this%nc_in_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nc_in_TAU'
         end if
         allocate(this%qr_in_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qr_in_TAU'
         end if
         allocate(this%nr_in_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nr_in_TAU'
         end if
         allocate(this%qctend_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qctend_TAU'
         end if
         allocate(this%nctend_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nctend_TAU'
         end if
         allocate(this%qrtend_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qrtend_TAU'
         end if
         allocate(this%nrtend_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nrtend_TAU'
         end if
         allocate(this%gmnnn_lmnnn_TAU(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%gmnnn_lmnnn_TAU'
         end if
         allocate(this%ML_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%ML_fixer'
         end if
         allocate(this%QC_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%QC_fixer'
         end if
         allocate(this%NC_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%NC_fixer'
         end if
         allocate(this%QR_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%QR_fixer'
         end if
         allocate(this%NR_fixer(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%NR_fixer'
         end if
      else if (warm_rain == 'sb2001') then
         ! Classic default (non-ML) microphysics
         allocate(this%qctend_SB2001(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%qctend_SB2001'
         end if
         allocate(this%nctend_SB2001(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nctend_SB2001'
         end if
         allocate(this%qrtend_SB2001(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%artend_SB2001'
         end if
         allocate(this%nrtend_SB2001(psetcols,nlev), stat=ierr)
         if (ierr /= 0) then
           errstring='Error allocating this%nrtend_SB2001'
         end if
      end if

      ! Variables which are needed by all code (Machine Learning and non-ML)
      allocate(this%qctend_KK2000(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%qctend_KK2000'
      end if
      allocate(this%nctend_KK2000(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nctend_KK2000'
      end if
      allocate(this%qrtend_KK2000(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%artend_KK2000'
      end if
      allocate(this%nrtend_KK2000(psetcols,nlev), stat=ierr)
      if (ierr /= 0) then
        errstring='Error allocating this%nrtend_KK2000'
      end if

   end subroutine proc_rates_allocate

   subroutine proc_rates_deallocate(this, warm_rain)
   !--------------------------------------------------------------
   ! Routine to deallocate the elements of the proc_rates DDT
   !--------------------------------------------------------------

      class(proc_rates_type) :: this
      character(len=16), intent(in) :: warm_rain            ! 'tau','emulated','sb2001' or 'kk2000'

      deallocate(this%prodsnow)
      deallocate(this%evapsnow)
      deallocate(this%qcsevap)
      deallocate(this%qisevap)
      deallocate(this%qvres)
      deallocate(this%cmeitot)
      deallocate(this%vtrmc)
      deallocate(this%vtrmi)
      deallocate(this%umr)
      deallocate(this%ums)
      deallocate(this%umg)
      deallocate(this%qgsedten)
      deallocate(this%qcsedten)
      deallocate(this%qisedten)
      deallocate(this%qrsedten)
      deallocate(this%qssedten)
      deallocate(this%pratot)
      deallocate(this%prctot)
      deallocate(this%mnuccctot)
      deallocate(this%mnuccttot)
      deallocate(this%msacwitot)
      deallocate(this%psacwstot)
      deallocate(this%bergstot)
      deallocate(this%vapdepstot)
      deallocate(this%bergtot)
      deallocate(this%melttot)
      deallocate(this%meltstot)
      deallocate(this%meltgtot)
      deallocate(this%homotot)
      deallocate(this%qcrestot)
      deallocate(this%prcitot)
      deallocate(this%praitot)
      deallocate(this%qirestot)
      deallocate(this%mnuccrtot)
      deallocate(this%mnudeptot)
      deallocate(this%mnuccritot)
      deallocate(this%pracstot)
      deallocate(this%meltsdttot)
      deallocate(this%frzrdttot)
      deallocate(this%mnuccdtot)
      deallocate(this%pracgtot)
      deallocate(this%psacwgtot)
      deallocate(this%pgsacwtot)
      deallocate(this%pgracstot)
      deallocate(this%prdgtot)
      deallocate(this%qmultgtot)
      deallocate(this%qmultrgtot)
      deallocate(this%psacrtot)
      deallocate(this%npracgtot)
      deallocate(this%nscngtot)
      deallocate(this%ngracstot)
      deallocate(this%nmultgtot)
      deallocate(this%nmultrgtot)
      deallocate(this%npsacwgtot)
      deallocate(this%nnuccctot)
      deallocate(this%nnuccttot)
      deallocate(this%nnuccdtot)
      deallocate(this%nnudeptot)
      deallocate(this%nhomotot)
      deallocate(this%nnuccrtot)
      deallocate(this%nnuccritot)
      deallocate(this%nsacwitot)
      deallocate(this%npratot)
      deallocate(this%npsacwstot)
      deallocate(this%npraitot)
      deallocate(this%npracstot)
      deallocate(this%nprctot)
      deallocate(this%nraggtot)
      deallocate(this%nprcitot)
      deallocate(this%ncsedten)
      deallocate(this%nisedten)
      deallocate(this%nrsedten)
      deallocate(this%nssedten)
      deallocate(this%ngsedten)
      deallocate(this%nmelttot)
      deallocate(this%nmeltstot)
      deallocate(this%nmeltgtot)

      deallocate(this%qctend_KK2000)
      deallocate(this%nctend_KK2000)
      deallocate(this%qrtend_KK2000)
      deallocate(this%nrtend_KK2000)

      deallocate(this%lamc_out)
      deallocate(this%lamr_out)
      deallocate(this%pgam_out)
      deallocate(this%n0r_out)

      if (trim(warm_rain) == 'tau' .or. trim(warm_rain) == 'emulated') then
         deallocate(this%scale_qc)
         deallocate(this%scale_nc)
         deallocate(this%scale_qr)
         deallocate(this%scale_nr)
         deallocate(this%amk_c)
         deallocate(this%ank_c)
         deallocate(this%amk_r)
         deallocate(this%ank_r)
         deallocate(this%amk)
         deallocate(this%ank)
         deallocate(this%amk_out)
         deallocate(this%ank_out)
         deallocate(this%qc_out_TAU)
         deallocate(this%nc_out_TAU)
         deallocate(this%qr_out_TAU)
         deallocate(this%nr_out_TAU)
         deallocate(this%qc_in_TAU)
         deallocate(this%nc_in_TAU)
         deallocate(this%qr_in_TAU)
         deallocate(this%nr_in_TAU)
         deallocate(this%qctend_TAU)
         deallocate(this%nctend_TAU)
         deallocate(this%qrtend_TAU)
         deallocate(this%nrtend_TAU)
         deallocate(this%gmnnn_lmnnn_TAU)
         deallocate(this%ML_fixer)
         deallocate(this%QC_fixer)
         deallocate(this%NC_fixer)
         deallocate(this%QR_fixer)
         deallocate(this%NR_fixer)
      else if (trim(warm_rain) == 'sb2001') then
         deallocate(this%qctend_SB2001)
         deallocate(this%nctend_SB2001)
         deallocate(this%qrtend_SB2001)
         deallocate(this%nrtend_SB2001)
      end if

   end subroutine proc_rates_deallocate

end module micro_pumas_diags

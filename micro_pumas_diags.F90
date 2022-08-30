module micro_pumas_diags

!----------------------------------------
! PUMAS diagnostics support package
!----------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8

  type, public :: proc_rates_type

    real(r8), allocatable :: prodsnow(:,:)     ! production of snow (1/s)
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
  real(r8), allocatable :: nprctot(:,:)          ! change n  due to Autoconversion of cloud water [to rain]
  real(r8), allocatable :: nprcitot(:,:)         ! change n  due to Autoconversion of cloud ice to snow
  real(r8), allocatable :: ncsedten(:,:)         ! change n  due to cloud liquid sedimentation
  real(r8), allocatable :: nisedten(:,:)         ! change n  due to cloud ice sedimentation
  real(r8), allocatable :: nrsedten(:,:)         ! change n  due to rain sedimentation
  real(r8), allocatable :: nssedten(:,:)         ! change n  due to snow sedimentation
  real(r8), allocatable :: ngsedten(:,:)         ! change n  due to graupel sedimentation
  real(r8), allocatable :: nmelttot(:,:)         ! change n  due to Melting of cloud ice
  real(r8), allocatable :: nmeltstot(:,:)        ! change n  due to Melting of snow
  real(r8), allocatable :: nmeltgtot(:,:)        ! change n  due to Melting of graupel

    contains
      procedure :: allocate => proc_rates_allocate
      procedure :: deallocate => proc_rates_deallocate
  end type proc_rates_type

contains

   subroutine proc_rates_allocate(this, psetcols, nlev)
   !--------------------------------------------------------------
   ! Routine to allocate the elements of the proc_rates DDT
   !--------------------------------------------------------------

      class(proc_rates_type) :: this

      integer, intent(in) :: psetcols, nlev

      integer :: ierr

!CAC -- NEED TO ADD ALLOCATE CHECKS

      allocate(this%prodsnow(psetcols,nlev), stat=ierr)
      allocate(this%evapsnow(psetcols,nlev), stat=ierr)
      allocate(this%qcsevap(psetcols,nlev), stat=ierr)
      allocate(this%qisevap(psetcols,nlev), stat=ierr)
      allocate(this%qvres(psetcols,nlev), stat=ierr)
      allocate(this%cmeitot(psetcols,nlev), stat=ierr)
      allocate(this%vtrmc(psetcols,nlev), stat=ierr)
      allocate(this%vtrmi(psetcols,nlev), stat=ierr)
      allocate(this%umr(psetcols,nlev), stat=ierr)
      allocate(this%ums(psetcols,nlev), stat=ierr)
      allocate(this%umg(psetcols,nlev), stat=ierr)
      allocate(this%qgsedten(psetcols,nlev), stat=ierr)
      allocate(this%qcsedten(psetcols,nlev), stat=ierr)
      allocate(this%qisedten(psetcols,nlev), stat=ierr)
      allocate(this%qrsedten(psetcols,nlev), stat=ierr)
      allocate(this%qssedten(psetcols,nlev), stat=ierr)
      allocate(this%pratot(psetcols,nlev), stat=ierr)
      allocate(this%prctot(psetcols,nlev), stat=ierr)
      allocate(this%mnuccctot(psetcols,nlev), stat=ierr)
      allocate(this%mnuccttot(psetcols,nlev), stat=ierr)
      allocate(this%msacwitot(psetcols,nlev), stat=ierr)
      allocate(this%psacwstot(psetcols,nlev), stat=ierr)
      allocate(this%bergstot(psetcols,nlev), stat=ierr)
      allocate(this%vapdepstot(psetcols,nlev), stat=ierr)
      allocate(this%bergtot(psetcols,nlev), stat=ierr)
      allocate(this%melttot(psetcols,nlev), stat=ierr)
      allocate(this%meltstot(psetcols,nlev), stat=ierr)
      allocate(this%meltgtot(psetcols,nlev), stat=ierr)
      allocate(this%homotot(psetcols,nlev), stat=ierr)
      allocate(this%qcrestot(psetcols,nlev), stat=ierr)
      allocate(this%prcitot(psetcols,nlev), stat=ierr)
      allocate(this%praitot(psetcols,nlev), stat=ierr)
      allocate(this%qirestot(psetcols,nlev), stat=ierr)
      allocate(this%mnuccrtot(psetcols,nlev), stat=ierr)
      allocate(this%mnudeptot(psetcols,nlev), stat=ierr)
      allocate(this%mnuccritot(psetcols,nlev), stat=ierr)
      allocate(this%pracstot(psetcols,nlev), stat=ierr)
      allocate(this%meltsdttot(psetcols,nlev), stat=ierr)
      allocate(this%frzrdttot(psetcols,nlev), stat=ierr)
      allocate(this%mnuccdtot(psetcols,nlev), stat=ierr)
      allocate(this%pracgtot(psetcols,nlev), stat=ierr)
      allocate(this%psacwgtot(psetcols,nlev), stat=ierr)
      allocate(this%pgsacwtot(psetcols,nlev), stat=ierr)
      allocate(this%pgracstot(psetcols,nlev), stat=ierr)
      allocate(this%prdgtot(psetcols,nlev), stat=ierr)
      allocate(this%qmultgtot(psetcols,nlev), stat=ierr)
      allocate(this%qmultrgtot(psetcols,nlev), stat=ierr)
      allocate(this%psacrtot(psetcols,nlev), stat=ierr)
      allocate(this%npracgtot(psetcols,nlev), stat=ierr)
      allocate(this%nscngtot(psetcols,nlev), stat=ierr)
      allocate(this%ngracstot(psetcols,nlev), stat=ierr)
      allocate(this%nmultgtot(psetcols,nlev), stat=ierr)
      allocate(this%nmultrgtot(psetcols,nlev), stat=ierr)
      allocate(this%npsacwgtot(psetcols,nlev), stat=ierr)

      allocate(this%nnuccctot(psetcols,nlev), stat=ierr)
      allocate(this%nnuccttot(psetcols,nlev), stat=ierr)
      allocate(this%nnuccdtot(psetcols,nlev), stat=ierr)
      allocate(this%nnudeptot(psetcols,nlev), stat=ierr)
      allocate(this%nhomotot(psetcols,nlev), stat=ierr)
      allocate(this%nnuccrtot(psetcols,nlev), stat=ierr)
      allocate(this%nnuccritot(psetcols,nlev), stat=ierr)
      allocate(this%nsacwitot(psetcols,nlev), stat=ierr)
      allocate(this%npratot(psetcols,nlev), stat=ierr)
      allocate(this%npsacwstot(psetcols,nlev), stat=ierr)
      allocate(this%npraitot(psetcols,nlev), stat=ierr)
      allocate(this%npracstot(psetcols,nlev), stat=ierr)
      allocate(this%nprctot(psetcols,nlev), stat=ierr)
      allocate(this%nprcitot(psetcols,nlev), stat=ierr)
      allocate(this%ncsedten(psetcols,nlev), stat=ierr)
      allocate(this%nisedten(psetcols,nlev), stat=ierr)
      allocate(this%nrsedten(psetcols,nlev), stat=ierr)
      allocate(this%nssedten(psetcols,nlev), stat=ierr)
      allocate(this%ngsedten(psetcols,nlev), stat=ierr)
      allocate(this%nmelttot(psetcols,nlev), stat=ierr)
      allocate(this%nmeltstot(psetcols,nlev), stat=ierr)
      allocate(this%nmeltgtot(psetcols,nlev), stat=ierr)
   end subroutine proc_rates_allocate

   subroutine proc_rates_deallocate(this)
   !--------------------------------------------------------------
   ! Routine to deallocate the elements of the proc_rates DDT
   !--------------------------------------------------------------

      class(proc_rates_type) :: this

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
      deallocate(this%nprcitot)
      deallocate(this%ncsedten)
      deallocate(this%nisedten)
      deallocate(this%nrsedten)
      deallocate(this%nssedten)
      deallocate(this%ngsedten)
      deallocate(this%nmelttot)
      deallocate(this%nmeltstot)
      deallocate(this%nmeltgtot)

   end subroutine proc_rates_deallocate

end module micro_pumas_diags

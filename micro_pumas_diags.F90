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

   end subroutine proc_rates_deallocate

end module micro_pumas_diags

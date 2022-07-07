module micro_pumas_diags

!----------------------------------------
! PUMAS diagnostics support package
!----------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8

  type, public :: proc_rates_type

    real(r8), allocatable :: prodsnow(:,:)
    real(r8), allocatable :: evapsnow(:,:)

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

   end subroutine proc_rates_allocate

   subroutine proc_rates_deallocate(this)
   !--------------------------------------------------------------
   ! Routine to deallocate the elements of the proc_rates DDT
   !--------------------------------------------------------------

      class(proc_rates_type) :: this

      deallocate(this%prodsnow)
      deallocate(this%evapsnow)

   end subroutine proc_rates_deallocate

end module micro_pumas_diags

!Define Fortran kinds for use specifically
!within PUMAS.  This allows PUMAS to control
!the precision in its internal routines without
!having to depend on a specific host model
!implemention.

module pumas_kinds

   implicit none
   private
   save

   integer, public, parameter :: kind_r8 = selected_real_kind(12) !8-byte real
   integer, public, parameter :: kind_i8 = selected_int_kind(18)  !8-byte integer

end module pumas_kinds

module tester
    use shr_kind_mod,   only: r8=>shr_kind_r8
    implicit none
    integer, parameter, public :: i8 = selected_int_kind(18)

contains
    subroutine write_test_values(filename, num_variables, nn_variables, batch_size)
        character(len=128), intent(in) :: filename
        integer(i8), intent(in) :: num_variables
        integer(i8), intent(in) :: batch_size
        real(r8), dimension(batch_size, num_variables), intent(in) :: nn_variables
        logical :: itexists 

        inquire(file=filename, exist=itexists)
        if (itexists) then
            open(12, file=filename, status="old", position="append", action="write")
                write(12, *) nn_variables
            close(12)
        else
            open(12, file=filename, status="new", position="append", action="write")
                write(12, *) nn_variables
            close(12)
        endif
    end subroutine write_test_values
end module tester

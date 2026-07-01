MODULE data_types
USE precision

IMPLICIT NONE

! linked list
type intll(na)
    integer(inp), len :: na
    integer(inp), dimension(na) :: list
    type(intll(na)), pointer :: next
end type intll

type reall(na)
    integer(inp), len :: na
    real(dp), dimension(na) :: list
    type(reall(na)), pointer :: next
end type reall

INTERFACE print_array
    MODULE PROCEDURE :: printa, printl!, printr
END INTERFACE

contains
! Print integer 2-D array.
SUBROUTINE printa(array)
    IMPLICIT NONE
    integer(inp), allocatable, intent(in) :: array(:,:)
    integer :: i
    do i = 1, size(array(:,1))
        print *, array(i,:)
    end do
    print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
END SUBROUTINE printa

! Print logical 2-D array.
SUBROUTINE printl(array)
    IMPLICIT NONE
    logical, allocatable, intent(in) :: array(:,:)
    integer :: i
    do i = 1, size(array(:,1))
        print *, array(i,:)
    end do
    print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
END SUBROUTINE printl

END MODULE data_types

! rsa_hash.f90
! This module is the hash function for the rsa analysis
module rsa_hash
    use logger

contains

! hash function to turn a flaxible size sorted ring list to a hash value,
! the function returns the bucket number.
pure function hash10(key, table_size) result(h)
    use iso_fortran_env, only: int64
    implicit none

    integer, intent(in) :: key(10)
    integer, intent(in) :: table_size
    integer :: h
    integer :: i
    integer(int64) :: x

    x = 1469598103934665603_int64

    do i = 1, 10
        x = ieor(x, int(key(i), int64))
        x = x * 1099511628211_int64
    end do

    h = int(modulo(x, int(table_size, int64))) + 1
end function hash10

end module rsa_hash

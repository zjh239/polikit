! rsa_hash.f90
! This module is the hash function for the rsa analysis
module rsa_hash
    use precision
    use logger
    use iso_fortran_env, only: int64

    integer :: table_size = 500
    integer :: load_bar = 500
    real(dp) :: load_factor = 0.8
    integer, allocatable :: hash_table(:)

contains

subroutine init_hash()
    if (.not. allocated(hash_table)) allocate(hash_table(table_size))
    hash_table = 0
end subroutine init_hash

! hash function to turn a flaxible size sorted ring list to a hash value,
! the function returns the bucket number.
pure integer(int64) function hash_ring(key)
    use iso_fortran_env, only: int64
    implicit none
    ! IN:
    integer, intent(in) :: key(:)

    integer :: h
    integer :: i
    integer(int64) :: x

    x = 1469598103934665603_int64
    h = size(key)

    do i = 1, h
        x = ieor(x, int(key(i), int64))
        x = x * 1099511628211_int64
    end do

    hash_ring = x

end function hash_ring

! This part is called when the loading ratio is over a value.
! A larger hash table is created and elements in the old table are rehashed and moved to the new table.
integer function rehash(tail, raw_hash)
    implicit none
    ! IN:
    integer(int64), intent(in) :: raw_hash(:)
    ! INOUT:
    integer, allocatable, intent(inout) :: tail(:)
    !
    integer, allocatable :: tmp(:)
    integer :: old_size
    integer :: i, j, k
    integer :: head, new_index

    old_size = table_size
    table_size = 2*old_size

    allocate(tmp(table_size))
    tmp = 0
    do i = 1, old_size
        head = hash_table(i)
        do while(head /= 0)
            new_index = int(modulo(raw_hash(head), int(table_size, int64))) + 1
            ! insert to new table
            j = tail(head)
            tail(head) = tmp(new_index)
            tmp(new_index) = head
            head = j

        end do
    end do

    deallocate(hash_table)
    call move_alloc(tmp, hash_table)

    rehash = table_size*load_factor

end function rehash

! From a raw ring, check if it already exist in the hash table
! if not, insert it
pure logical function check_rp_hash(b_index, ar, h_table, tail, sorted_list)
    implicit none
    ! IN:
    integer, intent(in) :: b_index
    integer, intent(in) :: ar(:)
    integer, intent(in) :: h_table(:), tail(:), sorted_list(:,:)
    !
    integer :: h_value, head
    integer :: i

    if (h_table(b_index) == 0) then
        check_rp_hash = .false.
        return
    else
        head = h_table(b_index)
        do while (head /= 0)
            if (same_ring(sorted_list(head, :), ar)) then
                check_rp_hash = .true.
                return
            else
                head = tail(head)
            end if
        end do
        check_rp_hash = .false.
        return
    end if

end function check_rp_hash

! check if the rings are exactly the same
pure logical function same_ring(a, b)
    implicit none

    integer, intent(in) :: a(:), b(:)
    same_ring = all(a == b)

end function same_ring

! list expand
subroutine r_list_expand(l1, l2, l3, l4, l5)
    implicit none
    ! inout:
    integer, allocatable, intent(inout) :: l1(:,:), l2(:,:), l3(:), l4(:)
    integer(int64), allocatable, intent(inout) :: l5(:)
    ! PRIV:
    integer, allocatable :: tmp_2d(:,:), tmp_1d(:)
    integer(int64), allocatable :: tmp_64(:)
    integer :: old_size, new_size, w

    w = size(l1(1,:))
    old_size = size(l3)
    new_size = old_size*2

    allocate(tmp_2d(new_size, w))
    tmp_2d(:old_size, :) = l1
    deallocate(l1)
    call move_alloc(tmp_2d, l1)

    allocate(tmp_2d(new_size, w))
    tmp_2d(:old_size, :) = l2
    deallocate(l2)
    call move_alloc(tmp_2d, l2)

    allocate(tmp_1d(new_size))
    tmp_1d(:old_size) = l3
    deallocate(l3)
    call move_alloc(tmp_1d, l3)

    allocate(tmp_1d(new_size))
    tmp_1d(:old_size) = l4
    deallocate(l4)
    call move_alloc(tmp_1d, l4)

    allocate(tmp_64(new_size))
    tmp_64(:old_size) = l5
    deallocate(l5)
    call move_alloc(tmp_64, l5)

    print '(a,i0,a)', ' Ring list expanded, space cost: ', 2.5*sizeof(l1)/1024, ' KB;'

end subroutine r_list_expand

subroutine clean_hash()
    if (allocated(hash_table)) deallocate(hash_table)
end subroutine clean_hash

end module rsa_hash

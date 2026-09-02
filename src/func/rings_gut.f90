! rings_gut.f90
! guttman rings
MODULE rings_gut
    use rings_v2
    IMPLICIT NONE

CONTAINS

SUBROUTINE rsa_gut(maxlvl)     ! Ring statistics analysis simple
    IMPLICIT NONE
    ! in:
    integer, intent(in) :: maxlvl  ! Max length, decided by ring size limit
    ! PRIVATE:
    integer :: atom     ! Center node index

    integer, allocatable, dimension(:,:) :: pathArray
    real(dp) :: mean_path_size

    ringlist_cap = 1000
    ringlist_size = 0

    crude_ring_num = 0

    ! time counters
    tneilist = 0.
    tfindring = 0.
    tcheckrepi = 0.
    tcheckpr = 0.
    taddring = 0.

    max_ring_size = maxlvl*2
    print *, info//' Max branch length:', maxlvl
    if(.not. allocated(path_list_size)) allocate(path_list_size(natom))

    allocate(ring_list_sorted(ringlist_cap, max_ring_size))
    ring_list_sorted = 0
    allocate(ring_list_raw(ringlist_cap, max_ring_size))
    ring_list_raw = 0
    allocate(ring_len_list(ringlist_cap))
    ring_len_list = 0
    allocate(next_in_bucket(ringlist_cap))
    next_in_bucket = 0
    allocate(raw_hash_value(ringlist_cap))
    raw_hash_value = 0_int64
    call init_hash()

    path_list_size = 0

    DO atom = 1, natom
!         print *, atom
        CALL create_path_list(atom, maxlvl, pathArray)

        CALL find_rings(pathArray)

        call progress_bar(atom, natom)
    END DO

    call print_ringno(ring_len_list(:ringlist_size))

    mean_path_size = sum(path_list_size(:))/size(path_list_size)
    print *, info//' Average path list array size is:', mean_path_size

!     do atom = 1, ringlist_size
!         print '(I0, *(2x,I0))', atom, ringList(atom)%element
!     end do

    print *, info//' Crude ring number found:', crude_ring_num

    call clean_rsa_v2()

END SUBROUTINE rsa_gut

! Guttman definition
subroutine mod_pr_gut(branch1, branch2, vis, pathlist)
    implicit none
    integer(inp), intent(in) :: pathlist(:,:)
    integer(inp), intent(in) :: branch1(:), branch2(:)
    logical, intent(inout) :: vis(:,:)
    integer(inp) :: k, m, a, i, j
!     logical :: isodd
    logical, allocatable :: mask1(:), mask2(:), tmp(:)

    k = size(branch1)
    m = size(vis(:,1))
    allocate(mask1(m), source = .true.)
    allocate(mask2(m), source = .true.)

    mask1 = pathlist(:,2) == branch1(2)
    mask2 = pathlist(:,2) == branch2(2)

    vis(trueloc(mask1), trueloc(mask2)) = .false.
    vis(trueloc(mask2), trueloc(mask1)) = .false.

    deallocate(mask1)
    deallocate(mask2)
end subroutine mod_pr_gut


END MODULE rings_gut

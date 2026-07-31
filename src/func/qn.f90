! qn.f90
! this module is about oxygen classification and Q_n analysis
module qn
    use precision
    use logger
    use data_input, only: o_type, natom, coord_data, ntype, type_name
    use neighbor_finder, only: neigh_list, print_hist, cn_by_type
    use params, only: ox_type_set
    use stdlib_array

    integer, allocatable :: q_n(:,:)
contains

! check the cation and report their oxygen CN distribution.
subroutine qn_analysis()
    implicit none
    ! IN:
    ! PRIV:
    integer :: i, j, k, id
    integer :: ox_type
    logical, allocatable :: bo(:)

    print *, info//' Performing Q_n analysis;'

    ox_type = ox_type_set
    if (.not. allocated(neigh_list)) stop error//' Performing Q_n analysis without building neighbor list first, quit.'
    if (ox_type /= o_type) print *, warn//' Caution: the manually set oxygen type is not the auto detected one;'

    ! note: q_n contains a total q_n and pair-wise q_n. The total one is stored at ntype=oxygen column
    if (.not. allocated(q_n)) allocate(q_n(natom, ntype))
    q_n = 0

    associate(n_neighbor => neigh_list%n_neighbor, &
                neighbors => neigh_list%neighbors, &
                ptype => coord_data%ptype)

    do j = 1, ntype
        if (j == ox_type) then
            bo = bo_check(ox_type, n_neighbor(:)-cn_by_type(:, j))
!             print *, debug//' Number of this cation:', count(ptype /= j)
        else
            bo = bo_check(ox_type, cn_by_type(:, j))
!             print *, debug//' Number of this cation:', count(ptype == j)
        end if
!         print *, debug//' Number of BO:', count(bo .eqv. .true.)
        do i = 1, natom
            if (j == ox_type) then
                if (ptype(i) == ox_type) cycle
            else
                if (ptype(i) /= j) cycle
            end if

            do k = 1, n_neighbor(i)
                id = neighbors(i,k)
                if (bo(id)) then
                    q_n(i, j) = q_n(i, j)+1
                end if
            end do
        end do
    end do

    end associate

    call print_qn(ox_type)

end subroutine qn_analysis

! check if the oxygen is bridging-oxygen
function bo_check(ox_type, list) result(bo)
    implicit none
    ! IN:
    integer, intent(in) :: ox_type
    integer, intent(in) :: list(:)
    ! OUT:
    logical, allocatable :: bo(:)
    ! PRIV:
    integer :: i

    if (size(list) /= natom) stop error//' CN count list has wrong size, quit.'
    if (.not. allocated(bo)) allocate(bo(natom))
    bo = .false.
    associate(ptype => coord_data%ptype)
    do i = 1, natom
        if (ptype(i) == ox_type .and. list(i) == 2) then
            bo(i) = .true.
        end if
    end do
    end associate
end function bo_check

! print qn data by cation type
subroutine print_qn(ox_type)
    implicit none
    ! IN:
    integer, intent(in) :: ox_type
    ! PRIV:
    integer :: i

    print *, info//' - Q_n distribution -'
    print *, info//' - - - - - - - - - - - - - - - - - - - - - -'
    print *, '|  N  |  0  |  1  |  2  |  3  |  4  |  5  |'

    do i = ntype, 1, -1
        if (i == ox_type) then
            call print_hist(q_n(trueloc(coord_data%ptype /= i), i) , ' | Cation |')
        else
            call print_hist(q_n(trueloc(coord_data%ptype == i), i) , ' |  '//type_name(i)//'  | ')
        end if
    end do
    print *, info//' - - - - - - - - - - - - - - - - - - - - - -'

end subroutine print_qn

subroutine clean_qn()
    if (allocated(q_n)) deallocate(q_n)
end subroutine clean_qn

end module

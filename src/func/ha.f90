! ha.f90
! Honeycutt-Anderson parameter calculation.
module ha
    use precision
    use neighbor_finder, only: neigh_list
    use data_input, only: coord_data, natom, o_type
    use stdlib_array
    use logger
    implicit none

contains

subroutine calculate_ha()
    implicit none
    integer :: id, i, check_id, totaln, m
    integer :: index2, index3
    integer, allocatable, dimension(:) :: ind23
    integer, allocatable, dimension(:) :: shared_list

    associate(n_neighbor => neigh_list%n_neighbor, neighbors => neigh_list%neighbors)

    totaln = sum(n_neighbor)
    allocate(ind23(totaln), source=0)

    m = 1
    do id = 1, natom
        do i = 1, n_neighbor(id)
            check_id = neighbors(id, i)
            if (id < check_id) then
                call common_elements(neighbors(id, :n_neighbor(id)), &
                    neighbors(check_id, :n_neighbor(check_id)), index2, shared_list)
                call bonds_num(shared_list, index3)
                ind23(m) = index2*10 + index3
                m = m+1
                if (allocated(shared_list)) deallocate(shared_list)
            end if
        end do
    end do
    end associate

    call print_ha(ind23)

end subroutine calculate_ha

subroutine common_elements(list1, list2, n, out_list)
    implicit none
    ! IN:
    integer, dimension(:), intent(in) :: list1, list2
    ! OUT:
    integer, intent(out) :: n
    integer, allocatable, dimension(:), intent(inout) :: out_list
    ! PRIVATE:
    integer :: n1, n2
    integer :: i, j
    integer, allocatable, dimension(:) :: tmp_list

    n1 = size(list1)
    n2 = size(list2)
    allocate(tmp_list(n1))

    i = 1
    j = 1
    n = 0

    do while (i <= n1 .and. j <= n2)
        if (list1(i) == list2(j)) then
            n = n + 1
            tmp_list(n) = list1(i)
            i = i + 1
            j = j + 1
        else if (list1(i) < list2(j)) then
            i = i + 1
        else
            j = j + 1
        end if
    end do

    allocate(out_list(n))
    out_list = tmp_list(:n)

!     print *, "Number of common neighbor:", n
end subroutine common_elements

subroutine bonds_num(idlist, n)
    implicit none
    ! IN:
    integer, dimension(:), intent(in) :: idlist
    ! OUT:
    integer, intent(out) :: n
    ! PRIVATE:
    integer :: i, j, k, id, check_id

    n = 0
    do i = 1, size(idlist)-1
        id = idlist(i)
        do j = i+1, size(idlist)
            check_id = idlist(j)
            k = neigh_list%n_neighbor(check_id)
            if (any(neigh_list%neighbors(check_id, :k) == id)) then
                n = n+1
            end if
        end do
    end do

end subroutine bonds_num

subroutine print_ha(list)
    implicit none
    ! IN:
    integer, dimension(:), intent(in) :: list
    !
    integer :: maxn, i
    integer, allocatable :: rank(:), amount(:)

    maxn = maxval(list)

    allocate(rank(0:maxn), amount(0:maxn))

    rank = [(i, i=0, maxn, 1)]
    do i = 0, maxn, 1
      amount(i) = count(list == i)
    end do

    print *, ' ### Honeycutt-Anderson Parameter'
    print *, '***************************'
    print 125, ' | 1**1 | ', pack(rank, amount(0:maxn) /= 0)
    print 125, ' | No.  | ', pack(amount(0:maxn), amount(0:maxn) /= 0)
    print *, '***************************'
125 format (a11,*(i10, ' | '))
end subroutine print_ha

subroutine calculate_qn()
    implicit none

    integer, allocatable :: qn(:)
    integer :: i, j, id

    if(.not. allocated(qn)) allocate(qn(natom))
    qn = 0

    do i = 1, natom
        if (coord_data%ptype(i) == o_type) continue
        do j = 1, neigh_list%n_neighbor(i)
            id  = neigh_list%neighbors(i,j)
            if (neigh_list%n_neighbor(id) == 2 .and. coord_data%ptype(id) == o_type) then
                qn(i) = qn(i) + 1
            end if
        end do
    end do

    j = maxval(qn(:))
    do i = 1, j
        print '(a,i0,a,i0)', ' '//info//' ',  count(qn==i), ' cations with Qn = ', i
    end do

end subroutine calculate_qn

end module ha

MODULE rings_v2
    USE precision
    use iso_fortran_env, only: int64
    USE data_input, only: natom, coord_data
    USE neighbor_finder, only: neigh_list
    USE stdlib_array
    use rsa_hash
    IMPLICIT NONE

    integer, allocatable, dimension(:,:) :: ring_list_sorted
    integer, allocatable, dimension(:,:) :: ring_list_raw

    integer, allocatable, dimension(:) :: ring_len_list
    integer, allocatable, dimension(:) :: next_in_bucket
    integer(int64), allocatable, dimension(:) :: raw_hash_value

    integer(int64), allocatable :: path_list_size(:)

    integer :: ringlist_cap, ringlist_size  ! Main ring list capacity and current size.

    real(dp) :: tstart, tcheck, tneilist, tfindring, tcheckrepi, tcheckpr, taddring
    integer :: crude_ring_num

    integer :: max_ring_size = 24

CONTAINS

SUBROUTINE rsa_v2(maxlvl)     ! Ring statistics analysis simple
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

!     maxlvl = 8
    print *, 'Max branch length:', maxlvl
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
    END DO

    call print_ringno(ring_len_list(:ringlist_size))

    mean_path_size = sum(path_list_size(:))/size(path_list_size)
    print *, 'Average path list array size is:', mean_path_size

!     do atom = 1, ringlist_size
!         print '(I0, *(2x,I0))', atom, ringList(atom)%element
!     end do

    print *, info//' Crude ring number found:', crude_ring_num

    call clean_rsa_v2()

END SUBROUTINE rsa_v2

! This subroutine creates the shortest paths list of a given center node(atom).
SUBROUTINE create_path_list(id_in, lvlim, pathArray)
    IMPLICIT NONE
    ! IN:
    integer, intent(in) :: id_in, lvlim
    ! OUT:
    integer, allocatable, dimension(:,:), intent(out) :: pathArray ! pathlist
    ! PRIVATE:
    integer, allocatable, dimension(:,:) :: tmp
    integer :: i, row, k !
    integer :: lvl  !

    call cpu_time(tstart)

    allocate(pathArray(1,1), source=id_in) !     lvl = 1

    associate(n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)
    do lvl= 1, lvlim
        if (n_neighbor(id_in)==0) exit
        ! level is the current distance to the center node.
        ! Expand the path length for one step
        allocate(tmp(size(pathArray(:,1)), lvl+1), source=0)
        tmp(:,:lvl) = pathArray
        deallocate(pathArray)
        call move_alloc(tmp, pathArray)

        row = 1   ! loop over the second last column
        do  ! j = 1, size(pathArray(:,lvl))
            ! end point of odd ring
            if (lvl > 1) then
                if (any(pathArray(:, lvl-1) == pathArray(row,lvl))) then
                    go to 102
                else if (pathArray(row,lvl) == 0) then
                    go to 102
                end if
            end if
            k = 1
            do ! k = 1, n_neighbor(pathArray(row,lvl))
                ! skip if the atom already exist in former level.
                if (lvl > 1) then
                    if (any(pathArray(:, lvl-1) == neighbor(pathArray(row,lvl), k))) then
                        go to 101
                    end if
                end if

                if (pathArray(row,lvl+1)==0) then
                    pathArray(row,lvl+1) = neighbor(pathArray(row,lvl), k)
                else
                    ! Expand the path list for new path
                    allocate(tmp(size(pathArray(:,1))+1, size(pathArray(1,:))))
                    tmp(:size(pathArray(:,1)),:) = pathArray
                    deallocate(pathArray)  ! deallocated?
                    call move_alloc(tmp, pathArray)

                    row = row+1
                    pathArray(row:,:) = eoshift(pathArray(row:,:), shift=-1)
                    pathArray(row,:) = pathArray(row-1,:)
                    pathArray(row,lvl+1) = neighbor(pathArray(row,lvl), k)
                    ! because moved row-1 to row, so here the neighbor list is still the same
                end if
101             if (k >= n_neighbor(pathArray(row,lvl))) exit
                k = k+1
            end do

102         if (row==size(pathArray(:,lvl))) exit
            row = row+1

        end do
!     call printa(pathArray)
    end do
    end associate

    path_list_size(id_in) = row

    call cpu_time(tcheck)
    tneilist = tneilist + (tcheck - tstart)
END SUBROUTINE create_path_list

! This subroutine finds all the possible rings around a center atom, given the
!   constructed shortest paths list. Push all the found rings to a data container.
SUBROUTINE find_rings(pathlist)
! natomring, ringatom, noprlist, numnopr, noprindex)
    IMPLICIT NONE
    ! IN:
    integer, dimension(:,:), intent(in) :: pathlist
    ! inOUT:
    ! PRIVATE:
!     type(ring), allocatable, dimension(:) :: ringlist
    integer :: mxrow
    logical, allocatable, dimension(:,:) :: vis
    logical, allocatable, dimension(:) :: bpoint
    integer, allocatable, dimension(:,:) :: vispl
!     integer(inp), dimension(:), allocatable :: smlst
    integer :: rma, rmb, id, q, t(1), lvl, n, mxlvl, l, j, id2, id0, n1, n2
    integer :: row, row_2
    real(dp) :: start, end

    call cpu_time(tstart)

    mxrow = size(pathlist(:,1))
    mxlvl = size(pathlist(1,:))
    allocate(bpoint(mxrow), source = .false.)

    id0 = pathlist(1,1)

    allocate(vis(mxrow,mxrow), source=.true.)
    forall (j = 1:mxrow)
        vis(j,j) = .false.
    end forall

    bpoint(1) = .true.
    rma = 1
    if (mxlvl > 1) id  = pathlist(1,2)
    do row = 1, mxrow
        if (mxlvl == 1) exit

        if (pathlist(row,2) /= id) then
            rmb = row-1
            vis(rma:rmb, rma:rmb) = .false.
            bpoint(row) = .true.

            id = pathlist(row,2)
            rma = row
        else if (row==mxrow) then
            vis(rma:row, rma:row) = .false.
        end if
    end do

    do lvl = 3, mxlvl
        id  = pathlist(1,lvl)
        do row = 1, mxrow
            if (pathlist(row, lvl) == 0) then
                vis(row,:) = .false.
                vis(:,row) = .false.
                cycle
            end if
            if (bpoint(row) &   !
                .or. pathlist(row, lvl) /= id) then !
                bpoint(row) = .true.      !
                id = pathlist(row,lvl)

                id2 = id
                do row_2= row, mxrow

                    if (bpoint(row_2) .or. pathlist(row_2, lvl) /= id2) then

                        if (vis(row, row_2)) then
                        id2 = pathlist(row_2,lvl)

                            ! check for odd ring
                            if   (pathlist(row_2, lvl-1)==pathlist(row, lvl)&
                            .and. pathlist(row_2, lvl)==pathlist(row, lvl-1)) then
                                call cpu_time(tcheck)
                                tfindring = tfindring + (tcheck - tstart)
                                crude_ring_num = crude_ring_num +1

                                call add_ring(pathlist(row,:lvl), pathlist(row_2,:lvl))

                                call cpu_time(tstart)
                                vis(row_2,:) = .false.
                                vis(:,row_2) = .false.
                                vis(row,:) = .false.
                                vis(:,row) = .false.
                            end if

                            ! check for even rings
                            if (pathlist(row_2, lvl)==pathlist(row, lvl)) then
                                call cpu_time(tcheck)
                                tfindring = tfindring + (tcheck - tstart)
                                crude_ring_num = crude_ring_num +1

                                call add_ring(pathlist(row,:lvl), pathlist(row_2,:lvl))

                                call mod_pr(pathlist(row,:lvl), pathlist(row_2,:lvl), vis, pathlist)
                                call cpu_time(tstart)
                            end if
                        end if
                    end if
                end do
            end if
    !     call printl(vb)
        end do

        if (all(vis .eqv. .false.)) then
    !         print *, 'Visibility array full of FALSE at level ', lvl,', quit.'
            exit
        end if
    end do

    deallocate(bpoint)

    call cpu_time(tcheck)
    tfindring = tfindring + (tcheck - tstart)

END SUBROUTINE find_rings

FUNCTION checkShortCut(rr) RESULT(ifpr)
    IMPLICIT NONE
    ! IN:
    integer, intent(in) :: rr(:)
    ! OUT:
    logical :: ifpr
    ! PRIVATE:
    integer :: src1, src2, src3   !
    logical :: isodd
    integer, allocatable, dimension(:) :: elem  !
    integer, allocatable, dimension(:) :: head1, head2  !
    integer, allocatable, dimension(:) :: last1, last2  !
    integer, allocatable, dimension(:) :: scndlast1, scndlast2
    integer, allocatable, dimension(:) :: tmp
    integer :: lvl, j, n, m, l, distance
    integer :: brlen, clen, mxlvl  !
    
    call cpu_time(tstart)
    l = size(rr)

    isodd = .false.
    if (mod(l,2)/=0) isodd = .true.

    ifpr = .true.
    allocate(elem(l*2))
    elem = [rr(:), rr(:)]   !

    brlen = ceiling((l+1)/2.)
    mxlvl = ceiling(brlen/2.)

    associate(n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)
    do m = 1, brlen !-1

        distance = 0

        src1 = elem(m)
        src2 = elem(m+brlen-1)
        if (isodd) src3 = elem(m+brlen)

        allocate(last1(1), source = 0)
        allocate(last2(1), source = 0)

        allocate(head1(1), source = src1)
        if (isodd) then
            allocate(head2(2), source = [src2, src3])
        else
            allocate(head2(1), source = src2)
        endif

        do lvl = 2, mxlvl
            call move_alloc(last1, scndlast1)
            call move_alloc(head1, last1)
            allocate(head1(1), source=0)

            n = 1
            do while (n <= size(last1))
                do j = 1, n_neighbor(last1(n))
                    if (any(scndlast1==neighbor(last1(n), j))) cycle
                    if (any(last1==neighbor(last1(n), j))) cycle
                    if (any(head1==neighbor(last1(n), j))) cycle

                    if (head1(1)==0) then
                        head1(1) = neighbor(last1(n), j)
                    else
                        allocate(tmp(size(head1)+1))
                        tmp(:size(head1)) = head1
                        tmp(size(head1)+1) = neighbor(last1(n), j)
                        call move_alloc(tmp, head1)
                    end if
                end do

                n = n+1
            end do
            do n = 1, size(head1)
                if (any(head2 == head1(n))) then
                    ifpr = .false.
                    return
                end if
            end do
            clen = lvl*2 - 1
            if (clen == brlen) then
                if (m == brlen) then
                    ifpr = .true.
                    return
                else
                    cycle
                end if
            end if

            call move_alloc(last2, scndlast2)
            call move_alloc(head2, last2)
            allocate(head2(1), source=0)
            n = 1
            do while(n <= size(last2))
                do j = 1, n_neighbor(last2(n))
                    if (any(scndlast2==neighbor(last2(n), j))) cycle
                    if (any(last2==neighbor(last2(n), j))) cycle
                    if (any(head2==neighbor(last2(n), j))) cycle

                    if (head2(1)==0) then
                        head2(1) = neighbor(last2(n), j)
                    else
                        allocate(tmp(size(head2)+1))
                        tmp(:size(head2)) = head2
                        tmp(size(head2)+1) = neighbor(last2(n), j)
                        call move_alloc(tmp, head2)
                    end if
                end do
                n = n+1
            end do
            clen = lvl*2
            do n = 1, size(head2)
                if (any(head1 == head2(n))) then
                    ifpr = .false.
                    return
                end if
            end do
            if (allocated(scndlast1)) deallocate(scndlast1)
            if (allocated(scndlast2)) deallocate(scndlast2)
        end do

    if (allocated(head1)) deallocate(head1)
    if (allocated(head2)) deallocate(head2)

    if (allocated(last1)) deallocate(last1)
    if (allocated(last2)) deallocate(last2)

    end do
    end associate
    deallocate(elem)

    call cpu_time(tcheck)
    tcheckpr = tcheckpr + (tcheck - tstart)
END FUNCTION checkShortCut

! This subroutine adds a new ring type element to the ring list. The inputs are two
!   integer type lists, means the two branches of a ring.
SUBROUTINE add_ring(branch1, branch2)
    IMPLICIT NONE
    ! IN:
    integer, intent(in) :: branch1(:), branch2(:)
    ! INOUT:
    ! PRIVATE:
    integer, allocatable, dimension(:,:) :: tmp_2d
    integer, allocatable, dimension(:) :: tmp_1d, element, sorted

    logical :: isodd, doexist, ispr, gofound
    integer :: k, i, t, l
    integer :: rpos, b_index
    integer(int64) :: h_value

    k = size(branch1)
    
    ! need allocation
    l = 0
    allocate(element(max_ring_size))
    allocate(sorted(max_ring_size))
    element = 0
    sorted = 0

    if (branch1(k) == branch2(k)) then
        isodd=.false.
    else
        isodd=.true.
    end if

    if (isodd) then
        t = 2*k-3
        l = t
        element(:t) = [branch1(:k-2), branch2(k:2:-1)]
    else
         t = 2*k-2
         l = t
        element(:t) = [branch1(:k-1), branch2(k:2:-1)]
    end if

    ! Sort the ring elements and check if it already exist in the list.
    sorted(:t) = element(:t)
    call bubble_sort(t, sorted)

    call cpu_time(tstart)

    ! Use sorted ring to calculate hash.
    h_value = hash_ring(sorted(:t))

    b_index = int(modulo(h_value, int(table_size, int64))) + 1

    if (check_rp_hash(b_index, sorted, hash_table, next_in_bucket, ring_list_sorted)) then
        return
    else
        if (checkShortCut(element(:t))) then
            ringlist_size = ringlist_size+1
            ring_list_raw(ringlist_size, :) = element
            ring_list_sorted(ringlist_size, :) = sorted
            ring_len_list(ringlist_size) = t
            raw_hash_value(ringlist_size) = h_value

            next_in_bucket(ringlist_size) = hash_table(b_index)
            hash_table(b_index) = ringlist_size
        end if
    end if

    if (ringlist_size >= load_bar) load_bar = rehash(next_in_bucket, raw_hash_value)
    if (ringlist_size == ringlist_cap) then
        call r_list_expand(ring_list_raw, &
                        ring_list_sorted, next_in_bucket, ring_len_list, raw_hash_value)
        ringlist_cap = ringlist_cap*2
    end if

    call cpu_time(tcheck)
    tcheckrepi = tcheckrepi + (tcheck - tstart)


END SUBROUTINE add_ring

! Modify the visibility array according to the primitive ring definition.
SUBROUTINE mod_pr(branch1, branch2, vis, pathlist)
    IMPLICIT NONE
    integer, intent(in) :: pathlist(:,:)
    integer, intent(in) :: branch1(:), branch2(:)
    logical, intent(inout) :: vis(:,:)
    integer :: k, m, a, i, j
    logical :: isodd
    logical, allocatable :: mask1(:), mask2(:), tmp(:)

    call cpu_time(tstart)

    k = size(branch1)
    m = size(vis(:,1))
    allocate(mask1(m), source = .true.)
    allocate(mask2(m), source = .true.)
    allocate(tmp(m), source = .true.)

    if (branch1(k) == branch2(k)) then
        isodd=.false.
    else
        isodd=.true.
    end if

    if (isodd) then
        do i = 2, k-1
            mask1 = .true.
            mask2 = .true.

            do j = 2,i
                tmp = pathlist(:,j)==branch1(j)
                mask1 = mask1 .and. tmp
            end do
            a = k+1-i
    !         branch2(:a)
            do j = 2, a
                tmp = pathlist(:,j) == branch2(j)
                mask2 = mask2 .and. tmp
            end do
            vis(trueloc(mask1), trueloc(mask2)) = .false.
            vis(trueloc(mask2), trueloc(mask1)) = .false.
        end do
    else !if (isodd .eqv. .false.) then
    ! even ring
        do i = 2, k
            mask1 = .true.
            mask2 = .true.

            do j = 2,i
                tmp = pathlist(:,j) == branch1(j)
                mask1 = mask1 .and. tmp
            end do
            a = k+2-i

            do j = 2, a
                tmp = pathlist(:,j) == branch2(j)
                mask2 = mask2 .and. tmp
            end do
            vis(trueloc(mask1), trueloc(mask2)) = .false.
            vis(trueloc(mask2), trueloc(mask1)) = .false.
        end do
    end if

    deallocate(mask1)
    deallocate(mask2)
    deallocate(tmp)

    call cpu_time(tcheck)
    taddring = taddring + (tcheck - tstart)
END SUBROUTINE mod_pr

subroutine print_ringno(ring_l)
    implicit none
    ! IN:
    integer, dimension(:), intent(in) :: ring_l
    !
    integer :: maxn, i
    integer, allocatable :: rank(:), amount(:)

    maxn = maxval(ring_l)
    print *, 'Maximum ring length is: ', maxn
!     if (maxn > 20) maxn = 20

    allocate(rank(0:maxn), amount(0:maxn))

    rank = [(i, i=0, maxn, 1)]
    do i = 0, maxn, 1
      amount(i) = count(ring_l == i)
    end do

    print *, ' ### RSA Size Distribution'
    print *, '***************************'
    print 107, ' | Size  | ',rank
    print 107, ' | Count | ',amount(0:maxn)
    print *, '***************************'
107 format (a11,*(i6, ' | '))

    print *, ' ### RSA Time Cost'
    print *, '_________________________________________________________________________________'
    print *, '| T(Path List)  | T(Find Ring)  | T(Chk&Insert) | T(PR Check)   | T(VA modif)   |'
    print 108, tneilist, tfindring, tcheckrepi, tcheckpr, taddring
    print *, '|_______________|_______________|_______________|_______________|_______________|'
108 format (' | ', *(f11.3,' s | '))

end subroutine print_ringno

subroutine clean_rsa_v2()

    if (allocated(ring_list_sorted))    deallocate(ring_list_sorted)
    if (allocated(ring_list_raw))       deallocate(ring_list_raw)

    if (allocated(ring_len_list))       deallocate(ring_len_list)
    if (allocated(next_in_bucket))      deallocate(next_in_bucket)
    if (allocated(raw_hash_value))      deallocate(raw_hash_value)

    call clean_hash()
end subroutine clean_rsa_v2

END MODULE rings_v2

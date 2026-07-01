MODULE rings_simple
    USE precision
    USE data_input, only: natom, coord_data
    USE neighbor_finder, only: neigh_list
    USE stdlib_array
    IMPLICIT NONE

    TYPE ring
        integer :: l
        integer :: element(24) = 0
        integer :: sorted(24) = 0
    END TYPE

    integer, allocatable :: path_list_size(:)

    integer :: ringlist_cap, ringlist_size  ! Main ring list capacity and current size.
    real(dp) :: tstart, tcheck, tneilist, tfindring, tcheckrepi, tcheckpr, taddring
    integer :: crude_ring_num

    type(ring) :: emptyring

CONTAINS

SUBROUTINE rsa_simple(maxlvl)     ! Ring statistics analysis simple
    IMPLICIT NONE
    ! in:
    integer, intent(in) :: maxlvl  ! Max length, decided by ring size limit
    ! PRIVATE:
    integer :: atom     ! Center node index

    type(ring), allocatable, dimension(:) :: ringList
    integer(inp), allocatable, dimension(:,:) :: pathArray

    real(dp) :: mean_path_size

    emptyring%l = 0
    emptyring%element = 0
    emptyring%sorted = 0

    ringlist_cap = 1000
    ringlist_size = 0

    crude_ring_num = 0

    tneilist = 0.
    tfindring = 0.
    tcheckrepi = 0.
    tcheckpr = 0.
    taddring = 0.

!     maxlvl = 8
    print *, 'Max branch length:', maxlvl
    if(.not. allocated(path_list_size)) allocate(path_list_size(natom))

    allocate(ringList(ringlist_cap))
    ringList%l = 0

    path_list_size = 0

    DO atom = 1, natom
!         print *, atom
        CALL create_path_list(atom, maxlvl, pathArray)

        CALL find_rings(pathArray, ringList)
    END DO

    call print_ringno(ringList(:ringlist_size)%l)

    mean_path_size = sum(path_list_size(:))/size(path_list_size)
    print *, 'Average path list array size is:', mean_path_size

!     do atom = 1, ringlist_size
!         print '(I0, *(2x,I0))', atom, ringList(atom)%element
!     end do

    print *, 'Crude ring number found:', crude_ring_num

END SUBROUTINE rsa_simple

! This subroutine creates the shortest paths list of a given center node(atom).
SUBROUTINE create_path_list(id_in, lvlim, pathArray)
    IMPLICIT NONE
    ! IN:
    integer(inp), intent(in) :: id_in, lvlim
    ! OUT:
    integer(inp), allocatable, dimension(:,:), intent(out) :: pathArray ! pathlist
    ! PRIVATE:
    integer(inp), allocatable, dimension(:,:) :: tmp
    integer(inp) :: i, row, k !
    integer(inp) :: lvl  !

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
SUBROUTINE find_rings(pathlist, mainringlist)
! natomring, ringatom, noprlist, numnopr, noprindex)
    IMPLICIT NONE
    ! IN:
    integer(inp), allocatable, dimension(:,:), intent(in) :: pathlist
    ! inOUT:
    type(ring), allocatable, dimension(:), intent(inout) :: mainringlist
    ! PRIVATE:
    type(ring), allocatable, dimension(:) :: ringlist
    integer(inp) :: mxrow
    logical, allocatable, dimension(:,:) :: vis
    logical, allocatable, dimension(:) :: bpoint
    integer(inp), allocatable, dimension(:,:) :: vispl
!     integer(inp), dimension(:), allocatable :: smlst
    integer(inp) :: rma, rmb, id, q, t(1), lvl, n, mxlvl, l, j, id2, id0, n1, n2
    integer(inp) :: row, row_2
    real(dp) :: start, end

    call cpu_time(tstart)

    allocate(ringList(1))
    ringList(1)%l = 0
    ringList(1)%element = 0

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
    ! This part garantees the two paths split at the beginning.
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
            if (bpoint(row) &   ! This means the second last element changed.
                .or. pathlist(row, lvl) /= id) then ! This means the current element changed.
                bpoint(row) = .true.      ! Modify this for the next lvl check.
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

                                call add_ring(pathlist(row,:lvl), pathlist(row_2,:lvl), mainringlist)

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

                                call add_ring(pathlist(row,:lvl), pathlist(row_2,:lvl), mainringlist)

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

!     call rm_not_pr(ringList)

!     if (ringList(1)%l /= 0) call add_ringlist(mainringlist, ringList)

    deallocate(bpoint)

    call cpu_time(tcheck)
    tfindring = tfindring + (tcheck - tstart)

END SUBROUTINE find_rings

FUNCTION checkShortCut(rr) RESULT(ifpr)
    IMPLICIT NONE
    ! IN:
    type(ring), intent(in) :: rr    ! Input ring
    ! OUT:
    logical :: ifpr
    ! PRIVATE:
    integer :: src1, src2, src3   ! check-node and its mid-node(even ring) or mid-nodes(odd ring).
    logical :: isodd
    integer, allocatable, dimension(:) :: elem  ! Stores two repeated ring list elements for a constant offset.
    integer, allocatable, dimension(:) :: head1, head2  ! Heads list of the wave.
    integer, allocatable, dimension(:) :: last1, last2  ! Heads list of last level.
    integer, allocatable, dimension(:) :: scndlast1, scndlast2
    integer, allocatable, dimension(:) :: tmp
    integer :: lvl, j, n, m, l, distance
    integer :: brlen, clen, mxlvl  ! branch length, current length
    
    call cpu_time(tstart)

    isodd = .false.
    if (mod(rr%l,2)/=0) isodd = .true.

    ifpr = .true.
    allocate(elem(rr%l*2))
    elem = [rr%element(:rr%l), rr%element(:rr%l)]   ! Avoid seg. fault.

    brlen = ceiling((rr%l+1)/2.)
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
            ! update head1.
            call move_alloc(last1, scndlast1)
            call move_alloc(head1, last1)
            allocate(head1(1), source=0)

            n = 1
            do while (n <= size(last1))
                do j = 1, n_neighbor(last1(n))
                    ! if this atom is already in the wave, cycle
                    if (any(scndlast1==neighbor(last1(n), j))) cycle
                    if (any(last1==neighbor(last1(n), j))) cycle
                    if (any(head1==neighbor(last1(n), j))) cycle

                    ! if not cycled, push it to the list
                    if (head1(1)==0) then
                        head1(1) = neighbor(last1(n), j)
                    else
                        allocate(tmp(size(head1)+1))
                        tmp(:size(head1)) = head1
                        tmp(size(head1)+1) = neighbor(last1(n), j)
!                         deallocate(head1)
                        call move_alloc(tmp, head1)
                    end if
                end do

!                 if (n >= size(last1)) exit
                n = n+1
            end do
            ! check if the branches meet.
            do n = 1, size(head1)
                if (any(head2 == head1(n))) then
                    ! print *, 'Short cut found, length is', 2*i-2, ', meet at ', tmp1(n)
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
            ! update head2.
            call move_alloc(last2, scndlast2)
            call move_alloc(head2, last2)
            allocate(head2(1), source=0)
            n = 1
            do while(n <= size(last2))
                do j = 1, n_neighbor(last2(n))
                    ! if this atom is already in the wave, cycle
                    if (any(scndlast2==neighbor(last2(n), j))) cycle
                    if (any(last2==neighbor(last2(n), j))) cycle
                    if (any(head2==neighbor(last2(n), j))) cycle

                    ! if not cycled, push it to the list
                    if (head2(1)==0) then
                        head2(1) = neighbor(last2(n), j)
                    else
                        allocate(tmp(size(head2)+1))
                        tmp(:size(head2)) = head2
                        tmp(size(head2)+1) = neighbor(last2(n), j)
!                         deallocate(head2)
                        call move_alloc(tmp, head2)
                    end if
                end do

!                 if (n >= size(last2)) exit
                n = n+1
            end do
            clen = lvl*2
            ! check if the branches meet.
            do n = 1, size(head2)
                if (any(head1 == head2(n))) then
                    ifpr = .false.
                    return
                end if
            end do
            ! deallocate second last lists.
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
SUBROUTINE add_ring(branch1, branch2, mainringlist)
    IMPLICIT NONE
    ! IN:
    integer(inp), intent(in) :: branch1(:), branch2(:)
    ! INOUT:
    type(ring), allocatable, intent(inout) :: mainringlist(:)
    ! PRIVATE:
    type(ring), allocatable :: tmp(:)
    type(ring) :: ar
    logical :: isodd, doexist, ispr, gofound
    integer(inp) :: k, i, t
    integer :: rpos

    k = size(branch1)
    
    ar%l = 0
    ar%element = 0
    ar%sorted = 0

    if (branch1(k) == branch2(k)) then
        isodd=.false.
    else
        isodd=.true.
    end if

    if (isodd) then
        t = 2*k-3
        ar%l = t
        ar%element(:t) = [branch1(:k-2), branch2(k:2:-1)]
    else
         t = 2*k-2
         ar%l = t
        ar%element(:t) = [branch1(:k-1), branch2(k:2:-1)]
    end if

    ar%sorted(:t) = ar%element(:t)
    call bubble_sort(t, ar%sorted)

    ! Sort the ring elements and check if it already exist in the list.
    call cpu_time(tstart)
    call new_check_rp(ar, mainringlist, rpos, doexist)


    call cpu_time(tcheck)
    tcheckrepi = tcheckrepi + (tcheck - tstart)

    if (doexist .eqv. .false.) then
        ispr = checkShortCut(ar)
        if (ispr) then
            call cpu_time(tstart)

            if (ringlist_size == 0) then
                mainringlist(1) = ar
                ringlist_size = 1
            else if (ringlist_size < ringlist_cap) then
                mainringlist(rpos:ringlist_size+1) = eoshift(mainringlist(rpos:ringlist_size+1), shift=-1, boundary=emptyring)
                mainringlist(rpos) = ar
                ringlist_size = ringlist_size + 1

            else if (ringlist_size == ringlist_cap) then
                ! Expand the ring list if full.
                ringlist_cap = ringlist_cap*2
                allocate(tmp(ringlist_cap))
                tmp(:ringlist_size) = mainringlist(:ringlist_size)
                deallocate(mainringlist)
                call move_alloc(tmp, mainringlist)
                print '(a,i0,a)', ' Main ring list expanded to size: ', sizeof(mainringlist)/1024, ' KB;'
            end if

            call cpu_time(tcheck)
            taddring = taddring + (tcheck - tstart)
        end if
    end if

END SUBROUTINE add_ring

subroutine new_check_rp(ar, mainringlist, pos, goal_found)
    implicit none
    ! IN:
    type(ring), intent(in) :: ar
    type(ring), allocatable, intent(in) :: mainringlist(:)
    ! OUT:
    integer, intent(out) :: pos
    logical, intent(out) :: goal_found
    ! Private:
    integer :: low, high, middle, level, goal, row, n_elem, last_low, last_high

    n_elem = ar%l
    low = 1
    high = ringlist_size

    level = 1
    pos = 1

    goal_found = .true.

    do while(level <= n_elem)
        goal = ar%sorted(level)

        last_low = low
        last_high = high
        DO WHILE(low <= high)! .AND. pos == -1)
            ! If item out of range, return
            if (goal < mainringlist(low)%sorted(level)) then
                pos = low
                goal_found = .false.
                return
            else if (goal > mainringlist(high)%sorted(level)) then
                pos = high+1
                goal_found = .false.
                return
            end if

            ! Now searching element middle.
            middle = (low + high)/2
            IF (goal == mainringlist(middle)%sorted(level)) THEN
                pos = middle
                exit
            ELSE IF (goal < mainringlist(middle)%sorted(level)) THEN
                high = middle-1
            ELSE
                low = middle+1
            END IF
        END DO

        ! Get the new range of the list.
        do row = pos, last_low, -1
            if(mainringlist(row)%sorted(level) == mainringlist(pos)%sorted(level)) then
                low = row
            else
                exit
            end if
        end do
        do row = pos, last_high
            if(mainringlist(row)%sorted(level) == mainringlist(pos)%sorted(level)) then
                high = row
            ELSE
                exit
            end if
        end do

        level = level+1
    end do

    if (ringlist_size == 0) then
        goal_found = .false.
        pos = 1
        return
    end if
end subroutine new_check_rp

! Modify the visibility array according to the primitive ring definition.
SUBROUTINE mod_pr(branch1, branch2, vis, pathlist)
    IMPLICIT NONE
    integer(inp), allocatable, intent(in) :: pathlist(:,:)
    integer(inp), intent(in) :: branch1(:), branch2(:)
    logical, intent(inout) :: vis(:,:)
    integer(inp) :: k, m, a, i, j
    logical :: isodd
    logical, allocatable :: mask1(:), mask2(:), tmp(:)

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
    print *, '| T(Path List)  | T(Find Ring)  | T(Rep. Check) | T(PR Check)   | T(Add Ring)   |'
    print 108, tneilist, tfindring, tcheckrepi, tcheckpr, taddring
    print *, '|_______________|_______________|_______________|_______________|_______________|'
108 format (' | ', *(f11.3,' s | '))

end subroutine print_ringno

pure function randomness(ringa)
    implicit none
    ! Input:
    type(ring), intent(in) :: ringa
    ! Output:
    real(dp) :: randomness

    integer :: i, id, t, r_count

    t = coord_data%ptype(ringa%element(1))

    do i = 2, ringa%l
        id = ringa%element(i)
        if (t /= coord_data%ptype(id)) r_count = r_count + 1
        t = coord_data%ptype(id)
    end do
    if (t /= coord_data%ptype(1)) r_count = r_count + 1

    randomness = r_count/ringa%l

end function

END MODULE rings_simple

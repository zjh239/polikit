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
    print *, info//' Starting RSA analysis with Guttman"s criteria;'
    print *, warn//' (Still working on this, results seem wrong)'
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

        CALL find_rings_gut(pathArray)

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

SUBROUTINE find_rings_gut(pathlist)
! natomring, ringatom, noprlist, numnopr, noprindex)
    IMPLICIT NONE
    ! IN:
    integer, dimension(:,:), intent(in) :: pathlist
    ! inOUT:
    ! PRIVATE:
!     type(ring), allocatable, dimension(:) :: ringlist
    integer :: mxrow
    logical, allocatable, dimension(:,:) :: vis, vis_copy
    logical, allocatable, dimension(:) :: bpoint
    integer, allocatable, dimension(:,:) :: vispl
!     integer(inp), dimension(:), allocatable :: smlst
    integer :: rma, rmb, id, q, t(1), lvl, n, mxlvl, l, j, id2, id0, n1, n2
    integer :: row, row_2
    real(dp) :: start, end

    logical :: add_success

    call cpu_time(tstart)

    mxrow = size(pathlist(:,1))
    mxlvl = size(pathlist(1,:))
    allocate(bpoint(mxrow), source = .false.)

    id0 = pathlist(1,1)

    allocate(vis(mxrow,mxrow), source=.true.)
    allocate(vis_copy(mxrow,mxrow), source=.true.)
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
        vis_copy = vis
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

                        if (vis_copy(row, row_2)) then
                        id2 = pathlist(row_2,lvl)

                            ! check for odd ring
                            if   (pathlist(row_2, lvl-1)==pathlist(row, lvl)&
                            .and. pathlist(row_2, lvl)==pathlist(row, lvl-1)) then
                                call cpu_time(tcheck)
                                tfindring = tfindring + (tcheck - tstart)
                                crude_ring_num = crude_ring_num +1

                                add_success = add_ring(pathlist(row,:lvl), pathlist(row_2,:lvl))

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

                                add_success = add_ring(pathlist(row,:lvl), pathlist(row_2,:lvl))

                                call mod_pr_gut(pathlist(row,:lvl), pathlist(row_2,:lvl), vis, pathlist)
                                call cpu_time(tstart)
                            end if
                        end if
                    end if
                end do
            end if
    !     call printl(vb)
        end do

!         if (all(vis .eqv. .false.)) then
!             exit
!         end if
    end do

    deallocate(bpoint)

    call cpu_time(tcheck)
    tfindring = tfindring + (tcheck - tstart)

END SUBROUTINE find_rings_gut

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

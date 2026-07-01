!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main program

program main
    use flags
    use parser
    use init
    use logger

    use neighbor_finder
    use poly_analysis
    use rdf
    use tct
    use bad
    use dynamic_data
    use rings_simple
    use data_input
    use ha
    use d2min
    use cluster

    use omp_test
    implicit none

    real :: start, finish
    call logger_init()
    call cpu_time(start)
! put test code here
!     call t1()

!
    call get_input_options() ! get the command line input strings

    if (flag_test .eqv. .true.) then
        print *, 'Performing test module ...'
        call test_run()

        call cpu_time(finish)
        print '("Wall time = ",f7.2," seconds.")', finish-start
    else if (static .eqv. .true.) then
        print *, warn//" Starting static analysis;"
        call static_analysis(0)

        call cpu_time(finish)
        print '(a,f7.2,a)', ' '//warn//" Wall time = ",(finish-start)/omp_get_max_threads()," seconds."
    else
        if (frame_interval == 0) stop error//' Frame interval missing for dynamic analysis!'
        call dynamic()

        call cpu_time(finish)
        if (finish-start>0.01) &
        print '(a,f7.2,a)', ' '//warn//" Wall time = ",(finish-start)/omp_get_max_threads()," seconds."
!         print *,
!         print *, (interval < 0), (flag_nc .eqv. .true.)
!         stop error//' Dynamic analysis report error'
    end if

contains

! this subroutine is for dynamic comparison, which means a constant interval between current
! frame and reference frame is given. It checks if the interval meet and the results should be compared and exported.
SUBROUTINE dynamic()
    IMPLICIT NONE
    integer :: fcounter

    do fcounter = skip_frame+1, fnumber
        if (fcounter == 1) print *, warn//" Starting dynamic analysis ..."

        file_name = trim(fnames(fcounter))
        print '(A, A, A, I0)', warn//' Performing dynamic analysis on ', trim(file_name), ', this is frame number ', fcounter

        call static_analysis(fcounter)
        call collect_data(fcounter)
        if (fcounter > frame_interval + skip_frame) call compare_data()
        if (fcounter > frame_interval + skip_frame) call static_post_d(fcounter)

        call mem_clean()

        print *, info//' ---------------------------End of Frame-----------------------------'
    end do

END SUBROUTINE dynamic

subroutine collect_data(cur_frame)
    implicit none
    ! IN:
    integer, intent(in) :: cur_frame

    if (flag_nc)        call collect_neighbor(cur_frame)
    if (flag_d2min)     call collect_xyz(cur_frame)

end subroutine collect_data

! D2min analysis need the position vectors of past to compare.
subroutine collect_xyz(cur_frame)
    implicit none
    ! IN:
    integer, intent(in) :: cur_frame
    !
    integer :: n

    if (frame_interval==0) then
    ! Ref is the first frame.
        n = 1

        if (.not. allocated(xyz_bf)) then
            allocate(xyz_bf(n,natom,3))
            xyz_bf = 0
        end if

        if (.not. allocated(box_bf)) then
            allocate(box_bf(n,3))
            box_bf = 0
        end if
        if (cur_frame==1) then
            xyz_bf(n,:,:) = coord_data%coord
            box_bf(n,:) = [coord_data%lx, coord_data%ly, coord_data%lz]
        end if
    else
    ! Ref is a dynamic frame.
        n = frame_interval+1

        if (.not. allocated(xyz_bf)) then
            allocate(xyz_bf(n,natom,3))
            xyz_bf = 0
        end if

        if (.not. allocated(box_bf)) then
            allocate(box_bf(n,3))
            box_bf = 0
        end if

        xyz_bf(:n-1,:,:) = xyz_bf(2:,:,:)
        xyz_bf(n,:,:) = coord_data%coord

        box_bf(:n-1,:) = box_bf(2:,:)
        box_bf(n,:) = [coord_data%lx, coord_data%ly, coord_data%lz]
    end if

end subroutine collect_xyz

!
subroutine collect_neighbor(cur_frame)
    implicit none
    !IN:
    integer, intent(in) :: cur_frame
    !
    integer :: n

    if (frame_interval == 0) then
    ! Ref is the first frame.
        n = 2

        if (.not. allocated(n_neighbor_bf)) then
            allocate(n_neighbor_bf(n,natom))
            n_neighbor_bf = 0
        end if

        if (.not. allocated(neighbors_bf)) then
            allocate(neighbors_bf(n,natom, n_cap))
            neighbors_bf = 0
            print *, 'Initializaing neighbor list data container ...'
        end if

        n_neighbor_bf(n,:) = neigh_list%n_neighbor   ! First column is n_neighbor
        neighbors_bf(n,:,:) = neigh_list%neighbors   ! Then the neighbors.

        if (cur_frame==1) then
            n_neighbor_bf(1,:) = neigh_list%n_neighbor
            neighbors_bf(1,:,:) = neigh_list%neighbors
        end if
    else
    ! Ref is a dynamic frame.
        n = frame_interval+1
!         print *, n

        if (.not. allocated(n_neighbor_bf)) then
            allocate(n_neighbor_bf(n,natom))
            n_neighbor_bf = 0
        end if

        if (.not. allocated(neighbors_bf)) then
            allocate(neighbors_bf(n,natom, n_cap))
            neighbors_bf = 0
            print *, 'Initializaing neighbor list data container ...'
        end if

!         print *, 'tets'

        n_neighbor_bf(:n-1,:) = n_neighbor_bf(2:,:)
        n_neighbor_bf(n,:) = neigh_list%n_neighbor   ! First column is n_neighbor
        neighbors_bf(:n-1,:,:) = neighbors_bf(2:,:,:)
        neighbors_bf(n,:,:) = neigh_list%neighbors   ! Then the neighbors.

    end if

end subroutine collect_neighbor

subroutine compare_data()
    implicit none
!     integer, intent(in) :: cur_frame

    if (flag_nc)  call compare_neighbor() ! neighbor_finder
    if (flag_d2min) call get_d2min()

!     if (flag_nf)  call neighbor_finder_old
!     if (flag_nfd) call compare_neighbor_d()
!     if (flag_poly) call compare_poly()
!     if (flag_bad) call compare_tct()

end subroutine compare_data

! Compare always the first and last column of the data.
subroutine compare_neighbor()
    implicit none
    integer :: i, j, n
    integer :: cn_increase, cn_decrease, cn_changed

    print *, 'Comparing coordinates ... Start'
    cn_increase = 0
    cn_decrease = 0
    cn_changed = 0

    n = size(neighbors_bf(:,1,1))
    associate(this_frame => neighbors_bf(n, :, :), &
            ref_frame => neighbors_bf(1,:,:), &
            this_n=>n_neighbor_bf(n, :), &
            ref_n=>n_neighbor_bf(1,:))

    do i = 1, natom
!         print *, this_frame(i,1), ref_frame(i,1)
        if (this_n(i) > ref_n(i)) then
            cn_increase = cn_increase + 1
        elseif (this_n(i) < ref_n(i)) then
            cn_decrease = cn_decrease + 1
        elseif (this_n(i) == ref_n(i)) then
            do j = 2, this_n(i)
                if (this_frame(i,j) /= ref_frame(i,j)) then
                    cn_changed = cn_changed + 1
                    exit
                end if
            end do
        end if
    end do
    print *, '******************************************'
    print *, '| cn_increase | cn_decrease | cn_changed |'
    end associate

    print *, 'z|', cn_increase, cn_decrease, cn_changed
    print *, '******************************************'
end subroutine compare_neighbor

subroutine test_run()
    implicit none
    ! This code is for abnormal memory cost test, modify if needed.
    integer :: fcounter

!     print *, fnumber
    do fcounter = 1, 100

        print *, 'Performing dynamic analysis on ', trim(file_name)

        call static_analysis(fcounter)

        if (frame_interval /= 0) then
            call collect_data(fcounter)
            call compare_data()
        end if

        call mem_clean()

        print *, '---------------------------End of Frame-----------------------------'
    end do

end subroutine test_run

subroutine mem_clean
    implicit none

    if (flag_nf .or. flag_nfd)    call clean_neighbor()

    if (flag_rdf .or. flag_wa) then
        call clean_rdf()
    end if

    if (flag_poly) call clean_poly()
    if (flag_cluster) call clean_cluster()

    call clean_xyz_data()

end subroutine

END PROGRAM main


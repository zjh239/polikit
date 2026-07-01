! cluster.f90
module cluster
    use precision
    use data_input, only: natom, coord_data
    use d2min, only: d2min_data, mean_d2min
    use neighbor_finder, only: neigh_list

    logical, allocatable :: selected(:), not_checked(:)
    integer, allocatable :: atomic_c_id(:), c_size(:)
    real(dp), allocatable :: c_pos(:,:)

    real(dp) :: ratio_th = 0.0

    contains

subroutine cluster_analysis()
    implicit none
    ! in:
!     integer, intent(in) :: c_frame
    ! private:
    integer :: i, t
    integer :: cluster_id, max_c_id, max_c_size

    print *, 'Performing cluster analysis ...'

    if (.not. allocated(selected)) allocate(selected(natom))
    if (.not. allocated(not_checked)) allocate(not_checked(natom))

    if (.not. allocated(atomic_c_id)) allocate(atomic_c_id(natom))

!     selected = .true.
    selected = d2min_data > mean_d2min*2

    cluster_id = 1
    not_checked = .true.

    do i = 1, natom
        if (selected(i) .and. not_checked(i)) then
            call set_cluster_id(i, cluster_id) ! set cluster id for all linked atoms.
            cluster_id = cluster_id+1
        end if
    end do

    max_c_id = cluster_id-1

    if (.not. allocated(c_size)) allocate(c_size(max_c_id))

    do i = 1, max_c_id
        c_size(i) = count(atomic_c_id == i)
!         print *, i, c_size(i)
    end do

    max_c_size = maxval(c_size)

    print *, ' ### Cluster Size Distribution'
    print *, '********************************'
    print *, 'c1|   Size       Count'

    do i = 1, max_c_size
        t = count(c_size == i)
        if (t /= 0) print *, i, t
    end do
    print *, 'c2| ****************************'

!     call cluster_pos(c_frame)

end subroutine cluster_analysis

! From one atom, recursively set its neighbor's cluster id to the same as the first.
recursive subroutine set_cluster_id(id, c_id)
    implicit none
    ! IN:
    integer, intent(in) :: id, c_id
    !
    integer :: checkid, i

    atomic_c_id(id) = c_id
    not_checked(id) = .false.
    do i = 1, neigh_list%n_neighbor(id)
        checkid = neigh_list%neighbors(id,i)
        if (selected(checkid) .and. not_checked(checkid)) then
            call set_cluster_id(checkid, c_id)
        end if
    end do

end subroutine set_cluster_id

! Loop over the cluster list and check one by one;
! for each cluster loop over its atoms and their cluster info in previous frame.
subroutine compare_cluster(ref_c_id, c_frame, interval)
    implicit none
    ! in:
    integer, intent(in) :: ref_c_id(:)
    integer, intent(in) :: c_frame, interval
    !
    integer, parameter :: size_threshold = 10
    integer :: k, i, n_clus
    character(len=5) :: frame_num

    ! calculate threshold ratio based on interval number.
    ratio_th = 0.9*(0.8/0.9)**(1/real(interval, dp))

    print *, 'Performing cluster inherit analysis with threshold ratio of ',ratio_th

    write(frame_num, '(i0)') c_frame

    open(unit = 276, file = "inherit_data_"//trim(frame_num)//".txt", status = "replace")
    write(276, *) "id  size  center_x  center_y  center_z  sub_cluster_size  parent_size  parent_id  inherit_type"

    n_clus = size(c_size)
    if (.not. allocated(c_pos)) allocate(c_pos(n_clus, 3))

    do i = 1, n_clus
        if (c_size(i) > size_threshold) then
            call check_inherit(i, c_size(i), ref_c_id(:))
        end if
    end do
    close(276)

    print *, 'Cluster inherit analysis ends ...'
end subroutine compare_cluster

! Input a cluster id, find its largest sub-cluster and the real size
! at reference frame.
subroutine check_inherit(target_id, t_size, ref_cluster_id)
    implicit none
    ! in:
    integer, intent(in) :: ref_cluster_id(:) ! complete atomic cluster id list at ref. frame
    integer, intent(in) :: target_id         ! the cluster id to be analyzed
    integer, intent(in) :: t_size
    !
    integer, allocatable :: ref_c_id_list(:)   ! cluster ids at ref. frame
    integer, allocatable :: id_list(:)       ! atom id list in target cluster
    integer :: tmp_size, max_size, tmp_c_id, max_c_id, max_ref_c
    integer :: i, k
    integer :: inherit_type
    logical :: ismerge, issplit, isinherit
    integer, parameter :: size_threshold = 10
!     real(dp), parameter :: ratio_th = 0.8

    max_size = 1
    ! get id list in target cluster
    allocate(id_list(t_size))
    allocate(ref_c_id_list(t_size))

    k = 1
    do i = 1, natom !t_size
        if (atomic_c_id(i) == target_id) then
            id_list(k) = i
            k = k+1
        end if
    end do

    ! get mass center.
    c_pos(target_id, :) = single_c_pos(id_list)

    ! get atomic cluster id at reference frame (ref_c_id_list)
    do i = 1, t_size
        ref_c_id_list(i) = ref_cluster_id(id_list(i))
    end do
    if (all(ref_c_id_list==0)) return

    call bubble_sort(t_size, ref_c_id_list)

    ! count id number in ref list, find the max one.
    tmp_c_id = 1
    do i = 1, t_size
        if (ref_c_id_list(i) /= tmp_c_id .and. ref_c_id_list(i) /= 0) then
            tmp_c_id = ref_c_id_list(i)
            tmp_size = count(ref_c_id_list == tmp_c_id)
            if (tmp_size > max_size) then
                max_size = tmp_size
                max_c_id = tmp_c_id
            end if
        end if
    end do

    tmp_c_id = 1
    ! suppose there are more than one sub-clusters with max size.
    do i = 1, t_size
        if (ref_c_id_list(i) /= tmp_c_id .and. ref_c_id_list(i) /= 0) then
            tmp_c_id = ref_c_id_list(i)
            tmp_size = count(ref_c_id_list == tmp_c_id) ! get ref. cluster size
            if (tmp_size == max_size) then
                max_ref_c = count(ref_cluster_id == max_c_id)
            ! determine if it is merge, split or inherit, by comparing
            ! max_ref_c, max_size, and t_size
                if (max_ref_c > size_threshold) then
                    ismerge = max_size > ratio_th*max_ref_c
                    issplit = max_size > ratio_th*t_size
                    isinherit = ismerge .and. issplit

                    if(ismerge .or. issplit .or. isinherit) then
                        if (isinherit) then
                            inherit_type = 3
                        else if (ismerge) then
                            inherit_type = 1
                        else if (issplit) then
                            inherit_type = 2
                        end if
                        write(276, *) target_id, t_size, c_pos(target_id, 1), c_pos(target_id, 2), &
                            c_pos(target_id, 3), max_size, max_ref_c, tmp_c_id, inherit_type

!                         print *, 'Is merge: ', ismerge, '; is split: ', issplit, '; is inherit: ', isinherit
                    end if
                end if
            end if
        end if
    end do

end subroutine check_inherit

! compute center of mass for one cluster.
function single_c_pos(list) result(tmp_pos)
    implicit none
    ! in:
    integer, INTENT(IN) :: list(:)
    ! private:
    integer :: k, i, tar_size
    real(dp) :: tmp_pos(3)

    tar_size = size(list)

    ! use id mask to get average position
    associate(xyz => coord_data%coord(list, :))
        tmp_pos(1) = sum(xyz(:,1))/tar_size
        tmp_pos(2) = sum(xyz(:,2))/tar_size
        tmp_pos(3) = sum(xyz(:,3))/tar_size
    end associate

end function single_c_pos

! compute center of mass for all clusters.
subroutine cluster_pos(c_frame)
    implicit none
    ! in:
    integer, INTENT(IN) :: c_frame
    ! private:
    integer, allocatable :: id_list(:)
    integer :: k, i, c_id, tar_size, n_clus, max_size
    real(dp) :: tmp_pos(3)

    print *, 'Performing cluster position analysis ...'

    n_clus = size(c_size)
    if (.not. allocated(c_pos)) allocate(c_pos(n_clus, 3))
    max_size = maxval(c_size)
    allocate(id_list(max_size))

    do c_id = 1, n_clus
        id_list = 0
        tar_size = c_size(c_id)

        ! get atom id list
        k = 1
        do i = 1, natom !t_size
            if (atomic_c_id(i) == c_id) then
                id_list(k) = i
                k = k+1
            end if
        end do

        ! use id mask to get average position
        associate(xyz => coord_data%coord(id_list(:tar_size), :))

            tmp_pos(1) = sum(xyz(:,1))/tar_size
            tmp_pos(2) = sum(xyz(:,2))/tar_size
            tmp_pos(3) = sum(xyz(:,3))/tar_size

            c_pos(c_id,:) = [tmp_pos(1), tmp_pos(2), tmp_pos(3)]
        end associate
    end do

!     do i = 1, 20
!         print *, c_pos(i, 1), c_pos(i, 2), c_pos(i, 3), c_size(i)
!     end do
    call draw_hist(c_frame)
    print *, 'Cluster position anlysis ends ...'

end subroutine cluster_pos

! export the cluster position, size to a file.
subroutine draw_hist(c_frame)
    implicit none
    ! in:
    integer, intent(in) :: c_frame
    ! private:
    real(dp), allocatable :: dis_to_center(:)
    real(dp) :: center(2), bin_size
    real(dp) :: hist_pos(200), avrg_size(200)
    integer :: i, k
    character(len = 4) :: frame_num

    center = [coord_data%lx/2 + coord_data%xmin, coord_data%ly/2 + coord_data%ymin]
    allocate(dis_to_center(size(c_size)))

    dis_to_center = sqrt((c_pos(:,1) - center(1))**2 + (c_pos(:,2) - center(2))**2)

    ! Note: sort it when there is a need of histogram, do not sort it when using '-lpse' and '-ci' together.
!     call dist_sort(dis_to_center, c_size)

    write(frame_num, '(i0)') c_frame

    open(unit = 275, file = "cluster_pos_"//trim(frame_num)//".txt", status = "replace")

    write(275, *) "r(angstrom)  size"
    do i = 1, size(c_size)
        if (c_size(i)>2) write(275, *) dis_to_center(i), c_size(i)
    end do

    close(275)

end subroutine draw_hist

! sort the distance list and corresponding cluster size list
subroutine dist_sort(array, acom)
    implicit none
    ! inout:
    real(dp), intent(inout) :: array(:)
    integer, intent(inout) :: acom(:)
    ! private:
    real(dp) :: temp
    integer :: last, i, j, k, tmp_int

    last = size(array)
    do i=last-1,1,-1
        do j=1,i
            if (array(j+1).lt.array(j)) then
                temp=array(j+1)
                array(j+1)=array(j)
                array(j)=temp

                tmp_int = acom(j+1)
                acom(j+1) = acom(j)
                acom(j) = tmp_int
            endif
        enddo
    enddo

end subroutine dist_sort

subroutine clean_cluster()
    implicit none

    if (allocated(atomic_c_id)) atomic_c_id = 0
    if (allocated(c_size)) deallocate(c_size)
    if (allocated(c_pos)) deallocate(c_pos)

end subroutine clean_cluster

end module cluster

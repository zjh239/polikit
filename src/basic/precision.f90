MODULE precision
    IMPLICIT NONE
    save
    public
    integer, parameter :: dp = selected_real_kind(6, 37)
    integer, parameter :: inp = selected_int_kind(9)

!     integer, parameter :: inp = selected_int_kind(18)
!     integer, parameter :: dp = selected_real_kind(15, 307)

    real(dp), parameter :: pi = 3.141592653589793 , r2d = 180.0/pi
    integer :: ierr
    character(len=80) :: emsg
!     real :: start, finish

    contains

subroutine bubble_sort(last,array)
    implicit none
    ! last is the position of last element in array.
    integer(4), intent(inout) :: array(:) !length should match the main subroutine
    integer(4), intent(in) :: last
    integer(4) :: temp, i, j, k

    do i=last-1,1,-1
        do j=1,i
            if (array(j+1).lt.array(j)) then
                temp=array(j+1)
                array(j+1)=array(j)
                array(j)=temp
            endif
        enddo
    enddo
    return
end subroutine bubble_sort

! This version maintains its own stack, to avoid needing to call
! itself recursively. By always pushing the larger "half" to the
! stack, and moving directly to calculate the smaller "half",
! it can guarantee that the stack needs no more than log_2(N)
! entries
subroutine quicksort_nr(array)
real(dp), intent(inout)::array(:)
real(dp) :: temp,pivot
integer(4) :: i,j,left,right,low,high
! If your compiler lacks storage_size(), replace
! storage_size(i) by 64
integer(4) :: stack(2, 64), stack_ptr

low=1
high=size(array)
stack_ptr=1

do
    if (high-low.lt.50) then ! use insertion sort on small arrays
        do i=low+1,high
            temp=array(i)
            do j=i-1,low,-1
            if (array(j).le.temp) exit
            array(j+1)=array(j)
            enddo
            array(j+1)=temp
        enddo
        ! now pop from stack
        if (stack_ptr.eq.1) return
        stack_ptr=stack_ptr-1
        low=stack(1,stack_ptr)
        high=stack(2,stack_ptr)
        cycle
    endif

    ! find median of three pivot
    ! and place sentinels at first and last elements
    temp=array((low+high)/2)
    array((low+high)/2)=array(low+1)
    if (temp.gt.array(high)) then
        array(low+1)=array(high)
        array(high)=temp
    else
        array(low+1)=temp
    endif
    if (array(low).gt.array(high)) then
        temp=array(low)
        array(low)=array(high)
        array(high)=temp
    endif
    if (array(low).gt.array(low+1)) then
        temp=array(low)
        array(low)=array(low+1)
        array(low+1)=temp
    endif
    pivot=array(low+1)

    left=low+2
    right=high-1
    do
        do while(array(left).lt.pivot)
            left=left+1
        enddo
        do while(array(right).gt.pivot)
            right=right-1
        enddo
        if (left.ge.right) exit
        temp=array(left)
        array(left)=array(right)
        array(right)=temp
        left=left+1
        right=right-1
    enddo
    if (left.eq.right) left=left+1
    !          call quicksort(array(1:left-1))
    !          call quicksort(array(left:))
    if (left.lt.(low+high)/2) then
        stack(1,stack_ptr)=left
        stack(2,stack_ptr)=high
        stack_ptr=stack_ptr+1
        high=left-1
    else
        stack(1,stack_ptr)=low
        stack(2,stack_ptr)=left-1
        stack_ptr=stack_ptr+1
        low=left
    endif

enddo
end subroutine quicksort_nr

! Sort two arrays, one is the driving array and the other is the idler array.
subroutine para_sort(d_array, i_array)
    implicit none
    real(dp), intent(inout) :: d_array(:)
    integer, intent(inout) :: i_array(:,:)

    real(dp) :: temp,pivot
    integer, allocatable :: i_temp(:)
    integer(4) :: i,j,left,right,low,high
    ! If your compiler lacks storage_size(), replace
    ! storage_size(i) by 64
    integer(4) :: stack(2, storage_size(i)), stack_ptr

    low=1
    high=size(d_array)
    stack_ptr=1

    allocate(i_temp(size(i_array(1,:))))

    do
        if (high-low.lt.50) then ! use insertion sort on small arrays
            do i=low+1,high
                temp=d_array(i)
                i_temp = i_array(i,:)
                do j=i-1,low,-1
                    if (d_array(j).le.temp) exit
                    d_array(j+1)=d_array(j)
                    i_array(j+1,:)=i_array(j,:)
                enddo
                d_array(j+1)=temp
                i_array(j+1,:)=i_temp
            enddo
            ! now pop from stack
            if (stack_ptr.eq.1) return
            stack_ptr=stack_ptr-1
            low=stack(1,stack_ptr)
            high=stack(2,stack_ptr)
            cycle
        endif

        ! find median of three pivot
        ! and place sentinels at first and last elements
        temp=d_array((low+high)/2)
        i_temp=i_array((low+high)/2,:)
        d_array((low+high)/2)=d_array(low+1)
        i_array((low+high)/2,:)=i_array(low+1,:)
        if (temp.gt.d_array(high)) then
            d_array(low+1)=d_array(high)
            i_array(low+1,:)=i_array(high,:)
            d_array(high)=temp
            i_array(high,:)=i_temp
        else
            d_array(low+1)=temp
            i_array(low+1,:)=i_temp
        endif
        if (d_array(low).gt.d_array(high)) then
            temp=d_array(low)
            i_temp=i_array(low,:)
            d_array(low)=d_array(high)
            i_array(low,:)=i_array(high,:)
            d_array(high)=temp
            i_array(high,:)=i_temp
        endif
        if (d_array(low).gt.d_array(low+1)) then
            temp=d_array(low)
            i_temp=i_array(low,:)
            d_array(low)=d_array(low+1)
            i_array(low,:)=i_array(low+1,:)
            d_array(low+1)=temp
            i_array(low+1,:)=i_temp
        endif
        pivot=d_array(low+1)

        left=low+2
        right=high-1
        do
            do while(d_array(left).lt.pivot)
                left=left+1
            enddo
            do while(d_array(right).gt.pivot)
                right=right-1
            enddo
            if (left.ge.right) exit
            temp=d_array(left)
            i_temp=i_array(left,:)
            d_array(left)=d_array(right)
            i_array(left,:)=i_array(right,:)
            d_array(right)=temp
            i_array(right,:)=i_temp
            left=left+1
            right=right-1
        enddo
        if (left.eq.right) left=left+1

        if (left.lt.(low+high)/2) then
            stack(1,stack_ptr)=left
            stack(2,stack_ptr)=high
            stack_ptr=stack_ptr+1
            high=left-1
        else
            stack(1,stack_ptr)=low
            stack(2,stack_ptr)=left-1
            stack_ptr=stack_ptr+1
            low=left
        endif
    enddo

end subroutine para_sort

END MODULE precision

module tct
  use precision
  use data_input, n_atom => natom
  use neighbor_finder, only: neigh_list, delta
  use parser, only: frame_interval

  implicit none

!   save
  type bond_length_data(atom_num, volumn, frame_num)
    integer, len :: atom_num, volumn, frame_num
    integer, dimension(atom_num, volumn) :: headers
    real(dp), dimension(atom_num, volumn, frame_num) :: length_data
  end type

  type bond_angle_data(atom_num, volumn, frame_num)
    integer, len :: atom_num, volumn, frame_num
    integer, dimension(atom_num, volumn, 2) :: vertices
    real(dp), dimension(atom_num, volumn, frame_num) :: angle_data
  end type

  type constraints(atom_num)
    integer, len :: atom_num
    real(dp), dimension(atom_num) :: bend_constr
    real(dp), dimension(atom_num) :: stretch_constr
  end type

  type(bond_length_data(atom_num = :, volumn = :, frame_num = :)), allocatable :: bl_data
  type(bond_angle_data(atom_num = :, volumn = :, frame_num = :)), allocatable :: ba_data
  type(constraints(atom_num = :)), allocatable :: constraint

  contains

! subroutine get_delta_data()
!   implicit none
!   integer :: fn  ! Frame number
!   integer :: length_list_size = 10, angle_list_size = 30
!
!   integer :: atom_index, neigh_index, data_index, angle_index
!   integer :: topa, topb, k
!   real(dp) :: angle, length
!   real(dp), parameter :: pi = 3.141592653 , r2d = 180.0/pi  ! 3.141592653589793
!
!   integer, allocatable :: counter(:,:,:) ! Number of occupied place, max = frame_num
!
!   fn = frame_interval
!
!   if (.not. allocated(bl_data)) allocate(bond_length_data(n_atom, length_list_size, fn) :: bl_data)
!   if (.not. allocated(ba_data)) allocate(bond_angle_data(n_atom, angle_list_size, fn) :: ba_data)
!
!   allocate(counter(n_atom, length_list_size, 2))
!
!   associate(n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors, &
!             length_header => bl_data%headers, bond_length => bl_data%length_data, &
!             angle_header => ba_data%vertices, bond_angle => ba_data%angle_data &
!             )
!
! ! Always store the lengths and angles data to the first column.
!   bond_length(:, :, 2:fn) = bond_length(:, :, :fn - 1)  ! Move lengths data one column afterward to make space.
!   bond_angle(:, : ,2:fn) = bond_angle(:, :, :fn - 1)    ! Move angles data one column afterward to make space.
!
!   do atom_index = 1, n_atom
!     do neigh_index = 1, n_neighbor(atom_index)
!
!       ! Put length data into the list
!       length = norm2(delta(atom_index, neigh_index, :))
!       do data_index = 1, length_list_size
!         if (neighbor(atom_index, neigh_index) == length_header(atom_index, data_index)) then
!           bond_length(atom_index, data_index, 1) = length
! !           counter(atom_index, data_index, 1)=counter(atom_index, data_index, 1)+1
! !           counter(atom_index, data_index, 2)=fn
!           exit  ! Exit with found header.
!         elseif (neighbor(atom_index, neigh_index) < length_header(atom_index, data_index)) then
!           length_header(atom_index, data_index+1:)=length_header(atom_index, data_index:-1)
!           length_header(atom_index, data_index)=neighbor(atom_index, neigh_index)
!           bond_length(atom_index, data_index+1:, :)=bond_length(atom_index, data_index:-1, :)
!           bond_length(atom_index, data_index, 1)=length
! !           counter(atom_index, data_index+1:, :)=counter(atom_index, data_index:-1, :)
! !           counter(atom_index, data_index, 1)=1
! !           counter(atom_index, data_index, 2)=fn
!           exit ! Exit with new header inserted.
!         elseif (length_header(atom_index, data_index) == 0) then
!           length_header(atom_index, data_index+1:)=length_header(atom_index, data_index:-1)
!           length_header(atom_index, data_index)=neighbor(atom_index, neigh_index)
!           bond_length(atom_index, data_index, 1)=length
! !           counter(atom_index, data_index, 1)=counter(atom_index, data_index, 1)+1
! !           counter(atom_index, data_index, 2)=fn
!           exit  ! Exit with new header inserted at the end.
!         endif
!         if (data_index == length_list_size) then
!           print *, "Can't find a place to insert bond length data."
!         endif
!       end do
!
!       ! Put angle data into the list
!       do data_index = neigh_index + 1, n_neighbor(atom_index)
!         topa = neighbor(atom_index,neigh_index)
!         topb = neighbor(atom_index,data_index)
!         angle = r2d*acos((delta(atom_index, neigh_index, 1) * delta(atom_index, data_index, 1) &
!               + delta(atom_index, neigh_index, 2)*delta(atom_index, data_index, 2) &
!               + delta(atom_index,neigh_index,3)*delta(atom_index,data_index,3)) &
!               /(norm2(delta(atom_index, neigh_index, :))*norm2(delta(atom_index, data_index, :)))) ! Calculate bond angle
!
!         do angle_index = 1, angle_list_size
!           if (topa==angle_header(atom_index,angle_index,1) .and. topb==angle_header(atom_index,angle_index,2)) then
!             bond_angle(atom_index,angle_index+1:,:)=bond_angle(atom_index,angle_index:-1,:)
!             bond_angle(atom_index,angle_index,1)=angle
!             exit
!           elseif (topa<angle_header(atom_index,angle_index,1)) then
!             angle_header(atom_index,angle_index+1:,:)=angle_header(atom_index,angle_index:-1,:)
!             angle_header(atom_index,angle_index,1)=topa
!             angle_header(atom_index,angle_index,2)=topb
!             bond_angle(atom_index,angle_index+1:,:)=bond_angle(atom_index,angle_index:-1,:)
!             bond_angle(atom_index,angle_index,1)=angle
!             exit
!           elseif (topa==angle_header(atom_index,angle_index,1) .and. topb<angle_header(atom_index,angle_index,2)) then
!             angle_header(atom_index,angle_index+1:,:)=angle_header(atom_index,angle_index:-1,:)
!             angle_header(atom_index,angle_index,1)=topa
!             angle_header(atom_index,angle_index,2)=topb
!             bond_angle(atom_index,angle_index+1:,:)=bond_angle(atom_index,angle_index:-1,:)
!             bond_angle(atom_index,angle_index,1)=angle
!             exit
!           elseif (topa > angle_header(atom_index,angle_index,1) .and. angle_index>1 .and. topa <= angle_header(atom_index,angle_index-1,1)) then
!             angle_header(atom_index,angle_index+1:,:)=angle_header(atom_index,angle_index:-1,:)
!             angle_header(atom_index,angle_index,1)=topa
!             angle_header(atom_index,angle_index,2)=topb
!             bond_angle(atom_index,angle_index+1:,:)=bond_angle(atom_index,angle_index:-1,:)
!             bond_angle(atom_index,angle_index,1)=angle
!             exit
!           elseif (angle_header(atom_index,angle_index,1)==0) then
!             angle_header(atom_index,angle_index+1:,:)=angle_header(atom_index,angle_index:-1,:)
!             angle_header(atom_index,angle_index,1)=topa
!             angle_header(atom_index,angle_index,2)=topb
!             bond_angle(atom_index,angle_index+1:,:)=bond_angle(atom_index,angle_index:-1,:)
!             bond_angle(atom_index,angle_index,1)=angle
!             exit
!           endif
!           if (angle_index == angle_list_size) then
!             print *, "Can't find a place to insert bond angle data."
!           endif
!         enddo
!       enddo
!     end do
!
!     ! Remove invalid data.
!     ! If a bond is break, the latest length becomes 0, remove the corresponding header and data.
!     data_index = 1
!     do while (data_index <= length_list_size)
!       if (bond_length() == 0) then
!         ! Move length
!         ! Move header
!         ! Move angle
!         data_index = data_index - 1
!       end if
!       data_index = data_index + 1
!     end do
!
!
!     ! Remove the invalid data.
!     do neigh_index = 1, length_list_size
!       if (counter(atom_index, neigh_index,2) /= fn) then
!         counter(atom_index, neigh_index,1)=0
!       endif
!       if (counter(atom_index, neigh_index,2) > 0) then
!         counter(atom_index, neigh_index,2)=counter(atom_index, neigh_index,2)-1
!       endif
!     enddo
!     neigh_index = 1
!     do ! Clear space for bond length
!       if (counter(atom_index, neigh_index,1)==0 .and. length_header(atom_index, neigh_index) > 0) then
!         bond_length(atom_index, neigh_index:, :) = eoshift(bond_length(atom_index, neigh_index:, :), shift = 1)
! !         bond_length(atom_index, neigh_index:-1,:)=bond_length(atom_index, neigh_index+1:,:)
! !         bond_length(atom_index, length_list_size,:) = 0
!         angle_index = 1
!         do !l = angle_index, 20
!           do k = 1,2
!             !clear space for bond angle
!             if (angle_header(atom_index, angle_index, k)==length_header(atom_index, neigh_index)) then
!               bond_angle(atom_index, angle_index:, :) = eoshift(bond_angle(atom_index, angle_index:, :), shift = 1)
! !               bond_angle(atom_index, angle_index:-1,:)=bond_angle(atom_index, angle_index+1:, :)
! !               bond_angle(atom_index, 30, :)=0
!               angle_header(atom_index, angle_index:, :) = eoshift(angle_header(atom_index, angle_index:, :), shift = 1)
! !               angle_header(atom_index, angle_index:-1, :)=angle_header(atom_index, angle_index+1:, :)
! !               angle_header(atom_index, 30, :)=0
!               if (angle_index < 29) then
!               angle_index = angle_index - 1
!               endif
!               exit
!             endif
!           enddo
!           angle_index = angle_index + 1
!           ! if (l==30) then
!           !   print *, i, j, l
!           ! endif
!           if (angle_index > 30) exit
!         enddo
!         length_header(atom_index, neigh_index:) = eoshift(length_header(atom_index, neigh_index:), shift = 1)
! !         length_header(atom_index, neigh_index:-1)=length_header(atom_index, neigh_index+1:)
! !         length_header(atom_index, length_list_size) = 0
!         counter(atom_index, neigh_index:-1,:)=counter(atom_index, neigh_index+1:,:)
!         counter(atom_index, length_list_size, :)=0
!         neigh_index = length_list_size - 1
!       end if
!       neigh_index = neigh_index + 1
!       if (neigh_index > 10) exit
!     enddo
!   end do
!   end associate
! !     deallocate(bond_length)
! !   deallocate(bond_angle)
! !   deallocate(length_header)
! !   deallocate(angle_header)
! !   deallocate(counter)
!
! end subroutine get_delta_data
!
! SUBROUTINE move_angle_data()
!   IMPLICIT NONE
!   dimension(:,:), intent(inout) :: array
!
!
!
! END SUBROUTINE move_angle_data

subroutine tct_calculate_old()
! Constraint is calculated as follow:
! 1. For each atom, 2 lists are created for bond lengths and bond angles correspond
!   to it. It is a n*ndev list for each atom, ndev is the number of frames will be
!   used to calculate standard deviation.
! 2. When having enough data, which means a bond is stable and unbroken for enough
!   time, standard deviation of its corralated angles are calculated.
! 3. If standard deviation is lower than given threshold, it is considered as an
!   active constraint.

  implicit none

  integer(inp) :: i,j,k,l,topa,topb,tm,tm_list(10)

  integer :: fn  ! Frame number
  integer :: length_list_size = 10, angle_list_size = 30

!   integer(inp) :: ndev=2  ! frame number to calculate a statistical result
  real(dp) :: thv=15.   ! threshold value for constraint determine
  real(dp) :: angle, length, std_dva
  real(dp), allocatable :: bond_length(:,:,:), bond_angle(:,:,:)
  integer(inp), allocatable :: length_header(:,:), angle_header(:,:,:)
  integer, allocatable :: counter(:,:,:) ! Number of occupied place, max = frame_num

  if (.not. allocated(bl_data)) allocate(bond_length_data(n_atom, length_list_size, fn) :: bl_data)
  if (.not. allocated(ba_data)) allocate(bond_angle_data(n_atom, angle_list_size, fn) :: ba_data)

  if (.not. allocated(counter)) allocate(counter(n_atom, length_list_size, 2))

  if (.not. allocated(constraint)) allocate(constraints(n_atom) :: constraint)

  associate(n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors, &
            length_header => bl_data%headers, bond_length => bl_data%length_data, &
            angle_header => ba_data%vertices, bond_angle => ba_data%angle_data, &
            bs_constr => constraint%stretch_constr, bb_constr => constraint%bend_constr)

  ! Always store the lengths and angles data from the first column.
  bond_length(:, :, 2:fn) = bond_length(:, :, :fn - 1) ! Move lengths data one column ahead to make space.
  bond_angle(:, : ,2:fn) = bond_angle(:, :, :fn - 1) ! Move angles data one column ahead to make space.

  do i =1, n_atom
    do j=1, length_list_size
      if (counter(i,j,2) /= fn) then
        counter(i,j,1)=0
      endif
      if (counter(i,j,2) > 0) then
        counter(i,j,2)=counter(i,j,2)-1
      endif
    enddo
    j = 1
    do ! Clear space for bond length
      if (counter(i,j,1)==0 .and. length_header(i,j) > 0) then
        bond_length(i,j:,:) = eoshift(bond_length(i,j:,:), shift=1)

        l=1
        do !l = 1,20
          do k = 1,2
            !clear space for bond angle
            if (angle_header(i,l,k)==length_header(i,j)) then
              bond_angle(i,l:,:) = eoshift(bond_angle(i,l:,:), shift = 1)
              angle_header(i,l:,:) = eoshift(angle_header(i,l:,:), shift = 1)
              if (l<29) then
              l=l-1
              endif
              exit
            endif
          enddo
          l=l+1
          if (l>30) exit
        enddo
        length_header(i,j:) = eoshift(length_header(i,j:), shift = 1)
        counter(i,j:,:) = eoshift(counter(i,j:,:), shift = 1)
        j=j-1
      end if
      j=j+1
      if (j > length_list_size) exit
    enddo

    ! Put data into the list
    do j = 1, n_neighbor(i)
      length = norm2(delta(i,j,:))
      do l = 1, length_list_size
        if (neighbor(i,j)==length_header(i,l)) then
          bond_length(i,l,1)=length
          counter(i,l,1)=counter(i,l,1)+1
          counter(i,l,2)=fn
          exit  ! Exit if this bond length is stored
        elseif (neighbor(i,j)<length_header(i,l)) then
          length_header(i,l:) = eoshift(length_header(i,l:), shift = -1, boundary = neighbor(i,j))
          bond_length(i,l:,:) = eoshift(bond_length(i,l:,:), shift = -1)
          bond_length(i,l,1)=length
          counter(i,l+1:,:)=counter(i,l:-1,:)
          counter(i,l,1)=1
          counter(i,l,2)=fn
          exit  ! Exit with header and length newly inserted.
        elseif (length_header(i,l) == 0) then
          length_header(i,l:) = eoshift(length_header(i,l:), shift = -1, boundary = neighbor(i,j))
          bond_length(i,l,1)=length
          counter(i,l,1)=counter(i,l,1)+1
          counter(i,l,2)=fn
          exit
        endif
        if (l == length_list_size) then
          print *, 'Can"t find a place to insert bond length data.'
        endif
      enddo

      do l = j+1,n_neighbor(i)
        topa=neighbor(i,j)
        topb=neighbor(i,l)
        angle= r2d*acos((delta(i,j,1)*delta(i,l,1)+delta(i,j,2)*delta(i,l,2)&
        &+delta(i,j,3)*delta(i,l,3))/(norm2(delta(i,j,:))*norm2(delta(i,l,:)))) !calculate bond angle

        do k = 1, angle_list_size
          if (topa==angle_header(i,k,1) .and. topb==angle_header(i,k,2)) then
            bond_angle(i,k+1:,:)=bond_angle(i,k:-1,:)
            bond_angle(i,k,1)=angle
            exit
          elseif (topa<angle_header(i,k,1)) then
            angle_header(i,k+1:,:)=angle_header(i,k:-1,:)
            angle_header(i,k,1)=topa
            angle_header(i,k,2)=topb
            bond_angle(i,k+1:,:)=bond_angle(i,k:-1,:)
            bond_angle(i,k,1)=angle
            exit
          elseif (topa==angle_header(i,k,1) .and. topb<angle_header(i,k,2)) then
            angle_header(i, k:, :) = eoshift(angle_header(i, k:, :), shift = -1)
            angle_header(i,k,1)=topa
            angle_header(i,k,2)=topb
            bond_angle(i,k:,:) = eoshift(bond_angle(i,k:,:), shift = -1)
            bond_angle(i,k,1)=angle
            exit
          elseif (topa>angle_header(i,k,1) .and. k>1 .and. topa<=angle_header(i,k-1,1)) then
            angle_header(i, k:, :) = eoshift(angle_header(i, k:, :), shift = -1)
            angle_header(i,k,1)=topa
            angle_header(i,k,2)=topb
            bond_angle(i,k:,:) = eoshift(bond_angle(i,k:,:), shift = -1)
            bond_angle(i,k,1)=angle
            exit
          elseif (angle_header(i,k,1)==0) then
            angle_header(i, k:, :) = eoshift(angle_header(i, k:, :), shift = -1)
            angle_header(i,k,1)=topa
            angle_header(i,k,2)=topb
            bond_angle(i,k:,:) = eoshift(bond_angle(i,k:,:), shift = -1)
            bond_angle(i,k,1)=angle
            exit
          endif
          if (k == angle_list_size) then
            print *, 'Can"t find a place to insert bond angle data.'
          endif
        enddo
      enddo
    enddo

    tm=0
    tm_list=0
    do j=1, length_list_size
      if (counter(i,j,1)>=fn .and. counter(i,j,2)==fn) then
        tm=tm+1
        tm_list(tm)=length_header(i,j)
        !stddva of bond length
        std_dva = stddev(bond_length(i,j,:))
        if (std_dva < 1.) then
          bs_constr(i) = bs_constr(i) + 1.
        endif
        do l=j+1, length_list_size
          if (counter(i,l,1)>=fn .and. counter(i,l,2)==fn) then
            do k=1, angle_list_size
              if (angle_header(i,k,1)==length_header(i,j) .and. &
              angle_header(i,k,2)==length_header(i,l)) then
                std_dva = stddev(bond_angle(i,k,:))
                if (std_dva < thv) then   !here is the threshold value that need to be modified
                  bb_constr(i) = bb_constr(i) + 1.
                endif
                exit
              endif
            enddo
          endif
        enddo
      endif
    enddo
    bs_constr(i) = 0.5 * bs_constr(i)
    if (bb_constr(i) >= 1.) then
      bb_constr(i) = angle2ncbb(bb_constr(i))
      ! sqrt(8.*constrain(i,2)+1.)-2.
    endif
  enddo

  ! print *, bond_length(886315,2,:)
  ! print *, delta(886315,2,:)
  end associate

end subroutine tct_calculate_old

function stddev(array)
    ! a function to calculate standard deviation from a 1D array.
    implicit none
    real(dp), intent(in) :: array(:)
    real(dp) :: bar, sm, stddev
    integer(inp) :: p, q
    sm=0.
    q=size(array)
    bar = sum(array)/q
    do p = 1,q
      sm = sm+(array(p)-bar)**2.
    enddo
    stddev=sqrt(sm/q)
    return
end function stddev

function angle2ncbb(angle_number)
    ! a step function to get ncbb from stable angle number, is a little bit ambigious.
    implicit none
    real(dp), intent(in) :: angle_number
    real(dp) :: angle2ncbb
    if (angle_number < 1.) then
      angle2ncbb = 0.
    else if (angle_number < 3.) then
      angle2ncbb = 1.
    else if (angle_number < 6.) then
      angle2ncbb = 3.
    else if (angle_number < 10.) then
      angle2ncbb = 5.
    else if (angle_number < 15.) then
      angle2ncbb = 7.
    else if (angle_number < 21.) then
      angle2ncbb = 9.
    else if (angle_number < 28.) then
      angle2ncbb = 11.
    else
      angle2ncbb = 13.
    end if
    return
end function angle2ncbb

subroutine poly_volume(n_atom,n_neigh,delta,ptype,ctype,n_ctype,hull_v,std_volume)
  implicit none
  integer(inp), intent(in) :: n_atom,n_neigh(n_atom),ptype(n_atom),ctype,n_ctype
  real(dp), intent(in) :: delta(n_atom,10,3),hull_v(n_ctype)
  real(dp), intent(out) :: std_volume(n_ctype)
  real(dp) :: edge(n_atom,20), vfactor,mean
  integer(inp) :: n,k,i,j,c
  vfactor=6.*sqrt(2.)
  c=1
  do n=1,n_atom
    if (ptype(n)==ctype .and. n_neigh(n)==4) then
      k=1
      do i=1,n_neigh(n)-1
        do j=i+1,n_neigh(n)
          edge(n,k)=norm2((delta(n,j,:)-delta(n,i,:)))
          k=k+1
        enddo
      enddo
      mean = sum(edge(n,:k-1))/(k-1)
      std_volume(c)=hull_v(c)/mean**3*vfactor
      c=c+1
    endif
  enddo
end subroutine poly_volume

end module tct

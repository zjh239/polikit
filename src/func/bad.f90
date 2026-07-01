! bad.f90
module bad  ! bond angle distribution
    use precision
    use neighbor_finder, only: neigh_list, delta
    use data_input, only: coord_data, natom, o_type, ntype
    contains

subroutine bond_angle()
  implicit none
  integer(inp) :: id, j, k, c, p
  integer :: bad_id
  integer :: bad_cap
  integer :: neigh_1, neigh_2

  real(dp) :: angle = 0.0
  real(dp), dimension(:), allocatable :: bad_angle, tmp
  integer, dimension(:,:), allocatable :: type_array, type_tmp

  print *, 'Performing bond angle distribution analysis...'

  bad_cap = 100

    allocate(bad_angle(bad_cap))
    allocate(type_array(bad_cap, 3))
  ! Given a A--B--C bond, the 3 columns store the type data of A, B, C respectively.

    bad_id = 1

    associate(n_neighbor => neigh_list%n_neighbor, ptype => coord_data%ptype)

    do id = 1, natom
      if (n_neighbor(id) < 2) cycle

      do j = 1, n_neighbor(id)-1
        neigh_1 = neigh_list%neighbors(id, j)

        do k = j+1, n_neighbor(id)
          neigh_2 = neigh_list%neighbors(id, k)

          angle= r2d*acos((delta(id,j,1)*delta(id,k,1)+delta(id,j,2)*delta(id,k,2)&
          +delta(id,j,3)*delta(id,k,3))/(norm2(delta(id,j,:))*norm2(delta(id,k,:))))

          bad_angle(bad_id) = angle
          type_array(bad_id, 1)=ptype(id)
          type_array(bad_id, 2)=ptype(neigh_1)
          type_array(bad_id, 3)=ptype(neigh_2)

          bad_id = bad_id + 1

          if (bad_id == bad_cap) then
            bad_cap = bad_cap*2

            call move_alloc(bad_angle, tmp)
            allocate(bad_angle(bad_cap))
            bad_angle(:bad_id) = tmp
            deallocate(tmp)

            call move_alloc(type_array, type_tmp)
            allocate(type_array(bad_cap,3))
            type_array(:bad_id, :) = type_tmp
            deallocate(type_tmp)
          end if

        enddo
      enddo
    enddo
  end associate

  c = bad_id-1

  call para_sort(bad_angle(:c), type_array(:c, :))

  call hist_of_bad(bad_angle(:c), type_array(:c,:))

end subroutine bond_angle

! Sort the bond angles and draw histogram.
subroutine hist_of_bad(array, type_array)
  implicit none
  ! INOUT:
  real(dp), intent(inout) :: array(:)
  integer, intent(inout) :: type_array(:,:)
  !
  integer :: cap, i, k, t1, t2, t3, n
  integer, allocatable :: raw_data(:,:,:,:), bad_data(:,:)
  logical :: mask(ntype, ntype, ntype)
  real(dp) :: binedge, bin_center
  character(len=8) :: str
  character(:), allocatable :: head

  n = ntype*ntype*(ntype+1)/2
  allocate(raw_data(180, ntype, ntype, ntype))
  allocate(bad_data(180, n))

  ! Build upper triangle mask.
  do i = 1, ntype
     do k = 1, ntype
     do t1 = 1, ntype
        mask(i,t1,k) = (i <= k)
     end do
     end do
  end do

  ! Construct header line.
  head = 'b2| theta     sum    '


  do i = 1, ntype
      do t1=1, ntype
      do k = i, ntype
          write (str, '(i0,"-",i0,"-",i0)') i,t1,k
          head  = head//str
      end do
      end do
  end do

  cap = size(array)

  binedge = 1.0
  bin_center = 0.5

  raw_data = 0

  k = 1
  i = 1
  do while(i <= cap)
    t1=type_array(i,1)
    t2=type_array(i,2)
    t3=type_array(i,3)

    if (array(i) < binedge) then
      raw_data(k,t1,t2,t3) = raw_data(k,t1,t2,t3)+1
      raw_data(k,t3,t2,t1) = raw_data(k,t3,t2,t1)+1
    else
      k = k+1
      binedge = binedge + 1.0
      i = i-1
      if (k>180) exit
    end if
    i = i+1
  end do

  forall (k=1:180)
     bad_data(k,:) = pack(raw_data(k,:,:,:), mask)
  end forall

  print *, 'b1| ============================================'
  print *, head
  do i = 1,180
    print "(f8.3, *(i8))", bin_center, sum(bad_data(i,:)), bad_data(i, :)
    bin_center = bin_center + 1.0
  end do
  print *, '================================================'

end subroutine hist_of_bad

end module bad

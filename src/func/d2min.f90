module d2min
  use precision
  use data_input, only: natom
  use dynamic_data, only: xyz_bf, box_bf
  use neighbor_finder, only: neigh_list, delta, n_cap
  use parser, only: frame_interval, pbcs
  implicit none

  real(dp), allocatable :: d2min_data(:)
  real(dp) :: tstart, tcheck, tref, ttensor, td2

  real(dp) :: mean_d2min

  contains

subroutine get_d2min()
  implicit none
  integer :: atom, n, i
  real(dp) :: tmp(3,3), v_inv(3,3), j(3,3), v(3,3), w(3,3)
  real(dp) :: d_tmp(n_cap, 3), d_test(n_cap, 3), d_ref(n_cap, 3)

  tref = 0.
  ttensor = 0.
  td2 = 0.

  if (.not. allocated(d2min_data)) allocate(d2min_data(natom))
  d2min_data = 0

  do atom = 1, natom
    n = neigh_list%n_neighbor(atom)
    associate(d_now=> delta(atom,:n,:))
      d_ref = get_ref_delta(atom, neigh_list%neighbors(atom, :n))

!   d2min = (d_now - d_ref*J)^2
!   J = (d_ref^T.d_ref)^{-1}.d_ref^T.d_now

      v=0.
      w=0.
      j=0.
      v_inv=0.
      tmp=0.

      v = matmul(transpose(d_ref(:n,:)), d_ref(:n,:))

!       if (matdet3(v) < 1e-4) then
!         d2min_data(atom) = 0
!       else

      v_inv = matinv3(v)
      w = matmul(transpose(d_ref(:n,:)), d_now)
      j = matmul(v_inv, w)

      d_tmp(:n,:) = matmul(d_ref(:n,:), j)
      d_test(:n,:) = d_now - d_tmp(:n,:)

      do i = 1, n
        d2min_data(atom) = d2min_data(atom) + norm2(d_test(i,:))**2
      end do

      if (d2min_data(atom) /= d2min_data(atom)) then
        ! If the value is 'NaN', V is probably singular. Set D2min to 0 directly.
        d2min_data(atom) = 0
      end if

    end associate
  end do

  mean_d2min = sum(d2min_data)/natom
  print *, 'd| ', mean_d2min, ' is the mean value of D2min.'

  call hist_of_d2min(d2min_data)

!   print *, 'T_ref   T_tensor   T_d2'
!   print *, tref, ttensor, td2, tstart-tcheck

end subroutine get_d2min

function get_ref_delta(cid, n_list) result(delta)
  implicit none
  ! IN:
  integer, intent(in) :: n_list(:), cid
  ! xyz_bf, box_bf from dynamic data.
  ! OUT:
  real(dp) :: delta(n_cap,3)
  ! PRIVATE:
  real(dp) :: x_tmp, y_tmp, z_tmp
  real(dp) :: delta_vec(3), min_delta(3)
  integer :: n, i, j, checkid

  delta = 0

  n = size(n_list)

  ! Reference coordinates (central atom)
  x_tmp = xyz_bf(1, cid, 1)
  y_tmp = xyz_bf(1, cid, 2)
  z_tmp = xyz_bf(1, cid, 3)

  ! Loop over neighbor atoms
  do i = 1, n
    checkid = n_list(i)
    delta_vec = [xyz_bf(1, checkid, 1) - x_tmp, &
      xyz_bf(1, checkid, 2) - y_tmp, xyz_bf(1, checkid, 3) - z_tmp]

    min_delta = delta_vec

    do j = 1, 3  ! Check each dimension (x, y, z)
      ! Only check neccessary parts.
      if (delta_vec(j) > 0) then
        delta_vec(j) = delta_vec(j) - box_bf(1,j)
      else
        delta_vec(j) = delta_vec(j) + box_bf(1,j)
      end if

      ! Check if the new displacement is shorter
      if (abs(delta_vec(j)) < abs(min_delta(j))) then
        min_delta(j) = delta_vec(j)
      end if
    end do
    delta(i, :) = min_delta
  end do

end function get_ref_delta

! inverse of 3x3 matrix, copied from fortran wiki.
pure function matinv3(A) result(B)
  !! Performs a direct calculation of the inverse of a 3×3 matrix.
  real(dp), intent(in) :: A(3,3)   !! Matrix
  real(dp)             :: B(3,3)   !! Inverse matrix
  real(dp)             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end function

! determinant of 3x3 matrix.
pure function matdet3(A) result(det)
    implicit none
    real(dp), intent(in) :: A(3,3)
    real(dp) :: det

    det = A(1,1) * (A(2,2) * A(3,3) - A(2,3) * A(3,2)) &
        - A(1,2) * (A(2,1) * A(3,3) - A(2,3) * A(3,1)) &
        + A(1,3) * (A(2,1) * A(3,2) - A(2,2) * A(3,1))
end function matdet3

! Sort the d2min value and draw histogram
subroutine hist_of_d2min(array_in)
  implicit none
  ! INOUT:
  real(dp), intent(in) :: array_in(:)
  !
  real(dp) :: array(natom)
  integer :: bin(400), i, k
  real(dp) :: binedge, bin_center, bin_size

  array = array_in

  bin_size = maxval(array)/400.

  binedge = bin_size
  bin_center = bin_size/2.

  call quicksort_nr(array)

  bin = 0

  k = 1
  i = 1
  do while(i <= natom)

    if (array(i) < binedge) then
      bin(k) = bin(k)+1

    else
      k = k+1
      binedge = binedge + bin_size
      i = i-1
      if (k == 400) then
        bin(k) = natom-i
        exit
      end if
    end if
    i = i+1
  end do

  print *, "d1| D2min distribution ===================="
  do i = 1, 400
    print "(f8.2, '    ', i0)", bin_center, bin(i)
    bin_center = bin_center + bin_size
  end do

end subroutine hist_of_d2min

end module d2min


! For cluster analysis, we need a mask on the atom list to select certain part of the atoms. Because a neighbor list has already been constructed, identifying of cluster can start with a cluster id. The id can be attributed to all the atoms that can be reached by recursively searching the neighbor list of the atoms in the cluster already.

! The next step of the cluster analysis is about calculating the mass center of the cluster. Its appearing, merging and movement during the plastic deformation can be investigated.

module neighbor_finder
  use omp_lib
  use precision
  use data_types
  use stdlib_array
  use parser, only: pbcs
  use data_input, only: coord_data, natom, ntype
  use logger

  implicit none
!   save
  ! neighbor list
  type neighbor_list(atom_number, capacity)
      integer, len :: atom_number, capacity
      integer, dimension(atom_number) :: n_neighbor = 0
      integer, dimension(atom_number, capacity) :: neighbors = 0
  end type

  ! bin type define, abolished due to stability issue.
  type bins(capacity)
      integer, len :: capacity   !capacity of the bin
      integer :: n = 0
      integer :: x_pbc = 0, y_pbc = 0, z_pbc = 0
      integer, dimension(capacity) :: ids = 0  !id of the atoms in bin
  end type

  real(dp), save, allocatable :: delta(:,:,:)

  type(neighbor_list(atom_number = :, capacity = :)), allocatable :: neigh_list

  integer, dimension(:,:), allocatable :: n_by_type

  integer, save :: n_cap

  contains

! Divide the simulation box into bins, and put atoms into their corresponding bins.
SUBROUTINE create_bins(rCut, cells_n, cells_xpbc, cells_ypbc, cells_zpbc, cells_ids, xbin_max, ybin_max, zbin_max)
  IMPLICIT NONE
  ! in:
  real(dp), intent(in) :: rCut
  ! inOUT:
  integer, intent(inout) :: xbin_max, ybin_max, zbin_max
  integer, dimension(:,:,:), allocatable, intent(inout) :: cells_n, cells_xpbc, cells_ypbc, cells_zpbc
  integer, dimension(:,:,:,:), allocatable, intent(inout) :: cells_ids
  ! PRIVATE:
  real(dp), allocatable, dimension(:,:) :: realxyz
  REAL(dp) :: x_min, y_min, z_min, bin_size_factor
  INTEGER :: xbin, ybin, zbin, atom, bincap
  integer :: i, amount, maxn

  bin_size_factor = 1.0

198 bincap = bin_size_factor*(rCut**3)
199 format (' ',a,' ', i0,a)
  print 199, info,bincap," is the bin capacity;"

  associate(xyz0 => coord_data%coord, lx => coord_data%lx, ly => coord_data%ly, lz => coord_data%lz)

    x_min = coord_data%xmin
    y_min = coord_data%ymin
    z_min = coord_data%zmin

    allocate(realxyz(natom, 3), STAT=ierr, ERRMSG=emsg)
    realxyz(:,1) = xyz0(:,1) - x_min
    realxyz(:,2) = xyz0(:,2) - y_min
    realxyz(:,3) = xyz0(:,3) - z_min

    xbin_max = CEILING(lx/rCut) - 1
    ybin_max = CEILING(ly/rCut) - 1
    zbin_max = CEILING(lz/rCut) - 1

    print 199, info, xbin_max, " bins on dimension X;"
    print 199, info, ybin_max, " bins on dimension Y;"
    print 199, info, zbin_max, " bins on dimension Z;"

  end associate

  allocate(cells_n(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))
  allocate(cells_xpbc(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))
  allocate(cells_ypbc(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))
  allocate(cells_zpbc(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))
  allocate(cells_ids(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1, bincap))

  print '(" ",a,i0,a)', info//' Deviding box into bins, memory cost of bins: ', sizeof(cells_ids)/1024, ' KB;'
  cells_n = 0
  cells_xpbc = 0
  cells_ypbc = 0
  cells_zpbc = 0
  cells_ids = 0

!!$omp parallel do default(private) shared(cells_n, cells_ids)
  do atom = 1, natom
    print *, 'This is thread: ', OMP_GET_THREAD_NUM()
    xbin = CEILING(realxyz(atom, 1)/rCut)
    ybin = CEILING(realxyz(atom, 2)/rCut)
    zbin = CEILING(realxyz(atom, 3)/rCut)

    if (xbin > xbin_max) xbin = xbin_max
    if (ybin > ybin_max) ybin = ybin_max
    if (zbin > zbin_max) zbin = zbin_max
    if (xbin == 0) xbin = 1
    if (ybin == 0) ybin = 1
    if (zbin == 0) zbin = 1
    
    if (atom == 3610) print *, xbin, ybin, zbin


    cells_n(xbin,ybin,zbin) = cells_n(xbin,ybin,zbin) + 1
    if (cells_n(xbin,ybin,zbin) > bincap) then
      bin_size_factor = bin_size_factor*2.0
      print *, warn//" Bin capacity full, size factor expanded to: ", bin_size_factor
      deallocate(cells_n)
      deallocate(cells_xpbc)
      deallocate(cells_ypbc)
      deallocate(cells_zpbc)
      deallocate(cells_ids)
      goto 198
      if (bin_size_factor > 7.0) then   ! For 1 million system, this limit the memory cost on ~GB level.
        stop error//"Bin capacity full due to unknown reason, stopping."
      end if
    end if
    cells_ids(xbin,ybin,zbin,cells_n(xbin,ybin,zbin)) = atom

  end do
!!$omp end parallel do

! Create periodic image on the right dimension.
!   if (pbc == 1) then
    do xbin = 0, xbin_max + 1
    do ybin = 0, ybin_max + 1
    do zbin = 0, zbin_max + 1

      associate(x_pbc => cells_xpbc(xbin, ybin, zbin), &
                y_pbc => cells_ypbc(xbin, ybin, zbin), &
                z_pbc => cells_zpbc(xbin, ybin, zbin))

        if (pbcs(1) == 1 .and. xbin == 0) x_pbc = 1
        if (pbcs(1) == 1 .and. xbin == xbin_max + 1) x_pbc = -1
        if (pbcs(2) == 1 .and. ybin == 0) y_pbc = 1
        if (pbcs(2) == 1 .and. ybin == ybin_max + 1) y_pbc = -1
        if (pbcs(3) == 1 .and. zbin == 0) z_pbc = 1
        if (pbcs(3) == 1 .and. zbin == zbin_max + 1) z_pbc = -1
        cells_n(xbin, ybin, zbin) = &
        cells_n(xbin + x_pbc*xbin_max, ybin + y_pbc*ybin_max, zbin + z_pbc*zbin_max)
        cells_ids(xbin, ybin, zbin, :) = &
        cells_ids(xbin + x_pbc*xbin_max, ybin + y_pbc*ybin_max, zbin + z_pbc*zbin_max, :)
      end associate
    end do
    end do
    end do
!   endif

    maxn = maxval(cells_n)

    if (maxn > 16) maxn = 16

    do i = 0, maxn
      amount = count(cells_n == i)
      print '(a,i0,a,i0,a)', ' '//info//' ', amount, ' bins with atom number ', i,';'
    end do

    print *, info//' Leaving bins construction subroutine ...'

END SUBROUTINE create_bins

SUBROUTINE find_neighbors(cutoffs)
  IMPLICIT NONE
  ! IN:
  real(dp) :: cutoffs(:)
  !
  integer :: xbin_max, ybin_max, zbin_max
  integer :: xbin, ybin, zbin, atom, atom2, i, p, q, o
  integer :: id, checkid, type1, type2
  real(dp) :: d, x_tmp, y_tmp, z_tmp

  integer, dimension(:,:,:), allocatable :: cells_n, cells_xpbc, cells_ypbc, cells_zpbc
  integer, dimension(:,:,:,:), allocatable :: cells_ids

  real(dp), allocatable :: r(:,:)

  if (.not. allocated(r)) allocate(r(ntype, ntype))

  if (size(cutoffs) == 1) then
    r = cutoffs(1)
  else
    i = 1
    do p = 1, ntype
      do q = p, ntype
        r(p,q) = cutoffs(i)
        r(q,p) = cutoffs(i)
        i = i+1
      end do
    end do
  end if

  n_cap = 2*maxval(r)**3
  print *, info//' Capacity of neighbor list is set to: ', n_cap

  if (.not. allocated(neigh_list)) allocate(neighbor_list(atom_number = natom, capacity = n_cap) :: neigh_list)
  print '(a,i0,a)', ' '//info//' Constructing neighbor list, memory cost: ', sizeof(neigh_list%neighbors)/1024, ' KB;'
  neigh_list%n_neighbor = 0
  neigh_list%neighbors = 0

  if (.not. allocated(n_by_type)) allocate(n_by_type(natom, ntype))
  n_by_type = 0

  d = maxval(r(:,:))
  r = r**2

  call create_bins(d, cells_n, cells_xpbc, cells_ypbc, cells_zpbc, cells_ids, xbin_max, ybin_max, zbin_max)

  associate(xyz => coord_data%coord, n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)

!$omp parallel do
  do xbin = 1, xbin_max
  do ybin = 1, ybin_max
  do zbin = 1, zbin_max

      if (cells_n(xbin, ybin, zbin) == 0) cycle
      do atom = 1, cells_n(xbin, ybin, zbin)
        id = cells_ids(xbin, ybin, zbin, atom)
        type1 = coord_data%ptype(id)

!         print *, id, xbin, ybin, zbin, xyz(id,1), xyz(id,2), xyz(id,3)

        do p = -1, 1
        do q = -1, 1
        do o = -1, 1

        associate(checked_n => cells_n(xbin+p, ybin+q, zbin+o),&
                  checked_ids => cells_ids(xbin+p, ybin+q, zbin+o, :),&
                  x_pbc => cells_xpbc(xbin+p, ybin+q, zbin+o), &
                  y_pbc => cells_ypbc(xbin+p, ybin+q, zbin+o), &
                  z_pbc => cells_zpbc(xbin+p, ybin+q, zbin+o))

          do atom2 = 1, checked_n
            checkid = checked_ids(atom2)
            type2 = coord_data%ptype(checkid)

            if (checkid < id) then   !to avoid repeat calculation

              x_tmp = xyz(checkid,1) - x_pbc*coord_data%lx
              y_tmp = xyz(checkid,2) - y_pbc*coord_data%ly
              z_tmp = xyz(checkid,3) - z_pbc*coord_data%lz

              d = (x_tmp - xyz(id,1))**2&     !x
                    +(y_tmp - xyz(id,2))**2&  !y
                    +(z_tmp - xyz(id,3))**2   !z

              if (d < r(type1, type2)) then
                n_neighbor(checkid) = n_neighbor(checkid) + 1
                n_neighbor(id) = n_neighbor(id) + 1
                neighbor(checkid, n_neighbor(checkid)) = id
                neighbor(id, n_neighbor(id)) = checkid

                n_by_type(id, coord_data%ptype(checkid)) = n_by_type(id, coord_data%ptype(checkid))+1
                n_by_type(checkid, coord_data%ptype(id)) = n_by_type(checkid, coord_data%ptype(id))+1

              endif
            endif
          end do

        end associate

        end do
        end do
        end do
      end do

  end do
  end do
  end do
!$omp end parallel do
  call print_cn


    do i = 1, natom
      call bubble_sort(n_neighbor(i), neighbor(i,:))
    enddo
  end associate

END SUBROUTINE find_neighbors

SUBROUTINE find_neighbors_d(cutoffs, flag_d2min)
  IMPLICIT NONE
  ! IN:
  real(dp) :: cutoffs(:)
  logical, intent(in) :: flag_d2min
        ! for d2min analysis, the neighbor list won't be sorted so that neighbor id and vector match.
  !
  integer :: xbin_max, ybin_max, zbin_max
  integer :: xbin, ybin, zbin, atom, atom2, i, p, q, o
  integer :: id, checkid, type1, type2
  real(dp) :: d, x_tmp, y_tmp, z_tmp

  integer, dimension(:,:,:), allocatable :: cells_n, cells_xpbc, cells_ypbc, cells_zpbc
  integer, dimension(:,:,:,:), allocatable :: cells_ids

  real(dp), allocatable :: r(:,:)

  print *, info//'Entering neighbor list constructing function ...'

  if (.not. allocated(r)) allocate(r(ntype, ntype))

  if (size(cutoffs) == 1) then
    r = cutoffs(1)
  else
    i = 1
    do p = 1, ntype
      do q = p, ntype
        r(p,q) = cutoffs(i)
        r(q,p) = cutoffs(i)
        i = i+1
      end do
    end do
  end if

  n_cap = 2*maxval(r)**3
  print *, info//' Capacity of neighbor list is set to: ', n_cap

  if (.not. allocated(neigh_list)) allocate(neighbor_list(atom_number = natom, capacity = n_cap) :: neigh_list)
  print '(a,i0,a)', info//' Constructing neighbor list, memory cost: ', sizeof(neigh_list%neighbors)/1024, ' KB;'
  neigh_list%n_neighbor = 0
  neigh_list%neighbors = 0

  if ((.not. allocated(delta))) then
    allocate(delta(natom, n_cap, 3), STAT=ierr, ERRMSG=emsg)
    delta = 0.
    print '(a,i0,a)', ' '//info//' Constructing delta array, memory cost: ', sizeof(delta)/1024, ' KB;'
  end if

  if (.not. allocated(n_by_type)) allocate(n_by_type(natom, ntype))
  n_by_type = 0

  d = maxval(r(:,:))
  r = r**2

  call create_bins(d, cells_n, cells_xpbc, cells_ypbc, cells_zpbc, cells_ids, xbin_max, ybin_max, zbin_max)

  associate(xyz => coord_data%coord, n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)

!!$omp parallel do
  do xbin = 1, xbin_max
  do ybin = 1, ybin_max
  do zbin = 1, zbin_max

      do atom = 1, cells_n(xbin, ybin, zbin)
        id = cells_ids(xbin, ybin, zbin, atom)
        type1 = coord_data%ptype(id)

        do p = -1, 1
        do q = -1, 1
        do o = -1, 1

        associate(checked_n => cells_n(xbin+p, ybin+q, zbin+o),&
                  checked_ids => cells_ids(xbin+p, ybin+q, zbin+o, :),&
                  x_pbc => cells_xpbc(xbin+p, ybin+q, zbin+o), &
                  y_pbc => cells_ypbc(xbin+p, ybin+q, zbin+o), &
                  z_pbc => cells_zpbc(xbin+p, ybin+q, zbin+o))

          do atom2 = 1, checked_n
            checkid = checked_ids(atom2)
            type2 = coord_data%ptype(checkid)

            if (checkid < id) then   !to avoid repeat calculation

              x_tmp = xyz(checkid,1) - x_pbc*coord_data%lx
              y_tmp = xyz(checkid,2) - y_pbc*coord_data%ly
              z_tmp = xyz(checkid,3) - z_pbc*coord_data%lz

              d = (x_tmp - xyz(id,1))**2&     !x
                    +(y_tmp - xyz(id,2))**2&  !y
                    +(z_tmp - xyz(id,3))**2   !z

              if (d < r(type1, type2)) then
                n_neighbor(checkid) = n_neighbor(checkid) + 1
                n_neighbor(id) = n_neighbor(id) + 1
                neighbor(checkid, n_neighbor(checkid)) = id
                neighbor(id, n_neighbor(id)) = checkid

                n_by_type(id, coord_data%ptype(checkid)) = n_by_type(id, coord_data%ptype(checkid))+1
                n_by_type(checkid, coord_data%ptype(id)) = n_by_type(checkid, coord_data%ptype(id))+1

                delta(id,n_neighbor(id),1) = x_tmp - xyz(id,1)
                delta(id,n_neighbor(id),2) = y_tmp - xyz(id,2)
                delta(id,n_neighbor(id),3) = z_tmp - xyz(id,3)

                delta(checkid,n_neighbor(checkid),1) = xyz(id,1) - x_tmp
                delta(checkid,n_neighbor(checkid),2) = xyz(id,2) - y_tmp
                delta(checkid,n_neighbor(checkid),3) = xyz(id,3) - z_tmp
              endif
            endif
          end do

        end associate

        end do
        end do
        end do
      end do

  end do
  end do
  end do
!!$omp end parallel do
  call print_cn

  if (.not. flag_d2min) then
    do i = 1, natom
      call bubble_sort(n_neighbor(i), neighbor(i,:))
    enddo
  end if
! Note that when dumping the delta, the neighbor list is not sorted.

  end associate

END SUBROUTINE find_neighbors_d

subroutine print_cn()
  implicit none
  integer :: i, j
  character(len=13) :: title

  print *, info//' ### Coordination Distribution'
  print *, info//'***************************'

  call print_hist(neigh_list%n_neighbor,'c| Count | ')

  do j = 1, ntype
    do i = 1, ntype
      write (title, "(A5,I1,A1,I1,A5)") ' |   ',j,'-',i,'   | '
      call print_hist(n_by_type(trueloc(coord_data%ptype==j),i) , title)
    end do
  end do
  print *, info//'***************************'

end subroutine print_cn

subroutine print_hist(alist, title)
  implicit none
  ! IN:
  integer, intent(in) :: alist(:)
  character(len=*), intent(in) :: title
  ! PRIV:
  integer :: maxn, i
  integer, allocatable :: rank(:), amount(:)

    maxn = maxval(alist)

    if (maxn > 16) maxn = 16

    allocate(rank(0:maxn), amount(0:maxn))

    rank = [(i, i=0, maxn)]
    do i = 0, maxn
      amount(i) = count(alist == i)
    end do

    if (verify('Count',title)==0) print 117, ' | CN    | ',rank
    print 117, title,amount

    117 format (a11,*(i8, ' | '))

end subroutine print_hist

subroutine clean_neighbor
  implicit none

  if (allocated(neigh_list)) then
    neigh_list%n_neighbor = 0
    neigh_list%neighbors = 0
  end if

  if (allocated(delta)) delta = 0.
  if (allocated(n_by_type)) n_by_type = 0.

end subroutine clean_neighbor

end module neighbor_finder

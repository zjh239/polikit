! dynamic_data.f90
module dynamic_data
    use precision

    implicit none

    real(dp), allocatable, dimension(:,:,:) :: xyz_bf
    real(dp), allocatable, dimension(:,:) :: box_bf

    integer, allocatable, dimension(:,:,:) :: neighbors_bf
    integer, allocatable, dimension(:,:) :: n_neighbor_bf

    integer, allocatable, dimension(:,:,:) :: poly_list_bf

    real(dp), allocatable, dimension(:,:,:,:) :: delta_bf

    integer, allocatable, dimension(:,:) :: cluster_id_bf

!     type(neighbor_list(atom_number = :, capacity = :)), allocatable, dimension(:) :: d_neighbors

!     type(polyhedron), allocatable, dimension(:) :: d_polys
!     type(length_data), allocatable, dimension(:) :: d_bondlengths
!     type(angle_data), allocatable, dimension(:) :: d_bondangles
    contains




end module dynamic_data

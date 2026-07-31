! params.f90
! this module is used to store the parameters so that either cli or file reader can access.
module params
  use precision
  implicit none
    ! basic params
    character(len=80) :: path, file_name, conf_file   ! path that contain files, or file names
    character(len=20) :: coption, doption
    character(len=40), allocatable :: fnames(:)
    logical :: static = .true.

    ! dynamic params
    integer :: fnumber, frame_interval
    integer :: skip_frame = 0

    ! analysis params
    character(len=80) :: cutoff_str, pbc_str
    real(dp) :: d2min_r(1,1)
    real(dp) :: rdf_r
    integer :: max_ring_lim
    integer :: pbcs(3)
    real(dp), dimension(:,:), allocatable :: cutoffs
    integer :: ox_type_set
end module

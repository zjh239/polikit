! params.f90
! this module is used to store the parameters so that either cli or file reader can access.
module params
  use precision
  implicit none
    ! basic params
    character(len=80) :: path, file_name, conf_file   ! path that contain files, or file names
    integer :: pbcs(3)

    character(len=20) :: coption, doption
    character(len=20), allocatable :: fnames(:)
    integer :: fnumber, frame_interval
    integer :: skip_frame = 0
    character(len=30) :: cutoff_str, pbc_str
    logical :: static
    ! dynamic params
    real(dp), allocatable :: cutoffs(:)
    ! analysis params
    real(dp) :: rdf_r, d2min_r
    integer :: max_ring_lim

end module

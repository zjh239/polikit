! conf_parser.f90
! This module handles the input and output of config file
! In CLI mode, export a config file; in file-mode, read the config file.
module config
    use precision
    use logger  ! colored text printing
    use tomlf
    use flags
    use params
    implicit none

contains

subroutine from_config(file_name)
  implicit none

  character(len=*), intent(in) :: file_name   ! path that contain files, or file names

  type(toml_table), allocatable :: table
  type(toml_error), allocatable :: error
  type(toml_table), pointer :: base_table, anal_table, exp_table, var_table
  type(toml_array), pointer :: tmp_array

  integer :: io, i

  ! Read the TOML file
  open(file=file_name, newunit=io, status="old")
  call toml_parse(table, io, error)
  close(io)

  if (allocated(error)) then
    print '(a)', "Error parsing config.toml: " // error%message
    stop 1
  end if

  ! --- Read basic values ---
  call get_value(table, "basic", base_table)
  if (.not. associated(base_table)) stop 'Basic parameters missing!'

  call get_value(base_table, 'static', static)
  if (static) print *, info//' Analaysis mode is static;'
  
  !call get_value(base_table, 'file', file_name)

  call get_value(base_table, 'pbc', tmp_array)
  do i = 1, 3
    call get_value(tmp_array, i, pbcs(i))
  end do

  ! check analysis part one by one.
  call get_value(table, 'analysis', anal_table)

  call get_value(anal_table, "nf",    flag_nf,    .false.)
  if (flag_nf) then
    call get_value(anal_table, 'nf_params', var_table)
    if (.not. associated(var_table)) stop 'Method parameters missing!'
 

  end if

  call get_value(anal_table, "nfd",   flag_nfd,   .false.)
  call get_value(anal_table, "poly",  flag_poly,  .false.)
  call get_value(anal_table, "bad",   flag_bad,   .false.)
  call get_value(anal_table, "rstat", flag_rstat, .false.)
  call get_value(anal_table, "rdf",   flag_rdf,   .false.)
  call get_value(anal_table, "wa",    flag_wa,    .false.)
  call get_value(anal_table, "ha",    flag_ha,    .false.)
  
  if (flag_nf) then
    call get_value(anal_table, 'nf_params', var_table, requested=.false.)
    if (.not. associated(var_table)) stop 'No variable provided for requisted analysis!'
    
    !call get_value(var_table, 'cutoff', cutoffs)
  end if

end subroutine from_config
end module config

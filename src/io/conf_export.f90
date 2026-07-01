! conf_export.f90
module conf_export
    use params
    use tomlf
    implicit none
    type(toml_table), allocatable :: vars
contains

subroutine write_conf()
     implicit none
! write basic part

! write analysis part

! write export part

end subroutine write_conf

end module

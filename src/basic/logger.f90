! logger.f90
! Colored log
module logger
    use iso_fortran_env, only: output_unit
    implicit none

    interface
        ! isatty is a C function often available in Fortran environments
        logical function isatty(unit) bind(C, name="isatty")
            import :: output_unit
            integer, value :: unit
        end function isatty
    end interface


    character(len=*), parameter :: white = char(27) // '[0m'
    ! red error
    character(len=*), parameter :: red = char(27) // '[31m'
    ! green info
    character(len=*), parameter :: green  = char(27) // '[32m'
    ! yellow warning
    character(len=*), parameter :: yellow = char(27) // '[33m'
    ! blue trace
    character(len=*), parameter :: blue = char(27) // '[34m'
    ! magenta other
    character(len=*), parameter :: magenta = char(27) // '[35m'
    ! cyan debug
    character(len=*), parameter :: cyan = char(27) // '[36m'

!!!!!!!!!!!! line init
    character(len=:), allocatable :: error
    ! red error
    character(len=:), allocatable :: info
    ! green info
    character(len=:), allocatable :: warn
    ! yellow warning
    character(len=:), allocatable :: trace
    ! blue trace
    character(len=:), allocatable :: other
    ! magenta other
    character(len=:), allocatable :: debug
    ! cyan debug

contains

subroutine logger_init()
    implicit none
    logical :: use_color

    ! Check if standard output (unit 6) is a terminal
    ! Note: isatty usually takes the file descriptor (0, 1, 2)
    ! For stdout, the descriptor is 1.
    use_color = isatty(1)

    if (use_color) then
        error = red//'[-error-]'//white
        info  = green//'[--info-]'//white
        warn  = yellow//'[--warn-]'//white
        trace = blue//'[-trace-]'//white
        other = magenta//'[-other-]'//white
        debug = cyan//'[-debug-]'//white
    else
        error = '[-error-]'
        info  = '[--info-]'
        warn  = '[--warn-]'
        trace = '[-trace-]'
        other = '[-other-]'
        debug = '[-debug-]'
    end if

    print *, warn//' Initializing logger;'
end subroutine
end module logger

! cli_parser.f90
! this module reads parameters from command line.
module cli_parser
    use parser
    implicit none

contains
  SUBROUTINE get_input_options()
    IMPLICIT NONE
!
    character(len=80) :: args
    integer(inp) :: n, i, k

    n = iargc()

    if(n == 0) then
        call cite_msg()
        call help_msg()
        call version_msg()
        stop
    end if

    do i = 1, n
        call get_command_argument(i, args)

        select case (args)

        case ('-c','--config')        !computing coption, now read from config
            call get_command_argument(i+1, args)
            conf_file = trim(args)
            !call from_config(conf_file)
            print *, info//' Reading config file ', trim(conf_file)
            exit
            ! stop error//' Computing option is deprecated!'
        case ("-p")         !check if pbc is applied
            call get_command_argument(i+1, pbc_str)
            call get_pbc(pbc_str)
            print '(a,i2,i2,i2)', ' '//info//' Periodic boundary conditions are: ', pbcs
        case ("-f")
            !data is a single file
            static = .true.
            print *, info//" Reading data from a single file;"
            path=''
            call get_command_argument(i+1, args)
            file_name = trim(args)
            print *, info//" File name is ", trim(file_name)
        case ("-d")
            !data is a directory
            static = .false.
            print *, info//" Input data is a directory ... True"
            call get_command_argument(i+1, args)
            path = trim(args)
            call read_file_names()
        case('-os') ! offset
            if (static .eqv. .true.) stop error//' Should not have offset in static analysis mode.'
            call get_command_argument(i+1, args)
            read (args, *, iostat=k) frame_interval
            print *, info//' Frame interval is ', frame_interval
        case('-rdf')
            flag_rdf = .true.
            call get_command_argument(i+1, args)
            read (args, *) rdf_r
            print *, info//' RDF cutoff value is:', rdf_r
        case('-d2min')
            flag_nfd = .true.
            flag_d2min = .true.
            call get_command_argument(i+1, args)
            read (args, *) d2min_r
            print *, info//' D2min cutoff value is:', d2min_r
        case('-bad')
            flag_nfd = .true.
            flag_bad = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
        case('-ring')
            flag_nf = .true.
            flag_rstat = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
            call get_command_argument(i+2, args)
            read (args, *, iostat=k) max_ring_lim
            if (k /= 0) then
                stop error//' 2nd argument after -ring should be the limit of max. path length.'
            end if
        case('-tct')
            flag_nfd = .true.
            flag_tct = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
        case('-poly')
            flag_nf = .true.
            flag_poly = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
        case('-wa')
            flag_wa = .true.
            call get_command_argument(i+1, args)
            read (args, *) rdf_r
            print *, info//' RDF cutoff value is:', rdf_r
        case('-ha')
            flag_nf = .true.
            flag_ha = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
        case('-cluster')
            flag_cluster = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
        case('-lpse')
            flag_nfd = .true.
            flag_d2min = .true.
            flag_cluster = .true.
            flag_lpse = .true.

            call get_command_argument(i+1, args)
            read (args, *) d2min_r
            print *, info//' D2min cutoff value is:', d2min_r

            call get_command_argument(i+2, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
        case('-ci')
            flag_nfd = .true.
            flag_d2min = .true.
            flag_cluster = .true.
            flag_ci = .true.

            call get_command_argument(i+1, args)
            read (args, *) d2min_r
            print *, info//' D2min cutoff value is:', d2min_r

            call get_command_argument(i+2, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
        case('-nc')
            flag_nc=.true.
            flag_nf = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
        case('-nf')
            flag_nf = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, info//' Cutoff values are:', cutoffs
            end if
        case('-skip')
            call get_command_argument(i+1, args)
            read (args, *, iostat=k) skip_frame
            print '(a, i4, a, i4)', info//' Skipping the first ', skip_frame, ' frames and starting from frame ', skip_frame+1
        case ("--help", "-h")     !help document
            call help_msg()
        case ("--version", "-v")
            CALL version_msg()
        case default
            if(verify('-', args) == 0) then
                print *, error//' Input contains unknown variable: ', args
                stop
            end if
        end select
    end do
    return
  end SUBROUTINE

end module

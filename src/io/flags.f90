module flags
    implicit none

    logical :: flag_test = .false.

    ! static analysis flags
    LOGICAL :: flag_nf = .false.    ! neighbor finder flag.
    LOGICAL :: flag_nfd = .false.
    logical :: flag_poly = .false.  ! poly_analysis flag.
    logical :: flag_bad = .false.   ! if bond angle distribution is analyzed.
    logical :: flag_rstat = .false. ! ring statistics analysis.
    logical :: flag_rdf = .false.   ! radial distribution function.
    logical :: flag_wa = .false.    ! Wendt-Abraham parameter from RDF.
    logical :: flag_ha = .false.    ! Honeycutt-Anderson parameters.

    ! dynamic analysis flags
    logical :: flag_d2min = .false.
    logical :: flag_cluster = .false.
    logical :: flag_tct = .false.   ! tct_analysis flag.
    logical :: flag_nc = .false.    ! if neighbor change is performed.
    logical :: flag_ci = .false.    ! cluster inheritance check.
    logical :: flag_lpse = .false.  ! LPSE analysis. This actually only controls dumping of the position.

end module

! This module runs the analysis based on flags and variables.
module init
    use precision
    use flags
    ! use parser
    use logger
    use cli_parser

    use neighbor_finder
    use poly_analysis
    use rdf
    use tct
    use bad
    use dynamic_data
    use rings_simple
    use data_input
    use ha
    use d2min
    use cluster
    implicit none

contains

subroutine initialize()

    call get_input_options()

end subroutine initialize

subroutine static_analysis(cur_frame)
    implicit none
    integer, intent(in) :: cur_frame

    call get_data_from_file(file_name, path)

    if (flag_nf)  call find_neighbors(cutoffs) ! if (flag_nf)  call neighbor_finder_old
    if (flag_nfd .and. .not. flag_d2min) call find_neighbors_d(cutoffs, flag_d2min)
    if (flag_nfd .and. flag_d2min) call find_neighbors_d([d2min_r], flag_d2min)
    if (flag_rdf) call calculate_rdf()
    if (flag_wa) call wa_parameter()
    if (flag_poly) call poly_neighbor()
    if (flag_bad) call bond_angle()
    if (flag_rstat) call rsa_simple(max_ring_lim)
    if (flag_ha) call calculate_ha()

    if (flag_tct) then
        stop error//" TCT analysis can not be static."
    end if
!     call calculate_qn()

end subroutine static_analysis

! Those static analysis that should be performed after some dynamic analysis.
subroutine static_post_d(cur_frame)
    implicit none
    integer, intent(in) :: cur_frame
    real(dp) :: tmp_r(1)

    if (flag_cluster) call find_neighbors(cutoffs)
    if (flag_cluster) call cluster_analysis()

    if (flag_lpse) call cluster_pos(cur_frame)

    if (flag_ci)   call collect_cluster()
    if (flag_ci .and. cur_frame > frame_interval + skip_frame + 1) &
        call compare_cluster(cluster_id_bf(1,:), cur_frame, frame_interval)

end subroutine static_post_d

subroutine collect_cluster()
    implicit none
    ! IN:
!     integer, intent(in) :: c_frame
    !
    integer :: n

    if (frame_interval==0) then
    ! Ref is the first frame.
        stop error//' Cluster tracking can not be performed by refering to the first frame.'
    else
    ! Ref is a dynamic frame.
        n = 2
        if (.not. allocated(cluster_id_bf)) then
            allocate(cluster_id_bf(n,natom))
            cluster_id_bf = 0
        end if
        cluster_id_bf(:n-1,:) = cluster_id_bf(2:,:)
        cluster_id_bf(n,:) = atomic_c_id(:)
    end if

end subroutine collect_cluster

subroutine init_msg()
    print *,    "Polikit - Atomistic Simulation Analysis Tool"
    print *,    "    V0.4"
    print *,    "Bug report: zjh239@foxmail.com"
    print *,    "    Please kindly cite:"
    print *,    " ---- "
    print *,    "    Room temperature plasticity in amorphous SiO2 and amorphous Al2O3: A computational"
    print *,    "   and topological study. Zhang, J., Frankberg, E. J., Kalikka, J. & Kuronen, A. Acta Mater."
    print *,    "   259, 119223. https://doi.org/10.1016/j.actamat.2023.119223 (2023)."
    print *,    " ---- "
end subroutine

end module

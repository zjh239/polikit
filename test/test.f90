program verlet_list
    use omp_lib
    implicit none

    integer, parameter :: n_atoms = 100000
    integer, parameter :: max_neighbors = 00
    real(8), parameter :: r_cut = 1.0d0, r_skin = 0.5d0
    real(8) :: dx, dy, dz, r2, t1, t2

    ! Atomic positions (assuming 3D)
    real(8), dimension(n_atoms,3) :: positions
    integer, dimension(n_atoms, max_neighbors) :: neighbor_list
    integer, dimension(n_atoms) :: num_neighbors

    integer :: i, j, n_count


    ! Initialize positions (for example purposes, you should load actual data)
    call omp_set_num_threads(16)
    call cpu_time(t1)
    call random_seed()
    call random_number(positions)
    positions = positions * 10.0d0  ! Scale positions to a 10x10x10 box

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(positions, neighbor_list, num_neighbors) PRIVATE(i, j, dx, dy, dz, r2, n_count)
    do i = 1, n_atoms
        n_count = 0
        do j = 1, n_atoms
            if (i /= j) then
                dx = positions(i,1) - positions(j,1)
                dy = positions(i,2) - positions(j,2)
                dz = positions(i,3) - positions(j,3)
                r2 = dx*dx + dy*dy + dz*dz

                if (r2 < (r_cut + r_skin)**2 .and. n_count < max_neighbors) then
                    n_count = n_count + 1
                    neighbor_list(i, n_count) = j
                end if
            end if
        end do
        num_neighbors(i) = n_count
    end do

    !$OMP END PARALLEL DO

    print *, "Verlet list construction completed!"
    call cpu_time(t2)
    print *, (t2-t1)/omp_get_max_threads()

end program verlet_list

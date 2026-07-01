! omp_test.f90
module omp_test
    use omp_lib
    use logger

contains

SUBROUTINE t1()
    implicit none
    call omp_set_num_threads(8)
!$OMP PARALLEL

    PRINT *, debug//"Hello from process: ", OMP_GET_THREAD_NUM()

!$OMP END PARALLEL

END SUBROUTINE t1

end module omp_test

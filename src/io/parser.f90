! parser.f90
! This module read data file with different formats and analysis parameters.
! When performing analysis, this module is called within specific functional parts.
MODULE parser
    use precision
    use flags
    use logger  ! colored text printing
    use params  ! parameter storage module.
    implicit none
    save

contains

! read file names in the directory and sort the files.
SUBROUTINE read_file_names()
    IMPLICIT NONE
    integer :: i,j,k,m, tmp
    integer, allocatable, dimension(:) :: t
    character(len=40) :: tmp_name

    call system("ls "//trim(path)//"> fname.tmp")

    fnumber = get_n_lines("fname.tmp")

    allocate(fnames(fnumber))
    allocate(t(fnumber))

    open(unit=30, file="fname.tmp", status='old', iostat=ierr, iomsg=emsg)

    do i = 1, fnumber
        read (30,*, iostat=ierr) fnames(i)
    end do
    close(unit=30)
    call system("rm fname.tmp")

    ! This part is to sort the files.
    do i = 1, fnumber
        call get_digits(fnames(i), t(i))
        do j = 1, i-1
            if (t(i) < t(j)) then
                tmp=t(j)
                t(j) = t(i)
                t(i) = tmp

                tmp_name = fnames(j)
                fnames(j) = fnames(i)
                fnames(i) = tmp_name
            endif
        enddo
    end do

!     ! Pring the sorted files.
!     do i = 1, fnumber
!         print *, fnames(i)
!     end do

    deallocate(t)
END SUBROUTINE

! get the number of lines(frame number) from a file
FUNCTION get_n_lines(filename) RESULT(nlines)
    use precision
    implicit none
    integer :: nlines
    character(len=*) :: filename

    open(unit=21, file= filename, status='old', iostat=ierr, iomsg=emsg)
    nlines = 0
    do
        read(21, *, iostat=ierr, iomsg=emsg)
        if(ierr/=0) exit    ! print emsg if need to debug here
        nlines = nlines+1
    enddo
    close(unit=21)
    print '(a,i0," files in directory.")', ' '//info//' Found ', nlines

END FUNCTION

! Get the number part of the file names.
subroutine get_digits(filename, num)
    implicit none
    ! IN:
    CHARACTER(LEN=*), INTENT(IN) :: filename
    ! OUT:
    INTEGER, INTENT(OUT) :: num
    ! Private:
    INTEGER :: i, start, end
    CHARACTER(LEN=256) :: num_str

    start = 0
    end = 0
    num_str = ''

    ! Find first digit in the filename
    DO i = 1, LEN_TRIM(filename)
        IF (filename(i:i) >= '0' .AND. filename(i:i) <= '9') THEN
            start = i
            EXIT
        END IF
    END DO

    ! Find last digit
    IF (start > 0) THEN
        DO i = start, LEN_TRIM(filename)
            IF (.NOT. (filename(i:i) >= '0' .AND. filename(i:i) <= '9')) THEN
                end = i - 1
                EXIT
            END IF
        END DO
        IF (end == 0) end = LEN_TRIM(filename)  ! If only numbers at the end
        num_str = filename(start:end)
    END IF

    ! Convert string to number
    READ (num_str, *, IOSTAT=ierr) num
    IF (ierr /= 0) print *, error//'No number found!'  ! Default if no number found

end subroutine get_digits

! Read pair-wise cutoffs.
function get_cutoff(str_in) result(r_list)
    implicit none
    ! IN:
    character(len=*), intent(in) :: str_in
    ! PRIV:
    integer :: p, k, i, n
    ! out:
    real(dp), allocatable :: r_list(:)

    if (allocated(r_list)) return

    if (index(str_in, ",") == 0) then
        allocate(r_list(1))
        read(str_in, *) r_list(1)
    else
        n = 0
        do i = 1, len_trim(str_in)
            if (str_in(i:i) == ',') then
                n=n + 1
            end if
        end do
        allocate(r_list(n+1))
        p = 1
        k = 0
        i = 1
        do while(index(str_in(p:), ",") /= 0)
            k = index(str_in(p:), ",") + k
            read(str_in(p:k-1),*) r_list(i)
            p = k+1
            i = i+1
        end do
        read(str_in(p:),*) r_list(i)
    end if
end function get_cutoff

! Read PBCs if more than 1 are given.
subroutine get_pbc(str_in)
    implicit none
    ! IN:
    character(len=*), intent(in) :: str_in
    ! PRIV:
    integer :: p, k, i

    if (index(str_in, ",") == 0) then
!         if (.not. allocated(pbcs)) allocate(pbcs(3))
        read(str_in, *) pbcs(1)
        pbcs = pbcs(1)
    else
!         if (.not. allocated(pbcs)) allocate(pbcs(3))
        p = 1
        k = 0
        i = 1
        do while(index(str_in(p:), ",") /= 0)
            k = index(str_in(p:), ",") + k
            read(str_in(p:k-1),*) pbcs(i)
            p = k+1
            i = i+1
        end do
        read(str_in(p:),*) pbcs(i)
    end if
end subroutine get_pbc

SUBROUTINE help_msg()

print *, info//" Example usage:"
print *, info//"  ./polikit -f ../test/ga2o3_test.xyz -p 1 -poly 2.3      (polyhedral analysis)"
print *, info//"  ./polikit -f ../test/ga2o3_test.xyz -p 1 -bad 2.3       (bond angle analysis)"
print *, info//"  ./polikit -f ../test/ga2o3_test.xyz -p 1 -rdf 10        (radial distribution)"
print *, info//"  ./polikit -f ../test/ga2o3_test.xyz -p 1 -wa 5          (Wendt-Abraham parameter)"
print *, info//"  ./polikit -f ../test/ga2o3_test.xyz -p 1 -ha 2.3        (Honeycutt-Anderson parameters)"
print *, info//"  ./polikit -f ../test/ga2o3_test.xyz -p 1 -ring 2.3 8    (ring statistics analysis)"
print *, info//"  ./polikit -d ../test/test_dir/ -os 3 -p 1 -nc 2.3       (dynamic neighbor change)"
print *, info//"  ./polikit -d ../test/test_dir/ -os 1 -p 1 -d2min 4.6    (LPSE inheritance analysis)"
print *, info//"  ./polikit -d ../test/test_dir/ -os 1 -p 1 -lpse 4.6 2.3 (LPSE analysis)"
print *, info//"  ./polikit -d ../test/test_dir/ -os 1 -p 1 -ci 4.6 2.3   (LPSE inheritance analysis)"
print *, info//" Variables:  "
print *, info//"   -f [string]     File name, supports .xyz, .dump, .data formats, incompatible with '-d' option."
print *, info//"                   Also be careful to match the data in each colume."
print *, info//"   -d [string]     Directory name and interval in dynamic analysis. Interval 0 will just perform  "
print *, info//"                   static analysis without comparing. Incompatible with '-f' option."
print *, info//"   -p [1 or 0]     Periodic boundary condition on all dimensions;"
print *, info//"      [1,1,1]      Or periodic boundary condition on each dimension."
print *, info//'   -os [int]       Offset for dynamic analysis, should be used with -d.'
print *, info//'   -skip [int]     Skipping the first N frames in dynamic analysis, starting from N+1 frame.'
print *, info//"   -[key] [float]  Analyzing options. 'poly' - polyhedral analysis; 'd2min' - non-affine displacement"
print *, info//"                 analysis; 'ring' - ring statistics analysis; 'bad' - bond angle distribution; 'rdf' - "
print *, info//"                 radial distribution function; 'wa' - WA parameter; 'ha' - HA parameter."
! print *, "   -o [string]        Atomic output options, won't dump atomic file if not set. 'n' - atomic coordination;"
! print *, "                        't' - tct results; 'p' - poly. neighbor; 'l' - linked state."

END SUBROUTINE

SUBROUTINE version_msg()

print *,    warn//" PoAM - Polyhedral Analysis Module"
print *,    warn//"    V0.4"
print *,    warn//" Bug report: zjh239@foxmail.com"

END SUBROUTINE

subroutine cite_msg()

    print *, warn//"    Please kindly cite:"
    print *, warn//"---"
    print *, warn//" Room temperature plasticity in amorphous SiO2 and amorphous Al2O3 : A &
    computational and topological study. Zhang, J., Frankberg, E. J., Kalikka, J. &
    Kuronen, A. Acta Mater. 259, (2023) 119223. https://doi.org/10.1016/j.actamat.2023.119223"
    print *, warn//"---"

end subroutine cite_msg

END MODULE parser

! This module read data file with different formats and analysis parameters.
! When performing analysis, this module is called within specific functional parts.
MODULE data_input
    use precision
    use data_types
    use params, only: cutoffs, pbcs
    use logger
    implicit none
    save

    ! xyz data type define
    type coordinates(atom_number)
        integer, len :: atom_number
        real(dp) :: lx = 0.0, ly = 0.0, lz = 0.0
        real(dp) :: xmin = 0.0, ymin = 0.0, zmin = 0.0
        integer, dimension(atom_number) :: ptype = 0
        real(dp), dimension(atom_number, 3) :: coord = 0.0
    end type

    integer(inp) :: natom, o_type
    integer :: ntype    ! Number of types
    real(dp), dimension(:), allocatable :: atom_frac

    type(coordinates(atom_number = :)), allocatable :: coord_data
contains

! read arguments include filename, pbc, cutoff, atomtype
SUBROUTINE get_data_from_file(file_name, path)
    IMPLICIT NONE
    ! in:
    character(len=*), intent(in)::file_name
    character(len=*), intent(in)::path
    !
    integer :: slength

!     print *, 'Reading xyz from file named ', file_name ! trim(dpath)//fname
!     slength = len(trim(fname))

    ! check the format
    if (verify('xyz', file_name) == 0) then
        print *, info//' File name: ', trim(file_name), '; XYZ format detected.'
        call read_xyz_file(file_name, path)
    else if (verify('dump', file_name) == 0) then
        print *, info//' File name: ', trim(file_name), '; LAMMPS dump format detected.'
        call read_dump_file(file_name, path)
    else if (verify('data', file_name) == 0) then
        print *, info//' File name: ', trim(file_name), '; LAMMPS data format detected.'
        call read_data_file(file_name, path)
    else
        print *, error//' Unknown format!'
        stop
    end if
    
    ! verify atom type and cutoffs match.
    if (size(cutoffs) /= 1) then
      slength = ntype*(ntype-1)/2
      if (size(cutoffs) /= slength) stop 'Using pair-wise cut-offs, but incorrect number!'
    end if
    print *, info//' Leaving get xyz subroutine ...'
END SUBROUTINE

SUBROUTINE read_xyz_file(file_name, path)
    IMPLICIT NONE
    ! in:
    character(len=*), intent(in) :: file_name
    character(len=*), intent(in) :: path
    !
    integer(inp) :: i
    character(len=10), allocatable :: typechar(:)

    open (unit=20, file=trim(path)//file_name, status='old', iostat=ierr, iomsg=emsg)

    if (ierr /=0) then
        print *, emsg
        stop
    end if

    read (20,*, iostat=ierr) natom
    print *, info, natom," atoms read from ", trim(file_name)
    read (20,*, iostat=ierr)

    if (.not. allocated(coord_data)) allocate(coordinates(natom) :: coord_data)
    allocate(typechar(natom), STAT=ierr, ERRMSG=emsg)

    associate(xyz => coord_data%coord, lx => coord_data%lx, ly => coord_data%ly, lz => coord_data%lz)
        do i=1, natom
            read (20,*, iostat=ierr) typechar(i), xyz(i,1), xyz(i,2), xyz(i,3)
        end do
        close(20)
        call type_convert(typechar)

        coord_data%xmin = MINVAL(xyz(:,1))
        coord_data%ymin = MINVAL(xyz(:,2))
        coord_data%zmin = MINVAL(xyz(:,3))

        lx = MAXVAL(xyz(:,1)) - coord_data%xmin
        ly = MAXVAL(xyz(:,2)) - coord_data%ymin
        lz = MAXVAL(xyz(:,3)) - coord_data%zmin

        PRINT *, info//" Box size on dimension X: ", coord_data%xmin, ", ", coord_data%xmin + lx
        PRINT *, info//" Box size on dimension Y: ", coord_data%ymin, ", ", coord_data%ymin + ly
        PRINT *, info//" Box size on dimension Z: ", coord_data%zmin, ", ", coord_data%zmin + lz

        deallocate(typechar)
    end associate

END SUBROUTINE

! this sub would be called if data is in .dump format(used for LAMMPS).
SUBROUTINE read_dump_file(file_name, path)
    IMPLICIT NONE
    ! in:
    character(len=*), intent(in) :: file_name
    character(len=*), intent(in) :: path
    !
    integer :: i, id
    real(dp) :: xlo, xhi, ylo, yhi, zlo, zhi
    real(dp) :: tmp_x, tmp_y, tmp_z
    character(len=10), allocatable :: typechar(:)
    character(len=10) :: tmptype

    open(unit=20, file=trim(path)//file_name, status='old', iostat=ierr, iomsg=emsg)

    read (20,*, iostat=ierr) ! ITEM: timestep
    read (20,*, iostat=ierr)
    read (20,*, iostat=ierr) ! ITEM: NUMBER OF ATOMS
    read (20,*, iostat=ierr) natom
    print '(" There are ", i0, " atoms in ", a)', natom, trim(file_name)

    if (.not. allocated(coord_data)) allocate(coordinates(natom) :: coord_data)
    allocate(typechar(natom))

    read (20,*, iostat=ierr) ! ITEM: BOX BOUNDS pp pp pp
    read (20,*, iostat=ierr) xlo, xhi
    read (20,*, iostat=ierr) ylo, yhi
    read (20,*, iostat=ierr) zlo, zhi
    read (20,*, iostat=ierr) ! ITEM: ATOMS

    coord_data%xmin = xlo
    coord_data%ymin = ylo
    coord_data%zmin = zlo

    coord_data%lx = xhi - xlo
    coord_data%ly = yhi - ylo
    coord_data%lz = zhi - zlo

    PRINT '(" Box size on dimension X: ", f10.4, ", ", f10.4, ";")', coord_data%xmin, coord_data%xmin + coord_data%lx
    PRINT '(" Box size on dimension Y: ", f10.4, ", ", f10.4, ";")', coord_data%ymin, coord_data%ymin + coord_data%ly
    PRINT '(" Box size on dimension Z: ", f10.4, ", ", f10.4, ";")', coord_data%zmin, coord_data%zmin + coord_data%lz

    associate(xyz => coord_data%coord)

        do i = 1, natom
            read (20,*, iostat=ierr) id, typechar(i), xyz(i,1), xyz(i,2), xyz(i,3)
        end do
    end associate

    close(20)
    call type_convert(typechar)
    deallocate(typechar)
END SUBROUTINE

! this sub would be called if data is in .data format(used for LAMMPS).
SUBROUTINE read_data_file(file_name, path)
    IMPLICIT NONE
    ! in:
    character(len=*), intent(in) :: file_name
    character(len=*), intent(in) :: path
    !
    integer :: i, id, ntype
    real(dp) :: xlo, xhi, ylo, yhi, zlo, zhi
    real(dp) :: tmp_x, tmp_y, tmp_z
    character(len=10), allocatable :: typechar(:)
    character(len=10) :: tmptype

    open(unit=20, file=trim(path)//file_name, status='old', iostat=ierr, iomsg=emsg)

    read (20,*, iostat=ierr) ! Comments
    read (20,*, iostat=ierr)
    read (20,*, iostat=ierr) natom  ! atoms

    print '(" There are ", i0, " atoms in ", a)', natom, trim(file_name)
    if (.not. allocated(coord_data)) allocate(coordinates(natom) :: coord_data)
    allocate(typechar(natom))

    read (20,*, iostat=ierr) ntype  ! atom types
    read (20,*, iostat=ierr)

    read (20,*, iostat=ierr) xlo, xhi   ! xlo xhi
    read (20,*, iostat=ierr) ylo, yhi   ! ylo yhi
    read (20,*, iostat=ierr) zlo, zhi   ! zlo zhi

!     do i = 1, 8
!         read (20,*, iostat=ierr)
!         ! skip the masses.
!     end do

    coord_data%xmin = xlo
    coord_data%ymin = ylo
    coord_data%zmin = zlo

    coord_data%lx = xhi - xlo
    coord_data%ly = yhi - ylo
    coord_data%lz = zhi - zlo

    PRINT '(" Box size on dimension X: ", f10.4, ", ", f10.4, ";")', coord_data%xmin, coord_data%xmin + coord_data%lx
    PRINT '(" Box size on dimension Y: ", f10.4, ", ", f10.4, ";")', coord_data%ymin, coord_data%ymin + coord_data%ly
    PRINT '(" Box size on dimension Z: ", f10.4, ", ", f10.4, ";")', coord_data%zmin, coord_data%zmin + coord_data%lz

    read (20,*, iostat=ierr)
    read (20,*, iostat=ierr)    ! Atoms  # atomic
    read (20,*, iostat=ierr)

    associate(xyz => coord_data%coord)
        do i = 1, natom
            read (20,*, iostat=ierr) id, typechar(i), xyz(i,1), xyz(i,2), xyz(i,3)
!             typechar(id) = tmptype
!             xyz(id,1) = tmp_x
!             xyz(id,2) = tmp_y
!             xyz(id,3) = tmp_z
        end do
    end associate

    close(20)
    call type_convert(typechar)
    deallocate(typechar)
END SUBROUTINE

! if type names are string in data file, we need to convert it to integer. Oxygen atom is
! automatically distinguished by assuming oxygen takes the largest part.
SUBROUTINE type_convert(charin)
    IMPLICIT NONE

    character(len=10), intent(in) :: charin(natom)
    character(len=10) :: typename(10)
    integer :: i, j

    integer :: tmp_number

!     allocate(ptype(natom), STAT=ierr, ERRMSG=emsg)

    associate(ptype => coord_data%ptype)
        ptype = 0
        ntype = 1
        typename(1) = charin(1)

        do i = 1, natom     !compare char type with each existing type, add it to typename if not exist
            do j = 1, ntype
                if(charin(i) == typename(j)) then
                    ptype(i) = j
                    exit
                else if(j==ntype) then
                    ntype = ntype + 1
                    typename(ntype) = charin(i)
                    ptype(i) = ntype
                end if
            end do

        end do
        print *, " ### Atom Type Number"
        print *, "****************************"
        print *, "Typename    Typeid    Number"

        o_type = 1
        if (.not. allocated(atom_frac)) allocate(atom_frac(ntype))

        do i = 1, ntype
            tmp_number = count(ptype == i)
            atom_frac(i) = tmp_number/real(natom)
            if (tmp_number > count(ptype==o_type)) then
                o_type = i
            end if
            print "('   ', a, i0, '         ', i0)", typename(i), i, tmp_number
        end do
        print *, "****************************"
        print *, 'Oxygen type is: ', o_type !, atom_frac
    end associate

END SUBROUTINE

subroutine clean_xyz_data
    implicit none

    if (allocated(coord_data)) then
        coord_data%lx = 0.0
        coord_data%ly = 0.0
        coord_data%lz = 0.0

        coord_data%xmin = 0.0
        coord_data%ymin = 0.0
        coord_data%zmin = 0.0

        coord_data%ptype = 0
        coord_data%coord = 0.0
    end if

end subroutine

END MODULE data_input

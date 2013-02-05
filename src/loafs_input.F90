module loafs_input

  implicit none
  private
  public :: configure_loafs 

contains

!===============================================================================
! CONFIGURE_LOAFS
!===============================================================================

  subroutine configure_loafs()
    
    use error,   only: fatal_error
    use global
    use output,  only: write_message
    use string,  only: lower_case
    use xml_data_loafs_t

    integer :: ng
    logical :: file_exists ! does loafs.xml exist?
    character(MAX_LINE_LEN) :: filename

    ! read loafs infput file
    filename = trim(path_input) // "loafs.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
    
      message = "No LOAFS XML file, '" // trim(filename) // "' does not exist!"
      call fatal_error()
      
    else

      ! tell user
      message = "Reading LOAFS XML file..."
      call write_message(5)

    end if

    ! parse loafs.xml file
    call read_xml_file_loafs_t(filename)

    ! Check number of particle sites
    if (sites_ == 0) then
      message = "Need to specify number of total starting particles."
      call fatal_error()
    end if
      
    n_particles = sites_
    
    
    if (associated(mesh_ % energy)) then
      ng = size(mesh_ % energy)
      if(.not.allocated(loafs%egrid)) allocate(loafs%egrid(ng))
      loafs%egrid = mesh_ % energy 
      loafs%n_egroups = ng - 1
    end if

  end subroutine configure_loafs

end module loafs_input

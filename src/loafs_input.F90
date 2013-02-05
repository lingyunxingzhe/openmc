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



  end subroutine configure_loafs

end module loafs_input

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

    ! read energy grid
    if (associated(mesh_ % energy)) then
      ng = size(mesh_ % energy)
      if(.not.allocated(loafs % egrid)) allocate(loafs % egrid(ng))
      loafs % egrid = mesh_ % energy 
      loafs % n_egroups = ng - 1
    else
      message = "Need to specify energy group structure."
      call fatal_error()
    end if
    
    ! read number of sites
    if (associated(mesh_ % sites)) then
      ng = size(mesh_ % sites)
      if (ng /= loafs % n_egroups) then
        message = "Size of storage sites must match the number of LOAFS zones."
        call fatal_error()
      end if
      if(.not.allocated(loafs % max_sites)) allocate(loafs % max_sites(ng))
      loafs % max_sites = mesh_ % sites
    else
      message = "Need to specify number of sites to store in each zone."
      call fatal_error()
    end if
    
    if (fission_bank_size_ >= 0) then
      loafs % fission_bank_size = fission_bank_size_
    else
      message = "Need to specify size of LOAFS fission bank."
      call fatal_error()
    end if

  end subroutine configure_loafs

end module loafs_input

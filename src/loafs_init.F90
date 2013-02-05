module loafs_init

  use error,            only: fatal_error
  use global
  use string,           only: to_str

  implicit none
  private
  public :: allocate_loafs_banks

contains

!===============================================================================
! ALLOCATE_LOAFS_BANKS
!===============================================================================

  subroutine allocate_loafs_banks()
  
    integer    :: i
    integer    :: alloc_err  ! allocation error code

    ! Allocate scatter banks
    allocate(loafs % scatter_banks(loafs % n_egroups), STAT=alloc_err)
    allocate(loafs % scatter_bank_idx(loafs % n_egroups), STAT=alloc_err)
    do i=1,loafs % n_egroups
      allocate(loafs % scatter_banks(i) % sites(n_particles), STAT=alloc_err)
      loafs % scatter_bank_idx(i) = 1
      
      ! Check for allocation errors 
      if (alloc_err /= 0) then
        message = "Failed to allocate LOAFS scatter bank for group " & 
                  // to_str(i) // "."
        call fatal_error()
      end if
    end do
    
  end subroutine allocate_loafs_banks

end module loafs_init

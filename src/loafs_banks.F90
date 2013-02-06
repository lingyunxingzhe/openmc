module loafs_banks

  use error,            only: fatal_error
  use global
  use string,           only: to_str

  implicit none
  private
  public :: allocate_loafs_banks, loafs_add_to_bank

contains

!===============================================================================
! ALLOCATE_LOAFS_BANKS
!===============================================================================

  subroutine allocate_loafs_banks()
  
    integer    :: i
    integer    :: alloc_err  ! allocation error code

    ! Allocate banks
    allocate(loafs % site_banks(loafs % n_egroups), STAT=alloc_err)
    allocate(loafs % site_bank_idx(loafs % n_egroups), STAT=alloc_err)
    do i=1,loafs % n_egroups
      allocate(loafs % site_banks(i) % sites(loafs % max_sites(i)), & 
                                                                 STAT=alloc_err)
      loafs % site_bank_idx(i) = 1
      
      ! Check for allocation errors 
      if (alloc_err /= 0) then
        message = "Failed to allocate LOAFS site bank for group " & 
                  // to_str(i) // "."
        call fatal_error()
      end if

    end do
    
  end subroutine allocate_loafs_banks

!===============================================================================
! LOAFS_ADD_TO_BANK
!===============================================================================

  subroutine loafs_add_to_bank(bin)

    integer :: bin
    
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % is_fission_site = .false.
    
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % wgt = p % wgt
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % xyz = p % coord0 % xyz
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % uvw = p % coord0 % uvw
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % E = p % E
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % mu = p % mu
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % last_xyz = p % last_xyz
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % last_wgt = p % last_wgt
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % last_E = p % last_E
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % absorb_wgt = p % absorb_wgt
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % event = p % event
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % event_nuclide = p % event_nuclide
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % event_MT = p % event_MT
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % n_bank = p % n_bank
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % wgt_bank = p % wgt_bank
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % surface = p % surface
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % cell_born = p % cell_born
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % material = p % material
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % last_material = p % last_material
    loafs % site_banks(bin) % sites(loafs % site_bank_idx(bin)) % n_collision = p % n_collision

    
  end subroutine loafs_add_to_bank

end module loafs_banks

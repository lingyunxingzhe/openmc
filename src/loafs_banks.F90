module loafs_banks

  use error,            only: fatal_error
  use global
  use particle_header,  only: deallocate_coord
  use string,           only: to_str

  implicit none

contains

!===============================================================================
! ALLOCATE_LOAFS_BANKS
!===============================================================================

  subroutine allocate_loafs_banks()
  
    integer    :: i
    integer    :: alloc_err  ! allocation error code

    allocate(loafs % extra_weights(loafs % n_egroups), STAT=alloc_err)

    ! Allocate banks
    allocate(loafs % site_banks(loafs % n_egroups), STAT=alloc_err)
    allocate(loafs % site_bank_idx(loafs % n_egroups), STAT=alloc_err)
    do i=1,loafs % n_egroups
      allocate(loafs % site_banks(i) % sites(loafs % max_sites(i)), & 
                                                                 STAT=alloc_err)
      loafs % site_bank_idx(i) = 0
      
      ! Check for allocation errors 
      if (alloc_err /= 0) then
        message = "Failed to allocate LOAFS site bank for group " & 
                  // to_str(i) // "."
        call fatal_error()
      end if

    end do

    allocate(fission_bank(loafs % fission_bank_size), STAT=alloc_err)
    
    if (alloc_err /= 0) then
      message = "Failed to allocate fission bank."
      call fatal_error()
    end if
    
  end subroutine allocate_loafs_banks

!===============================================================================
! LOAFS_PARTICLE_TO_BANK
!===============================================================================

  subroutine loafs_particle_to_bank(bin,i)

    integer, intent(in) :: bin
    integer, intent(in) :: i
    
    loafs % site_banks(bin) % sites(i) % id = p % id
    
    loafs % site_banks(bin) % sites(i) % wgt = p % wgt
    loafs % site_banks(bin) % sites(i) % xyz = p % coord0 % xyz
    loafs % site_banks(bin) % sites(i) % uvw = p % coord0 % uvw
    loafs % site_banks(bin) % sites(i) % E = p % E
    loafs % site_banks(bin) % sites(i) % mu = p % mu
    loafs % site_banks(bin) % sites(i) % last_xyz = p % last_xyz
    loafs % site_banks(bin) % sites(i) % last_wgt = p % last_wgt
    loafs % site_banks(bin) % sites(i) % last_E = p % last_E
    loafs % site_banks(bin) % sites(i) % absorb_wgt = p % absorb_wgt
    loafs % site_banks(bin) % sites(i) % event = p % event
    loafs % site_banks(bin) % sites(i) % event_nuclide = p % event_nuclide
    loafs % site_banks(bin) % sites(i) % event_MT = p % event_MT
    loafs % site_banks(bin) % sites(i) % n_bank = p % n_bank
    loafs % site_banks(bin) % sites(i) % wgt_bank = p % wgt_bank
    loafs % site_banks(bin) % sites(i) % surface = p % surface
    loafs % site_banks(bin) % sites(i) % cell_born = p % cell_born
    loafs % site_banks(bin) % sites(i) % material = p % material
    loafs % site_banks(bin) % sites(i) % last_material = p % last_material
    loafs % site_banks(bin) % sites(i) % n_collision = p % n_collision

    
  end subroutine loafs_particle_to_bank


!===============================================================================
! LOAFS_BANK_TO_PARTICLE
!===============================================================================

  subroutine loafs_bank_to_particle(bin, i)

    integer, intent(in) :: bin
    integer, intent(in) :: i
    
    call deallocate_coord(p % coord0 % next)
    p % coord => p % coord0
    
    p % id = loafs % site_banks(bin) % sites(i) % id
    p % wgt = loafs % site_banks(bin) % sites(i) % wgt
    p % coord0 % xyz = loafs % site_banks(bin) % sites(i) % xyz
    p % coord0 % uvw = loafs % site_banks(bin) % sites(i) % uvw
    p % E = loafs % site_banks(bin) % sites(i) % E
    p % mu = loafs % site_banks(bin) % sites(i) % mu
    p % last_xyz = loafs % site_banks(bin) % sites(i) % last_xyz
    p % last_wgt = loafs % site_banks(bin) % sites(i) % last_wgt
    p % last_E = loafs % site_banks(bin) % sites(i) % last_E
    p % absorb_wgt = loafs % site_banks(bin) % sites(i) % absorb_wgt
    p % event = loafs % site_banks(bin) % sites(i) % event
    p % event_nuclide = loafs % site_banks(bin) % sites(i) % event_nuclide
    p % event_MT = loafs % site_banks(bin) % sites(i) % event_MT
    p % n_bank = loafs % site_banks(bin) % sites(i) % n_bank
    p % wgt_bank = loafs % site_banks(bin) % sites(i) % wgt_bank
    p % surface = loafs % site_banks(bin) % sites(i) % surface
    p % cell_born = loafs % site_banks(bin) % sites(i) % cell_born
    p % material = loafs % site_banks(bin) % sites(i) % material
    p % last_material = loafs % site_banks(bin) % sites(i) % last_material
    p % n_collision = loafs % site_banks(bin) % sites(i) % n_collision

    p % alive = .true.
    
  end subroutine loafs_bank_to_particle

!===============================================================================
! DISTRIBUTE_EXTRA_WEIGHT
!===============================================================================

  subroutine distribute_extra_weight()

    integer :: bin, i
    real(8) :: sliver
    
    do bin = 1, loafs % n_egroups
    
      sliver = loafs % extra_weights(bin) / loafs % max_sites(bin)
      
      do i = 1, loafs % max_sites(bin)
        
        loafs % site_banks(bin) % sites(i) % wgt = &
                               loafs % site_banks(bin) % sites(i) % wgt + sliver
        
      end do
    
    end do

  end subroutine distribute_extra_weight

end module loafs_banks

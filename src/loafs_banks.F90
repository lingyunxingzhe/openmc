module loafs_banks

  use constants
  use error,            only: fatal_error
  use global
  use particle_header,  only: deallocate_coord
  use random_lcg,       only: prn
  use search,           only: binary_search
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
    allocate(loafs % source_banks(loafs % n_egroups), STAT=alloc_err)
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
      
      allocate(loafs % source_banks(i) % sites(loafs % max_sites(i)), & 
                                                                 STAT=alloc_err)
                                                                 
      ! Check for allocation errors 
      if (alloc_err /= 0) then
        message = "Failed to allocate LOAFS source bank for group " & 
                  // to_str(i) // "."
        call fatal_error()
      end if
      
    end do

    ! TODO: currently the fission bank isn't used

    allocate(fission_bank(loafs % fission_bank_size), STAT=alloc_err)
    
    if (alloc_err /= 0) then
      message = "Failed to allocate fission bank."
      call fatal_error()
    end if
    
  end subroutine allocate_loafs_banks


!===============================================================================
! LOAFS_TRANSPORT_CHECK
!===============================================================================

  subroutine loafs_transport_check()
  
    loafs_bin = binary_search(loafs % egrid, loafs % n_egroups + 1, p % E)
    
    if (loafs_site_gen) then
    
      if (loafs_bin /= loafs_last_bin) then
        if (loafs % site_bank_idx(loafs_bin) < loafs % max_sites(loafs_bin)) then
          loafs % site_bank_idx(loafs_bin) = loafs % site_bank_idx(loafs_bin) + 1
          call loafs_particle_to_bank(loafs_bin,loafs % site_bank_idx(loafs_bin))
        else
          loafs % extra_weights(loafs_bin) = loafs % extra_weights(loafs_bin) + p % wgt
        end if
      end if
      
      loafs_last_bin = loafs_bin
      
    else
    
      if (loafs_bin /= loafs_active_bin) p % alive = .false.
      
    end if
    
  end subroutine loafs_transport_check


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
    loafs % site_banks(bin) % sites(i) % n_collision = p % n_collision

  end subroutine loafs_particle_to_bank


!===============================================================================
! COPY_SITES_TO_SOURCE
!===============================================================================

  subroutine copy_sites_to_source()
  
    integer :: bin
    integer :: site_idx
    integer :: source_idx
    
    real(8) :: real_weight
    real(8) :: extra_weight
    real(8) :: ratio
    
    
    do bin = 1, loafs % n_egroups
      
      source_idx = 0  
      site_idx = 0
      extra_weight = 0.0_8
      real_weight = 0.0_8
      
      do while (source_idx < loafs % max_sites(bin))
      
        source_idx = source_idx + 1
      
        if (site_idx <= loafs % site_bank_idx(bin)) then
          ! copy from banked sites to the source
          site_idx = site_idx + 1
          real_weight = real_weight + &
                                 loafs % site_banks(bin) % sites(site_idx) % wgt
        else
          ! sample a site from the banked sites and note the extra weight
          site_idx = 0
          do while (site_idx == 0) 
            site_idx = prn()*loafs % max_sites(bin)
          end do
          extra_weight = extra_weight + &
                                 loafs % site_banks(bin) % sites(site_idx) % wgt
        end if
      
        loafs % source_banks(bin) % sites(source_idx) = &
                                       loafs % site_banks(bin) % sites(site_idx) 
      
      end do 
      
      ! reduce the weights if we had extra weights
      if (extra_weight > 0.0_8) then
      
        ratio = real_weight/(real_weight+extra_weight)
        do source_idx = 1,loafs % max_sites(bin)
          loafs % source_banks(bin) % sites(source_idx) % wgt = &
                     loafs % source_banks(bin) % sites(source_idx) % wgt * ratio
        end do
      
      end if
    
    end do
  

    
  end subroutine copy_sites_to_source

!===============================================================================
! LOAFS_SOURCE_TO_PARTICLE
!===============================================================================

  subroutine loafs_source_to_particle(bin, i)

    integer, intent(in) :: bin
    integer, intent(in) :: i
    
    call deallocate_coord(p % coord0 % next)
    p % coord => p % coord0
    p % coord % cell = NONE
    
    p % id = loafs % source_banks(bin) % sites(i) % id
    p % wgt = loafs % source_banks(bin) % sites(i) % wgt
    p % coord % xyz = loafs % source_banks(bin) % sites(i) % xyz
    p % coord % uvw = loafs % source_banks(bin) % sites(i) % uvw
    p % E = loafs % source_banks(bin) % sites(i) % E
    p % mu = loafs % source_banks(bin) % sites(i) % mu
    p % last_xyz = loafs % source_banks(bin) % sites(i) % last_xyz
    p % last_wgt = loafs % source_banks(bin) % sites(i) % last_wgt
    p % last_E = loafs % source_banks(bin) % sites(i) % last_E
    p % absorb_wgt = loafs % source_banks(bin) % sites(i) % absorb_wgt
    p % event = loafs % source_banks(bin) % sites(i) % event
    p % event_nuclide = loafs % source_banks(bin) % sites(i) % event_nuclide
    p % event_MT = loafs % source_banks(bin) % sites(i) % event_MT
    p % n_bank = loafs % source_banks(bin) % sites(i) % n_bank
    p % wgt_bank = loafs % source_banks(bin) % sites(i) % wgt_bank
    p % surface = loafs % source_banks(bin) % sites(i) % surface
    p % cell_born = loafs % source_banks(bin) % sites(i) % cell_born
    p % material = NONE
    p % last_material = NONE
    p % n_collision = loafs % source_banks(bin) % sites(i) % n_collision

    p % alive = .true.
    
  end subroutine loafs_source_to_particle

!===============================================================================
! DISTRIBUTE_EXTRA_WEIGHT
!===============================================================================

  subroutine distribute_extra_weight()

    integer :: bin, i
    real(8) :: sliver
    
    ! TODO: add a bin-wise total weight counter to the loafs object to properly
    ! re-weight with ratios instead of the additive stuff (incorrect when not
    ! all particles have the same weight)
    
    do bin = 1, loafs % n_egroups
    
      sliver = loafs % extra_weights(bin) / loafs % max_sites(bin)
      
      do i = 1, loafs % max_sites(bin)
        
        loafs % site_banks(bin) % sites(i) % wgt = &
                               loafs % site_banks(bin) % sites(i) % wgt + sliver
        
      end do
    
    end do

  end subroutine distribute_extra_weight

end module loafs_banks

module loafs_main

  use bank_header
  use global
  use output,       only: header
  use physics,      only: transport
  use random_lcg,   only: set_particle_seed
  use sort,         only: heap_sort
  use source,       only: initialize_particle, sample_external_source, &
                          copy_source_attributes

  type(Bank), pointer :: source_site => null()

contains

!===============================================================================
! RUN_LOAFS
!===============================================================================

  subroutine run_loafs()
  
    integer :: n_cycles
    integer :: c
  
    if (master) call header("LOAFS SIMULATION", level=1)
  
    call mc_create_sites()
    call sort_sites()
    
    n_cycles = 1
    
    do c = 1, n_cycles
    
!      call loo_process_inputs()
!      call loo_solve()
      
!      call mc_flush_tallies()
      call mc_fixed_source()
      
    end do
    
    stop
  
  end subroutine run_loafs


!===============================================================================
! CREATE_SITES
!===============================================================================

  subroutine mc_create_sites()
  
    integer :: i = 0
  
    ! Allocate particle
    allocate(p)
    allocate(source_site)
  
    ! set work arbitratily high so the site generation loop can run longer
    work = huge(0)
    
    loafs % total_weight = 0.0
    
    ! ====================================================================
    ! LOOP UNTIL ALL SITES GENERATED
    SITE_LOOP: do while (.not. all(loafs % site_bank_idx >= loafs % max_sites))

      p % id = i

      ! set random number seed
      call set_particle_seed(p % id)

      ! get new source particle
      call sample_source_particle()

      ! transport particle
      call transport()

      write(*,*)loafs % site_bank_idx

      i = i + 1

    end do SITE_LOOP
    
  end subroutine mc_create_sites

!===============================================================================
! SORT_SITES
!===============================================================================

  subroutine sort_sites()
  
    integer :: bin
    
    ! ====================================================================
    ! LOOP OVER LOAFS BINS - TODO: parallelize this loop!
    LOAFS_BIN_LOOP: do bin = 1, loafs % n_egroups

      call heap_sort(loafs % site_banks(bin) % sites)

    end do LOAFS_BIN_LOOP

    
  end subroutine sort_sites


!===============================================================================
! MC_FIXED_SOURCE
!===============================================================================

  subroutine mc_fixed_source()
  
    integer :: i, bin
    
    ! ====================================================================
    ! LOOP OVER LOAFS BINS - TODO: parallelize this loop!
    LOAFS_BIN_LOOP: do bin = 1, loafs % n_egroups

      ! ====================================================================
      ! LOOP OVER STARTING SITES IN THIS BIN
      PARTICLE_LOOP: do i = 1, loafs % max_sites(bin)

        write(*,*) loafs % site_banks(bin) % sites(i) % E

      end do PARTICLE_LOOP

    end do LOAFS_BIN_LOOP
  
    
  end subroutine mc_fixed_source

!===============================================================================
! SAMPLE_SOURCE_PARTICLE
!===============================================================================

  subroutine sample_source_particle()

    ! Set particle
    call initialize_particle()

    ! Sample the external source distribution
    call sample_external_source(source_site)

    ! Copy source attributes to the particle
    call copy_source_attributes(source_site)

  end subroutine sample_source_particle

end module loafs_main

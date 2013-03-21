module loafs_main

  use bank_header
  use cmfd_execute, only: cmfd_init_batch, execute_cmfd
  use global
  use loafs_banks,  only: loafs_source_to_particle, loafs_particle_to_bank, &
                          distribute_extra_weight, copy_sites_to_source
  use math,         only: t_percentile
  use output,       only: write_message, header, print_columns, &
                          print_batch_keff, print_generation
  use physics,      only: transport
  use random_lcg,   only: set_particle_seed, prn
  use search,       only: binary_search
  use sort,         only: heap_sort
  use source,       only: initialize_particle, sample_external_source, &
                          copy_source_attributes
  use string,       only: to_str
  use tally,        only: synchronize_tallies, setup_active_usertallies, &
                          reset_result

  type(Bank), pointer :: source_site => null()

contains

!===============================================================================
! RUN_LOAFS
!===============================================================================

  subroutine run_loafs()
  
    integer :: i
  
    if (master) call header("LOAFS SIMULATION", level=1)
  
    ! these normally done in input_xml
    n_batches = 10
    n_inactive = 2
    n_active = n_batches - n_inactive
    current_gen = 1
    gen_per_batch = 1
    allocate(k_batch(n_batches))
    allocate(entropy(n_batches*gen_per_batch))
    entropy = ZERO
    
    ! Allocate particle and source site
    allocate(p)
    allocate(source_site)
    
    ! fool the timers for now so we don't segfault
    call time_inactive % start()
    call time_active % start()
    call time_transport % start()
    call time_bank % start()
    

    ! ==========================================================================
    ! LOOP OVER BATCHES
    BATCH_LOOP: do current_batch = 1, n_batches

      call initialize_batch()

      ! ==========================================================================
      ! LOOP OVER LOAFS BINS - TODO: parallelize this loop
      LOAFS_BIN_LOOP: do loafs_active_bin = loafs % n_egroups, 1, -1

        call initialize_bin()

        ! ========================================================================
        ! RUN ALL PARTICLES IN BANK
        PARTICLE_LOOP: do i=loafs % source_bank_idx(loafs_active_bin), 1, -1

          call loafs_source_to_particle(loafs_active_bin, i)

          call transport()
          
        end do PARTICLE_LOOP

        call finalize_bin()

      end do LOAFS_BIN_LOOP

      call finalize_batch()
      
    end do BATCH_LOOP

    call time_inactive % stop()
    call time_active % stop()
    call time_transport % stop()
    call time_bank % stop()
    
  end subroutine run_loafs

!===============================================================================
! INITIALIZE_BATCH
!===============================================================================

  subroutine initialize_batch()

    message = "Simulating batch " // trim(to_str(current_batch)) // "..."
    call write_message(8)

    call initialize_source()
    
    call cmfd_init_batch()
    
    loafs % site_bank_idx = 0     ! clear loafs site banks
    loafs % total_weights = ZERO

    if (current_batch == n_inactive + 1) then
      ! Switch from inactive batch timer to active batch timer
      call time_inactive % stop()
      call time_active % start()

      ! Enable active batches (and tallies_on if it hasn't been enabled)
      active_batches = .true.
      tallies_on = .true.

      ! Add user tallies to active tallies list
      call setup_active_usertallies()
    end if


  end subroutine initialize_batch


!===============================================================================
! INITIALIZE_BIN
!===============================================================================

  subroutine initialize_bin()
    total_weight = ZERO
  end subroutine initialize_bin

!===============================================================================
! FINALIZE_BIN
!===============================================================================

  subroutine finalize_bin()
    loafs % total_weights(loafs_active_bin) = total_weight
  end subroutine finalize_bin

!===============================================================================
! INITIALIZE_SOURCE
!===============================================================================

  subroutine initialize_source()

    integer :: i, bin

    if (current_batch == 1) then ! sample random starting particles within bin
    
      do bin = 1,loafs % n_egroups
        do i = 1,loafs % max_sites(bin)
        
          ! sample xyz
          call sample_external_source(source_site)
        
          ! override starting energy
          source_site % E = loafs % egrid(bin) + &
                                 prn()*(loafs % egrid(bin+1)-loafs % egrid(bin))
          
          ! insert particle in source bank
          call initialize_particle()
          call copy_source_attributes(source_site)
          call loafs_particle_to_bank(bin,i)
          loafs % source_banks(bin) % sites(i) = &
                                              loafs % site_banks(bin) % sites(i)
          
        end do
        loafs % source_bank_idx(bin) = loafs % max_sites(bin)
      end do
    
    else
    
      call copy_sites_to_source()
    
    end if

  end subroutine initialize_source

!===============================================================================
! FINALIZE_BATCH handles synchronization and accumulation of tallies,
! calculation of Shannon entropy, getting single-batch estimate of keff, and
! turning on tallies when appropriate
!===============================================================================

  subroutine finalize_batch()

    integer :: i, j
    integer :: loafs_bin, filter_bin
    real(8) :: ratio
    
    type(TallyObject),    pointer :: t => null() ! pointer for tally object

    ! Collect tallies
    call time_tallies % start()
    call synchronize_tallies()
    call time_tallies % stop()

    call execute_cmfd()

    write(*,*) current_batch, cmfd % keff
    write(*,*) cmfd % totalxs


  end subroutine finalize_batch


end module loafs_main

module loafs_main

  use bank_header
  use global
  use loafs_banks,  only: loafs_bank_to_particle, loafs_particle_to_bank, &
                          distribute_extra_weight
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
  
    if (master) call header("LOAFS SIMULATION", level=1)
  
    ! these normally done in input_xml
    n_batches = 5
    n_inactive = 1
    n_active = n_batches - n_inactive
    gen_per_batch = 1
    allocate(k_batch(n_batches))
    allocate(entropy(n_batches*gen_per_batch))
    entropy = ZERO
    
    
    ! Allocate particle and source site
    allocate(p)
    allocate(source_site)
    
    ! Display column titles
    call print_columns()

    ! Turn on inactive timer
    call time_inactive % start()
    
    ! ==========================================================================
    ! LOOP OVER BATCHES
    BATCH_LOOP: do current_batch = 1, n_batches

      call initialize_batch()

      ! =======================================================================
      ! LOOP OVER GENERATIONS
      GENERATION_LOOP: do current_gen = 1, gen_per_batch

        call initialize_generation()

        ! ====================================================================
        ! RUN PARTICLES
        
        ! generate bin sites
        call time_transport % start()
        if (current_batch == 1) then
          call mc_create_sites(.false.)
        else
          call mc_create_sites(.true.)
        end if
        call time_transport % stop()
        
        ! sort by starting energy
        call sort_sites()
        
        ! run bin sites
        call time_transport % start()
        call mc_fixed_source()
        call time_transport % stop()
        
        ! Distribute fission bank across processors evenly
        call time_bank % start()
!        call synchronize_bank()
        call time_bank % stop()

        ! Calculate shannon entropy
!        if (entropy_on) call shannon_entropy() ! TODO

        ! Write generation output
        if (master .and. current_gen /= gen_per_batch) call print_generation()

      end do GENERATION_LOOP

      call finalize_batch()

    end do BATCH_LOOP

    call time_active % stop()

  
  end subroutine run_loafs


!===============================================================================
! INITIALIZE_BATCH
!===============================================================================

  subroutine initialize_batch()

    message = "Simulating batch " // trim(to_str(current_batch)) // "..."
    call write_message(8)

    ! Reset total starting particle weight used for normalizing tallies
    total_weight = ZERO

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
! INITIALIZE_GENERATION
!===============================================================================

  subroutine initialize_generation()

    ! TODO

  end subroutine initialize_generation

!===============================================================================
! FINALIZE_BATCH handles synchronization and accumulation of tallies,
! calculation of Shannon entropy, getting single-batch estimate of keff, and
! turning on tallies when appropriate
!===============================================================================

  subroutine finalize_batch()

    ! Collect tallies
    call time_tallies % start()
    call synchronize_tallies()
    call time_tallies % stop()

    ! Collect results and statistics
    call calculate_keff()

    ! Perform CMFD calculation if on
!    if (cmfd_on) call execute_cmfd() ! TODO

    ! Display output
    if (master) call print_batch_keff()
    
!    write(*,*)total_weight, loafs % total_weight
!    write(*,*)loafs % site_bank_idx

  end subroutine finalize_batch


!===============================================================================
! CREATE_SITES
!===============================================================================

  subroutine mc_create_sites(from_fission_bank)
  
    logical, intent(in) :: from_fission_bank
  
    integer :: i
  
    i = 0
  
    ! set work arbitratily high so the site generation loop can run freely
    work = huge(0)
    
    loafs_site_gen = .true.
    
    loafs % total_weight = ZERO
    loafs % extra_weights = ZERO
    
    loafs % site_bank_idx = 0
    
    ! ==========================================================================
    ! LOOP UNTIL ALL SITES GENERATED
    SITE_LOOP: do while (.not. all(loafs % site_bank_idx >= loafs % max_sites))
    
      p % id = i
      
      ! set random number seed
      call set_particle_seed(p % id)

      ! get new source particle
      if (from_fission_bank) then
        call sample_from_fission_bank()
      else
        call sample_source_particle()
      end if

      ! accumulate total starting weight of particles
      loafs % total_weight = loafs % total_weight + p % wgt

      ! transport particle
      loafs_last_ebin = -1 ! force the starting site to be banked
      call transport()

!      write(*,*)p%id,loafs % site_bank_idx

      i = i + 1

    end do SITE_LOOP
    
    loafs_site_gen = .false.

!    write(*,*)loafs % extra_weights
!    write(*,*)loafs % total_weight
    
    call distribute_extra_weight()
    
  end subroutine mc_create_sites

!===============================================================================
! SORT_SITES
!===============================================================================

  subroutine sort_sites()
  
    integer :: bin
    
    ! ==========================================================================
    ! LOOP OVER LOAFS BINS - TODO: parallelize this loop
    LOAFS_BIN_LOOP: do bin = 1, loafs % n_egroups

      call heap_sort(loafs % site_banks(bin) % sites, loafs % site_bank_idx(bin))

    end do LOAFS_BIN_LOOP

    
  end subroutine sort_sites


!===============================================================================
! MC_FIXED_SOURCE
!===============================================================================

  subroutine mc_fixed_source()
    
    integer :: i
    
    n_bank = 0
    
    total_weight = ZERO
    
    ! ==========================================================================
    ! LOOP OVER LOAFS BINS - TODO: parallelize this loop
    LOAFS_BIN_LOOP: do loafs_active_ebin = loafs % n_egroups, 1, -1

!      write(*,*) loafs_active_ebin

      ! ========================================================================
      ! RUN ALL PARTICLES IN BANK
      PARTICLE_LOOP: do i=loafs % site_bank_idx(loafs_active_ebin), 1, -1

        call loafs_bank_to_particle(loafs_active_ebin, i)
        
        call transport()
        
      end do PARTICLE_LOOP

    end do LOAFS_BIN_LOOP
  
  end subroutine mc_fixed_source

!===============================================================================
! MC_FIXED_SOURCE_ORDER - not working yet
!===============================================================================

  subroutine mc_fixed_source_order()
  
    integer :: i(1), ilast(1)
    integer :: bin
    integer :: n_finished
    
    real(8), allocatable :: energies(:)
    
    ! ==========================================================================
    ! LOOP OVER LOAFS BINS - TODO: parallelize this loop
    LOAFS_BIN_LOOP: do bin = loafs % n_egroups, 1, -1

      n_finished = 0
      
      ! TODO: replace the energy list and maxloc with something more efficient
      allocate(energies(loafs % max_sites(bin)))
      energies = loafs % site_banks(bin) % sites % E

      ilast = maxloc(energies)

      ! ========================================================================
      ! RUN ALL PARTICLES IN BANK
      PARTICLE_LOOP: do while (n_finished < loafs % max_sites(bin))

        i = maxloc(energies)

        if (i(1) /= ilast(1)) then
          call loafs_particle_to_bank(bin,ilast(1))
          call loafs_bank_to_particle(bin, i(1))
        end if
        
        call transport(one_collision=.true.)
        
!        write(*,*) p % event
        
        write(*,*) n_finished, i, energies(i(1)),p % E!, p % event

!        if (energies(i(1)) < p % E) then
        if (p % E < 0.06025e-6_8) then
          call transport()
        end if

!        write(*,*) p % alive, p % last_material, p % event

        if (p % alive) then
          energies(i(1)) = p % E
        else
          energies(i(1)) = 0.0_8
          n_finished = n_finished + 1
        end if

        ilast(1) = i(1)

      end do PARTICLE_LOOP

      deallocate(energies)

      stop

    end do LOAFS_BIN_LOOP
  
    
  end subroutine mc_fixed_source_order

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
  

!===============================================================================
! SAMPLE_FROM_FISSION_BANK
!===============================================================================

  subroutine sample_from_fission_bank()

    integer :: i

    ! Set particle
    call initialize_particle()

    i = int(prn()*dble(n_bank))
    do while (i == 0)
      i = int(prn()*dble(n_bank))
    end do

    source_site => fission_bank(i)

    ! Copy source attributes to the particle
    call copy_source_attributes(source_site)

  end subroutine sample_from_fission_bank


!===============================================================================
! CALCULATE_KEFF calculates the single batch estimate of keff as well as the
! mean and standard deviation of the mean for active batches
!===============================================================================

  subroutine calculate_keff()

    real(8) :: temp(2) ! used to reduce sum and sum_sq
    real(8) :: alpha   ! significance level for CI
    real(8) :: t_value ! t-value for confidence intervals

    message = "Calculate batch keff..."
    call write_message(8)

    ! =========================================================================
    ! SINGLE-BATCH ESTIMATE OF K-EFFECTIVE

#ifdef MPI
    if (.not. reduce_tallies) then
      ! Reduce value of k_batch if running in parallel
      if (master) then
        call MPI_REDUCE(MPI_IN_PLACE, k_batch(current_batch), 1, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
      else
        ! Receive buffer not significant at other processors
        call MPI_REDUCE(k_batch(current_batch), temp, 1, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
      end if
    end if
#endif

    ! Normalize single batch estimate of k
    if (master) then
      k_batch(current_batch) = k_batch(current_batch) / &
!           (n_particles * gen_per_batch)
!            loafs % total_weight
            total_weight
    end if

    ! =========================================================================
    ! ACCUMULATED K-EFFECTIVE AND ITS VARIANCE

    if (reduce_tallies) then
      ! In this case, global_tallies has already been reduced, so we don't
      ! need to perform any more reductions and just take the values from
      ! global_tallies directly

      ! Sample mean of keff
      keff = global_tallies(K_TRACKLENGTH) % sum / n_realizations

      if (n_realizations > 1) then
        if (confidence_intervals) then
          ! Calculate t-value for confidence intervals
          alpha = ONE - CONFIDENCE_LEVEL
          t_value = t_percentile(ONE - alpha/TWO, n_realizations - 1)
        else
          t_value = ONE
        end if

        ! Standard deviation of the sample mean of k
        keff_std = t_value * sqrt((global_tallies(K_TRACKLENGTH) % sum_sq / &
             n_realizations - keff * keff) / (n_realizations - 1))
      end if
    else
      ! In this case, no reduce was ever done on global_tallies. Thus, we need
      ! to reduce the values in sum and sum^2 to get the sample mean and its
      ! standard deviation

#ifdef MPI
      call MPI_REDUCE(global_tallies(K_TRACKLENGTH) % sum, temp, 2, &
           MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
#else
      temp(1) = global_tallies(K_TRACKLENGTH) % sum
      temp(2) = global_tallies(K_TRACKLENGTH) % sum_sq
#endif

      ! Sample mean of k
      keff = temp(1) / n_realizations

      if (n_realizations > 1) then
        if (confidence_intervals) then
          ! Calculate t-value for confidence intervals
          alpha = ONE - CONFIDENCE_LEVEL
          t_value = t_percentile(ONE - alpha/TWO, n_realizations - 1)
        else
          t_value = ONE
        end if

        ! Standard deviation of the sample mean of k
        keff_std = t_value * sqrt((temp(2)/n_realizations - keff*keff) / &
             (n_realizations - 1))
      end if
    end if

#ifdef MPI
    ! Broadcast new keff value to all processors
    call MPI_BCAST(keff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
#endif

    ! Reset global tally results
    if (.not. active_batches) then
      call reset_result(global_tallies)
      n_realizations = 0
    end if

  end subroutine calculate_keff

end module loafs_main

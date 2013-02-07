module loafs_main

  use bank_header
  use global
  use loafs_banks,  only: loafs_bank_to_particle, loafs_particle_to_bank, &
                          distribute_extra_weight
  use output,       only: header
  use physics,      only: transport
  use random_lcg,   only: set_particle_seed, prn
  use search,       only: binary_search
  use sort,         only: heap_sort
  use source,       only: initialize_particle, sample_external_source, &
                          copy_source_attributes
  use tally,        only: synchronize_tallies, setup_active_usertallies

  type(Bank), pointer :: source_site => null()

contains

!===============================================================================
! RUN_LOAFS
!===============================================================================

  subroutine run_loafs()
  
    integer :: n_cycles
    integer :: c
  
    if (master) call header("LOAFS SIMULATION", level=1)
  
    n_cycles = 5
    
    do c = 1, n_cycles

!      TODO: run some inactive cycles
    
      if (c == 1) then
        call mc_create_sites(.false.)
      else
        call mc_create_sites(.true.)
      end if
      call sort_sites() ! sort by starting energy
    
!      call loo_process_inputs()
!      call loo_solve()
!      call loo_reweight_sites()
      
!      TODO: conditionally flux tallies

      call mc_fixed_source()
      
      total_weight = loafs % total_weight
      call synchronize_tallies()
      
      write(*,*)global_tallies(K_TRACKLENGTH) % sum / n_realizations
      
    end do
  
  end subroutine run_loafs


!===============================================================================
! CREATE_SITES
!===============================================================================

  subroutine mc_create_sites(from_fission_bank)
  
    logical, intent(in) :: from_fission_bank
  
    integer :: i = 0
  
    ! Allocate particle
    allocate(p)
    allocate(source_site)
  
    ! set work arbitratily high so the site generation loop can run freely
    work = huge(0)
    
    loafs_site_gen = .true.
    
    loafs % total_weight = 0.0_8
    loafs % extra_weights = 0.0_8
    
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

      loafs_last_ebin = -1
      
      ! accumulate total starting weight of particles
      loafs % total_weight = loafs % total_weight + p % wgt

      ! transport particle
      call transport()

!      write(*,*)loafs % site_bank_idx

      i = i + 1

    end do SITE_LOOP
    
    loafs_site_gen = .false.


!    write(*,*)loafs % extra_weights
!    write(*,*)loafs % total_weight
!    total_weight = 0.0_8
    
    
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

      call heap_sort(loafs % site_banks(bin) % sites)

    end do LOAFS_BIN_LOOP

    
  end subroutine sort_sites


!===============================================================================
! MC_FIXED_SOURCE
!===============================================================================

  subroutine mc_fixed_source()
    
    integer :: i
    
    tallies_on = .true.
    call setup_active_usertallies()
    
    n_bank = 0
    
    ! ==========================================================================
    ! LOOP OVER LOAFS BINS - TODO: parallelize this loop
    LOAFS_BIN_LOOP: do loafs_active_ebin = loafs % n_egroups, 1, -1

!      write(*,*) loafs_active_ebin

      ! ========================================================================
      ! RUN ALL PARTICLES IN BANK
      PARTICLE_LOOP: do i=loafs % max_sites(loafs_active_ebin), 1, -1

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

end module loafs_main

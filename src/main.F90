program main

  use constants
  use global
  use initialize,      only: initialize_run
  use intercycle,      only: shannon_entropy, calculate_keff
  use mpi_routines,    only: synchronize_bank
  use output,          only: write_message, header, print_runtime
  use particle_header, only: Particle
  use plot,            only: run_plot
  use physics,         only: transport
  use random_lcg,      only: set_particle_seed
  use source,          only: get_source_particle
  use string,          only: int_to_str
  use tally,           only: synchronize_tallies, write_tallies, &
                             tally_statistics
  use timing,          only: timer_start, timer_stop

#ifdef MPI
  use mpi
#endif

  implicit none

  ! start timer for total run time
  call timer_start(time_total)

  ! set up problem
  call initialize_run()

  ! start problem
  if (plotting) then
     call run_plot()
  else
     call run_problem()

     ! Calculate statistics for tallies and write to tallies.out
     call tally_statistics()
     if (master) call write_tallies()

     ! show timing statistics
     call timer_stop(time_total)
     if (master) call print_runtime()
  end if

  ! deallocate arrays
  call free_memory()
  
contains

!===============================================================================
! RUN_PROBLEM encompasses all the main logic where iterations are performed over
! the cycles and histories.
!===============================================================================

  subroutine run_problem()

    integer                 :: i_cycle     ! cycle index
    integer(8)              :: i_particle  ! history index
    type(Particle), pointer :: p => null()

    if (master) call header("BEGIN SIMULATION", level=1)

    tallies_on = .false.
    call timer_start(time_inactive)

    ! Display column titles
    message = " Cycle   k(cycle)   Entropy         Average k"
    call write_message(1)
    message = " =====   ========   =======    ==================="
    call write_message(1)

    ! ==========================================================================
    ! LOOP OVER CYCLES
    CYCLE_LOOP: do i_cycle = 1, n_cycles

       ! Start timer for computation
       call timer_start(time_compute)

       message = "Simulating cycle " // trim(int_to_str(i_cycle)) // "..."
       call write_message(8)

       ! Set all tallies to zero
       n_bank = 0

       ! =======================================================================
       ! LOOP OVER HISTORIES
       HISTORY_LOOP: do

          ! grab source particle from bank
          p => get_source_particle()
          if ( .not. associated(p) ) then
             ! no particles left in source bank
             exit HISTORY_LOOP
          end if

          ! set random number seed
          i_particle = (i_cycle-1)*n_particles + p % id
          call set_particle_seed(i_particle)
          
          ! set particle trace
          trace = .false.
          if (i_cycle == trace_cycle .and. &
               p % id == trace_particle) trace = .true.

          ! transport particle
          call transport(p)

       end do HISTORY_LOOP

       ! Accumulate time for computation
       call timer_stop(time_compute)

       ! =======================================================================
       ! WRAP UP FISSION BANK AND COMPUTE TALLIES, KEFF, ETC

       ! Start timer for inter-cycle synchronization
       call timer_start(time_intercycle)

       ! Collect tallies
       if (tallies_on) then
          call timer_start(time_ic_tallies)
          call synchronize_tallies()
          call timer_stop(time_ic_tallies)
       end if

       ! Calculate shannon entropy
       if (entropy_on) call shannon_entropy()

       ! Distribute fission bank across processors evenly
       call synchronize_bank(i_cycle)

       ! Collect results and statistics
       call calculate_keff(i_cycle)

       ! print cycle information

       ! Turn tallies on once inactive cycles are complete
       if (i_cycle == n_inactive) then
          tallies_on = .true.
          call timer_stop(time_inactive)
          call timer_start(time_active)
       end if

       ! Stop timer for inter-cycle synchronization
       call timer_stop(time_intercycle)

    end do CYCLE_LOOP

    call timer_stop(time_active)

    ! ==========================================================================
    ! END OF RUN WRAPUP

    if (master) call header("SIMULATION FINISHED", level=1)

  end subroutine run_problem

end program main
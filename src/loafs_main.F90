module loafs_main

  use global
  use output,       only: header
  use physics,      only: transport
  use source,       only: get_source_particle

contains

!===============================================================================
! RUN_LOAFS
!===============================================================================

  subroutine run_loafs()
  
!    integer :: n_cycles
!    integer :: i
  
    if (master) call header("LOAFS SIMULATION", level=1)
  
    call mc_create_sites()
    
    ! communicate sites to nodes
    
!    do i=1,n_cycles
!    
!      call loo_process_inputs()
!      call loo_solve()
!      
!      call mc_flush_tallies()
!      call mc_fixed_source()
!      
!    end do
    
    stop
  
  end subroutine run_loafs


!===============================================================================
! CREATE_SITES
!===============================================================================

  subroutine mc_create_sites()
  
    integer(8) :: i
  
    write(*,*)loafs%egrid
    
    write(*,*)n_particles
    
    ! Allocate particle
    allocate(p)
  
    ! ====================================================================
    ! LOOP OVER PARTICLES
    PARTICLE_LOOP: do i = 1, work

      ! grab source particle from bank
      call get_source_particle(i)

      ! transport particle
      call transport()

      write(*,*)loafs % scatter_bank_idx,n_bank

    end do PARTICLE_LOOP
  
  
  end subroutine mc_create_sites

end module loafs_main

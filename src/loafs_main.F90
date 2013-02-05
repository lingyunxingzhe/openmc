module loafs_main

  use global
  use output,       only: header

contains

!===============================================================================
! RUN_LOAFS
!===============================================================================

  subroutine run_loafs()
  
    integer :: n_cycles
    integer :: i
  
    if (master) call header("LOAFS SIMULATION", level=1)
  
    call mc_create_sites()
    
    ! communicate sites to nodes
    
!    do i=1,n_cycles
!    
!      call loo_inputs()
!      call loo_solve()
!      
!      call mc_flush_tallies()
!      call mc_fixed_source()
!      
!    end do
    
  
  end subroutine run_loafs


!===============================================================================
! CREATE_SITES
!===============================================================================

  subroutine mc_create_sites()
  
  
  end subroutine mc_create_sites

end module loafs_main

module loafs_main

  use global
  use output,       only: header

contains

  subroutine run_loafs()
  
    if (master) call header("LOAFS SIMULATION", level=1)
  
  end subroutine run_loafs


end module loafs_main

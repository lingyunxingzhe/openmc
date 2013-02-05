module loafs_header 

  implicit none
  private

!===============================================================================
! LOAFSBANKSITE -- Used to hold the starting sites for loafs runs. Here we
! require more information than fission banks, since we're banking particles
! that just scattered from one energy group to another, rather than fission
! sites
!===============================================================================

  type, public :: LoafsBankSite

    sequence

    real(8)    :: wgt    ! weight of bank site
    real(8)    :: xyz(3) ! location of bank particle
    real(8)    :: uvw(3) ! diretional cosines
    real(8)    :: E      ! energy
    
    real(8)    :: mu            ! angle of scatter

    ! Pre-collision physical data
    real(8)    :: last_xyz(3)   ! previous coordinates
    real(8)    :: last_wgt      ! pre-collision particle weight
    real(8)    :: last_E        ! pre-collision energy
    real(8)    :: absorb_wgt    ! weight absorbed for survival biasing

    ! What event last took place
    integer    :: event         ! scatter, absorption, fission
    integer    :: event_nuclide ! index in nuclides array
    integer    :: event_MT      ! reaction MT

    ! Post-collision physical data
    integer    :: n_bank        ! number of fission sites banked
    real(8)    :: wgt_bank      ! weight of fission sites banked

    integer    :: cell_born     ! index for cell particle was born in

    ! Statistical data
    integer    :: n_collision   ! # of collisions
    
  end type LoafsBankSite

!===============================================================================
! LOAFSBANK
!===============================================================================

  type, public :: LoafsBank

    type(LoafsBankSite), allocatable :: sites(:)
    
  end type LoafsBank

!===============================================================================
! LOAFS_TYPE
!===============================================================================

  type, public :: loafs_type

    ! energy grid
    real(8), allocatable :: egrid(:)
    
    integer              :: n_egroups
    
    type(LoafsBank), allocatable :: scatter_banks(:)
    integer, allocatable         :: scatter_bank_idx(:)

  end type loafs_type

contains


end module loafs_header

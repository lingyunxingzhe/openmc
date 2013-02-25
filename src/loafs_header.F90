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


    integer(8) :: id

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

    ! Energy grid data
    integer    :: index_grid    ! index on unionized energy grid
    real(8)    :: interp        ! interpolation factor for energy grid

    ! Indices for various arrays
    integer    :: surface       ! index for surface particle is on
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

    real(8), allocatable          :: egrid(:)         ! energy grid

    integer                       :: n_egroups
    
    integer, allocatable          :: max_sites(:)
    
    type(LoafsBank), allocatable  :: site_banks(:)
    integer, allocatable          :: site_bank_idx(:)
    
    type(LoafsBank), allocatable  :: source_banks(:)
    
    real(8)                       :: total_weight
    real(8),allocatable           :: extra_weights(:)
    
    integer                       :: fission_bank_size
    
  end type loafs_type

contains


end module loafs_header

module sort

  use error,          only: fatal_error
  use global,         only: message
  use loafs_header,   only: LoafsBankSite

  ! TODO: implement smooth sort and quicksort for comparison

  interface heap_sort
    module procedure heap_sort_loafs_bank_sites, heap_sort_real8
  end interface heap_sort
  
  interface sift_down
    module procedure sift_down_loafs_bank_sites, sift_down_real8
  end interface sift_down

contains

!===============================================================================
! HEAP_SORT performs an in-place heap sort on an array, guaranteed O(n log n)
! and O(1) in memory.
!===============================================================================

  !=============================================
  ! HEAP_SORT_REAL8
  subroutine heap_sort_real8(array, length)

    real(8), intent(inout)        :: array(:)
    integer, intent(in), optional :: length

    integer                            :: start, n, end_
    real(8)                            :: tmp

    if (present(length)) then
      n = length
    else
      n = size(array)
    end if

    do start = (n - 2) / 2 + 1, 1, -1
      call sift_down(array, start, n)
    end do

    do end_ = n, 2, -1
      tmp = array(1)
      array(1) = array(end_)
      array(end_) = tmp
      call sift_down(array, 1, end_)
    end do
    
  end subroutine heap_sort_real8


  !=============================================
  ! HEAP_SORT_LOAFS_BANK_SITES
  subroutine heap_sort_loafs_bank_sites(array, length)

    type(LoafsBankSite), intent(inout) :: array(:)
    integer, intent(in), optional      :: length

    integer                            :: start, n, end_
    type(LoafsBankSite)                :: tmp

    if (present(length)) then
      n = length
    else
      n = size(array)
    end if

    do start = (n - 2) / 2 + 1, 1, -1
      call sift_down(array, start, n)
    end do

    do end_ = n, 2, -1
      tmp = array(1)
      array(1) = array(end_)
      array(end_) = tmp
      call sift_down(array, 1, end_)
    end do
    
  end subroutine heap_sort_loafs_bank_sites


!===============================================================================
! SIFT_DOWN
!===============================================================================

  !=============================================
  ! SIFT_DOWN_REAL8
  subroutine sift_down_real8(array, start, end_)
   
    real(8), intent(inout)    :: array(:)
    integer, intent(in)       :: start, end_
    
    integer                   :: child, root
    real(8)                   :: tmp
   
    root = start
    do while(root*2 < end_)
      child = root * 2
   
      if ((child + 1 < end_) .and. (array(child) < array(child+1))) then
        child = child + 1
      end if
   
      if (array(root) < array(child)) then
        tmp = array(child)
        array(child) = array(root)
        array(root) = tmp
        root = child
      else
        return
      end if
      
    end do      
   
  end subroutine sift_down_real8
  
  !========================================================
  ! SIFT_DOWN_LOAFS_BANK_SITES
  subroutine sift_down_loafs_bank_sites(array, start, end_)
   
    type(LoafsBankSite), intent(inout)    :: array(:)
    integer, intent(in)                   :: start, end_
    
    integer                   :: child, root
    type(LoafsBankSite)       :: tmp
   
    root = start
    do while(root*2 < end_)
      child = root * 2
   
      if ((child + 1 < end_) .and. (array(child) % E < array(child+1) % E)) then
        child = child + 1
      end if
   
      if (array(root) % E < array(child) % E) then
        tmp = array(child)
        array(child) = array(root)
        array(root) = tmp
        root = child
      else
        return
      end if
      
    end do      
   
  end subroutine sift_down_loafs_bank_sites

end module sort
